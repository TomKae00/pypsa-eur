#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Run in the PyPSA-Eur environment.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

year="${1:?Expected year 2024 or 2025}"
cores="${2:-8}"
grid="${3:-reviewed}"
configs=(config/hindcast/build.yaml "config/hindcast/build-$year.yaml")
case "$year" in
    2024) configs+=(config/hindcast/build-full-2024.yaml) ;;
    2025) ;;
    *) echo "Unsupported full-year hindcast: $year" >&2; exit 2 ;;
esac

native="hindcast-full-$year"
grid_config=config/hindcast/grid-reviewed-full-years.yaml
case "$grid" in
    reviewed) run_name="hindcast-reviewed-$year" ;;
    native) run_name="$native" ;;
    *) echo "Expected a reviewed or native grid, got: $grid" >&2; exit 2 ;;
esac

input="resources/$run_name/networks/composed_$year.nc"
output="results/$run_name/network.nc"
if [[ -e "$output" || -e "${output%.nc}.audit.json" ]]; then
    echo "Existing hindcast output: $output" >&2
    exit 2
fi

# Reuse a composed or solved native-grid network when available. All
# year-specific generation, load and cost inputs are then identical by construction.
native_input="resources/$native/networks/composed_$year.nc"
native_result="results/$native/network.nc"
grid_args=()
if [[ "$grid" == reviewed && ( -s "$native_input" || -s "$native_result" ) ]]; then
    input="$native_input"
    if [[ ! -s "$input" ]]; then
        input="$native_result"
    fi
    grid_args=(--grid-config "$grid_config")
else
    build_network=true
    if [[ "$grid" == reviewed ]]; then
        configs+=("$grid_config")
    fi
fi

if [[ "${build_network:-false}" == true ]]; then
    costs=data/hindcast/paper-procedure
    if [[ ! -s "$costs/fuel_${year}_hourly.csv" || ! -s "$costs/co2_${year}_hourly.csv" ]]; then
        markets="data/hindcast/market-data-$year"
        if [[ ! -s "$markets/manifest.json" ]]; then
            python scripts/retrieve_hindcast_market_data.py --years "$year" --directory "$markets"
        fi
        python scripts/prepare_hindcast_cost_inputs.py --years "$year" \
            --source-directory "$markets" --directory "$costs"
    fi
    # Pass the run name explicitly because the reviewed overlay serves two years.
    snakemake --executor local --cores "$cores" --resources mem_mb="${HINDCAST_BUILD_MEM_MB:-90000}" \
        --rerun-incomplete --configfile "${configs[@]}" --config "run={name: $run_name}" -- "$input"
fi

mapping=(python -m scripts.hindcast prepare-zonal "$input"
    --config config/hindcast/benchmark-full-years.yaml --year "$year")
"${mapping[@]}"
# The first call may have created a new mapping template.
"${mapping[@]}"
python -m scripts.hindcast solve "$input" --output "$output" --solver gurobi \
    --threads "$cores" "${grid_args[@]}"
if [[ -n "${HINDCAST_SUBMITTED_REVISION:-}" ]]; then
    printf '%s\n' "$HINDCAST_SUBMITTED_REVISION" > "$(dirname "$output")/submitted-revision.txt"
fi
