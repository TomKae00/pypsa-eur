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
reviewed="hindcast-reviewed-$year"
grid_config=config/hindcast/grid-reviewed-full-years.yaml
wind_args=()
case "$grid" in
    reviewed) run_name="$reviewed" ;;
    wind-calibrated)
        run_name="hindcast-wind-calibrated-$year"
        grid=reviewed
        wind_args=(--wind-config config/hindcast/wind-calibrated-full-years.yaml)
        ;;
    native) run_name="$native" ;;
    *) echo "Expected reviewed, native or wind-calibrated, got: $grid" >&2; exit 2 ;;
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
reviewed_result="results/$reviewed/network.nc"
grid_args=()
if [[ ${#wind_args[@]} -gt 0 && -s "$reviewed_result" && -s "${reviewed_result%.nc}.audit.json" ]]; then
    input="$reviewed_result"
elif [[ "$grid" == reviewed && ( -s "$native_input" || -s "$native_result" ) ]]; then
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
    --threads "$cores" "${grid_args[@]}" "${wind_args[@]}"
if [[ -n "${HINDCAST_SUBMITTED_REVISION:-}" ]]; then
    printf '%s\n' "$HINDCAST_SUBMITTED_REVISION" > "$(dirname "$output")/submitted-revision.txt"
fi
