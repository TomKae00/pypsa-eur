#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Run in the PyPSA-Eur environment.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

year="${1:?Expected year 2024 or 2025}"
cores="${2:-8}"
configs=(config/hindcast/build.yaml "config/hindcast/build-$year.yaml")
case "$year" in
    2024) configs+=(config/hindcast/build-full-2024.yaml) ;;
    2025) ;;
    *) echo "Unsupported full-year hindcast: $year" >&2; exit 2 ;;
esac

input="resources/hindcast-full-$year/networks/composed_$year.nc"
output="results/hindcast-full-$year/network.nc"
if [[ -e "$output" || -e "${output%.nc}.audit.json" ]]; then
    echo "Existing hindcast output: $output" >&2
    exit 2
fi

costs=data/hindcast/paper-procedure
if [[ ! -s "$costs/fuel_${year}_hourly.csv" || ! -s "$costs/co2_${year}_hourly.csv" ]]; then
    markets="data/hindcast/market-data-$year"
    if [[ ! -s "$markets/manifest.json" ]]; then
        python scripts/retrieve_hindcast_market_data.py --years "$year" --directory "$markets"
    fi
    python scripts/prepare_hindcast_cost_inputs.py --years "$year" \
        --source-directory "$markets" --directory "$costs"
fi

snakemake --executor local --cores "$cores" --resources mem_mb="${HINDCAST_BUILD_MEM_MB:-90000}" \
    --rerun-incomplete --configfile "${configs[@]}" -- "$input"

mapping=(python -m scripts.hindcast prepare-zonal "$input"
    --config config/hindcast/benchmark-full-years.yaml --year "$year")
"${mapping[@]}"
# The first call may have created a new mapping template.
"${mapping[@]}"
python -m scripts.hindcast solve "$input" --output "$output" --solver gurobi --threads "$cores"
