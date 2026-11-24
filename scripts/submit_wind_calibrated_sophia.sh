#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Submit from the hindcast checkout on Sophia.
set -euo pipefail

hindcast="$(git rev-parse --show-toplevel)"
cd "$hindcast"
bornholm="$(cd "${1:-$HOME/projects/pypsa-eur-bornholm}" && pwd -P)"
template="$bornholm/resources/bornholm_2035_final_v1/central_hybrid/networks/composed_2035.nc"
worker="$bornholm/scripts/bornholm_historical.sbatch"
for path in "$template" "$worker" \
    data/hindcast/observed/danish_wind_{2024,2025}.csv.gz; do
    if [[ ! -s "$path" ]]; then
        echo "Missing input: $path" >&2
        exit 1
    fi
done

revision="$(git rev-parse HEAD)"
options=(--cpus-per-task=16 --exclusive --mem=0 --time=48:00:00 --output=slurm-%x-%j.out)
previous=()
for year in 2024 2025; do
    source="$hindcast/results/hindcast-wind-calibrated-$year/network.nc"
    paired_dep=()
    if [[ -s "$source" && -s "${source%.nc}.audit.json" ]]; then
        printf '%s calibrated hindcast already saved at %s\n' "$year" "$source"
    else
        job="$(sbatch --parsable "${options[@]}" "${previous[@]}" \
            --job-name="hindcast-wind-calibrated-$year" \
            scripts/hindcast_sophia.sbatch full-year "$year" "$revision" wind-calibrated)"
        job="${job%%;*}"
        printf '%s calibrated hindcast job %s\n' "$year" "$job"
        paired_dep=(--dependency="afterok:$job" --kill-on-invalid-dep=yes)
        previous=(--dependency="afterany:$job")
    fi
    (
        cd "$bornholm"
        job="$(sbatch --parsable "${paired_dep[@]}" --job-name="bornholm-wind-calibrated-$year" \
            --output=slurm-%x-%j.out "$worker" "$source" --year "$year" --template "$template")"
        printf '%s paired Bornholm job %s\n' "$year" "${job%%;*}"
    )
done
