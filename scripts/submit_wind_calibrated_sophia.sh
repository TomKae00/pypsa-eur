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
previous_job=
for year in 2024 2025; do
    source="$hindcast/results/hindcast-wind-calibrated-$year/network.nc"
    paired_options=(--job-name="bornholm-wind-calibrated-$year" --output=slurm-%x-%j.out)
    if [[ -s "$source" && -s "${source%.nc}.audit.json" ]]; then
        printf '%s calibrated hindcast already saved at %s\n' "$year" "$source"
    else
        hindcast_options=("${options[@]}" --job-name="hindcast-wind-calibrated-$year")
        if [[ -n "$previous_job" ]]; then
            hindcast_options+=(--dependency="afterany:$previous_job")
        fi
        job="$(sbatch --parsable "${hindcast_options[@]}" \
            scripts/hindcast_sophia.sbatch full-year "$year" "$revision" wind-calibrated)"
        job="${job%%;*}"
        printf '%s calibrated hindcast job %s\n' "$year" "$job"
        paired_options+=(--dependency="afterok:$job" --kill-on-invalid-dep=yes)
        previous_job="$job"
    fi
    (
        cd "$bornholm"
        job="$(sbatch --parsable "${paired_options[@]}" \
            "$worker" "$source" --year "$year" --template "$template")"
        printf '%s paired Bornholm job %s\n' "$year" "${job%%;*}"
    )
done
