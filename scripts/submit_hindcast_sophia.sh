#!/usr/bin/env bash
# Human-invoked submission helper; uses host Git/Slurm and never starts Python
# or connects to another machine. Run from your Sophia terminal.
set -euo pipefail

# Reuse a still-known weather array after changing only solve/reporting code.
# Caller must first cancel the old downstream build/solve array.
existing_weather=""
if [[ "$#" -gt 0 ]]; then
    if [[ "$#" -ne 2 || "$1" != --weather-array || ! "$2" =~ ^[0-9]+$ ]]; then
        echo "Usage: bash scripts/submit_hindcast_sophia.sh [--weather-array JOBID]" >&2
        exit 2
    fi
    existing_weather="$2"
fi

repo_dir="$(git rev-parse --show-toplevel)"
cd "$repo_dir"
if [[ "$(git rev-parse --abbrev-ref HEAD)" != feature/historical-hindcast ]]; then
    echo "Use the feature/historical-hindcast checkout." >&2
    exit 1
fi
if ! git diff --quiet HEAD --; then
    echo "Commit tracked source/configuration changes before submitting." >&2
    exit 1
fi
if [[ ! -x "$HOME/bin/pixi-sophia" ]]; then
    echo "Missing executable wrapper: $HOME/bin/pixi-sophia" >&2
    exit 1
fi
if [[ ! -s resources/hindcast-build-2023/networks/composed_2023.nc ]]; then
    echo "Missing composed_2023.nc; this campaign expects the completed 2023 build." >&2
    exit 1
fi
if [[ -e results/hindcast-build-2023/network.nc ]]; then
    echo "2023 result already exists; use individual jobs documented in doc/hindcast-sophia.md." >&2
    exit 1
fi
mkdir -p logs tmp
revision="$(git rev-parse HEAD)"
worker="$repo_dir/scripts/hindcast_sophia.sbatch"
record="logs/hindcast-submissions-$(date +%Y%m%dT%H%M%S)-$$.txt"
printf 'revision=%s\n' "$revision" > "$record"

# No dependency: the 2023 solve starts from the already-composed file.
solve_job="$(sbatch --parsable --job-name=hindcast-solve-2023 \
    --cpus-per-task=8 --output=logs/hindcast-solve-2023-%j.log \
    "$worker" solve 2023 "$revision")"
solve_job="${solve_job%%;*}"
printf 'solve_2023=%s\n' "$solve_job" | tee -a "$record"

# A restart can attach new annual builds to the existing weather array without
# interrupting downloads. Otherwise submit the three independent cutout jobs.
if [[ -n "$existing_weather" ]]; then
    weather_job="$existing_weather"
else
    weather_job="$(sbatch --parsable --job-name=hindcast-weather \
        --cpus-per-task=4 --array=0-2%3 \
        --output=logs/hindcast-weather-%A_%a.log \
        "$worker" weather array "$revision")"
    weather_job="${weather_job%%;*}"
fi
printf 'weather_array=%s (0=2022, 1=2024, 2=2025)\n' "$weather_job" | tee -a "$record"

# aftercorr releases each build as soon as its OWN weather task succeeds.
# The 2022/2024 jobs do not wait for 2025 weather or for the 2023 solve.
annual_job="$(sbatch --parsable --job-name=hindcast-build-solve \
    --cpus-per-task=8 --array=0-1%2 \
    --dependency="aftercorr:$weather_job" \
    --output=logs/hindcast-build-solve-%A_%a.log \
    "$worker" build-solve array "$revision")"
annual_job="${annual_job%%;*}"
printf 'build_solve_array=%s (0=2022, 1=2024)\n' "$annual_job" | tee -a "$record"
printf '%s\n' "2025: weather only; annual model inputs still require review." | tee -a "$record"

printf 'Submission record: %s\n' "$record"
printf '%s\n' "Keep this checkout and environment unchanged while the jobs are queued/running."
printf '%s\n' "Inspect jobs with: squeue -u \"\$USER\""
