#!/usr/bin/env bash
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00

# Submitted by submit_hindcast_sophia.sh; no explicit --mem on Sophia workq.
set -euo pipefail

: "${SLURM_JOB_ID:?Run through sbatch, not on the login node}"
mode="${1:?Expected weather, build-solve or solve}"
year="${2:?Expected year or array}"
expected_revision="${3:?Expected submitted Git revision}"

if [[ "$year" == array ]]; then
    case "${SLURM_ARRAY_TASK_ID:?Missing array index}" in
        0) year=2022 ;;
        1) year=2024 ;;
        2) year=2025 ;;
        *) echo "Unsupported array index" >&2; exit 2 ;;
    esac
fi
case "$year" in 2022|2023|2024|2025) ;; *) exit 2 ;; esac
case "$mode" in weather|build-solve|solve) ;; *) exit 2 ;; esac

cd "${SLURM_SUBMIT_DIR:?Submit from the repository root}"
actual_revision="$(git rev-parse HEAD)"
if [[ "$actual_revision" != "$expected_revision" ]]; then
    echo "Repository changed since submission. Expected $expected_revision; got $actual_revision." >&2
    echo "Resubmit at the intended revision; do not pull while this campaign is running." >&2
    exit 1
fi
if ! git diff --quiet HEAD --; then
    echo "Tracked source/configuration changes detected; commit them before resubmitting." >&2
    exit 1
fi

cores="${SLURM_CPUS_PER_TASK:-8}"
export TMPDIR="$PWD/tmp/hindcast-${SLURM_JOB_ID}-$year"
mkdir -p "$TMPDIR"
export SINGULARITYENV_TMPDIR="$TMPDIR"
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export SINGULARITYENV_PYTHONUNBUFFERED=1
# Avoid each concurrently scheduled preprocessing task spawning a full BLAS pool.
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1
export SINGULARITYENV_NUMEXPR_NUM_THREADS=1

# Reuse the installed environment without concurrent package/lock-file changes.
run_pixi() {
    "$HOME/bin/pixi-sophia" run --as-is "$@"
}

network="resources/hindcast-build-$year/networks/composed_$year.nc"
output="results/hindcast-build-$year/network.nc"
mapping="data/hindcast/bus_zones_$year.csv"
configs=(config/hindcast/build.yaml "config/hindcast/build-$year.yaml")

printf 'Stage=%s year=%s job=%s node=%s CPUs=%s revision=%s\n' \
    "$mode" "$year" "$SLURM_JOB_ID" "$(hostname)" "$cores" "$actual_revision"

if [[ "$mode" == weather ]]; then
    # This native rule has no upstream inputs; annual cutout outputs are disjoint.
    # Keep Snakemake locks enabled, including when other years run concurrently.
    run_pixi snakemake "data/cutout/build/unknown/europe-$year-era5.nc" \
        --cores "$cores" --configfile "${configs[@]}" --rerun-incomplete
    exit 0
fi

if [[ "$mode" == build-solve ]]; then
    if [[ "$year" == 2025 ]]; then
        echo "2025 build is pending review of annual nuclear availability and other input coverage." >&2
        echo "The configured data/nuclear_p_max_pu.csv currently ends at 2024; weather preparation is allowed." >&2
        exit 1
    fi
    for kind in fuel co2; do
        input="data/hindcast/paper-procedure/${kind}_${year}_hourly.csv"
        if [[ ! -s "$input" ]]; then
            echo "Missing prepared annual input: $input" >&2
            exit 1
        fi
    done
    # Native workflow reuses completed common downloads and the preceding weather job.
    # Run directories differ by year. Any unexpected shared writer is blocked by
    # Snakemake rather than bypassed with --nolock.
    run_pixi snakemake "$network" \
        --cores "$cores" --configfile "${configs[@]}" --rerun-incomplete
fi

if [[ ! -s "$network" ]]; then
    echo "Missing composed network: $network" >&2
    exit 1
fi
if [[ -e "$output" || -e "${output%.nc}.audit.json" ]]; then
    echo "Existing solve output/audit found at $output; refusing to overwrite it." >&2
    exit 1
fi
if [[ ! -f "$mapping" ]]; then
    run_pixi python scripts/hindcast.py prepare-zonal "$network" --year "$year"
fi
# A second invocation checks a new template and fails on unknown/unreviewed IDs.
# Existing reviewed CSVs (including 2023) are preserved.
run_pixi python scripts/hindcast.py prepare-zonal "$network" --year "$year"

mkdir -p "results/hindcast-build-$year"
printf '%s\n' "$actual_revision" > "results/hindcast-build-$year/submitted-revision.txt"
run_pixi python scripts/hindcast.py solve "$network" \
    --output "$output" --solver gurobi --threads "$cores"
printf 'Completed year=%s; solved network=%s\n' "$year" "$output"
