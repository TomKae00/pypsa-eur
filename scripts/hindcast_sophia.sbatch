#!/usr/bin/env bash
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00

# Submitted by submit_hindcast_sophia.sh; no explicit --mem on Sophia workq.
set -euo pipefail

: "${SLURM_JOB_ID:?Run through sbatch, not on the login node}"
mode="${1:?Expected weather, nominal-build, build-solve, rebuild-solve or solve}"
year="${2:?Expected year or array}"
submitted_revision="${3:?Expected submitted Git revision}"

if [[ "$year" == array ]]; then
    case "${SLURM_ARRAY_TASK_ID:?Missing array index}" in
        0) year=2022 ;;
        1) year=2024 ;;
        2) year=2025 ;;
        *) echo "Unsupported array index" >&2; exit 2 ;;
    esac
fi
case "$year" in 2022|2023|2024|2025) ;; *) exit 2 ;; esac
case "$mode" in weather|nominal-build|build-solve|rebuild-solve|solve) ;; *) exit 2 ;; esac

cd "${SLURM_SUBMIT_DIR:?Submit from the repository root}"
# Host Git exists on the login node but is absent on Sophia compute nodes.
# The submit helper checks the branch/working tree and passes its revision.
# Record that provenance without requiring Git here. Keep source unchanged
# while queued/running; workers cannot independently verify the checkout.

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
rebuild_args=()
flow_args=()
if [[ "$mode" == rebuild-solve ]]; then
    output="results/hindcast-corrected-$year/network.nc"
    # Reconstruct spatial assignments and all dependent resources, retaining
    # the completed raw-data/weather cache. Never solve the stale composed file.
    rebuild_args=(--forcerun base_network compose_network)
    flow_args=(--flow-sources "data/hindcast/flow-sources-corrected-$year.yaml")
fi
mapping="data/hindcast/bus_zones_$year.csv"
configs=(config/hindcast/build.yaml "config/hindcast/build-$year.yaml")
if [[ "$mode" == nominal-build ]]; then
    # Isolate physical-grid preparation by year; no composed dispatch network,
    # weather preparation or solver is requested in this mode.
    nominal_run="hindcast-nominal-$year"
    nominal_config="$TMPDIR/nominal-run.yaml"
    printf 'run:\n  name: %s\n' "$nominal_run" > "$nominal_config"
    configs+=(config/hindcast/line-ratings.yaml "$nominal_config")
    nominal_network="resources/$nominal_run/networks/base.nc"
    nominal_audit="resources/$nominal_run/line_rating_audit.json"
    printf 'Stage=nominal-build year=%s job=%s revision=%s\n' "$year" "$SLURM_JOB_ID" "$submitted_revision"
    run_pixi snakemake "$nominal_network" "$nominal_audit" \
        --cores "$cores" --configfile "${configs[@]}" --rerun-incomplete
    printf 'Completed physical nominal ratings: %s\n' "$nominal_audit"
    exit 0
fi
if [[ "$mode" != weather && ( -e "$output" || -e "${output%.nc}.audit.json" ) ]]; then
    echo "Existing solve output/audit found at $output; refusing to overwrite it." >&2
    exit 1
fi

printf 'Stage=%s year=%s job=%s node=%s CPUs=%s revision=%s\n' \
    "$mode" "$year" "$SLURM_JOB_ID" "$(hostname)" "$cores" "$submitted_revision"

if [[ "$mode" == weather ]]; then
    # This native rule has no upstream inputs; annual cutout outputs are disjoint.
    # Keep Snakemake locks enabled, including when other years run concurrently.
    run_pixi snakemake "data/cutout/build/unknown/europe-$year-era5.nc" \
        --cores "$cores" --configfile "${configs[@]}" --rerun-incomplete
    exit 0
fi

if [[ "$mode" == build-solve || "$mode" == rebuild-solve ]]; then
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
        --cores "$cores" --configfile "${configs[@]}" --rerun-incomplete \
        "${rebuild_args[@]}"
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
    run_pixi python scripts/hindcast.py prepare-zonal "$network" --year "$year" "${flow_args[@]}"
fi
# A second invocation checks a new template and fails on unknown/unreviewed IDs.
# Existing reviewed CSVs (including 2023) are preserved.
run_pixi python scripts/hindcast.py prepare-zonal "$network" --year "$year" "${flow_args[@]}"

mkdir -p "$(dirname "$output")"
printf '%s\n' "$submitted_revision" > "$(dirname "$output")/submitted-revision.txt"
run_pixi python scripts/hindcast.py solve "$network" \
    --output "$output" --solver gurobi --threads "$cores"
printf 'Completed year=%s; solved network=%s\n' "$year" "$output"
