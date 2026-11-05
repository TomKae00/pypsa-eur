#!/usr/bin/env bash
#SBATCH --job-name=bornholm-hybrid-1h
#SBATCH --partition=fatq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=3-00:00:00
#SBATCH --output=slurm-bornholm-hybrid-1h-%j.out

# One node: build the hourly input locally, then solve only fixed DK2-DE Hybrid.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
[[ -n "${SLURM_JOB_ID:-}" && "${SLURM_JOB_PARTITION:-}" == fatq ]] || { echo 'Submit this job to fatq with sbatch.' >&2; exit 2; }
[[ "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || { echo 'Reserve 32 CPUs; solver threads are selected separately.' >&2; exit 2; }
threads="${1:-16}"
case "$threads" in 16|32) ;; *) echo 'Use 16 or 32 solver threads.' >&2; exit 2 ;; esac
memory_profile="${2:-standard}"
case "$memory_profile" in standard|retry|aggfill0) ;; *) echo 'Use standard, retry or aggfill0 as the memory profile.' >&2; exit 2 ;; esac
(( $# <= 3 )) || { echo 'Usage: sbatch scripts/bornholm_hourly_hybrid_fatq.sbatch [16|32] [standard|retry|aggfill0] [reference-job-id]' >&2; exit 2; }
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export SINGULARITYENV_SLURM_CPUS_PER_TASK="$SLURM_CPUS_PER_TASK"
export SINGULARITYENV_SLURM_JOB_ID="$SLURM_JOB_ID"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1

prefix=bornholm_1h_s1_w3_o3_reviewed_biomass30
scenario=bornholm_zonal_hybrid
target="resources/$prefix/$scenario/networks/composed_2035.nc"
directory="results/$prefix/$SLURM_JOB_ID/$scenario"
[[ ! -e "$directory" ]] || { echo "Output exists: $directory" >&2; exit 2; }
mkdir -p "$(dirname "$directory")"
echo "Node $(hostname); allocated CPUs=$SLURM_CPUS_PER_TASK; solver threads=$threads; memory profile=$memory_profile"
head -3 /proc/meminfo

reference_args=()
if [[ "$memory_profile" != standard ]]; then
    reference_job="${3:-}"
    [[ "$reference_job" =~ ^[0-9]+$ ]] || { echo 'Retry requires the previous numeric job ID as argument 3.' >&2; exit 2; }
    reference="results/$prefix/$reference_job/$scenario/summary.json"
    [[ -f "$target" && -f "$reference" ]] || { echo "Missing existing hourly input or reference summary: $reference" >&2; exit 2; }
    reference_args=(--reference-summary "$reference")
    echo "Reusing existing hourly input from job $reference_job; no Snakemake rebuild."
else
    [[ -z "${3:-}" ]] || { echo 'A reference job ID is only used with retry or aggfill0.' >&2; exit 2; }
    # No nested Slurm jobs. The last config layer isolates hourly processed inputs.
    build_memory_mb=$(awk '/^MemAvailable:/ {printf "%d", $2 * 0.75 * 1024 / 1000000}' /proc/meminfo)
    srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
        "$pixi_wrapper" run snakemake --executor local --cores 16 \
        --resources "mem_mb=$build_memory_mb" --rerun-incomplete --latency-wait 120 \
        --configfile config/config.bornholm.yaml config/config.bornholm-3h.yaml \
        config/config.bornholm-calibrated-base.yaml config/config.bornholm-hourly-hybrid.yaml -- "$target"
fi

# Deduct preprocessing from the 70-hour budget, leaving two hours for shutdown
# and export under the default 72-hour walltime. Respect shorter allocations
# when Slurm exposes the actual end timestamp.
time_limit=$(( 252000 - SECONDS ))
if [[ "${SLURM_JOB_END_TIME:-}" =~ ^[0-9]+$ ]]; then
    remaining=$(( SLURM_JOB_END_TIME - $(date +%s) - 7200 ))
    if (( remaining < time_limit )); then time_limit=$remaining; fi
fi
(( time_limit > 60 )) || { echo 'No solve time remains after preprocessing.' >&2; exit 2; }
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
    "$pixi_wrapper" run python -m scripts.run_bornholm_hourly \
    --threads "$threads" --time-limit "$time_limit" --memory-profile "$memory_profile" \
    "${reference_args[@]}" --output-dir "$directory"

