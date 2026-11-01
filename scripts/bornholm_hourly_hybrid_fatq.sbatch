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
case "$threads" in 16|32) ;; *) echo 'Usage: sbatch scripts/bornholm_hourly_hybrid_fatq.sbatch [16|32]' >&2; exit 2 ;; esac
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export SINGULARITYENV_SLURM_CPUS_PER_TASK="$SLURM_CPUS_PER_TASK"
export SINGULARITYENV_SLURM_JOB_ID="$SLURM_JOB_ID"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1

prefix=bornholm_1h_s1_w2_o3
scenario=bornholm_zonal_hybrid
target="resources/$prefix/$scenario/networks/composed_2035.nc"
directory="results/$prefix/$SLURM_JOB_ID/$scenario"
[[ ! -e "$directory" ]] || { echo "Output exists: $directory" >&2; exit 2; }
mkdir -p "$(dirname "$directory")"
echo "Node $(hostname); allocated CPUs=$SLURM_CPUS_PER_TASK; solver threads=$threads"
head -3 /proc/meminfo

# No nested Slurm jobs: use the existing container wrapper inside this allocation.
# Resource classes come from the 3h layer; the last layer removes temporal
# averaging and isolates processed inputs for the fixed Hybrid's fresh namespace.
build_memory_mb=$(awk '/^MemAvailable:/ {printf "%d", $2 * 0.75 * 1024 / 1000000}' /proc/meminfo)
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
    "$pixi_wrapper" run snakemake --executor local --cores 16 \
    --resources "mem_mb=$build_memory_mb" --rerun-incomplete --latency-wait 120 \
    --configfile config/config.bornholm.yaml config/config.bornholm-3h.yaml \
    config/config.bornholm-hourly-hybrid.yaml -- "$target"

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
    --threads "$threads" --time-limit "$time_limit" --output-dir "$directory"
