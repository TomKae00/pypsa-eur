#!/usr/bin/env bash
#SBATCH --job-name=bornholm-full-sector
#SBATCH --partition=fatq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=3-00:00:00
#SBATCH --output=slurm-bornholm-full-sector-%j.out

# Submit from a separate worktree so the active 15-case code hashes stay valid.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
[[ -n "${SLURM_JOB_ID:-}" && "${SLURM_JOB_PARTITION:-}" == fatq ]] || { echo 'Submit to fatq with sbatch.' >&2; exit 2; }
[[ "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || { echo 'Reserve 32 CPUs; the solver uses 16 threads.' >&2; exit 2; }
(( $# == 0 )) || { echo 'This script runs one central Hybrid case with no arguments.' >&2; exit 2; }
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export SINGULARITYENV_SLURM_CPUS_PER_TASK="$SLURM_CPUS_PER_TASK"
export SINGULARITYENV_SLURM_JOB_ID="$SLURM_JOB_ID"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1

prefix=bornholm_2035_full_sector_v1
scenario=central_hybrid
target="resources/$prefix/$scenario/networks/composed_2035.nc"
directory="results/$prefix/$SLURM_JOB_ID/$scenario"
[[ ! -e "$directory" ]] || { echo "Output exists: $directory" >&2; exit 2; }
echo "Node $(hostname); full sector central Hybrid; 3h full year; 16 solver threads"
head -3 /proc/meminfo

# Native preprocessing, local to this allocation; no nested Slurm submission.
build_memory_mb=$(awk '/^MemAvailable:/ {printf "%d", $2 * 0.75 * 1024 / 1000000}' /proc/meminfo)
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
    "$pixi_wrapper" run snakemake --executor local --cores 4 \
    --resources "mem_mb=$build_memory_mb" --rerun-incomplete --latency-wait 120 \
    --configfile config/config.bornholm.yaml config/config.bornholm-3h.yaml \
    config/config.bornholm-calibrated-base.yaml config/config.bornholm-2035.yaml \
    config/config.bornholm-full-sector.yaml -- "$target"

# Count build time against the 70h solve budget; reserve 2h for export/shutdown.
time_limit=$(( 252000 - SECONDS ))
if [[ "${SLURM_JOB_END_TIME:-}" =~ ^[0-9]+$ ]]; then
    remaining=$(( SLURM_JOB_END_TIME - $(date +%s) - 7200 ))
    if (( remaining < time_limit )); then time_limit=$remaining; fi
fi
(( time_limit > 60 )) || { echo 'No solve time remains after preprocessing.' >&2; exit 2; }
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
    "$pixi_wrapper" run python -m scripts.run_bornholm_full_sector \
    --time-limit "$time_limit" --output-dir "$directory"
