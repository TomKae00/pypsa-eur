#!/usr/bin/env bash
#SBATCH --job-name=bornholm-selected
#SBATCH --partition=workq
#SBATCH --array=0-7
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-selected-%A_%a.out

# Replay all eight existing annual 3h scenarios with nf1_aggfill_t16, the
# selected baseline configuration. Reserve the same exclusive 32-CPU node
# used by that candidate; Gurobi itself uses 16 threads.
# Request a 100 GB Gurobi cap, reduced to 75% of available node RAM by the runner.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
[[ "${SLURM_JOB_PARTITION:-}" == workq ]] || { echo 'Scenario replays require workq.' >&2; exit 2; }
[[ "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || { echo 'Allocate 32 CPUs for this array.' >&2; exit 2; }
index="${SLURM_ARRAY_TASK_ID:-}"
[[ "$index" =~ ^[0-7]$ ]] || { echo 'Submit as array tasks 0-7 with sbatch.' >&2; exit 2; }
scenarios=(
    bornholm_zonal_baseline
    bornholm_zonal_radial
    bornholm_zonal_hybrid
    bornholm_zonal_radial_de
    bornholm_zonal_optimised
    bornholm_zonal_open
    bornholm_zonal_baseline_expandable
    bornholm_zonal_open_expandable
)
scenario="${scenarios[$index]}"
campaign="${SLURM_ARRAY_JOB_ID:?Submit with sbatch}"
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export SINGULARITYENV_SLURM_CPUS_PER_TASK="$SLURM_CPUS_PER_TASK"
export SINGULARITYENV_SLURM_JOB_ID="${SLURM_JOB_ID:-}"
export SINGULARITYENV_SLURM_ARRAY_JOB_ID="$campaign"
export SINGULARITYENV_SLURM_ARRAY_TASK_ID="$index"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
    "$pixi_wrapper" run python -m scripts.run_bornholm_accuracy \
    --scenario "$scenario" --profile nf1_aggfill_t16 --threads 16 --soft-mem-limit-gb 100 \
    --output-dir "results/bornholm_accuracy/$campaign/$scenario"
