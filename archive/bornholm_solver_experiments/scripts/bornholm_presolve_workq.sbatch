#!/usr/bin/env bash
#SBATCH --job-name=bornholm-presolve
#SBATCH --partition=workq
#SBATCH --array=0-2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-presolve-%A_%a.out

# Existing annual 3h baseline, NF1, crossover off, objective constant included.
# Tasks 0/1 change PreDual/AggFill independently against the existing 16-thread
# reference. Task 2 matches task 0 with Gurobi Threads=0 (automatic, normally <=32).
# Reserve 32 CPUs for every task and bind the solver step to its allocation.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
[[ "${SLURM_JOB_PARTITION:-}" == workq ]] || { echo 'Presolve trials require workq.' >&2; exit 2; }
[[ "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || { echo 'Allocate 32 CPUs for this array.' >&2; exit 2; }
case "${SLURM_ARRAY_TASK_ID:-}" in
    0) profile=nf1_predual_t16; threads=16; label=tol8_nf1_predual_t16 ;;
    1) profile=nf1_aggfill_t16; threads=16; label=tol8_nf1_aggfill_t16 ;;
    2) profile=nf1_predual_auto; threads=0; label=tol8_nf1_predual_tauto ;;
    *) echo 'Submit as array tasks 0-2 with sbatch.' >&2; exit 2 ;;
esac
pair="${SLURM_ARRAY_JOB_ID:?Submit with sbatch}"
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export SINGULARITYENV_SLURM_CPUS_PER_TASK="$SLURM_CPUS_PER_TASK"
export SINGULARITYENV_SLURM_JOB_ID="${SLURM_JOB_ID:-}"
export SINGULARITYENV_SLURM_ARRAY_JOB_ID="$pair"
export SINGULARITYENV_SLURM_ARRAY_TASK_ID="$SLURM_ARRAY_TASK_ID"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
    "$pixi_wrapper" run python -m scripts.run_bornholm_accuracy \
    --profile "$profile" --threads "$threads" \
    --output-dir "results/bornholm_accuracy/$pair/$label"
