#!/usr/bin/env bash
#SBATCH --job-name=bornholm-objective-pair
#SBATCH --partition=workq
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-objective-pair-%A_%a.out

# Matched final NF1 candidates. Only include_objective_constant differs:
# task 0 includes it; task 1 omits it and retains the comparable cost correction.
# Both use 32 threads, AggFill=-1, PreDual=0, BarHomogeneous=0,
# crossover off and BarConvTol=1e-8 on the same composed annual 3h baseline.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
[[ "${SLURM_JOB_PARTITION:-}" == workq ]] || { echo 'Objective pair requires workq.' >&2; exit 2; }
[[ "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || { echo 'Allocate 32 CPUs for this array.' >&2; exit 2; }
case "${SLURM_ARRAY_TASK_ID:-}" in
    0) profile=nf1_final_conston_t32 ;;
    1) profile=nf1_final_t32 ;;
    *) echo 'Submit as array tasks 0-1 with sbatch.' >&2; exit 2 ;;
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
    --profile "$profile" --threads 32 \
    --output-dir "results/bornholm_accuracy/$pair/tol8_$profile"
