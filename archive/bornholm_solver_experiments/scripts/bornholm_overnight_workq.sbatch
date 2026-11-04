#!/usr/bin/env bash
#SBATCH --job-name=bornholm-overnight
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-overnight-%j.out

# One final NF1 candidate: automatic AggFill, PreDual=0 and standard barrier.
# The runner omits the objective-constant variable and records a comparable cost.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
[[ "${SLURM_JOB_PARTITION:-}" == workq ]] || { echo 'Overnight candidate requires workq.' >&2; exit 2; }
[[ "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || { echo 'Allocate 32 CPUs for this run.' >&2; exit 2; }
trial="${SLURM_JOB_ID:?Submit with sbatch}"
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export SINGULARITYENV_SLURM_CPUS_PER_TASK="$SLURM_CPUS_PER_TASK"
export SINGULARITYENV_SLURM_JOB_ID="$trial"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
    "$pixi_wrapper" run python -m scripts.run_bornholm_accuracy \
    --profile nf1_final_t32 --threads 32 \
    --output-dir "results/bornholm_accuracy/$trial/tol8_nf1_final_t32"
