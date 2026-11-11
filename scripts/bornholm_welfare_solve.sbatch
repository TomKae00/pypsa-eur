#!/usr/bin/env bash
#SBATCH --job-name=bornholm-welfare-solve
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-welfare-solve-%A_%a.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
[[ "${SLURM_JOB_PARTITION:-}" == workq ]] || exit 2
[[ "${SLURM_CPUS_PER_TASK:-0}" == 16 ]] || exit 2
index="${SLURM_ARRAY_TASK_ID:?Submit an array with tasks 0-11}"
[[ "$index" =~ ^([0-9]|1[01])$ ]] || exit 2
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export SINGULARITYENV_SLURM_CPUS_PER_TASK="$SLURM_CPUS_PER_TASK"
export SINGULARITYENV_SLURM_JOB_ID="$SLURM_JOB_ID"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1
srun --ntasks=1 --cpus-per-task=16 --cpu-bind=cores \
    "$pixi_wrapper" run python -m scripts.run_bornholm_welfare --index "$index"
