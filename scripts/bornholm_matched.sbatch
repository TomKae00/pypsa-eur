#!/usr/bin/env bash
#SBATCH --job-name=bornholm-matched
#SBATCH --array=0-2
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-matched-%A_%a.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
[[ "${SLURM_JOB_PARTITION:-}" == workq && "${SLURM_CPUS_PER_TASK:-0}" == 16 ]] || exit 2
index="${SLURM_ARRAY_TASK_ID:?Submit as an array}"
if [[ $# == 0 && "$index" =~ ^[0-2]$ ]]; then
    args=("$index")
elif [[ $# == 1 && "$1" == --layout-screen && "$index" =~ ^([0-9]|1[0-4])$ ]]; then
    args=(--layout-screen "$index")
elif [[ $# == 2 && "$1" == --layout-screen && "$2" == --check && "$index" == 0 ]]; then
    args=(--layout-screen --check)
else
    printf 'Use the default matched array or --array=0-14%%6 with --layout-screen.\n' >&2
    exit 2
fi
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export SINGULARITYENV_SLURM_CPUS_PER_TASK="$SLURM_CPUS_PER_TASK"
export SINGULARITYENV_SLURM_ARRAY_JOB_ID="$SLURM_ARRAY_JOB_ID"
export SINGULARITYENV_SLURM_JOB_ID="$SLURM_JOB_ID"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1
srun --ntasks=1 --cpus-per-task=16 --cpu-bind=cores \
    "$pixi_wrapper" run python -m scripts.run_bornholm_matched "${args[@]}"
