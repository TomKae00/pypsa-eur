#!/usr/bin/env bash
#SBATCH --job-name=bornholm-commitment
#SBATCH --array=0-8%6
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-commitment-%A_%a.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
[[ "${SLURM_JOB_PARTITION:-}" == workq && "${SLURM_CPUS_PER_TASK:-0}" == 16 ]] || exit 2
index="${SLURM_ARRAY_TASK_ID:?Submit array tasks 0-8}"
[[ "$index" =~ ^[0-8]$ ]] || exit 2
if [[ $# == 0 ]]; then
    args=(--index "$index")
elif [[ $# == 1 && "$1" == --check ]]; then
    args=(--check)
else
    printf 'Use the default array or pass --check with --array=0.\n' >&2
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
    "$pixi_wrapper" run python -m scripts.run_bornholm_commitment_stress "${args[@]}"
