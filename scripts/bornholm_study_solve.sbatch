#!/usr/bin/env bash
#SBATCH --job-name=bornholm-study-solve
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-study-solve-%A_%a.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
[[ "${SLURM_JOB_PARTITION:-}" == workq ]] || { echo 'Solves must run in workq.' >&2; exit 2; }
[[ "${SLURM_CPUS_PER_TASK:-0}" == 4 ]] || { echo 'Expected four allocated solver CPUs.' >&2; exit 2; }
index="${SLURM_ARRAY_TASK_ID:?Submit with --array=0-7}"
[[ "$index" =~ ^[0-7]$ ]] || exit 2
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1
"$pixi_wrapper" run python -m scripts.run_bornholm_study --index "$index"
