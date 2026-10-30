#!/usr/bin/env bash
#SBATCH --job-name=bornholm-solver-sweep
#SBATCH --partition=workq
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-solver-sweep-%A_%a.out

# Four independent annual replays of the same existing 3h baseline. No rebuild.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
[[ "${SLURM_JOB_PARTITION:-}" == workq ]] || { echo 'Solver trials require workq.' >&2; exit 2; }
[[ "${SLURM_CPUS_PER_TASK:-0}" == 8 ]] || { echo 'Allocate eight solver CPUs.' >&2; exit 2; }
case "${SLURM_ARRAY_TASK_ID:-}" in
    0) profile=nf1; label=tol8_nf1_t8 ;;
    1) profile=nf1_scale2; label=tol8_nf1_scale2_t8 ;;
    2) profile=nf1_crossover; label=tol8_nf1_crossover_t8 ;;
    3) profile=nf1_tol7; label=tol7_nf1_t8 ;;
    *) echo 'Submit as array tasks 0-3 with sbatch.' >&2; exit 2 ;;
esac
sweep="${SLURM_ARRAY_JOB_ID:?Submit with sbatch}"
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1
"$pixi_wrapper" run python -m scripts.run_bornholm_accuracy \
    --profile "$profile" --threads 8 \
    --output-dir "results/bornholm_accuracy/$sweep/$label"
