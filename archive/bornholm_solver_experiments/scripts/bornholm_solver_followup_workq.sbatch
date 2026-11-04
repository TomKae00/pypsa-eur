#!/usr/bin/env bash
#SBATCH --job-name=bornholm-solver-followup
#SBATCH --partition=workq
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-solver-followup-%A_%a.out

# Four additional solves; do not resubmit or modify the running eight-thread sweep.
# All tasks reserve 16 CPUs on an exclusive node. Tasks 2/3 use eight solver threads
# deliberately, to isolate crossover and NumericFocus from thread-count changes.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
[[ "${SLURM_JOB_PARTITION:-}" == workq ]] || { echo 'Solver trials require workq.' >&2; exit 2; }
[[ "${SLURM_CPUS_PER_TASK:-0}" == 16 ]] || { echo 'Allocate 16 CPUs for this array.' >&2; exit 2; }
case "${SLURM_ARRAY_TASK_ID:-}" in
    0) profile=nf1_t16; threads=16; label=tol8_nf1_t16 ;;
    1) profile=nf1_homogeneous_t16; threads=16; label=tol8_nf1_homogeneous_t16 ;;
    2) profile=nf1_homogeneous_crossover; threads=8; label=tol8_nf1_homogeneous_crossover_t8 ;;
    3) profile=nf0_homogeneous; threads=8; label=tol8_nf0_homogeneous_t8 ;;
    *) echo 'Submit as array tasks 0-3 with sbatch.' >&2; exit 2 ;;
esac
followup="${SLURM_ARRAY_JOB_ID:?Submit with sbatch}"
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1
"$pixi_wrapper" run python -m scripts.run_bornholm_accuracy \
    --profile "$profile" --threads "$threads" \
    --output-dir "results/bornholm_accuracy/$followup/$label"
