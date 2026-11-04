#!/usr/bin/env bash
#SBATCH --job-name=bornholm-fourway
#SBATCH --partition=workq
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-fourway-%A_%a.out

# Four combinations of homogeneous barrier (-1/1) and the PyPSA objective
# constant (included/omitted). All retain AggFill=-1, PreDual=0, NF1,
# BarConvTol=1e-8 and crossover off, on the same composed annual 3h baseline.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
[[ "${SLURM_JOB_PARTITION:-}" == workq ]] || { echo 'Four-way comparison requires workq.' >&2; exit 2; }
[[ "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || { echo 'Allocate 32 CPUs for this array.' >&2; exit 2; }
case "${SLURM_ARRAY_TASK_ID:-}" in
    0) profile=nf1_aggfill_homoff_conston_t32 ;;
    1) profile=nf1_aggfill_homoff_constoff_t32 ;;
    2) profile=nf1_aggfill_homon_conston_t32 ;;
    3) profile=nf1_aggfill_homon_constoff_t32 ;;
    *) echo 'Submit as array tasks 0-3 with sbatch.' >&2; exit 2 ;;
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
