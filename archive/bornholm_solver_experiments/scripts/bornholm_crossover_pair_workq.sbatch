#!/usr/bin/env bash
#SBATCH --job-name=bornholm-crossover-pair
#SBATCH --partition=workq
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-crossover-pair-%A_%a.out

# Same existing annual 3h baseline, NF1 and Crossover=2; only solver threads differ.
# Both exclusive-node tasks reserve 16 CPUs; task 0 uses eight Gurobi threads.
# The runner explicitly retains include_objective_constant=True for these profiles.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
[[ "${SLURM_JOB_PARTITION:-}" == workq ]] || { echo 'Crossover trials require workq.' >&2; exit 2; }
[[ "${SLURM_CPUS_PER_TASK:-0}" == 16 ]] || { echo 'Allocate 16 CPUs for this array.' >&2; exit 2; }
case "${SLURM_ARRAY_TASK_ID:-}" in
    0) profile=nf1_crossover2_t8; threads=8; label=tol8_nf1_crossover2_t8 ;;
    1) profile=nf1_crossover2_t16; threads=16; label=tol8_nf1_crossover2_t16 ;;
    *) echo 'Submit as array tasks 0-1 with sbatch.' >&2; exit 2 ;;
esac
pair="${SLURM_ARRAY_JOB_ID:?Submit with sbatch}"
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1
"$pixi_wrapper" run python -m scripts.run_bornholm_accuracy \
    --profile "$profile" --threads "$threads" \
    --output-dir "results/bornholm_accuracy/$pair/$label"
