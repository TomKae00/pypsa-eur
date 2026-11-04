#!/usr/bin/env bash
#SBATCH --job-name=bornholm-saved-audit
#SBATCH --partition=workq
#SBATCH --array=0-7
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:30:00
#SBATCH --output=slurm-bornholm-saved-audit-%A_%a.out

# Inspect eight saved candidates, or the original study with argument "original".
# No optimization model or solver is called.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
[[ "${SLURM_JOB_PARTITION:-}" == workq ]] || { echo 'Saved audits require workq.' >&2; exit 2; }
[[ "${SLURM_CPUS_PER_TASK:-0}" == 1 ]] || { echo 'Allocate one CPU per audit task.' >&2; exit 2; }
campaign="${1:-4930796}"
if [[ "$campaign" == original ]]; then
    input_root=results/bornholm_3h_s1_w2_o3
elif [[ "$campaign" =~ ^[0-9]+$ ]]; then
    input_root="results/bornholm_accuracy/$campaign"
else
    echo 'Pass a numeric solve array ID or original.' >&2
    exit 2
fi
index="${SLURM_ARRAY_TASK_ID:-}"
[[ "$index" =~ ^[0-7]$ ]] || { echo 'Submit as array tasks 0-7 with sbatch.' >&2; exit 2; }
scenarios=(
    bornholm_zonal_baseline
    bornholm_zonal_radial
    bornholm_zonal_hybrid
    bornholm_zonal_radial_de
    bornholm_zonal_optimised
    bornholm_zonal_open
    bornholm_zonal_baseline_expandable
    bornholm_zonal_open_expandable
)
scenario="${scenarios[$index]}"
audit_id="${SLURM_ARRAY_JOB_ID:?Submit with sbatch}"
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1
srun --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
    "$pixi_wrapper" run --as-is python -m scripts.audit_bornholm_candidate \
    "$input_root/$scenario" \
    --output "results/bornholm_saved_audit/$campaign/$audit_id/$scenario"
