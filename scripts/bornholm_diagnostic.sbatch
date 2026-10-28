#!/usr/bin/env bash
#SBATCH --job-name=bornholm-diagnostic
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=slurm-bornholm-diagnostic-%j.out

# Rebuild the corrected ORIGINAL baseline, then solve a short diagnostic.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
stage="${1:-all}"
hours="${2:-24}"
method="${3:-barrier}"
case "$stage" in
  all|solve|dry-run) ;;
  *) echo "Usage: sbatch scripts/bornholm_diagnostic.sbatch [all|solve|dry-run] [hours] [barrier|pdhg]" >&2; exit 2 ;;
esac
[[ "$hours" =~ ^[0-9]+$ ]] && (( hours > 0 && hours < 8760 )) || { echo "Invalid hours: $hours" >&2; exit 2; }
[[ "$method" == barrier || "$method" == pdhg ]] || { echo "Invalid method: $method" >&2; exit 2; }
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo "Submit with sbatch from the repository root." >&2; exit 2; }
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1

if [[ "$stage" == dry-run ]]; then
  bash scripts/bornholm_gas_capacity_rebuild.sbatch dry-run
  exit 0
fi
if [[ "$stage" == all ]]; then
  bash scripts/bornholm_gas_capacity_rebuild.sbatch build
fi
out="results/bornholm_diagnostic/${SLURM_JOB_ID}_${hours}h_${method}"
mkdir -p "$(dirname "$out")"
echo "Diagnostic: ${hours} consecutive hours, ${method}, all 52 zones and full baseline sectors."
echo "One-hour physical storage steps; prorated annual flow budgets and capital costs."
echo "Results: $out"
"$pixi_wrapper" run python -m scripts.run_bornholm_diagnostic \
  --hours "$hours" --method "$method" --threads "${SLURM_CPUS_PER_TASK:-4}" \
  --time-limit 3600 --output-dir "$out"
