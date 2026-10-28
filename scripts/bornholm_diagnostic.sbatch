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

# Rebuild the selected corrected baseline, then solve a short diagnostic.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
stage="${1:-all}"
hours="${2:-24}"
method="${3:-barrier}"
variant="${4:-original}"
case "$stage" in
  all|build|solve|dry-run) ;;
  *) echo "Usage: sbatch scripts/bornholm_diagnostic.sbatch [all|build|solve|dry-run] [hours] [barrier|pdhg] [original|simplified]" >&2; exit 2 ;;
esac
[[ "$variant" == original || "$variant" == simplified ]] || { echo "Unknown baseline: $variant" >&2; exit 2; }
[[ "$hours" =~ ^[0-9]+$ ]] && (( hours > 0 && hours < 8760 )) || { echo "Invalid hours: $hours" >&2; exit 2; }
[[ "$method" == barrier || "$method" == pdhg ]] || { echo "Invalid method: $method" >&2; exit 2; }
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo "Submit with sbatch from the repository root." >&2; exit 2; }
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1

if [[ "$stage" == dry-run ]]; then
  bash scripts/bornholm_gas_capacity_rebuild.sbatch dry-run original
  if [[ "$variant" == simplified ]]; then
    bash scripts/bornholm_gas_capacity_rebuild.sbatch dry-run simplified
  fi
  exit 0
fi
if [[ "$stage" == all || "$stage" == build ]]; then
  bash scripts/bornholm_gas_capacity_rebuild.sbatch build original
  if [[ "$variant" == simplified ]]; then
    bash scripts/bornholm_gas_capacity_rebuild.sbatch build simplified
  fi
fi
if [[ "$variant" == simplified ]]; then
  "$pixi_wrapper" run python scripts/check_bornholm_simplified.py --stage composed
fi
if [[ "$stage" == build ]]; then
  echo "Build and input checks complete. No solver was launched."
  exit 0
fi
out="results/bornholm_diagnostic/${SLURM_JOB_ID}_${variant}_${hours}h_${method}"
mkdir -p "$(dirname "$out")"
echo "Diagnostic: ${variant} baseline, ${hours} consecutive hours, ${method}, all 52 zones."
echo "One-hour physical storage steps; prorated annual flow budgets and capital costs."
echo "Results: $out"
"$pixi_wrapper" run python -m scripts.run_bornholm_diagnostic \
  --variant "$variant" --hours "$hours" --method "$method" --threads "${SLURM_CPUS_PER_TASK:-4}" \
  --time-limit 3600 --output-dir "$out"
