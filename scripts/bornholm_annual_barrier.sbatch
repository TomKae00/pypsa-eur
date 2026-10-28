#!/usr/bin/env bash
#SBATCH --job-name=bornholm-annual-barrier
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-annual-barrier-%j.out

# Full-year solve of an already rebuilt source, with independent outputs.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
variant="${1:-original}"
case "$variant" in
  original|simplified) ;;
  *) echo "Usage: sbatch scripts/bornholm_annual_barrier.sbatch [original|simplified]" >&2; exit 2 ;;
esac
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo "Submit with sbatch from the repository root." >&2; exit 2; }
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1
if [[ "$variant" == simplified ]]; then
  "$pixi_wrapper" run python scripts/check_bornholm_simplified.py --stage composed
fi
out="results/bornholm_annual_barrier/${SLURM_JOB_ID}_${variant}_8760h_barrier"
echo "Annual barrier: ${variant}, 8760 hours, four threads, no rebuild or annual-input scaling."
echo "24-hour job walltime; 20-hour solver limit; automatic Gurobi memory cap."
echo "Results: $out"
"$pixi_wrapper" run python -m scripts.run_bornholm_diagnostic \
  --variant "$variant" --annual --hours 8760 --method barrier \
  --threads "${SLURM_CPUS_PER_TASK:-4}" --time-limit 72000 --output-dir "$out"
