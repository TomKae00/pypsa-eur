#!/usr/bin/env bash
#SBATCH --job-name=bornholm-pdhg-relax
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-pdhg-relax-%j.out

# Solve-only test: reuse the completed, unchanged simplified composed network.
# NEVER recomposes the network or overwrites results for the ordinary simplified run.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

mode="${1:-solve}"
if [[ "$mode" != "solve" && "$mode" != "dry-run" ]]; then
  echo "Usage: sbatch scripts/bornholm_simplified_solver_experiment.sbatch [solve|dry-run]" >&2
  exit 2
fi

base="config/config.bornholm.yaml"
simplified="config/config.bornholm-simplified.yaml"
variant="config/config.bornholm-simplified-pdhg-relaxed.yaml"
for cfg in "$base" "$simplified" "$variant"; do test -s "$cfg"; done

source_nc="resources/bornholm_simplified_baseline/networks/composed_2035.nc"
test -s "$source_nc" || { echo "Missing composed simplified baseline: $source_nc" >&2; exit 1; }

case_name="bornholm_simplified_pdhg_relaxed"
experiment_nc="resources/${case_name}/networks/composed_2035.nc"
solved_nc="results/${case_name}/networks/solved_2035.nc"
mkdir -p "$(dirname "$experiment_nc")"
if [[ -e "$experiment_nc" && ! -L "$experiment_nc" ]]; then
  echo "Refusing to overwrite existing non-symlink $experiment_nc" >&2
  exit 1
fi
ln -sfn "$(realpath "$source_nc")" "$experiment_nc"

pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1

cores="${SLURM_CPUS_PER_TASK:-16}"
echo "Solver experiment: $case_name"
echo "Network: $source_nc"
echo "Output: $solved_nc"
echo "PDHGConvTol=1e-3; PDHGRelTol=1e-5; no crossover"
echo "Original simplified solve is unaffected."

cmd=("$pixi_wrapper" run snakemake
  --cores "$cores"
  --resources mem_mb=110000
  --rerun-incomplete
  --configfile "$base" "$simplified" "$variant"
  --allowed-rules solve_network)
if [[ "$mode" == "dry-run" ]]; then
  "${cmd[@]}" --dry-run -- "$solved_nc"
else
  "${cmd[@]}" -- "$solved_nc"
fi
