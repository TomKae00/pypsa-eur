#!/usr/bin/env bash
#SBATCH --job-name=bornholm-dual
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-dual-%j.out

# Solver-only dual-simplex run using the existing simplified 2035 network.
# No Gurobi time limit; 24h is the HPC scheduler walltime, as for the baseline.
# Reserve the node exclusively for its memory; dual simplex uses 1 CPU thread.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

mode="${1:-solve}"
if [[ "$mode" != "solve" && "$mode" != "dry-run" ]]; then
  echo "Usage: sbatch scripts/bornholm_simplified_dual_simplex.sbatch [solve|dry-run]" >&2
  exit 2
fi

base="config/config.bornholm.yaml"
simplified="config/config.bornholm-simplified.yaml"
variant="config/config.bornholm-simplified-dual-simplex.yaml"
for cfg in "$base" "$simplified" "$variant"; do
  test -s "$cfg" || { echo "Missing config: $cfg" >&2; exit 1; }
done

source_nc="resources/bornholm_simplified_baseline/networks/composed_2035.nc"
test -s "$source_nc" || { echo "Missing composed network: $source_nc" >&2; exit 1; }

run_name="bornholm_simplified_dual_simplex"
experiment_nc="resources/${run_name}/networks/composed_2035.nc"
solved_nc="results/${run_name}/networks/solved_2035.nc"
mkdir -p "$(dirname "$experiment_nc")"
if [[ -e "$experiment_nc" && ! -L "$experiment_nc" ]]; then
  echo "Refusing to overwrite non-symlink: $experiment_nc" >&2
  exit 1
fi
ln -sfn "$(realpath "$source_nc")" "$experiment_nc"

pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper" || { echo "Missing pixi wrapper: $pixi_wrapper" >&2; exit 1; }
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1

echo "Run: $run_name"
echo "Input: $source_nc"
echo "Output: $solved_nc"
echo "Gurobi: Method=1, Threads=1, no solver TimeLimit"
echo "Slurm walltime: 24h"

cmd=("$pixi_wrapper" run snakemake
  --cores 1
  --resources mem_mb=110000
  --rerun-incomplete
  --configfile "$base" "$simplified" "$variant"
  --allowed-rules solve_network)

if [[ "$mode" == "dry-run" ]]; then
  "${cmd[@]}" --dry-run -- "$solved_nc"
else
  "${cmd[@]}" -- "$solved_nc"
fi
