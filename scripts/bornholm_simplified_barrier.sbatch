#!/usr/bin/env bash
#SBATCH --job-name=bornholm-barrier
#SBATCH --partition=fatq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-barrier-%j.out

# Solver-only barrier comparison on exactly the same simplified 2035 network.
# fatq provides higher-memory nodes; --mem=0 requests all node memory.
# Does NOT rebuild the model or overwrite any other solver experiment.
# Barrier may still exceed available memory due to its matrix factorization.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

mode="${1:-solve}"
if [[ "$mode" != "solve" && "$mode" != "dry-run" ]]; then
  echo "Usage: sbatch scripts/bornholm_simplified_barrier.sbatch [solve|dry-run]" >&2
  exit 2
fi

base="config/config.bornholm.yaml"
simplified="config/config.bornholm-simplified.yaml"
variant="config/config.bornholm-simplified-barrier.yaml"
for cfg in "$base" "$simplified" "$variant"; do
  test -s "$cfg" || { echo "Missing config: $cfg" >&2; exit 1; }
done

source_nc="resources/bornholm_simplified_baseline/networks/composed_2035.nc"
test -s "$source_nc" || {
  echo "Missing simplified 2035 composed network: $source_nc" >&2
  exit 1
}

run_name="bornholm_simplified_barrier"
experiment_nc="resources/${run_name}/networks/composed_2035.nc"
solved_nc="results/${run_name}/networks/solved_2035.nc"
mkdir -p "$(dirname "$experiment_nc")"
if [[ -e "$experiment_nc" && ! -L "$experiment_nc" ]]; then
  echo "Refusing to overwrite non-symlink: $experiment_nc" >&2
  exit 1
fi
ln -sfn "$(realpath "$source_nc")" "$experiment_nc"

pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper" || { echo "Missing executable pixi wrapper: $pixi_wrapper" >&2; exit 1; }

export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1

echo "Solver experiment: $run_name"
echo "Input: $source_nc"
echo "Output: $solved_nc"
echo "Gurobi: Method=2, Threads=4, Crossover=0; no solver TimeLimit"
echo "Slurm: fatq, exclusive 256GB-class node, full memory, 24h walltime"
echo "Other network inputs and solver runs remain unchanged."

cmd=("$pixi_wrapper" run snakemake
  --cores 4
  --resources mem_mb=240000
  --rerun-incomplete
  --configfile "$base" "$simplified" "$variant"
  --allowed-rules solve_network)

if [[ "$mode" == "dry-run" ]]; then
  "${cmd[@]}" --dry-run -- "$solved_nc"
else
  "${cmd[@]}" -- "$solved_nc"
fi
