#!/usr/bin/env bash
#SBATCH --job-name=bornholm-pdhg-exp
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --time=06:00:00
#SBATCH --output=slurm-bornholm-pdhg-exp-%j.out

# Six-hour wall-clock budget per test. The Gurobi overlay imposes a shorter\n# optimizer time limit to leave room for Snakemake and result export.\n# Run a solver-only experiment against the already-built Bornholm baseline.
# The experiment gets its own results, solver logs and network input symlink.
# It NEVER recomposes or rewrites the baseline network.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

case "${1:-scale2}" in
  scale2)
    case_name="bornholm_solver_scale2"
    overlay="config/config.bornholm-solver-scale2.yaml"
    ;;
  relaxed)
    case_name="bornholm_solver_relaxed"
    overlay="config/config.bornholm-solver-relaxed.yaml"
    ;;
  *)
    echo "Usage: sbatch scripts/bornholm_solver_experiment.sbatch [scale2|relaxed]" >&2
    exit 2
    ;;
esac

pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
test -f "$overlay"

# Baseline build must have completed before this job is submitted.
source_nc="resources/bornholm_zonal_baseline/networks/composed_2035.nc"
test -s "$source_nc" || {
  echo "Missing completed baseline composed network: $source_nc" >&2
  exit 1
}

# Symlink only the immutable composed network; all outputs/logs are separate.
experiment_nc="resources/${case_name}/networks/composed_2035.nc"
mkdir -p "$(dirname "$experiment_nc")"
if [[ -e "$experiment_nc" && ! -L "$experiment_nc" ]]; then
  echo "Refusing to overwrite non-symlink $experiment_nc" >&2
  exit 1
fi
ln -sfn "$(realpath "$source_nc")" "$experiment_nc"

export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1

cores="${SLURM_CPUS_PER_TASK:-16}"
echo "Experiment: $case_name"
echo "Input: $source_nc"
echo "Output: results/${case_name}/networks/solved_2035.nc"
echo "Overlay: $overlay"

# Restrict to solve_network so Snakemake cannot rebuild or modify any
# upstream energy-system components; input is the completed baseline NetCDF.
"$pixi_wrapper" run snakemake \
  --cores "$cores" \
  --resources mem_mb=110000 \
  --rerun-incomplete \
  --configfile config/config.bornholm.yaml "$overlay" \
  --allowed-rules solve_network \
  -- "results/${case_name}/networks/solved_2035.nc"
