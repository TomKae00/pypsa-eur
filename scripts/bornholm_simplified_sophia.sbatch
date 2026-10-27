#!/usr/bin/env bash
#SBATCH --job-name=bornholm-simple
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-simplified-%j.out

# Independent build/solve of the simplified full-year Bornholm baseline.
# No change to resources/results for bornholm_zonal_baseline or solver experiments.
# Use two dependent jobs: build must finish successfully before solve starts.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

stage="${1:-build}"
case "$stage" in
  build|solve|all|check|dry-run) ;;
  *) echo "Usage: sbatch scripts/bornholm_simplified_sophia.sbatch [build|solve|all|check|dry-run]" >&2; exit 2 ;;
esac

base="config/config.bornholm.yaml"
overlay="config/config.bornholm-simplified.yaml"
for f in "$base" "$overlay"; do test -s "$f"; done

pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1

name="bornholm_simplified_baseline"
composed="resources/${name}/networks/composed_2035.nc"
solved="results/${name}/networks/solved_2035.nc"
cores="${SLURM_CPUS_PER_TASK:-16}"
smk=("$pixi_wrapper" run snakemake --cores "$cores" --resources mem_mb=110000 \
     --rerun-incomplete --configfile "$base" "$overlay")

echo "Stage: $stage"
echo "Input composed network: $composed"
echo "Solved network: $solved"
echo "Simplifications: EU CO2 pooling; no CO2 pipelines, methanation, H2 fuel cell, solar thermal"

if [[ "$stage" == dry-run ]]; then
  "${smk[@]}" --dry-run -- "$solved"
  exit 0
fi

# Validate the existing composed network without repeating the successful build.
if [[ "$stage" == check ]]; then
  "$pixi_wrapper" run python scripts/check_bornholm_simplified.py --stage composed
  exit 0
fi

if [[ "$stage" == build || "$stage" == all ]]; then
  "${smk[@]}" -- "$composed"
  "$pixi_wrapper" run python scripts/check_bornholm_simplified.py --stage composed
fi

if [[ "$stage" == solve || "$stage" == all ]]; then
  test -s "$composed" || { echo "Simplified build has not completed: $composed" >&2; exit 1; }
  "$pixi_wrapper" run python scripts/check_bornholm_simplified.py --stage composed
  # Use only solve_network to avoid accidental regeneration of any upstream files.
  "${smk[@]}" --allowed-rules solve_network -- "$solved"
  "$pixi_wrapper" run python scripts/check_bornholm_simplified.py --stage solved
fi
