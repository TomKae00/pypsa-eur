#!/usr/bin/env bash
#SBATCH --job-name=hindcast-reviewed-2023
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=slurm-hindcast-reviewed-2023-%j.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
stage="${1:-all}"
case "$stage" in all|build|solve|dry-run) ;; *) echo 'Usage: sbatch scripts/hindcast_reviewed_workq.sbatch [all|build|solve|dry-run]' >&2; exit 2 ;; esac
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
input="resources/hindcast-reviewed-2023/networks/composed_2023.nc"
output="results/hindcast-reviewed-2023/network.nc"
smk=("$pixi_wrapper" run snakemake --executor local --cores 16 --resources mem_mb=90000
    --configfile config/hindcast/build.yaml config/hindcast/build-2023.yaml)
if [[ "$stage" == dry-run ]]; then
    "${smk[@]}" --dry-run -- "$input"
    exit
fi
if [[ "$stage" == all || "$stage" == build ]]; then
    "${smk[@]}" -- "$input"
fi
if [[ "$stage" == all || "$stage" == solve ]]; then
    test -f "$input"
    [[ ! -e "$output" ]] || { echo "Accepted result already exists: $output" >&2; exit 2; }
    "$pixi_wrapper" run python -m scripts.hindcast solve "$input" --output "$output" --solver gurobi --threads 16
    "$pixi_wrapper" run python -m scripts.hindcast prepare-zonal "$output" --config config/hindcast/benchmark-reviewed-2023.yaml
    if [[ -f data/hindcast/prices_zonal_2023.csv ]]; then
        "$pixi_wrapper" run python -m scripts.hindcast benchmark --config config/hindcast/benchmark-reviewed-2023.yaml
    else
        echo 'Solve passed. Fetch observed prices before running the optional benchmark.'
    fi
fi
