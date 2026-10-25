#!/usr/bin/env bash
#SBATCH --job-name=bornholm-zonal
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-%j.out

# Sophia workq reports nominal memory (1 MB) to Slurm, so explicit --mem
# requests are rejected. Reserve the node exclusively for this large model.
# The Snakemake mem_mb budget below schedules work; it does not reserve RAM.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
test -f config/config.bornholm.yaml

# Uses the existing Sophia Singularity/Pixi wrapper and tested token-server licence.
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1

stage="${1:-build}"
case "$stage" in
    all|build|solve|dry-run) ;;
    *) echo "Usage: sbatch scripts/bornholm_sophia.sbatch [build|solve|all|dry-run] [short-scenario-name|all]" >&2; exit 2 ;;
esac
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "Run inside a Slurm allocation on a compute node, including dry-run; use srun or sbatch." >&2
    exit 2
fi

cores="${SLURM_CPUS_PER_TASK:-16}"
if (( cores < 16 )); then
    echo "The configured Gurobi preset requires at least 16 allocated CPUs." >&2
    exit 2
fi
smk=("$pixi_wrapper" run snakemake --cores "$cores" --resources mem_mb=110000 --rerun-incomplete
     --configfile config/config.bornholm.yaml)
composed=()
solved=()
available_text=$("$pixi_wrapper" run python -c 'import yaml; from pathlib import Path; print("\n".join(k.removeprefix("bornholm_zonal_") for k in yaml.safe_load(Path("config/scenarios_bornholm.yaml").read_text())))')
mapfile -t available <<< "$available_text"
selection="${2:-all}"
if [[ "$stage" == solve || "$stage" == all ]]; then
    selection="${2:-baseline}"
fi
selected=()
for scenario in "${available[@]}"; do
    if [[ "$selection" != all && "$selection" != "$scenario" ]]; then
        continue
    fi
    selected+=("$scenario")
    composed+=("resources/bornholm_zonal_${scenario}/networks/composed_2035.nc")
    solved+=("results/bornholm_zonal_${scenario}/networks/solved_2035.nc")
done
if (( ${#selected[@]} == 0 )); then
    echo "Unknown scenario '$selection'. Available: ${available[*]}" >&2
    exit 2
fi
if [[ "$selection" != baseline && "$selection" != all ]]; then
    composed=("resources/bornholm_zonal_baseline/networks/composed_2035.nc" "${composed[@]}")
fi
# Every comparison uses the same built baseline. Run the baseline pilot first.
if [[ ( "$stage" == solve || "$stage" == all ) && "$selection" != baseline && "$selection" != all && ! -f results/bornholm_zonal_baseline/networks/solved_2035.nc ]]; then
    echo "Solve the baseline pilot before other scenarios." >&2
    exit 2
fi
if [[ "$stage" == dry-run ]]; then
    "${smk[@]}" --dry-run -- "${solved[@]}"
    exit
fi
if [[ "$stage" == all || "$stage" == build ]]; then
    "${smk[@]}" -- "${composed[@]}"
fi
"$pixi_wrapper" run python scripts/check_bornholm_networks.py --scenarios "${selected[@]}"
if [[ "$stage" == all || "$stage" == solve ]]; then
    "${smk[@]}" -- "${solved[@]}"
    "$pixi_wrapper" run python scripts/check_bornholm_networks.py --solved --scenarios "${selected[@]}"
fi
