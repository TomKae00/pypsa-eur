#!/usr/bin/env bash
#SBATCH --job-name=bornholm-study-check
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=slurm-bornholm-study-check-%j.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
[[ "${SLURM_JOB_PARTITION:-}" == workq ]] || { echo 'Checks must run in workq.' >&2; exit 2; }
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1
"$pixi_wrapper" run python -m scripts.check_bornholm_networks \
    --overlay config/config.bornholm-3h.yaml --solved
