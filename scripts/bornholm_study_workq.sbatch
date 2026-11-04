#!/usr/bin/env bash
#SBATCH --job-name=bornholm-reviewed-build
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-reviewed-build-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
exec bash scripts/bornholm_calibrated_workq.sbatch "$@"
