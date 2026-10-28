#!/usr/bin/env bash
#SBATCH --job-name=bornholm-gas-rebuild
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --output=slurm-bornholm-gas-rebuild-%j.out

# Rebuild gas-dependent inputs and the selected composed baseline; no optimisation.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
mode="${1:-dry-run}"
variant="${2:-original}"
case "$mode" in
  dry-run|build) ;;
  *) echo "Usage: sbatch scripts/bornholm_gas_capacity_rebuild.sbatch [dry-run|build] [original|simplified]" >&2; exit 2 ;;
esac
case "$variant" in
  original) name="bornholm_zonal_baseline" ;;
  simplified) name="bornholm_simplified_baseline" ;;
  *) echo "Unknown baseline: $variant" >&2; exit 2 ;;
esac
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "Submit this script with sbatch from the project directory." >&2
  exit 2
fi
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
test -f config/config.bornholm.yaml
configs=(config/config.bornholm.yaml)
if [[ "$variant" == simplified ]]; then
  configs+=(config/config.bornholm-simplified.yaml)
fi
target="resources/${name}/networks/composed_2035.nc"

if [[ "$mode" == build ]]; then
  backup_dir="gas_capacity_backup_$(date +%Y%m%dT%H%M%S)_${SLURM_JOB_ID}_${variant}"
  mkdir -p "$backup_dir"
  git rev-parse HEAD > "$backup_dir/git_head.txt"
  git diff HEAD > "$backup_dir/code_changes.patch"
  for item in "$target" \
      "resources/${name}/gas_network.csv" \
      "resources/${name}/gas_network_clustered.csv" \
      "${configs[@]}"; do
    if [[ -f "$item" ]]; then
      cp --parents -p -- "$item" "$backup_dir/"
    fi
  done
  echo "Pre-rebuild inputs saved in $backup_dir"
fi

cmd=("$pixi_wrapper" run snakemake
  --cores "${SLURM_CPUS_PER_TASK:-4}"
  --resources mem_mb=110000
  --rerun-incomplete
  --configfile "${configs[@]}"
  --forcerun build_gas_network cluster_gas_network)
if [[ "$mode" == dry-run ]]; then
  "${cmd[@]}" --dry-run -- "$target"
else
  "${cmd[@]}" -- "$target"
fi
