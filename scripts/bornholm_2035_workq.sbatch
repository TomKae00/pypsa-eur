#!/usr/bin/env bash
#SBATCH --job-name=bornholm-2035-build
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-2035-build-%j.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
stage="${1:-all}"
case "$stage" in all|build|check|plan|dry-run) ;; *) echo 'Use all, build, check, plan or dry-run.' >&2; exit 2 ;; esac
[[ "${SLURM_JOB_PARTITION:-}" == workq ]] || { echo 'Submit this controller in workq.' >&2; exit 2; }
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export PIXI_SOPHIA="$pixi_wrapper"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
if [[ "$stage" == plan ]]; then
    "$pixi_wrapper" run python -m scripts.run_bornholm_2035 --plan
    exit
fi
if [[ "$stage" == check ]]; then
    "$pixi_wrapper" run python -m scripts.run_bornholm_2035 --prepare
    exit
fi
spool="$PWD/results/bornholm_2035_final_v1/submissions/${SLURM_JOB_ID:?}"
mkdir -p "$spool"
export BORNHOLM_WORKQ_BRIDGE="$spool"
export SINGULARITYENV_BORNHOLM_WORKQ_BRIDGE="$spool"
export SINGULARITYENV_PIXI_SOPHIA="$pixi_wrapper"
bridge_pid=""
workflow_pid=""
cleanup() {
    code=$?
    trap - EXIT
    if [[ -n "$workflow_pid" ]]; then kill "$workflow_pid" 2>/dev/null || true; fi
    touch "$spool/stop"
    if [[ -n "$bridge_pid" ]]; then wait "$bridge_pid" || true; fi
    if (( code != 0 )) && [[ -f "$spool/job_ids.txt" ]]; then
        mapfile -t ids < "$spool/job_ids.txt"
        if (( ${#ids[@]} )); then scancel "${ids[@]}" || true; fi
    fi
    exit "$code"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
"$pixi_wrapper" run python -m scripts.bornholm_provenance --output-dir "$spool"
targets_text=$("$pixi_wrapper" run python -m scripts.run_bornholm_2035 --list-build-targets)
mapfile -t targets <<< "$targets_text"
(( ${#targets[@]} == 15 )) || { echo 'Expected 15 composed inputs.' >&2; exit 2; }
smk=("$pixi_wrapper" run snakemake --profile profiles/bornholm-workq
    --configfile config/config.bornholm.yaml config/config.bornholm-3h.yaml
    config/config.bornholm-calibrated-base.yaml config/config.bornholm-2035.yaml)
if [[ "$stage" == dry-run ]]; then
    "${smk[@]}" --dry-run -- "${targets[@]}"
    exit
fi
bash scripts/bornholm_workq_bridge.sh "$spool" "$pixi_wrapper" &
bridge_pid=$!
"${smk[@]}" -- "${targets[@]}" &
workflow_pid=$!
wait "$workflow_pid"
workflow_pid=""
touch "$spool/stop"
wait "$bridge_pid"
bridge_pid=""
"$pixi_wrapper" run python -m scripts.run_bornholm_2035 --prepare
if [[ "$stage" == all ]]; then
    concurrency="${BORNHOLM_2035_CONCURRENCY:-6}"
    [[ "$concurrency" =~ ^([1-9]|1[0-5])$ ]] || { echo 'Concurrency must be 1-15.' >&2; exit 2; }
    array=$(sbatch --parsable --array="0-14%${concurrency}" scripts/bornholm_2035_solve.sbatch)
    printf '%s\n' "$array" > "$spool/solve_array.txt"
    if ! collector=$(sbatch --parsable --dependency="afterany:${array%%;*}" \
        --kill-on-invalid-dep=yes scripts/bornholm_2035_collect.sbatch); then
        scancel "${array%%;*}"
        exit 1
    fi
    printf '%s\n' "$collector" > "$spool/collect_job.txt"
    echo "2035 solve array: $array (15 cases, concurrency $concurrency)"
    echo "Collection after all tasks finish: $collector"
fi
echo "Build and input checks complete: $spool"
