#!/usr/bin/env bash
#SBATCH --job-name=bornholm-study-build
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm-bornholm-study-build-%j.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
stage="${1:-all}"
case "$stage" in all|build|dry-run) ;; *) echo 'Usage: sbatch scripts/bornholm_study_workq.sbatch [all|build|dry-run]' >&2; exit 2 ;; esac
[[ -n "${SLURM_JOB_ID:-}" && "${SLURM_JOB_PARTITION:-}" == workq ]] || { echo 'This controller must run in workq.' >&2; exit 2; }
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
for tool in sbatch sacct scancel; do command -v "$tool" >/dev/null; done
export PIXI_SOPHIA="$pixi_wrapper"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1
export SINGULARITYENV_GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
spool="$PWD/results/bornholm_3h_s1_w2_o3/submissions/${SLURM_JOB_ID}"
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
targets_text=$("$pixi_wrapper" run python -m scripts.run_bornholm_study --list-build-targets)
mapfile -t targets <<< "$targets_text"
(( ${#targets[@]} == 8 )) || { echo 'Expected exactly eight build targets.' >&2; exit 2; }
if [[ "$stage" != dry-run ]]; then
    bash scripts/bornholm_workq_bridge.sh "$spool" "$pixi_wrapper" &
    bridge_pid=$!
fi
smk=("$pixi_wrapper" run snakemake --profile profiles/bornholm-workq
    --configfile config/config.bornholm.yaml config/config.bornholm-3h.yaml)
if [[ "$stage" == dry-run ]]; then
    "${smk[@]}" --dry-run -- "${targets[@]}"
    exit
fi
"${smk[@]}" -- "${targets[@]}" &
workflow_pid=$!
wait "$workflow_pid"
workflow_pid=""
"$pixi_wrapper" run python -m scripts.check_bornholm_networks --overlay config/config.bornholm-3h.yaml
"$pixi_wrapper" run python -m scripts.run_bornholm_study --preflight-all
touch "$spool/stop"
wait "$bridge_pid"
bridge_pid=""
if [[ "$stage" == all ]]; then
    # Eight independent annual solves; no baseline-solve dependency.
    solve_array=$(sbatch --parsable --partition=workq --array=0-7 scripts/bornholm_study_solve.sbatch)
    printf '%s\n' "$solve_array" > "$spool/solve_array.txt"
    echo "Eight workq solves submitted: array $solve_array (tasks 0 through 7)."
    check_job=$(sbatch --parsable --partition=workq --dependency="afterok:${solve_array%%;*}" \
        --kill-on-invalid-dep=yes scripts/bornholm_study_check.sbatch)
    printf '%s\n' "$check_job" > "$spool/check_job.txt"
    echo "Comparison check after all eight solves: $check_job"
fi
echo "Build checks complete. Submission records and worker logs: $spool"
