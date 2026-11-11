#!/usr/bin/env bash
#SBATCH --job-name=bornholm-welfare-prepare
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=slurm-bornholm-welfare-prepare-%j.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
stage="${1:-all}"
case "$stage" in all|prepare|plan) ;; *) echo 'Usage: sbatch scripts/bornholm_welfare_workq.sbatch [all|prepare|plan]' >&2; exit 2 ;; esac
[[ "${SLURM_JOB_PARTITION:-}" == workq ]] || { echo 'Submit this controller in workq.' >&2; exit 2; }
pixi_wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
test -x "$pixi_wrapper"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_OPENBLAS_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1
if [[ "$stage" == plan ]]; then
    "$pixi_wrapper" run python -m scripts.run_bornholm_welfare --plan
    exit
fi
"$pixi_wrapper" run python -m scripts.run_bornholm_welfare --prepare
if [[ "$stage" == all ]]; then
    concurrency="${BORNHOLM_WELFARE_CONCURRENCY:-6}"
    [[ "$concurrency" =~ ^([1-9]|1[0-2])$ ]] || { echo 'Concurrency must be 1-12.' >&2; exit 2; }
    spool="results/bornholm_welfare_2035_v1/submissions/${SLURM_JOB_ID}"
    mkdir -p "$spool"
    array=$(sbatch --parsable --array="0-11%${concurrency}" scripts/bornholm_welfare_solve.sbatch)
    printf '%s\n' "$array" > "$spool/solve_array.txt"
    # Numerical flags retain candidates; collect diagnostics even if an array task fails its screens.
    if ! collector=$(sbatch --parsable --dependency="afterany:${array%%;*}" \
        --kill-on-invalid-dep=yes scripts/bornholm_welfare_collect.sbatch); then
        scancel "${array%%;*}"
        exit 1
    fi
    printf '%s\n' "$collector" > "$spool/collect_job.txt"
    echo "Welfare array: $array (12 cases; at most $concurrency concurrently)"
    echo "Collection job: $collector (after every task ends, including flagged candidates)"
fi
