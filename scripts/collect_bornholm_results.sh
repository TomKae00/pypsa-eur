#!/usr/bin/env bash
# Read-only collection of the eight completed 4921050 Bornholm solutions.
# Run from the cluster repository: bash collect_bornholm_results.sh
# Compression runs on workq. No rebuild, optimization, commit or push.
set -euo pipefail
prefix=bornholm_3h_s1_w2_o3
cases=(baseline radial hybrid radial_de optimised open baseline_expandable open_expandable)
review="${2:-bornholm_results_review_4921050}"
[[ "$review" != */* && "$review" != .* ]] || { echo 'Use a plain output directory name.' >&2; exit 2; }
[[ -f Snakefile && -f config/config.bornholm-3h.yaml ]] || { echo 'Run from /home/s253303/projects/pypsa-eur-bornholm.' >&2; exit 2; }
for case in "${cases[@]}"; do
    for path in "resources/$prefix/bornholm_zonal_$case/networks/composed_2035.nc" \
                "results/$prefix/bornholm_zonal_$case/networks/solved_2035.nc" \
                "results/$prefix/bornholm_zonal_$case/summary.json" \
                "results/$prefix/bornholm_zonal_$case/logs/solve_network/solver_2035.log"; do
        [[ -f "$path" ]] || { echo "Missing required file: $path" >&2; exit 2; }
    done
done
mkdir -p "$review/metadata"

if [[ "${1:-}" != --worker ]]; then
    # Git is available on the login node but can be absent on compute nodes.
    if command -v git >/dev/null 2>&1; then
        git rev-parse HEAD > "$review/metadata/git_head.txt"
        git status --short > "$review/metadata/git_status.txt"
        git diff HEAD -- config scripts rules profiles Snakefile pixi.toml pixi.lock \
            data/energy_islands data/transmission_projects data/custom_extra_functionality.py \
            > "$review/metadata/uncommitted_source_changes.patch"
    else
        printf '%s\n' 'Git unavailable; use recorded run hashes and archived source snapshot.' > "$review/metadata/git_unavailable.txt"
    fi
    pwd > "$review/metadata/project_path.txt"
    date -u '+%Y-%m-%dT%H:%M:%SZ' > "$review/metadata/collection_started_utc.txt"
    if [[ -z "${SLURM_JOB_ID:-}" ]]; then
        script=$(readlink -f "$0")
        srun --partition=workq --nodes=1 --ntasks=1 --cpus-per-task=1 \
            --time=01:00:00 bash "$script" --worker "$review"
        exit
    fi
fi

[[ "${SLURM_JOB_PARTITION:-}" == workq ]] || { echo 'Collection/compression must run on workq.' >&2; exit 2; }
printf '%s\n' "${SLURM_JOB_ID:-unrecorded}" > "$review/metadata/collection_job.txt"
hostname > "$review/metadata/collection_node.txt"
if command -v sacct >/dev/null 2>&1; then
    sacct -S 2026-10-28 -j 4921050,4921051 \
        --format=JobIDRaw,JobName%32,Partition,State,ExitCode,Elapsed,Timelimit,AllocCPUS,ReqMem,MaxRSS,NodeList -P \
        > "$review/metadata/job_accounting.txt" 2> "$review/metadata/job_accounting.stderr" || true
fi
wrapper="${PIXI_SOPHIA:-$HOME/bin/pixi-sophia}"
if [[ -x "$wrapper" ]]; then
    "$wrapper" run python - <<'PY' > "$review/metadata/current_environment.txt" 2>&1 || true
import importlib.metadata as metadata
import sys
print(sys.version)
for package in ('pypsa', 'linopy', 'gurobipy', 'highspy', 'pandas', 'numpy', 'xarray', 'snakemake'):
    try:
        print(package, metadata.version(package))
    except metadata.PackageNotFoundError:
        print(package, 'not installed')
PY
fi

for case in "${cases[@]}"; do
    result="results/$prefix/bornholm_zonal_$case"
    members=("resources/$prefix/bornholm_zonal_$case/networks/composed_2035.nc"
             "$result/networks/solved_2035.nc" "$result/summary.json" "$result/logs")
    [[ ! -d "$result/benchmarks" ]] || members+=("$result/benchmarks")
    printf '%s\n' "${members[@]}" > "$review/metadata/${case}_members.txt"
    echo "Packing $case ..."
    tar --dereference -cf - "${members[@]}" | gzip -1 > "$review/${case}.tar.gz.tmp"
    mv "$review/${case}.tar.gz.tmp" "$review/${case}.tar.gz"
done

manifest="$review/metadata/support_members.nul"
: > "$manifest"
: > "$review/metadata/missing_optional_files.txt"
add_support() {
    if [[ -e "$1" ]]; then
        printf '%s\0' "$1" >> "$manifest"
    else
        printf '%s\n' "$1" >> "$review/metadata/missing_optional_files.txt"
    fi
}
for path in config scripts rules profiles Snakefile pixi.toml pixi.lock envs/environment.yaml \
    data/energy_islands data/transmission_projects/tyndp2026 data/custom_extra_functionality.py \
    data/custom_costs.csv data/parameter_corrections.yaml data/unit_commitment.csv \
    resources/networks/base.nc resources/networks/base_extended.nc \
    resources/networks/simplified.nc resources/networks/clustered.nc \
    resources/transmission_projects "logs/$prefix" "results/$prefix/submissions"; do
    add_support "$path"
done
# Small shared tables include gas inputs, storage bounds and clustering maps.
while IFS= read -r -d '' path; do add_support "$path"; done \
    < <(find resources -maxdepth 1 -type f -name '*.csv' -print0)
for case in "${cases[@]}"; do
    while IFS= read -r -d '' path; do add_support "$path"; done \
        < <(find "resources/$prefix/bornholm_zonal_$case" -maxdepth 1 -type f -name '*.csv' -print0)
done
shopt -s nullglob
for path in doc/bornholm*.md logs/build_transmission_projects.log \
    logs/add_transmission_projects_and_dlr.log logs/cluster_network.log \
    logs/simplify_network.log logs/build_gas_network.log logs/cluster_gas_network.log \
    logs/base_network.log logs/time_aggregation_elec.log \
    slurm-bornholm-study-build-*.out slurm-bornholm-study-solve-4921050_*.out \
    slurm-bornholm-study-check-4921051.out; do
    add_support "$path"
done
add_support "$review/metadata"
echo 'Packing source, grid lineage, costs, provenance and accounting ...'
tar --dereference --exclude='__pycache__' --exclude='*.pyc' --null -T "$manifest" -cf - \
    | gzip -1 > "$review/support.tar.gz.tmp"
mv "$review/support.tar.gz.tmp" "$review/support.tar.gz"
(
    cd "$review"
    sha256sum ./*.tar.gz > SHA256SUMS
)
ls -lh "$review"/*.tar.gz
du -sh "$review"
echo "Ready: $PWD/$review"
echo 'Copy this directory to your laptop; upload the eight scenario archives and support.tar.gz for assessment.'
