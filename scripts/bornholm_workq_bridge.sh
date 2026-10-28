#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Run on the host, beside containerised Snakemake. No host Python/Git required.
set -uo pipefail
spool="$1"
wrapper="$2"
shopt -s nullglob
while [[ ! -e "$spool/stop" ]]; do
    for request in "$spool"/*.request; do
        mapfile -d '' -t fields < "$request"
        token="${request%.request}"
        operation="${fields[0]:-}"
        code=0
        case "$operation" in
          submit)
            threads="${fields[1]:-}"
            jobscript="${fields[2]:-}"
            if [[ ! "$threads" =~ ^[1-4]$ || ! -f "$jobscript" ]]; then
                printf '%s\n' 'Invalid build submission' > "$token.err"
                code=2
                : > "$token.out"
            else
                printf -v wrapped '%q run bash %q' "$wrapper" "$jobscript"
                sbatch --parsable --partition=workq --nodes=1 --ntasks=1 \
                    --cpus-per-task="$threads" --exclusive --mem=0 --time=04:00:00 \
                    --job-name=bornholm-build --chdir="$PWD" \
                    --output="$spool/build-%j.out" --wrap="$wrapped" \
                    > "$token.out" 2> "$token.err" || code=$?
                if (( code == 0 )); then
                    id=$(head -n 1 "$token.out")
                    printf '%s\n' "${id%%;*}" >> "$spool/job_ids.txt"
                fi
            fi
            ;;
          status)
            sacct -n -X -j "${fields[1]}" --format=JobIDRaw,State,ExitCode -P \
                > "$token.out" 2> "$token.err" || code=$?
            ;;
          cancel)
            scancel "${fields[@]:1}" > "$token.out" 2> "$token.err" || code=$?
            ;;
          *)
            printf '%s\n' 'Unknown bridge operation' > "$token.err"
            : > "$token.out"
            code=2
            ;;
        esac
        rm -f "$request"
        printf '%s\n' "$code" > "$token.exit.tmp"
        mv "$token.exit.tmp" "$token.exit"
    done
    sleep 0.2
done
