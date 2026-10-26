# Parallel hindcast jobs on Sophia

Run these commands yourself in the Sophia SSH terminal. This helper does not
connect to Sophia, start a remote agent or install an IDE backend. Submit from
the separate `feature/historical-hindcast` checkout.

## Rebuild and solve the corrected 2023 case

With no hindcast jobs active, the human-invoked command below rebuilds the
corrected network and immediately solves it in one allocation. No intermediate
manual review or separate geography-only solve is required for this run.

```bash
cd ~/projects/pypsa-eur-hindcast
git pull --ff-only
mkdir -p logs
sbatch --job-name=hindcast-corrected-2023 --cpus-per-task=8 \
  --output=logs/hindcast-corrected-2023-%j.log \
  scripts/hindcast_sophia.sbatch rebuild-solve 2023 "$(git rev-parse HEAD)"
squeue -r -u "$USER"
```

The native `--forcerun base_network compose_network` rebuilds spatial
assignments and dependent resources, including the composed network. The
completed ERA5 cutout and unchanged raw-data cache are reused. Generated
resources in `resources/hindcast-build-2023` are updated; the original solved
network and audit in `results/hindcast-build-2023` remain intact. The corrected
solve, audit and submitted revision go to `results/hindcast-corrected-2023`.
Flow-query products use a separate `flow-sources-corrected-2023` filename.
An existing corrected result/audit prevents a new submission from rebuilding or
overwriting it. Do not submit this mode twice or edit this checkout while active.

This uses the successful electricity-only 2023 case's `workq`, Gurobi options
and native dynamic costs. It does not submit the multi-year campaign. Shell
syntax and local dispatch checks cover the launcher; the new HPC run itself
has not been executed here.

## Campaign layout

| Work | CPUs per job | Dependency |
| --- | --- | --- |
| Solve the existing 2023 composed network | 8 | None |
| Build/reuse 2022 weather | 4 | None |
| Build/reuse 2024 weather | 4 | None |
| Build/reuse 2025 weather | 4 | None |
| Build, prepare mapping and solve 2022 | 8 | Successful 2022 weather only |
| Build, prepare mapping and solve 2024 | 8 | Successful 2024 weather only |

Each job requests one node and one task in `workq`, with a 24-hour limit.
Slurm chooses placement; jobs can share a physical node when resources permit.
Multiple independent years provide parallelism. A standard Gurobi solve uses
threads within one node; requesting extra nodes does not distribute that solve.
Eight solve threads is an initial resource choice, not a measured optimum or
a promise of linear acceleration. CDS server queues and throughput may limit
concurrent weather downloads.

Weather and build/solve stages are separate allocations. Slurm `aftercorr`
releases each annual build immediately after its matching weather array task
succeeds. The 2022/2024 builds do not wait for each other, for 2023, or for 2025.
Failed weather tasks leave their downstream builds dependent; inspect/cancel
those pending jobs before submitting replacements.

The native workflow is retained. Annual `run.name` values separate network,
profile, demand and other generated-resource paths. Common completed downloads
from the 2023 build are reused. Snakemake permits concurrent instances with
disjoint output sets; its locks remain **enabled**. If it identifies an
unexpected shared writer, inspect the reported collision and rerun the affected
year after that writer completes. Do not use `--nolock` or `--unlock` to bypass
an active job. No duplicate checkouts, copied cutouts or separate environments
are needed.

## Requirements

- Existing `~/bin/pixi-sophia` container wrapper and installed project Pixi
  environment. The helper uses `pixi run --as-is` so concurrent jobs do not
  install packages or update the lock file.
- Existing CDS credentials in `~/.cdsapirc`.
- Gurobi licence at `~/gurobi.lic`, or set host `GRB_LICENSE_FILE` before
  submission. The worker forwards its path into Singularity.
- The completed 2023 network and reviewed 2023 bus mapping.
- Prepared fuel and CO2 files for 2022 and 2024 under
  `data/hindcast/paper-procedure/`. Missing inputs stop the relevant annual job.

No explicit Slurm memory request is used: this Sophia `workq` configuration
advertises placeholder memory and rejected `--mem=32G`. Every job uses its own
on-disk project `tmp/hindcast-JOBID-YEAR` directory, avoiding RAM-backed `/tmp`
for weather/solver temporary files. Temporary directories are retained on
failure for diagnosis; completed weather preparation normally removes its GRIB
working files.

## Submit

On Sophia:

```bash
cd ~/projects/pypsa-eur-hindcast
git pull --ff-only
bash scripts/submit_hindcast_sophia.sh
squeue -u "$USER"
```

The submission helper uses ordinary host Git/Slurm commands on the login node.
Python, Snakemake and Gurobi execute inside Slurm allocations. It submits:

1. The independent 2023 solve.
2. Weather array indexes `0=2022`, `1=2024`, `2=2025`.
3. Build/solve array indexes `0=2022`, `1=2024`, dependent on the
   corresponding successful weather tasks.

The printed IDs are saved in `logs/hindcast-submissions-TIMESTAMP-PID.txt`.
Do not submit the campaign twice while its jobs remain active. Keep the checkout
and environment unchanged while the campaign is queued/running. The login-node
submission helper checks the branch and tracked working-tree edits and passes
the submitted revision to each worker for recording. Compute nodes do not
provide host Git, so workers do not invoke it or independently verify that the
checkout still matches that revision. To change code, finish/cancel affected
jobs and submit replacements at the new revision. Queued or mid-run source
edits cannot be prevented by this helper.

The annual worker creates a new mapping template only when none exists, then
runs `prepare-zonal` again to enforce complete assignments and native shared
market connections. Current exact aliases include the three model-only islands.
Unknown bus IDs still stop the job for review rather than receiving a guessed
country assignment. The 2023 mapping is reused.

The 2023 `solve` stage does **not** call Snakemake or rebuild the composed file.
Native dynamic fuel/CO2 costs are used directly: no `--marginal-costs` argument.

## Watch progress and restart individual jobs

Replace job IDs with those printed at submission:

```bash
tail -f logs/hindcast-solve-2023-JOBID.log
sacct -j JOBID --format=JobID,JobName,State,ExitCode,Elapsed
```

Weather log names are `logs/hindcast-weather-ARRAYID_INDEX.log`.
Annual build/solve logs are
`logs/hindcast-build-solve-ARRAYID_INDEX.log`.

If a year fails, inspect its log and fix the reported cause. Once any old job
for that year has stopped, submit only the needed stage, for example:

```bash
sbatch --job-name=hindcast-build-solve-2022 \
  --output=logs/hindcast-build-solve-2022-%j.log \
  scripts/hindcast_sophia.sbatch build-solve 2022 "$(git rev-parse HEAD)"
```

This reuses completed native inputs and builds any still-missing dependencies.
Use `solve 2023` instead of `build-solve 2022` to restart only the 2023 solve.
Existing solve outputs or audits are not overwritten automatically; inspect any
partial output before choosing a new destination or removing it yourself.

## Restart after suboptimal barrier termination

The first 2023 solve on Sophia reached Gurobi status 13 (suboptimal) after
barrier progress stalled. The helper refused to export a solved network.
This is not an accepted price benchmark result.

The Gurobi solve now uses `Method=2`, `Crossover=-1` (automatic crossover
strategy) and `NumericFocus=2`, retaining the requested thread count.
Crossover permits simplex cleanup of the barrier point; NumericFocus requests
more care in numerical computations and can increase runtime. No feasibility,
optimality or barrier convergence tolerances are relaxed. Only an
`ok/optimal` termination is exported. Success of this retry is not established
by the configuration change.

The helper also explicitly sets PyPSA's `include_objective_constant=False`.
It already rejects extendable assets, so the fixed-capacity dispatch model has
no investment-offset term to optimise. This avoids the changing-default warning;
it is not evidence that an objective constant caused the observed failure.
The result metadata records the solver options and termination status.

If weather jobs are still active, **do not resubmit the whole campaign**.
Cancel the old downstream build/solve array before updating source, then attach
the replacement jobs to the existing weather array:

```bash
# Example IDs from the first campaign; use your actual array IDs.
scancel 4907353
git pull --ff-only
bash scripts/submit_hindcast_sophia.sh --weather-array 4907352
squeue -r -u "$USER"
```

This retries the 2023 solve and submits replacement 2022/2024 build/solve jobs
using the existing weather-array dependencies. It does not create new weather
jobs. The old 2023 solve must have stopped and no successful result may already
exist. Do not cancel an annual job that has already completed successfully;
submit just the missing individual stages in that case.

Reusing a weather array across this update is appropriate because only
`hindcast.py`, the submission helper and this guide change. Weather rules,
`build_cutout.py`, shared helpers, build configurations and the installed
environment are unchanged. Do not apply this restart pattern to arbitrary
source changes during active builds. New solve/build workers record the new
submitted revision, while the already-running weather workers retain theirs.
The weather array must still be known to Slurm for `aftercorr` to resolve.
If Slurm has purged a completed array, submit a new campaign: completed cutouts
are reused by the native workflow rather than downloaded from scratch.

Sources: [Gurobi numerical guidance](https://docs.gurobi.com/projects/optimizer/en/current/concepts/numericguide/numeric_parameters.html),
[Gurobi crossover parameters](https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#parametercrossover).
No tests or cluster solves were run to establish the retry outcome.

## 2025 is weather-only pending annual inputs

The configured `data/nuclear_p_max_pu.csv` has year columns only through 2024.
The native attachment code deliberately rejects a missing selected year.
Neither 2024 nuclear availability nor 2024 demand/capacities is silently
substituted. Before submitting a 2025 network build/solve, establish the intended
2025 nuclear input and confirm annual demand and renewable-capacity coverage.
The worker explicitly blocks a 2025 `build-solve` until this guard is updated
with that resolved input decision. A prepared 2025 weather file remains reusable.
This campaign therefore does not claim to produce a complete 2025 result yet.

## Bring completed results back to the laptop

No other year needs to finish before copying a successful 2023 result.
After Slurm reports the 2023 job as `COMPLETED` with exit `0:0`, run on the
**laptop**, not Sophia:

```bash
cd ~/PycharmProjects/pypsa-eur
git pull --ff-only
mkdir -p results/hindcast-build-2023 data/hindcast
rsync -avP \
  s253303@sophia.dtu.dk:projects/pypsa-eur-hindcast/results/hindcast-build-2023/ \
  results/hindcast-build-2023/
rsync -avP \
  s253303@sophia.dtu.dk:projects/pypsa-eur-hindcast/data/hindcast/bus_zones_2023.csv \
  data/hindcast/
```

The result directory contains `network.nc`, `network.audit.json` and
`submitted-revision.txt`. Download completed years separately by replacing
2023 in the paths. Do not transfer ERA5 cutouts, all intermediates or the Pixi
environment for results inspection.

Historical price retrieval and `hindcast.py benchmark --year 2023` can then run
locally using the existing documented price workflow. The other HPC jobs
continue independently. Numerical reconstruction quality, observed-price
coverage and model feasibility/load shedding still need assessment; successful
job completion alone is not validation.

No tests, dry runs, model builds, solves or Sophia connections were performed
while preparing these scripts.
