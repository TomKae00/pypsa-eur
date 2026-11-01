# Hourly fixed DK2–Germany Hybrid on Sophia fatq

This job builds and solves **only `bornholm_zonal_hybrid`**, with the original
2013 hourly chronology (8,760 one-hour snapshots) and 2035 planning assumptions.
It keeps the selected three-hour study's spatial/sector setup: 52 mainland
bidding zones, solar/onshore/offshore resource classes 1/2/3, fixed mainland
transmission (`vfixed`), fixed 3,000 MW project wind, 1,200 MW export to DK2 and
2,000 MW export to Germany. Native internal bus IDs remain `DKE1` and `DE00`;
DK2 and DE-LU are their presentation labels. Export capacities do not impose an
annual energy split.

## Submit on Sophia

```bash
cd /home/s253303/projects/pypsa-eur-bornholm
git pull --ff-only origin feature/bornholm-energy-island
sbatch scripts/bornholm_hourly_hybrid_fatq.sbatch
```

The job requests one exclusive fatq node, **32 allocated CPUs, all node RAM
(`--mem=0`) and 72 hours**. Gurobi uses **16 solver threads** by default, as in
the selected `nf1_aggfill_t16` campaign. Reserving 32 CPUs does not force 32
Gurobi threads. The existing `~/bin/pixi-sophia` wrapper and Gurobi licence are
used. Snakemake runs locally inside this allocation with at most 16 build CPUs;
it submits no nested Slurm jobs. Building the fresh hourly input happens first,
so there is no manual dependency or separate baseline solve to submit.

Configuration order is essential:

1. `config/config.bornholm.yaml`: physical study assumptions.
2. `config/config.bornholm-3h.yaml`: selected resource classes.
3. `config/config.bornholm-hourly-hybrid.yaml`: hourly chronology, fresh namespace,
   only fixed Hybrid, selected solver preset.

The third layer switches averaging, segmentation and representative periods
off and disables shared processed resources for this single case. In particular,
hourly snapshot weights cannot replace the three-hour campaign's common cache.
This may require additional preprocessing on the first submission.
Preflight rejects a three-hour source, non-unit physical time weights,
stale resource classes, wrong project capacities/endpoints or expandable
mainland lines. Annual budgets and costs are not prorated. The native model
preparation and supplementary constraints remain in use.

## Memory retry after job 4937219

The first hourly attempt completed input building in 18m56s. Gurobi then stopped
with `Status=17` (`MEM_LIMIT`), before any barrier iteration. Its soft cap was
175.903 GB. Gurobi recorded 205.768 GB peak allocated memory; the process's
`ru_maxrss` was 164.778 GiB and Slurm recorded 153.21G. These metrics use different
accounting and sampling; they must not be added together. The node reported
261.547 GB total physical RAM and 234.537 GB available after model construction.
The lower final `MemUsed=37.706 GB` is not the peak allocation.

An explicit retry allows a larger budget while preserving the original default:

```bash
cd /home/s253303/projects/pypsa-eur-bornholm
git pull --ff-only origin feature/bornholm-energy-island
sbatch scripts/bornholm_hourly_hybrid_fatq.sbatch 16 retry 4937219
```

This mode skips Snakemake and reuses the existing hourly composed network. The
runner verifies its SHA256, the pinned native physical code/configuration/custom
constraints and all solver options except memory/time budgets against the
previous job's summary. A changed input, changed thread count or changed physical
formulation aborts the retry. Outputs go into the new job's results directory.

After model construction, `retry` sets `SoftMemLimit` to the smallest of:

- **210 decimal GB**;
- **90% of the then-available physical RAM**;
- **available RAM minus a 20 GB reserve**.

At the first job's measured availability this selects 210 GB, 34.1 GB above the
previous configured limit and 4.2 GB above the observed Gurobi peak. It may get
past the previous stop, but the subsequent factorisation/iterations may require
more. The reserve is a budgeting rule, not a guaranteed physical free-memory
floor: Gurobi can overshoot its soft cap between checks and its accounting
differs from process RSS. More memory or a separately tested method may still
be needed. Keep 16 threads so this experiment changes the memory budget only.

The optional third launcher argument is the numeric reference job ID. For a
later retry, use the relevant prior job's ID instead. The original
`sbatch scripts/bornholm_hourly_hybrid_fatq.sbatch` command still uses the standard
190 GB / 75% policy.

## Solver and memory settings

The starting preset preserves the selected campaign's settings:

| Parameter | Setting |
|---|---:|
| Method | 2: barrier |
| Threads | 16 |
| Crossover | 0 |
| BarConvTol | 1e-8 |
| NumericFocus | 1 |
| AggFill / PreDual | -1 / 0 |
| ScaleFlag / BarHomogeneous | -1 / -1 |
| Seed | 123 |
| Objective constant | included |

This is the best-supported starting preset from the existing study, not a
guarantee of convergence or the fastest hourly configuration. The inherited
`GURO_PAR_BARDENSETHRESH=200` setting is retained from the tested campaign.
There is no automatic fallback, crossover activation, tolerance relaxation or
technology removal.

After native Linopy model construction, the runner measures Linux
`MemAvailable`. It caps Gurobi at the smaller of **190 decimal GB** and **75% of
the RAM still available**. This reserves space for allocations outside Gurobi's
accounting. It is a soft Gurobi limit, not a process-wide memory guarantee;
Python/Linopy construction can exhaust memory before that limit is installed.
The solver receives the remaining part of a 70-hour total build/model/solve
budget, leaving two hours under the default walltime for export and shutdown.

Sophia's published fatq specification is 256 GB per node. Earlier hourly
Bornholm jobs failed around 222 GiB RSS, and an earlier simplified hourly
barrier reported roughly 1 TB of factorisation memory. These are different
inputs from this 1/2/3-class Hybrid, so they are warnings about feasibility on
fatq, not predictions for this job. A fatq allocation is not a 1 TB allocation.
The selected 16-thread three-hour solves peaked around 54–59 GiB; a 32-thread
trial stopped at its soft memory limit around 96.65 GiB. These runs were not a
controlled demonstration that doubling barrier threads doubles RAM.

To run an explicit 32-thread experiment instead of the default:

```bash
sbatch scripts/bornholm_hourly_hybrid_fatq.sbatch 32
```

It uses the same node allocation and memory protection, with 32 Gurobi threads.
Increasing threads does not increase the node's available RAM and is not
guaranteed to speed up the factorisation. Start with 16; inspect the actual
factorisation estimate and measured peak memory before choosing 32.

## Outputs and monitoring

```text
resources/bornholm_1h_s1_w2_o3/bornholm_zonal_hybrid/networks/composed_2035.nc
results/bornholm_1h_s1_w2_o3/<job-id>/bornholm_zonal_hybrid/
```

Each job gets a new results directory and refuses to overwrite it. It records
`summary.json`, `effective_config.yaml`, `solver.log`, `python.log`, `memory.log`,
input/code hashes, native Gurobi quality, barrier progress and coal-store
screening. A finite assigned solution is first saved as `candidate_2035.nc`.
It becomes `solved_2035.nc` only after native optimal termination and the
existing diagnostic quality screens pass. A `SUBOPTIMAL` or failed-screen
candidate needs review. A memory/time stop without an assigned solution may
leave logs and a summary only. A hard Slurm kill may leave the last checkpoint
marked running; check `sacct` as well. Screening is not an independent
optimality certificate. No solved network is written into the three-hour
campaign's directory.

Replace `JOBID` below with the number printed by `sbatch`:

```bash
squeue -j JOBID
tail -f slurm-bornholm-hybrid-1h-JOBID.out
tail -F results/bornholm_1h_s1_w2_o3/JOBID/bornholm_zonal_hybrid/solver.log
sacct -j JOBID --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,ExitCode
```

The solver log appears after preprocessing and model construction. Check
`Barrier statistics` / `Factor NZ` before interpreting runtime. If the factor
estimate exceeds the available RAM, more threads will not make it fit; the next
step is a larger-memory allocation or a separately validated alternative solve.

Official references:

- [Sophia partitions and node memory](https://dtu-sophia.github.io/docs/scheduler/)
- [Gurobi continuous-model parameter guidance](https://docs.gurobi.com/projects/optimizer/en/current/concepts/parameters/guidelines.html)
- [Gurobi parameters, including SoftMemLimit](https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html)


## Additional AggFill=0 trial

To compare presolve aggregation against the automatic setting while retaining
16-thread barrier and the guarded 210 GB memory policy:

```bash
git pull --ff-only origin feature/bornholm-energy-island
sbatch scripts/bornholm_hourly_hybrid_fatq.sbatch 16 aggfill0 4937219
```

This trial reuses the same existing hourly input and verifies its provenance.
It permits only `AggFill=0` in addition to the memory/time budget differences;
all other solver settings and the physical formulation remain pinned to the
reference. It writes to a separate results directory under its new job ID and
can run alongside the unchanged `retry` trial. Lower fill may reduce memory,
but the presolved model can retain more rows/columns, so savings are not assured.
This addition was pushed without running tests, as requested.
