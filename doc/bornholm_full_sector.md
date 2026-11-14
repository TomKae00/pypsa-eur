# One additional 2035 full-sector Hybrid test

This is a single memory/feasibility test alongside the existing 15-case campaign.
It uses the same central 2035 Hybrid assumptions: 2013 weather, three-hour
chronology, 3 GW Bornholm wind, 1.2 GW DK2 and 2 GW German exports, EV share 45%,
managed charging 50%, V2G 25%, industrial transition progress 45%, background
offshore factor 0.85 and dedicated Bornholm factor 0.8855. Background electricity
transmission remains fixed. The processed inputs and outputs have a separate
prefix, `bornholm_2035_full_sector_v1`, and only `central_hybrid` is built/solved.

## What is restored

| Native setting | Main campaign | Extra run |
| --- | --- | --- |
| `sector.H2_network` | false | true |
| `sector.H2_retrofit` | false | true |
| `sector.shipping` | false | true |
| `sector.aviation` | false | true |
| `sector.co2_spatial` | true | true |
| `sector.co2_network` | true | true |

Local hydrogen production, storage and use, carbon capture, DAC, regional CO2
storage, gas networks, heating, industry, biomass, road transport and agriculture
are already present. This run does not indiscriminately enable unrelated optional
technologies or external hydrogen imports.

The model uses PyPSA-Eur's native hydrogen pipeline candidates and coupled gas/H2
retrofit constraints. Existing gas pipelines become bounded conversion candidates;
their optimized remaining gas capacity is coupled to retrofitted hydrogen capacity.
The extra preflight accepts these native bounds rather than requiring fixed gas
pipelines. Shipping inherits the native 2035 fuel shares: 35% methanol, 65% oil,
and zero direct hydrogen. Aviation adds native kerosene demand. Enabling these
sectors does not force all their fuel to be renewable.

The relative 2035 CO2 ceiling remains 25% of the native 1990 reference. That
reference includes the enabled sectors, so restoring shipping/aviation changes
the absolute emissions allowance as well as demand. The effective constraint is
recorded in `summary.json`; compare it before interpreting price/cost changes.
An unpaired full-sector Hybrid run cannot establish Bornholm's incremental
welfare benefit. A matching full-sector baseline can be considered after the
memory test succeeds.

## Submit without changing the running campaign

The active campaign fingerprints all scripts/configuration, including newly added
files. Do not pull this branch into its active checkout. From that checkout on
Sophia, create a separate worktree:

```bash
git fetch origin feature/bornholm-2035-full-sector-test
git worktree add -b feature/bornholm-2035-full-sector-test ../pypsa-eur-bornholm-full-sector FETCH_HEAD
cd ../pypsa-eur-bornholm-full-sector
sbatch scripts/bornholm_full_sector_fatq.sbatch
```

Native cached HTTP downloads reuse the user's download cache automatically.
The new worktree rebuilds processed inputs, so preprocessing may take time.
If necessary, existing raw datasets can be copied into this worktree before
submission; do not share processed `resources` or mutable configuration with the
active campaign.

The job reserves one whole fatq node, 32 CPUs and all node RAM for up to 72h.
Preprocessing uses four cores; the single native annual solve uses 16 Gurobi
threads, barrier, crossover off and NumericFocus 1. No rolling horizon, artificial
load shedding or load sinks are introduced. Build/model time counts against a
70h solve budget, leaving two hours for export and shutdown. After Linopy model
construction the Gurobi SoftMemLimit is at most 210 decimal GB and at most 75% of
then-available memory. Model construction itself can still exhaust RAM.

The shared annual solver retains the existing hourly runner's default behavior.
The additional run supplies its own input namespace and sector preflight.

## Inspect the outcome

```bash
squeue -j JOBID
sacct -j JOBID --units=G --format=JobID%24,State,Elapsed,MaxRSS,ExitCode
cat results/bornholm_2035_full_sector_v1/JOBID/central_hybrid/summary.json
```

The result directory records effective configuration, input/code hashes, model
dimensions, phase checkpoints, Python peak RSS, a 30-second memory log, solver
log/quality and termination status. Slurm `MaxRSS` on the job steps additionally
covers preprocessing and hard kills. If the kernel kills the process, the last
checkpoint remains; Python cannot write a final exception report after SIGKILL.
A finite assigned candidate is retained before quality screening. Only a passed
optimal solve becomes `solved_2035.nc`; failed/review candidates are not labelled
as accepted results. This preparation does not itself submit a Sophia job.
