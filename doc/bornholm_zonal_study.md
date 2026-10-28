# Sector-coupled Bornholm comparison

## Current submission: all eight scenarios at three-hour resolution

Use `config/config.bornholm-3h.yaml` after the main study configuration. It sets
solar and tracking solar to one resource class, onshore wind to two, and offshore
AC/DC/floating wind to three. Floating offshore wind previously inherited one
class; it now also follows the requested three-class offshore setting. Dedicated
Bornholm project wind remains one fixed site. All 52 bidding zones, the full
2013 weather year, 2035 demand/cost assumptions and the original sector
formulation are retained. This is the eight-case study, not eight copies of the
separate simplified-baseline experiment. Spatial CO2, the CO2 network and capture
technologies retain the original configuration.

Native chronological averaging produces 2,920 snapshots, each with three-hour
objective, generator and storage weights. Annual costs, emissions budgets,
feedstock quantities and physical storage capacities are not prorated. Mean
renewable availability and demand smooth within-block variability, so these
results have a different temporal approximation from the hourly study.

```bash
cd /home/s253303/projects/pypsa-eur-bornholm
git pull --ff-only origin feature/bornholm-energy-island
sbatch scripts/bornholm_study_workq.sbatch
```

This single submission performs the following:

1. A workq controller runs Snakemake with `--executor cluster-generic` via the
   profile `profiles/bornholm-workq`. One DAG covers all eight composed networks.
   Up to eight rule jobs can run concurrently on separate exclusive workq nodes,
   using up to four CPUs each. Dependencies determine which jobs are ready.
2. Native `run.shared_resources.policy: base` shares eligible preprocessing
   outputs, including weather profiles, clustered geography and the time weights.
   Scenario-specific composition and horizon-specific resources retain their
   own paths. Processed costs and one unrelated postprocessing output are
   excluded from sharing to keep input/output/log wildcard handling consistent.
3. All composed networks pass comparison checks, including time weights, common
   gas/storage inputs and scenario definitions. Source preflight also checks
   undefined buses, positive fixed gas capacities, annual gas supply, resource
   classes and the study namespace. No solver is submitted if these checks fail.
4. A Slurm array `0-7` launches exactly eight independent native annual solves,
   each on an exclusive workq node, with four-thread Gurobi barrier, no crossover,
   `BarConvTol=1e-5`, 24-hour job walltime and a 20-hour solver limit. They have no
   dependency on a solved baseline. Availability of workq nodes and Gurobi licences
   determines how many actually start together.
5. A workq comparison job runs after all eight solves succeed. Slurm cancels this
   dependent check if its success dependency becomes impossible.

The controller, worker submission bridge, solve array and final check all specify
workq. The controller and solve/check scripts reject execution in any other
partition. Host Slurm commands are used through a filesystem bridge, while each
worker re-enters the existing `~/bin/pixi-sophia` container wrapper. Neither host
Git/Python nor container Slurm executables are required. The cluster-generic
plugin is already declared in the Pixi environment. All nodes must see the same
project filesystem and wrapper, as in the existing Sophia jobs.

The native solve entry point and its configured solve mode remain in use,
including supplementary constraints and the expandable-grid scenarios. Solver
memory is capped at the smaller of 70 decimal GB and 75% of available RAM before
native solving. This is Gurobi-accounted memory, not a whole-process cap; OS kills
and memory growth beyond graceful check points remain possible. The reduction
in variable count does not establish that factorisation fits workq. Three-hour
barrier runs can still stop for memory, time or infeasibility. No automatic
technology removal, solver fallback or change of partition is performed.

The overlay uses a separate namespace, leaving prior hourly composed/solved
networks and diagnostic results available:

```text
resources/bornholm_3h_s1_w2_o3/bornholm_zonal_<case>/networks/composed_2035.nc
results/bornholm_3h_s1_w2_o3/bornholm_zonal_<case>/networks/solved_2035.nc
results/bornholm_3h_s1_w2_o3/bornholm_zonal_<case>/summary.json
results/bornholm_3h_s1_w2_o3/bornholm_zonal_<case>/logs/solve_network/solver_2035.log
```

Array mapping: `0 baseline`, `1 radial`, `2 hybrid`, `3 radial_de`, `4 optimised`,
`5 open`, `6 baseline_expandable`, `7 open_expandable`. Solves refuse to overwrite
an existing solved network. A non-optimal output from the native entry point is
renamed `failed_2035.nc`, not accepted as a canonical solved result.

```bash
# Optional alternatives: build only, or inspect the build DAG without workers.
sbatch scripts/bornholm_study_workq.sbatch build
sbatch scripts/bornholm_study_workq.sbatch dry-run

# Live state, including the partition for every job.
watch -n 5 'squeue -u "$USER" -o "%.18i %.25j %.10P %.12T %.12M %.25R"'

# Follow baseline solver iterations after its solve starts.
tail -n 30 -F results/bornholm_3h_s1_w2_o3/bornholm_zonal_baseline/logs/solve_network/solver_2035.log
```

The controller log is `slurm-bornholm-study-build-<id>.out`; worker logs and job
IDs are under `results/bornholm_3h_s1_w2_o3/submissions/<controller-id>/`. That
directory also records the solve array and final comparison job IDs. Solve task
logs are `slurm-bornholm-study-solve-<array-id>_<index>.out`.

The sections below document the original hourly configuration and earlier
diagnostic submissions; use the new command above for this three-hour campaign.

Development validation: 69 targeted tests passed (two unrelated weather-profile
integration tests excluded). A native eight-scenario DAG dry-run expanded 216
build jobs from an empty resource cache, with eight distinct composition jobs,
one clustered network, one project profile and one time-aggregation scheme.
Remote HTTP metadata checks were skipped for this development dry-run. The real
cluster-generic executor also completed a small dependency graph against fake
Slurm commands: one shared task ran once, followed by eight dependent tasks;
every submission specified workq. Solver handoff, memory cap, failed-output
rejection and build/check failure propagation were tested. These tests do not
establish annual Gurobi memory use or full-dataset physical feasibility. No
Sophia jobs were submitted from development.

`config/config.bornholm.yaml` now defines the full study. It replaces the former
90-cluster, 3-hour electricity-only configuration. The separate smoke configuration
is retained for debugging; it is not the analysis run.

## Study definition

| Setting | Full comparison |
|---|---|
| Spatial footprint | Native 34-country European footprint, aggregated to bidding-zone polygons |
| Shared electricity markets | DE00 + LUG1; IE00 + GBNI, with native lossless `copper` links |
| Time resolution | Native default year: 8,760 hourly snapshots, 1 January–31 December 2013; no temporal reduction or rolling horizon |
| Planning and technology-cost year | 2035, overnight capacity expansion |
| Weather input | Native archived `europe-2013-sarah3-era5` cutout (SARAH-3/ERA5) |
| Electricity-load profile year | Native `load.fixed_year: false`, using the snapshot year (2013) |
| Sector demand base year | Native 2023 energy totals; native 2035 technology and end-use assumptions |
| Sectors | Electricity, heating, transport, industry, biomass, agriculture, hydrogen and gas infrastructure; shipping and aviation disabled in the current performance configuration |
| CO2 constraint | Native relative budget: 25% of the configured 1990 reference emissions, covering the enabled sectors |
| Main electricity grid | Native OSM 0.7 plus the covered TYNDP 2026 project subset for 2035; fixed background (`vfixed`) |
| Grid benchmarks | The same installed planned assets, with native further reinforcement (`vopt`) |
| Gas / hydrogen transmission | Existing gas pipelines fixed; native gas connectivity candidates enabled; hydrogen pipelines/retrofits disabled, in all eight cases |
| Losses | Native sector DC and distribution losses; AC loss approximation disabled |
| Renewable resource classes | Three for onshore wind, offshore AC/DC wind and solar |
| Project counterfactual | Project site reserved from generic offshore development in every case |

The comparison is entirely within the Bornholm study. It is independent of the
historical price/flow benchmark and does not require that branch's built network.
The weather and electricity-load settings match this checkout's
`config/config.default.yaml`; the custom 2023 ERA5 cutout has been removed.
Weather year **2013**, statistical energy-balance base year **2023**, and planning
year **2035** have different roles. Historical hourly weather supplies a coherent
wind/solar/temperature year for the future system; it does not turn the study into
a 2013 system reconstruction. Sector service demands use the native statistical
data and 2035 assumptions consistently across scenarios. This is a future scenario,
not a calibrated price forecast.

The native bidding-zone implementation aggregates the electrical network by zone
and retains AC network constraints between zones. It does not reproduce an NTC or
flow-based market-clearing model. DE/LU and Ireland/Northern Ireland keep separate
country buses for sector data while sharing an electricity price through their
unconstrained links. Sector buses increase the **total** bus count. Radial and
hybrid also add `BEI_HUB`; the onshore electricity bus IDs must be identical across
the Bornholm scenarios. The actual count is printed from the built networks,
not inferred from `n_clusters` (which is not a target in administrative mode).

### Common pipeline assumption

The Bornholm configuration retains the native natural-gas connectivity default:

```yaml
sector:
  gas_network: true
  gas_network_connectivity_upgrade: 1
  H2_network: false
  H2_retrofit: false
```

`gas_network: true` retains regional gas balances, existing pipelines and native
supply entry points. With `H2_retrofit: false`, existing pipeline capacities stay
fixed by native PyPSA-Eur. Connectivity level 1 generates additional investment
candidates to connect disconnected parts of the gas graph. Their capacities are
endogenous and carry the native investment costs. The candidate graph only counts
existing connections above 1,500 MW; a weaker existing corridor may therefore be
offered an additional candidate. This does not enable general reinforcement of
already connected, larger corridors.

The former Bornholm override of zero and its custom NetworkX guard have been
reverted. All eight cases share existing pipeline inputs and the same candidate
set and investment costs. Optimal new-gas investments may differ across cases.
Previously composed networks recording connectivity zero are stale and must be
rebuilt before solving or comparing cases.

`H2_network: false` still removes interregional hydrogen pipeline candidates;
regional hydrogen production, consumption and storage remain available.
`H2_retrofit: false` is already the upstream default. The hydrogen-network
restriction is a separate study assumption and is not required to keep existing
methane pipelines fixed.

The gas-data preprocessing correction preserves accepted source capacities and
lengths when applying selective replacements. Previously, accepted values became
NaN and clustering could sum missing capacities to zero. In the configured
SciGRID 1.1.2 dataset, the HU-RS source segment retains about 5,465 MW after this
correction. Restoring connectivity candidates does not replace correcting the
existing pipeline inputs. Rebuild gas-dependent resources before another solve.

### Shared geography and reservoir corrections

The Bornholm branch includes the generic preprocessing corrections from
hindcast commits `a477dca8249f2614ffddaf9aeaf6e65d35875bf1` and
`0b30b179d332232e380aa17c330808b4fb6456d1`:

- At administrative level `bz`, base bus regions use the actual native bidding
  zone polygons. Dissolving entire NUTS3 regions with majority-overlap labels
  previously moved boundaries, plants and transmission endpoints between zones.
  Corrected regions feed the existing country-constrained bus and plant mapping
  and the existing line aggregation; no plant-specific reassignment or manual
  interzonal transfer-capacity calibration is added.
- Missing reservoir durations are completed at individual plants before native
  capacity-weighted aggregation. Reported durations are retained and only the
  residual national storage estimate is allocated to plants with missing data.
  The erroneous overwrite of the first reservoir entry of each country is removed.
  Native inflows, cyclic storage and discharge constraints remain in use.

The assignment methods match the corrected hindcast branch. Numerical capacities
need not match a 2023 hindcast: Bornholm uses the planned 2035 grid, 2035 plant
eligibility and future investment assumptions. Its weather remains the native
2013 input. Cross-zone AC ratings retain native line-rating aggregation and
`s_max_pu`; they are not observed NTC or flow-based market limits.

The networks built before these corrections are stale. After pulling, run the
full `build` stage before any `solve`, so the changed base-network script and its
explicit bidding-zone input rebuild dependent regions, plants and networks.
Existing retrieved datasets and unchanged weather cutouts can be reused. The
comparison checker rejects the old pipeline settings and checks that every case
retains the same existing gas network and native gas candidate inputs, with no hydrogen pipeline candidates.

### Scenario comparisons

All run names have the prefix `bornholm_zonal_`:

| Suffix | Bornholm cables | Background grid |
|---|---|---|
| `baseline` | No project wind, hub or cables | Fixed planned grid |
| `radial` | DK2: fixed 3 GW | Fixed planned grid |
| `radial_de` | Germany: fixed 3 GW | Fixed planned grid |
| `hybrid` | DK2: fixed 1.2 GW; Germany: fixed 2 GW | Fixed planned grid |
| `optimised` | DK2 and Germany: independently 0–6 GW each | Fixed planned grid |
| `open` | DK2, Germany, SE4 and Poland: independently 0–6 GW each | Fixed planned grid |
| `baseline_expandable` | No project wind, hub or cables | Further reinforcement allowed |
| `open_expandable` | Same candidates as `open` | Further reinforcement allowed |

Every project case fixes wind investment at 3 GW. In `optimised` the **total and
allocation** are endogenous: there is no fixed-total 3.2 GW constraint. `open`
may choose zero Danish capacity or zero capacity to any other destination, so
separate forced DE–PL, DE–SE–PL and "no Denmark" cases are not duplicated here.
The 6 GW per-corridor upper bound is a research search limit, not an engineering
rating or a policy target. The post-solve checker warns if it binds; enlarge the
bound and rerun before claiming an unconstrained optimum. No minimum connection
or fixed construction charge is imposed on optional cables; native continuous
capacity investment is used.

Receiving buses are `DKE1`, `DE00`, `SE04` and `PL00`. The added SE4 and Polish
landing coordinates (14.13 E, 55.43 N and 15.57 E, 54.18 N) are illustrative coastal
route proxies near Ystad and Kołobrzeg, not approved landing sites. Their connection
costs use the same hub-to-landing distance × 1.25 convention as the existing routes.
Route feasibility, onshore landing reinforcement and fixed station costs would
need a separate engineering sensitivity before recommending construction.

All project cases use the same native, site-specific wind profile. Exports use
the native bidirectional DC-loss formulation: two non-negative directional flows,
synchronised capacity variables, one investment charge per physical cable.
Efficiency is `0.98 * 0.977 ** (length_km / 1000)` in either direction. Cable
lengths remain hub-to-landing distance times 1.25 estimates. See
[resource and cost accounting](bornholm_native_resources.md) for the site and
engineering approximations. The offshore hub receives no heat/hydrogen demand.

The common demand assumption is the native sector-coupled 2035 pathway, with
2023 energy totals, the default 2013 electricity-load profile and native sector shares/efficiency
assumptions across cases. It is not a blanket percentage increase in electricity
load. The checks compare all exogenous load series and static loads. Electricity
consumption by heat pumps, electrolysis and other endogenous conversion can
respond to the project. This supports a statement about optimal electrification
**by 2035**, not the timing of acceleration along a transition trajectory.

Compare each project case against `baseline` for its net system effect; compare
`hybrid` against both radial cases for the incremental hybrid effect; compare
`optimised` against `hybrid` for sizing/allocation and `open` against `optimised`
for destination choice. Compare `open_expandable - baseline_expandable` with
`open - baseline` for the effect of background grid flexibility. A cheaper system
objective is not itself country welfare or a consumer/producer surplus estimate.
The 70% German / 30% Danish financing arrangement motivated by expected power
flows is not a cable-capacity ratio and is not imposed on hourly flows. CfD
payments, country/consumer/producer welfare, budget exhaustion and the post-CfD
period remain subsequent financial analysis, using the saved dispatch and prices.

See [planned-grid selection, provenance and coverage limits](bornholm_2035_grid.md).

Scenario resources have separate paths from earlier smoke, 90-cluster and historical
outputs. The native retrieval rule obtains the default archived weather cutout;
this configuration does not request a new 2023 ERA5 download. After pulling this
correction, rerun the workflow build before solving so weather-dependent profiles
and composed networks are regenerated for 2013. The checker rejects networks with
the old 2023 snapshots. Native source datasets that remain valid can be reused.

## Run on Sophia

Use the existing `~/bin/pixi-sophia` Singularity/Pixi wrapper and the working
Gurobi token-server licence at `~/gurobi.lic`. The job script exports the licence
path into the container; it does not copy credentials. Run from the repository
root on the Bornholm branch after applying these changes. If historical jobs are
still running, use a separate Bornholm checkout or worktree; do not change the
branch or files in their active checkout.

The wrapper requires a Slurm allocation even for a dry-run. From the login node,
request an interactive allocation without an explicit memory request:

```bash
cd ~/projects/pypsa-eur-bornholm
srun --partition=workq --time=00:30:00 --nodes=1 --ntasks=1 \
  --cpus-per-task=16 --pty bash -l
```

Once the compute-node shell opens:

```bash
cd ~/projects/pypsa-eur-bornholm
bash scripts/bornholm_sophia.sbatch dry-run
```

The dry-run resolves the complete dependency graph, so it may need access to
remote source-data manifests. After a successful dry-run, use `exit` to release
the interactive allocation and return to the login node. Build first without
solving:

```bash
sbatch scripts/bornholm_sophia.sbatch build
```

After building and passing the within-study comparison checks, solve the
Bornholm baseline pilot. No historical-benchmark network or `HINDCAST_NETWORK`
environment variable is required.

```bash
sbatch --partition=fatq scripts/bornholm_sophia.sbatch solve baseline
```

After the baseline pilot completes, solve another case, for example:

```bash
sbatch --partition=fatq scripts/bornholm_sophia.sbatch solve hybrid
```

Use the suffixes in the table above for the remaining cases. `build` and `dry-run`
default to all eight scenarios; `solve` defaults to the baseline pilot. A second
argument `all` explicitly selects all cases. `all <suffix>` builds and solves that
selection; a non-baseline case requires the solved baseline pilot first.
A selected build also builds its comparison baseline. Building all cases first
is the recommended sequence. The full solve stops before optimisation if
the actual electricity bus IDs or country assignments differ. Investigate an
upstream shape/network version mismatch rather than suppressing this check.

The batch script defaults to one exclusive `workq` node for builds; the solve
commands above override this to `fatq`. Both use 24 hours and 16 CPU workers;
Gurobi uses 16 threads. Sophia's observed `workq` configuration
reports only 1 MB per node to Slurm, which rejects explicit requests such as
`--mem=16G` or `--mem=112G`. This is scheduler metadata, not the node's physical
RAM. Omit `--mem` and `--mem-per-cpu`; the exclusive allocation avoids sharing the
node with another job. Check physical memory with `free -h` on the compute node.
The first hourly baseline failed with Gurobi out of memory during presolve on a
125 GiB `workq` node. A subsequent hardware check found 243 GiB on `fatq`.
This motivates the larger-node retry after rebuilding, but is not evidence that
243 GiB will be sufficient to finish the solve.
Snakemake's `mem_mb=110000` is an internal scheduling budget, not a Slurm memory
reservation or a hard runtime limit. It admits one full solve at a time. The
selected cases run sequentially within an allocation. Use one full solve per job
initially, then size the allocation from the pilot. The 24-hour allocation is a starting budget,
not a measured completion time or memory guarantee. Large sector models may
need a larger allocation; a walltime-limited solve does not automatically retain
an optimiser checkpoint. Finished network targets are reusable on a subsequent
submission. Do not start a second job against the same output directories while
the first is still running.

After a failed or cancelled build, wait for that job to exit before updating the
checkout and submitting `build` again. The script uses `--rerun-incomplete` to
regenerate unfinished outputs; completed, up-to-date targets can be reused.

### Restart after missing land-cell geometries

Build `4911929` at `89bd671a` failed in `build_electricity_demand_base` for all
eight scenarios: `rasterstats` passed a missing polygon to Shapely. This matches
the missing clipped-cell failure fixed on the hindcast branch in `98c6ba9`.
The same generic correction is included here. Clipped onshore Voronoi cells
with no positive land area are omitted from region exports; their electrical
buses, generators and transmission connections remain in the network. Substation
demand allocation filters existing land regions instead of reindexing them and
recreating missing polygons. Remaining missing geometries or countries without
eligible land cells produce an explicit error rather than silently losing load.

The native JRC, GB-statistics and NUTS source priority is retained. Combined
weights are normalised by country so allocation conserves every national hourly
load, including when GB and Northern Ireland use different sources. The two
affected production functions match the hindcast fix. These changes do not
alter the selected bidding-zone boundaries, pipeline settings or weather year.

After the failed build has stopped, pull the update and use `rebuild` once:

```bash
cd ~/projects/pypsa-eur-bornholm
git pull --ff-only origin feature/bornholm-energy-island

build_job=$(sbatch --parsable scripts/bornholm_sophia.sbatch rebuild)
pilot_job=$(sbatch --parsable --partition=fatq \
  --dependency="afterok:${build_job}" --kill-on-invalid-dep=yes \
  scripts/bornholm_sophia.sbatch solve baseline)
echo "Build: $build_job"
echo "Solve: $pilot_job"
squeue -j "$build_job,$pilot_job"
```

`rebuild` forces the native `base_network` rule for the selected scenarios and
rebuilds its dependants, then runs the comparison checks. It defaults to all
eight cases, does not solve, and reuses retrieved inputs and unchanged cutouts.
A failed rebuild returns a non-zero status, preventing the dependent baseline
solve. An old solve dependency cannot be reused after its build has failed;
submit a new dependent job as above. Ordinary subsequent runs can use `build`.
Local regression checks cover geometry coverage, retention of electrical buses,
national hourly demand conservation and job-stage selection. A complete revised
sector build and baseline solve still require the Sophia run.

### Restart after missing reservoir capacities

The rebuild at `eecee9a` reached `compose_network` on 27 October 2026, then
failed the reservoir capacity validation before duration completion. The
hindcast branch's correction from `d3e54e2` is also applied here: reservoir
records with missing discharge capacity are excluded before duration completion
and aggregation, with source row IDs, names, countries, buses and values logged.
Their omitted MW and storage are unknown; this is a source-data limitation.
The source CSV is retained, and capacity is never invented or set to zero.
Known-capacity plants still preserve reported durations and share residual
national storage through native capacity-weighted aggregation. Negative or
infinite capacities continue to fail, now naming the offending records.

After pulling this correction, use `build` to resume the unfinished composition
steps and reuse the corrected spatial inputs from the previous rebuild:

```bash
cd ~/projects/pypsa-eur-bornholm
git pull --ff-only origin feature/bornholm-energy-island

build_job=$(sbatch --parsable scripts/bornholm_sophia.sbatch build)
pilot_job=$(sbatch --parsable --partition=fatq \
  --dependency="afterok:${build_job}" --kill-on-invalid-dep=yes \
  scripts/bornholm_sophia.sbatch solve baseline)
echo "Build: $build_job"
echo "Solve: $pilot_job"
squeue -j "$build_job,$pilot_job"
```

Review the logged excluded records after the build. A successful build and
baseline solve remain to be verified on Sophia.

The script refuses all execution outside a Slurm job, including dry-run. `PIXI_SOPHIA` and
`GRB_LICENSE_FILE` can override the wrapper and licence paths. Record the job ID
returned by `sbatch`; the combined log is `slurm-bornholm-<jobid>.out`.

## Checks and results

`scripts/check_bornholm_networks.py` runs separately from model construction,
before and after solving. It checks hourly snapshots and weights, matched
onshore electricity buses and countries, sector buses, retained market links,
the reserved-area configuration, equal generic offshore potentials, equal project
wind profiles, fixed project capacities, cable endpoints, losses and investment
accounting. It also checks equal exogenous demands, the same initial background
grid, fixed background capacities in the main cases and the bounds on optional
project links. After solving it checks successful optimal termination, finite
objectives, synchronised cable capacities and equal prices within the two
shared-market pairs. Checks raise an
error instead of continuing with a mismatched comparison.

The resulting networks are:

```text
results/bornholm_zonal_<suffix>/networks/solved_2035.nc
```

Their hourly dispatch, prices, flows, sector conversion, curtailment and investment
provide the inputs for subsequent analysis. Use the same cost scope and snapshot
weights when comparing annual objectives, and distinguish system cost from
project revenue. A successful solver exit alone is not evidence that the study
assumptions or future prices are calibrated.

## Local verification

Before the weather-default correction on 25 October 2026, the expanded regression
suite passed all 84 collected tests;
the final job-script rerun passed all eight cases, including two subsequently
added guards (86 distinct tests overall). Coverage includes the original
project/resource implementation, native composition and clustering helpers,
sector coupling, reverse export flow, single-counted cable costs, shared-market
prices, rejection of mismatched comparison networks, the planned-grid filter,
replacement retirement and independent DK/DE/SE/PL cable investment. Tests use
Python 3.12.14, PyPSA 1.2.4, atlite 0.6.1, linopy 0.8.0 and HiGHS 1.15.1 in a
development environment. This was not an exact reproduction of Sophia's Pixi
environment and does not validate full-model Gurobi runtime or memory use.

After restoring the native weather/load defaults, 18 targeted checks passed:
within-study comparisons, rejection of stale 2023 networks, consistency with
native default inputs, and job stages without a historical network. Native
workflow expansion also confirmed 2013 snapshots, the default archived cutout
for generic solar and project wind, `load.fixed_year: false`, the 2023 sector
statistical base year, and the 2035 planning horizon. The subsequent Sophia build
completed all eight full-hourly sector networks and their comparison checks at
`cb45ce41`; the baseline solve then failed with Gurobi out of memory.

```bash
pixi run python -m pytest --noconftest -q test/test_energy_islands.py test/test_energy_island_resources.py test/test_compose_network.py test/test_cluster_network.py test/test_add_electricity.py test/test_bornholm_zonal.py test/test_bornholm_grid.py test/test_bornholm_job.py
```

Workflow input/parameter expansion was checked for all eight scenarios, including
34 countries, sector mode, native bidding zones, hourly timing, the four explicit
project CSV dependencies, the 2035 grid year and 16 solver threads. The selected
input tables and audits reproduce byte-for-byte from the pinned converter
checkout. Shell syntax, Python lint/format, Snakemake formatting and diff
whitespace checks passed.

The grid-only dependency dry-run passed. A native build through
`base_extended.nc` checks the real OSM grid, project import and retirement step;
its mapping audit is discussed in the [grid note](bornholm_2035_grid.md).
This grid-only test is separate from the subsequent full build on Sophia.

For the fixed-pipeline and shared-preprocessing revision, 112 targeted tests
passed, including pipeline candidate switches, stale-network rejection, native
plant and AC/DC corridor assignment across a boundary inside a county, reservoir
energy conservation through aggregation, small-network optimisation, NetCDF
round trips, configuration/schema synchronisation and the 12 job-script cases.
Two existing sector/resource integration tests were excluded from this run after
their multiprocessing weather-profile fixture stalled in the development
executor; the full resource-profile suite was not rerun successfully. The local
runtime uses native CPU-affinity detection in Dask because its optional psutil
probe cannot resolve process IDs in the executor's PID namespace. No production
model functions were replaced for these tests.

The ten relevant geography, clustering, plant-mapping and reservoir functions
were compared structurally with hindcast commit `0b30b179` and are identical.
A check using the available native zone polygons and Swedish plant coordinates
assigns Ringhals 3/4 and Oskarshamn 3 to SE03 and keeps Forsmark 1/2/3 in SE03.
This checks the mapping method, not a completed 2035 network rebuild.

The eight sector networks must be rebuilt with the restored native gas connectivity
and corrected gas capacities, pass the updated comparison checker, and complete
a diagnostic baseline solve before a new full-year pilot on Sophia. Existing
geography and hydro corrections remain. No HPC job was submitted from this audit.

### Gas-storage correction and short solver diagnostic

Gas storage retains the native reported GWh-to-MWh conversion and 98th-percentile
clipping rule, but the percentile now includes **only modelled gas stores**,
including zero-capacity gas nodes. Previously the series was expanded to all
stores first; adding or removing CO2/H2/heat stores changed the clipping threshold
and hence mandatory gas-storage investments. This is a correction to the input
calculation and changes the old baseline's bounds. Rebuild both original and
simplified composed networks before comparing them. The comparison checks now
reject mismatched gas-storage bounds and costs.

The default diagnostic uses the **original baseline**, with existing capture options,
spatial CO2, CO2 transport, methanation, all 52 bidding zones and all other baseline
sector technologies retained. The previously corrected gas-capacity preprocessing
and restored native connectivity candidates are included. Existing gas pipelines
remain fixed. No emissions-budget relaxation, load shedding or invented imports
are introduced.

From the repository root on Sophia:

```bash
git pull --ff-only origin feature/bornholm-energy-island
sbatch scripts/bornholm_diagnostic.sbatch all 24 barrier
```

This requests one exclusive workq node, four CPUs and four hours of walltime.
It first backs up and rebuilds the gas inputs and original composed baseline using
`scripts/bornholm_gas_capacity_rebuild.sbatch`, then runs a separate 24-hour solve
with a one-hour Gurobi time limit. Barrier uses the native `gurobi-default` preset
with four threads and the diagnostic time limit. The main full-year solver preset
is unchanged. A failed rebuild prevents the solve.

The diagnostic preserves one-hour physical storage chronology and power/energy
capacity bounds. Annual flow budgets, generator energy quotas and capital costs
are prorated by hours/8760. Annual non-sequestered HVC feedstock is also prorated,
with mandatory depletion imposed at the diagnostic endpoint as in the native
annual construction. Geological storage's flow limit is added by native
`prepare_network` using the shortened horizon. Native supplementary constraints
are added by `create_optimization_model`; a configured external constraint hook
must exist and is not silently bypassed. This is a short-window feasibility and
solver diagnostic, not an annual investment result. Cyclic stores cycle over the
short window, so the test cannot establish seasonal or full-year feasibility.

Preflight checks reject undefined buses, conflicting capacity bounds, zero-capacity
existing gas pipelines, stale gas-storage bounds and stale sector flags. An
optimistic annual max-flow check also tests a necessary condition for fixed direct
and industrial gas demand: it deliberately relaxes upstream conversion inputs,
competition for fuels and variable gas consumption for heat/power. In particular,
unlimited methanation candidates are treated as unlimited supply in this check;
passing it does not establish energy-system feasibility.

Each submission writes a fresh directory under
`results/bornholm_diagnostic/<job-id>_<variant>_<hours>h_<method>/` containing `python.log`,
`solver.log`, `input.nc`, `summary.json`, and, only after optimal termination,
`solved.nc`. The summary records model size, elapsed time, process peak RSS, actual
solver options, source-network hash and the external constraint hook's hash when
present. Slurm accounting provides the job-wide memory observation. Failed
preflight/solves leave a summary and return a failing exit code.

After the first diagnostic passes, use the already rebuilt source for a week or
a PDHG comparison; these are subsequent submissions, not automatically launched:

```bash
sbatch scripts/bornholm_diagnostic.sbatch solve 168 barrier
sbatch scripts/bornholm_diagnostic.sbatch solve 24 pdhg
```

Do not infer full-year memory requirements or seasonal feasibility from an optimal
24-hour result. No cluster job was submitted from this development environment.

Development validation for this revision: 58 focused tests passed across the gas
preprocessing, gas/storage integration, scenario comparison, diagnostic accounting,
preflight and job-control suites. Two existing offshore weather integration cases
were deselected. Native solve-rule parameter expansion was checked, including the
resolved tracked `data/custom_extra_functionality.py` hook (a no-op in this
checkout). Lint, formatting, shell syntax and diff-whitespace checks passed.
Separately, a 24-hour smoke test of the supplied **old original** composed network
built 168,242 variables and reached optimal termination with HiGHS in 100.94 solver
seconds using the built-in sector constraints. This checks short-model
construction; it is not a rebuild of the corrected source or a Gurobi/Sophia
performance result.


### Paired original and simplified diagnostics

The diagnostic now accepts a fourth submission argument, `original` (default) or
`simplified`. The latter applies the existing
`config/config.bornholm-simplified.yaml` overlay, including pooled CO2, no CO2
network, no methanation, no hydrogen fuel cells and no solar thermal. DAC and SMR
with capture remain enabled in both. Both keep 52 electricity zones and the same
24-hour window, annual accounting adjustments, solver preset and time limit.
The simplified formulation is an alternative assumption; pooled CO2 also changes
the industrial capture electricity representation in this code revision.

A simplified `all` submission rebuilds the original reference and the simplified
network with corrected gas inputs, then runs the shared-input comparison checker
before solving. Backup directories include the variant name. A simplified `solve`
submission checks the existing pair without rebuilding. Variant-specific
preflight checks reject sources built for the wrong baseline.

Prepare both sources in one build job, then submit two solve-only jobs with the
same successful-build dependency. The two solves have no dependency on each other
and can run in parallel on separate allocated nodes. Scheduler availability
controls their actual start times. Each reads the prepared inputs and writes its
own result directory; neither solve rebuilds the networks.

```bash
git pull --ff-only origin feature/bornholm-energy-island
build_job=$(sbatch --parsable scripts/bornholm_diagnostic.sbatch build 24 barrier simplified)
sbatch --dependency=afterok:"$build_job" scripts/bornholm_diagnostic.sbatch solve 24 barrier original
sbatch --dependency=afterok:"$build_job" scripts/bornholm_diagnostic.sbatch solve 24 barrier simplified
```

Here `build ... simplified` prepares both the reference and simplified networks
and checks shared inputs, then exits without launching a solver. A failed build
or input comparison prevents both dependent solves from starting. Do not launch
another rebuilding job while either solve is reading the composed sources.
Each job has its own four-hour allocation; each solve has a one-hour solver limit.
Existing commands without the fourth argument still select the original
baseline. Neither job was submitted from development.

The `24` argument denotes 24 modelled hourly snapshots (the first day of the source
year), not a 24-hour job duration. `barrier` selects Gurobi's interior-point solver
(Method=2), rather than the Bornholm PDHG preset (Method=6). The model retains the
selected variant's technology set, with annual input accounting prorated as
specified above. This diagnostic tests feasibility and solver progress before
longer chronological/full-year tests.

Validation of paired selection: 24 targeted diagnostic, job-control and gas/storage
tests passed. Native solve-rule expansion selected the correct source network and
CO2/methanation flags for each variant; DAC and SMR with capture remained true in
both. Shell syntax, lint, format and whitespace checks passed. No full-year build
or solver performance claim is made by these checks.

Validation of the parallel submission path: all nine job-control tests passed,
including stopping after build/input checks without launching a solver for either
variant. Shell syntax, lint, formatting and whitespace checks passed. The native
solver/model formulation is unchanged by this scheduling update.

### Compute-node Git availability

Sophia build job 4920111 stopped during backup bookkeeping with exit 127 because
host `git` was absent from the batch PATH, before invoking the network workflow.
Revision bookkeeping now runs through the existing Pixi wrapper. The new
`scripts/bornholm_provenance.py` always hashes the scripts, workflow rules,
configuration and environment manifests, and records a Git commit/diff only when
Git is available. Missing/failing Git is recorded explicitly and does not stop a
model job. The diagnostic solver also uses this helper so it cannot fail later
on the same mandatory Git command. Model inputs, solver settings and dependencies
are unchanged by this fix.

Validation: 21 targeted provenance, diagnostic and job-control tests passed,
including execution of the real build shell script with Git removed from its
PATH and the provenance Python module running in the test runtime. The model
workflow was mocked in that regression test; no Sophia rebuild was performed.
Lint, formatting, shell syntax and whitespace checks passed.

### Full-year barrier pilots

After the successful 24-hour and 168-hour barrier diagnostics, full-year pilots
can use the already rebuilt original and simplified networks. The annual wrapper
uses explicit `--annual --hours 8760` mode in the same runner. It retains the
source annual constraints, costs, generator energy quotas, HVC inventories and
storage chronology; it does not apply the short-window prorating function or
rebuild source networks. Native solve preparation and supplementary constraints
remain in use. Source preflight, provenance and simplified/reference comparison
checks remain enabled.

```bash
git pull --ff-only origin feature/bornholm-energy-island
sbatch scripts/bornholm_annual_barrier.sbatch original
sbatch scripts/bornholm_annual_barrier.sbatch simplified
```

The two independent pilots can run in parallel, each on one exclusive workq node
with four CPUs, a 24-hour Slurm walltime and a 20-hour Gurobi solver limit. The
solver is the same native barrier preset used for the weekly reference, including
Method=2, Crossover=0 and BarConvTol=1e-5. The main Bornholm full-year PDHG preset is
not edited. Each job writes fresh results under
`results/bornholm_annual_barrier/<job-id>_<variant>_8760h_barrier/` and its launcher
log is `slurm-bornholm-annual-barrier-<job-id>.out`.

The runner records physical RAM and available RAM on the allocated node. After
Linopy model construction, Gurobi's `SoftMemLimit` is set to 75% of then-available
physical RAM in decimal GB, leaving headroom for other allocations. A CLI memory
cap can lower this value but cannot raise it. This parameter covers Gurobi's
tracked memory, not all Python/Linopy allocations; it is checked only at graceful
termination points, so an OS-level OOM remains possible. It changes the permitted
computational resources, not physical capacities or model assumptions.

The summary is checkpointed atomically at loading, preflight, preparation, export,
model construction and solving phases. A caught failure records the previous
phase and error. If the OS kills the process, the last checkpoint and Slurm
accounting are needed; a final summary cannot be guaranteed. Only an optimally
terminated run exports `solved.nc`. A solver time/memory limit is a stopped pilot,
not an accepted annual solution.

Memory failures often arise in presolve, ordering or factorisation, but annual
model construction/loading must happen first and memory can also grow later.
An early failure is not guaranteed. Inspect the phase log and `Barrier statistics`
before making runtime estimates or deciding on larger-memory hardware.

Development checks: 26 targeted annual/diagnostic/provenance/job-control tests
passed, including annual-input invariance, memory-budget units/headroom, a
simulated solver memory failure with preserved progress, and selection of the
correct annual variant without rebuilding. Lint, format, shell syntax and
whitespace checks passed. No full annual model was built in development and no
Sophia job was submitted here.
