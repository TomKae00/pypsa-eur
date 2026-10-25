# Hourly, sector-coupled Bornholm comparison

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
| Sectors | Electricity, heating, transport, industry, biomass, shipping, aviation, agriculture, hydrogen and gas infrastructure |
| CO2 constraint | Native relative budget: 25% of the configured 1990 reference emissions, covering the enabled sectors |
| Main electricity grid | Native OSM 0.7 plus the covered TYNDP 2026 project subset for 2035; fixed background (`vfixed`) |
| Grid benchmarks | The same installed planned assets, with native further reinforcement (`vopt`) |
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

```bash
cd ~/projects/pypsa-eur
bash scripts/bornholm_sophia.sbatch dry-run
```

The dry-run resolves the complete dependency graph, so it may need access to
remote source-data manifests. To build first without solving:

```bash
sbatch scripts/bornholm_sophia.sbatch build
```

After building and passing the within-study comparison checks, solve the
Bornholm baseline pilot. No historical-benchmark network or `HINDCAST_NETWORK`
environment variable is required.

```bash
sbatch scripts/bornholm_sophia.sbatch solve baseline
```

After the baseline pilot completes, solve another case, for example:

```bash
sbatch scripts/bornholm_sophia.sbatch solve hybrid
```

Use the suffixes in the table above for the remaining cases. `build` and `dry-run`
default to all eight scenarios; `solve` defaults to the baseline pilot. A second
argument `all` explicitly selects all cases. `all <suffix>` builds and solves that
selection; a non-baseline case requires the solved baseline pilot first.
A selected build also builds its comparison baseline. Building all cases first
is the recommended sequence. The full solve stops before optimisation if
the actual electricity bus IDs or country assignments differ. Investigate an
upstream shape/network version mismatch rather than suppressing this check.

The script requests one `workq` node, 16 CPUs, 112 GiB and 24 hours. Gurobi uses
16 threads; the Snakemake memory budget admits one full solve at a time. The
selected cases run sequentially within an allocation. Use one full solve per job
initially, then size the allocation from the pilot. The 24-hour allocation is a starting budget,
not a measured completion time or memory guarantee. Large sector models may
need a larger allocation; a walltime-limited solve does not automatically retain
an optimiser checkpoint. Finished network targets are reusable on a subsequent
submission. Do not start a second job against the same output directories while
the first is still running.

The script refuses non-dry execution outside a Slurm job. `PIXI_SOPHIA` and
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
statistical base year, and the 2035 planning horizon. Full European profiles and
sector networks have not been rebuilt for the corrected weather input here.

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
This is not a full European sector build or solve. The earlier full sector
dependency dry-run was stopped during DAG construction and is not reported as
passed. No cluster submission was performed here. Complete that dependency
dry-run, all composed-network comparison checks and the baseline pilot on Sophia.
