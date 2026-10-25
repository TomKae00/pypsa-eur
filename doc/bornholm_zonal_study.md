# Hourly, sector-coupled Bornholm comparison

`config/config.bornholm.yaml` now defines the full study. It replaces the former
90-cluster, 3-hour electricity-only configuration. The separate smoke configuration
is retained for debugging; it is not the analysis run.

## Study definition

| Setting | Full comparison |
|---|---|
| Spatial footprint | Same 34 countries and native bidding-zone polygons as `feature/historical-hindcast/config/hindcast/build-2023.yaml` |
| Shared electricity markets | DE00 + LUG1; IE00 + GBNI, with native lossless `copper` links |
| Time resolution | 8,760 hourly snapshots, 1 January–31 December 2023; no temporal reduction or rolling horizon |
| Planning and technology-cost year | 2035, overnight capacity expansion |
| Weather and electricity-load year | 2023; ERA5 at 0.3 degrees, with the same cutout extent as the historical branch |
| Sector demand base year | 2023 energy totals; native 2035 technology and end-use assumptions |
| Sectors | Electricity, heating, transport, industry, biomass, shipping, aviation, agriculture, hydrogen and gas infrastructure |
| CO2 constraint | Native relative budget: 25% of the configured 1990 reference emissions, covering the enabled sectors |
| Grid | Native optimisable transmission (`vopt`); native sector DC and distribution losses; AC loss approximation disabled |
| Renewable resource classes | Three for onshore wind, offshore AC/DC wind and solar |
| Project counterfactual | Project site reserved from generic offshore development in every case |

The historical configuration was checked on 24 October 2026 (GitHub file blob
`0c5d638d76524174b829e62eeaa6b971055f454a`). Its electricity bus layout is reused;
its fixed historical generation fleet, historical fuel/CO2 prices and fixed grid
are not transferred into the 2035 expansion study. Baseline-versus-project
differences are controlled comparisons within this study. Differences from the
historical branch also reflect the year, sectors, investment and policy assumptions.
This is a future scenario, not a calibrated price forecast.

The native bidding-zone implementation aggregates the electrical network by zone
and retains AC network constraints between zones. It does not reproduce an NTC or
flow-based market-clearing model. DE/LU and Ireland/Northern Ireland keep separate
country buses for sector data while sharing an electricity price through their
unconstrained links. Sector buses increase the **total** bus count. Radial and
hybrid also add `BEI_HUB`; the onshore electricity bus IDs must be identical to
the historical reference. The actual count is printed from the built networks,
not inferred from `n_clusters` (which is not a target in administrative mode).

| Run name | Project wind | Danish export to DKE1 | German export to DE00 |
|---|---:|---:|---:|
| `bornholm_zonal_baseline` | None | None | None |
| `bornholm_zonal_radial` | 3,000 MW | 3,000 MW | None |
| `bornholm_zonal_hybrid` | 3,000 MW | 1,200 MW | 2,000 MW |

Both project cases use the same native, site-specific wind profile. Exports use
the native bidirectional DC-loss formulation: two non-negative directional flows,
synchronised capacity variables, one investment charge per physical cable.
Efficiency is `0.98 * 0.977 ** (length_km / 1000)` in either direction. Cable
lengths remain hub-to-landing distance times 1.25 estimates. See
[resource and cost accounting](bornholm_native_resources.md) for the site and
engineering approximations. The offshore hub receives no heat/hydrogen demand.

Scenario resources have separate paths so earlier smoke, 90-cluster and historical
outputs cannot silently stand in for these networks. Native source datasets and
the annual weather cutout can still be reused. The full Europe 2023 cutout must
already exist or be built using the user's configured CDS access.

## Run on Sophia

Use the existing `~/bin/pixi-sophia` Singularity/Pixi wrapper and the working
Gurobi token-server licence at `~/gurobi.lic`. The job script exports the licence
path into the container; it does not copy credentials. Run from the repository
root on the Bornholm branch after applying these changes.

```bash
cd ~/projects/pypsa-eur
bash scripts/bornholm_sophia.sbatch dry-run
```

The dry-run resolves the complete dependency graph, so it may need access to
remote source-data manifests. To build first without solving:

```bash
sbatch scripts/bornholm_sophia.sbatch build
```

For the actual solve, point to the **built network from the historical branch**,
not its YAML file. If it remains at
`resources/hindcast-build-2023/networks/composed_2023.nc` in this checkout, the
script finds it automatically. Otherwise set its real absolute path; the path
below is a placeholder to replace. It must also be visible inside the container.

```bash
export HINDCAST_NETWORK=/absolute/path/to/historical/composed_2023.nc
sbatch scripts/bornholm_sophia.sbatch solve
```

Alternatively, after setting that reference, `sbatch scripts/bornholm_sophia.sbatch all`
builds and solves in one allocation. The full solve stops before optimisation if
the actual electricity bus IDs or country assignments differ. Investigate an
upstream shape/network version mismatch rather than suppressing this check.

The script requests one `workq` node, 16 CPUs, 112 GiB and 24 hours. Gurobi uses
16 threads; the Snakemake memory budget admits one full solve at a time. The
three cases run sequentially. The 24-hour allocation is a starting budget,
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
accounting. After solving it also checks finite objectives, optimised project
capacities and equal prices within the two shared-market pairs. Checks raise an
error instead of continuing with a mismatched comparison.

The resulting networks are:

```text
results/bornholm_zonal_baseline/networks/solved_2035.nc
results/bornholm_zonal_radial/networks/solved_2035.nc
results/bornholm_zonal_hybrid/networks/solved_2035.nc
```

Their hourly dispatch, prices, flows, sector conversion, curtailment and investment
provide the inputs for subsequent analysis. Use the same cost scope and snapshot
weights when comparing annual objectives, and distinguish system cost from
project revenue. A successful solver exit alone is not evidence that the study
assumptions or future prices are calibrated.

## Local verification

On 24 October 2026, all 68 tests in the focused regression suite passed. It covers
the original project/resource implementation,
native composition and clustering helpers, sector coupling, reverse export flow,
single-counted cable costs, shared-market prices and rejection of mismatched
comparison networks. Tests use PyPSA 1.2.4, atlite 0.6.1, linopy 0.8.0 and HiGHS;
they do not validate the full model's Gurobi runtime or memory requirements.

```bash
pixi run python -m pytest --noconftest -q test/test_energy_islands.py test/test_energy_island_resources.py test/test_compose_network.py test/test_cluster_network.py test/test_add_electricity.py test/test_bornholm_zonal.py
```

Workflow input/parameter expansion was checked for all three composed networks
and Gurobi solves, including 34 countries, sector mode, native bidding zones,
hourly timing and 16 solver threads. Two build-output paths were corrected to
accept the native string-valued dataset folders. Shell syntax and all four job
stages were checked using a mock wrapper, including reference/Slurm/CPU guards.
Python lint/format, Snakemake formatting and diff whitespace checks passed.
The full dependency dry-run was stopped after 60 seconds during DAG construction;
it is not reported as passed. No full European build or solve and no cluster
submission was performed in the development environment. Complete the dependency
dry-run with working remote source-data access on Sophia before submission.
