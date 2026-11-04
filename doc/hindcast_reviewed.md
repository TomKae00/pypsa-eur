# Reviewed 2023 hindcast

The `feature/historical-hindcast` branch now builds the reviewed case under
`hindcast-reviewed-2023`. It retains all 34 study countries, 8,760 hourly 2023
snapshots, fixed historical capacities, the native demand/weather series and
the existing historical fuel/CO2 cost procedure. Inputs must be recomposed;
re-solving an earlier network does not install these corrections.

## Changes

The shared [transmission specification](transmission_assumptions.md) adds dated
asset rules and operating limits to the physical network. Historical 2023 limits
are kept separate from 2035 assumptions. The Nordic annual maxima are upper
envelopes, not an hourly market reconstruction. Core MaxBEX diagnostics are not
installed as independent NTCs. Complete historical Core flow-based and Nordic
joint/hourly availability remain outstanding.

Danish onshore availability is estimated independently for DK1 and DK2 using
Energinet's `ProductionConsumptionSettlement`, `CapacityPerMunicipality` and
`Elspotprices` datasets. The retrieval validates every UTC hour, handles the
price areas separately and records raw-source/output hashes. Municipality codes
map the reported monthly fleet to the two areas. Monthly capacity is a proxy
for within-month commissioning, not an exact dated asset register.

Existing wind-class capacity fractions are preserved within each area. The
fleet is scaled to the maximum reported monthly capacity, while the monthly
capacity fraction enters hourly availability. A positive scalar is fitted to
observed energy in hours with prices at least EUR 5/MWh and native fleet CF at
least 0.05. At least 1,000 fit hours are required, and the factor must be within
0.3–1.5. The fit retains weather chronology, zero-output hours and class ordering;
there is no constant-CF replacement. Factors are computed from the newly built
profiles rather than copied from the Bornholm 0.60 scenario multiplier.

This is **empirical historical effective availability**. Observed production
includes curtailment and outages; filtering low-price hours reduces, but cannot
eliminate, that contamination. Validate simulated generation and curtailment
after solving, and do not interpret the fitted factors as independent estimates
of meteorological bias or transfer them directly to future potential sites.
Offshore profiles and the native 0.8855 offshore correction are unchanged.

A 2023 Swedish reactor guard relocates any nonzero fixed SE4 nuclear capacity
to SE3 without changing total capacity or other countries. Already-correct
SE3 inputs are unchanged. It does not reconstruct dated reactor outages.
No Bornholm coal-store assumption is imported into this electricity-only case.

## Run on Sophia

```bash
git pull --ff-only
sbatch scripts/hindcast_reviewed_workq.sbatch dry-run
sbatch scripts/hindcast_reviewed_workq.sbatch all
```

The controller reserves one exclusive workq node, requests 16 solver CPUs and
uses all node memory (`--exclusive --mem=0`). Sophia advertises placeholder
memory to Slurm, so a positive request such as `--mem=100G` is unsuitable.
The internal Snakemake build-concurrency budget remains 90,000 MB. It builds
`resources/hindcast-reviewed-2023/networks/composed_2023.nc`, then solves to
`results/hindcast-reviewed-2023/network.nc`. Required fuel/CO2 files from the
existing procedure and ERA5 credentials must already be available. The Danish
observations are retrieved automatically as an explicit workflow dependency.
Use `build` or `solve` instead of `all` to run a single stage. Existing accepted
outputs cannot be overwritten.

The solve requires optimal termination and checks corridor/joint limits, bus
balances, generator availability and chronological storage equations before
export. It retains source hashes and a physical-quality report in metadata.
The zonal mapping is prepared using the existing exact aliases. Review any
new mapping template or unresolved zones before interpreting the benchmark.

With the observed price file and reviewed bus mapping present:

```bash
pixi-sophia run python -m scripts.hindcast benchmark \
  --config config/hindcast/benchmark-reviewed-2023.yaml
```

Compare hourly/monthly prices, border flows, generation, installed capacities,
wind curtailment, hydro/storage operation and scarcity hours. A factor fitted
to historical production can still yield too little dispatched wind if the
model curtails it again; report that separately from weather error.

The original benchmark configurations/results remain reproducible. The 2022,
2024 and 2025 cases are not silently updated with 2023 transmission values.
Full dated plant outages, reserves, hydro foresight/cyclicity, CHP operation and
the Core flow-based feasible domain remain further validation tasks. This
revision improves the run foundation; it does not certify a market replica.

Local verification covers small solved systems, UTC source completeness and
workflow configuration/dependencies. No full annual Sophia result is claimed
by this code update.
