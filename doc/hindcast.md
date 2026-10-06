# Historical price and flow benchmark

## What is implemented

This is the first validation stage before Bornholm scenario runs. Start by
benchmarking the authors' published **Hindcast-dynamic 2023** network. It already
contains the hourly demand, renewable availability, fixed capacities, and dynamic
generator marginal costs. A second command re-solves those inputs with the local
PyPSA installation. A separate configuration starts a rebuild with this branch.

These are distinct results:

1. Benchmarking the archived solved network checks our analysis against the
   published model. It does not validate our own new network construction.
2. Re-solving the archived inputs checks reproducibility with our PyPSA/solver.
3. Building a new historical network and applying the same benchmark validates
   the construction workflow that will later support the Bornholm analysis.

The CLI implements the paper's MAE/RMSE/SMAPE comparisons at hourly, daily and
weekly resolution. The default six-country selection targets Northern Europe.
It additionally calculates price spreads, border flows, signed bias, correlation,
coverage and load shedding. No arbitrary pass/fail threshold is imposed.

## Sources and provenance

- Paper: https://arxiv.org/abs/2606.16486
- Published networks, archive version 1:
  https://data.dtu.dk/articles/media/31248568/1
- Reference prices and loads:
  https://github.com/marco-saretta/pypsa-network-analyzer
  at commit `fd282cc8b1705b68aad0edd286d6b8f1c612ac5d`.
- Additional zonal prices and retrieval:
  https://github.com/martavp/pv_wind_historical_timeseries/tree/main/data/entsoe

The downloader selects only `hindcast-dyn-wy2023.nc`, verifies the archive MD5,
then records local SHA256 hashes and URLs in `sources-2023.json`. It downloads
about 165 MB for one network plus the reference tables; no ENTSO-E key is needed
for this first price benchmark. Files are cached under `data/hindcast/reference`.

The archived 2023 file was inspected during implementation: it contains **35
buses and 8,760 hourly snapshots**, has no extendable components, and records
PyPSA 0.33.2. The paper text describes 39 nodes. We reproduce the actual archived
file and record that difference, rather than silently changing its topology.
The file contains no load-shedding generators; its metadata sets transmission
losses to 2. The re-solve preserves these choices (it does not add load shedding).

## First run on your laptop

Run from the root of your existing PyPSA-Eur checkout, using its Pixi environment:

```bash
pixi run python scripts/hindcast.py fetch-reference --year 2023
pixi run python scripts/hindcast.py benchmark --config config/hindcast/benchmark.yaml
```

This uses the published solved network and creates `results/hindcast-dynamic-2023/`:

- `scores.csv`: raw price and spread errors, coverage, hourly/daily/weekly results;
- `price_*.png`, `spread_*.png`: weekly profiles and duration curves;
- `*_paired.csv`: exactly the matched hourly observations used in each comparison;
- `model_prices_raw.csv`, `model_border_flows_mw.csv`, `load_shedding_mw.csv`;
- `audit.json`: model/input hashes, mappings, settings and validation status.

**Flow extraction is not flow validation.** Without a configured observed-flow
file, the program clearly reports `NOT PERFORMED` for flow validation. The price
benchmark can run independently.

To re-solve the complete year from the archived inputs:

```bash
pixi run python scripts/hindcast.py solve \
  data/hindcast/reference/hindcast-dyn-wy2023.nc \
  --output results/hindcast-resolved-2023/network.nc \
  --solver gurobi --threads 4

pixi run python scripts/hindcast.py benchmark \
  --network results/hindcast-resolved-2023/network.nc \
  --output results/hindcast-resolved-2023/benchmark
```

Use `--solver highs --threads 1` if Gurobi is unavailable. On Sophia replace
`pixi run` with your existing `~/bin/pixi-sophia run` inside an allocated job,
using no more threads than allocated. The script uses full-year perfect foresight,
not a rolling horizon. It refuses capacity-expansion, multi-investment-period,
unit-commitment, temporally aggregated and static-cost networks. It preserves
embedded constraints and reads transmission-loss settings from network metadata.
The original network is never overwritten; choose a new output for repeat solves.
Custom constraints defined only in external Python callbacks cannot be recovered
from NetCDF. The archived inputs and current PyPSA version alone are therefore
not a guarantee of an identical optimisation problem; inspect the re-solve check.

## Historical power flows

Retrieve directional data with your existing ENTSO-E credentials in the environment
(`ENTSOE_API_KEY`). Do not commit the token. `entsoe-py` is already included in
this branch's Pixi dependencies and is needed for retrieval only.

```bash
pixi run python scripts/retrieve_hindcast_flows.py --year 2023 --kind physical
```

Then set in `config/hindcast/benchmark.yaml`:

```yaml
flows: data/hindcast/flows_physical_2023.csv
flow_kind: physical
```

Re-run the benchmark command. Border mappings in `flow-sources.yaml` sum
DK1–Germany and DK2–Germany into DK–Germany, and the two relevant Danish–Swedish
borders into DK–Sweden. SE–PL uses SE4–PL. Both directional series are queried and
subtracted; a missing series raises an error and is not assumed to be zero.
Quarter-hourly or half-hourly data are averaged only for complete hours. Review
coverage in the resulting JSON before interpreting the errors. If a source
changes temporal resolution within the year, split and normalise it explicitly;
the strict helper does not silently fill the coarser intervals.

Alternatively, export existing pipeline data as a UTC hourly wide CSV:

```text
utc_timestamp,DK->DE,DK->SE,SE->PL
2023-01-01T00:00:00Z,120.0,-50.0,200.0
```

Values are MW, positive from the first named region to the second. Model flow
uses the actual first-region branch terminal (`p0` or `p1`), including losses.
Only AC lines and DC electricity links are included, never storage converters.
Physical observations can differ through losses, measurement conventions, loop
flows, redispatch and outages. Check terminal conventions before interpretation.
Commercial day-ahead exchanges can be retrieved separately with
`--kind scheduled_day_ahead`; they are a distinct market diagnostic and must not
be labelled physical-flow validation.

## Geographic mapping and missing data

The paper file has country-level buses, including one DK, SE and NO bus each.
Prices for countries with multiple zones are aggregated using historical zonal
loads. All declared zones and loads must be present; a missing zone invalidates
that hour. Single-zone regions do not need load weights. A country DK price is
**not a DK2 price**. The default DK–DE spread is labelled accordingly.

For your own bidding-zone network, supply an explicit mapping in a copied config:

```yaml
regions:
  DK2:
    buses: [YOUR_ACTUAL_DK2_BUS_NAME]
    price_zones: [DK_2]
  DE_LU:
    buses: [YOUR_ACTUAL_DE_BUS_NAME, YOUR_ACTUAL_LU_BUS_NAME]
    price_zones: [DE]
```

Replace placeholders with actual bus names. Inspect `buses` and the benchmark
price CSV headers before mapping. Multi-bus model prices use demand weighting;
this is an aggregation statistic, not an implementation of zonal market clearing.
Marta's individual price files have `utc_timestamp,price_eur_mwh`; pivot/rename
them into the wide zonal table before using them here.

All timestamps are converted to UTC. Naive network timestamps are treated as
UTC. Duplicate timestamps and non-hourly inputs fail. The model snapshot index
is authoritative: omission of 29 February never shifts March dates. No price
clipping, interpolation, forward filling or backward filling is applied to
benchmark observations. Daily/weekly scores use identical paired hours on both
sides and require complete coverage of modelled hours by default. Partial
first/last weeks remain labelled calendar bins. Correlation is undefined for a
constant series. SMAPE uses 0 for a zero/zero pair and is reported alongside
absolute errors, because percentage errors are sensitive near zero.

## Rebuild using this PyPSA-Eur branch

`config/hindcast/build-2023.yaml` is a separately labelled starting point. It fixes
existing capacities, disables Bornholm and planned grid projects, uses hourly
2023 demand/weather and annual 2023 renewable capacity estimates, and clusters by
bidding zones. It uses 2020 technology assumptions (the source does not provide a 2023 cost
table), then replaces fuel/CO2 dispatch costs with historical 2023 inputs.
It disables CO2 quantity caps and both default cost-price pipelines
so that verified historical total marginal costs can be supplied explicitly.

Build the **composed** network, not the ordinary solved target:

```bash
pixi run snakemake --cores 4 \
  --configfile config/hindcast/build-2023.yaml \
  resources/hindcast-build-2023/networks/composed_2023.nc
```

First perform the same command with `--dry-run`. Full construction requires the
normal upstream datasets, weather cutout and any CDS credentials. The config was
validated against this branch's JSON schema; the complete data-building DAG has
not been executed in the implementation environment.

This rebuild is NOT an exact reconstruction of the paper: current topology and
plant databases, 2023 conventional capacities, bidding-zone clustering, and ERA5
solar differ. Removing planned projects does not reconstruct commissioning dates
of every asset already present in today's grid database. Audit the historical
network before calling the new model validated. The paper used SARAH3 solar,
and held its conventional/hydro fleet constant from a 2020 database. The original
archive's embedded metadata can be exported with `hindcast.py inspect` to guide
reconciliation. Unpushed local Bornholm changes were not inspected or modified.

To supply total generator costs, use an hourly CSV with one column per exact
generator name. The helper builds it from hourly historical fuel prices,
CO2 prices, network efficiencies/emission intensities and explicit VOM:

```bash
pixi run python scripts/build_hindcast_costs.py \
  resources/hindcast-build-2023/networks/composed_2023.nc \
  --fuel-prices data/hindcast/fuel_2023_hourly.csv \
  --co2-prices data/hindcast/co2_2023_hourly.csv \
  --vom data/hindcast/vom.yaml \
  --output data/hindcast/marginal_costs_2023.csv

pixi run python scripts/hindcast.py solve \
  resources/hindcast-build-2023/networks/composed_2023.nc \
  --marginal-costs data/hindcast/marginal_costs_2023.csv \
  --output results/hindcast-build-2023/network.nc --solver gurobi --threads 4
```

Required fuel columns depend on the carriers present: gas, coal, lignite, oil,
biomass and waste. Units are EUR/MWh_th. CO2 column `co2` is EUR/tCO2; VOM YAML
maps carriers such as `CCGT: 3.0` in EUR/MWh_el (illustrative, not an assumption).
Costs are `VOM + (fuel + CO2_price * emissions_per_MWh_th) / efficiency`.
Non-fuel technologies retain their input costs. The CSV **replaces** total
marginal costs, so emissions charges are not added twice. Source fuel/CO2 data
are not fabricated or substituted with monthly defaults. Obtain the authors'
processed daily/weekly inputs for close methodological replication, document
hourly expansion and currency basis, and retain their provenance.

## Tests

```bash
pixi run -e test python -m pytest -q test/test_hindcast.py
```

Tests cover price spikes/negative prices, zero SMAPE denominators, missing zonal
weights, DST and leap-day alignment, both branch orientations with losses,
subhourly flow completeness, generator efficiency/CO2 accounting, historical
capacity guards and a complete HiGHS solve/export/benchmark on a small system.

## Recorded implementation check (2026-10-06)

All 11 tests passed using PyPSA 1.2.4 and HiGHS 1.15.1; Ruff checks passed.
The real archived 2023 network was benchmarked with the default configuration:

| Country | Matched hours | Hourly MAE (EUR/MWh) | Hourly bias (EUR/MWh) |
| --- | ---: | ---: | ---: |
| Germany | 8,760 | 29.80 | +16.54 |
| Denmark | 8,759 | 34.37 | +19.21 |
| Sweden | 8,759 | 49.03 | +43.36 |
| Norway | 8,759 | 44.91 | +41.75 |
| Poland | 8,760 | 26.51 | +13.58 |
| Finland | 8,760 | 50.94 | +37.48 |

Bias is model minus observation. These are raw-price results for the published
network with the aggregation rules above, not evidence that our rebuilt model
is validated. The Nordic overestimation warrants investigation before using
the model for Bornholm market-value conclusions. No model load shedding was found.

The full-year re-solve was attempted but the execution environment killed the
process (exit 137) before the solver returned a solution. No successful full-year
re-solve is claimed. Run the documented command in an adequately provisioned
job. Historical flow validation remains pending observed-flow inputs or an
ENTSO-E API key. The complete native network rebuild also remains to be run.
