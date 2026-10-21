# Historical price and flow benchmark

## Independent feature branch

This workflow lives on `feature/historical-hindcast`, based directly on `master`
(PyPSA-Eur configuration version v2026.08.0). Only the historical configuration,
helper scripts, documentation and existing test file were transferred from the
previous draft branch. It does not import energy-island scripts or modify the
standard composition rules. Bornholm's feature branch remains separate.

The market-price convention now uses the native `cluster_network` copperplate
option; its existing helper has guards for invalid bus groups. There is no
separate network-merging script.

The build config does not need an `energy_islands` switch: that feature is not
installed in this branch. No tests, simulations or benchmark validation were run
for this migration.

## What is implemented

The main workflow is **four independent European bidding-zone hindcasts for
2022, 2023, 2024 and 2025**, before Bornholm scenario runs. `build.yaml` contains
shared assumptions; `build-YEAR.yaml` selects each independent annual case.
Each case includes the default 34-country study footprint and hourly inputs. It does not inherit the four-country smoke-case selection.
This footprint does not literally include every European state: extensions to
Ukraine, Moldova or Turkey require corresponding historical network and input
data. Keep the build and benchmark country lists consistent if extending it.

`benchmark.yaml` is now the main zonal benchmark. Every electricity bus must be
mapped to its historical bidding zone. All mapped zones and all modelled
inter-zone AC/DC connections enter the reporting; missing observations are
reported as unscored, never used to remove a country from the optimisation.
Separate DK1/DK2, SE1–SE4 and NO1–NO5 are required. Italian and other market zones
must also use their historical boundaries and identifiers in the bus mapping.

The authors' archived 2023 solution is a separate reference exercise using
`benchmark-reference.yaml`. Its six-country selection affects only reported
metrics: the archived optimisation itself covers its full European network.
Those six selections did not truncate the simulated system or cause the price
errors. However, the archive has country-level DK/SE/NO buses and cannot be made
into a bidding-zone model by renaming or duplicating their prices.

Building and solving our historical zonal network remains necessary. Benchmarking
or re-solving an archived country-level network does not validate that new build.
The CLI reports raw hourly/daily/weekly MAE, RMSE, SMAPE, bias, correlation, coverage,
price spreads, border flows and load shedding without arbitrary pass thresholds.

## Native market-zone convention

The goal is **one hourly market price per zone**, not necessarily one physical
bus. Native bidding-zone clustering is retained with
`allow_ac_dc_mixing_in_bus_clusters: true`. Native `copperplate_regions` connects
the following expected clustered IDs:

| Market price area | Native cluster IDs | National bus information retained |
| --- | --- | --- |
| Germany–Luxembourg (`DE_LU`) | `DE00`, `LUG1` | DE and LU |
| Ireland/Northern Ireland (`IE_SEM`) | `IE00`, `GBNI` | IE and GB |

These IDs follow `build_bidding_zones.format_names` in this branch. They are
configuration expectations, not inspected outputs from an executed 2023 build.
If the actual shape/network combination uses different IDs, clustering raises
an error listing the unknown and available buses. Adjust the group IDs and exact
aliases together; a misspelt group must not silently leave a market separated.
Other islands or exceptional regions need explicit historical assignment;
country-prefix matching is not used.

Copperplating retains buses and their countries, removes their internal branches,
and connects them by fixed, infinite-capacity, bidirectional, lossless, zero-cost
`copper` links. Existing national input assignment is retained. Inter-zone AC
lines and DC links remain in the standard pipeline. These represent aggregated
electrical constraints, not a reconstruction of historical NTC or flow-based
market capacity allocation. Internal congestion and redispatch are not modelled.

The main benchmark requires buses sharing an observed price to be connected by
these native copperplate links. A copperplate group cannot be split across
different benchmark zones. After solving, bus-price spread within a zone must
be no larger than `price_atol + price_rtol * max(abs(bus prices))` per hour.
Defaults are 0.01 EUR/MWh absolute and 1e-6 relative. Disagreement is an error;
prices are never averaged to conceal it. The lexicographically first bus supplies
the reported price, and an hour is missing if any member's price is unavailable.
`audit.json` records the representative buses and maximum within-zone spreads.
The optional archived country-reference comparison keeps its explicitly labelled
demand-weighted aggregation.

`zone-aliases.yaml` supplies exact native-ID-to-price-area aliases for the
mapping template. It is a naming aid, not a certified historical shape catalogue.
Pin/review the source shapes and reconcile boundaries with 2023 before claiming
historical spatial fidelity. No tests or model construction were run for this
revision.

## Multi-year workflow: 2022–2025

Each year is solved independently, with its own snapshots, demand year, ERA5
cutout, standard plant-date filter and IRENA capacity year. Countries, technology
assumptions, model formulation and native copperplate settings are shared.
The native grid datasets are used as supplied; no year-specific reconstruction
of every transmission asset or observed-capacity calibration is introduced.
There is no investment or storage-state coupling between the annual cases.
Do not put all four years into one `planning_horizons` list in overnight mode.

### 1. Prepare observed electricity prices

```bash
pixi run python scripts/hindcast.py fetch-prices --years 2022 2023 2024 2025
```

This downloads only `day_ahead_price_*.csv` from Marta's pinned repository commit
`ad89e2d8e8c0c8b1c7e717dded8ebd43d52b3beb`. It does not run the ENTSO-E query
script or require an API key. Source files are cached once under
`data/hindcast/source/COMMIT/`. One wide UTC-hourly table and a provenance/coverage
JSON are written for each year:

- `data/hindcast/prices_zonal_2022.csv` and `.json`;
- `data/hindcast/prices_zonal_2023.csv` and `.json`;
- `data/hindcast/prices_zonal_2024.csv` and `.json`;
- `data/hindcast/prices_zonal_2025.csv` and `.json`.

No interpolation, clipping or neighbouring-zone substitution is performed.
Model zones absent from the source (for example IT_CALA or GB in this pinned
folder) remain all-missing columns. Each JSON reports source-file presence,
observed hours and coverage separately by zone and year. Existing annual outputs
are not overwritten; choose another `--directory` for a separately pinned refresh.
The importer uses already-hourly published data: underlying subhourly gaps cannot
be reconstructed from those files. In 2024 it retains all 8,784 observation hours;
the model currently omits 29 February and benchmarks against its own 8,760 timestamps.

These are observed day-ahead electricity prices, not generator fuel/CO2 inputs.
Retrieve the latter independently as described next.

### 1b. Retrieve historical fuel and CO2 source data

```bash
pixi run python scripts/retrieve_hindcast_market_data.py
```

The defaults in `config/hindcast/market-data.yaml` retrieve 2022–2025, plus
December 2021 so that a later, explicit carry-forward policy can cover the
beginning of 2022. This step does not build or solve a network. It uses existing
environment dependencies and public endpoints, without an ENTSO-E token,
market-data account, `yfinance` installation or Trading Economics subscription.
Live endpoint access is not guaranteed: Investing.com can return HTTP 403 from
some connections even when the same public request succeeds elsewhere.

| Input | Retrieved source | Quoted units and convention |
| --- | --- | --- |
| Gas | Yahoo Finance `TTF=F`, daily close | EUR/MWh as quoted; exchange-local trading date |
| Coal | Investing.com Rotterdam coal futures, instrument 1112838, daily close | USD/metric tonne |
| CO2 | Investing.com ICE carbon-emissions futures, instrument 8848, daily close | EUR/tCO2 |
| Oil products | Commission Weekly Oil Bulletin history, `Prices wo taxes`, all six EU aggregate products | EUR/tonne for fuel oils; EUR/1000 litres for petrol, diesel, heating oil and LPG |
| Currency | ECB euro foreign-exchange reference history | USD per EUR |

Coal uses the provider and market page cited in the paper. For gas, the analyzer
repository's plotting script uses `TTF=F`, while the paper cites `NG=F`; the exact
model-input contract and rolling convention remain unconfirmed. The CO2 source
is an **explicit substitution** for the paper's Trading Economics series. It is
an EUA futures price history, not a primary-market auction-price substitution.
Both futures histories retain the provider's continuous-series convention;
the workflow does not reconstruct individual delivery contracts or their rolls.

The output directory is `data/hindcast/market-data-2022-2025/`, containing:

- `gas_ttf_daily.csv`: `date,gas_eur_mwh_quoted`;
- `coal_rotterdam_daily.csv`: `date,coal_usd_t`;
- `co2_eua_daily.csv`: `date,co2_eur_tco2`;
- `oil_eu_weekly.csv`: dates and six separately named EU oil-product columns;
- `fx_daily.csv`: `date,usd_per_eur`;
- `raw/`: original response bytes and download receipts;
- `manifest.json`: configuration snapshot, instrument identity, source URLs,
  retrieval times, SHA256 hashes and per-year observation counts/date ranges.

Downloads and normalized tables are ignored by Git. Running the command in
another checkout retrieves that machine's own collection. Sources can revise
their histories; retain the raw files and manifest used for any published run.
An interrupted retrieval reuses downloaded sources with matching receipts and
hashes. A completed collection is not overwritten: choose a new `--directory`
for a refresh. `--years` can select another set of completed historical years.
Provider schema/access failures stop retrieval; no alternative ticker or provider
is substituted automatically.

For an existing source collection, an offline ZIP bundle can transfer the
original responses and their download receipts without querying any provider:

```bash
pixi run python scripts/retrieve_hindcast_market_data.py \
  --source-bundle "$HOME/Downloads/hindcast-market-sources-2022-2025.zip" \
  --directory data/hindcast/market-data-2022-2025-offline
```

The bundle contains `raw/` entries using the same URL-derived filenames and
SHA256 receipts as the source cache. Source URLs, sizes and hashes are checked
before import, and original retrieval dates are retained. Files are read by
their exact expected names, without extracting archive paths. All five inputs
must be present for the configured range; offline mode never falls back to live
requests or another provider. The manifest records the bundle SHA256 too.
Matching cached files are reused. If an existing cache has different source
bytes, choose a new `--directory` to preserve both versions.

The 403 error now reports this offline option directly instead of displaying a
Python traceback. The pandas date slice uses explicit timestamps for pandas 3
and later. These changes do not alter source prices or their units.

These are **source histories, not the prepared hourly CSVs consumed by
the native composition workflow**. Retrieval keeps missing source values and non-trading
days, without interpolating, filling, rounding CO2, converting currencies,
deflating prices or selecting an oil fuel. The `date` column is an observation
date, not an assertion that the closing quote was available at midnight.

### 1c. Apply the procedure recovered from the paper's 2023 archive

```bash
pixi run python scripts/prepare_hindcast_cost_inputs.py
```

The default in `config/hindcast/cost-inputs.yaml` now follows the implemented
fuel-price procedure inferred from `hindcast-dyn-wy2023.nc`. It replaces the
previous additional gas HHV/LHV adjustment, daily ECB coal conversion and
40.4 GJ/tonne oil conversion. Output goes to **`data/hindcast/paper-procedure/`**,
so earlier files under `data/hindcast/` are preserved and must not be mixed into
these revised runs. Use `--source-directory` for another source collection,
`--years` for a subset and `--directory` for a different output location.

| Input | Implemented procedure | Evidence and limits |
| --- | --- | --- |
| Gas | Unmodified Yahoo TTF close; multiplier 1.0 | All 8,760 hours match the recovered 2023 input; no additional division by 0.9 |
| Coal | Rotterdam USD/tonne quote rounded to two decimals, multiplied by 0.1141439205955335 | All 8,760 hours match the recovered 2023 input; this is one combined currency/energy factor |
| Oil | National untaxed low-sulphur fuel-oil prices aggregated with published 2023 consumption weights, divided by 10.83 | Matches 357 of 365 reference days, or 8,568 of 8,760 hours |
| Lignite | Converted coal price multiplied by 0.8911423821185138 | Constant ratio across every hour of the 2023 archive |
| CO2 | Investing.com EUA futures retained without scaling or integer rounding | Explicit provider substitution; the original Trading Economics series remains unresolved |

The gas match establishes equality to the quote, not an independently certified
HHV/LHV basis. The combined coal factor must not be interpreted as a standalone
calorific value, and **ECB rates are not applied in addition to it**. The raw ECB
history remains in the retrieved collection, but is unused by this procedure.
These choices reproduce input-cost behaviour; they are not fitted to electricity
prices and do not change native generation assets, efficiencies or intensities.

The oil calculation reads the **original Commission workbook** from the
collection's `raw/` directory. Using its already published EU aggregate would
not reproduce the identified weighting procedure. For each observation date:

1. Read national `*_price_wo_tax_fuel_oil_1` prices from `Prices wo taxes`.
2. Read the fixed 2023 `*_consumption_fuel_oil_1` weights from `Consumption`.
3. Over the configured EU-27 countries with positive reported consumption and
   an available price, compute `sum(price * consumption) / sum(consumption)`.
4. Divide that aggregate EUR/tonne quote by 10.83 and hold it until the next
   usable observation.

Countries with unreported consumption are excluded; weights are renormalized
when a national price is missing. The workbook's source hash, actual weights,
omitted countries and number of countries contributing to each hourly value
are recorded. No future national price is used to fill another country's gap.

The unresolved 2023 oil dates are **1 January and 27 November–3 December**.
They remain generated from the source-based procedure: no hand-picked patches
or values copied from network outputs are introduced. The data comparison
establishes the preceding match counts, but does not establish the reason for
those eight dates. Likewise, the recovered CO2 series contains whole-euro
values, but rounding or shifting the retrieved Investing.com series does not
reproduce it. The script therefore keeps that provider substitution explicit
instead of inventing a correction to force a match.

**Parameter evidence is currently limited to 2023.** The same recovered coal
and lignite factors, oil divisor and fixed 2023 oil-consumption weights are
applied to 2022, 2024 and 2025 as an explicit extension. Their exact equivalence
to the authors' other archived years has not been established, and no archived
2025 case exists. Every annual JSON records this distinction and sets
`exact_replica: false`, together with the CO2 source difference and unresolved
oil intervals. This is an alignment of the known fuel procedure, not a claim
that every paper input or parameter has been replicated.

Daily and weekly prices are expanded using only the latest observation with an
effective date on or before each model hour. Missing source rows are not new
observations. There is no interpolation or backward filling from a future date;
December observations seed the year start. Daily prices expire after 10 calendar
days and weekly oil after 21 days, including the end of the last allowed day.
Missing seed observations or excessively long gaps stop preparation.

This is a **retrospective full-foresight convention**: daily closes and weekly
observations are assigned to 00:00 UTC on their labelled observation dates
(`effective_lag_days: 0`). Those timestamps are not publication times; the
result must not be presented as a causal real-time backtest. No publication
dates are inferred. Explicit non-negative lags can move effective dates forward
for a separately defined timing convention.

Each selected year produces the following under `data/hindcast/paper-procedure/`:

- `fuel_YEAR_hourly.csv`: `utc_timestamp,gas,coal,lignite,oil`, nominal
  EUR/MWh thermal following the implemented archive price convention;
- `co2_YEAR_hourly.csv`: `utc_timestamp,co2`, nominal EUR/tCO2;
- `cost_source_dates_YEAR_hourly.csv`: observation/effective dates and the
  number of contributing oil-price countries;
- `cost-inputs_YEAR.json`: source/configuration hashes, conversion factors,
  national oil weights, holding ages, annual statistics, output hashes and
  the remaining alignment limitations.

The source CSVs and raw workbook must match their collection-manifest hashes.
There are 8,760 hours in 2022, 2023 and 2025, and 8,784 in 2024. Preparation
retains leap day; the native input loader subsequently selects the network's
actual timestamps, so omitting 29 February in the build causes no date shift.
Existing outputs are never overwritten, and every selected year is prepared
before export. These files supply the four traded fuels and CO2. PyPSA-Eur's
native composition code retains its own VOM, efficiencies, emission intensities
and cost-table fuel assumptions for other carriers; no separate VOM YAML or
biomass/waste price file is required for this build.

### 2. Construct the independent yearly systems

Both build configurations enable PyPSA-Eur's native dynamic fuel and emissions
pathways with these optional external input fields:

```yaml
conventional:
  dynamic_fuel_price: true
  fuel_price_file: data/hindcast/paper-procedure/fuel_{horizon}_hourly.csv
costs:
  emission_prices:
    enable: true
    dynamic: true
    co2: 0.0
    price_file: data/hindcast/paper-procedure/co2_{horizon}_hourly.csv
```

`{horizon}` resolves to each independent case's year (2022–2025). The files are
explicit inputs of `compose_network`, so they must exist before composition and
Snakemake tracks changes to them. Prepared inputs bypass the monthly World Bank
fuel builder and the Instrat CO2 builder, including their rolling averages. There
is no additional conversion, smoothing or price filling during composition.
With these optional fields unset, the standard PyPSA-Eur price pathways remain
available with their existing behaviour.

The loader converts timestamp labels to UTC, matches network snapshots exactly,
and rejects gaps, duplicate timestamps and non-finite prices. It selects the
actual 2024 snapshots without shifting dates when leap day is omitted. Missing
files or hours cause an error rather than substitution with another source.
The resolved paths and SHA256 hashes are stored in the composed network's
`meta.native_price_inputs` and retained by the hindcast solve.

Native generator attachment maps `gas` to CCGT/OCGT and calculates fuel cost
divided by each generator's efficiency plus its native VOM. Other carrier fuel
values come from the native processed technology costs. Native network
preparation then adds the hourly CO2 charge using carrier emission intensities
and generator efficiencies. Its dynamic/static branch applies CO2 once.
No marginal-cost formula is reimplemented in the hindcast build.

Use separate native workflow invocations, merging the shared configuration with
one annual configuration. This branch rejects scenario overrides that change
`planning_horizons`, so the annual cases cannot be run as entries in one native
scenario file. `scenarios.yaml` is retained as a reference catalogue only.
The 2023 configuration also works on its own.

```bash
for year in 2022 2023 2024 2025; do
  pixi run snakemake \
    "resources/hindcast-build-${year}/networks/composed_${year}.nc" \
    --cores 4 --configfile config/hindcast/build.yaml "config/hindcast/build-${year}.yaml" || break
done
```

These commands are for when inputs and assumptions are ready; they have not been
run here. Make shared changes in `build.yaml` and annual changes in the matching
`build-YEAR.yaml`. The standalone `build-2023.yaml` also contains shared settings,
so keep it aligned when changing those assumptions.

### Historical weather cutouts

The hindcast configurations explicitly set `data.cutout.source: build` (version
`latest`) and use PyPSA-Eur's native `build_cutout` rule with atlite. Annual names
such as `europe-2023-era5` describe the local cutout configuration; they are not
promises that matching files exist in the default download archive. Leaving the
source at its default `archive` causes a `MissingInputException` for the remote
`.../cutout/v1.0/europe-2023-era5.nc` URL before any build starts.

The native version registry currently resolves local build output to
`data/cutout/build/unknown/europe-YEAR-era5.nc`. Monthly ERA5 requests are enabled
and the time chunk size is 500. All ERA5 features are prepared: the hindcast does
not restrict preparation to the default example's temperature/height/runoff
subset, which would omit weather variables needed by wind and solar.

Building requires CDS API access configured on the machine running the workflow;
see [PyPSA-Eur's weather-data instructions](data-cutouts.md). The first full-Europe
annual build downloads and prepares weather data before generation profiles can
be built. Completed cutouts are reused by later workflow runs.

To prepare just the 2023 cutout before continuing the composed-network build:

```bash
pixi run snakemake data/cutout/build/unknown/europe-2023-era5.nc \
  --cores 4 --configfile config/hindcast/build-2023.yaml
```

No cutout download, weather preparation or model run was performed for this fix.

**2025 input availability is not assumed.** The branch's native
`data/nuclear_p_max_pu.csv` currently ends in 2024; its nuclear availability loader
selects the model year. A 2025 run therefore needs a supplied 2025 column from the
upstream data/source or a separately documented explicit assumption. No 2024
values have been silently substituted. Likewise, availability of 2025 demand,
IRENA capacities and processed fuel/CO2 prices must be established before solving.
This is an input-coverage requirement, not calibration to observed prices or
replacement of the native plant database. The paper's downloadable solved
reference networks cover 2020–2024, so `fetch-reference` deliberately has no 2025 case.

### 3. Prepare a market mapping for each year

```bash
for year in 2022 2023 2024 2025; do
  pixi run python scripts/hindcast.py prepare-zonal \
    "resources/hindcast-build-${year}/networks/composed_${year}.nc" --year "$year"
done
```

The first invocation for each year creates `data/hindcast/bus_zones_YEAR.csv`
using exact aliases in `zone-aliases.yaml`. Review every assignment and fill any
unknown IDs from that year's generated shapes/bus maps. Existing mappings are
never overwritten. Native shared-market groups must have the same benchmark zone:
DE/LU maps to DE_LU, and IE/Northern Ireland maps to IE_SEM. Re-run after review
to generate `data/hindcast/flow-sources-zonal-YEAR.yaml` for every modelled border.
The alias file accepts 2022–2025 as a naming convention; it does not certify that
source geometries reproduce the historical boundaries in every year.

### 4. Solve each case with its native composed costs

After composition and mapping review, solve each selected year's network. The
example uses 2023:

```bash
year=2023
pixi run python scripts/hindcast.py solve \
  "resources/hindcast-build-${year}/networks/composed_${year}.nc" \
  --output "results/hindcast-build-${year}/network.nc" --solver gurobi --threads 4
```

The composed network already contains historical fuel and CO2 costs calculated
by PyPSA-Eur. Do not pass `--marginal-costs` for this route: that option replaces
the total generator costs. `build_hindcast_costs.py` remains an optional manual
or archived-reference override and is not part of the main build.

Rebuild previously composed networks with the updated build configuration
before solving; older files may contain static or manually replaced costs.
No own-network solve is run by the price importer or benchmark. On Sophia use
an allocated job and your established Pixi/Gurobi wrapper as described below.

### 5. Compare annual prices and errors

```bash
pixi run python scripts/hindcast.py benchmark-years --years 2022 2023 2024 2025
```

The shared `benchmark.yaml` expands `{year}` in network, mapping, observation
and output paths. The command requires all selected annual inputs, compares each
year separately and writes `results/hindcast-europe-zonal/scores_all_years.csv`.
The table retains year, zone, resolution, MAE/RMSE/SMAPE, bias, correlation and
coverage, including explicitly unscored observations. Annual errors are not
pooled into one four-year score, so gaps or high-price years are not hidden.
All selected annual comparisons must finish before a combined table is written.

For a single case:

```bash
pixi run python scripts/hindcast.py benchmark --year 2024
```

The default for single-year commands remains 2023. Existing archived-reference
configs continue to work without a year override. The European aggregate and
generation-comparison extensions discussed separately remain pending.

Physical flows are optional and do not gate price comparison. Retrieval now
selects the matching `flow-sources-zonal-YEAR.yaml` automatically with `--year`.
To enable annual flow scoring set this template in `benchmark.yaml`:

```yaml
flows: data/hindcast/flows_physical_{year}.csv
flow_kind: physical
```

Missing configured flow inputs are errors. Leave `flows: null` for the price-only
work package. Unsupported directional API queries are not replaced with zero.

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

## Optional published-network reference comparison

Run from the root of your existing PyPSA-Eur checkout, using its Pixi environment:

```bash
pixi run python scripts/hindcast.py fetch-reference --year 2023
pixi run python scripts/hindcast.py benchmark --config config/hindcast/benchmark-reference.yaml
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
  --config config/hindcast/benchmark-reference.yaml \
  --network results/hindcast-resolved-2023/network.nc \
  --output results/hindcast-resolved-2023/benchmark
```

Gurobi is the default for the native build configurations (`gurobi-default`)
and for `hindcast.py solve`. The standalone helper uses barrier with crossover
disabled and accepts an explicit `--threads` value (default 1).
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

## Country-level reference flow comparison

Retrieve directional data with your existing ENTSO-E credentials in the environment
(`ENTSOE_API_KEY`). Do not commit the token. `entsoe-py` is already included in
this branch's Pixi dependencies and is needed for retrieval only.

```bash
pixi run python scripts/retrieve_hindcast_flows.py --year 2023 --kind physical \
  --config config/hindcast/flow-sources-reference.yaml \
  --output data/hindcast/flows_reference_2023.csv
```

Then set in `config/hindcast/benchmark-reference.yaml`:

```yaml
flows: data/hindcast/flows_reference_2023.csv
flow_kind: physical
```

Re-run the reference benchmark command. Border mappings in `flow-sources-reference.yaml` sum
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
**not a DK2 price**. The reference DK–DE spread is labelled accordingly.

The main zonal workflow uses the complete bus CSV above. For a separate, explicitly
limited diagnostic, the CLI also accepts inline regions as in this example:

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
price CSV headers before mapping. This separate diagnostic uses demand-weighted
multi-bus prices by default; it is not zonal market clearing. The main `bus_zones`
workflow always requires native copperplating and checks common prices.
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

`config/hindcast/build.yaml` and the annual `build-YEAR.yaml` files are the main starting point.
`build-2023.yaml` retains the earlier standalone example. It fixes
existing capacities, disables planned grid projects, uses hourly
2023 demand/weather and annual 2023 renewable capacity estimates, and clusters by
bidding zones. It uses 2020 technology assumptions (the source does not provide a 2023 cost
table), while historical 2023 fuel and CO2 inputs enter its native cost pathways.
It disables CO2 quantity caps. The prepared hourly files bypass the standard
price-source builders; generator cost calculation stays in PyPSA-Eur.

Build the **composed** network, not the ordinary solved target:

```bash
pixi run snakemake resources/hindcast-build-2023/networks/composed_2023.nc \
  --cores 4 --configfile config/hindcast/build-2023.yaml
```

First perform the same command with `--dry-run`. Full construction requires the
normal upstream datasets, weather cutout and any CDS credentials. The earlier configuration was checked against the source branch's JSON schema.
No schema check or complete data-building DAG has been run for this standalone
branch; construction remains pending.

This rebuild is NOT an exact reconstruction of the paper: current topology and
plant databases, 2023 conventional capacities, bidding-zone clustering, and ERA5
solar differ. Removing planned projects does not reconstruct commissioning dates
of every asset already present in today's grid database. Audit the historical
network before calling the new model validated. The paper used SARAH3 solar,
and held its conventional/hydro fleet constant from a 2020 database. The original
archive's embedded metadata can be exported with `hindcast.py inspect` to guide
reconciliation. The standalone workflow does not alter the Bornholm implementation.

The composed network already includes the historical hourly costs. Solve it
without a total marginal-cost override:

```bash
pixi run python scripts/hindcast.py solve \
  resources/hindcast-build-2023/networks/composed_2023.nc \
  --output results/hindcast-build-2023/network.nc --solver gurobi --threads 4
```

The native input handling described in section 2 also applies to this standalone
configuration. Its fuel CSV is in EUR/MWh thermal and its `co2` column is in
EUR/tCO2. VOM and other fuels use the native processed technology costs.
The optional `build_hindcast_costs.py` helper and `--marginal-costs` override are
reserved for workflows that explicitly supply total generator costs; they are
not needed for the native historical build.

The four traded-fuel and CO2 files follow section 1c, including its explicit
alignment limitations. Obtaining the authors' original CO2 series and resolving
the eight oil dates remains necessary before claiming exact source replication.

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
The real archived 2023 network was benchmarked with the six-country reference
configuration, now named `benchmark-reference.yaml`:

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

## Scope revision (2026-10-06)

The Europe-wide zonal configuration, complete bus mapping and all-border
reporting were added after the recorded implementation check above. No tests,
model solves or benchmark validation were run for this revision, as requested.
The earlier 11-test result does not describe the revised code.

## Native copperplate revision (2026-10-21)

Enabled native AC/DC cluster mixing and configured DE–LU and IE/NI copperplate
groups. Added exact native-ID aliases, structural guards for shared-market
mappings and representative-bus price extraction with a numerical agreement
check. No tests, schema checks, network builds, solves or benchmark validation
were run for this revision, as requested. Generated IDs and historical shape
boundaries still need review during construction; a completed historical price
reconstruction is not claimed.

## Multi-year revision (2026-10-21)

Added native annual scenarios for 2022–2025, a pinned Marta price importer,
year-expanded benchmark paths and a combined annual score table. No data
retrieval, tests, schema checks, network builds, solves or benchmark validation
were run for this revision. Completing these configuration/helpers does not
establish that all 2025 upstream model inputs are currently available.
