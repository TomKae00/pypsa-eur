# Bornholm within the native PyPSA-Eur workflow

This supersedes the proxy-profile/manual-allocation instructions in
[the first input patch](bornholm_input_patch.md) for the supplied Bornholm configs.
The integration is for **overnight** scenarios. Inspection remains a separate
command; no analysis or financial-accounting stage is added to the workflow.

## What is reused

- `determine_availability_matrix.py` and atlite's existing fine-resolution
  exclusion container remove the project footprint from **all offshore wind
  categories**, before regional potential and profiles are calculated. The
  exclusion is intersected with the ordinary siting restrictions, rather than
  subtracting 3000 MW or multiplying coarse weather-cell area fractions.
- Two project rules call the **same availability and renewable-profile scripts**.
  They use the native cutout, turbine, correction factor, CF-weighted layout and
  clipping settings of `energy_islands.profile_technology`, with one resource
  class per project. There is no second wind conversion algorithm.
- Hubs and project wind are attached after clustering and **before native time
  aggregation**. Export links are attached after grid preparation, preserving
  their scenario capacities and physical route lengths.
- Native prepared, annualised technology costs are used. A shared connection
  cost helper preserves the existing generic-offshore formula.

The project profile uses its planning area as an exogenous site. Generic depth,
shore-distance, Natura and shipping eligibility screens are **not reapplied to
this site**. Its whole footprint is reserved, and turbine positions are represented
by the native CF-weighted layout. This is a planning-envelope approximation, not
an engineering layout or a claim that every square kilometre is buildable.
The project's specified 3000 MW is independent of the generic 2 MW/km² density;
profile generation does not impose that generic density as a project capacity cap.

## Geography and provenance

`data/energy_islands/bornholm_area.geojson` dissolves the three polygons from the
Danish Energy Agency's `Havvindmoelleparkomraade_010224` dataset, reprojected from
ETRS89 / UTM zone 33N (EPSG:25833) to WGS84 (EPSG:4326), without simplification.
The feature name `BEI_WIND` matches the wind CSV's `generator_id`.

- [Official GIS archive](https://ens.dk/media/6117/download)
- [Adopted Plan for Program Energiø Bornholm, 13 May 2024](https://ens.dk/media/2101/download/), section 3.1
- [Agency plan and data page](https://ens.dk/energikilder/plan-og-miljoevurderingsproces-energioe-bornholm/plan-program-energioe-bornholm)
- Downloaded 23 October 2026. Archive SHA256:
  `6d56a9070275b193312d79b8b090bd2a7bdf403762de1c94e43c4354b87143fd`.

The GIS polygon areas in their original projected CRS are approximately 123.15,
117.93 and 407.18 km², consistent with the plan's rounded 123, 118 and 408 km².
This is the **2024 planning envelope**, not a certified final tender boundary or
final turbine layout. Updating that input requires a documented boundary review.
The weather cutout must cover the complete area; otherwise the project rule fails.
Coverage is checked against atlite's native rectangular raster bounds, including
the outer half-cells. Unioning individual cell polygons can leave floating-point
slivers on 0.3-degree grids and falsely report missing coverage. Genuine partial
or complete non-coverage still raises an error, with both bounds included.

## Scenario and cost accounting

| Item | Baseline | Radial | Hybrid |
|---|---|---|---|
| Planning area reserved from generic offshore development | Yes | Yes | Yes |
| Explicit project wind | None | 3000 MW | 3000 MW |
| Export to Danish receiving bus | None | 3000 MW | 1200 MW |
| Export to German receiving bus | None | None | 2000 MW |

This baseline means **reserved site, project not built**. Allowing generic wind to
use the site in the baseline would define a different counterfactual. Keep
`exclude_project_area` identical across the primary comparison. It is independent
of `enabled`, which controls only construction of explicit project assets.
`subtract_generic_potential: false` prevents additional MW subtraction; attempting
to combine geographic exclusion with manual allocation raises an error.
Native existing/minimum-capacity handling is unchanged.

Project wind carries turbine cost plus `offwind-ac-station` and the native AC
submarine connection cost over the CF-weighted weather-cell-to-hub distance.
Collection landfall length is zero, and collection length factor is 1.0. The
station term is the native AC station proxy, not a bottom-up Bornholm substation
estimate. Each export Link carries its HVDC cable cost and **one inverter-pair
cost**. No DC wind-to-mainland connection package is also added to the generator.
The hub itself has no additional equipment cost. Review these approximations
before treating project costs as engineering estimates.

`fixed` assets use equal extendable lower/upper bounds, so mandatory investment
enters the optimisation objective. For finite native transmission budgets, the
mandatory export-link cost or volume is added to the budget RHS: forced Bornholm
assets do not use up the ordinary grid expansion allowance. Optional expansion
still competes for that allowance. Project links remain ordinary `DC` carriers.

No duplicate Bornholm-labelled project was found in the committed TYNDP2020, NEP
or manual project tables. These native transmission inputs remain enabled; no
second project is inserted into those tables. Existing-ID collisions are errors.
That is not a geographic deduplication guarantee for future OSM/TYNDP releases:
review their Bornholm assets before combining a new release with this addition.
Do not remove the existing Bornholm–Sweden supply connection as if it were a new
energy-island export project.

## Run and validate

The missing smoke CSVs are now included. They use the buses observed in the user's
20-cluster DK/DE/SE/PL baseline: `DK0 0AC` and `DE1 1AC`. These are explicit nodal
bindings, not a bidding-zone price model. Recheck them if clustering or network
inputs change. The regular 90-cluster links CSV deliberately has **blank receiving
buses**: compose its baseline, inspect the mapping and fill those IDs independently.
The regular project case will fail clearly until that mapping is supplied.

From the repository root:

```bash
pixi run python -m pytest --noconftest -q test/test_energy_islands.py test/test_energy_island_resources.py
pixi run snakemake --cores 4 --configfile config/config.bornholm.smoke.yaml -- resources/smoke_baseline/networks/composed_2035.nc resources/smoke_radial/networks/composed_2035.nc resources/smoke_hybrid/networks/composed_2035.nc
pixi run python scripts/inspect_energy_islands.py resources/smoke_baseline/networks/composed_2035.nc --links data/energy_islands/bornholm_links.smoke.csv --wind data/energy_islands/bornholm_wind.smoke.csv
pixi run snakemake --cores 4 --configfile config/config.bornholm.smoke.yaml
```

The inventory's generic donor/profile suggestions are no longer needed for the
project-profile mode. It remains useful for receiving-bus review. Cable lengths
are still hub-to-landing distance × 1.25 estimates unless `length_km` is supplied;
explicit physical lengths are not multiplied again. The supplied 2035 commissioning
year and hub/landing coordinates are scenario assumptions retained from the branch.

New exclusion inputs and project profile dependencies are tracked by Snakemake;
rebuild affected availability matrices, profiles, composed and solved networks.
Do not reuse pre-exclusion solved networks for comparison. With custom shared
resource settings, keep exclusion settings identical for runs that share offshore
resource files, or give them separate resource paths.

Regression tests run the actual atlite and native profile code with synthetic
weather, exercise native temporal averaging, and solve tiny radial/hybrid cases.
They verify cost and capacity accounting; they do not replace a full European
smoke run with the actual cutouts, clustered buses and cost data.

### Verification of this change

On 24 October 2026, all **50 tests** in `test_energy_islands.py`,
`test_energy_island_resources.py` and the existing `test_compose_network.py`
passed with PyPSA 1.2.4, atlite 0.6.1 and linopy 0.8.0. The integration cases
export and reload the composed network before solving, matching the native
compose/solve rule boundary. Both transmission-volume and transmission-cost
budgets are tested, including objective accounting for mandatory project assets.
The real workflow's input/parameter expansion was also checked for all three
smoke scenarios, all three offshore resource categories and onshore wind.
Python lint/format, Snakemake formatting and diff whitespace checks passed.
The coverage correction adds four regressions using the committed Bornholm
polygon: successful availability calculation on 0.25- and 0.3-degree grids,
and rejection of partially and completely uncovered sites. The 0.3-degree
regression reproduces the reported error before the correction.

The execution environment required an external psutil process-ID workaround;
that workaround is not part of the repository. A complete dry-run did not finish
while checking remote source-data availability, so it is not reported as passed.
No full European build or solve was performed here.

The regular configuration still uses 90 clusters and 3-hour averaging, and the
smoke configuration uses 20 clusters and 6-hour averaging. Export links remain
lossless. This change completes the resource/profile integration; it does not
calibrate the future electricity system or introduce a bidding-zone market model.
