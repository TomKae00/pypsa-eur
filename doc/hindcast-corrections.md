<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Geographic and reservoir preprocessing corrections

## Scope and native integration

These corrections address two preprocessing failures found when inspecting the
fixed-capacity European 2023 hindcast. They use the existing PyPSA-Eur inputs,
spatial aggregation and hydro representation. No regional price observations
enter either correction.

| Stage | Existing capability | Correction |
| --- | --- | --- |
| Zone definitions | `build_bidding_zones` supplies native zone polygons | `base_network.build_admin_shapes` consumes these polygons for `bz`, including country overrides, instead of dissolving whole NUTS3 regions labelled by majority overlap |
| Grid aggregation | Native bus maps and PyPSA line aggregation | Rebuild the maps using the corrected geometry; retain line-rating sums and the configured `s_max_pu` |
| Plant placement | Native country-constrained spatial assignment and coastal nearest-region fallback | Consume the rebuilt regions; no named-plant reassignment table |
| Reservoir estimation | National storage estimates and plant-level reported durations | Complete missing durations before capacity-weighted aggregation; remove the subsequent first-entry overwrite |
| Reservoir operation | Native `StorageUnit` inflows, efficiency, dispatch limits and cyclic state of charge | Retained |
| Workflow | `compose_network` already lists `add_electricity.py` as a code dependency and supplies `hydro_capacities.csv` | Reuse that dependency and input; no new retrieval or composition rule |

The changes apply generically to the countries and years selected by the build
configuration. They do not tune Swedish transfer limits, Norwegian reservoir
volumes, generator costs or availability to match observed prices.

## Evidence from the original 2023 case

The supplied solved network records revision
`86dd9ad9fb5a24a2e0b31e78ad285b6d3b2415ed`, PyPSA 1.2.4, and 8,760 hourly
snapshots. The investigation used the corresponding `simplified.nc`, cluster
bus/line maps, individual `powerplants.csv`, administrative shapes and native
bidding-zone polygons. These generated files are not committed here.

### Geography

The original build mapped Ringhals 3/4 and Oskarshamn 3 to SE4. The native
administrative-region construction dissolved NUTS3 regions after assigning each
whole region to its largest-overlap bidding zone. Consequently, boundaries
inside a NUTS3 region were lost for grid and asset placement.

Reassigning the *same saved simplified network* using its native zone polygons
changes the SE3–SE4 AC corridor as follows:

| Quantity | Original mapping | Native polygon reassignment |
| --- | ---: | ---: |
| AC member lines | 9 | 5 |
| Sum of nominal ratings | 16.087 GW | 8.937 GW |
| Operating limit, `s_max_pu = 0.7` | 11.261 GW | 6.256 GW |

Eight old border lines become internal to SE3 and four other lines become border
lines. The 260 MW mainland–Gotland link becomes internal to SE3; the 1,200 MW
South-West Link remains an SE3–SE4 connection. Fourteen of the 134 Swedish buses
in that simplified network change zone across all Swedish boundaries.

These figures isolate reassignment of an existing network. They are not results
of a full rebuild, and do not establish historical market transfer capacities.
A rebuild can also change simplification, regions, demand and asset allocation.
The existing handling of areas without grid buses remains a separate limitation
to review for islands and other small regions.

### Clipped cells at coastlines and borders

The first corrected 2023 rebuild failed before solving, during demand
disaggregation. The supplied diagnostic identified 17 missing/empty geometries
in the base onshore regions and two surviving after simplification:
`virtual_way/691457479:b:1-300` (NO) and `way/693745438-220` (NL).
All 17 reported coordinates lie outside their country's native zone outlines
in the supplied polygon dataset. A local reproducer confirms that a bus mapped
to its nearest administrative region can have a Voronoi cell with no
intersection with that region. Native clipping then removes the cell, while
subsequent reindexing creates a row with missing geometry.

`process_onshore_regions` now exports only cells with positive intersection
area. It logs the affected bus IDs and retains all electrical buses and their
asset assignments. It does not expand the zone polygon or move a bus to another
country/zone. A bus without a land cell receives no area-based demand or
renewable potential directly; it remains part of the grid and subsequent
network aggregation. This handles the selected geometries consistently; it
does not establish that the underlying coastline or bus coordinates are exact.

Demand allocation now filters existing cells to eligible substations, rather
than reindexing and recreating missing shapes. A country with national demand
and grid buses must still have eligible land cells. Remaining malformed input
geometries are rejected with their IDs.

The existing JRC Energy Atlas, GB consumption and NUTS GDP/population keys keep
their priority. After combining them, relative weights are normalised once per
country to preserve the supplied national load at every hour. This generalises
the existing intended GB/NI renormalisation; the old GB membership test checked
the Series index, and mixing individually normalised source keys can otherwise
inflate national demand. Missing, negative, non-finite or zero-total weights
stop the build rather than triggering an arbitrary uniform allocation. Logs
record the combined weight sum and number of eligible cells per country.

Focused regressions check coverage of the original outline, retention of
electrical buses, singleton cells, substation selection, malformed geometries,
missing country coverage and unusable distribution data. Integration cases run
the actual allocation functions against small GeoTIFF, GB Excel/geometry and
NUTS fixtures, including a landless bus and GB/NI source mixing, and verify
hourly country-total conservation. Full corrected HPC rebuilds remain pending.

### Reservoirs

The Norwegian input contains 918 reservoir plants, 491 with missing durations,
representing 5,896.39 MW of discharge capacity. Reported durations account for
44.900471 TWh against the native national estimate of 84.4 TWh.

Two operations caused the failure:

1. Capacity-weighted aggregation summed over missing durations, leaving a
   positive regional duration even when some constituent plants had no duration.
   Subsequent missing-value detection could no longer identify those plants.
2. `~hydro.country.duplicated()` selected the first entry of every country,
   including countries with multiple reservoir entries. NO3 was selected and
   its existing 2.445429 TWh was overwritten by the entire 39.499529 TWh residual.

Replaying this arithmetic on the supplied individual plant file reproduces the
five saved Norwegian durations to floating-point precision. The resulting
national total is 81.954571 TWh. Related overwrite errors affect other countries,
including countries with a single aggregated reservoir entry.

Completing missing plant durations first, using the same national residual
heuristic and **the original plant-to-zone assignments**, gives:

| Zone | Original storage (TWh) | Completion before aggregation (TWh) |
| --- | ---: | ---: |
| NO1 | 4.441805 | 14.637181 |
| NO2 | 11.632508 | 17.951345 |
| NO3 | 39.499529 | 7.247761 |
| NO4 | 16.369579 | 19.208922 |
| NO5 | 10.011150 | 25.354791 |
| Total | 81.954571 | 84.400000 |

These are estimated allocations, not independent observations of regional
reservoir capacities. Their distribution may change again after geographic
reassignment. The geography and hydro corrections therefore need separate
comparisons before attributing changes in prices to either one.

## Reservoir completion convention

For `energy_capacity_totals_by_country`, let `K` contain plants with positive
reported durations and `M` contain plants with missing or zero durations.
The existing national target retains its 0.2 TWh floor where an estimate exists:

```text
E_reported[c] = sum(p_nom[i] * max_hours[i] for i in K[c])
E_remaining[c] = max(E_target[c] - E_reported[c], 0)
h_missing[c] = E_remaining[c] / sum(p_nom[i] for i in M[c])
```

Only missing durations receive `h_missing`. Applying the same duration to these
plants allocates the remaining energy proportionally to their discharge
capacities. Native capacity-weighted aggregation then preserves
`sum(p_nom * max_hours)` for each region.

The implementation explicitly handles these cases:

- Positive reported durations remain unchanged, including in single-plant
  countries. A complete dataset is not rescaled to a national estimate.
- If reported energy already exceeds the estimate, it is retained, missing
  plants receive zero residual storage, and the build logs a warning. Zero
  estimated duration is not reinterpreted as missing during component creation.
- If a positive residual has no positive-capacity missing plants to receive it,
  it remains unallocated with a warning. No arbitrary first plant is selected.
- Where a country estimate is unavailable or cannot yield a finite duration,
  missing durations use the existing six-hour fallback, with a warning.
- `estimate_by_large_installations` retains its existing national energy-to-power
  ratio for missing plants. It is now also applied before aggregation.
- Pumped storage continues to use `PHS_max_hours`; run-of-river is unaffected by
  reservoir-duration completion. Natural inflow allocation remains the native
  country-profile/capacity-share approach.

Build logs identify the heuristic, number of completed plant records, and
reported versus completed storage totals by country. These conventions preserve
the intended role of national estimates as missing-data heuristics; they do not
make them mandatory totals when reported data conflict or are already complete.

## Rebuild and evaluation

Preserve the original solved network, audit, submission revision and mappings.
The geometry-only implementation is available at revision
`a73c75bfb867b8bafb88318dd505e84f34782ebe`; the reservoir correction follows it.
Use separate checkouts/output locations for controlled comparisons. Do not
change a checkout underneath queued or running network builds.

The geometry correction requires a rebuild from the spatial stages. The hydro
correction requires recomposition from individual plants, followed by solving;
editing a solved network or merely re-solving the old composed file is not
sufficient. Annual ERA5 cutouts and unchanged retrieved datasets can be reused.

Compare baseline, geography-only, and geography-plus-hydro results using the
same observations, snapshot resolution, costs, solver settings and score
definitions. Report mean bias, MAE, correlation and price variability together.
A successful solve alone is not evidence of historical-price accuracy.

Focused regression cases cover boundaries cutting through NUTS3 regions,
per-country `bz` overrides, unaffected national clustering, mixed known/missing
durations, row-order invariance, preservation of energy through the actual plant
aggregation path, complete/conflicting country data, zero capacities and the
existing alternative/fallback heuristics. Component creation is checked to
preserve a completed zero duration and the existing inflow allocation.

The following focused selection passed all 14 cases with PyPSA 1.2.4,
pandas 3.0.5 and powerplantmatching 0.8.0:

```bash
python -m pytest -q test/test_add_electricity.py test/test_base_network.py \
  -k 'hydro or reservoir or bidding_zone_shapes or attach_load'
```

Ruff lint and formatting checks passed for the changed Python files. These are
preprocessing checks; a full rebuild and solved-price effects remain to be
evaluated on the HPC.

## Source references

- [Native reservoir implementation in PyPSA-Eur v2026.08.0](https://github.com/PyPSA/pypsa-eur/blob/v2026.08.0/scripts/add_electricity.py)
- [Native PyPSA 1.2.4 spatial aggregation](https://github.com/PyPSA/PyPSA/blob/v1.2.4/pypsa/clustering/spatial.py)
- [Zone polygon dataset version used by the build](https://github.com/electricitymaps/electricitymaps-contrib/blob/v1.238.0/web/geo/world.geojson)
- [Svenska kraftnät network development plan, including Gotland in SE3](https://www.svk.se/49011f/contentassets/d92ea45b663f4e1f839d025479182f13/svk_natutveckling_nup_2026-2035_eng.pdf)
