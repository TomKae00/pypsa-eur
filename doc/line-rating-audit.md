# Physical AC line ratings: audit and integration proposal

This work checks the physical ratings of the original transmission circuits.
It does not infer ratings from electricity prices or bidding-zone transfer
limits. Commissioning dates, outages, NTCs and flow-based market constraints
are outside this step. The tool never modifies a network.

## Flexible source selection

`scripts/audit_line_ratings.py` accepts a model year and one or more explicitly
selected source editions. Neither the algorithm nor the output directory is
tied to 2023. The source registry is
`config/hindcast/line-rating-sources.yaml`. A new edition is a new registry
entry, not a change to the matching algorithm.

Four JAO editions have been downloaded and their workbook layouts inspected:

| Source ID | Reference records | Approximate endpoints available | Important qualification |
|---|---:|---:|---|
| `jao-core-2022-09` | 2,190 | 0 | Separate TSO workbooks; one zero-byte member excluded; older Slovak duplicate excluded; no bundled map |
| `jao-core-2023-09` | 2,933 | 2,839 | Combined workbook; includes a malformed seasonal cell |
| `jao-core-2024-09` | 2,939 | 2,890 | Combined workbook; includes the same malformed seasonal cell |
| `jao-core-2025-09` | 2,975 | 2,928 | Publication workbook selected explicitly; `_versions_` drafts excluded |

Counts are source records, not unique physical circuits. The same cross-border
circuit can be reported by both TSOs. The year in an edition name does not
establish full-year applicability. In particular, the 2022 archive includes a
Belgian workbook labelled Q4 2021 and no usable TransnetBW workbook in this
download. These gaps are not filled from a later edition automatically.

The 2023 and 2024 workbook dates in the registry come from workbook filenames;
they are edition labels, not an independently verified asset-validity date.
URLs, exact archive hashes, selected/excluded members and input-network hashes
are retained in `audit.json`. Changing the source config also changes its
recorded hash. Source choice is always explicit; there is no nearest-year or
latest-edition fallback.

## Run the audit

Use an original-voltage `base.nc` or `base_extended.nc`, before simplification
to 380 kV. The normal PyPSA-Eur environment supplies pandas, numpy, PyYAML,
netCDF4 and the Excel reader. No solver or PyPSA optimization is used.

Example for a 2024 network and an explicitly selected 2024 edition:

```bash
python scripts/audit_line_ratings.py \
  resources/hindcast-build-2024/networks/base_extended.nc \
  --year 2024 \
  --source-id jao-core-2024-09 \
  --output results/line-ratings/2024
```

Add `--offline` to require cached inputs. Downloads are cached under
`data/hindcast/line-ratings/sources/`. The cache retains source edition IDs and
downloaded bytes are checked against the registry hashes.

To compare source editions against the **same** physical network, repeat
`--source-id`. This is a source-comparison exercise, not a hindcast of the later
years. The `--year` flag records the target model year; it does not shift
network timestamps, alter topology or claim that an edition applies to that
whole year.

To trace which original lines remain in the zonal network, also supply:

```bash
  --clustered-network results/hindcast-corrected-2024/network.nc \
  --linemap resources/hindcast-build-2024/linemap_cluster_network.csv
```

These two inputs must come from the same build as the original network. The
audit checks nominal sums and reports mismatches. A matching sum establishes
capacity arithmetic only, not electrical equivalence or correct market limits.

## Additional TSO sources

A second adapter accepts a normalized CSV. Add an explicit edition entry:

```yaml
sources:
  tso-reviewed-edition:
    format: normalized_csv
    path: data/hindcast/line-ratings/tso-reviewed-edition.csv
    provenance_url: https://the-tsos-actual-publication-url
    snapshot_date: ''
```

Replace the example URL with the real source. Required CSV columns are
`reference_id`, `name`, `voltage_kv`, `lon0`, `lat0`, `lon1`, `lat1`, and
`fixed_current_a`. A missing rating or coordinate is blank, never zero.
Optional columns are `tso`, `eic`, `station0`, `station1`,
`reference_length_km`, `period_1_a` through `period_6_a`, `dlr_min_a`,
`dlr_max_a`, `map_geometry_count`, and `source_workbook`.

Coordinates are longitude/latitude in WGS84; currents are amperes; lengths are
kilometres. Use a unique reference ID within each edition. The adapter
preserves periodic and DLR fields separately. It does not average them or
choose a seasonal value on the user's behalf. Invalid numeric source values
are retained in `source_parse_issues`, and the row is marked
`source_data_issue`.

For a source with no georeferencing, retain the source inventory and add
reviewed endpoint coordinates through a normalized CSV, with provenance for
those coordinates. Do not silently borrow a later year's line identities.

## What the output files mean

* `original_line_inventory.csv`: original OSM IDs, circuit counts, voltage,
  nominal capacity, assumed current per circuit, endpoints, type and source tags.
* `reference_ratings.csv`: source records and the original rating categories.
* `candidate_matches.csv`: geographically plausible pairs at similar voltage.
* `rating_review.csv`: reference details beside each candidate, with blank
  reviewer, validity, allocation and approved-capacity fields. All `approved`
  values start false.
* `interzone_rating_review.csv`: candidates whose original line is mapped to
  an AC line retained in the supplied zonal network. This is the main starting
  point for the present price study.
* `country_coverage.csv`: candidate coverage for all model countries. A border
  line is counted as incident to both countries; these counts must not be
  summed into a European total. Candidate coverage is not validated coverage.
* `reference_coverage.csv`: source-record counts by edition, TSO and rating kind.
* `aggregation_check.csv`: nominal-capacity sums, when lineage inputs are supplied.
* `audit.json`: selections, hashes, counts and qualifications.

Candidate matching checks both endpoint orientations with a configurable
distance tolerance (default 5 km). Voltage labels may differ by at most 6%,
allowing 220/225 and 380/400 candidates, but differing labels are flagged.
The JAO map is schematic: some parallel circuits have displaced endpoints.
Proximity does not prove circuit identity. Circuit names, EICs, OSM relation
IDs and circuit references still need review.

`priority_for_review` means a uniquely located, single-circuit, fixed-rating
candidate with matching voltage, both endpoints within 1 km, route length
within 0.8–1.25 of the reference and one map geometry. It is a transparent
shortlist, not an approval or statistical confidence estimate. In particular,
JAO tie-line lengths can cover only the reporting TSO's share up to the border;
their total-route length ratio can legitimately exceed this range.

## Conditions for a physical correction

Before approving a row, establish:

1. The physical circuit identity and voltage, including parallel circuits,
   branches, shared routes and source records for the two sides of a border.
2. Whether the current is a continuous, seasonal or temporary rating, and
   whether it describes the conductor, a limiting terminal or a whole circuit.
   A temporary emergency rating must not become an unrestricted annual limit.
3. The equipment and dates to which the value applies. `valid_from` and
   `valid_to` are decisions supported by evidence, not inferred from the filename.
4. How the rating maps to the model record. A source circuit split into several
   serial model segments must not have its rating divided among those segments.
   Two TSO records for the same circuit must not be added as parallel capacity.
   Multiple circuits within one model record require explicit reconciliation.
5. How different limiting ratings along a circuit are handled. Do not average
   unequal terminal or segment limits; identify the binding equipment and rating
   regime. Leave unresolved conflicts unapproved.

For a supported constant circuit current, the nominal apparent rating is
`sqrt(3) * voltage_kv * current_a / 1000` MVA, multiplied by the reconciled
number of parallel circuits. Keep the physical nominal rating separate from
`s_max_pu` and its security approximation. Do not encode a market NTC by
reducing physical conductor ratings.

## Reviewed bindings and native integration

The first circuit groups have now been reviewed and are stored in
`config/hindcast/line-rating-bindings.yaml`. This file binds original OSM
records to circuit EICs and documents the identity evidence. It contains no
fixed year or hard-coded current value.

The preparer `scripts/prepare_line_rating_corrections.py` combines those
identities with one explicitly selected source edition. It accepts a group only
when every EIC is reported by both TSOs, every report uses a positive fixed
current, all reported currents agree, and the number of source circuits equals
the original model's circuit count. Conflicts, incomplete data and seasonal/DLR
records are held. It never averages conflicting ratings or assumes the lower
of two values is necessarily a comparable continuous limit.

For Krško–Tumbri, the association of the two EICs with the two OSM records is
not resolved individually. Both source circuits have the same agreed current,
so applying that current to the two verified single-circuit records is
permutation invariant. The group is rejected if that equivalence changes.

The present policy is **fixed snapshot adopted as a representative annual
physical input**. It does not claim verified unchanged ampacity for every hour
of the year. This choice excludes commissioning dates and outages, as agreed
for this work package. Source edition, model/planning year, and weather year
remain distinct. A source snapshot may be selected for another model year
explicitly; there is no automatic nearest-year fallback.

For example:

```bash
python scripts/prepare_line_rating_corrections.py \
  resources/hindcast-build-2024/networks/base_extended.nc \
  --year 2024 \
  --source-id jao-core-2024-09 \
  --rating-policy fixed-snapshot-representative \
  --output data/hindcast/line-ratings/corrections_2024.yaml
```

The generated file records model year, source URL/hash, bindings hash, original
network hash, the accepted and held groups, and per-line expected signatures.
Its native correction payload is `Line: {s_nom: {index: {original_osm_id: MVA}}}`.
The network hash is provenance; runtime compatibility uses the selected lines'
original endpoints, voltage, circuit count, type and default nominal capacity,
not a byte-for-byte comparison between base and base-extended NetCDF files.

Enable the optional overlay `config/hindcast/line-ratings.yaml` alongside the
annual build config when constructing a **separate comparison run**. Its option
is:

```yaml
lines:
  rating_corrections: data/hindcast/line-ratings/corrections_{year}.yaml
```

`{year}` comes from the single planning horizon, not the weather snapshots.
The default is null, so existing cases retain their baseline behaviour.
Multiple planning horizons in one build are rejected for this input; use
independent annual hindcast builds. The optional correction file is a declared
Snakemake input, so a changed payload triggers the affected build stages.

The base-network hook checks the compiled policy, model year, provenance,
original IDs and signatures, and positive finite ratings **before** updating
anything. It restricts this file to nominal AC ratings, then reuses the existing
`_apply_parameter_corrections` function. The old GridKit correction file remains
separate. Missing or changed original lines fail explicitly instead of being
silently ignored. Fixed source corrections and active DLR cannot be combined
because DLR would replace the static loading limit.

The hook acts after standard `s_nom` calculation and before simplification.
Physical-rating provenance is preserved through simplification, clustering and
composition in `network.meta["physical_line_ratings"]`. No new optimization
constraints or alternative aggregation routines are introduced.

Correcting `s_nom` changes the limit, not the assumed impedance or conductor
type. Source 380-kV and model 400-kV labels are retained separately; the agreed
current is converted using the model nominal voltage. It is not presented as
the TSO's actual operating MW limit. During native simplification, the equivalent
parallel-circuit count changes to preserve the corrected nominal capacity.

The prepared examples use the supplied original network signatures. A rebuild
against another grid release or with different line types must regenerate and,
where necessary, re-review the bindings. Existing weather cutouts can be reused.
The original source audits remain candidate inventories; the compiled correction
files and their decision records identify the accepted subset.

## First reviewed subset

| Circuit group | Accepted fixed current | Edition behaviour |
|---|---:|---|
| Maribor–Kainachtal, two circuits | 1,920 A per circuit | Both TSOs agree in all four selected editions |
| Krško–Tumbri, two circuits | 1,920 A per circuit | Both TSOs agree in all four editions; equal-current group mapping |
| Divača–Pehlin, one circuit | 960 A | Both TSOs agree in all four editions |
| Sokolnice–Stupava, one circuit | 2,000 A | Accepted in 2023–2025; held in 2022 (1,530 versus 2,000 A) |
| Sokolnice–Križovany, one circuit | 2,000 A | Accepted in 2022; held in 2023–2024 for disagreement and in 2025 for incomplete reporting |

Each selected edition therefore produces five original model-record updates,
representing six physical circuits. This is a small supported subset, not a
claim that Europe has been fully corrected. Nordic circuits have not received
new physical-rating overrides in this subset. Market transfer limits, stability
constraints, and internal constraints removed by copperplating remain separate.

## Sources

* JAO Core Static Grid Model editions, workbooks, maps and handbooks:
  https://www.jao.eu/static-grid-model
* PyPSA-Eur line-type defaults:
  https://pypsa-eur.readthedocs.io/en/latest/configuration/
* Xiong et al. (2025), grid construction, validation and data limitations:
  https://doi.org/10.1038/s41597-025-04550-7
* 50Hertz static grid model, a possible additional German reference:
  https://www.50hertz.com/Transparency/GridData/Gridfigures/Staticgridmodel/
* NESO ETYS Appendix B, a candidate source for a separate British adapter;
  it has not yet been normalized or matched in this package:
  https://www.neso.energy/publications/electricity-ten-year-statement-etys/etys-documents-and-appendices

Nordic, Iberian, Italian, British, Irish and other countries without adequate
Core coverage still need additional circuit-rating sources. Existing Nordic
market-capacity publications are not substitutes for circuit ampacities.
