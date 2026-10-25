# Bornholm: common planned electricity grid for 2035

## Definition and scope

The main comparison uses the native OSM 0.7 physical grid plus a pinned selection
of TYNDP 2026 investments expected by 2035. It is a **planned-project scenario**,
not a guarantee that every project will be completed, and not the official
ENTSO-E 2035 NTC reference grid. Native bidding-zone aggregation retains AC
network physics. No zonal transfer-capacity table is substituted for physical
cable ratings or AC line parameters.

The input snapshot comes from
[p-glaum/tyndp2026-to-pypsa](https://github.com/p-glaum/tyndp2026-to-pypsa/tree/7c207622c12c47447907d67a3a4de4ae53d633b4),
commit `7c207622c12c47447907d67a3a4de4ae53d633b4`. Its converted data are reused;
its example electricity-only configuration and converter implementation are not.
The original data are ENTSO-E's
[TYNDP 2026 project portfolio](https://www.entsoe.eu/news/2026/04/02/tyndp-2026-project-portfolio-expands-further-to-199-transmission-and-69-storage-project/)
and project map, with technical attributes inherited by investment ID from
earlier TYNDP workbooks where the 2026 workbook omits them. This is not a claim
that all electrical parameters were freshly validated in 2026.

The [ENTSO-E 2026 grid methodology](https://2026.entsos-tyndp-scenarios.eu/scenarios-methodology-report/grid/)
uses additional maturity evidence for its reference grid. Our status filter
does not independently reproduce that project-by-project national-plan review.
In particular, `planned_not_yet_permitted` is an assumed delivery pipeline, not
a synonym for a final investment decision. TYNDP editions are biennial; there
is no separate electricity-project TYNDP 2025 edition to switch to.

## Selection and reproducibility

`scripts/prepare_bornholm_grid.py` is a small, offline selection script. Run it
against a clean checkout of the pinned converter commit:

```bash
pixi run python scripts/prepare_bornholm_grid.py /path/to/tyndp2026-to-pypsa
```

The committed `data/transmission_projects/tyndp2026/` tables include only:

- commissioning year at or before 2035;
- `under_construction`, `in_permitting` or `planned_not_yet_permitted` status;
- valid geometry endpoints/lengths and positive finite DC ratings;
- projects other than Bornholm Energy Island (project ID 1106).

The converter supplies 59 new AC entries, 19 new DC entries and 22 AC upgrades.
One additional DC entry recovers **Celtic, investment 810 / project 107**. Its
2026 converter audit gives 700 MW, 575 km, 2027 and `under_construction`, but no
route. The endpoints are copied from the existing native
`tyndp2020/new_links.csv` row `TYNDP2020_17` (Knockraha–La Martyre).
Only the endpoint geometry is reused; the older year and length are not.
The source identifier and geometry fallback are recorded in the cable tags.
The technical review below moves Harmony from DC to AC and excludes two other
rows. Final input: **99 branch entries** (60 new AC, 17 new DC, 22 AC upgrades)
before native mapping/deduplication, not 99 distinct interconnectors or guaranteed
additions to the final grid.

The `audit/` subdirectory records:

- `manifest.json`: upstream commit, input SHA256 hashes, filters and counts;
- `selection.csv`: the selection decision for each converted entry;
- `unconverted_2035.csv`: all 100 eligible entries omitted by the converter,
  including substations/PSTs, missing routes and short/unsupported branches.
  Celtic is marked recovered; the other rows remain explicit coverage limits.

Do not infer that omitted projects are absent from the native base grid. Conversely,
do not assume they are already present. Of the eligible omitted rows, 27 lack
route geometry, including Celtic before recovery. Examples still needing coverage
review include Nautilus and several national reinforcements. These are not
silently converted using guessed locations. The data therefore describe the
**covered planned-project subset**, not a complete reconstruction of all TYNDP
investments. Substations, transformers and PSTs require different treatment from
ordinary branch additions.

The converter assumes two parallel circuits for AC entries and derives ampacity
from standard line types; geometry and some electrical attributes are proxies.
These approximations matter when interpreting optimal Bornholm cable ratings.
Native OSM assets marked under construction remain included under the existing
`under_construction: keep` settings, assuming completion by 2035 unless a supplied
build year says otherwise. OSM itself does not provide a verified 2035 schedule.

### Reviewed corrections and exclusions

**Harmony Link (170/1034):** the converter inherits the former offshore HVDC
design (330 km, Darbenai–Zarnowiec) from 2024, despite the 2026 description and
map referring to the new onshore AC project. The [2026 Litgrid announcement](https://www.litgrid.eu/en/news/litgrid-seeks-eu-funding-to-strengthen-the-resilience-of-harmony-link-interconnection)
describes 700 MW at 220 kV. The [connection-point agreement](https://litgrid.eu/en/news/litgrid-and-pse-signed-connection-point-agreement-for-the-harmony-link-lithuania-poland-electricity-interconnector)
describes two AC circuits between Elk Bis and Gizai via Norki and Wigry, plus
refurbishment of the existing Ostroleka–Elk corridor.

The study therefore replaces the obsolete DC entry with a **two-circuit 220 kV
AC proxy**. Its start coordinate is the pinned OSM 0.7 Elk Bis station
(`way/303870256-400`); its remaining route uses the 2026 map vertices towards
Gizai. Length is the geodesic length of that partial route, recorded in the CSV.
The old Ostroleka–Elk route is not added again as a parallel new connection.
Refurbishment there and the new intermediate substations are not modelled
explicitly. Gizai maps to an existing Lithuanian bus before zonal aggregation.
The native 220 kV standard line type supplies approximate resistance, reactance
and thermal rating, even though the planned connection includes underground
cable. Neither that rating nor the general 0.7 loading factor is a calibrated
700 MW cross-border NTC. This approximation needs a grid sensitivity check
before treating optimal Bornholm capacities as engineering recommendations.

Two known incompatible DC rows are explicitly excluded rather than repaired
with guessed endpoints or silently accepted if the base network later changes:

- **SACOI3 (299/1458):** [Terna](https://www.terna.it/en/projects/projects-common-interest/sardinia-corsica-mainland-italy-electrical-interconnection)
  describes replacement of SACOI2 with a three-terminal 400 MW system. An
  independent additional two-terminal 400 MW branch would misrepresent it and
  could double-count the existing link. Existing native infrastructure remains.
- **GRITA2 (1112/1799):** the converter's landfall coordinates fail the native
  endpoint-distance check. [IPTO](https://www.admie.gr/en/erga/erga-diasyndeseis/GRITA2)
  identifies a new Thesprotia substation, rather than the existing GRITA1
  terminal. IPTO's project page gives 1 GW for the new link, consistent with the
  converter's rating. The new endpoint and associated onshore works need a
  separate review; it is not silently connected to the GRITA1 terminal.

These exclusions limit the common background; they are not assumptions that the
projects will be cancelled. They are recorded in the manifest and selection
audit. All Bornholm cases share them. A complete official 2035 grid would require
further project engineering and coverage work beyond this integration.

## Native integration and audit

The configuration registers `include.tyndp2026` in the validated schema; merely
adding an unknown YAML key would not reliably enable it. TYNDP2020, NEP and the
manual project bundle are disabled in this study to avoid overlapping additions.
`new_link_capacity: keep` installs the selected DC ratings. `skip: []` permits
AC upgrades. Native geometry-based duplicate detection remains active.

`reference_year: 2035` removes not-yet-built and expired electricity branches
**after applying replacement lifetimes, before clustering**. In particular,
an upgraded line replaces its old line rather than adding both capacities.
This option defines one static grid year and is intended for overnight studies;
leave it unset for a multi-period transmission model.

Every native project import writes
`resources/<run>/transmission_projects/project_mapping.csv`. It records accepted
additions, existing equivalents, upgrades, outside-area exclusions and failed
endpoint/upgrade mappings. Inspect this file alongside the converter audit.
An input row count alone is not evidence that all intended capacity survived
mapping. Two distant offshore endpoints receive distinct bus IDs; a branch that
maps to a single bus is dropped and logged. CSV content changes are explicit
Snakemake dependencies.

The native OSM 0.7 import on 25 October 2026 produced this audit:

| Input table | Added/replaced | Existing equivalent | Same bus | Unmatched upgrade |
|---|---:|---:|---:|---:|
| New AC lines | 51 | 5 | 4 | 0 |
| New DC links | 17 | 0 | 0 | 0 |
| AC upgrades | 21 | 0 | 0 | 1 |

The unmatched upgrade is **47/689, Vöhringen–Westtirol**; its original native grid
representation remains. The four same-bus entries are 1088/1755, 1088/1978,
1182/1997 and 1250/2056. These omissions require review if more detailed internal
grid representation is introduced. The mapping result is 68 additions and 21
replacements, not 99 additions. Harmony maps from the Elk Bis bus in Poland to
`way/258321412-330` in Lithuania; this distant native receiving-bus proxy is another
reason not to interpret the partial-route parameters as an engineering model.

The whole Bornholm project (ID 1106) is removed from this background input. Its
wind, hub and candidate cables are added only by the energy-island module.
The ordinary existing Bornholm–Sweden supply cable is not deliberately removed.

## Fixed and expandable backgrounds

`electricity.transmission_limit: vfixed` fixes prepared background AC and DC
capacities, without adding a global expansion budget. Bornholm links are added
after this step and retain their own fixed or optimisable bounds. Hydrogen and
gas infrastructure keep the same native investment options in all cases.

The two expandable benchmarks use `vopt`, the same installed planned assets,
and the native extension limits (20 GW per AC line, 30 GW per DC link). They
allow reinforcement of the represented network, not arbitrary new onshore
corridors. Background capacity already installed is treated as sunk in the
fixed cases. For comparisons across grid regimes, account consistently for any
fixed investment cost offsets included by the native expandable formulation;
prefer the difference between project and baseline within each regime.

See [the scenario and run guide](bornholm_zonal_study.md) for the eight runs.
