# Nordic physical nominal-rating investigation

The Nordic nominal ratings have **not** been corrected by this investigation.
It establishes the original circuit inventory, reviews public physical references,
and identifies the specific missing data needed for defensible replacements.
Neither a reduced border transfer limit nor a conductor bundle count alone is
used as an ampacity.

## Findings on the supplied original network

The supplied `base_extended.nc` and its two busmaps contain 877 original AC line
records incident to DK, FI, NO or SE, including 56 that connect different model
zones. There are also 19 DC branch records and 24 other/converter link records.
These counts describe this input topology, not a historical asset census.

| AC corridor | Original records | Circuits | Sum of nominal ratings, MVA |
|---|---:|---:|---:|
| SE1–SE2 | 4 | 4 | 7,149.91 |
| SE2–SE3 | 11 | 11 | 15,774.48 |
| SE3–SE4 | 5 | 5 | 8,937.38 |
| DK2–SE4 | 2 | 2 | 3,574.95 |
| NO1–SE3 | 2 | 2 | 3,753.70 |
| NO3–SE2 | 1 | 1 | 1,876.85 |
| NO4–SE2 | 1 | 1 | 491.56 |
| DK1–DE | 3 | 4 | 7,149.91 |

These are sums of `s_nom`, before any security factor. They are not simultaneously
deliverable active power, market capacities, or a proposal to set corridor limits.
Norwegian model identifiers are translated here to their usual bidding-zone names.

### SE2–SE3: uniform electrical assumptions hide different construction

The eight 400-kV records each have the native `Al/St 240/40 4-bundle 380.0`
type, corresponding to 2,580 A and 1,787.48 MVA per circuit at model voltage.
The three 220-kV records each have 1,290 A and 491.56 MVA.

[Svenska kraftnät's 2023 Munga–Hamra application](https://www.svk.se/4a7f0c/siteassets/2.utveckling-av-kraftsystemet/transmissionsnatet/utbyggnadsprojekt/munga-hamra/koncessionsansokan/koncessionsansokan.pdf)
independently identifies eight 400-kV and three 220-kV lines across Snitt 2.
Agreement in count does not validate the assumed conductor types or currents.

Current OSM identity observations on model route members include:

| Route identity | Observed bundle tags on members | Model type |
|---|---|---|
| CL1 S2, Midskog–Kättbo | double / 2 | four-bundle |
| CL22 S1–2 | double | four-bundle |
| CL3 S3 | double | four-bundle |
| CL4, CL5, CL7, CL26 and UL17 members | triple on inspected members | four-bundle |

The observations expose a construction-assumption mismatch worth investigating.
They do not establish conductor cross-section, continuous current, terminal
limits, or the historical construction of every member. Multi-voltage routes
also have mixed tags. **Do not scale the rating by 2/4 or 3/4.**
The review manifest records individual OSM URLs and original model signatures.

### Øresund: circuit-specific cable power is published

The joint TSO [SOA appendix 6, section 4.1.1](https://www.fingrid.fi/globalassets/dokumentit/fi/kantaverkko/kantaverkon-kehittaminen/appendix-1-9-to-soa-annex-os-v5.pdf)
lists 830 **MW per cable circuit** for FL23 and FL25 in each temperature column.
That is stronger circuit evidence than the overall border TTC. The model instead
uses 1,787.48 **MVA per circuit** from its generic overhead type.

This is an important unresolved discrepancy. The published active-power value
does not establish the nominal MVA or continuous current. No power-factor or
reactive-current allowance is silently inferred. The two model IDs are
`relation/5528438:b-400` and `relation/5528440:b-400`; their individual allocation
to FL23/FL25 is unresolved, though the two published powers are equal.

The source's 130-kV links and the existing Bornholm cable also lie below the
model's voltage cutoff. Their representation is a separate topology issue;
it is not repaired by changing these two 400-kV nominal ratings.

### Hasle–Borgvik: a specific conductor reference exists

[NVE's asset description](https://www.nve.no/om-nve/nves-listefoerte-kulturminner/kraftledninger/hasle-borgvik/)
identifies FeAl 2×326 duplex conductors. The original Norwegian model record
`way/585198893-420` uses the generic four-bundle type. This corroborates that
the blanket type is not a verified physical description. A continuous circuit
ampacity, its reference conditions and the limiting terminal equipment remain
necessary before replacing `s_nom`. The Norwegian record is only part of the
international circuit; a national segment rating cannot certify its other end.

### DC: keep installed power and reference points explicit

Sixteen of the 19 DC branch records match an explicitly recorded project or
installed reference value. This is a numerical comparison, not validation of
historical availability. Sources and their rating bases remain in the manifest.

Three records remain on hold: Konti-Skan and the two model Skagerrak aggregates.
The [2024 ENTSO-E HVDC report](https://eepublicdownloads.entsoe.eu/clean-documents/SOC%20documents/Nordic/2024/HVDC_Utilisation_and_Availability_Statistics_2024.pdf)
distinguishes sending/receiving conventions and reports the Konti-Skan convention
change in 2020. Its receiving-end Skagerrak values must not be substituted for
the manufacturer's historic pole ratings or later joint TSO nominal figures.

Fenno-Skan 1 is a useful counterexample: SOA appendix 2 explicitly identifies
500 MW installed power and separately explains operation at 400 MW due to
reduced voltage. The model's 500 MW matches the installed value. The restriction
is recorded without changing the physical nominal input in this work package.

## Reproduce the audit for any model year

```bash
python scripts/audit_nordic_nominal_ratings.py \
  resources/hindcast-build-2023/networks/base_extended.nc \
  --year 2023 \
  --simplify-busmap resources/hindcast-build-2023/busmap_simplify_network.csv \
  --cluster-busmap resources/hindcast-build-2023/busmap_cluster_network.csv \
  --output data/hindcast/line-ratings/nordic-review-2023
```

Change the year and corresponding input paths for another build. Both busmaps
are optional together; they are needed for corridor reports and data-request
tables. The original-voltage network is essential. A simplified or composed
network has already lost the individual voltage/type assumptions being examined.

The script reuses the original inventory functions and produces AC and DC
inventories, a separate converter inventory, AC corridor sums, a review JSON,
and `nordic_ac_rating_data_request.csv` containing all interzone AC records.
It performs no downloads, optimization or network export. Changed model
signatures and missing records are reported for re-review; matches are not
silently carried to another topology.

`config/hindcast/nordic-rating-review.yaml` records 31 AC identity/construction
reviews and 19 DC comparisons. It is **not** a correction input and is not added
to the automatic nominal-rating study. AC records without evidence remain
unreviewed. Existing reviewed corrections elsewhere in Europe are unaffected.

Model year and source edition remain independent. The public SOA bundle contains
appendices of different dates, including later revisions. The ENTSO-E report
describes 2024 and was published in 2025. Neither is silently relabelled as a
complete 2022/2023/2025 historical dataset. Source fingerprints identify the
reviewed bytes; they are not assertions that mutable web pages remain unchanged.

## Data needed to complete the corrections

The generated request table identifies the exact original circuits and their
assumptions. For each circuit, the remaining asset-owner information is:

* TSO circuit code and endpoint correspondence, including every parallel circuit.
* Continuous phase current in A, or continuous apparent power in MVA, with
  reference voltage and the applicable temperature/cooling assumptions.
* Conductor or cable construction and limiting terminal equipment.
* Source edition and validity dates, including rating-changing refurbishments.

Prioritize the eleven Snitt 2 circuits, FL23/FL25, the SE3–SE4 circuits, and the
Norway–Sweden circuits. The same table also includes other Nordic interzone lines.
This is a prepared data specification; no request has been sent to a TSO.

Once those fields are supported, add normalized current references and reviewed
bindings to the existing nominal-rating compiler. No new optimization formulation
or independent network aggregation is needed. Until then, the Nordic physical
rating correction remains incomplete and should be reported that way.
