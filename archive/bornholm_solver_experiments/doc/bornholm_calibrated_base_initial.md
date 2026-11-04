# Bornholm: calibrated common base

This opt-in revision combines three onshore resource classes, regional onshore
profile adjustments, a fixed coal buffer and directional electricity-transfer
ceilings. It applies the same assumptions to all eight Bornholm cases. The
2013 weather, 2035 system horizon, demand assumptions, technology costs and
native land-eligibility calculation remain the basis of the study.

## Assumptions overview

These are pragmatic scenario assumptions supported by a plausibility review,
not fitted estimates of ERA5 error. A reference-fleet available CF around
20–30% is the working guide for onshore wind; stronger individual classes can
exceed it. No constraint forces the optimised country mix to attain a target CF
or a prescribed wind/solar/offshore share.

The following preview uses the uploaded `candidate_2035.nc`, SHA256
`212a8310cb9238d933b59a98e7773721112174437b7d3862636d702778cca68d`.
Its two existing classes are held fixed solely to illustrate the multipliers.
The reference weights are fixed `p_nom` for nonextendable generators and
`p_nom_min` for extendable generators, never `p_nom_opt`. They are the model's
required background fleet, **not a reconstructed 2013 fleet**. Snapshot generator
weights define annual CF. The rebuilt three-class results must be checked again.

The full country preview is in
[bornholm_calibration_preview.csv](../../../doc/data/bornholm_calibration_preview.csv).

| Country | Reference available CF before | Multiplier | Preview after clipping | Decision |
|---|---:|---:|---:|---|
| Denmark | 42.1% | 0.60 | 25.3% | Reduce clear high outlier |
| Netherlands | 35.0% | 0.85 | 29.8% | Moderate reduction |
| Great Britain | 38.7% | 0.65 | 25.2% | Reduce excessive fleet/best-class levels |
| Ireland | 33.9% | 0.80 | 27.2% | Moderate reduction; retain stronger coastal sites |
| Austria | 16.9% | 1.35 | 22.3% | Raise low availability |
| Czechia | 17.5% | 1.25 | 21.7% | Rounded Central European scenario assumption |
| Switzerland | 6.1% | 4.00 | 19.6% | Provisional proxy; large uplift needs particular scrutiny |
| Germany | 22.7% | 1.00 | 22.7% | Retain |
| Poland | 25.1% | 1.00 | 25.1% | Retain |
| Sweden | 20.5% | 1.00 | 20.5% | Retain |
| Norway | 23.2% | 1.00 | 23.2% | Retain |
| Finland | 21.6% | 1.00 | 21.6% | Retain |
| Belgium / France | 29.0% / 23.5% | 1.00 | Unchanged | Retain |

The [IEA Wind TCP 2015 report, table 13](https://www.seai.ie/sites/default/files/publications/IEA_Wind_Annual_Report_2015.pdf)
provides historical context: Austria 24%, Switzerland about 20%, Ireland 30.5%
and Dutch onshore 22.3% in 2013. It also shows country and weather-year variation.
Its Danish figure combines onshore and offshore and is not an onshore target.
[CBS's onshore series](https://dashboards.cbs.nl/longreads/2026/hernieuwbare-energie-in-nederland-2025/H4_windenergie.html)
and [UK government load factors](https://www.gov.uk/government/publications/renewables-obligation-level-calculations-2026-to-2027/calculating-the-level-of-the-renewables-obligation-for-2026-to-2027)
support using broad fleet ranges. These sources do not determine the exact
rounded multipliers above, nor validate a future fleet using a historical year.
Czechia's uplift is an explicit modelling approximation rather than an independently
estimated national bias. Switzerland is the weakest assumption: scaling cannot
repair missing local terrain detail or prove that the spatial distribution is right.

Each affected country's three onshore classes receive the same multiplier on
their **hourly availability**, clipped to [0, 1] before temporal averaging. This
preserves chronology, zeros and within-country class distinctions; upward
clipping changes peaks. The audit records its extent (the old-network preview
clips about 1.7% of Austrian, 1.4% of Czech and 3.9% of Swiss profile entries).
The native class-specific finite potentials and cost assumptions are retained.
These remain geographical eligibility ceilings, not commitments that the whole
potential can be deployed by 2035. No additional guessed national build-rate
limit is introduced in this revision.
Three classes partition the native eligible resource rather than tripling it.
The extra class therefore requires rebuilding profiles from the cutout; it cannot
be created faithfully by editing the old solved network.

Balkan countries AL, BA, BG, GR, HR, ME, MK, RO, RS, SI and XK receive no regional
multiplier. The global class-count change can still alter their class partition.
Other countries retain factor 1 unless listed in the committed calibration CSV.

## Offshore, solar and the investment interpretation

The audited fixed/rooftop solar and tracking-solar ranges have a coherent
geographical pattern; no extra adjustment is applied. Offshore retains its
existing 0.8855 correction factor, native turbine assumptions and three classes.
High offshore classes deserve inspection but do not, by themselves, justify a
second blanket loss factor. The fixed Bornholm project's profile is untouched.

The old solved file has `suboptimal` status, so its investments are provisional.
It also **does build solar in Germany**: about 104.6 GW additional rooftop PV,
alongside about 36.2 GW additional onshore wind. Generic utility PV alone misses
the rooftop investment on low-voltage buses. Rooftop and utility PV have different
costs and grid value. Changing renewable availability can change the mix; it does
not guarantee offshore or utility PV investment. Diagnose cost, connection cost,
potential and curtailment effects after the common corrected solve.

## Coal buffer

Only coal fuel stores are fixed at `e_nom = 1000 MWh`,
`e_nom_extendable = false`, cyclic, zero initial inventory, zero standing loss
and ordinary 0–100% energy bounds. **MWh belongs to `e_nom`, not `p_nom`.**
The independent fuel-supply generator still supplies annual coal demand, so
1,000 MWh is a buffer, not an annual fuel quota or 1,000 MW supply limit.

The uploaded network had about 1.165 billion MWh optimised coal capacity but only
an 8.35 MWh annual inventory swing. Removing this redundant investment dimension
addresses that numerical pathology. It does not prove that every other numerical
problem has disappeared. Post-solve checks verify the buffer bounds and the
existing numerical acceptance screens. Gas and other stores retain their native
formulation. Fixed buffer cost is a common constant across the revised cases.

## Directional transfer assumptions

The [SvK transfer map](https://www.svk.se/siteassets/press-och-nyheter/temasidor-bilder-och-block/bilder/svk_karta_elomraden_overforingskapacitet.pdf)
is a published **2021 nominal-maximum reference**, not hourly historical NTC or
today's flow-based capacity domain. It supports the directional approximations
below. A 2023 row set is supplied for a static hindcast approximation; the 2035
study selects the 2035 rows independently of its 2013 weather year.

| Direction (reverse is opposite) | 2023 forward / reverse MW | 2035 forward / reverse MW | Treatment |
|---|---:|---:|---|
| SE1 → SE2 | 3300 / 3300 | 3300 / 3300 | Retain published maximum |
| SE2 → SE3 | 7300 / 7300 | 10350 / 10350 | SvK LMA2026 bidirectional 2035 stage |
| SE3 → SE4 | 6200 / 2800 | 6200 / 2800 | Common cap across parallel AC and DC |
| DK2 → SE4 | 1700 / 1300 | 1700 / 1300 | Retain directional maximum |

[SvK LMA2026, table 10, printed page 54](https://www.svk.se/49da32/siteassets/om-oss/rapporter/langsiktig-marknadsanalys-2026/langsiktig-marknadsanalys-2026.pdf)
assumes SE2↔SE3 stages of 8,100 MW in 2028, 9,600 MW in 2034, 10,350 MW in 2035
and 10,700 MW in 2040. The table uses a bidirectional arrow. This model therefore
uses **10,350 MW in both directions for 2035**, adopting the completed-2035-stage
scenario rather than claiming guaranteed hourly market availability. The earlier
9,600/7,300 MW asymmetric entry was unsupported and has been replaced.
[LMA2024, table 9](https://www.svk.se/49006b/siteassets/om-oss/rapporter/2024/lma_2024.pdf)
also used a bidirectional arrow for 9,600 MW from 2034, before the extra 2035 stage
was added in LMA2026. The generic NordSyd webpage is less specific than these
scenario tables. The other unchanged 2035 maxima are static assumptions, not
claims of verified 2035 adequacy.

The remaining directional pairs are also reproduced in the
[Nordic TSOs' maximum-NTC map, valid from January 2024](https://www.nordpoolgroup.com/globalassets/download-center/tso/max-ntc.pdf).
In particular, DK2→SE4 is 1,700 MW and SE4→DK2 is 1,300 MW; SE3→SE4 is 6,200 MW
and SE4→SE3 is 2,800 MW. This confirms the orientation, not invariant availability
throughout the year. [The wider transmission review](../../../doc/bornholm_transmission_review.md)
records additional Norway–Sweden, internal Norwegian, Danish and Finnish issues
from the earlier audit. Those are not all corrected by the four active rows.

These limits supplement native AC physics and installed TYNDP projects. They
are not substituted for line MVA ratings, do not receive another 0.7 multiplier,
and are not added to project capacities. Each constraint sums all directly
crossing AC lines and electricity DC links once, using the exporting terminal's
MW and native link efficiency. Reverse-link records are handled by their actual
terminals. The limits are a zonal transport approximation, not an N-1 or full
flow-based assessment. Nonreviewed borders retain their native assumptions.

For the two `vopt` cases, additional investment raises each ceiling proportionally
to its directional physical envelope above the composed minimum capacities:

`limit = base_limit + (base_limit / base_envelope) * added_envelope`.

The baseline envelope already includes the planned grid, so those projects are
not added twice. This is an explicit expansion proxy that preserves reinforcement
options; it is not a security assessment of expanded infrastructure. A `fixed`
policy is also supported. A nominal ceiling can remain above a corridor's smaller
physical envelope; native physical constraints then bind first.

## Build and verification

Use the overlay last, in a new namespace:

```bash
pixi run python -m scripts.run_bornholm_calibrated --list-build-targets
# Optional read-only check against the uploaded two-class candidate:
pixi run python -m scripts.run_bornholm_calibrated --preview-network /path/to/candidate_2035.nc --preview-output results/calibration_preview.json
sbatch scripts/bornholm_calibrated_workq.sbatch build
pixi run python -m scripts.run_bornholm_calibrated --preflight-all
# Build and then submit the eight annual solves:
sbatch scripts/bornholm_calibrated_workq.sbatch all
```

The controller uses the existing Sophia workq bridge; each solve reserves 16
CPUs. Builds and solves are separate jobs. The default remains chronological
three-hour averaging (2,920 snapshots, 8,760 annual weight hours). This revision
does not overwrite or reuse the previous `w2` composed/solved networks. Native
wind/PV profiles and class-region files are rebuilt once under the new
`resources/bornholm_3h_s1_w3_o3_calibrated/shared/` namespace, shared by the eight
cases. Existing weather cutouts and eligibility matrices remain reusable.

Composition records input-table hashes, per-asset before/after CF and clipping,
country reference-fleet CF and the exact corridor membership in network metadata.
The solver rejects missing/stale assumption metadata, and preflight reports the
rebuilt class results. Neither preflight nor a multiplier preview is a 2035 solve.
Inspect solved fleet-weighted available CF, realised CF, curtailment, binding
transfer limits and coal buffer utilisation after the full cluster runs.

Local tests exercise small optimisations, including asymmetric and parallel
corridors, lossy reverse links and reinforcement, plus profile calibration,
fixed coal supply and metadata round-trips. A full rebuild requires the weather
cutouts and preprocessing inputs on Sophia; a local test or old-network preview
does not replace that rebuild and annual solve.

The new assumptions, study-controller and accuracy suites passed 359 checks:

```bash
pixi run python -m pytest -q test/test_bornholm_calibrated.py test/test_bornholm_study.py test/test_bornholm_accuracy.py
```

The full uploaded network also passed the read-only multiplier/corridor preview.
The zonal, island-project, planned-grid and hourly suites passed another 64
checks, giving 423 passing tests across the relevant suites.
Ruff checks, formatting, Bash syntax checks and Git whitespace checks passed.
Local native preprocessing tests use an environment-only fallback from psutil
to the OS CPU-affinity reading, because sandboxed child PIDs cannot always be
resolved by psutil. That accommodation is not part of the model or this branch.
