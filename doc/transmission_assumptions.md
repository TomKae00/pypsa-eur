# Reviewed operating transmission assumptions

The reviewed cases add **operating MW constraints** to the physical network.
They retain the native AC impedances, MVA ratings, KVL and converter losses.
An AC rating inferred from voltage and line type is not an observed trading
capacity. The two quantities are not substituted for each other.

## Dated inputs

| Case | Operating reference | Interpretation |
|---|---|---|
| 2023 hindcast | 20 Nordic corridors from the NVE/RME 2023 report and cited TSO evidence | Static annual upper envelopes; hourly outages and joint Nordic restrictions remain incomplete |
| 2035 Bornholm | ENTSO-E ERAA 2025 TY2035 public NTC tables, with reviewed TSO/project overrides | Annual maximum NTC scenario approximation, including Core fallback NTCs; not a validated flow-based market domain |
| Poland 2035 | DE/CZ/SK legs and a shared 3,000 MW import/export profile | Conservative sum of positive directional border use; imports cannot cancel exports to obtain extra capacity |

`data/model_assumptions/reviewed_transfer_limits.csv` contains 20 historical
and 100 future corridor pairs, each with an explicit source and two directions.
References without a direct model crossing are recorded as skipped; no
physical asset is invented. Missing countries retain their native topology.
The CSV does not create a geographically reduced model.

The future reference comes from [ERAA 2025 modelling data](https://www.entsoe.eu/eraa/2025/modelling-data/),
`NTCs Consolidated TY2035.xlsx`. The historical reference is
[NVE/RME report 4/2024](https://publikasjoner.nve.no/rme_rapport/2024/rme_rapport2024_04.pdf),
including figure 41 and the 2023 maximum-trading-capacity tables. Annual maxima
are deliberately labelled as upper-envelope assumptions, not hourly available
capacity or evidence that every corridor can simultaneously reach its maximum.

## Specific corrections and timing assumptions

| Connection | Reviewed treatment |
|---|---|
| SE2–SE3 in 2035 | 10,350 MW in both directions, SvK LMA 2026 table 10. The 9,600 MW step belongs to 2034; the former 2035 value is superseded. |
| KontiSkan in 2035 | 1,000 MW replacement assumed delivered; later project delivery remains a sensitivity. |
| Kontek | DE→DK2 600 MW, DK2→DE 585 MW, using 50Hertz's explicit direction. No unconditional additional 400 MW is added for the wind-dependent KF CGS path. |
| EstLink in 2035 | Existing 1,000 MW retained. EstLink 3's later delivery is excluded from the firm base. |
| Fenno-Skan 1 | Engineering nameplate retained; 400 MW operating limit. Fenno-Skan 2 remains a separate 800 MW asset. |
| Italy–Montenegro | Existing first pole corrected to 600 MW; an added second pole must be counted separately. |
| Viking in 2023 | Zero before 30 December 00:00 UTC, then an initial 800 MW ceiling. This is a conservative first-full-day approximation after the 29 December opening, not a verified first commercial hour. |
| Shetland in 2023 | Zero availability; August 2024 energisation cannot be backdated. Island assets remain in the network and require their own adequacy review. |
| FAB, GridLink, Aquind | Identified proposals excluded from the firm 2035 background following the cited February 2026 CRE/Ofgem decision. |
| Bornholm project | The 2035 study remains an explicit early-delivery counterfactual relative to the later published project schedule. |

The asset identifiers, dates and sources are in
`data/model_assumptions/transmission_asset_rules.csv`. Rules are applied before
sector coupling creates reversed lossy converters and before temporal averaging.
Consequently, an unavailable hour cannot return as an unrestricted converter.

## Implementation and verification

`scripts/transmission_assumptions.py` installs corridor membership after final
composition, then adds Linopy constraints through both solve paths. Flows are
measured at the exporting zone's terminal. Parallel AC lines and both records
of a reversed lossy DC pair contribute to one net corridor expression.

An optional `hourly_limits` CSV has columns `snapshot,zone0,zone1,forward_MW,reverse_MW`.
For every overridden corridor it must cover every final snapshot exactly once;
timestamps are UTC. No forward-fill is performed. For a temporally averaged
case the supplied limits must already use that case's snapshot grid; the
implementation does not silently choose an hourly-to-block aggregation rule.

The historical policy is `fixed`. In Bornholm's `vopt` cases the individual
corridor ceilings scale with investment above the existing physical directional
envelope. This is a proportional planning assumption, not a TSO-certified
transfer benefit of reinforcement. Poland's joint 3,000 MW profiles remain fixed
in both grid regimes; relaxing those requires a separately justified scenario.
Within each grid regime, all eight cases use the same background assumptions.

Input hashes, the planning year, final snapshots, asset changes and the installed
constraints are retained in network metadata. Stale source tables, changed
snapshots or changed corridor topology/loss coefficients require recomposition.
The checks use the system horizon, never the 2013 Bornholm weather year.
Tests solve small systems for both directions, losses, reinforcement, hourly
availability and the shared Polish profile, and verify a NetCDF round trip.

## Remaining validation work

The historical Core flow-based feasible domain has **not** been reconstructed.
JAO MaxBEX values are conditional diagnostics and are rejected as input NTC
tables. The previous Core 2035 independent KPI comparisons did not reconcile;
they are not a basis for claiming validated Core flow-based restrictions.
Joint Nordic hourly restrictions, OL3-dependent FI–SE1 availability and complete
historical outages still require dated data. Physical cable ratings alone do
not supply those missing market constraints.

The 2022, 2024 and 2025 configurations remain separate cases awaiting their own
dated review. These changes do not silently apply 2023 limits to other years.
Hydro foresight/cyclicity, nuclear outages, reserves and CHP/heat assumptions
remain modelling limitations to examine against the first updated hindcast.
Passing numerical checks does not establish that every real-world assumption
is correct or provide an independent optimality certificate.
