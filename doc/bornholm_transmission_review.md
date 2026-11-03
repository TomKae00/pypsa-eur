# Bornholm: wider Nordic transmission review

Verification on 3 November 2026 of the previous Nordic nominal-rating audit and
the four active zonal transfer overrides. This review corrects the SE2–SE3
planning assumption and records the wider outstanding scope. It does not add
unverified physical ratings or automatically turn every reference below into a
2035 operational ceiling.

## SE2–SE3: verified planning stages and both directions

[SvK LMA2024, table 9, printed page 40](https://www.svk.se/49006b/siteassets/om-oss/rapporter/2024/lma_2024.pdf)
uses SE2↔SE3 for 8,100 MW in 2028, 9,600 MW in 2034 and 10,500 MW in 2040.
Thus 9,600 MW has a published scenario basis, but not the asymmetric interpretation
previously implemented here.

[SvK LMA2026, table 10, printed page 54](https://www.svk.se/49da32/siteassets/om-oss/rapporter/langsiktig-marknadsanalys-2026/langsiktig-marknadsanalys-2026.pdf)
updates the stages:

| Stage year | Bidirectional capacity assumption, MW |
|---|---:|
| 2028 | 8,100 |
| 2034 | 9,600 |
| 2035 | 10,350 |
| 2040 | 10,700 |

The revised base uses **10,350 MW in both directions in 2035**. This represents
the completed-2035-stage scenario. The original 9,600 MW southbound / 7,300 MW
northbound pair is withdrawn. A delay assumption would need to be stated as such;
the 9,600 MW stage itself is also bidirectional in the source. These scenario
values remain distinct from a guaranteed hourly market-capacity domain.

## Additional connections from the earlier audit

The earlier investigation is `nordic-nominal-ratings.md` and the accompanying
`nordic_ac_interzone_review.csv`, `nordic_ac_corridor_nominal_sums.csv`,
`nordic_ac_rating_data_request.csv` and `nordic_dc_nominal_review.csv`. It applied
zero Nordic nominal-rating corrections: missing continuous-current references
and differing rating conventions were explicitly left unresolved.

The table below compares those findings with the uploaded hybrid network,
SHA256 `212a8310cb9238d933b59a98e7773721112174437b7d3862636d702778cca68d`.
The model figures are sums of individual branch operating envelopes in the
linear model, after the existing loading factors. They are not independent
estimates of simultaneously deliverable transfer capacity. They come from the
input representation and do not depend on accepting the suboptimal dispatch.

| Additional corridor | Uploaded model envelope, GW | Comparison and outstanding treatment |
|---|---:|---|
| NO1–SE3 | 2.628 each way | Historical maximum NTC: NO1→SE3 2.145 and reverse 2.095 GW. Physical Hasle–Borgvik/Halden ratings also need review. |
| NO3–SE2 | 1.314 each way | Historical maximum NTC: NO3→SE2 0.600 and reverse 1.000 GW. Check 2035 reinforcement before setting a final override. |
| NO4–SE1 | 1.251 each way | Historical maximum NTC: NO4→SE1 0.700 and reverse 0.600 GW. Physical circuit rating remains unverified. |
| NO4–SE2 | 0.344 each way | Historical maximum NTC: NO4→SE2 0.250 and reverse 0.300 GW. |
| DK1–DE | 5.005 each way | LMA2026 assumes 3.500 GW bidirectionally from 2027. Existing AC nominal ratings and planned reinforcement require one consistent treatment. |
| SE1–FI | 5.005 each way | Strong review priority. Account for Aurora Line once and distinguish original circuit ratings, commissioning year, TTC and market margins. The joint SOA appendix 2 gives normal AC TTC of 2.100 GW each way, with further operating restrictions. |
| FI–SE3, Fenno-Skan | 1.300 each way | Model keeps 500 + 800 MW installed. SOA appendix 2 distinguishes 400 MW normal operation on Fenno-Skan 1 from its 500 MW installed rating, giving 1,200 MW normal aggregate operation. |
| NO1–NO2 | 7.367 each way | Historical maximum NTC is 2.200 GW NO1→NO2 and 3.500 GW reverse. Norwegian internal constraints were inventoried but not physically calibrated. |
| NO1–NO5 | 6.757 each way | Historical maximum NTC is 0.300 GW NO1→NO5 and 3.900 GW reverse. Material directional issue. |
| Other internal Norwegian borders | NO1–NO3 0.704; NO2–NO5 3.519; NO3–NO4 2.018; NO3–NO5 1.314 | Included in the earlier inventory but recorded as unreviewed native assumptions. Reconcile operational constraints and 2035 reinforcements before replacing them. |
| NO4–FI | 0.344 each way | Original AC connection has conditional operation and trading through neighbouring borders. A planned 150 MW back-to-back connection in 2032 is a separate LMA2026 assumption, not evidence for 344 MW unrestricted trade. |

Historical maximum NTC comparisons are from the
[Nordic TSO map valid from 1 January 2024](https://www.nordpoolgroup.com/globalassets/download-center/tso/max-ntc.pdf).
Its revision log also explains internal Norwegian directional changes. These
values are comparison references, not automatically adopted 2035 limits.
[Statnett's transport-channel study](https://www.statnett.no/contentassets/fa81f0dc3b074d309af63f90b73e2a61/analyse-av-transportkanaler-2021-2040-revidert-6.10.pdf)
separately distinguishes NO3→SE2 600 MW from SE2→NO3 1,000 MW.

## DC links: retain the rating reference point

The earlier DC audit matched 16 of 19 branch records to recorded installed/project
references. The three outstanding records are Konti-Skan and the two Skagerrak
aggregates. It did not establish that all DC links were incorrectly rated.

- **Konti-Skan, DK1–SE3:** model 680 MW. Joint SOA appendix 5 gives 715 MW at the
  receiving end in bipolar operation; the ENTSO-E report distinguishes sending-end
  capability. Review nominal capacity, link losses and operational limits together.
- **Skagerrak, DK1–NO2:** model 940 + 700 = 1,640 MW. Joint SOA appendix 4 gives
  different nominal conventions and distinguishes 1,680 MW day-ahead sending-end
  from 1,632 MW intraday receiving-end capacity. These are not interchangeable
  replacement values for `p_nom`.
- **Great Belt, DK1–DK2:** model installed power is 600 MW, which matches its
  nominal project rating. The same SOA gives 590 MW DK1→DK2 and 600 MW reverse
  trading capacity. Confirm the loss/reference-point convention before calling
  this a nominal-rating error.

References: [joint Nordic SOA appendices](https://www.fingrid.fi/globalassets/dokumentit/fi/kantaverkko/kantaverkon-kehittaminen/appendix-1-9-to-soa-annex-os-v5.pdf)
and [ENTSO-E's 2024 HVDC report](https://eepublicdownloads.entsoe.eu/clean-documents/SOC%20documents/Nordic/2024/HVDC_Utilisation_and_Availability_Statistics_2024.pdf).

## Future projects and scope

LMA2026 is useful for planning assumptions, but its complete project list must not
be copied blindly. The report explicitly says that the Konti-Skan Connect upgrade
remained in its scenarios despite being paused in May 2026, after the scenario
cutoff. Its 1,000 MW 2035 entry therefore needs a separate project decision.
It also assumes Bornholm's 1,200 MW DK2 and 2,000 MW German connections in 2037;
this study's 2035 Bornholm cases are an explicit earlier-commissioning scenario.

The earlier TYNDP input review separately corrected Harmony Link's obsolete DC
representation and excluded SACOI3 and GRITA2 pending compatible endpoint/
replacement modelling; see [the planned-grid input review](bornholm_2035_grid.md).
Those interventions are separate from zonal transfer ceilings.

Only the existing four active corridors are constrained by
`data/model_assumptions/zonal_transfer_limits.csv`. This verification changes the
SE2–SE3 2035 row and its source. The wider pending list above is not represented
as a completed Nordic calibration. Physical assets, TYNDP capacities and operating
ceilings must be reconciled corridor by corridor so projects are counted once.
