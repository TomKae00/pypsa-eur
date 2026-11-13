# Bornholm 2035 campaign for the final investment comparisons

This opt-in campaign retains the reviewed 2035 grid, 2013 weather, three-hour
chronology, EUR2025 biomass price of 30/MWh_fuel, native carbon/resource budgets
and the existing 16-thread Gurobi barrier preset. Earlier results and scenario
definitions are preserved. The new prefix is `bornholm_2035_final_v1`.

## Submit on Sophia

From the existing project checkout, after updating the branch:

```bash
cd ~/projects/pypsa-eur-bornholm
git checkout feature/bornholm-energy-island
git pull --ff-only
sbatch scripts/bornholm_2035_workq.sbatch all
```

The controller rebuilds the affected native inputs using the existing distributed
workq build profile. It validates all fifteen composed networks before submitting
any annual solve. The solve array uses at most six nodes simultaneously by default
(`BORNHOLM_2035_CONCURRENCY=1..15` can override this). A collector runs after all
tasks finish, including failed or numerically flagged tasks. Submission/job IDs are
recorded under `results/bornholm_2035_final_v1/submissions/<controller job>/`.

Other stages: `plan`, `dry-run`, `build` (build and check without solving), and
`check` (check existing composed inputs without rebuilding). Run these with the
same controller script. The walltime is a limit, not a completion-time prediction.

## Run matrix

| Tasks | Demand/transition world | Designs |
|---|---|---|
| 0–2 | Reference, EV45 | Baseline, hybrid, split radial |
| 3–5 | EV55 and faster industrial transition | Same three designs |
| 6–8 | EV65 and faster industrial transition | Same three designs |
| 9–10 | Reference only | 3 GW radial DK; 3 GW radial DE |
| 11–12 | Reference, V2G0 | Baseline, hybrid |
| 13–14 | Reference, V2G50 | Baseline, hybrid |

The first eleven cases use V2G25. Managed-charging participation remains 50%
in every case. The two single-country radial cases have no extra sensitivities.
Split radial retains exactly 1 GW of wind and cable to DK2 and 2 GW to Germany.
It uses two separate offshore buses, with the same Bornholm wind profile and
collection-cost basis. No offshore interconnector joins those buses. The hybrid
retains 1.2 GW to DK2 and 2 GW to Germany, so topology and total export capacity
both differ between these complete designs.

## Offshore availability

Background `offwind-ac`, `offwind-dc` and `offwind-float` correction factors fall
from 0.8855 to 0.85, a roughly 4.0% reduction in available output before clipping.
This uniform study derating affects existing and new generic offshore resources;
it is not a fleet-vintage calibration. Bornholm keeps its original 0.8855 factor
through a dedicated project-profile override. Its realised CF remains endogenous
because curtailment can change. Native existing wake-loss correction is not
applied again elsewhere.

## Native demand assumptions and additional fixed electricity

| Driver | Reference | EV55 world | EV65 world |
|---|---:|---:|---:|
| `sector.land_transport_electric_share[2035]` | 0.45 | 0.55 | 0.65 |
| Combustion share (fuel-cell share stays zero) | 0.55 | 0.45 | 0.35 |
| `industry.sector_ratios_fraction_future[2035]` | 0.45 | 0.55 | 0.65 |
| Additional DK data-centre load, TWh/y | 0 | 3 | 6 |
| Additional German data-centre load, TWh/y | 0 | 40 | 80 |

EV55/EV65 are explicit faster fleet-transition assumptions, not new-sales shares
and not reproduced official national scenarios. The native EV and industrial
settings apply across all modelled countries. Increased EV electricity use comes
with the associated fleet batteries and charging capacity. Industry progress
interpolates between current and future process/fuel intensities; it is not the
share of industrial energy supplied by electricity, nor general industrial growth.
Other route/recycling assumptions and material-production volumes are unchanged.

These assumptions represent faster electrification and industrial process
conversion than the reference, together with additional large data-centre demand.
The reference has no separately imposed future data-centre growth beyond what is
already embedded in the underlying electricity profiles. It is therefore the
conservative reference for these drivers, not a claim to reproduce a complete
official central 2035 forecast.

Data-centre additions are incremental study inputs. Danish additions are allocated
two-thirds to DK1 and one-third to DK2, close to the AF26 regional allocation.
Large centres connect to the regional AC bus and use a flat annual profile. No
flexibility, behind-the-meter generation or heat recovery is added for them.
AF26 reports 10.55086 TWh of total Danish data-centre demand in 2035
(`Elforbrug!L135`; DK1 7.13422 and DK2 3.41664 TWh, cells L133/L134).
This is a reference for scale, not the increment applied here. The German additions
are conditional low/high growth stresses, compared for scale with BNetzA's
2037 data-centre totals of 78.4/97.3/116.2 TWh. The existing profile already
contains some data-centre demand; none of those totals is simply added in full.

Sources: [AF26 dataset](https://ens.dk/media/8833/download),
[AF26 data-centre methodology](https://ens.dk/media/8838/download),
[BNetzA 2025 approval, Table 14, p.47](https://www.netzentwicklungsplan.de/sites/default/files/2025-05/250430_Genehmigung_Szenariorahmen_2025.pdf),
[native transport](https://pypsa-eur.readthedocs.io/en/latest/design/transport/),
[native industry](https://pypsa-eur.readthedocs.io/en/latest/design/industry/).

The old +10%/+20% fixed-load multiplication is not applied. The new worlds aim
at materially higher demand, approximately those orders of magnitude, but do not
promise exact total-consumption percentages. Native conversion choices, efficiency
and dispatch affect the outcome. Sector inputs are identical across designs
within each world and are never retuned separately to hit a demand percentage.

Local hydrogen production and industrial use remain available; pipelines and
external imports remain disabled. Electrolyser electricity is endogenous and is
not added again to fixed demand. The 2045 hydrogen/network/import design is deferred.

## Independent V2G participation

The optional `sector.bev_v2g_availability` is a share of all BEVs, bounded by the
managed-charging share. It partitions the native fleet into V2G-capable and
remaining cars, each with its own mobility demand, chargers and accessible battery
energy. Thus V2G25 cannot discharge the batteries of drivers permitting only
managed charging. Total mobility demand, total charger nameplate and total managed
battery energy stay unchanged across V2G0/V2G25/V2G50. Default null preserves the
previous pooled native model, so earlier configurations keep their behaviour.

## Checks and outputs

Existing onshore/PV ceilings from the central welfare campaign are reused,
including all resource classes and committed capacities. The existing numerical,
physical, biomass and reviewed-transmission checks are retained. An optimal result
flagged by the existing numerical screen is retained as `optimal_needs_review`;
it is not promoted to the accepted target and does not prevent other results being
collected. Code/config/source-network hashes protect prepared and queued inputs.

Outputs include `campaign_status.csv`, `input_demand_by_carrier.csv`,
`electricity_use_by_carrier.csv`, `electricity_use_totals.csv` and `transfer/`.
The electricity-use totals include fixed electricity, mobility including charging
losses, and actual electricity withdrawals for sector/fuel conversion. They exclude
electrical grid transfers and stationary battery cycling. EV driving is measured
from final battery demand, not gross V2G recharge. Electrolysis and other conversion
uses remain separate rows. These are explicitly defined model accounts, not a
replacement for official gross national electricity statistics.

Use `rsync -avhL` for the transfer folder so candidate-network symlinks are followed.
Do not call a completed collector evidence that every annual solve succeeded.

After this campaign, use two additional coherent weather years for fixed-capacity
dispatch of the central baseline/hybrid pair. Those four operational tests are
separate from this submission. Myopic investment and detailed country/CfD welfare
accounting are later stages.
