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

## Central matched-cable experiment

Run from the campaign checkout on Sophia with the two original unsolved networks
and their prepared.json files still under the final-v1 resources and results
prefixes. Submit from the repository root.

```bash
sbatch scripts/bornholm_matched.sbatch
```

The array solves split radial (1/2 GW), hybrid (1/2 GW), and hybrid (1.2/2 GW),
in that order. It verifies the original input hashes before use. Both directional
models of the Danish cable are resized in the middle case. The original per-MW
line and converter costs and losses apply. Dedicated and shared fixed-station
costs are not independently itemised by this model. All three use the same new
numerical settings (NumericFocus 2, BarConvTol 1e-9) and the original acceptance
thresholds. Exit code 3 retains a flagged candidate for review. A new submission
has its own results and resources under bornholm_2035_matched/ARRAY_JOB_ID.

After all three tasks have finished, copy the transfer folder **from your laptop**,
substituting the Slurm array ID and the destination path to your local model checkout.

```bash
rsync -avL s253303@sophia.dtu.dk:projects/pypsa-eur-bornholm/results/bornholm_2035_matched/ARRAY_JOB_ID/transfer/ ~/bornholm-matched-results/transfer/
```

Open `doc/bornholm_matched_central.ipynb` in the local checkout of this branch and
set DATA_DIR to the copied transfer directory (absolute path). It reports both
increments and their total, solver flags and settlement residuals. Repeatability
means comparing the sizes and signs of the *economic* differences as well as the
solver diagnostics. All figures represent one operating year. The pooled cable-cost
and landing-country allocations bracket two accounting assumptions; replace the
shared equipment assumption with itemised costs when available. The notebook
uses a two-sided actual-output proxy with a 70% German charge or credit and equal
rent shares.

From the original executed central notebooks, the treaty hybrid saved €18.7651m/year
against split radial, while its annual cable charge was about €24.7m higher
(€365.1m versus €340.4m). Other system cost therefore fell by about €43.5m/year.
Those two terms do not separately identify the operational and capacity effects.
Under pooled export costs at €75/MWh, the Danish, German and developer gains were
+€4.78m, −€9.23m and +€9.99m. With both designs charging their bundled cables to
their landing countries and keeping congestion rent pooled, the incremental German
cable cost is unchanged (same German route and rating). The €24.7m incremental
cable charge moves entirely to Denmark. This gives approximately −€7.6m for Denmark
and +€3.1m for Germany. This is an ownership sensitivity from rounded previous
results, before allowing for jointly owned equipment. All original cases triggered
numerical quality flags.

The [Bornholm assessment README](bornholm_assessment/README.md) sets out the
decision method and shows how to run the strike and funding-share screen.

## Five-layout screen across the three demand pathways

Run this screen to compare connection choices before adding separate demand
and offshore stresses. The wind build remains 3 GW. Each pathway uses its
archived final-v1 split-radial and treaty-hybrid inputs.

| Layout | Denmark / Germany rating | Comparison |
|---|---|---|
| Split radial | 1.0 / 2.0 GW | Modelled fallback |
| Treaty hybrid | 1.2 / 2.0 GW | Negotiated configuration |
| Matched hybrid | 1.0 / 2.0 GW | Compare hybrid routing at equal ratings |
| Hybrid, more DK capacity | 1.4 / 1.8 GW | Shift 200 MW toward Denmark at the same total rating |
| Hybrid, more DE capacity | 1.0 / 2.2 GW | Shift 200 MW toward Germany at the same total rating |

Tasks 0–4 run Central, 5–9 EV55, and 10–14 EV65, in table order within each
group. This re-solves the six original split-radial and treaty-hybrid cases
alongside the nine additional combinations under the same solver settings.
The shifted ratings are modelling sensitivities. The model scales existing
per-MW cable costs and keeps the original losses. Both flow directions use
the new rating. Any change in mandatory transmission allowance is applied
to its matching global constraint to preserve the background grid budget.

From the Sophia checkout, run the input check and queue the dependent array.
The check reads all six prepared input networks and their provenance records
and verifies matched demand, snapshots, wind availability and background
resource budgets after resizing.

```bash
check_id=$(sbatch --parsable --array=0 --time=00:30:00 \
  scripts/bornholm_matched.sbatch --layout-screen --check) &&
sbatch --array=0-14%6 --dependency=afterok:"${check_id%%;*}" \
  --kill-on-invalid-dep=yes scripts/bornholm_matched.sbatch --layout-screen
```

The input check does not start a solver. Outputs go to
`results/bornholm_2035_layout_screen/<job>/transfer/`. Each case includes the
candidate network, prepared record and summary. An optimal candidate with
solver screening flags is labelled `optimal_needs_review` and the task exits
with code 3. Use comparable, screened runs for the final welfare comparison.

## Commitment stress runs

This optional nine-task array retains the 3 GW Bornholm build. It crosses
split radial (1/2 GW), treaty hybrid (1.2/2 GW) and matched hybrid (1/2 GW)
with three stresses. The matched hybrid is an equal-rating modelling control.
All three stresses start from Central demand and use the same weather,
solver settings and other system inputs.

| Tasks | Stress | Intervention |
|---|---|---|
| 0–2 | More demand in DK2 | Add 5 TWh/year as a flat fixed load at the DK2 receiving bus. |
| 3–5 | More demand in Germany | Add the same 5 TWh/year and flat profile at the German receiving bus. |
| 6–8 | More wind competing with Bornholm | Raise native offshore capacity minima in DK2 and Germany by 1 GW per region above the original solved Central split-radial portfolio. Allocate the increment to eligible sites in order of correlation with the Bornholm availability profile. |

The offshore stress uses existing native profiles, potential bounds and
investment costs. It commits the reference portfolio in the two regions
and the additional capacity on top of it. Other investment and dispatch
can respond. Germany’s modelled offshore resources are aggregated at
national resolution. The doses are controlled sensitivities.

Run from the campaign checkout on Sophia. Keep the original Central hybrid
and split-radial composed inputs, their prepared records, and the original
Central split-radial solved network and summary. The input check verifies
the original hashes, the reference solution, capacity headroom and the
paired demand and snapshot weights before starting the solves.

```bash
check_id=$(sbatch --parsable --array=0 --time=00:30:00 \
  scripts/bornholm_commitment_stress.sbatch --check) &&
sbatch --dependency=afterok:"${check_id%%;*}" --kill-on-invalid-dep=yes \
  scripts/bornholm_commitment_stress.sbatch
```

The check does not start a solver. At most six solves can run concurrently.
Logs are written to `slurm-bornholm-commitment-<job>_<task>.out`.
Transfer outputs are under `results/bornholm_2035_commitment_stress/<job>/transfer/`.
Each case includes its candidate network, summary and prepared record. An
optimal candidate with numerical screening flags is labelled
`optimal_needs_review` and its Slurm task exits with code 3. Inspect those
flags before interpreting small differences.

For the joint funding comparison, pair the three layouts *within* each
stress. Compare each hybrid against split radial, and use the difference
between the two hybrids to measure the value of the extra Danish rating.
Test a common funding rule against the original demand cases and these
stress cases after the physical solves.

## Central follow-up

The follow-up uses the prepared final-v1 Central treaty case as its common
input. It changes the required Bornholm wind and cable capacities. Wind
availability, cable losses and per-MW costs retain the treaty assumptions.
Removing a landing removes both directed links and adjusts the project's
allowance in any native DC expansion budget. The one-country layouts are
physical counterfactuals at the Bornholm site. They do not assume a particular
national ownership or CfD agreement. The follow-up uses the original final-v1
Gurobi barrier preset with 16 threads, `NumericFocus=1` and `BarConvTol=1e-8`.

| Array task | Layout | DK / DE cable GW | Wind GW |
| --- | --- | --- | --- |
| 0 | Matched hybrid | 1.0 / 2.0 | 3.0 |
| 1 | Germany alone | 0 / 2.0 | 2.0 |
| 2 | Denmark alone | 1.2 / 0 | 1.2 |
| 3 | Smaller island | 1.2 / 2.0 | 2.0 |

Reuse the Central no-island, split radial, treaty hybrid and more-DE results
already available. Compare the smaller island directly with the treaty hybrid
to assess the incremental third GW at identical cable ratings. On Sophia,
from the Bornholm checkout:

```bash
git pull --ff-only
sbatch --array=0-3 scripts/bornholm_matched.sbatch --central-followup
```

For an input-only run, use `--array=0` and pass
`--central-followup --check` to the script. The prepared input is
`resources/bornholm_2035_final_v1/central_hybrid/networks/composed_2035.nc`
and its run-level `prepared.json` must be available. New candidates and
summaries are under `results/bornholm_2035_central_followup/<array-job>/`.

## Historical fixed-capacity comparison

Use the full-year electricity-only hindcasts for 2024 and 2025. The runner
checks complete hourly coverage (8,784 hours for 2024 and 8,760 for 2025),
fixed capacities and dynamic fuel costs. Both cases use the same historical
demand, fleet, grid and operating constraints. They optimise dispatch and
flows separately. The project case adds 3 GW of wind, a 1.2 GW cable to DK2
and a 2 GW cable to Germany. The wind profile uses capacity-weighted
availability from the existing DK2 offshore fleet in the same weather year
as a proxy for the project site. Cable losses and annualised project costs
come from the prepared 2035 Central treaty template.

From the same Bornholm checkout on Sophia:

```bash
for year in 2024 2025; do
  historical="$HOME/projects/pypsa-eur-hindcast/results/hindcast-full-$year/network.nc"
  if [[ -s "$historical" ]]; then
    sbatch --job-name="bornholm-hindcast-$year" \
      scripts/bornholm_historical.sbatch "$historical" --year "$year"
  else
    printf 'Waiting for hindcast output %s\n' "$historical"
  fi
done
```

For 2023, use `results/hindcast-reviewed-2023/network.nc` in the hindcast
checkout and pass `--year 2023`. The historical job also accepts `--check` for
input validation without dispatch. Results are written to
`results/bornholm_historical/<year>/<job>/`. Examine operating savings and load
shedding before comparing against annualised project cost. The summary shows
operating savings in the original price basis and an indicative net figure in
EUR2025. The historical operating-cost conversion uses [Eurostat euro-area
annual average HICP](https://ec.europa.eu/eurostat/statistics-explained/SEPDF/cache/4176.pdf)
(2.4% in 2024 and 2.1% in 2025). Use `--operating-to-eur2025` if the source
requires another conversion. This is a one-year counterfactual with the
historical fleet fixed.

## Interpretation of demand sensitivities

In final-v1, the treaty hybrid's incremental cost against the corresponding
no-island case rises from €165.7m/year (Central) to €201.3m/year (EV65).
Project costs are nearly constant and savings elsewhere fall from about
€1,111m to €1,075m/year. EV55 and EV65 jointly change EV uptake and available
flexibility, industrial conversion, and fixed data-centre demand. At Central
demand, the separate V2G test lowered the Bornholm premium when participation
increased. The bundled results leave the contribution of each driver open.

The earlier +20% DK/DE fixed-load sensitivity increased total costs and reduced
the Bornholm premium from €196.6m to €135.1m/year. In that run the hybrid
displaced about 3.49 GW more German onshore wind capacity than it did at the
starting demand. The slower-onshore ceiling remained above both solutions.
This avoided investment depends on how much alternative capacity would still
be open to change. Keep the fixed-fleet operating comparison alongside the
endogenous 2035 investment result.
