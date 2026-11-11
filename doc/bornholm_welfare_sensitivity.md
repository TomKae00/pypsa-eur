# Bornholm welfare sensitivity campaign — 11 November 2026

This campaign studies baseline versus the **fixed 3 GW hybrid**, with the 1.2 GW
DK2 and 2 GW German export connections, in six matched electricity-system
environments. It retains the calibrated three-hourly 2013 weather data and the
reviewed 2035 transmission inputs. It does not repeat preprocessing or overwrite
campaign 4970565. The original **unsolved composed** baseline/hybrid networks are
copied, validated and modified under `bornholm_welfare_2035_v1`.

## Run exactly these 12 cases

| Array tasks | Onshore deployment ceilings | Fixed DK/DE electricity demand | Project pair |
|---|---|---|---|
| 0, 1 | Central | Current (100%) | Baseline, hybrid |
| 2, 3 | Central | +10% (110%) | Baseline, hybrid |
| 4, 5 | Central | +20% (120%) | Baseline, hybrid |
| 6, 7 | Slower onshore | Current (100%) | Baseline, hybrid |
| 8, 9 | Slower onshore | +10% (110%) | Baseline, hybrid |
| 10, 11 | Slower onshore | +20% (120%) | Baseline, hybrid |

There is **no lower-demand case**. Even task indices are baseline; odd indices
are fixed hybrid. The controller submits at most six simultaneous 16-thread
solves, so the 12 tasks can run in two waves subject to scheduler availability.
The 24-hour per-task walltime is a limit, not a runtime prediction.
Set `BORNHOLM_WELFARE_CONCURRENCY=12` in the submission environment to permit
all 12 cases to run simultaneously if the scheduler permits it.

## Onshore wind ceilings: total installed capacity, not additional build

| Region | Central ceiling, GW | Slower-onshore ceiling, GW | Basis |
|---|---:|---:|---|
| DK1 | 6.14418 | 5.74418 | AF26 primo-2035; slower case halves additions after primo-2030 |
| DK2 | 0.92852 | 0.82852 | Same method, retaining AF26's geographical allocation |
| Denmark total | 7.07270 | 6.57270 | Sum of the two bidding zones |
| Germany | 150.000 | 123.285714 | Central: published EEG-2035 pathway; slower case is a study interpolation |

**Published inputs:** [AF26 workbook](https://ens.dk/media/8833/download),
`Vindmøller!N198` (DK1), `N191` (DK2), with new-build stocks `I195`, `N195`,
`I188`, `N188`. AF26 reports capacity at the **start** of the year, so these are
primo-2035, not the workbook's separate end-2035 totals. Its onshore methodology
explicitly accounts for land-use, municipal approval and neighbours, with a
200 MW/year national commercial onshore build ceiling after 2030.
[AF26 land renewables note, pp. 5–6](https://ens.dk/media/8845/download).

**Study judgement:** the slower Danish case reduces this to 100 MW/year during
2030–2035. Retirements and the AF26 geographical split are retained. This removes
400 MW in DK1 and 100 MW in DK2 from the central 2035 ceilings. It is a stress
test, not another official forecast.

**Germany:** [BNetzA approval, April 2025, p.52](https://www.bundesnetzagentur.de/DE/Fachthemen/ElektrizitaetundGas/NEP/DL_Szenariorahmen/2025/Genehm_SR_2025Strom.pdf?__blob=publicationFile&v=2)
records a 150 GW onshore target for 2035 and a lower A2037 endpoint of 126.6 GW
(Table 16, p.53). The slower ceiling is explicitly constructed as
`115 + (5/7)*(126.6 - 115) = 123.285714 GW`, using the published 2030 target
of 115 GW and the lower 2037 endpoint. This is an analyst interpolation; it is
not an official A2035 value or an internally reproduced NEP scenario.

These are **upper bounds**, not forced deployment targets or measured social
acceptance limits. Existing fixed capacity and extendable minimum capacity count
towards each ceiling. All resource classes are summed. Inputs fail preflight
if a ceiling is below committed/minimum capacity; assets are not silently deleted.
German limits may not bind: the previous baseline has about 94 GW onshore, below
both ceilings. The particularly material correction is DK2, previously about
4.11 GW against the AF26 0.929 GW reference ceiling.

## Solar ceilings, identical in the two onshore cases

Solar is constrained to a documented central envelope, not independently varied
in this first campaign. Rooftop and ground-mounted assets are separate in Denmark.
Tracker variants count towards ground-mounted solar.

| Region and type | Published AF26 MW AC | Applied model-nameplate MW |
|---|---:|---:|
| DK1 rooftop | 2,306.83 | 2,768.196 |
| DK2 rooftop | 1,213.32 | 1,455.984 |
| DK1 ground, including tracker/battery categories | 14,832.01 | 17,798.412 |
| DK2 ground, including tracker/battery categories | 6,581.39 | 7,897.668 |
| Germany, all PV | — | 308,142.857 |

AF26 `Solceller!M6:M12` and `M15:M21` report primo-2035 **AC connection** capacity.
To avoid treating this as module-peak capacity, a **study DC/AC ratio of 1.2** is
applied to the Danish ceilings. It is a declared unit-harmonisation assumption,
not an AF26 value or a site-specific inverter model; it can be varied later.
The German module-nameplate ceiling is a study interpolation of the published
2030 PV target 215 GW and B2037 345.4 GW:
`215 + (5/7)*(345.4-215) = 308.142857 GW`. Germany is not split into a fabricated
rooftop/ground target. National upper bounds do not reproduce official build mixes.

The AF26 source cells, downloaded-workbook SHA256 and locators are retained in
`data/model_assumptions/bornholm_welfare_sources.csv`.

## Demand: fixed electricity end use, not total future gross electricity

The explicit study factors are **1.0, 1.1 and 1.2**, applied identically to baseline
and hybrid in each pair. They multiply the existing hourly profiles of fixed
electricity Loads on **AC and low-voltage buses in Denmark and Germany**, including
the static industry/agriculture electricity Loads. Hourly shape is retained.
Other countries, final heat demand, mobility/EV battery demand, hydrogen demand
and endogenous electricity use by converters are not multiplied.

The saved campaign-4970565 baseline gives this check of the affected boundary:

| Fixed electricity Loads | Current, TWh | +10%, TWh | +20%, TWh |
|---|---:|---:|---:|
| Denmark | 34.4271 | 37.8698 | 41.3125 |
| Germany | 527.4942 | 580.2436 | 632.9930 |

These are **experimental demand stresses**, not imported AF26/NEP total-demand
forecasts. Total electricity consumption, including heat pumps and other
conversion technologies, will not necessarily rise by exactly 10% or 20%.
AF26's classical and gross demand and BNetzA's total net electricity demand have
different sector boundaries from this model; imposing them blindly would
double-count some sector-coupled consumption. AF26 2035 classical demand is
33.54584 TWh net / 35.64636 TWh gross (workbook `Elforbrug!L31:L32`), supplied as
context rather than an imposed replacement. In every paired comparison all final
demand time series, snapshot weights and carbon/resource budgets must match.

## Wind CAPEX and CfD settlement after the system solves

No extra annual solve is needed for changing the CAPEX of this **fixed** 3 GW
wind asset. The collector produces 18 hybrid cost observations: six system
environments times these three cost cases:

| Case | Wind-farm investment, EUR2025/kW | Source |
|---|---:|---|
| Current | 2,039.704 | Native technology-data v0.15.0, 2035 `offwind` investment |
| Intermediate | 2,760 | Study midpoint, rounded |
| Higher | 3,480 | BCG Bornholm wind component, Figure 3-2, p.25 |

[Native costs](https://github.com/PyPSA/technology-data/blob/v0.15.0/outputs/costs_2035.csv);
[BCG report](https://ens.dk/media/7826/download). BCG's 3,480 excludes its separate
1,050 EUR/kW electrical component. The higher case does not reproduce BCG's
whole project or financial business case and is not a forecast upper bound.

Hold the native **absolute annual O&M**, collection cost, dispatch and 30-year
wind lifetime fixed. Apply the existing 7% discount rate only to the incremental
wind CAPEX. Do not automatically increase O&M through a constant percentage of
higher CAPEX. For each hybrid:

`C_wind(case) = C_wind(native) + 3000 MW * 1000 * annuity(30, 7%) * (CAPEX_case - 2039.704)`.

Divide this annual cost by actual wind generation for wind-and-collection LCOE.
Exclude both mainland export links from that LCOE, while retaining their costs
in complete-project and national welfare accounts. Collection costs remain in
wind costs. The collector separately exports market revenue, available and
generated energy, and the annual net support needed to cover model costs.

The later CfD grid is 60, 65, 70, 75, 80, 90, 100, 110, 120 and 130 EUR2025/MWh,
with a contract-specific equal-remuneration strike. A bankable 20-year contract
needs separate lifetime cash-flow/financing work; annual cost recovery is not a
tender bid. Keep support funding DK30/DE70 and new-project congestion income
50/50 in the main assessment. Developer ownership remains a sensitivity, not
implicitly a fact that every receipt belongs to Denmark. Government-financed
mainland connections do not become part of the developer's strike requirement.

## Other assumptions retained

- Biomass 30 EUR2025/MWh of fuel Europe-wide; native other fuel/technology costs.
- Reviewed directional transmission limits and fixed 2035 background grid.
- Existing sector boundaries, carbon budget and calibrated renewable profiles.
- Original Method=2, crossover=0, threads=16, NumericFocus=1, BarConvTol=1e-8
  preset; Gurobi SoftMemLimit capped by available physical RAM as before.
- Wind cost assumptions do not alter European background offshore costs.
- No new weather years or outages tonight; those are subsequent robustness tests.

Carbon-cap scarcity rents still need a justified national-account treatment before
definitive welfare conclusions. This campaign does not silently assign them.

## Sophia submission

From the existing Bornholm checkout after pulling the update:

```bash
cd ~/projects/pypsa-eur-bornholm &&
git checkout feature/bornholm-energy-island &&
git pull --ff-only &&
sbatch scripts/bornholm_welfare_workq.sbatch
```

The controller first validates the original two composed inputs and prepares all
12 new inputs. Missing original inputs cause an explanatory failure **before**
any solve submissions. In that case build them with the existing
`sbatch scripts/bornholm_calibrated_workq.sbatch build`, then submit this campaign.
It never starts another eight-case rebuild unnecessarily.

Numerical/physical screening thresholds are unchanged. If a solve is optimal but
fails numerical screens, its candidate and diagnostics are retained with
`optimal_needs_review` and Slurm exit code 3; it is not promoted to an accepted
network. The collector depends on **afterany**, so one flagged case cannot block
retrieval of the other results. A successful collector is not an acceptance
certificate: read `campaign_status.csv`.

Outputs:

- `results/bornholm_welfare_2035_v1/campaign_status.csv`
- `results/bornholm_welfare_2035_v1/wind_cost_cases.csv`
- `results/bornholm_welfare_2035_v1/transfer/` (candidate/config/summary per case)

Monitor with `squeue -r -u "$USER" -o "%.18i %.32j %.12T %R"`. For transfer,
use `rsync -avhL` so candidate-network symlinks are followed.

Preparation stores source-network hashes, code/config fingerprints, regional
capacity minima/ceilings and before/after demand totals. Changed prepared inputs
or changed queued code are rejected. Original outputs are preserved.

## Local verification

Meaningful small-network tests check summed capacity limits including existing
assets using HiGHS, all resource bins, low-voltage and static-load demand scaling,
preservation of other end uses, paired final demand, CAPEX-only cost changes,
preparation-failure submission guards, afterany collection and array index 11.
The highest-demand/slower-onshore transformation was also checked against the
actual saved baseline network. Full annual Gurobi execution and cluster scheduling
are checked by the Sophia jobs, not claimed as locally completed.

The existing timestamp hash check also now verifies exact unit-conversion
precision explicitly, avoiding a removed pandas keyword. Stored hashes of
unchanged timestamps are preserved; the actual saved network still validates.
