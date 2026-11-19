# Bornholm welfare and project design assessment

## Research question

Which combination of grid design and wind-support terms creates the most value
while giving Denmark and Germany a fair share given their alternatives,
project-specific contributions and risks? How much is it worth to consider grid
design and support terms jointly before either country makes a binding grid
commitment?

Bornholm is the worked case for the special course. The master's thesis uses the
same method to compare portfolios of Baltic hubs, with interacting routes,
market effects and investment dates. [Decision record and thesis design](https://chatgpt.com/library?file_id=libfile_fa4966d695588191bf8706bc2fa14188&account_id=personal)

## 1. Define the decision and validate the networks

| Decision | Reference | Comparison |
|---|---|---|
| Should we build the island? | No island, with the system free to build alternative assets under the same demand, climate and adequacy requirements | Candidate hub versus an adapting system |
| How should we connect the same 3 GW at Bornholm? | Split radial with 1 GW to Denmark and 2 GW to Germany | Both hybrids against split radial |
| What is the extra Danish capacity worth? | Hybrid with 1 GW to Denmark and 2 GW to Germany | Treaty hybrid with 1.2 GW to Denmark and 2 GW to Germany versus the matched hybrid |

The three Sophia tasks implement the connection and sizing comparisons. Verify
solver status, objective reconstruction and welfare reconciliation. Keep wind
profiles, wind costs, demand, emissions limits, transfers and background
investment assumptions comparable. The routing and sizing differences must add
up to the treaty-hybrid versus split-radial difference. Identify any change in
market-zone exposure when interpreting price and welfare differences. The 1/2
GW hybrid is a modelling control pending confirmation of engineering feasibility.
[Campaign instructions and transfer layout](../bornholm_2035_campaign.md#central-matched-cable-experiment)

Keep four balances visible in the welfare ledger: Denmark, Germany, the wind
developer and the rest of the modelled system. Each country includes consumer
and producer surplus, network income, investment costs, and public transfers.
Record the developer's ownership separately if it affects national attribution.
Deduct resource costs and count support, congestion income and other transfers
exactly once. Physical gains and changes in who pays should reconcile across
the chosen accounting boundary.

## 2. Screen strikes and funding shares

**Implemented for the central matched campaign.** From the repository root,
use the campaign Python environment with PyPSA and Matplotlib. Supply a complete
`transfer` directory for one matched Slurm array.

```bash
pixi run --as-is python -m scripts.assess_bornholm_agreement \
  /path/to/transfer --output results/bornholm_assessment/central
```

| Output | What to inspect |
|---|---|
| `matched_runs.csv` | Solver flags, source hashes, total modelled project cable costs and rents, generation and capture prices |
| `welfare_changes.csv` | The existing paired welfare and reconciliation checks at the four illustrative strikes |
| `agreement_screen.csv` | Strike limits with a fixed 70% German share and with a negotiable share, the surplus available over the three participation floors, and a selected agreement where one exists |
| `strike_and_funding_share.png` | Feasible terms by hybrid design and cable cost allocation, with the 70% German share and the selected agreement marked |

Use `--allow-provisional` to diagnose numerically flagged solves or an annual
welfare reconciliation residual above €10,000. The table and figure label those
outputs provisional.

### Settlement and fairness

The split-radial reference keeps its illustrative €75/MWh strike and 70% German
share during the candidate sweep. For the candidate, the two-sided annual
actual-output proxy pays `generation × (strike − capture price)`. Use
TWh/year and EUR/MWh to get million real 2025 EUR/year. Each country pays its
share of positive support or receives its share of negative support. Changes
in support affect the government and developer accounts by equal and opposite
amounts. The sweep holds dispatch, investment and market prices fixed. The
actual capability-based terms require a separate implementation when the
payment formula is available.

Each country's gain must reach its participation floor and its positive net
CfD payment must fit any annual budget ceiling supplied. The developer must
meet its own floor. By default, all three match the financial outcome of the
split-radial project under the modelled reference contract. Treat the developer
floor as a stand-in until we can calculate revenues, costs and financing over
the whole project and establish whether radial investment is a genuine option.
Run a separate financeability test when that option is unavailable. The flags
`--dk-floor`, `--de-floor`, `--developer-floor`, `--dk-cap` and `--de-cap` use
million EUR/year. Annual payment ceilings are model inputs. Translate a
multi-year fiscal ceiling using the contract-year settlement and price basis.

At each strike, total CfD support fixes the sum of the countries' payments.
The national welfare floors and budgets put an upper bound on each payment.
Those limits determine the feasible range of German shares, including years
with repayments. A country can still fail the test even if the available
three-party surplus is positive. The available contract terms can restrict
who pays and who receives a credit.

The selected benchmark uses the minimum strike that meets the developer floor
and then maximises the smaller of the countries' gains above their floors. The
default gives the gains equal weight. Use `--dk-weight` for a declared
project-specific sensitivity. A weight of 2 targets twice the extra gain for
Denmark where the limits permit it, with Germany's weight fixed at 1. Agree and
document the weights before comparing designs. Base them on evidenced project
impacts or risks that merit different sharing. Report national gains per MWh of
reference electricity demand as a measure of relative impact. Treat local host
impacts and cancellation risk explicitly where feasible.

The current code holds cable congestion income at equal shares and shows two
cable cost sensitivities, equal pooling and assignment to each landing country.
The asset model bundles some shared equipment within the legs. Replace the
sensitivity with an itemised allocation when evidence supports it. Track common
assets, onshore reinforcement, grants, target credits and third-country impacts
separately before drawing a distributional conclusion.

### Check against the original rounded results

The inputs below come from the earlier executed local-results and national
welfare notebooks. They use pooled cable costs and a split-radial reference.

| 2035 scenario | Surplus above all three floors, €m/year | Feasible strike with 70% German funding | Feasible strike with variable funding |
|---|---:|---|---|
| Central | +5.54 | None | €74.29–74.68/MWh |
| EV55 | −4.27 | None | None |
| EV65 | −9.70 | None | None |

For central demand, the equal-gain benchmark gives each country **€2.77m/year**
at about **€74.29/MWh** and a **63.8% German share**, with the developer at its
assumed split-radial return. These results use the old candidates and the
annual proxy. Use the validated matched outputs for new central results.

## 3. Test one agreement across futures

Extend the screen once the central accounting closes. Use equivalent designs
and references for EV55, EV65, selected weather years, cost and delay scenarios.
Rerun dispatch for physical changes, such as new offshore build lowering
capture prices. A price haircut alone tests financial exposure at fixed output.
V2G results need a matching split-radial reference before they support that
comparison. Keep no island as a separate investment comparison.

For each design and each strike, calculate the feasible cost-share interval in
each future and intersect the intervals. An empty intersection means that
strike cannot work across all those futures. Find the range of common strikes
and report the scenario and participation or budget condition that sets each
boundary. If no pair of terms passes, show the minimum shortfall and which
countries bear it. Test alternative risk allocation only through rules that
can be written into the agreement in advance. Use multi-year physical outcomes
and cash flows before assessing lifetime bankability or the full public cap.

## 4. Decide when the countries should commit

Compare three processes in a small decision tree using the same equipment,
auction and market information at the same dates.

| Process | Decisions the model can still change |
|---|---|
| Grid first | Commit to a grid package early. Negotiate the best feasible wind-support terms and cost split for that grid once wind-market evidence improves. |
| Joint assessment before grid commitment | Compare the grid and support together. Keep a conditional supplier slot only if the terms were available. Set the technical design before binding bids and let the auction set the strike. |
| Wind terms first | Wait for credible wind costs or bids before placing a binding grid order. Use the delivery slot and construction timetable that would then be available. |

Represent each equipment reservation, order and wind award by its decision
date, available information, fees, remaining freedom to revise, and
cancellation exposure. Include an outside option to defer or stop. Use the
welfare screen to decide which packages can be agreed at each stage. Value
each viable outcome in discounted system terms, including project timing,
avoidable resource expenditure and costs of maintaining options. Show annual
and discounted results for each country, including exposure to delay,
cancellation and changing settlement costs.

When scenario probabilities can be supported, rank the feasible processes by
expected discounted system value subject to the participation and risk
conditions. Otherwise report the thresholds at which the ranking changes,
such as the maximum reservation fee or delay that the value of flexibility can
justify. Calculate the value of joint assessment as the best feasible discounted
outcome under that process minus the best feasible outcome with the grid
already fixed. Test it with a common delivery date first. Then allow each
process to use its evidenced delivery dates and include reservation fees and
cancellation costs. Ground those trade-offs in dated supplier and
procurement evidence. The PyPSA dispatch solutions supply the physical welfare
and price inputs for this calculation.

## Course handover and thesis boundary

The course should leave a reproducible chain from matched networks to
reconciled ledgers, feasible CfD terms and a short decision note for Bornholm.
Report system value, each country's gain, the developer hurdle and downside
against the relevant reference. State which decisions the evidence supports
and which depend on missing engineering, procurement or finance data.

For the thesis, evaluate candidate hubs jointly within the Baltic network.
Make location, wind capacity, cable rating, landing country and build date
explicit choices. Recompute the alternatives available to each country and
the value of combinations of hubs as the network expands. Reuse the Bornholm
welfare, participation and staged-commitment tests for these portfolios.

## Analysis of the original 15 runs

The [executed final-v1 notebook](bornholm_final_v1_analysis.ipynb) uses the saved
rounded tables to compare all 15 runs and evaluate the same central CfD terms
in EV55 and EV65. It records the source tables, screening flags and the eight
available national welfare comparisons.
