# Full-year 2024 and 2025 hindcasts

These electricity-only baselines hold background capacities fixed. The 2024
run includes 29 February and covers 8,784 hours. The 2025 run covers 8,760 hours.
Each run uses its annual load and ERA5 weather. Nuclear availability is held
constant by country within each year.

## Reviewed grid

The normal full-year launcher applies 29 directional corridor profiles and
dated asset rules. It can reoptimise an existing native composed or solved
network with exactly the same fleet, load, fuels and weather. If neither native
input is available, it composes the annual network with the grid overlay.
Solved results go to `results/hindcast-reviewed-YEAR/network.nc`. To replay
the native-grid version, pass `native` as the third argument to the launcher.
That version uses the existing `results/hindcast-full-YEAR/` directory.

- **2024 through 29 October local delivery.** Use the Nordic operators' maximum
  directional NTCs as fixed upper envelopes. The FI–NO4 line uses the operating
  agreement's physical TTC upper bounds because trade was arranged via Sweden.
- **From 30 October 2024 local delivery.** Use the published JAO Max Border Flow
  by direction as hourly upper envelopes. The interface covers both historical
  years after the switch. From the October 2025 quarter-hour change, each input
  hour uses the arithmetic mean of its four quarter-hour ceilings. This
  preserves the corresponding upper limit on hourly energy.
- **Individual assets.** The 2024 Shetland link stays unavailable through the
  reported energisation day. Fenno-Skan 1 retains its installed nameplate and
  uses the 400 MW operating limit. The first Italy–Montenegro pole is 600 MW.
  The Finnish–Swedish operating limits reflect the change around the November
  2025 Aurora Line commissioning. Existing physical branch limits still apply.

JAO publishes single-border extrema of the flow-based domain. These bounds
retain the observed time variation and leave the original PyPSA AC and DC
constraints in place. The model can reach combinations of border flows that
the joint market domain would exclude. NO4–FI, the two-way SE3–SE4 aggregate
and the external interconnectors use the reported real-zone border fields.
The separately reported virtual South-West Link fields are excluded from the
sum. Unreviewed borders use their physical branch constraints. For example,
the JAO DK2–Germany operating number covers the market interface while the
historical model can still be limited by the physical representation of that
interface. Assets missing from the zonal network and outages outside these
source series remain outside the historical grid inputs.

Sources: [Nordic maximum NTC map](https://www.nordpoolgroup.com/globalassets/download-center/tso/max-ntc.pdf),
[Fingrid/Statnett operating agreement](https://www.fingrid.fi/globalassets/dokumentit/fi/kantaverkko/sahkonsiirto/appendix-1-9-to-soa-annex-os-v3.pdf),
[Nordic RCC go-live](https://nordic-rcc.net/nordic-flow-based-market-coupling-project-confirms-go-live-for-29th-october-2024/),
[JAO Publication Tool](https://publicationtool.jao.eu/nordic/) and
[Fingrid on Aurora Line](https://www.fingrid.fi/en/news/news/2025/aurora-line-power-transmission-connection-to-be-commissioned-on-november-13).

The compressed inputs are committed with their source response hashes.
To refresh the data from the published API, use
`python -m scripts.prepare_hindcast_transmission --cache data/hindcast/jao-cache`.
The script checks coverage, duplicates and missing directions over both calendar
years before preparing the constraints.

## Prepare and solve

Run in the project environment with Gurobi available.

```bash
pixi run --as-is bash scripts/run_hindcast_full_year.sh 2024 16
pixi run --as-is bash scripts/run_hindcast_full_year.sh 2025 16
```

### Submit on Sophia

From the checkout containing this update, submit the two years in sequence so
builds can share the existing data cache. The 2025 job starts once the 2024 job
exits, regardless of its outcome.

```bash
(
  set -e
  revision=$(git rev-parse HEAD)
  options=(--cpus-per-task=16 --exclusive --mem=0 --time=48:00:00 --output=slurm-%x-%j.out)
  job_2024=$(sbatch --parsable "${options[@]}" --job-name=hindcast-reviewed-2024 \
    scripts/hindcast_sophia.sbatch full-year 2024 "$revision")
  printf '2024 job %s\n' "${job_2024%%;*}"
  sbatch "${options[@]}" --job-name=hindcast-reviewed-2025 --dependency="afterany:${job_2024%%;*}" \
    scripts/hindcast_sophia.sbatch full-year 2025 "$revision"
)
```

## Benchmark and plot

```bash
pixi run python -m scripts.hindcast benchmark --year 2025 \
  --config config/hindcast/benchmark-full-years.yaml \
  --network results/hindcast-reviewed-2025/network.nc \
  --output results/hindcast-reviewed-benchmark-2025
```

Fetch the observed price file with `python -m scripts.hindcast fetch-prices --years 2024 2025`
if it is missing. Use `scripts/retrieve_hindcast_flows.py` for an observed
border-flow comparison. For the 2025 and 2035 grid figures, pass the reviewed
historical NetCDF and the central Bornholm hybrid candidate on the other branch.
The plotter also accepts the matching prepared networks before solving.

```bash
pixi run python -m scripts.plot_bidding_zone_grid \
  --network-2025 results/hindcast-reviewed-2025/network.nc \
  --network-2035 /path/to/central_hybrid/candidate_2035.nc \
  --zones /path/to/bidding_zones.geojson \
  --output results/grid-capacity-comparison
```

The figure bundle displays the same zones with line widths, directly labelled
Nordic and Baltic corridors, and a difference map with both directions in a
table. A numeric CSV includes every represented cross-zone corridor and its
source. The historical summary is the 95th percentile of the hourly upper
bound after also applying physical ratings. The 2035 comparison uses the
fixed reviewed ceilings over the selected TYNDP 2026 grid. Width shows the
larger direction. The difference map includes changes in the source of the
transfer assumptions as well as physical investments.

## 2025 data provenance

- The February 2026 ENTSO-E and NESO archives supply national load.
  The native demand workflow fills gaps and allocates load within countries.
- The IRENASTAT 2026 H1 [source CSV](https://tubcloud.tu-berlin.de/s/qALFkR8eB73FWA9/download/IRENASTAT_capacities_2000-2025.csv)
  supplies 2025 PV, onshore and offshore capacity totals. The committed 2025
  table uses the country and technology conventions from powerplantmatching
  0.9.1, with capacities summed across grid connection categories.
  Source SHA-256 `bcf94df099d4d9da3382242a65a7268f6b39052635f6eed192e91b3460045bac`.
- The nuclear factors weight each reactor's reported 2025 annual energy
  availability by its reference unit MW. The calculation includes
  reported partial-year operation and uses the published 2025 reactor
  entries. Flamanville 3 has no 2025 factor in the report and is excluded
  from the French average. The capacity filter excludes Mochovce 4,
  which was under construction at the end of 2025.
  [IAEA operating experience 2026](https://www-pub.iaea.org/MTCD/publications/PDF/OPEX_2026_web.pdf)
  and [IAEA reactor status at year-end](https://www-pub.iaea.org/MTCD/publications/PDF/RDS-2_46web.pdf).
