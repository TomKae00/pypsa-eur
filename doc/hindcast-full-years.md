# Full-year 2024 and 2025 hindcasts

These electricity-only baselines hold background capacities fixed. The 2024
run includes 29 February and covers 8,784 hours. The 2025 run covers 8,760 hours.
Each run uses its annual load and ERA5 weather. The grid uses native PyPSA-Eur
topology and limits. The country-level nuclear factors hold availability
constant across each year.

## Prepare and solve the baselines

Run in the project environment with Gurobi available. Each command reuses
existing weather and market inputs, prepares missing hourly fuel and CO₂ files,
and writes a solved network under `results/hindcast-full-YEAR/`.

```bash
pixi run --as-is bash scripts/run_hindcast_full_year.sh 2024 16
pixi run --as-is bash scripts/run_hindcast_full_year.sh 2025 16
```

The 2024 output path is distinct from the earlier corrected 2024 archive.

### Submit on Sophia

From the `feature/historical-hindcast` checkout on the login node, submit the
two years in sequence so their Snakemake builds share the data cache safely.
The 2025 job can start once the 2024 job exits, regardless of its outcome.

```bash
(
  set -e
  git pull --ff-only origin feature/historical-hindcast
  revision=$(git rev-parse HEAD)
  options=(--cpus-per-task=16 --exclusive --mem=0 --time=48:00:00 --output=slurm-%x-%j.out)
  job_2024=$(sbatch --parsable "${options[@]}" --job-name=hindcast-full-2024 \
    scripts/hindcast_sophia.sbatch full-year 2024 "$revision")
  printf '2024 job %s\n' "${job_2024%%;*}"
  sbatch "${options[@]}" --job-name=hindcast-full-2025 --dependency="afterany:${job_2024%%;*}" \
    scripts/hindcast_sophia.sbatch full-year 2025 "$revision"
)
```

## Inspect the baselines

Set `YEAR` in [inspect_hindcast.ipynb](../notebooks/inspect_hindcast.ipynb)
to 2024 or 2025. The notebook uses the solved NetCDF, the mapping at
`data/hindcast/bus_zones_YEAR.csv`, and observed day-ahead prices.
Fetch prices for each year whose `data/hindcast/prices_zonal_YEAR.csv` is missing.

```bash
pixi run python -m scripts.hindcast fetch-prices --years 2024
pixi run python -m scripts.hindcast fetch-prices --years 2025
```

Compare hourly prices and spreads around DK2 before interpreting the
Bornholm capture price. For observed border flows, the benchmark can use
the output of `scripts/retrieve_hindcast_flows.py`.

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
