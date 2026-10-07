# Bornholm: first code patch

Base: `feature/bornholm-energy-island`, commit
`0194b59cb1e7b6244af1e70133bd82d5ef31b142`.

This patch corrects physical-input handling and adds a standalone input-inventory helper.
The helper is not part of the Snakemake model-building or solving workflow.
It is not a calibrated 2035/2040 scenario, a bidding-zone conversion, or a CfD model.
It intentionally requires the regional choices to be supplied before an island
case can compose. The baseline can still compose without those choices.

## Apply

From your repository root, after extracting this package somewhere outside the repository:

```bash
git apply --check /path/to/bornholm_physical_inputs.patch
git apply /path/to/bornholm_physical_inputs.patch
```

Run the second command only if the check succeeds. Do not overwrite conflicting
local changes. The package also contains complete replacement files under `files/`.
No changes have been pushed to GitHub.

## 1. Compose the baseline and inspect actual IDs

From the repository root:

```bash
pixi run snakemake resources/smoke_baseline/networks/composed_2035.nc --cores 4 --configfile config/config.bornholm.smoke.yaml
pixi run python scripts/inspect_energy_islands.py resources/smoke_baseline/networks/composed_2035.nc --links data/energy_islands/bornholm_links.smoke.csv --wind data/energy_islands/bornholm_wind.smoke.csv
```

The second command writes:

- `results/bornholm_input_inventory/landing_bus_candidates.csv`
- `results/bornholm_input_inventory/wind_profile_and_potential_candidates.csv`

Use a baseline COMPOSED network, not a solved Bornholm network whose potential
has already been reduced. If you changed `run.prefix`, use the corresponding
resource path. These commands do not require solving the European model.

## 2. Set receiving buses in the smoke links CSV

In `data/energy_islands/bornholm_links.smoke.csv`, fill `target_bus` for all three rows:

- `BEI_TO_DK2`: the bus representing the Danish receiving region.
- `BEI_TO_DE`: the German receiving bus.
- `BEI_TO_DK2_RADIAL`: the same Danish receiving bus as the hybrid.

Copy exact IDs from the inventory after reviewing their geography. The link
name `BEI_TO_DK2` is only a label. Country and distance alone do not prove DK2
membership. If the clustering merges relevant regions, change the clustering
before interpreting prices; do not solve that by renaming a link.

An explicit `length_km` is now treated as the full physical cable route and is
NOT multiplied by `length_factor` again. If it is blank, the code estimates the
route as hub-to-landing great-circle distance times `length_factor` (currently
1.25). This estimate is independent of the selected cluster's centroid and is
logged as an estimate. Actual engineering route lengths remain preferable.

## 3. Select the wind proxy explicitly

Fill `source_generator` in `data/energy_islands/bornholm_wind.smoke.csv` with
an exact generator ID from the inventory. It must have carrier `offwind-dc`
and a complete dynamic profile in [0, 1].

This makes the proxy choice visible; it does not create a site-specific Bornholm
profile. The source's connected bus is not the wind site's location. A dedicated
project-area profile and geometric resource exclusion remain a subsequent step.

## 4. Define the resource overlap rather than guessing it

Replace the empty mapping in `config/config.bornholm.smoke.yaml`. The following
is syntax only: replace the ID and MW after reviewing the regional resources.

```yaml
energy_islands:
  potential_allocation:
    BEI_WIND:
      "ACTUAL_GENERIC_GENERATOR_ID": 3000
  potential_allocation_notes:
    BEI_WIND: "Document which regional wind resource overlaps the project."
```

Edit the existing `energy_islands` section; do not add a second YAML section
with the same name. Multiple source IDs can be used, each with its own MW value.

For a fully represented 3 GW project the allocation sums to 3000 MW. If only
part of the project is represented in the existing resource pool, subtract only
that overlap and supply a note explaining it. An explicit empty project mapping
(`BEI_WIND: {}`) also requires a note explaining why no subtraction is needed,
such as a prior verified geographic exclusion.

Insufficient Danish potential is not evidence of partial overlap. The inventory
shows numerical bounds, not project-area overlap. Do not simply subtract all
remaining Danish potential, disable subtraction, or take the difference from
Sweden to get the model running. Revisit the resource geography when necessary.

The patch protects existing/minimum capacity, rejects foreign-bus donors in
this Danish bookkeeping approximation, and validates the complete allocation
before changing any bounds. Donors must currently use the same wind carrier.
It does not yet resolve overlapping potential across offshore technology types.

## 5. Recompose, then solve

```bash
pixi run snakemake resources/smoke_hybrid/networks/composed_2035.nc resources/smoke_radial/networks/composed_2035.nc --cores 4 --configfile config/config.bornholm.smoke.yaml
pixi run snakemake --cores 4 --configfile config/config.bornholm.smoke.yaml
```

Do not treat an old solved network as a result of the new inputs. These config
and script changes are tracked by the existing compose rule and require the
affected networks to be rebuilt.

The 90-cluster regular case uses the original `bornholm_links.csv` and
`bornholm_wind.csv`; configure those independently using its baseline inventory.
The new smoke copies prevent 20-cluster IDs from overwriting regular inputs.
Repeat the mapping review whenever clustering or the base network changes.

## Verification

```bash
pixi run python -m pytest --noconftest -q test/test_energy_islands.py
```

All 13 focused regression checks passed with PyPSA 1.2.4 and HiGHS, including a
small solved hybrid network, including a NetCDF round trip. A full European
Snakemake run and your local input mapping were not validated here.
The current patch has no analysis/CfD scripts and does not change your solver settings.

The next implementation steps are a dedicated project profile/resource mask,
an explicit hourly future-system/market configuration, and national financial
accounting. This patch leaves the optimisation objective and baseline/radial/
hybrid scenario definitions intact.
