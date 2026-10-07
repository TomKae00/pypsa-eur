# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
# SPDX-License-Identifier: MIT

"""Export the real bus/profile IDs needed to configure an energy-island case.

Run on the COMPOSED BASELINE for the same clustering and horizon as the case.
This inventory does not infer bidding zones or establish resource overlap.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa


def distance_km(x, y, x0, y0):
    lon, lat, lon0, lat0 = [np.radians(v) for v in (x, y, x0, y0)]
    a = (
        np.sin((lat - lat0) / 2) ** 2
        + np.cos(lat) * np.cos(lat0) * np.sin((lon - lon0) / 2) ** 2
    )
    return 6371 * 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def inventory(n, links, wind):
    hubs = set(wind.hub_id) | set(links.hub_id)
    if hubs.intersection(n.buses.index):
        raise ValueError(
            "Use the composed BASELINE: this network already contains an energy-island hub"
        )
    if isinstance(n.snapshots, pd.MultiIndex):
        raise ValueError("This first inventory supports a single overnight horizon")
    candidates = []
    for row in links.itertuples():
        buses = n.buses.loc[
            n.buses.carrier.eq("AC") & n.buses.country.eq(row.target_country)
        ].copy()
        columns = ["country", "x", "y"] + [
            c for c in ["bidding_zone", "zone", "location"] if c in buses
        ]
        buses = buses[columns]
        buses.index.name = "bus_id"
        buses["link_id"] = row.link_id
        buses["distance_to_landing_km"] = distance_km(
            buses.x, buses.y, row.target_x, row.target_y
        )
        candidates.append(buses.reset_index())
    buses = pd.concat(candidates, ignore_index=True).sort_values(
        ["link_id", "distance_to_landing_km"]
    )

    generators = n.generators.loc[n.generators.carrier.isin(wind.carrier)].copy()
    generators.index.name = "generator_id"
    generators["country"] = generators.bus.map(n.buses.country)
    generators["bus_x"] = generators.bus.map(n.buses.x)
    generators["bus_y"] = generators.bus.map(n.buses.y)
    floor = generators[["p_nom", "p_nom_min"]].max(axis=1)
    generators["reducible_capacity_mw"] = (generators.p_nom_max - floor).clip(lower=0)
    generators.loc[
        ~generators.p_nom_extendable | ~np.isfinite(generators.p_nom_max),
        "reducible_capacity_mw",
    ] = 0.0
    weights = n.snapshot_weightings.generators
    if not np.isfinite(weights).all() or weights.sum() <= 0:
        raise ValueError("Invalid snapshot weights")
    generators["has_valid_dynamic_profile"] = False
    generators["weighted_mean_availability"] = np.nan
    for name in generators.index.intersection(n.generators_t.p_max_pu.columns):
        profile = n.generators_t.p_max_pu[name].reindex(n.snapshots)
        valid = bool(np.isfinite(profile).all() and profile.between(0, 1).all())
        generators.at[name, "has_valid_dynamic_profile"] = valid
        if valid:
            generators.at[name, "weighted_mean_availability"] = (
                profile * weights
            ).sum() / weights.sum()
    columns = [
        "bus",
        "country",
        "carrier",
        "bus_x",
        "bus_y",
        "p_nom",
        "p_nom_min",
        "p_nom_max",
        "p_nom_extendable",
        "reducible_capacity_mw",
        "has_valid_dynamic_profile",
        "weighted_mean_availability",
    ]
    return buses, generators[columns].reset_index()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "network", help="Composed baseline .nc file, before island additions"
    )
    parser.add_argument("--links", default="data/energy_islands/bornholm_links.csv")
    parser.add_argument("--wind", default="data/energy_islands/bornholm_wind.csv")
    parser.add_argument("--output", default="results/bornholm_input_inventory")
    args = parser.parse_args()
    n = pypsa.Network(args.network)
    buses, generators = inventory(n, pd.read_csv(args.links), pd.read_csv(args.wind))
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    buses.to_csv(out / "landing_bus_candidates.csv", index=False)
    generators.to_csv(out / "wind_profile_and_potential_candidates.csv", index=False)
    print("Candidate buses (distance is to the clustered bus, not cable route length):")
    print(buses.groupby("link_id", sort=False).head(5).to_string(index=False))
    print(
        "\nOffshore generators (country refers to connected bus; inspect resource geography separately):"
    )
    print(generators.to_string(index=False))
    print(f"\nCSV inventories written to {out}")
    print(
        "Choose the receiving zone explicitly. Do not select a potential donor just because its bus is nearby."
    )


if __name__ == "__main__":
    main()
