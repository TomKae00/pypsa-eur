# SPDX-License-Identifier: MIT
"""
Build total hourly generator marginal costs from explicit historical inputs.

Fuel CSV: utc_timestamp,gas,coal,lignite,oil (EUR/MWh_th, same money basis as
observed prices). CO2 CSV: utc_timestamp,co2 (EUR/tCO2). VOM YAML: mapping from
thermal generator carrier to EUR/MWh_el. Every emitting generator needs a fuel
mapping and VOM; missing inputs cause an error. No fuel, CO2 or VOM is guessed.
Daily/weekly series must first be expanded to hourly with a documented policy.
"""

import argparse
from pathlib import Path

import numpy as np
import yaml
from hindcast import digest, hourly_network, read_hourly, write_json

FUEL_MAP = {
    "OCGT": "gas",
    "CCGT": "gas",
    "gas": "gas",
    "coal": "coal",
    "lignite": "lignite",
    "oil": "oil",
    "biomass": "biomass",
    "waste": "waste",
}


def build(n, fuel, co2, vom):
    idx = hourly_network(n)
    fuel, co2 = fuel.reindex(idx), co2.reindex(idx)
    result = n.get_switchable_as_dense("Generator", "marginal_cost").copy()
    result.index = idx
    for name, g in n.generators.iterrows():
        carrier = g.carrier
        intensity = n.carriers.at[carrier, "co2_emissions"]
        if carrier not in FUEL_MAP:
            if intensity > 0:
                raise ValueError(f"No fuel mapping for emitting carrier {carrier}")
            continue
        if carrier not in vom:
            raise ValueError(f"Explicit variable O&M required for {carrier}")
        if not np.isfinite(g.efficiency) or g.efficiency <= 0:
            raise ValueError(f"Invalid efficiency: {name}")
        if "efficiency" in n.generators_t and name in n.generators_t.efficiency:
            efficiency = n.generators_t.efficiency[name].set_axis(idx)
        else:
            efficiency = g.efficiency
        if not np.isfinite(efficiency).all() or (np.asarray(efficiency) <= 0).any():
            raise ValueError(f"Invalid hourly efficiency: {name}")
        if not np.isfinite(intensity):
            raise ValueError(f"Invalid CO2 intensity: {carrier}")
        result[name] = (
            fuel[FUEL_MAP[carrier]] + co2["co2"] * intensity
        ) / efficiency + vom[carrier]
    if not np.isfinite(result.to_numpy()).all():
        raise ValueError(
            "Missing or non-finite costs; check fuel/CO2 timestamps and efficiencies"
        )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("network")
    parser.add_argument("--fuel-prices", required=True)
    parser.add_argument("--co2-prices", required=True)
    parser.add_argument("--vom", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    import pypsa

    n = pypsa.Network(args.network)
    costs = build(
        n,
        read_hourly(args.fuel_prices),
        read_hourly(args.co2_prices),
        yaml.safe_load(Path(args.vom).read_text()),
    )
    output = Path(args.output)
    if output.exists():
        raise ValueError(f"{output} already exists")
    output.parent.mkdir(parents=True, exist_ok=True)
    costs.to_csv(output, index_label="utc_timestamp")
    write_json(
        output.with_suffix(".json"),
        {
            k: {"path": getattr(args, k), "sha256": digest(getattr(args, k))}
            for k in ["network", "fuel_prices", "co2_prices", "vom"]
        },
    )
    print(output)


if __name__ == "__main__":
    main()
