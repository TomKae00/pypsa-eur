# SPDX-License-Identifier: MIT
"""Check that the isolated simplified Bornholm baseline retains study essentials."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa

SIMPLIFIED = "bornholm_simplified_baseline"
REFERENCE = "bornholm_zonal_baseline"
HORIZON = 2035


def require(condition, message):
    if not condition:
        raise ValueError(message)


def count_carrier(network, component, carrier):
    return int(getattr(network, component).carrier.eq(carrier).sum())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("composed", "solved"), default="composed")
    args = parser.parse_args()

    root = "resources" if args.stage == "composed" else "results"
    target = Path(root) / SIMPLIFIED / "networks" / f"{args.stage}_{HORIZON}.nc"
    reference = Path("resources") / REFERENCE / "networks" / f"composed_{HORIZON}.nc"
    require(target.is_file(), f"Missing simplified network: {target}")
    require(reference.is_file(), f"Missing reference composed network: {reference}")

    n = pypsa.Network(str(target))
    b = pypsa.Network(str(reference))

    require(n.snapshots.equals(b.snapshots), "Time snapshots differ from reference")
    require(len(n.snapshots) == 8760, "Expected 8760 hourly snapshots")
    require(
        np.allclose(n.snapshot_weightings, b.snapshot_weightings),
        "Snapshot weights changed",
    )
    ac = n.buses.index[n.buses.carrier.eq("AC")]
    original_ac = b.buses.index[b.buses.carrier.eq("AC")]
    require(ac.equals(original_ac) and len(ac) == 52, "Electricity zones changed")
    require("BEI_HUB" not in n.buses.index, "Baseline contains Bornholm project hub")

    sector = n.meta.get("sector", {})
    expected = {
        "co2_spatial": False,
        "co2_network": False,
        "methanation": False,
        "hydrogen_fuel_cell": False,
        "solar_thermal": False,
        "v2g": True,
        "hydrogen_turbine": True,
        "tes": True,
        "gas_network": True,
        "gas_network_connectivity_upgrade": 1,
        "H2_network": False,
        "H2_retrofit": False,
    }
    for key, value in expected.items():
        require(sector.get(key) == value, f"Unexpected sector.{key}: {sector.get(key)}")
    require(
        sector.get("regional_co2_sequestration_potential", {}).get("enable") is False,
        "Regional sequestration potential is still enabled",
    )
    # NetCDF metadata round-tripping stringifies YAML integer dictionary keys.
    # Accept both representations, but still verify the actual 2035 cap.
    sequestration_limits = sector.get("co2_sequestration_potential", {})
    limit = sequestration_limits.get(str(HORIZON), sequestration_limits.get(HORIZON))
    require(
        limit is not None and np.isclose(float(limit), 100.0),
        f"Unexpected global CO2 sequestration limit in network metadata: "
        f"2035={limit!r}; available values={sequestration_limits!r}",
    )
    carrier = n.links.carrier
    require(not carrier.eq("CO2 pipeline").any(), "CO2 pipelines remain")
    require(
        not carrier.str.contains("H2 pipeline", na=False).any(), "H2 pipelines remain"
    )
    require(not carrier.eq("H2 Fuel Cell").any(), "Hydrogen fuel cells remain")
    require(not carrier.eq("Sabatier").any(), "Methanation remains")
    require(
        not n.generators.carrier.str.contains(
            "solar thermal", case=False, na=False
        ).any(),
        "Solar thermal remains",
    )
    require(carrier.eq("V2G").any(), "Vehicle-to-grid disappeared")
    require(carrier.eq("H2 turbine").any(), "Hydrogen turbine disappeared")
    require(carrier.eq("gas pipeline").any(), "Gas transport disappeared")
    require(n.buses.carrier.eq("H2").any(), "Hydrogen buses disappeared")
    require(n.buses.carrier.eq("co2").any(), "Carbon accounting bus disappeared")
    require(
        n.stores.carrier.eq("co2 sequestered").any(),
        "Geological CO2 storage disappeared",
    )

    gas_stores = n.stores.index[n.stores.carrier.eq("gas")]
    reference_gas_stores = b.stores.index[b.stores.carrier.eq("gas")]
    pd.testing.assert_frame_equal(
        n.stores.loc[
            gas_stores,
            ["bus", "e_nom_min", "e_nom_max", "e_nom_extendable", "capital_cost"],
        ].sort_index(),
        b.stores.loc[
            reference_gas_stores,
            ["bus", "e_nom_min", "e_nom_max", "e_nom_extendable", "capital_cost"],
        ].sort_index(),
        check_dtype=False,
        obj="shared gas-storage inputs (rebuild both networks after storage fix)",
    )

    electricity = n.loads.index[n.loads.carrier.eq("electricity")]
    ref_electricity = b.loads.index[b.loads.carrier.eq("electricity")]
    require(electricity.equals(ref_electricity), "Fixed electricity load set changed")
    pd.testing.assert_frame_equal(
        n.loads_t.p_set.reindex(columns=electricity),
        b.loads_t.p_set.reindex(columns=ref_electricity),
        check_dtype=False,
    )
    if args.stage == "solved":
        require(
            n.meta.get("solve_result")
            == {"status": "ok", "termination_condition": "optimal"},
            "Solver did not report optimal termination",
        )
        require(np.isfinite(n.objective), "Solved objective is not finite")
        require(
            np.isfinite(n.buses_t.marginal_price.loc[:, ac]).all().all(),
            "Missing finite electricity prices",
        )

    print(f"Simplified baseline check passed ({args.stage}):")
    print(f"  snapshots: {len(n.snapshots)}; electricity zones: {len(ac)}")
    for comp in ("buses", "generators", "links", "stores", "storage_units", "loads"):
        old, new = len(getattr(b, comp)), len(getattr(n, comp))
        print(f"  {comp}: {old} -> {new} ({new - old:+d})")
    print("  global 2035 CO2 storage limit: 100 MtCO2; transport: pooled")
    print(
        "  WARNING: pooled CO2 also changes the industrial capture electricity representation."
    )


if __name__ == "__main__":
    main()
