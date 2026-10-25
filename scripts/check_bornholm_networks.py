# SPDX-License-Identifier: MIT
"""Standalone checks of the full Bornholm comparison, outside model construction."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
import yaml


def require(condition, message):
    if not condition:
        raise ValueError(message)


def electricity_buses(n):
    return n.buses.index[n.buses.carrier.isin(["AC", "DC"])].difference(["BEI_HUB"])


def check_networks(
    networks, config, *, solved=False, reference=None, require_reference=False
):
    require(
        reference is not None or not require_reference,
        "Supply --reference with the other branch's built historical network",
    )
    baseline = networks["baseline"]
    expected = pd.date_range(**config["snapshots"], freq="h")
    buses = electricity_buses(baseline)
    require("BEI_HUB" not in baseline.buses.index, "Baseline contains the project hub")
    require(
        "BEI_WIND" not in baseline.generators.index, "Baseline contains project wind"
    )
    require(
        set(baseline.buses.loc[buses, "country"]) == set(config["countries"]),
        "Electricity buses do not cover the configured study countries",
    )
    for case, n in networks.items():
        require(n.snapshots.equals(expected), f"{case}: expected the full hourly year")
        require(np.allclose(n.snapshot_weightings, 1.0), f"{case}: non-hourly weights")
        require(
            electricity_buses(n).equals(buses), f"{case}: electricity bus set differs"
        )
        pd.testing.assert_series_equal(
            n.buses.loc[buses, "country"], baseline.buses.loc[buses, "country"]
        )
        require(n.buses.carrier.str.contains("heat").any(), f"{case}: no heat buses")
        require(n.buses.carrier.eq("H2").any(), f"{case}: no hydrogen buses")
        require(n.meta.get("sector", {}).get("enabled"), f"{case}: sector mode missing")
        require(
            n.meta.get("energy_islands", {}).get("exclude_project_area"),
            f"{case}: project-area exclusion is not recorded",
        )
        for group in config["clustering"]["copperplate_regions"]:
            require(set(group) <= set(buses), f"{case}: missing market buses {group}")
            copper = n.links.loc[
                n.links.carrier.eq("copper")
                & n.links.bus0.isin(group)
                & n.links.bus1.isin(group)
            ]
            require(
                len(copper) == 1, f"{case}: missing shared-price connection {group}"
            )
            require(
                np.isposinf(copper.p_nom).all()
                and not copper.p_nom_extendable.any()
                and np.allclose(copper.p_min_pu, -1)
                and np.allclose(copper.efficiency, 1)
                and np.allclose(copper.marginal_cost, 0),
                f"{case}: invalid shared-price connection {group}",
            )
            if solved:
                prices = n.buses_t.marginal_price[group]
                require(
                    np.isfinite(prices).all().all()
                    and np.allclose(
                        prices.iloc[:, 0], prices.iloc[:, 1], atol=0.01, rtol=1e-6
                    ),
                    f"{case}: shared-market prices differ for {group}",
                )
        if solved:
            require(np.isfinite(n.objective), f"{case}: no finite solved objective")
        print(
            f"{case}: {len(buses)} onshore electricity buses; {len(n.buses)} total buses; {len(n.snapshots)} hours"
        )
        if case == "baseline":
            continue
        require("BEI_WIND" in n.generators.index, f"{case}: project wind missing")
        wind = n.generators.loc["BEI_WIND"]
        require(
            wind.p_nom_extendable and wind.p_nom_min == wind.p_nom_max == 3000,
            f"{case}: Bornholm wind is not fixed at 3000 MW",
        )
        if solved:
            require(
                np.isclose(wind.p_nom_opt, 3000),
                f"{case}: incorrect solved wind capacity",
            )
        profile = n.generators_t.p_max_pu["BEI_WIND"]
        require(
            np.isfinite(profile).all()
            and profile.between(0, 1).all()
            and profile.max() > 0,
            f"{case}: invalid project wind profile",
        )
        require(
            set(n.buses.index[n.buses.location.eq("BEI_HUB")]) == {"BEI_HUB"},
            f"{case}: unintended sector buses were attached to the offshore hub",
        )
        expected_links = (
            {"BEI_TO_DK2_RADIAL": ("DKE1", 3000)}
            if case == "radial"
            else {"BEI_TO_DK2": ("DKE1", 1200), "BEI_TO_DE": ("DE00", 2000)}
        )
        for name, (target, capacity) in expected_links.items():
            forward, reverse = n.links.loc[name], n.links.loc[name + "-reversed"]
            require(
                forward.bus0 == reverse.bus1 == "BEI_HUB"
                and forward.bus1 == reverse.bus0 == target,
                f"{case}: incorrect endpoints for {name}",
            )
            require(
                forward.p_nom_min
                == forward.p_nom_max
                == reverse.p_nom_min
                == reverse.p_nom_max
                == capacity,
                f"{case}: incorrect export capacity for {name}",
            )
            require(
                0 < forward.efficiency < 1
                and np.isclose(forward.efficiency, reverse.efficiency)
                and forward.p_min_pu == reverse.p_min_pu == 0
                and reverse.capital_cost == reverse.length == 0,
                f"{case}: inconsistent losses or double-counted investment for {name}",
            )
            if solved:
                require(
                    np.allclose([forward.p_nom_opt, reverse.p_nom_opt], capacity),
                    f"{case}: incorrect solved export capacity for {name}",
                )
    pd.testing.assert_series_equal(
        networks["radial"].generators_t.p_max_pu["BEI_WIND"],
        networks["hybrid"].generators_t.p_max_pu["BEI_WIND"],
        check_exact=False,
        rtol=1e-12,
        atol=1e-12,
    )
    generic = baseline.generators.index[
        baseline.generators.carrier.str.startswith("offwind")
    ]
    for case in ["radial", "hybrid"]:
        pd.testing.assert_series_equal(
            baseline.generators.loc[generic, "p_nom_max"],
            networks[case].generators.loc[generic, "p_nom_max"],
            check_index_type=False,
        )
    if reference is not None:
        require(
            electricity_buses(reference).equals(buses),
            "Electricity bus IDs differ from the supplied historical reference network",
        )
        pd.testing.assert_series_equal(
            reference.buses.loc[buses, "country"], baseline.buses.loc[buses, "country"]
        )
        print("Electricity bus IDs match the historical reference network.")
    else:
        print("Actual historical network comparison not performed (no --reference).")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/config.bornholm.yaml")
    parser.add_argument("--solved", action="store_true")
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--require-reference", action="store_true")
    args = parser.parse_args()
    config = yaml.safe_load(Path(args.config).read_text())
    horizon = config["planning_horizons"][0]
    root, stage = ("results", "solved") if args.solved else ("resources", "composed")
    networks = {
        case: pypsa.Network(
            f"{root}/bornholm_zonal_{case}/networks/{stage}_{horizon}.nc"
        )
        for case in ["baseline", "radial", "hybrid"]
    }
    reference = pypsa.Network(args.reference) if args.reference else None
    check_networks(
        networks,
        config,
        solved=args.solved,
        reference=reference,
        require_reference=args.require_reference,
    )
    print("Bornholm comparison checks passed.")
