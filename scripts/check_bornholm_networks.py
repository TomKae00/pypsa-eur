# SPDX-License-Identifier: MIT
"""Standalone checks of the full Bornholm comparison, outside model construction."""

import argparse
import warnings
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


def scenarios(config):
    definitions = yaml.safe_load(Path(config["run"]["scenarios"]["file"]).read_text())
    return {
        name.removeprefix("bornholm_zonal_"): override
        for name, override in definitions.items()
    }


def project_links(n):
    return n.links.index[n.links.bus0.eq("BEI_HUB") | n.links.bus1.eq("BEI_HUB")]


def check_networks(networks, config, *, solved=False):
    baseline = networks["baseline"]
    definitions = scenarios(config)
    cables = pd.read_csv(config["energy_islands"]["links"])
    project_profiles = []
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
        spec = definitions[case]
        enabled = spec["energy_islands"]["enabled"]
        sector = n.meta.get("sector", {})
        for key, value in {
            "gas_network": True,
            "gas_network_connectivity_upgrade": 0,
            "H2_network": False,
            "H2_retrofit": False,
        }.items():
            require(
                sector.get(key) == value,
                f"{case}: stale or inconsistent sector setting {key}; rebuild networks",
            )
        gas_pipes = n.links.index[n.links.carrier.eq("gas pipeline")]
        require(
            not n.links.loc[gas_pipes, "p_nom_extendable"].any(),
            f"{case}: existing gas pipelines are expandable",
        )
        require(
            not n.links.carrier.eq("gas pipeline new").any(),
            f"{case}: new gas pipeline candidates are present",
        )
        require(
            not n.links.carrier.str.startswith("H2 pipeline").any(),
            f"{case}: hydrogen pipelines are present",
        )
        # The pipeline assumption is common to the electric-grid benchmarks too.
        baseline_gas = baseline.links.index[baseline.links.carrier.eq("gas pipeline")]
        pd.testing.assert_frame_equal(
            n.links.loc[gas_pipes, ["bus0", "bus1", "p_nom", "length"]].sort_index(),
            baseline.links.loc[
                baseline_gas, ["bus0", "bus1", "p_nom", "length"]
            ].sort_index(),
            check_dtype=False,
            check_index_type=False,
        )
        grid_policy = spec.get("electricity", {}).get(
            "transmission_limit", config["electricity"]["transmission_limit"]
        )
        background_dc = n.links.index[n.links.carrier.eq("DC")].difference(
            project_links(n)
        )
        if grid_policy == "vfixed":
            require(
                not n.lines.s_nom_extendable.any()
                and not n.links.loc[background_dc, "p_nom_extendable"].any(),
                f"{case}: background electricity grid is expandable",
            )
        # Capacity is endogenous in the expandable benchmark after solving.
        if not solved:
            pd.testing.assert_frame_equal(
                n.lines[["bus0", "bus1", "s_nom", "length"]].sort_index(),
                baseline.lines[["bus0", "bus1", "s_nom", "length"]].sort_index(),
                check_dtype=False,
                check_index_type=False,
            )
            baseline_dc = baseline.links.index[baseline.links.carrier.eq("DC")]
            pd.testing.assert_frame_equal(
                n.links.loc[
                    background_dc, ["bus0", "bus1", "p_nom", "length"]
                ].sort_index(),
                baseline.links.loc[
                    baseline_dc, ["bus0", "bus1", "p_nom", "length"]
                ].sort_index(),
                check_dtype=False,
                check_index_type=False,
            )
        # Fixed service/input demands; electricity used by conversion links can respond.
        pd.testing.assert_frame_equal(
            n.loads[["bus", "carrier", "p_set"]],
            baseline.loads[["bus", "carrier", "p_set"]],
            check_dtype=False,
            check_index_type=False,
        )
        pd.testing.assert_frame_equal(n.loads_t.p_set, baseline.loads_t.p_set)
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
            require(
                n.meta.get("solve_result")
                == {"status": "ok", "termination_condition": "optimal"},
                f"{case}: missing successful optimal solver termination",
            )
        print(
            f"{case}: {len(buses)} onshore electricity buses; {len(n.buses)} total buses; {len(n.snapshots)} hours"
        )
        if not enabled:
            require(
                "BEI_HUB" not in n.buses.index
                and "BEI_WIND" not in n.generators.index
                and project_links(n).empty,
                f"{case}: baseline contains project assets",
            )
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
        project_profiles.append(profile)
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
        expected_links = cables.loc[
            cables.case.eq(spec["energy_islands"]["case"])
        ].set_index("link_id")
        require(not expected_links.empty, f"{case}: no candidate cables configured")
        expected_names = set(expected_links.index) | {
            name + "-reversed" for name in expected_links.index
        }
        require(
            set(project_links(n)) == expected_names,
            f"{case}: unexpected/missing project cables",
        )
        for name, row in expected_links.iterrows():
            target, maximum = row.target_bus, row.capacity_mw
            minimum = maximum if row.build_mode == "fixed" else 0
            forward, reverse = n.links.loc[name], n.links.loc[name + "-reversed"]
            require(
                forward.bus0 == reverse.bus1 == "BEI_HUB"
                and forward.bus1 == reverse.bus0 == target,
                f"{case}: incorrect endpoints for {name}",
            )
            require(
                forward.p_nom_extendable
                and reverse.p_nom_extendable
                and forward.p_nom_min == reverse.p_nom_min == minimum
                and forward.p_nom_max == reverse.p_nom_max == maximum,
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
                    np.isclose(forward.p_nom_opt, reverse.p_nom_opt)
                    and minimum - 1e-3 <= forward.p_nom_opt <= maximum + 1e-3,
                    f"{case}: incorrect solved export capacity for {name}",
                )
                if minimum == 0 and np.isclose(forward.p_nom_opt, maximum, atol=1):
                    warnings.warn(
                        f"{case}: {name} reaches its {maximum} MW search bound; enlarge the bound before interpreting optimal sizing",
                        stacklevel=2,
                    )
    for profile in project_profiles[1:]:
        pd.testing.assert_series_equal(
            project_profiles[0], profile, check_exact=False, rtol=1e-12, atol=1e-12
        )
    generic = baseline.generators.index[
        baseline.generators.carrier.str.startswith("offwind")
    ]
    for case in networks:
        pd.testing.assert_series_equal(
            baseline.generators.loc[generic, "p_nom_max"],
            networks[case].generators.loc[generic, "p_nom_max"],
            check_index_type=False,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/config.bornholm.yaml")
    parser.add_argument("--solved", action="store_true")
    parser.add_argument(
        "--scenarios",
        nargs="+",
        help="Short scenario names; baseline is always included",
    )
    args = parser.parse_args()
    config = yaml.safe_load(Path(args.config).read_text())
    horizon = config["planning_horizons"][0]
    root, stage = ("results", "solved") if args.solved else ("resources", "composed")
    selected = args.scenarios or list(scenarios(config))
    unknown = set(selected).difference(scenarios(config))
    if unknown:
        parser.error(f"Unknown scenarios: {sorted(unknown)}")
    selected = list(dict.fromkeys(["baseline", *selected]))
    networks = {
        case: pypsa.Network(
            f"{root}/bornholm_zonal_{case}/networks/{stage}_{horizon}.nc"
        )
        for case in selected
    }
    check_networks(networks, config, solved=args.solved)
    print("Bornholm comparison checks passed.")
