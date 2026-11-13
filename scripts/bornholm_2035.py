# SPDX-License-Identifier: MIT
"""The agreed 2035 scenario matrix and checks, using native sector assumptions."""

from pathlib import Path

import numpy as np
import yaml

from scripts import bornholm_welfare_sensitivity as welfare

AUDIT_KEY = "bornholm_2035_audit"
SPEC = Path("config/bornholm_2035_campaign.yaml")
SCENARIOS = Path("config/scenarios_bornholm_2035.yaml")


def read_config():
    config = yaml.safe_load(SPEC.read_text())
    welfare.require(
        config["version"] == 2 and config["planning_year"] == 2035,
        "Expected the 2035 campaign specification",
    )
    capacities = yaml.safe_load(Path(config["capacity_assumptions"]).read_text())
    config.update(deployment=capacities["deployment"], solar=capacities["solar"])
    return config


def scenarios(config):
    rows = []

    def add(name, world, project, v2g):
        rows.append(
            dict(
                index=len(rows),
                name=name,
                world=world,
                project=project,
                deployment="central",
                v2g_share=float(v2g),
            )
        )

    for world in ("central", "ev55", "ev65"):
        for project in ("baseline", "hybrid", "split_radial"):
            add(f"{world}_{project}", world, project, config["central_v2g_share"])
    for project in ("radial_dk", "radial_de"):
        add(f"central_{project}", "central", project, config["central_v2g_share"])
    for label, share in (("v2g0", 0.0), ("v2g50", 0.5)):
        for project in ("baseline", "hybrid"):
            add(f"{label}_{project}", "central", project, share)
    return rows


def definitions(config):
    result = {}
    for case in scenarios(config):
        world = config["worlds"][case["world"]]
        settings = dict(
            enabled=True,
            scenario=case,
            additional_data_centres_TWh=world["additional_data_centres_TWh"],
        )
        result[case["name"]] = {
            "energy_islands": {
                "enabled": case["project"] != "baseline",
                "case": "radial" if case["project"] == "radial_dk" else case["project"],
            },
            "sector": {
                "land_transport_electric_share": {2035: world["ev_share"]},
                "land_transport_ice_share": {2035: 1.0 - world["ev_share"]},
                "land_transport_fuel_cell_share": {2035: 0.0},
                "bev_v2g_availability": case["v2g_share"],
            },
            "industry": {
                "sector_ratios_fraction_future": {2035: world["industrial_progress"]}
            },
            "bornholm_2035": settings,
            "bornholm_welfare_sensitivity": {
                "enabled": True,
                "scenario": case,
                "capacity_limits": welfare.capacity_limits(config, case),
            },
        }
    return result


def demand_report(n):
    """Final service demands by carrier; converter electricity is reported after solve."""
    ids = n.loads.index
    panel = welfare.demand_panel(n, ids, nonnegative=False)
    annual = panel.mul(n.snapshot_weightings.generators, axis=0).sum() / 1e6
    country = welfare.countries(n, n.loads.bus)
    rows = []
    for ct in ("DK", "DE"):
        selected = ids[country.eq(ct)]
        for carrier, names in n.loads.loc[selected].groupby("carrier").groups.items():
            rows.append(
                dict(
                    country=ct,
                    carrier=carrier,
                    final_demand_TWh=float(annual.loc[names].sum()),
                )
            )
    return rows


def electricity_use(n):
    """
    End-use/fuel-conversion electricity, excluding grid transfers and battery cycling.

    EV driving is measured at the battery and divided by charger efficiency.
    Thus V2G arbitrage does not inflate mobility demand. Other fuel-conversion
    withdrawals, including electrolysis, are explicit and remain identifiable.
    """
    rows = []
    weights = n.snapshot_weightings.generators
    electric = n.buses.carrier.isin(["AC", "low voltage"])
    countries = welfare.countries(n, n.loads.bus)
    demand = welfare.demand_panel(n, n.loads.index, nonnegative=False)
    for ct in ("DK", "DE"):
        for carrier, ids in (
            n.loads.loc[countries.eq(ct)].groupby("carrier").groups.items()
        ):
            if carrier == "land transport EV":
                amount = demand[ids].mul(weights, axis=0).sum().sum()
                amount /= n.meta["sector"]["bev_charge_efficiency"]
                label = "EV driving including charging losses"
            else:
                ids = ids[n.loads.loc[ids, "bus"].map(electric)]
                if len(ids) == 0:
                    continue
                amount = demand[ids].mul(weights, axis=0).sum().sum()
                label = carrier
            rows.append(
                dict(country=ct, carrier=label, electricity_TWh=float(amount / 1e6))
            )
    for column in n.links.columns[n.links.columns.str.fullmatch(r"bus\d+")]:
        port = column[3:]
        flows = n.links_t.get(f"p{port}")
        if flows is None or flows.empty:
            continue
        bus = n.links[column]
        country = welfare.countries(n, bus)
        excluded = n.links.carrier.isin(
            [
                "DC",
                "B2B",
                "copper",
                "electricity distribution grid",
                "BEV charger",
                "V2G",
            ]
        )
        excluded |= n.links.carrier.str.contains("battery", case=False)
        for ct in ("DK", "DE"):
            ids = n.links.index[
                country.eq(ct) & bus.map(electric).fillna(False) & ~excluded
            ]
            ids = ids.intersection(flows.columns)
            for carrier, names in n.links.loc[ids].groupby("carrier").groups.items():
                amount = flows[names].clip(lower=0).mul(weights, axis=0).sum().sum()
                if amount > 1e-6:
                    rows.append(
                        dict(
                            country=ct,
                            carrier=carrier,
                            electricity_TWh=float(amount / 1e6),
                        )
                    )
    return rows


def apply(n, settings):
    if not settings.get("enabled", False):
        return
    welfare.require(AUDIT_KEY not in n.meta, "2035 assumptions already applied")
    hours = float(n.snapshot_weightings.generators.sum())
    welfare.require(np.isclose(hours, 8760), "Expected one complete 2035 model year")
    buses = {"DK1": "DKW1", "DK2": "DKE1", "DE": "DE00"}
    for region, twh in settings["additional_data_centres_TWh"].items():
        welfare.require(
            np.isfinite(twh) and twh >= 0, "Invalid additional data-centre demand"
        )
        if twh == 0:
            continue
        bus = buses[region]
        welfare.require(
            bus in n.buses.index and n.buses.at[bus, "carrier"] == "AC",
            f"Missing receiving electricity bus {bus}",
        )
        n.add(
            "Load",
            f"{bus} additional data centres",
            bus=bus,
            carrier="additional data centres",
            p_set=float(twh) * 1e6 / hours,
        )
    panel = welfare.demand_panel(n, n.loads.index, nonnegative=False)
    n.meta[AUDIT_KEY] = dict(
        settings=settings,
        demand_sha256=welfare.panel_digest(panel),
        final_demand=demand_report(n),
    )


def project_check(n, project):
    names = n.generators.index[n.generators.index.str.startswith("BEI_WIND")]
    expected = (
        {"BEI_WIND_DK": 1000.0, "BEI_WIND_DE": 2000.0}
        if project == "split_radial"
        else {}
        if project == "baseline"
        else {"BEI_WIND": 3000.0}
    )
    welfare.require(set(names) == set(expected), "Wrong project wind assets")
    for name, capacity in expected.items():
        g = n.generators.loc[name]
        welfare.require(
            g.p_nom_extendable and g.p_nom_min == g.p_nom_max == capacity,
            f"Wrong mandatory capacity for {name}",
        )
    expected_links = {
        "baseline": {},
        "hybrid": {"BEI_TO_DK2": 1200, "BEI_TO_DE": 2000},
        "split_radial": {"BEI_TO_DK2_SPLIT": 1000, "BEI_TO_DE_SPLIT": 2000},
        "radial_dk": {"BEI_TO_DK2_RADIAL": 3000},
        "radial_de": {"BEI_TO_DE_RADIAL": 3000},
    }[project]
    for name, capacity in expected_links.items():
        welfare.require(name in n.links.index, f"Missing project export {name}")
        link = n.links.loc[name]
        welfare.require(
            link.p_nom_extendable and link.p_nom_min == link.p_nom_max == capacity,
            f"Wrong export capacity for {name}",
        )
    if project == "split_radial":
        welfare.require(
            n.generators.at["BEI_WIND_DK", "bus"]
            != n.generators.at["BEI_WIND_DE", "bus"],
            "Split radial cannot share an offshore bus",
        )
        offshore = {"BEI_HUB_DK", "BEI_HUB_DE"}
        links = n.links[n.links.bus0.isin(offshore) | n.links.bus1.isin(offshore)]
        welfare.require(
            not (links.bus0.isin(offshore) & links.bus1.isin(offshore)).any(),
            "Split radial contains an offshore interconnector",
        )
        a, b = n.generators_t.p_max_pu[["BEI_WIND_DK", "BEI_WIND_DE"]].T.to_numpy()
        welfare.require(
            np.array_equal(a, b), "Split radial must use the same project resource"
        )
    return dict(wind_MW=sum(expected.values()), export_MW=sum(expected_links.values()))


def validate(n, config, case, solved=False):
    settings = definitions(config)[case["name"]]["bornholm_2035"]
    audit = n.meta.get(AUDIT_KEY, {})
    welfare.require(audit.get("settings") == settings, "Wrong 2035 scenario settings")
    welfare.require(
        audit.get("demand_sha256")
        == welfare.panel_digest(
            welfare.demand_panel(n, n.loads.index, nonnegative=False)
        ),
        "Final demands changed",
    )
    sector = n.meta["sector"]
    world = config["worlds"][case["world"]]
    shares = sector["land_transport_electric_share"]
    welfare.require(
        np.isclose(shares.get(2035, shares.get("2035")), world["ev_share"]),
        "Wrong native EV share",
    )
    welfare.require(
        np.isclose(sector["bev_dsm_availability"], config["managed_charging_share"]),
        "Managed charging must remain fixed",
    )
    welfare.require(
        np.isclose(sector["bev_v2g_availability"], case["v2g_share"]),
        "Wrong independent V2G share",
    )
    for tech in ("offwind-ac", "offwind-dc", "offwind-float"):
        welfare.require(
            n.meta["renewable"][tech]["correction_factor"] == 0.85,
            "Wrong background offshore correction",
        )
    welfare.require(
        n.meta["energy_islands"]["project_correction_factor"] == 0.8855,
        "Wrong dedicated Bornholm correction",
    )
    project_check(n, case["project"])
    return welfare.capacity_report(
        n, welfare.capacity_limits(config, case), solved=solved
    )
