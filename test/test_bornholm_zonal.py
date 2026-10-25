# SPDX-License-Identifier: MIT
"""Bidding-zone, sector-loss and full-run verification regressions."""

from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
import pytest
import yaml

from scripts._helpers import sanitize_custom_columns
from scripts.add_electricity import set_transmission_costs
from scripts.add_energy_islands import main as add_energy_islands
from scripts.check_bornholm_networks import check_networks, scenarios
from scripts.cluster_network import copperplate_buses
from scripts.prepare_sector_network import (
    lossy_bidirectional_links,
    remove_elec_base_techs,
)
from scripts.solve_network import add_lossy_bidirectional_link_constraints
from test.test_energy_island_resources import project_profile, setup_project  # noqa: F401


@pytest.mark.parametrize("case", ["radial", "hybrid"])
def test_sector_project_losses_and_costs(project_profile, case, tmp_path):  # noqa: F811
    n, inputs, config, costs = setup_project(project_profile, case)
    n.add("Carrier", ["heat", "H2", "heat pump", "electrolysis"])
    n.add("Bus", "DK heat", carrier="heat", location="DK")
    n.add("Bus", "DK H2", carrier="H2", location="DK")
    n.add(
        "Link",
        "heat pump",
        bus0="DK",
        bus1="DK heat",
        carrier="heat pump",
        efficiency=3,
        p_nom=500,
    )
    n.add(
        "Link",
        "electrolysis",
        bus0="DK",
        bus1="DK H2",
        carrier="electrolysis",
        efficiency=0.7,
        p_nom=1000,
    )
    n.add("Load", "heat demand", bus="DK heat", p_set=400)
    n.add("Load", "hydrogen demand", bus="DK H2", p_set=300)
    n.add("Load", "electricity demand", bus="DK", p_set=1000)
    n.add("Generator", "DK gas", bus="DK", carrier="gas", p_nom=4000, marginal_cost=100)
    n.add("Generator", "DE gas", bus="DE", carrier="gas", p_nom=4000, marginal_cost=10)
    n.add(
        "Link",
        "existing",
        bus0="DK",
        bus1="DE",
        carrier="DC",
        p_nom=0,
        p_min_pu=-1,
        length=100,
        underwater_fraction=1,
    )
    efficiency = {"efficiency_static": 0.98, "efficiency_per_1000km": 0.977}
    lossy_bidirectional_links(n, "DC", efficiency)
    costs.loc["HVAC overhead"] = [10, 0, 1, 30]
    set_transmission_costs(n, costs)
    assert n.links.at["existing", "capital_cost"] == 2100
    assert n.links.at["existing-reversed", "capital_cost"] == 0
    existing = n.links.loc[["existing", "existing-reversed"]].copy()

    # Sector additions precede project generation; project cables follow grid prep.
    add_energy_islands(n, inputs, config, costs, 2035, stage="generation")
    n.generators_t.p_max_pu.loc[n.snapshots[0], "BEI_WIND"] = 0
    before = n.links.index.copy()
    add_energy_islands(n, inputs, config, costs, 2035, stage="links")
    exports = n.links.index.difference(before)
    investment = (
        n.links.loc[exports, "p_nom_min"] @ n.links.loc[exports, "capital_cost"]
    )
    lossy_bidirectional_links(n, "DC", efficiency, link_indices=exports)
    pd.testing.assert_frame_equal(
        n.links.loc[existing.index, existing.columns], existing, check_column_type=False
    )
    assert not n.links.index.str.contains("reversed-reversed").any()
    assert n.links.reversed.notna().all()
    for cable in exports:
        reverse = cable + "-reversed"
        assert n.links.at[reverse, "capital_cost"] == 0
        assert n.links.at[reverse, "length"] == 0
        assert n.links.at[cable, "efficiency"] == n.links.at[reverse, "efficiency"]
    assert n.links.loc[n.links.carrier.eq("DC"), "p_nom_min"] @ n.links.loc[
        n.links.carrier.eq("DC"), "capital_cost"
    ] == pytest.approx(investment)
    sanitize_custom_columns(n)
    path = tmp_path / "sector.nc"
    n.export_to_netcdf(path)
    n = pypsa.Network(path)
    status, condition = n.optimize(
        solver_name="highs",
        threads=1,
        extra_functionality=lambda n, snapshots: (
            add_lossy_bidirectional_link_constraints(n)
        ),
    )
    assert (status, condition) == ("ok", "optimal")
    assert "Link-bidirectional_sync" in n.model.constraints
    for cable in exports:
        reverse = cable + "-reversed"
        assert n.links.at[cable, "p_nom_opt"] == n.links.at[reverse, "p_nom_opt"]
        for direction in [cable, reverse]:
            np.testing.assert_allclose(
                -n.links_t.p1[direction],
                n.links_t.p0[direction] * n.links.at[direction, "efficiency"],
                atol=1e-7,
            )
    if case == "hybrid":
        assert n.links_t.p0.at[n.snapshots[0], "BEI_TO_DE-reversed"] > 0
    assert n.links_t.p0[exports].sum().sum() > 0
    np.testing.assert_allclose(-n.links_t.p1["heat pump"], 400)
    np.testing.assert_allclose(-n.links_t.p1["electrolysis"], 300)
    capex = (
        n.generators.p_nom_opt @ n.generators.capital_cost
        + n.links.p_nom_opt @ n.links.capital_cost
    )
    opex = n.snapshot_weightings.objective @ n.generators_t.p.mul(
        n.generators.marginal_cost
    ).sum(axis=1)
    assert n.objective == pytest.approx(capex + opex)


def test_copperplate_survives_sector_cleanup_and_shares_prices(tmp_path):
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2023-01-01", periods=2, freq="h"))
    n.add("Carrier", ["AC", "DC", "gas"])
    n.add("Bus", ["DE00", "LUG1"], carrier="AC", country=["DE", "LU"])
    n.add("Link", "old", bus0="DE00", bus1="LUG1", carrier="DC", p_nom=1)
    copperplate_buses(n, [["DE00", "LUG1"]])
    config = yaml.safe_load(Path("config/config.bornholm.yaml").read_text())
    remove_elec_base_techs(n, {"Link": config["pypsa_eur"]["Link"]})
    assert list(n.buses.country) == ["DE", "LU"]
    assert "old" not in n.links.index
    assert n.links.carrier.eq("copper").all() and len(n.links) == 1
    n.add(
        "Generator", "DE gas", bus="DE00", carrier="gas", p_nom=1000, marginal_cost=10
    )
    n.add(
        "Generator", "LU gas", bus="LUG1", carrier="gas", p_nom=1000, marginal_cost=100
    )
    n.add("Load", "LU load", bus="LUG1", p_set=100)
    path = tmp_path / "copper.nc"
    n.export_to_netcdf(path)
    n = pypsa.Network(path)
    assert n.optimize(solver_name="highs", threads=1) == ("ok", "optimal")
    np.testing.assert_allclose(n.buses_t.marginal_price[["DE00", "LUG1"]], 10)
    np.testing.assert_allclose(n.generators_t.p["DE gas"], 100)
    np.testing.assert_allclose(n.generators_t.p["LU gas"], 0)


@pytest.mark.parametrize(
    "groups",
    [
        [["DE00", "missing"]],
        [["DE00"]],
        [["DE00", "DE00"]],
        [["DE00", "LUG1"], ["LUG1", "DKE1"]],
        [["DE00", "heat"]],
    ],
)
def test_copperplate_rejects_invalid_groups(groups):
    n = pypsa.Network()
    n.add("Bus", ["DE00", "LUG1", "DKE1"], carrier="AC")
    n.add("Bus", "heat", carrier="heat")
    with pytest.raises(ValueError):
        copperplate_buses(n, groups)
    assert n.links.empty


@pytest.fixture
def comparison():
    """Full-hourly check fixture; no optimisation or external input downloads."""
    config = yaml.safe_load(Path("config/config.bornholm.yaml").read_text())
    config["countries"] = ["DK", "DE", "LU", "IE", "GB", "SE", "PL"]
    n = pypsa.Network()
    n.set_snapshots(pd.date_range(**config["snapshots"], freq="h"))
    n.add(
        "Bus",
        ["DKE1", "DE00", "LUG1", "IE00", "GBNI", "SE04", "PL00"],
        carrier="AC",
        country=config["countries"],
    )
    n.add("Bus", "DK heat", carrier="heat", location="DKE1")
    n.add("Bus", "DK H2", carrier="H2", location="DKE1")
    n.add("Load", "heat service", bus="DK heat", p_set=100)
    n.add(
        "Link",
        "existing grid",
        bus0="DKE1",
        bus1="DE00",
        carrier="DC",
        p_nom=100,
        p_min_pu=-1,
        length=100,
    )
    lossy_bidirectional_links(n, "DC", {"efficiency_static": 0.98})
    n.add(
        "Generator",
        "generic offshore",
        bus="DKE1",
        carrier="offwind-ac",
        p_nom_extendable=True,
        p_nom_max=10000,
    )
    copperplate_buses(n, config["clustering"]["copperplate_regions"])
    n.meta = {
        "sector": {"enabled": True},
        "energy_islands": {"exclude_project_area": True},
    }
    reference = n.copy()
    networks = {"baseline": n}
    cables = pd.read_csv(config["energy_islands"]["links"])
    for case, spec in scenarios(config).items():
        if case == "baseline":
            continue
        m = n.copy()
        if spec.get("electricity", {}).get("transmission_limit") == "vopt":
            dc = m.links.carrier.eq("DC")
            m.links.loc[dc, "p_nom_extendable"] = True
            m.links.loc[dc, "p_nom_min"] = m.links.loc[dc, "p_nom"]
        if not spec["energy_islands"]["enabled"]:
            networks[case] = m
            continue
        m.add("Bus", "BEI_HUB", carrier="AC", country="DK", location="BEI_HUB")
        m.add(
            "Generator",
            "BEI_WIND",
            bus="BEI_HUB",
            p_nom_extendable=True,
            p_nom_min=3000,
            p_nom_max=3000,
            p_max_pu=pd.Series(0.5, index=m.snapshots),
        )
        selected = cables.loc[cables.case.eq(spec["energy_islands"]["case"])]
        for row in selected.itertuples():
            m.add(
                "Link",
                row.link_id,
                bus0="BEI_HUB",
                bus1=row.target_bus,
                carrier="DC",
                p_min_pu=-1,
                p_nom_extendable=True,
                p_nom_min=row.capacity_mw if row.build_mode == "fixed" else 0,
                p_nom_max=row.capacity_mw,
                length=100,
                capital_cost=200,
            )
        lossy_bidirectional_links(
            m, "DC", {"efficiency_static": 0.98}, link_indices=selected.link_id
        )
        networks[case] = m
    return networks, config, reference


def test_comparison_accepts_matching_reference(comparison):
    networks, config, reference = comparison
    check_networks(networks, config, reference=reference, require_reference=True)


@pytest.mark.parametrize(
    "fault,message",
    [
        ("missing_reference", "Supply --reference"),
        ("different_reference", "bus IDs differ"),
        ("baseline_hub", "Baseline contains the project hub"),
        ("nonhourly", "full hourly year"),
        ("missing_copper", "missing shared-price connection"),
        ("double_cost", "double-counted investment"),
    ],
)
def test_comparison_rejects_invalid_runs(comparison, fault, message):
    networks, config, reference = comparison
    if fault == "missing_reference":
        reference = None
    elif fault == "different_reference":
        reference.add("Bus", "extra", carrier="AC")
    elif fault == "baseline_hub":
        networks["baseline"].add("Bus", "BEI_HUB", carrier="AC")
    elif fault == "nonhourly":
        networks["hybrid"].set_snapshots(networks["hybrid"].snapshots[::2])
    elif fault == "missing_copper":
        m = networks["hybrid"]
        m.remove("Link", m.links.index[m.links.carrier.eq("copper")][0])
    elif fault == "double_cost":
        networks["hybrid"].links.at["BEI_TO_DE-reversed", "capital_cost"] = 200
    with pytest.raises(ValueError, match=message):
        check_networks(networks, config, reference=reference, require_reference=True)


def test_comparison_rejects_changed_service_demand(comparison):
    networks, config, _ = comparison
    networks["open"].loads.at["heat service", "p_set"] = 200
    with pytest.raises(AssertionError):
        check_networks(networks, config)


def test_comparison_rejects_expandable_main_grid(comparison):
    networks, config, _ = comparison
    networks["optimised"].links.at["existing grid", "p_nom_extendable"] = True
    with pytest.raises(ValueError, match="background electricity grid is expandable"):
        check_networks(networks, config)
