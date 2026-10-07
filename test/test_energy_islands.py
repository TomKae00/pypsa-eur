"""Focused regression checks for the Bornholm physical-input patch."""

import numpy as np
import pandas as pd
import pypsa
import pytest

from scripts.add_energy_islands import (
    _add_hubs,
    _add_links,
    _haversine_km,
    _select_profile_generator,
    _select_target_bus,
    _subtract_generic_potential,
    main,
)
from scripts.inspect_energy_islands import inventory


@pytest.fixture
def network():
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2013-01-01", periods=4, freq="h"))
    n.snapshot_weightings.loc[:, :] = 2190.0
    for carrier in ["AC", "DC", "offwind-dc", "gas"]:
        n.add("Carrier", carrier)
    for name, country, x, y in [
        ("DK2_TEST", "DK", 12.0, 55.5),
        ("DE_TEST", "DE", 13.0, 54.0),
        ("SE_TEST", "SE", 14.0, 56.0),
    ]:
        n.add("Bus", name, country=country, x=x, y=y, carrier="AC", v_nom=380)
    for name, bus in [("DK offshore", "DK2_TEST"), ("SE offshore", "SE_TEST")]:
        n.add(
            "Generator",
            name,
            bus=bus,
            carrier="offwind-dc",
            p_nom_extendable=True,
            p_nom_max=5000,
            capital_cost=1e9,
            p_max_pu=[1.0, 0.0, 0.5, 1.0],
        )
    return n


@pytest.fixture
def tables():
    hubs = pd.read_csv("data/energy_islands/bornholm_hubs.csv")
    wind = pd.read_csv("data/energy_islands/bornholm_wind.csv")
    wind["source_generator"] = "DK offshore"
    links = pd.read_csv("data/energy_islands/bornholm_links.csv")
    links["target_bus"] = links.target_country.map({"DK": "DK2_TEST", "DE": "DE_TEST"})
    costs = pd.DataFrame(
        {
            "capital_cost": [100000.0, 10.0, 10.0, 100.0],
            "marginal_cost": 0.0,
            "efficiency": 1.0,
            "lifetime": 30.0,
        },
        index=["offwind", "HVDC overhead", "HVDC submarine", "HVDC inverter pair"],
    )
    return hubs, wind, links, costs


def subtract(n, allocation, note=""):
    _subtract_generic_potential(n, "offwind-dc", 3000, set(), allocation, "DK", note)


def test_missing_target_is_not_guessed(network, tables):
    row = tables[2].iloc[0].copy()
    row["target_bus"] = np.nan
    with pytest.raises(ValueError, match="set target_bus explicitly"):
        _select_target_bus(network, row)


def test_exact_target_country_is_checked(network, tables):
    row = tables[2].iloc[0].copy()
    row["target_bus"] = "SE_TEST"
    with pytest.raises(ValueError, match="not in DK"):
        _select_target_bus(network, row)


def test_hub_cannot_be_mainland_target(network, tables):
    _add_hubs(network, tables[0])
    row = tables[2].iloc[0].copy()
    row["target_bus"] = "BEI_HUB"
    with pytest.raises(ValueError, match="must not be an energy-island hub"):
        _select_target_bus(network, row)


def test_route_cost_is_independent_of_cluster_coordinates(network, tables):
    hubs, _, links, costs = tables
    _add_hubs(network, hubs)
    other = network.copy()
    other.buses.loc["DK2_TEST", ["x", "y"]] = [15.01, 55.12]
    _add_links(network, links.iloc[[0]], costs, 1.25)
    _add_links(other, links.iloc[[0]], costs, 1.25)
    expected = float(_haversine_km(15.0, 55.12, 12.27, 55.6)) * 1.25
    assert network.links.at["BEI_TO_DK2", "length"] == pytest.approx(expected)
    assert other.links.at["BEI_TO_DK2", "capital_cost"] == pytest.approx(
        network.links.at["BEI_TO_DK2", "capital_cost"]
    )


def test_explicit_route_is_not_multiplied_twice(network, tables):
    hubs, _, links, costs = tables
    _add_hubs(network, hubs)
    links = links.iloc[[0]].copy()
    links["length_km"] = 250.0
    _add_links(network, links, costs, 1.25)
    assert network.links.at["BEI_TO_DK2", "length"] == 250.0
    assert network.links.at["BEI_TO_DK2", "capital_cost"] == 250.0 * 10.0 + 100.0


def test_profile_must_be_explicit(network):
    with pytest.raises(ValueError, match="Set source_generator explicitly"):
        _select_profile_generator(network, "offwind-dc", "")


def test_profile_must_cover_snapshots(network):
    network.generators_t.p_max_pu.loc[network.snapshots[1], "DK offshore"] = np.nan
    with pytest.raises(ValueError, match="cover every snapshot"):
        _select_profile_generator(network, "offwind-dc", "DK offshore")


def test_no_automatic_potential_subtraction(network):
    before = network.generators.p_nom_max.copy()
    with pytest.raises(ValueError, match="No automatic subtraction"):
        subtract(network, None)
    pd.testing.assert_series_equal(before, network.generators.p_nom_max)


def test_potential_cannot_spill_to_sweden(network):
    before = network.generators.p_nom_max.copy()
    with pytest.raises(ValueError, match="outside project country"):
        subtract(network, {"DK offshore": 1000, "SE offshore": 2000})
    pd.testing.assert_series_equal(before, network.generators.p_nom_max)


def test_existing_capacity_is_preserved(network):
    network.generators.at["DK offshore", "p_nom"] = 2500
    with pytest.raises(ValueError, match="existing/minimum capacity"):
        subtract(network, {"DK offshore": 3000})
    assert network.generators.at["DK offshore", "p_nom_max"] == 5000


def test_partial_overlap_requires_documentation(network):
    with pytest.raises(ValueError, match="document why"):
        subtract(network, {"DK offshore": 973})
    subtract(
        network,
        {"DK offshore": 973},
        "Synthetic test: only 973 MW overlap this resource pool",
    )
    assert network.generators.at["DK offshore", "p_nom_max"] == 4027
    assert network.generators.at["SE offshore", "p_nom_max"] == 5000


def test_inventory_requires_baseline(network, tables):
    buses, generators = inventory(network, tables[2], tables[1])
    assert "DK2_TEST" in buses.bus_id.values
    assert generators.has_valid_dynamic_profile.all()
    _add_hubs(network, tables[0])
    with pytest.raises(ValueError, match="BASELINE"):
        inventory(network, tables[2], tables[1])


def test_composition_solve_and_roundtrip(network, tables, tmp_path):
    hubs, wind, links, costs = tables
    paths = {}
    for key, frame in [("hubs", hubs), ("wind", wind), ("links", links)]:
        paths[f"energy_island_{key}"] = tmp_path / f"{key}.csv"
        frame.to_csv(paths[f"energy_island_{key}"], index=False)
    config = {
        "enabled": True,
        "case": "hybrid",
        "potential_allocation": {"BEI_WIND": {"DK offshore": 3000}},
    }
    main(network, paths, config, costs, current_horizon=2035)
    assert network.generators.at["DK offshore", "p_nom_max"] == 2000
    assert network.generators.at["SE offshore", "p_nom_max"] == 5000
    assert network.generators.at["BEI_WIND", "p_nom_min"] == 3000
    assert network.generators.at["BEI_WIND", "p_nom_max"] == 3000
    for bus, load, cost in [("DK2_TEST", 500, 120), ("DE_TEST", 2500, 60)]:
        network.add("Load", bus, bus=bus, p_set=load)
        network.add(
            "Generator",
            f"gas {bus}",
            bus=bus,
            carrier="gas",
            p_nom=10000,
            marginal_cost=cost,
        )
    status, condition = network.optimize(solver_name="highs", threads=1)
    assert status == "ok" and condition == "optimal"
    capital = (network.generators.capital_cost * network.generators.p_nom_opt).sum() + (
        network.links.capital_cost * network.links.p_nom_opt
    ).sum()
    operating = float(
        network.snapshot_weightings.objective
        @ network.generators_t.p.mul(network.generators.marginal_cost).sum(axis=1)
    )
    assert network.objective == pytest.approx(capital + operating)
    network.export_to_netcdf(tmp_path / "solved.nc")
    restored = pypsa.Network(tmp_path / "solved.nc")
    assert restored.objective == pytest.approx(network.objective)
    pd.testing.assert_series_equal(restored.links.bus1, network.links.bus1)
    assert restored.generators.at["DK offshore", "p_nom_max"] == 2000
    assert restored.generators.at["SE offshore", "p_nom_max"] == 5000
