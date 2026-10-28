# SPDX-License-Identifier: MIT
"""Window accounting and preflight tests for the short Bornholm diagnostic."""

import copy

import numpy as np
import pandas as pd
import pypsa
import pytest

from scripts.run_bornholm_diagnostic import annual_gas_supply_check, short_network


def test_short_window_preserves_physical_storage_and_scales_annual_inputs():
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2013-01-01", periods=8760, freq="h"))
    n.add("Carrier", ["gas", "non-sequestered HVC"])
    n.add("Bus", ["gas", "hvc"], carrier=["gas", "non-sequestered HVC"])
    n.add(
        "Store",
        "gas storage",
        bus="gas",
        carrier="gas",
        e_nom_extendable=True,
        e_nom_min=1000,
        e_cyclic=True,
        capital_cost=8760,
    )
    n.add(
        "Store",
        "annual HVC",
        bus="hvc",
        carrier="non-sequestered HVC",
        e_nom=8760,
        e_initial=8760,
        e_max_pu=pd.Series(1.0, index=n.snapshots),
    )
    n.stores_t.e_max_pu.iloc[-1, 0] = 0
    n.add(
        "Generator",
        "biogas quota",
        bus="gas",
        carrier="gas",
        p_nom=10,
        e_sum_min=8760,
        e_sum_max=17520,
        capital_cost=4380,
    )
    n.add("GlobalConstraint", "CO2", type="co2_atmosphere", constant=8760)
    n.add("GlobalConstraint", "biomass quota", type="operational_limit", constant=17520)
    n.add(
        "GlobalConstraint",
        "grid volume",
        type="transmission_volume_expansion_limit",
        constant=500,
    )
    before = n.stores.copy()
    meta = copy.deepcopy(n.meta)
    small = short_network(n, 24, 48)
    assert len(small.snapshots) == 24 and small.snapshots[0] == n.snapshots[48]
    assert np.allclose(small.snapshot_weightings, 1)
    assert small.stores.at["gas storage", "e_nom_min"] == 1000
    assert small.stores.at["gas storage", "e_cyclic"]
    assert small.stores.at["gas storage", "capital_cost"] == 24
    assert small.generators.at["biogas quota", "e_sum_min"] == 24
    assert small.generators.at["biogas quota", "e_sum_max"] == 48
    assert small.global_constraints.at["CO2", "constant"] == 24
    assert small.global_constraints.at["biomass quota", "constant"] == 48
    assert small.global_constraints.at["grid volume", "constant"] == 500
    assert small.stores.at["annual HVC", "e_initial"] == 24
    assert small.stores.at["annual HVC", "e_nom"] == 24
    assert small.stores_t.e_max_pu.at[small.snapshots[-1], "annual HVC"] == 0
    assert small.stores_t.e_max_pu.at[small.snapshots[0], "annual HVC"] == 1
    pd.testing.assert_frame_equal(n.stores, before)
    assert n.meta == meta


def test_annual_gas_preflight_detects_zero_capacity_connection():
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2013-01-01", periods=24, freq="h"))
    n.add(
        "Bus",
        ["HU gas", "RS gas", "RS industry"],
        carrier=["gas", "gas", "gas for industry"],
    )
    n.add("Generator", "existing gas supply", bus="HU gas", p_nom=100)
    n.add("Link", "existing pipe", bus0="HU gas", bus1="RS gas", p_nom=0)
    n.add(
        "Link",
        "industrial supply",
        bus0="RS gas",
        bus1="RS industry",
        p_nom_extendable=True,
    )
    n.add("Load", "fixed gas demand", bus="RS industry", p_set=50)
    with pytest.raises(ValueError, match="Annual fixed-gas supply shortfall"):
        annual_gas_supply_check(n)
    n.links.at["existing pipe", "p_nom"] = 100
    assert annual_gas_supply_check(n)["optimistic_supply_shortfall_TWh"] == 0
    n.links.at["existing pipe", "p_nom"] = 0
    n.links.at["existing pipe", "p_nom_extendable"] = True
    assert annual_gas_supply_check(n)["optimistic_supply_shortfall_TWh"] == 0


def test_window_can_solve_with_native_extra_constraints():
    from types import SimpleNamespace

    import yaml

    from scripts.solve_network import create_optimization_model

    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2013-01-01", periods=8760, freq="h"))
    n.add("Carrier", "AC")
    n.add("Bus", "AC", carrier="AC")
    n.add(
        "Generator",
        "supply",
        bus="AC",
        p_nom_extendable=True,
        marginal_cost=1,
        capital_cost=8760,
    )
    n.add("Load", "demand", bus="AC", p_set=2)
    n = short_network(n, 24)
    config = yaml.safe_load(open("config/config.default.yaml"))
    create_optimization_model(
        n, config, SimpleNamespace(custom_extra_functionality=[]), {}, {}
    )
    assert n.optimize.solve_model(solver_name="highs", threads=1) == ("ok", "optimal")
    assert n.generators.at["supply", "p_nom_opt"] == 2
    assert n.objective == pytest.approx(96)


def preflight_source():
    from scripts.prepare_sector_network import gas_storage_minimum_capacities

    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2013-01-01", periods=8760, freq="h"))
    n.add("Bus", [f"zone-{i}" for i in range(52)], carrier="AC")
    n.add("Bus", ["A gas", "B gas"], carrier="gas")
    n.add("Generator", "supply", bus="A gas", p_nom=100)
    n.add(
        "Link",
        "existing pipe",
        carrier="gas pipeline",
        bus0="A gas",
        bus1="B gas",
        p_nom=100,
    )
    n.add("Load", "demand", bus="B gas", p_set=50)
    n.add(
        "Store",
        ["A gas Store", "B gas Store"],
        carrier="gas",
        bus=["A gas", "B gas"],
        e_nom_extendable=True,
    )
    gas_inputs = pd.DataFrame({"storage": [1.0, 10.0]}, index=["A", "B"])
    n.stores["e_nom_min"] = gas_storage_minimum_capacities(gas_inputs, n.stores.index)
    n.meta["sector"] = {
        "gas_network": True,
        "gas_network_connectivity_upgrade": 1,
        "H2_network": False,
        "H2_retrofit": False,
        "co2_spatial": True,
        "co2_network": True,
    }
    return n, gas_inputs


def test_preflight_accepts_rebuilt_source():
    from scripts.run_bornholm_diagnostic import check_source

    n, gas_inputs = preflight_source()
    assert (
        check_source(n, gas_inputs)["annual_gas_check"][
            "optimistic_supply_shortfall_TWh"
        ]
        == 0
    )


@pytest.mark.parametrize(
    "fault, message",
    [
        ("storage", "Stale gas-storage bounds"),
        ("zero", "zero-capacity"),
        ("bus", "Undefined buses"),
        ("bounds", "Conflicting bounds"),
    ],
)
def test_preflight_rejects_bad_source(fault, message):
    from scripts.run_bornholm_diagnostic import check_source

    n, gas_inputs = preflight_source()
    if fault == "storage":
        n.stores.at["A gas Store", "e_nom_min"] = 1
    elif fault == "zero":
        n.links.at["existing pipe", "p_nom"] = 0
    elif fault == "bus":
        n.loads.at["demand", "bus"] = "nonexistent gas bus"
    else:
        n.stores.at["A gas Store", "e_nom_max"] = 0
    with pytest.raises(ValueError, match=message):
        check_source(n, gas_inputs)
