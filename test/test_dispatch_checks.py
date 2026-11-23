# SPDX-License-Identifier: MIT
import pandas as pd
import pypsa
import pytest

from scripts.dispatch_checks import check, storage_check


def test_lossy_bus_ports_and_corrupted_result():
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2023-01-01", periods=3, freq="h"))
    n.add("Bus", "source")
    n.add("Bus", "sink")
    n.add("Link", "cable", bus0="source", bus1="sink", p_nom=100, efficiency=0.9)
    n.add("Generator", "supply", bus="source", p_nom=100, marginal_cost=1)
    n.add("Load", "demand", bus="sink", p_set=9)
    assert n.optimize(solver_name="highs", include_objective_constant=False) == (
        "ok",
        "optimal",
    )
    assert check(n)["maximum_absolute_bus_balance_MW"] < 1e-6
    n.links_t.p1.loc[n.snapshots[1], "cable"] += 1
    with pytest.raises(ValueError, match="bus imbalance"):
        check(n)


def test_storage_equation_includes_weighting_and_cyclic_boundary():
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2023-01-01", periods=3, freq="3h"))
    n.snapshot_weightings.stores = 3.0
    n.add("Bus", "b")
    n.add("Store", "inventory", bus="b", e_nom=100, e_cyclic=True)
    n.add("Generator", "supply", bus="b", p_nom=100)
    n.add("Load", "demand", bus="b", p_set=20)
    n.stores_t.e = pd.DataFrame({"inventory": [30.0, 60.0, 0.0]}, index=n.snapshots)
    n.stores_t.p = pd.DataFrame({"inventory": [-10.0, -10.0, 20.0]}, index=n.snapshots)
    n.generators_t.p["supply"] = [30.0, 30.0, 0.0]
    n.loads_t.p["demand"] = 20.0
    n.buses_t.marginal_price["b"] = 1.0
    assert storage_check(n) == 0

    n.stores_t.e.iloc[0, 0] += 0.046701
    report = check(n)
    assert report["maximum_storage_equation_residual_MWh"] == pytest.approx(0.046701)
    assert report["storage_tolerance_MWh"] == 0.1
    assert report["absolute_tolerance_MW"] == 0.01

    n.stores_t.e.iloc[0, 0] += 1
    with pytest.raises(ValueError, match="storage balance"):
        check(n)
