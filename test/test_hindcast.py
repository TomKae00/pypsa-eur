# SPDX-License-Identifier: MIT
"""Regression tests for benchmark alignment, physical flows and dispatch."""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
import pytest

spec = importlib.util.spec_from_file_location(
    "hindcast", Path(__file__).parents[1] / "scripts/hindcast.py"
)
h = importlib.util.module_from_spec(spec)
spec.loader.exec_module(h)


def tiny_network():
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2023-01-01", periods=24, freq="h"))
    n.add("Carrier", "AC")
    n.add("Carrier", "DC")
    n.add("Carrier", "gas", co2_emissions=0.2)
    n.add("Bus", "a", carrier="AC", country="DK")
    n.add("Bus", "b", carrier="AC", country="DE")
    n.add(
        "Generator",
        "cheap",
        bus="a",
        carrier="gas",
        p_nom=100,
        marginal_cost=pd.Series(np.arange(24) + 10, index=n.snapshots),
    )
    n.add("Generator", "expensive", bus="b", carrier="gas", p_nom=100, marginal_cost=80)
    n.add("Load", "demand", bus="b", p_set=60)
    n.add("Link", "cable", bus0="a", bus1="b", carrier="DC", p_nom=40, efficiency=0.9)
    return n


def test_raw_metrics_keep_spikes_zeros_negative_prices():
    idx = pd.date_range("2023-01-01", periods=4, freq="h", tz="UTC")
    obs = pd.Series([0, -20, 40, 2000], index=idx)
    sim = pd.Series([0, -10, 60, 3000], index=idx)
    m, _ = h.paired_metrics(sim, obs)
    assert m["mae"] == 257.5
    assert m["smape_percent"] == pytest.approx((0 + 200 / 3 + 40 + 40) / 4)


def test_daily_aggregation_uses_identical_hours():
    idx = pd.date_range("2023-01-01", periods=24, freq="h", tz="UTC")
    obs = pd.Series(1.0, index=idx)
    sim = pd.Series(1.0, index=idx)
    obs.iloc[0] = np.nan
    sim.iloc[0] = 10000
    m, _ = h.paired_metrics(sim, obs, "daily", min_coverage=0.9)
    assert m["mae"] == 0
    assert m["coverage"] == 23 / 24
    with pytest.raises(ValueError, match="No valid"):
        h.paired_metrics(sim, obs, "daily", min_coverage=1)


def test_leap_day_alignment_does_not_shift_march():
    idx = pd.date_range(
        "2024-02-28", "2024-03-02", freq="h", inclusive="left", tz="UTC"
    )
    obs = pd.Series(np.arange(len(idx)), index=idx)
    sim = obs.loc[~((idx.month == 2) & (idx.day == 29))]
    m, _ = h.paired_metrics(sim, obs)
    assert m["mae"] == 0 and m["n_matched_hours"] == 48


def test_timezone_dst_duplicates_and_hourly_guard():
    idx = pd.date_range("2023-10-29", periods=5, freq="h", tz="Europe/Berlin")
    assert not h.utc_index(idx).has_duplicates
    with pytest.raises(ValueError, match="Duplicate"):
        h.utc_index(["2023-01-01", "2023-01-01"])
    n = tiny_network()
    n.snapshot_weightings.loc[:, :] = 6.0
    with pytest.raises(ValueError, match="weights"):
        h.hourly_network(n)


def test_reference_requires_all_zone_weights():
    idx = pd.date_range("2023-01-01", periods=2, freq="h", tz="UTC")
    prices = pd.DataFrame({"DK_1": [10, 10], "DK_2": [40, 40]}, index=idx)
    loads = pd.DataFrame({"DK_1": [1, 1], "DK_2": [3, np.nan]}, index=idx)
    result = h.reference_prices(
        prices, loads, {"DK": {"price_zones": ["DK_1", "DK_2"]}}
    )
    assert result.DK.iloc[0] == 32.5
    assert np.isnan(result.DK.iloc[1])
    with pytest.raises(KeyError):
        h.reference_prices(
            prices,
            loads.drop(columns="DK_2"),
            {"DK": {"price_zones": ["DK_1", "DK_2"]}},
        )


def test_branch_end_losses_and_reversed_orientation():
    n = tiny_network()
    n.links_t.p0 = pd.DataFrame({"cable": 40.0}, index=n.snapshots)
    n.links_t.p1 = pd.DataFrame({"cable": -36.0}, index=n.snapshots)
    flow = h.border_flows(n, {"DK": ["a"], "DE": ["b"]}, [("DK", "DE"), ("DE", "DK")])
    assert (flow["DK->DE"] == 40).all()
    assert (flow["DE->DK"] == -36).all()
    with pytest.raises(ValueError, match="overlapping"):
        h.region_buses(n, {"x": {"buses": ["a"]}, "y": {"buses": ["a"]}})


def test_shedding_uses_generator_sign():
    n = tiny_network()
    n.add("Generator", "shed", bus="b", carrier="load", p_nom=1000, sign=0.001)
    n.generators_t.p = pd.DataFrame(
        {"cheap": 10.0, "expensive": 10.0, "shed": 1000.0}, index=n.snapshots
    )
    assert (h.load_shedding(n)["b"] == 1).all()


def test_fixed_capacity_solve_and_complete_benchmark(tmp_path):
    n = tiny_network()
    src = tmp_path / "input.nc"
    dst = tmp_path / "solved.nc"
    n.export_to_netcdf(src)
    h.solve(src, dst, "highs", 1)
    solved = pypsa.Network(dst)
    assert np.allclose(solved.links_t.p0.cable, 40)
    assert np.allclose(solved.buses_t.marginal_price.b, 80)
    idx = h.hourly_network(solved)
    pd.DataFrame({"DK": np.arange(24) + 10, "DE": 80.0}, index=idx).to_csv(
        tmp_path / "prices.csv"
    )
    pd.DataFrame({"DK->DE": 40.0}, index=idx).to_csv(tmp_path / "flows.csv")
    cfg = {
        "year": 2023,
        "prices": str(tmp_path / "prices.csv"),
        "regions": {
            "DK": {"country": "DK", "price_zones": ["DK"]},
            "DE": {"country": "DE", "price_zones": ["DE"]},
        },
        "spreads": [["DK", "DE"]],
        "borders": [["DK", "DE"]],
        "flows": str(tmp_path / "flows.csv"),
        "flow_kind": "physical",
    }
    h.benchmark(dst, cfg, tmp_path / "results")
    scores = pd.read_csv(tmp_path / "results/scores.csv")
    assert np.allclose(scores.mae, 0)
    assert (tmp_path / "results/flow_DK_to_DE.png").exists()
    with pytest.raises(ValueError, match="new output"):
        h.solve(src, dst, "highs", 1)


def test_reject_capacity_expansion_and_static_costs(tmp_path):
    n = tiny_network()
    n.generators.loc["cheap", "p_nom_extendable"] = True
    src = tmp_path / "input.nc"
    n.export_to_netcdf(src)
    with pytest.raises(ValueError, match="extendable"):
        h.solve(src, tmp_path / "output.nc", "highs", 1)
    n.generators.loc["cheap", "p_nom_extendable"] = False
    n.generators_t.marginal_cost = pd.DataFrame(index=n.snapshots)
    n.export_to_netcdf(src)
    with pytest.raises(ValueError, match="No dynamic"):
        h.solve(src, tmp_path / "output.nc", "highs", 1)


def load_script(name):
    import sys

    sys.path.insert(0, str(Path(__file__).parents[1] / "scripts"))
    spec = importlib.util.spec_from_file_location(
        name, Path(__file__).parents[1] / f"scripts/{name}.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_subhourly_flow_gaps_are_not_partial_hour_means():
    f = load_script("retrieve_hindcast_flows")
    start, end = (
        pd.Timestamp("2023-01-01", tz="UTC"),
        pd.Timestamp("2023-01-01 02:00", tz="UTC"),
    )
    idx = pd.date_range(start, end, freq="15min", inclusive="left")
    raw = pd.Series(100.0, index=idx).drop(idx[1])
    hourly = f.hourly_complete(raw, start, end)
    assert np.isnan(hourly.iloc[0]) and hourly.iloc[1] == 100


def test_cost_builder_efficiency_emissions_and_missing_data():
    c = load_script("build_hindcast_costs")
    n = tiny_network()
    n.generators.efficiency = 0.5
    idx = h.hourly_network(n)
    fuel = pd.DataFrame({"gas": 20.0}, index=idx)
    co2 = pd.DataFrame({"co2": 100.0}, index=idx)
    costs = c.build(n, fuel, co2, {"gas": 3.0})
    assert np.allclose(costs, 83.0)  # (20 + 100 * 0.2) / 0.5 + 3
    fuel.iloc[0] = np.nan
    with pytest.raises(ValueError, match="Missing"):
        c.build(n, fuel, co2, {"gas": 3.0})
    fuel.iloc[0] = 20.0
    n.generators_t.efficiency = pd.DataFrame(
        -0.5, index=n.snapshots, columns=n.generators.index
    )
    with pytest.raises(ValueError, match="Invalid hourly efficiency"):
        c.build(n, fuel, co2, {"gas": 3.0})
