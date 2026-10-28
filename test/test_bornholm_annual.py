# SPDX-License-Identifier: MIT
"""Annual mode preserves physical inputs and reports memory-limited failures."""

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pypsa
import pytest

from scripts import run_bornholm_diagnostic as runner


def test_annual_selection_retains_year_budgets_costs_and_inventory():
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2013-01-01", periods=8760, freq="h"))
    n.add("Bus", "fuel")
    n.add(
        "Generator",
        "annual fuel",
        bus="fuel",
        p_nom=100,
        e_sum_min=8760,
        e_sum_max=17520,
        capital_cost=500,
    )
    n.add(
        "Store",
        "feedstock",
        bus="fuel",
        carrier="non-sequestered HVC",
        e_nom=8760,
        e_initial=8760,
        e_max_pu=pd.Series(1.0, index=n.snapshots),
    )
    n.stores_t.e_max_pu.iloc[-1, 0] = 0
    n.add("GlobalConstraint", "carbon cap", type="co2_atmosphere", constant=1e8)
    before = {
        name: getattr(n, name).copy(deep=True)
        for name in [
            "generators",
            "stores",
            "global_constraints",
            "snapshot_weightings",
        ]
    }
    n = runner.select_network(n, 8760, annual=True)
    assert len(n.snapshots) == 8760
    for name, expected in before.items():
        pd.testing.assert_frame_equal(getattr(n, name), expected)
    assert n.stores_t.e_max_pu.iloc[-1, 0] == 0
    assert n.meta["diagnostic"]["annual_fraction"] == 1
    with pytest.raises(ValueError, match="Annual mode requires"):
        runner.select_network(n, 168, annual=True)


def test_memory_budget_uses_available_ram_and_decimal_units(tmp_path):
    info = tmp_path / "meminfo"
    info.write_text("MemTotal: 125000000 kB\nMemAvailable: 100000000 kB\n")
    snapshot = runner.memory_snapshot(info)
    assert snapshot == {"MemTotal": 128, "MemAvailable": 102.4}
    assert runner.soft_memory_limit(snapshot) == pytest.approx(76.8)
    assert runner.soft_memory_limit(snapshot, 30) == 30
    assert runner.soft_memory_limit(snapshot, 200) == pytest.approx(76.8)
    with pytest.raises(ValueError):
        runner.soft_memory_limit(snapshot, float("nan"))


def test_memory_error_retains_annual_progress_summary(tmp_path, monkeypatch):
    source = tmp_path / "resources/baseline/networks/composed_2035.nc"
    source.parent.mkdir(parents=True)
    source.write_text("source network hash fixture")
    (source.parent.parent / "gas_input_locations_simplified.csv").write_text(
        "node,storage\nA,1\n"
    )
    out = tmp_path / "results"
    captured = {}

    def solve(**kwargs):
        captured.update(kwargs)
        assert json.loads((out / "summary.json").read_text())["phase"] == "solving"
        raise MemoryError("simulated annual solver allocation failure")

    network = SimpleNamespace(
        snapshots=pd.date_range("2013-01-01", periods=8760, freq="h"),
        meta={},
        model=SimpleNamespace(nvars=100, ncons=200),
        export_to_netcdf=lambda path: Path(path).write_text("network"),
        optimize=SimpleNamespace(solve_model=solve),
    )
    context = SimpleNamespace(
        input=SimpleNamespace(network=str(source)),
        config={
            "solving": {
                "solver_options": {"gurobi-default": {"method": 2, "crossover": 0}}
            }
        },
        params=SimpleNamespace(
            custom_extra_functionality=[],
            solving={"options": {}},
            foresight="overnight",
            co2_sequestration_potential={2035: 100},
            sector={"limit_max_growth": {}},
        ),
    )
    monkeypatch.setattr(runner, "mock_snakemake", lambda *a, **kw: context)
    monkeypatch.setattr(runner, "set_scenario_config", lambda *a: None)
    monkeypatch.setattr(runner, "code_provenance", lambda: {"git_commit": None})
    monkeypatch.setattr(runner.pypsa, "Network", lambda *a: network)
    monkeypatch.setattr(runner, "check_source", lambda *a: {})
    monkeypatch.setattr(runner, "prepare_network", lambda *a, **kw: None)
    monkeypatch.setattr(runner, "collect_kwargs", lambda *a, **kw: ({}, {}))
    monkeypatch.setattr(runner, "create_optimization_model", lambda *a, **kw: None)
    monkeypatch.setattr(
        runner, "memory_snapshot", lambda: {"MemTotal": 128, "MemAvailable": 100}
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "runner",
            "--annual",
            "--hours",
            "8760",
            "--time-limit",
            "72000",
            "--output-dir",
            str(out),
        ],
    )
    with pytest.raises(MemoryError):
        runner.main()
    report = json.loads((out / "summary.json").read_text())
    assert report["status"] == "failed" and report["phase"] == "failed"
    assert report["failed_phase"] == "solving"
    assert report["diagnostic"]["annual_fraction"] == 1
    assert captured["solver_options"]["SoftMemLimit"] == 75
    assert captured["solver_options"]["TimeLimit"] == 72000
    assert report["solver_options"]["SoftMemLimit"] == 75
    assert not (out / "solved.nc").exists()
