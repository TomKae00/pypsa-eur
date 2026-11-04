# SPDX-License-Identifier: MIT
"""Hourly chronology, fixed project identity and fatq submission safeguards."""

import json
import os
import subprocess
import sys
from contextlib import nullcontext
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pypsa
import pytest
from snakemake.iocontainers import Log, OutputFiles, Params

from scripts import run_bornholm_hourly as hourly


@pytest.fixture
def hybrid():
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2013-01-01", periods=8760, freq="h"))
    for bus in ("BEI_HUB", "DKE1", "DE00"):
        n.add("Bus", bus)
    n.add(
        "Generator",
        "BEI_WIND",
        bus="BEI_HUB",
        p_nom_extendable=True,
        p_nom_min=3000,
        p_nom_max=3000,
    )
    for name, target, size in (
        ("BEI_TO_DK2", "DKE1", 1200),
        ("BEI_TO_DE", "DE00", 2000),
    ):
        n.add(
            "Link",
            name,
            bus0="BEI_HUB",
            bus1=target,
            p_nom_extendable=True,
            p_nom_min=size,
            p_nom_max=size,
        )
    n.meta = {
        "energy_islands": {"enabled": True, "case": "hybrid"},
        "electricity": {"transmission_limit": "vfixed"},
    }
    return n


def test_hourly_guard_accepts_fixed_hybrid(hybrid):
    assert hourly.check_hybrid(hybrid)["DK2_export_MW"] == 1200


@pytest.mark.parametrize(
    "change", ["3h", "weights", "dispatch", "expandable", "endpoint", "capacity"]
)
def test_hourly_guard_rejects_wrong_source(hybrid, change):
    if change == "3h":
        hybrid.set_snapshots(hybrid.snapshots[::3])
    elif change == "weights":
        hybrid.snapshot_weightings.loc[:, "stores"] = 3
    elif change == "dispatch":
        hybrid.meta["energy_islands"]["case"] = "optimised"
    elif change == "expandable":
        hybrid.meta["electricity"]["transmission_limit"] = "vopt"
    elif change == "endpoint":
        hybrid.links.at["BEI_TO_DK2", "bus1"] = "DE00"
    else:
        hybrid.links.at["BEI_TO_DK2", "p_nom_max"] = 6000
    with pytest.raises(ValueError):
        hourly.check_hybrid(hybrid)


def test_gurobi_budget_reserves_memory_after_model_construction():
    # 100 decimal GB left after model construction => at most 75 GB to Gurobi.
    assert hourly.memory_cap({"MemAvailable": 100}, 190) == 75
    assert hourly.memory_cap({"MemAvailable": 260}, 190) == 190
    with pytest.raises(ValueError):
        hourly.memory_cap({"MemAvailable": 100}, float("inf"))


def test_runtime_context_uses_hourly_namespace_and_selected_resource_classes():
    from scripts._helpers import mock_snakemake

    context = hourly.study.context_for(hourly.RUN, configfiles=hourly.CONFIGS)
    assert Path(context.input.network).resolve() == Path(hourly.TARGET).resolve()
    assert context.config["run"]["name"] == [hourly.RUN]
    assert context.config["clustering"]["temporal"] == {
        "averaging": False,
        "segmentation": False,
        "representative": False,
    }
    assert context.config["energy_islands"]["case"] == "hybrid"
    assert context.config["energy_islands"]["enabled"]
    assert context.config["electricity"]["transmission_limit"] == "vfixed"
    assert context.config["run"]["shared_resources"]["policy"] is False
    compose = mock_snakemake(
        "compose_network", configfiles=hourly.CONFIGS, run=hourly.RUN, horizon="2035"
    )
    expected = f"resources/{hourly.PREFIX}/{hourly.RUN}/snapshot_weightings.csv"
    assert Path(compose.input.snapshot_weightings).resolve() == Path(expected).resolve()
    # Constructing the hourly context must not mutate the ordinary 3h defaults.
    reference = hourly.study.context_for(hourly.RUN)
    assert "/bornholm_3h_s1_w2_o3/" in str(reference.input.network)
    assert reference.config["clustering"]["temporal"]["averaging"] == "3h"
    for technology, count in hourly.reviewed.CLASSES.items():
        assert context.config["renewable"][technology]["resource_classes"] == count


@pytest.mark.parametrize("threads", [16, 32])
def test_solver_threads_are_separate_from_allocated_cpus(
    monkeypatch, tmp_path, threads
):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "32")
    monkeypatch.setattr(os, "sched_getaffinity", lambda _: set(range(32)))
    solving = {
        "solver": {"name": "gurobi", "options": "test"},
        "solver_options": {"test": {"method": 2, "crossover": 0}},
    }
    context = SimpleNamespace(
        params=Params(fromdict={"solving": solving}),
        config={"solving": solving},
        output=OutputFiles(fromdict={"network": "old.nc"}),
        log=Log(fromdict={"solver": "old.log"}),
        threads=4,
    )
    options = hourly.configure(context, tmp_path, threads, 3600)
    assert options["threads"] == threads
    assert options["NumericFocus"] == 1
    assert options["BarConvTol"] == 1e-8
    assert options["AggFill"] == -1
    assert options["TimeLimit"] == 3600
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    with pytest.raises(ValueError, match="Insufficient"):
        hourly.configure(context, tmp_path, threads, 3600)


@pytest.mark.parametrize("build_fails", [False, True])
def test_launcher_stops_on_failed_build_and_never_submits_other_scenarios(
    tmp_path, build_fails
):
    root = Path.cwd()
    bins = tmp_path / "bin"
    bins.mkdir()
    calls = tmp_path / "calls.jsonl"
    wrapper = bins / "pixi-mock"
    wrapper.write_text(
        f"#!{sys.executable}\n"
        + """import json, os, sys
with open(os.environ['CALLS'], 'a') as f: f.write(json.dumps(sys.argv[1:])+'\\n')
if 'snakemake' in sys.argv and os.environ['BUILD_FAILS'] == '1': sys.exit(9)
"""
    )
    wrapper.chmod(0o755)
    srun = bins / "srun"
    srun.write_text('#!/bin/bash\nwhile [[ "$1" == --* ]]; do shift; done\nexec "$@"\n')
    srun.chmod(0o755)
    env = os.environ | {
        "PATH": f"{bins}:{os.environ['PATH']}",
        "SLURM_SUBMIT_DIR": str(tmp_path),
        "SLURM_JOB_ID": "1234",
        "SLURM_JOB_PARTITION": "fatq",
        "SLURM_CPUS_PER_TASK": "32",
        "PIXI_SOPHIA": str(wrapper),
        "CALLS": str(calls),
        "BUILD_FAILS": str(int(build_fails)),
    }
    env.pop("SLURM_JOB_END_TIME", None)
    result = subprocess.run(
        ["bash", str(root / "scripts/bornholm_hourly_hybrid_fatq.sbatch")],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    recorded = [json.loads(line) for line in calls.read_text().splitlines()]
    assert recorded[0][-1] == hourly.TARGET
    assert "config/config.bornholm-hourly-hybrid.yaml" in recorded[0]
    assert "local" in recorded[0]
    if build_fails:
        assert result.returncode == 9
        assert len(recorded) == 1
    else:
        assert result.returncode == 0, result.stderr
        assert len(recorded) == 2
        assert recorded[1][recorded[1].index("--threads") + 1] == "16"
        assert recorded[1][-1] == f"results/{hourly.PREFIX}/1234/{hourly.RUN}"


@pytest.mark.parametrize("optimal", [True, False])
def test_solve_caps_memory_after_model_and_preserves_rejected_candidate(
    monkeypatch, tmp_path, optimal
):
    """Exercise the driver; replace the expensive native model/solver only."""
    from scripts import _benchmark, _helpers, solve_network

    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "32")
    monkeypatch.setattr(os, "sched_getaffinity", lambda _: set(range(32)))
    source = Path(hourly.TARGET)
    source.parent.mkdir(parents=True)
    source.write_bytes(b"unchanged hourly test input")
    solving = {
        "solver": {"name": "gurobi", "options": "test"},
        "solver_options": {"test": {"method": 2, "crossover": 0, "SoftMemLimit": 190}},
        "options": {"rolling_horizon": False},
        "check_objective": {"enable": False},
    }
    context = SimpleNamespace(
        input=SimpleNamespace(network=str(source)),
        params=Params(
            fromdict={
                "solving": solving,
                "foresight": "overnight",
                "co2_sequestration_potential": {},
                "sector": {"limit_max_growth": {}},
                "custom_extra_functionality": None,
            }
        ),
        config={
            "solving": solving,
            "foresight": "overnight",
            "bornholm_model_assumptions": {"coal_store_MWh": 1000},
        },
        output=OutputFiles(fromdict={"network": "old.nc"}),
        log=Log(fromdict={"solver": "old.log"}),
        threads=4,
        wildcards={"run": hourly.RUN, "horizon": "2035"},
    )
    events = []

    def memory():
        events.append("memory")
        return {"MemTotal": 256, "MemAvailable": 100 if "model" in events else 200}

    monkeypatch.setitem(
        sys.modules,
        "scripts.run_bornholm_diagnostic",
        SimpleNamespace(memory_snapshot=memory),
    )
    monkeypatch.setattr(hourly.study, "context_for", lambda *a, **k: context)
    monkeypatch.setattr(hourly.study, "preflight", lambda *a, **k: {})
    monkeypatch.setattr(hourly, "check_hybrid", lambda n: {})
    from scripts import bornholm_model_assumptions, dispatch_checks

    monkeypatch.setattr(
        bornholm_model_assumptions, "validate_assumptions", lambda *a: {}
    )
    monkeypatch.setattr(
        bornholm_model_assumptions, "check_solved_assumptions", lambda *a: {}
    )
    monkeypatch.setattr(dispatch_checks, "check", lambda *a: {})
    monkeypatch.setattr(hourly, "code_provenance", dict)
    monkeypatch.setattr(_helpers, "configure_logging", lambda c: None)
    monkeypatch.setattr(_benchmark, "memory_logger", lambda **k: nullcontext())
    monkeypatch.setattr(solve_network, "prepare_network", lambda *a, **k: None)
    monkeypatch.setattr(
        solve_network,
        "create_optimization_model",
        lambda *a, **k: events.append("model"),
    )
    native = SimpleNamespace(
        getAttr=lambda key: (2 if optimal else 13) if key == "Status" else 0
    )

    def native_solve(**kwargs):
        events.append("solve")
        assert kwargs["solver_options"]["SoftMemLimit"] == 75
        Path(kwargs["log_fn"]).write_text(" 1 100 100 0 0 0 1s\n")
        return ("ok", "optimal") if optimal else ("warning", "suboptimal")

    n = SimpleNamespace(
        objective=100,
        meta={},
        model=SimpleNamespace(nvars=2, ncons=2, solver_model=native),
        optimize=SimpleNamespace(solve_model=native_solve),
        export_to_netcdf=lambda path: Path(path).write_bytes(b"assigned solution"),
    )
    monkeypatch.setattr(pypsa, "Network", lambda *a: n)
    monkeypatch.setattr(
        hourly.accuracy,
        "coal_check",
        lambda n, **kwargs: {"same_dispatch_avoidable_cost_eur": 0},
    )
    directory = Path("results") / hourly.PREFIX / "job" / hourly.RUN
    if optimal:
        hourly.solve(directory, time_limit=3600)
    else:
        with pytest.raises(ValueError, match="Candidate needs review"):
            hourly.solve(directory, time_limit=3600)
    assert events == ["memory", "model", "memory", "solve"]
    summary = json.loads((directory / "summary.json").read_text())
    assert summary["status"] == ("checks_passed" if optimal else "needs_review")
    assert (directory / "solved_2035.nc").exists() is optimal
    assert (directory / "candidate_2035.nc").exists() is not optimal
    # Reusing a job directory must leave its previous artifacts untouched.
    with pytest.raises(FileExistsError):
        hourly.solve(directory, time_limit=3600)
