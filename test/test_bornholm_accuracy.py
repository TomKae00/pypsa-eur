# SPDX-License-Identifier: MIT
"""Isolated outputs, numerical rejection and paired workq submission."""

import copy
import json
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from snakemake.iocontainers import (
    InputFiles,
    Log,
    OutputFiles,
    Params,
    ResourceList,
    Snakemake,
    Wildcards,
)

from scripts import run_bornholm_accuracy as accuracy


def context():
    solving = {
        "solver": {"name": "gurobi", "options": "test"},
        "solver_options": {
            "test": {
                "method": 2,
                "crossover": 0,
                "threads": 4,
                "BarConvTol": 1e-5,
                "TimeLimit": 72000,
                "SoftMemLimit": 100,
            }
        },
        "options": {"skip_iterations": True},
    }
    return SimpleNamespace(
        input=SimpleNamespace(network="original.nc"),
        params=Params(fromdict={"solving": solving}),
        config={"solving": copy.deepcopy(solving)},
        output=OutputFiles(fromdict={"network": "old/solved.nc"}),
        log=Log(fromdict={"solver": "old/solver.log"}),
        threads=4,
    )


def test_pair_changes_only_focus_and_uses_independent_outputs(tmp_path):
    first, second = context(), context()
    options = [
        accuracy.configure_context(c, tmp_path / str(n), n, 8)
        for c, n in [(first, 0), (second, 2)]
    ]
    assert options[0] | {"NumericFocus": 2} == options[1]
    assert options[0]["BarConvTol"] == 1e-8
    assert options[0]["TimeLimit"] == 72000
    assert options[0]["SoftMemLimit"] == 100
    assert first.input.network == second.input.network == "original.nc"
    assert first.output.network != second.output.network
    assert first.log.solver != second.log.solver
    assert first.params.solving is first.config["solving"]


def test_invalid_thread_count_rejected(tmp_path):
    with pytest.raises(ValueError, match="eight threads"):
        accuracy.configure_context(context(), tmp_path, 0, 16)


def test_effective_options_survive_real_snakemake_params_access(tmp_path):
    source = context()
    native = Snakemake(
        InputFiles(fromdict={"network": "original.nc"}),
        source.output,
        source.params,
        Wildcards(fromdict={"run": accuracy.RUN}),
        4,
        ResourceList(),
        source.log,
        source.config,
        "solve_network",
        None,
    )
    options = accuracy.configure_context(native, tmp_path, 2, 8)
    options["SoftMemLimit"] = 90
    actual = native.params.solving["solver_options"]["test"]
    assert actual == options
    assert actual["NumericFocus"] == 2
    assert actual["threads"] == 8
    assert actual["BarConvTol"] == 1e-8
    assert actual["SoftMemLimit"] == 90


def test_logged_gap_and_missing_metrics_are_not_accepted(tmp_path):
    log = tmp_path / "solver.log"
    log.write_text("185 5.15893949e+11 5.15300685e+11 1.09 3.62e-04 1.81e-01 25376s\n")
    final = accuracy.barrier_result(log)
    assert final["displayed_separation_eur"] == 593264000
    issues = accuracy.accuracy_issues(
        {"values": {}}, {"same_dispatch_avoidable_cost_eur": 588e6}, final
    )
    assert any("coal" in i for i in issues)
    assert any("5m" in i for i in issues)
    assert any("BoundVio" in i for i in issues)


def test_analytic_coal_slack_detects_original_artifact():
    n = SimpleNamespace(
        stores=pd.DataFrame(
            {
                "e_nom_extendable": [True],
                "e_cyclic": [True],
                "standing_loss": [0],
                "e_nom_min": [0],
                "e_min_pu": [0],
                "e_max_pu": [1],
                "e_nom_opt": [5.88e9],
                "capital_cost": [0.1],
            },
            index=["EU coal Store"],
        ),
        stores_t=SimpleNamespace(
            e=pd.DataFrame({"EU coal Store": [2.39e9, 2.39e9 + 3]}),
            e_min_pu=pd.DataFrame(),
            e_max_pu=pd.DataFrame(),
            standing_loss=pd.DataFrame(),
        ),
    )
    x = accuracy.coal_check(n)
    assert x["capacity_only_avoidable_cost_eur"] == pytest.approx(349e6 - 0.3)
    assert x["same_dispatch_avoidable_cost_eur"] == pytest.approx(588e6 - 0.3)


@pytest.mark.parametrize(
    "task,focus,label", [("0", "0", "tol8_nf0_t8"), ("1", "2", "tol8_nf2_t8")]
)
def test_workq_array_handoff(tmp_path, task, focus, label):
    wrapper = tmp_path / "pixi-sophia"
    wrapper.write_text(
        "#!/bin/bash\npython3 -c 'import sys,json; print(json.dumps(sys.argv[1:]))' \"$@\"\n"
    )
    wrapper.chmod(0o755)
    script = Path("scripts/bornholm_accuracy_workq.sbatch").resolve()
    env = os.environ | {
        "SLURM_SUBMIT_DIR": str(tmp_path),
        "SLURM_JOB_PARTITION": "workq",
        "SLURM_CPUS_PER_TASK": "8",
        "SLURM_ARRAY_TASK_ID": task,
        "SLURM_ARRAY_JOB_ID": "123",
        "PIXI_SOPHIA": str(wrapper),
    }
    result = subprocess.run(
        ["bash", str(script)], env=env, text=True, capture_output=True
    )
    assert result.returncode == 0, result.stderr
    args = json.loads(result.stdout)
    assert args[:4] == ["run", "python", "-m", "scripts.run_bornholm_accuracy"]
    assert args[args.index("--numeric-focus") + 1] == focus
    assert (
        args[args.index("--output-dir") + 1] == f"results/bornholm_accuracy/123/{label}"
    )


def test_workq_only(tmp_path):
    script = Path("scripts/bornholm_accuracy_workq.sbatch").resolve()
    result = subprocess.run(
        ["bash", str(script)],
        text=True,
        capture_output=True,
        env=os.environ
        | {"SLURM_SUBMIT_DIR": str(tmp_path), "SLURM_JOB_PARTITION": "fatq"},
    )
    assert result.returncode == 2 and "workq" in result.stderr


def test_live_pypsa_gurobi_quality_capture(tmp_path):
    """Actual native solver object retained after exporting the candidate."""
    import pypsa

    pytest.importorskip("gurobipy")
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2035-01-01", periods=3, freq="3h"))
    n.add("Bus", "coal")
    n.add("Carrier", "coal")
    n.add(
        "Generator",
        "fuel",
        bus="coal",
        carrier="coal",
        p_nom_extendable=True,
        marginal_cost=7.5147,
    )
    n.add("Load", "industry", bus="coal", p_set=100)
    n.add(
        "Store",
        "EU coal Store",
        bus="coal",
        carrier="coal",
        e_nom_extendable=True,
        e_cyclic=True,
        capital_cost=0.1,
    )
    status, condition = n.optimize(
        solver_name="gurobi", method=2, crossover=0, BarConvTol=1e-8, threads=1
    )
    assert str(condition) == "optimal"
    n.export_to_netcdf(tmp_path / "candidate.nc")
    q = accuracy.solver_quality(n.model.solver_model)
    assert q["values"]["ObjVal"] == pytest.approx(n.objective)
    assert np.isfinite(q["values"]["MaxMemUsed"])
    assert accuracy.coal_check(n)["same_dispatch_avoidable_cost_eur"] < 1


@pytest.mark.parametrize("excess,passes", [(0, True), (588e6, False)])
def test_native_handoff_retains_failed_candidate_and_original_outputs(
    tmp_path, monkeypatch, excess, passes
):
    from scripts import bornholm_provenance

    monkeypatch.chdir(tmp_path)
    original = tmp_path / "original.nc"
    original.write_text("unchanged source")
    reference_dir = tmp_path / "results" / accuracy.study.PREFIX / accuracy.RUN
    reference_dir.mkdir(parents=True)
    previous = reference_dir / "networks/solved_2035.nc"
    previous.parent.mkdir()
    previous.write_text("previous completed solution")
    reference = {
        "status": "complete",
        "source_sha256": accuracy.sha256(original),
        "code_files_sha256": {},
        "objective": 515e9,
    }
    (reference_dir / "summary.json").write_text(json.dumps(reference))
    ctx = context()
    ctx.input.network = str(original)
    monkeypatch.setattr(accuracy.study, "context_for", lambda _: ctx)
    monkeypatch.setattr(accuracy.study, "preflight", lambda _: {"snapshots": 2920})
    monkeypatch.setattr(accuracy, "PINNED_FILES", ())
    monkeypatch.setattr(bornholm_provenance, "code_provenance", lambda: {})
    monkeypatch.setitem(
        sys.modules,
        "scripts.run_bornholm_diagnostic",
        SimpleNamespace(
            memory_snapshot=lambda: {"MemTotal": 135, "MemAvailable": 120},
            soft_memory_limit=lambda memory, requested: min(
                requested, 0.75 * memory["MemAvailable"]
            ),
        ),
    )

    class Solver:
        def getAttr(self, key):
            return 123 if key == "ObjVal" else 0

    n = SimpleNamespace(
        meta={"solve_result": {"status": "ok", "termination_condition": "optimal"}},
        objective=123,
        model=SimpleNamespace(solver_model=Solver()),
    )
    monkeypatch.setattr(
        accuracy, "coal_check", lambda _: {"same_dispatch_avoidable_cost_eur": excess}
    )
    directory = tmp_path / "results/bornholm_accuracy/pair/tol8_nf0_t8"

    def native(path, *, run_name, init_globals):
        assert init_globals["snakemake"] is ctx
        assert ctx.params.solving["solver_options"]["test"]["SoftMemLimit"] == 90
        Path(ctx.output.network).write_text("new candidate")
        Path(ctx.log.solver).write_text("9 1.230e+02 1.230e+02 0 0 0 1s\n")
        return {"n": n}

    monkeypatch.setattr(accuracy.runpy, "run_path", native)
    if passes:
        accuracy.solve(0, 8, directory)
    else:
        with pytest.raises(ValueError, match="needs review"):
            accuracy.solve(0, 8, directory)
    report = json.loads((directory / "summary.json").read_text())
    assert report["status"] == ("checks_passed" if passes else "failed")
    assert (directory / "candidate_2035.nc").exists() != passes
    assert (directory / "solved_2035.nc").exists() == passes
    assert previous.read_text() == "previous completed solution"
    assert original.read_text() == "unchanged source"
    with pytest.raises(FileExistsError):
        accuracy.solve(0, 8, directory)
