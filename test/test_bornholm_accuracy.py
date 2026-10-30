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


def profile_settings(profile, default_focus=1):
    spec = accuracy.PROFILES[profile] if profile else {}
    return spec.get("NumericFocus", 1 if profile else default_focus), spec.get(
        "threads", 8
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
    with pytest.raises(ValueError, match="eight or 16 threads"):
        accuracy.configure_context(context(), tmp_path, 0, 4)


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


@pytest.mark.parametrize("profile", [None, *accuracy.PROFILES])
def test_live_pypsa_gurobi_quality_capture(tmp_path, profile):
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
    focus, threads = profile_settings(profile)
    options = accuracy.configure_context(
        context(), tmp_path, focus, threads, profile=profile
    )
    status, condition = n.optimize(solver_name="gurobi", **options)
    assert str(condition) == "optimal"
    n.export_to_netcdf(tmp_path / "candidate.nc")
    q = accuracy.solver_quality(n.model.solver_model)
    assert q["values"]["ObjVal"] == pytest.approx(n.objective)
    assert np.isfinite(q["values"]["MaxMemUsed"])
    assert accuracy.coal_check(n)["same_dispatch_avoidable_cost_eur"] < 1
    if options["crossover"] != 0:
        # A poor earlier log row must not reject a clean final basic solution.
        assert (
            accuracy.accuracy_issues(
                q,
                accuracy.coal_check(n),
                {"displayed_separation_eur": 1e9},
                crossover=options["crossover"],
            )
            == []
        )


@pytest.mark.parametrize(
    "excess,passes,native_failure",
    [(0, True, False), (588e6, False, False), (0, False, True)],
)
@pytest.mark.parametrize("profile", [None, *accuracy.PROFILES])
def test_native_handoff_retains_failed_candidate_and_original_outputs(
    tmp_path, monkeypatch, excess, passes, native_failure, profile
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
            if key == "Status":
                return 13 if native_failure else 2
            return 123 if key == "ObjVal" else 0

    n = SimpleNamespace(
        meta={"solve_result": {"status": "ok", "termination_condition": "optimal"}},
        objective=123,
        model=SimpleNamespace(
            solver_model=Solver(),
            status="warning" if native_failure else "ok",
            termination_condition="suboptimal" if native_failure else "optimal",
        ),
    )
    monkeypatch.setattr(
        accuracy, "coal_check", lambda _: {"same_dispatch_avoidable_cost_eur": excess}
    )
    directory = tmp_path / "results/bornholm_accuracy/pair/tol8_nf0_t8"

    def native(path, *, run_name, init_globals):
        assert init_globals["snakemake"] is ctx
        assert ctx.params.solving["solver_options"]["test"]["SoftMemLimit"] == 90
        # Simulate a poor pre-cleanup barrier row for the crossover profile.
        gap = (
            200e6 if accuracy.PROFILES.get(profile, {}).get("crossover", 0) != 0 else 0
        )
        Path(ctx.log.solver).write_text(f"9 {123 + gap:.8e} 1.23000000e+02 0 0 0 1s\n")
        if native_failure:
            exec(
                compile(
                    "raise RuntimeError(\"Solving status 'warning'. Discarding solution.\")",
                    "scripts/solve_network.py",
                    "exec",
                ),
                {"n": n},
            )
        Path(ctx.output.network).write_text("new candidate")
        return {"n": n}

    monkeypatch.setattr(accuracy.runpy, "run_path", native)
    focus, threads = profile_settings(profile, default_focus=0)
    if passes:
        accuracy.solve(focus, threads, directory, profile=profile)
    else:
        exception, match = (
            (RuntimeError, "Discarding solution")
            if native_failure
            else (ValueError, "needs review")
        )
        with pytest.raises(exception, match=match):
            accuracy.solve(focus, threads, directory, profile=profile)
    report = json.loads((directory / "summary.json").read_text())
    assert report["status"] == ("checks_passed" if passes else "failed")
    assert (directory / "candidate_2035.nc").exists() == (
        not passes and not native_failure
    )
    assert (directory / "solved_2035.nc").exists() == passes
    if native_failure:
        assert report["solver_quality"]["values"]["Status"] == 13
        assert report["native_termination_condition"] == "suboptimal"
        expected_gap = (
            200e6 if accuracy.PROFILES.get(profile, {}).get("crossover", 0) != 0 else 0
        )
        assert report["barrier_final"]["displayed_separation_eur"] == expected_gap
    if accuracy.PROFILES.get(profile, {}).get("crossover", 0) != 0:
        assert (
            report["screening_basis"] == "returned_solution_quality_crossover_enabled"
        )
    assert previous.read_text() == "previous completed solution"
    assert original.read_text() == "unchanged source"
    with pytest.raises(FileExistsError):
        accuracy.solve(focus, threads, directory, profile=profile)


@pytest.mark.parametrize(
    "profile,changed",
    [
        ("nf1", None),
        ("nf1_scale2", "ScaleFlag"),
        ("nf1_crossover", "crossover"),
        ("nf1_homogeneous", "BarHomogeneous"),
    ],
)
def test_sweep_changes_one_setting_and_preserves_physical_inputs(
    tmp_path, profile, changed
):
    base_context, trial_context = context(), context()
    base = accuracy.configure_context(
        base_context, tmp_path / "base", 1, 8, profile="nf1"
    )
    trial = accuracy.configure_context(
        trial_context, tmp_path / profile, 1, 8, profile=profile
    )
    differences = {key for key in base if base[key] != trial[key]}
    assert differences == ({changed} if changed else set())
    assert trial_context.input.network == "original.nc"
    assert trial_context.params.solving is trial_context.config["solving"]
    assert trial_context.params.solving["options"] == {"skip_iterations": True}
    assert trial["threads"] == 8 and trial["NumericFocus"] == 1
    assert trial["BarConvTol"] == 1e-8
    assert trial["TimeLimit"] == 72000 and trial["SoftMemLimit"] == 100


def test_profile_focus_mismatch_rejected(tmp_path):
    with pytest.raises(ValueError, match="NumericFocus does not match profile"):
        accuracy.configure_context(context(), tmp_path, 2, 8, profile="nf1")


@pytest.mark.parametrize("profile", accuracy.PROFILES)
def test_sweep_options_survive_real_snakemake_params_access(tmp_path, profile):
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
    focus, threads = profile_settings(profile)
    options = accuracy.configure_context(
        native, tmp_path, focus, threads, profile=profile
    )
    assert native.params.solving["solver_options"]["test"] == options
    for key, value in accuracy.PROFILES[profile].items():
        if key != "label":
            assert native.params.solving["solver_options"]["test"][key] == value


@pytest.mark.parametrize("task,profile", list(enumerate(accuracy.ORIGINAL_SWEEP)))
def test_sweep_array_handoff(tmp_path, task, profile):
    wrapper = tmp_path / "pixi-sophia"
    wrapper.write_text(
        "#!/bin/bash\npython3 -c 'import sys,json; print(json.dumps(sys.argv[1:]))' \"$@\"\n"
    )
    wrapper.chmod(0o755)
    script = Path("scripts/bornholm_solver_sweep_workq.sbatch").resolve()
    env = os.environ | {
        "SLURM_SUBMIT_DIR": str(tmp_path),
        "SLURM_JOB_PARTITION": "workq",
        "SLURM_CPUS_PER_TASK": "8",
        "SLURM_ARRAY_TASK_ID": str(task),
        "SLURM_ARRAY_JOB_ID": "456",
        "PIXI_SOPHIA": str(wrapper),
    }
    result = subprocess.run(
        ["bash", str(script)], env=env, text=True, capture_output=True
    )
    assert result.returncode == 0, result.stderr
    args = json.loads(result.stdout)
    assert args[:4] == ["run", "python", "-m", "scripts.run_bornholm_accuracy"]
    assert args[args.index("--profile") + 1] == profile
    assert (
        args[args.index("--output-dir") + 1]
        == "results/bornholm_accuracy/456/" + accuracy.PROFILES[profile]["label"]
    )


@pytest.mark.parametrize(
    "override",
    [
        {"SLURM_JOB_PARTITION": "fatq"},
        {"SLURM_CPUS_PER_TASK": "4"},
        {"SLURM_ARRAY_TASK_ID": "9"},
    ],
)
def test_sweep_rejects_wrong_allocation(tmp_path, override):
    script = Path("scripts/bornholm_solver_sweep_workq.sbatch").resolve()
    env = (
        os.environ
        | {
            "SLURM_SUBMIT_DIR": str(tmp_path),
            "SLURM_JOB_PARTITION": "workq",
            "SLURM_CPUS_PER_TASK": "8",
            "SLURM_ARRAY_TASK_ID": "0",
        }
        | override
    )
    result = subprocess.run(
        ["bash", str(script)], env=env, text=True, capture_output=True
    )
    assert result.returncode == 2


def test_compare_sweep_reports_options_failures_and_common_input(tmp_path, capsys):
    for index, (profile, spec) in enumerate(accuracy.PROFILES.items()):
        directory = tmp_path / spec["label"]
        directory.mkdir()
        focus, threads = profile_settings(profile)
        options = accuracy.configure_context(
            context(), directory, focus, threads, profile=profile
        )
        (directory / "summary.json").write_text(
            json.dumps(
                {
                    "status": "failed" if index == 0 else "checks_passed",
                    "profile": profile,
                    "source_sha256": "same",
                    "solver_options": options,
                    "error": "native solver rejected" if index == 0 else None,
                    "solver_quality": {
                        "values": {"Status": 13 if index == 0 else 2, "ObjVal": 123}
                    },
                }
            )
        )
    accuracy.compare(tmp_path)
    output = json.loads(capsys.readouterr().out)
    assert len(output["runs"]) == len(accuracy.PROFILES)
    assert output["runs"]["tol8_nf1_t8"]["Gurobi_status"] == 13
    homogeneous = output["runs"]["tol8_nf1_homogeneous_t8"]["solver_options"]
    assert homogeneous["BarConvTol"] == 1e-8
    assert homogeneous["BarHomogeneous"] == 1
    assert output["runs"]["tol8_nf1_t16"]["threads"] == 16
    assert (
        output["runs"]["tol8_nf0_homogeneous_t8"]["solver_options"]["NumericFocus"] == 0
    )
    path = tmp_path / "tol8_nf1_homogeneous_t8/summary.json"
    report = json.loads(path.read_text())
    report["source_sha256"] = "different"
    path.write_text(json.dumps(report))
    with pytest.raises(ValueError, match="inputs differ"):
        accuracy.compare(tmp_path)


def clean_quality():
    return {
        "values": {
            name: 0.0
            for name in (
                "BoundVio",
                "ConstrVio",
                "ConstrResidual",
                "DualVio",
                "DualResidual",
                "ComplVioSum",
            )
        }
        | {"Status": 2}
    }


def test_crossover_uses_returned_solution_not_last_barrier_row():
    quality = clean_quality()
    coal = {"same_dispatch_avoidable_cost_eur": 0}
    barrier = {"displayed_separation_eur": 500e6}
    assert accuracy.accuracy_issues(quality, coal, barrier, crossover=-1) == []
    assert accuracy.accuracy_issues(quality, coal, None, crossover=-1) == []
    assert any(
        "5m" in issue for issue in accuracy.accuracy_issues(quality, coal, barrier)
    )


@pytest.mark.parametrize(
    "metric,value",
    [
        ("Status", 13),
        ("BoundVio", 0.01),
        ("ConstrVio", 0.01),
        ("ConstrResidual", 0.01),
        ("DualVio", 1e-3),
        ("DualResidual", 1e-3),
        ("ComplVioSum", 5e6 + 1),
        ("ComplVioSum", None),
        ("ComplVioSum", float("nan")),
        ("ComplVioSum", float("inf")),
    ],
)
def test_crossover_still_rejects_poor_or_missing_final_quality(metric, value):
    quality = clean_quality()
    quality["values"][metric] = value
    issues = accuracy.accuracy_issues(
        quality,
        {"same_dispatch_avoidable_cost_eur": 0},
        {"displayed_separation_eur": 0},
        crossover=-1,
    )
    assert issues


@pytest.mark.parametrize("value", [2e6, None, float("nan"), float("inf")])
def test_crossover_retains_coal_cost_screen(value):
    assert accuracy.accuracy_issues(
        clean_quality(),
        {"same_dispatch_avoidable_cost_eur": value},
        {"displayed_separation_eur": 0},
        crossover=-1,
    )


@pytest.mark.parametrize(
    "trial,reference,changed",
    [
        ("nf1_t16", "nf1", "threads"),
        ("nf1_homogeneous_t16", "nf1_homogeneous", "threads"),
        ("nf1_homogeneous_crossover", "nf1_homogeneous", "crossover"),
        ("nf0_homogeneous", "nf1_homogeneous", "NumericFocus"),
    ],
)
def test_followup_changes_only_one_setting_from_existing_reference(
    tmp_path, trial, reference, changed
):
    first, second = context(), context()
    focus, threads = profile_settings(reference)
    base = accuracy.configure_context(
        first, tmp_path / reference, focus, threads, profile=reference
    )
    focus, threads = profile_settings(trial)
    options = accuracy.configure_context(
        second, tmp_path / trial, focus, threads, profile=trial
    )
    assert {key for key in options if options[key] != base[key]} == {changed}
    assert options["TimeLimit"] == base["TimeLimit"] == 72000
    assert options["SoftMemLimit"] == base["SoftMemLimit"] == 100
    assert first.input.network == second.input.network == "original.nc"
    assert second.threads == threads


def test_profile_rejects_mismatched_thread_count(tmp_path):
    with pytest.raises(ValueError, match="Thread count does not match profile"):
        accuracy.configure_context(context(), tmp_path, 1, 8, profile="nf1_t16")


@pytest.mark.parametrize("task,profile", list(enumerate(accuracy.FOLLOWUP_PROFILES)))
def test_followup_array_handoff(tmp_path, task, profile):
    wrapper = tmp_path / "pixi-sophia"
    wrapper.write_text(
        "#!/bin/bash\npython3 -c 'import sys,json; print(json.dumps(sys.argv[1:]))' \"$@\"\n"
    )
    wrapper.chmod(0o755)
    script = Path("scripts/bornholm_solver_followup_workq.sbatch").resolve()
    env = os.environ | {
        "SLURM_SUBMIT_DIR": str(tmp_path),
        "SLURM_JOB_PARTITION": "workq",
        "SLURM_CPUS_PER_TASK": "16",
        "SLURM_ARRAY_TASK_ID": str(task),
        "SLURM_ARRAY_JOB_ID": "789",
        "PIXI_SOPHIA": str(wrapper),
    }
    result = subprocess.run(
        ["bash", str(script)], env=env, text=True, capture_output=True
    )
    assert result.returncode == 0, result.stderr
    args = json.loads(result.stdout)
    assert args[args.index("--profile") + 1] == profile
    assert int(args[args.index("--threads") + 1]) == profile_settings(profile)[1]
    assert (
        args[args.index("--output-dir") + 1]
        == "results/bornholm_accuracy/789/" + accuracy.PROFILES[profile]["label"]
    )


@pytest.mark.parametrize(
    "override",
    [
        {"SLURM_JOB_PARTITION": "fatq"},
        {"SLURM_CPUS_PER_TASK": "8"},
        {"SLURM_ARRAY_TASK_ID": "4"},
    ],
)
def test_followup_rejects_wrong_allocation(tmp_path, override):
    script = Path("scripts/bornholm_solver_followup_workq.sbatch").resolve()
    env = (
        os.environ
        | {
            "SLURM_SUBMIT_DIR": str(tmp_path),
            "SLURM_JOB_PARTITION": "workq",
            "SLURM_CPUS_PER_TASK": "16",
            "SLURM_ARRAY_TASK_ID": "0",
        }
        | override
    )
    result = subprocess.run(
        ["bash", str(script)], env=env, text=True, capture_output=True
    )
    assert result.returncode == 2


@pytest.mark.parametrize("profile", accuracy.FOLLOWUP_PROFILES)
def test_cli_derives_focus_and_threads_from_followup_profile(
    tmp_path, monkeypatch, profile
):
    calls = []
    monkeypatch.setattr(
        accuracy, "solve", lambda *args, **kwargs: calls.append((args, kwargs))
    )
    monkeypatch.setattr(
        sys,
        "argv",
        ["run_bornholm_accuracy", "--profile", profile, "--output-dir", str(tmp_path)],
    )
    accuracy.main()
    assert calls == [((*profile_settings(profile), tmp_path), {"profile": profile})]


def test_compare_across_arrays_keeps_labels_and_allows_thread_experiment(
    tmp_path, capsys
):
    first, second = tmp_path / "4929796", tmp_path / "followup"
    for root, profiles in [
        (first, accuracy.ORIGINAL_SWEEP),
        (second, accuracy.FOLLOWUP_PROFILES),
    ]:
        for profile in profiles:
            spec = accuracy.PROFILES[profile]
            directory = root / spec["label"]
            directory.mkdir(parents=True)
            focus, threads = profile_settings(profile)
            options = accuracy.configure_context(
                context(), directory, focus, threads, profile=profile
            )
            (directory / "summary.json").write_text(
                json.dumps(
                    {
                        "status": "checks_passed",
                        "profile": profile,
                        "threads": threads,
                        "source_sha256": "same",
                        "solver_options": options,
                    }
                )
            )
    accuracy.compare([first, second])
    result = json.loads(capsys.readouterr().out)
    assert len(result["runs"]) == 8
    assert result["runs"]["4929796/tol8_nf1_t8"]["threads"] == 8
    assert result["runs"]["followup/tol8_nf1_t16"]["threads"] == 16
    path = second / "tol8_nf1_t16/summary.json"
    report = json.loads(path.read_text())
    report["solver_options"]["TimeLimit"] = 100
    path.write_text(json.dumps(report))
    with pytest.raises(ValueError, match="Other solver settings differ"):
        accuracy.compare([first, second])
