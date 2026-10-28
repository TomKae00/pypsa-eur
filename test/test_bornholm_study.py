# SPDX-License-Identifier: MIT
"""Weighted aggregation, scenario sharing and host/container job submission."""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from snakemake.utils import update_config

from scripts import run_bornholm_study as study
from scripts.bornholm_workq_client import slurm_status
from scripts.check_bornholm_networks import check_networks, expected_timing
from scripts.prepare_sector_network import set_temporal_aggregation
from test.test_bornholm_zonal import comparison  # noqa: F401


def test_native_three_hour_aggregation_preserves_annual_energy(comparison, tmp_path):  # noqa: F811
    networks, config = comparison
    import yaml

    update_config(config, yaml.safe_load(Path(study.CONFIGS[1]).read_text()))
    weights = expected_timing(config)
    frame = pd.DataFrame(
        {key: weights for key in ("objective", "stores", "generators")}
    )
    path = tmp_path / "weights.csv"
    frame.to_csv(path)
    for case, n in networks.items():
        demand = pd.Series(np.tile([10, 20, 90], 2920), index=n.snapshots)
        n.loads_t.p_set["heat service"] = demand
        n.generators.loc[:, "e_sum_max"] = 123456
        n.stores_t.e_max_pu["DK gas Store"] = 1.0
        n.stores_t.e_max_pu.loc[n.snapshots[-1], "DK gas Store"] = 0.0
        before = demand.sum()
        n = set_temporal_aggregation(n, config["clustering"]["temporal"], path)
        assert len(n.snapshots) == 2920
        assert (n.snapshot_weightings.sum() == 8760).all()
        assert (
            n.snapshot_weightings.generators @ n.loads_t.p_set["heat service"] == before
        )
        assert n.generators.e_sum_max.eq(123456).all()
        assert n.stores_t.e_max_pu.iloc[-1]["DK gas Store"] == 0
        assert n.stores.at["DK gas Store", "e_nom_min"] == 1000
        networks[case] = n
    check_networks(networks, config)
    networks["open"].snapshot_weightings.loc[:, "stores"] = 1
    with pytest.raises(ValueError, match="physical snapshot weights"):
        check_networks(networks, config)


def test_native_sharing_keeps_scenario_composition_independent():
    from scripts._helpers import mock_snakemake

    runs = study.study_runs()
    assert len(set(runs)) == 8
    contexts = [
        mock_snakemake(
            "compose_network", configfiles=study.CONFIGS, run=run, horizon="2035"
        )
        for run in (runs[0], runs[-1])
    ]
    baseline, expandable = contexts
    assert baseline.input.network == expandable.input.network
    assert baseline.input.profile_onwind == expandable.input.profile_onwind
    assert baseline.input.snapshot_weightings == expandable.input.snapshot_weightings
    assert "/bornholm_zonal_baseline/" in baseline.input.tech_costs
    assert "/bornholm_zonal_open_expandable/" in expandable.input.tech_costs
    assert baseline.output[0] != expandable.output[0]
    assert baseline.params.electricity["transmission_limit"] == "vfixed"
    assert expandable.params.electricity["transmission_limit"] == "vopt"
    assert not baseline.params.energy_islands["enabled"]
    assert expandable.params.energy_islands["case"] == "open"
    assert all(c.params.clustering_temporal["averaging"] == "3h" for c in contexts)


@pytest.mark.parametrize(
    "state,exit_code,expected",
    [
        ("PENDING", "0:0", "running"),
        ("RUNNING", "0:0", "running"),
        ("COMPLETED", "0:0", "success"),
        ("COMPLETED", "1:0", "failed"),
        ("OUT_OF_MEMORY", "0:125", "failed"),
        ("TIMEOUT", "0:15", "failed"),
        ("CANCELLED by 1000", "0:15", "failed"),
    ],
)
def test_slurm_accounting_status(state, exit_code, expected):
    assert slurm_status(f"123|{state}|{exit_code}\n", "123") == expected
    assert slurm_status("", "123") == "running"


def test_container_client_submits_through_host_bridge(tmp_path):
    """Exercise the real bridge and client; fake only external Slurm commands."""
    root = Path.cwd()
    spool = tmp_path / "spool"
    spool.mkdir()
    binaries = tmp_path / "bin"
    binaries.mkdir()
    tool = (
        f"#!{sys.executable}\n"
        + """import json, os, sys
from pathlib import Path
name = Path(sys.argv[0]).name
with open(os.environ['SLURM_CALLS'], 'a') as f: f.write(json.dumps([name, *sys.argv[1:]])+'\\n')
if name == 'sbatch': print('123;cluster')
elif name == 'sacct': print('123|OUT_OF_MEMORY|0:125')
"""
    )
    for name in ("sbatch", "sacct", "scancel"):
        path = binaries / name
        path.write_text(tool)
        path.chmod(0o755)
    env = os.environ | {
        "PATH": f"{binaries}:{os.environ['PATH']}",
        "BORNHOLM_WORKQ_BRIDGE": str(spool),
        "SLURM_CALLS": str(tmp_path / "calls.jsonl"),
    }
    bridge = subprocess.Popen(
        [
            "bash",
            "scripts/bornholm_workq_bridge.sh",
            str(spool),
            "/home/a user/bin/pixi-sophia",
        ],
        env=env,
    )
    client = [sys.executable, "scripts/bornholm_workq_client.py"]
    try:
        jobscript = tmp_path / "a worker.sh"
        jobscript.write_text(
            '#!/bin/sh\n# properties = {"rule": "compose_network", "threads": 4}\n'
        )
        submitted = subprocess.run(
            [*client, "submit", str(jobscript)],
            env=env,
            capture_output=True,
            text=True,
            check=True,
        )
        assert submitted.stdout.strip() == "123"
        status = subprocess.run(
            [*client, "status", "123"],
            env=env,
            capture_output=True,
            text=True,
            check=True,
        )
        assert status.stdout.strip() == "failed"
        subprocess.run([*client, "cancel", "123"], env=env, check=True)
    finally:
        (spool / "stop").touch()
        assert bridge.wait(timeout=10) == 0
    calls = [
        json.loads(line) for line in (tmp_path / "calls.jsonl").read_text().splitlines()
    ]
    submit = calls[0]
    assert (
        "--partition=workq" in submit
        and "--exclusive" in submit
        and "--mem=0" in submit
    )
    assert "--cpus-per-task=4" in submit
    wrapped = next(a.removeprefix("--wrap=") for a in submit if a.startswith("--wrap="))
    import shlex

    assert shlex.split(wrapped) == [
        "/home/a user/bin/pixi-sophia",
        "run",
        "bash",
        str(jobscript),
    ]
    assert (spool / "job_ids.txt").read_text().strip() == "123"
    assert calls[-1] == ["scancel", "123"]
    assert Path.cwd() == root


@pytest.mark.parametrize(
    "stage,failed_check,submits",
    [
        ("all", False, True),
        ("build", False, False),
        ("all", True, False),
    ],
)
def test_controller_submits_exactly_eight_solves_after_checks(
    tmp_path, stage, failed_check, submits
):
    binaries = tmp_path / "bin"
    binaries.mkdir()
    calls = tmp_path / "calls.jsonl"
    wrapper = binaries / "pixi-mock"
    wrapper.write_text(
        f"#!{sys.executable}\n"
        + """import json, os, sys
args = sys.argv[1:]
with open(os.environ['STUDY_CALLS'], 'a') as f: f.write(json.dumps(args)+'\\n')
if '--list-build-targets' in args:
    for i in range(8): print(f'resources/study/case{i}/networks/composed_2035.nc')
if 'scripts.check_bornholm_networks' in args and os.environ['FAILED_CHECK'] == '1': sys.exit(9)
"""
    )
    wrapper.chmod(0o755)
    tool = (
        f"#!{sys.executable}\n"
        + """import json, os, sys
from pathlib import Path
with open(os.environ['STUDY_CALLS'], 'a') as f: f.write(json.dumps([Path(sys.argv[0]).name,*sys.argv[1:]])+'\\n')
if Path(sys.argv[0]).name == 'sbatch': print('999')
"""
    )
    for name in ("sbatch", "sacct", "scancel"):
        path = binaries / name
        path.write_text(tool)
        path.chmod(0o755)
    # The controller references only this host shell script outside the wrapper.
    (tmp_path / "scripts").mkdir()
    import shutil

    shutil.copy("scripts/bornholm_workq_bridge.sh", tmp_path / "scripts")
    env = os.environ | {
        "PATH": f"{binaries}:{os.environ['PATH']}",
        "PIXI_SOPHIA": str(wrapper),
        "STUDY_CALLS": str(calls),
        "SLURM_JOB_ID": "test",
        "SLURM_JOB_PARTITION": "workq",
        "SLURM_SUBMIT_DIR": str(tmp_path),
        "FAILED_CHECK": str(int(failed_check)),
    }
    result = subprocess.run(
        ["bash", str(Path("scripts/bornholm_study_workq.sbatch").resolve()), stage],
        env=env,
        capture_output=True,
        text=True,
        timeout=20,
    )
    assert result.returncode == (9 if failed_check else 0), result.stderr
    commands = [json.loads(line) for line in calls.read_text().splitlines()]
    submissions = [cmd for cmd in commands if cmd[0] == "sbatch"]
    assert len(submissions) == 2 * int(submits)
    if submits:
        assert submissions[0] == [
            "sbatch",
            "--parsable",
            "--partition=workq",
            "--array=0-7",
            "scripts/bornholm_study_solve.sbatch",
        ]
        assert commands.index(submissions[0]) > next(
            i for i, c in enumerate(commands) if "--preflight-all" in c
        )
        assert submissions[1] == [
            "sbatch",
            "--parsable",
            "--partition=workq",
            "--dependency=afterok:999",
            "--kill-on-invalid-dep=yes",
            "scripts/bornholm_study_check.sbatch",
        ]
    build = next(cmd for cmd in commands if cmd[:2] == ["run", "snakemake"])
    assert len(build[build.index("--") + 1 :]) == 8


@pytest.mark.parametrize(
    "script", ["bornholm_study_workq.sbatch", "bornholm_study_solve.sbatch"]
)
def test_study_refuses_other_partitions(tmp_path, script):
    env = os.environ | {
        "SLURM_JOB_PARTITION": "fatq",
        "SLURM_JOB_ID": "test",
        "SLURM_SUBMIT_DIR": str(tmp_path),
    }
    result = subprocess.run(
        ["bash", str(Path("scripts") / script)], env=env, capture_output=True, text=True
    )
    assert result.returncode == 2 and "workq" in result.stderr


@pytest.mark.parametrize("optimal", [True, False])
def test_native_solve_handoff_caps_memory_and_rejects_failed_exports(
    tmp_path, monkeypatch, optimal
):
    from types import SimpleNamespace

    from scripts import bornholm_provenance, run_bornholm_diagnostic

    source = tmp_path / "source.nc"
    source.write_text("source")
    output = tmp_path / "networks/solved_2035.nc"
    options = {"method": 2, "threads": 4, "SoftMemLimit": 70, "TimeLimit": 72000}
    context = SimpleNamespace(
        input=SimpleNamespace(network=str(source)),
        output=SimpleNamespace(network=str(output)),
        log=[],
        params=SimpleNamespace(
            solving={"solver": {"options": "test"}, "solver_options": {"test": options}}
        ),
    )
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(study, "context_for", lambda run: context)
    monkeypatch.setattr(study, "preflight", lambda ctx: {"snapshots": 2920})
    monkeypatch.setattr(bornholm_provenance, "code_provenance", lambda: {})
    monkeypatch.setattr(
        run_bornholm_diagnostic,
        "memory_snapshot",
        lambda: {"MemTotal": 100, "MemAvailable": 50},
    )

    def native(path, *, run_name, init_globals):
        assert path == "scripts/solve_network.py" and run_name == "__main__"
        assert init_globals["snakemake"] is context
        assert options["SoftMemLimit"] == 37.5 and options["TimeLimit"] == 72000
        output.write_text("native export")
        result = (
            {"status": "ok", "termination_condition": "optimal"}
            if optimal
            else {"status": "error", "termination_condition": "internal_solver_error"}
        )
        return {"n": SimpleNamespace(meta={"solve_result": result}, objective=123)}

    monkeypatch.setattr(study.runpy, "run_path", native)
    if optimal:
        study.solve("case")
        assert output.exists()
    else:
        with pytest.raises(ValueError, match="not optimal"):
            study.solve("case")
        assert not output.exists() and output.with_name("failed_2035.nc").exists()
    report = json.loads(
        (tmp_path / "results" / study.PREFIX / "case/summary.json").read_text()
    )
    assert report["status"] == ("complete" if optimal else "failed")
