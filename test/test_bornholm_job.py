# SPDX-License-Identifier: MIT
"""Submission selection and guards without calling Slurm or starting a model."""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.fixture
def job(tmp_path):
    (tmp_path / "config").mkdir()
    for name in ["config.bornholm.yaml", "scenarios_bornholm.yaml"]:
        shutil.copy(Path("config") / name, tmp_path / "config" / name)
    wrapper = tmp_path / "pixi-mock"
    wrapper.write_text(
        f"#!{sys.executable}\n"
        "import json, os, subprocess, sys\n"
        "with open(os.environ['MOCK_CALLS'], 'a') as f: f.write(json.dumps(sys.argv[1:]) + '\\n')\n"
        "if sys.argv[1:4] == ['run', 'python', '-c']: sys.exit(subprocess.call([sys.executable, *sys.argv[3:]]))\n"
    )
    wrapper.chmod(0o755)
    env = os.environ | {
        "PIXI_SOPHIA": str(wrapper),
        "SLURM_JOB_ID": "test",
        "SLURM_SUBMIT_DIR": str(tmp_path),
        "SLURM_CPUS_PER_TASK": "16",
        "MOCK_CALLS": str(tmp_path / "calls.jsonl"),
    }
    env.pop("HINDCAST_NETWORK", None)
    script = Path("scripts/bornholm_sophia.sbatch").resolve()
    return tmp_path, script, env


@pytest.mark.parametrize(
    "stage,count", [("dry-run", 8), ("build", 8), ("solve", 1), ("all", 1)]
)
def test_job_defaults(job, stage, count):
    root, script, env = job
    result = subprocess.run(
        ["bash", str(script), stage], env=env, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    calls = [
        json.loads(line) for line in (root / "calls.jsonl").read_text().splitlines()
    ]
    smk = next(call for call in calls if call[:2] == ["run", "snakemake"])
    assert len(smk[smk.index("--") + 1 :]) == count
    if stage == "solve":
        assert smk[-1] == "results/bornholm_zonal_baseline/networks/solved_2035.nc"
    assert all(
        "--reference" not in call and "--require-reference" not in call
        for call in calls
    )


def test_job_rejects_unknown_case(job):
    _, script, env = job
    result = subprocess.run(
        ["bash", str(script), "dry-run", "typo"],
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert "Unknown scenario" in result.stderr


@pytest.mark.parametrize("stage", ["solve", "all"])
def test_job_requires_baseline_pilot(job, stage):
    _, script, env = job
    result = subprocess.run(
        ["bash", str(script), stage, "open"], env=env, capture_output=True, text=True
    )
    assert result.returncode == 2
    assert "baseline pilot" in result.stderr


def test_selected_build_includes_comparison_baseline(job):
    root, script, env = job
    result = subprocess.run(
        ["bash", str(script), "build", "open"], env=env, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    calls = [
        json.loads(line) for line in (root / "calls.jsonl").read_text().splitlines()
    ]
    smk = next(call for call in calls if call[:2] == ["run", "snakemake"])
    assert smk[smk.index("--") + 1 :] == [
        "resources/bornholm_zonal_baseline/networks/composed_2035.nc",
        "resources/bornholm_zonal_open/networks/composed_2035.nc",
    ]


def test_job_rejects_login_node_execution(job):
    _, script, env = job
    env.pop("SLURM_JOB_ID")
    result = subprocess.run(
        ["bash", str(script), "build"], env=env, capture_output=True, text=True
    )
    assert result.returncode == 2
    assert "compute node" in result.stderr
