# SPDX-License-Identifier: MIT
"""Diagnostic submission preserves build failure and passes solver/window choices."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.fixture
def diagnostic_job(tmp_path):
    (tmp_path / "scripts").mkdir()
    (tmp_path / "scripts/bornholm_gas_capacity_rebuild.sbatch").write_text(
        '#!/bin/bash\nprintf "%s\\n" "$1" > build_stage.txt\nprintf "%s\\n" "$*" >> build_calls.txt\nexit "${MOCK_BUILD_EXIT:-0}"\n'
    )
    wrapper = tmp_path / "pixi-mock"
    wrapper.write_text(
        f"#!{sys.executable}\nimport json,sys\n"
        'open("solve_args.json", "w").write(json.dumps(sys.argv[1:]))\n'
    )
    wrapper.chmod(0o755)
    env = os.environ | {
        "PIXI_SOPHIA": str(wrapper),
        "SLURM_JOB_ID": "test-job",
        "SLURM_SUBMIT_DIR": str(tmp_path),
        "SLURM_CPUS_PER_TASK": "4",
    }
    return tmp_path, Path("scripts/bornholm_diagnostic.sbatch").resolve(), env


def test_diagnostic_build_then_solves_default_window(diagnostic_job):
    root, script, env = diagnostic_job
    result = subprocess.run(
        ["bash", str(script)], env=env, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    assert (root / "build_stage.txt").read_text().strip() == "build"
    args = json.loads((root / "solve_args.json").read_text())
    assert args[:4] == ["run", "python", "-m", "scripts.run_bornholm_diagnostic"]
    assert args[args.index("--hours") + 1] == "24"
    assert args[args.index("--method") + 1] == "barrier"
    assert args[args.index("--time-limit") + 1] == "3600"


def test_diagnostic_does_not_solve_after_build_failure(diagnostic_job):
    root, script, env = diagnostic_job
    result = subprocess.run(
        ["bash", str(script)], env=env | {"MOCK_BUILD_EXIT": "7"}, capture_output=True
    )
    assert result.returncode == 7
    assert not (root / "solve_args.json").exists()


def test_diagnostic_solve_only_supports_pdhg_week(diagnostic_job):
    root, script, env = diagnostic_job
    result = subprocess.run(
        ["bash", str(script), "solve", "168", "pdhg"], env=env, capture_output=True
    )
    assert result.returncode == 0
    assert not (root / "build_stage.txt").exists()
    args = json.loads((root / "solve_args.json").read_text())
    assert args[args.index("--hours") + 1] == "168"
    assert args[args.index("--method") + 1] == "pdhg"


def test_diagnostic_dry_run_does_not_solve(diagnostic_job):
    root, script, env = diagnostic_job
    assert (
        subprocess.run(
            ["bash", str(script), "dry-run"], env=env, capture_output=True
        ).returncode
        == 0
    )
    assert (root / "build_stage.txt").read_text().strip() == "dry-run"
    assert not (root / "solve_args.json").exists()


def test_simplified_rebuilds_both_sources_and_selects_overlay(diagnostic_job):
    root, script, env = diagnostic_job
    result = subprocess.run(
        ["bash", str(script), "all", "24", "barrier", "simplified"],
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert (root / "build_calls.txt").read_text().splitlines() == [
        "build original",
        "build simplified",
    ]
    args = json.loads((root / "solve_args.json").read_text())
    assert args[args.index("--variant") + 1] == "simplified"
    assert args[args.index("--output-dir") + 1].endswith(
        "test-job_simplified_24h_barrier"
    )


@pytest.mark.parametrize(
    "variant, expected",
    [
        ("original", "bornholm_zonal_baseline"),
        ("simplified", "bornholm_simplified_baseline"),
    ],
)
def test_real_rebuild_helper_uses_selected_resources(diagnostic_job, variant, expected):
    root, _, env = diagnostic_job
    (root / "config").mkdir()
    (root / "config/config.bornholm.yaml").touch()
    (root / "config/config.bornholm-simplified.yaml").touch()
    script = Path("scripts/bornholm_gas_capacity_rebuild.sbatch").resolve()
    result = subprocess.run(
        ["bash", str(script), "dry-run", variant],
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    args = json.loads((root / "solve_args.json").read_text())
    assert args[:2] == ["run", "snakemake"]
    assert args[-1] == f"resources/{expected}/networks/composed_2035.nc"
    assert ("config/config.bornholm-simplified.yaml" in args) == (
        variant == "simplified"
    )
