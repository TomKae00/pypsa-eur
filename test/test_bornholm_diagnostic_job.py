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


@pytest.mark.parametrize("variant", ["original", "simplified"])
def test_build_stage_stops_before_solving(diagnostic_job, variant):
    root, script, env = diagnostic_job
    result = subprocess.run(
        ["bash", str(script), "build", "24", "barrier", variant],
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    expected = ["build original"]
    if variant == "simplified":
        expected.append("build simplified")
        # This invocation only checks the composed pair, without optimisation.
        assert json.loads((root / "solve_args.json").read_text()) == [
            "run",
            "python",
            "scripts/check_bornholm_simplified.py",
            "--stage",
            "composed",
        ]
    else:
        assert not (root / "solve_args.json").exists()
    assert (root / "build_calls.txt").read_text().splitlines() == expected
    assert not (root / "results/bornholm_diagnostic").exists()


def test_real_build_reaches_workflow_without_git_on_host(diagnostic_job):
    import shutil

    root, _, env = diagnostic_job
    (root / "config").mkdir()
    (root / "config/config.bornholm.yaml").touch()
    shutil.copy("scripts/bornholm_provenance.py", root / "scripts")
    # Keep the ordinary backup utilities, but remove Git from the batch PATH.
    binary_dir = root / "bin"
    binary_dir.mkdir()
    for name in ["date", "mkdir", "cp"]:
        (binary_dir / name).symlink_to(shutil.which(name))
    wrapper = Path(env["PIXI_SOPHIA"])
    wrapper.write_text(
        f"#!{sys.executable}\nimport json,subprocess,sys\n"
        'open("solve_args.json", "w").write(json.dumps(sys.argv[1:]))\n'
        'if sys.argv[1:3] == ["run", "python"]:\n'
        "    sys.exit(subprocess.call([sys.executable, *sys.argv[3:]]))\n"
    )
    script = Path("scripts/bornholm_gas_capacity_rebuild.sbatch").resolve()
    result = subprocess.run(
        [shutil.which("bash"), str(script), "build"],
        env=env | {"PATH": str(binary_dir)},
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    (report_path,) = root.glob("gas_capacity_backup_*/code_provenance.json")
    report = json.loads(report_path.read_text())
    assert report["git_commit"] is None
    assert "FileNotFoundError" in report["git_provenance_error"]
    assert "scripts/bornholm_provenance.py" in report["code_files_sha256"]
    assert json.loads((root / "solve_args.json").read_text())[:2] == [
        "run",
        "snakemake",
    ]


@pytest.mark.parametrize("variant", ["original", "simplified"])
def test_annual_job_selects_full_year_without_rebuild(diagnostic_job, variant):
    root, _, env = diagnostic_job
    script = Path("scripts/bornholm_annual_barrier.sbatch").resolve()
    result = subprocess.run(
        ["bash", str(script), variant], env=env, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    assert not (root / "build_calls.txt").exists()
    args = json.loads((root / "solve_args.json").read_text())
    assert "--annual" in args
    assert args[args.index("--hours") + 1] == "8760"
    assert args[args.index("--variant") + 1] == variant
    assert args[args.index("--time-limit") + 1] == "72000"
    assert (
        args[args.index("--output-dir") + 1]
        == f"results/bornholm_annual_barrier/test-job_{variant}_8760h_barrier"
    )
