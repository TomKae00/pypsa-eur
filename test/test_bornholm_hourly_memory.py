# SPDX-License-Identifier: MIT
"""Memory retry controls; run with python -m unittest test.test_bornholm_hourly_memory."""

import copy
import json
import os
import subprocess
import sys
import tempfile
import unittest
from contextlib import nullcontext
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from scripts import run_bornholm_hourly as hourly


class MemoryBudgetTests(unittest.TestCase):
    def test_observed_node_keeps_standard_limit_and_allows_210_retry(self):
        snapshot = {"MemTotal": 261.547442176, "MemAvailable": 234.537353216}
        self.assertAlmostEqual(hourly.memory_cap(snapshot, 190), 175.903014912)
        policy = hourly.MEMORY_PROFILES["retry"]
        cap = hourly.memory_cap(
            snapshot,
            policy["limit_gb"],
            fraction=policy["fraction"],
            reserve_gb=policy["reserve_gb"],
        )
        self.assertEqual(cap, 210)
        self.assertGreater(cap, 205.768039224)
        self.assertGreaterEqual(snapshot["MemAvailable"] - cap, 20)

    def test_retry_budget_shrinks_when_less_memory_is_available(self):
        self.assertEqual(
            hourly.memory_cap({"MemAvailable": 200}, 210, fraction=0.9, reserve_gb=20),
            180,
        )
        self.assertAlmostEqual(
            hourly.memory_cap({"MemAvailable": 50}, 210, fraction=0.9, reserve_gb=20),
            30,
        )

    def test_invalid_budgets_are_rejected(self):
        cases = [
            ({"MemAvailable": 100}, float("nan"), {}),
            ({"MemAvailable": float("inf")}, 210, {}),
            ({"MemAvailable": 100}, 0, {}),
            ({"MemAvailable": 100}, 210, {"fraction": 1}),
            ({"MemAvailable": 100}, 210, {"reserve_gb": -1}),
            ({"MemAvailable": 19}, 210, {"fraction": 0.9, "reserve_gb": 20}),
        ]
        for snapshot, requested, kwargs in cases:
            with (
                self.subTest(snapshot=snapshot, requested=requested, kwargs=kwargs),
                self.assertRaises(ValueError),
            ):
                hourly.memory_cap(snapshot, requested, **kwargs)

    def test_retry_requires_reference_before_loading_solver_dependencies(self):
        with self.assertRaisesRegex(ValueError, "requires --reference-summary"):
            hourly.solve("unused", memory_profile="retry")

    def test_driver_applies_profile_after_model_and_records_memory_failure(self):
        for profile, expected in (("standard", 175.903014912), ("retry", 210)):
            with (
                self.subTest(profile=profile),
                tempfile.TemporaryDirectory() as directory,
            ):
                root = Path(directory)
                output_path = root / "results" / hourly.PREFIX / "job" / hourly.RUN
                options = {"SoftMemLimit": 190, "TimeLimit": 3600}
                solving = {"options": {"rolling_horizon": False}}
                context = SimpleNamespace(
                    input=SimpleNamespace(network=hourly.TARGET),
                    params=SimpleNamespace(
                        solving=solving,
                        foresight="overnight",
                        co2_sequestration_potential={},
                        sector={"limit_max_growth": {}},
                        custom_extra_functionality=None,
                    ),
                    config={"solving": solving},
                    wildcards={"run": hourly.RUN},
                    log=SimpleNamespace(
                        solver=str(output_path / "solver.log"),
                        memory=str(output_path / "memory.log"),
                    ),
                )
                events = []

                def memory(events=events):
                    events.append("memory")
                    return {"MemTotal": 261.547442176, "MemAvailable": 234.537353216}

                def native_solve(
                    events=events, expected=expected, context=context, **kwargs
                ):
                    events.append("solve")
                    self.assertAlmostEqual(
                        kwargs["solver_options"]["SoftMemLimit"], expected
                    )
                    Path(context.log.solver).write_text("Memory limit reached\n")
                    return "error", "memory_limit"

                n = SimpleNamespace(
                    objective=None,
                    meta={},
                    model=SimpleNamespace(nvars=2, ncons=2, solver_model=object()),
                    optimize=SimpleNamespace(solve_model=native_solve),
                )
                modules = {
                    "numpy": SimpleNamespace(
                        random=SimpleNamespace(seed=lambda value: None)
                    ),
                    "pypsa": SimpleNamespace(
                        Network=lambda *a, network=n: network,
                        options=SimpleNamespace(
                            params=SimpleNamespace(optimize=SimpleNamespace())
                        ),
                    ),
                    "scripts._helpers": SimpleNamespace(
                        configure_logging=lambda c: None
                    ),
                    "scripts._benchmark": SimpleNamespace(
                        memory_logger=lambda **k: nullcontext()
                    ),
                    "scripts.run_bornholm_diagnostic": SimpleNamespace(
                        memory_snapshot=memory
                    ),
                    "scripts.solve_network": SimpleNamespace(
                        prepare_network=lambda *a, **k: None,
                        check_objective_value=lambda *a: None,
                        collect_kwargs=lambda *a, settings=options, **k: (
                            {},
                            {"solver_options": settings},
                        ),
                        create_optimization_model=lambda *a, events=events, **k: (
                            events.append("model")
                        ),
                    ),
                }
                with (
                    patch.dict(sys.modules, modules),
                    patch.object(hourly, "configure", return_value=options),
                    patch.object(hourly, "check_reference", return_value={}),
                    patch.object(hourly, "check_hybrid", return_value={}),
                    patch.object(hourly, "code_provenance", return_value={}),
                    patch.object(hourly.study, "context_for", return_value=context),
                    patch.object(hourly.study, "preflight", return_value={}),
                    patch.object(hourly.accuracy, "sha256", return_value="source-sha"),
                    patch.object(
                        hourly.accuracy,
                        "solver_quality",
                        return_value={"values": {"Status": 17}},
                    ),
                ):
                    # Keep the real path checks; isolate all outputs in a temp cwd.
                    previous = Path.cwd()
                    try:
                        os.chdir(root)
                        output = Path("results") / hourly.PREFIX / "job" / hourly.RUN
                        with self.assertRaisesRegex(
                            ValueError, "No finite assigned solution"
                        ):
                            hourly.solve(
                                output,
                                time_limit=3600,
                                memory_profile=profile,
                                reference_summary="reference.json"
                                if profile == "retry"
                                else None,
                            )
                        report = json.loads((output / "summary.json").read_text())
                    finally:
                        os.chdir(previous)
                self.assertEqual(events, ["memory", "model", "memory", "solve"])
                self.assertEqual(report["memory_profile"], profile)
                self.assertAlmostEqual(
                    report["memory_policy"]["effective_limit_gb"], expected
                )
                self.assertEqual(
                    report["solver_stop_reason"], "MEM_LIMIT (SoftMemLimit reached)"
                )
                self.assertEqual(report["status"], "failed")


class ReferenceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.options = {
            "threads": 16,
            "method": 2,
            "crossover": 0,
            "SoftMemLimit": 175.903014912,
            "TimeLimit": 250793.29,
        }
        self.reference = {
            "scenario": hourly.RUN,
            "source_sha256": "source-sha",
            "solver_options": self.options,
            "code_files_sha256": {},
        }
        self.physical_files = {}
        for name in dict.fromkeys([*hourly.accuracy.PINNED_FILES, *hourly.CONFIGS]):
            path = self.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("original formulation: " + name)
            self.physical_files[name] = path
            digest = hourly.accuracy.sha256(path)
            if name == "data/custom_extra_functionality.py":
                self.reference["custom_constraints_sha256"] = digest
            else:
                self.reference["code_files_sha256"][name] = digest
        self.path = self.root / "summary.json"

    def check(self, reference=None, source_hash="source-sha", options=None):
        self.path.write_text(
            json.dumps(self.reference if reference is None else reference)
        )
        original_sha = hourly.accuracy.sha256

        def physical_sha(name):
            return original_sha(self.physical_files[name])

        with patch.object(hourly.accuracy, "sha256", side_effect=physical_sha):
            return hourly.check_reference(
                self.path, source_hash, self.options if options is None else options
            )

    def test_identical_input_and_formulation_allow_only_new_budgets(self):
        options = self.options | {"SoftMemLimit": 210, "TimeLimit": 252000}
        result = self.check(options=options)
        self.assertEqual(result["reference_source_sha256"], "source-sha")

    def test_rebuilt_source_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Hourly input changed"):
            self.check(source_hash="different-source")

    def test_other_scenario_is_rejected(self):
        reference = copy.deepcopy(self.reference)
        reference["scenario"] = "bornholm_zonal_optimised"
        with self.assertRaisesRegex(ValueError, "fixed Hybrid"):
            self.check(reference=reference)

    def test_changed_threads_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "Only memory and time"):
            self.check(options=self.options | {"threads": 32})

    def test_missing_solver_options_are_rejected(self):
        reference = copy.deepcopy(self.reference)
        reference.pop("solver_options")
        with self.assertRaisesRegex(ValueError, "Only memory and time"):
            self.check(reference=reference)

    def test_changed_native_code_or_config_is_rejected(self):
        for name in self.physical_files:
            with self.subTest(name=name):
                path = self.physical_files[name]
                original = path.read_bytes()
                path.write_text("modified formulation")
                try:
                    with self.assertRaisesRegex(
                        ValueError, "physical formulation changed"
                    ):
                        self.check()
                finally:
                    path.write_bytes(original)

    def test_missing_physical_manifest_is_rejected(self):
        reference = copy.deepcopy(self.reference)
        reference.pop("code_files_sha256")
        with self.assertRaisesRegex(ValueError, "lacks provenance"):
            self.check(reference=reference)


class LauncherTests(unittest.TestCase):
    def run_launcher(self, args, *, build_fails=False, existing_source=True):
        launcher = Path(hourly.__file__).with_name("bornholm_hourly_hybrid_fatq.sbatch")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bins = root / "bin"
            bins.mkdir()
            calls = root / "calls.jsonl"
            wrapper = bins / "pixi-mock"
            wrapper.write_text(
                f"#!{sys.executable}\n"
                + """import json, os, sys
with open(os.environ['CALLS'], 'a') as stream:
    stream.write(json.dumps(sys.argv[1:])+'\\n')
if 'snakemake' in sys.argv and os.environ['BUILD_FAILS'] == '1': sys.exit(9)
"""
            )
            wrapper.chmod(0o755)
            srun = bins / "srun"
            srun.write_text(
                '#!/bin/bash\nwhile [[ "$1" == --* ]]; do shift; done\nexec "$@"\n'
            )
            srun.chmod(0o755)
            source = root / hourly.TARGET
            source.parent.mkdir(parents=True)
            if existing_source:
                source.write_bytes(b"original hourly input")
            reference = (
                root
                / "results"
                / hourly.PREFIX
                / "4937219"
                / hourly.RUN
                / "summary.json"
            )
            reference.parent.mkdir(parents=True)
            reference.write_text('{"original": true}')
            env = os.environ | {
                "PATH": f"{bins}:{os.environ['PATH']}",
                "CALLS": str(calls),
                "BUILD_FAILS": str(int(build_fails)),
                "SLURM_SUBMIT_DIR": str(root),
                "SLURM_JOB_ID": "5555555",
                "SLURM_JOB_PARTITION": "fatq",
                "SLURM_CPUS_PER_TASK": "32",
                "PIXI_SOPHIA": str(wrapper),
            }
            env.pop("SLURM_JOB_END_TIME", None)
            result = subprocess.run(
                ["bash", str(launcher), *args],
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )
            records = (
                [json.loads(line) for line in calls.read_text().splitlines()]
                if calls.exists()
                else []
            )
            self.assertEqual(reference.read_text(), '{"original": true}')
            if existing_source:
                self.assertEqual(source.read_bytes(), b"original hourly input")
            return result, records

    def test_standard_submission_still_builds_and_solves(self):
        result, calls = self.run_launcher([])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(calls), 2)
        self.assertIn("snakemake", calls[0])
        self.assertEqual(calls[0][-1], hourly.TARGET)
        self.assertEqual(calls[1][calls[1].index("--memory-profile") + 1], "standard")
        self.assertEqual(calls[1][calls[1].index("--threads") + 1], "16")

    def test_retry_skips_build_and_passes_reference_and_profile(self):
        result, calls = self.run_launcher(["16", "retry", "4937219"])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(calls), 1)
        self.assertIn("scripts.run_bornholm_hourly", calls[0])
        self.assertNotIn("snakemake", calls[0])
        self.assertEqual(calls[0][calls[0].index("--memory-profile") + 1], "retry")
        reference = calls[0][calls[0].index("--reference-summary") + 1]
        self.assertIn("/4937219/", reference)
        self.assertIn("/5555555/", calls[0][-1])

    def test_failed_standard_build_never_starts_solver(self):
        result, calls = self.run_launcher([], build_fails=True)
        self.assertEqual(result.returncode, 9)
        self.assertEqual(len(calls), 1)

    def test_retry_without_source_never_starts_build_or_solver(self):
        result, calls = self.run_launcher(
            ["16", "retry", "4937219"], existing_source=False
        )
        self.assertEqual(result.returncode, 2)
        self.assertEqual(calls, [])

    def test_retry_requires_a_valid_reference_id(self):
        for args in (["16", "retry"], ["16", "retry", "../other"], ["16", "unknown"]):
            with self.subTest(args=args):
                result, calls = self.run_launcher(args)
                self.assertEqual(result.returncode, 2)
                self.assertEqual(calls, [])


if __name__ == "__main__":
    unittest.main()
