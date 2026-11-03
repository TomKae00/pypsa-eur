# SPDX-License-Identifier: MIT
"""Rebuild, preflight and solve the eight cases with one corrected common base."""

import argparse
import json
import logging
import os
import resource
import runpy
import socket
import time
from pathlib import Path

from scripts import run_bornholm_accuracy as accuracy
from scripts import run_bornholm_study as study
from scripts.bornholm_model_assumptions import (
    AUDIT_KEY,
    apply_resource_and_coal_assumptions,
    check_solved_assumptions,
    install_transfer_assumptions,
    resource_report,
    validate_assumptions,
)

CONFIGS = [*study.CONFIGS, "config/config.bornholm-calibrated-base.yaml"]
PREFIX = "bornholm_3h_s1_w3_o3_calibrated"
CLASSES = {**study.CLASSES, "onwind": 3}


def preview(network_path):
    """Read an old network and preview factors in memory; never write a network."""
    import pypsa
    import yaml

    n = pypsa.Network(network_path)
    original_classes = (
        n.meta.get("renewable", {}).get("onwind", {}).get("resource_classes")
    )
    settings = yaml.safe_load(Path(CONFIGS[-1]).read_text())[
        "bornholm_model_assumptions"
    ]
    audit = apply_resource_and_coal_assumptions(n, settings)
    install_transfer_assumptions(n, settings, 2035)
    return {
        "source": str(network_path),
        "source_sha256": accuracy.sha256(network_path),
        "source_onshore_classes": original_classes,
        "note": "Old-profile multiplier preview only. The three classes require a full rebuild. "
        "Upward clipping of already averaged profiles differs from hourly clipping. "
        "No optimisation or network export is performed.",
        "assumptions": audit,
    }


def context_for(run):
    return study.context_for(run, configfiles=CONFIGS)


def preflight(context):
    import pypsa

    report = study.preflight(
        context, configfiles=CONFIGS, prefix=PREFIX, classes=CLASSES
    )
    n = pypsa.Network(context.input.network)
    settings = context.config["bornholm_model_assumptions"]
    audit = validate_assumptions(n, settings)
    study.require(
        audit["transfers"]["planning_year"] == 2035,
        "Transfer limits must use the planning year, not weather year",
    )
    report.update(
        assumption_fingerprint=audit["fingerprint"],
        onshore=resource_report(n),
        coal=audit["coal"],
        transfers=audit["transfers"],
    )
    return report


def solve(run):
    import numpy as np
    import pypsa

    from scripts.bornholm_provenance import code_provenance
    from scripts.run_bornholm_diagnostic import memory_snapshot, soft_memory_limit

    context = context_for(run)
    source, output = Path(context.input.network), Path(context.output.network)
    study.require(not output.exists(), f"Accepted output already exists: {output}")
    directory = Path("results") / PREFIX / run
    attempt = (
        directory / "attempts" / (os.environ.get("SLURM_JOB_ID") or str(time.time_ns()))
    )
    attempt.mkdir(parents=True, exist_ok=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    report = {
        "status": "running",
        "phase": "preflight",
        "scenario": run,
        "node": socket.gethostname(),
        "source": str(source),
        "attempt": str(attempt),
        **code_provenance(),
    }

    def checkpoint():
        report["elapsed_seconds"] = time.monotonic() - started
        report["peak_rss_GiB"] = (
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2
        )
        for path in (attempt / "summary.json", directory / "summary.json"):
            temporary = path.with_suffix(".json.tmp")
            temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
            temporary.replace(path)

    try:
        checkpoint()
        report["preflight"] = preflight(context)
        report["source_sha256"] = accuracy.sha256(source)
        solving = context.params.solving
        options = solving["solver_options"][solving["solver"]["options"]]
        study.require(
            options["threads"] == 16
            and options["method"] == 2
            and options["crossover"] == 0,
            "Expected 16-thread barrier preset",
        )
        allocated = int(os.environ.get("SLURM_CPUS_PER_TASK", "0"))
        affinity = (
            len(os.sched_getaffinity(0))
            if hasattr(os, "sched_getaffinity")
            else allocated
        )
        study.require(
            allocated >= 16 and affinity >= 16, "Allocate and bind 16 solver CPUs"
        )
        snapshot = memory_snapshot()
        options["SoftMemLimit"] = soft_memory_limit(snapshot, options["SoftMemLimit"])
        # Snakemake Params is reconstructed on access: retain the mutated mapping.
        context.config["solving"] = solving
        context.threads = 16
        context.output = type(context.output)(
            fromdict={
                "network": str(attempt / "candidate_2035.nc"),
                "config": str(attempt / "effective_config.yaml"),
            }
        )
        context.log = type(context.log)(
            fromdict={
                key: str(attempt / f"{key}.log")
                for key in ("python", "solver", "memory")
            }
        )
        # Retain the study's native objective accounting and all supplementary constraints.
        pypsa.options.params.optimize.include_objective_constant = True
        report.update(
            solver_options=dict(options),
            node_memory_before_solve_GB=snapshot,
            include_objective_constant=True,
            phase="native_solve",
        )
        checkpoint()
        namespace = runpy.run_path(
            "scripts/solve_network.py",
            run_name="__main__",
            init_globals={"snakemake": context},
        )
        n = namespace["n"]
        report.update(
            phase="quality_checks", solve_result=n.meta.get("solve_result", {})
        )
        study.require(
            report["solve_result"]
            == {"status": "ok", "termination_condition": "optimal"},
            f"Run was not optimal: {report['solve_result']}",
        )
        study.require(
            n.objective is not None and np.isfinite(n.objective), "Non-finite objective"
        )
        report["objective"] = float(n.objective)
        report["objective_constant"] = float(n.objective_constant)
        report["solver_quality"] = accuracy.solver_quality(n.model.solver_model)
        report["barrier_final"] = accuracy.barrier_result(attempt / "solver.log")
        settings = context.config["bornholm_model_assumptions"]
        report["coal"] = accuracy.coal_check(n, fixed_MWh=settings["coal_store_MWh"])
        report["assumption_checks"] = check_solved_assumptions(n, settings)
        report["screening_issues"] = accuracy.accuracy_issues(
            report["solver_quality"],
            report["coal"],
            report["barrier_final"],
            crossover=0,
        )
        checkpoint()
        study.require(
            not report["screening_issues"],
            f"Numerical result needs review: {report['screening_issues']}",
        )
        n.meta["calibrated_base_quality"] = {
            "screening_passed": True,
            "assumption_fingerprint": n.meta[AUDIT_KEY]["fingerprint"],
            "note": "Diagnostic checks; not an independent optimality certificate.",
        }
        # Write the accepted target only after all numerical/physical screens pass.
        n.export_to_netcdf(output)
        report.update(
            status="checks_passed", phase="complete", accepted_network=str(output)
        )
    except Exception as error:
        report.update(
            status="failed",
            failed_phase=report["phase"],
            phase="failed",
            error=f"{type(error).__name__}: {error}",
        )
        if output.exists():
            output.replace(attempt / "rejected_2035.nc")
        raise
    finally:
        checkpoint()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--list-build-targets", action="store_true")
    action.add_argument("--preflight-all", action="store_true")
    action.add_argument("--index", type=int, choices=range(8))
    action.add_argument("--preview-network", type=Path)
    parser.add_argument("--preview-output", type=Path)
    args = parser.parse_args()
    runs = study.study_runs()
    if args.preview_network:
        report = (
            json.dumps(preview(args.preview_network), indent=2, allow_nan=False) + "\n"
        )
        if args.preview_output:
            args.preview_output.parent.mkdir(parents=True, exist_ok=True)
            args.preview_output.write_text(report)
            print(args.preview_output)
        else:
            print(report, end="")
    elif args.list_build_targets:
        for run in runs:
            print(f"resources/{PREFIX}/{run}/networks/composed_2035.nc")
    elif args.preflight_all:
        for run in runs:
            context = context_for(run)
            report = preflight(context)
            directory = Path("results") / PREFIX / run
            directory.mkdir(parents=True, exist_ok=True)
            (directory / "preflight.json").write_text(
                json.dumps(report, indent=2) + "\n"
            )
            print(run, json.dumps(report))
    else:
        logging.basicConfig(level=logging.INFO)
        solve(runs[args.index])


if __name__ == "__main__":
    main()
