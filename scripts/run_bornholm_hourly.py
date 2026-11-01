# SPDX-License-Identifier: MIT
"""Build-target reporting and one full-year fixed-Hybrid solve on Sophia fatq."""

import argparse
import json
import logging
import math
import os
import resource
import socket
import time
from pathlib import Path

import yaml

from scripts import run_bornholm_accuracy as accuracy
from scripts import run_bornholm_study as study
from scripts.bornholm_provenance import code_provenance

RUN = "bornholm_zonal_hybrid"
PREFIX = "bornholm_1h_s1_w2_o3"
CONFIGS = [*study.CONFIGS, "config/config.bornholm-hourly-hybrid.yaml"]
TARGET = f"resources/{PREFIX}/{RUN}/networks/composed_2035.nc"
MEMORY_PROFILES = {
    "standard": {"limit_gb": None, "fraction": 0.75, "reserve_gb": 0.0},
    "retry": {"limit_gb": 210.0, "fraction": 0.90, "reserve_gb": 20.0},
    "aggfill0": {"limit_gb": 210.0, "fraction": 0.90, "reserve_gb": 20.0},
}


def check_hybrid(n):
    """Reject stale hourly inputs or an accidentally optimised project."""
    import numpy as np
    import pandas as pd

    study.require(
        n.snapshots.equals(pd.date_range("2013-01-01", periods=8760, freq="h"))
        and np.allclose(n.snapshot_weightings, 1),
        "Expected 8760 consecutive hours with unit objective/generator/store weights",
    )
    project = n.meta.get("energy_islands", {})
    study.require(
        project.get("enabled") is True and project.get("case") == "hybrid",
        "Expected the fixed Hybrid project",
    )
    study.require(
        n.meta.get("electricity", {}).get("transmission_limit") == "vfixed"
        and not n.lines.s_nom_extendable.any(),
        "This runner requires fixed mainland transmission and a single native solve",
    )
    for frame, name, capacity, target in (
        (n.generators, "BEI_WIND", 3000, "BEI_HUB"),
        (n.links, "BEI_TO_DK2", 1200, "DKE1"),
        (n.links, "BEI_TO_DE", 2000, "DE00"),
    ):
        study.require(name in frame.index, f"Missing fixed project asset: {name}")
        asset = frame.loc[name]
        fixed = (
            np.isclose(asset.p_nom_min, capacity)
            and np.isclose(asset.p_nom_max, capacity)
            if asset.p_nom_extendable
            else np.isclose(asset.p_nom, capacity)
        )
        study.require(fixed, f"Wrong fixed capacity for {name}; expected {capacity} MW")
        if name == "BEI_WIND":
            study.require(asset.bus == target, "Wind must connect to the offshore hub")
        else:
            study.require(
                asset.bus0 == "BEI_HUB" and asset.bus1 == target,
                f"Wrong endpoint for {name}",
            )
    return {"wind_MW": 3000, "DK2_export_MW": 1200, "DE_export_MW": 2000}


def configure(context, directory, threads, time_limit):
    study.require(threads in (16, 32), "Use 16 or 32 solver threads")
    allocated = int(os.environ.get("SLURM_CPUS_PER_TASK", "0"))
    affinity = (
        len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else allocated
    )
    study.require(
        allocated >= threads and affinity >= threads,
        "Insufficient allocated/bound CPUs",
    )
    study.require(0 < time_limit <= 252000, "Solver time limit must be within 70 hours")
    options = accuracy.configure_context(
        context, directory, 1, threads, profile=f"nf1_aggfill_t{threads}"
    )
    options["TimeLimit"] = time_limit
    return options


def memory_cap(snapshot, requested, *, fraction=0.75, reserve_gb=0.0):
    """Budget decimal GB from RAM available after Linopy model construction."""
    available = float(snapshot["MemAvailable"])
    study.require(
        math.isfinite(requested)
        and requested > 0
        and math.isfinite(available)
        and available > 0,
        "Positive finite memory budget required",
    )
    study.require(
        math.isfinite(fraction)
        and 0 < fraction <= 0.90
        and math.isfinite(reserve_gb)
        and reserve_gb >= 0,
        "Use at most 90% of available memory and a non-negative finite reserve",
    )
    cap = min(requested, fraction * available, available - reserve_gb)
    study.require(cap > 0, "Insufficient available RAM for the memory reserve")
    return cap


def check_reference(path, source_hash, options, *, aggfill_zero=False):
    """Preserve input/formulation; allow only the explicitly selected solver trial."""
    reference = json.loads(Path(path).read_text())
    study.require(reference.get("scenario") == RUN, "Reference must be fixed Hybrid")
    study.require(
        reference.get("source_sha256") == source_hash,
        "Hourly input changed since the reference job; do not rebuild for this retry",
    )
    excluded = {"TimeLimit", "SoftMemLimit"}
    previous = {
        k: v
        for k, v in reference.get("solver_options", {}).items()
        if k not in excluded
    }
    if aggfill_zero:
        study.require(options.get("AggFill") == 0, "AggFill trial requires AggFill=0")
        previous["AggFill"] = 0
    current = {k: v for k, v in options.items() if k not in excluded}
    study.require(
        previous == current,
        "Only memory and time budgets may differ from the reference"
        + ("; this trial also permits AggFill=0" if aggfill_zero else ""),
    )
    hashes = reference.get("code_files_sha256", {})
    for name in dict.fromkeys([*accuracy.PINNED_FILES, *CONFIGS]):
        expected = (
            reference.get("custom_constraints_sha256")
            if name == "data/custom_extra_functionality.py"
            else hashes.get(name)
        )
        study.require(
            expected is not None and accuracy.sha256(name) == expected,
            f"Reference physical formulation changed or lacks provenance: {name}",
        )
    return {
        "reference_summary": str(Path(path).resolve()),
        "reference_source_sha256": source_hash,
    }


def solve(
    directory,
    threads=16,
    time_limit=252000,
    *,
    memory_profile="standard",
    reference_summary=None,
):
    study.require(memory_profile in MEMORY_PROFILES, "Unknown memory profile")
    study.require(
        memory_profile == "standard" or reference_summary is not None,
        "Memory retry requires --reference-summary from the previous job",
    )
    import numpy as np
    import pypsa
    from scripts._benchmark import memory_logger
    from scripts._helpers import configure_logging
    from scripts.run_bornholm_diagnostic import memory_snapshot
    from scripts.solve_network import (
        check_objective_value,
        collect_kwargs,
        create_optimization_model,
        prepare_network,
    )

    directory = Path(directory).resolve()
    root = (Path("results") / PREFIX).resolve()
    study.require(
        directory.is_relative_to(root) and directory != root,
        "Use a new hourly results subdirectory",
    )
    directory.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    report = {
        "status": "running",
        "phase": "initialising",
        "scenario": RUN,
        "node": socket.gethostname(),
        "threads": threads,
        "memory_profile": memory_profile,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "slurm_cpus_per_task": os.environ.get("SLURM_CPUS_PER_TASK"),
        "include_objective_constant": True,
        "quality_screening_note": "Diagnostic screening, not an independent optimality certificate.",
        **code_provenance(),
    }

    def checkpoint(phase=None):
        if phase is not None:
            report["phase"] = phase
        report["elapsed_seconds"] = time.monotonic() - started
        report["peak_rss_GiB"] = (
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2
        )
        temporary = directory / "summary.json.tmp"
        temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        temporary.replace(directory / "summary.json")

    n = None
    try:
        checkpoint()
        context = study.context_for(RUN, configfiles=CONFIGS)
        source = Path(context.input.network)
        study.require(
            source.resolve() == Path(TARGET).resolve(), "Wrong hourly source path"
        )
        options = configure(context, directory, threads, time_limit)
        if memory_profile == "aggfill0":
            options["AggFill"] = 0
        configure_logging(context)
        # Preserve the cost convention used by the selected eight-case campaign.
        pypsa.options.params.optimize.include_objective_constant = True
        report.update(source=str(source), source_sha256=accuracy.sha256(source))
        if reference_summary is not None:
            report.update(
                check_reference(
                    reference_summary, report["source_sha256"], options,
                    aggfill_zero=memory_profile == "aggfill0",
                )
            )
        checkpoint("preflight")
        report["preflight"] = study.preflight(
            context, configfiles=CONFIGS, prefix=PREFIX, step_hours=1
        )
        n = pypsa.Network(source)
        report["fixed_hybrid"] = check_hybrid(n)
        solving = context.params.solving
        solve_opts = solving["options"]
        study.require(
            not solve_opts.get("rolling_horizon") and not solve_opts.get("nhours"),
            "Full-year chronological solve required",
        )
        study.require(
            not solve_opts.get("load_shedding", {}).get("enable")
            and not solve_opts.get("load_sinks", {}).get("enable"),
            "Do not introduce balancing slack",
        )
        hook = context.params.custom_extra_functionality
        if hook:
            study.require(Path(hook).is_file(), f"Missing custom constraints: {hook}")
            report["custom_constraints_sha256"] = accuracy.sha256(hook)
        report["solver_options"] = dict(options)
        report["node_memory_before_model_GB"] = memory_snapshot()
        (directory / "effective_config.yaml").write_text(yaml.safe_dump(context.config))
        checkpoint("building_model")
        np.random.seed(solve_opts.get("seed", 123))
        # Use the same native preparation, objective, model and supplementary
        # constraints as solve_network.py's fixed-grid single-pass branch.
        with memory_logger(filename=context.log.memory, interval=30):
            prepare_network(
                n,
                config=context.config,
                solve_opts=solve_opts,
                foresight=context.params.foresight,
                planning_horizons="2035",
                co2_sequestration_potential=context.params.co2_sequestration_potential,
                limit_max_growth=context.params.sector["limit_max_growth"],
                rolling_horizon=False,
            )
            model_kwargs, solve_kwargs = collect_kwargs(
                context.config,
                solving,
                "2035",
                log_fn=context.log.solver,
                mode="single",
            )
            create_optimization_model(
                n,
                config=context.config,
                params=context.params,
                model_kwargs=model_kwargs,
                solve_kwargs=solve_kwargs,
                planning_horizons="2035",
                snakemake=context,
            )
            # Model construction itself can be large. Measure AFTER it, rather
            # than offering Gurobi memory already occupied by Python/Linopy.
            memory = memory_snapshot()
            policy = MEMORY_PROFILES[memory_profile]
            requested = (
                policy["limit_gb"]
                if policy["limit_gb"] is not None
                else options["SoftMemLimit"]
            )
            cap = memory_cap(
                memory,
                requested,
                fraction=policy["fraction"],
                reserve_gb=policy["reserve_gb"],
            )
            solve_kwargs["solver_options"]["SoftMemLimit"] = cap
            options["SoftMemLimit"] = cap
            remaining = time_limit - (time.monotonic() - started)
            study.require(
                remaining > 60, "No solver time remains after model construction"
            )
            solve_kwargs["solver_options"]["TimeLimit"] = remaining
            options["TimeLimit"] = remaining
            report.update(
                solver_options=dict(solve_kwargs["solver_options"]),
                node_memory_before_solver_GB=memory,
                memory_policy={
                    **policy,
                    "requested_gb": requested,
                    "effective_limit_gb": cap,
                },
                variables=int(n.model.nvars),
                constraints=int(n.model.ncons),
            )
            (directory / "effective_config.yaml").write_text(
                yaml.safe_dump(context.config)
            )
            checkpoint("native_solve")
            print(
                f"Hourly fixed Hybrid: threads={threads}, memory profile={memory_profile}, "
                f"Gurobi cap={cap:.2f} GB; results={directory}",
                flush=True,
            )
            status, condition = n.optimize.solve_model(**solve_kwargs)
        report["solve_result"] = {
            "status": getattr(status, "value", status),
            "termination_condition": getattr(condition, "value", condition),
        }
        n.meta = dict(context.config, wildcards=dict(context.wildcards))
        n.meta["solve_result"] = report["solve_result"]
        report["solver_quality"] = accuracy.solver_quality(n.model.solver_model)
        report["barrier_final"] = accuracy.barrier_result(directory / "solver.log")
        if report["solver_quality"]["values"].get("Status") == 17:
            report["solver_stop_reason"] = "MEM_LIMIT (SoftMemLimit reached)"
        checkpoint("quality_checks")
        study.require(
            n.objective is not None and math.isfinite(float(n.objective)),
            "No finite assigned solution; inspect solver status/log",
        )
        report["objective"] = float(n.objective)
        # Preserve an assigned solution before any acceptance screen can fail.
        n.export_to_netcdf(directory / "candidate_2035.nc")
        report["coal"] = accuracy.coal_check(n)
        report["screening_issues"] = accuracy.accuracy_issues(
            report["solver_quality"],
            report["coal"],
            report["barrier_final"],
            crossover=0,
        )
        check_objective_value(n, solving, "2035")
        study.require(
            report["solve_result"]
            == {"status": "ok", "termination_condition": "optimal"}
            and not report["screening_issues"],
            "Candidate needs review; see summary.json and solver.log",
        )
        (directory / "candidate_2035.nc").replace(directory / "solved_2035.nc")
        report.update(status="checks_passed", phase="complete")
    except Exception as error:
        if n is not None and getattr(n, "model", None) is not None:
            report["solver_quality"] = accuracy.solver_quality(
                getattr(n.model, "solver_model", None)
            )
        report.update(
            status="needs_review"
            if (directory / "candidate_2035.nc").exists()
            else "failed",
            failed_phase=report["phase"],
            phase="failed",
            error=f"{type(error).__name__}: {error}",
        )
        raise
    finally:
        if (directory / "solver.log").is_file():
            report["solver_threads"] = accuracy.solver_thread_usage(
                directory / "solver.log"
            )
        checkpoint()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-target", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--threads", type=int, choices=(16, 32), default=16)
    parser.add_argument("--time-limit", type=int, default=252000)
    parser.add_argument(
        "--memory-profile", choices=tuple(MEMORY_PROFILES), default="standard"
    )
    parser.add_argument("--reference-summary", type=Path)
    args = parser.parse_args()
    if args.build_target:
        print(TARGET)
    else:
        if args.output_dir is None:
            parser.error("--output-dir is required")
        if args.memory_profile != "standard" and args.reference_summary is None:
            parser.error(f"--memory-profile {args.memory_profile} requires --reference-summary")
        logging.basicConfig(level=logging.INFO)
        solve(
            args.output_dir,
            args.threads,
            args.time_limit,
            memory_profile=args.memory_profile,
            reference_summary=args.reference_summary,
        )


if __name__ == "__main__":
    main()
