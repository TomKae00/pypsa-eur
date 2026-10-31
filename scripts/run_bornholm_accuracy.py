# SPDX-License-Identifier: MIT
"""Independent annual baseline replays for barrier accuracy; no rebuilding."""

import argparse
import hashlib
import json
import logging
import math
import os
import re
import resource
import runpy
import socket
import time
from pathlib import Path

from scripts import run_bornholm_study as study

RUN = "bornholm_zonal_baseline"
# Each variant changes one solver setting relative to nf1. No physical changes.
PROFILES = {
    "nf1": {
        "label": "tol8_nf1_t8",
        "BarConvTol": 1e-8,
        "ScaleFlag": -1,
        "crossover": 0,
        "BarHomogeneous": -1,
    },
    "nf1_scale2": {
        "label": "tol8_nf1_scale2_t8",
        "BarConvTol": 1e-8,
        "ScaleFlag": 2,
        "crossover": 0,
        "BarHomogeneous": -1,
    },
    "nf1_crossover": {
        "label": "tol8_nf1_crossover_t8",
        "BarConvTol": 1e-8,
        "ScaleFlag": -1,
        "crossover": -1,
        "BarHomogeneous": -1,
    },
    "nf1_homogeneous": {
        "label": "tol8_nf1_homogeneous_t8",
        "BarConvTol": 1e-8,
        "ScaleFlag": -1,
        "crossover": 0,
        "BarHomogeneous": 1,
    },
}
# Keep the running four-task sweep fixed when adding new trial profiles.
ORIGINAL_SWEEP = tuple(PROFILES)
FOLLOWUP_PROFILES = {
    "nf1_t16": PROFILES["nf1"] | {"label": "tol8_nf1_t16", "threads": 16},
    "nf1_homogeneous_t16": PROFILES["nf1_homogeneous"]
    | {"label": "tol8_nf1_homogeneous_t16", "threads": 16},
    "nf1_homogeneous_crossover": PROFILES["nf1_homogeneous"]
    | {"label": "tol8_nf1_homogeneous_crossover_t8", "crossover": -1},
    "nf0_homogeneous": PROFILES["nf1_homogeneous"]
    | {"label": "tol8_nf0_homogeneous_t8", "NumericFocus": 0},
}
PROFILES.update(FOLLOWUP_PROFILES)
# Matched dual-simplex crossover trials; keep both existing arrays unchanged.
CROSSOVER_PROFILES = {
    "nf1_crossover2_t8": PROFILES["nf1"]
    | {"label": "tol8_nf1_crossover2_t8", "crossover": 2, "threads": 8},
    "nf1_crossover2_t16": PROFILES["nf1"]
    | {"label": "tol8_nf1_crossover2_t16", "crossover": 2, "threads": 16},
}
PROFILES.update(CROSSOVER_PROFILES)
# Independent presolve trials, plus a matched automatic-thread PreDual trial.
PRESOLVE_PROFILES = {
    "nf1_predual_t16": PROFILES["nf1_t16"]
    | {"label": "tol8_nf1_predual_t16", "PreDual": -1, "AggFill": 0},
    "nf1_aggfill_t16": PROFILES["nf1_t16"]
    | {"label": "tol8_nf1_aggfill_t16", "PreDual": 0, "AggFill": -1},
    "nf1_predual_auto": PROFILES["nf1_t16"]
    | {
        "label": "tol8_nf1_predual_tauto",
        "PreDual": -1,
        "AggFill": 0,
        "threads": 0,
    },
}
PROFILES.update(PRESOLVE_PROFILES)
QUALITY = (
    "Status",
    "ObjVal",
    "ObjCon",
    "Runtime",
    "Work",
    "BarIterCount",
    "IterCount",
    "MemUsed",
    "MaxMemUsed",
    "MaxVio",
    "BoundVio",
    "BoundVioIndex",
    "BoundSVio",
    "ConstrVio",
    "ConstrVioIndex",
    "ConstrSVio",
    "ConstrResidual",
    "ConstrSResidual",
    "DualVio",
    "DualVioIndex",
    "DualSVio",
    "DualResidual",
    "DualResidualIndex",
    "DualSResidual",
    "ComplVio",
    "ComplVioSum",
    "NumVars",
    "NumConstrs",
    "IsMIP",
)
PINNED_FILES = (
    "scripts/solve_network.py",
    "scripts/_helpers.py",
    "config/config.default.yaml",
    "config/config.bornholm.yaml",
    "config/scenarios_bornholm.yaml",
    "data/custom_extra_functionality.py",
)
# The campaign's old source manifest omitted data/*.py; this is the separately
# audited no-op custom hook, not an assumption that the current hook is unchanged.
AUDITED_CUSTOM_HOOK_SHA256 = (
    "c61dcafdc27617f4e9a984368bc8a8737b64b4a65b6034dba876ed633243fc5c"
)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def configure_context(context, directory, numeric_focus, threads, *, profile=None):
    """Change solver options and output destinations, preserving source inputs."""
    study.require(numeric_focus in (0, 1, 2), "NumericFocus must be 0, 1 or 2")
    study.require(
        threads in (0, 8, 16),
        "Comparison runs require 0 (automatic), eight or 16 threads",
    )
    study.require(profile is None or profile in PROFILES, "Unknown solver profile")
    if profile is not None:
        spec = PROFILES[profile]
        study.require(
            numeric_focus == spec.get("NumericFocus", 1),
            "NumericFocus does not match profile",
        )
        study.require(
            threads == spec.get("threads", 8), "Thread count does not match profile"
        )
    workflow_threads = threads
    if threads == 0:
        # Gurobi's automatic setting normally uses up to 32 CPUs. Reserve those
        # CPUs explicitly; Snakemake still needs a positive resource count.
        allocated = os.environ.get("SLURM_CPUS_PER_TASK", "")
        study.require(
            allocated.isdigit() and int(allocated) >= 32,
            "Automatic Gurobi threads require SLURM_CPUS_PER_TASK >= 32",
        )
        workflow_threads = int(allocated)
    # Snakemake 9 reconstructs Params on each access. Mutate its stored dict:
    # assigning context.params.solving would change a temporary container only.
    solving = context.params.solving
    selected = solving["solver"]["options"]
    options = solving["solver_options"][selected]
    study.require(solving["solver"]["name"] == "gurobi", "Expected Gurobi")
    study.require(
        options["method"] == 2 and options["crossover"] == 0,
        "Expected the existing barrier / crossover-off preset",
    )
    options.update(threads=threads, BarConvTol=1e-8, NumericFocus=numeric_focus)
    if profile is not None:
        options.update(
            {key: value for key, value in PROFILES[profile].items() if key != "label"}
        )
    context.config["solving"] = solving
    context.threads = workflow_threads
    context.output = type(context.output)(
        fromdict={
            "network": str(directory / "candidate_2035.nc"),
            "config": str(directory / "effective_config.yaml"),
        }
    )
    context.log = type(context.log)(
        fromdict={
            "python": str(directory / "python.log"),
            "solver": str(directory / "solver.log"),
            "memory": str(directory / "memory.log"),
        }
    )
    return options


def solver_quality(model):
    values, unavailable = {}, {}
    if model is None:
        return {
            "values": {},
            "unavailable": {"model": "Native solver object unavailable"},
        }
    for name in QUALITY:
        try:
            value = model.getAttr(name)
            values[name] = (
                value if not isinstance(value, float) or math.isfinite(value) else None
            )
        except Exception as error:
            # Includes Gurobi's exception for attributes unavailable for an LP.
            unavailable[name] = f"{type(error).__name__}: {error}"
    try:
        import gurobipy

        values["gurobi_version"] = list(gurobipy.gurobi.version())
    except ImportError:
        pass
    return {"values": values, "unavailable": unavailable}


def barrier_result(path):
    last = None
    for line in Path(path).read_text().splitlines():
        fields = line.split()
        if len(fields) == 7 and re.fullmatch(r"\d+", fields[0]):
            try:
                last = dict(
                    zip(
                        (
                            "iteration",
                            "primal_eur",
                            "dual_eur",
                            "primal_residual",
                            "dual_residual",
                            "complementarity",
                            "seconds",
                        ),
                        [int(fields[0]), *[float(x.rstrip("s")) for x in fields[1:]]],
                    )
                )
            except ValueError:
                continue
    if last is not None:
        last["displayed_separation_eur"] = abs(last["primal_eur"] - last["dual_eur"])
        last["note"] = (
            "Rounded log iterates; dual feasibility is not certified. Not an error bound."
        )
    return last


def solver_thread_usage(path):
    """Retain reported solver limits, distinct from the requested Threads=0."""
    result = {}
    with Path(path).open() as stream:
        for line in stream:
            match = re.match(r"Thread count: .*using up to (\d+) threads", line)
            if match:
                result["reported_max_threads"] = int(match[1])
            match = re.match(r"\s*Threads\s*:\s*(\d+)\s*$", line)
            if match:
                result["barrier_threads"] = int(match[1])
    return result


def coal_check(n):
    name = "EU coal Store"
    study.require(name in n.stores.index, "Expected the audited coal Store")
    store = n.stores.loc[name]
    study.require(
        store.e_nom_extendable and store.e_cyclic and store.standing_loss == 0,
        "Coal storage formulation differs from the audited baseline",
    )
    study.require(
        store.e_nom_min == 0 and store.e_min_pu == 0 and store.e_max_pu == 1,
        "Coal storage bounds differ from the audited baseline",
    )
    for attribute, expected in (("e_min_pu", 0), ("e_max_pu", 1), ("standing_loss", 0)):
        series = getattr(n.stores_t, attribute)
        study.require(
            name not in series or series[name].eq(expected).all(),
            f"Time-varying coal {attribute} needs a separate audit",
        )
    energy = n.stores_t.e[name]
    capacity = float(store.e_nom_opt)
    cost = float(store.capital_cost)
    return {
        "capacity_MWh": capacity,
        "capital_cost_eur_per_MWh_year": cost,
        "inventory_min_MWh": float(energy.min()),
        "inventory_max_MWh": float(energy.max()),
        "inventory_range_MWh": float(energy.max() - energy.min()),
        "capacity_only_avoidable_cost_eur": max(0.0, capacity - float(energy.max()))
        * cost,
        "same_dispatch_avoidable_cost_eur": max(
            0.0, capacity - float(energy.max() - energy.min())
        )
        * cost,
    }


def accuracy_issues(quality, coal, barrier, *, crossover=0):
    """Diagnostic screening, not an independent primal/dual optimality certificate."""
    issues = []
    avoidable = coal["same_dispatch_avoidable_cost_eur"]
    if avoidable is None or not math.isfinite(avoidable) or avoidable > 1e6:
        issues.append(
            "Avoidable coal storage cost missing, non-finite or above EUR 1m/year"
        )
    values = quality["values"]
    if values.get("Status") != 2:
        issues.append("Native Gurobi status is not OPTIMAL (2)")
    for name, threshold in (
        ("BoundVio", 1e-3),
        ("ConstrVio", 1e-3),
        ("ConstrResidual", 1e-3),
        ("DualVio", 1e-5),
        ("DualResidual", 1e-5),
    ):
        value = values.get(name)
        if value is None or not math.isfinite(value) or value > threshold:
            issues.append(
                f"{name} missing or exceeds screening threshold {threshold:g}"
            )
    if crossover == 0:
        # The rounded log row is only an extra diagnostic screen. It need not
        # describe the returned candidate and is never a certified bound.
        separation = None if barrier is None else barrier["displayed_separation_eur"]
        if separation is None or not math.isfinite(separation) or separation > 5e6:
            issues.append(
                "Final displayed primal/dual separation missing or exceeds EUR 5m"
            )
    else:
        # Crossover can clean up a poor barrier iterate. Judge the returned LP
        # solution, not the earlier log row; retain feasibility/coal screens.
        # Complementarity sums are in objective units for this unscaled EUR LP.
        # This is a diagnostic threshold, not a certified objective-error bound.
        value = values.get("ComplVioSum")
        if value is None or not math.isfinite(value) or value > 5e6:
            issues.append(
                "Final ComplVioSum missing, non-finite or exceeds EUR 5m after crossover"
            )
    return issues


def solve(numeric_focus, threads, directory, *, profile=None):
    from scripts.bornholm_provenance import code_provenance
    from scripts.run_bornholm_diagnostic import memory_snapshot, soft_memory_limit

    directory = Path(directory).resolve()
    root = (Path("results") / "bornholm_accuracy").resolve()
    study.require(
        directory.is_relative_to(root) and directory != root,
        "Output must be a new subdirectory of results/bornholm_accuracy",
    )
    directory.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    report = {
        "status": "running",
        "phase": "initialising",
        "scenario": RUN,
        "numeric_focus": numeric_focus,
        "profile": profile,
        "threads": threads,
        "slurm_cpus_per_task": os.environ.get("SLURM_CPUS_PER_TASK"),
        "visible_cpu_count": os.cpu_count(),
        "cpu_affinity_count": (
            len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None
        ),
        "node": socket.gethostname(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "slurm_array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),
        "quality_screening_note": "Diagnostic thresholds only; no independent optimality certificate.",
    }

    def checkpoint():
        report["elapsed_seconds"] = time.monotonic() - started
        report["peak_rss_GiB"] = (
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2
        )
        temporary = directory / "summary.json.tmp"
        temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        temporary.replace(directory / "summary.json")

    try:
        if profile in CROSSOVER_PROFILES or profile in PRESOLVE_PROFILES:
            import pypsa

            # Retain the audited accounting convention even if an environment
            # or module-wide option overrides PyPSA's current default.
            pypsa.options.params.optimize.include_objective_constant = True
            report["include_objective_constant"] = True
        checkpoint()
        context = study.context_for(RUN)
        source = Path(context.input.network)
        reference_path = Path("results") / study.PREFIX / RUN / "summary.json"
        reference = json.loads(reference_path.read_text())
        study.require(
            reference.get("status") == "complete",
            "Missing completed reference baseline",
        )
        source_hash = sha256(source)
        study.require(
            source_hash == reference["source_sha256"],
            "Composed input changed since the completed baseline; do not rebuild for this comparison",
        )
        hashes = reference["code_files_sha256"]
        for path in PINNED_FILES:
            expected = hashes.get(path) or (
                reference.get("custom_constraints_sha256")
                if path.endswith("custom_extra_functionality.py")
                else None
            )
            if path.endswith("custom_extra_functionality.py") and expected is None:
                expected = AUDITED_CUSTOM_HOOK_SHA256
            study.require(
                expected is not None and sha256(path) == expected,
                f"Reference physical formulation changed: {path}",
            )
        report.update(
            source=str(source),
            source_sha256=source_hash,
            reference_summary=str(reference_path),
            reference_objective=reference["objective"],
            **code_provenance(),
        )
        report["phase"] = "preflight"
        checkpoint()
        report["preflight"] = study.preflight(context)
        options = configure_context(
            context, directory, numeric_focus, threads, profile=profile
        )
        memory = memory_snapshot()
        options["SoftMemLimit"] = soft_memory_limit(memory, options["SoftMemLimit"])
        report.update(
            solver_options=dict(options),
            workflow_threads=context.threads,
            node_memory_before_solve_GB=memory,
            screening_basis=(
                "returned_solution_quality_crossover_enabled"
                if options["crossover"] != 0
                else "returned_solution_quality_and_barrier_log_screen"
            ),
            phase="native_solve",
        )
        checkpoint()
        print(
            f"Baseline accuracy comparison: profile={profile or 'legacy'}, "
            f"BarConvTol={options['BarConvTol']:g}, NumericFocus={numeric_focus}, "
            f"ScaleFlag={options.get('ScaleFlag', -1)}, crossover={options['crossover']}, "
            f"BarHomogeneous={options.get('BarHomogeneous', -1)}, "
            f"PreDual={options.get('PreDual', -1)}, AggFill={options.get('AggFill', -1)}, "
            f"threads={threads}, allocated_cpus={report['slurm_cpus_per_task']}",
            flush=True,
        )
        print(f"Source SHA256: {source_hash}\nResults: {directory}", flush=True)
        namespace = runpy.run_path(
            "scripts/solve_network.py",
            run_name="__main__",
            init_globals={"snakemake": context},
        )
        n = namespace["n"]
        report["phase"] = "quality_checks"
        report["solve_result"] = n.meta.get("solve_result", {})
        objective = float(n.objective)
        study.require(math.isfinite(objective), "Non-finite objective")
        report["objective"] = objective
        report["solver_quality"] = solver_quality(
            getattr(n.model, "solver_model", None)
        )
        report["barrier_final"] = barrier_result(directory / "solver.log")
        report["coal"] = coal_check(n)
        report["screening_issues"] = accuracy_issues(
            report["solver_quality"],
            report["coal"],
            report["barrier_final"],
            crossover=options["crossover"],
        )
        checkpoint()
        study.require(
            report["solve_result"]
            == {"status": "ok", "termination_condition": "optimal"},
            f"Solver did not terminate optimally: {report['solve_result']}",
        )
        study.require(
            not report["screening_issues"],
            f"Numerical result needs review: {report['screening_issues']}; candidate retained",
        )
        (directory / "candidate_2035.nc").replace(directory / "solved_2035.nc")
        report.update(status="checks_passed", phase="complete")
    except Exception as error:
        # Native solve_network.py can raise before exporting a non-optimal result.
        # Preserve native solver diagnostics without weakening its rejection or
        # altering the pinned physical solve script.
        native_path = Path("scripts/solve_network.py").resolve()
        trace = error.__traceback__
        while trace is not None:
            if Path(trace.tb_frame.f_code.co_filename).resolve() == native_path:
                network = trace.tb_frame.f_globals.get("n")
                if network is not None:
                    try:
                        model = getattr(network, "model", None)
                        report["solver_quality"] = solver_quality(
                            getattr(model, "solver_model", None)
                        )
                        report["native_model_status"] = str(
                            getattr(model, "status", "unavailable")
                        )
                        report["native_termination_condition"] = str(
                            getattr(model, "termination_condition", "unavailable")
                        )
                    except Exception as diagnostic_error:
                        report["failure_diagnostics_error"] = (
                            f"{type(diagnostic_error).__name__}: {diagnostic_error}"
                        )
                    break
            trace = trace.tb_next
        if "barrier_final" not in report and (directory / "solver.log").is_file():
            report["barrier_final"] = barrier_result(directory / "solver.log")
        report.update(
            status="failed",
            failed_phase=report["phase"],
            phase="failed",
            error=f"{type(error).__name__}: {error}",
        )
        raise
    finally:
        if (directory / "solver.log").is_file():
            report["solver_threads"] = solver_thread_usage(directory / "solver.log")
        checkpoint()


def compare(directories):
    if isinstance(directories, (str, Path)):
        directories = [directories]
    roots = list(dict.fromkeys(Path(directory).resolve() for directory in directories))
    paths = sorted(path for root in roots for path in root.glob("*/summary.json"))
    study.require(paths, "No run summaries found under comparison directory")
    labels = [
        f"{path.parent.parent.name}/{path.parent.name}"
        if len(roots) > 1
        else path.parent.name
        for path in paths
    ]
    study.require(len(labels) == len(set(labels)), "Duplicate run labels in comparison")
    reports = [json.loads(path.read_text()) for path in paths]
    hashes = [r.get("source_sha256") for r in reports]
    study.require(
        all(hashes) and len(set(hashes)) == 1,
        "Comparison inputs differ or lack source hashes",
    )
    options = [
        {
            k: v
            for k, v in r.get("solver_options", {}).items()
            if k
            not in (
                "NumericFocus",
                "SoftMemLimit",
                "ScaleFlag",
                "crossover",
                "BarConvTol",
                "BarHomogeneous",
                "threads",
                "PreDual",
                "AggFill",
            )
        }
        for r in reports
    ]
    study.require(
        all(option == options[0] for option in options),
        "Other solver settings differ between the runs",
    )
    columns = []
    for report in reports:
        quality = report.get("solver_quality", {}).get("values", {})
        columns.append(
            {
                "status": report["status"],
                "profile": report.get("profile"),
                "threads": report.get(
                    "threads", report.get("solver_options", {}).get("threads")
                ),
                "node": report.get("node"),
                "slurm_cpus_per_task": report.get("slurm_cpus_per_task"),
                "cpu_affinity_count": report.get("cpu_affinity_count"),
                "solver_threads": report.get("solver_threads"),
                "solver_options": report.get("solver_options", {}),
                "error": report.get("error"),
                "hours": round(report.get("elapsed_seconds", 0) / 3600, 3),
                "peak_RSS_GiB": report.get("peak_rss_GiB"),
                "Gurobi_MaxMemUsed_GB": quality.get("MaxMemUsed"),
                "barrier_iterations": quality.get("BarIterCount"),
                "objective_EUR": report.get("objective", quality.get("ObjVal")),
                "Gurobi_status": quality.get("Status"),
                "native_termination_condition": report.get(
                    "native_termination_condition"
                ),
                "coal_avoidable_EUR_per_year": report.get("coal", {}).get(
                    "same_dispatch_avoidable_cost_eur"
                ),
                "displayed_separation_EUR": report.get("barrier_final", {}).get(
                    "displayed_separation_eur"
                )
                if report.get("barrier_final")
                else None,
                "MaxVio": quality.get("MaxVio"),
                "DualVio": quality.get("DualVio"),
                "DualResidual": quality.get("DualResidual"),
                "ComplVio": quality.get("ComplVio"),
                "ComplVioSum": quality.get("ComplVioSum"),
                "screening_basis": report.get("screening_basis"),
                "screening_issues": report.get("screening_issues"),
            }
        )
    print(
        json.dumps(
            {
                "same_source_sha256": hashes[0],
                "runs": {label: column for label, column in zip(labels, columns)},
                "note": "Screening only; logged dual values are not certified bounds. Check full solver options and final quality; runtime depends on node hardware. A failed run's ObjVal is an unaccepted candidate.",
            },
            indent=2,
            allow_nan=False,
        )
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--numeric-focus", type=int, choices=(0, 1, 2))
    action.add_argument("--profile", choices=tuple(PROFILES), help="Named solver trial")
    action.add_argument(
        "--compare",
        type=Path,
        nargs="+",
        help="Read summaries in one or more array output directories",
    )
    parser.add_argument(
        "--threads",
        type=int,
        choices=(0, 8, 16),
        help="Default: profile thread count, or 8; 0 is automatic and needs 32 allocated CPUs",
    )
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    if args.compare:
        compare(args.compare)
    else:
        if args.output_dir is None:
            parser.error("--output-dir is required for a solve")
        spec = PROFILES[args.profile] if args.profile else {}
        focus = spec.get("NumericFocus", 1) if args.profile else args.numeric_focus
        threads = args.threads if args.threads is not None else spec.get("threads", 8)
        solve(focus, threads, args.output_dir, profile=args.profile)


if __name__ == "__main__":
    main()
