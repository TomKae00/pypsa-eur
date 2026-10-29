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
QUALITY = (
    "Status",
    "ObjVal",
    "ObjCon",
    "Runtime",
    "Work",
    "BarIterCount",
    "MemUsed",
    "MaxMemUsed",
    "MaxVio",
    "BoundVio",
    "BoundSVio",
    "ConstrVio",
    "ConstrSVio",
    "ConstrResidual",
    "ConstrSResidual",
    "DualVio",
    "DualSVio",
    "DualResidual",
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


def configure_context(context, directory, numeric_focus, threads):
    """Change solver options and output destinations, preserving source inputs."""
    study.require(numeric_focus in (0, 2), "NumericFocus must be 0 or 2")
    study.require(threads == 8, "Both comparison runs must use eight threads")
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
    context.config["solving"] = solving
    context.threads = threads
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


def accuracy_issues(quality, coal, barrier):
    """Diagnostic screening, not an independent primal/dual optimality certificate."""
    issues = []
    if coal["same_dispatch_avoidable_cost_eur"] > 1e6:
        issues.append("More than EUR 1m/year of avoidable coal storage cost")
    values = quality["values"]
    for name, threshold in (
        ("BoundVio", 1e-3),
        ("ConstrVio", 1e-3),
        ("ConstrResidual", 1e-3),
        ("DualVio", 1e-5),
        ("DualResidual", 1e-5),
    ):
        value = values.get(name)
        if value is None or value > threshold:
            issues.append(
                f"{name} missing or exceeds screening threshold {threshold:g}"
            )
    if barrier is None or barrier["displayed_separation_eur"] > 5e6:
        issues.append(
            "Final displayed primal/dual separation missing or exceeds EUR 5m"
        )
    return issues


def solve(numeric_focus, threads, directory):
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
        "threads": threads,
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
        options = configure_context(context, directory, numeric_focus, threads)
        memory = memory_snapshot()
        options["SoftMemLimit"] = soft_memory_limit(memory, options["SoftMemLimit"])
        report.update(
            solver_options=dict(options),
            node_memory_before_solve_GB=memory,
            phase="native_solve",
        )
        checkpoint()
        print(
            f"Baseline accuracy comparison: BarConvTol=1e-8, NumericFocus={numeric_focus}, threads={threads}",
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
            report["solver_quality"], report["coal"], report["barrier_final"]
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
        report.update(
            status="failed",
            failed_phase=report["phase"],
            phase="failed",
            error=f"{type(error).__name__}: {error}",
        )
        raise
    finally:
        checkpoint()


def compare(directory):
    reports = [
        json.loads((directory / f"tol8_nf{focus}_t8/summary.json").read_text())
        for focus in (0, 2)
    ]
    hashes = [r.get("source_sha256") for r in reports]
    study.require(
        all(hashes) and hashes[0] == hashes[1],
        "Comparison inputs differ or lack source hashes",
    )
    options = [
        {
            k: v
            for k, v in r.get("solver_options", {}).items()
            if k not in ("NumericFocus", "SoftMemLimit")
        }
        for r in reports
    ]
    study.require(
        options[0] == options[1], "Other solver settings differ between the two runs"
    )
    columns = []
    for report in reports:
        quality = report.get("solver_quality", {}).get("values", {})
        columns.append(
            {
                "status": report["status"],
                "hours": round(report.get("elapsed_seconds", 0) / 3600, 3),
                "peak_RSS_GiB": report.get("peak_rss_GiB"),
                "Gurobi_MaxMemUsed_GB": quality.get("MaxMemUsed"),
                "barrier_iterations": quality.get("BarIterCount"),
                "objective_EUR": report.get("objective"),
                "coal_avoidable_EUR_per_year": report.get("coal", {}).get(
                    "same_dispatch_avoidable_cost_eur"
                ),
                "displayed_separation_EUR": report.get("barrier_final", {}).get(
                    "displayed_separation_eur"
                )
                if report.get("barrier_final")
                else None,
                "MaxVio": quality.get("MaxVio"),
                "DualResidual": quality.get("DualResidual"),
                "screening_issues": report.get("screening_issues"),
            }
        )
    print(
        json.dumps(
            {
                "same_source_sha256": hashes[0],
                "NumericFocus_0": columns[0],
                "NumericFocus_2": columns[1],
                "note": "Screening only; logged dual values are not certified bounds. Thread speedup versus prior four-thread run is confounded by tighter convergence.",
            },
            indent=2,
            allow_nan=False,
        )
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--numeric-focus", type=int, choices=(0, 2))
    action.add_argument(
        "--compare", type=Path, help="Read summaries in an array output directory"
    )
    parser.add_argument("--threads", type=int, default=8, choices=(8,))
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    if args.compare:
        compare(args.compare)
    else:
        if args.output_dir is None:
            parser.error("--output-dir is required for a solve")
        solve(args.numeric_focus, args.threads, args.output_dir)


if __name__ == "__main__":
    main()
