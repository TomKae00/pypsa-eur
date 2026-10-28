# SPDX-License-Identifier: MIT
"""Build targets, input preflight and native annual solves for the eight 3h cases."""

import argparse
import hashlib
import json
import logging
import resource
import runpy
import socket
import time
from pathlib import Path

import yaml

CONFIGS = ["config/config.bornholm.yaml", "config/config.bornholm-3h.yaml"]
PREFIX = "bornholm_3h_s1_w2_o3"
CLASSES = {
    "solar": 1,
    "solar-hsat": 1,
    "onwind": 2,
    "offwind-ac": 3,
    "offwind-dc": 3,
    "offwind-float": 3,
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def study_runs():
    definitions = yaml.safe_load(Path("config/scenarios_bornholm.yaml").read_text())
    require(len(definitions) == 8, "Expected exactly eight study scenarios")
    return list(definitions)


def context_for(run):
    from scripts._helpers import mock_snakemake, set_scenario_config

    context = mock_snakemake(
        "solve_network", configfiles=CONFIGS, run=run, horizon="2035"
    )
    set_scenario_config(context)
    return context


def preflight(context):
    import pandas as pd
    import pypsa

    from scripts._helpers import mock_snakemake
    from scripts.run_bornholm_diagnostic import check_source

    source_path = Path(context.input.network)
    n = pypsa.Network(source_path)
    for technology, classes in CLASSES.items():
        require(
            n.meta.get("renewable", {}).get(technology, {}).get("resource_classes")
            == classes,
            f"Stale {technology} resource classes; rebuild {source_path}",
        )
    require(n.meta.get("run", {}).get("prefix") == PREFIX, "Wrong study namespace")
    # Gas inputs are shared under the native base policy; use the actual rule
    # dependency rather than assuming a scenario-local CSV.
    compose = mock_snakemake(
        "compose_network",
        configfiles=CONFIGS,
        run=context.wildcards.run,
        horizon="2035",
    )
    candidates = [
        Path(p)
        for p in compose.input
        if str(p).endswith("gas_input_locations_simplified.csv")
    ]
    require(len(candidates) == 1, "Expected one shared gas-input dependency")
    inputs = candidates[0]
    report = check_source(n, pd.read_csv(inputs, index_col=0), step_hours=3)
    report.update(
        snapshots=len(n.snapshots),
        annual_weight_hours=float(n.snapshot_weightings.objective.sum()),
        components={c.list_name: len(c.static) for c in n.components},
    )
    return report


def solve(run):
    import numpy as np

    from scripts.bornholm_provenance import code_provenance
    from scripts.run_bornholm_diagnostic import memory_snapshot, soft_memory_limit

    context = context_for(run)
    source = Path(context.input.network)
    output = Path(context.output.network)
    require(
        not output.exists(),
        f"Solved output already exists: {output}; archive it before resubmitting",
    )
    directory = Path("results") / PREFIX / run
    directory.mkdir(parents=True, exist_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    for path in context.log:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    report = {
        "status": "running",
        "phase": "preflight",
        "scenario": run,
        "node": socket.gethostname(),
        "source": str(source),
        **code_provenance(),
    }

    def checkpoint():
        report["elapsed_seconds"] = time.monotonic() - started
        report["peak_rss_GiB"] = (
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2
        )
        temporary = directory / "summary.json.tmp"
        temporary.write_text(json.dumps(report, indent=2) + "\n")
        temporary.replace(directory / "summary.json")

    try:
        checkpoint()
        report["preflight"] = preflight(context)
        report["source_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
        solving = context.params.solving
        options = solving["solver_options"][solving["solver"]["options"]]
        snapshot = memory_snapshot()
        options["SoftMemLimit"] = soft_memory_limit(snapshot, options["SoftMemLimit"])
        require(
            options["method"] == 2 and options["threads"] == 4,
            "Expected four-thread barrier",
        )
        report.update(
            solver_options=dict(options),
            node_memory_before_solve_GB=snapshot,
            phase="native_solve",
        )
        checkpoint()
        logging.info(
            "Scenario %s; native annual 3h solve; memory cap %.2f GB",
            run,
            options["SoftMemLimit"],
        )
        # Execute the native entry point and respect its configured solve mode,
        # including all supplementary constraints and logging.
        namespace = runpy.run_path(
            "scripts/solve_network.py",
            run_name="__main__",
            init_globals={"snakemake": context},
        )
        n = namespace["n"]
        result = n.meta.get("solve_result", {})
        report.update(
            solver_status=result.get("status"),
            termination_condition=result.get("termination_condition"),
        )
        require(
            result == {"status": "ok", "termination_condition": "optimal"},
            f"Run was not optimal: {result}",
        )
        require(np.isfinite(n.objective), "Non-finite objective")
        report.update(status="complete", phase="complete", objective=float(n.objective))
    except Exception as error:
        report.update(
            status="failed",
            failed_phase=report["phase"],
            phase="failed",
            error=f"{type(error).__name__}: {error}",
        )
        # Native PyPSA-Eur can export after an error status. Do not leave an
        # unaccepted file at the canonical solved target.
        if output.exists():
            output.replace(output.with_name("failed_2035.nc"))
        raise
    finally:
        checkpoint()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--list-build-targets", action="store_true")
    action.add_argument("--preflight-all", action="store_true")
    action.add_argument("--index", type=int, choices=range(8))
    args = parser.parse_args()
    runs = study_runs()
    if args.list_build_targets:
        for run in runs:
            print(f"resources/{PREFIX}/{run}/networks/composed_2035.nc")
    elif args.preflight_all:
        for run in runs:
            print(run, json.dumps(preflight(context_for(run))))
    else:
        logging.basicConfig(level=logging.INFO)
        solve(runs[args.index])


if __name__ == "__main__":
    main()
