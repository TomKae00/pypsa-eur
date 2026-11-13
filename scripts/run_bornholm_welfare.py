# SPDX-License-Identifier: MIT
"""Prepare, solve and collect the 12 paired fixed-Bornholm welfare sensitivities."""

import argparse
import csv
import hashlib
import json
import logging
import os
import resource
import runpy
import socket
import time
from pathlib import Path

import numpy as np
import yaml

from scripts import bornholm_welfare_sensitivity as sensitivity

DEFAULT_CONFIG = Path("config/bornholm_welfare_sensitivity.yaml")


def digest(path):
    sha = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b""):
            sha.update(block)
    return sha.hexdigest()


def read_config(path=DEFAULT_CONFIG):
    config = yaml.safe_load(Path(path).read_text())
    sensitivity.require(
        config["version"] == 1 and config["planning_year"] == 2035,
        "Expected version 1, planning year 2035",
    )
    sensitivity.require(
        set(config["demand"]["factors"]) == {"central", "plus10", "plus20"},
        "Demand cases must be current, +10%, +20%",
    )
    sensitivity.require(
        list(config["demand"]["factors"].values()) == [1.0, 1.1, 1.2],
        "Unexpected demand factors",
    )
    sensitivity.require(
        config["demand"]["scope"] == "fixed_electricity_loads",
        "Unsupported demand scope",
    )
    return config


def paths(config, scenario):
    name = scenario["name"]
    return (
        Path("resources") / config["prefix"] / name / "networks/composed_2035.nc",
        Path("results") / config["prefix"] / name,
    )


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def code_fingerprint():
    from scripts.bornholm_provenance import code_provenance

    hashes = [
        code_provenance()["code_sha256"],
        digest("data/custom_extra_functionality.py"),
    ]
    return hashlib.sha256(json.dumps(hashes).encode()).hexdigest()


def plan(config):
    rows = sensitivity.scenarios(config)
    print(
        "index\tscenario\tdemand_factor\tDK1_onshore_MW\tDK2_onshore_MW\tDE_onshore_MW"
    )
    for row in rows:
        d = config["deployment"][row["deployment"]]
        print(
            f"{row['index']}\t{row['name']}\t{row['demand_factor']:.2f}\t"
            f"{d['onshore_DK1_MW']:.2f}\t{d['onshore_DK2_MW']:.2f}\t{d['onshore_DE_MW']:.2f}"
        )
    return rows


def check_pair(a, b):
    """The final-demand boundary must match within each project/no-project pair."""
    sensitivity.require(
        a.loads.index.equals(b.loads.index), "Paired load assets differ"
    )
    pa = sensitivity.demand_panel(a, a.loads.index, nonnegative=False)
    pb = sensitivity.demand_panel(b, b.loads.index, nonnegative=False)
    sensitivity.require(pa.equals(pb), "Paired final-demand time series differ")
    sensitivity.require(
        a.snapshot_weightings.equals(b.snapshot_weightings), "Paired weights differ"
    )
    ga = a.global_constraints.drop(columns="mu", errors="ignore")
    gb = b.global_constraints.drop(columns="mu", errors="ignore")
    sensitivity.require(ga.equals(gb), "Paired carbon/resource budgets differ")


def prepare(config, config_path):
    import pypsa

    from scripts import run_bornholm_calibrated as calibrated
    from scripts.bornholm_provenance import code_provenance

    sensitivity.require(
        config["source_prefix"] == calibrated.PREFIX, "Wrong source campaign"
    )
    provenance = code_provenance()
    code_sha = code_fingerprint()
    config_sha = digest(config_path)
    records = []
    pair_networks = []
    for scenario in sensitivity.scenarios(config):
        context = calibrated.context_for(scenario["source_run"])
        source = Path(context.input.network)
        sensitivity.require(
            source.is_file(),
            f"Missing composed source: {source}; run the existing calibrated build first",
        )
        source_sha = digest(source)
        target, directory = paths(config, scenario)
        marker = directory / "prepared.json"
        expected = dict(
            scenario=scenario,
            source=str(source),
            source_sha256=source_sha,
            campaign_config_sha256=config_sha,
            code_sha256=code_sha,
        )
        if target.exists() or marker.exists():
            sensitivity.require(
                target.exists() and marker.exists(),
                f"Incomplete prepared input: {target}",
            )
            record = json.loads(marker.read_text())
            sensitivity.require(
                all(record.get(k) == v for k, v in expected.items()),
                "Prepared campaign changed; use a new prefix rather than overwriting inputs",
            )
            sensitivity.require(
                digest(target) == record["prepared_sha256"],
                "Prepared input was modified",
            )
            n = pypsa.Network(target)
            sensitivity.validate(n, config, scenario)
        else:
            # Use the original native preflight before changing any demand or capacity constraints.
            source_preflight = calibrated.preflight(context)
            n = pypsa.Network(source)
            sensitivity.require(
                "solve_result" not in n.meta,
                "Use unsolved composed inputs, not a solved candidate",
            )
            sensitivity.require(
                not n.lines.s_nom_extendable.any(),
                "This campaign requires the fixed background grid",
            )
            if scenario["project"] == "hybrid":
                wind = n.generators.loc["BEI_WIND"]
                sensitivity.require(
                    wind.p_nom_extendable and wind.p_nom_min == wind.p_nom_max == 3000,
                    "Expected a fixed 3 GW Bornholm wind project",
                )
            else:
                sensitivity.require(
                    "BEI_WIND" not in n.generators.index,
                    "Baseline contains Bornholm wind",
                )
            audit = sensitivity.apply(n, config, scenario)
            settings = dict(
                enabled=True,
                scenario=scenario,
                capacity_limits=sensitivity.capacity_limits(config, scenario),
            )
            n.meta["bornholm_welfare_sensitivity"] = settings
            target.parent.mkdir(parents=True, exist_ok=True)
            temporary = target.with_name("composed_2035.tmp.nc")
            n.export_to_netcdf(temporary)
            temporary.replace(target)
            record = dict(
                expected,
                prepared_network=str(target),
                prepared_sha256=digest(target),
                source_preflight=source_preflight,
                assumptions=audit,
            )
            write_json(marker, record)
        pair_networks.append(n)
        if len(pair_networks) == 2:
            check_pair(*pair_networks)
            pair_networks.clear()
        records.append(record)
        print(f"Prepared {scenario['index']:02d}: {scenario['name']}", flush=True)
    manifest = dict(
        config=config,
        config_sha256=config_sha,
        code_sha256=code_sha,
        provenance=provenance,
        prepared=records,
    )
    write_json(Path("results") / config["prefix"] / "campaign.json", manifest)
    print("All 12 inputs and six paired demand/budget checks passed.")


def solve(
    config,
    config_path,
    index,
    *,
    scenario_table=None,
    context_factory=None,
    validate_fn=None,
):
    import pypsa

    from scripts import run_bornholm_accuracy as accuracy
    from scripts import run_bornholm_calibrated as calibrated
    from scripts import transmission_assumptions
    from scripts.bornholm_model_assumptions import (
        check_solved_assumptions,
        validate_assumptions,
    )
    from scripts.dispatch_checks import check as check_dispatch
    from scripts.run_bornholm_diagnostic import memory_snapshot, soft_memory_limit

    scenario = (
        sensitivity.scenarios(config) if scenario_table is None else scenario_table
    )[index]
    validate_network = sensitivity.validate if validate_fn is None else validate_fn
    source, directory = paths(config, scenario)
    prepared = json.loads((directory / "prepared.json").read_text())
    sensitivity.require(
        prepared["campaign_config_sha256"] == digest(config_path),
        "Campaign config changed after prepare",
    )
    sensitivity.require(
        prepared["code_sha256"] == code_fingerprint(), "Code changed after prepare"
    )
    sensitivity.require(
        prepared["prepared_sha256"] == digest(source), "Prepared network changed"
    )
    accepted = directory / "networks/solved_2035.nc"
    sensitivity.require(not accepted.exists(), f"Accepted output exists: {accepted}")
    attempt = (
        directory / "attempts" / (os.environ.get("SLURM_JOB_ID") or str(time.time_ns()))
    )
    attempt.mkdir(parents=True, exist_ok=False)
    context = (
        calibrated.context_for(scenario["source_run"])
        if context_factory is None
        else context_factory(scenario)
    )
    context.input = type(context.input)(fromdict={"network": str(source)})
    context.output = type(context.output)(
        fromdict={
            "network": str(attempt / "candidate_2035.nc"),
            "config": str(attempt / "effective_config.yaml"),
        }
    )
    context.log = type(context.log)(
        fromdict={k: str(attempt / f"{k}.log") for k in ("python", "solver", "memory")}
    )
    context.config["bornholm_welfare_sensitivity"] = dict(
        enabled=True,
        scenario=scenario,
        capacity_limits=sensitivity.capacity_limits(config, scenario),
    )
    context.config["run"]["prefix"] = config["prefix"]
    context.config["run"]["name"] = scenario["name"]
    context.threads = 16
    allocated = int(os.environ.get("SLURM_CPUS_PER_TASK", "0"))
    affinity = (
        len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else allocated
    )
    sensitivity.require(
        allocated >= 16 and affinity >= 16, "Allocate and bind 16 solver CPUs"
    )
    solving = context.params.solving
    options = solving["solver_options"][solving["solver"]["options"]]
    sensitivity.require(
        options["threads"] == 16
        and options["method"] == 2
        and options["crossover"] == 0,
        "Expected existing calibrated 16-thread barrier preset",
    )
    options["SoftMemLimit"] = soft_memory_limit(
        memory_snapshot(), options["SoftMemLimit"]
    )
    context.config["solving"] = solving
    started = time.monotonic()
    report = dict(
        status="running",
        phase="preflight",
        scenario=scenario,
        node=socket.gethostname(),
        source=str(source),
        source_sha256=prepared["prepared_sha256"],
        attempt=str(attempt),
        solver_options=dict(options),
        code_sha256=prepared["code_sha256"],
    )

    def checkpoint():
        report.update(
            elapsed_seconds=time.monotonic() - started,
            peak_rss_GiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2,
        )
        write_json(attempt / "summary.json", report)
        write_json(directory / "summary.json", report)

    try:
        n = pypsa.Network(source)
        validate_network(n, config, scenario)
        validate_assumptions(n, context.config["bornholm_model_assumptions"])
        transmission_assumptions.validate(n, context.config["transmission_assumptions"])
        calibrated.check_biomass_fuel(n)
        del n
        pypsa.options.params.optimize.include_objective_constant = True
        report["phase"] = "native_solve"
        checkpoint()
        namespace = runpy.run_path(
            "scripts/solve_network.py",
            run_name="__main__",
            init_globals={"snakemake": context},
        )
        n = namespace["n"]
        report.update(
            phase="quality_checks",
            solve_result=n.meta.get("solve_result", {}),
            candidate_network=str(attempt / "candidate_2035.nc"),
        )
        sensitivity.require(
            report["solve_result"]
            == {"status": "ok", "termination_condition": "optimal"},
            f"Non-optimal solve: {report['solve_result']}",
        )
        sensitivity.require(
            n.objective is not None and np.isfinite(n.objective), "Non-finite objective"
        )
        report["objective"] = float(n.objective)
        report["objective_constant"] = float(n.objective_constant)
        report["capacity_checks"] = validate_network(n, config, scenario, solved=True)
        report["assumption_checks"] = check_solved_assumptions(
            n, context.config["bornholm_model_assumptions"]
        )
        report["physical_quality"] = check_dispatch(n)
        report["solver_quality"] = accuracy.solver_quality(n.model.solver_model)
        report["barrier_final"] = accuracy.barrier_result(attempt / "solver.log")
        report["coal"] = accuracy.coal_check(
            n, fixed_MWh=context.config["bornholm_model_assumptions"]["coal_store_MWh"]
        )
        report["screening_issues"] = accuracy.accuracy_issues(
            report["solver_quality"],
            report["coal"],
            report["barrier_final"],
            crossover=0,
        )
        if report["screening_issues"]:
            report.update(status="optimal_needs_review", phase="complete")
            # Retain the optimal candidate and all diagnostics. Do not weaken the screens.
            return 3
        n.meta["bornholm_welfare_quality"] = dict(screening_passed=True)
        accepted.parent.mkdir(parents=True, exist_ok=True)
        n.export_to_netcdf(accepted)
        report.update(
            status="checks_passed", phase="complete", accepted_network=str(accepted)
        )
        return 0
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


def cost_cases(n, config):
    """Wind-only CAPEX sensitivity; fixed absolute O&M, collection and actual dispatch."""
    if "BEI_WIND" not in n.generators.index:
        return []
    g = n.generators.loc["BEI_WIND"]
    capacity = float(g.p_nom_opt)
    hours = n.snapshot_weightings.generators
    generation = float(n.generators_t.p.BEI_WIND @ hours)
    potential = float(n.generators_t.p_max_pu.BEI_WIND @ hours * capacity)
    operating = float(n.generators_t.p.BEI_WIND @ hours * g.marginal_cost)
    revenue = float(
        (n.generators_t.p.BEI_WIND * n.buses_t.marginal_price[g.bus])
        @ n.snapshot_weightings.objective
    )
    cfg = config["postprocessing"]
    native = cfg["wind_capex_EUR_kW"]["current"]
    rate = cfg["discount_rate"]
    life = cfg["wind_lifetime_years"]
    annuity = rate / (1 - (1 + rate) ** -life) if rate else 1 / life
    # capital_cost includes FOM. Only replace its annuitised wind CAPEX portion.
    rows = []
    for name, investment in cfg["wind_capex_EUR_kW"].items():
        annual = float(
            g.capital_cost * capacity
            + operating
            + (investment - native) * 1000 * capacity * annuity
        )
        rows.append(
            dict(
                cost_case=name,
                wind_capex_EUR_kW=investment,
                wind_capacity_MW=capacity,
                generation_MWh=generation,
                available_MWh=potential,
                market_revenue_EUR=revenue,
                annual_wind_collection_cost_EUR=annual,
                wind_collection_LCOE_EUR_MWh=annual / generation,
                cost_per_available_MWh=annual / potential,
                minimum_annual_net_support_EUR=annual - revenue,
                export_cost_included=False,
            )
        )
    return rows


def collect(config):
    import pypsa

    root = Path("results") / config["prefix"]
    rows, costs = [], []
    transfer = root / "transfer"
    transfer.mkdir(parents=True, exist_ok=True)
    for scenario in sensitivity.scenarios(config):
        _, directory = paths(config, scenario)
        summary = directory / "summary.json"
        report = (
            json.loads(summary.read_text())
            if summary.exists()
            else {"status": "missing"}
        )
        row = dict(
            scenario,
            status=report["status"],
            objective=report.get("objective"),
            elapsed_seconds=report.get("elapsed_seconds"),
            screening_issues="; ".join(report.get("screening_issues", [])),
            error=report.get("error", ""),
        )
        network = report.get("accepted_network") or report.get("candidate_network")
        row["network"] = network or ""
        target = transfer / scenario["name"]
        target.mkdir(exist_ok=True)
        link = target / "candidate_2035.nc"
        if link.is_symlink():
            link.unlink()
        else:
            sensitivity.require(
                not link.exists(),
                f"Refusing to replace a regular transfer file: {link}",
            )
        for name, source in (
            ("summary.json", summary),
            ("prepared.json", directory / "prepared.json"),
        ):
            if source.exists():
                (target / name).write_bytes(source.read_bytes())
        if network and Path(network).is_file():
            link.symlink_to(Path(network).resolve())
            n = pypsa.Network(network)
            if report.get("solve_result") == {
                "status": "ok",
                "termination_condition": "optimal",
            }:
                costs.extend(
                    dict(scenario=scenario["name"], status=report["status"], **r)
                    for r in cost_cases(n, config)
                )
            del n
            attempt = Path(report["attempt"])
            for name in ("effective_config.yaml", "solver.log"):
                if (attempt / name).is_file():
                    (target / name).write_bytes((attempt / name).read_bytes())
        rows.append(row)
    for name, data in (("campaign_status.csv", rows), ("wind_cost_cases.csv", costs)):
        if data:
            with (root / name).open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(data[0]))
                writer.writeheader()
                writer.writerows(data)
            (transfer / name).write_bytes((root / name).read_bytes())
    for name in ("campaign.json",):
        if (root / name).exists():
            (transfer / name).write_bytes((root / name).read_bytes())
    for row in rows:
        print(f"{row['index']:02d} {row['name']}: {row['status']}")
    print(f"Transfer directory (use rsync -L): {transfer}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument("--plan", action="store_true")
    actions.add_argument("--prepare", action="store_true")
    actions.add_argument("--index", type=int, choices=range(12))
    actions.add_argument("--collect", action="store_true")
    args = parser.parse_args()
    config = read_config(args.config)
    logging.basicConfig(level=logging.INFO)
    if args.plan:
        plan(config)
    elif args.prepare:
        prepare(config, args.config)
    elif args.collect:
        collect(config)
    else:
        raise SystemExit(solve(config, args.config, args.index))


if __name__ == "__main__":
    main()
