# SPDX-License-Identifier: MIT
"""
Short chronological Bornholm diagnostic using the native solve constraints.

This is a feasibility/solver test, not an annual investment result. Annual flow
budgets, annual HVC feedstock and capital costs are prorated to the window;
physical power/energy capacities and one-hour storage chronology are retained.
"""

import argparse
import copy
import hashlib
import json
import logging
import resource
import subprocess
import time
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import pypsa

from scripts._helpers import mock_snakemake, set_scenario_config
from scripts.prepare_sector_network import gas_storage_minimum_capacities
from scripts.solve_network import (
    collect_kwargs,
    create_optimization_model,
    prepare_network,
)

logger = logging.getLogger(__name__)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def annual_gas_supply_check(n):
    """
    Necessary, optimistic annual flow check for fixed direct/industrial gas.

    Ignores competition for biogas, electricity, hydrogen and captured CO2, and
    variable gas demand for heat/power. Passing cannot establish feasibility.
    Unlimited conversion candidates are deliberately treated as unlimited supply.
    """
    gas_buses = n.buses.index[n.buses.carrier.eq("gas")]
    demand = pd.Series(0.0, index=gas_buses)
    loads = (
        n.get_switchable_as_dense("Load", "p_set")
        .mul(n.snapshot_weightings.generators, axis=0)
        .sum()
    )
    for name, load in n.loads.iterrows():
        if load.bus in demand.index:
            demand[load.bus] += max(0, loads[name])
        elif n.buses.at[load.bus, "carrier"] == "gas for industry":
            routes = n.links[n.links.bus0.isin(gas_buses) & n.links.bus1.eq(load.bus)]
            require(len(routes) > 0, f"No gas supply route to {load.bus}")
            require(routes.bus0.nunique() == 1, f"Ambiguous gas supply to {load.bus}")
            demand[routes.bus0.iloc[0]] += max(0, loads[name]) / routes.efficiency.max()
    total = float(demand.sum())
    graph = nx.DiGraph()
    graph.add_nodes_from(["source", "sink", *gas_buses])
    supply = pd.Series(0.0, index=gas_buses)
    hours = float(n.snapshot_weightings.generators.sum())
    for _, gen in n.generators[n.generators.bus.isin(gas_buses)].iterrows():
        capacity = gen.p_nom_max if gen.p_nom_extendable else gen.p_nom
        supply[gen.bus] += max(0, capacity * hours)
    for _, link in n.links.iterrows():
        if link.bus1 not in gas_buses or link.efficiency <= 0:
            continue
        capacity = link.p_nom_max if link.p_nom_extendable else link.p_nom
        if link.bus0 in gas_buses:
            # An optimistic capacity bound; native hourly constraints still apply.
            amount = min(total, capacity * max(0, link.p_max_pu) * hours)
            previous = graph.get_edge_data(link.bus0, link.bus1, {}).get("capacity", 0)
            graph.add_edge(link.bus0, link.bus1, capacity=previous + amount)
            if link.p_min_pu < 0:
                previous = graph.get_edge_data(link.bus1, link.bus0, {}).get(
                    "capacity", 0
                )
                graph.add_edge(
                    link.bus1,
                    link.bus0,
                    capacity=previous + min(total, capacity * -link.p_min_pu * hours),
                )
        else:
            # Relax all upstream constraints, including energy limits on biogas.
            supply[link.bus1] += min(
                total, capacity * max(0, link.p_max_pu) * link.efficiency * hours
            )
    for bus in gas_buses:
        graph.add_edge("source", bus, capacity=min(total, float(supply[bus])))
        graph.add_edge(bus, "sink", capacity=float(demand[bus]))
    flow = nx.maximum_flow_value(graph, "source", "sink")
    shortage = max(0.0, total - flow)
    require(
        shortage <= max(1e-3, total * 1e-8),
        f"Annual fixed-gas supply shortfall: {shortage / 1e6:.6f} TWh",
    )
    return {
        "fixed_gas_demand_TWh": total / 1e6,
        "optimistic_supply_shortfall_TWh": shortage / 1e6,
        "scope": "necessary condition only; upstream conversion inputs and variable gas demands relaxed",
    }


def check_source(n, gas_inputs):
    require(
        len(n.snapshots) == 8760, "Diagnostic requires the rebuilt full-year source"
    )
    require(np.allclose(n.snapshot_weightings, 1), "Source must have hourly weights")
    require(
        n.snapshots.equals(pd.date_range(n.snapshots[0], periods=8760, freq="h")),
        "Source must have consecutive hourly snapshots",
    )
    require(
        n.buses.carrier.eq("AC").sum() == 52, "Expected 52 electricity bidding zones"
    )
    for key, expected in {
        "gas_network": True,
        "gas_network_connectivity_upgrade": 1,
        "H2_network": False,
        "H2_retrofit": False,
        "co2_spatial": True,
        "co2_network": True,
    }.items():
        require(
            n.meta.get("sector", {}).get(key) == expected,
            f"Stale source setting sector.{key}; rebuild baseline",
        )
    gas = n.links[n.links.carrier.eq("gas pipeline")]
    require(
        not gas.empty and (gas.p_nom > 0).all(),
        "Missing/zero-capacity existing gas pipelines; rebuild gas inputs",
    )
    require(not gas.p_nom_extendable.any(), "Existing gas pipelines must remain fixed")
    for component in n.components:
        table = component.static
        for column in table.columns[table.columns.str.match(r"^bus\d*$")]:
            refs = table[column]
            # Empty optional ports are native multi-port Link placeholders.
            bad = refs.notna() & refs.ne("") & ~refs.isin(n.buses.index)
            require(
                not bad.any(),
                f"Undefined buses in {component.name}.{column}: {refs[bad].to_dict()}",
            )
        for minimum in ["p_nom_min", "e_nom_min", "s_nom_min"]:
            maximum = minimum.replace("_min", "_max")
            if minimum in table and maximum in table:
                require(
                    (table[minimum] <= table[maximum]).all(),
                    f"Conflicting bounds in {component.name}",
                )
    gas_stores = n.stores.index[n.stores.carrier.eq("gas")]
    expected = gas_storage_minimum_capacities(gas_inputs, gas_stores)
    require(
        np.allclose(n.stores.loc[gas_stores, "e_nom_min"], expected),
        "Stale gas-storage bounds; rebuild composed baseline",
    )
    return {
        "gas_storage_minimum_TWh": float(expected.sum() / 1e6),
        "annual_gas_check": annual_gas_supply_check(n),
    }


def short_network(source, hours, start=0):
    require(0 < hours < len(source.snapshots), "Hours must be between 1 and 8759")
    require(
        0 <= start <= len(source.snapshots) - hours,
        "Diagnostic window lies outside source year",
    )
    n = source.copy(snapshots=source.snapshots[start : start + hours])
    n.meta = copy.deepcopy(source.meta)
    fraction = (
        n.snapshot_weightings.generators.sum()
        / source.snapshot_weightings.generators.sum()
    )
    # Physical storage steps remain one hour, rather than representing 8760/hours.
    require(
        np.allclose(n.snapshot_weightings, 1), "Diagnostic lost physical hourly weights"
    )
    annual = n.global_constraints.type.isin(
        ["operational_limit", "primary_energy", "co2_atmosphere"]
    )
    known = annual | n.global_constraints.type.isin(
        [
            "transmission_volume_expansion_limit",
            "transmission_expansion_cost_limit",
            "tech_capacity_expansion_limit",
        ]
    )
    require(known.all(), "Unrecognised global constraint: review diagnostic scaling")
    n.global_constraints.loc[annual, "constant"] *= fraction
    for component in n.components:
        table = component.static
        if "capital_cost" in table:
            table["capital_cost"] *= fraction
        for column in ["e_sum_min", "e_sum_max"]:
            if column in table:
                finite = np.isfinite(table[column])
                table.loc[finite, column] *= fraction
    # HVC is an annual waste feedstock inventory, not a physical energy-storage
    # capacity. Preserve native mandatory depletion at the diagnostic endpoint.
    hvc = n.stores.index[n.stores.carrier.eq("non-sequestered HVC")]
    n.stores.loc[hvc, ["e_nom", "e_initial"]] *= fraction
    if len(hvc):
        n.stores_t.e_max_pu.loc[n.snapshots[-1], hvc] = 0.0
    n.meta["diagnostic"] = {
        "hours": hours,
        "start_offset": start,
        "annual_fraction": float(fraction),
        "interpretation": "short-window feasibility and solver diagnostic; not an annual scenario result",
    }
    return n


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hours", type=int, default=24)
    parser.add_argument(
        "--start", type=int, default=0, help="Hourly offset into source year"
    )
    parser.add_argument("--method", choices=["barrier", "pdhg"], default="barrier")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--time-limit", type=int, default=3600)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    require(
        args.threads > 0 and args.time_limit > 0,
        "Positive threads and time limit required",
    )
    args.output_dir.mkdir(parents=True, exist_ok=False)
    logging.basicConfig(
        level=logging.INFO,
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(args.output_dir / "python.log"),
        ],
    )
    started = time.monotonic()
    report = {
        "arguments": {**vars(args), "output_dir": str(args.output_dir)},
        "status": "failed",
    }
    try:
        context = mock_snakemake(
            "solve_network",
            configfiles="config/config.bornholm.yaml",
            run="bornholm_zonal_baseline",
            horizon="2035",
        )
        set_scenario_config(context)
        config = context.config
        source_path = Path(context.input.network)
        source = pypsa.Network(source_path)
        inputs_path = source_path.parent.parent / "gas_input_locations_simplified.csv"
        report["preflight"] = check_source(
            source, pd.read_csv(inputs_path, index_col=0)
        )
        report["git_commit"] = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip()
        report["source"] = str(source_path)
        report["source_sha256"] = hashlib.sha256(source_path.read_bytes()).hexdigest()
        hook = context.params.custom_extra_functionality
        if hook:
            hook = Path(hook)
            require(
                hook.is_file(),
                f"Configured custom constraints file is missing: {hook}; do not silently disable it",
            )
            report["custom_constraints_sha256"] = hashlib.sha256(
                hook.read_bytes()
            ).hexdigest()
        options = dict(context.params.solving["options"])
        require(
            not options.get("nhours") and not options.get("rolling_horizon"),
            "Use the diagnostic window, not nhours/rolling-horizon overrides",
        )
        require(
            not options.get("load_shedding", {}).get("enable")
            and not options.get("load_sinks", {}).get("enable"),
            "Diagnostic must not introduce balancing slack",
        )
        n = short_network(source, args.hours, args.start)
        prepare_network(
            n,
            config=config,
            solve_opts=options,
            foresight=context.params.foresight,
            planning_horizons="2035",
            co2_sequestration_potential=context.params.co2_sequestration_potential,
            limit_max_growth=context.params.sector["limit_max_growth"],
        )
        n.export_to_netcdf(args.output_dir / "input.nc")
        solving = context.params.solving
        model_kwargs, solve_kwargs = collect_kwargs(
            config, solving, "2035", log_fn=str(args.output_dir / "solver.log")
        )
        if args.method == "barrier":
            solver_options = dict(config["solving"]["solver_options"]["gurobi-default"])
        else:
            solver_options = dict(
                solving["solver_options"][solving["solver"]["options"]]
            )
            solver_options.update(method=6, crossover=0)
        solver_options.update(threads=args.threads, TimeLimit=args.time_limit)
        solve_kwargs.update(solver_name="gurobi", solver_options=solver_options)
        report["solver_options"] = solver_options
        create_optimization_model(
            n,
            config=config,
            params=context.params,
            model_kwargs=model_kwargs,
            solve_kwargs=solve_kwargs,
            planning_horizons="2035",
            snakemake=context,
        )
        report.update(
            variables=int(n.model.nvars),
            constraints=int(n.model.ncons),
            diagnostic=n.meta["diagnostic"],
        )
        status, condition = n.optimize.solve_model(**solve_kwargs)
        report.update(status=str(status), termination_condition=str(condition))
        require(
            str(status) == "ok" and str(condition) == "optimal",
            f"Diagnostic did not converge optimally: {status}/{condition}; inspect solver.log",
        )
        report["objective"] = float(n.objective)
        n.meta["solve_result"] = {
            "status": str(status),
            "termination_condition": str(condition),
        }
        n.export_to_netcdf(args.output_dir / "solved.nc")
    except Exception as error:
        report["error"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        report["elapsed_seconds"] = time.monotonic() - started
        report["peak_rss_GiB"] = (
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2
        )
        (args.output_dir / "summary.json").write_text(
            json.dumps(report, indent=2) + "\n"
        )


if __name__ == "__main__":
    main()
