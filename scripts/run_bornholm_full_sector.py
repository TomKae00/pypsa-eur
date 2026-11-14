# SPDX-License-Identifier: MIT
"""One central 2035 Hybrid run with native H2 pipelines, shipping and aviation."""

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from scripts import bornholm_2035 as campaign
from scripts import run_bornholm_2035 as final_campaign
from scripts import run_bornholm_calibrated as calibrated
from scripts import run_bornholm_hourly as annual
from scripts import run_bornholm_study as study

RUN = "central_hybrid"
PREFIX = "bornholm_2035_full_sector_v1"
CONFIGS = [*final_campaign.CONFIGS, "config/config.bornholm-full-sector.yaml"]
TARGET = f"resources/{PREFIX}/{RUN}/networks/composed_2035.nc"
RESTORED = ("H2_network", "H2_retrofit", "shipping", "aviation")


def check_scope(n):
    """Reject stale simplified inputs while allowing native gas/H2 conversion."""
    study.require("solve_result" not in n.meta, "Use an unsolved composed input")
    study.require(n.meta.get("run", {}).get("prefix") == PREFIX, "Wrong namespace")
    study.require(
        n.snapshots.as_unit("ns").equals(
            pd.date_range("2013-01-01", periods=2920, freq="3h").as_unit("ns")
        )
        and np.allclose(n.snapshot_weightings, 3),
        "Expected the complete 2013 weather year at three-hour resolution",
    )
    sector = n.meta["sector"]
    for key in (*RESTORED, "co2_spatial", "co2_network", "gas_network"):
        study.require(sector.get(key) is True, f"Expected sector.{key}=true")
    study.require(
        n.meta["electricity"]["transmission_limit"] == "vfixed"
        and not n.lines.s_nom_extendable.any(),
        "The background electricity grid must remain fixed",
    )
    mainland = ~n.buses.index.str.startswith("BEI_")
    study.require(
        n.buses.loc[mainland, "carrier"].eq("AC").sum() == 52,
        "Expected the reviewed 52 mainland AC buses",
    )
    for technology, count in calibrated.CLASSES.items():
        study.require(
            n.meta["renewable"][technology]["resource_classes"] == count,
            f"Wrong {technology} resource classes",
        )
    for carrier in ("H2 pipeline", "H2 pipeline retrofitted", "CO2 pipeline"):
        pipes = n.links.loc[n.links.carrier.eq(carrier)]
        study.require(
            not pipes.empty and pipes.p_nom_extendable.all(),
            f"Missing native investment candidates: {carrier}",
        )
    gas = n.links.loc[n.links.carrier.eq("gas pipeline") & n.links.p_nom.gt(0)]
    study.require(
        not gas.empty
        and gas.p_nom_extendable.all()
        and np.allclose(gas.p_nom_min, 0)
        and np.allclose(gas.p_nom_max, gas.p_nom),
        "Existing gas pipes must retain the native bounded retrofit formulation",
    )
    energy = (
        n.get_switchable_as_dense("Load", "p_set")
        .mul(n.snapshot_weightings.generators, axis=0)
        .sum()
    )
    demand = energy.groupby(n.loads.carrier).sum().div(1e6)
    for carrier in ("shipping oil", "shipping methanol", "kerosene for aviation"):
        study.require(demand.get(carrier, 0) > 0, f"Missing native {carrier} demand")
    budgets = n.global_constraints.loc[n.global_constraints.type.eq("co2_atmosphere")]
    study.require(not budgets.empty, "Missing native carbon budget")
    return {
        "restored_switches": {key: sector[key] for key in RESTORED},
        "snapshots": len(n.snapshots),
        "annual_weight_hours": float(n.snapshot_weightings.objective.sum()),
        "links_by_carrier": n.links.carrier.value_counts().to_dict(),
        "final_demand_TWh_by_carrier": demand.to_dict(),
        "carbon_budget": budgets[
            ["type", "carrier_attribute", "sense", "constant"]
        ].to_dict("index"),
        "components": {c.list_name: len(c.static) for c in n.components},
    }


def preflight(context):
    import pypsa

    from scripts import transmission_assumptions
    from scripts.bornholm_model_assumptions import validate_assumptions

    n = pypsa.Network(context.input.network)
    report = check_scope(n)
    spec = campaign.read_config()
    case = next(case for case in campaign.scenarios(spec) if case["name"] == RUN)
    report["capacity_limits"] = campaign.validate(n, spec, case)
    report["biomass"] = calibrated.check_biomass_fuel(n)
    audit = validate_assumptions(n, context.config["bornholm_model_assumptions"])
    transmission_assumptions.validate(n, context.config["transmission_assumptions"])
    study.require(audit["transfers"]["planning_year"] == 2035, "Wrong grid year")
    report["assumption_fingerprint"] = audit["fingerprint"]
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-target", action="store_true")
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--time-limit", type=int, default=252000)
    args = parser.parse_args()
    if args.build_target:
        print(TARGET)
    elif args.preflight:
        print(
            json.dumps(preflight(study.context_for(RUN, configfiles=CONFIGS)), indent=2)
        )
    else:
        if args.output_dir is None:
            parser.error("--output-dir is required")
        logging.basicConfig(level=logging.INFO)
        annual.solve(
            args.output_dir,
            threads=16,
            time_limit=args.time_limit,
            annual_run=annual.AnnualRun(
                RUN,
                PREFIX,
                tuple(CONFIGS),
                TARGET,
                preflight,
                lambda n: campaign.project_check(n, "hybrid"),
            ),
        )


if __name__ == "__main__":
    main()
