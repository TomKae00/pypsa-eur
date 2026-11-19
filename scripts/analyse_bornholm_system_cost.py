# SPDX-License-Identifier: MIT
"""Attribute 2035 Bornholm system-cost changes directly from solved networks."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa

from scripts.audit_bornholm_candidate import Audit
from scripts.bornholm_welfare_sensitivity import countries, demand_panel, panel_digest

PATHWAYS = ("central", "ev55", "ev65")
DESIGNS = ("baseline", "hybrid")
GROUPS = ["scope", "region", "component", "carrier", "cost_type"]


def asset_regions(n, component, assets):
    ports = n.components[component].static.loc[assets].filter(regex=r"^bus\d*$")
    endpoint_countries = ports.apply(lambda buses: countries(n, buses))
    endpoint_countries = endpoint_countries.where(
        endpoint_countries.apply(
            lambda values: (
                values.str.fullmatch("[A-Z]{2}", na=False) & ~values.isin(["EU", "NA"])
            )
        )
    )

    def region(row):
        zones = row.dropna().unique()
        if len(zones) == 1:
            return zones[0]
        return "Cross-border" if len(zones) else "Unallocated"

    return endpoint_countries.apply(region, axis=1)


def read_case(path, design):
    n = pypsa.Network(path)
    if isinstance(n.snapshots, pd.MultiIndex) or not np.isclose(
        n.snapshot_weightings.objective.sum(), 8760
    ):
        raise ValueError(f"Expected one complete model year in {path}")
    if n.objective is None or not np.isfinite(n.objective):
        raise ValueError(f"Missing solved objective in {path}")

    costs, audit = Audit(n).costs()
    if abs(audit["reconstruction_error_eur"]) > 2:
        raise ValueError(f"Objective reconstruction failed for {path}: {audit}")

    costs["cost_eur"] = costs["gross_eur"]
    if n.objective_constant:
        costs["cost_eur"] -= costs["existing_constant_eur"]
    costs["cost_type"] = costs["cost_type"].where(
        costs["cost_type"].eq("investment"), "operation"
    )
    costs["scope"] = "Background"
    wind = costs.component.eq("Generator") & costs.asset.str.startswith("BEI_WIND")
    cable = costs.component.eq("Link") & costs.asset.str.startswith("BEI_TO_")
    if (wind.any(), cable.any()) != (design == "hybrid", design == "hybrid"):
        raise ValueError(f"Unexpected Bornholm wind or export assets in {path}")
    costs.loc[wind, "scope"] = "Bornholm wind"
    costs.loc[cable, "scope"] = "Bornholm cables"
    for component, rows in costs.groupby("component"):
        region = asset_regions(n, component, rows.asset.unique())
        costs.loc[rows.index, "region"] = rows.asset.map(region)

    load_ids = n.loads.index.sort_values()
    signature = (
        panel_digest(demand_panel(n, load_ids, nonnegative=False)),
        n.loads.loc[load_ids, ["bus", "carrier", "sign"]].copy(),
        n.snapshot_weightings.copy(),
    )
    return costs.groupby(GROUPS).cost_eur.sum(), audit, signature


def compare(root):
    cases = [(world, design) for world in PATHWAYS for design in DESIGNS]
    paths = {
        (world, design): Path(root) / f"{world}_{design}" / "candidate_2035.nc"
        for world, design in cases
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing solved networks:\n" + "\n".join(missing))

    ledger, records = {}, []
    for world in PATHWAYS:
        baseline, base_audit, base_signature = read_case(
            paths[world, "baseline"], "baseline"
        )
        hybrid, hybrid_audit, hybrid_signature = read_case(
            paths[world, "hybrid"], "hybrid"
        )
        if base_signature[0] != hybrid_signature[0] or any(
            not left.equals(right)
            for left, right in zip(base_signature[1:], hybrid_signature[1:])
        ):
            raise ValueError(f"Fixed loads or snapshot weights differ in {world}")
        ledger["baseline", world] = baseline
        ledger["hybrid impact", world] = hybrid.subtract(baseline, fill_value=0)
        scope = hybrid.groupby(level="scope").sum()
        project = scope.drop(index="Background").sum()
        no_island_cost = base_audit["saved_objective_eur"]
        hybrid_cost = hybrid_audit["saved_objective_eur"]
        records.append(
            {
                "pathway": world,
                "no_island_cost_eur": no_island_cost,
                "hybrid_cost_eur": hybrid_cost,
                "hybrid_minus_no_island_eur": hybrid_cost - no_island_cost,
                "hybrid_project_cost_eur": project,
                "background_cost_avoided_eur": project - (hybrid_cost - no_island_cost),
                "baseline_reconstruction_error_eur": base_audit[
                    "reconstruction_error_eur"
                ],
                "hybrid_reconstruction_error_eur": hybrid_audit[
                    "reconstruction_error_eur"
                ],
            }
        )

    summary = pd.DataFrame(records).set_index("pathway")
    summary["hybrid_impact_change_vs_central_eur"] = (
        summary.hybrid_minus_no_island_eur
        - summary.at["central", "hybrid_minus_no_island_eur"]
    )
    comparisons = {}
    for world in PATHWAYS:
        impact = ledger["hybrid impact", world]
        comparisons["hybrid minus no island", world] = impact
        if world != "central":
            comparisons["no island minus central no island", world] = ledger[
                "baseline", world
            ].subtract(ledger["baseline", "central"], fill_value=0)
            comparisons["hybrid impact minus central hybrid impact", world] = (
                impact.subtract(ledger["hybrid impact", "central"], fill_value=0)
            )
    drivers = (
        pd.concat(comparisons, names=["comparison", "pathway"])
        .rename("cost_change_eur")
        .reset_index()
    )
    return summary, drivers


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("transfer", type=Path, help="Final-v1 campaign transfer folder")
    parser.add_argument("--output", type=Path, required=True, help="Output directory")
    args = parser.parse_args()
    summary, drivers = compare(args.transfer)
    args.output.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output / "system_cost_summary.csv")
    drivers.to_csv(args.output / "system_cost_drivers.csv", index=False)
    shown = summary.drop(columns=list(summary.filter(like="reconstruction_error")))
    shown = shown.div(1e6).rename(columns=lambda name: name.removesuffix("_eur"))
    print("Annual system-cost bridge [million real 2025 EUR]")
    print(shown.round(2).to_string())
    for (comparison, world), rows in drivers.groupby(
        ["comparison", "pathway"], sort=False
    ):
        if comparison == "hybrid minus no island":
            continue
        top = rows.loc[rows.cost_change_eur.abs().nlargest(8).index]
        top = top.assign(cost_change_million_eur=top.cost_change_eur / 1e6)
        print(f"\n{world}: {comparison}")
        print(top[[*GROUPS, "cost_change_million_eur"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
