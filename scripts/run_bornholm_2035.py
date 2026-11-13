# SPDX-License-Identifier: MIT
"""Build checks, native annual solves and result collection for the 15-case campaign."""

import argparse
import json
import shutil
from pathlib import Path

import pandas as pd
import yaml

from scripts import bornholm_2035 as campaign
from scripts import run_bornholm_calibrated as calibrated
from scripts import run_bornholm_study as study
from scripts import run_bornholm_welfare as runner

CONFIGS = [*calibrated.CONFIGS, "config/config.bornholm-2035.yaml"]


def context_for(case):
    return study.context_for(case["name"], configfiles=CONFIGS)


def check_definitions(config):
    actual = yaml.safe_load(campaign.SCENARIOS.read_text())
    campaign.welfare.require(
        actual == campaign.definitions(config),
        "Scenario definitions differ from the campaign specification",
    )


def prepare(config):
    import pypsa

    code_sha = runner.code_fingerprint()
    config_sha = runner.digest(campaign.SPEC)
    records, demand_rows, references = [], [], {}
    for case in campaign.scenarios(config):
        context = context_for(case)
        source, directory = runner.paths(config, case)
        campaign.welfare.require(
            Path(context.input.network).resolve() == source.resolve(),
            "Native composition path does not match the campaign",
        )
        campaign.welfare.require(source.is_file(), f"Missing composed input: {source}")
        preflight = calibrated.preflight(
            context, configfiles=CONFIGS, prefix=config["prefix"]
        )
        n = pypsa.Network(source)
        campaign.welfare.require(
            "solve_result" not in n.meta, "Use unsolved composed inputs"
        )
        limits = campaign.validate(n, config, case)
        group = (case["world"], case["v2g_share"])
        if case["project"] == "baseline":
            # Save only demands/budgets, not several full European networks.
            references[group] = (
                n.loads.copy(),
                campaign.welfare.demand_panel(n, n.loads.index, nonnegative=False),
                n.snapshot_weightings.copy(),
                n.global_constraints.drop(columns="mu", errors="ignore").copy(),
            )
        else:
            loads, panel, weights, budgets = references[group]
            campaign.welfare.require(
                n.loads.index.equals(loads.index), "Paired load names differ"
            )
            campaign.welfare.require(
                campaign.welfare.demand_panel(
                    n, n.loads.index, nonnegative=False
                ).equals(panel),
                "Paired final service demand differs",
            )
            campaign.welfare.require(
                n.snapshot_weightings.equals(weights), "Paired weights differ"
            )
            campaign.welfare.require(
                n.global_constraints.drop(columns="mu", errors="ignore").equals(
                    budgets
                ),
                "Paired emissions/resource budgets differ",
            )
        record = dict(
            scenario=case,
            campaign_config_sha256=config_sha,
            code_sha256=code_sha,
            prepared_network=str(source),
            prepared_sha256=runner.digest(source),
            preflight=preflight,
            capacity_limits=limits,
            project=campaign.project_check(n, case["project"]),
        )
        marker = directory / "prepared.json"
        if marker.exists():
            old = json.loads(marker.read_text())
            campaign.welfare.require(
                all(
                    old.get(k) == record[k]
                    for k in (
                        "scenario",
                        "campaign_config_sha256",
                        "code_sha256",
                        "prepared_sha256",
                    )
                ),
                "Prepared inputs changed; use a new prefix rather than overwriting",
            )
        runner.write_json(marker, record)
        records.append(record)
        demand_rows.extend(
            dict(case=case["name"], **row) for row in campaign.demand_report(n)
        )
        print(f"Checked {case['index']:02d}: {case['name']}", flush=True)
        del n
    root = Path("results") / config["prefix"]
    pd.DataFrame(demand_rows).to_csv(root / "input_demand_by_carrier.csv", index=False)
    runner.write_json(root / "campaign.json", dict(config=config, prepared=records))
    print("All 15 composed inputs and matched final-demand/budget checks passed.")


def collect(config):
    """Retain flagged candidates and failed-case diagnostics; no CfD accounting here."""
    import pypsa

    root = Path("results") / config["prefix"]
    rows, energy, conversion_rows = [], [], []
    for case in campaign.scenarios(config):
        _, directory = runner.paths(config, case)
        summary_path = directory / "summary.json"
        summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
        rows.append(
            dict(
                index=case["index"],
                case=case["name"],
                status=summary.get("status", "missing"),
                objective=summary.get("objective"),
                error=summary.get("error", ""),
            )
        )
        transfer = root / "transfer" / case["name"]
        transfer.mkdir(parents=True, exist_ok=True)
        for path in (summary_path, directory / "prepared.json"):
            if path.exists():
                shutil.copy2(path, transfer / path.name)
        candidate = summary.get("accepted_network") or summary.get("candidate_network")
        if not candidate or not Path(candidate).is_file():
            continue
        candidate = Path(candidate)
        link = transfer / "candidate_2035.nc"
        if link.is_symlink():
            link.unlink()
        if not link.exists():
            link.symlink_to(candidate.resolve())
        effective = (
            Path(summary.get("attempt", candidate.parent)) / "effective_config.yaml"
        )
        if effective.is_file():
            shutil.copy2(effective, transfer / effective.name)
        n = pypsa.Network(candidate)
        conversion_rows.extend(
            dict(case=case["name"], **row) for row in campaign.electricity_use(n)
        )
        from scripts.bornholm_welfare_sensitivity import countries, demand_panel

        w = n.snapshot_weightings.generators
        panel = demand_panel(n, n.loads.index, nonnegative=False)
        for ct in ("DK", "DE"):
            country_loads = countries(n, n.loads.bus).eq(ct)
            electric = n.loads.bus.map(n.buses.carrier).isin(["AC", "low voltage"])
            fixed = float(
                panel.loc[:, country_loads & electric].mul(w, axis=0).sum().sum() / 1e6
            )
            energy.append(
                dict(case=case["name"], country=ct, fixed_electricity_TWh=fixed)
            )
        del n
    root.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(root / "campaign_status.csv", index=False)
    pd.DataFrame(energy).to_csv(root / "fixed_electricity_results.csv", index=False)
    if conversion_rows:
        detail = pd.DataFrame(conversion_rows)
        detail.to_csv(root / "electricity_use_by_carrier.csv", index=False)
        totals = detail.groupby(
            ["case", "country"], as_index=False
        ).electricity_TWh.sum()
        reference = (
            totals[totals.case.eq("central_baseline")]
            .set_index("country")
            .electricity_TWh
        )
        totals["change_from_central_baseline_percent"] = 100 * (
            totals.electricity_TWh / totals.country.map(reference) - 1
        )
        totals.to_csv(root / "electricity_use_totals.csv", index=False)
    for name in (
        "campaign.json",
        "input_demand_by_carrier.csv",
        "campaign_status.csv",
        "electricity_use_by_carrier.csv",
        "electricity_use_totals.csv",
    ):
        source = root / name
        if source.exists():
            shutil.copy2(source, root / "transfer" / name)
    print(pd.DataFrame(rows).to_string(index=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    for action in (
        "plan",
        "list-build-targets",
        "prepare",
        "collect",
        "write-scenarios",
    ):
        actions.add_argument(f"--{action}", action="store_true")
    actions.add_argument("--index", type=int, choices=range(15))
    args = parser.parse_args()
    config = campaign.read_config()
    cases = campaign.scenarios(config)
    if args.write_scenarios:
        campaign.SCENARIOS.write_text(
            yaml.safe_dump(campaign.definitions(config), sort_keys=False)
        )
        return 0
    check_definitions(config)
    if args.plan:
        for case in cases:
            world = config["worlds"][case["world"]]
            print(
                f"{case['index']:02d}\t{case['name']}\tEV={world['ev_share']:.0%}\t"
                f"industry={world['industrial_progress']:.0%}\tV2G={case['v2g_share']:.0%}\t"
                f"extraDC={world['additional_data_centres_TWh']}"
            )
    elif args.list_build_targets:
        for case in cases:
            print(runner.paths(config, case)[0])
    elif args.prepare:
        prepare(config)
    elif args.collect:
        collect(config)
    else:
        return runner.solve(
            config,
            campaign.SPEC,
            args.index,
            scenario_table=cases,
            context_factory=context_for,
            validate_fn=campaign.validate,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
