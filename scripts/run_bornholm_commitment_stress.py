# SPDX-License-Identifier: MIT
"""Run paired 2035 Bornholm commitments under demand-location and offshore stresses."""

import argparse
import json
import os
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd

from scripts import bornholm_2035 as campaign
from scripts import run_bornholm_2035 as original
from scripts import run_bornholm_welfare as runner
from scripts.run_bornholm_matched import project_legs, require, resize_danish_leg

ROOT = "bornholm_2035_commitment_stress"
LAYOUTS = (
    ("split_radial", "central_split_radial", 1000),
    ("treaty_hybrid", "central_hybrid", 1200),
    ("matched_hybrid", "central_hybrid", 1000),
)
STRESSES = ("demand_dk2", "demand_de", "competing_offshore")
ADDED_DEMAND_TWH = 5.0
ADDED_OFFSHORE_MW = 1000.0
OFFSHORE_CARRIERS = ("offwind-ac", "offwind-dc", "offwind-float")
CAPACITY_TOLERANCE = 0.01


def cases():
    return [
        dict(
            index=3 * i + j,
            name=f"{stress}_{layout}",
            stress=stress,
            layout=layout,
            parent=parent,
            dk_mw=dk_mw,
            de_mw=2000,
        )
        for i, stress in enumerate(STRESSES)
        for j, (layout, parent, dk_mw) in enumerate(LAYOUTS)
    ]


def load_original(config, source_case):
    from pypsa import Network

    name = source_case["name"]
    path, directory = runner.paths(config, source_case)
    record = json.loads((directory / "prepared.json").read_text())
    require(
        record["scenario"] == source_case
        and record["campaign_config_sha256"] == runner.digest(campaign.SPEC)
        and runner.digest(path) == record["prepared_sha256"],
        f"Original {name} settings or prepared network changed",
    )
    n = Network(path)
    require("solve_result" not in n.meta, f"Expected an unsolved network: {path}")
    campaign.validate(n, config, source_case)
    project_legs(n, source_case, 1000 if name == "central_split_radial" else 1200)
    return n, path, directory, record


def offshore_minimums(source, solved):
    """Reserve existing resource headroom above the same solved radial benchmark."""
    profile = source.generators_t.p_max_pu
    require(
        source.snapshots.equals(solved.snapshots), "Offshore reference snapshots differ"
    )
    project = profile["BEI_WIND_DK"]
    minimums, additions = {}, {}
    for region in ("DK2", "DE"):
        ids = campaign.welfare.selected_generators(
            source, {"region": region, "carriers": OFFSHORE_CARRIERS}
        )
        ids = ids.intersection(source.components["Generator"].active_assets)
        ids = ids[source.generators.loc[ids, "p_nom_extendable"]]
        ids = ids[~ids.str.startswith("BEI_WIND")]
        require(len(ids) > 0, f"No native extendable offshore sites in {region}")
        stock = source.generators.loc[ids, ["p_nom_min", "p_nom_max"]]
        require(
            ids.isin(solved.generators.index).all(),
            f"Missing reference sites in {region}",
        )
        reference = solved.generators.loc[ids, "p_nom_opt"]
        require(
            np.isfinite(reference).all()
            and reference.ge(stock.p_nom_min - CAPACITY_TOLERANCE).all()
            and reference.le(stock.p_nom_max + CAPACITY_TOLERANCE).all(),
            f"Invalid solved offshore capacity in {region}",
        )
        base = reference.clip(lower=stock.p_nom_min, upper=stock.p_nom_max)
        headroom = stock.p_nom_max - base
        require(
            np.isfinite(headroom).all(),
            f"Expected finite offshore resource bounds in {region}",
        )
        eligible = ids.intersection(profile.columns)
        score = profile[eligible].corrwith(project).dropna()
        score = score.sort_index().sort_values(ascending=False, kind="stable")
        require(
            headroom.reindex(score.index).sum() >= ADDED_OFFSHORE_MW,
            f"Insufficient profiled offshore headroom for {ADDED_OFFSHORE_MW:g} MW in {region}",
        )
        remaining = ADDED_OFFSHORE_MW
        for name in score.index:
            extra = min(float(headroom[name]), remaining)
            if extra > 0:
                base[name] += extra
                additions[name] = dict(
                    region=region, extra_mw=extra, correlation=float(score[name])
                )
                remaining -= extra
            if remaining <= 1e-8:
                break
        minimums.update({name: float(value) for name, value in base.items()})
    return dict(minimums_mw=minimums, added_sites=additions)


def read_offshore_plan(config, sources):
    from pypsa import Network

    n, _, directory, prepared = load_original(config, sources["central_split_radial"])
    summary = json.loads((directory / "summary.json").read_text())
    require(
        summary["status"] in {"checks_passed", "optimal_needs_review"}
        and summary["source_sha256"] == prepared["prepared_sha256"]
        and summary["scenario"] == sources["central_split_radial"],
        "Use the matching solved Central split-radial candidate for the offshore stress",
    )
    solved_path = Path(summary.get("accepted_network") or summary["candidate_network"])
    solved = Network(solved_path)
    require(
        solved.objective is not None
        and abs(solved.objective - summary["objective"]) <= 2,
        "Central split-radial objective differs from its summary",
    )
    plan = offshore_minimums(n, solved)
    plan.update(reference=str(solved_path), reference_status=summary["status"])
    return plan


def stressed_config(config, stress):
    changed = deepcopy(config)
    if stress in {"demand_dk2", "demand_de"}:
        region = "DK2" if stress == "demand_dk2" else "DE"
        loads = changed["worlds"]["central"]["additional_data_centres_TWh"]
        require(
            not any(loads.values()), "Expected zero added load in the Central parent"
        )
        loads[region] = ADDED_DEMAND_TWH
    return changed


def validate_stress(n, config, source_case, design, offshore, *, solved=False):
    project_legs(n, source_case, design["dk_mw"], solved=solved)
    report = campaign.validate(
        n, config, source_case, solved, hybrid_dk_mw=design["dk_mw"]
    )
    if offshore is not None:
        bounds = pd.Series(offshore["minimums_mw"], dtype=float)
        column = "p_nom_opt" if solved else "p_nom_min"
        capacities = n.generators.loc[bounds.index, column]
        valid = (
            capacities.ge(bounds - CAPACITY_TOLERANCE).all()
            if solved
            else np.allclose(capacities, bounds, rtol=0, atol=1e-8)
        )
        require(valid, "Competing offshore capacity differs from its commitment")
    return report


def apply_stress(n, config, source_case, design, offshore):
    """Modify the already checked native inputs, then repeat the existing audits."""
    if design["layout"] == "matched_hybrid":
        resize_danish_leg(n)
    if design["stress"] in {"demand_dk2", "demand_de"}:
        n.meta.pop(campaign.AUDIT_KEY)
        campaign.apply(
            n, campaign.definitions(config)[source_case["name"]]["bornholm_2035"]
        )
    if offshore is not None:
        bounds = pd.Series(offshore["minimums_mw"], dtype=float)
        g = n.generators.loc[bounds.index]
        require(
            g.p_nom_extendable.all()
            and bounds.ge(g.p_nom_min).all()
            and bounds.le(g.p_nom_max).all(),
            "Offshore commitments cannot fit the parent generator bounds",
        )
        n.generators.loc[bounds.index, "p_nom_min"] = bounds
    validate_stress(n, config, source_case, design, offshore)


def check(config, sources):
    """Stop before submitting expensive solves if archived inputs or bounds are wrong."""
    offshore = read_offshore_plan(config, sources)
    by_stress = {}
    for design in cases():
        n, _, _, _ = load_original(config, sources[design["parent"]])
        stress = design["stress"]
        apply_stress(
            n,
            stressed_config(config, stress),
            sources[design["parent"]],
            design,
            offshore if stress == "competing_offshore" else None,
        )
        signature = (
            n.meta[campaign.AUDIT_KEY]["demand_sha256"],
            n.snapshot_weightings.copy(),
        )
        if stress in by_stress:
            prior, weights = by_stress[stress]
            require(
                signature[0] == prior and signature[1].equals(weights),
                f"Demand or weights differ across layouts for {stress}",
            )
        else:
            by_stress[stress] = signature
        print(f"Checked {design['name']}", flush=True)
        del n
    print("All nine inputs passed. No optimisation was run.")


def run(config, sources, design, batch):
    source_case = sources[design["parent"]]
    n, source, _, prepared = load_original(config, source_case)
    stress = design["stress"]
    offshore = (
        read_offshore_plan(config, sources) if stress == "competing_offshore" else None
    )
    cfg = stressed_config(config, stress)
    apply_stress(n, cfg, source_case, design, offshore)
    trial = dict(source_case, **design)
    cfg["prefix"] = f"{ROOT}/{batch}"
    target, directory = runner.paths(cfg, trial)
    directory.mkdir(parents=True, exist_ok=False)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(".tmp.nc")
    n.export_to_netcdf(temporary)
    temporary.replace(target)
    del n

    def validate(network, cfg, case, solved=False):
        return validate_stress(
            network, cfg, source_case, design, offshore, solved=solved
        )

    def context_for(case):
        context = original.context_for(source_case)
        context.config["run"]["scenarios"]["enable"] = False
        context.config["bornholm_2035"] = campaign.definitions(cfg)[
            source_case["name"]
        ]["bornholm_2035"]
        options = context.params.solving["solver_options"]
        options[context.params.solving["solver"]["options"]].update(
            NumericFocus=2, BarConvTol=1e-9
        )
        context.config["bornholm_commitment_stress"] = dict(
            design=design, offshore=offshore
        )
        return context

    runner.write_json(
        directory / "prepared.json",
        dict(
            scenario=trial,
            original_input=str(source),
            offshore=offshore,
            source_sha256=prepared["prepared_sha256"],
            prepared_sha256=runner.digest(target),
            campaign_config_sha256=runner.digest(campaign.SPEC),
            code_sha256=runner.code_fingerprint(),
        ),
    )
    status = runner.solve(
        cfg,
        campaign.SPEC,
        0,
        scenario_table=[trial],
        context_factory=context_for,
        validate_fn=validate,
    )
    if status in (0, 3):
        summary = directory / "summary.json"
        record = json.loads(summary.read_text())
        transfer = directory.parent / "transfer" / trial["name"]
        transfer.mkdir(parents=True, exist_ok=False)
        for name, source_path in (
            ("candidate_2035.nc", Path(record["candidate_network"])),
            ("summary.json", summary),
            ("prepared.json", directory / "prepared.json"),
        ):
            (transfer / name).symlink_to(source_path.resolve())
    return status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--plan", action="store_true")
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--index", type=int, choices=range(len(cases())))
    args = parser.parse_args()
    matrix = cases()
    if args.plan:
        for design in matrix:
            print(
                f"{design['index']}\t{design['name']}\t{design['dk_mw']}/{design['de_mw']} MW"
            )
        return 0
    config = campaign.read_config()
    original.check_definitions(config)
    sources = {case["name"]: case for case in campaign.scenarios(config)}
    if args.check:
        check(config, sources)
        return 0
    batch = os.environ["SLURM_ARRAY_JOB_ID"]
    require(batch.isdecimal(), "Expected a numeric Slurm array job ID")
    return run(config, sources, matrix[args.index], batch)


if __name__ == "__main__":
    raise SystemExit(main())
