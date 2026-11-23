# SPDX-License-Identifier: MIT
"""Solve matched Bornholm layouts from the audited final-v1 campaign inputs."""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from scripts import bornholm_2035 as campaign
from scripts import run_bornholm_2035 as original
from scripts import run_bornholm_welfare as runner

DESIGNS = (
    ("central_split_radial", "central_split_radial", 1000),
    ("central_hybrid_matched", "central_hybrid", 1000),
    ("central_hybrid", "central_hybrid", 1200),
)
LAYOUTS = (
    ("split_radial", "split_radial", 1000, 2000),
    ("hybrid", "hybrid", 1200, 2000),
    ("hybrid_matched", "hybrid", 1000, 2000),
    ("hybrid_more_dk", "hybrid", 1400, 1800),
    ("hybrid_more_de", "hybrid", 1000, 2200),
)
LAYOUT_WORLDS = ("central", "ev55", "ev65")
CENTRAL_FOLLOWUPS = (
    ("central_hybrid_matched", 1000, 2000, 3000),
    ("central_germany_alone", 0, 2000, 2000),
    ("central_denmark_alone", 1200, 0, 1200),
    ("central_smaller_island", 1200, 2000, 2000),
)
BOUNDS = ["p_nom_min", "p_nom_max"]
ROOT = "bornholm_2035_matched"
LAYOUT_ROOT = "bornholm_2035_layout_screen"
require = campaign.welfare.require


def project_legs(n, source_case, dk_mw, *, de_mw=2000, solved=False):
    split = source_case["project"] == "split_radial"
    legs = []
    for country, onshore, capacity in (("DK2", "DKE1", dk_mw), ("DE", "DE00", de_mw)):
        if not capacity:
            continue
        name = f"BEI_TO_{country}" + ("_SPLIT" if split else "")
        hub = "BEI_HUB" + ("_DK" if country == "DK2" else "_DE") if split else "BEI_HUB"
        reverse = name + "-reversed"
        a, b = n.links.loc[name], n.links.loc[reverse]
        require(
            (a.bus0, a.bus1, b.bus0, b.bus1) == (hub, onshore, onshore, hub),
            f"Unexpected terminals for {name}",
        )
        require(
            a.capital_cost > 0 and b.capital_cost == 0 and b.reversed,
            f"Unexpected investment accounting for {name}",
        )
        for leg in (a, b):
            require(
                leg.p_nom_extendable
                and leg.p_nom == 0
                and leg.p_nom_min == leg.p_nom_max == capacity
                and leg.p_min_pu == 0
                and leg.p_max_pu > 0
                and 0 < leg.efficiency <= 1,
                f"Unexpected directional capacity or losses for {leg.name}",
            )
            if solved:
                require(
                    np.isclose(leg.p_nom_opt, capacity, rtol=0, atol=0.01),
                    f"Wrong solved rating for {leg.name}",
                )
        legs.extend([name, reverse])
    require(
        set(n.links.index[n.links.index.str.startswith("BEI_TO_")]) == set(legs),
        "Unexpected Bornholm export assets",
    )


def resize_hybrid_legs(n, dk_mw, de_mw):
    """Resize both directions without changing the background grid budget."""
    changes = {
        "BEI_TO_DK2": dk_mw - n.links.at["BEI_TO_DK2", "p_nom_min"],
        "BEI_TO_DE": de_mw - n.links.at["BEI_TO_DE", "p_nom_min"],
    }
    for name, row in n.global_constraints.iterrows():
        unit = {
            "transmission_volume_expansion_limit": "length",
            "transmission_expansion_cost_limit": "capital_cost",
        }.get(row.type)
        if unit and "DC" in str(row.carrier_attribute).replace(" ", "").split(","):
            allowance = row.get("energy_island_allowance", np.nan)
            require(np.isfinite(allowance), f"Missing forced-link allowance for {name}")
            # The reverse copies have zero length and capital cost in PyPSA.
            delta = sum(
                change * n.links.at[leg, unit] for leg, change in changes.items()
            )
            n.global_constraints.at[name, "constant"] += delta
            n.global_constraints.at[name, "energy_island_allowance"] += delta
    for leg, capacity in (("BEI_TO_DK2", dk_mw), ("BEI_TO_DE", de_mw)):
        directions = [leg, leg + "-reversed"]
        if capacity:
            n.links.loc[directions, BOUNDS] = capacity
        else:
            n.remove("Link", directions)


def resize_danish_leg(n):
    """Resize the original Central matched control to 1 GW on the Danish shore."""
    resize_hybrid_legs(n, 1000, 2000)


def design_table(layout_screen=False):
    if not layout_screen:
        return [(name, parent, dk_mw, 2000) for name, parent, dk_mw in DESIGNS]
    return [
        (f"{world}_{name}", f"{world}_{parent}", dk_mw, de_mw)
        for world in LAYOUT_WORLDS
        for name, parent, dk_mw, de_mw in LAYOUTS
    ]


def validate_design(
    n, config, source_case, dk_mw, de_mw, *, wind_mw=3000, solved=False
):
    project_legs(n, source_case, dk_mw, de_mw=de_mw, solved=solved)
    if solved and source_case["project"] == "hybrid":
        require(
            np.isclose(
                n.generators.at["BEI_WIND", "p_nom_opt"], wind_mw, rtol=0, atol=0.01
            ),
            "Wrong solved Bornholm wind capacity",
        )
    return campaign.validate(
        n,
        config,
        source_case,
        solved,
        hybrid_dk_mw=dk_mw,
        hybrid_de_mw=de_mw,
        hybrid_wind_mw=wind_mw,
    )


def load_design(config, source_case, dk_mw, de_mw, *, wind_mw=3000):
    import pypsa

    parent = source_case["name"]
    source, source_dir = runner.paths(config, source_case)
    provenance = json.loads((source_dir / "prepared.json").read_text())
    require(
        provenance["scenario"] == source_case
        and provenance["campaign_config_sha256"] == runner.digest(campaign.SPEC),
        f"Original scenario or campaign settings differ for {parent}",
    )
    require(
        runner.digest(source) == provenance["prepared_sha256"],
        f"Changed archived input: {source}",
    )
    n = pypsa.Network(source)
    require("solve_result" not in n.meta, "Expected unsolved composed network")
    campaign.validate(n, config, source_case)
    default_dk_mw = 1000 if source_case["project"] == "split_radial" else 1200
    project_legs(n, source_case, default_dk_mw)
    resize_legs = (dk_mw, de_mw) != (default_dk_mw, 2000)
    changed = resize_legs or wind_mw != 3000
    if changed:
        require(source_case["project"] == "hybrid", "Only rescale the hybrid")
        require(
            not n.meta.get("energy_islands", {}).get(
                "subtract_generic_potential", False
            ),
            "Resizing expects unchanged background wind potential",
        )
        if resize_legs:
            resize_hybrid_legs(n, dk_mw, de_mw)
        n.generators.loc["BEI_WIND", BOUNDS] = wind_mw
    validate_design(n, config, source_case, dk_mw, de_mw, wind_mw=wind_mw)
    return n, source, provenance, changed


def comparison_inputs(n, source_case):
    constraints = n.global_constraints.drop(columns="mu", errors="ignore").copy()
    if "energy_island_allowance" in constraints:
        constraints["constant"] -= constraints.pop("energy_island_allowance").fillna(0)
    wind = "BEI_WIND_DK" if source_case["project"] == "split_radial" else "BEI_WIND"
    return (
        n.meta[campaign.AUDIT_KEY]["demand_sha256"],
        n.snapshot_weightings.copy(),
        n.generators_t.p_max_pu[wind].copy(),
        constraints.sort_index(),
    )


def check_layout_inputs():
    config = campaign.read_config()
    original.check_definitions(config)
    sources = {c["name"]: c for c in campaign.scenarios(config)}
    reference = {}
    for name, parent, dk_mw, de_mw in design_table(layout_screen=True):
        source_case = sources[parent]
        n, *_ = load_design(config, source_case, dk_mw, de_mw)
        world = source_case["world"]
        inputs = comparison_inputs(n, source_case)
        if world in reference:
            base = reference[world]
            require(
                inputs[0] == base[0]
                and inputs[1].equals(base[1])
                and np.array_equal(inputs[2], base[2]),
                f"Unmatched demand, snapshots or wind resource in {name}",
            )
            pd.testing.assert_frame_equal(inputs[3], base[3], rtol=1e-10, atol=1e-6)
        else:
            reference[world] = inputs
        print(f"Checked {name}: {dk_mw}/{de_mw} MW", flush=True)
        del n
    print("All fifteen inputs passed. No optimisation was run.")


def run(index, batch, *, layout_screen=False, central_followup=False):
    config = campaign.read_config()
    original.check_definitions(config)
    if central_followup:
        name, dk_mw, de_mw, wind_mw = CENTRAL_FOLLOWUPS[index]
        parent = "central_hybrid"
    else:
        name, parent, dk_mw, de_mw = design_table(layout_screen)[index]
        wind_mw = 3000
    sources = {c["name"]: c for c in campaign.scenarios(config)}
    source_case = sources[parent]
    n, source, provenance, changed = load_design(
        config, source_case, dk_mw, de_mw, wind_mw=wind_mw
    )
    trial = dict(source_case, index=index, name=name)
    if layout_screen or central_followup:
        trial.update(dk_mw=dk_mw, de_mw=de_mw, wind_mw=wind_mw)
    root = LAYOUT_ROOT if layout_screen else ROOT
    if central_followup:
        root = "bornholm_2035_central_followup"
    config["prefix"] = f"{root}/{batch}"
    target, directory = runner.paths(config, trial)
    directory.mkdir(parents=True, exist_ok=False)
    target.parent.mkdir(parents=True, exist_ok=True)
    if changed:
        temp = target.with_suffix(".tmp.nc")
        n.export_to_netcdf(temp)
        temp.replace(target)
    else:
        target.symlink_to(source.resolve())
    del n

    def validate(network, cfg, case, solved=False):
        return validate_design(
            network, cfg, source_case, dk_mw, de_mw, wind_mw=wind_mw, solved=solved
        )

    def context_for(case):
        context = original.context_for(source_case)
        # Scenario overrides have already been loaded by original.context_for.
        context.config["run"]["scenarios"]["enable"] = False
        options = context.params.solving["solver_options"]
        options[context.params.solving["solver"]["options"]].update(
            NumericFocus=2, BarConvTol=1e-9
        )
        identity = dict(name=name, parent=parent, dk_mw=dk_mw)
        if layout_screen:
            identity["de_mw"] = de_mw
        if central_followup:
            identity.update(de_mw=de_mw, wind_mw=wind_mw)
        context.config["matched_design"] = identity
        return context

    runner.write_json(
        directory / "prepared.json",
        dict(
            scenario=trial,
            original_input=str(source),
            source_sha256=provenance["prepared_sha256"],
            original_code_sha256=provenance["code_sha256"],
            prepared_sha256=runner.digest(target),
            campaign_config_sha256=runner.digest(campaign.SPEC),
            code_sha256=runner.code_fingerprint(),
            dk_leg_mw=dk_mw,
            de_leg_mw=de_mw,
            wind_mw=wind_mw,
        ),
    )
    status = runner.solve(
        config,
        campaign.SPEC,
        0,
        scenario_table=[trial],
        context_factory=context_for,
        validate_fn=validate,
    )

    result = json.loads((directory / "summary.json").read_text())
    if status in (0, 3):
        transfer = directory.parent / "transfer" / name
        transfer.mkdir(parents=True, exist_ok=True)
        for filename, src in (
            ("candidate_2035.nc", Path(result["candidate_network"])),
            ("summary.json", directory / "summary.json"),
            ("prepared.json", directory / "prepared.json"),
        ):
            (transfer / filename).symlink_to(src.resolve(), target_is_directory=False)
    return status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("index", type=int, nargs="?")
    scope = parser.add_mutually_exclusive_group()
    scope.add_argument("--layout-screen", action="store_true")
    scope.add_argument("--central-followup", action="store_true")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--plan", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    designs = (
        CENTRAL_FOLLOWUPS if args.central_followup else design_table(args.layout_screen)
    )
    if args.plan or args.check:
        if args.index is not None or (
            args.check and not (args.layout_screen or args.central_followup)
        ):
            parser.error(
                "Use --layout-screen or --central-followup with --check and no index"
            )
        if args.plan:
            for i, (name, second, third, fourth) in enumerate(designs):
                if args.central_followup:
                    print(f"{i}\t{name}\t{second}/{third} MW\twind={fourth} MW")
                else:
                    print(f"{i}\t{name}\t{third}/{fourth} MW")
        elif args.central_followup:
            config = campaign.read_config()
            original.check_definitions(config)
            source = next(
                c for c in campaign.scenarios(config) if c["name"] == "central_hybrid"
            )
            for name, dk, de, wind in CENTRAL_FOLLOWUPS:
                n, *_ = load_design(config, source, dk, de, wind_mw=wind)
                del n
                print(f"Checked {name}", flush=True)
        else:
            check_layout_inputs()
        return 0
    if args.index not in range(len(designs)):
        parser.error(f"Expected a task index from 0 to {len(designs) - 1}")
    batch = os.environ["SLURM_ARRAY_JOB_ID"]
    require(batch.isdecimal(), "Expected a numeric Slurm array job ID")
    return run(
        args.index,
        batch,
        layout_screen=args.layout_screen,
        central_followup=args.central_followup,
    )


if __name__ == "__main__":
    raise SystemExit(main())
