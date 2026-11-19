# SPDX-License-Identifier: MIT
"""Re-solve the central split radial, matched hybrid and treaty hybrid from audited inputs."""

import argparse
import json
import os
from pathlib import Path

import numpy as np

from scripts import bornholm_2035 as campaign
from scripts import run_bornholm_2035 as original
from scripts import run_bornholm_welfare as runner

DESIGNS = (
    ("central_split_radial", "central_split_radial", 1000),
    ("central_hybrid_matched", "central_hybrid", 1000),
    ("central_hybrid", "central_hybrid", 1200),
)
BOUNDS = ["p_nom_min", "p_nom_max"]
ROOT = "bornholm_2035_matched"
require = campaign.welfare.require


def project_legs(n, source_case, dk_mw, *, solved=False):
    split = source_case["project"] == "split_radial"
    legs = []
    for country, onshore, capacity in (("DK2", "DKE1", dk_mw), ("DE", "DE00", 2000)):
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


def resize_danish_leg(n):
    """Preserve the background transmission budget if the native cap is enabled."""
    forward = "BEI_TO_DK2"
    rows = [forward, forward + "-reversed"]
    for name, row in n.global_constraints.iterrows():
        unit = {
            "transmission_volume_expansion_limit": "length",
            "transmission_expansion_cost_limit": "capital_cost",
        }.get(row.type)
        if unit and "DC" in str(row.carrier_attribute).replace(" ", "").split(","):
            allowance = row.get("energy_island_allowance", np.nan)
            require(np.isfinite(allowance), f"Missing forced-link allowance for {name}")
            delta = 200 * n.links.at[forward, unit]
            n.global_constraints.at[name, "constant"] -= delta
            n.global_constraints.at[name, "energy_island_allowance"] -= delta
    n.links.loc[rows, BOUNDS] = 1000.0


def run(index, batch):
    import pypsa

    config = campaign.read_config()
    original.check_definitions(config)
    name, parent, dk_mw = DESIGNS[index]
    cases = {c["name"]: c for c in campaign.scenarios(config)}
    source_case = cases[parent]
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
    trial = dict(source_case, index=index, name=name)
    config["prefix"] = f"{ROOT}/{batch}"
    target, directory = runner.paths(config, trial)
    directory.mkdir(parents=True, exist_ok=False)
    target.parent.mkdir(parents=True, exist_ok=True)
    n = pypsa.Network(source)
    require("solve_result" not in n.meta, "Expected unsolved composed network")
    campaign.validate(n, config, source_case)
    project_legs(
        n, source_case, 1000 if source_case["project"] == "split_radial" else 1200
    )
    if name == "central_hybrid_matched":
        resize_danish_leg(n)
        temp = target.with_suffix(".tmp.nc")
        n.export_to_netcdf(temp)
        temp.replace(target)
    else:
        target.symlink_to(source.resolve())
    del n

    def validate(network, cfg, case, solved=False):
        project_legs(network, source_case, dk_mw, solved=solved)
        return campaign.validate(network, cfg, source_case, solved, hybrid_dk_mw=dk_mw)

    def context_for(case):
        context = original.context_for(source_case)
        # Scenario overrides have already been loaded by original.context_for.
        context.config["run"]["scenarios"]["enable"] = False
        options = context.params.solving["solver_options"]
        options[context.params.solving["solver"]["options"]].update(
            NumericFocus=2, BarConvTol=1e-9
        )
        context.config["matched_design"] = dict(name=name, parent=parent, dk_mw=dk_mw)
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
        ):
            (transfer / filename).symlink_to(src.resolve(), target_is_directory=False)
    return status


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("index", choices=range(len(DESIGNS)), type=int)
    args = parser.parse_args()
    batch = os.environ["SLURM_ARRAY_JOB_ID"]
    if not batch.isdecimal():
        raise ValueError("Expected a numeric Slurm array job ID")
    raise SystemExit(run(args.index, batch))
