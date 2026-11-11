# SPDX-License-Identifier: MIT
"""Auditable demand perturbations and aggregate capacity ceilings for fixed BEI."""

import hashlib
import json

import numpy as np
import pandas as pd

AUDIT_KEY = "bornholm_welfare_sensitivity_audit"
SOLAR_GROUND = ["solar", "solar-hsat"]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def scenarios(config):
    """Stable adjacent baseline/hybrid pairs; IDs are also Slurm array indices."""
    result = []
    for deployment in ("central", "slow_onshore"):
        for demand in ("central", "plus10", "plus20"):
            pair = f"{deployment}_{demand}"
            for project in ("baseline", "hybrid"):
                result.append(
                    dict(
                        index=len(result),
                        name=f"{pair}_{project}",
                        pair=pair,
                        deployment=deployment,
                        demand=demand,
                        demand_factor=float(config["demand"]["factors"][demand]),
                        project=project,
                        source_run=f"bornholm_zonal_{project}",
                    )
                )
    return result


def capacity_limits(config, scenario):
    d = config["deployment"][scenario["deployment"]]
    pv = config["solar"]
    ratio = float(pv["DK_DC_AC_ratio"])
    require(np.isfinite(ratio) and ratio >= 1, "PV DC/AC ratio must be finite and >= 1")
    limits = [
        dict(
            id=f"onshore_{region}",
            region=region,
            carriers=["onwind"],
            maximum_MW=float(d[f"onshore_{region}_MW"]),
        )
        for region in ("DK1", "DK2", "DE")
    ]
    for region in ("DK1", "DK2"):
        for kind, carriers in (
            ("rooftop", ["solar rooftop"]),
            ("ground", SOLAR_GROUND),
        ):
            limits.append(
                dict(
                    id=f"solar_{region}_{kind}",
                    region=region,
                    carriers=carriers,
                    maximum_MW=float(pv[f"{region}_{kind}_MW_AC"]) * ratio,
                )
            )
    limits.append(
        dict(
            id="solar_DE",
            region="DE",
            carriers=[*SOLAR_GROUND, "solar rooftop"],
            maximum_MW=float(pv["DE_total_MW"]),
        )
    )
    return limits


def countries(n, buses):
    country = buses.map(n.buses.country)
    parent = buses.map(n.buses.location).map(n.buses.country)
    return country.where(country.fillna("").str.fullmatch("[A-Z]{2}"), parent)


def selected_generators(n, limit):
    bus = n.generators.bus
    location = bus.map(n.buses.location).fillna("")
    region = limit["region"]
    if region in ("DK1", "DK2"):
        prefix = {"DK1": "DKW1", "DK2": "DKE1"}[region]
        geographic = bus.str.startswith(prefix) | location.str.startswith(prefix)
    else:
        geographic = countries(n, bus).eq(region)
    return n.generators.index[geographic & n.generators.carrier.isin(limit["carriers"])]


def capacity_report(n, limits, solved=False):
    rows = []
    active = n.components["Generator"].active_assets
    for limit in limits:
        ids = selected_generators(n, limit).intersection(active)
        require(len(ids) > 0, f"No active generators selected for {limit['id']}")
        g = n.generators.loc[ids]
        required = g.p_nom.where(~g.p_nom_extendable, g.p_nom_min)
        require(
            np.isfinite(required).all() and required.ge(0).all(),
            "Invalid capacity minima",
        )
        maximum = float(limit["maximum_MW"])
        require(np.isfinite(maximum) and maximum >= 0, "Invalid capacity ceiling")
        require(
            required.sum() <= maximum + 1e-6,
            f"{limit['id']}: existing/minimum capacity {required.sum():.3f} exceeds cap {maximum:.3f}; inspect retirements",
        )
        row = dict(
            limit,
            minimum_MW=float(required.sum()),
            assets=len(ids),
            extendable_assets=int(g.p_nom_extendable.sum()),
        )
        if solved:
            total = float(g.p_nom.where(~g.p_nom_extendable, g.p_nom_opt).sum())
            require(
                total <= maximum + max(0.01, maximum * 1e-7),
                f"{limit['id']}: solved capacity {total} exceeds {maximum}",
            )
            row.update(solved_MW=total, headroom_MW=maximum - total)
        rows.append(row)
    return rows


def add_capacity_constraints(n, settings):
    """One sum per region/technology, including every resource bin and existing asset."""
    limits = settings["capacity_limits"]
    capacity_report(n, limits)
    variable = n.model["Generator-p_nom"]
    dimension = variable.dims[0]
    names = variable.coords[dimension].to_index()
    for limit in limits:
        ids = selected_generators(n, limit).intersection(
            n.components["Generator"].active_assets
        )
        g = n.generators.loc[ids]
        fixed = float(g.loc[~g.p_nom_extendable, "p_nom"].sum())
        extendable = ids.intersection(names)
        if len(extendable):
            n.model.add_constraints(
                variable.sel({dimension: extendable}).sum()
                <= limit["maximum_MW"] - fixed,
                name=f"BEI-welfare-cap-{limit['id']}",
            )


def selected_loads(n, config):
    electric = n.loads.bus.map(n.buses.carrier).isin(["AC", "low voltage"])
    country = countries(n, n.loads.bus)
    return n.loads.index[electric & country.isin(config["demand"]["countries"])]


def demand_panel(n, ids, nonnegative=True):
    panel = pd.DataFrame(
        np.broadcast_to(
            n.loads.loc[ids, "p_set"].to_numpy(), (len(n.snapshots), len(ids))
        ).copy(),
        index=n.snapshots,
        columns=ids,
    )
    columns = ids.intersection(n.loads_t.p_set.columns)
    panel.loc[:, columns] = n.loads_t.p_set.loc[n.snapshots, columns]
    require(np.isfinite(panel).all().all(), "Non-finite demand")
    if nonnegative:
        require(
            panel.ge(0).all().all(), "Selected electricity demand must be non-negative"
        )
    return panel


def panel_digest(panel):
    # NetCDF can change datetime storage units without changing instants.
    panel = panel.astype("float64").copy()
    panel.index = pd.DatetimeIndex(panel.index).as_unit("ns")
    return hashlib.sha256(
        pd.util.hash_pandas_object(panel, index=True).values.tobytes()
        + json.dumps(list(panel.columns)).encode()
    ).hexdigest()


def apply(n, config, scenario):
    require(
        AUDIT_KEY not in n.meta,
        "Sensitivity already applied; use the original composed input",
    )
    limits = capacity_limits(config, scenario)
    report = capacity_report(n, limits)
    ids = selected_loads(n, config)
    require(len(ids) > 0, "No DK/DE fixed electricity demand selected")
    original = demand_panel(n, ids)
    factor = scenario["demand_factor"]
    require(np.isfinite(factor) and factor > 0, "Positive demand factor required")
    n.loads.loc[ids, "p_set"] *= factor
    columns = ids.intersection(n.loads_t.p_set.columns)
    n.loads_t.p_set.loc[:, columns] *= factor
    modified = demand_panel(n, ids)
    require(
        np.allclose(modified, original * factor, rtol=1e-12, atol=1e-10),
        "Demand scaling did not produce the requested factor",
    )
    weights = n.snapshot_weightings.generators
    country = countries(n, n.loads.loc[ids, "bus"])
    demand = {
        c: dict(
            before_MWh=float(
                original.loc[:, country.eq(c)].mul(weights, axis=0).sum().sum()
            ),
            after_MWh=float(
                modified.loc[:, country.eq(c)].mul(weights, axis=0).sum().sum()
            ),
        )
        for c in config["demand"]["countries"]
    }
    audit = dict(
        scenario=scenario,
        capacity_limits=report,
        demand=demand,
        demand_loads=list(ids),
        demand_sha256=panel_digest(modified),
        note="Fixed AC/low-voltage electricity Loads scaled; heat, transport, electrolysis and other endogenous conversion demand unchanged",
    )
    n.meta[AUDIT_KEY] = audit
    return audit


def validate(n, config, scenario, solved=False):
    audit = n.meta.get(AUDIT_KEY, {})
    require(audit.get("scenario") == scenario, "Wrong sensitivity scenario")
    require(
        panel_digest(demand_panel(n, selected_loads(n, config)))
        == audit.get("demand_sha256"),
        "Demand changed after preparation",
    )
    return capacity_report(n, capacity_limits(config, scenario), solved=solved)
