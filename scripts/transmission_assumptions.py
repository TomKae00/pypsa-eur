# SPDX-License-Identifier: MIT
"""
Dated operating MW constraints, separate from physical AC MVA ratings.

This is a zonal NTC approximation, not a reconstructed flow-based market.
JAO MaxBEX diagnostics cannot be installed as simultaneous bilateral limits.
"""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

AUDIT_KEY = "transmission_assumptions_audit"
ALIASES = {
    "SE01": "SE1",
    "SE02": "SE2",
    "SE03": "SE3",
    "SE04": "SE4",
    "DKW1": "DK1",
    "DKE1": "DK2",
    "NOS1": "NO1",
    "NOS2": "NO2",
    "NOM1": "NO3",
    "NON1": "NO4",
    "NOS5": "NO5",
    "FI00": "FI",
    "EE00": "EE",
    "LV00": "LV",
    "LT00": "LT",
    "DE00": "DE-LU",
    "LUG1": "DE-LU",
    "IE00": "IE/NI",
    "GBNI": "IE/NI",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def paths(settings, year):
    if not settings.get("enabled", False):
        return {}
    return {
        key: str(value).format(year=int(year), horizon=int(year))
        for key in ("limits", "hourly_limits", "joint_limits", "asset_rules")
        if (value := settings.get(key))
    }


def fingerprint(settings, year):
    hashes = {
        key: hashlib.sha256(Path(p).read_bytes()).hexdigest()
        for key, p in paths(settings, year).items()
    }
    payload = {"version": 1, "year": int(year), "settings": settings, "sha256": hashes}
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode()
    ).hexdigest(), hashes


def zones(n):
    return pd.Series(n.buses.index, index=n.buses.index).replace(ALIASES)


def dense(n, component, attr, names):
    frame = getattr(n, component.lower() + "s")
    result = pd.DataFrame(
        np.broadcast_to(
            frame.loc[names, attr].to_numpy(), (len(n.snapshots), len(names))
        ).copy(),
        index=n.snapshots,
        columns=names,
    )
    dynamic = getattr(getattr(n, component.lower() + "s_t"), attr)
    shared = dynamic.columns.intersection(names)
    if len(shared):
        result.loc[:, shared] = dynamic.reindex(n.snapshots)[shared].to_numpy()
    require(np.isfinite(result.to_numpy()).all(), f"Incomplete {component} {attr}")
    return result


def branches(n, zone0, zone1):
    """Coefficients measure exports at the sending zone's boundary, including DC losses."""
    z = zones(n)
    result = []
    for component, frame in [("Line", n.lines), ("Link", n.links)]:
        z0, z1 = frame.bus0.map(z), frame.bus1.map(z)
        selected = frame.index[
            ((z0 == zone0) & (z1 == zone1)) | ((z0 == zone1) & (z1 == zone0))
        ]
        if component == "Link":
            selected = selected[
                frame.loc[selected, "carrier"].isin(["AC", "DC", "B2B"])
            ]
        for name in selected:
            row = frame.loc[name]
            nominal = "s_nom" if component == "Line" else "p_nom"
            if component == "Line":
                coefficients = [1.0, -1.0]
                envelope = [
                    float(dense(n, component, "s_max_pu", [name])[name].max())
                ] * 2
            else:
                require(
                    all(not row.get(f"bus{i}", "") for i in range(2, 5)),
                    f"Multiterminal transfer needs explicit coefficients: {name}",
                )
                efficiency = dense(n, component, "efficiency", [name])[name]
                require(
                    efficiency.eq(efficiency.iloc[0]).all(),
                    f"Time-varying DC efficiency needs explicit coefficients: {name}",
                )
                eta = float(efficiency.iloc[0])
                require(0 < eta <= 1, f"Invalid efficiency: {name}")
                coefficients = [1.0, -eta]
                envelope = [
                    max(
                        0.0, float(dense(n, component, "p_max_pu", [name])[name].max())
                    ),
                    max(
                        0.0, -float(dense(n, component, "p_min_pu", [name])[name].min())
                    )
                    * eta,
                ]
            if z0[name] != zone0:
                coefficients.reverse()
                envelope.reverse()
            baseline = float(
                row[nominal + "_min"] if row[nominal + "_extendable"] else row[nominal]
            )
            require(
                np.isfinite(baseline) and baseline >= 0, f"Invalid capacity: {name}"
            )
            result.append(
                {
                    "component": component,
                    "name": str(name),
                    "bus0": str(row.bus0),
                    "bus1": str(row.bus1),
                    "baseline_MW": baseline,
                    "flow": coefficients,
                    "envelope": envelope,
                }
            )
    return result


def _utc(index):
    return pd.DatetimeIndex(pd.to_datetime(index, utc=True))


def apply_asset_rules(n, settings, year):
    """Apply before sector losses split bidirectional converters; never apply twice."""
    if not settings.get("enabled", False) or not settings.get("asset_rules"):
        return []
    path = Path(paths(settings, year)["asset_rules"])
    source_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    previous = n.meta.get("transmission_asset_rules")
    if previous:
        require(
            previous.get("year") == int(year) and previous.get("sha256") == source_hash,
            "Asset rules changed on a prepared network; rebuild",
        )
        return previous["assets"]
    table = pd.read_csv(path).fillna("")
    table = table[
        (table.year_min.astype(int) <= year) & (table.year_max.astype(int) >= year)
    ]
    audit = []
    sns = _utc(n.snapshots)
    for r in table.to_dict("records"):
        names = n.links.index[
            n.links.index.to_series()
            .str.replace("-reversed$", "", regex=True)
            .eq(r["asset"])
        ]
        if not len(names):
            audit.append(
                {
                    "asset": r["asset"],
                    "not_present": True,
                    "basis": r.get("basis", ""),
                    "source": r.get("source", ""),
                    "note": "Asset absent after upstream filtering or aggregation; no converter invented",
                }
            )
        for name in names:
            old = float(n.links.at[name, "p_nom"])
            if r.get("nominal_MW") != "":
                value = float(r["nominal_MW"])
                require(
                    np.isfinite(value) and value >= 0,
                    f"Invalid converter rating: {name}",
                )
                if value == 0:
                    n.remove("Link", name)
                    audit.append(
                        {
                            "asset": str(name),
                            "removed": True,
                            "basis": r.get("basis", ""),
                            "source": r.get("source", ""),
                        }
                    )
                    continue
                n.links.at[name, "p_nom"] = value
                if n.links.at[name, "p_nom_extendable"]:
                    n.links.at[name, "p_nom_min"] = value
                    require(
                        n.links.at[name, "p_nom_max"] >= value,
                        f"Inconsistent capacity maximum: {name}",
                    )
                elif n.links.at[name, "p_nom_max"] == old:
                    n.links.at[name, "p_nom_max"] = value
            nominal = float(n.links.at[name, "p_nom"])
            multiplier = pd.Series(1.0, index=n.snapshots)
            if r.get("available_from_utc"):
                multiplier.loc[sns < pd.Timestamp(r["available_from_utc"])] = 0.0
            if r.get("operating_MW") != "":
                require(
                    nominal > 0 and float(r["operating_MW"]) >= 0,
                    f"Invalid operating rating: {name}",
                )
                multiplier *= min(1.0, float(r["operating_MW"]) / nominal)
            for attr in ["p_min_pu", "p_max_pu"]:
                getattr(n.links_t, attr)[name] = (
                    dense(n, "Link", attr, [name])[name] * multiplier
                )
            audit.append(
                {
                    "asset": str(name),
                    "original_nominal_MW": old,
                    "nominal_MW": nominal,
                    "basis": r.get("basis", ""),
                    "source": r.get("source", ""),
                    "available_from_utc": r.get("available_from_utc", ""),
                }
            )
    n.meta["transmission_asset_rules"] = {
        "year": int(year),
        "sha256": source_hash,
        "assets": audit,
    }
    return audit


def install(n, settings, year, *, ac_losses=0):
    if not settings.get("enabled", False):
        return None
    require(
        ac_losses == 0,
        "Operating corridor approximation currently requires lossless AC",
    )
    require(
        n.investment_periods.empty,
        "Operating bounds currently support a single horizon",
    )
    digest, hashes = fingerprint(settings, year)
    policy = settings.get("expansion_policy", "fixed")
    require(policy in {"fixed", "proportional"}, "Unknown transfer expansion policy")
    input_paths = paths(settings, year)
    table = pd.read_csv(input_paths["limits"])
    table = table[table.planning_year.astype(int).eq(int(year))]
    require(not table.empty, f"No operating limits for planning year {year}")
    require(
        not table.apply(lambda r: tuple(sorted([r.zone0, r.zone1])), axis=1)
        .duplicated()
        .any(),
        "Duplicate operating corridor",
    )
    require(
        not table.get("basis", pd.Series("", index=table.index))
        .str.contains("MaxBEX", case=False)
        .any(),
        "MaxBEX is a diagnostic, not a simultaneous NTC table",
    )
    hourly = None
    if "hourly_limits" in input_paths:
        hourly = pd.read_csv(input_paths["hourly_limits"])
        hourly["snapshot"] = pd.to_datetime(hourly.snapshot, utc=True)
        require(
            not hourly.duplicated(["snapshot", "zone0", "zone1"]).any(),
            "Duplicate hourly transfer limits",
        )
    records, skipped = [], []
    n_zones = set(zones(n))
    for r in table.to_dict("records"):
        if r["zone0"] not in n_zones or r["zone1"] not in n_zones:
            continue
        crossed = branches(n, r["zone0"], r["zone1"])
        if not crossed:
            skipped.append(
                {
                    "zone0": r["zone0"],
                    "zone1": r["zone1"],
                    "reason": "No direct model crossing; no asset invented",
                }
            )
            continue
        limits = [float(r["forward_MW"]), float(r["reverse_MW"])]
        require(
            np.isfinite(limits).all() and min(limits) >= 0,
            "Finite nonnegative operating MW required",
        )
        envelope = [
            sum(b["baseline_MW"] * b["envelope"][i] for b in crossed) for i in range(2)
        ]
        require(
            policy != "proportional" or min(envelope) > 0,
            f"Empty expansion baseline: {r['zone0']}–{r['zone1']}",
        )
        record = {
            "zone0": r["zone0"],
            "zone1": r["zone1"],
            "limits_MW": limits,
            "baseline_envelope_MW": envelope,
            "branches": crossed,
            "basis": r.get("basis", ""),
            "source": r.get("source", ""),
        }
        if hourly is not None:
            h = hourly[(hourly.zone0 == r["zone0"]) & (hourly.zone1 == r["zone1"])]
            if len(h):
                h = h.set_index("snapshot").reindex(_utc(n.snapshots))
                values = h[["forward_MW", "reverse_MW"]].to_numpy(float)
                require(
                    np.isfinite(values).all() and (values >= 0).all(),
                    "Hourly limits must cover all snapshots",
                )
                record["hourly_limits_MW"] = values.tolist()
        records.append(record)
    joint = []
    if "joint_limits" in input_paths:
        j = pd.read_csv(input_paths["joint_limits"])
        j = j[j.planning_year.astype(int).eq(int(year))]
        for name, group in j.groupby("constraint"):
            require(
                group.limit_MW.nunique() == 1
                and np.isfinite(group.limit_MW).all()
                and group.limit_MW.ge(0).all(),
                "Invalid joint RHS",
            )
            terms = []
            for r in group.to_dict("records"):
                crossed = branches(n, r["zone0"], r["zone1"])
                require(
                    np.isfinite(float(r["coefficient"])), "Invalid joint coefficient"
                )
                if crossed:
                    terms.append(
                        {"coefficient": float(r["coefficient"]), "branches": crossed}
                    )
            if terms:
                positive = str(
                    group.get("positive_only", pd.Series([False])).iloc[0]
                ).lower() in {"1", "true"}
                joint.append(
                    {
                        "name": str(name),
                        "limit_MW": float(group.limit_MW.iloc[0]),
                        "terms": terms,
                        "positive_only": positive,
                    }
                )
    audit = {
        "version": 1,
        "year": int(year),
        "fingerprint": digest,
        "input_sha256": hashes,
        "settings": settings,
        "expansion_policy": policy,
        "corridors": records,
        "joint_limits": joint,
        "skipped_references": skipped,
        "snapshots_sha256": hashlib.sha256(
            _utc(n.snapshots).asi8.tobytes()
        ).hexdigest(),
        "interpretation": "Additional zonal NTC approximation; not a reconstructed flow-based market domain",
    }
    n.meta[AUDIT_KEY] = audit
    return audit


def validate(n, settings=None):
    audit = n.meta.get(AUDIT_KEY)
    if not audit:
        require(
            not settings or not settings.get("enabled", False),
            "Missing transmission assumptions; rebuild composition",
        )
        return None
    actual = hashlib.sha256(_utc(n.snapshots).asi8.tobytes()).hexdigest()
    require(
        actual == audit["snapshots_sha256"], "Transmission snapshots changed; recompose"
    )
    if settings is not None:
        require(
            fingerprint(settings, audit["year"])[0] == audit["fingerprint"],
            "Transmission inputs changed; recompose",
        )
    for corridor in audit["corridors"]:
        now = branches(n, corridor["zone0"], corridor["zone1"])

        def signature(rows):
            return [
                (
                    b["component"],
                    b["name"],
                    b["bus0"],
                    b["bus1"],
                    b["flow"],
                    b["envelope"],
                )
                for b in rows
            ]

        require(
            signature(now) == signature(corridor["branches"]),
            "Operating corridor topology or losses changed; recompose",
        )
    return audit


def _expression(n, branch_list, direction):
    total = 0
    for b in branch_list:
        component, name = b["component"], b["name"]
        total += b["flow"][direction] * n.model[
            f"{component}-{'s' if component == 'Line' else 'p'}"
        ].sel(name=name, drop=True)
    return total


def add_constraints(n, snapshots, settings=None):
    audit = validate(n, settings)
    if audit is None:
        return
    import xarray as xr

    for c in audit["corridors"]:
        for i in range(2):
            lhs = _expression(n, c["branches"], i)
            if "hourly_limits_MW" in c:
                limits = pd.Series(
                    np.asarray(c["hourly_limits_MW"])[:, i], index=n.snapshots
                ).reindex(snapshots)
                rhs = xr.DataArray(
                    limits.to_numpy(), coords={"snapshot": snapshots}, dims="snapshot"
                )
            else:
                rhs = c["limits_MW"][i]
            if audit["expansion_policy"] == "proportional":
                added = 0
                for b in c["branches"]:
                    component, name = b["component"], b["name"]
                    frame = getattr(n, component.lower() + "s")
                    nominal = "s_nom" if component == "Line" else "p_nom"
                    capacity = (
                        n.model[f"{component}-{nominal}"].sel(name=name, drop=True)
                        if frame.at[name, nominal + "_extendable"]
                        else float(frame.at[name, nominal])
                    )
                    added += b["envelope"][i] * (capacity - b["baseline_MW"])
                lhs -= rhs / c["baseline_envelope_MW"][i] * added
            a, b = (c["zone0"], c["zone1"]) if i == 0 else (c["zone1"], c["zone0"])
            n.model.add_constraints(lhs <= rhs, name=f"Operating-transfer-{a}-to-{b}")
    for c in audit["joint_limits"]:
        if c.get("positive_only", False):
            term_index = pd.Index(range(len(c["terms"])), name="term")
            positive = n.model.add_variables(
                lower=0,
                coords=[snapshots, term_index],
                name="Operating-positive-" + c["name"],
            )
            for i, t in enumerate(c["terms"]):
                expression = t["coefficient"] * _expression(n, t["branches"], 0)
                n.model.add_constraints(
                    positive.sel(term=i, drop=True) >= expression,
                    name=f"Operating-positive-{c['name']}-{i}",
                )
            lhs = positive.sum("term")
        else:
            lhs = sum(
                t["coefficient"] * _expression(n, t["branches"], 0) for t in c["terms"]
            )
        n.model.add_constraints(
            lhs <= c["limit_MW"], name="Operating-joint-" + c["name"]
        )


def _flows(n, branches_, direction=0):
    flow = pd.Series(0.0, index=n.snapshots)
    for b in branches_:
        flow += (
            b["flow"][direction]
            * getattr(getattr(n, b["component"].lower() + "s_t"), "p0")[b["name"]]
        )
    return flow


def check_solution(n, settings=None, tolerance_MW=1e-3):
    audit = validate(n, settings)
    if audit is None:
        return []
    out = []
    for c in audit["corridors"]:
        for i in range(2):
            flow = _flows(n, c["branches"], i)
            limit = (
                pd.Series(np.asarray(c["hourly_limits_MW"])[:, i], index=n.snapshots)
                if "hourly_limits_MW" in c
                else pd.Series(c["limits_MW"][i], index=n.snapshots)
            )
            if audit["expansion_policy"] == "proportional":
                added = 0.0
                for b in c["branches"]:
                    frame = getattr(n, b["component"].lower() + "s")
                    nom = "s_nom" if b["component"] == "Line" else "p_nom"
                    value = float(
                        frame.at[b["name"], nom + "_opt"]
                        if frame.at[b["name"], nom + "_extendable"]
                        else frame.at[b["name"], nom]
                    )
                    added += b["envelope"][i] * (value - b["baseline_MW"])
                limit *= 1 + added / c["baseline_envelope_MW"][i]
            violation = float((flow - limit).max())
            require(
                np.isfinite(flow).all() and violation <= tolerance_MW,
                f"Operating limit violated on {c['zone0']}–{c['zone1']}",
            )
            out.append(
                {
                    "zone0": c["zone0"],
                    "zone1": c["zone1"],
                    "direction": i,
                    "max_flow_MW": float(flow.max()),
                    "max_violation_MW": max(0.0, violation),
                }
            )
    for c in audit["joint_limits"]:
        flow = pd.Series(0.0, index=n.snapshots)
        for t in c["terms"]:
            value = t["coefficient"] * _flows(n, t["branches"])
            flow += value.clip(lower=0) if c.get("positive_only", False) else value
        require(
            np.isfinite(flow).all()
            and float((flow - c["limit_MW"]).max()) <= tolerance_MW,
            "Joint operating limit violated: " + c["name"],
        )
    return out
