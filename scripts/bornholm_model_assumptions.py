# SPDX-License-Identifier: MIT
"""Opt-in, auditable assumptions for the common Bornholm study base."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

AUDIT_KEY = "bornholm_assumptions_audit"
VERSION = 1
ZONE_ALIASES = {
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
    "DE00": "DE-LU",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def input_paths(settings):
    """Return explicit workflow dependencies, including disabled-mode support."""
    if not settings.get("enabled", False):
        return {}
    return {
        "onshore_calibration": settings["onshore_calibration"],
        "zonal_transfer_limits": settings["transfer_limits"],
    }


def fingerprint(settings):
    """Pin assumption settings and table contents, independent of output paths."""
    paths = input_paths(settings)
    hashes = {
        key: hashlib.sha256(Path(path).read_bytes()).hexdigest()
        for key, path in paths.items()
    }
    payload = {"version": VERSION, "settings": settings, "table_sha256": hashes}
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    return digest, hashes


def dense(n, component, attribute, names):
    """Combine PyPSA's static defaults with explicit time-dependent columns."""
    frame = getattr(n, component)
    result = pd.DataFrame(
        np.broadcast_to(
            frame.loc[names, attribute].to_numpy(), (len(n.snapshots), len(names))
        ).copy(),
        index=n.snapshots,
        columns=names,
    )
    dynamic = getattr(getattr(n, component + "_t"), attribute)
    columns = dynamic.columns.intersection(names)
    result.loc[:, columns] = dynamic.reindex(n.snapshots)[columns]
    return result


def generator_countries(n, names):
    buses = n.generators.loc[names, "bus"]
    countries = buses.map(n.buses.country)
    valid = countries.fillna("").str.fullmatch("[A-Z]{2}")
    parents = buses.map(n.buses.location).map(n.buses.country)
    countries = countries.where(valid, parents)
    require(
        countries.fillna("").str.fullmatch("[A-Z]{2}").all(),
        "Onshore generator country is unresolved; inspect bus country/location",
    )
    return countries


def reference_cf(n, names, cf):
    weights = n.generators.loc[names, "p_nom"].where(
        ~n.generators.loc[names, "p_nom_extendable"],
        n.generators.loc[names, "p_nom_min"],
    )
    require(
        np.isfinite(weights).all() and weights.ge(0).all(),
        "Invalid required background capacities",
    )
    return float(cf.loc[names] @ weights / weights.sum()) if weights.sum() > 0 else None


def fix_coal_buffer(n, capacity):
    require(np.isfinite(capacity) and capacity > 0, "Positive finite coal MWh required")
    names = n.stores.index[n.stores.carrier.eq("coal")]
    require(len(names) == 1, "Expected the single audited coal fuel Store")
    name = names[0]
    store = n.stores.loc[name]
    require(
        store.e_cyclic
        and store.standing_loss == 0
        and store.e_min_pu == 0
        and store.e_max_pu == 1,
        "Coal inventory formulation differs from the audited fuel buffer",
    )
    for attr, expected in (("e_min_pu", 0), ("e_max_pu", 1), ("standing_loss", 0)):
        series = getattr(n.stores_t, attr)
        require(
            name not in series or series[name].eq(expected).all(),
            f"Unexpected time-varying coal {attr}",
        )
    supply = n.generators.loc[
        n.generators.bus.eq(store.bus) & n.generators.carrier.eq("coal")
    ]
    require(
        not supply.empty and (supply.p_nom_extendable | supply.p_nom.gt(0)).any(),
        "Coal buffer requires independent fuel supply",
    )
    n.stores.loc[name, ["e_nom", "e_nom_min", "e_nom_max", "e_initial"]] = [
        capacity,
        0,
        capacity,
        0,
    ]
    n.stores.loc[name, "e_nom_extendable"] = False
    return {"store": str(name), "fixed_MWh": float(capacity)}


def apply_resource_and_coal_assumptions(n, settings):
    """Scale hourly profiles once, before temporal aggregation, and fix coal."""
    if not settings.get("enabled", False):
        return
    require(n.investment_periods.empty, "Bornholm calibration supports overnight only")
    digest, hashes = fingerprint(settings)
    previous = n.meta.get(AUDIT_KEY)
    if previous is not None:
        require(
            previous.get("fingerprint") == digest,
            "Assumptions changed on an already calibrated network; recompose",
        )
        return previous

    table = pd.read_csv(settings["onshore_calibration"])
    require(
        {"country", "carrier", "factor"}.issubset(table),
        "Missing onshore calibration columns",
    )
    require(
        not table.duplicated(["country", "carrier"]).any(),
        "Duplicate country/carrier calibration",
    )
    require(
        table.country.str.fullmatch("[A-Z]{2}").all()
        and table.carrier.eq("onwind").all(),
        "Only explicit country onwind rows supported",
    )
    require(
        np.isfinite(table.factor).all() and table.factor.gt(0).all(),
        "Calibration factors must be positive and finite",
    )
    require(
        not set(table.country) & set(settings.get("excluded_countries", [])),
        "An excluded country appears in the calibration table",
    )

    names = n.generators.index[n.generators.carrier.eq("onwind")]
    require(len(names) > 0, "Missing onshore generators")
    countries = generator_countries(n, names)
    # A subset model may omit a table country; retain the decision in provenance.
    factors = countries.map(table.set_index("country").factor).fillna(1.0)
    profiles = dense(n, "generators", "p_max_pu", names)
    require(
        np.isfinite(profiles).all().all()
        and profiles.ge(-1e-8).all().all()
        and profiles.le(1 + 1e-8).all().all(),
        "Invalid input onshore availability",
    )
    weights = n.snapshot_weightings.generators
    require(
        np.isfinite(weights).all() and weights.gt(0).all(),
        "Positive generator snapshot weights required",
    )
    scaled = profiles.mul(factors, axis=1)
    corrected = scaled.clip(0, 1)
    before = profiles.mul(weights, axis=0).sum() / weights.sum()
    after = corrected.mul(weights, axis=0).sum() / weights.sum()
    assets = {
        str(name): {
            "country": str(countries[name]),
            "factor": float(factors[name]),
            "before_cf": float(before[name]),
            "after_cf": float(after[name]),
            "clipped_fraction": float((scaled[name] > 1).mean()),
        }
        for name in names
    }
    regions = {}
    for country in sorted(countries.unique()):
        subset = countries.index[countries.eq(country)]
        regions[str(country)] = {
            "reference_before_cf": reference_cf(n, subset, before),
            "reference_after_cf": reference_cf(n, subset, after),
            "class_min_cf": float(after[subset].min()),
            "class_max_cf": float(after[subset].max()),
            "factor": float(factors[subset].iloc[0]),
        }

    # Validate coal before mutating the availability matrix.
    coal = fix_coal_buffer(n, float(settings["coal_store_MWh"]))
    changed = names[factors.ne(1)]
    for name in changed:
        n.generators_t.p_max_pu[name] = corrected[name]
    audit = {
        "version": VERSION,
        "fingerprint": digest,
        "table_sha256": hashes,
        "resources": {"assets": assets, "countries": regions},
        "coal": coal,
    }
    n.meta[AUDIT_KEY] = audit
    return audit


def corridor_branches(n, zone0, zone1):
    """Describe each direct AC/DC crossing once, with terminal MW coefficients."""
    records = []
    zones = pd.Series(n.buses.index, index=n.buses.index).replace(ZONE_ALIASES)
    for component, frame, nominal, attrs in (
        ("Line", n.lines, "s_nom", ["s_max_pu"]),
        (
            "Link",
            n.links[n.links.carrier.isin(["DC", "AC"])],
            "p_nom",
            ["p_min_pu", "p_max_pu", "efficiency"],
        ),
    ):
        z0, z1 = frame.bus0.map(zones), frame.bus1.map(zones)
        selected = frame.index[
            ((z0 == zone0) & (z1 == zone1)) | ((z0 == zone1) & (z1 == zone0))
        ]
        for attr in attrs:
            dynamic = getattr(getattr(n, component.lower() + "s_t"), attr)
            require(
                dynamic.columns.intersection(selected).empty,
                f"Dynamic {component} {attr} is unsupported on calibrated corridors",
            )
        for name in selected:
            row = frame.loc[name]
            forward = z0[name] == zone0
            if component == "Line":
                flow = [1.0, -1.0] if forward else [-1.0, 1.0]
                envelope = [float(row.s_max_pu)] * 2
            else:
                require(
                    all(not row.get(f"bus{i}", "") for i in range(2, 5)),
                    f"Multiterminal crossing needs a separate audit: {name}",
                )
                require(
                    np.isfinite(row.efficiency) and 0 < row.efficiency <= 1,
                    f"Invalid DC transmission efficiency: {name}",
                )
                outgoing = max(float(row.p_max_pu), 0)
                incoming = max(-float(row.p_min_pu), 0) * float(row.efficiency)
                flow = [1.0, -float(row.efficiency)]
                envelope = [outgoing, incoming]
                if not forward:
                    flow.reverse()
                    envelope.reverse()
            baseline = float(
                row[nominal + "_min"] if row[nominal + "_extendable"] else row[nominal]
            )
            require(
                np.isfinite(baseline)
                and baseline >= 0
                and np.isfinite(envelope).all()
                and min(envelope) >= 0,
                f"Invalid corridor envelope: {name}",
            )
            records.append(
                {
                    "component": component,
                    "name": str(name),
                    "bus0": str(row.bus0),
                    "bus1": str(row.bus1),
                    "baseline_MW": baseline,
                    "flow": flow,
                    "envelope": envelope,
                }
            )
    require(records, f"No physical crossing for {zone0}--{zone1}")
    return records


def install_transfer_assumptions(n, settings, planning_year, *, ac_losses=0):
    """Pin corridor membership after native transmission and project preparation."""
    if not settings.get("enabled", False):
        return
    require(ac_losses == 0, "Calibrated corridor formulation requires lossless AC")
    audit = n.meta.get(AUDIT_KEY, {})
    require(
        audit.get("fingerprint") == fingerprint(settings)[0],
        "Apply resource/coal assumptions before installing transfer limits",
    )
    require(
        settings.get("expansion_policy") in ("fixed", "proportional"),
        "Unknown corridor expansion policy",
    )
    table = pd.read_csv(settings["transfer_limits"])
    columns = {"planning_year", "zone0", "zone1", "forward_MW", "reverse_MW"}
    require(columns.issubset(table), "Missing zonal transfer-limit columns")
    table = table.loc[table.planning_year.eq(int(planning_year))]
    require(
        not table.empty, f"No transfer assumptions for planning year {planning_year}"
    )
    pairs = table.apply(lambda r: tuple(sorted([r.zone0, r.zone1])), axis=1)
    require(
        not pairs.duplicated().any() and table.zone0.ne(table.zone1).all(),
        "Duplicate or invalid corridor",
    )
    require(
        np.isfinite(table[["forward_MW", "reverse_MW"]]).all().all()
        and table[["forward_MW", "reverse_MW"]].gt(0).all().all(),
        "Positive finite directional MW limits required",
    )
    corridors = []
    for row in table.itertuples(index=False):
        branches = corridor_branches(n, row.zone0, row.zone1)
        envelope = [
            sum(b["baseline_MW"] * b["envelope"][i] for b in branches) for i in range(2)
        ]
        require(
            min(envelope) > 0,
            f"Missing directional physical envelope: {row.zone0}--{row.zone1}",
        )
        corridors.append(
            {
                "zone0": row.zone0,
                "zone1": row.zone1,
                "limits_MW": [float(row.forward_MW), float(row.reverse_MW)],
                "baseline_envelope_MW": envelope,
                "branches": branches,
            }
        )
    audit["transfers"] = {
        "planning_year": int(planning_year),
        "expansion_policy": settings["expansion_policy"],
        "corridors": corridors,
    }
    return audit["transfers"]


def validate_assumptions(n, settings):
    """Reject stale composition rather than silently calibrating at solve time."""
    if not settings.get("enabled", False):
        require(
            AUDIT_KEY not in n.meta,
            "Calibrated network requires its assumption configuration",
        )
        return
    audit = n.meta.get(AUDIT_KEY, {})
    require(
        audit.get("fingerprint") == fingerprint(settings)[0]
        and "resources" in audit
        and "transfers" in audit,
        "Missing or stale Bornholm assumptions; rebuild composed network",
    )
    store = n.stores.loc[audit["coal"]["store"]]
    require(
        not store.e_nom_extendable
        and np.isclose(store.e_nom, settings["coal_store_MWh"])
        and store.e_cyclic
        and store.standing_loss == 0
        and store.e_initial == 0
        and store.e_min_pu == 0
        and store.e_max_pu == 1,
        "Coal buffer no longer matches the composed assumptions",
    )
    for corridor in audit["transfers"]["corridors"]:
        current = corridor_branches(n, corridor["zone0"], corridor["zone1"])

        # Capacity may change in native iterative reinforcement; topology,
        # ratings per unit and terminal efficiencies may not change silently.
        def signatures(records):
            return [{k: v for k, v in r.items() if k != "baseline_MW"} for r in records]

        require(
            signatures(current) == signatures(corridor["branches"]),
            "Calibrated corridor topology or coefficients changed; recompose",
        )
    return audit


def add_zonal_transfer_constraints(n, snapshots, settings):
    """Supplement native KVL and thermal constraints with directional ceilings."""
    audit = validate_assumptions(n, settings)
    if audit is None:
        return
    transfers = audit["transfers"]
    for corridor in transfers["corridors"]:
        for direction in range(2):
            flow, added_envelope = 0, 0
            for branch in corridor["branches"]:
                component, name = branch["component"], branch["name"]
                flow_var = "s" if component == "Line" else "p"
                flow += branch["flow"][direction] * n.model[
                    f"{component}-{flow_var}"
                ].sel(
                    name=name,
                    snapshot=snapshots,
                )
                nominal = "s_nom" if component == "Line" else "p_nom"
                frame = getattr(n, component.lower() + "s")
                if frame.at[name, nominal + "_extendable"]:
                    capacity = n.model[f"{component}-{nominal}"].sel(name=name)
                else:
                    capacity = float(frame.at[name, nominal])
                added_envelope += branch["envelope"][direction] * (
                    capacity - branch["baseline_MW"]
                )
            limit = corridor["limits_MW"][direction]
            if transfers["expansion_policy"] == "proportional":
                flow -= (
                    limit / corridor["baseline_envelope_MW"][direction] * added_envelope
                )
            source, target = (corridor["zone0"], corridor["zone1"])
            if direction == 1:
                source, target = target, source
            n.model.add_constraints(
                flow <= limit, name=f"Bornholm-transfer-{source}-to-{target}"
            )


def resource_report(n, *, solved=False):
    """Report reference and returned fleet CF without using dispatch to calibrate."""
    names = n.generators.index[n.generators.carrier.eq("onwind")]
    countries = generator_countries(n, names)
    profiles = dense(n, "generators", "p_max_pu", names)
    weights = n.snapshot_weightings.generators
    cf = profiles.mul(weights, axis=0).sum() / weights.sum()
    report = {}
    for country in sorted(countries.unique()):
        subset = names[countries.eq(country)]
        record = {
            "reference_available_cf": reference_cf(n, subset, cf),
            "class_min_cf": float(cf[subset].min()),
            "class_max_cf": float(cf[subset].max()),
        }
        if solved:
            capacity = n.generators.loc[subset, "p_nom_opt"]
            total = float(capacity.sum())
            record.update(
                capacity_GW=total / 1000,
                fleet_available_cf=float(cf[subset] @ capacity / total)
                if total > 0
                else None,
                fleet_realised_cf=float(
                    n.generators_t.p[subset].mul(weights, axis=0).sum().sum()
                    / (total * weights.sum())
                )
                if total > 0
                else None,
            )
        report[str(country)] = record
    return report


def check_solved_assumptions(n, settings, *, tolerance_MW=1e-3):
    """Check returned inventory and combined directional flows against the model."""
    audit = validate_assumptions(n, settings)
    require(audit is not None, "Solved checks require enabled assumptions")
    coal = audit["coal"]
    inventory = n.stores_t.e[coal["store"]]
    require(
        np.isfinite(inventory).all()
        and inventory.min() >= -tolerance_MW
        and inventory.max() <= coal["fixed_MWh"] + tolerance_MW,
        "Returned coal inventory exceeds fixed buffer bounds",
    )
    results = []
    for corridor in audit["transfers"]["corridors"]:
        for direction in range(2):
            flow = pd.Series(0.0, index=n.snapshots)
            added_envelope = 0.0
            for branch in corridor["branches"]:
                component, name = branch["component"], branch["name"]
                frame = getattr(n, component.lower() + "s")
                dynamic = getattr(n, component.lower() + "s_t")
                flow += branch["flow"][direction] * dynamic.p0[name]
                nominal = "s_nom" if component == "Line" else "p_nom"
                # Iterative AC reinforcement can fix the final nominal capacity.
                capacity = float(frame.at[name, nominal + "_opt"])
                added_envelope += branch["envelope"][direction] * (
                    capacity - branch["baseline_MW"]
                )
            limit = corridor["limits_MW"][direction]
            if audit["transfers"]["expansion_policy"] == "proportional":
                limit *= (
                    1 + added_envelope / corridor["baseline_envelope_MW"][direction]
                )
            require(
                np.isfinite(flow).all() and np.isfinite(limit),
                "Non-finite solved corridor flow or limit",
            )
            excess = float((flow - limit).max())
            require(
                excess <= tolerance_MW,
                f"Returned flow exceeds {corridor['zone0']}--{corridor['zone1']} limit",
            )
            source, target = corridor["zone0"], corridor["zone1"]
            if direction == 1:
                source, target = target, source
            results.append(
                {
                    "direction": f"{source}->{target}",
                    "limit_MW": float(limit),
                    "maximum_net_export_MW": float(flow.max()),
                    "maximum_excess_MW": max(0.0, excess),
                    "binding_snapshot_count": int(
                        (flow - limit).abs().le(tolerance_MW).sum()
                    ),
                }
            )
    return {
        "coal_min_MWh": float(inventory.min()),
        "coal_max_MWh": float(inventory.max()),
        "corridors": results,
        "onshore": resource_report(n, solved=True),
    }
