#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only audit of a deterministic, single-period PyPSA candidate.

No optimization model is built and no solver is called. Tested with PyPSA 1.2.4.
Run: python -m scripts.audit_bornholm_candidate CANDIDATE_DIRECTORY --output audit_output
Supports accuracy candidates and original study networks/solved_2035.nc
with summary.json. Configuration comes from effective_config.yaml when present,
otherwise from the configuration embedded in the saved network.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
import yaml
from pypsa.optimization.optimize import lookup

NOMINAL = {
    "Generator": "p_nom",
    "Link": "p_nom",
    "Line": "s_nom",
    "Transformer": "s_nom",
    "Store": "e_nom",
    "StorageUnit": "p_nom",
}
MODEL_VARIABLES = {
    "Generator": {"p": "p"},
    "Link": {"p": "p0"},
    "Line": {"s": "p0"},
    "Transformer": {"s": "p0"},
    "Store": {"e": "e"},
    "StorageUnit": {
        "p_dispatch": "p_dispatch",
        "p_store": "p_store",
        "state_of_charge": "state_of_charge",
    },
}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024**2), b""):
            h.update(block)
    return h.hexdigest()


def safe_product(a, b):
    """Native PyPSA convention: a zero cost/bound times infinity is zero."""
    with np.errstate(invalid="ignore"):
        out = np.asarray(a) * np.asarray(b)
    return np.where((np.asarray(a) == 0) | (np.asarray(b) == 0), 0, out)


def storage_residual(energy, power, hours, standing_loss, cyclic, initial):
    previous = np.roll(energy, 1, axis=0).copy()
    previous[0, ~cyclic] = initial[~cyclic]
    # Initial values enter the RHS without a standing-loss multiplier in PyPSA.
    retention = (1 - standing_loss) ** hours[:, None]
    retention[0, ~cyclic] = 1
    return energy - retention * previous + hours[:, None] * power


class Audit:
    def __init__(self, n):
        self.n = n
        self.sns = n.snapshots
        self.metrics = []
        self.worst = []
        self.scope = []
        self.ids = {k: n.components[k].active_assets for k in [*NOMINAL, "Load"]}

    def dense(self, component, attr, ids=None):
        c = self.n.components[component]
        ids = self.ids[component] if ids is None else pd.Index(ids)
        base = c.static.loc[ids, attr].to_numpy(dtype=float)
        out = np.broadcast_to(base, (len(self.sns), len(ids))).copy()
        dynamic = c.dynamic.get(attr)
        if dynamic is not None and not dynamic.empty:
            cols = ids.intersection(dynamic.columns)
            out[:, ids.get_indexer(cols)] = dynamic.loc[self.sns, cols].to_numpy()
        return out

    def result(self, component, attr, ids=None):
        ids = self.ids[component] if ids is None else pd.Index(ids)
        c = self.n.components[component]
        frame = c.dynamic[attr]
        missing = ids.difference(frame.columns)
        if not missing.empty:
            # PyPSA omits an output table when every value equals its schema
            # default (usually zero). Restore that documented storage convention.
            default = float(c.defaults.at[attr, "default"])
            if not np.isfinite(default):
                raise ValueError(f"Missing non-default result: {component}.{attr}")
            frame = frame.reindex(index=self.sns, columns=ids, fill_value=default)
            note = (
                f"Restored omitted {component}.{attr} columns using the PyPSA "
                f"schema default {default:g}."
            )
            if note not in self.scope:
                self.scope.append(note)
        out = frame.loc[self.sns, ids].to_numpy()
        if not np.isfinite(out).all():
            raise ValueError(f"Non-finite or missing result: {component}.{attr}")
        return out

    def record(
        self,
        name,
        residual,
        columns,
        unit,
        *,
        component="",
        timestamps=True,
        inequality=False,
        scale=None,
    ):
        r = np.asarray(residual, dtype=float)
        columns = pd.Index(columns)
        if r.size == 0:
            return
        if np.isnan(r).any():
            raise ValueError(f"NaN in evaluated check {name}")
        # A negative infinite residual is a satisfied unbounded inequality.
        error = np.maximum(r, 0) if inequality else np.abs(r)
        if not np.isfinite(error).all():
            raise ValueError(f"Non-finite violation in check {name}")
        pos = int(error.argmax())
        if r.ndim == 2:
            ti, ci = np.unravel_index(pos, r.shape)
            timestamp = str(self.sns[ti]) if timestamps else ""
        else:
            ci, ti, timestamp = pos, 0, ""
        row = {
            "check": name,
            "component": component,
            "unit": unit,
            "equations_or_bounds": int(r.size),
            "max_violation": float(error.flat[pos]),
            "asset_or_bus": str(columns[ci]),
            "timestamp": timestamp,
            "signed_residual_at_max": float(r.flat[pos]),
            "count_above_1e-3": int((error > 1e-3).sum()),
            "count_above_1e-6": int((error > 1e-6).sum()),
        }
        if scale is not None:
            activity = np.broadcast_to(scale, r.shape)
            row["activity_at_max"] = float(activity.flat[pos])
            row["relative_at_max"] = float(
                error.flat[pos] / max(1, abs(activity.flat[pos]))
            )
        self.metrics.append(row)
        # Retain several large residuals, with their original units and location.
        count = min(10, error.size)
        top = np.argpartition(error.ravel(), -count)[-count:]
        for p in top[np.argsort(error.ravel()[top])[::-1]]:
            if error.flat[p] <= 1e-9:
                continue
            if r.ndim == 2:
                t, j = np.unravel_index(p, r.shape)
                date = str(self.sns[t]) if timestamps else ""
            else:
                j, date = p, ""
            self.worst.append(
                {
                    "check": name,
                    "component": component,
                    "asset_or_bus": str(columns[j]),
                    "timestamp": date,
                    "violation": float(error.flat[p]),
                    "unit": unit,
                }
            )

    def costs(self):
        rows = []
        for component, nominal in NOMINAL.items():
            c = self.n.components[component]
            ids = self.ids[component]
            if ids.empty:
                continue
            static = c.static.loc[ids]
            extendable = static[nominal + "_extendable"].to_numpy(bool)
            periodic = c.periodized_cost.sel(name=ids).to_numpy()
            # Only extendables enter the native investment objective.
            gross = safe_product(periodic, static[nominal + "_opt"].to_numpy())
            constant = safe_product(periodic, static[nominal].to_numpy())
            for j, asset in enumerate(ids):
                rows.append(
                    {
                        "component": component,
                        "asset": str(asset),
                        "carrier": str(static.loc[asset, "carrier"]),
                        "cost_type": "investment",
                        "extendable": bool(extendable[j]),
                        "periodized_unit_cost": float(periodic[j]),
                        "nominal": float(static.loc[asset, nominal]),
                        "nominal_opt": float(static.loc[asset, nominal + "_opt"]),
                        "gross_eur": float(gross[j]) if extendable[j] else 0,
                        "existing_constant_eur": float(constant[j])
                        if extendable[j]
                        else 0,
                        "native_objective_eur": float(gross[j] - constant[j])
                        if extendable[j]
                        else 0,
                    }
                )
            if static.get("committable", pd.Series(False, index=ids)).any():
                raise ValueError("Unit commitment costs need a separate audit")
            if (
                static.get("marginal_cost_quadratic", pd.Series(0, index=ids)) != 0
            ).any():
                raise ValueError("Quadratic costs need a separate audit")

        weighting = self.n.snapshot_weightings.objective.to_numpy()
        for cost_type in ["marginal_cost", "marginal_cost_storage", "spill_cost"]:
            for component, variable in lookup.query(cost_type).index:
                if component not in self.ids or self.ids[component].empty:
                    continue
                ids = self.ids[component]
                cost = self.dense(component, cost_type, ids)
                if not np.any(cost):
                    continue
                attribute = (
                    "p0" if component == "Link" and variable == "p" else variable
                )
                operation = self.result(component, attribute, ids)
                amounts = (operation * cost * weighting[:, None]).sum(axis=0)
                for j in np.flatnonzero(np.any(cost != 0, axis=0)):
                    rows.append(
                        {
                            "component": component,
                            "asset": str(ids[j]),
                            "carrier": self.n.components[component].static.loc[
                                ids[j], "carrier"
                            ],
                            "cost_type": cost_type,
                            "extendable": False,
                            "gross_eur": float(amounts[j]),
                            "existing_constant_eur": 0,
                            "native_objective_eur": float(amounts[j]),
                        }
                    )
        table = pd.DataFrame(rows)
        totals = table.groupby("cost_type")[
            ["gross_eur", "existing_constant_eur", "native_objective_eur"]
        ].sum()
        constant = float(totals.existing_constant_eur.sum())
        saved_constant = float(self.n.objective_constant)
        # include_objective_constant=False is also supported: saved constant is zero.
        use_constant = constant if saved_constant != 0 else 0
        if saved_constant != 0 and abs(constant - saved_constant) > 1:
            raise ValueError(
                f"Constant reconstruction mismatch: {constant} vs {saved_constant}"
            )
        reconstructed = float(table.gross_eur.sum() - use_constant)
        result = {
            "gross_investment_eur": float(totals.loc["investment", "gross_eur"]),
            "existing_capacity_constant_eur": saved_constant,
            "operating_eur": float(totals.drop("investment").gross_eur.sum()),
            "reconstructed_objective_eur": reconstructed,
            "saved_objective_eur": float(self.n.objective),
            "reconstruction_error_eur": reconstructed - float(self.n.objective),
        }
        return table, result

    def balance(self):
        n = self.n
        ids = n.buses.index
        residual = np.zeros((len(self.sns), len(ids)))
        activity = np.zeros_like(residual)

        def accumulate(values, bus):
            bus = pd.Index(bus)
            mapped = ids.get_indexer(bus)
            if (mapped < 0).any():
                raise ValueError(f"Undefined buses: {list(bus[mapped < 0])}")
            # Transpose permits indexed accumulation when multiple assets share a bus.
            np.add.at(residual.T, mapped, values.T)
            np.add.at(activity.T, mapped, np.abs(values.T))

        for component in ["Generator", "Store", "StorageUnit", "Load"]:
            ci = self.ids[component]
            static = n.components[component].static.loc[ci]
            values = self.result(component, "p") * static.sign.to_numpy()[None, :]
            accumulate(values, static.bus)
            if component == "Load":
                self.record(
                    "load_equals_p_set",
                    self.result(component, "p") - self.dense(component, "p_set"),
                    ci,
                    "bus-flow units",
                    component=component,
                )
            if component == "StorageUnit":
                self.record(
                    "storage_unit_net_dispatch",
                    self.result(component, "p")
                    - self.result(component, "p_dispatch")
                    + self.result(component, "p_store"),
                    ci,
                    "MW",
                    component=component,
                )
        for component in ["Link", "Line", "Transformer"]:
            c = n.components[component]
            ci = self.ids[component]
            for col in c.static.columns:
                if not col.startswith("bus") or not col[3:].isdigit():
                    continue
                valid = c.static.loc[ci, col].fillna("").ne("")
                subset = ci[valid]
                if subset.empty:
                    continue
                accumulate(
                    -self.result(component, "p" + col[3:], subset),
                    c.static.loc[subset, col],
                )
        self.record(
            "nodal_balance",
            residual,
            ids,
            "bus-flow units",
            component="Bus",
            scale=activity,
        )
        bus_table = pd.DataFrame(
            {
                "bus": ids,
                "carrier": n.buses.carrier,
                "unit": n.buses.unit,
                "max_balance_error": np.abs(residual).max(axis=0),
                "max_activity": activity.max(axis=0),
            }
        )
        grouped = bus_table.groupby(["carrier", "unit"], dropna=False).agg(
            max_balance_error=("max_balance_error", "max"), buses=("bus", "size")
        )
        return bus_table, grouped.reset_index()

    def limits(self):
        for component, nominal in NOMINAL.items():
            ids = self.ids[component]
            if ids.empty:
                continue
            c = self.n.components[component]
            static = c.static.loc[ids]
            ext = static[nominal + "_extendable"].to_numpy(bool)
            cap = np.where(ext, static[nominal + "_opt"], static[nominal])
            nominal_unit = (
                "MWh or carrier-energy units"
                if component == "Store"
                else "MW or carrier-flow units"
            )
            for side in ["min", "max"]:
                bound = static[nominal + "_" + side].to_numpy()[ext]
                actual = static[nominal + "_opt"].to_numpy()[ext]
                r = bound - actual if side == "min" else actual - bound
                self.record(
                    "nominal_" + side,
                    r,
                    ids[ext],
                    nominal_unit,
                    component=component,
                    timestamps=False,
                    inequality=True,
                )
            fixed = ~ext & np.isfinite(cap)
            self.record(
                "fixed_nominal_retained",
                static[nominal + "_opt"].to_numpy()[fixed] - cap[fixed],
                ids[fixed],
                nominal_unit,
                component=component,
                timestamps=False,
            )
            if nominal + "_set" in static:
                selected = ext & static[nominal + "_set"].notna().to_numpy()
                self.record(
                    "nominal_set",
                    static[nominal + "_opt"].to_numpy()[selected]
                    - static[nominal + "_set"].to_numpy()[selected],
                    ids[selected],
                    nominal_unit,
                    component=component,
                    timestamps=False,
                )

            for variable, saved in MODEL_VARIABLES[component].items():
                values = self.result(component, saved)
                # Use the same public bounds method as native PyPSA constraints.
                lower, upper = c.get_bounds_pu(attr=variable)

                def expand(bound, ids=ids, shape=values.shape):
                    bound = bound.sel(name=ids)
                    if "snapshot" in bound.dims:
                        return (
                            bound.sel(snapshot=self.sns)
                            .transpose("snapshot", "name")
                            .to_numpy()
                        )
                    return np.broadcast_to(bound.to_numpy(), shape)

                low = safe_product(expand(lower), cap[None, :])
                high = safe_product(expand(upper), cap[None, :])
                unit = (
                    "carrier-energy units"
                    if saved in ["e", "state_of_charge"]
                    else "carrier-flow units"
                )
                self.record(
                    variable + "_lower",
                    low - values,
                    ids,
                    unit,
                    component=component,
                    inequality=True,
                )
                self.record(
                    variable + "_upper",
                    values - high,
                    ids,
                    unit,
                    component=component,
                    inequality=True,
                )
                set_attr = variable + "_set"
                if set_attr in static:
                    target = self.dense(component, set_attr)
                    valid = np.isfinite(target)
                    if valid.any():
                        self.record(
                            variable + "_set",
                            np.where(valid, values - target, 0),
                            ids,
                            unit,
                            component=component,
                        )
            for ramp in ["ramp_limit_up", "ramp_limit_down"]:
                if ramp in static and (
                    static[ramp].notna().any() or not c.dynamic[ramp].empty
                ):
                    self.scope.append(
                        f"NOT CHECKED: {component} {ramp} constraints are configured"
                    )

        ids = self.ids["StorageUnit"]
        spill = self.result("StorageUnit", "spill")
        inflow = self.dense("StorageUnit", "inflow")
        self.record(
            "spill_lower", -spill, ids, "MW", component="StorageUnit", inequality=True
        )
        self.record(
            "spill_upper",
            spill - inflow,
            ids,
            "MW",
            component="StorageUnit",
            inequality=True,
        )

    def temporal(self):
        hours = self.n.snapshot_weightings.stores.to_numpy()
        ids = self.ids["Store"]
        static = self.n.stores.loc[ids]
        e = self.result("Store", "e")
        p = self.result("Store", "p")
        res = storage_residual(
            e,
            p,
            hours,
            self.dense("Store", "standing_loss"),
            static.e_cyclic.to_numpy(bool),
            static.e_initial.to_numpy(),
        )
        self.record(
            "energy_balance", res, ids, "carrier-energy units", component="Store"
        )
        ids = self.ids["StorageUnit"]
        static = self.n.storage_units.loc[ids]
        e = self.result("StorageUnit", "state_of_charge")
        change = (
            self.result("StorageUnit", "p_dispatch")
            / self.dense("StorageUnit", "efficiency_dispatch")
            - self.result("StorageUnit", "p_store")
            * self.dense("StorageUnit", "efficiency_store")
            - self.dense("StorageUnit", "inflow")
            + self.result("StorageUnit", "spill")
        )
        res = storage_residual(
            e,
            change,
            hours,
            self.dense("StorageUnit", "standing_loss"),
            static.cyclic_state_of_charge.to_numpy(bool),
            static.state_of_charge_initial.to_numpy(),
        )
        self.record("energy_balance", res, ids, "MWh", component="StorageUnit")

    def link_equations(self):
        c = self.n.components["Link"]
        ids = self.ids["Link"]
        for col in c.static.columns:
            if not col.startswith("bus") or not col[3:].isdigit() or col == "bus0":
                continue
            selected = ids[c.static.loc[ids, col].fillna("").ne("")]
            port = col[3:]
            efficiency = "efficiency" if port == "1" else "efficiency" + port
            r = self.result("Link", "p" + port, selected) + self.result(
                "Link", "p0", selected
            ) * self.dense("Link", efficiency, selected)
            self.record(
                "port_" + port + "_conversion",
                r,
                selected,
                "output-carrier-flow units",
                component="Link",
            )
        if "p" in c.dynamic and not c.dynamic["p"].empty:
            self.record(
                "p_equals_p0",
                self.result("Link", "p") - self.result("Link", "p0"),
                ids,
                "input-carrier-flow units",
                component="Link",
            )
        for component in ["Line", "Transformer"]:
            ids = self.ids[component]
            if not ids.empty:
                self.record(
                    "lossless_flow_conservation",
                    self.result(component, "p0") + self.result(component, "p1"),
                    ids,
                    "MW",
                    component=component,
                )

    def kvl(self):
        self.n.calculate_dependent_values()
        cycles = self.n.cycle_matrix(apply_weights=True)
        if cycles.empty:
            return
        total = np.zeros((len(self.sns), len(cycles.columns)))
        for component in cycles.index.unique("type"):
            branch = cycles.loc[component]
            total += self.result(component, "p0", branch.index) @ branch.to_numpy()
        self.record(
            "Kirchhoff_voltage_law_native_scaling",
            total * 1e5,
            cycles.columns,
            "cycle-equation units (native 1e5 scaling)",
            component="AC cycle",
        )

    def global_budgets(self):
        rows = []
        w = self.n.snapshot_weightings.generators.to_numpy()
        for name, glc in self.n.global_constraints.iterrows():
            carrier = glc.carrier_attribute
            if glc.type == "operational_limit":
                gen = self.ids["Generator"][
                    self.n.generators.loc[self.ids["Generator"], "carrier"].eq(carrier)
                ]
                value = float((self.result("Generator", "p", gen) * w[:, None]).sum())
                for component, cyclic, initial, result in [
                    ("Store", "e_cyclic", "e_initial", "e"),
                    (
                        "StorageUnit",
                        "cyclic_state_of_charge",
                        "state_of_charge_initial",
                        "state_of_charge",
                    ),
                ]:
                    c = self.n.components[component]
                    ci = self.ids[component]
                    selected = ci[
                        c.static.loc[ci, "carrier"].eq(carrier)
                        & ~c.static.loc[ci, cyclic]
                    ]
                    if not selected.empty:
                        value += float(
                            c.static.loc[selected, initial].sum()
                            - self.result(component, result, selected)[-1].sum()
                        )
            elif glc.type == "co2_atmosphere":
                emitting = self.n.carriers.index[self.n.carriers[carrier].ne(0)]
                bus_carrier = self.n.stores.bus.map(self.n.buses.carrier)
                selected = self.n.stores.index[
                    bus_carrier.isin(emitting) & ~self.n.stores.e_cyclic
                ]
                value = float(self.result("Store", "e", selected)[-1].sum())
            else:
                self.scope.append(
                    f"NOT CHECKED: global constraint {name}, type {glc.type}"
                )
                continue
            if glc.sense in ["==", "="]:
                violation = abs(value - glc.constant)
            elif glc.sense == "<=":
                violation = max(0, value - glc.constant)
            elif glc.sense == ">=":
                violation = max(0, glc.constant - value)
            else:
                raise ValueError(f"Unsupported sense {glc.sense}")
            rows.append(
                {
                    "constraint": name,
                    "type": glc.type,
                    "carrier_attribute": carrier,
                    "sense": glc.sense,
                    "lhs": value,
                    "rhs": float(glc.constant),
                    "violation": float(violation),
                    "saved_mu": float(glc.mu),
                }
            )
        ci = self.ids["Generator"]
        annual = (self.result("Generator", "p") * w[:, None]).sum(axis=0)
        for side in ["min", "max"]:
            limit = self.n.generators.loc[ci, "e_sum_" + side].to_numpy()
            selected = np.isfinite(limit)
            r = (
                limit[selected] - annual[selected]
                if side == "min"
                else annual[selected] - limit[selected]
            )
            self.record(
                "annual_energy_" + side,
                r,
                ci[selected],
                "MWh",
                component="Generator",
                inequality=True,
                timestamps=False,
            )
        return pd.DataFrame(rows)

    def shared_capacities(self, config):
        links, stores, gen = self.n.links, self.n.stores, self.n.generators
        for label, token in [
            ("battery_capacity_ratio", "battery"),
            (
                "TES_charger_ratio",
                "water tanks|water pits|aquifer thermal energy storage",
            ),
        ]:
            # Parentheses are essential when the pattern contains alternatives.
            chargers = links.index[
                links.index.str.contains(r"(?:" + token + r") charger")
                & links.p_nom_extendable
            ]
            dischargers = chargers.str.replace(" charger", " discharger", regex=False)
            if len(chargers):
                r = (
                    links.loc[chargers, "p_nom_opt"].to_numpy()
                    - links.loc[dischargers, "p_nom_opt"].to_numpy()
                    * links.loc[dischargers, "efficiency"].to_numpy()
                )
                self.record(
                    label,
                    r,
                    chargers,
                    "MW or carrier-flow units",
                    component="Link",
                    timestamps=False,
                )
        chargers = links.index[
            links.index.str.contains(r"water tanks charger|water pits charger")
            & links.p_nom_extendable
        ]
        store_ids = chargers.str.replace(" charger", "", regex=False)
        r = (
            stores.loc[store_ids, "e_nom_opt"].to_numpy()
            - links.loc[chargers, "energy to power ratio"].to_numpy()
            * links.loc[chargers, "p_nom_opt"].to_numpy()
        )
        self.record(
            "TES_energy_to_power_ratio",
            r,
            store_ids,
            "MWh",
            component="Store/Link",
            timestamps=False,
        )
        backwards = links.index[links.get("reversed", False) & links.p_nom_extendable]
        forwards = backwards.str.replace("-reversed", "", regex=False)
        r = (
            links.loc[backwards, "p_nom_opt"].to_numpy()
            - links.loc[forwards, "p_nom_opt"].to_numpy()
        )
        self.record(
            "bidirectional_capacity_sync",
            r,
            backwards,
            "MW or carrier-flow units",
            component="Link",
            timestamps=False,
        )
        factor = (
            config["renewable"]["solar"]["capacity_per_sqkm"]
            / config["renewable"]["solar-hsat"]["capacity_per_sqkm"]
        )
        solar = gen[gen.carrier.isin(["solar", "solar-hsat"]) & gen.p_nom_extendable]
        today = gen[(gen.carrier == "solar") & gen.p_nom_extendable]
        tracking = gen[gen.carrier == "solar-hsat"]
        rhs = (
            today.groupby("bus").p_nom_max.sum()
            - tracking.groupby("bus").p_nom.sum() * factor
        ).clip(lower=0)
        land = pd.Series(
            np.where(solar.carrier == "solar-hsat", factor, 1), index=solar.index
        )
        lhs = (solar.p_nom_opt * land).groupby(solar.bus).sum()
        valid = rhs.notna() & lhs.reindex(rhs.index).notna()
        self.record(
            "solar_shared_land",
            (lhs.reindex(rhs.index) - rhs)[valid].to_numpy(),
            rhs.index[valid],
            "solar-equivalent MW",
            component="Generator",
            inequality=True,
            timestamps=False,
        )


def coal_audit(n):
    ids = n.stores.index[n.stores.carrier == "coal"]
    if len(ids) != 1:
        raise ValueError("Expected one coal store")
    name = ids[0]
    store = n.stores.loc[name]
    e = n.stores_t.e[name]
    conditions = {
        "extendable": bool(store.e_nom_extendable),
        "cyclic": bool(store.e_cyclic),
        "zero_standing_loss": bool(store.standing_loss == 0),
        "zero_minimum_capacity": bool(store.e_nom_min == 0),
        "no_fixed_capacity": bool(pd.isna(store.get("e_nom_set", np.nan))),
        "no_fixed_energy_level": bool(pd.isna(store.get("e_set", np.nan))),
        "zero_minimum_energy_pu": bool(store.e_min_pu == 0),
        "maximum_energy_pu_one": bool(store.e_max_pu == 1),
        "zero_inventory_holding_cost": bool(
            store.get("marginal_cost_storage", 0) == 0
            and (
                name
                not in n.stores_t.get("marginal_cost_storage", pd.DataFrame()).columns
                or n.stores_t["marginal_cost_storage"][name].eq(0).all()
            )
        ),
        "no_time_varying_energy_limits": all(
            name not in n.stores_t[a].columns
            for a in ["e_min_pu", "e_max_pu", "standing_loss", "e_set"]
        ),
    }
    if not all(conditions.values()):
        raise ValueError(f"Coal shift conditions not established: {conditions}")
    cost = float(n.c.stores.periodized_cost.sel(name=name))
    capacity = float(store.e_nom_opt)
    maximum, minimum = float(e.max()), float(e.min())
    return {
        "store": name,
        "conditions_for_lossless_cyclic_shift": conditions,
        "periodized_cost_eur_per_MWh_year": cost,
        "capacity_MWh": capacity,
        "minimum_inventory_MWh": minimum,
        "maximum_inventory_MWh": maximum,
        "inventory_range_MWh": maximum - minimum,
        "current_capacity_cost_eur": cost * capacity,
        "capacity_only_avoidable_cost_eur": cost * max(0, capacity - maximum),
        "same_dispatch_shift_avoidable_cost_eur": cost
        * max(0, capacity - maximum + minimum),
        "shifted_capacity_cost_eur": cost * (maximum - minimum),
    }


def price_table(n):
    rows = []
    w = n.snapshot_weightings.objective.to_numpy()
    buses = n.buses.index[n.buses.carrier.eq("AC")]
    for bus in buses:
        price = n.buses_t.marginal_price[bus].to_numpy()
        if not np.isfinite(price).all():
            raise ValueError(f"Non-finite price {bus}")
        rows.append(
            {
                "bus": bus,
                "mean_eur_MWh": float(np.average(price, weights=w)),
                "min_eur_MWh": float(price.min()),
                "max_eur_MWh": float(price.max()),
                "p05_eur_MWh": float(np.quantile(price, 0.05)),
                "median_eur_MWh": float(np.median(price)),
                "p95_eur_MWh": float(np.quantile(price, 0.95)),
                "hours_price_above_200": float(w[price > 200].sum()),
            }
        )
    return pd.DataFrame(rows)


def run(root, output):
    output.mkdir(parents=True, exist_ok=False)
    metadata = json.loads((root / "summary.json").read_text())
    candidates = [
        root / "candidate_2035.nc",
        root / "solved_2035.nc",
        root / "networks" / "solved_2035.nc",
    ]
    path = next((path for path in candidates if path.is_file()), None)
    if path is None:
        raise FileNotFoundError(f"No saved solved network found under {root}")
    print("Loading candidate ...", flush=True)
    n = pypsa.Network(path)
    config_path = root / "effective_config.yaml"
    if config_path.is_file():
        config = yaml.safe_load(config_path.read_text())
        config_source = str(config_path)
    else:
        config = n.meta
        config_source = "saved_network_meta"
    if not isinstance(config, dict) or not {"solving", "renewable"} <= config.keys():
        raise ValueError("Saved network has no usable audit configuration")
    if n.objective is None or not np.isfinite(n.objective):
        raise ValueError("Saved candidate has no finite objective")
    if n._multi_invest or n.has_scenarios or isinstance(n.snapshots, pd.MultiIndex):
        raise ValueError("This audit supports deterministic single-period networks")
    a = Audit(n)
    for component, ids in a.ids.items():
        if len(ids) != len(n.components[component].static):
            raise ValueError("Time-varying/inactive assets require additional masking")
    if config["solving"]["options"].get("transmission_losses", 0):
        raise ValueError("Passive-branch transmission losses require another audit")
    print("Reconstructing objective ...", flush=True)
    costs, objective = a.costs()
    print("Checking nodal balances and operating limits ...", flush=True)
    buses, carriers = a.balance()
    a.limits()
    print("Checking storage, conversion, KVL and shared capacities ...", flush=True)
    a.temporal()
    a.link_equations()
    a.kvl()
    a.shared_capacities(config)
    globals_table = a.global_budgets()
    coal = coal_audit(n)
    prices = price_table(n)
    metrics = pd.DataFrame(a.metrics)
    worst = pd.DataFrame(a.worst).sort_values("violation", ascending=False)
    quality = metadata.get("solver_quality", {}).get("values", {})
    report = {
        "audit_type": "read_only_saved_candidate",
        "optimization_called": False,
        "pypsa_version": pypsa.__version__,
        "input_network": str(path),
        "input_sha256": sha256(path),
        "audit_code_sha256": sha256(Path(__file__)),
        "run_summary_sha256": sha256(root / "summary.json"),
        "effective_config_source": config_source,
        "effective_config_sha256": sha256(config_path)
        if config_path.is_file()
        else None,
        "source_sha256_recorded_by_run": metadata.get("source_sha256"),
        "scenario": metadata.get("scenario"),
        "snapshots": len(n.snapshots),
        "hours": float(n.snapshot_weightings.objective.sum()),
        "weights": {
            c: sorted(n.snapshot_weightings[c].unique().tolist())
            for c in n.snapshot_weightings
        },
        "objective": objective,
        "coal": coal,
        "reference_objective_eur": metadata.get("reference_objective"),
        "native_solver_quality": quality,
        "original_screening_issues": metadata.get("screening_issues"),
        "barrier_final": metadata.get("barrier_final"),
        "original_run_status": metadata.get("status"),
        "original_solver_options": metadata.get("solver_options"),
        "original_termination_condition": metadata.get("solve_result", {}).get(
            "termination_condition", metadata.get("termination_condition")
        ),
        "no_full_dual_certificate": True,
        "scope_notes": a.scope
        + [
            "NetCDF does not contain the full Linopy/native Gurobi constraint dual vector or reduced costs.",
            "Native row/column names and full dual vectors were not archived; quality indices alone do not identify the offending physical equation or variable.",
            "Physical checks do not establish LP optimality or certify welfare estimates.",
            "No result status or validation threshold was changed.",
        ],
    }
    report["coal_shifted_candidate_objective_eur"] = (
        objective["saved_objective_eur"]
        - coal["same_dispatch_shift_avoidable_cost_eur"]
    )
    costs.to_csv(output / "costs_by_asset.csv", index=False)
    costs.groupby(["component", "carrier", "cost_type"])[
        ["gross_eur", "existing_constant_eur", "native_objective_eur"]
    ].sum().to_csv(output / "costs_by_carrier.csv")
    metrics.to_csv(output / "constraint_checks.csv", index=False)
    worst.to_csv(output / "largest_residuals.csv", index=False)
    buses.sort_values("max_balance_error", ascending=False).to_csv(
        output / "bus_balances.csv", index=False
    )
    carriers.to_csv(output / "carrier_balances.csv", index=False)
    globals_table.to_csv(output / "global_budgets.csv", index=False)
    prices.to_csv(output / "zonal_prices.csv", index=False)
    # Keep full electricity prices and weights for later matched-series comparisons.
    electric_buses = n.buses.index[n.buses.carrier.eq("AC")]
    n.buses_t.marginal_price.loc[n.snapshots, electric_buses].to_csv(
        output / "electricity_prices.csv.gz", index_label="snapshot"
    )
    n.snapshot_weightings.to_csv(
        output / "snapshot_weightings.csv", index_label="snapshot"
    )
    (output / "audit.json").write_text(json.dumps(report, indent=2, allow_nan=False))
    print(
        json.dumps(
            {
                "objective": objective,
                "coal": coal,
                "largest_checks": metrics.nlargest(12, "max_violation").to_dict(
                    "records"
                ),
                "globals": globals_table.to_dict("records"),
                "scope_notes": a.scope,
            },
            indent=2,
        )
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate_directory", type=Path)
    parser.add_argument("--output", type=Path, default=Path("audit_output"))
    args = parser.parse_args()
    run(args.candidate_directory, args.output)
