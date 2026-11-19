# SPDX-License-Identifier: MIT
"""Central 2035 matched-design welfare decomposition, in real EUR 2025."""

import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa

from scripts.audit_bornholm_candidate import Audit
from scripts.run_bornholm_matched import DESIGNS, project_legs


def extract(n):
    if isinstance(n.snapshots, pd.MultiIndex):
        raise ValueError("Expected a single planning year")
    w = n.snapshot_weightings.objective.to_numpy(float)
    if not (
        np.isclose(w.sum(), 8760) and np.allclose(w, n.snapshot_weightings.generators)
    ):
        raise ValueError("Check annual settlement weights")
    reader = Audit(n)
    costs, audit = reader.costs()
    if abs(audit["reconstruction_error_eur"]) > 2:
        raise ValueError(f"Objective reconstruction mismatch for this network: {audit}")
    costs["net"] = costs.gross_eur - (
        costs.existing_constant_eur if n.objective_constant else 0.0
    )
    amounts = costs.groupby(["component", "asset"]).net.sum()
    prices = n.buses_t.marginal_price
    bus_country = n.buses.country.copy()
    valid = bus_country.str.fullmatch("[A-Z]{2}", na=False) & ~bus_country.isin(
        ["EU", "NA"]
    )
    bus_country = bus_country.where(valid).fillna(
        n.buses.location.map(bus_country.where(valid))
    )
    electric = set(
        n.buses.index[n.buses.carrier.isin(["AC", "DC", "low voltage", "electricity"])]
    )
    ledger = defaultdict(float)
    wind_mwh = wind_revenue = cash_sum = 0.0
    cable_cost = defaultdict(float)
    cable_rent = defaultdict(float)
    for component, ids in reader.ids.items():
        if not len(ids):
            continue
        table = n.components[component].static.loc[ids]
        single = "bus" in table.columns
        ports = (
            ["bus"]
            if single
            else table.columns[table.columns.str.fullmatch(r"bus\d+")].tolist()
        )
        cash = pd.Series(0.0, index=ids)
        for port in ports:
            names = table.index[table[port].isin(n.buses.index)]
            if len(names) == 0:
                continue
            p = reader.result(component, "p" if single else "p" + port[3:], names)
            bp = prices.loc[:, table.loc[names, port].tolist()].to_numpy()
            if not np.isfinite(bp).all():
                raise ValueError(f"Missing price in {component} {port}")
            sign = table.loc[names, "sign"].to_numpy() if single else -1.0
            cash.loc[names] += (w @ (p * bp)) * sign
        cash_sum += cash.sum()
        for asset, row in table.iterrows():
            cashflow = cash[asset]
            cost = amounts.get((component, str(asset)), 0.0)
            if component == "Generator" and str(asset).startswith("BEI_WIND"):
                ledger["Developer", "Wind operating and investment surplus"] += (
                    cashflow - cost
                )
                wind_mwh += w @ reader.result(component, "p", pd.Index([asset]))[:, 0]
                wind_revenue += cashflow
                continue
            if component == "Link" and str(asset).startswith("BEI_TO_"):
                landing = "DK" if str(asset).startswith("BEI_TO_DK2") else "DE"
                cable_cost[landing] += cost
                cable_rent[landing] += cashflow
                continue
            countries = sorted(
                {
                    bus_country.get(row[col])
                    for col in ports
                    if pd.notna(bus_country.get(row[col]))
                }
            ) or ["Unallocated"]
            if component == "Store" and str(asset) == "co2 atmosphere":
                ledger["Unallocated", "CO2-cap rent after costs"] += cashflow - cost
                continue
            network = component in ("Line", "Transformer") or (
                component == "Link"
                and (
                    (row.bus0 in electric and row.bus1 in electric)
                    or "pipeline" in row.carrier
                )
            )
            if component == "Load":
                category = (
                    "Fixed electricity consumers"
                    if row.bus in electric
                    else "Other fixed demand"
                )
            elif network:
                category = "Background network net income"
            else:
                category = "Producer and storage surplus after costs"
            for ct in countries:
                ledger[ct, category] += (cashflow - cost) / len(countries)
    load_ids = reader.ids["Load"].sort_values()
    h = hashlib.sha256()
    h.update(
        pd.util.hash_pandas_object(n.snapshot_weightings, index=True).values.tobytes()
    )
    load_table = n.loads.loc[load_ids, ["bus", "carrier", "sign"]]
    h.update(pd.util.hash_pandas_object(load_table, index=True).values.tobytes())
    h.update(reader.dense("Load", "p_set", load_ids).astype("float64").tobytes())
    wind_ids = sorted(n.generators.index[n.generators.index.str.startswith("BEI_WIND")])
    wind_profiles = reader.dense("Generator", "p_max_pu", wind_ids)
    if not np.all(wind_profiles == wind_profiles[:, [0]]):
        raise ValueError("Different wind profiles across project generators")
    h.update(wind_profiles[:, 0].astype("float64").tobytes())
    h.update(pd.util.hash_pandas_object(n.snapshot_weightings.index).values.tobytes())
    budgets = n.global_constraints.drop(
        columns=["mu", "energy_island_allowance"], errors="ignore"
    ).copy()
    if "energy_island_allowance" in n.global_constraints:
        budgets["constant"] -= n.global_constraints.energy_island_allowance.fillna(0)
    h.update(
        pd.util.hash_pandas_object(budgets.sort_index(), index=True).values.tobytes()
    )
    error = (
        sum(ledger.values())
        + sum(cable_rent.values())
        - sum(cable_cost.values())
        + n.objective
        - cash_sum
    )
    if abs(error) > 3:
        raise ValueError(f"Unreconciled costs for this network: {error:,.1f} EUR")
    return dict(
        ledger=pd.Series(ledger, dtype=float),
        energy=wind_mwh,
        revenue=wind_revenue,
        cable_cost=cable_cost,
        cable_rent=cable_rent,
        signature=h.hexdigest(),
        objective=float(n.objective),
        cash=cash_sum,
    )


def accounts(d, strike, pooled_cost_share):
    """Equal rents; each shore funds its dedicated leg, shared cost is split equally."""
    values = d["ledger"].copy()
    support = strike * d["energy"] - d["revenue"]
    for ct, support_share in (("DK", 0.3), ("DE", 0.7)):
        values.at[ct, "Project cable rent"] = sum(d["cable_rent"].values()) / 2
        values.at[ct, "Project cable cost"] = -(
            (1 - pooled_cost_share) * d["cable_cost"][ct]
            + pooled_cost_share * sum(d["cable_cost"].values()) / 2
        )
        values.at[ct, "Net CfD transfer"] = -support_share * support
    values.at["Developer", "Net CfD transfer"] = support
    return values


def compare(root):
    """Read a single Slurm array's transferred outputs, never mix revisions."""
    data, runs = {}, []
    for name, parent, capacity in DESIGNS:
        directory = Path(root) / name
        report = json.loads((directory / "summary.json").read_text())
        if report["status"] not in {"checks_passed", "optimal_needs_review"}:
            raise ValueError(f"{name}: {report['status']}")
        n = pypsa.Network(directory / "candidate_2035.nc")
        if (
            report["scenario"]["name"] != name
            or n.meta.get("matched_design")
            != dict(name=name, parent=parent, dk_mw=capacity)
            or abs(report["objective"] - n.objective) > 2
        ):
            raise ValueError(f"Mismatched summary or network for {name}")
        project_legs(
            n, {"project": parent.removeprefix("central_")}, capacity, solved=True
        )
        d = data[name] = extract(n)
        runs.append(
            dict(
                design=name,
                status=report["status"],
                issues=report["screening_issues"],
                objective_eur=d["objective"],
                generation_twh=d["energy"] / 1e6,
                capture_price=d["revenue"] / d["energy"],
                cable_cost_million=sum(d["cable_cost"].values()) / 1e6,
                cable_rent_million=sum(d["cable_rent"].values()) / 1e6,
                source_hash=report["source_sha256"],
                code_hash=report["code_sha256"],
            )
        )
        del n
    if (
        len({d["signature"] for d in data.values()}) != 1
        or len({r["code_hash"] for r in runs}) != 1
    ):
        raise ValueError(
            "Demands, budgets, wind profiles, snapshot weights or code differ"
        )
    a, b, c = data
    rows = []
    for ref, design in ((a, b), (b, c), (a, c)):
        for strike in (75, 82.5, 90, 97.5):
            for pooled in (0.0, 1.0):
                delta = (
                    accounts(data[design], strike, pooled)
                    .subtract(accounts(data[ref], strike, pooled), fill_value=0)
                    .groupby(level=0)
                    .sum()
                    / 1e6
                )
                saving = (data[ref]["objective"] - data[design]["objective"]) / 1e6
                rows.append(
                    dict(
                        reference=ref,
                        design=design,
                        strike=strike,
                        pooled_cost_share=pooled,
                        system=saving,
                        DK=delta.DK,
                        DE=delta.DE,
                        developer=delta.Developer,
                        others=delta.drop(["DK", "DE", "Developer"]).sum(),
                        residual=delta.sum() - saving,
                    )
                )
    return pd.DataFrame(runs), pd.DataFrame(rows)
