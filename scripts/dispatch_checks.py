# SPDX-License-Identifier: MIT
"""Post-solve MW balances, renewable bounds and chronological storage equations."""

import numpy as np
import pandas as pd


def _dense(n, component, attr):
    frame = getattr(n, component)
    result = pd.DataFrame(
        np.broadcast_to(frame[attr].to_numpy(), (len(n.snapshots), len(frame))).copy(),
        index=n.snapshots,
        columns=frame.index,
    )
    dynamic = getattr(getattr(n, component + "_t"), attr)
    shared = dynamic.columns.intersection(frame.index)
    result.loc[:, shared] = dynamic.reindex(n.snapshots)[shared]
    return result


def _finite(values, label):
    if not np.isfinite(values.to_numpy()).all():
        raise ValueError("Missing or nonfinite solved values: " + label)


def storage_check(n, tolerance_MWh=0.01, relative_tolerance=1e-7):
    maximum = 0.0
    hours = n.snapshot_weightings.stores
    for comp, state_attr, initial_attr, cyclic_attr in [
        ("stores", "e", "e_initial", "e_cyclic"),
        (
            "storage_units",
            "state_of_charge",
            "state_of_charge_initial",
            "cyclic_state_of_charge",
        ),
    ]:
        frame = getattr(n, comp)
        dynamic = getattr(n, comp + "_t")
        if frame.empty:
            continue
        state = getattr(dynamic, state_attr).reindex(
            index=n.snapshots, columns=frame.index
        )
        _finite(state, comp + " state")
        previous = state.shift(1)
        previous.iloc[0] = frame[initial_attr].where(
            ~frame[cyclic_attr], state.iloc[-1]
        )
        retained = previous * (1 - _dense(n, comp, "standing_loss")).pow(hours, axis=0)
        if comp == "stores":
            expected = retained - dynamic.p.reindex(columns=frame.index).mul(
                hours, axis=0
            )
        else:
            charge = dynamic.p_store.reindex(columns=frame.index) * _dense(
                n, comp, "efficiency_store"
            )
            discharge = dynamic.p_dispatch.reindex(columns=frame.index) / _dense(
                n, comp, "efficiency_dispatch"
            )
            spill = dynamic.spill.reindex(
                index=n.snapshots, columns=frame.index, fill_value=0
            )
            expected = retained + (
                charge - discharge + _dense(n, comp, "inflow") - spill
            ).mul(hours, axis=0)
        _finite(expected, comp + " storage equation")
        error = (state - expected).abs()
        maximum = max(maximum, float(error.to_numpy().max(initial=0)))
        if (error > tolerance_MWh + relative_tolerance * state.abs()).any().any():
            raise ValueError("Chronological storage balance violated: " + comp)
    return maximum


def check(n, absolute_tolerance_MW=0.01, relative_tolerance=1e-7):
    balance = pd.DataFrame(0.0, index=n.snapshots, columns=n.buses.index)
    activity = balance.copy()

    def add(values, buses):
        nonlocal balance, activity
        _finite(values, "bus-port dispatch")
        grouped = values.T.groupby(buses.reindex(values.columns)).sum().T
        absolute = values.abs().T.groupby(buses.reindex(values.columns)).sum().T
        balance = balance.add(grouped, fill_value=0)
        activity = activity.add(absolute, fill_value=0)

    for component in ["generators", "loads", "stores", "storage_units"]:
        frame = getattr(n, component)
        if frame.empty:
            continue
        p = getattr(n, component + "_t").p.reindex(
            index=n.snapshots, columns=frame.index
        )
        signs = frame.get("sign", pd.Series(1.0, index=frame.index))
        add(p.mul(signs, axis=1), frame.bus)
    for component in ["lines", "transformers", "links"]:
        frame = getattr(n, component)
        dynamic = getattr(n, component + "_t")
        for i in range(5 if component == "links" else 2):
            column = f"bus{i}"
            if column not in frame:
                continue
            names = frame.index[frame[column].ne("")]
            if len(names):
                add(
                    -getattr(dynamic, f"p{i}").reindex(
                        index=n.snapshots, columns=names
                    ),
                    frame[column],
                )
    _finite(balance, "bus balance")
    excess = balance.abs() - absolute_tolerance_MW - relative_tolerance * activity
    if excess.to_numpy().max(initial=0) > 0:
        row, column = np.unravel_index(np.nanargmax(excess.to_numpy()), excess.shape)
        raise ValueError(
            f"Solved bus imbalance at {balance.columns[column]}, {balance.index[row]}: {balance.iloc[row, column]:.6g} MW"
        )
    upper_violation = 0.0
    if not n.generators.empty:
        cap = n.generators.p_nom.where(
            ~n.generators.p_nom_extendable, n.generators.p_nom_opt
        )
        bounds = _dense(n, "generators", "p_max_pu").mul(cap, axis=1)
        violation = n.generators_t.p - bounds
        _finite(violation, "generator availability")
        upper_violation = float(violation.to_numpy().max(initial=0))
        if (
            (violation > absolute_tolerance_MW + relative_tolerance * bounds.abs())
            .any()
            .any()
        ):
            raise ValueError("Solved generator exceeds its available MW")
    electric = n.buses.index[n.buses.carrier.isin(["AC", "DC", ""])]
    if len(electric):
        _finite(
            n.buses_t.marginal_price.reindex(index=n.snapshots, columns=electric),
            "electricity marginal prices",
        )
    storage_residual = storage_check(n, absolute_tolerance_MW, relative_tolerance)
    return {
        "maximum_absolute_bus_balance_MW": float(
            balance.abs().to_numpy().max(initial=0)
        ),
        "maximum_generator_upper_violation_MW": upper_violation,
        "maximum_storage_equation_residual_MWh": storage_residual,
        "absolute_tolerance_MW": absolute_tolerance_MW,
        "relative_tolerance": relative_tolerance,
    }
