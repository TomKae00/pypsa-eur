# SPDX-License-Identifier: MIT
"""
Danish historical effective availability fitted separately in DK1 and DK2.

Observed production contains outages and curtailment. This fit is an empirical
hindcast adjustment, not an estimate of unconstrained meteorological potential.
Weather chronology, zero-output hours and resource-class ordering are retained.
"""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

KEY = "hindcast_danish_wind_audit"


def observations_path(settings, year):
    return str(settings["observations"]).format(year=int(year))


def _payload(settings, year):
    source = Path(observations_path(settings, year))
    payload = {
        "version": 1,
        "year": int(year),
        "settings": settings,
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    }
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    return payload, digest


def apply(n, settings, year):
    if not settings.get("enabled", False):
        return None
    if n.generators.p_nom_extendable.any() or not n.investment_periods.empty:
        raise ValueError("Historical Danish wind needs a fixed historical fleet")
    sns = pd.DatetimeIndex(pd.to_datetime(n.snapshots, utc=True))
    if not (sns.year == year).all() or not n.snapshot_weightings.generators.eq(1).all():
        raise ValueError(
            "Apply historical Danish wind at hourly native resolution in its actual year"
        )
    payload, digest = _payload(settings, year)
    previous = n.meta.get(KEY)
    if previous:
        if previous.get("fingerprint") != digest:
            raise ValueError(
                "Danish wind correction changed on a prepared network; recompose"
            )
        return previous
    data = pd.read_csv(observations_path(settings, year))
    data["snapshot"] = pd.to_datetime(data.snapshot, utc=True)
    if data.duplicated(["snapshot", "zone", "carrier"]).any():
        raise ValueError("Duplicate Danish wind observation")
    rows = []
    for area, bus in [("DK1", "DKW1"), ("DK2", "DKE1")]:
        for technology in settings.get("carriers", ["onwind"]):
            if technology != "onwind":
                raise ValueError(
                    "Only onshore calibration is reviewed; offshore native losses must not be fitted twice"
                )
            selected = n.generators.index[
                n.generators.bus.isin([area, bus]) & n.generators.carrier.eq(technology)
            ]
            if selected.empty or n.generators.loc[selected, "p_nom"].sum() <= 0:
                raise ValueError(f"Missing historical wind fleet: {area} {technology}")
            observed = (
                data[data.zone.eq(area) & data.carrier.eq(technology)]
                .set_index("snapshot")
                .reindex(sns)
            )
            obs_values = observed[
                ["generation_MW", "capacity_MW", "spot_price_EUR_MWh"]
            ].to_numpy(float)
            if (
                not np.isfinite(obs_values).all()
                or (obs_values[:, :2] < 0).any()
                or (obs_values[:, 1] <= 0).any()
            ):
                raise ValueError(f"Incomplete {area} {technology} observations")
            capacities = n.generators.loc[selected, "p_nom"].to_numpy(float)
            fractions = capacities / capacities.sum()
            raw = pd.DataFrame(
                np.broadcast_to(
                    n.generators.loc[selected, "p_max_pu"].to_numpy(),
                    (len(n.snapshots), len(selected)),
                ).copy(),
                index=n.snapshots,
                columns=selected,
            )
            shared = n.generators_t.p_max_pu.columns.intersection(selected)
            raw.loc[:, shared] = n.generators_t.p_max_pu.reindex(n.snapshots)[
                shared
            ].to_numpy()
            values = raw.to_numpy(float)
            if (
                not np.isfinite(values).all()
                or values.min() < -1e-8
                or values.max() > 1 + 1e-8
            ):
                raise ValueError("Invalid native wind availability")
            values = np.clip(values, 0, 1)
            actual_capacity = observed.capacity_MW.to_numpy(float)
            reference_capacity = actual_capacity.max()
            potential = (values @ fractions) * actual_capacity
            fit = (
                observed.spot_price_EUR_MWh.to_numpy()
                >= settings.get("minimum_fit_price_EUR_MWh", 5.0)
            ) & ((values @ fractions) >= 0.05)
            if fit.sum() < settings.get("minimum_fit_hours", 1000):
                raise ValueError(
                    f"Too few positive-price fit hours: {area} {technology}"
                )
            target = float(observed.generation_MW.to_numpy()[fit].sum())
            lower, upper = settings.get("factor_bounds", [0.3, 1.5])

            def energy(factor):
                return float(
                    ((np.clip(values * factor, 0, 1) @ fractions) * actual_capacity)[
                        fit
                    ].sum()
                )

            if not energy(lower) <= target <= energy(upper) or target <= 0:
                raise ValueError(
                    f"Danish wind fit outside reviewed factor bounds: {area} {technology}"
                )
            for _ in range(55):
                middle = (lower + upper) / 2
                if energy(middle) < target:
                    lower = middle
                else:
                    upper = middle
            factor = (lower + upper) / 2
            corrected = (
                np.clip(values * factor, 0, 1)
                * (actual_capacity / reference_capacity)[:, None]
            )
            n.generators.loc[selected, "p_nom"] = fractions * reference_capacity
            n.generators.loc[selected, "p_nom_min"] = fractions * reference_capacity
            for i, name in enumerate(selected):
                n.generators_t.p_max_pu[name] = corrected[:, i]
            after = (corrected @ fractions) * reference_capacity
            rows.append(
                {
                    "zone": area,
                    "carrier": technology,
                    "factor": factor,
                    "original_capacity_MW": float(capacities.sum()),
                    "reference_capacity_MW": float(reference_capacity),
                    "fit_hours": int(fit.sum()),
                    "observed_generation_MWh": float(observed.generation_MW.sum()),
                    "native_potential_MWh_with_actual_capacity": float(potential.sum()),
                    "corrected_effective_available_MWh": float(after.sum()),
                    "fit_relative_energy_error": abs(energy(factor) / target - 1),
                    "hourly_correlation": float(
                        np.corrcoef(after, observed.generation_MW.to_numpy())[0, 1]
                    ),
                }
            )
    audit = {
        **payload,
        "fingerprint": digest,
        "series": rows,
        "interpretation": "Empirical historical effective availability. Positive-price fitting reduces, but does not eliminate, curtailment contamination. Validate post-solve curtailment and generation.",
        "offshore_policy": "Unchanged; do not duplicate the native 0.8855 offshore correction.",
    }
    n.meta[KEY] = audit
    return audit


def validate(n, settings, year):
    if not settings.get("enabled", False):
        return
    audit = n.meta.get(KEY)
    if not audit or audit.get("year") != int(year):
        raise ValueError(
            "Historical Danish wind was not corrected at composition; rebuild"
        )
    if _payload(settings, year)[1] != audit["fingerprint"]:
        raise ValueError("Danish wind settings or observations changed; recompose")
