# SPDX-License-Identifier: MIT
"""Dated guard for existing Swedish reactor locations; no capacity change."""

KEY = "historical_asset_corrections_audit"


def apply(n, settings, year):
    if not settings.get("enabled", False):
        return []
    if year != 2023:
        raise ValueError("Swedish reactor mapping is reviewed for 2023 only")
    previous = n.meta.get(KEY)
    if previous:
        if previous["year"] != year or previous["settings"] != settings:
            raise ValueError("Historical asset assumptions changed; rebuild")
        return previous["corrections"]
    target = [b for b in n.buses.index if b in {"SE3", "SE03"}]
    source = [b for b in n.buses.index if b in {"SE4", "SE04"}]
    changes = []
    generators = n.generators
    selected = generators.index[
        generators.bus.isin(source)
        & generators.carrier.eq("nuclear")
        & generators.p_nom.gt(0)
    ]
    if len(selected):
        if len(target) != 1 or generators.loc[selected, "p_nom_extendable"].any():
            raise ValueError("Historical correction requires SE3 and fixed reactors")
        for name in selected:
            changes.append(
                {
                    "generator": str(name),
                    "old_bus": str(generators.at[name, "bus"]),
                    "new_bus": target[0],
                    "capacity_MW": float(generators.at[name, "p_nom"]),
                    "basis": "Active Swedish 2023 reactors at Forsmark, Ringhals and Oskarshamn belong to SE3; no operating reactor in SE4",
                }
            )
            generators.at[name, "bus"] = target[0]
    n.meta[KEY] = {"year": year, "settings": settings, "corrections": changes}
    return changes
