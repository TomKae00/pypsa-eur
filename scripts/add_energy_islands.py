# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Add explicitly specified offshore energy islands to a PyPSA-Eur network."""

from __future__ import annotations

import logging
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pypsa

logger = logging.getLogger(__name__)


def _read_csv(path: str | Path, required: set[str]) -> pd.DataFrame:
    """Read and validate one energy-island input table."""
    frame = pd.read_csv(path)
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
    return frame


def _optional_string(value: Any, default: str = "") -> str:
    if pd.isna(value):
        return default
    value = str(value).strip()
    return value if value else default


def _as_bool(value: Any, default: bool = False) -> bool:
    if pd.isna(value):
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    value = str(value).strip().lower()
    if value in {"true", "yes", "y", "1"}:
        return True
    if value in {"false", "no", "n", "0"}:
        return False
    raise ValueError(f"Cannot interpret {value!r} as a boolean")


def _active_in_horizon(frame: pd.DataFrame, horizon: int | None) -> pd.DataFrame:
    """Filter assets that are commissioned by the current planning horizon."""
    if horizon is None or "commissioning_year" not in frame:
        return frame

    commissioning = pd.to_numeric(frame["commissioning_year"], errors="coerce")
    active = commissioning.isna() | commissioning.le(horizon)

    if "decommissioning_year" in frame:
        decommissioning = pd.to_numeric(
            frame["decommissioning_year"], errors="coerce"
        )
        active &= decommissioning.isna() | decommissioning.gt(horizon)

    return frame.loc[active].copy()


def _haversine_km(x0: Any, y0: Any, x1: Any, y1: Any) -> np.ndarray:
    """Return great-circle distances in kilometres for lon/lat coordinates."""
    lon0 = np.radians(np.asarray(x0, dtype=float))
    lat0 = np.radians(np.asarray(y0, dtype=float))
    lon1 = np.radians(np.asarray(x1, dtype=float))
    lat1 = np.radians(np.asarray(y1, dtype=float))
    dlon = lon1 - lon0
    dlat = lat1 - lat0
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat0) * np.cos(lat1) * np.sin(
        dlon / 2.0
    ) ** 2
    return 6371.0 * 2.0 * np.arcsin(np.sqrt(a))


def _capacity_attributes(capacity_mw: float, build_mode: str) -> dict[str, Any]:
    """Translate a build mode into PyPSA nominal-capacity attributes."""
    if capacity_mw <= 0:
        raise ValueError(f"capacity_mw must be positive, got {capacity_mw}")

    mode = build_mode.strip().lower()
    if mode == "fixed":
        # A forced investment: bounds are equal so that capital costs enter the
        # objective while the installed capacity remains exogenous.
        return {
            "p_nom_extendable": True,
            "p_nom_min": capacity_mw,
            "p_nom_max": capacity_mw,
        }
    if mode == "extendable":
        return {
            "p_nom_extendable": True,
            "p_nom_min": 0.0,
            "p_nom_max": capacity_mw,
        }
    if mode == "existing":
        # Existing assets are treated as sunk and therefore do not add capital
        # expenditure to the optimisation objective.
        return {"p_nom_extendable": False, "p_nom": capacity_mw}

    raise ValueError(
        f"Unknown build_mode {build_mode!r}; use fixed, extendable or existing"
    )


def _require_unique(frame: pd.DataFrame, column: str, table: str) -> None:
    duplicated = frame.loc[frame[column].duplicated(), column].tolist()
    if duplicated:
        raise ValueError(f"Duplicate {column} values in {table}: {duplicated}")


def _select_target_bus(n: pypsa.Network, row: pd.Series) -> str:
    """Select an exact bus or the closest eligible clustered AC bus."""
    exact = _optional_string(row.get("target_bus"))
    if exact:
        if exact not in n.buses.index:
            raise KeyError(f"Configured target_bus {exact!r} does not exist")
        return exact

    candidates = n.buses.copy()
    if "carrier" in candidates:
        candidates = candidates.loc[candidates.carrier.eq("AC")]

    country = _optional_string(row.get("target_country"))
    if country:
        if "country" not in candidates:
            raise KeyError("Network buses have no country column")
        candidates = candidates.loc[candidates.country.eq(country)]

    prefix = _optional_string(row.get("target_bus_prefix"))
    if prefix:
        candidates = candidates.loc[candidates.index.str.startswith(prefix)]

    candidates = candidates.dropna(subset=["x", "y"])
    if candidates.empty:
        raise ValueError(
            "No target bus matches "
            f"country={country!r}, prefix={prefix!r}"
        )

    target_x = float(row["target_x"])
    target_y = float(row["target_y"])
    distance = pd.Series(
        _haversine_km(candidates.x, candidates.y, target_x, target_y),
        index=candidates.index,
    )
    selected = str(distance.idxmin())
    logger.info(
        "Selected mainland bus %s for %s (%.1f km from target coordinate)",
        selected,
        row["link_id"],
        distance.loc[selected],
    )
    return selected


def _select_profile_generator(
    n: pypsa.Network,
    carrier: str,
    x: float,
    y: float,
    configured_source: str,
) -> str:
    """Select the configured or nearest existing generator profile."""
    if configured_source:
        if configured_source not in n.generators.index:
            raise KeyError(
                f"Configured source_generator {configured_source!r} does not exist"
            )
        return configured_source

    candidates = n.generators.loc[n.generators.carrier.eq(carrier)].copy()
    dynamic_profiles = set(n.generators_t.p_max_pu.columns)
    candidates = candidates.loc[candidates.index.isin(dynamic_profiles)]
    if candidates.empty:
        raise ValueError(
            f"No existing {carrier!r} generator has a dynamic p_max_pu profile"
        )

    coordinates = n.buses.reindex(candidates.bus)[["x", "y"]]
    coordinates.index = candidates.index
    valid = coordinates.notna().all(axis=1)
    candidates = candidates.loc[valid]
    coordinates = coordinates.loc[valid]
    if candidates.empty:
        raise ValueError(f"No {carrier!r} profile source has valid coordinates")

    distance = pd.Series(
        _haversine_km(coordinates.x, coordinates.y, x, y),
        index=candidates.index,
    )
    selected = str(distance.idxmin())
    logger.info(
        "Using profile from %s (associated bus %.1f km from hub)",
        selected,
        distance.loc[selected],
    )
    return selected


def _subtract_generic_potential(
    n: pypsa.Network,
    carrier: str,
    capacity_mw: float,
    x: float,
    y: float,
    protected_generators: set[str],
) -> None:
    """Remove dedicated project capacity from nearby generic wind potential."""
    candidates = n.generators.loc[
        n.generators.carrier.eq(carrier)
        & n.generators.p_nom_extendable.fillna(False).astype(bool)
        & ~n.generators.index.isin(protected_generators)
    ].copy()
    p_nom_max = pd.to_numeric(candidates.p_nom_max, errors="coerce")
    candidates = candidates.loc[p_nom_max.notna() & np.isfinite(p_nom_max)]
    if candidates.empty:
        raise ValueError(f"No finite generic {carrier!r} potential can be reduced")

    coordinates = n.buses.reindex(candidates.bus)[["x", "y"]]
    coordinates.index = candidates.index
    valid = coordinates.notna().all(axis=1)
    candidates = candidates.loc[valid]
    coordinates = coordinates.loc[valid]
    distance = pd.Series(
        _haversine_km(coordinates.x, coordinates.y, x, y),
        index=candidates.index,
    ).sort_values()

    remaining = float(capacity_mw)
    for generator in distance.index:
        upper = float(n.generators.at[generator, "p_nom_max"])
        lower_raw = n.generators.at[generator, "p_nom_min"]
        lower = 0.0 if pd.isna(lower_raw) else float(lower_raw)
        reducible = max(0.0, upper - lower)
        reduction = min(remaining, reducible)
        if reduction:
            n.generators.at[generator, "p_nom_max"] = upper - reduction
            remaining -= reduction
            logger.info(
                "Reduced generic potential %s by %.1f MW", generator, reduction
            )
        if remaining <= 1e-6:
            break

    if remaining > 1e-6:
        raise ValueError(
            f"Could not subtract {capacity_mw:.1f} MW of {carrier} potential; "
            f"{remaining:.1f} MW remains"
        )


def _add_hubs(n: pypsa.Network, hubs: pd.DataFrame) -> None:
    for _, row in hubs.iterrows():
        hub_id = str(row["hub_id"])
        if hub_id in n.buses.index:
            raise ValueError(f"Bus {hub_id!r} already exists")

        n.add(
            "Bus",
            hub_id,
            carrier=_optional_string(row.get("carrier"), "AC"),
            x=float(row["x"]),
            y=float(row["y"]),
            country=str(row["country"]),
            location=hub_id,
        )
        n.buses.loc[hub_id, "energy_island"] = hub_id


def _add_wind(
    n: pypsa.Network,
    wind: pd.DataFrame,
    hubs: pd.DataFrame,
    costs: pd.DataFrame,
    subtract_potential: bool,
) -> None:
    hub_coordinates = hubs.set_index("hub_id")[["x", "y"]]
    added_generators: set[str] = set()

    for _, row in wind.iterrows():
        generator_id = str(row["generator_id"])
        hub_id = str(row["hub_id"])
        if generator_id in n.generators.index:
            raise ValueError(f"Generator {generator_id!r} already exists")
        if hub_id not in n.buses.index:
            raise KeyError(f"Unknown hub_id {hub_id!r} for {generator_id}")

        x = float(hub_coordinates.at[hub_id, "x"])
        y = float(hub_coordinates.at[hub_id, "y"])
        carrier = str(row["carrier"])
        capacity_mw = float(row["capacity_mw"])
        source = _select_profile_generator(
            n,
            carrier,
            x,
            y,
            _optional_string(row.get("source_generator")),
        )
        profile = n.generators_t.p_max_pu[source].copy()

        if subtract_potential and _as_bool(row.get("subtract_potential"), True):
            _subtract_generic_potential(
                n,
                carrier,
                capacity_mw,
                x,
                y,
                protected_generators=added_generators,
            )

        cost_technology = _optional_string(
            row.get("capital_cost_technology"), "offwind"
        )
        operating_technology = carrier.split("-", 1)[0]
        capacity_attrs = _capacity_attributes(
            capacity_mw,
            _optional_string(row.get("build_mode"), "fixed"),
        )

        n.add(
            "Generator",
            generator_id,
            bus=hub_id,
            carrier=carrier,
            capital_cost=float(costs.at[cost_technology, "capital_cost"]),
            marginal_cost=float(costs.at[operating_technology, "marginal_cost"]),
            efficiency=float(costs.at[operating_technology, "efficiency"]),
            lifetime=float(costs.at[operating_technology, "lifetime"]),
            p_max_pu=profile,
            **capacity_attrs,
        )
        n.generators.loc[generator_id, "energy_island"] = hub_id
        n.generators.loc[generator_id, "profile_source"] = source
        added_generators.add(generator_id)


def _link_capital_cost(
    costs: pd.DataFrame,
    length_km: float,
    underwater_fraction: float,
    length_factor: float,
) -> float:
    cable_cost = length_km * length_factor * (
        (1.0 - underwater_fraction)
        * float(costs.at["HVDC overhead", "capital_cost"])
        + underwater_fraction
        * float(costs.at["HVDC submarine", "capital_cost"])
    )
    converter_cost = float(costs.at["HVDC inverter pair", "capital_cost"])
    return cable_cost + converter_cost


def _add_links(
    n: pypsa.Network,
    links: pd.DataFrame,
    costs: pd.DataFrame,
    length_factor: float,
) -> None:
    for _, row in links.iterrows():
        link_id = str(row["link_id"])
        hub_id = str(row["hub_id"])
        if link_id in n.links.index:
            raise ValueError(f"Link {link_id!r} already exists")
        if hub_id not in n.buses.index:
            raise KeyError(f"Unknown hub_id {hub_id!r} for {link_id}")

        target_bus = _select_target_bus(n, row)
        hub = n.buses.loc[hub_id]
        target = n.buses.loc[target_bus]

        configured_length = pd.to_numeric(row.get("length_km"), errors="coerce")
        if pd.isna(configured_length):
            length_km = float(
                _haversine_km(hub.x, hub.y, target.x, target.y)
            )
        else:
            length_km = float(configured_length)

        underwater_fraction = float(row.get("underwater_fraction", 1.0))
        if not 0.0 <= underwater_fraction <= 1.0:
            raise ValueError(
                f"underwater_fraction for {link_id} must lie between zero and one"
            )

        bidirectional = _as_bool(row.get("bidirectional"), True)
        efficiency = float(row.get("efficiency", 1.0))
        if bidirectional and not np.isclose(efficiency, 1.0):
            raise ValueError(
                f"Bidirectional Link {link_id} must use efficiency=1.0. "
                "A single PyPSA Link cannot represent symmetric losses in both directions."
            )

        capacity_mw = float(row["capacity_mw"])
        capacity_attrs = _capacity_attributes(
            capacity_mw,
            _optional_string(row.get("build_mode"), "fixed"),
        )
        capital_cost = _link_capital_cost(
            costs,
            length_km,
            underwater_fraction,
            length_factor,
        )

        n.add(
            "Link",
            link_id,
            bus0=hub_id,
            bus1=target_bus,
            carrier="DC",
            length=length_km,
            capital_cost=capital_cost,
            p_min_pu=-1.0 if bidirectional else 0.0,
            efficiency=efficiency,
            **capacity_attrs,
        )
        n.links.loc[link_id, "underwater_fraction"] = underwater_fraction
        n.links.loc[link_id, "energy_island"] = hub_id
        n.links.loc[link_id, "target_coordinate"] = (
            f"{float(row['target_x']):.5f},{float(row['target_y']):.5f}"
        )


def main(
    n: pypsa.Network,
    inputs: Mapping[str, Any],
    config: Mapping[str, Any],
    costs: pd.DataFrame,
    current_horizon: int | None = None,
) -> pypsa.Network:
    """Add the configured energy-island case to ``n`` in place."""
    if not config.get("enabled", False):
        return n

    hubs = _read_csv(
        inputs["energy_island_hubs"],
        {"hub_id", "country", "x", "y"},
    )
    wind = _read_csv(
        inputs["energy_island_wind"],
        {"generator_id", "hub_id", "carrier", "capacity_mw"},
    )
    links = _read_csv(
        inputs["energy_island_links"],
        {
            "case",
            "link_id",
            "hub_id",
            "target_country",
            "target_bus_prefix",
            "target_x",
            "target_y",
            "capacity_mw",
        },
    )

    _require_unique(hubs, "hub_id", "hubs")
    _require_unique(wind, "generator_id", "wind")

    case = str(config.get("case", "hybrid"))
    links = links.loc[links["case"].isin([case, "all"])].copy()
    _require_unique(links, "link_id", f"links for case {case}")

    wind = _active_in_horizon(wind, current_horizon)
    links = _active_in_horizon(links, current_horizon)
    active_hubs = set(wind.hub_id.astype(str)) | set(links.hub_id.astype(str))
    hubs = hubs.loc[hubs.hub_id.astype(str).isin(active_hubs)].copy()

    if hubs.empty:
        logger.info(
            "No energy-island assets are active in planning horizon %s",
            current_horizon,
        )
        return n

    unknown_wind_hubs = set(wind.hub_id.astype(str)).difference(hubs.hub_id.astype(str))
    unknown_link_hubs = set(links.hub_id.astype(str)).difference(hubs.hub_id.astype(str))
    if unknown_wind_hubs or unknown_link_hubs:
        raise ValueError(
            "Unknown hub references: "
            f"wind={sorted(unknown_wind_hubs)}, links={sorted(unknown_link_hubs)}"
        )

    logger.info(
        "Adding energy-island case %s for horizon %s: %d hub(s), %d wind farm(s), %d link(s)",
        case,
        current_horizon,
        len(hubs),
        len(wind),
        len(links),
    )
    _add_hubs(n, hubs)
    _add_wind(
        n,
        wind,
        hubs,
        costs,
        subtract_potential=bool(config.get("subtract_generic_potential", True)),
    )
    _add_links(
        n,
        links,
        costs,
        length_factor=float(config.get("length_factor", 1.25)),
    )
    return n
