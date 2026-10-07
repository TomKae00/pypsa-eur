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
        decommissioning = pd.to_numeric(frame["decommissioning_year"], errors="coerce")
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
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat0) * np.cos(lat1) * np.sin(dlon / 2.0) ** 2
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
    """Require an explicit receiving bus and validate its basic attributes.

    A country and a nearest cluster are not a bidding-zone/landing mapping.
    Use inspect_energy_islands.py on a composed baseline to find bus IDs.
    """
    exact = _optional_string(row.get("target_bus"))
    if not exact:
        raise ValueError(
            f"{row['link_id']}: set target_bus explicitly in the links CSV. "
            "Run scripts/inspect_energy_islands.py on the composed baseline. "
            "A link name containing DK2 does not identify the receiving zone."
        )
    if exact not in n.buses.index:
        raise KeyError(f"Configured target_bus {exact!r} does not exist")
    target = n.buses.loc[exact]
    if exact == str(row["hub_id"]) or _optional_string(target.get("energy_island")):
        raise ValueError(
            f"{row['link_id']}: target_bus must not be an energy-island hub"
        )
    if target.get("carrier") != "AC":
        raise ValueError(f"{row['link_id']}: target_bus {exact!r} must be an AC bus")
    country = _optional_string(row.get("target_country"))
    if country and target.get("country") != country:
        raise ValueError(f"{row['link_id']}: {exact!r} is not in {country}")
    prefix = _optional_string(row.get("target_bus_prefix"))
    if prefix and not exact.startswith(prefix):
        raise ValueError(
            f"{row['link_id']}: {exact!r} does not match prefix {prefix!r}"
        )
    return exact


def _select_profile_generator(
    n: pypsa.Network,
    carrier: str,
    configured_source: str,
) -> str:
    """Use an explicit proxy profile; bus proximity is not wind-site proximity."""
    if not configured_source:
        raise ValueError(
            "Set source_generator explicitly in the wind CSV. "
            "The inventory command lists available dynamic profiles. "
            "This is a documented proxy, not a dedicated Bornholm wind profile."
        )
    if configured_source not in n.generators.index:
        raise KeyError(
            f"Configured source_generator {configured_source!r} does not exist"
        )
    if n.generators.at[configured_source, "carrier"] != carrier:
        raise ValueError(
            f"Profile source {configured_source!r} does not have carrier {carrier!r}"
        )
    if configured_source not in n.generators_t.p_max_pu:
        raise ValueError(f"Profile source {configured_source!r} has no dynamic profile")
    profile = n.generators_t.p_max_pu[configured_source].reindex(n.snapshots)
    if not np.isfinite(profile).all() or not profile.between(0.0, 1.0).all():
        raise ValueError(
            f"Profile {configured_source!r} must cover every snapshot and lie in [0, 1]"
        )
    logger.warning("Using explicitly selected proxy profile %s", configured_source)
    return configured_source


def _subtract_generic_potential(
    n: pypsa.Network,
    carrier: str,
    capacity_mw: float,
    protected_generators: set[str],
    allocation: Mapping[str, float] | None = None,
    country: str | None = None,
    allocation_note: str = "",
) -> None:
    """Apply a reviewed allocation; never spill into a neighbouring country.

    This is capacity bookkeeping, not a geographic project-area exclusion.
    All rows are validated before any potential is changed.
    """
    if not isinstance(allocation, Mapping):
        raise ValueError(
            "Provide energy_islands.potential_allocation.<generator_id> as "
            "a mapping of generic generator IDs to MW. No automatic subtraction "
            "is performed. Do not disable subtraction to bypass this check."
        )
    reductions = {str(key): float(value) for key, value in allocation.items()}
    if any(not np.isfinite(value) or value <= 0 for value in reductions.values()):
        raise ValueError(
            "Every potential allocation must be finite and strictly positive"
        )
    removed = sum(reductions.values())
    if removed > capacity_mw + 1e-6:
        raise ValueError(
            f"Potential allocation exceeds project capacity {capacity_mw:g} MW"
        )
    if (
        not np.isclose(removed, capacity_mw, rtol=0, atol=1e-6)
        and not allocation_note.strip()
    ):
        raise ValueError(
            f"Allocation removes {removed:g} of {capacity_mw:g} MW. "
            "Set potential_allocation_notes.<generator_id> to document why only "
            "this overlap is represented (or why it has already been excluded). "
            "Do not fill the difference with unrelated potential."
        )
    for generator, reduction in reductions.items():
        if generator not in n.generators.index or generator in protected_generators:
            raise ValueError(f"Invalid generic potential source {generator!r}")
        source = n.generators.loc[generator]
        if source.carrier != carrier or not bool(source.p_nom_extendable):
            raise ValueError(f"{generator!r} must be an extendable {carrier} generator")
        if country and n.buses.at[source.bus, "country"] != country:
            raise ValueError(
                f"{generator!r} is outside project country {country}; review the resource mapping"
            )
        upper = float(source.p_nom_max)
        floor = max(float(source.p_nom_min), float(source.p_nom))
        if (
            not np.isfinite(upper)
            or not np.isfinite(floor)
            or upper - reduction < floor - 1e-6
        ):
            raise ValueError(
                f"{generator!r}: cannot remove {reduction:g} MW from p_nom_max={upper:g}; "
                f"existing/minimum capacity is {floor:g} MW. Review the resource representation."
            )
    for generator, reduction in reductions.items():
        n.generators.at[generator, "p_nom_max"] -= reduction
        logger.info(
            "Reduced explicitly mapped potential %s by %.1f MW", generator, reduction
        )
    if allocation_note:
        logger.info("Potential allocation basis: %s", allocation_note)


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
    potential_allocation: Mapping[str, Mapping[str, float]] | None = None,
    potential_allocation_notes: Mapping[str, str] | None = None,
) -> None:
    added_generators: set[str] = set()

    for _, row in wind.iterrows():
        generator_id = str(row["generator_id"])
        hub_id = str(row["hub_id"])
        if generator_id in n.generators.index:
            raise ValueError(f"Generator {generator_id!r} already exists")
        if hub_id not in n.buses.index:
            raise KeyError(f"Unknown hub_id {hub_id!r} for {generator_id}")

        carrier = str(row["carrier"])
        capacity_mw = float(row["capacity_mw"])
        source = _select_profile_generator(
            n,
            carrier,
            _optional_string(row.get("source_generator")),
        )
        profile = n.generators_t.p_max_pu[source].copy()

        if subtract_potential and _as_bool(row.get("subtract_potential"), True):
            _subtract_generic_potential(
                n,
                carrier,
                capacity_mw,
                protected_generators=added_generators,
                allocation=(potential_allocation or {}).get(generator_id),
                country=str(hubs.set_index("hub_id").at[hub_id, "country"]),
                allocation_note=_optional_string(
                    (potential_allocation_notes or {}).get(generator_id)
                ),
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
    cable_cost = (
        length_km
        * length_factor
        * (
            (1.0 - underwater_fraction)
            * float(costs.at["HVDC overhead", "capital_cost"])
            + underwater_fraction * float(costs.at["HVDC submarine", "capital_cost"])
        )
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
        configured_length = pd.to_numeric(row.get("length_km"), errors="coerce")
        if pd.isna(configured_length):
            # Estimate from the physical landing point, never a cluster centroid.
            coordinates = np.array(
                [hub.x, hub.y, row["target_x"], row["target_y"]], dtype=float
            )
            if (
                not np.isfinite(coordinates).all()
                or not np.isfinite(length_factor)
                or length_factor < 1
            ):
                raise ValueError(
                    f"{link_id}: finite coordinates and length_factor >= 1 are required"
                )
            length_km = float(_haversine_km(*coordinates)) * length_factor
            length_source = "landing_distance_times_factor"
            logger.warning(
                "%s: estimated cable route %.1f km; supply length_km for an actual route",
                link_id,
                length_km,
            )
        else:
            # An explicit length is the full route length: do not multiply again.
            length_km = float(configured_length)
            length_source = "configured_route"
        if not np.isfinite(length_km) or length_km <= 0:
            raise ValueError(f"{link_id}: route length must be finite and positive")

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
            1.0,
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
        n.links.loc[link_id, "length_source"] = length_source
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
    unknown_link_hubs = set(links.hub_id.astype(str)).difference(
        hubs.hub_id.astype(str)
    )
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
        potential_allocation=config.get("potential_allocation", {}),
        potential_allocation_notes=config.get("potential_allocation_notes", {}),
    )
    _add_links(
        n,
        links,
        costs,
        length_factor=float(config.get("length_factor", 1.25)),
    )
    return n