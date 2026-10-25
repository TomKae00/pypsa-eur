# SPDX-License-Identifier: MIT
"""
Reproduce the committed 2035 project subset from a pinned converter checkout.

Run from the repository root:
python scripts/prepare_bornholm_grid.py /path/to/tyndp2026-to-pypsa
No network requests or automatic updates are performed.
"""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import Geod

SOURCE_URL = "https://github.com/p-glaum/tyndp2026-to-pypsa"
SOURCE_COMMIT = "7c207622c12c47447907d67a3a4de4ae53d633b4"
STATUSES = {"under_construction", "in_permitting", "planned_not_yet_permitted"}
TABLES = ("new_lines", "new_links", "upgraded_lines", "upgraded_links")


def correct_harmony(tables, source):
    """Replace obsolete offshore DC attributes with an explicit AC corridor proxy."""
    name = "TYNDP2026_170_1034"
    old = tables["new_links"].loc[name]
    route_path = source / "tyndp/entsoe/2026/tyndp2026_maps_all.geojson"
    features = json.loads(route_path.read_text())["features"]
    feature = next(f for f in features if f["properties"].get("investmentId") == 1034)
    if "HVAC" not in feature["properties"]["investmentName"]:
        raise ValueError("Review the Harmony Link project definition")
    # The complete map includes the existing Ostroleka-Elk corridor, which is
    # refurbished, not a second new circuit. Start at the OSM 0.7 Elk Bis station
    # and retain the map's remaining route towards Gizai. These are proxies for
    # a zonal study, not engineering coordinates/impedances of the planned cable.
    route = [(22.301858, 53.757587), *feature["geometry"]["coordinates"][2:]]
    x, y = zip(*route)
    length = round(Geod(ellps="WGS84").line_length(x, y) / 1000, 1)
    row = dict(
        project_status=old.project_status,
        build_year=old.build_year,
        length=length,
        underground=True,
        v_nom=220.0,
        num_parallel=2,
        type="Al/St 240/40 2-bundle 220.0",
        x0=x[0],
        y0=y[0],
        x1=x[-1],
        y1=y[-1],
        tags='"origin"=>"tyndp", "version"=>"2026", "project_id"=>"170", "investment_id"=>"1034", "location0"=>"Elk Bis", "location1"=>"Gizai", "attribute_source"=>"tyndp2026_description", "geometry_source"=>"osm0.7_ElkBis_and_entsoe2026_map", "length_source"=>"partial_route_geodesic", "electrical_parameters"=>"native_220kV_two_circuit_proxy", "url"=>"https://tyndp2026.entsoe.eu/projects-map/Transmission/170"',
    )
    tables["new_links"] = tables["new_links"].drop(name)
    tables["new_lines"] = pd.concat(
        [tables["new_lines"], pd.DataFrame([row], index=[name])]
    )
    return route_path


def select_projects(frame):
    """Return selected rows and a reason for every excluded converted row."""
    reasons = pd.Series("included", index=frame.index)
    reasons.loc[~frame.project_status.isin(STATUSES)] = "status_not_selected"
    years = pd.to_numeric(frame.build_year, errors="coerce")
    reasons.loc[years.isna() | years.gt(2035)] = "missing_year_or_after_2035"
    # Entire project 1106, independent of its commissioning dates or ratings.
    reasons.loc[frame.index.str.startswith("TYNDP2026_1106_")] = (
        "bornholm_modelled_separately"
    )
    selected = frame.loc[reasons.eq("included")].copy()
    for field in ["length", "x0", "y0", "x1", "y1"]:
        if not np.isfinite(selected[field]).all():
            raise ValueError(f"Missing/nonfinite {field} in selected projects")
    if (selected.length <= 0).any():
        raise ValueError("Selected branch lengths must be positive")
    if "p_nom" in selected and (
        not np.isfinite(selected.p_nom).all() or selected.p_nom.le(0).any()
    ):
        raise ValueError("A fixed planned DC project needs a positive finite rating")
    return selected, reasons


def prepare(source, destination):
    commit = subprocess.check_output(
        ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
    ).strip()
    if commit != SOURCE_COMMIT:
        raise ValueError(f"Expected converter commit {SOURCE_COMMIT}, got {commit}")
    # A pinned HEAD alone does not protect against edits in the checkout.
    dirty = subprocess.check_output(
        ["git", "-C", str(source), "status", "--porcelain", "--untracked-files=no"],
        text=True,
    )
    if dirty:
        raise ValueError("Use a clean converter checkout")
    destination.mkdir(parents=True, exist_ok=True)
    audit = destination / "audit"
    audit.mkdir(exist_ok=True)
    manifest = {
        "source_url": SOURCE_URL,
        "source_commit": commit,
        "horizon": 2035,
        "statuses": sorted(STATUSES),
        "excluded_project_ids": [1106],
        "input_sha256": {},
        "selected_counts": {},
    }
    decisions = []
    output_tables = {}
    paths = [Path("transmission_projects/tyndp2026") / f"{t}.csv" for t in TABLES]
    paths += [Path("audit/dropped_investments.csv")]
    for relative in paths:
        path = source / relative
        manifest["input_sha256"][str(relative)] = hashlib.sha256(
            path.read_bytes()
        ).hexdigest()
        frame = pd.read_csv(path, index_col=0 if relative.stem in TABLES else None)
        if relative.stem in TABLES:
            selected, reasons = select_projects(frame)
            output_tables[relative.stem] = selected
            decisions.append(
                reasons.rename("decision").to_frame().assign(table=relative.stem)
            )
        else:
            # Converter status spelling differs before export. Retain omissions,
            # including non-branch assets, so coverage is never implied complete.
            eligible = frame.build_year.le(2035) & frame.project_status.isin(
                STATUSES | {"in_planning"}
            )
            unconverted = frame.loc[eligible].copy()
    # A material cross-border omission: Celtic has no route in the converter's
    # 2026 map snapshot. Reuse the reviewed endpoint pair already in this repo's
    # TYNDP2020 table, retaining the 2026 status/year and current technical data.
    legacy_path = (
        Path(__file__).resolve().parents[1]
        / "data/transmission_projects/tyndp2020/new_links.csv"
    )
    legacy = pd.read_csv(legacy_path, index_col=0).loc["TYNDP2020_17"]
    if '"tyndp2020_invest_id"=>"810"' not in legacy.tags:
        raise ValueError("Celtic fallback no longer identifies investment 810")
    celtic = unconverted.loc[unconverted.investment_id.eq(810)].squeeze()
    if celtic.technology != "DC" or celtic.p_nom != legacy.p_nom:
        raise ValueError("Review Celtic technology/rating before using old endpoints")
    name = f"TYNDP2026_{int(celtic.project_id)}_810"
    if name in output_tables["new_links"].index:
        raise ValueError("Celtic is already converted; remove the geometry fallback")
    row = {key: legacy[key] for key in ["x0", "y0", "x1", "y1"]}
    row.update(
        project_status=celtic.project_status,
        length=celtic.length,
        build_year=celtic.build_year,
        p_nom=celtic.p_nom,
        underground=True,
        tags='"origin"=>"tyndp", "version"=>"2026", "project_id"=>"107", "investment_id"=>"810", "geometry_source"=>"pypsa-eur_tyndp2020_TYNDP2020_17", "attribute_source"=>"tyndp2024", "url"=>"https://tyndp2026.entsoe.eu/projects-map/Transmission/107"',
    )
    output_tables["new_links"] = pd.concat(
        [output_tables["new_links"], pd.DataFrame([row], index=[name])]
    )
    decisions.append(
        pd.DataFrame(
            {"decision": ["included_geometry_fallback"], "table": ["new_links"]},
            index=[name],
        )
    )
    route_path = correct_harmony(output_tables, source)
    manifest["input_sha256"][str(route_path.relative_to(source))] = hashlib.sha256(
        route_path.read_bytes()
    ).hexdigest()
    manifest["corrected_investment_ids"] = [1034]
    for frame in decisions:
        if "TYNDP2026_170_1034" in frame.index:
            frame.loc["TYNDP2026_170_1034", ["decision", "table"]] = [
                "included_ac_corridor_proxy",
                "new_lines",
            ]
    # Known technical/representation conflicts must not silently become new DC
    # additions if endpoint matching changes with a later native base network.
    excluded = {
        "TYNDP2026_299_1458": "excluded_multiterminal_replacement_needs_review",
        "TYNDP2026_1112_1799": "excluded_new_substation_endpoint_needs_review",
    }
    output_tables["new_links"] = output_tables["new_links"].drop(list(excluded))
    for frame in decisions:
        for name, reason in excluded.items():
            if name in frame.index:
                frame.loc[name, "decision"] = reason
    manifest["technical_review_exclusions"] = excluded
    unconverted["resolution"] = "not_included"
    unconverted.loc[unconverted.investment_id.eq(810), "resolution"] = (
        "recovered_from_native_2020_endpoints"
    )
    unconverted.to_csv(audit / "unconverted_2035.csv", index=False)
    manifest["unconverted_2035_rows"] = len(unconverted)
    manifest["recovered_investment_ids"] = [810]
    manifest["input_sha256"]["native/tyndp2020/new_links.csv"] = hashlib.sha256(
        legacy_path.read_bytes()
    ).hexdigest()
    for table, frame in output_tables.items():
        frame.sort_index().to_csv(destination / f"{table}.csv")
        manifest["selected_counts"][table] = len(frame)
    pd.concat(decisions).rename_axis("project").to_csv(audit / "selection.csv")
    (audit / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest["selected_counts"], indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument(
        "--output", type=Path, default=Path("data/transmission_projects/tyndp2026")
    )
    args = parser.parse_args()
    prepare(args.source, args.output)
