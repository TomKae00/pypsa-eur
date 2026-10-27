#!/usr/bin/env python3
"""Prepare physical AC rating comparisons. Never changes a network.

Sources are selected explicitly, independently of model year. The JAO adapter
preserves every current field; it does not turn seasonal/DLR values into a
constant rating. A normalized CSV adapter allows other TSO sources to use the
same matching and reporting code. Map-based matches are review candidates.
"""

import argparse
import hashlib
import html
import io
import json
import re
from pathlib import Path
from urllib.request import Request, urlopen
from zipfile import ZipFile, is_zipfile

import numpy as np
import pandas as pd
import yaml


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_component(ds, component, columns):
    return pd.DataFrame(
        {c: ds[f"{component}_{c}"][:] for c in columns
         if f"{component}_{c}" in ds.variables},
        index=pd.Index(ds[f"{component}_i"][:], name={"lines": "line_id", "links": "link_id"}.get(component, "bus")),
    )


def inventory(path):
    from netCDF4 import Dataset

    with Dataset(path) as ds:
        buses = read_component(ds, "buses", ["country", "x", "y"])
        lines = read_component(ds, "lines", [
            "bus0", "bus1", "v_nom", "s_nom", "num_parallel", "type",
            "length", "s_max_pu", "underground", "tags",
        ])
    return annotate_inventory(lines, buses)


def annotate_inventory(lines, buses):
    """Use identical nominal-input checks for saved and in-construction networks."""
    lines = lines.copy()
    if lines.index.has_duplicates:
        raise ValueError("Original line identifiers must be unique")
    for end in [0, 1]:
        for attr in ["country", "x", "y"]:
            lines[f"{attr}{end}"] = lines[f"bus{end}"].map(buses[attr])
    lines["model_current_a_per_circuit"] = (
        1000 * lines.s_nom / (np.sqrt(3) * lines.v_nom * lines.num_parallel)
    )
    # These are review flags, never automatic voltage/type/capacity replacements.
    flags = pd.Series("", index=lines.index, dtype=object)
    checks = {
        "nonstandard_nominal_voltage": ~lines.v_nom.isin([220, 225, 275, 300, 330, 380, 400, 420, 500, 750]),
        "invalid_nominal_capacity": ~np.isfinite(lines.s_nom) | lines.s_nom.le(0),
        "invalid_circuit_count": ~np.isfinite(lines.num_parallel) | lines.num_parallel.le(0),
        "fractional_original_circuit_count": ~np.isclose(lines.num_parallel, lines.num_parallel.round()),
    }
    if "underground" in lines:
        checks["generic_overhead_type_on_underground_asset"] = (
            lines.underground.eq(True) & lines.type.str.contains("Al/St", regex=False, na=False)
        )
    for label, mask in checks.items():
        flags.loc[mask] = flags.loc[mask] + label + ";"
    lines["nominal_input_flags"] = flags.str.rstrip(";")
    return lines


def nominal_link_inventory(path):
    """Keep physical DC p_nom separate from AC MVA and model connectors."""
    from netCDF4 import Dataset

    with Dataset(path) as ds:
        buses = read_component(ds, "buses", ["country"])
        if "links_i" not in ds.variables:
            return pd.DataFrame()
        links = read_component(ds, "links", ["bus0", "bus1", "carrier", "p_nom", "voltage", "tags"])
    if "carrier" not in links:
        links["carrier"] = ""
    if "p_nom" not in links:
        links["p_nom"] = np.nan
    for end in [0, 1]:
        links[f"country{end}"] = links[f"bus{end}"].map(buses.country)
    links["inventory_role"] = np.where(links.carrier.eq("DC"), "dc_branch_record", "model_connector_or_other_link")
    links = links.rename(columns={"p_nom": "nominal_power_mw"})
    links["reference_status"] = "not_reviewed_against_physical_project_reference"
    return links


def cached_source(spec, cache, offline):
    if "path" in spec:
        path = Path(spec["path"])
    else:
        path = cache / (spec["id"] + ".zip")
        if not path.exists():
            if offline:
                raise FileNotFoundError(f"Offline source is missing: {path}")
            cache.mkdir(parents=True, exist_ok=True)
            request = Request(spec["url"], headers={"User-Agent": "PyPSA-Eur physical-rating audit"})
            with urlopen(request, timeout=60) as response:
                content = response.read()
            temporary = path.with_suffix(".download")
            temporary.write_bytes(content)
            temporary.replace(path)
    digest = sha256(path)
    if spec.get("sha256") and spec["sha256"] != digest:
        raise ValueError(f"Source checksum differs for {spec['id']}")
    return path, digest


def map_coordinates(archive):
    """Read approximate endpoints from a bundled JAO Leaflet map, if present.

    JAO displaces some parallel lines for display. Coordinates are candidate
    evidence, not proof of circuit identity. No external map requests are made.
    """
    maps = []
    for name in archive.namelist():
        if name.lower().endswith(".html"):
            maps.append(archive.read(name))
        elif name.lower().endswith(".zip") and not name.startswith("__MACOSX"):
            nested_bytes = io.BytesIO(archive.read(name))
            if not is_zipfile(nested_bytes):
                continue  # Optional map member; absence remains visible as unlocated rows.
            with ZipFile(nested_bytes) as nested:
                maps.extend(nested.read(n) for n in nested.namelist() if n.lower().endswith(".html"))
    rows = []
    for raw in maps:
        for payload in re.findall(r'<script[^>]*type="application/json"[^>]*>(.*?)</script>', raw.decode("utf-8"), re.S):
            data = json.loads(payload)
            for call in data.get("x", {}).get("calls", []):
                if call.get("method") != "addPolylines":
                    continue
                args = call["args"]
                popups = args[6] if len(args) > 6 and isinstance(args[6], list) else args[4]
                if not isinstance(popups, list):
                    continue
                for geometry, popup in zip(args[0], popups):
                    name = re.search(r"NE name:\s*(.*?)</b>", popup)
                    eic = re.search(r"EIC Code:<b>\s*([^<]+)", popup)
                    def points(obj):
                        if isinstance(obj, dict) and {"lng", "lat"} <= obj.keys():
                            return list(zip(obj["lng"], obj["lat"]))
                        return [p for item in obj for p in points(item)] if isinstance(obj, list) else []
                    coords = points(geometry)
                    if name and eic and len(coords) >= 2:
                        rows.append({"name": html.unescape(name[1]).strip(), "eic": eic[1].strip(),
                                     "lon0": coords[0][0], "lat0": coords[0][1],
                                     "lon1": coords[-1][0], "lat1": coords[-1][1]})
    return pd.DataFrame(rows, columns=["name", "eic", "lon0", "lat0", "lon1", "lat1"]).drop_duplicates()


def read_jao(path, spec):
    rows = []
    with ZipFile(path) as archive:
        available = [n for n in archive.namelist() if n.lower().endswith(".xlsx") and not n.startswith("__MACOSX")]
        workbooks = spec.get("workbooks", [n for n in available if n not in spec.get("exclude_workbooks", {})])
        if not workbooks or set(workbooks) - set(available):
            raise ValueError("Selected JAO workbooks are absent")
        for name in workbooks:
            excel = pd.ExcelFile(io.BytesIO(archive.read(name)))
            for sheet in ["Lines", "Tielines"]:
                frame = pd.read_excel(excel, sheet_name=sheet, header=1)
                frame.columns = frame.columns.str.strip()
                frame = frame.dropna(how="all")
                expected = ["NE_name", "EIC_Code", "TSO", "Full_name", "Full_name.1", "Voltage_level(kV)"]
                if not set(expected) <= set(frame.columns):
                    raise ValueError(f"Unrecognized JAO identity columns in {name}/{sheet}")
                required = ["Period 1", "Period 2", "Period 3", "Period 4", "Period 5", "Period 6", "Fixed", "DLRmin(A)", "DLRmax(A)"]
                if not set(required) <= set(frame.columns):
                    raise ValueError(f"Unrecognized JAO rating columns in {name}/{sheet}")
                for idx, r in frame.iterrows():
                    row = {"reference_id": f"{name}:{sheet}:{idx + 3}", "sheet": sheet, "excel_row": idx + 3,
                           "name": r["NE_name"], "eic": r["EIC_Code"], "tso": r["TSO"],
                           "station0": r["Full_name"], "station1": r["Full_name.1"],
                           "voltage_kv": r["Voltage_level(kV)"], "reference_length_km": r["Length_(km)"],
                           "fixed_current_a": r["Fixed"], "dlr_min_a": r["DLRmin(A)"], "dlr_max_a": r["DLRmax(A)"],
                           "source_workbook": name}
                    row.update({f"period_{i}_a": r[f"Period {i}"] for i in range(1, 7)})
                    rows.append(row)
        coords = map_coordinates(archive)
    ref = pd.DataFrame(rows)
    ref["eic"] = ref.eic.fillna("").str.strip()
    ref["name"] = ref.name.fillna("").str.strip()
    counts = coords.groupby(["name", "eic"]).size().rename("map_geometry_count")
    coords = coords.drop_duplicates(["name", "eic"])
    ref = ref.merge(coords, on=["name", "eic"], how="left").merge(counts, on=["name", "eic"], how="left")
    return ref


def normalize(spec, path):
    if spec["format"] == "jao_core_sgm_zip":
        ref = read_jao(path, spec)
    elif spec["format"] == "normalized_csv":
        ref = pd.read_csv(path)
        required = {"reference_id", "name", "voltage_kv", "lon0", "lat0", "lon1", "lat1", "fixed_current_a"}
        if not required <= set(ref):
            raise ValueError(f"Normalized source lacks {sorted(required - set(ref))}")
        for c in ["reference_length_km", "dlr_min_a", "dlr_max_a", "map_geometry_count"] + [f"period_{i}_a" for i in range(1, 7)]:
            if c not in ref:
                ref[c] = np.nan
        for c in ["tso", "eic", "station0", "station1", "source_workbook"]:
            if c not in ref:
                ref[c] = ""
    else:
        raise ValueError(f"Unsupported source format: {spec['format']}")
    if ref.reference_id.isna().any() or ref.reference_id.duplicated().any():
        raise ValueError("Reference IDs must be present and unique within an edition")
    numeric = ["voltage_kv", "fixed_current_a", "dlr_min_a", "dlr_max_a", "reference_length_km",
               "lon0", "lat0", "lon1", "lat1"] + [f"period_{i}_a" for i in range(1, 7)]
    parse_issues = {idx: {} for idx in ref.index}
    for c in numeric:
        # The JAO workbook uses a dash for unavailable data, never for zero.
        def clean_number(value):
            if not isinstance(value, str):
                return value
            value = value.strip()
            if value in {"-", ""}:
                return np.nan
            if re.fullmatch(r"[+-]?\d+,\d+", value):
                value = value.replace(",", ".")
            return value
        cleaned = ref[c].map(clean_number)
        parsed = pd.to_numeric(cleaned, errors="coerce")
        invalid = cleaned.notna() & (parsed.isna() | ~np.isfinite(parsed))
        for idx in ref.index[invalid]:
            parse_issues[idx][c] = str(ref.at[idx, c])
        ref[c] = parsed.replace([np.inf, -np.inf], np.nan)
    ref["source_parse_issues"] = [json.dumps(parse_issues[idx], ensure_ascii=False) for idx in ref.index]
    currents = ref[["fixed_current_a", "dlr_min_a", "dlr_max_a"] + [f"period_{i}_a" for i in range(1, 7)]]
    if (currents.dropna(how="all") < 0).any().any():
        raise ValueError("Negative current in source")
    ref["source_id"] = spec["id"]
    # Other TSOs need not publish EICs. Preserve their actual circuit codes or
    # explicit reference IDs instead of inventing an EIC for a local identifier.
    if "circuit_id" not in ref:
        eic = ref["eic"].fillna("").astype(str).str.strip()
        ref["circuit_id"] = eic.where(eic.ne(""), ref.reference_id.astype(str))
    if ref.circuit_id.isna().any() or ref.circuit_id.astype(str).str.strip().eq("").any():
        raise ValueError("Circuit identity must be explicit in every reference row")
    ref["circuit_id"] = ref.circuit_id.astype(str)
    ref["source_url"] = spec.get("url", spec.get("provenance_url", ""))
    ref["snapshot_date"] = spec.get("snapshot_date", "")
    ref["temporal_scope"] = "snapshot; not evidence of availability or whole-year validity"
    has_period = ref[[f"period_{i}_a" for i in range(1, 7)]].notna().any(axis=1)
    has_dlr = ref[["dlr_min_a", "dlr_max_a"]].notna().any(axis=1)
    ref["rating_kind"] = np.select(
        [ref.source_parse_issues.ne("{}"), ref.fixed_current_a.gt(0) & ~has_period & ~has_dlr, has_dlr, has_period],
        ["source_data_issue", "fixed", "dynamic_or_mixed", "periodic_or_mixed"], default="missing_or_zero")
    return ref


def distance_km(lon, lat, x, y):
    lon, lat, x, y = [np.radians(np.asarray(v, dtype=float)) for v in [lon, lat, x, y]]
    a = np.sin((lat-y)/2)**2 + np.cos(lat)*np.cos(y)*np.sin((lon-x)/2)**2
    return 6371.0088 * 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def candidates(lines, reference, tolerance_km):
    located = reference.dropna(subset=["lon0", "lat0", "lon1", "lat1"])
    result = []
    for line_id, line in lines.iterrows():
        # Accept 220/225 and 380/400 labels as candidates, but report mismatch.
        ref = located.loc[(located.voltage_kv-line.v_nom).abs() <= 0.06*line.v_nom]
        if ref.empty:
            continue
        forward = np.maximum(distance_km(ref.lon0, ref.lat0, line.x0, line.y0),
                             distance_km(ref.lon1, ref.lat1, line.x1, line.y1))
        reverse = np.maximum(distance_km(ref.lon1, ref.lat1, line.x0, line.y0),
                             distance_km(ref.lon0, ref.lat0, line.x1, line.y1))
        d = np.asarray(np.minimum(forward, reverse))
        for idx in ref.index[d <= tolerance_km]:
            r = ref.loc[idx]
            result.append({"line_id": line_id, "source_id": r.source_id, "reference_id": r.reference_id,
                           "endpoint_max_distance_km": d[ref.index.get_loc(idx)],
                           "voltage_label_differs": line.v_nom != r.voltage_kv,
                           "model_circuits": line.num_parallel,
                           "model_current_a_per_circuit": line.model_current_a_per_circuit,
                           "reference_fixed_current_a": r.fixed_current_a,
                           "fixed_to_model_current_ratio": r.fixed_current_a/line.model_current_a_per_circuit,
                           "model_to_reference_length_ratio": line.length/r.reference_length_km if r.reference_length_km > 0 else np.nan,
                           "rating_kind": r.rating_kind,
                           "status": "endpoint_candidate_requires_circuit_review"})
    columns = ["line_id", "source_id", "reference_id", "endpoint_max_distance_km", "voltage_label_differs",
               "model_circuits", "model_current_a_per_circuit", "reference_fixed_current_a",
               "fixed_to_model_current_ratio", "model_to_reference_length_ratio", "rating_kind", "status"]
    result = pd.DataFrame(result, columns=columns)
    if not result.empty:
        result["candidate_count_for_line_and_source"] = result.groupby(["line_id", "source_id"]).line_id.transform("size")
        result = result.sort_values(["line_id", "source_id", "endpoint_max_distance_km", "reference_id"])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("network", type=Path, help="Original base/base_extended network, before 380-kV simplification")
    parser.add_argument("--year", required=True, type=int, help="Model year; independent of source edition")
    parser.add_argument("--source-config", type=Path, default=Path("config/hindcast/line-rating-sources.yaml"))
    parser.add_argument("--source-id", action="append", required=True, help="Explicit edition ID; repeat to compare editions")
    parser.add_argument("--cache", type=Path, default=Path("data/hindcast/line-ratings/sources"))
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--endpoint-tolerance-km", type=float, default=5.0)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--clustered-network", type=Path)
    parser.add_argument("--linemap", type=Path)
    args = parser.parse_args()
    if args.endpoint_tolerance_km <= 0 or not np.isfinite(args.endpoint_tolerance_km):
        parser.error("Endpoint tolerance must be finite and positive")
    if len(set(args.source_id)) != len(args.source_id):
        parser.error("Select each source edition only once")
    if bool(args.clustered_network) != bool(args.linemap):
        parser.error("Supply --clustered-network and --linemap together")
    registry = yaml.safe_load(args.source_config.read_text())["sources"]
    lines = inventory(args.network)
    links = nominal_link_inventory(args.network)
    if lines.v_nom.nunique() < 2:
        raise ValueError("Network appears simplified; provide the original-voltage network")
    args.output.mkdir(parents=True, exist_ok=True)
    references, sources = [], []
    for source_id in args.source_id:
        spec = dict(registry[source_id], id=source_id)
        path, digest = cached_source(spec, args.cache, args.offline)
        ref = normalize(spec, path)
        references.append(ref)
        sources.append({**spec, "cached_path": str(path), "download_sha256": digest,
                        "rows": len(ref), "located_rows": int(ref.lon0.notna().sum())})
    reference = pd.concat(references, ignore_index=True)
    match = candidates(lines, reference, args.endpoint_tolerance_km)
    lines.to_csv(args.output / "original_line_inventory.csv")
    links.to_csv(args.output / "original_link_nominal_inventory.csv")
    reference.to_csv(args.output / "reference_ratings.csv", index=False)
    match.to_csv(args.output / "candidate_matches.csv", index=False)
    aggregation = None
    if args.clustered_network:
        from netCDF4 import Dataset

        with Dataset(args.clustered_network) as ds:
            clustered = read_component(ds, "lines", ["s_nom", "bus0", "bus1"])
        linemap = pd.read_csv(args.linemap, index_col=0).iloc[:, 0].astype(str)
        lines["clustered_line_id"] = linemap.reindex(lines.index).where(
            linemap.reindex(lines.index).isin(clustered.index.astype(str)))
        aggregates = []
        for line_id, row in clustered.iterrows():
            members = linemap.index[linemap.eq(str(line_id))]
            missing = members.difference(lines.index)
            original_sum = lines.s_nom.reindex(members).sum() if not len(missing) else np.nan
            aggregates.append({"clustered_line_id": line_id, "source_records": len(members),
                               "missing_source_records": len(missing), "original_sum_mva": original_sum,
                               "clustered_s_nom_mva": row.s_nom,
                               "matches": bool(len(members) and not len(missing) and np.isclose(original_sum, row.s_nom, rtol=1e-8, atol=1e-6))})
        aggregate = pd.DataFrame(aggregates)
        aggregate.to_csv(args.output / "aggregation_check.csv", index=False)
        aggregation = {"lines": len(aggregate), "matching_nominal_sums": int(aggregate.matches.sum()),
                       "network_sha256": sha256(args.clustered_network), "linemap_sha256": sha256(args.linemap)}
    coverage = []
    matched = set(match.line_id)
    for country in sorted(set(lines.country0) | set(lines.country1)):
        subset = lines.loc[lines.country0.eq(country) | lines.country1.eq(country)]
        coverage.append({"country": country, "incident_original_line_records": len(subset),
                         "records_with_candidates": len(subset.index.intersection(matched)),
                         "records_without_candidates": len(subset.index.difference(matched)),
                         "approved_corrections": 0})
    pd.DataFrame(coverage).to_csv(args.output / "country_coverage.csv", index=False)
    lines.to_csv(args.output / "original_line_inventory.csv")
    reference.groupby(["source_id", "tso", "rating_kind"], dropna=False).size().rename("reference_records").reset_index().to_csv(
        args.output / "reference_coverage.csv", index=False)
    summary = {"model_year": args.year, "network": str(args.network), "network_sha256": sha256(args.network),
               "source_config_sha256": sha256(args.source_config), "sources": sources,
               "original_lines": len(lines), "countries": sorted(set(lines.country0) | set(lines.country1)),
               "link_inventory_role_counts": {str(k): int(v) for k, v in links.get("inventory_role", pd.Series(dtype=str)).value_counts().items()},
               "candidate_rows": len(match), "lines_with_candidates": int(match.line_id.nunique()),
               "approved_corrections": 0, "network_modified": False,
               "endpoint_tolerance_km": args.endpoint_tolerance_km,
               "source_selection": "explicit; source edition does not imply whole-year applicability",
               "scope": "physical AC rating comparison; no market NTC, commissioning-date or outage adjustments",
               "aggregation_check": aggregation}
    (args.output / "audit.json").write_text(json.dumps(summary, indent=2) + "\n")
    review = match.merge(reference.drop(columns="rating_kind"), on=["source_id", "reference_id"])
    review = review.merge(lines[[c for c in ["country0", "country1", "clustered_line_id"] if c in lines]],
                          left_on="line_id", right_index=True, how="left")
    review["priority_for_review"] = (
        review.rating_kind.eq("fixed") & review.model_circuits.eq(1)
        & review.candidate_count_for_line_and_source.eq(1)
        & ~review.voltage_label_differs
        & review.endpoint_max_distance_km.le(1)
        & review.model_to_reference_length_ratio.between(0.8, 1.25)
        & review.map_geometry_count.eq(1)
    )
    review["approved"] = False
    for column in ["reviewer", "review_note", "valid_from", "valid_to", "approved_s_nom_mva", "allocation_basis", "rating_basis"]:
        review[column] = ""
    review.to_csv(args.output / "rating_review.csv", index=False)
    if "clustered_line_id" in review:
        review.loc[review.clustered_line_id.notna()].to_csv(args.output / "interzone_rating_review.csv", index=False)
    print(f"Inventoried {len(lines)} lines; {len(match)} candidate pairs for {match.line_id.nunique()} lines.")
    print("No corrections approved or applied. Inspect rating_review.csv and audit.json.")


if __name__ == "__main__":
    main()
