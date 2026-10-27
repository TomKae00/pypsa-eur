#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Report Nordic physical nominal inputs and reviewed source evidence.

This is a read-only companion to audit_line_ratings, not another correction
mechanism. AC MVA, DC MW, published cable MW and market TTC remain distinct.
Source editions are recorded in the review manifest independently of model year.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

if __package__:
    from .audit_line_ratings import inventory, nominal_link_inventory, sha256
else:
    from audit_line_ratings import inventory, nominal_link_inventory, sha256


def read_busmap(path):
    frame = pd.read_csv(path, index_col=0, dtype=str, keep_default_na=False)
    if len(frame.columns) != 1 or frame.index.has_duplicates:
        raise ValueError(f"Expected one mapping column and unique source buses: {path}")
    return frame.iloc[:, 0]


def review_inputs(frame, component, manifest):
    """Compare reviewed identities before attaching any source finding."""
    result = frame.copy()
    result["asset_name"] = result.index.astype(str)
    result["review_status"] = "native_assumption_unreviewed"
    for field in ["source_ids", "source_urls", "identity_urls", "rating_basis", "review_note",
                  "observed_route_refs", "observed_bundle_tags"]:
        result[field] = ""
    result["reference_nominal_mw"] = np.nan
    result["published_cable_capacity_mw"] = np.nan
    decisions = []
    for review in manifest["reviews"]:
        if review["component"] != component:
            continue
        key = review["id"]
        decision = {"component": component, "id": key, "name": review["name"]}
        if key not in result.index:
            decisions.append({**decision, "status": "model_record_absent"})
            continue
        result.at[key, "asset_name"] = review["name"]
        expected = review["expected"]
        required = ({"bus0", "bus1", "type", "v_nom", "num_parallel"}
                    if component == "Line" else {"bus0", "bus1", "carrier", "voltage"})
        if set(expected) != required:
            raise ValueError(f"Incomplete identity signature: {component} {key}")
        changed = []
        for attr, value in expected.items():
            actual = result.at[key, attr]
            equal = (np.isfinite(actual) and np.isclose(actual, value, rtol=1e-8, atol=1e-6)
                     if isinstance(value, (int, float)) else actual == value)
            if not equal:
                changed.append(attr)
        sources = review["source_ids"]
        result.at[key, "source_ids"] = ";".join(sources)
        result.at[key, "source_urls"] = " | ".join(manifest["sources"][s]["url"] for s in sources)
        result.at[key, "identity_urls"] = " | ".join(review.get("identity_urls", []))
        result.at[key, "rating_basis"] = review.get("rating_basis", "continuous_AC_ampacity_not_established")
        result.at[key, "review_note"] = review["note"]
        if changed:
            status = "model_signature_changed"
            result.at[key, "review_note"] = f"Re-review changed fields {changed}. " + review["note"]
        else:
            status = review["finding"]
            if "reference_nominal_mw" in review:
                value = review["reference_nominal_mw"]
                if component != "Link" or isinstance(value, bool) or not np.isfinite(value) or value <= 0:
                    raise ValueError(f"Invalid DC reference power: {key}")
                result.at[key, "reference_nominal_mw"] = value
                status = ("reference_value_matches" if np.isclose(result.at[key, "nominal_power_mw"], value,
                                                                 rtol=1e-8, atol=1e-6)
                          else "reference_value_differs")
            if "published_cable_capacity_mw" in review:
                value = review["published_cable_capacity_mw"]
                if component != "Line" or isinstance(value, bool) or not np.isfinite(value) or value <= 0:
                    raise ValueError(f"Invalid published cable MW: {key}")
                # Deliberately no comparison with, or conversion to, s_nom in MVA.
                result.at[key, "published_cable_capacity_mw"] = value
            for attr in ["observed_route_refs", "observed_bundle_tags"]:
                result.at[key, attr] = " | ".join(review.get(attr, []))
        result.at[key, "review_status"] = status
        decisions.append({**decision, "status": status})
    return result, decisions


def build_report(network, review_path, year, output, simplify_map=None, cluster_map=None):
    manifest = yaml.safe_load(Path(review_path).read_text())
    if manifest.get("schema_version") != 1:
        raise ValueError("Unsupported Nordic review manifest")
    keys = [(r["component"], r["id"]) for r in manifest["reviews"]]
    if len(set(keys)) != len(keys) or any(c not in {"Line", "Link"} for c, _ in keys):
        raise ValueError("Review records must have unique Line/Link identities")
    if bool(simplify_map) != bool(cluster_map):
        raise ValueError("Supply both busmaps or neither")
    lines = inventory(network)
    if lines.v_nom.nunique() < 2:
        raise ValueError("Use the original-voltage base/base_extended network")
    links = nominal_link_inventory(network)
    countries = manifest["countries"]
    lines = lines.loc[lines.country0.isin(countries) | lines.country1.isin(countries)].copy()
    links = links.loc[links.country0.isin(countries) | links.country1.isin(countries)].copy()
    fingerprints = {"network_sha256": sha256(network), "review_manifest_sha256": sha256(review_path)}
    if simplify_map:
        simple, cluster = read_busmap(simplify_map), read_busmap(cluster_map)
        # Removed buses may map through an identity rather than a listed merge.
        mapping = simple.map(cluster).fillna(simple.where(simple.isin(cluster.values)))
        for end in [0, 1]:
            lines[f"zone{end}"] = lines[f"bus{end}"].map(mapping)
        lines["zone_mapping_complete"] = lines[["zone0", "zone1"]].notna().all(axis=1)
        lines["interzone"] = lines.zone_mapping_complete & lines.zone0.ne(lines.zone1)
        fingerprints.update(simplify_busmap_sha256=sha256(simplify_map),
                            cluster_busmap_sha256=sha256(cluster_map))
    lines, ac_decisions = review_inputs(lines, "Line", manifest)
    dc = links.loc[links.inventory_role.eq("dc_branch_record")].copy()
    dc, dc_decisions = review_inputs(dc, "Link", manifest)
    connectors = links.loc[~links.inventory_role.eq("dc_branch_record")].copy()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    lines.to_csv(output / "nordic_ac_nominal_inventory.csv")
    dc.to_csv(output / "nordic_dc_nominal_review.csv")
    connectors.to_csv(output / "nordic_model_connectors.csv")
    if simplify_map:
        borders = lines.loc[lines.interzone].copy()
        borders["zone_pair"] = ["--".join(sorted([x, y])) for x, y in zip(borders.zone0, borders.zone1)]
        corridor = borders.groupby("zone_pair").agg(
            original_records=("s_nom", "size"), circuits=("num_parallel", "sum"),
            sum_nominal_mva=("s_nom", "sum"))
        corridor.to_csv(output / "nordic_ac_corridor_nominal_sums.csv")
        borders.to_csv(output / "nordic_ac_interzone_review.csv")
        # Specific missing fields, ready to complete from an asset-owner source.
        request = borders[["asset_name", "country0", "country1", "zone0", "zone1", "bus0", "bus1",
                           "v_nom", "num_parallel", "type", "s_nom", "model_current_a_per_circuit",
                           "observed_route_refs", "observed_bundle_tags", "review_status"]].copy()
        for field in ["tso_circuit_id", "continuous_current_a_per_circuit", "continuous_rating_mva",
                      "rating_conditions", "limiting_component", "valid_from", "valid_to", "reference_url"]:
            request[field] = ""
        request.to_csv(output / "nordic_ac_rating_data_request.csv")
    report = {"model_year": year, "reviewed_on": manifest["reviewed_on"], **fingerprints,
              "scope": manifest["scope"], "temporal_policy": manifest["temporal_policy"],
              "network_modified": False, "nominal_corrections_applied": 0,
              "ac_records": len(lines), "dc_branch_records": len(dc),
              "other_link_or_converter_records": len(connectors),
              "ac_review_status_counts": {str(k): int(v) for k, v in lines.review_status.value_counts().items()},
              "dc_review_status_counts": {str(k): int(v) for k, v in dc.review_status.value_counts().items()},
              "sources": manifest["sources"], "decisions": ac_decisions + dc_decisions,
              "interpretation": "Matching DC reference values is not certification of historical availability; AC corridor sums are MVA, not market transfer capacities."}
    if simplify_map:
        report.update(ac_interzone_records=int(lines.interzone.sum()),
                      ac_unmapped_records=int((~lines.zone_mapping_complete).sum()))
    (output / "nordic_nominal_review.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("network", type=Path)
    parser.add_argument("--year", required=True, type=int, help="Reporting/model year, independent of source edition")
    parser.add_argument("--review-config", type=Path, default=Path("config/hindcast/nordic-rating-review.yaml"))
    parser.add_argument("--simplify-busmap", type=Path)
    parser.add_argument("--cluster-busmap", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = build_report(args.network, args.review_config, args.year, args.output,
                          args.simplify_busmap, args.cluster_busmap)
    print(f"Reviewed {result['ac_records']} Nordic AC records and {result['dc_branch_records']} DC branch records.")
    print("DC reference comparisons:", result["dc_review_status_counts"])
    print("No nominal ratings changed. Results:", args.output)


if __name__ == "__main__":
    main()
