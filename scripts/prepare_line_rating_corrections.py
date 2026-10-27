#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compile reviewed circuit bindings and a selected source into native overrides.

No network is changed. A fixed source snapshot is explicitly adopted as a
representative annual rating; availability and commissioning dates are outside
this policy. Only reviewed, fixed-current EIC groups are accepted. Border
circuits require agreeing reports from both TSOs; internal circuits use their owner.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

if __package__:
    from .audit_line_ratings import cached_source, inventory, normalize, sha256
    from .line_rating_study import resolve_study
else:
    from audit_line_ratings import cached_source, inventory, normalize, sha256
    from line_rating_study import resolve_study


def compile_corrections(lines, reference, bindings):
    guards, updates, decisions = {}, {}, []
    used_lines, used_eics = set(), set()
    for binding_id, binding in bindings.items():
        ids = list(binding["lines"])
        if bool(binding.get("eics")) == bool(binding.get("circuit_ids")):
            raise ValueError(f"Binding {binding_id} must select either EICs or source circuit IDs")
        identity_column = "eic" if binding.get("eics") else "circuit_id"
        eics = binding.get("eics", binding.get("circuit_ids"))
        if any(not isinstance(key, str) or not key.strip() for key in eics):
            raise ValueError(f"Physical circuit IDs in {binding_id} must be nonempty strings")
        if len(eics) != len(set(eics)):
            raise ValueError(f"Duplicated physical circuit identities in {binding_id}")
        if used_lines.intersection(ids) or used_eics.intersection(eics):
            raise ValueError(f"Overlapping circuit identity bindings: {binding_id}")
        used_lines.update(ids)
        used_eics.update(eics)
        record = {"binding": binding_id, "line_ids": ids,
                  "eics" if identity_column == "eic" else "circuit_ids": eics,
                  "evidence": binding["evidence"], "evidence_url": binding["evidence_url"]}
        missing = set(ids) - set(lines.index)
        if missing:
            record.update(status="hold", reason=f"Original model IDs absent: {sorted(missing)}")
            decisions.append(record)
            continue
        # A circuit binding must be re-reviewed when original endpoints or counts change.
        for line_id, expected in binding["lines"].items():
            for attr, value in expected.items():
                actual = lines.at[line_id, attr]
                equal = np.isclose(actual, value, rtol=1e-9, atol=1e-9) if isinstance(value, (int, float)) else actual == value
                if not equal:
                    raise ValueError(f"Binding {binding_id}: {line_id} {attr} changed ({actual!r} != {value!r})")
        if not np.isclose(lines.loc[ids, "num_parallel"].sum(), len(set(eics))):
            raise ValueError(f"Binding {binding_id} does not reconcile the number of circuits")
        if lines.loc[ids, "v_nom"].nunique() != 1:
            raise ValueError(f"Binding {binding_id} must represent parallel circuits at one nominal voltage")
        selected = reference.loc[reference[identity_column].isin(eics)]
        record["reference_ids"] = selected.reference_id.tolist()
        record["reported_currents_a"] = sorted(selected.fixed_current_a.dropna().unique().tolist())
        record["reporting_tsos"] = sorted(selected.tso.dropna().unique().tolist())
        record["reference_circuit_names"] = selected["name"].tolist()
        record["reference_station_pairs"] = selected[["station0", "station1"]].fillna("").values.tolist()
        owner = binding.get("owning_tso")
        if owner:
            country = binding.get("country")
            if not country or not (
                lines.loc[ids, "country0"].eq(country).all()
                and lines.loc[ids, "country1"].eq(country).all()
            ):
                raise ValueError(f"Single-owner binding {binding_id} must be internal to its declared country")
            complete = all(
                set(selected.loc[selected[identity_column].eq(eic), "tso"].dropna().str.strip().str.upper())
                == {owner.strip().upper()} for eic in eics
            )
            record["rating_authority"] = f"Internal circuit; owning TSO {owner}"
        else:
            complete = all(selected.loc[selected[identity_column].eq(eic), "tso"].nunique() >= 2 for eic in eics)
            record["rating_authority"] = "Cross-border circuit; agreeing reports from both TSOs"
        if not complete:
            record.update(status="hold", reason=(
                "Each internal circuit needs its declared owner and no unexpected reporting TSO"
                if owner else "Each border circuit needs reports from both TSOs"
            ))
        elif selected.voltage_kv.isna().any() or not (
            (selected.voltage_kv - lines.loc[ids, "v_nom"].iloc[0]).abs()
            <= 0.06 * lines.loc[ids, "v_nom"].iloc[0]
        ).all():
            record.update(status="hold", reason="Source nominal voltage differs from the reviewed physical circuit")
        elif not selected.rating_kind.eq("fixed").all():
            record.update(status="hold", reason="Source has a seasonal, dynamic, missing or malformed rating")
        elif selected.fixed_current_a.isna().any() or not selected.fixed_current_a.gt(0).all():
            record.update(status="hold", reason="Source has missing or non-positive fixed current")
        elif not np.allclose(selected.fixed_current_a, selected.fixed_current_a.iloc[0], rtol=0, atol=1e-6):
            record.update(status="hold", reason="Fixed currents differ between TSOs or parallel circuits; no averaging/minimum inferred")
        else:
            current = float(selected.fixed_current_a.iloc[0])
            for line_id in ids:
                old = lines.loc[line_id]
                updates[line_id] = float(np.sqrt(3) * old.v_nom * current / 1000 * old.num_parallel)
                guards[line_id] = {k: str(old[k]) for k in ["bus0", "bus1", "type"]}
                guards[line_id].update({k: float(old[k]) for k in ["v_nom", "num_parallel", "s_nom"]})
            record.update(status="accepted", current_a_per_circuit=current,
                          reason=("Reviewed internal circuit; positive fixed rating from its owning TSO"
                                  if owner else "Reviewed border group; both TSOs agree on one positive fixed current"))
        decisions.append(record)
    return updates, guards, decisions


def nominal_coverage(lines, updates, decisions):
    """Account for every original line, including all retained assumptions."""
    rows = lines[["country0", "country1", "v_nom", "num_parallel", "s_nom"]].copy()
    if "nominal_input_flags" in lines:
        rows["nominal_input_flags"] = lines.nominal_input_flags
    rows = rows.rename(columns={"s_nom": "original_s_nom_mva"})
    rows["revised_s_nom_mva"] = rows.index.to_series().map(updates).fillna(rows.original_s_nom_mva)
    rows["status"] = "native_assumption_unreviewed"
    rows["binding"] = ""
    rows["reason"] = "No reviewed physical reference; retain native nominal rating"
    for decision in decisions:
        ids = rows.index.intersection(decision["line_ids"])
        rows.loc[ids, "binding"] = decision["binding"]
        rows.loc[ids, "reason"] = decision["reason"]
        if decision["status"] == "accepted":
            changed = ~np.isclose(rows.loc[ids, "original_s_nom_mva"], rows.loc[ids, "revised_s_nom_mva"], rtol=1e-8, atol=1e-6)
            rows.loc[ids[changed], "status"] = "nominal_rating_corrected"
            rows.loc[ids[~changed], "status"] = "nominal_rating_verified_unchanged"
        else:
            rows.loc[ids, "status"] = "reviewed_binding_on_hold"
    countries = []
    for country in sorted(set(rows.country0) | set(rows.country1)):
        incident = rows.loc[rows.country0.eq(country) | rows.country1.eq(country)]
        counts = incident.status.value_counts().to_dict()
        countries.append({"country": country, "incident_original_records": len(incident),
                          **{key: int(counts.get(key, 0)) for key in [
                              "nominal_rating_corrected", "nominal_rating_verified_unchanged",
                              "reviewed_binding_on_hold", "native_assumption_unreviewed"]}})
    return rows, pd.DataFrame(countries)


def nominal_flag_counts(lines):
    """Review indicators, not evidence that an individual rating is incorrect."""
    if "nominal_input_flags" not in lines:
        return {}
    flags = lines.nominal_input_flags.str.split(";").explode()
    return {str(k): int(v) for k, v in flags[flags.ne("")].value_counts().items()}


def payload_from_study(lines, study_path, planning_horizons, source_files):
    """Compile within native base construction, avoiding a baseline-file cycle."""
    selection = resolve_study(study_path, planning_horizons)
    if isinstance(source_files, (str, Path)):
        source_files = [source_files]
    if len(source_files) != len(selection["sources"]):
        raise ValueError("Declared source inputs do not match the selected editions")
    frames, sources = [], []
    for spec, path in zip(selection["sources"], source_files):
        digest = sha256(path)
        if spec.get("sha256") and digest != spec["sha256"]:
            raise ValueError(f"Physical-rating checksum differs: {spec['id']}")
        frames.append(normalize(spec, path))
        sources.append({"id": spec["id"], "sha256": digest,
                        "url": spec.get("url", spec.get("provenance_url", ""))})
    reference = pd.concat(frames, ignore_index=True)
    bindings = yaml.safe_load(Path(selection["bindings_path"]).read_text())["bindings"]
    updates, guards, decisions = compile_corrections(lines, reference, bindings)
    if not updates:
        raise ValueError("The selected study has no supported physical nominal ratings")
    coverage, countries = nominal_coverage(lines, updates, decisions)
    source_digest = (sources[0]["sha256"] if len(sources) == 1 else
                     hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest())
    metadata = {"schema_version": 1, "model_year": selection["model_year"],
                "source_id": "+".join(s["id"] for s in sources),
                "source_url": " | ".join(s["url"] for s in sources), "source_sha256": source_digest,
                "source_hash_kind": "archive_sha256" if len(sources) == 1 else "source_manifest_sha256",
                "sources": sources, "bindings_sha256": sha256(selection["bindings_path"]),
                "study_sha256": sha256(study_path), "registry_sha256": sha256(selection["registry_path"]),
                "rating_policy": selection["policy"],
                "coverage": "partial; all unreviewed lines retain native nominal assumptions",
                "original_line_status_counts": {str(k): int(v) for k, v in coverage.status.value_counts().items()},
                "nominal_input_review_flags": nominal_flag_counts(lines),
                "country_coverage": countries.to_dict("records"),
                "scope": "physical nominal ratings only; no operating or market capacity adjustments"}
    return {"metadata": metadata, "expected_lines": guards,
            "corrections": {"Line": {"s_nom": {"index": updates}}}, "decisions": decisions}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("network", type=Path, help="Original-voltage base/base_extended network")
    parser.add_argument("--year", type=int, required=True, help="Model/planning year, independent of weather or source edition")
    parser.add_argument("--source-id", required=True)
    parser.add_argument("--rating-policy", required=True, choices=["fixed-snapshot-representative"],
                        help="Explicitly adopt the selected fixed ratings as representative annual physical inputs")
    parser.add_argument("--source-config", type=Path, default=Path("config/hindcast/line-rating-sources.yaml"))
    parser.add_argument("--bindings", type=Path, default=Path("config/hindcast/line-rating-bindings.yaml"))
    parser.add_argument("--cache", type=Path, default=Path("data/hindcast/line-ratings/sources"))
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    registry = yaml.safe_load(args.source_config.read_text())["sources"]
    spec = dict(registry[args.source_id], id=args.source_id)
    source_path, digest = cached_source(spec, args.cache, args.offline)
    reference = normalize(spec, source_path)
    lines = inventory(args.network)
    bindings = yaml.safe_load(args.bindings.read_text())["bindings"]
    updates, guards, decisions = compile_corrections(lines, reference, bindings)
    coverage, countries = nominal_coverage(lines, updates, decisions)
    if not updates:
        raise ValueError("No supported fixed ratings; no correction file produced")
    metadata = {"schema_version": 1, "model_year": args.year, "source_id": args.source_id,
                "source_url": spec.get("url", spec.get("provenance_url", "")),
                "source_sha256": digest, "bindings_sha256": sha256(args.bindings),
                "audit_network_sha256": sha256(args.network), "rating_policy": args.rating_policy,
                "coverage": "partial; unreviewed lines retain native assumptions",
                "original_line_status_counts": {str(k): int(v) for k, v in coverage.status.value_counts().items()},
                "nominal_input_review_flags": nominal_flag_counts(lines),
                "country_coverage": countries.to_dict("records"),
                "scope": "physical current limits; no market transfer limits or availability corrections"}
    payload = {"metadata": metadata, "expected_lines": guards,
               "corrections": {"Line": {"s_nom": {"index": updates}}}, "decisions": decisions}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(yaml.safe_dump(payload, sort_keys=False, allow_unicode=True))
    args.output.with_suffix(".audit.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    coverage.to_csv(args.output.with_suffix(".nominal_coverage.csv"), index_label="line_id")
    countries.to_csv(args.output.with_suffix(".country_coverage.csv"), index=False)
    print(f"Prepared {len(updates)} original-line corrections for model year {args.year} from {args.source_id}.")
    for decision in decisions:
        print(f"{decision['binding']}: {decision['status']} — {decision['reason']}")


if __name__ == "__main__":
    main()
