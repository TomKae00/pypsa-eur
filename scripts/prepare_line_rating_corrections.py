#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compile reviewed circuit bindings and a selected source into native overrides.

No network is changed. A fixed source snapshot is explicitly adopted as a
representative annual rating; availability and commissioning dates are outside
this policy. Only complete, agreeing, fixed-current EIC groups are accepted.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from audit_line_ratings import cached_source, inventory, normalize, sha256


def compile_corrections(lines, reference, bindings):
    guards, updates, decisions = {}, {}, []
    used_lines, used_eics = set(), set()
    for binding_id, binding in bindings.items():
        ids = list(binding["lines"])
        eics = binding["eics"]
        if used_lines.intersection(ids) or used_eics.intersection(eics):
            raise ValueError(f"Overlapping circuit identity bindings: {binding_id}")
        used_lines.update(ids)
        used_eics.update(eics)
        record = {"binding": binding_id, "line_ids": ids, "eics": eics,
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
        selected = reference.loc[reference.eic.isin(eics)]
        record["reference_ids"] = selected.reference_id.tolist()
        record["reported_currents_a"] = sorted(selected.fixed_current_a.dropna().unique().tolist())
        record["reporting_tsos"] = sorted(selected.tso.dropna().unique().tolist())
        complete = all(selected.loc[selected.eic.eq(eic), "tso"].nunique() >= 2 for eic in eics)
        if not complete:
            record.update(status="hold", reason="Each circuit needs reports from both TSOs")
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
                          reason="Reviewed circuit group; both TSOs agree on one positive fixed current")
        decisions.append(record)
    return updates, guards, decisions


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
    if not updates:
        raise ValueError("No supported fixed ratings; no correction file produced")
    metadata = {"schema_version": 1, "model_year": args.year, "source_id": args.source_id,
                "source_url": spec.get("url", spec.get("provenance_url", "")),
                "source_sha256": digest, "bindings_sha256": sha256(args.bindings),
                "audit_network_sha256": sha256(args.network), "rating_policy": args.rating_policy,
                "coverage": "partial; unreviewed lines retain native assumptions",
                "scope": "physical current limits; no market transfer limits or availability corrections"}
    payload = {"metadata": metadata, "expected_lines": guards,
               "corrections": {"Line": {"s_nom": {"index": updates}}}, "decisions": decisions}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(yaml.safe_dump(payload, sort_keys=False, allow_unicode=True))
    args.output.with_suffix(".audit.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    print(f"Prepared {len(updates)} original-line corrections for model year {args.year} from {args.source_id}.")
    for decision in decisions:
        print(f"{decision['binding']}: {decision['status']} — {decision['reason']}")


if __name__ == "__main__":
    main()
