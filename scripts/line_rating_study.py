# SPDX-License-Identifier: MIT
"""Resolve explicit physical-rating sources without building a network."""

from pathlib import Path

import yaml


def resolve_study(path, planning_horizons):
    years = planning_horizons if isinstance(planning_horizons, list) else [planning_horizons]
    if len(years) != 1:
        raise ValueError("A physical nominal-rating study needs one planning year")
    year = int(years[0])
    study = yaml.safe_load(Path(path).read_text())
    if study.get("policy") != "fixed-snapshot-representative":
        raise ValueError("Only fixed physical nominal-rating snapshots are supported")
    choices = study.get("sources_by_year", {})
    ids = choices.get(year, choices.get(str(year)))
    if not isinstance(ids, list) or not ids or len(ids) != len(set(ids)):
        raise ValueError(f"Select distinct source editions explicitly for model year {year}")
    registry = yaml.safe_load(Path(study["registry"]).read_text())["sources"]
    specs = []
    for source_id in ids:
        if source_id not in registry:
            raise ValueError(f"Unregistered physical-rating source: {source_id}")
        spec = dict(registry[source_id], id=source_id)
        spec["local_path"] = spec.get("path", f"data/hindcast/line-ratings/sources/{source_id}.zip")
        specs.append(spec)
    return {"model_year": year, "study_path": str(path), "registry_path": study["registry"],
            "bindings_path": study["bindings"], "policy": study["policy"], "sources": specs}
