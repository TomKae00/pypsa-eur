# SPDX-License-Identifier: MIT
"""Retrieve one explicitly registered physical nominal-rating source."""

from pathlib import Path

import yaml

from scripts.audit_line_ratings import cached_source


if __name__ == "__main__":
    registry = yaml.safe_load(Path(snakemake.input.registry).read_text())["sources"]
    source_id = str(snakemake.wildcards.source_id)
    spec = dict(registry[source_id], id=source_id)
    if "path" in spec:
        raise ValueError("Local source files must be supplied directly, not retrieved as archives")
    target = Path(snakemake.output.archive)
    path, digest = cached_source(spec, target.parent, offline=False)
    if path != target:
        raise ValueError("Source cache path differs from declared Snakemake output")
    print(f"Physical nominal-rating source {source_id}: sha256={digest}")
