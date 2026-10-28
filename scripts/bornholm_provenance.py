# SPDX-License-Identifier: MIT
"""Record Bornholm source provenance without requiring Git on compute nodes."""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


def code_provenance(root=Path(".")):
    """Record code hashes and, when available, the repository commit."""
    root = Path(root)
    paths = [root / name for name in ["Snakefile", "pixi.toml", "pixi.lock"]]
    for directory in ["scripts", "rules", "config", "profiles"]:
        paths.extend(
            p
            for p in (root / directory).rglob("*")
            if p.suffix in [".py", ".smk", ".yaml", ".yml", ".json", ".sbatch", ".sh"]
        )
    hashes = {
        p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(paths)
        if p.is_file()
    }
    report = {
        "code_files_sha256": hashes,
        "code_sha256": hashlib.sha256(
            json.dumps(hashes, sort_keys=True).encode()
        ).hexdigest(),
        "git_commit": None,
    }
    try:
        report["git_commit"] = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            text=True,
            stderr=subprocess.PIPE,
            timeout=10,
        ).strip()
    except (OSError, subprocess.SubprocessError) as error:
        report["git_provenance_error"] = f"{type(error).__name__}: {error}"
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report = code_provenance()
    (args.output_dir / "git_head.txt").write_text(
        (report["git_commit"] or "unavailable") + "\n"
    )
    if report["git_commit"] is not None:
        try:
            patch = subprocess.check_output(
                ["git", "diff", "HEAD"], text=True, stderr=subprocess.PIPE, timeout=10
            )
            (args.output_dir / "code_changes.patch").write_text(patch)
        except (OSError, subprocess.SubprocessError) as error:
            report["git_diff_error"] = f"{type(error).__name__}: {error}"
    (args.output_dir / "code_provenance.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    if report["git_commit"] is None:
        print(
            "Git revision unavailable; recorded source-file hashes in code_provenance.json."
        )


if __name__ == "__main__":
    main()
