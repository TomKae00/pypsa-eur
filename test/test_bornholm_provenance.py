# SPDX-License-Identifier: MIT
"""Revision bookkeeping must work without compute-node Git."""

import subprocess

from scripts.bornholm_provenance import code_provenance


def test_source_hashes_survive_missing_git_and_track_edits(tmp_path, monkeypatch):
    (tmp_path / "scripts").mkdir()
    source = tmp_path / "scripts/model.py"
    source.write_text("x = 1\n")

    def missing_git(*args, **kwargs):
        raise FileNotFoundError("git")

    monkeypatch.setattr(subprocess, "check_output", missing_git)
    before = code_provenance(tmp_path)
    assert before["git_commit"] is None
    assert "FileNotFoundError" in before["git_provenance_error"]
    assert set(before["code_files_sha256"]) == {"scripts/model.py"}
    source.write_text("x = 2\n")
    assert code_provenance(tmp_path)["code_sha256"] != before["code_sha256"]


def test_git_failure_is_optional_provenance(tmp_path, monkeypatch):
    def failed_git(*args, **kwargs):
        raise subprocess.CalledProcessError(128, ["git", "rev-parse", "HEAD"])

    monkeypatch.setattr(subprocess, "check_output", failed_git)
    report = code_provenance(tmp_path)
    assert report["git_commit"] is None
    assert "CalledProcessError" in report["git_provenance_error"]
