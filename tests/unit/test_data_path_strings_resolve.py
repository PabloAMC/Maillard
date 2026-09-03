"""Every repo path a curated data file names must exist.

Found 2026-09-03: six intake references silently vanished from the literature backlog
because their ``runtime_artifacts[].path`` pointed at ``data/lit/slr_incorporation_matrix.json``,
which had moved to ``results/literature/`` the day before; two more had pointed at a
benchmark quarantined weeks earlier. The registry stores file locations as data and
nothing checked them. This test does, for every string under ``data/`` that looks like a
repository path.
"""
from __future__ import annotations

import glob
import re
from pathlib import Path

from src import data_access, data_paths

PATH_RE = re.compile(r"^(data|results|docs|src|scripts|tests)/[A-Za-z0-9_./-]+\.[A-Za-z0-9]+$")
# Paths that are cited as history (a retired file named in a provenance note) rather than
# as a live reference. Keep this list short and dated.
HISTORICAL = {
    "data/benchmarks/maillard_validation_benchmarks.md",  # retired record, kept at its path
}


def _strings(node, out):
    if isinstance(node, dict):
        for value in node.values():
            _strings(value, out)
    elif isinstance(node, list):
        for value in node:
            _strings(value, out)
    elif isinstance(node, str) and PATH_RE.match(node.strip()):
        out.add(node.strip())


def test_every_repo_path_named_in_a_curated_data_file_exists():
    files = sorted(
        glob.glob(str(data_paths.DATA_ROOT / "**" / "*.json"), recursive=True)
        + glob.glob(str(data_paths.DATA_ROOT / "**" / "*.yml"), recursive=True)
        + glob.glob(str(data_paths.DATA_ROOT / "**" / "*.yaml"), recursive=True)
    )
    missing = {}
    for path in files:
        payload = data_access.load_mapping(path, missing_ok=True)
        if payload is None:
            continue
        found: set[str] = set()
        _strings(payload, found)
        for rel in sorted(found):
            if rel in HISTORICAL:
                continue
            if not (data_paths.REPO_ROOT / rel).exists():
                missing.setdefault(data_paths.rel(Path(path)), []).append(rel)
    assert not missing, (
        "curated data files name repository paths that do not exist "
        "(update the reference or, for a retired file, add it to HISTORICAL with a date):\n"
        + "\n".join(f"  {f}: {p}" for f, ps in missing.items() for p in ps)
    )
