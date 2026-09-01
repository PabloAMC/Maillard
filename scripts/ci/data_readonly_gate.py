#!/usr/bin/env python
"""Data read-only gate.

``data/`` is curated input. Nothing that runs -- a test tier, a generator, a CLI
smoke run -- may leave a modified, added or deleted file under it. Until
2026-09-01 the integration suite wrote a calibration-history JSON into ``data/``
on every run (102 of them had accumulated), a generator could silently rewrite the
frozen hold-out bundles, and the ranker committed its own output as a protocol.
None of that was visible to CI because nothing looked.

This gate is run AFTER a test tier and fails if ``git status`` reports anything
under ``data/``. It is stdlib-only and offline like the other gates in this tier.

Usage:
    python scripts/ci/data_readonly_gate.py            # check data/
    python scripts/ci/data_readonly_gate.py data docs  # check several trees
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TREES = ("data",)


def _git_lines(*args: str) -> list[str]:
    completed = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=True
    )
    return [line for line in completed.stdout.splitlines() if line.strip()]


def dirty_paths(trees: tuple[str, ...]) -> list[str]:
    """Paths under ``trees`` changed since the index, plus untracked files.

    Deliberately compares the WORKING TREE against the INDEX (not HEAD): a developer
    with staged edits under data/ is not what this gate is for. What it catches is a
    file that something wrote, rewrote or deleted while running.
    """
    modified = [f"M {p}" for p in _git_lines("diff", "--name-only", "--", *trees)]
    untracked = [
        f"? {p}"
        for p in _git_lines("ls-files", "--others", "--exclude-standard", "--", *trees)
    ]
    return modified + untracked


def main(argv: list[str] | None = None) -> int:
    trees = tuple(argv) if argv else DEFAULT_TREES
    dirty = dirty_paths(trees)
    if not dirty:
        print(f"data_readonly_gate: PASS ({', '.join(trees)} untouched)")
        return 0
    print(
        f"data_readonly_gate: FAIL -- {len(dirty)} path(s) under {', '.join(trees)} were "
        "modified, added or deleted by something that ran. Curated inputs are read-only at "
        "runtime; generated artifacts belong under results/.",
        file=sys.stderr,
    )
    for line in dirty:
        print(f"  {line}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
