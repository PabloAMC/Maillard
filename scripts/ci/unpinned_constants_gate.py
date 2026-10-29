#!/usr/bin/env python
"""
unpinned_constants_gate -- no unpinned constant that matters is silently left out of the uncertainty (2026-10-29).

Reads the generated sulfur reaction evidence table (``results/validation/sulfur_reaction_evidence.json``): every fitted,
unidentified constant on a step that matters whose RATE the envelope does not sample must be listed, with the sub-model
that will test it, in ``docs/reference/unpinned_constants.yml``; and every listed constant must still be one (the list
only shrinks). Structural and offline (JSON and YAML only).

    python scripts/ci/unpinned_constants_gate.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
TABLE = ROOT / "results" / "validation" / "sulfur_reaction_evidence.json"
REGISTER = ROOT / "docs" / "reference" / "unpinned_constants.yml"


def check() -> list:
    summary = json.loads(TABLE.read_text(encoding="utf-8"))["summary"]
    must = set(summary["unpinned_and_rate_not_sampled"])
    unpinned = set(summary["unpinned_constants_that_matter"])
    listed = {}
    problems = []
    for entry in (yaml.safe_load(REGISTER.read_text(encoding="utf-8")) or {}).get("entries") or []:
        key = entry.get("key")
        if key in listed:
            problems.append(f"{key}: listed twice")
        listed[key] = entry
        if not str(entry.get("planned_test") or "").strip():
            problems.append(f"{key}: no planned_test")
    for key in sorted(must - set(listed)):
        problems.append(f"{key}: unpinned, matters, its rate is not sampled by the envelope, and it is not in {REGISTER.name}")
    for key in sorted(set(listed) - unpinned):
        problems.append(f"{key}: listed but no longer an unpinned constant that matters; remove it (the list only shrinks)")
    return problems


def main() -> int:
    problems = check()
    if problems:
        print("unpinned_constants_gate: FAIL")
        for p in problems:
            print("  - " + p)
        return 1
    n = len((yaml.safe_load(REGISTER.read_text(encoding="utf-8")) or {}).get("entries") or [])
    print(f"unpinned_constants_gate: PASS ({n} unpinned constants that matter, each with its planned test)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
