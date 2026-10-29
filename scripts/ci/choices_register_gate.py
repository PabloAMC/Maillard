"""
choices_register_gate -- the register of non-measured choices is complete and points at real code (2026-10-28).

Checks `docs/reference/choices_register.yml`: every entry has its fields; `kind` and `sensitivity.status` are
from the vocabulary; an `alternatives` or `cutoff` entry carries at least two alternatives; a `fact` or
`policy` entry carries none and is `not_applicable`; and the `where.symbol` string occurs in `where.file`, so
an entry cannot outlive the code it describes. Structural and offline: it does not judge a choice.

    python scripts/ci/choices_register_gate.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
REGISTER = ROOT / "docs" / "reference" / "choices_register.yml"
KINDS = {"fact", "alternatives", "cutoff", "policy"}
STATUSES = {"untested", "harmless", "moves_headline", "not_applicable"}
FIELDS = ("id", "what", "where", "kind", "current", "basis", "alternatives", "decided_by", "sensitivity")


def check() -> list:
    problems = []
    entries = (yaml.safe_load(REGISTER.read_text(encoding="utf-8")) or {}).get("entries") or []
    seen = set()
    for entry in entries:
        key = entry.get("id", "?")
        missing = [f for f in FIELDS if f not in entry]
        if missing:
            problems.append(f"{key}: missing {missing}")
            continue
        if key in seen:
            problems.append(f"{key}: duplicate id")
        seen.add(key)
        kind, status = entry["kind"], (entry["sensitivity"] or {}).get("status")
        if kind not in KINDS:
            problems.append(f"{key}: kind {kind!r} not in {sorted(KINDS)}")
        if status not in STATUSES:
            problems.append(f"{key}: sensitivity.status {status!r} not in {sorted(STATUSES)}")
        n_alt = len(entry["alternatives"] or [])
        if kind in ("alternatives", "cutoff") and n_alt < 2:
            problems.append(f"{key}: a {kind} entry needs at least two alternatives")
        if kind in ("fact", "policy") and (n_alt or status != "not_applicable"):
            problems.append(f"{key}: a {kind} entry carries no alternatives and is not_applicable")
        where = entry["where"] or {}
        path = ROOT / str(where.get("file", ""))
        if not path.is_file():
            problems.append(f"{key}: {where.get('file')!r} does not exist")
        elif str(where.get("symbol", "")) not in path.read_text(encoding="utf-8"):
            problems.append(f"{key}: {where.get('symbol')!r} not found in {where.get('file')}")
    if not entries:
        problems.append("the register is empty")
    return problems


def main() -> int:
    problems = check()
    if problems:
        print("choices_register_gate: FAIL\n  " + "\n  ".join(problems))
        return 1
    print(f"choices_register_gate: PASS ({REGISTER.relative_to(ROOT)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
