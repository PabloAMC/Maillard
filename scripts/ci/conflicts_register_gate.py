#!/usr/bin/env python
"""
conflicts_register_gate -- the register of conflicts is well-formed and points at real sources (2026-10-29).

Checks `docs/reference/conflicts_register.yml`: every entry has its fields; `kind`, `status` and `scope` are
from the vocabulary; every source names a dossier that exists; every ledger id in `affects` exists in a
measurement ledger, and every {file, symbol} pair names a file that contains the symbol; a resolved entry
(anything but `open`) says how it was resolved. Structural and offline: it does not judge an entry.

    python scripts/ci/conflicts_register_gate.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
REGISTER = ROOT / "docs" / "reference" / "conflicts_register.yml"
DOSSIERS = ROOT / "data" / "lit" / "extraction_dossiers"
LEDGERS = ROOT / "data" / "lit" / "measurement_ledger"
KINDS = {"conflict", "duplicate", "condition_error", "code_mismatch", "self_contradiction"}
STATUSES = {"open", "explained", "duplicate_confirmed", "error_confirmed", "fixed", "not_an_error"}
SCOPES = {"sulfur_rows", "other_lane", "corpus"}
FIELDS = ("id", "kind", "scope", "sources", "quantity", "observation", "explanations", "tests", "status",
          "remedy", "affects")


def _ledger_ids() -> set:
    ids = set()
    for path in sorted(LEDGERS.glob("*.yml")):
        for entry in (yaml.safe_load(path.read_text(encoding="utf-8")) or {}).get("entries") or []:
            ids.add(entry["id"])
    return ids


def check() -> list:
    problems = []
    entries = (yaml.safe_load(REGISTER.read_text(encoding="utf-8")) or {}).get("entries") or []
    ledger_ids = _ledger_ids()
    dossier_names = {p.name for p in DOSSIERS.iterdir()}
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
        for field, vocabulary in (("kind", KINDS), ("status", STATUSES), ("scope", SCOPES)):
            if entry[field] not in vocabulary:
                problems.append(f"{key}: {field} {entry[field]!r} not in {sorted(vocabulary)}")
        if not entry["sources"]:
            problems.append(f"{key}: no sources")
        for stem in entry["sources"] or []:
            # exact names, not Path.exists(): macOS disks are case-insensitive, CI's Linux is not (2026-10-29,
            # 'xin2026b' passed locally and failed in CI against Xin2026b_extraction.md)
            if f"{stem}_extraction.md" not in dossier_names and f"{stem}.md" not in dossier_names:
                problems.append(f"{key}: source {stem!r} has no dossier (names are case-sensitive)")
        for item in entry["affects"] or []:
            if isinstance(item, str):
                if item not in ledger_ids:
                    problems.append(f"{key}: ledger id {item!r} not in any measurement ledger")
            elif isinstance(item, dict) and {"file", "symbol"} <= set(item):
                path = ROOT / item["file"]
                if not path.exists():
                    problems.append(f"{key}: {item['file']} does not exist")
                elif str(item["symbol"]) not in path.read_text(encoding="utf-8"):
                    problems.append(f"{key}: symbol {item['symbol']!r} not found in {item['file']}")
            else:
                problems.append(f"{key}: affects item {item!r} is neither a ledger id nor {{file, symbol}}")
        if entry["status"] != "open" and not str(entry["remedy"] or "").strip() and not entry["tests"]:
            problems.append(f"{key}: status {entry['status']!r} but neither tests nor remedy say how")
    return problems


def main() -> int:
    problems = check()
    if problems:
        print("conflicts_register_gate: FAIL")
        for p in problems:
            print("  - " + p)
        return 1
    n = len((yaml.safe_load(REGISTER.read_text(encoding="utf-8")) or {}).get("entries") or [])
    print(f"conflicts_register_gate: PASS ({n} entries, docs/reference/conflicts_register.yml)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
