#!/usr/bin/env python
"""
source_verification_gate -- nothing enters a fit unverified against its paper (2026-10-29).

Checks `data/lit/source_verification.yml`: every entry is well-formed, names a ledger id that
exists, and records all five checks; and every row of a keyed fit table (`data/lit/fit_rows/*.yml`) points at a
ledger id verified there with every check passed. It also REPORTS, without failing, how many ledger entries that
already feed the model (`in_model_fit`) are unverified: the debt the older waves carry.

    python scripts/ci/source_verification_gate.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
LEDGERS = ROOT / "data" / "lit" / "measurement_ledger"
VERIFICATION = ROOT / "data" / "lit" / "source_verification.yml"
VERIFICATION_DIR = ROOT / "data" / "lit" / "source_verification"   # one file per sub-model (2026-10-29)
FIT_ROWS = ROOT / "data" / "lit" / "fit_rows"
CHECKS = ("conditions", "units", "quantification", "duplicates", "values")
FIELDS = ("ledger_id", "source_file", "read_by", "second_reader", "date", "checks", "quote")
VERDICTS = {"pass", "corrected", "fail"}


def _ledger_entries() -> dict:
    out = {}
    for path in sorted(LEDGERS.glob("*.yml")):
        for entry in (yaml.safe_load(path.read_text(encoding="utf-8")) or {}).get("entries") or []:
            out[entry["id"]] = entry
    return out


def check() -> tuple:
    problems, notes = [], []
    ledger = _ledger_entries()
    verified = {}
    files = [VERIFICATION, *sorted(VERIFICATION_DIR.glob("*.yml"))] if VERIFICATION_DIR.exists() else [VERIFICATION]
    for entry in (e for f in files for e in (yaml.safe_load(f.read_text(encoding="utf-8")) or {}).get("entries") or []):
        key = entry.get("ledger_id", "?")
        missing = [f for f in FIELDS if f not in entry]
        if missing:
            problems.append(f"{key}: missing {missing}")
            continue
        if key not in ledger:
            problems.append(f"{key}: not in any measurement ledger")
        checks = entry["checks"] or {}
        absent = [c for c in CHECKS if c not in checks]
        if absent:
            problems.append(f"{key}: checks missing {absent}")
        bad = {c: v for c, v in checks.items() if v not in VERDICTS}
        if bad:
            problems.append(f"{key}: verdicts not in {sorted(VERDICTS)}: {bad}")
        if not str(entry["quote"] or "").strip():
            problems.append(f"{key}: no quote from the source")
        if str(entry["second_reader"] or "").strip() == str(entry["read_by"] or "").strip():
            problems.append(f"{key}: the second reader must differ from the first")
        # a check by one reader is not a verification: "pending" (the second reader has not run) is unverified
        second = str(entry["second_reader"] or "").strip()
        verified[key] = all(checks.get(c) in ("pass", "corrected") for c in CHECKS) and second not in ("", "pending")
    rows_checked, drafts = 0, []
    for path in sorted(FIT_ROWS.glob("*.yml")) if FIT_ROWS.exists() else []:
        table = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        if table.get("status", "draft") != "ready":
            # a table still being read and verified: reported, not failed; fit_rows.load refuses it for fitting
            drafts.append(path.name)
            continue
        cited = {pot.get("ledger_id") for pot in (table.get("pots") or {}).values()}
        rows_checked += len(table.get("rows") or [])
        for key in sorted(k for k in cited if k):
            if not verified.get(key, False):
                problems.append(f"{path.name} is ready but its source {key!r} is not verified by two readers")
    debt = sorted(k for k, e in ledger.items() if e.get("status") == "in_model_fit" and not verified.get(k))
    notes.append(f"draft tables (not yet fit-ready): {drafts or 'none'}; {sum(verified.values())} ledger entries verified by two readers ({len(verified)} read); {rows_checked} keyed fit rows checked; "
                 f"{len(debt)} in_model_fit entries of earlier waves unverified (debt, reported, not failed)")
    return problems, notes


def main() -> int:
    problems, notes = check()
    if problems:
        print("source_verification_gate: FAIL")
        for p in problems:
            print("  - " + p)
        return 1
    print("source_verification_gate: PASS (" + "; ".join(notes) + ")")
    return 0


if __name__ == "__main__":
    sys.exit(main())
