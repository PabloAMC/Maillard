#!/usr/bin/env python
"""
reconcile_source_readings.py -- two independent readings of the same sources become one verified record (2026-10-29).

A source enters a fit only when two readers verified it against the paper (data/lit/source_verification*). The readers
work in parallel and without seeing each other (double entry): the first writes data/lit/source_verification/<name>.yml
(sub-model 1's first reader: data/lit/source_verification.yml) and data/lit/fit_rows/<name>.yml; the second writes the
same names under second/. This script compares them, per ledger id:
  * inventory   -- a data set only one reader chose is listed, not promoted;
  * verdicts    -- both readers' five checks must be pass / corrected;
  * pots        -- matched across readers by their numeric conditions (the two name pots differently): temperature,
                   time, initial mmol/L, label pH must agree to 1 %;
  * rows        -- matched by (pot, species, time); values must agree within the larger of the two rows' errors
                   (dex for amounts, pH units for pH), or to 1 % where both transcribe a printed number exactly.
It writes results/validation/source_reading_reconciliation_<name>.md. With --promote it records the second reader in
the first reader's verification entry for every ledger id that agrees on all four, which makes that source verified;
every disagreement stays for a person to settle against the paper.

    python scripts/reconcile_source_readings.py ph_acid_submodel [--promote]
"""
from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402
from src.kinetic_core import fit_rows  # noqa: E402

CHECKS = ("conditions", "units", "quantification", "duplicates", "values")
REL = 0.01


def _first_verification(name: str) -> Path:
    own = data_paths.SOURCE_VERIFICATION_DIR / f"{name}.yml"
    return own if own.exists() else data_paths.SOURCE_VERIFICATION


def _entries(path: Path) -> Dict[str, dict]:
    if not path.exists():
        return {}
    return {e["ledger_id"]: e for e in (yaml.safe_load(path.read_text(encoding="utf-8")) or {}).get("entries") or []}


def _table(path: Path) -> Tuple[dict, list]:
    if not path.exists():
        return {}, []
    t = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    return t.get("pots") or {}, t.get("rows") or []


def _close(a, b) -> bool:
    if a is None or b is None:
        return a is None and b is None
    a, b = float(a), float(b)
    return math.isclose(a, b, rel_tol=REL, abs_tol=1e-9)


def _pot_signature_match(p: dict, q: dict) -> bool:
    if p.get("ledger_id") != q.get("ledger_id"):
        return False
    if not (_close(p.get("t_c"), q.get("t_c")) and _close(p.get("minutes"), q.get("minutes"))
            and _close(p.get("label_ph"), q.get("label_ph"))):
        return False
    a, b = p.get("initial_mmol_L") or {}, q.get("initial_mmol_L") or {}
    keys = {k for k, v in a.items() if v} | {k for k, v in b.items() if v}
    return all(_close(a.get(k, 0.0), b.get(k, 0.0)) for k in keys)


def _row_error(row: dict) -> float:
    try:
        return fit_rows._sigma(row)
    except Exception:  # noqa: BLE001 - a row without an error is reported by the loader, not here
        return 0.0


def _rows_agree(r: dict, s: dict) -> bool:
    if r.get("unit") != s.get("unit"):
        return False
    v, w = float(r["value"]), float(s["value"])
    if math.isclose(v, w, rel_tol=REL, abs_tol=1e-12):
        return True
    tol = max(_row_error(r), _row_error(s))
    diff = abs(v - w) if r["unit"] == "pH" else abs(math.log10(max(v, 1e-30) / max(w, 1e-30)))
    return diff <= tol


def reconcile(name: str) -> dict:
    first_v = _entries(_first_verification(name))
    second_v = _entries(data_paths.SOURCE_VERIFICATION_DIR / "second" / f"{name}.yml")
    pots1, rows1 = _table(fit_rows.FIT_ROWS_DIR / f"{name}.yml")
    pots2, rows2 = _table(fit_rows.FIT_ROWS_DIR / "second" / f"{name}.yml")
    out = {"only_first": sorted(set(first_v) - set(second_v)), "only_second": sorted(set(second_v) - set(first_v)),
           "agree": [], "disagree": {}}
    for lid in sorted(set(first_v) & set(second_v)):
        problems: List[str] = []
        for reader, entry in (("first", first_v[lid]), ("second", second_v[lid])):
            bad = [c for c in CHECKS if (entry.get("checks") or {}).get(c) not in ("pass", "corrected")]
            if bad:
                problems.append(f"{reader} reader: {', '.join(f'{c}={(entry.get('checks') or {}).get(c)}' for c in bad)}")
        mine = {k: p for k, p in pots1.items() if p.get("ledger_id") == lid}
        theirs = {k: p for k, p in pots2.items() if p.get("ledger_id") == lid}
        mapping = {}
        for k, p in mine.items():
            match = [k2 for k2, q in theirs.items() if _pot_signature_match(p, q) and k2 not in mapping.values()]
            if match:
                mapping[k] = match[0]
            else:
                problems.append(f"pot {k} ({p.get('t_c')} C, {p.get('minutes')} min, label pH {p.get('label_ph')}, "
                                f"{p.get('initial_mmol_L')}) has no matching pot in the second reading")
        for k2 in set(theirs) - set(mapping.values()):
            q = theirs[k2]
            problems.append(f"second reader's pot {k2} ({q.get('t_c')} C, {q.get('minutes')} min, {q.get('initial_mmol_L')}) "
                            f"has no matching pot in the first reading")
        index2 = {(r["pot"], str(r["species"]), float(r["t_min"])): r for r in rows2}
        for r in rows1:
            if r.get("pot") not in mapping:
                continue
            key = (mapping[r["pot"]], str(r["species"]), float(r["t_min"]))
            s = index2.pop(key, None)
            if s is None:
                problems.append(f"row {r['pot']}/{r['species']}@{r['t_min']} = {r['value']} {r['unit']}: not in the second reading")
            elif not _rows_agree(r, s):
                problems.append(f"row {r['pot']}/{r['species']}@{r['t_min']}: {r['value']} vs {s['value']} {r['unit']}")
        for (pot, species, t), s in index2.items():
            if pot in mapping.values():
                problems.append(f"second reader's row {pot}/{species}@{t} = {s['value']} {s['unit']}: not in the first reading")
        if problems:
            out["disagree"][lid] = problems
        else:
            out["agree"].append(lid)
    return out


def promote(name: str, agreed: List[str], second_v: Dict[str, dict]) -> int:
    path = _first_verification(name)
    text = path.read_text(encoding="utf-8")
    n = 0
    for lid in agreed:
        reader = str(second_v[lid].get("read_by") or "second reader")
        pattern = re.compile(rf"(- ledger_id: {re.escape(lid)}\n(?:(?!  - ledger_id:).*\n)*?\s+second_reader: )pending")
        text, k = pattern.subn(rf"\g<1>{reader}", text, count=1)
        n += k
    path.write_text(text, encoding="utf-8")
    return n


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("name")
    parser.add_argument("--promote", action="store_true")
    args = parser.parse_args(argv)
    result = reconcile(args.name)
    lines = [f"# Source readings reconciled: {args.name}", "",
             f"*Generated by `scripts/reconcile_source_readings.py {args.name}`. Two readers, independent and in parallel; "
             "a ledger id is promoted to verified only when inventory, verdicts, pots and rows all agree.*", "",
             f"**Agree ({len(result['agree'])}):** " + (", ".join(f"`{x}`" for x in result["agree"]) or "none"), "",
             f"**Chosen by the first reader only ({len(result['only_first'])}):** " + (", ".join(f"`{x}`" for x in result["only_first"]) or "none"), "",
             f"**Chosen by the second reader only ({len(result['only_second'])}):** " + (", ".join(f"`{x}`" for x in result["only_second"]) or "none"), "",
             f"## Disagreements ({len(result['disagree'])})", ""]
    for lid, problems in result["disagree"].items():
        lines.append(f"- `{lid}`")
        lines += [f"  - {p}" for p in problems[:12]] + ([f"  - ... {len(problems) - 12} more"] if len(problems) > 12 else [])
    report = data_paths.VALIDATION_DIR / f"source_reading_reconciliation_{args.name}.md"
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"{args.name}: {len(result['agree'])} agree, {len(result['disagree'])} disagree, "
          f"{len(result['only_first'])} first-only, {len(result['only_second'])} second-only -> {data_paths.rel(report)}")
    if args.promote and result["agree"]:
        second_v = _entries(data_paths.SOURCE_VERIFICATION_DIR / "second" / f"{args.name}.yml")
        print(f"promoted {promote(args.name, result['agree'], second_v)} ledger ids to verified by two readers")
    return 0


if __name__ == "__main__":
    sys.exit(main())
