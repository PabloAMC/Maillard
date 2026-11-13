#!/usr/bin/env python
"""Compare the two blind readers' PREREG A6 field checks (docs/reference/lever_table/field_checks/, STEP8_FIELDS_SPEC
sec. 2): F1 exact, F2 within 0.1 pH unit, F3 on compound classes (not names). Prints a markdown report of agreement
and the rows the resolver must re-read; exit status 0 always (a report, not a gate).

    python scripts/tools/compare_field_checks.py [<stem> ...]      # stems like xu2007 or zhou2023.F3; default: all
"""
from __future__ import annotations

import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.generators.build_lever_field_readings import CHECKS_DIR, classes  # noqa: E402

TOL_PH = 0.1


def _load(stem: str, reader: str):
    f = CHECKS_DIR / f"{stem}.reader_{reader}.yml"
    return yaml.safe_load(f.read_text()) if f.exists() else None


def _cls(names):
    out = set()
    for n in names or []:
        out |= classes(n)
    return sorted(out)


def compare(stem: str):
    a, b = _load(stem, "a"), _load(stem, "b")
    rows, tally = [], {"F1": [0, 0], "F2": [0, 0], "F3": [0, 0]}
    if a is None or b is None:
        return [f"| {stem} | - | - | a reader file is missing |"], tally
    fields = set(a.get("fields") or []) & set(b.get("fields") or [])
    if "F1" in fields:
        ok = a.get("ph_control") == b.get("ph_control") and (a.get("ph_control_overrides") or {}) == (
            b.get("ph_control_overrides") or {})
        tally["F1"][0] += ok
        tally["F1"][1] += 1
        if not ok:
            rows.append(f"| {stem} | F1 | (source) | a {a.get('ph_control')} / b {b.get('ph_control')} |")
    A, B = a.get("arms") or {}, b.get("arms") or {}
    for label in sorted(set(A) | set(B)):
        ra, rb = A.get(label) or {}, B.get(label) or {}
        if "F2" in fields:
            fa, fb = ra.get("final_ph"), rb.get("final_ph")
            ok = (fa is None and fb is None) or (fa is not None and fb is not None and abs(float(fa) - float(fb)) <= TOL_PH)
            tally["F2"][0] += ok
            tally["F2"][1] += 1
            if not ok:
                rows.append(f"| {stem} | F2 | {label} | a {fa} / b {fb} |")
        if "F3" in fields:
            ca, cb = _cls(ra.get("present")), _cls(rb.get("present"))
            ok = ca == cb and label in A and label in B
            tally["F3"][0] += ok
            tally["F3"][1] += 1
            if not ok:
                rows.append(f"| {stem} | F3 | {label} | a {'+'.join(ca)} {ra.get('present')} / "
                            f"b {'+'.join(cb)} {rb.get('present')} |")
    return rows, tally


def main(argv=None) -> int:
    stems = list(argv if argv is not None else sys.argv[1:]) or sorted(
        {p.name.split(".reader_")[0] for p in CHECKS_DIR.glob("*.reader_*.yml")})
    total = {"F1": [0, 0], "F2": [0, 0], "F3": [0, 0]}
    out = ["| source | field | arm | readers |", "|---|---|---|---|"]
    for s in stems:
        rows, t = compare(s)
        out += rows
        for k in total:
            total[k][0] += t[k][0]
            total[k][1] += t[k][1]
    print("Agreement: " + "; ".join(f"{k} {v[0]}/{v[1]}" for k, v in total.items() if v[1]))
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
