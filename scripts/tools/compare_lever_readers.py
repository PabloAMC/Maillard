#!/usr/bin/env python
"""Compare the two blind readers' keyings of one lever-table source (docs/reference/lever_table/PREREG.md sec. 3).

Only the fields that define the cell, the lever step and the band are compared: odorant, lever class, the lever
values, unit type and the arm values (a difference above 0.05 dex). Contrasts are matched by (resolved odorant,
lever class, location, sorted lever values). Prints a markdown report of the rows the resolver must re-read; exit
status 0 always (it is a report, not a gate).

    python scripts/tools/compare_lever_readers.py <source> [<source> ...]
"""
from __future__ import annotations

import math
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple

import yaml

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import compound_keys, data_paths  # noqa: E402

TOL_DEX = 0.05


def _num(x) -> Optional[float]:
    try:
        return float(str(x).strip().split()[0].replace(",", ""))
    except (ValueError, IndexError, TypeError):
        return None


def _odorant(s: str) -> str:
    k = compound_keys.resolve(str(s))
    return k.display_name if k else f"UNRESOLVED:{s}"


def _loc(s: str) -> str:
    m = re.search(r"(table|fig\w*)\s*\.?\s*([0-9]+[a-z]?)", str(s).lower())
    return f"{m.group(1)[:3]} {m.group(2)}" if m else str(s).strip().lower()


def _lv(a: Mapping[str, Any]) -> str:
    v = _num(a.get("lever_value"))
    return f"{v:g}" if v is not None else str(a.get("lever_value")).strip().lower()


def _key(c: Mapping[str, Any]) -> Tuple:
    return (_odorant(c.get("odorant")), c.get("lever_class"), _loc(c.get("location", "")),
            tuple(sorted(_lv(a) for a in c.get("arms") or [])))


def _arm_val(a):
    cen = a.get("censor") or "none"
    v = _num(a.get("limit")) if cen in ("nd_with_limit", "below_value") else _num(a.get("value"))
    return cen, v


def compare(source: str) -> List[str]:
    d = data_paths.LEVER_TABLE_CONTRASTS_DIR
    ra = yaml.safe_load((d / f"{source}.reader_a.yml").read_text()) or {}
    rb = yaml.safe_load((d / f"{source}.reader_b.yml").read_text()) or {}
    from collections import defaultdict

    GA, GB = defaultdict(list), defaultdict(list)
    for c in ra.get("contrasts") or []:
        GA[_key(c)].append(c)
    for c in rb.get("contrasts") or []:
        GB[_key(c)].append(c)

    def dist(ca, cb):
        vb = {_lv(a): _arm_val(a)[1] for a in cb.get("arms") or []}
        d = 0.0
        for a in ca.get("arms") or []:
            x, y = _arm_val(a)[1], vb.get(_lv(a))
            d += abs(math.log10(x) - math.log10(y)) if x and y and x > 0 and y > 0 else (0.0 if x == y else 1.0)
        return d

    matched, only_a, only_b = [], [], []
    for k in set(GA) | set(GB):
        la, lb = list(GA.get(k, [])), list(GB.get(k, []))
        while la and lb:
            i, j = min(((i, j) for i in range(len(la)) for j in range(len(lb))), key=lambda t: dist(la[t[0]], lb[t[1]]))
            matched.append((k, la.pop(i), lb.pop(j)))
        only_a += [(k, c) for c in la]
        only_b += [(k, c) for c in lb]
    n_a, n_b = sum(map(len, GA.values())), sum(map(len, GB.values()))
    out = [f"## {source}", "", f"reader a: {n_a} contrasts; reader b: {n_b}; matched {len(matched)}", ""]
    for k, c in sorted(only_a, key=lambda t: str(t[0])):
        out.append(f"- ONLY A: {c.get('id')} {k[0]} / {k[1]} / {k[2]} / levels {list(k[3])}")
    for k, c in sorted(only_b, key=lambda t: str(t[0])):
        out.append(f"- ONLY B: {c.get('id')} {k[0]} / {k[1]} / {k[2]} / levels {list(k[3])}")
    for k, ca, cb in sorted(matched, key=lambda t: str(t[0])):
        diffs = []
        if (ca.get("quant_type") or "") != (cb.get("quant_type") or ""):
            diffs.append(f"quant_type {ca.get('quant_type')!r} vs {cb.get('quant_type')!r}")
        arms_b = {_lv(a): a for a in cb.get("arms") or []}
        for aa in ca.get("arms") or []:
            ab = arms_b.get(_lv(aa))
            if ab is None:
                continue
            (cea, va), (ceb, vb) = _arm_val(aa), _arm_val(ab)
            if cea != ceb:
                diffs.append(f"arm {_lv(aa)}: censor {cea} vs {ceb}")
            elif va is not None and vb is not None and va > 0 and vb > 0:
                if abs(math.log10(va) - math.log10(vb)) > TOL_DEX:
                    diffs.append(f"arm {_lv(aa)}: {va:g} vs {vb:g}")
            elif (va is None) != (vb is None):
                diffs.append(f"arm {_lv(aa)}: {va} vs {vb}")
        if diffs:
            out.append(f"- DIFFER {ca.get('id')} / {cb.get('id')} ({k[0]}, {k[1]}, {k[2]}): " + "; ".join(diffs))
    out.append("")
    return out


def main(argv: List[str]) -> int:
    for s in argv:
        print("\n".join(compare(s)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
