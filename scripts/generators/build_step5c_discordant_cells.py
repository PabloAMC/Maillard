#!/usr/bin/env python
"""
The engine on the lever table's discordant cells (STATUS_AND_PLAN sec. 0, step 5c).

Question: for each cell the disagreement audit classes "conflict between laboratories", does the engine as shipped,
run on each laboratory's own keyed conditions, reproduce the difference between the laboratories?

Rules (frozen in this file before ``--score`` was first run; a rule changed afterwards is recorded in the plan):

* Cells: the audit's (``lever_disagreement_audit.json``) cells of class "conflict between laboratories".
* Items: every keyed contrast of those cells, both partitions (DEVELOPMENT and TEST), built exactly as step 5a builds
  its items (``build_step5a_display_check.cell_items``), and mapped to engine runs by step 5a's ``classify`` under
  its three readings (front door; added precursors; added precursors with the pH held). The front door is primary.
* Per laboratory: the mean measured log ratio per step over its items, and the engine's mean over the items it ran.
* A cell is **explained by conditions** when the engine runs at least one item of each of its two laboratories and
  its between-laboratory difference has the measured difference's sign AND misses it by less than the measured
  difference itself (it beats "no difference between the laboratories"). It is **not explained** when the engine
  runs both and fails that, and **not runnable** otherwise (the reasons are counted, as in step 5a).
* This is a diagnostic: no display, no credit, no constant moves. It is the second engine look at TEST contrasts
  (sec. 0 rules allow three before 01-23).

    ./scripts/docker_maillard.sh run "python scripts/generators/build_step5c_discordant_cells.py --items"
    ./scripts/docker_maillard.sh run "python scripts/generators/build_step5c_discordant_cells.py --score"
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT, ROOT / "scripts" / "generators"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import build_lever_table as LT  # noqa: E402
import build_step5a_display_check as S5A  # noqa: E402
from src import artifact_io, data_paths, provenance  # noqa: E402

AUDIT = data_paths.VALIDATION_DIR / "lever_disagreement_audit.json"
ITEMS_OUT = data_paths.VALIDATION_DIR / "step5c_discordant_cells_items.json"
SCORE_OUT = data_paths.VALIDATION_DIR / "step5c_discordant_cells.json"
CONFLICT = "conflict between laboratories"


def conflict_cells() -> List[Dict[str, Any]]:
    return [c for c in json.loads(AUDIT.read_text())["cells"] if c["class"] == CONFLICT]


def build() -> Dict[str, Any]:
    """The conflict cells, their items (both partitions) and each item's mapping under every reading."""
    labs, _ = LT.lab_partition()
    obs = S5A.keyed_observations(LT.load_sources(), labs)
    by_cell = defaultdict(list)
    for o in obs:
        by_cell[(o["compound"], o["lever_class"], o["band"])].append(o)
    cells = []
    for c in conflict_cells():
        key = (c["compound"], c["lever_class"], c["band"])
        if key in by_cell:          # a numeric cell: the audit's band is the table's
            items = S5A.cell_items(key, by_cell[key])
        else:                       # a categorical pair: the audit names it 'b vs a', levels in name order
            items = []
            for k, rows in by_cell.items():
                if k[:2] == key[:2]:
                    items += [it for it in (_oriented(i) for i in S5A.cell_items(k, rows)) if it["pair"] == key[2]]
        cells.append({"cell": list(key), "lab_means_audit": c["lab_means"], "items": items})
    return {"cells": cells}


def _oriented(it: Dict[str, Any]) -> Dict[str, Any]:
    """A categorical item turned to the audit's orientation (log of the later level in name order over the earlier,
    ``build_lever_disagreement_audit.items``); a numeric item unchanged."""
    if it["kind"] != "categorical":
        return {**it, "pair": it["band"]}
    a, b = sorted((it["level"], it["against"]))
    if (it["level"], it["against"]) == (b, a):
        return {**it, "pair": f"{b} vs {a}"}
    return {**it, "pair": f"{b} vs {a}", "m": -it["m"], "level": b, "against": a, "arms": list(reversed(it["arms"]))}


def _ident(it: Mapping[str, Any]) -> Dict[str, Any]:
    return {k: it.get(k) for k in ("lab", "source", "contrast", "level", "against", "kind", "m")}


def items_payload() -> Dict[str, Any]:
    data = build()
    out_cells = []
    for c in data["cells"]:
        rows = []
        for it in c["items"]:
            rows.append({**_ident(it), "mapping": {r: S5A.classify(it, r)[0] for r in S5A.READINGS}})
        out_cells.append({"cell": c["cell"], "laboratories": sorted({it["lab"] for it in c["items"]}),
                          "items": rows})
    counts = {r: dict(Counter(S5A._why_class(row["mapping"][r]) for c in out_cells for row in c["items"])
                      .most_common()) for r in S5A.READINGS}
    return {"artifact": ITEMS_OUT.stem, "status": "diagnostic (step 5c); no engine output",
            "cells": out_cells, "mapping_counts": counts}


def _mean(xs: Sequence[float]) -> Optional[float]:
    return sum(xs) / len(xs) if xs else None


def judge(lab_rows: Mapping[str, List[Mapping[str, Any]]]) -> Dict[str, Any]:
    """The cell's verdict from each laboratory's measured and engine means (rules in the module docstring)."""
    labs = sorted(lab_rows)
    per_lab = {lab: {"measured": _mean([r["m"] for r in rows]),
                     "engine": _mean([r["engine"] for r in rows if r.get("engine") is not None]),
                     "items": len(rows), "engine_items": sum(1 for r in rows if r.get("engine") is not None)}
               for lab, rows in lab_rows.items()}
    if len(labs) != 2 or any(per_lab[lab]["engine"] is None for lab in labs):
        return {"per_lab": per_lab, "verdict": "not runnable"}
    a, b = labs
    meas = per_lab[a]["measured"] - per_lab[b]["measured"]
    pred = per_lab[a]["engine"] - per_lab[b]["engine"]
    ok = pred != 0 and math.copysign(1, pred) == math.copysign(1, meas) and abs(pred - meas) < abs(meas)
    return {"per_lab": per_lab, "measured_difference": meas, "engine_difference": pred,
            "verdict": "explained by conditions" if ok else "not explained"}


def score_payload() -> Dict[str, Any]:
    data = build()
    runner = S5A.Runner()
    readings = {}
    for reading in S5A.READINGS:
        cells, reasons = [], Counter()
        for c in data["cells"]:
            lab_rows: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
            for it in c["items"]:
                state, specs = S5A.classify(it, reading)
                row = {**_ident(it), "engine": None, "why": None}
                if specs is None:
                    row["why"] = state
                else:
                    p, why = S5A.engine_prediction(it, specs, runner)
                    row["engine"], row["why"] = p, (None if p is not None else "the engine refuses the ratio: " + why)
                if row["why"]:
                    reasons[S5A._why_class(row["why"])] += 1
                lab_rows[it["lab"]].append(row)
            cells.append({"cell": c["cell"], **judge(lab_rows), "rows": [r for rs in lab_rows.values() for r in rs]})
        readings[reading] = {"cells": cells, "verdicts": dict(Counter(c["verdict"] for c in cells)),
                             "not_run": dict(reasons.most_common())}
    return {"artifact": SCORE_OUT.stem, "status": "diagnostic (step 5c); no display, no credit",
            "primary_reading": S5A.PRIMARY, "readings": readings, "engine_integrations": runner.integrations}


def render_items(p: Mapping[str, Any]) -> str:
    lines = ["# Step 5c items: the conflict cells and how each contrast maps to the engine", "",
             "Generated by `scripts/generators/build_step5c_discordant_cells.py --items`. No engine output.", "",
             "| reading | mapping (contrasts or pairs) |", "|---|---|"]
    for r, counts in p["mapping_counts"].items():
        lines.append(f"| {r} | " + "; ".join(f"{k}: {n}" for k, n in counts.items()) + " |")
    lines += ["", "| cell | laboratories | items |", "|---|---|---|"]
    for c in p["cells"]:
        lines.append(f"| {' x '.join(c['cell'])} | {', '.join(c['laboratories'])} | {len(c['items'])} |")
    return "\n".join(lines) + "\n"


def _f(x: Optional[float]) -> str:
    return "--" if x is None else f"{x:+.2f}"


def render_score(p: Mapping[str, Any]) -> str:
    lines = ["# Step 5c: the engine on the discordant cells", "",
             "Generated by `scripts/generators/build_step5c_discordant_cells.py --score`; rules in its docstring. "
             "A diagnostic: nothing is displayed or credited. Log ratio per lever step, dex.", ""]
    for r, rd in p["readings"].items():
        tag = " (primary)" if r == p["primary_reading"] else ""
        lines += [f"## Reading: {r}{tag}", "",
                  "Verdicts: " + ", ".join(f"{k} {n}" for k, n in rd["verdicts"].items()) + ".", "",
                  "Not run: " + ("; ".join(f"{k} ({n})" for k, n in rd["not_run"].items()) or "none") + ".", "",
                  "| cell | laboratory: measured / engine (items run of items) | verdict |", "|---|---|---|"]
        for c in rd["cells"]:
            labs = "; ".join(f"{lab}: {_f(v['measured'])} / {_f(v['engine'])} ({v['engine_items']} of {v['items']})"
                             for lab, v in sorted(c["per_lab"].items()))
            lines.append(f"| {' x '.join(c['cell'])} | {labs} | {c['verdict']} |")
        lines.append("")
    lines.append(f"Engine integrations: {p['engine_integrations']}.")
    return "\n".join(lines) + "\n"


def _inputs() -> List[Path]:
    return [AUDIT, Path(S5A.__file__), Path(__file__)]


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--items", action="store_true")
    mode.add_argument("--score", action="store_true")
    args = ap.parse_args(argv)
    payload, out, render = ((items_payload(), ITEMS_OUT, render_items) if args.items
                            else (score_payload(), SCORE_OUT, render_score))
    payload["provenance"] = provenance.provenance_block(
        out.stem, generated_by="scripts/generators/build_step5c_discordant_cells.py",
        inputs=_inputs() + ([ITEMS_OUT] if args.score else []))
    artifact_io.write_artifact(payload, out, render=render)
    print(f"wrote {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
