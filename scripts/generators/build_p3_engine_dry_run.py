"""P3 engine dry run (STATUS_AND_PLAN step 16e; FACTORIAL_LEVER_STUDY.md, 'Premortem').

Before the P3 study is sent: on how many of the study's lever x odorant cells can the shipped engine make a
prediction at all? Two laboratories are the default only because step 5 (the engine judged on held-out contrasts
from two TEST laboratories) needs two; if the engine predicts few of the cells, that reason is void and the design
falls back to one laboratory (2^3 x 4).

For each of the 8 corners (sugar glucose / ribose 30 mM x sulfur source cysteine / thiamine 10 mM x pH 6.0 / 7.0;
140 C, 20 min) and each analyte of the design, the shipped engine is asked through its own front door
(``declare_envelope``, then ``compare``) under three readings:

* ``pea_isolate``: 15 % w/w pea isolate through ``src.ingredients`` (protein sites only: no dossier gives its free
  sugars or amino acids; this is the placeholder composition), no buffer salt declared (the isolate is the buffer);
* ``pea_isolate_ph_held``: the same with the pH held at the set value (the engine's clamped buffer), a sensitivity
  for pH drift;
* ``added_precursors``: the added sugar and sulfur source in water, no isolate (what step 5a could run).

A lever's prediction is A4's main effect: the mean, over the 4 matched pairs (the other two factors fixed), of the
engine's log10 ratio high / low. A cell is **predicted** when all 4 pairs return a ratio, **partial** when some do.
The buffered pair (glucose + cysteine, 0.1 M phosphate vs none, pH 6.0 and 7.0) is reported beside: the engine's pH
slope with and without phosphate.

A diagnostic for a design decision: it changes no constant, ships nothing and reads no study data. Integrations are
cached per (formulation, target).
"""

from __future__ import annotations

import argparse
import itertools
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "generators"))

import build_lever_design_preposterior as PP  # noqa: E402
import build_step5a_display_check as S5  # noqa: E402
from src import artifact_io, data_paths, ingredients, provenance  # noqa: E402
from src.kinetic_core import engine  # noqa: E402
from src.kinetic_core.ph_state import BufferSpec  # noqa: E402

ARTIFACT = "p3_engine_dry_run"
OUT_DIR = ROOT / "results" / "validation"

TEMP_C, MINUTES = 140.0, 20.0
ISOLATE_PERCENT = 15.0
SUGARS = ("glucose", "ribose")          # reference first
SULFUR = ("cysteine", "thiamine")
PH = (6.0, 7.0)
SUGAR_MM, SULFUR_MM = 30.0, 10.0
PHOSPHATE_M = 0.1
LEVERS = ("sugar identity", "sulfur source", "pH")
READINGS = ("pea_isolate", "pea_isolate_ph_held", "added_precursors")
PRIMARY = "pea_isolate"
ANALYTES = list(PP.ANALYTES)


def corner_spec(sugar: str, sulfur: str, ph: float, reading: str, *, phosphate: bool = False
                ) -> engine.FormulationSpec:
    """The engine run of one P3 pot under ``reading``."""
    precursors = {sugar: SUGAR_MM, sulfur: SULFUR_MM}
    matrix, protein = "water", None
    if reading.startswith("pea_isolate"):
        resolved, _ = ingredients.expand({"ingredients": {"pea_isolate": ISOLATE_PERCENT},
                                          "precursors": precursors})
        precursors = dict(resolved["precursors"])
        matrix, protein = resolved["matrix"], resolved["protein_g_per_l"]
    if phosphate:
        buffer = BufferSpec(kind="phosphate", phosphate_mol_l=PHOSPHATE_M, declared=True,
                            source="P3 buffered pair: 0.1 M phosphate")
    elif reading.endswith("ph_held"):
        buffer = BufferSpec(kind="clamped", declared=True, source="P3 dry run: pH held at the set value")
    else:
        buffer = BufferSpec(kind="none", declared=True, source="P3: no buffer salt (the isolate is the buffer)")
    process = engine.ProcessSpec(thermal=engine.ThermalProgram.isothermal(TEMP_C, MINUTES), ph=ph, buffer=buffer,
                                 matrix=matrix, protein_g_per_l=protein)
    name = f"{reading}:{sugar}/{sulfur}/pH{ph}{'/phosphate' if phosphate else ''}"
    return engine.FormulationSpec(name=name, precursors=precursors, process=process)


def _pairs(lever: str) -> List[Tuple[Tuple[str, str, float], Tuple[str, str, float]]]:
    """The 4 matched (low, high) corner pairs of one lever."""
    out = []
    for a, b in itertools.product((0, 1), repeat=2):
        if lever == "sugar identity":
            out.append(((SUGARS[0], SULFUR[a], PH[b]), (SUGARS[1], SULFUR[a], PH[b])))
        elif lever == "sulfur source":
            out.append(((SUGARS[a], SULFUR[0], PH[b]), (SUGARS[a], SULFUR[1], PH[b])))
        else:
            out.append(((SUGARS[a], SULFUR[b], PH[0]), (SUGARS[a], SULFUR[b], PH[1])))
    return out


def _refusal(spec: engine.FormulationSpec, target: str) -> Optional[str]:
    decl = engine.declare_envelope(spec, [target])
    if decl.is_answerable:
        return None
    return "the engine refuses the arm: " + (decl.reasons[0] if decl.reasons else decl.state)[:160]


def _why_class(why: str) -> str:
    return why.split(":")[0].strip()


def cell(compound: str, lever: str, reading: str, runner: S5.Runner) -> Dict[str, Any]:
    target = S5.target_alias(compound)
    row: Dict[str, Any] = {"compound": compound, "lever": lever, "reading": reading, "target": target}
    if target is None:
        return {**row, "status": "not predicted", "why": ["no engine target for this odorant"]}
    logs, whys = [], []
    for lo, hi in _pairs(lever):
        a, b = corner_spec(*hi, reading), corner_spec(*lo, reading)
        why = _refusal(a, target) or _refusal(b, target)
        if why is None:
            r, why = runner.log_ratio(a, b, target)
            if r is not None:
                logs.append(r)
                continue
            why = "the engine refuses the ratio: " + why
        whys.append(why)
    status = "predicted" if len(logs) == 4 else ("partial" if logs else "not predicted")
    row.update(status=status, pairs_predicted=len(logs), why=sorted(set(whys)))
    if logs:
        row["main_effect_dex"] = sum(logs) / len(logs)
        row["pair_log_ratios"] = logs
        row["corner_signs_agree"] = len({x > 0 for x in logs if abs(x) > S5.FLOOR_DEX}) <= 1
    return row


def buffered_pair(compound: str, runner: S5.Runner) -> Dict[str, Any]:
    """The engine's pH slope (dex per pH unit) at glucose + cysteine in the isolate, with and without phosphate."""
    target = S5.target_alias(compound)
    out: Dict[str, Any] = {"compound": compound}
    if target is None:
        return {**out, "why": "no engine target for this odorant"}
    for label, phos in (("matrix_only", False), ("phosphate", True)):
        lo = corner_spec(SUGARS[0], SULFUR[0], PH[0], PRIMARY, phosphate=phos)
        hi = corner_spec(SUGARS[0], SULFUR[0], PH[1], PRIMARY, phosphate=phos)
        why = _refusal(hi, target) or _refusal(lo, target)
        r = None
        if why is None:
            r, why = runner.log_ratio(hi, lo, target)
        out[label] = (r / (PH[1] - PH[0])) if r is not None else None
        if r is None:
            out[f"{label}_why"] = why
    if out.get("matrix_only") is not None and out.get("phosphate") is not None:
        out["difference"] = out["phosphate"] - out["matrix_only"]
    return out


def payload() -> Dict[str, Any]:
    runner = S5.Runner()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cells = [cell(c, lv, rd, runner) for rd in READINGS for lv in LEVERS for c in ANALYTES]
        pair = [buffered_pair(c, runner) for c in ANALYTES]
    summary: Dict[str, Any] = {}
    for rd in READINGS:
        rows = [r for r in cells if r["reading"] == rd]
        reasons: Dict[str, int] = {}
        for r in rows:
            for w in r["why"]:
                reasons[_why_class(w)] = reasons.get(_why_class(w), 0) + 1
        summary[rd] = {"cells": len(rows),
                       "predicted": sum(r["status"] == "predicted" for r in rows),
                       "partial": sum(r["status"] == "partial" for r in rows),
                       "by_lever": {lv: sum(r["status"] == "predicted" for r in rows if r["lever"] == lv)
                                    for lv in LEVERS},
                       "reasons": dict(sorted(reasons.items(), key=lambda kv: -kv[1]))}
    return {"artifact": ARTIFACT,
            "design": {"temp_C": TEMP_C, "minutes": MINUTES, "isolate_percent_w_w": ISOLATE_PERCENT,
                       "sugars_mM": {s: SUGAR_MM for s in SUGARS}, "sulfur_mM": {s: SULFUR_MM for s in SULFUR},
                       "ph": list(PH), "phosphate_M": PHOSPHATE_M},
            "engine": {"fingerprint": getattr(engine, "ENGINE_FINGERPRINT", None)},
            "integrations": runner.integrations, "primary_reading": PRIMARY,
            "summary": summary, "cells": cells, "buffered_pair": pair}


def _f(x: Optional[float]) -> str:
    return "-" if x is None else f"{x:+.2f}"


def render(p: Mapping[str, Any]) -> str:
    s = p["summary"]
    prim = s[PRIMARY]
    lines = [f"# P3 engine dry run: which of the study's cells can the shipped engine predict? (step 16e)", "",
             "*Generated by `scripts/generators/build_p3_engine_dry_run.py` (STATUS_AND_PLAN step 16e; "
             "`docs/guides/FACTORIAL_LEVER_STUDY.md`, 'Premortem'). The shipped engine on the 8 corners of the P3 "
             "design through its own front door; A4 main effects. A design diagnostic: it changes no constant, "
             "ships nothing and reads no study data.*", "",
             f"**Result.** In 15 % pea isolate (protein sites only, the placeholder composition), the engine predicts "
             f"**{prim['predicted']} of {prim['cells']}** lever x odorant cells (all 4 matched pairs), and "
             f"{prim['partial']} more in part. {p['integrations']} integrations.", "",
             "| reading | cells | predicted | partial | sugar identity | sulfur source | pH | main reasons not predicted |",
             "|---|---|---|---|---|---|---|---|"]
    for rd in READINGS:
        r = s[rd]
        reasons = "; ".join(f"{k} ({v} pairs)" for k, v in list(r["reasons"].items())[:3]) or "-"
        lines.append(f"| {rd} | {r['cells']} | {r['predicted']} | {r['partial']} | "
                     f"{r['by_lever']['sugar identity']} | {r['by_lever']['sulfur source']} | {r['by_lever']['pH']} "
                     f"| {reasons} |")
    lines += ["", "Readings: `pea_isolate` = 15 % pea isolate via `src/ingredients.py`, no buffer salt; "
              "`pea_isolate_ph_held` = the same, pH held at the set value; `added_precursors` = the added sugar "
              "and sulfur source in water. 'Reasons' count matched pairs, not cells.", "",
              f"## Cells, reading `{PRIMARY}`", "",
              "| odorant | engine target | sugar: ribose vs glucose | sulfur: thiamine vs cysteine | pH: 7.0 vs 6.0 |",
              "|---|---|---|---|---|"]
    by = {(r["compound"], r["lever"]): r for r in p["cells"] if r["reading"] == PRIMARY}
    for c in ANALYTES:
        def show(lv: str) -> str:
            r = by[(c, lv)]
            if r["status"] == "not predicted":
                return "not predicted: " + (_why_class(r["why"][0]) if r["why"] else "")
            tag = "" if r["status"] == "predicted" else f" ({r['pairs_predicted']}/4 pairs)"
            flag = "" if r.get("corner_signs_agree", True) else ", corners disagree"
            return f"{_f(r['main_effect_dex'])}{tag}{flag}"
        lines.append(f"| {c} | {by[(c, LEVERS[0])]['target'] or '-'} | {show(LEVERS[0])} | {show(LEVERS[1])} | "
                     f"{show(LEVERS[2])} |")
    lines += ["", "## The buffered pair: the engine's pH slope (dex per pH unit), glucose + cysteine in the isolate", "",
              "| odorant | matrix only | 0.1 M phosphate | difference |", "|---|---|---|---|"]
    for r in p["buffered_pair"]:
        if "why" in r:
            continue
        lines.append(f"| {r['compound']} | {_f(r.get('matrix_only'))} | {_f(r.get('phosphate'))} | "
                     f"{_f(r.get('difference'))} |")
    lines += ["", "Values are the engine's log10 ratios (dex), the A4 main effect (mean of the matched pairs). "
              "'Corners disagree': the matched pairs differ in sign (the A7 corner check would flag the cell).", ""]
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    p = payload()
    p["provenance"] = provenance.provenance_block(
        ARTIFACT, generated_by="scripts/generators/build_p3_engine_dry_run.py",
        inputs=[data_paths.INGREDIENTS])
    artifact_io.write_artifact(p, OUT_DIR / f"{ARTIFACT}.json", render=render)
    print(f"wrote {ARTIFACT}: {p['summary'][PRIMARY]['predicted']} of {p['summary'][PRIMARY]['cells']} predicted")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
