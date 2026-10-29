#!/usr/bin/env python
"""
Does the mechanism pay for itself? A NULL-MODEL BASELINE against the kinetic core (2026-10-28, north star C1).

`tasks/data_restructure_plan.md` sec. 7, "North star", target 6: on the same rows the scorecard scores, the engine
must beat the simplest model that knows no chemistry: a compound's TYPICAL YIELD. For each row, the median of
y = log10(measured ppb / limiting precursor mM) over the same compound in the OTHER studies (leave-one-study-out:
every bundle from the held-out paper leaves together), times the row's limiting precursor. Deliberately simple
(owner, 2026-10-28: "let's not overcomplicate the statistical analysis"); a log-linear regression on T, pH and
time was tried first and did worse than this.

A row the baseline cannot answer (its compound appears in no other study) is counted separately, as the
engine's refusals are. The comparison is on the rows all three answer, which are all ppb rows; ratio rows are
excluded and counted. The baseline is given MORE than the engine: it trains on hold-out rows of other studies,
which the engine never read. That asymmetry favours the baseline, so a win for the engine is conservative.

    python scripts/generators/build_null_model_baselines.py
    # writes results/validation/null_model_baselines.{json,md}
"""
from __future__ import annotations

import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import artifact_io, compound_keys, data_paths, paper_keys, provenance  # noqa: E402
from src.kinetic_core.panel import limiting_precursor_molar, load_bundle  # noqa: E402

SCORES = data_paths.VALIDATION_DIR / "core_panel_scores.json"
OUT = data_paths.VALIDATION_DIR / "null_model_baselines.json"
BAND = 3.0


def _study(bench: Dict[str, Any], bench_file: str) -> str:
    doi = bench.get("source_doi")
    if doi:
        key = paper_keys.for_doi(doi)
        return key.paper_id if key else (paper_keys.normalise_doi(doi) or doi)
    citation = (bench.get("source_metadata") or {}).get("citation", "")
    return citation[:40] or Path(bench_file).stem


def _rows(scores: Optional[Dict[str, Any]] = None) -> Tuple[List[Dict[str, Any]], Dict[str, int]]:
    """The scored ppb rows of ``scores`` (default: the tracked scorecard)."""
    if scores is None:
        scores = json.loads(SCORES.read_text(encoding="utf-8"))
    rows, excluded = [], defaultdict(int)
    for b in scores["benchmarks"]:
        if not b["compounds"]:
            continue
        bench = load_bundle(ROOT / b["bench_file"])
        _, limiting = limiting_precursor_molar(bench)
        c = b["conditions"]
        literature = b["signal_origin"] == "external_literature" and b["evidence_role"] != "fit_recovery"
        for r in b["compounds"]:
            if r["target_unit"] != "ppb":
                excluded["ratio_or_non_ppb_rows"] += 1
                continue
            if not limiting or not r["measured"] or r["measured"] <= 0 or None in (c.get("temp_C"), c.get("time_min"), c.get("ph")):
                excluded["no_limiting_precursor_or_conditions"] += 1
                continue
            rows.append({
                "benchmark": b["benchmark_id"], "study": _study(bench, b["bench_file"]),
                "compound": compound_keys.canonical_name(r["compound"]), "lane": r["lane"],
                "T": float(c["temp_C"]), "pH": float(c["ph"]), "logt": math.log10(float(c["time_min"])),
                "limiting_mM": float(limiting) * 1000.0,  # limiting_precursor_molar returns mol/L
                "measured": float(r["measured"]), "engine": r["predicted"],
                "out_of_sample": bool(literature and not r["in_core_fit"]), "literature": bool(literature),
            })
    return rows, dict(excluded)


def _y(r: Dict[str, Any]) -> float:
    return math.log10(r["measured"] / r["limiting_mM"])


def _constant(train: List[Dict[str, Any]], r: Dict[str, Any]) -> Optional[float]:
    ys = [_y(t) for t in train if t["compound"] == r["compound"]]
    return 10.0 ** statistics.median(ys) * r["limiting_mM"] if ys else None



def _score(rows: List[Dict[str, Any]], key: str) -> Dict[str, Any]:
    logs = [abs(math.log10(r[key] / r["measured"])) for r in rows]
    return {"rows": len(rows), "within_3x": sum(v <= math.log10(BAND) for v in logs),
            "median_fold": 10.0 ** statistics.median(logs) if logs else None,
            "median_abs_log10": statistics.median(logs) if logs else None}



def main() -> int:
    rows, excluded = _rows()
    for r in rows:
        train = [t for t in rows if t["study"] != r["study"]]
        r["constant"] = _constant(train, r)
    answered = [r for r in rows if r["engine"] and r["engine"] > 0 and r["constant"]]
    unanswered = len(rows) - len(answered)
    out: Dict[str, Any] = {}
    for label, subset in (("all_ppb_rows", answered),
                          ("out_of_sample", [r for r in answered if r["out_of_sample"]])):
        out[label] = {m: _score(subset, m) for m in ("engine", "constant")}
    by_lane: Dict[str, Any] = {}
    for lane in sorted({r["lane"] for r in answered}):
        sub = [r for r in answered if r["lane"] == lane]
        by_lane[lane] = {m: _score(sub, m) for m in ("engine", "constant")}
    payload = {
        "artifact": "null_model_baselines",
        "provenance": provenance.provenance_block(
            "null_model_baselines", generated_by="scripts/generators/build_null_model_baselines.py", inputs=[SCORES]),
        "definition": __doc__.split("    python")[0].strip(),
        "band": BAND, "studies": len({r["study"] for r in rows}),
        "excluded": excluded, "rows_no_baseline_answer": unanswered,
        "scores": out, "by_lane": by_lane,
        "rows": [{k: r[k] for k in ("benchmark", "study", "compound", "lane", "measured", "engine", "constant",
                                     "out_of_sample")} for r in rows],
    }
    artifact_io.write_artifact(payload, OUT, render=render)
    o = out["out_of_sample"]
    print(f"out of sample: engine {o['engine']['within_3x']}/{o['engine']['rows']} (median {o['engine']['median_fold']:.2f}x) | "
          f"typical yield {o['constant']['within_3x']} ({o['constant']['median_fold']:.2f}x)")
    return 0


def render(payload: Dict[str, Any]) -> str:
    lines = ["# Null-model baselines against the kinetic core", "",
             f"*Generated by `scripts/generators/build_null_model_baselines.py` from `core_panel_scores.json`. "
             f"Leave-one-study-out over {payload['studies']} studies; within {payload['band']:g}x. "
             f"Excluded: {payload['excluded']}; rows no baseline can answer: {payload['rows_no_baseline_answer']}.*", "",
             payload["definition"], ""]
    for label, s in payload["scores"].items():
        lines += [f"## {label.replace('_', ' ')}", "", "| model | rows | within 3x | median fold |", "| --- | --- | --- | --- |"]
        for m in ("engine", "constant"):
            v = s[m]
            lines.append(f"| {'typical yield' if m == 'constant' else m} | {v['rows']} | {v['within_3x']} | {v['median_fold']:.2f}x |")
        lines.append("")
    lines += ["## By lane (all ppb rows)", "", "| lane | rows | engine within 3x (median) | typical yield |", "| --- | --- | --- | --- |"]
    for lane, s in payload["by_lane"].items():
        cell = lambda v: f"{v['within_3x']} ({v['median_fold']:.1f}x)"  # noqa: E731
        lines.append(f"| {lane} | {s['engine']['rows']} | {cell(s['engine'])} | {cell(s['constant'])} |")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
