"""
Sweep the registered choices and report whether any headline depends on them (2026-10-28).

For each entry of `docs/reference/choices_register.yml` that can be swept on the shipped engine, the headlines
are recomputed under every alternative: the out-of-sample rows within the band and their median fold error,
the same per lane, and the within-study contrast diagnostic. A choice MOVES A HEADLINE when any of those
integer counts differs between its alternatives; no tolerance is needed, so no new number is introduced.
Entries that need a refit (the thermal readings, the shared barrier) are listed as pending with the wave that
sweeps them. Writes results/validation/choice_robustness.{json,md}.

    python scripts/generators/build_choice_robustness.py
"""
from __future__ import annotations

import dataclasses
import json
import math
import statistics
import sys
from pathlib import Path
from typing import Any, Dict, List

import yaml

ROOT = Path(__file__).resolve().parents[2]
for p in (ROOT, ROOT / "scripts" / "generators"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import build_contrast_scores as CONTRASTS  # noqa: E402
import build_null_model_baselines as NULL  # noqa: E402
from src import artifact_io, data_paths, provenance  # noqa: E402
from src.kinetic_core import parameters_lipid  # noqa: E402

REGISTER = ROOT / "docs" / "reference" / "choices_register.yml"
OUT = data_paths.VALIDATION_DIR / "choice_robustness.json"
PENDING = {"CR-03": "B64 (both readings refitted)", "CR-04": "B64 sensitivity (refit constants, rescored)",
           "CR-07": "B63"}


def _oos_rows(scores: Dict[str, Any]) -> List[tuple]:
    """(lane, fold error) of every out-of-sample row, by the scorecard's own definition (scoring.score_panel)."""
    out = []
    for b in scores["benchmarks"]:
        literature = b["signal_origin"] == "external_literature" and b["evidence_role"] != "fit_recovery"
        for r in b["compounds"]:
            if literature and not r["in_core_fit"] and r.get("fold_error") is not None:
                out.append((str(r["lane"]), float(r["fold_error"])))
    return out


def _headlines(scores: Dict[str, Any], band: float) -> Dict[str, Any]:
    oos = _oos_rows(scores)
    by_lane: Dict[str, int] = {}
    for lane, fold in oos:
        by_lane[lane] = by_lane.get(lane, 0) + int(fold < band)
    pairs = [p for p in CONTRASTS._pairs(NULL._rows(scores)[0]) if p["out_of_sample"]]
    contrast = CONTRASTS._score(pairs)
    return {"oos_rows": len(oos), "oos_within": sum(f < band for _, f in oos),
            "oos_median_fold": statistics.median(f for _, f in oos) if oos else None, "by_lane_within": by_lane,
            "contrast_median_dex": contrast.get("engine", {}).get("median_error_dex"),
            "contrast_signs": contrast.get("engine", {}).get("sign_agreement_on_resolved"),
            "contrast_resolved": contrast.get("resolved_pairs")}


def _moves(results: Dict[str, Dict[str, Any]]) -> bool:
    keys = ("oos_within", "by_lane_within", "contrast_signs")
    return any(len({repr(v[k]) for v in results.values()}) > 1 for k in keys)


def _rescored() -> Dict[str, Any]:
    """The scorecard recomputed live on the current engine state."""
    from src.kinetic_core.scoring import score_panel

    return score_panel()


def main() -> int:
    register = yaml.safe_load(REGISTER.read_text(encoding="utf-8"))["entries"]
    by_id = {e["id"]: e for e in register}
    shipped = json.loads(NULL.SCORES.read_text(encoding="utf-8"))
    report: Dict[str, Any] = {}

    band_alts = by_id["CR-05"]["alternatives"]
    report["CR-05"] = {a["name"]: _headlines(shipped, float(a["value"])) for a in band_alts}

    noise = {}
    saved_noise = CONTRASTS.NOISE_DEX
    for a in by_id["CR-06"]["alternatives"]:
        CONTRASTS.NOISE_DEX = math.log10(float(a["value"]))
        noise[a["name"]] = _headlines(shipped, 3.0)
    CONTRASTS.NOISE_DEX = saved_noise
    report["CR-06"] = noise

    q10 = {}
    saved_q10 = parameters_lipid.Q10_ASSUMPTION
    try:
        for a in by_id["CR-02"]["alternatives"]:
            parameters_lipid.Q10_ASSUMPTION = dataclasses.replace(saved_q10, default=float(a["value"]))
            q10[a["name"]] = _headlines(_rescored(), 3.0)
    finally:
        parameters_lipid.Q10_ASSUMPTION = saved_q10
    report["CR-02"] = q10

    verdicts = {}
    for entry in register:
        key = entry["id"]
        if key in report:
            verdicts[key] = "moves_headline" if _moves(report[key]) else "harmless"
        elif key in PENDING:
            verdicts[key] = f"pending: {PENDING[key]}"
        else:
            verdicts[key] = entry["sensitivity"]["status"]
    payload = {
        "artifact": "choice_robustness",
        "provenance": provenance.provenance_block(
            "choice_robustness", generated_by="scripts/generators/build_choice_robustness.py",
            inputs=[REGISTER, NULL.SCORES]),
        "definition": __doc__.split("    python")[0].strip(),
        "verdicts": verdicts, "sweeps": report,
    }
    artifact_io.write_artifact(payload, OUT, render=render)
    for key, v in verdicts.items():
        print(key, v)
    return 0


def render(payload: Dict[str, Any]) -> str:
    lines = ["# Choice robustness", "", "*Generated by `scripts/generators/build_choice_robustness.py` from "
             "`docs/reference/choices_register.yml`.*", "", payload["definition"], "", "## Verdicts", "",
             "| entry | verdict |", "| --- | --- |"]
    lines += [f"| {k} | {v} |" for k, v in payload["verdicts"].items()]
    for key, sweep in payload["sweeps"].items():
        lines += ["", f"## {key}", "", "| alternative | out of sample within band | median fold | by lane | "
                  "contrast median (dex) | contrast signs |", "| --- | --- | --- | --- | --- | --- |"]
        for name, h in sweep.items():
            lanes = ", ".join(f"{k} {v}" for k, v in sorted(h["by_lane_within"].items()))
            med = f"{h['oos_median_fold']:.2f}x" if h["oos_median_fold"] else "-"
            cm = f"{h['contrast_median_dex']:.2f}" if h["contrast_median_dex"] is not None else "-"
            lines.append(f"| {name} | {h['oos_within']}/{h['oos_rows']} | {med} | {lanes} | {cm} | "
                         f"{h['contrast_signs']}/{h['contrast_resolved']} |")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
