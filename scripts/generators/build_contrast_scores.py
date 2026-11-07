"""
Within-study contrasts: how well the engine predicts the RATIO between two pots of one study (D4, 2026-10-28).

Owner, 2026-10-28: "when we make predictions, those should be lab independent ... we want general
predictions". A ratio between two pots of one study cancels everything the study holds fixed, its
calibration and its unreported conditions included, so it tests the chemistry with the laboratory
removed. It is a DIAGNOSTIC: nothing is fitted on it and no laboratory term exists anywhere.

For every pair of scored ppb rows of the same compound in the same study: measured m = log10(a / b) and
predicted p = log10(engine a / engine b). Reported: the median |p - m| (dex), the share within 3x, the
sign agreement on pairs whose measured effect clears the measurement noise (|m| > log10 1.4, the
directional panel's MEASUREMENT_NOISE_FOLD), and the same numbers for the no-change baseline (p = 0).
Pairs from one study share rows, so they are not independent; the study count is printed beside them.

    python scripts/generators/build_contrast_scores.py
"""
from __future__ import annotations

import itertools
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Tuple

ROOT = Path(__file__).resolve().parents[2]
for p in (ROOT, ROOT / "scripts" / "generators"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import build_null_model_baselines as NULL  # noqa: E402  (the same rows and the same study key)
from src import artifact_io, data_paths, provenance  # noqa: E402

OUT = data_paths.VALIDATION_DIR / "contrast_scores.json"
BAND_DEX = math.log10(3.0)
#: pairs refused by S3 (the pH balance; filled by ``_pairs``)
BALANCE_REFUSED_PAIRS: List[Dict[str, Any]] = []
NOISE_DEX = math.log10(1.4)


def _absolute(r: Dict[str, Any]) -> bool:
    """True for a row with a shipped absolute level; False for a relative-only row (centred log level only)."""
    return bool(r.get("engine")) and r["engine"] > 0


def _balance_refused(a: Dict[str, Any], b: Dict[str, Any], engine_dex: float) -> bool:
    """2026-11-07 (S3 application of M1b part (i)): the pair's ratio is refused when B0 and the carried pH-balance
    rivals disagree on its sign or spread it beyond omega (``ph_balance_refusal.offsets_ratio_refusal``)."""
    from src.kinetic_core.ph_balance_refusal import offsets_ratio_refusal

    return offsets_ratio_refusal(engine_dex, a.get("balance_s3_offsets_dex"), b.get("balance_s3_offsets_dex")) is not None


def _pairs(rows: List[Dict[str, Any]]) -> Tuple[List[Dict[str, Any]], int]:
    """The within-study pairs, and the number of MIXED pairs skipped.

    A mixed pair holds one row with an absolute level and one with only a centred log level (2026-11-05, H0
    follow-up), e.g. MGO measured in a trunk pot and in a sulfur pot of one study. The centred value is relative
    to the mean of the study's refused rows only, so no ratio with an absolute level exists: the pair is skipped
    and counted, never scored."""
    BALANCE_REFUSED_PAIRS.clear()
    groups: Dict[tuple, List[Dict[str, Any]]] = defaultdict(list)
    for r in rows:
        if _absolute(r) or r.get("engine_centred_log10") is not None:
            groups[(r["study"], r["compound"])].append(r)
    out, mixed = [], 0
    for (study, compound), members in sorted(groups.items()):
        for a, b in itertools.combinations(sorted(members, key=lambda r: r["benchmark"]), 2):
            if a["benchmark"] == b["benchmark"]:
                continue
            if _absolute(a) != _absolute(b):
                mixed += 1
                continue
            # 2026-11-05 (H0 / EX-059): a thiol-lane row carries its centred log level, not a ppb level; the
            # difference of two centred values IS the log ratio (one group, one mean)
            engine_dex = (math.log10(a["engine"] / b["engine"]) if _absolute(a)
                          else a["engine_centred_log10"] - b["engine_centred_log10"])
            if _balance_refused(a, b, engine_dex):
                BALANCE_REFUSED_PAIRS.append({"study": study, "compound": compound, "a": a["benchmark"],
                                              "b": b["benchmark"]})
                continue
            out.append({
                "study": study, "compound": compound, "lane": a["lane"], "a": a["benchmark"], "b": b["benchmark"],
                "measured_dex": math.log10(a["measured"] / b["measured"]),
                "engine_dex": engine_dex,
                "out_of_sample": a["out_of_sample"] and b["out_of_sample"],
            })
    return out, mixed


def _score(pairs: List[Dict[str, Any]]) -> Dict[str, Any]:
    if not pairs:
        return {"pairs": 0, "studies": 0}
    resolved = [p for p in pairs if abs(p["measured_dex"]) > NOISE_DEX]

    def stats(predict):
        errors = [abs(predict(p) - p["measured_dex"]) for p in pairs]
        signs = [math.copysign(1.0, predict(p)) == math.copysign(1.0, p["measured_dex"]) and predict(p) != 0.0
                 for p in resolved]
        return {"median_error_dex": statistics.median(errors),
                "within_3x": sum(e < BAND_DEX for e in errors),
                "sign_agreement_on_resolved": sum(signs)}

    return {"pairs": len(pairs), "studies": len({p["study"] for p in pairs}), "resolved_pairs": len(resolved),
            "engine": stats(lambda p: p["engine_dex"]), "no_change": stats(lambda p: 0.0)}


def main() -> int:
    rows, _excluded = NULL._rows(include_relative_only=True)   # 2026-11-05: ratios of the thiol lane ship (H0 / EX-059)
    pairs, mixed = _pairs(rows)
    scores = {"all_pairs": _score(pairs), "out_of_sample": _score([p for p in pairs if p["out_of_sample"]])}
    by_lane = {lane: _score([p for p in pairs if p["lane"] == lane]) for lane in sorted({p["lane"] for p in pairs})}
    payload = {
        "artifact": "contrast_scores",
        "provenance": provenance.provenance_block(
            "contrast_scores", generated_by="scripts/generators/build_contrast_scores.py", inputs=[NULL.SCORES]),
        "definition": __doc__.split("    python")[0].strip(),
        "scores": scores, "by_lane": by_lane,
        # one absolute row and one relative-only row: no ratio exists, skipped (2026-11-05, H0 follow-up)
        "mixed_pairs_skipped": mixed,
        # 2026-11-07 (S3, ph_balance_refusal): pairs whose ratio the falsified pH balance moves beyond omega or in sign
        "balance_refused_pairs": list(BALANCE_REFUSED_PAIRS),
        "pairs": pairs,
    }
    artifact_io.write_artifact(payload, OUT, render=render)
    o = scores["out_of_sample"]
    if o["pairs"]:
        print(f"out of sample: {o['pairs']} pairs in {o['studies']} studies; engine median {o['engine']['median_error_dex']:.2f} dex, "
              f"no-change {o['no_change']['median_error_dex']:.2f} dex; signs {o['engine']['sign_agreement_on_resolved']}"
              f"/{o['resolved_pairs']}")
    return 0


def render(payload: Dict[str, Any]) -> str:
    lines = ["# Within-study contrasts", "",
             "*Generated by `scripts/generators/build_contrast_scores.py` from `core_panel_scores.json`. "
             "A diagnostic: nothing is fitted on it.*", "", payload["definition"], "",
             "| subset | pairs (studies) | resolved | engine median error | engine within 3x | engine signs | "
             "no-change median error | no-change within 3x |", "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for label, s in list(payload["scores"].items()) + [(f"lane: {k}", v) for k, v in payload["by_lane"].items()]:
        if not s["pairs"]:
            lines.append(f"| {label} | 0 | | | | | | |")
            continue
        e, n = s["engine"], s["no_change"]
        lines.append(f"| {label.replace('_', ' ')} | {s['pairs']} ({s['studies']}) | {s['resolved_pairs']} | "
                     f"{e['median_error_dex']:.2f} dex | {e['within_3x']} | {e['sign_agreement_on_resolved']}/{s['resolved_pairs']} | "
                     f"{n['median_error_dex']:.2f} dex | {n['within_3x']} |")
    lines += ["", "Mixed pairs skipped (one row with an absolute level, one relative-only; no ratio exists): "
              f"{payload.get('mixed_pairs_skipped', 0)}."]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
