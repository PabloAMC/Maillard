#!/usr/bin/env python3
"""Freeze the CURRENT model's predictions on the MAILLARD-PATH hold-out.

WHY THIS EXISTS
---------------
Wave S1 (2026-08-27) established that all four bundles in the repo's external
hold-out run the ``matrix_only`` execution path: they read a lipid-oxidation
load off a matrix profile and never call ``Recommender.predict_from_steps``.
Three consecutive waves of reaction-network changes -- an additive flux
propagator, a pH/water-activity rewiring, and a shipped barrier revert -- moved
in-panel predictions and left all eight hold-out points BIT-IDENTICAL. The
reaction network, which is the scientific content of this repository, had never
been tested out of sample.

This script scores the FREE-PRECURSOR hold-out bundles under
``data/benchmarks/external_validation/maillard_path/`` -- systems whose
predictions are produced by ``predict_from_steps`` -- and writes the result to
``results/validation/maillard_path_holdout_frozen_predictions.{json,md}``.

THE ORDERING IS THE POINT. This artifact is a PRE-REGISTRATION. It records the
git HEAD it was generated from, BEFORE any rate-calibration wave has seen these
points. A later calibration wave must be scored against THIS file; regenerating
it after a calibration and comparing to itself proves nothing.

UNITS
-----
Papers that report a yield in mol % are kept in mol %. The repo learned in Wave
S2b that converting a literature mol % into ppb requires a molar basis, and that
an assumed basis becomes a free multiplicative parameter hiding underneath the
score (the 342/200 lesson). So:

  * ``target_unit == "ppb"``   -> scored directly against the model's ppb output.
  * ``target_unit == "mol_percent"`` -> the model's ppb output is converted to
    mol % by the model's OWN declared projection identity and nothing else:

        ppb = molar_conc(mol/L) * MW(g/mol) * ppb_conversion_factor(=1e6)

    (``src/projection.py``: ``DEFAULT_PROJECTION_STRATEGY.ppb_conversion_factor
    = 1.0e6``, i.e. the ppb basis is ug/L read as ug/kg at 1 kg/L), hence

        predicted_mol_percent = 100 * (ppb / 1e6 / MW) / limiting_precursor_M

    where ``limiting_precursor_M`` is the LOWEST precursor molarity declared in
    the bundle -- the same "limiting precursor" rule
    ``src/projection._estimate_projection_budget`` uses -- and the bundle must
    declare, in ``molar_basis``, that this matches the denominator the PAPER
    quotes. No other basis is introduced. If a bundle's paper defines its mol %
    against something the bundle cannot reproduce, the bundle sets
    ``unit_commensurable: false`` and the pair is scored ORDINALLY ONLY.

NOTHING IN THIS SCRIPT MAY EVER BE FITTED. The bundles are
``evidence_class: external_validation_only`` and are asserted so by
``scripts/ci/holdout_guard.py``.

Usage:
    python scripts/generators/generate_maillard_path_holdout_frozen_predictions.py
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import subprocess
import sys
from datetime import date
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

MAILLARD_PATH_HOLDOUT_DIR = ROOT / "data" / "benchmarks" / "external_validation" / "maillard_path"
DEFAULT_BASENAME = "maillard_path_holdout_frozen_predictions"

# The projection identity, quoted from src/projection.py so a drift is visible here.
PPB_CONVERSION_FACTOR = 1.0e6

# Basis-free ratio units a paper may report a yield in, and the multiplier that
# turns a bare mole fraction (mol product / mol limiting precursor) into them.
# These are the ONLY non-ppb target units this artifact scores quantitatively,
# precisely because they need no assumed molar basis (the 342/200 lesson).
RATIO_UNIT_FACTORS = {
    "mol_percent": 100.0,
    "umol_per_mol_limiting_precursor": 1.0e6,
}


def _git_head() -> Dict[str, str]:
    def _run(*args: str) -> str:
        try:
            return subprocess.run(
                ["git", *args], cwd=ROOT, capture_output=True, text=True, check=True
            ).stdout.strip()
        except Exception:
            return "unknown"

    return {
        "commit": _run("rev-parse", "HEAD"),
        "short": _run("rev-parse", "--short", "HEAD"),
        "branch": _run("rev-parse", "--abbrev-ref", "HEAD"),
        "dirty": "yes" if _run("status", "--porcelain") else "no",
    }


def _molecular_weight(smiles: str) -> Optional[float]:
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except Exception:  # pragma: no cover - rdkit is a hard dep of the model anyway
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return float(Descriptors.MolWt(mol))


def _limiting_precursor_molar(bench: Dict[str, Any]) -> Tuple[Optional[str], Optional[float]]:
    """Lowest declared precursor molarity, in mol/L.

    This mirrors ``src/projection._estimate_projection_budget``'s limiting-precursor
    rule (``min`` over positive concentrations) so the mol % denominator the score
    uses is the same one the model's own budget uses.
    """
    items = [
        (str(name), float(data.get("concentration_mM", 0.0)) / 1000.0)
        for name, data in (bench.get("precursors") or {}).items()
        if float(data.get("concentration_mM", 0.0)) > 0.0
    ]
    if not items:
        return None, None
    return min(items, key=lambda kv: kv[1])


def _load_bundles(directory: Path) -> List[Path]:
    return sorted(directory.glob("*.json"))


def build_frozen_predictions(
    *,
    directory: Path = MAILLARD_PATH_HOLDOUT_DIR,
    target_tag: str = "meaty",
) -> Dict[str, Any]:
    from src.benchmark_validation import evaluate_benchmark, load_benchmark

    bundles = _load_bundles(directory)
    if not bundles:
        raise SystemExit(f"No maillard-path hold-out bundles found under {directory}")

    rows: List[Dict[str, Any]] = []
    commensurable_folds: List[float] = []
    commensurable_logs: List[float] = []

    for path in bundles:
        bench = load_benchmark(path)
        evaluation = evaluate_benchmark(path, target_tag=target_tag)
        limiting_name, limiting_molar = _limiting_precursor_molar(bench)
        targets = bench.get("holdout_targets") or {}

        compound_rows: List[Dict[str, Any]] = []
        for comparison in evaluation.comparisons:
            spec = targets.get(comparison.compound, {})
            unit = str(spec.get("target_unit", "ppb"))
            commensurable = bool(spec.get("unit_commensurable", unit == "ppb"))
            predicted_ppb = float(comparison.predicted_ppb)

            predicted_native: Optional[float] = None
            if unit == "ppb":
                predicted_native = predicted_ppb
            elif unit in RATIO_UNIT_FACTORS:
                smiles = spec.get("smiles")
                mw = _molecular_weight(str(smiles)) if smiles else None
                if mw and limiting_molar:
                    predicted_native = (
                        RATIO_UNIT_FACTORS[unit]
                        * (predicted_ppb / PPB_CONVERSION_FACTOR / mw)
                        / limiting_molar
                    )
                else:
                    commensurable = False
            else:
                commensurable = False

            target_value = spec.get("target_value")
            if target_value is None:
                target_value = float((bench.get("reference_volatiles") or {}).get(
                    comparison.compound, {}
                ).get("conc_ppb", 0.0))

            fold_error: Optional[float] = None
            abs_log10: Optional[float] = None
            if (
                commensurable
                and predicted_native is not None
                and predicted_native > 0.0
                and float(target_value) > 0.0
            ):
                ratio = predicted_native / float(target_value)
                fold_error = max(ratio, 1.0 / ratio)
                abs_log10 = abs(math.log10(ratio))
                commensurable_folds.append(fold_error)
                commensurable_logs.append(abs_log10)

            compound_rows.append(
                {
                    "compound": comparison.compound,
                    "matched_prediction_key": comparison.matched_name,
                    "target_unit": unit,
                    "target_value": None if target_value is None else float(target_value),
                    "target_quote": spec.get("target_quote"),
                    "molar_basis": spec.get("molar_basis"),
                    "unit_commensurable": commensurable,
                    "predicted_ppb": predicted_ppb,
                    "predicted_native_unit": predicted_native,
                    "fold_error": fold_error,
                    "abs_log10_error": abs_log10,
                    "scoring": _scoring_label(fold_error, commensurable, predicted_ppb, target_value),
                }
            )

        # Ordinal check: does the model rank the measured compounds in the paper's order?
        ordinal = _score_ordinal(compound_rows)

        rows.append(
            {
                "benchmark_id": bench.get("benchmark_id"),
                "bench_file": str(path.relative_to(ROOT)),
                "source_doi": bench.get("source_doi"),
                "citation": (bench.get("source_metadata") or {}).get("citation"),
                "access_route": (bench.get("content_verification") or {}).get("access_route"),
                "access_level": (bench.get("content_verification") or {}).get("access_level"),
                "retrieval_date": (bench.get("content_verification") or {}).get("retrieval_date"),
                "execution_path": (bench.get("metadata") or {}).get("execution_path"),
                "precursors": bench.get("precursors"),
                "conditions": bench.get("conditions"),
                "limiting_precursor": limiting_name,
                "limiting_precursor_molar": limiting_molar,
                "compounds": compound_rows,
                "ordinal": ordinal,
            }
        )

    series = _score_series(rows)
    within_10x = sum(1 for f in commensurable_folds if f <= 10.0)
    summary = {
        "bundle_count": len(rows),
        "target_count": sum(len(r["compounds"]) for r in rows),
        "quantitatively_scored_count": len(commensurable_folds),
        "ordinal_only_count": sum(
            1 for r in rows for c in r["compounds"] if c["scoring"] == "ordinal_only"
        ),
        "structural_zero_count": sum(
            1 for r in rows for c in r["compounds"] if c["scoring"] == "structural_zero"
        ),
        "median_fold_error": statistics.median(commensurable_folds) if commensurable_folds else None,
        "worst_fold_error": max(commensurable_folds) if commensurable_folds else None,
        "best_fold_error": min(commensurable_folds) if commensurable_folds else None,
        "median_abs_log10_error": statistics.median(commensurable_logs) if commensurable_logs else None,
        "within_10x": within_10x,
        "within_10x_rate": (within_10x / len(commensurable_folds)) if commensurable_folds else None,
        "ordinal_pairs_correct": sum(r["ordinal"]["correct_pairs"] for r in rows),
        "ordinal_pairs_total": sum(r["ordinal"]["total_pairs"] for r in rows),
        "series_directions_correct": sum(s["directions_correct"] for s in series),
        "series_directions_total": sum(s["directions_total"] for s in series),
    }

    return {
        "artifact": DEFAULT_BASENAME,
        "generated_on": date.today().isoformat(),
        "git": _git_head(),
        "target_tag": target_tag,
        "pre_registration": (
            "These predictions were generated BEFORE any rate-calibration wave saw these "
            "points. A calibration wave must be scored against THIS file at THIS git commit. "
            "Regenerating this artifact after a calibration and comparing it to itself proves "
            "nothing."
        ),
        "benchmarks": rows,
        "series": series,
        "summary": summary,
    }


# --- Cross-bundle series ---------------------------------------------------
# A per-point fold error cannot see whether the model gets the SHAPE right. Three
# of the harvested sources vary exactly one axis at a time and hold everything
# else fixed, so their bundles form series whose measured ratio can be compared
# with the model's predicted ratio. This is where an uncalibrated branch shows
# what it actually got wrong, and it is deliberately computed here rather than
# left to prose.
DECLARED_SERIES = (
    {
        "series_id": "sulfur_temperature_series_Yiltirak2026",
        "axis": "temperature up / time down: 100 C 4 h -> 130 C 0.5 h",
        "low": "mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026",
        "high": "mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026",
        "compounds": ("2-Methyl-3-furanthiol (MFT)", "2-Furfurylthiol (FFT)"),
    },
    {
        "series_id": "acrylamide_time_series_Chang2021",
        "axis": "time at fixed 180 C: 10 min -> 30 min",
        "low": "mp_holdout_glucose_asparagine_180C_10min_Chang2021",
        "high": "mp_holdout_glucose_asparagine_180C_30min_Chang2021",
        "compounds": ("Acrylamide",),
    },
    {
        "series_id": "ph_pair_Schibilsky2019",
        "axis": "pH at fixed 130 C / 2 h: pH 5 -> pH 8",
        "low": "mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019",
        "high": "mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019",
        "compounds": ("Furfural", "DMHF", "5-Hydroxymethylfurfural (HMF)"),
    },
)


def _score_series(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    by_id = {str(row["benchmark_id"]): row for row in rows}
    out: List[Dict[str, Any]] = []
    for spec in DECLARED_SERIES:
        low, high = by_id.get(spec["low"]), by_id.get(spec["high"])
        if low is None or high is None:
            continue
        entries: List[Dict[str, Any]] = []
        for compound in spec["compounds"]:
            lo = next((c for c in low["compounds"] if c["compound"] == compound), None)
            hi = next((c for c in high["compounds"] if c["compound"] == compound), None)
            if lo is None or hi is None:
                continue
            m_lo, m_hi = lo["target_value"], hi["target_value"]
            p_lo, p_hi = lo["predicted_ppb"], hi["predicted_ppb"]
            if not (m_lo and m_hi and p_lo and p_hi):
                continue
            measured_ratio = float(m_hi) / float(m_lo)
            predicted_ratio = float(p_hi) / float(p_lo)
            direction_ok = (measured_ratio > 1.0) == (predicted_ratio > 1.0)
            entries.append(
                {
                    "compound": compound,
                    "measured_ratio_high_over_low": measured_ratio,
                    "predicted_ratio_high_over_low": predicted_ratio,
                    "direction_correct": direction_ok,
                    "response_ratio": predicted_ratio / measured_ratio,
                }
            )
        out.append({**{k: v for k, v in spec.items() if k != "compounds"}, "entries": entries,
                    "directions_correct": sum(1 for e in entries if e["direction_correct"]),
                    "directions_total": len(entries)})
    return out


def _scoring_label(
    fold_error: Optional[float],
    commensurable: bool,
    predicted_ppb: float,
    target_value: Optional[float],
) -> str:
    """Name what kind of score a target got, so no failure hides inside a median.

    ``structural_zero`` is the case the fold-error median would otherwise swallow:
    the model returned literally nothing for a compound the source measured. It has
    no finite fold error, so it cannot enter the median -- which is exactly why it
    has to be counted and named separately rather than dropped as "not evaluable".
    """
    if fold_error is not None:
        return "quantitative"
    if commensurable and predicted_ppb <= 0.0 and target_value and float(target_value) > 0.0:
        return "structural_zero"
    return "ordinal_only"


def _score_ordinal(compound_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Pairwise ordinal agreement between measured and predicted magnitudes.

    Only pairs whose measured values share a unit are compared, so a mol %
    target is never ranked against a ppb target.
    """
    pairs: List[Dict[str, Any]] = []
    correct = 0
    for i in range(len(compound_rows)):
        for j in range(i + 1, len(compound_rows)):
            a, b = compound_rows[i], compound_rows[j]
            if a["target_unit"] != b["target_unit"]:
                continue
            if a["target_value"] is None or b["target_value"] is None:
                continue
            if a["predicted_ppb"] <= 0.0 and b["predicted_ppb"] <= 0.0:
                continue
            measured_order = a["target_value"] > b["target_value"]
            predicted_order = a["predicted_ppb"] > b["predicted_ppb"]
            hit = measured_order == predicted_order
            correct += int(hit)
            pairs.append(
                {
                    "a": a["compound"],
                    "b": b["compound"],
                    "measured_a_gt_b": measured_order,
                    "predicted_a_gt_b": predicted_order,
                    "hit": hit,
                }
            )
    return {"pairs": pairs, "correct_pairs": correct, "total_pairs": len(pairs)}


def _fmt(value: Optional[float], digits: int = 4) -> str:
    if value is None:
        return "n/a"
    if value == 0:
        return "0"
    if abs(value) >= 1e5 or abs(value) < 1e-3:
        return f"{value:.3e}"
    return f"{value:.{digits}g}"


def render_markdown(payload: Dict[str, Any]) -> str:
    git = payload["git"]
    summary = payload["summary"]
    lines: List[str] = []
    lines.append(
        "<!-- Auto-generated by "
        "scripts/generators/generate_maillard_path_holdout_frozen_predictions.py. "
        "This file is a PRE-REGISTRATION: do not regenerate it to make a later "
        "calibration look good. -->"
    )
    lines.append("")
    lines.append("# Maillard-path hold-out — FROZEN predictions")
    lines.append("")
    lines.append(
        f"**Generated from git HEAD `{git['short']}` ({git['commit']}), branch "
        f"`{git['branch']}`, working tree dirty: {git['dirty']}, on {payload['generated_on']}.**"
    )
    lines.append("")
    lines.append(
        "> **A later rate-calibration wave MUST be scored against this file, at this "
        "commit.** These are the predictions the model made on these systems *before* any "
        "calibration wave had seen them. This is the repository's first out-of-sample test "
        "of the reaction network itself: every bundle below runs the `free_precursor` "
        "execution path, i.e. `Recommender.predict_from_steps`. The pre-existing eight-point "
        "external hold-out runs `matrix_only` and never reaches the network (Wave S1)."
    )
    lines.append("")
    lines.append(
        "> **These points were not fitted, and may not be fitted.** Every bundle carries "
        "`evidence_class: external_validation_only` and lives under "
        "`data/benchmarks/external_validation/`, which `scripts/ci/holdout_guard.py` "
        "asserts statically. They are deliberately kept in the `maillard_path/` "
        "subdirectory so they do NOT enter `external_validation_report`'s median: that "
        "median is a matrix-only lipid-oxidation number and averaging two different "
        "execution paths into one figure would be meaningless."
    )
    lines.append("")
    lines.append("## Baseline scores (the number S3 has to beat)")
    lines.append("")
    lines.append("| Metric | Value |")
    lines.append("|---|---|")
    lines.append(f"| Bundles | {summary['bundle_count']} |")
    lines.append(f"| Targets | {summary['target_count']} |")
    lines.append(f"| Quantitatively scored (unit-commensurable) | {summary['quantitatively_scored_count']} |")
    lines.append(f"| Ordinal-only (unit-incommensurable) | {summary['ordinal_only_count']} |")
    lines.append(
        f"| **Structural zeroes** (model returned nothing; excluded from every fold statistic) | "
        f"**{summary['structural_zero_count']}** |"
    )
    lines.append(f"| **Median fold error** | **{_fmt(summary['median_fold_error'])}x** |")
    lines.append(f"| Worst fold error | {_fmt(summary['worst_fold_error'])}x |")
    lines.append(f"| Best fold error | {_fmt(summary['best_fold_error'])}x |")
    lines.append(f"| Median abs log10 error | {_fmt(summary['median_abs_log10_error'])} dex |")
    lines.append(
        f"| Within 10x | {summary['within_10x']}/{summary['quantitatively_scored_count']} |"
    )
    lines.append(
        f"| Ordinal pairs correct | {summary['ordinal_pairs_correct']}/{summary['ordinal_pairs_total']} |"
    )
    lines.append("")

    if payload.get("series"):
        lines.append("## Cross-bundle series — does the model get the SHAPE right?")
        lines.append("")
        lines.append(
            "Three of the harvested sources vary exactly one axis and hold everything else fixed. "
            "A per-point fold error cannot see whether the model's RESPONSE to that axis is right; "
            "these rows can. `response ratio` is predicted-change divided by measured-change: 1.0 is "
            "a perfect response, below 1 is under-response, above 1 is over-response, and a negative "
            "direction verdict means the model moved the opposite way to the measurement."
        )
        lines.append("")
        lines.append("| Series | Axis | Compound | Measured x | Predicted x | Response ratio | Direction |")
        lines.append("|---|---|---|---|---|---|---|")
        for group in payload["series"]:
            for entry in group["entries"]:
                verdict = "OK" if entry["direction_correct"] else "**WRONG WAY**"
                lines.append(
                    f"| `{group['series_id']}` | {group['axis']} | {entry['compound']} | "
                    f"{_fmt(entry['measured_ratio_high_over_low'], 3)}x | "
                    f"{_fmt(entry['predicted_ratio_high_over_low'], 3)}x | "
                    f"{_fmt(entry['response_ratio'], 3)} | {verdict} |"
                )
        lines.append("")

    lines.append("## Per-point record")
    lines.append("")
    for row in payload["benchmarks"]:
        lines.append(f"### `{row['benchmark_id']}`")
        lines.append("")
        lines.append(f"* **Citation:** {row.get('citation')}")
        lines.append(f"* **DOI:** `{row.get('source_doi')}`")
        lines.append(
            f"* **Access:** {row.get('access_level')} via {row.get('access_route')} "
            f"(retrieved {row.get('retrieval_date')})"
        )
        precursors = ", ".join(
            f"{name} {data.get('concentration_mM')} mM"
            for name, data in (row.get("precursors") or {}).items()
        )
        conditions = row.get("conditions") or {}
        lines.append(f"* **System:** {precursors}")
        lines.append(
            f"* **Conditions:** {conditions.get('temp_C')} C, {conditions.get('time_min')} min, "
            f"pH {conditions.get('ph')}, aw {conditions.get('water_activity')}"
        )
        lines.append(
            f"* **Limiting precursor (mol % denominator):** {row.get('limiting_precursor')} "
            f"= {_fmt(row.get('limiting_precursor_molar'))} mol/L"
        )
        lines.append("")
        lines.append("| Compound | Target (native) | Unit | Predicted (native) | Predicted ppb | Fold error | Scoring |")
        lines.append("|---|---|---|---|---|---|---|")
        for compound in row["compounds"]:
            fold = f"{_fmt(compound['fold_error'])}x" if compound["fold_error"] is not None else "—"
            lines.append(
                f"| {compound['compound']} | {_fmt(compound['target_value'])} | "
                f"{compound['target_unit']} | {_fmt(compound['predicted_native_unit'])} | "
                f"{_fmt(compound['predicted_ppb'])} | {fold} | {compound['scoring']} |"
            )
        lines.append("")
        for compound in row["compounds"]:
            if compound.get("target_quote"):
                lines.append(f"> **{compound['compound']} — quoted source:** {compound['target_quote']}")
                lines.append("")
        ordinal = row["ordinal"]
        if ordinal["total_pairs"]:
            lines.append(
                f"**Ordinal:** {ordinal['correct_pairs']}/{ordinal['total_pairs']} pairs correct."
            )
            for pair in ordinal["pairs"]:
                verdict = "OK" if pair["hit"] else "MISS"
                lines.append(
                    f"  * {verdict}: measured {pair['a']} "
                    f"{'>' if pair['measured_a_gt_b'] else '<'} {pair['b']}; "
                    f"predicted {'>' if pair['predicted_a_gt_b'] else '<'}"
                )
            lines.append("")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="results/validation")
    parser.add_argument("--basename", default=DEFAULT_BASENAME)
    parser.add_argument("--target-tag", default="meaty")
    args = parser.parse_args()

    payload = build_frozen_predictions(target_tag=args.target_tag)
    out_dir = ROOT / args.output_dir if not Path(args.output_dir).is_absolute() else Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / f"{args.basename}.json"
    md_path = out_dir / f"{args.basename}.md"
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    md_path.write_text(render_markdown(payload), encoding="utf-8")

    summary = payload["summary"]
    print(json.dumps({"json": str(json_path), "md": str(md_path)}, indent=2))
    print(
        f"Maillard-path hold-out: {summary['bundle_count']} bundles, "
        f"{summary['quantitatively_scored_count']} scored, median "
        f"{_fmt(summary['median_fold_error'])}x, within-10x "
        f"{summary['within_10x']}/{summary['quantitatively_scored_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
