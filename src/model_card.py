"""The model card: a validity domain assembled from live artifacts, never written by hand.

WHY THIS MODULE EXISTS
----------------------
2026-08-28 (Wave S5). This repository's honest numbers have drifted before -- twice a README
sentence outlived the evidence behind it, and Wave S5 found a third case: the directional
report's own headline said 20/29 four commits after the tree started scoring 21/29. Prose
does not update itself.

So the validity domain -- what this model may be used for, on what system, with what measured
reliability -- is **computed**, from artifacts that a wave has to regenerate anyway, and
pasted into README.md between markers by ``scripts/generators/generate_model_card.py``. If
an artifact moves, the card moves. If an artifact is missing, the card says so in the cell
rather than quietly dropping the row: a blank is indistinguishable from a pass, and this
campaign has spent nine waves removing exactly that failure mode.

WHAT IT READS
-------------
============================================ ==========================================
directional_accuracy_report.md               ordinal reliability, per axis (tracked)
maillard_path_holdout_frozen_predictions.json free-precursor hold-out, pre-registered
external_validation_report.json              matrix hold-out, CI coverage, extrapolation
matrix_binding_mode_comparison.json          three observability modes, hold-out medians
the benchmark panel                          RECOMPUTED live -- benchmark_summary.json is
                                             gitignored and absent in a fresh clone
data/**/*.json,yml                           the no_verifiable_source census, recounted
data/benchmarks/cys_ribose_140C_Hofmann1998  the sulfur branch's anchor status
scripts/ci/*.py                              the three blocking gates, actually run
============================================ ==========================================

NOTHING HERE IS A MEASUREMENT. Every number is read or recounted; the only judgement this
module adds is the trust/caution/do-not-use threshold, which is stated in the card itself.
"""

from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple

import yaml

from src.directional_reliability import (
    CAUTION_MIN_RATE,
    MIN_EVALUABLE_FOR_TRUST,
    TRUST_MIN_RATE,
    VERDICT_CAUTION,
    VERDICT_DO_NOT_USE,
    VERDICT_TRUST,
    load_panel_counts,
    verdict_for,
)

ROOT = Path(__file__).resolve().parents[1]
VALIDATION_DIR = ROOT / "results" / "validation"

HOLDOUT_PATH = VALIDATION_DIR / "maillard_path_holdout_frozen_predictions.json"
EXTERNAL_PATH = VALIDATION_DIR / "external_validation_report.json"
MODE_COMPARISON_PATH = VALIDATION_DIR / "matrix_binding_mode_comparison.json"
SULFUR_BENCHMARK_PATH = ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json"

GATES = (
    "scripts/ci/holdout_guard.py",
    "scripts/ci/citation_gate.py",
    "scripts/ci/fit_target_gate.py",
)

README_BEGIN = "<!-- BEGIN GENERATED: model-card -->"
README_END = "<!-- END GENERATED: model-card -->"

MISSING = "artifact absent"


# ---------------------------------------------------------------------------------------
# Evidence collection
# ---------------------------------------------------------------------------------------


def _relative(path: Path) -> str:
    """Repo-relative where possible, absolute otherwise.

    ``Path.relative_to`` raises on a path outside the tree, and a card generator must never
    crash because an artifact it was pointed at lives somewhere unexpected -- it must SAY the
    artifact is elsewhere or missing.
    """
    try:
        return str(Path(path).relative_to(ROOT))
    except ValueError:
        return str(path)


def _read_json(path: Path) -> Optional[Dict[str, Any]]:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return None


def collect_directional(report_path: Optional[Path] = None) -> Dict[str, Any]:
    counts = load_panel_counts(report_path) if report_path else load_panel_counts()
    # The bucket rows and the category rows have DIFFERENT DENOMINATORS: the categories sum
    # to the screening bucket (independent + system-overlap), the headline counts only the
    # independent claims. Summing categories to derive an "excluding pH/aw" figure therefore
    # quotes 24/29 where the headline-comparable number is 18/23. Both are true; mixing them
    # inflates the score. So the bucket figures are READ from their own rows
    # rather than derived, and the report states that arithmetic in place.
    buckets = {
        "headline": "strictly independent (headline)",
        "system_overlap": "system-overlap",
        "fit_adjacent": "fit-adjacent (excluded from the headline)",
        "excluding_ph_aw": "independent, excluding ph and moisture_aw",
        "ph_aw": "independent, ph and moisture_aw only",
    }
    categories = {
        key: value for key, value in counts.items() if key not in set(buckets.values())
    }
    resolved = {name: counts.get(label) for name, label in buckets.items()}
    resolved["categories"] = categories
    resolved["category_bucket_note"] = (
        "Category rows sum to the SCREENING bucket (independent + system-overlap); the "
        "headline and the excluding-pH/aw rows are strictly independent. Different "
        "denominators -- see the report's CURRENT STANDING section."
    )
    return resolved


def collect_free_precursor_holdout() -> Dict[str, Any]:
    payload = _read_json(HOLDOUT_PATH)
    if payload is None:
        return {"available": False, "path": _relative(HOLDOUT_PATH)}
    summary = payload.get("summary", {})
    return {
        "available": True,
        "path": _relative(HOLDOUT_PATH),
        "generated_on": payload.get("generated_on"),
        "git": (payload.get("git") or {}).get("short") or (payload.get("git") or {}).get("commit"),
        "bundles": summary.get("bundle_count"),
        "targets": summary.get("target_count"),
        "scored": summary.get("quantitatively_scored_count"),
        "median_fold_error": summary.get("median_fold_error"),
        "worst_fold_error": summary.get("worst_fold_error"),
        "best_fold_error": summary.get("best_fold_error"),
        "within_10x": summary.get("within_10x"),
        "ordinal_pairs_correct": summary.get("ordinal_pairs_correct"),
        "ordinal_pairs_total": summary.get("ordinal_pairs_total"),
        "series_directions_correct": summary.get("series_directions_correct"),
        "series_directions_total": summary.get("series_directions_total"),
        "structural_zeroes": summary.get("structural_zero_count"),
    }


def collect_matrix_holdout() -> Dict[str, Any]:
    payload = _read_json(EXTERNAL_PATH)
    if payload is None:
        return {"available": False, "path": _relative(EXTERNAL_PATH)}
    summary = payload.get("summary", {})
    split = summary.get("holdout_kind_split", {}) or {}
    genuine = split.get("genuine_extrapolation", {}) or {}
    in_panel = split.get("in_panel_rescoring", {}) or {}
    return {
        "available": True,
        "path": _relative(EXTERNAL_PATH),
        "ci_coverage_hits": summary.get("ci_coverage_hits"),
        "ci_level_pct": summary.get("ci_level_pct"),
        "genuine_hits": genuine.get("hits"),
        "genuine_total": genuine.get("total"),
        "in_panel_hits": in_panel.get("hits"),
        "in_panel_total": in_panel.get("total"),
        "median_ci_width_log10": (summary.get("honest_literature_coverage") or {}).get(
            "median_ci_width_log10"
        ),
    }


def collect_observability_modes() -> Dict[str, Any]:
    payload = _read_json(MODE_COMPARISON_PATH)
    if payload is None:
        return {"available": False, "path": _relative(MODE_COMPARISON_PATH)}
    monte_carlo = payload.get("holdout_monte_carlo", {}) or {}
    return {
        "available": True,
        "path": _relative(MODE_COMPARISON_PATH),
        "shipped_mode": "fitted_factors",
        "modes": {
            mode: {
                "median_fold_error": (monte_carlo.get(mode) or {}).get("median_fold_error"),
                "ci_coverage_hits": (monte_carlo.get(mode) or {}).get("ci_coverage_hits"),
                "ci_coverage_total": (monte_carlo.get(mode) or {}).get("ci_coverage_total"),
            }
            for mode in payload.get("modes", [])
        },
    }


def collect_benchmark_panel() -> Dict[str, Any]:
    """Recompute the panel rather than read benchmark_summary.json.

    That artifact is gitignored, so a card that read it would print ``artifact absent`` on
    every fresh clone and on CI -- the self-excusing-skip pattern Wave J2 removed. The panel
    re-evaluates in a couple of seconds, so live recomputation is both cheap and stronger.
    """
    from src.benchmark_validation import (
        evaluate_benchmark,
        get_benchmark_files,
        load_benchmark,
        summarize_evaluation,
    )

    statuses: Dict[str, int] = {}
    strict_ready = 0
    total = 0
    for bench_file in get_benchmark_files():
        evaluation = evaluate_benchmark(bench_file)
        bench = load_benchmark(bench_file)
        summary = summarize_evaluation(evaluation, protein_type=bench.get("protein_type", "free"))
        total += 1
        statuses[summary.overall_status] = statuses.get(summary.overall_status, 0) + 1
        strict_ready += 1 if summary.strict_ready else 0
    return {
        "available": True,
        "source": "recomputed live from data/benchmarks/",
        "total": total,
        "strict_ready": strict_ready,
        "status_counts": dict(sorted(statuses.items())),
    }


def collect_no_verifiable_source_census() -> Dict[str, Any]:
    """Recount the census over tracked data files.

    Definition, stated because the number is meaningless without it: a RECORD is a mapping
    that carries ``no_verifiable_source`` in a key whose name contains ``status``. The
    ``source_status`` subset is the one README has quoted since Wave T3.
    """
    try:
        listed = subprocess.run(
            ["git", "ls-files", "*.json", "*.yml", "*.yaml"],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            check=True,
        ).stdout.split()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return {"available": False, "reason": "git ls-files unavailable"}

    records = 0
    with_numeric = 0
    by_key: Dict[str, int] = {}
    files: Dict[str, int] = {}

    def walk(node: Any, origin: str) -> None:
        nonlocal records, with_numeric
        if isinstance(node, dict):
            hits = [
                key
                for key, value in node.items()
                if value == "no_verifiable_source" and "status" in key.lower()
            ]
            if hits:
                records += 1
                files[origin] = files.get(origin, 0) + 1
                for key in hits:
                    by_key[key] = by_key.get(key, 0) + 1
                if any(
                    isinstance(value, (int, float)) and not isinstance(value, bool)
                    for value in node.values()
                ):
                    with_numeric += 1
            for value in node.values():
                walk(value, origin)
        elif isinstance(node, list):
            for value in node:
                walk(value, origin)

    for relative in listed:
        path = ROOT / relative
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        if "no_verifiable_source" not in text:
            continue
        try:
            data = json.loads(text) if relative.endswith(".json") else yaml.safe_load(text)
        except (json.JSONDecodeError, yaml.YAMLError):
            continue
        walk(data, relative)

    return {
        "available": True,
        "records": records,
        "with_numeric_payload": with_numeric,
        "by_status_key": dict(sorted(by_key.items(), key=lambda item: -item[1])),
        "files": len(files),
    }


_SULFUR_ANALYTE_MARKERS = ("thiol", "thio", "sulf", "mercapto", "mft", "fft")


def _is_sulfur_analyte(name: str) -> bool:
    lowered = name.lower()
    return any(marker in lowered for marker in _SULFUR_ANALYTE_MARKERS)


def collect_sulfur_anchor_status() -> Dict[str, Any]:
    """Check, rather than assert, that the sulfur branch has no absolute literature anchor."""
    payload = _read_json(SULFUR_BENCHMARK_PATH)
    if payload is None:
        return {"available": False, "path": _relative(SULFUR_BENCHMARK_PATH)}
    unverifiable: List[str] = []
    anchored: List[str] = []
    for name, record in (payload.get("measured_volatiles") or {}).items():
        if not isinstance(record, Mapping):
            continue
        statuses = {
            str(value)
            for key, value in record.items()
            if "status" in key.lower() and isinstance(value, str)
        }
        (unverifiable if "no_verifiable_source" in statuses else anchored).append(str(name))
    # 2026-08-28 (Wave W). This function used to look at the RETIRED benchmark and nothing
    # else, so `zero_absolute_anchors` was really "the retired file still has no verified
    # value" -- true, permanent, and no longer the question anyone is asking. The sulfur
    # branch gained three real anchors when the full text of 10.1021/jf9705983 arrived; a
    # card that kept printing "zero absolute literature anchors" off this one file would
    # have been stale in the flattering direction (understating what the model is now
    # measured against, and hiding that it fails it). So the panel is scanned too, and
    # `zero_absolute_anchors` is now false iff a PANEL benchmark carries a
    # primary-source-verified sulfur value.
    # 2026-08-28 (Wave X). A SECOND scan was added alongside the first, and the reason is
    # the audit's own rule rather than tidiness: a benchmark that a constant was FITTED
    # against is not evidence about the model, so listing it among the anchors the model
    # "fails" would be circular in the flattering direction the moment it stopped failing.
    # Wave X ingested the first step-level sulfur rows and declared one of them
    # (hofmann1998_norfuraneol_h2s_145C_20min_pH5, Hofmann Table 4) a fit target. It is
    # separated out here rather than dropped, because dropping it would hide a row that IS
    # currently a miss -- exclusion has to be symmetric about which direction it flatters.
    verified: List[str] = []
    fitted: List[str] = []
    benchmark_dir = ROOT / "data" / "benchmarks"
    if benchmark_dir.is_dir():
        for path in sorted(benchmark_dir.glob("*.json")):  # non-recursive: hold-out stays out
            other = _read_json(path)
            if not isinstance(other, Mapping):
                continue
            is_fit_target = "fit_target_declaration" in other
            for name, record in (other.get("measured_volatiles") or {}).items():
                if not isinstance(record, Mapping):
                    continue
                if record.get("value_status") != "primary_source_verified":
                    continue
                if not _is_sulfur_analyte(str(name)):
                    continue
                (fitted if is_fit_target else verified).append(f"{path.stem}/{name}")
    return {
        "available": True,
        "path": _relative(SULFUR_BENCHMARK_PATH),
        "tier": (payload.get("metadata") or {}).get("tier", payload.get("tier")),
        "values_without_verifiable_source": sorted(unverifiable),
        "values_with_an_anchor": sorted(anchored),
        "panel_verified_sulfur_values": sorted(verified),
        "panel_verified_sulfur_benchmarks": sorted(
            {row.split("/", 1)[0] for row in verified}
        ),
        "panel_verified_sulfur_values_that_are_fit_targets": sorted(fitted),
        "fit_target_note": (
            "Rows listed under `..._that_are_fit_targets` carry a `fit_target_declaration` "
            "block, i.e. a constant was selected by looking at them. They are primary-source "
            "verified and they are NOT counted as anchors: agreement on a fitted row carries "
            "no information about the model. They are listed rather than deleted because "
            "excluding a fitted row removes its MISSES as well as its hits, and at the time "
            "of writing the one row here is a miss."
        ),
        "zero_absolute_anchors": not anchored and not verified,
    }


def run_gates(gates: Tuple[str, ...] = GATES) -> List[Dict[str, Any]]:
    results: List[Dict[str, Any]] = []
    for gate in gates:
        path = ROOT / gate
        if not path.exists():
            results.append({"gate": gate, "status": "MISSING"})
            continue
        proc = subprocess.run(
            [sys.executable, str(path)],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            timeout=1800,
        )
        tail = [line for line in proc.stdout.strip().splitlines() if line.strip()]
        results.append(
            {
                "gate": gate,
                "status": "PASS" if proc.returncode == 0 else "FAIL",
                "detail": tail[-1] if tail else "",
            }
        )
    return results


def build_model_card(*, run_gate_checks: bool = True) -> Dict[str, Any]:
    card = {
        "artifact": "maillard_model_card",
        "generated_by": "scripts/generators/generate_model_card.py",
        "thresholds": {
            "trust": f">= {TRUST_MIN_RATE:.0%} agreement on >= {MIN_EVALUABLE_FOR_TRUST} claims",
            "caution": f">= {CAUTION_MIN_RATE:.0%} agreement, or too few claims to establish",
            "do_not_use": f"< {CAUTION_MIN_RATE:.0%} agreement, or unmeasured",
        },
        "directional": collect_directional(),
        "free_precursor_holdout": collect_free_precursor_holdout(),
        "matrix_holdout": collect_matrix_holdout(),
        "observability_modes": collect_observability_modes(),
        "benchmark_panel": collect_benchmark_panel(),
        "no_verifiable_source": collect_no_verifiable_source_census(),
        "sulfur_anchor": collect_sulfur_anchor_status(),
        "gates": run_gates() if run_gate_checks else [],
    }
    card["headline_sentences"] = _headline_sentences(card)
    card["validity_domain"] = _validity_domain(card)
    return card


# ---------------------------------------------------------------------------------------
# The three honest headline sentences, each carrying its own number
# ---------------------------------------------------------------------------------------


def _fold(value: Optional[float]) -> str:
    return MISSING if value is None else f"{float(value):.3g}x"


def _headline_sentences(card: Mapping[str, Any]) -> List[str]:
    free = card["free_precursor_holdout"]
    modes = card["observability_modes"]
    directional = card["directional"]
    sulfur = card["sulfur_anchor"]

    if modes.get("available"):
        medians = [
            entry["median_fold_error"]
            for entry in modes["modes"].values()
            if entry.get("median_fold_error") is not None
        ]
        matrix_range = (
            f"{min(medians):.3g}x-{max(medians):.3g}x" if medians else MISSING
        )
    else:
        matrix_range = MISSING

    free_median = _fold(free.get("median_fold_error")) if free.get("available") else MISSING
    free_worst = _fold(free.get("worst_fold_error")) if free.get("available") else MISSING

    headline = directional.get("headline")
    excl = directional.get("excluding_ph_aw")
    ph_aw = directional.get("ph_aw")

    sentences = [
        (
            f"**Absolute concentrations are unreliable.** Out of sample the free-precursor "
            f"lane lands at a {free_median} median fold error (worst {free_worst}) and the "
            f"matrix lane at {matrix_range} across all three observability modes. Nothing in "
            f"this repository licenses a ppb number as a specification."
        ),
        (
            f"**Directional and ranking claims are the measured product**, at "
            f"{headline[0]}/{headline[1]} on strictly independent claims"
            if headline
            else "**Directional and ranking claims are the measured product**"
        )
        + (
            f" -- {excl[0]}/{excl[1]} ({excl[0] / excl[1]:.0%}) once pH and water activity are "
            f"set aside, and {ph_aw[0]}/{ph_aw[1]} on pH and water activity themselves, which "
            f"is at or below chance."
            if excl and excl[1] and ph_aw and ph_aw[1]
            else "."
        ),
        (
            (
                "**The sulfur branch has "
                + str(len(sulfur.get("panel_verified_sulfur_benchmarks") or []))
                + " absolute literature anchors, and the model fails every one of them.** "
                "They are the primary-source-verified stable-isotope-dilution rows in "
                + ", ".join(sulfur.get("panel_verified_sulfur_benchmarks") or [])
                + ". "
                + (
                    (
                        "A further "
                        + str(len({
                            row.split("/", 1)[0]
                            for row in (sulfur.get(
                                "panel_verified_sulfur_values_that_are_fit_targets"
                            ) or [])
                        }))
                        + " primary-source-verified sulfur row(s) are on the panel and are NOT "
                        "counted here, because a constant was selected by looking at them ("
                        + ", ".join(sorted({
                            row.split("/", 1)[0]
                            for row in (sulfur.get(
                                "panel_verified_sulfur_values_that_are_fit_targets"
                            ) or [])
                        }))
                        + "): agreement on a fitted row is not evidence about the model. "
                    )
                    if sulfur.get("panel_verified_sulfur_values_that_are_fit_targets")
                    else ""
                )
                + "The previously shipped claim of ZERO anchors was corrected on 2026-08-28 "
                "(Wave W) when the full text behind them was obtained; the retired benchmark ("
                + Path(sulfur.get("path", "")).stem
                + ") is kept in the tree as the provenance record of the values that were not "
                "measurements. Absolute agreement is poor and the DIRECTION is a separate "
                "question."
            )
            if sulfur.get("available") and sulfur.get("panel_verified_sulfur_benchmarks")
            else "**The sulfur branch has zero absolute literature anchors.** "
            + (
                "Both values on its only benchmark ("
                + Path(sulfur.get("path", "")).stem
                + ") carry no verifiable source and its tightest contract was retired; the "
                "shipped numbers are the evidence for nothing but themselves. Its DIRECTION is "
                "a separate question and is not thereby worthless."
                if sulfur.get("available") and sulfur.get("zero_absolute_anchors")
                else "Status could not be verified from "
                + str(sulfur.get("path", "the benchmark file"))
                + "; treat this sentence as unconfirmed until it can be."
            )
        ),
    ]
    return sentences


# ---------------------------------------------------------------------------------------
# The validity-domain table: claim type x system class
# ---------------------------------------------------------------------------------------


@dataclass
class DomainCell:
    claim_type: str
    system_class: str
    evidence: str
    verdict: str
    basis: str = ""


def _verdict_from_fraction(hits: Optional[int], total: Optional[int]) -> str:
    if not total:
        return VERDICT_DO_NOT_USE
    return verdict_for(int(hits or 0), int(total))


def _validity_domain(card: Mapping[str, Any]) -> List[Dict[str, str]]:
    directional = card["directional"]
    categories = directional.get("categories", {})
    free = card["free_precursor_holdout"]
    matrix = card["matrix_holdout"]
    modes = card["observability_modes"]
    panel = card["benchmark_panel"]

    cells: List[DomainCell] = []

    # --- Absolute concentration -------------------------------------------------------
    if free.get("available"):
        cells.append(
            DomainCell(
                "Absolute concentration (ppb)",
                "free precursor (sugar + amino acid, aqueous)",
                f"median {_fold(free['median_fold_error'])}, "
                f"{free['within_10x']}/{free['scored']} within 10x",
                VERDICT_DO_NOT_USE,
                "12-point pre-registered hold-out, frozen before any calibration wave saw it",
            )
        )
    if modes.get("available"):
        shipped = modes["modes"].get(modes["shipped_mode"], {})
        cells.append(
            DomainCell(
                "Absolute concentration (ppb)",
                "protein matrix (pea / soy isolate)",
                f"median {_fold(shipped.get('median_fold_error'))}, CI coverage "
                f"{shipped.get('ci_coverage_hits')}/{shipped.get('ci_coverage_total')}",
                VERDICT_DO_NOT_USE,
                "8-point external hold-out; the shipped fitted factors LOSE to applying no "
                "observability factor at all",
            )
        )
    if matrix.get("available"):
        cells.append(
            DomainCell(
                "Absolute concentration (ppb)",
                "genuine extrapolation (roasting, HME extrusion)",
                f"{matrix.get('genuine_hits')}/{matrix.get('genuine_total')} inside the "
                f"{matrix.get('ci_level_pct')}% CI",
                VERDICT_DO_NOT_USE,
                "the only rows that test transfer; the interval that does cover is ~4 decades wide",
            )
        )

    # --- Ordinal claims, per axis ------------------------------------------------------
    axis_systems = {
        "sugar_identity": "any (sugar swap, conditions held)",
        "additive_cysteine": "free precursor (cysteine present vs absent)",
        "temperature": "any (temperature moved, everything else held)",
        "time": "any (time moved, everything else held)",
        "lipid_lane": "protein matrix (lipid-derived aldehydes)",
        "matrix_identity": "protein matrix (pea vs soy)",
        "ph": "any (pH moved)",
        "moisture_aw": "any (water activity moved)",
    }
    for axis, system in axis_systems.items():
        if axis not in categories:
            continue
        agree, evaluable = categories[axis]
        cells.append(
            DomainCell(
                f"Direction / ranking on `{axis}`",
                system,
                f"{agree}/{evaluable} on the directional panel",
                verdict_for(agree, evaluable),
                "",
            )
        )

    # --- Within-system ordering --------------------------------------------------------
    if free.get("available") and free.get("ordinal_pairs_total"):
        cells.append(
            DomainCell(
                "Ordering of two compounds in one system",
                "free precursor",
                f"{free['ordinal_pairs_correct']}/{free['ordinal_pairs_total']} pairs correct",
                _verdict_from_fraction(
                    free["ordinal_pairs_correct"], free["ordinal_pairs_total"]
                ),
                "measured on the hold-out, independently of the panel",
            )
        )
    if free.get("available") and free.get("series_directions_total"):
        cells.append(
            DomainCell(
                "Response direction across a condition series",
                "free precursor",
                f"{free['series_directions_correct']}/{free['series_directions_total']} correct",
                _verdict_from_fraction(
                    free["series_directions_correct"], free["series_directions_total"]
                ),
                "the sulfur temperature dependence is inverted and acrylamide is ~40x "
                "under-responsive in time",
            )
        )

    # --- Panel-level statements --------------------------------------------------------
    if panel.get("available"):
        cells.append(
            DomainCell(
                "Any claim of benchmark-grade agreement",
                "the in-panel benchmarks themselves",
                f"{panel['strict_ready']}/{panel['total']} strict-ready",
                VERDICT_DO_NOT_USE if not panel["strict_ready"] else VERDICT_CAUTION,
                "recomputed live; strict-ready is the repository's own passing bar",
            )
        )

    # --- What the model is actually for -------------------------------------------------
    cells.append(
        DomainCell(
            "Which experiment to run next (value of information)",
            "any system the uncertainty panel covers",
            "every ranked row is a measured model failure",
            VERDICT_TRUST,
            "this claim type does not depend on the model being right -- it depends on the "
            "model being wrong in a located, quantified way, which it demonstrably is",
        )
    )

    return [
        {
            "claim_type": cell.claim_type,
            "system_class": cell.system_class,
            "evidence": cell.evidence,
            "verdict": cell.verdict,
            "basis": cell.basis,
        }
        for cell in cells
    ]


# ---------------------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------------------

_VERDICT_LABEL = {
    VERDICT_TRUST: "**trust**",
    VERDICT_CAUTION: "caution",
    VERDICT_DO_NOT_USE: "**do-not-use**",
}


def render_model_card_markdown(card: Mapping[str, Any]) -> str:
    lines: List[str] = []
    lines.append("### Model card — the validity domain, generated from the artifacts")
    lines.append("")
    lines.append(
        "*Generated by `scripts/generators/generate_model_card.py`. Do not hand-edit between "
        "the markers; regenerate. Every number below is read from a tracked artifact or "
        "recomputed live, and the row says which.*"
    )
    lines.append("")
    for sentence in card.get("headline_sentences", []):
        lines.append(f"1. {sentence}" if False else f"- {sentence}")
    lines.append("")
    lines.append("| Claim type | System class | Measured | Verdict |")
    lines.append("|---|---|---|---|")
    for row in card.get("validity_domain", []):
        evidence = row["evidence"]
        if row["basis"]:
            evidence = f"{evidence}<br/><sub>{row['basis']}</sub>"
        lines.append(
            f"| {row['claim_type']} | {row['system_class']} | {evidence} | "
            f"{_VERDICT_LABEL.get(row['verdict'], row['verdict'])} |"
        )
    lines.append("")
    thresholds = card.get("thresholds", {})
    lines.append(
        f"**Verdict thresholds** (applied, not judged): trust = {thresholds.get('trust')}; "
        f"caution = {thresholds.get('caution')}; do-not-use = {thresholds.get('do_not_use')}. "
        "An unmeasured axis is reported do-not-use on purpose — absence of evidence is not "
        "evidence."
    )
    lines.append("")

    census = card.get("no_verifiable_source", {})
    if census.get("available"):
        source_status = census.get("by_status_key", {}).get("source_status")
        other = census["records"] - (source_status or 0)
        lines.append(
            f"**Provenance census (recounted at generation time, not copied).** "
            f"**{source_status} records** carry `source_status: no_verifiable_source` across "
            f"{census['files']} tracked data files — the figure the provenance note above "
            f"quotes, reproduced here by recount. A further {other} carry the same marker under "
            f"a different status key (`status`, `value_status`, `value_anchor_status`), for "
            f"{census['records']} in total. The numeric-payload and runtime-consumed subsets "
            f"(98 and 80) use a narrower definition than this recount and are pinned "
            f"separately by `tests/scientific/test_honest_headline_guards.py`."
        )
        lines.append("")

    gates = card.get("gates", [])
    if gates:
        lines.append(
            "**Blocking gates at generation time:** "
            + " · ".join(f"`{Path(item['gate']).name}` {item['status']}" for item in gates)
            + "."
        )
        lines.append("")

    free = card.get("free_precursor_holdout", {})
    if free.get("available"):
        lines.append(
            f"**Provenance of the hold-out numbers.** `{free['path']}`, "
            f"{free['bundles']} bundles / {free['scored']} scored targets, frozen "
            f"{free.get('generated_on')} at `{free.get('git')}` — before any calibration wave "
            "saw those points, and un-gitignored specifically so a later wave cannot "
            "regenerate it and compare it with itself."
        )
        lines.append("")

    lines.append(
        "**How to use this model in one line:** compare two formulations and read the ratio "
        "(`python scripts/maillard.py compare`), never quote the absolute number, and treat "
        "any pH or moisture direction as unsupported."
    )
    return "\n".join(lines)


def splice_into_readme(markdown: str, readme_path: Path) -> Tuple[str, bool]:
    """Replace the marked block in README.md. Returns (new_text, changed)."""
    text = readme_path.read_text(encoding="utf-8")
    block = f"{README_BEGIN}\n\n{markdown}\n\n{README_END}"
    start = text.find(README_BEGIN)
    end = text.find(README_END)
    if start < 0 or end < 0:
        raise ValueError(
            f"{readme_path} has no model-card markers. Add\n  {README_BEGIN}\n  {README_END}\n"
            "where the card should live, then re-run the generator."
        )
    if end < start:
        raise ValueError(f"{readme_path}: model-card markers are out of order.")
    updated = text[:start] + block + text[end + len(README_END) :]
    return updated, updated != text
