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

WHAT IT READS (2026-09-03, retirement step B5: the legacy lane's artifacts are gone)
-------------
============================================ ==========================================
the union panel                              RECOMPUTED live on the kinetic core
                                             (src.kinetic_core.scoring.score_panel, ~15 s)
core_prediction_uncertainty.json             the core's Monte-Carlo envelope (tracked)
data/**/*.json,yml                           the no_verifiable_source census, recounted
data/benchmarks/cys_ribose_140C_Hofmann1998  the sulfur branch's anchor status
scripts/ci/*.py                              the three blocking gates, actually run
============================================ ==========================================

core_directional_scores.json                 the core on the directional claims panel (B7)

NOT READ ANY MORE: the retired lane's directional-accuracy report, hold-out reports and
observability-mode comparison. They were measured on the retired screening lane.

NOTHING HERE IS A MEASUREMENT. Every number is read or recounted; the only judgement this
module adds is the trust/caution/do-not-use threshold, which is stated in the card itself.
"""

from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple

import yaml

from src import data_paths
from src.directional_reliability import (
    CAUTION_MIN_RATE,
    MIN_EVALUABLE_FOR_TRUST,
    TRUST_MIN_RATE,
    VERDICT_CAUTION,
    VERDICT_DO_NOT_USE,
    VERDICT_TRUST,
    verdict_for,
)

ENVELOPE_PATH = data_paths.CORE_PREDICTION_UNCERTAINTY
DIRECTIONAL_PATH = data_paths.CORE_DIRECTIONAL_SCORES

ROOT = data_paths.REPO_ROOT

# Quarantined 2026-09-04 (a repo-internal derivation, not a measurement); still read here as the
# record of the sulfur branch's anchor status.
SULFUR_BENCHMARK_PATH = data_paths.QUARANTINED_BENCHMARKS_DIR / "cys_ribose_140C_Hofmann1998.json"

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


def collect_directional() -> Dict[str, Any]:
    """The core on the directional claims panel (step B7), read from the tracked scorecard."""
    payload = _read_json(DIRECTIONAL_PATH)
    if payload is None or payload.get("artifact") != "core_directional_scores":
        return {"available": False, "path": _relative(DIRECTIONAL_PATH)}
    s = payload["summary"]
    ind = s["independent"]
    return {
        "available": True,
        "path": _relative(DIRECTIONAL_PATH),
        "claims": payload["panel"]["claims"],
        "headline": tuple(s["headline"]),
        "excluding_ph_aw": (ind["excluding_ph_aw"]["agree"], ind["excluding_ph_aw"]["evaluable"]),
        "ph_aw": (ind["ph_aw"]["agree"], ind["ph_aw"]["evaluable"]),
        "not_evaluable": ind["total"]["not_evaluable"],
        "mechanism_absent_misses": ind["total"]["mechanism_absent"],
        "categories": {k: (v["agree"], v["evaluable"]) for k, v in ind["by_category"].items()},
        "misses": {k: list(v["misses"]) for k, v in ind["by_category"].items() if v["misses"]},
    }


def collect_core_envelope() -> Dict[str, Any]:
    """The kinetic core's Monte-Carlo envelope (retirement step B2), read from the tracked
    artifact: n=200 takes ~40 min, so this one is not recomputed."""
    payload = _read_json(ENVELOPE_PATH)
    if payload is None:
        return {"available": False, "path": _relative(ENVELOPE_PATH)}
    s = payload.get("summary", {})
    lit = s.get("honest_literature_coverage", {}) or {}
    oos = s.get("out_of_sample_literature_coverage", {}) or {}
    return {
        "available": True,
        "path": _relative(ENVELOPE_PATH),
        "n_samples": s.get("n_samples"),
        "seed": s.get("seed"),
        "ci_level_pct": s.get("ci_level_pct"),
        "literature_hits": lit.get("hits"),
        "literature_total": lit.get("total"),
        "literature_not_evaluable": lit.get("not_evaluable"),
        "median_ci_width_log10": lit.get("median_ci_width_log10"),
        "out_of_sample_hits": oos.get("hits"),
        "out_of_sample_total": oos.get("total"),
        "unsampled_lanes": list(s.get("unsampled_lanes") or []),
    }


def collect_core_panel() -> Dict[str, Any]:
    """The kinetic core scored on the union panel (retirement step B3), recomputed live.

    ``results/validation/core_panel_scores.json`` is tracked, but the card recomputes (about
    15 s) for the same reason ``collect_benchmark_panel`` does: a number the card prints must
    be one the card can stand behind on a fresh clone.
    """
    from src.kinetic_core.scoring import score_panel

    payload = score_panel()
    s = payload["summary"]
    return {
        "available": True,
        "source": "recomputed live: src.kinetic_core.scoring.score_panel over the union panel",
        "by_lane": {
            lane: {"within_band": b["within_band"], "rows": b["rows"],
                   "median_fold_error": b["fold_summary"]["median_fold_error"]}
            for lane, b in s["by_lane"].items()
        },
        "median_fold_error": s["fold_summary"]["median_fold_error"],
        "worst_fold_error": s["fold_summary"]["worst_fold_error"],
        "out_of_sample_median_fold_error": s["out_of_sample"]["fold_summary"]["median_fold_error"],
        "total": s["panel_benchmark_count"],
        "scored": s["scored_benchmark_count"],
        "rows": s["matched_compound_count"],
        "within_band": s["within_band_count"],
        "out_of_sample_within_band": s["out_of_sample"]["within_band"],
        "out_of_sample_rows": s["out_of_sample"]["rows"],
        "strict_ready": len(s["strict_ready"]),
        "strict_ready_ids": list(s["strict_ready"]),
        "pass_band_level": payload["pass_band_level"],
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
    benchmark_dir = data_paths.BENCHMARKS_DIR
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
        "core_panel": collect_core_panel(),
        "core_envelope": collect_core_envelope(),
        "directional": collect_directional(),
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


# ---------------------------------------------------------------------------------------
# The validity-domain table: claim type x system class
# ---------------------------------------------------------------------------------------


def _headline_sentences(card: Mapping[str, Any]) -> List[str]:
    core = card["core_panel"]
    env = card["core_envelope"]
    sulfur = card["sulfur_anchor"]

    if core.get("available"):
        first = (
            f"**Absolute concentrations are unreliable.** On the union panel the kinetic core "
            f"lands {core['within_band']}/{core['rows']} rows within {core['pass_band_level']:.0f}x "
            f"(median fold error {_fold(core['median_fold_error'])}, worst "
            f"{_fold(core['worst_fold_error'])}); out of sample -- every row a core fit read "
            f"removed -- {core['out_of_sample_within_band']}/{core['out_of_sample_rows']} "
            f"(median {_fold(core['out_of_sample_median_fold_error'])}). Nothing in this "
            f"repository licenses a ppb number as a specification."
        )
    else:
        first = "**Absolute concentrations are unreliable.** The core scorecard could not be computed."
    if env.get("available"):
        first += (
            f" The core's {env['ci_level_pct']}% Monte-Carlo interval covers "
            f"{env['literature_hits']}/{env['literature_total']} evaluable literature rows "
            f"({env['literature_not_evaluable']} not evaluable: the "
            f"{', '.join(env['unsampled_lanes']) or 'no'} lane carries no sampled uncertainty), "
            f"{env['out_of_sample_hits']}/{env['out_of_sample_total']} out of sample."
        )
    directional = card.get("directional") or {}
    if directional.get("available"):
        h = directional["headline"]; ex = directional["excluding_ph_aw"]; pa = directional["ph_aw"]
        second = (
            f"**Directional and ranking claims are the product, and on the kinetic core they score "
            f"{h[0]}/{h[1]} on strictly independent literature claims** "
            f"({directional['not_evaluable']} independent claims not evaluable: refused arms, "
            f"prose-only claims, observables the core does not represent) -- "
            f"{ex[0]}/{ex[1]} once pH and water activity are set aside, and {pa[0]}/{pa[1]} on pH and "
            f"water activity themselves, {directional['mechanism_absent_misses']} of the misses being "
            f"identical predictions across an axis the lane carries no term for. A coin scores ~0.5 "
            f"on binary orderings; read the per-axis rows below, not the aggregate."
        )
    else:
        second = (
            "**Directional and ranking claims are NOT YET MEASURED on the kinetic core.** The "
            "directional scorecard is absent; every ranking claim is reported do-not-use by the "
            "card's own rule, not because it is known to be wrong."
        )
    third = (
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
                        for row in (sulfur.get("panel_verified_sulfur_values_that_are_fit_targets") or [])
                    }))
                    + " primary-source-verified sulfur row(s) are on the panel and are NOT "
                    "counted here, because a constant was selected by looking at them ("
                    + ", ".join(sorted({
                        row.split("/", 1)[0]
                        for row in (sulfur.get("panel_verified_sulfur_values_that_are_fit_targets") or [])
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
    )
    return [first, second, third]


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


_LANE_SYSTEMS = {
    "trunk": "free precursor, sugar + amine browning / furanics",
    "sulfur": "free precursor, cysteine / ribose meaty thiols",
    "acrylamide": "free precursor, asparagine + reducing sugar",
    "lipid": "protein matrix, lipid-derived aldehydes",
}


def _validity_domain(card: Mapping[str, Any]) -> List[Dict[str, str]]:
    core = card["core_panel"]
    env = card["core_envelope"]
    cells: List[DomainCell] = []

    # --- Absolute concentration, per core lane ---------------------------------------
    if core.get("available"):
        for lane, entry in sorted(core["by_lane"].items()):
            if lane == "None":
                continue
            cells.append(
                DomainCell(
                    "Absolute concentration (ppb)",
                    _LANE_SYSTEMS.get(lane, lane) + f" [{lane} lane]",
                    f"{entry['within_band']}/{entry['rows']} rows within "
                    f"{core['pass_band_level']:.0f}x, median {_fold(entry['median_fold_error'])}",
                    VERDICT_DO_NOT_USE,
                    "recomputed live on the union panel; an absolute is never trust by rule",
                )
            )
    if env.get("available"):
        cells.append(
            DomainCell(
                "Absolute concentration interval (90% CI)",
                "every lane with sampled uncertainty",
                f"{env['literature_hits']}/{env['literature_total']} evaluable literature rows "
                f"inside; {env['out_of_sample_hits']}/{env['out_of_sample_total']} out of sample; "
                f"{env['literature_not_evaluable']} not evaluable",
                VERDICT_DO_NOT_USE,
                f"{env['path']} (n={env['n_samples']}); the "
                f"{', '.join(env['unsampled_lanes']) or 'no'} lane has no sampled uncertainty",
            )
        )

    # --- Ordinal claims, per axis, from the core's directional scorecard -----------------
    directional = card.get("directional") or {}
    axis_systems = {
        "sugar_identity": "any (sugar swap, conditions held)",
        "additive_cysteine": "free precursor (cysteine present vs absent)",
        "temperature": "any (temperature moved, everything else held)",
        "time": "any (time moved, everything else held)",
        "lipid_lane": "protein matrix (lipid-derived aldehydes)",
        "matrix_identity": "protein matrix (pea vs soy)",
        "ph": "any (pH moved)",
        "moisture_aw": "any (water activity moved)",
        "ranking": "several compounds ordered in one system",
        "process_heating": "processed vs raw",
    }
    if directional.get("available"):
        for axis, system in axis_systems.items():
            agree, evaluable = directional["categories"].get(axis, (0, 0))
            misses = directional["misses"].get(axis, [])
            cells.append(
                DomainCell(
                    f"Direction / ranking on `{axis}`",
                    system,
                    f"{agree}/{evaluable} on the directional panel (independent claims)"
                    if evaluable else "no evaluable independent claim on the core",
                    verdict_for(agree, evaluable),
                    ("misses: " + ", ".join(misses)) if misses else "",
                )
            )
    else:
        cells.append(
            DomainCell(
                "Direction / ranking on any axis",
                "any",
                "not measured: the core directional scorecard is absent",
                VERDICT_DO_NOT_USE,
                "regenerate results/validation/core_directional_scores.json",
            )
        )

    # --- Panel-level statements --------------------------------------------------------
    if core.get("available"):
        ids = ", ".join(core.get("strict_ready_ids") or []) or "none"
        cells.append(
            DomainCell(
                "Any claim of benchmark-grade agreement",
                "the union panel: trust loop + hold-outs + matrix bundles",
                f"{core['strict_ready']}/{core['total']} strict-ready ({ids}); "
                f"{core['within_band']}/{core['rows']} rows within {core['pass_band_level']:.0f}x, "
                f"out-of-sample {core['out_of_sample_within_band']}/{core['out_of_sample_rows']}",
                VERDICT_DO_NOT_USE if not core["strict_ready"] else VERDICT_CAUTION,
                "recomputed live; strict-ready is the repository's own passing bar",
            )
        )

    # --- What the model is actually for -------------------------------------------------
    cells.append(
        DomainCell(
            "Which experiment to run next (value of information)",
            "any system the core envelope covers",
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
            f"(65 and 65) use a narrower definition than this recount and are pinned "
            f"separately by the headline guards under `tests/scientific/`."
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
