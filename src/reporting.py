"""
src/reporting.py

Consolidated report generation for Maillard framework simulations.
Outputs machine-readable JSON and human-readable Markdown.
"""

import datetime
import hashlib
import io
import json
import math
import platform
import re
import shlex
import subprocess
import sys
from collections import defaultdict
from contextlib import redirect_stdout
from typing import Any, Dict, List, Optional, Sequence, Mapping
from functools import lru_cache
from pathlib import Path
from typing import List, Dict, Any, Optional, Iterable

from src import data_paths
from src.input_normalization import resolve_condition_float, resolve_condition_value
from src.pipeline import FormulationResult, UncertaintyEnvelope
from src.literature_learning_loop import build_literature_learning_loop_payload
from src.family_lane_sensitivity import build_family_lane_sensitivity_payload
from src.literature_family_registry import build_family_payload_coverage_artifact, resolve_family_descriptor
from src.projection_metadata import ProjectionMetadataMap, normalize_projection_metadata_row
from src.safety import build_safety_reference_context
from src.literature_runtime import build_flavor_reference_policy_summary
from src.usability_reports import DomainWarning
from src.projection_utils import build_projection_rows, build_artifact_provenance
from src.presentation import (
    render_decision_summary_cli,
    render_flavor_axis_markdown,
    render_projection_rows_markdown,
    render_provenance_markdown,
)

SCHEMA_VERSION = "2026-03-18"


def _md_cell(value: Any) -> str:
    """Render ``value`` safely inside a GitHub-flavoured Markdown table cell.

    Several emitted labels are themselves pipe-joined (notably
    ``observable_assumption_summary``, e.g.
    ``"static_class_profile | class_level | standard_matrix_support"``), which
    used to split an 8-column row into 10 fields and corrupt the table. Pipes
    are escaped and newlines collapsed so one value stays one cell.
    """
    text = "" if value is None else str(value)
    return text.replace("\\", "\\\\").replace("|", "\\|").replace("\n", " ").strip()


def _get_uncertainty_envelope(result: FormulationResult, compound: str) -> Optional[UncertaintyEnvelope]:
    envelopes = getattr(result, "uncertainty_envelopes", {}) or {}
    return envelopes.get(compound)


def _format_compound_prediction(result: FormulationResult, compound: str, observable_ppb: float) -> str:
    envelope = _get_uncertainty_envelope(result, compound)
    if envelope is None:
        return f"{float(observable_ppb):.3g} ppb"
    return (
        f"{float(envelope.predicted_p50):.3g} ppb "
        f"[{float(envelope.predicted_p5):.3g}-{float(envelope.predicted_p95):.3g}, {int(envelope.ci_level_pct)}% CI]"
    )


def _serialize_uncertainty_envelopes(result: FormulationResult) -> Dict[str, Dict[str, Any]]:
    envelopes = getattr(result, "uncertainty_envelopes", {}) or {}
    return {
        compound: {
            "compound": envelope.compound,
            "predicted_ppb": float(envelope.predicted_ppb),
            "predicted_p5": float(envelope.predicted_p5),
            "predicted_p50": float(envelope.predicted_p50),
            "predicted_p95": float(envelope.predicted_p95),
            "ci_level_pct": int(envelope.ci_level_pct),
            "support_count": int(envelope.support_count),
            "envelope_source": str(envelope.envelope_source),
        }
        for compound, envelope in envelopes.items()
    }


def _sorted_projection_metadata(result: FormulationResult) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for raw_row in (result.projection_metadata or {}).values():
        normalized = normalize_projection_metadata_row(
            raw_row,
            compound_fallback=str(raw_row.get("compound", "unknown")),
            observable_ppb_fallback=float(raw_row.get("observable_ppb", 0.0) or 0.0),
        )
        rows.append(dict(normalized))
    rows.sort(key=lambda row: float(row.get("observable_ppb", 0.0) or 0.0), reverse=True)
    return rows


def _evidence_ladder_flags(meta: Dict[str, Any]) -> Dict[str, bool]:
    evidence_state = str(meta.get("evidence_state", "")).lower()
    source = str(meta.get("calibration_source", "")).lower()
    strength = str(meta.get("calibration_evidence_strength", "")).lower()
    fallback = str(meta.get("calibration_fallback_mode", "")).lower()
    notes_blob = " ".join(
        [
            source,
            str(meta.get("calibration_notes", "")).lower(),
            str(meta.get("retention_runtime_mode", "")).lower(),
        ]
    )

    direct_anchor = evidence_state in {"externally_benchmarked", "internally_benchmarked"} or (
        strength == "literature_anchored" and fallback == "compound_specific"
    )
    transferred_prior = evidence_state == "transferred_prior" or (
        strength in {"conditional_literature_anchored", "process_state_mismatch"}
        or fallback in {"nearest_process_state", "compound_specific_process_state"}
        or "transfer" in source
        or "carryover" in source
        or "ratio" in source
    )
    computational_refinement = any(token in notes_blob for token in ["dft", "xtb", "qm", "semiempirical", "computational", "refinement"])
    mechanistic_surrogate = (
        not direct_anchor and not transferred_prior and not computational_refinement
    ) or evidence_state == "still_missing" or strength == "heuristic" or float(meta.get("melanoidin_trapping_factor", 1.0) or 1.0) < 1.0

    return {
        "direct_anchor": bool(direct_anchor),
        "transferred_prior": bool(transferred_prior),
        "mechanistic_surrogate": bool(mechanistic_surrogate),
        "computational_refinement": bool(computational_refinement),
    }


def _support_origin(meta: Dict[str, Any]) -> str:
    process_state = str(meta.get("process_state", "unknown"))
    calibration_state = str(meta.get("calibration_process_state", "unknown"))
    fallback = str(meta.get("calibration_fallback_mode", "class_level"))
    source = str(meta.get("calibration_source", "class_fallback")).lower()
    strength = str(meta.get("calibration_evidence_strength", "heuristic")).lower()
    extrusion_state = process_state in {"aqueous_pre_extrusion_model", "extrusion_structured"}

    if extrusion_state:
        if fallback == "nearest_process_state" or (calibration_state not in {"unknown", process_state}):
            return "lower_regime_transfer"
        if strength == "heuristic" or source == "class_fallback":
            return "extrusion_extrapolation"
        return "extrusion_specific_support"

    return "standard_matrix_support"


def _build_scope_assessment_payload(result: FormulationResult) -> Dict[str, Any]:
    """Lane E (sprint 2026-05-10b): expose the calibration scope verdict.

    Returns ``{"in_scope": True, ...}`` when the formulation lives inside the
    calibration convex hull. When out of scope, the payload carries the
    human-readable reasons and the nearest calibrated (protein_type,
    process_state) so the report can render a banner and downgrade tiers.
    """
    confidence = result.confidence_metadata or {}
    raw = confidence.get("scope_assessment") or {}
    if not isinstance(raw, dict):
        raw = {}
    in_scope = bool(raw.get("in_scope", True))
    reasons = list(raw.get("reasons", []))
    nearest = dict(raw.get("nearest_calibrated", {}))
    return {
        "in_scope": in_scope,
        "reasons": [str(r) for r in reasons],
        "nearest_calibrated": {str(k): str(v) for k, v in nearest.items()},
    }


def _build_compound_evidence_ladder(result: FormulationResult, *, top_n: int = 8) -> List[Dict[str, Any]]:
    scope = _build_scope_assessment_payload(result)
    out_of_scope = not scope["in_scope"]
    ladder_rows: List[Dict[str, Any]] = []
    for meta in _sorted_projection_metadata(result)[:top_n]:
        flags = _evidence_ladder_flags(meta)
        original_strength = str(meta.get("calibration_evidence_strength", "heuristic"))
        # Lane E: visible tier label is downgraded one notch when we're outside
        # the calibration convex hull. The original tier is preserved in
        # `scope_demoted_from` so the audit trail is intact.
        if out_of_scope:
            displayed_strength = _demote_evidence_strength(original_strength)
            scope_demoted_from: Optional[str] = original_strength if displayed_strength != original_strength else None
        else:
            displayed_strength = original_strength
            scope_demoted_from = None
        ladder_rows.append(
            {
                "compound": str(meta.get("compound", "unknown")),
                "observable_ppb": float(meta.get("observable_ppb", 0.0) or 0.0),
                "direct_anchor": flags["direct_anchor"],
                "transferred_prior": flags["transferred_prior"],
                "mechanistic_surrogate": flags["mechanistic_surrogate"],
                "computational_refinement": flags["computational_refinement"],
                "evidence_state": str(meta.get("evidence_state", "still_missing")),
                "chemistry_family": str(meta.get("chemistry_family", "") or ""),
                "target_class": str(meta.get("target_class", "unknown")),
                "decision_panel_source": str(meta.get("decision_panel_source", "")),
                "support_origin": _support_origin(meta),
                "calibration_source": str(meta.get("calibration_source", "class_fallback")),
                "calibration_evidence_strength": displayed_strength,
                "calibration_fallback_mode": str(meta.get("calibration_fallback_mode", "class_level")),
                "scope_demoted_from": scope_demoted_from,
            }
        )
    return ladder_rows


_SCOPE_DEMOTION_MAP: Dict[str, str] = {
    "literature_anchored": "transferred_literature",
    "conditional_literature_anchored": "transferred_literature",
    "class_anchored": "surrogate_family",
    "directional_transferred": "surrogate_family",
    "transferred_literature": "surrogate_family",
    "bounded_calibration": "surrogate_family",
    "surrogate_family": "surrogate_only",
    "process_state_mismatch": "surrogate_only",
    "heuristic": "surrogate_only",
}


def _demote_evidence_strength(strength: str) -> str:
    return _SCOPE_DEMOTION_MAP.get(strength, "surrogate_only")



def _build_calibration_summary(result: FormulationResult, *, top_n: int = 5) -> List[Dict[str, Any]]:
    grouped: Dict[str, Dict[str, Any]] = {}
    for meta in _sorted_projection_metadata(result):
        source = str(meta.get("calibration_source", "class_fallback"))
        if not source:
            continue
        support_origin = _support_origin(meta)
        bucket = grouped.setdefault(
            f"{source}::{support_origin}",
            {
                "source": source,
                "support_origin": support_origin,
                "compounds": [],
                "observable_ppb_total": 0.0,
                "evidence_strength": str(meta.get("calibration_evidence_strength", "heuristic")),
                "fallback_mode": str(meta.get("calibration_fallback_mode", "class_level")),
            },
        )
        compound = str(meta.get("compound", "unknown"))
        if compound not in bucket["compounds"]:
            bucket["compounds"].append(compound)
        bucket["observable_ppb_total"] += float(meta.get("observable_ppb", 0.0) or 0.0)

    rows = sorted(grouped.values(), key=lambda row: float(row["observable_ppb_total"]), reverse=True)
    return rows[:top_n]


def _flatten_condition_values(value: Any) -> Iterable[str]:
    if isinstance(value, dict):
        for nested in value.values():
            yield from _flatten_condition_values(nested)
        return
    if isinstance(value, (list, tuple, set)):
        for nested in value:
            yield from _flatten_condition_values(nested)
        return
    if value is None:
        return
    yield str(value)


def _build_benchmark_neighborhood_summary(
    result: FormulationResult,
    conditions_dict: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    confidence = result.confidence_metadata or {}
    neighborhood = str(confidence.get("benchmark_neighborhood", "unknown"))
    prediction_mode = str(confidence.get("prediction_mode", "unknown"))
    process_regime = str(confidence.get("process_regime", "unknown"))
    condition_tokens = " ".join(
        token.lower() for token in _flatten_condition_values(conditions_dict or {})
    )
    hydrolysate_proxy = any(token in condition_tokens for token in ["hydrolysate", "peptide", "glutathione"])

    if hydrolysate_proxy:
        category = "hydrolysate_proxy"
        summary = "Run depends on hydrolysate/peptide-like inputs that sit outside the primary free-precursor benchmark surface."
    elif process_regime in {"extrusion_like", "extrusion_heavy"}:
        category = "extrusion_regime"
        summary = str(confidence.get("process_regime_summary", "Run is being interpreted through an extrusion-like regime without direct benchmark closure."))
    elif neighborhood in {"primary_free_precursor", "free_precursor_partial_analogy"}:
        category = "free_system_anchor"
        summary = "Run is anchored primarily by the free-system benchmark family, with varying degrees of analogy completeness."
    elif neighborhood == "matrix_intake_only":
        category = "matrix_transfer"
        summary = "Run uses matrix intake/headspace support and transferred accessibility priors rather than direct ranking benchmarks."
    else:
        category = "structural_gap"
        summary = "Run still sits beyond the strongest benchmark neighborhood and should be treated as a structural-gap extrapolation."

    return {
        "benchmark_neighborhood": neighborhood,
        "category": category,
        "prediction_mode": prediction_mode,
        "process_regime": process_regime,
        "summary": summary,
    }


def _build_missing_data_summary(result: FormulationResult) -> Dict[str, Any]:
    items: List[str] = []
    hypothesis_only: List[str] = []
    structurally_unsupported: List[str] = []

    for meta in _sorted_projection_metadata(result)[:8]:
        compound = str(meta.get("compound", "unknown"))
        source = str(meta.get("calibration_source", "class_fallback"))
        strength = str(meta.get("calibration_evidence_strength", "heuristic"))
        if source == "class_fallback" or strength == "heuristic":
            structurally_unsupported.append(compound)
            items.append(f"{compound}: still relies on class-level fallback or heuristic calibration.")
        elif str(meta.get("calibration_fallback_mode", "")) == "nearest_process_state":
            hypothesis_only.append(compound)
            items.append(f"{compound}: uses nearest-process-state transfer rather than a direct benchmark-condition anchor.")

    for compound in result.flavor_axis_summary.get("furanone_missing", []) or []:
        hypothesis_only.append(str(compound))
        items.append(f"{compound}: mechanistically expected but still unobserved in the current prediction surface.")

    benchmark_neighborhood = str((result.confidence_metadata or {}).get("benchmark_neighborhood", "unknown"))
    if benchmark_neighborhood in {"matrix_intake_only", "exploratory_matrix", "sparse_precursor_analogy"}:
        items.append(
            "Benchmark neighborhood remains extrapolative relative to the primary free-precursor validation envelope."
        )

    process_regime = str((result.confidence_metadata or {}).get("process_regime", "unknown"))
    extrusion_panel = (result.confidence_metadata or {}).get("extrusion_observable_panel", {})
    if process_regime in {"extrusion_like", "extrusion_heavy"} and not bool(extrusion_panel.get("minimum_panel_ready", False)):
        missing_categories = [
            category.replace("_", " ")
            for category in ("meaty_positive", "off_notes", "severity_markers")
            if not extrusion_panel.get(category, {}).get("present")
        ]
        items.append(
            "Extrusion observable panel remains incomplete: missing support for " + ", ".join(missing_categories) + "."
        )

    deduped_items = list(dict.fromkeys(items))[:8]
    return {
        "items": deduped_items,
        "hypothesis_only_compounds": list(dict.fromkeys(hypothesis_only)),
        "structurally_unsupported_compounds": list(dict.fromkeys(structurally_unsupported)),
    }


def _strecker_support_marker(result: FormulationResult) -> str:
    score = float(result.strecker_balance_score)
    if score >= 0.75:
        return "strong"
    if score >= 0.4:
        return "moderate"
    return "weak"


def _sulfur_trapping_summary(result: FormulationResult) -> Dict[str, Any]:
    sulfur_rows = []
    for meta in _sorted_projection_metadata(result):
        volatile_class = str(meta.get("volatile_class", "")).lower()
        compound = str(meta.get("compound", "")).lower()
        if volatile_class == "sulfur" or any(token in compound for token in ["thiol", "sulfide", "sulfur", "methional", "thiazole", "thiophene"]):
            sulfur_rows.append(meta)

    if not sulfur_rows:
        return {
            "status": "not_applicable",
            "avg_trapping_factor": 1.0,
            "summary": "No sulfur-focused observable rows were present in this run.",
        }

    avg_trapping = sum(float(row.get("melanoidin_trapping_factor", 1.0) or 1.0) for row in sulfur_rows) / len(sulfur_rows)
    if avg_trapping < 0.55:
        status = "severe"
    elif avg_trapping < 0.85:
        status = "moderate"
    else:
        status = "mild"
    return {
        "status": status,
        "avg_trapping_factor": float(avg_trapping),
        "summary": f"Average sulfur trapping factor is {avg_trapping:.2f} across {len(sulfur_rows)} sulfur-linked observable rows.",
    }


def _build_safety_reference_summary(result: FormulationResult) -> Dict[str, Any]:
    analyte = "acrylamide" if "Acrylamide" in (result.flagged_toxics or []) or float(result.safety_score) > 0.0 else "acrylamide"
    return build_safety_reference_context(analyte=analyte)


def _build_flavor_reference_policy(result: FormulationResult) -> List[Dict[str, Any]]:
    policy_rows = result.flavor_axis_summary.get("flavor_reference_policy") if getattr(result, "flavor_axis_summary", None) else None
    if isinstance(policy_rows, list) and policy_rows:
        return [dict(row) for row in policy_rows]
    return build_flavor_reference_policy_summary()


def _build_family_evidence_ladder(result: FormulationResult) -> List[Dict[str, Any]]:
    grouped: Dict[str, Dict[str, Any]] = {}
    family_lane_summary = ((result.flavor_axis_summary or {}).get("family_lane_summary", {}) or {})
    for row in _build_compound_evidence_ladder(result, top_n=24):
        chemistry_family = str(row.get("chemistry_family") or row.get("target_class") or "unknown")
        descriptor = resolve_family_descriptor(chemistry_family)
        bucket = grouped.setdefault(
            chemistry_family,
            {
                "chemistry_family": chemistry_family,
                "slr_family": str(descriptor.get("slr_family", "")).zfill(2) if descriptor else "",
                "display_name": str(descriptor.get("display_name", chemistry_family)) if descriptor else chemistry_family,
                "direct_anchor_count": 0,
                "transferred_prior_count": 0,
                "mechanistic_surrogate_count": 0,
                "computational_refinement_count": 0,
                "compounds": [],
            },
        )
        bucket["direct_anchor_count"] += int(bool(row.get("direct_anchor", False)))
        bucket["transferred_prior_count"] += int(bool(row.get("transferred_prior", False)))
        bucket["mechanistic_surrogate_count"] += int(bool(row.get("mechanistic_surrogate", False)))
        bucket["computational_refinement_count"] += int(bool(row.get("computational_refinement", False)))
        bucket["compounds"].append(str(row.get("compound", "unknown")))

    rows: List[Dict[str, Any]] = []
    for chemistry_family, bucket in grouped.items():
        posture = "structural_gap_extrapolation"
        if chemistry_family == "amino_acid_sugar_core" and bucket["direct_anchor_count"] > 0:
            posture = "core_benchmarked_chemistry"
        elif bucket["direct_anchor_count"] > 0:
            posture = "calibration_grade_family_payloads"
        elif bucket["transferred_prior_count"] > 0:
            posture = "directional_priors"
        rows.append(
            {
                **bucket,
                "evidence_posture": posture,
                "active_runtime_lane": any(str(lane.get("family_id", "")) == chemistry_family for lane in family_lane_summary.values()),
                "compounds": sorted(dict.fromkeys(bucket["compounds"])),
            }
        )
    rows.sort(key=lambda row: (row.get("slr_family", "99"), row.get("chemistry_family", "unknown")))
    return rows


def _build_family_runtime_support_summary(result: FormulationResult) -> Dict[str, Any]:
    family_lane_summary = ((result.flavor_axis_summary or {}).get("family_lane_summary", {}) or {})
    prior_bundle = ((result.flavor_axis_summary or {}).get("family_prior_bundle", {}) or {})
    evidence_ladder = _build_family_evidence_ladder(result)
    evidence_by_family = {str(row.get("chemistry_family", "unknown")): row for row in evidence_ladder}
    payload_coverage = build_family_payload_coverage_artifact()
    coverage_by_family = {str(row.get("family_id", "unknown")): row for row in payload_coverage.get("families", [])}

    rows: List[Dict[str, Any]] = []
    for slr_family, lane in sorted(family_lane_summary.items()):
        family_id = str(lane.get("family_id", "unknown"))
        evidence = evidence_by_family.get(family_id, {})
        coverage = coverage_by_family.get(family_id, {})
        prior_count = len(prior_bundle.get(family_id, []))
        evidence_posture = str(evidence.get("evidence_posture", "structural_gap_extrapolation"))
        rows.append(
            {
                "slr_family": str(slr_family),
                "family_id": family_id,
                "display_name": str(lane.get("display_name", family_id)),
                "active": bool(lane.get("active", False)),
                "strategic_posture": str(lane.get("strategic_posture", "unknown")),
                "evidence_posture": evidence_posture,
                "primary_payload_count": int(coverage.get("total_primary_payload_count", 0)),
                "supporting_payload_count": int(coverage.get("total_supporting_payload_count", 0)),
                "prior_count": int(prior_count),
                "summary": str(lane.get("summary", "")),
            }
        )
    return {
        "active_family_lane_count": sum(1 for row in rows if row.get("active")),
        "family_lanes": rows,
    }


def _build_family_specific_open_gaps(result: FormulationResult) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    support_summary = _build_family_runtime_support_summary(result)
    for row in support_summary.get("family_lanes", []):
        if not row.get("active"):
            continue
        if row.get("evidence_posture") == "structural_gap_extrapolation" or int(row.get("primary_payload_count", 0)) == 0:
            rows.append(
                {
                    "slr_family": row.get("slr_family", ""),
                    "family_id": row.get("family_id", "unknown"),
                    "display_name": row.get("display_name", "unknown"),
                    "gap_reason": "active_runtime_lane_without_direct_benchmark_or_calibration_grade_payload_closure",
                }
            )
    for compound in (result.flavor_axis_summary or {}).get("furanone_missing", []) or []:
        rows.append(
            {
                "slr_family": "01",
                "family_id": "amino_acid_sugar_core",
                "display_name": "Amino acid plus sugar core Maillard chemistry",
                "gap_reason": f"expected_marker_missing:{compound}",
            }
        )
    deduped: List[Dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for row in rows:
        key = (str(row.get("family_id", "unknown")), str(row.get("gap_reason", "unknown")))
        if key in seen:
            continue
        seen.add(key)
        deduped.append(row)
    return deduped


def _repo_root() -> Path:
    return data_paths.REPO_ROOT


def _validation_artifact(root: Path, name: str) -> Path:
    """``results/validation/<name>`` under ``root`` (a caller-supplied repo root)."""
    return root / data_paths.rel(data_paths.VALIDATION_DIR / name)


def _to_repo_relative(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path)


def _normalize_compound_key(name: str) -> str:
    if not name:
        return ""
    stripped = re.sub(r"\([^)]*\)", " ", str(name))
    return re.sub(r"[^a-z0-9]+", " ", stripped.lower()).strip()


def _configure_report_plot_style() -> None:
    from src.plot_style import configure_science_plot_style

    configure_science_plot_style()


@lru_cache(maxsize=4)
def _load_external_validation_compounds(root_str: str) -> set[str]:
    path = Path(root_str) / data_paths.rel(data_paths.EXTERNAL_VALIDATION_REPORT)
    if not path.exists():
        return set()
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return set()

    compounds: set[str] = set()
    for benchmark in payload.get("benchmarks", []) or []:
        for row in benchmark.get("compounds", []) or []:
            compound_key = _normalize_compound_key(str(row.get("compound", "")))
            if compound_key:
                compounds.add(compound_key)
    return compounds


@lru_cache(maxsize=4)
def _load_external_failing_compounds(root_str: str) -> set[str]:
    """Lane F (sprint 2026-05-10b): compounds whose mean |log10 error| exceeds
    the failing threshold on the external hold-out. Returned as normalized keys.
    """
    path = Path(root_str) / data_paths.rel(data_paths.EXTERNAL_FAILING_COMPOUNDS)
    if not path.exists():
        return set()
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return set()
    compounds: set[str] = set()
    for row in payload.get("external_failing_compounds", []) or []:
        key = _normalize_compound_key(str(row.get("compound", "")))
        if key:
            compounds.add(key)
    return compounds


def _resolve_plot_evidence_class(meta: Dict[str, Any], *, root: Path) -> str:
    compound_key = _normalize_compound_key(str(meta.get("compound", "")))
    # Lane F: external_failing wins over calibration/external/surrogate so the
    # report is honest about hexanal/nonanal-class over-prediction.
    if compound_key in _load_external_failing_compounds(str(root)):
        return "external_failing"
    flags = _evidence_ladder_flags(meta)
    if flags.get("direct_anchor", False):
        return "calibration"
    if compound_key in _load_external_validation_compounds(str(root)):
        return "external"
    return "surrogate"


def _build_compound_overlay_rows(
    result: FormulationResult,
    *,
    top_n: int = 8,
    root: Optional[Path] = None,
) -> List[Dict[str, Any]]:
    repo_root = root or _repo_root()
    confidence_rows = {
        _normalize_compound_key(str(row.get("compound", ""))): row
        for row in (result.confidence_metadata or {}).get("compound_confidence", [])
    }

    rows: List[Dict[str, Any]] = []
    for meta in _sorted_projection_metadata(result)[:top_n]:
        compound = str(meta.get("compound", "unknown"))
        envelope = _get_uncertainty_envelope(result, compound)
        observable_ppb = float(meta.get("observable_ppb", 0.0) or 0.0)
        predicted_p50 = float(envelope.predicted_p50) if envelope is not None else observable_ppb
        predicted_p5 = float(envelope.predicted_p5) if envelope is not None else predicted_p50
        predicted_p95 = float(envelope.predicted_p95) if envelope is not None else predicted_p50
        if predicted_p50 <= 0.0:
            continue
        confidence_row = confidence_rows.get(_normalize_compound_key(compound), {})
        rows.append(
            {
                "compound": compound,
                "predicted_p5": predicted_p5,
                "predicted_p50": predicted_p50,
                "predicted_p95": predicted_p95,
                "tier": str(confidence_row.get("tier", "unknown")),
                "evidence_class": _resolve_plot_evidence_class(meta, root=repo_root),
            }
        )

    if rows:
        return rows

    for compound, envelope in sorted(
        (getattr(result, "uncertainty_envelopes", {}) or {}).items(),
        key=lambda item: float(item[1].predicted_p50),
        reverse=True,
    )[:top_n]:
        rows.append(
            {
                "compound": compound,
                "predicted_p5": float(envelope.predicted_p5),
                "predicted_p50": float(envelope.predicted_p50),
                "predicted_p95": float(envelope.predicted_p95),
                "tier": "unknown",
                "evidence_class": "external_failing"
                if _normalize_compound_key(compound) in _load_external_failing_compounds(str(repo_root))
                else (
                    "external"
                    if _normalize_compound_key(compound) in _load_external_validation_compounds(str(repo_root))
                    else "surrogate"
                ),
            }
        )
    return rows


def _write_compound_confidence_overlay(
    result: FormulationResult,
    *,
    output_dir: Path,
    root: Optional[Path] = None,
) -> Optional[Path]:
    if not (getattr(result, "uncertainty_envelopes", {}) or {}):
        return None

    rows = _build_compound_overlay_rows(result, root=root)
    if not rows:
        return None

    import matplotlib

    matplotlib.use("Agg")
    _configure_report_plot_style()

    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    palette = {
        "calibration": "#2f855a",
        "external": "#d97706",
        "external_failing": "#dc2626",
        "surrogate": "#4b5563",
    }
    rows = list(reversed(rows))
    fig, ax = plt.subplots(figsize=(8.4, max(3.5, 0.65 * len(rows) + 1.6)))

    min_x = min(max(float(row["predicted_p5"]), 1.0e-6) for row in rows)
    max_x = max(max(float(row["predicted_p95"]), float(row["predicted_p50"])) for row in rows)

    for idx, row in enumerate(rows):
        color = palette.get(str(row.get("evidence_class", "surrogate")), palette["surrogate"])
        p5 = float(row["predicted_p5"])
        p50 = float(row["predicted_p50"])
        p95 = float(row["predicted_p95"])
        lower = max(p50 - p5, 0.0)
        upper = max(p95 - p50, 0.0)
        ax.errorbar(
            p50,
            idx,
            xerr=[[lower], [upper]],
            fmt="o",
            color=color,
            ecolor=color,
            capsize=3,
            elinewidth=1.8,
            markersize=6,
        )
        ax.text(max(p95, p50) * 1.08, idx, str(row.get("tier", "unknown")).upper(), va="center", fontsize=8)

    ax.set_xscale("log")
    ax.set_xlim(min_x * 0.8, max_x * 1.6)
    ax.set_yticks(list(range(len(rows))))
    ax.set_yticklabels([str(row["compound"]) for row in rows])
    ax.set_xlabel("Predicted observable ppb (P50 with 90% CI)")
    ax.set_title(f"Compound Confidence Overlay: {result.name}")
    ax.grid(True, axis="x", alpha=0.25)
    ax.legend(
        handles=[
            Patch(facecolor=color, edgecolor=color, label=label.title())
            for label, color in palette.items()
        ],
        loc="lower right",
        frameon=True,
    )
    fig.tight_layout()

    output_path = output_dir / "compound_confidence_overlay.png"
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output_path


def _format_intervention_class_label(raw_value: str) -> str:
    normalized = str(raw_value or "other").strip().lower()
    aliases = {
        "sulfur": "Sulfur / Thiols",
        "aldehyde": "Aldehydes",
        "alcohol": "Alcohols",
        "furan": "Furans",
        "pyrazine": "Pyrazines",
        "severity_markers": "Severity Markers",
        "adverse_lipid_markers": "Adverse Lipid Markers",
    }
    return aliases.get(normalized, normalized.replace("_", " ").title())


def _aggregate_intervention_class_totals(result: FormulationResult) -> Dict[str, float]:
    totals: Dict[str, float] = defaultdict(float)
    for meta in _sorted_projection_metadata(result):
        compound = str(meta.get("compound", "unknown"))
        envelope = _get_uncertainty_envelope(result, compound)
        signal_ppb = float(envelope.predicted_p50) if envelope is not None else float(meta.get("observable_ppb", 0.0) or 0.0)
        if signal_ppb <= 0.0:
            continue
        class_key = str(meta.get("volatile_class") or meta.get("target_class") or "other")
        totals[class_key] += signal_ppb
    return dict(totals)


def _normalize_precursor_key(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(name or "").lower()).strip()


def _extract_precursor_ratios(formulation_inputs: Optional[Mapping[str, Any]]) -> Dict[str, Dict[str, Any]]:
    if not isinstance(formulation_inputs, Mapping):
        return {}

    extracted: Dict[str, Dict[str, Any]] = {}
    raw_ratios = formulation_inputs.get("molar_ratios", {}) or {}
    if isinstance(raw_ratios, Mapping):
        for raw_name, raw_ratio in raw_ratios.items():
            key = _normalize_precursor_key(str(raw_name))
            if not key:
                continue
            try:
                ratio = float(raw_ratio)
            except (TypeError, ValueError):
                continue
            extracted[key] = {"precursor": str(raw_name), "ratio": ratio}

    for family_key in ("sugars", "amino_acids", "additives", "lipids"):
        values = formulation_inputs.get(family_key, ()) or ()
        if isinstance(values, str):
            values = [values]
        for raw_name in values:
            key = _normalize_precursor_key(str(raw_name))
            if not key or key in extracted:
                continue
            extracted[key] = {"precursor": str(raw_name), "ratio": 1.0}
    return extracted


def _build_intervention_compound_rows(result: FormulationResult) -> Dict[str, Dict[str, Any]]:
    rows: Dict[str, Dict[str, Any]] = {}
    for meta in _sorted_projection_metadata(result):
        compound = str(meta.get("compound", "unknown"))
        key = _normalize_compound_key(compound)
        if not key:
            continue
        envelope = _get_uncertainty_envelope(result, compound)
        predicted_p50 = float(envelope.predicted_p50) if envelope is not None else float(meta.get("observable_ppb", 0.0) or 0.0)
        rows[key] = {
            "compound": compound,
            "volatile_class": str(meta.get("volatile_class") or meta.get("target_class") or "other"),
            "predicted_p50": predicted_p50,
        }

    if rows:
        return rows

    for compound, envelope in (getattr(result, "uncertainty_envelopes", {}) or {}).items():
        key = _normalize_compound_key(compound)
        if not key:
            continue
        rows[key] = {
            "compound": str(compound),
            "volatile_class": "other",
            "predicted_p50": float(envelope.predicted_p50),
        }
    return rows


def _build_precursor_intervention_summary(
    current: FormulationResult,
    baseline: FormulationResult,
    *,
    current_inputs: Optional[Mapping[str, Any]] = None,
    baseline_inputs: Optional[Mapping[str, Any]] = None,
) -> Optional[Dict[str, Any]]:
    current_precursors = _extract_precursor_ratios(current_inputs)
    baseline_precursors = _extract_precursor_ratios(baseline_inputs)
    changed_precursors: List[Dict[str, Any]] = []

    for precursor_key in sorted(set(current_precursors) | set(baseline_precursors)):
        baseline_row = baseline_precursors.get(precursor_key, {})
        current_row = current_precursors.get(precursor_key, {})
        baseline_ratio = float(baseline_row.get("ratio", 0.0) or 0.0)
        current_ratio = float(current_row.get("ratio", 0.0) or 0.0)
        if math.isclose(baseline_ratio, current_ratio, rel_tol=1e-9, abs_tol=1e-9):
            continue
        changed_precursors.append(
            {
                "precursor": str(current_row.get("precursor") or baseline_row.get("precursor") or precursor_key),
                "baseline_ratio": baseline_ratio,
                "current_ratio": current_ratio,
                "delta_ratio": current_ratio - baseline_ratio,
            }
        )

    if not changed_precursors:
        return None

    current_compounds = _build_intervention_compound_rows(current)
    baseline_compounds = _build_intervention_compound_rows(baseline)
    compound_rows: List[Dict[str, Any]] = []
    precursor_totals: Dict[str, Dict[str, Any]] = {
        row["precursor"]: {
            **dict(row),
            "attributed_delta_ppb": 0.0,
            "compound_count": 0,
            "top_compounds": [],
        }
        for row in changed_precursors
    }

    if len(changed_precursors) == 1:
        weighted_precursors = [(changed_precursors[0], 1.0)]
        attribution_mode = "single_changed_precursor"
    else:
        total_shift = sum(abs(float(row["delta_ratio"])) for row in changed_precursors)
        if total_shift <= 0.0:
            return None
        weighted_precursors = [
            (row, abs(float(row["delta_ratio"])) / total_shift)
            for row in changed_precursors
        ]
        attribution_mode = "proportional_ratio_shift"

    for compound_key in sorted(set(current_compounds) | set(baseline_compounds)):
        baseline_row = baseline_compounds.get(compound_key, {})
        current_row = current_compounds.get(compound_key, {})
        baseline_ppb = float(baseline_row.get("predicted_p50", 0.0) or 0.0)
        current_ppb = float(current_row.get("predicted_p50", 0.0) or 0.0)
        delta_ppb = current_ppb - baseline_ppb
        if math.isclose(delta_ppb, 0.0, rel_tol=1e-9, abs_tol=1e-9):
            continue
        compound_name = str(current_row.get("compound") or baseline_row.get("compound") or compound_key)
        volatile_class = str(current_row.get("volatile_class") or baseline_row.get("volatile_class") or "other")

        for precursor_row, attribution_weight in weighted_precursors:
            attributed_delta_ppb = delta_ppb * attribution_weight
            compound_rows.append(
                {
                    "precursor": str(precursor_row["precursor"]),
                    "compound": compound_name,
                    "volatile_class": volatile_class,
                    "baseline_ppb": baseline_ppb,
                    "current_ppb": current_ppb,
                    "delta_ppb": delta_ppb,
                    "attribution_weight": attribution_weight,
                    "attributed_delta_ppb": attributed_delta_ppb,
                }
            )
            precursor_total = precursor_totals[str(precursor_row["precursor"])]
            precursor_total["attributed_delta_ppb"] += attributed_delta_ppb
            precursor_total["compound_count"] += 1
            precursor_total["top_compounds"].append(
                {
                    "compound": compound_name,
                    "attributed_delta_ppb": attributed_delta_ppb,
                }
            )

    if not compound_rows:
        return None

    compound_rows.sort(key=lambda row: abs(float(row["attributed_delta_ppb"])), reverse=True)
    precursor_summary_rows: List[Dict[str, Any]] = []
    for row in changed_precursors:
        precursor_name = str(row["precursor"])
        total_row = precursor_totals[precursor_name]
        top_compounds = sorted(
            list(total_row["top_compounds"]),
            key=lambda item: abs(float(item["attributed_delta_ppb"])),
            reverse=True,
        )[:3]
        precursor_summary_rows.append(
            {
                **dict(row),
                "attributed_delta_ppb": float(total_row["attributed_delta_ppb"]),
                "compound_count": int(total_row["compound_count"]),
                "top_compounds": top_compounds,
            }
        )

    total_compound_delta = sum(
        float(current_compounds.get(key, {}).get("predicted_p50", 0.0) or 0.0)
        - float(baseline_compounds.get(key, {}).get("predicted_p50", 0.0) or 0.0)
        for key in set(current_compounds) | set(baseline_compounds)
    )

    return {
        "attribution_mode": attribution_mode,
        "changed_precursors": changed_precursors,
        "precursor_totals": precursor_summary_rows,
        "rows": compound_rows,
        "total_compound_delta_ppb": float(total_compound_delta),
        "attributed_total_ppb": float(sum(float(row["attributed_delta_ppb"]) for row in compound_rows)),
    }


def _build_intervention_waterfall_summary(
    current: FormulationResult,
    baseline: FormulationResult,
    *,
    current_inputs: Optional[Mapping[str, Any]] = None,
    baseline_inputs: Optional[Mapping[str, Any]] = None,
    max_steps: int = 6,
) -> Optional[Dict[str, Any]]:
    current_totals = _aggregate_intervention_class_totals(current)
    baseline_totals = _aggregate_intervention_class_totals(baseline)
    step_rows: List[Dict[str, Any]] = []

    for class_key in sorted(set(current_totals) | set(baseline_totals)):
        baseline_ppb = float(baseline_totals.get(class_key, 0.0))
        current_ppb = float(current_totals.get(class_key, 0.0))
        delta_ppb = current_ppb - baseline_ppb
        if abs(delta_ppb) < 1.0e-9:
            continue
        step_rows.append(
            {
                "class_key": class_key,
                "display_name": _format_intervention_class_label(class_key),
                "baseline_ppb": baseline_ppb,
                "current_ppb": current_ppb,
                "delta_ppb": delta_ppb,
                "delta_pct": None if baseline_ppb <= 0.0 else (100.0 * delta_ppb / baseline_ppb),
            }
        )

    if not step_rows:
        return None

    step_rows.sort(key=lambda row: abs(float(row["delta_ppb"])), reverse=True)
    if len(step_rows) > max_steps:
        kept_rows = step_rows[: max_steps - 1]
        remainder_rows = step_rows[max_steps - 1 :]
        remainder_baseline = sum(float(row["baseline_ppb"]) for row in remainder_rows)
        remainder_current = sum(float(row["current_ppb"]) for row in remainder_rows)
        remainder_delta = remainder_current - remainder_baseline
        kept_rows.append(
            {
                "class_key": "other",
                "display_name": "Other",
                "baseline_ppb": remainder_baseline,
                "current_ppb": remainder_current,
                "delta_ppb": remainder_delta,
                "delta_pct": None if remainder_baseline <= 0.0 else (100.0 * remainder_delta / remainder_baseline),
            }
        )
        step_rows = kept_rows

    return {
        "baseline_name": baseline.name,
        "current_name": current.name,
        # Sum by sorted key: float addition is not associative, so an
        # order-dependent sum makes two identical runs differ by ~1 ULP.
        "baseline_total_ppb": sum(float(baseline_totals[key]) for key in sorted(baseline_totals)),
        "current_total_ppb": sum(float(current_totals[key]) for key in sorted(current_totals)),
        "steps": step_rows,
        "precursor_attribution": _build_precursor_intervention_summary(
            current,
            baseline,
            current_inputs=current_inputs,
            baseline_inputs=baseline_inputs,
        ),
    }


def _write_intervention_waterfall(
    summary: Mapping[str, Any],
    *,
    output_dir: Path,
    filename: str,
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    _configure_report_plot_style()

    import matplotlib.pyplot as plt

    steps = list(summary.get("steps", []) or [])
    baseline_total = float(summary.get("baseline_total_ppb", 0.0) or 0.0)
    current_total = float(summary.get("current_total_ppb", 0.0) or 0.0)
    labels = [str(summary.get("baseline_name", "Baseline"))] + [str(step.get("display_name", "unknown")) for step in steps] + [str(summary.get("current_name", "Current"))]

    bottoms = [0.0]
    heights = [baseline_total]
    colors = ["#6b7280"]
    annotations = [f"{baseline_total:.2f} ppb"]
    cumulative = baseline_total

    for step in steps:
        delta_ppb = float(step.get("delta_ppb", 0.0) or 0.0)
        next_total = cumulative + delta_ppb
        bottoms.append(min(cumulative, next_total))
        heights.append(abs(delta_ppb))
        colors.append("#2f855a" if delta_ppb >= 0.0 else "#c53030")
        delta_pct = step.get("delta_pct")
        pct_label = "" if delta_pct is None else f" ({float(delta_pct):+.0f}%)"
        annotations.append(f"{delta_ppb:+.2f} ppb{pct_label}")
        cumulative = next_total

    bottoms.append(0.0)
    heights.append(current_total)
    colors.append("#1f2937")
    annotations.append(f"{current_total:.2f} ppb")

    fig, ax = plt.subplots(figsize=(max(7.5, 1.15 * len(labels)), 4.6))
    max_height = max([baseline_total, current_total] + [abs(float(step.get("delta_ppb", 0.0) or 0.0)) for step in steps] + [1.0])
    for idx, (bottom, height, color, annotation) in enumerate(zip(bottoms, heights, colors, annotations)):
        ax.bar(idx, height, bottom=bottom, color=color, width=0.72)
        ax.text(idx, bottom + height + (0.04 * max_height), annotation, ha="center", va="bottom", fontsize=8)

    ax.set_xticks(list(range(len(labels))))
    ax.set_xticklabels(labels, rotation=22, ha="right")
    ax.set_ylabel("Observable ppb contribution")
    ax.set_title(f"Intervention Waterfall: {summary.get('baseline_name', 'Baseline')} -> {summary.get('current_name', 'Current')}")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()

    output_path = output_dir / filename
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output_path


def _safe_git_output(root: Path, args: List[str]) -> Optional[str]:
    try:
        completed = subprocess.run(
            ["git", *args],
            cwd=root,
            capture_output=True,
            text=True,
            check=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None
    return completed.stdout.strip()


def _build_scientific_surface(root: Path) -> Dict[str, str]:
    references = {
        "scientific_reference": root / data_paths.rel(data_paths.SCIENTIFIC_REFERENCE_MD),
        "benchmark_summary": _validation_artifact(root, "benchmark_summary.md"),
        "validated_envelope": _validation_artifact(root, "validated_envelope.md"),
        "validation_overview": _validation_artifact(root, "validation_overview.md"),
        "matrix_decision_panel": root / data_paths.rel(data_paths.MATRIX_DECISION_PANEL),
        "chemistry_family_scope_registry": root / data_paths.rel(data_paths.CHEMISTRY_FAMILY_SCOPE_REGISTRY),
        "family_ingestion_plan_registry": root / data_paths.rel(data_paths.FAMILY_INGESTION_PLAN),
        "family_identifier_contract": _validation_artifact(root, "family_identifier_contract.md"),
        "family_identifier_contract_json": _validation_artifact(root, "family_identifier_contract.json"),
        "family_strategy_policy": _validation_artifact(root, "family_strategy_policy.md"),
        "family_strategy_policy_json": _validation_artifact(root, "family_strategy_policy.json"),
        "family_payload_coverage": _validation_artifact(root, "family_payload_coverage.md"),
        "family_payload_coverage_json": _validation_artifact(root, "family_payload_coverage.json"),
        "matrix_family_coverage_registry": root / data_paths.rel(data_paths.MATRIX_FAMILY_COVERAGE_REGISTRY),
        "benchmark_intake_registry": root / data_paths.rel(data_paths.BENCHMARK_INTAKE_REGISTRY),
        "computational_priors": root / data_paths.rel(data_paths.COMPUTATIONAL_PRIORS),
        "slr_incorporation_matrix": root / data_paths.rel(data_paths.SLR_INCORPORATION_MATRIX),
        "flavor_reference_payloads": root / data_paths.rel(data_paths.FLAVOR_REFERENCE_PAYLOADS),
        "process_state_calibrations": root / data_paths.rel(data_paths.PROCESS_STATE_CALIBRATIONS),
        "retention_reference_payloads": root / data_paths.rel(data_paths.RETENTION_REFERENCE_PAYLOADS),
        "process_gap_registry": root / data_paths.rel(data_paths.PROCESS_GAP_REGISTRY),
        "safety_reference_payloads": root / data_paths.rel(data_paths.SAFETY_REFERENCE_PAYLOADS),
        "primary_benchmark_protocol": root / data_paths.rel(data_paths.PRIMARY_BENCHMARK_PROTOCOL_MD),
        "primary_benchmark_contract": root / data_paths.rel(data_paths.PPI_SPI_PRIMARY_BENCHMARK_CONTRACT),
        "literature_learning_loop": _validation_artifact(root, "literature_learning_loop.md"),
        "literature_learning_loop_json": _validation_artifact(root, "literature_learning_loop.json"),
        "literature_backlog": _validation_artifact(root, "literature_backlog.md"),
        "literature_backlog_json": _validation_artifact(root, "literature_backlog.json"),
        "deep_research_runtime_queue": _validation_artifact(root, "deep_research_runtime_queue.md"),
        "deep_research_runtime_queue_json": _validation_artifact(root, "deep_research_runtime_queue.json"),
        "literature_runtime_templates": _validation_artifact(root, "literature_runtime_templates.json"),
        "family_ingestion_plan": _validation_artifact(root, "family_ingestion_plan.md"),
        "family_ingestion_plan_json": _validation_artifact(root, "family_ingestion_plan.json"),
        "matrix_target_status": _validation_artifact(root, "matrix_target_status.md"),
        "matrix_target_status_json": _validation_artifact(root, "matrix_target_status.json"),
        "chemistry_family_scope": _validation_artifact(root, "chemistry_family_scope.md"),
        "chemistry_family_scope_json": _validation_artifact(root, "chemistry_family_scope.json"),
        "matrix_family_coverage": _validation_artifact(root, "matrix_family_coverage.md"),
        "matrix_family_coverage_json": _validation_artifact(root, "matrix_family_coverage.json"),
        "family_sensitivity": _validation_artifact(root, "family_sensitivity.md"),
        "family_sensitivity_json": _validation_artifact(root, "family_sensitivity.json"),
        "family_lane_validation": _validation_artifact(root, "family_lane_validation.md"),
        "family_lane_validation_json": _validation_artifact(root, "family_lane_validation.json"),
        "refinement_surrogate_patches": root / data_paths.rel(data_paths.REFINEMENT_SURROGATE_PATCHES),
    }
    payload: Dict[str, str] = {}
    for key, path in references.items():
        payload[key] = _to_repo_relative(path, root)
    return payload


def _build_literature_evidence_summary(root: Optional[Path] = None) -> Dict[str, Any]:
    repo_root = root or _repo_root()
    intake_path = repo_root / data_paths.rel(data_paths.BENCHMARK_INTAKE_REGISTRY)
    if not intake_path.exists():
        return {}

    with open(intake_path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)

    eligible = payload.get("eligible_references", []) or []
    structural_gaps = payload.get("structural_gaps", []) or []
    ready_refs = [entry for entry in eligible if str(entry.get("status", "")).startswith("ready_for_")]
    no_primary_data_refs = [entry for entry in eligible if not bool(entry.get("requires_primary_data", False))]

    modules: Dict[str, int] = defaultdict(int)
    matrix_families: Dict[str, int] = defaultdict(int)
    kinds: Dict[str, int] = defaultdict(int)
    for entry in ready_refs:
        kinds[str(entry.get("kind", "unknown"))] += 1
        matrix_families[str(entry.get("matrix_family", "unknown"))] += 1
        for module in entry.get("target_modules", []) or []:
            modules[str(module)] += 1

    return {
        "source": _to_repo_relative(intake_path, repo_root),
        "eligible_reference_count": len(eligible),
        "ready_reference_count": len(ready_refs),
        "closable_without_primary_data_count": len(no_primary_data_refs),
        "structural_gap_count": len(structural_gaps),
        "ready_reference_ids": [str(entry.get("id", "unknown")) for entry in ready_refs[:8]],
        "ready_by_kind": dict(sorted(kinds.items())),
        "ready_by_module": dict(sorted(modules.items())),
        "ready_by_matrix_family": dict(sorted(matrix_families.items())),
        "structural_gap_ids": [str(entry.get("id", entry.get("gap_id", "unknown"))) for entry in structural_gaps[:8]],
    }


def _build_literature_learning_loop_summary(root: Optional[Path] = None) -> Dict[str, Any]:
    repo_root = root or _repo_root()
    payload = build_literature_learning_loop_payload(repo_root)
    return dict(payload.get("summary", {}))


def _render_glossary_markdown() -> str:
    """Plain-language glossary footer for scientist-facing reports.

    Keeps a single source of truth for the evidence-tier vocabulary that appears
    inline in `report.md` and `comparison.md`, so a reader of the markdown alone
    never has to leave the file to decode a label.
    """
    return (
        "## 6. Glossary\n"
        "Plain-language meaning of the labels used above. The model is honest about *how* it knows what it claims; this section names that vocabulary.\n\n"
        "**Three different tier vocabularies appear in this report. They are not interchangeable.** Each answers a different question:\n\n"
        "**1. `tier` — how well benchmark-supported is *this run*?** Emitted on the run-level *Confidence & Support* block, on every *Compound Confidence* row, and on every *Aggregate Sensory Confidence* row. It is a band on a 0-100 confidence score, and it always travels with a `prediction_mode`:\n\n"
        "| `tier` | score band | paired `prediction_mode` | what it licenses |\n"
        "| :--- | :--- | :--- | :--- |\n"
        "| `high` | >= 85 | `benchmark_supported_quantitative` | quantitative prioritisation before wet-lab confirmation |\n"
        "| `medium` | 65-84 | `ranking_supported` | ranking and triage; verify absolute levels experimentally |\n"
        "| `low` | 45-64 | `directional_only` | direction only; absolute ppb is provisional |\n"
        "| `exploratory` | < 45 | `hypothesis_only` | hypothesis generation, not decision-grade |\n\n"
        "**2. `calibration_evidence_strength` (shown as *Evidence* in the Calibration Summary) — what kind of anchor stands behind a compound's projection?** Strongest to weakest, and only these values are emitted at compound level:\n"
        "- `literature_anchored` — a published measurement on a directly comparable system backs this compound's retention/partition treatment.\n"
        "- `conditional_literature_anchored` — literature-anchored, but only under stated conditions (pH / process-state caveats attached).\n"
        "- `class_anchored` — anchored at compound-class level (e.g. \"sulfur volatiles\"), not for this molecule specifically.\n"
        "- `directional_transferred` — a prior transferred from an adjacent matrix or process state; direction is meaningful, magnitude is not.\n"
        "- `process_state_mismatch` — an anchor exists, but for a different process state than the one you asked about; the nearest state was substituted.\n"
        "- `heuristic` — no anchor at all; a built-in class default. Ranking use only.\n\n"
        "  When the run is out of calibration scope every one of these is demoted one notch (see the scope banner below), and the pre-demotion value is preserved as `scope_demoted_from` in `report.json`.\n\n"
        "**3. `confidence_tier` — how strong is the *literature prior* behind a chemistry lane?** This one comes from the curated literature registries (`data/lit/`), not from your run, and uses a five-point scale: `high`, `medium_high`, `medium`, `medium_low`, `low`. It grades the source, not the prediction: a `high` `confidence_tier` prior can still feed an `exploratory` `tier` prediction, because your formulation may sit far from where that prior was measured.\n\n"
        "  *Name collision, stated plainly:* the campaign/comparison JSON also carries a key spelled `confidence_tier` that holds the run-level `tier` value (vocabulary 1), kept as a legacy alias. Prefer the `tier` key alongside it; only `confidence_tier` inside `data/lit/` payloads means vocabulary 3.\n\n"
        "**Scope banner.** A `⚠️ Out of calibration scope` banner at the top of the report means the matrix or process you asked about lies outside the convex hull of formulations the model has been calibrated against. The predictions are still emitted, but every compound's evidence strength is demoted one notch.\n\n"
        "**Reachability** (`reachability_status`). `chemically_reachable` — the compound is produced by an enumerated, barrier-scored pathway from your precursors. `conditionally_reachable` — reachable only under an assumption the run had to make. `merely_plausible` — no enumerated route; the number is a class-level projection.\n\n"
        "**Observable assumption** (`observable_assumption_summary`). A pipe-joined triple: retention runtime mode, calibration fallback mode, support origin — e.g. `static_class_profile | class_level | standard_matrix_support` means the volatile's matrix retention came from a static class profile, its calibration fell back to class level, and no special matrix-support route was used.\n\n"
        "**Confidence envelope.** `0.038 ppb [0.012-0.089, 90% CI]` means the p50 (median) Monte-Carlo prediction is `0.038 ppb`, with the central 90 % of samples between `0.012` and `0.089 ppb`. A compound printed without an interval had no envelope sampled. Wide envelopes make coverage cheap — read coverage and width together.\n\n"
        "**Intervention waterfall.** When two formulations are compared, the per-compound delta is decomposed into class-level (e.g. \"thiols\", \"aldehydes\") and per-precursor (e.g. \"cysteine\", \"glutathione\") contributions. Per-precursor attribution sums to the compound delta and is explicit about attribution mode.\n\n"
        "Full machine-readable trust artifacts: `results/validation/`. Per-compound 90 % envelope: `results/validation/prediction_uncertainty.md`. External hold-out: `results/validation/external_validation_report.md`.\n\n"
    )


def _render_next_experiment_markdown(
    *,
    ranking_path: Optional[Path] = None,
    prediction_path: Optional[Path] = None,
    top_n: int = 3,
    matrix_filter: Optional[Sequence[str]] = None,
) -> str:
    """Top-N value-of-information experiment recommendations.

    Sourced from `results/validation/experiment_value_ranking.json` when
    available, otherwise re-ranked on the fly via `src.experiment_value`. Emits
    a graceful stub when neither source is available — never raises.

    When ``matrix_filter`` is supplied, the table is restricted to candidates
    whose inferred ``matrix_family`` falls in the filter set. If no candidate
    matches, we fall back to the global top-N with an explicit one-line note
    so the scientist still sees actionable rows.
    """
    repo_root = _repo_root()
    default_ranking = data_paths.EXPERIMENT_VALUE_RANKING
    default_prediction = data_paths.PREDICTION_UNCERTAINTY
    ranking_path = Path(ranking_path) if ranking_path else default_ranking
    prediction_path = Path(prediction_path) if prediction_path else default_prediction

    payload: Optional[Dict[str, Any]] = None
    source_label = ""

    def _display(p: Path) -> str:
        try:
            return str(p.relative_to(repo_root))
        except ValueError:
            return str(p)

    if ranking_path.exists():
        try:
            payload = json.loads(ranking_path.read_text(encoding="utf-8"))
            source_label = _display(ranking_path)
        except (OSError, ValueError):
            payload = None

    if payload is None and prediction_path.exists():
        try:
            from src.experiment_value import (
                build_ranking_payload,
                load_prediction_payload,
                rank_experiments,
            )

            prediction_payload = load_prediction_payload(prediction_path)
            candidates = rank_experiments(prediction_payload, top_n=None)
            payload = build_ranking_payload(candidates, source_path=prediction_path)
            source_label = f"{_display(prediction_path)} (re-ranked at report time)"
        except Exception:
            payload = None

    header = "## 7. Recommended next experiment\n"
    how_to = (
        "_How to use this: run `./scripts/docker_maillard.sh next-experiment --top "
        f"{top_n}` to materialise pre-filled intake YAMLs and protocol Markdown "
        "for each row. Ingest the resulting measurement via "
        "`./scripts/docker_maillard.sh ingest --file results.csv ...`._\n\n"
    )

    if payload is None:
        return (
            header
            + "_No value-of-information ranking is available yet. Generate one with_ "
            "`./scripts/docker_maillard.sh experiment-value-ranking` _(requires "
            "`results/validation/prediction_uncertainty.json`)._\n\n"
        )

    candidates = list(payload.get("candidates", []) or [])
    matrix_note = ""
    wanted = (
        {str(m).strip().lower() for m in matrix_filter if str(m).strip()}
        if matrix_filter
        else set()
    )
    if wanted:
        scoped = [
            c for c in candidates
            if str(c.get("matrix_family", "unknown")).lower() in wanted
        ]
        if scoped:
            candidates = scoped
            matrix_note = (
                f"Scoped to matrix `{', '.join(sorted(wanted))}` (filtered from the global ranking).\n\n"
            )
        else:
            matrix_note = (
                f"_No `{', '.join(sorted(wanted))}` candidates currently above the VoI floor; "
                "showing the global top-N instead._\n\n"
            )
    candidates = candidates[: max(int(top_n), 0)]
    if not candidates:
        return (
            header
            + "_The current model has no candidate (benchmark, compound) pairs above "
            "the value-of-information floor. The trust loop is closed for now._\n\n"
        )

    lines = [
        header,
        "Ranked by value-of-information (envelope miss × CI width × ODT-anchored decision relevance). "
        f"Source: `{source_label}`.\n\n",
    ]
    if matrix_note:
        lines.append(matrix_note)
    lines.extend([
        "| Rank | VoI | Benchmark | Matrix | Compound | DoE template | Why this one |\n",
        "| ---: | ---: | --- | --- | --- | --- | --- |\n",
    ])
    for cand in candidates:
        rationale = str(cand.get("rationale", "")).replace("|", "\\|").strip()
        if len(rationale) > 220:
            rationale = rationale[:217].rstrip() + "..."
        lines.append(
            f"| {cand.get('rank')} | {float(cand.get('voi_score', 0.0)):.2f} | "
            f"`{cand.get('benchmark_id')}` | `{cand.get('matrix_family', 'unknown')}` | "
            f"{cand.get('compound')} | "
            f"`{cand.get('suggested_doe_template')}` | {rationale} |\n"
        )
    lines.append("\n")
    lines.append(how_to)
    return "".join(lines)


def _format_oav(value: Optional[float]) -> str:
    if value is None:
        return "n/a"
    if value <= 0.0:
        return "0"
    if value >= 100.0:
        return f"{value:.0f}"
    if value >= 10.0:
        return f"{value:.1f}"
    if value >= 1.0:
        return f"{value:.2f}"
    if value >= 0.01:
        return f"{value:.3f}"
    return f"{value:.2e}"


def _render_sensory_readout_markdown(
    result: FormulationResult,
    *,
    heading: str = "## 8. Sensory readout",
) -> str:
    """`## 8. Sensory readout` — per-axis OAV summary + per-compound table.

    OAV = predicted_ppb / odour_threshold_ug_per_kg using curated thresholds
    from ``data/species/{desirable_targets,off_flavour_targets}.yml``. Rows
    inherit the kernel's 90 % CI envelope (``predicted_p5/p50/p95``). Compounds
    without a curated ODT are surfaced but excluded from axis roll-ups.
    Never raises — falls back to a one-line stub on any internal error.

    Pass a different ``heading`` (e.g. ``"### Soy iso 1.0% — sensory readout"``)
    when embedding inside a comparison report so the heading hierarchy stays
    sane.
    """
    header = heading.rstrip() + "\n"
    explanation = (
        "Per-compound odour activity value (OAV = predicted ppb ÷ curated odour threshold). "
        "An axis is *above threshold* when at least one of its compounds reaches OAV ≥ 1. "
        "Compounds without a curated odour threshold are listed but excluded from axis roll-ups. "
        "Source thresholds: `data/species/desirable_targets.yml`, `data/species/off_flavour_targets.yml`.\n\n"
    )
    try:
        from src.sensory_readout import (
            AXIS_MEATY,
            AXIS_OFF_NOTE,
            AXIS_SAFETY,
            build_sensory_readout,
        )

        readout = build_sensory_readout(result)
    except Exception as exc:  # pragma: no cover - defensive
        return (
            header
            + f"_Sensory readout unavailable for this run ({exc.__class__.__name__})._\n\n"
        )

    if not readout.per_compound:
        return (
            header
            + "_FormulationResult has no predicted ppb — nothing to score._\n\n"
        )

    lines: List[str] = [header, explanation]
    lines.append(f"**Headline:** {readout.headline}.\n\n")

    lines.append("### Axis roll-up\n")
    lines.append("| Axis | Compounds (with ODT) | Above threshold | Max OAV | Top contributor |\n")
    lines.append("| --- | ---: | ---: | ---: | --- |\n")
    axis_labels = {
        AXIS_MEATY: "meaty",
        AXIS_OFF_NOTE: "off-note",
        AXIS_SAFETY: "safety",
    }
    for axis_key, label in axis_labels.items():
        rollup = readout.axes.get(axis_key)
        if rollup is None or rollup.compounds_with_odt == 0:
            lines.append(
                f"| {label} | 0 | 0 | n/a | _no compound with curated ODT in this run_ |\n"
            )
            continue
        lines.append(
            f"| {label} | {rollup.compounds_with_odt} | "
            f"{rollup.above_threshold_count} | {_format_oav(rollup.max_oav)} | "
            f"{rollup.top_contributor or '-'} |\n"
        )
    lines.append("\n")

    lines.append("### Per-compound OAV (90 % CI)\n")
    lines.append("| Compound | Axis | ODT (μg/kg) | Predicted ppb (p50) | OAV (p50) | OAV p5 | OAV p95 | ≥1? |\n")
    lines.append("| --- | --- | ---: | ---: | ---: | ---: | ---: | :---: |\n")
    sorted_rows = sorted(
        readout.per_compound,
        key=lambda r: (r.oav is None, -(r.oav or 0.0)),
    )
    for row in sorted_rows:
        odt_str = (
            f"{row.odour_threshold_ug_per_kg:.4g}"
            if row.odour_threshold_ug_per_kg is not None
            else "n/a"
        )
        ppb_basis = row.predicted_p50 if row.predicted_p50 is not None else row.predicted_ppb
        ppb_str = f"{ppb_basis:.3g}" if ppb_basis is not None else "n/a"
        flag = "✅" if row.above_threshold else ("·" if row.oav is not None else "—")
        lines.append(
            f"| {row.compound} | {row.axis} | {odt_str} | {ppb_str} | "
            f"{_format_oav(row.oav)} | {_format_oav(row.oav_p5)} | "
            f"{_format_oav(row.oav_p95)} | {flag} |\n"
        )
    lines.append("\n")

    for note in readout.notes:
        lines.append(f"_{note}_\n")
    if readout.notes:
        lines.append("\n")

    return "".join(lines)


def generate_report(
    result: FormulationResult, 
    warnings: List[DomainWarning], 
    conditions_dict: Dict[str, Any],
    output_dir: Optional[Path] = None,
    campaign_metadata: Optional[Dict[str, Any]] = None,
    baseline_result: Optional[FormulationResult] = None,
) -> Path:
    """
    Generates a consolidated report (JSON + MD) for a formulation result.
    If output_dir is None, creates a timestamped folder in 'reports/'.
    Returns the path to the report directory.
    """
    if output_dir is None:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = Path(f"reports/run_{timestamp}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    provenance = build_artifact_provenance(
        artifact_kind="single_run_report",
        output_dir=output_dir,
        inputs=conditions_dict,
        campaign_metadata=campaign_metadata,
    )
    compound_evidence_ladder = _build_compound_evidence_ladder(result)
    calibration_summary = _build_calibration_summary(result)
    missing_data_summary = _build_missing_data_summary(result)
    benchmark_neighborhood_summary = _build_benchmark_neighborhood_summary(result, conditions_dict)
    safety_reference_summary = _build_safety_reference_summary(result)
    flavor_reference_policy = _build_flavor_reference_policy(result)
    literature_evidence_summary = _build_literature_evidence_summary()
    literature_learning_loop_summary = _build_literature_learning_loop_summary()
    family_evidence_ladder = _build_family_evidence_ladder(result)
    family_runtime_support_summary = _build_family_runtime_support_summary(result)
    family_specific_open_gaps = _build_family_specific_open_gaps(result)
    family_lane_sensitivity = build_family_lane_sensitivity_payload(result.flavor_axis_summary or {})
    generated_assets: Dict[str, str] = {}
    compound_overlay_path = _write_compound_confidence_overlay(result, output_dir=output_dir, root=_repo_root())
    if compound_overlay_path is not None:
        generated_assets["compound_confidence_overlay_png"] = str(compound_overlay_path)
    intervention_waterfall = None
    intervention_waterfall_path = None
    if baseline_result is not None:
        intervention_waterfall = _build_intervention_waterfall_summary(result, baseline_result)
        if intervention_waterfall is not None:
            intervention_waterfall_path = _write_intervention_waterfall(
                intervention_waterfall,
                output_dir=output_dir,
                filename="intervention_waterfall.png",
            )
            generated_assets["intervention_waterfall_png"] = str(intervention_waterfall_path)
    
    # 1. Save JSON Report
    json_path = output_dir / "report.json"
    report_data = {
        "timestamp": datetime.datetime.now().isoformat(),
        "schema_version": SCHEMA_VERSION,
        "provenance": provenance,
        "inputs": conditions_dict,
        "results": {
            "name": result.name,
            "target_score": float(result.target_score),
            "off_flavour_risk": float(result.off_flavour_risk),
            "safety_score": float(result.safety_score),
            "lysine_budget": float(result.lysine_budget),
            "trapping_efficiency": float(result.trapping_efficiency),
            "mft_to_furfural_ratio": float(result.mft_to_furfural_ratio),
            "meaty_quality_penalty": float(result.meaty_quality_penalty),
            "strecker_balance_score": float(result.strecker_balance_score),
            "strecker_gap_penalty": float(result.strecker_gap_penalty),
            "pyrazine_propensity": float(result.pyrazine_propensity),
            "pyrazine_burden": float(result.pyrazine_burden),
            "pyrazine_penalty": float(result.pyrazine_penalty),
            "furanone_penalty": float(result.furanone_penalty),
            "flagged_toxics": result.flagged_toxics,
            "radar": {k: float(v[0]) for k, v in result.radar.items()},
            "matrix_explainability": result.matrix_explainability,
            "confidence_metadata": result.confidence_metadata,
            "scope_assessment": _build_scope_assessment_payload(result),
            "compound_evidence_ladder": compound_evidence_ladder,
            "calibration_summary": calibration_summary,
            "missing_data_summary": missing_data_summary,
            "benchmark_neighborhood_summary": benchmark_neighborhood_summary,
            "safety_reference_summary": safety_reference_summary,
            "flavor_reference_policy": flavor_reference_policy,
            "literature_evidence_summary": literature_evidence_summary,
            "literature_learning_loop_summary": literature_learning_loop_summary,
            "family_evidence_ladder": family_evidence_ladder,
            "family_runtime_support_summary": family_runtime_support_summary,
            "family_specific_open_gaps": family_specific_open_gaps,
            "family_lane_sensitivity": family_lane_sensitivity,
            "projection_metadata": dict(result.projection_metadata),
            "uncertainty_envelopes": _serialize_uncertainty_envelopes(result),
            "generated_assets": generated_assets,
            "intervention_waterfall": {
                **dict(intervention_waterfall),
                "png": str(intervention_waterfall_path),
            }
            if intervention_waterfall is not None and intervention_waterfall_path is not None
            else None,
            "flavor_axis_summary": result.flavor_axis_summary,
            "predicted_ppb": {k: float(v) for k, v in result.predicted_ppb.items()},
            "detected_targets": result.detected_targets,
            "detected_minimize": result.detected_minimize,
            "skipped_formulations": [
                dict(row) for row in (getattr(result, "skipped_formulations", []) or [])
            ],
        },
        "domain_warnings": [
            {"category": w.category, "level": w.level, "message": w.message}
            for w in warnings
        ]
    }
    
    with open(json_path, "w") as f:
        json.dump(report_data, f, indent=4, default=str)
        
    # 2. Save Markdown Report
    md_path = output_dir / "report.md"
    with open(md_path, "w") as f:
        f.write(f"# Maillard Simulation Report - {result.name}\n\n")
        f.write(f"**Date:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        scope_payload = _build_scope_assessment_payload(result)
        if not scope_payload["in_scope"]:
            nearest = scope_payload.get("nearest_calibrated", {}) or {}
            nearest_label = (
                f"{nearest.get('protein_type', 'unknown')} / {nearest.get('process_state', 'unknown')}"
                if nearest else "no calibrated neighbor"
            )
            f.write("> ⚠️ **Out of calibration scope.** This formulation lies outside the convex hull of matrices we have calibrated against, so all per-compound tier labels below are downgraded one notch (see `scope_demoted_from` in the JSON). Treat the predicted ppb values as directional only.\n>\n")
            for reason in scope_payload.get("reasons", []):
                f.write(f"> - {reason}\n")
            f.write(f">\n> Nearest calibrated regime: **{nearest_label}**.\n\n")

        f.write("## 1. Input Formulation & Conditions\n")
        f.write("| Parameter | Value |\n")
        f.write("| :--- | :--- |\n")
        for k, v in conditions_dict.items():
            f.write(f"| {_md_cell(k)} | {_md_cell(v)} |\n")
        f.write("\n")

        skipped_rows = list(getattr(result, "skipped_formulations", []) or [])
        if skipped_rows:
            f.write("> ⚠️ **Formulations dropped from this evaluation.** "
                    "The following candidates were removed before scoring because at least one "
                    "precursor could not be resolved. They are absent from every ranking below.\n>\n")
            for row in skipped_rows:
                f.write(
                    f"> - **{_md_cell(row.get('name', 'unknown'))}** — "
                    f"unresolved: `{_md_cell(row.get('unresolved_precursors', 'unknown'))}` "
                    f"({_md_cell(row.get('reason', 'unknown reason'))})\n"
                )
            f.write("\n")


        f.write("## 2. Decision Summary\n")
        f.write("```text\n")
        with io.StringIO() as buf, redirect_stdout(buf):
            render_decision_summary_cli(result, warnings)
            f.write(buf.getvalue())
        f.write("```\n\n")
        
        f.write("## 3. Detailed Results\n")
        f.write(f"- **Target Score:** {result.target_score:.2f}\n")
        f.write(f"- **Off-Flavour Risk:** {result.off_flavour_risk:.2f}\n")
        # 2026-08-27 (Wave I): unlabelled, this reads like a grade and was being read as
        # "higher is better". It is 2x the band-relative risk (src/pipeline.py s_penalty).
        f.write(
            f"- **Safety Score:** {result.safety_score:.2f} "
            "*(2\u00d7 band-relative risk, range [0, 2]; **higher is worse**; "
            ">1.0 = above the action band)*\n\n"
        )
        f.write(f"- **MFT/Furfural Ratio:** {result.mft_to_furfural_ratio:.4f}\n")
        f.write(f"- **Meaty Quality Penalty:** {result.meaty_quality_penalty:.2f}\n\n")
        f.write(f"- **Strecker Balance Score:** {result.strecker_balance_score:.2f}\n")
        f.write(f"- **Strecker Gap Penalty:** {result.strecker_gap_penalty:.2f}\n")
        f.write(f"- **Pyrazine Propensity:** {result.pyrazine_propensity:.2f}\n")
        f.write(f"- **Pyrazine Burden:** {result.pyrazine_burden:.2f}\n")
        f.write(f"- **Pyrazine Penalty:** {result.pyrazine_penalty:.2f}\n")
        f.write(f"- **Furanone Penalty:** {result.furanone_penalty:.2f}\n\n")

        calibration = {}

        if result.confidence_metadata:
            f.write("### Confidence & Support\n")
            f.write(f"- **tier:** {result.confidence_metadata.get('tier', 'unknown')}\n")
            f.write(f"- **score:** {result.confidence_metadata.get('score', 0):.1f}\n")
            f.write(f"- **benchmark_neighborhood:** {result.confidence_metadata.get('benchmark_neighborhood', 'unknown')}\n")
            f.write(f"- **decision_mode:** {result.confidence_metadata.get('decision_mode', 'directional_hypothesis')}\n")
            f.write(f"- **prediction_mode:** {result.confidence_metadata.get('prediction_mode', 'unknown')}\n")
            f.write(f"- **recommended_posture:** {result.confidence_metadata.get('recommended_posture', '')}\n")
            for factor in result.confidence_metadata.get("dominant_factors", []):
                f.write(f"- **factor:** {factor}\n")
            f.write("\n")

            calibration = result.confidence_metadata.get("calibration_diagnostics", {})
            if calibration:
                f.write("### Calibration Diagnostics\n")
                f.write(f"- **supported_envelope:** {calibration.get('supported_envelope', False)}\n")
                f.write(f"- **summary:** {calibration.get('summary', '')}\n")
                axes = calibration.get("extrapolation_axes", [])
                if axes:
                    f.write(f"- **extrapolation_axes:** {', '.join(str(axis) for axis in axes)}\n")
                f.write("\n")

            f.write("### Benchmark Neighborhood\n")
            f.write(f"- **benchmark_neighborhood:** {benchmark_neighborhood_summary.get('benchmark_neighborhood', 'unknown')}\n")
            f.write(f"- **category:** {benchmark_neighborhood_summary.get('category', 'unknown')}\n")
            f.write(f"- **prediction_mode:** {benchmark_neighborhood_summary.get('prediction_mode', 'unknown')}\n")
            f.write(f"- **summary:** {benchmark_neighborhood_summary.get('summary', '')}\n\n")

            if calibration_summary:
                f.write("### Calibration Summary\n")
                f.write("| Source | Support Origin | Evidence | Fallback | Compounds | Observable ppb |\n")
                f.write("| :--- | :--- | :--- | :--- | :--- | ---: |\n")
                for row in calibration_summary:
                    f.write(
                        f"| {_md_cell(row.get('source', 'unknown'))} | {_md_cell(row.get('support_origin', 'standard_matrix_support'))} | {_md_cell(row.get('evidence_strength', 'unknown'))} | {_md_cell(row.get('fallback_mode', 'unknown'))} | {_md_cell(', '.join(str(item) for item in row.get('compounds', [])))} | {float(row.get('observable_ppb_total', 0.0)):.2f} |\n"
                    )
                f.write("\n")

            extrusion_panel = result.confidence_metadata.get("extrusion_observable_panel", {})
            if result.confidence_metadata.get("process_regime") in {"extrusion_like", "extrusion_heavy"} and extrusion_panel:
                f.write("### Extrusion Observable Panel\n")
                f.write("| Category | Present | Missing | Ready |\n")
                f.write("| :--- | :--- | :--- | :---: |\n")
                for category in ("meaty_positive", "off_notes", "severity_markers"):
                    row = extrusion_panel.get(category, {})
                    f.write(
                        f"| {category.replace('_', ' ')} | {', '.join(row.get('present', [])) or '-'} | {', '.join(row.get('missing', [])) or '-'} | {'yes' if row.get('present') else '-'} |\n"
                    )
                f.write(f"\n- **minimum_panel_ready:** {extrusion_panel.get('minimum_panel_ready', False)}\n\n")

            extrusion_process = result.confidence_metadata.get("extrusion_process", {})
            if extrusion_process:
                total_damage = extrusion_process.get("total_damage_load", {})
                f.write("### Extrusion Process Model\n")
                f.write(f"- **model:** {extrusion_process.get('model', 'unknown')}\n")
                f.write(f"- **moisture_regime:** {extrusion_process.get('moisture_regime', 'unknown')}\n")
                f.write(f"- **jacket_temperature_celsius:** {float(extrusion_process.get('jacket_temperature_celsius', 0.0)):.1f}\n")
                f.write(f"- **effective_temperature_celsius:** {float(extrusion_process.get('effective_temperature_celsius', 0.0)):.1f}\n")
                f.write(f"- **die_exit_temperature_celsius:** {float(extrusion_process.get('die_exit_temperature_celsius', 0.0)):.1f}\n")
                f.write(f"- **sme_kj_per_kg:** {float(extrusion_process.get('sme_kj_per_kg', 0.0)):.1f}\n")
                f.write(f"- **mean_residence_time_seconds:** {float(extrusion_process.get('mean_residence_time_seconds', 0.0)):.1f}\n")
                f.write(f"- **rtd_decision:** {extrusion_process.get('rtd_assessment', {}).get('decision', 'unknown')}\n")
                f.write(f"- **furosine_mg_per_kg:** {float(total_damage.get('furosine_mg_per_kg', 0.0)):.1f}\n")
                f.write(f"- **lal_mg_per_kg:** {float(total_damage.get('lal_mg_per_kg', 0.0)):.1f}\n\n")

                transport_panel = extrusion_process.get('volatile_transport', {}).get('panel', {})
                if transport_panel:
                    f.write("| Extrusion volatile transport | Class | Shear-release | Die stripping | Combined factor |\n")
                    f.write("| :--- | :--- | ---: | ---: | ---: |\n")
                    for compound, row in transport_panel.items():
                        f.write(
                            f"| {compound} | {row.get('compound_class', 'unknown')} | {float(row.get('shear_release_factor', 0.0)):.2f} | {float(row.get('die_stripping_fraction', 0.0)):.2f} | {float(row.get('combined_headspace_factor', 0.0)):.2f} |\n"
                        )
                    f.write("\n")

            compound_rows = result.confidence_metadata.get("compound_confidence", [])
            if compound_rows:
                f.write("### Compound Confidence\n")
                f.write("| Compound | Predicted | Tier | Score | Mode | Reachability | Calibration Source | Observable Assumption |\n")
                f.write("| :--- | :--- | :---: | ---: | :--- | :--- | :--- | :--- |\n")
                for row in compound_rows:
                    compound_name = str(row.get("compound", "unknown"))
                    f.write(
                        f"| {_md_cell(compound_name)} "
                        f"| {_md_cell(_format_compound_prediction(result, compound_name, float(row.get('observable_ppb', 0.0))))} "
                        f"| {_md_cell(row.get('tier', 'unknown'))} "
                        f"| {float(row.get('score', 0.0)):.1f} "
                        f"| {_md_cell(row.get('prediction_mode', 'unknown'))} "
                        f"| {_md_cell(row.get('reachability_status', 'merely_plausible'))} "
                        f"| {_md_cell(row.get('calibration_source', 'unknown'))} "
                        f"| {_md_cell(row.get('observable_assumption_summary', 'unknown'))} |\n"
                    )
                f.write("\n")

            if compound_overlay_path is not None:
                f.write("### Compound Confidence Overlay\n")
                f.write(f"- **figure:** {compound_overlay_path.name}\n")
                f.write("- **evidence classes:** calibration = green, external = amber, surrogate = slate.\n\n")

            aggregate_rows = result.confidence_metadata.get("aggregate_confidence", {})
            if aggregate_rows:
                f.write("### Aggregate Sensory Confidence\n")
                f.write("| Tag | Sensory Score | Supporting Compounds | Tier | Mode |\n")
                f.write("| :--- | ---: | ---: | :---: | :--- |\n")
                for tag, row in aggregate_rows.items():
                    f.write(
                        f"| {tag} | {float(row.get('score', 0.0)):.2f} | {int(row.get('support_count', 0))} | {row.get('tier', 'unknown')} | {row.get('prediction_mode', 'unknown')} |\n"
                    )
                f.write("\n")

        if intervention_waterfall is not None and intervention_waterfall_path is not None:
            f.write("### Intervention Waterfall\n")
            f.write(f"- **baseline:** {intervention_waterfall.get('baseline_name', 'Baseline')}\n")
            f.write(f"- **current:** {intervention_waterfall.get('current_name', result.name)}\n")
            for step in (intervention_waterfall.get("steps", []) or [])[:3]:
                delta_ppb = float(step.get("delta_ppb", 0.0) or 0.0)
                direction = "raised" if delta_ppb >= 0.0 else "lowered"
                delta_pct = step.get("delta_pct")
                pct_text = "" if delta_pct is None else f" ({float(delta_pct):+.0f}%)"
                f.write(
                    f"- **{step.get('display_name', 'unknown')}:** {direction} {abs(delta_ppb):.2f} ppb{pct_text}\n"
                )
            f.write(f"- **figure:** {intervention_waterfall_path.name}\n\n")

        if compound_evidence_ladder:
            f.write("### Compound Evidence Ladder\n")
            f.write("| Compound | Class | Evidence State | Direct Anchor | Transferred Prior | Mechanistic Surrogate | Computational Refinement | Support Origin | Source |\n")
            f.write("| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :--- | :--- |\n")
            for row in compound_evidence_ladder:
                f.write(
                    f"| {_md_cell(row.get('compound', 'unknown'))} | {_md_cell(row.get('target_class', 'unknown'))} | {_md_cell(row.get('evidence_state', 'still_missing'))} | {'yes' if row.get('direct_anchor') else '-'} | {'yes' if row.get('transferred_prior') else '-'} | {'yes' if row.get('mechanistic_surrogate') else '-'} | {'yes' if row.get('computational_refinement') else '-'} | {_md_cell(row.get('support_origin', 'standard_matrix_support'))} | {_md_cell(row.get('calibration_source', 'unknown'))} |\n"
                )
            f.write("\n")

        f.write("### Missing Data\n")
        missing_items = missing_data_summary.get("items", [])
        if missing_items:
            for item in missing_items:
                f.write(f"- {item}\n")
        else:
            f.write("- No high-priority missing-data flags were triggered for the top reported compounds.\n")
        f.write("\n")

        safety_defaults = safety_reference_summary.get("default_entries", [])
        if safety_defaults:
            f.write("### Safety Reference Context\n")
            f.write("| Visibility | Kind | Source | Summary |\n")
            f.write("| :--- | :--- | :--- | :--- |\n")
            for row in safety_defaults:
                f.write(
                    f"| {row.get('report_visibility', 'default')} | {row.get('kind', 'unknown')} | {row.get('source_citation', 'unknown')} | {row.get('summary', '')} |\n"
                )
            extended_count = len(safety_reference_summary.get("extended_entries", []))
            if extended_count:
                f.write(f"\n- Extended safety provenance entries available in JSON: {extended_count}\n")
            # Evidence that could not be matched to this analyte is stated, not
            # silently dropped (Wave G2, 2026-08-27).
            exclusion_note = str(safety_reference_summary.get("exclusion_note", "") or "")
            if exclusion_note:
                f.write(f"\n- {exclusion_note}\n")
            f.write("\n")

        if flavor_reference_policy:
            f.write("### Flavor Reference Policy\n")
            f.write("| Compound | Pipeline Role | Benchmark Role | Source |\n")
            f.write("| :--- | :--- | :--- | :--- |\n")
            for row in flavor_reference_policy:
                f.write(
                    f"| {row.get('compound', 'unknown')} | {row.get('pipeline_role', 'reference_only')} | {row.get('benchmark_role', 'unknown')} | {row.get('source_citation', 'unknown')} |\n"
                )
            f.write("\n")

        if literature_evidence_summary:
            f.write("### Literature Evidence Summary\n")
            f.write(f"- **source:** {literature_evidence_summary.get('source', 'unknown')}\n")
            f.write(f"- **eligible_reference_count:** {literature_evidence_summary.get('eligible_reference_count', 0)}\n")
            f.write(f"- **ready_reference_count:** {literature_evidence_summary.get('ready_reference_count', 0)}\n")
            f.write(f"- **closable_without_primary_data_count:** {literature_evidence_summary.get('closable_without_primary_data_count', 0)}\n")
            f.write(f"- **structural_gap_count:** {literature_evidence_summary.get('structural_gap_count', 0)}\n")
            if literature_evidence_summary.get('ready_reference_ids'):
                f.write(f"- **ready_reference_ids:** {', '.join(str(item) for item in literature_evidence_summary.get('ready_reference_ids', []))}\n")
            if literature_evidence_summary.get('structural_gap_ids'):
                f.write(f"- **structural_gap_ids:** {', '.join(str(item) for item in literature_evidence_summary.get('structural_gap_ids', []))}\n")
            f.write("\n")

        if literature_learning_loop_summary:
            f.write("### Literature Learning Loop Summary\n")
            for key in [
                "ready_reference_count",
                "encoded_runtime_reference_count",
                "template_queue_count",
                "matrix_family_count",
                "intake_structural_gap_count",
                "process_gap_count",
            ]:
                f.write(f"- **{key}:** {literature_learning_loop_summary.get(key, 0)}\n")
            if literature_learning_loop_summary.get("matrix_prior_families"):
                f.write(f"- **matrix_prior_families:** {', '.join(str(item) for item in literature_learning_loop_summary.get('matrix_prior_families', []))}\n")
            f.write("\n")

        if family_runtime_support_summary.get("family_lanes"):
            f.write("### Family Runtime Support Summary\n")
            f.write("| SLR | Family | Active | Posture | Evidence Posture | Primary Payloads | Supporting Payloads | Priors |\n")
            f.write("| :--- | :--- | :---: | :--- | :--- | ---: | ---: | ---: |\n")
            for row in family_runtime_support_summary.get("family_lanes", []):
                f.write(
                    f"| {row.get('slr_family', '') or 'n/a'} | {row.get('display_name', 'unknown')} | {'yes' if row.get('active') else '-'} | {row.get('strategic_posture', 'unknown')} | {row.get('evidence_posture', 'structural_gap_extrapolation')} | {int(row.get('primary_payload_count', 0))} | {int(row.get('supporting_payload_count', 0))} | {int(row.get('prior_count', 0))} |\n"
                )
            f.write("\n")

        if family_evidence_ladder:
            f.write("### Family Evidence Ladder\n")
            f.write("| SLR | Family | Evidence Posture | Direct Anchors | Transferred Priors | Surrogates | Compounds |\n")
            f.write("| :--- | :--- | :--- | ---: | ---: | ---: | :--- |\n")
            for row in family_evidence_ladder:
                f.write(
                    f"| {row.get('slr_family', '') or 'n/a'} | {row.get('display_name', row.get('chemistry_family', 'unknown'))} | {row.get('evidence_posture', 'structural_gap_extrapolation')} | {int(row.get('direct_anchor_count', 0))} | {int(row.get('transferred_prior_count', 0))} | {int(row.get('mechanistic_surrogate_count', 0))} | {', '.join(str(item) for item in row.get('compounds', [])) or '-'} |\n"
                )
            f.write("\n")

        if family_specific_open_gaps:
            f.write("### Family Specific Open Gaps\n")
            for row in family_specific_open_gaps:
                f.write(
                    f"- **{row.get('display_name', row.get('family_id', 'unknown'))}:** {row.get('gap_reason', 'unknown')}\n"
                )
            f.write("\n")

        if family_lane_sensitivity.get("family_lanes"):
            f.write("### Family Lane Sensitivity\n")
            f.write("| SLR | Family Lane | Active | Target Δ | Closure Δ | Off-flavour Δ | Toggle Magnitude |\n")
            f.write("| :--- | :--- | :---: | ---: | ---: | ---: | ---: |\n")
            for row in family_lane_sensitivity.get("family_lanes", [])[:10]:
                f.write(
                    f"| {row.get('slr_family', '99')} | {row.get('display_name', 'unknown')} | {'yes' if row.get('active') else '-'} | {float(row.get('target_score_delta', 0.0)):+.2f} | {float(row.get('maillard_closure_delta', 0.0)):+.2f} | {float(row.get('off_flavour_risk_delta', 0.0)):+.2f} | {float(row.get('toggle_magnitude', 0.0)):.2f} |\n"
                )
            f.write("\n")

        if result.confidence_metadata:
            sensitivity = result.confidence_metadata.get("sensitivity_summary", {})
            if sensitivity:
                f.write("### Sensitivity Summary\n")
                f.write(f"- **mode:** {sensitivity.get('mode', 'unknown')}\n")
                f.write(f"- **evaluated_perturbations:** {sensitivity.get('evaluated_perturbations', 0)}\n")
                for item in sensitivity.get("ranking_drivers", [])[:3]:
                    f.write(
                        f"- **ranking_driver:** {item.get('input', 'unknown')} | {item.get('perturbation', 'unknown')} | "
                        f"Δdecision {float(item.get('decision_delta', 0.0)):+.2f} | Δsafety {float(item.get('safety_delta', 0.0)):+.2f}\n"
                    )
                for item in sensitivity.get("safety_drivers", [])[:3]:
                    f.write(
                        f"- **safety_driver:** {item.get('input', 'unknown')} | {item.get('perturbation', 'unknown')} | "
                        f"Δsafety {float(item.get('safety_delta', 0.0)):+.2f}\n"
                    )
                f.write("\n")
        
        if result.detected_targets:
            f.write("### Predicted Desirable Targets\n")
            f.write("| Compound | Barrier |\n")
            f.write("| :--- | :--- |\n")
            for t in result.detected_targets[:10]:
                f.write(f"| {t} | - |\n")
            f.write("\n")

        if getattr(result, "projection_metadata", None):
            projection_rows = build_projection_rows(result)
            f.write(render_projection_rows_markdown(projection_rows, heading="### Projection Calibration", variant="compact"))

            if projection_rows:
                f.write("### Trust Surface\n")
                f.write(f"- **decision_mode:** {result.confidence_metadata.get('decision_mode', 'directional_hypothesis')}\n")
                f.write(f"- **benchmark_neighborhood:** {benchmark_neighborhood_summary.get('benchmark_neighborhood', 'unknown')}\n")
                f.write(f"- **extrapolation_axes:** {', '.join(str(axis) for axis in calibration.get('extrapolation_axes', [])) if calibration else 'none'}\n")
                top_row = projection_rows[0]
                f.write(f"- **top_calibration_source:** {top_row.get('calibration_source', 'unknown')}\n")
                f.write(f"- **top_observable_assumption:** {top_row.get('observable_assumption_summary', 'unknown')}\n")
                f.write(f"- **top_reachability_status:** {top_row.get('reachability_status', 'merely_plausible')}\n\n")

        if getattr(result, "flavor_axis_summary", None):
            f.write(render_flavor_axis_markdown(result.flavor_axis_summary, heading="### Flavor Axis Diagnostics", variant="detailed"))

        f.write("## 4. Analytical Metadata\n")
        f.write("### Matrix Explainability\n")
        for k, v in result.matrix_explainability.items():
            f.write(f"- **{k}:** {v}\n")
        f.write("\n")
        f.write(render_provenance_markdown(provenance))
        f.write(_render_glossary_markdown())
        formulation_matrix = (
            str(conditions_dict.get("protein_type")).strip()
            if conditions_dict.get("protein_type")
            else None
        )
        f.write(
            _render_next_experiment_markdown(
                matrix_filter=[formulation_matrix] if formulation_matrix else None
            )
        )
        f.write(_render_sensory_readout_markdown(result))

    return output_dir

def generate_comparison_report(
    results: List[FormulationResult],
    conditions_list: List[Dict[str, Any]],
    warnings_list: Optional[List[List[DomainWarning]]] = None,
    output_dir: Optional[Path] = None,
    campaign_metadata: Optional[Dict[str, Any]] = None,
) -> Path:
    """
    Generates a side-by-side comparison report for multiple formulation results.
    """
    if output_dir is None:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = Path(f"reports/comparison_{timestamp}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    provenance = build_artifact_provenance(
        artifact_kind="formulation_comparison",
        output_dir=output_dir,
        inputs=conditions_list,
        campaign_metadata=campaign_metadata,
    )
    
    # 1. Save JSON Comparison
    json_path = output_dir / "comparison.json"
    comparison_data = {
        "timestamp": datetime.datetime.now().isoformat(),
        "schema_version": SCHEMA_VERSION,
        "provenance": provenance,
        "campaign": campaign_metadata or {},
        "runs": []
    }
    
    for res, cond in zip(results, conditions_list):
        compound_evidence_ladder = _build_compound_evidence_ladder(res)
        calibration_summary = _build_calibration_summary(res)
        missing_data_summary = _build_missing_data_summary(res)
        benchmark_neighborhood_summary = _build_benchmark_neighborhood_summary(res, cond)
        sulfur_trapping_summary = _sulfur_trapping_summary(res)
        safety_reference_summary = _build_safety_reference_summary(res)
        flavor_reference_policy = _build_flavor_reference_policy(res)
        literature_evidence_summary = _build_literature_evidence_summary()
        literature_learning_loop_summary = _build_literature_learning_loop_summary()
        family_runtime_support_summary = _build_family_runtime_support_summary(res)
        family_specific_open_gaps = _build_family_specific_open_gaps(res)
        family_lane_sensitivity = build_family_lane_sensitivity_payload(res.flavor_axis_summary or {})
        comparison_data["runs"].append({
            "name": res.name,
            "inputs": cond,
            "metrics": {
                "target_score": float(res.target_score),
                "off_flavour_risk": float(res.off_flavour_risk),
                "safety_score": float(res.safety_score),
                "lysine_budget": float(res.lysine_budget),
                "trapping_efficiency": float(res.trapping_efficiency),
                "mft_to_furfural_ratio": float(res.mft_to_furfural_ratio),
                "meaty_quality_penalty": float(res.meaty_quality_penalty),
                "strecker_balance_score": float(res.strecker_balance_score),
                "strecker_gap_penalty": float(res.strecker_gap_penalty),
                "pyrazine_propensity": float(res.pyrazine_propensity),
                "pyrazine_burden": float(res.pyrazine_burden),
                "pyrazine_penalty": float(res.pyrazine_penalty),
                "furanone_penalty": float(res.furanone_penalty),
                # `tier` is the canonical name for the run-level confidence band
                # (high/medium/low/exploratory). `confidence_tier` is kept as a
                # legacy alias of the same value — note it collides with the
                # literature-prior `confidence_tier` scale in data/lit/, which is
                # a different five-point vocabulary. See report §6.
                "tier": res.confidence_metadata.get("tier", "unknown"),
                "confidence_tier": res.confidence_metadata.get("tier", "unknown"),
                "confidence_score": float(res.confidence_metadata.get("score", 0.0)),
                "benchmark_neighborhood": res.confidence_metadata.get("benchmark_neighborhood", "unknown"),
                "prediction_mode": res.confidence_metadata.get("prediction_mode", "unknown"),
                "strecker_support_marker": _strecker_support_marker(res),
            },
            "compound_evidence_ladder": compound_evidence_ladder,
            "calibration_summary": calibration_summary,
            "missing_data_summary": missing_data_summary,
            "benchmark_neighborhood_summary": benchmark_neighborhood_summary,
            "sulfur_trapping_summary": sulfur_trapping_summary,
            "safety_reference_summary": safety_reference_summary,
            "flavor_reference_policy": flavor_reference_policy,
            "literature_evidence_summary": literature_evidence_summary,
            "literature_learning_loop_summary": literature_learning_loop_summary,
            "family_runtime_support_summary": family_runtime_support_summary,
            "family_specific_open_gaps": family_specific_open_gaps,
            "family_lane_sensitivity": family_lane_sensitivity,
            "flavor_axis_summary": res.flavor_axis_summary,
        })

    comparison_waterfall = None
    comparison_waterfall_path = None
    if len(results) == 2:
        comparison_waterfall = _build_intervention_waterfall_summary(
            results[1],
            results[0],
            current_inputs=conditions_list[1],
            baseline_inputs=conditions_list[0],
        )
        if comparison_waterfall is not None:
            comparison_waterfall_path = _write_intervention_waterfall(
                comparison_waterfall,
                output_dir=output_dir,
                filename="comparison_intervention_waterfall.png",
            )
            comparison_data["intervention_waterfall"] = {
                **dict(comparison_waterfall),
                "png": str(comparison_waterfall_path),
            }
    
    with open(json_path, "w") as f:
        json.dump(comparison_data, f, indent=4, default=str)
        
    # 2. Save Markdown Comparison
    md_path = output_dir / "comparison.md"
    with open(md_path, "w") as f:
        f.write("# Maillard Formulation Comparison Report\n\n")
        f.write(f"**Date:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## 1. Metric Overview\n")
        f.write("| Metric | " + " | ".join([res.name for res in results]) + " |\n")
        f.write("| :--- | " + " | ".join([":---:"] * len(results)) + " |\n")
        
        f.write("| **Target Score** | " + " | ".join([f"{res.target_score:.2f}" for res in results]) + " |\n")
        f.write("| **Off-Flavour Risk** | " + " | ".join([f"{res.off_flavour_risk:.2f}" for res in results]) + " |\n")
        # 2026-08-27 (Wave I): see the single-run writer above -- higher is worse.
        f.write("| **Safety Score** (2\u00d7 band risk, higher is worse) | " + " | ".join([f"{res.safety_score:.2f}" for res in results]) + " |\n")
        f.write("| **Lysine Budget** | " + " | ".join([f"{res.lysine_budget:.1f}%" for res in results]) + " |\n")
        f.write("| **Trapping Eff.** | " + " | ".join([f"{res.trapping_efficiency:.1f}%" for res in results]) + " |\n")
        f.write("| **MFT/Furfural Ratio** | " + " | ".join([f"{res.mft_to_furfural_ratio:.4f}" for res in results]) + " |\n")
        f.write("| **Meaty Quality Penalty** | " + " | ".join([f"{res.meaty_quality_penalty:.2f}" for res in results]) + " |\n")
        f.write("| **Strecker Balance** | " + " | ".join([f"{res.strecker_balance_score:.2f}" for res in results]) + " |\n")
        f.write("| **Strecker Penalty** | " + " | ".join([f"{res.strecker_gap_penalty:.2f}" for res in results]) + " |\n")
        f.write("| **Pyrazine Burden** | " + " | ".join([f"{res.pyrazine_burden:.2f}" for res in results]) + " |\n")
        f.write("| **Pyrazine Penalty** | " + " | ".join([f"{res.pyrazine_penalty:.2f}" for res in results]) + " |\n")
        f.write("| **Furanone Penalty** | " + " | ".join([f"{res.furanone_penalty:.2f}" for res in results]) + " |\n")
        f.write("| **Confidence** | " + " | ".join([f"{res.confidence_metadata.get('tier', 'unknown')} ({float(res.confidence_metadata.get('score', 0.0)):.0f})" for res in results]) + " |\n")
        f.write("| **Prediction Mode** | " + " | ".join([str(res.confidence_metadata.get('prediction_mode', 'unknown')) for res in results]) + " |\n")
        f.write("\n")
        
        f.write("## 2. Competitive Highlights\n")
        best_target = max(results, key=lambda x: x.target_score)
        best_safety = min(results, key=lambda x: x.safety_score)
        best_risk = min(results, key=lambda x: x.off_flavour_risk)
        
        f.write(f"- 🏆 **Highest Target Score:** {best_target.name} ({best_target.target_score:.2f})\n")
        f.write(f"- 🛡️ **Safest Formulation:** {best_safety.name} ({best_safety.safety_score:.2f})\n")
        f.write(f"- 🍃 **Lowest Off-Flavour Risk:** {best_risk.name} ({best_risk.off_flavour_risk:.2f})\n\n")

        if comparison_waterfall is not None and comparison_waterfall_path is not None:
            f.write("## Intervention Waterfall\n")
            f.write(f"- **baseline:** {comparison_waterfall.get('baseline_name', results[0].name)}\n")
            f.write(f"- **current:** {comparison_waterfall.get('current_name', results[1].name)}\n")
            for step in (comparison_waterfall.get("steps", []) or [])[:3]:
                delta_ppb = float(step.get("delta_ppb", 0.0) or 0.0)
                direction = "raised" if delta_ppb >= 0.0 else "lowered"
                delta_pct = step.get("delta_pct")
                pct_text = "" if delta_pct is None else f" ({float(delta_pct):+.0f}%)"
                f.write(f"- **{step.get('display_name', 'unknown')}:** {direction} {abs(delta_ppb):.2f} ppb{pct_text}\n")
            f.write(f"- **figure:** {comparison_waterfall_path.name}\n\n")

            precursor_attribution = comparison_waterfall.get("precursor_attribution") or {}
            precursor_totals = list(precursor_attribution.get("precursor_totals", []) or [])
            if precursor_totals:
                f.write("### Per-precursor intervention deltas\n")
                f.write(f"- **attribution mode:** {precursor_attribution.get('attribution_mode', 'unknown')}\n\n")
                f.write("| Precursor | Baseline Ratio | Current Ratio | Δ Ratio | Attributed Δ ppb | Top Compounds |\n")
                f.write("| :--- | ---: | ---: | ---: | ---: | :--- |\n")
                for row in precursor_totals:
                    top_compounds = ", ".join(
                        f"{item.get('compound', 'unknown')} ({float(item.get('attributed_delta_ppb', 0.0) or 0.0):+.2f})"
                        for item in (row.get("top_compounds", []) or [])[:3]
                    ) or "-"
                    f.write(
                        f"| {row.get('precursor', 'unknown')} | {float(row.get('baseline_ratio', 0.0) or 0.0):.2f} | {float(row.get('current_ratio', 0.0) or 0.0):.2f} | {float(row.get('delta_ratio', 0.0) or 0.0):+.2f} | {float(row.get('attributed_delta_ppb', 0.0) or 0.0):+.2f} | {top_compounds} |\n"
                    )
                f.write("\n")

                precursor_rows = list(precursor_attribution.get("rows", []) or [])
                if precursor_rows:
                    f.write("| Precursor | Compound | Class | Baseline ppb | Current ppb | Δ ppb | Attributed Δ ppb |\n")
                    f.write("| :--- | :--- | :--- | ---: | ---: | ---: | ---: |\n")
                    for row in precursor_rows[:12]:
                        f.write(
                            f"| {row.get('precursor', 'unknown')} | {row.get('compound', 'unknown')} | {_format_intervention_class_label(str(row.get('volatile_class', 'other')))} | {float(row.get('baseline_ppb', 0.0) or 0.0):.2f} | {float(row.get('current_ppb', 0.0) or 0.0):.2f} | {float(row.get('delta_ppb', 0.0) or 0.0):+.2f} | {float(row.get('attributed_delta_ppb', 0.0) or 0.0):+.2f} |\n"
                        )
                    f.write("\n")

        f.write("## 3. Cross-Marker Context\n")
        f.write("| Formulation | Strecker Balance | Strecker Support | Pyrazine Burden | Sulfur Trapping | Furanone Penalty | Benchmark Neighborhood | Thiamine Pathway | Thiamine Source | Expected Furanones | Missing Furanones |\n")
        f.write("| :--- | ---: | :---: | ---: | :--- | ---: | :--- | :---: | :--- | :--- | :--- |\n")
        for res in results:
            flavor_axis = res.flavor_axis_summary or {}
            trapping = _sulfur_trapping_summary(res)
            benchmark_summary = _build_benchmark_neighborhood_summary(res)
            f.write(
                f"| {res.name} | {res.strecker_balance_score:.2f} | {_strecker_support_marker(res)} | {res.pyrazine_burden:.2f} | {trapping.get('status', 'n/a')} ({float(trapping.get('avg_trapping_factor', 1.0)):.2f}) | {res.furanone_penalty:.2f} | {benchmark_summary.get('category', 'unknown')} | {str(flavor_axis.get('thiamine_pathway_active', False))} | {str(flavor_axis.get('thiamine_availability_source', 'unknown'))} | {', '.join(str(item) for item in flavor_axis.get('furanone_expected', [])) or '-'} | {', '.join(str(item) for item in flavor_axis.get('furanone_missing', [])) or '-'} |\n"
            )
        f.write("\n")

        f.write("## 4. Calibration Contrast\n")
        f.write("| Formulation | Decision Mode | Benchmark Neighborhood | Top Calibration Source | Observable Assumption | Extrapolation Axes | Missing Data Flags | Benchmark Summary |\n")
        f.write("| :--- | :--- | :--- | :--- | :--- | :--- | ---: | :--- |\n")
        for res, cond in zip(results, conditions_list):
            calibration_summary = _build_calibration_summary(res)
            top_calibration = calibration_summary[0] if calibration_summary else {
                "source": "class_fallback",
                "evidence_strength": "heuristic",
            }
            missing_data = _build_missing_data_summary(res)
            benchmark_summary = _build_benchmark_neighborhood_summary(res, cond)
            projection_rows = build_projection_rows(res)
            top_projection = projection_rows[0] if projection_rows else {}
            diagnostics = (res.confidence_metadata or {}).get("calibration_diagnostics", {})
            f.write(
                f"| {res.name} | {res.confidence_metadata.get('decision_mode', 'directional_hypothesis')} | {benchmark_summary.get('benchmark_neighborhood', 'unknown')} | {top_calibration.get('source', 'class_fallback')} | {top_projection.get('observable_assumption_summary', 'unknown')} | {', '.join(str(axis) for axis in diagnostics.get('extrapolation_axes', [])) or 'none'} | {len(missing_data.get('items', []))} | {benchmark_summary.get('summary', '')} |\n"
            )
        f.write("\n")

        f.write("## Trust Surface\n")
        f.write("| Formulation | Prediction Mode | Decision Mode | Top Reachability | Support Origin | Recommended Use |\n")
        f.write("| :--- | :--- | :--- | :--- | :--- | :--- |\n")
        for res in results:
            projection_rows = build_projection_rows(res)
            top_projection = projection_rows[0] if projection_rows else {}
            f.write(
                f"| {res.name} | {res.confidence_metadata.get('prediction_mode', 'unknown')} | {res.confidence_metadata.get('decision_mode', 'directional_hypothesis')} | {top_projection.get('reachability_status', 'merely_plausible')} | {top_projection.get('support_origin', 'standard_matrix_support')} | {res.confidence_metadata.get('recommended_posture', '')} |\n"
            )
        f.write("\n")

        for index, res in enumerate(results):
            f.write(f"### {res.name}\n")
            f.write("```text\n")
            with io.StringIO() as buf, redirect_stdout(buf):
                item_warnings = warnings_list[index] if warnings_list and index < len(warnings_list) else []
                render_decision_summary_cli(res, item_warnings)
                f.write(buf.getvalue())
            f.write("```\n\n")
        f.write(render_provenance_markdown(provenance))
        f.write(_render_glossary_markdown())
        comparison_matrices = sorted({
            str(cond.get("protein_type")).strip()
            for cond in conditions_list
            if cond.get("protein_type")
        })
        # Only narrow when every formulation in the comparison shares one matrix;
        # otherwise show the global ranking so cross-matrix studies stay unbiased.
        comparison_filter = comparison_matrices if len(comparison_matrices) == 1 else None
        f.write(
            _render_next_experiment_markdown(matrix_filter=comparison_filter)
        )
        for res in results:
            f.write(
                _render_sensory_readout_markdown(
                    res,
                    heading=f"### {res.name} — sensory readout",
                )
            )

    return output_dir


def build_campaign_leaderboard(
    results: Sequence[FormulationResult],
    conditions_list: Sequence[Any],
    shared_conditions: Optional[Mapping[str, Any]] = None,
    effective_conditions: Optional[Mapping[str, Any]] = None,
) -> List[Dict[str, Any]]:
    """Build the campaign roll-up leaderboard, ranked by target score.

    2026-08-27 (audit remediation A): the roll-up used to read
    ``conditions.get("ph")`` / ``.get("temp")`` straight off the per-run
    formulation dict and default to ``0.0``. Campaign formulations come from
    ``data/formulation_grid.yml``, which carries no process conditions at all,
    so EVERY leaderboard row printed pH 0.00 / Temp 0.0 (and Protein "free")
    while the per-run reports showed the real values. Conditions are now
    resolved through ``src.input_normalization.resolve_condition_value`` over
    the same fallback chain ``src.pipeline.evaluate_all`` uses — per-formulation
    override first, then the campaign's shared conditions, then the global
    conditions the run actually used — so alternative key spellings
    (``pH``/``temperature``) resolve too.
    """
    shared_conditions = shared_conditions or {}
    leaderboard: List[Dict[str, Any]] = []
    for result, conditions in zip(results, conditions_list):
        condition_sources = (conditions, shared_conditions, effective_conditions)
        leaderboard.append(
            {
                "name": result.name,
                "protein_type": resolve_condition_value(
                    "protein_type", condition_sources, default="free"
                ),
                # None (not 0.0) when a condition genuinely cannot be resolved:
                # an unknown pH must read as unknown, not as pH 0.
                "ph": resolve_condition_float("ph", condition_sources),
                "temp": resolve_condition_float("temp", condition_sources),
                "target_score": float(result.target_score),
                "off_flavour_risk": float(result.off_flavour_risk),
                "safety_score": float(result.safety_score),
                # See the note in generate_comparison_report: `tier` is canonical,
                # `confidence_tier` is a legacy alias of the same run-level value.
                "tier": result.confidence_metadata.get("tier", "unknown"),
                "confidence_tier": result.confidence_metadata.get("tier", "unknown"),
                "prediction_mode": result.confidence_metadata.get("prediction_mode", "unknown"),
            }
        )
    leaderboard.sort(key=lambda item: item["target_score"], reverse=True)
    return leaderboard


def generate_campaign_report(
    campaign_spec: Dict[str, Any],
    results: List[FormulationResult],
    conditions_list: List[Dict[str, Any]],
    run_artifacts: List[Dict[str, Any]],
    warnings_list: Optional[List[List[DomainWarning]]] = None,
    output_dir: Optional[Path] = None,
    effective_conditions: Optional[Mapping[str, Any]] = None,
) -> Path:
    """Roll a set of per-formulation runs up into one campaign artifact.

    `effective_conditions` is the resolved global condition set the pipeline
    actually ran with (run_campaign passes its `ReactionConditions` values). It
    is the last fallback for a condition a per-run formulation does not override,
    ahead of nothing: the leaderboard must print what was simulated.
    """
    campaign_metadata = campaign_spec.get("campaign", {})
    comparison_dir = generate_comparison_report(
        results=results,
        conditions_list=conditions_list,
        warnings_list=warnings_list,
        output_dir=output_dir,
        campaign_metadata=campaign_metadata,
    )

    leaderboard = build_campaign_leaderboard(
        results=results,
        conditions_list=conditions_list,
        shared_conditions=campaign_spec.get("shared_conditions", {}) or {},
        effective_conditions=effective_conditions,
    )

    provenance = build_artifact_provenance(
        artifact_kind="campaign_screen",
        output_dir=comparison_dir,
        inputs=campaign_spec,
        campaign_metadata=campaign_metadata,
    )

    campaign_payload = {
        "timestamp": datetime.datetime.now().isoformat(),
        "schema_version": SCHEMA_VERSION,
        "campaign": campaign_metadata,
        "shared_conditions": campaign_spec.get("shared_conditions", {}),
        "leaderboard": leaderboard,
        "run_artifacts": run_artifacts,
        "comparison_artifacts": {
            "markdown": str(comparison_dir / "comparison.md"),
            "json": str(comparison_dir / "comparison.json"),
        },
        "provenance": provenance,
    }

    json_path = comparison_dir / "campaign.json"
    with open(json_path, "w") as handle:
        json.dump(campaign_payload, handle, indent=4, default=str)

    md_path = comparison_dir / "campaign.md"
    with open(md_path, "w") as handle:
        handle.write("# Maillard Campaign Report\n\n")
        handle.write(f"**Campaign:** {campaign_metadata.get('name', 'Unnamed campaign')}\n\n")
        handle.write(f"**Objective:** {campaign_metadata.get('objective', campaign_metadata.get('description', ''))}\n\n")
        if campaign_metadata.get("audience"):
            handle.write(f"**Audience:** {campaign_metadata.get('audience')}\n\n")
        handle.write("## 1. Shared Conditions\n")
        handle.write("| Parameter | Value |\n")
        handle.write("| :--- | :--- |\n")
        for key, value in campaign_spec.get("shared_conditions", {}).items():
            handle.write(f"| {key} | {value} |\n")
        handle.write("\n")
        handle.write("## 2. Leaderboard\n")
        handle.write("| Formulation | Protein | pH | Temp C | Target | Off-flavour | Safety | Tier | Mode |\n")
        handle.write("| :--- | :--- | ---: | ---: | ---: | ---: | ---: | :---: | :--- |\n")
        for row in leaderboard:
            ph_cell = f"{row['ph']:.2f}" if row["ph"] is not None else "n/a"
            temp_cell = f"{row['temp']:.1f}" if row["temp"] is not None else "n/a"
            handle.write(
                f"| {_md_cell(row['name'])} | {_md_cell(row['protein_type'])} | {ph_cell} | {temp_cell} | {row['target_score']:.2f} | {row['off_flavour_risk']:.2f} | {row['safety_score']:.2f} | {_md_cell(row['tier'])} | {_md_cell(row['prediction_mode'])} |\n"
            )
        handle.write("\n")
        handle.write("## 3. Generated Artifacts\n")
        handle.write(f"- comparison markdown: {comparison_dir / 'comparison.md'}\n")
        handle.write(f"- comparison json: {comparison_dir / 'comparison.json'}\n")
        handle.write(f"- campaign json: {json_path}\n")
        for artifact in run_artifacts:
            handle.write(
                f"- run artifact: {artifact.get('name', 'unknown')} -> {artifact.get('directory', '')}\n"
            )
        handle.write("\n")
        handle.write(render_provenance_markdown(provenance))

    return comparison_dir
