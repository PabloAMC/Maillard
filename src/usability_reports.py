from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional, Tuple, TYPE_CHECKING

from src.benchmark_validation import summarize_benchmarks
from src.conditions import ReactionConditions
from src.projection_utils import build_projection_rows
from src.pipeline import FormulationResult, MaillardPipeline


@dataclass(frozen=True)
class ValidatedEnvelopeReport:
    target_tag: str
    strict_ready_benchmarks: List[str]
    supported_benchmarks: int
    total_benchmarks: int
    matrix_only_benchmarks: List[str]
    warnings: List[str]
    next_priorities: List[str]


def build_validated_envelope_report(target_tag: str = "meaty") -> ValidatedEnvelopeReport:
    summaries = summarize_benchmarks(target_tag=target_tag)
    strict_ready = [summary.benchmark_id for summary in summaries if summary.strict_ready]
    matrix_only = [summary.benchmark_id for summary in summaries if summary.execution_path == "matrix_only"]
    warnings = []
    if matrix_only:
        warnings.append(
            "Matrix benchmarks are executable intake/headspace checks, but remain outside the strict gate and target snapshots."
        )
    warnings.append(
        "Benchmark-facing concentrations still use the FAST observable projection; Cantera remains diagnostic-reference-only."
    )
    warnings.append(
        "Peptide-bound and intact-protein accessibility remain outside the validated precursor envelope."
    )
    warnings.append(
        "The validated plant-matrix scope is currently limited to pea/soy matrix-only systems and not yet a broad release gate."
    )
    next_priorities = [
        "Expose matrix-state and projection explainability in user-facing artifacts.",
        "Promote domain-of-applicability warnings into routine CLI/report outputs.",
        "Replace bulk matrix retention with broader compound-aware observability across plant-matrix systems.",
    ]
    return ValidatedEnvelopeReport(
        target_tag=target_tag,
        strict_ready_benchmarks=strict_ready,
        supported_benchmarks=sum(1 for summary in summaries if summary.supported),
        total_benchmarks=len(summaries),
        matrix_only_benchmarks=matrix_only,
        warnings=warnings,
        next_priorities=next_priorities,
    )


def build_formulation_explainability_payload(
    formulation: Dict[str, object],
    result: FormulationResult,
    *,
    target_tag: str,
    minimize_tag: str,
) -> Dict[str, object]:
    all_rows = build_projection_rows(result)
    top_rows = all_rows[:8] 

    return {
        "formulation_name": formulation.get("name", "unknown"),
        "target_tag": target_tag,
        "minimize_tag": minimize_tag,
        "protein_type": formulation.get("protein_type", "free"),
        "matrix_explainability": result.matrix_explainability,
        "scores": {
            "target_score": result.target_score,
            "off_flavour_risk": result.off_flavour_risk,
            "safety_score": result.safety_score,
            "texture_risk": result.texture_risk,
            "lysine_budget": result.lysine_budget,
            "trapping_efficiency": result.trapping_efficiency,
            "strecker_balance_score": result.strecker_balance_score,
            "pyrazine_burden": result.pyrazine_burden,
        },
        "flavor_axis_summary": result.flavor_axis_summary,
        "top_projection_rows": top_rows,
        "detected_targets": result.detected_targets,
        "detected_minimize": result.detected_minimize,
    }


@dataclass(frozen=True)
class DomainWarning:
    level: str  # "WARNING" or "CAUTION"
    category: str
    message: str


@dataclass(frozen=True)
class ConfidenceAssessment:
    tier: str
    score: float
    benchmark_neighborhood: str
    support_ratio: float
    avg_uncertainty_kcal: float
    dominant_factors: List[str]
    recommended_posture: str


def _confidence_tier_from_score(score: float) -> str:
    if score >= 85.0:
        return "high"
    if score >= 65.0:
        return "medium"
    if score >= 45.0:
        return "low"
    return "exploratory"


def _prediction_mode_from_tier(tier: str) -> str:
    if tier == "high":
        return "benchmark_supported_quantitative"
    if tier == "medium":
        return "ranking_supported"
    if tier == "low":
        return "directional_only"
    return "hypothesis_only"


def _decision_mode_from_prediction_mode(prediction_mode: str) -> str:
    if str(prediction_mode) == "benchmark_supported_quantitative":
        return "quantitative_recommendation"
    return "directional_hypothesis"


def _clamp_confidence_score(score: float) -> float:
    return max(5.0, min(100.0, float(score)))


def _recommended_posture_from_tier(tier: str) -> str:
    if tier == "high":
        return "Suitable for quantitative prioritization before wet-lab confirmation."
    if tier == "medium":
        return "Reliable for ranking and formulation triage; verify absolute levels experimentally."
    if tier == "low":
        return "Use directionally only; absolute concentrations should be treated as provisional."
    return "Use for hypothesis generation, not decision-grade prioritization."


def _classify_process_regime(
    *,
    protein_type: str,
    temp_c: Optional[float],
    aw: Optional[float],
) -> Dict[str, object]:
    protein = str(protein_type or "free").strip().lower()
    temperature = None if temp_c is None else float(temp_c)
    water_activity = None if aw is None else float(aw)

    if protein == "free":
        return {
            "process_regime": "free_aqueous",
            "process_neighborhood": "in_domain",
            "prediction_mode_override": None,
            "summary": "Free-system conditions stay inside the baseline aqueous process neighborhood.",
        }

    if temperature is None or water_activity is None:
        return {
            "process_regime": "matrix_hydrated",
            "process_neighborhood": "unspecified",
            "prediction_mode_override": None,
            "summary": "Process severity is not fully specified, so the run defaults to the hydrated matrix neighborhood.",
        }

    if temperature >= 160.0 and water_activity <= 0.45:
        return {
            "process_regime": "extrusion_heavy",
            "process_neighborhood": "out_of_domain",
            "prediction_mode_override": "hypothesis_only",
            "summary": "Temperature and water activity fall into an extrusion-heavy regime that currently sits outside the validated matrix benchmark neighborhood.",
        }

    if temperature >= 140.0 and water_activity <= 0.65:
        return {
            "process_regime": "extrusion_like",
            "process_neighborhood": "near_domain",
            "prediction_mode_override": "directional_only",
            "summary": "Process severity is extrusion-like, so support is transferred from lower-severity matrix benchmarks rather than directly benchmarked here.",
        }

    return {
        "process_regime": "matrix_hydrated",
        "process_neighborhood": "in_domain",
        "prediction_mode_override": None,
        "summary": "Run remains in the hydrated matrix neighborhood currently supported by pea/soy accessibility and observability assumptions.",
    }


def _apply_prediction_mode_override(
    assessment: ConfidenceAssessment,
    *,
    prediction_mode_override: Optional[str],
    dominant_reason: Optional[str] = None,
) -> ConfidenceAssessment:
    if prediction_mode_override not in {"directional_only", "hypothesis_only"}:
        return assessment

    current_mode = _prediction_mode_from_tier(assessment.tier)
    if prediction_mode_override == current_mode:
        return assessment

    adjusted_score = min(float(assessment.score), 58.0 if prediction_mode_override == "directional_only" else 40.0)
    adjusted_tier = _confidence_tier_from_score(adjusted_score)
    dominant_factors = list(assessment.dominant_factors)
    if dominant_reason and dominant_reason not in dominant_factors:
        dominant_factors.insert(0, dominant_reason)

    return ConfidenceAssessment(
        tier=adjusted_tier,
        score=adjusted_score,
        benchmark_neighborhood=assessment.benchmark_neighborhood,
        support_ratio=assessment.support_ratio,
        avg_uncertainty_kcal=assessment.avg_uncertainty_kcal,
        dominant_factors=dominant_factors[:3],
        recommended_posture=_recommended_posture_from_tier(adjusted_tier),
    )


def _build_extrusion_observable_panel(result: "FormulationResult") -> Dict[str, object]:
    rows = build_projection_rows(result)
    signal_map = {
        str(row.get("compound", "")).strip().lower(): float(row.get("observable_ppb", 0.0) or 0.0)
        for row in rows
        if str(row.get("compound", "")).strip()
    }
    panel_targets = {
        "meaty_positive": [
            "2-Methyl-3-furanthiol (MFT)",
            "2-Furfurylthiol (FFT)",
            "Methional",
            "2,5-Dimethylpyrazine",
        ],
        "off_notes": [
            "Hexanal",
            "2-Pentylfuran",
            "Nonanal",
            "1-Hexanol",
        ],
        "severity_markers": [
            "Furfural",
            "5-Hydroxymethylfurfural (HMF)",
            "2-Acetylfuran",
        ],
    }

    panel: Dict[str, object] = {}
    minimum_panel_ready = True
    for category, compounds in panel_targets.items():
        present = [compound for compound in compounds if signal_map.get(compound.strip().lower(), 0.0) > 0.0]
        missing = [compound for compound in compounds if compound not in present]
        minimum_panel_ready = minimum_panel_ready and bool(present)
        panel[category] = {
            "required_count": len(compounds),
            "present_count": len(present),
            "present": present,
            "missing": missing,
        }

    panel["minimum_panel_ready"] = minimum_panel_ready
    return panel


def _decision_proxy(result: "FormulationResult") -> float:
    return float(result.target_score) - float(result.safety_score) - float(result.off_flavour_risk) - 0.01 * float(result.texture_risk)


def _iter_precursor_inputs(formulation: Dict[str, object]) -> List[Tuple[str, str]]:
    ordered = []
    for key in ("sugars", "amino_acids", "additives", "lipids"):
        for item in formulation.get(key, []) or []:
            if isinstance(item, str) and item.strip():
                ordered.append((key, item.strip()))
    return ordered


def _build_conditions_for_variant(
    formulation: Dict[str, object],
    baseline_conditions: Optional[ReactionConditions],
) -> ReactionConditions:
    baseline = baseline_conditions or ReactionConditions()
    return ReactionConditions(
        pH=float(formulation.get("ph", getattr(baseline, "pH", 6.0))),
        temperature_celsius=float(formulation.get("temp", getattr(baseline, "temperature_celsius", 120.0))),
        water_activity=float(formulation.get("aw", getattr(baseline, "water_activity", 0.8))),
        fat_fraction=float(getattr(baseline, "fat_fraction", 0.0)),
        protein_fraction=float(getattr(baseline, "protein_fraction", 1.0)),
        dielectric_constant=float(getattr(baseline, "dielectric_constant", 78.4)),
        solvent_name=str(getattr(baseline, "solvent_name", "water")),
        matrix_fiber=float(getattr(baseline, "matrix_fiber", 0.0)),
        metal_catalyst=getattr(baseline, "metal_catalyst", None),
        protein_type=str(formulation.get("protein_type", getattr(baseline, "protein_type", "free"))),
        sme_kj_per_kg=float(formulation.get("sme_kj_per_kg", getattr(baseline, "sme_kj_per_kg", 0.0))),
        moisture_regime=formulation.get("moisture_regime", getattr(baseline, "moisture_regime", None)),
        sterilization_temperature_celsius=formulation.get(
            "sterilization_temperature_celsius",
            getattr(baseline, "sterilization_temperature_celsius", None),
        ),
        sterilization_time_minutes=float(
            formulation.get(
                "sterilization_time_minutes",
                getattr(baseline, "sterilization_time_minutes", 0.0),
            )
        ),
        barrel_zone_temperatures=formulation.get(
            "barrel_zone_temperatures",
            getattr(baseline, "barrel_zone_temperatures", None),
        ),
        barrel_zone_time_fractions=formulation.get(
            "barrel_zone_time_fractions",
            getattr(baseline, "barrel_zone_time_fractions", None),
        ),
    )


def _build_calibration_diagnostics(
    assessment: ConfidenceAssessment,
    result: "FormulationResult",
    warnings: List[DomainWarning],
) -> Dict[str, object]:
    axes: List[str] = []
    categories = sorted({warning.category for warning in warnings})
    if assessment.benchmark_neighborhood != "primary_free_precursor":
        axes.append("benchmark_neighborhood")
    if categories:
        axes.extend(category.lower() for category in categories)
    if float(assessment.avg_uncertainty_kcal) > 5.0:
        axes.append("pathway_uncertainty")
    if float(getattr(result, "bottleneck_severity", 0.0) or 0.0) > 0.6:
        axes.append("projection_severity")
    axes = sorted(set(axes))

    if not axes:
        summary = "Run stays inside the most trusted benchmark/process envelope."
    else:
        summary = "Recommendation extrapolates beyond the strongest support on: " + ", ".join(axes) + "."

    return {
        "supported_envelope": len(axes) == 0,
        "prediction_mode": _prediction_mode_from_tier(assessment.tier),
        "extrapolation_axes": axes,
        "warning_categories": categories,
        "summary": summary,
    }


def _build_compound_confidence_rows(
    result: "FormulationResult",
    assessment: ConfidenceAssessment,
    *,
    top_n: int = 5,
) -> List[Dict[str, object]]:
    derived_projection_rows = {
        str(row.get("compound", "")).strip().lower(): row
        for row in build_projection_rows(result)
        if str(row.get("compound", "")).strip()
    }
    target_map = {
        str(target.get("name", "")).strip().lower(): target
        for target in getattr(result, "targets", [])
        if isinstance(target, dict) and target.get("name")
    }
    ranked_metadata = sorted(
        result.projection_metadata.items(),
        key=lambda item: float(item[1].get("observable_ppb", 0.0)),
        reverse=True,
    )
    rows: List[Dict[str, object]] = []
    for canon, meta in ranked_metadata[:top_n]:
        compound = str(meta.get("compound", canon))
        derived_meta = derived_projection_rows.get(compound.strip().lower(), {})
        target = target_map.get(compound.strip().lower(), {})
        score = float(assessment.score)
        dominant_factors: List[str] = []

        span_uncertainty = float(target.get("span_uncertainty", getattr(result, "avg_uncertainty", 5.0)) or getattr(result, "avg_uncertainty", 5.0))
        if span_uncertainty > 8.0:
            score -= 20.0
            dominant_factors.append(f"High pathway uncertainty ({span_uncertainty:.1f} kcal/mol).")
        elif span_uncertainty > 5.0:
            score -= 10.0
            dominant_factors.append(f"Moderate-high pathway uncertainty ({span_uncertainty:.1f} kcal/mol).")
        elif span_uncertainty > 3.0:
            score -= 4.0

        observable_ratio = float(meta.get("proxy_to_observable_ratio", 1.0) or 1.0)
        if observable_ratio < 0.2:
            score -= 18.0
            dominant_factors.append("Most of the proxy signal is lost after matrix/headspace projection.")
        elif observable_ratio < 0.5:
            score -= 10.0
            dominant_factors.append("Observable yield is materially lower than proxy yield.")

        matrix_factor = float(meta.get("matrix_factor", 1.0) or 1.0)
        headspace_factor = float(meta.get("headspace_factor", 1.0) or 1.0)
        accessibility_dominated = min(matrix_factor, headspace_factor) < 0.35
        if min(matrix_factor, headspace_factor) < 0.15:
            score -= 12.0
            dominant_factors.append("Severe physical suppression limits observability.")
        elif min(matrix_factor, headspace_factor) < 0.35:
            score -= 6.0
            dominant_factors.append("Accessibility (not chemistry) is the main source of uncertainty here.")

        score = _clamp_confidence_score(score)
        tier = _confidence_tier_from_score(score)
        rows.append({
            "compound": compound,
            "observable_ppb": float(meta.get("observable_ppb", 0.0)),
            "proxy_ppb": float(meta.get("proxy_ppb", 0.0)),
            "span_uncertainty_kcal": span_uncertainty,
            "tier": tier,
            "score": score,
            "prediction_mode": _prediction_mode_from_tier(tier),
            "dominant_factors": dominant_factors[:2],
            "calibration_source": str(meta.get("calibration_source", "unknown")),
            "calibration_evidence_strength": str(meta.get("calibration_evidence_strength", "heuristic")),
            "evidence_state": str(meta.get("evidence_state", "still_missing")),
            "target_class": str(meta.get("target_class", "unknown")),
            "matrix_factor": matrix_factor,
            "headspace_factor": headspace_factor,
            "accessibility_dominated": accessibility_dominated,
            "support_origin": str(derived_meta.get("support_origin", meta.get("support_origin", "standard_matrix_support"))),
            "observable_assumption_summary": str(
                derived_meta.get("observable_assumption_summary", meta.get("observable_assumption_summary", "unknown"))
            ),
            "chemically_reachable": bool(derived_meta.get("chemically_reachable", meta.get("chemically_reachable", False))),
            "reachability_status": str(
                derived_meta.get("reachability_status", meta.get("reachability_status", "merely_plausible"))
            ),
            "reachability_basis": str(
                derived_meta.get("reachability_basis", meta.get("reachability_basis", "mechanistic_surrogate_only"))
            ),
        })
    return rows


def _build_aggregate_confidence_rows(
    result: "FormulationResult",
    assessment: ConfidenceAssessment,
) -> Dict[str, Dict[str, object]]:
    rows: Dict[str, Dict[str, object]] = {}
    radar = getattr(result, "radar", {}) or {}
    for tag, raw in radar.items():
        score_value, support_count = raw
        score_value = float(score_value)
        support_count = int(support_count)
        if score_value <= 0.0:
            continue
        confidence_score = float(assessment.score)
        dominant_factors: List[str] = []
        if support_count <= 1:
            confidence_score -= 12.0
            dominant_factors.append("Aggregate sensory score is carried by a single compound signal.")
        elif support_count <= 3:
            confidence_score -= 5.0
            dominant_factors.append("Aggregate sensory score has limited compound support.")
        if score_value < 1.0:
            confidence_score -= 4.0
            dominant_factors.append("Aggregate score is low-intensity and may be less robust experimentally.")
        confidence_score = _clamp_confidence_score(confidence_score)
        tier = _confidence_tier_from_score(confidence_score)
        rows[tag] = {
            "score": score_value,
            "support_count": support_count,
            "tier": tier,
            "confidence_score": confidence_score,
            "prediction_mode": _prediction_mode_from_tier(tier),
            "dominant_factors": dominant_factors[:2],
        }
    return rows


def _fallback_sensitivity_summary(result: "FormulationResult") -> Dict[str, object]:
    precursor_attribution = getattr(result, "precursor_contributions", {}) or {}
    ranking_drivers = []
    total = sum(float(value) for value in precursor_attribution.values())
    for name, value in sorted(precursor_attribution.items(), key=lambda item: item[1], reverse=True)[:3]:
        share = 100.0 * float(value) / total if total > 0 else 0.0
        ranking_drivers.append({
            "input": name,
            "perturbation": "heuristic_attribution",
            "decision_delta": 0.0,
            "target_delta": 0.0,
            "safety_delta": 0.0,
            "off_flavour_delta": 0.0,
            "interpretation": f"Contributes {share:.1f}% of the observable aroma yield.",
        })

    safety_drivers = []
    for toxic in getattr(result, "flagged_toxics", [])[:3]:
        safety_drivers.append({
            "input": toxic,
            "perturbation": "detected_toxic_marker",
            "decision_delta": 0.0,
            "target_delta": 0.0,
            "safety_delta": float(result.safety_score),
            "off_flavour_delta": 0.0,
            "interpretation": "Detected toxic marker dominates the safety posture.",
        })

    return {
        "baseline_decision_score": _decision_proxy(result),
        "evaluated_perturbations": 0,
        "ranking_drivers": ranking_drivers,
        "safety_drivers": safety_drivers,
        "mode": "heuristic",
    }


def _analyze_formulation_sensitivity(
    result: "FormulationResult",
    *,
    formulation: Optional[Dict[str, object]] = None,
    baseline_conditions: Optional[ReactionConditions] = None,
    designer: Optional[MaillardPipeline] = None,
) -> Dict[str, object]:
    if not formulation or designer is None:
        return _fallback_sensitivity_summary(result)

    baseline_score = _decision_proxy(result)
    perturbations: List[Dict[str, object]] = []

    precursor_priority = {str(name).lower(): float(value) for name, value in (getattr(result, "precursor_contributions", {}) or {}).items()}
    precursor_inputs = _iter_precursor_inputs(formulation)
    precursor_inputs.sort(key=lambda item: precursor_priority.get(item[1].lower(), 0.0), reverse=True)

    for category, precursor in precursor_inputs[:6]:
        variant = deepcopy(formulation)
        variant_list = [item for item in (variant.get(category, []) or []) if str(item).strip() != precursor]
        if len(variant_list) == len(variant.get(category, []) or []):
            continue
        variant[category] = variant_list
        variant_conditions = _build_conditions_for_variant(variant, baseline_conditions)
        variant_result = designer.evaluate_single(variant, variant_conditions)
        perturbations.append({
            "input": precursor,
            "perturbation": f"remove_{category}",
            "decision_delta": _decision_proxy(variant_result) - baseline_score,
            "target_delta": float(variant_result.target_score) - float(result.target_score),
            "safety_delta": float(variant_result.safety_score) - float(result.safety_score),
            "off_flavour_delta": float(variant_result.off_flavour_risk) - float(result.off_flavour_risk),
            "interpretation": "Removal test against the baseline formulation.",
        })

    process_candidates = [
        ("temp", float(formulation.get("temp", getattr(baseline_conditions or ReactionConditions(), "temperature_celsius", 120.0))) + 10.0, "increase_temp"),
        ("temp", max(60.0, float(formulation.get("temp", getattr(baseline_conditions or ReactionConditions(), "temperature_celsius", 120.0))) - 10.0), "decrease_temp"),
        ("ph", min(9.5, float(formulation.get("ph", getattr(baseline_conditions or ReactionConditions(), "pH", 6.0))) + 0.5), "increase_ph"),
        ("ph", max(2.5, float(formulation.get("ph", getattr(baseline_conditions or ReactionConditions(), "pH", 6.0))) - 0.5), "decrease_ph"),
    ]
    if "time_minutes" in formulation:
        baseline_time = float(formulation.get("time_minutes", 60.0))
        process_candidates.extend([
            ("time_minutes", baseline_time * 1.25, "increase_time"),
            ("time_minutes", max(5.0, baseline_time * 0.8), "decrease_time"),
        ])

    for key, value, label in process_candidates:
        variant = deepcopy(formulation)
        variant[key] = value
        variant_conditions = _build_conditions_for_variant(variant, baseline_conditions)
        variant_result = designer.evaluate_single(variant, variant_conditions)
        perturbations.append({
            "input": key,
            "perturbation": label,
            "decision_delta": _decision_proxy(variant_result) - baseline_score,
            "target_delta": float(variant_result.target_score) - float(result.target_score),
            "safety_delta": float(variant_result.safety_score) - float(result.safety_score),
            "off_flavour_delta": float(variant_result.off_flavour_risk) - float(result.off_flavour_risk),
            "interpretation": f"Local process perturbation to {value:.2f}.",
        })

    ranking_drivers = sorted(perturbations, key=lambda item: abs(float(item["decision_delta"])), reverse=True)[:3]
    safety_drivers = sorted(perturbations, key=lambda item: abs(float(item["safety_delta"])), reverse=True)[:3]
    return {
        "baseline_decision_score": baseline_score,
        "evaluated_perturbations": len(perturbations),
        "ranking_drivers": ranking_drivers,
        "safety_drivers": safety_drivers,
        "mode": "local_oat",
    }


def build_confidence_package(
    result: "FormulationResult",
    warnings: List[DomainWarning],
    *,
    precursor_names: List[str],
    protein_type: str,
    formulation: Optional[Dict[str, object]] = None,
    baseline_conditions: Optional[ReactionConditions] = None,
    designer: Optional[MaillardPipeline] = None,
) -> Dict[str, object]:
    matrix_explainability = getattr(result, "matrix_explainability", {}) or {}
    temp_c = matrix_explainability.get("temperature_celsius")
    if temp_c is None and formulation is not None:
        temp_c = formulation.get("temp")
    if temp_c is None and baseline_conditions is not None:
        temp_c = getattr(baseline_conditions, "temperature_celsius", None)

    aw = None
    if formulation is not None:
        aw = formulation.get("aw")
    if aw is None and baseline_conditions is not None:
        aw = getattr(baseline_conditions, "water_activity", None)

    process_regime = _classify_process_regime(
        protein_type=protein_type,
        temp_c=None if temp_c is None else float(temp_c),
        aw=None if aw is None else float(aw),
    )
    assessment = assess_formulation_confidence(
        result,
        warnings,
        precursor_names=precursor_names,
        protein_type=protein_type,
    )
    assessment = _apply_prediction_mode_override(
        assessment,
        prediction_mode_override=process_regime.get("prediction_mode_override"),
        dominant_reason=str(process_regime.get("summary", "")) or None,
    )
    payload = asdict(assessment)
    payload["prediction_mode"] = _prediction_mode_from_tier(assessment.tier)
    payload["decision_mode"] = _decision_mode_from_prediction_mode(payload["prediction_mode"])
    payload["process_regime"] = process_regime.get("process_regime", "unknown")
    payload["process_neighborhood"] = process_regime.get("process_neighborhood", "unknown")
    payload["process_regime_summary"] = process_regime.get("summary", "")
    if baseline_conditions is not None and hasattr(baseline_conditions, "extrusion_profile"):
        extrusion_process = getattr(baseline_conditions, "extrusion_profile") or {}
        if extrusion_process.get("active"):
            payload["extrusion_process"] = extrusion_process
    payload["calibration_diagnostics"] = _build_calibration_diagnostics(assessment, result, warnings)
    payload["compound_confidence"] = _build_compound_confidence_rows(result, assessment)
    payload["aggregate_confidence"] = _build_aggregate_confidence_rows(result, assessment)
    payload["extrusion_observable_panel"] = _build_extrusion_observable_panel(result)

    if payload["process_regime"] in {"extrusion_like", "extrusion_heavy"}:
        panel = payload["extrusion_observable_panel"]
        if not bool(panel.get("minimum_panel_ready", False)):
            missing_categories = [
                category.replace("_", " ")
                for category in ("meaty_positive", "off_notes", "severity_markers")
                if not panel.get(category, {}).get("present")
            ]
            factor = "Extrusion observable panel is incomplete; missing live support for " + ", ".join(missing_categories) + "."
            dominant_factors = list(payload.get("dominant_factors", []))
            if factor not in dominant_factors:
                dominant_factors.append(factor)
            payload["dominant_factors"] = dominant_factors[:3]

            diagnostics = dict(payload.get("calibration_diagnostics", {}))
            axes = list(diagnostics.get("extrapolation_axes", []))
            if "extrusion_observable_panel" not in axes:
                axes.append("extrusion_observable_panel")
            diagnostics["extrapolation_axes"] = sorted(set(axes))
            diagnostics["supported_envelope"] = False
            diagnostics["summary"] = "Recommendation extrapolates beyond the strongest support on: " + ", ".join(diagnostics["extrapolation_axes"]) + "."
            payload["calibration_diagnostics"] = diagnostics

    payload["sensitivity_summary"] = _analyze_formulation_sensitivity(
        result,
        formulation=formulation,
        baseline_conditions=baseline_conditions,
        designer=designer,
    )
    return payload


def prepare_cli_confidence(
    result: "FormulationResult",
    *,
    target_tag: str,
    precursor_names: List[str],
    protein_type: str,
    temp_c: float,
    ph: float,
    aw: float,
    formulation: Dict[str, object],
    baseline_conditions: Optional[ReactionConditions] = None,
    designer: Optional[MaillardPipeline] = None,
) -> List[DomainWarning]:
    cleaned_precursors = [str(name).strip() for name in precursor_names if str(name).strip()]
    checker = DomainOfValidityChecker(target_tag)
    warnings = checker.check(
        precursor_names=cleaned_precursors,
        protein_type=protein_type,
        temp_c=temp_c,
        ph=ph,
        aw=aw,
        matrix_explainability=result.matrix_explainability,
    )
    result.confidence_metadata = build_confidence_package(
        result,
        warnings,
        precursor_names=cleaned_precursors,
        protein_type=protein_type,
        formulation=formulation,
        baseline_conditions=baseline_conditions,
        designer=designer,
    )
    return warnings


class DomainOfValidityChecker:
    """
    Checks if a formulation and conditions are within the "trusted envelope" 
    of validated benchmarks and physical assumptions.
    """
    
    TRUSTED_PRECURSORS = {
        "cysteine", "l-cysteine", "ribose", "d-ribose", 
        "glucose", "d-glucose", "asparagine", "l-asparagine"
    }
    
    MAX_TRUSTED_TEMP_C = 180.0
    MIN_TRUSTED_PH = 3.0
    MAX_TRUSTED_PH = 9.0

    def __init__(self, target_tag: str = "meaty"):
        self.target_tag = target_tag

    def check(
        self, 
        precursor_names: List[str], 
        protein_type: str, 
        temp_c: float, 
        ph: float,
        aw: Optional[float] = None,
        matrix_explainability: Optional[Dict[str, object]] = None,
    ) -> List[DomainWarning]:
        warnings = []
        matrix_expl = matrix_explainability if isinstance(matrix_explainability, dict) else {}
        
        # 1. Matrix Check
        if protein_type != "free":
            warnings.append(DomainWarning(
                level="CAUTION",
                category="MATRIX",
                message=f"Matrix '{protein_type}' uses speculative accessibility scaling; PRIMARY benchmarks are free-precursor only."
            ))

        if bool(matrix_expl.get("accessibility_warning", False)):
            accessibility_profile = str(matrix_expl.get("accessibility_profile", "unknown"))
            dominant_source = str(matrix_expl.get("accessibility_dominant_source", "unknown"))
            warnings.append(DomainWarning(
                level="CAUTION",
                category="ACCESSIBILITY",
                message=(
                    f"Accessibility state '{accessibility_profile}' is still controlling the prediction; "
                    f"reactive-site exposure is inferred via {dominant_source} rather than benchmarked release behavior."
                ),
            ))

        # 2. Precursor Support Check
        unsupported_precursors = []
        for name in precursor_names:
            norm_name = name.strip().lower()
            if norm_name not in self.TRUSTED_PRECURSORS:
                unsupported_precursors.append(name)
        
        if unsupported_precursors:
            warnings.append(DomainWarning(
                level="WARNING",
                category="PRECURSORS",
                message=f"Sparse benchmark analogies: {', '.join(unsupported_precursors)} lack PRIMARY quantitative validation."
            ))

        # 3. Process Condition Check
        if temp_c > self.MAX_TRUSTED_TEMP_C:
            warnings.append(DomainWarning(
                level="CAUTION",
                category="SEVERITY",
                message=f"Temperature {temp_c}°C exceeds the validated Arrhenius envelope (> {self.MAX_TRUSTED_TEMP_C}°C)."
            ))
            
        if ph < self.MIN_TRUSTED_PH or ph > self.MAX_TRUSTED_PH:
            warnings.append(DomainWarning(
                level="WARNING",
                category="PHYSICOCHEMICAL",
                message=f"pH {ph} is outside the trusted kinetic range ({self.MIN_TRUSTED_PH}-{self.MAX_TRUSTED_PH})."
            ))

        process_regime = _classify_process_regime(
            protein_type=protein_type,
            temp_c=temp_c,
            aw=aw,
        )
        override = process_regime.get("prediction_mode_override")
        if override == "directional_only":
            warnings.append(DomainWarning(
                level="CAUTION",
                category="EXTRUSION",
                message=str(process_regime.get("summary", "Extrusion-like processing reduces confidence to directional support.")),
            ))
        elif override == "hypothesis_only":
            warnings.append(DomainWarning(
                level="WARNING",
                category="EXTRUSION",
                message=str(process_regime.get("summary", "Extrusion-heavy processing is out of the validated neighborhood.")),
            ))

        return warnings


def assess_formulation_confidence(
    result: "FormulationResult",
    warnings: List[DomainWarning],
    *,
    precursor_names: List[str],
    protein_type: str,
) -> ConfidenceAssessment:
    trusted_precursors = DomainOfValidityChecker.TRUSTED_PRECURSORS
    normalized_precursors = [name.strip().lower() for name in precursor_names if isinstance(name, str) and name.strip()]
    supported_precursors = sum(1 for name in normalized_precursors if name in trusted_precursors)
    support_ratio = supported_precursors / len(normalized_precursors) if normalized_precursors else 1.0
    matrix_explainability = getattr(result, "matrix_explainability", {}) or {}

    score = 100.0
    dominant_factors: List[str] = []

    if protein_type == "free":
        if support_ratio >= 0.999:
            benchmark_neighborhood = "primary_free_precursor"
            dominant_factors.append("Matches the free-precursor benchmark family used by the PRIMARY quantitative envelope.")
        elif support_ratio >= 0.5:
            benchmark_neighborhood = "free_precursor_partial_analogy"
            score -= 18.0
            dominant_factors.append("Only part of the precursor set is covered by PRIMARY benchmark analogies.")
        else:
            benchmark_neighborhood = "sparse_precursor_analogy"
            score -= 32.0
            dominant_factors.append("Most precursors are outside the validated PRIMARY benchmark family.")
    elif protein_type in {"pea_iso", "soy_iso", "pea_conc", "soy_conc"}:
        benchmark_neighborhood = "matrix_intake_only"
        score -= 28.0
        dominant_factors.append("Plant-matrix support is still intake/headspace validated rather than target-ranking validated.")
    else:
        benchmark_neighborhood = "exploratory_matrix"
        score -= 38.0
        dominant_factors.append("This matrix family remains exploratory relative to the current validated envelope.")

    avg_uncertainty = float(getattr(result, "avg_uncertainty", 5.0) or 5.0)
    if avg_uncertainty > 8.0:
        score -= 24.0
        dominant_factors.append(f"Average pathway uncertainty is high ({avg_uncertainty:.1f} kcal/mol).")
    elif avg_uncertainty > 5.0:
        score -= 14.0
        dominant_factors.append(f"Average pathway uncertainty is moderate-high ({avg_uncertainty:.1f} kcal/mol).")
    elif avg_uncertainty > 3.0:
        score -= 6.0
        dominant_factors.append(f"Average pathway uncertainty is moderate ({avg_uncertainty:.1f} kcal/mol).")

    if bool(matrix_explainability.get("accessibility_warning", False)):
        score -= 10.0
        dominant_factors.append(
            "Prediction is dominated by accessibility assumptions rather than directly benchmarked release behavior."
        )

    warning_penalty = 0.0
    for warning in warnings:
        warning_penalty += 12.0 if warning.level == "WARNING" else 7.0
    if warning_penalty:
        score -= min(30.0, warning_penalty)
        dominant_factors.append("Domain-of-validity warnings indicate extrapolation beyond the most trusted scientific envelope.")

    score = _clamp_confidence_score(score)

    tier = _confidence_tier_from_score(score)
    recommended_posture = _recommended_posture_from_tier(tier)

    return ConfidenceAssessment(
        tier=tier,
        score=score,
        benchmark_neighborhood=benchmark_neighborhood,
        support_ratio=support_ratio,
        avg_uncertainty_kcal=avg_uncertainty,
        dominant_factors=dominant_factors[:3],
        recommended_posture=recommended_posture,
    )