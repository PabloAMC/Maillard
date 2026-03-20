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


def _clamp_confidence_score(score: float) -> float:
    return max(5.0, min(100.0, float(score)))


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
            "matrix_factor": matrix_factor,
            "headspace_factor": headspace_factor,
            "accessibility_dominated": accessibility_dominated,
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
    assessment = assess_formulation_confidence(
        result,
        warnings,
        precursor_names=precursor_names,
        protein_type=protein_type,
    )
    payload = asdict(assessment)
    payload["prediction_mode"] = _prediction_mode_from_tier(assessment.tier)
    payload["calibration_diagnostics"] = _build_calibration_diagnostics(assessment, result, warnings)
    payload["compound_confidence"] = _build_compound_confidence_rows(result, assessment)
    payload["aggregate_confidence"] = _build_aggregate_confidence_rows(result, assessment)
    payload["sensitivity_summary"] = _analyze_formulation_sensitivity(
        result,
        formulation=formulation,
        baseline_conditions=baseline_conditions,
        designer=designer,
    )
    return payload


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
        aw: float = 1.0
    ) -> List[DomainWarning]:
        warnings = []
        
        # 1. Matrix Check
        if protein_type != "free":
            warnings.append(DomainWarning(
                level="CAUTION",
                category="MATRIX",
                message=f"Matrix '{protein_type}' uses speculative accessibility scaling; PRIMARY benchmarks are free-precursor only."
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

    warning_penalty = 0.0
    for warning in warnings:
        warning_penalty += 12.0 if warning.level == "WARNING" else 7.0
    if warning_penalty:
        score -= min(30.0, warning_penalty)
        dominant_factors.append("Domain-of-validity warnings indicate extrapolation beyond the most trusted scientific envelope.")

    score = _clamp_confidence_score(score)

    tier = _confidence_tier_from_score(score)
    if tier == "high":
        recommended_posture = "Suitable for quantitative prioritization before wet-lab confirmation."
    elif tier == "medium":
        recommended_posture = "Reliable for ranking and formulation triage; verify absolute levels experimentally."
    elif tier == "low":
        recommended_posture = "Use directionally only; absolute concentrations should be treated as provisional."
    else:
        recommended_posture = "Use for hypothesis generation, not decision-grade prioritization."

    return ConfidenceAssessment(
        tier=tier,
        score=score,
        benchmark_neighborhood=benchmark_neighborhood,
        support_ratio=support_ratio,
        avg_uncertainty_kcal=avg_uncertainty,
        dominant_factors=dominant_factors[:3],
        recommended_posture=recommended_posture,
    )