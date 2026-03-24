from __future__ import annotations

import json
import math
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple

from src.mlp_external_benchmarks import build_external_mlp_evidence_index
from src.mlp_adoption_contract import MLPAdoptionDecision, MLPModelCandidate, load_mlp_candidates
from src.reaction_benchmark import ReactionBenchmarkEntry, build_reaction_benchmark_alias_index, load_reaction_benchmark_entries


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _normalize_key(value: str) -> str:
    text = str(value).strip().lower()
    for token in ("-", " ", "/"):
        text = text.replace(token, "_")
    return text


def _load_json(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _extract_predicted_barrier(raw_value: Any) -> Optional[float]:
    if raw_value is None:
        return None
    if isinstance(raw_value, (int, float)):
        return float(raw_value)
    if isinstance(raw_value, Mapping):
        candidate = raw_value.get("estimated_barrier_kcal")
        if candidate is None:
            return None
        return float(candidate)
    return None


def _rank_mapping(values: Mapping[str, float]) -> Dict[str, int]:
    ordered = sorted(values.items(), key=lambda item: (float(item[1]), item[0]))
    return {key: index + 1 for index, (key, _) in enumerate(ordered)}


def _spearman_rank_correlation(reference: Mapping[str, float], predicted: Mapping[str, float]) -> Optional[float]:
    comparable_keys = sorted(set(reference) & set(predicted))
    if len(comparable_keys) < 3:
        return None
    ref_ranks = _rank_mapping({key: float(reference[key]) for key in comparable_keys})
    pred_ranks = _rank_mapping({key: float(predicted[key]) for key in comparable_keys})
    n_terms = len(comparable_keys)
    diff_sq = sum((ref_ranks[key] - pred_ranks[key]) ** 2 for key in comparable_keys)
    return 1.0 - (6.0 * diff_sq) / float(n_terms * (n_terms**2 - 1))


def _load_candidate_predictions(candidate: MLPModelCandidate) -> Dict[str, float]:
    if not candidate.benchmark_results_path:
        return {}
    raw_path = Path(candidate.benchmark_results_path)
    path = raw_path if raw_path.is_absolute() else _repo_root() / raw_path
    if not path.exists():
        return {}
    payload = _load_json(path)
    predictions: Dict[str, float] = {}
    for key, value in payload.items():
        predicted = _extract_predicted_barrier(value)
        if predicted is not None:
            predictions[_normalize_key(str(key))] = float(predicted)
    return predictions


def _load_geometry_evidence(candidate: MLPModelCandidate) -> Dict[str, Any]:
    if not candidate.geometry_benchmark_path:
        return {"available": False}
    raw_path = Path(candidate.geometry_benchmark_path)
    path = raw_path if raw_path.is_absolute() else _repo_root() / raw_path
    if not path.exists():
        return {"available": False}
    payload = _load_json(path)
    if "candidate_assessments" in payload:
        for row in payload.get("candidate_assessments", []):
            if str(row.get("candidate_id", "")) == candidate.candidate_id:
                result = dict(row)
                result["available"] = bool(row.get("available", False))
                return result
        return {"available": False}
    payload["available"] = True
    return payload


def _nonphysical_threshold(entry: ReactionBenchmarkEntry) -> float:
    return max(120.0, float(entry.reference_barrier_kcal) + max(20.0, 5.0 * float(entry.barrier_uncertainty_kcal)))


def _external_prior_for(candidate: MLPModelCandidate) -> Dict[str, Any]:
    evidence_index = build_external_mlp_evidence_index()
    evidence = evidence_index.get(str(candidate.external_evidence_id or ""))
    if evidence is None:
        return {
            "external_evidence_id": candidate.external_evidence_id,
            "external_prior_strength": "unknown",
            "external_domain_relevance": "unknown",
            "external_selection_priority": "unknown",
            "external_supported_roles": [],
            "external_unproven_roles": [],
            "external_benchmark_scope": "unknown",
        }

    if candidate.proposed_role in evidence.supported_roles and evidence.prior_strength == "strong":
        selection_priority = "high"
    elif candidate.proposed_role in evidence.supported_roles and evidence.prior_strength == "moderate":
        selection_priority = "medium"
    elif evidence.prior_strength in {"strong", "moderate"}:
        selection_priority = "watch"
    else:
        selection_priority = "low"

    return {
        "external_evidence_id": evidence.evidence_id,
        "external_prior_strength": evidence.prior_strength,
        "external_domain_relevance": evidence.domain_relevance,
        "external_selection_priority": selection_priority,
        "external_supported_roles": list(evidence.supported_roles),
        "external_unproven_roles": list(evidence.unproven_roles),
        "external_benchmark_scope": evidence.benchmark_scope,
    }


def evaluate_candidate_against_reaction_benchmark(
    candidate: MLPModelCandidate,
    benchmark_entries: Iterable[ReactionBenchmarkEntry],
) -> Dict[str, Any]:
    benchmark_list = [entry for entry in benchmark_entries if not candidate.target_motif_families or entry.motif_family in candidate.target_motif_families]
    alias_index = build_reaction_benchmark_alias_index(benchmark_list)
    predictions = _load_candidate_predictions(candidate)
    geometry_evidence = _load_geometry_evidence(candidate)
    external_prior = _external_prior_for(candidate)

    comparable_rows: List[Dict[str, Any]] = []
    missing_reactions: List[str] = []
    nonphysical_reactions: List[str] = []
    reference_for_ordering: Dict[str, float] = {}
    predicted_for_ordering: Dict[str, float] = {}

    for entry in benchmark_list:
        predicted_barrier = None
        matched_alias = None
        for alias in entry.aliases:
            if alias in predictions:
                matched_alias = alias
                predicted_barrier = float(predictions[alias])
                break
        if predicted_barrier is None:
            missing_reactions.append(entry.reaction_family)
            continue

        nonphysical = (not math.isfinite(predicted_barrier)) or predicted_barrier < 0.0 or predicted_barrier > _nonphysical_threshold(entry)
        if nonphysical:
            nonphysical_reactions.append(entry.reaction_family)
        else:
            reference_for_ordering[entry.reaction_family] = float(entry.reference_barrier_kcal)
            predicted_for_ordering[entry.reaction_family] = float(predicted_barrier)

        comparable_rows.append(
            {
                "reaction_family": entry.reaction_family,
                "motif_family": entry.motif_family,
                "matched_alias": matched_alias,
                "reference_barrier_kcal": float(entry.reference_barrier_kcal),
                "predicted_barrier_kcal": float(predicted_barrier),
                "absolute_error_kcal": abs(float(predicted_barrier) - float(entry.reference_barrier_kcal)),
                "nonphysical": bool(nonphysical),
                "benchmark_visible_gap": entry.benchmark_visible_gap,
            }
        )

    comparable_valid = [row for row in comparable_rows if not row["nonphysical"]]
    coverage_ratio = len(comparable_rows) / max(len(benchmark_list), 1)
    mean_abs_error = None
    max_abs_error = None
    if comparable_valid:
        errors = [float(row["absolute_error_kcal"]) for row in comparable_valid]
        mean_abs_error = sum(errors) / len(errors)
        max_abs_error = max(errors)
    rank_correlation = _spearman_rank_correlation(reference_for_ordering, predicted_for_ordering)

    stop_reasons: List[str] = []
    hard_fail = False
    if candidate.materials_first:
        stop_reasons.append("materials_first_models_not_allowed_before_chemistry_domain_shortlist")
        hard_fail = True
    if nonphysical_reactions:
        stop_reasons.append("nonphysical_barrier_predictions")
        hard_fail = True
    if candidate.proposed_role == "local_barrier_surrogate":
        if coverage_ratio < 0.5:
            stop_reasons.append("insufficient_reaction_benchmark_coverage")
        if rank_correlation is None or rank_correlation < 0.8:
            stop_reasons.append("barrier_ordering_not_preserved")
            hard_fail = True
        if mean_abs_error is None or mean_abs_error > 6.0:
            stop_reasons.append("barrier_error_above_adoption_window")
            hard_fail = True
    if candidate.proposed_role in {"geom_preopt", "ts_initialization", "conformer_screening"}:
        if not geometry_evidence.get("available", False):
            stop_reasons.append("missing_geometry_benchmark_evidence")
        elif float(geometry_evidence.get("max_rmsd_angstrom", 999.0)) > 0.35:
            stop_reasons.append("geometry_drift_above_threshold")
            hard_fail = True

    if hard_fail:
        decision = "quarantine"
    elif stop_reasons:
        decision = "defer"
    else:
        decision = "adopt_offline"

    if decision == "adopt_offline":
        rationale = "Candidate stays inside the P4 chemistry benchmark window and is suitable only for the proposed offline accelerator role."
    elif decision == "defer":
        rationale = "Candidate remains unapproved until missing benchmark evidence is collected for the proposed offline role."
    else:
        rationale = "Candidate is quarantined because it fails core P4 stop rules and would reduce trust if integrated."

    if external_prior["external_prior_strength"] in {"strong", "moderate"}:
        rationale += (
            " External community benchmarks justify shortlist priority, but they do not replace Maillard-specific validation "
            "for proton-transfer, sulfur chemistry, TS behavior, or benchmark-visible barrier ordering."
        )

    return {
        "candidate_id": candidate.candidate_id,
        "model_family": candidate.model_family,
        "model_name": candidate.model_name,
        "proposed_role": candidate.proposed_role,
        "decision": decision,
        "coverage_ratio": float(coverage_ratio),
        "rank_correlation": None if rank_correlation is None else float(rank_correlation),
        "mean_abs_error_kcal": None if mean_abs_error is None else float(mean_abs_error),
        "max_abs_error_kcal": None if max_abs_error is None else float(max_abs_error),
        "nonphysical_reaction_count": len(nonphysical_reactions),
        "stop_reasons": stop_reasons,
        "rationale": rationale,
        "fallback_comparator": candidate.fallback_comparator,
        "benchmark_visible_gap": candidate.benchmark_visible_gap,
        "approved_for_default": False,
        "comparable_reactions": comparable_rows,
        "missing_reactions": missing_reactions,
        "likely_failure_modes": list(candidate.likely_failure_modes),
        "status": candidate.status,
        **external_prior,
    }


def build_mlp_assessment_artifact() -> Dict[str, Any]:
    benchmark_entries = load_reaction_benchmark_entries()
    candidates = load_mlp_candidates()
    candidate_rows = [evaluate_candidate_against_reaction_benchmark(candidate, benchmark_entries) for candidate in candidates]
    decision_counts: Dict[str, int] = {}
    for row in candidate_rows:
        decision_counts[row["decision"]] = decision_counts.get(row["decision"], 0) + 1

    return {
        "summary": {
            "benchmark_set_id": "maillard_reaction_benchmark_v1",
            "candidate_count": len(candidate_rows),
            "decision_counts": decision_counts,
            "default_policy": "no_default_mlp_until_reaction_benchmark_passes",
        },
        "candidates": candidate_rows,
    }


def build_adoption_decisions_from_assessment(payload: Mapping[str, Any]) -> List[MLPAdoptionDecision]:
    benchmark_set_id = str(payload.get("summary", {}).get("benchmark_set_id", "unknown"))
    decisions: List[MLPAdoptionDecision] = []
    for row in payload.get("candidates", []):
        decisions.append(
            MLPAdoptionDecision(
                candidate_id=str(row.get("candidate_id", "unknown")),
                model_family=str(row.get("model_family", "unknown")),
                model_name=str(row.get("model_name", "unknown")),
                proposed_role=str(row.get("proposed_role", "unknown")),
                decision=str(row.get("decision", "defer")),
                benchmark_set_id=benchmark_set_id,
                coverage_ratio=float(row.get("coverage_ratio", 0.0)),
                rank_correlation=None if row.get("rank_correlation") is None else float(row.get("rank_correlation")),
                mean_abs_error_kcal=None if row.get("mean_abs_error_kcal") is None else float(row.get("mean_abs_error_kcal")),
                max_abs_error_kcal=None if row.get("max_abs_error_kcal") is None else float(row.get("max_abs_error_kcal")),
                stop_reasons=[str(item) for item in row.get("stop_reasons", [])],
                rationale=str(row.get("rationale", "")),
                fallback_comparator=str(row.get("fallback_comparator", "unknown")),
                benchmark_visible_gap=str(row.get("benchmark_visible_gap", "unknown")),
                approved_for_default=bool(row.get("approved_for_default", False)),
            )
        )
    return decisions


def render_mlp_assessment_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# P4 MLP Assessment",
        "",
        "| Candidate | Role | Decision | External Prior | Priority | Coverage | Rank Correlation | MAE (kcal/mol) | Nonphysical Reactions | Stop Reasons |",
        "| --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in payload.get("candidates", []):
        correlation = row.get("rank_correlation")
        mae = row.get("mean_abs_error_kcal")
        lines.append(
            f"| {row.get('candidate_id', 'unknown')} | {row.get('proposed_role', 'unknown')} | {row.get('decision', 'unknown')} | "
            f"{row.get('external_prior_strength', 'unknown')} | {row.get('external_selection_priority', 'unknown')} | "
            f"{100.0 * float(row.get('coverage_ratio', 0.0)):.0f}% | "
            f"{'' if correlation is None else f'{float(correlation):.2f}'} | "
            f"{'' if mae is None else f'{float(mae):.2f}'} | {int(row.get('nonphysical_reaction_count', 0))} | "
            f"{', '.join(str(item) for item in row.get('stop_reasons', [])) or 'none'} |"
        )
    summary = payload.get("summary", {})
    lines.extend(
        [
            "",
            f"Benchmark set: {summary.get('benchmark_set_id', 'unknown')}",
            f"Candidates evaluated: {int(summary.get('candidate_count', 0))}",
            f"Default policy: {summary.get('default_policy', 'unknown')}",
        ]
    )
    return "\n".join(lines) + "\n"