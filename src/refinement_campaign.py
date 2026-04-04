from __future__ import annotations

import json
import os
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Optional, Tuple

from src.benchmark_validation import (
    DEFAULT_TARGET_TAG,
    _best_prediction_match,
    benchmark_to_formulation,
    build_matrix_target_status_artifact,
    evaluate_benchmark,
    get_benchmark_files,
    get_benchmark_metadata,
    get_matrix_ranking_contract,
    load_benchmark,
    summarize_evaluation,
)
from src.dft_refinement_contract import DFTTargetCandidate, build_offline_dft_job, triage_dft_candidate
from src.family_sensitivity import build_family_sensitivity_artifact
from src.refinement_watchlist import build_refinement_watchlist


DEFAULT_PROCESS_SCENARIOS: Tuple[Tuple[str, str, float], ...] = (
    ("temperature_celsius", "up", 10.0),
    ("temperature_celsius", "down", -10.0),
    ("pH", "up", 0.4),
    ("pH", "down", -0.4),
    ("time_minutes", "up", 15.0),
    ("time_minutes", "down", -15.0),
    ("water_activity", "up", 0.05),
    ("water_activity", "down", -0.05),
)

DEFAULT_FORMULATION_SCENARIOS: Tuple[Tuple[str, str, float], ...] = (
    ("sugars", "up", 1.25),
    ("sugars", "down", 0.75),
    ("amino_acids", "up", 1.25),
    ("amino_acids", "down", 0.75),
)
CHEAP_SCREENING_MAGNITUDES: Tuple[float, ...] = (1.0, 2.0, 3.0)


def _status_score(status: str) -> int:
    return {
        "pass": 0,
        "pass-no-ranking": 1,
        "coverage-gap": 2,
        "ranking-gap": 3,
        "scale-gap": 4,
        "unsupported": 5,
    }.get(str(status), 4)


def _improves_benchmark_decision(row: Mapping[str, Any]) -> bool:
    baseline_status = _status_score(str(row.get("baseline_status", row.get("scenario_status", "scale-gap"))))
    scenario_status = _status_score(str(row.get("scenario_status", "scale-gap")))
    if scenario_status < baseline_status:
        return True
    baseline_ranking = str(row.get("baseline_ranking_contract_status", "n/a"))
    scenario_ranking = str(row.get("scenario_ranking_contract_status", "n/a"))
    return baseline_ranking != "pass" and scenario_ranking == "pass"


def _execution_weight(execution_path: str) -> float:
    normalized = str(execution_path).strip().lower()
    if normalized == "matrix_precursor_augmented":
        return 1.0
    if normalized == "free_precursor":
        return 0.7
    return 0.5


def _requires_explicit_solvent(family: str) -> bool:
    normalized = str(family).strip().lower()
    return any(token in normalized for token in ["hydrolysis", "schiff", "rearrangement", "proton"])


def _requires_irc(family: str) -> bool:
    normalized = str(family).strip().lower()
    return any(token in normalized for token in ["rearrangement", "elimination", "strecker", "thiol_addition"])


@contextmanager
def _barrier_offsets(offsets: Mapping[str, float]) -> Iterator[None]:
    previous = os.environ.get("BARRIER_OFFSETS")
    try:
        os.environ["BARRIER_OFFSETS"] = json.dumps(dict(offsets))
        yield
    finally:
        if previous is None:
            os.environ.pop("BARRIER_OFFSETS", None)
        else:
            os.environ["BARRIER_OFFSETS"] = previous


def _evaluate_benchmark_payload(bench: Mapping[str, Any], *, target_tag: str) -> Tuple[Any, Any]:
    with tempfile.TemporaryDirectory() as td:
        bench_path = Path(td) / f"{bench.get('benchmark_id', 'temp_benchmark')}.json"
        bench_path.write_text(json.dumps(dict(bench), indent=2), encoding="utf-8")
        evaluation = evaluate_benchmark(bench_path, target_tag=target_tag)
        summary = summarize_evaluation(evaluation, protein_type=str(bench.get("protein_type", "free")))
        return evaluation, summary


def _build_benchmark_contexts(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> List[Dict[str, Any]]:
    bench_files = list(benchmark_files) if benchmark_files is not None else get_benchmark_files()
    contexts: List[Dict[str, Any]] = []
    for bench_file in bench_files:
        bench = load_benchmark(bench_file)
        metadata = get_benchmark_metadata(bench)
        if metadata.execution_path not in {"free_precursor", "matrix_precursor_augmented"}:
            continue
        evaluation = evaluate_benchmark(bench_file, target_tag=target_tag)
        summary = summarize_evaluation(evaluation, protein_type=str(bench.get("protein_type", "free")))
        contexts.append(
            {
                "bench": bench,
                "bench_file": Path(bench_file),
                "benchmark_id": summary.benchmark_id,
                "execution_path": summary.execution_path,
                "protein_type": str(bench.get("protein_type", "free")),
                "evaluation": evaluation,
                "summary": summary,
                "weight": _execution_weight(summary.execution_path),
            }
        )
    return contexts


def _observable_names(bench: Mapping[str, Any], evaluation: Any) -> List[str]:
    contract = get_matrix_ranking_contract(dict(bench))
    contract_targets = [str(item.get("name", "")).strip() for item in contract.get("observable_targets", []) if item.get("name")]
    if contract_targets:
        return contract_targets
    return [comparison.compound for comparison in evaluation.comparisons if comparison.matched_name is not None]


def _predicted_value_for_target(target_name: str, predicted_ppb: Mapping[str, float]) -> Optional[float]:
    matched_name, predicted_value, _score = _best_prediction_match(target_name, dict(predicted_ppb))
    if matched_name is None:
        return None
    return float(predicted_value)


def _panel_shift_score(bench: Mapping[str, Any], baseline_eval: Any, scenario_eval: Any) -> Tuple[float, int]:
    targets = _observable_names(bench, baseline_eval)
    if not targets:
        return 0.0, 0
    total_shift = 0.0
    terms = 0
    for target in targets:
        baseline_value = _predicted_value_for_target(target, baseline_eval.predicted_ppb)
        scenario_value = _predicted_value_for_target(target, scenario_eval.predicted_ppb)
        if baseline_value is None or scenario_value is None:
            continue
        total_shift += abs(scenario_value - baseline_value) / max(abs(baseline_value), 1.0e-12)
        terms += 1
    return total_shift / max(terms, 1), terms


def _benchmark_objective(summary: Any) -> float:
    score = 3.0 * _status_score(summary.overall_status)
    if str(summary.ranking_contract_status) not in {"pass", "n/a"}:
        score += 1.0
    if summary.mean_abs_log10_error is not None:
        score += min(float(summary.mean_abs_log10_error), 2.0)
    if summary.max_ratio is not None:
        score += max(0.0, float(summary.max_ratio) - 1.0)
    return float(score)


def _apply_process_scenario(bench: Mapping[str, Any], axis: str, magnitude: float) -> Dict[str, Any]:
    payload = json.loads(json.dumps(dict(bench)))
    conditions = payload.setdefault("conditions", {})
    if axis == "temperature_celsius":
        conditions["temp_C"] = max(20.0, float(conditions.get("temp_C", 25.0)) + magnitude)
    elif axis == "pH":
        conditions["ph"] = min(8.5, max(3.0, float(conditions.get("ph", 7.0)) + magnitude))
    elif axis == "time_minutes":
        conditions["time_min"] = max(1.0, float(conditions.get("time_min", 10.0)) + magnitude)
    elif axis == "water_activity":
        conditions["water_activity"] = min(0.99, max(0.2, float(conditions.get("water_activity", 0.9)) + magnitude))
    payload["benchmark_id"] = f"{payload.get('benchmark_id', 'benchmark')}__process__{axis}__{magnitude:+.2f}"
    return payload


def _apply_formulation_scenario(bench: Mapping[str, Any], bucket: str, factor: float) -> Dict[str, Any]:
    payload = json.loads(json.dumps(dict(bench)))
    formulation = benchmark_to_formulation(payload)
    bucket_names = list(formulation.get(bucket, []))
    if not bucket_names:
        return payload
    for name in bucket_names:
        if name in payload.get("precursors", {}):
            current = float(payload["precursors"][name].get("concentration_mM", 0.0))
            payload["precursors"][name]["concentration_mM"] = max(1.0e-6, current * factor)
    payload["benchmark_id"] = f"{payload.get('benchmark_id', 'benchmark')}__formulation__{bucket}__x{factor:.2f}"
    return payload


def _aggregate_axis_sensitivity(contexts: List[Dict[str, Any]], *, target_tag: str) -> Dict[str, List[Dict[str, Any]]]:
    process_rows: List[Dict[str, Any]] = []
    formulation_rows: List[Dict[str, Any]] = []

    for axis, direction, magnitude in DEFAULT_PROCESS_SCENARIOS:
        scenario_details: List[Dict[str, Any]] = []
        weighted_drift = 0.0
        ranking_changes = 0
        status_changes = 0
        touched = 0
        for context in contexts:
            scenario_bench = _apply_process_scenario(context["bench"], axis, magnitude)
            scenario_eval, scenario_summary = _evaluate_benchmark_payload(scenario_bench, target_tag=target_tag)
            panel_shift, terms = _panel_shift_score(context["bench"], context["evaluation"], scenario_eval)
            weighted_drift += context["weight"] * panel_shift
            if scenario_summary.ranking_contract_status != context["summary"].ranking_contract_status:
                ranking_changes += 1
            if scenario_summary.overall_status != context["summary"].overall_status:
                status_changes += 1
            if panel_shift > 0.05 or scenario_summary.overall_status != context["summary"].overall_status:
                touched += 1
            scenario_details.append(
                {
                    "benchmark_id": str(context["benchmark_id"]),
                    "panel_shift_score": float(panel_shift),
                    "matched_targets": int(terms),
                    "baseline_status": str(context["summary"].overall_status),
                    "scenario_status": str(scenario_summary.overall_status),
                    "baseline_ranking_contract_status": str(context["summary"].ranking_contract_status),
                    "scenario_ranking_contract_status": str(scenario_summary.ranking_contract_status),
                }
            )
        process_rows.append(
            {
                "axis": axis,
                "direction": direction,
                "magnitude": float(magnitude),
                "weighted_drift_score": float(weighted_drift),
                "ranking_change_count": int(ranking_changes),
                "status_change_count": int(status_changes),
                "touched_benchmark_count": int(touched),
                "benchmarks": scenario_details,
            }
        )

    for bucket, direction, factor in DEFAULT_FORMULATION_SCENARIOS:
        scenario_details = []
        weighted_drift = 0.0
        ranking_changes = 0
        status_changes = 0
        touched = 0
        for context in contexts:
            scenario_bench = _apply_formulation_scenario(context["bench"], bucket, factor)
            scenario_eval, scenario_summary = _evaluate_benchmark_payload(scenario_bench, target_tag=target_tag)
            panel_shift, terms = _panel_shift_score(context["bench"], context["evaluation"], scenario_eval)
            weighted_drift += context["weight"] * panel_shift
            if scenario_summary.ranking_contract_status != context["summary"].ranking_contract_status:
                ranking_changes += 1
            if scenario_summary.overall_status != context["summary"].overall_status:
                status_changes += 1
            if panel_shift > 0.05 or scenario_summary.overall_status != context["summary"].overall_status:
                touched += 1
            scenario_details.append(
                {
                    "benchmark_id": str(context["benchmark_id"]),
                    "panel_shift_score": float(panel_shift),
                    "matched_targets": int(terms),
                    "baseline_status": str(context["summary"].overall_status),
                    "scenario_status": str(scenario_summary.overall_status),
                    "baseline_ranking_contract_status": str(context["summary"].ranking_contract_status),
                    "scenario_ranking_contract_status": str(scenario_summary.ranking_contract_status),
                }
            )
        formulation_rows.append(
            {
                "bucket": bucket,
                "direction": direction,
                "factor": float(factor),
                "weighted_drift_score": float(weighted_drift),
                "ranking_change_count": int(ranking_changes),
                "status_change_count": int(status_changes),
                "touched_benchmark_count": int(touched),
                "benchmarks": scenario_details,
            }
        )

    process_rows.sort(key=lambda row: (-float(row["weighted_drift_score"]), row["axis"], row["direction"]))
    formulation_rows.sort(key=lambda row: (-float(row["weighted_drift_score"]), row["bucket"], row["direction"]))
    return {"process_axes": process_rows, "formulation_axes": formulation_rows}


def build_global_sensitivity_artifact(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> Dict[str, Any]:
    contexts = _build_benchmark_contexts(benchmark_files, target_tag=target_tag)
    barrier_payload = build_family_sensitivity_artifact(benchmark_files, target_tag=target_tag)
    axis_payload = _aggregate_axis_sensitivity(contexts, target_tag=target_tag)
    return {
        "summary": {
            "evaluated_benchmark_count": len(contexts),
            "barrier_family_count": int(barrier_payload.get("summary", {}).get("family_count", 0)),
            "process_axis_count": len(axis_payload["process_axes"]),
            "formulation_axis_count": len(axis_payload["formulation_axes"]),
            "top_barrier_family": (barrier_payload.get("families") or [{}])[0].get("reaction_family"),
            "top_process_axis": (axis_payload["process_axes"] or [{}])[0].get("axis"),
            "top_formulation_axis": (axis_payload["formulation_axes"] or [{}])[0].get("bucket"),
        },
        "barrier_families": list(barrier_payload.get("families", [])),
        "process_axes": axis_payload["process_axes"],
        "formulation_axes": axis_payload["formulation_axes"],
    }


def _evaluate_offsets(contexts: List[Dict[str, Any]], offsets: Mapping[str, float], *, target_tag: str) -> Dict[str, Any]:
    benchmark_rows: List[Dict[str, Any]] = []
    total_score = 0.0
    total_improvement = 0.0
    with _barrier_offsets(offsets):
        for context in contexts:
            evaluation = evaluate_benchmark(context["bench_file"], target_tag=target_tag)
            summary = summarize_evaluation(evaluation, protein_type=context["protein_type"])
            panel_shift, terms = _panel_shift_score(context["bench"], context["evaluation"], evaluation)
            baseline_score = _benchmark_objective(context["summary"])
            scenario_score = _benchmark_objective(summary)
            improvement = baseline_score - scenario_score
            total_score += context["weight"] * scenario_score
            total_improvement += context["weight"] * improvement
            benchmark_rows.append(
                {
                    "benchmark_id": str(context["benchmark_id"]),
                    "execution_path": str(context["execution_path"]),
                    "baseline_status": str(context["summary"].overall_status),
                    "scenario_status": str(summary.overall_status),
                    "baseline_ranking_contract_status": str(context["summary"].ranking_contract_status),
                    "scenario_ranking_contract_status": str(summary.ranking_contract_status),
                    "baseline_score": float(baseline_score),
                    "scenario_score": float(scenario_score),
                    "score_improvement": float(improvement),
                    "panel_shift_score": float(panel_shift),
                    "matched_targets": int(terms),
                }
            )
    return {
        "offsets": dict(offsets),
        "total_score": float(total_score),
        "total_improvement": float(total_improvement),
        "benchmarks": benchmark_rows,
    }


def build_cheap_screening_artifact(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> Dict[str, Any]:
    contexts = _build_benchmark_contexts(benchmark_files, target_tag=target_tag)
    global_payload = build_global_sensitivity_artifact(benchmark_files, target_tag=target_tag)
    matrix_status = build_matrix_target_status_artifact(benchmark_files, target_tag=target_tag)
    watchlist = build_refinement_watchlist(benchmark_files, target_tag=target_tag)
    mechanistic_benchmarks = {
        str(row.get("benchmark_id"))
        for row in matrix_status.get("benchmarks", [])
        if bool(row.get("mechanistic_priority_ready"))
    }
    watchlist_lookup = {str(row.get("reaction_family", "")).strip().lower(): row for row in watchlist.get("candidates", [])}
    baseline_total_score = sum(context["weight"] * _benchmark_objective(context["summary"]) for context in contexts)

    candidate_rows: List[Dict[str, Any]] = []
    family_rows = list(global_payload.get("barrier_families", []))
    max_family_impact = max((float(row.get("max_weighted_impact_score", 0.0)) for row in family_rows), default=1.0)

    for row in family_rows:
        family = str(row.get("reaction_family", "")).strip()
        if not family:
            continue
        scenarios = {str(item.get("direction")): item for item in row.get("scenarios", [])}
        scenario_trials: List[Tuple[str, Dict[str, float], Dict[str, Any]]] = []
        offset_key = str(row.get("offset_key", ""))
        if not offset_key:
            continue
        for direction_label, sign in (("down", -1.0), ("up", 1.0)):
            for magnitude in CHEAP_SCREENING_MAGNITUDES:
                recommended_offsets = {offset_key: sign * float(magnitude)}
                performance = _evaluate_offsets(contexts, recommended_offsets, target_tag=target_tag)
                scenario_trials.append((f"{direction_label}_{magnitude:.1f}", recommended_offsets, performance))
        if not scenario_trials:
            continue
        best_direction, recommended_offsets, performance = max(
            scenario_trials,
            key=lambda item: (float(item[2].get("total_improvement", 0.0)), -float(item[2].get("total_score", 0.0))),
        )
        linked_mechanistic = sorted(mechanistic_benchmarks.intersection(set(row.get("affected_benchmark_ids", []))))
        touched_benchmarks = [item for item in performance["benchmarks"] if item["panel_shift_score"] > 0.05 or item["score_improvement"] != 0.0]
        linked_benchmark_rows = [item for item in performance["benchmarks"] if item["benchmark_id"] in set(linked_mechanistic)]
        linked_decision_change = any(_improves_benchmark_decision(item) for item in linked_benchmark_rows)
        linked_activity = any(item["panel_shift_score"] > 0.05 or item["score_improvement"] != 0.0 for item in linked_benchmark_rows)
        if not touched_benchmarks:
            screen_decision = "reject"
            no_escalation_reason = "cheap-first perturbation does not move the target panel or benchmark diagnostics"
        elif linked_mechanistic and linked_decision_change:
            screen_decision = "advance"
            no_escalation_reason = "none"
        elif linked_mechanistic and linked_activity:
            screen_decision = "defer"
            no_escalation_reason = "cheap-first perturbation moves mechanistic-priority surrogates but does not change the mechanistic-priority benchmark decisions"
        elif linked_mechanistic and performance["total_improvement"] > 0.05:
            screen_decision = "defer"
            no_escalation_reason = "cheap-first perturbation improves non-priority diagnostics but does not move the mechanistic-priority benchmark decisions"
        elif performance["total_improvement"] > 0.05:
            screen_decision = "advance"
            no_escalation_reason = "none"
        else:
            screen_decision = "defer"
            no_escalation_reason = "cheap-first perturbation is sensitivity-visible but does not improve benchmark-visible decisions"
        watch = watchlist_lookup.get(family.lower(), {})
        candidate_rows.append(
            {
                "reaction_family": family,
                "offset_key": str(row.get("offset_key", "")),
                "recommended_offsets": recommended_offsets,
                "dominant_direction": str(row.get("dominant_direction", "none")),
                "best_screen_direction": best_direction,
                "max_weighted_impact_score": float(row.get("max_weighted_impact_score", 0.0)),
                "normalized_sensitivity_score": float(row.get("max_weighted_impact_score", 0.0)) / max(max_family_impact, 1.0e-9),
                "affected_benchmark_ids": list(row.get("affected_benchmark_ids", [])),
                "linked_mechanistic_benchmarks": linked_mechanistic,
                "baseline_total_score": float(baseline_total_score),
                "scenario_total_score": float(performance["total_score"]),
                "total_improvement": float(performance["total_improvement"]),
                "screen_decision": screen_decision,
                "no_escalation_reason": no_escalation_reason,
                "screening_method": "barrier_offset_surrogate",
                "uncertainty_posture": "directional_until_qm",
                "current_barrier_source": str(watch.get("current_barrier_source", "unknown")),
                "benchmark_impact_score": float(watch.get("benchmark_impact_score", min(1.0, len(row.get("affected_benchmark_ids", [])) / max(len(contexts), 1)))),
                "evidence_gap_score": float(watch.get("evidence_gap_score", 0.8)),
                "performance": performance,
            }
        )

    candidate_rows.sort(key=lambda item: ({"advance": 0, "defer": 1, "reject": 2}.get(item["screen_decision"], 3), -float(item["total_improvement"]), -float(item["max_weighted_impact_score"]), item["reaction_family"]))
    return {
        "summary": {
            "candidate_count": len(candidate_rows),
            "advance": sum(1 for row in candidate_rows if row["screen_decision"] == "advance"),
            "defer": sum(1 for row in candidate_rows if row["screen_decision"] == "defer"),
            "reject": sum(1 for row in candidate_rows if row["screen_decision"] == "reject"),
            "mechanistic_priority_benchmark_count": len(mechanistic_benchmarks),
        },
        "candidates": candidate_rows,
    }


def build_selective_dft_plan(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> Dict[str, Any]:
    cheap_payload = build_cheap_screening_artifact(benchmark_files, target_tag=target_tag)
    rows: List[Dict[str, Any]] = []
    offline_jobs: List[Dict[str, Any]] = []
    for candidate_row in cheap_payload.get("candidates", []):
        family = str(candidate_row.get("reaction_family", ""))
        if candidate_row.get("screen_decision") == "advance":
            dft_candidate = DFTTargetCandidate(
                reaction_id=f"family::{family.lower()}",
                reaction_family=family,
                sensitivity_score=float(candidate_row.get("normalized_sensitivity_score", 0.0)),
                benchmark_impact_score=float(candidate_row.get("benchmark_impact_score", 0.0)),
                evidence_gap_score=float(candidate_row.get("evidence_gap_score", 0.0)),
                current_barrier_source=str(candidate_row.get("current_barrier_source", "unknown")),
                rationale=(
                    f"Cheap-first screening found a benchmark-visible improvement of {float(candidate_row.get('total_improvement', 0.0)):.2f} "
                    f"for {family} with offsets {candidate_row.get('recommended_offsets', {})}."
                ),
                requires_explicit_solvent=_requires_explicit_solvent(family),
                requires_irc=_requires_irc(family),
                current_confidence="heuristic" if str(candidate_row.get("current_barrier_source", "")).lower().startswith("heuristic") else "mixed",
            )
            decision = triage_dft_candidate(dft_candidate)
            row = {
                "reaction_family": family,
                "reaction_id": dft_candidate.reaction_id,
                "recommended_offsets": dict(candidate_row.get("recommended_offsets", {})),
                "cheap_screen_decision": str(candidate_row.get("screen_decision", "unknown")),
                "total_improvement": float(candidate_row.get("total_improvement", 0.0)),
                "decision": decision.decision,
                "priority_tier": decision.priority_tier,
                "rationale": decision.rationale,
                "required_outputs": list(decision.required_outputs),
                "surrogate_plan": decision.surrogate_plan,
            }
        else:
            row = {
                "reaction_family": family,
                "reaction_id": f"family::{family.lower()}",
                "recommended_offsets": dict(candidate_row.get("recommended_offsets", {})),
                "cheap_screen_decision": str(candidate_row.get("screen_decision", "unknown")),
                "total_improvement": float(candidate_row.get("total_improvement", 0.0)),
                "decision": "defer" if candidate_row.get("screen_decision") == "defer" else "reject",
                "priority_tier": "cheap_first_hold" if candidate_row.get("screen_decision") == "defer" else "not_applicable",
                "rationale": str(candidate_row.get("no_escalation_reason", "cheap-first screening did not support escalation")),
                "required_outputs": [],
                "surrogate_plan": "Preserve current surrogate barrier until a cheap-first screen produces benchmark-visible improvement.",
            }
        rows.append(row)
        if row["decision"] == "run_now":
            job = build_offline_dft_job(dft_candidate, artifact_prefix=f"p3_{family.lower()}")
            offline_jobs.append(
                {
                    "reaction_id": job.reaction_id,
                    "priority_tier": job.priority_tier,
                    "solvent_name": job.solvent_name,
                    "temperature_k": job.temperature_k,
                    "charge": job.charge,
                    "spin": job.spin,
                    "use_explicit_solvent": job.use_explicit_solvent,
                    "n_water": job.n_water,
                    "generate_irc": job.generate_irc,
                    "expected_outputs": list(job.expected_outputs),
                    "artifact_prefix": job.artifact_prefix,
                }
            )
    rows.sort(key=lambda item: ({"run_now": 0, "defer": 1, "reject": 2}.get(item["decision"], 3), -float(item["total_improvement"]), item["reaction_family"]))
    return {
        "summary": {
            "candidate_count": len(rows),
            "run_now": sum(1 for row in rows if row["decision"] == "run_now"),
            "defer": sum(1 for row in rows if row["decision"] == "defer"),
            "reject": sum(1 for row in rows if row["decision"] == "reject"),
        },
        "candidates": rows,
        "offline_jobs": offline_jobs,
    }


def _has_severe_regression(baseline_rows: List[Dict[str, Any]], scenario_rows: List[Dict[str, Any]]) -> bool:
    baseline_lookup = {row["benchmark_id"]: row for row in baseline_rows}
    for row in scenario_rows:
        baseline = baseline_lookup.get(row["benchmark_id"])
        if baseline is None:
            continue
        baseline_status = _status_score(str(baseline.get("baseline_status", baseline.get("scenario_status", "scale-gap"))))
        scenario_status = _status_score(str(row.get("scenario_status", "scale-gap")))
        if scenario_status - baseline_status >= 2:
            return True
        if str(baseline.get("baseline_ranking_contract_status", "n/a")) == "pass" and str(row.get("scenario_ranking_contract_status", "n/a")) != "pass":
            return True
    return False


def build_refinement_impact_artifact(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> Dict[str, Any]:
    contexts = _build_benchmark_contexts(benchmark_files, target_tag=target_tag)
    cheap_payload = build_cheap_screening_artifact(benchmark_files, target_tag=target_tag)
    baseline = _evaluate_offsets(contexts, {}, target_tag=target_tag)
    accepted_offsets: Dict[str, float] = {}
    accepted_candidates: List[Dict[str, Any]] = []
    current = baseline

    for candidate in cheap_payload.get("candidates", []):
        if candidate.get("screen_decision") != "advance":
            continue
        trial_offsets = dict(accepted_offsets)
        trial_offsets.update({str(key): float(value) for key, value in candidate.get("recommended_offsets", {}).items()})
        trial = _evaluate_offsets(contexts, trial_offsets, target_tag=target_tag)
        if trial["total_score"] >= current["total_score"] - 0.05:
            continue
        if _has_severe_regression(current["benchmarks"], trial["benchmarks"]):
            continue
        accepted_offsets = trial_offsets
        accepted_candidates.append(
            {
                "reaction_family": str(candidate.get("reaction_family", "unknown")),
                "recommended_offsets": dict(candidate.get("recommended_offsets", {})),
                "screen_decision": str(candidate.get("screen_decision", "unknown")),
                "total_improvement": float(candidate.get("total_improvement", 0.0)),
            }
        )
        current = trial

    patch_payload = {
        "schema_version": "1.0",
        "description": "Accepted surrogate barrier offsets selected only when they improve benchmark-visible diagnostics without severe regressions.",
        "accepted_offsets": accepted_offsets,
        "accepted_candidates": accepted_candidates,
        "baseline_total_score": float(baseline["total_score"]),
        "patched_total_score": float(current["total_score"]),
    }
    return {
        "patch": patch_payload,
        "baseline": baseline,
        "patched": current,
        "summary": {
            "accepted_candidate_count": len(accepted_candidates),
            "baseline_total_score": float(baseline["total_score"]),
            "patched_total_score": float(current["total_score"]),
            "total_improvement": float(baseline["total_score"] - current["total_score"]),
        },
    }


def render_global_sensitivity_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# P3 Global Sensitivity",
        "",
        "## Barrier Families",
        "",
        "| Reaction Family | Dominant Direction | Max Impact | Affected Benchmarks |",
        "| --- | --- | ---: | --- |",
    ]
    for row in payload.get("barrier_families", []):
        lines.append(
            f"| {row['reaction_family']} | {row.get('dominant_direction', 'none')} | {float(row.get('max_weighted_impact_score', 0.0)):.2f} | {', '.join(row.get('affected_benchmark_ids', [])[:4]) or 'none'} |"
        )
    lines.extend([
        "",
        "## Process Axes",
        "",
        "| Axis | Direction | Magnitude | Weighted Drift | Touched Benchmarks |",
        "| --- | --- | ---: | ---: | ---: |",
    ])
    for row in payload.get("process_axes", [])[:8]:
        lines.append(
            f"| {row['axis']} | {row['direction']} | {float(row['magnitude']):.2f} | {float(row['weighted_drift_score']):.2f} | {int(row['touched_benchmark_count'])} |"
        )
    lines.extend([
        "",
        "## Formulation Axes",
        "",
        "| Bucket | Direction | Factor | Weighted Drift | Touched Benchmarks |",
        "| --- | --- | ---: | ---: | ---: |",
    ])
    for row in payload.get("formulation_axes", [])[:8]:
        lines.append(
            f"| {row['bucket']} | {row['direction']} | {float(row['factor']):.2f} | {float(row['weighted_drift_score']):.2f} | {int(row['touched_benchmark_count'])} |"
        )
    summary = payload.get("summary", {})
    lines.extend([
        "",
        f"Benchmarks evaluated: {int(summary.get('evaluated_benchmark_count', 0))}",
        f"Top barrier family: {summary.get('top_barrier_family', 'n/a')}",
        f"Top process axis: {summary.get('top_process_axis', 'n/a')}",
        f"Top formulation axis: {summary.get('top_formulation_axis', 'n/a')}",
    ])
    return "\n".join(lines) + "\n"


def render_cheap_screening_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Cheap Refinement Screening",
        "",
        "| Reaction Family | Screen Decision | Improvement | Direction | Offsets | Mechanistic Benchmarks |",
        "| --- | --- | ---: | --- | --- | --- |",
    ]
    for row in payload.get("candidates", []):
        lines.append(
            f"| {row['reaction_family']} | {row['screen_decision']} | {float(row['total_improvement']):.2f} | {row.get('best_screen_direction', row['dominant_direction'])} | {json.dumps(row['recommended_offsets'], sort_keys=True)} | {', '.join(row.get('linked_mechanistic_benchmarks', [])[:4]) or 'none'} |"
        )
    summary = payload.get("summary", {})
    lines.extend([
        "",
        f"Advance candidates: {int(summary.get('advance', 0))}",
        f"Deferred candidates: {int(summary.get('defer', 0))}",
        f"Rejected candidates: {int(summary.get('reject', 0))}",
    ])
    return "\n".join(lines) + "\n"


def render_selective_dft_plan_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Selective DFT Plan",
        "",
        "| Reaction Family | DFT Decision | Tier | Cheap Improvement | Offsets |",
        "| --- | --- | --- | ---: | --- |",
    ]
    for row in payload.get("candidates", []):
        lines.append(
            f"| {row['reaction_family']} | {row['decision']} | {row['priority_tier']} | {float(row['total_improvement']):.2f} | {json.dumps(row['recommended_offsets'], sort_keys=True)} |"
        )
    summary = payload.get("summary", {})
    lines.extend([
        "",
        f"Run-now jobs: {int(summary.get('run_now', 0))}",
        f"Deferred jobs: {int(summary.get('defer', 0))}",
        f"Rejected jobs: {int(summary.get('reject', 0))}",
    ])
    return "\n".join(lines) + "\n"


def render_refinement_impact_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Refinement Impact",
        "",
        "| Benchmark | Baseline Status | Patched Status | Baseline Ranking | Patched Ranking | Score Improvement |",
        "| --- | --- | --- | --- | --- | ---: |",
    ]
    patched_lookup = {row['benchmark_id']: row for row in payload.get('patched', {}).get('benchmarks', [])}
    for baseline_row in payload.get("baseline", {}).get("benchmarks", []):
        patched_row = patched_lookup.get(baseline_row["benchmark_id"], baseline_row)
        lines.append(
            f"| {baseline_row['benchmark_id']} | {baseline_row['scenario_status']} | {patched_row['scenario_status']} | {baseline_row['scenario_ranking_contract_status']} | {patched_row['scenario_ranking_contract_status']} | {float(patched_row.get('score_improvement', 0.0)):.2f} |"
        )
    summary = payload.get("summary", {})
    lines.extend([
        "",
        f"Accepted candidates: {int(summary.get('accepted_candidate_count', 0))}",
        f"Total improvement: {float(summary.get('total_improvement', 0.0)):.2f}",
        f"Baseline total score: {float(summary.get('baseline_total_score', 0.0)):.2f}",
        f"Patched total score: {float(summary.get('patched_total_score', 0.0)):.2f}",
    ])
    return "\n".join(lines) + "\n"