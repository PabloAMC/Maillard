"""
src/benchmark_reporting.py — Benchmark summary, artifact, and reporting functions.

Extracted from src/benchmark_validation.py (Priority 2 decomposition).
Contains all functions that produce BenchmarkSummary objects, build complex
JSON artifacts (evidence audits, promotion contracts, target-status payloads,
family-lane validation), and render markdown tables.  No IO, metadata
resolution, or raw prediction logic lives here.
"""
from __future__ import annotations

from copy import deepcopy
import math
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional

from src.benchmark_types import (
    BenchmarkEvaluation,
    BenchmarkIndexEntry,
    BenchmarkSummary,
    BenchmarkTargetSnapshot,
    BenchmarkThresholds,
    CompoundComparison,
    MatrixBenchmarkBranchDelta,
    MatrixBenchmarkDelta,
    MatrixBenchmarkEvidence,
    MatrixPromotionFamilyStatus,
    ThermodynamicGatingAudit,
)
from src.benchmark_registry import (
    DEFAULT_TARGET_TAG,
    get_benchmark_files,
    get_benchmark_metadata,
    get_matrix_ranking_contract,
    load_benchmark,
)


def _benchmark_cache_key(benchmark_files: Optional[Iterable[Path | str]]) -> tuple[str, ...]:
    files = list(benchmark_files) if benchmark_files is not None else get_benchmark_files()
    return tuple(str(Path(file_path).resolve()) for file_path in files)
from src.benchmark_evaluator import (
    evaluate_benchmark,
    _best_prediction_match,
    _mean_abs_log10_error,
    _resolve_scale_thresholds,
    _run_benchmark_recommendation,
    MATRIX_BENCHMARK_PROFILES,
)
from src.text_utils import normalize_compound_name as _normalize_name
from src.matrix_targets import get_compound_panel_entry
from src.literature_family_registry import resolve_family_descriptor
from src.validation_contract import DEFAULT_VALIDATION_CONTRACT

DEFAULT_BENCHMARK_THRESHOLDS = DEFAULT_VALIDATION_CONTRACT.thresholds
THERMODYNAMIC_GATING_MIN_ABSOLUTE_MAE_IMPROVEMENT_PPB = 5.0
THERMODYNAMIC_GATING_MIN_RELATIVE_MAE_IMPROVEMENT = 0.05
THERMODYNAMIC_GATING_MIN_RATIO_IMPROVEMENT = 0.05


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _projection_metadata_for_match(
    evaluation: BenchmarkEvaluation,
    comparison: CompoundComparison,
) -> Dict[str, Any]:
    if comparison.matched_name and comparison.matched_name in evaluation.projection_metadata:
        return evaluation.projection_metadata[comparison.matched_name]
    normalized = _normalize_name(comparison.compound)
    for key, meta in evaluation.projection_metadata.items():
        if (
            _normalize_name(str(key)) == normalized
            or _normalize_name(str(meta.get("compound", ""))) == normalized
        ):
            return meta
    return {}


def _evaluate_matrix_ranking_contract(
    bench: dict,
    predicted_ppb: Dict[str, float],
) -> Dict[str, object]:
    contract = get_matrix_ranking_contract(bench)
    observable_targets = contract.get("observable_targets", [])
    if not observable_targets:
        return {
            "status": "missing_contract",
            "ranked_observable_targets": [],
            "adverse_markers": contract.get("adverse_markers", []),
        }

    expected = sorted(
        observable_targets,
        key=lambda item: int(item.get("expected_rank", 999)),
    )
    predicted_rows = []
    missing_targets = []
    for item in expected:
        matched_name, predicted_value, _score = _best_prediction_match(
            str(item.get("name", "")), predicted_ppb
        )
        if matched_name is None:
            missing_targets.append(str(item.get("name", "")))
            continue
        predicted_rows.append((str(item.get("name", "")), float(predicted_value)))

    if missing_targets:
        status = "missing_targets"
    else:
        predicted_order = [
            name
            for name, _value in sorted(predicted_rows, key=lambda row: row[1], reverse=True)
        ]
        expected_order = [str(item.get("name", "")) for item in expected]
        status = "pass" if predicted_order == expected_order else "order_mismatch"

    return {
        "status": status,
        "ranked_observable_targets": [str(item.get("name", "")) for item in expected],
        "adverse_markers": contract.get("adverse_markers", []),
    }


# ---------------------------------------------------------------------------
# Summary builders
# ---------------------------------------------------------------------------

def summarize_evaluation_for_benchmark(
    evaluation: BenchmarkEvaluation,
    bench: dict,
    *,
    protein_type: str = "free",
    thresholds: BenchmarkThresholds = DEFAULT_BENCHMARK_THRESHOLDS,
) -> BenchmarkSummary:
    matched = [c for c in evaluation.comparisons if c.matched_name is not None]
    ratios = [c.ratio for c in matched if math.isfinite(c.ratio)]
    max_ratio = max(ratios) if ratios else None
    mean_ratio = sum(ratios) / len(ratios) if ratios else None
    mean_abs_log10_error = _mean_abs_log10_error(matched)

    scale_thresholds = _resolve_scale_thresholds(bench, protein_type=protein_type, thresholds=thresholds)
    ratio_threshold = scale_thresholds["max_ratio"]
    log_error_threshold = scale_thresholds["mean_abs_log10_error"]
    metadata = get_benchmark_metadata(bench)
    conditions = bench.get("conditions", {})
    ranking_contract = (
        _evaluate_matrix_ranking_contract(bench, evaluation.predicted_ppb)
        if metadata.execution_path in {"matrix_only", "matrix_precursor_augmented"}
        else {"status": "n/a", "ranked_observable_targets": [], "adverse_markers": []}
    )
    matrix_contract = (
        get_matrix_ranking_contract(bench)
        if metadata.execution_path in {"matrix_only", "matrix_precursor_augmented"}
        else {}
    )
    process_state = matrix_contract.get("process_state")

    if not evaluation.supported:
        return BenchmarkSummary(
            benchmark_id=evaluation.benchmark_id,
            bench_file=evaluation.bench_file,
            tier=metadata.tier,
            family=metadata.family,
            execution_path=metadata.execution_path,
            benchmark_engine=metadata.benchmark_engine,
            comparator_signal=metadata.comparator_signal,
            cantera_role=metadata.cantera_role,
            target_snapshot_policy=metadata.target_snapshot_policy,
            thermodynamic_gating_policy=metadata.thermodynamic_gating_policy,
            supported=False,
            reason=evaluation.reason,
            protein_type=protein_type,
            coverage=0.0,
            matched_compounds=0,
            total_compounds=0,
            pearson_r=None,
            mae_ppb=None,
            max_ratio=None,
            mean_ratio=None,
            ranking_status="unsupported",
            scale_status="unsupported",
            overall_status="unsupported",
            strict_ready=False,
            blocking_issues=[evaluation.reason or "unsupported"],
            conditions=conditions,
            process_state=process_state,
            ranked_observable_targets=list(ranking_contract.get("ranked_observable_targets", [])),
            adverse_markers=list(ranking_contract.get("adverse_markers", [])),
            ranking_contract_status=str(ranking_contract.get("status", "n/a")),
            calibration_mode=matrix_contract.get("calibration_mode"),
            reference_signal_origin=evaluation.reference_signal_origin,
            mean_abs_log10_error=None,
        )

    if len(matched) >= thresholds.min_matched_for_ranking and evaluation.pearson_r is not None:
        ranking_status = "pass" if evaluation.pearson_r >= thresholds.ranking_threshold else "fail"
    elif len(matched) > 0:
        ranking_status = "n/a"
    else:
        ranking_status = "fail"

    if max_ratio is None:
        scale_status = "fail"
    elif max_ratio <= ratio_threshold and (
        mean_abs_log10_error is None or mean_abs_log10_error <= log_error_threshold
    ):
        scale_status = "pass"
    else:
        scale_status = "fail"

    overall_status = "pass"
    if evaluation.coverage < thresholds.full_coverage_threshold:
        overall_status = "coverage-gap"
    elif ranking_status == "fail":
        overall_status = "ranking-gap"
    elif scale_status == "fail":
        overall_status = "scale-gap"
    elif ranking_status == "n/a":
        overall_status = "pass-no-ranking"

    blocking_issues: List[str] = []
    if evaluation.coverage < thresholds.full_coverage_threshold:
        blocking_issues.append(
            f"coverage {evaluation.coverage:.1%} < {thresholds.full_coverage_threshold:.0%}"
        )
    if ranking_status == "fail":
        pearson = "n/a" if evaluation.pearson_r is None else f"{evaluation.pearson_r:.3f}"
        blocking_issues.append(f"ranking {pearson} < {thresholds.ranking_threshold:.2f}")
    if scale_status == "fail":
        ratio_value = "n/a" if max_ratio is None else f"{max_ratio:.3f}"
        blocking_issues.append(f"max ratio {ratio_value} > {ratio_threshold:.2f}")
        if mean_abs_log10_error is not None and mean_abs_log10_error > log_error_threshold:
            blocking_issues.append(
                f"mean |log10 ratio| {mean_abs_log10_error:.3f} > {log_error_threshold:.3f}"
            )

    strict_ready = (
        evaluation.coverage >= thresholds.full_coverage_threshold
        and ranking_status != "fail"
        and scale_status == "pass"
    )
    if not DEFAULT_VALIDATION_CONTRACT.is_strict_gate_eligible(
        tier=metadata.tier, execution_path=metadata.execution_path
    ):
        strict_ready = False
        if not blocking_issues:
            blocking_issues.append(
                "matrix-only intake path is executable but not yet in the strict release gate"
            )

    if (
        metadata.execution_path in {"matrix_only", "matrix_precursor_augmented"}
        and ranking_contract.get("status") not in {"pass", "n/a"}
    ):
        blocking_issues.append(f"matrix ranking contract: {ranking_contract.get('status')}")

    return BenchmarkSummary(
        benchmark_id=evaluation.benchmark_id,
        bench_file=evaluation.bench_file,
        tier=metadata.tier,
        family=metadata.family,
        execution_path=metadata.execution_path,
        benchmark_engine=metadata.benchmark_engine,
        comparator_signal=metadata.comparator_signal,
        cantera_role=metadata.cantera_role,
        target_snapshot_policy=metadata.target_snapshot_policy,
        thermodynamic_gating_policy=metadata.thermodynamic_gating_policy,
        supported=True,
        reason=None,
        protein_type=protein_type,
        coverage=evaluation.coverage,
        matched_compounds=len(matched),
        total_compounds=len(evaluation.comparisons),
        pearson_r=evaluation.pearson_r,
        mae_ppb=evaluation.mae_ppb,
        max_ratio=max_ratio,
        mean_ratio=mean_ratio,
        ranking_status=ranking_status,
        scale_status=scale_status,
        overall_status=overall_status,
        strict_ready=strict_ready,
        blocking_issues=blocking_issues,
        conditions=conditions,
        process_state=process_state,
        ranked_observable_targets=list(ranking_contract.get("ranked_observable_targets", [])),
        adverse_markers=list(ranking_contract.get("adverse_markers", [])),
        ranking_contract_status=str(ranking_contract.get("status", "n/a")),
        calibration_mode=matrix_contract.get("calibration_mode"),
        reference_signal_origin=evaluation.reference_signal_origin,
        mean_abs_log10_error=mean_abs_log10_error,
    )


def summarize_evaluation(
    evaluation: BenchmarkEvaluation,
    *,
    protein_type: str = "free",
    thresholds: BenchmarkThresholds = DEFAULT_BENCHMARK_THRESHOLDS,
) -> BenchmarkSummary:
    bench = load_benchmark(evaluation.bench_file)
    return summarize_evaluation_for_benchmark(
        evaluation, bench, protein_type=protein_type, thresholds=thresholds
    )


# ---------------------------------------------------------------------------
# Matrix evidence helpers
# ---------------------------------------------------------------------------

def _matrix_source_origin(bench: dict) -> str:
    source_metadata = bench.get("source_metadata") or {}
    origin = str(source_metadata.get("origin", "")).strip()
    if origin:
        return origin
    if bench.get("source_doi"):
        return "external_literature"
    return "unspecified"


def _matrix_source_reference(bench: dict) -> str:
    source_metadata = bench.get("source_metadata") or {}
    if bench.get("source_doi"):
        return str(bench["source_doi"])
    generator = str(source_metadata.get("generator", "")).strip()
    origin = str(source_metadata.get("origin", "")).strip()
    if origin and generator:
        return f"{origin}:{generator}"
    return origin or generator or "unspecified"


def _matrix_target_profile(bench: dict) -> str:
    contract = get_matrix_ranking_contract(bench)
    roles = {str(item.get("role", "")).strip().lower() for item in contract.get("observable_targets", [])}
    has_desirable = "desirable_marker" in roles
    has_adverse = bool(contract.get("adverse_markers"))
    if has_desirable and has_adverse:
        return "mixed"
    if has_desirable:
        return "meaty_positive"
    if has_adverse:
        return "adverse_only"
    return "untyped"


def _matrix_external_data_status(bench: dict) -> str:
    source_origin = _matrix_source_origin(bench)
    has_measured = bool(bench.get("measured_volatiles"))
    if has_measured and (bench.get("source_doi") or source_origin.startswith("external")):
        return "external_quantitative"
    if has_measured:
        return "quantitative_unspecified_origin"
    if bench.get("reference_volatiles"):
        return "internal_reference_only"
    return "no_comparator_signal"


def assess_matrix_benchmark_evidence(bench: Any) -> MatrixBenchmarkEvidence:
    from src.benchmark_registry import ROOT as _ROOT
    if isinstance(bench, (Path, str)):
        bench_path = Path(bench)
        bench = load_benchmark(bench_path)
    else:
        bench_path = _ROOT / "data" / "benchmarks" / f"{bench.get('benchmark_id', 'unknown')}.json"

    metadata = get_benchmark_metadata(bench)
    process_state = get_matrix_ranking_contract(bench).get("process_state")
    reference_signal_origin = (
        "measured_volatiles" if bench.get("measured_volatiles") else "reference_volatiles"
    )
    source_origin = _matrix_source_origin(bench)
    target_profile = _matrix_target_profile(bench)
    external_data_status = _matrix_external_data_status(bench)
    promotable = (
        metadata.execution_path in {"matrix_only", "matrix_precursor_augmented"}
        and target_profile in {"meaty_positive", "mixed"}
        and external_data_status == "external_quantitative"
        and reference_signal_origin == "measured_volatiles"
    )

    if target_profile == "adverse_only":
        blocker = "benchmark only anchors adverse/off-flavour markers; no external meaty-positive targets are present"
    elif external_data_status != "external_quantitative":
        blocker = "missing external quantitative matrix evidence for meaty-positive targets"
    elif reference_signal_origin != "measured_volatiles":
        blocker = "comparator signal is not wet-lab measured_volatiles"
    else:
        blocker = ""

    return MatrixBenchmarkEvidence(
        benchmark_id=str(bench.get("benchmark_id", bench_path.stem)),
        bench_file=bench_path,
        protein_type=str(bench.get("protein_type", "free")),
        execution_path=metadata.execution_path,
        process_state=process_state,
        reference_signal_origin=reference_signal_origin,
        source_origin=source_origin,
        source_reference=_matrix_source_reference(bench),
        target_profile=target_profile,
        external_data_status=external_data_status,
        promotable=promotable,
        promotion_blocker=blocker,
    )


def build_matrix_benchmark_evidence_audit(
    benchmark_files: Optional[Iterable[Path | str]] = None,
) -> List[MatrixBenchmarkEvidence]:
    bench_files = list(benchmark_files) if benchmark_files is not None else get_benchmark_files()
    rows: List[MatrixBenchmarkEvidence] = []
    for bench_file in bench_files:
        bench = load_benchmark(bench_file)
        metadata = get_benchmark_metadata(bench)
        if metadata.execution_path not in {"matrix_only", "matrix_precursor_augmented"}:
            continue
        rows.append(assess_matrix_benchmark_evidence(Path(bench_file)))
    return rows


# ---------------------------------------------------------------------------
# Delta builders
# ---------------------------------------------------------------------------

def build_matrix_benchmark_deltas(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> List[MatrixBenchmarkDelta]:
    bench_files = list(benchmark_files) if benchmark_files is not None else get_benchmark_files()
    rows: List[MatrixBenchmarkDelta] = []
    for bench_file in bench_files:
        bench = load_benchmark(bench_file)
        metadata = get_benchmark_metadata(bench)
        if metadata.execution_path not in {"matrix_only", "matrix_precursor_augmented"}:
            continue
        evaluation = evaluate_benchmark(bench_file, target_tag=target_tag)
        summary = summarize_evaluation(
            evaluation, protein_type=bench.get("protein_type", "free")
        )
        role_lookup = {
            str(item.get("name", "")).strip().lower(): str(item.get("role", "target"))
            for item in get_matrix_ranking_contract(bench).get("observable_targets", [])
        }
        adverse = {str(item).strip().lower() for item in summary.adverse_markers}
        for comparison in evaluation.comparisons:
            meta = _projection_metadata_for_match(evaluation, comparison)
            compound_key = comparison.compound.strip().lower()
            role = role_lookup.get(
                compound_key, "adverse_marker" if compound_key in adverse else "reference"
            )
            abs_delta = abs(float(comparison.measured_ppb) - float(comparison.predicted_ppb))
            pct_delta = (
                None
                if comparison.measured_ppb <= 0.0
                else abs_delta / float(comparison.measured_ppb)
            )
            rows.append(
                MatrixBenchmarkDelta(
                    benchmark_id=evaluation.benchmark_id,
                    bench_file=evaluation.bench_file,
                    protein_type=bench.get("protein_type", "free"),
                    execution_path=metadata.execution_path,
                    process_state=summary.process_state,
                    reference_signal_origin=evaluation.reference_signal_origin,
                    ranking_contract_status=summary.ranking_contract_status,
                    compound=comparison.compound,
                    role=role,
                    reference_ppb=float(comparison.measured_ppb),
                    predicted_ppb=float(comparison.predicted_ppb),
                    abs_delta_ppb=abs_delta,
                    pct_delta=pct_delta,
                    ratio=float(comparison.ratio),
                    calibration_source=str(meta.get("calibration_source", "class_fallback")),
                    calibration_evidence_strength=str(
                        meta.get("calibration_evidence_strength", "heuristic")
                    ),
                    calibration_fallback_mode=str(meta.get("calibration_fallback_mode", "class_level")),
                )
            )
    return rows


# ---------------------------------------------------------------------------
# Benchmark collection-level summarizer
# ---------------------------------------------------------------------------

def _benchmark_compound_names(bench: dict, summary: BenchmarkSummary) -> List[str]:
    names: List[str] = []
    measured = bench.get("measured_volatiles", {}) or {}
    if isinstance(measured, dict):
        names.extend(str(name) for name in measured.keys())
    reference = bench.get("reference_volatiles", {}) or {}
    if isinstance(reference, dict):
        names.extend(str(name) for name in reference.keys())
    names.extend(str(name) for name in summary.ranked_observable_targets)
    names.extend(str(name) for name in summary.adverse_markers)
    deduped: List[str] = []
    seen: set[str] = set()
    for name in names:
        normalized = _normalize_name(name)
        if not normalized or normalized in seen:
            continue
        seen.add(normalized)
        deduped.append(str(name))
    return deduped


def _payload_role_from_evidence_state(evidence_state: str) -> str:
    normalized = str(evidence_state).strip().lower()
    if normalized in {"externally_benchmarked", "internally_benchmarked"}:
        return "benchmark_payload"
    if normalized == "conditional_calibration":
        return "calibration_payload"
    if normalized in {"transferred_prior", "safety_reference"}:
        return "directional_prior"
    return "structural_gap_extrapolation"


def _enrich_benchmark_summary_family_metadata(
    summary: BenchmarkSummary, bench: dict
) -> BenchmarkSummary:
    chemistry_families: List[str] = []
    slr_families: List[str] = []
    family_lane_names: List[str] = []
    payload_roles: List[str] = []

    for compound_name in _benchmark_compound_names(bench, summary):
        panel_entry = get_compound_panel_entry(compound_name) or {}
        chemistry_family = str(panel_entry.get("chemistry_family", "")).strip()
        if not chemistry_family:
            continue
        descriptor = resolve_family_descriptor(chemistry_family)
        chemistry_families.append(chemistry_family)
        if descriptor:
            slr_families.append(str(descriptor.get("slr_family", "")).zfill(2))
            family_lane_names.append(str(descriptor.get("display_name", chemistry_family)))
        payload_roles.append(
            _payload_role_from_evidence_state(
                str(panel_entry.get("evidence_state", "still_missing"))
            )
        )

    return BenchmarkSummary(
        **{k: v for k, v in summary.__dict__.items() if k not in {
            "chemistry_families", "slr_families", "payload_roles", "family_lane_names"
        }},
        chemistry_families=sorted(dict.fromkeys(chemistry_families)),
        slr_families=sorted(dict.fromkeys(slr_families)),
        payload_roles=sorted(dict.fromkeys(payload_roles)),
        family_lane_names=sorted(dict.fromkeys(family_lane_names)),
    )


def summarize_benchmarks(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    target_tag: str = DEFAULT_TARGET_TAG,
    thresholds: BenchmarkThresholds = DEFAULT_BENCHMARK_THRESHOLDS,
) -> List[BenchmarkSummary]:
    bench_files = list(benchmark_files) if benchmark_files is not None else get_benchmark_files()
    summaries: List[BenchmarkSummary] = []
    for bench_file in bench_files:
        bench_path = Path(bench_file)
        bench = load_benchmark(bench_path)
        evaluation = evaluate_benchmark(bench_path, target_tag=target_tag)
        summary = summarize_evaluation(
            evaluation, protein_type=bench.get("protein_type", "free"), thresholds=thresholds
        )
        summaries.append(_enrich_benchmark_summary_family_metadata(summary, bench))
    return summaries


# ---------------------------------------------------------------------------
# Branch delta comparison
# ---------------------------------------------------------------------------

def compare_matrix_benchmark_delta_sets(
    current_rows: Iterable[MatrixBenchmarkDelta],
    baseline_rows: Iterable[MatrixBenchmarkDelta],
    *,
    current_evidence: Optional[Iterable[MatrixBenchmarkEvidence]] = None,
    baseline_evidence: Optional[Iterable[MatrixBenchmarkEvidence]] = None,
) -> List[MatrixBenchmarkBranchDelta]:
    current_lookup = {(row.benchmark_id, _normalize_name(row.compound)): row for row in current_rows}
    baseline_lookup = {(row.benchmark_id, _normalize_name(row.compound)): row for row in baseline_rows}
    current_evidence_lookup = {row.benchmark_id: row for row in (current_evidence or [])}
    baseline_evidence_lookup = {row.benchmark_id: row for row in (baseline_evidence or [])}

    rows: List[MatrixBenchmarkBranchDelta] = []
    for key in sorted(set(current_lookup) | set(baseline_lookup)):
        current = current_lookup.get(key)
        baseline = baseline_lookup.get(key)
        current_meta = current_evidence_lookup.get(key[0])
        baseline_meta = baseline_evidence_lookup.get(key[0])

        if current is None:
            change_type = "removed"
        elif baseline is None:
            change_type = "added"
        else:
            predicted_changed = not math.isclose(
                current.predicted_ppb, baseline.predicted_ppb, rel_tol=1.0e-9, abs_tol=1.0e-12
            )
            ratio_changed = not math.isclose(
                current.ratio, baseline.ratio, rel_tol=1.0e-9, abs_tol=1.0e-12
            )
            metadata_changed = (
                current.execution_path != baseline.execution_path
                or current.reference_signal_origin != baseline.reference_signal_origin
                or (current_meta.source_origin if current_meta else "n/a")
                != (baseline_meta.source_origin if baseline_meta else "n/a")
                or (current_meta.external_data_status if current_meta else "n/a")
                != (baseline_meta.external_data_status if baseline_meta else "n/a")
            )
            if predicted_changed or ratio_changed:
                change_type = "modified"
            elif metadata_changed:
                change_type = "metadata_changed"
            else:
                continue

        predicted_delta = (
            None if current is None or baseline is None else current.predicted_ppb - baseline.predicted_ppb
        )
        ratio_delta = (
            None if current is None or baseline is None else current.ratio - baseline.ratio
        )
        rows.append(
            MatrixBenchmarkBranchDelta(
                benchmark_id=key[0],
                compound=current.compound if current is not None else baseline.compound,
                change_type=change_type,
                current_present=current is not None,
                baseline_present=baseline is not None,
                current_execution_path=current.execution_path if current is not None else "n/a",
                baseline_execution_path=baseline.execution_path if baseline is not None else "n/a",
                current_reference_signal_origin=(
                    current.reference_signal_origin if current is not None else "n/a"
                ),
                baseline_reference_signal_origin=(
                    baseline.reference_signal_origin if baseline is not None else "n/a"
                ),
                current_source_origin=(
                    current_meta.source_origin if current_meta is not None else "n/a"
                ),
                baseline_source_origin=(
                    baseline_meta.source_origin if baseline_meta is not None else "n/a"
                ),
                current_external_data_status=(
                    current_meta.external_data_status if current_meta is not None else "n/a"
                ),
                baseline_external_data_status=(
                    baseline_meta.external_data_status if baseline_meta is not None else "n/a"
                ),
                current_predicted_ppb=current.predicted_ppb if current is not None else None,
                baseline_predicted_ppb=baseline.predicted_ppb if baseline is not None else None,
                predicted_delta_ppb=predicted_delta,
                current_ratio=current.ratio if current is not None else None,
                baseline_ratio=baseline.ratio if baseline is not None else None,
                ratio_delta=ratio_delta,
            )
        )
    return rows


def render_matrix_branch_deltas_markdown(
    rows: Iterable[MatrixBenchmarkBranchDelta],
    *,
    base_ref: str,
) -> str:
    from src.benchmark_markdown import render_matrix_branch_deltas_markdown as _render_matrix_branch_deltas_markdown

    return _render_matrix_branch_deltas_markdown(rows, base_ref=base_ref)


# ---------------------------------------------------------------------------
# Family lane validation
# ---------------------------------------------------------------------------

def build_family_lane_validation_artifact(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
    thresholds: BenchmarkThresholds = DEFAULT_BENCHMARK_THRESHOLDS,
) -> Dict[str, Any]:
    summaries = summarize_benchmarks(benchmark_files, target_tag=target_tag, thresholds=thresholds)
    enriched_summaries: List[BenchmarkSummary] = []
    for summary in summaries:
        bench = load_benchmark(summary.bench_file)
        enriched_summaries.append(_enrich_benchmark_summary_family_metadata(summary, bench))

    family_groups: Dict[str, Dict[str, Any]] = {}
    lane_groups: Dict[str, Dict[str, Any]] = defaultdict(lambda: {
        "execution_path": "unknown",
        "benchmark_count": 0,
        "strict_ready_count": 0,
        "supported_count": 0,
        "status_counts": defaultdict(int),
        "chemistry_families": set(),
        "payload_roles": set(),
    })

    for summary in enriched_summaries:
        lane_bucket = lane_groups[summary.execution_path]
        lane_bucket["execution_path"] = summary.execution_path
        lane_bucket["benchmark_count"] += 1
        lane_bucket["strict_ready_count"] += int(bool(summary.strict_ready))
        lane_bucket["supported_count"] += int(bool(summary.supported))
        lane_bucket["status_counts"][summary.overall_status] += 1
        lane_bucket["chemistry_families"].update(summary.chemistry_families)
        lane_bucket["payload_roles"].update(summary.payload_roles)

        for chemistry_family in summary.chemistry_families:
            descriptor = resolve_family_descriptor(chemistry_family)
            family_bucket = family_groups.setdefault(
                chemistry_family,
                {
                    "chemistry_family": chemistry_family,
                    "slr_family": str(descriptor.get("slr_family", "")).zfill(2) if descriptor else "",
                    "display_name": str(descriptor.get("display_name", chemistry_family)) if descriptor else chemistry_family,
                    "strategic_posture": str(descriptor.get("strategic_posture", "unknown")) if descriptor else "unknown",
                    "benchmark_count": 0,
                    "strict_ready_count": 0,
                    "supported_count": 0,
                    "status_counts": defaultdict(int),
                    "execution_paths": defaultdict(int),
                    "payload_roles": set(),
                    "benchmark_ids": [],
                },
            )
            family_bucket["benchmark_count"] += 1
            family_bucket["strict_ready_count"] += int(bool(summary.strict_ready))
            family_bucket["supported_count"] += int(bool(summary.supported))
            family_bucket["status_counts"][summary.overall_status] += 1
            family_bucket["execution_paths"][summary.execution_path] += 1
            family_bucket["payload_roles"].update(summary.payload_roles)
            family_bucket["benchmark_ids"].append(summary.benchmark_id)

    family_rows = []
    for chemistry_family, row in sorted(
        family_groups.items(), key=lambda item: (item[1]["slr_family"], item[0])
    ):
        family_rows.append({
            "chemistry_family": chemistry_family,
            "slr_family": row["slr_family"],
            "display_name": row["display_name"],
            "strategic_posture": row["strategic_posture"],
            "benchmark_count": int(row["benchmark_count"]),
            "strict_ready_count": int(row["strict_ready_count"]),
            "supported_count": int(row["supported_count"]),
            "status_counts": dict(sorted(row["status_counts"].items())),
            "execution_paths": dict(sorted(row["execution_paths"].items())),
            "payload_roles": sorted(row["payload_roles"]),
            "benchmark_ids": sorted(dict.fromkeys(row["benchmark_ids"])),
        })

    lane_rows = []
    for execution_path, row in sorted(lane_groups.items()):
        lane_rows.append({
            "execution_path": execution_path,
            "benchmark_count": int(row["benchmark_count"]),
            "strict_ready_count": int(row["strict_ready_count"]),
            "supported_count": int(row["supported_count"]),
            "status_counts": dict(sorted(row["status_counts"].items())),
            "chemistry_families": sorted(row["chemistry_families"]),
            "payload_roles": sorted(row["payload_roles"]),
        })

    return {
        "summary": {
            "benchmark_count": len(enriched_summaries),
            "family_count": len(family_rows),
            "lane_count": len(lane_rows),
            "tracked_execution_paths": [row["execution_path"] for row in lane_rows],
        },
        "families": family_rows,
        "lanes": lane_rows,
        "benchmarks": [
            {
                "benchmark_id": summary.benchmark_id,
                "execution_path": summary.execution_path,
                "overall_status": summary.overall_status,
                "strict_ready": bool(summary.strict_ready),
                "chemistry_families": list(summary.chemistry_families),
                "slr_families": list(summary.slr_families),
                "payload_roles": list(summary.payload_roles),
                "family_lane_names": list(summary.family_lane_names),
            }
            for summary in enriched_summaries
        ],
    }


# ---------------------------------------------------------------------------
# Thermodynamic gating audit
# ---------------------------------------------------------------------------

def _benchmark_status_score(summary: BenchmarkSummary) -> int:
    ranking = {
        "unsupported": 0,
        "coverage-gap": 1,
        "ranking-gap": 2,
        "scale-gap": 3,
        "partial-pass": 4,
        "pass-no-ranking": 5,
        "pass": 5,
    }
    return ranking.get(summary.overall_status, 0)


def thermodynamic_gating_materially_improves(
    baseline: BenchmarkSummary, gated: BenchmarkSummary
) -> bool:
    if not baseline.supported or not gated.supported:
        return False
    if gated.coverage + 1.0e-12 < baseline.coverage:
        return False
    if _benchmark_status_score(gated) < _benchmark_status_score(baseline):
        return False

    baseline_mae = baseline.mae_ppb if baseline.mae_ppb is not None else float("inf")
    gated_mae = gated.mae_ppb if gated.mae_ppb is not None else float("inf")
    baseline_ratio = baseline.max_ratio if baseline.max_ratio is not None else float("inf")
    gated_ratio = gated.max_ratio if gated.max_ratio is not None else float("inf")

    mae_improvement = baseline_mae - gated_mae
    ratio_improvement = baseline_ratio - gated_ratio
    mae_threshold = max(
        THERMODYNAMIC_GATING_MIN_ABSOLUTE_MAE_IMPROVEMENT_PPB,
        (
            baseline_mae * THERMODYNAMIC_GATING_MIN_RELATIVE_MAE_IMPROVEMENT
            if math.isfinite(baseline_mae)
            else float("inf")
        ),
    )
    return (
        (math.isfinite(mae_improvement) and mae_improvement >= mae_threshold)
        or (math.isfinite(ratio_improvement) and ratio_improvement >= THERMODYNAMIC_GATING_MIN_RATIO_IMPROVEMENT)
    )


def audit_thermodynamic_gating(
    bench_file: Path | str,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
    thresholds: BenchmarkThresholds = DEFAULT_BENCHMARK_THRESHOLDS,
) -> ThermodynamicGatingAudit:
    bench_path = Path(bench_file)
    bench = load_benchmark(bench_path)
    metadata = get_benchmark_metadata(bench)
    protein_type = bench.get("protein_type", "free")

    if metadata.execution_path != "free_precursor" or metadata.family == "safety":
        return ThermodynamicGatingAudit(
            benchmark_id=bench["benchmark_id"],
            bench_file=bench_path,
            execution_path=metadata.execution_path,
            applicable=False,
            baseline_overall_status="n/a",
            gated_overall_status="n/a",
            baseline_mae_ppb=None,
            gated_mae_ppb=None,
            baseline_max_ratio=None,
            gated_max_ratio=None,
            delta_mae_ppb=None,
            delta_max_ratio=None,
            material_improvement=False,
            recommended_policy="diagnostic_only",
            notes="Thermodynamic gating audit is currently only meaningful for non-safety free-precursor FAST benchmarks.",
        )

    baseline_eval = evaluate_benchmark(bench_path, target_tag=target_tag, thermodynamic_gating="off")
    gated_eval = evaluate_benchmark(bench_path, target_tag=target_tag, thermodynamic_gating="on")
    baseline_summary = summarize_evaluation(baseline_eval, protein_type=protein_type, thresholds=thresholds)
    gated_summary = summarize_evaluation(gated_eval, protein_type=protein_type, thresholds=thresholds)

    delta_mae = (
        None
        if baseline_summary.mae_ppb is None or gated_summary.mae_ppb is None
        else baseline_summary.mae_ppb - gated_summary.mae_ppb
    )
    delta_ratio = (
        None
        if baseline_summary.max_ratio is None or gated_summary.max_ratio is None
        else baseline_summary.max_ratio - gated_summary.max_ratio
    )
    material = thermodynamic_gating_materially_improves(baseline_summary, gated_summary)

    return ThermodynamicGatingAudit(
        benchmark_id=bench["benchmark_id"],
        bench_file=bench_path,
        execution_path=metadata.execution_path,
        applicable=True,
        baseline_overall_status=baseline_summary.overall_status,
        gated_overall_status=gated_summary.overall_status,
        baseline_mae_ppb=baseline_summary.mae_ppb,
        gated_mae_ppb=gated_summary.mae_ppb,
        baseline_max_ratio=baseline_summary.max_ratio,
        gated_max_ratio=gated_summary.max_ratio,
        delta_mae_ppb=delta_mae,
        delta_max_ratio=delta_ratio,
        material_improvement=material,
        recommended_policy="benchmark_facing_candidate" if material else "diagnostic_only",
        notes=(
            "Thermodynamic gating materially improves benchmark error without degrading supported coverage/status."
            if material
            else "Thermodynamic gating does not materially improve benchmark error under the current threshold and remains diagnostic-only."
        ),
    )


def audit_all_thermodynamic_gating(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
    thresholds: BenchmarkThresholds = DEFAULT_BENCHMARK_THRESHOLDS,
) -> List[ThermodynamicGatingAudit]:
    bench_files = list(benchmark_files) if benchmark_files is not None else get_benchmark_files()
    return [
        audit_thermodynamic_gating(bf, target_tag=target_tag, thresholds=thresholds)
        for bf in bench_files
    ]


# ---------------------------------------------------------------------------
# Snapshot and index
# ---------------------------------------------------------------------------

def snapshot_benchmark_targets(
    bench_file: Path | str,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> List[BenchmarkTargetSnapshot]:
    from src.benchmark_registry import is_matrix_only_benchmark
    bench_path = Path(bench_file)
    bench = load_benchmark(bench_path)
    if is_matrix_only_benchmark(bench):
        return []
    from src.benchmark_evaluator import _is_supported_formulation, benchmark_to_formulation
    formulation = benchmark_to_formulation(bench)
    supported, _reason = _is_supported_formulation(formulation)
    if not supported:
        return []
    rec_result = _run_benchmark_recommendation(bench, target_tag=target_tag)

    snapshots: List[BenchmarkTargetSnapshot] = []
    for target in sorted(
        rec_result.get("targets", []),
        key=lambda row: (row.get("type", ""), -float(row.get("concentration", 0.0)), row.get("name", "")),
    ):
        snapshots.append(
            BenchmarkTargetSnapshot(
                benchmark_id=bench["benchmark_id"],
                bench_file=bench_path,
                target_name=str(target.get("name", "")),
                target_type=str(target.get("type", "unknown")),
                roles=list(target.get("roles", [target.get("type", "unknown")])),
                target_class=str(target.get("projection", {}).get("target_class", "unknown")),
                evidence_state=str(target.get("projection", {}).get("evidence_state", "still_missing")),
                predicted_ppb=float(target.get("concentration", 0.0)),
                proxy_ppb=float(target.get("proxy_concentration", target.get("concentration", 0.0))),
                observable_ratio=float(target.get("projection", {}).get("proxy_to_observable_ratio", 1.0)),
                weighted_flux=float(target.get("weighted_flux", 0.0)),
                span=float(target.get("span", math.inf)),
                depth=int(target.get("depth", 0)),
                volatile_class=str(target.get("projection", {}).get("volatile_class", "other")),
                matrix_factor=float(target.get("projection", {}).get("matrix_factor", 1.0)),
                headspace_factor=float(target.get("projection", {}).get("headspace_factor", 1.0)),
                headspace_observable=bool(target.get("headspace_observable", True)),
                headspace_class=str(target.get("headspace_class", "assumed_observable")),
                henry_kaw_25c=(
                    float(target["henry_kaw_25c"]) if target.get("henry_kaw_25c") is not None else None
                ),
                henry_source_name=target.get("henry_source_name"),
            )
        )
    return snapshots


def snapshot_all_benchmark_targets(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> List[BenchmarkTargetSnapshot]:
    bench_files = list(benchmark_files) if benchmark_files is not None else get_benchmark_files()
    snapshots: List[BenchmarkTargetSnapshot] = []
    for bench_file in bench_files:
        snapshots.extend(snapshot_benchmark_targets(bench_file, target_tag=target_tag))
    return snapshots


def build_benchmark_index(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> List[BenchmarkIndexEntry]:
    summaries = summarize_benchmarks(benchmark_files=benchmark_files, target_tag=target_tag)
    return [
        BenchmarkIndexEntry(
            benchmark_id=summary.benchmark_id,
            bench_file=summary.bench_file,
            tier=summary.tier,
            family=summary.family,
            protein_type=summary.protein_type,
            execution_path=summary.execution_path,
            benchmark_engine=summary.benchmark_engine,
            cantera_role=summary.cantera_role,
            thermodynamic_gating_policy=summary.thermodynamic_gating_policy,
            supported=summary.supported,
            reason=summary.reason,
            status=summary.overall_status,
            strict_ready=summary.strict_ready,
            process_state=summary.process_state,
            ranking_contract_status=summary.ranking_contract_status,
        )
        for summary in summaries
    ]


# ---------------------------------------------------------------------------
# Matrix promotion family status
# ---------------------------------------------------------------------------

def build_matrix_promotion_family_status(
    benchmark_files: Optional[Iterable[Path | str]] = None,
) -> List[MatrixPromotionFamilyStatus]:
    evidence_rows = build_matrix_benchmark_evidence_audit(benchmark_files)
    proteins = sorted({row.protein_type for row in evidence_rows})
    rows: List[MatrixPromotionFamilyStatus] = []
    for protein_type in proteins:
        subset = [row for row in evidence_rows if row.protein_type == protein_type]
        off_flavour_anchor_count = sum(
            1 for row in subset if row.target_profile == "adverse_only" and row.external_data_status == "external_quantitative"
        )
        meaty_candidate_count = sum(
            1 for row in subset if row.target_profile in {"meaty_positive", "mixed"}
        )
        external_meaty_anchor_count = sum(
            1
            for row in subset
            if row.target_profile in {"meaty_positive", "mixed"}
            and row.external_data_status == "external_quantitative"
            and row.reference_signal_origin == "measured_volatiles"
        )
        candidate_set_ready = off_flavour_anchor_count > 0 and meaty_candidate_count > 0
        external_assessment_unlocked = off_flavour_anchor_count > 0 and external_meaty_anchor_count > 0
        if off_flavour_anchor_count == 0:
            blocker = "missing external off-flavour anchor"
        elif meaty_candidate_count == 0:
            blocker = "missing meaty-positive benchmark candidate"
        elif external_meaty_anchor_count == 0:
            blocker = "missing external meaty-positive benchmark"
        else:
            blocker = "none"
        rows.append(
            MatrixPromotionFamilyStatus(
                protein_type=protein_type,
                off_flavour_anchor_count=off_flavour_anchor_count,
                meaty_candidate_count=meaty_candidate_count,
                external_meaty_anchor_count=external_meaty_anchor_count,
                candidate_set_ready=candidate_set_ready,
                external_assessment_unlocked=external_assessment_unlocked,
                blocker=blocker,
            )
        )
    return rows


# ---------------------------------------------------------------------------
# Matrix compound support status + target status artifact
# ---------------------------------------------------------------------------

def _matrix_compound_support_status(
    *,
    evidence_state: str,
    calibration_evidence_strength: str,
    reference_signal_origin: str,
    source_origin: str,
) -> str:
    evidence = str(evidence_state).strip().lower()
    strength = str(calibration_evidence_strength).strip().lower()
    signal_origin = str(reference_signal_origin).strip().lower()
    origin = str(source_origin).strip().lower()

    if evidence == "externally_benchmarked" and signal_origin == "measured_volatiles" and origin.startswith("external"):
        return "quantitative_closed"
    if evidence in {"internally_benchmarked", "conditional_calibration"} or signal_origin == "reference_volatiles":
        return "internal_candidate"
    if evidence in {"transferred_prior", "safety_reference"} or strength in {"class_anchored", "directional_transferred"}:
        return "directional_support"
    return "open_gap"


def _matrix_observable_or_calibration_priority(
    *,
    target_profile: str,
    ranking_contract_status: str,
    reference_signal_origin: str,
    external_data_status: str,
    support_counts: Mapping[str, int],
    compounds: Iterable[Mapping[str, Any]],
) -> bool:
    if target_profile not in {"mixed", "meaty_positive"}:
        return False
    if ranking_contract_status != "pass":
        return True
    if reference_signal_origin != "measured_volatiles":
        return True
    if external_data_status != "external_quantitative":
        return True
    if int(support_counts.get("quantitative_closed", 0)) < 2:
        return True
    if int(support_counts.get("open_gap", 0)) > 0:
        return True
    for compound in compounds:
        strength = str(compound.get("calibration_evidence_strength", "heuristic")).strip().lower()
        if strength in {"process_state_mismatch", "heuristic", "conditional_literature_anchored"}:
            return True
    return False


def _matrix_promotion_blocker(
    *,
    target_profile: str,
    ranking_contract_status: str,
    reference_signal_origin: str,
    external_data_status: str,
    support_counts: Mapping[str, int],
    compounds: Iterable[Mapping[str, Any]],
) -> str:
    if target_profile not in {"mixed", "meaty_positive"}:
        return "benchmark lacks meaty-positive targets"
    if ranking_contract_status != "pass":
        return "ranking contract not yet passing"
    if reference_signal_origin != "measured_volatiles":
        return "comparator is not measured_volatiles"
    if external_data_status != "external_quantitative":
        return "source is not external quantitative evidence"
    if int(support_counts.get("quantitative_closed", 0)) < 2:
        return "insufficient externally measured target closure"
    if int(support_counts.get("open_gap", 0)) > 0:
        return "named targets remain unresolved"
    for compound in compounds:
        strength = str(compound.get("calibration_evidence_strength", "heuristic")).strip().lower()
        if strength == "process_state_mismatch":
            return "observable calibration still mismatched to process state"
        if strength in {"heuristic", "conditional_literature_anchored"}:
            return "observable calibration still depends on non-closed evidence"
    if int(support_counts.get("internal_candidate", 0)) > 0 or int(support_counts.get("directional_support", 0)) > 0:
        return "depends on internal or transferred support"
    return "none"


@lru_cache(maxsize=8)
def _build_matrix_target_status_artifact_cached(
    benchmark_file_keys: tuple[str, ...],
    target_tag: str,
) -> Dict[str, Any]:
    bench_files = [Path(file_path) for file_path in benchmark_file_keys]
    benchmark_rows: List[Dict[str, Any]] = []
    support_totals = {
        "quantitative_closed": 0,
        "internal_candidate": 0,
        "directional_support": 0,
        "open_gap": 0,
    }

    for bench_file in bench_files:
        bench = load_benchmark(bench_file)
        metadata = get_benchmark_metadata(bench)
        if metadata.execution_path not in {"matrix_only", "matrix_precursor_augmented"}:
            continue

        evaluation = evaluate_benchmark(bench_file, target_tag=target_tag)
        summary = summarize_evaluation(evaluation, protein_type=str(bench.get("protein_type", "free")))
        evidence = assess_matrix_benchmark_evidence(bench_file)
        contract = get_matrix_ranking_contract(bench)
        adverse_markers = {str(item).strip().lower() for item in contract.get("adverse_markers", [])}
        compounds: List[Dict[str, Any]] = []
        benchmark_counts = {
            "quantitative_closed": 0,
            "internal_candidate": 0,
            "directional_support": 0,
            "open_gap": 0,
        }

        for item in contract.get("observable_targets", []):
            compound_name = str(item.get("name", "")).strip()
            if not compound_name:
                continue
            meta = _projection_metadata_for_match(
                evaluation,
                CompoundComparison(
                    compound=compound_name,
                    measured_ppb=0.0,
                    predicted_ppb=0.0,
                    matched_name=None,
                    uncertainty_pct=None,
                ),
            )
            role = str(item.get("role", "adverse_marker" if compound_name.lower() in adverse_markers else "desirable_marker"))
            support_status = _matrix_compound_support_status(
                evidence_state=str(meta.get("evidence_state", item.get("evidence_state", "still_missing"))),
                calibration_evidence_strength=str(meta.get("calibration_evidence_strength", "heuristic")),
                reference_signal_origin=summary.reference_signal_origin,
                source_origin=evidence.source_origin,
            )
            benchmark_counts[support_status] += 1
            support_totals[support_status] += 1
            compounds.append(
                {
                    "compound": compound_name,
                    "role": role,
                    "target_class": str(meta.get("target_class", item.get("target_class", "unknown"))),
                    "evidence_state": str(meta.get("evidence_state", item.get("evidence_state", "still_missing"))),
                    "calibration_source": str(meta.get("calibration_source", "unknown")),
                    "calibration_evidence_strength": str(meta.get("calibration_evidence_strength", "heuristic")),
                    "support_status": support_status,
                }
            )

        promotion_blocker = _matrix_promotion_blocker(
            target_profile=evidence.target_profile,
            ranking_contract_status=summary.ranking_contract_status,
            reference_signal_origin=summary.reference_signal_origin,
            external_data_status=evidence.external_data_status,
            support_counts=benchmark_counts,
            compounds=compounds,
        )
        external_decision_ready = (
            evidence.target_profile in {"mixed", "meaty_positive"}
            and summary.ranking_contract_status == "pass"
            and summary.reference_signal_origin == "measured_volatiles"
            and evidence.external_data_status == "external_quantitative"
            and benchmark_counts["quantitative_closed"] >= 2
            and benchmark_counts["internal_candidate"] == 0
            and benchmark_counts["directional_support"] == 0
            and benchmark_counts["open_gap"] == 0
        )
        evidence_or_calibration_priority_ready = _matrix_observable_or_calibration_priority(
            target_profile=evidence.target_profile,
            ranking_contract_status=summary.ranking_contract_status,
            reference_signal_origin=summary.reference_signal_origin,
            external_data_status=evidence.external_data_status,
            support_counts=benchmark_counts,
            compounds=compounds,
        )
        mechanistic_priority_ready = (
            evidence.target_profile in {"mixed", "meaty_positive"}
            and summary.ranking_contract_status == "pass"
            and not external_decision_ready
            and not evidence_or_calibration_priority_ready
            and (benchmark_counts["internal_candidate"] + benchmark_counts["directional_support"]) >= 1
        )

        if external_decision_ready:
            blocker_class = "externally_decision_ready"
            promotion_claim_posture = "external_decision_ready"
            next_best_action = "use_for_external_decision"
            best_computational_action = "none"
        elif mechanistic_priority_ready:
            blocker_class = "mechanistic_blocker"
            promotion_claim_posture = "mechanistic_triage_lane"
            next_best_action = "prioritize_mechanistic_refinement"
            best_computational_action = "named_mechanistic_refinement_only"
        elif evidence.target_profile in {"mixed", "meaty_positive"}:
            blocker_class = "observable_or_calibration_blocker"
            promotion_claim_posture = (
                "observable_or_calibration_hold"
                if summary.reference_signal_origin == "measured_volatiles" and evidence.external_data_status == "external_quantitative"
                else "not_a_promotion_lane"
            )
            next_best_action = "improve_observable_or_calibration"
            best_computational_action = "improve_observable_or_calibration_before_qm"
        else:
            blocker_class = "not_a_promotion_lane"
            promotion_claim_posture = "not_a_promotion_lane"
            next_best_action = "retain_as_adverse_anchor"
            best_computational_action = "none"

        benchmark_rows.append(
            {
                "benchmark_id": summary.benchmark_id,
                "bench_file": str(summary.bench_file),
                "protein_type": summary.protein_type,
                "execution_path": summary.execution_path,
                "process_state": summary.process_state,
                "target_profile": evidence.target_profile,
                "reference_signal_origin": summary.reference_signal_origin,
                "source_origin": evidence.source_origin,
                "ranking_contract_status": summary.ranking_contract_status,
                "support_counts": benchmark_counts,
                "quantitative_support_ready": benchmark_counts["quantitative_closed"] > 0,
                "promotion_ready": external_decision_ready,
                "evidence_or_calibration_priority_ready": evidence_or_calibration_priority_ready,
                "mechanistic_priority_ready": mechanistic_priority_ready,
                "promotion_blocker": promotion_blocker,
                "blocker_class": blocker_class,
                "promotion_claim_posture": promotion_claim_posture,
                "next_best_action": next_best_action,
                "best_computational_action": best_computational_action,
                "compounds": compounds,
            }
        )

    return {
        "schema_version": "1.0",
        "description": "Matrix target support status artifact distinguishing external decision-ready support from mechanistic-priority candidates and unresolved external gaps.",
        "benchmarks": benchmark_rows,
        "summary": {
            "total_benchmarks": len(benchmark_rows),
            "quantitative_support_ready": sum(1 for row in benchmark_rows if row["quantitative_support_ready"]),
            "promotion_ready": sum(1 for row in benchmark_rows if row["promotion_ready"]),
            "evidence_or_calibration_priority_ready": sum(
                1 for row in benchmark_rows if row["evidence_or_calibration_priority_ready"]
            ),
            "mechanistic_priority_ready": sum(1 for row in benchmark_rows if row["mechanistic_priority_ready"]),
            **support_totals,
        },
    }


def build_matrix_target_status_artifact(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> Dict[str, Any]:
    return deepcopy(_build_matrix_target_status_artifact_cached(_benchmark_cache_key(benchmark_files), target_tag))


# ---------------------------------------------------------------------------
# Matrix promotion contract artifact
# ---------------------------------------------------------------------------

def _matrix_promotion_requirement_rows(
    benchmark_row: Mapping[str, Any],
    evidence_row: MatrixBenchmarkEvidence,
) -> List[Dict[str, Any]]:
    support_counts = benchmark_row.get("support_counts", {})
    return [
        {
            "key": "meaty_positive_targets_present",
            "label": "Target profile includes meaty-positive compounds",
            "passed": benchmark_row.get("target_profile") in {"mixed", "meaty_positive"},
            "detail": str(benchmark_row.get("target_profile", "unknown")),
        },
        {
            "key": "ranking_contract_passes",
            "label": "Ranking contract passes",
            "passed": benchmark_row.get("ranking_contract_status") == "pass",
            "detail": str(benchmark_row.get("ranking_contract_status", "unknown")),
        },
        {
            "key": "comparator_is_measured_volatiles",
            "label": "Comparator signal is wet-lab measured_volatiles",
            "passed": benchmark_row.get("reference_signal_origin") == "measured_volatiles",
            "detail": str(benchmark_row.get("reference_signal_origin", "unknown")),
        },
        {
            "key": "external_quantitative_origin",
            "label": "Source is externally quantitative",
            "passed": evidence_row.external_data_status == "external_quantitative",
            "detail": evidence_row.external_data_status,
        },
        {
            "key": "minimum_quantitative_closed_targets",
            "label": "At least two compounds are quantitatively closed",
            "passed": int(support_counts.get("quantitative_closed", 0)) >= 2,
            "detail": str(int(support_counts.get("quantitative_closed", 0))),
        },
        {
            "key": "no_internal_or_directional_dependencies",
            "label": "No internal-candidate or directional dependencies remain",
            "passed": int(support_counts.get("internal_candidate", 0)) == 0 and int(support_counts.get("directional_support", 0)) == 0,
            "detail": f"internal={int(support_counts.get('internal_candidate', 0))}; directional={int(support_counts.get('directional_support', 0))}",
        },
    ]


def _select_matrix_promotion_target(
    benchmark_rows: Iterable[Mapping[str, Any]],
    evidence_rows: Iterable[MatrixBenchmarkEvidence],
) -> Optional[Dict[str, Any]]:
    benchmark_list = list(benchmark_rows)
    evidence_list = list(evidence_rows)
    candidates = [
        row for row in benchmark_list
        if row.get("target_profile") in {"mixed", "meaty_positive"}
        and row.get("ranking_contract_status") == "pass"
    ]
    if not candidates:
        return None

    external_anchor_counts: Dict[str, int] = defaultdict(int)
    distinct_external_states: Dict[str, set] = defaultdict(set)
    for row in evidence_list:
        if row.external_data_status != "external_quantitative":
            continue
        external_anchor_counts[row.protein_type] += 1
        if row.process_state:
            distinct_external_states[row.protein_type].add(str(row.process_state))

    def rank_tuple(row: Mapping[str, Any]) -> tuple:
        protein_type = str(row.get("protein_type", "free"))
        counts = row.get("support_counts", {})
        return (
            int(external_anchor_counts.get(protein_type, 0)),
            len(distinct_external_states.get(protein_type, set())),
            int(counts.get("quantitative_closed", 0)),
            -1 * (int(counts.get("internal_candidate", 0)) + int(counts.get("directional_support", 0)) + int(counts.get("open_gap", 0))),
            "1" if row.get("target_profile") == "mixed" else "0",
            str(row.get("benchmark_id", "unknown")),
        )

    selected = sorted(candidates, key=rank_tuple, reverse=True)[0]
    protein_type = str(selected.get("protein_type", "free"))
    rationale = [
        f"same_protein_external_anchor_count={int(external_anchor_counts.get(protein_type, 0))}",
        f"same_protein_external_process_states={len(distinct_external_states.get(protein_type, set()))}",
        f"quantitative_closed={int(selected.get('support_counts', {}).get('quantitative_closed', 0))}",
        f"internal_candidate={int(selected.get('support_counts', {}).get('internal_candidate', 0))}",
    ]
    return {
        "benchmark_id": selected.get("benchmark_id"),
        "protein_type": protein_type,
        "process_state": selected.get("process_state"),
        "target_profile": selected.get("target_profile"),
        "rationale": rationale,
        "selection_policy": "prefer_mixed_matrix_lanes_with_the_broadest_same_protein_external_anchor_span_before_mechanistic_escalation",
    }


def build_matrix_promotion_contract_artifact(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> Dict[str, Any]:
    status_payload = build_matrix_target_status_artifact(benchmark_files, target_tag=target_tag)
    evidence_rows = build_matrix_benchmark_evidence_audit(benchmark_files)
    evidence_lookup = {row.benchmark_id: row for row in evidence_rows}

    benchmark_assessments: List[Dict[str, Any]] = []
    for row in status_payload.get("benchmarks", []):
        benchmark_id = str(row.get("benchmark_id", "unknown"))
        evidence = evidence_lookup.get(benchmark_id)
        if evidence is None:
            continue
        requirements = _matrix_promotion_requirement_rows(row, evidence)
        benchmark_assessments.append(
            {
                "benchmark_id": benchmark_id,
                "protein_type": row.get("protein_type"),
                "process_state": row.get("process_state"),
                "target_profile": row.get("target_profile"),
                "promotion_ready": bool(row.get("promotion_ready", False)),
                "promotion_blocker": row.get("promotion_blocker", "unknown"),
                "requirements": requirements,
            }
        )

    selected_target = _select_matrix_promotion_target(status_payload.get("benchmarks", []), evidence_rows)
    summary = status_payload.get("summary", {})
    return {
        "schema_version": "1.0",
        "description": "Explicit matrix promotion contract defining how a matrix benchmark moves from internal candidate support to external decision readiness.",
        "promotion_rule": {
            "contract_id": "matrix_external_decision_ready_v1",
            "minimum_quantitative_closed_targets": 2,
            "disallow_internal_candidate_support": True,
            "disallow_directional_support": True,
            "requires_measured_volatiles": True,
            "requires_external_quantitative_origin": True,
            "requires_mixed_or_meaty_positive_target_profile": True,
            "requires_passing_ranking_contract": True,
            "notes": [
                "External decision readiness is a benchmark-level promotion state, not a generic matrix-family claim.",
                "Internal reproducibility candidates and transferred priors can strengthen triage, but they do not by themselves unlock promotion.",
                "Mechanistic refinement stays secondary until the observable audit says the remaining blocker is no longer external evidence or transfer dependence.",
            ],
        },
        "selected_promotion_target": selected_target,
        "benchmarks": benchmark_assessments,
        "summary": {
            "benchmarks_assessed": len(benchmark_assessments),
            "promotion_ready": int(summary.get("promotion_ready", 0)),
            "mechanistic_priority_ready": int(summary.get("mechanistic_priority_ready", 0)),
        },
    }


# ---------------------------------------------------------------------------
# Matrix observable closure audit
# ---------------------------------------------------------------------------

def _matrix_closure_action(
    *,
    compound_row: Mapping[str, Any],
    benchmark_row: Mapping[str, Any],
) -> str:
    support_status = str(compound_row.get("support_status", "open_gap"))
    calibration_strength = str(compound_row.get("calibration_evidence_strength", "heuristic"))

    if support_status == "quantitative_closed":
        return "already_closed"
    if support_status in {"internal_candidate", "directional_support"} and (
        calibration_strength in {"class_anchored", "directional_transferred"}
    ):
        return "class_level_transfer_acceptable"
    if calibration_strength in {"literature_anchored", "conditional_literature_anchored"}:
        return "literature_anchor_available"
    if calibration_strength in {"process_state_mismatch", "heuristic"}:
        return "evidence_or_calibration_blocker"
    if benchmark_row.get("mechanistic_priority_ready") and support_status == "internal_candidate":
        return "mechanistic_blocker"
    if benchmark_row.get("evidence_or_calibration_priority_ready"):
        return "evidence_or_calibration_blocker"
    return "external_data_blocker"


def _mechanistic_refinement_expected_change(compound_rows: Iterable[Mapping[str, Any]]) -> str:
    rows = list(compound_rows)
    roles = {str(row.get("role", "unknown")) for row in rows}
    if roles == {"adverse_marker"}:
        return "clarify whether the lane remains meaty-positive once named adverse-marker closure is resolved"
    if roles == {"desirable_marker"}:
        return "clarify whether named desirable-marker closure can move the lane toward external decision readiness"
    return "clarify whether named mechanistic blockers materially change benchmark readiness before broader retuning"


@lru_cache(maxsize=8)
def _build_matrix_observable_closure_audit_cached(
    benchmark_file_keys: tuple[str, ...],
    target_tag: str,
) -> Dict[str, Any]:
    benchmark_files = [Path(file_path) for file_path in benchmark_file_keys]
    status_payload = _build_matrix_target_status_artifact_cached(benchmark_file_keys, target_tag)
    evidence_rows = build_matrix_benchmark_evidence_audit(benchmark_files)
    evidence_lookup = {row.benchmark_id: row for row in evidence_rows}
    selected_target = _select_matrix_promotion_target(status_payload.get("benchmarks", []), evidence_rows)

    audit_rows: List[Dict[str, Any]] = []
    action_counts: Dict[str, int] = defaultdict(int)
    mechanistic_watchlist: List[Dict[str, Any]] = []
    for benchmark_row in status_payload.get("benchmarks", []):
        if benchmark_row.get("target_profile") not in {"mixed", "meaty_positive"}:
            continue
        benchmark_id = str(benchmark_row.get("benchmark_id", "unknown"))
        evidence = evidence_lookup.get(benchmark_id)
        compound_rows: List[Dict[str, Any]] = []
        benchmark_action_counts: Dict[str, int] = defaultdict(int)
        for compound_row in benchmark_row.get("compounds", []):
            closure_action = _matrix_closure_action(compound_row=compound_row, benchmark_row=benchmark_row)
            benchmark_action_counts[closure_action] += 1
            action_counts[closure_action] += 1
            compound_rows.append({**compound_row, "closure_action": closure_action})
        audit_rows.append(
            {
                "benchmark_id": benchmark_id,
                "protein_type": benchmark_row.get("protein_type"),
                "process_state": benchmark_row.get("process_state"),
                "target_profile": benchmark_row.get("target_profile"),
                "promotion_blocker": benchmark_row.get("promotion_blocker"),
                "source_origin": evidence.source_origin if evidence is not None else "unknown",
                "compounds": compound_rows,
                "closure_action_counts": dict(sorted(benchmark_action_counts.items())),
            }
        )

        if benchmark_row.get("mechanistic_priority_ready"):
            blocker_rows = [
                row for row in compound_rows
                if str(row.get("closure_action", "unknown")) == "mechanistic_blocker"
            ]
            if blocker_rows:
                mechanistic_watchlist.append(
                    {
                        "benchmark_id": benchmark_id,
                        "protein_type": benchmark_row.get("protein_type"),
                        "process_state": benchmark_row.get("process_state"),
                        "promotion_blocker": benchmark_row.get("promotion_blocker"),
                        "target_compounds": [str(row.get("compound", "unknown")) for row in blocker_rows],
                        "target_roles": sorted({str(row.get("role", "unknown")) for row in blocker_rows}),
                        "expected_decision_change": _mechanistic_refinement_expected_change(blocker_rows),
                        "allowed_scope": "named_compound_refinement_only",
                        "offline_compute_gate": "escalate_only_if_named_compounds_change_benchmark_visible_decision_readiness",
                    }
                )

    return {
        "schema_version": "1.0",
        "description": "Compound-level observable closure audit for mixed matrix lanes, labeling the next closure action needed for each decision-driving compound.",
        "selected_promotion_target": selected_target,
        "benchmarks": audit_rows,
        "mechanistic_refinement_watchlist": mechanistic_watchlist,
        "summary": {
            "benchmarks_audited": len(audit_rows),
            "mechanistic_watchlist_count": len(mechanistic_watchlist),
            "closure_action_counts": dict(sorted(action_counts.items())),
        },
    }


def build_matrix_observable_closure_audit(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> Dict[str, Any]:
    return deepcopy(_build_matrix_observable_closure_audit_cached(_benchmark_cache_key(benchmark_files), target_tag))
