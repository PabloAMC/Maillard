from __future__ import annotations

import json
import math
import os
import statistics
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple

from src import data_paths
from src.benchmark_validation import evaluate_benchmark, get_benchmark_files, get_benchmark_metadata, load_benchmark
from src.cross_validation import compute_leverage
from src.uncertainty_propagation import propagate_benchmarks


_RUNTIME_MULTIPLIER_ENV = "MAILLARD_MATRIX_CALIBRATION_MULTIPLIERS"
_MATRIX_EXECUTION_PATHS = ("matrix_only", "matrix_precursor_augmented")
_TRUST_LOOP_EXECUTION_PATHS = ("free_precursor", "matrix_precursor_augmented")
_DEFAULT_GRID = (0.55, 0.7, 0.85, 1.0, 1.15, 1.35, 1.6, 1.9)


def _matrix_benchmark_groups() -> Dict[str, List[Path]]:
    grouped: Dict[str, List[Path]] = {}
    for bench_file in get_benchmark_files():
        bench = load_benchmark(bench_file)
        metadata = get_benchmark_metadata(bench)
        if metadata.execution_path not in _MATRIX_EXECUTION_PATHS:
            continue
        protein_type = str(bench.get("protein_type", "")).strip()
        if not protein_type or protein_type == "free":
            continue
        grouped.setdefault(protein_type, []).append(Path(bench_file))
    return grouped


@contextmanager
def matrix_observable_multiplier_overrides(overrides: Mapping[str, float]) -> Iterator[None]:
    previous = os.environ.get(_RUNTIME_MULTIPLIER_ENV)
    try:
        os.environ[_RUNTIME_MULTIPLIER_ENV] = json.dumps({str(key): float(value) for key, value in overrides.items()})
        yield
    finally:
        if previous is None:
            os.environ.pop(_RUNTIME_MULTIPLIER_ENV, None)
        else:
            os.environ[_RUNTIME_MULTIPLIER_ENV] = previous


def _collect_abs_log10_errors(benchmark_files: Sequence[Path]) -> List[float]:
    errors: List[float] = []
    for bench_file in benchmark_files:
        evaluation = evaluate_benchmark(bench_file)
        for comparison in evaluation.comparisons:
            if comparison.matched_name is None:
                continue
            measured = float(comparison.measured_ppb)
            predicted = float(comparison.predicted_ppb)
            if measured <= 0.0 or predicted <= 0.0:
                continue
            errors.append(abs(math.log10(predicted / measured)))
    return errors


def _mean_error_for_multiplier(protein_type: str, benchmark_files: Sequence[Path], multiplier: float) -> Optional[float]:
    with matrix_observable_multiplier_overrides({protein_type: multiplier}):
        errors = _collect_abs_log10_errors(benchmark_files)
    if not errors:
        return None
    return statistics.fmean(errors)


def select_best_multiplier(
    protein_type: str,
    benchmark_files: Sequence[Path],
    *,
    grid: Sequence[float] = _DEFAULT_GRID,
) -> Dict[str, Any]:
    candidates: List[Dict[str, Any]] = []
    for multiplier in grid:
        mean_error = _mean_error_for_multiplier(protein_type, benchmark_files, float(multiplier))
        candidates.append(
            {
                "protein_type": protein_type,
                "multiplier": float(multiplier),
                "mean_abs_log10_error": mean_error,
            }
        )
    usable = [row for row in candidates if row["mean_abs_log10_error"] is not None]
    if not usable:
        return {
            "protein_type": protein_type,
            "selected_multiplier": 1.0,
            "baseline_mean_abs_log10_error": None,
            "selected_mean_abs_log10_error": None,
            "grid": candidates,
        }

    selected = min(usable, key=lambda row: (float(row["mean_abs_log10_error"]), abs(float(row["multiplier"]) - 1.0)))
    baseline = next(row for row in usable if math.isclose(float(row["multiplier"]), 1.0, rel_tol=1e-9, abs_tol=1e-9))
    baseline_error = float(baseline["mean_abs_log10_error"])
    selected_error = float(selected["mean_abs_log10_error"])
    improvement_pct = ((baseline_error - selected_error) / baseline_error * 100.0) if baseline_error > 0.0 else 0.0
    return {
        "protein_type": protein_type,
        "selected_multiplier": float(selected["multiplier"]),
        "baseline_mean_abs_log10_error": baseline_error,
        "selected_mean_abs_log10_error": selected_error,
        "improvement_pct": improvement_pct,
        "grid": candidates,
    }


def _median_benchmark_residual(leverage_report: Mapping[str, Any]) -> Optional[float]:
    residuals = [
        float(row.get("mean_abs_log10_residual", 0.0))
        for row in leverage_report.get("benchmarks", []) or []
        if row.get("mean_abs_log10_residual") is not None
    ]
    if not residuals:
        return None
    return statistics.median(residuals)


def _inside_ci_lookup(payload: Mapping[str, Any]) -> Dict[Tuple[str, str], bool]:
    lookup: Dict[Tuple[str, str], bool] = {}
    for benchmark in payload.get("benchmarks", []) or []:
        benchmark_id = str(benchmark.get("benchmark_id", ""))
        for compound in benchmark.get("compounds", []) or []:
            compound_name = str(compound.get("compound", ""))
            lookup[(benchmark_id, compound_name)] = bool(compound.get("inside_ci"))
    return lookup


def _count_ci_regressions(baseline_payload: Mapping[str, Any], candidate_payload: Mapping[str, Any]) -> int:
    baseline_lookup = _inside_ci_lookup(baseline_payload)
    candidate_lookup = _inside_ci_lookup(candidate_payload)
    regressions = 0
    for key, baseline_inside in baseline_lookup.items():
        if baseline_inside and not candidate_lookup.get(key, False):
            regressions += 1
    return regressions


def assess_refit_candidate(
    baseline_uncertainty: Mapping[str, Any],
    candidate_uncertainty: Mapping[str, Any],
    baseline_trust: Mapping[str, Any],
    candidate_trust: Mapping[str, Any],
) -> Dict[str, Any]:
    baseline_loo = compute_leverage(baseline_uncertainty)
    candidate_loo = compute_leverage(candidate_uncertainty)
    baseline_median = _median_benchmark_residual(baseline_loo)
    candidate_median = _median_benchmark_residual(candidate_loo)
    ci_regressions = _count_ci_regressions(baseline_uncertainty, candidate_uncertainty)
    baseline_trust_rate = baseline_trust.get("summary", {}).get("ci_coverage_rate")
    candidate_trust_rate = candidate_trust.get("summary", {}).get("ci_coverage_rate")

    accepted = (
        baseline_median is not None
        and candidate_median is not None
        and candidate_median < baseline_median
        and ci_regressions == 0
        and (baseline_trust_rate is None or candidate_trust_rate is None or float(candidate_trust_rate) >= float(baseline_trust_rate))
    )
    return {
        "accepted": accepted,
        "baseline_median_abs_log10_residual": baseline_median,
        "candidate_median_abs_log10_residual": candidate_median,
        "ci_regressions": ci_regressions,
        "baseline_trust_rate": baseline_trust_rate,
        "candidate_trust_rate": candidate_trust_rate,
        "baseline_loo": baseline_loo,
        "candidate_loo": candidate_loo,
    }


def _write_multiplier_file(overrides: Mapping[str, float], *, decision_summary: Mapping[str, Any]) -> Path:
    payload = {
        "version": 1,
        "entries": {
            protein_type: {
                "observable_factor_multiplier": float(multiplier),
                "source": "matrix_recalibration_refit",
                "provenance": str(decision_summary.get("decision_id", "matrix_recalibration_refit")),
            }
            for protein_type, multiplier in overrides.items()
        },
    }
    data_paths.MATRIX_CALIBRATION_OFFSETS.parent.mkdir(parents=True, exist_ok=True)
    data_paths.MATRIX_CALIBRATION_OFFSETS.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return data_paths.MATRIX_CALIBRATION_OFFSETS


def render_refit_decision_markdown(report: Mapping[str, Any]) -> str:
    lines = [
        "# Matrix Calibration Refit Decision",
        "",
        f"- Decision: {report.get('decision', 'unknown')}",
        f"- Decision id: {report.get('decision_id', 'unknown')}",
        f"- Baseline matrix-panel median |log10 error|: {float(report.get('baseline_median_abs_log10_residual', 0.0) or 0.0):.3f}",
        f"- Candidate matrix-panel median |log10 error|: {float(report.get('candidate_median_abs_log10_residual', 0.0) or 0.0):.3f}",
        f"- CI regressions: {int(report.get('ci_regressions', 0) or 0)}",
        "",
        "## Protein Multipliers",
        "",
        "| Protein Type | Selected Multiplier | Baseline Mean | Candidate Mean | Improvement % |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for row in report.get("protein_sweeps", []) or []:
        lines.append(
            f"| {row.get('protein_type', 'unknown')} | {float(row.get('selected_multiplier', 1.0) or 1.0):.3f} | "
            f"{float(row.get('baseline_mean_abs_log10_error', 0.0) or 0.0):.3f} | "
            f"{float(row.get('selected_mean_abs_log10_error', 0.0) or 0.0):.3f} | "
            f"{float(row.get('improvement_pct', 0.0) or 0.0):+.1f} |"
        )
    lines.extend(
        [
            "",
            "## Trust Loop",
            "",
            f"- Before: {report.get('trust_loop_before', {}).get('ci_coverage_hits', 0)}/{report.get('trust_loop_before', {}).get('matched_compound_count', 0)}",
            f"- After: {report.get('trust_loop_after', {}).get('ci_coverage_hits', 0)}/{report.get('trust_loop_after', {}).get('matched_compound_count', 0)}",
        ]
    )
    if report.get("persistence_path"):
        lines.append(f"- Persistence path: {report.get('persistence_path')}")
    lines.append("")
    return "\n".join(lines)


def write_refit_artifacts(report: Mapping[str, Any], *, output_dir: Path | str = data_paths.VALIDATION_DIR) -> Dict[str, Path]:
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    json_path = destination / "calibration_refit_decision.json"
    md_path = destination / "calibration_refit_decision.md"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    md_path.write_text(render_refit_decision_markdown(report), encoding="utf-8")
    return {"json": json_path, "md": md_path}


def run_matrix_recalibration(
    *,
    write_changes: bool = False,
    uncertainty_samples: int = 60,
    uncertainty_seed: int = 0,
) -> Dict[str, Any]:
    grouped = _matrix_benchmark_groups()
    protein_sweeps = [select_best_multiplier(protein_type, bench_files) for protein_type, bench_files in sorted(grouped.items())]
    candidate_overrides = {
        row["protein_type"]: float(row["selected_multiplier"])
        for row in protein_sweeps
        if row.get("selected_multiplier") is not None and not math.isclose(float(row["selected_multiplier"]), 1.0, rel_tol=1e-9, abs_tol=1e-9)
    }

    baseline_uncertainty = propagate_benchmarks(
        execution_paths=_MATRIX_EXECUTION_PATHS,
        n_samples=uncertainty_samples,
        seed=uncertainty_seed,
    )
    baseline_trust = propagate_benchmarks(
        execution_paths=_TRUST_LOOP_EXECUTION_PATHS,
        n_samples=uncertainty_samples,
        seed=uncertainty_seed,
    )
    with matrix_observable_multiplier_overrides(candidate_overrides):
        candidate_uncertainty = propagate_benchmarks(
            execution_paths=_MATRIX_EXECUTION_PATHS,
            n_samples=uncertainty_samples,
            seed=uncertainty_seed,
        )
        candidate_trust = propagate_benchmarks(
            execution_paths=_TRUST_LOOP_EXECUTION_PATHS,
            n_samples=uncertainty_samples,
            seed=uncertainty_seed,
        )

    assessment = assess_refit_candidate(
        baseline_uncertainty,
        candidate_uncertainty,
        baseline_trust,
        candidate_trust,
    )
    report: Dict[str, Any] = {
        "decision_id": "refit_2026-05-10",
        "decision": "accepted" if assessment["accepted"] else "rejected",
        "protein_sweeps": protein_sweeps,
        "candidate_overrides": candidate_overrides,
        "uncertainty_samples": uncertainty_samples,
        "uncertainty_seed": uncertainty_seed,
        "baseline_median_abs_log10_residual": assessment["baseline_median_abs_log10_residual"],
        "candidate_median_abs_log10_residual": assessment["candidate_median_abs_log10_residual"],
        "ci_regressions": assessment["ci_regressions"],
        "matrix_panel_before": baseline_uncertainty.get("summary", {}),
        "matrix_panel_after": candidate_uncertainty.get("summary", {}),
        "trust_loop_before": baseline_trust.get("summary", {}),
        "trust_loop_after": candidate_trust.get("summary", {}),
    }

    persistence_path: Optional[Path] = None
    if assessment["accepted"] and write_changes and candidate_overrides:
        persistence_path = _write_multiplier_file(candidate_overrides, decision_summary=report)
    report["persistence_path"] = str(persistence_path) if persistence_path is not None else None
    return report