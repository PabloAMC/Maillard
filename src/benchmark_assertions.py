"""
src/benchmark_assertions.py — Structured matrix benchmark assertions.

Extracted from src/benchmark_validation.py (Priority 2 decomposition).
Contains the logic for building and rendering per-benchmark assertion rows
that test whether a matrix benchmark meets its ranking contract, coverage,
and scale requirements.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional

from src.benchmark_types import (
    BenchmarkThresholds,
    CompoundComparison,
    MatrixBenchmarkAssertion,
)
from src.benchmark_registry import (
    DEFAULT_TARGET_TAG,
    get_benchmark_files,
    get_benchmark_metadata,
    get_matrix_ranking_contract,
    load_benchmark,
)
from src.benchmark_evaluator import (
    evaluate_benchmark,
    _best_prediction_match,
)
from src.benchmark_reporting import (
    summarize_evaluation,
    assess_matrix_benchmark_evidence,
    _projection_metadata_for_match,
)
from src.validation_contract import DEFAULT_VALIDATION_CONTRACT


# ---------------------------------------------------------------------------
# Threshold helpers
# ---------------------------------------------------------------------------

def _matrix_assertion_thresholds(
    bench: dict,
    *,
    protein_type: str,
    thresholds: BenchmarkThresholds,
) -> Dict[str, float]:
    contract = bench.get("matrix_ranking_contract") or {}
    configured = (
        contract.get("assertion_thresholds")
        or contract.get("validation_contract", {}).get("scale_thresholds")
        or {}
    )
    observable_targets = get_matrix_ranking_contract(bench).get("observable_targets", [])
    return {
        "min_coverage": float(
            configured.get("min_coverage", thresholds.full_coverage_threshold)
        ),
        "top_k": float(configured.get("top_k", len(observable_targets))),
        "max_ratio": float(
            configured.get("max_ratio", thresholds.ratio_threshold_for(protein_type))
        ),
    }


def _predicted_order_lookup(evaluation) -> List[str]:
    from src.text_utils import normalize_compound_name as _normalize_name
    normalized_rows = []
    seen: set[str] = set()
    for name, value in sorted(
        evaluation.predicted_ppb.items(), key=lambda item: item[1], reverse=True
    ):
        normalized = _normalize_name(name)
        if normalized in seen:
            continue
        seen.add(normalized)
        normalized_rows.append(name)
    return normalized_rows


def _matched_contract_prediction_rows(
    evaluation, contract_names: Iterable[str]
) -> List[tuple[str, float]]:
    rows: List[tuple[str, float]] = []
    for contract_name in contract_names:
        matched_name, predicted_value, _score = _best_prediction_match(
            str(contract_name), evaluation.predicted_ppb
        )
        if matched_name is None:
            continue
        rows.append((str(contract_name), float(predicted_value)))
    return rows


# ---------------------------------------------------------------------------
# Assertion builder
# ---------------------------------------------------------------------------

def build_matrix_benchmark_assertions(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
    thresholds: BenchmarkThresholds = None,
) -> List[MatrixBenchmarkAssertion]:
    from src.validation_contract import DEFAULT_VALIDATION_CONTRACT
    if thresholds is None:
        thresholds = DEFAULT_VALIDATION_CONTRACT.thresholds

    bench_files = list(benchmark_files) if benchmark_files is not None else get_benchmark_files()
    rows: List[MatrixBenchmarkAssertion] = []
    for bench_file in bench_files:
        bench = load_benchmark(bench_file)
        metadata = get_benchmark_metadata(bench)
        if metadata.execution_path not in {"matrix_only", "matrix_precursor_augmented"}:
            continue

        evaluation = evaluate_benchmark(bench_file, target_tag=target_tag)
        summary = summarize_evaluation(
            evaluation,
            protein_type=bench.get("protein_type", "free"),
            thresholds=thresholds,
        )
        evidence = assess_matrix_benchmark_evidence(bench_file)
        assertion_thresholds = _matrix_assertion_thresholds(
            bench,
            protein_type=str(bench.get("protein_type", "free")),
            thresholds=thresholds,
        )
        top_k = max(0, int(assertion_thresholds["top_k"]))
        matched_ranked_targets = _matched_contract_prediction_rows(
            evaluation,
            summary.ranked_observable_targets,
        )
        ranked_targets = [name for name, _value in matched_ranked_targets[:top_k]]
        predicted_top_k = {
            name
            for name, _value in sorted(
                matched_ranked_targets, key=lambda row: row[1], reverse=True
            )[:top_k]
        }
        top_k_hits = sum(1 for name in ranked_targets if name in predicted_top_k)
        if top_k == 0:
            top_k_status = "n/a"
        elif top_k_hits == top_k and summary.ranking_contract_status == "pass":
            top_k_status = "pass"
        else:
            top_k_status = "fail"

        adverse_rows = _matched_contract_prediction_rows(
            evaluation, summary.adverse_markers
        )
        adverse_markers = list(summary.adverse_markers)
        if len(adverse_markers) <= 1:
            adverse_order_status = "n/a"
        else:
            predicted_adverse = [
                name
                for name, _value in sorted(adverse_rows, key=lambda row: row[1], reverse=True)
            ]
            if len(predicted_adverse) != len(adverse_markers):
                adverse_order_status = "missing"
            else:
                adverse_order_status = (
                    "pass" if predicted_adverse == adverse_markers else "fail"
                )

        min_coverage = float(assertion_thresholds["min_coverage"])
        ratio_tolerance = float(assertion_thresholds["max_ratio"])
        if summary.max_ratio is None:
            ratio_status = "fail"
        elif summary.max_ratio <= ratio_tolerance:
            ratio_status = "pass"
        else:
            ratio_status = "fail"

        overall_pass = (
            summary.coverage >= min_coverage
            and top_k_status in {"pass", "n/a"}
            and adverse_order_status in {"pass", "n/a"}
            and ratio_status == "pass"
            and summary.ranking_contract_status == "pass"
        )
        blocker = evidence.promotion_blocker or "matrix strict gate remains disabled by contract"
        rows.append(
            MatrixBenchmarkAssertion(
                benchmark_id=summary.benchmark_id,
                bench_file=summary.bench_file,
                protein_type=summary.protein_type,
                execution_path=summary.execution_path,
                process_state=summary.process_state,
                target_profile=evidence.target_profile,
                ranking_contract_status=summary.ranking_contract_status,
                coverage=summary.coverage,
                min_coverage=min_coverage,
                top_k=top_k,
                top_k_hits=top_k_hits,
                top_k_status=top_k_status,
                adverse_order_status=adverse_order_status,
                max_ratio=summary.max_ratio,
                ratio_tolerance=ratio_tolerance,
                ratio_status=ratio_status,
                overall_status="pass" if overall_pass else "fail",
                strict_gate_blocked=True,
                blocker=blocker,
            )
        )
    return rows


# ---------------------------------------------------------------------------
# Markdown renderer
# ---------------------------------------------------------------------------

def render_matrix_benchmark_assertions_markdown(
    rows: Iterable[MatrixBenchmarkAssertion],
) -> str:
    assertion_rows = list(rows)
    lines = [
        "# Matrix Benchmark Assertions",
        "",
        "| Benchmark | Protein | Path | Process State | Target Profile | Ranking Contract | Coverage | Min Coverage | Top-k | Top-k Hits | Top-k Status | Adverse Order | Max Ratio | Ratio Tol. | Ratio Status | Overall | Strict Gate Blocked | Blocker |",
        "| --- | --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- | ---: | ---: | --- | --- | --- | --- |",
    ]
    for row in assertion_rows:
        max_ratio = f"{row.max_ratio:.3f}" if row.max_ratio is not None else "n/a"
        lines.append(
            f"| {row.benchmark_id} | {row.protein_type} | {row.execution_path}"
            f" | {row.process_state or 'n/a'} | {row.target_profile}"
            f" | {row.ranking_contract_status} | {row.coverage:.3f}"
            f" | {row.min_coverage:.3f} | {row.top_k} | {row.top_k_hits}"
            f" | {row.top_k_status} | {row.adverse_order_status}"
            f" | {max_ratio} | {row.ratio_tolerance:.3f} | {row.ratio_status}"
            f" | {row.overall_status} | {'yes' if row.strict_gate_blocked else 'no'}"
            f" | {row.blocker} |"
        )
    lines.extend(
        [
            "",
            f"Benchmarks asserted: {len(assertion_rows)}",
            f"Assertion passes: {sum(1 for row in assertion_rows if row.overall_status == 'pass')}",
            "Strict gate remains blocked for all matrix benchmarks by contract until external evidence exists.",
        ]
    )
    return "\n".join(lines) + "\n"
