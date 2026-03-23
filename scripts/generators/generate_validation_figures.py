#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scienceplots

plt.style.use(["science", "no-latex"])

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.generate_benchmark_coverage_gaps import _build_rows as build_coverage_gap_rows
from src.benchmark_validation import build_matrix_promotion_family_status, evaluate_benchmark, summarize_benchmarks


def _status_color(status: str) -> str:
    if status == "pass":
        return "#2a7f62"
    if status == "pass-no-ranking":
        return "#4c956c"
    if status == "partial-pass":
        return "#cc8b3b"
    return "#b03a2e"


def _path_marker(execution_path: str) -> str:
    if execution_path == "matrix_only":
        return "^"
    if execution_path == "matrix_precursor_augmented":
        return "s"
    return "o"


def _path_hatch(execution_path: str) -> str:
    if execution_path == "matrix_only":
        return "///"
    if execution_path == "matrix_precursor_augmented":
        return "xx"
    return ""


def _iter_quantitative_benchmarks(summaries):
    for row in summaries:
        if row.supported and row.matched_compounds > 0:
            yield row


def _build_payload() -> dict[str, object]:
    summaries = summarize_benchmarks()
    readiness = build_matrix_promotion_family_status()
    coverage_rows = build_coverage_gap_rows()
    authoritative_benchmarks = [row for row in summaries if row.execution_path == "free_precursor"]
    quantitative_benchmarks = list(_iter_quantitative_benchmarks(summaries))
    authoritative_points = []
    quantitative_points = []

    for summary in quantitative_benchmarks:
        evaluation = evaluate_benchmark(summary.bench_file)
        for comparison in evaluation.comparisons:
            if comparison.matched_name is None:
                continue
            point = {
                "benchmark_id": summary.benchmark_id,
                "compound": comparison.compound,
                "measured_ppb": comparison.measured_ppb,
                "predicted_ppb": comparison.predicted_ppb,
                "uncertainty_pct": comparison.uncertainty_pct,
                "max_ratio": comparison.ratio,
                "overall_status": summary.overall_status,
                "strict_ready": summary.strict_ready,
                "execution_path": summary.execution_path,
            }
            quantitative_points.append(point)
            if summary.execution_path == "free_precursor":
                authoritative_points.append(point)

    return {
        "benchmark_count": len(summaries),
        "strict_ready_count": sum(1 for row in summaries if row.strict_ready),
        "status_counts": dict(Counter(row.overall_status for row in summaries)),
        "quantitative_benchmarks": [
            {
                "benchmark_id": row.benchmark_id,
                "execution_path": row.execution_path,
                "overall_status": row.overall_status,
                "strict_ready": row.strict_ready,
                "coverage": row.coverage,
                "pearson_r": row.pearson_r,
                "max_ratio": row.max_ratio,
                "mae_ppb": row.mae_ppb,
            }
            for row in quantitative_benchmarks
        ],
        "authoritative_benchmarks": [
            {
                "benchmark_id": row.benchmark_id,
                "overall_status": row.overall_status,
                "strict_ready": row.strict_ready,
                "coverage": row.coverage,
                "pearson_r": row.pearson_r,
                "max_ratio": row.max_ratio,
                "mae_ppb": row.mae_ppb,
            }
            for row in authoritative_benchmarks
        ],
        "authoritative_points": authoritative_points,
        "quantitative_points": quantitative_points,
        "matrix_readiness": [
            {
                "protein_type": row.protein_type,
                "off_flavour_anchor_count": row.off_flavour_anchor_count,
                "meaty_candidate_count": row.meaty_candidate_count,
                "external_meaty_anchor_count": row.external_meaty_anchor_count,
                "candidate_set_ready": row.candidate_set_ready,
                "external_assessment_unlocked": row.external_assessment_unlocked,
                "blocker": row.blocker,
            }
            for row in readiness
        ],
        "coverage_gaps": coverage_rows,
    }


def _render_markdown(payload: dict[str, object]) -> str:
    coverage_rows = payload["coverage_gaps"]
    authoritative_benchmarks = payload["authoritative_benchmarks"]
    authoritative_points = payload["authoritative_points"]
    quantitative_benchmarks = payload["quantitative_benchmarks"]
    quantitative_points = payload["quantitative_points"]
    gap_count = sum(1 for row in coverage_rows if row["status"] == "gap")
    median_ratio = np.median([row["max_ratio"] for row in authoritative_points]) if authoritative_points else float("nan")
    matrix_quantitative = [row for row in quantitative_benchmarks if row["execution_path"] == "matrix_only"]
    augmented_quantitative = [row for row in quantitative_benchmarks if row["execution_path"] == "matrix_precursor_augmented"]

    lines = [
        "# Validation Overview",
        "",
        "This artifact shows the full quantitative benchmark surface while still separating strict-gate free-precursor trust from matrix-only and matrix-augmented evidence.",
        "",
        f"- Benchmarks summarized: {payload['benchmark_count']}",
        f"- Strict-ready benchmarks: {payload['strict_ready_count']}",
        f"- Quantitative benchmarks plotted: {len(quantitative_benchmarks)}",
        f"- Quantitative matched compounds plotted: {len(quantitative_points)}",
        f"- Authoritative free-precursor benchmarks: {len(authoritative_benchmarks)}",
        f"- Quantitative matrix-only benchmarks: {len(matrix_quantitative)}",
        f"- Quantitative matrix-augmented benchmarks: {len(augmented_quantitative)}",
        f"- Authoritative matched compounds: {len(authoritative_points)}",
        f"- Median compound max-ratio in authoritative set: {median_ratio:.3f}",
        f"- Coverage gaps still open: {gap_count}",
        "",
        "How to read the PNG:",
        "",
        "- left: quantitative parity for all numeric benchmarks; circles are free-precursor, triangles are matrix-only, squares are matrix-augmented",
        "- right: per-benchmark quantitative error; the 1.5x line is the strict free-precursor gate and the 2.0x line is the matrix benchmark tolerance",
        "",
        "Matrix readiness, benchmark evidence, and coverage gaps remain in their own dedicated artifacts, but quantitative matrix benchmarks are no longer hidden from the overview.",
    ]
    return "\n".join(lines) + "\n"


def _render_figure(payload: dict[str, object], output_path: Path) -> None:
    quantitative_benchmarks = payload["quantitative_benchmarks"]
    quantitative_points = payload["quantitative_points"]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.8))
    fig.suptitle("Maillard Validation Overview", fontsize=15)

    ax = axes[0]
    measured = np.array([row["measured_ppb"] for row in quantitative_points], dtype=float)
    predicted = np.array([row["predicted_ppb"] for row in quantitative_points], dtype=float)
    if len(measured) == 0:
        raise RuntimeError("No quantitative benchmark points available for validation overview plot.")
    lower = min(np.min(measured[measured > 0.0]), np.min(predicted[predicted > 0.0])) * 0.7
    upper = max(np.max(measured), np.max(predicted)) * 1.4
    reference_line = np.geomspace(lower, upper, 200)
    ax.fill_between(reference_line, reference_line / 1.5, reference_line * 1.5, color="#d9ead3", alpha=0.7, label="within 1.5x")
    ax.plot(reference_line, reference_line, color="#2f4858", linewidth=1.5, label="ideal parity")
    palette = ["#2a7f62", "#4c956c", "#cc8b3b", "#b03a2e", "#2f4858", "#bc4749"]
    benchmark_ids = sorted({row["benchmark_id"] for row in quantitative_points})
    for idx, benchmark_id in enumerate(benchmark_ids):
        subset = [row for row in quantitative_points if row["benchmark_id"] == benchmark_id]
        execution_path = subset[0]["execution_path"]
        ax.errorbar(
            [row["measured_ppb"] for row in subset],
            [row["predicted_ppb"] for row in subset],
            xerr=[row["measured_ppb"] * ((row["uncertainty_pct"] or 0.0) / 100.0) for row in subset],
            fmt=_path_marker(execution_path),
            color=palette[idx % len(palette)],
            ecolor=palette[idx % len(palette)],
            capsize=3,
            markersize=7,
            label=f"{benchmark_id} [{execution_path}]",
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lower, upper)
    ax.set_ylim(lower, upper)
    ax.set_xlabel("Measured concentration (ppb)")
    ax.set_ylabel("Predicted concentration (ppb)")
    ax.set_title("Quantitative Benchmark Parity")
    ax.legend(fontsize=7, loc="lower right")

    ax = axes[1]
    ranked_benchmarks = sorted(
        quantitative_benchmarks,
        key=lambda row: float(row["max_ratio"] or 0.0),
        reverse=True,
    )
    benchmark_labels = [row["benchmark_id"] for row in ranked_benchmarks]
    max_ratios = [row["max_ratio"] if row["max_ratio"] is not None else 0.0 for row in ranked_benchmarks]
    bars = ax.barh(benchmark_labels, max_ratios, color=[_status_color(row["overall_status"]) for row in ranked_benchmarks])
    for bar, row in zip(bars, ranked_benchmarks):
        hatch = _path_hatch(row["execution_path"])
        if hatch:
            bar.set_hatch(hatch)
            bar.set_edgecolor("#2f4858")
            bar.set_linewidth(1.0)
    ax.axvline(1.5, color="#8c1c13", linestyle="--", linewidth=1.5, label="strict threshold")
    ax.axvline(2.0, color="#6c757d", linestyle=":", linewidth=1.5, label="matrix tolerance")
    ax.set_title("Per-Benchmark Quantitative Error")
    ax.set_xlabel("Max measured/predicted ratio")
    ax.legend(loc="lower right")

    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    payload = _build_payload()
    markdown = _render_markdown(payload)

    png_path = output_dir / "validation_overview.png"
    md_path = output_dir / "validation_overview.md"
    json_path = output_dir / "validation_overview.json"

    _render_figure(payload, png_path)
    md_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(markdown)
    print(f"Wrote {png_path}")
    print(f"Wrote {md_path}")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
