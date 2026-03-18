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


def _build_payload() -> dict[str, object]:
    summaries = summarize_benchmarks()
    readiness = build_matrix_promotion_family_status()
    coverage_rows = build_coverage_gap_rows()
    authoritative_benchmarks = [row for row in summaries if row.execution_path == "free_precursor"]
    authoritative_points = []

    for summary in authoritative_benchmarks:
        evaluation = evaluate_benchmark(summary.bench_file)
        for comparison in evaluation.comparisons:
            if comparison.matched_name is None:
                continue
            authoritative_points.append(
                {
                    "benchmark_id": summary.benchmark_id,
                    "compound": comparison.compound,
                    "measured_ppb": comparison.measured_ppb,
                    "predicted_ppb": comparison.predicted_ppb,
                    "uncertainty_pct": comparison.uncertainty_pct,
                    "max_ratio": comparison.ratio,
                    "overall_status": summary.overall_status,
                    "strict_ready": summary.strict_ready,
                }
            )

    return {
        "benchmark_count": len(summaries),
        "strict_ready_count": sum(1 for row in summaries if row.strict_ready),
        "status_counts": dict(Counter(row.overall_status for row in summaries)),
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
    gap_count = sum(1 for row in coverage_rows if row["status"] == "gap")
    median_ratio = np.median([row["max_ratio"] for row in authoritative_points]) if authoritative_points else float("nan")

    lines = [
        "# Validation Overview",
        "",
        "This artifact now focuses only on the two first-pass trust panels: quantitative parity and per-benchmark error.",
        "",
        f"- Benchmarks summarized: {payload['benchmark_count']}",
        f"- Strict-ready benchmarks: {payload['strict_ready_count']}",
        f"- Authoritative free-precursor benchmarks: {len(authoritative_benchmarks)}",
        f"- Authoritative matched compounds: {len(authoritative_points)}",
        f"- Median compound max-ratio in authoritative set: {median_ratio:.3f}",
        f"- Coverage gaps still open: {gap_count}",
        "",
        "How to read the PNG:",
        "",
        "- left: proof surface for the validated free-precursor envelope",
        "- right: per-benchmark quantitative error against the 1.5x contract",
        "",
        "Matrix readiness, benchmark evidence, and coverage gaps remain in their own dedicated artifacts.",
    ]
    return "\n".join(lines) + "\n"


def _render_figure(payload: dict[str, object], output_path: Path) -> None:
    authoritative_benchmarks = payload["authoritative_benchmarks"]
    authoritative_points = payload["authoritative_points"]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.8))
    fig.suptitle("Maillard Validation Overview", fontsize=15)

    ax = axes[0]
    measured = np.array([row["measured_ppb"] for row in authoritative_points], dtype=float)
    predicted = np.array([row["predicted_ppb"] for row in authoritative_points], dtype=float)
    uncertainties = np.array([
        row["uncertainty_pct"] if row["uncertainty_pct"] is not None else 0.0
        for row in authoritative_points
    ], dtype=float)
    measured_err = measured * uncertainties / 100.0
    lower = min(np.min(measured[measured > 0.0]), np.min(predicted[predicted > 0.0])) * 0.7
    upper = max(np.max(measured), np.max(predicted)) * 1.4
    reference_line = np.geomspace(lower, upper, 200)
    ax.fill_between(reference_line, reference_line / 1.5, reference_line * 1.5, color="#d9ead3", alpha=0.7, label="within 1.5x")
    ax.plot(reference_line, reference_line, color="#2f4858", linewidth=1.5, label="ideal parity")
    palette = ["#2a7f62", "#4c956c", "#cc8b3b", "#b03a2e"]
    benchmark_ids = sorted({row["benchmark_id"] for row in authoritative_points})
    for idx, benchmark_id in enumerate(benchmark_ids):
        subset = [row for row in authoritative_points if row["benchmark_id"] == benchmark_id]
        ax.errorbar(
            [row["measured_ppb"] for row in subset],
            [row["predicted_ppb"] for row in subset],
            xerr=[row["measured_ppb"] * ((row["uncertainty_pct"] or 0.0) / 100.0) for row in subset],
            fmt="o",
            color=palette[idx % len(palette)],
            ecolor=palette[idx % len(palette)],
            capsize=3,
            label=benchmark_id,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lower, upper)
    ax.set_ylim(lower, upper)
    ax.set_xlabel("Measured concentration (ppb)")
    ax.set_ylabel("Predicted concentration (ppb)")
    ax.set_title("Authoritative Free-Precursor Parity")
    ax.legend(fontsize=8, loc="lower right")

    ax = axes[1]
    benchmark_labels = [row["benchmark_id"] for row in authoritative_benchmarks]
    max_ratios = [row["max_ratio"] if row["max_ratio"] is not None else 0.0 for row in authoritative_benchmarks]
    ax.barh(benchmark_labels, max_ratios, color=[_status_color(row["overall_status"]) for row in authoritative_benchmarks])
    ax.axvline(1.5, color="#8c1c13", linestyle="--", linewidth=1.5, label="strict threshold")
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
    ax.set_xlabel("Max measured/predicted ratio")
    ax.legend(loc="lower right")

    ax = axes[1, 0]
    proteins = [row["protein_type"] for row in readiness]
    off_counts = [row["off_flavour_anchor_count"] for row in readiness]
    meaty_counts = [row["meaty_candidate_count"] for row in readiness]
    external_counts = [row["external_meaty_anchor_count"] for row in readiness]
    xpos = np.arange(len(proteins))
    ax.bar(xpos - 0.25, off_counts, width=0.25, label="off-flavour anchors", color="#487a6a")
    ax.bar(xpos, meaty_counts, width=0.25, label="meaty candidates", color="#c98f3d")
    ax.bar(xpos + 0.25, external_counts, width=0.25, label="external meaty anchors", color="#8c1c13")
    ax.set_xticks(xpos)
    ax.set_xticklabels(proteins)
    ax.set_ylabel("Benchmarks")
    ax.set_title("Matrix Readiness")
    ax.legend()

    ax = axes[1, 1]
    ax.axis("off")
    gap_rows = [row for row in coverage_rows if row["status"] == "gap"]
    notes = [
        "Open limits:",
        "",
    ]
    notes.extend(f"- {row['dimension']}: {row['category']}" for row in gap_rows)
    ax.text(0.0, 1.0, "\n".join(notes), va="top", ha="left", fontsize=10)
    ax.set_title("What Still Blocks Broader Trust")

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
