#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.generate_benchmark_coverage_gaps import _build_rows as build_coverage_gap_rows
from src.benchmark_validation import build_matrix_promotion_family_status, summarize_benchmarks


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

    return {
        "benchmark_count": len(summaries),
        "strict_ready_count": sum(1 for row in summaries if row.strict_ready),
        "status_counts": dict(Counter(row.overall_status for row in summaries)),
        "benchmarks": [
            {
                "benchmark_id": row.benchmark_id,
                "execution_path": row.execution_path,
                "status": row.overall_status,
                "strict_ready": row.strict_ready,
                "coverage": row.coverage,
                "pearson_r": row.pearson_r,
                "max_ratio": row.max_ratio,
                "mae_ppb": row.mae_ppb,
            }
            for row in summaries
        ],
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
    status_counts = payload["status_counts"]
    readiness = payload["matrix_readiness"]
    coverage_rows = payload["coverage_gaps"]
    gap_count = sum(1 for row in coverage_rows if row["status"] == "gap")
    unlocked = sum(1 for row in readiness if row["external_assessment_unlocked"])
    lines = [
        "# Validation Overview",
        "",
        "This artifact summarizes current reliability and current limits in a single place.",
        "",
        f"- Benchmarks summarized: {payload['benchmark_count']}",
        f"- Strict-ready benchmarks: {payload['strict_ready_count']}",
        f"- Status distribution: {status_counts}",
        f"- Matrix families externally unlocked: {unlocked}",
        f"- Coverage gaps still open: {gap_count}",
        "",
        "The accompanying PNG combines four views:",
        "",
        "- benchmark max-ratio behavior versus the strict free-precursor threshold",
        "- benchmark status distribution",
        "- matrix readiness by protein family",
        "- open coverage gaps that still limit scientific scope",
        "",
        "Interpretation:",
        "",
        "- free-precursor systems form the current high-trust envelope",
        "- matrix candidate sets exist for pea and soy",
        "- external meaty-positive matrix evidence is still the main blocker",
    ]
    return "\n".join(lines) + "\n"


def _render_figure(payload: dict[str, object], output_path: Path) -> None:
    benchmarks = payload["benchmarks"]
    readiness = payload["matrix_readiness"]
    coverage_rows = payload["coverage_gaps"]

    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle("Maillard Validation Overview", fontsize=16)

    ax = axes[0, 0]
    bench_labels = [row["benchmark_id"] for row in benchmarks]
    max_ratios = [row["max_ratio"] if row["max_ratio"] is not None else 0.0 for row in benchmarks]
    colors = [_status_color(row["status"]) for row in benchmarks]
    ax.barh(bench_labels, max_ratios, color=colors)
    ax.axvline(1.5, color="#8c1c13", linestyle="--", linewidth=1.5, label="strict free-AA threshold")
    ax.set_title("Benchmark Max Ratio")
    ax.set_xlabel("Measured / predicted max ratio")
    ax.legend(loc="lower right")

    ax = axes[0, 1]
    status_counts = Counter(row["status"] for row in benchmarks)
    labels = list(status_counts.keys())
    values = [status_counts[label] for label in labels]
    ax.bar(labels, values, color=[_status_color(label) for label in labels])
    ax.set_title("Benchmark Status Distribution")
    ax.set_ylabel("Benchmarks")
    ax.tick_params(axis="x", rotation=20)

    ax = axes[1, 0]
    proteins = [row["protein_type"] for row in readiness]
    off_counts = [row["off_flavour_anchor_count"] for row in readiness]
    meaty_counts = [row["meaty_candidate_count"] for row in readiness]
    external_counts = [row["external_meaty_anchor_count"] for row in readiness]
    xpos = range(len(proteins))
    ax.bar([x - 0.25 for x in xpos], off_counts, width=0.25, label="off-flavour anchors", color="#487a6a")
    ax.bar(xpos, meaty_counts, width=0.25, label="meaty candidates", color="#c98f3d")
    ax.bar([x + 0.25 for x in xpos], external_counts, width=0.25, label="external meaty anchors", color="#8c1c13")
    ax.set_xticks(list(xpos), proteins)
    ax.set_title("Matrix Readiness By Protein Family")
    ax.set_ylabel("Benchmarks")
    ax.legend()

    ax = axes[1, 1]
    gap_rows = [row for row in coverage_rows if row["status"] == "gap"]
    gap_labels = [f"{row['dimension']}:{row['category']}" for row in gap_rows]
    gap_values = [row["benchmark_count"] for row in gap_rows]
    ax.barh(gap_labels, gap_values, color="#b03a2e")
    ax.set_title("Open Coverage Gaps")
    ax.set_xlabel("Benchmarks present")

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
