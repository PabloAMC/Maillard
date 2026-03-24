#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scienceplots  # noqa: F401

matplotlib.use("Agg")
plt.style.use(["science"])  # full LaTeX rendering

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.generators.generate_benchmark_coverage_gaps import _build_rows as build_coverage_gap_rows
from src.benchmark_validation import build_matrix_promotion_family_status, evaluate_benchmark, summarize_benchmarks

# Human-readable benchmark labels (LaTeX-safe)
_BENCH_LABELS: dict[str, str] = {
    "acrylamide_asparagine_glucose_Parker2012":                  r"Asn + glucose acrylamide (Parker 2012)",
    "cys_glucose_150C_Farmer1999":                               r"Cys + glucose, $150\,^\circ$C (Farmer 1999)",
    "cys_ribose_140C_Hofmann1998":                               r"Cys + ribose, $140\,^\circ$C (Hofmann 1998)",
    "cys_ribose_150C_Mottram1994":                               r"Cys + ribose, $150\,^\circ$C (Mottram 1994)",
    "pea_isolate_40C_PratapSingh2021":                           r"Pea isolate, $40\,^\circ$C (Pratap Singh 2021)",
    "pea_isolate_ribose_cysteine_100C_45min_Internal2026":       r"Pea iso.\ + Rib/Cys, $100\,^\circ$C (internal)",
    "pea_isolate_uht_140C_Trikusuma2019":                        r"Pea isolate UHT, $140\,^\circ$C (Trikusuma 2019)",
    "soy_isolate_40C_PratapSingh2021":                           r"Soy isolate, $40\,^\circ$C (Pratap Singh 2021)",
    "soy_isolate_ribose_cysteine_100C_45min_Internal2026":       r"Soy iso.\ + Rib/Cys, $100\,^\circ$C (internal)",
}

def _bench_label(benchmark_id: str) -> str:
    return _BENCH_LABELS.get(benchmark_id, benchmark_id.replace("_", r"\_"))


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

    fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))
    fig.suptitle(r"\textbf{Maillard model --- benchmark-level quantitative accuracy}", fontsize=12, y=0.99)
    fig.subplots_adjust(left=0.06, right=0.97, top=0.93, bottom=0.10, wspace=0.30)

    ax = axes[0]
    measured = np.array([row["measured_ppb"] for row in quantitative_points], dtype=float)
    predicted = np.array([row["predicted_ppb"] for row in quantitative_points], dtype=float)
    if len(measured) == 0:
        raise RuntimeError("No quantitative benchmark points available for validation overview plot.")
    lower = min(np.min(measured[measured > 0.0]), np.min(predicted[predicted > 0.0])) * 0.5
    upper = max(np.max(measured), np.max(predicted)) * 2.0
    reference_line = np.geomspace(lower, upper, 200)
    ax.fill_between(reference_line, reference_line / 2.0, reference_line * 2.0, color="#FFF3CD", alpha=0.55, label=r"Within $2\times$")
    ax.fill_between(reference_line, reference_line / 1.5, reference_line * 1.5, color="#D9EAD3", alpha=0.70, label=r"Within $1.5\times$")
    ax.plot(reference_line, reference_line, color="#2f4858", linewidth=1.5, label=r"Ideal parity ($y = x$)")
    palette = ["#0077BB", "#EE7733", "#009988", "#CC3311", "#AA4499", "#DDCC77", "#882255", "#44AA99", "#332288"]
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
            capsize=2.5,
            markersize=6.5,
            linewidth=0.8,
            label=_bench_label(benchmark_id),
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lower, upper)
    ax.set_ylim(lower, upper)
    ax.set_xlabel(r"Measured concentration (ppb)")
    ax.set_ylabel(r"Predicted concentration (ppb)")
    ax.set_title(r"Quantitative parity --- all benchmark compounds")
    # Place legend inside the left panel without overlapping the parity markers.
    ax.legend(fontsize=6.0, loc="upper left", bbox_to_anchor=(0.02, 0.98),
              bbox_transform=ax.transAxes, framealpha=0.85, edgecolor="0.75")

    ax = axes[1]
    ranked_benchmarks = sorted(
        quantitative_benchmarks,
        key=lambda row: float(row["max_ratio"] or 0.0),
        reverse=True,
    )
    benchmark_labels = [_bench_label(row["benchmark_id"]) for row in ranked_benchmarks]
    max_ratios = [row["max_ratio"] if row["max_ratio"] is not None else 0.0 for row in ranked_benchmarks]
    bars = ax.barh(
        benchmark_labels, max_ratios,
        color=[_status_color(row["overall_status"]) for row in ranked_benchmarks],
        edgecolor="white", linewidth=0.5,
    )
    for bar, row in zip(bars, ranked_benchmarks):
        hatch = _path_hatch(row["execution_path"])
        if hatch:
            bar.set_hatch(hatch)
            bar.set_edgecolor("#2f4858")
            bar.set_linewidth(0.8)
        ax.text(
            (row["max_ratio"] or 0.0) + 0.04,
            bar.get_y() + bar.get_height() / 2.0,
            f"{row['max_ratio']:.2f}$\\times$" if row["max_ratio"] else "",
            va="center", fontsize=6.5,
        )
    ax.axvline(1.5, color="#8c1c13", linestyle="--", linewidth=1.4, label=r"Strict gate ($1.5\times$)")
    ax.axvline(2.0, color="#555555", linestyle=":",  linewidth=1.4, label=r"Matrix tolerance ($2\times$)")
    ax.set_title(r"Worst-case ratio per benchmark")
    ax.set_xlabel(r"Max ratio (predicted / measured)")
    # Move the right-panel legend inside the axes in the top-right corner.
    ax.legend(loc="upper right", fontsize=6.5, bbox_to_anchor=(0.98, 0.98),
              bbox_transform=ax.transAxes, framealpha=0.85, edgecolor="0.75")

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
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
