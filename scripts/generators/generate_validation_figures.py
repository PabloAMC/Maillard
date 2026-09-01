#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import shutil
import sys
from collections import Counter
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths
from src.benchmark_coverage_gaps import build_benchmark_coverage_gap_rows as build_coverage_gap_rows
from src.benchmark_labels import benchmark_label
from src.benchmark_validation import build_matrix_promotion_family_status, evaluate_benchmark, summarize_benchmarks
from src.family_validation_overview import build_family_validation_overview_artifact
from src.plot_style import configure_science_plot_style

_FAMILY_LABELS: dict[str, str] = {
    "amino_acid_sugar_core": r"F01 Amino acid + sugar",
    "lipid_oxidation_and_carbonylic_crosstalk": r"F02 Lipid oxidation",
    "thiamine_fragmentation_support": r"F03 Thiamine",
    "nucleotide_and_ribose_support": r"F04 Nucleotide support",
    "glutathione_and_peptide_support": r"F05 Peptide support",
    "alternative_protein_matrix_scope": r"F06 Matrix scope",
    "carbonyl_donor_hierarchy": r"F07 Donor hierarchy",
    "off_note_and_maillard_suppression": r"F08 Off-notes",
    "carbohydrate_pyrolysis_and_caramelization": r"F09 Caramelisation",
    "fermentation_pretreatment": r"F10 Pretreatment",
    "lipid_maillard_crosstalk": r"F11 Lipid crosstalk",
    "protein_damage_markers": r"F12 Damage markers",
}

_FAMILY_COLORS: dict[str, str] = {
    "amino_acid_sugar_core": "#0077BB",
    "lipid_oxidation_and_carbonylic_crosstalk": "#EE7733",
    "thiamine_fragmentation_support": "#228833",
    "nucleotide_and_ribose_support": "#AA4499",
    "glutathione_and_peptide_support": "#66A61E",
    "alternative_protein_matrix_scope": "#8C564B",
    "carbonyl_donor_hierarchy": "#DDCC77",
    "off_note_and_maillard_suppression": "#CC3311",
    "carbohydrate_pyrolysis_and_caramelization": "#009988",
    "fermentation_pretreatment": "#CC79A7",
    "lipid_maillard_crosstalk": "#4477AA",
    "protein_damage_markers": "#882255",
}
_FAMILY_COLOR_NONE = "#BBBBBB"

def _bench_label(benchmark_id: str) -> str:
    return benchmark_label(benchmark_id, style="latex")


def _compact_point_summary(point: dict[str, object] | None) -> dict[str, object] | None:
    if not point:
        return None
    return {
        "benchmark_id": point["benchmark_id"],
        "benchmark_label": _bench_label(str(point["benchmark_id"])),
        "compound": point["compound"],
        "measured_ppb": point["measured_ppb"],
        "predicted_ppb": point["predicted_ppb"],
        "max_ratio": point["max_ratio"],
        "reference_signal_origin": point["reference_signal_origin"],
    }


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
        if not (row.supported and row.matched_compounds > 0):
            continue
        try:
            benchmark_payload = json.loads(Path(row.bench_file).read_text(encoding="utf-8"))
        except (FileNotFoundError, json.JSONDecodeError, OSError, TypeError):
            benchmark_payload = {}
        tier = str(benchmark_payload.get("metadata", {}).get("tier", "")).strip().upper()
        if tier == "SECONDARY":
            continue
        yield row


def _build_payload() -> dict[str, object]:
    summaries = summarize_benchmarks()
    readiness = build_matrix_promotion_family_status()
    coverage_rows = build_coverage_gap_rows()
    family_overview = build_family_validation_overview_artifact()
    authoritative_benchmarks = [row for row in summaries if row.execution_path == "free_precursor"]
    quantitative_benchmarks = list(_iter_quantitative_benchmarks(summaries))
    experimental_quantitative_benchmarks = [row for row in quantitative_benchmarks if row.reference_signal_origin == "measured_volatiles"]
    reference_quantitative_benchmarks = [row for row in quantitative_benchmarks if row.reference_signal_origin != "measured_volatiles"]
    inside_1_5x_benchmark_count = sum(1 for row in quantitative_benchmarks if float(row.max_ratio or 0.0) <= 1.5)
    inside_2x_benchmark_count = sum(1 for row in quantitative_benchmarks if float(row.max_ratio or 0.0) <= 2.0)
    experimental_inside_1_5x_benchmark_count = sum(1 for row in experimental_quantitative_benchmarks if float(row.max_ratio or 0.0) <= 1.5)
    experimental_inside_2x_benchmark_count = sum(1 for row in experimental_quantitative_benchmarks if float(row.max_ratio or 0.0) <= 2.0)
    authoritative_points = []
    quantitative_points = []

    for summary in quantitative_benchmarks:
        evaluation = evaluate_benchmark(summary.bench_file)
        for comparison in evaluation.comparisons:
            if comparison.matched_name is None:
                continue
            point = {
                "benchmark_id": summary.benchmark_id,
                "benchmark_label": _bench_label(summary.benchmark_id),
                "compound": comparison.compound,
                "measured_ppb": comparison.measured_ppb,
                "predicted_ppb": comparison.predicted_ppb,
                "uncertainty_pct": comparison.uncertainty_pct,
                "max_ratio": comparison.ratio,
                "overall_status": summary.overall_status,
                "strict_ready": summary.strict_ready,
                "execution_path": summary.execution_path,
                "reference_signal_origin": summary.reference_signal_origin,
            }
            quantitative_points.append(point)
            if summary.execution_path == "free_precursor":
                authoritative_points.append(point)

    points_by_benchmark = {}
    for point in quantitative_points:
        points_by_benchmark.setdefault(str(point["benchmark_id"]), []).append(point)

    def _select_representative_point(summary_row, points):
        if summary_row is None:
            return None
        benchmark_points = points_by_benchmark.get(str(summary_row.benchmark_id), [])
        return max(benchmark_points, key=lambda row: float(row["max_ratio"]), default=None)

    # CONSISTENCY FIX (2026-08-27, audit remediation Part 2a): `worst_quantitative_ratio`
    # is a max over ALL quantitative benchmarks, but `worst_quantitative_point` used to be
    # selected from `reference_quantitative_benchmarks` (the reference_volatiles subset)
    # only. The two happened to coincide until the Bolton1994 benchmark entered the panel,
    # at which point the artifact would report a "worst quantitative benchmark ratio" of
    # ~163x next to a "worst quantitative point" of Cerny at ~1.6x. The pair is now drawn
    # from the same population (all quantitative benchmarks); the reference-only view it
    # used to provide is preserved below under an explicitly-named key.
    worst_quantitative_benchmark = max(quantitative_benchmarks, key=lambda row: float(row.max_ratio or 0.0), default=None)
    worst_quantitative_point = _select_representative_point(worst_quantitative_benchmark, quantitative_points)
    reference_quantitative_points = [
        row for row in quantitative_points if row["reference_signal_origin"] != "measured_volatiles"
    ]
    reference_worst_benchmark = max(
        reference_quantitative_benchmarks,
        key=lambda row: float(row.max_ratio or 0.0),
        default=None,
    )
    reference_worst_quantitative_point = _select_representative_point(
        reference_worst_benchmark, reference_quantitative_points
    )
    experimental_quantitative_points = [
        row for row in quantitative_points if row["reference_signal_origin"] == "measured_volatiles"
    ]
    experimental_worst_benchmark = max(
        experimental_quantitative_benchmarks,
        key=lambda row: float(row.max_ratio or 0.0),
        default=None,
    )
    experimental_worst_quantitative_point = _select_representative_point(experimental_worst_benchmark, experimental_quantitative_points)

    return {
        "benchmark_count": len(summaries),
        "strict_ready_count": sum(1 for row in summaries if row.strict_ready),
        "status_counts": dict(Counter(row.overall_status for row in summaries)),
        "inside_1_5x_benchmark_count": inside_1_5x_benchmark_count,
        "inside_2x_benchmark_count": inside_2x_benchmark_count,
        "outside_1_5x_benchmark_count": max(0, len(quantitative_benchmarks) - inside_1_5x_benchmark_count),
        "outside_2x_benchmark_count": max(0, len(quantitative_benchmarks) - inside_2x_benchmark_count),
        "worst_quantitative_ratio": max((float(row.max_ratio or 0.0) for row in quantitative_benchmarks), default=0.0),
        "experimental_quantitative_benchmark_count": len(experimental_quantitative_benchmarks),
        "reference_quantitative_benchmark_count": len(reference_quantitative_benchmarks),
        "experimental_inside_1_5x_benchmark_count": experimental_inside_1_5x_benchmark_count,
        "experimental_outside_1_5x_benchmark_count": max(0, len(experimental_quantitative_benchmarks) - experimental_inside_1_5x_benchmark_count),
        "experimental_inside_2x_benchmark_count": experimental_inside_2x_benchmark_count,
        "experimental_outside_2x_benchmark_count": max(0, len(experimental_quantitative_benchmarks) - experimental_inside_2x_benchmark_count),
        "experimental_worst_quantitative_ratio": max((float(row.max_ratio or 0.0) for row in experimental_quantitative_benchmarks), default=0.0),
        "worst_quantitative_point": _compact_point_summary(worst_quantitative_point),
        "worst_quantitative_population": "all_quantitative_benchmarks",
        "reference_worst_quantitative_ratio": max((float(row.max_ratio or 0.0) for row in reference_quantitative_benchmarks), default=0.0),
        "reference_worst_quantitative_point": _compact_point_summary(reference_worst_quantitative_point),
        "experimental_worst_quantitative_point": _compact_point_summary(experimental_worst_quantitative_point),
        "quantitative_benchmarks": [
            {
                "benchmark_id": row.benchmark_id,
                "execution_path": row.execution_path,
                "overall_status": row.overall_status,
                "strict_ready": row.strict_ready,
                "reference_signal_origin": row.reference_signal_origin,
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
        "family_overview": family_overview,
        "integrated_family_count": int(family_overview.get("summary", {}).get("integrated_family_count", 0)),
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
    worst_point = payload.get("worst_quantitative_point")
    reference_worst_point = payload.get("reference_worst_quantitative_point")
    experimental_worst_point = payload.get("experimental_worst_quantitative_point")

    lines = [
        '<!-- Auto-regenerated by ./scripts/docker_maillard.sh run "python scripts/generators/generate_validation_figures.py". Manual edits will be overwritten. -->',
        "",
        "# Validation Overview",
        "",
        "This artifact shows the full quantitative benchmark surface while still separating strict-gate free-precursor trust from matrix-only and matrix-augmented evidence.",
        "",
        f"- Benchmarks summarized: {payload['benchmark_count']}",
        f"- Strict-ready benchmarks: {payload['strict_ready_count']}",
        f"- Quantitative benchmarks plotted: {len(quantitative_benchmarks)}",
        f"- Experimental quantitative benchmarks: {payload['experimental_quantitative_benchmark_count']}",
        f"- Reference-only quantitative anchors: {payload['reference_quantitative_benchmark_count']}",
        f"- Experimental benchmarks inside 1.5x: {payload['experimental_inside_1_5x_benchmark_count']}",
        f"- Experimental benchmarks outside 1.5x: {payload['experimental_outside_1_5x_benchmark_count']}",
        f"- Experimental benchmarks outside 2x: {payload['experimental_outside_2x_benchmark_count']}",
        f"- Worst experimental benchmark ratio: {float(payload['experimental_worst_quantitative_ratio']):.3f}x",
        f"- Worst experimental point: {experimental_worst_point['benchmark_label']} / {experimental_worst_point['compound']} ({float(experimental_worst_point['max_ratio']):.3f}x)" if experimental_worst_point else "- Worst experimental point: n/a",
        f"- Quantitative benchmarks inside 1.5x: {payload['inside_1_5x_benchmark_count']}",
        f"- Quantitative benchmarks outside 1.5x: {payload['outside_1_5x_benchmark_count']}",
        f"- Quantitative benchmarks outside 2x: {payload['outside_2x_benchmark_count']}",
        f"- Worst quantitative benchmark ratio (all quantitative benchmarks): {float(payload['worst_quantitative_ratio']):.3f}x",
        f"- Worst quantitative point (all quantitative benchmarks): {worst_point['benchmark_label']} / {worst_point['compound']} ({float(worst_point['max_ratio']):.3f}x; {worst_point['reference_signal_origin']})" if worst_point else "- Worst quantitative point (all quantitative benchmarks): n/a",
        f"- Worst reference-only benchmark ratio: {float(payload['reference_worst_quantitative_ratio']):.3f}x",
        f"- Worst reference-only point: {reference_worst_point['benchmark_label']} / {reference_worst_point['compound']} ({float(reference_worst_point['max_ratio']):.3f}x; {reference_worst_point['reference_signal_origin']})" if reference_worst_point else "- Worst reference-only point: n/a",
        f"- Quantitative matched compounds plotted: {len(quantitative_points)}",
        f"- Authoritative free-precursor benchmarks: {len(authoritative_benchmarks)}",
        f"- Quantitative matrix-only benchmarks: {len(matrix_quantitative)}",
        f"- Quantitative matrix-augmented benchmarks: {len(augmented_quantitative)}",
        f"- Integrated runtime families tracked in the overview artifact: {int(payload.get('integrated_family_count', 0))}",
        f"- Authoritative matched compounds: {len(authoritative_points)}",
        f"- Median compound max-ratio in authoritative set: {median_ratio:.3f}",
        f"- Coverage gaps still open: {gap_count}",
        "",
        "How to read the PNG:",
        "",
        "- single panel: quantitative parity for all numeric benchmarks",
        "- colour: benchmark / literature system, using formatted study references rather than raw benchmark ids",
        "- marker shape: circles are free-precursor, triangles are matrix-only, squares are matrix-augmented",
        "- filled markers denote wet-lab measured comparators; hollow markers denote reference-only anchors",
        "- green/yellow bands denote the 1.5x and 2x tolerance envelopes around ideal parity",
        "- families without direct matched numeric compounds are tracked in the family overview artifact, not forced into this benchmark-level scatter",
        "",
        "Matrix readiness, benchmark evidence, family integration tiers, and per-benchmark worst-case error remain in their own dedicated artifacts, while this overview stays single-panel for README-safe embedding.",
    ]
    return "\n".join(lines) + "\n"


def _render_figure(payload: dict[str, object], output_path: Path) -> None:
    configure_science_plot_style()
    quantitative_benchmarks = payload["quantitative_benchmarks"]
    quantitative_points = payload["quantitative_points"]

    fig, ax = plt.subplots(figsize=(10.5, 6.5))

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
        reference_only = subset[0]["reference_signal_origin"] != "measured_volatiles"
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
            markerfacecolor="white" if reference_only else palette[idx % len(palette)],
            markeredgecolor=palette[idx % len(palette)],
            markeredgewidth=1.0,
            alpha=0.95 if reference_only else 0.90,
            label=_bench_label(benchmark_id),
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lower, upper)
    ax.set_ylim(lower, upper)
    ax.set_xlabel(r"Measured concentration (ppb)")
    ax.set_ylabel(r"Predicted concentration (ppb)")
    ax.set_title(r"Predicted vs. measured concentration --- all benchmark compounds")

    execution_handles = [
        plt.Line2D([0], [0], marker="o", color="0.35", linestyle="None", markersize=6, label=r"Free precursor"),
        plt.Line2D([0], [0], marker="^", color="0.35", linestyle="None", markersize=6, label=r"Matrix-only"),
        plt.Line2D([0], [0], marker="s", color="0.35", linestyle="None", markersize=6, label=r"Matrix + precursor"),
        plt.Line2D([0], [0], marker="o", color="0.35", markerfacecolor="white", markeredgecolor="0.35", linestyle="None", markersize=6, label=r"Reference-only anchor"),
    ]
    band_handles = [
        mpatches.Patch(facecolor="#D9EAD3", alpha=0.70, label=r"Within $1.5\times$"),
        mpatches.Patch(facecolor="#FFF3CD", alpha=0.55, label=r"Within $2\times$"),
        plt.Line2D([0], [0], color="#2f4858", linewidth=1.5, label=r"Ideal parity ($y=x$)"),
    ]
    benchmark_handles = [
        plt.Line2D(
            [0],
            [0],
            marker=_path_marker(next(row["execution_path"] for row in quantitative_points if row["benchmark_id"] == benchmark_id)),
            color=palette[idx % len(palette)],
            markerfacecolor=(
                "white"
                if next(row["reference_signal_origin"] for row in quantitative_points if row["benchmark_id"] == benchmark_id) != "measured_volatiles"
                else palette[idx % len(palette)]
            ),
            markeredgecolor=palette[idx % len(palette)],
            markeredgewidth=1.0,
            linestyle="None",
            markersize=6,
            label=_bench_label(benchmark_id),
        )
        for idx, benchmark_id in enumerate(benchmark_ids)
    ]
    helper_legend = ax.legend(
        handles=band_handles + execution_handles,
        fontsize=6.0,
        loc="upper left",
        bbox_to_anchor=(0.02, 0.98),
        bbox_transform=ax.transAxes,
        framealpha=0.90,
        edgecolor="0.75",
        title=r"\textit{Guide}",
        title_fontsize=6.0,
    )
    ax.add_artist(helper_legend)
    ax.legend(
        handles=benchmark_handles,
        fontsize=5.7,
        loc="center left",
        bbox_to_anchor=(1.02, 0.50),
        framealpha=0.90,
        edgecolor="0.75",
        title=r"\textit{Benchmarks / references}",
        title_fontsize=6.0,
        labelspacing=0.18,
        handletextpad=0.30,
        borderpad=0.25,
        ncol=1,
    )

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default=data_paths.rel(data_paths.VALIDATION_DIR))
    parser.add_argument("--docs-asset-dir", default=data_paths.rel(data_paths.DOCS_ASSETS_DIR))
    parser.add_argument(
        "--skip-figures",
        action="store_true",
        help=(
            "write validation_overview.md/.json only. Added 2026-08-27 (Wave H) so the "
            "artifact can be regenerated on a machine without dvipng; the PNG stays stale "
            "and the staleness is stated rather than hidden."
        ),
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    docs_asset_dir = ROOT / args.docs_asset_dir
    docs_asset_dir.mkdir(parents=True, exist_ok=True)

    payload = _build_payload()
    markdown = _render_markdown(payload)

    png_path = output_dir / "validation_overview.png"
    md_path = output_dir / "validation_overview.md"
    json_path = output_dir / "validation_overview.json"

    md_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    if args.skip_figures:
        print(markdown)
        print(f"Wrote {md_path}")
        print(f"Wrote {json_path}")
        print("Figure skipped (--skip-figures).")
        return 0

    _render_figure(payload, png_path)
    shutil.copyfile(png_path, docs_asset_dir / "validation_overview.png")

    print(markdown)
    print(f"Wrote {png_path}")
    print(f"Copied validation_overview.png to {docs_asset_dir / 'validation_overview.png'}")
    print(f"Wrote {md_path}")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
