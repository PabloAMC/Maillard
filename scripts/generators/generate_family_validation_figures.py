#!/usr/bin/env python3
"""Generates three standalone validation figures for the family overview.

Figures produced
----------------
family_parity.png
    Concentration parity scatter — all matched compounds, coloured by chemistry
    family, marker shape indicating execution lane.

family_benchmark_accuracy.png
    Per-benchmark worst-case ratio (predicted / measured) as a horizontal bar
    chart with human-readable study labels.

family_coverage.png
    Coverage census across all tracked chemistry families showing how many
    quantitative compound points each family currently has.

Each figure is saved independently at (8 x 6) inches so it can be placed in
documentation or publications without cropping.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import summarize_benchmarks  # noqa: E402
from src.benchmark_labels import benchmark_label  # noqa: E402
from src.plot_style import configure_science_plot_style  # noqa: E402
from src.family_validation_overview import (  # noqa: E402
    build_family_validation_overview_artifact,
    render_family_validation_overview_markdown,
)

configure_science_plot_style()

# ---------------------------------------------------------------------------
# Human-readable label tables (LaTeX-safe: & → \&, % → \%)
# ---------------------------------------------------------------------------

# Selected short labels for tracked SLR families used in the coverage figure
_FAMILY_LABELS: dict[str, str] = {
    "amino_acid_sugar_core":                     r"F01 \ Amino acid + sugar Maillard core",
    "lipid_oxidation_and_carbonylic_crosstalk":  r"F02 \ Lipid oxidation \& carbonylic crosstalk",
    "thiamine_fragmentation_support":            r"F03 \ Thiamine degradation",
    "nucleotide_and_ribose_support":             r"F04 \ Nucleotide \& ribose support",
    "glutathione_and_peptide_support":           r"F05 \ Glutathione \& peptide support",
    "alternative_protein_matrix_scope":          r"F06 \ Alternative protein matrix scope",
    "carbonyl_donor_hierarchy":                  r"F07 \ Carbonyl donor hierarchy",
    "off_note_and_maillard_suppression":         r"F08 \ Off-notes \& Maillard suppression",
    "carbohydrate_pyrolysis_and_caramelization": r"F09 \ Caramelisation \& pyrolysis",
    "fermentation_pretreatment":                 r"F10 \ Fermentation pre-treatment",
    "lipid_maillard_crosstalk":                  r"F11 \ Maillard-lipid crosstalk",
    "protein_damage_markers":                    r"F12 \ Protein damage markers",
    "polyphenol_amino_capping":                  r"F13 \ Polyphenol-amino capping",
    "ascorbic_acid_maillard":                    r"F14 \ Ascorbic acid Maillard",
    "phospholipid_amine_sink":                   r"F15 \ Phospholipid-amine sink",
    "melanoidin_polymerization":                 r"F16 \ Melanoidin polymerisation",
}

# Short labels for scatter plot legend (space-constrained)
_FAMILY_LEGEND_LABELS: dict[str, str] = {
    "amino_acid_sugar_core":                     r"F01 \ Amino acid + sugar",
    "lipid_oxidation_and_carbonylic_crosstalk":  r"F02 \ Lipid oxidation",
    "thiamine_fragmentation_support":            r"F03 \ Thiamine",
    "nucleotide_and_ribose_support":             r"F04 \ Nucleotide support",
    "glutathione_and_peptide_support":           r"F05 \ Peptide support",
    "alternative_protein_matrix_scope":          r"F06 \ Matrix scope",
    "carbonyl_donor_hierarchy":                  r"F07 \ Donor hierarchy",
    "off_note_and_maillard_suppression":         r"F08 \ Off-notes \& suppression",
    "carbohydrate_pyrolysis_and_caramelization": r"F09 \ Caramelisation",
    "fermentation_pretreatment":                 r"F10 \ Pretreatment",
    "lipid_maillard_crosstalk":                  r"F11 \ Lipid crosstalk",
    "protein_damage_markers":                    r"F12 \ Damage markers",
}

def _bench_label(benchmark_id: str) -> str:
    return benchmark_label(benchmark_id, style="latex")


# ---------------------------------------------------------------------------
# Visual constants
# ---------------------------------------------------------------------------

_FAMILY_COLORS: dict[str, str] = {
    "amino_acid_sugar_core":                     "#0077BB",
    "lipid_oxidation_and_carbonylic_crosstalk":  "#EE7733",
    "thiamine_fragmentation_support":            "#228833",
    "nucleotide_and_ribose_support":             "#AA4499",
    "glutathione_and_peptide_support":           "#66A61E",
    "alternative_protein_matrix_scope":          "#8C564B",
    "carbonyl_donor_hierarchy":                  "#DDCC77",
    "off_note_and_maillard_suppression":         "#CC3311",
    "carbohydrate_pyrolysis_and_caramelization": "#009988",
    "fermentation_pretreatment":                 "#CC79A7",
    "lipid_maillard_crosstalk":                  "#4477AA",
    "protein_damage_markers":                    "#882255",
    "polyphenol_amino_capping":                  "#999999",
    "ascorbic_acid_maillard":                    "#BBBBBB",
    "phospholipid_amine_sink":                   "#A6A6A6",
    "melanoidin_polymerization":                 "#C7C7C7",
}
_FAMILY_COLOR_NONE = "#BBBBBB"

_STATUS_COLORS: dict[str, str] = {
    "pass":            "#2a7f62",
    "pass-no-ranking": "#4c956c",
    "partial-pass":    "#cc8b3b",
    "fail":            "#b03a2e",
}

_PATH_MARKERS: dict[str, str] = {
    "free_precursor":              "o",
    "matrix_only":                 "^",
    "matrix_precursor_augmented":  "s",
}
_PATH_LABELS: dict[str, str] = {
    "free_precursor":              r"Free precursor (strict gate)",
    "matrix_only":                 r"Matrix-only lane",
    "matrix_precursor_augmented":  r"Matrix + precursor augmented",
}

_FS = 12  # base font size for all figures


def _path_marker(execution_path: str) -> str:
    return _PATH_MARKERS.get(execution_path, "o")


def _family_palette(family_ids: list[str]) -> dict[str, str]:
    """Backward-compat palette accessor kept for callers outside this module."""
    return {fid: _FAMILY_COLORS.get(fid, _FAMILY_COLOR_NONE) for fid in family_ids}


# ---------------------------------------------------------------------------
# Figure 1: Concentration parity scatter
# ---------------------------------------------------------------------------

def _render_parity(payload: dict[str, object], output_path: Path) -> None:
    """Standalone parity scatter: predicted vs. measured ppb for all compounds."""
    quantitative_points = payload["quantitative_points"]
    if not quantitative_points:
        raise RuntimeError("No quantitative points available for parity figure.")

    measured = np.array([r["measured_ppb"] for r in quantitative_points], dtype=float)
    predicted = np.array([r["predicted_ppb"] for r in quantitative_points], dtype=float)
    lower = min(measured[measured > 0].min(), predicted[predicted > 0].min()) * 0.4
    upper = max(measured.max(), predicted.max()) * 2.5
    ref = np.geomspace(lower, upper, 300)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.fill_between(ref, ref / 2.0, ref * 2.0, color="#FFF3CD", alpha=0.55, zorder=1)
    ax.fill_between(ref, ref / 1.5, ref * 1.5, color="#D9EAD3", alpha=0.70, zorder=2)
    ax.plot(ref, ref, color="#2F4858", linewidth=1.5, zorder=3)

    family_ids_with_data = sorted({str(r["chemistry_family"]) for r in quantitative_points})
    execution_paths_seen: set[str] = set()

    for fid in family_ids_with_data:
        color = _FAMILY_COLORS.get(fid, _FAMILY_COLOR_NONE)
        subset = [r for r in quantitative_points if r["chemistry_family"] == fid]
        for ep in sorted({r["execution_path"] for r in subset}):
            lane = [r for r in subset if r["execution_path"] == ep]
            x = np.array([r["measured_ppb"] for r in lane], dtype=float)
            y = np.array([r["predicted_ppb"] for r in lane], dtype=float)
            xerr = np.array(
                [r["measured_ppb"] * ((r.get("uncertainty_pct") or 0.0) / 100.0) for r in lane],
                dtype=float,
            )
            ax.errorbar(
                x, y, xerr=xerr,
                fmt=_path_marker(ep), color=color, ecolor=color,
                capsize=3, markersize=7, linewidth=0.8, zorder=4,
            )
            execution_paths_seen.add(ep)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lower, upper)
    ax.set_ylim(lower, upper)
    ax.set_xlabel(r"Measured concentration (ppb)", fontsize=_FS)
    ax.set_ylabel(r"Predicted concentration (ppb)", fontsize=_FS)
    ax.set_title(r"Predicted vs.\ measured concentration --- all matched compounds", fontsize=_FS)
    ax.tick_params(labelsize=_FS - 1)

    # Build legend handles
    family_handles = [
        mpatches.Patch(
            facecolor=_FAMILY_COLORS.get(fid, _FAMILY_COLOR_NONE),
            label=_FAMILY_LEGEND_LABELS.get(fid, _FAMILY_LABELS.get(fid, fid)),
        )
        for fid in family_ids_with_data
    ]
    band_handles = [
        mpatches.Patch(facecolor="#D9EAD3", alpha=0.9, label=r"Within $1.5\times$"),
        mpatches.Patch(facecolor="#FFF3CD", alpha=0.85, label=r"Within $2\times$"),
        plt.Line2D([0], [0], color="#2F4858", linewidth=1.5, label=r"Ideal parity ($y = x$)"),
    ]
    path_handles = [
        plt.Line2D(
            [0], [0], marker=_path_marker(ep), color="0.4",
            linestyle="None", markersize=7,
            label=_PATH_LABELS.get(ep, ep),
        )
        for ep in sorted(execution_paths_seen)
    ]
    guide_legend = ax.legend(
        handles=band_handles + path_handles,
        fontsize=_FS - 4,
        loc="upper left",
        bbox_to_anchor=(0.04, 0.96),
        bbox_transform=ax.transAxes,
        framealpha=0.9,
        edgecolor="0.75",
        title=r"\textit{Guide}",
        title_fontsize=_FS - 4,
        labelspacing=0.12,
        handletextpad=0.3,
        borderpad=0.2,
        handlelength=1.0,
        columnspacing=0.6,
        ncol=1,
    )
    ax.add_artist(guide_legend)
    ax.legend(
        handles=family_handles,
        fontsize=_FS - 4,
        loc="lower right",
        bbox_to_anchor=(0.98, 0.04),
        bbox_transform=ax.transAxes,
        framealpha=0.9,
        edgecolor="0.75",
        title=r"\textit{Chemistry families}",
        title_fontsize=_FS - 4,
        labelspacing=0.12,
        handletextpad=0.3,
        borderpad=0.2,
    )

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2: Per-benchmark worst-case ratio
# ---------------------------------------------------------------------------

def _render_benchmark_accuracy(payload: dict[str, object], output_path: Path) -> None:
    """Standalone bar chart: worst-case predicted/measured ratio per benchmark."""
    bench_summaries = payload["bench_summaries"]
    quantitative = [s for s in bench_summaries if s.supported and s.matched_compounds > 0]
    ranked = sorted(quantitative, key=lambda s: (s.max_ratio or 0.0), reverse=True)

    labels = [_bench_label(s.benchmark_id) for s in ranked]
    ratios = [float(s.max_ratio or 0.0) for s in ranked]
    colors = [_STATUS_COLORS.get(s.overall_status, "#607d8b") for s in ranked]

    fig, ax = plt.subplots(figsize=(8, 6))
    y_pos = np.arange(len(ranked))
    bars = ax.barh(y_pos, ratios, color=colors, edgecolor="white", linewidth=0.6, zorder=3)
    for bar, ratio in zip(bars, ratios):
        ax.text(
            bar.get_width() + 0.05,
            bar.get_y() + bar.get_height() / 2.0,
            f"{ratio:.2f}$\\times$",
            va="center", fontsize=_FS - 2,
        )

    ax.axvline(1.5, color="#8C1C13", linestyle="--", linewidth=1.5, zorder=4)
    ax.axvline(2.0, color="#555555", linestyle=":",  linewidth=1.5, zorder=4)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=_FS - 1)
    ax.set_xlabel(r"Worst-case ratio \ (predicted / measured)", fontsize=_FS)
    ax.set_title(r"Per-benchmark quantitative accuracy", fontsize=_FS)
    ax.tick_params(axis="x", labelsize=_FS - 1)

    legend_handles = [
        mpatches.Patch(facecolor=_STATUS_COLORS["pass"],         label=r"Pass \ ($\leq 1.5\times$)"),
        mpatches.Patch(facecolor=_STATUS_COLORS["partial-pass"], label=r"Partial \ ($\leq 2\times$)"),
        mpatches.Patch(facecolor=_STATUS_COLORS["fail"],         label=r"Fail \ ($> 2\times$)"),
        plt.Line2D([0], [0], color="#8C1C13", linestyle="--", linewidth=1.5, label=r"Strict gate \ (1.5$\times$)"),
        plt.Line2D([0], [0], color="#555555", linestyle=":",  linewidth=1.5, label=r"Matrix tolerance \ (2$\times$)"),
    ]
    # Put benchmark status legend inside the plot area (top-right) with smaller
    # font so it doesn't overlap the bars.
    ax.legend(
        handles=legend_handles,
        fontsize=_FS - 3,
        loc="upper right",
        bbox_to_anchor=(0.98, 0.98),
        bbox_transform=ax.transAxes,
        framealpha=0.9,
        edgecolor="0.75",
        labelspacing=0.2,
    )

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 3: family benchmark coverage census
# ---------------------------------------------------------------------------

def _render_family_coverage(payload: dict[str, object], output_path: Path) -> None:
    """Standalone horizontal bar chart: quantitative compound points per family."""
    families = payload["families"]
    plotted = list(reversed(families))
    family_count = len(families)

    pt_counts = [int(r["quantitative_point_count"]) for r in plotted]
    bar_colors = [_FAMILY_COLORS.get(str(r["chemistry_family"]), _FAMILY_COLOR_NONE) for r in plotted]
    tick_labels = []
    for r in plotted:
        fid = str(r["chemistry_family"])
        pts = int(r["quantitative_point_count"])
        label = _FAMILY_LABELS.get(fid, fid.replace("_", " "))
        if pts > 0:
            tick_labels.append(rf"{label} \ ({pts} compounds)")
        else:
            tick_labels.append(rf"{label}")

    fig, ax = plt.subplots(figsize=(8, 6))
    y_pos = np.arange(len(plotted))
    ax.barh(y_pos, pt_counts, color=bar_colors, edgecolor="white", linewidth=0.6, zorder=3)

    # Annotate gaps explicitly on the bar row
    for idx, r in enumerate(plotted):
        if int(r["quantitative_point_count"]) == 0:
            status_text = r"\textit{Active in simulation}" if bool(r.get("has_runtime_support", False)) else r"\textit{Not yet modeled}"
            ax.text(0.12, idx, status_text, va="center",
                    fontsize=_FS - 3, color="#777777", style="italic")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(tick_labels, fontsize=_FS - 2)
    ax.set_xlabel(r"Number of matched compound points", fontsize=_FS)
    ax.set_title(rf"Benchmark coverage across all {family_count} tracked chemistry families", fontsize=_FS)
    ax.tick_params(axis="x", labelsize=_FS - 1)
    ax.set_xlim(left=0)

    legend_handles = [
        mpatches.Patch(facecolor="#555555", label=r"Reaction family with matched data"),
        mpatches.Patch(facecolor=_FAMILY_COLOR_NONE, label=r"No quantitative data points yet"),
    ]
    ax.legend(
        handles=legend_handles,
        fontsize=_FS - 3,
        loc="lower right",
        bbox_to_anchor=(0.98, 0.02),
        bbox_transform=ax.transAxes,
        framealpha=0.9,
        edgecolor="0.75",
        labelspacing=0.2,
    )



    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Backward-compat wrapper (kept so existing callers don't break)
# ---------------------------------------------------------------------------

def _render_figure(payload: dict[str, object], output_path: Path) -> None:
    """Alias for the parity scatter; kept for any external callers."""
    _render_parity(payload, output_path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    parser.add_argument("--docs-asset-dir", default="docs/assets")
    args = parser.parse_args()

    output_dir = ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    docs_asset_dir = ROOT / args.docs_asset_dir
    docs_asset_dir.mkdir(parents=True, exist_ok=True)

    payload = build_family_validation_overview_artifact()
    payload["bench_summaries"] = summarize_benchmarks()
    markdown = render_family_validation_overview_markdown(payload)

    md_path = output_dir / "family_validation_overview.md"
    json_path = output_dir / "family_validation_overview.json"
    json_payload = {k: v for k, v in payload.items() if k != "bench_summaries"}
    md_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(json_payload, indent=2), encoding="utf-8")

    figures = [
        ("family_parity.png",             _render_parity),
        ("family_benchmark_accuracy.png", _render_benchmark_accuracy),
        ("family_coverage.png",           _render_family_coverage),
    ]
    for fname, render_fn in figures:
        png_path = output_dir / fname
        render_fn(payload, png_path)
        shutil.copyfile(png_path, docs_asset_dir / fname)
        print(f"Wrote {png_path}")
        print(f"Copied {fname} to {docs_asset_dir / fname}")

    # Keep backward-compat copy for any existing README/doc reference
    shutil.copyfile(output_dir / "family_parity.png", output_dir / "family_validation_overview.png")
    shutil.copyfile(output_dir / "family_parity.png", docs_asset_dir / "family_validation_overview.png")

    print(f"Wrote {md_path}")
    print(f"Wrote {json_path}")
    print(markdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())