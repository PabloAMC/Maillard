#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scienceplots

plt.style.use(["science", "no-latex"])

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.benchmark_validation import evaluate_benchmark, summarize_evaluation


def _infer_protein_type(benchmark_id: str) -> str:
    if "pea_isolate" in benchmark_id:
        return "pea_iso"
    if "soy_isolate" in benchmark_id:
        return "soy_iso"
    return "free"


def compare(lit_path: str) -> int:
    evaluation = evaluate_benchmark(lit_path)
    if not evaluation.supported:
        print(f"Benchmark not supported yet: {evaluation.reason}")
        return 1

    summary = summarize_evaluation(
        evaluation,
        protein_type=_infer_protein_type(evaluation.benchmark_id),
    )

    compounds = [comparison.compound for comparison in evaluation.comparisons]
    measured = np.array([comparison.measured_ppb for comparison in evaluation.comparisons], dtype=float)
    predicted = np.array([comparison.predicted_ppb for comparison in evaluation.comparisons], dtype=float)
    uncertainty_pct = np.array([
        comparison.uncertainty_pct if comparison.uncertainty_pct is not None else 0.0
        for comparison in evaluation.comparisons
    ], dtype=float)
    measured_err = measured * uncertainty_pct / 100.0
    fold_change = np.divide(predicted, measured, out=np.full_like(predicted, np.nan), where=measured > 0.0)
    acceptable_low = 1.0 / 1.5
    acceptable_high = 1.5

    print("\nComparison Results:")
    print("-" * 84)
    for comparison, fold in zip(evaluation.comparisons, fold_change):
        exp_display = f"{comparison.measured_ppb:.4g}"
        pred_display = f"{comparison.predicted_ppb:.4g}"
        fold_display = f"{fold:.3f}" if np.isfinite(fold) else "n/a"
        print(
            f"{comparison.compound:30} | "
            f"Exp: {exp_display:>10} ppb | "
            f"Pred: {pred_display:>10} ppb | "
            f"Pred/Exp: {fold_display:>6} | "
            f"Match: {comparison.matched_name or '-'}"
        )

    print("-" * 84)
    print(f"Coverage:      {evaluation.coverage:.2%}")
    print(
        f"Pearson R:     {evaluation.pearson_r:.4f}"
        if evaluation.pearson_r is not None
        else "Pearson R:     n/a (requires at least 3 matched compounds)"
    )
    print(f"MAE:           {evaluation.mae_ppb:.2f} ppb" if evaluation.mae_ppb is not None else "MAE:           n/a")
    print(f"Max ratio:     {summary.max_ratio:.3f}" if summary.max_ratio is not None else "Max ratio:     n/a")
    print(f"Overall status:{summary.overall_status:>16}")
    print(f"Strict-ready:  {'yes' if summary.strict_ready else 'no'}")

    out_dir = Path("results/validation")
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.5))
    fig.suptitle(f"Benchmark Comparison: {evaluation.benchmark_id}", fontsize=14)

    ax = axes[0]
    positive_min = min(np.min(measured[measured > 0.0]), np.min(predicted[predicted > 0.0]))
    positive_max = max(np.max(measured), np.max(predicted))
    lower = positive_min * 0.7
    upper = positive_max * 1.4
    reference_line = np.geomspace(lower, upper, 200)
    ax.fill_between(reference_line, reference_line * acceptable_low, reference_line * acceptable_high, color="#d9ead3", alpha=0.7, label="within 1.5x")
    ax.plot(reference_line, reference_line, color="#2f4858", linewidth=1.5, label="ideal parity")
    ax.errorbar(measured, predicted, xerr=measured_err, fmt="o", color="#b03a2e", ecolor="#6b7280", capsize=3)
    for compound, exp_value, pred_value in zip(compounds, measured, predicted):
        ax.annotate(compound, (exp_value, pred_value), textcoords="offset points", xytext=(4, 4), fontsize=8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lower, upper)
    ax.set_ylim(lower, upper)
    ax.set_xlabel("Measured concentration (ppb)")
    ax.set_ylabel("Predicted concentration (ppb)")
    ax.set_title("Parity Plot")
    ax.legend(loc="lower right")

    ax = axes[1]
    positions = np.arange(len(compounds))
    width = 0.36
    ax.bar(positions - width / 2, measured, width, label="Measured", color="#2a7f62")
    ax.bar(positions + width / 2, predicted, width, label="Predicted", color="#b03a2e")
    ax.errorbar(positions - width / 2, measured, yerr=measured_err, fmt="none", ecolor="#1f2933", capsize=3)
    ax.set_yscale("log")
    ax.set_ylabel("Concentration (ppb, log scale)")
    ax.set_title("Absolute Yields")
    ax.set_xticks(positions)
    ax.set_xticklabels(compounds, rotation=35, ha="right")
    ax.legend()

    ax = axes[2]
    ax.axis("off")
    pearson_text = f"{evaluation.pearson_r:.3f}" if evaluation.pearson_r is not None else "n/a"
    mae_text = f"{evaluation.mae_ppb:.2f} ppb" if evaluation.mae_ppb is not None else "n/a"
    max_ratio_text = f"{summary.max_ratio:.3f}" if summary.max_ratio is not None else "n/a"
    lines = [
        f"Tier: {summary.tier}",
        f"Execution path: {summary.execution_path}",
        f"Reference origin: {summary.reference_signal_origin}",
        f"Coverage: {summary.coverage:.1%}",
        f"Pearson R: {pearson_text}",
        f"MAE: {mae_text}",
        f"Max ratio: {max_ratio_text}",
        f"Overall status: {summary.overall_status}",
        f"Strict-ready: {'yes' if summary.strict_ready else 'no'}",
    ]
    if summary.blocking_issues:
        lines.append("")
        lines.append("Blocking issues:")
        lines.extend(f"- {issue}" for issue in summary.blocking_issues)
    ax.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", fontsize=10)

    fig.tight_layout()
    plot_path = out_dir / f"{evaluation.benchmark_id}_comparison.png"
    markdown_path = out_dir / f"{evaluation.benchmark_id}_comparison.md"
    json_path = out_dir / f"{evaluation.benchmark_id}_comparison.json"
    fig.savefig(plot_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    markdown_path.write_text(
        "\n".join([
            f"# Benchmark Comparison: {evaluation.benchmark_id}",
            "",
            f"- Coverage: {summary.coverage:.1%}",
            f"- Pearson R: {pearson_text}",
            f"- MAE: {mae_text}",
            f"- Max ratio: {max_ratio_text}",
            f"- Overall status: {summary.overall_status}",
            f"- Strict-ready: {'yes' if summary.strict_ready else 'no'}",
        ]) + "\n",
        encoding="utf-8",
    )
    json_path.write_text(
        json.dumps(
            {
                "benchmark_id": evaluation.benchmark_id,
                "coverage": summary.coverage,
                "pearson_r": evaluation.pearson_r,
                "mae_ppb": evaluation.mae_ppb,
                "max_ratio": summary.max_ratio,
                "overall_status": summary.overall_status,
                "strict_ready": summary.strict_ready,
                "comparisons": [
                    {
                        "compound": comparison.compound,
                        "measured_ppb": comparison.measured_ppb,
                        "predicted_ppb": comparison.predicted_ppb,
                        "uncertainty_pct": comparison.uncertainty_pct,
                        "predicted_over_measured": fold,
                    }
                    for comparison, fold in zip(evaluation.comparisons, fold_change)
                ],
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    print(f"\nPlot saved to {plot_path}")
    print(f"Summary saved to {markdown_path}")
    print(f"Data saved to {json_path}")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--lit", default="data/benchmarks/cys_ribose_140C_Hofmann1998.json")
    args = parser.parse_args()

    raise SystemExit(compare(args.lit))
