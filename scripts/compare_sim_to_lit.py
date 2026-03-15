import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import evaluate_benchmark


def compare(lit_path: str):
    evaluation = evaluate_benchmark(lit_path)
    if not evaluation.supported:
        print(f"Benchmark not supported yet: {evaluation.reason}")
        return 1

    compounds = [comparison.compound for comparison in evaluation.comparisons]
    y_exp = np.array([comparison.measured_ppb for comparison in evaluation.comparisons], dtype=float)
    y_sim = np.array([comparison.predicted_ppb for comparison in evaluation.comparisons], dtype=float)

    print("\nComparison Results:")
    print("-" * 72)
    for comparison in evaluation.comparisons:
        exp_display = f"{comparison.measured_ppb:.4g}"
        pred_display = f"{comparison.predicted_ppb:.4g}"
        print(
            f"{comparison.compound:30} | "
            f"Exp: {exp_display:>10} ppb | "
            f"Pred: {pred_display:>10} ppb | "
            f"Match: {comparison.matched_name or '-'}"
        )

    print("-" * 72)
    print(f"Coverage:  {evaluation.coverage:.2%}")
    print(
        f"Pearson R: {evaluation.pearson_r:.4f}"
        if evaluation.pearson_r is not None
        else "Pearson R: n/a (requires at least 3 matched compounds)"
    )
    print(f"MAE:       {evaluation.mae_ppb:.2f} ppb" if evaluation.mae_ppb is not None else "MAE: n/a")

    out_dir = Path("results/validation")
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    x = np.arange(len(compounds))
    width = 0.35

    ax1.bar(x - width / 2, y_exp, width, label="Exp (Lit)")
    ax1.bar(x + width / 2, y_sim, width, label="Pred (Framework)")
    ax1.set_ylabel("Concentration (ppb)")
    title = f"Absolute Yields ({evaluation.benchmark_id})"
    if evaluation.pearson_r is not None:
        title += f"\nR = {evaluation.pearson_r:.3f}"
    else:
        title += "\nR = n/a (< 3 matches)"
    ax1.set_title(title)
    ax1.set_xticks(x)
    ax1.set_xticklabels(compounds, rotation=45, ha="right")
    ax1.legend()

    y_exp_norm = y_exp / np.max(y_exp) if np.max(y_exp) > 0 else y_exp
    y_sim_norm = y_sim / np.max(y_sim) if np.max(y_sim) > 0 else y_sim
    ax2.bar(x - width / 2, y_exp_norm, width, label="Exp (Normalized)")
    ax2.bar(x + width / 2, y_sim_norm, width, label="Pred (Normalized)")
    ax2.set_ylabel("Relative abundance")
    ax2.set_title("Normalized abundance profile")
    ax2.set_xticks(x)
    ax2.set_xticklabels(compounds, rotation=45, ha="right")
    ax2.legend()

    plt.tight_layout()
    plot_path = out_dir / f"{evaluation.benchmark_id}_comparison.png"
    plt.savefig(plot_path)
    print(f"\nPlot saved to {plot_path}")
    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--lit", default="data/benchmarks/cys_ribose_140C_Hofmann1998.json")
    args = parser.parse_args()

    raise SystemExit(compare(args.lit))
