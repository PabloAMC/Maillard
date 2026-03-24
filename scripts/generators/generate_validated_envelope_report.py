#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import scienceplots

plt.style.use(["science", "no-latex"])

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.usability_reports import build_validated_envelope_report
from src.presentation import render_validated_envelope_markdown


def _render_validated_envelope_figure(report, output_path: Path) -> None:
    strict_count = len(report.strict_ready_benchmarks)
    matrix_only_count = len(report.matrix_only_benchmarks)
    directional_count = max(report.supported_benchmarks - strict_count - matrix_only_count, 0)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6))
    fig.suptitle("Maillard Validated Envelope", fontsize=15)

    ax = axes[0]
    categories = ["strict-ready\nfree-precursor", "directional\nsupported", "matrix-only\nexecutable"]
    values = [strict_count, directional_count, matrix_only_count]
    colors = ["#2a7f62", "#c98f3d", "#8c1c13"]
    bars = ax.bar(categories, values, color=colors)
    ax.set_ylabel("Benchmark count")
    ax.set_title("Current Trust Envelope")
    ax.set_ylim(0, max(max(values), 1) + 1)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.05, str(value), ha="center", va="bottom", fontsize=10)

    ax = axes[1]
    ax.axis("off")
    lines = [
        f"Target tag: {report.target_tag}",
        f"Supported benchmarks: {report.supported_benchmarks}/{report.total_benchmarks}",
        "",
        "Strict-ready benchmarks:",
    ]
    lines.extend(f"- {name}" for name in report.strict_ready_benchmarks)
    if report.matrix_only_benchmarks:
        lines.append("")
        lines.append("Matrix-only executable benchmarks:")
        lines.extend(f"- {name}" for name in report.matrix_only_benchmarks)
    lines.append("")
    lines.append("Main caveats:")
    lines.extend(f"- {warning}" for warning in report.warnings[:3])
    ax.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", fontsize=10)
    ax.set_title("Boundary Conditions")

    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    parser.add_argument("--target-tag", default="meaty")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    report = build_validated_envelope_report(target_tag=args.target_tag)
    markdown = render_validated_envelope_markdown(report)

    markdown_path = output_dir / "validated_envelope.md"
    json_path = output_dir / "validated_envelope.json"
    png_path = output_dir / "validated_envelope.png"
    markdown_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(asdict(report), indent=2), encoding="utf-8")
    _render_validated_envelope_figure(report, png_path)

    print(markdown)
    print(f"Wrote {png_path}")
    print(f"Wrote {markdown_path}")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())