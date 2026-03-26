#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.objective_progress import (  # noqa: E402
    build_objective_progress_artifact,
    render_objective_progress_markdown,
)


def _render_objective_closure_panel(payload: dict[str, object], output_path: Path) -> None:
    objectives = payload.get("objectives", [])

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    fig.subplots_adjust(left=0.28, right=0.96, top=0.88, bottom=0.18)

    labels = [str(row.get("label", "unknown")) for row in objectives]
    closed = np.array([int(row.get("closed_count", 0)) for row in objectives], dtype=float)
    target = np.array([int(row.get("target_count", 0)) for row in objectives], dtype=float)
    remaining = np.maximum(target - closed, 0.0)
    positions = np.arange(len(labels))

    ax.barh(positions, closed, color="#2a7f62", label="closed")
    ax.barh(positions, remaining, left=closed, color="#d9d9d9", label="remaining")
    ax.set_yticks(positions)
    ax.set_yticklabels(labels)
    ax.set_xlabel("tracked units")
    ax.set_title("Objective closure")
    ax.legend(frameon=False, fontsize=8)
    for idx, (closed_value, target_value) in enumerate(zip(closed, target)):
        ax.text(max(target_value, 0.0) + 0.05, idx, f"{int(closed_value)}/{int(target_value)}", va="center", fontsize=8)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def _render_hexanal_nonanal_panel(payload: dict[str, object], output_path: Path) -> None:
    prediction_rows = payload.get("hexanal_nonanal_prediction_change", [])

    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    fig.subplots_adjust(left=0.12, right=0.97, top=0.88, bottom=0.28)

    ratio_positions = np.arange(len(prediction_rows))
    ratio_values = np.array([float(row.get("ratio", 0.0) or 0.0) for row in prediction_rows], dtype=float)
    ratio_labels = [f"{row.get('protein_type', 'unknown')} {row.get('compound', 'unknown')}" for row in prediction_rows]
    ax.axhspan(0.5, 2.0, color="#eaf4ea", alpha=0.8, label="accepted band")
    ax.axhline(1.0, color="#2f4858", linestyle="--", linewidth=1.2, label="ideal parity")
    ax.bar(ratio_positions, ratio_values, color="#3f88c5", width=0.6)
    ax.set_xticks(ratio_positions)
    ax.set_xticklabels(ratio_labels, rotation=20, ha="right")
    ax.set_ylabel("ProtocolPilot2026 / Internal2026")
    ax.set_title("Hexanal/Nonanal prediction change")
    ax.set_ylim(0.0, max(2.2, float(np.max(ratio_values)) + 0.2 if len(ratio_values) else 2.2))
    ax.legend(frameon=False, fontsize=8)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def _render_figure(payload: dict[str, object], output_path: Path) -> None:
    objectives = payload.get("objectives", [])
    prediction_rows = payload.get("hexanal_nonanal_prediction_change", [])

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    fig.subplots_adjust(left=0.08, right=0.97, top=0.90, bottom=0.16, wspace=0.32)
    fig.suptitle("Objective Progress and Prediction Impact", fontsize=12)

    labels = [str(row.get("label", "unknown")) for row in objectives]
    closed = np.array([int(row.get("closed_count", 0)) for row in objectives], dtype=float)
    target = np.array([int(row.get("target_count", 0)) for row in objectives], dtype=float)
    remaining = np.maximum(target - closed, 0.0)
    positions = np.arange(len(labels))

    axes[0].barh(positions, closed, color="#2a7f62", label="closed")
    axes[0].barh(positions, remaining, left=closed, color="#d9d9d9", label="remaining")
    axes[0].set_yticks(positions)
    axes[0].set_yticklabels(labels)
    axes[0].set_xlabel("tracked units")
    axes[0].set_title("Objective closure")
    axes[0].legend(frameon=False, fontsize=8)
    for idx, (closed_value, target_value) in enumerate(zip(closed, target)):
        axes[0].text(max(target_value, 0.0) + 0.05, idx, f"{int(closed_value)}/{int(target_value)}", va="center", fontsize=8)

    ratio_positions = np.arange(len(prediction_rows))
    ratio_values = np.array([float(row.get("ratio", 0.0) or 0.0) for row in prediction_rows], dtype=float)
    ratio_labels = [f"{row.get('protein_type', 'unknown')} {row.get('compound', 'unknown')}" for row in prediction_rows]
    axes[1].axhspan(0.5, 2.0, color="#eaf4ea", alpha=0.8, label="accepted band")
    axes[1].axhline(1.0, color="#2f4858", linestyle="--", linewidth=1.2, label="ideal parity")
    axes[1].bar(ratio_positions, ratio_values, color="#3f88c5", width=0.6)
    axes[1].set_xticks(ratio_positions)
    axes[1].set_xticklabels(ratio_labels, rotation=20, ha="right")
    axes[1].set_ylabel("ProtocolPilot2026 / Internal2026")
    axes[1].set_title("Hexanal/Nonanal prediction change")
    axes[1].set_ylim(0.0, max(2.2, float(np.max(ratio_values)) + 0.2 if len(ratio_values) else 2.2))
    axes[1].legend(frameon=False, fontsize=8)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    parser.add_argument("--docs-asset-dir", default="docs/assets")
    args = parser.parse_args()

    payload = build_objective_progress_artifact()
    output_dir = ROOT / args.output_dir
    docs_asset_dir = ROOT / args.docs_asset_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    docs_asset_dir.mkdir(parents=True, exist_ok=True)

    (output_dir / "objective_progress.md").write_text(
        render_objective_progress_markdown(payload),
        encoding="utf-8",
    )
    (output_dir / "objective_progress.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    _render_figure(payload, docs_asset_dir / "objective_progress.png")
    _render_objective_closure_panel(payload, docs_asset_dir / "objective_closure.png")
    _render_hexanal_nonanal_panel(payload, docs_asset_dir / "hexanal_nonanal_prediction_change.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())