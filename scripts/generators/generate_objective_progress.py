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

from src import data_paths  # noqa: E402
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


def _render_figure(payload: dict[str, object], output_path: Path) -> None:
    """Single-panel figure (the Pilot/Internal ratio panel was removed 2026-09-01: it plotted
    a benchmark against a byte-identical copy of itself)."""
    _render_objective_closure_panel(payload, output_path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default=data_paths.rel(data_paths.VALIDATION_DIR))
    parser.add_argument("--docs-asset-dir", default=data_paths.rel(data_paths.DOCS_ASSETS_DIR))
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())