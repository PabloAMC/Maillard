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

from src.family_barrier_progress import (  # noqa: E402
    build_family_barrier_progress_artifact,
    render_family_barrier_progress_markdown,
)


def _render_figure(payload: dict[str, object], output_path: Path) -> None:
    rows = payload.get("families", [])
    labels = [str(row.get("slr_family", "99")) for row in rows]
    display_names = [str(row.get("display_name", "unknown")) for row in rows]
    explicit_fast = np.array([int(row.get("explicit_fast_barrier_count", 0)) for row in rows], dtype=float)
    arrhenius = np.array([int(row.get("arrhenius_anchor_count", 0)) for row in rows], dtype=float)
    quantitative = np.array([int(row.get("quantitative_point_count", 0)) for row in rows], dtype=float)
    max_ratios = np.array([
        float(row.get("max_ratio")) if row.get("max_ratio") is not None else np.nan for row in rows
    ], dtype=float)
    pos = np.arange(len(rows))

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.subplots_adjust(left=0.15, right=0.95, top=0.92, bottom=0.15)
    fig.suptitle("Family Barrier Implementation Progress", fontsize=12)

    ax.barh(pos, explicit_fast, color="#8fb996", label="explicit FAST lanes")
    ax.barh(pos, arrhenius, color="#2f6690", alpha=0.6, label="Arrhenius anchors")
    ax.set_yticks(pos)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("lane count")
    ax.set_title("Kinetic barrier support by family")
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    
    for idx, name in enumerate(display_names):
        max_val = max(explicit_fast[idx], arrhenius[idx])
        ax.text(max_val + 0.2, idx, name, va="center", fontsize=8, color="#333333")
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    parser.add_argument("--docs-asset-dir", default="docs/assets")
    args = parser.parse_args()

    payload = build_family_barrier_progress_artifact()
    output_dir = ROOT / args.output_dir
    docs_asset_dir = ROOT / args.docs_asset_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    docs_asset_dir.mkdir(parents=True, exist_ok=True)

    (output_dir / "family_barrier_progress.md").write_text(
        render_family_barrier_progress_markdown(payload),
        encoding="utf-8",
    )
    (output_dir / "family_barrier_progress.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    _render_figure(payload, docs_asset_dir / "family_barrier_progress.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())