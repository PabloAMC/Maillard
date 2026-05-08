"""S22.3 — render a benchmark × compound gap heatmap.

Reads `results/validation/experiment_value_ranking.json` (S22.1 output) and
draws a PNG heatmap where:

* rows are benchmarks (sorted by max VoI desc),
* columns are compounds (sorted by max VoI desc),
* cell colour = VoI score (yellow→red),
* cell annotation = "*" if measured value sits outside the 90% CI envelope.

The artifact lives at `results/validation/gap_heatmap.png` and is embedded
into the README "Where the next experiments matter most" section.

Run inside the Docker container:

    ./scripts/docker_maillard.sh run "python scripts/generators/generate_gap_heatmap.py"
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "results" / "validation" / "experiment_value_ranking.json"
DEFAULT_OUTPUT = ROOT / "results" / "validation" / "gap_heatmap.png"


def _short_compound(name: str) -> str:
    name = name.split("(")[0].strip()
    if len(name) > 32:
        name = name[:30] + "…"
    return name


def _short_benchmark(name: str) -> str:
    if len(name) > 38:
        name = name[:36] + "…"
    return name


def build_grid(payload: dict) -> Tuple[List[str], List[str], np.ndarray, np.ndarray]:
    """Return (benchmark_labels, compound_labels, voi_grid, outside_grid)."""
    candidates = payload.get("candidates", [])
    if not candidates:
        raise ValueError("experiment_value_ranking.json has no candidates")

    bench_max: Dict[str, float] = {}
    comp_max: Dict[str, float] = {}
    for cand in candidates:
        b = str(cand["benchmark_id"])
        c = _short_compound(str(cand["compound"]))
        v = float(cand.get("voi_score", 0.0))
        bench_max[b] = max(bench_max.get(b, 0.0), v)
        comp_max[c] = max(comp_max.get(c, 0.0), v)

    bench_order = [b for b, _ in sorted(bench_max.items(), key=lambda kv: -kv[1])]
    comp_order = [c for c, _ in sorted(comp_max.items(), key=lambda kv: -kv[1])]

    voi = np.full((len(bench_order), len(comp_order)), np.nan)
    outside = np.zeros_like(voi, dtype=bool)
    bidx = {b: i for i, b in enumerate(bench_order)}
    cidx = {c: i for i, c in enumerate(comp_order)}
    for cand in candidates:
        b = str(cand["benchmark_id"])
        c = _short_compound(str(cand["compound"]))
        i, j = bidx[b], cidx[c]
        voi[i, j] = float(cand.get("voi_score", 0.0))
        if not bool(cand.get("inside_ci", True)):
            outside[i, j] = True
    return bench_order, comp_order, voi, outside


def render(
    bench_labels: List[str],
    comp_labels: List[str],
    voi: np.ndarray,
    outside: np.ndarray,
    output: Path,
) -> Path:
    n_rows, n_cols = voi.shape
    fig_w = max(7.5, 0.45 * n_cols + 4.5)
    fig_h = max(4.5, 0.32 * n_rows + 2.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    cmap = plt.get_cmap("YlOrRd")
    cmap.set_bad(color="#eeeeee")
    masked = np.ma.masked_invalid(voi)
    vmax = float(np.nanmax(voi)) if np.isfinite(np.nanmax(voi)) else 1.0
    im = ax.imshow(masked, cmap=cmap, vmin=0.0, vmax=vmax, aspect="auto")

    ax.set_xticks(np.arange(n_cols))
    ax.set_yticks(np.arange(n_rows))
    ax.set_xticklabels([_short_compound(c) for c in comp_labels], rotation=55, ha="right", fontsize=8)
    ax.set_yticklabels([_short_benchmark(b) for b in bench_labels], fontsize=8)

    for i in range(n_rows):
        for j in range(n_cols):
            if outside[i, j]:
                ax.text(j, i, "*", ha="center", va="center",
                        color="black", fontsize=9, fontweight="bold")

    ax.set_title(
        "Benchmark × compound experiment-value gaps\n"
        "Cell colour: VoI score (yellow → red).  '*' marks measured value outside 90% CI.",
        fontsize=10,
    )
    cbar = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label("VoI score", fontsize=9)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160)
    plt.close(fig)
    return output


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)

    payload = json.loads(args.input.read_text(encoding="utf-8"))
    bench_labels, comp_labels, voi, outside = build_grid(payload)
    output = render(bench_labels, comp_labels, voi, outside, args.output)
    finite_voi = voi[np.isfinite(voi)]
    print(
        f"wrote {output.relative_to(ROOT)} | "
        f"{len(bench_labels)} benchmarks × {len(comp_labels)} compounds | "
        f"max VoI={float(np.nanmax(finite_voi)):.2f} | "
        f"outside-CI cells={int(outside.sum())}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
