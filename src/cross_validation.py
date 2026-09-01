"""S20.2 — Leave-one-benchmark-out diagnostics.

We do not (yet) have a programmatic "refit priors on N-1 benchmarks" loop,
so a literal LOO over a learning step is not what we report. Instead, this
module surfaces the **leverage** each benchmark has on the headline trust
metric (90% CI coverage rate from `src.uncertainty_propagation`):

* `coverage_full` = hit-rate across all matched compounds.
* `coverage_minus_b` = hit-rate when benchmark `b` is removed.
* `leverage` = `coverage_full - coverage_minus_b`. Positive ⇒ removing `b`
  *lowers* coverage (so `b` is currently load-bearing — its compounds are
  inside the envelope and removing them would expose model weakness).
  Negative ⇒ removing `b` *raises* coverage, i.e. `b` is the residual
  outlier dragging the panel down — usually the right place to spend the
  next experiment.

Per benchmark we also report `mean_abs_log10_residual` of `predicted_p50`
vs `measured_ppb` so the scientist can rank residual magnitude
independently of whether the wide envelope happens to contain the point.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

from src import data_paths

PREDICTION_UNCERTAINTY_PATH = data_paths.PREDICTION_UNCERTAINTY


@dataclass(frozen=True)
class BenchmarkLeverage:
    benchmark_id: str
    matched_compounds: int
    inside_ci: int
    coverage_self: Optional[float]
    coverage_minus_self: Optional[float]
    leverage: Optional[float]
    mean_abs_log10_residual: Optional[float]
    max_abs_log10_residual: Optional[float]


def _abs_log10_residual(measured: float, p50: float) -> Optional[float]:
    if measured <= 0.0 or p50 <= 0.0:
        return None
    if not (math.isfinite(measured) and math.isfinite(p50)):
        return None
    return abs(math.log10(p50 / measured))


def _aggregate(
    benchmarks: Sequence[Mapping[str, Any]],
    *,
    skip_id: Optional[str] = None,
) -> Dict[str, Any]:
    matched = 0
    inside = 0
    residuals: List[float] = []
    for bench in benchmarks:
        if skip_id is not None and str(bench.get("benchmark_id")) == skip_id:
            continue
        for compound in bench.get("compounds", []) or []:
            matched += 1
            if compound.get("inside_ci"):
                inside += 1
            res = _abs_log10_residual(
                float(compound.get("measured_ppb", 0.0) or 0.0),
                float(compound.get("predicted_p50", 0.0) or 0.0),
            )
            if res is not None:
                residuals.append(res)
    coverage = inside / matched if matched else None
    mean_res = sum(residuals) / len(residuals) if residuals else None
    max_res = max(residuals) if residuals else None
    return {
        "matched": matched,
        "inside_ci": inside,
        "coverage": coverage,
        "mean_abs_log10_residual": mean_res,
        "max_abs_log10_residual": max_res,
    }


def compute_leverage(payload: Mapping[str, Any]) -> Dict[str, Any]:
    """Return panel-level + per-benchmark leverage diagnostics."""
    benchmarks: List[Mapping[str, Any]] = list(payload.get("benchmarks", []) or [])

    full = _aggregate(benchmarks)
    leverage_rows: List[BenchmarkLeverage] = []

    for bench in benchmarks:
        bench_id = str(bench.get("benchmark_id"))
        self_stats = _aggregate([bench])
        rest_stats = _aggregate(benchmarks, skip_id=bench_id)
        leverage_value: Optional[float] = None
        if full["coverage"] is not None and rest_stats["coverage"] is not None:
            leverage_value = float(full["coverage"] - rest_stats["coverage"])
        leverage_rows.append(
            BenchmarkLeverage(
                benchmark_id=bench_id,
                matched_compounds=int(self_stats["matched"]),
                inside_ci=int(self_stats["inside_ci"]),
                coverage_self=self_stats["coverage"],
                coverage_minus_self=rest_stats["coverage"],
                leverage=leverage_value,
                mean_abs_log10_residual=self_stats["mean_abs_log10_residual"],
                max_abs_log10_residual=self_stats["max_abs_log10_residual"],
            )
        )

    leverage_rows.sort(
        key=lambda r: (r.leverage is None, -(r.leverage or 0.0))
    )

    return {
        "panel": {
            "matched_compound_count": full["matched"],
            "inside_ci_count": full["inside_ci"],
            "coverage": full["coverage"],
            "mean_abs_log10_residual": full["mean_abs_log10_residual"],
            "max_abs_log10_residual": full["max_abs_log10_residual"],
        },
        "benchmarks": [
            {
                "benchmark_id": r.benchmark_id,
                "matched_compounds": r.matched_compounds,
                "inside_ci": r.inside_ci,
                "coverage_self": r.coverage_self,
                "coverage_minus_self": r.coverage_minus_self,
                "leverage": r.leverage,
                "mean_abs_log10_residual": r.mean_abs_log10_residual,
                "max_abs_log10_residual": r.max_abs_log10_residual,
            }
            for r in leverage_rows
        ],
    }


def render_markdown(report: Mapping[str, Any]) -> str:
    panel = report.get("panel", {})
    coverage = panel.get("coverage")
    coverage_str = f"{coverage * 100:.1f}%" if coverage is not None else "n/a"
    lines = [
        '<!-- Auto-regenerated by ./scripts/docker_maillard.sh run "python scripts/generators/generate_loo_leverage.py". Manual edits will be overwritten. -->',
        "",
        "# Leave-One-Benchmark-Out Leverage",
        "",
        "_Per-benchmark contribution to the panel's 90% CI coverage rate._",
        "",
        f"Panel coverage: **{coverage_str}** ({panel.get('inside_ci_count', 0)} / {panel.get('matched_compound_count', 0)}).",
        f"Panel mean |log₁₀(P50/measured)| = "
        + (
            f"{float(panel.get('mean_abs_log10_residual', 0.0)):.2f} dex"
            if panel.get("mean_abs_log10_residual") is not None
            else "n/a"
        )
        + ".",
        "",
        "Sorted by descending leverage (positive ⇒ benchmark currently *carries* coverage; negative ⇒ benchmark *drags* it down and is the natural next experiment).",
        "",
        "| Benchmark | Matched | Inside CI | Self coverage | Coverage minus self | Leverage | Mean |Δlog₁₀| | Max |Δlog₁₀| |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for entry in report.get("benchmarks", []) or []:
        cov_self = entry.get("coverage_self")
        cov_rest = entry.get("coverage_minus_self")
        lev = entry.get("leverage")
        mean_res = entry.get("mean_abs_log10_residual")
        max_res = entry.get("max_abs_log10_residual")
        lines.append(
            "| `{bid}` | {matched} | {inside} | {cov_self} | {cov_rest} | {lev} | {mean_res} | {max_res} |".format(
                bid=entry.get("benchmark_id"),
                matched=entry.get("matched_compounds", 0),
                inside=entry.get("inside_ci", 0),
                cov_self=f"{float(cov_self) * 100:.1f}%" if cov_self is not None else "n/a",
                cov_rest=f"{float(cov_rest) * 100:.1f}%" if cov_rest is not None else "n/a",
                lev=f"{float(lev):+.3f}" if lev is not None else "n/a",
                mean_res=f"{float(mean_res):.2f}" if mean_res is not None else "n/a",
                max_res=f"{float(max_res):.2f}" if max_res is not None else "n/a",
            )
        )
    lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def write_artifact(
    report: Mapping[str, Any],
    *,
    output_dir: Path | str = data_paths.VALIDATION_DIR,
    basename: str = "loo_leverage",
) -> Dict[str, Path]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / f"{basename}.json"
    md_path = output_dir / f"{basename}.md"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    md_path.write_text(render_markdown(report), encoding="utf-8")
    return {"json": json_path, "md": md_path}


def load_prediction_payload(
    path: Path | str = PREDICTION_UNCERTAINTY_PATH,
) -> Dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))
