#!/usr/bin/env python
"""
SLICE PROFILES OF THE B8 SULFUR OBJECTIVE AROUND ITS FROZEN OPTIMUM (step 4, 2026-09-03).

The Laplace covariance (``generate_kinetic_core_b8_laplace.py``) is a local Gaussian at an
optimum that sits on two declared bounds. This script asks, coordinate by coordinate,
whether that Gaussian is a fair picture: it evaluates B8's own cost
(``0.5 * r . r`` over the 62 sigma-weighted rows) along each free coordinate at
``x_opt + k * sigma`` for k in -3..3, the other 22 coordinates held at the optimum, and
compares the measured rise with the quadratic the Jacobian predicts.

It is a SLICE (conditional) profile, not a full profile likelihood: no re-optimisation of
the other coordinates at each point. That is the honest scope at ~4 s per evaluation
(207 evaluations, parallelised over coordinates); a re-optimising profile is ~500x more
and is the backlog's next step if a slice looks non-quadratic.

Per coordinate the artifact reports: the step sigma used (the Laplace sigma when
identified, else a tenth of the declared band), the cost rise at each k, the quadratic
prediction ``0.5 * (J^T J)_ii * (k sigma)^2``, the left/right asymmetry at 2 sigma, the
location of the slice minimum, how many points were clipped by a bound, and a verdict:

    quadratic      rise within a factor 2 of the prediction on both sides
    asymmetric     one side rises > 2x the other (a skewed or one-sided constraint)
    flat           rise below 0.5 (one sigma) out to 3 sigma on both sides
    bound_limited  the band ends before 2 sigma on a side
    shifted        the slice minimum is more than one sigma from the frozen value

Nothing here refits; the frozen report is not modified.

Usage:
    python scripts/generators/generate_kinetic_core_b8_profile.py [--workers 6]
"""
from __future__ import annotations

import argparse
import json
import multiprocessing
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "generators"))

from generate_kinetic_core_b8_laplace import frozen_vector, wave_module  # noqa: E402
from src import data_paths  # noqa: E402

B8 = None  # bound in main() to the wave's fit module
OUT = data_paths.VALIDATION_DIR / "kinetic_core_b8_profile.json"
KS = (-3, -2, -1, 0, 1, 2, 3)


def _cost(x: np.ndarray) -> float:
    r = B8.residual_vector(x, True)
    return 0.5 * float(r @ r)


def _slice(args):
    j, i, key, kind, x_opt, lo, hi, sigma, jtj_ii = args
    points: List[Dict[str, Any]] = []
    for k in KS:
        x = x_opt.copy()
        value = float(x_opt[i]) + k * sigma
        clipped = value < lo or value > hi
        value = min(max(value, lo), hi)
        x[i] = value
        points.append({"k": k, "value": value, "clipped": bool(clipped), "cost": _cost(x)})
    return j, key, kind, sigma, points


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--output", default=None)
    parser.add_argument("--wave", default="b9")
    args = parser.parse_args(argv)
    global B8
    B8 = wave_module(args.wave)
    if not Path(B8.OUT_FIT_REPORT).exists() and args.wave != "b8":
        B8 = wave_module("b8")
    wave = B8.WAVE.lower() if hasattr(B8, "WAVE") else "b8"
    out = Path(args.output) if args.output else data_paths.VALIDATION_DIR / f"kinetic_core_{wave}_profile.json"
    lap = json.loads((data_paths.VALIDATION_DIR / f"kinetic_core_{wave}_laplace_covariance.json").read_text(encoding="utf-8"))
    report = json.loads(Path(B8.OUT_FIT_REPORT).read_text(encoding="utf-8"))
    x_opt = frozen_vector(report)
    lower, upper = B8.full_bounds()
    idx = list(B8.FREE_INDEX)
    keys = list(B8.FREE_KEYS)
    # J^T J diagonal from the Laplace artifact: |J_j|^2 is exactly (JtJ)_jj.
    jtj_diag = [float(v) ** 2 for v in lap["jacobian_column_norms"]]

    t0 = time.time()
    cost0 = _cost(x_opt)
    tasks = []
    for j, i in enumerate(idx):
        lo, hi = float(lower[i]), float(upper[i])
        sigma = float(lap["sigma"][j])
        if not (np.isfinite(sigma) and lap["identified"][j]):
            sigma = 0.1 * (hi - lo)
        tasks.append((j, i, keys[j], lap["coordinates"][j]["kind"], x_opt, lo, hi, sigma, jtj_diag[j]))
    if args.workers > 1:
        with multiprocessing.get_context("fork").Pool(args.workers) as pool:
            results = pool.map(_slice, tasks, chunksize=1)
    else:
        results = [_slice(t) for t in tasks]

    rows = []
    verdict_counts: Dict[str, int] = {}
    for (j, key, kind, sigma, points), jj in zip(sorted(results), jtj_diag):
        by_k = {p["k"]: p for p in points}
        rise = {k: by_k[k]["cost"] - cost0 for k in KS}
        quad = {k: 0.5 * jj * (k * sigma) ** 2 for k in KS}
        left2, right2 = rise[-2], rise[2]
        clipped_left = any(by_k[k]["clipped"] for k in (-1, -2))
        clipped_right = any(by_k[k]["clipped"] for k in (1, 2))
        min_k = min(KS, key=lambda k: by_k[k]["cost"])
        if clipped_left or clipped_right:
            verdict = "bound_limited"
        elif abs(min_k) >= 1 and by_k[min_k]["cost"] < cost0 - 0.05:
            verdict = "shifted"
        elif max(rise[-3], rise[3]) < 0.5:
            verdict = "flat"
        elif (left2 > 0 and right2 > 0) and (max(left2, right2) / max(min(left2, right2), 1e-9) > 2.0):
            verdict = "asymmetric"
        else:
            ratio = max(rise[2], 1e-9) / max(quad[2], 1e-9)
            verdict = "quadratic" if 0.5 <= ratio <= 2.0 else "asymmetric"
        verdict_counts[verdict] = verdict_counts.get(verdict, 0) + 1
        rows.append(
            {
                "key": key, "kind": kind, "sigma_used": sigma,
                "identified_in_laplace": bool(lap["identified"][j]),
                "cost_rise": {str(k): rise[k] for k in KS},
                "quadratic_prediction": {str(k): quad[k] for k in KS},
                "clipped": {str(k): by_k[k]["clipped"] for k in KS},
                "slice_min_k": min_k,
                "asymmetry_2sigma": (max(left2, right2) / max(min(left2, right2), 1e-9)) if min(left2, right2) > 0 else None,
                "verdict": verdict,
            }
        )
    payload = {
        "artifact": "kinetic_core_b8_profile",
        "generated_on": time.strftime("%Y-%m-%d"),
        "method": "slice profile (other coordinates fixed at the optimum), quick pH fixed point, k in -3..3",
        "cost_at_optimum": cost0,
        "n_evaluations": len(KS) * len(keys),
        "wall_seconds": round(time.time() - t0, 1),
        "verdict_counts": dict(sorted(verdict_counts.items())),
        "rows": rows,
        "reading": (
            "A 'quadratic' slice supports the Laplace picture on that coordinate; 'asymmetric' "
            "means the Gaussian over-states one side; 'flat' means the data do not constrain it "
            "within 3 sigma; 'bound_limited' means the declared band, not the data, ends the slice; "
            "'shifted' means the frozen vector is not at the slice minimum (quick vs careful "
            "integration, or a converged-on-a-bound optimum)."
        ),
    }
    payload["wave"] = wave
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    md = [
        "# B8 sulfur objective: slice profiles around the frozen optimum", "",
        f"cost at optimum {cost0:.3f}; {payload['n_evaluations']} evaluations; verdicts {payload['verdict_counts']}", "",
        "| coordinate | kind | sigma used | rise -2s | rise +2s | quadratic +2s | asymmetry | slice min k | verdict |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        md.append(
            f"| {r['key']} | {r['kind']} | {r['sigma_used']:.3g} | {r['cost_rise']['-2']:.3f} | {r['cost_rise']['2']:.3f} | "
            f"{r['quadratic_prediction']['2']:.3f} | {'-' if r['asymmetry_2sigma'] is None else f'{r['asymmetry_2sigma']:.2f}'} | "
            f"{r['slice_min_k']} | **{r['verdict']}** |"
        )
    out.with_suffix(".md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"wrote {out}: {payload['verdict_counts']} wall={payload['wall_seconds']}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
