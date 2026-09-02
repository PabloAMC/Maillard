#!/usr/bin/env python
"""
A LAPLACE (GAUSS-NEWTON) COVARIANCE AT THE B8 SULFUR OPTIMUM (step B8, 2026-09-03).

The B8 fit report (``results/validation/kinetic_core_b8_fit_report.json``) freezes
23 free sulfur coordinates but carries no parameter uncertainty, so the core's
Monte-Carlo envelope (``src/kinetic_core/uncertainty.py``) samples NOTHING on the
sulfur lane and 24 of its 49 panel rows are "not evaluable". This script computes
what B2.3's ``intervals()`` computed for the earlier sulfur waves and B8 did not:

    J      = d r / d x   (central finite differences over the 23 FREE coordinates,
                          one-sided at a bound; r is B8's own sigma-weighted residual
                          vector over its 62 rows, at the frozen optimum)
    chi2   = (r . r) / (n_rows - n_free)
    Sigma  = pinv(J^T J) * chi2

and writes ``results/validation/kinetic_core_b8_laplace_covariance.json``: the
free keys, their fit-report coordinates, the optimum, the step sizes, Sigma, the
marginal sigmas, which directions are IDENTIFIED (sigma below the same thresholds
B2.3 used: 3 dex for a log10 k, 60 kJ/mol for an Ea, 0.5 for the acid yield, 1.5
for the pKa), the eigen-spectrum, and the residual norm the frozen vector
actually scores (quick and careful integration). The envelope samples the
IDENTIFIED sub-space from Sigma and clips every draw to the declared bounds;
unidentified directions stay at the optimum and are listed as such -- a flat
direction is not a wide prior, it is no information.

Nothing here refits anything: the optimum is read from the frozen report and the
report is not modified.

Usage:
    python scripts/generators/generate_kinetic_core_b8_laplace.py [--careful]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "generators"))

import generate_kinetic_core_b2_3_fit as B23  # noqa: E402
import generate_kinetic_core_b8_fit as B8  # noqa: E402
from src import data_paths  # noqa: E402

OUT = data_paths.KINETIC_CORE_B8_LAPLACE
THRESHOLDS = {"Ea": 60.0, "log10k": 3.0, "acid_yield": 0.5, "pka": 1.5}
STEPS = {"Ea": 2.0, "log10k": 0.05, "acid_yield": 0.01, "pka": 0.1}


def coordinate_of(key: str) -> Dict[str, str]:
    """Where a free key lives in the fit report's ``frozen_parameters``."""
    if key in B23.PARAM_ORDER:
        return {"block": "log10_k_ref_at_145C", "key": key, "kind": "log10k"}
    if key == "Ea_lumped_formation":
        return {"block": "lumped_formation_Ea_kJ_mol", "key": "", "kind": "Ea"}
    if key.startswith("Ea_decay_"):
        return {"block": "decay_Ea_kJ_mol", "key": key[len("Ea_decay_"):], "kind": "Ea"}
    if key == "ph_acid_yield_per_sink_event":
        return {"block": "ph_drift", "key": "acid_yield_per_sink_event", "kind": "acid_yield"}
    if key == "ph_arp_secondary_ammonium_pKa":
        return {"block": "ph_drift", "key": "arp_secondary_ammonium_pKa", "kind": "pka"}
    raise KeyError(key)


def frozen_vector(report: Dict[str, Any]) -> np.ndarray:
    fr = report["frozen_parameters"]
    return np.array(
        [fr["log10_k_ref_at_145C"][k] for k in B23.PARAM_ORDER]
        + [fr["lumped_formation_Ea_kJ_mol"]]
        + [fr["decay_Ea_kJ_mol"][f] for f in B23.DECAY_FAMILY_ORDER]
        + [fr["ph_drift"]["acid_yield_per_sink_event"], fr["ph_drift"]["arp_secondary_ammonium_pKa"]],
        dtype=float,
    )


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--careful", action="store_true", help="careful pH fixed point (slower) for the Jacobian")
    parser.add_argument("--output", default=str(OUT))
    args = parser.parse_args(argv)
    quick = not args.careful

    report = json.loads(Path(B8.OUT_FIT_REPORT).read_text(encoding="utf-8"))
    x_opt = frozen_vector(report)
    lower, upper = B8.full_bounds()
    idx = np.array(B8.FREE_INDEX, dtype=int)
    keys = list(B8.FREE_KEYS)
    assert len(keys) == len(idx) == 23
    coords = [coordinate_of(k) for k in keys]

    t0 = time.time()
    r0 = B8.residual_vector(x_opt, quick)
    n_rows, n_free = len(r0), len(idx)
    print(f"[laplace] rows={n_rows} free={n_free} cost_at_optimum={0.5 * float(r0 @ r0):.4f} ({time.time() - t0:.1f}s/eval, quick={quick})", flush=True)

    jac = np.zeros((n_rows, n_free))
    steps: List[float] = []
    sides: List[str] = []
    for j, (i, coord) in enumerate(zip(idx, coords)):
        h = STEPS[coord["kind"]]
        lo, hi = float(lower[i]), float(upper[i])
        x0 = float(x_opt[i])
        up_ok, down_ok = x0 + h <= hi, x0 - h >= lo
        if up_ok and down_ok:
            xp, xm = x_opt.copy(), x_opt.copy(); xp[i] += h; xm[i] -= h
            jac[:, j] = (B8.residual_vector(xp, quick) - B8.residual_vector(xm, quick)) / (2 * h); sides.append("central")
        elif up_ok:
            xp = x_opt.copy(); xp[i] += h
            jac[:, j] = (B8.residual_vector(xp, quick) - r0) / h; sides.append("forward (at lower bound)")
        else:
            xm = x_opt.copy(); xm[i] -= h
            jac[:, j] = (r0 - B8.residual_vector(xm, quick)) / h; sides.append("backward (at upper bound)")
        steps.append(h)
        print(f"[laplace] {j + 1:2d}/{n_free} {keys[j]:34s} |dr/dx| = {np.linalg.norm(jac[:, j]):.3g} ({sides[-1]})", flush=True)

    dof = max(n_rows - n_free, 1)
    chi2_red = float(r0 @ r0) / dof
    jtj = jac.T @ jac
    cov = np.linalg.pinv(jtj) * chi2_red
    sigmas = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    eig = np.linalg.eigvalsh(jtj)
    # A coordinate is IDENTIFIED when its sigma is below the kind's threshold AND the
    # Gaussian is a usable local picture: an optimum sitting ON a bound with a sigma that
    # spans more than a tenth of the declared band is a one-sided, poorly constrained
    # direction (the acid yield at 4e-6 with sigma 0.29 over [0, 1], the thiol-sink Ea at
    # its 102 kJ/mol ceiling with sigma 27) -- sampling a Gaussian there would invent mass
    # the fit never saw. Such coordinates stay at the optimum and are flagged `at_bound`.
    at_bound = []
    identified = []
    for j, (i, s, c) in enumerate(zip(idx, sigmas, coords)):
        lo, hi = float(lower[i]), float(upper[i])
        width = hi - lo
        x0 = float(x_opt[i])
        on_bound = (x0 - lo) <= 1e-3 * width or (hi - x0) <= 1e-3 * width
        at_bound.append(bool(on_bound))
        ok = bool(np.isfinite(s) and s <= THRESHOLDS[c["kind"]])
        if ok and on_bound and s > 0.1 * width:
            ok = False
        identified.append(ok)
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = cov / np.outer(sigmas, sigmas)
        corr[~np.isfinite(corr)] = 0.0

    payload = {
        "artifact": "kinetic_core_b8_laplace_covariance",
        "generated_on": time.strftime("%Y-%m-%d"),
        "fit_report": data_paths.rel(Path(B8.OUT_FIT_REPORT)),
        "method": (
            "Gauss-Newton / Laplace at the frozen B8 optimum: Sigma = pinv(J^T J) * chi2_red, "
            "J by finite differences of B8's sigma-weighted residual vector over its 62 rows; "
            f"{'quick' if quick else 'careful'} pH fixed point"
        ),
        "n_rows": int(n_rows),
        "n_free": int(n_free),
        "dof": int(dof),
        "cost_at_optimum": 0.5 * float(r0 @ r0),
        "reduced_chi_square": chi2_red,
        "keys": keys,
        "coordinates": coords,
        "optimum": [float(x_opt[i]) for i in idx],
        "bounds": [[float(lower[i]), float(upper[i])] for i in idx],
        "steps": steps,
        "difference_sides": sides,
        "sigma": [float(s) for s in sigmas],
        "ci95_halfwidth": [None if not ok else 1.96 * float(s) for s, ok in zip(sigmas, identified)],
        "identified": identified,
        "at_bound": at_bound,
        "identified_count": int(sum(identified)),
        "identification_rule": (
            "sigma <= threshold for the coordinate's kind, AND NOT (optimum on a declared bound "
            "with sigma > 10 % of the band width)"
        ),
        "identification_thresholds": THRESHOLDS,
        "covariance": cov.tolist(),
        "correlation": corr.tolist(),
        "jtj_eigenvalues": [float(v) for v in eig],
        "jacobian_column_norms": [float(np.linalg.norm(jac[:, j])) for j in range(n_free)],
        "residuals_at_optimum": [float(v) for v in r0],
        "row_ids": [row["id"] for row in B23.ACTIVE_FIT_ROWS],
        "wall_seconds": round(time.time() - t0, 1),
        "reading": (
            "Only IDENTIFIED directions carry information; the envelope samples those from Sigma "
            "and clips to the declared bounds. An unidentified coordinate (sigma above its "
            "threshold, or a flat J column) stays at the optimum and is reported, not widened."
        ),
    }
    Path(args.output).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.output}: identified {payload['identified_count']}/{n_free}, chi2_red={chi2_red:.3f}, wall={payload['wall_seconds']}s")
    for k, s, ok, c, ab in zip(keys, sigmas, identified, coords, at_bound):
        print(f"  {k:34s} sigma={s:9.3g} {'IDENTIFIED' if ok else 'not sampled'} ({c['kind']}{', at bound' if ab else ''})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
