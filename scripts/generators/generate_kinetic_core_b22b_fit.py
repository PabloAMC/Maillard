#!/usr/bin/env python
"""
Wave B22b -- THE ROUTE DENG'S OWN EXPERIMENT NAMES (2026-09-10).

`results/validation/kinetic_core_b22b_prereg.md`. Two free coordinates on five rows: Deng 2022
Table 1's methional time course from the FED methionine-glucose Amadori compound at 120 C, which
rises to 120 minutes and then falls.

B22 built methional as free dicarbonyl times methionine and two laboratories refuted it. Deng's own
comparison names the successor: the fed Amadori compound gives 1.4 to 2.6 times MORE methional than
methionine plus glucose, so the dicarbonyl arrives inside the molecule.

The binary Met + Glc arm of the same table is a CHECK and not a fit row: it depends on the trunk's
Amadori route charged as glycine by declaration, which is the substitution B22 already found
overshoots. If it is wrong the check says so.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths, provenance  # noqa: E402
from src.kinetic_core import operative_parameters, trunk_conditions  # noqa: E402
from src.kinetic_core.engine import b1_fitted  # noqa: E402
from src.kinetic_core.integrate import integrate  # noqa: E402
from src.kinetic_core.parameters_methionine import (  # noqa: E402
    DENG_ARP_METHIONAL_UMOL_L, DENG_ARP_MMOL_L, DENG_BINARY_METHIONAL_UMOL_L, DENG_PH, DENG_SOURCE,
    DENG_T_C, METHIONINE_B22B_COORDINATES, with_fitted_methionine_b22b,
)

WAVE = "B22b"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b22b_prereg.md"
OUT_JSON = data_paths.VALIDATION_DIR / "kinetic_core_b22b_fit_report.json"
OUT_MD = data_paths.VALIDATION_DIR / "kinetic_core_b22b_fit_report.md"
CELSIUS = 273.15
SEED = 20260910
SIGMA_LOG = 0.30
KEYS: Tuple[str, ...] = METHIONINE_B22B_COORDINATES
#: The Amadori compound's decomposition starts near "half of it in an hour"; its competing loss
#: three times slower, which is what a rise-then-fall with a peak at 120 min roughly implies.
PRIOR = np.array([math.log10(0.012), math.log10(0.004)], dtype=float)
LOWER, UPPER = PRIOR - 4.0, PRIOR + 4.0
BASE = dict(operative_parameters(b1_fitted()))
TIMES = tuple(sorted(DENG_ARP_METHIONAL_UMOL_L))


def parameters_for(x):
    p = dict(BASE)
    p.update(with_fitted_methionine_b22b(float(x[0]), float(x[1])))
    return p


def _series(x, initial) -> Dict[float, float]:
    params, _ = trunk_conditions.apply(parameters_for(x), SimpleNamespace(ph=DENG_PH, water_activity=None))
    grid = np.array((0.0,) + TIMES, dtype=float)
    run = integrate(params, DENG_T_C + CELSIUS, initial, grid, rtol=1e-8, atol=1e-16)
    mtal = run.series("MTAL")
    # mmol/L -> umol/L
    return {t: 1.0e3 * float(mtal[i + 1]) for i, t in enumerate(TIMES)}


def predictions(x) -> Dict[float, float]:
    return _series(x, {"MARP": DENG_ARP_MMOL_L})


def residuals(x) -> np.ndarray:
    pred = predictions(x)
    return np.array([
        max(-1e3, min(1e3, math.log10((pred[t] + 1e-30) / DENG_ARP_METHIONAL_UMOL_L[t]) / SIGMA_LOG))
        for t in TIMES
    ])


def fit_member(start: int, max_nfev: int) -> Dict[str, Any]:
    from scipy.optimize import least_squares

    if start == 0:
        x0 = PRIOR.copy()
    else:
        rng = np.random.default_rng(SEED + start)
        x0 = np.clip(PRIOR + rng.normal(0.0, 0.8, size=PRIOR.shape), LOWER, UPPER)
    sol = least_squares(residuals, x0, bounds=(LOWER, UPPER), method="trf", max_nfev=max_nfev,
                        xtol=1e-10, ftol=1e-10)
    r = residuals(sol.x)
    return {"start": start, "x0": x0.tolist(), "x": sol.x.tolist(), "cost": float(np.sum(r * r)),
            "nfev": int(sol.nfev),
            "residuals_dex": {str(t): float(v) * SIGMA_LOG for t, v in zip(TIMES, r)}}


def laplace(x, r) -> Dict[str, Any]:
    h = 0.02
    jac = np.zeros((len(TIMES), len(KEYS)))
    for j in range(len(KEYS)):
        xp, xm = x.copy(), x.copy()
        xp[j] += h
        xm[j] -= h
        jac[:, j] = (residuals(xp) - residuals(xm)) / (2 * h)
    dof = max(len(TIMES) - len(KEYS), 1)
    chi2_red = float(r @ r) / dof
    cov = np.linalg.pinv(jac.T @ jac) * max(chi2_red, 1.0)
    sig = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    on_bound = [bool(abs(x[i] - LOWER[i]) < 1e-3 or abs(UPPER[i] - x[i]) < 1e-3) for i in range(len(KEYS))]
    return {"sigma": dict(zip(KEYS, sig.tolist())), "chi2_reduced": chi2_red, "dof": dof,
            "identified": {k: bool(s < 1.0 and not b) for k, s, b in zip(KEYS, sig, on_bound)},
            "on_bound": dict(zip(KEYS, on_bound))}


def two_arm_check(x) -> Dict[str, Any]:
    """
    The Met + Glc arm, charged the DECLARED way and never fitted. Deng measures the fed Amadori arm
    ABOVE the binary arm at every time; this asks whether the model reproduces that direction.
    """
    binary = _series(x, {"MET": 200.0, "Gly": 200.0, "Glc": 200.0})
    arp = predictions(x)
    rows = {}
    for t in TIMES:
        rows[str(t)] = {
            "printed_binary_umol_l": DENG_BINARY_METHIONAL_UMOL_L[t],
            "printed_arp_umol_l": DENG_ARP_METHIONAL_UMOL_L[t],
            "printed_arp_over_binary": DENG_ARP_METHIONAL_UMOL_L[t] / DENG_BINARY_METHIONAL_UMOL_L[t],
            "model_binary_umol_l": binary[t],
            "model_arp_umol_l": arp[t],
            "model_arp_over_binary": (arp[t] + 1e-30) / (binary[t] + 1e-30),
        }
    direction_ok = all(r["model_arp_over_binary"] > 1.0 for r in rows.values())
    return {"rows": rows, "direction_holds_at_every_time": bool(direction_ok),
            "note": ("The binary arm rests on charging methionine as GLYCINE for the Amadori "
                     "chemistry, a declared substitution B22 already found overshoots Deng's pot by "
                     "1.5 to 2.9 decades. This check is where that shows.")}


def peak_time(x) -> float:
    pred = predictions(x)
    return float(max(pred, key=lambda t: pred[t]))


def build(max_nfev: int) -> Dict[str, Any]:
    members = [fit_member(s, max_nfev) for s in range(3)]
    best = min(range(len(members)), key=lambda i: members[i]["cost"])
    x = np.array(members[best]["x"], dtype=float)
    r = residuals(x)
    pred = predictions(x)
    return {
        "wave": WAVE, "artifact": "kinetic_core_b22b_fit_report",
        "provenance": provenance.provenance_block(
            "kinetic_core_b22b_fit_report",
            generated_by="scripts/generators/generate_kinetic_core_b22b_fit.py", inputs=[]),
        "prereg": data_paths.rel(PREREG),
        "source": DENG_SOURCE,
        "objective": {"times_min": list(TIMES), "targets_umol_l": dict(DENG_ARP_METHIONAL_UMOL_L),
                      "sigma_log": SIGMA_LOG, "n_rows": len(TIMES), "n_free": len(KEYS),
                      "final_cost": float(np.sum(r * r))},
        "bounds": {k: [float(LOWER[i]), float(UPPER[i])] for i, k in enumerate(KEYS)},
        "members": members, "best_start": best,
        "frozen_parameters": {"methionine_b22b": dict(zip(KEYS, x.tolist()))},
        "predicted_umol_l": {str(t): v for t, v in pred.items()},
        "residual_by_row_dex": {str(t): float(v) * SIGMA_LOG for t, v in zip(TIMES, r)},
        "peak_time_min": peak_time(x),
        "laplace": laplace(x, r),
        "two_arm_check": two_arm_check(x),
        "reference_temperature_K": DENG_T_C + CELSIUS, "reference_ph": DENG_PH,
    }


def render(p: Dict[str, Any]) -> str:
    L = [f"# Wave {WAVE} fit report", "",
         f"Cost {p['objective']['final_cost']:.3f} on {p['objective']['n_rows']} rows, "
         f"{p['objective']['n_free']} free. Model peaks at {p['peak_time_min']:.0f} min; the source "
         f"peaks at 120.", "",
         "| minutes | printed umol/L | model | residual (dex) |", "|---:|---:|---:|---:|"]
    for t in p["objective"]["times_min"]:
        k = str(t)
        L.append(f"| {t:.0f} | {p['objective']['targets_umol_l'][k] if k in p['objective']['targets_umol_l'] else p['objective']['targets_umol_l'][t]:.4g} | "
                 f"{p['predicted_umol_l'][k]:.4g} | {p['residual_by_row_dex'][k]:+.3f} |")
    lap = p["laplace"]
    L += ["", f"Laplace: sigma {json.dumps({k: round(v, 3) for k, v in lap['sigma'].items()})}; "
              f"identified {lap['identified']}; on bound {lap['on_bound']}.", "",
          "## The two-arm check, not fitted", "",
          "| minutes | printed ARP / binary | model ARP / binary |", "|---:|---:|---:|"]
    for t, row in sorted(p["two_arm_check"]["rows"].items(), key=lambda kv: float(kv[0])):
        L.append(f"| {float(t):.0f} | {row['printed_arp_over_binary']:.3g} | "
                 f"{row['model_arp_over_binary']:.3g} |")
    L += ["", f"> {p['two_arm_check']['note']}", ""]
    return "\n".join(L)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-nfev", type=int, default=200)
    args = parser.parse_args(argv)
    assert PREREG.exists(), "B22b is pre-registered; write the prereg before running the fit"
    payload = build(args.max_nfev)
    OUT_JSON.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    OUT_MD.write_text(render(payload), encoding="utf-8")
    print(f"cost {payload['objective']['final_cost']:.3f}; peak {payload['peak_time_min']:.0f} min; "
          f"frozen {json.dumps({k: round(v, 3) for k, v in payload['frozen_parameters']['methionine_b22b'].items()})}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
