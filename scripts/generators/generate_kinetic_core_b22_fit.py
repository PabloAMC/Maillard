#!/usr/bin/env python
"""
Build Wave B22 -- THE METHIONINE CHAIN ON THE SUGAR PATH (2026-09-09).

Pre-registered in ``results/validation/kinetic_core_b22_prereg.md`` before any B22 number existed.
Four steps entered the trunk network (network.METHIONINE_REACTIONS). WHAT IS FITTED: four
coordinates (the methionine-to-glycine identity ratio on B18's two Strecker constants; log10
k_mtal_msh and its barrier; log10 k_msh_dmds with a declared barrier) on Pan 2025 Table 2's nine
zero-order constants, each modelled as the mean formation rate over 30-600 s in Pan's pot
(methionine 0.268 + fructose 111 + glucose 83 mmol/L, pH 6.2; sucrose omitted, declared), the unit
read as micromoles per litre per second (inferred, declared). Two starts, scipy least_squares (trf)
on log10 residuals, Laplace at the optimum. Diagnostics: Deng 2022's methional levels at 120 C
(methionine 200 + glucose 200 mmol/L, pH 7.5) and Chin & Lindsay 1994's methanethiol half-life.

Usage:
    python scripts/generators/generate_kinetic_core_b22_fit.py
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from datetime import date
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402
from src.kinetic_core import operative_parameters, trunk_conditions  # noqa: E402
from src.kinetic_core.engine import b1_fitted  # noqa: E402
from src.kinetic_core.integrate import integrate  # noqa: E402
from src.kinetic_core.parameters_methionine import (  # noqa: E402
    EA_MTAL_MSH_BAND_KJ_MOL, FROZEN_B22, METHIONINE_COORDINATES, METHIONINE_FIT_PH, with_fitted_methionine,
)

WAVE = "B22"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b22_prereg.md"
OUT_JSON = data_paths.VALIDATION_DIR / "kinetic_core_b22_fit_report.json"
OUT_MD = data_paths.VALIDATION_DIR / "kinetic_core_b22_fit_report.md"
CELSIUS = 273.15
SEED = 20260909
PAN_INITIAL = {"Fru": 111.0, "Glc": 83.0, "MET": 0.268, "Gly": 0.268}     # glycine = methionine's amine, declared
PAN_WINDOW_S = (30.0, 600.0)
PAN_ANCHOR = ("Pan et al. 2025 Table 2 (unit inferred as umol L-1 s-1 from the printed 626.31 ug/L methional and 26.43 ug/L "
              "methanethiol at 140 C / 600 s; pan2025_extraction.md sec. 4)")
#: Pan's zero-order constants, converted to umol L-1 min-1 (x 60), per species and temperature.
PAN_RATES_UMOL_L_MIN: Dict[str, Dict[float, float]] = {
    "MTAL": {100.0: 1.823e-4 * 60, 120.0: 16.8e-4 * 60, 140.0: 89.9e-4 * 60},
    "MSH": {100.0: 1.579e-4 * 60, 120.0: 2.342e-4 * 60, 140.0: 6.335e-4 * 60},
    "DMDS": {100.0: 0.015e-4 * 60, 120.0: 0.070e-4 * 60, 140.0: 0.187e-4 * 60},
}
SIGMA = {"MTAL": 0.15, "MSH": 0.15, "DMDS": 0.25}
DENG_PRINTED_UG_L = {30.0: 18.52, 60.0: 25.27, 120.0: 75.76, 180.0: 87.77}
MW_MTAL = 104.17


def rows() -> Tuple[Dict[str, Any], ...]:
    out: List[Dict[str, Any]] = []
    for sp, by_t in PAN_RATES_UMOL_L_MIN.items():
        for t_c, k in by_t.items():
            out.append(dict(id=f"pan_{sp}_rate_{int(t_c)}C", species=sp, t_c=t_c, target=k, sigma_log=SIGMA[sp], anchor=PAN_ANCHOR,
                            kind="mean_rate_umol_l_min", decisive=(sp != "DMDS")))
    return tuple(out)


ROWS = rows()
KEYS: Tuple[str, ...] = METHIONINE_COORDINATES
PRIOR = np.array([FROZEN_B22[k] for k in KEYS], dtype=float)
LOWER = np.array([-2.0, PRIOR[1] - 3.0, EA_MTAL_MSH_BAND_KJ_MOL[0], PRIOR[3] - 3.0])
UPPER = np.array([2.0, PRIOR[1] + 3.0, EA_MTAL_MSH_BAND_KJ_MOL[1], PRIOR[3] + 3.0])
BASE = dict(operative_parameters(b1_fitted()))


def parameters_for(x: np.ndarray) -> Dict[str, Any]:
    p = dict(BASE)
    p.update(with_fitted_methionine(*[float(v) for v in x]))
    return p


def _run(x, initial, t_c, minutes, ph, grid=None):
    params, _ = trunk_conditions.apply(parameters_for(x), SimpleNamespace(ph=float(ph), water_activity=None))
    grid = np.array([0.0, PAN_WINDOW_S[0] / 60.0, PAN_WINDOW_S[1] / 60.0]) if grid is None else np.asarray(grid, dtype=float)
    return integrate(params, t_c + CELSIUS, initial, grid, rtol=1e-8, atol=1e-16)


def predictions(x: np.ndarray) -> Dict[str, float]:
    pred: Dict[str, float] = {}
    runs = {t: _run(x, PAN_INITIAL, t, PAN_WINDOW_S[1] / 60.0, METHIONINE_FIT_PH) for t in (100.0, 120.0, 140.0)}
    span_min = (PAN_WINDOW_S[1] - PAN_WINDOW_S[0]) / 60.0
    for r in ROWS:
        c = runs[r["t_c"]].series(r["species"])
        pred[r["id"]] = float(c[-1] - c[-2]) / span_min * 1000.0
    return pred


def residuals(x: np.ndarray) -> np.ndarray:
    pred = predictions(x)
    out = np.empty(len(ROWS))
    for i, r in enumerate(ROWS):
        p, t = pred[r["id"]], float(r["target"])
        # no clipping short of the numerical floor: a residual of forty sigma must still carry its gradient,
        # or an optimiser started four decades away never moves (the first run of this generator did not)
        v = math.log10((p + 1e-24) / t) / float(r["sigma_log"]) if (math.isfinite(p) and p > 0) else 1.0e3
        out[i] = max(-1.0e3, min(1.0e3, v))
    return out


def start_vector(start: int) -> np.ndarray:
    if start == 0:
        return PRIOR.copy()
    rng = np.random.default_rng(SEED + start)
    return np.clip(PRIOR + rng.normal(0.0, 0.5, size=PRIOR.shape) * np.array([1.0, 1.0, 30.0, 1.0]), LOWER, UPPER)


def fit_member(start: int, max_nfev: int) -> Dict[str, Any]:
    from scipy.optimize import least_squares

    x0 = start_vector(start)
    sol = least_squares(residuals, x0, bounds=(LOWER, UPPER), method="trf", max_nfev=max_nfev, xtol=1e-10, ftol=1e-10,
                        x_scale=np.array([1.0, 1.0, 30.0, 1.0]))
    r = residuals(sol.x)
    return {"start": start, "x0": x0.tolist(), "x": sol.x.tolist(), "cost": float(np.sum(r * r)), "nfev": int(sol.nfev), "status": int(sol.status),
            "residuals_dex": {row["id"]: float(v) * float(row["sigma_log"]) for row, v in zip(ROWS, r)}}


def laplace(x: np.ndarray, r: np.ndarray) -> Dict[str, Any]:
    steps = np.array([0.02, 0.02, 1.0, 0.02])
    jac = np.zeros((len(ROWS), len(KEYS)))
    for j in range(len(KEYS)):
        xp, xm = x.copy(), x.copy()
        xp[j] = min(xp[j] + steps[j], UPPER[j])
        xm[j] = max(xm[j] - steps[j], LOWER[j])
        jac[:, j] = (residuals(xp) - residuals(xm)) / (xp[j] - xm[j])
    dof = max(len(ROWS) - len(KEYS), 1)
    chi2_red = float(r @ r) / dof
    cov = np.linalg.pinv(jac.T @ jac) * max(chi2_red, 1.0)
    sig = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    on_bound = [bool(abs(x[i] - LOWER[i]) < 1e-3 * (UPPER[i] - LOWER[i]) or abs(UPPER[i] - x[i]) < 1e-3 * (UPPER[i] - LOWER[i])) for i in range(len(KEYS))]
    thresh = [1.0, 1.0, 60.0, 1.0]
    return {"sigma": dict(zip(KEYS, sig.tolist())), "chi2_reduced": chi2_red, "dof": dof,
            "identified": {k: bool(s < t and not b) for k, s, t, b in zip(KEYS, sig, thresh, on_bound)}, "on_bound": dict(zip(KEYS, on_bound))}


def diagnostics(x: np.ndarray) -> Dict[str, Any]:
    grid = [0.0, 30.0, 60.0, 120.0, 180.0]
    run = _run(x, {"MET": 200.0, "Glc": 200.0, "Gly": 200.0}, 120.0, 180.0, 7.5, grid=grid)
    mtal = run.series("MTAL")
    deng = {f"{int(t)}min": {"model_ug_per_l": float(mtal[i]) * MW_MTAL * 1000.0, "printed_ug_per_l": DENG_PRINTED_UG_L[t],
                              "dex": math.log10(float(mtal[i]) * MW_MTAL * 1000.0 / DENG_PRINTED_UG_L[t]) if mtal[i] > 0 else float("-inf")}
            for i, t in enumerate(grid) if t in DENG_PRINTED_UG_L}
    deng["rising_30_to_120"] = bool(mtal[3] > mtal[1])
    p = parameters_for(x)
    k30 = p["k_msh_dmds"].k_at(303.15)      # L/(mmol min)
    half_life_min = math.log(2.0) / (k30 * 0.0416) if k30 > 0 else float("inf")   # second order at 41.6 uM, approximate
    pan140 = _run(x, PAN_INITIAL, 140.0, 10.0, METHIONINE_FIT_PH)
    supply = {k: float(pan140.series(k)[-1]) for k in ("GO", "MGO", "MTAL", "MSH", "DMDS", "MET")}
    return {"deng2022_methional_120C": deng,
            "chin1994_methanethiol_half_life_30C": {"model_min": half_life_min, "chin_min_with_cu": 17.0, "note": "Chin's is copper 1 ppm in air; the model's constant is the apparent one from Pan's pot"},
            "pan_pot_140C_10min_mmol_l": supply}


def _git_head() -> Dict[str, str]:
    try:
        return {"head": subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(ROOT), capture_output=True, text=True, timeout=10).stdout.strip()}
    except Exception:  # noqa: BLE001
        return {"head": "unknown"}


def build(max_nfev: int) -> Dict[str, Any]:
    members = [fit_member(s, max_nfev) for s in (0, 1)]
    best = min(members, key=lambda m: m["cost"])
    x = np.array(best["x"], dtype=float)
    r = residuals(x)
    lap = laplace(x, r)
    pred = predictions(x)
    return {"artifact": "kinetic_core_b22_fit_report", "wave": f"{WAVE} -- the methionine chain on the sugar path",
            "generated_on": date.today().isoformat(), "generated_by": "scripts/generators/generate_kinetic_core_b22_fit.py", "git": _git_head(),
            "prereg": data_paths.rel(PREREG), "declaration": "docs/reference/FIT_HOLDOUT_DECLARATION.md Amendment 32",
            "objective": {"form": "nine zero-order rates (Pan 2025, 100 / 120 / 140 C) as mean formation rates over 30-600 s in Pan's pot; four free",
                          "n_rows": len(ROWS), "n_free_parameters": len(KEYS), "final_cost": best["cost"], "best_start": best["start"], "reduced_chi2": lap["chi2_reduced"]},
            "rows": [dict(rr, predicted=pred[rr["id"]], residual_dex=best["residuals_dex"][rr["id"]]) for rr in ROWS],
            "members": members, "frozen_parameters": {"methionine": dict(zip(KEYS, x.tolist()))},
            "bands": {"lower": dict(zip(KEYS, LOWER.tolist())), "upper": dict(zip(KEYS, UPPER.tolist()))},
            "laplace": lap, "diagnostics": diagnostics(x), "declared": {"pan_initial_mmol_l": PAN_INITIAL, "unit": "umol L-1 s-1 inferred", "sucrose": "omitted"}}


def render(p: Dict[str, Any]) -> str:
    L = [f"# {p['wave']}", "", f"*Generated {p['generated_on']} by `{p['generated_by']}`; pre-registration `{p['prereg']}`; {p['declaration']}.*", "",
         f"Objective: {p['objective']['form']}. Cost {p['objective']['final_cost']:.2f} on {p['objective']['n_rows']} rows, reduced chi-square {p['objective']['reduced_chi2']:.1f}.", "",
         "## The fitted coordinates", "", "| coordinate | value | sigma | identified | on bound | band |", "|---|---|---|---|---|---|"]
    for k in KEYS:
        L.append(f"| {k} | {p['frozen_parameters']['methionine'][k]:.3f} | {p['laplace']['sigma'][k]:.3g} | {p['laplace']['identified'][k]} | {p['laplace']['on_bound'][k]} | [{p['bands']['lower'][k]:.1f}, {p['bands']['upper'][k]:.1f}] |")
    L += ["", "## The rows (mean formation rate, umol L-1 min-1)", "", "| row | printed | model | residual (dex) | decisive |", "|---|---|---|---|---|"]
    for r in p["rows"]:
        L.append(f"| {r['id']} | {r['target']:.3g} | {r['predicted']:.3g} | {r['residual_dex']:+.2f} | {r['decisive']} |")
    d = p["diagnostics"]
    L += ["", "## Diagnostics", ""]
    for t, v in d["deng2022_methional_120C"].items():
        if isinstance(v, dict):
            L.append(f"- Deng 2022, methional at 120 C, {t}: model {v['model_ug_per_l']:.3g} ug/L vs printed {v['printed_ug_per_l']} ({v['dex']:+.2f} dex)")
    L.append(f"- Deng 2022, rising 30 -> 120 min: {d['deng2022_methional_120C']['rising_30_to_120']}")
    c = d["chin1994_methanethiol_half_life_30C"]
    L.append(f"- Chin & Lindsay 1994, methanethiol half-life at 30 C: model {c['model_min']:.3g} min vs {c['chin_min_with_cu']} min with copper")
    L.append(f"- Pan's pot at 140 C / 10 min, mmol/L: {json.dumps({k: float(f'{v:.3g}') for k, v in d['pan_pot_140C_10min_mmol_l'].items()})}")
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-nfev", type=int, default=150)
    args = parser.parse_args(argv)
    assert PREREG.exists(), "B22 is pre-registered; write the prereg before running the fit"
    payload = build(args.max_nfev)
    OUT_JSON.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    OUT_MD.write_text(render(payload), encoding="utf-8")
    print(f"wrote {OUT_JSON}: cost {payload['objective']['final_cost']:.2f}; frozen {json.dumps(payload['frozen_parameters']['methionine'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
