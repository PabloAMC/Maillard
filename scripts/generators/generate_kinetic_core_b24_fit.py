#!/usr/bin/env python
"""
Build Wave B24 -- 2-ACETYL-1-PYRROLINE FROM PROLINE (2026-09-09).

Pre-registered in ``results/validation/kinetic_core_b24_prereg.md`` before any B24 number existed.
Two steps entered the trunk network (network.PROLINE_REACTIONS). WHAT IS FITTED: log10 k_pyrl_ap
and log10 k_mgo_pro at 100 C, on Hofmann & Schieberle 1998b's 30-minute yields at 100 C, pH 7,
0.5 M phosphate: Table 7 experiments 1 and 2 (fed 1-pyrroline 2 mmol/L + methylglyoxal 10 / 2) and
Table 9 (proline 400 mmol/L + methylglyoxal 4 / 40 / 400), each pot integrated with the fed species
charged directly and the yield read at 30 min. Barriers declared. Two starts, least_squares (trf),
Laplace at the optimum. Diagnostics: Hofmann's experiment 3 (excess pyrroline), the apparent barrier
of the whole cascade in a glucose + proline pot against Chan & Reineccius 1994's 60.2 kJ/mol, and
the B18 test pot's pyrazines with and without proline.

Usage:
    python scripts/generators/generate_kinetic_core_b24_fit.py
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
from src.kinetic_core.parameters_proline import FROZEN_B24, PROLINE_COORDINATES, PROLINE_FIT_PH, with_fitted_proline  # noqa: E402

WAVE = "B24"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b24_prereg.md"
OUT_JSON = data_paths.VALIDATION_DIR / "kinetic_core_b24_fit_report.json"
OUT_MD = data_paths.VALIDATION_DIR / "kinetic_core_b24_fit_report.md"
CELSIUS = 273.15
SEED = 20260909
HOF = ("Hofmann & Schieberle 1998b, JAFC 46:2270, Tables 7 and 9 (0.5 M phosphate pH 7, 100 C, 30 min); hofmann1998b_extraction.md sec. 4")
#: id, initial (mmol/L), yield basis species, printed mol %, sigma, decisive
ROWS_RAW = [
    ("hof_t7_e1_pyrl2_mgo10", {"PYRL": 2.0, "MGO": 10.0}, "PYRL", 28.7, 0.15, True),
    ("hof_t7_e2_pyrl2_mgo2", {"PYRL": 2.0, "MGO": 2.0}, "PYRL", 5.3, 0.15, True),
    ("hof_t9_pro400_mgo4", {"PRO": 400.0, "Gly": 400.0, "MGO": 4.0}, "PRO", 0.0058, 0.25, True),
    ("hof_t9_pro400_mgo40", {"PRO": 400.0, "Gly": 400.0, "MGO": 40.0}, "PRO", 0.0125, 0.25, True),
    ("hof_t9_pro400_mgo400", {"PRO": 400.0, "Gly": 400.0, "MGO": 400.0}, "PRO", 0.0179, 0.25, True),
]
MINUTES = 30.0
KEYS: Tuple[str, ...] = PROLINE_COORDINATES
PRIOR = np.array([FROZEN_B24[k] for k in KEYS], dtype=float)
LOWER, UPPER = PRIOR - 3.0, PRIOR + 3.0
BASE = dict(operative_parameters(b1_fitted()))


def rows() -> Tuple[Dict[str, Any], ...]:
    return tuple(dict(id=r[0], initial=r[1], basis=r[2], target_mol_pct=r[3], sigma_log=r[4], decisive=r[5], anchor=HOF, kind="fed_yield_mol_pct")
                 for r in ROWS_RAW)


ROWS = rows()


def parameters_for(x):
    p = dict(BASE)
    p.update(with_fitted_proline(float(x[0]), float(x[1])))
    return p


def _run(x, initial, t_c, minutes, ph, grid=None):
    params, _ = trunk_conditions.apply(parameters_for(x), SimpleNamespace(ph=float(ph), water_activity=None))
    grid = np.array([0.0, minutes]) if grid is None else np.asarray(grid, dtype=float)
    return integrate(params, t_c + CELSIUS, initial, grid, rtol=1e-8, atol=1e-16)


def predictions(x) -> Dict[str, float]:
    out = {}
    for r in ROWS:
        run = _run(x, r["initial"], 100.0, MINUTES, PROLINE_FIT_PH)
        out[r["id"]] = 100.0 * float(run.series("AP")[-1]) / float(r["initial"][r["basis"]])
    return out


def residuals(x) -> np.ndarray:
    pred = predictions(x)
    return np.array([max(-1e3, min(1e3, math.log10((pred[r["id"]] + 1e-24) / r["target_mol_pct"]) / r["sigma_log"])) for r in ROWS])


def start_vector(start: int) -> np.ndarray:
    if start == 0:
        return PRIOR.copy()
    rng = np.random.default_rng(SEED + start)
    return np.clip(PRIOR + rng.normal(0.0, 0.5, size=PRIOR.shape), LOWER, UPPER)


def fit_member(start: int, max_nfev: int) -> Dict[str, Any]:
    from scipy.optimize import least_squares

    x0 = start_vector(start)
    sol = least_squares(residuals, x0, bounds=(LOWER, UPPER), method="trf", max_nfev=max_nfev, xtol=1e-10, ftol=1e-10)
    r = residuals(sol.x)
    return {"start": start, "x0": x0.tolist(), "x": sol.x.tolist(), "cost": float(np.sum(r * r)), "nfev": int(sol.nfev), "status": int(sol.status),
            "residuals_dex": {row["id"]: float(v) * float(row["sigma_log"]) for row, v in zip(ROWS, r)}}


def laplace(x, r) -> Dict[str, Any]:
    h = 0.02
    jac = np.zeros((len(ROWS), len(KEYS)))
    for j in range(len(KEYS)):
        xp, xm = x.copy(), x.copy()
        xp[j] += h
        xm[j] -= h
        jac[:, j] = (residuals(xp) - residuals(xm)) / (2 * h)
    dof = max(len(ROWS) - len(KEYS), 1)
    chi2_red = float(r @ r) / dof
    cov = np.linalg.pinv(jac.T @ jac) * max(chi2_red, 1.0)
    sig = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    on_bound = [bool(abs(x[i] - LOWER[i]) < 1e-3 or abs(UPPER[i] - x[i]) < 1e-3) for i in range(len(KEYS))]
    return {"sigma": dict(zip(KEYS, sig.tolist())), "chi2_reduced": chi2_red, "dof": dof,
            "identified": {k: bool(s < 1.0 and not b) for k, s, b in zip(KEYS, sig, on_bound)}, "on_bound": dict(zip(KEYS, on_bound))}


def diagnostics(x) -> Dict[str, Any]:
    e3 = _run(x, {"PYRL": 10.0, "MGO": 2.0}, 100.0, MINUTES, PROLINE_FIT_PH)
    per_mgo = 100.0 * float(e3.series("AP")[-1]) / 2.0
    temps = (75.0, 95.0, 115.0)
    rates = {}
    for t in temps:
        run = _run(x, {"Glc": 100.0, "PRO": 100.0, "Gly": 100.0}, t, 60.0, 7.0, grid=[0.0, 30.0, 60.0])
        rates[t] = max(float(run.series("AP")[-1]) / 60.0, 1e-30)
    xs = np.array([1.0 / (t + CELSIUS) for t in temps])
    ys = np.array([math.log(rates[t]) for t in temps])
    slope = float(np.polyfit(xs, ys, 1)[0])
    ea_apparent = -slope * 8.314e-3
    b18 = {}
    for tag, init in (("without_proline", {"Glc": 100.0, "Gly": 100.0}), ("with_proline_10mM", {"Glc": 100.0, "Gly": 110.0, "PRO": 10.0})):
        run = _run(x, init, 120.0, 60.0, 6.8, grid=[0.0, 60.0])
        b18[tag] = {k: float(run.series(k)[-1]) for k in ("PZ", "MPZ", "DMP", "AP")}
    return {"hofmann_expt3_pyrl10_mgo2": {"model_mol_pct_of_mgo": per_mgo, "printed_mol_pct_of_mgo": 0.33, "dex": math.log10(per_mgo / 0.33) if per_mgo > 0 else float("-inf")},
            "apparent_ea_glucose_proline_pot": {"model_kj_mol": ea_apparent, "chan1994b_kj_mol": 60.2, "temps_c": list(temps), "rates_mmol_l_min": rates},
            "b18_pot_120C_60min_mmol_l": b18}


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
    return {"artifact": "kinetic_core_b24_fit_report", "wave": f"{WAVE} -- 2-acetyl-1-pyrroline from proline",
            "generated_on": date.today().isoformat(), "generated_by": "scripts/generators/generate_kinetic_core_b24_fit.py", "git": _git_head(),
            "prereg": data_paths.rel(PREREG), "declaration": "docs/reference/FIT_HOLDOUT_DECLARATION.md Amendment 33",
            "objective": {"form": "five 30-minute yields (Hofmann 1998b Tables 7 and 9, 100 C, pH 7), log10 residuals; two log10 constants at 100 C free, barriers declared",
                          "n_rows": len(ROWS), "n_free_parameters": len(KEYS), "final_cost": best["cost"], "best_start": best["start"], "reduced_chi2": lap["chi2_reduced"]},
            "rows": [dict({k: v for k, v in rr.items()}, predicted_mol_pct=pred[rr["id"]], residual_dex=best["residuals_dex"][rr["id"]]) for rr in ROWS],
            "members": members, "frozen_parameters": {"proline": dict(zip(KEYS, x.tolist()))},
            "bands_log10": {"lower": dict(zip(KEYS, LOWER.tolist())), "upper": dict(zip(KEYS, UPPER.tolist()))},
            "laplace": lap, "diagnostics": diagnostics(x)}


def render(p) -> str:
    L = [f"# {p['wave']}", "", f"*Generated {p['generated_on']} by `{p['generated_by']}`; pre-registration `{p['prereg']}`; {p['declaration']}.*", "",
         f"Objective: {p['objective']['form']}. Cost {p['objective']['final_cost']:.3f} on {p['objective']['n_rows']} rows, reduced chi-square {p['objective']['reduced_chi2']:.2f}.", "",
         "## The fitted constants (log10 at 100 C, L/(mmol min); barriers declared)", "", "| coordinate | value | sigma (dex) | identified |", "|---|---|---|---|"]
    for k in KEYS:
        L.append(f"| {k} | {p['frozen_parameters']['proline'][k]:.4f} | {p['laplace']['sigma'][k]:.3f} | {p['laplace']['identified'][k]} |")
    L += ["", "## The rows (yield at 30 min, mol % of the basis species)", "", "| row | printed | model | residual (dex) |", "|---|---|---|---|"]
    for r in p["rows"]:
        L.append(f"| {r['id']} | {r['target_mol_pct']:.4g} | {r['predicted_mol_pct']:.4g} | {r['residual_dex']:+.2f} |")
    d = p["diagnostics"]
    e3, ea, b = d["hofmann_expt3_pyrl10_mgo2"], d["apparent_ea_glucose_proline_pot"], d["b18_pot_120C_60min_mmol_l"]
    L += ["", "## Diagnostics", "",
          f"- Hofmann experiment 3 (1-pyrroline 10 + methylglyoxal 2 mmol/L): model {e3['model_mol_pct_of_mgo']:.3g} mol % of the methylglyoxal vs printed 0.33 ({e3['dex']:+.2f} dex): the source's suppression by excess pyrroline is not written",
          f"- apparent barrier of the whole cascade, glucose 100 + proline 100 mmol/L, pH 7, 75-115 C: model {ea['model_kj_mol']:.0f} kJ/mol vs Chan & Reineccius 1994's 60.2",
          f"- the B18 pot at 120 C / 60 min, mmol/L: without proline {json.dumps({k: float(f'{v:.3g}') for k, v in b['without_proline'].items()})}; with 10 mmol/L proline {json.dumps({k: float(f'{v:.3g}') for k, v in b['with_proline_10mM'].items()})}"]
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-nfev", type=int, default=150)
    args = parser.parse_args(argv)
    assert PREREG.exists(), "B24 is pre-registered; write the prereg before running the fit"
    payload = build(args.max_nfev)
    OUT_JSON.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    OUT_MD.write_text(render(payload), encoding="utf-8")
    print(f"wrote {OUT_JSON}: cost {payload['objective']['final_cost']:.3f}; frozen {json.dumps(payload['frozen_parameters']['proline'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
