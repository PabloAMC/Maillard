#!/usr/bin/env python
"""
Build Wave B20 -- THE GLYCATION ARM (2026-09-09).

Pre-registered in ``results/validation/kinetic_core_b20_prereg.md`` before any B20 number existed.
Five steps entered the trunk network (network.GLYCATION_REACTIONS) on protein-bound lysine: the
sugar glycates it to the bound Amadori compound, which oxidises to CML, goes to CEL, or decays back
to the sugar path returning the lysine; CML is lost into the melanoidin pools.

WHAT IS FITTED: five coordinates, the log10 of each constant at the trunk's 100 C reference, on
the ten rate constants Nguyen 2016 printed for casein-bound lysine + glucose in water at 120 and
130 C (nguyen2016_extraction.md Table 1, system M1), each weighted by its printed 95 % interval.
The barriers are DECLARED (parameters_glycation), so a row is the fitted constant evaluated at the
row's temperature against the printed value; no integration enters the objective. Two starts, scipy
least_squares (trf) on log10 residuals, then the Laplace covariance at the optimum. Beside the fit:
Nguyen's pot integrated for 30 min (levels), and the comparators the prereg names.

Usage:
    python scripts/generators/generate_kinetic_core_b20_fit.py            # both starts, Laplace, report
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
from src.kinetic_core.parameters import T_REF_K  # noqa: E402
from src.kinetic_core.parameters_glycation import (  # noqa: E402
    EA_CML_LOSS_KJ_MOL, EA_FLP_CEL_KJ_MOL, EA_FLP_CML_KJ_MOL, EA_FLP_DECAY_KJ_MOL, EA_GLYC_KJ_MOL, FROZEN_B20,
    GLYCATION_COORDINATES, with_fitted_glycation,
)

WAVE = "B20"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b20_prereg.md"
OUT_JSON = data_paths.VALIDATION_DIR / "kinetic_core_b20_fit_report.json"
OUT_MD = data_paths.VALIDATION_DIR / "kinetic_core_b20_fit_report.md"
CELSIUS = 273.15
R_KJ = 8.314e-3
SEED = 20260909

# ===========================================================================
# 1. THE ROWS: Nguyen 2016 Table 1, system M1 (value, half-width of the 95 % HPD), per temperature
# ===========================================================================
NGUYEN_ANCHOR = ("Nguyen, van der Fels-Klerx & van Boekel 2016, Food Chem. 192:125, Table 1 system M1: sodium caseinate 30 g/L "
                 "(about 16 mmol/L lysine residues) + glucose 150 mmol/L, 0.1 M phosphate pH 6.8, 120 / 130 C, sealed with air; "
                 "nguyen2016_extraction.md sec. 4")
#: key -> {T_C: (value, half-width)}; k3 in L/(mmol min), the rest per minute
NGUYEN_M1: Dict[str, Dict[float, Tuple[float, float]]] = {
    "k_glyc": {120.0: (1.5e-4, 3.0e-5), 130.0: (1.6e-4, 2.2e-5)},
    "k_flp_cml": {120.0: (8.8e-3, 6.6e-3), 130.0: (6.0e-3, 1.6e-3)},
    "k_flp_cel": {120.0: (2.3e-3, 2.7e-3), 130.0: (2.0e-3, 3.0e-4)},
    "k_flp_decay": {120.0: (5.2e-2, 3.3e-2), 130.0: (1.5e-1, 3.4e-2)},
    "k_cml_loss": {120.0: (2.9e-1, 2.7e-1), 130.0: (7.7e-2, 3.3e-2)},
}
EA_OF: Dict[str, float] = {"k_glyc": EA_GLYC_KJ_MOL, "k_flp_cml": EA_FLP_CML_KJ_MOL, "k_flp_cel": EA_FLP_CEL_KJ_MOL,
                           "k_flp_decay": EA_FLP_DECAY_KJ_MOL, "k_cml_loss": EA_CML_LOSS_KJ_MOL}
KEY_OF_COORD: Dict[str, str] = {"log10_k_glyc_100C": "k_glyc", "log10_k_flp_cml_100C": "k_flp_cml", "log10_k_flp_cel_100C": "k_flp_cel",
                                "log10_k_flp_decay_100C": "k_flp_decay", "log10_k_cml_loss_100C": "k_cml_loss"}
SIGMA_SPANS_ZERO = 0.5


def rows() -> Tuple[Dict[str, Any], ...]:
    out: List[Dict[str, Any]] = []
    for key, by_t in NGUYEN_M1.items():
        for t_c, (value, half) in by_t.items():
            spans_zero = half >= value
            sigma = SIGMA_SPANS_ZERO if spans_zero else math.log10(1.0 + half / value)
            out.append(dict(id=f"nguyen_{key}_{int(t_c)}C", key=key, t_c=t_c, target=value, half_width=half,
                            sigma_log=max(sigma, 0.05), spans_zero=spans_zero, anchor=NGUYEN_ANCHOR,
                            kind="rate_constant", decisive=not spans_zero))
    return tuple(out)


ROWS = rows()
assert len(ROWS) == 10

# ===========================================================================
# 2. THE VECTOR
# ===========================================================================
KEYS: Tuple[str, ...] = GLYCATION_COORDINATES
PRIOR = np.array([FROZEN_B20[k] for k in KEYS], dtype=float)
LOWER = PRIOR - 2.0
UPPER = PRIOR + 2.0
BASE_PARAMETERS = dict(operative_parameters(b1_fitted()))


def k_at(log10_k100: float, ea: float, t_c: float) -> float:
    t_k = t_c + CELSIUS
    return 10.0 ** log10_k100 * math.exp(-ea / R_KJ * (1.0 / t_k - 1.0 / T_REF_K))


def predictions(x: np.ndarray) -> Dict[str, float]:
    coord = dict(zip(KEYS, (float(v) for v in x)))
    pred: Dict[str, float] = {}
    for r in ROWS:
        c = next(cc for cc, kk in KEY_OF_COORD.items() if kk == r["key"])
        pred[r["id"]] = k_at(coord[c], EA_OF[r["key"]], r["t_c"])
    return pred


def residuals(x: np.ndarray) -> np.ndarray:
    pred = predictions(x)
    out = np.empty(len(ROWS))
    for i, r in enumerate(ROWS):
        out[i] = math.log10(pred[r["id"]] / float(r["target"])) / float(r["sigma_log"])
    return out


def cost(x: np.ndarray) -> float:
    r = residuals(x)
    return float(np.sum(r * r))


def parameters_for(x: np.ndarray) -> Dict[str, Any]:
    p = dict(BASE_PARAMETERS)
    p.update(with_fitted_glycation(*[float(v) for v in x]))
    return p


# ===========================================================================
# 3. THE FIT
# ===========================================================================
def start_vector(start: int) -> np.ndarray:
    if start == 0:
        return PRIOR.copy()
    rng = np.random.default_rng(SEED + start)
    return np.clip(PRIOR + rng.normal(0.0, 0.5, size=PRIOR.shape), LOWER, UPPER)


def fit_member(start: int, max_nfev: int = 200) -> Dict[str, Any]:
    from scipy.optimize import least_squares

    x0 = start_vector(start)
    sol = least_squares(residuals, x0, bounds=(LOWER, UPPER), method="trf", max_nfev=max_nfev, xtol=1e-12, ftol=1e-12)
    r = residuals(sol.x)
    return {"start": start, "x0": x0.tolist(), "x": sol.x.tolist(), "cost": float(np.sum(r * r)), "nfev": int(sol.nfev),
            "status": int(sol.status), "residuals": {row["id"]: float(v) for row, v in zip(ROWS, r)},
            "residuals_dex": {row["id"]: float(v) * float(row["sigma_log"]) for row, v in zip(ROWS, r)}}


def laplace(x: np.ndarray, r: np.ndarray) -> Dict[str, Any]:
    h = 1e-4
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
            "identified": {k: bool(s < 1.0 and not b) for k, s, b in zip(KEYS, sig, on_bound)}, "on_bound": dict(zip(KEYS, on_bound)),
            "note": "Gauss-Newton at the optimum; sigma scaled by the reduced chi-square when it exceeds one (the rows' intervals are the sigmas)"}


# ===========================================================================
# 4. DIAGNOSTICS: Nguyen's pot integrated, and the comparators
# ===========================================================================
def nguyen_pot(x: np.ndarray, t_c: float, minutes: float = 30.0) -> Dict[str, float]:
    p = parameters_for(x)
    proc = SimpleNamespace(ph=6.8, water_activity=None)
    params, _notes = trunk_conditions.apply(p, proc)
    grid = np.linspace(0.0, minutes, 7)
    run = integrate(params, t_c + CELSIUS, {"Glc": 150.0, "Gly": 0.0, "LYSP": 16.0}, grid, rtol=1e-8, atol=1e-14)
    end = {k: float(run.series(k)[-1]) for k in ("LYSP", "FLP", "CML", "CEL", "TDG")}
    end["lysine_lost_fraction"] = 1.0 - end["LYSP"] / 16.0
    return end


def diagnostics(x: np.ndarray) -> Dict[str, Any]:
    coord = dict(zip(KEYS, (float(v) for v in x)))
    out: Dict[str, Any] = {"nguyen_pot_30min": {f"{int(t)}C": nguyen_pot(x, t) for t in (120.0, 130.0)}}
    out["nguyen_printed_cml_range_mmol_l"] = [0.025, 0.135]
    k_cml_180 = k_at(coord["log10_k_flp_cml_100C"], EA_FLP_CML_KJ_MOL, 180.0)
    out["berk2021_flp_to_cml_180C"] = {"model_per_min": k_cml_180, "berk_per_min": 5.54e-3, "dex": math.log10(k_cml_180 / 5.54e-3),
                                        "note": "sesame seed, dry, 180 C; Berk's own barrier (113) is the declared one, so this compares the 100 C anchor"}
    ham = {110.0: 1.7e-4, 120.0: 3.6e-3, 130.0: 1.4e-3, 140.0: 1.1e-3}
    out["hamzalioglu2026_laclys_to_cml"] = {f"{int(t)}C": {"model_per_min": k_at(coord["log10_k_flp_cml_100C"], EA_FLP_CML_KJ_MOL, t), "milk_per_min": v,
                                                           "dex": math.log10(k_at(coord["log10_k_flp_cml_100C"], EA_FLP_CML_KJ_MOL, t) / v)}
                                            for t, v in ham.items()}
    out["troise2015_direction"] = ("expanded soybean, 110 C, 60 min: about 25 % of the bound lysine lost (a moist solid with sucrose, "
                                   "which the engine cannot charge; direction only)")
    return out


def _git_head() -> Dict[str, str]:
    try:
        sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(ROOT), capture_output=True, text=True, timeout=10).stdout.strip()
    except Exception:  # noqa: BLE001
        sha = "unknown"
    return {"head": sha}


def build(max_nfev: int) -> Dict[str, Any]:
    members = [fit_member(s, max_nfev) for s in (0, 1)]
    best = min(members, key=lambda m: m["cost"])
    x = np.array(best["x"], dtype=float)
    r = residuals(x)
    lap = laplace(x, r)
    frozen = dict(zip(KEYS, x.tolist()))
    pred = predictions(x)
    return {
        "artifact": "kinetic_core_b20_fit_report",
        "wave": f"{WAVE} -- the glycation arm: protein-bound lysine as a reactant on the trunk lane",
        "generated_on": date.today().isoformat(),
        "generated_by": "scripts/generators/generate_kinetic_core_b20_fit.py",
        "git": _git_head(),
        "prereg": data_paths.rel(PREREG),
        "declaration": "docs/reference/FIT_HOLDOUT_DECLARATION.md Amendment 30",
        "objective": {"form": "ten printed rate constants (Nguyen 2016 M1, 120 / 130 C), log10 residuals over the printed interval; "
                              "five log10 constants at 100 C free, five barriers declared",
                      "n_rows": len(ROWS), "n_free_parameters": len(KEYS), "final_cost": best["cost"], "best_start": best["start"],
                      "reduced_chi2": lap["chi2_reduced"]},
        "rows": [dict(r, predicted=pred[r["id"]], residual_dex=best["residuals_dex"][r["id"]]) for r in ROWS],
        "members": members,
        "frozen_parameters": {"glycation": frozen,
                              "declared_barriers_kj_mol": EA_OF,
                              "constants_at_120C": {k: k_at(frozen[c], EA_OF[k], 120.0) for c, k in KEY_OF_COORD.items()},
                              "constants_at_130C": {k: k_at(frozen[c], EA_OF[k], 130.0) for c, k in KEY_OF_COORD.items()}},
        "bands_log10": {"lower": dict(zip(KEYS, LOWER.tolist())), "upper": dict(zip(KEYS, UPPER.tolist()))},
        "laplace": lap,
        "diagnostics": diagnostics(x),
        "start_vectors": {str(m["start"]): m["x0"] for m in members},
    }


def render(p: Dict[str, Any]) -> str:
    L = [f"# {p['wave']}", "", f"*Generated {p['generated_on']} by `{p['generated_by']}`; pre-registration `{p['prereg']}`; {p['declaration']}.*", "",
         f"Objective: {p['objective']['form']}. Cost {p['objective']['final_cost']:.3f} on {p['objective']['n_rows']} rows, reduced chi-square "
         f"{p['objective']['reduced_chi2']:.2f}; best start {p['objective']['best_start']}.", "",
         "## The fitted constants (log10 at 100 C; barriers declared)", "", "| coordinate | value | sigma (dex) | identified | band |", "|---|---|---|---|---|"]
    for k in KEYS:
        L.append(f"| {k} | {p['frozen_parameters']['glycation'][k]:.4f} | {p['laplace']['sigma'][k]:.3f} | {p['laplace']['identified'][k]} | "
                 f"[{p['bands_log10']['lower'][k]:.2f}, {p['bands_log10']['upper'][k]:.2f}] |")
    L += ["", "## The rows", "", "| row | printed | model | residual (dex) | decisive |", "|---|---|---|---|---|"]
    for r in p["rows"]:
        L.append(f"| {r['id']} | {r['target']:.3g} +/- {r['half_width']:.2g} | {r['predicted']:.3g} | {r['residual_dex']:+.2f} | {r['decisive']} |")
    d = p["diagnostics"]
    L += ["", "## Nguyen's pot, integrated (LYSP 16, glucose 150 mmol/L, pH 6.8, 30 min), mmol/L", ""]
    for t, v in d["nguyen_pot_30min"].items():
        L.append(f"- {t}: CML {v['CML']:.4f} (printed range {d['nguyen_printed_cml_range_mmol_l']}), CEL {v['CEL']:.4f}, fructosyl-lysine {v['FLP']:.3f}, "
                 f"lysine lost {100 * v['lysine_lost_fraction']:.1f} %")
    b = d["berk2021_flp_to_cml_180C"]
    L += ["", f"- Berk 2021, fructosyl-lysine -> CML at 180 C (dry sesame): model {b['model_per_min']:.3g} vs {b['berk_per_min']:.3g} per minute ({b['dex']:+.2f} dex)"]
    for t, v in d["hamzalioglu2026_laclys_to_cml"].items():
        L.append(f"- Hamzalioglu 2026, lactulosyl-lysine -> CML at {t} (milk): model {v['model_per_min']:.3g} vs {v['milk_per_min']:.3g} ({v['dex']:+.2f} dex)")
    L.append(f"- {d['troise2015_direction']}")
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-nfev", type=int, default=200)
    args = parser.parse_args(argv)
    assert PREREG.exists(), "B20 is pre-registered; write the prereg before running the fit"
    payload = build(args.max_nfev)
    OUT_JSON.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    OUT_MD.write_text(render(payload), encoding="utf-8")
    print(f"wrote {OUT_JSON} and {OUT_MD}: cost {payload['objective']['final_cost']:.3f}; "
          f"frozen {json.dumps(payload['frozen_parameters']['glycation'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
