#!/usr/bin/env python
"""
Wave B24b -- THE BRANCH THAT REFUSED B24 (2026-09-10).

`results/validation/kinetic_core_b24b_prereg.md`. Two free coordinates on four rows:

  the SHAPE   Hofmann & Schieberle 1998b Table 9's three AP : ATHP molar ratios, 0.16 / 0.51 / 12.8
              across a hundredfold methylglyoxal ladder against 400 mmol/L proline. A ratio inside
              one analysis, so the response factor and the extraction cancel -- which is the whole
              reason the target is a ratio and not a level.
  the SINK    Table 7 experiment 3: 1-pyrroline in fivefold excess over methylglyoxal gives 0.33
              mol % of the methylglyoxal, where B24 made 41.

B24's two constants are HELD at their frozen optimum. The question is whether adding the branch
fixes the shape, and refitting them would let the fit buy the shape with the constants that already
worked. Their two fed rows stay in as CHECKS: the new sinks change the 1-pyrroline pool, so their
predictions move even though their constants do not.

The pH ladder is a CHECK too, not a fit row -- see the pre-registration's section 3b. Fitting it
would need a piecewise slope of its own on two independent ratios, which is a coordinate per data
point. The step takes B18's pH term instead and the ladder scores that transfer out of sample.
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
from src.kinetic_core.parameters_proline import (  # noqa: E402
    AP_ATHP_SWITCH, AP_ATHP_SWITCH_SOURCE, ATHP_PH_LADDER_CENSORED, ATHP_PH_LADDER_UG,
    ATHP_PH_SOURCE, FROZEN_B24, PROLINE_B24B_COORDINATES, PROLINE_COORDINATES, PROLINE_FIT_PH,
    PYRL_EXCESS_ROW, with_fitted_proline, with_fitted_proline_b24b,
)

WAVE = "B24b"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b24b_prereg.md"
OUT_JSON = data_paths.VALIDATION_DIR / "kinetic_core_b24b_fit_report.json"
OUT_MD = data_paths.VALIDATION_DIR / "kinetic_core_b24b_fit_report.md"
CELSIUS = 273.15
SEED = 20260910
MINUTES = 30.0
PRO_MMOL_L = 400.0
KEYS: Tuple[str, ...] = PROLINE_B24B_COORDINATES
#: Prior centres. The branch constant starts at the acylation's fitted value, because the two steps
#: are condensations of the same 1-pyrroline with a small carbonyl in the same pot; the loss starts
#: three decades below it, which is roughly what "a slow drain" means against a 30-minute cook.
PRIOR = np.array([FROZEN_B24["log10_k_pyrl_ap_100C"], FROZEN_B24["log10_k_pyrl_ap_100C"] - 3.0], dtype=float)
LOWER, UPPER = PRIOR - 4.0, PRIOR + 4.0
BASE = dict(operative_parameters(b1_fitted()))
B24_FROZEN_ARGS = [FROZEN_B24[k] for k in PROLINE_COORDINATES]


def parameters_for(x):
    p = dict(BASE)
    p.update(with_fitted_proline(*B24_FROZEN_ARGS))          # HELD, not refitted
    p.update(with_fitted_proline_b24b(float(x[0]), float(x[1])))
    return p


def _run(x, initial, ph=PROLINE_FIT_PH, t_c=100.0, minutes=MINUTES):
    params, _ = trunk_conditions.apply(parameters_for(x), SimpleNamespace(ph=float(ph), water_activity=None))
    return integrate(params, t_c + CELSIUS, initial, np.array([0.0, float(minutes)]), rtol=1e-8, atol=1e-16)


def _switch_rows(x) -> Dict[str, float]:
    """The AP : ATHP molar ratio at each rung of Table 9's methylglyoxal ladder."""
    out = {}
    for mgo, _target in sorted(AP_ATHP_SWITCH.items()):
        run = _run(x, {"PRO": PRO_MMOL_L, "Gly": PRO_MMOL_L, "MGO": float(mgo)})
        ap = float(run.series("AP")[-1])
        athp = float(run.series("ATHP")[-1])
        out[f"switch_mgo{mgo:g}"] = (ap + 1e-30) / (athp + 1e-30)
    return out


def _excess_row(x) -> float:
    r = PYRL_EXCESS_ROW
    run = _run(x, {"PYRL": r["pyrl_mmol_l"], "MGO": r["mgo_mmol_l"]})
    return 100.0 * float(run.series("AP")[-1]) / r["mgo_mmol_l"]


def rows() -> Tuple[Dict[str, Any], ...]:
    out = []
    for mgo, target in sorted(AP_ATHP_SWITCH.items()):
        out.append(dict(id=f"switch_mgo{mgo:g}", kind="within_study_ratio", target=float(target),
                        sigma_log=0.20, anchor=AP_ATHP_SWITCH_SOURCE))
    out.append(dict(id="pyrl_excess_ap_molpct", kind="fed_yield_mol_pct",
                    target=float(PYRL_EXCESS_ROW["ap_molpct_of_mgo"]), sigma_log=0.30,
                    anchor="Hofmann & Schieberle 1998b Table 7 experiment 3"))
    return tuple(out)


ROWS = rows()


def predictions(x) -> Dict[str, float]:
    out = dict(_switch_rows(x))
    out["pyrl_excess_ap_molpct"] = _excess_row(x)
    return out


def residuals(x) -> np.ndarray:
    pred = predictions(x)
    return np.array([
        max(-1e3, min(1e3, math.log10((pred[r["id"]] + 1e-30) / r["target"]) / r["sigma_log"]))
        for r in ROWS
    ])


def fit_member(start: int, max_nfev: int) -> Dict[str, Any]:
    from scipy.optimize import least_squares

    if start == 0:
        x0 = PRIOR.copy()
    else:
        rng = np.random.default_rng(SEED + start)
        x0 = np.clip(PRIOR + rng.normal(0.0, 0.7, size=PRIOR.shape), LOWER, UPPER)
    sol = least_squares(residuals, x0, bounds=(LOWER, UPPER), method="trf", max_nfev=max_nfev,
                        xtol=1e-10, ftol=1e-10)
    r = residuals(sol.x)
    return {"start": start, "x0": x0.tolist(), "x": sol.x.tolist(), "cost": float(np.sum(r * r)),
            "nfev": int(sol.nfev),
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
            "identified": {k: bool(s < 1.0 and not b) for k, s, b in zip(KEYS, sig, on_bound)},
            "on_bound": dict(zip(KEYS, on_bound))}


def checks(x) -> Dict[str, Any]:
    """Everything scored and NOT fitted: B24's two fed rows, and the pH transfer."""
    fed = {}
    for rid, initial, target in (("hof_t7_e1_pyrl2_mgo10", {"PYRL": 2.0, "MGO": 10.0}, 28.7),
                                 ("hof_t7_e2_pyrl2_mgo2", {"PYRL": 2.0, "MGO": 2.0}, 5.3)):
        run = _run(x, initial)
        pred = 100.0 * float(run.series("AP")[-1]) / initial["PYRL"]
        fed[rid] = {"target_mol_pct": target, "predicted_mol_pct": pred,
                    "dex": math.log10((pred + 1e-30) / target)}
    # the pH ladder, as ratios to the pH-7 rung, against B18's transferred term
    ref = None
    ladder = {}
    for ph, ug in sorted(ATHP_PH_LADDER_UG.items(), key=lambda kv: float(kv[0])):
        run = _run(x, {"PYRL": 2.0, "ACETOL": 2.0}, ph=float(ph))
        ladder[ph] = {"printed_ug": ug, "model_athp_mmol_l": float(run.series("ATHP")[-1])}
    ref = ladder.get("7.0", {}).get("model_athp_mmol_l") or 0.0
    ref_ug = ATHP_PH_LADDER_UG["7.0"]
    for ph, row in ladder.items():
        row["printed_ratio_to_ph7"] = row["printed_ug"] / ref_ug
        row["model_ratio_to_ph7"] = (row["model_athp_mmol_l"] / ref) if ref > 0 else float("nan")
        row["dex"] = (math.log10((row["model_ratio_to_ph7"] + 1e-30) / row["printed_ratio_to_ph7"])
                      if row["printed_ratio_to_ph7"] > 0 else float("nan"))
    return {"b24_fed_rows": fed, "ph_ladder_vs_b18_transfer": ladder,
            "ph_ladder_censored_rung": dict(ATHP_PH_LADDER_CENSORED),
            "ph_ladder_source": ATHP_PH_SOURCE,
            "note": ("Neither block is in the objective. The fed rows check that adding the branch did "
                     "not break what B24 already fitted; the ladder checks B18's pH slopes on a "
                     "chemistry they were not fitted on.")}


def build(max_nfev: int) -> Dict[str, Any]:
    members = [fit_member(s, max_nfev) for s in range(3)]
    best = min(range(len(members)), key=lambda i: members[i]["cost"])
    x = np.array(members[best]["x"], dtype=float)
    r = residuals(x)
    pred = predictions(x)
    return {
        "wave": WAVE,
        "artifact": "kinetic_core_b24b_fit_report",
        "provenance": provenance.provenance_block(
            "kinetic_core_b24b_fit_report",
            generated_by="scripts/generators/generate_kinetic_core_b24b_fit.py", inputs=[]),
        "prereg": data_paths.rel(PREREG),
        "held_not_refitted": {"reason": "B24's two constants; refitting them would let the fit buy "
                                        "the shape with the constants that already worked",
                              "values": dict(zip(PROLINE_COORDINATES, B24_FROZEN_ARGS))},
        "objective": {"rows": [dict(r) for r in ROWS], "n_rows": len(ROWS),
                      "final_cost": float(np.sum(r * r)), "n_free": len(KEYS)},
        "bounds": {k: [float(LOWER[i]), float(UPPER[i])] for i, k in enumerate(KEYS)},
        "members": members, "best_start": best,
        "frozen_parameters": {"proline_b24b": dict(zip(KEYS, x.tolist()))},
        "predicted_by_row": pred,
        "residual_by_row_dex": {row["id"]: float(v) * float(row["sigma_log"]) for row, v in zip(ROWS, r)},
        "laplace": laplace(x, r),
        "checks": checks(x),
        "reference_temperature_K": 100.0 + CELSIUS,
        "reference_ph": PROLINE_FIT_PH,
    }


def render(p: Dict[str, Any]) -> str:
    L = [f"# Wave {WAVE} fit report", "",
         f"Cost {p['objective']['final_cost']:.3f} on {p['objective']['n_rows']} rows, "
         f"{p['objective']['n_free']} free coordinates. B24's two constants held, not refitted.", "",
         "| row | target | predicted | residual (dex) |", "|---|---:|---:|---:|"]
    for row in p["objective"]["rows"]:
        rid = row["id"]
        L.append(f"| {rid} | {row['target']:.4g} | {p['predicted_by_row'][rid]:.4g} | "
                 f"{p['residual_by_row_dex'][rid]:+.3f} |")
    lap = p["laplace"]
    L += ["", f"Laplace: sigma {json.dumps({k: round(v, 3) for k, v in lap['sigma'].items()})}; "
              f"identified {lap['identified']}; on bound {lap['on_bound']}.", "",
          "## Checks, not fitted", "",
          "| check | printed | model | dex |", "|---|---:|---:|---:|"]
    for rid, row in p["checks"]["b24_fed_rows"].items():
        L.append(f"| {rid} | {row['target_mol_pct']:.4g} mol % | {row['predicted_mol_pct']:.4g} | "
                 f"{row['dex']:+.3f} |")
    for ph, row in sorted(p["checks"]["ph_ladder_vs_b18_transfer"].items()):
        L.append(f"| ATHP ratio to pH 7 at pH {ph} | {row['printed_ratio_to_ph7']:.4g} | "
                 f"{row['model_ratio_to_ph7']:.4g} | {row['dex']:+.3f} |")
    L += ["", f"> {p['checks']['note']}", ""]
    return "\n".join(L)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-nfev", type=int, default=200)
    args = parser.parse_args(argv)
    assert PREREG.exists(), "B24b is pre-registered; write the prereg before running the fit"
    payload = build(args.max_nfev)
    OUT_JSON.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    OUT_MD.write_text(render(payload), encoding="utf-8")
    print(f"cost {payload['objective']['final_cost']:.3f}; frozen "
          f"{json.dumps({k: round(v, 3) for k, v in payload['frozen_parameters']['proline_b24b'].items()})}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
