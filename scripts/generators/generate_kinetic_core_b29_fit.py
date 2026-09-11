#!/usr/bin/env python
"""
Wave B29 -- THE FIRST OXYGEN AXIS ON THE TRUNK (2026-09-10).

`results/validation/kinetic_core_b29_prereg.md`. Two coordinates on four within-study ratios:
Hofmann & Schieberle 2000b Table 2's Strecker aldehyde under argon, air and air + copper, from the
Amadori compound and from the sugar pot.

THE TEST IS THAT ONE MULTIPLIER HAS TO EXPLAIN TWO DIFFERENT RATIOS. The ARP pot gives 9.2x and the
sugar pot 3.5x, and both share `f(argon)`. They can only differ if the two pots already weight the
Amadori route and the direct sugar route to glucosone differently, by about the right amount -- and
nothing here was tuned to make that so, because those constants come from a hazelnut paper and a
milk paper. So this is close to an out-of-sample test of the trunk's own branching.

The lever goes on the two oxidative entries and nowhere else, because the source's COMPANION paper
feeds the dicarbonyls directly and finds the Strecker aldehyde oxygen-independent. The sensitivity
is upstream of the dicarbonyl by measurement, not by assumption.
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

WAVE = "B29"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b29_prereg.md"
OUT_JSON = data_paths.VALIDATION_DIR / "kinetic_core_b29_fit_report.json"
OUT_MD = data_paths.VALIDATION_DIR / "kinetic_core_b29_fit_report.md"
CELSIUS = 273.15
SEED = 20260910
T_C, MINUTES, PH = 100.0, 120.0, 7.0
#: 1 mmol of each precursor in 10 mL = 100 mmol/L.
CHARGE_MMOL_L = 100.0
KEYS: Tuple[str, ...] = ("log10_f_argon", "log10_f_air_cu")
PRIOR = np.array([math.log10(1.0 / 9.2), math.log10(2.5)], dtype=float)
LOWER, UPPER = np.array([-3.0, 0.0]), np.array([0.0, 2.0])
BASE = dict(operative_parameters(b1_fitted()))
#: The printed within-study ratios. The Strecker aldehyde only; the ACID is recorded and not fitted,
#: because the model has no Strecker acid and the companion paper shows it is a different step with
#: its own oxygen term.
TARGETS: Dict[str, float] = {
    "arp_air_over_argon": 9.2, "glc_air_over_argon": 3.5,
    "arp_aircu_over_air": 2.5, "glc_aircu_over_air": 1.9,
}
SIGMA_LOG = 0.15
SOURCE = ("Hofmann & Schieberle 2000b, J. Agric. Food Chem. 48:4301, Table 2; "
          "hofmann2000b_extraction.md")
#: THE OBSERVABLE, CORRECTED 2026-09-10 ON REVIEW. The first run used "AKG" alone, the Strecker
#: product of GLYOXAL. That is one dicarbonyl's Strecker, not the Strecker aldehyde: methylglyoxal's
#: Strecker gives AKM by a different route, and Hofmann's phenylacetaldehyde is the aldehyde
#: whichever dicarbonyl did the Strecker. Measuring AKG alone made the aldehyde look 100 % oxidative
#: -- because glyoxal comes only through glucosone -- and the wave concluded a non-oxidative route
#: was missing. It is not missing: AMA -> 1-deoxyosone -> methylglyoxal -> AKM is that route, and
#: with the observable summed over both Strecker products the non-oxidative share is 45 % in the
#: Amadori pot and 36 % in the sugar pot. The first conclusion was an artefact of the observable.
OBSERVABLE: Tuple[str, ...] = ("AKG", "AKM")


def _factors(x) -> Dict[str, float]:
    return {"air": 1.0, "argon": 10.0 ** float(x[0]), "air_cu": 10.0 ** float(x[1])}


def _run(x, initial, atmosphere) -> float:
    process = SimpleNamespace(ph=PH, water_activity=None, atmosphere=atmosphere)
    params, _ = trunk_conditions.apply(dict(BASE), process, atmosphere_factors=_factors(x))
    run = integrate(params, T_C + CELSIUS, initial, np.array([0.0, MINUTES]), rtol=1e-8, atol=1e-16)
    return sum(float(run.series(k)[-1]) for k in OBSERVABLE)


def predictions(x) -> Dict[str, float]:
    arp = {"AMA": CHARGE_MMOL_L, "Gly": CHARGE_MMOL_L}
    glc = {"Glc": CHARGE_MMOL_L, "Gly": CHARGE_MMOL_L}
    out = {}
    for tag, initial in (("arp", arp), ("glc", glc)):
        argon, air, aircu = (_run(x, initial, a) for a in ("argon", "air", "air_cu"))
        out[f"{tag}_air_over_argon"] = (air + 1e-30) / (argon + 1e-30)
        out[f"{tag}_aircu_over_air"] = (aircu + 1e-30) / (air + 1e-30)
    return out


def residuals(x) -> np.ndarray:
    pred = predictions(x)
    return np.array([
        max(-1e3, min(1e3, math.log10((pred[k] + 1e-30) / TARGETS[k]) / SIGMA_LOG))
        for k in TARGETS
    ])


def fit_member(start: int, max_nfev: int) -> Dict[str, Any]:
    from scipy.optimize import least_squares

    if start == 0:
        x0 = PRIOR.copy()
    else:
        rng = np.random.default_rng(SEED + start)
        x0 = np.clip(PRIOR + rng.normal(0.0, 0.4, size=PRIOR.shape), LOWER, UPPER)
    sol = least_squares(residuals, x0, bounds=(LOWER, UPPER), method="trf", max_nfev=max_nfev,
                        xtol=1e-12, ftol=1e-12)
    r = residuals(sol.x)
    return {"start": start, "x0": x0.tolist(), "x": sol.x.tolist(), "cost": float(np.sum(r * r)),
            "nfev": int(sol.nfev),
            "residuals_dex": {k: float(v) * SIGMA_LOG for k, v in zip(TARGETS, r)}}


def laplace(x, r) -> Dict[str, Any]:
    h = 0.02
    jac = np.zeros((len(TARGETS), len(KEYS)))
    for j in range(len(KEYS)):
        xp, xm = x.copy(), x.copy()
        xp[j] += h
        xm[j] -= h
        jac[:, j] = (residuals(xp) - residuals(xm)) / (2 * h)
    dof = max(len(TARGETS) - len(KEYS), 1)
    chi2_red = float(r @ r) / dof
    cov = np.linalg.pinv(jac.T @ jac) * max(chi2_red, 1.0)
    sig = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    on_bound = [bool(abs(x[i] - LOWER[i]) < 1e-3 or abs(UPPER[i] - x[i]) < 1e-3) for i in range(len(KEYS))]
    return {"sigma": dict(zip(KEYS, sig.tolist())), "chi2_reduced": chi2_red, "dof": dof,
            "identified": {k: bool(s < 1.0 and not b) for k, s, b in zip(KEYS, sig, on_bound)},
            "on_bound": dict(zip(KEYS, on_bound))}


def air_is_exactly_one() -> Dict[str, Any]:
    """T3: a pot that declares nothing, and a pot that declares air, must be bit-for-bit identical."""
    glc = {"Glc": CHARGE_MMOL_L, "Gly": CHARGE_MMOL_L}
    none_p = SimpleNamespace(ph=PH, water_activity=None)
    air_p = SimpleNamespace(ph=PH, water_activity=None, atmosphere="air")
    a, _ = trunk_conditions.apply(dict(BASE), none_p)
    b, _ = trunk_conditions.apply(dict(BASE), air_p, atmosphere_factors={"air": 1.0, "argon": 0.1})
    same = all(a[k].k_ref == b[k].k_ref for k in trunk_conditions.OXIDATIVE_ENTRY_STEPS)
    ra = integrate(a, T_C + CELSIUS, glc, np.array([0.0, MINUTES]), rtol=1e-8, atol=1e-16)
    rb = integrate(b, T_C + CELSIUS, glc, np.array([0.0, MINUTES]), rtol=1e-8, atol=1e-16)
    oa = sum(float(ra.series(k)[-1]) for k in OBSERVABLE)
    ob = sum(float(rb.series(k)[-1]) for k in OBSERVABLE)
    return {"parameters_identical": bool(same),
            "observable_identical": bool(oa == ob),
            "pass": bool(same)}


def build(max_nfev: int) -> Dict[str, Any]:
    members = [fit_member(s, max_nfev) for s in range(3)]
    best = min(range(len(members)), key=lambda i: members[i]["cost"])
    x = np.array(members[best]["x"], dtype=float)
    r = residuals(x)
    return {
        "wave": WAVE, "artifact": "kinetic_core_b29_fit_report",
        "provenance": provenance.provenance_block(
            "kinetic_core_b29_fit_report",
            generated_by="scripts/generators/generate_kinetic_core_b29_fit.py", inputs=[]),
        "prereg": data_paths.rel(PREREG), "source": SOURCE,
        "observable": list(OBSERVABLE),
        "lever_on": list(trunk_conditions.OXIDATIVE_ENTRY_STEPS),
        "objective": {"targets": dict(TARGETS), "sigma_log": SIGMA_LOG, "n_rows": len(TARGETS),
                      "n_free": len(KEYS), "final_cost": float(np.sum(r * r))},
        "bounds": {k: [float(LOWER[i]), float(UPPER[i])] for i, k in enumerate(KEYS)},
        "members": members, "best_start": best,
        "frozen_parameters": {"atmosphere": dict(zip(KEYS, x.tolist()))},
        "atmosphere_factors": _factors(x),
        "predicted": predictions(x),
        "residual_by_row_dex": {k: float(v) * SIGMA_LOG for k, v in zip(TARGETS, r)},
        "laplace": laplace(x, r),
        "air_is_exactly_one": air_is_exactly_one(),
        "not_fitted": {"strecker_acid_air_over_argon": {"arp": 10.0, "glc": 9.0},
                       "why": ("the model has no Strecker acid, and the companion paper shows the "
                               "ACID is oxygen-dependent even from FED dicarbonyls while the "
                               "aldehyde is not -- so it is a different step with its own oxygen "
                               "term, not this one")},
        "reference_temperature_K": T_C + CELSIUS, "reference_ph": PH,
    }


def render(p: Dict[str, Any]) -> str:
    L = [f"# Wave {WAVE} fit report", "",
         f"Cost {p['objective']['final_cost']:.3f} on {p['objective']['n_rows']} ratios, "
         f"{p['objective']['n_free']} free. Lever on {', '.join(p['lever_on'])}; observable "
         f"{' + '.join(p['observable'])} (the total Strecker flux, not one dicarbonyl's).", "",
         f"Fitted factors: {json.dumps({k: round(v, 4) for k, v in p['atmosphere_factors'].items()})}", "",
         "| ratio | printed | model | residual (dex) |", "|---|---:|---:|---:|"]
    for k, target in p["objective"]["targets"].items():
        L.append(f"| {k} | {target:.3g} | {p['predicted'][k]:.4g} | {p['residual_by_row_dex'][k]:+.3f} |")
    lap = p["laplace"]
    L += ["", f"Laplace: sigma {json.dumps({k: round(v, 3) for k, v in lap['sigma'].items()})}; "
              f"identified {lap['identified']}; on bound {lap['on_bound']}.", "",
          f"Air is exactly 1: parameters identical {p['air_is_exactly_one']['parameters_identical']}, "
          f"observable identical {p['air_is_exactly_one']['observable_identical']}.", ""]
    return "\n".join(L)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-nfev", type=int, default=200)
    args = parser.parse_args(argv)
    assert PREREG.exists(), "B29 is pre-registered; write the prereg before running the fit"
    payload = build(args.max_nfev)
    OUT_JSON.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    OUT_MD.write_text(render(payload), encoding="utf-8")
    print(f"cost {payload['objective']['final_cost']:.3f}; factors "
          f"{json.dumps({k: round(v, 3) for k, v in payload['atmosphere_factors'].items()})}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
