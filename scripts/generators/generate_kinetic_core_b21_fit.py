#!/usr/bin/env python
"""
Build Wave B21 -- THE AQUEOUS GLUCOSONE ROUTE TO GLYOXAL (2026-09-09).

Pre-registered in ``results/validation/kinetic_core_b21_prereg.md`` before any B21 number existed.
One step entered the trunk network (network.AQUEOUS_GLYOXAL_REACTIONS: the Amadori compound ->
glucosone + glycine) and the glucosone -> glyoxal constant received an aqueous value that replaces
the B13 glass value in the operative set.

WHAT IS FITTED: two coordinates, log10 k_ama_g and log10 k_g_go at the trunk's 100 C reference, on
Hamzalioglu 2026's first-order constants (whole milk, lactulosyl-lysine, 110-140 C): four for the
glucosone formation, two determinate ones for glucosone -> glyoxal; barriers declared (75.9 and 4.2
kJ/mol, that laboratory's measured values). No integration in the objective. Diagnostics: Quan
2020's glyoxal level, Xia 2022's glyoxal-to-methylglyoxal ordering, Leahy 1989's total pyrazine,
and the B1 browning hold-out re-scored through the frozen B1 hold-out generator, before and after.

Usage:
    python scripts/generators/generate_kinetic_core_b21_fit.py
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from dataclasses import replace
from datetime import date
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts" / "generators") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))

from src import data_paths  # noqa: E402
from src.kinetic_core import operative_parameters, trunk_conditions  # noqa: E402
from src.kinetic_core import parameters_dicarbonyl as PD  # noqa: E402
from src.kinetic_core.engine import b1_fitted  # noqa: E402
from src.kinetic_core.integrate import integrate  # noqa: E402
from src.kinetic_core.parameters import T_REF_K  # noqa: E402

WAVE = "B21"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b21_prereg.md"
OUT_JSON = data_paths.VALIDATION_DIR / "kinetic_core_b21_fit_report.json"
OUT_MD = data_paths.VALIDATION_DIR / "kinetic_core_b21_fit_report.md"
B18_SHIP = data_paths.VALIDATION_DIR / "kinetic_core_b18_ship_rule.json"
CELSIUS = 273.15
R_KJ = 8.314e-3
SEED = 20260909

# ===========================================================================
# 1. THE ROWS: Hamzalioglu 2026 Table 1 (value, half-width) per temperature, per minute
# ===========================================================================
HAM_ANCHOR = ("Hamzalioglu, Kocadagli & Gokmen 2026 (whole milk, 110-140 C, multiresponse fit), Table 1 steps 4 and 10, Table 2; "
              "hamzalioglu2026_extraction.md sec. 4")
ROWS_RAW = [
    # id, key, T_C, value, half-width (None = not printed / narrow), decisive
    ("ham_k_ama_g_110C", "k_ama_g", 110.0, 1.9e-2, None, True),
    ("ham_k_ama_g_120C", "k_ama_g", 120.0, 2.1e-2, None, True),
    ("ham_k_ama_g_130C", "k_ama_g", 130.0, 2.5e-2, 7.1e-2, False),
    ("ham_k_ama_g_140C", "k_ama_g", 140.0, 1.5e-1, None, True),
    ("ham_k_g_go_120C", "k_g_go", 120.0, 3.3e-1, None, True),
    ("ham_k_g_go_130C", "k_g_go", 130.0, 3.5e-1, 3.7e-1, False),
]
SIGMA_DEFAULT = 0.15     # the paper prints no interval for these four; a five-temperature, 30 K fit
SIGMA_WIDE = 0.5
EA_OF = {"k_ama_g": PD.EA_AMA_G_KJ_MOL, "k_g_go": PD.EA_G_GO_AQUEOUS_KJ_MOL}
KEYS: Tuple[str, ...] = PD.AQUEOUS_GLYOXAL_COORDINATES
KEY_OF_COORD = {"log10_k_ama_g_100C": "k_ama_g", "log10_k_g_go_aqueous_100C": "k_g_go"}


def rows() -> Tuple[Dict[str, Any], ...]:
    out: List[Dict[str, Any]] = []
    for rid, key, t_c, value, half, decisive in ROWS_RAW:
        sigma = SIGMA_WIDE if half is not None and half >= value else SIGMA_DEFAULT
        out.append(dict(id=rid, key=key, t_c=t_c, target=value, half_width=half, sigma_log=sigma, decisive=decisive, anchor=HAM_ANCHOR,
                        kind="rate_constant"))
    return tuple(out)


ROWS = rows()
PRIOR = np.array([PD.FROZEN_B21[k] for k in KEYS], dtype=float)
LOWER, UPPER = PRIOR - 2.0, PRIOR + 2.0


def k_at(log10_k100: float, ea: float, t_c: float) -> float:
    t_k = t_c + CELSIUS
    return 10.0 ** log10_k100 * math.exp(-ea / R_KJ * (1.0 / t_k - 1.0 / T_REF_K))


def predictions(x: np.ndarray) -> Dict[str, float]:
    coord = dict(zip(KEYS, (float(v) for v in x)))
    return {r["id"]: k_at(coord[next(c for c, k in KEY_OF_COORD.items() if k == r["key"])], EA_OF[r["key"]], r["t_c"]) for r in ROWS}


def residuals(x: np.ndarray) -> np.ndarray:
    pred = predictions(x)
    return np.array([math.log10(pred[r["id"]] / float(r["target"])) / float(r["sigma_log"]) for r in ROWS])


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
            "identified": {k: bool(s < 1.0 and not b) for k, s, b in zip(KEYS, sig, on_bound)}, "on_bound": dict(zip(KEYS, on_bound))}


# ===========================================================================
# 2. DIAGNOSTICS, before (glass) and after (the candidate)
# ===========================================================================
def _install(x) -> None:
    """Point the operative set at the candidate (x) or at the glass 'before' (x is None)."""
    if x is None:
        before = dict(PD.with_aqueous_glyoxal(-30.0, -30.0))
        before["k_ama_g"] = replace(before["k_ama_g"], k_ref=0.0)
        before["k_g_go"] = PD.GLASS_K_G_GO
        PD.AQUEOUS_GLYOXAL_PARAMETERS = before
    else:
        PD.AQUEOUS_GLYOXAL_PARAMETERS = PD.with_aqueous_glyoxal(float(x[0]), float(x[1]))


def _pot(initial, t_c, minutes, ph, points=8):
    params, _ = trunk_conditions.apply(dict(operative_parameters(b1_fitted())), SimpleNamespace(ph=float(ph), water_activity=None))
    grid = np.linspace(0.0, minutes, points)
    return integrate(params, t_c + CELSIUS, initial, grid, rtol=1e-8, atol=1e-14)


def _end(run, key):
    return float(run.series(key)[-1])


def diagnostics_for(x) -> Dict[str, Any]:
    _install(x)
    out: Dict[str, Any] = {}
    quan = {}
    for t_c, lo, hi in ((100.0, 0.052, 0.127), (130.0, 0.144, 0.605)):
        run = _pot({"Glc": 100.0, "Gly": 30.0}, t_c, 21.0, 7.0)
        go = _end(run, "GO")
        quan[f"{int(t_c)}C"] = {"model_go_mmol_l": go, "printed_range": [lo, hi],
                                "within_range_widened_0.5_dex": bool(lo / 10 ** 0.5 <= go <= hi * 10 ** 0.5),
                                "dex_from_range": (0.0 if lo <= go <= hi else (math.log10(go / hi) if go > hi else math.log10(go / lo)))}
    out["quan2020_glyoxal"] = quan
    xia = _pot({"Glc": 200.0, "Gly": 200.0}, 130.0, 80.0, 7.5)
    out["xia2022_ordering_130C_80min"] = {"go_mmol_l": _end(xia, "GO"), "mgo_mmol_l": _end(xia, "MGO"), "go_above_mgo": bool(_end(xia, "GO") > _end(xia, "MGO"))}
    leahy = _pot({"Glc": 100.0, "Gly": 100.0}, 95.0, 120.0, 9.0, points=13)
    mw = {"PZ": 80.09, "MPZ": 94.12, "DMP": 108.14}
    total = sum(_end(leahy, k) * mw[k] * 1000.0 for k in mw)        # ug/L
    out["leahy1989_total_pyrazine_95C_2h"] = {"model_ug_per_l": total, "leahy_ug_per_l": 13100.0,
                                             "dex": math.log10(total / 13100.0) if total > 0 else float("-inf"),
                                             "b18_dex": json.loads(B18_SHIP.read_text())["leahy"]["T4_total"]["dex"] if B18_SHIP.exists() else None}
    import generate_kinetic_core_b1_holdout as B1H

    fit_payload = json.loads((data_paths.VALIDATION_DIR / "kinetic_core_b1_fit_report.json").read_text(encoding="utf-8"))
    frozen = fit_payload["frozen_parameters"]["variant_B_reactant_side_sink"]
    scored = B1H.score_variant("B (reactant-side sink, out-of-sample)", frozen, B1H.load_holdout_series())
    out["b1_browning_holdout"] = {k: scored[k] for k in ("median_fold_error", "max_fold_error", "fraction_within_2x", "fraction_within_3x", "n_holdout_points")}
    return out


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
    before = diagnostics_for(None)
    after = diagnostics_for(x)
    _install(x)
    frozen = dict(zip(KEYS, x.tolist()))
    return {
        "artifact": "kinetic_core_b21_fit_report",
        "wave": f"{WAVE} -- the aqueous glucosone route to glyoxal",
        "generated_on": date.today().isoformat(), "generated_by": "scripts/generators/generate_kinetic_core_b21_fit.py", "git": _git_head(),
        "prereg": data_paths.rel(PREREG), "declaration": "docs/reference/FIT_HOLDOUT_DECLARATION.md Amendment 31",
        "objective": {"form": "six printed first-order constants (Hamzalioglu 2026, 110-140 C), log10 residuals; two log10 constants at 100 C free, barriers declared",
                      "n_rows": len(ROWS), "n_free_parameters": len(KEYS), "final_cost": best["cost"], "best_start": best["start"], "reduced_chi2": lap["chi2_reduced"]},
        "rows": [dict(r, predicted=pred[r["id"]], residual_dex=best["residuals_dex"][r["id"]]) for r in ROWS],
        "members": members,
        "frozen_parameters": {"aqueous_glyoxal": frozen, "declared_barriers_kj_mol": EA_OF,
                              "glass_k_g_go_at_100C": PD.GLASS_K_G_GO.k_ref, "glass_k_g_go_ea": PD.GLASS_K_G_GO.ea_kj_mol,
                              "constants_at_120C": {k: k_at(frozen[c], EA_OF[k], 120.0) for c, k in KEY_OF_COORD.items()}},
        "bands_log10": {"lower": dict(zip(KEYS, LOWER.tolist())), "upper": dict(zip(KEYS, UPPER.tolist()))},
        "laplace": lap,
        "diagnostics": {"before_glass": before, "after_candidate": after},
    }


def render(p: Dict[str, Any]) -> str:
    L = [f"# {p['wave']}", "", f"*Generated {p['generated_on']} by `{p['generated_by']}`; pre-registration `{p['prereg']}`; {p['declaration']}.*", "",
         f"Objective: {p['objective']['form']}. Cost {p['objective']['final_cost']:.3f} on {p['objective']['n_rows']} rows, reduced chi-square {p['objective']['reduced_chi2']:.2f}.", "",
         "## The fitted constants (log10 at 100 C; barriers declared)", "", "| coordinate | value | sigma (dex) | identified |", "|---|---|---|---|"]
    for k in KEYS:
        L.append(f"| {k} | {p['frozen_parameters']['aqueous_glyoxal'][k]:.4f} | {p['laplace']['sigma'][k]:.3f} | {p['laplace']['identified'][k]} |")
    L += [f"", f"The glass value of k_g_go at 100 C was {p['frozen_parameters']['glass_k_g_go_at_100C']:.3g} per minute (Ea {p['frozen_parameters']['glass_k_g_go_ea']}); "
          f"the aqueous value is {10 ** p['frozen_parameters']['aqueous_glyoxal']['log10_k_g_go_aqueous_100C']:.3g}.", "",
          "## The rows", "", "| row | printed | model | residual (dex) | decisive |", "|---|---|---|---|---|"]
    for r in p["rows"]:
        L.append(f"| {r['id']} | {r['target']:.3g} | {r['predicted']:.3g} | {r['residual_dex']:+.2f} | {r['decisive']} |")
    for tag in ("before_glass", "after_candidate"):
        d = p["diagnostics"][tag]
        L += ["", f"## Diagnostics, {tag.replace('_', ' ')}", ""]
        for t, v in d["quan2020_glyoxal"].items():
            L.append(f"- Quan 2020 glyoxal at {t}, 21 min: model {v['model_go_mmol_l']:.4f} mmol/L, printed {v['printed_range']} ({v['dex_from_range']:+.2f} dex from the range)")
        x = d["xia2022_ordering_130C_80min"]
        L.append(f"- Xia 2022 at 130 C, 80 min: glyoxal {x['go_mmol_l']:.3f} vs methylglyoxal {x['mgo_mmol_l']:.3f} mmol/L; glyoxal above: {x['go_above_mgo']}")
        le = d["leahy1989_total_pyrazine_95C_2h"]
        L.append(f"- Leahy 1989 total pyrazine, 95 C, 2 h: model {le['model_ug_per_l']:.3g} ug/L vs 13100 ({le['dex']:+.2f} dex; B18 recorded {le['b18_dex']:+.2f})")
        b = d["b1_browning_holdout"]
        L.append(f"- B1 browning hold-out: median fold {b['median_fold_error']:.2f}, max {b['max_fold_error']:.2f}, within 2x {b['fraction_within_2x']:.2f}, within 3x {b['fraction_within_3x']:.2f}")
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-nfev", type=int, default=200)
    args = parser.parse_args(argv)
    assert PREREG.exists(), "B21 is pre-registered; write the prereg before running the fit"
    payload = build(args.max_nfev)
    OUT_JSON.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    OUT_MD.write_text(render(payload), encoding="utf-8")
    print(f"wrote {OUT_JSON}: cost {payload['objective']['final_cost']:.3f}; frozen {json.dumps(payload['frozen_parameters']['aqueous_glyoxal'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
