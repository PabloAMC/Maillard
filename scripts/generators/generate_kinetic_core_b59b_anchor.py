"""
WAVE B59b V2 / V3 (kinetic_core_b59b_prereg.md sec. 3 item 4).

V2 ANCHOR: the one unmeasured constant of de Roos's route, the TTCA -> cysteine Amadori formation rate at 90 C, set by
reproducing Cao 2017's blank pot in the engine (0.4 M Cys + Xyl, 0.2 M phosphate pH 5.5, 90 C, 0.5-5 h; Cys-Amadori
0.006 / 0.063 / 0.090 / 0.106 / 0.109 / 0.133 mol/L): one-dimensional least squares in log10 k on the Cys-Amadori
series alone. TTCA, free cysteine and the conserved cysteine pool are then CHECKS. V3: a pot of the cysteine
Amadori compound alone at pH 5 and 6, 100 C, must decay at de Roos's 4.5e-3 / 1.4e-2 1/min (clamped pH).
Writes results/validation/kinetic_core_b59b_anchor.json.

    python scripts/generators/generate_kinetic_core_b59b_anchor.py
"""
from __future__ import annotations

import dataclasses
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402
from src.kinetic_core import engine as E  # noqa: E402
from src.kinetic_core import parameters_sulfur as ps  # noqa: E402
from src.kinetic_core import ph_state as PH  # noqa: E402
from src.kinetic_core import sulfur as S  # noqa: E402

OUT = data_paths.VALIDATION_DIR / "kinetic_core_b59b_anchor.json"
TIMES = [30.0, 60.0, 120.0, 180.0, 240.0, 300.0]
CAO_BLANK = {  # mol/L, cao2017_extraction.md digitised Figs. 4-5
    "CAMA": [0.006, 0.063, 0.090, 0.106, 0.109, 0.133],
    "TTCA": [0.356, 0.338, 0.290, 0.268, 0.254, 0.236],
    "Cys": [0.076, 0.068, 0.066, 0.064, 0.062, 0.056],
}
BUFFER = PH.BufferSpec(kind="phosphate", phosphate_mol_l=0.2, declared=True, source="Cao 2017 blank pot")


def cao_run(parameters, drift, log10_k):
    ps.CYS_AMADORI_ROUTE = dataclasses.replace(ps.CYS_AMADORI_ROUTE, enabled=True, k_form_ref_per_min=10.0 ** log10_k)
    run = S.integrate_sulfur(parameters, 363.15, {"Cys": 400.0, "PENT": 400.0}, [0.0] + TIMES, ph=5.5,
                             buffer_spec=BUFFER, ph_drift=drift, rtol=1e-8, atol=1e-12)
    return {k: [float(run.series(k)[i + 1]) / 1000.0 for i in range(len(TIMES))] for k in ("CAMA", "TTCA", "Cys")}, run


def main() -> int:
    S.PH_COUPLING = "inline"
    parameters = E.core_parameters(E.SULFUR)
    drift = E.core_ph_drift()
    saved = ps.CYS_AMADORI_ROUTE

    def loss(log10_k):
        pred, _ = cao_run(parameters, drift, log10_k)
        return sum((math.log10(max(p, 1e-9)) - math.log10(m)) ** 2 for p, m in zip(pred["CAMA"][1:], CAO_BLANK["CAMA"][1:]))

    # the 0.5 h point (0.006) is the noisiest read and the steepest; anchored on the five later points, reported in full
    best = minimize_scalar(loss, bounds=(-6.0, -1.0), method="bounded", options={"xatol": 1e-3})
    pred, run = cao_run(parameters, drift, best.x)
    pool = [pred["CAMA"][i] + pred["TTCA"][i] + pred["Cys"][i] for i in range(len(TIMES))]
    measured_pool = [sum(CAO_BLANK[k][i] for k in CAO_BLANK) for i in range(len(TIMES))]
    checks = {k: [round(math.log10(max(p, 1e-9) / m), 3) for p, m in zip(pred[k], CAO_BLANK[k])] for k in CAO_BLANK}
    # V3: de Roos's decay, CAMA alone, clamped pH, 100 C
    v3 = {}
    for ph in (5.0, 6.0):
        r = S.integrate_sulfur(parameters, 373.15, {"CAMA": 10.0}, [0.0, 10.0], ph=ph,
                               buffer_spec=PH.BufferSpec(kind="clamped"), rtol=1e-10, atol=1e-14)
        v3[str(ph)] = -math.log(r.final("CAMA") / 10.0) / 10.0
    ps.CYS_AMADORI_ROUTE = saved
    report = {
        "artifact": "kinetic_core_b59b_anchor", "prereg": data_paths.rel(data_paths.VALIDATION_DIR / "kinetic_core_b59b_prereg.md"),
        "k_form_ref_per_min_90C": 10.0 ** best.x, "log10_k": best.x, "anchor_loss_sum_sq_dex": best.fun,
        "predicted_mol_l": pred, "measured_mol_l": CAO_BLANK, "log10_pred_over_meas": checks,
        "cysteine_pool_pred": pool, "cysteine_pool_measured": measured_pool,
        "final_ph_in_situ": float(run.ph_final),
        "v3_decay_per_min": v3, "v3_measured": {"5.0": 4.5e-3, "6.0": 1.4e-2},
        "v3_pass": all(abs(math.log10(v3[k] / m)) < 0.02 for k, m in (("5.0", 4.5e-3), ("6.0", 1.4e-2))),
    }
    OUT.write_text(json.dumps(report, indent=2))
    print(json.dumps({k: report[k] for k in ("k_form_ref_per_min_90C", "anchor_loss_sum_sq_dex", "log10_pred_over_meas",
                                              "cysteine_pool_pred", "v3_decay_per_min", "v3_pass")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
