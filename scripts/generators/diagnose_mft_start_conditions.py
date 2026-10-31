#!/usr/bin/env python
"""
diagnose_mft_start_conditions.py -- dry check T1 of the 100 C MFT lag (thiol prereg sec. 14, 2026-10-31). No fit.

Does Schieberle 2000's early MFT excess come from the start conditions? Variants on one pot (ribose 100 + Cys 33 mM,
pH 5, 100 C): as fitted (isothermal, all Cys free); a 13-min linear come-up 20 -> 100 C; TTCA pre-equilibrated at
t = 0 (Cao 2017's K); both -- at R3En's optimum and at sec. 13's Cerny-passing point (Ea 160, DDP route +3).
Prints the MFT residuals (dex); the numbers are recorded in the prereg.

    python scripts/generators/diagnose_mft_start_conditions.py
"""
from __future__ import annotations

import functools
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
for entry in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

from src import data_paths  # noqa: E402
from src.kinetic_core import fit_rows, parameters_sulfur, sulfur  # noqa: E402

import build_submodel_probe as probe  # noqa: E402
import fit_joint_deoxyosone_thiol as joint  # noqa: E402

POT = "schieberle2000_rib_cys_ph5_100C"
COME_UP = ([0.0, 13.0, 720.0], [293.15, 373.15, 373.15])   # ~9 min dead time in the data; the autoclave's ramp


def ttca_equilibrated(charge):
    """Cysteine and pentose moved into TTCA at Cao 2017's apparent K, x / ((c - x)(s - x)) = K."""
    k = parameters_sulfur.TTCA_FORMATION.k_l_per_mmol if parameters_sulfur.TTCA_FORMATION.enabled else 0.0
    c, s = charge.get("Cys", 0.0), charge.get("PENT", 0.0)
    b = k * (c + s) + 1.0
    x = (b - math.sqrt(b * b - 4.0 * k * k * c * s)) / (2.0 * k) if k else 0.0
    out = dict(charge, Cys=c - x, PENT=s - x)
    out["TTCA"] = out.get("TTCA", 0.0) + x
    return out, x


def run(job0, rows, theta, ramp, ttca):
    keys = tuple(theta)
    j = dict(job0, **joint.trial([theta[k] for k in keys], keys), strecker_share=0.042)
    bound = 0.0
    if ttca:
        j["charge"], bound = ttca_equilibrated(job0["charge"])
    plain = probe.integrate_sulfur
    probe.integrate_sulfur = functools.partial(plain, temperature_profile=COME_UP) if ramp else plain
    try:
        res = joint.run_task(j)
    finally:
        probe.integrate_sulfur = plain
    if "failed" in res:
        return res["failed"], bound
    return {(r.species, r.t_min): math.log10(probe.predict(r, j, res) / r.value) for r in rows}, bound


def main() -> int:
    sulfur.PH_COUPLING = "inline"
    joint.DDP_LAW_STATIC = {"law": "none"}   # R3En's law
    fit = json.loads((data_paths.VALIDATION_DIR / "kinetic_core_joint_deoxyosone_thiol_fit_s042_r3en.json").read_text())
    theta0 = dict(fit["fit"]["theta"])
    pots, rows = fit_rows.load("thiol_submodel")
    job0 = next(j for j in probe.plan(pots, rows) if j["pot"] == POT)
    mft = sorted((r for r in rows if r.pot == POT and r.species == "MFT"), key=lambda r: r.t_min)
    meas = {r.t_min: r.value for r in mft}
    print("measured MFT log10 ratios: 60/30 %+.2f, 720/360 %+.2f, 60/720 %+.2f" % (
        math.log10(meas[60.0] / meas[30.0]), math.log10(meas[720.0] / meas[360.0]), math.log10(meas[60.0] / meas[720.0])))
    points = {"R3En optimum": theta0,
              "sec. 13 Cerny-passing (Ea 160, DDP +3)": dict(theta0, ea_k_nf_mft=160.0,
                                                              offset_ddp_route=theta0["offset_ddp_route"] + 3.0)}
    for name, theta in points.items():
        for ramp in (False, True):
            for ttca in (False, True):
                vals, bound = run(job0, mft, theta, ramp, ttca)
                if isinstance(vals, str):
                    print(f"{name} | ramp={ramp} ttca={ttca} | FAILED {vals}")
                    continue
                resid = " ".join(f"{t:.0f}:{vals[('MFT', t)]:+.2f}" for t in (30.0, 60.0, 360.0, 720.0))
                print(f"{name} | ramp={ramp!s:5} ttca={ttca!s:5} | MFT residuals {resid}"
                      + (f" | TTCA at t=0 {bound:.1f} mM" if ttca else ""), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
