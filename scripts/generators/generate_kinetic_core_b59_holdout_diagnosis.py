"""
WAVES B59a-r and B59b: the pre-declared hold-out diagnosis (not gating; kinetic_core_b59a_prereg.md B59a-r and
kinetic_core_b59b_prereg.md sec. 4), at every optimum of one B64 thermal reading with its own switches.

Arms: the reading's baseline optimum (both switches off), B59a-r at each share (Strecker sulfur on), B59b (de Roos's
route on) without and with the Strecker sulfur at each share; an arm is skipped until its member exists.
Hold-outs, none of them a row of the objective:
  - Hofmann 1998 ribose + Cys pH 3 / 5 / 7 (145 C along the reading's profile): MFT, FFT log10 pred/meas; H2S made.
  - Schieberle 2000 (same pot, 100 C): MFT 100 C-6 h / 145 C factor (measured 7.9-13) and the MFT shape over 30 min.
  - US 10,076,126 buffer arm (Xyl + Cys 0.25 M, pH 5.5, 110 C): MFT at 0.5 / 1 / 2 / 4 h over 4 h
    (measured 0.008 / 0.067 / 0.31 / 1). The buffer is not printed; 0.2 M phosphate is the sibling patent's (declared).
  - Li 2011 (Glc + Cys) is NOT evaluated: the sulfur lane carries no hexose species, so the declared hexose check has
    no pot to run in; recorded as such in the report.
  - Zhang 2025b blank (Xyl + Cys 40 mM, 0.2 M phosphate pH 6.5, 120 C): TTCA and Cys-Amadori at 10-50 min.
Every pot carries the Hofmann systems' ambient dissolved oxygen (closed vessel, air headspace).
Writes results/validation/kinetic_core_b59_holdout_diagnosis_<reading>.json.

    python scripts/generators/generate_kinetic_core_b59_holdout_diagnosis.py --reading isothermal
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for p in (ROOT, ROOT / "scripts" / "generators"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import generate_kinetic_core_b64_fit as B64  # noqa: E402  (inline pH coupling, B60's objective)
from generate_kinetic_core_b59a_fit import best_member  # noqa: E402
from src import data_paths  # noqa: E402
from src.kinetic_core import parameters_sulfur as ps  # noqa: E402
from src.kinetic_core import ph_state, sulfur  # noqa: E402

CELSIUS = 273.15
OX = {k: v for k, v in B64.B23.SYSTEMS["hofmann_pentose_pH5"]["initial"].items() if k == "OX"}
HOFMANN = {"PENT": 100.0, "Cys": 33.0, **OX}
HOFMANN_BUFFER = ph_state.BufferSpec(kind="phosphate", phosphate_mol_l=0.5, declared=True, source="Hofmann 1998 pot")
HOFMANN_PPB = {3.0: {"MFT": 553.0, "FFT": 229.0}, 5.0: {"MFT": 198.0, "FFT": 121.0}, 7.0: {"MFT": 25.0, "FFT": 12.0}}
SCHIEBERLE_MFT_UG = {30.0: 4.5, 60.0: 13.8, 360.0: 156.0, 720.0: 179.0}
NESTLE_MFT_NG_G = {30.0: 15.0, 60.0: 124.0, 120.0: 575.0, 240.0: 1855.0}
NESTLE_BUFFER = ph_state.BufferSpec(kind="phosphate", phosphate_mol_l=0.2, declared=True,
                                    source="US 8,920,862's buffer (US 10,076,126 does not print it)")
ZHANG_BUFFER = ph_state.BufferSpec(kind="phosphate", phosphate_mol_l=0.2, declared=True, source="Zhang 2025b blank")
ZHANG_MM = {"TTCA": {10.0: 11.2, 20.0: 15.1, 30.0: 17.7, 50.0: 18.2},
            "CAMA": {10.0: 0.07, 20.0: 1.38, 30.0: 1.60, 50.0: 1.61}}
MW_MFT = 114.17


def _vector(path: Path) -> np.ndarray:
    payload = json.loads(path.read_text())
    return np.array(payload["x_full"] if "x_full" in payload else B64.frozen_vector(payload), dtype=float)


def _hot_profile(reading: str):
    points = B64.READINGS[reading].get("hofmann_pentose_pH5", ((0.0, 145.0), (20.0, 145.0)))
    return np.array([t for t, _ in points]), np.array([c for _, c in points]) + CELSIUS


def _run(parameters, drift, t_c, charge, times, ph, buffer, profile=None):
    return sulfur.integrate_sulfur(parameters, t_c + CELSIUS, charge, [0.0] + list(times), ph=ph, buffer_spec=buffer,
                                   ph_drift=drift, rtol=1e-8, atol=1e-14, temperature_profile=profile)


def _log(p, m):
    return math.log10(max(p, 1e-15) / m)


def diagnose(x: np.ndarray, reading: str) -> dict:
    _f, _e, _d, drift = B64.B23.unpack(x)
    parameters = B64.B23.build_parameters(x)
    prof = _hot_profile(reading)
    out = {"hofmann_pH_sweep_log10_pred_over_meas": {}, "hofmann_h2s_made_mM": {}}
    hot = {}
    for ph, meas in HOFMANN_PPB.items():
        run = _run(parameters, drift, 145.0, HOFMANN, [float(prof[0][-1])], ph, HOFMANN_BUFFER, profile=prof)
        hot[ph] = run.final("MFT")
        out["hofmann_pH_sweep_log10_pred_over_meas"][str(ph)] = {
            s: _log(run.final(s) * MW_MFT * 1000.0, meas[s]) for s in ("MFT", "FFT")}
        out["hofmann_h2s_made_mM"][str(ph)] = float(run.final("H2S"))
    cool = _run(parameters, drift, 100.0, HOFMANN, list(SCHIEBERLE_MFT_UG), 5.0, HOFMANN_BUFFER)
    mft = {t: float(cool.series("MFT")[i + 1]) for i, t in enumerate(SCHIEBERLE_MFT_UG)}
    out["schieberle_mft_100C_6h_over_145C"] = mft[360.0] / hot[5.0]
    out["schieberle_mft_shape_over_30min"] = {str(t): mft[t] / mft[30.0] for t in (60.0, 360.0, 720.0)}
    nes = _run(parameters, drift, 110.0, {"PENT": 250.0, "Cys": 250.0, **OX}, list(NESTLE_MFT_NG_G), 5.5, NESTLE_BUFFER)
    nm = {t: float(nes.series("MFT")[i + 1]) for i, t in enumerate(NESTLE_MFT_NG_G)}
    out["nestle_mft_shape_over_4h"] = {str(t): nm[t] / nm[240.0] for t in (30.0, 60.0, 120.0)}
    zr = _run(parameters, drift, 120.0, {"PENT": 40.0, "Cys": 40.0, **OX}, list(ZHANG_MM["TTCA"]), 6.5, ZHANG_BUFFER)
    # with the route off the engine carries no TTCA or Cys-Amadori at all (r_ttca_form is inactive): None, not a miss
    out["zhang2025b_blank_log10_pred_over_meas"] = {
        s: {str(t): (_log(float(zr.series(s)[i + 1]), m) if float(zr.series(s)[i + 1]) > 1e-9 else None)
            for i, (t, m) in enumerate(ZHANG_MM[s].items())} for s in ZHANG_MM}
    out["zhang2025b_free_cys_fraction_50min"] = float(zr.series("Cys")[-1]) / 40.0
    return out


def arms(reading: str):
    """(label, member path, switch settings) for every arm whose member exists."""
    yield "baseline (switches off)", best_member(reading), {}
    candidates = [(f"B59a-r share {s:g}", f"kinetic_core_b59a_{reading}_share{s:g}_members", {"strecker": s})
                  for s in ps.STRECKER_SULFUR.alternatives]
    candidates.append(("B59b", f"kinetic_core_b59b_{reading}_members", {"route": True}))
    candidates += [(f"B59b + Strecker {s:g}", f"kinetic_core_b59b_{reading}_strecker{s:g}_members",
                    {"route": True, "strecker": s}) for s in ps.STRECKER_SULFUR.alternatives]
    for label, dirname, switches in candidates:
        files = sorted((data_paths.VALIDATION_DIR / dirname).glob("*.json"))
        if files:
            yield label, min(files, key=lambda f: json.loads(f.read_text())["cost"]), switches


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reading", choices=sorted(B64.READINGS), required=True)
    args = parser.parse_args(argv)
    B64.install_reading(args.reading)
    saved = (ps.STRECKER_SULFUR, ps.CYS_AMADORI_ROUTE)
    report = {"artifact": "kinetic_core_b59_holdout_diagnosis", "reading": args.reading,
              "measured": {"hofmann_ppb": HOFMANN_PPB, "schieberle_factor": [7.9, 13.0],
                           "schieberle_mft_ug": SCHIEBERLE_MFT_UG, "nestle_mft_ng_g": NESTLE_MFT_NG_G,
                           "li2011": "not evaluable: no hexose species in the sulfur lane", "zhang2025b_mM": ZHANG_MM, "zhang2025b_free_cys_fraction_50min": 0.024},
              "arms": {}}
    for label, path, switches in arms(args.reading):
        if "strecker" in switches:
            ps.STRECKER_SULFUR = dataclasses.replace(saved[0], enabled=True, h2s_share=switches["strecker"])
        if switches.get("route"):
            ps.CYS_AMADORI_ROUTE = dataclasses.replace(saved[1], enabled=True)
        try:
            member = json.loads(path.read_text())
            if "cost" not in member:  # B60's fit report (the isothermal baseline): its cost on this objective
                r = B64.B60.residual_vector(_vector(path), True)
                member = {"cost": float(0.5 * np.sum(r ** 2)), "status": "B60 shipped optimum", "converged": None}
            report["arms"][label] = {"member": data_paths.rel(path), "cost": member.get("cost"),
                                     "optimiser_status": member.get("status"), "converged": member.get("converged"),
                                     **diagnose(_vector(path), args.reading)}
        finally:
            ps.STRECKER_SULFUR, ps.CYS_AMADORI_ROUTE = saved
        arm = report["arms"][label]
        print(args.reading, label, "cost", round(arm["cost"], 2) if arm["cost"] else None,
              "factor", round(arm["schieberle_mft_100C_6h_over_145C"], 2),
              "nestle", {k: round(v, 3) for k, v in arm["nestle_mft_shape_over_4h"].items()}, flush=True)
    out = data_paths.VALIDATION_DIR / f"kinetic_core_b59_holdout_diagnosis_{args.reading}.json"
    out.write_text(json.dumps(report, indent=2, default=str))
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
