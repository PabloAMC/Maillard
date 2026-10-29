"""
WAVE B64, amendment 1: the measurement that bears directly on the thermal reading (2026-10-28).

Schieberle 2000 (Table IV and p. 141) runs the reference pot of Hofmann 1998 (ribose 100 mM + cysteine 33 mM,
100 mL, 0.5 M phosphate pH 5.0, the same laboratory autoclave) at 100 C for 30-720 min, and reports MFT "by a
factor of 13 higher" after 6 h at 100 C than after the "20 min at 145 C" pot. A ratio from one laboratory and
one pot cancels its calibration, and it is not a row of B60's objective. Under each reading's own optimum this
generator predicts that factor and the 100 C time shape (MFT and FFT at 60 / 360 / 720 min over 30 min), with
the 145 C pot integrated along that reading's profile. The 100 C pot is run isothermal and, as a sensitivity,
after a 20 -> 100 C come-up at the 145 C ramp's own rate (6.25 K/min, 12.8 min); neither is fitted.
Writes results/validation/kinetic_core_b64_discriminator.json.

    python scripts/generators/generate_kinetic_core_b64_discriminator.py
"""
from __future__ import annotations

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
from src import data_paths  # noqa: E402
from src.kinetic_core import ph_state, sulfur  # noqa: E402

OUT = data_paths.VALIDATION_DIR / "kinetic_core_b64_discriminator.json"
CHARGE = {"PENT": 100.0, "Cys": 33.0}
BUFFER = ph_state.BufferSpec(kind="phosphate", phosphate_mol_l=0.5, declared=True, source="Hofmann 1998 pot")
MEASURED = {  # Schieberle 2000 Table IV, ug per 100 mL pot, and p. 141
    "MFT": {30.0: 4.5, 60.0: 13.8, 360.0: 156.0, 720.0: 179.0},
    "FFT": {30.0: 2.0, 60.0: 3.1, 360.0: 110.0, 720.0: 132.0},
    "MFT_100C_6h_over_145C": 13.0,
}
PROFILES_145 = {
    "ramp_only": ((0.0, 20.0), (20.0, 145.0)),
    "ramp_plus_hold": ((0.0, 20.0), (20.0, 145.0), (40.0, 145.0)),
    "isothermal": None,
    **{name: B64.READINGS[name]["hofmann_pentose_pH5"] for name in B64.READINGS if "newtonian" in name},
}


def _members():
    out = {}
    for reading in sorted(B64.READINGS):
        files = sorted((data_paths.VALIDATION_DIR / f"kinetic_core_b64_{reading}_members").glob("*.json"))
        best = min((json.loads(f.read_text()) for f in files), key=lambda m: m["cost"], default=None)
        if best is not None:
            out[reading] = best
    return out


def _run(parameters, drift, temperature_c, minutes, profile=None, grid=None):
    times = np.array(grid if grid is not None else [0.0, minutes], dtype=float)
    prof = None
    if profile is not None:
        prof = (np.array([t for t, _ in profile]), np.array([c for _, c in profile]) + 273.15)
    return sulfur.integrate_sulfur(parameters, temperature_c + 273.15, CHARGE, times, ph=5.0, buffer_spec=BUFFER,
                                   ph_drift=drift, rtol=1e-8, atol=1e-14, temperature_profile=prof)


def evaluate(x: np.ndarray) -> dict:
    _f, _e, _d, drift = B64.B23.unpack(x)
    parameters = B64.B23.build_parameters(x)
    out = {}
    times = [30.0, 60.0, 360.0, 720.0]
    for label, come_up in (("isothermal_100C", None), ("come_up_100C", ((0.0, 20.0), (12.8, 100.0)))):
        if come_up is None:
            run = _run(parameters, drift, 100.0, 720.0, grid=[0.0] + times)
            series = {s: {t: float(run.series(s)[i + 1]) for i, t in enumerate(times)} for s in ("MFT", "FFT")}
        else:
            series = {s: {} for s in ("MFT", "FFT")}
            for t in times:
                profile = come_up + ((12.8 + t, 100.0),)
                run = _run(parameters, drift, 100.0, 12.8 + t, profile=profile)
                for s in ("MFT", "FFT"):
                    series[s][t] = run.final(s)
        out[label] = series
    hot = {}
    for reading, profile in PROFILES_145.items():
        minutes = 20.0 if profile is None else profile[-1][0]
        hot[reading] = _run(parameters, drift, 145.0, minutes, profile=profile).final("MFT")
    out["mft_145C_by_reading"] = hot
    return out


def main() -> int:
    members = _members()
    report = {"artifact": "kinetic_core_b64_discriminator", "prereg": data_paths.rel(B64.PREREG),
              "measured": {k: v for k, v in MEASURED.items()}, "by_optimum": {}}
    for reading, member in members.items():
        res = evaluate(np.array(member["x_full"], dtype=float))
        own = res["mft_145C_by_reading"][reading]
        factors = {kind: res[kind]["MFT"][360.0] / own for kind in ("isothermal_100C", "come_up_100C")}
        shape = {kind: {s: {str(t): res[kind][s][t] / res[kind][s][30.0] for t in (60.0, 360.0, 720.0)}
                        for s in ("MFT", "FFT")} for kind in ("isothermal_100C", "come_up_100C")}
        report["by_optimum"][reading] = {
            "member_cost": member["cost"],
            "mft_100C_6h_over_145C": factors,
            "log10_factor_vs_measured": {k: math.log10(v / MEASURED["MFT_100C_6h_over_145C"]) for k, v in factors.items()},
            "shape_over_30min": shape,
            "raw_mM": res,
        }
        print(reading, "cost", round(member["cost"], 2), "factor", {k: round(v, 1) for k, v in factors.items()},
              "(measured 13)", flush=True)
    OUT.write_text(json.dumps(report, indent=2, default=str))
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
