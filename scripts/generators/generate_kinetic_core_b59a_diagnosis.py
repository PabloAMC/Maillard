"""
WAVE B59a recorded diagnosis (prereg sec. 4; not gating): the Hofmann ribose pH 3 / 5 / 7 shape and the
Schieberle 2000 100 C factor, at the B64 baseline's optimum (switch off) and at each B59a optimum (switch on).
Writes results/validation/kinetic_core_b59a_diagnosis.json.

    python scripts/generators/generate_kinetic_core_b59a_diagnosis.py --baseline ramp_plus_hold
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

import generate_kinetic_core_b64_fit as B64  # noqa: E402
from generate_kinetic_core_b59a_fit import best_member  # noqa: E402
from src import data_paths  # noqa: E402
from src.kinetic_core import parameters_sulfur as ps  # noqa: E402
from src.kinetic_core import ph_state, sulfur  # noqa: E402

BUFFER = ph_state.BufferSpec(kind="phosphate", phosphate_mol_l=0.5, declared=True, source="Hofmann 1998 pot")
CHARGE = {"PENT": 100.0, "Cys": 33.0}
MEASURED_PPB = {3.0: {"MFT": 553.0, "FFT": 229.0}, 5.0: {"MFT": 198.0, "FFT": 121.0}, 7.0: {"MFT": 25.0, "FFT": 12.0}}
MW = {"MFT": 114.17, "FFT": 114.17}


def _profile(points):
    return np.array([t for t, _ in points]), np.array([c for _, c in points]) + 273.15


def diagnose(x, profile_points):
    _f, _e, _d, drift = B64.B23.unpack(x)
    parameters = B64.B23.build_parameters(x)
    prof = _profile(profile_points)
    out = {"pH_sweep_log10_pred_over_meas": {}, "h2s_made_mM": {}}
    for ph, meas in MEASURED_PPB.items():
        run = sulfur.integrate_sulfur(parameters, 418.15, CHARGE, [0.0, float(prof[0][-1])], ph=ph, buffer_spec=BUFFER,
                                      ph_drift=drift, rtol=1e-8, atol=1e-14, temperature_profile=prof)
        out["pH_sweep_log10_pred_over_meas"][str(ph)] = {
            s: math.log10(max(run.final(s) * MW[s] * 1000.0, 1e-12) / meas[s]) for s in ("MFT", "FFT")}
        out["h2s_made_mM"][str(ph)] = float(run.final("H2S"))
    hot = sulfur.integrate_sulfur(parameters, 418.15, CHARGE, [0.0, float(prof[0][-1])], ph=5.0, buffer_spec=BUFFER,
                                  ph_drift=drift, rtol=1e-8, atol=1e-14, temperature_profile=prof).final("MFT")
    cool = sulfur.integrate_sulfur(parameters, 373.15, CHARGE, [0.0, 360.0], ph=5.0, buffer_spec=BUFFER,
                                   ph_drift=drift, rtol=1e-8, atol=1e-14).final("MFT")
    out["mft_100C_6h_over_145C"] = cool / hot
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", choices=sorted(B64.READINGS), required=True)
    args = parser.parse_args(argv)
    points = B64.READINGS[args.baseline]["hofmann_pentose_pH5"]
    report = {"artifact": "kinetic_core_b59a_diagnosis", "baseline": args.baseline, "measured_ppb": MEASURED_PPB,
              "measured_factor": 13.0, "arms": {}}
    base = json.loads(best_member(args.baseline).read_text())
    report["arms"]["baseline (switch off)"] = {"cost": base["cost"], **diagnose(np.array(base["x_full"]), points)}
    saved = ps.STRECKER_SULFUR
    for share in saved.alternatives:
        files = sorted((data_paths.VALIDATION_DIR / f"kinetic_core_b59a_{args.baseline}_share{share:g}_members").glob("*.json"))
        if not files:
            continue
        member = min((json.loads(f.read_text()) for f in files), key=lambda m: m["cost"])
        ps.STRECKER_SULFUR = dataclasses.replace(saved, enabled=True, h2s_share=share)
        try:
            report["arms"][f"B59a share {share:g}"] = {"cost": member["cost"], **diagnose(np.array(member["x_full"]), points)}
        finally:
            ps.STRECKER_SULFUR = saved
    out = data_paths.VALIDATION_DIR / f"kinetic_core_b59a_diagnosis_{args.baseline}.json"
    out.write_text(json.dumps(report, indent=2))
    for arm, r in report["arms"].items():
        sweep = {ph: {s: round(v, 2) for s, v in d.items()} for ph, d in r["pH_sweep_log10_pred_over_meas"].items()}
        print(f"{arm:24s} cost {r['cost']:.2f}  log10(pred/meas) {sweep}  H2S pH7 {r['h2s_made_mM']['7.0']:.3g} mM  "
              f"100C/145C {r['mft_100C_6h_over_145C']:.2f} (13)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
