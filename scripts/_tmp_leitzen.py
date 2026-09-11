"""Leitzen 2021 Table 4 (121 C / 18 min, 10 % glucose in water, AMINE-FREE) against the trunk.
Five observables in one pot, and every one of them is a HOLD-OUT level: this VALIDATES, it does not fit."""
import sys, json
from pathlib import Path
from dataclasses import replace
import numpy as np
sys.path.insert(0, "/workspace")
from src.kinetic_core import panel
from src.kinetic_core.engine import TRUNK, ACRYLAMIDE, core_parameters, predict

MW = {"GO": 58.04, "MGO": 72.06, "G": 178.14, "TDG": 162.14, "HMF": 126.11, "DDG": 144.13}
# Table 4, scheme A, 121 C / 18 min, ug/mL of solution (= mg/L)
MEASURED_UG_PER_ML = {"GO": 5.6, "MGO": 2.6, "G": 7.5, "TDG": 52.2, "HMF": 17.4, "DDG": 55.5}
SD = {"GO": 1.3, "MGO": 0.2, "G": 1.4, "TDG": 4.0, "HMF": 3.9, "DDG": 1.7}
LABEL = {"GO": "glyoxal", "MGO": "methylglyoxal", "G": "glucosone", "TDG": "3-deoxyglucosone", "HMF": "5-HMF", "DDG": "3,4-dideoxyglucosone"}
DIRECT = ("k_glc_tdg", "k_fru_odg", "k_fru_int")

bench = json.loads(Path("data/benchmarks/external_validation/maillard_path/"
                        "mp_holdout_glucose_only_autoclave_121C_Steinhagen2021.json").read_text())
print("pot:", bench["precursors"], "|", bench["conditions"]["temp_C"], "C,", bench["conditions"]["time_min"], "min, pH", bench["conditions"]["ph"])
spec = panel.core_spec(bench)

def scaled(params, mult):
    out = dict(params)
    for k in DIRECT:
        if k in out and out[k].k_ref:
            out[k] = replace(out[k], k_ref=out[k].k_ref * mult)
    return out

base = core_parameters(TRUNK)
for mult in (1.0,):
    run = predict(spec, [], parameters=scaled(base, mult))
    st = run.run_metadata.get("final_state") or {}
    if not st:
        # fall back: integrate the trunk directly
        from src.kinetic_core.integrate import integrate
        from src.kinetic_core import trunk_conditions
        from types import SimpleNamespace
        p, _ = trunk_conditions.apply(scaled(base, mult), SimpleNamespace(ph=bench["conditions"]["ph"], water_activity=None, atmosphere=None))
        r = integrate(p, bench["conditions"]["temp_C"] + 273.15,
                      {"Glc": float(bench["precursors"]["D-Glucose"]["concentration_mM"])},
                      np.array([0.0, float(bench["conditions"]["time_min"])]), rtol=1e-8, atol=1e-16)
        st = {k: float(r.series(k)[-1]) for k in MW}
    print(f"\n--- amine-free entries x{mult:g} ---")
    print(f"{'species':20s} {'measured ug/mL':>15} {'model ug/mL':>13} {'fold':>8}")
    for k in ("TDG", "DDG", "HMF", "MGO", "G", "GO"):
        model = float(st.get(k, 0.0)) * MW[k]   # mmol/L * g/mol = mg/L = ug/mL
        m = MEASURED_UG_PER_ML[k]
        fold = max(model / m, m / model) if model > 0 else float("inf")
        print(f"{LABEL[k]:20s} {m:10.1f} +/-{SD[k]:<4.1f} {model:13.4g} {fold:8.2f}")
