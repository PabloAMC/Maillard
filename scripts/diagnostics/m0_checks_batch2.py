#!/usr/bin/env python
"""
m0_checks_batch2.py -- the numbers of the misfit inventory's M0 checks batch 2 (2026-11-06): CK-23, CK-34, CK-42,
CK-42b and CK-56. Read-only diagnostics: the shipped (B60) engine, scripted parameter variants in memory, arithmetic on
printed values. Nothing is fitted, nothing is written under data/ or src/. The pass rules are in the record
(results/validation/m0_checks_batch2.md, committed before this script existed).

Each engine part runs in its own subprocess (the B60 objective's import sets its quick-setting integrator patch; the
Schieberle and Munch pots run on the plain engine path, as M1.0 ran them):

    python scripts/diagnostics/m0_checks_batch2.py            # all parts, writes the diagnostics JSON
    python scripts/diagnostics/m0_checks_batch2.py --part ck23_schieberle --out /tmp/x.json   # one part

Not computed under any variant (open cases' decisive pots): Hofmann 1998 NF + Cys (`fed_nf_cys`) and the 145 C ribose
+ Cys pot (M2 Rt-H), Cerny 2003 pot D (M2), the Zhai 2023 JAFC 140 C pots (M2 T-140), Whitfield 1999 (M1 D6), the Zhou
2023 / Zhang 2024 pots (M1 D3 - D5; M1b guard), the US 10,076,126 pot (B59b hold-out), van Seeventer 2001 Table 1.
Output: results/validation/m0_checks_batch2_diagnostics.json.
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import math
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence

ROOT = Path(__file__).resolve().parents[2]
for _p in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

R_GAS = 8.314462618e-3  # kJ/(mol K)
CELSIUS = 273.15
T_REF_K = 418.15        # B2.3's reference temperature for k_ref (145 C)
OUT = ROOT / "results" / "validation" / "m0_checks_batch2_diagnostics.json"

# ---------------------------------------------------------------------------------------------------------------
# shared helpers
# ---------------------------------------------------------------------------------------------------------------


def arrhenius(k_ref: float, ea_kj: float, t_k: float, t_ref_k: float = T_REF_K) -> float:
    return float(k_ref) * math.exp(-float(ea_kj) / R_GAS * (1.0 / t_k - 1.0 / t_ref_k))


def scaled(params: Mapping[str, Any], factors: Mapping[str, float], ea: Optional[Mapping[str, float]] = None
           ) -> Dict[str, Any]:
    """A copy of the operative parameter dict with k_ref multiplied (and optionally Ea replaced)."""
    out = dict(params)
    for key, f in factors.items():
        p = out[key]
        changes = {"k_ref": float(p.k_ref) * float(f)}
        if ea and key in ea:
            changes["ea_kj_mol"] = float(ea[key])
        out[key] = dataclasses.replace(p, **changes)
    for key, e in (ea or {}).items():
        if key not in factors:
            out[key] = dataclasses.replace(out[key], ea_kj_mol=float(e))
    return out


def whitened(pred: float, row: Mapping[str, Any]) -> float:
    """B2.3's residual (residuals()) for one row."""
    floor = 1e-12
    t, sigma, kind = float(row["target"]), float(row["sigma_log"]), row["kind"]
    if not math.isfinite(pred):
        return 25.0
    if kind == "ph_endpoint":
        r = (pred - t) / sigma
    elif kind in ("ceiling", "peak_fraction_ceiling"):
        r = max(0.0, math.log10((pred + floor) / (t + floor))) / sigma
    elif kind == "floor":
        r = min(0.0, math.log10((pred + floor) / (t + floor))) / sigma
    else:
        r = math.log10((pred + floor) / (t + floor)) / sigma
    return max(-25.0, min(25.0, r))


# ---------------------------------------------------------------------------------------------------------------
# B60 in-sample guard (subset of B60's 52 rows; quick setting, as B60's residual_by_row)
# ---------------------------------------------------------------------------------------------------------------

#: systems whose rows serve as in-sample guards. Excluded: fed_nf_cys and hofmann_pentose_pH5 (M2 Rt-H), whitfield_*
#: (M1 D6), zhou_* and zhang_fig1_* (M1 D3 - D5; M1b guard), kumazawa_* (v1 TEST source), cerny_ternary (needs the
#: flux budget; not an MFT-sink row), hofmann_glucose / fructose (no active row).
GUARD_SYSTEMS_CK23 = ("fed_ribose_h2s", "fed_tdp_h2s", "fed_furfural_h2s", "fed_nf_h2s", "fed_c2c3", "fed_c2c3_pH3",
                      "fed_c2c3_pH7", "fed_thiamine", "fed_mgo_h2s_1to1", "fed_mgo_h2s_1to2", "kang_ttca_100",
                      "kang_ttca_120", "kang_freecys_100", "kang_freecys_120", "vanseeventer_130C", "feng_arp_100",
                      "feng_arp_120")
GUARD_SYSTEMS_CK34 = ("fed_ribose_h2s", "fed_tdp_h2s", "fed_furfural_h2s", "fed_nf_h2s", "fed_mgo_h2s_1to1",
                      "fed_mgo_h2s_1to2")


def _b60():
    import generate_kinetic_core_b60_fit as B60  # noqa: F401  (B47's rows; B51's pH state; B60's quick setting)
    from generate_kinetic_core_b8_laplace import frozen_vector
    B23 = B60.B23
    report = json.loads((ROOT / "results" / "validation" / "kinetic_core_b60_fit_report.json").read_text())
    x = frozen_vector(report)
    params = B23.build_parameters(x)
    drift = B23.unpack(x)[3]
    b60_resid = report["members"][report["best_start"]]["residual_by_row"]
    return B60, B23, params, drift, b60_resid


def guard_rows(B23, params, drift, systems: Sequence[str]) -> Dict[str, Dict[str, float]]:
    saved_sys, saved_rows = B23.SYSTEMS, B23.ACTIVE_FIT_ROWS
    keep = set(systems)
    rows = [r for r in saved_rows if r["system"] in keep and r.get("system_b", r["system"]) in keep]
    try:
        B23.SYSTEMS = {k: saved_sys[k] for k in systems}
        B23.ACTIVE_FIT_ROWS = rows
        pred, _ = B23.observables(params, True, drift)
    finally:
        B23.SYSTEMS, B23.ACTIVE_FIT_ROWS = saved_sys, saved_rows
    return {r["id"]: {"pred": float(pred[r["id"]]), "target": float(r["target"]), "sigma_log": float(r["sigma_log"]),
                      "r": whitened(float(pred[r["id"]]), r)} for r in rows}


def moves_out(base: Mapping[str, Mapping[str, float]], var: Mapping[str, Mapping[str, float]]) -> List[str]:
    return [k for k in var if abs(var[k]["r"]) > 2.0 and abs(var[k]["r"]) - abs(base[k]["r"]) > 1.0]


# ---------------------------------------------------------------------------------------------------------------
# CK-23
# ---------------------------------------------------------------------------------------------------------------

F_MIN_BAND = (1e-6, 1e-2)       # declared: the co-heated standard's surviving fraction a SIM quantitation tolerates
PATENT_T_C, PATENT_MIN, PATENT_PH = 110.0, 240.0, 5.5
GIGL_EA = {"min": 7.2, "central": 60.2, "max": 102.4}   # Gigl 2021 (gigl2021 dossier sec. 6d), kJ/mol
GROSCH_LAMBDA_100C = {"thiamine": (0.044, 0.035, 0.053), "thiamine_cys": (0.034, 0.026, 0.041)}  # dossier sec. 3 [Z]


def lambda_max(f_min: float, minutes: float = PATENT_MIN) -> float:
    return math.log(1.0 / f_min) / minutes


def standard_survival(params, drift, t_c: float, ph: float, buffer, minutes: float) -> Dict[str, float]:
    """The engine's first-order loss of a trace MFT / FFT 'standard' alone in buffer (no precursors)."""
    import m1_preflight_checks as M10
    c0 = 1e-3
    run, _check, _meta = M10.dense_run(params, drift, {"MFT": c0, "FFT": c0}, t_c, minutes, ph, buffer,
                                       check_species=("MFT", "FFT"))
    out = {}
    for s in ("MFT", "FFT"):
        surv = float(run.final(s)) / c0
        out[f"{s}_survival"] = surv
        out[f"{s}_lambda_eff_per_min"] = -math.log(max(surv, 1e-300)) / minutes
    out["ph_in_situ"] = float(run.metadata.get("ph_final_in_situ", float("nan")))
    return out


def part_ck23_schieberle() -> Dict[str, Any]:
    import m1_preflight_checks as M10
    from src.kinetic_core.engine import SULFUR, core_parameters, core_ph_drift
    from src.kinetic_core.ph_state import BufferSpec
    params, drift = core_parameters(SULFUR), core_ph_drift()
    km, kf = params["k_mft_decay"], params["k_fft_decay"]
    patent_buffer = BufferSpec(kind="phosphate", phosphate_mol_l=0.2, declared=True,
                               source="US 10,076,126: 'an aqueous phosphate buffer at pH 5.5', strength not printed; "
                                      "0.2 M declared (CF-013 analysis)")
    grosch_buffer = BufferSpec(kind="phosphate", phosphate_mol_l=0.2, declared=True,
                               source="stand-in for Grosch 1992's 0.2 M pyrophosphate pH 5.7 (thiolate share only)")
    out: Dict[str, Any] = {"shipped": {"k_mft_decay_k_ref_145C": km.k_ref, "k_fft_decay_k_ref_145C": kf.k_ref,
                                       "Ea_decay_thiol_sink": km.ea_kj_mol}}
    t110 = PATENT_T_C + CELSIUS
    out["arithmetic_k_mft_decay"] = {f"{t}C": arrhenius(km.k_ref, km.ea_kj_mol, t + CELSIUS) for t in (100, 110, 145)}
    out["arithmetic_k_fft_decay"] = {f"{t}C": arrhenius(kf.k_ref, kf.ea_kj_mol, t + CELSIUS) for t in (100, 110, 145)}
    std110 = standard_survival(params, drift, PATENT_T_C, PATENT_PH, patent_buffer, PATENT_MIN)
    std100 = standard_survival(params, drift, 100.0, 5.7, grosch_buffer, 120.0)
    out["standard_only_110C_pH5_5_240min"] = std110
    out["standard_only_100C_pH5_7_120min"] = std100
    bound = {}
    for f in F_MIN_BAND:
        lm = lambda_max(f)
        bound[f"f_min_{f:g}"] = {
            "lambda_max_110C": lm,
            "engine_MFT_over_bound": std110["MFT_lambda_eff_per_min"] / lm,
            "engine_FFT_over_bound": std110["FFT_lambda_eff_per_min"] / lm,
            "lambda_max_100C_by_Ea": {k: lm * math.exp(-ea / R_GAS * (1 / 373.15 - 1 / t110)) for k, ea in GIGL_EA.items()},
            "lambda_max_100C_at_shipped_Ea": lm * math.exp(-km.ea_kj_mol / R_GAS * (1 / 373.15 - 1 / t110)),
        }
    out["bound"] = bound
    # Ea that would let Grosch's 100 C lambda (thiamine + Cys central) sit under the bound at 110 C
    lam100 = GROSCH_LAMBDA_100C["thiamine_cys"][0]
    out["grosch_crosscheck"] = {
        "grosch_lambda_100C": GROSCH_LAMBDA_100C,
        "engine_lambda_eff_100C": std100["MFT_lambda_eff_per_min"],
        "max_Ea_for_grosch_central_under_bound": {
            f"f_min_{f:g}": (R_GAS * math.log(lambda_max(f) / lam100) / (1 / 373.15 - 1 / t110)
                             if lambda_max(f) > lam100 else None) for f in F_MIN_BAND},
    }
    # forward Schieberle T4 (M1.0's runner): baseline, the two bound corners (Ea held), k_mft_decay = 0, and the
    # Gigl-central Ea with k(110 C) at the loose corner
    other = std110["MFT_lambda_eff_per_min"] - arrhenius(km.k_ref, km.ea_kj_mol, t110)
    k110 = arrhenius(km.k_ref, km.ea_kj_mol, t110)
    variants = {"shipped": ({}, None)}
    for f in F_MIN_BAND:
        target = lambda_max(f) - other
        variants[f"bound_f_min_{f:g}"] = ({"k_mft_decay": max(target, 0.0) / k110}, None)
    variants["k_mft_decay_zero"] = ({"k_mft_decay": 0.0}, None)
    # Ea at Gigl's central value, k_mft_decay(110 C) at the loose bound: k_ref(145) = k110_target * exp(...)
    target = lambda_max(F_MIN_BAND[0]) - other
    kref_gigl = target / math.exp(-GIGL_EA["central"] / R_GAS * (1 / t110 - 1 / T_REF_K))
    variants["bound_f_min_1e-06_Ea_gigl_central"] = ({"k_mft_decay": kref_gigl / km.k_ref},
                                                     {"k_mft_decay": GIGL_EA["central"]})
    out["other_MFT_loss_110C_per_min"] = other
    out["schieberle"] = {}
    for name, (fac, ea) in variants.items():
        p = scaled(params, fac, ea) if (fac or ea) else params
        d1 = M10.d1(p, drift)
        mft = d1["species"]["MFT"]
        fft = d1["species"]["FFT"]
        out["schieberle"][name] = {
            "factor_on_k_mft_decay_k_ref": fac.get("k_mft_decay", 1.0),
            "Ea_k_mft_decay": (ea or {}).get("k_mft_decay", km.ea_kj_mol),
            "k_mft_decay_100C": arrhenius(p["k_mft_decay"].k_ref, p["k_mft_decay"].ea_kj_mol, 373.15),
            "k_mft_decay_145C": arrhenius(p["k_mft_decay"].k_ref, p["k_mft_decay"].ea_kj_mol, 418.15),
            "MFT_ug_l": mft["R0_ug_l"], "MFT_peak_min": mft["R0_peak_min"],
            "MFT_360_over_30_resid_dex": mft["ratios"]["360_over_30"]["residual_dex"],
            "MFT_720_over_30_resid_dex": mft["ratios"]["720_over_30"]["residual_dex"],
            "MFT_measured_dex": {k: v["measured_dex"] for k, v in mft["ratios"].items()},
            "FFT_360_over_30_resid_dex": fft["ratios"]["360_over_30"]["residual_dex"],
            "FFT_720_over_30_resid_dex": fft["ratios"]["720_over_30"]["residual_dex"],
            "_params": {"factor": fac.get("k_mft_decay", 1.0), "ea": (ea or {}).get("k_mft_decay")},
        }
    return out


def part_ck23_guard(schieberle_part: Mapping[str, Any]) -> Dict[str, Any]:
    _B60, B23, params, drift, b60_resid = _b60()
    base = guard_rows(B23, params, drift, GUARD_SYSTEMS_CK23)
    repro = {k: base[k]["r"] - b60_resid[k] for k in base}
    out = {"baseline": base, "reproduces_B60_max_abs_dr": max(abs(v) for v in repro.values()), "variants": {}}
    for name, rec in schieberle_part["schieberle"].items():
        if name == "shipped":
            continue
        fac, ea = rec["_params"]["factor"], rec["_params"]["ea"]
        p = scaled(params, {"k_mft_decay": fac}, {"k_mft_decay": ea} if ea else None)
        var = guard_rows(B23, p, drift, GUARD_SYSTEMS_CK23)
        out["variants"][name] = {"rows": {k: {"r": var[k]["r"], "dr": var[k]["r"] - base[k]["r"]} for k in var},
                                 "moved_out": moves_out(base, var)}
    return out


# ---------------------------------------------------------------------------------------------------------------
# CK-34
# ---------------------------------------------------------------------------------------------------------------

KAW_25C, DH_SOL_KJ = 0.403, -17.5   # henry_constants.yml, Hydrogen Sulfide (Sander 5.0.0)
H2S_REACTANT_KEYS = ("k_fur_fft", "k_fur_fft_hs", "k_nf_mft", "k_nf_mp3p", "k_nf_shet", "k_ddp_mft", "k_ddp_mft_hs",
                     "k_mgo_mp", "k_dmhf_h2s", "k_redox_mp3p", "k_h2s_loss")
MUNCH = {  # Munch 1997 Table 5 (furfural 1 mmol; FFT ug per pot); expt 3 = Hofmann 1998 T3 (CF-101)
    "expt1_cys_1mmol": 25.3, "expt2_h2s_0.25mmol": 130.1, "expt3_h2s_1mmol": 550.8}


def kaw(t_k: float) -> float:
    b = -DH_SOL_KJ / R_GAS  # K
    return KAW_25C * (298.15 / t_k) * math.exp(-b * (1.0 / t_k - 1.0 / 298.15))


def phi_dissolved(t_k: float, ph_in_situ: float, vg_over_vl: float) -> Dict[str, float]:
    from src.kinetic_core.parameters_sulfur import h2s_pk1_suleimenov_seward
    pk1 = h2s_pk1_suleimenov_seward(t_k)
    alpha0 = 1.0 / (1.0 + 10.0 ** (ph_in_situ - pk1))
    k = kaw(t_k)
    return {"Kaw": k, "pK1": pk1, "alpha0_neutral": alpha0, "Vg_over_Vl": vg_over_vl,
            "phi": 1.0 / (1.0 + k * alpha0 * vg_over_vl)}


def _munch_runs(params, drift, scale_l: float) -> Dict[str, float]:
    """Munch T5's three furfural pots, FFT final (mmol/L); scale_l 1 = 50 mL (20 mM), 0.5 = 100 mL (10 mM)."""
    import numpy as np
    from generate_kinetic_core_b2_3_fit import BUFFER_HOFMANN
    from src.kinetic_core.sulfur import integrate_sulfur
    c = 20.0 * scale_l
    pots = {"expt1_cys_1mmol": {"FUR": c, "Cys": c}, "expt2_h2s_0.25mmol": {"FUR": c, "H2S": c / 4.0},
            "expt3_h2s_1mmol": {"FUR": c, "H2S": c}}
    out = {}
    for name, init in pots.items():
        run = integrate_sulfur(params, 145.0 + CELSIUS, dict(init), np.array([0.0, 20.0]), ph=5.0,
                               buffer_spec=BUFFER_HOFMANN, ph_drift=drift, rtol=1e-8, atol=1e-14)
        out[name] = float(run.final("FFT"))
        out[name + "_ph_in_situ"] = float(run.metadata["ph_final_in_situ"])
        out[name + "_FUR_left_frac"] = float(run.final("FUR")) / c
    out["contrast_cys_over_h2s"] = out["expt1_cys_1mmol"] / out["expt3_h2s_1mmol"]
    out["dose_ratio_1_over_0.25"] = out["expt3_h2s_1mmol"] / out["expt2_h2s_0.25mmol"]
    return out


def part_ck34() -> Dict[str, Any]:
    _B60, B23, params, drift, b60_resid = _b60()
    import numpy as np
    from src.kinetic_core.sulfur import integrate_sulfur
    missing = [k for k in H2S_REACTANT_KEYS if k not in params]
    keys = [k for k in H2S_REACTANT_KEYS if k in params]
    spec = B23.SYSTEMS["fed_furfural_h2s"]
    run = integrate_sulfur(params, spec["t_c"] + CELSIUS, spec["initial"], np.array([0.0, 20.0]), ph=5.0,
                           buffer_spec=spec["buffer"], ph_drift=drift, rtol=1e-8, atol=1e-14)
    ph_is = float(run.metadata["ph_final_in_situ"])
    t_k = 145.0 + CELSIUS
    meas_contrast = MUNCH["expt1_cys_1mmol"] / MUNCH["expt3_h2s_1mmol"]
    meas_dose = MUNCH["expt3_h2s_1mmol"] / MUNCH["expt2_h2s_0.25mmol"]
    out: Dict[str, Any] = {"missing_keys": missing, "scaled_keys": keys, "ph_in_situ_fed_furfural": ph_is,
                           "measured": {"contrast": meas_contrast, "dose_ratio": meas_dose}}
    base_guard = guard_rows(B23, params, drift, GUARD_SYSTEMS_CK34)
    out["guard_baseline"] = base_guard
    out["guard_reproduces_B60_max_abs_dr"] = max(abs(base_guard[k]["r"] - b60_resid[k]) for k in base_guard)
    readings = {"50mL_VgVl3": (1.0, 3.0), "100mL_VgVl1": (0.5, 1.0)}
    out["munch"] = {}
    for rname, (scale_l, vg) in readings.items():
        ph_op = phi_dissolved(t_k, ph_is, vg)
        base = _munch_runs(params, drift, scale_l)
        arms = {"no_partition": base}
        for loss_scaled in (True, False):
            fac = {k: ph_op["phi"] for k in keys if loss_scaled or k != "k_h2s_loss"}
            arms["partition" + ("" if loss_scaled else "_h2s_loss_unscaled")] = _munch_runs(scaled(params, fac), drift,
                                                                                           scale_l)
        for a in arms.values():
            a["contrast_resid_dex"] = math.log10(a["contrast_cys_over_h2s"] / meas_contrast)
            a["dose_resid_dex"] = math.log10(a["dose_ratio_1_over_0.25"] / meas_dose)
        out["munch"][rname] = {"operator": ph_op, "arms": arms}
    # in-sample guard and the bias (Hofmann pots: 50 mL in the 200 mL autoclave, Vg/Vl 3)
    op3 = phi_dissolved(t_k, ph_is, 3.0)
    p_part = scaled(params, {k: op3["phi"] for k in keys})
    var_guard = guard_rows(B23, p_part, drift, GUARD_SYSTEMS_CK34)
    out["guard_partition"] = {"rows": {k: {"r": var_guard[k]["r"], "dr": var_guard[k]["r"] - base_guard[k]["r"],
                                           "pred": var_guard[k]["pred"]} for k in var_guard},
                              "moved_out": moves_out(base_guard, var_guard)}

    def bisect_restore(row_id: str, system: str, extra_keys: Sequence[str]) -> Dict[str, Any]:
        # the response is not monotone in the constant (fast formation then the 145 C sink empties the pot), so
        # scan up from 0 dex for the FIRST crossing, then bisect inside that step
        target = base_guard[row_id]["pred"]

        def pred_at(dex: float) -> float:
            fac = {k: op3["phi"] for k in keys}
            for k in extra_keys:
                fac[k] = fac.get(k, 1.0) * 10.0 ** dex
            return guard_rows(B23, scaled(params, fac), drift, (system,))[row_id]["pred"]

        scan = []
        lo = hi = None
        prev = 0.0
        for i in range(0, 41):
            dex = 0.05 * i
            p = pred_at(dex)
            scan.append((dex, p))
            if p >= target:
                lo, hi = prev, dex
                break
            prev = dex
        if lo is None:
            return {"restoring_factor_dex": None, "target_pred_B60": target, "scan": scan,
                    "note": "no crossing within 0 - 2 dex: the row cannot be restored by this constant alone"}
        for _ in range(12):
            mid = 0.5 * (lo + hi)
            if pred_at(mid) < target:
                lo = mid
            else:
                hi = mid
        return {"restoring_factor_dex": 0.5 * (lo + hi), "target_pred_B60": target,
                "pred_at_restore": pred_at(0.5 * (lo + hi)), "peak_scan_max": max(p for _d, p in scan)}

    out["bias"] = {
        "k_fur_fft_and_hs_for_fed_furfural_h2s_FFT": bisect_restore("fed_furfural_h2s_FFT", "fed_furfural_h2s",
                                                                    ("k_fur_fft", "k_fur_fft_hs")),
        "k_nf_mft_for_fed_nf_h2s_MFT": bisect_restore("fed_nf_h2s_MFT", "fed_nf_h2s", ("k_nf_mft",)),
    }
    # the Munch contrast at the bias-restored furfural constants (50 mL reading)
    b = out["bias"]["k_fur_fft_and_hs_for_fed_furfural_h2s_FFT"]["restoring_factor_dex"] or 0.0
    fac = {k: op3["phi"] for k in keys}
    for k in ("k_fur_fft", "k_fur_fft_hs"):
        fac[k] *= 10.0 ** b
    restored = _munch_runs(scaled(params, fac), drift, 1.0)
    restored["contrast_resid_dex"] = math.log10(restored["contrast_cys_over_h2s"] / meas_contrast)
    restored["dose_resid_dex"] = math.log10(restored["dose_ratio_1_over_0.25"] / meas_dose)
    out["munch_partition_with_restored_furfural_constants_50mL"] = restored
    return out


# ---------------------------------------------------------------------------------------------------------------
# CK-42 / CK-42b (engine-free; only k_cys_h2s read from the shipped parameters)
# ---------------------------------------------------------------------------------------------------------------

GROSCH = {  # Grosch & Zeiler-Hilgart 1992 (grosch1992 dossier sec. 3), ug MFT in 500 mL, Method I
    "T2": {"no1_thiamine": (4.7, 1.2, 3), "no3_thi_rib": (5.2, 0.4, 2), "no4_thi_cys": (33.0, 3.7, 4),
           "no6_thi_cys_rib": (19.0, 2.5, 3), "no7_thi_gsh_rib": (23.0, 1.2, 3), "no8_thi_gsh_cys_rib": (9.0, 2.8, 4)},
    "T3_cys_umol": (0, 85, 170, 255, 325, 340, 680, 1360, 2720),
    "T3_thi30": (5.2, 14, 18, 18, 19, 14, 10, 4.8, 2.0),
    "T3_thi6": {170: 2.1, 340: 2.8, 680: 1.6, 1360: 0.5},
    "T6_h2s": {"thiamine_plus_h2s": (15.0, 7.0, 4), "h2s_charge_mmol": 26.8},
    "volume_L": 0.5, "thiamine_umol": 30.0, "gsh_umol": 230.0, "cys_T2_umol": 325.0,
}
ZHENG_HO_EA_BAND = (123.0, 135.0)   # Zheng & Ho 1994, the Cys -> H2S class (M2 prereg sec. 2)
PYROPHOSPHATE = {"conc_M": 0.2, "pKa3": 6.6, "pH": 5.7}  # pKa3 of H4P2O7 ~6.6 (declared textbook value)


def part_ck42() -> Dict[str, Any]:
    from src.kinetic_core.engine import SULFUR, core_parameters
    p = core_parameters(SULFUR)["k_cys_h2s"]
    out: Dict[str, Any] = {"k_cys_h2s": {"k_ref_145C": p.k_ref, "Ea": p.ea_kj_mol, "source": p.source_anchor}}
    vol = GROSCH["volume_L"]
    cys_mm = GROSCH["cys_T2_umol"] / 1000.0 / vol          # 0.65 mM
    rows = {}
    for ea in (ZHENG_HO_EA_BAND[0], p.ea_kj_mol, ZHENG_HO_EA_BAND[1]):
        # anchored at the shipped k at 145 C; the barrier only moves the 100 C extrapolation
        k100 = arrhenius(p.k_ref, ea, 373.15)
        frac = 1.0 - math.exp(-k100 * 120.0)
        rows[f"Ea_{ea:g}"] = {"k_cys_h2s_100C_per_min": k100, "fraction_converted_120min": frac,
                              "h2s_released_uM_at_0.65mM_Cys": 1000.0 * cys_mm * frac}
    out["cys_derived_sulfide"] = rows
    h2s_charge_mm = GROSCH["T6_h2s"]["h2s_charge_mmol"] / vol   # 53.6 mM if all dissolved
    retained = (1e-3, 1e-2, 1e-1, 1.0)
    worst = max(r["h2s_released_uM_at_0.65mM_Cys"] for r in rows.values())
    out["table_vi_arm"] = {
        "h2s_charge_mM_if_all_dissolved": h2s_charge_mm,
        "dissolved_uM_by_retained_fraction": {f"{f:g}": 1000.0 * h2s_charge_mm * f for f in retained},
        "rate_factor_needed_by_retained_fraction": {f"{f:g}": 1000.0 * h2s_charge_mm * f / worst for f in retained},
        "response_measured": "thiamine + H2S 15 +- 7 vs thiamine 4.7 +- 1.2 ug (x3.2); thiamine + Cys 33 +- 3.7 (x7.0)",
    }
    # falling limb: MFT ~ F / (lambda0 + k_c c^n) with lambda0 = Grosch's Cys-free 0.044 /min (Table VIII)
    lam0 = GROSCH_LAMBDA_100C["thiamine"][0]
    c325, c2720 = 325 / 1000 / vol, 2720 / 1000 / vol
    m325, m2720 = 19.0, 2.0
    fall = m325 / m2720
    upper = GROSCH_LAMBDA_100C["thiamine_cys"][2]

    def capture_law(n: int, lam0_: float, fall_: float) -> Dict[str, Any]:
        # lam(c) = lam0 + k c^n ; lam(c2720) / lam(c325) = fall (constant formation F over the 2 h)
        den = c2720 ** n - fall_ * c325 ** n
        if den <= 0.0:   # the law cannot make this fall at any k >= 0 (its ceiling is (c2720 / c325)^n)
            return {"feasible": False, "ceiling_fall": (c2720 / c325) ** n, "lambda_at_0.65mM": None,
                    "refuted": True}
        k = lam0_ * (fall_ - 1.0) / den
        lam065 = lam0_ + k * c325 ** n
        return {"feasible": True, "k_c": k, "lambda_at_0.65mM": lam065, "lambda_at_5.44mM": lam0_ + k * c2720 ** n,
                "refuted": lam065 > upper}

    # the rule's central reading, then the noise bands: lam0 over Table VIII's Cys-free band; the fall over 2 SD of
    # its two Table III values (Table II's relative SDs transferred: 13 % at 19 ug, 26 % at the 4.7 ug level)
    sd_fall_dex = math.hypot(0.13 / math.log(10), 0.26 / math.log(10))
    falls = [fall * 10 ** (s * 2.0 * sd_fall_dex) for s in (-1.0, -0.5, 0.0, 0.5, 1.0)]
    lam0s = [GROSCH_LAMBDA_100C["thiamine"][1], lam0, GROSCH_LAMBDA_100C["thiamine"][2]]
    cap = {}
    for n in (1, 2):
        grid = [dict(capture_law(n, l0, fl), lambda0=l0, fall=fl) for l0 in lam0s for fl in falls]
        cap[f"n_{n}"] = {"central": capture_law(n, lam0, fall), "band_scan": grid,
                         "refuted_everywhere_in_bands": all(g["refuted"] for g in grid),
                         "survives_at": [(g["lambda0"], g["fall"]) for g in grid if not g["refuted"]]}
    out["capture"] = {"lambda0_cys_free": lam0, "fall_325_to_2720": fall, "fall_2sd_dex": 2.0 * sd_fall_dex,
                      "table_viii_band_with_cys": GROSCH_LAMBDA_100C["thiamine_cys"], "laws": cap}
    # the thiamine charge caps any thiamine-derived dicarbonyl-assisted release of H2S from Cys in no. 4 (no ribose)
    out["dicarbonyl_cap_no4_uM"] = 1000.0 * GROSCH["thiamine_umol"] / 1000.0 / vol
    t3c, t3m = GROSCH["T3_cys_umol"], GROSCH["T3_thi30"]
    out["high_dose_log_slope"] = {
        "thi30_680_to_2720": math.log(t3m[-1] / t3m[6]) / math.log(t3c[-1] / t3c[6]),
        "thi6_680_to_1360": math.log(0.5 / 1.6) / math.log(1360 / 680),
    }
    # CK-42b: ratio vs concentration
    thi30_opt = [c for c, m in zip(t3c, t3m) if m >= 0.9 * max(t3m)]
    thi6 = GROSCH["T3_thi6"]
    thi6_opt = max(thi6, key=thi6.get)
    out["ck42b_ratio_vs_concentration"] = {
        "thi30_optimum_umol_within_10pct_of_max": thi30_opt,
        "thi6_printed_points": thi6, "thi6_argmax_umol": thi6_opt,
        "ratio_law_predicted_thi6_optimum_umol": [round(0.2 * c) for c in (170, 325)],
        "concentration_law_predicted_thi6_optimum_umol": [170, 340],
        "thi6_series_lowest_printed_point_umol": min(thi6),
    }
    # total-thiol concentration law from Table III (thiamine 30 umol), log-log interpolation in thiol mM
    xs = [c / 1000 / vol for c in t3c[1:]]
    ys = list(t3m[1:])
    # Table III has two near-duplicate doses 325 / 340 with 19 / 14: use their mean at 332.5 umol (declared)
    pts = sorted(zip(xs, ys))
    merged: List[List[float]] = []
    for x, y in pts:
        if merged and abs(x - merged[-1][0]) / x < 0.06:
            merged[-1] = [0.5 * (merged[-1][0] + x), 0.5 * (merged[-1][1] + y)]
        else:
            merged.append([x, y])

    def interp(c_mm: float) -> float:
        lx = math.log(c_mm)
        for (x0, y0), (x1, y1) in zip(merged, merged[1:]):
            if x0 <= c_mm <= x1:
                w = (lx - math.log(x0)) / (math.log(x1) - math.log(x0))
                return math.exp(math.log(y0) + w * (math.log(y1) - math.log(y0)))
        raise ValueError(c_mm)

    preds = {}
    for name, thiol_umol in (("no7_thi_gsh_rib", GROSCH["gsh_umol"]),
                             ("no8_thi_gsh_cys_rib", GROSCH["gsh_umol"] + GROSCH["cys_T2_umol"]),
                             ("no6_thi_cys_rib", GROSCH["cys_T2_umol"])):
        mean, sd, n = GROSCH["T2"][name]
        c = thiol_umol / 1000 / vol
        pred = interp(c)
        tol = 2.0 * math.sqrt((sd / mean / math.log(10)) ** 2 + 0.05 ** 2)
        resid = math.log10(pred / mean)
        preds[name] = {"total_thiol_mM": c, "pred_ug": pred, "meas_ug": mean, "sd": sd, "n": n,
                       "resid_dex": resid, "tol_dex": tol, "pass": abs(resid) <= tol}
    out["ck42b_total_thiol_law"] = {"interpolation_nodes_mM_ug": merged, "predictions": preds}
    # pH artefact: Cys.HCl acid load in 0.2 M pyrophosphate (Henderson-Hasselbalch on the third proton)
    c_buf, pka, ph0 = PYROPHOSPHATE["conc_M"], PYROPHOSPHATE["pKa3"], PYROPHOSPHATE["pH"]
    hcl = 2720e-6 / vol
    a_base = c_buf / (1.0 + 10 ** (pka - ph0))      # HP2O7(3-) form
    a_acid = c_buf - a_base
    ph1 = pka + math.log10((a_base - hcl) / (a_acid + hcl))
    out["ck42b_ph_artefact"] = {"HCl_M": hcl, "pH_before": ph0, "pH_after": ph1, "delta_pH": ph1 - ph0,
                                "max_effect_dex_at_1_dex_per_pH": abs(ph1 - ph0),
                                "retired": abs(ph1 - ph0) * 1.0 < 0.5}
    return out


# ---------------------------------------------------------------------------------------------------------------
# CK-56 (engine-free)
# ---------------------------------------------------------------------------------------------------------------

def part_ck56() -> Dict[str, Any]:
    kcal = 4.184

    def case1(g: float):
        e = 32.9 + 1.94 * g - 1.32 * g * g        # pH_c = a_w,c = 0
        k100 = 0.0124 + 0.0043 * g
        return e * kcal, k100

    def case2(g: float):
        e = 32.9 + 1.10 * g ** 3 - 1.32 * g * g
        k100 = 0.0122 + 0.00175 * g ** 3
        return e * kcal, k100

    def k_at(e_kj: float, k100: float, t_k: float) -> float:
        return k100 * math.exp(-(e_kj / R_GAS) * (1.0 / t_k - 1.0 / 373.15))

    t105, t125 = 378.15, 398.15
    inv = 1.0 / t105 - 1.0 / t125
    k105_meas = (0.024, 0.026)
    sd_lnk_e = 2.0 * kcal / R_GAS * abs(1 / t125 - 1 / 373.15)   # E's residual SD moves ln k125 by this
    out: Dict[str, Any] = {"engine_Ea": 96.8, "wolf_E_residual_sd_kJ": 2.0 * kcal, "cases": {}}
    for cname, fn in (("case1", case1), ("case2", case2)):
        for g in (1.68, 2.0):
            e, k100 = fn(g)
            k125 = k_at(e, k100, t125)
            k105_surface = k_at(e, k100, t105)
            sd_lnk125 = math.sqrt((0.003 / k100) ** 2 + sd_lnk_e ** 2)
            es = [R_GAS * math.log(k125 / k) / inv for k in k105_meas]
            e_mid = R_GAS * math.log(k125 / 0.025) / inv
            sd_e = R_GAS / inv * math.sqrt(sd_lnk125 ** 2 + (math.log(0.026 / 0.024) / 2) ** 2)
            band = max(2.0 * kcal, sd_e)
            k80 = k_at(e, k100, 353.15)   # the plateau-prone end point (reported, not used by the rule)
            e_80_105 = R_GAS * math.log(0.025 / k80) / (1.0 / 353.15 - 1.0 / t105)
            out["cases"][f"{cname}_G{g:g}"] = {
                "surface_k80": k80, "E_80_105_kJ_with_surface_k80": e_80_105,
                "surface_E_kJ": e, "surface_k100": k100, "surface_k105": k105_surface, "surface_k125": k125,
                "E_105_125_kJ": e_mid, "E_105_125_range_kJ": es, "propagated_sd_kJ": sd_e, "band_kJ": band,
                "explains": abs(e_mid - 96.8) <= band, "refutes_a": (e_mid - 96.8) > 2.0 * band,
            }
    # side reading: the companion's 125 C / 2 min pass 1 (9.6 -> 5.0 mg/g)
    k125_comp = math.log(9.6 / 5.0) / 2.0
    out["companion_side_reading"] = {"k125_lower_bound": k125_comp,
                                     "E_105_125_kJ": R_GAS * math.log(k125_comp / 0.025) / inv,
                                     "note": "composition of the companion's pass not matched to the 105 C run"}
    return out


# ---------------------------------------------------------------------------------------------------------------

PARTS = ("ck23_schieberle", "ck23_guard", "ck34", "ck42", "ck56")


def run_part(name: str, inputs: Optional[Path]) -> Dict[str, Any]:
    t0 = time.time()
    if name == "ck23_schieberle":
        res = part_ck23_schieberle()
    elif name == "ck23_guard":
        res = part_ck23_guard(json.loads(inputs.read_text()))
    elif name == "ck34":
        res = part_ck34()
    elif name == "ck42":
        res = part_ck42()
    elif name == "ck56":
        res = part_ck56()
    else:
        raise SystemExit(f"unknown part {name}")
    res["_seconds"] = time.time() - t0
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--part", choices=PARTS)
    ap.add_argument("--inputs", type=Path)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args(argv)
    if args.part:
        res = run_part(args.part, args.inputs)
        (args.out or Path("/dev/stdout")).write_text(json.dumps(res, indent=1, default=float))
        return 0
    results: Dict[str, Any] = {}
    with tempfile.TemporaryDirectory() as tmp:
        for part in PARTS:
            out = Path(tmp) / f"{part}.json"
            cmd = [sys.executable, str(Path(__file__).resolve()), "--part", part, "--out", str(out)]
            if part == "ck23_guard":
                cmd += ["--inputs", str(Path(tmp) / "ck23_schieberle.json")]
            print("running", part, flush=True)
            subprocess.run(cmd, check=True, cwd=str(ROOT))
            results[part] = json.loads(out.read_text())
    from src import provenance
    payload = {
        "artifact": "m0_checks_batch2_diagnostics",
        "generated_by": "scripts/diagnostics/m0_checks_batch2.py",
        "record": "results/validation/m0_checks_batch2.md",
        "provenance": provenance.provenance_block(
            "m0_checks_batch2_diagnostics", generated_by="scripts/diagnostics/m0_checks_batch2.py", wave="M0",
            inputs=[ROOT / "results" / "validation" / "kinetic_core_b60_fit_report.json"]),
        "parts": results,
    }
    OUT.write_text(json.dumps(payload, indent=1, default=float) + "\n")
    print("wrote", OUT.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
