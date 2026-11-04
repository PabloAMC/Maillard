#!/usr/bin/env python
"""
Plan step M1.0 ("no fit", 2026-11-04): what the SHIPPED engine (R0 = B60, today's default configuration) predicts on
M1's decisive data, before any rival (C1 - C3) is built. No constant and no switch is changed; the process asserts
that its configuration is the shipped one before and after every check.

Case: docs/reference/cases/M1_thiol_loss_oxidant_ph_drift.md (secs. 3 - 5); plan: docs/guides/STATUS_AND_PLAN.md
sec. 5, W2, row M1; conflict CF-107.

  G1  B27's gate at B60. In Zhou 2023's and Zhang 2024's ambient pots (B23.SYSTEMS, the systems B60 was fitted on):
      the fraction of the charged oxidant OX the consumers use over the cook, and the relative rise in that pool if the
      whole mercaptoketone flux delivered one OX each (phi = 1). Reported beside the B9 numbers of
      kinetic_core_b27_ship_rule.json and, for attribution only, the B9 vector under today's numerics.
  D1  Schieberle 2000 Table IV (100 C, 0.5 M phosphate pH 5, Rib 100 + Cys 33 mM): R0's MFT and FFT ratios 360/30,
      720/360 and 720/30 min against the measured ones (the reference pot's own scale cancels).
  D2  Free cysteine and pentose over time in Cao 2017 (90 C), Zhang 2025b (120 C) and Zhai 2019 (70 - 130 C): R0's
      time to half of the charge, free and with the TTCA / cysteine-Amadori reservoir counted, against the measured
      points (fit_rows route_submodel.yml).

Blinding: every row read is DEVELOPMENT in contrast_split_v1 and in split_v2_per_change (checked here, below); the
frozen Maillard-path hold-out bundle of Zhang 2025b (Table 2 volatiles) is not read and no volatile is computed or
reported for that pot.

    python scripts/generators/m1_preflight_checks.py
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
for p in (ROOT, ROOT / "scripts" / "generators"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import generate_kinetic_core_b2_3_fit as B23  # noqa: E402  (the fit systems; no global is set at import)
import generate_kinetic_core_b27_ship_rule as B27SR  # noqa: E402  (G1's pot / consumer lists, _predict, _operative_and_drift)
from src import artifact_io, data_paths, provenance  # noqa: E402
from src.kinetic_core import config, engine, fit_rows  # noqa: E402
from src.kinetic_core.engine import SULFUR, ProcessSpec, ThermalProgram, core_parameters, core_ph_drift  # noqa: E402
from src.kinetic_core.parameters_sulfur import OX_AMBIENT_MMOL_L  # noqa: E402
from src.kinetic_core.ph_state import BufferSpec  # noqa: E402
from src.kinetic_core.sulfur import FULL_REACTIONS, flux_closure, integrate_sulfur, reaction_flux_integrals  # noqa: E402

CELSIUS = 273.15
V = data_paths.VALIDATION_DIR
OUT = V / "kinetic_core_m1_preflight_checks.json"
B60_REPORT = V / "kinetic_core_b60_fit_report.json"
B9_REPORT = V / "kinetic_core_b9_fit_report.json"
B27_SHIP = V / "kinetic_core_b27_ship_rule.json"
S042 = V / "kinetic_core_joint_deoxyosone_thiol_fit_s042.json"
CASE = data_paths.REPO_ROOT / "docs" / "reference" / "cases" / "M1_thiol_loss_oxidant_ph_drift.md"
SPLIT_V1 = data_paths.REPO_ROOT / "docs" / "reference" / "contrast_benchmark" / "split.yml"
SPLIT_V2 = data_paths.REPO_ROOT / "docs" / "reference" / "contrast_benchmark" / "split_v2_per_change.yml"
THIOL_TABLE, ROUTE_TABLE = "thiol_submodel", "route_submodel"

#: the studies whose rows this artifact reads, with their laboratory (contrast_benchmark/lab_ids.yml)
STUDIES = {"schieberle2000": "dfa_garching", "cao2017": "btbu_beijing", "zhang2025b": "btbu_beijing",
           "zhai2019": "ho_network", "zhou2023": "ho_network", "Zhang2024": "btbu_beijing"}
#: held out and therefore NOT read, with the reason (listed in the report)
NOT_READ = [
    {"item": "mp_holdout_xylose_cysteine_phosphate_120C_50min_pH65_Zhang2025b (Table 2 blank: MFT, FFT, MFT disulfide, "
             "2-acetylthiazole at 50 min)",
     "why": "frozen Maillard-path hold-out bundle (FIT_HOLDOUT_DECLARATION.md Amendment 45); D2 integrates that pot but "
            "reads and reports only cysteine, xylose, TTCA and cysteine-Amadori"},
    {"item": "Kumazawa 2003 (F28#3), Esterbauer 1975 (F28#17), Stack 2018 (F28#18)", "why": "contrast_split_v1 TEST; not used"},
    {"item": "Ehrenberg 1989 (F28#8), Xu 2011 (F28#20), Wang 2022 (F28#25), Liu 2023b (F28#26)",
     "why": "split_v2_per_change E12 TEST contrasts; not used"},
]

D1_POT = "schieberle2000_rib_cys_ph5_100C"
D1_MINUTES = (30.0, 60.0, 360.0, 720.0)
D1_RATIOS = ((360.0, 30.0), (720.0, 360.0), (720.0, 30.0))
SCHIEBERLE_BUFFER = BufferSpec(kind="phosphate", phosphate_mol_l=0.5, declared=True,
                               source="Schieberle & Hofmann 2000 (thiol_submodel.yml pot): phosphate 0.5 mol/L, pH 5.0")

PHOSPHATE_02 = lambda src: BufferSpec(kind="phosphate", phosphate_mol_l=0.2, declared=True, source=src)  # noqa: E731
#: D2 pots: (title, measured pot ids in route_submodel.yml, the integrated pot, buffer, charged Cys / pentose mmol/L)
D2_POTS = {
    "cao2017_90C": dict(pots=("cao2017_blank_fig4", "cao2017_blank_fig5"), run="cao2017_blank_fig4",
                        buffer=PHOSPHATE_02("Cao 2017: 0.2 M sodium phosphate, pH 5.5"), charge={"Cys": 400.0, "PENT": 400.0}),
    "zhang2025b_120C": dict(pots=("zhang2025b_blank_fig2gh", "zhang2025b_blank"), run="zhang2025b_blank",
                            buffer=PHOSPHATE_02("Zhang 2025b: 0.2 mol/L phosphate, pH 6.5"), charge={"Cys": 40.0, "PENT": 40.0}),
}
for _t in (70, 90, 110, 130):
    D2_POTS[f"zhai2019_{_t}C"] = dict(
        pots=(f"zhai2019_{_t}c",), run=f"zhai2019_{_t}c",
        buffer=BufferSpec(kind="none", declared=True, source="Zhai 2019: deionized water, pH 7.4 set with NaOH"),
        charge={"Cys": 82.6, "PENT": 165.2})
RESERVOIR = ("TTCA", "CAMA")  # each holds one cysteine and one pentose
D2_SPECIES = ("Cys", "PENT", "TTCA", "CAMA")


# ---------------------------------------------------------------------------------------------------------------
# guards
# ---------------------------------------------------------------------------------------------------------------

def assert_shipped() -> str:
    """The process runs the shipped configuration (no switch moved by an import or a check); returns its digest."""
    here = config.Config.from_globals()
    diff = here.diff(config.Config.shipped())
    if diff:
        raise SystemExit(f"not the shipped configuration: {sorted(diff)}")
    return here.digest()


def split_status() -> Dict[str, Any]:
    """Every study read is DEVELOPMENT in both frozen splits; raise otherwise."""
    v1 = yaml.safe_load(SPLIT_V1.read_text(encoding="utf-8"))
    v2 = yaml.safe_load(SPLIT_V2.read_text(encoding="utf-8"))
    out, problems = {}, []
    v1_contrasts = {c["id"]: c for c in v1.get("contrasts") or []}
    for study, lab in STUDIES.items():
        rec: Dict[str, Any] = {"lab": lab}
        rec["v1_lab"] = ("development" if lab in v1["split"]["development_labs"]
                         else "test" if lab in v1["split"]["test_labs"] else "absent")
        rec["v1_contrasts"] = {cid: c["partition"] for cid, c in v1_contrasts.items() if cid.endswith("@" + study)}
        rec["v2_test_contrasts"] = sorted(f"{e}:{c}" for e, blk in v2["split"]["errors"].items()
                                          for c in (blk.get("test_contrasts") or []) if c.endswith("@" + study))
        rec["v2_calibration_input"] = study in (v2["split"].get("calibration_input_papers") or [])
        if rec["v1_lab"] == "test" or "test" in rec["v1_contrasts"].values() or rec["v2_test_contrasts"]:
            problems.append(study)
        out[study] = rec
    if problems:
        raise SystemExit(f"TEST rows would be read: {problems}")
    out["_v2_E10_development_labs"] = v2["split"]["errors"]["E10"]["development_labs"]
    out["_v2_E12_development_labs"] = v2["split"]["errors"]["E12"]["development_labs"]
    return out


# ---------------------------------------------------------------------------------------------------------------
# G1
# ---------------------------------------------------------------------------------------------------------------

OX_CONSUMING = tuple((r.key, r.reactants["OX"]) for r in FULL_REACTIONS if r.reactants.get("OX"))
OX_PRODUCING = tuple((r.key, r.products["OX"]) for r in FULL_REACTIONS if r.products.get("OX"))


def g1_pot(params: Mapping[str, Any], drift, name: str) -> Dict[str, Any]:
    """B27's G1 probe on one fit system, on B60's closed flux budget (sulfur.sulfur_flux_budget's grid), with the run
    kept so that the OX budget can be checked against the state."""
    s = B23.SYSTEMS[name]
    dur = float(s["minutes"])
    grid = np.unique(np.concatenate([np.linspace(0.0, dur, 201), np.geomspace(dur * 1e-7, dur, 201)]))
    run = integrate_sulfur(params, float(s["t_c"]) + CELSIUS, s["initial"], grid, ph=float(s["ph"]),
                           buffer_spec=s.get("buffer"), ph_drift=drift)
    fl = reaction_flux_integrals(run)
    ox0 = float(s["initial"].get("OX", 0.0))
    used_b27 = sum(float(fl.get(k, 0.0)) for k in B27SR.OX_CONSUMERS)
    used_all = sum(nu * float(fl.get(k, 0.0)) for k, nu in OX_CONSUMING)
    made = sum(nu * float(fl.get(k, 0.0)) for k, nu in OX_PRODUCING)
    supply_b27 = sum(float(fl.get(k, 0.0)) for k in B27SR.MERCAPTOKETONE)
    supply_nf = float(fl.get("r_nf_mp3p", 0.0))
    closure = flux_closure(run, fl, ("OX", "MFT", "FFT", "Cys"))
    return {
        "system": name, "anchor": s.get("anchor"), "t_c": s["t_c"], "minutes": dur, "ph_label": s["ph"],
        "ox_charged_mmol_l": ox0, "ox_final_mmol_l": run.final("OX"),
        "ox_consumed_b27_consumers": used_b27, "ox_consumed_all_reactions": used_all, "ox_made": made,
        "consumed_fraction": used_all / ox0 if ox0 else None,
        "consumed_fraction_b27_list": used_b27 / ox0 if ox0 else None,
        "consumers": {k: float(fl.get(k, 0.0)) for k, _ in OX_CONSUMING if fl.get(k, 0.0)},
        "mercaptoketone_flux_b27_sum_mmol_l": supply_b27, "r_nf_mp3p_flux_mmol_l": supply_nf,
        "rise_at_phi_1_b27_sum": supply_b27 / ox0 if ox0 else None,
        "rise_at_phi_1_r_nf_mp3p_only": supply_nf / ox0 if ox0 else None,
        "dimer_flux_mmol_l": float(fl.get("ch_dimer_mft", 0.0)) + float(fl.get("ch_dimer_fft", 0.0)),
        "flux_closure": closure,
        "ph_final_in_situ": run.ph_final,
        "ph_fixed_point_converged": run.metadata.get("ph_fixed_point_converged"),
    }


def g1(params, drift, params_b9, drift_b9) -> Dict[str, Any]:
    b9_rec = json.loads(B27_SHIP.read_text(encoding="utf-8"))["G1"]["pots"]
    pots = {}
    for name in B27SR.AMBIENT_POTS:
        b60 = g1_pot(params, drift, name)
        b9_now = g1_pot(params_b9, drift_b9, name)
        rec = b9_rec[name]
        pots[name] = {"B60": b60,
                      "B9_as_recorded": {"consumed_fraction": rec["consumed_fraction"],
                                         "rise_at_phi_1": rec["supply_at_phi_1_over_pool"],
                                         "dimer_flux_mmol_l": rec["dimer_flux"],
                                         "basis": "kinetic_core_b27_ship_rule.json (2026-09-11: B9 vector, pre-B51 pH state, "
                                                  "pre-B60 linear flux budget)"},
                      "B9_vector_today_numerics": {k: b9_now[k] for k in (
                          "consumed_fraction", "consumed_fraction_b27_list", "rise_at_phi_1_b27_sum",
                          "rise_at_phi_1_r_nf_mp3p_only", "dimer_flux_mmol_l")}}
    best = max(v["B60"]["rise_at_phi_1_b27_sum"] or 0.0 for v in pots.values())
    worst_used = max(v["B60"]["consumed_fraction"] or 0.0 for v in pots.values())
    holds = bool(worst_used < 0.01 and best < 0.01)
    return {"pots": pots, "needed_rise_for_one_decade": 9.0, "max_consumed_fraction_B60": worst_used,
            "max_rise_at_phi_1_B60": best, "holds": holds,
            "criterion": "G1 holds when the consumers use < 1 % of the charged OX in every ambient pot and the whole "
                         "mercaptoketone flux at phi = 1 would raise the pool by < 1 % (B27: < 1 % and < 0.5 %; a decade on "
                         "the share needs ~9x)"}


# ---------------------------------------------------------------------------------------------------------------
# D1
# ---------------------------------------------------------------------------------------------------------------

def _first_below(t: np.ndarray, f: np.ndarray, level: float = 0.5) -> Optional[float]:
    """First time f falls to ``level`` (f(t) linear in log t between the bracketing points); None if never."""
    idx = np.nonzero(f <= level)[0]
    if idx.size == 0:
        return None
    i = int(idx[0])
    if i == 0:
        return float(t[0])
    t0, t1, f0, f1 = float(t[i - 1]), float(t[i]), float(f[i - 1]), float(f[i])
    if t0 <= 0.0:
        return t1
    w = (f0 - level) / (f0 - f1) if f0 != f1 else 1.0
    return float(math.exp(math.log(t0) + w * (math.log(t1) - math.log(t0))))


def dense_run(params, drift, initial: Mapping[str, float], t_c: float, minutes: float, ph: float, buffer: BufferSpec,
              extra_times: Sequence[float] = (), check_species: Sequence[str] = ("Cys", "PENT")):
    """The engine's sulfur path (engine._integrate_program's state and arguments: ambient OX, the declared headspace
    reservoir, rtol 1e-8) on a dense log-and-linear output grid."""
    process = ProcessSpec(ThermalProgram.isothermal(float(t_c), float(minutes)), ph=float(ph), buffer=buffer)
    state = dict(initial)
    state.setdefault("OX", OX_AMBIENT_MMOL_L)
    state.setdefault("OXR", engine.oxygen_reservoir_units(process)[0])
    state.setdefault("OXV", 0.0)
    grid = np.unique(np.concatenate([[0.0], np.geomspace(minutes * 1e-9, minutes, 601), np.linspace(0.0, minutes, 301),
                                     np.asarray(extra_times, dtype=float)]))
    run = integrate_sulfur(params, float(t_c) + CELSIUS, state, grid, ph=float(ph), buffer_spec=buffer, ph_drift=drift,
                           rtol=1e-8, atol=1e-14)
    # the endpoint against the engine's own path (one segment, grid [0, minutes])
    final, meta = engine._integrate_program(SULFUR, params, dict(initial), process, ph_drift=drift)
    check = {k: abs(math.log10((run.final(k) + 1e-30) / (final[k] + 1e-30))) for k in check_species}
    return run, check, meta


def d1(params, drift) -> Dict[str, Any]:
    _pots, rows = fit_rows.load(THIOL_TABLE)
    meas = {(r.species, r.t_min): r for r in rows if r.pot == D1_POT}
    pred = {m: B27SR._predict({"D-Ribose": 100.0, "L-Cysteine": 33.0}, 100.0, m, 5.0, SCHIEBERLE_BUFFER, params, drift,
                              ("MFT", "FFT")) for m in D1_MINUTES}
    run, check, _meta = dense_run(params, drift, {"PENT": 100.0, "Cys": 33.0}, 100.0, 720.0, 5.0, SCHIEBERLE_BUFFER,
                                  extra_times=D1_MINUTES, check_species=("Cys", "PENT", "MFT", "FFT"))
    b9 = json.loads(B27_SHIP.read_text(encoding="utf-8"))["b9_reference"]["T1"]
    out: Dict[str, Any] = {"pot": D1_POT, "species": {}}
    for s in ("MFT", "FFT"):
        series = run.series(s)
        imax = int(np.argmax(series))
        rec: Dict[str, Any] = {
            "measured_ug_l": {int(m): meas[(s, m)].value for m in D1_MINUTES},
            "measured_read_error_dex": {int(m): meas[(s, m)].sigma for m in D1_MINUTES},
            "R0_ug_l": {int(m): pred[m][s] for m in D1_MINUTES},
            "R0_peak_min": float(run.times_min[imax]),
            "R0_peak_over_720": float(series[imax] / series[-1]) if series[-1] > 0 else None,
            "R0_still_rising_6_to_12_h": bool(pred[720.0][s] > pred[360.0][s]),
            "measured_still_rising_6_to_12_h": bool(meas[(s, 720.0)].value > meas[(s, 360.0)].value),
            "ratios": {},
        }
        for a, b in D1_RATIOS:
            m_ratio = meas[(s, a)].value / meas[(s, b)].value
            p_ratio = pred[a][s] / pred[b][s] if pred[b][s] else float("inf")
            err = math.hypot(meas[(s, a)].sigma, meas[(s, b)].sigma)
            resid = math.log10(p_ratio / m_ratio) if p_ratio > 0 and math.isfinite(p_ratio) else float("inf")
            key = f"{int(a)}_over_{int(b)}"
            b9_key = f"{s}_{int(a)}_over_30"
            rec["ratios"][key] = {"measured": m_ratio, "measured_dex": math.log10(m_ratio), "R0": p_ratio,
                                  "R0_dex": math.log10(p_ratio) if p_ratio > 0 else None,
                                  "residual_dex": resid, "read_error_dex": err,
                                  "residual_over_read_error": resid / err if err else None,
                                  "B9_residual_dex": b9["ratios_dex"].get(b9_key) if int(b) == 30 else None}
        out["species"][s] = rec
    out["dense_vs_engine_endpoint_dex"] = check
    s042 = json.loads(S042.read_text(encoding="utf-8")).get("series_extra_noise", {})
    out["noise"] = {
        "source": "Table IV prints no SD and no n; the only stated error is the transcription error (half the last printed "
                  "digit), 0.001 - 0.011 dex per point (thiol_submodel.yml); a ratio's error is the two in quadrature.",
        "measurement_model_series_noise_dex": {k.split("|")[-1]: v for k, v in s042.items() if D1_POT in k},
        "measurement_model_basis": "the per-series extra noise the joint thiol fit s042 estimated from its own residual scatter "
                                   "(fitting.estimate_series_noise): a property of that fit's misfit, not a measured "
                                   "replicate error; quoted as the case's sec. 5 says, not used to score",
    }
    out["B9_as_recorded"] = {"ratios_dex": b9["ratios_dex"], "still_rising_6_to_12_h": b9["still_rising_6_to_12_h"],
                             "basis": "kinetic_core_b27_ship_rule.json b9_reference.T1 (residuals log10 model/measured)"}
    return out


# ---------------------------------------------------------------------------------------------------------------
# D2
# ---------------------------------------------------------------------------------------------------------------

def _measured_half(times: Sequence[float], fracs: Sequence[float]) -> Dict[str, Any]:
    t, f = np.asarray(times, dtype=float), np.asarray(fracs, dtype=float)
    order = np.argsort(t)
    t, f = t[order], f[order]
    if f.size == 0:
        return {"t_half_min": None, "statement": "not measured"}
    if f[0] <= 0.5:
        return {"t_half_min": None, "upper_bound_min": float(t[0]),
                "statement": f"< {t[0]:g} min (already {f[0]:.2f} at the first measured time)"}
    th = _first_below(t, f)
    if th is None:
        return {"t_half_min": None, "lower_bound_min": float(t[-1]),
                "statement": f"> {t[-1]:g} min ({f[-1]:.2f} left at the last measured time)"}
    return {"t_half_min": th, "statement": f"~{th:.3g} min (log-time interpolation between measured points)"}


def _predicted_half(t: np.ndarray, f: np.ndarray, minutes: float) -> Dict[str, Any]:
    th = _first_below(t, f)
    if th is None:
        return {"t_half_min": None, "statement": f"> {minutes:g} min ({f[-1]:.2f} left at the end)"}
    return {"t_half_min": th, "t_half_s": th * 60.0, "statement": f"{th * 60.0:.3g} s" if th < 1.0 else f"{th:.3g} min"}


def d2(params, drift) -> Dict[str, Any]:
    pots, rows = fit_rows.load(ROUTE_TABLE)
    out = {}
    for title, spec in D2_POTS.items():
        pot = pots[spec["run"]]
        initial = {k: float(v) for k, v in pot["initial_mmol_L"].items()}
        charge = spec["charge"]
        my_rows = [r for r in rows if r.pot in spec["pots"] and r.species in D2_SPECIES]
        mtimes = sorted({r.t_min for r in my_rows})
        run, check, meta = dense_run(params, drift, initial, float(pot["t_c"]), float(pot["minutes"]),
                                     float(pot["label_ph"]), spec["buffer"], extra_times=mtimes)
        t = run.times_min
        series = {k: run.series(k) for k in D2_SPECIES}
        res = series["TTCA"] + series["CAMA"]
        pred_frac = {"Cys_free": series["Cys"] / charge["Cys"], "PENT_free": series["PENT"] / charge["PENT"],
                     "Cys_free_plus_reservoir": (series["Cys"] + res) / charge["Cys"],
                     "PENT_free_plus_reservoir": (series["PENT"] + res) / charge["PENT"]}
        meas = {(r.species, r.t_min): r for r in my_rows}

        def m_frac(species_free: str, with_reservoir: bool):
            ts, fs = [], []
            for tm in mtimes:
                parts = [(species_free, tm)] + ([("TTCA", tm), ("CAMA", tm)] if with_reservoir else [])
                if all(p in meas for p in parts):
                    ts.append(tm)
                    fs.append(sum(meas[p].value for p in parts) / charge["Cys" if species_free == "Cys" else "PENT"])
            return ts, fs

        def at(times, frac):
            return {f"{tm:g}": float(np.interp(tm, t, frac)) for tm in times}

        rec: Dict[str, Any] = {"pot": spec["run"], "measured_pots": list(spec["pots"]), "t_c": pot["t_c"],
                               "minutes": pot["minutes"], "label_ph": pot["label_ph"], "initial_mmol_l": initial,
                               "charge_mmol_l": charge, "buffer": spec["buffer"].source,
                               "R0_final_ph_in_situ": run.ph_final,
                               "R0_ph_fixed_point_converged": run.metadata.get("ph_fixed_point_converged"),
                               "dense_vs_engine_endpoint_dex": check,
                               "pools": {}}
        for pool, free, with_res in (("cysteine_free", "Cys", False), ("pentose_free", "PENT", False),
                                     ("cysteine_free_plus_reservoir", "Cys", True),
                                     ("pentose_free_plus_reservoir", "PENT", True)):
            key = {"cysteine_free": "Cys_free", "pentose_free": "PENT_free",
                   "cysteine_free_plus_reservoir": "Cys_free_plus_reservoir",
                   "pentose_free_plus_reservoir": "PENT_free_plus_reservoir"}[pool]
            ts, fs = m_frac(free, with_res)
            entry = {"measured_fraction_of_charge": {f"{a:g}": b for a, b in zip(ts, fs)},
                     "measured_half": _measured_half(ts, fs),
                     "R0_fraction_at_measured_times": at(ts, pred_frac[key]) if ts else {},
                     "R0_half": _predicted_half(t, pred_frac[key], float(pot["minutes"]))}
            if ts:
                entry["log10_R0_over_measured"] = {f"{a:g}": (math.log10(max(float(np.interp(a, t, pred_frac[key])), 1e-30) / b)
                                                             if b > 0 else None) for a, b in zip(ts, fs)}
            rec["pools"][pool] = entry
        # Zhai 2019 measures only the reservoir: a LOWER bound on the cysteine not destroyed (free Cys >= 0)
        res_ts = [tm for tm in mtimes if ("TTCA", tm) in meas]
        rec["reservoir_measured_mmol_l"] = {f"{tm:g}": meas[("TTCA", tm)].value + (meas[("CAMA", tm)].value if ("CAMA", tm) in meas else 0.0)
                                            for tm in res_ts}
        rec["reservoir_R0_mmol_l"] = {f"{tm:g}": float(np.interp(tm, t, res)) for tm in res_ts}
        rec["reservoir_lower_bound_on_cysteine_fraction"] = {f"{tm:g}": v / charge["Cys"] for tm, v in
                                                             zip(res_ts, rec["reservoir_measured_mmol_l"].values())}
        rec["R0_cysteine_free_plus_reservoir_fraction"] = at(res_ts, pred_frac["Cys_free_plus_reservoir"])
        rec["R0_pentose_free_plus_reservoir_fraction"] = at(res_ts, pred_frac["PENT_free_plus_reservoir"])
        rec["measurement_noise"] = sorted({f"{r.species}: sigma {r.sigma:.3f} dex"
                                           + (f" (sd {r.raw.get('sd')})" if r.raw.get("sd") is not None else " (read error only)")
                                           for r in my_rows if r.t_min == (mtimes[-1] if mtimes else None)})
        out[title] = rec
    return out


# ---------------------------------------------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------------------------------------------

def _pct(x: Optional[float], nd: int = 3) -> str:
    return "n/a" if x is None else f"{100.0 * x:.{nd}f} %"


def _d1_reading(d1r: Mapping[str, Any]) -> str:
    m, f = d1r["species"]["MFT"], d1r["species"]["FFT"]
    nz = d1r["noise"]["measurement_model_series_noise_dex"]
    mr, fr = m["ratios"], f["ratios"]
    return (
        f"**Reading.** MFT: R0 rises {mr['360_over_30']['R0']:.2g}x from 30 min to 6 h where the pot rises "
        f"{mr['360_over_30']['measured']:.2g}x, peaks at {m['R0_peak_min'] / 60:.1f} h and falls from 6 to 12 h "
        f"(x{mr['720_over_360']['R0']:.2f}; measured x{mr['720_over_360']['measured']:.2f}, still rising). Residuals "
        f"{mr['360_over_30']['residual_dex']:+.2f} / {mr['720_over_360']['residual_dex']:+.2f} / "
        f"{mr['720_over_30']['residual_dex']:+.2f} dex (360/30, 720/360, 720/30), against B9's "
        f"{mr['360_over_30']['B9_residual_dex']:+.2f} on 360/30; the s042 series noise is {nz.get('MFT', float('nan')):.2f} dex. "
        f"FFT: R0 rises FASTER than the pot from 30 min to 6 h ({fr['360_over_30']['residual_dex']:+.2f} dex), then is flat "
        f"(x{fr['720_over_360']['R0']:.2f} 6 -> 12 h; measured x{fr['720_over_360']['measured']:.2f}), peak at "
        f"{f['R0_peak_min'] / 60:.1f} h; both FFT residuals lie inside the s042 FFT series noise "
        f"({nz.get('FFT', float('nan')):.2f} dex). Every residual is far outside the printed read error, which is the "
        "only error the source states. Net: R0 at B60 still does not reproduce the reference pot's late MFT rise; its FFT "
        "shape is closer than at B9 but plateaus where the pot keeps rising.")


def _d2_reading(d2r: Mapping[str, Any]) -> str:
    cao, zh = d2r["cao2017_90C"]["pools"], d2r["zhang2025b_120C"]["pools"]
    zhai = {k: v for k, v in d2r.items() if k.startswith("zhai")}
    zt = lambda pools, key: pools[key]["R0_half"]["statement"]  # noqa: E731
    cao_last = list(cao["cysteine_free"]["R0_fraction_at_measured_times"])[-1]
    zhai_lines = "; ".join(
        f"{k.split('_')[1]}: free Cys halves in {v['pools']['cysteine_free']['R0_half']['statement']}, pentose in "
        f"{v['pools']['pentose_free']['R0_half']['statement']}; at 160 min R0 holds "
        f"{list(v['R0_cysteine_free_plus_reservoir_fraction'].values())[-1]:.2g} of the cysteine (free + reservoir) "
        f"where the reservoir alone shows >= {list(v['reservoir_lower_bound_on_cysteine_fraction'].values())[-1]:.2f}"
        for k, v in zhai.items())
    return (
        f"**Reading.** The pentose: R0 consumes it in seconds in every pot (Cao 90 C {zt(cao, 'pentose_free')}, Zhang 120 C "
        f"{zt(zh, 'pentose_free')}), where Cao keeps {cao['pentose_free_plus_reservoir']['measured_fraction_of_charge'][cao_last]:.2f} "
        f"of it (free + reservoir) at 5 h and Zhang {list(zh['pentose_free_plus_reservoir']['measured_fraction_of_charge'].values())[-1]:.2f} "
        f"at 50 min. The cysteine: in Zhang's pot R0 consumes it in {zt(zh, 'cysteine_free')} (measured free half-life "
        f"{zh['cysteine_free']['measured_half']['statement']}; {list(zh['cysteine_free_plus_reservoir']['measured_fraction_of_charge'].values())[-1]:.2f} "
        f"still accounted at 50 min). In Cao's equimolar 0.4 M pot R0 does NOT consume it in seconds: the pentose runs out "
        f"first and {cao['cysteine_free']['R0_fraction_at_measured_times'][cao_last]:.2f} of the cysteine is left free at 5 h, "
        f"i.e. a third destroyed where the pot destroys none measurably ({cao['cysteine_free_plus_reservoir']['measured_fraction_of_charge'][cao_last]:.2f} "
        f"accounted); the pot's FREE cysteine is lower than R0's ({cao['cysteine_free']['measured_fraction_of_charge'][cao_last]:.2f}) because "
        f"it sits in TTCA, which R0 does not form. Zhai 2019 ({zhai_lines}). R0's reservoir there is only the measured "
        "t = 0 TTCA (14.5 mM) dissociating; it never forms more. At 130 C the measured reservoir also empties by 160 min, so "
        "the bound says little late in that pot.")


def render(p: Mapping[str, Any]) -> str:
    g, d1r, d2r = p["G1"], p["D1"], p["D2"]
    L: List[str] = [
        "# M1.0: the shipped engine (R0 = B60) on M1's decisive data, no fit", "",
        f"*Plan step M1.0 (docs/guides/STATUS_AND_PLAN.md sec. 5, W2, row M1); case `{p['case']}` secs. 3 - 5; CF-107. "
        f"Configuration: the shipped one, asserted before and after (digest `{p['config_digest'][:16]}`). B60 vector from "
        f"`{p['b60_report']}`. No constant, switch or row was changed. Wall time {p['wall_seconds']:.0f} s.*", "",
        "These are R0's numbers only. Which rival (C1 - C3) the data favour is not stated here; that is M1's comparison.", "",
        "## G1: B27's gate at B60", "",
        "In each ambient pot (the fit systems B60 was fitted on, `B23.SYSTEMS`), the fraction of the charged oxidant OX "
        "(1.0 mM ambient; 62.4 mM cystine in Zhang's GCys arm) that the OX-consuming steps use over the cook, and the "
        "relative rise in the pool if the whole mercaptoketone flux delivered one OX each (phi = 1; B27's sum "
        "r_nf_mp3p + r_mgo_mp + ch_redox_mp3p, and r_nf_mp3p alone). The rise is an upper bound: it books the whole "
        "cook's supply as if present from t = 0.", "",
        "| pot | OX used, B9 (recorded) | OX used, B60 | rise at phi = 1, B9 (recorded) | rise at phi = 1, B60 (sum / r_nf_mp3p only) | B9 vector, today's numerics: used / rise |",
        "|---|---|---|---|---|---|",
    ]
    for name, v in g["pots"].items():
        b60, b9, b9n = v["B60"], v["B9_as_recorded"], v["B9_vector_today_numerics"]
        L.append(f"| {name} | {_pct(b9['consumed_fraction'])} | {_pct(b60['consumed_fraction'])} | {_pct(b9['rise_at_phi_1'], 2)} | "
                 f"{_pct(b60['rise_at_phi_1_b27_sum'], 2)} / {_pct(b60['rise_at_phi_1_r_nf_mp3p_only'], 2)} | "
                 f"{_pct(b9n['consumed_fraction'])} / {_pct(b9n['rise_at_phi_1_b27_sum'], 2)} |")
    L += ["",
          f"**G1 {'HOLDS' if g['holds'] else 'DOES NOT HOLD'} at B60.** Largest use of the pool {_pct(g['max_consumed_fraction_B60'])}; "
          f"largest rise at phi = 1 {_pct(g['max_rise_at_phi_1_B60'], 2)}, against the ~9x (900 %) a decade on the disulfide "
          "share would need. " + g["criterion"] + ".", "",
          "OX consumers carrying flux at B60 (mmol/L over the cook): "
          + "; ".join(f"{n}: " + ", ".join(f"{k} {val:.3g}" for k, val in v['B60']['consumers'].items()) for n, v in g["pots"].items())
          + ". OX budget closure (|change - net flux| / scale): "
          + "; ".join(f"{n} {v['B60']['flux_closure']['OX']:.1e}" for n, v in g["pots"].items()) + ".", "",
          "**Reading.** As at B9, the ambient pots are not oxidant-limited: the dimer steps use a negligible share of a "
          "pool that is a constant multiplier on their rate, and no oxidant SOURCE of the mercaptoketone's size can move "
          "that multiplier. The B9 numbers were computed before B51's pH state and B60's closed flux budget; the "
          "'B9 vector, today's numerics' column separates the vector's change from the numerics': under today's numerics "
          "the B9 vector uses up to "
          + _pct(max(v["B9_vector_today_numerics"]["consumed_fraction"] or 0.0 for v in g["pots"].values()))
          + " of the pool, so the fall to B60's numbers comes from the B60 vector, not from the numerics (which constant "
          "carries it was not resolved here).", "",
          "## D1: Schieberle 2000 Table IV, the reference pot's own ratios", "",
          "100 C, 0.5 M phosphate pH 5, ribose 100 + cysteine 33 mM, one pot read at 30 / 60 / 360 / 720 min. "
          "Residual = log10(R0 ratio / measured ratio).", "",
          "| species | ratio | measured | R0 | residual (dex) | read error (dex) | B9 residual (recorded) |",
          "|---|---|---|---|---|---|---|"]
    for s, rec in d1r["species"].items():
        for key, r in rec["ratios"].items():
            b9 = "--" if r["B9_residual_dex"] is None else f"{r['B9_residual_dex']:+.2f}"
            L.append(f"| {s} | {key.replace('_over_', ' / ')} min | {r['measured']:.3g} ({r['measured_dex']:+.2f} dex) | "
                     f"{r['R0']:.3g} ({r['R0_dex']:+.2f} dex) | {r['residual_dex']:+.2f} | {r['read_error_dex']:.4f} | {b9} |")
    L += [""]
    for s, rec in d1r["species"].items():
        L.append(f"- {s}: measured {', '.join(f'{v:g}' for v in rec['measured_ug_l'].values())} ug/L; R0 "
                 f"{', '.join(f'{v:.3g}' for v in rec['R0_ug_l'].values())} ug/L at 30 / 60 / 360 / 720 min. R0 peaks at "
                 f"{rec['R0_peak_min']:.0f} min ({rec['R0_peak_over_720']:.3g}x its 12 h value); still rising 6 -> 12 h: "
                 f"R0 {rec['R0_still_rising_6_to_12_h']}, measured {rec['measured_still_rising_6_to_12_h']}.")
    nz = d1r["noise"]
    L += ["", f"**Noise.** {nz['source']} Measurement-model series noise (s042): "
          + ", ".join(f"{k} {v:.2f} dex" for k, v in nz["measurement_model_series_noise_dex"].items())
          + f" -- {nz['measurement_model_basis']}.", "",
          _d1_reading(d1r), "",
          "## D2: cysteine and pentose over time", "",
          "Fractions of the charge. 'Free' is the free amino acid or sugar; '+ reservoir' adds TTCA and the cysteine "
          "Amadori compound (each holds one cysteine and one pentose), i.e. the material not yet destroyed. R0 forms "
          "neither (TTCA formation and the cysteine-Amadori route are off in the shipped configuration).", "",
          "| pot | pool | measured: time to half | R0: time to half | measured at last time | R0 at last time |",
          "|---|---|---|---|---|---|"]
    for title, rec in d2r.items():
        for pool, e in rec["pools"].items():
            if not e["measured_fraction_of_charge"] and not title.startswith("zhai"):
                continue
            mf = e["measured_fraction_of_charge"]
            last = list(mf)[-1] if mf else None
            pf = e["R0_fraction_at_measured_times"]
            L.append(f"| {title} | {pool.replace('_', ' ')} | {e['measured_half']['statement']} | {e['R0_half']['statement']} | "
                     f"{'--' if last is None else f'{mf[last]:.3f} at {last} min'} | "
                     f"{'--' if last is None else f'{pf[last]:.2e}'} |")
    L += ["", "Zhai 2019 measures TTCA and the cysteine Amadori compound only (no free cysteine or xylose). Their sum is a "
          "LOWER bound on the cysteine not destroyed (fraction of the 82.6 mM charge), against R0's reservoir:", "",
          "| pot | measured reservoir / charge (lower bound on cysteine left; first, last time) | R0 reservoir (mmol/L; first, last time) | R0 cysteine free + reservoir / charge (first, last time) |",
          "|---|---|---|---|"]
    for title, rec in d2r.items():
        if not title.startswith("zhai"):
            continue
        lb, r0 = rec["reservoir_lower_bound_on_cysteine_fraction"], rec["reservoir_R0_mmol_l"]
        r0c = rec["R0_cysteine_free_plus_reservoir_fraction"]
        ks = list(lb)
        L.append(f"| {title} | {lb[ks[0]]:.2f} at {ks[0]} min, {lb[ks[-1]]:.2f} at {ks[-1]} min | "
                 f"{r0[ks[0]]:.3g}, {r0[ks[-1]]:.3g} | {r0c[ks[0]]:.3g}, {r0c[ks[-1]]:.3g} |")
    L += ["", "Endpoint check, dense grid vs the engine's own one-segment path (max |dex| over Cys, PENT): "
          + "; ".join(f"{k} {max(v['dense_vs_engine_endpoint_dex'].values()):.1e}" for k, v in d2r.items())
          + f"; Schieberle pot (Cys, PENT, MFT, FFT) {max(d1r['dense_vs_engine_endpoint_dex'].values()):.1e}.", "",
          "**Noise.** Cao 2017 and Zhang 2025b: printed SD (n = 3) where an error bar is drawn, plus the digitisation error "
          "(route_submodel.yml); Zhang 2025b's cysteine and xylose calibrations print identical (CF-063) and its free "
          "series are read as ln(C/C0) with C0 taken as the charge. Zhai 2019: raster read error only (no SD printed).", "",
          _d2_reading(d2r), "",
          "## Rows excluded as held out (not read)", ""]
    for item in p["not_read"]:
        L.append(f"- {item['item']}: {item['why']}.")
    L += ["", "Every row read is DEVELOPMENT in `contrast_split_v1` and has no TEST contrast in `split_v2_per_change` "
          "(laboratories " + ", ".join(sorted({v['lab'] for k, v in p['split_status'].items() if not k.startswith('_')}))
          + "; checked by this generator).", ""]
    return "\n".join(L) + "\n"


def main() -> int:
    t0 = time.time()
    digest = assert_shipped()
    splits = split_status()
    params, drift = core_parameters(SULFUR), core_ph_drift()
    b9_report = json.loads(B9_REPORT.read_text(encoding="utf-8"))
    params_b9, drift_b9 = B27SR._operative_and_drift(b9_report, None)
    G1 = g1(params, drift, params_b9, drift_b9)
    assert_shipped()
    D1 = d1(params, drift)
    assert_shipped()
    D2 = d2(params, drift)
    assert assert_shipped() == digest
    payload = {
        "artifact": "kinetic_core_m1_preflight_checks",
        "provenance": provenance.provenance_block(
            "kinetic_core_m1_preflight_checks", generated_by="scripts/generators/m1_preflight_checks.py",
            wave="M1.0",
            inputs=[B60_REPORT, B9_REPORT, B27_SHIP, S042, CASE, SPLIT_V1, SPLIT_V2,
                    fit_rows.FIT_ROWS_DIR / f"{THIOL_TABLE}.yml", fit_rows.FIT_ROWS_DIR / f"{ROUTE_TABLE}.yml"]),
        "case": data_paths.rel(CASE),
        "b60_report": data_paths.rel(B60_REPORT),
        "config_digest": digest,
        "fitted": False,
        "split_status": splits,
        "not_read": NOT_READ,
        "G1": G1, "D1": D1, "D2": D2,
        "wall_seconds": time.time() - t0,
    }
    artifact_io.write_artifact(payload, OUT, render=render)
    print(f"wrote {data_paths.rel(OUT)} | G1 holds {G1['holds']} (max used {G1['max_consumed_fraction_B60']:.2e}, "
          f"max rise {G1['max_rise_at_phi_1_B60']:.2e}) | D1 MFT 360/30 resid "
          f"{D1['species']['MFT']['ratios']['360_over_30']['residual_dex']:+.2f} dex | {payload['wall_seconds']:.0f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
