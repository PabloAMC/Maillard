"""
CYSTEINE-PENTOSE ROUTE CALIBRATION, stage 1 (2026-10-29): de Roos's route, the TTCA ring opening and cysteine's
own loss, calibrated on the measured time courses of their own species (TTCA, cysteine Amadori compound, free
cysteine, free pentose) in Cao 2017, Zhai 2019 Fig. 3A, Zhai 2021 and Kang 2026 Figs. S3-S4. Flavour levels are
not data here. Everything downstream stays at B60's isothermal optimum. Pre-registration:
results/validation/kinetic_core_cysteine_route_calibration_prereg.md. Writes
results/validation/kinetic_core_cysteine_route_calibration/<arm>[_without_<source>].json.

    python scripts/generators/generate_kinetic_core_cysteine_route_calibration.py --baseline
    python scripts/generators/generate_kinetic_core_cysteine_route_calibration.py --arm extend_last_slope
    python scripts/generators/generate_kinetic_core_cysteine_route_calibration.py --arm flat --leave-out zhai2019
    python scripts/generators/generate_kinetic_core_cysteine_route_calibration.py --checks extend_last_slope
"""
from __future__ import annotations

import argparse
import contextlib
import dataclasses
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for p in (ROOT, ROOT / "scripts" / "generators"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import generate_kinetic_core_b64_fit as B64  # noqa: E402  (inline pH coupling, B60's objective)
import generate_kinetic_core_b59_holdout_diagnosis as HOLDOUT  # noqa: E402
from generate_kinetic_core_b59a_fit import best_member  # noqa: E402
from src import data_paths  # noqa: E402
from src.kinetic_core import fitting  # noqa: E402
from src.kinetic_core import parameters_sulfur as ps  # noqa: E402
from src.kinetic_core import ph_state, sulfur  # noqa: E402

B23 = B64.B23
CELSIUS = 273.15
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_cysteine_route_calibration_prereg.md"
OUT_DIR = data_paths.VALIDATION_DIR / "kinetic_core_cysteine_route_calibration"
ARMS = ("extend_last_slope", "flat")
OX = {"OX": ps.OX_AMBIENT_MMOL_L}
WATER = B23.BUFFER_NONE
CAO_BUFFER = ph_state.BufferSpec(kind="phosphate", phosphate_mol_l=0.2, declared=True, source="Cao 2017 blank pot")
RTOL, ATOL = 1e-7, 1e-12

# ---------------------------------------------------------------------------------------------------------------
# The pots and their rows: (source, pot, species, minutes, measured mmol/L, read-error floor mmol/L)
# ---------------------------------------------------------------------------------------------------------------
POTS = {
    "cao_blank_90C": dict(t_c=90.0, initial={"Cys": 400.0, "PENT": 400.0, **OX}, ph=5.5, buffer=CAO_BUFFER),
    **{f"zhai2019_{t}C": dict(t_c=float(t), initial={"TTCA": 14.5, "Cys": 68.1, "PENT": 150.7, **OX}, ph=7.4,
                              buffer=WATER) for t in (70, 90, 110, 130)},
    **{f"zhai2021_ttca_{t}C": dict(t_c=float(t), initial={"TTCA": 10.0, **OX}, ph=7.0, buffer=WATER)
       for t in (100, 120, 140)},
    **{f"kang_ttca_{t}C": dict(t_c=float(t), initial={"TTCA": 10.0, **OX}, ph=7.0, buffer=WATER) for t in (100, 120, 140)},
    **{f"kang_freecys_{t}C": dict(t_c=float(t), initial={"Cys": 10.0, **OX}, ph=7.0, buffer=WATER) for t in (100, 120)},
}


def _rows():
    rows = []
    cao_t = (30.0, 60.0, 120.0, 180.0, 240.0, 300.0)
    cao = {  # cao2017_extraction.md, Fig. 5 (TTCA, Cys-Amadori; +-5 mM) and Fig. 4 (Cys, Xyl; +-4 mM), mol/L
        "TTCA": ((0.356, 0.338, 0.290, 0.268, 0.254, 0.236), 5.0),
        "CAMA": ((0.006, 0.063, 0.090, 0.106, 0.109, 0.133), 5.0),
        "Cys": ((0.076, 0.068, 0.066, 0.064, 0.062, 0.056), 4.0),
        "PENT": ((0.080, 0.077, 0.073, 0.071, 0.070, 0.069), 4.0),
    }
    for sp, (vals, floor) in cao.items():
        rows += [("cao2017", "cao_blank_90C", sp, t, 1000.0 * v, floor) for t, v in zip(cao_t, vals)]
    zt = (20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0)
    zhai = {  # zhai2019_extraction.md, digitised Fig. 3A (TTCA +-0.7, Cys-ARP +-0.3 mM)
        70: ((41.1, 50.5, 55.6, 54.0, 54.1, 53.1, 54.0, 54.3), (0.2, 0.5, 1.0, 1.3, 1.6, 1.9, 2.1, 2.8)),
        90: ((48.3, 59.2, 57.7, 59.0, 58.3, 59.1, 61.9, 63.2), (0.9, 1.4, 1.7, 1.9, 2.6, 3.0, 3.4, 3.7)),
        110: ((58.5, 55.9, 53.9, 52.0, 52.8, 48.7, 42.6, 39.0), (2.3, 3.7, 5.6, 5.3, 4.8, 4.6, 4.2, 3.8)),
        130: ((45.7, 42.5, 34.5, 32.2, 21.2, 15.8, 11.3, 6.9), (4.3, 4.8, 3.9, 3.1, 1.6, 1.2, 0.0, 0.0)),
    }
    for t_c, (ttca, cama) in zhai.items():
        rows += [("zhai2019", f"zhai2019_{t_c}C", "TTCA", t, v, 0.7) for t, v in zip(zt, ttca)]
        rows += [("zhai2019", f"zhai2019_{t_c}C", "CAMA", t, v, 0.3) for t, v in zip(zt, cama)]
    for t_c, (c0, k) in {100: (10.331, 0.0271), 120: (9.9718, 0.0651), 140: (9.3375, 0.0813)}.items():
        rows.append(("zhai2021", f"zhai2021_ttca_{t_c}C", "TTCA", 60.0, c0 - k * 60.0, 0.0))  # zhai2021_extraction.md
    kt = (20.0, 40.0, 60.0, 80.0, 100.0, 120.0)
    s3 = {100: (0.121, 0.242, 0.693, 0.872, 1.30, 1.233), 120: (0.205, 0.989, 1.219, 1.112, 0.988, 0.874),
          140: (0.500, 1.630, 1.574, 1.619, 1.30, 1.088)}  # kang2026_SI_extraction.md sec. 6c, +-0.03 mM
    for t_c, vals in s3.items():
        rows += [("kang_s3", f"kang_ttca_{t_c}C", "Cys", t, v, 0.03) for t, v in zip(kt, vals)]
    s4 = {100: (9.77, 9.61, 9.41, 9.12, 8.64, 8.38), 120: (8.63, 7.81, 7.11, 6.90, 6.55, 6.13)}  # sec. 6b, +-0.03
    for t_c, vals in s4.items():
        rows += [("kang_s4", f"kang_freecys_{t_c}C", "Cys", t, v, 0.03) for t, v in zip(kt, vals)]
    return tuple(dict(source=s, pot=p, species=sp, minutes=t, measured=m, floor=f) for s, p, sp, t, m, f in rows)


ROWS = _rows()
SOURCES = tuple(sorted({r["source"] for r in ROWS}))

# ---------------------------------------------------------------------------------------------------------------
# The 8 free coordinates (prereg sec. 2)
# ---------------------------------------------------------------------------------------------------------------
_X = HOLDOUT._vector(best_member("isothermal"))
_B60_PARAMS = B23.build_parameters(_X)
_DRIFT = B23.unpack(_X)[3]
_ROUTE0 = dataclasses.replace(ps.CYS_AMADORI_ROUTE, enabled=True)
_TTCA0 = ps.TTCA_FORMATION
KEYS = ("log10_k_ttca_cys_145C", "ea_ttca_ring_kj_mol", "log10_k_ttca_deg_145C", "log10_K_ttca_l_per_mmol",
        "log10_k_cama_form_90C", "ea_cama_form_kj_mol", "ea_cama_decay_kj_mol", "log10_k_cys_thermal_145C")
START = np.array([
    math.log10(_B60_PARAMS["k_ttca_cys"].k_ref), _B60_PARAMS["k_ttca_cys"].ea_kj_mol,
    math.log10(_B60_PARAMS["k_ttca_deg"].k_ref), math.log10(_TTCA0.k_l_per_mmol),
    math.log10(_ROUTE0.k_form_ref_per_min), _ROUTE0.ea_form_kj_mol, _ROUTE0.ea_decay_kj_mol,
    math.log10(_B60_PARAMS["k_cys_thermal"].k_ref)])
EA_BAND = ps.FORMATION_EA_BOUNDS_BY_ROUTE["sugar_trunk"]
LOWER = np.array([START[0] - 3, EA_BAND[0], START[2] - 3, math.log10(0.019), START[4] - 3, EA_BAND[0], EA_BAND[0],
                  START[7] - 3])
UPPER = np.array([START[0] + 3, EA_BAND[1], START[2] + 3, math.log10(0.075), START[4] + 3, EA_BAND[1], EA_BAND[1],
                  START[7] + 3])


@contextlib.contextmanager
def installed(theta: np.ndarray, arm: str):
    """The route on with ``theta``'s constants, restored afterwards (module switches are process-global)."""
    saved = (ps.CYS_AMADORI_ROUTE, ps.TTCA_FORMATION)
    ps.CYS_AMADORI_ROUTE = dataclasses.replace(_ROUTE0, k_form_ref_per_min=10.0 ** theta[4], ea_form_kj_mol=theta[5],
                                               ea_decay_kj_mol=theta[6], above_measured_ph=arm)
    ps.TTCA_FORMATION = dataclasses.replace(_TTCA0, k_l_per_mmol=10.0 ** theta[3])
    try:
        yield
    finally:
        ps.CYS_AMADORI_ROUTE, ps.TTCA_FORMATION = saved


def overridden(parameters: dict, theta: np.ndarray) -> dict:
    out = dict(parameters)
    out["k_ttca_cys"] = dataclasses.replace(out["k_ttca_cys"], k_ref=10.0 ** theta[0], ea_kj_mol=theta[1])
    out["k_ttca_deg"] = dataclasses.replace(out["k_ttca_deg"], k_ref=10.0 ** theta[2], ea_kj_mol=theta[1])
    out["k_cys_thermal"] = dataclasses.replace(out["k_cys_thermal"], k_ref=10.0 ** theta[7])
    return out


def run_pot(name: str, parameters: dict, times) -> "sulfur.SulfurRun":
    spec = POTS[name]
    return sulfur.integrate_sulfur(parameters, spec["t_c"] + CELSIUS, spec["initial"], [0.0] + sorted(times),
                                   ph=spec["ph"], buffer_spec=spec["buffer"], ph_drift=_DRIFT, rtol=RTOL, atol=ATOL)


def predictions(theta: np.ndarray, arm: str, rows=ROWS) -> np.ndarray:
    parameters = overridden(_B60_PARAMS, theta)
    by_pot = {}
    for r in rows:
        by_pot.setdefault(r["pot"], set()).add(r["minutes"])
    with installed(theta, arm):
        runs = {}
        for name, times in by_pot.items():
            run = run_pot(name, parameters, times)
            grid = [0.0] + sorted(times)
            runs[name] = {s: dict(zip(grid, run.series(s))) for s in ("TTCA", "CAMA", "Cys", "PENT")}
    return np.array([float(runs[r["pot"]][r["species"]][r["minutes"]]) for r in rows])


def residuals(theta: np.ndarray, arm: str, rows=ROWS) -> np.ndarray:
    """log10 (pred + floor) / (meas + floor), in dex; a failed integration is handled by ``fitting.Problem``."""
    pred = predictions(theta, arm, rows)
    meas = np.array([r["measured"] for r in rows])
    floor = np.array([r["floor"] for r in rows])
    return np.log10((np.maximum(pred, 0.0) + floor + 1e-12) / (meas + floor + 1e-12))


def by_source(r: np.ndarray, rows=ROWS) -> dict:
    out = {}
    for s in sorted({row["source"] for row in rows}):
        idx = [i for i, row in enumerate(rows) if row["source"] == s]
        out[s] = {"n": len(idx), "sum_sq": float(np.sum(r[idx] ** 2)), "rms_dex": float(np.sqrt(np.mean(r[idx] ** 2)))}
    return out


def problem(arm: str, upper: np.ndarray = UPPER) -> fitting.Problem:
    """The calibration as a ``fitting.Problem`` (2026-10-29: the shared module replaced this generator's own
    optimiser call and identification code, both of which had been wrong; run log 4)."""
    return fitting.Problem(keys=KEYS, kinds=tuple("Ea" if k.startswith("ea_") else "log10k" for k in KEYS),
                           lower=LOWER, upper=np.asarray(upper, dtype=float), rows=ROWS,
                           residuals=lambda theta, rows: residuals(theta, arm, rows),
                           group_of=lambda row: row["source"], failed_residual=5.0)


def member_record(prob: fitting.Problem, f_: fitting.Fit, lap: dict, arm: str, leave_out: str | None) -> dict:
    rows = tuple(r for r in ROWS if r["source"] != leave_out)
    return {"arm": arm, "leave_out": leave_out, "n_rows": len(rows), "keys": KEYS, "fitting": "src/kinetic_core/fitting.py",
            "start": f_.start.tolist(), "theta": f_.theta.tolist(), "start_cost": f_.start_cost, "cost": f_.cost,
            "status": f_.status, "budget_exhausted": f_.budget_exhausted, "evaluations": f_.evaluations,
            "seconds": f_.seconds, "by_source": by_source(f_.residuals, rows), "laplace": lap,
            "lower": prob.lower.tolist(), "upper": prob.upper.tolist(),
            "summary": fitting.summary_lines(prob, f_, lap),
            "residuals": [dict(row, residual_dex=float(v)) for row, v in zip(rows, f_.residuals)]}


# ---------------------------------------------------------------------------------------------------------------
# Pre-declared checks (prereg sec. 4), not fitted
# ---------------------------------------------------------------------------------------------------------------
def checks(theta: np.ndarray, arm: str) -> dict:
    parameters = overridden(_B60_PARAMS, theta)
    out = {}
    with installed(theta, arm):
        z = sulfur.integrate_sulfur(parameters, 120.0 + CELSIUS, {"PENT": 40.0, "Cys": 40.0, **OX}, [0, 10, 20, 30, 50],
                                    ph=6.5, buffer_spec=HOLDOUT.ZHANG_BUFFER, ph_drift=_DRIFT, rtol=RTOL, atol=ATOL)
        out["zhang2025b_log10_pred_over_meas"] = {
            s: {str(t): math.log10(max(float(z.series(s)[i + 1]), 1e-12) / m) for i, (t, m) in enumerate(HOLDOUT.ZHANG_MM[s].items())}
            for s in ("TTCA", "CAMA")}
        d = sulfur.integrate_sulfur(parameters, 120.0 + CELSIUS, {"TTCA": 5.0, **OX}, [0, 10, 20, 30, 40, 50, 60], ph=7.4,
                                    buffer_spec=WATER, ph_drift=_DRIFT, rtol=RTOL, atol=ATOL)
        meas_um = {"DPO": (42, 68, 76, 89, 91, 98), "TDP": (89, 164, 199, 205, 228, 249)}  # 1-DX, 3-DX (Zhai 2019 Fig. 3B)
        out["zhai2019_fig3b_log10_pred_over_meas"] = {
            s: [math.log10(max(1000.0 * float(d.series(s)[i + 1]), 1e-9) / m) for i, m in enumerate(v)] for s, v in meas_um.items()}
        out["zhai2019_fig3b_3dx_over_1dx_pred_60min"] = float(d.series("TDP")[-1] / max(d.series("DPO")[-1], 1e-15))
        cao_ph = {}
        for t, m in zip((30, 60, 120, 180, 240, 300), (4.67, 4.56, 4.49, 4.33, 4.26, 4.15)):
            run = run_pot("cao_blank_90C", parameters, [float(t)])
            cao_ph[str(t)] = {"pred_cooled": float(run.metadata["ph_final_cooled"]), "measured": m}
        out["cao_ph_trajectory"] = cao_ph
        kang_ph = {}
        for ph0, m in ((5.5, 5.1), (7.0, 4.9), (8.0, 4.7)):
            run = sulfur.integrate_sulfur(parameters, 120.0 + CELSIUS, {"TTCA": 10.0, **OX}, [0.0, 120.0], ph=ph0,
                                          buffer_spec=WATER, ph_drift=_DRIFT, rtol=RTOL, atol=ATOL)
            kang_ph[str(ph0)] = {"pred_cooled": float(run.metadata["ph_final_cooled"]), "measured": m}
        out["kang_final_ph"] = kang_ph
        cao = run_pot("cao_blank_90C", parameters, [300.0])
        pool = [float(cao.series("Cys")[i] + cao.series("TTCA")[i] + cao.series("CAMA")[i]) for i in (0, 1)]
        out["zhao2017_cysteine_pool_loss_per_min"] = -math.log(pool[1] / pool[0]) / 300.0
    out["flavour_level_diagnosis"] = flavour_diagnosis(theta, arm)
    return out


def flavour_diagnosis(theta: np.ndarray, arm: str) -> dict:
    """B60's 52 rows and the B59 hold-outs at the stage-1 constants (everything downstream at B60)."""
    original = B23.build_parameters
    B23.build_parameters = lambda x: overridden(original(x), theta)
    try:
        with installed(theta, arm):
            B64.install_reading("isothermal")
            r = B64.B60.residual_vector(_X, True)
            out = {"b60_rows_cost": float(0.5 * np.sum(r ** 2)),
                   "b60_rows_worst": sorted(((row["id"], float(v)) for row, v in zip(B23.ACTIVE_FIT_ROWS, r)),
                                            key=lambda kv: -abs(kv[1]))[:8]}
            p = B23.build_parameters(_X)
            hot = {}
            for reading in ("isothermal", "ramp_plus_hold", "ramp_only_newtonian_tau4", "ramp_only_newtonian_tau6.7", "ramp_only"):
                prof = HOLDOUT._hot_profile(reading)
                hot[reading] = HOLDOUT._run(p, _DRIFT, 145.0, HOLDOUT.HOFMANN, [float(prof[0][-1])], 5.0,
                                            HOLDOUT.HOFMANN_BUFFER, profile=prof).final("MFT")
            iso = HOLDOUT.diagnose(_X, "isothermal")
            cool = iso["schieberle_mft_100C_6h_over_145C"] * hot["isothermal"]
            out.update({"schieberle_factor_by_145C_reading": {k: cool / v for k, v in hot.items()},
                        "schieberle_shape_over_30min": iso["schieberle_mft_shape_over_30min"],
                        "nestle_shape_over_4h": iso["nestle_mft_shape_over_4h"],
                        "hofmann_ph_sweep_log10": iso["hofmann_pH_sweep_log10_pred_over_meas"]})
    finally:
        B23.build_parameters = original
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm", choices=ARMS)
    parser.add_argument("--leave-out", dest="leave_out", choices=SOURCES)
    parser.add_argument("--budget", type=int, default=1200, help="objective evaluations, Jacobian columns included")
    parser.add_argument("--baseline", action="store_true", help="V3: the rows at the start vector, both arms, no fit")
    parser.add_argument("--checks", choices=ARMS, help="the pre-declared checks at that arm's full fit")
    parser.add_argument("--leave-one-out", dest="leave_one_out",
                        help="leave-one-source-out from this stored full fit (a file stem), with the predictive check")
    parser.add_argument("--tag", default="", help="suffix for the output file")
    parser.add_argument("--reflag", action="store_true", help="re-apply the bound / identification rule to stored fits")
    parser.add_argument("--k-upper", dest="k_upper", type=float,
                        help="sensitivity: K's upper bound in L/mmol (run log 2), instead of Cao's band top 0.075")
    parser.add_argument("--start-from", dest="start_from",
                        help="second start: the theta of this member (a file stem in the output folder)")
    args = parser.parse_args(argv)
    assert PREREG.exists(), "pre-register before running"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if args.baseline:
        report = {}
        for arm in ARMS:
            r = residuals(START, arm)
            report[arm] = {"cost": float(0.5 * np.sum(r ** 2)), "by_source": by_source(r),
                           "residuals": [dict(row, residual_dex=float(v)) for row, v in zip(ROWS, r)]}
            print(arm, "start cost", round(report[arm]["cost"], 2), {s: round(v["rms_dex"], 2) for s, v in report[arm]["by_source"].items()})
        (OUT_DIR / "baseline_start_vector.json").write_text(json.dumps(report, indent=1))
        return 0
    if args.reflag:
        for path in sorted(OUT_DIR.glob("*.json")):
            member = json.loads(path.read_text())
            if "laplace" not in member:
                continue
            upper = np.array(member.get("upper", UPPER))
            if "k_upper" in path.stem and "upper" not in member:
                upper = UPPER.copy()
                upper[3] = math.log10(float(path.stem.split("k_upper_")[1]))
            sig = [member["laplace"]["sigma"][k] for k in KEYS]
            member["laplace"].update(fitting.identification(problem(member["arm"], upper), member["theta"], sig))
            path.write_text(json.dumps(member, indent=1, default=str))
            print(path.stem, "on bound:", [k for k, v in member["laplace"]["at_bound"].items() if v],
                  "unidentified:", [k for k, v in member["laplace"]["identified"].items() if not v])
        return 0
    if args.leave_one_out:
        path = OUT_DIR / f"{args.leave_one_out}.json"
        member = json.loads(path.read_text())
        prob = problem(member["arm"], np.array(member["upper"]))
        out = fitting.leave_one_group_out(prob, np.array(member["theta"]), member["laplace"], args.budget)
        report = {g: {**{k: v for k, v in r.items() if k != "fit"}, "theta": r["fit"].theta.tolist(),
                      "cost": r["fit"].cost, "status": r["fit"].status} for g, r in out.items()}
        (OUT_DIR / f"{args.leave_one_out}_leave_one_source_out.json").write_text(json.dumps(report, indent=1, default=str))
        for g, r in report.items():
            pr = r["predictive"]
            print(f"without {g:10s} beyond 2 sigma {r['beyond_2_sigma']} on bound {r['on_bound']}  held-out rms "
                  f"{pr['rms_held_out']:.3f} vs in-sample {pr['rms_in_sample']:.3f}  F {pr['F']:.2f} p {pr['p']:.3g}"
                  f"{'  TENSION' if pr['tension'] else ''}")
        return 0
    if args.checks:
        path = OUT_DIR / f"{args.checks}.json"
        member = json.loads(path.read_text())
        member["checks"] = checks(np.array(member["theta"]), args.checks)
        path.write_text(json.dumps(member, indent=1, default=str))
        print(json.dumps({k: v for k, v in member["checks"].items() if k != "flavour_level_diagnosis"}, indent=1)[:3000])
        print(json.dumps(member["checks"]["flavour_level_diagnosis"], indent=1, default=str)[:2500])
        return 0
    start = START if args.start_from is None else np.array(json.loads((OUT_DIR / f"{args.start_from}.json").read_text())["theta"])
    upper = UPPER.copy()
    if args.k_upper:
        upper[3] = math.log10(args.k_upper)
    prob = problem(args.arm, upper)
    rows = tuple(r for r in ROWS if r["source"] != args.leave_out)
    f_ = fitting.fit(prob, start, args.budget, rows=rows)
    lap = fitting.laplace(prob, f_.theta, rows=rows)
    member = member_record(prob, f_, lap, args.arm, args.leave_out)
    tag = (args.arm + (f"_without_{args.leave_out}" if args.leave_out else "")
           + ("_second_start" if args.start_from else "") + (f"_k_upper_{args.k_upper:g}" if args.k_upper else "")
           + (f"_{args.tag}" if args.tag else ""))
    (OUT_DIR / f"{tag}.json").write_text(json.dumps(member, indent=1, default=str))
    print(tag)
    print("\n".join(member["summary"]))
    print("  by source", {s_: round(v["rms_dex"], 3) for s_, v in member["by_source"].items()})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
