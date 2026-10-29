#!/usr/bin/env python
"""
fit_ph_acid_submodel.py -- sub-model 1's fit: pH and acid production (results/validation/kinetic_core_ph_acid_submodel_prereg.md
sec. 4a and sec. 5, the design review's fixes; 2026-10-30).

The probe found the sulfur lane's acid supply blind to pH (the trunk's measured pH terms never applied there) and a
fitted acetic-acid sink that two data sets contradict and one requires. With ``sulfur.TRUNK_PH_TERMS`` and
``ph_state.TRUNK_SOLUTES_IN_BALANCE`` on, this fits:

  * the pH slopes (decades per pH unit) of the two observed Amadori-decay steps, k_ama_tdg and k_ama_odg (k_ama_mgo,
    which no row here sees, tied to k_ama_odg), and of the 3-DG -> formic exit k_tdg_fa. FLAT in [0, 1.5]: Martins
    2003 Table 3, the source of the shipped slopes, is fitted to the Martins 2003b rows in this table (review B1);
  * a log10 offset of each of those steps' tabulated constant (k_ama_mgo tied to k_ama_odg): Martins 2005's
    constants are averages over the drift from pH 6.8 (Martins 2005, p. 267: the pH effect "is concealed in the rate
    constants"), so with the slopes applied along the drift the constant at the reference pH is a free quantity,
    not the tabulated one (review B2);
  * the pentose acid yield per sink event and the ARP secondary-ammonium pKa (B2.2's two constants; no prior);
  * arm A only: log10 k_ref and Ea of ``k_aa_frag`` (review S4: the Ea free, the data span 80-120 C). Arm B: the
    step off.

Rows: the READY table less the casein systems (Brands 2001, Nguyen 2016), Davidek's NaOH uptake (it counts acids the
lane cannot make: NaOH / (FA+AA) 1.1-2.2, review B5), Zhao 2017b's aldehyde pots (the aldehyde binds cysteine, they
are not blanks: B6), pots the lane cannot run, and t = 0 pH rows. Errors: each row's own plus an extra noise per
data set x observable, estimated from the scatter WITHIN each series around its own mean (a series' bias is misfit,
not noise: review B7) and re-estimated inside each leave-one-laboratory-out fold on the kept rows only. Two starts
(S6). Done-tests on RAW residuals against tolerances declared in the prereg before this ran.

The pH coupling is ``inline`` (checked against the shipped fixed point on this table, with these switches, before
the fit: prereg sec. 5); every Jacobian point's pots run in one pool (``fitting.Problem.evaluate_many``).

    python scripts/generators/fit_ph_acid_submodel.py --arm A|B [--budget N] [--rounds N] [--lolo]
Writes results/validation/kinetic_core_ph_acid_submodel_fit_<arm>.{json,md}.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for entry in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

from src import data_paths, provenance  # noqa: E402
from src.kinetic_core import engine, fit_rows, fitting, parallel, ph_state, sulfur  # noqa: E402
from src.kinetic_core.ph_state import ACID_YIELD_BOUNDS, ARP_AMINE_PKA_BOUNDS, PhDrift  # noqa: E402

import build_submodel_probe as probe  # noqa: E402

TABLE = "ph_acid_submodel"
#: diagnostic only (prereg sec. 4a finding 5; review B5, B6)
DIAGNOSTIC_DATA_SETS = ("brands2001_fig3_ph_titratable_acid", "nguyen2016_figs1_ph_mass_balance")
DIAGNOSTIC_OBSERVABLES = ("NaOH_uptake",)
NOT_BLANKS = ("zhao2017b_hexanal", "zhao2017b_decadienal")
#: sloped steps: coordinate suffix -> the parameter keys it moves
SLOPED = {"ama_tdg": ("k_ama_tdg",), "ama_odg": ("k_ama_odg", "k_ama_mgo"), "tdg_fa": ("k_tdg_fa",)}
#: done-test tolerances on raw residuals (prereg sec. 5, declared before the fit): RMS per data set x observable
TOLERANCE = {"pH": 0.30, "pH_6x_dilution_25C": 0.30, "amount": 0.30}
#: the second start (S6): the other side of every band from the first
START_B = {"slope_": 0.40, "offset_": -0.30, "acid_yield": 0.30, "arp_amine_pka": 8.0, "log10_k_aa_frag": -3.0,
           "ea_k_aa_frag": 80.0}


def coordinates(arm: str):
    drift = engine.core_ph_drift()
    keys, kinds, lower, upper, start = [], [], [], [], []
    for name in SLOPED:
        keys += [f"slope_{name}", f"offset_{name}"]
        kinds += ["ph_slope", "log10_offset"]
        lower += [0.0, -1.5]
        upper += [1.5, 1.5]
        first = sulfur._trunk_ph_exponents()[SLOPED[name][0]]
        start += [first, 0.0]
    keys += ["acid_yield", "arp_amine_pka"]
    kinds += ["acid_yield", "pka"]
    lower += [ACID_YIELD_BOUNDS[0], ARP_AMINE_PKA_BOUNDS[0]]
    upper += [ACID_YIELD_BOUNDS[1], ARP_AMINE_PKA_BOUNDS[1]]
    start += [drift.acid_yield, drift.arp_amine_pka]
    if arm == "A":
        k_aa = engine.core_parameters(engine.SULFUR)["k_aa_frag"]
        keys += ["log10_k_aa_frag", "ea_k_aa_frag"]
        kinds += ["log10k", "Ea"]
        lower += [-6.0, 20.0]
        upper += [0.0, 160.0]
        start += [math.log10(k_aa.k_ref), max(float(k_aa.ea_kj_mol), 20.0)]
    return tuple(keys), tuple(kinds), np.array(lower), np.array(upper), np.array(start, dtype=float)


def second_start(keys, lower, upper) -> np.ndarray:
    out = []
    for k, lo, hi in zip(keys, lower, upper):
        v = next((val for prefix, val in START_B.items() if k.startswith(prefix)), None)
        out.append(min(max(v, lo), hi))
    return np.array(out, dtype=float)


def trial(theta: Sequence[float], keys) -> Dict[str, Any]:
    """The overrides one trial vector makes: slopes for the lane, k_ref values and the drift for the pot runner."""
    t = dict(zip(keys, map(float, theta)))
    base = engine.core_parameters(engine.SULFUR)
    slopes, k_ref = {}, {}
    for name, steps in SLOPED.items():
        for key in steps:
            slopes[key] = t[f"slope_{name}"]
            k_ref[key] = base[key].k_ref * 10.0 ** t[f"offset_{name}"]
    ea = {}
    if "log10_k_aa_frag" in t:
        k_ref["k_aa_frag"] = 10.0 ** t["log10_k_aa_frag"]
        ea["k_aa_frag"] = t["ea_k_aa_frag"]
    return {"slopes": slopes, "k_ref": k_ref, "ea": ea,
            "without": () if "log10_k_aa_frag" in t else ("k_aa_frag",),
            "drift": PhDrift(acid_yield=t["acid_yield"], arp_amine_pka=t["arp_amine_pka"])}


def run_task(job: Dict[str, Any]):
    """One pot at one trial vector, in whichever process: the lane's switches and slopes are set here."""
    sulfur.TRUNK_PH_TERMS = True
    ph_state.TRUNK_SOLUTES_IN_BALANCE = True
    sulfur.TRUNK_PH_EXPONENT_OVERRIDE = dict(job["slopes"])
    try:
        return probe.run_pot(job)
    except Exception as exc:  # noqa: BLE001 - a failed pot is counted and its rows become failed residuals
        return {"failed": f"{type(exc).__name__}: {exc}"}


FAILURES: List[str] = []


def make_raw_many(jobs_by_pot: Dict[str, dict], keys):
    """raw_many(points, rows): each point's residuals in the rows' own units, every (point, pot) in one pool."""

    def raw_many(points, rows):
        pots = sorted({r.pot for r in rows})
        trials = [trial(x, keys) for x in points]
        tasks = [dict(jobs_by_pot[p], **tr) for tr in trials for p in pots]
        results = parallel.map_ordered(run_task, tasks)
        out = []
        for i, tr in enumerate(trials):
            res = dict(zip(pots, results[i * len(pots):(i + 1) * len(pots)]))
            jobs = {p: dict(jobs_by_pot[p], **tr) for p in pots}
            r_i = []
            for r in rows:
                if "failed" in res[r.pot]:
                    FAILURES.append(f"{r.pot}: {res[r.pot]['failed']}")
                    r_i.append(float("nan"))
                else:
                    r_i.append(fit_rows.residual(r, probe.predict(r, jobs[r.pot], res[r.pot])))
            out.append(np.array(r_i, dtype=float))
        return out

    return raw_many


def pooled_noise(model, raw_residuals, rows, n_constants=0) -> Dict[str, float]:
    """Extra noise per data set x observable from the scatter WITHIN each series around its own mean (review B7: a
    series' bias is model misfit and must stay in the cost, not become noise), assigned to every series of the group."""
    r = np.asarray(raw_residuals, dtype=float)
    groups: Dict[tuple, Dict[str, List[int]]] = {}
    for i, row in enumerate(rows):
        groups.setdefault((row.ledger_id, row.species), {}).setdefault(model.series_of(row), []).append(i)
    out: Dict[str, float] = {}
    n_total = max(len(rows), 1)
    for series in groups.values():
        idx_all = [i for idx in series.values() for i in idx]
        ss = sum(float(np.sum((r[idx] - np.mean(r[idx])) ** 2)) for idx in series.values())
        dof = len(idx_all) - len(series) - n_constants * len(idx_all) / n_total
        if dof < 2:
            continue
        s2 = ss / dof - float(np.mean([model.sigma_of(rows[i]) ** 2 for i in idx_all]))
        if s2 > 0:
            for name in series:
                out[name] = math.sqrt(s2)
    return out


def attributes(pots):
    def sugar(r):
        charge = pots[r.pot]["initial_mmol_L"]
        return next((k for k in ("Glc", "Fru", "AMA", "PENT", "ARP") if charge.get(k)), "-")

    def amine(r):
        charge = pots[r.pot]["initial_mmol_L"]
        return "+".join(k for k in ("Gly", "Cys") if charge.get(k)) or "none"

    def buffer(r):
        return probe.buffer_of(pots[r.pot], float(pots[r.pot]["label_ph"])).kind

    return {"laboratory": lambda r: r.laboratory, "data_set": lambda r: r.ledger_id, "observable": lambda r: r.species,
            "temperature_C": lambda r: float(pots[r.pot]["t_c"]), "label_pH": lambda r: float(pots[r.pot]["label_ph"]),
            "time_min": lambda r: r.t_min, "sugar": sugar, "amine": amine, "buffer": buffer}


def raw_table(rows, raw) -> Dict[str, Dict[str, Any]]:
    """Per data set x observable: n, mean and RMS of the raw residual, and the declared tolerance's verdict."""
    groups: Dict[str, List[float]] = {}
    for row, v in zip(rows, raw):
        groups.setdefault(f"{row.ledger_id} | {row.species}", []).append(float(v))
    out = {}
    for key, vals in sorted(groups.items()):
        species = key.split(" | ")[1]
        tol = TOLERANCE.get(species, TOLERANCE["amount"])
        rms = float(np.sqrt(np.mean(np.square(vals))))
        out[key] = {"n": len(vals), "mean": float(np.mean(vals)), "rms": rms, "tolerance": tol, "within": rms <= tol}
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--arm", choices=("A", "B"), required=True)
    parser.add_argument("--budget", type=int, default=480)
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--lolo", action="store_true")
    parser.add_argument("--ph-coupling", choices=("fixed_point_nodes", "inline"), default="inline")
    args = parser.parse_args(argv)
    t0 = time.time()
    sulfur.PH_COUPLING = args.ph_coupling
    sulfur.TRUNK_PH_TERMS = True
    ph_state.TRUNK_SOLUTES_IN_BALANCE = True
    pots, all_rows = fit_rows.load(TABLE)
    jobs = {j["pot"]: j for j in probe.plan(pots, all_rows) if "skip" not in j and j["pot"] not in NOT_BLANKS}
    data = tuple(r for r in all_rows
                 if r.pot in jobs and r.ledger_id not in DIAGNOSTIC_DATA_SETS and r.species not in DIAGNOSTIC_OBSERVABLES
                 and not (r.species == "pH" and r.t_min == 0.0))
    diagnostic = tuple(r for r in all_rows if r.pot in jobs and r not in data and not (r.species == "pH" and r.t_min == 0.0))
    keys, kinds, lower, upper, start = coordinates(args.arm)
    raw_many = make_raw_many(jobs, keys)

    def raw(theta, rows):
        return raw_many([theta], rows)[0]

    model = fit_rows.measurement_model()

    def make_problem(m, rows=data):
        return fitting.Problem(keys=keys, kinds=kinds, lower=lower, upper=upper, rows=rows,
                               residuals=lambda th, rs: m.whiten(raw(th, rs), rs),
                               evaluate_many=lambda pts, rs: [m.whiten(v, rs) for v in raw_many(pts, rs)],
                               group_of=lambda r: r.laboratory, nuisance_dof=m.nuisance_dof)

    raw_start = raw(start, data)
    best, final_model = fitting.fit_with_series_noise(make_problem, model, raw, start, budget=args.budget,
                                                      rounds=args.rounds, estimator=pooled_noise)
    problem = make_problem(final_model)
    starts = fitting.two_starts(problem, best.theta, second_start(keys, lower, upper), budget=args.budget)
    best = starts["best"]
    lap = starts["laplace"]
    raw_best = raw(best.theta, data)
    attr = attributes(pots)
    report: Dict[str, Any] = {
        "artifact": f"kinetic_core_ph_acid_submodel_fit_{args.arm}",
        "provenance": provenance.provenance_block(
            f"kinetic_core_ph_acid_submodel_fit_{args.arm}", generated_by="scripts/generators/fit_ph_acid_submodel.py",
            inputs=[fit_rows.FIT_ROWS_DIR / f"{TABLE}.yml", ROOT / "src/kinetic_core/sulfur.py",
                    ROOT / "src/kinetic_core/ph_state.py", ROOT / "src/kinetic_core/observations.py",
                    ROOT / "src/kinetic_core/trunk_conditions.py", ROOT / "src/kinetic_core/parameters.py",
                    ROOT / "src/kinetic_core/fitting.py", ROOT / "scripts/generators/build_submodel_probe.py"]),
        "arm": args.arm, "ph_coupling": args.ph_coupling, "rows": len(data), "pots": len({r.pot for r in data}),
        "start": dict(zip(keys, map(float, start))),
        "fit": best.as_dict(keys),
        "two_starts": {"converged": starts["converged"], "cost_agreement": starts["cost_agreement"],
                       "disagreeing": starts["disagreeing"], "unidentified": starts["unidentified"],
                       "costs": [f.cost for f in starts["fits"]],
                       "thetas": [dict(zip(keys, map(float, f.theta))) for f in starts["fits"]]},
        "laplace": lap,
        "series_extra_noise": dict(final_model.extra_noise),
        "failed_pot_runs": len(FAILURES), "failures_sample": FAILURES[:10],
        "raw_at_start": raw_table(data, raw_start),
        "raw_at_optimum": raw_table(data, raw_best),
        "diagnostic_rows_at_optimum": raw_table(diagnostic, raw(best.theta, diagnostic)) if diagnostic else {},
        "diagnostics_whitened": fitting.diagnostics(problem, best.theta, attr),
        "diagnostics_raw": fitting.diagnostics(
            fitting.Problem(keys=keys, kinds=kinds, lower=lower, upper=upper, rows=data,
                            residuals=lambda th, rs: raw(th, rs)), best.theta, attr),
    }
    if args.ph_coupling != "fixed_point_nodes":
        sulfur.PH_COUPLING = "fixed_point_nodes"
        shipped = raw(best.theta, data)
        sulfur.PH_COUPLING = args.ph_coupling
        report["shipped_coupling_max_abs_difference"] = float(np.nanmax(np.abs(shipped - raw_best)))
    if args.lolo:
        report["leave_one_laboratory_out"] = leave_one_lab_out(make_problem, model, raw, data, best.theta, keys, lap,
                                                               args.budget)
    report["wall_seconds"] = round(time.time() - t0, 1)
    out = data_paths.VALIDATION_DIR / f"kinetic_core_ph_acid_submodel_fit_{args.arm}.json"
    out.write_text(json.dumps(report, indent=1, default=float) + "\n", encoding="utf-8")
    out.with_suffix(".md").write_text(markdown(report, problem, best, lap), encoding="utf-8")
    print(f"arm {args.arm}: cost {best.start_cost:.4g} -> {best.cost:.4g}; two starts converged {starts['converged']}; "
          f"{len(FAILURES)} failed pot runs; {report['wall_seconds']} s -> {data_paths.rel(out.with_suffix('.md'))}")
    return 0


def leave_one_lab_out(make_problem, model, raw, data, theta, keys, lap, budget) -> Dict[str, Any]:
    """Each laboratory left out: the noise re-estimated and the constants refitted on the kept rows only (review B7);
    the held-out rows' RAW residuals per data set against the declared tolerances; constant shifts in full-fit sigmas."""
    out = {}
    sig = np.array([lap["sigma"][k] for k in keys])
    for lab in sorted({r.laboratory for r in data}):
        kept = tuple(r for r in data if r.laboratory != lab)
        held = tuple(r for r in data if r.laboratory == lab)
        f_, m = fitting.fit_with_series_noise(lambda mm: make_problem(mm, kept), model, raw, theta,
                                              budget=budget // 2, rounds=2, estimator=pooled_noise)
        shift = (f_.theta - np.asarray(theta)) / sig
        out[lab] = {"rows_kept": len(kept), "rows_held": len(held), "theta": dict(zip(keys, map(float, f_.theta))),
                    "shift_sigma": dict(zip(keys, map(float, shift))),
                    "beyond_2_sigma": [k for k, s in zip(keys, shift) if np.isfinite(s) and abs(s) > 2.0],
                    "held_out_raw": raw_table(held, raw(f_.theta, held))}
    return out


def markdown(report: Dict[str, Any], problem, best, lap) -> str:
    arm = "k_aa_frag free (rate and Ea)" if report["arm"] == "A" else "k_aa_frag off"
    lines = [f"# Sub-model 1 fit, arm {report['arm']} ({arm})", "",
             f"*Generated by `scripts/generators/fit_ph_acid_submodel.py --arm {report['arm']}`; pre-registration "
             "`kinetic_core_ph_acid_submodel_prereg.md` sec. 4a and 5 (the design review's fixes). "
             f"{report['rows']} rows on {report['pots']} pots; pH coupling `{report['ph_coupling']}`. Raw residual = "
             "predicted minus measured (pH units; dex for amounts).*", "", "## Result", ""]
    lines += [f"- {line}" for line in fitting.summary_lines(problem, best, lap)]
    ts = report["two_starts"]
    lines += [f"- two starts: converged {ts['converged']} (costs {', '.join(f'{c:.4g}' for c in ts['costs'])}; "
              f"disagreeing {ts['disagreeing'] or 'none'})",
              f"- failed pot runs: {report['failed_pot_runs']}"]
    if "shipped_coupling_max_abs_difference" in report:
        lines.append(f"- largest raw-residual difference, inline vs the shipped fixed point, at the optimum: "
                     f"{report['shipped_coupling_max_abs_difference']:.2g}")
    lines += ["", "## Done-test 1 on raw residuals (tolerances declared before the fit)", "",
              "| data set \\| observable | n | start mean / RMS | optimum mean / RMS | tolerance | within |", "|---|---|---|---|---|---|"]
    for key, v in report["raw_at_optimum"].items():
        s = report["raw_at_start"][key]
        lines.append(f"| {key} | {v['n']} | {s['mean']:+.2f} / {s['rms']:.2f} | {v['mean']:+.2f} / {v['rms']:.2f} | "
                     f"{v['tolerance']:.2f} | {'yes' if v['within'] else 'NO'} |")
    if report.get("diagnostic_rows_at_optimum"):
        lines += ["", "## Diagnostic rows (not fitted) at the optimum", "", "| data set \\| observable | n | mean | RMS |",
                  "|---|---|---|---|"]
        lines += [f"| {k} | {v['n']} | {v['mean']:+.2f} | {v['rms']:.2f} |" for k, v in report["diagnostic_rows_at_optimum"].items()]
    if "leave_one_laboratory_out" in report:
        lines += ["", "## Leave one laboratory out (noise re-estimated on the kept rows)", ""]
        for lab, res in report["leave_one_laboratory_out"].items():
            within = sum(v["within"] for v in res["held_out_raw"].values())
            lines.append(f"- **{lab}** ({res['rows_held']} rows held out): {within}/{len(res['held_out_raw'])} data set x "
                         f"observable within tolerance; constants beyond 2 sigma: {', '.join(res['beyond_2_sigma']) or 'none'}")
            lines += [f"  - {k}: mean {v['mean']:+.2f}, RMS {v['rms']:.2f} ({'within' if v['within'] else 'OUTSIDE'})"
                      for k, v in res["held_out_raw"].items()]
    if report.get("series_extra_noise"):
        lines += ["", "## Extra noise per data set x observable (within-series scatter)", ""]
        seen = {}
        for s, v in report["series_extra_noise"].items():
            seen.setdefault("|".join(s.split("|")[::2]), v)
        lines += [f"- `{k}`: {v:.3g}" for k, v in sorted(seen.items(), key=lambda kv: -kv[1])]
    lines += ["", fitting.diagnostics_markdown(report["diagnostics_raw"], "Residual diagnostics at the optimum (raw)"),
              fitting.diagnostics_markdown(report["diagnostics_whitened"], "Residual diagnostics at the optimum (whitened)")]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    sys.exit(main())
