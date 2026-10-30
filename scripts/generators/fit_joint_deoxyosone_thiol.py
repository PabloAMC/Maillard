#!/usr/bin/env python
"""
fit_joint_deoxyosone_thiol.py -- the joint sub-model 3 + 4 fit (results/validation/kinetic_core_thiol_submodel_prereg.md
sec. 4, and its design review; 2026-10-30).

Sub-model 3's parent arm (k_dpo_ddp slowed, k_dpo_c2c3 raised) closes the calibrated methylglyoxal but starves the
H2S that MFT (from norfuraneol) and FFT (from furfural) need: the cysteine Strecker step releases H2S (B59a). The
constants are fitted together on the deoxyosone and thiol tables, with the H2S loss free. Engine state as sub-model 3
(route switches B58 / B59b / B59a / B61 on, sub-model 1's off, inline coupling).

Coordinates: log10 multipliers on the shipped k_ref of k_dpo_ddp, k_dpo_c2c3, k_dpo_nf, k_nf_mft, k_fur_fft and
k_h2s_loss, and the norfuraneol -> MGO share. Rows: both tables less glyoxal, the hexose + cysteine pots (MFT is
structurally zero until the hexose -> C5 entry exists) and the hold-out (Zhai 2023 JAFC's 140 C pots, both tables).

    python scripts/generators/fit_joint_deoxyosone_thiol.py [--budget N] [--rounds N] [--lopo]
Writes results/validation/kinetic_core_joint_deoxyosone_thiol_fit.{json,md}.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for entry in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

from src import data_paths, provenance  # noqa: E402
from src.kinetic_core import engine, fit_rows, fitting, parallel, sulfur  # noqa: E402

import build_submodel_probe as probe  # noqa: E402
import fit_deoxyosone_submodel as sm3  # noqa: E402
from fit_ph_acid_submodel import pooled_noise, raw_table  # noqa: E402

TABLES = ("deoxyosone_submodel", "thiol_submodel")
MULTIPLIED = ("k_dpo_ddp", "k_dpo_c2c3", "k_dpo_nf", "k_nf_mft", "k_fur_fft", "k_h2s_loss")
START = {"k_dpo_ddp": -1.0, "k_dpo_c2c3": 1.5, "k_dpo_nf": 0.0, "k_nf_mft": 0.0, "k_fur_fft": 0.0, "k_h2s_loss": 0.0}
HEXOSES = ("Glc", "Fru")
OUT = data_paths.VALIDATION_DIR / "kinetic_core_joint_deoxyosone_thiol_fit"
POT_SECONDS = 120
FAILURES: List[str] = []


def coordinates():
    keys = ["mgo_share_nf"] + [f"offset_{k}" for k in MULTIPLIED]
    kinds = ["acid_yield"] + ["log10_offset"] * len(MULTIPLIED)
    lower = [0.0] + [-4.0] * len(MULTIPLIED)
    upper = [1.0] + [3.0] * len(MULTIPLIED)
    start = [0.0142] + [START[k] for k in MULTIPLIED]
    return tuple(keys), tuple(kinds), np.array(lower), np.array(upper), np.array(start, dtype=float)


def trial(theta, keys) -> Dict[str, Any]:
    t = dict(zip(keys, map(float, theta)))
    base = engine.core_parameters(engine.SULFUR)
    k_ref = {k: base[k].k_ref * 10.0 ** t[f"offset_{k}"] for k in MULTIPLIED}
    routes = {"k_tdp_mgo_100c_ph5": 0.0, "ea_tdp_mgo_kj_mol": 100.0, "k_pent_mgo_100c_ph5": 0.0, "ea_pent_mgo_kj_mol": 100.0}
    return {"k_ref": k_ref, "routes": routes, "mgo_share": t["mgo_share_nf"]}


def run_task(job):
    sm3.install(job)
    try:
        return probe.run_pot_limited(job, POT_SECONDS)
    except Exception as exc:  # noqa: BLE001 - counted; its rows become failed residuals
        return {"failed": f"{type(exc).__name__}: {exc}"}


def key_of(row) -> str:
    return f"{row.table}:{row.pot}"


def make_raw_many(jobs, keys):
    def raw_many(points, rows):
        pots = sorted({key_of(r) for r in rows})
        trials = [trial(x, keys) for x in points]
        results = parallel.map_ordered(run_task, [dict(jobs[p], **tr) for tr in trials for p in pots])
        out = []
        for i, tr in enumerate(trials):
            res = dict(zip(pots, results[i * len(pots):(i + 1) * len(pots)]))
            r_i = []
            for r in rows:
                got = res[key_of(r)]
                if "failed" in got:
                    FAILURES.append(f"{key_of(r)}: {got['failed']}")
                    r_i.append(float("nan"))
                    continue
                predicted = probe.predict(r, dict(jobs[key_of(r)], **tr), got)
                r_i.append(fit_rows.residual(r, predicted) if predicted > 0 else float("nan"))
            out.append(np.array(r_i, dtype=float))
        return out

    return raw_many


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--budget", type=int, default=200)
    parser.add_argument("--rounds", type=int, default=2)
    parser.add_argument("--lopo", action="store_true")
    args = parser.parse_args(argv)
    t0 = time.time()
    sulfur.PH_COUPLING = "inline"
    predictable = set(sulfur.SULFUR_STATE_KEYS)
    jobs, rows_all, pots_all = {}, [], {}
    for table in TABLES:
        pots, rows = fit_rows.load(table)
        for j in probe.plan(pots, rows):
            if "skip" in j or any(pots[j["pot"]]["initial_mmol_L"].get(h) for h in HEXOSES):
                continue
            jobs[f"{table}:{j['pot']}"] = j
            pots_all[f"{table}:{j['pot']}"] = pots[j["pot"]]
        rows_all += [r for r in rows if r.species in predictable]
    holdout = tuple(r for r in rows_all if key_of(r) in jobs and r.ledger_id.startswith(sm3.HOLDOUT_LEDGER_PREFIX)
                    and sm3.HOLDOUT_POTS_SUBSTRING in r.pot)
    data = tuple(r for r in rows_all if key_of(r) in jobs and r not in holdout)
    keys, kinds, lower, upper, start = coordinates()
    raw_many = make_raw_many(jobs, keys)

    def raw(theta, rows):
        return raw_many([theta], rows)[0]

    model = fitting.MeasurementModel(
        series_of=lambda r: f"{r.ledger_id}|{r.species}" if sm3.method_of(r) in ("semi_figure", "zhou_free") else r.series,
        sigma_of=lambda r: r.sigma, method_of=sm3.method_of, scale_tau=sm3.SCALE_TAU)

    def make_problem(m, rows=data):
        return fitting.Problem(keys=keys, kinds=kinds, lower=lower, upper=upper, rows=rows,
                               residuals=lambda th, rs: m.whiten(raw(th, rs), rs),
                               evaluate_many=lambda pts, rs: [m.whiten(v, rs) for v in raw_many(pts, rs)],
                               group_of=sm3.paper_of, nuisance_dof=m.nuisance_dof)

    raw_start = raw(start, data)
    best, final_model = fitting.fit_with_series_noise(make_problem, model, raw, start, budget=args.budget,
                                                      rounds=args.rounds, estimator=pooled_noise)
    problem = make_problem(final_model)
    lap = fitting.laplace(problem, best.theta)
    raw_best = raw(best.theta, data)
    report: Dict[str, Any] = {
        "artifact": OUT.name,
        "provenance": provenance.provenance_block(
            OUT.name, generated_by="scripts/generators/fit_joint_deoxyosone_thiol.py",
            inputs=[fit_rows.FIT_ROWS_DIR / f"{t}.yml" for t in TABLES]
                   + [ROOT / "src/kinetic_core/sulfur.py", ROOT / "src/kinetic_core/parameters_sulfur.py",
                      ROOT / "src/kinetic_core/fitting.py", ROOT / "scripts/generators/build_submodel_probe.py"]),
        "engine_state": "B60 + B58/B59b/B59a/B61 on; sub-model 1 switches off; inline coupling",
        "rows": len(data), "pots": len({key_of(r) for r in data}), "start": dict(zip(keys, map(float, start))),
        "fit": best.as_dict(keys), "laplace": lap, "series_extra_noise": dict(final_model.extra_noise),
        "failed_pot_runs": len(FAILURES), "failures_sample": FAILURES[:10],
        "raw_at_start": raw_table(data, raw_start), "raw_at_optimum": raw_table(data, raw_best),
        "holdout_raw_at_optimum": raw_table(holdout, raw(best.theta, holdout)) if holdout else {},
        "diagnostics_whitened": fitting.diagnostics(problem, best.theta, {
            "paper": sm3.paper_of, "data_set": lambda r: r.ledger_id, "observable": lambda r: r.species,
            "temperature_C": lambda r: float(pots_all[key_of(r)]["t_c"]),
            "label_pH": lambda r: float(pots_all[key_of(r)]["label_ph"]), "time_min": lambda r: r.t_min}),
    }
    sulfur.PH_COUPLING = "fixed_point_nodes"
    report["shipped_coupling_max_abs_difference"] = float(np.nanmax(np.abs(raw(best.theta, data) - raw_best)))
    sulfur.PH_COUPLING = "inline"
    if args.lopo:
        out, sig = {}, np.array([lap["sigma"][k] for k in keys])
        for paper in sorted({sm3.paper_of(r) for r in data}):
            kept = tuple(r for r in data if sm3.paper_of(r) != paper)
            held = tuple(r for r in data if sm3.paper_of(r) == paper)
            f_, _m = fitting.fit_with_series_noise(lambda mm: make_problem(mm, kept), model, raw, best.theta,
                                                   budget=args.budget // 2, rounds=1, estimator=pooled_noise)
            out[paper] = {"rows_held": len(held), "theta": dict(zip(keys, map(float, f_.theta))),
                          "shift_sigma": dict(zip(keys, map(float, (f_.theta - best.theta) / sig))),
                          "held_out_raw": raw_table(held, raw(f_.theta, held))}
        report["leave_one_paper_out"] = out
    report["wall_seconds"] = round(time.time() - t0, 1)
    OUT.with_suffix(".json").write_text(json.dumps(report, indent=1, default=float) + "\n", encoding="utf-8")
    lines = ["# Joint sub-model 3 + 4 fit (deoxyosones + thiols)", "",
             "*Generated by `scripts/generators/fit_joint_deoxyosone_thiol.py`; thiol prereg sec. 4. "
             f"{report['rows']} rows on {report['pots']} pots; engine state: {report['engine_state']}. Raw residual = "
             "predicted minus measured (dex); semi-quantitative rows after their figure scale.*", "", "## Result", ""]
    lines += [f"- {line}" for line in fitting.summary_lines(problem, best, lap)]
    lines += [f"- failed pot runs: {len(FAILURES)}",
              f"- largest raw-residual difference, inline vs the shipped fixed point, at the optimum: "
              f"{report['shipped_coupling_max_abs_difference']:.2g}", "",
              "## Raw residuals per data set x observable", "",
              "| data set \\| observable | n | start mean / RMS | optimum mean / RMS |", "|---|---|---|---|"]
    lines += [f"| {k} | {v['n']} | {report['raw_at_start'][k]['mean']:+.2f} / {report['raw_at_start'][k]['rms']:.2f} | "
              f"{v['mean']:+.2f} / {v['rms']:.2f} |" for k, v in report["raw_at_optimum"].items()]
    if report["holdout_raw_at_optimum"]:
        lines += ["", "## Hold-out (Zhai 2023 JAFC, 140 C), at the optimum", ""]
        lines += [f"- {k}: n {v['n']}, mean {v['mean']:+.2f}, RMS {v['rms']:.2f}" for k, v in report["holdout_raw_at_optimum"].items()]
    if "leave_one_paper_out" in report:
        lines += ["", "## Leave one paper out", ""]
        for paper, res in report["leave_one_paper_out"].items():
            lines.append(f"- **{paper}** ({res['rows_held']} rows): " + "; ".join(
                f"{k} {v['mean']:+.2f} / {v['rms']:.2f}" for k, v in res["held_out_raw"].items()))
    lines += ["", fitting.diagnostics_markdown(report["diagnostics_whitened"], "Residual diagnostics at the optimum (whitened)")]
    OUT.with_suffix(".md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"joint: cost {best.start_cost:.4g} -> {best.cost:.4g}; {len(FAILURES)} failed; {report['wall_seconds']} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
