#!/usr/bin/env python
"""
fit_deoxyosone_submodel.py -- sub-model 3's arms for methylglyoxal in pentose systems
(results/validation/kinetic_core_deoxyosone_submodel_prereg.md sec. 5, the revised design; 2026-10-30).

Engine state (named, frozen for every arm): the shipped B60 constants with the four pre-registered route steps on
(B58 TTCA equilibrium, B59b cysteine Amadori route, B59a Strecker sulfur, B61 norfuraneol hydroxide decay); sub-model
1's switches off (the winner is re-scored after). Every arm frees the methylglyoxal share of B61's decay
(NF_HYDROXIDE_DECAY.mgo_share, pinned by Mikami's MGO rows) and:

  A1  log10 multiplier on k_dpo_c2c3 (1-deoxypentosone -> MGO + glycolaldehyde)
  A2  A1 + log10 multiplier on k_dpo_ddp (the cysteine Strecker step: DPO's main exit, the parent deficit)
  A3  3-deoxypentosone -> MGO + glycolaldehyde (k at 100 C / pH 5, Ea; base-favoured)
  A4  pentose -> MGO + glycolaldehyde, amine-free (k at 100 C / pH 5, Ea; base-favoured)

Measurement model: the semi-quantitative 1-DX / 3-DX share one scale per figure (ledger id x species) with a prior of
0.3 dex (quinoxalines against a structurally close internal standard); each Zhou paper's MGO gets a free scale
(CF-060); Zhai's MGO is absolute; extra noise per data set x observable from within-series scatter. Rows the lane
cannot predict (glyoxal) and the declared hold-out (Zhai 2023 JAFC's 140 C pots) are left out. Leave-one-PAPER-out
(132 of 135 MGO rows are one laboratory's: leave-one-laboratory-out cannot be run, declared).

    python scripts/generators/fit_deoxyosone_submodel.py --arm A1|A2|A3|A4 [--budget N] [--rounds N] [--lopo]
Writes results/validation/kinetic_core_deoxyosone_submodel_fit_<arm>.{json,md}.
"""
from __future__ import annotations

import argparse
import dataclasses
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
from src.kinetic_core import engine, fit_rows, fitting, parallel, parameters_sulfur, sulfur  # noqa: E402

import build_submodel_probe as probe  # noqa: E402
from fit_ph_acid_submodel import pooled_noise, raw_table  # noqa: E402

TABLE = "deoxyosone_submodel"
ROUTE_SWITCHES = ("TTCA_FORMATION", "CYS_AMADORI_ROUTE", "STRECKER_SULFUR", "NF_HYDROXIDE_DECAY")
HOLDOUT_POTS_SUBSTRING = "140"   # Zhai 2023 JAFC's 140 C pots (prereg sec. 5)
HOLDOUT_LEDGER_PREFIX = "zhai2023jafc"
FIGURE_TAU = 0.3
SCALE_TAU = {**fitting.SCALE_TAU, "semi_figure": FIGURE_TAU, "zhou_free": float("inf")}


def coordinates(arm: str):
    keys = ["mgo_share_nf"]
    kinds = ["acid_yield"]           # a share in [0, 1]: the acid_yield kind's bounds and step
    lower, upper, start = [0.0], [1.0], [0.02]
    if arm in ("A1", "A2"):
        keys.append("offset_dpo_c2c3"); kinds.append("log10_offset"); lower.append(-1.0); upper.append(5.0); start.append(1.0)
    if arm == "A2":
        keys.append("offset_dpo_ddp"); kinds.append("log10_offset"); lower.append(-4.0); upper.append(1.0); start.append(-1.0)
    if arm == "A3":
        keys += ["log10_k_tdp_mgo", "ea_tdp_mgo"]; kinds += ["log10k", "Ea"]
        lower += [-8.0, 40.0]; upper += [0.0, 160.0]; start += [-4.0, 100.0]
    if arm == "A4":
        keys += ["log10_k_pent_mgo", "ea_pent_mgo"]; kinds += ["log10k", "Ea"]
        lower += [-9.0, 40.0]; upper += [-1.0, 160.0]; start += [-5.0, 100.0]
    return tuple(keys), tuple(kinds), np.array(lower, dtype=float), np.array(upper, dtype=float), np.array(start, dtype=float)


def trial(theta, keys) -> Dict[str, Any]:
    t = dict(zip(keys, map(float, theta)))
    base = engine.core_parameters(engine.SULFUR)
    k_ref = {}
    if "offset_dpo_c2c3" in t:
        k_ref["k_dpo_c2c3"] = base["k_dpo_c2c3"].k_ref * 10.0 ** t["offset_dpo_c2c3"]
    if "offset_dpo_ddp" in t:
        k_ref["k_dpo_ddp"] = base["k_dpo_ddp"].k_ref * 10.0 ** t["offset_dpo_ddp"]
    routes = {"k_tdp_mgo_100c_ph5": 10.0 ** t["log10_k_tdp_mgo"] if "log10_k_tdp_mgo" in t else 0.0,
              "ea_tdp_mgo_kj_mol": t.get("ea_tdp_mgo", 100.0),
              "k_pent_mgo_100c_ph5": 10.0 ** t["log10_k_pent_mgo"] if "log10_k_pent_mgo" in t else 0.0,
              "ea_pent_mgo_kj_mol": t.get("ea_pent_mgo", 100.0)}
    return {"k_ref": k_ref, "routes": routes, "mgo_share": t["mgo_share_nf"]}


def install(job: Dict[str, Any]) -> None:
    """The engine state and this trial's switches, set in whichever process runs the pot."""
    for name in ROUTE_SWITCHES:
        setattr(parameters_sulfur, name, dataclasses.replace(getattr(parameters_sulfur, name), enabled=True))
    parameters_sulfur.NF_HYDROXIDE_DECAY = dataclasses.replace(parameters_sulfur.NF_HYDROXIDE_DECAY,
                                                               mgo_share=float(job["mgo_share"]))
    parameters_sulfur.PENTOSE_MGO_ROUTES = parameters_sulfur.PentoseMgoRoutes(**job["routes"])


def run_task(job: Dict[str, Any]):
    install(job)
    try:
        return probe.run_pot_limited(job, POT_SECONDS)
    except Exception as exc:  # noqa: BLE001 - counted; its rows become failed residuals
        return {"failed": f"{type(exc).__name__}: {exc}"}


FAILURES: List[str] = []
#: wall-clock limit per pot integration (normal pots take 1-3 s; a stiff corner of parameter space can hang)
POT_SECONDS = 120


def make_raw_many(jobs_by_pot, keys):
    def raw_many(points, rows):
        pots = sorted({r.pot for r in rows})
        trials = [trial(x, keys) for x in points]
        results = parallel.map_ordered(run_task, [dict(jobs_by_pot[p], **tr) for tr in trials for p in pots])
        out = []
        for i, tr in enumerate(trials):
            res = dict(zip(pots, results[i * len(pots):(i + 1) * len(pots)]))
            r_i = []
            for r in rows:
                if "failed" in res[r.pot]:
                    FAILURES.append(f"{r.pot}: {res[r.pot]['failed']}")
                    r_i.append(float("nan"))
                    continue
                predicted = probe.predict(r, dict(jobs_by_pot[r.pot], **tr), res[r.pot])
                r_i.append(fit_rows.residual(r, predicted) if predicted > 0 else float("nan"))
            out.append(np.array(r_i, dtype=float))
        return out

    return raw_many


def method_of(row) -> str:
    if row.method in ("semi_quantitative", "relative"):
        return "semi_figure"
    if row.species == "MGO" and row.ledger_id.startswith("zhou"):
        return "zhou_free"
    return row.method


def paper_of(row) -> str:
    return row.ledger_id.split("_")[0]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--arm", choices=("A1", "A2", "A3", "A4"), required=True)
    parser.add_argument("--budget", type=int, default=240)
    parser.add_argument("--rounds", type=int, default=2)
    parser.add_argument("--lopo", action="store_true")
    parser.add_argument("--ph-coupling", choices=("fixed_point_nodes", "inline"), default="inline")
    args = parser.parse_args(argv)
    t0 = time.time()
    sulfur.PH_COUPLING = args.ph_coupling
    pots, all_rows = fit_rows.load(TABLE)
    jobs = {j["pot"]: j for j in probe.plan(pots, all_rows) if "skip" not in j}
    predictable = set(sulfur.SULFUR_STATE_KEYS)
    holdout = tuple(r for r in all_rows if r.ledger_id.startswith(HOLDOUT_LEDGER_PREFIX) and HOLDOUT_POTS_SUBSTRING in r.pot
                    and r.pot in jobs and r.species in predictable)
    data = tuple(r for r in all_rows if r.pot in jobs and r.species in predictable and r not in holdout)
    keys, kinds, lower, upper, start = coordinates(args.arm)
    raw_many = make_raw_many(jobs, keys)

    def raw(theta, rows):
        return raw_many([theta], rows)[0]

    model = fitting.MeasurementModel(
        series_of=lambda r: f"{r.ledger_id}|{r.species}" if method_of(r) in ("semi_figure", "zhou_free") else r.series,
        sigma_of=lambda r: r.sigma, method_of=method_of, scale_tau=SCALE_TAU)

    def make_problem(m, rows=data):
        return fitting.Problem(keys=keys, kinds=kinds, lower=lower, upper=upper, rows=rows,
                               residuals=lambda th, rs: m.whiten(raw(th, rs), rs),
                               evaluate_many=lambda pts, rs: [m.whiten(v, rs) for v in raw_many(pts, rs)],
                               group_of=paper_of, nuisance_dof=m.nuisance_dof)

    raw_start = raw(start, data)
    best, final_model = fitting.fit_with_series_noise(make_problem, model, raw, start, budget=args.budget,
                                                      rounds=args.rounds, estimator=pooled_noise)
    problem = make_problem(final_model)
    lap = fitting.laplace(problem, best.theta)
    raw_best = raw(best.theta, data)
    report: Dict[str, Any] = {
        "artifact": f"kinetic_core_deoxyosone_submodel_fit_{args.arm}",
        "provenance": provenance.provenance_block(
            f"kinetic_core_deoxyosone_submodel_fit_{args.arm}", generated_by="scripts/generators/fit_deoxyosone_submodel.py",
            inputs=[fit_rows.FIT_ROWS_DIR / f"{TABLE}.yml", ROOT / "src/kinetic_core/sulfur.py",
                    ROOT / "src/kinetic_core/parameters_sulfur.py", ROOT / "src/kinetic_core/fitting.py",
                    ROOT / "scripts/generators/build_submodel_probe.py"]),
        "arm": args.arm, "engine_state": "B60 + B58/B59b/B59a/B61 on; sub-model 1 switches off",
        "ph_coupling": args.ph_coupling, "rows": len(data), "pots": len({r.pot for r in data}),
        "start": dict(zip(keys, map(float, start))), "fit": best.as_dict(keys), "laplace": lap,
        "series_extra_noise": dict(final_model.extra_noise), "failed_pot_runs": len(FAILURES),
        "raw_at_start": raw_table(data, raw_start), "raw_at_optimum": raw_table(data, raw_best),
        "holdout_raw_at_optimum": raw_table(holdout, raw(best.theta, holdout)) if holdout else {},
        "diagnostics_whitened": fitting.diagnostics(problem, best.theta, {
            "paper": paper_of, "data_set": lambda r: r.ledger_id, "observable": lambda r: r.species,
            "temperature_C": lambda r: float(pots[r.pot]["t_c"]), "label_pH": lambda r: float(pots[r.pot]["label_ph"]),
            "time_min": lambda r: r.t_min}),
    }
    if args.lopo:
        out = {}
        sig = np.array([lap["sigma"][k] for k in keys])
        for paper in sorted({paper_of(r) for r in data}):
            kept = tuple(r for r in data if paper_of(r) != paper)
            held = tuple(r for r in data if paper_of(r) == paper)
            f_, _m = fitting.fit_with_series_noise(lambda mm: make_problem(mm, kept), model, raw, best.theta,
                                                   budget=args.budget // 2, rounds=1, estimator=pooled_noise)
            out[paper] = {"rows_held": len(held), "theta": dict(zip(keys, map(float, f_.theta))),
                          "shift_sigma": dict(zip(keys, map(float, (f_.theta - best.theta) / sig))),
                          "held_out_raw": raw_table(held, raw(f_.theta, held))}
        report["leave_one_paper_out"] = out
    report["wall_seconds"] = round(time.time() - t0, 1)
    path = data_paths.VALIDATION_DIR / f"kinetic_core_deoxyosone_submodel_fit_{args.arm}.json"
    path.write_text(provenance.dumps_compact(report, default=float) + "\n", encoding="utf-8")
    lines = [f"# Sub-model 3 fit, arm {args.arm}", "",
             f"*Generated by `scripts/generators/fit_deoxyosone_submodel.py --arm {args.arm}`; prereg sec. 5. "
             f"{report['rows']} rows on {report['pots']} pots; engine state: {report['engine_state']}. Raw residual = "
             "predicted minus measured (dex); semi-quantitative rows after their figure scale.*", "", "## Result", ""]
    lines += [f"- {line}" for line in fitting.summary_lines(problem, best, lap)] + [f"- failed pot runs: {len(FAILURES)}", ""]
    lines += ["## Raw residuals per data set x observable", "", "| data set \\| observable | n | start mean / RMS | optimum mean / RMS |",
              "|---|---|---|---|"]
    lines += [f"| {k} | {v['n']} | {report['raw_at_start'][k]['mean']:+.2f} / {report['raw_at_start'][k]['rms']:.2f} | "
              f"{v['mean']:+.2f} / {v['rms']:.2f} |" for k, v in report["raw_at_optimum"].items()]
    if report["holdout_raw_at_optimum"]:
        lines += ["", "## Hold-out (Zhai 2023 JAFC, 140 C), at the optimum", ""]
        lines += [f"- {k}: n {v['n']}, mean {v['mean']:+.2f}, RMS {v['rms']:.2f}" for k, v in report["holdout_raw_at_optimum"].items()]
    lines += ["", fitting.diagnostics_markdown(report["diagnostics_whitened"], "Residual diagnostics at the optimum (whitened)")]
    path.with_suffix(".md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"arm {args.arm}: cost {best.start_cost:.4g} -> {best.cost:.4g}; {len(FAILURES)} failed; {report['wall_seconds']} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
