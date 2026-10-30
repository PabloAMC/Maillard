#!/usr/bin/env python
"""
fit_joint_deoxyosone_thiol.py -- the joint sub-model 3 + 4 fit (results/validation/kinetic_core_thiol_submodel_prereg.md
sec. 4, and its design review; 2026-10-30).

Sub-model 3's parent arm (k_dpo_ddp slowed, k_dpo_c2c3 raised) closes the calibrated methylglyoxal but starves the
H2S that MFT (from norfuraneol) and FFT (from furfural) need: the cysteine Strecker step releases H2S (B59a). The
constants are fitted together on the deoxyosone and thiol tables, with the H2S loss free. Engine state as sub-model 3
(route switches B58 / B59b / B59a / B61 on, sub-model 1's off, inline coupling).

After its design review (thiol prereg sec. 5): coordinates are log10 multipliers on the shipped k_ref of k_dpo_ddp
(down to 10^-8), k_dpo_c2c3, k_dpo_nf, k_nf_mft, BOTH FFT branches together (k_fur_fft, k_fur_fft_hs), k_h2s_loss and
k_mft_decay (``--mft-decay-cap``: the CF-013-capped rival, at most 1/5 of the shipped loss), and the norfuraneol ->
MGO share. The B59a Strecker H2S share is held at a registered CR-11 alternative (``--strecker-share``).
Rows: both tables, deduplicated; OUT as declared hold-outs (FIT_HOLDOUT_DECLARATION): every Zhai 140 C pot, Yiltirak,
Hofmann T1 and T2 pH 3 / 7; SCORE-ONLY (no coordinate moves them): Kumazawa survival, Hofmann T8 and the T10
thiamine pot; dropped: glyoxal and the hexose + cysteine pots. Thiols semi-quantitative per figure with tau = inf
(CR-15), Zhou's thiols and MGO a per-paper scale (CF-014, CF-060), the quinoxalines tau 0.3 (CR-16). Extra noise
per data set x observable from within-series scatter, with a declared 0.10 dex floor where it cannot be estimated.
2-Acetylthiazole in Zhou 2023's ARP + Cys pot (11.70 ug/L) is scored as a falsifier, not fitted.

    python scripts/generators/fit_joint_deoxyosone_thiol.py --strecker-share 0.042|0.080 [--mft-decay-cap] [--lopo]
Writes results/validation/kinetic_core_joint_deoxyosone_thiol_fit_s<share>[_capped].{json,md}.
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
#: coordinate -> the parameter keys its multiplier moves
MULTIPLIED = {"k_dpo_ddp": ("k_dpo_ddp",), "k_dpo_c2c3": ("k_dpo_c2c3",), "k_dpo_nf": ("k_dpo_nf",),
              "k_nf_mft": ("k_nf_mft",), "fft_formation": ("k_fur_fft", "k_fur_fft_hs"), "k_h2s_loss": ("k_h2s_loss",),
              "k_mft_decay": ("k_mft_decay",)}
START = {"k_dpo_ddp": -1.0, "k_dpo_c2c3": 1.5, "k_dpo_nf": 0.0, "k_nf_mft": 0.0, "fft_formation": 0.0,
         "k_h2s_loss": 0.0, "k_mft_decay": -0.3}
LOWER = {"k_dpo_ddp": -8.0}
#: declared hold-outs (FIT_HOLDOUT_DECLARATION.md) and score-only rows (design review B3, S3)
HOLDOUT_LEDGER_PREFIXES = ("yiltirak2026", "hofmann1998_t1_", )
HOLDOUT_POTS = ("hofmann1998_cys_ribose_ph3", "hofmann1998_cys_ribose_ph7")
SCORE_ONLY_LEDGER_PREFIXES = ("kumazawa2003", "hofmann1998_t8_")
SCORE_ONLY_POTS = ("hofmann1998_thiamin_ph5",)
NOISE_FLOOR_DEX = 0.10
#: the falsifier (design review B4): Zhou 2023 Table 1, ARP + Cys 20 + 20 mM, 120 C, 60 min, pH 7, unbuffered
ACTZ_FALSIFIER = {"pot": "zhou2023_arp_cys_ph7_actz", "ledger_id": "zhou2023 T1 (B60 fit row zhou_pH7_ACTZ)",
                  "t_c": 120.0, "label_ph": 7.0, "charge": {"ARP": 20.0, "Cys": 20.0}, "dropped": [], "times": [0.0, 60.0],
                  "buffer": "none (water)", "notes": [], "measured_ug_l": 11.70}
HEXOSES = ("Glc", "Fru")
OUT = data_paths.VALIDATION_DIR / "kinetic_core_joint_deoxyosone_thiol_fit"
POT_SECONDS = 120
FAILURES: List[str] = []


def coordinates(mft_decay_cap: bool):
    keys = ["mgo_share_nf"] + [f"offset_{k}" for k in MULTIPLIED]
    kinds = ["acid_yield"] + ["log10_offset"] * len(MULTIPLIED)
    lower = [0.0] + [LOWER.get(k, -4.0) for k in MULTIPLIED]
    upper = [1.0] + [math.log10(0.2) if (k == "k_mft_decay" and mft_decay_cap) else 3.0 for k in MULTIPLIED]
    start = [0.0142] + [min(START[k], u) for k, u in zip(MULTIPLIED, upper[1:])]
    return tuple(keys), tuple(kinds), np.array(lower), np.array(upper), np.array(start, dtype=float)


def trial(theta, keys) -> Dict[str, Any]:
    t = dict(zip(keys, map(float, theta)))
    base = engine.core_parameters(engine.SULFUR)
    k_ref = {key: base[key].k_ref * 10.0 ** t[f"offset_{name}"] for name, keys_ in MULTIPLIED.items() for key in keys_}
    routes = {"k_tdp_mgo_100c_ph5": 0.0, "ea_tdp_mgo_kj_mol": 100.0, "k_pent_mgo_100c_ph5": 0.0, "ea_pent_mgo_kj_mol": 100.0}
    return {"k_ref": k_ref, "routes": routes, "mgo_share": t["mgo_share_nf"]}


STRECKER_SHARE = 0.042


def run_task(job):
    sm3.install(job)
    import dataclasses as _dc
    from src.kinetic_core import parameters_sulfur as _ps
    _ps.STRECKER_SULFUR = _dc.replace(_ps.STRECKER_SULFUR, h2s_share=float(job.get("strecker_share", STRECKER_SHARE)))
    try:
        return probe.run_pot_limited(job, POT_SECONDS)
    except Exception as exc:  # noqa: BLE001 - counted; its rows become failed residuals
        return {"failed": f"{type(exc).__name__}: {exc}"}


def key_of(row) -> str:
    return f"{row.table}:{row.pot}"


def method_of(row) -> str:
    """The quinoxalines (deoxyosone table) keep sub-model 3's per-figure tau 0.3 (CR-16); the thiol table's
    semi-quantitative rows are per-figure tau = inf (CR-15, CF-065); Zhou's thiols and MGO a per-paper scale."""
    if row.ledger_id.startswith("zhou") and row.species in ("MFT", "FFT", "MFTD", "FFTD", "FUR", "MGO"):
        return "zhou_free"
    if row.method in ("semi_quantitative", "relative"):
        return "semi_figure" if row.table == "deoxyosone_submodel" else "zhou_free"
    return row.method


def noise_with_floor(model, raw_residuals, rows, n_constants=0):
    """pooled_noise, and NOISE_FLOOR_DEX (declared; CF-016 / CF-055 between-run scatter) for every series whose
    noise cannot be estimated (design review S1: Garching's 10 % sd would otherwise dominate the cost)."""
    out = pooled_noise(model, raw_residuals, rows, n_constants)
    for row in rows:
        out.setdefault(model.series_of(row), NOISE_FLOOR_DEX)
    return out


def dedupe(rows):
    """Rows printed in both tables (design review B7) counted once: by (ledger id, species, time, value)."""
    seen, out = set(), []
    for r in rows:
        k = (r.ledger_id, r.species, float(r.t_min), round(float(r.value), 6))
        if k not in seen:
            seen.add(k)
            out.append(r)
    return tuple(out)


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


def actz_falsifier(theta, start, keys) -> Dict[str, float]:
    """2-Acetylthiazole in Zhou 2023's ARP + Cys pot at the optimum and at the start (design review B4)."""
    from src.kinetic_core.observations import amount_in

    out = {}
    for label, x in (("start", start), ("optimum", theta)):
        job = dict(ACTZ_FALSIFIER, **trial(x, keys), strecker_share=STRECKER_SHARE)
        res = run_task(job)
        out[f"{label}_ug_l"] = float("nan") if "failed" in res else amount_in("ug_L", "ACTZ", res["states"][-1]["ACTZ"])
    out["log10_ratio"] = math.log10(max(out["optimum_ug_l"], 1e-30) / ACTZ_FALSIFIER["measured_ug_l"])
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--budget", type=int, default=200)
    parser.add_argument("--rounds", type=int, default=2)
    parser.add_argument("--lopo", action="store_true")
    parser.add_argument("--strecker-share", type=float, choices=(0.042, 0.080), default=0.042)
    parser.add_argument("--mft-decay-cap", action="store_true")
    args = parser.parse_args(argv)
    global STRECKER_SHARE, OUT
    STRECKER_SHARE = float(args.strecker_share)
    OUT = data_paths.VALIDATION_DIR / (f"kinetic_core_joint_deoxyosone_thiol_fit_s{args.strecker_share:.3f}"
                                       + ("_capped" if args.mft_decay_cap else ""))
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
    rows_all = list(dedupe(r for r in rows_all if key_of(r) in jobs))

    def is_holdout(r):
        return ((r.ledger_id.startswith("zhai2023") and float(pots_all[key_of(r)]["t_c"]) >= 140.0)
                or r.ledger_id.startswith(HOLDOUT_LEDGER_PREFIXES) or r.pot in HOLDOUT_POTS)

    def is_score_only(r):
        return r.ledger_id.startswith(SCORE_ONLY_LEDGER_PREFIXES) or r.pot in SCORE_ONLY_POTS

    holdout = tuple(r for r in rows_all if is_holdout(r))
    score_only = tuple(r for r in rows_all if not is_holdout(r) and is_score_only(r))
    data = tuple(r for r in rows_all if not is_holdout(r) and not is_score_only(r))
    keys, kinds, lower, upper, start = coordinates(args.mft_decay_cap)
    raw_many = make_raw_many(jobs, keys)

    def raw(theta, rows):
        return raw_many([theta], rows)[0]

    model = fitting.MeasurementModel(
        series_of=lambda r: f"{r.ledger_id}|{r.species}" if method_of(r) in ("semi_figure", "zhou_free") else r.series,
        sigma_of=lambda r: r.sigma, method_of=method_of, scale_tau=sm3.SCALE_TAU)

    def make_problem(m, rows=data):
        return fitting.Problem(keys=keys, kinds=kinds, lower=lower, upper=upper, rows=rows,
                               residuals=lambda th, rs: m.whiten(raw(th, rs), rs),
                               evaluate_many=lambda pts, rs: [m.whiten(v, rs) for v in raw_many(pts, rs)],
                               group_of=sm3.paper_of, nuisance_dof=m.nuisance_dof)

    raw_start = raw(start, data)
    best, final_model = fitting.fit_with_series_noise(make_problem, model, raw, start, budget=args.budget,
                                                      rounds=args.rounds, estimator=noise_with_floor)
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
        "strecker_share": STRECKER_SHARE, "mft_decay_cap": bool(args.mft_decay_cap),
        "holdout_raw_at_optimum": raw_table(holdout, raw(best.theta, holdout)) if holdout else {},
        "holdout_raw_at_start": raw_table(holdout, raw(start, holdout)) if holdout else {},
        "score_only_raw_at_optimum": raw_table(score_only, raw(best.theta, score_only)) if score_only else {},
        "actz_falsifier": actz_falsifier(best.theta, start, keys),
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
                                                   budget=args.budget // 2, rounds=1, estimator=noise_with_floor)
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
    fa = report["actz_falsifier"]
    lines += ["", "## Falsifier: 2-acetylthiazole, Zhou 2023 ARP + Cys pH 7 (measured 11.70 ug/L; not fitted)", "",
              f"- predicted at the start {fa['start_ug_l']:.4g} ug/L, at the optimum {fa['optimum_ug_l']:.4g} ug/L "
              f"(log10 ratio at the optimum {fa['log10_ratio']:+.2f})"]
    if report["holdout_raw_at_optimum"]:
        lines += ["", "## Declared hold-outs (every Zhai 140 C pot, Yiltirak, Hofmann T1 / T2 pH 3 / 7), scored once", "",
                  "| data set \\| observable | n | start mean / RMS | optimum mean / RMS |", "|---|---|---|---|"]
        lines += [f"| {k} | {v['n']} | {report['holdout_raw_at_start'][k]['mean']:+.2f} / {report['holdout_raw_at_start'][k]['rms']:.2f} | "
                  f"{v['mean']:+.2f} / {v['rms']:.2f} |" for k, v in report["holdout_raw_at_optimum"].items()]
    if report["score_only_raw_at_optimum"]:
        lines += ["", "## Score-only rows at the optimum", ""]
        lines += [f"- {k}: n {v['n']}, mean {v['mean']:+.2f}, RMS {v['rms']:.2f}" for k, v in report["score_only_raw_at_optimum"].items()]
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
