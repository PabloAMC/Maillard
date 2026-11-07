#!/usr/bin/env python
"""
m6_reduced_power_check.py -- case M6 (pyrazines) secs. 4 - 5: the falsification arithmetic on spent / seen data and the
reduced-model power check of the designed decider under the ADOPTED sec. 7 item 5 / M1.1 judge
(src/kinetic_core/predictive_densities.py: UnitNoiseModel, build_null_pool, judge_from_draws / judge_pair).

The mini-model is src/m6_reduced_power.py; the case is docs/reference/cases/M6_pyrazines.md. Writes
results/validation/kinetic_core_m6_reduced_power.{json,md}. SYNTHETIC data only for job (c); job (a) is arithmetic on
the printed spent / seen values the case quotes. No engine run, no engine fit, no decisive, reserved or TEST value.

    python scripts/generators/m6_reduced_power_check.py --checks      # build checks only (fast)
    python scripts/generators/m6_reduced_power_check.py                # jobs (a) and (c); the record

Rules enforced here: a memory guard (>= 35 % free before the heavy part; --no-memory-guard for tests); threads pinned
to 1 by the caller (OMP_NUM_THREADS=1).
"""
from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
import time
import zlib
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import m6_reduced_power as M  # noqa: E402
from src import provenance  # noqa: E402
from src.kinetic_core import predictive_densities as pdn  # noqa: E402

STEM = "kinetic_core_m6_reduced_power"
CASE = "docs/reference/cases/M6_pyrazines.md (branch m6-case-ready)"
MEMORY_MIN_FREE_PCT = 35
PREC = 100.0      # the declared precision requirement: ESS >= max(100, PREC) (the shared floor)
RETRY_MULT = 3    # the re-estimation ensembles are this many times the scoring ensembles

#: (pair, designs): each design a tuple of units; "full" first, then the leave-one-arm-out and two-arm variants
PART_C = (("C1", "C0"), ("C1", "C3"), ("C0", "C3"))
PART_S = (("S1", "S0"), ("S1", "S5"), ("S5", "S0"))


def memory_free_pct() -> Optional[int]:
    try:
        out = subprocess.run(["memory_pressure"], capture_output=True, text=True, timeout=60).stdout
    except Exception:
        return None
    m = re.search(r"System-wide memory free percentage:\s*(\d+)%", out)
    return int(m.group(1)) if m else None


def memory_guard(enabled: bool) -> Dict[str, Any]:
    free = memory_free_pct()
    if enabled and free is not None and free < MEMORY_MIN_FREE_PCT:
        raise SystemExit(f"memory guard: {free} % free < {MEMORY_MIN_FREE_PCT} %; refusing to start")
    return {"free_pct_at_start": free, "guard": enabled}


#: the designs per part (name -> units, the n_rep values run); the first is the primary design
DESIGNS_C = {
    "M": (("mgo_step",), (2, 3)),
    "D": (("dilution_step",), (2,)),
    "M+D": (("mgo_step", "dilution_step"), (2,)),
    "M+D+K": (("mgo_step", "dilution_step", "ak_trap"), (2,)),
}
DESIGNS_S = {
    "plain+phosphate": (("budget_plain", "budget_phosphate"), (2, 3)),
    "plain only": (("budget_plain",), (2,)),
    "phosphate only": (("budget_phosphate",), (2,)),
}


def _ens(d: pdn.RivalDraws) -> Dict[str, Dict[str, pdn.RivalDraws]]:
    return {pdn.BASE_FORM: {pdn.BASE_FORM: d}}


def run_pair(a: str, b: str, units: Sequence[str], spec: M.NoiseSpec, n_sets: int, n_draws: int, n_pool: int,
             ppc_reps: int, seed: int, extra_factors: Sequence[float]) -> Dict[str, Any]:
    """Power and false credit of the designed units for the pair (a, b), each rival as truth; then the misspecified
    false credit (the generator adds ``extra_factor`` x the per-point SD beyond the scoring noise)."""
    rng = np.random.default_rng(seed)
    noise = M.noise_models(units, spec)
    design = pdn.MeasurementDesign(noise)
    ea, eb = (M.rival_draws(r, units, n_draws, rng) for r in (a, b))
    ta, tb = (M.rival_draws(r, units, n_draws, rng, defensive=False) for r in (a, b))
    sa, sb = (M.rival_draws(r, units, n_draws, rng, defensive=False) for r in (a, b))
    t0 = time.time()
    # the re-estimation set (M1.1, PR #105): larger, independent ensembles, used for every unresolved replicate and
    # for an observed set whose densities the first ensembles do not resolve
    ra, rb = (M.rival_draws(r, units, RETRY_MULT * n_draws, rng) for r in (a, b))
    pa = pdn.build_null_pool("a", {pdn.BASE_FORM: ta}, design, _ens(ea), _ens(eb), n_pool, rng, precision_ess=PREC,
                             retry_draws=(_ens(ra), _ens(rb)))
    pb = pdn.build_null_pool("b", {pdn.BASE_FORM: tb}, design, _ens(ea), _ens(eb), n_pool, rng, precision_ess=PREC,
                             retry_draws=(_ens(ra), _ens(rb)))
    pool_s = time.time() - t0
    retried = [0]
    # The M1.1 gate is one number per pair: a set is judged once with a supplied, passing placeholder, and each gate
    # scenario then turns every credit into "undecided" when its false credit fails (exactly what judge_pair does).
    passing = (0, max(n_sets, pdn.power_check_min_n()))

    def tally(truth_draws: pdn.RivalDraws, extra: Optional[Mapping[str, float]]) -> Dict[str, int]:
        out: Dict[str, int] = {}
        for _ in range(n_sets):
            j = int(rng.integers(truth_draws.n))
            y = {}
            for u in units:
                f = truth_draws.pred(u)[j][None, :]
                yu = noise[u].sample(f, rng)[0]
                if extra:
                    yu = yu + rng.standard_normal(yu.size) * extra[u]
                y[u] = yu
            v = pdn.judge_from_draws(y, design, _ens(ea), _ens(eb), pa, pb, rng, ppc_reps=ppc_reps,
                                     precision_ess=PREC, misspecified_false_credit=passing)
            if v.outcome == "not_computable" and any("unresolved density" in r for r in v.reasons):
                retried[0] += 1
                v = pdn.judge_from_draws(y, design, _ens(ra), _ens(rb), pa, pb, rng, ppc_reps=ppc_reps,
                                         precision_ess=PREC, misspecified_false_credit=passing)
            key = f"credited_{v.credited}" if v.credited else v.outcome
            out[key] = out.get(key, 0) + 1
        return out

    res = {"a_true": tally(sa, None), "b_true": tally(sb, None)}
    k_spec = res["a_true"].get("credited_b", 0) + res["b_true"].get("credited_a", 0)
    gates = {"0.0": pdn.misspecified_false_credit_gate(k_spec, 2 * n_sets)}
    mis = {}
    for fac in extra_factors:
        if fac <= 0:
            continue
        extra = M.misspecified_extra_sd(units, spec, fac)
        t = {"a_true": tally(sa, extra), "b_true": tally(sb, extra)}
        k = t["a_true"].get("credited_b", 0) + t["b_true"].get("credited_a", 0)
        mis[str(fac)] = {"outcomes": t, "k": k, "n": 2 * n_sets, "clopper_pearson95": M.clopper_pearson(k, 2 * n_sets)}
        gates[str(fac)] = pdn.misspecified_false_credit_gate(k, 2 * n_sets)

    def rate(t: Mapping[str, int], key: str) -> Dict[str, Any]:
        k = int(t.get(key, 0))
        return {"k": k, "n": n_sets, "rate": k / n_sets, "clopper_pearson95": M.clopper_pearson(k, n_sets),
                "wilson95": pdn.wilson_interval(k, n_sets)}
    power_under_gate = {f: {"gate_passes": g["passes"],
                            "power_a": (res["a_true"].get("credited_a", 0) / n_sets) if g["passes"] else 0.0,
                            "power_b": (res["b_true"].get("credited_b", 0) / n_sets) if g["passes"] else 0.0}
                        for f, g in gates.items()}
    return {"pair": [a, b], "units": list(units), "spec": spec.__dict__, "n_sets": n_sets, "n_draws": n_draws,
            "n_pool": n_pool, "ppc_reps": ppc_reps, "pool_seconds": pool_s,
            "pool_unresolved": {"a": float(np.mean(~pa.resolved)), "b": float(np.mean(~pb.resolved))},
            "pool_reestimation": {"a": pa.extra, "b": pb.extra}, "observed_sets_reestimated": retried[0],
            "power_a": rate(res["a_true"], "credited_a"), "power_b": rate(res["b_true"], "credited_b"),
            "false_credit_b_when_a": rate(res["a_true"], "credited_b"),
            "false_credit_a_when_b": rate(res["b_true"], "credited_a"),
            "outcomes": res, "gates": gates, "misspecified": mis, "power_under_gate": power_under_gate,
            "note": "power_a / power_b and false credit are BEFORE the gate (placeholder passing); power_under_gate "
                    "applies each scenario's gate: '0.0' = the correctly specified generator (a run-level "
                    "reproducibility term cancels in every within-run unit), '0.5' / '1.0' = an unmodelled per-point "
                    "term of that multiple of the declared per-point SD"}


def job_c(args) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    seed = args.seed
    plan = []
    if "C" in args.parts:
        plan += [(pair, DESIGNS_C) for pair in PART_C if not args.pairs or f"{pair[0]}_vs_{pair[1]}" in args.pairs]
    if "S" in args.parts:
        plan += [(pair, DESIGNS_S) for pair in PART_S if not args.pairs or f"{pair[0]}_vs_{pair[1]}" in args.pairs]
    for (a, b), designs in plan:
        for dname, (units, n_reps) in designs.items():
            for n_rep in n_reps:
                key = f"{a}_vs_{b}|{dname}|n{n_rep}"
                seed = args.seed + zlib.crc32(key.encode()) % 100000       # per key: splits reproduce
                if args.only and args.only not in key:
                    continue
                spec = M.NoiseSpec(cv_dex=args.cv_dex, n_rep=n_rep)
                t0 = time.time()
                out[key] = run_pair(a, b, units, spec, args.n_sets, args.n_draws, args.n_pool, args.ppc_reps,
                                    seed, args.extra_factors)
                out[key]["seconds"] = time.time() - t0
                r = out[key]
                print(f"{key}: P(a) {r['power_a']['rate']:.2f} P(b) {r['power_b']['rate']:.2f} "
                      f"FC {r['false_credit_b_when_a']['rate']:.2f}/{r['false_credit_a_when_b']['rate']:.2f} "
                      f"gates {[(f, g['k'], g['passes']) for f, g in r['gates'].items()]} unres {r['pool_unresolved']} "
                      f"({r['seconds']:.0f} s)", flush=True)
                if args.partial:
                    Path(args.partial).write_text(json.dumps(out, default=_jsonable))
    return out


def checks() -> Dict[str, Any]:
    rng = np.random.default_rng(0)
    c = {}
    for r in ("C0", "C1", "C3"):
        p = M.predict_c(M.draw(r, 64, rng)[0], r)
        c[r] = {u: [float(np.nanmin(v)), float(np.nanmax(v))] for u, v in p.items()}
        c[r]["finite"] = all(bool(np.all(np.isfinite(v))) for v in p.values())
    c["C0_mgo_step_is_log3"] = abs(c["C0"]["mgo_step"][0] - math.log10(3)) < 0.06
    c["C3_dilution_is_log_1_256"] = abs(c["C3"]["dilution_step"][0] - math.log10(1 / 256)) < 0.12
    c["pass"] = all(c[r]["finite"] for r in ("C0", "C1", "C3")) and c["C0_mgo_step_is_log3"] and \
        c["C3_dilution_is_log_1_256"]
    return c


def summarise(jc: Mapping[str, Any]) -> List[Dict[str, Any]]:
    rows = []
    for k, r in jc.items():
        pair, design, n = k.split("|")
        nc = sum(v for t in r["outcomes"].values() for kk, v in t.items() if kk == "not_computable") / (2 * r["n_sets"])
        rows.append({"pair": pair, "design": design, "n_rep": int(n[1:]), "units": r["units"],
                     "P_a": r["power_a"]["rate"], "P_a_cp95": r["power_a"]["clopper_pearson95"],
                     "P_b": r["power_b"]["rate"], "P_b_cp95": r["power_b"]["clopper_pearson95"],
                     "false_credit": max(r["false_credit_b_when_a"]["rate"], r["false_credit_a_when_b"]["rate"]),
                     "false_credit_cp95_upper": max(r["false_credit_b_when_a"]["clopper_pearson95"][1],
                                                    r["false_credit_a_when_b"]["clopper_pearson95"][1]),
                     "gates": {f: {"k": g["k"], "n": g["n"], "wilson95": g["wilson95"], "passes": g["passes"]}
                               for f, g in r["gates"].items()},
                     "P_separates_correct_spec": min(r["power_under_gate"]["0.0"]["power_a"],
                                                     r["power_under_gate"]["0.0"]["power_b"]),
                     "not_computable": nc, "pool_unresolved": r["pool_unresolved"]})
    return rows


def render(o: Mapping[str, Any]) -> str:
    a = o["job_a"]
    L = [f"# Case M6: falsification arithmetic and reduced-model power check ({o['provenance']['generated_on']})", "",
         f"Case: `{CASE}`. Model: `src/m6_reduced_power.py`. SYNTHETIC data for job (c); job (a) is arithmetic on the "
         "spent / seen printed values. No engine run or fit; no decisive, reserved or TEST value read. Judge: the "
         "adopted sec. 7 item 5 / M1.1 (`src/kinetic_core/predictive_densities.py`).", "",
         "## Job (a): falsifications on data on disk (falsify-only; spent data never credit)", "",
         "**FX-C1, Zhou 2024's six zero-order R2 vs a free-aminoketone pool with Elfsson's class and no first-order "
         "exit.** p = max over noise SD and form of P(R2 >= printed):", "",
         "| unit | tau (h) | R2 noiseless | p max | p Holm | k2 needed for a 10-min lag / Elfsson's fastest corner (dex) |",
         "|---|---:|---:|---:|---:|---:|"]
    for u, v in a["onset"]["no_exit"].items():
        h = a["onset"]["no_exit_holm"][u]
        L.append(f"| {u} | {v['tau_h']:.1f} | {v['r2_noiseless']:.3f} | {v['p_max']:.4f} | {h['p_holm']:.4f} | "
                 f"{v['gap_dex']:.2f} |")
    imp = a["onset"].get("implied_at_ks_min") or {}
    L += ["", f"Falsified: **{a['onset']['no_exit_falsified']}**. Smallest first-order AK exit that survives (Holm, "
          f"every unit): **ks >= {a['onset']['ks_min_per_min']:.3g} /min** at 20 mM MGO. At that edge (100 C, pH 8, "
          f"Elfsson's fastest corner): AK {imp.get('AK_steady_uM', float('nan')):.0f} uM, Strecker flux "
          f"{imp.get('strecker_flux_uM_min', float('nan')):.2f} uM/min, kS {imp.get('kS_M_min', float('nan')):.2e} "
          f"M^-1 min^-1 = B18 + {imp.get('kS_over_B18_dex', float('nan')):.2f} dex = Hofmann-implied "
          f"{imp.get('kS_over_hofmann_implied_dex', float('nan')):+.2f} dex; DMP per Strecker event "
          f"{imp.get('yield_DMP_per_strecker_event', float('nan')):.4f}; exit product "
          f"{imp.get('pyrrole_uM_at_120min', float('nan')):.0f} uM by 120 min.", ""]
    s = a["strecker_only_exit"]
    L += ["**FX-S1, the incumbent's 'Strecker is the only amine exit'.**", "",
          f"- Hofmann 2000 Table 3, within one pot: (PA + PAA) / Phe consumed {s['hofmann']['rows'][0]['products_over_consumed']:.3f} "
          f"at 10 min ({-s['hofmann']['deficit_10min_dex']:.2f} dex) to {s['hofmann']['rows'][-1]['products_over_consumed']:.3f} "
          f"at 500 min; one-sided p under omega on each value {s['hofmann']['p_one_sided_at_omega']:.3f}; break-even SD "
          f"{s['hofmann']['break_even_sd_dex']:.2f} dex. Flag: {s['hofmann']['integrity_flag']}.",
          f"- Piloty 1979 Ala loss vs Zhou's DMP under yield 1: gap {s['piloty_vs_zhou']['rows']['pH6.97']['gap_dex']:.2f} dex "
          f"(pH 6.97; an upper-bound MGO exposure); p {s['piloty_vs_zhou']['p_one_sided_no_buffer']:.3f} (no buffer "
          f"term) / {s['piloty_vs_zhou']['p_one_sided_buffer_max']:.3f} (phosphate at M5's kappa edge), cross-laboratory "
          f"SD {s['piloty_vs_zhou']['sd_cross_lab_tau_abs']:.2f} dex (tau_ABS and omega, the largest listed); break-even "
          f"between-laboratory SD {s['piloty_vs_zhou']['break_even_between_lab_sd_no_buffer']:.2f} / "
          f"{s['piloty_vs_zhou']['break_even_between_lab_sd_buffer_max']:.2f} dex.",
          f"- Holm: {json.dumps({k: round(v['p_holm'], 3) for k, v in s['holm'].items()})}; falsified: **{s['falsified']}**.",
          "", f"**FX-S2, 'the amine-dependent MGO loss is a 1:1 amine-consuming adduct' (Kanzler 2026, one study):** "
          f"MGO lost {a['kanzler_one_to_one']['mgo_lost_mM']:.1f} mM by 120 min with Ala, at most "
          f"{a['kanzler_one_to_one']['alone_path_max_mM']:.1f} mM by the MGO-alone path, so >= "
          f"{a['kanzler_one_to_one']['amine_dependent_min_mM']:.1f} mM amine-dependent; a 1:1 adduct leaves <= "
          f"{a['kanzler_one_to_one']['ala_left_predicted_max_pct']:.1f} % Ala, observed "
          f"{a['kanzler_one_to_one']['ala_left_observed_pct']:.0f} % (SD used {a['kanzler_one_to_one']['sd_used']:.0f}); "
          f"t(2 df) {a['kanzler_one_to_one']['t_2df']:.1f}, p {a['kanzler_one_to_one']['p_two_sided']:.4f}; falsified: "
          f"**{a['kanzler_one_to_one']['falsified']}**. The Piloty-sized sink on the same course consumes "
          f"{a['kanzler_one_to_one']['consistent_check_piloty_sized_sink']['ala_consumed_pct_at_lower_bound']:.1f} % Ala "
          "(consistent with ~0 %).",
          "", f"**D2, Adams 2008 pH 4 vs 7 under R1:** {a['adams']['reading']} (AK flux ratio "
          f"{a['adams']['AK_flux_ratio_pH4_over_pH7']:.1f}).", "",
          "## Job (c): the designed decider under the adopted judge (synthetic)", "",
          f"Sets per truth {o['args']['n_sets']}, scoring draws {o['args']['n_draws']}, null replicates "
          f"{o['args']['n_pool']} per rival, adequacy replicates {o['args']['ppc_reps']}, per-reading CV "
          f"{o['args']['cv_dex']} dex. P(a) / P(b): P(credit the generating rival) with the gate passing; 95 % "
          "Clopper-Pearson in brackets; false credit: the worse direction. Gate: the M1.1 misspecified false credit "
          "k / n with the Wilson lower end (must not exceed 0.05), for the correctly specified generator (x0: a "
          "run-level reproducibility term cancels in a within-run unit) and an unmodelled per-point term of 0.5 / 1 x "
          "the declared per-point SD. A failing gate turns every credit of that pair into 'undecided'.",
          "", "| pair | design | n | P(a) [CP 95] | P(b) [CP 95] | false credit (CP upper) | gate k/n (Wilson low): x0 / x0.5 / x1 | not computable | null unresolved a / b |",
          "|---|---|---:|---|---|---|---|---:|---|"]
    for r in o["summary"]:
        gs = " / ".join(f"{g['k']}/{g['n']} ({g['wilson95'][0]:.3f}) {'pass' if g['passes'] else 'FAIL'}"
                        for _, g in sorted(r["gates"].items(), key=lambda kv: float(kv[0])))
        L.append(f"| {r['pair']} | {r['design']} | {r['n_rep']} | {r['P_a']:.2f} [{r['P_a_cp95'][0]:.2f}, "
                 f"{r['P_a_cp95'][1]:.2f}] | {r['P_b']:.2f} [{r['P_b_cp95'][0]:.2f}, {r['P_b_cp95'][1]:.2f}] | "
                 f"{r['false_credit']:.2f} ({r['false_credit_cp95_upper']:.2f}) | {gs} | {r['not_computable']:.2f} | "
                 f"{r['pool_unresolved']['a']:.3f} / {r['pool_unresolved']['b']:.3f} |")
    L += ["", "## Declared limits", ""] + [f"- {n}" for n in NOTES]
    return "\n".join(L)


NOTES = (
    "The rivals' structures are exact in the mini-model, so every power is an upper bound for real, misspecified data.",
    "MGO is held at its charge (the design measures it and forces the model with it); forced inputs are therefore not "
    "redrawn by the null (case_inputs: forced_inputs False). With a measured MGO course the null must redraw it.",
    "The spent rows enter only through R1's band (ks >= job (a)'s ks_min; kS near Hofmann's implied rate), not as "
    "conditioning rows of the density: a declared simplification of p(y_dec | y_spent).",
    "The heat-up is a measured design input (1 min, every rival alike); an unmeasured heat-up mimics R1's lag.",
    "The AK trap reading is (recovery x AK + blank) / reporting limit (0.5 uM), recovery 0.1 - 1 uncalibrated, blank "
    "0.2 - 1 x the limit: a non-detect is a blank-level reading, not a censored one.",
    "R1's exit (ks) and the ALA -> aminoacetone transfer (k2 +- 1 dex) are declared, unsourced widths; R3 and C0 "
    "predict the step units structurally (constant-free).",
    "Scoring ensembles come from a defensive proposal (0.9 prior + 0.1 edge-weighted; log p / q passed to the judge); "
    "unresolved bootstrap replicates and observed sets are re-estimated once on independent ensembles 3x larger.",
    "The precision requirement is the shared floor (ESS >= 100), declared; the misspecified generator's extra term is "
    "a multiple of the per-point SD because the pyrazine lane has no within-laboratory term.",
)


def run(args) -> Dict[str, Any]:
    c = checks()
    if args.checks or not c["pass"]:
        return {"checks": c}
    mem = memory_guard(not args.no_memory_guard)
    rng = np.random.default_rng(args.seed)
    job_a = {"onset": M.onset_falsification(rng, reps=args.r2_reps), "strecker_only_exit": M.strecker_only_amine_exit(),
             "kanzler_one_to_one": M.kanzler_one_to_one(),
             "adams": M.adams_ph_budget()}
    if args.merge:
        jc = {}
        for f in args.merge:
            jc.update(json.loads(Path(f).read_text()))
    else:
        jc = job_c(args)
    return {"checks": c, "memory": mem, "job_a": job_a, "job_c": jc, "summary": summarise(jc),
            "args": {k: (list(v) if isinstance(v, tuple) else v) for k, v in vars(args).items()}}


def _jsonable(x):
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (set, frozenset, tuple)):
        return list(x)
    return str(x)


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--checks", action="store_true")
    ap.add_argument("--n-sets", type=int, default=50)
    ap.add_argument("--n-draws", type=int, default=20000)
    ap.add_argument("--n-pool", type=int, default=1000)
    ap.add_argument("--ppc-reps", type=int, default=100)
    ap.add_argument("--cv-dex", type=float, default=0.08)
    ap.add_argument("--extra-factors", type=float, nargs="*", default=[0.5, 1.0])
    ap.add_argument("--r2-reps", type=int, default=20000)
    ap.add_argument("--parts", default="CS", help="C, S or CS")
    ap.add_argument("--only", default="", help="run only the keys containing this text")
    ap.add_argument("--pairs", nargs="*", default=[], help="run only these pairs, e.g. C1_vs_C0")
    ap.add_argument("--partial", default="", help="write job (c) to this JSON as it runs (no record)")
    ap.add_argument("--merge", nargs="*", default=[], help="assemble the record from --partial files")
    ap.add_argument("--seed", type=int, default=20261107)
    ap.add_argument("--no-memory-guard", action="store_true")
    ap.add_argument("--out-dir", default=str(ROOT / "results" / "validation"))
    ap.add_argument("--no-write", action="store_true")
    args = ap.parse_args(argv)
    o = run(args)
    if args.partial and not args.merge:
        return 0
    if args.checks or "job_c" not in o:
        print(json.dumps(o["checks"], default=_jsonable, indent=1))
        return 0 if o["checks"]["pass"] else 1
    payload = {"artifact": STEM,
               "provenance": provenance.provenance_block(
                   STEM, generated_by=f"scripts/generators/{Path(__file__).name}",
                   inputs=[ROOT / "src" / "m6_reduced_power.py", Path(__file__).resolve()]),
               "design": CASE + " secs. 4 - 5", "notes": list(NOTES), **o}
    md = render(payload)
    if args.no_write:
        print(md)
        return 0
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = out_dir / STEM
    stem.with_suffix(".json").write_text(provenance.dumps_compact(payload, default=_jsonable) + "\n")
    stem.with_suffix(".md").write_text(md + "\n")
    print(f"wrote {stem.with_suffix('.md')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
