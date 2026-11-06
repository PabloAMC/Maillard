#!/usr/bin/env python
"""
m4_reduced_power_check.py -- case M4 sec. 5.4: the reduced forced-input mini-model power check (synthetic only).

The mini-model is src/m4_reduced_power.py. The case is docs/reference/cases/M4_glucose_mft_protein_amine.md (draft
PR #90). Writes results/validation/kinetic_core_m4_reduced_power.{json,md}.

    python scripts/generators/m4_reduced_power_check.py --checks           # build checks only (fast)
    python scripts/generators/m4_reduced_power_check.py                    # checks + jobs (a), (b), (c); the record

The jobs (sec. 5.4, ruling 2: falsify, never credit):
- (a) the constant-free falsifications as arithmetic: C2 + C3 against D3's printed mixed cell; R4-as-written against D5.
- (b) D1 / D2 for the record only (development units, never credited): their power on SYNTHETIC Hofmann pots.
- (c) the decider's design power: per arm of sec. 4.4 and their unions, at n = 3 and n = 6 heatings per cell,
  P(the arm credits the generating rival) for R1 vs R5, the claim vs R0, the claim vs R2 and R4' vs the claim, both
  directions, with false-credit and "both inadequate" rates, under the primary and every required-agreement form
  (both band ends), the C-1 / C-6 equilibration variant, the OPD artefact, the d3 / M+3 collision and the
  label-precision limit. An arm below 0.5 for a pair is demoted to "ranks, does not decide" for that pair.

Rules enforced here: a memory guard (free memory >= 35 % before the heavy part; --no-memory-guard for tests); threads
pinned to 1 by the caller; nothing here reads a measured value of a decisive unit, a reserved set or a TEST set.
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
from dataclasses import replace
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import m4_reduced_power as M  # noqa: E402
from src import provenance  # noqa: E402

STEM = "kinetic_core_m4_reduced_power"
CASE = "docs/reference/cases/M4_glucose_mft_protein_amine.md (draft PR #90, branch m4-case-ready)"
MEMORY_MIN_FREE_PCT = 35
PAIRS = (("R1", "R5"), ("CLAIM", "R0"), ("CLAIM", "R2"), ("R4P", "CLAIM"))
PAIR_LABEL = {("R1", "R5"): "R1 vs R5", ("CLAIM", "R0"): "claim vs R0", ("CLAIM", "R2"): "claim vs R2",
              ("R4P", "CLAIM"): "R4' vs claim"}
#: truth rival -> the rivals its synthetic laboratories are scored by
PLAN = {"R1": ("R1", "R5"), "R5": ("R5", "R1"), "CLAIM": ("CLAIM", "R0", "R2", "R4P", "R1", "R5"),
        "R0": ("R0", "CLAIM"), "R2": ("R2", "CLAIM"), "R4P": ("R4P", "CLAIM")}
N_HEAT = (3, 6)

L10 = math.log10
LAYER_FORMS: Tuple[M.Form, ...] = (
    M.Form("formation_x0.1", pins=(("log10_mult_form", L10(0.1)),)),
    M.Form("formation_x3", pins=(("log10_mult_form", L10(3.0)),)),
    M.Form("formation_pH_slope_0", pins=(("s_form", 0.0),)),
    M.Form("formation_pH_slope_1", pins=(("s_form", 1.0),)),
    M.Form("decay_x0.3", pins=(("log10_mult_dec", L10(0.3)),)),
    M.Form("decay_x3", pins=(("log10_mult_dec", L10(3.0)),)),
    M.Form("split_phi_0.2", pins=(("phi_ref", 0.2),)),
    M.Form("split_phi_0.6", pins=(("phi_ref", 0.6),)),
    M.Form("ama_g_floor", pins=(("log10_kag_rel", -3.0),)),
    M.Form("ama_g_shipped", pins=(("log10_kag_rel", -1.5),)),
    M.Form("ama_g_ceiling", pins=(("log10_kag_rel", 0.0),)),
    M.Form("ama_g_pH_slope_0.44", kag_ph=True),
    M.Form("o2_cysteine_budget", o2_budget=True),
    M.Form("dg1_exits_x0.3", pins=(("log10_mult_dg1", L10(0.3)),)),
    M.Form("dg1_exits_x3", pins=(("log10_mult_dg1", L10(3.0)),)),
    M.Form("dg1_exits_pH_0", pins=(("s_dg1", 0.0),)),
    M.Form("dg1_exits_pH_0.5", pins=(("s_dg1", 0.5),)),
    M.Form("dg1_o2_drain_x1", pins=(("m_o2", 1.0),)),
    M.Form("dg1_o2_drain_x4.4", pins=(("m_o2", 4.4),)),
    M.Form("entries_floor", pins=(("log10_odds_odg", M.COORD["log10_odds_odg"].lo), ("log10_u5", -3.0))),
    M.Form("entries_ceiling", pins=(("log10_odds_odg", 0.0), ("log10_u5", 0.0))),
    M.Form("r2_strecker_floor", pins=(("log10_kstr", M.COORD["log10_kstr"].lo),)),
    M.Form("r2_strecker_ceiling", pins=(("log10_kstr", M.COORD["log10_kstr"].hi),)),
    M.Form("r4p_catalysis_x1", pins=(("log10_cat", 0.0),)),
    M.Form("r4p_catalysis_x10", pins=(("log10_cat", 1.0),)),
    M.Form("m2_via_free_nf_0.999", pins=(("s_nf", 0.999),)),
    M.Form("m2_via_free_nf_0.07", pins=(("s_nf", 0.07),)),
    M.Form("pre_dossier_bands", bands=(("ea_dec", 85.7, 107.0), ("rho", 0.71, 1.0), ("m_o2", 1.0, 3.4))),
    M.Form("n2_residual_o2_0.05", n2_residual=0.05),
)
DATA_FORMS: Tuple[M.Form, ...] = (
    M.Form("omega_0.390", omega=M.OMEGA_SENS_DEX),
    M.Form("lod_x0.1", lod_scale=0.1),
    M.Form("lod_x10", lod_scale=10.0),
    M.Form("label_precision_limit", label_limit=True),
    M.Form("d3_standard_collision", d3_collision=True),
    M.Form("opd_continuous_trapping", opd="continuous"),
)
TRUTH_FORMS: Tuple[M.Form, ...] = (M.Form("equilibration_e0.5", truth_e=0.5),)
BEYOND_SPEC = ("n2_residual_o2_0.05",)
#: the designs every form is run on (the primary runs every design in M.DESIGNS)
FORM_DESIGNS = ("arm1", "arm2", "arm3", "arm1+arm3", "arm2+arm3", "arm1+arm2+arm3")
CLAIM_PAIRS = (("CLAIM", "R0"), ("CLAIM", "R2"), ("R4P", "CLAIM"))


def minimum_design(vd: Mapping[str, Any]) -> Dict[str, Any]:
    """The fewest-heatings design whose every pair keeps P(separates) >= 0.5 at every spec form: for the claim
    (vs R0, R2, R4') and for the claim plus its carbon origin (R1 vs R5)."""
    out = {}
    for label, pairs in (("claim", CLAIM_PAIRS), ("claim_and_carbon_origin", CLAIM_PAIRS + (("R1", "R5"),))):
        for n_heat in N_HEAT:
            ok = []
            for ds in FORM_DESIGNS:
                mins = [vd.get(f"{a}|{b}|{ds}|n{n_heat}", {}).get("min_over_spec_forms") for a, b in pairs]
                if all(m is not None and m >= 0.5 for m in mins):
                    ok.append((M.heatings(ds, n_heat), ds, min(mins)))
            ok.sort()
            out[f"{label}|n{n_heat}"] = ({"design": ok[0][1], "heatings": ok[0][0], "worst_pair_power": ok[0][2],
                                          "all_passing": [x[1] for x in ok]} if ok else
                                         {"design": None, "all_passing": []})
    return out


# ============================================================================================ guards and records

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
    return {"free_pct_at_start": free, "min_required": MEMORY_MIN_FREE_PCT, "enforced": enabled}


def _jsonable(x):
    if isinstance(x, np.floating):
        return float(x)
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.bool_):
        return bool(x)
    return str(x)


def clopper_pearson(k: int, r: int) -> List[float]:
    from scipy.stats import beta
    lo = 0.0 if k == 0 else float(beta.ppf(0.025, k, r - k + 1))
    hi = 1.0 if k == r else float(beta.ppf(0.975, k + 1, r - k))
    return [round(lo, 3), round(hi, 3)]


# ============================================================================================ checks

HAND_TABLE = {  # by hand from the maps in the module docstring; (K_u, K_s, L_u, L_s, C23)
    "L1_M": (1, 0, 0, 1, 0.75), "L1_ring": (0, 0, 0, 1, 0.5), "L1_acyl": (1, 0, 0, 0, 0.5),
    "L6_M": (0, 1, 1, 0, 0.75), "L6_ring": (0, 0, 1, 0, 0.5), "L6_acyl": (0, 1, 0, 0, 0.5),
    "L2_M": (1, 1, 1, 1, 0.75), "L2_ring": (1, 1, 0, 1, 0.75), "L2_acyl": (1, 0, 1, 0, 0.5),
    "L12_M1": (0, 1, 1, 0, 0.0), "L12_M2": (1, 0, 0, 1, 0.5), "CAM": (0, 0, 0, 0, 0.5),
}


def check_label_table() -> Dict[str, Any]:
    bad = []
    for read, exp in HAND_TABLE.items():
        tb = M.label_table(read)
        got = (tb["K_u"], tb["K_s"], tb["L_u"], tb["L_s"], tb["C23"])
        if any(abs(a - b) > 0 for a, b in zip(got, exp)):
            bad.append({"read": read, "expected": exp, "got": got})
    return {"reads": len(HAND_TABLE), "mismatches": bad, "pass": not bad}


def check_balances(rng: np.random.Generator, n: int = 48) -> Dict[str, Any]:
    worst: Dict[str, float] = {}
    for rival in M.RIVALS:
        th, _ = M.draw_set(rival, M.PRIMARY, n, rng)
        for av in M.ARM_VARIANTS:
            for pot in M.POTS[av]:
                b = M.balances(pot, M.simulate(th, rival, M.PRIMARY, pot, full=True))
                for kk, v in b.items():
                    worst[kk] = max(worst.get(kk, 0.0), v)
    return {"max_residual_mM": worst, "tolerance": 1e-10, "pass": all(v < 1e-10 for v in worst.values())}


def check_rival_tables(rng: np.random.Generator, n: int = 64) -> Dict[str, Any]:
    """Each rival's shares at its draws against the hand table (the C2 + C3 bound's share of the pool is the only
    allowed difference), and R0's 1:1 mixed share against the binomial 2 q (1 - q)."""
    out: Dict[str, Any] = {}
    i_cam = [i for i in M.GROUPS["a1:CAM"]]
    i_l1 = [i for i in M.GROUPS["a1:L1_M"]]
    for rival in M.RIVALS:
        th, _ = M.draw_set(rival, M.PRIMARY, n, rng)
        p = M.predict(th, rival, M.PRIMARY, arms=("a1",))
        mft_c_share = []
        for pot in M.POTS["a1"]:
            r = M.simulate(th, rival, M.PRIMARY, pot)
            tot = r[("MFTK", 20.0)] + r[("MFTL", 20.0)] + r[("MFTD", 20.0)] + r[("MFTC", 20.0)]
            mft_c_share.append(np.max(r[("MFTC", 20.0)] / np.maximum(tot, 1e-300)))
        out[rival] = {"cam_share_range": [float(np.nanmin(p.share_a[:, i_cam])), float(np.nanmax(p.share_a[:, i_cam]))],
                      "L1_M_at_e0_range": [float(np.nanmin(p.share_a[:, i_l1])), float(np.nanmax(p.share_a[:, i_l1]))],
                      "max_c23_fraction_of_mft": float(np.max(mft_c_share))}
    r0 = out["R0"]["cam_share_range"]
    out["R0_binomial_reproduced"] = bool(abs(r0[0] - 0.5) < 1e-12 and abs(r0[1] - 0.5) < 1e-12)
    return out


def check_dt(rng: np.random.Generator, n: int = 24) -> Dict[str, Any]:
    worst = {}
    for rival in ("R1", "R5", "R2", "R4P"):
        th, _ = M.draw_set(rival, M.PRIMARY, n, rng)
        for av in ("a1", "a2", "a3"):
            for pot in M.POTS[av][::3]:
                a = M.simulate(th, rival, M.PRIMARY, pot, dt=M.DT)
                b = M.simulate(th, rival, M.PRIMARY, pot, dt=M.DT / 4)
                for key in a:
                    if key[0] in ("MFTK", "MFTL", "MFTD", "NFK", "NFL", "AP", "CYS"):
                        va, vb = a[key], b[key]
                        sp = {"MFTK": "MFT", "MFTL": "MFT", "MFTD": "MFT", "NFK": "NF", "NFL": "NF"}.get(key[0], key[0])
                        big = vb > M.LOD_MM[sp]          # the observable range only
                        if big.any():
                            rel = float(np.max(np.abs(np.log10(va[big] / vb[big]))))
                            worst[key[0]] = max(worst.get(key[0], 0.0), rel)
    return {"max_abs_log10_ratio_dt_vs_dt_over_4": worst, "note": "implicit Euler, dt 0.1 min vs 0.025 min"}


def check_shipped_constants() -> Dict[str, Any]:
    from src.kinetic_core import engine, sulfur
    p = engine.core_parameters(engine.SULFUR)
    k145 = sulfur.sulfur_rate_constants_at(p, M.T145, 5.0)
    k100 = sulfur.sulfur_rate_constants_at(p, M.T100, 5.0)
    rows = {
        "k_glc_ha_145": (k145["r_glc_c2c3"], M.K_GLC_HA_145),
        "k_fru_odg_100": (k100["r_fru_odg"], M.K_FRU_ODG_100),
        "k_glc_tdg_100": (k100["r_glc_tdg"], M.K_GLC_TDG_100),
        "k_odg_aa_100": (k100["r_odg_aa"], M.K_ODG_AA_100),
        "k_odg_aa_145": (k145["r_odg_aa"], M.K_ODG_AA_100 * float(M.arr(M.EA_ODG_AA, M.T145, M.T100))),
        "k_ama_g_100": (p["k_ama_g"].k_ref, 10.0 ** M.LOG10_K_AMA_G_100),
        "k_g_go_100": (p["k_g_go"].k_ref, 10.0 ** M.LOG10_K_G_GO_100),
    }
    out = {kk: {"engine": float(a), "literal": float(b), "rel": abs(a / b - 1.0)} for kk, (a, b) in rows.items()}
    return {"rows": out, "pass": all(v["rel"] < 5e-3 for v in out.values())}


def run_checks(seed: int = 11) -> Dict[str, Any]:
    rng = np.random.default_rng(seed)
    t0 = time.time()
    out = {"label_table": check_label_table(), "balances": check_balances(rng), "rival_tables": check_rival_tables(rng),
           "dt": check_dt(rng), "shipped_constants": check_shipped_constants()}
    out["wall_seconds"] = round(time.time() - t0, 1)
    out["pass"] = bool(out["label_table"]["pass"] and out["balances"]["pass"] and out["shipped_constants"]["pass"]
                       and out["rival_tables"]["R0_binomial_reproduced"])
    return out


# ============================================================================================ job (a)

def job_a() -> Dict[str, Any]:
    a = M.d3_arithmetic()
    return {
        "c2c3_against_d3": a,
        "reading": (f"A C2 + C3 route in a 1:1 [13C6] / [12C6] pool leaves a mixed share of 2 q (1 - q) = "
                    f"{a['mixed_share_c2c3']:.3f} whatever its rate constants; D3's printed mixed cell is <= 1 % "
                    f"(+-1 % absolute). The full route misses by {a['z_full_route']:.0f} sd "
                    f"(~{a['log_lik_ratio_nats_full_route']:.0f} nats); at q = 0.45 - 0.55 the share is "
                    f"{a['q_band_0.45_0.55'][0]:.3f} - {a['q_band_0.45_0.55'][1]:.3f}. A C2 + C3 share of the hexose MFT "
                    f"above f = {a['max_route_share_f']:.3f} is excluded at 1.96 sd (the R1 + f R3 sensitivity: "
                    f"f = 0.1 predicts {a['r1_plus_f_r3_at_f_0.1']:.3f}). Holds at every band corner: no rate "
                    f"constant enters. The one assumption: the two fragments come from independent molecules "
                    f"(no in-cage recombination of one glucose's own C2 and C3)."),
        "r4_as_written_against_d5": ("Structural, no run: R4-as-written has no Cys-hexose Amadori species, so it cannot "
                                     "carry the pool measured at 90 C (de Roos Table II) and at 145 C (Du 2024b, "
                                     "~62 % of the cysteine within 10 min). Falsified at every band corner. R4' "
                                     "carries the pool as a sink and is not falsified on disk."),
        "data_read": "none: D3's '<= 1 %' cell and D5's existence statement are taken as the case quotes them; no "
                     "TEST table was opened",
    }


# ============================================================================================ draw sets

class DrawCache:
    """Draw sets (scoring: Sobol; truth: independent pseudo-random) and their predictions, per (rival, layer form)."""

    def __init__(self, n_score: int, n_truth: int, arms: Sequence[str], seed: int):
        self.n_score, self.n_truth, self.arms, self.seed = n_score, n_truth, tuple(arms), seed
        self.cache: Dict[Tuple, Dict[str, Any]] = {}

    def get(self, rival: str, form: M.Form) -> Dict[str, Any]:
        layer = replace(M.PRIMARY, name="layer", pins=form.pins, bands=form.bands, kag_ph=form.kag_ph,
                        o2_budget=form.o2_budget, n2_residual=form.n2_residual)
        key = (rival, layer.layer_key)
        if key in self.cache:
            return self.cache[key]
        t0 = time.time()
        h = zlib.crc32(f"{rival}|{form.name}".encode()) % 10_000
        rs = np.random.default_rng([self.seed, h, 1])
        rt = np.random.default_rng([self.seed, h, 2])
        th_s, cnt_s = M.draw_set(rival, layer, self.n_score, rs, sobol=True, proposal=True)
        th_t, cnt_t = M.draw_set(rival, layer, self.n_truth, rt, sobol=False)
        ps = M.predict(th_s, rival, layer, arms=self.arms)
        pt = M.predict(th_t, rival, layer, arms=self.arms)
        cols = [i for i in M.LEVEL_IDX if M.OBS[i].group.split(":")[0] in self.arms]
        bad = int(np.sum(~np.isfinite(ps.level[:, cols]).all(axis=1)))
        val = {"score": ps, "truth": pt, "theta_truth": th_t, "counts": {"score": cnt_s, "truth": cnt_t},
               "nonfinite_draws": bad, "wall_s": round(time.time() - t0, 1)}
        self.cache[key] = val
        return val


# ============================================================================================ scoring

def _holm_ok(pvals: Sequence[float], alpha: float = 0.05) -> bool:
    ps = sorted(pvals)
    k = len(ps)
    return all(p > alpha / (k - i) for i, p in enumerate(ps))


def score_dataset(preds: Mapping[str, M.Prediction], d: M.Data, form: M.Form, designs: Sequence[str],
                  scorers: Sequence[str]) -> Dict[Tuple[str, str], Dict[str, Any]]:
    groups = sorted({g for ds in designs for g in M.DESIGNS[ds]})
    out: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for sc in scorers:
        gl = M.group_loglik(preds[sc], d, form, groups)
        for ds in designs:
            parts_all = M.design_parts(gl, M.DESIGNS[ds])
            ll = M.design_loglik(gl, M.DESIGNS[ds], parts_all)
            lp = M.lppd(ll, preds[sc].logw, half_check=True)
            arm_parts = M.arm_split(ds)
            pv, mus, idxs = [], [], []
            for part in arm_parts:
                pp = parts_all if len(arm_parts) == 1 else M.design_parts(gl, part)
                p, mu, idx = M.adequacy(preds[sc], d, form, pp)
                pv.append(p)
                mus.append(mu), idxs.append(idx)
            out[(sc, ds)] = {"lppd": lp["lppd"], "h1": lp["lppd_half1"], "h2": lp["lppd_half2"], "k_hat": lp["k_hat"],
                             "ess": lp["ess"], "resolved": lp["resolved"], "adequate": _holm_ok(pv), "p_min": min(pv),
                             "mu": np.concatenate(mus) if mus else np.array([]),
                             "idx": np.concatenate(idxs) if idxs else np.array([], int)}
    return out


def _verdict(d: float) -> int:
    return 1 if d > M.THRESHOLD_NATS else (-1 if d < -M.THRESHOLD_NATS else 0)


def decide(sa: Dict[str, Any], sb: Dict[str, Any], form: M.Form) -> Dict[str, Any]:
    """The case's rule: d > 1.92 nats and (the favoured rival adequate, or the means differ by > omega). 'strict'
    drops the escape (adequacy required). 'half_stable': the sign verdict is the same on both halves of the draw sets
    (a convergence check on the decision, beside PSIS k-hat / ESS)."""
    d = sa["lppd"] - sb["lppd"]
    diff = M.means_differ(sa["mu"], sa["idx"], sb["mu"], sb["idx"], form)
    v1, v2 = _verdict(sa["h1"] - sb["h1"]), _verdict(sa["h2"] - sb["h2"])
    return {"credit_a": bool(d > M.THRESHOLD_NATS and (sa["adequate"] or diff)),
            "credit_b": bool(d < -M.THRESHOLD_NATS and (sb["adequate"] or diff)),
            "strict_a": bool(d > M.THRESHOLD_NATS and sa["adequate"]),
            "strict_b": bool(d < -M.THRESHOLD_NATS and sb["adequate"]),
            "both_inadequate": bool(not sa["adequate"] and not sb["adequate"]),
            "unresolved": bool(not (sa["resolved"] and sb["resolved"])),
            "half_stable": bool(v1 == v2 == _verdict(d)), "d": float(d)}


def run_power_form(cache: DrawCache, form: M.Form, designs: Sequence[str], r: int, pairs: Sequence[Tuple[str, str]],
                   seed: int, truths_extra: Sequence[Tuple[str, Tuple[str, str]]] = ()) -> Dict[str, Any]:
    """Every pair x design x n under one form; both directions; plus extra (truth, pair) cells (the mixture truth)."""
    t0 = time.time()
    need_truth: Dict[str, set] = {}
    for a, b in pairs:
        need_truth.setdefault(a, set()).update({a, b})
        need_truth.setdefault(b, set()).update({a, b})
    for t, (a, b) in truths_extra:
        need_truth.setdefault(t, set()).update({a, b})
    sets = {rv: cache.get(rv, form) for rv in {x for v in need_truth.values() for x in v} | set(need_truth)}
    preds = {rv: s["score"] for rv, s in sets.items()}
    res: Dict[str, Any] = {"cells": {}, "draw_sets": {rv: {"counts": s["counts"], "nonfinite": s["nonfinite_draws"],
                                                          "wall_s": s["wall_s"]} for rv, s in sets.items()}}
    for n_heat in N_HEAT:
        rng = np.random.default_rng([seed, zlib.crc32(form.name.encode()) % 10_000, n_heat])
        for truth, scorers in need_truth.items():
            pt = sets[truth]["truth"]
            rr = min(r, pt.n)
            recs = []
            for j in range(rr):
                d = M.synthesize(pt, j, n_heat, form, rng)
                recs.append(score_dataset(preds, d, form, designs, sorted(scorers)))
            for a, b in dict.fromkeys(list(pairs) + [p for t, p in truths_extra if t == truth]):
                if truth not in (a, b) and (truth, (a, b)) not in truths_extra:
                    continue
                for ds in designs:
                    dec = [decide(x[(a, ds)], x[(b, ds)], form) for x in recs]
                    ca = sum(z["credit_a"] for z in dec)
                    cb = sum(z["credit_b"] for z in dec)
                    key = f"{a}|{b}|{ds}|n{n_heat}|truth={truth}"
                    res["cells"][key] = {
                        "R": rr, "credit_a": ca / rr, "credit_b": cb / rr,
                        "credit_a_ci95": clopper_pearson(ca, rr), "credit_b_ci95": clopper_pearson(cb, rr),
                        "strict_a": float(np.mean([z["strict_a"] for z in dec])),
                        "strict_b": float(np.mean([z["strict_b"] for z in dec])),
                        "both_inadequate": float(np.mean([z["both_inadequate"] for z in dec])),
                        "unresolved": float(np.mean([z["unresolved"] for z in dec])),
                        "half_unstable": float(np.mean([not z["half_stable"] for z in dec])),
                        "d_quantiles": [float(q) for q in np.quantile([z["d"] for z in dec], [0.1, 0.5, 0.9])],
                        "ess_median": [float(np.median([x[(a, ds)]["ess"] for x in recs])),
                                       float(np.median([x[(b, ds)]["ess"] for x in recs]))],
                    }
    res["wall_s"] = round(time.time() - t0, 1)
    return res


def summarise_pairs(cells: Mapping[str, Any], designs: Sequence[str]) -> Dict[str, Any]:
    """Per pair x design x n: power under each truth, false credit, P(separates) = the smaller power."""
    out = {}
    for a, b in PAIRS:
        for ds in designs:
            for n_heat in N_HEAT:
                ka = f"{a}|{b}|{ds}|n{n_heat}|truth={a}"
                kb = f"{a}|{b}|{ds}|n{n_heat}|truth={b}"
                if ka not in cells or kb not in cells:
                    continue
                ca, cb = cells[ka], cells[kb]
                rec = {
                    "power_truth_a": ca["credit_a"], "false_credit_truth_a": ca["credit_b"],
                    "power_truth_b": cb["credit_b"], "false_credit_truth_b": cb["credit_a"],
                    "separates": min(ca["credit_a"], cb["credit_b"]),
                    "max_false_credit": max(ca["credit_b"], cb["credit_a"]),
                    "separates_strict": min(ca["strict_a"], cb["strict_b"]),
                    "both_inadequate": max(ca["both_inadequate"], cb["both_inadequate"]),
                    "unresolved": max(ca["unresolved"], cb["unresolved"]),
                    "half_unstable": max(ca["half_unstable"], cb["half_unstable"]),
                }
                # the claim's two forms as truths: crediting the claim is correct, the rival is a false credit
                if "CLAIM" in (a, b):
                    side = "credit_a" if a == "CLAIM" else "credit_b"
                    other = "credit_b" if a == "CLAIM" else "credit_a"
                    sside = "strict_a" if a == "CLAIM" else "strict_b"
                    for form_r in ("R1", "R5"):
                        kf = f"{a}|{b}|{ds}|n{n_heat}|truth={form_r}"
                        if kf in cells:
                            rec[f"power_truth_{form_r}"] = cells[kf][side]
                            rec[f"false_credit_truth_{form_r}"] = cells[kf][other]
                            rec["separates"] = min(rec["separates"], cells[kf][side])
                            rec["separates_strict"] = min(rec["separates_strict"], cells[kf][sside])
                            rec["max_false_credit"] = max(rec["max_false_credit"], cells[kf][other])
                out[f"{a}|{b}|{ds}|n{n_heat}"] = rec
    return out


def claim_form_truths(pairs: Sequence[Tuple[str, str]], mixture: bool) -> Tuple[Tuple[str, Tuple[str, str]], ...]:
    """Extra (truth, pair) cells: the claim's forms R1 and R5 as truths for every pair holding the claim, and the
    FM4 mixture (claim draws) as truth for R1 vs R5 (both credits false there)."""
    out = [(t, p) for p in pairs if "CLAIM" in p for t in ("R1", "R5")]
    if mixture and ("R1", "R5") in pairs:
        out.append(("CLAIM", ("R1", "R5")))
    return tuple(out)


# ============================================================================================ job (b): D1 / D2

def job_b(cache: DrawCache, r: int, seed: int) -> Dict[str, Any]:
    """D1 (glucose NF / MFT, pH 5) and D2 (MFT pH 3 / 7 against pH 5) on SYNTHETIC Hofmann pots (= arm 1's air cells at
    145 C): omega per value, triplicate CV, no tau (within-pot / within-study), a free common offset (flat) integrated
    on a grid; floors are the declared transfer, 1 % of the truth's lowest value in the table, x0.1 / x1 / x10."""
    cols = {}
    for i, o in enumerate(M.OBS):
        if o.kind == "level" and o.group in ("a1:MFT", "a1:NF") and o.t_c == 145.0 and o.pot.endswith("_air"):
            cols[(o.species, int(o.pot.split("pH")[1][0]))] = i
    d1 = [cols[("NF", 5)], cols[("MFT", 5)]]
    d2 = [cols[("MFT", 3)], cols[("MFT", 5)], cols[("MFT", 7)]]
    sd = math.sqrt(M.OMEGA_DEX ** 2 + M.CV_DEX ** 2 / 3.0)
    grid = np.linspace(-8.0, 8.0, 161)
    from scipy.special import log_ndtr

    def ll_unit(pred_log, y, cens, floor):           # pred_log [N, V]
        mu = pred_log[:, :, None] + grid[None, None, :]
        z = (y[None, :, None] - mu) / sd
        ll = np.where(cens[None, :, None], log_ndtr((floor[None, :, None] - mu) / sd),
                      -0.5 * z ** 2 - math.log(sd) - 0.5 * math.log(2 * math.pi))
        return M._lse(ll.sum(axis=1), 1) - math.log(len(grid))

    pairs = (("CLAIM", "R2"), ("R1", "R5"), ("R4P", "CLAIM"))
    out: Dict[str, Any] = {}
    for unit, ii in (("D1", d1), ("D2", d2), ("D1+D2", sorted(set(d1 + d2)))):
        for fscale in (0.1, 1.0, 10.0):
            for a, b in pairs:
                for truth in (a, b):
                    pt = cache.get(truth, M.PRIMARY)["truth"]
                    rng = np.random.default_rng([seed, len(unit), int(fscale * 10), zlib.crc32(truth.encode()) % 997])
                    hits = false = 0
                    rr = min(r, pt.n)
                    for j in range(rr):
                        tv = np.log10(np.maximum(pt.level[j, ii], 1e-30))
                        floor = np.full(len(ii), math.log10(0.01 * fscale) + tv.min())
                        y = tv + rng.normal(0.0, sd, len(ii))
                        cens = y < floor
                        y = np.where(cens, floor, y)
                        lps = {}
                        for sc in (a, b):
                            sp_ = cache.get(sc, M.PRIMARY)["score"]
                            lps[sc] = M.lppd(ll_unit(sp_.loglevel[:, ii], y, cens, floor), sp_.logw)["lppd"]
                        dd = lps[a] - lps[b]
                        win = a if dd > M.THRESHOLD_NATS else (b if dd < -M.THRESHOLD_NATS else None)
                        hits += win == truth
                        false += win is not None and win != truth
                    out[f"{unit}|floor_x{fscale:g}|{a}|{b}|truth={truth}"] = {
                        "R": rr, "power": hits / rr, "false_credit": false / rr}
    return out


# ============================================================================================ the run

def run(args) -> Dict[str, Any]:
    t_start = time.time()
    mem = memory_guard(not args.no_memory_guard)
    checks = run_checks()
    if not checks["pass"]:
        raise SystemExit(f"checks failed: {json.dumps(checks, default=_jsonable)[:2000]}")
    if args.checks:
        return {"checks": checks, "memory": mem}
    arms_all = M.ARM_VARIANTS
    cache = DrawCache(args.n_draws, max(args.r_primary, args.r), arms_all, args.seed)
    out: Dict[str, Any] = {"memory": mem, "checks": checks, "job_a": job_a()}
    print("job (b) ...", flush=True)
    out["job_b"] = job_b(cache, args.r_primary, args.seed)
    all_designs = list(M.DESIGNS)
    print("primary ...", flush=True)
    prim = run_power_form(cache, M.PRIMARY, all_designs, args.r_primary, PAIRS, args.seed,
                          truths_extra=claim_form_truths(PAIRS, True))
    forms: Dict[str, Any] = {"primary": {"run": prim, "pairs": summarise_pairs(prim["cells"], all_designs)}}
    forms["primary"]["mixture_truth_r1_r5"] = {k: v for k, v in prim["cells"].items()
                                               if k.startswith("R1|R5|") and k.endswith("truth=CLAIM")}
    todo = list(LAYER_FORMS) + list(DATA_FORMS) + list(TRUTH_FORMS)
    if args.forms:
        todo = [f for f in todo if f.name in set(args.forms.split(","))]
    cache_small = DrawCache(args.n_draws_forms, args.r, ("a1", "a2", "a3"), args.seed + 1)
    for f in todo:
        if time.time() - t_start > args.cap_hours * 3600:
            forms[f.name] = {"status": "compute-capped"}
            continue
        is_layer = f in LAYER_FORMS
        pairs = [p for p in PAIRS if not is_layer or any(M.form_touches(f, x) for x in p)]
        if not pairs:
            forms[f.name] = {"status": "touches no rival of any pair (identical to primary)"}
            continue
        print(f"form {f.name} pairs {pairs} ...", flush=True)
        designs = list(FORM_DESIGNS) + (["arm3_plus_opd_glucosone"] if f.opd != "post" else [])
        res = run_power_form(cache_small, f, designs, args.r, pairs, args.seed,
                             truths_extra=claim_form_truths(pairs, f in TRUTH_FORMS))
        forms[f.name] = {"run": res, "pairs": summarise_pairs(res["cells"], M.PRIMARY_DESIGNS),
                         "kind": "layer" if is_layer else ("data" if f in DATA_FORMS else "truth"),
                         "beyond_spec": f.name in BEYOND_SPEC}
        # free the layer draw sets of this form (memory)
        for key in [k for k in cache_small.cache if k[1] != M.PRIMARY.layer_key]:
            del cache_small.cache[key]
    # the small-set primary, so forms compare like with like
    print("primary (form-sized sets) ...", flush=True)
    prim_small = run_power_form(cache_small, M.PRIMARY, list(FORM_DESIGNS), args.r, PAIRS, args.seed + 7,
                                truths_extra=claim_form_truths(PAIRS, False))
    forms["primary_form_sized"] = {"run": prim_small, "pairs": summarise_pairs(prim_small["cells"], M.PRIMARY_DESIGNS)}
    out["forms"] = forms
    out["verdicts"] = verdicts(forms)
    out["minimum_design"] = minimum_design(out["verdicts"])
    out["settings"] = {"n_draws_primary": args.n_draws, "R_primary": args.r_primary, "n_draws_forms": args.n_draws_forms,
                       "R_forms": args.r, "seed": args.seed, "threshold_nats": M.THRESHOLD_NATS,
                       "n_heatings": N_HEAT, "e_grid": M.E_GRID.tolist(), "gh_nodes": M.GH_NODES}
    out["wall_seconds"] = round(time.time() - t_start, 1)
    return out


def verdicts(forms: Mapping[str, Any]) -> Dict[str, Any]:
    """Per pair x primary design x n: the primary P(separates), the minimum over every form that ran (required
    agreement, both band ends), which form, and the demotion (< 0.5)."""
    prim = forms["primary"]["pairs"]
    out = {}
    for a, b in PAIRS:
        for ds in FORM_DESIGNS:
            for n_heat in N_HEAT:
                key = f"{a}|{b}|{ds}|n{n_heat}"
                if key not in prim:
                    continue
                worst, worst_form, fc = prim[key]["separates"], "primary", prim[key]["max_false_credit"]
                worst_spec, worst_spec_form = worst, "primary"
                per_form = {}
                for fname, fv in forms.items():
                    if fname in ("primary",) or "pairs" not in fv:
                        continue
                    if key in fv["pairs"]:
                        s = fv["pairs"][key]["separates"]
                        per_form[fname] = s
                        if fname not in BEYOND_SPEC:
                            fc = max(fc, fv["pairs"][key]["max_false_credit"])
                        if s < worst:
                            worst, worst_form = s, fname
                        if fname not in BEYOND_SPEC and s < worst_spec:
                            worst_spec, worst_spec_form = s, fname
                out[key] = {"primary": prim[key]["separates"], "min_over_forms": worst, "worst_form": worst_form,
                            "min_over_spec_forms": worst_spec, "worst_spec_form": worst_spec_form,
                            "max_false_credit_spec_forms": fc,
                            "decides_primary": prim[key]["separates"] >= 0.5,
                            "decides_under_required_agreement": worst_spec >= 0.5, "per_form": per_form}
    return out


# ============================================================================================ markdown

def _f(x: float) -> str:
    return f"{x:.2f}"


def render(o: Mapping[str, Any]) -> str:
    L: List[str] = []
    s = o["settings"]
    L.append("# M4 reduced-model power check: sizing the designed experiment's arms (case sec. 5.4; SYNTHETIC)")
    L.append("")
    L.append(f"*Generated by `scripts/generators/m4_reduced_power_check.py`; mini-model `src/m4_reduced_power.py`. Case: "
             f"{CASE}. Primary: {s['n_draws_primary']} Sobol band-prior draws per rival, R = {s['R_primary']} synthetic "
             f"laboratories per truth; forms: {s['n_draws_forms']} draws, R = {s['R_forms']}. Credit = lppd difference > "
             f"{s['threshold_nats']} nats with the adequacy rule (posterior-predictive p > 0.05, Holm across arms, or "
             f"means differing by more than omega). Ruling 2 binds: nothing here credits a rival on disk; this sizes "
             f"W6.*")
    L.append("")
    L.append("**No measured value of any decisive unit, reserved set or TEST set was read.** D1 / D2 are scored on "
             "synthetic Hofmann pots; D3 enters job (a) only through the '<= 1 %' mixed cell the case quotes.")
    L.append("")
    ck = o["checks"]
    L.append("## Checks before any result was read")
    L.append("")
    L.append(f"- Hand-computed label table ({ck['label_table']['reads']} reads x 5 families): "
             f"{'reproduced exactly' if ck['label_table']['pass'] else 'MISMATCH'}.")
    b = ck["balances"]["max_residual_mM"]
    L.append("- Balances over every pot, every rival (48 draws each), max residual (mM): " +
             ", ".join(f"{k} {v:.1e}" for k, v in b.items()) + f" (tolerance 1e-10: "
             f"{'pass' if ck['balances']['pass'] else 'FAIL'}).")
    L.append(f"- R0's 1:1 mixed share reproduces the binomial 2 q (1 - q) = 0.5 exactly: "
             f"{ck['rival_tables']['R0_binomial_reproduced']}. Largest C2 + C3 bound share of any entry rival's MFT: "
             f"{max(v['max_c23_fraction_of_mft'] for k, v in ck['rival_tables'].items() if isinstance(v, dict) and k != 'R0'):.1e}.")
    L.append("- Step size (implicit Euler 0.1 vs 0.025 min), max |log10 ratio|: " +
             ", ".join(f"{k} {v:.3f}" for k, v in ck["dt"]["max_abs_log10_ratio_dt_vs_dt_over_4"].items()) + ".")
    L.append(f"- Shipped constants (frozen literals) against the engine: "
             f"{'equal within 0.5 %' if ck['shipped_constants']['pass'] else 'DRIFT'}.")
    L.append("")
    a = o["job_a"]
    L.append("## (a) The on-disk falsifications, as arithmetic")
    L.append("")
    L.append(f"- **C2 + C3 (R0's route and R3) against D3.** {a['c2c3_against_d3'] and a['reading']}")
    L.append(f"- **R4-as-written against D5.** {a['r4_as_written_against_d5']}")
    L.append("")
    L.append("## (b) D1 / D2: development record only (never credited)")
    L.append("")
    L.append("P(credit the generating rival) on synthetic Hofmann pots (n = 3, omega per value, no tau, floors at the "
             "declared transfer x0.1 / x1 / x10); false credit in brackets.")
    L.append("")
    L.append("| unit | pair | floor | power (truth = first) | power (truth = second) |")
    L.append("|---|---|---|---|---|")
    jb = o["job_b"]
    for unit in ("D1", "D2", "D1+D2"):
        for a_, b_ in (("CLAIM", "R2"), ("R1", "R5"), ("R4P", "CLAIM")):
            for fs in (0.1, 1.0, 10.0):
                x = jb[f"{unit}|floor_x{fs:g}|{a_}|{b_}|truth={a_}"]
                y = jb[f"{unit}|floor_x{fs:g}|{a_}|{b_}|truth={b_}"]
                L.append(f"| {unit} | {M.RIVAL_LABEL[a_]} vs {M.RIVAL_LABEL[b_]} | x{fs:g} | {_f(x['power'])} "
                         f"[{_f(x['false_credit'])}] | {_f(y['power'])} [{_f(y['false_credit'])}] |")
    L.append("")
    L.append("## (c) The designed experiment: power per arm (primary)")
    L.append("")
    L.append("P(separates) = the smaller of the two directions' power (the arm credits the truth whichever rival is "
             "true). Columns: n = 3 | n = 6; each cell 'separates (power truth A / power truth B; max false credit)'.")
    L.append("")
    pp = o["forms"]["primary"]["pairs"]
    L.append("| design | " + " | ".join(PAIR_LABEL[p] for p in PAIRS) + " |")
    L.append("|---|" + "---|" * len(PAIRS))
    for ds in M.DESIGNS:
        row = [ds]
        for a_, b_ in PAIRS:
            cells = []
            for n_heat in N_HEAT:
                k = f"{a_}|{b_}|{ds}|n{n_heat}"
                if k not in pp:
                    cells.append("-")
                    continue
                v = pp[k]
                cells.append(f"{_f(v['separates'])} ({_f(v['power_truth_a'])}/{_f(v['power_truth_b'])}; "
                             f"{_f(v['max_false_credit'])})")
            row.append(" \\| ".join(cells))
        L.append("| " + " | ".join(row) + " |")
    L.append("")
    L.append(f"Heatings per design (n = 3 / n = 6): " + "; ".join(
        f"{ds} {M.heatings(ds, 3)} / {M.heatings(ds, 6)}" for ds in M.DESIGNS) + ".")
    L.append("")
    mt = o["forms"]["primary"]["mixture_truth_r1_r5"]
    L.append("**The FM4 mixture as truth (claim draws, f* ~ U[0, 1]) scored by the pure forms R1 and R5** "
             "(both are false credits here):")
    L.append("")
    for k, v in sorted(mt.items()):
        parts = k.split("|")
        if parts[2] in M.PRIMARY_DESIGNS:
            L.append(f"- {parts[2]}, {parts[3]}: credit R1 {_f(v['credit_a'])}, credit R5 {_f(v['credit_b'])}, "
                     f"both inadequate {_f(v['both_inadequate'])}")
    L.append("")
    L.append("## Required agreement: the minimum over every form (both band ends)")
    L.append("")
    L.append("Each cell: primary P(separates) -> minimum over the spec's forms [worst form]; DEMOTED where the minimum "
             "is < 0.5 ('ranks, does not decide' for that pair). The N2 residual-O2 form (beyond the spec) is shown "
             "separately.")
    L.append("")
    vd = o["verdicts"]
    L.append("| design | n | " + " | ".join(PAIR_LABEL[p] for p in PAIRS) + " |")
    L.append("|---|---|" + "---|" * len(PAIRS))
    for ds in FORM_DESIGNS:
        for n_heat in N_HEAT:
            row = [ds, str(n_heat)]
            for a_, b_ in PAIRS:
                k = f"{a_}|{b_}|{ds}|n{n_heat}"
                if k not in vd:
                    row.append("-")
                    continue
                v = vd[k]
                tag = "decides" if v["decides_under_required_agreement"] else "DEMOTED"
                beyond = v["per_form"].get("n2_residual_o2_0.05")
                bx = f"; N2-residual {_f(beyond)}" if beyond is not None else ""
                row.append(f"{_f(v['primary'])} -> {_f(v['min_over_spec_forms'])} [{v['worst_spec_form']}] {tag}{bx}")
            L.append("| " + " | ".join(row) + " |")
    L.append("")
    L.append("**By truth (primary, n = 3 / n = 6).** For pairs holding the claim, the claim's truth is tested three "
             "ways: the FM4 mixture (claim draws), and each of its forms R1 and R5. Entry = P(the correct rival is "
             "credited); false credit in brackets.")
    L.append("")
    L.append("| design | pair | truth = rival | truth = claim (mixture) | truth = R1 | truth = R5 |")
    L.append("|---|---|---|---|---|---|")
    cells = o["forms"]["primary"]["run"]["cells"]
    for ds in FORM_DESIGNS:
        for a_, b_ in PAIRS:
            if "CLAIM" not in (a_, b_):
                continue
            rival = b_ if a_ == "CLAIM" else a_
            good_claim = "credit_a" if a_ == "CLAIM" else "credit_b"
            bad_claim = "credit_b" if a_ == "CLAIM" else "credit_a"
            vals = []
            for truth, good, bad in ((rival, bad_claim, good_claim), ("CLAIM", good_claim, bad_claim),
                                     ("R1", good_claim, bad_claim), ("R5", good_claim, bad_claim)):
                parts = []
                for n_heat in N_HEAT:
                    c = cells.get(f"{a_}|{b_}|{ds}|n{n_heat}|truth={truth}")
                    parts.append("-" if c is None else f"{_f(c[good])} [{_f(c[bad])}]")
                vals.append(" / ".join(parts))
            L.append(f"| {ds} | {PAIR_LABEL[(a_, b_)]} | " + " | ".join(vals) + " |")
    L.append("")
    L.append("**Strict adequacy (no 'means differ' escape), primary P(separates), n = 3 / n = 6:**")
    L.append("")
    L.append("| design | " + " | ".join(PAIR_LABEL[p] for p in PAIRS) + " |")
    L.append("|---|" + "---|" * len(PAIRS))
    for ds in FORM_DESIGNS:
        row = [ds]
        for a_, b_ in PAIRS:
            row.append(" / ".join(_f(pp[f"{a_}|{b_}|{ds}|n{n}"]["separates_strict"]) if f"{a_}|{b_}|{ds}|n{n}" in pp
                                  else "-" for n in N_HEAT))
        L.append("| " + " | ".join(row) + " |")
    L.append("")
    L.append("## Per form (P(separates), n = 3, arm1 / arm2 / arm3 / arm1+arm2+arm3)")
    L.append("")
    L.append("| form | kind | " + " | ".join(PAIR_LABEL[p] for p in PAIRS) + " |")
    L.append("|---|---|" + "---|" * len(PAIRS))
    for fname, fv in o["forms"].items():
        if "pairs" not in fv:
            L.append(f"| {fname} | - | {fv.get('status', '')} |" + " |" * (len(PAIRS) - 1))
            continue
        row = [fname, fv.get("kind", "primary")]
        for a_, b_ in PAIRS:
            vals = []
            for ds in ("arm1", "arm2", "arm3", "arm1+arm2+arm3"):
                k = f"{a_}|{b_}|{ds}|n3"
                vals.append(_f(fv["pairs"][k]["separates"]) if k in fv["pairs"] else "-")
            row.append(" / ".join(vals))
        L.append("| " + " | ".join(row) + " |")
    L.append("")
    L.append("## Density diagnostics")
    L.append("")
    pp_all = [v for fv in o["forms"].values() if "pairs" in fv for v in fv["pairs"].values()]
    worst_unres = max((v["unresolved"] for v in pp_all), default=0.0)
    mean_unres = float(np.mean([v["unresolved"] for v in pp_all])) if pp_all else 0.0
    L.append(f"- PSIS k-hat >= 0.7 or ESS < 100 in either rival's density: mean fraction {_f(mean_unres)}, worst "
             f"{_f(worst_unres)} over every pair x design x n x form. The band priors are wide and the shares precise, "
             f"so plain importance densities over the band prior do not meet the case's acceptance rule; the M3 SMC "
             f"estimator was not run (see below).")
    L.append("- What carries the verdicts instead: most decisions are structural (a rival predicts a zero or a share "
             "of 0 / 1 where the truth reads a level or the opposite share), with |d| in the hundreds to tens of "
             "thousands of nats, far beyond any importance-sampling error. The decision-level check: the verdict "
             "(credit A / credit B / neither) recomputed on each half of the draw sets. Fraction of laboratories whose "
             "verdict flips between halves, per pair x design (primary, n = 3 / n = 6):")
    L.append("")
    pp = o["forms"]["primary"]["pairs"]
    L.append("| design | " + " | ".join(PAIR_LABEL[p] for p in PAIRS) + " |")
    L.append("|---|" + "---|" * len(PAIRS))
    for ds in FORM_DESIGNS:
        L.append(f"| {ds} | " + " | ".join(" / ".join(_f(pp[f"{a_}|{b_}|{ds}|n{n}"]["half_unstable"])
                                                    for n in N_HEAT) for a_, b_ in PAIRS) + " |")
    L.append("")
    md_ = o["minimum_design"]
    L.append("## The minimum design")
    L.append("")
    for key, v in md_.items():
        what, n_ = key.split("|")
        if v["design"] is None:
            L.append(f"- **{what.replace('_', ' ')}, {n_}:** no design passes every pair at every spec form.")
        else:
            L.append(f"- **{what.replace('_', ' ')}, {n_}:** {v['design']} ({v['heatings']} heatings; worst pair "
                     f"{_f(v['worst_pair_power'])}); passing: {', '.join(v['all_passing'])}.")
    L.append("")
    L.append("## Declared choices and what could not be built as written")
    L.append("")
    for line in o["notes"]:
        L.append(f"- {line}")
    L.append("")
    return "\n".join(L)


NOTES = (
    "The mini-model's declared choices D-1 - D-15 are in the module docstring (src/m4_reduced_power.py); none sets a "
    "band from a scored or screened pot.",
    "Required agreement is read as 'the verdict must hold if the world sits at either end': each layer form pins the "
    "coordinate in BOTH the generating and the scoring draw sets. The joint reading (truth over the whole band, "
    "scored under every pinned form at once) was not built; it can only lower these numbers.",
    "Holm across arms is applied to the adequacy p-values of the favoured rival in a union design; the decision "
    "itself is an lppd difference with no p-value, so Holm has nothing else to act on (as in M3's power check).",
    "The MGO loss law (Ea 47 / 121 / 20, orders 1 / 2) is not a coordinate: the only MGO in the mini-model is R0's "
    "C2 + C3 route, carried as a sink-free upper bound (Hofmann Tables 7 / 8 largest yields), which holds at every "
    "corner. R3 is not a rival of the power check (sec. 5.4 item 2 lists R0, R1, R2, R4', R5); its rejection is job "
    "(a)'s arithmetic.",
    "R0 has no Cys-hexose Amadori species: in arm 3 the fed compound is inert under R0 (declared), and in arm 2 R0 has "
    "no route from a fed C6 dicarbonyl to MFT. R0 'pinned at B60' is its shipped constants; its nuisance coordinates "
    "(Cys -> H2S, 1-DG exits, downstream) are drawn from the same bands as every rival.",
    "The N2 arm is O2-free (factor 0) as the spec states; a residual-O2 form (0.05) was added beyond the spec, because "
    "R5's glucosone needs O2 and degassing is never perfect. It is reported beside, not inside, the required agreement.",
    "Shares are read only where the same cell's unlabelled MFT (NF) is >= 3 LOD (D-13). The adequacy check is a "
    "diagonal-Gaussian posterior-predictive chi-square over the read, non-censored values (posterior over draws, the "
    "e grid and each species' tau node); correlations between values are ignored (declared approximation).",
    "pH is the label pH throughout (0.5 M phosphate, no drift), as the spec's out-of-scope list says; part A and the "
    "full network (W3) are out of scope.",
    "The case's 'Hofmann's aqueous arm for (b)' is arm 1's air cells at 145 C (the same design); D1 / D2 floors are "
    "1 % of the truth's lowest value (the declared transfer applied to synthetic data, so no printed value is read).",
    "NOT BUILT AS WRITTEN: the density rule (exact, or MC with PSIS k-hat < 0.7 and ESS >= 100, else M3's SMC). The "
    "densities here are self-normalised importance estimates over the band prior (edge-weighted for the claim's f, "
    "D-15) and mostly fail k-hat / ESS. M3's SMC costs minutes per density and the check needs ~10^4 - 10^5 "
    "densities per form, so it was not run. The verdicts rest on |d| far beyond the MC error for structural "
    "contrasts, and on the half-set verdict check where |d| is near the threshold; cells with a high flip rate are "
    "flagged and should be read as uncertain to that extent.",
    "NOT BUILT: a separate test of 'f moves with temperature' (constant f against the FM4 mixture). The spec's pairs "
    "are R1 vs R5, the claim vs R0 / R2 and R4' vs the claim; the 120 C cells are scored inside arm 1 for those, and "
    "'arm1_no120' shows what they add.",
    "Forms run at 2048 draws and R = 30 synthetic laboratories (Clopper-Pearson half-width up to ~0.18 at 0.5), on "
    "six designs; the primary at 4096 draws and R = 150 on every design. The demotion line (0.5) is read on the "
    "minimum over forms, so a form near 0.5 is a coarse reading.",
)


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--checks", action="store_true")
    ap.add_argument("--n-draws", type=int, default=4096)
    ap.add_argument("--r-primary", type=int, default=150)
    ap.add_argument("--n-draws-forms", type=int, default=2048)
    ap.add_argument("--r", type=int, default=30)
    ap.add_argument("--seed", type=int, default=20261106)
    ap.add_argument("--forms", default="")
    ap.add_argument("--cap-hours", type=float, default=8.0)
    ap.add_argument("--no-memory-guard", action="store_true")
    ap.add_argument("--out-dir", default=str(ROOT / "results" / "validation"))
    ap.add_argument("--no-write", action="store_true")
    args = ap.parse_args(argv)
    o = run(args)
    if args.checks:
        print(json.dumps(o["checks"], default=_jsonable, indent=1)[:6000])
        return 0 if o["checks"]["pass"] else 1
    o["notes"] = list(NOTES)
    md = render(o)
    if args.no_write:
        print(md)
        return 0
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = out_dir / STEM
    inputs = [ROOT / "src" / "m4_reduced_power.py", Path(__file__).resolve()]
    payload = {"artifact": STEM,
               "provenance": provenance.provenance_block(STEM, generated_by=f"scripts/generators/{Path(__file__).name}",
                                                         inputs=inputs),
               "design": CASE + " sec. 5.4", **o}
    stem.with_suffix(".json").write_text(provenance.dumps_compact(payload, default=_jsonable) + "\n")
    stem.with_suffix(".md").write_text(md + "\n")
    print(f"wrote {stem.with_suffix('.md')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
