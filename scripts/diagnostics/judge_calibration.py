#!/usr/bin/env python
"""
judge_calibration.py -- why the decision procedure cannot credit a true rival: planted-truth ablations of the shared
judge (M1.1's measurement model and the sec. 7 rules) on M1b part (i) and M3 phase 1 (2026-11-06).

A DIAGNOSTIC. It changes no prereg, no rule and no src/kinetic_core/ file, and ships nothing. Its only inputs are the
two cases' own draw ensembles (engine predictions already on disk), SYNTHETIC data generated from a planted truth, and
development / spent rows used only to pick a realistic truth. It never reads a case's decisive data on real values
(M1b's C-Mey / C-Hou / C-Zhou-ARP controls; M3's claim units) and never the blind hofmann2000b Figs. 2 - 3 (the M3
layer refuses the keyed table; this script never names it outside this sentence).

WHY. Every power check of the last two weeks found ~0 power to credit a rival even when the synthetic data were
generated FROM that rival (M3 phase 1: results/validation/kinetic_core_m3_phase1_power.md; M1b part (i): PR #94,
kinetic_core_m1b_part_i_power.md, 0 / 200 for every pair). A judge that cannot credit a planted truth is deaf. Five
shared components are suspected; each is ablated ONE AT A TIME, then all together, in a full 2^5 factorial:

  (i)   omega  : H0's within-laboratory 0.419 dex (estimated on thiol concentrations) in every unit's covariance (M1b:
                 as exact +-omega per-species re-solves, delta-pH 0.4 - 3.8 pH units; M3: per point) and as the
                 adequacy gate's |dmu| alternative. Ablation: the measurement's OWN reported noise only (sec. 7 item 5);
                 the |dmu| alternative becomes the unit's own measurement sd. Within one laboratory, so no tau.
                 For M1b this IS the ablation an independent review asked for (2026-11-06: the +-omega corner term is a
                 mis-transfer: anti-scaled corners break stoichiometry, a joint worst-case bound is used as a 1-sd
                 scale, and it flattens B1 / B13 ~4x structurally): "per-point measurement noise: the printed
                 replicate SD, else max(read error, meter 0.05 pH)" (the prereg's sec. 3.3 measurement column; the
                 censored C-Mey at 2 r^2), with delta-pH set to 0.
  (ii)  nuisance: coordinates integrated over their full prior bands. Ablation: conditioned on the case's development
                 / spent rows by importance weighting of the same draws (lppd(u | dev) = log E_{p(theta | dev)} L_u).
                 M1b: the three Zhou 2023 ARP + Cys rows (C-Zhou-Cys, B60 fit rows that calibrated the balance's two
                 constants), dropped from the guards when used. M3: the prereg already conditions on the calibrations
                 and prequentially on the fed units; the only spent rows outside it are G-T3 (hofmann2000 Table 3), added.
  (iii) conjunction: M1b: C-Mey (Cys) AND C-Hou must both favour, C-Zhou-ARP not oppose, guards not oppose after Holm;
                 M3: the primary favours AND no claim unit favours the other side (each at its own one-sided 0.05).
                 Ablation: Holm across units. M1b: credit iff Holm (family-wise 0.05, one-sided p = Phi(-d / sd_rep))
                 rejects for >= 1 decisive unit whose reading is the rival (adequacy kept), and no decisive or guard unit
                 favours the incumbent after Holm. M3: the primary as before; the claim-unit veto only after Holm.
  (iv)  forms  : required agreement over every form (M1b: omega 0.5x / 2x, r 0.02, the unforced C-Hou, B0-pinned; M3:
                 o_Ar 0.1, the Cu law, 20 mL headspace, and for transfer the (F) / (T) conditioning). Ablation: the
                 primary form only.
  (v)   unresolved: a density whose PSIS k-hat >= 0.7 or ESS < max(100, the tie-band precision rule) reads "tie".
                 Ablation: the estimate is read as estimated.

DESIGN (declared before any result; the numbers are the records' own unless stated)
- Truth = one rival's CENTRAL draw conditioned on development / spent rows: the draw nearest (standardised Euclidean
  distance, weighted sd) to the importance-weighted posterior mean of the rival's draws, weighted by the spent rows'
  likelihood at their OWN reported noise. M1b: C-Zhou-Cys (real, spent). M3: the units outside the pair's claim set
  (the prereg's own truth reweighting: route U1l; transfer and decay law U4, U1e, U1l; real, seen under the layer).
- Synthetic data, EVERY unit including the development rows, from the truth at each dataset's OWN reported noise
  (world "own"). A second world "omega" adds H0's omega to the generator as the original power checks did (data that
  carry an H0-sized model discrepancy): it is the misspecification check on false credit.
- Scoring: the cases' own functions (M1b: fit_m1b_ph_balance.py from PR #94's branch, imported read-only; M3:
  fit_m3_strecker_layer.py and src/m3_strecker_layer.py on main), the same draw ensembles (M1b: the screen's; M3: the
  v2.1 power check's 4 096-feasible sets, reloaded from its draw-set cache and checked against the record's counts).
- Tie bands: once per pair, truth, variant, form and unit, at the truth's expected data (M1b's declared reading,
  4 000 measurement-noise replicates; M3: 200, against the prereg's 40 per synthetic set; checked on a subset).
- R = 200 synthetic sets per truth (M1b); M3 per the --m3-reps argument (recorded). Power = P(credit rival a | truth
  a); false credit = P(credit a | truth b). Wilson 95 % intervals; power vs false credit by a one-sided Fisher test.
- Reproduction check (M1b): the baseline judge on the omega world at the power record's own truth draws (read by
  index from kinetic_core_m1b_part_i_power.json, never re-derived from decisive data) must give ~0 power, as recorded.

ATTRIBUTION AND CONTROL RULES (declared before any result)
- A component CARRIES the deafness if (a) ablating it alone raises power above the baseline (one-sided Fisher p < 0.05,
  200 vs 200), or (b) restoring it alone inside the all-together set lowers power below the all-together set (same
  test). Main effects over the 32 configurations are printed beside.
- False-credit control: the all-together set keeps false credit at or below nominal (0.05, the per-unit one-sided
  level and Holm's family-wise level) unless an exact one-sided binomial test rejects rate <= 0.05 at 0.05; both
  worlds are tested, and a failure is stated plainly.

RUN (one core, niced, threads pinned; each case in its own process because the M1b generator imports PR #94's src):
    OMP_NUM_THREADS=1 nice -n 10 ./scripts/docker_maillard.sh run "python scripts/diagnostics/judge_calibration.py \
        --case m1b --m1b-root <PR #94 worktree> --part-out <scratch>/m1b.json"
    ... --case m3 --m3-cache <draw-set cache> --part-out <scratch>/m3.json
    ... --case m1b --revised ... / --case m3 --revised ...   # the confirmation run of the REVISED rule (see that block)
    ... --merge <scratch>/m1b.json <scratch>/m3.json [<revised parts>]   # writes results/validation/judge_calibration.*

WHAT THIS IS NOT. Not a rule change. Any revision it suggests must be frozen by principle (sec. 7 item 5: "the
measurement's own reported noise") and confirmed independently BEFORE it is applied to any case, and the amended
judge's power must then be accepted whatever it is.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import dataclasses  # noqa: E402
import hashlib  # noqa: E402
import itertools  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import pickle  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
STEM = "judge_calibration"
GENERATED_BY = "scripts/diagnostics/judge_calibration.py"
COMPONENTS = ("omega", "nuisance", "conjunction", "forms", "unresolved")
COMPONENT_LABEL = {"omega": "(i) omega -> own per-point measurement noise", "nuisance": "(ii) nuisance conditioned on dev rows",
                   "conjunction": "(iii) Holm instead of the conjunction", "forms": "(iv) primary form only",
                   "unresolved": "(v) unresolved density read as estimated"}
Z90 = 1.6448536269514722
NOMINAL = 0.05
M1B_R = 200
M1B_TIE_R = 2000          # the prereg's 4 000 halved for compute (sd_rep relative SE ~1.6 %); declared
M1B_REPRO_R = int(os.environ.get("JUDGE_REPRO_R", "50"))
M1B_SEED = 20261107
M3_TIE_R = 200
M3_SEED = 20261108
M3_SCORE_FEASIBLE, M3_TRUTH_FEASIBLE, M3_SET_SEED = 4096, 2048, 20261105     # the v2.1 power check's sets


# ================================================================================================ statistics

def wilson(k: int, n: int, z: float = 1.959963984540054) -> List[float]:
    if n == 0:
        return [float("nan"), float("nan")]
    p = k / n
    den = 1.0 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [round(max(0.0, c - h), 4), round(min(1.0, c + h), 4)]


def fisher_greater(k1: int, k2: int, n1: int, n2: int) -> float:
    """One-sided p that rate 1 > rate 2."""
    from scipy.stats import fisher_exact
    return float(fisher_exact([[k1, n1 - k1], [k2, n2 - k2]], alternative="greater")[1])


def binom_above(k: int, n: int, p0: float = NOMINAL) -> float:
    """One-sided exact p of k or more successes when the rate is p0 (H0: rate <= p0)."""
    from scipy.stats import binom
    return float(binom.sf(k - 1, n, p0))


def holm_reject(p: Mapping[str, float], alpha: float = NOMINAL) -> Dict[str, bool]:
    items = sorted(p.items(), key=lambda kv: (kv[1], kv[0]))
    m, out, still = len(items), {}, True
    for i, (g, pv) in enumerate(items):
        still = still and pv <= alpha / (m - i)
        out[g] = bool(still)
    return out


def phi(x: float) -> float:
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def p_side(d: float, sd: float, side: str) -> float:
    """One-sided p that a unit favours ``side`` ('a': d > 0) by chance: Phi(-d / sd_rep) for a, Phi(d / sd_rep) for b
    (the M1b generator's guard convention). An infeasible rival (d = +-inf) gives p 0 for the other side."""
    if not math.isfinite(d):
        if math.isnan(d):
            return 1.0
        return 0.0 if (d > 0) == (side == "a") else 1.0
    if not (math.isfinite(sd) and sd > 0):
        return 0.5
    return phi(-d / sd) if side == "a" else phi(d / sd)


def lse(x: np.ndarray, axis: Optional[int] = None):
    x = np.asarray(x, float)
    m = np.max(x, axis=axis, keepdims=True)
    m = np.where(np.isfinite(m), m, 0.0)
    out = np.log(np.sum(np.exp(x - m), axis=axis, keepdims=True)) + m
    return float(out.squeeze()) if axis is None else np.squeeze(out, axis=axis)


def configs() -> List[Dict[str, bool]]:
    return [dict(zip(COMPONENTS, bits)) for bits in itertools.product((False, True), repeat=len(COMPONENTS))]


def config_name(c: Mapping[str, bool]) -> str:
    on = [k for k in COMPONENTS if c[k]]
    if not on:
        return "baseline"
    if len(on) == len(COMPONENTS):
        return "all"
    return "+".join(on)


def central_index(x: np.ndarray, log_w: np.ndarray) -> Tuple[int, Dict[str, float]]:
    """The draw nearest the weighted posterior mean (columns standardised by the weighted sd, floored at 1e-3 of the
    prior sd). Returns the index and the weights' ESS."""
    lw = np.where(np.isfinite(log_w), log_w, -np.inf)
    w = np.exp(lw - np.max(lw[np.isfinite(lw)]))
    w = w / w.sum()
    mu = w @ x
    sd_w = np.sqrt(np.maximum(w @ (x - mu) ** 2, 0.0))
    sd_p = np.std(x, axis=0) + 1e-12
    sd = np.maximum(sd_w, 1e-3 * sd_p)
    dist = np.sum(((x - mu) / sd) ** 2, axis=1)
    j = int(np.argmin(dist))
    ess = float(1.0 / np.sum(w ** 2))
    return j, {"ess_dev_weights": round(ess, 2), "n_draws": int(x.shape[0]), "weight_of_truth": float(w[j])}


def tally(rows: Sequence[Mapping[str, Any]], truth_a: str, truth_b: str, getter) -> Dict[str, Any]:
    ka = sum(bool(getter(r)) for r in rows if r["truth"] == truth_a)
    kb = sum(bool(getter(r)) for r in rows if r["truth"] == truth_b)
    na = sum(1 for r in rows if r["truth"] == truth_a)
    nb = sum(1 for r in rows if r["truth"] == truth_b)
    return {"power": ka / na if na else float("nan"), "power_k": ka, "n_a": na, "power_wilson95": wilson(ka, na),
            "false_credit": kb / nb if nb else float("nan"), "false_k": kb, "n_b": nb,
            "false_wilson95": wilson(kb, nb), "fisher_p_power_gt_false": fisher_greater(ka, kb, na, nb),
            "p_false_above_nominal": binom_above(kb, nb)}


def _jsonable(x):
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, float) and not math.isfinite(x):
        return None if math.isnan(x) else ("inf" if x > 0 else "-inf")
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    return x


def sha256_file(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def memory_free_pct() -> Optional[int]:
    try:
        out = subprocess.run(["memory_pressure"], capture_output=True, text=True, timeout=30).stdout
        for line in out.splitlines():
            if "free percentage" in line:
                return int(line.split(":")[1].strip().rstrip("%"))
    except Exception:
        return None
    return None


# ================================================================================================ M1b part (i)

def run_m1b(args) -> Dict[str, Any]:
    m1b_root = Path(args.m1b_root).resolve()
    for p in (str(m1b_root / "scripts" / "generators"), str(m1b_root)):
        if p not in sys.path:
            sys.path.insert(0, p)
    import fit_m1b_ph_balance as G  # the case's own generator, read-only (PR #94's branch)
    from src.kinetic_core import predictive_densities as pdens
    if Path(G.__file__).resolve().parents[2] != m1b_root or Path(pdens.__file__).resolve().parents[2] != m1b_root:
        raise SystemExit("the M1b generator or its src did not import from --m1b-root")
    draws_path = m1b_root / "results" / "validation" / "kinetic_core_m1b_part_i_screen_draws.json.gz"
    power_path = m1b_root / "results" / "validation" / "kinetic_core_m1b_part_i_power.json"
    ens = G.load_screen_ensembles(draws_path)
    t0 = time.time()

    # ---- noise templates: the units' reported noise only; every decisive y is dropped here and never used
    raw = G.measured_units()
    tmpl: Dict[str, Any] = {}
    for u, d in raw.items():
        tmpl[u] = (dataclasses.replace(d, within=(True,) * len(d.within)) if d.censored
                   else dataclasses.replace(d, y=tuple(float("nan") for _ in d.y)))
    zc_real = raw["C-Zhou-Cys"]                      # SPENT rows (B60 fit rows; the balance's two constants)
    del raw
    R_METER = G.METER_R

    def own_scale(u: str, r: float) -> float:
        d = tmpl[u]
        return math.sqrt(2.0) * r if d.censored else float(np.mean(d.sd))

    def forms_for(noise: str, b: str):
        fs = [f for f in G.LIKELIHOOD_FORMS if not (noise == "own" and f.name in ("omega_half", "omega_double"))]
        if b == "B0":
            fs.append(G.RuleForm("vs_b0_pinned", b_pinned=True))
        return fs

    def arrays(E, unit: str, form, noise: str) -> Tuple[np.ndarray, np.ndarray]:
        f = E.f[unit]
        dl = E.delta[unit][form.omega] if noise == "omega" else np.zeros_like(f)
        return f, dl

    def with_data(u: str, obs: Mapping[str, Any]):
        d = tmpl[u]
        if d.censored:
            return dataclasses.replace(d, within=tuple(map(bool, obs["within"])))
        return dataclasses.replace(d, y=tuple(map(float, obs["y"])))

    def dev_logw(E, zc_obs, form, noise: str) -> Optional[np.ndarray]:
        f, dl = arrays(E, "C-Zhou-Cys", form, noise)
        if f.shape[0] == 1:
            return None
        return G.unit_loglik(with_data("C-Zhou-Cys", zc_obs), f, dl, form.r_meter)

    def score(u: str, obs, f, dl, lw, r, prec, noise):
        dat = with_data(u, obs)
        ll = G.unit_loglik(dat, f, dl, r)
        if not np.all(np.isfinite(ll)):
            raise pdens.MeasurementModelRefusal(f"non-finite likelihood on {u}")
        if ll.size == 1:
            lppd, acc, khat = float(ll[0]), True, None
        elif lw is None:
            est = pdens.importance_lppd(ll, np.zeros_like(ll), precision_ess=prec)
            lppd, acc, khat = float(est.lppd), bool(est.accepted), est.k_hat
        else:
            en = pdens.importance_lppd(lw + ll, np.zeros_like(ll), precision_ess=prec)
            ed = pdens.importance_lppd(lw, np.zeros_like(ll))
            lppd, acc, khat = float(en.lppd - ed.lppd), bool(en.accepted and ed.accepted), en.k_hat
        post = ll + (0.0 if lw is None else lw)
        w = np.exp(post - np.max(post))
        w = w / w.sum()
        mean = w @ f
        scale = float(w @ np.mean(dl, axis=1)) if noise == "omega" else own_scale(u, r)
        if dat.censored:
            p = float(np.sum(w * np.exp(ll)))
        else:
            yy = np.asarray(dat.y, float)
            var = np.asarray(dat.sd, float) ** 2 + dl ** 2
            cov = np.diag(w @ var) + (f - mean).T @ ((f - mean) * w[:, None])
            p = pdens.posterior_predictive_p(yy, mean, cov)
        return lppd, acc, khat, mean, p, scale

    def band(u: str, exp_obs, fa, dla, lwa, fb, dlb, lwb, r, rng) -> float:
        rep = G.replicate_data(with_data(u, exp_obs), rng, M1B_TIE_R, r)
        dat = tmpl[u]
        la0 = 0.0 if lwa is None else lwa
        lb0 = 0.0 if lwb is None else lwb
        na = lse(np.asarray(la0) + np.zeros(fa.shape[0]))
        nb = lse(np.asarray(lb0) + np.zeros(fb.shape[0]))
        ds = []
        for i in range(0, M1B_TIE_R, 500):
            sl = slice(i, min(i + 500, M1B_TIE_R))
            kw = ({"within": rep["within"][sl], "edge": rep["edge"][sl]} if dat.censored else {"y": rep["y"][sl]})
            la = lse(la0 + G.unit_loglik(dat, fa, dla, r, **kw), axis=1) - na
            lb = lse(lb0 + G.unit_loglik(dat, fb, dlb, r, **kw), axis=1) - nb
            ds.append(la - lb)
        dd = np.concatenate(ds)
        return float(np.std(dd, ddof=1))

    def synthetic(T, s: int, rng, units, world: str) -> Dict[str, Dict[str, np.ndarray]]:
        out = {}
        for u in units:
            f = T.f[u][s]
            dl = T.delta[u]["primary"][s] if world == "omega" else np.zeros_like(f)
            if tmpl[u].censored:
                d_obs = f + rng.normal(size=f.size) * np.sqrt(2.0 * R_METER ** 2 + dl ** 2)
                out[u] = {"within": np.abs(d_obs) <= G.MEYNIER_BOUND}
            else:
                out[u] = {"y": f + rng.normal(size=f.size) * np.sqrt(np.asarray(tmpl[u].sd) ** 2 + dl ** 2)}
        out["C-Hou-unforced"] = out["C-Hou"]
        return out

    def expected(T, s: int, units) -> Dict[str, Dict[str, np.ndarray]]:
        out = {}
        for u in units:
            f = T.f[u][s]
            out[u] = {"within": np.abs(f) <= G.MEYNIER_BOUND} if tmpl[u].censored else {"y": f.copy()}
        out["C-Hou-unforced"] = out["C-Hou"]
        return out

    # ---- truths: the central development-conditioned draw of each rival (own noise on the spent rows)
    truths: Dict[str, Dict[str, Any]] = {}
    for rv in G.RIVAL_NAMES:
        E = ens[rv]
        f = E.f["C-Zhou-Cys"]
        lw = G.unit_loglik(zc_real, f, np.zeros_like(f), R_METER)
        cols = []
        for c in G.coordinates(rv):
            v = np.array([d[c.name] for d in E.draws], float)
            cols.append(np.log10(v) if c.scale == "log10" else v)
        j, info = central_index(np.column_stack(cols), lw)
        truths[rv] = {"draw": j, "coordinates": E.draws[j], **info,
                      "zhou_cys_rows_pred_ph": [round(x, 3) for x in f[j].tolist()],
                      "zhou_cys_rows_meas_ph": list(zc_real.y)}
    recorded = json.loads(power_path.read_text())["pairs"]

    variants = [("omega", "bands"), ("own", "bands"), ("omega", "dev"), ("own", "dev")]
    BAND_CACHE: Dict[Tuple[str, str, str, int], Dict[Tuple[str, str, str, str], float]] = {}

    def pair_units(a: str, nuis: str) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
        guards = tuple(g for g in G.guards_of(a) if not (nuis == "dev" and g == "C-Zhou-Cys"))
        return G.DECISIVE, guards

    def evaluate(a: str, b: str, truth: str, s: int, world: str, r_sets: int, rng, primary_only_variants=False,
                 use_variants=None):
        """Per synthetic set: every variant's per-form per-unit readings."""
        vlist = variants if use_variants is None else use_variants
        T = ens[truth]
        all_units = tuple(dict.fromkeys(G.DECISIVE + G.guards_of(a) + ("C-Zhou-Cys",)))
        exp_obs = expected(T, s, all_units)
        # bands once per (variant, form, unit), at the truth's expected data (shared by both worlds)
        bands: Dict[Tuple[str, str, str, str], float] = BAND_CACHE.setdefault((a, b, truth, s), {})
        for noise, nuis in vlist:
            for form in forms_for(noise, b):
                if primary_only_variants and form.name != "primary":
                    continue
                Eb = ens["B0_pinned"] if form.b_pinned else ens[b]
                lwa = dev_logw(ens[a], exp_obs["C-Zhou-Cys"], form, noise) if nuis == "dev" else None
                lwb = dev_logw(Eb, exp_obs["C-Zhou-Cys"], form, noise) if nuis == "dev" else None
                dec, guards = pair_units(a, nuis)
                for u in dec + guards:
                    if (noise, nuis, form.name, u) in bands:
                        continue
                    unit = form.hou_unit if u == "C-Hou" else u
                    fa, dla = arrays(ens[a], unit, form, noise)
                    fb, dlb = arrays(Eb, unit, form, noise)
                    bands[(noise, nuis, form.name, u)] = band(unit, exp_obs[unit], fa, dla, lwa, fb, dlb, lwb,
                                                              form.r_meter, rng)
        rows = []
        for _ in range(r_sets):
            syn = synthetic(T, s, rng, all_units, world)
            rec: Dict[str, Any] = {"truth": truth}
            for noise, nuis in vlist:
                vk = f"{noise}|{nuis}"
                rec[vk] = {}
                for form in forms_for(noise, b):
                    if primary_only_variants and form.name != "primary":
                        continue
                    Eb = ens["B0_pinned"] if form.b_pinned else ens[b]
                    lwa = dev_logw(ens[a], syn["C-Zhou-Cys"], form, noise) if nuis == "dev" else None
                    lwb = dev_logw(Eb, syn["C-Zhou-Cys"], form, noise) if nuis == "dev" else None
                    dec, guards = pair_units(a, nuis)
                    units = {}
                    for u in dec + guards:
                        unit = form.hou_unit if u == "C-Hou" else u
                        fa, dla = arrays(ens[a], unit, form, noise)
                        fb, dlb = arrays(Eb, unit, form, noise)
                        sd_rep = bands[(noise, nuis, form.name, u)]
                        prec = pdens.tie_band_precision_ess(sd_rep)
                        la, acc_a, ka, ma, pa, sa = score(unit, syn[unit], fa, dla, lwa, form.r_meter, prec, noise)
                        lb, acc_b, kb, mb, pb, sb = score(unit, syn[unit], fb, dlb, lwb, form.r_meter, prec, noise)
                        d = la - lb
                        reading = pdens.unit_reading(d, Z90 * sd_rep, pa, pb, ma, mb, max(sa, sb), G.ADEQUACY_P)
                        resolved = acc_a and acc_b
                        units[u] = {"d": d, "sd": sd_rep, "raw": reading, "rule": reading if resolved else "tie",
                                    "resolved": resolved}
                    rec[vk][form.name] = units
            rows.append(rec)
        return rows

    def favoured_one(units: Mapping[str, Mapping[str, Any]], a: str, nuis: str, holm: bool, key: str) -> bool:
        dec, guards = pair_units(a, nuis)
        if not holm:
            gp = {g: p_side(units[g]["d"], units[g]["sd"], "b") for g in guards}
            rej = holm_reject(gp)
            veto = any(rej[g] and units[g][key] == "b" for g in guards)
            return (units["C-Mey-Cys"][key] == "a" and units["C-Hou"][key] == "a"
                    and units["C-Zhou-ARP"][key] != "b" and not veto)
        rej_a = holm_reject({u: p_side(units[u]["d"], units[u]["sd"], "a") for u in dec})
        credit = any(rej_a[u] and units[u][key] == "a" for u in dec)
        rej_b = holm_reject({u: p_side(units[u]["d"], units[u]["sd"], "b") for u in dec + guards})
        veto = any(rej_b[u] and units[u][key] == "b" for u in dec + guards)
        return credit and not veto

    def pooled_one(units, a, nuis, key) -> bool:
        dec, guards = pair_units(a, nuis)
        dsum = sum(units[u]["d"] for u in dec)
        sdsum = math.sqrt(sum(units[u]["sd"] ** 2 for u in dec))
        ok = key == "raw" or all(units[u]["resolved"] for u in dec)
        rej_b = holm_reject({u: p_side(units[u]["d"], units[u]["sd"], "b") for u in dec + guards})
        veto = any(rej_b[u] and units[u][key] == "b" for u in dec + guards)
        return bool(ok and dsum > Z90 * sdsum and not veto)

    def credited(rec, a: str, cfg: Mapping[str, bool], pooled: bool = False) -> bool:
        noise = "own" if cfg["omega"] else "omega"
        nuis = "dev" if cfg["nuisance"] else "bands"
        key = "raw" if cfg["unresolved"] else "rule"
        per = rec[f"{noise}|{nuis}"]
        names = ["primary"] if cfg["forms"] else list(per)
        if pooled:
            return all(pooled_one(per[n], a, nuis, key) for n in names)
        return all(favoured_one(per[n], a, nuis, cfg["conjunction"], key) for n in names)

    def diagnostics(rows, a: str, truth: str, vk: str) -> Dict[str, Any]:
        out = {}
        sel = [r for r in rows if r["truth"] == truth]
        for u in sel[0][vk]["primary"]:
            ds = np.array([r[vk]["primary"][u]["d"] for r in sel], float)
            out[u] = {"d_median": float(np.median(ds)), "band": float(Z90 * sel[0][vk]["primary"][u]["sd"]),
                      "resolved_rate": float(np.mean([r[vk]["primary"][u]["resolved"] for r in sel])),
                      "reads_a_raw": float(np.mean([r[vk]["primary"][u]["raw"] == "a" for r in sel])),
                      "reads_b_raw": float(np.mean([r[vk]["primary"][u]["raw"] == "b" for r in sel])),
                      "both_inadequate_raw": float(np.mean([r[vk]["primary"][u]["raw"] == "both inadequate"
                                                            for r in sel]))}
        return out

    rng = np.random.default_rng(M1B_SEED)
    out_pairs: Dict[str, Any] = {}
    want = [x for x in args.m1b_pairs.split(",") if x]
    for a, b in G.PAIRS:
        key = f"{a}_vs_{b}"
        if want and key not in want:
            continue
        t1 = time.time()
        res: Dict[str, Any] = {"a": a, "b": b, "truth_a": truths[a], "truth_b": truths[b], "worlds": {}}
        for world in ("own", "omega"):
            rows = (evaluate(a, b, a, truths[a]["draw"], world, args.m1b_reps, rng)
                    + evaluate(a, b, b, truths[b]["draw"], world, args.m1b_reps, rng))
            tab = {}
            for cfg in configs():
                tab[config_name(cfg)] = {**tally(rows, a, b, lambda r, c=cfg: credited(r, a, c)), "config": cfg}
            allcfg = {k: True for k in COMPONENTS}
            tab["all|pooled_lppd"] = {**tally(rows, a, b, lambda r: credited(r, a, allcfg, pooled=True)),
                                      "config": allcfg, "note": "all, with the pooled decisive lppd in place of Holm"}
            res["worlds"][world] = {"table": tab,
                                    "unit_diagnostics_baseline": {t: diagnostics(rows, a, t, "omega|bands")
                                                                  for t in (a, b)},
                                    "unit_diagnostics_all": {t: diagnostics(rows, a, t, "own|dev") for t in (a, b)}}
            print(f"m1b {key} world {world}: baseline power {tab['baseline']['power']:.3f} false "
                  f"{tab['baseline']['false_credit']:.3f} | all power {tab['all']['power']:.3f} false "
                  f"{tab['all']['false_credit']:.3f} ({time.time() - t1:.0f} s)", flush=True)
        # the reproduction: omega world, the power record's own truth draws, the baseline judge
        rv = [("omega", "bands"), ("own", "dev")]
        rows = (evaluate(a, b, a, int(recorded[key]["power"]["truth_draw"]), "omega", M1B_REPRO_R, rng, use_variants=rv)
                + evaluate(a, b, b, int(recorded[key]["false_credit"]["truth_draw"]), "omega", M1B_REPRO_R, rng,
                           use_variants=rv))
        base = {k: False for k in COMPONENTS}
        raw_only = {**base, "unresolved": True}
        res["reproduction"] = {
            "truth_draws": [int(recorded[key]["power"]["truth_draw"]), int(recorded[key]["false_credit"]["truth_draw"])],
            "recorded_power_rate": recorded[key]["power"]["rate"],
            "recorded_power_rate_unresolved_as_tie": recorded[key]["power"]["rate_unresolved_density_as_tie"],
            "baseline": tally(rows, a, b, lambda r: credited(r, a, base)),
            "baseline_unresolved_read_as_estimated": tally(rows, a, b, lambda r: credited(r, a, raw_only)),
            "all": tally(rows, a, b, lambda r: credited(r, a, {k: True for k in COMPONENTS}))}
        res["seconds"] = round(time.time() - t1, 1)
        out_pairs[key] = res
    return {"case": "M1b part (i)", "pairs": out_pairs, "truths": truths, "R": args.m1b_reps, "tie_R": M1B_TIE_R,
            "repro_R": M1B_REPRO_R, "seconds": round(time.time() - t0, 1),
            "inputs": {"generator": str(Path(G.__file__).resolve()), "generator_sha256": sha256_file(Path(G.__file__)),
                       "draws": str(draws_path), "draws_sha256": sha256_file(draws_path),
                       "power_record_sha256": sha256_file(power_path),
                       "predictive_densities_sha256": sha256_file(Path(pdens.__file__))},
            "variants": [f"{n}|{c}" for n, c in variants],
            "forms": {"omega": [f.name for f in forms_for("omega", "B0")], "own": [f.name for f in forms_for("own", "B0")]},
            "spent_rows_for_conditioning": "C-Zhou-Cys (Zhou 2023 Fig. 2 ARP + Cys, the three B60 fit rows)"}


# ================================================================================================ M3 phase 1

def run_m3(args) -> Dict[str, Any]:
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    gen_dir = str(ROOT / "scripts" / "generators")
    if gen_dir not in sys.path:
        sys.path.insert(0, gen_dir)
    import fit_m3_strecker_layer as G
    from src import m3_strecker_layer as L
    from src.kinetic_core import predictive_densities as pd
    from dataclasses import replace

    data = L.load_data()
    real = G.real_obs(data)                       # seen units only (U3 is not in it: synthetic only)
    t0 = time.time()

    # ---- draw sets: the v2.1 power check's, from its cache, matched by (rival, form, seed) and checked against the
    #      record's feasible / drawn counts
    rec_sets = json.loads((ROOT / "results" / "validation" / "kinetic_core_m3_phase1_power.json").read_text())["draw_sets"]
    spec_args = argparse.Namespace(score_feasible=M3_SCORE_FEASIBLE, truth_feasible=M3_TRUTH_FEASIBLE, seed=M3_SET_SEED)
    specs = G.power_set_specs(spec_args, tuple(G.POWER_RIVALS))

    class _U(pickle.Unpickler):
        def find_class(self, module, name):
            if name == "DrawSet" and module in ("__main__", "fit_m3_strecker_layer"):
                return G.DrawSet
            return super().find_class(module, name)

    loaded = []
    for p in sorted(Path(args.m3_cache).glob("drawset_*.pkl")):
        with p.open("rb") as fh:
            ds = _U(fh).load()
        loaded.append((p, ds))
    sets: Dict[Tuple[str, str], Any] = {}
    set_log = {}
    for sp in specs:
        hit = [(p, ds) for p, ds in loaded if ds.rival == sp["rival"] and ds.seed == sp["seed"] and ds.form == sp["form"]]
        if len(hit) != 1:
            raise SystemExit(f"draw set {sp['key']} not found exactly once in {args.m3_cache} ({len(hit)})")
        p, ds = hit[0]
        name, fname = sp["key"]
        want = rec_sets[f"{name}|{fname}"]
        if (ds.n, ds.n_drawn, ds.n_failed) != (want["feasible"], want["drawn"], want["failed"]):
            raise SystemExit(f"draw set {sp['key']}: {ds.n}/{ds.n_drawn} differs from the record {want}")
        sets[(name, fname)] = ds
        set_log[f"{name}|{fname}"] = {"feasible": ds.n, "drawn": ds.n_drawn, "file": p.name,
                                      "sha256": sha256_file(p)[:16]}
    print(f"m3: {len(sets)} draw sets loaded and matched to the record ({time.time() - t0:.0f} s)", flush=True)

    FORMS = G.POWER_FORMS
    FED = ("U4", "U1e", "U1l")
    PAIRS = {"route": ("R1", "R2"), "transfer": ("R1", "R3"), "decay_law": ("R4", "R1drift")}
    REWEIGHT = {"route": ("U1l",), "transfer": ("U4", "U1e", "U1l"), "decay_law": ("U4", "U1e", "U1l")}
    CLAIMS = {c: G.CONTRASTS[c]["claim"] for c in PAIRS}

    def form_v(fname: str, noise: str):
        f = FORMS[fname]
        return replace(f, omega_scale=0.0) if noise == "own" else f

    def own_scale(u) -> float:
        return float(math.sqrt(np.median(np.diag(L.unit_covariance(u, L.PRIMARY, measurement_only=True)))))

    # ---- synthetic data
    def synth_unit(unit, f_points, rng, world):
        y = unit.design @ f_points
        cov = L.unit_covariance(unit, L.PRIMARY, measurement_only=(world == "own"))
        return replace(unit, y=y + rng.multivariate_normal(np.zeros(len(y)), cov))

    def trapped(tset, j, rng, world, noiseless=False):
        f3 = tset.points["U3"][j]
        if noiseless:
            y_log = f3.copy()
        else:
            cov_p = L.unit_covariance(L.u3_unit(10.0 ** f3), L.PRIMARY, measurement_only=(world == "own"))
            y_log = f3 + rng.multivariate_normal(np.zeros(len(f3)), cov_p)
        u3 = L.u3_unit(10.0 ** y_log)
        y_lin = np.where(u3.censored, 0.5 * L.U3_CENSOR, 10.0 ** u3.y)
        err = data.arp_t2_err["argon"][0]
        sd = math.sqrt(err ** 2 + L.S_REP_DEX ** 2 + L.SIGMA_RUN_DEX ** 2
                       + (L.OMEGA_DEX ** 2 if world == "omega" else 0.0))
        pa = 10.0 ** (tset.points["U4"][j][0] + (0.0 if noiseless else sd * rng.standard_normal()))
        return {"U3": u3, "U3a": L.u3a_unit(y_lin, pa, err), "U3b": L.u3b_unit(y_lin)}

    def uarp_noiseless(k_1min):
        t = np.arange(0.0, 501.0, 1.0)
        ts, vs = data.forced_inputs["PHEARP"]
        a_ = np.interp(t, ts, vs)
        f = np.maximum(np.gradient(a_, t) + k_1min * a_, 0.0)
        arp = np.zeros_like(a_)
        for i in range(1, t.size):
            k = k_1min[i - 1]
            ex = math.exp(-k)
            arp[i] = arp[i - 1] * ex + (f[i - 1] / k * (1 - ex) if k > 0 else f[i - 1])
        obs = {tt: max(round(arp[int(tt)] * 10.0), 1) / 10.0 for tt in sorted(data.arp_glc)}
        return np.array([[math.log(obs[a0] / obs[b0]) for a0, b0 in L.ARP_INTERVALS]])

    UNITS_NEEDED = ("U4", "U1e", "U1l", "GT3", "U5a", "U5b")

    def make_obs(tset, j, rng, world, noiseless=False):
        units = dict(real.units)
        for u in UNITS_NEEDED:
            if noiseless:
                units[u] = replace(real.units[u], y=real.units[u].design @ tset.points[u][j])
            else:
                units[u] = synth_unit(real.units[u], tset.points[u][j], rng, world)
        units.update(trapped(tset, j, rng, world, noiseless))
        uy = uarp_noiseless(tset.k_glc[j]) if noiseless else G.synth_uarp(data, tset.k_glc[j], rng)
        return G.Obs(units, uy)

    # ---- scoring with a per-observation loglik cache
    class Scorer:
        def __init__(self, obs):
            self.obs = obs
            self.cache: Dict[Tuple[int, str, Any], np.ndarray] = {}

        def ll(self, ds, unit, fv, obs=None):
            if obs is not None:
                return G.loglik(ds, unit, obs, fv, data)
            k = (id(ds), unit, fv)
            if k not in self.cache:
                self.cache[k] = G.loglik(ds, unit, self.obs, fv, data)
            return self.cache[k]

        def cond(self, ds, units, fv):
            out = np.zeros(ds.n)
            for u in units:
                out = out + self.ll(ds, u, fv)
            return out

    def sd_rep_item(sc: Scorer, a, b, item, fv, rng) -> float:
        kind, units, cond = item
        ca, cb = sc.cond(a, cond, fv), sc.cond(b, cond, fv)
        nca, ncb = G.lse(ca), G.lse(cb)
        reps = []
        for _ in range(M3_TIE_R):
            o2 = sc.obs
            for u in units:
                o2 = G.redraw(o2, u, data, rng, fv)
            la = sum(G.loglik(a, u, o2, fv, data) for u in units)
            lb = sum(G.loglik(b, u, o2, fv, data) for u in units)
            reps.append(G.lse(ca + la) - nca - G.lse(cb + lb) + ncb)
        reps = np.array([x for x in reps if np.isfinite(x)])
        return float(np.std(reps, ddof=1)) if reps.size > 2 else float("nan")

    def verdict(sc: Scorer, a, b, item, fv, noise, sd_rep) -> Dict[str, Any]:
        kind, units, cond = item
        if a.n == 0 or b.n == 0:
            v = G._one_sided(a.n > 0, b.n > 0)
            return {"d": v["d"], "sd": float("nan"), "raw": v["reading"], "rule": v["reading"], "resolved": True}
        ca, cb = sc.cond(a, cond, fv), sc.cond(b, cond, fv)
        la = sum(sc.ll(a, u, fv) for u in units)
        lb = sum(sc.ll(b, u, fv) for u in units)
        ra, rb = G._representable(a, ca, la), G._representable(b, cb, lb)
        if not (ra and rb):
            v = G._one_sided(ra, rb)
            return {"d": v["d"], "sd": float("nan"), "raw": v["reading"], "rule": v["reading"], "resolved": True}
        ea, eb = G.conditional_lppd(ca, la, True), G.conditional_lppd(cb, lb, True)
        d = ea["lppd"] - eb["lppd"]
        band = Z90 * sd_rep
        resolved = bool(ea.get("resolved", True) and eb.get("resolved", True))
        if kind == "block" or units[0] == "UARP":
            raw = "tie" if not (np.isfinite(d) and abs(d) > band) else ("a" if d > 0 else "b")
        else:
            u = sc.obs.units[units[0]]
            mu_a, cov_a = G.predictive_moments(a, u, a.points[units[0]], ca, fv)
            mu_b, cov_b = G.predictive_moments(b, u, b.points[units[0]], cb, fv)
            p_a, p_b = G.adequacy_p(u, mu_a, cov_a), G.adequacy_p(u, mu_b, cov_b)
            scale = (L.OMEGA_DEX * fv.omega_scale * math.sqrt(u.c)) if noise == "omega" else own_scale(u)
            raw = pd.unit_reading(d, band, p_a, p_b, mu_a, mu_b, scale)
        return {"d": d, "sd": sd_rep, "raw": raw, "rule": raw if resolved else "tie", "resolved": resolved}

    def items_for(cname: str, dev: bool) -> Dict[str, Any]:
        """The contrast's scored items (kind, units, cond) by role, as fit_m3_strecker_layer._contrast_on builds them;
        ``dev`` prepends G-T3 to every conditioning set."""
        pre = ("GT3",) if dev else ()
        if cname == "route":
            cond = {"U4": (), "U1e": ("U4",), "U3": FED, "U3a": FED, "U3b": FED}
            per = {u: ("unit", (u,), pre + cond[u]) for u in CLAIMS["route"]}
            return {"F": {"primary": ("block", ("U3a", "U3b"), pre + FED), "per": per}}
        if cname == "decay_law":
            block = ("U3", "U5a", "U5b", "UARP")
            per, prev = {}, FED
            for u in block:
                per[u] = ("unit", (u,), pre + prev)
                prev = prev + (u,)
            return {"F": {"primary": ("block", block, pre + FED), "per": per}}
        out = {}
        for tag, cond0 in (("F", FED + ("U3",)), ("T", FED)):
            per, prev = {}, cond0
            for u in G.GLC_ORDER:
                per[u] = ("unit", (u,), pre + prev)
                prev = prev + (u,)
            out[tag] = {"primary": ("block", G.GLC_ORDER, pre + cond0), "per": per}
        return out

    def credit_side(prim, per, key, holm: bool) -> str:
        for side, other in (("a", "b"), ("b", "a")):
            if prim[key] != side:
                continue
            if not holm:
                veto = any(v[key] == other for v in per.values())
            else:
                rej = holm_reject({u: p_side(v["d"], v["sd"], other) for u, v in per.items()})
                veto = any(rej[u] and v[key] == other for u, v in per.items())
            return side if not veto else "none"
        return "none"

    def credit(rec, cname, cfg, direction: str) -> bool:
        noise = "own" if cfg["omega"] else "omega"
        vk = f"{noise}|{'dev' if cfg['nuisance'] else 'prereg'}"
        key = "raw" if cfg["unresolved"] else "rule"
        fnames = ["primary"] if cfg["forms"] else list(FORMS)
        sides = []
        for fname in fnames:
            r = rec[vk][fname]
            tags = ["F"] if (cfg["forms"] or cname != "transfer") else ["F", "T"]
            per_tag = [credit_side(r[t]["primary"], r[t]["per"], key, cfg["conjunction"]) for t in tags]
            sides.append(per_tag[0] if len(set(per_tag)) == 1 else "none")
        return all(s == direction for s in sides)

    variants = [("omega", False), ("own", False), ("omega", True), ("own", True)]

    def run_pair(cname: str, world: str, reps: int, forms_subset: Sequence[str], rng) -> Dict[str, Any]:
        a_name, b_name = PAIRS[cname]
        rows, truth_info = [], {}
        for truth in (a_name, b_name):
            tset = sets[(truth, "truth")]
            lw = np.zeros(tset.n)
            for u in REWEIGHT[cname]:
                lw = lw + G.loglik(tset, u, real, replace(L.PRIMARY, omega_scale=0.0), data)
            j, info = central_index(tset.u, lw)
            truth_info[truth] = {"draw": j, **info,
                                 "theta": {k: float(np.asarray(v)[j]) for k, v in tset.theta.items()
                                           if np.ndim(v) == 1 and len(v) == tset.n}}
            # bands at the truth's expected data
            exp_sc = Scorer(make_obs(tset, j, rng, world, noiseless=True))
            bands = {}
            for noise, dev in variants:
                for fname in forms_subset:
                    fv = form_v(fname, noise)
                    a, b = sets[(a_name, fname)], sets[(b_name, fname)]
                    for tag, its in items_for(cname, dev).items():
                        for role, item in [("primary", its["primary"])] + list(its["per"].items()):
                            k = (noise, dev, fname, item)
                            if k not in bands and a.n > 0 and b.n > 0:
                                bands[k] = sd_rep_item(exp_sc, a, b, item, fv, rng)
            print(f"m3 {cname} {world} truth {truth}: {len(bands)} bands ({time.time() - t0:.0f} s)", flush=True)
            for i in range(reps):
                sc = Scorer(make_obs(tset, j, rng, world))
                rec: Dict[str, Any] = {"truth": truth}
                for noise, dev in variants:
                    vk = f"{noise}|{'dev' if dev else 'prereg'}"
                    rec[vk] = {}
                    for fname in forms_subset:
                        fv = form_v(fname, noise)
                        a, b = sets[(a_name, fname)], sets[(b_name, fname)]
                        rec[vk][fname] = {}
                        for tag, its in items_for(cname, dev).items():
                            prim = verdict(sc, a, b, its["primary"], fv, noise,
                                           bands.get((noise, dev, fname, its["primary"]), float("nan")))
                            per = {u: verdict(sc, a, b, it, fv, noise, bands.get((noise, dev, fname, it), float("nan")))
                                   for u, it in its["per"].items()}
                            rec[vk][fname][tag] = {"primary": prim, "per": per}
                rows.append(rec)
                if (i + 1) % 50 == 0:
                    print(f"m3 {cname} {world} truth {truth}: {i + 1} / {reps} ({time.time() - t0:.0f} s)", flush=True)
        primary_only = list(forms_subset) == ["primary"]
        tab = {}
        for cfg in configs():
            if primary_only and not cfg["forms"]:
                continue
            for direction, hit, other in (("a", a_name, b_name), ("b", b_name, a_name)):
                tab[f"{config_name(cfg)}|credit_{hit}"] = {
                    **tally(rows, hit, other, lambda r, c=cfg, dd=direction: credit(r, cname, c, dd)), "config": cfg}
        diag = {}
        for truth in (a_name, b_name):
            sel = [r for r in rows if r["truth"] == truth]
            for vk in ("omega|prereg", "own|dev"):
                prim = [r[vk]["primary"]["F"]["primary"] for r in sel]
                diag[f"{truth}|{vk}"] = {
                    "primary_resolved_rate": float(np.mean([p["resolved"] for p in prim])),
                    "primary_d_median": float(np.median([p["d"] for p in prim])),
                    "primary_band": float(Z90 * prim[0]["sd"]) if math.isfinite(prim[0]["sd"]) else None,
                    "primary_reads_a_raw": float(np.mean([p["raw"] == "a" for p in prim])),
                    "primary_reads_b_raw": float(np.mean([p["raw"] == "b" for p in prim])),
                    "claim_unit_vetoes_raw": {u: float(np.mean([r[vk]["primary"]["F"]["per"][u]["raw"] == (
                        "b" if truth == a_name else "a") for r in sel])) for u in sel[0][vk]["primary"]["F"]["per"]}}
        return {"pair": [a_name, b_name], "truths": truth_info, "R": reps, "forms": list(forms_subset), "table": tab,
                "diagnostics_primary_form": diag}

    def band_check(cname: str, n_sets: int, rng) -> Dict[str, Any]:
        """The fixed-band approximation against the prereg's per-set bands (40 redraws), baseline variant, primary
        form, own world, truth a: agreement of the primary's and the claim units' raw readings."""
        a_name, b_name = PAIRS[cname]
        tset = sets[(a_name, "truth")]
        lw = sum(G.loglik(tset, u, real, replace(L.PRIMARY, omega_scale=0.0), data) for u in REWEIGHT[cname])
        j, _ = central_index(tset.u, lw)
        fv = form_v("primary", "omega")
        a, b = sets[(a_name, "primary")], sets[(b_name, "primary")]
        its = items_for(cname, False)["F"]
        exp_sc = Scorer(make_obs(tset, j, rng, "own", noiseless=True))
        fixed = {it: sd_rep_item(exp_sc, a, b, it, fv, rng) for it in [its["primary"]] + list(its["per"].values())}
        agree = total = 0
        ratio = []
        for _ in range(n_sets):
            sc = Scorer(make_obs(tset, j, rng, "own"))
            for it in [its["primary"]] + list(its["per"].values()):
                v_fixed = verdict(sc, a, b, it, fv, "omega", fixed[it])
                sd_set = _per_set_sd(sc, a, b, it, fv, rng, 40)
                v_set = verdict(sc, a, b, it, fv, "omega", sd_set)
                agree += int(v_fixed["raw"] == v_set["raw"])
                total += 1
                if math.isfinite(sd_set) and fixed[it] > 0:
                    ratio.append(sd_set / fixed[it])
        return {"contrast": cname, "n_sets": n_sets, "readings_compared": total, "agreement": agree / total,
                "sd_ratio_per_set_over_fixed_quantiles": [float(x) for x in np.quantile(ratio, [0.1, 0.5, 0.9])]}

    def _per_set_sd(sc, a, b, item, fv, rng, n):
        kind, units, cond = item
        ca, cb = sc.cond(a, cond, fv), sc.cond(b, cond, fv)
        nca, ncb = G.lse(ca), G.lse(cb)
        reps = []
        for _ in range(n):
            o2 = sc.obs
            for u in units:
                o2 = G.redraw(o2, u, data, rng, fv)
            reps.append(G.lse(ca + sum(G.loglik(a, u, o2, fv, data) for u in units)) - nca
                        - G.lse(cb + sum(G.loglik(b, u, o2, fv, data) for u in units)) + ncb)
        reps = np.array([x for x in reps if np.isfinite(x)])
        return float(np.std(reps, ddof=1)) if reps.size > 2 else float("nan")

    rng = np.random.default_rng(M3_SEED)
    out: Dict[str, Any] = {"case": "M3 phase 1", "pairs": {}, "draw_sets": set_log, "tie_R": M3_TIE_R}
    contrasts = [c for c in args.m3_contrasts.split(",") if c]
    for cname in contrasts:
        out["pairs"][cname] = {"worlds": {
            "own": run_pair(cname, "own", args.m3_reps, list(FORMS), rng),
            "omega": run_pair(cname, "omega", args.m3_reps_omega, ["primary"], rng)}}
        for w, v in out["pairs"][cname]["worlds"].items():
            k0 = [k for k in v["table"] if k.startswith("all|")][0]
            print(f"m3 {cname} world {w}: {k0} power {v['table'][k0]['power']:.3f} false "
                  f"{v['table'][k0]['false_credit']:.3f}", flush=True)
    if args.m3_band_check:
        out["band_check"] = [band_check(c, args.m3_band_check, rng) for c in contrasts]
    out["seconds"] = round(time.time() - t0, 1)
    out["R_own"], out["R_omega"] = args.m3_reps, args.m3_reps_omega
    out["inputs"] = {"layer_sha256": sha256_file(ROOT / "src" / "m3_strecker_layer.py"),
                     "generator_sha256": sha256_file(ROOT / "scripts" / "generators" / "fit_m3_strecker_layer.py"),
                     "draw_set_cache": "the v2.1 power check's cache (4096 / 2048 feasible, seed 20261105), counts "
                                       "matched to kinetic_core_m3_phase1_power.json"}
    out["spent_rows_for_conditioning"] = ("the prereg's own: calibrations (exact inversions) and the fed units "
                                          "prequentially; ablation (ii) adds G-T3 (hofmann2000 Table 3)")
    return out


# ================================================================================================ the REVISED rule
#
# The confirmation run (coordinator, 2026-11-06, after an independent review of this record's proposal: "adopt with
# fixes"). It scores the revised sec. 7 item 5 and M1.1 EXACTLY as worded in REVISED_RULE_TEXT, on synthetic data in
# three generator worlds, and fixes this record's own flaws: truths come from the null rival's spent-conditioned
# posterior OUTSIDE the scoring ensemble; thresholds are computed per synthetic set; an unresolved estimate is never
# read as resolved (the set is "not computable", with the cost of re-estimating it).
#
# DECLARED READINGS (where the wording needs a reading to be computed; each is printed in the record)
# - The null distribution of T for one synthetic set is the parametric bootstrap from the null rival's posterior given
#   THAT set's spent rows. It is computed by importance-reweighting one bootstrap pool per (pair, null rival): pool
#   element j has a truth theta_j drawn uniformly from the null rival's held-out draws and its own redrawn data (decisive,
#   spent, replicate SDs, noise components, discrete nuisance index), with T_j computed exactly as on the set; a set's
#   null is the pool weighted by p(y_spent(set) | theta_j). This is exact up to Monte Carlo error, which is reported per
#   set (weighted-quantile SE from the weights' ESS).
# - "Outside the ensemble": M1b's screen ensembles are split by draw parity (even: truths; odd: scoring); M3's truths
#   come from its separate truth sets (other Sobol seeds). Forced inputs are NOT redrawn (each redraw is a new engine
#   run per draw: not computable here; stated per case).
# - Required agreement: M1b, the forced / unforced C-Hou (a modelling choice); M3, o_Ar 0.1, the Cu law and the 20 mL
#   headspace (the prereg's "unsourced choices"), each form's T against its own null, data from the primary form.
#   Discrete nuisance: M1b's meter repeatability r in {0.02, 0.05, 0.10} pH (A1's set), equal prior weights, mixed in the
#   density. Omega forms dropped. B0-pinned (a separate pairing) is not run.
# - M1b noise (revised item 5): C-Zhou rows: the printed SD of n = 3 independent runs, the mean's scale s / sqrt 3,
#   Student-t with n - 1 = 2 df, floored at max(read error, r). C-Hou: no printed SD, so A1's same-observable run-to-run
#   reproducibility (run offset 0.447 pH, chi-square 1 df; per-point shape 0.162 pH, 3 df), each integrated over its
#   estimation uncertainty (scaled inverse chi-square, 4-point quantile quadrature), floored per point at max(read
#   0.05, r). C-Mey: a within-run drift (the offset cancels), sqrt(2) r. No same-observable model-discrepancy estimate
#   exists, so absolute levels are falsification-only: T credits on C-Mey's drift, C-Hou's within-run shape (y_t -
#   y_20) and the signs of C-Zhou-ARP's within-study contrasts (7 - 6, 8 - 7; each as a t_2 difference, the two
#   treated as independent); C-Hou's and C-Zhou-ARP's absolute values enter only the adequacy (falsification) test.
#   Spent rows: C-Zhou-Cys.
# - M3 noise: the measurement's own terms plus the lane's within-lab run scatter (sigma_run 0.075 dex as the run
#   offset, and per point: the "between points" form), no omega; sigma_run's estimation uncertainty is not on disk, so
#   it is not integrated (stated). No discrepancy estimate: T credits on within-run ratios and shapes only (route: U4,
#   U1e, U3b; transfer: U5a / U5b PA shapes, U-ARPglc; decay law: U3b, the PA shapes, U-ARPglc); absolute levels (U3,
#   U3a, U5a, U5b) enter the adequacy test only. Spent rows: route U1l; transfer and decay law U4, U1e, U1l.
# - Adequacy: posterior-predictive p of each decisive unit (absolute form included) under the rival's posterior given
#   the spent and decisive data (M1b: Monte Carlo on -log p_pred, 200 replicates; M3: the layer's Gaussian Mahalanobis
#   p; U-ARPglc, censored, has none), Holm at 0.05 across the units; computed only where T crosses its threshold.
# - Leave-one-unit-out: where removing unit u drops T below the null quantile of T without u, the credit rests on u and
#   stands only if T_u alone exceeds the null quantile of T_u (item 4 for that unit alone).
# - Generator worlds: "own" (measurement noise only, r = 0.05), "reproducibility" (adds the same-observable run-to-run
#   term at its point estimate: M1b A1's offset and shape on the unbuffered C-Hou and C-Zhou rows; M3 sigma_run between
#   points), "omega" (adds H0's omega as before). The scorer is the same in every world.

REVISED_RULE_TEXT = (
    "REVISED item 5: noise at the level the compared value varies: SD of independent runs where printed, integrated "
    "over n-1 df; else the same lab class's measured run-to-run reproducibility of the same observable on a repeated "
    "design (run offset + per-point term where shown); instrument repeatability / reading error is a floor, never the "
    "whole noise of a reaction-system observable. A discrepancy term is a variance component (per point, or shared "
    "within a run), estimated on the same observable and lab class from rows on which no shared constant was fitted, "
    "integrated over its estimation uncertainty, the same for every rival, redrawn in calibration; where none exists, "
    "zero is NOT assumed: only values in which a run-level term cancels (within-run ratios, shapes, signs of "
    "within-study contrasts) can credit; absolute levels falsification-only. Across labs: the measured between-lab "
    "spread of the same observable and basis, else not decisive. No transfer from another observable; no bound as an "
    "SD. REVISED M1.1: one statistic T = pooled lppd difference over decisive units, each density p(y_dec | y_spent). "
    "a is credited iff T > the 95th percentile of T under b as truth AND a is adequate on every decisive unit "
    "(posterior-predictive p >= 0.05 after Holm); b by the mirror rule; else undecided, both carried (S3). A rival "
    "falsified under S3 is neither credited nor a null. Each null: parametric bootstrap from the null rival's "
    "spent-conditioned posterior, truths outside the ensemble, redrawing decisive data, conditioning rows, forced "
    "inputs, replicate SDs, the discrepancy component and the form index; >= 1000 replicates where affordable (else "
    "state the count and its MC error). Forms: uncertain experimental facts as discrete nuisance coordinates; "
    "unsourced modelling choices keep required agreement; omega hedges dropped. Print unit contributions and "
    "leave-one-unit-out verdicts; a credit resting on one unit stands only if that unit alone meets item 4.")

REV_WORLDS = ("own", "reproducibility", "omega")
REV_POOL = int(os.environ.get("JUDGE_REV_POOL", "2048"))
REV_R = int(os.environ.get("JUDGE_REV_R", "200"))
REV_PPC = 200
A1_OFF, A1_OFF_DF = 0.447, 1          # M1b prereg v3.7 sec. 20.2: one pair of Zhou 2024 / 2025b runs
A1_SH, A1_SH_DF = 0.162, 3            # the sd (n - 1 = 3) of the four differences, / sqrt 2
R_METER_SET = (0.02, 0.05, 0.10)      # A1's r, its sensitivity 0.02 and the confirmer's 0.10
ZHOU_N = 3                            # ph_acid_submodel.yml sd_basis: "mean +- SD, n = 3"
QUAD_K = 4


def invchi2_quad(s: float, nu: int, k: int = QUAD_K) -> np.ndarray:
    """K variances at the (i - 1/2) / K quantiles of the scaled inverse chi-square posterior of sigma^2 (estimate s,
    nu df): the declared quadrature for 'integrated over its estimation uncertainty'."""
    from scipy.stats import chi2
    u = (np.arange(k) + 0.5) / k
    return nu * s ** 2 / chi2.ppf(1.0 - u, nu)


def invchi2_draw(s: float, nu: int, rng, size=None):
    return nu * s ** 2 / rng.chisquare(nu, size)


def weighted_quantile(x: np.ndarray, w: np.ndarray, q: float) -> Tuple[float, float, float]:
    """(quantile, its Monte Carlo SE, ESS) of a weighted sample; SE = sqrt(q (1 - q) / ESS) / f(quantile), f a weighted
    Gaussian kernel density (Silverman bandwidth)."""
    ok = np.isfinite(x) & (w > 0)
    x, w = x[ok], w[ok]
    if x.size == 0:
        return float("nan"), float("nan"), 0.0
    w = w / w.sum()
    o = np.argsort(x)
    xs, cw = x[o], np.cumsum(w[o])
    qv = float(xs[min(np.searchsorted(cw, q), xs.size - 1)])
    ess = float(1.0 / np.sum(w ** 2))
    mu = float(w @ x)
    sd = float(np.sqrt(max(w @ (x - mu) ** 2, 1e-300)))
    h = 1.06 * sd * max(ess, 2.0) ** (-0.2)
    dens = float(w @ np.exp(-0.5 * ((x - qv) / h) ** 2)) / (h * math.sqrt(2 * math.pi)) if h > 0 else float("nan")
    se = math.sqrt(q * (1 - q) / max(ess, 1.0)) / dens if dens > 0 else float("inf")
    return qv, float(se), ess


LAST_ESS: List[float] = []


def resolved_rows(log_w: np.ndarray, k_max: float = 0.7, ess_min: float = 100.0) -> Tuple[np.ndarray, np.ndarray]:
    """Per row of a (n, S) matrix of log importance weights: PSIS k-hat < k_max and ESS >= ess_min (M1.1's rule,
    predictive_densities.K_HAT_MAX / ESS_MIN). Returns (resolved, k-hat)."""
    from src.kinetic_core import predictive_densities as pdens
    lw = np.atleast_2d(log_w)
    res, ks = np.zeros(lw.shape[0], bool), np.zeros(lw.shape[0])
    for i, row in enumerate(lw):
        fin = np.isfinite(row)
        if fin.sum() < 10:
            ks[i] = float("inf")
            continue
        sm, k = pdens.psis(row[fin] - np.max(row[fin]))
        w = np.exp(sm - np.max(sm))
        ess = float(w.sum() ** 2 / np.sum(w ** 2))
        ks[i] = k
        LAST_ESS.append(ess)
        res[i] = (k < k_max) and (ess >= ess_min)
    return res, ks


def _verdicts(T: Mapping[str, Mapping[str, float]], null_a, null_b, units: Sequence[str], forms: Sequence[str],
              adequate_a, adequate_b) -> Dict[str, Any]:
    """The revised M1.1 verdict of one set. T[form]['full' | '-u' | 'u'] the observed statistics; null_x(form, key)
    -> the (quantile, se, ess) of that key's null when x is the truth (for crediting a: the 95th percentile under b;
    for b: the 5th percentile under a). adequate_x(form) -> bool (Holm-adjusted). Returns the credit, the
    leave-one-unit-out readings and the threshold MC errors."""
    out: Dict[str, Any] = {"credit": "none", "loo": {}, "thr_se": []}
    for side in ("a", "b"):
        ok_all = True
        loo_side = {}
        for f in forms:
            if side == "a":
                q, se, ess = null_b(f, "full", 0.95)
                hit = T[f]["full"] > q
            else:
                q, se, ess = null_a(f, "full", 0.05)
                hit = T[f]["full"] < q
            out["thr_se"].append(se)
            if not (np.isfinite(q) and hit):
                ok_all = False
                break
            if not (adequate_a(f) if side == "a" else adequate_b(f)):
                ok_all = False
                out[f"inadequate_{side}"] = True
                break
            for u in units:
                if side == "a":
                    qm = null_b(f, f"-{u}", 0.95)[0]
                    rests = not (T[f][f"-{u}"] > qm)
                    alone = T[f][u] > null_b(f, u, 0.95)[0]
                else:
                    qm = null_a(f, f"-{u}", 0.05)[0]
                    rests = not (T[f][f"-{u}"] < qm)
                    alone = T[f][u] < null_a(f, u, 0.05)[0]
                loo_side[f"{f}|{u}"] = {"rests_on": bool(rests), "alone_credits": bool(alone)}
                if rests and not alone:
                    ok_all = False
        out["loo"][side] = loo_side
        if ok_all:
            out["credit"] = side
            break
    return out


# ------------------------------------------------------------------------------------------------ M1b, revised

def run_revised_m1b(args) -> Dict[str, Any]:
    m1b_root = Path(args.m1b_root).resolve()
    for p in (str(m1b_root / "scripts" / "generators"), str(m1b_root)):
        if p not in sys.path:
            sys.path.insert(0, p)
    import fit_m1b_ph_balance as G
    from src.kinetic_core import predictive_densities as pdens  # noqa: F401  (resolved_rows imports it too)
    from scipy.stats import t as tdist
    from scipy.special import log_ndtr
    if Path(G.__file__).resolve().parents[2] != m1b_root:
        raise SystemExit("the M1b generator did not import from --m1b-root")
    draws_path = m1b_root / "results" / "validation" / "kinetic_core_m1b_part_i_screen_draws.json.gz"
    ens = G.load_screen_ensembles(draws_path)
    t0 = time.time()
    raw = G.measured_units()
    zc_real_y = np.array(raw["C-Zhou-Cys"].y, float)                       # SPENT rows only
    zmeas = G.zhou_measured()
    zc_sd = np.array([zmeas[f"zhou_cys_pH{x:g}"]["sd"] for x in G.ZHOU_LABELS], float)
    za_sd = np.array([zmeas[f"zhou_arp_pH{x:g}"]["sd"] for x in G.ZHOU_LABELS], float)
    zread = 0.02
    hou_read = float(G.HOU["read_ph"])
    del raw
    sh_var = invchi2_quad(A1_SH, A1_SH_DF)
    off_var = invchi2_quad(A1_OFF, A1_OFF_DF)
    r_set = np.array(R_METER_SET)
    lw_r = np.log(np.full(len(r_set), 1.0 / len(r_set)))
    FORMS = ("forced", "unforced")
    HOU_UNIT = {"forced": "C-Hou", "unforced": "C-Hou-unforced"}
    UNITS = ("mey", "hou_shape", "zhou_signs")
    A_SH = np.zeros((5, 6))
    for k in range(5):
        A_SH[k, k + 1], A_SH[k, 0] = 1.0, -1.0

    # ---- scorer densities, vectorised: data with a leading set axis (n, ...), predictions (S, ...) -> (n, S)
    def ll_mey(within: np.ndarray, f: np.ndarray) -> np.ndarray:
        comps = []
        for r in r_set:
            s = math.sqrt(2.0) * r
            hi, lo = (G.MEYNIER_BOUND - f) / s, (-G.MEYNIER_BOUND - f) / s
            lin = G._log_inside(lo, hi)
            lout = G._log_outside(lo, hi)
            w = within.astype(float)
            comps.append(w @ lin.T + (1.0 - w) @ lout.T)
        return lse(np.stack(comps, 0) + lw_r[:, None, None], axis=0)

    def ll_hou_shape(y: np.ndarray, f: np.ndarray) -> np.ndarray:
        z = y @ A_SH.T
        fz = f @ A_SH.T
        res = z[:, None, :] - fz[None, :, :]
        s1, s2 = np.sum(res ** 2, axis=2), np.sum(res, axis=2) ** 2
        q = s1 - s2 / 6.0
        comps = []
        for ri, r in enumerate(r_set):
            meas2 = max(hou_read, r) ** 2
            for v_sh in sh_var:
                v = meas2 + v_sh
                comps.append(-0.5 * (q / v + 5 * math.log(v) + math.log(6.0) + 5 * math.log(2 * math.pi))
                             + lw_r[ri] - math.log(len(sh_var)))
        return lse(np.stack(comps, 0), axis=0)

    def ll_hou_abs(y: np.ndarray, f: np.ndarray) -> np.ndarray:
        res = y[:, None, :] - f[None, :, :]
        s1, s2 = np.sum(res ** 2, axis=2), np.sum(res, axis=2) ** 2
        comps = []
        for ri, r in enumerate(r_set):
            meas2 = max(hou_read, r) ** 2
            for v_sh in sh_var:
                for o in off_var:
                    v = meas2 + v_sh
                    q = (s1 - o / (v + 6 * o) * s2) / v
                    comps.append(-0.5 * (q + 5 * math.log(v) + math.log(v + 6 * o) + 6 * math.log(2 * math.pi))
                                 + lw_r[ri] - math.log(len(sh_var) * len(off_var)))
        return lse(np.stack(comps, 0), axis=0)

    def zscale(s_rep: np.ndarray, r: float) -> np.ndarray:
        return np.maximum(np.maximum(s_rep / math.sqrt(ZHOU_N), zread), r)

    def ll_zhou_abs(y: np.ndarray, s_rep: np.ndarray, f: np.ndarray) -> np.ndarray:
        comps = []
        for ri, r in enumerate(r_set):
            c = zscale(s_rep, r)                                          # (n, 3)
            x = (y[:, None, :] - f[None, :, :]) / c[:, None, :]
            comps.append(np.sum(tdist.logpdf(x, ZHOU_N - 1) - np.log(c)[:, None, :], axis=2) + lw_r[ri])
        return lse(np.stack(comps, 0), axis=0)

    PAIRS_C = ((1, 0), (2, 1))

    def ll_zhou_signs(y: np.ndarray, s_rep: np.ndarray, f: np.ndarray) -> np.ndarray:
        comps = []
        for ri, r in enumerate(r_set):
            c = zscale(s_rep, r)
            tot = 0.0
            for hi_, lo_ in PAIRS_C:
                sg = np.sign(y[:, hi_] - y[:, lo_])[:, None]
                sd = np.sqrt(c[:, hi_] ** 2 + c[:, lo_] ** 2)[:, None]
                df_ = (f[:, hi_] - f[:, lo_])[None, :]
                tot = tot + tdist.logcdf(sg * df_ / sd, ZHOU_N - 1)
            comps.append(tot + lw_r[ri])
        return lse(np.stack(comps, 0), axis=0)

    # ---- data generators: world "own" / "reproducibility" / "omega", and "scorer" (the bootstrap's own model)
    def gen(E, idx: np.ndarray, world: str, rng) -> Dict[str, np.ndarray]:
        n = idx.size
        out: Dict[str, np.ndarray] = {}
        if world == "scorer":
            r = rng.choice(r_set, size=n)
            v_sh = invchi2_draw(A1_SH, A1_SH_DF, rng, n)
            v_off = invchi2_draw(A1_OFF, A1_OFF_DF, rng, n)
        else:
            r = np.full(n, 0.05)
        # C-Mey
        f = E.f["C-Mey-Cys"][idx]
        var = 2.0 * r[:, None] ** 2 + (E.delta["C-Mey-Cys"]["primary"][idx] ** 2 if world == "omega" else 0.0)
        out["within"] = np.abs(f + rng.standard_normal(f.shape) * np.sqrt(var)) <= G.MEYNIER_BOUND
        # C-Hou (forced run's predictions; the unforced form is scored on the same data)
        f = E.f["C-Hou"][idx]
        meas = np.maximum(hou_read, r)[:, None]
        y = f + rng.standard_normal(f.shape) * meas
        if world == "scorer":
            y = y + rng.standard_normal(f.shape) * np.sqrt(v_sh)[:, None] + rng.standard_normal((n, 1)) * np.sqrt(v_off)[:, None]
        elif world == "reproducibility":
            y = y + rng.standard_normal(f.shape) * A1_SH + rng.standard_normal((n, 1)) * A1_OFF
        elif world == "omega":
            y = y + rng.standard_normal(f.shape) * E.delta["C-Hou"]["primary"][idx]
        out["hou"] = y
        # C-Zhou rows: the mean of n = 3 independent runs; the printed SD redrawn
        for key, unit, sd_print in (("za", "C-Zhou-ARP", za_sd), ("zc", "C-Zhou-Cys", zc_sd)):
            f = E.f[unit][idx]
            if world == "scorer":
                sig2 = invchi2_draw(1.0, ZHOU_N - 1, rng, f.shape) * sd_print[None, :] ** 2
            else:
                sig2 = np.broadcast_to(sd_print[None, :] ** 2, f.shape)
            y = f + rng.standard_normal(f.shape) * np.sqrt(sig2 / ZHOU_N)
            if world == "reproducibility":
                y = y + rng.standard_normal(f.shape) * A1_SH + rng.standard_normal(f.shape) * A1_OFF
            elif world == "omega":
                y = y + rng.standard_normal(f.shape) * E.delta[unit]["primary"][idx]
            out[key] = y
            out[key + "_s"] = np.sqrt(sig2 * rng.chisquare(ZHOU_N - 1, f.shape) / (ZHOU_N - 1))
        return out

    def take(d: Mapping[str, np.ndarray], i) -> Dict[str, np.ndarray]:
        return {k: v[i] for k, v in d.items()}

    def stats(dat: Mapping[str, np.ndarray], Ea, Eb) -> Dict[str, Any]:
        """T and its leave-one-out / alone variants per form, and the resolution of every density, for n sets."""
        lls = {}
        for side, E in (("a", Ea), ("b", Eb)):
            spent = ll_zhou_abs(dat["zc"], dat["zc_s"], E.f["C-Zhou-Cys"])
            per = {"mey": ll_mey(dat["within"], E.f["C-Mey-Cys"]),
                   "zhou_signs": ll_zhou_signs(dat["za"], dat["za_s"], E.f["C-Zhou-ARP"])}
            for form in FORMS:
                per[f"hou_shape|{form}"] = ll_hou_shape(dat["hou"], E.f[HOU_UNIT[form]])
            lls[side] = (spent, per)
        out: Dict[str, Any] = {"T": {}, "resolved": None, "contrib": {}, "diag": {}}
        res_all = None
        for form in FORMS:
            tf = {}
            for side in ("a", "b"):
                spent, per = lls[side]
                u_ll = {"mey": per["mey"], "hou_shape": per[f"hou_shape|{form}"], "zhou_signs": per["zhou_signs"]}
                den = lse(spent, axis=1)
                tot = spent + sum(u_ll.values())
                full = lse(tot, axis=1) - den
                vals = {"full": full}
                for u in UNITS:
                    vals[f"-{u}"] = lse(tot - u_ll[u], axis=1) - den
                    vals[u] = lse(spent + u_ll[u], axis=1) - den
                del LAST_ESS[:]
                rs, ks_ = resolved_rows(tot)
                ess_tot = list(LAST_ESS)
                del LAST_ESS[:]
                rd, kd_ = resolved_rows(spent)
                ess_den = list(LAST_ESS)
                out["diag"][f"{form}|{side}"] = {"k_tot_median": float(np.median(ks_)), "ess_tot_median":
                                                 float(np.median(ess_tot)), "k_spent_median": float(np.median(kd_)),
                                                 "ess_spent_median": float(np.median(ess_den)),
                                                 "resolved_fraction": float(np.mean(rs & rd))}
                r_ = rs & rd
                res_all = r_ if res_all is None else (res_all & r_)
                tf[side] = vals
            out["T"][form] = {k: tf["a"][k] - tf["b"][k] for k in tf["a"]}
        out["resolved"] = res_all
        out["lls"] = lls
        return out

    def spent_ll_over(E, idx_pool: np.ndarray, zc: np.ndarray, zc_s: np.ndarray) -> np.ndarray:
        """log p(spent | theta) for one set's spent rows over the truth draws idx_pool."""
        return ll_zhou_abs(zc[None, :], zc_s[None, :], E.f["C-Zhou-Cys"][idx_pool])[0]

    def adequacy(side_E, dat_one: Mapping[str, np.ndarray], lls_side, i: int, form: str, rng) -> Tuple[bool, Dict]:
        """Holm-adjusted posterior-predictive adequacy on the decisive units (absolute forms included), Monte Carlo
        on -log p_pred with REV_PPC replicates, under the rival's posterior given the spent and decisive data."""
        spent, per = lls_side
        logw = spent[i] + per["mey"][i] + per[f"hou_shape|{form}"][i] + per["zhou_signs"][i]
        w = np.exp(logw - np.max(logw))
        w = w / w.sum()
        S = w.size
        draws = rng.choice(S, size=REV_PPC, p=w)
        sel = {"mey": ("C-Mey-Cys",), "hou": (HOU_UNIT[form],), "za": ("C-Zhou-ARP",)}
        ps = {}
        # replicate data from the scorer's model at the posterior draws
        fake = type("E", (), {})()
        fake.f = {k: side_E.f[k] for k in ("C-Mey-Cys", "C-Hou", "C-Zhou-ARP", "C-Zhou-Cys")}
        fake.f["C-Hou"] = side_E.f[HOU_UNIT[form]]
        fake.delta = side_E.delta
        rep = gen(fake, draws, "scorer", rng)
        lw = np.log(np.maximum(w, 1e-300))
        for u in ("mey", "hou", "za"):
            if u == "mey":
                obs = lse(lw[None, :] + ll_mey(dat_one["within"][None, :], side_E.f["C-Mey-Cys"]), axis=1)
                r_ = lse(lw[None, :] + ll_mey(rep["within"], side_E.f["C-Mey-Cys"]), axis=1)
            elif u == "hou":
                obs = lse(lw[None, :] + ll_hou_abs(dat_one["hou"][None, :], side_E.f[HOU_UNIT[form]]), axis=1)
                r_ = lse(lw[None, :] + ll_hou_abs(rep["hou"], side_E.f[HOU_UNIT[form]]), axis=1)
            else:
                obs = lse(lw[None, :] + ll_zhou_abs(dat_one["za"][None, :], dat_one["za_s"][None, :],
                                                    side_E.f["C-Zhou-ARP"]), axis=1)
                r_ = lse(lw[None, :] + ll_zhou_abs(rep["za"], rep["za_s"], side_E.f["C-Zhou-ARP"]), axis=1)
            ps[u] = float(np.mean(r_ <= obs[0]))         # P(log p_pred(y_rep) <= log p_pred(y_obs))
        rej = holm_reject(ps)
        return (not any(rej.values())), ps

    rng = np.random.default_rng(M1B_SEED + 101)
    halves = {}
    for rv in G.RIVAL_NAMES:
        S = ens[rv].size
        halves[rv] = {"truth": np.arange(0, S, 2), "score": np.arange(1, S, 2)}
    score_E = {rv: ens[rv].subset(halves[rv]["score"]) for rv in G.RIVAL_NAMES}
    full_E = ens

    # real spent rows -> each rival's evaluation-truth weights over its held-out draws
    def truth_weights(rv: str) -> np.ndarray:
        lw = spent_ll_over(full_E[rv], halves[rv]["truth"], zc_real_y, zc_sd)
        w = np.exp(lw - np.max(lw))
        return w / w.sum()

    results: Dict[str, Any] = {}
    want = [x for x in args.m1b_pairs.split(",") if x]
    for a, b in G.PAIRS:
        key = f"{a}_vs_{b}"
        if want and key not in want:
            continue
        t1 = time.time()
        Ea, Eb = score_E[a], score_E[b]
        # ---- bootstrap pools (the scorer's own model), one per null rival
        pools = {}
        for null in (a, b):
            tidx = halves[null]["truth"]
            pick = rng.integers(0, tidx.size, REV_POOL)
            dat = gen(full_E[null], tidx[pick], "scorer", rng)
            st = stats(dat, Ea, Eb)
            pools[null] = {"pick": pick, "T": st["T"], "resolved": st["resolved"]}
        cost_pool = time.time() - t1

        def null_fn(null: str, spent_lw_over_truth: np.ndarray, resolved_only: bool = True):
            P = pools[null]
            w = np.exp(spent_lw_over_truth - np.max(spent_lw_over_truth))[P["pick"]]
            if resolved_only:
                w = w * P["resolved"]

            def f(form, k, q):
                return weighted_quantile(P["T"][form][k], w, q)
            return f

        world_out = {}
        diag_w: Dict[str, Any] = {}
        for world in REV_WORLDS:
            rows = []
            for truth in (a, b):
                tw = truth_weights(truth)
                tidx = halves[truth]["truth"][rng.choice(halves[truth]["truth"].size, size=REV_R, p=tw)]
                dat = gen(full_E[truth], tidx, world, rng)
                st = stats(dat, Ea, Eb)
                diag_w[f"{world}|truth_{truth}"] = st["diag"]
                print(f"  resolution {world} truth {truth}: {json.dumps(st['diag'])}", flush=True)
                for i in range(REV_R):
                    row: Dict[str, Any] = {"truth": truth}
                    T_i = {form: {k: float(v[i]) for k, v in st["T"][form].items()} for form in FORMS}
                    la = spent_ll_over(full_E[a], halves[a]["truth"], dat["zc"][i], dat["zc_s"][i])
                    lb = spent_ll_over(full_E[b], halves[b]["truth"], dat["zc"][i], dat["zc_s"][i])
                    one = take(dat, i)
                    ad_a = lambda form, one=one, i=i: adequacy(Ea, one, st["lls"]["a"], i, form, rng)[0]  # noqa: E731
                    ad_b = lambda form, one=one, i=i: adequacy(Eb, one, st["lls"]["b"], i, form, rng)[0]  # noqa: E731
                    v_est = _verdicts(T_i, null_fn(a, la, False), null_fn(b, lb, False), UNITS, FORMS, ad_a, ad_b)
                    row["credit_est"] = v_est["credit"]
                    row["loo_est"] = v_est["loo"]
                    row["contrib_all"] = {u: T_i["forced"][u] for u in UNITS}
                    if not st["resolved"][i]:
                        row["credit"] = "not computable"
                        rows.append(row)
                        continue
                    v = _verdicts(T_i, null_fn(a, la), null_fn(b, lb), UNITS, FORMS, ad_a, ad_b)
                    row.update({"credit": v["credit"], "T": T_i["forced"]["full"],
                                "contrib": {u: T_i["forced"][u] for u in UNITS},
                                "loo": v["loo"], "thr_se": [x for x in v["thr_se"] if np.isfinite(x)],
                                "thr_ess": weighted_quantile(pools[b]["T"]["forced"]["full"],
                                                             np.exp(lb - lb.max())[pools[b]["pick"]], 0.95)[2]})
                    rows.append(row)
            world_out[world] = summarise_revised(rows, a, b, UNITS)
            world_out[world]["estimated_diagnostic"] = summarise_revised(rows, a, b, UNITS, key="credit_est")
            print(f"m1b revised {key} {world}: credit {a} power {world_out[world]['credit_a']['power']:.3f} false "
                  f"{world_out[world]['credit_a']['false_credit']:.3f}; not computable "
                  f"{world_out[world]['not_computable']} ({time.time() - t1:.0f} s)", flush=True)
        results[key] = {"a": a, "b": b, "worlds": world_out, "resolution": diag_w,
                        "pool": {n: {"size": REV_POOL, "resolved_fraction": float(np.mean(p["resolved"]))}
                                 for n, p in pools.items()},
                        "seconds": round(time.time() - t1, 1), "pool_seconds": round(cost_pool, 1)}
    return {"case": "M1b part (i) REVISED", "pairs": results, "R": REV_R, "pool": REV_POOL, "ppc": REV_PPC,
            "seconds": round(time.time() - t0, 1),
            "split": "screen ensembles by draw parity: even = truth pool, odd = scoring ensemble",
            "forced_inputs": ("Hou's measured ARP(t) forcing is NOT redrawn: each redraw is a re-integration of the "
                              "forced pot for every scoring draw (~640 draws x 2 rivals x ~1 - 2 s per bootstrap "
                              "replicate, i.e. ~0.5 h per replicate): not computable here"),
            "noise": {"A1_offset": [A1_OFF, A1_OFF_DF], "A1_shape": [A1_SH, A1_SH_DF], "quadrature_K": QUAD_K,
                      "r_set": list(R_METER_SET), "zhou_n": ZHOU_N},
            "inputs": {"draws_sha256": sha256_file(draws_path), "generator_sha256": sha256_file(Path(G.__file__))}}


def summarise_revised(rows: Sequence[Mapping[str, Any]], a: str, b: str, units: Sequence[str],
                      key: str = "credit") -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for side, hit, other in (("a", a, b), ("b", b, a)):
        ka = sum(1 for r in rows if r["truth"] == hit and r.get(key) == side)
        kb = sum(1 for r in rows if r["truth"] == other and r.get(key) == side)
        na = sum(1 for r in rows if r["truth"] == hit)
        nb = sum(1 for r in rows if r["truth"] == other)
        out[f"credit_{side}"] = {"rival": hit, "power": ka / na, "power_k": ka, "n": na, "power_wilson95": wilson(ka, na),
                                 "false_credit": kb / nb, "false_k": kb, "false_wilson95": wilson(kb, nb),
                                 "p_false_above_nominal": binom_above(kb, nb)}
    out["not_computable"] = sum(1 for r in rows if r.get(key) == "not computable")
    loo_key, con_key = ("loo", "contrib") if key == "credit" else ("loo_est", "contrib_all")
    se = [x for r in rows for x in r.get("thr_se", [])]
    out["threshold_mc_se"] = {"median": float(np.median(se)) if se else None,
                              "q90": float(np.quantile(se, 0.9)) if se else None}
    ess = [r["thr_ess"] for r in rows if "thr_ess" in r]
    out["threshold_ess_median"] = float(np.median(ess)) if ess else None
    loo = {}
    for side in ("a", "b"):
        for u in units:
            # every set where T crossed the side's threshold and the unit-removal was examined
            ex = [r for r in rows if r.get(loo_key) and r[loo_key].get(side)]
            rests = [r for r in ex if any(v["rests_on"] for k, v in r[loo_key][side].items() if k.endswith(f"|{u}"))]
            alone = [r for r in rests if all(v["alone_credits"] for k, v in r[loo_key][side].items()
                                             if k.endswith(f"|{u}") and v["rests_on"])]
            loo[f"{side}|{u}"] = {"examined": len(ex), "rests_on_it": len(rests), "of_which_it_alone_credits": len(alone)}
    out["leave_one_unit_out"] = loo
    for t in (a, b):
        sel = [r for r in rows if r["truth"] == t and r.get(con_key)]
        out[f"contrib_median_truth_{t}"] = ({u: float(np.median([r[con_key][u] for r in sel])) for u in units}
                                            if sel else {})
    return out


# ------------------------------------------------------------------------------------------------ M3, revised

def _load_m3(args):
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    gen_dir = str(ROOT / "scripts" / "generators")
    if gen_dir not in sys.path:
        sys.path.insert(0, gen_dir)
    import fit_m3_strecker_layer as G
    from src import m3_strecker_layer as L
    data = L.load_data()
    real = G.real_obs(data)
    rec_sets = json.loads((ROOT / "results" / "validation" / "kinetic_core_m3_phase1_power.json").read_text())["draw_sets"]
    spec_args = argparse.Namespace(score_feasible=M3_SCORE_FEASIBLE, truth_feasible=M3_TRUTH_FEASIBLE, seed=M3_SET_SEED)
    specs = G.power_set_specs(spec_args, tuple(G.POWER_RIVALS))

    class _U(pickle.Unpickler):
        def find_class(self, module, name):
            if name == "DrawSet" and module in ("__main__", "fit_m3_strecker_layer"):
                return G.DrawSet
            return super().find_class(module, name)

    loaded = []
    for p in sorted(Path(args.m3_cache).glob("drawset_*.pkl")):
        with p.open("rb") as fh:
            loaded.append(_U(fh).load())
    sets = {}
    for sp in specs:
        hit = [ds for ds in loaded if ds.rival == sp["rival"] and ds.seed == sp["seed"] and ds.form == sp["form"]]
        name, fname = sp["key"]
        want = rec_sets[f"{name}|{fname}"]
        if len(hit) != 1 or (hit[0].n, hit[0].n_drawn) != (want["feasible"], want["drawn"]):
            raise SystemExit(f"draw set {sp['key']} does not match the record")
        sets[(name, fname)] = hit[0]
    return G, L, data, real, sets


def run_revised_m3(args) -> Dict[str, Any]:
    from dataclasses import replace
    G, L, data, real, sets = _load_m3(args)
    t0 = time.time()
    FORMS = tuple(G.POWER_FORMS)
    PAIRS = {"route": ("R1", "R2"), "transfer": ("R1", "R3"), "decay_law": ("R4", "R1drift")}
    SPENT = {"route": ("U1l",), "transfer": ("U4", "U1e", "U1l"), "decay_law": ("U4", "U1e", "U1l")}
    CRED = {"route": ("U4", "U1e", "U3b"), "transfer": ("U5a_PA_shape", "U5b_PA_shape", "UARP"),
            "decay_law": ("U3b", "U5a_PA_shape", "U5b_PA_shape", "UARP")}
    ADEQ = {"route": ("U4", "U1e", "U3b", "U3a"), "transfer": ("U5a", "U5b", "U5a_PA_shape", "U5b_PA_shape"),
            "decay_law": ("U3", "U5a", "U5b", "U3b")}
    pool_n = int(os.environ.get("JUDGE_REV_POOL_M3", "1000"))

    def scorer_form(fname: str):
        return replace(G.POWER_FORMS[fname], omega_scale=0.0, run_between_points=True)

    WORLD_FORM = {"own": (replace(L.PRIMARY, run_between_points=False), True),
                  "reproducibility": (replace(L.PRIMARY, run_between_points=True), True),
                  "omega": (L.PRIMARY, False)}

    def point_cov(unit, world):
        fw, mo = WORLD_FORM[world]
        return L.unit_covariance(replace(unit, design=np.eye(unit.design.shape[1])), fw, measurement_only=mo)

    # consistency of shared points (U1e / U1l are U4's Fig. 1 PA points; the PA shapes are U5's PA points)
    ts = sets[("R1", "truth")]
    assert np.allclose(ts.points["U1e"], ts.points["U4"][:, 6:13]), "U1 points are not U4's Fig. 1 PA points"
    for u5, sh in (("U5a", "U5a_PA_shape"), ("U5b", "U5b_PA_shape")):
        assert np.allclose(ts.points[sh], ts.points[u5][:, :ts.points[sh].shape[1]]), f"{sh} is not {u5}'s PA"

    def uarp_synth(k_1min, rng, world):
        y = G.synth_uarp(data, k_1min, rng)
        if world == "reproducibility":
            y = y + L.SIGMA_RUN_DEX * math.log(10.0) * math.sqrt(2.0) * rng.standard_normal(y.shape)
        return y

    def make_obs(tset, j, rng, world):
        units = dict(real.units)
        u4 = real.units["U4"]
        p4 = rng.multivariate_normal(tset.points["U4"][j], point_cov(u4, world))
        units["U4"] = replace(u4, y=u4.design @ p4)
        for u in ("U1e", "U1l"):
            units[u] = replace(real.units[u], y=real.units[u].design @ p4[6:13])
        for u5, sh in (("U5a", "U5a_PA_shape"), ("U5b", "U5b_PA_shape")):
            uu = real.units[u5]
            p5 = rng.multivariate_normal(tset.points[u5][j], point_cov(uu, world))
            units[u5] = replace(uu, y=p5)
            n_t = real.units[sh].design.shape[1]
            units[sh] = replace(real.units[sh], y=real.units[sh].design @ p5[:n_t])
        f3 = tset.points["U3"][j]
        fw, mo = WORLD_FORM[world]
        cov3 = L.unit_covariance(L.u3_unit(10.0 ** f3), fw, measurement_only=mo)
        u3 = L.u3_unit(10.0 ** (f3 + rng.multivariate_normal(np.zeros(f3.size), cov3)))
        y_lin = np.where(u3.censored, 0.5 * L.U3_CENSOR, 10.0 ** u3.y)
        err = data.arp_t2_err["argon"][0]
        sd = math.sqrt(err ** 2 + L.S_REP_DEX ** 2 + L.SIGMA_RUN_DEX ** 2
                       + (L.OMEGA_DEX ** 2 if world == "omega" else 0.0)
                       + (L.SIGMA_RUN_DEX ** 2 if world == "reproducibility" else 0.0))
        pa = 10.0 ** (tset.points["U4"][j][0] + sd * rng.standard_normal())
        units.update({"U3": u3, "U3a": L.u3a_unit(y_lin, pa, err), "U3b": L.u3b_unit(y_lin)})
        return G.Obs(units, uarp_synth(tset.k_glc[j], rng, world))

    def ll(ds, unit, obs, fv):
        return G.loglik(ds, unit, obs, fv, data)

    def stats_one(cname, a_name, b_name, obs) -> Dict[str, Any]:
        out: Dict[str, Any] = {"T": {}, "resolved": True, "skipped_forms": [], "lw": {}, "diag": {}}
        for fname in FORMS:
            fv = scorer_form(fname)
            A, B = sets[(a_name, fname)], sets[(b_name, fname)]
            if A.n == 0 or B.n == 0:
                out["skipped_forms"].append(fname)
                continue
            tf, keep = {}, {}
            for side, ds in (("a", A), ("b", B)):
                spent = sum(ll(ds, u, obs, fv) for u in SPENT[cname])
                lu = {u: ll(ds, u, obs, fv) for u in CRED[cname]}
                tot = spent + sum(lu.values())
                den = G.lse(spent)
                vals = {"full": G.lse(tot) - den}
                for u in CRED[cname]:
                    vals[f"-{u}"] = G.lse(tot - lu[u]) - den
                    vals[u] = G.lse(spent + lu[u]) - den
                del LAST_ESS[:]
                r1, k1 = resolved_rows(tot[None, :])
                r2, k2 = resolved_rows(spent[None, :])
                e1, e2 = (LAST_ESS + [float("nan"), float("nan")])[:2]
                out["diag"][f"{fname}|{side}"] = [round(float(k1[0]), 3), round(e1, 1), round(float(k2[0]), 3),
                                                  round(e2, 1)]
                out["resolved"] = bool(out["resolved"] and r1[0] and r2[0])
                tf[side] = vals
                keep[side] = tot
            out["T"][fname] = {k: tf["a"][k] - tf["b"][k] for k in tf["a"]}
            out["lw"][fname] = keep
        return out

    def adequate(cname, ds, obs, fname, lw_post) -> bool:
        fv = scorer_form(fname)
        ps = {}
        for u in ADEQ[cname]:
            uo = obs.units[u]
            mu, cov = G.predictive_moments(ds, uo, ds.points[u], lw_post, fv)
            ps[u] = G.adequacy_p(uo, mu, cov)
        return not any(holm_reject(ps).values())

    rng = np.random.default_rng(M3_SEED + 101)
    out: Dict[str, Any] = {}
    contrasts = [c for c in args.m3_contrasts.split(",") if c]
    for cname in contrasts:
        t1 = time.time()
        a_name, b_name = PAIRS[cname]
        fv0 = scorer_form("primary")
        pools = {}
        for null in (a_name, b_name):
            tset = sets[(null, "truth")]
            pick = rng.integers(0, tset.n, pool_n)
            Ts, res = [], []
            for j in pick:
                st = stats_one(cname, a_name, b_name, make_obs(tset, int(j), rng, "reproducibility"))
                Ts.append(st["T"])
                res.append(st["resolved"])
            pools[null] = {"pick": pick, "resolved": np.array(res, bool),
                           "T": {f: {k: np.array([t[f][k] for t in Ts]) for k in Ts[0][f]} for f in Ts[0]}}
            print(f"m3 revised {cname}: pool {null} {pool_n} ({time.time() - t1:.0f} s), resolved "
                  f"{np.mean(res):.3f}", flush=True)
        pool_s = time.time() - t1

        def null_fn(null, obs, resolved_only=True):
            tset = sets[(null, "truth")]
            lw = sum(ll(tset, u, obs, fv0) for u in SPENT[cname])
            P = pools[null]
            w = np.exp(lw - np.max(lw))[P["pick"]]
            if resolved_only:
                w = w * P["resolved"]

            def f(form, k, q):
                return weighted_quantile(P["T"][form][k], w, q)
            return f

        world_out = {}
        for world in REV_WORLDS:
            rows = []
            for truth in (a_name, b_name):
                tset = sets[(truth, "truth")]
                lw_real = sum(ll(tset, u, real, fv0) for u in SPENT[cname])
                tw = np.exp(lw_real - np.max(lw_real))
                tw = tw / tw.sum()
                for i in range(args.m3_reps):
                    j = int(rng.choice(tset.n, p=tw))
                    obs = make_obs(tset, j, rng, world)
                    st = stats_one(cname, a_name, b_name, obs)
                    row: Dict[str, Any] = {"truth": truth, "diag": st["diag"]}
                    if i < 2:
                        print(f"  diag {cname} {world} {truth} [k_tot, ess_tot, k_spent, ess_spent]: "
                              f"{json.dumps(st['diag'])}", flush=True)
                    forms = list(st["T"])
                    ad_a = lambda f, st=st, obs=obs: adequate(cname, sets[(a_name, f)], obs, f, st["lw"][f]["a"])  # noqa
                    ad_b = lambda f, st=st, obs=obs: adequate(cname, sets[(b_name, f)], obs, f, st["lw"][f]["b"])  # noqa
                    v_est = _verdicts(st["T"], null_fn(a_name, obs, False), null_fn(b_name, obs, False), CRED[cname],
                                      forms, ad_a, ad_b)
                    row["credit_est"] = v_est["credit"]
                    row["loo_est"] = v_est["loo"]
                    row["contrib_all"] = {u: st["T"]["primary"][u] for u in CRED[cname]} if "primary" in st["T"] else {}
                    if not st["resolved"]:
                        row["credit"] = "not computable"
                        row["contrib_est"] = {u: st["T"]["primary"][u] for u in CRED[cname]} if "primary" in st["T"] else {}
                        rows.append(row)
                        continue
                    na_, nb_ = null_fn(a_name, obs), null_fn(b_name, obs)
                    v = _verdicts(st["T"], na_, nb_, CRED[cname], forms, ad_a, ad_b)
                    row.update({"credit": v["credit"], "T": st["T"]["primary"]["full"],
                                "contrib": {u: st["T"]["primary"][u] for u in CRED[cname]}, "loo": v["loo"],
                                "thr_se": [x for x in v["thr_se"] if np.isfinite(x)],
                                "skipped_forms": st["skipped_forms"]})
                    rows.append(row)
            world_out[world] = summarise_revised(rows, a_name, b_name, CRED[cname])
            world_out[world]["estimated_diagnostic"] = summarise_revised(rows, a_name, b_name, CRED[cname],
                                                                         key="credit_est")
            world_out[world]["skipped_forms"] = sorted({f for r in rows for f in r.get("skipped_forms", [])})
            print(f"m3 revised {cname} {world}: credit {a_name} power {world_out[world]['credit_a']['power']:.3f} "
                  f"false {world_out[world]['credit_a']['false_credit']:.3f}; credit {b_name} power "
                  f"{world_out[world]['credit_b']['power']:.3f} false {world_out[world]['credit_b']['false_credit']:.3f};"
                  f" not computable {world_out[world]['not_computable']} ({time.time() - t1:.0f} s)", flush=True)
        out[cname] = {"pair": [a_name, b_name], "worlds": world_out, "seconds": round(time.time() - t1, 1),
                      "pool_seconds": round(pool_s, 1),
                      "pool": {n: {"size": pool_n, "resolved_fraction": float(np.mean(p["resolved"]))}
                               for n, p in pools.items()},
                      "units": {"spent": SPENT[cname], "credit": CRED[cname], "adequacy": ADEQ[cname]}}
    return {"case": "M3 phase 1 REVISED", "pairs": out, "R": args.m3_reps, "pool": pool_n,
            "seconds": round(time.time() - t0, 1),
            "forced_inputs": ("the measured dicarbonyl / ARP-Phe courses forcing the layer are NOT redrawn: each redraw "
                              "re-evaluates the layer for every scoring draw (~2.5 min per 4 096-draw set and form, "
                              "16 sets per replicate: ~40 min per bootstrap replicate): not computable here"),
            "noise": ("measurement + sigma_run 0.075 dex as run offset and per point (between-points form), no omega; "
                      "sigma_run not integrated (its estimation uncertainty is not on disk)")}


# ================================================================================================ the record

def _row(t: Mapping[str, Any]) -> str:
    return (f"{t['power']:.3f} [{t['power_wilson95'][0]:.3f}, {t['power_wilson95'][1]:.3f}] | "
            f"{t['false_credit']:.3f} [{t['false_wilson95'][0]:.3f}, {t['false_wilson95'][1]:.3f}] | "
            f"{t['fisher_p_power_gt_false']:.2g}")


SINGLE = ["baseline"] + list(COMPONENTS) + ["all"]


def row_label(cfg: str) -> str:
    if cfg in COMPONENTS:
        return COMPONENT_LABEL[cfg]
    for c in COMPONENTS:
        if cfg == leave_one_in(c):
            return f"all but {COMPONENT_LABEL[c].split(' ', 1)[0]} (restores {c})"
    return cfg


def leave_one_in(c: str) -> str:
    return "+".join(k for k in COMPONENTS if k != c)


def attribution(tab: Mapping[str, Any], suffix: str = "") -> Dict[str, Any]:
    """Sec. ATTRIBUTION: (a) single ablation vs baseline, (b) restoring the component inside 'all', and main effects
    (mean power with the component ablated minus without, over the configurations present)."""
    base = tab.get("baseline" + suffix)
    allc = tab.get("all" + suffix)
    out = {}
    for c in COMPONENTS:
        e: Dict[str, Any] = {}
        single = tab.get(c + suffix)
        loi = tab.get(leave_one_in(c) + suffix)
        if base and single:
            e["single_power"] = single["power"]
            e["single_p"] = fisher_greater(single["power_k"], base["power_k"], single["n_a"], base["n_a"])
        if allc and loi:
            e["restored_power"] = loi["power"]
            e["restored_p"] = fisher_greater(allc["power_k"], loi["power_k"], allc["n_a"], loi["n_a"])
        on = [v["power"] for k, v in tab.items() if k.endswith(suffix) and "config" in v and v["config"][c]
              and k.count("|") == suffix.count("|")]
        off = [v["power"] for k, v in tab.items() if k.endswith(suffix) and "config" in v and not v["config"][c]
               and k.count("|") == suffix.count("|")]
        if on and off:
            e["main_effect"] = float(np.mean(on) - np.mean(off))
        e["carries"] = bool(e.get("single_p", 1.0) < NOMINAL or e.get("restored_p", 1.0) < NOMINAL)
        out[c] = e
    return out


def merge(paths: Sequence[str], out_dir: Optional[str] = None) -> int:
    sys.path.insert(0, str(ROOT))
    from src import data_paths, provenance
    parts = [json.loads(Path(p).read_text()) for p in paths]
    rec: Dict[str, Any] = {"artifact": STEM, "generated_by": GENERATED_BY, "kind": "diagnostic (no rule, prereg or "
                           "src/kinetic_core change; nothing ships)", "nominal_false_credit": NOMINAL,
                           "components": COMPONENT_LABEL, "cases": {}}
    for p in parts:
        key = ("m1b" if p["case"].startswith("M1b") else "m3") + ("_revised" if p["case"].endswith("REVISED") else "")
        rec["cases"][key] = p
    # attribution and control
    findings: Dict[str, Any] = {}
    if "m1b" in rec["cases"]:
        for pair, v in rec["cases"]["m1b"]["pairs"].items():  # first-round parts only (the revised are rendered apart)
            tab = v["worlds"]["own"]["table"]
            om = v["worlds"]["omega"]["table"]
            findings[f"m1b|{pair}"] = {
                "attribution": attribution(tab),
                "all_false_credit_ok_own": tab["all"]["p_false_above_nominal"] >= NOMINAL,
                "all_false_credit_ok_omega_world": om["all"]["p_false_above_nominal"] >= NOMINAL}
    if "m3" in rec["cases"]:
        for cname, v in rec["cases"]["m3"]["pairs"].items():
            a_name, b_name = v["worlds"]["own"]["pair"]  # first round
            for hit in (a_name, b_name):
                sfx = f"|credit_{hit}"
                tab = {k: x for k, x in v["worlds"]["own"]["table"].items() if k.endswith(sfx)}
                om = {k: x for k, x in v["worlds"]["omega"]["table"].items() if k.endswith(sfx)}
                findings[f"m3|{cname}|credit_{hit}"] = {
                    "attribution": attribution(tab, sfx),
                    "all_false_credit_ok_own": tab["all" + sfx]["p_false_above_nominal"] >= NOMINAL,
                    "all_false_credit_ok_omega_world": om["all" + sfx]["p_false_above_nominal"] >= NOMINAL}
    rec["findings"] = findings
    rec["reading"] = READING
    rec["provenance"] = provenance.provenance_block(STEM, generated_by=GENERATED_BY)
    out = (Path(out_dir) if out_dir else data_paths.VALIDATION_DIR) / f"{STEM}.json"
    out.write_text(provenance.dumps_compact(_jsonable(rec), default=float) + "\n", encoding="utf-8")
    out.with_suffix(".md").write_text(render(_jsonable(rec)), encoding="utf-8")
    print(f"wrote {out} and .md")
    return 0


#: the confirmation run's reading, written AFTER its numbers (not blind)
REVISED_READING = """\
*Written after the numbers (not blind).*

**1. As worded, the revised rule is NOT COMPUTABLE on the cases' draw ensembles.** M1b: 400 of 400 sets per pair and
world are not computable in four of five pairs; B3 vs B0 is computable in 33 - 155 of 400, with rule power 0.000 - 0.010
and false credit 0.000. M3: 400 of 400 in every pair and world. The bootstrap pools are 0 - 20 % resolved. The cause is
the one this record found in round one, now without omega to hide it. At the noise item 5 prescribes, the decisive
likelihoods are sharp relative to the prior draw ensembles: joint ESS 2 - 25 of 320 - 4 096 draws (M3 route: spent ESS
~2 000, joint 4 - 13). In M1b, conditioning on the three spent C-Zhou-Cys rows already leaves B1 / B13 / B2 an ESS of
3 - 90. The gate is right to refuse: at ESS ~4 the estimates are biased, not just noisy. For example, the TRUE B1 loses
its own C-Hou shape (median contribution -0.84 to -1.25 nats).

**2. Cost of re-estimating (the rule's alternative to "not computable").**
- Per set, every rival, form, numerator and denominator needs a density with ESS >= 100.
- M3: the layer is cheap, but tuned per-set densities (the reusable SMC, 4.5 - 10 min each, 16 per set) cost ~1 - 3 h
  per set, i.e. ~4 000 - 8 000 core-hours per pair at R = 200 plus pools of 1 000.
- M3, the affordable route: one-time, ~10 - 20x larger draw sets per rival and form (~40 000 - 80 000 feasible draws;
  ~2.5 min per 4 096) cost ~10 - 20 core-hours for all 20 rival-form sets. That would lift joint ESS to ~50 - 250, so
  most sets should become computable (not shown here).
- M1b: engine-bound. The screen cost ~29 core-s per draw. A ~10x larger ensemble is ~250 core-hours, and B1 / B13
  need ~30x for their spent conditioning alone. Spent-conditioned proposals would need new engine runs per synthetic
  set (~5 core-hours each). On one machine this is not computable for the confirmation.
- Forced inputs (Hou's ARP(t); M3's measured courses) are not redrawn in either case: each redraw re-runs the engine or
  the layer for every scoring draw.

**3. The same rule with every estimate read as estimated (a DIAGNOSTIC, not the rule; it credits nothing).**
- **False credit is controlled everywhere.** Reproducibility world: max 0.020 (M1b) and 0.010 (M3). Own world: max
  0.035. Omega world: <= 0.005. Across all worlds every exact one-sided p >= 0.88 against 0.05. The per-set
  calibrated threshold with the Holm adequacy gate holds false credit at or below 0.05 + MC error, even under
  H0-sized misspecification, where every round-one configuration failed (0.155 - 0.545).
- **Power** (reproducibility world):
  - M3 route: 0.70 (crediting R1) and 0.76 (crediting R2); own world 0.86 / 0.90.
  - M3 transfer: <= 0.04. M3 decay law: <= 0.005.
  - M1b: <= 0.06 in every pair (own world <= 0.17).
  - Omega world: <= 0.025 everywhere. The rule stays silent rather than crediting. Which gate silences it (the
    threshold or the adequacy test) is not separated here.
- **Threshold MC error** (rule-computable sets only): ~0.07 - 0.09 nats in M1b B3 vs B0 (weights' ESS ~850 - 1 400,
  pool 2 048); 0.71 nats in B2 vs B0's three computable omega-world sets.
- **Leave-one-unit-out.**
  - M1b: C-Zhou-ARP's dominance is gone. Restricted by item 5 to the signs of its within-study contrasts, it
    contributes ~0.01 nats. Almost every crossing rests on C-Hou's within-run shape; C-Mey's drift adds
    |T_u| <= 0.4 nats (reproducibility world).
  - M3 route: U4 carries 25 - 35 % of the crossings, and in 69 - 90 % of those U4 alone crosses its own null. U1e
    carries 9 - 21 %, U3b 2 - 7 % (own and reproducibility worlds).

**4. What this means for the proposal.** The revised judge is calibrated: its false credit holds in every world tested,
including misspecification. Its power is real where the data are informative (M3 route ~0.7 - 0.9) and honest where
they are not (M1b, M3 transfer and decay law ~0). But no verdict can be issued under it until the density estimator
resolves at item 5's noise. For M3 that is a bounded, affordable step (bigger draw sets). For M1b it is not affordable
on one machine, and the M1b confirmation stands at "not computable", not at a power figure. The diagnostic numbers
bound what a resolved estimator should give; they must be re-confirmed with resolved densities before any use. Nothing
here reopens M1b part (i) (A1) or applies to any case.
"""

#: written AFTER the numbers (not blind); the record's interpretation and the proposed principled revision
READING = """\
*Written after the numbers (not blind). Own world unless stated; "all" = the five ablations together with Holm.*

**1. The judge is deaf by construction, through two universal blockers, each sufficient on its own.** No single
ablation gives power in any pair except M3 route crediting R2 (the primary form only: 0.95); every other pair that gains
power needs (i) AND (v) lifted together.
- **(i) omega.** In M1b the +-omega per-species corner term (delta-pH 0.4 - 3.8 pH units) makes the TRUE rival's
  density lower than B0's on its own data: truth B1, C-Hou d = -4.1 nats against a 0.09 band, C-Zhou-ARP -2.9; truth
  B13, -4.05 / -2.6; every decisive unit reads B0 in 100 % of the sets. This is the independent review's finding seen
  from the data side: the corner term grows with the produced titratable species, so it penalises the richer rival
  structurally. Replaced by the per-point measurement noise (printed replicate SD, else max(read error, meter 0.05
  pH)), the true B1 wins C-Hou (d +7.5, band 2.4) and C-Zhou-ARP (d +28.6, band 11.8). In M3 omega makes the
  unconditioned U4 veto the true R1 in 98.5 % of route sets.
- **(v) "unresolved density reads as a tie".** The resolution rate is 0.00 for every M1b unit in every configuration,
  and for M3's route and transfer primaries: the precision rule ESS >= 2 / (0.15 sd_rep)^2 asks for more effective
  draws than the 640 - 4 096 prior draws hold whenever the band is narrow, and at the data's own noise the prior draws
  are too few for ESS >= 100. In a power check this is an estimator limit, not evidence; the score was to use tuned
  proposals, so the power check never described the scorer.
- **Case-specific blockers.** M1b: the conjunction (C-Mey (Cys) AND C-Hou): C-Mey's censored statement barely
  separates any pair (|d| <= 0.12 nats), so requiring it to favour the rival removes most power (all but (iii): 0.00 -
  0.12). Nuisance conditioning matters for B1 vs B0 (all but (ii) 0.00 -> all 0.97), B13 vs B0 (0.12 -> 1.00) and
  B13 vs B3 (0.10 -> 0.98); required agreement for B1 (0.01 -> 0.97), B13 (0.76 -> 1.00) and B13 vs B3 (0.81 -> 0.98).
  M3: required
  agreement carries route-crediting-R2 (R2 has no feasible draw under the 20 mL headspace form, so four-form agreement
  can never credit R2) and transfer (with the (F) / (T) agreement). M3's decay law has zero power in every
  configuration: its block tie band is heavy-tailed (sd_rep up to ~2e4 nats, from the censored U-ARPglc and U3
  replicates); not resolved here.

**2. False credit is NOT controlled by the all-together set.** At the data's own noise "all" exceeds nominal in M1b B2
vs B0 (0.165) and B13 vs B3 (0.195), and in M3 transfer crediting R1 (0.120). Two causes: Holm's "any unit" lets one
unit carry the credit, and every tie band is a noise-only band at FIXED predictive distributions, which ignores the
sampling variability of the conditioning rows and, for the censored C-Mey unit, of the outcome itself (truth B3: C-Mey
reads B13 in 51 % of sets at d 0.03 against a 0.03 band). The pooled variant (one decisive statistic: the summed d over
the decisive units with one band) keeps false credit at or below nominal in all five M1b pairs (max 0.055, exact p
0.42) at power 0.975 / 1.00 / 0.38 / 1.00 / 0.975. **Under misspecification (the omega world) every configuration
that drops omega gives false credit far above nominal**: M1b 0.155 - 0.545, M3 transfer 0.18 - 0.235 (route 0.06 /
0.015). M1b's omega world uses the corner delta-pH itself as the generator's discrepancy (an extreme stress); M3's is
the prereg's per-point 0.419 dex. So removing omega is safe only where the data carry no model discrepancy of that
size, and that has to be measured, not assumed.

**3. What follows (a proposal; nothing here is applied).** The deafness is not "the data cannot decide": it is the
judge. The cure is not deleting terms but calibrating the judge. Proposed principled revision, to be frozen by
principle and confirmed independently before any use:
- *Sec. 7 item 5 (wording):* "Tolerances from the data. Each comparison uses the measurement's own reported noise: the
  printed replicate SD, else the reading error combined with the instrument's stated repeatability. A model-discrepancy
  term enters only if it was estimated on the same observable and laboratory class, on rows the case lists as spent,
  and only as a per-point variance (never as a bound used as a standard deviation). Across laboratories, the
  between-laboratory spread (H0), on absolute levels only. No term estimated on another observable is transferred."
- *M1.1 (decision statistic):* "A contrast is decided by ONE pre-declared statistic, the pooled lppd difference over its
  decisive units, with the nuisance coordinates conditioned on the case's spent rows. Its threshold is the 95th
  percentile of that statistic under planted-truth data from the incumbent (a parametric bootstrap that redraws every
  input, the conditioning rows included), so false credit is 0.05 by construction. Required-agreement forms and
  adequacy are reported beside the verdict and veto only where a form is a declared transfer of the decisive
  mechanism. The power check uses the scorer's own density estimator; an unresolved density is re-estimated, never read
  as a tie."
- *Expected power:* this record bounds it; it does not establish it. At own noise with the pooled statistic: M1b B1 /
  B2 / B13 vs B0 and B13 vs B3 ~0.97 - 1.00, B3 vs B0 ~0.38; M3 route ~0.87 (R1) / ~0.98 (R2) with the Holm variant at
  0 false credit; M3 transfer ~0.8 at 0.12 false credit here (must fall once the threshold is calibrated); M3 decay law
  ~0. A calibrated threshold lowers these where false credit was above nominal, and a non-zero estimated discrepancy
  lowers them further. Whatever the confirmed power is, it is accepted, "undecidable on disk" included.

**4. Limits of this record.** M1b scores the truth's own draw inside the ensembles (no refit: power overstated); tie
bands are fixed at the truth's expected data (M3's check against the prereg's per-set bands: route 0.89, transfer 0.68
of readings agree, so M3 transfer's numbers are the least reliable); truths conditioned on spent rows at own noise have
small dev-weight ESS where the real rows fit poorly (M1b B13 5.3, B1 11; M3 transfer / decay law ~1: the truth is the
posterior mode, the argon CF candidate's misfit); R = 200 per truth; tie replicates 2 000 (M1b) and 200 (M3). The
reproduction (baseline judge, omega world, the power record's own truth draws) gives 0 power in all five M1b pairs, as
recorded.
"""

PRINCIPLE = """\
**What this is not.** A diagnostic: no prereg, rule or `src/kinetic_core/` file changes and nothing ships. Any revision
it suggests must be frozen BY PRINCIPLE (plan sec. 7 item 5: "each comparison uses the measurement's own reported
noise; across laboratories, the between-laboratory spread (H0)"), confirmed independently (a reviewer who re-derives
the judge's calibration without this record's conclusion) BEFORE it is applied to any case, and the amended judge's
power must then be accepted whatever it is, including "undecidable on disk"."""


def render_revised(rec: Mapping[str, Any]) -> List[str]:
    """The confirmation run's section of the record."""
    cases = [(k, rec["cases"][k]) for k in ("m1b_revised", "m3_revised") if k in rec["cases"]]
    if not cases:
        return []
    L_ = ["## Confirmation run: the REVISED item 5 / M1.1, scored as worded", "",
          "*Coordinator's brief, 2026-11-06, after an independent review (\"adopt with fixes\"). The rule scored, "
          "verbatim:*", "", "> " + REVISED_RULE_TEXT, "",
          "Readings needed to compute it are declared in the script (block \"the REVISED rule\"). Truths come from "
          "the null rival's spent-conditioned posterior OUTSIDE the scoring ensemble. Thresholds are per synthetic "
          "set (the bootstrap pool reweighted by that set's spent rows). A set with any unresolved density (PSIS "
          "k-hat >= 0.7 or ESS < 100, on any rival, form, numerator or denominator) is **not computable** under the "
          "rule. The \"estimated\" columns are a labelled DIAGNOSTIC: the same rule with every importance estimate "
          "read as estimated. They are NOT the rule and credit nothing.", ""]
    for key, c in cases:
        L_ += [f"### {c['case']}: R = {c['R']} sets per truth, bootstrap pool {c['pool']} per null rival, "
               f"{c['seconds'] / 3600:.2f} h on one core", ""]
        L_ += ["| pair | world | rule: credit a power / false | rule: credit b power / false | not computable (of "
               "2R) | estimated (diagnostic): credit a power / false [Wilson 95 %] | credit b power / false | "
               "diagnostic false credit <= 0.05 (exact test, both sides) | threshold MC SE, nats, rule-computable "
               "sets (ESS) |", "|---|---|---|---|---|---|---|---|---|"]
        for pair, v in c["pairs"].items():
            for world in REV_WORLDS:
                w = v["worlds"][world]
                ra, rb = w["credit_a"], w["credit_b"]
                ea, eb = w["estimated_diagnostic"]["credit_a"], w["estimated_diagnostic"]["credit_b"]
                se = w["threshold_mc_se"]["median"]
                L_.append(f"| {pair} | {world} | {ra['power']:.3f} / {ra['false_credit']:.3f} | {rb['power']:.3f} / "
                          f"{rb['false_credit']:.3f} | {w['not_computable']} | {ea['power']:.3f} / "
                          f"{ea['false_credit']:.3f} [{ea['false_wilson95'][0]:.3f}, {ea['false_wilson95'][1]:.3f}] | "
                          f"{eb['power']:.3f} / {eb['false_credit']:.3f} | "
                          f"{ea['p_false_above_nominal'] >= NOMINAL and eb['p_false_above_nominal'] >= NOMINAL} | "
                          f"{'--' if se is None else f'{se:.3g}'} "
                          f"({'--' if w.get('threshold_ess_median') is None else round(w['threshold_ess_median'])}) |")
        L_ += ["", "Leave-one-unit-out and unit contributions (estimated diagnostic, reproducibility world): for each "
               "unit, among the sets where T crossed its threshold, how often the crossing rests on that unit, and of "
               "those how often the unit alone crosses its own null (item 4 for that unit alone); median T_u per "
               "truth.", ""]
        for pair, v in c["pairs"].items():
            w = v["worlds"]["reproducibility"]["estimated_diagnostic"]
            loo = "; ".join(f"{k}: {x['rests_on_it']} / {x['examined']} rest on it, {x['of_which_it_alone_credits']} "
                            f"alone" for k, x in w["leave_one_unit_out"].items() if x["examined"])
            con = "; ".join(f"truth {k.split('_truth_')[1]}: " + ", ".join(f"{u} {x:+.2f}" for u, x in m.items())
                            for k, m in w.items() if k.startswith("contrib_median_truth_") and m)
            L_.append(f"- **{pair}**: {loo or 'no threshold crossings'}. Median contributions: {con}.")
            res = v.get("resolution") or {}
            pool = v.get("pool", {})
            L_.append(f"  - bootstrap pools resolved: " + ", ".join(f"{n} {x['resolved_fraction']:.3f}"
                                                                   for n, x in pool.items())
                      + f"; compute {v['seconds']:.0f} s (pools {v['pool_seconds']:.0f} s)")
        L_ += ["", f"Forced inputs: {c['forced_inputs']}", ""]
    L_ += ["### Reading of the confirmation run", "", REVISED_READING, ""]
    return L_


def render(rec: Mapping[str, Any]) -> str:
    L_ = ["# Judge calibration: why the decision procedure cannot credit a true rival", "",
          f"*Generated by `{GENERATED_BY}` (`--merge`). Planted-truth ablations of the shared judge (M1.1's measurement "
          "model, `src/kinetic_core/predictive_densities.py`, and the sec. 7 / prereg rules) on M1b part (i) (PR #94's "
          "screen ensembles) and M3 phase 1 (the v2.1 power check's draw sets). Synthetic data only, plus spent / "
          "development rows to choose a realistic truth; no decisive or blind data read.*", "", PRINCIPLE, "",
          "## Components ablated (one at a time, then all; the full 2^5 factorial is in the JSON)", ""]
    for k in COMPONENTS:
        L_.append(f"- **{COMPONENT_LABEL[k]}**")
    L_ += ["", "Truth = the rival's central draw conditioned on spent rows (own noise). World **own**: every unit "
           "synthetic at its own reported noise. World **omega**: the generator also adds H0's omega (the original "
           "power checks' generator: data with an H0-sized model discrepancy), the false-credit check under "
           "misspecification. Power = P(credit a | truth a); false credit = P(credit a | truth b). Cells: rate "
           "[Wilson 95 %] | false credit [Wilson 95 %] | one-sided Fisher p (power > false credit).", ""]
    m1b = rec["cases"].get("m1b")
    if m1b:
        L_ += ["## M1b part (i): a rival vs B0 (or B3), R = %d per truth, tie bands from %d replicates" %
               (m1b["R"], m1b["tie_R"]), ""]
        L_ += ["### Truths (central development-conditioned draws; C-Zhou-Cys at its own noise)", "",
               "| rival | draw | ESS of the dev weights | coordinates |", "|---|---|---|---|"]
        for rv, t in m1b["truths"].items():
            co = ", ".join(f"{k} {float(x):.3g}" for k, x in t["coordinates"].items())
            L_.append(f"| {rv} | {t['draw']} / {t['n_draws']} | {t['ess_dev_weights']} | {co} |")
        L_ += ["", "### Reproduction (omega world, the power record's own truth draws, R = %d)" % m1b["repro_R"], "",
               "| pair | recorded power | baseline here | baseline, unresolved read as estimated | all here |",
               "|---|---|---|---|---|"]
        for pair, v in m1b["pairs"].items():
            r = v["reproduction"]
            L_.append(f"| {pair} | {r['recorded_power_rate']:.3f} | {r['baseline']['power']:.3f} | "
                      f"{r['baseline_unresolved_read_as_estimated']['power']:.3f} | {r['all']['power']:.3f} "
                      f"(false {r['all']['false_credit']:.3f}) |")
        for world in ("own", "omega"):
            L_ += ["", f"### World {world}: power | false credit", "",
                   "| configuration | " + " | ".join(m1b["pairs"]) + " |", "|---|" + "---|" * len(m1b["pairs"])]
            rows = SINGLE + [leave_one_in(c) for c in COMPONENTS] + ["all|pooled_lppd"]
            for cfg in rows:
                label = row_label(cfg)
                cells = []
                for pair, v in m1b["pairs"].items():
                    t = v["worlds"][world]["table"][cfg]
                    cells.append(f"{t['power']:.3f} / {t['false_credit']:.3f}")
                L_.append(f"| {label} | " + " | ".join(cells) + " |")
        L_ += ["", "### Detail, world own (rate [Wilson 95 %] | false credit [Wilson 95 %] | Fisher p)", ""]
        for pair, v in m1b["pairs"].items():
            L_ += [f"**{pair}**", "", "| configuration | power | false credit | Fisher p |", "|---|---|---|---|"]
            for cfg in SINGLE + ["all|pooled_lppd"]:
                L_.append(f"| {cfg} | " + _row(v["worlds"]["own"]["table"][cfg]) + " |")
            L_.append("")
            for tname, dg in (("baseline", v["worlds"]["own"]["unit_diagnostics_baseline"]),
                              ("all", v["worlds"]["own"]["unit_diagnostics_all"])):
                for truth, units in dg.items():
                    L_.append(f"- {tname}, truth {truth}, primary form: " + "; ".join(
                        f"{u} d {x['d_median']:.2f} vs band {x['band']:.2f}, resolved {x['resolved_rate']:.2f}, "
                        f"raw a/b/inadequate {x['reads_a_raw']:.2f}/{x['reads_b_raw']:.2f}/"
                        f"{x['both_inadequate_raw']:.2f}" for u, x in units.items()))
            L_.append("")
    m3 = rec["cases"].get("m3")
    if m3:
        L_ += ["## M3 phase 1: the three pairs, world own R = %d, world omega R = %d (primary-form configurations "
               "only), tie bands from %d redraws at the truth's expected data" % (m3["R_own"], m3["R_omega"],
                                                                                 m3["tie_R"]), ""]
        for cname, v in m3["pairs"].items():
            w = v["worlds"]["own"]
            a_name, b_name = w["pair"]
            L_ += [f"### {cname} ({a_name} vs {b_name})", "",
                   "truths: " + "; ".join(f"{t} draw {x['draw']} / {x['n_draws']}, ESS {x['ess_dev_weights']}"
                                          for t, x in w["truths"].items()), "",
                   f"| configuration | credit {a_name}: power | false | credit {b_name}: power | false | omega world: "
                   f"credit {a_name} power / false | credit {b_name} power / false |", "|---|---|---|---|---|---|---|"]
            om = v["worlds"]["omega"]["table"]
            for cfg in SINGLE + [leave_one_in(c) for c in COMPONENTS]:
                ta, tb = w["table"].get(f"{cfg}|credit_{a_name}"), w["table"].get(f"{cfg}|credit_{b_name}")
                oa, ob = om.get(f"{cfg}|credit_{a_name}"), om.get(f"{cfg}|credit_{b_name}")
                f_ = lambda t: "--" if t is None else f"{t['power']:.3f} / {t['false_credit']:.3f}"  # noqa: E731
                L_.append(f"| {row_label(cfg)} | {ta['power']:.3f} | {ta['false_credit']:.3f} | {tb['power']:.3f} | "
                          f"{tb['false_credit']:.3f} | {f_(oa)} | {f_(ob)} |")
            L_ += ["", "diagnostics (primary form, (F)): " + json.dumps(w["diagnostics_primary_form"]), ""]
        if m3.get("band_check"):
            L_ += ["Fixed-band check (baseline, primary, own world, truth a): " + "; ".join(
                f"{b['contrast']}: {b['agreement']:.3f} of {b['readings_compared']} readings agree, per-set / fixed "
                f"sd quantiles {b['sd_ratio_per_set_over_fixed_quantiles']}" for b in m3["band_check"]), ""]
    L_ += ["## Attribution and false-credit control (rules declared in the script before the run)", "",
           "| case | pair | " + " | ".join(COMPONENT_LABEL[c] for c in COMPONENTS) + " | all: false credit <= "
           "nominal (own / omega world) |", "|---|---|" + "---|" * len(COMPONENTS) + "---|"]
    for k, f in rec["findings"].items():
        case, pair = k.split("|", 1)
        cells = []
        for c in COMPONENTS:
            e = f["attribution"][c]
            cells.append(("**carries**" if e["carries"] else "no") + (f" (single {e['single_power']:.2f}, p "
                         f"{e['single_p']:.2g}; restored {e['restored_power']:.2f}, p {e['restored_p']:.2g})"
                         if "single_p" in e and "restored_p" in e else ""))
        L_.append(f"| {case} | {pair} | " + " | ".join(cells) + f" | {f['all_false_credit_ok_own']} / "
                  f"{f['all_false_credit_ok_omega_world']} |")
    L_ += ["", "## Reading", "", rec.get("reading", "(see the PR body; the reading is written after the numbers)"),
           ""]
    L_ += render_revised(rec)
    return "\n".join(L_)


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--case", choices=("m1b", "m3"))
    ap.add_argument("--m1b-root", help="a checkout of PR #94's branch (m1b-part-i-build), read-only")
    ap.add_argument("--m1b-reps", type=int, default=M1B_R)
    ap.add_argument("--m1b-pairs", default="", help="comma list like B1_vs_B0 (default: every pair)")
    ap.add_argument("--m3-cache", help="the M3 v2.1 power check's draw-set cache (pickles written by its generator)")
    ap.add_argument("--m3-reps", type=int, default=200)
    ap.add_argument("--m3-reps-omega", type=int, default=200)
    ap.add_argument("--m3-contrasts", default="route,transfer,decay_law")
    ap.add_argument("--m3-band-check", type=int, default=0, help="sets for the fixed-band check (0: skip)")
    ap.add_argument("--part-out", help="where a case writes its part (scratch; never under data/)")
    ap.add_argument("--merge", nargs="+", help="merge parts into results/validation/judge_calibration.{json,md}")
    ap.add_argument("--revised", action="store_true", help="score the REVISED item 5 / M1.1 (the confirmation run)")
    ap.add_argument("--out-dir", help="--merge: write the record here instead (a dry run; default results/validation)")
    args = ap.parse_args(argv)
    if args.merge:
        return merge(args.merge, args.out_dir)
    mem = memory_free_pct()
    t0 = time.time()
    if args.revised:
        part = run_revised_m1b(args) if args.case == "m1b" else run_revised_m3(args)
    else:
        part = run_m1b(args) if args.case == "m1b" else run_m3(args)
    part["memory_free_pct_at_launch"] = mem
    part["wall_s"] = round(time.time() - t0, 1)
    Path(args.part_out).write_text(json.dumps(_jsonable(part)) + "\n")
    print(f"wrote part {args.part_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
