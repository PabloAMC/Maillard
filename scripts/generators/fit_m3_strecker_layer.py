#!/usr/bin/env python
"""
fit_m3_strecker_layer.py -- M3 phase 1: the forced-input Strecker layer's checks, mechanism screen, power checks, and
(built, tested, NOT run before the freeze) its decisive fits and scoring.

Design: results/validation/kinetic_core_m3_strecker_layer_prereg.md (v2 DRAFT). The layer: src/m3_strecker_layer.py.

    python scripts/generators/fit_m3_strecker_layer.py --checks   # V-acid, V-time, V1/V2/V4/V5, V-density, PSIS pilot,
                                                                  # the blind-order (V-blind) stub
    python scripts/generators/fit_m3_strecker_layer.py --screen   # sec. 8: >= 600 Sobol draws per rival, per form
    python scripts/generators/fit_m3_strecker_layer.py --power    # sec. 12.1 (U3, synthetic) and 12.2 (three pairs)
    python scripts/generators/fit_m3_strecker_layer.py --fit      # REFUSED until the prereg is frozen (sec. 20)
    python scripts/generators/fit_m3_strecker_layer.py --score    # REFUSED until frozen AND V-blind (i)-(iii) pass

Writes results/validation/kinetic_core_m3_phase1_{checks,screen,power}.{json,md}.

Rules enforced here:
- memory guard: a compute mode refuses to start unless free memory >= 40 % (``memory_pressure``) [--no-memory-guard
  only for tests];
- compute cap: --cap-hours (default 24 for --fit / --score, 6 otherwise); a capped part reads "compute-capped";
- failed runs: every draw set counts feasible draws whose predictions are not finite; a part that used a set with a
  failure is REFUSED in the record;
- run-once: --fit and --score refuse when their record exists;
- freeze checks: --fit / --score need the prereg's sec. 20 freeze commit, an ancestor of HEAD, with the design text
  above sec. 20 unchanged since it;
- V-blind (sec. 7): --score needs (i) the freeze commit an ancestor of the commit that ADDED the keyed file, (ii) that
  commit an ancestor of HEAD, (iii) the keyed file's hash now equal to its hash at that commit (and no later commit
  changed it), before any real-data trapped-pot number exists (iv).
Nothing here reads hofmann2000b Figs. 2 - 3 or the keyed file outside ``v_blind`` and the post-V-blind scorer.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths, provenance  # noqa: E402
from src import m3_strecker_layer as L  # noqa: E402
from src.kinetic_core import fit_rows  # noqa: E402
from src.kinetic_core import predictive_densities as pd  # noqa: E402

PREREG = "results/validation/kinetic_core_m3_strecker_layer_prereg.md"
KEYED = f"data/lit/fit_rows/{L.KEYED_TRAPPED_TABLE}.yml"
STEM = "kinetic_core_m3_phase1"
MEMORY_MIN_FREE_PCT = 40
Z975 = 1.959963984540054
FED_ORDER = ("U4", "U1e", "U1l", "U3")
GLC_ORDER = ("U5a", "U5b", "UARP")
CONTRASTS = {
    "route": {"a": "R1", "b": "R2", "claim": ("U3", "U1e", "U4"), "primary": "U3", "guards": ("U1l", "GT3")},
    "transfer": {"a": "R1", "b": "R3", "claim": ("U5a", "U5b", "UARP"), "primary": "block",
                 "guards": ("U1e", "U1l", "U4", "GT3")},
    "decay_law": {"a": "R4", "b": "R1", "claim": ("U3", "U5a", "U5b", "UARP"), "primary": "block",
                  "guards": ("U1e", "U1l", "U4", "GT3")},
}
POWER_FORMS = {"primary": L.PRIMARY, "o_ar_0.1": L.Form(o_ar=0.1), "cu_linear": L.Form(cu_law="linear"),
               "headspace_20mL": L.Form(headspace_ml=20.0)}


# =================================================================================================== guards

class Refusal(RuntimeError):
    """A rule of the prereg or of the run (memory, freeze, blinding, run-once, cap) refuses the requested part."""


def memory_free_pct() -> Optional[int]:
    try:
        out = subprocess.run(["memory_pressure"], capture_output=True, text=True, timeout=60).stdout
    except Exception:
        return None
    m = re.search(r"System-wide memory free percentage:\s*(\d+)%", out)
    return int(m.group(1)) if m else None


def swap_used_mb() -> Optional[float]:
    try:
        out = subprocess.run(["sysctl", "vm.swapusage"], capture_output=True, text=True, timeout=30).stdout
    except Exception:
        return None
    m = re.search(r"used = ([\d.]+)M", out)
    return float(m.group(1)) if m else None


def memory_guard(enabled: bool) -> Dict[str, Any]:
    free, swap = memory_free_pct(), swap_used_mb()
    rec = {"free_pct_at_launch": free, "swap_used_mb_at_launch": swap, "threshold_pct": MEMORY_MIN_FREE_PCT,
           "guard": enabled}
    if enabled and (free is None or free < MEMORY_MIN_FREE_PCT):
        raise Refusal(f"memory guard: free {free} % < {MEMORY_MIN_FREE_PCT} %; wait (check at most every 10 min)")
    return rec


@dataclass
class Budget:
    """The compute cap (wall clock). ``check`` raises once the cap is passed; parts record "compute-capped"."""
    hours: float
    t0: float = field(default_factory=time.time)

    def elapsed_h(self) -> float:
        return (time.time() - self.t0) / 3600.0

    def check(self, what: str) -> None:
        if self.elapsed_h() > self.hours:
            raise Refusal(f"compute-capped at {self.hours} h wall during {what}")


def _git(*args: str, cwd: Path = ROOT, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=str(cwd), capture_output=True, text=True, check=check)


def is_ancestor(a: str, b: str, cwd: Path = ROOT) -> bool:
    return _git("merge-base", "--is-ancestor", a, b, cwd=cwd, check=False).returncode == 0


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def freeze_commit(prereg_text: str) -> Optional[str]:
    """The freeze commit entered in the prereg's sec. 20 ("Freeze commit: <40 hex>"), or None."""
    part = prereg_text.split("## 20.", 1)
    if len(part) < 2:
        return None
    m = re.search(r"[Ff]reeze commit[^0-9a-f]*([0-9a-f]{40})", part[1])
    return m.group(1) if m else None


def freeze_check(cwd: Path = ROOT, prereg: str = PREREG) -> Dict[str, Any]:
    """The prereg is frozen: sec. 20 names a freeze commit that is an ancestor of HEAD, and the design text above
    sec. 20 is unchanged since it. Raises Refusal otherwise."""
    text = (cwd / prereg).read_text()
    sha = freeze_commit(text)
    if sha is None:
        raise Refusal("prereg not frozen: sec. 20 names no freeze commit (v2 is a DRAFT)")
    if not is_ancestor(sha, "HEAD", cwd):
        raise Refusal(f"freeze commit {sha} is not an ancestor of HEAD")
    then = _git("show", f"{sha}:{prereg}", cwd=cwd).stdout
    if then.split("## 20.", 1)[0] != text.split("## 20.", 1)[0]:
        raise Refusal("the design text above sec. 20 changed after the freeze commit")
    return {"freeze_commit": sha, "design_unchanged": True}


def v_blind(freeze_sha: str, keyed: str = KEYED, cwd: Path = ROOT) -> Dict[str, Any]:
    """Prereg sec. 7, V-blind (i)-(iii). (i) the frozen design commit is an ANCESTOR of the commit that added the
    keyed file; (ii) that commit is an ancestor of HEAD; (iii) the keyed file's hash at scoring equals its hash at
    that commit, and no later commit changed it. (iv) is the caller's: nothing real is predicted before this passes.
    Raises Refusal with the failing condition."""
    added = _git("log", "--diff-filter=A", "--format=%H", "--", keyed, cwd=cwd, check=False).stdout.split()
    if not added:
        raise Refusal(f"V-blind: {keyed} is not committed (keying comes after the freeze)")
    keyed_commit = added[-1]                       # the first commit that added it (log is newest first)
    if not is_ancestor(freeze_sha, keyed_commit, cwd):
        raise Refusal(f"V-blind (i): freeze commit {freeze_sha} is not an ancestor of the keyed commit {keyed_commit}")
    if not is_ancestor(keyed_commit, "HEAD", cwd):
        raise Refusal(f"V-blind (ii): the keyed commit {keyed_commit} is not an ancestor of HEAD")
    at_commit = subprocess.run(["git", "show", f"{keyed_commit}:{keyed}"], cwd=str(cwd), capture_output=True,
                               check=True).stdout
    path = cwd / keyed
    now = path.read_bytes() if path.exists() else b""
    changed = _git("log", "--format=%H", f"{keyed_commit}..HEAD", "--", keyed, cwd=cwd, check=False).stdout.split()
    if sha256_bytes(now) != sha256_bytes(at_commit) or changed:
        raise Refusal("V-blind (iii): the keyed file's hash changed since its commit")
    return {"freeze_commit": freeze_sha, "keyed_commit": keyed_commit, "keyed_sha256": sha256_bytes(now),
            "conditions": ["i", "ii", "iii"], "passed": True}


def run_once(path: Path) -> None:
    if path.exists():
        raise Refusal(f"run-once: {path.name} exists; the decisive part runs once")


# ============================================================================================ draw sets

@dataclass
class DrawSet:
    """One rival's FEASIBLE draws under one layer form, with every unit's predicted points (sec. 12.2's precomputed
    draw sets). Infeasible draws have weight 0 under the calibration (sec. 4.4); their share is recorded."""
    rival: str
    form: L.Form
    seed: int
    n_drawn: int
    u: np.ndarray
    points: Dict[str, np.ndarray]
    uarp_c: Optional[np.ndarray]
    cv: Optional[np.ndarray]
    f_iso: Optional[np.ndarray]
    k_glc: Optional[np.ndarray]           # (n, 501) in-pot exit rate on a 1-min grid (truth sets only)
    theta: Dict[str, np.ndarray]
    diag: Dict[str, np.ndarray]
    n_failed: int
    wall_s: float
    reasons: Dict[str, int]

    @property
    def n(self) -> int:
        return self.u.shape[0]


UNIT_NAMES = ("U4", "U1e", "U1l", "U3", "U5a", "U5b", "U5a_PA_shape", "U5b_PA_shape", "GT3")


def cached_draw_set(cache: Optional[Path], rival: str, form: L.Form, data: L.Data, target_feasible: int, seed: int,
                    **kw) -> DrawSet:
    """``make_draw_set`` with an on-disk cache (pickle, written by this generator only) keyed by everything that
    determines the set, the layer's hash and make_draw_set's source, so a restarted or second process reuses it."""
    if cache is None:
        return make_draw_set(rival, form, data, target_feasible, seed, **kw)
    import inspect
    import pickle
    h = hashlib.sha256()
    h.update((ROOT / "src" / "m3_strecker_layer.py").read_bytes())
    h.update(inspect.getsource(make_draw_set).encode())
    h.update(repr((rival, form, target_feasible, seed, sorted((k, v) for k, v in kw.items() if k != "budget")))
             .encode())
    path = cache / f"drawset_{rival}_{h.hexdigest()[:16]}.pkl"
    if path.exists():
        with path.open("rb") as fh:
            return pickle.load(fh)
    ds = make_draw_set(rival, form, data, target_feasible, seed, **kw)
    cache.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    with tmp.open("wb") as fh:
        pickle.dump(ds, fh)
    tmp.rename(path)
    return ds


def make_draw_set(rival: str, form: L.Form, data: L.Data, target_feasible: int, seed: int, batch: int = 256,
                  max_draws: int = 12288, keep_k: bool = False, budget: Optional[Budget] = None,
                  trapped: bool = True, min_draws: int = 0) -> DrawSet:
    """Sobol draws in batches until ``target_feasible`` feasible draws (or ``max_draws``); at least ``min_draws``."""
    d = len(L.coords(rival, form))
    total = max(target_feasible, min_draws)
    u_all = L.sobol_unit(max_draws, d, seed)
    t0 = time.time()
    keep: Dict[str, List[np.ndarray]] = {}
    us, ths, diags, cs, cvs, fis, ks = [], [], [], [], [], [], []
    n_drawn, n_feas, n_failed = 0, 0, 0
    reasons: Dict[str, int] = {}
    while n_drawn < max_draws and (n_feas < target_feasible or n_drawn < min_draws):
        if budget is not None:
            budget.check(f"draw set {rival} {form.label()}")
        u = u_all[n_drawn:n_drawn + batch]
        n_drawn += u.shape[0]
        ev = L.evaluate(u, rival, form, data, trapped=trapped and rival != "R0")
        for k, v in ev.reasons.items():
            reasons[k] = reasons.get(k, 0) + int(np.sum(v))
        f = ev.feasible & ~ev.diag.get("failed_mask", np.zeros(len(u), bool))
        n_failed += ev.failed
        n_feas += int(f.sum())
        us.append(u[f])
        ths.append({k: np.asarray(v)[f] for k, v in ev.theta.items()})
        for unit in UNIT_NAMES:
            try:
                pts = L.points_for(unit, ev)
            except KeyError:
                continue
            keep.setdefault(unit, []).append(pts[f])
        if rival != "R0":
            cs.append(ev.pred["UARP_c"][f])
            cvs.append(ev.theta["cv_arp"][f])
            fis.append(ev.theta["f_iso"][f] if "f_iso" in ev.theta else np.zeros(int(f.sum())))
            diags.append({k: np.asarray(v)[f] for k, v in ev.diag.items() if np.ndim(v) == 1 and len(v) == len(u)})
            if keep_k:
                ks.append(ev.diag["k_eff_1min"][f])
    cat = lambda xs: np.concatenate(xs) if xs else None  # noqa: E731
    th = {k: np.concatenate([t[k] for t in ths]) for k in ths[0]} if ths else {}
    dg = {k: np.concatenate([t[k] for t in diags]) for k in diags[0]} if diags else {}
    return DrawSet(rival, form, seed, n_drawn, np.concatenate(us), {k: np.concatenate(v) for k, v in keep.items()},
                   cat(cs), cat(cvs), cat(fis), cat(ks), th, dg, n_failed, round(time.time() - t0, 1), reasons)


# ===================================================================================== log likelihoods

@dataclass
class Obs:
    """A data set for scoring: the units (real or synthetic) and U-ARPglc's measured ln ratios."""
    units: Dict[str, L.Unit]
    uarp_y: np.ndarray


def real_obs(data: L.Data) -> Obs:
    return Obs(L.build_units(data), L.arp_interval_measured(data))


def loglik(ds: DrawSet, unit: str, obs: Obs, form: L.Form, data: L.Data) -> np.ndarray:
    if unit == "UARP":
        meas = L.arp_interval_measured(data, ds.f_iso) if form.isomers else obs.uarp_y
        cv = ds.cv if form.cv_arp is None else np.full(ds.n, form.cv_arp)
        return L.uarp_loglik(data, ds.uarp_c, cv, form, measured=meas)
    u = obs.units[unit]
    return L.unit_loglik(u, ds.points[unit], form)


def lse(x: np.ndarray) -> float:
    x = np.asarray(x, float)
    fin = x[np.isfinite(x)]
    if fin.size == 0:
        return -math.inf
    m = float(np.max(fin))
    return m + math.log(float(np.sum(np.exp(fin - m))))


def conditional_lppd(ll_cond: np.ndarray, ll_unit: np.ndarray, psis_check: bool = True) -> Dict[str, Any]:
    """lppd(u | cond) by reweighting the draw set (proposal = the band prior conditioned on the calibrations):
    log sum exp(cond + u) - log sum exp(cond). With ``psis_check`` the numerator's PSIS k-hat and ESS are computed
    (predictive_densities.importance_lppd); "resolved" iff the prereg's acceptance rule holds."""
    num, den = ll_cond + ll_unit, ll_cond
    out = {"lppd": lse(num) - lse(den)}
    if psis_check:
        z = np.zeros_like(num)
        e_num = pd.importance_lppd(np.where(np.isfinite(num), num, -np.inf), z)
        e_den = pd.importance_lppd(np.where(np.isfinite(den), den, -np.inf), z)
        out.update({"k_hat": e_num.k_hat, "ess": e_num.ess, "resolved": bool(e_num.accepted and e_den.accepted),
                    "lppd_psis": e_num.lppd - e_den.lppd})
    return out


def predictive_moments(ds: DrawSet, unit: L.Unit, pts: np.ndarray, w_log: np.ndarray, form: L.Form
                       ) -> Tuple[np.ndarray, np.ndarray]:
    """The conditioned predictive mean and covariance of a unit's values (Sigma_u + the predictive spread)."""
    w = np.exp(w_log - np.max(w_log[np.isfinite(w_log)]))
    w = np.where(np.isfinite(w), w, 0.0)
    w = w / w.sum()
    f = pts @ unit.design.T
    mu = w @ f
    dev = f - mu
    cov = L.unit_covariance(unit, form) + (dev * w[:, None]).T @ dev
    return mu, cov


def adequacy_p(unit: L.Unit, mu: np.ndarray, cov: np.ndarray) -> float:
    """Posterior-predictive p of the Mahalanobis discrepancy; censored values are left out (declared)."""
    keep = ~unit.censored if unit.censored is not None else np.ones(len(unit.y), bool)
    if not keep.any():
        return 1.0
    return pd.posterior_predictive_p(unit.y[keep], mu[keep], cov[np.ix_(keep, keep)])


# ===================================================================================== synthetic data

def synth_unit(unit: L.Unit, f_points: np.ndarray, rng: np.random.Generator, form: L.Form = L.PRIMARY) -> L.Unit:
    """A synthetic copy of a unit from a truth's points, with noise from the unit's FULL Sigma_u (omega included)."""
    if unit.name == "U3":
        cov_p = L.unit_covariance(L.u3_unit(10.0 ** f_points), form)
        y_log = f_points + rng.multivariate_normal(np.zeros(len(f_points)), cov_p)
        return L.u3_unit(10.0 ** y_log)
    y = unit.design @ f_points
    cov = L.unit_covariance(unit, form)
    return replace(unit, y=y + rng.multivariate_normal(np.zeros(len(y)), cov))


def synth_uarp(data: L.Data, k_1min: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """U-ARPglc's generative model (sec. 12.2, M8): the formation F(t) >= 0 of minimum norm that reproduces the measured
    interpolant under the truth's k(t), clipped at 0; ARP' = F - k ARP; integers (umol) with the analyser CV drawn
    from its band. Returns the ln ratios per interval."""
    t = np.arange(0.0, 501.0, 1.0)
    ts, vs = data.forced_inputs["PHEARP"]
    a = np.interp(t, ts, vs)
    da = np.gradient(a, t)
    f = np.maximum(da + k_1min * a, 0.0)
    arp = np.zeros_like(a)
    for i in range(1, t.size):
        k = k_1min[i - 1]
        ex = math.exp(-k)
        arp[i] = arp[i - 1] * ex + (f[i - 1] / k * (1 - ex) if k > 0 else f[i - 1])
    cv = rng.uniform(*L.CV_ARP_BAND)
    times = sorted(data.arp_glc)
    obs = {}
    for tt in times:
        v = arp[int(tt)] * math.exp(cv * rng.standard_normal())
        obs[tt] = max(round(v * 10.0), 1) / 10.0
    return np.array([[math.log(obs[a0] / obs[b0]) for a0, b0 in L.ARP_INTERVALS]])


def weights_from(ds: DrawSet, units: Sequence[str], obs: Obs, form: L.Form, data: L.Data) -> np.ndarray:
    ll = np.zeros(ds.n)
    for u in units:
        ll = ll + loglik(ds, u, obs, form, data)
    return ll


def sample_truths(ds: DrawSet, log_w: np.ndarray, r: int, rng: np.random.Generator) -> np.ndarray:
    w = np.exp(log_w - np.max(log_w[np.isfinite(log_w)]))
    w = np.where(np.isfinite(w), w, 0.0)
    return rng.choice(ds.n, size=r, replace=True, p=w / w.sum())


# ============================================================================================= decisions

def _representable(ds: DrawSet, cond_ll: np.ndarray, unit_ll: Optional[np.ndarray] = None) -> bool:
    """A rival represents the conditioning (and unit) data iff it has a feasible draw with finite weight."""
    if ds.n == 0:
        return False
    w = cond_ll if unit_ll is None else cond_ll + unit_ll
    return bool(np.any(np.isfinite(w)))


def _one_sided(a_ok: bool, b_ok: bool) -> Dict[str, Any]:
    """A rival whose calibration is infeasible under a form (no feasible draw: sigma-stoich / sigma-O2 there) has
    lppd = -inf; the other side is read as favoured (flagged; its adequacy is not tested: optimistic)."""
    side = "a" if a_ok and not b_ok else "b" if b_ok and not a_ok else "tie"
    d = math.inf if side == "a" else -math.inf if side == "b" else math.nan
    return {"d": d, "band": math.nan, "sd_rep": math.nan, "reading": side, "reading_raw": side, "resolved": True,
            "by_infeasibility": True, "p_a": math.nan, "p_b": math.nan, "k_hat_a": None, "k_hat_b": None,
            "ess_a": None, "ess_b": None}


def unit_verdict(a: DrawSet, b: DrawSet, unit: str, cond: Sequence[str], obs: Obs, data: L.Data, form: L.Form,
                 rng: np.random.Generator, n_tie: int, psis_rule: bool = True) -> Dict[str, Any]:
    """One unit's reading between rivals a and b: d = lppd_a - lppd_b (prequential: conditioned on ``cond``), the tie
    band from measurement-noise redraws (sec. 9.2), adequacy (sec. 9.4), "density not resolved" -> tie."""
    if a.n == 0 or b.n == 0:
        return _one_sided(a.n > 0, b.n > 0)
    ca = weights_from(a, cond, obs, form, data)
    cb = weights_from(b, cond, obs, form, data)
    la, lb = loglik(a, unit, obs, form, data), loglik(b, unit, obs, form, data)
    ra, rb = _representable(a, ca, la), _representable(b, cb, lb)
    if not (ra and rb):
        return _one_sided(ra, rb)
    ea, eb = conditional_lppd(ca, la, psis_rule), conditional_lppd(cb, lb, psis_rule)
    d = ea["lppd"] - eb["lppd"]
    reps = []
    for _ in range(n_tie):
        o2 = redraw(obs, unit, data, rng, form)
        reps.append(lse(ca + loglik(a, unit, o2, form, data)) - lse(ca)
                    - lse(cb + loglik(b, unit, o2, form, data)) + lse(cb))
    reps = np.array([x for x in reps if np.isfinite(x)])
    sd_rep = float(np.std(reps, ddof=1)) if reps.size > 2 else float("nan")
    band = pd.Z_90 * sd_rep
    resolved = (ea.get("resolved", True) and eb.get("resolved", True)) if psis_rule else True
    if unit == "UARP":
        reading = "tie" if not (np.isfinite(d) and abs(d) > band) else ("a" if d > 0 else "b")
        p_a = p_b = float("nan")
    else:
        u = obs.units[unit]
        mu_a, cov_a = predictive_moments(a, u, a.points[unit], ca, form)
        mu_b, cov_b = predictive_moments(b, u, b.points[unit], cb, form)
        p_a, p_b = adequacy_p(u, mu_a, cov_a), adequacy_p(u, mu_b, cov_b)
        reading = pd.unit_reading(d, band, p_a, p_b, mu_a, mu_b, L.OMEGA_DEX * form.omega_scale * math.sqrt(u.c))
    reading_raw = reading
    if not resolved:
        reading = "tie"
    return {"d": d, "band": band, "sd_rep": sd_rep, "reading": reading, "reading_raw": reading_raw,
            "resolved": bool(resolved),
            "p_a": p_a, "p_b": p_b, "k_hat_a": ea.get("k_hat"), "k_hat_b": eb.get("k_hat"),
            "ess_a": ea.get("ess"), "ess_b": eb.get("ess")}


def redraw(obs: Obs, unit: str, data: L.Data, rng: np.random.Generator, form: L.Form) -> Obs:
    """A replicate of one unit's data from MEASUREMENT noise only (read error, rounding, s_rep, sigma_run,
    calibration coordinates, the analyser CV at its band midpoint); never omega (sec. 9.2)."""
    if unit == "UARP":
        cv = 0.5 * sum(L.CV_ARP_BAND)
        ts = sorted(data.arp_glc)
        umol = np.array([data.arp_glc[t] * 10.0 for t in ts])
        sd = np.sqrt(cv ** 2 + ((0.5 / math.sqrt(3.0)) / umol) ** 2)
        s = np.array([math.hypot(sd[ts.index(a)], sd[ts.index(b)]) for a, b in L.ARP_INTERVALS])
        return Obs(obs.units, obs.uarp_y + s * rng.standard_normal(s.size))
    u = obs.units[unit]
    cov = L.unit_covariance(u, form, measurement_only=True)
    if u.name == "U3":
        y_lin = 10.0 ** u.y
        y_lin = np.where(u.censored, 0.5 * L.U3_CENSOR, y_lin)
        y_log = np.log10(y_lin) + rng.multivariate_normal(np.zeros(len(u.y)), cov)
        new = L.u3_unit(10.0 ** y_log)
    else:
        new = replace(u, y=u.y + rng.multivariate_normal(np.zeros(len(u.y)), cov))
    return Obs({**obs.units, unit: new}, obs.uarp_y)


def block_verdict(a: DrawSet, b: DrawSet, block: Sequence[str], cond: Sequence[str], obs: Obs, data: L.Data,
                  form: L.Form, rng: np.random.Generator, n_tie: int, psis_rule: bool = True) -> Dict[str, Any]:
    """The summed block's reading (transfer / decay law primary): d = lppd_a(block | cond) - lppd_b(block | cond),
    with its own tie band (all of the block's data redrawn)."""
    if a.n == 0 or b.n == 0:
        return _one_sided(a.n > 0, b.n > 0)
    ca, cb = weights_from(a, cond, obs, form, data), weights_from(b, cond, obs, form, data)
    ba, bb = weights_from(a, block, obs, form, data), weights_from(b, block, obs, form, data)
    ra, rb = _representable(a, ca, ba), _representable(b, cb, bb)
    if not (ra and rb):
        return _one_sided(ra, rb)
    ea, eb = conditional_lppd(ca, ba, psis_rule), conditional_lppd(cb, bb, psis_rule)
    d = ea["lppd"] - eb["lppd"]
    reps = []
    for _ in range(n_tie):
        o2 = obs
        for u in block:
            o2 = redraw(o2, u, data, rng, form)
        reps.append(lse(ca + weights_from(a, block, o2, form, data)) - lse(ca)
                    - lse(cb + weights_from(b, block, o2, form, data)) + lse(cb))
    reps = np.array([x for x in reps if np.isfinite(x)])
    sd_rep = float(np.std(reps, ddof=1)) if reps.size > 2 else float("nan")
    band = pd.Z_90 * sd_rep
    resolved = (ea.get("resolved", True) and eb.get("resolved", True)) if psis_rule else True
    reading = "tie" if (not resolved or not (np.isfinite(d) and abs(d) > band)) else ("a" if d > 0 else "b")
    return {"d": d, "band": band, "reading": reading, "resolved": bool(resolved)}


def contrast_credit(primary: Dict[str, Any], per_unit: Mapping[str, Dict[str, Any]]) -> str:
    """Sec. 10.1's rule on the claim units (Holm's guard step is NOT simulated in the power check: optimistic):
    credit a iff the primary favours a and no claim unit favours b; symmetric for b; else none."""
    if primary["reading"] == "a" and not any(v["reading"] == "b" for v in per_unit.values()):
        return "a"
    if primary["reading"] == "b" and not any(v["reading"] == "a" for v in per_unit.values()):
        return "b"
    return "none"


# ============================================================================================ the checks

def _record(name: str, payload: Dict[str, Any], md: str, out_dir: Path, inputs: Sequence[Path]) -> None:
    stem = out_dir / f"{STEM}_{name}"
    payload = {"artifact": stem.name,
               "provenance": provenance.provenance_block(stem.name, generated_by=f"scripts/generators/{Path(__file__).name}"
                                                         f" --{name}", inputs=list(inputs)),
               "design": f"{PREREG} (v2 DRAFT; blob {prereg_blob()})", **payload}
    out_dir.mkdir(parents=True, exist_ok=True)
    stem.with_suffix(".json").write_text(provenance.dumps_compact(payload, default=_jsonable) + "\n")
    stem.with_suffix(".md").write_text(md + "\n")
    print(f"wrote {stem.with_suffix('.md')}")


def _jsonable(x):
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.bool_):
        return bool(x)
    return str(x)


def prereg_blob() -> str:
    r = _git("hash-object", PREREG, check=False)
    return r.stdout.strip()[:12]


def layer_inputs() -> List[Path]:
    return [ROOT / "src" / "m3_strecker_layer.py", Path(__file__).resolve()] + \
           [fit_rows.FIT_ROWS_DIR / f"{t}.yml" for t in L.TABLES]


def v1_ex049(data: L.Data) -> Dict[str, Any]:
    """V1: route A's calibration and the glucose pot's dicarbonyl part bit-exact against EX-049 at its structure
    limit (the same theta)."""
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))
    import feasibility_strecker_forced_input as fs  # noqa: E402
    rng = np.random.default_rng(5)
    th = fs.draw(64, rng)
    a = fs.calibrate(th, fs.load_data())
    b = L.calibrate_route_a(th, data)
    same = all(np.array_equal(a[key][d], b[key][d]) for key in ("k", "phi") for d in L.DONORS)
    same &= all(np.array_equal(a["omega"][c], b["omega"][c]) for c in ("C23", "C6"))
    same &= np.array_equal(a["t3_total"], b["t3_total"]) and np.array_equal(a["t3_pa"], b["t3_pa"])
    cfg = fs.config("DC")
    g_old = fs.glucose_pot(th, a, fs.interpolate(fs.load_data().forced_inputs), cfg)
    # the new glucose pot at the structure limit: no ARP-made donors, no route B
    u = L.sobol_unit(64, len(L.coords("R1")), 7)
    thn = L.from_unit(u, "R1")
    for k in ("log10_c_mgo", "nu", "log10_k_loss", "off_kappa_go", "off_kappa_g", "off_kappa_tdg", "cal_mgo",
              "cal_go", "cal_tdg", "cal_g", "cal_arp"):
        thn[k] = th[k]
    thn["o2_reflux"] = th["o2"]
    p = L.derive(thn, "R1", L.PRIMARY, data)
    p.sigma_ox = np.zeros(64)
    p.a = (np.zeros(64), np.zeros(64), np.ones(64))
    g_new = L.glucose_pot(p, data)
    same_glc = np.array_equal(g_old["PA_dc"], g_new["PA_dc"]) and np.array_equal(g_old["PAA_dc"], g_new["PAA_dc"])
    return {"route_a_calibration_bit_exact": bool(same), "glucose_pot_dicarbonyl_bit_exact": bool(same_glc),
            "draws": 64}


def v2_balance(data: L.Data, n: int = 32) -> Dict[str, Any]:
    """V2: Phe-moiety, sugar-moiety and O2 balance in every ARP pot, relative error."""
    out = {}
    for rival in ("R1", "R2", "R4"):
        for form in (L.PRIMARY, L.Form(headspace_ml=20.0, nu_o2=0.5)):
            u = L.sobol_unit(n, len(L.coords(rival, form)), 11)
            p = L.derive(L.from_unit(u, rival, form), rival, form, data)
            sig, cu = np.full(n, 0.3), np.full(n, 3.0)
            worst = {"phe": 0.0, "sugar": 0.0, "o2": 0.0}
            for pot in list(L.POT_T2.values()) + [L.POT_FIG1] + list(L.POT_TRAP.values()):
                g = L.solve_arp_pot(p, pot, sig, cu, ledger=True)["ledger"]
                phe = g["arp"] + g["phe"] + g["pa"] + g["pa_lost"] + g["paa"] + g["phe_sink"] + g["phe_bound"]
                worst["phe"] = max(worst["phe"], float(np.max(np.abs(phe - g["arp0"] - g["phe0"]) / g["arp0"])))
                sugar = g["donors"] + g["trapped"] + g["donor_sink"] + g["pa"] + g["pa_lost"] + g["paa"] + \
                    g["sugar_frag"]
                worst["sugar"] = max(worst["sugar"], float(np.max(np.abs(sugar - (g["arp0"] - g["arp"])) / g["arp0"])))
                if g["cap"] is not None:
                    ok = g["o2_clipped"] == 0
                    if ok.any():
                        worst["o2"] = max(worst["o2"], float(np.max(np.abs(g["o2"] + g["o2_used"] - g["cap"])[ok]
                                                                    / g["cap"])))
            out[f"{rival}|{form.label()}"] = worst
    tol = 1e-9
    return {"per_rival_form": out, "tolerance": tol,
            "pass": all(v <= tol for w in out.values() for v in w.values())}


def v_time(data: L.Data, n: int = 256) -> Dict[str, Any]:
    out = {}
    for rival in L.RIVALS:
        u = L.sobol_unit(n, len(L.coords(rival)), 3)
        t = time.time()
        ev = L.evaluate(u, rival, L.PRIMARY, data)
        dt = time.time() - t
        out[rival] = {"ms_per_draw": round(1000 * dt / n, 2), "batch": n, "feasible_share": float(ev.feasible.mean()),
                      "failed": ev.failed, "coordinates": len(L.coords(rival))}
    return out


def v_density_toy(data: L.Data, rng: np.random.Generator) -> Dict[str, Any]:
    """V-density: on a two-coordinate R3 toy (log10 k_d,7 and log10 k_d,7,glc free; every other coordinate at its
    band centre), U-ARPglc's lppd by a fine trapezoid grid, by adaptive Gauss-Hermite and by PSIS from a Student-t
    proposal must agree."""
    cs = L.coords("R3")
    names = [c.name for c in cs]
    i1, i2 = names.index("log10_kd7"), names.index("log10_kd7_glc")
    g = 41
    grid = np.linspace(0.0, 1.0, g)
    uu = np.full((g * g, len(cs)), 0.5)
    m1, m2 = np.meshgrid(grid, grid, indexing="ij")
    uu[:, i1], uu[:, i2] = m1.ravel(), m2.ravel()
    ev = L.evaluate(uu, "R3", L.PRIMARY, data, trapped=False)
    ll = L.uarp_loglik(data, ev.pred["UARP_c"], ev.theta["cv_arp"], L.PRIMARY, measured=L.arp_interval_measured(data))
    ll = np.where(ev.feasible, ll, -np.inf).reshape(g, g)
    lmax = np.max(ll[np.isfinite(ll)])
    w = np.exp(ll - lmax)
    trap = np.trapezoid(np.trapezoid(w, grid, axis=1), grid)
    lppd_grid = math.log(trap) + lmax
    # PSIS with a Student-t proposal at the grid mode, scale from the integrand's width (lesson 2026-11-05)
    wn = w / w.sum()
    mean = np.array([np.sum(wn * m1), np.sum(wn * m2)])
    cov = np.cov(np.vstack([m1.ravel(), m2.ravel()]), aweights=wn.ravel()) * 2.0
    x, lq = pd.multivariate_t_draws(mean, cov, 2048, rng)
    inside = np.all((x >= 0) & (x <= 1), axis=1)
    up = np.full((x.shape[0], len(cs)), 0.5)
    up[:, i1], up[:, i2] = np.clip(x[:, 0], 0, 1), np.clip(x[:, 1], 0, 1)
    ev2 = L.evaluate(up, "R3", L.PRIMARY, data, trapped=False)
    l2 = L.uarp_loglik(data, ev2.pred["UARP_c"], ev2.theta["cv_arp"], L.PRIMARY,
                       measured=L.arp_interval_measured(data))
    li = np.where(inside & ev2.feasible, l2, -np.inf)            # uniform prior on the unit square: log p = 0
    est = pd.importance_lppd(li, lq)

    def gh(pts):
        upp = np.full((pts.shape[0], len(cs)), 0.5)
        upp[:, i1], upp[:, i2] = np.clip(pts[:, 0], 0, 1), np.clip(pts[:, 1], 0, 1)
        e = L.evaluate(upp, "R3", L.PRIMARY, data, trapped=False)
        v = L.uarp_loglik(data, e.pred["UARP_c"], e.theta["cv_arp"], L.PRIMARY, measured=L.arp_interval_measured(data))
        ok = np.all((pts >= 0) & (pts <= 1), axis=1) & e.feasible
        return np.where(ok, v, -1e300)

    q = pd.quadrature_lppd(gh, mean, cov / 2.0, n_nodes=20)
    return {"lppd_trapezoid_41x41": round(lppd_grid, 4), "lppd_psis": round(est.lppd, 4), "k_hat": est.k_hat,
            "ess": est.ess, "mc_se": est.mc_se, "lppd_gauss_hermite_20": round(q.lppd, 4),
            "psis_minus_grid": round(est.lppd - lppd_grid, 4), "gh_minus_grid": round(q.lppd - lppd_grid, 4),
            "pass": bool(abs(est.lppd - lppd_grid) < max(0.05, 3 * (est.mc_se or 0)))}


def psis_pilot(data: L.Data, rng: np.random.Generator, budget: Budget, evals: int = 300) -> Dict[str, Any]:
    """The one-unit PSIS pilot at full dimension (sec. 13 item 2): U3 under R2 on SYNTHETIC data (a truth drawn from
    R2's prior reweighted by the non-U3 fed units), conditioned on U4, U1e, U1l. The mode of log(L x p) is searched
    from three starts with fitting.fit (central polish, require_converged='flag'), <= ``evals`` evaluations each; a
    Student-t (5 df) proposal at the best mode, scaled from the integrand's width; S = 256 -> 4096. Records k-hat,
    ESS, wall time and the projection to the full PSIS plan (385 cells)."""
    from src.kinetic_core import fitting
    t_start = time.time()
    rival = "R2"
    cs = L.coords(rival)
    d = len(cs)
    obs = real_obs(data)
    prior = make_draw_set(rival, L.PRIMARY, data, 256, seed=901, budget=budget)
    lw = weights_from(prior, ("U4", "U1e", "U1l"), obs, L.PRIMARY, data)
    j = int(sample_truths(prior, lw, 1, rng)[0])
    truth_u3 = prior.points["U3"][j]
    u3 = synth_unit(L.u3_unit(10.0 ** truth_u3), truth_u3, rng)
    obs_syn = Obs({**obs.units, "U3": u3}, obs.uarp_y)
    units = ("U4", "U1e", "U1l", "U3")
    covs = {u: L.unit_covariance(obs_syn.units[u]) for u in units}
    chol = {u: np.linalg.cholesky(covs[u]) for u in units if u != "U3"}

    def residuals_many(points_u: np.ndarray) -> np.ndarray:
        ev = L.evaluate(np.clip(points_u, 0.0, 1.0), rival, L.PRIMARY, data)
        res = []
        for u in units:
            pts = L.points_for(u, ev)
            if u == "U3":
                ll = L.unit_loglik(obs_syn.units[u], pts)
                res.append(np.sqrt(np.maximum(-2.0 * (ll - _ll_max_u3), 0.0))[:, None])
            else:
                un = obs_syn.units[u]
                r = (un.y[None, :] - pts @ un.design.T)
                res.append(np.linalg.solve(chol[u], r.T).T)
        r = np.concatenate(res, axis=1)
        return np.where(ev.feasible[:, None] & np.isfinite(r), r, 25.0)

    _ll_max_u3 = float(np.max(L.unit_loglik(u3, prior.points["U3"])))
    k_res = residuals_many(np.full((1, d), 0.5)).shape[1]
    prob = fitting.Problem(keys=tuple(c.name for c in cs), kinds=("acid_yield",) * d, lower=np.zeros(d),
                           upper=np.ones(d), rows=tuple(range(k_res)),
                           residuals=lambda th, rows: residuals_many(th[None, :])[0],
                           evaluate_many=lambda pts, rows: list(residuals_many(np.array(pts))))
    starts = [np.full(d, 0.5), prior.u[int(np.argmax(lw))], prior.u[j]]
    fits = []
    for s in starts:
        budget.check("PSIS pilot mode search")
        f = fitting.fit(prob, s, budget=evals, central_polish=1, polish_budget=max(2 * d + 2, 80),
                        require_converged="flag")
        fits.append(f)
    best = min(fits, key=lambda f: f.cost)
    t_mode = time.time() - t_start
    # proposal: scale from the integrand's width (the prior draws' weighted spread near the mode), never only the
    # point curvature (lesson 2026-11-05)
    lw_all = lw + L.unit_loglik(u3, prior.points["U3"])
    wn = np.exp(lw_all - np.max(lw_all))
    wn = wn / wn.sum()
    spread = np.cov(prior.u.T, aweights=wn) + 1e-4 * np.eye(d)
    out_s = {}
    draws_u, lq_all, li_all = [], [], []
    for s_target in (256, 1024, 4096):
        budget.check("PSIS pilot sampling")
        need = s_target - sum(len(x) for x in draws_u)
        x, lq = pd.multivariate_t_draws(best.theta, spread, need, rng)
        inside = np.all((x >= 0) & (x <= 1), axis=1)
        ev = L.evaluate(np.clip(x, 0, 1), rival, L.PRIMARY, data)
        li = np.zeros(need)
        for u in units:
            li = li + L.unit_loglik(obs_syn.units[u], L.points_for(u, ev))
        li = np.where(inside & ev.feasible, li, -np.inf)
        draws_u.append(x)
        lq_all.append(lq)
        li_all.append(li)
        est = pd.importance_lppd(np.concatenate(li_all), np.concatenate(lq_all))
        out_s[s_target] = {"k_hat": est.k_hat, "ess": est.ess, "accepted": est.accepted, "lppd": est.lppd,
                           "mc_se": est.mc_se, "reason": est.reason, "nonfinite": est.n_nonfinite}
        if est.accepted:
            break
    wall = time.time() - t_start
    n_eval = sum(f.evaluations for f in fits) + sum(len(x) for x in draws_u)
    per_cell_h = wall / 3600.0
    return {"rival": rival, "dimension": d, "unit": "U3 | U4, U1e, U1l (synthetic U3)",
            "mode_search": [{"cost": f.cost, "status": f.status, "evaluations": f.evaluations,
                             "converged": f.converged, "polish_flagged": (f.polish or {}).get("flagged")}
                            for f in fits],
            "proposal_scale": "the prior draws' U3 x fed-weighted spread (integrand width), Student-t 5 df",
            "psis_by_S": out_s, "mode_search_wall_s": round(t_mode, 1), "wall_s": round(wall, 1),
            "evaluations": n_eval, "projection_385_cells_h_serial": round(385 * per_cell_h, 1),
            "projection_385_cells_h_wall_2_workers": round(385 * per_cell_h / 1.8, 1),
            "exceeds_24h_cap": bool(385 * per_cell_h / 1.8 > 24.0)}


def blind_order_stub() -> Dict[str, Any]:
    """The V-blind check at build time: nothing is keyed or frozen yet, so the scorer must refuse."""
    out = {"keyed_file": KEYED, "keyed_file_present_on_disk": (ROOT / KEYED).exists()}
    try:
        sha = freeze_check()["freeze_commit"]
        out["freeze"] = sha
        out["v_blind"] = v_blind(sha)
    except Refusal as e:
        out["refusal"] = str(e)
        out["score_refused"] = True
    return out


def r1_drift_diagnostic(data: L.Data) -> Dict[str, Any]:
    """Sec. 3.4: R1's pH-drift effect in the 0.1 M pots, printed beside the units (diagnostic, central draws)."""
    u = L.sobol_unit(128, len(L.coords("R1")), 13)
    a = L.evaluate(u, "R1", L.PRIMARY, data, trapped=False)
    b = L.evaluate(u, "R1", L.Form(ph_drift_fed=True), data, trapped=False)
    f = a.feasible & b.feasible
    d1 = (b.pred["U1"] - a.pred["U1"])[f]
    d4 = (b.pred["U4"] - a.pred["U4"])[f]
    return {"draws": int(f.sum()), "U1_points_dex_median_abs": float(np.median(np.abs(d1))),
            "U1_points_dex_max_abs": float(np.max(np.abs(d1))), "U4_points_dex_median_abs": float(np.median(np.abs(d4))),
            "fig1_pH_500_min_median": float(np.median(b.diag["fig1_PH_500"][f]))}


def run_checks(args) -> int:
    budget = Budget(args.cap_hours or 6.0)
    mem = memory_guard(not args.no_memory_guard)
    data = L.load_data()
    rng = np.random.default_rng(args.seed)
    out: Dict[str, Any] = {"memory": mem, "ph_scale": L.PH_SCALE}
    va = L.v_acid_table()
    diffs = [x for r in va for x in r.get("abs_diff", [])]
    out["v_acid"] = {"rows": va, "max_abs_diff_pH": max(diffs), "tolerance": 0.02,
                     "pass": bool(max(diffs) <= 0.02 and all(r.get("pass", True) for r in va))}
    out["v1_ex049"] = v1_ex049(data)
    out["v2_balance"] = v2_balance(data)
    out["v4_r0_pinned"] = L.r0_pinned_constants()
    out["v5_ph_scale"] = {"scale": pd.require_ph_scale(L.PH_SCALE, L.PH_SCALE),
                          "note": "every dpH is on the 25 C label scale from label pH 7.0 (rule 13)"}
    out["v_time"] = v_time(data)
    out["r1_drift_in_0.1M_pots"] = r1_drift_diagnostic(data)
    out["v_density"] = v_density_toy(data, rng)
    out["blind_order_stub"] = blind_order_stub()
    try:
        out["psis_pilot"] = psis_pilot(data, rng, budget)
    except Refusal as e:
        out["psis_pilot"] = {"refused": str(e)}
    out["compute_wall_h"] = round(budget.elapsed_h(), 3)
    out["declared_choices"] = DECLARED
    _record("checks", out, render_checks(out), args.out_dir, layer_inputs())
    return 0


DECLARED = [line.strip() for line in (L.__doc__ or "").split("DECLARED CHOICES")[1].splitlines()
            if re.match(r"^\s*D-[a-z]", line)] if "DECLARED CHOICES" in (L.__doc__ or "") else []


def render_checks(o: Dict[str, Any]) -> str:
    L_ = ["# M3 phase 1: build checks (V-acid, V1, V2, V4, V5, V-time, V-density, PSIS pilot, V-blind stub)", "",
          f"*Generated by `scripts/generators/fit_m3_strecker_layer.py --checks`. Design: {PREREG} (v2 DRAFT, not "
          "frozen). Layer: `src/m3_strecker_layer.py`. Nothing here scores a decisive unit; the trapped pot ran on "
          "synthetic inputs only.*", ""]
    va = o["v_acid"]
    L_ += [f"## V-acid: {'PASS' if va['pass'] else 'FAIL'} (max |layer - prereg| {va['max_abs_diff_pH']:.4f} pH, "
           "tolerance 0.02)", "", "| pot | acid mM | prereg drop | layer drop |", "|---|---|---|---|"]
    for r in va["rows"]:
        L_.append(f"| {r['pot']} | {r['acid_mM']} | {r['printed_drop']} | "
                  f"{r.get('layer_drop', r.get('layer_ph'))} |")
    v1, v2 = o["v1_ex049"], o["v2_balance"]
    L_ += ["", f"## V1 (EX-049 at its structure limit): route-A calibration bit-exact {v1['route_a_calibration_bit_exact']}"
           f"; glucose-pot dicarbonyl part bit-exact {v1['glucose_pot_dicarbonyl_bit_exact']}", "",
           f"## V2 (Phe, sugar and O2 balance, relative): {'PASS' if v2['pass'] else 'FAIL'} at {v2['tolerance']:g}", "",
           "```", json.dumps(v2["per_rival_form"], indent=1), "```", "",
           f"## V4 (R0-pinned read from the code): {json.dumps(o['v4_r0_pinned'])}", "",
           f"## V5: pH scale {o['v5_ph_scale']['scale']} ({o['v5_ph_scale']['note']})", "",
           "## V-time (one process, nice'd, 256-draw batches)", "", "| rival | coordinates | ms / draw | feasible share | "
           "failed |", "|---|---|---|---|---|"]
    for r, v in o["v_time"].items():
        L_.append(f"| {r} | {v['coordinates']} | {v['ms_per_draw']} | {v['feasible_share']:.3f} | {v['failed']} |")
    dr = o["r1_drift_in_0.1M_pots"]
    L_ += ["", f"## R1's pH drift in the 0.1 M pots (diagnostic): U1 points median |d| {dr['U1_points_dex_median_abs']:.3f}"
           f" dex (max {dr['U1_points_dex_max_abs']:.3f}); U4 points median {dr['U4_points_dex_median_abs']:.3f}; Fig. 1"
           f" pot pH at 500 min {dr['fig1_pH_500_min_median']:.2f}", ""]
    vd = o["v_density"]
    L_ += [f"## V-density (two-coordinate R3 toy, U-ARPglc): {'PASS' if vd['pass'] else 'FAIL'}", "",
           f"trapezoid 41x41 {vd['lppd_trapezoid_41x41']}; PSIS {vd['lppd_psis']} (k-hat {vd['k_hat']:.3f}, ESS "
           f"{vd['ess']:.0f}); Gauss-Hermite 20x20 {vd['lppd_gauss_hermite_20']}", ""]
    pp = o["psis_pilot"]
    L_ += ["## One-unit PSIS pilot at full dimension (U3 under R2, synthetic U3)", "", "```", json.dumps(pp, indent=1,
                                                                                                       default=str),
           "```", "", "## V-blind (blind-order) stub", "", "```", json.dumps(o["blind_order_stub"], indent=1), "```", "",
           "## Declared choices where the prereg is silent", ""] + [f"- {d}" for d in o["declared_choices"]]
    L_ += ["", f"Memory at launch: {o['memory']}; compute {o['compute_wall_h']} h wall."]
    return "\n".join(L_)


# ============================================================================================ the screen

SIG_BEARS = {"R1": ("stoich", "ARP", "cross", "early"), "R2": ("stoich", "O2", "lag", "fp", "cross", "early"),
             "R3": ("stoich", "cross", "early"), "R4": ("stoich", "cross", "early")}
SCREEN_FORMS = {"primary": L.PRIMARY, "cf073": L.Form(cf073=True), "voigt": L.Form(voigt=True),
                "o_ar_0.1": L.Form(o_ar=0.1), "headspace_100mL": L.Form(headspace_ml=100.0),
                "unbuffered_trapped": L.Form(unbuffered_trapped=True, ph_drift_trapped=True),
                "kd_floor": L.Form(kd_floor=True), "isomers": L.Form(isomers=True),
                "cf070": L.Form(cf070=True), "early_step": L.Form(early_step=True), "t0_5": L.Form(t0=5.0)}
O2_FORMS = {f"headspace_{int(h)}mL_nu{nu}": L.Form(headspace_ml=h, nu_o2=nu) for h in L.HEADSPACES_ML
            for nu in L.NU_O2}


def screen_signatures(ev: L.Evaluation, data: L.Data, form: L.Form) -> Dict[str, np.ndarray]:
    """Per draw, pass / fail of each sec. 8 signature (sigma-trap excluded: after keying)."""
    lg = math.log10
    n = len(ev.feasible)
    out = {"stoich": ev.feasible.copy()}
    om2 = (L.OMEGA_DEX * form.omega_scale) ** 2
    # sigma-ARP: every Table 3 interval within its noise tolerance, under the noise forms
    y = L.arp_interval_measured(data, ev.theta.get("f_iso") if form.isomers else None)
    ts = sorted(data.arp_glc)
    umol = np.array([data.arp_glc[t] * 10.0 for t in ts])
    for label, cv, between in (("ARP_cv2", 0.02, False), ("ARP_cv10", 0.10, False), ("ARP_cv10_between", 0.10, True)):
        var = cv ** 2 + ((0.5 / math.sqrt(3.0)) / umol) ** 2 + ((L.SIGMA_RUN_DEX * math.log(10)) ** 2 if between else 0)
        s = np.array([math.sqrt(var[ts.index(a)] + var[ts.index(b)]) for a, b in L.ARP_INTERVALS])
        out[label] = ev.feasible & np.all(ev.pred["UARP_c"] - y >= -Z975 * s[None, :], axis=1)
    out["ARP"] = out["ARP_cv10_between"]
    # sigma-lag: Fig. 1 log PA(10)/PA(30), -0.42 measured; tolerance z sqrt(var + 2 omega^2) (vacuous: declared)
    f1 = data.fig1_arp["PHEAL"]
    meas_lag = lg(f1[10.0][0] / f1[30.0][0])
    var = f1[10.0][1] ** 2 + f1[30.0][1] ** 2
    lag = ev.pred["U1"][:, 0] - ev.pred["U1"][:, 1]
    out["lag"] = ev.feasible & (np.abs(lag - meas_lag) <= Z975 * math.sqrt(var + 2 * om2))
    # sigma-fp: Table 2 air + Cu PA / PAA 1.8
    a, b = data.arp_t2["air_cu5"]
    ea, eb = data.arp_t2_err["air_cu5"]
    fp = ev.pred["U4"][:, 2] - ev.pred["U4"][:, 5]
    out["fp"] = ev.feasible & (np.abs(fp - lg(a / b)) <= Z975 * math.sqrt(ea ** 2 + eb ** 2 + 2 * om2))
    # sigma-cross: PAA overtakes PA between 60 and 120 min (Fig. 2); tolerance on log PA/PAA at 60 and 120
    pa = ev.pred["U5_PA"]
    paa = ev.pred["U5_PAA"]
    ok = ev.feasible.copy()
    for t in (60.0, 120.0):
        j = L.GLC_T.index(t)
        m = lg(data.fig2["PHEAL"][t][0] / data.fig2["PHEAC"][t][0])
        v = data.fig2["PHEAL"][t][1] ** 2 + data.fig2["PHEAC"][t][1] ** 2
        ok &= np.abs((pa[:, j] - paa[:, j]) - m) <= Z975 * math.sqrt(v + 2 * om2)
    out["cross"] = ok
    cr = ev.diag.get("glc_crossover", np.full(n, np.inf))
    out["cross_in_window_binary"] = ev.feasible & (cr > 60.0) & (cr <= 120.0)
    # sigma-early: Fig. 2 at 30 min, log total(30)/total(120)
    tot = np.log10(10.0 ** pa + 10.0 ** paa)
    m = lg((data.fig2["PHEAL"][30.0][0] + data.fig2["PHEAC"][30.0][0]) /
           (data.fig2["PHEAL"][120.0][0] + data.fig2["PHEAC"][120.0][0]))
    out["early"] = ev.feasible & (np.abs(tot[:, L.GLC_T.index(30.0)] - tot[:, L.GLC_T.index(120.0)] - m)
                                  <= Z975 * math.sqrt(0.05 ** 2 + 2 * om2))
    return out


def run_screen(args) -> int:
    budget = Budget(args.cap_hours or 6.0)
    mem = memory_guard(not args.no_memory_guard)
    data = L.load_data()
    n = max(600, args.draws)
    res: Dict[str, Any] = {}
    failed_total = 0
    forms_for = {"R1": ("primary", "kd_floor", "isomers", "cf070", "early_step", "t0_5", "o_ar_0.1",
                        "headspace_100mL"),
                 "R2": ("primary", "cf073", "voigt", "o_ar_0.1", "headspace_100mL", "unbuffered_trapped",
                        "kd_floor", "t0_5"),
                 "R3": ("primary", "kd_floor"), "R4": ("primary", "kd_floor", "isomers")}
    for rival in ("R1", "R2", "R3", "R4"):
        res[rival] = {}
        forms = {k: SCREEN_FORMS[k] for k in forms_for[rival]}
        if rival in ("R1", "R2"):
            forms.update(O2_FORMS)
        for fname, form in forms.items():
            budget.check(f"screen {rival} {fname}")
            t = time.time()
            u = L.sobol_unit(n, len(L.coords(rival, form)), args.seed + 31 * RIVAL_INDEX[rival])
            sigs, feas, failed, diags = {}, [], 0, []
            for s0 in range(0, n, 300):
                ev = L.evaluate(u[s0:s0 + 300], rival, form, data, trapped=False)
                for k, v in screen_signatures(ev, data, form).items():
                    sigs.setdefault(k, []).append(v)
                failed += ev.failed
                diags.append({k: ev.diag[k] for k in ("sigma_ox", "c_cu5", "glc_ARP_share_PA_120", "k_glc_eff_45")})
                for k, v in ev.reasons.items():
                    feas.append((k, int(np.sum(v))))
            sigs = {k: np.concatenate(v) for k, v in sigs.items()}
            bears = SIG_BEARS[rival] if fname not in O2_FORMS else ("stoich",)
            joint = np.all([sigs[k] for k in bears], axis=0)
            reasons: Dict[str, int] = {}
            for k, v in feas:
                reasons[k] = reasons.get(k, 0) + v
            dg = {k: np.concatenate([d[k] for d in diags]) for k in diags[0]}
            res[rival][fname] = {"draws": n, "pass": {k: int(v.sum()) for k, v in sigs.items()},
                                 "bears_on": list(bears), "joint_pass": int(joint.sum()),
                                 "coverage_reported": bool(joint.sum() < 10), "failed_runs": failed,
                                 "infeasible": reasons, "wall_s": round(time.time() - t, 1),
                                 "sigma_ox_q10_50_90": _q(dg["sigma_ox"]), "c_cu5_q10_50_90": _q(dg["c_cu5"]),
                                 "glc_route_B_share_PA_120_q": _q(dg["glc_ARP_share_PA_120"][joint]) if joint.any()
                                 else None, "k_glc_eff_45_q": _q(dg["k_glc_eff_45"][joint]) if joint.any() else None}
            failed_total += failed
            print(rival, fname, res[rival][fname]["joint_pass"], "/", n, flush=True)
    # R0-pinned: its glucose-pot PA printed (sigma-R0)
    u0 = L.sobol_unit(n, len(L.coords("R0")), args.seed)
    ev0 = L.evaluate(u0, "R0", L.PRIMARY, data)
    meas = {t: data.fig2["PHEAL"][t][0] for t in L.FIG2_T}
    res["R0"] = {"PA_fig2_times_median_mM": {t: float(np.median(10.0 ** ev0.pred["U5a_PA_shape"][:, i]))
                                             for i, t in enumerate(L.FIG2_T)}, "measured_PA_mM": meas}
    kills, demotions = screen_kills(res)
    out = {"memory": mem, "draws_per_rival_form": n, "results": res, "kills": kills, "implied_demotions": demotions,
           "failed_runs_total": failed_total, "refused": failed_total > 0, "compute_wall_h": round(budget.elapsed_h(), 3),
           "excluded": "sigma-trap (runs after the keying commit, demotion only)"}
    _record("screen", out, render_screen(out), args.out_dir, layer_inputs())
    return 0


RIVAL_INDEX = {r: i for i, r in enumerate(L.RIVALS)}


def _q(x) -> Optional[List[float]]:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return [round(float(v), 4) for v in np.quantile(x, [0.1, 0.5, 0.9])] if x.size else None


def screen_kills(res: Dict[str, Any]) -> Tuple[Dict[str, Any], List[str]]:
    """Sec. 8's kill rule: no draw reproduces every signature bearing on the rival, under EVERY form that bears on it
    (R2: CF-073, Voigt, o_Ar 0.1, 100 mL, the unbuffered trapped pot; and every finite-O2 form for sigma-O2). A kill
    found by sampling still needs EX-049's bounded violation search before it is final (flagged)."""
    kills, dem = {}, []
    for rival in ("R1", "R2", "R3", "R4"):
        r = res[rival]
        main = [f for f in r if not f.startswith("headspace_") or f == "headspace_100mL"]
        killed = all(r[f]["joint_pass"] == 0 for f in main)
        o2_killed = all(r[f]["joint_pass"] == 0 for f in r if f.startswith("headspace_") and f != "headspace_100mL")
        kills[rival] = {"killed_by_sampling_under_every_form": bool(killed),
                        "forms_with_zero_joint_pass": [f for f in r if r[f]["joint_pass"] == 0],
                        "needs_violation_search_confirmation": bool(killed),
                        "sigma_O2_kill_every_budget_form": bool(o2_killed) if rival in ("R1", "R2") else None,
                        "feasible_share_primary": round(1 - sum(r["primary"]["infeasible"].values()) /
                                                        r["primary"]["draws"], 3),
                        "precondition_10pct": bool(1 - sum(r["primary"]["infeasible"].values()) /
                                                   r["primary"]["draws"] >= 0.10)}
        if not kills[rival]["precondition_10pct"]:
            dem.append(f"{rival}: feasible calibration in < 10 % of draws (sec. 9.5 precondition): every contrast "
                       f"with {rival} reads 'not separated on disk' unless more draws restore it")
    if kills["R2"]["killed_by_sampling_under_every_form"]:
        dem.append("route: R2 killed by the screen; credits R1 only if R1 passes adequacy on U3, U1e, U4 (M3)")
    if kills["R1"]["killed_by_sampling_under_every_form"]:
        dem.append("route / transfer / decay law: R1 killed; outcome rows 1 - 3 unreachable")
    return kills, dem


def render_screen(o: Dict[str, Any]) -> str:
    L_ = ["# M3 phase 1: the mechanism screen (sec. 8; sigma-trap excluded)", "",
          f"*Generated by `scripts/generators/fit_m3_strecker_layer.py --screen`; {o['draws_per_rival_form']} scrambled "
          "Sobol draws per rival and form, every coordinate drawn, sigma_ox and c_Cu5 inverted per draw. Shape "
          "signatures carry omega (near-vacuous, declared); feasibility signatures are hard. SINGLE LABORATORY.*", "",
          f"**Failed runs: {o['failed_runs_total']}{' -- REFUSED' if o['refused'] else ''}.**", "",
          "| rival | form | stoich | O2 budget | ARP (cv 2 / 10 / 10+between) | lag | fp | cross (tol / binary) | early |"
          " joint (bearing) | infeasible |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for rival in ("R1", "R2", "R3", "R4"):
        for f, v in o["results"][rival].items():
            p = v["pass"]
            L_.append(f"| {rival} | {f} | {p['stoich']} | {'-' if not f.startswith('headspace') else p['stoich']} | "
                      f"{p['ARP_cv2']} / {p['ARP_cv10']} / {p['ARP_cv10_between']} | {p['lag']} | {p['fp']} | "
                      f"{p['cross']} / {p['cross_in_window_binary']} | {p['early']} | {v['joint_pass']} / {v['draws']}"
                      f" ({', '.join(v['bears_on'])}) | {v['infeasible']} |")
    L_ += ["", "## Kills (sampling) and implied demotions", "", "```", json.dumps(o["kills"], indent=1), "```", ""]
    L_ += [f"- {d}" for d in o["implied_demotions"]] or ["- none implied"]
    L_ += ["", f"R0-pinned glucose-pot PA (median, mM) at Fig. 2's times: {json.dumps(o['results']['R0'])}", "",
           f"Memory at launch: {o['memory']}; compute {o['compute_wall_h']} h wall. {o['excluded']}."]
    return "\n".join(L_)


# ============================================================================================ the power checks

def binomial_decision(hits: int, false_hits: int, r: int, alpha: float = 0.05) -> Dict[str, Any]:
    """Sec. 12.2 (M8): power must exceed the false-credit rate by a one-sided exact two-proportion test (Fisher) at
    alpha; Clopper-Pearson 95 % intervals; power < 0.5 is 'low power'."""
    from scipy.stats import beta, fisher_exact

    def cp(k):
        lo = 0.0 if k == 0 else float(beta.ppf(0.025, k, r - k + 1))
        hi = 1.0 if k == r else float(beta.ppf(0.975, k + 1, r - k))
        return [round(lo, 3), round(hi, 3)]

    p = float(fisher_exact([[hits, r - hits], [false_hits, r - false_hits]], alternative="greater")[1])
    return {"power": hits / r, "power_ci95": cp(hits), "false_credit": false_hits / r, "false_credit_ci95": cp(false_hits),
            "fisher_p_one_sided": p, "decidable": bool(p < alpha), "low_power": bool(hits / r < 0.5), "R": r}


def power_u3(data: L.Data, sets: Mapping[Tuple[str, str], DrawSet], truths: Mapping[str, DrawSet], r: int,
             rng: np.random.Generator, n_tie: int, budget: Budget) -> Dict[str, Any]:
    """Sec. 12.1: U3 alone, route rule (R1 vs R2), on SYNTHETIC data at the design times; truths from each rival's
    prior reweighted by the non-U3 fed units (U4, U1e, U1l; Table 2 is the calibration); scored under the primary
    and the Cu-law form (required agreement), conditioned on the real non-U3 fed units."""
    obs = real_obs(data)
    out: Dict[str, Any] = {}
    readings = {}
    for truth in ("R1", "R2"):
        tset = truths[truth]
        lw = weights_from(tset, ("U4", "U1e", "U1l"), obs, L.PRIMARY, data)
        picks = sample_truths(tset, lw, r, rng)
        per = []
        for i, j in enumerate(picks):
            budget.check("U3 power")
            u3 = synth_unit(L.u3_unit(10.0 ** tset.points["U3"][j]), tset.points["U3"][j], rng)
            o = Obs({**obs.units, "U3": u3}, obs.uarp_y)
            rec = {}
            for fname in ("primary", "cu_linear"):
                v = unit_verdict(sets[("R1", fname)], sets[("R2", fname)], "U3", ("U4", "U1e", "U1l"), o, data,
                                 POWER_FORMS[fname], rng, n_tie)
                rec[fname] = v
            per.append(rec)
        readings[truth] = per
    for fname_set, label in ((("primary",), "primary_only"), (("primary", "cu_linear"), "primary_and_cu_law")):
        for psis_rule in (True, False):
            key = f"{label}{'' if psis_rule else '_ignoring_resolution'}"
            res = {}
            for credit in ("a", "b"):
                def credited(rec):
                    rs = [_reading(rec[f], psis_rule) for f in fname_set]
                    return all(x == credit for x in rs)
                hits = sum(credited(rec) for rec in readings["R1" if credit == "a" else "R2"])
                false = sum(credited(rec) for rec in readings["R2" if credit == "a" else "R1"])
                res[f"credit_{'R1' if credit == 'a' else 'R2'}"] = binomial_decision(hits, false, r)
            out[key] = res
    out["resolution_rate"] = {t: float(np.mean([rec["primary"]["resolved"] for rec in readings[t]])) for t in readings}
    out["d_quantiles"] = {t: _q([rec["primary"]["d"] for rec in readings[t]]) for t in readings}
    out["band_quantiles"] = {t: _q([rec["primary"]["band"] for rec in readings[t]]) for t in readings}
    out["k_hat_quantiles"] = {t: _q([max(rec["primary"]["k_hat_a"] or 0, rec["primary"]["k_hat_b"] or 0)
                                     for rec in readings[t]]) for t in readings}
    out["truth_reweighting"] = "U4, U1e, U1l on real data (Table 2 enters through the calibrations); never U3"
    return out


def _reading(v: Dict[str, Any], psis_rule: bool) -> str:
    if psis_rule:
        return v["reading"]
    return v["reading_raw"]


def power_pairs(data: L.Data, sets: Mapping[Tuple[str, str], DrawSet], truths: Mapping[str, DrawSet], r: int,
                rng: np.random.Generator, n_tie: int, budget: Budget,
                contrasts: Sequence[str] = tuple(CONTRASTS)) -> Dict[str, Any]:
    """Sec. 12.2: the three pairs. Synthetic claim units from a truth (truths reweighted by the real units outside the
    pair's claim set); the non-claim conditioning units real, except U3 (synthetic from the same truth: it is blind).
    Scored under the primary and the three forms with precomputed draw sets (o_Ar 0.1, Cu law, 20 mL headspace);
    the other forms are 'primary-only, optimistic'. Holm's guard step is not simulated (optimistic)."""
    obs = real_obs(data)
    out: Dict[str, Any] = {}
    pair_sets = {"route": ("R1", "R2"), "transfer": ("R1", "R3"), "decay_law": ("R4", "R1drift")}
    reweight = {"route": ("U1l",), "transfer": ("U4", "U1e", "U1l"), "decay_law": ("U4", "U1e", "U1l")}
    for cname, spec in CONTRASTS.items():
        if cname not in contrasts:
            continue
        a_name, b_name = pair_sets[cname]
        readings: Dict[str, List[Dict[str, str]]] = {}
        for truth in (a_name, b_name):
            tset = truths[truth]
            lw = weights_from(tset, reweight[cname], obs, L.PRIMARY, data)
            picks = sample_truths(tset, lw, r, rng)
            per = []
            for j in picks:
                budget.check(f"power {cname}")
                o = obs
                syn_units = dict(o.units)
                for u in spec["claim"]:
                    if u in ("U3",):
                        syn_units[u] = synth_unit(L.u3_unit(10.0 ** tset.points["U3"][j]), tset.points["U3"][j], rng)
                    elif u != "UARP":
                        syn_units[u] = synth_unit(o.units[u], tset.points[u][j], rng)
                if "U3" not in spec["claim"]:
                    syn_units["U3"] = synth_unit(L.u3_unit(10.0 ** tset.points["U3"][j]), tset.points["U3"][j], rng)
                uarp_y = synth_uarp(data, tset.k_glc[j], rng) if "UARP" in spec["claim"] else o.uarp_y
                os_ = Obs(syn_units, uarp_y)
                rec = {}
                for fname, form in POWER_FORMS.items():
                    a, b = sets[(a_name, fname)], sets[(b_name, fname)]
                    rec[fname] = _contrast_on(cname, spec, a, b, os_, data, form, rng, n_tie)
                per.append(rec)
            readings[truth] = per
        res = {}
        for label, forms in (("primary_only", ("primary",)), ("four_forms_required_agreement", tuple(POWER_FORMS))):
            for credit, truth_hit, truth_false in (("a", a_name, b_name), ("b", b_name, a_name)):
                hits = sum(all(rec[f] == credit for f in forms) for rec in readings[truth_hit])
                false = sum(all(rec[f] == credit for f in forms) for rec in readings[truth_false])
                res[f"{label}|credit_{truth_hit}"] = binomial_decision(hits, false, r)
        out[cname] = {"pair": [a_name, b_name], "claim_units": list(spec["claim"]), "truth_reweighting": reweight[cname],
                      "results": res,
                      "decidable_any_direction": any(v["decidable"] for k, v in res.items()
                                                     if k.startswith("four_forms")),
                      "demoted": not any(v["decidable"] for k, v in res.items() if k.startswith("four_forms"))}
    return out


def _contrast_on(cname: str, spec, a: DrawSet, b: DrawSet, obs: Obs, data: L.Data, form: L.Form,
                 rng: np.random.Generator, n_tie: int) -> str:
    claim = spec["claim"]
    if cname == "route":
        cond = {"U4": (), "U1e": ("U4",), "U3": ("U4", "U1e", "U1l")}
        per = {u: unit_verdict(a, b, u, cond[u], obs, data, form, rng, n_tie) for u in claim}
        return contrast_credit(per["U3"], per)
    fed = ("U4", "U1e", "U1l", "U3") if cname == "decay_law" else ("U4", "U1e", "U1l")
    if cname == "decay_law":
        block, cond0 = ("U3", "U5a", "U5b", "UARP"), ("U4", "U1e", "U1l")
        per = {}
        prev: Tuple[str, ...] = cond0
        for u in block:
            per[u] = unit_verdict(a, b, u, prev, obs, data, form, rng, n_tie)
            prev = prev + (u,)
        prim = block_verdict(a, b, block, cond0, obs, data, form, rng, n_tie)
        return contrast_credit(prim, per)
    # transfer: (F) with U3 in the conditioning, (T) without; required agreement
    sides = []
    for cond0 in (fed + ("U3",), fed):
        per = {}
        prev = cond0
        for u in GLC_ORDER:
            per[u] = unit_verdict(a, b, u, prev, obs, data, form, rng, n_tie)
            prev = prev + (u,)
        prim = block_verdict(a, b, GLC_ORDER, cond0, obs, data, form, rng, n_tie)
        sides.append(contrast_credit(prim, per))
    return sides[0] if sides[0] == sides[1] else "none"


def run_power(args) -> int:
    budget = Budget(args.cap_hours or 6.0)
    mem = memory_guard(not args.no_memory_guard)
    data = L.load_data()
    rng = np.random.default_rng(args.seed + 1)
    t0 = time.time()
    sets: Dict[Tuple[str, str], DrawSet] = {}
    truths: Dict[str, DrawSet] = {}
    set_log = {}
    drift = L.Form(ph_drift_trapped=True)
    rivals = {"R1": L.PRIMARY, "R2": L.PRIMARY, "R3": L.PRIMARY, "R4": L.PRIMARY, "R1drift": drift}
    part = args.power_part
    need = {"route": ("R1", "R2"), "transfer": ("R1", "R3"), "decay_law": ("R4", "R1drift")}
    want = ("R1", "R2") if part == "u3" else tuple(dict.fromkeys(
        (["R1", "R2"] if part == "all" else []) + [r for c in args.contrasts.split(",") for r in need[c]]))
    forms = tuple(POWER_FORMS)          # U3's check uses primary and cu_linear; the pairs all four
    for name in want:
        i = list(rivals).index(name)
        base = rivals[name]
        rival = "R1" if name == "R1drift" else name
        for fname in forms:
            k = list(POWER_FORMS).index(fname)
            f = replace(POWER_FORMS[fname], ph_drift_trapped=base.ph_drift_trapped)
            ds = cached_draw_set(args.cache_dir, rival, f, data, args.score_feasible, seed=args.seed + 1000 * i + 10 * k,
                                 budget=budget)
            sets[(name, fname)] = ds
            set_log[f"{name}|{fname}"] = {"feasible": ds.n, "drawn": ds.n_drawn, "failed": ds.n_failed,
                                          "wall_s": ds.wall_s, "infeasible": ds.reasons}
            print("set", name, fname, ds.n, "/", ds.n_drawn, "failed", ds.n_failed, f"{ds.wall_s}s", flush=True)
        tf = replace(L.PRIMARY, ph_drift_trapped=base.ph_drift_trapped)
        ts = cached_draw_set(args.cache_dir, rival, tf, data, args.truth_feasible, seed=args.seed + 7000 + 100 * i,
                             budget=budget, keep_k=True)
        truths[name] = ts
        set_log[f"{name}|truth"] = {"feasible": ts.n, "drawn": ts.n_drawn, "failed": ts.n_failed, "wall_s": ts.wall_s}
        print("truth", name, ts.n, "/", ts.n_drawn, flush=True)
    failed = sum(v["failed"] for v in set_log.values())
    out: Dict[str, Any] = {"memory": mem, "R": args.reps, "tie_redraws": args.tie_redraws, "draw_sets": set_log,
                           "failed_runs_total": failed, "sets_wall_s": round(time.time() - t0, 1)}
    if failed:
        out["refused"] = f"{failed} failed runs in the draw sets: every power result is refused"
    else:
        try:
            if part in ("u3", "all"):
                out["u3_synthetic_12_1"] = power_u3(data, sets, truths, args.reps, rng, args.tie_redraws, budget)
            if part in ("pairs", "all") and args.contrasts:
                out["pairs_12_2"] = power_pairs(data, sets, truths, args.reps, rng, args.tie_redraws, budget,
                                                contrasts=tuple(args.contrasts.split(",")))
        except Refusal as e:
            out["compute_capped"] = str(e)
    out["compute_wall_h"] = round(budget.elapsed_h(), 3)
    out["approximations"] = [
        "draw sets are scrambled-Sobol prior draws conditioned on the calibrations (feasible only), reweighted per "
        "synthetic set (sec. 12.2's 'by reweighting'); a density whose PSIS k-hat >= 0.7 or ESS < 100 reads 'not "
        "resolved' (a tie); results ignoring that rule are reported beside",
        f"tie bands from {args.tie_redraws} measurement-noise redraws per unit and set (the prereg's 4 000 for scoring)",
        "censored U3 values: product of conditional marginals (exact for one); adequacy leaves them out",
        "Holm's guard step is not simulated: power is optimistic on that account",
        "synthetic sets are generated from the truth's primary-form predictions and scored under each form",
        "a rival with no feasible draw under a form (its calibration infeasible there) has lppd -inf: the other side is "
        "read as favoured under that form without its adequacy test (flagged by_infeasibility; optimistic)",
    ]
    name = "power" if part != "u3" else "power_u3"
    if args.record_suffix:
        name += f"_{args.record_suffix}"
    _record(name, out, render_power(out), args.out_dir, layer_inputs())
    return 0


def render_power(o: Dict[str, Any]) -> str:
    L_ = ["# M3 phase 1: power checks (sec. 12.1 U3 synthetic; sec. 12.2 three pairs)", "",
          f"*Generated by `scripts/generators/fit_m3_strecker_layer.py --power`. R = {o['R']} synthetic sets per truth; "
          "binomial demotion rule: one-sided Fisher exact test of power against false credit at 0.05. U3 on SYNTHETIC "
          "data at 20 / 40 / 80 / 160 / 320 min (the keyed Figs. 2 - 3 were not opened). SINGLE LABORATORY.*", ""]
    if o.get("refused"):
        L_ += [f"**REFUSED: {o['refused']}**", ""]
    if o.get("compute_capped"):
        L_ += [f"**COMPUTE-CAPPED: {o['compute_capped']}**", ""]
    L_ += ["## Draw sets", "", "| set | feasible / drawn | failed | wall s |", "|---|---|---|---|"]
    for k, v in o["draw_sets"].items():
        L_.append(f"| {k} | {v['feasible']} / {v['drawn']} | {v['failed']} | {v['wall_s']} |")
    if "u3_synthetic_12_1" in o:
        u3 = o["u3_synthetic_12_1"]
        L_ += ["", "## Sec. 12.1: U3 alone (route rule, R1 vs R2)", "",
               "| scoring | credit | power [95 % CI] | false credit [95 % CI] | Fisher p | decidable | low power |",
               "|---|---|---|---|---|---|---|"]
        for key, res in u3.items():
            if not isinstance(res, dict) or "credit_R1" not in res:
                continue
            for c, v in res.items():
                L_.append(f"| {key} | {c} | {v['power']:.3f} {v['power_ci95']} | {v['false_credit']:.3f} "
                          f"{v['false_credit_ci95']} | {v['fisher_p_one_sided']:.3g} | {v['decidable']} | "
                          f"{v['low_power']} |")
        L_ += ["", f"resolution rate (primary): {u3['resolution_rate']}; d quantiles {u3['d_quantiles']}; tie band "
               f"quantiles {u3['band_quantiles']}; max k-hat quantiles {u3['k_hat_quantiles']}"]
    if "pairs_12_2" in o:
        L_ += ["", "## Sec. 12.2: the three pairs", "",
               "| contrast | scoring and credit | power [95 % CI] | false credit [95 % CI] | Fisher p | decidable |",
               "|---|---|---|---|---|---|"]
        for c, v in o["pairs_12_2"].items():
            for k, r in v["results"].items():
                L_.append(f"| {c} | {k} | {r['power']:.3f} {r['power_ci95']} | {r['false_credit']:.3f} "
                          f"{r['false_credit_ci95']} | {r['fisher_p_one_sided']:.3g} | {r['decidable']} |")
            L_.append(f"| {c} | **demoted: {v['demoted']}** | | | | |")
    L_ += ["", "## Approximations (declared)", ""] + [f"- {a}" for a in o["approximations"]]
    L_ += ["", f"Memory at launch: {o['memory']}; compute {o['compute_wall_h']} h wall."]
    return "\n".join(L_)


# =========================================================================== the decisive fits and scoring

def fit_cell(evaluate_points: Callable[[np.ndarray], Tuple[np.ndarray, np.ndarray]], d: int, starts: Sequence[np.ndarray],
             rng: np.random.Generator, evals: int = 300, s_schedule: Sequence[int] = (256, 1024, 4096),
             spread: Optional[np.ndarray] = None) -> Dict[str, Any]:
    """One (unit set, form, rival) cell of the decisive fits (sec. 9.1): the mode of log(L x p) on the unit cube from
    the three starts (fitting.fit, central polish, require_converged='flag'; <= ``evals`` evaluations each), then a
    Student-t (5 df) proposal at the best mode whose scale is the integrand's width (``spread``) or, failing it, the
    Laplace covariance widened x 4; S runs through ``s_schedule`` until PSIS accepts. ``evaluate_points`` maps
    (n, d) unit-cube points to (log L per draw, whitened residual matrix (n, k)): the residuals drive the mode
    search, the log L the density. Returns the draws, log q, log L and the acceptance."""
    from src.kinetic_core import fitting

    def res_many(pts):
        _, r = evaluate_points(np.clip(np.asarray(pts, float), 0.0, 1.0))
        return [row for row in r]

    k_rows = evaluate_points(np.full((1, d), 0.5))[1].shape[1]
    prob = fitting.Problem(keys=tuple(f"u{i}" for i in range(d)), kinds=("acid_yield",) * d, lower=np.zeros(d),
                           upper=np.ones(d), rows=tuple(range(k_rows)),
                           residuals=lambda th, rows: res_many(th[None, :])[0],
                           evaluate_many=lambda pts, rows: res_many(pts))
    fits = [fitting.fit(prob, s, budget=evals, central_polish=1, polish_budget=max(2 * d + 2, 20),
                        require_converged="flag") for s in starts]
    best = min(fits, key=lambda f: f.cost)
    if spread is None:
        # the unscaled Gauss-Newton covariance (J^T J)^-1 of the whitened residuals, widened x 4 (a kink at the mode
        # makes a point curvature too narrow: lesson 2026-11-05; PSIS's k-hat is the check)
        jac = fitting.central_jacobian(prob, best.theta, prob.lower, prob.upper)
        cov = np.linalg.pinv(jac.T @ jac)
        spread = 4.0 * np.minimum(cov, 1.0) + 1e-4 * np.eye(d)
    xs, lqs, lls = [], [], []
    est = None
    for s in s_schedule:
        need = s - sum(len(x) for x in xs)
        x, lq = pd.multivariate_t_draws(best.theta, spread, need, rng)
        inside = np.all((x >= 0) & (x <= 1), axis=1)
        ll, _ = evaluate_points(np.clip(x, 0, 1))
        xs.append(x)
        lqs.append(lq)
        lls.append(np.where(inside, ll, -np.inf))
        est = pd.importance_lppd(np.concatenate(lls), np.concatenate(lqs))
        if est.accepted:
            break
    return {"theta_mode": best.theta, "fits": [f.as_dict(prob.keys) for f in fits], "draws": np.concatenate(xs),
            "log_q": np.concatenate(lqs), "log_lik": np.concatenate(lls), "estimate": est,
            "not_resolved": not est.accepted}


def score_contrast(d_by_unit: Mapping[str, Mapping[str, float]], band_by_unit: Mapping[str, Mapping[str, float]],
                   adequate: Mapping[str, Mapping[str, bool]], primary: str, holm_p: Mapping[str, float]
                   ) -> Dict[str, Any]:
    """Sec. 10.1's rule on scored units: ``d_by_unit[unit][form]`` and the bands per form; required agreement over
    forms (predictive_densities.required_agreement), adequacy of the favoured side on every favouring unit, the
    primary claim unit (or block) favouring a, no claim unit favouring b, and Holm's guard step."""
    sides = {u: pd.required_agreement(d_by_unit[u], band_by_unit[u])["side"] for u in d_by_unit}
    for u, s in sides.items():
        if s in ("a", "b") and not adequate.get(u, {}).get(s, False):
            sides[u] = "both inadequate"
    holm = pd.holm(holm_p) if holm_p else {}
    guards_ok = all(v["rejected"] is False for v in holm.values()) if holm else True
    credit = "none"
    if sides.get(primary) == "a" and not any(s == "b" for s in sides.values()) and guards_ok:
        credit = "a"
    elif sides.get(primary) == "b" and not any(s == "a" for s in sides.values()) and guards_ok:
        credit = "b"
    return {"sides": sides, "holm": holm, "guards_ok": guards_ok, "credit": credit}


def outcome_row(route: str, transfer: str, decay: str, r1_beats_nulls_and_r0: bool) -> Tuple[str, str]:
    """Sec. 10.3, exclusive and in order. route: credit 'R1' / 'R2' / 'none'; transfer: 'R1' / 'R3' / 'none'
    (tie, inadequate, confounded, demoted); decay: 'R4' / 'none'."""
    if route == "R1" and transfer == "R1" and decay != "R4" and r1_beats_nulls_and_r0:
        return "1", "claim supported (single laboratory, PROVISIONAL)"
    if route == "R1" and decay == "R4":
        return "1b", "two routes; the decay moves with order / pH"
    if route == "R1" and transfer == "R3":
        return "2", "two routes; the fed decay does not transfer"
    if route == "R1":
        return "3", "two routes; transfer undecided"
    if route == "R2":
        return "4", "dicarbonyl only"
    return "5", "not separated on disk"


def run_fit(args) -> int:
    run_once(args.out_dir / f"{STEM}_fit.json")
    freeze = freeze_check()
    mem = memory_guard(not args.no_memory_guard)
    raise Refusal(f"--fit is built and tested but deliberately not wired to run in this build (freeze {freeze}, "
                  f"memory {mem}): the coordinator launches the decisive fits after the freeze")


def run_score(args) -> int:
    run_once(args.out_dir / f"{STEM}_score.json")
    freeze = freeze_check()
    vb = v_blind(freeze["freeze_commit"])
    raise Refusal(f"--score: V-blind passed ({vb}); scoring is launched by the coordinator after the fits")


# ================================================================================================== main

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    mode = ap.add_mutually_exclusive_group(required=True)
    for m in ("checks", "screen", "power", "fit", "score"):
        mode.add_argument(f"--{m}", action="store_true")
    ap.add_argument("--seed", type=int, default=20261105)
    ap.add_argument("--draws", type=int, default=600, help="screen draws per rival and form (>= 600)")
    ap.add_argument("--reps", type=int, default=200, help="synthetic sets per truth (R)")
    ap.add_argument("--tie-redraws", type=int, default=40)
    ap.add_argument("--score-feasible", type=int, default=1024)
    ap.add_argument("--truth-feasible", type=int, default=512)
    ap.add_argument("--power-part", choices=("u3", "pairs", "all"), default="all")
    ap.add_argument("--contrasts", default="route,transfer,decay_law", help="power: the 12.2 pairs to run")
    ap.add_argument("--record-suffix", default="", help="power: a part's record name suffix (merged later)")
    ap.add_argument("--cache-dir", type=Path, default=None, help="draw-set cache (scratch; never under data/)")
    ap.add_argument("--cap-hours", type=float, default=None)
    ap.add_argument("--no-memory-guard", action="store_true")
    ap.add_argument("--out-dir", type=Path, default=data_paths.VALIDATION_DIR)
    args = ap.parse_args(argv)
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    try:
        if args.checks:
            return run_checks(args)
        if args.screen:
            return run_screen(args)
        if args.power:
            return run_power(args)
        if args.fit:
            args.cap_hours = args.cap_hours or 24.0
            return run_fit(args)
        args.cap_hours = args.cap_hours or 24.0
        return run_score(args)
    except Refusal as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
