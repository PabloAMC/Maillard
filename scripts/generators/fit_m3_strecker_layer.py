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
    # v2.1 amendment: U3a / U3b added; the route primary is the {U3a, U3b} block (U3 stays a claim unit)
    "route": {"a": "R1", "b": "R2", "claim": ("U3", "U3a", "U3b", "U1e", "U4"), "primary": "U3a+U3b",
              "guards": ("U1l", "GT3")},
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


def memory_guard(enabled: bool, min_free: int = MEMORY_MIN_FREE_PCT) -> Dict[str, Any]:
    free, swap = memory_free_pct(), swap_used_mb()
    rec = {"free_pct_at_launch": free, "swap_used_mb_at_launch": swap, "threshold_pct": min_free,
           "guard": enabled}
    if enabled and (free is None or free < min_free):
        raise Refusal(f"memory guard: free {free} % < {min_free} %; wait (check at most every 10 min)")
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


UNIT_NAMES = ("U4", "U1e", "U1l", "U3", "U3a", "U3b", "U5a", "U5b", "U5a_PA_shape", "U5b_PA_shape", "GT3")


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
    while n_drawn < max_draws and (n_feas < target_feasible or n_drawn < min_draws) and \
            not (n_feas == 0 and n_drawn >= 12288):
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
    raw = "tie" if not (np.isfinite(d) and abs(d) > band) else ("a" if d > 0 else "b")
    reading = raw if resolved else "tie"
    return {"d": d, "band": band, "reading": reading, "reading_raw": raw, "resolved": bool(resolved)}


def contrast_credit(primary: Dict[str, Any], per_unit: Mapping[str, Dict[str, Any]], key: str = "reading") -> str:
    """Sec. 10.1's rule on the claim units (Holm's guard step is NOT simulated in the power check: optimistic):
    credit a iff the primary favours a and no claim unit favours b; symmetric for b; else none. ``key`` =
    "reading_raw" reads the units ignoring the density-resolution rule (a diagnostic)."""
    if primary[key] == "a" and not any(v[key] == "b" for v in per_unit.values()):
        return "a"
    if primary[key] == "b" and not any(v[key] == "a" for v in per_unit.values()):
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
        if rival != "R0":    # OPEN DIAGNOSTIC (not a scored unit): Table 2's argon PA + PAA, predicted vs measured
            f = ev.feasible
            out[rival]["t2_argon_total_mM_q10_50_90"] = _q(ev.diag["t2_total"][f, 0])
            out[rival]["t2_argon_total_measured_mM"] = round(sum(data.arp_t2["argon"]), 3)
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

    def gh(zpts):                                   # Gauss-Hermite in logit space, with the Jacobian
        uz = _sigmoid(zpts)
        upp = np.full((zpts.shape[0], len(cs)), 0.5)
        upp[:, i1], upp[:, i2] = uz[:, 0], uz[:, 1]
        e = L.evaluate(upp, "R3", L.PRIMARY, data, trapped=False)
        v = L.uarp_loglik(data, e.pred["UARP_c"], e.theta["cv_arp"], L.PRIMARY, measured=L.arp_interval_measured(data))
        jac = np.sum(-np.abs(zpts) - 2.0 * np.log1p(np.exp(-np.abs(zpts))), axis=1)
        return np.where(e.feasible, v + jac, -1e300)

    zg = _logit(np.column_stack([m1.ravel(), m2.ravel()]))
    zmean = wn.ravel() @ zg
    zcov = np.cov(zg.T, aweights=wn.ravel() + 1e-300) + 1e-3 * np.eye(2)
    q = pd.quadrature_lppd(gh, zmean, zcov, n_nodes=20)
    return {"lppd_trapezoid_41x41": round(lppd_grid, 4), "lppd_psis": round(est.lppd, 4), "k_hat": est.k_hat,
            "ess": est.ess, "mc_se": est.mc_se, "lppd_gauss_hermite_20_logit": round(q.lppd, 4),
            "psis_minus_grid": round(est.lppd - lppd_grid, 4), "gh_minus_grid": round(q.lppd - lppd_grid, 4),
            "pass": bool(abs(est.lppd - lppd_grid) < max(0.05, 3 * (est.mc_se or 0))),
            "gh_pass": bool(abs(q.lppd - lppd_grid) < 0.05)}


def psis_pilot(data: L.Data, rng: np.random.Generator, budget: Budget, evals: int = 300) -> Dict[str, Any]:
    """The one-unit PSIS pilot at full dimension (sec. 13 item 2): U3 under R2 on SYNTHETIC data (a truth drawn from
    R2's prior reweighted by the non-U3 fed units), conditioned on U4, U1e, U1l: the numerator integral
    int L(U4, U1e, U1l, U3) p by ``fit_cell`` (the decisive fits' estimator) from three starts (the prior draws with
    the largest weights). Records k-hat, ESS, wall time and the projection to the PSIS plan (385 cells)."""
    t_start = time.time()
    rival = "R2"
    d = len(L.coords(rival))
    obs = real_obs(data)
    prior = make_draw_set(rival, L.PRIMARY, data, 256, seed=901, budget=budget)
    lw = weights_from(prior, ("U4", "U1e", "U1l"), obs, L.PRIMARY, data)
    j = int(sample_truths(prior, lw, 1, rng)[0])
    truth_u3 = prior.points["U3"][j]
    u3 = synth_unit(L.u3_unit(10.0 ** truth_u3), truth_u3, rng)
    obs_syn = Obs({**obs.units, "U3": u3}, obs.uarp_y)
    units = ("U4", "U1e", "U1l", "U3")
    chol = {u: np.linalg.cholesky(L.unit_covariance(obs_syn.units[u])) for u in units if u != "U3"}
    ll_ref = float(np.max(L.unit_loglik(u3, prior.points["U3"])))
    n_ev = {"n": 0, "infeasible": 0}

    def evaluate_points(pts: np.ndarray):
        ev = L.evaluate(pts, rival, L.PRIMARY, data)
        n_ev["n"] += len(pts)
        n_ev["infeasible"] += int(np.sum(~ev.feasible))
        ll = np.zeros(len(pts))
        res = []
        for u in units:
            p_ = L.points_for(u, ev)
            l_ = L.unit_loglik(obs_syn.units[u], p_)
            ll = ll + l_
            if u == "U3":
                res.append(np.sqrt(np.maximum(-2.0 * (l_ - ll_ref), 0.0))[:, None])
            else:
                un = obs_syn.units[u]
                res.append(np.linalg.solve(chol[u], (un.y[None, :] - p_ @ un.design.T).T).T)
        r = np.concatenate(res, axis=1)
        ok = ev.feasible & ~ev.diag["failed_mask"]
        return np.where(ok, ll, -np.inf), np.where(ok[:, None] & np.isfinite(r), r, 25.0)

    w_all = lw + L.unit_loglik(u3, prior.points["U3"])
    starts = [prior.u[i] for i in np.argsort(-w_all)[:3]]
    cell = fit_cell(evaluate_points, d, starts, rng, evals=evals, budget=budget)
    est = cell["estimate"]
    wall = time.time() - t_start
    per_cell_h = wall / 3600.0
    return {"rival": rival, "dimension": d, "unit": "U3 | U4, U1e, U1l (synthetic U3)",
            "mode_search": [{"cost": f["cost"], "status": f["status"], "evaluations": f["evaluations"],
                             "converged": f.get("converged"), "polish_flagged": (f.get("polish") or {}).get("flagged")}
                            for f in cell["fits"]],
            "proposal": "Student-t 5 df in logit space at the best mode; Gauss-Newton covariance x 4, eigenvalues "
                        "clipped to [0.01, pi^2 / 3]",
            "psis": {"S": est.n_draws, "k_hat": est.k_hat, "ess": est.ess, "accepted": est.accepted,
                     "lppd": est.lppd, "mc_se": est.mc_se, "reason": est.reason,
                     "zero_weight_draws": est.n_nonfinite},
            "evaluations": n_ev["n"], "infeasible_evaluations": n_ev["infeasible"],
            "wall_s": round(wall, 1), "projection_385_cells_h_one_worker": round(385 * per_cell_h, 1),
            "exceeds_24h_cap": bool(385 * per_cell_h > 24.0)}


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
    mem = memory_guard(not args.no_memory_guard, args.memory_min_free)
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
    out["v_smc_grid"] = v_smc_grid(data, rng)
    out["v_smc_gaussian_32"] = v_smc_gaussian_32(rng)
    out["smc_pilot"] = {}
    for rival in ("R2", "R3", "R4"):
        try:
            out["smc_pilot"][rival] = smc_pilot(data, rival, rng, budget)
        except Refusal as e:
            out["smc_pilot"][rival] = {"refused": str(e)}
    out["v2_psis_pilot"] = ("superseded (v2.1 amendment): the v2 record at commit 5a6a30aa shows the Student-t PSIS "
                            "pilot unresolved (R2: k-hat 2.5, ESS 1.2 at S = 4096)")
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
           f"{vd['ess']:.0f}); Gauss-Hermite 20x20 in logit space {vd['lppd_gauss_hermite_20_logit']} (agrees: "
           f"{vd['gh_pass']})", ""]
    L_ += ["## Open diagnostic: Table 2's argon PA + PAA (not a scored unit)", "",
           "| rival | predicted q10 / q50 / q90 (mM) | measured (mM) |", "|---|---|---|"]
    for r, v in o["v_time"].items():
        if "t2_argon_total_mM_q10_50_90" in v:
            L_.append(f"| {r} | {v['t2_argon_total_mM_q10_50_90']} | {v['t2_argon_total_measured_mM']} |")
    vs, vg = o["v_smc_grid"], o["v_smc_gaussian_32"]
    L_ += ["", "## v2.1 estimator: adaptive-tempering SMC", "",
           f"- V-SMC (a), two-coordinate R3 toy against the grid: {'PASS' if vs['pass'] else 'FAIL'} (grid "
           f"{vs['lppd_grid']}, SMC {vs['smc']['lppd']} +- {vs['smc']['se']}, steps {vs['smc']['steps']}, max k-hat "
           f"{vs['smc']['max_k_hat']}, resolved {vs['smc']['resolved']})",
           f"- V-SMC (b), 32-dim Gaussian with a known answer: {'PASS' if vg['pass'] else 'FAIL'} (exact "
           f"{vg['exact_lppd']}, SMC {vg['smc']['lppd']:.4f} +- {vg['smc']['se_replicates']:.4f} over "
           f"{vg['replicates']} replicates {vg['smc']['replicates']}, max k-hat {vg['smc']['max_k_hat']:.3f}, min step "
           f"ESS {vg['smc']['min_ess']:.0f}; plain reweighting of 8 192 draws: {vg['reweighting_8192']})", "",
           "## SMC pilot at full dimension (U3 | U4, U1e, U1l, synthetic U3)", "",
           "| rival | coordinates | resolved | lppd +- se | steps | max k-hat | evaluations | wall s | reweighting "
           "(1 024) | projection 385 cells, h (1 / 4 workers) |", "|---|---|---|---|---|---|---|---|---|---|"]
    for r, v in o["smc_pilot"].items():
        if "smc" not in v:
            L_.append(f"| {r} | refused: {v.get('refused')} | | | | | | | | |")
            continue
        sm = v["smc"]
        L_.append(f"| {r} | {v['dimension']} | {sm['resolved']} ({'; '.join(sm['reasons']) or 'ok'}) | "
                  f"{sm['lppd']:.2f} +- {sm['se_replicates']:.2f} {sm['replicates']} | {sm['steps']} | "
                  f"{sm['max_k_hat']:.2f} | {sm['evaluations']} | {v['wall_s']} | "
                  f"k-hat {v['reweighting_1024']['k_hat']:.2f}, ESS {v['reweighting_1024']['ess']:.1f} | "
                  f"{v['projection_385_cells_h_one_worker']} / {v['projection_385_cells_h_four_workers']} |")
    L_ += ["", f"v2 pilot: {o['v2_psis_pilot']}", "",
           "## V-blind (blind-order) stub", "", "```", json.dumps(o["blind_order_stub"], indent=1), "```", "",
           "## Declared choices where the prereg is silent", ""] + [f"- {d}" for d in o["declared_choices"]]
    L_ += ["", f"Memory at launch: {o['memory']}; compute {o['compute_wall_h']} h wall."]
    return "\n".join(L_)


# ============================================================================================ the screen

SIG_BEARS = {"R1": ("stoich", "ARP", "cross", "early"), "R2": ("stoich", "lag", "fp", "cross", "early"),  # sigma-O2: the budget forms
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
    mem = memory_guard(not args.no_memory_guard, args.memory_min_free)
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


def synth_trapped(tset: DrawSet, j: int, data: L.Data, rng: np.random.Generator) -> Dict[str, L.Unit]:
    """The trapped pot's SYNTHETIC units from one truth (v2.1): U3 (full Sigma_u), and U3a / U3b built from the SAME
    synthetic trapped data; U3a's Table 2 argon PA is synthetic too (the truth's prediction with print error, s_rep,
    the vial's run offset and omega)."""
    f3 = tset.points["U3"][j]
    u3 = synth_unit(L.u3_unit(10.0 ** f3), f3, rng)
    y_lin = 10.0 ** u3.y
    y_lin = np.where(u3.censored, 0.5 * L.U3_CENSOR, y_lin)
    err = data.arp_t2_err["argon"][0]
    sd = math.sqrt(err ** 2 + L.S_REP_DEX ** 2 + L.SIGMA_RUN_DEX ** 2 + L.OMEGA_DEX ** 2)
    pa = 10.0 ** (tset.points["U4"][j][0] + sd * rng.standard_normal())
    return {"U3": u3, "U3a": L.u3a_unit(y_lin, pa, err), "U3b": L.u3b_unit(y_lin)}


TRAPPED_CHECKS = {"U3": ("U3",), "U3a": ("U3a",), "U3b": ("U3b",), "U3a+U3b": ("U3a", "U3b")}


def power_trapped(data: L.Data, sets: Mapping[Tuple[str, str], DrawSet], truths: Mapping[str, DrawSet], r: int,
                  rng: np.random.Generator, n_tie: int, budget: Budget) -> Dict[str, Any]:
    """Sec. 12.1 and the v2.1 amendment: the trapped pot's units, each ALONE by the route rule (R1 vs R2), and the
    {U3a, U3b} block (the amended route primary), on SYNTHETIC data at the design times. Truths from each rival's
    prior reweighted by the non-U3 fed units (U4, U1e, U1l on real data; Table 2 through the calibrations); scored
    under the primary and the Cu-law form (required agreement), conditioned on the real fed units."""
    obs = real_obs(data)
    fed = ("U4", "U1e", "U1l")
    readings: Dict[str, List[Dict[str, Dict[str, Any]]]] = {}
    for truth in ("R1", "R2"):
        tset = truths[truth]
        lw = weights_from(tset, fed, obs, L.PRIMARY, data)
        per = []
        for j in sample_truths(tset, lw, r, rng):
            budget.check("trapped-pot power")
            o = Obs({**obs.units, **synth_trapped(tset, j, data, rng)}, obs.uarp_y)
            rec: Dict[str, Dict[str, Any]] = {}
            for fname in ("primary", "cu_linear"):
                a, b = sets[("R1", fname)], sets[("R2", fname)]
                form = POWER_FORMS[fname]
                for check, units in TRAPPED_CHECKS.items():
                    v = (unit_verdict(a, b, units[0], fed, o, data, form, rng, n_tie) if len(units) == 1 else
                         block_verdict(a, b, units, fed, o, data, form, rng, n_tie))
                    rec[f"{check}|{fname}"] = v
            per.append(rec)
            if len(per) % 20 == 0:
                print(f"trapped power truth {truth}: {len(per)} / {r}, {budget.elapsed_h():.2f} h", flush=True)
        readings[truth] = per
    out: Dict[str, Any] = {}
    for check in TRAPPED_CHECKS:
        res = {}
        for forms, label in ((("primary",), "primary_only"), (("primary", "cu_linear"), "primary_and_cu_law")):
            for key, suffix in (("reading", ""), ("reading_raw", "_ignoring_resolution")):
                for credit, hit, false_t in (("a", "R1", "R2"), ("b", "R2", "R1")):
                    def credited(rec):
                        return all(rec[f"{check}|{f}"][key] == credit for f in forms)
                    res[f"{label}{suffix}|credit_{hit}"] = binomial_decision(
                        sum(credited(x) for x in readings[hit]), sum(credited(x) for x in readings[false_t]), r)
        prim = [f"{check}|primary"]
        out[check] = {"results": res,
                      "resolution_rate": {t: float(np.mean([x[prim[0]]["resolved"] for x in readings[t]]))
                                          for t in readings},
                      "d_quantiles": {t: _q([x[prim[0]]["d"] for x in readings[t]]) for t in readings},
                      "band_quantiles": {t: _q([x[prim[0]]["band"] for x in readings[t]]) for t in readings},
                      "decidable_under_the_rule": any(v["decidable"] for k, v in res.items()
                                                      if k.startswith("primary_and_cu_law|")),
                      "reaches_power_bar": any(v["decidable"] and not v["low_power"] for k, v in res.items()
                                               if k.startswith("primary_and_cu_law|"))}
    out["truth_reweighting"] = "U4, U1e, U1l on real data (Table 2 enters through the calibrations); never U3"
    return out


def power_pairs(data: L.Data, sets: Mapping[Tuple[str, str], DrawSet], truths: Mapping[str, DrawSet], r: int,
                rng: np.random.Generator, n_tie: int, budget: Budget,
                contrasts: Sequence[str] = tuple(CONTRASTS)) -> Dict[str, Any]:
    """Sec. 12.2: the three pairs. Synthetic claim units from a truth (truths reweighted by the real units outside the
    pair's claim set); the non-claim conditioning units real, except U3 (synthetic from the same truth: it is blind).
    Scored under the primary and the three forms with precomputed draw sets (o_Ar 0.1, Cu law, 20 mL headspace);
    the other forms are 'primary-only, optimistic'. Holm's guard step is not simulated (optimistic). The route
    contrast is read with the v2.1 primary (the {U3a, U3b} block) and, beside it, the literal v2 primary (U3)."""
    obs = real_obs(data)
    out: Dict[str, Any] = {}
    pair_sets = {"route": ("R1", "R2"), "transfer": ("R1", "R3"), "decay_law": ("R4", "R1drift")}
    reweight = {"route": ("U1l",), "transfer": ("U4", "U1e", "U1l"), "decay_law": ("U4", "U1e", "U1l")}
    for cname, spec in CONTRASTS.items():
        if cname not in contrasts:
            continue
        a_name, b_name = pair_sets[cname]
        readings: Dict[str, List[Dict[str, Dict[str, str]]]] = {}
        for truth in (a_name, b_name):
            tset = truths[truth]
            lw = weights_from(tset, reweight[cname], obs, L.PRIMARY, data)
            per = []
            for j in sample_truths(tset, lw, r, rng):
                budget.check(f"power {cname}")
                syn = dict(obs.units)
                syn.update(synth_trapped(tset, j, data, rng))
                for u in spec["claim"]:
                    if u not in ("U3", "U3a", "U3b", "UARP"):
                        syn[u] = synth_unit(obs.units[u], tset.points[u][j], rng)
                uarp_y = synth_uarp(data, tset.k_glc[j], rng) if "UARP" in spec["claim"] else obs.uarp_y
                os_ = Obs(syn, uarp_y)
                per.append({fname: _contrast_on(cname, spec, sets[(a_name, fname)], sets[(b_name, fname)], os_, data,
                                                form, rng, n_tie) for fname, form in POWER_FORMS.items()})
                if len(per) % 20 == 0:
                    print(f"power {cname} truth {truth}: {len(per)} / {r} sets, {budget.elapsed_h():.2f} h", flush=True)
            readings[truth] = per
        keys = list(readings[a_name][0]["primary"])
        res = {}
        for key in keys:
            suffix = {"rule": "", "raw": "|ignoring_resolution", "rule_literal": "|literal_U3_primary",
                      "raw_literal": "|literal_U3_primary|ignoring_resolution"}[key]
            for label, forms in (("primary_only", ("primary",)),
                                 ("four_forms_required_agreement", tuple(POWER_FORMS))):
                for credit, truth_hit, truth_false in (("a", a_name, b_name), ("b", b_name, a_name)):
                    hits = sum(all(rec[f][key] == credit for f in forms) for rec in readings[truth_hit])
                    false = sum(all(rec[f][key] == credit for f in forms) for rec in readings[truth_false])
                    res[f"{label}{suffix}|credit_{truth_hit}"] = binomial_decision(hits, false, r)

        def decided(pred):
            return any(v["decidable"] for k, v in res.items() if k.startswith("four_forms") and pred(k))

        strict = lambda k: "ignoring" not in k and "literal" not in k  # noqa: E731
        out[cname] = {"pair": [a_name, b_name], "claim_units": list(spec["claim"]), "truth_reweighting": reweight[cname],
                      "results": res, "demoted": not decided(strict),
                      "reaches_power_bar": any(v["decidable"] and not v["low_power"] for k, v in res.items()
                                               if k.startswith("four_forms") and strict(k)),
                      "decidable_ignoring_resolution": decided(lambda k: "ignoring" in k and "literal" not in k)}
        if cname == "route":
            out[cname]["decidable_literal_U3_primary"] = decided(lambda k: "literal" in k and "ignoring" not in k)
    return out


def _contrast_on(cname: str, spec, a: DrawSet, b: DrawSet, obs: Obs, data: L.Data, form: L.Form,
                 rng: np.random.Generator, n_tie: int) -> Dict[str, str]:
    """One synthetic set's contrast credit under one form: with the density-resolution rule ("rule") and ignoring
    it ("raw", a diagnostic); the route contrast also under its literal v2 primary (U3)."""

    def both(prim, per, tag=""):
        return {f"rule{tag}": contrast_credit(prim, per), f"raw{tag}": contrast_credit(prim, per, "reading_raw")}

    fed = ("U4", "U1e", "U1l")
    if cname == "route":
        cond = {"U4": (), "U1e": ("U4",), "U3": fed, "U3a": fed, "U3b": fed}
        per = {u: unit_verdict(a, b, u, cond[u], obs, data, form, rng, n_tie) for u in spec["claim"]}
        block = block_verdict(a, b, ("U3a", "U3b"), fed, obs, data, form, rng, n_tie)
        return {**both(block, per), **both(per["U3"], per, "_literal")}
    if cname == "decay_law":
        block = ("U3", "U5a", "U5b", "UARP")
        per = {}
        prev: Tuple[str, ...] = fed
        for u in block:
            per[u] = unit_verdict(a, b, u, prev, obs, data, form, rng, n_tie)
            prev = prev + (u,)
        return both(block_verdict(a, b, block, fed, obs, data, form, rng, n_tie), per)
    sides = []
    for cond0 in (fed + ("U3",), fed):
        per = {}
        prev = cond0
        for u in GLC_ORDER:
            per[u] = unit_verdict(a, b, u, prev, obs, data, form, rng, n_tie)
            prev = prev + (u,)
        sides.append(both(block_verdict(a, b, GLC_ORDER, cond0, obs, data, form, rng, n_tie), per))
    return {k: (sides[0][k] if sides[0][k] == sides[1][k] else "none") for k in ("rule", "raw")}


POWER_RIVALS = {"R1": L.PRIMARY, "R2": L.PRIMARY, "R3": L.PRIMARY, "R4": L.PRIMARY,
                "R1drift": L.Form(ph_drift_trapped=True)}


def power_set_specs(args, want: Sequence[str]) -> List[Dict[str, Any]]:
    """Every draw set the power parts need, with a stable seed (identical across processes: the cache key)."""
    specs = []
    for name in want:
        i = list(POWER_RIVALS).index(name)
        base = POWER_RIVALS[name]
        rival = "R1" if name == "R1drift" else name
        for fname, form in POWER_FORMS.items():
            k = list(POWER_FORMS).index(fname)
            specs.append({"key": (name, fname), "rival": rival,
                          "form": replace(form, ph_drift_trapped=base.ph_drift_trapped),
                          "target": args.score_feasible, "seed": args.seed + 1000 * i + 10 * k, "keep_k": False})
        specs.append({"key": (name, "truth"), "rival": rival,
                      "form": replace(L.PRIMARY, ph_drift_trapped=base.ph_drift_trapped),
                      "target": args.truth_feasible, "seed": args.seed + 7000 + 100 * i, "keep_k": True})
    return specs


def run_power(args) -> int:
    budget = Budget(args.cap_hours or 6.0)
    mem = memory_guard(not args.no_memory_guard, args.memory_min_free)
    data = L.load_data()
    t0 = time.time()
    part = args.power_part
    if part == "merge":
        return merge_power(args, mem)
    need = {"route": ("R1", "R2"), "transfer": ("R1", "R3"), "decay_law": ("R4", "R1drift")}
    contrasts = tuple(c for c in args.contrasts.split(",") if c)
    if part in ("u3", "trapped"):
        want: Tuple[str, ...] = ("R1", "R2")
    elif part == "sets":
        want = tuple(POWER_RIVALS)
    else:
        want = tuple(dict.fromkeys((["R1", "R2"] if part == "all" else []) + [r for c in contrasts for r in need[c]]))
    specs = power_set_specs(args, want)
    if part == "sets":
        i_shard, n_shard = (int(x) for x in args.shard.split("/"))
        specs = specs[i_shard::n_shard]
    sets: Dict[Tuple[str, str], DrawSet] = {}
    truths: Dict[str, DrawSet] = {}
    set_log = {}
    for sp in specs:
        ds = cached_draw_set(args.cache_dir, sp["rival"], sp["form"], data, sp["target"], seed=sp["seed"],
                             max_draws=args.max_draws, keep_k=sp["keep_k"], budget=budget)
        name, fname = sp["key"]
        if fname == "truth":
            truths[name] = ds
        else:
            sets[(name, fname)] = ds
        set_log[f"{name}|{fname}"] = {"feasible": ds.n, "drawn": ds.n_drawn, "failed": ds.n_failed,
                                      "wall_s": ds.wall_s, "infeasible": ds.reasons}
        print("set", name, fname, ds.n, "/", ds.n_drawn, "failed", ds.n_failed, f"{ds.wall_s}s", flush=True)
    if part == "sets":
        print(f"sets shard {args.shard} done in {time.time() - t0:.0f} s", flush=True)
        return 0
    rng = np.random.default_rng(args.seed + 1 + 101 * ("trapped", "pairs", "all", "u3").index(part)
                                + sum(ord(c) for c in ",".join(contrasts)))
    failed = sum(v["failed"] for v in set_log.values())
    out: Dict[str, Any] = {"memory": mem, "R": args.reps, "tie_redraws": args.tie_redraws, "draw_sets": set_log,
                           "failed_runs_total": failed, "sets_wall_s": round(time.time() - t0, 1)}
    if failed:
        out["refused"] = f"{failed} failed runs in the draw sets: every power result is refused"
    else:
        try:
            if part in ("u3", "trapped", "all"):
                out["trapped_12_1_and_v2_1"] = power_trapped(data, sets, truths, args.reps, rng, args.tie_redraws,
                                                             budget)
            if part in ("pairs", "all") and contrasts:
                out["pairs_12_2"] = power_pairs(data, sets, truths, args.reps, rng, args.tie_redraws, budget,
                                                contrasts=contrasts)
        except Refusal as e:
            out["compute_capped"] = str(e)
    out["compute_wall_h"] = round(budget.elapsed_h(), 3)
    out["approximations"] = POWER_APPROXIMATIONS(args)
    if args.record_suffix:      # a part: written to the cache dir (scratch) and merged into the record later
        part_dir = (args.cache_dir or args.out_dir) / "parts"
        part_dir.mkdir(parents=True, exist_ok=True)
        (part_dir / f"power_{args.record_suffix}.json").write_text(
            provenance.dumps_compact(out, default=_jsonable) + "\n")
        print(f"wrote part {args.record_suffix}", flush=True)
        return 0
    _record("power", out, render_power(out), args.out_dir, layer_inputs())
    return 0


def POWER_APPROXIMATIONS(args) -> List[str]:  # noqa: N802
    return [
        "draw sets are scrambled-Sobol prior draws conditioned on the calibrations (feasible only), reweighted per "
        "synthetic set (sec. 12.2's 'by reweighting'); a density whose PSIS k-hat >= 0.7 or ESS < 100 reads 'not "
        "resolved' (a tie, unchanged rule); results ignoring that rule are reported beside",
        f"draw sets enlarged to {args.score_feasible} feasible draws per rival and form ({args.truth_feasible} for the "
        "truths) under the v2.1 amendment; the SMC fallback (v2.1) is NOT run per synthetic set (its measured cost: "
        "checks record), so an unresolved density stays a tie here: power is conservative on that account",
        f"tie bands from {args.tie_redraws} measurement-noise redraws per unit and set (the prereg's 4 000 for scoring)",
        "censored U3 values: product of conditional marginals (exact for one); adequacy leaves them out",
        "Holm's guard step is not simulated: power is optimistic on that account",
        "synthetic sets are generated from the truth's primary-form predictions and scored under each form",
        "a rival with no feasible draw under a form (its calibration infeasible there) has lppd -inf: the other side is "
        "read as favoured under that form without its adequacy test (flagged by_infeasibility; optimistic)",
        "U3a / U3b are built from the same synthetic trapped data as U3; U3a's Table 2 argon PA is synthetic from the "
        "truth (its print error, s_rep, run offset and omega)",
    ]


def merge_power(args, mem) -> int:
    """Merge the power parts (written by --record-suffix runs into <cache>/parts) into the one record; refuses if a
    part refused or lost runs."""
    part_dir = (args.cache_dir or args.out_dir) / "parts"
    parts = sorted(part_dir.glob("power_*.json"))
    if not parts:
        raise Refusal(f"no power parts under {part_dir}")
    out: Dict[str, Any] = {"memory": mem, "parts": [p.name for p in parts], "draw_sets": {}, "failed_runs_total": 0,
                           "compute_wall_h_parts": {}}
    for p in parts:
        d = json.loads(p.read_text())
        out["R"], out["tie_redraws"] = d["R"], d["tie_redraws"]
        out["draw_sets"].update(d["draw_sets"])
        out["failed_runs_total"] += d["failed_runs_total"]
        out["compute_wall_h_parts"][p.stem] = d["compute_wall_h"]
        out["part_memory"] = {**out.get("part_memory", {}), p.stem: d["memory"]}
        for k in ("refused", "compute_capped"):
            if d.get(k):
                out[k] = f"{out.get(k, '')} {p.stem}: {d[k]}".strip()
        if "trapped_12_1_and_v2_1" in d:
            out["trapped_12_1_and_v2_1"] = d["trapped_12_1_and_v2_1"]
        for c, v in (d.get("pairs_12_2") or {}).items():
            out.setdefault("pairs_12_2", {})[c] = v
        out["approximations"] = d["approximations"]
    if out["failed_runs_total"]:
        out["refused"] = f"{out['failed_runs_total']} failed runs in the parts: every power result is refused"
    out["compute_wall_h"] = round(sum(out["compute_wall_h_parts"].values()), 3)
    _record("power", out, render_power(out), args.out_dir, layer_inputs())
    return 0


def render_power(o: Dict[str, Any]) -> str:
    L_ = ["# M3 phase 1: power checks (v2.1; sec. 12.1 trapped pot, synthetic; sec. 12.2 three pairs)", "",
          f"*Generated by `scripts/generators/fit_m3_strecker_layer.py --power`. R = {o['R']} synthetic sets per truth; "
          "binomial demotion rule: one-sided Fisher exact test of power against false credit at 0.05; power < 0.5 is "
          "'low power'. The trapped pot on SYNTHETIC data at 20 / 40 / 80 / 160 / 320 min (the keyed Figs. 2 - 3 were "
          "not opened). Design: prereg v2 + the v2.1 amendment (sec. 22). SINGLE LABORATORY.*", ""]
    if o.get("refused"):
        L_ += [f"**REFUSED: {o['refused']}**", ""]
    if o.get("compute_capped"):
        L_ += [f"**COMPUTE-CAPPED: {o['compute_capped']}**", ""]
    L_ += ["## Draw sets", "", "| set | feasible / drawn | failed | wall s |", "|---|---|---|---|"]
    for k, v in o["draw_sets"].items():
        L_.append(f"| {k} | {v['feasible']} / {v['drawn']} | {v['failed']} | {v['wall_s']} |")
    hdr = ["| scoring and credit | power [95 % CI] | false credit [95 % CI] | Fisher p | decidable | low power |",
           "|---|---|---|---|---|---|"]
    if "trapped_12_1_and_v2_1" in o:
        L_ += ["", "## The trapped pot's units alone and the {U3a, U3b} block (route rule, R1 vs R2)"]
        for check, v in o["trapped_12_1_and_v2_1"].items():
            if not isinstance(v, dict) or "results" not in v:
                continue
            L_ += ["", f"### {check}: decidable under the rule (primary and Cu law) {v['decidable_under_the_rule']}; "
                   f"reaches the power bar {v['reaches_power_bar']}", "",
                   f"resolution rate {v['resolution_rate']}; d quantiles {v['d_quantiles']}; tie band quantiles "
                   f"{v['band_quantiles']}", ""] + hdr
            for k, r in v["results"].items():
                L_.append(f"| {k} | {r['power']:.3f} {r['power_ci95']} | {r['false_credit']:.3f} "
                          f"{r['false_credit_ci95']} | {r['fisher_p_one_sided']:.3g} | {r['decidable']} | "
                          f"{r['low_power']} |")
    if "pairs_12_2" in o:
        L_ += ["", "## Sec. 12.2: the three pairs"]
        for c, v in o["pairs_12_2"].items():
            L_ += ["", f"### {c} ({' vs '.join(v['pair'])}): demoted {v['demoted']}; reaches the power bar "
                   f"{v['reaches_power_bar']}; decidable ignoring the resolution rule "
                   f"{v['decidable_ignoring_resolution']}"
                   + (f"; decidable under the literal U3 primary {v['decidable_literal_U3_primary']}"
                      if "decidable_literal_U3_primary" in v else ""), ""] + hdr
            for k, r in v["results"].items():
                L_.append(f"| {k} | {r['power']:.3f} {r['power_ci95']} | {r['false_credit']:.3f} "
                          f"{r['false_credit_ci95']} | {r['fisher_p_one_sided']:.3g} | {r['decidable']} | "
                          f"{r['low_power']} |")
    L_ += ["", "## Approximations (declared)", ""] + [f"- {a}" for a in o["approximations"]]
    L_ += ["", f"Memory at launch: {o['memory']}; compute {o['compute_wall_h']} h wall (sum over parts)."]
    return "\n".join(L_)


# =========================================================================== the decisive fits and scoring

def _logit(u):
    u = np.clip(u, 1e-9, 1 - 1e-9)
    return np.log(u) - np.log1p(-u)


def _sigmoid(z):
    return 1.0 / (1.0 + np.exp(-z))


def _clip_eigen(cov: np.ndarray, lo: float, hi: float) -> np.ndarray:
    w, v = np.linalg.eigh(0.5 * (cov + cov.T))
    w = np.clip(np.where(np.isfinite(w), w, hi), lo, hi)
    return (v * w) @ v.T


LOGIT_PRIOR_VAR = math.pi ** 2 / 3.0        # a uniform coordinate's variance in logit space


def fit_cell(evaluate_points: Callable[[np.ndarray], Tuple[np.ndarray, np.ndarray]], d: int, starts: Sequence[np.ndarray],
             rng: np.random.Generator, evals: int = 300, s_schedule: Sequence[int] = (256, 1024, 4096),
             spread_z: Optional[np.ndarray] = None, widen: float = 4.0,
             budget: Optional[Budget] = None) -> Dict[str, Any]:
    """One (unit set, form, rival) cell of the decisive fits (sec. 9.1). The mode of log(L x p) on the unit cube from
    the three starts (fitting.fit, central polish, require_converged='flag'; each start <= ``evals`` evaluations,
    the polish included); then a Student-t (5 df) proposal at the best mode IN LOGIT SPACE (every draw stays inside
    the bands; the density is carried back with the Jacobian), scaled by the Gauss-Newton covariance widened x
    ``widen`` with its eigenvalues clipped to [0.01, the uniform prior's logit variance] (a point curvature can be
    too narrow at a kink: lesson 2026-11-05; PSIS's k-hat is the check), or by ``spread_z``. S runs through
    ``s_schedule`` until PSIS accepts. ``evaluate_points`` maps (n, d) unit-cube points to (log L per draw, whitened
    residuals (n, k)): the residuals drive the mode search, log L the density (log p = 0 on the cube)."""
    from src.kinetic_core import fitting

    def res_many(pts):
        _, r = evaluate_points(np.clip(np.asarray(pts, float), 0.0, 1.0))
        return [row for row in r]

    k_rows = evaluate_points(np.full((1, d), 0.5))[1].shape[1]
    prob = fitting.Problem(keys=tuple(f"u{i}" for i in range(d)), kinds=("acid_yield",) * d, lower=np.zeros(d),
                           upper=np.ones(d), rows=tuple(range(k_rows)),
                           residuals=lambda th, rows: res_many(th[None, :])[0],
                           evaluate_many=lambda pts, rows: res_many(pts))
    polish = 2 * d + 3
    fits = []
    for s0 in starts:
        if budget is not None:
            budget.check("mode search")
        fits.append(fitting.fit(prob, s0, budget=max(evals - polish, d + 2), central_polish=1, polish_budget=polish,
                                require_converged="flag"))
    best = min(fits, key=lambda f: f.cost)
    u0 = np.clip(best.theta, 1e-6, 1 - 1e-6)
    z0 = _logit(u0)
    if spread_z is None:
        jac = fitting.central_jacobian(prob, u0, np.full(d, 1e-9), np.full(d, 1 - 1e-9)) * (u0 * (1 - u0))[None, :]
        spread_z = _clip_eigen(widen * np.linalg.pinv(jac.T @ jac), 1e-2, LOGIT_PRIOR_VAR)
    xs, lqs, lls = [], [], []
    est = None
    for s_n in s_schedule:
        if budget is not None:
            budget.check("PSIS sampling")
        need = s_n - sum(len(x) for x in xs)
        z, lqz = pd.multivariate_t_draws(z0, spread_z, need, rng)
        u = np.clip(_sigmoid(z), 1e-12, 1 - 1e-12)
        log_jac = -np.abs(z) - 2.0 * np.log1p(np.exp(-np.abs(z)))      # log u (1 - u), stable for any z
        lq = lqz - np.sum(log_jac, axis=1)                             # q in u-space: q_z / |du/dz|
        ll, _ = evaluate_points(u)
        xs.append(u)
        lqs.append(lq)
        lls.append(np.where(np.isfinite(ll), ll, -np.inf))
        est = pd.importance_lppd(np.concatenate(lls), np.concatenate(lqs))
        if est.accepted:
            break
    return {"theta_mode": best.theta, "fits": [f.as_dict(prob.keys) for f in fits], "draws": np.concatenate(xs),
            "log_q": np.concatenate(lqs), "log_lik": np.concatenate(lls), "estimate": est,
            "not_resolved": not est.accepted, "evaluations": sum(f.evaluations for f in fits)}


# ===================================================================== the v2.1 estimator: adaptive-tempering SMC

@dataclass
class SMCResult:
    lppd: float
    resolved: bool
    reason: str
    steps: int
    evaluations: int
    max_k_hat: float
    min_ess: float
    final_acceptance: float
    se: float
    betas: List[float]
    particles_u: np.ndarray
    log_cond: np.ndarray
    log_unit: np.ndarray
    extra: Dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, Any]:
        return {"lppd": self.lppd, "resolved": self.resolved, "reason": self.reason, "steps": self.steps,
                "evaluations": self.evaluations, "max_k_hat": self.max_k_hat, "min_ess": self.min_ess,
                "final_acceptance": self.final_acceptance, "se": self.se, **self.extra}


def _systematic(w: np.ndarray, n: int, rng: np.random.Generator) -> np.ndarray:
    c = np.cumsum(w / w.sum())
    c[-1] = 1.0
    return np.searchsorted(c, (rng.random() + np.arange(n)) / n)


def _log_jac_z(z: np.ndarray) -> np.ndarray:
    return np.sum(-np.abs(z) - 2.0 * np.log1p(np.exp(-np.abs(z))), axis=1)


def smc_lppd(eval_fn: Callable[[np.ndarray], Tuple[np.ndarray, np.ndarray]], u_init: np.ndarray,
             rng: np.random.Generator, n_moves: int = 2, ess_frac: float = 0.5, max_steps: int = 60,
             k_max: float = pd.K_HAT_MAX, ess_min: float = pd.ESS_MIN, acc_min: float = 0.05,
             budget: Optional[Budget] = None) -> SMCResult:
    """lppd(u | cond) = log E_{p(theta | cond)} L_u(theta) by adaptive-tempering sequential Monte Carlo (Del Moral,
    Doucet & Jasra 2006; Chopin 2002): particles from p(theta | cond) (``u_init``, unit-cube points), the targets
    p(theta | cond) L_u^beta with beta 0 -> 1 chosen so that each step's incremental-weight ESS is ``ess_frac`` N,
    systematic resampling, and ``n_moves`` independence Metropolis-Hastings moves per step in LOGIT space (a
    Student-t fitted to the particles, covariance x 1.5; the uniform band prior carried by the Jacobian).
    log Z = sum_t log mean_i W_i exp(d_beta l_i). ``eval_fn`` maps (n, d) unit-cube points to (log L_cond,
    log L_u), -inf for an infeasible draw. Resolved iff every step's incremental weights have PSIS k-hat < k_max,
    the particle ESS stays >= max(ess_min, ...) and the last moves' acceptance >= acc_min; else "density not
    resolved" (a tie, unchanged rule). The SE is the delta-method sum over steps (independent-steps approximation,
    declared)."""
    u = np.clip(np.asarray(u_init, float), 1e-9, 1 - 1e-9)
    n, d = u.shape
    lc, lu = eval_fn(u)
    n_eval = n
    beta, log_z, var_z = 0.0, 0.0, 0.0
    betas, khats, esss = [0.0], [], []
    scale = 1.5                      # the independence proposal's covariance inflation
    acc = float("nan")
    reason = ""

    def moves(u, lc, lu, beta, scale):
        """Independence Metropolis-Hastings: a Student-t (5 df) fitted to the particle cloud in logit space (its
        covariance x ``scale``); one accepted move is a fresh draw, so mixing does not degrade with d as a random
        walk's does (South, Pettitt & Drovandi 2019)."""
        from scipy.stats import multivariate_t
        nonlocal n_eval
        z = np.log(u) - np.log1p(-u)
        mu = z.mean(axis=0)
        cov = _clip_eigen(np.cov(z.T) * scale + 1e-6 * np.eye(d), 1e-6, LOGIT_PRIOR_VAR * 4)
        q = multivariate_t(loc=mu, shape=cov, df=5)
        lq_z = np.atleast_1d(q.logpdf(z))
        accs = []
        for _ in range(n_moves):
            zp = np.atleast_2d(q.rvs(size=n, random_state=rng)).reshape(n, d)
            lq_p = np.atleast_1d(q.logpdf(zp))
            up = np.clip(1.0 / (1.0 + np.exp(-zp)), 1e-12, 1 - 1e-12)
            lcp, lup = eval_fn(up)
            n_eval += n
            cur = lc + beta * lu + _log_jac_z(z) - lq_z
            new = lcp + beta * lup + _log_jac_z(zp) - lq_p
            with np.errstate(invalid="ignore"):
                acc_ = np.isfinite(new) & (np.log(rng.random(n)) < new - cur)
            z, u = np.where(acc_[:, None], zp, z), np.where(acc_[:, None], up, u)
            lc, lu = np.where(acc_, lcp, lc), np.where(acc_, lup, lu)
            lq_z = np.where(acc_, lq_p, lq_z)
            accs.append(float(acc_.mean()))
        return u, lc, lu, scale, float(np.mean(accs))

    if not np.all(np.isfinite(lc) & np.isfinite(lu)):
        keep = np.isfinite(lc) & np.isfinite(lu)
        if not keep.any():
            return SMCResult(-math.inf, False, "no finite initial particle", 0, n_eval, math.inf, 0.0, 0.0,
                             math.inf, betas, u, lc, lu)
        idx = rng.choice(np.flatnonzero(keep), size=n, replace=True)
        u, lc, lu = u[idx], lc[idx], lu[idx]
    u, lc, lu, scale, acc = moves(u, lc, lu, 0.0, scale)     # rejuvenate the (resampled) start at beta = 0
    for step in range(max_steps):
        if budget is not None:
            budget.check("SMC")

        def ess_at(b):
            w = np.exp((b - beta) * (lu - np.max(lu)))
            return w.sum() ** 2 / (w ** 2).sum()

        if ess_at(1.0) >= ess_frac * n:
            nb = 1.0
        else:
            lo, hi = beta, 1.0
            for _ in range(50):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if ess_at(mid) >= ess_frac * n else (lo, mid)
            nb = max(lo, beta + 1e-6)
        inc = (nb - beta) * lu
        m = float(np.max(inc))
        w = np.exp(inc - m)
        mean_w = float(w.mean())
        log_z += m + math.log(mean_w)
        var_z += float(np.var(w, ddof=1) / (n * mean_w ** 2))
        _, k = pd.psis(inc)
        khats.append(float(k))
        esss.append(float(w.sum() ** 2 / (w ** 2).sum()))
        idx = _systematic(w, n, rng)
        u, lc, lu = u[idx], lc[idx], lu[idx]
        beta = nb
        betas.append(beta)
        u, lc, lu, scale, acc = moves(u, lc, lu, beta, scale)
        if beta >= 1.0:
            break
    else:
        reason = f"beta reached only {beta:.3g} in {max_steps} steps"
    max_k = max(khats) if khats else math.inf
    min_ess = min(esss) if esss else 0.0
    reasons = [reason] if reason else []
    if not max_k < k_max:
        reasons.append(f"max step k-hat {max_k:.3g} >= {k_max}")
    if min_ess < ess_min:
        reasons.append(f"step ESS {min_ess:.1f} < {ess_min}")
    if not acc >= acc_min:
        reasons.append(f"final acceptance {acc:.3f} < {acc_min}")
    return SMCResult(log_z, not reasons, "; ".join(reasons), len(khats), n_eval, max_k, min_ess, acc,
                     math.sqrt(var_z), betas, u, lc, lu)


def smc_replicated(eval_fn, start_fn: Callable[[np.random.Generator], np.ndarray], rng: np.random.Generator,
                   k: int = 2, budget: Optional[Budget] = None, **kw) -> Dict[str, Any]:
    """``k`` independent SMC runs (each from its own SIR start): lppd = log of the mean of their Z (Z is unbiased,
    log Z is not), the stated error from the replicates' spread (the single-run delta-method SE is optimistic by
    ~1.5x at 32 dimensions: V-SMC (b)); resolved iff every replicate is."""
    runs = [smc_lppd(eval_fn, start_fn(rng), rng, budget=budget, **kw) for _ in range(k)]
    lz = np.array([r.lppd for r in runs])
    m = float(np.max(lz))
    lppd = m + math.log(float(np.mean(np.exp(lz - m)))) if np.isfinite(m) else -math.inf
    se = float(np.std(lz, ddof=1) / math.sqrt(k)) if k > 1 else runs[0].se
    return {"lppd": lppd, "se_replicates": se, "replicates": [round(float(x), 4) for x in lz],
            "resolved": all(r.resolved for r in runs), "reasons": [r.reason for r in runs if r.reason],
            "steps": [r.steps for r in runs], "evaluations": int(sum(r.evaluations for r in runs)),
            "max_k_hat": float(max(r.max_k_hat for r in runs)), "min_ess": float(min(r.min_ess for r in runs)),
            "final_acceptance": [round(r.final_acceptance, 3) for r in runs], "runs": runs}


def smc_sir_start(draws_u: np.ndarray, log_cond: np.ndarray, n: int, rng: np.random.Generator) -> np.ndarray:
    """Particles from p(theta | cond) by sampling-importance-resampling a prior draw set by its conditioning weights
    (the fed units' ESS is ~40 - 60 % of the set: checks record); the SMC's beta = 0 moves then rejuvenate them."""
    w = np.exp(log_cond - np.max(log_cond[np.isfinite(log_cond)]))
    w = np.where(np.isfinite(w), w, 0.0)
    return draws_u[_systematic(w, n, rng)]


def gaussian_toy_lppd(y: np.ndarray, s: float, m0: float, s0: float) -> float:
    """The 32-dim validation's exact answer: u_k ~ U[0, 1], cond N(u_k; m0, s0^2) (as a likelihood), unit
    N(y_k; u_k, s^2): lppd = sum_k log [int_0^1 N(y_k; u, s^2) N(u; m0, s0^2) du / int_0^1 N(u; m0, s0^2) du]."""
    from scipy.stats import norm
    out = 0.0
    for yk in np.atleast_1d(y):
        v = 1.0 / (1.0 / s ** 2 + 1.0 / s0 ** 2)
        mu = v * (yk / s ** 2 + m0 / s0 ** 2)
        lnorm = norm.logpdf(yk, m0, math.sqrt(s ** 2 + s0 ** 2))
        mass = norm.cdf(1.0, mu, math.sqrt(v)) - norm.cdf(0.0, mu, math.sqrt(v))
        den = norm.cdf(1.0, m0, s0) - norm.cdf(0.0, m0, s0)
        out += lnorm + math.log(mass) - math.log(den)
    return float(out)


def v_smc_gaussian_32(rng: np.random.Generator, d: int = 32, n: int = 512, k: int = 4) -> Dict[str, Any]:
    """V-SMC (b): a synthetic 32-dim case with a known answer (``gaussian_toy_lppd``). Plain reweighting of 8 192
    conditioned draws (the power check's estimator) is shown beside the replicated SMC (k runs of N particles, each
    from SIR of those draws). Pass: every replicate resolved and |lppd - exact| <= 3 x the replicates' SE."""
    m0, s0, s = 0.5, 0.25, 0.03
    y = rng.uniform(0.25, 0.75, d)
    exact = gaussian_toy_lppd(y, s, m0, s0)

    def eval_fn(u):
        lc = np.sum(-0.5 * ((u - m0) / s0) ** 2, axis=1)                 # unnormalised: cancels in the ratio
        lu = np.sum(-0.5 * ((y - u) / s) ** 2 - math.log(s * math.sqrt(2 * math.pi)), axis=1)
        return lc, lu

    prior = rng.random((8192, d))
    lc0, lu0 = eval_fn(prior)
    rw = conditional_lppd(lc0, lu0)
    res = smc_replicated(eval_fn, lambda r: smc_sir_start(prior, lc0, n, r), rng, k=k)
    err = res["lppd"] - exact
    out = {k_: v for k_, v in res.items() if k_ != "runs"}
    return {"dimension": d, "particles": n, "replicates": k, "exact_lppd": round(exact, 4), "smc": out,
            "smc_minus_exact": round(err, 4), "smc_within_3se": bool(abs(err) <= 3 * res["se_replicates"]),
            "reweighting_8192": {"lppd": round(rw["lppd"], 4), "k_hat": rw["k_hat"], "ess": rw["ess"],
                                 "resolved": rw["resolved"]},
            "pass": bool(res["resolved"] and abs(err) <= 3 * res["se_replicates"])}


def v_smc_grid(data: L.Data, rng: np.random.Generator, n: int = 256) -> Dict[str, Any]:
    """V-SMC (a): the V-density R3 toy (two coordinates, U-ARPglc) by SMC from the uniform prior, against the
    41 x 41 trapezoid grid."""
    cs = L.coords("R3")
    names = [c.name for c in cs]
    i1, i2 = names.index("log10_kd7"), names.index("log10_kd7_glc")
    meas = L.arp_interval_measured(data)

    def eval_fn(u2):
        up = np.full((u2.shape[0], len(cs)), 0.5)
        up[:, i1], up[:, i2] = u2[:, 0], u2[:, 1]
        e = L.evaluate(up, "R3", L.PRIMARY, data, trapped=False)
        v = L.uarp_loglik(data, e.pred["UARP_c"], e.theta["cv_arp"], L.PRIMARY, measured=meas)
        return np.where(e.feasible, 0.0, -np.inf), np.where(e.feasible, v, -np.inf)

    g = 41
    grid = np.linspace(0.0, 1.0, g)
    m1, m2 = np.meshgrid(grid, grid, indexing="ij")
    lc_g, lu_g = eval_fn(np.column_stack([m1.ravel(), m2.ravel()]))
    ll = np.where(np.isfinite(lc_g), lu_g, -np.inf).reshape(g, g)
    lmax = np.max(ll[np.isfinite(ll)])
    lppd_grid = math.log(np.trapezoid(np.trapezoid(np.exp(ll - lmax), grid, axis=1), grid)) + lmax
    t = time.time()
    res = smc_lppd(eval_fn, rng.random((n, 2)), rng)
    return {"lppd_grid": round(lppd_grid, 4), "smc": {k: (round(v, 4) if isinstance(v, float) else v)
                                                      for k, v in res.as_dict().items()},
            "smc_minus_grid": round(res.lppd - lppd_grid, 4), "wall_s": round(time.time() - t, 1),
            "pass": bool(res.resolved and abs(res.lppd - lppd_grid) <= max(0.05, 3 * res.se))}


def smc_pilot(data: L.Data, rival: str, rng: np.random.Generator, budget: Budget, n: int = 256) -> Dict[str, Any]:
    """The v2.1 PSIS-pilot replacement: U3 | U4, U1e, U1l under ``rival`` at full dimension on SYNTHETIC U3 (a truth
    from the rival's fed-reweighted prior), by SMC from SIR particles. Records resolution, steps, evaluations, wall
    time, and the projection to the decisive plan's cells."""
    t0 = time.time()
    obs = real_obs(data)
    prior = make_draw_set(rival, L.PRIMARY, data, 1024, seed=903 + RIVAL_INDEX[rival], budget=budget)
    lw = weights_from(prior, ("U4", "U1e", "U1l"), obs, L.PRIMARY, data)
    j = int(sample_truths(prior, lw, 1, rng)[0])
    u3 = synth_unit(L.u3_unit(10.0 ** prior.points["U3"][j]), prior.points["U3"][j], rng)
    obs_syn = Obs({**obs.units, "U3": u3}, obs.uarp_y)
    rw = conditional_lppd(lw, L.unit_loglik(u3, prior.points["U3"]))

    def eval_fn(u):
        ev = L.evaluate(u, rival, L.PRIMARY, data)
        lc = np.zeros(len(u))
        for un in ("U4", "U1e", "U1l"):
            lc = lc + L.unit_loglik(obs_syn.units[un], L.points_for(un, ev))
        lu = L.unit_loglik(u3, L.points_for("U3", ev))
        ok = ev.feasible & ~ev.diag["failed_mask"]
        return np.where(ok, lc, -np.inf), np.where(ok, lu, -np.inf)

    res = smc_replicated(eval_fn, lambda r: smc_sir_start(prior.u, lw, n, r), rng, k=2, budget=budget)
    res = {k: v for k, v in res.items() if k != "runs"}
    wall = time.time() - t0
    return {"rival": rival, "dimension": len(L.coords(rival)), "particles": n,
            "unit": "U3 | U4, U1e, U1l (synthetic U3)", "reweighting_1024": {k: rw[k] for k in
                                                                            ("lppd", "k_hat", "ess", "resolved")},
            "smc": res, "wall_s": round(wall, 1),
            "wall_per_density_s": round(wall, 1),
            "projection_385_cells_h_one_worker": round(385 * wall / 3600.0, 1),
            "projection_385_cells_h_four_workers": round(385 * wall / 3600.0 / 4.0, 1)}


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
    mem = memory_guard(not args.no_memory_guard, args.memory_min_free)
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
    ap.add_argument("--power-part", choices=("u3", "trapped", "pairs", "all", "sets", "merge"), default="all")
    ap.add_argument("--shard", default="0/1", help="power --power-part sets: build shard i of k")
    ap.add_argument("--contrasts", default="route,transfer,decay_law", help="power: the 12.2 pairs to run")
    ap.add_argument("--record-suffix", default="", help="power: a part's record name suffix (merged later)")
    ap.add_argument("--cache-dir", type=Path, default=None, help="draw-set cache (scratch; never under data/)")
    ap.add_argument("--max-draws", type=int, default=12288, help="power: draws per set at most")
    ap.add_argument("--cap-hours", type=float, default=None)
    ap.add_argument("--no-memory-guard", action="store_true")
    ap.add_argument("--memory-min-free", type=int, default=MEMORY_MIN_FREE_PCT,
                    help="memory guard threshold, % free (coordinator 2026-11-05: 30 with one worker)")
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
