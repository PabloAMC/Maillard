"""Probe: evaluations per carried rival on SUG-01's two arms.

Env: PROBE_BUDGET, PROBE_WORKERS, PROBE_PERTURB, PROBE_ONLY (substring of 'arm|point'),
PROBE_MODE: baseline | coldph (history-free pH solve) | bdf | radau | rtol8 (rtol 1e-8).
"""
import json
import os
import sys
import time

sys.path.insert(0, os.getcwd())
import numpy as np
import scipy.integrate as si

from src.kinetic_core import directional as d, ph_balance_refusal as s3, ph_balance_rivals as pbr, parallel, config
from src.kinetic_core import ph_state
from src.comparative_cli import spec_to_core

MODE = os.environ.get("PROBE_MODE", "baseline")
pbr.RIVAL_RHS_BUDGET = int(os.environ.get("PROBE_BUDGET", "200000"))
_orig = pbr.integrate_inline
_LOG = []
_TRACE = {}

if MODE in ("tightph", "coldtight"):
    ph_state.IONIC_STRENGTH_TOL_M = 1e-15
    _c0, _w0 = ph_state.solve_ph_and_ionic_strength_compiled, ph_state.solve_ph_and_ionic_strength_warm
    ph_state.solve_ph_and_ionic_strength_compiled = lambda T, sid, inv, g, i=0.0, **k: _c0(T, sid, inv, g, i, tol=1e-14)
    ph_state.solve_ph_and_ionic_strength_warm = lambda T, sid, inv, g, i=0.0, **k: _w0(T, sid, inv, g, i, tol=1e-14)
_FROZEN = {}
if MODE == "boundary45":
    from src.kinetic_core import parameters_sulfur as _psu
    _psu.ZHENG_MECHANISM_BOUNDARY_PH = 4.5
if MODE.startswith("snap"):
    from src.kinetic_core import sulfur as _S
    _kat = _S.sulfur_rate_constants_at
    _DEC = int(MODE[4:] or 9)

    def _kat_snap(parameters, temperature_k, ph, ph_label=None, **kw):
        return _kat(parameters, temperature_k, ph, ph_label=None if ph_label is None else round(float(ph_label), _DEC), **kw)

    _S.sulfur_rate_constants_at = _kat_snap
if MODE == "frozenph":
    _c1 = ph_state.solve_ph_and_ionic_strength_compiled

    def _frozen(T, sid, inv, g, i=0.0, **k):
        key = round(T, 6)
        if key not in _FROZEN:
            _FROZEN[key] = _c1(T, sid, inv, g, i, **k)
        return _FROZEN[key]

    ph_state.solve_ph_and_ionic_strength_compiled = _frozen
if MODE in ("coldph", "coldtight"):
    _c, _w = ph_state.solve_ph_and_ionic_strength_compiled, ph_state.solve_ph_and_ionic_strength_warm
    ph_state.solve_ph_and_ionic_strength_compiled = lambda T, sid, inv, g, i=0.0, **k: _c(T, sid, inv, 7.0, 0.0, **k)
    ph_state.solve_ph_and_ionic_strength_warm = lambda T, sid, inv, g, i=0.0, **k: _w(T, sid, inv, 7.0, 0.0, **k)

_solve_ivp = si.solve_ivp


def _traced_solve_ivp(fun, t_span, y0, **kw):
    tr = _TRACE
    tr.clear()
    tr.update(last_t=0.0, bins=np.zeros(20, int), jac_t=[], max_t=0.0)
    T = float(t_span[1])

    def f(t, y):
        tr["last_t"] = float(t)
        tr["max_t"] = max(tr["max_t"], float(t))
        tr["bins"][min(int(20 * t / T), 19)] += 1
        return fun(t, y)

    if "jac" in kw:
        j0 = kw["jac"]

        def j(t, y):
            if len(tr["jac_t"]) < 100000:
                tr["jac_t"].append(float(t))
            return j0(t, y)

        kw["jac"] = j
    if MODE in ("bdf", "radau"):
        kw["method"] = MODE.upper() if MODE == "bdf" else "Radau"
    if MODE.startswith("first"):
        kw["first_step"] = 10.0 ** (-int(MODE[5:]))
    if MODE == "atol12":
        kw["atol"] = 1e-12
    if MODE == "rtol8":
        kw["rtol"] = 1e-8
    return _solve_ivp(f, t_span, y0, **kw)


si.solve_ivp = _traced_solve_ivp

import scipy.integrate._ivp.lsoda as _L
_step0 = _L.LSODA._step_impl


def _step_traced(self):
    out = _step0(self)
    it = self._lsoda_solver._integrator
    st = _TRACE.setdefault("steps", [])
    st.append((float(self.t), float(it.rwork[10]), int(it.iwork[18]), int(it.iwork[13])))
    if len(st) > 4000:
        del st[1000:2000]
    return out


_L.LSODA._step_impl = _step_traced


def _trace_summary():
    jt = _TRACE.get("jac_t", [])
    st = _TRACE.get("steps", [])
    meth = [m for _, _, m, _ in st]
    return {"n_steps_kept": len(st), "adams_steps": meth.count(1), "bdf_steps": meth.count(2),
            "first_steps": st[:12], "last_steps": st[-12:],
            "last_t": _TRACE.get("last_t"), "max_t": _TRACE.get("max_t"),
            "bins": [int(v) for v in _TRACE.get("bins", [])], "n_jac": len(jt),
            "first_jac_t": jt[0] if jt else None}


def _wrapped(*a, **k):
    t = time.time()
    try:
        out = _orig(*a, **k)
    except pbr.RivalBudgetExceeded:
        _LOG.append({"budget_exceeded": True, "seconds": round(time.time() - t, 1), **_trace_summary()})
        raise
    diag = out[-1]
    _LOG.append({"evaluations": diag.get("evaluations"), "rhs": diag.get("rhs_calls"),
                 "method": diag.get("method"), "fallback": diag.get("solver_fallback"),
                 "seconds": round(time.time() - t, 2), **_trace_summary()})
    return out


pbr.integrate_inline = _wrapped


def job(item):
    name, spec, drift, cfg = item
    _LOG.clear()
    _FROZEN.clear()
    try:
        with config.applied(cfg):
            vals, ph = s3._rival_values((spec, ("MFT", "FFT"), drift))
        return name, dict(_LOG[-1]) if _LOG else {}, vals, ph, None
    except Exception as exc:  # noqa: BLE001
        return name, dict(_LOG[-1]) if _LOG else {}, None, None, f"{type(exc).__name__}: {exc}"[:200]


p = d.load_panel(d.PANEL_PATH)["panel"]
arms = d._arms(p[0])
perturb = float(os.environ.get("PROBE_PERTURB", "0"))
reps = int(os.environ.get("PROBE_REPS", "1"))
arm_filter = os.environ.get("PROBE_ARM", "")
base = config.Config.from_globals()
items = []
for rep in range(reps):
    for label, system in arms:
        if arm_filter and arm_filter not in label:
            continue
        sp = d._spec_from_system(label, system)
        eps = perturb * rep if reps > 1 else perturb
        if eps:
            sp["precursors"] = {k: v * (1 + eps) for k, v in sp["precursors"].items()}
        spec = spec_to_core(sp)
        for point, (changes, drift) in s3.rival_points().items():
            items.append((f"{label}|{point}|rep{rep}", spec, drift, base.replace(changes)))
only = os.environ.get("PROBE_ONLY")
if only:
    items = [it for it in items if any(o in it[0] for o in only.split(","))]
t0 = time.time()
res = parallel.map_ordered(job, items, n_workers=int(os.environ.get("PROBE_WORKERS", "1")))
print("MODE", MODE, "budget", pbr.RIVAL_RHS_BUDGET, "perturb", perturb)
for name, diag, vals, ph, err in res:
    print("ROW " + json.dumps({"name": name, "diag": diag, "vals": vals, "ph": ph, "err": err}), flush=True)
import collections
print("ERRKINDS", collections.Counter((r[4] or "")[:60] for r in res if r[4]))
print("EVALS median/max", float(np.median([r[1].get("evaluations") or 0 for r in res if not r[4]])), max((r[1].get("evaluations") or 0) for r in res))
print("FALLBACKS", sum(1 for r in res if r[1].get("fallback")), [r[0] for r in res if r[1].get("fallback")])
print("STALLS", MODE, sum(1 for r in res if r[4]), "of", len(res))
print("wall", round(time.time() - t0, 1), "max evals", max((r[1].get("evaluations") or 0) for r in res),
      "failures", [r[0] for r in res if r[4]])
