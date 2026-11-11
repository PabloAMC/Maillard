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

if MODE == "coldph":
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
    if MODE == "rtol8":
        kw["rtol"] = 1e-8
    return _solve_ivp(f, t_span, y0, **kw)


si.solve_ivp = _traced_solve_ivp


def _trace_summary():
    jt = _TRACE.get("jac_t", [])
    return {"last_t": _TRACE.get("last_t"), "max_t": _TRACE.get("max_t"),
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
                 "seconds": round(time.time() - t, 2), **_trace_summary()})
    return out


pbr.integrate_inline = _wrapped


def job(item):
    name, spec, drift, cfg = item
    _LOG.clear()
    try:
        with config.applied(cfg):
            vals, ph = s3._rival_values((spec, ("MFT", "FFT"), drift))
        return name, dict(_LOG[-1]) if _LOG else {}, vals, ph, None
    except Exception as exc:  # noqa: BLE001
        return name, dict(_LOG[-1]) if _LOG else {}, None, None, f"{type(exc).__name__}: {exc}"[:200]


p = d.load_panel(d.PANEL_PATH)["panel"]
arms = d._arms(p[0])
perturb = float(os.environ.get("PROBE_PERTURB", "0"))
base = config.Config.from_globals()
items = []
for label, system in arms:
    sp = d._spec_from_system(label, system)
    if perturb:
        sp["precursors"] = {k: v * (1 + perturb) for k, v in sp["precursors"].items()}
    spec = spec_to_core(sp)
    for point, (changes, drift) in s3.rival_points().items():
        items.append((f"{label}|{point}", spec, drift, base.replace(changes)))
only = os.environ.get("PROBE_ONLY")
if only:
    items = [it for it in items if any(o in it[0] for o in only.split(","))]
t0 = time.time()
res = parallel.map_ordered(job, items, n_workers=int(os.environ.get("PROBE_WORKERS", "1")))
print("MODE", MODE, "budget", pbr.RIVAL_RHS_BUDGET, "perturb", perturb)
for name, diag, vals, ph, err in res:
    print("ROW " + json.dumps({"name": name, "diag": diag, "vals": vals, "ph": ph, "err": err}), flush=True)
print("wall", round(time.time() - t0, 1), "max evals", max((r[1].get("evaluations") or 0) for r in res),
      "failures", [r[0] for r in res if r[4]])
