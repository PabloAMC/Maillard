"""At the stalled state of a rival integration: eigenvalues of the analytic Jacobian, and the smoothness of f.

Env: PROBE_MODE (coldph | baseline), PROBE_ONLY ('arm|point'), STALL_AT (rhs call count to analyse at).
"""
import json
import os
import sys

sys.path.insert(0, os.getcwd())
import numpy as np
import scipy.integrate as si

from src.kinetic_core import directional as d, ph_balance_refusal as s3, ph_balance_rivals as pbr, config, ph_state
from src.kinetic_core import sulfur as S
from src.comparative_cli import spec_to_core

MODE = os.environ.get("PROBE_MODE", "coldph")
STALL_AT = int(os.environ.get("STALL_AT", "20000"))
pbr.RIVAL_RHS_BUDGET = 10 ** 9
if MODE == "coldph":
    _c = ph_state.solve_ph_and_ionic_strength_compiled
    ph_state.solve_ph_and_ionic_strength_compiled = lambda T, sid, inv, g, i=0.0, **k: _c(T, sid, inv, 7.0, 0.0, **k)


class Done(Exception):
    pass


_solve = si.solve_ivp


def analyse(fun, jac, t, y, rtol, atol):
    names = list(S.SULFUR_INDEX) + list(pbr.TRACKERS)
    w = 1.0 / (rtol * np.abs(y) + atol)
    f0 = fun(t, y)
    J = jac(t, y)
    ev = np.linalg.eigvals(J)
    order = np.argsort(ev.real)[:5]
    print("t", t, "n", y.size)
    print("most negative eigenvalues (analytic J, pH frozen):", [complex(round(e.real, 3), round(e.imag, 3)) for e in ev[order]])
    print("max Re eig", float(ev.real.max()))
    # numerical J by central differences, per column with step relative to the weight scale
    Jn = np.zeros_like(J)
    for i in range(y.size):
        h = 1e-6 / w[i]
        yp, ym = y.copy(), y.copy()
        yp[i] += h
        ym[i] -= h
        Jn[:, i] = (fun(t, yp) - fun(t, ym)) / (2 * h)
    evn = np.linalg.eigvals(Jn)
    on = np.argsort(evn.real)[:5]
    print("most negative eigenvalues (numerical J, pH live):", [complex(round(e.real, 3), round(e.imag, 3)) for e in evn[on]])
    print("max Re eig numerical", float(evn.real.max()))
    # smoothness along a random direction scaled to the error weights
    rng = np.random.default_rng(0)
    v = rng.standard_normal(y.size) / w
    print("directional derivative ||(f(y+s v)-f(y))/s||_w vs s (v ~ 1 error-unit):")
    for s in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]:
        g = (fun(t, y + s * v) - f0) / s
        print(f"  s={s:.0e}  ||g||_w={np.sqrt(np.mean((g * w) ** 2)):.4e}")
    # noise: repeat f at identical (t, y) after evaluating elsewhere (history dependence)
    a = fun(t, y)
    fun(t, y + 1e-3 * v)
    b = fun(t, y)
    print("history: max |f(y) - f(y) after a detour| weighted:", float(np.max(np.abs(a - b) * w)))
    # species near zero / negative
    neg = [(names[i], float(y[i])) for i in np.flatnonzero(y < 0)]
    print("negative species:", neg[:20])
    top = np.argsort(-np.abs(f0 * w))[:8]
    print("largest weighted |f|:", [(names[i], float(f0[i]), float(y[i])) for i in top])
    # the dominant fast eigenvector's species
    i0 = order[0]
    vals, vecs = np.linalg.eig(J)
    k = np.argmin(vals.real)
    vv = np.abs(vecs[:, k])
    print("fastest mode", vals[k], "species:", [(names[i], round(float(vv[i]), 3)) for i in np.argsort(-vv)[:6]])


def solve_ivp(fun, t_span, y0, **kw):
    calls = {"n": 0}
    jac = kw.get("jac")

    def f(t, y):
        calls["n"] += 1
        if calls["n"] == STALL_AT:
            analyse(fun, jac, t, np.array(y, float), kw["rtol"], kw["atol"])
            raise Done()
        return fun(t, y)

    return _solve(f, t_span, y0, **kw)


si.solve_ivp = solve_ivp

p = d.load_panel(d.PANEL_PATH)["panel"]
arm_key, point = os.environ.get("PROBE_ONLY", "glucose|B1_d446").split("|")
label, system = next(a for a in d._arms(p[0]) if arm_key in a[0])
spec = spec_to_core(d._spec_from_system(label, system))
changes, drift = s3.rival_points()[point]
try:
    with config.applied(config.Config.from_globals().replace(changes)):
        s3._rival_values((spec, ("MFT", "FFT"), drift))
    print("finished before STALL_AT")
except Done:
    pass
