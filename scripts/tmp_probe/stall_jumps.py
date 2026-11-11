"""Record (t, y, f) over consecutive RHS calls in a stalled rival integration; which components of f jump."""
import os, sys
sys.path.insert(0, os.getcwd())
import numpy as np
import scipy.integrate as si
from src.kinetic_core import directional as d, ph_balance_refusal as s3, ph_balance_rivals as pbr, config
from src.kinetic_core import sulfur as S
from src.comparative_cli import spec_to_core

pbr.RIVAL_RHS_BUDGET = 10 ** 9
names = list(S.SULFUR_INDEX) + list(pbr.TRACKERS)
START = int(os.environ.get("STALL_AT", "15000"))
N = 40
rec = []


class Done(Exception):
    pass


_solve = si.solve_ivp


def solve_ivp(fun, t_span, y0, **kw):
    c = {"n": 0}

    def f(t, y):
        c["n"] += 1
        out = fun(t, y)
        if c["n"] >= START:
            rec.append((float(t), np.array(y, float), np.array(out, float)))
            if len(rec) >= N:
                yy = np.array(y, float)
                w = 1.0 / (kw["rtol"] * np.abs(yy) + kw["atol"])
                J = kw["jac"](t, yy)
                Jn = np.zeros_like(J)
                for i in range(yy.size):
                    h = max(1e-7 * abs(yy[i]), 1e-20)
                    yp = yy.copy(); yp[i] += h
                    Jn[:, i] = (fun(t, yp) - out) / h
                ic, ifr = names.index("Cys"), names.index("FRAG_C")
                for rel in (1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13):
                    for sgn in (1, -1):
                        yp = yy.copy(); yp[ic] *= (1 + sgn * rel)
                        df = fun(t, yp) - out
                        print(f"  Cys rel {sgn*rel:+.0e}: dF[FRAG_C]={df[ifr]:+.4e}  dF/dy={df[ifr]/(yp[ic]-yy[ic]):+.4e}  max|w dF|={np.max(np.abs(w*df)):.3e} at {names[int(np.argmax(np.abs(w*df)))]}")
                fun(t, yy); A = dict(LAST)
                yp = yy.copy(); yp[ic] *= (1 + 1e-6); fun(t, yp); B = dict(LAST)
                print("pH", A["ph"], B["ph"], A["kw"], B["kw"])
                fA = fun(t, yy); fB = fun(t, yp)
                for tau in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
                    fun(t, yy + tau * fA); la = LAST["kw"]["ph_label"]
                    fun(t, yp + tau * fB); lb = LAST["kw"]["ph_label"]
                    print(f"  tau={tau:.0e}: label pH along below-branch flow {la - A['kw']['ph_label']:+.3e}, along above-branch flow {lb - B['kw']['ph_label']:+.3e}")
                dk = sorted(((abs(B["k"][q] - A["k"][q]) / max(abs(A["k"][q]), 1e-300), q, A["k"][q], B["k"][q]) for q in A["k"]), reverse=True)[:6]
                print("largest relative k changes:", dk)
                dr = np.abs(B["r"] - A["r"]); top = np.argsort(-dr)[:6]
                print("largest rate changes:", [(S.FULL_REACTION_KEYS[q], float(A["r"][q]), float(B["r"][q])) for q in top])
                for nm, M in (("analytic", J), ("numerical", Jn)):
                    W = (w[:, None] * M) / w[None, :]
                    print(nm, "spectral radius", float(np.max(np.abs(np.linalg.eigvals(M)))),
                          "weighted 2-norm", float(np.linalg.norm(W, 2)), "max weighted entry", float(np.abs(W).max()))
                    i, j = np.unravel_index(np.argmax(np.abs(W)), W.shape)
                    print("   largest weighted entry: d f[%s] / d y[%s] = %.3e (w_i %.1e, w_j %.1e)" % (
                        names[i], names[j], M[i, j], w[i], w[j]))
                raise Done((kw["rtol"], kw["atol"]))
        return out

    return _solve(f, t_span, y0, **kw)


si.solve_ivp = solve_ivp
LAST = {}
_kat, _rv = S.sulfur_rate_constants_at, S._rates_vectorised


def _kat2(parameters, temperature_k, ph, **kw):
    out = _kat(parameters, temperature_k, ph, **kw)
    LAST.update(ph=ph, kw=dict(kw), k=dict(out))
    return out


def _rv2(y, kv):
    r = _rv(y, kv)
    LAST["r"] = r.copy()
    return r


S.sulfur_rate_constants_at, S._rates_vectorised = _kat2, _rv2
p = d.load_panel(d.PANEL_PATH)["panel"]
arm_key, point = os.environ["PROBE_ONLY"].split("|")
eps = float(os.environ.get("PROBE_PERTURB", "0"))
label, system = next(a for a in d._arms(p[0]) if arm_key in a[0])
sp = d._spec_from_system(label, system)
sp["precursors"] = {k: v * (1 + eps) for k, v in sp["precursors"].items()}
spec = spec_to_core(sp)
changes, drift = s3.rival_points()[point]
try:
    with config.applied(config.Config.from_globals().replace(changes)):
        s3._rival_values((spec, ("MFT", "FFT"), drift))
    print("finished")
    sys.exit()
except Done as e:
    rtol, atol = e.args[0]
names = list(S.SULFUR_INDEX) + list(pbr.TRACKERS)
ts = np.array([r[0] for r in rec]); Y = np.array([r[1] for r in rec]); F = np.array([r[2] for r in rec])
w = 1.0 / (rtol * np.abs(Y[0]) + atol)
print("t range", ts.min(), ts.max(), "dt between calls", np.diff(ts)[:10])
dF = np.abs(np.diff(F, axis=0)) * w
score = dF.max(axis=0)
for i in np.argsort(-score)[:8]:
    print(f"{names[i]:>12s} y={Y[-1, i]: .3e} w*max|dF|={score[i]:.3e} f:", np.array2string(F[:12, i], precision=4))
neg = np.flatnonzero((Y < 0).any(axis=0))
print("species negative at some call:", [(names[i], float(Y[:, i].min())) for i in neg])
