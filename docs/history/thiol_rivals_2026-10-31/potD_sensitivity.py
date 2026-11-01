"""At FIXED route supplies (same nuisance draw), how far does each loss law move Cerny pot D's NF share relative to no loss?"""
import sys, numpy as np
sys.path.insert(0, sys.argv[1])
import mft_loss as m
rng = np.random.default_rng(7); n = 2000
t = np.arange(0, 240.5, 0.5)
sr = m.supply(rng, n, t, (5, 20), (0.0, 1.5), (120, 1800)); tn = m.logu(rng, 60, 2000, n)[:, None]
sn = t[None, :] * np.exp(-t[None, :] / tn); rho = m.logu(rng, 0.003, 1.0, n)[:, None]
sr = sr / np.trapezoid(sr, t, axis=1)[:, None]; sn = rho * sn / np.trapezoid(sn, t, axis=1)[:, None]
cys = m.cys_series(rng, n, t, (10.0, 97.0), (100.0, 2000.0))
base = None
for rv in m.RIVALS:
    d = m.draw(rv, rng, n); lm = m.lam(d, 95.0, t, cys)
    a, b = m.integrate(sn, lm, 0.5)[:, -1], m.integrate(sr, lm, 0.5)[:, -1]; sh = a / (a + b)
    if rv == "NCa": base = sh
    print(rv, "median share", round(float(np.median(sh)), 4), "median lambda(95C)", float(np.median(lm[:, -1])))
    if base is not None: pass
    globals()["sh_" + rv] = sh
for rv in m.RIVALS[:-1]:
    r = np.log10(globals()["sh_" + rv] / base)
    print(rv, "log10 share/no-loss: p05 %.3f median %.3f p95 %.3f" % tuple(np.percentile(r, [5, 50, 95])))
print("nuisance rho spread of share (no loss): p05 %.4f p95 %.4f" % tuple(np.percentile(base, [5, 95])))
