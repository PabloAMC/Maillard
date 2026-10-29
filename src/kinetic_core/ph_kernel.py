"""
src/kinetic_core/ph_kernel.py -- THE CHARGE BALANCE, COMPILED (C5, 2026-10-28).

The inline-coupled sulfur integration (``sulfur.PH_COUPLING = "inline"``) solves the pot's pH and ionic
strength at every right-hand-side evaluation, and in pure Python that solve was ~80 % of a run. This module
is the same arithmetic as ``ph_state._residual_at_fixed_conditions`` and ``ph_state.ionic_strength_of`` --
the polyprotic distribution in log form, Davies activity, the conserved strong-ion difference -- on arrays,
compiled with numba when it is installed and run as plain Python when it is not, so the answer never depends
on the package being there. The root is a bisection inside a bracket grown from the previous call's pH,
to ``tol`` pH; the ionic strength is iterated to ``IONIC_STRENGTH_TOL_M`` as in ``ph_state``.

Nothing here is a parameter. The pK values at the pot's temperature (van 't Hoff with dCp) and the Davies A
are computed by ``ph_state`` and passed in, so the physics has one definition.
"""
from __future__ import annotations

import math

import numpy as np

try:  # pragma: no cover - which branch runs depends on the environment
    from numba import njit

    COMPILED = True
except ImportError:  # pragma: no cover
    COMPILED = False

    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda f: f


@njit(cache=True)
def _davies(ionic_m, a, on):
    if not on or ionic_m <= 0.0:
        return 0.0
    s = math.sqrt(ionic_m)
    return a * (s / (1.0 + s) - 0.3 * ionic_m)


@njit(cache=True)
def _weights(ph, pk, n, zmax, af, out):
    """The n + 1 species fractions of one group at ``ph``, into ``out``; pK' = pK + A f(I) (2 (z_max - i) - 1)."""
    out[0] = 0.0
    running = 0.0
    top = 0.0
    for i in range(n):
        running += ph - (pk[i] + af * (2.0 * (zmax - i) - 1.0))
        out[i + 1] = running
        if running > top:
            top = running
    total = 0.0
    for i in range(n + 1):
        out[i] = 10.0 ** (out[i] - top)
        total += out[i]
    for i in range(n + 1):
        out[i] /= total


@njit(cache=True)
def _residual(ph, pk, npk, zmax, conc, sid, pkw, a, on, ionic, work):
    af = _davies(ionic, a, on)
    inv_gamma = 10.0 ** af
    total = sid + 1000.0 * 10.0 ** (-ph) * inv_gamma - 1000.0 * 10.0 ** (ph - pkw) * inv_gamma
    for g in range(conc.shape[0]):
        n = npk[g]
        _weights(ph, pk[g], n, zmax[g], af, work)
        nbar = 0.0
        for i in range(n + 1):
            nbar += i * work[i]
        total += conc[g] * (zmax[g] - nbar)
    return total


@njit(cache=True)
def _ionic(ph, pk, npk, zmax, conc, sid, pkw, a, on, ionic, work):
    af = _davies(ionic, a, on)
    inv_gamma = 10.0 ** af
    two_i = abs(sid) + 1000.0 * 10.0 ** (-ph) * inv_gamma + 1000.0 * 10.0 ** (ph - pkw) * inv_gamma
    for g in range(conc.shape[0]):
        n = npk[g]
        _weights(ph, pk[g], n, zmax[g], af, work)
        contribution = 0.0
        for i in range(n + 1):
            z = zmax[g] - i
            contribution += work[i] * z * z
        two_i += conc[g] * contribution
    return 0.5 * two_i / 1000.0


@njit(cache=True)
def _root(pk, npk, zmax, conc, sid, pkw, a, on, ionic, guess, tol, work):
    step = 0.05
    lo = max(guess - step, -1.0)
    hi = min(guess + step, 15.0)
    f_lo = _residual(lo, pk, npk, zmax, conc, sid, pkw, a, on, ionic, work)
    f_hi = _residual(hi, pk, npk, zmax, conc, sid, pkw, a, on, ionic, work)
    while f_lo < 0.0 and lo > -1.0:
        step *= 4.0
        lo = max(lo - step, -1.0)
        f_lo = _residual(lo, pk, npk, zmax, conc, sid, pkw, a, on, ionic, work)
    while f_hi > 0.0 and hi < 15.0:
        step *= 4.0
        hi = min(hi + step, 15.0)
        f_hi = _residual(hi, pk, npk, zmax, conc, sid, pkw, a, on, ionic, work)
    if f_lo < 0.0:
        return lo
    if f_hi > 0.0:
        return hi
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if _residual(mid, pk, npk, zmax, conc, sid, pkw, a, on, ionic, work) > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


@njit(cache=True)
def solve(pk, npk, zmax, conc, sid, pkw, a, on, ph_guess, ionic_guess, tol, ionic_tol, sweeps):
    """(pH, I in mol/L) at the conserved SID: the root of the charge residual, I iterated to its fixed point."""
    work = np.empty(pk.shape[1] + 1)
    if not on:
        ph = _root(pk, npk, zmax, conc, sid, pkw, a, on, 0.0, ph_guess, tol, work)
        return ph, _ionic(ph, pk, npk, zmax, conc, sid, pkw, a, on, 0.0, work)
    ionic = ionic_guess
    ph = _root(pk, npk, zmax, conc, sid, pkw, a, on, ionic, ph_guess, tol, work)
    for _ in range(sweeps):
        fresh = _ionic(ph, pk, npk, zmax, conc, sid, pkw, a, on, ionic, work)
        ph = _root(pk, npk, zmax, conc, sid, pkw, a, on, fresh, ph, tol, work)
        moved = abs(fresh - ionic)
        ionic = fresh
        if moved < ionic_tol:
            break
    return ph, _ionic(ph, pk, npk, zmax, conc, sid, pkw, a, on, ionic, work)
