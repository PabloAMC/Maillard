"""
src/kinetic_core_jax/densities.py -- the bridge from the JAX network to M1.1's exact predictive densities
(``src.kinetic_core.predictive_densities``; A2 contract v2.2; 2026-11-05).

``predictive_densities`` is pure: every estimator (``importance_lppd``, ``quadrature_lppd``, ``plug_in_lppd``) takes
the log integrand log L_u(theta_s) + log p(theta_s) of each parameter draw, and the caller supplies the predictions
behind L. M1b / M1c / M3 get them from the engine, one LSODA run per draw and pot. Where the JAX network covers a
unit's pots, this module supplies them instead: ``Predictor`` maps a (D, p) array of draws of a coordinate vector to
the (D, n) predictions of a unit's rows, every draw of a pot solved in one compiled batch (``jax.vmap`` over draws,
chunked so the memory stays bounded), and ``log_likelihoods`` turns them into log L per draw under the unit's
covariance -- the same number ``predictive_densities.gaussian_log_density`` gives row by row (one Cholesky for all
draws). ``log_integrand`` composes both with a log prior for ``quadrature_lppd``'s callback.

FAILED DRAWS (an integration that did not finish) are NOT silently scored: ``Predictor`` returns NaN predictions with a
mask, and ``log_likelihoods`` refuses them unless the caller declares the policy (``failed="impossible"``: -inf, which
``importance_lppd`` counts in S and leaves out of the tail fit, as its docstring requires).
"""
from __future__ import annotations

import math
import time
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import jax
import jax.numpy as jnp
import numpy as np

from .fit import LOG_FLOOR, Observation
from .network import Network, integrate, n_rows, pot_arrays, pot_mode, theta_from_vector, theta_vector

FAILED_POLICIES = ("refuse", "impossible")


class Predictor:
    """Per-draw predictions of ``observations`` (``fit.Observation``: pot, row, observable, log10; value and sigma are
    not read) on ``pots`` for draws of the coordinates ``keys`` (the others at ``theta`` / the shipped point).
    ``predictor(draws)`` -> (predictions (D, n), failed (D,)). One compile per pot shape and chunk size."""

    def __init__(self, net: Network, pots: Sequence[Mapping[str, Any]], observations: Sequence[Observation],
                 keys: Sequence[Tuple[str, str]], *, theta=None, forcing=None, chunk: int = 512,
                 rtol: float = 1e-8, atol: float = 1e-14, solver=None, max_steps: Optional[int] = None):
        self.net, self.pots, self.observations = net, list(pots), list(observations)
        self.keys = [tuple(k) for k in keys]
        where = {c: i for i, c in enumerate(net.coordinates)}
        missing = [k for k in self.keys if k not in where]
        if missing:
            raise ValueError(f"not coordinates of the network: {missing}")
        base = jnp.asarray(theta_vector(net, net.theta0() if theta is None else theta))
        free = jnp.asarray([where[k] for k in self.keys], dtype=jnp.int64)
        self.chunk = int(chunk)
        table = {k: f for m in net.modules for k, f in m.observables.items()}
        unknown = sorted({o.observable for o in self.observations} - set(table))
        if unknown:
            raise ValueError(f"unknown observables {unknown}")
        if max_steps is None and forcing is not None:
            max_steps = 100_000                     # integrate picks its own default (pH-free 1e5, inline pH 2e5)
        solve_kw = dict(rtol=rtol, atol=atol, solver=solver, max_steps=max_steps, throw=False)

        def batch_for(mode):
            def one(x, pot):
                th = theta_from_vector(net, base.at[free].set(x))
                if forcing is None:
                    run = integrate(net, th, pot, mode=mode, **solve_kw)
                else:
                    from .forcing import integrate_forced

                    run = integrate_forced(net, th, pot, forcing, **solve_kw)
                return run.ys, jnp.any(run.failed)

            return jax.jit(jax.vmap(one, in_axes=(0, None)))

        # (v2.3) an inline-pH pot's static mode is a compile key; its numeric arrays are the argument. A row is a
        # segment end, or a pH node of an inline-pH pot; one outside the pot's output would be CLAMPED by JAX's
        # indexing (a silently wrong prediction), so it is refused here
        n_out = [n_rows(net, p) for p in self.pots]
        stray = [(i, o.pot, o.row, n_out[o.pot]) for i, o in enumerate(self.observations)
                 if not 0 <= int(o.row) < n_out[o.pot]]
        if stray:
            raise ValueError(f"observation rows outside their pot's output (observation, pot, row, n_rows): {stray[:5]}")
        self._modes = [pot_mode(net, p) for p in self.pots]
        self._arrays = [pot_arrays(net, p) for p in self.pots]
        self._batch = {mode: batch_for(mode) for mode in dict.fromkeys(self._modes)}
        # per pot: the rows it holds, evaluated on the (D, S, n) batch at once
        self._rows: Dict[int, List[Tuple[int, int, Callable, bool]]] = {}
        for i, o in enumerate(self.observations):
            self._rows.setdefault(o.pot, []).append((i, o.row, table[o.observable], o.log10))
        self._obs = {}
        for p, rows in self._rows.items():
            fns = [(row, fn, log10) for _, row, fn, log10 in rows]

            def evaluate(ys, fns=fns):
                cols = []
                for row, fn, log10 in fns:
                    v = jax.vmap(lambda y: fn(y, net.index))(ys[:, row])
                    cols.append(jnp.log10(jnp.maximum(v, LOG_FLOOR)) if log10 else v)
                return jnp.stack(cols, axis=1)

            self._obs[p] = jax.jit(evaluate)
        self.seconds: Dict[str, float] = {}

    def __call__(self, draws) -> Tuple[np.ndarray, np.ndarray]:
        draws = np.atleast_2d(np.asarray(draws, dtype=float))
        if draws.shape[1] != len(self.keys):
            raise ValueError(f"draws have {draws.shape[1]} columns for {len(self.keys)} coordinates")
        d = draws.shape[0]
        out = np.full((d, len(self.observations)), np.nan)
        failed = np.zeros(d, dtype=bool)
        t0 = time.perf_counter()
        for p, rows in self._rows.items():
            pot, batch = self._arrays[p], self._batch[self._modes[p]]
            cols = [i for i, *_ in rows]
            for a in range(0, d, self.chunk):
                block = draws[a: a + self.chunk]
                n = block.shape[0]
                if n < self.chunk and d > self.chunk:      # pad the last chunk: one compile per chunk size
                    block = np.concatenate([block, np.repeat(block[-1:], self.chunk - n, axis=0)])
                ys, bad = batch(jnp.asarray(block), pot)
                values = np.asarray(self._obs[p](ys))[:n]
                bad = np.asarray(bad)[:n]
                out[a: a + n, cols] = values
                failed[a: a + n] |= bad
        out[failed] = np.nan
        self.seconds["last_call"] = time.perf_counter() - t0
        return out, failed


def log_likelihoods(predictions: np.ndarray, y: Sequence[float], cov: np.ndarray, *, failed: Optional[np.ndarray] = None,
                    on_failed: str = "refuse") -> np.ndarray:
    """log N(y; f_s, cov) for every draw's prediction row f_s, nats, in the basis y / f / cov share; equal to
    ``predictive_densities.gaussian_log_density(y, f_s, cov)`` for each s. A failed draw (NaN row, or ``failed``) is
    refused unless ``on_failed="impossible"`` (-inf)."""
    from src.kinetic_core.predictive_densities import MeasurementModelRefusal

    if on_failed not in FAILED_POLICIES:
        raise ValueError(f"on_failed must be one of {FAILED_POLICIES}")
    f = np.atleast_2d(np.asarray(predictions, dtype=float))
    y = np.asarray(y, dtype=float)
    cov = np.atleast_2d(np.asarray(cov, dtype=float))
    bad = ~np.all(np.isfinite(f), axis=1)
    if failed is not None:
        bad |= np.asarray(failed, dtype=bool)
    if bad.any() and on_failed == "refuse":
        raise MeasurementModelRefusal(f"{int(bad.sum())} draws failed or are non-finite: declare the policy "
                                      "(on_failed='impossible') before scoring")
    try:
        chol = np.linalg.cholesky(cov)
    except np.linalg.LinAlgError as exc:
        raise MeasurementModelRefusal("the covariance is not positive definite") from exc
    logdet = 2.0 * float(np.sum(np.log(np.diag(chol))))
    r = np.where(bad[:, None], 0.0, y[None, :] - f)
    from scipy.linalg import solve_triangular

    z = solve_triangular(chol, r.T, lower=True)
    out = -0.5 * (np.sum(z * z, axis=0) + logdet + len(y) * math.log(2.0 * math.pi))
    return np.where(bad, -np.inf, out)


def log_integrand(predictor: Predictor, y: Sequence[float], cov: np.ndarray,
                  log_prior: Callable[[np.ndarray], np.ndarray], *, on_failed: str = "refuse"
                  ) -> Callable[[np.ndarray], np.ndarray]:
    """``g(theta (N, p)) -> log L + log p`` per row: the callback ``predictive_densities.quadrature_lppd`` takes, and
    the ``log_integrand`` array ``importance_lppd`` takes (``g(draws)``)."""

    def g(theta):
        theta = np.atleast_2d(np.asarray(theta, dtype=float))
        pred, failed = predictor(theta)
        return log_likelihoods(pred, y, cov, failed=failed, on_failed=on_failed) + np.asarray(log_prior(theta))

    return g
