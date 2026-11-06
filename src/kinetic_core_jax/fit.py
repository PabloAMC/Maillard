"""
src/kinetic_core_jax/fit.py -- the JAX fit driver (A2 contract v2.2 sec. 2.7 and 6.4; 2026-11-05): the counterpart of
the engine's ``src.kinetic_core.fitting.fit`` on a composed ``Network``.

WHAT IT DOES. A ``JaxProblem`` is a coordinate vector (a subset of ``Network.coordinates``, with bounds and kinds like
the engine's ``Problem``) over many pots, and weighted residuals r_i = (f_i - y_i) / sigma_i of named observables
(log10 of the level by default) at segment ends. ``fit`` runs scipy's ``least_squares`` with the engine's own settings
(``fitting.LEAST_SQUARES_SETTINGS``: trf, x_scale = "jac", ftol 1e-6, xtol 1e-8) on the JAX residual and its EXACT
forward-mode Jacobian (``jax.jacfwd`` through ``diffrax.ForwardMode()``: one JVP per coordinate, contract sec. 2.7),
inside an evaluation budget, and asserts convergence (``require_converged``, the engine's semantics: status <= 0 is
not convergence; also refused: a pot that failed to integrate at the optimum).

WHY NO CENTRAL POLISH. M1.1's ``central_polish`` exists because the engine's forward-difference Jacobian (B8's
absolute steps) is up to 7 % off its own central difference on the stiff pots (contract sec. 2.7), which stops the
optimiser short of the optimum. The Jacobian here is exact to the solver's tolerance, so the first trf run ends where
the polished engine fit ends (``build_jax_network_infra.py``'s equivalence stage measures it).

BATCHES. Pots are grouped by the shapes of their arrays (segment count, forcing); each group runs under ``jax.vmap``
(``batch="vmap"``: one compile per group, members in lockstep) or one pot at a time (``batch="loop"``: one compile per
shape, no slowest-member wait). On CPU vmap removes dispatch overhead only (contract sec. 2.8); both are offered and
timed. Pot data are ARGUMENTS, never traced constants: one compile serves every data set and start of a structure.

FAILURES. Pots integrate with ``throw=False``; a failed pot's residuals are ``failed_residual`` (25, the engine's
``Problem.failed_residual``) with a zero Jacobian row, and ``JaxFit.failed_pots`` names it.
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import jax
import jax.numpy as jnp
import numpy as np

from .network import Network, integrate, n_rows, pot_arrays, pot_mode, theta_from_vector, theta_vector

#: coordinate kinds: the engine's (``fitting.THRESHOLDS`` / ``EXTRA_THRESHOLDS``) plus "linear" (a k_ref the engine
#: carries at zero, ``ArrheniusBlock(allow_zero=True)``; bounded below by 0) and "share" (a fraction in [0, 1])
ENGINE_KINDS = ("Ea", "log10k", "acid_yield", "pka", "ph_slope", "log10_offset")
KINDS = ENGINE_KINDS + ("linear", "share", "q10")
FAILED_RESIDUAL = 25.0
#: a log10 observable of a non-positive prediction reads this floor (the residual is then large, not NaN)
LOG_FLOOR = 1e-300
#: the contract's fit-tolerance rule (sec. 2.6): fits run at rtol 1e-6 once the panel's residuals move < 1e-3 dex
#: between rtol 1e-6 and 1e-8
FIT_RTOL_CHECK_DEX = 1e-3


def kind_of(coordinate: str) -> Optional[str]:
    """The kind a module coordinate's name implies (None: declare it)."""
    head = coordinate.partition(":")[0]
    if head == "log10_k" or coordinate.startswith("log10_k"):
        return "log10k"
    if head == "ea" or coordinate.startswith("ea_"):
        return "Ea"
    if head == "k_ref":
        return "linear"
    if coordinate.startswith(("ph_exponent", "pyrazine_ph_slope")):
        return "ph_slope"
    if coordinate.startswith("share") or coordinate == "donor":
        return "share"
    if coordinate == "q10":
        return "q10"
    return None


@dataclass(frozen=True)
class Observation:
    """One residual: observable ``observable`` (a Network observable name, e.g. "mmol_per_L:HMF") of pot ``pot`` at
    output row ``row``, against ``value`` with ``sigma``, both in the residual's basis (log10 of the observable when
    ``log10``). A row is a segment end for a plain or forced pot, and a pH NODE for an inline-pH pot (S * N rows,
    segment-major; ``network.segment_end_rows`` gives the segment ends); ``network.n_rows`` bounds it, checked when the
    problem compiles."""

    pot: int
    row: int
    observable: str
    value: float
    sigma: float
    log10: bool = True


@dataclass(frozen=True)
class JaxProblem:
    net: Network
    keys: Tuple[Tuple[str, str], ...]               # (module, coordinate), the free vector's order
    kinds: Tuple[str, ...]
    lower: np.ndarray
    upper: np.ndarray
    pots: Tuple[Mapping[str, Any], ...]
    observations: Tuple[Observation, ...]
    #: the fixed coordinates (a full theta pytree); default the network's shipped point
    theta: Optional[Mapping[str, Mapping[str, Any]]] = None
    #: a ``forcing.Forcing`` when the pots are forced pots (``forcing.forced_pot(...).pot``)
    forcing: Any = None
    solve: Mapping[str, Any] = field(default_factory=lambda: {"rtol": 1e-8, "atol": 1e-14})
    batch: str = "vmap"
    failed_residual: float = FAILED_RESIDUAL

    def __post_init__(self):
        assert len(self.keys) == len(self.kinds) == len(self.lower) == len(self.upper), "keys, kinds, bounds align"
        bad = [k for k in self.kinds if k not in KINDS]
        assert not bad, f"unknown coordinate kinds {bad} (declared: {KINDS})"
        assert np.all(np.asarray(self.lower) < np.asarray(self.upper)), "every band must have positive width"
        coords = set(self.net.coordinates)
        missing = [k for k in self.keys if tuple(k) not in coords]
        assert not missing, f"not coordinates of the network: {missing}"
        for k, kind, lo in zip(self.keys, self.kinds, self.lower):
            if kind == "linear":
                assert lo >= 0.0, f"{k}: a linear rate constant is bounded below by 0"
        assert self.batch in ("vmap", "loop"), "batch is 'vmap' or 'loop'"
        for o in self.observations:
            assert 0 <= o.pot < len(self.pots) and o.sigma > 0, o

    def base_vector(self) -> np.ndarray:
        theta = self.net.theta0() if self.theta is None else self.theta
        return np.asarray(theta_vector(self.net, theta), dtype=float)

    def free_index(self) -> np.ndarray:
        where = {c: i for i, c in enumerate(self.net.coordinates)}
        return np.asarray([where[tuple(k)] for k in self.keys], dtype=np.int64)


def _signature(pot: Mapping[str, Any]) -> Tuple:
    return tuple(sorted((k, np.shape(v)) for k, v in pot.items()))


def _stack(pots: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    return {k: jnp.stack([jnp.asarray(p[k], dtype=jnp.float64) for p in pots]) for k in pots[0]}


class Compiled:
    """A problem's residual and exact Jacobian. The SOLVE is compiled once per pot shape, inline-pH mode and batch
    size: it returns every output row's state (segment ends; every pH node of an inline-pH pot) and, for the
    Jacobian, d state / d x (rows, n, p) by forward mode; the observables (cheap functions of one state) are applied
    after it with their exact gradients, so one compile serves every pot of a shape whatever is observed on it.
    ``residual(x)`` -> (r, failed pots); ``jacobian(x)`` -> (r, J, failed pots). Counts calls and records the first
    call's wall time (compile included).

    SWEEP COUNTS (inline pH). Every call records each pot's fixed-point sweep counts (``sweep_log``: pot -> one tuple
    per call). A derivative is taken at a fixed count; where the count changes between evaluations the residual jumps
    and the Jacobian does not see it, so ``sweep_changes()`` names those pots (``JaxFit.sweep_changes``)."""

    def __init__(self, problem: JaxProblem, **solve_override):
        import diffrax

        self.problem = problem
        net = problem.net
        solve = dict(problem.solve, **solve_override)
        self.solve = solve
        table = {k: f for m in net.modules for k, f in m.observables.items()}
        unknown = sorted({o.observable for o in problem.observations} - set(table))
        if unknown:
            raise ValueError(f"unknown observables {unknown}")
        # (v2.3) an inline-pH pot carries static entries (its "dynamic" / "clamped" mode, node minutes): the solve
        # takes its numeric arrays, and the mode is part of the compile key
        arrays = [pot_arrays(net, p) for p in problem.pots]
        modes = [pot_mode(net, p) for p in problem.pots]
        rows = [n_rows(net, p) for p in problem.pots]
        stray = [(i, o.pot, o.row, rows[o.pot]) for i, o in enumerate(problem.observations)
                 if not 0 <= int(o.row) < rows[o.pot]]
        if stray:
            raise ValueError("observation rows outside their pot's output (observation, pot, row, n_rows): "
                             f"{stray[:5]}")
        groups: Dict[Tuple, List[int]] = {}
        for i, pot in enumerate(arrays):
            groups.setdefault((_signature(pot), modes[i]) if problem.batch == "vmap" else (i,), []).append(i)
        self.groups = list(groups.values())
        self.modes = [modes[g[0]] for g in self.groups]
        if problem.batch == "vmap":
            self.data = [_stack([arrays[i] for i in g]) for g in self.groups]
        else:
            self.data = [arrays[g[0]] for g in self.groups]
        self.observables = {name: jax.jit(lambda y, fn=table[name]: (fn(y, net.index),
                                                                     jax.grad(lambda z: fn(z, net.index))(y)))
                            for name in sorted({o.observable for o in problem.observations})}
        base = jnp.asarray(problem.base_vector())
        free = jnp.asarray(problem.free_index())
        forcing = problem.forcing

        def solve_one(x, pot, adjoint, mode):
            """(ys, failed, sweeps): failed includes an unconverged pH fixed point (``network.integrate``); sweeps
            are the fixed point's counts per segment (zeros for a pH-free network)."""
            theta = theta_from_vector(net, base.at[free].set(x))
            if forcing is None:
                run = integrate(net, theta, pot, adjoint=adjoint, throw=False, mode=mode, **solve)
            else:
                from .forcing import integrate_forced

                run = integrate_forced(net, theta, pot, forcing, adjoint=adjoint, throw=False, **solve)
            sweeps = run.ph.sweeps if getattr(run, "ph", None) is not None else jnp.zeros(run.failed.shape, jnp.int32)
            return run.ys, jnp.any(run.failed), sweeps

        def compile_for(mode):
            def states(x, pot):
                return solve_one(x, pot, diffrax.RecursiveCheckpointAdjoint(), mode)

            def states_and_jac(x, pot):
                def f(x_):
                    ys, failed, sweeps = solve_one(x_, pot, diffrax.ForwardMode(), mode)
                    return ys, (ys, failed, sweeps)

                dys, (ys, failed, sweeps) = jax.jacfwd(f, has_aux=True)(x)
                return ys, dys, failed, sweeps

            if problem.batch == "vmap":
                return jax.jit(jax.vmap(states, in_axes=(None, 0))), jax.jit(jax.vmap(states_and_jac,
                                                                                     in_axes=(None, 0)))
            return jax.jit(states), jax.jit(states_and_jac)

        self._compiled = {mode: compile_for(mode) for mode in dict.fromkeys(self.modes)}
        self.calls = {"residual": 0, "jacobian": 0}
        self.compile_seconds: Dict[str, Optional[float]] = {"residual": None, "jacobian": None}
        self.sweep_log: Dict[int, List[Tuple[int, ...]]] = {i: [] for i in range(len(problem.pots))}

    def _per_pot(self, with_jac: bool, x):
        """{pot: (ys, dys or None, failed)} for every pot; logs each pot's sweep counts."""
        out = {}
        for g, d, mode in zip(self.groups, self.data, self.modes):
            states, states_jac = self._compiled[mode]
            if with_jac:
                ys, dys, failed, sweeps = states_jac(x, d)
            else:
                (ys, failed, sweeps), dys = states(x, d), None
            ys, failed, sweeps = np.asarray(ys), np.asarray(failed), np.asarray(sweeps)
            dys = None if dys is None else np.asarray(dys)
            if self.problem.batch == "vmap":
                for j, p in enumerate(g):
                    out[p] = (ys[j], None if dys is None else dys[j], bool(failed[j]))
                    self.sweep_log[p].append(tuple(int(s) for s in sweeps[j]))
            else:
                out[g[0]] = (ys, dys, bool(failed))
                self.sweep_log[g[0]].append(tuple(int(s) for s in sweeps))
        return out

    def sweep_changes(self) -> Dict[int, List[Tuple[int, ...]]]:
        """Pots whose fixed-point sweep counts differed between the calls so far -> the distinct counts in the order
        seen (empty for pH-free networks): where a fit's residual jumped and its Jacobian could not see it."""
        out = {}
        for p, seen in self.sweep_log.items():
            distinct = list(dict.fromkeys(seen))
            if len(distinct) > 1:
                out[p] = distinct
        return out

    def _assemble(self, states, with_jac: bool):
        obs = self.problem.observations
        p = len(self.problem.keys)
        r = np.zeros(len(obs))
        jac = np.zeros((len(obs), p)) if with_jac else None
        fail_value = float(self.problem.failed_residual)
        for i, o in enumerate(obs):
            ys, dys, failed = states[o.pot]
            if failed:
                r[i] = fail_value
                continue
            f, g = self.observables[o.observable](jnp.asarray(ys[o.row]))
            f, g = float(f), np.asarray(g)
            if o.log10:
                fl = max(f, LOG_FLOOR)
                r[i] = (math.log10(fl) - o.value) / o.sigma
                if with_jac:
                    jac[i] = (g @ dys[o.row]) / (fl * math.log(10.0)) / o.sigma if f > LOG_FLOOR else 0.0
            else:
                r[i] = (f - o.value) / o.sigma
                if with_jac:
                    jac[i] = (g @ dys[o.row]) / o.sigma
        failed = sorted(k for k, (_, _, bad) in states.items() if bad)
        return r, jac, failed

    def residual(self, x):
        t0 = time.perf_counter()
        r, _, failed = self._assemble(self._per_pot(False, jnp.asarray(x, dtype=jnp.float64)), False)
        if self.compile_seconds["residual"] is None:
            self.compile_seconds["residual"] = time.perf_counter() - t0
        self.calls["residual"] += 1
        return r, failed

    def jacobian(self, x):
        t0 = time.perf_counter()
        r, jac, failed = self._assemble(self._per_pot(True, jnp.asarray(x, dtype=jnp.float64)), True)
        if self.compile_seconds["jacobian"] is None:
            self.compile_seconds["jacobian"] = time.perf_counter() - t0
        self.calls["jacobian"] += 1
        return r, jac, failed




@dataclass
class JaxFit:
    theta: np.ndarray
    cost: float
    start: np.ndarray
    start_cost: float
    status: int
    message: str
    nfev: int
    njev: int
    budget: int
    seconds: float
    compile_seconds: Dict[str, Optional[float]]
    residuals: np.ndarray
    jacobian: np.ndarray
    failed_pots: List[int]
    converged: Optional[bool] = None
    #: (inline pH) pots whose fixed-point sweep counts changed between the fit's evaluations -> the distinct counts:
    #: there the residual jumped and the exact Jacobian could not see it (``Compiled.sweep_changes``)
    sweep_changes: Dict[int, List[Tuple[int, ...]]] = field(default_factory=dict)

    def as_dict(self, keys: Sequence[Tuple[str, str]]) -> Dict[str, Any]:
        names = [f"{m}.{k}" for m, k in keys]
        return {"theta": dict(zip(names, map(float, self.theta))), "cost": self.cost, "start_cost": self.start_cost,
                "start": dict(zip(names, map(float, self.start))), "status": self.status, "message": self.message,
                "nfev": self.nfev, "njev": self.njev, "budget": self.budget, "seconds": self.seconds,
                "compile_seconds": self.compile_seconds, "failed_pots": self.failed_pots,
                "converged": self.converged,
                "sweep_changes": {str(k): [list(t) for t in v] for k, v in self.sweep_changes.items()}}


def fit(problem: JaxProblem, start: Sequence[float], budget: int = 200, *, compiled: Optional[Compiled] = None,
        ftol: Optional[float] = None, require_converged: Optional[str] = "raise") -> JaxFit:
    """One least-squares fit with the engine's settings on the exact Jacobian. ``budget`` is least_squares' max_nfev
    (residual evaluations; each iteration adds one Jacobian, counted in ``njev``). ``require_converged``: "raise"
    (default) raises ``fitting.FitNotConverged`` when the status is <= 0 or a pot failed at the optimum; "flag" sets
    ``JaxFit.converged``; None does neither. ``compiled`` reuses a ``Compiled`` (one compile for many starts)."""
    from scipy.optimize import least_squares

    from src.kinetic_core.fitting import LEAST_SQUARES_SETTINGS, FitNotConverged

    if require_converged not in (None, "raise", "flag"):
        raise ValueError(f"require_converged must be None, 'raise' or 'flag', got {require_converged!r}")
    comp = compiled or Compiled(problem)
    lower, upper = np.asarray(problem.lower, dtype=float), np.asarray(problem.upper, dtype=float)
    x0 = np.clip(np.asarray(start, dtype=float), lower, upper)
    cache: Dict[bytes, Tuple[np.ndarray, np.ndarray]] = {}
    count = {"nfev": 0, "njev": 0}

    def fun(x):
        count["nfev"] += 1
        return comp.residual(x)[0]

    def jac(x):
        count["njev"] += 1
        key = np.asarray(x, dtype=float).tobytes()
        if key not in cache:
            cache.clear()
            r, j, _ = comp.jacobian(x)
            cache[key] = (r, j)
        return cache[key][1]

    t0 = time.perf_counter()
    r0 = comp.residual(x0)[0]
    settings = {**LEAST_SQUARES_SETTINGS, **({} if ftol is None else {"ftol": ftol})}
    result = least_squares(fun, x0, jac=jac, bounds=(lower, upper), max_nfev=int(budget), **settings)
    r, j, failed = comp.jacobian(result.x)
    out = JaxFit(theta=np.asarray(result.x, dtype=float), cost=float(0.5 * r @ r), start=x0,
                 start_cost=float(0.5 * r0 @ r0), status=int(result.status), message=str(result.message),
                 nfev=count["nfev"], njev=count["njev"], budget=int(budget),
                 seconds=round(time.perf_counter() - t0, 3), compile_seconds=dict(comp.compile_seconds),
                 residuals=r, jacobian=j, failed_pots=failed, sweep_changes=comp.sweep_changes())
    if require_converged is not None:
        ok = out.status > 0 and not failed
        if require_converged == "raise" and not ok:
            raise FitNotConverged(f"JAX fit not converged: status {out.status} ({out.message}); failed pots {failed}",
                                  out)
        if require_converged == "flag":
            out.converged = bool(ok)
    return out


def laplace(problem: JaxProblem, fit_or_jac, residuals: Optional[np.ndarray] = None, *,
            scale_by_chi2: bool = False) -> Dict[str, Any]:
    """Sigmas from the exact Jacobian at an optimum: J^T J inverted on its eigenbasis (eigenvalues below 1e-12 of the
    largest dropped; a coordinate with weight > 1e-6 in the dropped subspace gets an infinite sigma, as the engine's
    ``fitting.laplace``), optionally scaled by the reduced chi-square."""
    from src.kinetic_core.fitting import NULL_WEIGHT_TOL

    if isinstance(fit_or_jac, JaxFit):
        jac, residuals = fit_or_jac.jacobian, fit_or_jac.residuals
    else:
        jac = np.asarray(fit_or_jac, dtype=float)
    jtj = jac.T @ jac
    w, v = np.linalg.eigh(jtj)
    keep = w > 1e-12 * max(float(w.max()), 1e-300)
    cov = (v[:, keep] / w[keep]) @ v[:, keep].T
    null_weight = np.sum(v[:, ~keep] ** 2, axis=1) if (~keep).any() else np.zeros(len(w))
    sigma = np.sqrt(np.maximum(np.diag(cov), 0.0))
    sigma = np.where(null_weight > NULL_WEIGHT_TOL, np.inf, sigma)
    out = {"sigma": sigma, "covariance": cov, "eigenvalues": w, "dropped": int((~keep).sum())}
    if residuals is not None:
        dof = max(1, len(residuals) - len(problem.keys))
        chi2 = float(np.asarray(residuals) @ np.asarray(residuals)) / dof
        out["reduced_chi2"] = chi2
        if scale_by_chi2:
            out["sigma"] = sigma * math.sqrt(chi2)
            out["covariance"] = cov * chi2
    return out


def fd_check(problem: JaxProblem, x: Sequence[float], compiled: Optional[Compiled] = None,
             move: float = 1e-4) -> Dict[str, Any]:
    """The exact Jacobian against central differences of the same residual, with a step in each COLUMN's own scale
    (tasks/lessons.md 2026-11-05): h_j = move / max_i |J_ij|, so the worst residual moves by ``move`` (a fixed
    absolute step is wrong for a linear k_ref coordinate with no natural unit). One-sided inward at a bound."""
    comp = compiled or Compiled(problem)
    x = np.asarray(x, dtype=float)
    lower, upper = np.asarray(problem.lower, dtype=float), np.asarray(problem.upper, dtype=float)
    _, jac, _ = comp.jacobian(x)
    fd = np.zeros_like(jac)
    steps = []
    for j in range(len(x)):
        h = move / max(float(np.max(np.abs(jac[:, j]))), 1e-12)
        up, down = x[j] + h <= upper[j], x[j] - h >= lower[j]
        xp, xm = x.copy(), x.copy()
        xp[j] += h
        xm[j] -= h
        if up and down:
            fd[:, j] = (comp.residual(xp)[0] - comp.residual(xm)[0]) / (2 * h)
        elif up:
            fd[:, j] = (comp.residual(xp)[0] - comp.residual(x)[0]) / h
        else:
            fd[:, j] = (comp.residual(x)[0] - comp.residual(xm)[0]) / h
        steps.append(h)
    scale = np.maximum(np.max(np.abs(jac), axis=0), 1e-300)
    col = np.max(np.abs(jac - fd), axis=0) / scale
    return {"max_relative_column_error": float(col.max()), "per_column": col.tolist(), "steps": steps,
            "sweep_changes": comp.sweep_changes()}


def residual_tolerance_check(problem: JaxProblem, x: Sequence[float], loose: float = 1e-6,
                             tight: float = 1e-8) -> Dict[str, Any]:
    """Contract sec. 2.6: how far the residuals move (in the observables' own units, dex for log10 rows) between
    rtol ``loose`` and ``tight``; fits may run at ``loose`` when this is below ``FIT_RTOL_CHECK_DEX``."""
    sig = np.asarray([o.sigma for o in problem.observations])
    out = {}
    for name, rtol in (("loose", loose), ("tight", tight)):
        comp = Compiled(problem, rtol=rtol)
        out[name] = comp.residual(np.asarray(x, dtype=float))[0] * sig
    move = float(np.max(np.abs(out["loose"] - out["tight"])))
    return {"rtol_loose": loose, "rtol_tight": tight, "max_move": move, "passes": move < FIT_RTOL_CHECK_DEX}
