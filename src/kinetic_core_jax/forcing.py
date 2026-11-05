"""
src/kinetic_core_jax/forcing.py -- the forced-input driver of the A2 contract (v2.2, 2026-11-05; contract sec. 2.3,
2.4 and the sec. 4 rows "a measured pH trajectory as input" and "the forced-input driver").

A FORCED INPUT is a measured course imposed on a pot instead of being computed by the network:

* a pot SCALAR a rate function reads through ``RateCtx`` -- the temperature ``t_c``, the declared (label-scale) pH
  ``ph_label``, the water activity ``aw``. Forcing ``ph_label`` is the TRAJECTORY form of pH mode "prescribed"
  (M1c's measured pH courses, the measured-trajectory pH operator); forcing ``t_c`` is a thermal profile T(t) with
  the rate constants at T(t) (contract sec. 2.5's profile mode, for a pH-free network);
* a SPECIES level (M1b's ARP-Phe course, the measured 1-DG / 3-DG, GO / MGO courses of M3 / M4 / M6). A forced
  species LEAVES the state: its value in every mass-action term is the interpolant, its derivative is zero, and what
  the network would have produced into (or drawn from) it is integrated into a FORCED-EXCHANGE LEDGER (one extra state
  per forced species), so element closure still reports: E_free . dy_free + E_forced . ledger = 0 exactly.

THE COURSE is data, never a constant traced into the compiled code: knots t_0 = 0 < t_1 < ... < t_K (minutes) and the
values at the knots, interpolated per name either "linear" (piecewise-linear) or "previous" (a staircase: the value at
the interval's left knot). Each interval between knots is one SEGMENT of ``network.integrate``'s design: integrated on
a unit interval s in [0, 1] with the right-hand side scaled by its duration, the inputs at s being a + s (b - a). So a
knot is an integrator restart (``jump_ts`` by construction), one compile serves every course with the same number of
intervals, zero-duration intervals pad courses to a common length, and the whole thing runs under ``vmap``.

DERIVATIVES are exact through the solve (forward mode or the discrete adjoint, as ``network.jacobian``), in theta and
also in the knot values (they are arguments). A module that declares a BREAKPOINT of an input (``Module.breakpoints``:
a value where its rate function jumps or kinks -- the lipid breakage window at 100 C, B18's pyrazine knot at pH 7) gets
a knot inserted where a linear course crosses it (``forced_pot``), so the jump sits on a segment boundary: the solve
is exact there and, because the crossing time is a function of the DATA only, the theta-Jacobian has no missing
crossing term (contract sec. 2.4's one-sided caveat applies to an inline pH that depends on theta, not to a forced
one). A derivative in the KNOT VALUES does miss the crossing term at such a knot; ``ForcedPot.crossings`` lists them.

Refused (``UnsupportedConfiguration`` / ValueError, never a different model in silence): a network with an inline pH
module (the sulfur port's charge balance; forcing pH there is a different mode, "clamped"), forcing an input no module
reads, an in-situ pH scale (every ported module reads the label-scale declared pH), a forced species the network does
not carry, a course that is not finite or whose knots do not increase.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Mapping, NamedTuple, Optional, Sequence, Tuple

import jax
import jax.numpy as jnp
import numpy as np

from .lane import UnsupportedConfiguration
from .network import CELSIUS, Network, RateCtx, mass_action, theta_from_vector

#: the pot scalars a course may force (``RateCtx`` fields); "t_c" is always read
FORCEABLE = ("t_c", "ph_label", "aw")
INTERPOLATIONS = ("linear", "previous")
#: predictive_densities.PH_SCALES; "in_situ" is refused (no ported module reads an in-situ pH)
PH_SCALES = ("cooled", "label_25C", "in_situ")


@dataclass(frozen=True)
class Forcing:
    """The STATIC part of a forcing (what changes the compiled solve): which scalars and species are forced and how
    each is interpolated. The course itself (knots and values) is pot data (``forced_pot``).

    ``interpolation`` is a tuple of (name, "linear" | "previous") pairs; a name not listed is "linear".
    ``ph_scale`` declares the scale of a forced pH (M1b sec. 3.3, V5): "cooled" or "label_25C", both read as the
    label-scale declared pH ``RateCtx.ph_label``; "in_situ" is refused."""

    inputs: Tuple[str, ...] = ("t_c",)
    species: Tuple[str, ...] = ()
    interpolation: Tuple[Tuple[str, str], ...] = ()
    ph_scale: Optional[str] = None

    def __post_init__(self):
        if "t_c" not in self.inputs:
            raise ValueError("a forcing must carry the temperature course 't_c' (a thermal program is a 'previous' "
                             "course: ``program_course``)")
        bad = [x for x in self.inputs if x not in FORCEABLE]
        if bad:
            raise ValueError(f"cannot force {bad}: forceable scalars are {FORCEABLE}")
        if len(set(self.names)) != len(self.names):
            raise ValueError(f"a name is forced twice: {self.names}")
        for name, how in self.interpolation:
            if how not in INTERPOLATIONS:
                raise ValueError(f"{name}: interpolation {how!r} is not one of {INTERPOLATIONS}")
            if name not in self.names:
                raise ValueError(f"interpolation given for {name!r}, which is not forced")
        if "ph_label" in self.inputs:
            if self.ph_scale is None or self.ph_scale not in PH_SCALES:
                raise ValueError(f"a forced pH declares its scale, one of {PH_SCALES} (got {self.ph_scale!r})")
            if self.ph_scale == "in_situ":
                raise UnsupportedConfiguration("an in-situ pH course: every ported module reads the label-scale "
                                               "declared pH; an in-situ pH enters with the sulfur port's inline mode")

    @property
    def names(self) -> Tuple[str, ...]:
        return tuple(self.inputs) + tuple(self.species)

    def how(self, name: str) -> str:
        return dict(self.interpolation).get(name, "linear")

    @property
    def inputs_vary_within_segments(self) -> bool:
        """True when some forced scalar is linear: the rate constants are then evaluated in the right-hand side."""
        return any(self.how(x) == "linear" for x in self.inputs)


def ph_mode(net: Network, forcing: Optional[Forcing]) -> str:
    """``Network.ph_mode`` under a forcing: "prescribed (trajectory)" when the pH is a forced course."""
    if forcing is not None and "ph_label" in forcing.inputs:
        return "prescribed (trajectory)"
    return net.ph_mode


def check(net: Network, forcing: Forcing) -> None:
    """Refuse a forcing this network cannot take (module docstring)."""
    if net.needs_ph:
        raise UnsupportedConfiguration("an inline-pH network (contract sec. 2.4) is not ported; a forced pH there is "
                                       "the 'clamped' mode, which comes with the sulfur port")
    unread = [x for x in forcing.inputs if x != "t_c" and x not in net.pot_inputs]
    if unread:
        raise ValueError(f"forcing {unread}, which no module of the network reads (pot inputs {net.pot_inputs})")
    unknown = [s for s in forcing.species if s not in net.index]
    if unknown:
        raise UnsupportedConfiguration(f"forced species {unknown} are not carried by the network")


def breakpoints(net: Network, name: str) -> Tuple[float, ...]:
    """Every module's declared breakpoints of input ``name`` (sorted, unique)."""
    values = set()
    for m in net.modules:
        values.update(float(v) for v in dict(getattr(m, "breakpoints", {}) or {}).get(name, ()))
    return tuple(sorted(values))


# --------------------------------------------------------------------------- courses and pots
def program_course(segments: Sequence[Tuple[float, float]]) -> Tuple[np.ndarray, np.ndarray]:
    """A thermal program ``((duration_min, T_C), ...)`` as a 'previous' temperature course: (knots, values). The
    Forcing must declare ``("t_c", "previous")``: read as linear, the same knots are a ramp between the segments."""
    rows = [(float(d), float(t)) for d, t in segments]
    knots = np.concatenate([[0.0], np.cumsum([d for d, _ in rows])])
    values = np.asarray([t for _, t in rows] + [rows[-1][1]], dtype=float)
    return knots, values


def _value_at(knots: np.ndarray, values: np.ndarray, t: float, how: str) -> float:
    if how == "linear":
        return float(np.interp(t, knots, values))
    i = int(np.searchsorted(knots, t, side="right") - 1)
    return float(values[min(max(i, 0), len(values) - 1)])


class ForcedPot(NamedTuple):
    pot: Dict[str, Any]                 # the arrays ``integrate_forced`` takes
    rows: List[int]                     # rows[i]: the output row holding the state at save_times[i]
    knots: np.ndarray                   # the refined knots (minutes), K + 1 before padding
    crossings: List[Tuple[str, float, float]]   # (input, breakpoint, time) where a knot was inserted at a crossing


def forced_pot(net: Network, forcing: Forcing, y0: Sequence[float], courses: Mapping[str, Tuple[Sequence[float],
               Sequence[float]]], *, scalars: Optional[Mapping[str, float]] = None, save_times: Sequence[float] = (),
               n_segments: Optional[int] = None) -> ForcedPot:
    """A forced pot. ``courses[name] = (knots_min, values)`` for every name in ``forcing.names``; each course may have
    its own knots (they are merged; a course is held at its last value past its last knot and must start at 0).
    ``scalars`` carries the network's pot inputs that are NOT forced (NaN for the engine's None). Every save time and
    every crossing of a declared breakpoint by a linear course becomes a knot. ``n_segments`` pads with zero-duration
    intervals (the identity) so courses of different lengths share one compile."""
    check(net, forcing)
    missing = [n for n in forcing.names if n not in courses]
    if missing:
        raise ValueError(f"no course for forced {missing}")
    table: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for name in forcing.names:
        t, v = (np.asarray(x, dtype=float) for x in courses[name])
        if t.ndim != 1 or t.shape != v.shape or t.size < 1:
            raise ValueError(f"{name}: knots and values must be equal-length vectors")
        if not (np.all(np.isfinite(t)) and np.all(np.isfinite(v))):
            raise ValueError(f"{name}: a course must be finite (NaN is a pot scalar's None, not a course value)")
        if t[0] != 0.0 or np.any(np.diff(t) <= 0):
            raise ValueError(f"{name}: knots must start at 0 and increase strictly")
        if name in forcing.species and np.any(v < 0):
            raise ValueError(f"{name}: a forced species level must be >= 0")
        table[name] = (t, v)
    if any(m.name == "acrylamide" for m in net.modules):
        from .acrylamide import check_programme      # CF-106: the engine refuses a peak inside its window

        check_programme([(1.0, float(np.max(table["t_c"][1])))])
    end = max(float(t[-1]) for t, _ in table.values())
    cuts = {0.0, end}
    for t, _ in table.values():
        cuts.update(float(x) for x in t)
    save = [float(s) for s in save_times]
    if save and (min(save) <= 0.0 or max(save) > end * (1 + 1e-12)):
        raise ValueError(f"save times must lie in (0, {end}]")
    cuts.update(save)
    crossings: List[Tuple[str, float, float]] = []
    for name in forcing.inputs:
        if forcing.how(name) != "linear":
            continue
        t, v = table[name]
        for b in breakpoints(net, name):
            for (t0, t1, v0, v1) in zip(t[:-1], t[1:], v[:-1], v[1:]):
                if (v0 - b) * (v1 - b) < 0.0:
                    tc = float(t0 + (b - v0) / (v1 - v0) * (t1 - t0))
                    crossings.append((name, b, tc))
                    cuts.add(tc)
    knots = np.asarray(sorted(cuts), dtype=float)
    keep = np.concatenate([[True], np.diff(knots) > 1e-12 * max(1.0, end)])
    knots = knots[keep]
    if len(knots) < 2:
        raise ValueError("a forced pot needs a positive duration")
    k = len(knots) - 1
    start = np.zeros((k, len(forcing.names)))
    stop = np.zeros((k, len(forcing.names)))
    for j, name in enumerate(forcing.names):
        t, v = table[name]
        how = forcing.how(name)
        for i in range(k):
            a = _value_at(t, v, knots[i], how)
            start[i, j] = a
            stop[i, j] = _value_at(t, v, knots[i + 1], how) if how == "linear" else a
    durations = np.diff(knots)
    if n_segments is not None:
        if n_segments < k:
            raise ValueError(f"{k} intervals do not fit in {n_segments}")
        pad = n_segments - k
        durations = np.concatenate([durations, np.zeros(pad)])
        last = stop[-1:].repeat(pad, axis=0)
        start, stop = np.concatenate([start, last]), np.concatenate([stop, last])
    rows = [int(np.argmin(np.abs(knots[1:] - s))) for s in save]
    pot: Dict[str, Any] = {"y0": np.asarray(y0, dtype=float), "durations": durations, "start": start, "stop": stop}
    for x in net.pot_inputs:
        if x in forcing.inputs:
            continue
        if scalars is None or x not in scalars:
            raise ValueError(f"the pot lacks the scalar {x!r}, which the network reads and the forcing does not carry")
        pot[x] = float(scalars[x])
    return ForcedPot(pot=pot, rows=rows, knots=knots, crossings=crossings)


# --------------------------------------------------------------------------- the solve
class ForcedRun(NamedTuple):
    ys: Any          # (S, n_species): the state at the END of each interval; forced species read their course
    exchange: Any    # (S, n_forced): the forced-exchange ledger (what the network produced into each forced species)
    num_steps: Any   # (S,)
    failed: Any      # (S,)


def integrate_forced(net: Network, theta, pot: Mapping[str, Any], forcing: Forcing, *, solver=None,
                     rtol: float = 1e-8, atol: float = 1e-14, adjoint=None, max_steps: int = 100_000,
                     throw: bool = True) -> ForcedRun:
    """One forced pot (``forced_pot(...).pot``) through its course. With every forced scalar a staircase the rate
    constants are evaluated once per interval (as ``network.integrate``); with a linear one, in the right-hand side."""
    import diffrax

    check(net, forcing)
    n = net.n_species
    m_in = len(forcing.inputs)
    pos = {x: i for i, x in enumerate(forcing.inputs)}
    fidx = np.asarray([net.index[s] for s in forcing.species], dtype=np.int64)
    n_f = len(fidx)
    fixed = {x: jnp.asarray(pot[x], dtype=jnp.float64) for x in net.pot_inputs if x not in pos}
    stoich = jnp.asarray(net.stoichiometry)
    varying = forcing.inputs_vary_within_segments

    def ctx_of(u):
        kw = dict(fixed)
        for x in ("ph_label", "aw"):
            if x in pos:
                kw[x] = u[pos[x]]
        t_c = u[pos["t_c"]]
        return RateCtx(t_c=t_c, t_k=t_c + CELSIUS, **{x: kw[x] for x in net.pot_inputs})

    def rhs(s, z, args):
        k, duration, a, b = args
        u = a + s * (b - a)
        if varying:
            k = net.rate_constants(theta, ctx_of(u[:m_in]))
        y = z[:n]
        if n_f:
            y = y.at[fidx].set(u[m_in:])
        rate = stoich @ (k * mass_action(y, net))
        if n_f:
            return duration * jnp.concatenate([rate.at[fidx].set(0.0), rate[fidx]])
        return duration * rate

    term = diffrax.ODETerm(rhs)
    solver = solver or getattr(diffrax, net.solver)()
    adjoint = adjoint or diffrax.RecursiveCheckpointAdjoint()
    controller = diffrax.PIDController(rtol=rtol, atol=atol)

    def segment(z, seg):
        duration, a, b = seg
        k = jnp.zeros(len(net.reactions)) if varying else net.rate_constants(theta, ctx_of(a[:m_in]))
        sol = diffrax.diffeqsolve(term, solver, t0=0.0, t1=1.0, dt0=None, y0=z, args=(k, duration, a, b),
                                  saveat=diffrax.SaveAt(t1=True), stepsize_controller=controller,
                                  max_steps=max_steps, adjoint=adjoint, throw=throw)
        z1 = sol.ys[-1]
        return z1, (z1, sol.stats["num_steps"], sol.result != diffrax.RESULTS.successful)

    start = jnp.asarray(pot["start"], dtype=jnp.float64)
    stop = jnp.asarray(pot["stop"], dtype=jnp.float64)
    y0 = jnp.asarray(pot["y0"], dtype=jnp.float64)
    if n_f:
        y0 = y0.at[fidx].set(start[0, m_in:])
        y0 = jnp.concatenate([y0, jnp.zeros(n_f)])
    _, (zs, steps, failed) = jax.lax.scan(
        segment, y0, (jnp.asarray(pot["durations"], dtype=jnp.float64), start, stop))
    ys = zs[:, :n]
    if n_f:
        ys = ys.at[:, fidx].set(stop[:, m_in:])
        exchange = zs[:, n:]
    else:
        exchange = jnp.zeros((zs.shape[0], 0))
    return ForcedRun(ys=ys, exchange=exchange, num_steps=steps, failed=failed)


def closure(net: Network, forcing: Forcing, y0, run: ForcedRun) -> np.ndarray:
    """Per interval end, E_free . (y - y0)_free + E_forced . ledger, relative to E . |y0| per element (C, N, S):
    zero to the solver's tolerance when the forced-exchange ledger accounts for every forced flux."""
    y0 = np.asarray(y0, dtype=float)
    ys, ledger = np.asarray(run.ys), np.asarray(run.exchange)
    e = np.asarray(net.element_matrix)
    fidx = [net.index[s] for s in forcing.species]
    free = np.ones(net.n_species, dtype=bool)
    free[fidx] = False
    delta = (ys - y0)[:, free] @ e[:, free].T
    if fidx:
        delta = delta + ledger @ e[:, fidx].T
    scale = np.maximum(e @ np.abs(y0), 1e-300)
    return delta / scale


def jacobian_forced(net: Network, forcing: Forcing, observe, **solve_kw):
    """``f(theta_vector, pot) -> (observables, d observables / d theta)``: ``network.jacobian`` for a forced pot."""
    import diffrax

    solve_kw = dict(solve_kw, adjoint=diffrax.ForwardMode())

    def obs(vector, pot):
        value = observe(integrate_forced(net, theta_from_vector(net, vector), pot, forcing, **solve_kw).ys)
        return value, value

    def both(vector, pot):
        jac, value = jax.jacfwd(obs, has_aux=True)(vector, pot)
        return value, jac

    return jax.jit(both)


def batched_forced(net: Network, forcing: Forcing, *, theta_axis: Optional[int] = 0, pot_axis: Optional[int] = 0,
                   **solve_kw):
    """``network.batched`` for forced pots (one compile per (network, forcing, S, B))."""

    def one(theta, pot):
        return integrate_forced(net, theta, pot, forcing, **solve_kw).ys

    return jax.jit(jax.vmap(one, in_axes=(theta_axis, pot_axis)))


def staircase_course(knots: Sequence[float], values: Sequence[float], n_sub: int,
                     rule: str = "midpoint") -> Tuple[np.ndarray, np.ndarray]:
    """A linear course as a staircase of ``n_sub`` steps per interval, each at the course's value at the step's
    midpoint (or left end): the engine-runnable approximation of a linear course (it converges as n_sub^-2 for the
    midpoint rule). Returns (knots, values) for a 'previous' course."""
    t, v = np.asarray(knots, dtype=float), np.asarray(values, dtype=float)
    out_t: List[float] = [0.0]
    out_v: List[float] = []
    for t0, t1 in zip(t[:-1], t[1:]):
        edges = np.linspace(t0, t1, n_sub + 1)
        for a, b in zip(edges[:-1], edges[1:]):
            at = 0.5 * (a + b) if rule == "midpoint" else a
            out_v.append(float(np.interp(at, t, v)))
            out_t.append(float(b))
    out_v.append(out_v[-1])
    return np.asarray(out_t), np.asarray(out_v)
