"""
src/kinetic_core_jax/ph_coupling.py -- the network's OWN pH solve (A2 contract v2.3 sec. 2.4; plan step A2 + H3,
2026-11-06): the merged titratable table as arrays (``ChargeTable``), the charge balance and its (pH, I) root with the
implicit-function derivative (``PhSolver``), the ONE charge-closure check over the merged table (``closure``), the pH
pot (``pot_from_program``: node grids, temperature knots, the jump grid) and the integrator of an inline-pH network
(``integrate``: the engine's two couplings, "fixed_point_nodes" and "inline", and the clamped pH).

``network.integrate`` calls ``integrate`` here for any network with a ``needs_ph`` module; nothing here knows which
modules are composed. A module brings to the solve what it declares on ``Module`` (contract sec. 2.2):

* ``titratables`` -- ``network.Titratable`` pools: pK / dH / dCp at 25 C per protonation step, the charge of each
  state, mmol of pool per mmol of each species (``weights``), a YIELD coupling (``yield_weights`` scaled by one of the
  module's coordinates: the lumped acid pool x ``acid_yield``), a pK that IS a coordinate (``pka_coordinates``: the
  ARP amine pK), and pot-supplied pools (``source="pot"``: the buffer, mmol/L read from the pot's ``buffer`` vector);
* ``ph_model`` -- the solvent (pKw, water density, the Davies constant's inputs, the label and pK reference
  temperatures, R) and the coupling's numerics (``PhModel``), declared by exactly one module (or equal ones);
* ``centre_ledger`` -- per reaction, the titratable centres it moves and why (the engine's ``CENTRE_LEDGER``), which
  ``closure`` checks against the movement the merged table itself implies.

Moved here from ``sulfur_lane.py`` (2026-11-05, PR #83), where it was the sulfur lane's private driver because the
titratable table could not carry the inventory's dH / dCp / acid yield. The arithmetic is unchanged, operation for
operation, so the sulfur lane's parity (tests/scientific/test_jax_sulfur_parity.py, 1e-4 dex) carries over.

DIFFERENTIATION. The root's derivative is the implicit function (a 2 x 2 solve; forward and reverse mode). The fixed
point's sweep loop is chosen by the adjoint: under ``diffrax.ForwardMode()`` a ``lax.while_loop`` (``jax.jacfwd``,
``vmap``); under any reverse-capable adjoint (the default, ``RecursiveCheckpointAdjoint``) equinox's CHECKPOINTED
while loop, which ``jax.grad`` differentiates and which runs only the sweeps it needs (PR #83's blocker 4: ``jax.grad``
needed ``reverse=True``; it no longer does). ``reverse=True`` keeps PR #83's scan over ``max_sweeps`` (finished
sweeps skipped by ``lax.cond``). The sweep count is the engine's discrete stopping rule: a perturbation that changes
it is a jump no derivative sees (``PhRun.sweeps`` reports it; contract sec. 2.4).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, NamedTuple, Optional, Sequence, Tuple

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)

from .lane import UnsupportedConfiguration  # noqa: E402

CELSIUS = 273.15
LN10 = math.log(10.0)
PAD_PK = 1.0e3          # a padded protonation step that never ionises
GROUP_WIDTH = 4
BISECTION_STEPS = 14    # (-1, 15) to 1e-3 pH, then Newton (the bake-off's solve)
NEWTON_STEPS = 6
SID_SWEEPS = 60
COUPLINGS = ("fixed_point_nodes", "inline")
NODE_GRIDS = ("uniform", "log_and_linear")
FIXED_POINTS = ("relaxed_converged", "fixed_sweeps")
#: integrate_sulfur's stopping threshold under "fixed_sweeps" (a literal in its loop; pinned through its AST)
FIXED_SWEEPS_STOP = 1e-4
#: the numeric pot keys an inline-pH network reads (``pot["mode"]``, "dynamic" | "clamped", is static)
POT_KEYS = ("y0", "label_ph", "buffer", "durations", "node_s", "node_tc", "knot_s", "knot_tc", "label_nudge", "jump_s")
MODES = ("dynamic", "clamped")
#: nodes closer than this (as a fraction of the segment) are one kink of the rate vector for the step controller
JUMP_MERGE = 1e-12
#: the closure check's tolerance on a centre movement (pool weights are shares, e.g. 0.5 + 0.5)
CLOSURE_TOL = 1e-9


# --------------------------------------------------------------------------- the declarations
@dataclass(frozen=True)
class PhModel:
    """The solvent and the coupling's numerics of an inline-pH network (one per network). The solvent tables are the
    engine's (``ph_state``: R, the pK reference and label temperatures, pKw(T), the water density); the numerics are
    the Config's (``ph.coupling``, ``ph.node_grid``, ``ph.fixed_point*``, ``sulfur.profile_node_step_k``) and
    ``integrate_sulfur``'s own defaults (``ph_nodes``, ``ph_iterations``)."""

    coupling: str
    node_grid: str
    fixed_point: str
    relaxation: float
    tolerance: float
    max_sweeps: int
    flag: float
    profile_node_step_k: float
    ph_nodes: int
    ph_iterations: int
    r_kj: float
    pka_ref_t_k: float
    label_t_k: float
    pkw_table_k: Tuple[float, ...]
    pkw_table: Tuple[float, ...]
    rho_table_c: Tuple[float, ...]
    rho_table: Tuple[float, ...]
    basis: str = ""

    def __post_init__(self):
        for value, allowed in ((self.coupling, COUPLINGS), (self.node_grid, NODE_GRIDS),
                               (self.fixed_point, FIXED_POINTS)):
            if value not in allowed:
                raise UnsupportedConfiguration(f"{value!r}: not one of {allowed}")

    @property
    def const(self) -> Dict[str, Any]:
        """The solvent constants as the charge balance's functions read them."""
        return {"R_KJ": self.r_kj, "pka_ref_t_k": self.pka_ref_t_k, "label_t_k": self.label_t_k,
                "pkw_table_k": self.pkw_table_k, "pkw_table": self.pkw_table, "rho_table_c": self.rho_table_c,
                "rho_table": self.rho_table}


@dataclass(frozen=True)
class ChargeTable:
    """The merged titratable table as arrays (built by ``network.compose`` from every module's ``titratables``, in
    module order). Per pool (row): pK / dH / dCp at 25 C padded to GROUP_WIDTH steps, z_max (the charge of the fully
    protonated state), mmol/L of pool per mmol/L of each species (``weights``), the yield couplings (coordinate ->
    per-unit weights), the pK coordinates (coordinate -> (row, step) mask), the pot-supplied rows, and the
    titratable centres per unit of pool (``centres``: kinds x pools) the closure check reads."""

    names: Tuple[str, ...]
    pk25: np.ndarray
    dh: np.ndarray
    dcp: np.ndarray
    zmax: np.ndarray
    weights: np.ndarray
    yields: Tuple[Tuple[Tuple[str, str], np.ndarray], ...]
    pka_slots: Tuple[Tuple[Tuple[str, str], np.ndarray], ...]
    buffer_rows: Mapping[str, int]
    kinds: Tuple[str, ...]
    centres: np.ndarray
    owners: Tuple[str, ...]

    @property
    def n_pools(self) -> int:
        return len(self.names)

    def buffer_vector(self, spec) -> np.ndarray:
        """The pot's ``buffer`` vector (mmol/L per pool) from a buffer spec's ``groups()`` ((group, mmol/L) pairs whose
        ``group.key`` names a pot-supplied pool)."""
        out = np.zeros(self.n_pools)
        for group, mmol in spec.groups():
            if group.key not in self.buffer_rows:
                raise UnsupportedConfiguration(f"buffer group {group.key!r}: no pot-supplied pool declares it")
            out[self.buffer_rows[group.key]] += float(mmol)
        return out

    def structural_weights(self) -> np.ndarray:
        """mmol of pool per mmol of species with every yield coordinate at 1 (the ledger's view: a lumped acid
        equivalent is one centre, whatever share of it the fitted yield lets titrate)."""
        w = np.array(self.weights, dtype=float)
        for _, wy in self.yields:
            w = w + wy
        return w


def build_table(modules: Sequence[Any], index: Mapping[str, int]) -> Optional[ChargeTable]:
    """The ChargeTable of composed modules (None when no module declares a pool). ``network.compose`` has already
    checked reach, duplicates and charge lists."""
    pools = [(m, t) for m in modules for t in m.titratables]
    if not pools:
        return None
    n, g = len(index), len(pools)
    pk = np.full((g, GROUP_WIDTH), PAD_PK)
    dh = np.zeros((g, GROUP_WIDTH))
    dcp = np.zeros((g, GROUP_WIDTH))
    zmax = np.zeros(g)
    w = np.zeros((g, n))
    yields: Dict[Tuple[str, str], np.ndarray] = {}
    slots: Dict[Tuple[str, str], np.ndarray] = {}
    buffer_rows: Dict[str, int] = {}
    kinds = sorted({k for _, t in pools for k in (t.centres or {"titratable": 1.0})})
    centres = np.zeros((len(kinds), g))
    for r, (m, t) in enumerate(pools):
        k = len(t.pka_25c)
        if k > GROUP_WIDTH:
            raise UnsupportedConfiguration(f"pool {t.species}: {k} protonation steps (at most {GROUP_WIDTH})")
        pk[r, :k] = t.pka_25c
        dh[r, :k] = t.dh_kj_mol or (0.0,) * k
        dcp[r, :k] = t.dcp_j_mol_k or (0.0,) * k
        zmax[r] = float(t.charges[0])
        for s, x in t.pool_weights().items():
            w[r, index[s]] += float(x)
        if t.yield_coordinate is not None:
            key = (m.name, t.yield_coordinate)
            wy = yields.setdefault(key, np.zeros((g, n)))
            for s, x in (t.yield_weights or {}).items():
                wy[r, index[s]] += float(x)
        for step, coord in enumerate(t.pka_coordinates or ()):
            if coord is not None:
                slots.setdefault((m.name, coord), np.zeros((g, GROUP_WIDTH), dtype=bool))[r, step] = True
        if t.source == "pot":
            buffer_rows[t.pot_key] = r
        for kind, count in (t.centres or {"titratable": 1.0}).items():
            centres[kinds.index(kind), r] = float(count)
    return ChargeTable(names=tuple(t.species for _, t in pools), pk25=pk, dh=dh, dcp=dcp, zmax=zmax, weights=w,
                       yields=tuple(yields.items()), pka_slots=tuple(slots.items()), buffer_rows=buffer_rows,
                       kinds=tuple(kinds), centres=centres, owners=tuple(m.name for m, _ in pools))


# --------------------------------------------------------------------------- the charge balance (JAX)
def pk_at(pk25, dh, dcp, t_k, r_kj, t0):
    """ph_state.pk_at with B51's heat-capacity term."""
    out = pk25 + dh / (r_kj * LN10) * (1.0 / t_k - 1.0 / t0)
    return out - (1e-3 * dcp / (r_kj * LN10)) * (jnp.log(t_k / t0) + t0 / t_k - 1.0)


def davies_a(t_k, c):
    t_c = t_k - CELSIUS
    eps = 87.740 - 0.40008 * t_c + 9.398e-4 * t_c ** 2 - 1.410e-6 * t_c ** 3
    rho = jnp.interp(t_c, jnp.asarray(c["rho_table_c"]), jnp.asarray(c["rho_table"]))
    return 1.82483e6 * jnp.sqrt(rho) / (eps * t_k) ** 1.5


def davies(ionic_m, a):
    safe = jnp.maximum(ionic_m, 1e-300)
    s = jnp.sqrt(safe)
    return jnp.where(ionic_m > 0.0, a * (s / (1.0 + s) - 0.3 * safe), 0.0)


def pkw_at(t_k, c):
    x = 1.0 / jnp.asarray(c["pkw_table_k"])[::-1]
    return jnp.interp(1.0 / t_k, x, jnp.asarray(c["pkw_table"])[::-1])


class PhSolver:
    """(pH, I) of a network state at (T, SID): bisection at I = 0, one ionic sweep, Newton on the coupled (pH, I); the
    derivative by the implicit function (a 2 x 2 solve), so forward and reverse mode both see the exact root. ``theta``
    is the network's full pytree ({module: {coordinate: value}}): the yield and pK coordinates are read from their
    declaring module's entry."""

    def __init__(self, table: ChargeTable, model: PhModel):
        self.b = table
        # numpy only: a jnp array made here (inside a trace) would be a tracer the custom JVP closes over
        self.zmax = np.asarray(table.zmax, dtype=float)
        c = model.const

        @jax.custom_jvp
        def solve_tt(tt, sid, conc):
            return self._forward(jax.lax.stop_gradient(tt), jax.lax.stop_gradient(sid), jax.lax.stop_gradient(conc))

        @solve_tt.defjvp
        def _jvp(primals, tangents):
            tt, sid, conc = primals
            u = solve_tt(tt, sid, conc)
            jac_u = jax.jacfwd(lambda uu: self._equations(uu, tt, sid, conc))(u)
            _, d_f = jax.jvp(lambda t, s, cc: self._equations(u, t, s, cc), (tt, sid, conc), tangents)
            return u, -jnp.linalg.solve(jac_u, d_f)

        self._solve_tt = solve_tt
        self.c = c

    def pk25(self, theta):
        out = jnp.asarray(self.b.pk25)
        for (module, coord), mask in self.b.pka_slots:
            out = jnp.where(jnp.asarray(mask), theta[module][coord], out)
        return out

    def temperature_terms(self, t_k, pk25):
        c = self.c
        pk0 = pk_at(pk25, jnp.asarray(self.b.dh), jnp.asarray(self.b.dcp), t_k, c["R_KJ"], c["pka_ref_t_k"])
        return {"pk0": pk0, "pkw": pkw_at(t_k, c), "a": davies_a(t_k, c)}

    def _group_pk(self, ionic_m, tt):
        af = davies(ionic_m, tt["a"])
        return tt["pk0"] + af * (2.0 * (jnp.asarray(self.zmax)[:, None] - jnp.arange(float(GROUP_WIDTH))) - 1.0), af

    def group_terms(self, ph, pk):
        """Per pool: its mean charge and its ionic-strength term at ``ph`` (pK already at T and I)."""
        log_w = jnp.concatenate([jnp.zeros((pk.shape[0], 1)), jnp.cumsum(ph - pk, axis=1)], axis=1)
        w = 10.0 ** (log_w - jnp.max(log_w, axis=1, keepdims=True))
        w = w / jnp.sum(w, axis=1, keepdims=True)
        i = jnp.arange(GROUP_WIDTH + 1.0)
        z = jnp.asarray(self.zmax)
        return z - jnp.sum(i * w, axis=1), jnp.sum(w * (z[:, None] - i) ** 2, axis=1)

    def _residual_given(self, ph, pk, af, pkw, sid, conc):
        charge, ionic = self.group_terms(ph, pk)
        inv_g = 10.0 ** af
        h = 1000.0 * 10.0 ** (-ph) * inv_g
        oh = 1000.0 * 10.0 ** (ph - pkw) * inv_g
        res = sid + h - oh + jnp.sum(conc * charge)
        i_new = 0.5 * (jnp.abs(sid) + h + oh + jnp.sum(conc * ionic)) / 1000.0
        return res, i_new

    def residual(self, ph, ionic_m, tt, sid, conc):
        """SID + [H+] - [OH-] + SUM c z (mmol/L) and the ionic strength it implies (mol/L)."""
        pk, af = self._group_pk(ionic_m, tt)
        return self._residual_given(ph, pk, af, tt["pkw"], sid, conc)

    def _equations(self, u, tt, sid, conc):
        res, i_new = self.residual(u[0], u[1], tt, sid, conc)
        return jnp.stack([res, u[1] - i_new])

    def _forward(self, tt, sid, conc):
        pk, af = self._group_pk(0.0, tt)

        def bisect(_, bounds):
            lo, hi = bounds
            mid = 0.5 * (lo + hi)
            up = self._residual_given(mid, pk, af, tt["pkw"], sid, conc)[0] > 0.0
            return jnp.where(up, mid, lo), jnp.where(up, hi, mid)

        lo, hi = jax.lax.fori_loop(0, BISECTION_STEPS, bisect, (jnp.asarray(-1.0), jnp.asarray(15.0)))
        ph0 = 0.5 * (lo + hi)
        u = jnp.stack([ph0, self.residual(ph0, 0.0, tt, sid, conc)[1]])

        def newton(_, u):
            f = self._equations(u, tt, sid, conc)
            jac = jax.jacfwd(lambda uu: self._equations(uu, tt, sid, conc))(u)
            step = jnp.linalg.solve(jac, f)
            return jnp.stack([u[0] - jnp.clip(step[0], -1.0, 1.0), jnp.maximum(u[1] - step[1], 0.1 * u[1])])

        return jax.lax.fori_loop(0, NEWTON_STEPS, newton, u)

    def solve(self, t_k, sid, conc, pk25):
        """(pH, I) at ``t_k``."""
        return self._solve_tt(self.temperature_terms(t_k, pk25), sid, conc)

    def pool_weights(self, theta):
        w = jnp.asarray(self.b.weights)
        for (module, coord), wy in self.b.yields:
            w = w + theta[module][coord] * jnp.asarray(wy)
        return w

    def inventory(self, y, theta, buffer):
        return self.pool_weights(theta) @ jnp.maximum(y, 0.0) + buffer

    def state(self, y, t_k, sid, theta, buffer):
        """(in-situ pH, I, cooled label pH) of a state: the engine's ``_node_state``."""
        conc = self.inventory(y, theta, buffer)
        pk25 = self.pk25(theta)
        hot = self.solve(t_k, sid, conc, pk25)
        cool = self.solve(jnp.asarray(self.c["label_t_k"]), sid, conc, pk25)
        return hot[0], hot[1], cool[0]

    def strong_ion_difference(self, label_ph, y0, theta, buffer):
        """ph_state.strong_ion_difference: the SID that balances the declared label pH at 25 C, I solved with it."""
        conc = self.inventory(y0, theta, buffer)
        tt = self.temperature_terms(jnp.asarray(self.c["label_t_k"]), self.pk25(theta))
        sid = -self.residual(label_ph, 0.0, tt, 0.0, conc)[0]

        def body(_, state):
            sid, _ionic = state
            fresh = self.residual(label_ph, _ionic, tt, sid, conc)[1]
            return -self.residual(label_ph, fresh, tt, 0.0, conc)[0], fresh

        return jax.lax.fori_loop(0, SID_SWEEPS, body, (sid, jnp.asarray(0.0)))[0]


def solver_of(net) -> PhSolver:
    if net.charge is None or net.ph_model is None:
        raise UnsupportedConfiguration(
            "an inline-pH network needs a titratable table and one PhModel (Module.titratables / Module.ph_model); "
            f"this one has {0 if net.charge is None else net.charge.n_pools} pools and "
            f"{'a' if net.ph_model is not None else 'no'} PhModel")
    return PhSolver(net.charge, net.ph_model)


# --------------------------------------------------------------------------- the ONE charge-closure check
def closure(net, theta=None, *, ph: Optional[float] = None, t_c: float = 25.0, strict: bool = False) -> Dict[str, Any]:
    """The charge-closure check over the merged titratable table (contract sec. 2.3; v1's Gly / AMA / SB / PAMA gap).

    The engine keeps two tables for one job -- ``ph_state.titratable_inventory`` (what the charge balance titrates)
    and ``TITRATABLE_CENTRES`` (what ``validate_charge_closure`` checks each step against) -- so a solute could titrate
    in the balance while no step that moves it is checked. Here ONE table does both: each reaction's movement of
    titratable centres is computed from the very pool weights the solver uses (``ChargeTable.structural_weights`` x
    ``centres``), and checked against the modules' declared ``centre_ledger`` with the engine's three failure modes:
    a step that moves centres and is not declared; a declared movement that differs from the table's; a declaration
    for a step that moves nothing. With ``ph`` it also reports the charge each step moves at that pH (25 C label by
    default, I = 0): the protons it releases into the balance per event. ``strict`` raises ValueError on a problem."""
    table = net.charge
    out: Dict[str, Any] = {"pools": 0 if table is None else table.n_pools, "kinds": (), "moved": {}, "problems": []}
    ledger: Dict[str, Tuple[str, Mapping[str, Any]]] = {}
    owner_of = {r: m.name for m in net.modules for r in m.reactions}
    for m in net.modules:
        for r, entry in (m.centre_ledger or {}).items():
            if r not in m.reactions:
                out["problems"].append(f"{m.name}: centre_ledger declares {r!r}, which is not one of its reactions")
            ledger[r] = (m.name, entry)
    if table is None:
        out["problems"] += [f"{r}: declared in {o}'s centre_ledger, but the network titrates nothing"
                            for r, (o, _) in ledger.items()]
    else:
        species_centres = table.structural_weights().T @ table.centres.T          # (n, K)
        delta = net.stoichiometry.T @ species_centres                             # (R, K)
        out["kinds"] = table.kinds
        for j, r in enumerate(net.reactions):
            moved = {k: float(delta[j, i]) for i, k in enumerate(table.kinds) if abs(delta[j, i]) > CLOSURE_TOL}
            entry = ledger.get(r)
            if not moved:
                if entry is not None:
                    out["problems"].append(f"{r}: {entry[0]} declares a centre movement but the table's "
                                           "stoichiometry moves none (a stale entry)")
                continue
            out["moved"][r] = moved
            if entry is None:
                out["problems"].append(f"{r} ({owner_of[r]}): moves titratable centres {moved} and is not declared "
                                       "in any module's centre_ledger")
                continue
            declared = {k: float(entry[1].get(k, 0)) for k in table.kinds}
            if any(abs(declared[k] - moved.get(k, 0.0)) > CLOSURE_TOL for k in table.kinds):
                out["problems"].append(f"{r}: {entry[0]} declares {declared}, the table gives {moved}")
            if not str(entry[1].get("basis", "")).strip():
                out["problems"].append(f"{r}: its centre_ledger entry has no basis")
        if ph is not None:
            theta = net.theta0() if theta is None else theta
            solver = PhSolver(table, net.ph_model) if net.ph_model is not None else None
            if solver is None:
                raise UnsupportedConfiguration("the charge moved at a pH needs the network's PhModel")
            tt = solver.temperature_terms(jnp.asarray(t_c + CELSIUS), solver.pk25(theta))
            pk, _ = solver._group_pk(0.0, tt)
            charge, _ = solver.group_terms(jnp.asarray(float(ph)), pk)
            q = np.asarray(solver.pool_weights(theta)).T @ np.asarray(charge)    # (n,): charge per mmol species
            moved_q = net.stoichiometry.T @ q
            out["charge_moved"] = {"ph": float(ph), "t_c": float(t_c),
                                   "per_reaction": {r: float(v) for r, v in zip(net.reactions, moved_q) if abs(v) > 1e-12}}
    out["closes"] = not out["problems"]
    if strict and out["problems"]:
        raise ValueError("charge closure fails:\n" + "\n".join(out["problems"]))
    return out


# --------------------------------------------------------------------------- the pot
def node_times(model: PhModel, duration: float, profile: Optional[Tuple[np.ndarray, np.ndarray]] = None):
    """sulfur.ph_node_times (and, along a profile, profile_node_times at the model's step)."""
    n = max(int(model.ph_nodes), 5)
    linear = np.linspace(0.0, duration, n)
    nodes = linear if model.node_grid == "uniform" else np.unique(np.concatenate(
        [linear, np.geomspace(duration * 1e-5, duration, n)]))
    if profile is not None:
        times, temps = profile
        out = [np.array([0.0, duration])]
        for t0, t1, c0, c1 in zip(times[:-1], times[1:], temps[:-1], temps[1:]):
            k = max(1, int(np.ceil(abs(c1 - c0) / model.profile_node_step_k - 1e-12)))
            out.append(np.linspace(t0, t1, k + 1))
        nodes = np.unique(np.concatenate([nodes, np.unique(np.clip(np.concatenate(out), 0.0, duration))]))
    return nodes


def jump_grid(node_s: np.ndarray) -> np.ndarray:
    """The nodes as diffrax ``jump_ts``, with near-coincident nodes merged onto their predecessor. ``np.unique`` keeps
    two nodes 1.8e-15 min apart where a linspace and a geomspace node meet in floating point (9.0 and
    9.000000000000002 on a 90 min run); a jump that narrow stalls the step clipping. The rate vector between them
    changes by nothing, so the merged grid is the same kink set; the nodes themselves (saved states, interpolation)
    are untouched."""
    out = np.asarray(node_s, dtype=float).copy()
    for i in range(1, len(out)):
        if out[i] - out[i - 1] < JUMP_MERGE:
            out[i] = out[i - 1]
    return out


def pot_from_program(net, y0: np.ndarray, label_ph: float, buffer: np.ndarray, *,
                     segments: Optional[Sequence[Tuple[float, float]]] = None,
                     profile: Optional[Sequence[Tuple[float, float]]] = None, mode: str = "dynamic",
                     duration: Optional[float] = None) -> Dict[str, Any]:
    """An inline-pH pot: the state, the declared label pH, the buffer vector (``ChargeTable.buffer_vector``) and a
    staircase (``segments``: ((minutes, T_C), ...)) or one recorded profile (``profile``: ((minute, T_C), ...) from
    minute 0, one integration over ``duration``, default the last knot). Arrays: durations (S,), node_s / node_tc
    (S, N) (the pH nodes as fractions of each segment and their
    temperatures), knot_s / knot_tc (S, K) (the temperature course inside a segment, for the inline coupling),
    label_nudge (S, N; zero), jump_s (S, N); plus ``mode`` (static) and ``node_minutes``."""
    model = net.ph_model
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}")
    if (segments is None) == (profile is None):
        raise ValueError("give a staircase (segments) or a profile, not both")
    if profile is not None:
        times = np.array([t for t, _ in profile], dtype=float)
        temps = np.array([c for _, c in profile], dtype=float)
        duration = float(times[-1]) if duration is None else float(duration)
        if mode == "clamped":
            raise UnsupportedConfiguration("a clamped pH along a temperature profile (the engine's node-interpolated "
                                           "constants at T(t)) is not ported")
        nodes = node_times(model, duration, (times, temps + CELSIUS))
        rows = [(duration, nodes, np.interp(nodes, times, temps), times, temps)]
    else:
        rows = []
        for d, t_c in segments:
            nodes = node_times(model, float(d))
            rows.append((float(d), nodes, np.full(len(nodes), float(t_c)), np.array([0.0, float(d)]),
                         np.array([float(t_c), float(t_c)])))
    widths = {len(s[1]) for s in rows}
    if len(widths) != 1:
        raise UnsupportedConfiguration("segments with different node counts (a profile inside a staircase)")
    k_width = max(len(s[3]) for s in rows)

    def pad(a, n):
        return np.concatenate([a, np.full(n - len(a), a[-1])])

    return {"y0": np.asarray(y0, dtype=float), "label_ph": float(label_ph), "buffer": np.asarray(buffer, dtype=float),
            "mode": mode, "durations": np.array([s[0] for s in rows]),
            "node_s": np.stack([s[1] / s[0] for s in rows]), "node_tc": np.stack([s[2] for s in rows]),
            "knot_s": np.stack([pad(s[3] / s[0], k_width) for s in rows]),
            "knot_tc": np.stack([pad(s[4], k_width) for s in rows]),
            "label_nudge": np.zeros((len(rows), len(rows[0][1]))),
            "jump_s": np.stack([jump_grid(s[1] / s[0]) for s in rows]), "node_minutes": [s[1] for s in rows]}


def pot_arrays(pot: Mapping[str, Any]) -> Dict[str, Any]:
    """The numeric part of an inline-pH pot (the compiled functions' argument). ``label_nudge`` (S, N; zero) is added
    to the nodes' cooled pH where the measured tables read it -- only a parity check sets it, to put a node whose
    cooled pH sits ON a table's jump (a pot declared at exactly label pH 5.0, Zheng & Ho's boundary) on the engine's
    side of the jump, which there is decided by the root-finders' last digits."""
    missing = [k for k in POT_KEYS if k not in pot]
    if missing:
        raise ValueError(f"not an inline-pH pot (ph_coupling.pot_from_program): it lacks {missing}")
    return {k: jnp.asarray(pot[k], dtype=jnp.float64) for k in POT_KEYS}


def pot_mode(pot: Mapping[str, Any]) -> str:
    mode = pot.get("mode", "dynamic")
    if not isinstance(mode, str) or mode not in MODES:
        raise ValueError(f"an inline-pH pot's mode is static, one of {MODES}; got {mode!r}")
    return mode


# --------------------------------------------------------------------------- integration
class PhRun(NamedTuple):
    ys: Any          # (S, N, n): the state at every pH node of every segment (clipped at 0, as the engine reports)
    ph: Any          # (S, N): the in-situ pH at the nodes (fixed point: the relaxed node values the run ended on)
    ionic: Any       # (S, N): the ionic strength at the nodes, mol/L
    ph_cool: Any     # (S, N): the cooled (25 C label) pH of the node states
    sweeps: Any      # (S,): fixed-point sweeps (0 for inline / clamped)
    residual: Any    # (S,): the last sweep's max |fresh - node pH|
    num_steps: Any   # (S,): integrator steps in the segment (all sweeps)
    failed: Any      # (S,)
    sid: Any         # the strong-ion difference, mmol/L


def _interp_rows(s, grid, rows):
    """rows (N, R) linear in s between grid (N,) nodes, as sulfur.rate_vector_at."""
    n = grid.shape[0]
    lo = jnp.clip(jnp.searchsorted(grid, s, side="right") - 1, 0, n - 2)
    span = grid[lo + 1] - grid[lo]
    w = jnp.where(span > 0.0, jnp.clip((s - grid[lo]) / jnp.where(span > 0.0, span, 1.0), 0.0, 1.0), 0.0)
    return rows[lo] * (1.0 - w) + rows[lo + 1] * w


def _solve(net, rhs_k, data, y0, duration, save_s, jumps, *, rtol, atol, solver, adjoint, max_steps, throw):
    """One solve on s in [0, 1] (the right-hand side scaled by ``duration``), the rate constants ``rhs_k(s, y, data)``.
    Every traced number reaches the vector field through diffrax's ``args`` (``data``, ``duration``), never by
    closure, so the discrete adjoint can differentiate it."""
    import diffrax

    from .network import mass_action

    stoich = jnp.asarray(net.stoichiometry)

    def f(s, y, args):
        data_, duration_ = args
        return duration_ * (stoich @ (rhs_k(s, y, data_) * mass_action(y, net)))

    controller = diffrax.ClipStepSizeController(diffrax.PIDController(rtol=rtol, atol=atol), jump_ts=jumps)
    sol = diffrax.diffeqsolve(diffrax.ODETerm(f), solver, t0=0.0, t1=1.0, dt0=None, y0=y0, args=(data, duration),
                              saveat=diffrax.SaveAt(ts=save_s), stepsize_controller=controller,
                              max_steps=max_steps, adjoint=adjoint, throw=throw)
    return sol.ys, sol.stats["num_steps"], sol.result != diffrax.RESULTS.successful


def loop_kind(adjoint, reverse: Optional[bool]) -> str:
    """How the fixed point's sweeps run: "scan" (``reverse=True``, PR #83's form), "while" (forward mode) or
    "checkpointed" (any reverse-capable adjoint: ``jax.grad`` works with no flag)."""
    import diffrax

    if reverse:
        return "scan"
    return "while" if isinstance(adjoint, diffrax.ForwardMode) else "checkpointed"


def integrate(net, theta, pot: Mapping[str, Any], *, rtol: float = 1e-8, atol: float = 1e-14, solver=None,
              adjoint=None, reverse: Optional[bool] = None, max_steps: int = 200_000, throw: bool = True,
              force_sweeps: Optional[int] = None, mode: Optional[str] = None) -> PhRun:
    """One pot through its programme with the network's pH coupling (``mode`` "clamped" holds the declared pH; default
    ``pot["mode"]``). ``adjoint`` defaults to diffrax's discrete adjoint (``diffrax.ForwardMode()`` for
    ``jax.jacfwd``); the sweep loop follows it (``loop_kind``). ``force_sweeps`` replaces the stopping rule by a fixed
    sweep count (every segment): for a parity check on a pot whose last residual sits at the tolerance, where solver
    noise alone decides one sweep more or less."""
    import diffrax

    from .network import RateCtx

    bal, st = solver_of(net), net.ph_model
    solver = solver or getattr(diffrax, net.solver)()
    adjoint = adjoint or diffrax.RecursiveCheckpointAdjoint()
    kind = loop_kind(adjoint, reverse)
    a = pot_arrays(pot)
    mode = pot_mode(pot) if mode is None else mode
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}")
    kw = dict(rtol=rtol, atol=atol, solver=solver, adjoint=adjoint, max_steps=max_steps, throw=throw)
    y0 = a["y0"]
    sid = bal.strong_ion_difference(a["label_ph"], y0, theta, a["buffer"])

    def rates(t_k, ph, ionic, cool):
        return net.rate_constants(theta, RateCtx(t_c=t_k - CELSIUS, t_k=t_k, ph=ph, ph_label=cool, ionic_m=ionic))

    def node_states(ys, node_tk):
        return jax.vmap(lambda y, t: bal.state(y, t, sid, theta, a["buffer"]))(ys, node_tk)

    min_sweeps = max(int(st.ph_iterations), 1)
    relaxed = st.fixed_point == "relaxed_converged"
    max_sweeps = max(min_sweeps, int(st.max_sweeps)) if relaxed else min_sweeps
    if force_sweeps is not None:
        if not 1 <= int(force_sweeps) <= max_sweeps:
            raise ValueError(f"force_sweeps must lie in [1, {max_sweeps}]")
        min_sweeps = max_sweeps = int(force_sweeps)

    def segment(y_start, seg):
        duration, node_s, node_tc, knot_s, knot_tc, nudge, jump_s = seg
        node_tk = node_tc + CELSIUS
        if mode == "clamped":
            t_k = knot_tc[0] + CELSIUS
            k0 = rates(t_k, a["label_ph"], 0.0, jnp.nan)
            ys, steps, bad = _solve(net, lambda s, y, k: k, k0, y_start, duration, node_s, jump_s, **kw)
            ys = jnp.maximum(ys, 0.0)
            n = node_s.shape[0]
            return ys[-1], (ys, jnp.full(n, a["label_ph"]), jnp.zeros(n), jnp.full(n, a["label_ph"]),
                            jnp.asarray(0), jnp.asarray(0.0), steps, bad)
        if st.coupling == "inline":
            def rhs_k(s, y, d):
                theta_, sid_, buffer_, knot_s_, knot_tc_ = d
                t_k = jnp.interp(s, knot_s_, knot_tc_) + CELSIUS
                ph, ionic, cool = bal.state(y, t_k, sid_, theta_, buffer_)
                return net.rate_constants(theta_, RateCtx(t_c=t_k - CELSIUS, t_k=t_k, ph=ph, ph_label=cool,
                                                          ionic_m=ionic))

            ys, steps, bad = _solve(net, rhs_k, (theta, sid, a["buffer"], knot_s, knot_tc), y_start, duration,
                                    node_s, knot_s, **kw)
            ys = jnp.maximum(ys, 0.0)
            ph, ionic, cool = node_states(ys, node_tk)
            return ys[-1], (ys, ph, ionic, cool, jnp.asarray(0), jnp.asarray(0.0), steps, bad)
        # fixed_point_nodes (shipped): B64's per-node starts (identical at every node of an isothermal segment)
        p0, i0, c0 = node_states(jnp.broadcast_to(y_start, (node_s.shape[0],) + y_start.shape), node_tk)

        # every traced number enters the sweep as an OPERAND, never by closure: lax.cond's reverse-mode lowering
        # cannot carry a closed-over tracer (the branch would hold it as a constant)
        data = (theta, sid, a["buffer"], y_start, duration, node_s, node_tk, nudge, jump_s)

        def sweep(carry, d):
            theta_, sid_, buffer_, y_start_, duration_, node_s_, node_tk_, nudge_, jump_s_ = d
            count, node_ph, node_i, node_cool, _ys, _res, _done, steps, bad = carry
            k_nodes = jax.vmap(lambda t, p, i, c: net.rate_constants(
                theta_, RateCtx(t_c=t - CELSIUS, t_k=t, ph=p, ph_label=c, ionic_m=i)))(node_tk_, node_ph, node_i,
                                                                                         node_cool + nudge_)
            ys, n_steps, failed = _solve(net, lambda s, y, dd: _interp_rows(s, dd[0], dd[1]), (node_s_, k_nodes),
                                         y_start_, duration_, node_s_, jump_s_, **kw)
            ys = jnp.maximum(ys, 0.0)
            fresh, fi, fc = jax.vmap(lambda y, t: bal.state(y, t, sid_, theta_, buffer_))(ys, node_tk_)
            res = jnp.max(jnp.abs(fresh - node_ph))
            count = count + 1
            if relaxed:
                new_ph = node_ph + st.relaxation * (fresh - node_ph)
                done = ((res < st.tolerance) & (count >= min_sweeps)) | (count >= max_sweeps)
            else:
                new_ph = fresh
                done = (res < FIXED_SWEEPS_STOP) | (count >= max_sweeps)
            if force_sweeps is not None:
                done = count >= max_sweeps
            return (count, new_ph, fi, fc, ys, res, done, steps + n_steps, bad | failed)

        n = node_s.shape[0]
        init = (jnp.asarray(0), p0, i0, c0, jnp.zeros((n,) + y_start.shape), jnp.asarray(jnp.inf),
                jnp.asarray(False), jnp.asarray(0), jnp.asarray(False))
        if kind == "scan":
            def body(carry, _):
                return jax.lax.cond(carry[6], lambda c, d: c, sweep, carry, data), None

            out, _ = jax.lax.scan(body, init, None, length=max_sweeps)
        elif kind == "while":
            out = jax.lax.while_loop(lambda c: ~c[6], lambda c: sweep(c, data), init)
        else:
            import equinox.internal as eqxi

            # the operands travel in the carry (the checkpointed loop's custom VJP differentiates its carry)
            out, _ = eqxi.while_loop(lambda cd: ~cd[0][6], lambda cd: (sweep(cd[0], cd[1]), cd[1]), (init, data),
                                     max_steps=max_sweeps, kind="checkpointed", checkpoints=max_sweeps)
        count, node_ph, node_i, node_cool, ys, res, _, steps, bad = out
        return ys[-1], (ys, node_ph, node_i, node_cool, count, res, steps, bad)

    segs = (a["durations"], a["node_s"], a["node_tc"], a["knot_s"], a["knot_tc"], a["label_nudge"], a["jump_s"])
    _, (ys, ph, ionic, cool, count, res, steps, bad) = jax.lax.scan(segment, y0, segs)
    return PhRun(ys=ys, ph=ph, ionic=ionic, ph_cool=cool, sweeps=count, residual=res, num_steps=steps,
                 failed=bad, sid=sid)
