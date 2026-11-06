"""
src/kinetic_core_jax/network.py -- the A2 module contract v2 in JAX (docs/guides/A2_MODULE_CONTRACT.md sec. 2;
2026-11-05): ``Module`` (one chemistry module as static data plus one pure rate function), ``compose`` (modules ->
one ``Network``, refusing duplicate keys and unbalanced reactions), ``integrate`` (one pot through a piecewise-constant
thermal program, diffrax), the batch API (``batched`` / ``grid``: many parameter draws and pots in lockstep under
``jax.vmap``) and the exact Jacobian (``jacobian``, forward mode).

NOT part of the engine: nothing in ``src/kinetic_core`` imports it, it needs the separate ``maillard-jax`` env
(environment-jax.yml), and no shipped number comes from it. Modules are BUILT from the engine's own definitions
(``src.kinetic_core``, read-only); see ``lipid.py`` for the first one.

The three splits the contract makes:

* STATIC vs DYNAMIC. Everything structural -- species, reactions, stoichiometry, mass-action layout, which switches are
  on -- is numpy, frozen when the module is built from a ``Config`` (``src.kinetic_core.config``). The free numbers
  are a pytree ``theta = {module: {coordinate: value}}`` that JAX differentiates.
* NETWORK vs POT. A pot's numbers (initial state, thermal program; later SID and buffer) are ARGUMENTS of the compiled
  functions, never constants traced into them, so one compilation serves every pot with the same number of segments.
* TIME. Each segment is integrated on a unit interval s in [0, 1] with the right-hand side scaled by its duration, so
  pots whose segments differ in length share one compiled solve and run in lockstep under ``vmap``; a zero-duration
  segment (the padding that equalises segment counts) is the identity.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Mapping, NamedTuple, Optional, Sequence, Tuple

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)

from .lane import UnsupportedConfiguration  # noqa: E402  (one refusal type for every port)

CELSIUS = 273.15
#: the elements every reaction must balance (the engine's carbon / nitrogen / sulfur invariants)
ELEMENTS = ("C", "N", "S")
BALANCE_TOL = 1e-9


class RateCtx(NamedTuple):
    """What a module's rate function may read besides its own coordinates: the segment's temperature (both scales,
    T_K = T_C + 273.15 exactly as the engine computes it) and, for pH-coupled networks, the in-situ pH, the cooled
    label pH and the ionic strength (None in a pH-free network).

    In pH mode "prescribed" (contract v2 sec. 2.4) ``ph_label`` is the pot's DECLARED pH, read from the pot data
    (``pot["ph_label"]``) and constant over the programme -- the engine's declared-pH factors (trunk B12 / B40,
    acrylamide B15) -- and ``ph`` stays None. ``aw`` is the pot's declared water activity (``pot["aw"]``). Both are
    NaN where the engine's value is None; a module reads them only if it lists them in ``Module.pot_inputs``."""

    t_c: Any
    t_k: Any
    ph: Any = None
    ph_label: Any = None
    ionic_m: Any = None
    aw: Any = None


@dataclass(frozen=True)
class Module:
    """One chemistry module: static data plus one pure rate function (contract v2 sec. 2).

    ``stoichiometry`` and ``reactants`` are keyed by species NAME, so a module may touch the species of a module it
    ``requires``; ``compose`` maps names to rows. ``rate_constants(theta_m, ctx)`` returns this module's rate
    constants in ``reactions`` order, from its own coordinates ``theta_m`` (a dict with the keys of ``params``)."""

    name: str
    species: Tuple[str, ...]
    elements: Mapping[str, Tuple[int, int, int]]                 # owned species -> (C, N, S)
    reactions: Tuple[str, ...]
    stoichiometry: Mapping[str, Mapping[str, float]]            # reaction -> {species: net coefficient}
    reactants: Mapping[str, Mapping[str, int]]                  # reaction -> {species: mass-action order}
    params: Mapping[str, float]                                 # free coordinates at their SHIPPED values
    bands: Mapping[str, Optional[Tuple[float, float]]]
    basis: Mapping[str, str]                                    # coordinate -> where its value and band come from
    rate_constants: Callable[[Mapping[str, Any], RateCtx], Any]
    observables: Mapping[str, Callable[[Any, Mapping[str, int]], Any]] = field(default_factory=dict)
    requires: Tuple[str, ...] = ()
    needs_ph: bool = False
    config_digest: str = ""
    notes: Tuple[str, ...] = ()
    #: pot scalars the rate function reads through ``RateCtx`` (a subset of ``POT_INPUTS``); ``integrate`` refuses a
    #: pot that lacks one. "ph_label" makes the network's pH mode "prescribed" (a declared, constant pH).
    pot_inputs: Tuple[str, ...] = ()
    #: (v2.2) values of a RateCtx input ("t_c", "ph_label", "aw") where the rate function jumps or kinks: a forced
    #: linear course gets a knot where it crosses one (``forcing.forced_pot``), so the jump sits on a segment boundary
    breakpoints: Mapping[str, Tuple[float, ...]] = field(default_factory=dict)
    #: (v2.2) species owned by ANOTHER module that this module's reactions consume, declared (contract sec. 4,
    #: "declared competition"); ``compose(competition="refuse")`` refuses an undeclared one
    competes: Tuple[str, ...] = ()
    #: (v2.2) this module's titratable pools (a typed stub: merged and checked by ``compose``; the ONE charge-closure
    #: check needs the sulfur port's inventory, ``charge_closure``)
    titratables: Tuple["Titratable", ...] = ()


class Titratable(NamedTuple):
    """One titratable pool for the merged table (contract sec. 4 / 6.3): its pKa values at 25 C, the charge of each
    protonation state (most protonated first, one more than the pKa values), and where they come from."""

    species: str
    pka_25c: Tuple[float, ...]
    charges: Tuple[int, ...]
    basis: str = ""


#: the pot scalars a module may read (``pot[name]``, NaN for the engine's None)
POT_INPUTS = ("ph_label", "aw")
#: ``compose``'s competition modes: "report" lists undeclared shared consumption on the Network, "refuse" raises
COMPETITION_MODES = ("report", "refuse")


@dataclass(frozen=True)
class Network:
    """Modules composed into one recipe: frozen arrays plus the modules' rate functions in order."""

    modules: Tuple[Module, ...]
    species: Tuple[str, ...]
    reactions: Tuple[str, ...]
    index: Mapping[str, int]
    stoichiometry: np.ndarray        # (n_species, n_reactions)
    reactant_index: np.ndarray       # (n_reactions, width), padded with n_species (a 1.0 appended to the state)
    reactant_power: np.ndarray       # same shape, 0 for padding
    element_matrix: np.ndarray       # (3, n_species)
    coordinates: Tuple[Tuple[str, str], ...]   # (module, coordinate), the order of theta's flat vector
    needs_ph: bool
    #: the declared integrator (contract v2 sec. 2.6): chosen per composed network from its own benchmark, because
    #: neither dominates -- Kvaerno4 on the sulfur pot (no rejections, ~25 % faster), Kvaerno5 on the lipid pots
    #: (5x fewer steps once a pool is depleted far below atol: 158 vs 794 on the 160 C / 30 min pot)
    solver: str = "Kvaerno5"
    #: (v2.2) the merged titratable table (one entry per species), in module order
    titratables: Tuple[Titratable, ...] = ()
    #: (v2.2) species consumed by a module that does not own it -> (owner, consumers...), and the undeclared pairs
    competition: Mapping[str, Tuple[str, ...]] = field(default_factory=dict)
    undeclared_competition: Tuple[Tuple[str, str], ...] = ()

    @property
    def n_species(self) -> int:
        return len(self.species)

    @property
    def pot_inputs(self) -> Tuple[str, ...]:
        return tuple(sorted({x for m in self.modules for x in m.pot_inputs}))

    @property
    def ph_mode(self) -> str:
        """Contract v2 sec. 2.4: "inline" (a charge balance in the RHS), "prescribed" (a declared pH read from the
        pot; constant over the programme in every module ported so far) or "none"."""
        if self.needs_ph:
            return "inline"
        return "prescribed" if "ph_label" in self.pot_inputs else "none"

    def theta0(self) -> Dict[str, Dict[str, Any]]:
        """The shipped coordinates as the pytree the solvers take."""
        return {m.name: {k: jnp.asarray(float(v)) for k, v in m.params.items()} for m in self.modules}

    def rate_constants(self, theta: Mapping[str, Mapping[str, Any]], ctx: RateCtx):
        return jnp.concatenate([jnp.atleast_1d(m.rate_constants(theta[m.name], ctx)) for m in self.modules])

    def observe(self, y, names: Sequence[str]):
        """The named observables of one state (each module's ``observables``)."""
        table = {k: f for m in self.modules for k, f in m.observables.items()}
        return jnp.stack([table[n](y, self.index) for n in names])


def compose(modules: Sequence[Module], *, layout: Optional[Sequence[str]] = None,
            solver: str = "Kvaerno5", competition: str = "report") -> Network:
    """One Network from modules. Refuses (ValueError): a duplicate module name, a missing ``requires``, a species or
    reaction key defined twice, a reaction touching a species no module owns, and any reaction that does not balance
    C / N / S to 1e-9. ``layout`` pins a species order (a permutation of the owned species); ``solver`` names the
    diffrax integrator ``integrate`` uses unless given one.

    (v2.2) Declared competition: a species CONSUMED (a mass-action reactant with a negative net coefficient) by a
    module that does not own it must be in that module's ``competes``. ``competition="report"`` (the default: the
    ported trunk-lane / acrylamide modules predate the declaration) lists the undeclared pairs on
    ``Network.undeclared_competition``; ``"refuse"`` raises. The merged titratable table refuses a pool declared twice
    or by a module that neither owns nor requires its species."""
    if solver not in ("Kvaerno3", "Kvaerno4", "Kvaerno5"):
        raise ValueError(f"solver {solver!r} is not a declared choice")
    if competition not in COMPETITION_MODES:
        raise ValueError(f"competition must be one of {COMPETITION_MODES}")
    names = [m.name for m in modules]
    if len(set(names)) != len(names):
        raise ValueError(f"duplicate module names: {names}")
    for m in modules:
        missing = [r for r in m.requires if r not in names]
        if missing:
            raise ValueError(f"module {m.name} requires {missing}")
        unknown = [x for x in m.pot_inputs if x not in POT_INPUTS]
        if unknown:
            raise ValueError(f"module {m.name} reads unknown pot inputs {unknown} (declared: {POT_INPUTS})")
    owner: Dict[str, str] = {}
    for m in modules:
        for s in m.species:
            if s in owner:
                raise ValueError(f"species {s!r} is defined by both {owner[s]} and {m.name}")
            owner[s] = m.name
    species = tuple(s for m in modules for s in m.species)
    if layout is not None:
        if sorted(layout) != sorted(species):
            raise ValueError("layout must be a permutation of the composed species")
        species = tuple(layout)
    index = {s: i for i, s in enumerate(species)}
    reactions: List[str] = []
    for m in modules:
        for r in m.reactions:
            if r in reactions:
                raise ValueError(f"reaction {r!r} is defined twice")
            reactions.append(r)
    n_s, n_r = len(species), len(reactions)
    stoich = np.zeros((n_s, n_r))
    layout_rows: List[List[Tuple[int, int]]] = []
    j = 0
    for m in modules:
        allowed = set(m.species) | {s for r in m.requires for mm in modules if mm.name == r for s in mm.species}
        for r in m.reactions:
            for s, c in m.stoichiometry[r].items():
                if s not in allowed:
                    raise ValueError(f"{m.name}.{r} touches {s!r}, which neither it nor a required module owns")
                stoich[index[s], j] += float(c)
            layout_rows.append([(index[s], int(p)) for s, p in m.reactants[r].items()])
            j += 1
    width = max(1, max(len(x) for x in layout_rows))
    r_idx = np.full((n_r, width), n_s, dtype=np.int64)
    r_pow = np.zeros((n_r, width), dtype=np.int64)
    for j, row in enumerate(layout_rows):
        for q, (i, p) in enumerate(row):
            if p not in (1, 2):
                raise ValueError(f"{reactions[j]}: mass-action order {p} is not supported")
            r_idx[j, q], r_pow[j, q] = i, p
    elements = np.zeros((len(ELEMENTS), n_s))
    for m in modules:
        for s in m.species:
            elements[:, index[s]] = m.elements[s]
    imbalance = elements @ stoich
    for j in range(n_r):
        scale = max(1.0, float(np.max(np.abs(elements * np.abs(stoich[:, j])))))
        bad = np.abs(imbalance[:, j]) > BALANCE_TOL * scale
        if bad.any():
            raise ValueError(f"reaction {reactions[j]} does not balance: "
                             + ", ".join(f"{e} {imbalance[i, j]:+.3g}" for i, e in enumerate(ELEMENTS) if bad[i]))
    coordinates = tuple((m.name, k) for m in modules for k in m.params)
    shared, undeclared = _competition(modules, owner)
    if undeclared and competition == "refuse":
        raise ValueError("undeclared shared consumption (declare it in Module.competes): "
                         + ", ".join(f"{m} consumes {s} (owned by {owner[s]})" for s, m in undeclared))
    return Network(modules=tuple(modules), species=species, reactions=tuple(reactions), index=index,
                   stoichiometry=stoich, reactant_index=r_idx, reactant_power=r_pow, element_matrix=elements,
                   coordinates=coordinates, needs_ph=any(m.needs_ph for m in modules), solver=solver,
                   titratables=_merge_titratables(modules, owner), competition=shared,
                   undeclared_competition=undeclared)


def _competition(modules: Sequence[Module], owner: Mapping[str, str]):
    consumers: Dict[str, List[str]] = {}
    for m in modules:
        for r in m.reactions:
            for s in m.reactants[r]:
                if float(m.stoichiometry[r].get(s, 0.0)) < 0.0 and owner.get(s) != m.name:
                    if m.name not in consumers.setdefault(s, []):
                        consumers[s].append(m.name)
    declared = {m.name: set(m.competes) for m in modules}
    for m in modules:
        stray = sorted(set(m.competes) - set(owner))
        if stray:
            raise ValueError(f"module {m.name} declares competition for {stray}, which no module owns")
    shared = {s: (owner[s],) + tuple(c) for s, c in consumers.items()}
    undeclared = tuple((s, c) for s, cs in consumers.items() for c in cs if s not in declared[c])
    return shared, undeclared


def _merge_titratables(modules: Sequence[Module], owner: Mapping[str, str]) -> Tuple[Titratable, ...]:
    table: Dict[str, Titratable] = {}
    names = {m.name: m for m in modules}
    for m in modules:
        reach = set(m.species) | {s for r in m.requires for s in names[r].species}
        for t in m.titratables:
            if t.species in table:
                raise ValueError(f"titratable {t.species!r} is declared twice")
            if t.species not in reach:
                raise ValueError(f"{m.name} declares titratable {t.species!r}, which it neither owns nor requires")
            if len(t.charges) != len(t.pka_25c) + 1:
                raise ValueError(f"titratable {t.species!r}: {len(t.pka_25c)} pKa values need "
                                 f"{len(t.pka_25c) + 1} charge states")
            table[t.species] = t
    return tuple(table.values())


def charge_closure(net: Network):
    """The ONE charge-closure check over the merged titratable table (contract sec. 4 / 6.3; it closes v1's Gly /
    AMA / SB / PAMA gap). A TYPED STUB: it needs, per reaction, the charge each titratable reactant and product
    carries at the network's pH, i.e. the inline charge balance and the titratable inventory that only the sulfur
    lane carries today (``ph_state.titratable_inventory`` / ``TITRATABLE_CENTRES``, the bake-off's
    ``model.make_ph_solver``). It lands with the sulfur port."""
    raise UnsupportedConfiguration(
        "charge closure over the merged titratable table needs the sulfur port's inline charge balance and "
        f"titratable inventory ({len(net.titratables)} pools declared); not ported")


# --------------------------------------------------------------------------- the parameter vector
def theta_vector(net: Network, theta: Mapping[str, Mapping[str, Any]]):
    """theta's flat vector in ``net.coordinates`` order (for Jacobians and least squares)."""
    return jnp.stack([jnp.asarray(theta[m][k], dtype=jnp.float64) for m, k in net.coordinates])


def theta_from_vector(net: Network, vector) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {m.name: {} for m in net.modules}
    for i, (m, k) in enumerate(net.coordinates):
        out[m][k] = vector[i]
    return out


def stack_theta(thetas: Sequence[Mapping[str, Mapping[str, Any]]]):
    """A batch of parameter draws as ONE pytree with a leading axis (what ``batched`` maps over)."""
    return jax.tree_util.tree_map(lambda *xs: jnp.stack([jnp.asarray(x, dtype=jnp.float64) for x in xs]), *thetas)


# --------------------------------------------------------------------------- the thermal program and the pot
def program_segments(segments: Sequence[Tuple[float, float]], n_segments: Optional[int] = None) -> np.ndarray:
    """``((duration_min, T_C), ...)`` -- the engine's ThermalProgram.segments -- as an (S, 2) array, padded to
    ``n_segments`` with zero-duration segments at the last temperature (the identity), so pots with fewer segments
    share a compile with longer ones."""
    rows = [(float(d), float(t)) for d, t in segments]
    if not rows:
        raise ValueError("a thermal program needs at least one segment")
    if n_segments is not None:
        if n_segments < len(rows):
            raise ValueError(f"{len(rows)} segments do not fit in {n_segments}")
        rows += [(0.0, rows[-1][1])] * (n_segments - len(rows))
    return np.asarray(rows, dtype=float)


def refine_segments(segments: Sequence[Tuple[float, float]], save_times: Sequence[float]):
    """Split segments so every save time is a segment boundary. Returns (segments, rows): ``rows[i]`` is the row of
    ``integrate``'s output that holds the state at ``save_times[i]``. A split changes nothing physical: each piece
    keeps its segment's temperature."""
    edges = np.cumsum([0.0] + [float(d) for d, _ in segments])
    cuts = sorted({float(t) for t in save_times})
    if cuts and (cuts[0] <= 0.0 or cuts[-1] > edges[-1] * (1 + 1e-12)):
        raise ValueError("save times must lie in (0, total duration]")
    out: List[Tuple[float, float]] = []
    ends: List[float] = []
    for (d, temp), t0, t1 in zip(segments, edges[:-1], edges[1:]):
        inner = [c for c in cuts if t0 < c < t1 and not math.isclose(c, t1, rel_tol=1e-12)]
        points = [t0] + inner + [t1]
        for a, b in zip(points[:-1], points[1:]):
            out.append((b - a, float(temp)))
            ends.append(b)
    rows = [int(np.argmin(np.abs(np.asarray(ends) - c))) for c in save_times]
    return tuple(out), rows


def mass_action(y, net: Network):
    yc = jnp.append(jnp.maximum(y, 0.0), 1.0)[jnp.asarray(net.reactant_index)]
    p = jnp.asarray(net.reactant_power)
    return jnp.prod(jnp.where(p == 2, yc * yc, jnp.where(p == 1, yc, 1.0)), axis=1)


class Run(NamedTuple):
    ys: Any          # (S, n_species): the state at the END of each segment
    num_steps: Any   # (S,): accepted + rejected steps per segment
    failed: Any      # (S,): True where diffrax did not succeed (only reachable with ``throw=False``)


def _solver_of(net: Network):
    import diffrax

    return getattr(diffrax, net.solver)()


def integrate(net: Network, theta, pot: Mapping[str, Any], *, solver=None, rtol: float = 1e-8, atol: float = 1e-14,
              adjoint=None, max_steps: int = 100_000, throw: bool = True) -> Run:
    """One pot (``{"y0": (n,), "segments": (S, 2)}``, plus the scalars in ``net.pot_inputs``) through its thermal
    program.

    Rate constants are evaluated once per segment (temperature is constant on it, and so are the declared pot
    scalars; a network without inline pH has no state dependence in its constants). ``throw=False`` flags a failed
    pot (``Run.failed``) instead of raising, which is what a fit over many pots needs (XLA cannot be interrupted by
    a signal)."""
    import diffrax

    if net.needs_ph:
        raise UnsupportedConfiguration(
            "the composed network's pH coupling (charge balance inside the RHS, contract v2 sec. 2.4) is not ported "
            "yet; the sulfur lane's module carries it in its own driver (sulfur_lane.integrate)")
    missing = [x for x in net.pot_inputs if x not in pot]
    if missing:
        raise ValueError(f"the pot lacks {missing}, which the network's modules read (NaN stands for None)")
    scalars = {x: jnp.asarray(pot[x], dtype=jnp.float64) for x in net.pot_inputs}
    stoich = jnp.asarray(net.stoichiometry)
    term = diffrax.ODETerm(lambda s, y, args: args[1] * (stoich @ (args[0] * mass_action(y, net))))
    solver = solver or _solver_of(net)
    adjoint = adjoint or diffrax.RecursiveCheckpointAdjoint()
    controller = diffrax.PIDController(rtol=rtol, atol=atol)

    def segment(y, seg):
        duration, t_c = seg[0], seg[1]
        k = net.rate_constants(theta, RateCtx(t_c=t_c, t_k=t_c + CELSIUS, **scalars))
        sol = diffrax.diffeqsolve(term, solver, t0=0.0, t1=1.0, dt0=None, y0=y, args=(k, duration),
                                  saveat=diffrax.SaveAt(t1=True), stepsize_controller=controller,
                                  max_steps=max_steps, adjoint=adjoint, throw=throw)
        y1 = sol.ys[-1]
        return y1, (y1, sol.stats["num_steps"], sol.result != diffrax.RESULTS.successful)

    _, (ys, steps, failed) = jax.lax.scan(segment, jnp.asarray(pot["y0"], dtype=jnp.float64),
                                          jnp.asarray(pot["segments"], dtype=jnp.float64))
    return Run(ys=ys, num_steps=steps, failed=failed)


# --------------------------------------------------------------------------- batches and Jacobians
def batched(net: Network, *, theta_axis: Optional[int] = 0, pot_axis: Optional[int] = 0, **solve_kw):
    """``f(theta_batch, pot_batch) -> ys (B, S, n)``, compiled once: ``jax.vmap`` of ``integrate``. ``theta_axis`` /
    ``pot_axis`` = None broadcasts one theta (or one pot) across the batch. Members keep their own adaptive steps;
    the batched loop runs until the slowest member is done (diffrax's vmapped while loop). The state is NOT stacked
    into one ODE: an implicit solver would then factor one dense (B n) x (B n) Jacobian instead of B small ones."""

    def one(theta, pot):
        return integrate(net, theta, pot, **solve_kw).ys

    return jax.jit(jax.vmap(one, in_axes=(theta_axis, pot_axis)))


def grid(net: Network, **solve_kw):
    """``f(theta_draws, pots) -> ys (P, D, S, n)``: every draw on every pot, one compile (the posterior's shape)."""

    def one(theta, pot):
        return integrate(net, theta, pot, **solve_kw).ys

    return jax.jit(jax.vmap(jax.vmap(one, in_axes=(0, None)), in_axes=(None, 0)))


def jacobian(net: Network, observe: Callable[[Any], Any], **solve_kw):
    """``f(theta_vector, pot) -> (observables, d observables / d theta)``, exact, forward mode (one JVP per
    coordinate; the right mode up to ~50 coordinates, contract v2 sec. 4). ``observe(ys)`` maps the (S, n) output to
    the observable vector."""
    import diffrax

    solve_kw = dict(solve_kw, adjoint=diffrax.ForwardMode())

    def obs(vector, pot):
        value = observe(integrate(net, theta_from_vector(net, vector), pot, **solve_kw).ys)
        return value, value

    def both(vector, pot):
        jac, value = jax.jacfwd(obs, has_aux=True)(vector, pot)
        return value, jac

    return jax.jit(both)
