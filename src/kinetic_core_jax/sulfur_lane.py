"""
src/kinetic_core_jax/sulfur_lane.py -- the sulfur lane (``sulfur.FULL_REACTIONS``: the trunk's 26 steps, the fed
3-deoxy triangle and the 103 sulfur steps; 62 species) as an A2 contract module, with its dynamic pH (plan step
A2 + H3, docs/guides/A2_MODULE_CONTRACT.md v2.2 sec. 6.2.2; 2026-11-05). The last lane of A2 + H3.

    lane = sulfur_lane.build(cfg)                    # Module + Network + charge balance + the Config's pH numerics
    pot = sulfur_lane.make_pot(lane, initial, process)   # the engine's sulfur prelude, nodes, SID inputs
    run = sulfur_lane.integrate(lane, theta, pot)        # the engine's pH coupling, in JAX

ONE MODULE, NOT THREE. The engine's sulfur lane is one recipe with its OWN copy of the trunk (its 26 steps run on
``core_parameters("sulfur")``, under ``sulfur_rate_constants_at``'s generic rule, beside the lane's pH machinery), so
the port is one module, ``sulfur``, owning all 62 species. Composing it with ``trunk`` / ``fed_3deoxy`` (one owner per
species: Glc, Gly, MGO, ..., DGAL) is A2b's union network, which changes numbers and is judged by evidence
(contract sec. 6.3), not by identity with this lane.

DERIVED, NEVER RETYPED (``sulfur_structure.derive`` and ``build``; an edit to the engine's data carries over):
species and C / N / S, the 132 reactions in FULL_REACTIONS order, stoichiometry and mass-action orders, each
reaction's parameter key and pH-factor kind (``REACTION_PH_FACTOR``), every parameter's k_ref, Ea (None = held), T_ref
(``KineticParameter``: ``parameters.T_REF_K``; ``SulfurParameter``: its ``t_ref_k``) and temperature cap
(``parameters_furanic.TEMPERATURE_CAP_K``), the gas constant, the speciation pK tables (``PKA_GKL_CYSTEINE``,
``PKA_REFERENCE_T_K``, ``PH_NORMALISATION``), Zheng & Ho's cysteine set and its mechanism boundary, the glutathione
refit, the thiamine share points and floor, Stack's adduct equilibrium and dH, the Amadori split, the titratable
groups (pK, dH, dCp, z_max) of ``ph_state`` and the buffer groups (from ``BufferSpec.groups``), the water density and
pKw tables, the label temperature, the pot prelude's constants (``OX_AMBIENT_MMOL_L``, ``oxygen_reservoir_units``,
the matrix sites), the drift calibration (``engine.core_ph_drift``) and ``integrate_sulfur``'s node and sweep
defaults (read off its signature).

HAND-PORTED (each pinned: ``sulfur_structure.SOURCE_PINS`` hashes the engine function's AST, and
tests/scientific/test_jax_sulfur_parity.py compares the numbers):

* the generic rate rule ``k_at(min(T, cap)) x ph_factor(kind, pH, T, I)`` and the six pH-factor kinds with B51's
  speciation (Suleimenov & Seward's H2S pK1 polynomial, GKL's cysteine pair with dCp, Davies on the numerator only);
* the derived constants (``sulfur_structure.HAND_PORTED``): Zheng & Ho's cysteine and glutathione thermolysis on the
  label-scale pH, the thiamine H2S share (B50), the thioether release k_forward / K_adduct(T), the Amadori pin; and
  the 21 steps the shipped configuration holds at zero (their routes are off; the port refuses them on);
* the charge balance (``ph_state.titratable_inventory`` as a linear map of the state; ``_charge_residual``,
  ``ionic_strength_of``, ``strong_ion_difference``), solved for (pH, I) by bisection then Newton, differentiated by
  the implicit-function theorem (the bake-off's solver, ``model.make_ph_solver``, with the drift's pK a coordinate);
* the pH COUPLING of ``integrate_sulfur``: "fixed_point_nodes" (shipped: rate vectors at ~81 nodes, piecewise linear
  in time between them, the pH recomputed from the states at the nodes and relaxed, sweep after sweep, with the
  engine's stopping rule) and "inline" (the charge balance inside the right-hand side), the node grids
  (``ph_node_times``, ``profile_node_times``), the clamped pH (a ``BufferSpec("clamped")``), the SID carried across
  segments (B51), and a recorded temperature PROFILE (B64: one integration, nodes every ``profile_node_step_k``);
* the engine's sulfur prelude of ``_integrate_program`` (ambient OX, the protein disulfide pool, the headspace
  reservoir, OXV, MET / PRO / PYRL recorded and dropped).

THE COORDINATES (``Module.params``): ``log10_k:<key>`` (log10 of k_ref at its T_ref) or, where the engine carries the
key at exactly zero, a LINEAR ``k_ref:<key>``; ``ea:<key>`` where the key has an Ea (a held, Ea-less key has none);
and the drift's ``acid_yield`` and ``arp_amine_pka`` (``engine.core_ph_drift``), which enter the charge balance.

SWITCHES (``sulfur_structure.CONFIG_HANDLING``): the pH numerics are ported at every value; every rival route, B51
off, the trunk pH terms and the trunk solutes in the balance are HELD at their shipped values and refused otherwise
(they stay NumPy-only while M1b(i) / M1c reshape them). A pot-level ``ph_drift`` and a prescribed pH line (no drift)
are refused at the pot.

DIFFERENTIATION. The fixed point's sweep loop is a ``lax.while_loop`` (forward mode and ``vmap``) or, with
``reverse=True``, a ``lax.scan`` over ``max_sweeps`` whose finished sweeps are skipped by ``lax.cond`` (reverse mode,
``jax.grad``). The sweep count is the engine's discrete stopping rule; a perturbation that changes it is a jump the
derivative does not see (the run reports the count). ``network.integrate`` still refuses a pH network: the merged
titratable table (``Module.titratables``) cannot yet carry this inventory's dH / dCp / acid yield, so the pH solver
lives here until the charge-closure step moves it there.
"""
from __future__ import annotations

import inspect
import math
from dataclasses import dataclass
from typing import Any, Dict, List, Mapping, NamedTuple, Optional, Sequence, Tuple

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)

from . import sulfur_structure as structure  # noqa: E402
from .lane import UnsupportedConfiguration  # noqa: E402
from .network import CELSIUS, Module, Network, RateCtx, compose, mass_action  # noqa: E402
from .trunk import observables_for  # noqa: E402

NAME = structure.NAME
#: parameters_sulfur.ph_factor's kinds as integer codes (0 = no factor)
KIND_CODES = {"": 0, "acid": 1, "base": 2, "neutral_h2s": 3, "hs_anion": 4, "thiolate": 5, "hydroxide": 6}
DRIFT_COORDINATES = ("acid_yield", "arp_amine_pka")
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


# --------------------------------------------------------------------------- the Config
@dataclass(frozen=True)
class Settings:
    """The pH numerics a Config sets (every value ported), and integrate_sulfur's own defaults."""

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


def config_values(cfg=None) -> Tuple[Settings, str]:
    """The Config's pH numerics, after refusing every held key away from its shipped value (in the Config AND in the
    process, whose engine tables and ``core_parameters`` the build reads)."""
    from src.kinetic_core import config, sulfur

    cfg = config.Config.from_globals() if cfg is None else cfg
    shipped, here = config.Config.shipped(), config.Config.from_globals()
    for key in structure.HELD:
        how = structure.CONFIG_HANDLING[key][1]
        if cfg[key] != shipped[key]:
            raise UnsupportedConfiguration(f"{key} = {cfg[key]!r}: the JAX sulfur lane implements the shipped value "
                                           f"only ({how}); the rival stays NumPy-only")
        if here[key] != shipped[key]:
            raise UnsupportedConfiguration(f"{key}: the process has {here[key]!r} installed, and the build reads the "
                                           "engine's tables and core_parameters in-process; build in a shipped process")
    defaults = inspect.signature(sulfur.integrate_sulfur).parameters
    s = Settings(coupling=str(cfg["ph.coupling"]), node_grid=str(cfg["ph.node_grid"]),
                 fixed_point=str(cfg["ph.fixed_point"]), relaxation=float(cfg["ph.fixed_point.relaxation"]),
                 tolerance=float(cfg["ph.fixed_point.tolerance"]), max_sweeps=int(cfg["ph.fixed_point.max_sweeps"]),
                 flag=float(cfg["ph.fixed_point.flag"]), profile_node_step_k=float(cfg["sulfur.profile_node_step_k"]),
                 ph_nodes=int(defaults["ph_nodes"].default), ph_iterations=int(defaults["ph_iterations"].default))
    for value, allowed in ((s.coupling, COUPLINGS), (s.node_grid, NODE_GRIDS), (s.fixed_point, FIXED_POINTS)):
        if value not in allowed:
            raise UnsupportedConfiguration(f"{value!r}: not one of {allowed}")
    return s, cfg.digest()


# --------------------------------------------------------------------------- the charge balance (static half)
@dataclass(frozen=True)
class ChargeBalance:
    """``ph_state.titratable_inventory`` as a linear map of the state: per titratable group (row), its pK / dH / dCp
    at 25 C (padded to GROUP_WIDTH steps), z_max, its mmol/L per mmol/L of each species (``weights``; the lumped acid
    pool's share in ``acid_weights``, scaled by the ``acid_yield`` coordinate), where the ``arp_amine_pka`` coordinate
    sits (``arp_slot``), and the buffer's own rows (``buffer_rows``: group key -> row; the pot's mmol/L)."""

    names: Tuple[str, ...]
    pk25: np.ndarray
    dh: np.ndarray
    dcp: np.ndarray
    zmax: np.ndarray
    weights: np.ndarray
    acid_weights: np.ndarray
    arp_slot: np.ndarray
    buffer_rows: Mapping[str, int]
    const: Mapping[str, Any]

    def buffer_vector(self, spec) -> np.ndarray:
        out = np.zeros(len(self.names))
        for group, mmol in spec.groups():
            out[self.buffer_rows[group.key]] += float(mmol)
        return out


def _charge_balance(index: Mapping[str, int]) -> ChargeBalance:
    from src.kinetic_core import ph_state as ph

    n = len(index)
    rows: Dict[Any, List[Any]] = {}

    def row(key, group):
        if key not in rows:
            rows[key] = [group, np.zeros(n), np.zeros(n)]
        return rows[key]

    def add(group, species, weight=0.0, acid=0.0):
        r = row(("solute", group), group)
        r[1][index[species]] += weight
        r[2][index[species]] += acid

    for group, share in ph.ORGANIC_ACID_MIX:           # the lumped pool (x acid yield) and the carried carboxyls
        add(group, "ACID", acid=float(share))
        add(group, "CBX", weight=float(share))
    add(ph.FORMIC, "FA", 1.0)
    add(ph.ACETIC, "AA", 1.0)
    add(ph.cysteine_group(), "Cys", 1.0)
    template = ph.PhDrift(acid_yield=0.5, arp_amine_pka=7.0).arp_group()
    other = ph.PhDrift(acid_yield=0.5, arp_amine_pka=9.0).arp_group()
    slots = [i for i, (a, b) in enumerate(zip(template.pkas, other.pkas)) if a != b]
    if len(slots) != 1 or any(getattr(template, f) != getattr(other, f) for f in ("delta_h_kj_mol", "z_max")):
        raise UnsupportedConfiguration("PhDrift.arp_group moves more than one pK with arp_amine_pka: not ported")
    add(template, "ARP", 1.0)
    add(template, "CAMA", 1.0)
    add(ph.TTCA_GROUP, "TTCA", 1.0)
    buffer_rows: Dict[str, int] = {}
    kinds = [k for k in ph.BufferSpec._KINDS if k not in ("none", "clamped")]
    for kind in kinds:
        spec = ph.BufferSpec(kind=kind, phosphate_mol_l=1.0, citrate_mol_l=1.0, acetate_mol_l=1.0,
                             pyrophosphate_mol_l=1.0, phthalate_mol_l=1.0)
        for group, _ in spec.groups():
            row(("buffer", group.key), group)
    keys = list(rows)
    names, g = [], len(keys)
    pk = np.full((g, GROUP_WIDTH), PAD_PK)
    dh = np.zeros((g, GROUP_WIDTH))
    dcp = np.zeros((g, GROUP_WIDTH))
    zmax = np.zeros(g)
    w = np.zeros((g, n))
    wa = np.zeros((g, n))
    arp = np.zeros((g, GROUP_WIDTH), dtype=bool)
    for r, key in enumerate(keys):
        group, weights, acid = rows[key]
        k = len(group.pkas)
        if k > GROUP_WIDTH:
            raise UnsupportedConfiguration(f"group {group.key} has {k} protonation steps")
        pk[r, :k] = group.pkas
        dh[r, :k] = group.delta_h_kj_mol
        dcp[r, :k] = group.delta_cp_j_mol_k or (0.0,) * k
        zmax[r] = float(group.z_max)
        w[r], wa[r] = weights, acid
        if key[0] == "buffer":
            buffer_rows[key[1]] = r
            names.append(f"buffer_{key[1]}")
        else:
            names.append(group.key)
            if group.key == template.key:
                arp[r, slots[0]] = True
    const = {"R_KJ": float(ph.R_KJ), "pka_ref_t_k": float(ph.PKA_REFERENCE_T_K),
             "label_t_k": float(ph.LABEL_TEMPERATURE_K),
             "pkw_table_k": tuple(c + CELSIUS for c, _ in ph.PKW_TABLE_C), "pkw_table": tuple(p for _, p in ph.PKW_TABLE_C),
             "rho_table_c": tuple(c for c, _ in ph._WATER_DENSITY_C), "rho_table": tuple(r for _, r in ph._WATER_DENSITY_C)}
    return ChargeBalance(names=tuple(names), pk25=pk, dh=dh, dcp=dcp, zmax=zmax, weights=w, acid_weights=wa,
                         arp_slot=arp, buffer_rows=buffer_rows, const=const)


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
    """(pH, I) of an inventory at (T, SID): bisection at I = 0, one ionic sweep, Newton on the coupled (pH, I); the
    derivative by the implicit function (a 2 x 2 solve), so forward and reverse mode both see the exact root."""

    def __init__(self, balance: ChargeBalance):
        self.b = balance
        # numpy only: a jnp array made here (inside a trace) would be a tracer the custom JVP closes over
        self.zmax = np.asarray(balance.zmax, dtype=float)
        c = balance.const

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
        return jnp.where(jnp.asarray(self.b.arp_slot), theta["arp_amine_pka"], jnp.asarray(self.b.pk25))

    def temperature_terms(self, t_k, pk25):
        c = self.c
        pk0 = pk_at(pk25, jnp.asarray(self.b.dh), jnp.asarray(self.b.dcp), t_k, c["R_KJ"], c["pka_ref_t_k"])
        return {"pk0": pk0, "pkw": pkw_at(t_k, c), "a": davies_a(t_k, c)}

    def _group_pk(self, ionic_m, tt):
        af = davies(ionic_m, tt["a"])
        return tt["pk0"] + af * (2.0 * (jnp.asarray(self.zmax)[:, None] - jnp.arange(float(GROUP_WIDTH))) - 1.0), af

    def _group_terms(self, ph, pk):
        log_w = jnp.concatenate([jnp.zeros((pk.shape[0], 1)), jnp.cumsum(ph - pk, axis=1)], axis=1)
        w = 10.0 ** (log_w - jnp.max(log_w, axis=1, keepdims=True))
        w = w / jnp.sum(w, axis=1, keepdims=True)
        i = jnp.arange(GROUP_WIDTH + 1.0)
        z = jnp.asarray(self.zmax)
        return z - jnp.sum(i * w, axis=1), jnp.sum(w * (z[:, None] - i) ** 2, axis=1)

    def _residual_given(self, ph, pk, af, pkw, sid, conc):
        charge, ionic = self._group_terms(ph, pk)
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

    def inventory(self, y, theta, buffer):
        w = jnp.asarray(self.b.weights) + theta["acid_yield"] * jnp.asarray(self.b.acid_weights)
        return w @ jnp.maximum(y, 0.0) + buffer

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


# --------------------------------------------------------------------------- the rate constants
def _rate_static(derived: Mapping[str, Any], parameters: Mapping[str, Any]) -> Dict[str, Any]:
    """Per generic parameter key: k_ref, Ea (None = held), 1 / T_ref, its cap; per reaction: its parameter's index
    (-1 for a hand-ported step) and its pH-factor kind."""
    from src.kinetic_core import parameters as trunk_parameters
    from src.kinetic_core.parameters import KineticParameter
    from src.kinetic_core.parameters_furanic import TEMPERATURE_CAP_K

    keys = derived["generic_parameter_keys"]
    k_ref, ea, inv_t_ref, cap, linear = [], [], [], [], []
    for key in keys:
        p = parameters[key]
        if p.k_ref is None:
            raise UnsupportedConfiguration(f"{key}: unpopulated in the engine (it refuses to integrate too)")
        if isinstance(p, KineticParameter) and p.ea_kj_mol is None:
            raise UnsupportedConfiguration(f"{key}: a trunk parameter without Ea (the engine raises)")
        if float(p.k_ref) < 0.0:
            raise UnsupportedConfiguration(f"{key}: negative k_ref")
        k_ref.append(float(p.k_ref))
        ea.append(None if p.ea_kj_mol is None else float(p.ea_kj_mol))
        t_ref = float(trunk_parameters.T_REF_K) if isinstance(p, KineticParameter) else float(p.t_ref_k)
        inv_t_ref.append(1.0 / t_ref)
        cap.append(float(TEMPERATURE_CAP_K[key]) if key in TEMPERATURE_CAP_K else np.inf)
        linear.append(float(p.k_ref) == 0.0)
    pidx = np.array([keys.index(derived["parameter_key"][r]) if r not in structure.HAND_PORTED else -1
                     for r in derived["reactions"]])
    kind = np.array([KIND_CODES[derived["kind"][r]] for r in derived["reactions"]])
    return {"keys": tuple(keys), "k_ref": k_ref, "ea": ea, "inv_t_ref": np.asarray(inv_t_ref),
            "cap": np.asarray(cap), "linear": tuple(linear), "pidx": pidx, "kind": kind,
            "has_ea": np.asarray([e is not None for e in ea])}


def _term_constants() -> Dict[str, Any]:
    """The numbers the hand-ported terms read, from the engine."""
    from src.kinetic_core import parameters_sulfur as ps
    from src.kinetic_core.parameters import SCHIFF_AMADORI_SPLIT

    zheng = ps.ZHENG_CYSTEINE_THERMOLYSIS
    points = tuple(float(x) for x in zheng["ph_points"])
    if len(points) != 4 or float(ps.ZHENG_MECHANISM_BOUNDARY_PH) != points[1]:
        raise UnsupportedConfiguration("Zheng & Ho's table is no longer four pH points with the mechanism boundary "
                                       "at the second: the port's interpolation does not hold")
    if tuple(sorted(ps.GLUTATHIONE_ARRHENIUS_REFIT)) != (3.0, 5.0, 7.0, 9.0):
        raise UnsupportedConfiguration("the glutathione refit is no longer at pH 3 / 5 / 7 / 9")
    return {
        "R_KJ": float(ps.R_KJ), "ph_norm": float(ps.PH_NORMALISATION), "pka_ref_t_k": float(ps.PKA_REFERENCE_T_K),
        "pk_thiol": tuple(ps.PKA_GKL_CYSTEINE["thiol"]), "pk_amine": tuple(ps.PKA_GKL_CYSTEINE["cysteine_amine"]),
        "zheng_ph": points, "zheng_log10_a": tuple(math.log10(a * 60.0) for a in zheng["prefactor_per_second"]),
        "zheng_ea": tuple(float(e) for e in zheng["ea_kj_mol"]), "zheng_boundary": float(ps.ZHENG_MECHANISM_BOUNDARY_PH),
        "gsh_pairs": tuple(ps.GLUTATHIONE_ARRHENIUS_REFIT[p] for p in (3.0, 5.0, 7.0, 9.0)),
        "thi_share_points": tuple(ps.THIAMINE_H2S_SHARE_POINTS), "thi_share_floor": float(ps.THIAMINE_H2S_SHARE_FLOOR_PH),
        "adduct_k_ref": float(ps.THIOL_QUINONE_EQUILIBRIUM_USED_L_PER_MMOL),
        "adduct_dh": float(ps.STACK_DELTA_H_ADDUCT_KJ_MOL), "adduct_t_k": float(ps.STACK_MEASUREMENT_T_K),
        "amadori_factor": float(SCHIFF_AMADORI_SPLIT["ratio_amadori_over_schiff_pseudo_first_order"])
                          * float(SCHIFF_AMADORI_SPLIT["amine_loading_mmol_L_for_the_ratio"]),
    }


def _speciation_pk(name, t_k, ionic_m, c, a):
    """parameters_sulfur.pka_at under B51 (Suleimenov H2S, GKL cysteine pair), minus Davies."""
    if name == "h2s_1":
        pk = -(782.43945 + 0.361261 * t_k - 1.6722e-4 * t_k * t_k - 20565.7315 / t_k - 142.741722 * jnp.log(t_k))
    else:
        pk25, dh, dcp = c["pk_thiol"] if name == "thiol" else c["pk_amine"]
        pk = pk_at(pk25, dh, dcp, t_k, c["R_KJ"], c["pka_ref_t_k"])
    return pk - davies(ionic_m, a)


def ph_factors(ph, t_k, ionic_m, c, a):
    """parameters_sulfur.ph_factor for every kind (indexed by KIND_CODES); the normaliser at pH 5 and I = 0."""
    p5 = c["ph_norm"]

    def amine(x, i):
        return 1.0 / (1.0 + 10.0 ** (_speciation_pk("cysteine_amine", t_k, i, c, a) - x))

    def neutral(x, i):
        return 1.0 / (1.0 + 10.0 ** (x - _speciation_pk("h2s_1", t_k, i, c, a)))

    def thiolate(x, i):
        return 1.0 / (1.0 + 10.0 ** (_speciation_pk("thiol", t_k, i, c, a) - x))

    return jnp.stack([
        jnp.asarray(1.0), 10.0 ** (p5 - ph), amine(ph, ionic_m) / amine(p5, 0.0),
        neutral(ph, ionic_m) / neutral(p5, 0.0), (1.0 - neutral(ph, ionic_m)) / (1.0 - neutral(p5, 0.0)),
        thiolate(ph, ionic_m) / thiolate(p5, 0.0), 10.0 ** (ph - p5)])


def cysteine_thermolysis_k(ph_label, t_k, c):
    pts, la, ea = c["zheng_ph"], c["zheng_log10_a"], c["zheng_ea"]
    x = jnp.minimum(jnp.maximum(ph_label, pts[1]), pts[-1])
    lo = x <= pts[2]
    w = jnp.where(lo, (x - pts[1]) / (pts[2] - pts[1]), (x - pts[2]) / (pts[3] - pts[2]))
    log_a = jnp.where(lo, la[1] + w * (la[2] - la[1]), la[2] + w * (la[3] - la[2]))
    e = jnp.where(lo, ea[1] + w * (ea[2] - ea[1]), ea[2] + w * (ea[3] - ea[2]))
    below = ph_label < c["zheng_boundary"]
    log_a, e = jnp.where(below, la[0], log_a), jnp.where(below, ea[0], e)
    return 10.0 ** log_a * jnp.exp(-e / (c["R_KJ"] * t_k))


def glutathione_thermolysis_k(ph_label, t_k, c):
    (l3, e3), (l5, e5), (l7, e7), (l9, e9) = c["gsh_pairs"]
    x = jnp.minimum(ph_label, 9.0)
    lo = x <= 7.0
    w = jnp.where(lo, (x - 5.0) / 2.0, (x - 7.0) / 2.0)
    ln_a = jnp.where(lo, l5 + w * (l7 - l5), l7 + w * (l9 - l7))
    e = jnp.where(lo, e5 + w * (e7 - e5), e7 + w * (e9 - e7))
    below = ph_label < c["zheng_boundary"]
    ln_a, e = jnp.where(below, l3, ln_a), jnp.where(below, e3, e)
    return jnp.exp(ln_a) * jnp.exp(-e / (c["R_KJ"] * t_k))


def thiamine_h2s_share(ph_label, c):
    (p0, f0), (p1, f1) = c["thi_share_points"]
    mid = 10.0 ** (math.log10(f0) + (ph_label - p0) / (p1 - p0) * (math.log10(f1) - math.log10(f0)))
    return jnp.where(ph_label < c["thi_share_floor"], 0.0, jnp.where(ph_label >= p1, f1, mid))


def _make_rate_constants(static, c, balance_const, reaction_index):
    keys, linear = static["keys"], static["linear"]
    has_ea = jnp.asarray(static["has_ea"])
    inv_t_ref = jnp.asarray(static["inv_t_ref"])
    cap = jnp.asarray(static["cap"])
    pidx = static["pidx"]
    safe = jnp.asarray(np.where(pidx >= 0, pidx, 0))
    generic = jnp.asarray(pidx >= 0)
    kind = jnp.asarray(static["kind"])
    ix = reaction_index
    off = jnp.asarray([ix[r] for r in structure.OFF])
    r_kj = c["R_KJ"]

    def rate_constants(theta: Mapping[str, Any], ctx: RateCtx):
        t_k = jnp.asarray(ctx.t_k, dtype=jnp.float64)
        ph = jnp.asarray(ctx.ph, dtype=jnp.float64)
        ionic = jnp.asarray(ctx.ionic_m, dtype=jnp.float64)
        label = jnp.where(jnp.isnan(ctx.ph_label), ph, ctx.ph_label)     # B51: tables read the cooled pH
        k_ref = jnp.stack([jnp.asarray(theta[f"k_ref:{k}"], dtype=jnp.float64) if lin
                           else 10.0 ** jnp.asarray(theta[f"log10_k:{k}"], dtype=jnp.float64)
                           for k, lin in zip(keys, linear)])
        ea = jnp.stack([jnp.asarray(theta[f"ea:{k}"], dtype=jnp.float64) if e is not None else jnp.asarray(0.0)
                        for k, e in zip(keys, static["ea"])])
        t_eff = jnp.minimum(t_k, cap)
        k_par = jnp.where(has_ea, k_ref * jnp.exp(-(ea / r_kj) * (1.0 / t_eff - inv_t_ref)), k_ref)
        a = davies_a(t_k, balance_const)
        k = jnp.where(generic, k_par[safe], 0.0) * ph_factors(ph, t_k, ionic, c, a)[kind]
        share = thiamine_h2s_share(label, c)
        k_adduct = c["adduct_k_ref"] * jnp.exp(-(c["adduct_dh"] / r_kj) * (1.0 / t_k - 1.0 / c["adduct_t_k"]))
        derived = {
            "r_cys_h2s": cysteine_thermolysis_k(label, t_k, c),
            "r_gsh_h2s": glutathione_thermolysis_k(label, t_k, c),
            "r_thi_h2s": jnp.where(share > 0.0, share / (1.0 - share), 0.0) * (k[ix["r_thi_hmp"]]
                                                                               + k[ix["r_thi_mesh"]]),
            "ch_thioether_release": k[ix["ch_thioether_mft"]] / k_adduct,
            "ch_thioether_release_fft": k[ix["ch_thioether_fft"]] / k_adduct,
            "r_amadori": c["amadori_factor"] * k[ix["r_schiff"]],
        }
        names = sorted(derived)
        k = k.at[jnp.asarray([ix[n] for n in names])].set(jnp.stack([derived[n] for n in names]))
        return k.at[off].set(0.0)

    return rate_constants


# --------------------------------------------------------------------------- the module
def _bands(module_params: Mapping[str, float], parameters: Mapping[str, Any]):
    """Declared bands: the shipped wave's fit bounds (the Laplace file's ``bounds``) for the fitted log10 k_ref at
    145 C and the drift; an Ea's printed CI95 where it has one; None otherwise (no band invented here)."""
    from src.kinetic_core import ph_state as ph
    from src.kinetic_core import uncertainty

    bands: Dict[str, Optional[Tuple[float, float]]] = {k: None for k in module_params}
    basis: Dict[str, str] = {}
    lap = uncertainty._laplace() or {"coordinates": [], "bounds": []}
    fitted = {c["key"]: tuple(b) for c, b in zip(lap["coordinates"], lap["bounds"])
              if c.get("block") == "log10_k_ref_at_145C"}
    for name in module_params:
        kind, _, key = name.partition(":")
        p = parameters.get(key)
        src = f"{getattr(p, 'evidence_class', '')}; {str(getattr(p, 'source_anchor', ''))[:80]}" if p else ""
        if kind == "log10_k":
            t_ref = getattr(p, "t_ref_k", None)
            if key in fitted and t_ref is not None and abs(float(t_ref) - 418.15) < 1e-9:
                bands[name] = fitted[key]
                basis[name] = f"log10 k_ref of {key} at 145 C; band = the shipped wave's fit bounds ({src})"
            else:
                basis[name] = f"log10 k_ref of {key} at its T_ref; no band declared ({src})"
        elif kind == "k_ref":
            basis[name] = f"k_ref of {key}, LINEAR because the engine carries it at exactly zero ({src})"
        elif kind == "ea":
            ci = getattr(p, "ea_ci95_kj_mol", None)
            if ci is not None:
                bands[name] = (module_params[name] - float(ci), module_params[name] + float(ci))
                basis[name] = f"Ea of {key}; band = its printed 95 % CI ({src})"
            else:
                basis[name] = f"Ea of {key}; no band declared ({src})"
    bands["acid_yield"] = tuple(float(x) for x in ph.ACID_YIELD_BOUNDS)
    bands["arp_amine_pka"] = tuple(float(x) for x in ph.ARP_AMINE_PKA_BOUNDS)
    basis["acid_yield"] = "engine.core_ph_drift() acid yield per sink event (B2.2 frozen); band ph_state.ACID_YIELD_BOUNDS"
    basis["arp_amine_pka"] = ("engine.core_ph_drift() ARP secondary-ammonium pK (B2.2 frozen); band "
                              "ph_state.ARP_AMINE_PKA_BOUNDS")
    return bands, basis


@dataclass(frozen=True)
class SulfurLane:
    """The sulfur Module, its one-module Network, the charge balance and the Config's pH numerics."""

    module: Module
    network: Network
    balance: ChargeBalance
    settings: Settings
    parameters: Mapping[str, Any]

    @property
    def species(self):
        return self.network.species

    @property
    def index(self):
        return self.network.index

    def theta0(self):
        return self.network.theta0()


def build_module(cfg=None) -> Tuple[Module, ChargeBalance, Settings, Dict[str, Any]]:
    """The sulfur Module for an engine Config (shipped values only on the held keys)."""
    from src.kinetic_core import engine

    settings, digest = config_values(cfg)
    derived = structure.derive()
    problems = [p for p in structure.structure_problems(numeric=False)]
    if problems:
        raise UnsupportedConfiguration("the NumPy sulfur lane moved under the port:\n" + "\n".join(problems))
    parameters = dict(engine.core_parameters(engine.SULFUR))
    static = _rate_static(derived, parameters)
    c = _term_constants()
    index = {s: i for i, s in enumerate(derived["species"])}
    balance = _charge_balance(index)
    rix = {r: j for j, r in enumerate(derived["reactions"])}
    rate_constants = _make_rate_constants(static, c, balance.const, rix)
    params: Dict[str, float] = {}
    for key, lin, e in zip(static["keys"], static["linear"], static["ea"]):
        if lin:
            params[f"k_ref:{key}"] = 0.0
        else:
            params[f"log10_k:{key}"] = math.log10(float(parameters[key].k_ref))
        if e is not None:
            params[f"ea:{key}"] = float(e)
    drift = engine.core_ph_drift()
    params["acid_yield"] = float(drift.acid_yield)
    params["arp_amine_pka"] = float(drift.arp_amine_pka)
    bands, basis = _bands(params, parameters)
    hand = sorted(structure.HAND_PORTED)
    module = Module(
        name=NAME, species=derived["species"], elements=derived["elements"], reactions=derived["reactions"],
        stoichiometry=derived["stoichiometry"], reactants=derived["reactants"], params=params, bands=bands,
        basis=basis, rate_constants=rate_constants, observables=observables_for(derived["species"]),
        needs_ph=True, config_digest=digest,
        notes=(f"pH coupling {settings.coupling}; node grid {settings.node_grid}; fixed point {settings.fixed_point} "
               f"(relaxation {settings.relaxation}, tolerance {settings.tolerance}, max sweeps {settings.max_sweeps})",
               f"hand-ported reactions ({len(hand)}): {hand}",
               f"held at zero (route off at the shipped configuration): {list(structure.OFF)}",
               f"linear (zero k_ref) parameters: {[k for k, lin in zip(static['keys'], static['linear']) if lin]}",
               f"titratable groups: {list(balance.names)}"),
        breakpoints={"ph_label": (float(c["zheng_boundary"]), float(c["zheng_ph"][2]), float(c["zheng_ph"][3]), 7.0,
                                  9.0, float(c["thi_share_floor"]), float(c["thi_share_points"][1][0]))})
    return module, balance, settings, parameters


def build(cfg=None, *, solver: str = "Kvaerno5") -> SulfurLane:
    module, balance, settings, parameters = build_module(cfg)
    return SulfurLane(module=module, network=compose([module], solver=solver), balance=balance, settings=settings,
                      parameters=parameters)


def theta_from(lane: SulfurLane, parameters: Optional[Mapping[str, Any]] = None, drift=None) -> Dict[str, Any]:
    """theta for an engine parameter dict (``k_ref``, ``ea_kj_mol``) and a PhDrift, the rest at shipped values."""
    out = dict(lane.module.params)
    for name in out:
        kind, _, key = name.partition(":")
        if parameters is None or key not in parameters:
            continue
        p = parameters[key]
        if kind == "log10_k":
            if not float(p.k_ref) > 0.0:
                raise UnsupportedConfiguration(f"{key}: k_ref {p.k_ref!r} on a log10 coordinate")
            out[name] = math.log10(float(p.k_ref))
        elif kind == "k_ref":
            out[name] = float(p.k_ref)
        elif kind == "ea":
            out[name] = float(p.ea_kj_mol)
    if drift is not None:
        out["acid_yield"], out["arp_amine_pka"] = float(drift.acid_yield), float(drift.arp_amine_pka)
    return {NAME: {k: jnp.asarray(float(v)) for k, v in out.items()}}


# --------------------------------------------------------------------------- the pot
def node_times(settings: Settings, duration: float, profile: Optional[Tuple[np.ndarray, np.ndarray]] = None):
    """sulfur.ph_node_times (and, along a profile, profile_node_times at the Config's step)."""
    n = max(int(settings.ph_nodes), 5)
    linear = np.linspace(0.0, duration, n)
    nodes = linear if settings.node_grid == "uniform" else np.unique(np.concatenate(
        [linear, np.geomspace(duration * 1e-5, duration, n)]))
    if profile is not None:
        times, temps = profile
        out = [np.array([0.0, duration])]
        for t0, t1, c0, c1 in zip(times[:-1], times[1:], temps[:-1], temps[1:]):
            k = max(1, int(np.ceil(abs(c1 - c0) / settings.profile_node_step_k - 1e-12)))
            out.append(np.linspace(t0, t1, k + 1))
        nodes = np.unique(np.concatenate([nodes, np.unique(np.clip(np.concatenate(out), 0.0, duration))]))
    return nodes


def engine_prelude(initial: Mapping[str, float], process) -> Tuple[Dict[str, float], Dict[str, float]]:
    """``engine._integrate_program``'s sulfur prelude: ambient OX, PROT_SS from the matrix sites, the headspace
    reservoir OXR, OXV; MET / PRO / PYRL recorded and dropped. Returns (state, dropped)."""
    from src.kinetic_core import engine
    from src.kinetic_core.matrix_sites import resolve

    state = {k: float(v) for k, v in initial.items()}
    state.setdefault("OX", float(engine.OX_AMBIENT_MMOL_L))
    if "PROT_SS" not in state:
        charged, _ = resolve(process)
        if charged is not None and charged.disulfide > 0:
            state["PROT_SS"] = float(charged.disulfide)
    if "OXR" not in state:
        state["OXR"] = float(engine.oxygen_reservoir_units(process)[0])
    state.setdefault("OXV", 0.0)
    dropped = {k: state.pop(k) for k in ("MET", "PRO", "PYRL") if k in state}
    return state, dropped


def make_pot(lane: SulfurLane, initial: Mapping[str, float], process, *, prelude: bool = True) -> Dict[str, Any]:
    """A pot from the engine's ``ProcessSpec`` (thermal programme or profile, declared label pH, buffer) and its
    charged species. Arrays: y0 (n,), label_ph, buffer (G,), durations (S,), node_s / node_tc (S, N) (the pH nodes as
    fractions of each segment and their temperatures), knot_s / knot_tc (S, K) (the temperature course inside a
    segment, for the inline coupling), and ``mode`` ("dynamic" | "clamped", static)."""
    from src.kinetic_core import ph_state as ph

    if getattr(process, "ph_drift", None) is not None:
        raise UnsupportedConfiguration("a pot-level ph_drift: set the module's acid_yield / arp_amine_pka instead")
    state, dropped = engine_prelude(initial, process) if prelude else ({k: float(v) for k, v in initial.items()}, {})
    unknown = sorted(k for k in state if k not in lane.index)
    if unknown:
        raise UnsupportedConfiguration(f"the pot charges {unknown}, which the sulfur lane does not carry")
    y0 = np.zeros(lane.network.n_species)
    for k, v in state.items():
        if float(v) < 0.0:
            raise ValueError(f"negative initial concentration for {k!r}")
        y0[lane.index[k]] = float(v)
    spec = process.buffer if getattr(process, "buffer", None) is not None else ph.DEFAULT_BUFFER
    mode = "clamped" if spec.is_clamped else "dynamic"
    thermal = process.thermal
    if mode == "clamped" and thermal.profile is not None:
        raise UnsupportedConfiguration("a clamped pH along a temperature profile (the engine's node-interpolated "
                                       "constants at T(t)) is not ported")
    if mode == "clamped" and process.ph_final is not None and float(process.ph_final) != float(process.ph):
        raise UnsupportedConfiguration("a prescribed pH line between two measured endpoints (B2.1) is not ported")
    if thermal.profile is not None:
        times = np.array([t for t, _ in thermal.profile], dtype=float)
        temps = np.array([c for _, c in thermal.profile], dtype=float)
        duration = float(thermal.total_minutes)
        nodes = node_times(lane.settings, duration, (times, temps + CELSIUS))
        segments = [(duration, nodes, np.interp(nodes, times, temps), times, temps)]
    else:
        segments = []
        for d, t_c in thermal.segments:
            nodes = node_times(lane.settings, float(d))
            segments.append((float(d), nodes, np.full(len(nodes), float(t_c)), np.array([0.0, float(d)]),
                             np.array([float(t_c), float(t_c)])))
    widths = {len(s[1]) for s in segments}
    if len(widths) != 1:
        raise UnsupportedConfiguration("segments with different node counts (a profile inside a staircase)")
    k_width = max(len(s[3]) for s in segments)

    def pad(a, n):
        return np.concatenate([a, np.full(n - len(a), a[-1])])

    return {"y0": y0, "label_ph": float(process.ph), "buffer": lane.balance.buffer_vector(spec), "mode": mode,
            "durations": np.array([s[0] for s in segments]),
            "node_s": np.stack([s[1] / s[0] for s in segments]), "node_tc": np.stack([s[2] for s in segments]),
            "knot_s": np.stack([pad(s[3] / s[0], k_width) for s in segments]),
            "knot_tc": np.stack([pad(s[4], k_width) for s in segments]),
            "label_nudge": np.zeros((len(segments), len(segments[0][1]))),
            "node_minutes": [s[1] for s in segments], "dropped": dropped, "buffer_spec": spec}


def pot_arrays(pot: Mapping[str, Any]) -> Dict[str, Any]:
    """The numeric part of a pot (the compiled functions' argument). ``label_nudge`` (S, N; zero) is added to the
    nodes' cooled pH where the measured tables read it -- only a parity check sets it, to put a node whose cooled pH
    sits ON a table's jump (a pot declared at exactly label pH 5.0, Zheng & Ho's boundary) on the engine's side of
    the jump, which there is decided by the root-finders' last digits."""
    return {k: jnp.asarray(pot[k], dtype=jnp.float64) for k in ("y0", "label_ph", "buffer", "durations", "node_s",
                                                                  "node_tc", "knot_s", "knot_tc", "label_nudge")}


# --------------------------------------------------------------------------- integration
class SulfurRun(NamedTuple):
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


def _solve(net: Network, rhs_k, data, y0, duration, save_s, jumps, *, rtol, atol, solver, adjoint, max_steps,
           throw):
    """One solve on s in [0, 1] (the right-hand side scaled by ``duration``), the rate constants ``rhs_k(s, y, data)``.
    Every traced number reaches the vector field through diffrax's ``args`` (``data``, ``duration``), never by
    closure, so the discrete adjoint can differentiate it."""
    import diffrax

    stoich = jnp.asarray(net.stoichiometry)

    def f(s, y, args):
        data_, duration_ = args
        return duration_ * (stoich @ (rhs_k(s, y, data_) * mass_action(y, net)))

    controller = diffrax.ClipStepSizeController(diffrax.PIDController(rtol=rtol, atol=atol), jump_ts=jumps)
    sol = diffrax.diffeqsolve(diffrax.ODETerm(f), solver, t0=0.0, t1=1.0, dt0=None, y0=y0, args=(data, duration),
                              saveat=diffrax.SaveAt(ts=save_s), stepsize_controller=controller,
                              max_steps=max_steps, adjoint=adjoint, throw=throw)
    return sol.ys, sol.stats["num_steps"], sol.result != diffrax.RESULTS.successful


def integrate(lane: SulfurLane, theta, pot: Mapping[str, Any], *, rtol: float = 1e-8, atol: float = 1e-14,
              solver=None, adjoint=None, reverse: bool = False, max_steps: int = 200_000,
              throw: bool = True, force_sweeps: Optional[int] = None) -> SulfurRun:
    """One pot through its programme with the Config's pH coupling (``pot["mode"]`` "clamped" holds the declared pH).
    ``reverse=True`` runs the fixed point's sweeps as a scan (``jax.grad``); otherwise a while loop (forward mode,
    ``vmap``). ``adjoint`` defaults to diffrax's discrete adjoint (``diffrax.ForwardMode()`` for ``jax.jacfwd``).
    ``force_sweeps`` replaces the stopping rule by a fixed sweep count (every segment): for a parity check on a pot
    whose last residual sits at the tolerance, where solver noise alone decides one sweep more or less."""
    import diffrax

    net, bal, st = lane.network, PhSolver(lane.balance), lane.settings
    solver = solver or getattr(diffrax, net.solver)()
    adjoint = adjoint or diffrax.RecursiveCheckpointAdjoint()
    a = pot_arrays(pot)
    th = theta[NAME]
    mode = pot.get("mode", "dynamic")
    kw = dict(rtol=rtol, atol=atol, solver=solver, adjoint=adjoint, max_steps=max_steps, throw=throw)
    y0 = a["y0"]
    sid = bal.strong_ion_difference(a["label_ph"], y0, th, a["buffer"])

    def rates(t_k, ph, ionic, cool):
        return net.rate_constants(theta, RateCtx(t_c=t_k - CELSIUS, t_k=t_k, ph=ph, ph_label=cool, ionic_m=ionic))

    def node_states(ys, node_tk):
        return jax.vmap(lambda y, t: bal.state(y, t, sid, th, a["buffer"]))(ys, node_tk)

    min_sweeps = max(int(st.ph_iterations), 1)
    relaxed = st.fixed_point == "relaxed_converged"
    max_sweeps = max(min_sweeps, int(st.max_sweeps)) if relaxed else min_sweeps
    if force_sweeps is not None:
        if not 1 <= int(force_sweeps) <= max_sweeps:
            raise ValueError(f"force_sweeps must lie in [1, {max_sweeps}]")
        min_sweeps = max_sweeps = int(force_sweeps)

    def segment(y_start, seg):
        duration, node_s, node_tc, knot_s, knot_tc, nudge = seg
        node_tk = node_tc + CELSIUS
        if mode == "clamped":
            t_k = knot_tc[0] + CELSIUS
            k0 = rates(t_k, a["label_ph"], 0.0, jnp.nan)
            ys, steps, bad = _solve(net, lambda s, y, k: k, k0, y_start, duration, node_s, node_s, **kw)
            ys = jnp.maximum(ys, 0.0)
            n = node_s.shape[0]
            return ys[-1], (ys, jnp.full(n, a["label_ph"]), jnp.zeros(n), jnp.full(n, a["label_ph"]),
                            jnp.asarray(0), jnp.asarray(0.0), steps, bad)
        if st.coupling == "inline":
            def rhs_k(s, y, d):
                theta_, sid_, buffer_, knot_s_, knot_tc_ = d
                t_k = jnp.interp(s, knot_s_, knot_tc_) + CELSIUS
                ph, ionic, cool = bal.state(y, t_k, sid_, theta_[NAME], buffer_)
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
        data = (theta, sid, a["buffer"], y_start, duration, node_s, node_tk, nudge)

        def sweep(carry, d):
            theta_, sid_, buffer_, y_start_, duration_, node_s_, node_tk_, nudge_ = d
            count, node_ph, node_i, node_cool, _ys, _res, _done, steps, bad = carry
            k_nodes = jax.vmap(lambda t, p, i, c: net.rate_constants(
                theta_, RateCtx(t_c=t - CELSIUS, t_k=t, ph=p, ph_label=c, ionic_m=i)))(node_tk_, node_ph, node_i,
                                                                                         node_cool + nudge_)
            ys, n_steps, failed = _solve(net, lambda s, y, dd: _interp_rows(s, dd[0], dd[1]), (node_s_, k_nodes),
                                         y_start_, duration_, node_s_, node_s_, **kw)
            ys = jnp.maximum(ys, 0.0)
            fresh, fi, fc = jax.vmap(lambda y, t: bal.state(y, t, sid_, theta_[NAME], buffer_))(ys, node_tk_)
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
        if reverse:
            def body(carry, _):
                return jax.lax.cond(carry[6], lambda c, d: c, sweep, carry, data), None

            out, _ = jax.lax.scan(body, init, None, length=max_sweeps)
        else:
            out = jax.lax.while_loop(lambda c: ~c[6], lambda c: sweep(c, data), init)
        count, node_ph, node_i, node_cool, ys, res, _, steps, bad = out
        return ys[-1], (ys, node_ph, node_i, node_cool, count, res, steps, bad)

    segs = (a["durations"], a["node_s"], a["node_tc"], a["knot_s"], a["knot_tc"], a["label_nudge"])
    _, (ys, ph, ionic, cool, count, res, steps, bad) = jax.lax.scan(segment, y0, segs)
    return SulfurRun(ys=ys, ph=ph, ionic=ionic, ph_cool=cool, sweeps=count, residual=res, num_steps=steps,
                     failed=bad, sid=sid)


def compiled(lane: SulfurLane, pot_template: Mapping[str, Any], **kw):
    """``f(theta, arrays) -> SulfurRun``, jit-compiled once per lane, mode and pot shape."""
    mode = pot_template.get("mode", "dynamic")
    return jax.jit(lambda theta, arrays: integrate(lane, theta, dict(arrays, mode=mode), **kw))


def observe(lane: SulfurLane, run: SulfurRun, species: Sequence[str], rows: Sequence[int], *,
            floor: float = 1e-300):
    """log10 mmol/L of ``species`` at node ``rows`` of the LAST segment, then the cooled pH at the end."""
    ys = run.ys[-1]
    idx = [lane.index[s] for s in species]
    vals = [jnp.log10(jnp.maximum(ys[r, i], floor)) for i in idx for r in rows]
    vals.append(run.ph_cool[-1, -1])
    return jnp.stack(vals)


def jacobian(lane: SulfurLane, pot_template, coordinates: Sequence[str], species, rows, **kw):
    """``f(theta, arrays) -> (observables, d observables / d coordinates)``, exact, forward mode."""
    import diffrax

    mode = pot_template.get("mode", "dynamic")
    kw = dict(kw, adjoint=diffrax.ForwardMode())

    def obs(x, theta, arrays):
        th = {NAME: dict(theta[NAME], **{c: x[i] for i, c in enumerate(coordinates)})}
        v = observe(lane, integrate(lane, th, dict(arrays, mode=mode), **kw), species, rows)
        return v, v

    def both(theta, arrays):
        x = jnp.stack([theta[NAME][c] for c in coordinates])
        jac, value = jax.jacfwd(obs, has_aux=True)(x, theta, arrays)
        return value, jac

    return jax.jit(both)
