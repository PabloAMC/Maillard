"""
src/kinetic_core_jax/sulfur_lane.py -- the sulfur lane (``sulfur.FULL_REACTIONS``: the trunk's 26 steps, the fed
3-deoxy triangle and the 103 sulfur steps; 62 species) as A2 contract modules, with its dynamic pH (plan step
A2 + H3, docs/guides/A2_MODULE_CONTRACT.md v2.3 sec. 6.2.2; 2026-11-05, onto the network's own pH solve 2026-11-06).

    lane = sulfur_lane.build(cfg)                    # ONE module "sulfur" + its Network (the parity form)
    lane = sulfur_lane.build_composed(cfg)           # the SAME lane as three composed modules (sec. 6.2.2)
    pot = sulfur_lane.make_pot(lane, initial, process)   # the engine's sulfur prelude, nodes, SID inputs
    run = sulfur_lane.integrate(lane, theta, pot)        # = network.integrate(lane.network, ...).ph

TWO FORMS, ONE CODE PATH. The engine's sulfur lane is one recipe with its OWN copy of the trunk (its 26 steps run on
``core_parameters("sulfur")``, under ``sulfur_rate_constants_at``'s generic rule, beside the lane's pH machinery).
``build`` ports it as one module, ``sulfur``, owning all 62 species (PR #83's form; the parity test runs it).
``build_composed`` builds the same reactions as three modules -- ``sulfur_trunk`` (the trunk's 26 steps under the
lane's rule; owns Glc, Gly, FA, AA, ...), ``sulfur_fed_3deoxy`` (3 steps; owns DGAL) and ``sulfur`` (103 steps; owns
the rest, declares its consumption of the trunk's species and the lane's whole titratable table) -- composed by
``network.compose`` with declared competition refused. Both forms come from ``_block_module``, so their rate constants
are the same arithmetic. Composing ``sulfur`` with the DECLARED-pH trunk (``trunk.py``, the acrylamide lane's) is
A2b's union network, which changes numbers; ``compose`` refuses it (two pH conventions in one ``RateCtx.ph_label``).

DERIVED, NEVER RETYPED (``sulfur_structure.derive`` and ``build``; an edit to the engine's data carries over):
species and C / N / S, the 132 reactions in FULL_REACTIONS order, stoichiometry and mass-action orders, each
reaction's parameter key and pH-factor kind (``REACTION_PH_FACTOR``), every parameter's k_ref, Ea (None = held), T_ref
(``KineticParameter``: ``parameters.T_REF_K``; ``SulfurParameter``: its ``t_ref_k``) and temperature cap
(``parameters_furanic.TEMPERATURE_CAP_K``), the gas constant, the speciation pK tables (``PKA_GKL_CYSTEINE``,
``PKA_REFERENCE_T_K``, ``PH_NORMALISATION``), Zheng & Ho's cysteine set and its mechanism boundary, the glutathione
refit, the thiamine share points and floor, Stack's adduct equilibrium and dH, the Amadori split, the titratable
groups (pK, dH, dCp, z_max) of ``ph_state`` and the buffer groups (from ``BufferSpec.groups``), their centres
(``ph_state.TITRATABLE_CENTRES``) and the lane's ``sulfur.CENTRE_LEDGER``, the water density and pKw tables, the label
temperature, the pot prelude's constants (``OX_AMBIENT_MMOL_L``, ``oxygen_reservoir_units``, the matrix sites), the
drift calibration (``engine.core_ph_drift``) and ``integrate_sulfur``'s node and sweep defaults (read off its
signature).

HAND-PORTED (each pinned: ``sulfur_structure.SOURCE_PINS`` hashes the engine function's AST, and
tests/scientific/test_jax_sulfur_parity.py compares the numbers):

* the generic rate rule ``k_at(min(T, cap)) x ph_factor(kind, pH, T, I)`` and the six pH-factor kinds with B51's
  speciation (Suleimenov & Seward's H2S pK1 polynomial, GKL's cysteine pair with dCp, Davies on the numerator only);
* the derived constants (``sulfur_structure.HAND_PORTED``): Zheng & Ho's cysteine and glutathione thermolysis on the
  label-scale pH, the thiamine H2S share (B50), the thioether release k_forward / K_adduct(T), the Amadori pin; and
  the 21 steps the shipped configuration holds at zero (their routes are off; the port refuses them on);
* the titratable inventory (``ph_state.titratable_inventory``) as ``network.Titratable`` pools: the organic-acid mix
  fed by CBX at full strength and by ACID x the ``acid_yield`` coordinate, FA, AA, cysteine, the ARP class (ARP and
  CAMA, its amine pK the ``arp_amine_pka`` coordinate), TTCA, and the buffer groups as pot-supplied pools. The charge
  balance, its (pH, I) root, the couplings ("fixed_point_nodes", shipped, and "inline"), the clamped pH and the node
  grids are the NETWORK's (``ph_coupling``), moved there from this file unchanged;
* the engine's sulfur prelude of ``_integrate_program`` (ambient OX, the protein disulfide pool, the headspace
  reservoir, OXV, MET / PRO / PYRL recorded and dropped).

THE COORDINATES (``Module.params``): ``log10_k:<key>`` (log10 of k_ref at its T_ref) or, where the engine carries the
key at exactly zero, a LINEAR ``k_ref:<key>``; ``ea:<key>`` where the key has an Ea (a held, Ea-less key has none);
and the drift's ``acid_yield`` and ``arp_amine_pka`` (``engine.core_ph_drift``), which enter the charge balance only.

SWITCHES (``sulfur_structure.CONFIG_HANDLING``): the pH numerics are ported at every value; every rival route, B51
off, the trunk pH terms and the trunk solutes in the balance are HELD at their shipped values and refused otherwise
(they stay NumPy-only while M1b(i) / M1c reshape them). A pot-level ``ph_drift`` and a prescribed pH line (no drift)
are refused at the pot.
"""
from __future__ import annotations

import inspect
import math
from dataclasses import dataclass
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)

from . import ph_coupling  # noqa: E402
from . import sulfur_structure as structure  # noqa: E402
from .lane import UnsupportedConfiguration  # noqa: E402
from .network import CELSIUS, Module, Network, RateCtx, Titratable, compose  # noqa: E402
from .ph_coupling import (COUPLINGS, FIXED_POINTS, FIXED_SWEEPS_STOP, NODE_GRIDS, PhModel, PhRun,  # noqa: E402,F401
                          PhSolver, davies, davies_a, jump_grid, node_times, pk_at)
from .trunk import observables_for  # noqa: E402

NAME = structure.NAME
#: the composed form's modules, in order (each requires the ones before it)
COMPOSED = ("sulfur_trunk", "sulfur_fed_3deoxy", NAME)
#: parameters_sulfur.ph_factor's kinds as integer codes (0 = no factor)
KIND_CODES = {"": 0, "acid": 1, "base": 2, "neutral_h2s": 3, "hs_anion": 4, "thiolate": 5, "hydroxide": 6}
DRIFT_COORDINATES = ("acid_yield", "arp_amine_pka")
#: the hand-ported terms that read another reaction's constant: term -> the reactions it reads (same module only)
TERM_READS = {"r_thi_h2s": ("r_thi_hmp", "r_thi_mesh"), "ch_thioether_release": ("ch_thioether_mft",),
              "ch_thioether_release_fft": ("ch_thioether_fft",), "r_amadori": ("r_schiff",),
              "r_cys_h2s": (), "r_gsh_h2s": ()}
#: PR #83's name for the pH numerics (now the network's PhModel)
Settings = PhModel
SulfurRun = PhRun


# --------------------------------------------------------------------------- the Config
def config_values(cfg=None) -> Tuple[Dict[str, Any], str]:
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
    s = dict(coupling=str(cfg["ph.coupling"]), node_grid=str(cfg["ph.node_grid"]),
             fixed_point=str(cfg["ph.fixed_point"]), relaxation=float(cfg["ph.fixed_point.relaxation"]),
             tolerance=float(cfg["ph.fixed_point.tolerance"]), max_sweeps=int(cfg["ph.fixed_point.max_sweeps"]),
             flag=float(cfg["ph.fixed_point.flag"]), profile_node_step_k=float(cfg["sulfur.profile_node_step_k"]),
             ph_nodes=int(defaults["ph_nodes"].default), ph_iterations=int(defaults["ph_iterations"].default))
    for value, allowed in ((s["coupling"], COUPLINGS), (s["node_grid"], NODE_GRIDS), (s["fixed_point"], FIXED_POINTS)):
        if value not in allowed:
            raise UnsupportedConfiguration(f"{value!r}: not one of {allowed}")
    return s, cfg.digest()


def ph_model(cfg=None) -> PhModel:
    """The network's PhModel for this lane: the Config's numerics and ``ph_state``'s solvent tables."""
    from src.kinetic_core import ph_state as ph

    numerics, _ = config_values(cfg)
    return PhModel(**numerics, r_kj=float(ph.R_KJ), pka_ref_t_k=float(ph.PKA_REFERENCE_T_K),
                   label_t_k=float(ph.LABEL_TEMPERATURE_K),
                   pkw_table_k=tuple(c + CELSIUS for c, _ in ph.PKW_TABLE_C),
                   pkw_table=tuple(p for _, p in ph.PKW_TABLE_C),
                   rho_table_c=tuple(c for c, _ in ph._WATER_DENSITY_C),
                   rho_table=tuple(r for _, r in ph._WATER_DENSITY_C),
                   basis="ph_state (R, PKA_REFERENCE_T_K, LABEL_TEMPERATURE_K, PKW_TABLE_C, _WATER_DENSITY_C); the "
                         "Config's ph.* keys; integrate_sulfur's ph_nodes / ph_iterations defaults")


# --------------------------------------------------------------------------- the titratable table
def titratables() -> Tuple[Titratable, ...]:
    """``ph_state.titratable_inventory`` at the shipped configuration as pools, in the order the engine sums them
    (the organic-acid mix, formic and acetic, cysteine, the ARP class, TTCA; then every buffer group a
    ``BufferSpec`` can carry, pot-supplied). The ARP class's amine pK is the ``arp_amine_pka`` coordinate; the lumped
    acid pool enters at ``acid_yield`` per mmol. Centres per unit of pool from ``ph_state.TITRATABLE_CENTRES``."""
    from src.kinetic_core import ph_state as ph

    pools: Dict[Any, Dict[str, Any]] = {}

    def pool(group, key=None, **extra):
        k = key or group
        if k not in pools:
            pools[k] = {"group": group, "weights": {}, "yield": {}, **extra}
        return pools[k]

    def add(group, species, weight=0.0, acid=0.0):
        p = pool(group)
        if weight:
            p["weights"][species] = p["weights"].get(species, 0.0) + float(weight)
        if acid:
            p["yield"][species] = p["yield"].get(species, 0.0) + float(acid)
        p.setdefault("species", []).append(species)

    for group, share in ph.ORGANIC_ACID_MIX:           # the lumped pool (x acid yield) and the carried carboxyls
        add(group, "ACID", acid=share)
        add(group, "CBX", weight=share)
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
    out: List[Titratable] = []
    for p in pools.values():
        g = p["group"]
        centres = [dict(ph.TITRATABLE_CENTRES.get(s, {})) for s in p["species"]]
        if any(c != centres[0] for c in centres):
            raise UnsupportedConfiguration(f"pool {g.key}: its species carry different centres {centres}")
        k = len(g.pkas)
        out.append(Titratable(
            species=g.key, pka_25c=tuple(float(x) for x in g.pkas),
            charges=tuple(float(g.z_max) - i for i in range(k + 1)), basis=f"ph_state {g.key}: {g.source}"[:200],
            dh_kj_mol=tuple(float(x) for x in g.delta_h_kj_mol),
            dcp_j_mol_k=tuple(float(x) for x in (g.delta_cp_j_mol_k or (0.0,) * k)),
            weights=dict(p["weights"]), yield_weights=dict(p["yield"]) or None,
            yield_coordinate="acid_yield" if p["yield"] else None,
            pka_coordinates=tuple("arp_amine_pka" if (g.key == template.key and i == slots[0]) else None
                                  for i in range(k)),
            centres={kk: float(v) for kk, v in centres[0].items()}))
    seen = set()
    for kind in [k for k in ph.BufferSpec._KINDS if k not in ("none", "clamped")]:
        spec = ph.BufferSpec(kind=kind, phosphate_mol_l=1.0, citrate_mol_l=1.0, acetate_mol_l=1.0,
                             pyrophosphate_mol_l=1.0, phthalate_mol_l=1.0)
        for g, _ in spec.groups():
            if g.key in seen:
                continue
            seen.add(g.key)
            k = len(g.pkas)
            out.append(Titratable(
                species=f"buffer_{g.key}", pka_25c=tuple(float(x) for x in g.pkas),
                charges=tuple(float(g.z_max) - i for i in range(k + 1)), basis=f"buffer {g.key}: {g.source}"[:200],
                dh_kj_mol=tuple(float(x) for x in g.delta_h_kj_mol),
                dcp_j_mol_k=tuple(float(x) for x in (g.delta_cp_j_mol_k or (0.0,) * k)),
                source="pot", pot_key=g.key, centres={}))
    return tuple(out)


# --------------------------------------------------------------------------- the rate constants
def _rate_static(derived: Mapping[str, Any], parameters: Mapping[str, Any], reactions: Sequence[str]) -> Dict[str, Any]:
    """For a block of reactions: per generic parameter key (first appearance in the block): k_ref, Ea (None = held),
    1 / T_ref, its cap; per reaction: its parameter's index (-1 for a hand-ported step) and its pH-factor kind."""
    from src.kinetic_core import parameters as trunk_parameters
    from src.kinetic_core.parameters import KineticParameter
    from src.kinetic_core.parameters_furanic import TEMPERATURE_CAP_K

    keys = tuple(dict.fromkeys(derived["parameter_key"][r] for r in reactions if r not in structure.HAND_PORTED))
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
                     for r in reactions])
    kind = np.array([KIND_CODES[derived["kind"][r]] for r in reactions])
    return {"keys": keys, "k_ref": k_ref, "ea": ea, "inv_t_ref": np.asarray(inv_t_ref),
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


def _make_rate_constants(static, c, solvent, reactions: Sequence[str]):
    """The rate function of a block of reactions (the whole lane, or one of the composed modules): the generic rule on
    the block's own coordinates, then the hand-ported terms the block holds (each reading only its block's
    constants, ``TERM_READS``), then its held-at-zero steps."""
    keys, linear = static["keys"], static["linear"]
    has_ea = jnp.asarray(static["has_ea"])
    inv_t_ref = jnp.asarray(static["inv_t_ref"])
    cap = jnp.asarray(static["cap"])
    pidx = static["pidx"]
    safe = jnp.asarray(np.where(pidx >= 0, pidx, 0))
    generic = jnp.asarray(pidx >= 0)
    kind = jnp.asarray(static["kind"])
    ix = {r: j for j, r in enumerate(reactions)}
    terms = sorted(t for t in TERM_READS if t in ix)
    for t in terms:
        stray = [r for r in TERM_READS[t] if r not in ix]
        if stray:
            raise UnsupportedConfiguration(f"{t} reads {stray}, which another module owns: keep them in one module")
    unknown = sorted(r for r in reactions if r in structure.HAND_PORTED and r not in TERM_READS
                     and r not in structure.OFF)
    if unknown:
        raise UnsupportedConfiguration(f"hand-ported steps with no JAX term: {unknown}")
    off = [ix[r] for r in structure.OFF if r in ix]
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
        a = davies_a(t_k, solvent)
        k = jnp.where(generic, k_par[safe], 0.0) * ph_factors(ph, t_k, ionic, c, a)[kind]
        def k_adduct():
            return c["adduct_k_ref"] * jnp.exp(-(c["adduct_dh"] / r_kj) * (1.0 / t_k - 1.0 / c["adduct_t_k"]))

        def thiamine():
            share = thiamine_h2s_share(label, c)
            return jnp.where(share > 0.0, share / (1.0 - share), 0.0) * (k[ix["r_thi_hmp"]] + k[ix["r_thi_mesh"]])

        # every hand-ported term, evaluated only where this block holds the reaction (a term reads its own block's
        # constants, TERM_READS); tests/unit/test_jax_sulfur_structure.py reads this table's keys from the source
        derived = {
            "r_cys_h2s": lambda: cysteine_thermolysis_k(label, t_k, c),
            "r_gsh_h2s": lambda: glutathione_thermolysis_k(label, t_k, c),
            "r_thi_h2s": thiamine,
            "ch_thioether_release": lambda: k[ix["ch_thioether_mft"]] / k_adduct(),
            "ch_thioether_release_fft": lambda: k[ix["ch_thioether_fft"]] / k_adduct(),
            "r_amadori": lambda: c["amadori_factor"] * k[ix["r_schiff"]],
        }
        if terms:
            k = k.at[jnp.asarray([ix[n] for n in terms])].set(jnp.stack([derived[n]() for n in terms]))
        return k.at[jnp.asarray(off, dtype=jnp.int32)].set(0.0) if off else k

    return rate_constants


# --------------------------------------------------------------------------- the modules
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
    if "acid_yield" in module_params:
        bands["acid_yield"] = tuple(float(x) for x in ph.ACID_YIELD_BOUNDS)
        bands["arp_amine_pka"] = tuple(float(x) for x in ph.ARP_AMINE_PKA_BOUNDS)
        basis["acid_yield"] = ("engine.core_ph_drift() acid yield per sink event (B2.2 frozen); band "
                               "ph_state.ACID_YIELD_BOUNDS")
        basis["arp_amine_pka"] = ("engine.core_ph_drift() ARP secondary-ammonium pK (B2.2 frozen); band "
                                  "ph_state.ARP_AMINE_PKA_BOUNDS")
    return bands, basis


def _block_module(name: str, reactions: Sequence[str], species: Sequence[str], requires: Tuple[str, ...], *,
                  derived, parameters, c, model: PhModel, digest: str, charge: bool, competes: Tuple[str, ...] = ()
                  ) -> Module:
    """One module of the lane: ``reactions`` (in FULL_REACTIONS order), the ``species`` it owns; ``charge`` makes it
    the module that declares the lane's titratable table, the drift's coordinates and the PhModel."""
    from src.kinetic_core import engine, sulfur

    static = _rate_static(derived, parameters, reactions)
    params: Dict[str, float] = {}
    for key, lin, e in zip(static["keys"], static["linear"], static["ea"]):
        if lin:
            params[f"k_ref:{key}"] = 0.0
        else:
            params[f"log10_k:{key}"] = math.log10(float(parameters[key].k_ref))
        if e is not None:
            params[f"ea:{key}"] = float(e)
    if charge:
        drift = engine.core_ph_drift()
        params["acid_yield"] = float(drift.acid_yield)
        params["arp_amine_pka"] = float(drift.arp_amine_pka)
    bands, basis = _bands(params, parameters)
    hand = sorted(r for r in reactions if r in structure.HAND_PORTED)
    pools = titratables() if charge else ()
    ledger = {r: dict(e) for r, e in sulfur.CENTRE_LEDGER.items() if r in reactions}
    rset = set(reactions)
    return Module(
        name=name, species=tuple(species), elements={k: derived["elements"][k] for k in species},
        reactions=tuple(reactions), stoichiometry={r: derived["stoichiometry"][r] for r in reactions},
        reactants={r: derived["reactants"][r] for r in reactions}, params=params, bands=bands, basis=basis,
        rate_constants=_make_rate_constants(static, c, model.const, reactions),
        observables=observables_for(species), requires=requires, needs_ph=True, config_digest=digest,
        notes=(f"pH coupling {model.coupling}; node grid {model.node_grid}; fixed point {model.fixed_point} "
               f"(relaxation {model.relaxation}, tolerance {model.tolerance}, max sweeps {model.max_sweeps})",
               f"hand-ported reactions ({len(hand)}): {hand}",
               f"held at zero (route off at the shipped configuration): {[r for r in structure.OFF if r in rset]}",
               f"linear (zero k_ref) parameters: {[k for k, lin in zip(static['keys'], static['linear']) if lin]}",
               f"titratable pools: {[t.species for t in pools]}"),
        breakpoints=({"ph_label": (float(c["zheng_boundary"]), float(c["zheng_ph"][2]), float(c["zheng_ph"][3]), 7.0,
                                   9.0, float(c["thi_share_floor"]), float(c["thi_share_points"][1][0]))}
                     if rset & {"r_cys_h2s", "r_gsh_h2s", "r_thi_h2s"} else {}),
        competes=competes, titratables=pools, ph_model=model if charge else None, centre_ledger=ledger)


@dataclass(frozen=True)
class SulfurLane:
    """The sulfur lane's module(s), their Network and the PhModel (``settings``, PR #83's name)."""

    modules: Tuple[Module, ...]
    network: Network
    settings: PhModel
    parameters: Mapping[str, Any]

    @property
    def module(self) -> Module:
        """The one module of the ``build`` form."""
        if len(self.modules) != 1:
            raise ValueError(f"the composed lane has {len(self.modules)} modules: use .modules")
        return self.modules[0]

    @property
    def balance(self):
        """The merged titratable table (``ph_coupling.ChargeTable``)."""
        return self.network.charge

    @property
    def species(self):
        return self.network.species

    @property
    def index(self):
        return self.network.index

    def theta0(self):
        return self.network.theta0()


def _inputs(cfg):
    from src.kinetic_core import engine

    model = ph_model(cfg)
    _, digest = config_values(cfg)
    derived = structure.derive()
    problems = [p for p in structure.structure_problems(numeric=False)]
    if problems:
        raise UnsupportedConfiguration("the NumPy sulfur lane moved under the port:\n" + "\n".join(problems))
    return model, digest, derived, dict(engine.core_parameters(engine.SULFUR)), _term_constants()


def build_module(cfg=None) -> Tuple[Module, PhModel, Dict[str, Any]]:
    """The ONE-module sulfur lane for an engine Config (shipped values only on the held keys)."""
    model, digest, derived, parameters, c = _inputs(cfg)
    module = _block_module(NAME, derived["reactions"], derived["species"], (), derived=derived, parameters=parameters,
                           c=c, model=model, digest=digest, charge=True)
    return module, model, parameters


def build(cfg=None, *, solver: str = "Kvaerno5") -> SulfurLane:
    module, model, parameters = build_module(cfg)
    return SulfurLane(modules=(module,), network=compose([module], solver=solver), settings=model,
                      parameters=parameters)


def composed_blocks(derived) -> List[Tuple[str, Tuple[str, ...], Tuple[str, ...]]]:
    """(module, reactions, owned species) of the composed form: the engine's own blocks of FULL_REACTIONS
    (``network.REACTIONS``, ``network.FED_3DEOXY_REACTIONS``, the sulfur block), each species owned by the first block
    whose reactions touch it (an untouched one by the sulfur module)."""
    from src.kinetic_core import network as engine_network

    trunk = tuple(r.key for r in engine_network.REACTIONS)
    fed = tuple(r.key for r in engine_network.FED_3DEOXY_REACTIONS)
    full = tuple(derived["reactions"])
    if full[:len(trunk)] != trunk or full[len(trunk):len(trunk) + len(fed)] != fed:
        raise UnsupportedConfiguration("FULL_REACTIONS is no longer REACTIONS + FED_3DEOXY_REACTIONS + the sulfur "
                                       "block: the composed form's blocks do not hold")
    blocks = [(COMPOSED[0], trunk), (COMPOSED[1], fed), (COMPOSED[2], full[len(trunk) + len(fed):])]
    owner: Dict[str, str] = {}
    for name, rs in blocks:
        for r in rs:
            for s in list(derived["reactants"][r]) + list(derived["stoichiometry"][r]):
                owner.setdefault(s, name)
    out = []
    for name, rs in blocks:
        own = tuple(s for s in derived["species"] if owner.get(s, COMPOSED[2]) == name)
        out.append((name, rs, own))
    return out


def build_composed(cfg=None, *, solver: str = "Kvaerno5") -> SulfurLane:
    """The same lane as three modules composed by ``network.compose`` (``competition="refuse"``; the species order
    pinned to the engine's ``SULFUR_STATE_KEYS``). ``sulfur`` declares the titratable table (it reads FA / AA, which
    ``sulfur_trunk`` owns) and its consumption of the trunk's species. Coordinates: each block's own keys (no
    parameter key is shared between blocks: checked), the drift's in ``sulfur``."""
    model, digest, derived, parameters, c = _inputs(cfg)
    blocks = composed_blocks(derived)
    keys = [set(_rate_static(derived, parameters, rs)["keys"]) for _, rs, _ in blocks]
    shared = [(blocks[i][0], blocks[j][0], sorted(keys[i] & keys[j])) for i in range(3) for j in range(i + 1, 3)
              if keys[i] & keys[j]]
    if shared:
        raise UnsupportedConfiguration(f"parameter keys shared between blocks: {shared}")
    modules: List[Module] = []
    owned: Dict[str, str] = {s: name for name, _, own in blocks for s in own}
    for i, (name, rs, own) in enumerate(blocks):
        foreign = tuple(dict.fromkeys(s for r in rs for s in derived["reactants"][r]
                                      if owned[s] != name and derived["stoichiometry"][r].get(s, 0.0) < 0.0))
        modules.append(_block_module(name, rs, own, tuple(COMPOSED[:i]), derived=derived, parameters=parameters, c=c,
                                     model=model, digest=digest, charge=(name == NAME), competes=foreign))
    net = compose(modules, layout=derived["species"], solver=solver, competition="refuse")
    return SulfurLane(modules=tuple(modules), network=net, settings=model, parameters=parameters)


def theta_from(lane: SulfurLane, parameters: Optional[Mapping[str, Any]] = None, drift=None) -> Dict[str, Any]:
    """theta for an engine parameter dict (``k_ref``, ``ea_kj_mol``) and a PhDrift, the rest at shipped values."""
    out: Dict[str, Any] = {}
    for m in lane.modules:
        th = dict(m.params)
        for name in th:
            kind, _, key = name.partition(":")
            if parameters is None or key not in parameters:
                continue
            p = parameters[key]
            if kind == "log10_k":
                if not float(p.k_ref) > 0.0:
                    raise UnsupportedConfiguration(f"{key}: k_ref {p.k_ref!r} on a log10 coordinate")
                th[name] = math.log10(float(p.k_ref))
            elif kind == "k_ref":
                th[name] = float(p.k_ref)
            elif kind == "ea":
                th[name] = float(p.ea_kj_mol)
        if drift is not None and "acid_yield" in th:
            th["acid_yield"], th["arp_amine_pka"] = float(drift.acid_yield), float(drift.arp_amine_pka)
        out[m.name] = {k: jnp.asarray(float(v)) for k, v in th.items()}
    return out


# --------------------------------------------------------------------------- the pot
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
    charged species: ``ph_coupling.pot_from_program``'s arrays and ``mode`` (static), plus ``dropped`` and
    ``buffer_spec`` (what the engine's path needs to reproduce it)."""
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
    if mode == "clamped" and process.ph_final is not None and float(process.ph_final) != float(process.ph):
        raise UnsupportedConfiguration("a prescribed pH line between two measured endpoints (B2.1) is not ported")
    if thermal.profile is not None:
        pot = ph_coupling.pot_from_program(lane.network, y0, float(process.ph), lane.balance.buffer_vector(spec),
                                           profile=[(float(t), float(c)) for t, c in thermal.profile], mode=mode,
                                           duration=float(thermal.total_minutes))
    else:
        pot = ph_coupling.pot_from_program(lane.network, y0, float(process.ph), lane.balance.buffer_vector(spec),
                                           segments=[(float(d), float(t)) for d, t in thermal.segments], mode=mode)
    pot.update({"dropped": dropped, "buffer_spec": spec})
    return pot


pot_arrays = ph_coupling.pot_arrays


# --------------------------------------------------------------------------- integration
def integrate(lane: SulfurLane, theta, pot: Mapping[str, Any], *, rtol: float = 1e-8, atol: float = 1e-14,
              solver=None, adjoint=None, reverse: Optional[bool] = None, max_steps: int = 200_000,
              throw: bool = True, force_sweeps: Optional[int] = None) -> PhRun:
    """One pot through its programme: ``network.integrate`` on the lane's network (the network's own pH solve), as
    the node-resolved ``PhRun``. ``reverse=True`` keeps PR #83's scan over the sweeps; by default the sweep loop
    follows the adjoint (``ph_coupling.loop_kind``: ``jax.grad`` needs no flag)."""
    from . import network

    return network.integrate(lane.network, theta, pot, rtol=rtol, atol=atol, solver=solver, adjoint=adjoint,
                             max_steps=max_steps, throw=throw, reverse=reverse, force_sweeps=force_sweeps).ph


def compiled(lane: SulfurLane, pot_template: Mapping[str, Any], **kw):
    """``f(theta, arrays) -> PhRun``, jit-compiled once per lane, mode and pot shape."""
    mode = pot_template.get("mode", "dynamic")
    return jax.jit(lambda theta, arrays: integrate(lane, theta, dict(arrays, mode=mode), **kw))


def observe(lane: SulfurLane, run: PhRun, species: Sequence[str], rows: Sequence[int], *,
            floor: float = 1e-300):
    """log10 mmol/L of ``species`` at node ``rows`` of the LAST segment, then the cooled pH at the end."""
    ys = run.ys[-1]
    idx = [lane.index[s] for s in species]
    vals = [jnp.log10(jnp.maximum(ys[r, i], floor)) for i in idx for r in rows]
    vals.append(run.ph_cool[-1, -1])
    return jnp.stack(vals)


def jacobian(lane: SulfurLane, pot_template, coordinates: Sequence[str], species, rows, **kw):
    """``f(theta, arrays) -> (observables, d observables / d coordinates)``, exact, forward mode (coordinates of the
    one-module form, ``theta[NAME]``)."""
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
