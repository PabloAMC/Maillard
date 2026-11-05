"""
src/kinetic_core_jax/trunk.py -- the trunk as a module under the A2 contract v2 (2026-11-05): ``network.REACTIONS``,
the 26 steps the ACRYLAMIDE lane runs (B1's fifteen Martins steps and B7's eleven furanic steps), with their rate laws.

WHICH TRUNK (contract v1 sec. 1 found three, all built on this one):

* ``network.REACTIONS`` (26): the core every lane shares. The acrylamide lane runs exactly it
  (``acrylamide.FULL_ACRYLAMIDE_REACTIONS = REACTIONS + ACRYLAMIDE_REACTIONS``). THIS module.
* ``network.TRUNK_REACTIONS`` (55): the trunk LANE, ``REACTIONS`` plus seven blocks -- dicarbonyl (5), pyrazine (5),
  glycation (5), aqueous glyoxal (1), methionine (6), proline (4) and the fed 3-deoxy triangle (3). Those blocks are
  later modules that ``require`` this one; they are not ported here.
* sulfur's ``REACTIONS + FED_3DEOXY_REACTIONS + SULFUR_REACTIONS`` (132): this core, the fed 3-deoxy block and the
  sulfur steps; its copy of the trunk keeps the sulfur lane's own pH machinery. Ported with the sulfur module.

RATE LAW (read from the engine, never retyped): k_j = k_ref,j exp(-(Ea_j / R)(1/T - 1/T_ref)) with T_ref = 100 C
(``parameters.T_REF_K``) and R = ``parameters.R_KJ``; ``r_amadori`` is DERIVED, ratio x k(r_schiff) x loading
(``SCHIFF_AMADORI_SPLIT``), as both engine paths compute it. Values from ``engine.core_parameters(lane)``.

TWO ENGINE CONTEXTS, chosen by ``build(..., lane=)``:

* ``lane="acrylamide"`` (default): the trunk as the acrylamide lane runs it -- ``acrylamide_rate_constants_at``, NO
  condition terms (that lane applies ``acrylamide_conditions`` to its own steps only). pH mode "none".
* ``lane="trunk"``: the trunk LANE's declared condition terms (``trunk_conditions.apply``, B12 / B40-41) on the
  steps of this module: the a_w multiplier on ``k_schiff`` (Pereyra Gonzales 2010, scaled by ``aw_scale``), the
  pH factor 10^(e (pH - 6.8)) on the three Amadori-decay steps, and -- while ``trunk.three_deoxy_exit_ph_term`` is on
  -- each ``trunk.three_deoxy_exit_ph`` exit term (k_tdg_fa shipped). They scale k_ref at the pot's DECLARED pH and
  a_w (no dynamic pH): pH mode "prescribed", with the declared pH constant over the programme. Accepted at the rate
  level only: the trunk LANE's trajectories also run the seven blocks above, which this module does not carry. The
  pyrazine / methionine / proline pH terms and the atmosphere factor act on steps outside this module; a non-air
  atmosphere is refused at the pot (the engine raises there too).

Free coordinates: ``log10_k:<parameter>`` (log10 of k_ref at T_ref) and ``ea:<parameter>`` (kJ/mol) for each of the
25 parameters; with ``lane="trunk"`` also ``aw_scale``, ``ph_exponent`` and ``ph_exponent_exit:<parameter>``.
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import jax.numpy as jnp
import numpy as np

from .lane import UnsupportedConfiguration
from .network import Module, RateCtx

NAME = "trunk"
LANES = ("acrylamide", "trunk")


def _engine():
    from src.kinetic_core import config, engine, network, parameters, trunk_conditions

    return config, engine, network, parameters, trunk_conditions


# --------------------------------------------------------------------------- a block of Arrhenius steps
class ArrheniusBlock:
    """The static layout of a list of engine ``Reaction``s whose constants are ``parameter.k_at(T)``: per reaction
    the index of its parameter (or -1 for a derived step), per parameter its T_ref. ``k(theta, factors, t_k)`` returns
    the parameters' constants, ``k_ref x factor x exp(-(Ea / R)(1/T - 1/T_ref))``, in parameter order."""

    def __init__(self, reactions, parameters: Mapping[str, Any], t_ref_of, r_kj: float):
        self.reactions = tuple(reactions)
        self.keys: Tuple[str, ...] = tuple(dict.fromkeys(r.parameter_key for r in self.reactions
                                                         if r.parameter_key is not None))
        for key in self.keys:
            p = parameters.get(key)
            if p is None:
                raise KeyError(f"no engine parameter {key!r}")
            if p.k_ref is None or p.ea_kj_mol is None:
                raise UnsupportedConfiguration(f"{key}: unpopulated in the engine (it refuses to integrate too)")
            if not float(p.k_ref) > 0.0:
                raise UnsupportedConfiguration(f"{key}: k_ref {p.k_ref!r} has no log10 coordinate")
        self.k_ref = {k: float(parameters[k].k_ref) for k in self.keys}
        self.ea = {k: float(parameters[k].ea_kj_mol) for k in self.keys}
        self.parameters = {k: parameters[k] for k in self.keys}
        self.inv_t_ref = jnp.asarray([1.0 / float(t_ref_of(parameters[k])) for k in self.keys])
        self.r_kj = float(r_kj)
        self.pidx = np.array([self.keys.index(r.parameter_key) if r.parameter_key is not None else -1
                              for r in self.reactions])
        self.rxn_of = {r.key: j for j, r in enumerate(self.reactions)}

    def coordinates(self) -> Dict[str, float]:
        out: Dict[str, float] = {}
        for k in self.keys:
            out[f"log10_k:{k}"] = math.log10(self.k_ref[k])
            out[f"ea:{k}"] = self.ea[k]
        return out

    def k_params(self, theta: Mapping[str, Any], factor, t_k):
        log10_k = jnp.stack([jnp.asarray(theta[f"log10_k:{k}"], dtype=jnp.float64) for k in self.keys])
        ea = jnp.stack([jnp.asarray(theta[f"ea:{k}"], dtype=jnp.float64) for k in self.keys])
        return (10.0 ** log10_k) * factor * jnp.exp(-(ea / self.r_kj) * (1.0 / t_k - self.inv_t_ref))

    def k_reactions(self, k_params):
        """Per reaction, its parameter's constant (0 for a derived step, filled in by the caller)."""
        safe = np.where(self.pidx >= 0, self.pidx, 0)
        return jnp.where(jnp.asarray(self.pidx >= 0), k_params[jnp.asarray(safe)], 0.0)

    def stoichiometry(self) -> Tuple[Dict[str, Dict[str, float]], Dict[str, Dict[str, int]]]:
        stoich, reactants = {}, {}
        for r in self.reactions:
            net: Dict[str, float] = {}
            for s, c in r.reactants.items():
                net[s] = net.get(s, 0.0) - float(c)
            for s, c in r.products.items():
                net[s] = net.get(s, 0.0) + float(c)
            stoich[r.key] = {s: c for s, c in net.items() if c != 0.0}
            reactants[r.key] = {s: int(c) for s, c in r.reactants.items()}
        return stoich, reactants

    def bands_and_basis(self, ea_band_override: Optional[Mapping[str, Tuple[float, float]]] = None,
                        log10_k_bounds: Optional[Mapping[str, Tuple[float, float]]] = None):
        bands: Dict[str, Optional[Tuple[float, float]]] = {}
        basis: Dict[str, str] = {}
        for k in self.keys:
            p = self.parameters[k]
            src = f"{p.evidence_class}; {str(p.source_anchor)[:90]}"
            bounds = (log10_k_bounds or {}).get(k)
            bands[f"log10_k:{k}"] = tuple(bounds) if bounds is not None else None
            basis[f"log10_k:{k}"] = f"log10 k_ref of engine parameter {k} ({src})" + (
                "; band = the fit's search bounds" if bounds is not None else "")
            ci = getattr(p, "ea_ci95_kj_mol", None)
            if ea_band_override and k in ea_band_override:
                bands[f"ea:{k}"] = tuple(ea_band_override[k])
                basis[f"ea:{k}"] = f"Ea of {k}; band declared by the engine's envelope ({src})"
            elif ci is not None:
                bands[f"ea:{k}"] = (self.ea[k] - float(ci), self.ea[k] + float(ci))
                basis[f"ea:{k}"] = f"Ea of {k}; band = its printed 95 % CI ({src})"
            else:
                bands[f"ea:{k}"] = None
                basis[f"ea:{k}"] = f"Ea of {k}; no band declared ({src})"
        return bands, basis


def element_table(keys: Sequence[str]) -> Dict[str, Tuple[int, int, int]]:
    from src.kinetic_core.species_acrylamide import ACRYLAMIDE_BY_KEY

    return {k: (int(ACRYLAMIDE_BY_KEY[k].carbon), int(ACRYLAMIDE_BY_KEY[k].nitrogen), int(ACRYLAMIDE_BY_KEY[k].sulfur))
            for k in keys}


def species_of(reactions) -> Tuple[str, ...]:
    """The species a reaction list touches, in the engine's acrylamide state order (a trunk state is its prefix)."""
    from src.kinetic_core.species_acrylamide import ACRYLAMIDE_STATE_KEYS

    used = {s for r in reactions for s in list(r.reactants) + list(r.products)}
    return tuple(k for k in ACRYLAMIDE_STATE_KEYS if k in used)


def observables_for(keys: Sequence[str]) -> Dict[str, Any]:
    """``ug_per_L:<key>`` for every key with a molecular weight (``species_sulfur.mmol_per_litre_to_ug_per_litre``,
    the engine's reporting unit), and ``mmol_per_L:<key>`` for every key."""
    from src.kinetic_core.species_sulfur import MOLECULAR_WEIGHT_G_PER_MOL, mmol_per_litre_to_ug_per_litre

    out: Dict[str, Any] = {}
    for key in keys:
        out[f"mmol_per_L:{key}"] = (lambda y, index, key=key: y[index[key]])
        if key in MOLECULAR_WEIGHT_G_PER_MOL:
            per = float(mmol_per_litre_to_ug_per_litre(key, 1.0))
            out[f"ug_per_L:{key}"] = (lambda y, index, key=key, per=per: y[index[key]] * per)
    return out


def nan_or(value: Optional[float]) -> float:
    """A pot scalar: the engine's None as NaN."""
    return float("nan") if value is None else float(value)


# --------------------------------------------------------------------------- the module
def config_values(cfg=None, lane: str = "acrylamide") -> Dict[str, Any]:
    """The keys this module reads from an engine Config (default: the process's). On the acrylamide lane the trunk
    reads none (``trunk_conditions`` is not applied there); on the trunk lane, the two 3-DG exit keys."""
    config, *_ = _engine()
    cfg = config.Config.from_globals() if cfg is None else cfg
    out: Dict[str, Any] = {"digest": cfg.digest()}
    if lane == "trunk":
        out["exit_term"] = bool(cfg["trunk.three_deoxy_exit_ph_term"])
        out["exit_ph"] = dict(cfg["trunk.three_deoxy_exit_ph"])
    return out


def build(cfg=None, *, lane: str = "acrylamide") -> Module:
    """The trunk Module (``network.REACTIONS``) for an engine Config, in the context of ``lane``."""
    if lane not in LANES:
        raise ValueError(f"lane must be one of {LANES}")
    _, engine, network, parameters_mod, tc = _engine()
    values = config_values(cfg, lane)
    reactions = tuple(network.REACTIONS)
    for r in reactions:
        if r.statistical_factor != 1.0:
            # acrylamide_rate_constants_at ignores the factor and network.rate_constants_at applies it: the two engine
            # paths would disagree, and the port would have to pick one in silence
            raise UnsupportedConfiguration(f"{r.key}: a statistical factor {r.statistical_factor} on a trunk step")
    parameters = engine.core_parameters(engine.ACRYLAMIDE if lane == "acrylamide" else engine.TRUNK)
    t_ref = float(parameters_mod.T_REF_K)
    block = ArrheniusBlock(reactions, parameters, lambda p: t_ref, parameters_mod.R_KJ)
    split = parameters_mod.SCHIFF_AMADORI_SPLIT
    amadori_per_schiff = (float(split["ratio_amadori_over_schiff_pseudo_first_order"]),
                          float(split["amine_loading_mmol_L_for_the_ratio"]))
    j_amadori, j_schiff = block.rxn_of["r_amadori"], block.rxn_of["r_schiff"]
    derived = [r.key for r in reactions if r.parameter_key is None]
    if derived != ["r_amadori"]:
        raise UnsupportedConfiguration(f"derived trunk steps {derived}: only r_amadori is ported")

    params = block.coordinates()
    bands, basis = block.bands_and_basis()
    notes: List[str] = [f"lane context: {lane}"]
    pot_inputs: Tuple[str, ...] = ()
    n_par = len(block.keys)
    if lane == "trunk":
        exits = {k: v for k, v in values["exit_ph"].items()} if values["exit_term"] else {}
        outside = sorted(set(exits) - set(block.keys))
        if outside:
            raise UnsupportedConfiguration(
                f"trunk.three_deoxy_exit_ph names {outside}, steps outside network.REACTIONS (a trunk-lane block "
                "this port does not carry)")
        aw_table = np.asarray([a for a, _ in tc.AW_MULTIPLIER_TABLE])
        aw_mult = np.asarray([m for _, m in tc.AW_MULTIPLIER_TABLE])
        aw_rows = np.asarray([k in tc.AW_STEPS for k in block.keys])
        ph_rows = np.asarray([k in tc.PH_STEPS for k in block.keys])
        if aw_rows.sum() != len(tc.AW_STEPS) or ph_rows.sum() != len(tc.PH_STEPS):
            raise UnsupportedConfiguration("a trunk_conditions step is not a step of network.REACTIONS")
        exit_keys = tuple(exits)
        exit_rows = np.asarray([[k == e for k in block.keys] for e in exit_keys], dtype=float).reshape(
            len(exit_keys), n_par)
        params.update({"aw_scale": 1.0, "ph_exponent": float(tc.PH_EXPONENT_DECADES_PER_UNIT)})
        bands.update({"aw_scale": tuple(tc.AW_SCALE_BAND), "ph_exponent": tuple(tc.PH_EXPONENT_BAND)})
        basis.update({"aw_scale": f"trunk_conditions.AW_SCALE_BAND (scale on m - 1; {tc.AW_SOURCE[:70]})",
                      "ph_exponent": f"trunk_conditions.PH_EXPONENT_BAND ({tc.PH_SOURCE[:70]})"})
        for e, (exponent, band, source) in exits.items():
            params[f"ph_exponent_exit:{e}"] = float(exponent)
            bands[f"ph_exponent_exit:{e}"] = tuple(band)
            basis[f"ph_exponent_exit:{e}"] = f"trunk.three_deoxy_exit_ph ({source[:70]})"
        ref_ph, ref_aw, floor = float(tc.REFERENCE_PH), float(tc.REFERENCE_AW), float(tc.AW_TABLE_FLOOR)
        pot_inputs = ("ph_label", "aw")
        notes += [f"trunk_conditions: a_w on {tc.AW_STEPS}, pH on {tc.PH_STEPS}, exit terms on {exit_keys}",
                  "pH mode prescribed (the declared pH, constant over the programme)"]

        def factor(theta, ctx: RateCtx):
            aw_ok = ~jnp.isnan(ctx.aw)
            aw = jnp.where(aw_ok, ctx.aw, ref_aw)
            m = jnp.where(aw >= ref_aw, 1.0, jnp.where(aw <= floor, aw_mult[0], jnp.interp(aw, aw_table, aw_mult)))
            m_eff = jnp.where(aw_ok, 1.0 + (m - 1.0) * theta["aw_scale"], 1.0)
            ph_ok = ~jnp.isnan(ctx.ph_label)
            dph = jnp.where(ph_ok, ctx.ph_label - ref_ph, 0.0)
            f = 10.0 ** (theta["ph_exponent"] * dph)
            applied = jnp.abs(f - 1.0) > 1e-12      # trunk_conditions.apply scales (and adds the exits) only then
            out = jnp.where(jnp.asarray(aw_rows), m_eff, 1.0) * jnp.where(jnp.asarray(ph_rows) & applied, f, 1.0)
            for i, e in enumerate(exit_keys):
                fx = 10.0 ** (theta[f"ph_exponent_exit:{e}"] * dph)
                out = out * jnp.where((jnp.asarray(exit_rows[i]) > 0) & applied, fx, 1.0)
            return out
    else:
        def factor(theta, ctx: RateCtx):
            return jnp.ones(n_par)

    def rate_constants(theta: Mapping[str, Any], ctx: RateCtx):
        k = block.k_reactions(block.k_params(theta, factor(theta, ctx), ctx.t_k))
        return k.at[j_amadori].set(amadori_per_schiff[0] * k[j_schiff] * amadori_per_schiff[1])

    stoich, reactants = block.stoichiometry()
    species = species_of(reactions)
    return Module(
        name=NAME, species=species, elements=element_table(species), reactions=tuple(r.key for r in reactions),
        stoichiometry=stoich, reactants=reactants, params=params, bands=bands, basis=basis,
        rate_constants=rate_constants, observables=observables_for(species), requires=(), needs_ph=False,
        config_digest=values["digest"], notes=tuple(notes), pot_inputs=pot_inputs)


def theta_from_parameters(module: Module, parameters: Mapping[str, Any]) -> Dict[str, float]:
    """A module's coordinates read off an engine operative parameter dict (``k_ref``, ``ea_kj_mol``), the others at
    their shipped values: how a reference's perturbed parameter set enters the port."""
    out = dict(module.params)
    for name in module.params:
        kind, _, key = name.partition(":")
        if kind == "log10_k" and key in parameters:
            out[name] = math.log10(float(parameters[key].k_ref))
        elif kind == "ea" and key in parameters:
            out[name] = float(parameters[key].ea_kj_mol)
    return out


def check_atmosphere(atmosphere: Optional[str]) -> None:
    """The trunk lane's B29 axis: only air (or no declaration) is answered by the engine without a fitted factor; the
    factor scales steps outside this module, so any other atmosphere is refused here."""
    if atmosphere is not None and str(atmosphere) != "air":
        raise UnsupportedConfiguration(f"atmosphere {atmosphere!r}: the oxidative entries are outside this module")
