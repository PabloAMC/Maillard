"""
src/kinetic_core_jax/trunk_lane.py -- the trunk LANE's seven blocks under the A2 contract v2.1 (2026-11-05). Composed
after ``trunk.build(cfg, lane="trunk")`` they are the engine's trunk lane, ``network.TRUNK_REACTIONS`` (55 steps):

    compose([trunk.build(cfg, lane="trunk"), trunk_lane.build(cfg), trunk_lane.build_fed_3deoxy(cfg)])
    == REACTIONS (26) + DICARBONYL (5) + PYRAZINE (5) + GLYCATION (5) + AQUEOUS_GLYOXAL (1) + METHIONINE (6)
       + PROLINE (4) + FED_3DEOXY (3), in the engine's order        (``compose_lane``)

TWO MODULES, NOT SEVEN. ``trunk_lane`` carries the six TRUNK-ONLY blocks (dicarbonyl, pyrazine, glycation, aqueous
glyoxal, methionine, proline: 26 steps); ``fed_3deoxy`` carries the B39 triangle (3 steps, owns DGAL). The split
follows the contract's two rules on ownership and coordinates: (i) the six trunk-only blocks share one condition
coordinate pair -- B18's two pyrazine pH slopes scale the pyrazine, methionine AND proline Strecker steps
(``trunk_conditions.apply``) -- and a rate function reads only its own module's theta, so those blocks are one
module; (ii) the fed 3-deoxy triangle is ALSO in the sulfur network (``sulfur.FULL_REACTIONS = REACTIONS +
FED_3DEOXY_REACTIONS + ...``), so it is its own module, which the sulfur port can compose instead of redefining its
three reaction keys. Neither module owns a species another block also owns: GO, G, DA (dicarbonyl), AKG / AKM and the
pyrazines, LYSP and its products, MET and its chain, PRO and its chain belong to ``trunk_lane``; DGAL to
``fed_3deoxy``; Glc, Gly, MGO, TDG, DDG, FRAG_C, MEL_* to the trunk.

RATE LAW (read from the engine): ``KineticParameter.k_at`` (k_ref at T_ref = 100 C, ``parameters.T_REF_K``,
R = ``parameters.R_KJ``), times the step's ``statistical_factor`` (CF-077's 2 on the cross pair ``r_ak_mpz``), as
``network.rate_constants_at`` computes it. No derived step.

THE PARAMETER SET, FROM THE CONFIG WITHOUT INSTALLING IT (``lane_parameters``): ``engine.core_parameters("trunk")``
with the Config's values laid over the three blocks its keys decide, exactly as ``operative_parameters`` lays them
(the blocks are disjoint, which the build checks):

* ``dicarbonyl.aqueous_glyoxal_parameters`` (rate phase, a mapping): its ``k_ama_g`` / ``k_g_go`` replace the
  engine's, any value ported;
* ``methionine.shipped_b22`` / ``proline.shipped_b24`` (build phase): the block is the engine's own
  ``with_fitted_methionine(FROZEN_B22 if shipped else INERT_B22)`` (likewise ``with_fitted_proline``), so BOTH values
  are ported in any process (the engine itself would need a re-import to switch them);
* ``methionine.shipped_b22b`` / ``proline.shipped_b24b`` (build phase): REFUSED at True -- ``FROZEN_B22B`` /
  ``FROZEN_B24B`` are empty, so the engine cannot build that value either (KeyError at import);
* ``dicarbonyl.shipped_b39`` (build phase): ported at both values; it moves no number (``with_fed_3deoxy`` reads it
  for the flags, evidence class and source anchor only; ``uncertainty`` reads it for the envelope's draws).

Steps the engine carries at EXACTLY zero (the unshipped B22 / B22b / B24 / B24b steps, ``k_da_sink``) have a LINEAR
``k_ref:`` coordinate (``trunk.ArrheniusBlock(allow_zero=True)``); every other parameter ``log10_k:`` and ``ea:``.

CONDITION TERMS on these steps (``trunk_conditions.apply``, read from the pot's DECLARED pH; pH mode "prescribed"):

* B18's two-slope pyrazine pH term (knot 7, reference 6.8) on ``k_go_ak`` / ``k_mgo_ak`` and, B22 / B24, on the
  methionine (``k_go_met``, ``k_mgo_met``) and proline (``k_mgo_pro``, ``k_ha_athp``) Strecker steps, applied only
  where the factor differs from 1 by more than 1e-12 (as ``apply``); coordinates ``pyrazine_ph_slope_above_knot`` /
  ``pyrazine_ph_slope_below_knot`` at ``PYRAZINE_PH_SLOPES``;
* B29's atmosphere factor on ``k_glc_g`` / ``k_ama_g``: air (or none declared) is exactly 1; any other atmosphere is
  refused at the pot (``trunk.check_atmosphere``), as the engine refuses it while ``ATMOSPHERE_FACTORS`` holds air
  alone -- and the build refuses if that table ever holds more;
* no a_w term and no B12 / B40 pH term touches these steps (they act on the trunk module's steps); an exit term
  naming one of them is refused by ``trunk.build``.

THE POT (``make_pot``): methionine and proline are charged AS GLYCINE TOO (B22 / B24's declared convention in
``engine._integrate_program``: Gly += PRO, Gly += MET, then each keeps its own pool), the bound lysine LYSP is the
caller's (``lysp_charge`` computes the engine's, from the protein loading and the matrix's site density times the
available-fraction centre), and a non-air atmosphere or a charged species the network does not carry is refused.
"""
from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import jax.numpy as jnp
import numpy as np

from .lane import UnsupportedConfiguration
from .network import Module, RateCtx, compose, program_segments
from .trunk import ArrheniusBlock, check_atmosphere, nan_or, observables_for

NAME = "trunk_lane"
FED_NAME = "fed_3deoxy"
#: the six trunk-only blocks, in the engine's order (network.TRUNK_REACTIONS)
BLOCKS = ("DICARBONYL_REACTIONS", "PYRAZINE_REACTIONS", "GLYCATION_REACTIONS", "AQUEOUS_GLYOXAL_REACTIONS",
          "METHIONINE_REACTIONS", "PROLINE_REACTIONS")
SLOPE_COORDINATES = ("pyrazine_ph_slope_above_knot", "pyrazine_ph_slope_below_knot")


def _engine():
    from src.kinetic_core import config, engine, network, parameters, trunk_conditions

    return config, engine, network, parameters, trunk_conditions


def config_values(cfg=None) -> Dict[str, Any]:
    """The keys the trunk lane's blocks read from an engine Config (default: the process's)."""
    config, *_ = _engine()
    cfg = config.Config.from_globals() if cfg is None else cfg
    return {"digest": cfg.digest(),
            "aqueous_glyoxal": dict(cfg["dicarbonyl.aqueous_glyoxal_parameters"]),
            "shipped_b39": bool(cfg["dicarbonyl.shipped_b39"]),
            "methionine_b22": bool(cfg["methionine.shipped_b22"]),
            "methionine_b22b": bool(cfg["methionine.shipped_b22b"]),
            "proline_b24": bool(cfg["proline.shipped_b24"]),
            "proline_b24b": bool(cfg["proline.shipped_b24b"])}


def lane_parameters(cfg=None) -> Dict[str, Any]:
    """The trunk lane's operative parameter dict for a Config, built from the engine's own functions and never by
    installing the Config (module docstring)."""
    from src.kinetic_core import parameters_methionine as pm
    from src.kinetic_core import parameters_proline as pp

    _, engine, *_ = _engine()
    values = config_values(cfg)
    if values["methionine_b22b"]:
        raise UnsupportedConfiguration("methionine.shipped_b22b = True: FROZEN_B22B is empty, so the engine has no "
                                       "shipped B22b values (it would raise at import)")
    if values["proline_b24b"]:
        raise UnsupportedConfiguration("proline.shipped_b24b = True: FROZEN_B24B is empty, so the engine has no "
                                       "shipped B24b values (it would raise at import)")
    out = dict(engine.core_parameters(engine.TRUNK))
    met = pm.with_fitted_methionine(*[(pm.FROZEN_B22 if values["methionine_b22"] else pm.INERT_B22)[k]
                                      for k in pm.METHIONINE_COORDINATES])
    pro = pp.with_fitted_proline(*[(pp.FROZEN_B24 if values["proline_b24"] else pp.INERT_B24)[k]
                                   for k in pp.PROLINE_COORDINATES])
    aqueous = values["aqueous_glyoxal"]
    from src.kinetic_core.parameters_dicarbonyl import FED_3DEOXY_PARAMETERS

    blocks = (set(aqueous), set(met), set(pro), set(FED_3DEOXY_PARAMETERS))
    if sum(len(b) for b in blocks) != len(set().union(*blocks)):
        raise AssertionError("the Config-decided parameter blocks overlap; laying them over core_parameters would "
                             "depend on their order")
    missing = sorted(set().union(*blocks[:3]) - set(out))
    if missing:
        raise AssertionError(f"{missing}: not in engine.core_parameters('trunk')")
    out.update(aqueous)
    out.update(met)
    out.update(pro)
    return out


def _reactions(names: Sequence[str]):
    _, _, network, *_ = _engine()
    return tuple(r for name in names for r in getattr(network, name))


def _species(reactions, exclude: Sequence[str]) -> Tuple[str, ...]:
    from src.kinetic_core.species import SPECIES_KEYS

    used = {s for r in reactions for s in list(r.reactants) + list(r.products)}
    return tuple(k for k in SPECIES_KEYS if k in used and k not in set(exclude))


def _elements(keys: Sequence[str]) -> Dict[str, Tuple[int, int, int]]:
    from src.kinetic_core.species import BY_KEY

    return {k: (int(BY_KEY[k].carbon), int(BY_KEY[k].nitrogen), int(BY_KEY[k].sulfur)) for k in keys}


def _trunk_species() -> Tuple[str, ...]:
    from .trunk import species_of

    _, _, network, *_ = _engine()
    return species_of(network.REACTIONS)


def _check_exit_terms(cfg, keys: Sequence[str]) -> None:
    config, *_ = _engine()
    cfg = config.Config.from_globals() if cfg is None else cfg
    if bool(cfg["trunk.three_deoxy_exit_ph_term"]):
        named = sorted(set(dict(cfg["trunk.three_deoxy_exit_ph"])) & set(keys))
        if named:
            raise UnsupportedConfiguration(f"trunk.three_deoxy_exit_ph names {named}, steps of the trunk lane's "
                                           "blocks; the port carries the exit terms on network.REACTIONS only")


def _block_module(name: str, reactions, parameters, cfg, *, exclude: Sequence[str], factor_rows=None,
                  slopes=None, notes: Sequence[str] = (), digest: str = "") -> Module:
    _, _, _, parameters_mod, tc = _engine()
    for r in reactions:
        if r.parameter_key is None:
            raise UnsupportedConfiguration(f"{r.key}: a derived step in the trunk lane's blocks is not ported")
    t_ref = float(parameters_mod.T_REF_K)
    block = ArrheniusBlock(reactions, parameters, lambda p: t_ref, parameters_mod.R_KJ, allow_zero=True)
    _check_exit_terms(cfg, block.keys)
    sf = jnp.asarray([float(r.statistical_factor) for r in reactions])
    params = block.coordinates()
    bands, basis = block.bands_and_basis()
    pot_inputs: Tuple[str, ...] = ()
    breaks: Dict[str, Tuple[float, ...]] = {}
    if factor_rows is not None:
        from src.kinetic_core import parameters_pyrazine as pz

        knot, ref = float(pz.PYRAZINE_PH_KNOT), float(pz.PYRAZINE_REFERENCE_PH)
        breaks = {"ph_label": (knot,)}          # B18's two-slope term kinks at its knot
        rows = jnp.asarray(factor_rows)
        params.update(dict(zip(SLOPE_COORDINATES, (float(s) for s in slopes))))
        bands.update({c: None for c in SLOPE_COORDINATES})
        basis.update({c: f"parameters_pyrazine.PYRAZINE_PH_SLOPES (B18 fit; {pz.PYRAZINE_PH_SOURCE[:70]}); bands "
                         "= the B18 report's bounds, not read here" for c in SLOPE_COORDINATES})
        pot_inputs = ("ph_label",)

        def factor(theta, ctx: RateCtx):
            s_hi, s_lo = theta[SLOPE_COORDINATES[0]], theta[SLOPE_COORDINATES[1]]

            def shape(x):
                return jnp.where(x >= knot, s_hi * (x - knot), s_lo * (x - knot))

            ok = ~jnp.isnan(ctx.ph_label)
            ph = jnp.where(ok, ctx.ph_label, ref)
            fp = jnp.where(ok, 10.0 ** (shape(ph) - shape(jnp.asarray(ref))), 1.0)
            applied = jnp.abs(fp - 1.0) > 1e-12        # trunk_conditions.apply scales only then
            return jnp.where(rows & applied, fp, 1.0)
    else:
        n = len(block.keys)

        def factor(theta, ctx: RateCtx):
            return jnp.ones(n)

    def rate_constants(theta: Mapping[str, Any], ctx: RateCtx):
        return block.k_reactions(block.k_params(theta, factor(theta, ctx), ctx.t_k)) * sf

    stoich, reactants = block.stoichiometry()
    species = _species(reactions, exclude)
    return Module(
        name=name, species=species, elements=_elements(species), reactions=tuple(r.key for r in reactions),
        stoichiometry=stoich, reactants=reactants, params=params, bands=bands, basis=basis,
        rate_constants=rate_constants, observables=observables_for(species), requires=("trunk",),
        needs_ph=False, config_digest=digest, notes=tuple(notes) + (f"zero-k (linear) steps: {block.linear}",),
        pot_inputs=pot_inputs, breakpoints=breaks)


def build(cfg=None) -> Module:
    """The six trunk-only blocks (26 steps) for an engine Config, as one module (``trunk_lane``)."""
    from src.kinetic_core import parameters_methionine as pm
    from src.kinetic_core import parameters_proline as pp
    from src.kinetic_core import parameters_pyrazine as pz

    _, _, _, _, tc = _engine()
    values = config_values(cfg)
    parameters = lane_parameters(cfg)
    if dict(tc.ATMOSPHERE_FACTORS) != {tc.ATMOSPHERE_REFERENCE: 1.0}:
        raise UnsupportedConfiguration(f"trunk_conditions.ATMOSPHERE_FACTORS = {dict(tc.ATMOSPHERE_FACTORS)}: a "
                                       "fitted non-air factor is not ported")
    reactions = _reactions(BLOCKS)
    keys = list(dict.fromkeys(r.parameter_key for r in reactions))
    ph_steps = tuple(pz.PYRAZINE_PH_STEPS) + tuple(pm.METHIONINE_PH_STEPS) + tuple(pp.PROLINE_PH_STEPS)
    if not set(ph_steps) <= set(keys) or not set(tc.OXIDATIVE_ENTRY_STEPS) <= set(keys):
        raise UnsupportedConfiguration("a trunk_conditions step of the pyrazine / atmosphere terms is not a step of "
                                       "the trunk lane's blocks")
    touched = set(tc.AW_STEPS) | set(tc.PH_STEPS)
    if touched & set(keys):
        raise UnsupportedConfiguration(f"{sorted(touched & set(keys))}: a B12 term on a trunk-lane block")
    rows = np.asarray([k in ph_steps for k in keys])
    per_block = {b: len(getattr(_engine()[2], b)) for b in BLOCKS}
    return _block_module(
        NAME, reactions, parameters, cfg, exclude=_trunk_species(), factor_rows=rows, slopes=pz.PYRAZINE_PH_SLOPES,
        digest=values["digest"],
        notes=(f"blocks (steps): {per_block}",
               f"pyrazine pH term on {ph_steps}; atmosphere: air only (refused otherwise)",
               f"methionine.shipped_b22 = {values['methionine_b22']}, proline.shipped_b24 = {values['proline_b24']}",
               "pH mode prescribed (the declared pH, constant over the programme)"))


def build_fed_3deoxy(cfg=None) -> Module:
    """The B39 fed 3-deoxyglucosone triangle (3 steps, owns DGAL), on the trunk lane: no condition term touches it."""
    _, _, network, _, tc = _engine()
    values = config_values(cfg)
    reactions = tuple(network.FED_3DEOXY_REACTIONS)
    keys = {r.parameter_key for r in reactions}
    touched = (set(tc.AW_STEPS) | set(tc.PH_STEPS) | set(tc.OXIDATIVE_ENTRY_STEPS)) & keys
    if touched:
        raise UnsupportedConfiguration(f"{sorted(touched)}: a trunk_conditions term on the fed 3-deoxy block")
    return _block_module(FED_NAME, reactions, lane_parameters(cfg), cfg, exclude=_trunk_species(),
                         digest=values["digest"],
                         notes=(f"dicarbonyl.shipped_b39 = {values['shipped_b39']} (moves no number)",
                                "trunk-lane context: no condition term; the sulfur lane's context is its port's"))


def compose_lane(cfg=None, *, solver: str = "Kvaerno5"):
    """The engine's trunk lane, ``network.TRUNK_REACTIONS`` in its order, as one Network."""
    from . import trunk

    return compose([trunk.build(cfg, lane="trunk"), build(cfg), build_fed_3deoxy(cfg)], solver=solver)


def fold_amines(initial: Mapping[str, float]) -> Dict[str, float]:
    """``engine._integrate_program``'s trunk-lane convention (B24, then B22): proline and methionine are charged as
    glycine too, at the same molarity, and keep their own pools."""
    state = {k: float(v) for k, v in initial.items()}
    if state.get("PRO", 0.0) > 0.0:
        state["Gly"] = state.get("Gly", 0.0) + state["PRO"]
    if state.get("MET", 0.0) > 0.0:
        state["Gly"] = state.get("Gly", 0.0) + state["MET"]
    return state


def lysp_charge(process) -> Optional[float]:
    """The bound-lysine charge ``engine._integrate_program`` adds on the trunk lane (B20): the matrix's amine sites
    from the protein loading, times the available-fraction centre; None without a loading."""
    from src.kinetic_core.matrix_sites import resolve
    from src.kinetic_core.parameters_glycation import available_fraction

    try:
        charged, _ = resolve(process)
    except Exception:  # noqa: BLE001 - the engine's own guard: a malformed loading was refused by the declaration
        return None
    if charged is not None and charged.amine > 0:
        return float(charged.amine) * available_fraction(charged.amine_band)
    return None


def make_pot(net, initial: Mapping[str, float], segments: Sequence[Tuple[float, float]], *, ph: Optional[float],
             water_activity: Optional[float], n_segments: Optional[int] = None,
             atmosphere: Optional[str] = None, lysp: Optional[float] = None) -> Dict[str, Any]:
    """A trunk-lane pot: the DECLARED charge (``{species: mmol/L}``, before the amine folding, which this applies),
    the programme padded to ``n_segments``, the declared pH and a_w (None -> NaN). ``lysp`` is the engine's bound-
    lysine charge (``lysp_charge``), used only when ``initial`` has no LYSP, as the engine does."""
    check_atmosphere(atmosphere)
    state = fold_amines(initial)
    if "LYSP" not in state and lysp is not None:
        state["LYSP"] = float(lysp)
    unknown = sorted(k for k, v in state.items() if float(v) != 0.0 and k not in net.index)
    if unknown:
        raise UnsupportedConfiguration(f"the pot charges {unknown}, which the network does not carry")
    y0 = np.zeros(net.n_species)
    for k, v in state.items():
        if k in net.index:
            y0[net.index[k]] = float(v)
    return {"y0": y0, "segments": program_segments(segments, n_segments), "ph_label": nan_or(ph),
            "aw": nan_or(water_activity)}
