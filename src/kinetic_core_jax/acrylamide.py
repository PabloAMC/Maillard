"""
src/kinetic_core_jax/acrylamide.py -- the acrylamide lane's own block (``acrylamide.ACRYLAMIDE_REACTIONS``, 16 steps)
as a module under the A2 contract v2 (2026-11-05). Composed after ``trunk.build(lane="acrylamide")`` it is the
engine's acrylamide lane: ``FULL_ACRYLAMIDE_REACTIONS = network.REACTIONS + ACRYLAMIDE_REACTIONS``, 42 steps.

RATE LAW (read from the engine): ``AcrylamideParameter.k_at``, k_ref at T_ref = 160 C (each parameter's
``t_ref_k``), R = ``parameters_acrylamide.R_KJ``; ``a_sb_int1`` DERIVED from ``a_asn_glc_sb`` by
``ASN_SCHIFF_AMADORI_SPLIT`` (``acrylamide_rate_constants_at``). Values from ``engine.core_parameters("acrylamide")``
(the measured backbone and the frozen B3 fit).

THE SHIPPED CONDITION TERMS (``acrylamide_conditions.apply``, B14 / B15, and the envelope of CF-039 / CF-106,
decided 2026-11-03), every one ported, read from the pot's declared a_w and pH (pH mode "prescribed": the declared
INITIAL pH, constant over the programme, as the engine's factor is by construction):

* formation a_w multiplier on ``k_int1_acr``: ``aw_scale`` at any declared a_w (1.0 shipped; band AW_SCALE_BAND,
  drawn inside and outside the window since CF-104), exactly 1 at a_w None;
* elimination a_w multiplier on ``k_acr_dp``: 1 - ``aw_elimination_scale`` (1 - m(a_w)), m piecewise-linear through
  AW_ELIMINATION_TABLE inside the window 0.34-0.99 (held at the first point below its first knot, 1 above 0.88), 1
  outside the window and at None;
* pH factors 10^(e (pH - 6.8)) on ``k_asn_glc`` (``ph_exponent_formation``) and ``k_acr_dp``
  (``ph_exponent_elimination``), pH HELD at the window edge outside 4-8, 1 at None;
* the CF-039 Ea draws are the ``ea:k_int1_acr`` / ``ea:k_asn_glc`` coordinates (k(160 C) held, as ``apply``'s
  ``ea_overrides`` does), with the envelope's bands;
* CF-106 (low temperature): a programme whose PEAK lies in (60, 120) C is refused at the pot
  (``check_programme``), as the engine refuses its declaration; at or below 60 C the answer stands.

The lane has no other switch: it reads no ConfigKey (``acrylamide_conditions`` holds constants, not keys) and no pH
state. The thermal history is pot data, so both values of ``panel.thermal_history`` are ported.

SHARED SPECIES. The module owns its 13 species beyond the trunk's, among them ``Cys``, ``CBX``, ``FRAG_N`` and
``FRAG_S``, which the engine declares in the sulfur block of its state. ``compose`` therefore refuses acrylamide +
sulfur (``Cys`` defined twice), as ``engine.resolve_lanes`` refuses the two lanes together; declared competition for
cysteine is A2b's job (contract v2 sec. 6). Lysine is ``Lys`` here and ``LYSP`` on the trunk lane (two keys).
"""
from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import jax.numpy as jnp
import numpy as np

from .lane import UnsupportedConfiguration
from .network import Module, RateCtx
from .trunk import ArrheniusBlock, element_table, nan_or, observables_for, species_of

NAME = "acrylamide"
CONDITION_COORDINATES = ("aw_scale", "aw_elimination_scale", "ph_exponent_formation", "ph_exponent_elimination")


def _engine():
    from src.kinetic_core import acrylamide, acrylamide_conditions, config, engine, parameters_acrylamide

    return acrylamide, acrylamide_conditions, config, engine, parameters_acrylamide


def build(cfg=None) -> Module:
    """The acrylamide block for an engine Config (default: the process's). The lane reads no ConfigKey; the digest is
    recorded so a composed network names the configuration it was built in."""
    lane, ac, config, engine, pa = _engine()
    cfg = config.Config.from_globals() if cfg is None else cfg
    reactions = tuple(lane.ACRYLAMIDE_REACTIONS)
    parameters = engine.core_parameters(engine.ACRYLAMIDE)
    block = ArrheniusBlock(reactions, parameters, lambda p: float(p.t_ref_k), pa.R_KJ)
    derived = [r.key for r in reactions if r.parameter_key is None]
    if derived != ["a_sb_int1"]:
        raise UnsupportedConfiguration(f"derived acrylamide steps {derived}: only a_sb_int1 is ported")
    split = pa.ASN_SCHIFF_AMADORI_SPLIT
    ratio, loading = (float(split["ratio_amadori_over_schiff_pseudo_first_order"]),
                      float(split["amine_loading_mmol_L_for_the_ratio"]))
    j_int1, j_sb = block.rxn_of["a_sb_int1"], block.rxn_of["a_asn_glc_sb"]

    def rows(keys: Sequence[str]) -> np.ndarray:
        missing = [k for k in keys if k not in block.keys]
        if missing:
            raise UnsupportedConfiguration(f"acrylamide_conditions scales {missing}, not steps of this module")
        return np.asarray([k in keys for k in block.keys])

    f_aw_rows, e_aw_rows = rows(ac.SCALED_KEYS), rows(ac.ELIMINATION_KEYS)
    ph_f_rows, ph_e_rows = rows(ac.PH_FORMATION_KEYS), rows(ac.PH_ELIMINATION_KEYS)
    table_a = np.asarray([a for a, _ in ac.AW_ELIMINATION_TABLE])
    table_m = np.asarray([m for _, m in ac.AW_ELIMINATION_TABLE])
    w_lo, w_hi = (float(x) for x in ac.AW_WINDOW)
    p_lo, p_hi = (float(x) for x in ac.PH_WINDOW)
    ref_ph = float(ac.REFERENCE_PH)
    ref_aw = float(ac.REFERENCE_AW)

    def factor(theta, ctx: RateCtx):
        aw_given = ~jnp.isnan(ctx.aw)
        aw = jnp.where(aw_given, ctx.aw, ref_aw)
        f_aw = jnp.where(aw_given, theta["aw_scale"], 1.0)
        inside = aw_given & (aw >= w_lo - 1e-9) & (aw <= w_hi + 1e-9)
        m = jnp.where(aw <= table_a[0], table_m[0], jnp.where(aw >= table_a[-1], 1.0, jnp.interp(aw, table_a, table_m)))
        e_aw = jnp.where(inside, 1.0 - theta["aw_elimination_scale"] * (1.0 - m), 1.0)
        ph_given = ~jnp.isnan(ctx.ph_label)
        p = jnp.clip(jnp.where(ph_given, ctx.ph_label, ref_ph), p_lo, p_hi)
        f_ph = jnp.where(ph_given, 10.0 ** (theta["ph_exponent_formation"] * (p - ref_ph)), 1.0)
        e_ph = jnp.where(ph_given, 10.0 ** (theta["ph_exponent_elimination"] * (p - ref_ph)), 1.0)
        out = jnp.where(jnp.asarray(f_aw_rows), f_aw, 1.0)
        out = out * jnp.where(jnp.asarray(e_aw_rows), e_aw, 1.0)
        out = out * jnp.where(jnp.asarray(ph_f_rows), f_ph, 1.0)
        return out * jnp.where(jnp.asarray(ph_e_rows), e_ph, 1.0)

    def rate_constants(theta: Mapping[str, Any], ctx: RateCtx):
        k = block.k_reactions(block.k_params(theta, factor(theta, ctx), ctx.t_k))
        return k.at[j_int1].set(ratio * k[j_sb] * loading)

    params = block.coordinates()
    ea_bands = {"k_int1_acr": tuple(float(x) for x in ac.EA_INT1_BAND_KJ_MOL)}
    bands, basis = block.bands_and_basis(ea_band_override=ea_bands, log10_k_bounds=pa.FITTED_ACRYLAMIDE_BOUNDS_LOG10K)
    params.update({"aw_scale": float(ac.AW_MULTIPLIER), "aw_elimination_scale": 1.0,
                   "ph_exponent_formation": float(ac.PH_EXPONENT_FORMATION),
                   "ph_exponent_elimination": float(ac.PH_EXPONENT_ELIMINATION)})
    bands.update({"aw_scale": tuple(ac.AW_SCALE_BAND), "aw_elimination_scale": tuple(ac.AW_ELIMINATION_SCALE_BAND),
                  "ph_exponent_formation": tuple(ac.PH_EXPONENT_FORMATION_BAND),
                  "ph_exponent_elimination": tuple(ac.PH_EXPONENT_ELIMINATION_BAND)})
    basis.update({
        "aw_scale": f"acrylamide.aw_multiplier (B14 / CF-104): {ac.AW_SOURCE[:80]}",
        "aw_elimination_scale": "acrylamide.aw_elimination_deficit_scale (B15): De Vleeschouwer 2007 k_E shape",
        "ph_exponent_formation": f"acrylamide.ph_exponent_formation_decades_per_unit (B15): {ac.PH_SOURCE[:70]}",
        "ph_exponent_elimination": f"acrylamide.ph_exponent_elimination_decades_per_unit (B15): {ac.PH_SOURCE[:60]}",
        "ea:k_int1_acr": f"acrylamide.k_int1_acr.ea_kj_mol (CF-039, k(160 C) held): {ac.EA_INT1_SOURCE[:70]}"})
    stoich, reactants = block.stoichiometry()
    from src.kinetic_core import network as engine_network
    from src.kinetic_core.species_acrylamide import acrylamide_ppb

    trunk_species = set(species_of(engine_network.REACTIONS))
    species = tuple(s for s in species_of(reactions) if s not in trunk_species)
    obs = observables_for(species)

    per_ppb = float(acrylamide_ppb(1.0))
    obs["ppb:ACR"] = (lambda y, index: y[index["ACR"]] * per_ppb)
    return Module(
        name=NAME, species=species, elements=element_table(species), reactions=tuple(r.key for r in reactions),
        stoichiometry=stoich, reactants=reactants, params=params, bands=bands, basis=basis,
        rate_constants=rate_constants, observables=obs, requires=("trunk",), needs_ph=False,
        config_digest=cfg.digest(), pot_inputs=("ph_label", "aw"),
        notes=("acrylamide_conditions B14 / B15 / CF-039 / CF-104 ported; CF-106 refused at the pot",
               "owns Cys, CBX, FRAG_N, FRAG_S: refuses to compose with a sulfur module (engine.resolve_lanes)"))


def check_programme(segments: Sequence[Tuple[float, float]]) -> None:
    """CF-106: refuse a programme whose peak lies inside (60, 120) C (``engine.declare_envelope``'s rule, on
    ``ThermalProgram.peak_temperature_c``, the max over its segments)."""
    _, ac, *_ = _engine()
    low, high = ac.LOW_T_REFUSAL_C
    peak = max(float(t) for _, t in segments)
    if low < peak < high:
        raise UnsupportedConfiguration(
            f"ACRYLAMIDE BELOW ITS FIT WINDOW (CF-106): the programme peaks at {peak:.0f} C, inside {low:.0f}-"
            f"{high:.0f} C; the engine refuses it rather than extrapolate")


def make_pot(net, initial: Mapping[str, float], segments: Sequence[Tuple[float, float]], *,
             ph: Optional[float], water_activity: Optional[float], n_segments: Optional[int] = None,
             atmosphere: Optional[str] = None) -> Dict[str, Any]:
    """A pot for a network that holds this module: the initial state from ``{species: mmol/L}`` (a charged species
    the network does not carry is refused, never dropped), the programme (CF-106 checked; padded to ``n_segments``)
    and the declared pH and a_w (None -> NaN)."""
    from .network import program_segments
    from .trunk import check_atmosphere

    if any(m.name == NAME for m in net.modules):
        check_programme(segments)
    if any(m.name == "trunk" and "lane context: trunk" in m.notes for m in net.modules):
        check_atmosphere(atmosphere)
    unknown = sorted(k for k, v in initial.items() if float(v) != 0.0 and k not in net.index)
    if unknown:
        raise UnsupportedConfiguration(f"the pot charges {unknown}, which the network does not carry")
    y0 = np.zeros(net.n_species)
    for k, v in initial.items():
        if k in net.index:
            y0[net.index[k]] = float(v)
    return {"y0": y0, "segments": program_segments(segments, n_segments), "ph_label": nan_or(ph),
            "aw": nan_or(water_activity)}
