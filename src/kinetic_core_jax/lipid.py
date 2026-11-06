"""
src/kinetic_core_jax/lipid.py -- the lipid lane as the first module under the A2 contract v2 (2026-11-05).

The engine (``src/kinetic_core/lipid.py``) solves its lane ANALYTICALLY per segment: five hydroperoxide pools decay
first order with one constant k(T) (Schroen's 25 C anchor times the declared Q10), the decomposed linoleate becomes
the measured six-product slate at Frankel's frozen branch shares (scaled to Schroen's 1.2 % hexanal), 2-pentylfuran
follows hexanal by Frankel 1981's ratio, nonanal follows the oleate pool on the declared molar anchor, and above
100 C 2,4-decadienal breaks first order with a fixed share becoming hexanal (B48); the carbon no product names goes to
LIPID_FRAG_C. That arithmetic is the exact solution of a linear mass-action network, which is what this module
declares, so the lane composes with any other module and is integrated by the same diffrax solve.

STRUCTURE (static stoichiometry; every free number in a rate constant):

* ``lipid_decomp_<pool>``: LOOH -> 19 LIPID_FRAG_C, k(T), for the four linoleate cells and LOOH_OL;
* ``lipid_slate_<pool>_<product>``: a YIELD CHANNEL, LOOH -> LOOH + product - C(product) LIPID_FRAG_C, at
  k(T) x y(pool, product; donor). Each linoleate cell carries its OWN simplex (Frankel's (position, geometry) cells);
  the engine applies the pool-averaged yields to the summed decomposition, which is the same number because every
  cell decays with one k and the engine charges the cells in the composition's proportions;
* ``lipid_pentylfuran_<pool>``: the same channel for 2-pentylfuran, at the cell's hexanal yield x the 1981 ratio;
* ``lipid_nonanal_LOOH_OL``: the oleate pool's nonanal channel, at share x anchor x the named-slate total of the
  default pool composition (``engine.core_lipid_model``), as ``lipid._b28_yields`` computes it;
* ``lipid_decadienal_loss`` / ``lipid_decadienal_hexanal``: decadienal -> 10 LIPID_FRAG_C at k_brk(T), and its hexanal
  channel at share_hexanal x k_brk(T); k_brk is 0 below the breakage window (100 C), as the engine switches it off.

Everything is read from the engine at build time: the species and their carbon, the frozen branch model and pool
composition, the per-cell yields (``lipid.slate_yields`` on a single-cell composition), the 1981 ratios
(``lipid._b28_yields`` on unit inputs), the anchor, the Q10 and the breakage constants from the Config. The one
literal re-typed is the gas constant inside ``DecadienalBreakage.k_loss_per_min`` (CODATA 2018), pinned by the rate
test to 1e-12.

Free coordinates (the pytree; shipped values from the Config): ``q10``; ``log10_k_looh`` (a log10 multiplier on the
25 C anchor, 0 shipped); ``donor`` (the hydrogen-donor suppression of the homolytic channel, 0 shipped: the engine
stores no value); ``log10_k_decadienal`` (multiplier on k_brk at its reference), ``ea_decadienal_kj_mol`` and
``share_hexanal_from_decadienal``. Refused: an enabled covalent sink (the engine raises too).
"""
from __future__ import annotations

from typing import Any, Dict, Mapping, Optional, Tuple

import jax.numpy as jnp
import numpy as np

from .lane import UnsupportedConfiguration
from .network import CELSIUS, Module, RateCtx

NAME = "lipid"
#: DecadienalBreakage.k_loss_per_min's literal gas constant, J / (mol K) (CODATA 2018); the rate test pins it
R_J_PER_MOL_K = 8.314462618
COORDINATES = ("q10", "log10_k_looh", "donor", "log10_k_decadienal", "ea_decadienal_kj_mol",
               "share_hexanal_from_decadienal")


def _engine():
    from src.kinetic_core import config, engine
    from src.kinetic_core import lipid as lane
    from src.kinetic_core import parameters_lipid as pl
    from src.kinetic_core import species_lipid as sl

    return config, engine, lane, pl, sl


def config_values(cfg=None) -> Dict[str, Any]:
    """The lipid keys of an engine Config (default: the process's, ``Config.from_globals``), refused where the port
    does not implement them."""
    config, *_ = _engine()
    cfg = config.Config.from_globals() if cfg is None else cfg
    sink = cfg["lipid.covalent_sink"]
    if sink.enabled:
        raise UnsupportedConfiguration("lipid.covalent_sink enabled: the engine raises NotImplementedError and the "
                                       "port has no storage clock either")
    return {"q10": cfg["lipid.q10_assumption"], "breakage": cfg["lipid.decadienal_breakage"], "sink": sink,
            "digest": cfg.digest()}


def cell_yields(donor: float = 0.0) -> Dict[str, Dict[str, float]]:
    """The ENGINE's yields per mole of each linoleate cell decomposed (``lipid.slate_yields`` on a composition that is
    all that cell), at ``donor``: the per-cell stoichiometry this module uses."""
    _, engine, lane, _, sl = _engine()
    branch, _ = engine.core_lipid_model()
    out = {}
    for pool in sl.LOOH_POOLS:
        cells = {k: (1.0 if k == pool else 0.0) for k in sl.LOOH_POOLS}
        comp = lane.LOOHComposition(f13_ct=cells["LOOH_13_ct"], f13_tt=cells["LOOH_13_tt"],
                                    f9_ct=cells["LOOH_9_ct"], f9_tt=cells["LOOH_9_tt"])
        out[pool] = lane.slate_yields(branch, comp, donor=donor)
    return out


def build(cfg=None) -> Module:
    """The lipid Module for an engine Config (default: the process's)."""
    _, engine, lane, pl, sl = _engine()
    values = config_values(cfg)
    q10, brk = values["q10"], values["breakage"]
    branch, composition = engine.core_lipid_model()
    carbon = {s.key: int(s.carbon) for s in sl.LIPID_SPECIES}
    elements = {s.key: (int(s.carbon), int(s.nitrogen), int(s.sulfur)) for s in sl.LIPID_SPECIES}
    pools = tuple(sl.LOOH_POOLS)
    base = cell_yields(0.0)
    homolytic = {p: sl.CLEAVAGE_MECHANISM[p] == "homolytic" for p in sl.FRANKEL_SLATE}
    pentylfuran_ratio = float(lane._b28_yields({"HEXANAL": 1.0})["PENTYLFURAN"])
    nonanal_per_named = float(lane._b28_yields({"_": 1.0})["NONANAL_PER_OLEATE_LOOH"])
    weights = composition.as_dict()

    reactions, stoich, reactants = [], {}, {}
    kinds = []   # per reaction: ("decomp",) | ("slate", pool, product) | ("pentylfuran", pool) | ("nonanal",) | ...

    def add(key, coeffs, reactant, kind):
        reactions.append(key)
        stoich[key] = coeffs
        reactants[key] = {reactant: 1}
        kinds.append(kind)

    frag = "LIPID_FRAG_C"
    for pool in pools + ("LOOH_OL",):
        add(f"lipid_decomp_{pool}", {pool: -1.0, frag: float(carbon[pool])}, pool, ("decomp",))
    for pool, (position, _geometry) in sl.LOOH_POOLS.items():
        for product in sl.POSITION_PRODUCTS[position]:
            add(f"lipid_slate_{pool}_{product}", {product: 1.0, frag: -float(carbon[product])}, pool,
                ("slate", pool, product))
        if "HEXANAL" in sl.POSITION_PRODUCTS[position]:
            add(f"lipid_pentylfuran_{pool}", {"PENTYLFURAN": 1.0, frag: -float(carbon["PENTYLFURAN"])}, pool,
                ("pentylfuran", pool))
    add("lipid_nonanal_LOOH_OL", {"NONANAL": 1.0, frag: -float(carbon["NONANAL"])}, "LOOH_OL", ("nonanal",))
    if brk.enabled:
        add("lipid_decadienal_loss", {"DECADIENAL": -1.0, frag: float(carbon["DECADIENAL"])}, "DECADIENAL",
            ("brk_loss",))
        add("lipid_decadienal_hexanal", {"HEXANAL": 1.0, frag: -float(carbon["HEXANAL"])}, "DECADIENAL",
            ("brk_hexanal",))

    # static arrays of the rate function: k_j = k_looh x (a_j + b_j (1 - donor)) for the LOOH-driven reactions,
    # k_brk x (1 | share) for the two decadienal ones
    n = len(reactions)
    a = np.zeros(n)          # donor-independent part of the yield factor
    b = np.zeros(n)          # homolytic part (multiplied by 1 - donor)
    looh_driven = np.zeros(n, dtype=bool)
    for j, kind in enumerate(kinds):
        if kind[0] == "decomp":
            a[j], looh_driven[j] = 1.0, True
        elif kind[0] in ("slate", "pentylfuran"):
            pool = kind[1]
            product = kind[2] if kind[0] == "slate" else "HEXANAL"
            factor = 1.0 if kind[0] == "slate" else pentylfuran_ratio
            y = float(base[pool][product]) * factor
            if homolytic[product]:
                b[j] = y
            else:
                a[j] = y
            looh_driven[j] = True
        elif kind[0] == "nonanal":   # share x anchor x SUM_p yields_p(default composition, donor)
            for pool in pools:
                for product, y in base[pool].items():
                    w = float(weights[pool]) * float(y) * nonanal_per_named
                    if homolytic[product]:
                        b[j] += w
                    else:
                        a[j] += w
            looh_driven[j] = True
    brk_loss = np.array([k[0] == "brk_loss" for k in kinds])
    brk_hex = np.array([k[0] == "brk_hexanal" for k in kinds])
    anchor_per_hour = float(pl.K_LOOH_DECOMP_ANCHOR.value)
    t_ref_q10 = float(q10.reference_temperature_c)
    window = float(brk.window_min_c)
    k_brk_ref = float(brk.k_loss_per_min_at_ref)
    inv_t_brk_ref = 1.0 / (float(brk.t_ref_c) + CELSIUS)
    a_j, b_j = jnp.asarray(a), jnp.asarray(b)
    driven, loss, hexanal = jnp.asarray(looh_driven), jnp.asarray(brk_loss), jnp.asarray(brk_hex)

    def k_looh(theta, ctx: RateCtx):
        """parameters_lipid.k_looh_decomp_per_min, times the log10 multiplier."""
        per_hour = anchor_per_hour * theta["q10"] ** ((ctx.t_c - t_ref_q10) / 10.0)
        return per_hour / 60.0 * 10.0 ** theta["log10_k_looh"]

    def k_brk(theta, ctx: RateCtx):
        """DecadienalBreakage.k_loss_per_min, zero below its window (lipid.integrate_lipid's switch)."""
        k = (k_brk_ref * 10.0 ** theta["log10_k_decadienal"]
             * jnp.exp(-theta["ea_decadienal_kj_mol"] * 1000.0 / R_J_PER_MOL_K * (1.0 / ctx.t_k - inv_t_brk_ref)))
        return jnp.where(ctx.t_c >= window, k, 0.0)

    def rate_constants(theta: Mapping[str, Any], ctx: RateCtx):
        kl = k_looh(theta, ctx)
        kb = k_brk(theta, ctx) if brk.enabled else 0.0
        looh = kl * (a_j + b_j * (1.0 - theta["donor"]))
        brk_part = jnp.where(loss, kb, jnp.where(hexanal, kb * theta["share_hexanal_from_decadienal"], 0.0))
        return jnp.where(driven, looh, brk_part)

    mw = dict(sl.MOLECULAR_WEIGHT_G_PER_MOL)
    observables = {
        f"ug_per_L:{key}": (lambda y, index, key=key: y[index[key]] * mw[key] * 1.0e3) for key in mw
    }
    params = {"q10": float(q10.default), "log10_k_looh": 0.0, "donor": 0.0, "log10_k_decadienal": 0.0,
              "ea_decadienal_kj_mol": float(brk.ea_loss_kj_mol),
              "share_hexanal_from_decadienal": float(brk.share_hexanal)}
    q10_band = dict(_key("lipid.q10_assumption").bands).get("default", (float(q10.lo), float(q10.hi)))
    bands: Dict[str, Optional[Tuple[float, float]]] = {
        "q10": tuple(q10_band), "log10_k_looh": None, "donor": (0.0, 1.0), "log10_k_decadienal": None,
        "ea_decadienal_kj_mol": None, "share_hexanal_from_decadienal": (0.0, 1.0)}
    basis = {
        "q10": "lipid.q10_assumption (Schroen & Berton-Carabin 2022; band Tazi 2009 / Frankel 1993, config._BANDS)",
        "log10_k_looh": "a multiplier on K_LOOH_DECOMP_ANCHOR (Schroen 2022 k4, hand-fitted, no SE): no band declared",
        "donor": "BranchModel.slate_flux's hydrogen-donor suppression; d in [0, 1); no stored value (0 shipped)",
        "log10_k_decadienal": "a multiplier on DecadienalBreakage.k_loss_per_min_at_ref (Zamora 2015b Fig. 5C read)",
        "ea_decadienal_kj_mol": "DecadienalBreakage.ea_loss_kj_mol (Zamora 2015b, ~25-27 kJ/mol); no band declared",
        "share_hexanal_from_decadienal": "DecadienalBreakage.share_hexanal (Zamora 2015b sec. 3.5, 11.5 %)"}
    assert tuple(params) == COORDINATES
    return Module(
        name=NAME, species=tuple(sl.LIPID_KEYS), elements=elements, reactions=tuple(reactions),
        stoichiometry=stoich, reactants=reactants, params=params, bands=bands, basis=basis,
        rate_constants=rate_constants, observables=observables, requires=(), needs_ph=False,
        config_digest=values["digest"],
        breakpoints={"t_c": (window,)} if brk.enabled else {},     # k_brk jumps from 0 at the window (a jump)
        notes=(f"breakage={'on' if brk.enabled else 'off'} (window {window:g} C)",
               f"composition={ {k: round(float(v), 6) for k, v in weights.items()} }",
               "covalent sink off (required)"))


def _key(name: str):
    from src.kinetic_core import config

    return {k.name: k for k in config.KEYS}[name]


# --------------------------------------------------------------------------- pots
def initial_state(carriers: Mapping[str, Mapping[str, Optional[float]]], species: Tuple[str, ...]) -> np.ndarray:
    """The summed initial state of a pot's carriers, ``{carrier_key: {"lipid_fraction": f | None,
    "peroxide_value_meq_per_kg": pv | None}}`` (None = the carrier's declared centre), built by the engine's own
    ``lipid.charge_from_carrier`` and ``species_lipid.initial_lipid_state`` with the default composition. The engine
    integrates each carrier separately and sums the states; the network is linear, so one integration of the summed
    state is the same answer."""
    _, engine, lane, pl, sl = _engine()
    _, composition = engine.core_lipid_model()
    total = np.zeros(len(species))
    for key, over in carriers.items():
        charge = lane.charge_from_carrier(pl.LIPID_CARRIERS[key], composition,
                                          lipid_fraction=over.get("lipid_fraction"),
                                          peroxide_value_meq_per_kg=over.get("peroxide_value_meq_per_kg"))
        state = sl.initial_lipid_state({
            "LOOH_13_ct": charge.looh_linoleate_mmol_l * charge.composition.f13_ct,
            "LOOH_13_tt": charge.looh_linoleate_mmol_l * charge.composition.f13_tt,
            "LOOH_9_ct": charge.looh_linoleate_mmol_l * charge.composition.f9_ct,
            "LOOH_9_tt": charge.looh_linoleate_mmol_l * charge.composition.f9_tt,
            "LOOH_OL": charge.looh_oleate_mmol_l})
        total += np.array([float(state.get(s, 0.0)) for s in species])
    return total
