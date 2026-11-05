"""
src/kinetic_core_jax/lane.py -- the STATIC half of the port: everything the JAX model needs that does not depend on
the state or on the free parameters, read from the engine (``src.kinetic_core``) as data and frozen into numpy
arrays. Plain numpy, no JAX.

What is read, not re-derived: the species order, the stoichiometry matrix, the mass-action layout, every rate
constant's k_ref / Ea / reference temperature / temperature cap / pH-factor kind, the derived-rate constants of the
switches in force, the pKa tables of the charge balance, the buffer, the pH-drift calibration. What is re-implemented
(in ``model.py``): the arithmetic that turns those into rates -- Arrhenius, the pH factors, the measured pH tables,
the derived rates, the charge balance.

``build`` reads the engine configuration INSTALLED at call time (``src.kinetic_core.config``) and refuses any switch
setting the port does not implement (``UnsupportedConfiguration``), so it can never compute a different model in
silence.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, Mapping, Sequence, Tuple

import numpy as np

#: pH-factor kinds, as integer codes (parameters_sulfur.ph_factor's six kinds plus "none")
KIND_CODES = {"": 0, "acid": 1, "base": 2, "neutral_h2s": 3, "hs_anion": 4, "thiolate": 5, "hydroxide": 6}
#: free coordinates the port maps (the R3En joint-fit set and its sub-models' cousins); offsets are log10 multipliers
#: on the shipped k_ref of the parameter keys listed (fit_joint_deoxyosone_thiol.MULTIPLIED / ARM_MULTIPLIED)
OFFSET_KEYS = {
    "offset_k_dpo_ddp": ("k_dpo_ddp",), "offset_k_dpo_c2c3": ("k_dpo_c2c3",), "offset_k_dpo_nf": ("k_dpo_nf",),
    "offset_k_nf_mft": ("k_nf_mft",), "offset_fft_formation": ("k_fur_fft", "k_fur_fft_hs"),
    "offset_k_h2s_loss": ("k_h2s_loss",), "offset_k_mft_decay": ("k_mft_decay",),
    "offset_k_cys_actz": ("k_cys_actz",), "offset_k_aa_frag": ("k_aa_frag",), "offset_k_ttca_deg": ("k_ttca_deg",),
    "offset_k_pent_caramel": ("k_pent_caramel",),
}
#: the DDP route's product coordinate: k_ddp_mft / k_ddp_mft_hs move by (offset_ddp_route - offset_k_dpo_ddp)
DDP_ROUTE = ("offset_ddp_route", ("k_ddp_mft", "k_ddp_mft_hs"))
#: coordinates handled by name in model.rate_constants
NAMED_KEYS = ("mgo_share_nf", "log10_k_pent_amine", "entry_ph_slope", "log10_k_dpo_aa", "aa_ph_slope", "amine_return")
#: reaction keys whose constant the engine derives after the generic loop (sulfur.sulfur_rate_constants_at)
DERIVED = ("r_amadori", "r_cys_h2s", "r_gsh_h2s", "r_thi_h2s", "r_ttca_form", "r_cama_form", "r_cama_dpo",
           "r_cama_tdp", "r_dpo_ddp_h2s", "r_nf_decay_oh", "r_nf_decay_oh_mgo", "r_tdp_mgo", "r_pent_mgo",
           "r_mgo_cys_strecker", "r_mgo_cys_h2s", "r_tdp_cys_trap", "r_odg_dpo", "r_pent_amine", "r_pama_dpo",
           "r_pama_tdp", "r_pama_dpo_n", "r_pama_tdp_n", "r_dpo_aa", "ch_thioether_release",
           "ch_thioether_release_fft", "ch_oligomer_mft", "ch_oligomer_fft")
PAD_PK = 1.0e3   # a padded protonation step that never ionises


class UnsupportedConfiguration(ValueError):
    """The installed engine configuration has a switch the port does not implement."""


@dataclass(frozen=True)
class Lane:
    """The frozen arrays of one lane configuration, one pot's buffer and drift, and the free coordinates."""

    species: Tuple[str, ...]
    reactions: Tuple[str, ...]
    stoichiometry: np.ndarray          # (n_species, n_reactions)
    reactant_index: np.ndarray         # (n_reactions, max_reactants); padded with n_species (a 1.0 appended)
    reactant_power: np.ndarray         # same shape; 0 for padding
    k_ref: np.ndarray                  # shipped k_ref (0 for a derived reaction)
    ea_over_r: np.ndarray              # Ea / R in K (0 when no Ea)
    inv_t_ref: np.ndarray              # 1 / T_ref (unused when no Ea)
    cap_k: np.ndarray                  # temperature ceiling (inf when none)
    kind: np.ndarray                   # pH-factor code per reaction
    coordinate_keys: Tuple[str, ...]
    offset_matrix: np.ndarray          # (n_reactions, n_coordinates): d log10 k / d coordinate for offsets
    index: Mapping[str, int]
    const: Mapping[str, Any]           # scalar constants of the derived rates and the pH physics
    # the charge balance: titratable groups (rows) as a linear map of the state
    group_names: Tuple[str, ...]
    group_pk25: np.ndarray             # (n_groups, 4), padded with PAD_PK
    group_dh: np.ndarray
    group_dcp: np.ndarray
    group_zmax: np.ndarray
    group_weights: np.ndarray          # (n_groups, n_species): mmol/L of the group per mmol/L of the species
    group_offset: np.ndarray           # (n_groups,): the buffer's own mmol/L
    notes: Tuple[str, ...] = field(default_factory=tuple)

    @property
    def n_species(self) -> int:
        return len(self.species)


def _require(condition: bool, what: str) -> None:
    if not condition:
        raise UnsupportedConfiguration(what)


def _check_configuration() -> Dict[str, Any]:
    """The switches in force, refused unless the port implements them; returns the values model.py reads."""
    from src.kinetic_core import parameters_sulfur as ps
    from src.kinetic_core import ph_state, sulfur

    b51 = ph_state.B51
    _require(all(b51.as_dict().values()), f"B51 must be all on (shipped), got {b51.as_dict()}")
    _require(sulfur.PH_COUPLING == "inline", f"ph.coupling must be 'inline', got {sulfur.PH_COUPLING!r}")
    _require(not sulfur.TRUNK_PH_TERMS, "sulfur.TRUNK_PH_TERMS (sub-model 1) is not ported")
    _require(ps.CYS_AMADORI_ROUTE.enabled and ps.CYS_AMADORI_ROUTE.split_law == "none",
             "the port needs B59b on with split_law 'none'")
    _require(ps.CYS_AMADORI_ROUTE.above_measured_ph == "extend_last_slope", "CR-14 'flat' is not ported")
    _require(ps.TTCA_FORMATION.enabled, "the port needs B58 on")
    _require(ps.STRECKER_SULFUR.enabled, "the port needs B59a on")
    _require(ps.NF_HYDROXIDE_DECAY.enabled, "the port needs B61 on")
    _require(not ps.TDP_OWN_DECAY.enabled, "B7a R-T (TDP_OWN_DECAY) is not ported")
    _require(ps.DDP_FORMATION_PH.law in ("none", "base"), f"DDP law {ps.DDP_FORMATION_PH.law!r} is not ported")
    _require(ps.PENTOSE_MGO_ROUTES.k_tdp_mgo_100c_ph5 == 0.0 and ps.PENTOSE_MGO_ROUTES.k_pent_mgo_100c_ph5 == 0.0,
             "the pentose MGO routes (sub-model 3 A3 / A4) are not ported")
    _require(ps.MGO_CYSTEINE_STRECKER.k_l_per_mmol_min_100c_ph5 == 0.0, "R2a's MGO Strecker is not ported")
    _require(ps.TDP_CYSTEINE_TRAP.k_l_per_mmol_min_100c_ph5 == 0.0, "the 3-DX cysteine trap is not ported")
    _require(ps.HEXOSE_C5_ENTRY.k_per_min_100c == 0.0, "the hexose C5 entry is not ported")
    _require(ps.PENTOSE_AMINE_ENTRY.enabled and ps.PENTOSE_AMINE_ENTRY.aa_law == "slope",
             "the port needs the pentose + amine entry on with aa_law 'slope'")
    _require(ps.oligomerisation_rate() == 0.0, "k_oligomer is expected to be 0")
    return {"ddp_law": ps.DDP_FORMATION_PH.law, "trunk_solutes": bool(ph_state.TRUNK_SOLUTES_IN_BALANCE)}


def _charge_balance_groups(species_index, drift, buffer_spec, trunk_solutes):
    """titratable_inventory as a linear map of the state (ph_state.titratable_inventory, group by group)."""
    from src.kinetic_core import ph_state as ph

    n = len(species_index)
    rows: Dict[str, Tuple[Any, np.ndarray, float]] = {}

    def add(group, species=None, weight=1.0, offset=0.0, name=None):
        name = name or group.key
        if name not in rows:
            rows[name] = (group, np.zeros(n), 0.0)
        g, w, b = rows[name]
        if species is not None:
            w[species_index[species]] += weight
        rows[name] = (g, w, b + offset)

    y = float(drift.acid_yield)
    for group, share in ph.ORGANIC_ACID_MIX:   # the lumped pool (scaled) and the carried carboxyls (full strength)
        add(group, "ACID", share * y)
        add(group, "CBX", share)
    add(ph.FORMIC, "FA")
    add(ph.ACETIC, "AA")
    add(ph.cysteine_group(), "Cys")
    arp = drift.arp_group()
    add(arp, "ARP", name="arp")
    add(arp, "CAMA", name="arp")
    add(ph.TTCA_GROUP, "TTCA")
    if trunk_solutes:
        add(ph.GLYCINE, "Gly")
        for key in ("AMA", "SB", "PAMA"):
            add(ph.DFG_GROUP, key, name="dfg")
    for group, mmol in buffer_spec.groups():
        add(group, None, offset=float(mmol), name=f"buffer_{group.key}")
    names = tuple(rows)
    width = 4
    pk = np.full((len(names), width), PAD_PK)
    dh = np.zeros((len(names), width))
    dcp = np.zeros((len(names), width))
    zmax = np.zeros(len(names))
    weights = np.zeros((len(names), n))
    offset = np.zeros(len(names))
    for r, name in enumerate(names):
        group, w, b = rows[name]
        k = len(group.pkas)
        _require(k <= width, f"group {group.key} has {k} steps")
        pk[r, :k] = group.pkas
        dh[r, :k] = group.delta_h_kj_mol
        dcp[r, :k] = group.delta_cp_j_mol_k or (0.0,) * k
        zmax[r] = group.z_max
        weights[r] = w
        offset[r] = b
    return names, pk, dh, dcp, zmax, weights, offset


def build(parameters: Mapping[str, Any], coordinate_keys: Sequence[str], buffer_spec, drift) -> Lane:
    """The Lane for the installed configuration, the SHIPPED-k_ref ``parameters`` (engine.core_parameters(SULFUR);
    a trial's multipliers are the coordinates, applied in JAX), one pot's buffer and drift."""
    from src.kinetic_core import parameters as trunk_parameters
    from src.kinetic_core import parameters_sulfur as ps
    from src.kinetic_core import ph_state as ph
    from src.kinetic_core import sulfur
    from src.kinetic_core.parameters import KineticParameter
    from src.kinetic_core.parameters_furanic import TEMPERATURE_CAP_K

    switches = _check_configuration()
    species = tuple(sulfur.SULFUR_STATE_KEYS)
    reactions = tuple(r.key for r in sulfur.FULL_REACTIONS)
    n_r = len(reactions)
    index = {k: j for j, k in enumerate(reactions)}
    unknown = [k for k in coordinate_keys if k not in OFFSET_KEYS and k != DDP_ROUTE[0] and k not in NAMED_KEYS]
    _require(not unknown, f"coordinates not ported: {unknown}")
    _require(DDP_ROUTE[0] not in coordinate_keys or "offset_k_dpo_ddp" in coordinate_keys,
             "the DDP route coordinate needs offset_k_dpo_ddp")

    k_ref = np.zeros(n_r)
    ea_over_r = np.zeros(n_r)
    inv_t_ref = np.zeros(n_r)
    cap_k = np.full(n_r, np.inf)
    kind = np.zeros(n_r, dtype=np.int64)
    offset_matrix = np.zeros((n_r, len(coordinate_keys)))
    col = {k: c for c, k in enumerate(coordinate_keys)}
    for j, reaction in enumerate(sulfur.FULL_REACTIONS):
        key = reaction.parameter_key
        kind[j] = KIND_CODES[sulfur.REACTION_PH_FACTOR.get(reaction.key, "")]
        if key is None or reaction.key in DERIVED or key in ("k_cys_h2s", "k_gsh_h2s", "k_thi_h2s", "k_ttca_form",
                                                             "k_dpo_ddp_h2s", "k_oligomer"):
            continue
        p = parameters[key]
        _require(p.k_ref is not None, f"{key} unpopulated")
        k_ref[j] = float(p.k_ref)
        if p.ea_kj_mol is not None:
            ea_over_r[j] = float(p.ea_kj_mol) / ps.R_KJ
            inv_t_ref[j] = 1.0 / (float(trunk_parameters.T_REF_K) if isinstance(p, KineticParameter)
                                  else float(p.t_ref_k))
        if key in TEMPERATURE_CAP_K:
            cap_k[j] = float(TEMPERATURE_CAP_K[key])
        for coord, keys in OFFSET_KEYS.items():
            if coord in col and key in keys:
                offset_matrix[j, col[coord]] += 1.0
        if DDP_ROUTE[0] in col and key in DDP_ROUTE[1]:
            offset_matrix[j, col[DDP_ROUTE[0]]] += 1.0
            offset_matrix[j, col["offset_k_dpo_ddp"]] -= 1.0
    if switches["ddp_law"] == "none":   # sulfur_rate_constants_at: r_dpo_ddp / base x 1
        kind[index["r_dpo_ddp"]] = KIND_CODES[""]

    layout = [[(sulfur.SULFUR_INDEX[s], int(c)) for s, c in r.reactants.items()] for r in sulfur.FULL_REACTIONS]
    width = max(len(x) for x in layout)
    r_idx = np.full((n_r, width), len(species), dtype=np.int64)
    r_pow = np.zeros((n_r, width), dtype=np.int64)
    for j, lay in enumerate(layout):
        for q, (i, c) in enumerate(lay):
            _require(c in (1, 2), f"{reactions[j]}: reactant order {c} not ported")
            r_idx[j, q], r_pow[j, q] = i, c

    zheng = ps.ZHENG_CYSTEINE_THERMOLYSIS
    entry = ps.PENTOSE_AMINE_ENTRY
    route = ps.CYS_AMADORI_ROUTE
    nf = ps.NF_HYDROXIDE_DECAY
    from src.kinetic_core.parameters import SCHIFF_AMADORI_SPLIT

    const: Dict[str, Any] = {
        "R_KJ": ps.R_KJ,
        "ph_norm": ps.PH_NORMALISATION,
        # pKa of the speciation factors (B51: Suleimenov H2S, GKL cysteine pair)
        "pk_thiol": tuple(ps.PKA_GKL_CYSTEINE["thiol"]),
        "pk_amine": tuple(ps.PKA_GKL_CYSTEINE["cysteine_amine"]),
        "pka_ref_t_k": ph.PKA_REFERENCE_T_K,
        # Zheng & Ho cysteine thermolysis (log10 A in 1/min; Ea), the mechanism boundary
        "zheng_ph": tuple(zheng["ph_points"]),
        "zheng_log10_a": tuple(math.log10(a * 60.0) for a in zheng["prefactor_per_second"]),
        "zheng_ea": tuple(zheng["ea_kj_mol"]),
        "zheng_boundary": ps.ZHENG_MECHANISM_BOUNDARY_PH,
        "gsh_pairs": tuple(ps.GLUTATHIONE_ARRHENIUS_REFIT[p] for p in (3.0, 5.0, 7.0, 9.0)),
        "thi_share_points": tuple(ps.THIAMINE_H2S_SHARE_POINTS),
        "thi_share_floor": ps.THIAMINE_H2S_SHARE_FLOOR_PH,
        "ttca_k": float(ps.TTCA_FORMATION.k_l_per_mmol),
        "cama_k_form_ref": route.k_form_ref_per_min, "cama_t_form_ref_k": route.t_form_ref_c + 273.15,
        "cama_ea_form": route.ea_form_kj_mol,
        "cama_decay_ph": tuple(p for p, _ in route.decay_ph_table),
        "cama_decay_log10k": tuple(math.log10(k) for _, k in route.decay_ph_table),
        "cama_t_decay_ref_k": route.t_decay_ref_c + 273.15, "cama_ea_decay": route.ea_decay_kj_mol,
        "cama_dpo_share": float(route.dpo_share),
        "strecker_share": float(ps.STRECKER_SULFUR.h2s_share),
        "nf_log10_k_ph5": nf.log10_k_anchor_ph5, "nf_ea": nf.ea_kj_mol, "nf_t_anchor_k": nf.t_anchor_c + 273.15,
        "entry_ea": entry.ea_entry_kj_mol, "aa_ea": entry.ea_dpo_aa_kj_mol,
        "pama_dpo_factor": entry.pama_dpo_factor, "pama_tdp_factor": entry.pama_tdp_factor,
        "adduct_k_ref": ps.THIOL_QUINONE_EQUILIBRIUM_USED_L_PER_MMOL, "adduct_dh": ps.STACK_DELTA_H_ADDUCT_KJ_MOL,
        "adduct_t_k": ps.STACK_MEASUREMENT_T_K,
        "amadori_factor": float(SCHIFF_AMADORI_SPLIT["ratio_amadori_over_schiff_pseudo_first_order"])
                          * float(SCHIFF_AMADORI_SPLIT["amine_loading_mmol_L_for_the_ratio"]),
        # the charge balance's water and activity tables
        "pkw_table_k": tuple(c + 273.15 for c, _ in ph.PKW_TABLE_C), "pkw_table": tuple(p for _, p in ph.PKW_TABLE_C),
        "rho_table_c": tuple(c for c, _ in ph._WATER_DENSITY_C), "rho_table": tuple(r for _, r in ph._WATER_DENSITY_C),
        "label_t_k": ph.LABEL_TEMPERATURE_K,
        # the entry's / NF split's fixed values when a coordinate does not carry them
        "fixed": {"mgo_share_nf": float(nf.mgo_share),
                  "log10_k_pent_amine": math.log10(entry.k_entry_l_per_mmol_min_100c_ph6),
                  "entry_ph_slope": float(entry.entry_ph_slope),
                  "log10_k_dpo_aa": (math.log10(entry.k_dpo_aa_per_min_100c_ph6)
                                     if entry.k_dpo_aa_per_min_100c_ph6 > 0 else -300.0),
                  "aa_ph_slope": float(entry.aa_ph_slope), "amine_return": float(entry.amine_return)},
    }
    names, pk, dh, dcp, zmax, weights, offset = _charge_balance_groups(
        sulfur.SULFUR_INDEX, drift, buffer_spec, switches["trunk_solutes"])
    stoich = np.asarray(sulfur.SULFUR_STOICHIOMETRY, dtype=float)
    return Lane(species=species, reactions=reactions, stoichiometry=stoich, reactant_index=r_idx,
                reactant_power=r_pow, k_ref=k_ref, ea_over_r=ea_over_r, inv_t_ref=inv_t_ref, cap_k=cap_k, kind=kind,
                coordinate_keys=tuple(coordinate_keys), offset_matrix=offset_matrix, index=index, const=const,
                group_names=names, group_pk25=pk, group_dh=dh, group_dcp=dcp, group_zmax=zmax,
                group_weights=weights, group_offset=offset,
                notes=(f"ddp_law={switches['ddp_law']}", f"trunk_solutes={switches['trunk_solutes']}"))
