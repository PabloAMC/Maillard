"""
src/kinetic_core_jax/sulfur_structure.py -- the sulfur lane's STRUCTURE as the JAX port reads it, and the check that
the NumPy lane has not moved under the port (plan step A2 + H3; docs/guides/A2_MODULE_CONTRACT.md sec. 6.2.2;
2026-11-05). Plain numpy and the standard library, NO JAX: it is imported by the core env's unit tier
(tests/unit/test_jax_sulfur_structure.py), so a change to ``src/kinetic_core/sulfur.py`` that the port does not follow
fails CI in the env every PR runs, not only in ``maillard-jax``.

WHAT IS DERIVED (read from the engine at build time, never retyped; ``derive``): the 62 species and their C / N / S
counts, the 132 reactions in ``sulfur.FULL_REACTIONS`` order, their stoichiometry and mass-action layout, every
parameter's k_ref / Ea / T_ref / temperature cap, the pH-factor kind of each reaction (``sulfur.REACTION_PH_FACTOR``),
the measured tables the hand-ported terms read (Zheng & Ho, the glutathione refit, the thiamine share, Stack's adduct
equilibrium, the Amadori split, every pK / dH / dCp of the charge balance, the pKw and water-density tables).

WHAT IS HAND-PORTED (``HAND_PORTED``: the reactions whose constant ``sulfur.sulfur_rate_constants_at`` does NOT
compute by the generic rule ``parameter.k_at(min(T, cap)) x ph_factor(kind, pH, T, I)``), and the engine logic the
port re-implements (``SOURCE_PINS``). Both are CHECKED against the engine (``structure_problems``):

1. the set of hand-ported reactions equals the set the engine computes outside the generic rule, found NUMERICALLY
   (the engine's own ``sulfur_rate_constants_at`` against the engine's own ``k_at`` and ``ph_factor`` on a probe of
   (T, pH, label pH, I)), and each hand-ported reaction reads exactly the parameter keys declared (each key's k_ref
   and Ea perturbed in turn);
2. every reaction of ``FULL_REACTIONS`` is generic or hand-ported, and every parameter key a generic reaction names is
   in ``engine.core_parameters("sulfur")``;
3. every engine Config key bound to ``sulfur`` / ``parameters_sulfur`` / ``ph_state`` is handled (``CONFIG_HANDLING``:
   ported at every value, held at its shipped value and refused otherwise, or not read by the kinetics), and the
   shipped defaults of the held keys are the ones the port was written against (``SHIPPED_HELD_DIGEST``);
4. the AST (docstrings and comments stripped) of every engine function whose logic the port re-implements hashes to
   its pin, so ANY logic edit of those functions fails until the JAX lane is updated and re-pinned in the same PR.

A failure names the file to update: ``UPDATE`` below. The numerical parity of the two lanes (rates, panel
trajectories, random parameter draws) is tests/scientific/test_jax_sulfur_parity.py, in the ``maillard-jax`` env.
"""
from __future__ import annotations

import ast
import hashlib

from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
NAME = "sulfur"
JAX_FILE = "src/kinetic_core_jax/sulfur_lane.py"
UPDATE = (f"update {JAX_FILE} (and its registry in src/kinetic_core_jax/sulfur_structure.py) in the same PR, then "
          "run tests/scientific/test_jax_sulfur_parity.py in the maillard-jax env")

# --------------------------------------------------------------------------- the hand-ported terms (shipped config)
#: reaction -> (the JAX term that computes it, the engine parameter keys its constant reads). "off" terms are zero at
#: the shipped configuration because their route switch is off (the port refuses the switch on: those routes stay
#: NumPy-only while M1b(i) / M1c reshape them).
HAND_PORTED: Mapping[str, Tuple[str, Tuple[str, ...]]] = {
    "r_cys_h2s": ("zheng_cysteine_thermolysis (label-scale pH table, 1/min)", ()),
    "r_gsh_h2s": ("zheng_glutathione_thermolysis (label-scale pH table, 1/min)", ()),
    "r_thi_h2s": ("thiamine_h2s_share x (k(r_thi_hmp) + k(r_thi_mesh)) (B50)", ("k_thi_hmp", "k_thi_mesh")),
    "ch_thioether_release": ("k(ch_thioether_mft) / K_adduct(T) (Stack 2018)", ("k_thioether",)),
    "ch_thioether_release_fft": ("k(ch_thioether_fft) / K_adduct(T) (Stack 2018)", ("k_thioether",)),
    "r_amadori": ("SCHIFF_AMADORI_SPLIT ratio x loading x k(r_schiff) (B1)", ("k_schiff",)),
    "r_ttca_form": ("off: TTCA_FORMATION and CYS_AMADORI_ROUTE disabled (B58 / B59b)", ()),
    "r_cama_form": ("off: CYS_AMADORI_ROUTE disabled (B59b)", ()),
    "r_cama_dpo": ("off: CYS_AMADORI_ROUTE disabled (B59b)", ()),
    "r_cama_tdp": ("off: CYS_AMADORI_ROUTE disabled (B59b)", ()),
    "r_dpo_ddp_h2s": ("off: STRECKER_SULFUR disabled (B59a)", ()),
    "r_nf_decay_oh": ("off: NF_HYDROXIDE_DECAY disabled (B61)", ()),
    "r_nf_decay_oh_mgo": ("off: NF_HYDROXIDE_DECAY disabled, mgo_share 0 (sub-model 3)", ()),
    "r_tdp_mgo": ("off: PENTOSE_MGO_ROUTES k_tdp_mgo 0 (sub-model 3)", ()),
    "r_pent_mgo": ("off: PENTOSE_MGO_ROUTES k_pent_mgo 0 (sub-model 3)", ()),
    "r_mgo_cys_strecker": ("off: MGO_CYSTEINE_STRECKER k 0", ()),
    "r_mgo_cys_h2s": ("off: MGO_CYSTEINE_STRECKER k 0", ()),
    "r_tdp_cys_trap": ("off: TDP_CYSTEINE_TRAP k 0 (B7)", ()),
    "r_odg_dpo": ("off: HEXOSE_C5_ENTRY k 0 (sub-model 4)", ()),
    "r_pent_amine": ("off: PENTOSE_AMINE_ENTRY disabled", ()),
    "r_pama_dpo": ("off: PENTOSE_AMINE_ENTRY disabled", ()),
    "r_pama_tdp": ("off: PENTOSE_AMINE_ENTRY disabled", ()),
    "r_pama_dpo_n": ("off: PENTOSE_AMINE_ENTRY disabled", ()),
    "r_pama_tdp_n": ("off: PENTOSE_AMINE_ENTRY disabled", ()),
    "r_dpo_aa": ("off: PENTOSE_AMINE_ENTRY k_dpo_aa 0", ()),
    "ch_oligomer_mft": ("off: oligomerisation_rate() is 0", ()),
    "ch_oligomer_fft": ("off: oligomerisation_rate() is 0", ()),
}
OFF = tuple(k for k, (term, _) in HAND_PORTED.items() if term.startswith("off:"))

#: Every engine Config key bound to the sulfur lane's modules, and how the port handles it.
#: "ported": every value the key can take is implemented; "held": the port implements the shipped value only and
#: refuses any other (UnsupportedConfiguration); "not_read": the kinetics never read it.
CONFIG_HANDLING: Mapping[str, Tuple[str, str]] = {
    "ph.coupling": ("ported", "both values: 'fixed_point_nodes' (shipped; the node fixed point) and 'inline' (the "
                              "charge balance in the right-hand side)"),
    "ph.node_grid": ("ported", "both values: 'log_and_linear' (shipped) and 'uniform'"),
    "ph.fixed_point": ("ported", "both values: 'relaxed_converged' (shipped) and 'fixed_sweeps'"),
    "ph.fixed_point.relaxation": ("ported", "numeric, read per build"),
    "ph.fixed_point.tolerance": ("ported", "numeric, read per build"),
    "ph.fixed_point.max_sweeps": ("ported", "numeric, read per build (a static loop bound)"),
    "ph.fixed_point.flag": ("ported", "numeric, read per build (the convergence flag)"),
    "sulfur.profile_node_step_k": ("ported", "numeric: the profile's node spacing, read per pot"),
    "ph.b51": ("held", "B51 all on; the rivals (any correction off) stay NumPy-only"),
    "ph.trunk_solutes_in_balance": ("held", "False (sub-model 1's rival, off)"),
    "ph.trunk_terms": ("held", "False (sub-model 1's rival, off)"),
    "ph.trunk_exponent_override": ("held", "empty"),
    "sulfur.ttca_formation": ("held", "disabled (B58)"),
    "sulfur.nf_hydroxide_decay": ("held", "disabled (B61)"),
    "sulfur.pentose_mgo_routes": ("held", "both constants 0 (sub-model 3)"),
    "sulfur.pentose_amine_entry": ("held", "disabled, k_dpo_aa 0"),
    "sulfur.ddp_formation_ph": ("held", "law 'base' (the generic base factor)"),
    "sulfur.mgo_cysteine_strecker": ("held", "k 0"),
    "sulfur.tdp_cysteine_trap": ("held", "k 0"),
    "sulfur.hexose_c5_entry": ("held", "k 0"),
    "sulfur.strecker_sulfur": ("held", "disabled (B59a)"),
    "sulfur.cys_amadori_route": ("held", "disabled (B59b)"),
    "sulfur.tdp_own_decay": ("held", "disabled (B7a R-T)"),
    "sulfur.furfural_barriers": ("held", "disabled (B62); moves parameter values through core_parameters"),
    "sulfur.single_temperature_barriers": ("held", "disabled (B63); moves parameter values through core_parameters"),
    "ph.ionic_strength.sweeps": ("not_read", "the engine's I fixed-point iteration count: the port solves (pH, I) "
                                             "jointly by Newton to machine precision, the limit those sweeps approach"),
    "ph.ionic_strength.tol_m": ("not_read", "the engine's I fixed-point tolerance (as above)"),
    "sulfur.flux_budget": ("not_read", "the share rows' flux-budget post-processing, not the kinetics"),
    "fit.dimer_encounter_ceiling": ("not_read", "read by fit code only, never by the kinetics"),
    "fit.mft_decay_prior": ("not_read", "read by fit code only, never by the kinetics"),
}
HELD = tuple(k for k, (how, _) in CONFIG_HANDLING.items() if how == "held")
LANE_MODULES = ("sulfur", "parameters_sulfur", "ph_state")

#: sha256 (16 hex) of the shipped defaults of the HELD keys, as the port was written against them (2026-11-05)
SHIPPED_HELD_DIGEST = "e40256046969dff0"

#: (file, qualified name) -> sha256 (16 hex) of the AST without docstrings: the engine logic the port re-implements
SOURCE_PINS: Mapping[Tuple[str, str], str] = {
    ('src/kinetic_core/sulfur.py', 'sulfur_rate_constants_at'): '29db02f032049288',
    ('src/kinetic_core/sulfur.py', 'integrate_sulfur'): 'b856bb5a3fd80c82',
    ('src/kinetic_core/sulfur.py', '_integrate_inline_ph'): '3593b6f152b1dabb',
    ('src/kinetic_core/sulfur.py', 'ph_node_times'): '2c9d370551e51d33',
    ('src/kinetic_core/sulfur.py', 'profile_node_times'): '4fafa180d7d9e132',
    ('src/kinetic_core/sulfur.py', 'rate_vector_at'): '5574d5fec2ec287e',
    ('src/kinetic_core/sulfur.py', '_solve_on_time_nodes'): '3fac39ddfd47825e',
    ('src/kinetic_core/sulfur.py', '_rates_vectorised'): '75be8610ed97ed03',
    ('src/kinetic_core/sulfur.py', 'sulfur_stoichiometric_matrix'): 'b7568d16378e5d6b',
    ('src/kinetic_core/parameters_sulfur.py', 'SulfurParameter.k_at'): 'c693066311e4e7a1',
    ('src/kinetic_core/parameters_sulfur.py', 'cysteine_thermolysis_arrhenius'): '457463d4ead1d33d',
    ('src/kinetic_core/parameters_sulfur.py', 'cysteine_thermolysis_k'): 'c773761729026629',
    ('src/kinetic_core/parameters_sulfur.py', 'glutathione_thermolysis_k'): '74902265bd11db19',
    ('src/kinetic_core/parameters_sulfur.py', 'thiamine_h2s_share'): 'e3726b8b2acd9939',
    ('src/kinetic_core/parameters_sulfur.py', 'h2s_pk1_suleimenov_seward'): 'df1c9d6cfb4b0ff0',
    ('src/kinetic_core/parameters_sulfur.py', 'pka_at'): 'f4e23bb454b0e3ce',
    ('src/kinetic_core/parameters_sulfur.py', 'neutral_h2s_fraction'): 'd2cf0a127435f88c',
    ('src/kinetic_core/parameters_sulfur.py', 'hydrosulfide_fraction'): 'fec85dcf52e69e79',
    ('src/kinetic_core/parameters_sulfur.py', 'free_amine_fraction'): '7a514097386f0110',
    ('src/kinetic_core/parameters_sulfur.py', 'thiolate_fraction'): 'f0d479edfef69669',
    ('src/kinetic_core/parameters_sulfur.py', 'ph_factor'): 'feebe87205943f86',
    ('src/kinetic_core/parameters_sulfur.py', 'thiol_adduct_equilibrium_l_per_mmol'): '94032cde1a095c0a',
    ('src/kinetic_core/parameters_sulfur.py', 'oligomerisation_rate'): '75dd085efa2ebd10',
    ('src/kinetic_core/parameters.py', 'KineticParameter.k_at'): 'd3ccfb47a61953a0',
    ('src/kinetic_core/ph_state.py', 'pk_at'): 'c8809b2fc2f9c7b2',
    ('src/kinetic_core/ph_state.py', 'davies_a'): 'a232e97423ac08d0',
    ('src/kinetic_core/ph_state.py', 'davies_term'): 'dd38839a0b6736b3',
    ('src/kinetic_core/ph_state.py', '_interp'): 'ee387cfa74d9b1a8',
    ('src/kinetic_core/ph_state.py', 'pkw_at'): 'e18f5b2301a9cc07',
    ('src/kinetic_core/ph_state.py', 'TitratableGroup.pkas_at'): 'a0d7d239f54d5b1d',
    ('src/kinetic_core/ph_state.py', 'TitratableGroup._weights'): 'ea0ed89484096838',
    ('src/kinetic_core/ph_state.py', 'TitratableGroup.average_charge'): 'eb318e6c9324d406',
    ('src/kinetic_core/ph_state.py', 'TitratableGroup.ionic_strength_contribution'): '67409835d6df60e8',
    ('src/kinetic_core/ph_state.py', 'BufferSpec.groups'): 'c4526eeafd11a1ac',
    ('src/kinetic_core/ph_state.py', 'PhDrift.arp_group'): '9a98e3c5c0ccd2f5',
    ('src/kinetic_core/ph_state.py', 'cysteine_group'): 'b9d8ad380c89c68f',
    ('src/kinetic_core/ph_state.py', 'titratable_inventory'): '998c74d159d2950a',
    ('src/kinetic_core/ph_state.py', '_charge_residual'): '25fd1a380a396e57',
    ('src/kinetic_core/ph_state.py', 'ionic_strength_of'): 'f08e695bd92fc433',
    ('src/kinetic_core/ph_state.py', 'solve_ph_and_ionic_strength'): 'b759286fbe32928b',
    ('src/kinetic_core/ph_state.py', 'strong_ion_difference'): 'c84701bd17961d7f',
    ('src/kinetic_core/ph_state.py', 'ph_of_state'): '1561f4b9067f1523',
    ('src/kinetic_core/ph_state.py', 'ph_and_ionic_strength_of_state'): 'a4f7d5ae751936f9',
}


# --------------------------------------------------------------------------- source pins
def _strip_docstrings(node: ast.AST) -> ast.AST:
    for sub in ast.walk(node):
        body = getattr(sub, "body", None)
        if isinstance(sub, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef, ast.Module)) and body:
            first = body[0]
            if isinstance(first, ast.Expr) and isinstance(getattr(first, "value", None), ast.Constant) \
                    and isinstance(first.value.value, str):
                sub.body = body[1:] or [ast.Pass()]
    return node


def source_hash(path: str, qualname: str) -> str:
    """sha256 (16 hex) of one function's (or method's) AST, docstrings removed, comments ignored."""
    tree = ast.parse((ROOT / path).read_text(encoding="utf-8"))
    node: Any = tree
    for part in qualname.split("."):
        found = [n for n in getattr(node, "body", []) if isinstance(n, (ast.FunctionDef, ast.ClassDef))
                 and n.name == part]
        if not found:
            raise KeyError(f"{path}: no {qualname!r}")
        node = found[-1]
    node = _strip_docstrings(ast.parse(ast.unparse(node)))
    return hashlib.sha256(ast.dump(node, annotate_fields=True, include_attributes=False).encode()).hexdigest()[:16]


#: the engine functions whose LOGIC the port re-implements (everything else it reads as data)
PINNED_FUNCTIONS: Tuple[Tuple[str, str], ...] = (
    ("src/kinetic_core/sulfur.py", "sulfur_rate_constants_at"),
    ("src/kinetic_core/sulfur.py", "integrate_sulfur"),
    ("src/kinetic_core/sulfur.py", "_integrate_inline_ph"),
    ("src/kinetic_core/sulfur.py", "ph_node_times"),
    ("src/kinetic_core/sulfur.py", "profile_node_times"),
    ("src/kinetic_core/sulfur.py", "rate_vector_at"),
    ("src/kinetic_core/sulfur.py", "_solve_on_time_nodes"),
    ("src/kinetic_core/sulfur.py", "_rates_vectorised"),
    ("src/kinetic_core/sulfur.py", "sulfur_stoichiometric_matrix"),
    ("src/kinetic_core/parameters_sulfur.py", "SulfurParameter.k_at"),
    ("src/kinetic_core/parameters_sulfur.py", "cysteine_thermolysis_arrhenius"),
    ("src/kinetic_core/parameters_sulfur.py", "cysteine_thermolysis_k"),
    ("src/kinetic_core/parameters_sulfur.py", "glutathione_thermolysis_k"),
    ("src/kinetic_core/parameters_sulfur.py", "thiamine_h2s_share"),
    ("src/kinetic_core/parameters_sulfur.py", "h2s_pk1_suleimenov_seward"),
    ("src/kinetic_core/parameters_sulfur.py", "pka_at"),
    ("src/kinetic_core/parameters_sulfur.py", "neutral_h2s_fraction"),
    ("src/kinetic_core/parameters_sulfur.py", "hydrosulfide_fraction"),
    ("src/kinetic_core/parameters_sulfur.py", "free_amine_fraction"),
    ("src/kinetic_core/parameters_sulfur.py", "thiolate_fraction"),
    ("src/kinetic_core/parameters_sulfur.py", "ph_factor"),
    ("src/kinetic_core/parameters_sulfur.py", "thiol_adduct_equilibrium_l_per_mmol"),
    ("src/kinetic_core/parameters_sulfur.py", "oligomerisation_rate"),
    ("src/kinetic_core/parameters.py", "KineticParameter.k_at"),
    ("src/kinetic_core/ph_state.py", "pk_at"),
    ("src/kinetic_core/ph_state.py", "davies_a"),
    ("src/kinetic_core/ph_state.py", "davies_term"),
    ("src/kinetic_core/ph_state.py", "_interp"),
    ("src/kinetic_core/ph_state.py", "pkw_at"),
    ("src/kinetic_core/ph_state.py", "TitratableGroup.pkas_at"),
    ("src/kinetic_core/ph_state.py", "TitratableGroup._weights"),
    ("src/kinetic_core/ph_state.py", "TitratableGroup.average_charge"),
    ("src/kinetic_core/ph_state.py", "TitratableGroup.ionic_strength_contribution"),
    ("src/kinetic_core/ph_state.py", "BufferSpec.groups"),
    ("src/kinetic_core/ph_state.py", "PhDrift.arp_group"),
    ("src/kinetic_core/ph_state.py", "cysteine_group"),
    ("src/kinetic_core/ph_state.py", "titratable_inventory"),
    ("src/kinetic_core/ph_state.py", "_charge_residual"),
    ("src/kinetic_core/ph_state.py", "ionic_strength_of"),
    ("src/kinetic_core/ph_state.py", "solve_ph_and_ionic_strength"),
    ("src/kinetic_core/ph_state.py", "strong_ion_difference"),
    ("src/kinetic_core/ph_state.py", "ph_of_state"),
    ("src/kinetic_core/ph_state.py", "ph_and_ionic_strength_of_state"),
)


def current_pins() -> Dict[Tuple[str, str], str]:
    return {key: source_hash(*key) for key in PINNED_FUNCTIONS}


def held_digest(defaults: Mapping[str, Any]) -> str:
    text = "\n".join(f"{k}={defaults[k]!r}" for k in sorted(HELD))
    return hashlib.sha256(text.encode()).hexdigest()[:16]


# --------------------------------------------------------------------------- derivation from the engine
def _engine():
    from src.kinetic_core import config, engine, ph_state, sulfur
    from src.kinetic_core import parameters_sulfur as ps

    return config, engine, ph_state, sulfur, ps


def derive() -> Dict[str, Any]:
    """The lane's static structure, read from the engine: species, elements, reactions (FULL_REACTIONS order),
    stoichiometry and reactant orders keyed by name, each reaction's parameter key and pH-factor kind, and the
    parameter keys the generic reactions name."""
    _, _, _, sulfur, _ = _engine()
    from src.kinetic_core.species_sulfur import SULFUR_BY_KEY, SULFUR_STATE_KEYS

    reactions = tuple(sulfur.FULL_REACTIONS)
    stoich: Dict[str, Dict[str, float]] = {}
    reactants: Dict[str, Dict[str, int]] = {}
    for r in reactions:
        net: Dict[str, float] = {}
        for s, c in r.reactants.items():
            net[s] = net.get(s, 0.0) - float(c)
        for s, c in r.products.items():
            net[s] = net.get(s, 0.0) + float(c)
        stoich[r.key] = {s: c for s, c in net.items() if c != 0.0}
        reactants[r.key] = {s: int(c) for s, c in r.reactants.items()}
    generic = [r for r in reactions if r.key not in HAND_PORTED]
    return {
        "species": tuple(SULFUR_STATE_KEYS),
        "elements": {k: (int(SULFUR_BY_KEY[k].carbon), int(SULFUR_BY_KEY[k].nitrogen), int(SULFUR_BY_KEY[k].sulfur))
                     for k in SULFUR_STATE_KEYS},
        "reactions": tuple(r.key for r in reactions),
        "parameter_key": {r.key: r.parameter_key for r in reactions},
        "statistical_factor": {r.key: float(r.statistical_factor) for r in reactions},
        "kind": {r.key: sulfur.REACTION_PH_FACTOR.get(r.key, "") for r in reactions},
        "stoichiometry": stoich, "reactants": reactants,
        "generic_parameter_keys": tuple(dict.fromkeys(r.parameter_key for r in generic)),
    }


#: the probe on which the engine's rate function is compared with the generic rule: (T_C, in-situ pH, label pH, I)
PROBE = ((60.0, 5.0, 5.0, 0.0), (100.0, 3.4, 3.9, 0.02), (120.0, 6.3, 6.6, 0.4), (145.0, 7.6, 8.1, 0.9),
         (180.0, 4.6, 5.6, 0.15), (25.0, 9.3, 9.6, 0.05))


def generic_constant(parameters: Mapping[str, Any], reaction, temperature_k: float, ph: float,
                     ionic_m: float) -> float:
    """The GENERIC rule, from the engine's own functions: k_at at the capped temperature times the step's pH factor."""
    _, _, _, sulfur, ps = _engine()
    from src.kinetic_core.parameters_furanic import TEMPERATURE_CAP_K

    p = parameters[reaction.parameter_key]
    cap = TEMPERATURE_CAP_K.get(reaction.parameter_key)
    t = temperature_k if cap is None else min(temperature_k, float(cap))
    return p.k_at(t) * ps.ph_factor(sulfur.REACTION_PH_FACTOR.get(reaction.key, ""), ph, temperature_k, ionic_m)


def _engine_constants(parameters, probe=PROBE) -> np.ndarray:
    _, _, _, sulfur, _ = _engine()
    rows = []
    for t_c, ph, label, ionic in probe:
        k = sulfur.sulfur_rate_constants_at(parameters, t_c + 273.15, ph, ph_label=label, ionic_strength_m=ionic)
        rows.append([float(k[r.key]) for r in sulfur.FULL_REACTIONS])
    return np.asarray(rows)


def classify(parameters=None) -> Dict[str, Any]:
    """Which reactions the engine computes outside the generic rule (numerically), and which parameter keys each of
    those reads (each key's k_ref x 1.7 and, where it has one, Ea + 7 kJ/mol, in turn)."""
    _, engine, _, sulfur, _ = _engine()
    import dataclasses

    parameters = dict(engine.core_parameters(engine.SULFUR) if parameters is None else parameters)
    base = _engine_constants(parameters)
    custom = set()
    for j, r in enumerate(sulfur.FULL_REACTIONS):
        if r.parameter_key is None or r.parameter_key not in parameters:
            custom.add(r.key)
            continue
        gen = np.array([generic_constant(parameters, r, t + 273.15, ph, i) for t, ph, _, i in PROBE])
        if not np.allclose(base[:, j], gen, rtol=1e-13, atol=0.0):
            custom.add(r.key)
    reads: Dict[str, set] = {k: set() for k in custom}
    index = {r.key: j for j, r in enumerate(sulfur.FULL_REACTIONS)}
    for key, p in parameters.items():
        moved = []
        if p.k_ref is not None:
            moved.append(dataclasses.replace(p, k_ref=(float(p.k_ref) * 1.7 if float(p.k_ref) != 0.0 else 1e-3)))
        if getattr(p, "ea_kj_mol", None) is not None:
            moved.append(dataclasses.replace(p, ea_kj_mol=float(p.ea_kj_mol) + 7.0))
        for q in moved:
            k = _engine_constants({**parameters, key: q})
            for name in custom:
                if not np.allclose(k[:, index[name]], base[:, index[name]], rtol=1e-14, atol=0.0):
                    reads[name].add(key)
    zero = {k for k in custom if not np.any(base[:, index[k]])}
    return {"custom": custom, "reads": {k: tuple(sorted(v)) for k, v in reads.items()}, "zero": zero}


def lane_config_keys() -> Tuple[str, ...]:
    config, *_ = _engine()
    return tuple(k.name for k in config.KEYS if {m for m, _ in k.bindings} & set(LANE_MODULES))


def structure_problems(*, numeric: bool = True) -> List[str]:
    """Everything that says the NumPy lane moved under the port, each with ``UPDATE``. Empty when in step."""
    config, engine, _, sulfur, _ = _engine()
    out: List[str] = []
    keys = {r.key for r in sulfur.FULL_REACTIONS}
    stray = sorted(set(HAND_PORTED) - keys)
    if stray:
        out.append(f"HAND_PORTED names reactions sulfur.FULL_REACTIONS no longer has: {stray}")
    parameters = engine.core_parameters(engine.SULFUR)
    missing = sorted({r.parameter_key for r in sulfur.FULL_REACTIONS if r.key not in HAND_PORTED}
                     - set(parameters))
    if missing:
        out.append(f"generic reactions name parameters engine.core_parameters('sulfur') lacks: {missing}")
    factors = sorted(r.key for r in sulfur.FULL_REACTIONS if float(r.statistical_factor) != 1.0)
    if factors:
        out.append(f"reactions with a statistical factor (sulfur_rate_constants_at ignores it; the port would have "
                   f"to choose): {factors}")
    if numeric:
        found = classify(parameters)
        new = sorted(found["custom"] - set(HAND_PORTED))
        if new:
            out.append(f"sulfur.sulfur_rate_constants_at computes {new} outside the generic rule "
                       "(k_at x ph_factor) and the JAX lane has no hand-ported term for them")
        gone = sorted(set(HAND_PORTED) - found["custom"])
        if gone:
            out.append(f"HAND_PORTED lists {gone}, which the engine now computes by the generic rule")
        for name in sorted(found["custom"] & set(HAND_PORTED)):
            if tuple(sorted(HAND_PORTED[name][1])) != found["reads"][name]:
                out.append(f"{name} reads parameters {found['reads'][name]} in the engine, "
                           f"{tuple(sorted(HAND_PORTED[name][1]))} in the JAX lane")
        off_nonzero = sorted(set(OFF) - found["zero"])
        if off_nonzero:
            out.append(f"reactions the JAX lane holds at zero are non-zero in the engine: {off_nonzero}")
    bound = set(lane_config_keys())
    unhandled = sorted(bound - set(CONFIG_HANDLING))
    if unhandled:
        out.append(f"engine Config keys bound to the sulfur lane that the JAX lane neither ports nor refuses: "
                   f"{unhandled}")
    gone_keys = sorted(set(CONFIG_HANDLING) - bound)
    if gone_keys:
        out.append(f"CONFIG_HANDLING lists keys the engine no longer binds to the lane: {gone_keys}")
    shipped = config.Config.shipped()
    digest = held_digest({k: shipped[k] for k in HELD if k in bound})
    if digest != SHIPPED_HELD_DIGEST:
        out.append(f"the shipped defaults of the held keys moved (digest {digest}, port written against "
                   f"{SHIPPED_HELD_DIGEST}): {HELD}")
    for key, pin in current_pins().items():
        if SOURCE_PINS.get(key) != pin:
            out.append(f"{key[0]}::{key[1]} changed (AST {pin}, pinned {SOURCE_PINS.get(key)})")
    return [f"{p} -- {UPDATE}" for p in out]


def nan_safe(value) -> float:
    return float("nan") if value is None else float(value)
