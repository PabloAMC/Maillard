"""
src/kinetic_core/config.py -- THE ENGINE'S SWITCHES AS ONE FROZEN VALUE (A3 phase 1, 2026-11-01).

WHY THIS EXISTS
===============
The engine's structural and solver choices are module globals (``sulfur.PH_COUPLING``, ``ph_state.B51``, thirteen
``parameters_sulfur`` route dataclasses, the lipid / dicarbonyl / methionine / proline / panel / uncertainty flags).
Scripts set them by assignment, about 97 times in 39 files, and a job's configuration exists only as the state of the
process that runs it. ``parallel.map_ordered`` is correct only because it FORKS: a spawn-started or long-lived worker
would compute with whatever its globals hold, without an error (docs/guides/A2_MODULE_CONTRACT.md sec. 1).

This module gives that state a name and a value. ``KEYS`` declares one ``ConfigKey`` per switch (namespaced, with the
module attribute(s) it stands for, its shipped default, its phase and its basis); a ``Config`` is a frozen, hashable
assignment of a value to EVERY key, which can travel with a job and be installed by whichever process runs it.

PHASE 1 IS A BRIDGE, NOT A REWRITE
==================================
The engine still reads its globals; nothing in ``src/`` reads a ``Config``. ``install`` / ``applied`` set those
globals from a Config, ``Config.from_globals`` reads them back, so a Config reproduces exactly what assigning the
globals by hand did, bit for bit. Two consequences, stated so they are not mistaken for features:

* a ``"build"`` key is also read at IMPORT by some table (``parameters_methionine.METHIONINE_PARAMETERS`` reads
  ``METHIONINE_SHIPPED``; ``uncertainty.CORE_PRIORS`` reads ``B55_SINGLE_TEMPERATURE_EA_PRIOR``). Installing it at
  run time changes only the call-time readers, exactly as setting the global did;
* a switch re-imported by name into another module (``lipid`` imports ``DECADIENAL_BREAKAGE`` from
  ``parameters_lipid``) has one binding per module, and the engine reads different ones in different places. A key
  lists every binding and ``install`` sets them all; ``from_globals`` refuses a state in which they disagree, because
  no single value describes it.

``shipped()`` is built from the declarations below, not from the modules' current state, so it is the shipped
configuration whatever a script has set before importing this module. ``tests/unit/test_engine_config.py`` checks it
against the modules as first imported, checks that every switch in the engine has a key (a new switch without one
fails), and that a SPAWN-started worker given a Config reproduces a forked one with every switch flipped.

PHASE 2 (2026-11-01): THE OWNER'S DECISIONS AND THE REVIEW OF PHASE 1
=====================================================================
* ``"build"`` keys are READ-ONLY FACTS. Their tables are built at import, so installing one at run time would change
  the call-time readers only and leave a half-switched engine. ``install`` / ``applied`` refuse a Config whose build
  key differs from the shipped default, naming the key, and never set a build binding. A Config may still HOLD a
  non-default build value (``from_globals`` reports what a process holds); it cannot be installed.
* A key validates its value on construction: the type of the shipped default (any number stands, and is stored, as
  a float where the default is a float; a bool is never a number), ``choices``, a scalar ``band``, and FIELD-LEVEL
  ``bands`` ``{field: (lo, hi)}`` for the dataclass fields that fits and feasibility draws sample (``_BANDS``, with
  its basis). The fields of a dataclass value are type-checked against the shipped instance's.
* ``Config.digest()`` is a sha256 over a canonical serialisation (numbers as exact float hex, dataclasses with their
  qualified type and every field), stable across processes, unlike ``hash()``.
* ``update(changes)`` is what a script calls instead of assigning a global: the process's configuration with some
  keys replaced, validated, installed on EVERY binding (a script that assigned one binding of a re-imported switch
  left the other behind: build_choice_robustness's Q10, fixed 2026-11-01), and returned, so the script hands the
  same value to its jobs (``parallel.map_ordered(..., config=cfg)``).
"""
from __future__ import annotations

import contextlib
import dataclasses
import hashlib
import importlib
import json
import math
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Dict, Iterator, Literal, Mapping, Optional, Tuple

from . import parameters_dicarbonyl as _pd
from . import parameters_lipid as _pl
from . import parameters_sulfur as _ps
from . import ph_balance_rivals as _pbr
from . import ph_state as _ph

Phase = Literal["build", "rate", "solver"]
_PACKAGE = "src.kinetic_core"


def _is_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _same_kind(shipped: Any, value: Any) -> bool:
    """``value`` may stand where ``shipped`` stands: the same type, except that any number (never a bool) stands for
    a float, any tuple for a tuple, and None on either side (an Optional field)."""
    if shipped is None or value is None:
        return True
    if isinstance(shipped, bool) or isinstance(value, bool):
        return type(shipped) is type(value)
    if isinstance(shipped, float):
        return _is_number(value)
    if isinstance(shipped, int):
        return isinstance(value, int)
    if isinstance(shipped, tuple):
        return isinstance(value, tuple)
    return type(value) is type(shipped)


def _in_band(value: Any, band: Tuple[float, float]) -> bool:
    return _is_number(value) and not math.isnan(value) and band[0] <= value <= band[1]


@dataclass(frozen=True)
class ConfigKey:
    """One switch. ``bindings`` are ``(module, attribute)`` pairs, module relative to ``src.kinetic_core`` ("" for the
    package itself); the first is the one the engine's kinetics read. ``mapping`` keys hold a dict, kept in a Config
    as a tuple of its items (hashable, order kept) and installed as a fresh dict. ``band`` bounds a scalar value,
    ``bands`` named FIELDS of a dataclass value (``{field: (lo, hi)}``, inclusive)."""

    name: str
    bindings: Tuple[Tuple[str, str], ...]
    default: Any
    phase: Phase
    basis: str
    band: Optional[Tuple[Any, Any]] = None
    choices: Optional[Tuple[Any, ...]] = None
    mapping: bool = False
    bands: Mapping[str, Tuple[float, float]] = field(default_factory=lambda: MappingProxyType({}), hash=False)

    def freeze(self, value: Any) -> Any:
        if self.mapping:
            return tuple((k, v) for k, v in dict(value).items())
        if isinstance(self.default, float) and _is_number(value):
            return float(value)
        return value

    def thaw(self, value: Any) -> Any:
        return dict(value) if self.mapping else value

    def validate(self, value: Any) -> None:
        """Raise unless ``value`` (frozen) is one this key may take: its type, choices, band and field bands."""
        if self.mapping:
            if not isinstance(value, tuple) or not all(
                    isinstance(item, tuple) and len(item) == 2 and isinstance(item[0], str) for item in value):
                raise TypeError(f"{self.name}: a mapping from names, got {value!r}")
        elif not _same_kind(self.default, value):
            raise TypeError(f"{self.name}: a {type(self.default).__name__}, got {type(value).__name__} {value!r}")
        if self.choices is not None and value not in self.choices:
            raise ValueError(f"{self.name} = {value!r}: not one of {self.choices}")
        if self.band is not None and not _in_band(value, self.band):
            raise ValueError(f"{self.name} = {value!r}: outside its band {self.band}")
        if dataclasses.is_dataclass(self.default):
            for f in dataclasses.fields(self.default):
                got, shipped = getattr(value, f.name), getattr(self.default, f.name)
                if not _same_kind(shipped, got):
                    raise TypeError(f"{self.name}.{f.name}: a {type(shipped).__name__}, got {got!r}")
        for name, band in self.bands.items():
            got = getattr(value, name)
            if not _in_band(got, band):
                raise ValueError(f"{self.name}.{name} = {got!r}: outside its band {band} (config._BANDS)")


def _k(name: str, module: str, attr: str, default: Any, phase: Phase, basis: str, *, also: Tuple[str, ...] = (),
       band=None, choices=None, mapping: bool = False) -> ConfigKey:
    bindings = ((module, attr),) + tuple((m, attr) for m in also)
    return ConfigKey(name, bindings, default, phase, basis, band=band, choices=choices, mapping=mapping,
                     bands=MappingProxyType(dict(_BANDS.get(name, {}))))


#: FIELD BANDS (``ConfigKey.bands``, the owner's decision 6 of 2026-11-01): for every dataclass field a fit or a
#: feasibility draw samples, the union of the boxes the scripts sample, in the field's own units (a log10 coordinate
#: in [a, b] gives the band [0, 10**b]: 0 is the route off). Basis: the joint fit's coordinates
#: (fit_joint_deoxyosone_thiol.py: ENTRY_COORDINATES, AA_PK_COORDINATE, SPLIT_COORDINATES, TDP_MGO_COORDINATE,
#: TDP_CYS_TRAP_COORDINATE, MGO_STRECKER_BOUNDS, DDP_SLOPE_COORDINATE), sub-model 3's arms
#: (fit_deoxyosone_submodel.coordinates: k_tdp_mgo to 1, k_pent_mgo to 0.1, both Ea 40 - 160) and the feasibility
#: draws (feasibility_deoxyosone.BANDS, TDP_EXIT_BAND, and its +-1 dex about the fitted log10 k_dpo_aa, so 10**1.5);
#: shares are [0, 1]; the Q10 is Q10Assumption's own lo / hi (CR-02's three alternatives). A job outside a band is
#: refused: widen the band here, with its basis, before sampling outside it.
_SHARE = (0.0, 1.0)
_BANDS: Mapping[str, Mapping[str, Tuple[float, float]]] = {
    "sulfur.nf_hydroxide_decay": {"mgo_share": _SHARE},
    "sulfur.pentose_mgo_routes": {"k_tdp_mgo_100c_ph5": (0.0, 1.0), "ea_tdp_mgo_kj_mol": (40.0, 160.0),
                                  "k_pent_mgo_100c_ph5": (0.0, 0.1), "ea_pent_mgo_kj_mol": (40.0, 160.0)},
    "sulfur.pentose_amine_entry": {"k_entry_l_per_mmol_min_100c_ph6": (0.0, 1e-2), "entry_ph_slope": (-0.5, 1.5),
                                   "k_dpo_aa_per_min_100c_ph6": (0.0, 10.0 ** 1.5), "aa_ph_slope": (-0.5, 1.5),
                                   "aa_pk": (5.0, 7.5), "amine_return": _SHARE},
    "sulfur.ddp_formation_ph": {"slope": (-1.0, 2.0)},
    "sulfur.mgo_cysteine_strecker": {"k_l_per_mmol_min_100c_ph5": (0.0, 1.0)},
    "sulfur.tdp_cysteine_trap": {"k_l_per_mmol_min_100c_ph5": (0.0, 1.0)},
    "sulfur.strecker_sulfur": {"h2s_share": _SHARE},
    # M1 (2026-11-04; kinetic_core_m1_supply_loss_combinations_prereg.md secs. 2.2, 3): ea_form FITTED in 40 - 135
    # kJ/mol (the sugar_trunk class band, as stage 1), ea_decay drawn in the same band by the feasibility probe; the
    # route calibration samples both in FORMATION_EA_BOUNDS_BY_ROUTE["sugar_trunk"] (40, 135) too. k_form: the union of
    # the B59b anchor's search (log10 -6 .. -1), the route calibration's +-3 dex about 3.29e-3 and M1's +-0.44 dex.
    "sulfur.cys_amadori_route": {"dpo_share": _SHARE, "split_ph_c": (5.0, 10.0), "ea_form_kj_mol": (40.0, 135.0),
                                 "ea_decay_kj_mol": (40.0, 135.0), "k_form_ref_per_min": (1e-6, 3.3)},
    # K: the route calibration's box (0.019 .. its k_upper 0.6 sensitivity, run log 2); M1's probe draws 0.045 - 0.075
    "sulfur.ttca_formation": {"k_l_per_mmol": (0.019, 0.6)},
    "sulfur.tdp_own_decay": {"k_ref_per_min": (0.0, 0.1), "ea_kj_mol": (30.0, 200.0)},
    "lipid.q10_assumption": {"default": (1.6, 3.3)},
    # M1b part (i) (kinetic_core_m1b_ph_balance_pentose_steps_prereg.md sec. 3.1): fractions [0, 1]; k_lag from Zhai
    # 2023's 0.0097 /min to the band's 10 /min, widened to 1e6 for V2's fast-end reproduction; Ea the carbohydrate
    # class range (B55: 30 - 237 kJ/mol); k_iso [0.025, 0.24] /min at 120 C, widened to [1e-3, 10] for the checks
    "ph.amine_survival": {"phi_cys": _SHARE, "phi_arp": _SHARE},
    "ph.acid_lag": {"f_ret": _SHARE, "k_lag_per_min_120c": (1e-3, 1e6), "ea_lag_kj_mol": (0.0, 240.0)},
    "ph.acid_supply": {"y5": _SHARE, "y6": _SHARE, "y_hi": _SHARE, "k_iso_per_min_120c": (1e-3, 10.0),
                       "ea_iso_kj_mol": (0.0, 240.0)},
}


#: Every module-level switch of the engine. The default is the SHIPPED value, written here and checked against the
#: module as first imported (tests/unit/test_engine_config.py); a change of a shipped default changes both.
KEYS: Tuple[ConfigKey, ...] = (
    # --- the dynamic pH: coupling, nodes, fixed point, charge balance (sulfur.py, ph_state.py)
    _k("ph.coupling", "sulfur", "PH_COUPLING", "fixed_point_nodes", "solver",
       "C5 (2026-10-28); B64 runs on inline (kinetic_core_b64_prereg.md)", choices=("fixed_point_nodes", "inline")),
    _k("ph.node_grid", "sulfur", "PH_NODE_GRID", "log_and_linear", "solver",
       "B60 (kinetic_core_b60_prereg.md); frozen waves set uniform", choices=("uniform", "log_and_linear")),
    _k("ph.fixed_point", "sulfur", "PH_FIXED_POINT", "relaxed_converged", "solver",
       "B60 amendment 1; frozen waves set fixed_sweeps", choices=("fixed_sweeps", "relaxed_converged")),
    _k("ph.fixed_point.relaxation", "sulfur", "PH_FIXED_POINT_RELAXATION", 0.7, "solver",
       "B60 amendment 1: omega 1.0 / 0.7 / 0.5 / 0.3 reach the same answer", band=(0.3, 1.0)),
    _k("ph.fixed_point.tolerance", "sulfur", "PH_FIXED_POINT_TOLERANCE", 1e-3, "solver",
       "B60 amendment 1: the residual floors at 3e-5 - 6e-4 pH from the ODE tolerances"),
    _k("ph.fixed_point.max_sweeps", "sulfur", "PH_FIXED_POINT_MAX_SWEEPS", 40, "solver", "B60 amendment 1"),
    _k("ph.fixed_point.flag", "sulfur", "PH_FIXED_POINT_FLAG", 1e-2, "solver",
       "B60 amendment 1 / 2026-10-28: a run above it is marked unconverged and warns"),
    _k("ph.b51", "ph_state", "B51", _ph.B51_ALL_ON, "rate", "B51 pH-state corrections; all off reproduces B50"),
    _k("ph.trunk_solutes_in_balance", "ph_state", "TRUNK_SOLUTES_IN_BALANCE", False, "rate",
       "sub-model 1 Amendment 1 (Gly / AMA / SB / PAMA titrate in the charge balance)"),
    _k("ph.ionic_strength.sweeps", "ph_state", "IONIC_STRENGTH_SWEEPS", 30, "solver", "ph_state ionic-strength loop"),
    _k("ph.ionic_strength.tol_m", "ph_state", "IONIC_STRENGTH_TOL_M", 1e-9, "solver", "ph_state ionic-strength loop"),
    _k("ph.trunk_terms", "sulfur", "TRUNK_PH_TERMS", False, "rate",
       "sub-model 1 (2026-10-29): Martins & van Boekel 2003 Part II Table 3 slopes in the sulfur lane"),
    _k("ph.trunk_exponent_override", "sulfur", "TRUNK_PH_EXPONENT_OVERRIDE", (), "rate",
       "sub-model 1 fit slopes by parameter key (fit_ph_acid_submodel.py); empty otherwise", mapping=True),
    # --- M1b part (i)'s charge-balance rivals (ph_balance_rivals.py): every one OFF; inline dynamic pH only
    _k("ph.amine_survival", "ph_balance_rivals", "AMINE_SURVIVAL", _pbr.AmineSurvival(), "rate",
       "M1b prereg sec. 3.1, rival B1 (phi_Cys, phi_ARP)"),
    _k("ph.acid_lag", "ph_balance_rivals", "ACID_LAG", _pbr.AcidLag(), "rate",
       "M1b prereg sec. 3.1, rival B2 (f_ret, the latent carboxyl's lag)"),
    _k("ph.acid_supply", "ph_balance_rivals", "ACID_SUPPLY", _pbr.AcidSupply(), "rate",
       "M1b prereg secs. 3.1 - 3.1.1, rival B3 (Davidek 2006a's acetic share through a latent step; both branch "
       "forms)"),
    _k("ph.arp_forcing", "ph_balance_rivals", "ARP_FORCING", _pbr.ArpForcing(), "rate",
       "M1b prereg sec. 3.3, C-Hou: Hou 2026's measured ARP(t) as an input"),
    _k("ph.balance_refusal", "ph_balance_refusal", "ENABLED", True, "rate",
       "the S3 application of M1b part (i) (2026-11-07): engine.predict refuses the outputs that depend on the "
       "falsified pH balance. Off only where the model's own numbers are computed, not shipped: a fit's objective "
       "(build_frozen_fit_reproducibility.py) and the S3 sensitivity run"),
    # --- the sulfur lane's integration and flux budget (sulfur.py)
    _k("sulfur.flux_budget", "sulfur", "FLUX_BUDGET", "closed", "solver",
       "B60; frozen waves set linear_rederived", choices=("closed", "linear_rederived")),
    _k("sulfur.profile_node_step_k", "sulfur", "PROFILE_NODE_STEP_K", 2.0, "solver",
       "B64: largest temperature step between two rate nodes along a recorded profile"),
    _k("engine.ramp_step_k", "engine", "RAMP_STEP_K", 2.0, "solver",
       "ThermalProgram.from_profile's default segment step (a call-time read found by the phase 2 coverage scan)"),
    # --- the sulfur lane's declared routes (parameters_sulfur.py), each OFF or at its shipped value
    _k("sulfur.ttca_formation", "parameters_sulfur", "TTCA_FORMATION", _ps.TtcaFormation(), "rate",
       "B58 (kinetic_core_b58_prereg.md)"),
    _k("sulfur.nf_hydroxide_decay", "parameters_sulfur", "NF_HYDROXIDE_DECAY", _ps.NfHydroxideDecay(), "rate",
       "B61 (kinetic_core_b61_prereg.md)"),
    _k("sulfur.pentose_mgo_routes", "parameters_sulfur", "PENTOSE_MGO_ROUTES", _ps.PentoseMgoRoutes(), "rate",
       "sub-model 3 (2026-10-30), rival methylglyoxal arms"),
    _k("sulfur.pentose_amine_entry", "parameters_sulfur", "PENTOSE_AMINE_ENTRY", _ps.PentoseAmineEntry(), "rate",
       "kinetic_core_pentose_amine_entry_prereg.md (2026-10-30)"),
    _k("sulfur.ddp_formation_ph", "parameters_sulfur", "DDP_FORMATION_PH", _ps.DdpFormationPh(), "rate",
       "round 3, thiol prereg secs. 9-10 (2026-10-30)"),
    _k("sulfur.mgo_cysteine_strecker", "parameters_sulfur", "MGO_CYSTEINE_STRECKER", _ps.MgoCysteineStrecker(),
       "rate", "sub-models 3 + 4, joint-fit design review B5 (2026-10-30)"),
    _k("sulfur.tdp_cysteine_trap", "parameters_sulfur", "TDP_CYSTEINE_TRAP", _ps.TdpCysteineTrap(), "rate",
       "B7 hub feasibility, kinetic_core_dpo_partition_prereg.md sec. 20"),
    _k("sulfur.hexose_c5_entry", "parameters_sulfur", "HEXOSE_C5_ENTRY", _ps.HexoseC5Entry(), "rate",
       "sub-model 4 (2026-10-30)"),
    _k("sulfur.strecker_sulfur", "parameters_sulfur", "STRECKER_SULFUR", _ps.StreckerSulfur(), "rate",
       "B59a (2026-10-28), Kobayasi & Fujimaki 1965"),
    _k("sulfur.cys_amadori_route", "parameters_sulfur", "CYS_AMADORI_ROUTE", _ps.CysAmadoriRoute(), "rate",
       "B59b (2026-10-29), de Roos 2005; B7a split law"),
    _k("sulfur.tdp_own_decay", "parameters_sulfur", "TDP_OWN_DECAY", _ps.TdpOwnDecay(), "rate",
       "B7a arm R-T, kinetic_core_dpo_partition_prereg.md sec. 8"),
    _k("sulfur.furfural_barriers", "parameters_sulfur", "FURFURAL_BARRIERS", _ps.FurfuralBarriers(), "rate",
       "B62 (kinetic_core_b62_prereg.md)"),
    _k("sulfur.single_temperature_barriers", "parameters_sulfur", "SINGLE_TEMPERATURE_BARRIERS",
       _ps.SingleTemperatureBarriers(), "rate", "B63 route-class median barriers"),
    # --- fit-side declarations of plan step M1 (parameters_sulfur.py): read by fit code only, never by the kinetics
    _k("fit.dimer_encounter_ceiling", "parameters_sulfur", "DIMER_ENCOUNTER_CEILING", _ps.DimerEncounterCeiling(),
       "rate", "M1 prereg sec. 2.5: the encounter-limit ceiling on k_dimer_* replaces the numerical bound in a fit"),
    _k("fit.mft_decay_prior", "parameters_sulfur", "MFT_DECAY_PRIOR", _ps.MftDecayPrior(), "rate",
       "M1 prereg sec. 2.6: a weak Gaussian prior on log10 lambda_mft(100 C) as one residual of a fit's objective"),
    # --- the trunk (trunk_conditions.py)
    _k("trunk.three_deoxy_exit_ph_term", "trunk_conditions", "THREE_DEOXY_EXIT_PH_TERM", True, "rate",
       "B40: the measured 3-DG exit pH terms"),
    _k("trunk.three_deoxy_exit_ph", "trunk_conditions", "THREE_DEOXY_EXIT_PH", (
        ("k_tdg_fa", (0.77, (0.65, 0.89), "Martins 2003 Table 3 k6: 1.9e-3 vs 2.74e-2 (100 C), 4.30e-2 vs 3.04e-1 "
                                          "(120 C), pH 5.5 vs 6.8")),), "rate",
       "B41: the formic-acid exit only (B40's k_tdg_mgo term rejected by the Leitzen methylglyoxal row)", mapping=True),
    # --- the aqueous glyoxal supply block (parameters_dicarbonyl.py): B21 sets trial values through it
    _k("dicarbonyl.aqueous_glyoxal_parameters", "parameters_dicarbonyl", "AQUEOUS_GLYOXAL_PARAMETERS",
       tuple(_pd.with_aqueous_glyoxal(_pd.FROZEN_B21["log10_k_ama_g_100C"],
                                      _pd.FROZEN_B21["log10_k_g_go_aqueous_100C"]).items()), "rate",
       "B21 (kinetic_core_b21_prereg.md): FROZEN_B21's optimum", mapping=True),
    # --- the lipid lane (parameters_lipid.py, re-imported by lipid.py and the package)
    _k("lipid.q10_assumption", "parameters_lipid", "Q10_ASSUMPTION", _pl.Q10Assumption(), "rate",
       "Schroen & Berton-Carabin 2022: a declared Q10 with a band and a mandatory warning", also=("lipid", "")),
    _k("lipid.decadienal_breakage", "lipid", "DECADIENAL_BREAKAGE", _pl.DecadienalBreakage(), "rate",
       "Zamora 2015b: 2,4-decadienal thermal breakage", also=("parameters_lipid",)),
    _k("lipid.covalent_sink", "parameters_lipid", "COVALENT_SINK", _pl.CovalentSinkCeiling(), "rate",
       "covalent aldehyde sink ceiling, off", also=("lipid", "")),
    # --- shipped / unshipped fitted blocks (read at import by their parameter tables: build phase)
    _k("dicarbonyl.shipped_b39", "parameters_dicarbonyl", "SHIPPED_B39", True, "build",
       "B39 fed 3-deoxy block (results/validation/kinetic_core_b41_fit_report.json)"),
    _k("methionine.shipped_b22", "parameters_methionine", "METHIONINE_SHIPPED", False, "build",
       "B22, METHIONINE_NOT_SHIPPED_REASON"),
    _k("methionine.shipped_b22b", "parameters_methionine", "B22B_SHIPPED", False, "build", "B22b"),
    _k("proline.shipped_b24", "parameters_proline", "PROLINE_SHIPPED", False, "build",
       "B24, PROLINE_NOT_SHIPPED_REASON"),
    _k("proline.shipped_b24b", "parameters_proline", "B24B_SHIPPED", False, "build", "B24b"),
    # --- the panel and the uncertainty envelope
    _k("panel.thermal_history", "panel", "THERMAL_HISTORY", "isothermal", "rate",
       "B64: the bundle's recorded heating program, or an isothermal hold", choices=("isothermal", "as_recorded")),
    _k("uncertainty.b55_single_temperature_ea_prior", "uncertainty", "B55_SINGLE_TEMPERATURE_EA_PRIOR", True, "build",
       "B55: the measured-class Ea prior on single-temperature rows"),
    _k("uncertainty.b56_structural_discrepancy", "uncertainty", "B56_STRUCTURAL_DISCREPANCY", False, "build",
       "B56 structural discrepancy term, off"),
)

KEY_BY_NAME: Mapping[str, ConfigKey] = MappingProxyType({key.name: key for key in KEYS})
if len(KEY_BY_NAME) != len(KEYS):
    raise RuntimeError("duplicate ConfigKey names")
if set(_BANDS) - set(KEY_BY_NAME):
    raise RuntimeError(f"field bands for unknown keys: {sorted(set(_BANDS) - set(KEY_BY_NAME))}")
for _key in KEYS:   # a shipped default is a value its own key accepts
    _key.validate(_key.freeze(_key.default))
BUILD_KEYS: Tuple[str, ...] = tuple(key.name for key in KEYS if key.phase == "build")
#: ``(module, attribute) -> key name`` for every binding: how a script that named a global finds its key
KEY_OF_BINDING: Mapping[Tuple[str, str], str] = MappingProxyType(
    {binding: key.name for key in KEYS for binding in key.bindings})


def _module(name: str):
    return importlib.import_module(f"{_PACKAGE}.{name}" if name else _PACKAGE)


def _canonical(value: Any) -> Any:
    """A JSON-serialisable form of a Config value that two equal values share, whatever process made them: numbers
    as exact float hex (so 40 and 40.0, which are equal, agree), dataclasses tagged with their qualified type and
    every field in declaration order, sequences in order, mappings sorted by their canonical key."""
    if value is None or isinstance(value, (bool, str)):
        return value
    if _is_number(value):
        return {"num": float(value).hex()}
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {"type": f"{type(value).__module__}.{type(value).__qualname__}",
                "fields": [[f.name, _canonical(getattr(value, f.name))] for f in dataclasses.fields(value)]}
    if isinstance(value, (tuple, list)):
        return {"seq": [_canonical(v) for v in value]}
    if isinstance(value, Mapping):
        items = [[_canonical(k), _canonical(v)] for k, v in value.items()]
        return {"map": sorted(items, key=lambda kv: json.dumps(kv[0], sort_keys=True))}
    if isinstance(value, (frozenset, set)):
        return {"set": sorted((_canonical(v) for v in value), key=lambda c: json.dumps(c, sort_keys=True))}
    raise TypeError(f"no canonical form for a {type(value).__name__} ({value!r})")


@dataclass(frozen=True)
class Config:
    """A value for EVERY ``ConfigKey``: frozen, hashable, picklable, validated (type, choices, band, field bands) on
    construction. Build one with ``shipped()``, ``from_globals()`` or ``replace``; never by hand."""

    items: Tuple[Tuple[str, Any], ...]

    def __post_init__(self):
        names = [name for name, _ in self.items]
        if names != sorted(KEY_BY_NAME):
            missing = sorted(set(KEY_BY_NAME) - set(names))
            extra = sorted(set(names) - set(KEY_BY_NAME))
            raise ValueError(f"a Config sets every key once, in name order (missing {missing}, unknown {extra})")
        for name, value in self.items:
            KEY_BY_NAME[name].validate(value)
            try:
                hash(value)
            except TypeError as exc:
                raise TypeError(f"{name}: a Config value must be hashable ({exc})") from None

    def digest(self) -> str:
        """A sha256 of the canonical serialisation: the same for equal Configs in any process (``hash()`` is not:
        it is salted per process for strings). For records and cache keys."""
        blob = json.dumps([[name, _canonical(value)] for name, value in self.items], sort_keys=True,
                          separators=(",", ":"))
        return hashlib.sha256(blob.encode("utf-8")).hexdigest()

    @classmethod
    def _of(cls, values: Mapping[str, Any]) -> "Config":
        return cls(tuple((name, KEY_BY_NAME[name].freeze(values[name])) for name in sorted(values)))

    @property
    def values(self) -> Mapping[str, Any]:
        return MappingProxyType(dict(self.items))

    def __getitem__(self, name: str) -> Any:
        return KEY_BY_NAME[name].thaw(self.values[name])

    @classmethod
    def shipped(cls) -> "Config":
        """The shipped configuration, from the declarations (whatever the modules hold now)."""
        return cls._of({key.name: key.default for key in KEYS})

    @classmethod
    def from_globals(cls) -> "Config":
        """The configuration the engine's globals hold now (A3 migration bridge; deleted at the end of A2). Refuses
        a switch whose bindings disagree."""
        values = {}
        for key in KEYS:
            seen = [key.freeze(getattr(_module(m), a)) for m, a in key.bindings]
            if any(v != seen[0] for v in seen[1:]):
                where = ", ".join(f"{m or _PACKAGE}.{a}" for m, a in key.bindings)
                raise ValueError(f"{key.name}: its bindings disagree ({where}); set them together")
            values[key.name] = seen[0]
        return cls._of(values)

    def replace(self, changes: Optional[Mapping[str, Any]] = None, **kw: Any) -> "Config":
        """A copy with some keys changed: ``cfg.replace({"ph.coupling": "inline"})``, or as keywords with ``__`` for
        ``.`` (``cfg.replace(ph__coupling="inline")``). An unknown key is an error."""
        updates = dict(changes or {})
        updates.update({name.replace("__", "."): value for name, value in kw.items()})
        unknown = sorted(set(updates) - set(KEY_BY_NAME))
        if unknown:
            raise KeyError(f"unknown config keys {unknown}")
        values = dict(self.items)
        values.update({name: KEY_BY_NAME[name].freeze(value) for name, value in updates.items()})
        return Config(tuple(sorted(values.items())))

    def diff(self, other: Optional["Config"] = None) -> Dict[str, Tuple[Any, Any]]:
        """``{key: (other's value, this value)}`` for every key that differs (from ``shipped()`` by default)."""
        base = dict((other or Config.shipped()).items)
        return {name: (base[name], value) for name, value in self.items if base[name] != value}


def _snapshot() -> Dict[Tuple[str, str], Any]:
    return {(m, a): getattr(_module(m), a) for key in KEYS for m, a in key.bindings}


def _restore(snapshot: Mapping[Tuple[str, str], Any]) -> None:
    for (m, a), value in snapshot.items():
        setattr(_module(m), a, value)


class BuildPhaseKeyError(ValueError):
    """A Config that changes a build-phase key was to be installed (the owner's decision 1, 2026-11-01)."""


def _refuse_build_changes(config: Config) -> None:
    changed = [name for name in BUILD_KEYS if config.values[name] != KEY_BY_NAME[name].freeze(KEY_BY_NAME[name].default)]
    if changed:
        detail = "; ".join(f"{name} = {config.values[name]!r} (shipped {KEY_BY_NAME[name].default!r}, "
                           f"{', '.join(f'{m or _PACKAGE}.{a}' for m, a in KEY_BY_NAME[name].bindings)})"
                           for name in changed)
        raise BuildPhaseKeyError(
            f"refusing to install build-phase key(s) at a non-shipped value: {detail}. Their tables are built at "
            "import, so setting them at run time would switch only the call-time readers; change the shipped "
            "default (in the module and here) or run in a process that imports them so")


def install(config: Config) -> None:
    """Set every run-time switch (every binding of every ``"rate"`` / ``"solver"`` key) from ``config``, and leave
    it set: what a worker does at the start of a task. A build-phase key is a read-only fact: never set, and a
    Config that changes one is refused (``BuildPhaseKeyError``), before anything is set."""
    if not isinstance(config, Config):
        raise TypeError(f"install takes a Config, not {type(config).__name__}")
    _refuse_build_changes(config)
    values = dict(config.items)
    for key in KEYS:
        if key.phase == "build":
            continue
        value = key.thaw(values[key.name])
        for m, a in key.bindings:
            setattr(_module(m), a, value)


def update(changes: Optional[Mapping[str, Any]] = None, **kw: Any) -> Config:
    """What a script calls instead of assigning an engine global (A3 phase 2): the process's configuration with
    ``changes`` replaced (validated; every binding; build keys refused), installed and returned, so the same value
    can travel with the script's jobs (``parallel.map_ordered(fn, items, config=cfg)``)."""
    config = Config.from_globals().replace(changes, **kw)
    install(config)
    return config


@contextlib.contextmanager
def applied(config: Config) -> Iterator[Config]:
    """Run a block under ``config``; every binding gets back the very object it held, also on an exception."""
    snapshot = _snapshot()
    try:
        install(config)
        yield config
    finally:
        _restore(snapshot)


__all__ = ["BUILD_KEYS", "BuildPhaseKeyError", "Config", "ConfigKey", "KEYS", "KEY_BY_NAME", "KEY_OF_BINDING",
           "applied", "install", "update"]
