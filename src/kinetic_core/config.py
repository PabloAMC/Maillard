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

Next (A3 phase 2): the scripts build ``Config.replace(...)`` instead of assigning globals; then the persistent pool.
"""
from __future__ import annotations

import contextlib
import dataclasses
import importlib
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, Dict, Iterator, Literal, Mapping, Optional, Tuple

from . import parameters_dicarbonyl as _pd
from . import parameters_lipid as _pl
from . import parameters_sulfur as _ps
from . import ph_state as _ph

Phase = Literal["build", "rate", "solver"]
_PACKAGE = "src.kinetic_core"


@dataclass(frozen=True)
class ConfigKey:
    """One switch. ``bindings`` are ``(module, attribute)`` pairs, module relative to ``src.kinetic_core`` ("" for the
    package itself); the first is the one the engine's kinetics read. ``mapping`` keys hold a dict, kept in a Config
    as a tuple of its items (hashable, order kept) and installed as a fresh dict."""

    name: str
    bindings: Tuple[Tuple[str, str], ...]
    default: Any
    phase: Phase
    basis: str
    band: Optional[Tuple[Any, Any]] = None
    choices: Optional[Tuple[Any, ...]] = None
    mapping: bool = False

    def freeze(self, value: Any) -> Any:
        if self.mapping:
            return tuple((k, v) for k, v in dict(value).items())
        return value

    def thaw(self, value: Any) -> Any:
        return dict(value) if self.mapping else value


def _k(name: str, module: str, attr: str, default: Any, phase: Phase, basis: str, *, also: Tuple[str, ...] = (),
       band=None, choices=None, mapping: bool = False) -> ConfigKey:
    bindings = ((module, attr),) + tuple((m, attr) for m in also)
    return ConfigKey(name, bindings, default, phase, basis, band=band, choices=choices, mapping=mapping)


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
    # --- the sulfur lane's integration and flux budget (sulfur.py)
    _k("sulfur.flux_budget", "sulfur", "FLUX_BUDGET", "closed", "solver",
       "B60; frozen waves set linear_rederived", choices=("closed", "linear_rederived")),
    _k("sulfur.profile_node_step_k", "sulfur", "PROFILE_NODE_STEP_K", 2.0, "solver",
       "B64: largest temperature step between two rate nodes along a recorded profile"),
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


def _module(name: str):
    return importlib.import_module(f"{_PACKAGE}.{name}" if name else _PACKAGE)


@dataclass(frozen=True)
class Config:
    """A value for EVERY ``ConfigKey``: frozen, hashable, picklable. Build one with ``shipped()``, ``from_globals()``
    or ``replace``; never by hand."""

    items: Tuple[Tuple[str, Any], ...]

    def __post_init__(self):
        names = [name for name, _ in self.items]
        if names != sorted(KEY_BY_NAME):
            missing = sorted(set(KEY_BY_NAME) - set(names))
            extra = sorted(set(names) - set(KEY_BY_NAME))
            raise ValueError(f"a Config sets every key once, in name order (missing {missing}, unknown {extra})")
        for name, value in self.items:
            key = KEY_BY_NAME[name]
            if key.choices is not None and value not in key.choices:
                raise ValueError(f"{name} = {value!r}: not one of {key.choices}")
            try:
                hash(value)
            except TypeError as exc:
                raise TypeError(f"{name}: a Config value must be hashable ({exc})") from None

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


def install(config: Config) -> None:
    """Set every switch (every binding of every key) from ``config``, and leave it set: what a worker does at the
    start of a task, and what a script's ``install(job)`` does in the process that runs the job."""
    if not isinstance(config, Config):
        raise TypeError(f"install takes a Config, not {type(config).__name__}")
    values = dict(config.items)
    for key in KEYS:
        value = key.thaw(values[key.name])
        for m, a in key.bindings:
            setattr(_module(m), a, value)


@contextlib.contextmanager
def applied(config: Config) -> Iterator[Config]:
    """Run a block under ``config``; every binding gets back the very object it held, also on an exception."""
    snapshot = _snapshot()
    try:
        install(config)
        yield config
    finally:
        _restore(snapshot)


__all__ = ["Config", "ConfigKey", "KEYS", "KEY_BY_NAME", "applied", "install"]
