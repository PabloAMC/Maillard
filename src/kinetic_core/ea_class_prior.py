"""
src/kinetic_core/ea_class_prior.py -- WAVE B55: an activation-energy prior for the sulfur constants
known at ONE temperature only (pre-registered in results/validation/kinetic_core_b55_b56_prereg.md).

THE GAP. 31 of the shipped sulfur wave's 43 rate constants are quoted as log10 k at 145 C with a
barrier nobody measured: 30 share the lumped formation barrier (64.08 kJ/mol, UNIDENTIFIED in B2.3,
frozen since B9) and `k_mmft` has none (held at 145 C, policy 2). The envelope sampled none of the 31,
so every interval away from 145 C asserted a temperature slope as if it were known.

THE PRIOR. For each such constant, the MEASURED activation energies of its reaction class in this
corpus -- every value read from the registry that already carries it, never re-typed (the four Chan &
Reineccius 1994 values live only in a source string; ``CHAN_1994_TABLE_I_KJ_MOL`` transcribes them once
and a unit test asserts the string still prints them). Two classes, keyed by the route table the sulfur
module already declares (``FORMATION_ROUTE_OF``, B10):

  * ``carbohydrate`` (route ``sugar_trunk``): steps whose reactant is a sugar, an Amadori / Schiff
    product or a sugar-derived carbonyl, with no sulfur-bearing reactant;
  * ``sulfur`` (route ``thiol_assembly``, and the route-less thiol condensation ``k_mmft``): steps with
    a sulfur-bearing reactant.

A draw is the class's EMPIRICAL distribution (inverse empirical CDF, linear interpolation, at a uniform
quantile): it never leaves the measured range and assumes no shape. k(145 C) is held where the fit put
it, so a 145 C prediction is unchanged by construction; only the slope away from it is uncertain.

Nothing here moves a point prediction. The engine reads a drawn barrier only through the override
block ``single_temperature_ea_kj_mol`` a Monte-Carlo draw carries.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from types import MappingProxyType
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

#: The classes and which target route each serves.
CARBOHYDRATE = "carbohydrate"
SULFUR_CLASS = "sulfur"
CLASS_OF_ROUTE: Mapping[Optional[str], str] = {
    "sugar_trunk": CARBOHYDRATE,
    "thiol_assembly": SULFUR_CLASS,
    # `k_mmft` (thiol + methanethiol condensation) carries no route: no barrier at all (policy 2).
    None: SULFUR_CLASS,
}

#: Chan & Reineccius 1994, ACS Symp. Ser. 564 ch. 10, Table I (replicate means): methional at pH 6 / 7
#: / 8 and 2-acetylthiophene. Printed in parameters_sulfur.FORMATION_EA_PRIOR_SOURCE["thiol_assembly"];
#: tests/unit/test_kinetic_core_b55_b56.py asserts that string still prints each of them.
CHAN_1994_TABLE_I_KJ_MOL: Tuple[Tuple[str, float], ...] = (
    ("methional, pH 6", 102.1),
    ("methional, pH 7", 106.6),
    ("methional, pH 8", 81.4),
    ("2-acetylthiophene", 132.9),
)

#: Registry parameters whose printed barrier is a measurement: (class, lane, parameter key).
#: The lane names the registry the key is read from; the value, the evidence class and the source
#: are read from it at call time.
_REGISTRY_ROWS: Tuple[Tuple[str, str, str], ...] = (
    # Martins & van Boekel 2005 Table 2 (glucose / glycine, 80-120 C)
    (CARBOHYDRATE, "trunk", "k_glc_fru"), (CARBOHYDRATE, "trunk", "k_fru_glc"),
    (CARBOHYDRATE, "trunk", "k_schiff"), (CARBOHYDRATE, "trunk", "k_ama_tdg"),
    (CARBOHYDRATE, "trunk", "k_ama_odg"), (CARBOHYDRATE, "trunk", "k_ama_mgo"),
    (CARBOHYDRATE, "trunk", "k_fru_acids"), (CARBOHYDRATE, "trunk", "k_odg_aa"),
    (CARBOHYDRATE, "trunk", "k_tdg_fa"), (CARBOHYDRATE, "trunk", "k_tdg_mel"),
    # Kocadagli & Gokmen 2016 Table 2 (glucose system)
    (CARBOHYDRATE, "trunk", "k_fru_int"), (CARBOHYDRATE, "trunk", "k_fru_odg"),
    (CARBOHYDRATE, "trunk", "k_glc_g"), (CARBOHYDRATE, "trunk", "k_glc_tdg"),
    (CARBOHYDRATE, "trunk", "k_int_hmf"), (CARBOHYDRATE, "trunk", "k_odg_da"),
    (CARBOHYDRATE, "trunk", "k_tdg_mgo"),
    # De Vleeschouwer et al. 2009 Part I Table 3 (glucose column)
    (CARBOHYDRATE, "acrylamide", "k_asn_glc"), (CARBOHYDRATE, "acrylamide", "k_int1_acr"),
    # Zheng & Ho 1994 Tables I and III
    (SULFUR_CLASS, "sulfur", "k_cys_h2s"), (SULFUR_CLASS, "sulfur", "k_gsh_h2s"),
    # Hamzalioglu & Gokmen 2018 Table 1 (HMF + cysteine)
    (SULFUR_CLASS, "trunk", "k_hmf_cys"),
    # De Vleeschouwer, Van der Plancken, Van Loey & Hendrickx 2009 (precursors; cysteine steps)
    (SULFUR_CLASS, "acrylamide", "k_acr_cys"), (SULFUR_CLASS, "acrylamide", "k_cys_sink"),
    (SULFUR_CLASS, "acrylamide", "k_cys_glc"),
)

#: Module constants of parameters_sulfur that are measured barriers: (class, constant name, what).
_CONSTANT_ROWS: Tuple[Tuple[str, str, str], ...] = (
    (CARBOHYDRATE, "ZHANG_EA_CYS_AMADORI_TO_ALPHA_DC_KJ_MOL",
     "Zhang 2026 k16: Cys-Amadori -> alpha-dicarbonyl (R^2 1.000)"),
    (SULFUR_CLASS, "ZHANG_EA_THIOL_TO_DISULFIDE_KJ_MOL", "Zhang 2026 k17: thiol -> disulfide (R^2 0.971)"),
    (SULFUR_CLASS, "ZHANG_EA_CYS_TO_CYS_AMADORI_KJ_MOL", "Zhang 2026: cysteine -> Cys-Amadori"),
    (SULFUR_CLASS, "KANG_EA_FREE_CYS_DEPLETION_KJ_MOL", "Kang 2026: free-cysteine depletion"),
    (SULFUR_CLASS, "ZHAI_EA_FREE_CYS_DEPLETION_PRIMARY_KJ_MOL",
     "Zhai 2021: free-cysteine depletion (a second laboratory)"),
    (SULFUR_CLASS, "GIGL_EA_THIOL_COVALENT_CAPTURE_KJ_MOL",
     "Gigl 2021: thiol covalent capture (central value of the printed 7-102 kJ/mol range)"),
)

#: Evidence classes whose barrier counts as measured.
MEASURED_EVIDENCE = frozenset({"measured_rate", "measured_activation_energy"})


@dataclass(frozen=True)
class MeasuredEa:
    """One measured barrier in the class prior."""

    klass: str
    ea_kj_mol: float
    label: str
    source: str


@lru_cache(maxsize=None)
def measured_activation_energies() -> Tuple[MeasuredEa, ...]:
    """Every measured barrier in the two classes, resolved from the registries now."""
    from . import parameters_sulfur as ps
    from .engine import core_parameters

    registries: Dict[str, Mapping] = {}
    out: List[MeasuredEa] = []
    for klass, lane, key in _REGISTRY_ROWS:
        if lane not in registries:
            registries[lane] = core_parameters(lane)
        param = registries[lane][key]
        ec = str(getattr(param, "evidence_class", ""))
        ea = getattr(param, "ea_kj_mol", None)
        if ec not in MEASURED_EVIDENCE or ea is None or float(ea) <= 0.0:
            raise ValueError(
                f"{lane}.{key} is listed as a measured barrier but carries evidence_class {ec!r} and "
                f"Ea {ea!r}; the B55 prior takes measurements only"
            )
        out.append(MeasuredEa(klass, float(ea), f"{lane}.{key}", str(getattr(param, "source_anchor", ""))[:160]))
    for klass, name, what in _CONSTANT_ROWS:
        out.append(MeasuredEa(klass, float(getattr(ps, name)), f"parameters_sulfur.{name}", what))
    source = ps.FORMATION_EA_PRIOR_SOURCE["thiol_assembly"]
    for what, value in CHAN_1994_TABLE_I_KJ_MOL:
        out.append(MeasuredEa(SULFUR_CLASS, float(value), f"Chan & Reineccius 1994 Table I, {what}", source[:160]))
    return tuple(out)


@lru_cache(maxsize=None)
def class_values() -> Mapping[str, Tuple[float, ...]]:
    """``{class: sorted measured barriers}``."""
    out: Dict[str, List[float]] = {CARBOHYDRATE: [], SULFUR_CLASS: []}
    for row in measured_activation_energies():
        out[row.klass].append(row.ea_kj_mol)
    return MappingProxyType({k: tuple(sorted(v)) for k, v in out.items()})


def empirical_quantile(values: Sequence[float], u: float) -> float:
    """The inverse empirical CDF at ``u`` in [0, 1], linear interpolation between order statistics."""
    return float(np.quantile(np.asarray(values, dtype=float), float(u), method="linear"))


@lru_cache(maxsize=None)
def single_temperature_targets() -> Tuple[Tuple[str, str, Optional[float]], ...]:
    """
    ``((key, class, shipped barrier or None), ...)``: every key of the SHIPPED sulfur vector whose barrier
    is neither a measurement (MEASURED_EA_OVERRIDES) nor a fitted decay-family barrier. Computed from
    the shipped fit report, never listed by hand, so a later wave that measures a barrier drops it here.
    """
    from . import engine
    from . import parameters_sulfur as ps

    frozen = json.loads(engine._B2_FIT_REPORT.read_text(encoding="utf-8"))["frozen_parameters"]
    lumped = float(frozen["lumped_formation_Ea_kJ_mol"])
    out = []
    for key in frozen["log10_k_ref_at_145C"]:
        if key in ps.measured_ea_overrides() or ps.decay_family_of(key) is not None:
            continue
        route = ps.FORMATION_ROUTE_OF.get(key)
        shipped = None if key in ps.NO_EA_KEYS else lumped
        if shipped is not None and ps.SINGLE_TEMPERATURE_BARRIERS.enabled and route in ps.ROUTE_CLASS_MEDIAN_EA_KJ_MOL:
            shipped = float(ps.ROUTE_CLASS_MEDIAN_EA_KJ_MOL[route])  # B63: the point sits at its class's median
        out.append((key, CLASS_OF_ROUTE[route], shipped))
    return tuple(out)


def source_line() -> str:
    """One line naming the classes, their sizes and ranges, for a prior row's ``source``."""
    vals = class_values()
    parts = [
        f"{k}: n={len(v)}, median {np.median(v):.1f}, range {min(v):.1f}-{max(v):.1f} kJ/mol"
        for k, v in vals.items()
    ]
    return "B55 measured-barrier class prior (src/kinetic_core/ea_class_prior.py; " + "; ".join(parts) + ")"


__all__ = [
    "CARBOHYDRATE",
    "CHAN_1994_TABLE_I_KJ_MOL",
    "CLASS_OF_ROUTE",
    "MeasuredEa",
    "SULFUR_CLASS",
    "class_values",
    "empirical_quantile",
    "measured_activation_energies",
    "single_temperature_targets",
    "source_line",
]

