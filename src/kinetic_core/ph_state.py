"""
src/kinetic_core/ph_state.py

BUILD WAVE B2.2 -- THE pH TRAJECTORY AS A COMPUTED STATE, NOT A PRESCRIBED ONE.
=============================================================================

WHAT B2.1 HAD, AND WHY IT WAS NOT A STATE
-----------------------------------------
B2.1's ``integrate_sulfur(..., ph_final=...)`` moved the pH LINEARLY IN TIME
between a measured initial and a measured FINAL pH. That is interpolation
between two measurements, not a model: it can only be used where somebody has
already measured the answer, it cannot be asked what a buffer would do, and it
cannot be run at all on a system whose endpoint nobody published. The
declaration's sec.5 decision 1 scores Zhou's pH-6 column DIAGNOSTIC ONLY
"until the model carries a pH-trajectory state". A prescribed endpoint is not
that state.

WHAT THIS MODULE DOES INSTEAD
-----------------------------
The pH is an ALGEBRAIC FUNCTION OF THE STATE, solved at every point on the
trajectory from a charge balance:

    SID + [H+] - Kw/[H+] + SUM_i C_i * z_i(pH, T) = 0

  * ``SID`` is the STRONG ION DIFFERENCE (net strong cation, e.g. the NaOH used
    to set the initial pH, or the Na+ of a phosphate buffer). It is CONSERVED:
    no reaction in the network creates or destroys a strong ion, so SID is
    computed ONCE at t = 0 from the declared initial pH and the initial
    composition, and never moves again. This is the quantity that makes the
    three Zhou runs different from each other: an unbuffered pot taken to
    pH 8 with NaOH carries far more titratable base into the hold than the same
    pot taken to pH 6, and that -- not the pH label -- is what survives the
    heating.

  * ``C_i`` are the TRACKED pools: the organic-acid pool the network itself
    produces, plus whatever titratable solutes were charged (cysteine, the
    Amadori compound), plus the declared buffer.

  * ``z_i(pH, T)`` is the average charge of a polyprotic group at reaction
    temperature, from its pKa ladder corrected by van 't Hoff -- the SAME
    correction ``parameters_sulfur.pka_at`` already applies to the speciation
    factors, so the pH the rate law sees and the pH the charge balance solves
    are on ONE scale.

Because protonation equilibria are many orders of magnitude faster than any
step in this network, treating [H+] as an instantaneous function of the state
rather than as an ODE variable is the correct quasi-equilibrium reduction, and
it is also what keeps the system non-stiff.

WHAT IS FITTED HERE
-------------------
TWO constants, and they are calibrated on the THREE anchors
`FIT_HOLDOUT_DECLARATION.md` Amendment 7 declares FIT for exactly this purpose
(Zhou 2023 Fig. 2, ARP + Cys, initial 6/7/8 -> final 3.22/3.42/5.07):

  1. ``ACID_YIELD_PER_SINK_EVENT`` -- the titratable-acid yield of the lumped
     deoxyosone sink. The sink is a LUMP over everything the corpus never
     identified; Martins 2005 measures that PART of the hexose deoxyosone flux
     terminates as formic and acetic acid (steps 5 and 8) while another part
     terminates as browning polymer (step 9), so the acid yield of the pentose
     analogue is a fraction between 0 and 1 and nothing in the corpus measures
     it.
  2. ``ARP_AMINE_PKA`` -- the pKa of the Amadori compound's secondary
     ammonium. It sets how much NaOH went in to reach the declared initial pH,
     i.e. it sets SID. No value for it exists in this corpus or in the
     dossiers.

EVERY OTHER NUMBER IN THIS FILE IS TEXTBOOK AND IMMOVABLE, and each carries its
citation inline. None of them is a rate.

TWO pH SCALES, AND WHY BOTH ARE NEEDED
--------------------------------------
A pH meter reads a COOLED sample. Every pH number printed in this corpus --
Zhou's 6/7/8 labels, Zhou's 3.22/3.42/5.07 endpoints, Kumazawa's cold-pH grid,
Hofmann's "pH 5.0" -- is therefore a 25 C reading, while the rate law needs the
IN-SITU pH at 100-145 C. The two differ, and NOT by a constant: cysteine's
thiol pKa falls from 8.33 to 7.23 between 25 and 120 C, so a solution titrated
to a cold pH 7 has ~1 mmol/L of thiolate at the bench and ~8 mmol/L in the
autoclave. Conflating the scales puts a factor of eight into the base
inventory of exactly the systems the drift model is calibrated on, and this
module was measured doing so before the distinction was made explicit.

So the module carries both, and the conversion between them is free of any
fitted quantity:

  * the DECLARED LABEL pH is read at ``label_temperature_k`` (25 C by default)
    and is used ONCE, to compute the conserved SID;
  * the IN-SITU pH is solved from the same SID at the reaction temperature and
    is what every rate constant sees;
  * the COOLED pH is solved from the same SID at 25 C from the FINAL state,
    and is the only quantity comparable with a published endpoint.

At the acidic endpoints the two scales nearly coincide -- the ionisation
enthalpies of formic and acetic acid are ~0, so their pKa are
temperature-invariant -- which is why the difference bites at the START of a
run and not at its end.

UNITS: mmol/L throughout (the network's unit), Kelvin, kJ/mol.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

R_KJ = 8.314462618e-3

PKA_REFERENCE_T_K = 298.15


# ===========================================================================
# 1. THE TEXTBOOK ACID-BASE CONSTANTS
# ===========================================================================
# pKa at 25 C and standard ionisation enthalpy, kJ/mol. SOURCES, stated per
# block because this is the one place the module reaches outside its own
# dossiers, exactly as parameters_sulfur.PKA_25C does:
#
#   [CRC]  CRC Handbook of Chemistry and Physics, 95th ed. (2014), tables
#          "Dissociation Constants of Organic Acids and Bases" and
#          "Ionization Constants of Inorganic Acids and Bases".
#   [GKL]  Goldberg, Kishore & Lennen, "Thermodynamic Quantities for the
#          Ionization Reactions of Buffers", J. Phys. Chem. Ref. Data 31 (2002)
#          231-370 -- the standard compilation of ionisation enthalpies.
#   [MS]   Martell & Smith, Critical Stability Constants, vol. 3 (citrate).
#
# NONE of these is a rate, none is fitted, and none comes from a Maillard
# dossier. They are equilibrium constants of named small molecules.


@dataclass(frozen=True)
class TitratableGroup:
    """
    One polyprotic acid-base centre.

    ``pkas`` are the successive pKa of the FULLY PROTONATED form, in order.
    ``z_max`` is the charge of that fully protonated form (+1 for an amino
    acid's H3A+, 0 for a neutral carboxylic acid, 0 for H3PO4).
    """

    key: str
    label: str
    pkas: Tuple[float, ...]
    delta_h_kj_mol: Tuple[float, ...]
    z_max: float
    source: str

    def __post_init__(self) -> None:
        if len(self.pkas) != len(self.delta_h_kj_mol):
            raise ValueError(f"{self.key}: pKa and dH lists differ in length")
        if any(
            self.pkas[i] > self.pkas[i + 1] for i in range(len(self.pkas) - 1)
        ):
            raise ValueError(
                f"{self.key}: successive pKa must be non-decreasing; got "
                f"{self.pkas}"
            )

    def pkas_at(self, temperature_k: float) -> Tuple[float, ...]:
        """van 't Hoff, the same correction parameters_sulfur.pka_at applies."""
        scale = 1.0 / float(temperature_k) - 1.0 / PKA_REFERENCE_T_K
        return tuple(
            pka + (dh / (R_KJ * math.log(10.0))) * scale
            for pka, dh in zip(self.pkas, self.delta_h_kj_mol)
        )

    def average_charge(self, ph: float, temperature_k: float) -> float:
        """
        Mean charge of the group at ``ph``: ``z_max - n_bar``, where ``n_bar``
        is the mean number of protons lost.

        Computed from the standard polyprotic distribution in a form that
        cannot overflow: the un-normalised weights are carried in log10.
        """
        pkas = self.pkas_at(temperature_k)
        log_w = [0.0]
        running = 0.0
        for pka in pkas:
            running += float(ph) - pka
            log_w.append(running)
        top = max(log_w)
        weights = [10.0 ** (value - top) for value in log_w]
        total = sum(weights)
        n_bar = sum(i * w for i, w in enumerate(weights)) / total
        return float(self.z_max) - n_bar


#: pKw(T), TABULATED rather than van 't Hoff'd.
#:
#: A constant-enthalpy van 't Hoff extrapolation from 25 C is WRONG BY 0.23
#: UNITS AT 100 C and by 0.49 at 145 C -- water's ionisation has a large
#: negative dCp and the two-point form does not know it. Since every system in
#: this corpus is run between 100 and 145 C, the tabulated values are used and
#: interpolated linearly in 1/T.
#:
#: SOURCE: Marshall & Franck, "Ion product of water substance, 0-1000 C,
#: 1-10000 bars", J. Phys. Chem. Ref. Data 10 (1981) 295-304; the same values
#: are reproduced in the CRC Handbook's "Ionization Constant of Water" table.
PKW_TABLE_C: Tuple[Tuple[float, float], ...] = (
    (0.0, 14.938), (25.0, 13.995), (50.0, 13.275), (75.0, 12.700),
    (100.0, 12.265), (125.0, 11.911), (150.0, 11.638), (175.0, 11.431),
    (200.0, 11.289),
)
PKW_25C = 13.995


def pkw_at(temperature_k: float) -> float:
    """pKw at ``temperature_k``, interpolated linearly in 1/T on PKW_TABLE_C."""
    t = float(temperature_k)
    points = [(c + 273.15, p) for c, p in PKW_TABLE_C]
    if t <= points[0][0]:
        return points[0][1]
    if t >= points[-1][0]:
        return points[-1][1]
    for (t_lo, p_lo), (t_hi, p_hi) in zip(points, points[1:]):
        if t_lo <= t <= t_hi:
            x, x_lo, x_hi = 1.0 / t, 1.0 / t_lo, 1.0 / t_hi
            w = (x - x_lo) / (x_hi - x_lo)
            return p_lo + w * (p_hi - p_lo)
    raise AssertionError("unreachable")


#: --- the acids the NETWORK produces -------------------------------------
FORMIC = TitratableGroup(
    "formic", "formic acid (HCOOH)", (3.75,), (-0.13,), 0.0,
    "pKa 3.751 [CRC]; dH_ion -0.13 kJ/mol [GKL]")
ACETIC = TitratableGroup(
    "acetic", "acetic acid (CH3COOH)", (4.756,), (-0.41,), 0.0,
    "pKa 4.756 [CRC]; dH_ion -0.41 kJ/mol [GKL]")

#: The lumped organic-acid pool of the pentose lane. Its composition is not
#: measured; Martins' two named acids are formic and acetic, so the pool is
#: carried as a 1:1 mixture of the two and the resulting effective pKa (4.25)
#: is a DECLARED MIXTURE, not a fitted number. A pure-formic or pure-acetic
#: pool moves the predicted endpoints by ~0.4 pH units, which is reported as
#: the model's own composition sensitivity rather than tuned away.
ORGANIC_ACID_MIX: Tuple[Tuple[TitratableGroup, float], ...] = (
    (FORMIC, 0.5),
    (ACETIC, 0.5),
)

#: --- the titratable SOLUTES the panels actually charge -------------------
#: Cysteine, as the triprotic H3A+ it is.
CYSTEINE = TitratableGroup(
    "cysteine", "L-cysteine (COOH / SH / NH3+)",
    (1.92, 8.33, 10.28), (0.0, 26.0, 43.0), 1.0,
    "pKa 1.92 / 8.33 / 10.28 [CRC]; dH_ion 0 (carboxyl, |dH| < 2 kJ/mol) / "
    "26 / 43 kJ/mol [GKL]. The thiol and ammonium values are the SAME pair "
    "parameters_sulfur.PKA_25C already carries, so the charge balance and the "
    "speciation factors cannot disagree.")

#: The Amadori compound. Its carboxyl is an ordinary alpha-carboxyl; its
#: secondary ammonium is the CALIBRATED constant (see ARP_AMINE_PKA).
ARP_CARBOXYL_PKA = 3.00
ARP_CARBOXYL_DH = 0.0

#: TTCA -- 2-(tetrahydroxybutyl)thiazolidine-4-carboxylic acid. A thiazolidine
#: carboxylic acid: one carboxyl and one ring nitrogen. Values are those of the
#: parent thiazolidine-4-carboxylic acid.
TTCA_GROUP = TitratableGroup(
    "ttca", "thiazolidine-4-carboxylic acid (COOH / ring NH+)",
    (1.80, 6.24), (0.0, 40.0), 1.0,
    "pKa 1.80 / 6.24, thiazolidine-4-carboxylic acid [CRC, amino-acid table]; "
    "dH_ion 0 (carboxyl) / 40 kJ/mol (secondary ammonium) [GKL, amine class]. "
    "FLAG: the parent acid, not the tetrahydroxybutyl derivative -- the "
    "substituent is remote from both centres but no measurement of the "
    "derivative exists.")

#: --- the BUFFERS ----------------------------------------------------------
PHOSPHATE = TitratableGroup(
    "phosphate", "orthophosphate (H3PO4)",
    (2.148, 7.198, 12.35), (-7.9, 3.6, 16.0), 0.0,
    "pKa 2.148 / 7.198 / 12.35 [CRC]; dH_ion -7.9 / +3.6 / +16.0 kJ/mol [GKL]")
CITRATE = TitratableGroup(
    "citrate", "citric acid",
    (3.128, 4.761, 6.396), (4.07, 2.23, -3.38), 0.0,
    "pKa 3.128 / 4.761 / 6.396 [CRC]; dH_ion +4.07 / +2.23 / -3.38 kJ/mol "
    "[GKL, MS]")


# ===========================================================================
# 2. THE BUFFER SPEC -- an INPUT, with a declared default
# ===========================================================================


@dataclass(frozen=True)
class BufferSpec:
    """
    The declared buffer of one pot.

    ``kind``:
      * ``"none"``      -- unbuffered. Water autoprotolysis and the charged
                           solutes are the only base capacity. THIS IS THE
                           DECLARED DEFAULT when a caller supplies no spec, and
                           it raises the ``ph_trajectory_extrapolated`` warning,
                           because a pot whose buffer nobody recorded is a pot
                           whose pH trajectory is being extrapolated.
      * ``"phosphate"`` -- total orthophosphate, mol/L.
      * ``"citrate_phosphate"`` -- McIlvaine: citric acid AND disodium
                           hydrogen phosphate, both molarities required.
      * ``"clamped"``   -- the pH is HELD at the declared value. Reserved for
                           a source that states a pH-stat. It is NOT the
                           default and it is NOT applied to a source that
                           merely says "buffer": a buffer has a finite
                           capacity and this module computes it.
    """

    kind: str = "none"
    phosphate_mol_l: float = 0.0
    citrate_mol_l: float = 0.0
    declared: bool = False
    source: str = ""

    _KINDS = ("none", "phosphate", "citrate_phosphate", "clamped")

    def __post_init__(self) -> None:
        if self.kind not in self._KINDS:
            raise ValueError(
                f"unknown buffer kind {self.kind!r}; expected one of {self._KINDS}"
            )
        if self.kind == "phosphate" and self.phosphate_mol_l <= 0.0:
            raise ValueError("a phosphate buffer needs a positive molarity")
        if self.kind == "citrate_phosphate" and (
            self.phosphate_mol_l <= 0.0 or self.citrate_mol_l <= 0.0
        ):
            raise ValueError("a McIlvaine buffer needs both molarities")

    @property
    def is_clamped(self) -> bool:
        return self.kind == "clamped"

    def groups(self) -> Tuple[Tuple[TitratableGroup, float], ...]:
        """(group, mmol/L) pairs contributed by the buffer itself."""
        if self.kind == "phosphate":
            return ((PHOSPHATE, 1000.0 * float(self.phosphate_mol_l)),)
        if self.kind == "citrate_phosphate":
            return (
                (PHOSPHATE, 1000.0 * float(self.phosphate_mol_l)),
                (CITRATE, 1000.0 * float(self.citrate_mol_l)),
            )
        return ()

    def as_dict(self) -> Dict[str, Any]:
        return {
            "kind": self.kind,
            "phosphate_mol_L": self.phosphate_mol_l,
            "citrate_mol_L": self.citrate_mol_l,
            "declared_by_source": self.declared,
            "source": self.source,
        }


#: THE DECLARED DEFAULT. A caller who supplies no buffer gets an UNBUFFERED
#: pot and a warning that says so, rather than a silent clamp.
DEFAULT_BUFFER = BufferSpec(
    kind="none", declared=False,
    source="DEFAULT: no buffer was declared for this system.")

BUFFER_ABSENT_WARNING = (
    "no buffer was declared for this system, so the pH TRAJECTORY is "
    "EXTRAPOLATED: it is computed from water autoprotolysis and the charged "
    "solutes alone. If the experiment was in fact buffered, every pH-dependent "
    "rate in this run is wrong in the direction of too much drift."
)


# ===========================================================================
# 3. THE TWO CALIBRATED CONSTANTS
# ===========================================================================
# Both are calibrated ONLY on FIT_HOLDOUT_DECLARATION.md Amendment 7's three
# declared drift anchors. Their VALUES are written by the B2.2 fit and frozen
# into the fit report; the numbers below are the SEARCH BOUNDS and the
# placeholder defaults used when no fitted value is supplied.

ACID_YIELD_BOUNDS: Tuple[float, float] = (0.0, 1.0)
ARP_AMINE_PKA_BOUNDS: Tuple[float, float] = (5.0, 11.0)

#: Placeholders. `sulfur.integrate_sulfur` REQUIRES an explicit PhDrift when a
#: dynamic trajectory is asked for, so these are never silently operative in a
#: scored run; they exist so that the module can be exercised standalone.
ACID_YIELD_PLACEHOLDER = 0.5
ARP_AMINE_PKA_PLACEHOLDER = 8.0


@dataclass(frozen=True)
class PhDrift:
    """The two calibrated drift constants, carried together and immutably."""

    acid_yield: float = ACID_YIELD_PLACEHOLDER
    arp_amine_pka: float = ARP_AMINE_PKA_PLACEHOLDER

    def __post_init__(self) -> None:
        lo, hi = ACID_YIELD_BOUNDS
        if not lo <= float(self.acid_yield) <= hi:
            raise ValueError(f"acid_yield {self.acid_yield} outside {ACID_YIELD_BOUNDS}")
        lo, hi = ARP_AMINE_PKA_BOUNDS
        if not lo <= float(self.arp_amine_pka) <= hi:
            raise ValueError(
                f"arp_amine_pka {self.arp_amine_pka} outside {ARP_AMINE_PKA_BOUNDS}"
            )

    def arp_group(self) -> TitratableGroup:
        return TitratableGroup(
            "arp", "pentose Amadori compound (COOH / secondary NH2+)",
            (ARP_CARBOXYL_PKA, float(self.arp_amine_pka)),
            (ARP_CARBOXYL_DH, 45.0), 1.0,
            "carboxyl pKa 3.00 (alpha-carboxyl of an N-substituted amino acid, "
            "[CRC] class value, dH_ion ~ 0); secondary-ammonium pKa CALIBRATED "
            "on Amendment 7's three Zhou drift endpoints, dH_ion 45 kJ/mol "
            "[GKL, secondary-amine class]")

    def as_dict(self) -> Dict[str, float]:
        return {
            "acid_yield_per_sink_event": float(self.acid_yield),
            "arp_secondary_ammonium_pKa": float(self.arp_amine_pka),
        }


# ===========================================================================
# 4. THE CHARGE BALANCE
# ===========================================================================
#: Which state species carries which titratable group, and how many equivalents
#: of it. ``ACID`` is the network's lumped organic-acid pool; ``FA`` and ``AA``
#: are the trunk's own measured formic and acetic pools (Martins 2005 steps 5,
#: 8 and 10), so a hexose run gets its acid load from MEASURED kinetics and a
#: pentose run from the lumped sink.
_ACID_POOL_KEYS = ("ACID",)


def titratable_inventory(
    concentrations: Mapping[str, float], drift: PhDrift, buffer_spec: BufferSpec
) -> Tuple[Tuple[TitratableGroup, float], ...]:
    """
    Every titratable centre present, as ``(group, mmol/L)``.

    The organic-acid pools are scaled by the calibrated acid yield; the solutes
    and the buffer are not scaled by anything.
    """
    out: list = []
    yield_ = float(drift.acid_yield)

    # -- the acid the network made ---------------------------------------
    lumped = yield_ * float(concentrations.get("ACID", 0.0))
    for group, share in ORGANIC_ACID_MIX:
        if lumped > 0.0:
            out.append((group, share * lumped))
    # the trunk's own MEASURED acids enter at full strength: they are named
    # molecules with measured formation kinetics, not a lump.
    formic = float(concentrations.get("FA", 0.0))
    acetic = float(concentrations.get("AA", 0.0))
    if formic > 0.0:
        out.append((FORMIC, formic))
    if acetic > 0.0:
        out.append((ACETIC, acetic))

    # -- the charged solutes ---------------------------------------------
    cys = float(concentrations.get("Cys", 0.0))
    if cys > 0.0:
        out.append((CYSTEINE, cys))
    arp = float(concentrations.get("ARP", 0.0))
    if arp > 0.0:
        out.append((drift.arp_group(), arp))
    ttca = float(concentrations.get("TTCA", 0.0))
    if ttca > 0.0:
        out.append((TTCA_GROUP, ttca))

    # -- the buffer -------------------------------------------------------
    out.extend(buffer_spec.groups())
    return tuple(out)


def _charge_residual(
    ph: float,
    temperature_k: float,
    sid_mmol_l: float,
    inventory: Sequence[Tuple[TitratableGroup, float]],
) -> float:
    """
    ``SID + [H+] - [OH-] + SUM C_i z_i``, in mmol/L. Strictly DECREASING in pH,
    which is what makes the bisection below unconditionally safe.
    """
    h = 1000.0 * 10.0 ** (-float(ph))
    oh = 1000.0 * 10.0 ** (float(ph) - pkw_at(temperature_k))
    total = float(sid_mmol_l) + h - oh
    for group, concentration in inventory:
        if concentration <= 0.0:
            continue
        total += concentration * group.average_charge(ph, temperature_k)
    return total


def solve_ph(
    temperature_k: float,
    sid_mmol_l: float,
    inventory: Sequence[Tuple[TitratableGroup, float]],
    *,
    lo: float = -1.0,
    hi: float = 15.0,
    tol: float = 1e-9,
    max_iter: int = 200,
) -> float:
    """Bisection on the (monotone) charge residual. Returns pH."""
    f_lo = _charge_residual(lo, temperature_k, sid_mmol_l, inventory)
    f_hi = _charge_residual(hi, temperature_k, sid_mmol_l, inventory)
    if f_lo < 0.0:
        return lo
    if f_hi > 0.0:
        return hi
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        if hi - lo < tol:
            return mid
        if _charge_residual(mid, temperature_k, sid_mmol_l, inventory) > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


#: The temperature a pH LABEL is read at. Every pH number in this corpus is a
#: cooled-sample reading; none is an in-situ measurement.
LABEL_TEMPERATURE_K = 298.15


def strong_ion_difference(
    initial_ph_label: float,
    initial_concentrations: Mapping[str, float],
    drift: PhDrift,
    buffer_spec: BufferSpec,
    *,
    label_temperature_k: float = LABEL_TEMPERATURE_K,
) -> float:
    """
    The CONSERVED strong-ion difference, mmol/L, implied by the declared
    initial pH LABEL.

    This is the model's single most important quantity and it is not fitted:
    it is arithmetic on the declared initial pH and the declared charge, done
    AT THE TEMPERATURE THE LABEL WAS READ AT. It is what encodes "this pot was
    taken to pH 8 with NaOH and that pot was not", and it is the only place the
    three Zhou runs differ from one another at t = 0.
    """
    inventory = titratable_inventory(initial_concentrations, drift, buffer_spec)
    residual_without_sid = _charge_residual(
        initial_ph_label, float(label_temperature_k), 0.0, inventory
    )
    return -residual_without_sid


def ph_of_state(
    concentrations: Mapping[str, float],
    temperature_k: float,
    sid_mmol_l: float,
    drift: PhDrift,
    buffer_spec: BufferSpec,
) -> float:
    """The pH of one state, at one temperature, at the conserved SID."""
    inventory = titratable_inventory(concentrations, drift, buffer_spec)
    return solve_ph(temperature_k, sid_mmol_l, inventory)


def buffer_capacity_mmol_per_ph(
    concentrations: Mapping[str, float],
    temperature_k: float,
    ph: float,
    drift: PhDrift,
    buffer_spec: BufferSpec,
    *,
    delta: float = 0.01,
) -> float:
    """
    The van Slyke buffer capacity, mmol/L of strong base per pH unit, computed
    by finite difference on the charge balance. Reported, never fitted.
    """
    inventory = titratable_inventory(concentrations, drift, buffer_spec)
    up = _charge_residual(ph + delta, temperature_k, 0.0, inventory)
    down = _charge_residual(ph - delta, temperature_k, 0.0, inventory)
    return float((down - up) / (2.0 * delta))


PH_STATE_PROVENANCE: Dict[str, Any] = {
    "module": "B2.2",
    "fitted_constants": 2,
    "fitted_constant_names": [
        "acid_yield_per_sink_event", "arp_secondary_ammonium_pKa",
    ],
    "calibration_anchors": (
        "FIT_HOLDOUT_DECLARATION.md Amendment 7: Zhou 2023 Fig. 2, ARP + Cys, "
        "unbuffered, 120 C, initial pH 6/7/8 -> final 3.22 / 3.42 / 5.07. "
        "THREE anchors, TWO constants."
    ),
    "textbook_constants_not_fitted": [
        "pKw(T) tabulated 0-200 C [Marshall & Franck 1981]",
        "formic 3.751 / -0.13", "acetic 4.756 / -0.41",
        "cysteine 1.92, 8.33, 10.28 / 0, 26, 43",
        "phosphate 2.148, 7.198, 12.35 / -7.9, +3.6, +16.0",
        "citrate 3.128, 4.761, 6.396 / +4.07, +2.23, -3.38",
        "thiazolidine-4-carboxylic acid 1.80, 6.24 / 0, 40",
    ],
    "two_ph_scales": (
        "Every pH number in this corpus is a COOLED-SAMPLE reading. The label "
        "is read at 25 C and used once, to set the conserved SID; the IN-SITU "
        "pH that every rate constant sees is solved from that SID at reaction "
        "temperature; the COOLED pH that a published endpoint is compared "
        "against is solved from the same SID at 25 C. No fitted quantity "
        "converts between them. B2 and B2.1 used ONE scale, which put a factor "
        "of ~8 into the thiolate inventory of an unbuffered pot at a cold "
        "pH 7 (cysteine's thiol pKa falls 8.33 -> 7.23 between 25 and 120 C)."
    ),
    "declared_approximations": [
        "The lumped organic-acid pool is carried as a 1:1 formic/acetic "
        "mixture -- Martins' two named acids. A pure-formic or pure-acetic "
        "pool moves the predicted endpoints by ~0.4 pH units.",
        "TTCA's pKa are those of the unsubstituted parent acid.",
        "The ARP carboxyl pKa is a class value (3.00), not a measurement.",
    ],
    "conserved_quantity": (
        "the strong-ion difference. No reaction in the network creates or "
        "destroys a strong ion, so SID is computed once from the declared "
        "initial pH and is invariant along the trajectory -- which is the "
        "charge-conservation property tests/unit/test_kinetic_core_b2_2.py "
        "checks directly."
    ),
}
