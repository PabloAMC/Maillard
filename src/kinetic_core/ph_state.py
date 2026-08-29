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
      * ``"acetate"`` -- total acetate (acetic acid titrated with a strong
                           base), mol/L. B2.3, and it exists because a real
                           bundle needs it: Chang 2021's scored arm is 1%
                           acetic acid back-titrated to pH 6.0, which IS a
                           buffer and must not be recorded as water. It costs
                           no new constant -- the ``ACETIC`` group above is
                           already in this module, from the CRC, and is
                           already half of ``ORGANIC_ACID_MIX``.
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
    acetate_mol_l: float = 0.0
    declared: bool = False
    source: str = ""

    _KINDS = ("none", "phosphate", "citrate_phosphate", "acetate", "clamped")

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
        if self.kind == "acetate" and self.acetate_mol_l <= 0.0:
            raise ValueError("an acetate buffer needs a positive molarity")

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
        if self.kind == "acetate":
            return ((ACETIC, 1000.0 * float(self.acetate_mol_l)),)
        return ()

    def as_dict(self) -> Dict[str, Any]:
        return {
            "kind": self.kind,
            "phosphate_mol_L": self.phosphate_mol_l,
            "citrate_mol_L": self.citrate_mol_l,
            "acetate_mol_L": self.acetate_mol_l,
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

    # -- B2.3: the carboxyl centres the sinks CARRY OUT of a consumed --------
    #    molecule. FULL STRENGTH: a carboxyl carried out of a molecule that
    #    demonstrably had one is one equivalent, not a fraction of anything.
    #    Only `ACID` -- the lumped deoxyosone sink, whose acid yield nobody
    #    measured -- is scaled by the fitted yield.
    carried_acid = float(concentrations.get("CBX", 0.0))
    if carried_acid > 0.0:
        for group, share in ORGANIC_ACID_MIX:
            out.append((group, share * carried_acid))

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


# ===========================================================================
# 5. B2.3 -- THE NET-SOLUTE-CHARGE AUDIT
# ===========================================================================
# WHAT WENT WRONG, IN ONE SENTENCE. The charge balance conserves the STRONG ION
# DIFFERENCE, which is right -- no reaction makes or destroys sodium -- but the
# NETWORK was free to delete a titratable solute without saying where its
# titratable centre went, and the sodium that had titrated it was then left
# balancing nothing, so the pot drifted alkaline for a bookkeeping reason.
# B2.2's diagnosis sec. 3 measured it: Kang's 120 C pot at a predicted pH 11.4
# against a measured 4.9.
#
# WHY A "CONSERVATION" TEST IS NOT ENOUGH, AND WHAT IS ENFORCED INSTEAD.
# Titratable centres are NOT conserved quantities and a checker that demanded
# they were would be wrong in both directions:
#
#   * MAKING one is ordinary chemistry. `r_tdg_fa` turns a neutral deoxyosone
#     into neutral formic acid; the molecule that appears is uncharged and its
#     dissociation is the equilibrium's business, not the stoichiometry's.
#   * DESTROYING one is sometimes real too. `r_cys_actz` buries the cysteine
#     nitrogen in a thiazole ring whose conjugate acid sits near pKa 2.5. That
#     centre is genuinely gone, and the step genuinely releases a proton.
#
# So what is enforced is DECLARATION, not conservation: every reaction whose
# net titratable-centre count is non-zero must appear in the network's
# ``CENTRE_LEDGER`` with the exact delta and a written basis, and every ledger
# entry must match the stoichiometry it claims to describe. A step cannot
# quietly delete a carboxylate ever again -- it either carries it, or it says
# in the ledger that it destroyed it and why.
#
# THE TABLE IS DERIVED FROM `titratable_inventory` ABOVE, not written twice:
# a species has a centre here if and only if the charge balance gives it one,
# and `tests/unit/test_kinetic_core_b2_3.py` pins that correspondence. A
# species the charge balance cannot see contributes NOTHING here, which is why
# `UNTRACKED_TITRATABLE` exists directly below and is itself pinned.

CENTRE_KINDS: Tuple[str, ...] = ("carboxyl", "amine")

#: Titratable centres per unit of each species, AS THE CHARGE BALANCE MODELS
#: THEM. One entry per group `titratable_inventory` can put into the charge
#: balance, and nothing else.
TITRATABLE_CENTRES: Mapping[str, Mapping[str, int]] = {
    "Cys": {"carboxyl": 1, "amine": 1},   # CYSTEINE: pKa 1.92 / 8.33 / 10.28
    "ARP": {"carboxyl": 1, "amine": 1},   # drift.arp_group(): 3.00 / calibrated
    "TTCA": {"carboxyl": 1, "amine": 1},  # TTCA_GROUP: 1.80 / 6.24 (ring NH+)
    "FA": {"carboxyl": 1},                # FORMIC
    "AA": {"carboxyl": 1},                # ACETIC
    "ACID": {"carboxyl": 1},              # ORGANIC_ACID_MIX, scaled by the yield
    "CBX": {"carboxyl": 1},               # ORGANIC_ACID_MIX, full strength
    #: B3 lane only. S-(2-carbamoylethyl)cysteine IS an alpha-amino acid: the
    #: acrylamide adds across the thiol and leaves the carboxyl and the
    #: ammonium untouched. It is listed so that `a_acr_cys` nets to zero
    #: HONESTLY -- because the centres really do survive into the adduct -- and
    #: not because the acrylamide lane happens to have no pH state. The lane's
    #: charge balance is not wired up, so this entry moves no prediction today.
    "ACRCYS": {"carboxyl": 1, "amine": 1},
}

#: SPECIES WITH REAL TITRATABLE GROUPS THAT THE CHARGE BALANCE DELIBERATELY
#: DOES NOT MODEL. This is a DECLARED GAP LIST, not an exemption list: nothing
#: in it is allowed to satisfy a ledger entry, and the unit test pins the list
#: verbatim so that a new titratable species cannot join the network silently.
#: Each entry says why the omission is tolerable, or admits that it is not.
UNTRACKED_TITRATABLE: Mapping[str, str] = {
    "H2S": (
        "sulfide, pKa1 7.05 at 25 C, present at mmol/L. THE LARGEST OF THESE "
        "GAPS AND THE ONLY ONE THAT IS NOT OBVIOUSLY SMALL. It is left out "
        "because closing it is a MODELLING ADDITION (a new titratable solute "
        "in the balance), not a conservation fix, and B2.3's licence is the "
        "conservation fix. Sized, not guessed: at the acidic endpoints these "
        "pots reach, sulfide is fully protonated and carries no charge, so "
        "the omission bites at the START of an alkaline run and not at its "
        "end. It is reported as a residual defect, not as closed."
    ),
    "Gly": (
        "glycine, the B1 trunk's amine. The trunk lane carries NO pH state at "
        "all -- `integrate_sulfur`'s dynamic mode is sulfur-lane only -- so "
        "the trunk's zwitterion is invisible to the charge balance from both "
        "ends and its steps net to zero here by construction. If a pH state "
        "is ever given to the trunk, this entry becomes a defect."
    ),
    "SB": (
        "the Schiff base still carries glycine's carboxyl and its imine "
        "nitrogen; both are invisible to the balance for the same trunk-lane "
        "reason as Gly, and become defects the day the trunk gets a pH state."
    ),
    "AMA": (
        "the Amadori compound still carries glycine's carboxyl and its "
        "secondary ammonium; same trunk-lane reason as Gly. NOTE the "
        "asymmetry with the SULFUR lane's ARP, which IS in the balance and "
        "whose ammonium pKa is a calibrated constant."
    ),
    "MEL_N": (
        "melanoidins are polyanionic and retain amine nitrogen, so the "
        "terminal polymer really does hold titratable capacity. It is carried "
        "ELEMENTALLY (mmol of atom, not of molecule), so no per-unit centre "
        "count is definable. Same trunk-lane exemption as Gly."
    ),
    "THI": (
        "thiamine. Its thiazolium is a PERMANENT CATION and its pyrimidinium "
        "sits near pKa 4.8, so a thiamine pot carries a charge the balance "
        "does not see. THE SIZE IS NOT SMALL (Zhang's pot charges 44.5 "
        "mmol/L) and this is recorded as a genuine open defect, not as a "
        "rounding argument."
    ),
    "MFT": (
        "a furan-3-thiol, pKa near 7. Present at ug/L, i.e. <= 1e-3 mmol/L "
        "against an acid pool of tens of mmol/L -- four to five orders below "
        "the charge scale. Negligible on arithmetic, not on hope."
    ),
    "FFT": (
        "a furfuryl thiol, same pKa class as MFT and present at the same ug/L "
        "scale. Same arithmetic argument, same conclusion: four to five orders "
        "below the charge scale, negligible by size and not by hope."
    ),
    "MESH": (
        "methanethiol, pKa 10.3 and therefore barely dissociated anywhere "
        "these pots go, at a concentration lower still than MFT's. It is the "
        "least consequential entry on this list and it is here so that the "
        "list is exhaustive rather than convenient."
    ),
}


#: THE STRONGEST SINGLE ASSUMPTION IN THE pH MODEL, WRITTEN ONCE AND SHARED BY
#: EVERY LEDGER ENTRY THAT RELIES ON IT.
#:
#: B2.3 CARRIES THE CARBOXYL AND DECLARES THE AMINE DESTROYED, and that
#: asymmetry is not a convenience -- it is what the corpus's own DECLARED FIT
#: anchors force. Zhou 2023's system (Amendment 7) charges 20 mmol/L of Amadori
#: compound and 20 mmol/L of cysteine, i.e. 40 mmol/L of amino nitrogen, into
#: an unbuffered pot, and finishes at a MEASURED pH of 3.22 / 3.42. Forty
#: mmol/L of surviving ammonium (pKa 9.25, hence essentially fully protonated
#: anywhere below pH 7) makes those endpoints arithmetically unreachable at any
#: acid load this network can produce: the model was measured predicting
#: 5.04 / 5.11 with the amine carried, against 2.94 / 3.02 with it destroyed,
#: at B2.2's OWN frozen parameters (`scratch/b23_encoding_probe.py`, run before
#: this choice was made and before any refit).
#:
#: The chemistry says the same thing and says it more generally. Maillard
#: systems ACIDIFY, and the textbook reason is exactly this: the amino group is
#: consumed -- into pyrazines, into other N-heterocycles, into melanoidin
#: nitrogen, none of which is a base in this window -- while carboxylic acids
#: accumulate. A model in which liberated amino nitrogen persists as free
#: ammonia cannot reproduce the single most robust process observable in the
#: whole corpus.
#:
#: WHAT THIS IS NOT. It is not a claim that no ammonia is ever released --
#: Nedvidek 1992 Scheme 3 writes NH3 explicitly, and it is surely released. It
#: is a claim that it does not SURVIVE as a titratable base on the timescale of
#: these holds, and that the module, having no step for its consumption, must
#: account for the centre's loss at the point of release rather than pretend
#: the nitrogen is inert. This is the claim `FRAG_N` has silently encoded since
#: B1; B2.3's contribution is to make it explicit, name its evidence, and stop
#: it from concealing the CARBOXYL half, which does not vanish and which B2.2
#: was measured deleting.
#:
#: IF IT IS WRONG, IT IS WRONG IN A NAMED DIRECTION: the pots would be modelled
#: too acidic, most in the systems with the highest amino-nitrogen loading
#: (Zhou's ARP + Cys, Cerny's ternary), and least in Kang's TTCA pot.
#: WAVE B2.4 -- THE CITATION IS CORRECTED, AND THE EVIDENCE RE-DERIVED.
#:
#: Through B2.3 the last sentence of this declaration read "See
#: ph_state.AMINE_FATE_BASIS for ... the pre-freeze probe that sized it". Two
#: defects: the sentence sits INSIDE `AMINE_FATE_BASIS` and so pointed at
#: itself, and the probe it pointed at -- `scratch/b23_encoding_probe.py` --
#: DOES NOT RUN on the shipped tree. D1 sec. 7 found it raises `KeyError:
#: 'AMN'`: there is no amine pool in `src/kinetic_core/` to patch, so two of
#: that probe's three "encodings" collapsed onto one code path and the axis
#: could not be re-probed at all.
#:
#: B2.4 rebuilt the probe against the CURRENT species set --
#: `scripts/generators/probe_amine_fate_b2_4.py`, which derives the released
#: amino nitrogen as (Cys + ARP + TTCA at t=0) minus (Cys + ARP + TTCA at t)
#: and adds it back as ammonium at pKa 9.25, with no new species and no new
#: parameter. IT REPRODUCES THE B2.3 PRE-REGISTRATION'S PUBLISHED TABLE IN
#: EVERY CELL, TO TWO DECIMALS. So the defect was in the script, not in the
#: evidence, and the declaration keeps its basis rather than being weakened.
#:
#: The cell-by-cell reproduction check is NOT reprinted here on purpose -- a
#: number copied into a source comment is a number that can go stale silently.
#: It lives in the frozen artifact
#: `results/validation/kinetic_core_b2_4_amine_fate_probe.json`, under
#: `verdict.reproduction_cell_by_cell`. The one comparison this declaration
#: rests on is stated in the declaration itself, below.
#:
#: On the two anchors that actually acidify, the shipped encoding's mean |miss|
#: is 0.34 pH units and "carry both" is 1.76: the amine cannot be carried and
#: still reach those endpoints.
AMINE_FATE_BASIS = (
    "THE LIBERATED AMINO NITROGEN IS DECLARED NON-TITRATABLE AT THE POINT OF "
    "RELEASE. Basis: (i) Amendment 7's three DECLARED FIT drift endpoints -- "
    "40 mmol/L of amino nitrogen charged, measured finish at pH 3.22 / 3.42 -- "
    "which 40 mmol/L of surviving ammonium (pKa 9.25) makes arithmetically "
    "unreachable; (ii) the textbook reason Maillard systems acidify at all, "
    "which is that the amino group is consumed into pyrazines, N-heterocycles "
    "and melanoidin nitrogen while carboxylic acids accumulate. The module has "
    "no step for that consumption, so the centre's loss is accounted where the "
    "nitrogen leaves. (iii) A RUNNABLE probe, "
    "scripts/generators/probe_amine_fate_b2_4.py, which realises the "
    "counterfactual encoding against the CURRENT species set and puts Zhou's "
    "three cooled endpoints at 5.04 / 5.11 / 5.89 when the amine is carried, "
    "against 2.94 / 3.02 / 3.47 when it is not. The comment block above this "
    "declaration carries the full statement, that probe's frozen result, and "
    "the direction of the error if the declaration is wrong."
)


def centre_delta(
    reactants: Mapping[str, int], products: Mapping[str, int]
) -> Dict[str, int]:
    """
    Net change in each titratable centre across one step, products minus
    reactants. Zero for every kind means the step moves no titratable capacity.
    """
    out = {kind: 0 for kind in CENTRE_KINDS}
    for mapping, sign in ((reactants, -1), (products, +1)):
        for key, coefficient in mapping.items():
            for kind, count in TITRATABLE_CENTRES.get(key, {}).items():
                out[kind] += sign * int(coefficient) * int(count)
    return out


def validate_charge_closure(reactions, ledger: Mapping[str, Mapping[str, Any]]) -> None:
    """
    Raise unless every step's titratable-centre movement is DECLARED.

    Three failure modes, and each is its own message because they mean
    different things:

      1. a step moves a centre and is NOT in the ledger -- the B2.2 defect,
         and the one this function exists to make impossible;
      2. a step is in the ledger with the wrong delta -- the declaration and
         the stoichiometry have drifted apart, which is how a fix rots;
      3. a ledger entry names a step that no longer moves anything -- dead
         paperwork, which is how a ledger stops being read.
    """
    seen = set()
    for reaction in reactions:
        delta = centre_delta(reaction.reactants, reaction.products)
        moved = {k: v for k, v in delta.items() if v != 0}
        entry = ledger.get(reaction.key)
        if not moved:
            if entry is not None:
                raise ValueError(
                    f"{reaction.key}: CENTRE_LEDGER declares a movement "
                    f"{ {k: v for k, v in entry.items() if k in CENTRE_KINDS} } "
                    f"but the stoichiometry moves no titratable centre. Remove "
                    f"the stale entry; a ledger nobody can trust is worse than "
                    f"no ledger."
                )
            continue
        seen.add(reaction.key)
        if entry is None:
            raise ValueError(
                f"{reaction.key}: moves titratable centres {moved} and is NOT "
                f"declared in CENTRE_LEDGER. A step that deletes a carboxylate "
                f"or an amine without saying where it went makes the strong "
                f"ion that titrated it balance nothing, and the pot drifts for "
                f"a bookkeeping reason -- this is exactly the defect Build Wave "
                f"B2.2's diagnosis sec. 3 measured (Kang's pot at a predicted "
                f"pH 11.4 against a measured 4.9). Either CARRY the centre "
                f"into CBX, or declare the destruction and its basis."
            )
        declared = {k: int(entry.get(k, 0)) for k in CENTRE_KINDS}
        if declared != delta:
            raise ValueError(
                f"{reaction.key}: CENTRE_LEDGER declares {declared} but the "
                f"stoichiometry gives {delta}."
            )
        if not str(entry.get("basis", "")).strip():
            raise ValueError(
                f"{reaction.key}: its CENTRE_LEDGER entry has no `basis`. Every "
                f"declared centre movement must say, in words, what happens to "
                f"the group."
            )
    stale = sorted(set(ledger) - seen)
    if stale:
        raise ValueError(
            f"CENTRE_LEDGER declares steps that move nothing or do not exist: "
            f"{stale}"
        )


def net_solute_charge(
    concentrations: Mapping[str, float],
    ph: float,
    temperature_k: float,
    drift: "PhDrift",
    buffer_spec: "BufferSpec",
) -> float:
    """
    The net charge carried by every modelled solute at ``ph``, mmol/L.

    The reporting half of the audit: ``validate_charge_closure`` polices the
    NETWORK at import, this polices a TRAJECTORY at run time. Along a run the
    quantity ``SID + net_solute_charge + [H+] - [OH-]`` must stay at zero,
    which is what `_charge_residual` returns.
    """
    inventory = titratable_inventory(concentrations, drift, buffer_spec)
    return float(sum(
        c * g.average_charge(float(ph), float(temperature_k))
        for g, c in inventory if c > 0.0
    ))


PH_STATE_PROVENANCE: Dict[str, Any] = {
    "module": "B2.3",
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
        "B2.3: CBX carries the mixture's 4.25 for a heterogeneous set of real "
        "leaving groups (pyruvic 2.49, CO2/carbonic 6.35, acetic 4.76). No CO2 "
        "partition model exists in this corpus, so the mixture is declared "
        "rather than a Henry constant guessed.",
        "B2.3: the LIBERATED AMINO NITROGEN is declared NON-TITRATABLE at "
        "the point of release -- see AMINE_FATE_BASIS. That is a chemistry "
        "claim with a measured basis, and it is the strongest single "
        "assumption in the pH model.",
    ],
    "conserved_quantity": (
        "the strong-ion difference. No reaction in the network creates or "
        "destroys a strong ion, so SID is computed once from the declared "
        "initial pH and is invariant along the trajectory -- which is the "
        "charge-conservation property tests/unit/test_kinetic_core_b2_2.py "
        "checks directly."
    ),
    "b2_3_declared_invariant": (
        "DECLARATION, not conservation. Titratable centres are not conserved "
        "quantities -- making a neutral acid out of a neutral sugar is "
        "ordinary chemistry and burying an amine in a thiazole ring really "
        "does release a proton -- so `validate_charge_closure` requires every "
        "step whose net centre count moves to be declared in the network's "
        "CENTRE_LEDGER, with the exact delta and a written basis, and refuses "
        "the import otherwise. The B2.2 defect (a sink deleting a carboxylate "
        "in silence) is not merely fixed but unreachable."
    ),
    "b2_3_untracked_titratable_species": sorted(UNTRACKED_TITRATABLE),
    "b2_3_amine_fate": AMINE_FATE_BASIS,
}
