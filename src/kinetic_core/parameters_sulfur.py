"""
src/kinetic_core/parameters_sulfur.py

THE PARAMETER REGISTRY OF THE SULFUR MODULE (Build Wave B2).
=============================================================

Modules 1 (sulfur formation) and 2 (thiol consumption) of
``docs/reference/FIT_HOLDOUT_DECLARATION.md``. Same contract as B1's
``parameters.py``: every number carries its value, unit, source anchor down to
the table, the conditions it was measured under, its ``ph_of_measurement``, an
``evidence_class`` and a ``rate_transfer`` licence, and all of it travels into
runtime metadata.

FIVE STANDING POLICIES, ENFORCED HERE
-------------------------------------
1. **NO DFT.** Same as B1. ``assert_no_dft_sulfur()`` runs at import. Nothing
   from ``data/qm/`` is read, imported or referenced.

2. **NO SINGLE ARRHENIUS FOR THIOL CONSUMPTION.** K1 declared the gap and K3
   turned it into a structural finding: four papers at four temperatures
   measure four DIFFERENT dominant sinks, and each excludes the others'
   (inventory sec. B.7). Every consumption channel here is therefore a NAMED
   CHANNEL with its own declared temperature window and ``ea_kj_mol = None``.
   Pairing two channels' rates to extract an Ea is recorded in
   ``PROHIBITED_DERIVATIONS`` and there is no code path that can do it: the
   rate evaluator refuses to Arrhenius-scale a parameter whose Ea is None and
   emits an extrapolation warning instead.

3. **NO FIXED BRANCH FRACTIONS.** There is not one branch-fraction constant in
   this file, and there is no place to put one: the registry contains only rate
   constants, and every observable split (thiamine:sugar MFT, dimer:monomer,
   NF:intact-skeleton, FFT:furfural) is a RATIO OF MASS-ACTION FLUXES computed
   at run time. Cerny 2007 Table 5 measures a 2x precursor change moving the
   xylose share of MFT from 15% to 46%; a fixed fraction fails that row by
   construction. ``tests/unit/test_kinetic_core_b2.py`` pins the property.

4. **pH IS STRUCTURAL, NOT FITTED.** B1 has no pH term because the trunk corpus
   licenses none. The sulfur corpus is different: it contains MEASURED pH
   dependence (Zheng & Ho's four-pH cysteine thermolysis set) and two
   acid-base equilibria whose constants are physical rather than fitted. Those
   three things are the ONLY pH dependence in this module. There is no fitted
   pH polynomial, no pH-shape parameter, and no per-pH rate constant.
   ``PH_TERM_PROVENANCE`` states exactly what enters and what does not.

5. **EVERY SPECIES HAS A CONSUMPTION TERM.** Where the corpus supplies no rate
   for one, the constant is fitted under bounds that ALLOW ~ZERO, so the data
   may reject the channel -- B1's pattern. Where the corpus supplies no rate
   AND no FIT row can constrain it, the constant is declared in
   ``NO_MEASURED_RATE`` and left at zero with a flag, rather than invented.

UNITS: mmol/L, minutes, Kelvin, kJ/mol.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, replace
from typing import Any, Dict, Mapping, Optional, Tuple

from .parameters import EVIDENCE_CLASSES, R_KJ

#: Reference temperature of the sulfur module: 145 C, the temperature of the
#: Hofmann 1998 SIDA anchors and of the Cerny isotope experiments -- i.e. the
#: temperature at which most of the FIT panel was measured. B1's reference is
#: 100 C; the two modules therefore quote k_ref at different temperatures and
#: each says so on every row.
T_REF_S_K = 418.15

UNIT_BY_ORDER: Mapping[int, str] = {
    1: "1/min",
    2: "L/(mmol*min)",
    3: "L^2/(mmol^2*min)",
}


@dataclass(frozen=True)
class SulfurParameter:
    """
    One sulfur rate constant, with everything needed to defend it.

    Differs from B1's ``KineticParameter`` in exactly three ways, each forced
    by the sulfur corpus:

      * ``order`` may be 3 (the oxidative dimerisation is second order in thiol
        and first order in oxidant);
      * ``ea_kj_mol`` may be None, meaning NO ACTIVATION ENERGY EXISTS. Such a
        parameter is held at its measured temperature and any evaluation off
        that temperature raises an extrapolation warning. This is the encoding
        of inventory sec. C.1 / sec. B.7;
      * ``channel`` names the mechanism, so that the four thiol sinks are four
        separate objects rather than four values of one lumped constant.
    """

    key: str
    transformation: str
    k_ref: Optional[float]
    ea_kj_mol: Optional[float]
    unit: str
    order: int
    evidence_class: str
    source_anchor: str
    dossier_anchor: str
    conditions: str
    ph_of_measurement: Optional[float]
    #: the temperature this k_ref is quoted AT, K
    t_ref_k: float
    temperature_range_c: Tuple[float, float]
    rate_transfer: str
    #: named mechanism, for the consumption channels; "" for formation steps
    channel: str = ""
    #: which pH factor multiplies this step: "", "acid", "base", "neutral_h2s"
    ph_factor: str = ""
    k_ref_ci95: Optional[float] = None
    ea_ci95_kj_mol: Optional[float] = None
    flags: Tuple[str, ...] = ()
    note: str = ""

    def __post_init__(self) -> None:
        if self.evidence_class not in EVIDENCE_CLASSES:
            raise ValueError(
                f"{self.key}: evidence_class {self.evidence_class!r} not in "
                f"{EVIDENCE_CLASSES}. Computational (DFT-derived) barriers are "
                f"refused by owner policy."
            )
        if self.order not in UNIT_BY_ORDER:
            raise ValueError(f"{self.key}: order must be 1, 2 or 3, got {self.order}")
        if self.unit != UNIT_BY_ORDER[self.order]:
            raise ValueError(
                f"{self.key}: order {self.order} requires unit "
                f"{UNIT_BY_ORDER[self.order]!r}, got {self.unit!r}"
            )
        if self.ph_factor not in ("", "acid", "base", "neutral_h2s"):
            raise ValueError(f"{self.key}: unknown ph_factor {self.ph_factor!r}")

    @property
    def has_activation_energy(self) -> bool:
        return self.ea_kj_mol is not None

    def k_at(self, temperature_k: float) -> float:
        """
        k(T). With an Ea: the reparameterised Arrhenius form about ``t_ref_k``.
        WITHOUT an Ea: the measured value, UNSCALED.

        The unscaled branch is deliberate and is the whole of policy 2. It is
        not a silent default: ``extrapolation_warnings`` in ``sulfur.py``
        reports every channel evaluated off its measured temperature, on every
        run, with the temperature gap in degrees.
        """
        if self.k_ref is None:
            raise ValueError(f"{self.key}: parameter is not numerically populated")
        if self.ea_kj_mol is None:
            return float(self.k_ref)
        return float(self.k_ref) * math.exp(
            -(float(self.ea_kj_mol) / R_KJ)
            * (1.0 / float(temperature_k) - 1.0 / float(self.t_ref_k))
        )

    def as_metadata(self) -> Dict[str, Any]:
        payload = asdict(self)
        payload["flags"] = list(self.flags)
        payload["temperature_range_c"] = list(self.temperature_range_c)
        payload["dft_derived"] = False
        payload["has_activation_energy"] = self.has_activation_energy
        if not self.has_activation_energy:
            payload["temperature_policy"] = (
                "HELD FIXED at the measured temperature. No activation energy "
                "exists for this step in the corpus; see PROHIBITED_DERIVATIONS."
            )
        return payload


# ===========================================================================
# 1. THE MEASURED pH DEPENDENCE -- the only pH content in this module
# ===========================================================================
# THE SHAPE IS NOT FITTED. It comes from three measured things and nothing
# else. This block is the answer to the build brief's requirement that pH enter
# STRUCTURALLY rather than as a fitted polynomial.

#: (i) Zheng & Ho 1994, "Kinetics of the Release of Hydrogen Sulfide from
#: Cysteine and Glutathione During Thermal Treatment", ACS Symp. Ser. 564
#: (Sulfur Compounds in Foods) 138-146, doi 10.1021/bk-1994-0564.ch011.
#: 0.1 M aqueous L-cysteine, pH 3/5/7/9, 80-110 C, FIRST-ORDER H2S release.
#:
#: THE (Ea, A) PAIRS ARE THE SOURCE-CONSISTENT MATCHED ONES, from the full
#: twice-verified transcription in
#: `data/lit/extraction_dossiers/zheng1994_extraction.md` sec. 5b: an ordinary
#: least-squares Arrhenius refit of the paper's OWN Table I, which reproduces
#: every printed Ea to +/-0.1 kcal/mol AND every printed R^2 exactly, and which
#: yields the pre-exponentials the paper never prints.
#:
#: THIS SUPERSEDES BOTH EARLIER VALUES IN THE REPOSITORY:
#:   * the shipped `arrhenius_params.yml` pair (Ea 130.4 kJ/mol, A 1.0e14 1/s)
#:     is ~15x TOO FAST against its own source at pH 7 / 100 C -- it returns
#:     5.6e-5 1/s where Table I measures 3.8e-6 -- because the Ea was averaged
#:     across pH while A was re-fitted free. Averaging one member of a
#:     correlated pair is exactly how that discrepancy arose.
#:   * the K3 inventory sec. F row 4 rounding (9.8e11/1.9e12/2.4e13/1.0e12 with
#:     Ea 131.0/133.1/134.7/123.0) is the same fit to fewer figures. The
#:     inventory's pH-7 Ea of 134.7 came from the ABSTRACT's 32.2 kcal/mol;
#:     Table V prints 32.3 and the refit gives 32.37, so the citation is
#:     re-pointed at Table V p. 145 and the value is 135.5 kJ/mol.
#:
#: pH 3.0 IS A SEPARATE MECHANISTIC CHANNEL AND IS NOT INTERPOLATED THROUGH.
#: At 100 C, K(pH 3.0) = 4.15e-5 EXCEEDS K(pH 5.0) = 2.65e-5, breaking the
#: monotone pH ordering that holds at 80, 90 and 110 C -- an inversion that is
#: present in the printed half-lives and visible in Figure 3, and that the
#: authors flag themselves: "The nonconformity of reactions at pH 3.0 may mean
#: that a different reaction mechanism is involved." Their own pH regressions
#: (Table II) are fitted over pH 5.0-9.0 ONLY, and the dossier's instruction is
#: explicit: "Do not extrapolate these eight equations to pH 3, and do not use
#: them below pH 5." This module therefore interpolates over pH 5-9 and uses the
#: pH-3.0 pair AS MEASURED below pH 5, with no blending across the mechanistic
#: boundary.
#:
#: CROSS-CHECK, and it passes: at pH 5 and 145 C this set gives k = 2.80e-3
#: /min, so 20 min consumes 5.4% of a 33 mM cysteine charge = 1.80 mmol/L H2S =
#: 0.18 mmol in 100 mL, against the K3 inventory's independently derived
#: "~0.17 mmol in 100 mL" (sec. A.3.5). Agreement to 6%. That derived number is
#: what makes the fed-precursor experiments a 6-12x H2S SURPLUS, which is the
#: quantitative mechanism behind the Hofmann/Cerny disjointness.
ZHENG_CYSTEINE_THERMOLYSIS: Dict[str, Any] = {
    "ph_points": (3.0, 5.0, 7.0, 9.0),
    "prefactor_per_second": (9.79e11, 1.93e12, 2.36e13, 1.04e12),
    "ea_kj_mol": (131.2, 133.0, 135.5, 123.3),
    "ea_kcal_mol_as_printed": (31.3, 31.8, 32.3, 29.4),
    "ea_kcal_mol_refit": (31.36, 31.78, 32.37, 29.47),
    "r_squared_printed": (0.940, 0.975, 0.942, 0.967),
    "order": 1,
    "temperature_range_c": (80.0, 110.0),
    "source_anchor": (
        "Zheng & Ho 1994, ACS Symp. Ser. 564:138-146, "
        "doi 10.1021/bk-1994-0564.ch011, Table I (16 rate constants) and "
        "Table V p. 145 (4 activation energies), aqueous 0.1 M L-cysteine"
    ),
    "dossier_anchor": (
        "data/lit/extraction_dossiers/zheng1994_extraction.md sec. 5b (the "
        "matched (Ea, A) refit pairs), sec. 3 (Table I), sec. 9 (the caveats), "
        "sec. 10a (the consequences for arrhenius_params.yml); and "
        "k3_final_parameter_inventory.md sec. A.3.5 for the derived in-situ "
        "H2S budget this set is cross-checked against"
    ),
    "prefactors_are_not_the_papers": (
        "K0 is never printed anywhere in Zheng & Ho. All four pre-exponentials "
        "are the dossier's least-squares refit of the paper's own Table I, "
        "labelled [Z] there and flagged here."
    ),
    "interpolation": (
        "log10(A) and Ea are interpolated LINEARLY IN pH over pH 5-9 ONLY -- "
        "the same window the authors fitted their own pH regressions over. "
        "Below pH 5 the pH-3.0 pair is used AS MEASURED, because the authors "
        "identify pH 3.0 as a DIFFERENT MECHANISM and their data invert the pH "
        "ordering there at 100 C. Nothing is blended across that boundary. "
        "Above pH 9 the pH-9 pair is clamped."
    ),
    "what_this_buys": (
        "H2S availability RISES ~17x from pH 5 to pH 9 at 100 C, MEASURED "
        "(Table I: 2.65e-5 -> 4.53e-4 1/min). That is inventory sec. B2.8's "
        "first half, and it is one of the two factors of the two-factor law."
    ),
    "caveats_that_apply_to_every_use": (
        "(1) INITIAL-RATE CONSTANTS ONLY -- maximum observed conversion in the "
        "entire study is ~12% (110 C, pH 9, 180 min), so the first-order form "
        "is untested past ~15% consumption, and this module drives it to far "
        "higher conversion at 145 C. "
        "(2) NO MAILLARD CO-REACTANT: neat thermolysis, no sugar, no "
        "dicarbonyl, no lipid, N2-flushed and sealed. These constants are an "
        "UPPER BOUND on the free-H2S source term in any real system, where H2S "
        "is simultaneously consumed by carbonyls. "
        "(3) NOMINAL ROOM-TEMPERATURE pH: phosphate shifts ~-0.3 to -0.5 units "
        "between 25 and 110 C and the paper applies no correction, so the "
        "effective hot pH at the 'pH 9' label is nearer 8.4-8.6. "
        "(4) 0.1 M SUBSTRATE, 10-100x above typical food free-cysteine levels, "
        "so any bimolecular channel (cystine formation) is exaggerated. "
        "(5) NO ERROR BARS ANYWHERE; R^2 = 0.940 on four points is a weak "
        "constraint, and the dossier estimates a standing +/-2-4 kcal/mol on "
        "the Ea that the paper never quotes. Any sigma this repo attaches is "
        "the repo's own assumption. "
        "(6) THE ELECTRODE MEASURES TOTAL SULFIDE (S2- + HS- + H2S), not "
        "molecular H2S -- which is why this module's H2S state variable is the "
        "TOTAL sulfide pool and the neutral-H2S speciation is applied at "
        "reaction time rather than being folded into the source term. "
        "(7) NOTHING ABOVE 110 C, and every system in this module runs at "
        "115-145 C. "
        "(8) H2S RELEASE ONLY -- the paper contains no consumption kinetics."
    ),
}

#: Zheng & Ho Table II, the authors' own explicit pH law, carried as a
#: CROSS-CHECK and NOT operative: K = slope * log10[OH-] + intercept, K in
#: 1/min, fitted over pH 5.0-9.0 only. The 110 C intercept is the DOSSIER'S
#: CORRECTION -- the paper prints 1.55e-4, which returns a NEGATIVE rate
#: constant at pH 9 and is wrong by one decade; the three-point refit of Table I
#: gives 1.5509e-3 while reproducing the printed slope exactly
#: (zheng1994_extraction.md sec. 6a, "PRINTING ERROR #1").
#:
#: It is not operative because it is isothermal at four temperatures and
#: carries no barrier, so it cannot be evaluated at this module's 115-145 C.
#: The Arrhenius pairs above can, and the two agree by construction inside
#: 80-110 C.
ZHENG_TABLE_II_PH_LAW: Dict[str, Any] = {
    "form": "K(1/min) = slope * log10([OH-]) + intercept, pH 5.0-9.0 only",
    "80C": {"slope": 6.87e-6, "intercept": 6.19e-5, "r_squared": 0.939},
    "90C": {"slope": 3.01e-5, "intercept": 2.78e-4, "r_squared": 0.985},
    "100C": {"slope": 1.07e-4, "intercept": 9.88e-4, "r_squared": 0.999},
    "110C": {"slope": 1.66e-4, "intercept": 1.55e-3, "r_squared": 0.982,
             "intercept_as_printed": 1.55e-4,
             "correction": "printed exponent wrong by one decade; the printed "
                           "equation returns a NEGATIVE K at pH 9"},
    "operative": False,
    "do_not_use_below_ph_5": True,
}

#: Zheng & Ho Table I, the measured constants themselves, 1/min. Carried so the
#: Arrhenius pairs above can be checked against the data they were fitted to
#: rather than against each other.
ZHENG_TABLE_I_K_PER_MIN: Dict[float, Dict[float, float]] = {
    80.0: {3.0: 1.87e-6, 5.0: 2.02e-6, 7.0: 9.71e-6, 9.0: 2.95e-5},
    90.0: {3.0: 8.30e-6, 5.0: 1.24e-5, 7.0: 5.99e-5, 9.0: 1.33e-4},
    100.0: {3.0: 4.15e-5, 5.0: 2.65e-5, 7.0: 2.30e-4, 9.0: 4.53e-4},
    110.0: {3.0: 5.23e-5, 5.0: 7.97e-5, 7.0: 3.35e-4, 9.0: 7.45e-4},
}

#: Below this pH the pH-3.0 pair is used as measured: the authors identify
#: pH 3.0 as a different mechanism and fit their own pH law over 5.0-9.0 only.
ZHENG_MECHANISM_BOUNDARY_PH = 5.0


def cysteine_thermolysis_arrhenius(ph: float) -> Tuple[float, float]:
    """
    (A in 1/min, Ea in kJ/mol) for cysteine -> H2S at ``ph``.

    Interpolates log10 A and Ea linearly in pH over pH 5-9, the window the
    authors fitted their own pH regressions over. BELOW pH 5 the pH-3.0 pair is
    returned as measured, with no blending: Zheng & Ho identify pH 3.0 as a
    different mechanism ("The nonconformity of reactions at pH 3.0 may mean
    that a different reaction mechanism is involved") and their own data invert
    the pH ordering there at 100 C. Above pH 9 the pH-9 pair is clamped.
    """
    points = ZHENG_CYSTEINE_THERMOLYSIS["ph_points"]
    log_a = [
        math.log10(a * 60.0)
        for a in ZHENG_CYSTEINE_THERMOLYSIS["prefactor_per_second"]
    ]
    eas = ZHENG_CYSTEINE_THERMOLYSIS["ea_kj_mol"]
    x = float(ph)
    if x < ZHENG_MECHANISM_BOUNDARY_PH:
        return 10.0 ** log_a[0], eas[0]
    x = min(x, points[-1])
    for i in (1, 2):
        lo, hi = points[i], points[i + 1]
        if lo <= x <= hi:
            w = (x - lo) / (hi - lo)
            return 10.0 ** (log_a[i] + w * (log_a[i + 1] - log_a[i])), (
                eas[i] + w * (eas[i + 1] - eas[i])
            )
    raise AssertionError("unreachable: pH clamped into the measured range")


def cysteine_thermolysis_k(ph: float, temperature_k: float) -> float:
    """k for Cys -> H2S, 1/min, at ``ph`` and ``temperature_k``."""
    a, ea = cysteine_thermolysis_arrhenius(ph)
    return a * math.exp(-ea / (R_KJ * float(temperature_k)))


#: (ii) AQUEOUS ACID-BASE SPECIATION. Two equilibrium constants, both physical.
#: NEITHER IS FITTED, and neither is a rate.
#:
#: PROVENANCE, stated plainly because it is the one place this module reaches
#: outside its own dossiers: these are standard aqueous pKa values at 25 C, not
#: numbers extracted from the sulfur corpus. They carry the flag
#: `constant_not_from_corpus_dossier`. Their TEMPERATURE dependence is NOT
#: modelled (pKa1 of H2S falls by roughly half a unit between 25 and 100 C),
#: which is a declared approximation and is reported.
PKA_H2S_1 = 7.05          # H2S <-> HS- + H+
PKA_CYSTEINE_AMINE = 10.28  # the alpha-ammonium of cysteine
PH_NORMALISATION = 5.0    # every pH factor equals 1 here, so k_ref is "k at pH 5"


def neutral_h2s_fraction(ph: float) -> float:
    """Fraction of the sulfide pool present as NEUTRAL H2S."""
    return 1.0 / (1.0 + 10.0 ** (float(ph) - PKA_H2S_1))


def free_amine_fraction(ph: float) -> float:
    """Fraction of cysteine carrying a FREE (unprotonated) alpha-amino group."""
    return 1.0 / (1.0 + 10.0 ** (PKA_CYSTEINE_AMINE - float(ph)))


def ph_factor(kind: str, ph: float) -> float:
    """
    The multiplicative pH factor of a step, normalised to 1 at pH 5.

    THREE KINDS, and the assignment of a step to a kind is a STRUCTURAL claim
    about its mechanism, not a fitted choice:

      * ``"acid"``  -- 1,2-enolisation / dehydration to the 3-deoxyosone, and
        everything downstream of it (furfural, and hence FFT). Acid-catalysed,
        so the factor is [H+] relative to pH 5. THIS IS WHAT MAKES FFT FALL
        MONOTONICALLY WITH pH, and it is measured three independent ways:
        Hofmann 1998 T2 ribose FFT 229/121/12 ppb at pH 3/5/7; Zhou 2023 T1 FFT
        814/758/325 ug/L at pH 6/7/8; Cerny 2007 T2 FFT 431/368/364/185/0 peak
        areas at pH 4.0-7.0, falling to EXACT ZERO.

      * ``"base"``  -- 2,3-enolisation to the 1-deoxyosone, and the Nedvidek
        1,4-dideoxyosone step. Amine-catalysed, so the factor is the free-amine
        fraction relative to pH 5.

      * ``"neutral_h2s"`` -- any step whose nucleophile is the sulfide. The
        reactive species is taken to be NEUTRAL H2S, so the factor falls
        through pKa1 = 7.05.

    THE POINT: [1-deoxypentosone] rises with pH while [neutral H2S] falls
    through 7.05, so their PRODUCT -- which is what the MFT-forming flux is --
    has an interior maximum near the sulfide pKa. The pH-7 MFT maximum that
    Zhou 2023 measures is therefore a consequence of two measured shapes
    multiplying, not of a fitted pH term. There is no fitted pH parameter
    anywhere in this module.
    """
    if kind == "":
        return 1.0
    if kind == "acid":
        return 10.0 ** (PH_NORMALISATION - float(ph))
    if kind == "base":
        return free_amine_fraction(ph) / free_amine_fraction(PH_NORMALISATION)
    if kind == "neutral_h2s":
        return neutral_h2s_fraction(ph) / neutral_h2s_fraction(PH_NORMALISATION)
    raise ValueError(f"unknown pH factor kind {kind!r}")


PH_TERM_PROVENANCE: Dict[str, Any] = {
    "fitted_ph_parameters": 0,
    "what_enters": [
        "Zheng & Ho 1994's MEASURED four-pH cysteine thermolysis Arrhenius set "
        "(the source-consistent matched (Ea, A) pairs of "
        "zheng1994_extraction.md sec. 5b, interpolated over pH 5-9 only)",
        "the aqueous pKa1 of H2S (7.05), giving the reactive neutral-sulfide "
        "fraction",
        "the aqueous alpha-ammonium pKa of cysteine (10.28), giving the free-amine "
        "fraction that catalyses 2,3-enolisation",
        "the assignment of each enolisation branch to acid or base catalysis, "
        "which is a mechanism claim and carries no number",
    ],
    "what_does_not_enter": [
        "any fitted pH polynomial",
        "any per-pH rate constant",
        "any pH-shape parameter",
        "the temperature dependence of either pKa (declared approximation)",
    ],
    "where_the_product_law_puts_its_maximum": (
        "DERIVED, not fitted, and it is a sharper result than the build brief "
        "assumed. The product of the two measured shapes is PEAKED for every "
        "combination of exponents, and the peak's position is fixed by the "
        "exponents alone: base^1 x sulfide^2 peaks at pH 7.05 (the SULFIDE "
        "pKa, exactly); base^1 x sulfide^1 at 8.66; base^2 x sulfide^2 at "
        "8.67; base^2 x sulfide^1 at pH 10.28 (the CYSTEINE AMINE pKa, "
        "exactly). "
        "THE NETWORK'S OWN MFT FLUX IS base^2 x sulfide^1 -- 2,3-enolisation "
        "and the Nedvidek step are both base-catalysed and only the "
        "thiol-forming step takes the sulfide -- so its ALGEBRAIC maximum sits "
        "at 10.28, three pH units above the maximum Zhou measures. The product "
        "law alone therefore does NOT deliver the pH-7 maximum for this "
        "topology. What can is the COMPETITION FOR CYSTEINE: thermolysis to "
        "H2S accelerates with pH (Zheng & Ho, measured) and burns the same "
        "pool the carbon skeleton needs as a REAGENT (Nedvidek), so the two "
        "factors are coupled through a shared depleting reactant rather than "
        "merely multiplied. Whether that is enough is an out-of-sample "
        "question and is scored as one."
    ),
    "known_to_be_unsatisfiable": (
        "INVENTORY sec. B2.5 RECORDS A SIGN CONFLICT THIS MODULE CANNOT "
        "RESOLVE. Hofmann 1998 (buffered, free-sugar-fed) has MFT FALLING 19.8 "
        "-> 2.5 ug from pH 5 to 7; Zhou 2023 (unbuffered, Amadori-fed) has MFT "
        "RISING 697 -> 1589 ug/L from pH 6 to 7, with the absolute levels 64x "
        "apart at pH 7. The inventory's own recommendation is to ingest them as "
        "SEPARATE constraint records and NOT to merge them into one pH response "
        "curve. A single pH-free-parameter model necessarily gets at most one of "
        "them right. This module gets Zhou's shape and misses Hofmann's, and the "
        "hold-out report says so rather than splitting the difference."
    ),
    "ph_trajectory": (
        "Zhou's pH labels are INITIAL pH of an UNBUFFERED system (inventory "
        "sec. C.11): the pH-7 run ENDS at 3.42 and the pH-6 run at 3.22, i.e. "
        "within 0.2 units of each other. ``sulfur.integrate_sulfur`` therefore "
        "accepts a pH TRAJECTORY (initial -> measured final, linear in time) "
        "and re-evaluates every pH-dependent constant inside the right-hand "
        "side. BOTH endpoints are measured (sec. A.3.3, 'final pH' rows). This "
        "is what FIT_HOLDOUT_DECLARATION.md sec.5 decision 1 asks for before "
        "the pH-6 column can be anything but diagnostic."
    ),
}


# ===========================================================================
# 2. THE MEASURED THIOL SINK -- Module 2's calibration backbone
# ===========================================================================
# THE CROSS-VALIDATION (inventory sec. A.4, "THE CROSS-VALIDATION"): two labs,
# two methods, two matrices, agreement within 25%.
#     Hofmann 2002 Fig. 6   9.4e-4 /s   30 C, SIDA, 12.5 g/L coffee melanoidin
#     Hofmann 2002 Table 2  9.8e-4 /s   30 C, SIDA, CROSSPY model
#     Charles-Bernard 2005  >7.7e-4 /s  25 C, headspace ratio (a LOWER BOUND)
# => the inventory's own recommendation: adopt k(FFT sink) ~ 9e-4 /s at 30 C,
#    ~10 g/L coffee solids, pH 5-6.
#
# THE BIMOLECULAR RECAST is what makes this a mechanism rather than a fudge.
# Charles-Bernard titrated the matrix ELECTROPHILE SITE DENSITY at 8-10 mmol/g
# dry solids -- the only published value. At Hofmann's 12.5 g/L melanoidin that
# is 100-125 mmol/L of sites, so
#     k2 = 9.4e-4 /s / 0.1125 mol/L = 8.35e-3 L/(mol*s)
# which reproduces the inventory's independently derived k2(FFT) >~ 8.6e-3
# L/(mol*s) to 3%. In this module's units: 5.01e-4 L/(mmol*min).
#
# AND IT IS DEPLETABLE, which is the point. Hofmann's 80 C real-coffee brew is
# SLOWER than his own 30 C model systems, and the dossier's stated reason is
# that the brew's electrophile pool was partly consumed during extraction. A
# pseudo-first-order sink cannot express that; a bimolecular sink on a finite
# pool can. That row is a declared HOLD-OUT and this is the structure that will
# be scored against it.
MELANOIDIN_SITE_DENSITY_MMOL_PER_G: Tuple[float, float] = (8.0, 10.0)
MELANOIDIN_SITE_DENSITY_USED_MMOL_PER_G = 9.0  # the midpoint, stated not hidden

#: THE CONJUGATION IS REVERSIBLE, AND THAT CHANGES THE LONG-TIME PREDICTION.
#: Stack et al. 2018 (Chem. Res. Toxicol. 31:81) measure thiol-quinone
#: conjugation as an EQUILIBRIUM with K ~ 1e2-1e3 M^-1, i.e. a SEQUESTRATION
#: equilibrium rather than an irreversible sink -- which is the mechanism behind
#: Sun 2019's observation that ADDING cysteine RELEASES bound 2-furfurylthiol.
#: An irreversible-only sink over-predicts permanent aroma loss on long
#: timescales, and this module would have made that error.
#:
#: THE REVERSE RATE IS DERIVED, NOT FITTED: k_reverse = k_forward / K, so the
#: release step introduces NO free parameter. K is carried at the geometric
#: midpoint of the measured decade, with the decade itself recorded.
THIOL_QUINONE_EQUILIBRIUM_M_INV: Tuple[float, float] = (1.0e2, 1.0e3)
THIOL_QUINONE_EQUILIBRIUM_USED_L_PER_MMOL = 10.0 ** 2.5 / 1000.0  # 0.3162 L/mmol
THIOL_QUINONE_EQUILIBRIUM_ANCHOR = (
    "Stack et al. 2018, Chem. Res. Toxicol. 31:81 (thiol-quinone conjugation "
    "equilibrium, K ~ 1e2-1e3 M^-1); mechanism corroborated by Sun 2019, in "
    "which added cysteine RELEASES bound FFT. Geometric midpoint 10^2.5 M^-1 = "
    "0.3162 L/mmol is used and the decade is carried."
)

#: CONTEXT for k_thioether's magnitude, carried and NOT operative. Quinone +
#: thiol second-order constants from the biomedical literature are
#: k2 ~ 5e5-7e5 M^-1 s^-1 for cysteine-type thiols, against ~4.5-17 M^-1 s^-1
#: for amines -- a FIVE-ORDER selectivity gap. At those constants a NANOMOLAR
#: electrophile pool already reproduces the coffee-matrix loss rate, which is
#: the independent reason to believe the ~9e-4 1/s constant is pseudo-first-
#: order and matrix-borne rather than intrinsic to the thiol. Not operative
#: because it is not measured in a Maillard matrix and the corpus's own
#: bimolecular recast (8.6e-3 L/(mol*s)) is what the fit panel constrains.
QUINONE_THIOL_K2_CONTEXT_M_INV_S_INV: Tuple[float, float] = (5.0e5, 7.0e5)
QUINONE_AMINE_K2_CONTEXT_M_INV_S_INV: Tuple[float, float] = (4.5, 17.0)

_THIOETHER_CONDITIONS = (
    "pH 6.0, 0.1 M phosphate, 12.5 g/L coffee melanoidin (MW > 3000), "
    "FFT 438 uM, 30 C; and reconstituted coffee at 1% total solids, 0.01 M "
    "acetate pH 5.2, aerobic, 25 C"
)


def _sulfur_parameter(
    key: str,
    transformation: str,
    order: int,
    *,
    k_ref: Optional[float],
    ea: Optional[float],
    evidence_class: str,
    source_anchor: str,
    dossier_anchor: str,
    conditions: str,
    ph: Optional[float],
    t_ref_k: float,
    t_range: Tuple[float, float],
    rate_transfer: str,
    channel: str = "",
    ph_factor_kind: str = "",
    flags: Tuple[str, ...] = (),
    note: str = "",
) -> SulfurParameter:
    return SulfurParameter(
        key=key,
        transformation=transformation,
        k_ref=k_ref,
        ea_kj_mol=ea,
        unit=UNIT_BY_ORDER[order],
        order=order,
        evidence_class=evidence_class,
        source_anchor=source_anchor,
        dossier_anchor=dossier_anchor,
        conditions=conditions,
        ph_of_measurement=ph,
        t_ref_k=t_ref_k,
        temperature_range_c=t_range,
        rate_transfer=rate_transfer,
        channel=channel,
        ph_factor=ph_factor_kind,
        flags=flags,
        note=note,
    )


MEASURED_SULFUR: Mapping[str, SulfurParameter] = {
    "k_thioether": _sulfur_parameter(
        "k_thioether",
        "R-SH + matrix electrophile site -> matrix-bound thioether",
        2,
        k_ref=5.01e-4,
        ea=None,
        evidence_class="measured_rate",
        source_anchor=(
            "Hofmann & Schieberle 2002 Fig. 6 (9.4e-4 /s, 30 C, SIDA) and "
            "Table 2 (9.8e-4 /s); Charles-Bernard, Roberts & Kraehenbuehl 2005 "
            "Table 2 (>7.7e-4 /s, 25 C) as the corroborating lower bound; recast "
            "bimolecular on Charles-Bernard's titrated site density of 8-10 "
            "mmol per g dry solids (p. 4428)"
        ),
        dossier_anchor=(
            "k3_final_parameter_inventory.md sec. A.4 rows at 25/30 C and the "
            "'THE CROSS-VALIDATION' block; k1_kinetic_parameters.md secs. 1a, 1c"
        ),
        conditions=_THIOETHER_CONDITIONS,
        ph=5.6,
        t_ref_k=303.15,
        t_range=(25.0, 30.0),
        rate_transfer="not_licensed",
        channel="covalent_addition_to_matrix_electrophiles",
        flags=(
            "no_Ea_available",
            "temperature_held_fixed",
            "bimolecular_recast_from_pseudo_first_order",
            "site_density_midpoint_9_mmol_per_g",
            "charles_bernard_table2_unit_erratum_values_are_per_second",
        ),
        note=(
            "THE UNIT ERRATUM travels with the number: Charles-Bernard's Table 2 "
            "header prints '[mol^-1 s^-1]' and the values are s^-1, proven "
            "against five of the paper's own half-lives (inventory sec. F row 6). "
            "THE DISULFIDE BRANCH IS DEAD AT 30 C: Hofmann 2002 measures it at "
            "<1.5% of the thiol flux, which fixes the stoichiometry as 1:1 "
            "thioether addition, FIRST order in thiol. THE 25 C LADDER spans "
            "755x across ten thiols (FFT >7.7e-4 down to 2M2P 1.02e-6 /s) and "
            "FFT is the fastest member; MFT is not in the ladder, so MFT is "
            "assigned the same channel constant and that assignment is a "
            "declared assumption, not a measurement."
        ),
    ),
    "k_cys_h2s": _sulfur_parameter(
        "k_cys_h2s",
        "L-cysteine -> H2S (+ unmeasured C3/N1 residue -> FRAG_C, FRAG_N)",
        1,
        k_ref=cysteine_thermolysis_k(5.0, T_REF_S_K),
        ea=133.0,
        evidence_class="measured_activation_energy",
        source_anchor=ZHENG_CYSTEINE_THERMOLYSIS["source_anchor"],
        dossier_anchor=ZHENG_CYSTEINE_THERMOLYSIS["dossier_anchor"],
        conditions=(
            "aqueous 0.1 M L-cysteine, phosphate-buffered to a NOMINAL "
            "room-temperature pH of 3/5/7/9, N2-flushed and sealed, 80-110 C, "
            "first-order H2S release measured by sulfide/silver ISE (TOTAL "
            "sulfide). No sugar, no dicarbonyl, no lipid: neat thermolysis."
        ),
        ph=5.0,
        t_ref_k=T_REF_S_K,
        t_range=(80.0, 110.0),
        rate_transfer="licensed_over_ph_5_to_9_ph3_is_a_separate_mechanism",
        channel="cysteine_thermolysis",
        flags=(
            "matched_Ea_A_pair_from_zheng1994_extraction_sec_5b",
            "supersedes_repo_130.4_with_1e14_which_is_15x_too_fast",
            "prefactor_is_a_refit_never_printed_by_the_authors",
            "ph_resolved_measured_set",
            "ph3_is_a_separate_mechanism_not_interpolated_through",
            "initial_rate_only_max_12_percent_conversion_in_the_source",
            "upper_bound_on_the_free_H2S_source_term_no_maillard_partner",
            "nominal_room_temperature_ph_no_hot_ph_correction",
            "no_error_bars_published",
            "extrapolated_above_110C_in_every_fit_system",
        ),
        note=(
            "THE ONLY MEASURED pH DEPENDENCE IN THE MODULE, and the only sulfur "
            "step with a real activation energy. The k_ref printed here is the "
            "pH-5 value at 145 C; at run time the constant is re-evaluated from "
            "the pH-resolved set at the system's own pH, so the pH-5 number is a "
            "reporting convenience and not the operative constant. "
            "THE PAIR IS THE SOURCE-CONSISTENT ONE (zheng1994_extraction.md "
            "sec. 5b): the repository's shipped 130.4 kJ/mol with A = 1.0e14 "
            "1/s runs ~15x FASTER than Zheng's own Table I at pH 7 / 100 C, "
            "because the Ea was averaged across pH while A was re-fitted free. "
            "EXTRAPOLATION, stated twice because it matters twice: every system "
            "in the FIT panel runs at 115-145 C and Zheng measured 80-110 C, so "
            "this barrier is extrapolated 5-35 C beyond its window on EVERY "
            "row; and Zheng's constants are INITIAL-RATE constants validated "
            "only to ~12% conversion, while this module drives the same "
            "first-order law to far higher conversion at 145 C. "
            "DIRECTION OF THE RESULTING BIAS IS KNOWN: with no Maillard "
            "co-reactant present, Zheng's system has nothing consuming the H2S, "
            "so these constants are an UPPER BOUND on the free-sulfide source "
            "term. A model that over-predicts sulfide will over-predict every "
            "thiol, and that is the first place to look if it does."
        ),
    ),
}


# ===========================================================================
# 3. THE NAMED CONSUMPTION CHANNELS -- Module 2's architecture
# ===========================================================================
# Four papers, four temperatures, four dominant sinks, and each excludes the
# others' (inventory sec. B.7). They are four objects here, not four values of
# one constant.
THIOL_CHANNELS: Tuple[Dict[str, Any], ...] = (
    {
        "channel": "covalent_addition_to_matrix_electrophiles",
        "reactions": ("ch_thioether_mft", "ch_thioether_fft"),
        "dominant_at_c": (25.0, 30.0),
        "order": "second (thiol x electrophile SITE)",
        "parameter": "k_thioether (MEASURED, no Ea)",
        "what_excludes_the_neighbour": (
            "the disulfide branch is <1.5% of the flux at 30 C -- Hofmann 2002, "
            "measured (<6 ug disulfide from 400 ug FFT consumed)"
        ),
        "role": "FIT (Charles-Bernard 25 C ladder + T3 matrix; Hofmann 2002 30 C)",
        "reversible": (
            "YES, and it is implemented as such. Stack 2018 measures the "
            "thiol-quinone conjugation as an equilibrium (K ~ 1e2-1e3 M^-1), "
            "and Sun 2019 shows added cysteine RELEASES bound FFT. The reverse "
            "step ch_thioether_release carries k_forward / K and adds no free "
            "parameter. An irreversible-only sink over-predicts permanent "
            "aroma loss on long timescales."
        ),
        "electrophile_pool_is_a_lump": (
            "MELE lumps at least two measured, chemically distinct channels: "
            "melanoidin-bound pyrazinium (CROSSPY, Hofmann 2002) and the "
            "hydroxyhydroquinone-derived o-quinone, a discrete small molecule "
            "that Muller & Hofmann 2007 identify as the DOMINANT thiol-trapping "
            "electrophile in real brew. No result here separates them."
        ),
    },
    {
        "channel": "acid_catalysed_C5_oligomerisation",
        "reactions": ("ch_oligomer_mft", "ch_oligomer_fft"),
        "dominant_at_c": (50.0, 50.0),
        "order": "ZERO in thiol (59% of initial per DAY for MFT, 28% for FFT)",
        "parameter": "NONE -- see NO_MEASURED_RATE",
        "what_excludes_the_neighbour": (
            "air ~ argon (so it is NOT oxidative), and the MFT mass balance "
            "FAILS to close as thiol + disulfide while MB, FFT and DMFT all "
            "close. Mechanism evidence: 85% H-D exchange at C-5 against 10% at "
            "C-4 -- van Seeventer 2001"
        ),
        "role": (
            "HOLD-OUT (van Seeventer Table 1). The channel is therefore declared "
            "STRUCTURALLY and left at zero rate: there is no FIT row anywhere in "
            "the corpus that constrains it, and giving it a rate would mean "
            "reading the hold-out."
        ),
    },
    {
        "channel": "oxidative_dimerisation",
        "reactions": ("ch_dimer_mft", "ch_dimer_fft"),
        "dominant_at_c": (115.0, 120.0),
        "order": "second in thiol, first in OXIDANT EQUIVALENTS",
        "parameter": "k_dimer_mft, k_dimer_fft (FITTED on Zhang 2024 Fig. 1)",
        "what_excludes_the_neighbour": (
            "the dimer carries up to 49% of the MFT pool at 115 C -- the very "
            "channel the 30 C system rules out at <1.5%. The two cannot be the "
            "same reaction and cannot share an Arrhenius line."
        ),
        "role": "FIT (Zhang 2024 Fig. 1 redox series) / HOLD-OUT (Figs. 2d-f, Zhou 120 C)",
        "why_it_is_gated_on_an_oxidant_and_not_always_on": (
            "METAL-FREE thiol autoxidation in water is NEGLIGIBLE -- "
            "<= 2e-6 1/s (Ngamchuea 2016) -- while Cu-catalysed oxidation is "
            "measured and much faster. So a model system with no oxidant and no "
            "metal SHOULD show essentially no dimerisation, and this module "
            "reproduces that by making the channel first order in an explicit "
            "oxidant pool rather than always-on. The absence of an "
            "uncatalysed-oxidation channel is therefore MATRIX-DEPENDENT and "
            "defensible, not an omission -- but it does mean the module cannot "
            "predict dimerisation in a metal-containing matrix it is not told "
            "about."
        ),
    },
    {
        "channel": "radical_coupling_to_methanethiol",
        "reactions": ("ch_mmft",),
        "dominant_at_c": (115.0, 115.0),
        "order": "second (thiol x methanethiol)",
        "parameter": "k_mmft (FITTED on Zhang 2024 Fig. 1's MMFT column)",
        "what_excludes_the_neighbour": (
            "it requires a methanethiol partner, which the 30 C and 50 C systems "
            "do not have; and its slope BREAKS at 90 min, which no one-stage "
            "rate law can express (Zhang sec. 3.1)"
        ),
        "role": "FIT for Fig. 1; the Fig. 2 dose-response fractions are HOLD-OUT",
    },
)

#: A FULLY PARAMETERISED SET THE MODULE DOES NOT YET HAVE A SPECIES FOR.
#: Carried so the next wave inherits it rather than re-deriving it, and so that
#: the reason for deferring it is on the record.
CYSTINE_SULFIDE_EXCHANGE_DEFERRED: Dict[str, Any] = {
    "reaction": "cystine + sulfide exchange (forward and reverse both measured)",
    "k1_M_inv_min_inv_at_25C": 3.7,
    "ea_kj_mol": 66.1,
    "prefactor_M_inv_min_inv": 4.7e11,
    "source_anchor": "Liu & Chang 1987, Can. J. Chem. 65:770",
    "ph_of_measurement": "ALKALINE (the paper's own condition; not this module's pH 4.5-7)",
    "rate_transfer": "not_licensed",
    "operative": False,
    "why_deferred": (
        "THIS MODULE HAS NO CYSTINE SPECIES. Cystine enters only as the OXIDANT "
        "pool of Zhang 2024's Fig. 1 redox arm, where what is identifiable is "
        "the product k_dimer * [OX] and not the cystine chemistry itself. "
        "Adding a real cystine species with a sulfide-exchange lane would (a) "
        "change the Zhang Fig. 1 FIT row, which would require re-running the "
        "fit, and (b) import an ALKALINE rate into a pH 4.5-7 module, which the "
        "B1 convention does not license. Deferred to a later wave WITH its "
        "numbers intact, rather than half-adopted."
    ),
}

#: Steps that exist in the network, are mechanistically required, and have NO
#: rate obtainable from any DECLARED FIT ROW. They are carried at zero with a
#: flag rather than invented or borrowed. Each names what it would take.
NO_MEASURED_RATE: Mapping[str, str] = {
    "k_oligomer": (
        "van Seeventer 2001's 50 C zero-order MFT/FFT loss is the ONLY "
        "measurement of this channel and it is a declared HOLD-OUT "
        "(FIT_HOLDOUT_DECLARATION.md D.4). Fitting it here would read the "
        "hold-out. Held at 0.0. CONSEQUENCE, pre-registered: the model predicts "
        "NO oligomerisation loss at 50 C and will therefore FAIL the van "
        "Seeventer row. That failure is informative and is reported as such -- "
        "it localises a missing channel, not a wrong barrier."
    ),
}

#: Derivations that are FORBIDDEN, with the reason each is forbidden. Every one
#: is recorded in a dossier as a prohibited derivation; none is reachable by any
#: code path in this module.
PROHIBITED_DERIVATIONS: Tuple[Dict[str, str], ...] = (
    {
        "derivation": "Ea(thiol consumption) from the 30 C and 115 C rates",
        "why_forbidden": (
            "'Pairing Zhang's 115 C rate with K1's 30 C rate to extract an Ea "
            "would repeat exactly the two-point, two-lab, two-matrix, "
            "two-method error that K1 named and refused -- and it would be "
            "worse, because the two rates are not even the same reaction.' "
            "(Zhang2024_extraction.md sec. K3.7, quoted in inventory sec. C.1)"
        ),
        "enforced_how": (
            "every consumption channel carries ea_kj_mol=None and k_at() refuses "
            "to Arrhenius-scale it; there is no fitted Ea in the consumption "
            "block and no code path that pairs two channels"
        ),
    },
    {
        "derivation": "Ea(thiol consumption) from the 30 C model and the 80 C brew",
        "why_forbidden": (
            "'The same treatment applied to FFT gives a NEGATIVE Ea (the 80 C "
            "brew is slower than the 30 C models, because the real brew's "
            "electrophile pool was partly consumed during extraction).' "
            "(k1_kinetic_parameters.md sec. 1d)"
        ),
        "enforced_how": "as above; and the 80 C brew is a declared HOLD-OUT row",
    },
    {
        "derivation": "any fixed thiamine:sugar MFT branch fraction",
        "why_forbidden": (
            "'A 2x change in precursor loading moves the xylose share of MFT "
            "from 15% to 46% -- a 3.1x change in the branch fraction, one lab, "
            "one method, one pH, one temperature. Any model that carries a "
            "fixed thiamine:sugar MFT split is falsified by this one row pair.' "
            "(inventory sec. B10.1, Cerny 2007 Table 5)"
        ),
        "enforced_how": (
            "the registry contains no fraction of any kind; every split is a "
            "ratio of mass-action fluxes computed at run time, and a property "
            "test perturbs a precursor 2x and asserts the share MOVES"
        ),
    },
    {
        "derivation": "van Seeventer's 59 %/day bolted on as THE MFT sink",
        "why_forbidden": (
            "'Do not bolt van Seeventer's 59 %/day onto the network as the MFT "
            "sink.' ... 'the processing sink and the storage sink are DIFFERENT "
            "CHANNELS.' (vanseeventer2001_z3_addendum.md, inventory sec. C.17)"
        ),
        "enforced_how": "k_oligomer is 0.0 and is listed in NO_MEASURED_RATE",
    },
    {
        "derivation": "Zhang 2024's 0.0028 / 0.0031 zero-order constants as numbers",
        "why_forbidden": (
            "'VERDICT ON 0.0028 / 0.0031: DO NOT INGEST AS NUMBERS ... What IS "
            "ingestable is the RATIO 0.0031/0.0028 = 1.107 and the two-stage "
            "structure.' (inventory sec. C.12)"
        ),
        "enforced_how": "neither number appears in this module; only the ratio is scored",
    },
    {
        "derivation": "Bornhorst's norfuraneol k read as a DEGRADATION rate",
        "why_forbidden": (
            "'THIS IS AN APPROACH-TO-PLATEAU FORMATION LAW, NOT A DEGRADATION "
            "LAW ... Anyone reading these k values as norfuraneol degradation "
            "rates would be inverting the paper.' k read as a norfuraneol "
            "DEGRADATION rate: FORBIDDEN. (inventory sec. C.13)"
        ),
        "enforced_how": (
            "no Bornhorst number is operative in this module at all; the "
            "norfuraneol Ea is alkaline (pH 8.4-9.5) and carries "
            "rate_transfer: not_licensed"
        ),
    },
)

#: Alkaline-source numbers carried as UNVALIDATED PRIORS, never operative.
#: FIT_HOLDOUT_DECLARATION.md sec.5 decision 3.
ALKALINE_PRIORS: Tuple[Dict[str, Any], ...] = (
    {
        "quantity": "Ea, norfuraneol (M-2) net accumulation",
        "value_set_kj_mol": [121.1, 122.3, 104.9],
        "ci95_kj_mol": [8.1, 19.5, 8.9],
        "source_anchor": "Bornhorst et al. 2017b, LWT, Table 2 (1_R0.5_L / 1_R1_L / 2_R2_L)",
        "dossier_anchor": "k3_final_parameter_inventory.md sec. A.3.6(i)",
        "conditions": (
            "mashed-potato model gel with gellan gum and CaCl2, 80/90/100 C, "
            "come-up time 1.75 min excluded, pH 8.4-9.5"
        ),
        "ph_of_measurement": 8.95,
        "rate_transfer": "not_licensed",
        "operative": False,
        "why_not_operative": (
            "FOUR MANDATORY QUALIFICATIONS travel with it: (1) it is an "
            "APPARENT, LUMPED approach-to-plateau rate of NET accumulation with "
            "no destruction term at all, not an elementary step; (2) pH 8.4-9.5 "
            "against this module's pH 4.5-7 -- the alkaline-pH wall, sec. C.13; "
            "(3) three temperatures per fit, one of them imported from another "
            "paper; (4) CaCl2 in every row and the gellan matrix data withheld "
            "by the authors. And it is not a constant: across three studies the "
            "norfuraneol Ea spans 64-122 kJ/mol as a function of precursor "
            "loading and matrix (sec. B10.7). NEVER call it 'the norfuraneol "
            "barrier'."
        ),
    },
    {
        "quantity": "Ea, disappearance of furfural in the presence of cysteine",
        "value_kj_mol": 46.50,
        "ci95_kj_mol": 1.0,
        "source_anchor": "Yaghmur et al. 2005, sec. 3.7 p. 232 (the WATER arm)",
        "dossier_anchor": "k3_final_parameter_inventory.md sec. A.3.6(iii), sec. C.14",
        "conditions": "cysteine 6 mmol + furfural 3 mmol in 10 g phosphate 0.5 M, pH 5.0, 40-70 C",
        "ph_of_measurement": 5.0,
        "rate_transfer": "not_licensed",
        "operative": False,
        "why_not_operative": (
            "A LUMP OVER >=98.8% NON-FFT FLUX. The FFT branch is <1.2% of the "
            "furfural flux, so this Ea is overwhelmingly an Ea for something "
            "else, and it must NEVER be used under the step name 'furfural -> "
            "FFT'. It is also at 40-70 C, far below this module's window, and an "
            "Arrhenius Ea is not an Eyring barrier. Filed as an audit_flag, "
            "value unchanged (sec. F row 29). WHAT IS USED from Yaghmur is the "
            "ONE-SIDED CEILING -- the FFT share of the furfural flux is <=1.2% "
            "-- which is a declared FIT row and does constrain this network."
        ),
    },
)

#: The Yaghmur ceiling, the one Yaghmur quantity that IS a FIT row.
YAGHMUR_FFT_SHARE_CEILING = 0.012
YAGHMUR_CEILING_ANCHOR = (
    "Yaghmur et al. 2005 sec. 3.1 p. 226 (<1% FFT conversion of the furfural "
    "charged, 65 C / 15 h) against 84.9% furfural conversion => FFT is <=1.2% "
    "of the furfural flux. FIT_HOLDOUT_DECLARATION.md D.6 Module 9: 'Yaghmur "
    "FFT-share bound -- FIT (as a ceiling)'."
)


# ===========================================================================
# 4. THE FITTED STEPS -- placeholders, populated by the fit
# ===========================================================================
# Every one of these has NO literature rate constant. They are estimated by
# least squares on the DECLARED FIT ROWS ONLY and written back as
# evidence_class='derived_from_fit_data'. Until then k_ref is None and the
# network REFUSES to integrate.
#
# ACTIVATION ENERGIES: the sulfur FIT panel is almost entirely SINGLE
# TEMPERATURE per system (Hofmann 145 C, Whitfield 140 C, van Seeventer 130 C,
# Zhou 120 C, Zhang 115 C) and the systems differ in feedstock as well as in T,
# so a per-step Ea is confounded with the system. ONE LUMPED Ea is fitted for
# the whole FORMATION lane, flagged `lumped_cross_system_weakly_identified`,
# and its interval is reported honestly. The CONSUMPTION lane gets NO Ea at
# all, by policy 2.

_FITTED_SULFUR: Mapping[str, Tuple[str, int, str, str]] = {
    # key: (transformation, order, why it exists, pH factor)
    "k_pent_dpo": (
        "aldopentose + Cys -> 1-deoxypentosone + Cys (2,3-enolisation)", 2,
        "Amine-catalysed 2,3-enolisation. The amine is regenerated, so this is "
        "catalysis and not consumption. BASE-catalysed: the free-amine fraction "
        "is the catalyst concentration.", "base"),
    "k_pent_tdp": (
        "aldopentose + Cys -> 3-deoxypentosone + Cys (1,2-enolisation)", 2,
        "The competing enolisation. ACID-catalysed, which is what makes "
        "everything downstream of it -- furfural, and hence FFT -- fall "
        "monotonically with pH in all three papers that measure it.", "acid"),
    "k_dpo_c2c3": (
        "1-deoxypentosone -> hydroxyacetaldehyde + methylglyoxal (retro-aldol)", 1,
        "THE C2+C3 SPLIT, as one step producing BOTH fragments, so the C2 and "
        "the C3 arrive in the 1:1 stoichiometry the route needs rather than "
        "from two independent constants. Hofmann 1998 T10 measures the C2+C3 "
        "route as the STRONGEST MFT source at 145 C (0.24 mol%) and Cerny 2004 "
        "measures it as 'not relevant' at 95 C -- so the lane is "
        "temperature-scoped and must COMPETE, not be assigned a share.", ""),
    "k_arp_dpo": (
        "pentose Amadori compound -> 1-deoxypentosone + amine (2,3-enolisation)", 1,
        "Zhou 2023 charges the ARP directly. It partitions between the two "
        "enolisations exactly as the free sugar does, which is why Zhou's own "
        "Table 1 has FFT falling monotonically with pH while MFT peaks: the two "
        "products come off OPPOSITE branches of one fed precursor.", "base"),
    "k_arp_tdp": (
        "pentose Amadori compound -> 3-deoxypentosone + amine (1,2-enolisation)", 1,
        "", "acid"),
    "k_dpo_nf": ("1-deoxypentosone -> norfuraneol", 1,
                 "Ribose makes 54 530 ug norfuraneol against 19.8 ug MFT in the "
                 "same pot (a factor of ~2750) while glucose's ratio is ~10. "
                 "That is the authors' own reason for doubting NF as THE key "
                 "MFT intermediate, and it is why NF is a NODE with its own "
                 "competing fates rather than a funnel.", ""),
    "k_dpo_ptr": ("1-deoxypentosone -> 2,3,4-pentanetrione", 1,
                  "The OTHER half of the Nedvidek 1992 Scheme 2 partition, "
                  "verified by two negative controls. It exists so that "
                  "funnelling the whole 1-deoxyosone pool into norfuraneol is "
                  "impossible by construction.", ""),
    "k_dpo_ddp": (
        "1-deoxypentosone + Cys -> 1,4-dideoxypentosone + RCHO + CO2 + NH3", 2,
        "Nedvidek 1992 Scheme 3, balance verified C7H13NO6 both sides for "
        "glycine. THE AMINO ACID IS A REAGENT, NOT A SPECTATOR -- this step "
        "CONSUMES cysteine, which is what puts the sulfur donor and the carbon "
        "skeleton in direct competition for the same pool.", "base"),
    "k_tdp_fur": ("3-deoxypentosone -> 2-furaldehyde", 1, "", ""),
    "k_ddp_mft": ("1,4-dideoxypentosone + H2S -> MFT", 2,
                  "THE INTACT-SKELETON ROUTE. Cerny 2003 T2: 49% unlabelled / "
                  "46% 13C5 with no fragment pattern, 'pathways via ribose "
                  "fragmentation were not relevant' => ~93% of MFT carries the "
                  "intact pentose skeleton at 95 C.", "neutral_h2s"),
    "k_fur_fft": ("2-furaldehyde + H2S -> FFT", 2,
                  "Furfural is 60x ribose and 7x the 3-deoxyosone as an FFT "
                  "source (Hofmann 1998 T3).", "neutral_h2s"),
    "k_nf_mft": ("norfuraneol + H2S -> MFT", 2, "", "neutral_h2s"),
    "k_nf_mp3p": ("norfuraneol + H2S -> 2-mercapto-3-pentanone", 2,
                  "The NF route's DOMINANT fate: mercaptoketones : MFT = 16.3 : "
                  "1 from fed NF, and MFT is only 2.6% of everything fed NF "
                  "produces (Whitfield 1999). Four falsifiable ratio tests, "
                  "immune to the response-factor caveat.", "neutral_h2s"),
    "k_mgo_mp": ("methylglyoxal + H2S -> 1-mercapto-2-propanone", 2,
                 "Hofmann 1998 T7: 1.8 mol% at 1:1 and 4.0 mol% at 1:2. "
                 "Doubling H2S gives 2.2x product -- SUPER-linear, which a "
                 "first-order-in-H2S mass action law under-predicts. The "
                 "residual is reported rather than absorbed into an invented "
                 "reaction order.", "neutral_h2s"),
    "k_ha_mp_mft": ("hydroxyacetaldehyde + 1-mercapto-2-propanone -> MFT", 2,
                    "THE C2+C3 ROUTE, Hofmann 1998 T10's single most effective "
                    "measured MFT route (268.1 ug, 0.24 mol%).", ""),
    "k_glc_ha": ("glucose -> hydroxyacetaldehyde + methylglyoxal + C1 residue", 1,
                 "The hexose entry into the C2+C3 route. Hexoses reach MFT "
                 "ONLY through fragmentation, pentoses through fragmentation "
                 "AND the intact skeleton -- which is the structural reason "
                 "pentose beats hexose 10.4x on MFT in aqueous systems, with no "
                 "fitted sugar-reactivity factor anywhere.", ""),
    "k_thi_hmp": ("thiamine -> 5-hydroxy-3-mercapto-2-pentanone + residue", 1, "", ""),
    "k_thi_mesh": ("thiamine -> methanethiol + residue", 1,
                   "The MeSH source in Zhang's Fig. 1 system, which contains no "
                   "methionine.", ""),
    "k_hmp_mft": ("5-hydroxy-3-mercapto-2-pentanone -> MFT", 1, "", ""),
    "k_hmp_mp2p": ("5-hydroxy-3-mercapto-2-pentanone -> 3-mercapto-2-pentanone", 1,
                   "The THIAMINE-diagnostic isomer. Cerny 2007 T4 has it "
                   "77-90% thiamine-derived while its isomer 2-mercapto-3-"
                   "pentanone is 94->95% xylose-derived, in the same pot.", ""),
    "k_cys_actz": ("Cys + methylglyoxal -> 2-acetylthiazole + C1 residue", 2, "", ""),
    # ---- consumption -------------------------------------------------------
    "k_dimer_mft": ("2 MFT + oxidant -> bis(2-methyl-3-furyl) disulfide", 3, "", ""),
    "k_dimer_fft": ("2 FFT + oxidant -> bis(2-furfuryl) disulfide", 3, "", ""),
    "k_mmft": ("MFT + methanethiol -> MMFT", 2, "", ""),
    "k_mft_decay": ("MFT -> fragments (unassigned sink)", 1,
                    "Bounded to allow ~zero so the data may reject it.", ""),
    "k_fft_decay": ("FFT -> fragments (unassigned sink)", 1,
                    "Bounded to allow ~zero.", ""),
    "k_dimer_decay": ("disulfide dimer -> fragments", 1, "Bounded to allow ~zero.", ""),
    "k_nf_decay": ("norfuraneol -> fragments", 1, "", ""),
    "k_fur_decay": ("2-furaldehyde -> fragments (the large unidentified sink)", 1,
                    "Yaghmur's whole point: the FFT branch is <=1.2% of the "
                    "furfural flux, so almost all of the furfural that "
                    "disappears goes somewhere this corpus never identifies. "
                    "That somewhere is this step, and the ceiling is what "
                    "constrains the ratio between the two.", ""),
    "k_h2s_loss": ("H2S -> unassigned sulfur (volatilisation / non-thiol sinks)", 1, "", ""),
    "k_osone_decay": ("deoxyosone -> fragments (shared by DPO, TDP, DDP)", 1,
                      "ONE constant shared by three species. They are the same "
                      "compound class and nothing in the corpus resolves them "
                      "separately; three free constants would be three ways to "
                      "fit the same residual.", ""),
    "k_thiol_decay": ("minor thiols -> fragments (shared by MP, HMP, MP3P, MP2P)", 1,
                      "One shared constant, same reasoning.", ""),
    "k_pent_caramel": (
        "aldopentose -> 3-deoxypentosone (AMINE-INDEPENDENT 1,2-enolisation)", 1,
        "REQUIRED BY THE FIT ROWS, not optional. Hofmann 1998 T3/T6 feed ribose "
        "+ H2S with NO AMINE and still measure FFT and MFT; a network in which "
        "every enolisation is amine-catalysed predicts zero for both. Its "
        "existence is independently sized by Ajandouz 2008 sec. 3.4 (25-80% of "
        "A294, 7-55% of A420 amine-independent, share RISING with T) but no "
        "rate for it exists anywhere. Bounded to allow ~zero.", "acid"),
    "k_pent_thermal": (
        "aldopentose -> 1-deoxypentosone (AMINE-INDEPENDENT 2,3-enolisation)", 1,
        "The partner of k_pent_caramel, required by the same two rows. Bounded "
        "to allow ~zero.", ""),
    "k_glc_fur": (
        "glucose -> 2-furaldehyde + C1 residue", 1,
        "REQUIRED BY THE FIT ROWS. FFT's only precursor here is furfural, a C5, "
        "and Hofmann measures FFT from both hexoses (28 and 32 ppb), so a "
        "hexose must shed a carbon to reach it. Without this step the model "
        "predicts exactly zero hexose FFT. Bounded to allow ~zero.", "acid"),
}

FITTED_SULFUR_KEYS: Tuple[str, ...] = tuple(_FITTED_SULFUR)

#: Which fitted steps belong to the FORMATION lane (they share the lumped Ea)
#: and which to the CONSUMPTION lane (no Ea at all, by policy 2).
CONSUMPTION_KEYS: Tuple[str, ...] = (
    "k_dimer_mft", "k_dimer_fft", "k_mmft", "k_mft_decay", "k_fft_decay",
    "k_dimer_decay", "k_nf_decay", "k_fur_decay", "k_h2s_loss",
    "k_osone_decay", "k_thiol_decay",
)
FORMATION_KEYS: Tuple[str, ...] = tuple(
    k for k in FITTED_SULFUR_KEYS if k not in CONSUMPTION_KEYS
)

#: Search bounds on log10(k_ref at 145 C). Deliberately wide, and the fit is
#: started from RANDOM points inside them, so any agreement with a literature
#: number is a result and not an initialisation.
#: The UPPER bound is ~3 /min (first order) or ~3 L/(mmol*min) (higher order)
#: at 145 C. On the fit panel's loadings that is a half-life of well under a
#: second, i.e. instantaneous against a 20-minute hold, so a step at the bound
#: is already "as fast as the data can distinguish" and raising the bound
#: changes no observable -- it only makes the ODE stiffer and the integration
#: slower. This is a NUMERICAL bound and it is stated as one; it was tightened
#: from 10-100 after the optimiser was measured spending its whole budget
#: inside integrations at k > 10 that moved nothing.
#:
#: The LOWER bound of 1e-10 is effectively "this step does not happen", so
#: EVERY channel can be rejected by the data. That direction is the one that
#: matters for honesty and it is deliberately left wide open.
FITTED_SULFUR_BOUNDS_LOG10K: Mapping[str, Tuple[float, float]] = {
    key: (-10.0, 0.5) for key in _FITTED_SULFUR
}

#: The single lumped formation Ea, and its bounds.
LUMPED_FORMATION_EA_BOUNDS: Tuple[float, float] = (20.0, 250.0)


def sulfur_placeholders() -> Dict[str, SulfurParameter]:
    """The fitted sulfur steps, unpopulated. Integration refuses them as-is."""
    out: Dict[str, SulfurParameter] = {}
    for key, (transformation, order, why, ph_kind) in _FITTED_SULFUR.items():
        consumption = key in CONSUMPTION_KEYS
        out[key] = _sulfur_parameter(
            key,
            transformation,
            order,
            k_ref=None,
            ea=None,
            evidence_class="derived_from_fit_data",
            source_anchor=(
                "ESTIMATED IN THIS MODULE by least squares on the declared FIT "
                "rows of FIT_HOLDOUT_DECLARATION.md D.3/D.4 only. No literature "
                "rate constant exists for this step."
            ),
            dossier_anchor=(
                "docs/reference/FIT_HOLDOUT_DECLARATION.md D.3 (Module 1) and "
                "D.4 (Module 2); k3_final_parameter_inventory.md secs. A.3, A.4"
            ),
            conditions=(
                "aqueous, pH 4.5-7, 115-145 C -- the union of the FIT panel's "
                "systems. See the fit report for the per-row conditions."
            ),
            ph=5.0,
            t_ref_k=T_REF_S_K,
            t_range=(115.0, 145.0),
            rate_transfer="not_licensed",
            channel="" if not consumption else "fitted_consumption",
            ph_factor_kind=ph_kind,
            flags=(
                ("fitted_here", "no_literature_value")
                + (("no_Ea_available", "temperature_held_fixed")
                   if consumption else ("shares_lumped_formation_Ea",))
            ),
            note=why,
        )
    return out


def with_fitted_sulfur(
    fitted_log10k: Mapping[str, float], lumped_formation_ea: float
) -> Dict[str, SulfurParameter]:
    """
    Populate the fitted steps from ``{key: log10 k_ref at 145 C}``.

    Formation steps receive ``lumped_formation_ea``; consumption steps receive
    None, and are therefore held at 145 C with an extrapolation warning
    wherever they are evaluated. That asymmetry is policy 2, in code.
    """
    out = sulfur_placeholders()
    for key, log10k in fitted_log10k.items():
        if key not in out:
            raise KeyError(f"{key!r} is not a fitted sulfur step")
        ea = None if key in CONSUMPTION_KEYS else float(lumped_formation_ea)
        out[key] = replace(out[key], k_ref=10.0 ** float(log10k), ea_kj_mol=ea)
    # The one channel that is declared but deliberately unpopulated stays at
    # zero rather than at None, because zero is a PREDICTION (no oligomerisation)
    # and None would merely refuse to run.
    return out


def oligomerisation_rate() -> float:
    """The C-5 oligomerisation channel's rate constant: 0.0, and why."""
    return 0.0


# ===========================================================================
# 5. Policy checks
# ===========================================================================


def assert_no_dft_sulfur(
    registry: Mapping[str, SulfurParameter] = MEASURED_SULFUR
) -> None:
    """Owner policy, pinned: no DFT-derived barrier may enter this module."""
    banned = ("dft", "b3lyp", "wb97", "m06", "cbs-qb3", "g4", "computed", "calculated")
    for key, param in registry.items():
        haystack = " ".join(
            [param.source_anchor, param.dossier_anchor, param.note, param.evidence_class]
        ).lower()
        for token in banned:
            if token in haystack:
                raise ValueError(
                    f"{key}: parameter provenance mentions {token!r}. "
                    f"DFT-derived barriers are refused by owner policy."
                )


assert_no_dft_sulfur()


def sulfur_registry_metadata(
    parameters: Mapping[str, SulfurParameter]
) -> Dict[str, Any]:
    """Full runtime metadata block for a sulfur parameter set."""
    return {
        "reference_temperature_K": T_REF_S_K,
        "ph_term_present": True,
        "ph_term_is_fitted": False,
        "ph_term_provenance": dict(PH_TERM_PROVENANCE),
        "branch_fractions_present": False,
        "branch_fraction_policy": (
            "There is no branch-fraction constant in this module. Every split is "
            "a ratio of mass-action fluxes evaluated at run time."
        ),
        "thiol_channels": [dict(c) for c in THIOL_CHANNELS],
        "cysteine_thermolysis_source_set": dict(ZHENG_CYSTEINE_THERMOLYSIS),
        "cysteine_thermolysis_ph_law_carried_not_operative": dict(
            ZHENG_TABLE_II_PH_LAW
        ),
        "no_measured_rate": dict(NO_MEASURED_RATE),
        "prohibited_derivations": [dict(d) for d in PROHIBITED_DERIVATIONS],
        "alkaline_priors_carried_not_operative": [dict(p) for p in ALKALINE_PRIORS],
        "melanoidin_site_density_mmol_per_g": list(MELANOIDIN_SITE_DENSITY_MMOL_PER_G),
        "thiol_conjugation_is_reversible": {
            "equilibrium_constant_M_inv": list(THIOL_QUINONE_EQUILIBRIUM_M_INV),
            "used_L_per_mmol": THIOL_QUINONE_EQUILIBRIUM_USED_L_PER_MMOL,
            "source_anchor": THIOL_QUINONE_EQUILIBRIUM_ANCHOR,
            "reverse_rate_is_derived_not_fitted": True,
        },
        "quinone_thiol_k2_context_M_inv_s_inv": list(
            QUINONE_THIOL_K2_CONTEXT_M_INV_S_INV
        ),
        "quinone_amine_k2_context_M_inv_s_inv": list(
            QUINONE_AMINE_K2_CONTEXT_M_INV_S_INV
        ),
        "deferred_with_numbers_intact": dict(CYSTINE_SULFIDE_EXCHANGE_DEFERRED),
        "dft_free": True,
        "parameters": {k: p.as_metadata() for k, p in sorted(parameters.items())},
    }
