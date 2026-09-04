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
1. **NO DFT.** Same as B1. ``assert_no_dft_sulfur()`` runs at import. No computed
   barrier is read, imported or referenced (the ``data/qm/`` fixtures were deleted
   with the QM/DFT lane on 2026-09-01).

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
        if self.ph_factor not in (
            "", "acid", "base", "neutral_h2s", "hs_anion", "thiolate"
        ):
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


#: (ii) AQUEOUS ACID-BASE SPECIATION, **AT REACTION TEMPERATURE** (B2.1).
#:
#: ===================================================================
#: WHAT B2 GOT WRONG HERE, AND WHY IT WAS THE WAVE'S DOMINANT DEFECT
#: ===================================================================
#: B2's hold-out report named one defect behind most of its 20 failures: "the
#: structural pH term uses acid-catalysis proportional to [H+] and base-
#: catalysis proportional to the free-amine fraction, i.e. ONE DECADE PER pH
#: UNIT ... far too much SLOPE". Two independent things were wrong:
#:
#:   1. **THE pKa VALUES WERE 25 C VALUES USED AT 115-145 C.** B2 recorded this
#:      as "a declared approximation"; it is not a small one. Ionisation of a
#:      neutral acid is endothermic, so every pKa in this module FALLS with
#:      temperature, and the two that matter fall past the operating pH window:
#:      the cysteine alpha-ammonium moves from 10.28 to ~8.1 at 145 C and the
#:      sulfide's pKa1 from 7.05 to ~5.9. A speciation curve evaluated three
#:      units off its own midpoint is in its ASYMPTOTIC region, where the slope
#:      is exactly one decade per pH unit -- which is precisely the slope B2
#:      measured itself having. Correcting the pKa to reaction temperature moves
#:      the base lane's midpoint INTO the operating window and its effective
#:      slope from 3.00 decades to 1.83 decades over pH 5-8. That is a MEASURED
#:      correction (van 't Hoff on published ionisation enthalpies), not a
#:      fitted one.
#:
#:   2. **EVERY pH-SENSITIVE LANE WAS 100% CATALYSED.** If a product has only
#:      one route and that route is fully rate-limited by a catalyst whose
#:      concentration is a decade per pH unit, then the product is a decade per
#:      pH unit -- by construction, with no chemistry left over. Real enolisation
#:      is a SUM of a catalysed and an uncatalysed route, and the sulfide
#:      nucleophile is a SUM over its two protonation states. B2.1 therefore
#:      makes each pH-sensitive lane a PAIR OF PARALLEL ELEMENTARY STEPS (see
#:      sulfur.SULFUR_REACTIONS: r_arp_tdp / r_arp_tdp_thermal, r_ddp_mft /
#:      r_ddp_mft_hs, r_fur_fft / r_fur_fft_hs, r_pent_tdp / r_pent_caramel).
#:      The effective pH response is then the ratio of two ORDINARY RATE
#:      CONSTANTS, fitted on declared FIT rows like every other rate constant in
#:      the module -- it is NOT a pH exponent, a pH polynomial or a pH-shape
#:      parameter, and the count of fitted pH parameters is still ZERO.
#:
#: PROVENANCE OF THE NUMBERS, stated plainly because this is the one place the
#: module reaches outside its own dossiers: the four (pKa at 25 C, ionisation
#: enthalpy) pairs below are standard aqueous physical chemistry, not values
#: extracted from the sulfur corpus. Every one carries the flag
#: `constant_not_from_corpus_dossier` in PH_TERM_PROVENANCE. They are
#: equilibrium constants, not rates, and none of them is fitted.
PKA_25C: Mapping[str, float] = {
    "h2s_1": 7.05,          # H2S <-> HS- + H+
    "cysteine_amine": 10.28,  # the alpha-ammonium of cysteine
    "thiol": 8.33,          # R-SH <-> R-S- + H+, cysteine's own thiol as the
                            # generic reference (see THIOL_PKA_CAVEAT)
}

#: Standard ionisation enthalpies, kJ/mol, all POSITIVE (ionisation of a
#: neutral acid is endothermic), which is what makes every pKa fall with
#: temperature. Same provenance class as the pKa themselves.
DELTA_H_IONISATION_KJ_MOL: Mapping[str, float] = {
    "h2s_1": 22.0,
    "cysteine_amine": 43.0,
    "thiol": 26.0,
}

#: The temperature the pKa values are quoted at.
PKA_REFERENCE_T_K = 298.15

PH_NORMALISATION = 5.0    # every pH factor equals 1 here, so k_ref is "k at pH 5"

THIOL_PKA_CAVEAT = (
    "The thiolate factor uses CYSTEINE's thiol pKa (8.33) as a generic thiol "
    "reference. No pKa for 2-methyl-3-furanthiol or 2-furfurylthiol appears "
    "anywhere in this corpus. The two are not the same compound class -- FFT is "
    "a benzylic-type alkanethiol (pKa nearer 9.4 by analogy with benzyl "
    "mercaptan) and MFT is a heteroaryl thiol (nearer 7 by analogy with "
    "thiophenol) -- so a single reference value is an approximation that must "
    "be reported, and the direction of its error differs between the two "
    "thiols. What the factor is being asked to deliver is a SHAPE that is "
    "sub-decade over pH 3-7 at reaction temperature, and every candidate pKa in "
    "that range delivers one; the level is carried by the fitted rate constant."
)


def pka_at(name: str, temperature_k: float) -> float:
    """
    The pKa of ``name`` at ``temperature_k``, by van 't Hoff.

    pKa(T) = pKa(298.15) + (dH / (R ln10)) * (1/T - 1/298.15)

    Endothermic ionisation (dH > 0) therefore LOWERS the pKa as temperature
    rises, which is the correction B2 omitted. At this module's 418.15 K
    reference the three pKa move 7.05 -> 5.94, 10.28 -> 8.09 and 8.33 -> 6.99.
    """
    if name not in PKA_25C:
        raise KeyError(f"no pKa registered for {name!r}")
    shift = (
        float(DELTA_H_IONISATION_KJ_MOL[name]) / (R_KJ * math.log(10.0))
    ) * (1.0 / float(temperature_k) - 1.0 / PKA_REFERENCE_T_K)
    return float(PKA_25C[name]) + shift


def neutral_h2s_fraction(ph: float, temperature_k: float = T_REF_S_K) -> float:
    """Fraction of the sulfide pool present as NEUTRAL H2S, at temperature."""
    return 1.0 / (1.0 + 10.0 ** (float(ph) - pka_at("h2s_1", temperature_k)))


def hydrosulfide_fraction(ph: float, temperature_k: float = T_REF_S_K) -> float:
    """Fraction present as the HYDROSULFIDE ANION HS-, at temperature."""
    return 1.0 - neutral_h2s_fraction(ph, temperature_k)


def free_amine_fraction(ph: float, temperature_k: float = T_REF_S_K) -> float:
    """Fraction of cysteine carrying a FREE alpha-amino group, at temperature."""
    return 1.0 / (
        1.0 + 10.0 ** (pka_at("cysteine_amine", temperature_k) - float(ph))
    )


def thiolate_fraction(ph: float, temperature_k: float = T_REF_S_K) -> float:
    """Fraction of a thiol present as the THIOLATE ANION, at temperature."""
    return 1.0 / (1.0 + 10.0 ** (pka_at("thiol", temperature_k) - float(ph)))


def ph_factor(kind: str, ph: float, temperature_k: float = T_REF_S_K) -> float:
    """
    The multiplicative pH factor of a step, normalised to 1 at pH 5 AND at the
    SAME temperature -- so that a fitted k_ref always means "k at pH 5", with no
    temperature dependence smuggled in through the normalisation.

    FIVE KINDS. The assignment of a step to a kind is a STRUCTURAL claim about
    its mechanism and carries no number:

      * ``"acid"`` -- 1,2-enolisation / dehydration to the 3-deoxyosone. Acid-
        catalysed, so the factor is [H+] relative to pH 5. IT IS NO LONGER THE
        ONLY ROUTE TO ITS PRODUCT: every acid-catalysed step now runs in
        parallel with an uncatalysed thermal twin, so the OBSERVABLE pH slope is
        sub-decade and is set by the ratio of two fitted rate constants.

      * ``"base"`` -- 2,3-enolisation to the 1-deoxyosone, and the Nedvidek
        1,4-dideoxyosone step. Amine-catalysed: the factor is the free-amine
        fraction AT REACTION TEMPERATURE, whose pKa has moved to ~8.1 at 145 C.

      * ``"neutral_h2s"`` -- sulfide addition proceeding through the NEUTRAL
        molecule.

      * ``"hs_anion"`` -- sulfide addition proceeding through the HYDROSULFIDE
        ANION. B2 had only the neutral branch, which is a mechanism claim that
        runs against the elementary chemistry: HS- is the vastly better
        nucleophile, and the neutral molecule dominates only by abundance. A
        nucleophile with two protonation states adds through both, and carrying
        both branches is what turns a decade-per-pH sulfide term into a
        sub-decade one. The two branches are separate REACTIONS with separate
        fitted rate constants, so their crossover is a rate ratio and not a pH
        parameter.

      * ``"thiolate"`` -- thiol oxidation and disulfide formation, which proceed
        through the THIOLATE anion (textbook, and it is why Kumazawa 2003
        measures 2-furfurylthiol survival collapsing from 99.5% to <0.5% between
        pH 3 and pH 7 at 121 C while its DISULFIDE grows monotonically over the
        same grid). B2 had no pH dependence on either lane, which is why it
        loaded the whole observed FFT-versus-pH slope onto FORMATION.

    WHERE THE MAXIMUM COMES FROM. [1-deoxypentosone] rises with pH while
    [neutral H2S] falls, so their product is peaked; at reaction temperature
    both midpoints sit inside pH 5-9 rather than three units outside it, so the
    peak is broad instead of knife-edged. Nothing here is fitted.
    """
    t_k = float(temperature_k)
    if kind == "":
        return 1.0
    if kind == "acid":
        return 10.0 ** (PH_NORMALISATION - float(ph))
    if kind == "base":
        return (
            free_amine_fraction(ph, t_k)
            / free_amine_fraction(PH_NORMALISATION, t_k)
        )
    if kind == "neutral_h2s":
        return (
            neutral_h2s_fraction(ph, t_k)
            / neutral_h2s_fraction(PH_NORMALISATION, t_k)
        )
    if kind == "hs_anion":
        return (
            hydrosulfide_fraction(ph, t_k)
            / hydrosulfide_fraction(PH_NORMALISATION, t_k)
        )
    if kind == "thiolate":
        return (
            thiolate_fraction(ph, t_k) / thiolate_fraction(PH_NORMALISATION, t_k)
        )
    raise ValueError(f"unknown pH factor kind {kind!r}")


PH_TERM_PROVENANCE: Dict[str, Any] = {
    "fitted_ph_parameters": 0,
    "revision": (
        "B2.1. B2's structural pH slope was ~1 decade per pH unit and its own "
        "hold-out report named that as the wave's dominant defect. TWO changes, "
        "both structural, neither adding a fitted pH parameter: (1) every pKa is "
        "now evaluated AT REACTION TEMPERATURE by van 't Hoff on published "
        "ionisation enthalpies, which moves the cysteine amine from 10.28 to "
        "8.09 and the sulfide from 7.05 to 5.94 at 145 C and takes the base "
        "lane's slope over pH 5-8 from 3.00 decades to 1.83; (2) no pH-sensitive "
        "lane is 100% catalysed any more -- each runs in parallel with an "
        "uncatalysed or oppositely-charged twin REACTION, so the observable "
        "slope is a ratio of two ordinary fitted rate constants rather than a pH "
        "exponent."
    ),
    "what_enters": [
        "Zheng & Ho 1994's MEASURED four-pH cysteine thermolysis Arrhenius set "
        "(the source-consistent matched (Ea, A) pairs of "
        "zheng1994_extraction.md sec. 5b, interpolated over pH 5-9 only)",
        "the aqueous pKa1 of H2S (7.05 at 25 C) with its ionisation enthalpy "
        "(22.0 kJ/mol), giving BOTH the neutral-sulfide and the hydrosulfide-"
        "anion fraction at reaction temperature",
        "the aqueous alpha-ammonium pKa of cysteine (10.28 at 25 C, 43.0 kJ/mol), "
        "giving the free-amine fraction that catalyses 2,3-enolisation",
        "a generic thiol pKa (8.33 at 25 C, 26.0 kJ/mol) giving the thiolate "
        "fraction that carries thiol oxidation and disulfide formation -- see "
        "THIOL_PKA_CAVEAT for why one value stands for two different thiols",
        "the assignment of each step to acid, base, neutral-sulfide, "
        "hydrosulfide or thiolate, which is a mechanism claim and carries no "
        "number",
    ],
    "what_does_not_enter": [
        "any fitted pH polynomial",
        "any per-pH rate constant",
        "any pH-shape parameter or pH exponent",
        "the ACTIVITY-COEFFICIENT correction to the pKa (ionic strength is "
        "0.1-0.5 M in the buffered systems; a declared approximation)",
        "any hot-pH correction to the LABELS the sources print, which are "
        "room-temperature meter readings in every paper in the corpus "
        "(Kumazawa says so explicitly; Zheng's phosphate shifts -0.3 to -0.5 "
        "units). The pKa move with temperature and the pH LABEL does not, "
        "which is a real and unquantified asymmetry",
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
#: ===================================================================
#: B2.1 -- A CONTRADICTION FOUND, AND THE MEASURED PAIR THAT REPLACES IT
#: ===================================================================
#: B2 wrote "Stack 2018 measures ... K ~ 1e2-1e3 M^-1" and used the geometric
#: midpoint 10^2.5 = 316 M^-1. **NOTHING IN stack2018_extraction.md SUPPORTS
#: THAT RANGE.** The dossier's Figure 4C transcription gives the forward AND
#: reverse constants directly at four temperatures, and sec. 3.1 derives
#: K = k1/k-1 from them: 10.57 / 9.29 / 7.59 / 5.64 M^-1 at 4.6 / 9.6 / 14.5 /
#: 19.4 C. The measured equilibrium constant is therefore **56x SMALLER** than
#: the value B2 used, and B2's number is not traceable to the source. Reported
#: as a contradiction, corrected here, and the direction matters: a 56x smaller
#: K means the adduct is far more readily released than B2 assumed.
#:
#: AND THE SIGN OF ITS TEMPERATURE DEPENDENCE IS THE FOOD-RELEVANT RESULT.
#: K FALLS with temperature -- dossier sec. 3.2, van 't Hoff on the four
#: measured K, dH = -6.81 kcal/mol = -28.5 kJ/mol, R^2 0.961, and independently
#: confirmed by Ea(forward) - Ea(reverse) = 2.58 - 9.39 = -6.81 kcal/mol EXACTLY
#: from the two separate Arrhenius refits. The adduct is exothermic, so heating
#: pushes the equilibrium back toward FREE THIOL. That is the mechanism B2 was
#: missing when it lost FFT 18x too fast in Hofmann's 80 C brew.
#:
#: THE CORRECTED ACTIVATION PARAMETERS ARE THE DOSSIER'S REFIT, NOT THE PAPER'S.
#: FIT_HOLDOUT_DECLARATION.md Amendment 4 is explicit: "Stack 2018 NAC arm
#: forward+reverse constants (CORRECTED set -- the paper's printed activation
#: parameters carry a spurious ln(10); dossier's refit is canonical)". The
#: paper's printed dH-doubledagger values are 2.303x too large across all four
#: series and its printed dG-doubledagger regenerates the measured k 126-5870x
#: too fast. None of the twelve printed Figure 4D values is used here.
STACK_NAC_FORWARD_K_M_INV_S_INV_AT_19_4C = 496.0
STACK_NAC_REVERSE_K_S_INV_AT_19_4C = 88.0
STACK_MEASUREMENT_T_K = 292.55  # 19.4 C
STACK_EA_FORWARD_KJ_MOL = 10.8   # 2.58 kcal/mol, R^2 0.9405 [Z, dossier sec. 6.4]
STACK_EA_REVERSE_KJ_MOL = 39.3   # 9.39 kcal/mol, R^2 0.9938 [Z, dossier sec. 6.4]
STACK_DELTA_H_ADDUCT_KJ_MOL = -28.5  # -6.81 kcal/mol, van 't Hoff, R^2 0.961

#: The measured equilibrium constant, in this module's units.
#: 496 / 88 = 5.636 M^-1 = 5.636e-3 L/mmol at 19.4 C.
THIOL_QUINONE_EQUILIBRIUM_M_INV: Tuple[float, float] = (5.64, 10.57)
THIOL_QUINONE_EQUILIBRIUM_USED_L_PER_MMOL = (
    STACK_NAC_FORWARD_K_M_INV_S_INV_AT_19_4C
    / STACK_NAC_REVERSE_K_S_INV_AT_19_4C
) / 1000.0
THIOL_QUINONE_EQUILIBRIUM_ANCHOR = (
    "Stack et al. 2018, Chem. Res. Toxicol. 31:81, Figure 4C (p. 85): the NAC "
    "arm's forward k1 and reverse k-1 are both MEASURED, at four temperatures. "
    "K = k1/k-1 = 496/88 = 5.64 M^-1 at 19.4 C, rising to 10.57 M^-1 at 4.6 C. "
    "Temperature dependence from the dossier's van 't Hoff analysis (sec. 3.2), "
    "dH = -28.5 kJ/mol, R^2 0.961, cross-validated to the last digit by the two "
    "independent Arrhenius refits (sec. 6.2). Mechanism corroborated by Sun "
    "2019, in which added cysteine RELEASES bound FFT. "
    "SUPERSEDES B2's 'K ~ 1e2-1e3 M^-1', which is not traceable to this source "
    "and is 56x too large."
)
STACK_TRANSFER_CAVEATS = (
    "FOUR, all from stack2018_extraction.md and all carried: (1) the measured "
    "system is a MODEL benzoquinone (BPAQ) with N-acetylcysteine at pH 6.2, not "
    "a Maillard melanoidin with a furanthiol; (2) the forward Arrhenius refit "
    "spans 14.8 K with the rate moving only 1.27x, so Ea(forward) is weakly "
    "determined -- the dossier's own estimate is +/-30-50% -- and the reverse is "
    "the better-determined of the two (R^2 0.9938); (3) the dossier states "
    "plainly that these barriers are 'not adequate to extrapolate to process "
    "temperatures', and this module evaluates them at 30-145 C; (4) the "
    "ABSOLUTE forward constant does NOT transfer -- 496 M^-1 s^-1 against the "
    "titrated matrix site density would predict a thiol half-life of "
    "milliseconds, three orders faster than Hofmann measures -- so what is taken "
    "from Stack is the EQUILIBRIUM and its temperature dependence, while the "
    "forward MAGNITUDE stays Hofmann's and Charles-Bernard's matrix recast."
)


def thiol_adduct_equilibrium_l_per_mmol(temperature_k: float) -> float:
    """
    K(T) for the thiol-electrophile adduct, L/mmol, from Stack's MEASURED pair.

    van 't Hoff about the 19.4 C measurement, with the MEASURED dH of -28.5
    kJ/mol. The adduct is exothermic, so K falls as temperature rises: 5.64e-3
    L/mmol at 19.4 C, 3.74e-3 at 30 C, 7.55e-4 at 80 C, 1.5e-4 at 145 C.
    Heating releases bound thiol, and that is a measurement, not a modelling
    choice.
    """
    return THIOL_QUINONE_EQUILIBRIUM_USED_L_PER_MMOL * math.exp(
        -(STACK_DELTA_H_ADDUCT_KJ_MOL / R_KJ)
        * (1.0 / float(temperature_k) - 1.0 / STACK_MEASUREMENT_T_K)
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


# ===========================================================================
# 2b. THE PROTEIN-DISULFIDE THIOL SINK -- B2.1's new channel
# ===========================================================================
# THIS CHANNEL HAS NO MEASURED RATE AND NO FIT ROW, AND IT IS CARRIED ANYWAY.
# FIT_HOLDOUT_DECLARATION.md Amendment 5 makes it a required channel; the
# source (anantharamkrishnan2020b) reports a saturating dose, a 24 h binary
# endpoint and no y-axis, so nothing in it is a rate. What it does report is
# ONE timescale bracket (sec. 5c), and that is what the constant below is.
#
# THE DERIVATION, in full, because the number is weak enough that it must be
# checkable:
#   * thiol charged 12 g/L of propanethiol = 158 mmol/L (dossier sec. 3);
#   * the 1:1 exchange adduct is already present at 6 h  => the FIRST disulfide's
#     pseudo-first-order half-life is <= 360 min at that loading, so
#     k2 >= ln2 / (360 min * 158 mmol/L) = 1.22e-5 L/(mmol*min);
#   * the 2:1 adduct is "nearly none" at 6 h and "clearly seen" at 24 h => the
#     SECOND exchange's half-life is between 360 and 1440 min, so
#     k2 = 3.0e-6 to 1.2e-5 L/(mmol*min).
# The operative value is the geometric midpoint of that bracket. It is an
# ORDER OF MAGNITUDE, and the bracket travels with it everywhere.
#
# WHY IT IS SAFE TO CARRY AN UNMEASURED CONSTANT HERE, AND ONLY HERE: the
# partner is PROT_SS, a titrated protein-disulfide site pool that is ZERO in
# every system in both the fit panel and the hold-out panel. Kang's TTCA
# solution, Hofmann's phosphate buffers, Zhou's deionised water, Whitfield's,
# Cerny's and Zhang's model solutions and Kumazawa's citrate/phosphate all
# contain no protein. The channel therefore carries EXACTLY ZERO FLUX in every
# scored row, by mass action rather than by a switch, and it can move no number
# in either report. It exists so that the module stops being silently wrong in
# a protein matrix, which is the only place it fires.
#
# THE BOUND'S OWN VERDICT ON ITSELF, from the dossier: at 3 wt% protein the
# implied half-life is 37-760 DAYS, so "the covalent channel CANNOT explain
# ambient losses and bites only at process temperature" -- and there is no
# activation energy anywhere in the source, so this module cannot take it to
# process temperature either. That gap is declared, not papered over.
PROTEIN_DISULFIDE_EXCHANGE_BRACKET_L_PER_MMOL_MIN: Tuple[float, float] = (
    3.0e-6, 1.22e-5,
)
PROTEIN_DISULFIDE_EXCHANGE_USED_L_PER_MMOL_MIN = math.sqrt(
    PROTEIN_DISULFIDE_EXCHANGE_BRACKET_L_PER_MMOL_MIN[0]
    * PROTEIN_DISULFIDE_EXCHANGE_BRACKET_L_PER_MMOL_MIN[1]
)
#: 2 disulfides per beta-lactoglobulin, 10% (w/v), MW 18 400 => 1.09 mmol/L.
PROTEIN_DISULFIDE_SITES_BLG_10PCT_MMOL_L = 1.09
PROTEIN_DISULFIDE_DERIVATION = (
    "NOT A MEASURED RATE. The source reports a saturating dose, a 24 h binary "
    "endpoint and no y-axis on any adduct figure; it contains no rate constant "
    "of any kind. This value is the geometric midpoint of the ONLY timescale "
    "bracket in the paper (sec. 5c: the propanethiol 2:1 adduct is 'nearly "
    "none' at 6 h and 'clearly seen' at 24 h at 158 mmol/L thiol), i.e. "
    "3.0e-6 to 1.22e-5 L/(mmol*min), and the bracket travels with it. "
    "IT IS NOT FITTED AND CANNOT BE: no declared FIT row contains protein, so "
    "the objective has no gradient on it. "
    "IT IS ALSO INERT IN BOTH SCORED PANELS: PROT_SS is zero in every fit and "
    "hold-out system, so the channel's flux is identically zero there and it "
    "moves no number in either report. "
    "TWO THINGS IT IS EXPLICITLY NOT LICENSED FOR: (1) any temperature other "
    "than ambient -- there is no Ea in the source, and the dossier's own "
    "derived ambient bound (t1/2 37-760 days at 3 wt%) says the channel 'CANNOT "
    "explain ambient losses and bites only at process temperature', which is "
    "exactly the regime this constant cannot reach; (2) MFT -- the measurement "
    "is on FFT and on n-propanethiol, and extending it to the furan-3-thiol "
    "class is the dossier's own flagged inference, not a measurement. "
    "AND THE PAPER CARRIES A NEGATIVE GATE THAT BOUNDS ANY GENERIC VERSION OF "
    "THIS TERM: 32 of 47 flavour compounds formed NO adduct at all, which "
    "falsifies any blanket protein-binding loss applied to the whole volatile "
    "set."
)


MEASURED_SULFUR: Mapping[str, SulfurParameter] = {
    "k_thioether": _sulfur_parameter(
        "k_thioether",
        "R-SH + matrix electrophile site -> matrix-bound thioether",
        2,
        k_ref=5.01e-4,
        ea=STACK_EA_FORWARD_KJ_MOL,
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
            "'THE CROSS-VALIDATION' block; k1_kinetic_parameters.md secs. 1a, 1c; "
            "hofmann2002_extraction.md (full text read 2026-09-04: Table 2 gives "
            "9.8e-4 /s exactly, Fig. 6 supports 9.4e-4 within the read-off, "
            "disulfide <6 ug of 400 ug; the model reactions ran at pH 6.0)"
        ),
        conditions=_THIOETHER_CONDITIONS,
        ph=5.6,
        t_ref_k=303.15,
        t_range=(25.0, 30.0),
        rate_transfer="not_licensed",
        channel="covalent_addition_to_matrix_electrophiles",
        flags=(
            "forward_Ea_is_stack2018_dossier_refit_weakly_determined",
            "stack_forward_Ea_spans_14.8K_with_1.27x_rate_change",
            "stack_dossier_says_not_licensed_above_its_window",
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
    "k_protein_ss": _sulfur_parameter(
        "k_protein_ss",
        "R-SH + protein disulfide site -> protein-bound thiol (PRB)",
        2,
        k_ref=PROTEIN_DISULFIDE_EXCHANGE_USED_L_PER_MMOL_MIN,
        ea=None,
        evidence_class="bounded_from_a_timescale_bracket",
        source_anchor=(
            "Anantharamkrishnan & Reineccius 2020b (the parent adduct survey), "
            "Fig. 7c: 2-furfurylthiol forms BOTH 1:1 (+114 Da) and 2:1 (+228 Da) "
            "adducts with beta-lactoglobulin, 'in close analogy to the reactions "
            "with n-PrSH'; p. 12 names the mechanism as 'cleaving one or both of "
            "the (soft) electrophilic disulfide linkages in BLG'. The ONLY two "
            "timescale brackets in the paper (sec. 5c): the propanethiol 1:1 "
            "adduct is present at 6 h and the 2:1 adduct is 'nearly none' at 6 h "
            "but 'clearly seen' at 24 h, at 158 mmol/L thiol, ambient."
        ),
        dossier_anchor=(
            "data/lit/extraction_dossiers/anantharamkrishnan2020b_extraction.md "
            "secs. 5b, 5c, 5d, 6a and gate G-3; "
            "docs/reference/FIT_HOLDOUT_DECLARATION.md Amendment 5 ('New "
            "required channel for the sulfur branch: FFT is consumed by protein "
            "disulfides')"
        ),
        conditions=(
            "10% (w/v) beta-lactoglobulin in doubly distilled water, flavour "
            "dosed at 12 g/L (SATURATING), 24 h endpoint, AMBIENT temperature, "
            "pH unbuffered. Disulfide site inventory 2 per BLG = 1.09 mmol/L."
        ),
        ph=None,
        t_ref_k=298.15,
        t_range=(20.0, 25.0),
        rate_transfer="not_licensed",
        channel="thiol_disulfide_exchange_with_matrix_protein",
        flags=(
            "no_measured_rate_only_a_derived_bracket",
            "no_fit_row_constrains_it",
            "rate_bounded_not_fitted",
            "zero_flux_in_every_fit_and_holdout_system_no_protein_present",
            "saturating_dose_no_dose_response",
            "reversibility_untested_in_either_direction",
            "no_Ea_available",
            "temperature_held_fixed",
            "extrapolated_from_ambient_to_100_145C_in_any_protein_system",
        ),
        note=PROTEIN_DISULFIDE_DERIVATION,
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
    # ---- B2.1 adds two ---------------------------------------------------
    {
        "channel": "thiolate_mediated_oxidative_loss",
        "reactions": ("ch_thiolate_loss_mft", "ch_thiolate_loss_fft"),
        "dominant_at_c": (121.0, 121.0),
        "order": "first in thiol, first in THIOLATE FRACTION",
        "parameter": "k_thiolate_loss (FITTED on Kumazawa 2003's 121 C pH grid)",
        "what_excludes_the_neighbour": (
            "it is measured with NO FORMATION PRESENT -- 1 ppm 2-furfurylthiol "
            "dosed into buffer and heated -- so it cannot be confounded with "
            "the formation lane the way every in-situ pH series can. And it is "
            "the NON-disulfide half specifically: Kumazawa's own mass balance "
            "fails to close ('the total amount of the volatile degradation "
            "products ... was less than the amount of 2-furfurylthiol lost'), "
            "which is what this channel is."
        ),
        "role": "FIT (Kumazawa 2003 Fig. 3, six of seven pH levels)",
        "why_it_is_new": (
            "B2 had NO pH dependence on any consumption channel and was "
            "therefore forced to express the whole observed pH response of FFT "
            "through FORMATION, at a decade per pH unit. Its own hold-out "
            "report named that as the wave's dominant defect."
        ),
        "caveats": (
            "See KUMAZAWA_CAVEATS below, carried in full into runtime metadata "
            "under 'kumazawa2003_ph_grid'. The two that bite hardest: the pH "
            "labels are ROOM-TEMPERATURE meter readings on a citrate/phosphate "
            "buffer heated to 121 C, while this module's pKa DO move with "
            "temperature -- an unquantified asymmetry -- and there are no error "
            "bars, no SD and no n anywhere in the paper."
        ),
    },
    {
        "channel": "thiol_disulfide_exchange_with_matrix_protein",
        "reactions": ("ch_protein_ss_fft", "ch_protein_ss_mft"),
        "dominant_at_c": (20.0, 25.0),
        "order": "second (thiol x protein DISULFIDE SITE)",
        "parameter": (
            "k_protein_ss -- BOUNDED from a timescale bracket, never fitted, "
            "and inert in every scored row because no system in either panel "
            "contains protein"
        ),
        "what_excludes_the_neighbour": (
            "it requires a protein disulfide, which no model solution in the "
            "corpus has; and it is STOICHIOMETRIC rather than catalytic, so its "
            "capacity is the protein's own disulfide inventory"
        ),
        "role": (
            "NEITHER. anantharamkrishnan2020b is declared mechanism_reference: "
            "saturating dose, 24 h binary endpoints, no rates, no y-axes. The "
            "channel is required by FIT_HOLDOUT_DECLARATION.md Amendment 5 and "
            "carried rate-bounded."
        ),
        "known_gap": (
            "NO ACTIVATION ENERGY EXISTS, and the dossier's own derived ambient "
            "bound (t1/2 37-760 days at 3 wt%) says the channel 'CANNOT explain "
            "ambient losses and bites only at process temperature' -- which is "
            "precisely the regime this constant cannot reach. Declared."
        ),
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

# ===========================================================================
# 3a-bis. B8 -- THE K6a / K6b MEASURED BARRIERS, AND THE PROVENANCE CORRECTION
# ===========================================================================
# FIT_HOLDOUT_DECLARATION.md Amendment 17, ratified 2026-08-30. Four things
# happen here and each is a MEASUREMENT replacing either a fitted lump or an
# absence, never a taste:
#
#   (1) the 248 kJ/mol thiol sink is REFUTED and replaced by a PRIOR BAND;
#   (2) the disulfide channels gain the corpus's only measured barrier;
#   (3) the fed-Amadori enolisations gain theirs;
#   (4) the whole Kang ladder is re-pointed at its primary source and
#       re-banded from Tier A to Tier B semi-quant.
#
# ---------------------------------------------------------------------------
# (1) THE COVALENT-CAPTURE BAND -- Gigl, Hofmann & Frank 2021
# ---------------------------------------------------------------------------
#: B2.2's fit put `Ea_decay_thiol_sink` at 248.0 kJ/mol against a 250 ceiling;
#: B2.3 got 216.1 and B2.4 212.9-218.1. NONE of those is compatible with the
#: only direct measurement of a thiol sink's temperature structure that exists.
#: Gigl 2021 (qHNMR, coffee melanoidins, 279 / 300 / 333 K) measures the
#: covalent capture of 2-furfurylthiol at k = 1.17e-5 / 2.57e-4 / 7.83e-4 /min.
#:
#: THE FALSIFICATION IS ARITHMETIC, NOT RHETORICAL. 248 kJ/mol predicts
#: k(333)/k(279) = exp(248000 * 5.81226e-4 / 8.314) = 3.4e7. The MEASURED ratio
#: is 67.2. The prediction is wrong by 5.1e5x. No plateau choice and no rate law
#: recovers a 248 kJ/mol barrier from a 67-fold change over 54 K.
GIGL_EA_THIOL_COVALENT_CAPTURE_KJ_MOL = 60.2
#: The DEFENSIBLE RANGE across all six of the dossier's derivation routes. It is
#: a PRIOR BAND, not a numerical convenience: B8 searches inside it, and a fit
#: that presses either edge is reporting a residual conflict rather than a
#: barrier. (Three-point Arrhenius gives 58.5 with R^2 0.886.)
GIGL_EA_COVALENT_CAPTURE_RANGE_KJ_MOL: Tuple[float, float] = (7.0, 102.0)
GIGL_COVALENT_ANCHOR = (
    "Gigl, Hofmann & Frank 2021, J. Agric. Food Chem. 69:15334-15344, doi "
    "10.1021/acs.jafc.1c06163, Figure 4A (gigl2021_extraction.md secs. 5b, 6d, "
    "7b): 2-furfurylthiol + isolated coffee HMW melanoidins, qHNMR with ERETIC "
    "2 / PULCON absolute quantification, 279 / 300 / 333 K. Ea = "
    "R ln(67.2) / (1/279 - 1/333) = 60.2 kJ/mol; 3-point fit 58.5, R^2 0.886."
)
GIGL_COVALENT_CAVEATS = (
    "(1) A 3-5 DECADE EXTRAPOLATION. Measured at 6-60 C; used at 100-180 C. "
    "Transfer risk is HIGH for magnitude, LOW for sign, MODERATE for order of "
    "magnitude, and this module uses it ONLY as a band on a barrier -- never as "
    "a rate. (2) n = 1, no SDs, no LOD anywhere in the paper. (3) The ligand : "
    "site ratio is ~1e5 above an authentic brew. (4) The binding partner is a "
    "PREFORMED, purified melanoidin re-added to solution, not a melanoidin "
    "forming in situ. (5) The dossier's own calibrated verdict is that this "
    "paper is a REFUTATIONAL and STRUCTURAL anchor, not a calibration source: "
    "P(this paper alone suffices to refit a repo parameter) < 0.05, which is "
    "why what it supplies here is a BAND and a prior centre rather than a value."
)
#: THE MODEL-FORM GAP THIS BAND DOES NOT CLOSE, DECLARED RATHER THAN PATCHED.
#: Gigl's FFT is the only one of 20 odorants in "Scenario IV -- covalent AND
#: pi-pi". Two controls dissect it: N-acetylcysteine (thiol, no ring) binds
#: covalently only; methylated FFT (ring, no free -SH) pi-stacks only. The
#: pi-stacked pool is a REVERSIBLE reservoir that FEEDS the covalent sink, it
#: equilibrates in under five minutes, and its enthalpy is NEGATIVE:
GIGL_DELTA_H_PI_STACKING_KJ_MOL = -19.5
GIGL_PI_STACKING_DECLARATION = (
    "NOT IMPLEMENTED, DECLARED. At 60 min Gigl's FREE thiol fraction RISES with "
    "temperature (24.1% at 279 K, 32.3% at 300 K, 55.0% at 333 K): more thiol "
    "survives HOT than COLD. NO single first-order sink with a positive Ea can "
    "produce that at ANY prefactor -- it is a model-form falsification, not a "
    "fit-quality problem, and it is the mechanism by which an optimiser asked "
    "to represent two opposite-sign channels with one knob drives that knob to "
    "a rail. This module has NO reversible negative-enthalpy sequestration "
    "term. B8 declines to add one: adding a species after reading a hold-out "
    "scorecard is a refit in disguise. The enthalpy is on the record "
    "(-19.5 kJ/mol, R^2 0.979, confirmed independently on a pyrazine at -13.4) "
    "so a later wave can implement it as a pre-registered structural change."
)

# ---------------------------------------------------------------------------
# (2) and (3) THE ZHANG 2026 BARRIERS
# ---------------------------------------------------------------------------
#: Zhang et al. 2026, Food Chem. doi 10.1016/j.foodchem.2026.148681, Table 2 --
#: the ONLY paper in wave K6a that publishes rate constants AND activation
#: energies, and ALL NINETEEN of its Ea reproduce from its own k values (mean
#: |deviation| 2.1%, median 1.7%). The MAGNITUDES are refused (no units printed,
#: mixed reaction orders, 8 of 19 on arbitrary-scale nodes); the BARRIERS
#: survive, because a constant scale cancels in ln(k2/k1).
#:
#: k17, Cys-Gly -> di-(Cys-Gly): THE CORPUS'S ONLY MEASURED BARRIER FOR THIOL
#: OXIDATION TO A DISULFIDE. Before B8 the two dimerisation channels carried no
#: activation energy at all and were held at their 145 C value everywhere --
#: which is not conservatism but a strong claim (that a thermal oxidation runs
#: exactly as fast at 100 C as at 145 C) with a known direction of error.
ZHANG_EA_THIOL_TO_DISULFIDE_KJ_MOL = 122.2
ZHANG_DISULFIDE_ANCHOR = (
    "Zhang et al. 2026, Food Chem., Table 2 step k17 (zhang2026_extraction.md "
    "secs. 4b, 6): k = 0.0008 / 0.0015 / 0.0045 at 130 / 140 / 150 C, "
    "glucose-glutathione, pH 6.5, 0.2 M phosphate. Published Ea 122.9; this "
    "wave's independent refit 122.2, R^2 0.971, deviation -0.6%. Legs 87.1 "
    "(130->140) and 159.7 (140->150) -- the spread is recorded, and it is the "
    "reason the value is carried as a barrier and never as a rate."
)
#: k16, Cys-Amadori -> alpha-dicarbonyl: R^2 = 1.000 with legs 86.7 / 84.7, a
#: 2% spread. It is the FLAT-LEG member of a family (k4, k8, k11, k14) whose
#: other four members saturate on the 140->150 C leg and whose published
#: 53-68 kJ/mol are PLATEAU ARTEFACTS the authors noticed, wrote down, and
#: published a single Arrhenius Ea for anyway. k16 does NOT show it, which is
#: what makes it ingestible where they are not.
ZHANG_EA_CYS_AMADORI_TO_ALPHA_DC_KJ_MOL = 85.7
ZHANG_AMADORI_ANCHOR = (
    "Zhang et al. 2026, Food Chem., Table 2 step k16 (zhang2026_extraction.md "
    "secs. 4b, 6.2): k = 0.0263 / 0.0492 / 0.0881 at 130 / 140 / 150 C. "
    "Published Ea 87.9; refit 85.7, R^2 1.000, legs 86.7 / 84.7. The product "
    "node is a UV absorbance, so the Ea survives and the k does not; only the "
    "Ea is taken. Applied to the two AMINE-CATALYSED fed-Amadori enolisations "
    "(k_arp_dpo, k_arp_tdp), which are what Zhang's buffered pH-6.5 system "
    "measures; the UNCATALYSED partners (k_arp_tdp_th, k_arp_dpo_th) are a "
    "different mechanism and keep the lumped formation barrier."
)
#: k15, Cys + glucose -> Cys-Amadori, Ea = 100.9, R^2 1.000, legs 98.9 / 103.0 --
#: the best-conditioned constant of the whole K6a wave. IT HAS NO SHIPPED SITE.
#: This module's Amadori is FED (Zhou 2023 charges purified ARP at 20 mmol/L);
#: there is no `Cys + sugar -> ARP` edge anywhere in `_FITTED_SULFUR`. Recording
#: it as a declared prior with no site is the honest option; forcing it onto
#: `k_ttca_cys` (a retro-condensation, not an Amadori rearrangement) would be
#: assigning a measurement to a step that is not the one measured.
ZHANG_EA_CYS_TO_CYS_AMADORI_KJ_MOL = 100.9
ZHANG_CYS_AMADORI_FORMATION_NO_SITE = (
    "DECLARED, NO SHIPPED SITE. Zhang 2026 k15 (Cys + glucose -> Cys-Amadori, "
    "Ea 100.9, R^2 1.000, legs 98.9/103.0) is the best-conditioned barrier in "
    "wave K6a and this network has no edge for it: the Amadori compound is a "
    "FED precursor here, not a formed one. Ingested as a prior on the record so "
    "a later wave that adds the formation edge inherits the measurement rather "
    "than fitting one."
)

# ---------------------------------------------------------------------------
# (4) THE PROVENANCE CORRECTION -- the ladder is ONE experiment, not four
# ---------------------------------------------------------------------------
#: FIT_HOLDOUT_DECLARATION.md Amendment 17 clause 1. Kang et al. 2026's Table S4
#: -- the repo's only sulfur temperature ladder -- IS Zhai et al. 2023 (Food
#: Chem. 404:134420, online 28 Sep 2022) Table 1, TTCA columns, published 3.5
#: years earlier. Over the 102 shared cells (34 compounds x 3 temperatures) the
#: two agree cell for cell to the last printed decimal in 101 of 102, all six
#: class subtotals and all three grand totals agree EXACTLY -- including an
#: arithmetic error at 120 C reproduced identically -- the printed SDs are
#: identical, and Kang's corresponding author is Zhai's first author.
#: The one discrepant cell is diagnosed: 2-methylthiophene at 120 C reads
#: 1.128 +/- 0.126b in BOTH Zhai columns (a copy-paste), and Kang's 0.000 is the
#: correction. No shipped row uses that cell.
ZHAI_FOODCHEM_PRIMARY_SOURCE = (
    "Zhai, Y. et al. 2023, Food Chemistry 404:134420, doi "
    "10.1016/j.foodchem.2022.134420, Table 1 (online 28 Sep 2022) -- THE "
    "PRIMARY SOURCE. Kang et al. 2026, Sust. Food Technol. 4:3239, Table S4, "
    "which this module originally cited, is a RE-PUBLICATION of the same 102 "
    "cells (zhai2023foodchem_extraction.md sec. 8; k6a_sulfur_ladders_"
    "synthesis.md sec. 2). Any weighting that counted Kang S4 and S5 as two "
    "observations was wrong by 2x."
)
#: THE TIER A LABEL IS WITHDRAWN. `kang2026_SI_extraction.md` sec. 4d assigned
#: Tier A (+/-15%) on the strength of Kang's Table S3 calibration curves. The
#: PRIMARY source states the opposite basis in a printed equation --
#: `Wi = f' * (Ai*ms)/(As*V)` with "the correction factor was 1" -- and
#: identical numbers cannot come from two quantification routes.
ZHAI_LADDER_QUANTIFICATION = (
    "SINGLE-IS SEMI-QUANT, f' = 1, n = 1 pot (1,2-dichlorobenzene, 0.054 ug). "
    "TIER B: x/ 2 to x3 on absolute MAGNITUDE. Ea and SHAPE are unaffected, "
    "because a constant unknown response factor cancels in a ratio and "
    "therefore in an Arrhenius slope. THIS LICENSES within-study shape and "
    "ratio use and FORBIDS absolute yields, cross-compound magnitude "
    "comparison, cross-paper magnitude comparison, and class subtotals read as "
    "physical totals. The six shipped level rows already carry sigma_log = 0.4 "
    "= log10(2.51), which is INSIDE the x/ 2-3 interval, so this re-band is a "
    "label correction and costs no value change; widening a sigma to "
    "manufacture a visible edit would be worse, not better."
)
#: The switch-on the ladder was thought to show does not exist, and the failure
#: is recorded here so no later wave re-fits a two-regime barrier to it.
ZHAI_SWITCH_ON_RETIRED = (
    "RETIRED, NOT FAILED (Amendment 17 clause 2). The MFT/FFT 'switch-on "
    "between 120 and 140 C' is a HOLD-TIME ARTEFACT of a single 120-minute "
    "slice through curves that peak and turn over. It vanishes at 40, 80 and "
    "180 min in the same pot; it vanishes when an equimolar sugar is added "
    "(MFT 1.44x/1.30x, FFT 1.35x/0.96x in the +Xyl arm of the same table); and "
    "Wang 2026's independent 5-rung ladder has both thiols PEAK MID-LADDER AND "
    "FALL. A fixed-time slice across a family of curves whose maxima move left "
    "as temperature rises measures where each curve sits relative to its own "
    "peak, not an activation energy. What survives, and is better supported: "
    "above ~120-140 C thiol consumption outruns thiol formation."
)

# ===========================================================================
# 3b. B2.1's NEW MEASURED CONSTRAINTS -- Kang 2026 SI and Kumazawa 2003
# ===========================================================================
#: Kang 2026 SI sec. 6c: peak free cysteine is 1.63 mmol/L (140 C, 40 min)
#: against 10 mmol/L TTCA charged, so at most 16.3 mol% of the cysteine moiety
#: EVER appears as free cysteine. FIT_HOLDOUT_DECLARATION.md Amendment 5 makes
#: this a ONE-SIDED FIT bound. It is a mass-balance ceiling, not a level: it is
#: exactly the kind of constraint a network can violate silently.
KANG_TTCA_FREE_CYS_YIELD_CEILING_MOL_PCT = 16.3
KANG_TTCA_CEILING_ANCHOR = (
    "Kang 2026 SI Fig. S3 (digitised, kang2026_SI_extraction.md sec. 6c): peak "
    "free Cys 1.63 mmol/L at 140 C / 40 min against 10 mmol/L TTCA loaded. "
    "The pool RISES THEN FALLS at every temperature (t_max ~100 / ~57 / ~40 min "
    "at 100 / 120 / 140 C), so free cysteine is a transient intermediate and "
    "not an accumulating product."
)

#: Kang 2026 SI sec. 6b: the FIRST measured activation energy for cysteine
#: consumption under Maillard conditions. FIT per Amendment 5.
#: Digitised from Fig. S4 at 400 dpi with exact axis calibration; all three
#: curves start at exactly 10.00 mmol/L at t = 0, which validates the axis.
KANG_EA_FREE_CYS_DEPLETION_KJ_MOL = 55.1
#: B8 / Amendment 17: THE PRIMARY SOURCE OF THIS NUMBER IS ZHAI, NOT KANG, and
#: the two are not independent. Zhai Food Chem p. 6 PRINTS the conversions in a
#: sentence -- 16.2% / 38.8% / 62.2% at 100 / 120 / 140 C, Cys 10 mmol/L, pH 7,
#: 120 min -- which as pseudo-first order give k = 1.473e-3 / 4.092e-3 /
#: 8.107e-3 /min and Ea = 54.7 kJ/mol (R^2 0.996, legs 62.3 / 46.2).
#: `kang2026_SI_extraction.md` sec. 6b derived 55.1 by DIGITISING Kang's Fig. S4
#: -- the same experiment, one publication downstream, read off a figure instead
#: of a sentence. The 0.7% agreement is a VALIDATION OF THE DIGITISATION, not an
#: independent confirmation, and `kinetic_core_b2_2_diagnosis.md` K-4's phrase
#: "Kang's own immovable 55.1" should read "Zhai's".
#: B8 does NOT move the shipped value: 54.7 and 55.1 are two readings of one
#: measurement and swapping them would be churn, not a correction.
ZHAI_EA_FREE_CYS_DEPLETION_PRIMARY_KJ_MOL = 54.7
ZHAI_FREE_CYS_DEPLETION_NOTE = (
    "PRIMARY-SOURCE VALUE 54.7 kJ/mol (Zhai Food Chem 2023 p. 6, printed "
    "conversions 16.2 / 38.8 / 62.2% at 100 / 120 / 140 C; R^2 0.996). The "
    "shipped 55.1 is the digitised reading of the SAME experiment in Kang 2026 "
    "Fig. S4 and is kept. NOTE THE CURVATURE: legs 62.3 then 46.2 -- the "
    "cysteine sink DECELERATES across the very interval where the thiols were "
    "said to accelerate, so a single 55.1 is a 40 K average of a mildly curved "
    "line whose curvature has the OPPOSITE SIGN to the retired switch-on."
)
#: FIT ROWS ONLY. The 140 C rung of both ladders is the declared GATING
#: HOLD-OUT (Amendment 5) and is deliberately ABSENT from this file and from
#: everything under src/. tests/unit/test_kinetic_core_b2_1.py greps for its
#: literals.
KANG_CYS_DEPLETION_K_PER_MIN: Mapping[float, float] = {
    100.0: 1.472e-3, 120.0: 4.078e-3,
}
KANG_CYS_CONVERSION_AT_120_MIN: Mapping[float, float] = {
    100.0: 0.162, 120.0: 0.387,
}
KANG_CYS_ANCHOR = (
    "PRIMARY SOURCE Zhai et al. 2023 Food Chem. 404:134420 p. 6 (printed "
    "conversions); re-published as -- and originally cited here from -- "
    "Kang 2026 SI Fig. S4, digitised (kang2026_SI_extraction.md sec. 6b): free "
    "cysteine 10 mmol/L, pH 7, sealed, 0-120 min at 100/120/140 C. Arrhenius "
    "over all three points gives Ea = 55.1 kJ/mol, A = 8.0e4 /min, R^2 0.994, "
    "and refitting on the 0-20 min INITIAL rates instead gives 51 kJ/mol, so "
    "the value is robust to the fitting window at the +/-10% level."
)
KANG_CYS_CAVEATS = (
    "(1) THE DECAYS ARE NOT CLEANLY FIRST ORDER: 120 and 140 C decelerate (the "
    "0-20 min k is 1.6-1.8x the 120-min effective value) and 100 C mildly "
    "accelerates, so the endpoint k is an average over a non-exponential curve. "
    "(2) system_identity: free_Cys at 85% confidence -- the caption says only "
    "'Cys compound', and the reading that Fig. S4 tracks the cysteine moiety "
    "still bound in TTCA cannot be excluded from the published record. "
    "(3) DIGITISED, not printed. (4) n is not stated for this figure. "
    "(5) 140 C is NOT used in this module's fit: only the 100 and 120 C rows "
    "are FIT, because 140 C is the declared gating hold-out."
)

#: Kumazawa & Masuda 2003, J. Agric. Food Chem. 51:2674, Figure 3 and Table 3:
#: 2-furfurylthiol SURVIVAL after 121 C / 10 min in citrate/phosphate buffer,
#: seven pH levels. FIT per FIT_HOLDOUT_DECLARATION.md Amendment 4.
#: These are RESIDUAL RATIOS (heated / unheated half of the SAME solution), so
#: they are ratio-scale and immune to the absolute-calibration defects that
#: qualify most of the corpus's level rows.
KUMAZAWA_FFT_RESIDUAL_FRACTION: Mapping[float, float] = {
    3.0: 0.995, 4.0: 0.962, 5.0: 0.891, 5.4: 0.795, 6.0: 0.451, 6.4: 0.11,
}
KUMAZAWA_ANCHOR = (
    "Kumazawa & Masuda 2003, J. Agric. Food Chem. 51:2674-2678, Figure 3 "
    "(p. 2677), 900 dpi digitisation cross-validated against the paper's own "
    "Table 3 peak-area ratios at all seven pH levels (kumazawa2003_extraction.md "
    "secs. 5, 6.1). 1 ppm 2-furfurylthiol in mixed citric acid / Na2HPO4 "
    "buffer, canned WITHOUT deoxidisation, 121 C / 10 min."
)
KUMAZAWA_CAVEATS = (
    "(1) THE pH LABELS ARE ROOM-TEMPERATURE METER READINGS. The paper says so "
    "and applies no hot-pH correction; citrate/phosphate shifts with "
    "temperature, so the pH during the 121 C hold is unknown. This module's pKa "
    "DO move with temperature while the label does not, which is a real and "
    "unquantified asymmetry -- recorded, not corrected. "
    "(2) NO ERROR BARS, NO SD, NO n ANYWHERE in Table 2, Table 3 or Figure 3. "
    "(3) The pH 7.0 point is left-censored (<0.01 area units, ~0.1% residual) "
    "and is therefore EXCLUDED from the objective rather than fitted as a level. "
    "(4) The paper's own mass balance does NOT close -- the volatile degradation "
    "products recovered are less than the 2-furfurylthiol lost -- which is why "
    "the module needs a non-disulfide sink (k_thiolate_loss) and not only the "
    "dimerisation channel. "
    "(5) Table 2 and Figure 3 disagree at pH 6.0 (48.3% vs 45.1%); the "
    "digitised Figure 3 value is used, and the 3.2-point gap is the honest "
    "measure of this row's precision."
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
    # ---- B2.1: the parallel routes that break the decade-per-pH slope ------
    "k_arp_tdp_th": (
        "pentose Amadori compound -> 3-deoxypentosone (UNCATALYSED)", 1,
        "B2.1, NEW, AND IT IS THE pH FIX RATHER THAN A NEW CHEMICAL CLAIM. B2 "
        "gave the fed Amadori exactly one route to the 3-deoxyosone and made it "
        "fully acid-catalysed, so its product (furfural, and hence FFT) was a "
        "decade per pH unit BY CONSTRUCTION. The uncatalysed 1,2-enolisation is "
        "the same reaction without the proton, it is what keeps caramelisation "
        "running in neutral and alkaline solution, and the module already "
        "carries its free-sugar analogue (k_pent_caramel). With both routes "
        "present the OBSERVABLE pH slope is sub-decade and its crossover is the "
        "ratio of two ordinary rate constants fitted on FIT rows -- not a pH "
        "parameter. Bounded to allow ~zero, so the data may reject it.", ""),
    "k_arp_dpo_th": (
        "pentose Amadori compound -> 1-deoxypentosone (UNCATALYSED)", 1,
        "The base-lane partner of k_arp_tdp_th, for the same reason and with the "
        "same freedom to be rejected.", ""),
    "k_ddp_mft_hs": (
        "1,4-dideoxypentosone + HS- -> MFT", 2,
        "B2.1, NEW. THE SULFIDE NUCLEOPHILE HAS TWO PROTONATION STATES AND B2 "
        "CARRIED ONLY ONE. Taking the reactive species to be the neutral "
        "molecule alone is a mechanism claim that runs against the elementary "
        "chemistry -- HS- is by far the better nucleophile and H2S dominates "
        "only by abundance below the pKa -- and it is half of why B2's MFT and "
        "FFT collapsed at high pH. The two branches are separate steps with "
        "separate constants, so their crossover is a rate ratio and not a pH "
        "parameter. Bounded to allow ~zero.", "hs_anion"),
    "k_fur_fft_hs": (
        "2-furaldehyde + HS- -> FFT", 2,
        "The FFT partner of k_ddp_mft_hs.", "hs_anion"),
    # ---- B2.1: Kang 2026's fed TTCA intermediate ---------------------------
    "k_ttca_cys": (
        "TTCA -> free cysteine + aldopentose (ring-opening / retro-condensation)",
        1,
        "B2.1, NEW. Kang 2026 Fig. S3 measures free cysteine RISING THEN FALLING "
        "at all three temperatures with the maximum moving earlier as T rises "
        "(t_max ~100 / ~57 / ~40 min at 100 / 120 / 140 C) -- textbook "
        "consecutive A -> B -> C kinetics, and the dossier's explicit warning is "
        "that 'any model treating released Cys as an accumulating pool is "
        "wrong'. This is the release half; the module's existing cysteine sinks "
        "are the consumption half, so the transient emerges rather than being "
        "imposed.", ""),
    "k_ttca_deg": (
        "TTCA -> 1-deoxypentosone + fragments, WITHOUT liberating free cysteine",
        1,
        "B2.1, NEW, and it is what the 16.3 mol% ceiling is made of. Kang's peak "
        "free cysteine is 1.63 mmol/L against 10 mmol/L TTCA charged, so at the "
        "most favourable temperature and time at most 16.3% of the cysteine "
        "moiety EVER appears as free cysteine and >=84% leaves by a route that "
        "does not pass through it. A model with only the release route would "
        "have to violate that bound. The ceiling is a declared one-sided FIT "
        "row and this step is how the network can satisfy it.", ""),
    "k_cys_thermal": (
        "L-cysteine -> fragments, WITHOUT releasing sulfide", 1,
        "B2.1, NEW, and the one fitted step in this module with a MEASURED "
        "activation energy of its own (Kang 2026, 55.1 kJ/mol; see "
        "MEASURED_EA_OVERRIDES). Zheng & Ho's H2S-release constants account for "
        "only ~2.7% cysteine conversion at 100 C / 120 min against Kang's "
        "measured 16.2%, so most of what consumes cysteine in a sealed aqueous "
        "system releases no sulfide. This is that remainder. Bounded to allow "
        "~zero.", ""),
    # ---- B2.1: the pH-dependent thiol sink Kumazawa measures ----------------
    "k_thiolate_loss": (
        "R-SH (as thiolate) -> fragments; the non-disulfide oxidative sink", 1,
        "B2.1, NEW, and it is where the FFT-versus-pH slope actually lives. "
        "Kumazawa 2003 Fig. 3 / Table 3 measure 2-furfurylthiol SURVIVAL across "
        "seven pH levels at 121 C / 10 min -- 99.5 / 96.2 / 89.1 / 79.5 / 45.1 / "
        "11 / <0.5 % at pH 3.0 to 7.0 -- with the difurfuryl disulfide growing "
        "monotonically over the same grid AND a mass-balance shortfall the "
        "authors attribute to a Fenton-type route. B2 had no pH dependence on "
        "ANY consumption channel, so it had to load the entire observed pH "
        "response of FFT onto FORMATION, at a decade per pH unit. Half of that "
        "response is measurably a CONSUMPTION effect. Thiol oxidation proceeds "
        "through the thiolate, so the factor is the thiolate fraction at "
        "reaction temperature and carries no fitted pH number. Bounded to allow "
        "~zero.", "thiolate"),
}

FITTED_SULFUR_KEYS: Tuple[str, ...] = tuple(_FITTED_SULFUR)

#: ===================================================================
#: B2.1 -- POLICY 2 IS NARROWED TO THE FOUR NAMED CHANNELS IT WAS ABOUT
#: ===================================================================
#: B2 gave EVERY consumption constant ea_kj_mol = None and held it at its 145 C
#: value at all temperatures, and reported that as conservatism. It is not
#: conservative; it is a strong claim, and it is the claim that broke every
#: low-temperature hold-out B2 had. `k_fft_decay` was fitted at 0.188 /min at
#: 145 C and then evaluated UNCHANGED in Hofmann's 80 C brew, where it alone
#: predicts eight times the measured loss -- which is most of that row's 18x
#: failure -- and in van Seeventer's 50 C storage, where the same constant
#: destroys 99% of the thiol in a day against a measured 59% and produced a
#: "pass" that B2's own report refused to count.
#:
#: WHAT POLICY 2 ACTUALLY SAYS. Inventory sec. B.7's finding is about FOUR
#: NAMED CHANNELS: four papers, four temperatures, four DIFFERENT dominant
#: sinks, each excluding the others'. Deriving one Arrhenius line through two of
#: THOSE is the prohibited derivation, and it stays prohibited: k_thioether,
#: k_oligomer, k_dimer_*, k_mmft and k_protein_ss still carry no activation
#: energy and no code path can pair them.
#:
#: The residual `*_decay` lumps are not those channels. They are unassigned
#: sinks fitted to close a mass balance at one temperature, and nobody has ever
#: measured them at any temperature, so there is no cross-lab pairing to
#: prohibit. Asserting that they are TEMPERATURE-INDEPENDENT -- that a thermal
#: degradation runs exactly as fast at 50 C as at 145 C -- is a claim with no
#: evidence behind it and with a known direction of error. B2.1 therefore gives
#: them the lumped formation Ea, which says they are thermal chemistry like
#: everything else in the network. That is a weaker claim than B2 made, not a
#: stronger one, and it adds NO new parameter: it reuses the single Ea already
#: being fitted.
NAMED_CHANNEL_KEYS: Tuple[str, ...] = ("k_dimer_mft", "k_dimer_fft", "k_mmft")

#: B2.1: fitted steps whose BARRIER is measured even though their RATE is not.
#: These do NOT receive the lumped formation Ea; they receive the measurement,
#: and the fit cannot move it. This is the mechanism by which "the network's
#: temperature dependence flows through measured Ea where they exist".
#:
#: B8 (Amendment 17 clause 5) adds FOUR entries, on two measurements. Each is a
#: single paper's single step measured at THREE temperatures -- the same class
#: of object as Kang/Zhai's 55.1, and NOT the derivation policy 2 prohibits
#: (one Arrhenius line drawn through two DIFFERENT named channels from two
#: DIFFERENT papers at two DIFFERENT temperatures). `k_mmft`, `k_thioether`,
#: `k_oligomer` and `k_protein_ss` still carry no activation energy and no code
#: path may pair them; that prohibition is untouched.
MEASURED_EA_OVERRIDES: Mapping[str, float] = {
    "k_cys_thermal": KANG_EA_FREE_CYS_DEPLETION_KJ_MOL,
    # B8: the oxidative decay family. Zhang 2026 k17, R^2 0.971.
    "k_dimer_mft": ZHANG_EA_THIOL_TO_DISULFIDE_KJ_MOL,
    "k_dimer_fft": ZHANG_EA_THIOL_TO_DISULFIDE_KJ_MOL,
    # B8: the fed-Amadori enolisations. Zhang 2026 k16, R^2 1.000.
    "k_arp_dpo": ZHANG_EA_CYS_AMADORI_TO_ALPHA_DC_KJ_MOL,
    "k_arp_tdp": ZHANG_EA_CYS_AMADORI_TO_ALPHA_DC_KJ_MOL,
}

#: Which source each override came from, so a value can never be lifted out of
#: this module without its conditions travelling with it.
MEASURED_EA_OVERRIDE_ANCHORS: Mapping[str, str] = {
    "k_cys_thermal": KANG_CYS_ANCHOR,
    "k_dimer_mft": ZHANG_DISULFIDE_ANCHOR,
    "k_dimer_fft": ZHANG_DISULFIDE_ANCHOR,
    "k_arp_dpo": ZHANG_AMADORI_ANCHOR,
    "k_arp_tdp": ZHANG_AMADORI_ANCHOR,
}

#: Measured barriers with NO shipped site in this network. Recorded so the
#: measurement is not lost and is not forced onto a step it does not describe.
MEASURED_EA_WITHOUT_A_SITE: Mapping[str, Tuple[float, str]] = {
    "cys_plus_sugar_to_cys_amadori": (
        ZHANG_EA_CYS_TO_CYS_AMADORI_KJ_MOL, ZHANG_CYS_AMADORI_FORMATION_NO_SITE,
    ),
}

#: The unassigned sinks: fitted lumps, no measurement at any temperature, now
#: sharing the one lumped Ea rather than being pinned to 145 C.
UNASSIGNED_SINK_KEYS: Tuple[str, ...] = (
    "k_mft_decay", "k_fft_decay", "k_dimer_decay", "k_nf_decay", "k_fur_decay",
    "k_h2s_loss", "k_osone_decay", "k_thiol_decay", "k_thiolate_loss",
)

#: ===================================================================
#: B2.2 -- THE DECAY LUMPS GET THEIR OWN BARRIERS, IN TWO NAMED FAMILIES
#: ===================================================================
#: B2.1's own diagnosis (results/validation/kinetic_core_b2_1_diagnosis.md sec.
#: 2a) localised its 140 C failure precisely: giving the `*_decay` lumps the
#: LUMPED FORMATION Ea fixed eleven low-temperature rows and broke the top of
#: the ladder, because a sink that accelerates exactly as fast as its source
#: keeps running at full speed once the pot is precursor-exhausted. Its own
#: closing sentence names the fix -- "the `*_decay` lumps should not share the
#: FORMATION barrier. They need their own".
#:
#: THE CONSTRAINT ON DOING IT: a barrier can only be given to a family the
#: DECLARED FIT ROWS can identify, and identification needs the same sink
#: measured at two temperatures. Exactly two families qualify, and the third
#: group is deliberately left alone.
#:
#:   FAMILY A -- FREE-THIOL SINKS. `k_fft_decay`, `k_mft_decay`,
#:   `k_thiol_decay` (the shared minor-thiol lump) and `k_thiolate_loss`. ONE
#:   barrier, shared, because the chemistry is shared: every one of them is the
#:   loss of a free -SH group, and the corpus resolves them at no temperature
#:   separately. WHAT IDENTIFIES IT: Kumazawa 2003's pH grid measures
#:   2-furfurylthiol SURVIVAL at 121 C WITH NO FORMATION PRESENT -- a direct
#:   read of the sink -- while Hofmann 1998's pH-5 SIDA anchors pin the same
#:   thiols' NET level at 145 C and Kang's declared 100/120 C rungs pin them in
#:   a third pot. Three temperatures, one of them formation-free.
#:
#:   FAMILY B -- CARBONYL-POOL SINKS. `k_fur_decay`, `k_osone_decay` and
#:   `k_nf_decay`. ONE barrier, shared, because all three are the thermal
#:   degradation of a reactive carbonyl in the same pot and no fit row resolves
#:   them apart. WHAT IDENTIFIES IT: Kang's FURFURAL ladder -- the dossier's
#:   own "best-behaved Arrhenius series in the table" -- at the declared 100
#:   and 120 C rungs, against Hofmann 1998 Table 5's in-situ furfural and
#:   norfuraneol levels at 145 C.
#:
#:   NOT GIVEN A BARRIER: `k_dimer_decay` and `k_h2s_loss`. They keep the
#:   lumped formation Ea, i.e. B2.1's setting, UNCHANGED. Nothing in the fit
#:   panel measures a disulfide's or a sulfide's decay at two temperatures --
#:   the dimer rows are single-temperature RATIOS and no H2S is measured
#:   anywhere in this corpus (kang2026_SI_extraction.md sec. 6a) -- so a free
#:   barrier for either would be a parameter identified by nothing. Leaving
#:   them where they are is the conservative direction and it costs no claim.
#:
#: NET COST: TWO new fitted numbers. Not nine, and not one per lump.
#:
#: ===================================================================
#: B8 -- THE T-STRUCTURE. THE THIOL SINK IS TWO CHANNELS, NOT ONE.
#: ===================================================================
#: Amendment 17 clause 4. Wave K6a produced four independent measurements of
#: thiol decay and NONE is compatible with 248 kJ/mol. What they show is not one
#: sink with a better number but TWO channels of different chemistry, and B8
#: gives each its own barrier from its own measurement:
#:
#:   COVALENT CAPTURE -- the FAMILY A lumps above (`k_mft_decay`,
#:   `k_fft_decay`, `k_thiol_decay`, `k_thiolate_loss`). STILL FITTED, but now
#:   inside Gigl's measured band (7, 102) with a prior centre of 60.2 instead of
#:   the arbitrary (20, 250). See `DECAY_EA_BOUNDS_BY_FAMILY`.
#:
#:   OXIDATIVE DIMERISATION -- `k_dimer_mft` and `k_dimer_fft`. NOT FITTED AT
#:   ALL: they take Zhang 2026 k17's MEASURED 122.2 kJ/mol through
#:   `MEASURED_EA_OVERRIDES`, and the fit cannot move it. Before B8 they carried
#:   no barrier and were held at their 145 C value at every temperature.
#:
#: THE RESULTING PICTURE, and it is coherent where the B2.2 fit's was not:
#:   thiol FORMATION, aqueous      ~30-75   (Zhai, Feng, Meng)
#:   thiol DECAY, covalent capture ~60      (range 7-102, Gigl)
#:   thiol DECAY, oxidative dimer  ~122     (Zhang k17)
#:   precursor supply (Amadori)    ~86-106  (Zhang k15/k16 legs)
#: Decay barriers are the SAME ORDER as formation -- somewhat higher, not higher
#: by 170 kJ/mol. That is precisely the configuration that produces a
#: PEAK-AND-DECLINE: at low temperature formation wins, at high temperature
#: decay catches up, and above ~140 C precursor supply saturates while
#: consumption keeps accelerating. The corpus measures that turnover in three
#: laboratories and the pre-B8 core could not produce it at all.
#:
#: WHY A LUMPED SINK INFLATES ITS OWN FITTED BARRIER, on the record because it
#: is the failure mode that produced the 248: Zhang's k18 (alpha-DC ->
#: melanoidin) fits Ea = 118.0 from a FIFTH-ORDER rate law, against 28.3 and
#: 30.6 measured directly for browning by two other labs -- a 4x inflation
#: fully explained by four unmodelled co-reactant concentrations riding inside
#: one "barrier".
DECAY_FAMILY_THIOL_SINK: Tuple[str, ...] = (
    "k_mft_decay", "k_fft_decay", "k_thiol_decay", "k_thiolate_loss",
)
DECAY_FAMILY_CARBONYL_SINK: Tuple[str, ...] = (
    "k_fur_decay", "k_osone_decay", "k_nf_decay",
)
#: The lumps that keep the lumped FORMATION Ea, exactly as in B2.1.
DECAY_KEYS_ON_FORMATION_EA: Tuple[str, ...] = ("k_dimer_decay", "k_h2s_loss")

DECAY_FAMILIES: Mapping[str, Tuple[str, ...]] = {
    "thiol_sink": DECAY_FAMILY_THIOL_SINK,
    "carbonyl_sink": DECAY_FAMILY_CARBONYL_SINK,
}

DECAY_FAMILY_IDENTIFYING_ROWS: Mapping[str, str] = {
    "thiol_sink": (
        "Kumazawa 2003 Fig. 3 (2-furfurylthiol survival at 121 C, formation "
        "absent, six declared-FIT pH levels) + Hofmann 1998 T1/T2 pH-5 SIDA "
        "anchors at 145 C + Kang 2026 SI T-S4 MFT/FFT at the declared 100 and "
        "120 C rungs."
    ),
    "carbonyl_sink": (
        "Kang 2026 SI T-S4 furfural at the declared 100 and 120 C rungs (the "
        "dossier's own best-behaved Arrhenius series) + Hofmann 1998 T5 "
        "in-situ furfural and norfuraneol at 145 C + Hofmann 1998 T3/T4 fed "
        "furfural and fed norfuraneol."
    ),
}

#: Search bounds on the two decay barriers. Deliberately as wide as the
#: formation barrier's, and STARTED FROM RANDOM POINTS INSIDE THEM, so that
#: landing above or below the formation Ea is a RESULT. The direction matters:
#: a decay barrier BELOW the formation barrier makes the sinks relatively
#: slower as temperature rises, which is what B2.1's diagnosis says the 140 C
#: rung needs; a decay barrier ABOVE it makes them relatively faster, which is
#: what a flat 100->120 C yield ladder could instead demand. The fit is not
#: told which, and neither bound is tightened toward either answer.
DECAY_EA_BOUNDS: Tuple[float, float] = (20.0, 250.0)

#: B8 -- PER-FAMILY BANDS. `DECAY_EA_BOUNDS` above is B2.2's arbitrary-wide pair
#: and is KEPT UNCHANGED so that every B2.2/B2.3/B2.4 artefact stays exactly
#: reproducible and the tests that pin those runs against it keep meaning what
#: they meant. From B8 the fit uses the per-family bands below.
#:
#: `thiol_sink` (7, 102) is NOT a tightened search bound. It is Gigl 2021's
#: DEFENSIBLE RANGE for the covalent-capture barrier, and the fact that its
#: CEILING sits 2.1x below B2.3's landed 216.1 is the point: the incumbent value
#: is refuted by measurement. A B8 fit that presses 102 is reporting a residual
#: conflict between the objective and the measurement, and the report says so
#: rather than widening the band to hide it.
#:
#: `carbonyl_sink` is UNCHANGED at (20, 250). K6a measured no carbonyl-sink
#: barrier, so narrowing it would be inventing one. It landed at 249.9 in B2.3
#: -- pressed against its ceiling -- and that remains an open, reported defect
#: of exactly the kind Zhang's k18 explains: a lumped multi-order sink inflates
#: its own fitted barrier. B8 does not touch it.
DECAY_EA_BOUNDS_BY_FAMILY: Mapping[str, Tuple[float, float]] = {
    "thiol_sink": GIGL_EA_COVALENT_CAPTURE_RANGE_KJ_MOL,
    "carbonyl_sink": DECAY_EA_BOUNDS,
}

#: The prior centre each family's band is organised around, where one exists.
DECAY_EA_PRIOR_CENTRE: Mapping[str, Optional[float]] = {
    "thiol_sink": GIGL_EA_THIOL_COVALENT_CAPTURE_KJ_MOL,
    "carbonyl_sink": None,
}


def decay_ea_bounds(family: str) -> Tuple[float, float]:
    """The declared search band for one decay family (B8 and later)."""
    return DECAY_EA_BOUNDS_BY_FAMILY.get(family, DECAY_EA_BOUNDS)


#: Kept for the B2 API: every fitted step that consumes rather than forms.
CONSUMPTION_KEYS: Tuple[str, ...] = NAMED_CHANNEL_KEYS + UNASSIGNED_SINK_KEYS

#: The steps that get NO activation energy at all -- policy 2, narrowed AGAIN
#: in B8. `k_dimer_mft` and `k_dimer_fft` leave this tuple because Zhang 2026
#: k17 MEASURES their barrier (one paper, one step, three temperatures, R^2
#: 0.971). `k_mmft` stays: nothing measures a thiol + methanethiol condensation
#: at any second temperature, so it is still held at 145 C.
#: NAMED_CHANNEL_KEYS itself is UNCHANGED -- policy 2's prohibition on pairing
#: the four named channels into one Arrhenius line is about a DERIVATION, not
#: about membership of this tuple, and it still stands.
NO_EA_KEYS: Tuple[str, ...] = tuple(
    k for k in NAMED_CHANNEL_KEYS if k not in MEASURED_EA_OVERRIDES
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
        consumption = key in NO_EA_KEYS
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
            # B2.1: the fit panel now reaches down to 100 C (Kang's TTCA ladder)
            # and through 121 C (Kumazawa's pH grid), so the declared window is
            # widened to the union the fit actually spans.
            t_range=(100.0, 145.0),
            rate_transfer="not_licensed",
            channel=(
                "fitted_consumption" if key in CONSUMPTION_KEYS else ""
            ),
            ph_factor_kind=ph_kind,
            flags=(
                ("fitted_here", "no_literature_value")
                # B8: the barrier flag now says which of THREE regimes a step is
                # in, because there are three. Before B8 there were two and the
                # tuple was chosen by NO_EA_KEYS membership alone, which would
                # now mis-label the two dimer channels as sharing the lumped
                # formation Ea when they carry a measurement.
                + (("no_Ea_available", "temperature_held_fixed")
                   if consumption
                   else (("measured_Ea_override", "fit_cannot_move_this_barrier")
                         if key in MEASURED_EA_OVERRIDES
                         else ("shares_lumped_formation_Ea",)))
                + ((f"b2_2_decay_family:{decay_family_of(key)}",)
                   if decay_family_of(key) else ())
            ),
            note=why,
        )
    return out


def decay_family_of(key: str) -> Optional[str]:
    """Which B2.2 decay family a fitted key belongs to, or None."""
    for family, keys in DECAY_FAMILIES.items():
        if key in keys:
            return family
    return None


def with_fitted_sulfur(
    fitted_log10k: Mapping[str, float],
    lumped_formation_ea: float,
    decay_ea: Optional[Mapping[str, float]] = None,
) -> Dict[str, SulfurParameter]:
    """
    Populate the fitted steps from ``{key: log10 k_ref at 145 C}``.

    Formation steps receive ``lumped_formation_ea``. A step with a MEASURED
    barrier receives the measurement and the fit cannot move it. A step in one
    of the two B2.2 decay families receives that family's own fitted barrier;
    if ``decay_ea`` is not supplied the families fall back to
    ``lumped_formation_ea``, which is exactly B2.1's behaviour, so every B2.1
    artefact stays reproducible. The three NAMED CHANNELS still receive None
    and are held at 145 C -- that asymmetry is policy 2, unchanged.
    """
    out = sulfur_placeholders()
    families = dict(decay_ea or {})
    for key, log10k in fitted_log10k.items():
        if key not in out:
            raise KeyError(f"{key!r} is not a fitted sulfur step")
        family = decay_family_of(key)
        if key in NO_EA_KEYS:
            ea = None
        elif key in MEASURED_EA_OVERRIDES:
            ea = float(MEASURED_EA_OVERRIDES[key])
        elif family is not None and family in families:
            ea = float(families[family])
        else:
            ea = float(lumped_formation_ea)
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
            "equilibrium_constant_M_inv_measured_range": list(
                THIOL_QUINONE_EQUILIBRIUM_M_INV
            ),
            "used_L_per_mmol_at_19_4C": THIOL_QUINONE_EQUILIBRIUM_USED_L_PER_MMOL,
            "delta_H_kJ_mol": STACK_DELTA_H_ADDUCT_KJ_MOL,
            "K_is_temperature_dependent_and_falls_with_T": True,
            "source_anchor": THIOL_QUINONE_EQUILIBRIUM_ANCHOR,
            "transfer_caveats": STACK_TRANSFER_CAVEATS,
            "reverse_rate_is_derived_not_fitted": True,
            "b2_contradiction_corrected": (
                "B2 used K = 10^2.5 = 316 M^-1 on the stated authority of Stack "
                "2018. Stack 2018 measures 5.64-10.57 M^-1. B2's value is 56x "
                "too large and is not traceable to the source."
            ),
        },
        "ph_speciation_at_reaction_temperature": {
            "pKa_25C": dict(PKA_25C),
            "delta_H_ionisation_kJ_mol": dict(DELTA_H_IONISATION_KJ_MOL),
            "pKa_at_the_module_reference_T": {
                name: pka_at(name, T_REF_S_K) for name in PKA_25C
            },
            "thiol_pKa_caveat": THIOL_PKA_CAVEAT,
            "provenance": "constant_not_from_corpus_dossier",
        },
        "protein_disulfide_channel": {
            "bracket_L_per_mmol_min": list(
                PROTEIN_DISULFIDE_EXCHANGE_BRACKET_L_PER_MMOL_MIN
            ),
            "used_L_per_mmol_min": PROTEIN_DISULFIDE_EXCHANGE_USED_L_PER_MMOL_MIN,
            "site_inventory_blg_10pct_mmol_L": (
                PROTEIN_DISULFIDE_SITES_BLG_10PCT_MMOL_L
            ),
            "derivation": PROTEIN_DISULFIDE_DERIVATION,
            "flux_is_identically_zero_in_every_scored_row": True,
        },
        "kang2026_constraints": {
            "free_cys_yield_ceiling_mol_pct": (
                KANG_TTCA_FREE_CYS_YIELD_CEILING_MOL_PCT
            ),
            "ceiling_anchor": KANG_TTCA_CEILING_ANCHOR,
            "ea_free_cys_depletion_kJ_mol": KANG_EA_FREE_CYS_DEPLETION_KJ_MOL,
            "cys_conversion_at_120_min": dict(KANG_CYS_CONVERSION_AT_120_MIN),
            "anchor": KANG_CYS_ANCHOR,
            "caveats": KANG_CYS_CAVEATS,
        },
        "kumazawa2003_ph_grid": {
            "fft_residual_fraction_121C_10min": dict(
                KUMAZAWA_FFT_RESIDUAL_FRACTION
            ),
            "anchor": KUMAZAWA_ANCHOR,
            "caveats": KUMAZAWA_CAVEATS,
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
