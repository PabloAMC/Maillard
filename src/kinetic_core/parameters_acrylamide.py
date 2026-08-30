"""
src/kinetic_core/parameters_acrylamide.py

THE PARAMETER REGISTRY OF THE ACRYLAMIDE / SAFETY MODULE (Build Wave B3).
=========================================================================

Module 3 of ``docs/reference/FIT_HOLDOUT_DECLARATION.md``. Same contract as
B1's ``parameters.py`` and B2's ``parameters_sulfur.py``: every number carries
its value, unit, source anchor down to the table, the conditions it was
measured under, its ``ph_of_measurement`` AND its ``aw_of_measurement``, an
``evidence_class`` and a ``rate_transfer`` licence, and all of it travels into
runtime metadata.

The K3 inventory calls this "the best-parameterised block in the corpus", and
that is true -- seven of the sixteen steps in ``acrylamide.py`` carry a
literature rate constant with a published interval. It is also the block with
the corpus's worst fabrication, and both facts are encoded here.

SIX STANDING POLICIES, ENFORCED HERE
------------------------------------
1. **NO DFT.** Same as B1 and B2. ``assert_no_dft_acrylamide()`` runs at
   import.

2. **THE 129 kJ/mol BARRIER DOES NOT EXIST AND CANNOT ENTER.** The shipped
   ``Ea = 129`` for acrylamide (``src/barrier_constants.py``, ``src/safety.py``)
   is unsourceable: the inventory checked all three Knol papers and found
   maxima of 102 (2005, upper-95 % 116), none at all (2009) and 93 +/- 12
   (2010). Its companion ``A_f = 1.6e13`` carries a fit-circularity signature.
   Both are listed in ``REFUSED_PARAMETERS`` with the reason, and
   ``assert_no_fabricated_barrier()`` fails the import if either value ever
   appears in this registry.

3. **NO pH TERM, AND IN PARTICULAR NO ``_acrylamide_ph_factor``.** The shipped
   FAST lane multiplies its acrylamide rate by a pH factor worth ~2000x at
   pH 5.5, which is the whole of the 480x two-lane contradiction (inventory
   line 206, Z0 #13). Nothing in the Module 3 fit corpus varies pH: Claeys is
   0.05 M citrate at pH 6, the extrusion benchmark is pH 6, De Vleeschouwer is
   a freeze-dried powder with no pH at all, and Knol is unbuffered aqueous.
   One pH, therefore no pH term, therefore no place to put a 2000x factor.

4. **NO WATER-ACTIVITY TERM EITHER -- BUT a_w IS CARRIED ON EVERY ROW.** This
   is the one axis on which the module's own fit corpus is genuinely split:
   Claeys is dilute aqueous (a_w ~ 1.0), De Vleeschouwer is a freeze-dried
   powder at a_w 0.92, and the extrusion benchmark is a_w 0.35. NOTHING in the
   corpus measures the same step at two water activities, so an a_w term would
   be invented. What the module does instead is carry ``aw_of_measurement`` on
   every parameter and emit a WARNING on every run evaluated more than 0.1
   units away from it. The a_w gap is thereby visible in the metadata of every
   prediction rather than absorbed into a fitted constant.

5. **EVERY SPECIES HAS A CONSUMPTION TERM, AND ACRYLAMIDE ESPECIALLY.** The
   old FAST lane's acrylamide had no elimination at all, which is why it was
   ~40x under-responsive to time. Here acrylamide has three sinks: the
   first-order degradation (fitted to three independent measurements of it),
   the MEASURED bimolecular cysteine adduct, and the fitted amine adducts of
   the competitor systems.

6. **COMPETITION IS A MASS-ACTION CHANNEL, NEVER A MULTIPLIER.** There is not
   one "competitor factor", "yield modifier" or per-amino-acid scaling constant
   in this file, and there is no place to put one: the registry contains only
   rate constants on named, balanced reactions. A competing amino acid changes
   the acrylamide yield in exactly two ways, both of them measured mechanisms:
   it consumes the SHARED GLUCOSE POOL (De Vleeschouwer II measures this
   channel directly for cysteine, k_INT2 = 0.26 M^-1 min^-1) and it SCAVENGES
   ACRYLAMIDE as a Michael acceptor (measured for cysteine as k_E2, fitted for
   the amines under the family licence the inventory grants).

UNITS: mmol/L, minutes, Kelvin, kJ/mol.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, replace
from typing import Any, Dict, Mapping, Optional, Tuple

from .parameters import EVIDENCE_CLASSES, R_KJ

#: Reference temperature of the reparameterised Arrhenius form: 160 C.
#: This is not a choice. It is the T_ref that Claeys 2005 Table 2, De
#: Vleeschouwer 2009 Part I Table 3 and De Vleeschouwer 2009 Part II Table 3
#: are ALL printed at, and it sits in the middle of the 120-200 C window every
#: one of them was measured over.
T_REF_A_K = 433.15

UNIT_BY_ORDER: Mapping[int, str] = {1: "1/min", 2: "L/(mmol*min)"}


def per_molar_to_per_mmol(value: float) -> float:
    """
    M^-1 min^-1  ->  L/(mmol*min).

    Every second-order constant in the acrylamide corpus is printed per MOLAR
    and this module runs per MILLIMOLAR, so this factor of 1000 is applied
    exactly once, here, and is unit-tested. It is called out because a
    mM/M confusion is one of the defects the audit found in the shipped lane.
    """
    return float(value) * 1.0e-3


@dataclass(frozen=True)
class AcrylamideParameter:
    """
    One acrylamide-lane rate constant, with everything needed to defend it.

    Differs from B2's ``SulfurParameter`` in two ways, both forced by this
    corpus rather than chosen:

      * ``ea_kj_mol`` is REQUIRED. Unlike the thiol sinks, every step in this
        module has a published activation energy from at least one of three
        labs, so a None here would be an omission rather than a declared gap.

      * ``aw_of_measurement`` is required and is checked at run time, because
        a_w is the axis this module's own fit corpus is split on (policy 4).
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
    aw_of_measurement: Optional[float]
    #: the temperature this k_ref is quoted AT, K
    t_ref_k: float
    temperature_range_c: Tuple[float, float]
    rate_transfer: str
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
            raise ValueError(f"{self.key}: order must be 1 or 2, got {self.order}")
        if self.unit != UNIT_BY_ORDER[self.order]:
            raise ValueError(
                f"{self.key}: order {self.order} requires unit "
                f"{UNIT_BY_ORDER[self.order]!r}, got {self.unit!r}"
            )
        if self.k_ref is not None and self.ea_kj_mol is None:
            raise ValueError(
                f"{self.key}: a POPULATED acrylamide parameter must carry an "
                f"activation energy. Every step in this module has a published "
                f"Ea from at least one of Claeys 2005, De Vleeschouwer 2009 or "
                f"Knol 2005; a missing one would be an omission, not a declared "
                f"gap (contrast B2's thiol sinks, sec. C.1)."
            )

    def k_at(self, temperature_k: float) -> float:
        """k(T) = k_ref * exp(-(Ea/R)(1/T - 1/T_ref)), same unit as k_ref."""
        if self.k_ref is None or self.ea_kj_mol is None:
            raise ValueError(f"{self.key}: parameter is not numerically populated")
        return float(self.k_ref) * math.exp(
            -(float(self.ea_kj_mol) / R_KJ)
            * (1.0 / float(temperature_k) - 1.0 / float(self.t_ref_k))
        )

    def as_metadata(self) -> Dict[str, Any]:
        payload = asdict(self)
        payload["flags"] = list(self.flags)
        payload["temperature_range_c"] = list(self.temperature_range_c)
        payload["dft_derived"] = False
        return payload


# ===========================================================================
# THE MEASURED BACKBONE
# ===========================================================================
# Three sources, all printed at T_ref = 160 C, all in the 120-200 C window.
#
#   Claeys, W., De Vleeschouwer, K. & Hendrickx, M. (2005) "Kinetics of
#   acrylamide formation and elimination in equimolar asparagine/glucose model
#   systems", Biotechnol. Prog. 21:1525-1530, Table 2 p. 1528.
#   FIT (control AND the four competitor systems), declaration D.5.
#
#   De Vleeschouwer, K. et al. (2009) Part I, Food Chem. 114:116-126, Table 3
#   p. 124. FIT for the GLUCOSE, non-italic subset ONLY.
#
#   De Vleeschouwer, K. et al. (2009) Part II, J. Agric. Food Chem. 57:539-546,
#   Table 3 p. 542. FIT for the CYSTEINE system ONLY.
#
# THE COLUMNS THAT ARE HOLD-OUT ARE NOT IN THIS FILE, and the firewall test
# checks that mechanically. See the disclosure block at the bottom.

_DV1_SOURCE = (
    "De Vleeschouwer et al. 2009 Part I, Food Chem 114:116-126, Table 3 p. 124, "
    "GLUCOSE column (non-italic subset)"
)
_DV1_CONDITIONS = (
    "equimolar asparagine + glucose, freeze-dried powder equilibrated to "
    "a_w 0.92 at 4 C (moisture 12.6-17.8 %), ~3 mol/L reactant concentrations, "
    "hermetically sealed inox tubes, 120/140/160/180/200 C, 4 s thermocouple "
    "log, multiresponse fit by the determinant criterion"
)
_DV2_SOURCE = (
    "De Vleeschouwer et al. 2009 Part II, J Agric Food Chem 57:539-546, "
    "Table 3 p. 542, CYSTEINE column"
)
_DV2_CONDITIONS = (
    "as Part I plus equimolar L-cysteine; a_w 0.92 freeze-dried powder "
    "(moisture 19.65 % for the cysteine system), ~3 mol/L, 120-200 C. "
    "[Cys] is measured as cysteine + cystine"
)


def _measured(
    key: str,
    transformation: str,
    k_ref: float,
    k_ci: Optional[float],
    ea: float,
    ea_ci: Optional[float],
    order: int,
    *,
    source: str,
    conditions: str,
    dossier: str,
    ph: Optional[float],
    aw: Optional[float],
    t_range: Tuple[float, float] = (120.0, 200.0),
    rate_transfer: str = "licensed_at_measurement_aw_only",
    flags: Tuple[str, ...] = (),
    note: str = "",
) -> AcrylamideParameter:
    return AcrylamideParameter(
        key=key,
        transformation=transformation,
        k_ref=k_ref,
        ea_kj_mol=ea,
        unit=UNIT_BY_ORDER[order],
        order=order,
        evidence_class="measured_rate",
        source_anchor=source,
        dossier_anchor=dossier,
        conditions=conditions,
        ph_of_measurement=ph,
        aw_of_measurement=aw,
        t_ref_k=T_REF_A_K,
        temperature_range_c=t_range,
        rate_transfer=rate_transfer,
        k_ref_ci95=k_ci,
        ea_ci95_kj_mol=ea_ci,
        flags=flags,
        note=note,
    )


MEASURED_ACRYLAMIDE: Mapping[str, AcrylamideParameter] = {
    # ---- the asparagine trunk branch, De Vleeschouwer I glucose column ------
    "k_asn_glc": _measured(
        "k_asn_glc", "Asn + Glc -> Schiff base (bimolecular initiation)",
        per_molar_to_per_mmol(1.70), per_molar_to_per_mmol(1.05),
        117.5, 25.2, 2,
        source=_DV1_SOURCE + ", k_INTg / Ea_INTg",
        conditions=_DV1_CONDITIONS,
        dossier="k3_final_parameter_inventory.md sec. A.2 line 192",
        ph=None, aw=0.92,
        flags=("the_only_bimolecular_initiation_constant_in_the_corpus",),
        note=(
            "THE reason this module exists in mass-action form. It is the only "
            "genuinely second-order Maillard-initiation constant anywhere in the "
            "corpus, which is what lets the network respond to precursor "
            "CONCENTRATION instead of carrying a fixed yield -- the exact defect "
            "the audit's Z1 sec. E diagnosed as saturation. Printed per MOLAR; "
            "converted here by per_molar_to_per_mmol, once."
        ),
    ),
    "k_int1_acr": _measured(
        "k_int1_acr", "Int1 -> acrylamide (+ unmeasured C7 N1 residue)",
        3.57e-3, 1.38e-3, 159.2, 29.5, 1,
        source=_DV1_SOURCE + ", k_Fref / Ea_F",
        conditions=_DV1_CONDITIONS,
        dossier="k3_final_parameter_inventory.md sec. A.2 line 193",
        ph=None, aw=0.92,
        flags=("printed_unit_is_a_typo_10e-3_mm-1_means_min-1",),
        note=(
            "The source prints the unit as '(10^-3 mm^-1)', which the inventory "
            "identifies as a typo for min^-1 (K1 sec. 2b). Taken as min^-1. The "
            "acrylamide-forming step is FIRST order in the intermediate, which "
            "is why the acrylamide yield is set by the INT1 PARTITION and not by "
            "this constant alone."
        ),
    ),
    "k_asn_asp": _measured(
        "k_asn_asp", "Asn -> Asp + NH3 (deamidation)",
        26.43e-3, 5.76e-3, 105.4, 10.6, 1,
        source=_DV1_SOURCE + ", k_Aspref / Ea_Asp",
        conditions=_DV1_CONDITIONS,
        dossier="k3_final_parameter_inventory.md sec. A.2 lines 195-196",
        ph=None, aw=0.92,
        flags=("sugar_independent_across_all_three_sugars",),
        note=(
            "The competing fate of asparagine, and the most transferable number "
            "in the block: Ea_Asp is 105.4 / 108.3 / 109.4 kJ/mol across "
            "glucose, fructose and sucrose. Only the GLUCOSE column is used "
            "here; the other two are HOLD-OUT and the sugar-independence is "
            "noted from the inventory's verdict text, not from their values."
        ),
    ),
    # ---- the cysteine channels, De Vleeschouwer II cysteine column ---------
    "k_acr_cys": _measured(
        "k_acr_cys", "acrylamide + Cys -> S-(2-carbamoylethyl)cysteine",
        per_molar_to_per_mmol(49.36), per_molar_to_per_mmol(1.18),
        51.3, 1.5, 2,
        source=_DV2_SOURCE + ", k_E2ref / Ea_E2",
        conditions=_DV2_CONDITIONS,
        dossier="k3_final_parameter_inventory.md sec. A.2 lines 190-191",
        rate_transfer="licensed_to_the_thiol_michael_acceptor_family",
        ph=None, aw=0.92,
        flags=("tightest_parameter_in_the_corpus_2.4pct_rse",
               "order_assumed_never_tested_by_the_source"),
        note=(
            "The corpus's tightest parameter (2.4 % RSE) and its only measured "
            "acrylamide-scavenging constant. Ea_E2 = 51.3 is HALF of the "
            "elimination barrier, so the two channels CROSS OVER in temperature "
            "-- a structural fact this network reproduces by construction rather "
            "than by a fitted switch. The inventory's verdict grants the "
            "generalisation licence used for the fitted amine adducts: "
            "'generalises to any thiol/Michael-acceptor pair'. CAVEAT CARRIED: "
            "the second order in this step is ASSUMED by the source and never "
            "tested (K1 sec. 2c)."
        ),
    ),
    "k_cys_sink": _measured(
        "k_cys_sink", "Cys -> unidentified products (the De Vleeschouwer II Cys sink)",
        0.35, 0.01, 110.5, 8.5, 1,
        source=_DV2_SOURCE + ", k_Yref / Ea_Y",
        conditions=_DV2_CONDITIONS,
        dossier="k3_final_parameter_inventory.md sec. A.2 line 199",
        ph=None, aw=0.92,
        flags=("product_not_identified_by_the_source",),
        note=(
            "USE-Q in the inventory: the constant is well determined but 'the "
            "product is not identified'. It is operative here because it is the "
            "cysteine sink measured IN THE SYSTEM THIS MODULE FITS; B2's four "
            "thiol-consumption channels are measured in coffee brew and in "
            "thiamine/xylose pots at 25-120 C and are NOT composed with this "
            "one. Composing them would spend the same cysteine twice."
        ),
    ),
    "k_cys_glc": _measured(
        "k_cys_glc", "Cys + Glc -> melanoidin (LUMPED competitor consumption of sugar)",
        per_molar_to_per_mmol(0.26), per_molar_to_per_mmol(0.02),
        30.3, 1.6, 2,
        source=_DV2_SOURCE + ", k_INT2ref / Ea_INT2",
        conditions=_DV2_CONDITIONS,
        dossier="k3_final_parameter_inventory.md sec. A.2, K1 sec. 2c row k_INT2ref",
        rate_transfer="not_licensed",
        ph=None, aw=0.92,
        flags=("apparent_only_authors_say_not_comparable_to_k_INT",),
        note=(
            "THE MEASURED COMPETITION CHANNEL, and the template for the three "
            "fitted ones. A competing amino acid removes glucose from the "
            "asparagine lane; this is that reaction, measured. The authors' own "
            "qualification travels with it: k_INT2 is APPARENT and 'not "
            "comparable to k_INT', so it is used as the cysteine competitor's "
            "own sugar-consumption rate and never as a second initiation "
            "constant."
        ),
    ),
    "k_asp_sink": _measured(
        "k_asp_sink", "Asp -> unidentified products",
        0.04, 0.01, 97.2, 8.3, 1,
        source=_DV2_SOURCE + ", k_Xref / Ea_X (cysteine column)",
        conditions=_DV2_CONDITIONS,
        dossier="k3_final_parameter_inventory.md sec. A.2, K1 sec. 2c row k_Xref",
        ph=None, aw=0.92,
        flags=("the_only_usable_k_X_in_the_corpus",),
        note=(
            "Aspartic acid's consumption term. Every OTHER k_X in the corpus is "
            "refused: Part I's glucose value is 'Indeterminate' with Ea_X = "
            "668.9 kJ/mol (physically absurd) and the glutamine column inherits "
            "that fixed value. The cysteine column's 0.04 +/- 0.01 with Ea "
            "97.2 +/- 8.3 is the one that was actually estimated."
        ),
    ),
}


# ===========================================================================
# THE PINNED SCHIFF / AMADORI SPLIT FOR THE ASPARAGINE LANE
# ===========================================================================
# Exactly B1's device, for exactly B1's reason, and with B1's number.
#
# De Vleeschouwer's Scheme 4 does not resolve a Schiff base: everything between
# the reactants and the acrylamide-forming step is one lumped 'Int1'. The
# network carries two states anyway, because a later module needs somewhere to
# attach -- but the SPLIT IS NOT PARAMETERISED. The second step's constant is
# DERIVED from the first by a pinned ratio, and because the second step is fast,
# irreversible and sink-free, every condensed molecule reaches Int1 and the pair
# reduces EXACTLY to the measured one-step composite. The ratio's numeric value
# therefore changes nothing that is scored; what it changes is a sub-minute time
# lag at 160 C. That is why it is legitimate to reuse B1's glycine-derived
# ratio here, and why the reuse is flagged `not_licensed` rather than presented
# as an asparagine measurement.
ASN_SCHIFF_AMADORI_SPLIT: Dict[str, Any] = {
    "ratio_amadori_over_schiff_pseudo_first_order": 44.9,
    "amine_loading_mmol_L_for_the_ratio": 200.0,
    "evidence_class": "structural_constant",
    "source_anchor": (
        "REUSED from src/kinetic_core/parameters.SCHIFF_AMADORI_SPLIT, which was "
        "fitted to the Martins GLYCINE system, not to asparagine"
    ),
    "rate_transfer": "not_licensed",
    "why_it_does_not_matter": (
        "The rearrangement is the only fate of the Schiff base in this network, "
        "so the pair's throughput equals the condensation's rate for ANY ratio "
        "above ~1. The ratio sets a time lag (0.07 min at 160 C in the Claeys "
        "system), not a yield. If a later wave gives the Schiff base a sink, "
        "this constant stops being inert and must be re-derived for asparagine."
    ),
    "source_verdict_on_the_split": (
        "REFUSED, twice over. Martins 2005 T1 found that removing the Schiff "
        "intermediate 'fitted the data equally well'; De Vleeschouwer 2009 I "
        "never resolves it at all and calls the lump Int1."
    ),
}


# ===========================================================================
# THE FITTED EXTENSION
# ===========================================================================
# Five named mechanisms with no usable literature rate, in three groups.
#
# GROUP 1 -- the Int1 partition. De Vleeschouwer measures k_M (Int1 -> Int2) at
# 1.23 min^-1, and its own authors mark it "NO PHYSICAL MEANING". It is refused
# (inventory line 200). But Int1 MUST have a non-acrylamide fate: without one
# every molecule that condenses becomes acrylamide and the yield is ~100 %
# instead of the ~0.4-3 % the corpus measures. So the sink exists, is named
# (browning: Int1 into the melanoidin mass pool), and is FITTED.
#
# GROUP 2 -- acrylamide's first-order elimination. THREE independent labs
# measure this same constant and they agree to 1.26x on the rate at 160 C
# (0.1111 Claeys, 0.100 De Vleeschouwer, 0.0881 Knol) while disagreeing by 2x
# on the barrier (167.2 / 113.2 / 85.1 kJ/mol). All three are FIT rows. Fitting
# ONE (k, Ea) pair to all three is how three measurements of one quantity are
# combined; the residuals then show the barrier conflict instead of hiding it
# behind a choice of source.
#
# GROUP 3 -- the competition channels of the three amine competitors. Cysteine's
# two channels are measured (k_cys_glc, k_acr_cys); glutamine's, lysine's and
# alanine's are not, so each gets the SAME TWO NAMED CHANNELS with fitted
# constants and bounds that ALLOW ~ZERO. The acrylamide-scavenging channels
# share the MEASURED Ea_E2 = 51.3 kJ/mol under the inventory's own family
# licence ("generalises to any thiol/Michael-acceptor pair"), so the fit
# determines three pre-exponentials and no new barrier.

#: {key: (transformation, order, why, ea_policy)}
_FITTED_TEMPLATE: Dict[str, Tuple[str, int, str, str]] = {
    "k_int1_mel": (
        "Int1 -> melanoidin (the non-acrylamide fate of the intermediate)", 1,
        "THE PARTITION THAT SETS THE ACRYLAMIDE YIELD. Its literature value "
        "exists and is REFUSED: De Vleeschouwer's k_M = 1.23 min^-1 is one of "
        "the parameters its own authors mark 'NO PHYSICAL MEANING'. Fitted "
        "here against Claeys' dilute-aqueous lumped constants, which is the "
        "only place in the corpus where the partition is observable.",
        "fitted",
    ),
    "k_acr_dp": (
        "acrylamide -> degradation products (first order)", 1,
        "Acrylamide's basic elimination channel. THREE FIT measurements of it "
        "exist (Claeys k_E, De Vleeschouwer k_E, Knol k6 at five "
        "temperatures); one (k, Ea) pair is fitted to all three rather than "
        "one source being picked. The products are NOT identified by any of "
        "them -- Knol's own caveat is that 'the model was not restrained by "
        "experimental data for the products formed in the degradation "
        "reaction' -- so the carbon and nitrogen are routed to the fragment "
        "pools rather than to an invented species.",
        "fitted",
    ),
    "k_gln_glc": (
        "Gln + Glc -> melanoidin (competitor consumption of the shared sugar)", 2,
        "The glutamine analogue of the MEASURED k_cys_glc. Fitted; bounds "
        "allow ~zero.",
        "shared_competitor_ea",
    ),
    "k_lys_glc": (
        "Lys + Glc -> melanoidin (competitor consumption of the shared sugar)", 2,
        "The lysine analogue of the MEASURED k_cys_glc. Fitted; bounds allow "
        "~zero.",
        "shared_competitor_ea",
    ),
    "k_ala_glc": (
        "Ala + Glc -> melanoidin (competitor consumption of the shared sugar)", 2,
        "The alanine analogue of the MEASURED k_cys_glc. Claeys' alanine row "
        "is statistically indistinguishable from the control, so this constant "
        "is expected to come out at the bottom of its range; the bounds allow "
        "that and the fit report says whether it did.",
        "shared_competitor_ea",
    ),
    "k_acr_gln": (
        "acrylamide + Gln -> adduct (unidentified)", 2,
        "The glutamine analogue of the MEASURED k_acr_cys, under the "
        "inventory's family licence. Its BARRIER is not fitted: it is fixed at "
        "the measured Ea_E2 = 51.3 kJ/mol.",
        "measured_ea_e2",
    ),
    "k_acr_lys": (
        "acrylamide + Lys -> adduct (unidentified)", 2,
        "The lysine analogue. Lysine's epsilon-amino group is the strongest "
        "Michael nucleophile among the competitors and Claeys measures the "
        "largest elimination effect for it, so this is the channel the "
        "competitor panel most sharply identifies.",
        "measured_ea_e2",
    ),
    "k_acr_ala": (
        "acrylamide + Ala -> adduct (unidentified)", 2,
        "The alanine analogue. Expected ~zero; the bounds allow it.",
        "measured_ea_e2",
    ),
}

FITTED_ACRYLAMIDE_KEYS: Tuple[str, ...] = tuple(_FITTED_TEMPLATE)

#: Search bounds on log10(k_ref at 160 C). Deliberately wide, and every fit
#: start is random inside them, so agreement with a literature value is a
#: RESULT and not an initialisation. Every competitor channel's lower bound is
#: low enough to be indistinguishable from zero on the fit panel (1e-9
#: L/(mmol*min) against a 10 mmol/L competitor is 1e-8 min^-1, i.e. nothing in
#: 45 minutes), which is what "the data may reject the channel" means.
FITTED_ACRYLAMIDE_BOUNDS_LOG10K: Mapping[str, Tuple[float, float]] = {
    "k_int1_mel": (-6.0, 2.0),
    "k_acr_dp": (-4.0, 1.0),
    "k_gln_glc": (-9.0, -1.0),
    "k_lys_glc": (-9.0, -1.0),
    "k_ala_glc": (-9.0, -1.0),
    "k_acr_gln": (-9.0, 0.0),
    "k_acr_lys": (-9.0, 0.0),
    "k_acr_ala": (-9.0, 0.0),
}

#: The two fitted barriers, and the one that is NOT fitted.
FITTED_ACRYLAMIDE_EA_KEYS: Tuple[str, ...] = (
    "Ea_int1_mel", "Ea_acr_dp", "Ea_competitor_sugar",
)
FITTED_ACRYLAMIDE_EA_BOUNDS: Tuple[float, float] = (20.0, 260.0)

#: The barrier the three fitted acrylamide-scavenging channels are HELD at.
#: Measured, not fitted: De Vleeschouwer II's Ea_E2.
EA_E2_MEASURED_KJ_MOL = 51.3

#: Which fitted Ea each fitted rate uses.
_EA_ASSIGNMENT: Mapping[str, str] = {
    "k_int1_mel": "Ea_int1_mel",
    "k_acr_dp": "Ea_acr_dp",
    "k_gln_glc": "Ea_competitor_sugar",
    "k_lys_glc": "Ea_competitor_sugar",
    "k_ala_glc": "Ea_competitor_sugar",
    "k_acr_gln": "MEASURED_Ea_E2",
    "k_acr_lys": "MEASURED_Ea_E2",
    "k_acr_ala": "MEASURED_Ea_E2",
}

#: The a_w the fitted steps are declared at. The fit panel spans a_w 0.35-1.0
#: and NOTHING in it measures one step at two water activities, so the fitted
#: constants are cross-a_w LUMPS. Declaring them at the Claeys value (the panel
#: rows that actually identify them) with `rate_transfer: not_licensed` is the
#: honest encoding; the alternative -- an a_w exponent -- would be invented.
FITTED_AW_OF_DECLARATION = 1.0


def acrylamide_placeholders() -> Dict[str, AcrylamideParameter]:
    """The eight fitted steps, unpopulated. Integration refuses them as-is."""
    out: Dict[str, AcrylamideParameter] = {}
    for key, (transformation, order, why, ea_policy) in _FITTED_TEMPLATE.items():
        out[key] = AcrylamideParameter(
            key=key,
            transformation=transformation,
            k_ref=None,
            ea_kj_mol=None,
            unit=UNIT_BY_ORDER[order],
            order=order,
            evidence_class="derived_from_fit_data",
            source_anchor=(
                "ESTIMATED IN THIS MODULE by least squares on the declared FIT "
                "rows of docs/reference/FIT_HOLDOUT_DECLARATION.md D.5. No "
                "usable literature value exists (see the note)."
            ),
            dossier_anchor=(
                "docs/reference/FIT_HOLDOUT_DECLARATION.md D.5 (Module 3); "
                "k3_final_parameter_inventory.md sec. A.2"
            ),
            conditions=(
                "fitted across the Module 3 FIT panel: Claeys 2005 dilute "
                "aqueous pH 6 (140-200 C), De Vleeschouwer 2009 I/II a_w 0.92 "
                "powder (120-200 C), Knol 2005 aqueous (120-200 C)"
            ),
            ph_of_measurement=6.0,
            aw_of_measurement=FITTED_AW_OF_DECLARATION,
            t_ref_k=T_REF_A_K,
            temperature_range_c=(120.0, 200.0),
            rate_transfer="not_licensed",
            flags=("fitted_here", f"ea_policy:{ea_policy}"),
            note=why,
        )
    return out


def with_fitted_acrylamide(
    log10_k_ref: Mapping[str, float],
    fitted_ea: Mapping[str, float],
) -> Dict[str, AcrylamideParameter]:
    """
    Populate the fitted steps.

    ``log10_k_ref`` is ``{key: log10 k at 160 C}``; ``fitted_ea`` is
    ``{Ea_int1_mel, Ea_acr_dp, Ea_competitor_sugar}`` in kJ/mol. The three
    acrylamide-scavenging channels do NOT read ``fitted_ea``: they are held at
    the MEASURED Ea_E2, and there is no argument by which a caller can override
    that.
    """
    out = acrylamide_placeholders()
    for key, value in log10_k_ref.items():
        if key not in out:
            raise KeyError(f"{key!r} is not a fitted step; it is measured or unknown")
        ea_key = _EA_ASSIGNMENT[key]
        if ea_key == "MEASURED_Ea_E2":
            ea = EA_E2_MEASURED_KJ_MOL
        else:
            if ea_key not in fitted_ea:
                raise KeyError(f"{key!r} needs fitted activation energy {ea_key!r}")
            ea = float(fitted_ea[ea_key])
        out[key] = replace(
            out[key], k_ref=float(10.0 ** float(value)), ea_kj_mol=float(ea)
        )
    return out


# ===========================================================================
# WHAT IS REFUSED, AND WHY -- carried in the registry, not in a comment
# ===========================================================================

REFUSED_PARAMETERS: Tuple[Dict[str, Any], ...] = (
    {
        "quantity": "acrylamide activation energy, shipped",
        "value": 129.0,
        "unit": "kJ/mol",
        "where_it_lives": "src/barrier_constants.py:274, :493; src/safety.py:790-795",
        "verdict": "REFUSE -- FABRICATION-CLASS, TRIPLE-CONFIRMED, UN-RE-POINTABLE",
        "why": (
            "Attributed to Knol, and absent from all three Knol papers: Knol "
            "2005's maximum is 102 (upper-95 % 116), Knol 2009 publishes no Ea "
            "at all, Knol 2010's maximum is 93 +/- 12. It cannot be re-pointed "
            "to any source because there is no source. It is NOT carried into "
            "this module in any form, and assert_no_fabricated_barrier() fails "
            "the import if it ever appears."
        ),
        "dossier_anchor": "k3_final_parameter_inventory.md sec. A.2 line 205",
    },
    {
        "quantity": "acrylamide Arrhenius prefactor, shipped",
        "value": 1.6e13,
        "unit": "(unstated)",
        "where_it_lives": "src/barrier_constants.py",
        "verdict": "REFUSE -- fit-circularity signature",
        "why": (
            "A_f*exp(-129000/RT) at 160 C = 4.4e-3 against Knol's own k2 = "
            "4.5e-3, i.e. the prefactor was back-solved to reproduce Knol's "
            "rate at Knol's own average temperature using a barrier Knol never "
            "published. The ~2000x that separates the two shipped lanes is "
            "_acrylamide_ph_factor(5.5), not chemistry."
        ),
        "dossier_anchor": "k3_final_parameter_inventory.md sec. A.2 line 206",
    },
    {
        "quantity": "De Vleeschouwer's italicised parameters (k_M, k_B, k_C, k_X and their Ea)",
        "value": None,
        "unit": "various",
        "where_it_lives": "De Vleeschouwer 2009 I Table 3, italic rows",
        "verdict": "REFUSE -- the authors mark them 'NO PHYSICAL MEANING'",
        "why": (
            "Includes Ea_Cg = -6.7 kJ/mol (negative) and Ea_X = 668.9 kJ/mol "
            "(absurd). k_M is the Int1 partition, so refusing it is what forces "
            "k_int1_mel to be FITTED rather than carried. The one exception is "
            "the CYSTEINE column's k_X, which Part II actually estimated "
            "(0.04 +/- 0.01, Ea 97.2 +/- 8.3), and that one is operative."
        ),
        "dossier_anchor": "k3_final_parameter_inventory.md sec. A.2 line 200",
    },
    {
        "quantity": "acrylamide DEGRADATION parameters as constants",
        "value": None,
        "unit": "various",
        "where_it_lives": "Knol 2005 + Knol 2009 + Knol 2010",
        "verdict": "REFUSE AS CONSTANTS -- UNIDENTIFIABLE, three times over",
        "why": (
            "Knol 2005: fitted but 'the model was not restrained by "
            "experimental data for the products formed in the degradation "
            "reaction'. Knol 2009: EVERY degradation parameter has SD >= "
            "estimate (tc2 = 0.60 +/- 52, k2 = 3.5 +/- 8, tau = 22 +/- 15). "
            "Knol 2010: the Ea went negative and was deleted. THIS IS WHY "
            "k_acr_dp IS FITTED HERE rather than transcribed -- but note the "
            "distinction the inventory draws (sec. C.6): the degradation is "
            "unidentifiable in ITS OWN papers' multiresponse fits, not "
            "mis-valued. The three labs' RATES agree to 1.26x."
        ),
        "dossier_anchor": "k3_final_parameter_inventory.md sec. C.6",
    },
    {
        "quantity": "Quan 2020, all 88 rate constants",
        "value": None,
        "unit": "NONE PRINTED",
        "where_it_lives": "Quan et al. 2020 Table 1",
        "verdict": "REFUSE -- four disqualifying defects",
        "why": (
            "No units on any k; no rate laws (the SI figure defining them is "
            "absent, so every reaction ORDER is unknown); mislabelled "
            "continuation headers; no uncertainties despite a footnote claiming "
            "them. The authors themselves write that 'specific activation "
            "energies cannot be estimated'. Declared 'neither' by the "
            "declaration."
        ),
        "dossier_anchor": "k3_final_parameter_inventory.md sec. C.5",
    },
    {
        "quantity": "_acrylamide_ph_factor",
        "value": None,
        "unit": "dimensionless",
        "where_it_lives": "the shipped FAST lane",
        "verdict": "REFUSE -- no pH dependency is licensed for this module",
        "why": (
            "Worth ~2000x at pH 5.5 and responsible for the whole of the 480x "
            "two-lane acrylamide contradiction. Nothing in the Module 3 fit "
            "corpus varies pH, so there is no data from which such a factor "
            "could be estimated, and there is no place in this registry to put "
            "one: it holds rate constants on balanced reactions and nothing "
            "else."
        ),
        "dossier_anchor": "k3_final_parameter_inventory.md sec. A.2 line 206 (Z0 #13)",
    },
)

#: Numeric values that must never appear as a populated barrier in this module.
_FABRICATED_BARRIERS_KJ_MOL: Tuple[float, ...] = (129.0,)


def assert_no_fabricated_barrier(
    registry: Mapping[str, AcrylamideParameter] = MEASURED_ACRYLAMIDE,
) -> None:
    """Refuse the import if the unsourceable 129 kJ/mol barrier ever appears."""
    for key, parameter in registry.items():
        ea = getattr(parameter, "ea_kj_mol", None)
        if ea is None:
            continue
        for banned in _FABRICATED_BARRIERS_KJ_MOL:
            if abs(float(ea) - banned) < 1e-9:
                raise ValueError(
                    f"{key}: activation energy {banned} kJ/mol is the shipped "
                    f"acrylamide barrier that appears in NO Knol paper "
                    f"(k3_final_parameter_inventory.md sec. A.2 line 205). It "
                    f"may not enter the kinetic core."
                )


def assert_no_dft_acrylamide(
    registry: Mapping[str, AcrylamideParameter] = MEASURED_ACRYLAMIDE,
) -> None:
    """Owner policy, pinned: no DFT-derived barrier may enter this module."""
    banned = ("dft", "b3lyp", "wb97", "m06", "cbs-qb3", "g4", "computed", "calculated")
    for key, parameter in registry.items():
        haystack = " ".join([
            parameter.source_anchor, parameter.dossier_anchor,
            parameter.note, parameter.evidence_class,
        ]).lower()
        for token in banned:
            if token in haystack:
                raise ValueError(
                    f"{key}: parameter provenance mentions {token!r}. "
                    f"DFT-derived barriers are refused by owner policy."
                )


assert_no_dft_acrylamide()
assert_no_fabricated_barrier()


# ===========================================================================
# CROSS-LAB CONFLICTS -- reported, never averaged into a single number
# ===========================================================================

CROSS_LAB_CONFLICTS: Tuple[Dict[str, Any], ...] = (
    {
        "quantity": "acrylamide ELIMINATION rate constant at 160 C",
        "values": {
            "Claeys 2005 T2 control, aqueous pH 6": 111.1e-3,
            "De Vleeschouwer 2009 I glucose, a_w 0.92": 0.10,
            "Knol 2005 T1 k6, aqueous": 88.1e-3,
        },
        "unit": "1/min",
        "spread": "1.26x across three labs, two water activities and two decades",
        "how_this_module_handles_it": (
            "All three are FIT rows. ONE (k, Ea) pair is fitted to all three "
            "rather than one source being adopted; the row residuals then show "
            "the disagreement."
        ),
        "verdict": (
            "THE STRONGEST CROSS-LAB AGREEMENT IN THE MODULE, and it is on the "
            "very constant the inventory calls unidentifiable (sec. C.6). Both "
            "statements are true: the RATE is reproducible, the DEGRADATION "
            "MECHANISM behind it is not constrained by any of the three."
        ),
    },
    {
        "quantity": "acrylamide ELIMINATION activation energy",
        "values": {
            "Claeys 2005 T2 control": 167.21,
            "De Vleeschouwer 2009 I glucose": 113.2,
            "Knol 2005 T1 k6": 85.1,
        },
        "unit": "kJ/mol",
        "spread": "2.0x, intervals do not all overlap (167.21 +/- 4.30 vs 85.1 +/- 14)",
        "how_this_module_handles_it": (
            "Fitted to a single compromise value whose residual against each "
            "source is printed. NOT averaged silently, and NOT resolved."
        ),
        "verdict": (
            "A STANDING OWNER QUESTION. The same three labs agree to 1.26x on "
            "the RATE at 160 C while spanning 2x on the BARRIER, which is the "
            "signature of a reference temperature sitting in the middle of "
            "every data set: k(T_ref) is well determined and the slope is not."
        ),
    },
    {
        "quantity": "acrylamide FORMATION activation energy",
        "values": {
            "Claeys 2005 T2 control, aqueous": 168.25,
            "De Vleeschouwer 2009 I glucose, a_w 0.92": 159.2,
            "Knol 2005 T1 k4, aqueous": 94.4,
        },
        "unit": "kJ/mol",
        "spread": "74 kJ/mol between Claeys and Knol; intervals (168.25 +/- 3.80 "
                  "and 94.4 +/- 11) are nowhere near overlapping",
        "how_this_module_handles_it": (
            "Both are scored FIT rows and both residuals are printed. The model "
            "cannot satisfy both and no attempt is made to; which one it lands "
            "near is a reported outcome."
        ),
        "verdict": (
            "Claeys and De Vleeschouwer (same lab, different a_w) agree to "
            "9 kJ/mol; Knol is the outlier by 65-74. This is the largest "
            "unresolved numerical conflict in Module 3."
        ),
    },
    {
        "quantity": "condensation barrier, Glc + Asn -> Schiff base",
        "values": {
            "Knol 2005 T1 k1": 57.6,
            "Martins 2005 step 1 (glycine, B1's operative value)": 96.8,
            "De Vleeschouwer 2009 I k_INTg (this module's operative value)": 117.5,
        },
        "unit": "kJ/mol",
        "spread": "60 kJ/mol; no two intervals overlap",
        "how_this_module_handles_it": (
            "De Vleeschouwer's is operative because it is the only one measured "
            "as a genuinely BIMOLECULAR constant, which is what a mass-action "
            "network needs. The other two are carried here and in B1's "
            "CROSS_LAB_COMPARATORS."
        ),
        "verdict": "Owner question, already flagged by the inventory as Z2 #17.",
    },
)


# ===========================================================================
# DELIBERATE UNDER-FITS -- stated BEFORE the fit, not after
# ===========================================================================
# A model that cannot reproduce a FIT row it was given is normally a defect. In
# one case here it is a DECISION, and recording it in the registry is what
# stops it being presented later as a discovery or quietly patched with a term.

DELIBERATE_UNDERFITS: Tuple[Dict[str, Any], ...] = (
    {
        "row": "Claeys 2005 T2, the GLUTAMINE competitor system's formation constant",
        "what_the_model_will_do": (
            "UNDER-PREDICT it by roughly the measured promotion factor. Every "
            "competition channel in this module can only REMOVE precursor or "
            "REMOVE acrylamide, so the model's glutamine system can never form "
            "MORE acrylamide than its control."
        ),
        "why_no_promotion_TERM_IS_ADDED": (
            "Because the shape belongs to a HOLD-OUT. Inventory sec. B5.5: "
            "glutamine's acrylamide promotion GROWS with temperature in the "
            "liquid system and SHRINKS with temperature at a_w 0.92 -- same "
            "lab, same amino acid, and neither paper remarks on it. The liquid "
            "half is Claeys (FIT); the a_w 0.92 half is De Vleeschouwer II's "
            "glutamine system (HOLD-OUT). A promotion term fitted to the FIT "
            "half would be a term built toward a hold-out's shape, and the "
            "build brief for this wave forbids exactly that."
        ),
        "what_would_license_one": (
            "A measured mechanism. Glutamine deamidates to glutamic acid and "
            "the corpus has a rate for it -- in the HOLD-OUT column. Until a "
            "FIT-side measurement of a glutamine promotion route exists, the "
            "honest model is one that misses the row."
        ),
    },
)


# ===========================================================================
# THE HOLD-OUT DISCLOSURE
# ===========================================================================
# Recorded in the runtime registry, not only in the report, because an exposure
# that lives only in a generated file is one `rm` away from being invisible.

HOLDOUT_EXPOSURE_DISCLOSURE: Dict[str, Any] = {
    "what_was_seen": (
        "The build brief directed this wave to read "
        "k3_final_parameter_inventory.md sec. A.2 and "
        "k1_kinetic_parameters.md sec. 2, and BOTH print the declared HOLD-OUT "
        "columns beside the FIT ones: De Vleeschouwer 2009 Part I's FRUCTOSE "
        "and SUCROSE columns (sec. 2b prints all three sugars in one table, and "
        "the inventory prints the sucrose hydrolysis and isomerisation "
        "constants in its own sec. A.2), and De Vleeschouwer 2009 Part II's "
        "GLUTAMINE column (sec. 2c prints Gln and Cys side by side). Inventory "
        "sec. B5.5 additionally prints the glutamine promotion percentages that "
        "are the shape of the Part II hold-out."
    ),
    "what_was_not_seen": (
        "Knol 2010 Table 2's seven steps are NOT transcribed anywhere in the "
        "dossiers -- only the three organic-acid/isomerisation barriers, which "
        "orchestrator decision 2 assigns to Module 4 as FIT. The frozen "
        "external_validation bundles were never opened."
    ),
    "what_was_done_about_it": (
        "Nothing held out entered a parameter, a bound, an initialisation or a "
        "fit row. tests/unit/test_kinetic_core_b3.py enforces it mechanically "
        "with a literal-grep firewall over the executable code of every runtime "
        "file and the fit script, in the pattern Wave B2 established."
    ),
    "why_it_is_recorded": (
        "A hold-out whose exposure is undisclosed is worth less than one whose "
        "exposure is stated."
    ),
}


# ===========================================================================
# Metadata
# ===========================================================================

def acrylamide_registry_metadata(
    parameters: Mapping[str, AcrylamideParameter]
) -> Dict[str, Any]:
    """Full runtime metadata block for the acrylamide parameter set."""
    return {
        "acrylamide_reference_temperature_K": T_REF_A_K,
        "ph_term_present": False,
        "aw_term_present": False,
        "ph_term_absent_because": (
            "nothing in the Module 3 fit corpus varies pH: Claeys is 0.05 M "
            "citrate at pH 6, the extrusion benchmark is pH 6, De Vleeschouwer "
            "is a freeze-dried powder with no pH, and Knol is unbuffered "
            "aqueous. The shipped _acrylamide_ph_factor (~2000x at pH 5.5) is "
            "the whole of the 480x two-lane contradiction and is REFUSED."
        ),
        "aw_term_absent_because": (
            "the fit panel spans a_w 0.35-1.0 but NOTHING in it measures one "
            "step at two water activities, so an a_w exponent would be fitted "
            "to a between-system difference and would absorb every other "
            "difference between those systems as well. Instead every parameter "
            "carries aw_of_measurement and every run reports the gap."
        ),
        "competition_is_a_multiplier": False,
        "competition_mechanism": (
            "two named mass-action channels per competitor: consumption of the "
            "SHARED GLUCOSE pool (measured for cysteine as k_INT2), and "
            "Michael-acceptor scavenging of acrylamide (measured for cysteine "
            "as k_E2). There is no per-amino-acid yield multiplier in the "
            "registry and no place to put one."
        ),
        "acrylamide_has_elimination": True,
        "acrylamide_elimination_channels": [
            "k_acr_dp (first order, fitted to three independent measurements)",
            "k_acr_cys (bimolecular, MEASURED, De Vleeschouwer II k_E2)",
            "k_acr_gln / k_acr_lys / k_acr_ala (bimolecular, fitted, "
            "Ea held at the measured Ea_E2)",
        ],
        "dft_free": True,
        "fabricated_barrier_129_present": False,
        "schiff_amadori_split": dict(ASN_SCHIFF_AMADORI_SPLIT),
        "refused_parameters": [dict(row) for row in REFUSED_PARAMETERS],
        "cross_lab_conflicts": [dict(row) for row in CROSS_LAB_CONFLICTS],
        "deliberate_underfits": [dict(row) for row in DELIBERATE_UNDERFITS],
        "holdout_exposure_disclosure": dict(HOLDOUT_EXPOSURE_DISCLOSURE),
        "parameters": {k: p.as_metadata() for k, p in sorted(parameters.items())},
    }
