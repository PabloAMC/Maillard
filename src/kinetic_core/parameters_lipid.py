"""
src/kinetic_core/parameters_lipid.py

THE LIPID-OXIDATION PARAMETER REGISTRY (Build Wave B6, 2026-08-29).
=============================================================================

THE ONE THING TO READ FIRST
---------------------------
**The branch DISTRIBUTION is measured. The absolute RATE is not.**

That asymmetry is the whole design of this module, and it is a declared gap,
not an oversight:

  ``data/lit/extraction_dossiers/research_round3_channels.md`` §F.3 --
  "No measured rate constant exists for linoleate hydroperoxide -> hexanal at
  cooking temperatures in aqueous or protein systems. The only measured
  ``k_hexanal`` is at 25 C and was hand-fitted without standard errors.
  Extrapolating it to the repo's window using its authors' own Q10 = 2-3 spans
  11.5 decades of 10 C steps (~3e3-8e5x) and is **the repo's assumption, not
  the authors'**."

  ``k3`` §C.9 on Frankel 1989 -- "NOT a yield source ... no absolute yield and
  no Ea exist in it (one temperature, 180 C)."

So this registry carries:

  * the branch topology and the FIT column (Frankel 1989's three zero-additive
    columns) as MEASURED composition data;
  * ``K_LOOH_DECOMP_ANCHOR`` -- the Schroen 25 C constant, as a BOUNDED INPUT
    with every disqualifying flag its own authors printed;
  * ``Q10_ASSUMPTION`` -- a user-visible parameter with a declared default and a
    band, NOT a number baked into a rate constant. Nothing in this file
    multiplies the anchor by a temperature factor; that happens at call time,
    with the assumption travelling alongside the number.

THREE STANDING POLICIES, INHERITED
----------------------------------
1. NO DFT. ``assert_no_dft_lipid()`` is called at import.
2. NO INVENTED BRANCH FRACTIONS. ``PROHIBITED_DERIVATIONS`` names each one that
   was asked for and refused, with the reason.
3. Every declared assumption carries a BAND, and the band is propagated into
   the prediction interval rather than hidden inside a point value.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, Mapping, Optional, Tuple

# ---------------------------------------------------------------------------
# 0. Evidence classes -- B1's four, plus the two this module actually needs
# ---------------------------------------------------------------------------

LIPID_EVIDENCE_CLASSES: Tuple[str, ...] = (
    "measured_composition",       # a distribution printed in a paper
    "hand_fitted_no_uncertainty",  # a value the AUTHORS state was eyeballed
    "declared_assumption",        # not measured, not fitted: declared + bounded
    "structural_constant",        # a stoichiometry or a topology gate
    "derived_from_fit_data",      # estimated here, on FIT rows only
)


@dataclass(frozen=True)
class LipidParameter:
    """One lipid constant, with everything needed to defend or refuse it."""

    key: str
    what: str
    value: Optional[float]
    unit: str
    evidence_class: str
    source_anchor: str
    dossier_anchor: str
    conditions: str
    temperature_of_measurement_c: Optional[float]
    ph_of_measurement: Optional[float]
    #: "licensed" | "not_licensed" | "licensed_as_an_anchor_only"
    rate_transfer: str
    lo: Optional[float] = None
    hi: Optional[float] = None
    flags: Tuple[str, ...] = ()
    note: str = ""

    def __post_init__(self) -> None:
        if self.evidence_class not in LIPID_EVIDENCE_CLASSES:
            raise ValueError(
                f"{self.key}: evidence_class {self.evidence_class!r} not in "
                f"{LIPID_EVIDENCE_CLASSES}. Computational barriers are refused "
                f"by owner policy."
            )
        if self.lo is not None and self.hi is not None and self.lo > self.hi:
            raise ValueError(f"{self.key}: lo > hi")

    def as_metadata(self) -> Dict[str, Any]:
        return {
            "key": self.key,
            "what": self.what,
            "value": self.value,
            "unit": self.unit,
            "band": (None if self.lo is None else [self.lo, self.hi]),
            "evidence_class": self.evidence_class,
            "source": self.source_anchor,
            "dossier": self.dossier_anchor,
            "conditions": self.conditions,
            "T_of_measurement_C": self.temperature_of_measurement_c,
            "pH_of_measurement": self.ph_of_measurement,
            "rate_transfer": self.rate_transfer,
            "flags": list(self.flags),
            "note": self.note,
        }


# ===========================================================================
# 1. THE FIT COLUMN -- Frankel 1989's THREE zero-additive columns, and nothing
#    else from that paper
# ===========================================================================
# Relative percent of the sum of the SIX major volatiles, exactly as printed.
# Table 1 column "0"; Table 2 cis,trans-13 column "0"; Table 2 trans,trans
# 9+13 column "0". The tocopherol and 1,4-cyclohexadiene columns are HOLD-OUT
# and appear NOWHERE in this package -- a literal-grep firewall test asserts it.
#
# Provenance: data/articles/frankel1989.pdf, E.N. Frankel & H.W. Gardner,
# Lipids 24:603-608 (1989), thermolysis in a GC injector port at 180 C,
# relative to methyl hexanoate internal standard, RSD of duplicates 3.9-4.8 %.

FRANKEL_ZERO_ADDITIVE: Mapping[str, Mapping[str, float]] = {
    # Table 1, column "0": the AUTOXIDATION mixture (cis,trans AND trans,trans
    # 9- and 13-hydroperoxides). This is the composition a real, autoxidised
    # food lipid most nearly resembles.
    "mixed_ct_tt_9_13": {
        "PENTANE": 16.0,
        "HEXANAL": 11.0,
        "ME_OCTANOATE": 17.0,
        "DECADIENAL": 23.0,
        "ME_9_OXONONANOATE": 13.0,
        "ME_13_OXO_TRIDECADIENOATE": 20.0,
    },
    # Table 2, left block, column "0": PURE methyl cis,trans-13-hydroperoxide,
    # from soy lipoxygenase. The 9-derived products are NOT zero here, which
    # is what identifies the preparation's isomeric purity.
    "pure_ct_13": {
        "PENTANE": 21.0,
        "HEXANAL": 20.0,
        "ME_OCTANOATE": 5.0,
        "DECADIENAL": 3.5,
        "ME_9_OXONONANOATE": 4.3,
        "ME_13_OXO_TRIDECADIENOATE": 46.0,
    },
    # Table 2, right block, column "0": the trans,trans 9+13 pair.
    "tt_9_13": {
        "PENTANE": 1.5,
        "HEXANAL": 13.0,
        "ME_OCTANOATE": 13.0,
        "DECADIENAL": 30.0,
        "ME_9_OXONONANOATE": 26.0,
        "ME_13_OXO_TRIDECADIENOATE": 16.0,
    },
}

#: Which geometry each FIT system's hydroperoxide pool carries. The mixed
#: system carries BOTH, with the split a fitted parameter.
FRANKEL_SYSTEM_GEOMETRY: Mapping[str, str] = {
    "mixed_ct_tt_9_13": "mixed",
    "pure_ct_13": "ct",
    "tt_9_13": "tt",
}

#: Total peak area relative to the methyl hexanoate internal standard, at zero
#: additive, for each FIT system. NOT fitted (one temperature, no absolute
#: yield, and the three preparations were injected at different loadings); it
#: is carried so the fit report can state what it declined to use.
FRANKEL_ZERO_ADDITIVE_TOTAL_AREA: Mapping[str, float] = {
    "mixed_ct_tt_9_13": 16.0,
    "pure_ct_13": 2.2,
    "tt_9_13": 2.3,
}

FRANKEL_ANCHOR = (
    "Frankel & Gardner (1989) Lipids 24:603-608, Table 1 col. '0' and Table 2 "
    "cols. '0' (both blocks); data/articles/frankel1989.pdf, read directly by "
    "Wave B6"
)
FRANKEL_DOSSIER = (
    "data/lit/extraction_dossiers/k3_final_parameter_inventory.md sec. A.5 and "
    "sec. C.9; docs/reference/FIT_HOLDOUT_DECLARATION.md D.6 Module 5"
)

#: The refutation this FIT column exists to deliver, stated as data.
SHIPPED_VALUES_REFUTED: Mapping[str, str] = {
    "hexanal 0.37": (
        "data/lit/lipid_oxidation_calibration.json branching_ratios both paths, "
        "used at src/lipid_oxidation.py:446 and :512. Frankel's measured hexanal "
        "share spans 11-20 % across the three zero-additive columns and 11-26 % "
        "across all 26 columns: 0.37 sits ABOVE the paper's entire measured "
        "range. B6 does not edit that file (it belongs to the FAST lane); it "
        "declares its own registry and reports the contradiction."
    ),
    "nonanal 0.15": (
        "same file, hexanal_path. Refuted STRUCTURALLY, not numerically: nonanal "
        "is the C9 fragment of the OLEATE double bond and cannot come from a "
        "linoleate hydroperoxide at all. See species_lipid.NONANAL."
    ),
}


# ===========================================================================
# 2. THE RATE -- a bounded INPUT, never a fitted barrier
# ===========================================================================

#: The ONLY measured decomposition constant for a lipid hydroperoxide pool in a
#: protein-containing aqueous system anywhere in the corpus. Every one of its
#: flags is the AUTHORS' own statement, quoted in the dossier.
K_LOOH_DECOMP_ANCHOR = LipidParameter(
    key="k_LOOH_decomp",
    what="first-order decomposition of the lipid hydroperoxide pool to ALL "
         "secondary oxidation products (the authors' k4, an explicit LUMP)",
    value=6.0e-3,
    unit="1/h",
    evidence_class="hand_fitted_no_uncertainty",
    source_anchor="Schroen & Berton-Carabin (2022) Food Res. Int. 160:111621, "
                  "Table 1, k4,CD -- identical across all five emulsifiers "
                  "(beta-lactoglobulin, BSA, beta-casein, Tween 20, Tween 80)",
    dossier_anchor="research_round3_channels.md sec. D.1; "
                   "data/articles/schroen2022_fulltext.txt",
    conditions="rapeseed-oil O/W emulsion, droplets 1.4-1.8 um, pH 6.7 buffer, "
               "rotative agitation, controlled oxygen-to-lipid ratio",
    temperature_of_measurement_c=25.0,
    ph_of_measurement=6.7,
    rate_transfer="licensed_as_an_anchor_only",
    lo=None, hi=None,
    flags=(
        "hand_fitted_by_visual_agreement",
        "no_standard_error_anywhere_in_the_source",
        "lumped_over_all_secondary_products",
        "storage_temperature_25C_only",
        "temperature_dependence_UNMEASURED",
    ),
    note="Verbatim from the source: 'The parameter values (k1-k5) were "
         "determined based on visual agreement of fit.' and 'k4 is a lumped "
         "reaction rate, ultimately leading to formation of secondary "
         "oxidation products'. It is a FLOOR ANCHOR for the module's rate, "
         "never 'the hexanal barrier'.",
)

#: The companion hexanal-specific constant, carried for the RATIO it implies
#: rather than for its level. 6e-5 / 6e-3 = 1.2 %, which is the number the
#: dossier calls the robust export.
K_HEXANAL_SCHROEN = LipidParameter(
    key="k_hexanal_schroen",
    what="first-order hexanal formation from the lipid hydroperoxide pool",
    value=6.0e-5,
    unit="1/h",
    evidence_class="hand_fitted_no_uncertainty",
    source_anchor="Schroen & Berton-Carabin (2022) Table 1, k_hexanal -- "
                  "identical across all five emulsifiers",
    dossier_anchor="research_round3_channels.md sec. D.1",
    conditions="as k_LOOH_decomp",
    temperature_of_measurement_c=25.0,
    ph_of_measurement=6.7,
    rate_transfer="not_licensed",
    flags=("hand_fitted_by_visual_agreement", "no_standard_error",
           "25C_storage_only"),
    note="Its LEVEL does not transfer. Its RATIO to k4 does: the authors state "
         "'The corresponding constants in Table 1 can be interpreted as a "
         "percentage of the total aldehydes formed related to k4 (7 % and "
         "1.2 %, for propanal and hexanal, respectively).' That 1.2 % is a "
         "25 C, whole-oil (C18:1 618.5 / C18:2 191.7 / C18:3 92.2 mg/g) "
         "figure and is NOT the same quantity as Frankel's within-slate share, "
         "which is a fraction of six measured peaks from a PURE linoleate "
         "hydroperoxide feed. The two are reported side by side in the fit "
         "report and never averaged.",
)

#: The propanal : hexanal = 5.8 : 1 ratio, the dossier's other robust export.
PROPANAL_HEXANAL_RATIO = LipidParameter(
    key="propanal_over_hexanal",
    what="ratio of the two named secondary-product constants",
    value=7.0 / 1.2,
    unit="dimensionless",
    evidence_class="measured_composition",
    source_anchor="Schroen & Berton-Carabin (2022), 7 % and 1.2 % of k4",
    dossier_anchor="research_round3_channels.md sec. D.1 export 2",
    conditions="rapeseed oil, so the propanal reflects that oil's 92.2 mg/g "
               "alpha-linolenate. It is a PROPERTY OF THAT OIL, not a "
               "transferable branch ratio.",
    temperature_of_measurement_c=25.0,
    ph_of_measurement=6.7,
    rate_transfer="not_licensed",
    flags=("oil_composition_specific",),
    note="Carried as a CROSS-CHECK on any future linolenate channel. B6 does "
         "NOT implement a linolenate channel: Frankel 1989 fed linoleate only, "
         "so no propanal branch fraction is measured in the FIT column, and "
         "importing this ratio would be importing rapeseed oil's fatty-acid "
         "profile as if it were chemistry.",
)


@dataclass(frozen=True)
class Q10Assumption:
    """
    THE MODULE'S CENTRAL ASSUMPTION, made a first-class object so that it
    cannot be a number hiding inside a rate constant.

    ``default`` is the geometric mean of the authors' own stated band. It is a
    DECLARED ASSUMPTION with class ``declared_assumption`` -- not a measurement,
    not a fit. Every prediction that uses it carries ``warning`` verbatim.
    """

    default: float = math.sqrt(2.0 * 3.0)   # 2.449...
    lo: float = 2.0
    hi: float = 3.0
    reference_temperature_c: float = 25.0
    source: str = (
        "Schroen & Berton-Carabin (2022) sec. 3.2.2 and sec. 4.2, twice, "
        "verbatim: 'In general the reaction rate constant would need to be "
        "decreased by a factor of 2-3 for every 10 C temperature difference' "
        "and 'parameters will need to be adjusted according to their "
        "temperature-sensitivity, which typically is in the order of a factor "
        "of 2-3 per 10 C temperature difference'."
    )
    licensed_span_c: Tuple[float, float] = (15.0, 40.0)
    warning: str = (
        "EXTRAPOLATION WARNING -- THE LIPID LANE'S RATE IS AN ASSUMPTION, NOT A "
        "MEASUREMENT. The hydroperoxide decomposition constant is anchored at "
        "25 C (Schroen & Berton-Carabin 2022, k4 = 6e-3 /h), hand-fitted by "
        "visual agreement, with NO standard error anywhere in the source, in a "
        "rapeseed O/W emulsion at pH 6.7, and it is an explicit LUMP over all "
        "secondary products. Its TEMPERATURE DEPENDENCE IS MEASURED NOWHERE "
        "(declared gap: research_round3_channels.md sec. F.3, re-affirming k3 "
        "sec. C.9). The Q10 applied here is the authors' own stated 2-3, but "
        "they licensed ADJUSTMENT, not an extrapolation across ~11.5 decades of "
        "10 C steps (a factor of 3e3-8e5). The BRANCH DISTRIBUTION is measured; "
        "the ABSOLUTE RATE is not. Ratios between formulations at a common rate "
        "assumption are first-class; absolute ppb inherits this band."
    )

    def factor(self, temperature_c: float, q10: Optional[float] = None) -> float:
        """The multiplicative factor on the 25 C anchor. No default is hidden."""
        value = self.default if q10 is None else float(q10)
        return float(value) ** ((float(temperature_c) - self.reference_temperature_c) / 10.0)

    def decades_of_extrapolation(self, temperature_c: float) -> float:
        """How many 10 C steps beyond the anchor. The honesty metric."""
        return (float(temperature_c) - self.reference_temperature_c) / 10.0

    def as_metadata(self) -> Dict[str, Any]:
        return {
            "key": "Q10",
            "default": self.default,
            "band": [self.lo, self.hi],
            "evidence_class": "declared_assumption",
            "source": self.source,
            "licensed_span_C": list(self.licensed_span_c),
            "warning": self.warning,
        }


Q10_ASSUMPTION = Q10Assumption()


def k_looh_decomp_per_min(
    temperature_c: float, q10: Optional[float] = None
) -> float:
    """
    The operative decomposition constant, in 1/min.

    THE ONLY PLACE the anchor and the assumption meet. Callers that want the
    band call this three times, with ``q10`` = lo / default / hi, which is
    exactly what the interval reporting does -- the band is propagated by
    RE-INTEGRATION, not by adding a nominal width to the answer.
    """
    per_hour = float(K_LOOH_DECOMP_ANCHOR.value) * Q10_ASSUMPTION.factor(
        temperature_c, q10
    )
    return per_hour / 60.0


# ===========================================================================
# 3. THE OXIDATION-STATE PROXY AND THE LIPID CARRIERS
# ===========================================================================
# The absolute scale of every lipid prediction is set by the SIZE OF THE
# HYDROPEROXIDE POOL AT THE START OF THE THERMAL STEP, and no source in the fit
# corpus measures it for any of the exam's matrices. It is therefore a declared,
# BOUNDED input with a standard food-analysis unit (peroxide value), so that a
# user who has measured it can supply it and collapse the largest single band in
# the module.

#: Peroxide value is printed in milliequivalents of active oxygen per kg fat.
#: 1 meq active oxygen = 0.5 mmol hydroperoxide. Arithmetic, not a parameter.
MEQ_O2_PER_MMOL_LOOH = 2.0


@dataclass(frozen=True)
class LipidCarrier:
    """
    A food matrix, described as far as the corpus and general food composition
    allow, with every unmeasured quantity BANDED and flagged.
    """

    key: str
    display: str
    lipid_mass_fraction: float
    lipid_lo: float
    lipid_hi: float
    #: fatty-acid profile, % of total fatty acids
    linoleic_pct: float
    linolenic_pct: float
    oleic_pct: float
    #: oxidation-state proxy, meq active oxygen per kg fat
    peroxide_value_meq_per_kg: float
    pv_lo: float
    pv_hi: float
    #: bulk density used to turn a mass fraction into mmol/L; 1.0 kg/L is the
    #: bundles' own stated aqueous basis and is not a fitted quantity.
    density_kg_per_l: float
    evidence_class: str
    source_anchor: str
    note: str = ""


#: DECLARED ASSUMPTIONS, every one of them. The precedent is
#: parameters_matrix.MATRIX_LOADING['soy_paste_hong'], which carries the same
#: honesty label for the same reason: the papers do not report the composition.
LIPID_CARRIERS: Mapping[str, LipidCarrier] = {
    "pea_protein_isolate": LipidCarrier(
        "pea_protein_isolate", "pea protein isolate",
        0.025, 0.010, 0.060,
        50.0, 12.0, 22.0,
        10.0, 2.0, 40.0,
        1.0,
        "declared_assumption",
        "Fatty-acid profile and residual-lipid fraction carried forward "
        "VERBATIM from data/lit/lipid_oxidation_calibration.json "
        "matrix_lipid_profiles['pea_iso'], whose own provenance block says "
        "'literature-typical values already present in the codebase'. They are "
        "NOT corpus values and NOT in the FIT column.",
        "THE LARGEST ASSUMPTION IN THIS MODULE, together with Q10. The "
        "peroxide value is not reported by any source for any of the exam's "
        "matrices; the band [2, 40] meq/kg spans an unoxidised isolate to a "
        "badly rancid one. It is propagated into every interval, and a user "
        "who measures PV collapses it.",
    ),
    "soy_protein_isolate": LipidCarrier(
        "soy_protein_isolate", "soy protein isolate",
        0.020, 0.008, 0.050,
        53.0, 8.0, 23.0,
        10.0, 2.0, 40.0,
        1.0,
        "declared_assumption",
        "As pea, from matrix_lipid_profiles['soy_iso'].",
        "Same standing as the pea entry.",
    ),
    "frankel_pure_hydroperoxide": LipidCarrier(
        "frankel_pure_hydroperoxide",
        "fed methyl linoleate hydroperoxide (Frankel's reaction chromatography)",
        1.0, 1.0, 1.0,
        100.0, 0.0, 0.0,
        # A fed hydroperoxide is 100 % hydroperoxide by construction: PV is not
        # a proxy here, it is the definition. 1 mol LOOH/kg = 2000 meq/kg.
        2000.0, 2000.0, 2000.0,
        1.0,
        "structural_constant",
        "Frankel 1989 Methods: 13.1 mg methyl linoleate hydroperoxides injected "
        "from a 200 uL hexane solution, thermolysed in the injector port.",
        "Present so the FIT and the nonanal hold-out are run through the SAME "
        "code path as a real matrix, with no oleate at all.",
    ),
}


def looh_charge_mmol_per_l(
    carrier: LipidCarrier,
    *,
    lipid_fraction: Optional[float] = None,
    peroxide_value_meq_per_kg: Optional[float] = None,
) -> float:
    """
    The hydroperoxide pool at the start of the thermal step, in mmol/L.

    ``PV / 2`` mmol LOOH per kg fat, times kg fat per litre. Pure arithmetic.
    """
    fraction = carrier.lipid_mass_fraction if lipid_fraction is None else float(lipid_fraction)
    pv = (carrier.peroxide_value_meq_per_kg
          if peroxide_value_meq_per_kg is None else float(peroxide_value_meq_per_kg))
    kg_fat_per_l = fraction * carrier.density_kg_per_l
    return (pv / MEQ_O2_PER_MMOL_LOOH) * kg_fat_per_l


def oleate_fraction(carrier: LipidCarrier) -> float:
    """
    Share of the hydroperoxide pool that is OLEATE-derived.

    Weighted by fatty-acid abundance ONLY -- NOT by relative oxidisability.
    Oleate is roughly an order of magnitude less oxidisable than linoleate, so
    this OVERSTATES the oleate pool, which is the conservative direction: it
    makes the module refuse nonanal more often, never less.
    """
    total = carrier.linoleic_pct + carrier.linolenic_pct + carrier.oleic_pct
    if total <= 0.0:
        return 0.0
    return carrier.oleic_pct / total


# ===========================================================================
# 4. THE CONSUMPTION SIDE
# ===========================================================================
# Formation without consumption is the defect B1 was built to avoid, so every
# product here has a consumption term. What the corpus licenses for it is a
# CEILING and nothing else.

@dataclass(frozen=True)
class CovalentSinkCeiling:
    """
    The aldehyde -> lysine covalent channel, as a bounded ceiling.

    Amendment 6 ruling 1 brackets it two-sided across two labs, two decades,
    two protein families and two techniques. Ruling 2 computes that it supplies
    ~0.06 % of the hexanal matrix log-shift and does NOT close the gap; B4's
    ``covalent_channel_state`` therefore reports
    ``contribution_to_point_prediction: 0.0``, and this module agrees with B4
    rather than disagreeing with it silently.
    """

    k2_m_per_s_at_20c: float = 2.5e-5
    #: NOT the same bracket as parameters_matrix.COVALENT_CEILING's (37, 760).
    #: This is meynier2004 sec. 1(c)'s OVERLAP of two independent methods; that
    #: one is anantharamkrishnan2020b's single-method MS adduct-counting range.
    #: Wave B8 found the disagreement, verified both are as their sources print
    #: them, and CHANGED NEITHER: they answer different questions and picking
    #: one after reading a scorecard is not a repair. Reported, not resolved.
    ambient_half_life_days: Tuple[float, float] = (37.0, 74.0)
    activation_ea_threshold_kj_mol: float = 70.0
    #: WAVE B8 (Amendment 17 clause 6): MEASURED, 15-23 kJ/mol, central 20.0.
    #: Shepelev & Reineccius 2024 Fig. 3, 14-C on beta-lactoglobulin at 25/45/
    #: 65 C. The channel stays DISABLED -- but for the opposite reason from
    #: before: not because the deciding number is missing, but because it has
    #: been measured and it decides against.
    measured_ea_kj_mol: Optional[float] = 20.0
    measured_ea_range_kj_mol: Tuple[float, float] = (15.0, 23.0)
    enabled: bool = False
    why_disabled: str = (
        "INERT BY RULING (FIT_HOLDOUT_DECLARATION Amendment 6 ruling 2), AND "
        "NOW ALSO BY MEASUREMENT (Amendment 17 clause 6, Wave B8). Through B7 "
        "this field read: 'the channel's Ea on food proteins is UNMEASURED in "
        "every corpus source -- a NAMED WET-LAB GAP -- so evaluating it at "
        "process temperature would require inventing the one number that "
        "decides whether it matters at all.' WAVE K6b MEASURED THAT NUMBER: "
        "Ea = 15-23 kJ/mol against the >= 70 threshold, i.e. missed by "
        "3.5-4.7x, so at process temperature the channel removes 0.005%-0.21% "
        "of an aldehyde dose rather than the 1%-91% the assumption implied. "
        "Two further measurements disable it without needing the Arrhenius "
        "arithmetic at all: the sink's CAPACITY FALLS with temperature (day-28 "
        "binding 25.6 > 20.1 > 17.4 mg/g at 25/45/65 C, the 65 C series "
        "peaking at day 7 then declining 26%), and the EQUILIBRIUM UNWINDS on "
        "heating (Ea_rev 52 > Ea_fwd 44 => K_eq falls 3.0x from 25 to 180 C). "
        "It is an AMBIENT-STORAGE channel and no amount of heat makes it a "
        "process channel. IT REMAINS DISABLED AND ITS FLUX REMAINS EXACTLY "
        "0.0 -- the value did not move, the reason did. It is reported as a "
        "ceiling on every lipid prediction, and the lane-coupling guard in "
        "lipid.py asserts that it is still zero before permitting "
        "co-integration with a Maillard lane."
    )

    def rate_per_min_ea_zero_bound(
        self, aldehyde_mmol_l: float, lysine_mmol_l: float
    ) -> float:
        """
        The Ea = 0 LOWER bound on the sink, in mmol/(L min).

        Offered for sensitivity work only. It is a lower bound because a real
        Ea is positive; there is no upper bound without the missing Ea, which
        is precisely why the channel is off.
        """
        k_l_per_mmol_min = self.k2_m_per_s_at_20c * 60.0 / 1000.0
        return k_l_per_mmol_min * float(aldehyde_mmol_l) * float(lysine_mmol_l)


COVALENT_SINK = CovalentSinkCeiling()


#: The LOOH pool's own FORMATION during the thermal step: NOT MODELLED.
LOOH_FORMATION_GAP = (
    "NOT MODELLED, DECLARED. The hydroperoxide pool is an INPUT (an "
    "oxidation-state proxy), not a state variable that grows during heating. "
    "Schroen's initiation constants k1 (1.4e-5 to 6.5e-5 /h) and k5 are "
    "emulsifier-specific, 25 C, hand-fitted, and vary 4.6x between "
    "beta-lactoglobulin and beta-casein -- i.e. the ONE part of that paper's "
    "scheme the authors say is system-dependent. Extrapolating a "
    "system-dependent initiation term across 115 C would add a second "
    "unmeasured factor on top of the one this module already declares. The "
    "consequence is stated plainly: for a long, hot process the module "
    "UNDER-predicts, because it cannot make new hydroperoxide."
)

#: The other consumption channel Frankel's own introduction names, and the
#: reason it is off.
DECADIENAL_TO_HEXANAL_GAP = (
    "NOT MODELLED, DECLARED. Frankel 1989's introduction cites Schieberle & "
    "Grosch (1981) for 'hexanal may form ... after oxidative decomposition of "
    "2,4-decadienal'. If that route is live, hexanal and 2,4-decadienal are "
    "coupled and the fitted 13-OOH/9-OOH split is partly an artefact. No rate "
    "for it exists in the corpus, so it is a structurally enabled edge with no "
    "constant: reported as a named risk to the FIT, not silently assumed absent."
)


# ===========================================================================
# 5. WHAT WAS ASKED FOR AND REFUSED
# ===========================================================================

PROHIBITED_DERIVATIONS: Mapping[str, str] = {
    "Ea for LOOH decomposition from Frankel 1989": (
        "ONE TEMPERATURE (180 C). k3 sec. C.9, verbatim: 'NOT a yield source ... "
        "no absolute yield and no Ea exist in it'. Re-affirmed by "
        "research_round3_channels.md sec. F.3 after a further search round."
    ),
    "Q10 baked into a rate constant": (
        "The Q10 is the REPO'S assumption, not the authors'. It is exposed as "
        "Q10_ASSUMPTION with a band and a mandatory warning, and no stored "
        "constant is ever pre-multiplied by it."
    ),
    "oleate -> nonanal branch fraction": (
        "Measured nowhere in the fit corpus. Frankel 1989 fed linoleate only. "
        "The FAST lane's shipped 0.15 has no source. Requests for absolute "
        "nonanal in an oleate-bearing matrix are REFUSED."
    ),
    "linoleate -> 2-pentylfuran branch fraction": (
        "Not in Frankel's six-product slate and measured nowhere else in the "
        "corpus. The FAST lane's shipped 0.08 has no source."
    ),
    "aldehyde -> alcohol reduction (hexanal -> 1-hexanol)": (
        "No reduction step is measured anywhere in the corpus, and in a "
        "thermally processed extrudate the reductant pool is not even "
        "identified. Refused."
    ),
    "propanal branch from the linolenate channel": (
        "Schroen's 7 % propanal is a property of RAPESEED OIL's 92.2 mg/g "
        "alpha-linolenate, not a branch fraction of a linoleate hydroperoxide. "
        "Frankel's slate contains no propanal. Importing it would import an "
        "oil's composition as if it were chemistry."
    ),
    "a mass yield of hexanal per hydroperoxide": (
        "Frankel's shares are fractions of SIX MEASURED PEAKS, and two of the "
        "scissions' partners (2-nonenal, methyl 12-oxo-10-dodecenoate) are "
        "named in his introduction and quantified in none of his tables. The "
        "share is therefore a WITHIN-SLATE share, and this module carries the "
        "difference in LIPID_FRAG_C rather than pretending the slate closes."
    ),
}


#: Contradictions this wave found and is REPORTING, not resolving.
LIPID_SOURCE_CONTRADICTIONS: Mapping[str, str] = {
    "schroen_1.2pct_vs_frankel_11_20pct": (
        "Schroen's hexanal is 1.2 % of the secondary-product flux at 25 C in "
        "whole rapeseed oil; Frankel's is 11-20 % of a six-peak slate at 180 C "
        "from a PURE linoleate hydroperoxide. These differ by ~10-17x and they "
        "are NOT the same quantity: different denominator (all secondary "
        "products vs six peaks), different feed (whole oil, 62 % oleate, vs "
        "pure linoleate hydroperoxide), different temperature (25 vs 180 C) and "
        "different mechanism mix (Frankel's heterolytic Hock route is "
        "heat-promoted and cannot be active at 25 C). Reported side by side; "
        "never averaged, never reconciled by a fitted factor."
    ),
    "pentane_vs_me13oxo_pairing": (
        "Pentane and methyl 13-oxo-9,11-tridecadienoate are the two halves of "
        "ONE beta-scission and should be 1:1. Their measured ratio across the "
        "three zero-additive FIT columns spans 8.5x (0.80 / 0.46 / 0.094). "
        "Either the GC recovery of a C5 alkane and a C14 dienoate differ "
        "enormously, or the pairing is wrong. This module does NOT impose the "
        "pairing and does NOT absorb the discrepancy into a response factor; "
        "it fits free shares inside each isomer's simplex and reports the "
        "falsified pairing."
    ),
}


def assert_no_dft_lipid() -> None:
    """No barrier in this module is computational. Called at import."""
    for parameter in (K_LOOH_DECOMP_ANCHOR, K_HEXANAL_SCHROEN, PROPANAL_HEXANAL_RATIO):
        if parameter.evidence_class not in LIPID_EVIDENCE_CLASSES:
            raise AssertionError(f"{parameter.key}: illegal evidence class")
        if "dft" in parameter.source_anchor.lower():
            raise AssertionError(f"{parameter.key}: DFT source")


def assert_fit_column_is_zero_additive_only() -> None:
    """
    THE FIREWALL. The FIT column must be exactly three systems x six products,
    and every value must be a zero-additive value.

    A tocopherol column entering this table would show up as a fourth system or
    a seventh product; a substituted value would break the simplex sum. Both are
    checked here, at import, so the failure is loud.
    """
    from .species_lipid import FRANKEL_SLATE

    if set(FRANKEL_ZERO_ADDITIVE) != {"mixed_ct_tt_9_13", "pure_ct_13", "tt_9_13"}:
        raise AssertionError(
            "the FIT column must contain EXACTLY the three zero-additive "
            "systems. A fourth system is a hold-out column leaking in."
        )
    for system, shares in FRANKEL_ZERO_ADDITIVE.items():
        if tuple(shares) != FRANKEL_SLATE:
            raise AssertionError(
                f"{system}: the product slate must be Frankel's six, in order"
            )
        total = sum(shares.values())
        if not 99.0 <= total <= 101.0:
            raise AssertionError(
                f"{system}: relative percents sum to {total}, not ~100"
            )


def lipid_registry_metadata() -> Dict[str, Any]:
    """Everything a report needs to defend this module's numbers."""
    return {
        "module": "src/kinetic_core/parameters_lipid.py",
        "wave": "B6 -- the lipid-oxidation module",
        "fit_column": {
            "source": FRANKEL_ANCHOR,
            "dossier": FRANKEL_DOSSIER,
            "systems": list(FRANKEL_ZERO_ADDITIVE),
            "n_values": sum(len(v) for v in FRANKEL_ZERO_ADDITIVE.values()),
            "role": "FIT (declaration D.6 Module 5)",
        },
        "rate_is_an_assumption": True,
        "rate_anchor": K_LOOH_DECOMP_ANCHOR.as_metadata(),
        "hexanal_anchor": K_HEXANAL_SCHROEN.as_metadata(),
        "q10": Q10_ASSUMPTION.as_metadata(),
        "carriers": {
            k: {
                "lipid_fraction": [c.lipid_lo, c.lipid_mass_fraction, c.lipid_hi],
                "peroxide_value_meq_per_kg": [c.pv_lo, c.peroxide_value_meq_per_kg, c.pv_hi],
                "evidence_class": c.evidence_class,
                "source": c.source_anchor,
                "note": c.note,
            }
            for k, c in LIPID_CARRIERS.items()
        },
        "covalent_sink": {
            "enabled": COVALENT_SINK.enabled,
            "k2_M_per_s_at_20C_ceiling": COVALENT_SINK.k2_m_per_s_at_20c,
            "measured_ea_kj_mol": COVALENT_SINK.measured_ea_kj_mol,
            "measured_ea_range_kj_mol": list(COVALENT_SINK.measured_ea_range_kj_mol),
            "decision_threshold_ea_kj_mol": COVALENT_SINK.activation_ea_threshold_kj_mol,
            "ea_status": "MEASURED (Wave B8, Amendment 17 clause 6)",
            "why_disabled": COVALENT_SINK.why_disabled,
        },
        "declared_gaps": {
            "LOOH_formation": LOOH_FORMATION_GAP,
            "decadienal_to_hexanal": DECADIENAL_TO_HEXANAL_GAP,
        },
        "prohibited_derivations": dict(PROHIBITED_DERIVATIONS),
        "source_contradictions": dict(LIPID_SOURCE_CONTRADICTIONS),
        "shipped_values_refuted": dict(SHIPPED_VALUES_REFUTED),
        "no_dft": True,
    }


assert_no_dft_lipid()
assert_fit_column_is_zero_additive_only()


__all__ = [
    "COVALENT_SINK",
    "DECADIENAL_TO_HEXANAL_GAP",
    "FRANKEL_ANCHOR",
    "FRANKEL_DOSSIER",
    "FRANKEL_SYSTEM_GEOMETRY",
    "FRANKEL_ZERO_ADDITIVE",
    "FRANKEL_ZERO_ADDITIVE_TOTAL_AREA",
    "K_HEXANAL_SCHROEN",
    "K_LOOH_DECOMP_ANCHOR",
    "LIPID_CARRIERS",
    "LIPID_SOURCE_CONTRADICTIONS",
    "LOOH_FORMATION_GAP",
    "LipidCarrier",
    "LipidParameter",
    "PROHIBITED_DERIVATIONS",
    "PROPANAL_HEXANAL_RATIO",
    "Q10_ASSUMPTION",
    "Q10Assumption",
    "SHIPPED_VALUES_REFUTED",
    "assert_no_dft_lipid",
    "assert_fit_column_is_zero_additive_only",
    "k_looh_decomp_per_min",
    "lipid_registry_metadata",
    "looh_charge_mmol_per_l",
    "oleate_fraction",
]
