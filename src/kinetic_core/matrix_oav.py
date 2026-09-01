"""
src/kinetic_core/matrix_oav.py

THE MATRIX / OAV OUTPUT LAYER (Build Wave B4).
==============================================

The last module before the propagator cutover: it turns the network's predicted
headspace-relevant concentrations into something a plant-protein scientist can
act on.

THE DESIGN, AND WHY IT IS THIS AND NOT SOMETHING MORE CONFIDENT
---------------------------------------------------------------
**RATIOS LEAD.** The primary outputs are formulation-vs-formulation ratios and
rankings, each with a validity bound. The evidence forces this:

  * A general matrix correction factor is refuted. Cross-study ratios span
    2 000x with a 1-sigma band of 27-41x; on a SAME-PANEL set the band is still
    4.7x, and inversions occur (1 in 10 in Hong 2020, 3 in 21 in Leksrisompong
    2010). A uniform factor misplaces the two extreme compounds by 10x and 28x
    IN OPPOSITE DIRECTIONS -- worse than doing nothing, because it manufactures
    confident wrong answers where the repo currently has an admitted gap.
    (k2_matrix_and_thresholds.md sec. D.1, D.4.)
  * Partition and reversible binding are refuted as the mechanism THREE
    independent ways on matched samples: Hong's log P correlation is r = -0.05;
    Leksrisompong's partition model errs up to 30.6x and gets 2 of 3 signs
    wrong; Meynier's two-phase model fails in both directions on two structural
    isomers. (k4b sec. B; Amendment 6.)
  * Perception is LESS sensitive to matrix binding than headspace is, measured
    within-subject on the same mouthfuls (Baek: 1.5x n.s. in concentration
    against 3.0x P<0.001 in perception), so a headspace-calibrated free-fraction
    OVER-predicts the threshold shift and is still 2-3 orders short.

**ABSOLUTES ARE INTERVALS, ALWAYS.** Where an absolute ug/L is emitted it
carries the measured reliability band -- HS-SPME same-sample dispersion 10-23x,
plus +/-0.5 decades on the air/water partition constant -- as an explicit
interval. ``AbsoluteConcentration`` has no float coercion, so a bare point
cannot leave this module by accident.

**THREE NAMED SHIFT TERMS, THEN AN HONEST RESIDUAL.** Reversible binding (capped
at ~25 % of a log-shift), the alpha,beta-unsaturation penalty, and the covalent
ceiling (inert: it contributes zero and only ever REPORTS). Whatever an
observation shows beyond those three is reported as UNEXPLAINED RESIDUAL, in
decades, per compound. The output the layer is proud of looks like: "measured
shift 132x, reversible binding explains 2.6x, unsaturation 1.0x, covalent 1.0x
(inert), UNEXPLAINED RESIDUAL 51x".

**THRESHOLDS ARE INPUTS, NEVER PREDICTIONS.** Scoring a threshold would be a
category error (declaration D.6 Module 7). The table carries water values and
MATRIX-CONDITIONAL values only where a measurement exists, with an explicit
``no_measured_threshold_for_this_matrix`` state everywhere else. Nothing is
borrowed across matrices.

WHAT THIS MODULE DOES NOT TOUCH
-------------------------------
``src/matrix_correction.py`` and ``src/headspace.py`` are the OLD lane. They are
not imported and not modified. The ``protein_source_registry.json``-derived
protein differentiation they once carried (declared ``no_verifiable_source``;
file WITHDRAWN 2026-09-01) is not reproduced here in any form.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union

from .parameters_matrix import (
    ADDUCT_POSITIVE_CLASSES,
    CHAIN_LENGTH_SLOPE_PER_CH2,
    COMPOUND_STRUCTURE,
    COVALENT_CEILING,
    COVALENT_CEILING_RETIREMENT,
    HS_SPME_SAME_SAMPLE_DISPERSION,
    K_AW_UNCERTAINTY_DECADES,
    MATRIX_LOADING,
    NO_ADDUCT_GATE_CLASSES,
    NO_ADDUCT_GATE_COMPOUNDS,
    PH_ADDUCT_GATE_BELOW,
    PH_ADDUCT_GATE_UNCERTAIN_BELOW,
    REVERSIBLE_BINDING,
    REVERSIBLE_LOG_SHIFT_CEILING,
    ALPHA_BETA_UNSATURATION_OBSERVATIONS,
    MatrixParameter,
)

# ===========================================================================
# 1. THE THRESHOLD TABLE -- a LOOKUP, not a formula
# ===========================================================================

#: The sentinel returned when no measured threshold exists for a
#: (compound, matrix) pair. It is deliberately NOT a number and NOT None:
#: a caller that tries to divide by it gets a TypeError, not a silent zero.
class NoMeasuredThreshold:
    __slots__ = ("compound", "matrix", "reason")

    def __init__(self, compound: str, matrix: str, reason: str) -> None:
        self.compound = compound
        self.matrix = matrix
        self.reason = reason

    def __repr__(self) -> str:  # pragma: no cover - trivial
        return (f"NoMeasuredThreshold(compound={self.compound!r}, "
                f"matrix={self.matrix!r}, reason={self.reason!r})")

    def as_dict(self) -> Dict[str, str]:
        return {"state": "no_measured_threshold_for_this_matrix",
                "compound": self.compound, "matrix": self.matrix,
                "reason": self.reason}


@dataclass(frozen=True)
class ThresholdRecord:
    """One odour threshold, with the provenance the ratios actually need."""

    compound: str
    matrix: str
    value_ug_per_l: float
    temperature_c: Optional[float]
    #: '50%_forced_choice' | '75%_uncorrected' | 'BET_3AFC_ASTM_E679' |
    #: 'triangle_ASTM_E1432' | 'not_stated'
    criterion: str
    method: str
    source: str
    #: Model INPUT, never a scored prediction (declaration D.6 Module 7).
    role: str = "reference_input"
    thermal_step_after_dosing: Optional[bool] = None
    concentration_verified: Optional[bool] = None
    cross_study_cross_method: Optional[bool] = None
    aqueous_reference_source: Optional[str] = None
    dispersion_scale: str = "not_stated"
    provenance_flag: str = ""
    notes: str = ""


#: --- AQUEOUS ------------------------------------------------------------
#: Zhou 2023 SI Table S2. PROVENANCE UNKNOWN BY THE SOURCE'S OWN ADMISSION:
#: the SI prints no reference column and no footnote for any of them. Ingested
#: `basis_declared: true, provenance: unknown` exactly as K3 sec. A.7.1
#: requires, and carried with an uncited flag on every record.
_ZHOU_S2 = {
    "MFTD": 0.00032, "MFT": 0.005, "FFT": 0.006, "ACTZ": 10.0,
    "2_4_dimethylthiazole": 18.0, "thiazole": 38.0, "methylpyrazine": 60.0,
    "2_5_dimethylpyrazine": 80.0, "thiophene": 84.0,
    "2_6_dimethylpyrazine": 400.0, "FUR": 3000.0,
    "2_thiophenecarboxaldehyde": 5000.0, "furan": 4500.0,
    "2_acetylfuran": 10000.0, "pyrazine": 75000.0,
}

WATER_THRESHOLDS: Tuple[ThresholdRecord, ...] = tuple(
    ThresholdRecord(
        compound=key, matrix="water", value_ug_per_l=value, temperature_c=None,
        criterion="not_stated", method="not_stated",
        source="Zhou 2023 SI Table S2, via k3_final_parameter_inventory.md sec. A.7.1",
        cross_study_cross_method=True,
        provenance_flag="basis_declared_true__provenance_UNCITED",
        notes="The SI gives NO citation, no reference column and no footnote for "
              "any threshold in this table. Reference input only.",
    )
    for key, value in _ZHOU_S2.items()
) + (
    # Guadagni 1963/72, held SECOND-HAND through Vega 1994. Every ratio in
    # k2 sec. A.8 divides by these, and Vega mis-cites the 1972 venue.
    ThresholdRecord("pentanal", "water", 12.0, None, "not_stated", "classical dilution",
                    "Guadagni 1963 via vega1994 (reference only, D.6 Module 7)",
                    cross_study_cross_method=True,
                    provenance_flag="quoted_second_hand_by_vega1994"),
    ThresholdRecord("hexanal", "water", 4.5, None, "not_stated", "classical dilution",
                    "Guadagni 1972 via vega1994 (reference only)",
                    cross_study_cross_method=True,
                    provenance_flag="quoted_second_hand_by_vega1994"),
    ThresholdRecord("heptanal", "water", 3.0, None, "not_stated", "classical dilution",
                    "Guadagni 1972 via vega1994 (reference only)",
                    cross_study_cross_method=True,
                    provenance_flag="quoted_second_hand_by_vega1994"),
    ThresholdRecord("t_2_hexenal", "water", 3.0, None, "not_stated", "classical dilution",
                    "Guadagni 1963 via vega1994 (reference only)",
                    cross_study_cross_method=True,
                    provenance_flag="quoted_second_hand_by_vega1994"),
    ThresholdRecord("t_2_octenal", "water", 3.0, None, "not_stated", "classical dilution",
                    "Guadagni 1972 via vega1994 (reference only)",
                    cross_study_cross_method=True,
                    provenance_flag="quoted_second_hand_by_vega1994"),
    ThresholdRecord("tt_2_4_decadienal", "water", 0.07, None, "not_stated",
                    "classical dilution",
                    "Guadagni 1972 via vega1994 (reference only)",
                    cross_study_cross_method=True,
                    provenance_flag="quoted_second_hand_by_vega1994"),
    # Recovered by K2/K3 from the Xin papers' own OAV arithmetic.
    ThresholdRecord("hexanal", "water", 5.0, None, "not_stated", "not_stated",
                    "recovered from Xin 2026's own OAV arithmetic (k2 sec. A.6)",
                    cross_study_cross_method=True,
                    provenance_flag="derived_by_back_solving_a_published_OAV",
                    notes="A SECOND water value for hexanal, 1.11x from Guadagni's "
                          "4.5. Both are carried; the spread is reported, never averaged away."),
    ThresholdRecord("nonanal", "water", 1.1, None, "not_stated", "not_stated",
                    "recovered from Xin 2026's own OAV arithmetic (k2 sec. A.6)",
                    cross_study_cross_method=True,
                    provenance_flag="derived_by_back_solving_a_published_OAV"),
    ThresholdRecord("2_pentylfuran", "water", 5.800, None, "not_stated", "not_stated",
                    "recovered TWICE, from four independent numbers across Xin 2026 "
                    "and Xin 2026b (k3 sec. A.7.2)",
                    cross_study_cross_method=True,
                    provenance_flag="derived_but_doubly_corroborated"),
    ThresholdRecord("1_octen_3_ol", "water", 1.5, None, "not_stated", "not_stated",
                    "printed VERBATIM in Xin 2026b p. 9 (no citation given), and "
                    "independently back-solved to 1.500 from the same paper's OAVs",
                    cross_study_cross_method=True,
                    provenance_flag="printed_and_independently_reproduced"),
    ThresholdRecord("2_heptanone", "water", 140.0, None, "not_stated", "not_stated",
                    "back-solved from Xin 2026b's own OAV range (k3 sec. A.7.2)",
                    cross_study_cross_method=True,
                    provenance_flag="derived_by_back_solving_a_published_OAV"),
    ThresholdRecord("3_sulfanylhexan_1_ol", "water", 0.022, None, "50%_forced_choice",
                    "3-AFC, ascending, RETRONASAL (30 mL held in mouth 5 s)",
                    "Starkenmann 2008 p. 9577", cross_study_cross_method=False,
                    provenance_flag="ALREADY_AN_IN_MOUTH_NUMBER",
                    notes="A saliva correction applied on top of this would DOUBLE-COUNT."),
)

#: --- MATRIX-CONDITIONAL -------------------------------------------------
#: ONLY where a measurement exists. Vega 1994's 3 % gelatin ladder is the
#: cleanest matrix threshold set in the corpus -- lipid-free, dosed at 22 C
#: with NO thermal step after dosing -- and is declaration D.6 Module 7 FIT
#: (as lookup-table entries, which is exactly how it is used here).
_VEGA_GELATIN = {
    "pentanal": {4.0: 47.0, 22.0: 41.0, 37.0: 34.0, 60.0: 22.0},
    "hexanal": {4.0: 90.0, 22.0: 58.0, 37.0: 34.0, 60.0: 38.0},
    "heptanal": {4.0: 108.0, 22.0: 79.0, 37.0: 62.0, 60.0: 50.0},
    "t_2_hexenal": {4.0: 170.0, 22.0: 109.0, 37.0: 79.0, 60.0: 60.0},
    "t_2_octenal": {4.0: 140.0, 22.0: 109.0, 37.0: 105.0, 60.0: 81.0},
    "tt_2_4_decadienal": {4.0: 112.0, 22.0: 64.0, 37.0: 89.0, 60.0: 64.0},
}

MATRIX_THRESHOLDS: Tuple[ThresholdRecord, ...] = tuple(
    ThresholdRecord(
        compound=compound, matrix="gelatin_3pct", value_ug_per_l=value,
        temperature_c=temperature,
        criterion="75%_uncorrected",
        method="single sample vs a gelatin control, ascending series, LINEAR fit, "
               "16 panellists, 3 replicates",
        source="Vega & Brewer 1994 Table 1, via k2_matrix_and_thresholds.md sec. A.2",
        thermal_step_after_dosing=False, concentration_verified=False,
        cross_study_cross_method=False,
        aqueous_reference_source="Guadagni 1963/72 (second-hand, cross-method)",
        dispersion_scale="not_stated",
        notes="A 75 % UNCORRECTED criterion yields a systematically HIGHER "
              "threshold than a 50 % forced-choice one; sign known, magnitude not "
              "(~2x is a plausible size). The gelatin/water RATIO is therefore "
              "cross-method even though the ladder itself is one study.",
    )
    for compound, ladder in _VEGA_GELATIN.items()
    for temperature, value in ladder.items()
) + (
    # Paraffin oil, Guadagni 1972 quoted second-hand by Vega. Reference only.
    ThresholdRecord("hexanal", "paraffin_oil", 120.0, 22.0, "not_stated",
                    "classical dilution", "Guadagni 1972 via vega1994 pp. 234-235",
                    cross_study_cross_method=True,
                    provenance_flag="quoted_second_hand_by_vega1994"),
    ThresholdRecord("heptanal", "paraffin_oil", 250.0, 22.0, "not_stated",
                    "classical dilution", "Guadagni 1972 via vega1994",
                    cross_study_cross_method=True,
                    provenance_flag="quoted_second_hand_by_vega1994"),
    ThresholdRecord("tt_2_4_decadienal", "paraffin_oil", 135.0, 22.0, "not_stated",
                    "classical dilution", "Guadagni 1972 via vega1994",
                    cross_study_cross_method=True,
                    provenance_flag="quoted_second_hand_by_vega1994"),
)

#: Matrices the layer KNOWS about but for which it deliberately carries NO
#: threshold entries, with the governance reason. A caller asking for one of
#: these gets ``NoMeasuredThreshold`` carrying this reason, not a borrowed
#: number from a different matrix.
SEALED_OR_REFUSED_MATRICES: Mapping[str, str] = {
    "soy_paste_hong": (
        "Hong 2020's 10 paired soy/water thresholds are the GATING HOLD-OUT of "
        "declaration Amendment 4. Ingesting them as table entries would spend "
        "the hold-out. The layer PREDICTS this matrix; it does not look it up."),
    "caseinate_1pct": (
        "Leksrisompong 2010's 24 BETs are not fit-eligible: Amendment 4 makes "
        "only its K RATIOS FIT, and k4b sec. H.2 proposes the BETs as hold-out. "
        "Its K ratios are used as binding constants; its thresholds are not "
        "used at all."),
    "emulsion_10pct_fat": "as caseinate_1pct; and fat is confounded with "
                          "aqueous-phase protein by the paper's own design.",
    "emulsion_20pct_fat": "as emulsion_10pct_fat.",
    "soybean_oil": "as caseinate_1pct.",
    "cooked_beef": (
        "Brewer 1995 is declaration D.6 Module 7 HOLD-OUT and RECLASSIFIED "
        "`dose_added_pre_cook`: its numbers are doses added to RAW beef before a "
        "70 C cook, not concentrations present at the moment of perception. They "
        "are not thresholds in this repo's sense and are not tabulated here."),
    "milk_tian": (
        "UNIT RESOLVED 2026-08-30 (Wave B8, FIT_HOLDOUT_DECLARATION.md "
        "Amendment 17 clause 6); STILL NOT TABULATED, for a DIFFERENT and now "
        "stated reason. The block was: `Tian 2020's units cell prints a literal "
        "`?`, verified at 900 dpi. A factor-of-1000 basis risk; declared "
        "`neither, until the unit is settled`.` THE UNIT IS SETTLED BY "
        "ARITHMETIC, not by inference: Tian 2020 Table 1's own footnote says "
        "its concentration column is the Y6 sample of Tian et al. 2019, and "
        "SEVEN OF SEVEN values reproduce that column digit for digit (347; "
        "7,030; 1,719 vs 1,720; 29 vs 29.3; 198 vs 199; 244; 1,001 vs 1,000) -- "
        "a column its source heads `Concentration (ug/kg)` in plain type. The "
        "same `?/kg` notation heads the THRESHOLD column of the same table, so "
        "`?/kg` = `ug/kg` and the factor-of-1000 basis risk is CLOSED "
        "(k3 sec. C.18 RESOLVED). Two independent confirmations: k2 sec. A.4 "
        "recovered nonanal's aqueous threshold as ~1.1 ug/kg from Xin2026b and "
        "computed ~1,450x, against 1,455x from Tian's own 1.10 -- 0.3% apart; "
        "and tian2020b's yogurt thresholds (5,430-29,000 ug/kg) span the same "
        "range as the seven milk rows read as ug/kg, while ng/kg would be "
        "30-300x LOWER and mg/kg 1,000x higher and physically absurd. "
        "WHY THEY ARE STILL NOT TABULATED, and this is a role decision and not "
        "a defect: k6b sec. 7b proposes the seven newly-usable rows as a "
        "HOLD-OUT, not as fit data -- 'they were previously unusable, so they "
        "are newly valuable as a hold-out and should NOT be spent as fit data' "
        "-- and no amendment has assigned them a column. THREE CAVEATS MUST "
        "TRAVEL WITH THEM WHEREVER THEY LAND: (1) `same_matrix: FALSE` -- "
        "'milk fan' is a CHEESE and the thresholds were dosed into an "
        "unspecified 'fresh milk solution'; (2) the reference thresholds are "
        "UNSOURCED and mislabelled 'in air' (they are almost certainly "
        "aqueous, which is why nonanal's 1.10 matches); (3) the printed SDs "
        "are NOT threshold uncertainty -- 0.16-2.1% RSD on a 15-panellist "
        "sensory threshold is impossible and five of the seven values are "
        "exact powers-of-two dilution steps (1,024 = 2^10; 25,600 = 25 x "
        "1,024; 51,200 = 50 x 1,024). NO VALUE CHANGES IN THIS WAVE: nothing "
        "was ever tabulated, so unblocking is a provenance edit."),
    "saliva": (
        "Starkenmann's thiol binding is STRANDED (no basis, no stoichiometry, "
        "mechanism unresolved) and Baek 1999 EXCLUDES mucosal binding for a "
        "neutral ester. A single saliva factor applied to all odourants is "
        "refuted by the pair."),
}


def select_threshold(
    compound: str,
    matrix: str = "water",
    temperature_c: Optional[float] = None,
) -> Union[ThresholdRecord, NoMeasuredThreshold]:
    """
    The threshold for ``compound`` in ``matrix``, or an explicit
    ``NoMeasuredThreshold``. NOTHING IS EVER BORROWED FROM ANOTHER MATRIX.

    Matrix selection is by DECLARED VALIDITY DOMAIN: a record is eligible only
    if its ``matrix`` matches exactly. Where the domain is also temperature-
    resolved (the Vega ladder), the nearest measured temperature is used and the
    distance is reported by ``select_threshold_verbose``.
    """
    result, _ = select_threshold_verbose(compound, matrix, temperature_c)
    return result


def select_threshold_verbose(
    compound: str,
    matrix: str = "water",
    temperature_c: Optional[float] = None,
) -> Tuple[Union[ThresholdRecord, NoMeasuredThreshold], Dict[str, object]]:
    """``select_threshold`` plus the selection diagnostics."""
    pool = WATER_THRESHOLDS if matrix == "water" else MATRIX_THRESHOLDS
    candidates = [r for r in pool if r.compound == compound and r.matrix == matrix]
    if not candidates:
        reason = SEALED_OR_REFUSED_MATRICES.get(
            matrix,
            f"no measured threshold for {compound!r} in {matrix!r} anywhere in the corpus")
        return NoMeasuredThreshold(compound, matrix, reason), {
            "n_candidates": 0, "borrowed": False}
    if temperature_c is not None:
        with_t = [r for r in candidates if r.temperature_c is not None]
        if with_t:
            chosen = min(with_t, key=lambda r: abs(r.temperature_c - temperature_c))
            return chosen, {
                "n_candidates": len(candidates),
                "temperature_requested_c": temperature_c,
                "temperature_used_c": chosen.temperature_c,
                "temperature_gap_c": abs(chosen.temperature_c - temperature_c),
                "borrowed": False,
            }
    # Preference order when several sources exist for the same cell: a printed
    # and independently reproduced value, then a doubly corroborated derivation,
    # then anything else -- and the SPREAD across all candidates is reported.
    order = {"printed_and_independently_reproduced": 0,
             "derived_but_doubly_corroborated": 1,
             "ALREADY_AN_IN_MOUTH_NUMBER": 2}
    chosen = min(candidates, key=lambda r: order.get(r.provenance_flag, 3))
    values = [r.value_ug_per_l for r in candidates]
    return chosen, {
        "n_candidates": len(candidates),
        "alternates_ug_per_l": values,
        "spread_x": (max(values) / min(values)) if min(values) > 0 else float("inf"),
        "borrowed": False,
    }


# ===========================================================================
# 2. ABSOLUTES ARE INTERVALS
# ===========================================================================

@dataclass(frozen=True)
class AbsoluteConcentration:
    """
    An absolute concentration in ug/L, as an INTERVAL. There is deliberately no
    ``__float__``: an absolute cannot leave this module as a bare point.

    The band is the MEASURED reliability of the quantity, not a fitted error:
      * HS-SPME same-sample dispersion 10-23x (Xin 2026 vs Xin 2026b on the same
        samples; declaration D.6 Module 8 calls this a CALIBRATION FACT and
        forbids scoring it). A single measurement therefore sits within
        x/ sqrt(23) = 4.80 of the truth, two-sided: 0.681 decades.
      * +/-0.5 decades on the air/water partition constant, added in quadrature
        only when the value passed through a partition step (k4b sec. D.2:
        the total literature spread on hexanal K_aw is 9.5x, and the ruling is
        explicitly to CARRY the band, not to swap the constant).
    """

    point_ug_per_l: float
    lo_ug_per_l: float
    hi_ug_per_l: float
    halfwidth_decades: float
    components: Mapping[str, float]
    via_partition: bool
    provenance: str

    @property
    def band_x(self) -> float:
        """The full multiplicative width of the interval."""
        return self.hi_ug_per_l / self.lo_ug_per_l if self.lo_ug_per_l > 0 else float("inf")

    def as_dict(self) -> Dict[str, object]:
        return {
            "point_ug_per_L": self.point_ug_per_l,
            "interval_ug_per_L": [self.lo_ug_per_l, self.hi_ug_per_l],
            "halfwidth_decades": self.halfwidth_decades,
            "band_x": self.band_x,
            "components_decades": dict(self.components),
            "via_partition": self.via_partition,
            "provenance": self.provenance,
            "warning": "NEVER quote the point without this interval.",
        }


def absolute_concentration(
    point_ug_per_l: float,
    via_partition: bool = True,
    extra_decades: float = 0.0,
    provenance: str = "",
) -> AbsoluteConcentration:
    """Wrap a point value in its measured reliability band."""
    if point_ug_per_l < 0:
        raise ValueError("concentration must be non-negative")
    dispersion_decades = math.log10(math.sqrt(HS_SPME_SAME_SAMPLE_DISPERSION[1]))
    components = {"hs_spme_same_sample_dispersion": dispersion_decades}
    total_sq = dispersion_decades ** 2
    if via_partition:
        components["air_water_partition_constant"] = K_AW_UNCERTAINTY_DECADES
        total_sq += K_AW_UNCERTAINTY_DECADES ** 2
    if extra_decades:
        components["caller_supplied"] = extra_decades
        total_sq += extra_decades ** 2
    halfwidth = math.sqrt(total_sq)
    factor = 10.0 ** halfwidth
    return AbsoluteConcentration(
        point_ug_per_l=point_ug_per_l,
        lo_ug_per_l=point_ug_per_l / factor,
        hi_ug_per_l=point_ug_per_l * factor,
        halfwidth_decades=halfwidth,
        components=components,
        via_partition=via_partition,
        provenance=provenance or "B4 measured reliability band",
    )


# ===========================================================================
# 3. THE OAV LAYER
# ===========================================================================

@dataclass(frozen=True)
class OAVResult:
    compound: str
    matrix: str
    threshold_ug_per_l: Optional[float]
    threshold_source: Optional[str]
    threshold_state: str
    oav_point: Optional[float]
    oav_lo: Optional[float]
    oav_hi: Optional[float]
    concentration: Optional[AbsoluteConcentration]
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, object]:
        return {
            "compound": self.compound, "matrix": self.matrix,
            "threshold_ug_per_L": self.threshold_ug_per_l,
            "threshold_source": self.threshold_source,
            "threshold_state": self.threshold_state,
            "OAV_point": self.oav_point,
            "OAV_interval": (None if self.oav_lo is None else [self.oav_lo, self.oav_hi]),
            "concentration": (None if self.concentration is None
                              else self.concentration.as_dict()),
            "diagnostics": dict(self.diagnostics),
        }


def odour_activity(
    compound: str,
    concentration: Union[AbsoluteConcentration, float],
    matrix: str = "water",
    temperature_c: Optional[float] = None,
    via_partition: bool = True,
) -> OAVResult:
    """
    OAV = concentration / threshold, per species, as an INTERVAL.

    A bare float is accepted and immediately wrapped in its reliability band --
    it is never divided as a point.
    """
    if not isinstance(concentration, AbsoluteConcentration):
        concentration = absolute_concentration(
            float(concentration), via_partition=via_partition,
            provenance=f"auto-wrapped for {compound}")
    record, diagnostics = select_threshold_verbose(compound, matrix, temperature_c)
    if isinstance(record, NoMeasuredThreshold):
        return OAVResult(
            compound=compound, matrix=matrix, threshold_ug_per_l=None,
            threshold_source=None,
            threshold_state="no_measured_threshold_for_this_matrix",
            oav_point=None, oav_lo=None, oav_hi=None, concentration=concentration,
            diagnostics={**diagnostics, "reason": record.reason,
                         "borrowed_from_another_matrix": False})
    threshold = record.value_ug_per_l
    return OAVResult(
        compound=compound, matrix=matrix, threshold_ug_per_l=threshold,
        threshold_source=record.source, threshold_state="measured",
        oav_point=concentration.point_ug_per_l / threshold,
        oav_lo=concentration.lo_ug_per_l / threshold,
        oav_hi=concentration.hi_ug_per_l / threshold,
        concentration=concentration,
        diagnostics={**diagnostics,
                     "criterion": record.criterion,
                     "provenance_flag": record.provenance_flag,
                     "cross_study_cross_method": record.cross_study_cross_method,
                     "role": record.role})


def oav_table(
    concentrations_ug_per_l: Mapping[str, Union[AbsoluteConcentration, float]],
    matrix: str = "water",
    temperature_c: Optional[float] = None,
) -> Dict[str, object]:
    """
    OAV for every species supplied, INCLUDING the MFT dimer at its own
    threshold.

    THE DIMER IS THE POINT. bis(2-methyl-3-furyl) disulfide has a water
    threshold of 0.00032 ug/L against MFT's 0.005 -- it is 15.6x MORE POTENT
    than its own monomer. Only 6.5-9.6 % of MFT-equivalents sit in the dimer,
    but 0.096/0.064 > 1, so the dimer's OAV MATCHES or EXCEEDS the monomer's
    (Zhou SI Table S2 prints 3.21e5 against MFT's 3.18e5 at pH 7). Mass lost to
    dimerisation is NOT aroma lost, and any objective that scores the
    dimerisation channel as a pure loss is wrong by roughly the threshold ratio.
    """
    results: Dict[str, object] = {}
    for compound, concentration in concentrations_ug_per_l.items():
        results[compound] = odour_activity(
            compound, concentration, matrix, temperature_c).as_dict()
    scored = {k: v for k, v in results.items()
              if isinstance(v, dict) and v.get("OAV_point") is not None}
    ranked = sorted(scored.items(), key=lambda kv: kv[1]["OAV_point"], reverse=True)
    out: Dict[str, object] = {
        "matrix": matrix,
        "temperature_c": temperature_c,
        "per_species": results,
        "ranking_by_OAV": [k for k, _ in ranked],
        "n_without_threshold": len(results) - len(scored),
        "dimer_potency_ratio_over_monomer": (
            _ZHOU_S2["MFT"] / _ZHOU_S2["MFTD"]),
        "dimer_note": (
            "The MFT dimer is tracked as its own species at its own threshold. "
            "Dimerisation is NOT aroma loss."),
        "threshold_policy": (
            "Thresholds are model INPUTS, never scored predictions "
            "(declaration D.6 Module 7). No value is borrowed across matrices."),
    }
    if "MFT" in scored and "MFTD" in scored:
        out["dimer_over_monomer_OAV"] = (
            scored["MFTD"]["OAV_point"] / scored["MFT"]["OAV_point"]
            if scored["MFT"]["OAV_point"] > 0 else float("inf"))
    return out


# ===========================================================================
# 4. RATIOS LEAD -- formulation comparison, which is the layer's headline
# ===========================================================================

@dataclass(frozen=True)
class FormulationComparison:
    compound: str
    ratio: float
    direction: str
    within_reliability_band: bool
    band_x: float
    note: str

    def as_dict(self) -> Dict[str, object]:
        return {"compound": self.compound, "ratio_a_over_b": self.ratio,
                "direction": self.direction,
                "within_reliability_band": self.within_reliability_band,
                "reliability_band_x": self.band_x, "note": self.note}


def compare_formulations(
    formulation_a: Mapping[str, float],
    formulation_b: Mapping[str, float],
    label_a: str = "A",
    label_b: str = "B",
) -> Dict[str, object]:
    """
    The layer's PRIMARY output: per-compound ratios between two formulations,
    with an explicit validity bound.

    A ratio is the right unit because the two dominant error sources -- the
    HS-SPME calibration offset and the air/water partition constant -- are
    SHARED between the two arms and cancel exactly in a within-run ratio. That
    is the same argument that makes Meynier's and Leksrisompong's ratios
    fit-eligible while their absolutes are declared `absolute_scale_suspect`.

    A ratio inside the same-sample dispersion band is reported as NOT
    RESOLVED rather than as a small effect.
    """
    band_x = math.sqrt(HS_SPME_SAME_SAMPLE_DISPERSION[1])
    rows: List[FormulationComparison] = []
    shared = [k for k in formulation_a if k in formulation_b]
    for compound in shared:
        a, b = float(formulation_a[compound]), float(formulation_b[compound])
        if b <= 0:
            rows.append(FormulationComparison(
                compound, float("inf"), "undefined", False, band_x,
                f"{label_b} is zero; a ratio is undefined"))
            continue
        ratio = a / b
        resolved = ratio >= band_x or ratio <= 1.0 / band_x
        rows.append(FormulationComparison(
            compound, ratio,
            "higher_in_" + label_a if ratio > 1 else
            ("higher_in_" + label_b if ratio < 1 else "equal"),
            not resolved, band_x,
            "NOT RESOLVED: inside the measured same-sample dispersion band"
            if not resolved else "resolved above the same-sample dispersion band"))
    # Q1: an UNDEFINED ratio resolves nothing and must not be counted as
    # resolved. It was, because the undefined row is constructed with
    # ``within_reliability_band=False`` -- true in the literal sense (an
    # undefined ratio is not inside the dispersion band) but read by every
    # consumer as "this row resolved". "6 of 6 ratios resolve" then reads as six
    # claims when two of them are not ratios at all. Both the CLI and the HTML
    # report had grown their own compensating text; the count is fixed here
    # instead, in the layer that owns it, and ``n_undefined`` is published so a
    # renderer never has to re-derive the predicate.
    undefined_rows = [r for r in rows if r.direction == "undefined"]
    resolved_rows = [
        r for r in rows
        if r.direction != "undefined" and not r.within_reliability_band
    ]
    return {
        "formulation_a": label_a, "formulation_b": label_b,
        "n_compared": len(rows), "n_resolved": len(resolved_rows),
        "n_undefined": len(undefined_rows),
        "reliability_band_x": band_x,
        "rows": [r.as_dict() for r in rows],
        "ranking_by_ratio": [r.compound for r in sorted(
            rows, key=lambda r: r.ratio, reverse=True)],
        "policy": ("RATIOS LEAD. Shared multiplicative error cancels in a "
                   "within-run ratio; absolutes carry the full band and are "
                   "never emitted as bare points."),
    }


# ===========================================================================
# 5. THE THREE-TERM MATRIX SHIFT MODEL, AND THE HONEST RESIDUAL
# ===========================================================================

#: The per-class per-gram binding constants, derived from FIT ROWS ONLY.
#: This IS the fit: it is a deterministic, reproducible function of the frozen
#: registry, so freezing the registry freezes the fit.
def fit_class_binding_constants() -> Dict[str, Dict[str, object]]:
    """
    Pool the FIT-row per-gram constants by binding class, WITHIN METHOD FAMILY.

    Aldehyde constants are never pooled across headspace and dialysis (k2
    sec. B.3). Classes with no FIT-row constant get ``None`` -- an explicit gap,
    not an imputed default.
    """
    by_class: Dict[str, List[MatrixParameter]] = {}
    for parameter in REVERSIBLE_BINDING:
        structure = COMPOUND_STRUCTURE.get(parameter.compound)
        if structure is None:
            continue
        if parameter.provenance.get("quarantined_as_binding"):
            continue
        # Headspace family only: it is the family that measures what a
        # headspace prediction needs, and the two families must not be pooled.
        if parameter.method not in ("headspace_depletion", "static_headspace_partition"):
            continue
        by_class.setdefault(structure.binding_class, []).append(parameter)

    out: Dict[str, Dict[str, object]] = {}
    for chem_class, parameters in by_class.items():
        values = [p.value for p in parameters if p.value is not None]
        positives = [v for v in values if v > 0]
        if positives and len(positives) == len(values):
            # geometric mean, the right average for a multiplicative constant
            pooled = math.exp(sum(math.log(v) for v in positives) / len(positives))
            sign = "suppression"
        else:
            # at least one measured ENHANCEMENT: arithmetic mean, and the class
            # is flagged, because a geometric mean is undefined through zero.
            pooled = sum(values) / len(values)
            sign = "enhancement_measured"
        out[chem_class] = {
            "k_g_l_per_g": pooled,
            "n_fit_rows": len(values),
            "members": [p.compound for p in parameters],
            "reference_loading_g_per_l": max(
                MATRIX_LOADING[p.medium].protein_g_per_l or 0.0
                if p.medium in MATRIX_LOADING else 0.0 for p in parameters),
            "method_family": "headspace",
            "sign": sign,
            "sources": sorted({p.source_anchor.split(" sec.")[0] for p in parameters}),
        }
    # Chain-length surrogate: branched C5 alkanals from the C6 n-alkanal.
    if "n_alkanal" in out:
        out["branched_alkanal"] = {
            "k_g_l_per_g": out["n_alkanal"]["k_g_l_per_g"] / CHAIN_LENGTH_SLOPE_PER_CH2,
            "n_fit_rows": 0,
            "members": [],
            "reference_loading_g_per_l": out["n_alkanal"]["reference_loading_g_per_l"],
            "method_family": "headspace",
            "sign": "suppression",
            "surrogate": (
                "CHAIN-LENGTH EXTRAPOLATED from the C6 n-alkanal by the measured "
                f"{CHAIN_LENGTH_SLOPE_PER_CH2}x/CH2 slope (Andriot 2.72, Damodaran "
                "2.9, agreeing to 6 %). The slope is measured on BINDING "
                "CONSTANTS. It is NOT applied to any threshold shift, because k2 "
                "sec. D.3 shows the chain-length structure of the shift does not "
                "transfer between matrices. Branch position is unmeasured "
                "anywhere, so 2-methyl and 3-methyl butanal get the SAME value "
                "and the layer cannot distinguish them."),
            "sources": out["n_alkanal"]["sources"],
        }
    return out


def fit_unsaturation_penalty() -> Dict[str, object]:
    """
    The alpha,beta-unsaturation penalty, from FIT ROWS ONLY.

    Two observations survive the declaration: Vega's gelatin ladder (2.81x) and
    Meynier's skim-milk headspace contrast (4.95x). Brewer's beef 2.01x is
    EXCLUDED because Brewer is HOLD-OUT and reclassified `dose_added_pre_cook`.
    """
    values = [p.value for p in ALPHA_BETA_UNSATURATION_OBSERVATIONS if p.value]
    pooled = math.exp(sum(math.log(v) for v in values) / len(values))
    return {
        "penalty_x": pooled,
        "n_fit_rows": len(values),
        "observations": {p.key: p.value for p in ALPHA_BETA_UNSATURATION_OBSERVATIONS},
        "excluded": {"unsat_penalty_beef": 2.01},
        "spread_x": max(values) / min(values),
        "caveat": (
            "k2 sec. D.3 states the penalty as ~2-3x. This FIT-only estimate is "
            "ABOVE that band, and it is above it BECAUSE the observation that "
            "would pull it down (beef, 2.01x) is hold-out and excluded. k4b "
            "records the same disagreement: 'the term's magnitude is "
            "matrix-dependent'. Reported, not smoothed."),
        "applies_to": "alpha,beta-unsaturated CARBONYLS only (Michael acceptors). "
                      "A conjugated C=C without a carbonyl -- 4-vinyl phenol -- "
                      "does NOT qualify.",
    }


@dataclass(frozen=True)
class ShiftPrediction:
    compound: str
    matrix: str
    predicted_ratio: float
    predicted_lo: float
    predicted_hi: float
    terms: Mapping[str, object]
    state: str
    warnings: Tuple[str, ...]

    def as_dict(self) -> Dict[str, object]:
        return {
            "compound": self.compound, "matrix": self.matrix,
            "predicted_matrix_over_water_ratio": self.predicted_ratio,
            "predicted_interval": [self.predicted_lo, self.predicted_hi],
            "predicted_sign": ("elevated" if self.predicted_ratio > 1.0 else
                               ("inverted" if self.predicted_ratio < 1.0 else "flat")),
            "terms": dict(self.terms), "state": self.state,
            "warnings": list(self.warnings),
        }


def covalent_channel_state(compound: str, ph: Optional[float] = None) -> Dict[str, object]:
    """
    Whether the covalent channel is even STRUCTURALLY ALLOWED for ``compound``,
    and -- since it is a ceiling and contributes zero -- what it would report.

    Two structural gates, both binding:
      * the 32-of-47 no-adduct negative gate. A model that applies a generic
        protein-binding loss term across all volatile classes is falsified by
        32 counter-examples at a saturating (12 ppth) dose.
      * the pH-3 adduct gate: carbonyl-lysine adduct formation is ABOLISHED at
        pH 3, so any aldehyde-loss-to-protein term must go to zero at acid pH.
    """
    structure = COMPOUND_STRUCTURE.get(compound)
    display = structure.display.lower() if structure else compound.lower()
    gated_by_name = any(display == g.lower() or display.replace("-", " ") == g.lower()
                        or display.replace(" ", "") == g.lower().replace(" ", "")
                        or display.replace("-", "") == g.lower().replace("-", "")
                        for g in NO_ADDUCT_GATE_COMPOUNDS)
    gated_by_class = (structure is not None
                      and structure.binding_class in NO_ADDUCT_GATE_CLASSES)
    allowed_class = (structure is not None
                     and structure.binding_class in ADDUCT_POSITIVE_CLASSES)
    state: str
    if gated_by_name:
        state = "blocked_by_negative_gate_named_compound"
    elif gated_by_class and not allowed_class:
        state = "blocked_by_negative_gate_class"
    elif allowed_class:
        state = "structurally_allowed"
    else:
        state = "unknown_class_no_gate_evidence"
    ph_state = "not_gated"
    if ph is not None:
        if ph <= PH_ADDUCT_GATE_BELOW:
            state = "blocked_by_pH_gate"
            ph_state = f"pH {ph} <= {PH_ADDUCT_GATE_BELOW}: adduct formation ABOLISHED"
        elif ph < PH_ADDUCT_GATE_UNCERTAIN_BELOW:
            ph_state = (f"pH {ph} is between the abolished point ({PH_ADDUCT_GATE_BELOW}) "
                        f"and the demonstrated point ({PH_ADDUCT_GATE_UNCERTAIN_BELOW}); "
                        "no source measures the intermediate, so this is reported as "
                        "uncertain rather than interpolated")
    out: Dict[str, object] = {
        "state": state,
        "ph_state": ph_state,
        "contribution_to_point_prediction": 0.0,
        "why_zero": (
            "INERT BY RULING. Amendment 6 ruling 2: the covalent channel supplies "
            "~0.06 % of the 1 304x hexanal log-shift on a threshold-panel "
            "timescale and does NOT close the matrix gap. It is carried as a "
            "ceiling, not a term."),
        "ceiling_k2_M_per_s_at_20C": COVALENT_CEILING.k2_m_per_s_at_20c,
        "ambient_half_life_days": list(COVALENT_CEILING.ambient_half_life_days),
    }
    if state == "structurally_allowed":
        # WAVE B8, Amendment 17 clause 6. This report used to read "Could matter
        # at process temperature IF Ea >= 70 kJ/mol ... THAT Ea IS UNMEASURED IN
        # EVERY CORPUS SOURCE." It is measured now, and it missed the threshold.
        #
        # The text is a PLAIN STRING constant carried on the parameter object
        # rather than an f-string assembled here, deliberately: on Python 3.12 an
        # f-string tokenizes into FSTRING_MIDDLE parts that are NOT
        # tokenize.STRING, so every digit inside one leaks into what the B4
        # hold-out firewall reads as executable code. Interpolating this
        # paragraph would have tripped the firewall on a temperature in a
        # citation. Keeping the prose in the registry keeps the guard sharp.
        out["process_temperature_report"] = (
            COVALENT_CEILING.process_temperature_verdict)
        out["ea_kJ_per_mol"] = COVALENT_CEILING.ea_kj_mol
        out["ea_range_kJ_per_mol"] = list(COVALENT_CEILING.ea_range_kj_mol or ())
        out["ea_status"] = COVALENT_CEILING.ea_status
        out["retirement"] = COVALENT_CEILING_RETIREMENT
    else:
        out["process_temperature_report"] = (
            "Not applicable: the compound is outside the structurally allowed "
            "set, so there is no covalent channel to activate at any temperature.")
    if compound == "dimethyl_disulfide":
        out["source_contradiction"] = (
            "anantharamkrishnan2020b Table 2 says DMDS forms NO adduct (it is one "
            "of the 32); its own Results text says DMDS DOES, and names +46 Da "
            "(BLG-CysSSMe). Resolved CONSERVATIVELY here in favour of the table. "
            "Reported on every DMDS prediction.")
    return out


def predict_matrix_shift(
    compound: str,
    matrix: str,
    ph: Optional[float] = None,
    class_constants: Optional[Mapping[str, Dict[str, object]]] = None,
    unsaturation: Optional[Mapping[str, object]] = None,
) -> ShiftPrediction:
    """
    Predict the matrix/water threshold-shift ratio for ``compound`` from the
    three named terms, with an interval and an explicit state.

    The state is as important as the number. ``no_binding_constant_for_class``
    means the corpus supplies no term at all for this compound -- the layer
    reports that instead of manufacturing a plausible-looking 1.0.
    """
    constants = dict(class_constants or fit_class_binding_constants())
    penalty = dict(unsaturation or fit_unsaturation_penalty())
    structure = COMPOUND_STRUCTURE.get(compound)
    loading = MATRIX_LOADING.get(matrix)
    warnings: List[str] = []
    terms: Dict[str, object] = {}

    if structure is None:
        return ShiftPrediction(compound, matrix, float("nan"), float("nan"),
                               float("nan"), {}, "unknown_compound",
                               ("no structural record for this compound",))
    if loading is None:
        return ShiftPrediction(compound, matrix, float("nan"), float("nan"),
                               float("nan"), {}, "unknown_matrix",
                               ("no loading record for this matrix",))

    entry = constants.get(structure.binding_class)
    c_mid = loading.protein_g_per_l or 0.0
    c_lo = loading.protein_lo_g_per_l or 0.0
    c_hi = loading.protein_hi_g_per_l or 0.0

    # ---- TERM 1: reversible binding ------------------------------------
    if entry is None:
        terms["reversible_binding"] = {
            "state": "no_measured_constant_for_this_class",
            "chem_class": structure.binding_class,
            "factor_x": None,
            "why": (f"No per-gram binding constant exists anywhere in the FIT "
                    f"column for the class {structure.binding_class!r}. The "
                    f"corpus measures binding for esters, n-alkanals, methyl "
                    f"ketones, diketones, lactones and furanones and for nothing "
                    f"else. Imputing one would be inventing a parameter."),
        }
        state = "no_binding_constant_for_class"
        bind_mid = bind_lo = bind_hi = 1.0
    else:
        k_g = float(entry["k_g_l_per_g"])
        reference_loading = float(entry.get("reference_loading_g_per_l") or 0.0)
        if k_g < 0:
            # A measured ENHANCEMENT. The linear/Langmuir form 1 + K*C is not
            # licensed for it (it goes negative), and no source measures an
            # enhancement at another loading, so the layer REFUSES to
            # extrapolate and holds the shift at the measured value.
            measured_ratio = 1.0 + k_g * reference_loading
            bind_mid = bind_lo = bind_hi = max(measured_ratio, 1e-3)
            warnings.append(
                f"MEASURED ENHANCEMENT for class {structure.binding_class!r} "
                f"(K_g = {k_g:.3g} L/g < 0): the matrix makes this compound MORE "
                f"volatile. The linear binding form is not licensed for a negative "
                f"constant, so the shift is HELD at the measured value at "
                f"{reference_loading:g} g/L and NOT extrapolated to "
                f"{c_mid:g} g/L.")
            terms["reversible_binding"] = {
                "state": "measured_enhancement_extrapolation_refused",
                "chem_class": structure.binding_class, "k_g_l_per_g": k_g,
                "factor_x": bind_mid,
                "reference_loading_g_per_l": reference_loading,
                "sources": entry.get("sources"),
            }
        else:
            bind_mid = 1.0 + k_g * c_mid
            bind_lo = 1.0 + k_g * c_lo
            bind_hi = 1.0 + k_g * c_hi
            if reference_loading > 0 and c_mid > 3.0 * reference_loading:
                warnings.append(
                    f"EXTRAPOLATION WARNING: the protein loading used "
                    f"({c_mid:g} g/L) is {c_mid / reference_loading:.1f}x the "
                    f"loading the constant was measured at ({reference_loading:g} "
                    f"g/L). The linear binding form has been validated only to "
                    f"~2.5x (Andriot's own 4 % -> 10 % extrapolation).")
            terms["reversible_binding"] = {
                "state": "measured", "chem_class": structure.binding_class,
                "k_g_l_per_g": k_g, "factor_x": bind_mid,
                "factor_interval": [bind_lo, bind_hi],
                "protein_g_per_l": c_mid,
                "protein_band_g_per_l": [c_lo, c_hi],
                "reference_loading_g_per_l": reference_loading,
                "n_fit_rows": entry.get("n_fit_rows"),
                "surrogate": entry.get("surrogate"),
                "sources": entry.get("sources"),
            }
        state = "predicted"
        if entry.get("surrogate"):
            state = "predicted_with_chain_length_surrogate"

    # ---- TERM 2: the alpha,beta-unsaturation penalty --------------------
    if structure.alpha_beta_unsaturated_carbonyl:
        unsat = float(penalty["penalty_x"])
        terms["alpha_beta_unsaturation"] = {
            "state": "applied", "factor_x": unsat,
            "n_fit_rows": penalty["n_fit_rows"],
            "spread_x": penalty["spread_x"], "caveat": penalty["caveat"]}
    else:
        unsat = 1.0
        terms["alpha_beta_unsaturation"] = {
            "state": "not_applicable", "factor_x": 1.0,
            "why": ("not an alpha,beta-unsaturated carbonyl. "
                    + (structure.notes or ""))}

    # ---- TERM 3: the covalent ceiling (INERT) ---------------------------
    covalent = covalent_channel_state(compound, ph)
    terms["covalent_ceiling"] = covalent

    # A compound can have NO binding constant and still receive the
    # unsaturation term (every alkenal is in exactly that position: Meynier's
    # t-2-hexenal constant is quarantined as a binding constant because ~22-33 %
    # of it is irreversible chemistry). Say so precisely rather than reporting
    # "no term" when one term did apply.
    if state == "no_binding_constant_for_class" and unsat != 1.0:
        state = "partial_unsaturation_term_only"

    predicted = bind_mid * unsat
    predicted_lo = min(bind_lo, bind_hi) * unsat
    predicted_hi = max(bind_lo, bind_hi) * unsat

    if loading.evidence_class == "declared_assumption":
        warnings.append(
            f"The protein loading for {matrix!r} is a DECLARED ASSUMPTION, not a "
            f"corpus value: {loading.notes.splitlines()[0] if loading.notes else ''} "
            f"Band {c_lo:g}-{c_hi:g} g/L propagated into the interval.")
    if ph is None and matrix != "water":
        warnings.append(
            "pH not supplied and not stated by the matrix's source. The covalent "
            "gate is pH-dependent (abolished at pH 3), so its state is reported "
            "as ungated rather than assumed.")

    return ShiftPrediction(compound, matrix, predicted, predicted_lo, predicted_hi,
                           terms, state, tuple(warnings))


@dataclass(frozen=True)
class ResidualDecomposition:
    compound: str
    measured_ratio: float
    predicted_ratio: float
    measured_decades: float
    explained_decades: float
    residual_decades: float
    residual_x: float
    per_term_decades: Mapping[str, float]
    fold_error: float
    sign_measured: str
    sign_predicted: str
    sign_correct: bool
    flags: Tuple[str, ...]

    def as_dict(self) -> Dict[str, object]:
        return {
            "compound": self.compound,
            "measured_ratio": self.measured_ratio,
            "predicted_ratio": self.predicted_ratio,
            "measured_log10": self.measured_decades,
            "explained_log10": self.explained_decades,
            "residual_log10": self.residual_decades,
            "unexplained_residual_x": self.residual_x,
            "per_term_log10": dict(self.per_term_decades),
            "fold_error": self.fold_error,
            "sign_measured": self.sign_measured,
            "sign_predicted": self.sign_predicted,
            "sign_correct": self.sign_correct,
            "flags": list(self.flags),
        }


def decompose_residual(
    prediction: ShiftPrediction,
    measured_ratio: float,
) -> ResidualDecomposition:
    """
    "Measured shift 132x, model explains Nx, residual unexplained Mx" -- per
    compound, in decades, with the terms summing EXACTLY to the explained part.

    The self-check that matters: Amendment 6 ruling 2 caps reversible binding at
    ~25 % of an observed log-shift. If the reversible term claims more than that,
    the layer FLAGS it, because a term that outruns its own evidence ceiling is
    a bug in the layer, not a discovery about the matrix.
    """
    per_term: Dict[str, float] = {}
    for name in ("reversible_binding", "alpha_beta_unsaturation"):
        term = prediction.terms.get(name, {})
        factor = term.get("factor_x") if isinstance(term, dict) else None
        per_term[name] = math.log10(factor) if factor and factor > 0 else 0.0
    per_term["covalent_ceiling"] = 0.0  # inert, by ruling

    explained = sum(per_term.values())
    measured_decades = math.log10(measured_ratio) if measured_ratio > 0 else float("nan")
    residual = measured_decades - explained
    flags: List[str] = []

    if measured_decades != 0 and abs(measured_decades) > 1e-9:
        share = per_term["reversible_binding"] / measured_decades
        if share > REVERSIBLE_LOG_SHIFT_CEILING and measured_decades > 0:
            flags.append(
                f"REVERSIBLE TERM EXCEEDS ITS EVIDENCE CEILING: it claims "
                f"{share:.1%} of the observed log-shift against Amendment 6's "
                f"~{REVERSIBLE_LOG_SHIFT_CEILING:.0%} cap.")
    if prediction.state == "no_binding_constant_for_class":
        flags.append(
            "NO TERM AVAILABLE: the entire measured shift is unexplained "
            "residual, because the corpus supplies no constant for this class.")
    elif prediction.state == "partial_unsaturation_term_only":
        flags.append(
            "PARTIAL: the unsaturation term applied but no binding constant "
            "exists for this class, so the binding contribution is zero by "
            "absence of evidence rather than by measurement.")

    sign_measured = "elevated" if measured_ratio > 1.0 else (
        "inverted" if measured_ratio < 1.0 else "flat")
    sign_predicted = "elevated" if prediction.predicted_ratio > 1.0 else (
        "inverted" if prediction.predicted_ratio < 1.0 else "flat")
    if prediction.predicted_ratio == 1.0:
        sign_predicted = "no_prediction"
        flags.append("SIGN NOT PREDICTED: the model emits exactly 1.0, which is "
                     "the absence of a prediction, not a prediction of no shift.")

    fold = (max(measured_ratio, prediction.predicted_ratio)
            / min(measured_ratio, prediction.predicted_ratio)
            if min(measured_ratio, prediction.predicted_ratio) > 0 else float("inf"))

    return ResidualDecomposition(
        compound=prediction.compound,
        measured_ratio=measured_ratio,
        predicted_ratio=prediction.predicted_ratio,
        measured_decades=measured_decades,
        explained_decades=explained,
        residual_decades=residual,
        residual_x=10.0 ** residual,
        per_term_decades=per_term,
        fold_error=fold,
        sign_measured=sign_measured,
        sign_predicted=sign_predicted,
        sign_correct=(sign_measured == sign_predicted),
        flags=tuple(flags),
    )


# ===========================================================================
# 6. THE NAMED UNIMPLEMENTED TERM
# ===========================================================================
#: Registered here rather than implemented, so that the residual has a named
#: leading candidate instead of an anonymous gap -- and so that adding it later
#: is a governed decision rather than a quiet tuning step.
UNIMPLEMENTED_CANDIDATE_TERMS: Mapping[str, str] = {
    "lipid_phase_partition": (
        "Leksrisompong's oil arm (FIT: K ratios) measures a 170x suppression for "
        "a log P 3.4 lactone and a 1.85-1.94x ENHANCEMENT for two hydrophilic "
        "odourants over the same oil. A lipid term is therefore real, "
        "sign-switching, and fittable on FIT rows. It is NOT implemented here "
        "for two reasons: the build spec fixes this layer at THREE named terms, "
        "and the fat content of the one matrix it would be applied to (Hong's "
        "soy paste) is unreported, so it would be fitted on FIT rows and then "
        "evaluated at an invented loading. Leading candidate for part of the "
        "residual on the lipophilic members of the panel."),
    "background_odour_masking": (
        "The only mechanism proposed anywhere in the corpus that is "
        "compound-specific in a way that IGNORES hydrophobicity, and the only "
        "one that can produce an INVERSION. It is not implementable: Hong "
        "measured no background volatile concentrations in the soy matrix, so "
        "there is nothing to compute a masking term from. This is the single "
        "largest named reason the layer is expected to miss an inversion."),
    "delivery_kinetics": (
        "Baek 1999: perceived intensity tracks the RATE of release "
        "(Imax/2(t75-t25), r = 0.968) better than the maximum concentration "
        "(r = 0.860). No equilibrium partition model computes a time derivative "
        "at all. Would require a release-kinetics state this layer does not have."),
    "criterion_bias": (
        "Vega's 75 % UNCORRECTED criterion yields a systematically HIGHER "
        "threshold than a 50 % forced-choice one. Sign known, magnitude not "
        "(~2x is a plausible size, k2 sec. D.2(i)). Not implemented because no "
        "source measures it; it is a known contributor to every cross-method "
        "ratio in the corpus."),
}


def layer_metadata() -> Dict[str, object]:
    """Everything a report needs about the layer as built."""
    return {
        "wave": "B4 -- matrix / OAV output layer",
        "primary_output": "formulation-vs-formulation ratios and rankings",
        "n_water_thresholds": len(WATER_THRESHOLDS),
        "n_matrix_thresholds": len(MATRIX_THRESHOLDS),
        "matrices_with_thresholds": sorted({r.matrix for r in MATRIX_THRESHOLDS}),
        "matrices_sealed_or_refused": dict(SEALED_OR_REFUSED_MATRICES),
        "class_binding_constants": fit_class_binding_constants(),
        "unsaturation_penalty": fit_unsaturation_penalty(),
        "named_terms": ["reversible_binding", "alpha_beta_unsaturation",
                        "covalent_ceiling (INERT, contributes 0)"],
        "unimplemented_candidate_terms": dict(UNIMPLEMENTED_CANDIDATE_TERMS),
        "absolute_policy": (
            "every absolute ug/L is an AbsoluteConcentration interval; the class "
            "has no __float__, so a bare point cannot leave this module"),
    }
