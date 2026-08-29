"""
src/kinetic_core/parameters_matrix.py

THE PARAMETER REGISTRY OF THE MATRIX / OAV OUTPUT LAYER (Build Wave B4).
========================================================================

Module 7 of the kinetic-core rebuild: the layer that turns predicted
headspace-relevant concentrations into something a plant-protein scientist can
act on -- formulation-vs-formulation RATIOS and RANKINGS first, odour-activity
values second, and absolute ppb only ever as an interval.

WHAT IS IN HERE
---------------
Three, and only three, named matrix-shift terms, each capped by the evidence
that licenses it:

  1. REVERSIBLE BINDING -- per-gram constants, with ``method`` and ``pH`` as
     first-class fields (k2_matrix_and_thresholds.md sec. B.3: headspace-derived
     and dialysis-derived ALDEHYDE constants must never be pooled). Capped at
     ~25 % of an observed log-shift by Amendment 6 ruling 2.
  2. The ALPHA,BETA-UNSATURATION PENALTY, ~2-5x, the only parametric matrix
     term k2 sec. D.4.4 judged defensible, fitted here on FIT rows only.
  3. The COVALENT CEILING -- an INERT BOUND. It contributes ZERO to every point
     prediction and appears only in REPORTING. Through B7 it reported "this
     could matter at process temperature if Ea >= 70 kJ/mol, and that Ea is
     unmeasured in every corpus source" (Amendment 6 ruling 2, a named wet-lab
     gap). WAVE B8 RETIRES THAT: the Ea is MEASURED at 15-23 kJ/mol
     (Amendment 17 clause 6), the 70 threshold is missed by 3.5-4.7x, and the
     channel is now reported as MEASURED NEGLIGIBLE at process temperature and
     REAL over ambient storage. The term's numerical contribution was zero
     before and is zero after; only the reason changed.

Anything an observation shows beyond these three is reported as UNEXPLAINED
RESIDUAL, quantified per compound. That is the honest output and it is the
design the evidence forces: k2 sec. D.1 (ratios span 2000x, 1-sigma band 27-41x
cross-study, 4.7x same-panel), k4b sec. B (partition and binding refuted THREE
independent ways on matched samples), Amendment 6 (covalent supplies ~0.06 % of
the hexanal log-shift on a threshold-panel timescale).

WHAT IS DELIBERATELY *NOT* IN HERE
----------------------------------
* **No general matrix correction factor.** k2 sec. D.1 computed what one would
  cost: a uniform 33x misplaces the two extreme compounds by 10x and 28x in
  OPPOSITE directions. There is a lookup table with an explicit
  ``no_measured_threshold`` state instead.
* **No partition-derived threshold.** Amendment 6 and k4b sec. B: refuted three
  ways on matched samples (Hong r = -0.05 vs log P; Leksrisompong's partition
  model errs up to 30.6x with 2 of 3 signs wrong; Meynier's two-phase model
  fails in both directions on two structural isomers).
* **No log P term of any kind.** k4b hold-out guard #4: "any shipped matrix term
  that is a monotone function of log P is refuted".
* **No `protein_source_registry.json` values.** That file is declared
  ``no_verifiable_source`` and the OLD lane (``src/matrix_correction.py``,
  ``src/headspace.py``) raises RuntimeWarnings about it. Nothing here imports
  it, reproduces it, or re-derives its mocked protein differentiation.
  ``assert_no_mocked_protein_source()`` enforces that at import.
* **No DFT.** Standing owner policy. ``assert_no_dft_matrix()`` at import.
* **The Barallat-Perez lupin and mucin constants.** They are Module 6
  **HOLD-OUT** (declaration D.6). Their KEYS are registered so a caller can see
  the gap; their VALUES are not carried in this file.

Units: per-gram binding constants in L/g protein; protein loadings in g/L;
thresholds in ug/L (== ug/kg for these dilute aqueous matrices, and the sources
print both); temperatures in degrees C.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Mapping, Optional, Tuple

# ---------------------------------------------------------------------------
# Evidence classes. As in parameters.py there is NO computational class.
# ---------------------------------------------------------------------------
MATRIX_EVIDENCE_CLASSES: Tuple[str, ...] = (
    "measured_ratio",              # a within-study ratio printed or computed from printed cells
    "derived_from_fit_rows",       # estimated here from FIT rows only
    "measured_bound",              # a one-sided ceiling or floor that was measured
    "structural_gate",             # a binary/ordinal presence-absence constraint
    "declared_assumption",         # NOT from the corpus; stated, bounded, propagated
    "reference_input",             # a model input (e.g. a threshold), never a scored target
)


@dataclass(frozen=True)
class MatrixParameter:
    """One matrix-layer constant, with everything needed to defend or refuse it."""

    key: str
    quantity: str
    value: Optional[float]
    unit: str
    compound: str
    medium: str
    #: 'headspace_depletion' | 'equilibrium_dialysis' | 'gel_filtration' |
    #: 'static_headspace_partition' | 'sensory_BET' | 'not_applicable'
    method: str
    ph_of_measurement: Optional[float]
    temperature_c: Optional[float]
    source_anchor: str
    evidence_class: str
    declaration_role: str
    notes: str = ""
    #: The provenance fields Amendment 4 and k4b sec. H.4 make mandatory here.
    provenance: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.evidence_class not in MATRIX_EVIDENCE_CLASSES:
            raise ValueError(
                f"{self.key}: evidence_class {self.evidence_class!r} is not one of "
                f"{MATRIX_EVIDENCE_CLASSES}. There is no computational class, by policy."
            )


# ===========================================================================
# 0. THE COMPOUND STRUCTURAL REGISTRY
# ===========================================================================
# Structure only -- functional groups, chain length, saturation. These are
# facts about the molecule, derivable from its name, and they are what the
# three terms are gated on. NOTHING here is a measured matrix outcome.

@dataclass(frozen=True)
class CompoundStructure:
    key: str
    display: str
    chem_class: str
    n_carbon: int
    #: True only for an alpha,beta-unsaturated CARBONYL (a Michael acceptor).
    #: 4-vinyl phenol has a conjugated C=C and NO carbonyl -> False.
    alpha_beta_unsaturated_carbonyl: bool
    #: n-alkanal / branched_alkanal / alkenal / methyl_ketone / ester /
    #: lactone / furanone / diketone / alkylpyrazine / phenol /
    #: carboxylic_acid / alkylfuran / disulfide / thiol
    binding_class: str
    notes: str = ""


COMPOUND_STRUCTURE: Mapping[str, CompoundStructure] = {
    # -- the Vega / Brewer / Meynier aldehyde ladder ------------------------
    "pentanal": CompoundStructure("pentanal", "pentanal", "aliphatic aldehyde", 5, False, "n_alkanal"),
    "hexanal": CompoundStructure("hexanal", "hexanal", "aliphatic aldehyde", 6, False, "n_alkanal"),
    "heptanal": CompoundStructure("heptanal", "heptanal", "aliphatic aldehyde", 7, False, "n_alkanal"),
    "nonanal": CompoundStructure("nonanal", "nonanal", "aliphatic aldehyde", 9, False, "n_alkanal"),
    "t_2_hexenal": CompoundStructure(
        "t_2_hexenal", "trans-2-hexenal", "alpha,beta-unsaturated aldehyde", 6, True, "alkenal",
        "Michael acceptor for lysine/cysteine. Meynier: 6.88x skim-milk retention "
        "against hexanal's 1.39x, attributed by the authors to covalent chemistry."),
    "t_2_octenal": CompoundStructure(
        "t_2_octenal", "trans-2-octenal", "alpha,beta-unsaturated aldehyde", 8, True, "alkenal"),
    "tt_2_4_decadienal": CompoundStructure(
        "tt_2_4_decadienal", "trans,trans-2,4-decadienal",
        "doubly alpha,beta-unsaturated aldehyde", 10, True, "alkenal",
        "The extreme outlier in all three media in k2 sec. A.8. Doubly conjugated."),
    # -- Meynier's esters ---------------------------------------------------
    "isoamyl_acetate": CompoundStructure("isoamyl_acetate", "isoamyl acetate", "ester", 7, False, "ester"),
    "amyl_acetate": CompoundStructure("amyl_acetate", "amyl acetate", "ester", 7, False, "ester"),
    "ethyl_pentanoate": CompoundStructure("ethyl_pentanoate", "ethyl pentanoate", "ester", 7, False, "ester"),
    # -- Leksrisompong's three ---------------------------------------------
    "diacetyl": CompoundStructure("diacetyl", "2,3-butanedione", "alpha-diketone", 4, False, "diketone"),
    "delta_decalactone": CompoundStructure(
        "delta_decalactone", "delta-decalactone", "lactone", 10, False, "lactone"),
    "furaneol": CompoundStructure("furaneol", "furaneol (HDMF)", "hydroxyfuranone", 6, False, "furanone"),
    # -- the ketone series --------------------------------------------------
    "2_heptanone": CompoundStructure("2_heptanone", "2-heptanone", "methyl ketone", 7, False, "methyl_ketone"),
    "2_octanone": CompoundStructure("2_octanone", "2-octanone", "methyl ketone", 8, False, "methyl_ketone"),
    "2_nonanone": CompoundStructure("2_nonanone", "2-nonanone", "methyl ketone", 9, False, "methyl_ketone"),
    # -- the sulfur products the B2 network emits ---------------------------
    "MFT": CompoundStructure("MFT", "2-methyl-3-furanthiol", "furanthiol", 5, False, "thiol"),
    "FFT": CompoundStructure("FFT", "2-furfurylthiol", "furanmethanethiol", 5, False, "thiol"),
    "MFTD": CompoundStructure(
        "MFTD", "bis(2-methyl-3-furyl) disulfide", "disulfide", 10, False, "disulfide",
        "The MFT DIMER. 15.6x MORE potent than its own monomer. Carried as a "
        "species with its own threshold, never scored as aroma loss."),
    "ACTZ": CompoundStructure("ACTZ", "2-acetylthiazole", "acetylazole", 5, False, "methyl_ketone"),
    "FUR": CompoundStructure("FUR", "furfural", "furanaldehyde", 5, False, "n_alkanal"),
    # -- the Hong 2020 panel: STRUCTURE ONLY, from the public manifest ------
    # (results/validation/holdout_frozen/hong2020_public_manifest.json, which
    #  is firewalled to carry no threshold, ratio, sign or direction).
    "2_3_dimethylpyrazine": CompoundStructure(
        "2_3_dimethylpyrazine", "2,3-dimethyl pyrazine", "alkylpyrazine", 6, False, "alkylpyrazine"),
    "ethyl_4_methylpentanoate": CompoundStructure(
        "ethyl_4_methylpentanoate", "ethyl-4-methylpentanoate", "ester", 8, False, "ester",
        "Branched C6 acyl chain esterified to ethanol; the closest structural "
        "homologue in the corpus is Meynier's ethyl pentanoate (same ester, "
        "unbranched, one carbon shorter)."),
    "2_pentylfuran": CompoundStructure(
        "2_pentylfuran", "2-pentyl furan", "alkylfuran", 9, False, "alkylfuran"),
    "4_vinylphenol": CompoundStructure(
        "4_vinylphenol", "4-vinyl phenol", "vinyl phenol", 8, False, "phenol",
        "Conjugated C=C but NO carbonyl -> NOT an alpha,beta-unsaturated "
        "carbonyl, so the unsaturation penalty does NOT apply."),
    "3_methylbutanal": CompoundStructure(
        "3_methylbutanal", "3-methyl butanal", "branched aldehyde", 5, False, "branched_alkanal"),
    "2_methylbutanal": CompoundStructure(
        "2_methylbutanal", "2-methyl butanal", "branched aldehyde", 5, False, "branched_alkanal"),
    "butyric_acid": CompoundStructure(
        "butyric_acid", "butyric acid", "carboxylic acid", 4, False, "carboxylic_acid",
        "The only ionisable member of the Hong panel; pH-dependent speciation, "
        "and the soy matrix pH is unreported."),
    "4_ethylphenol": CompoundStructure(
        "4_ethylphenol", "4-ethyl phenol", "alkyl phenol", 8, False, "phenol"),
    "dimethyl_disulfide": CompoundStructure(
        "dimethyl_disulfide", "dimethyl disulfide", "disulfide", 2, False, "disulfide",
        "NOT a free thiol. Carries the anantharamkrishnan2020b Table-2-vs-text "
        "contradiction (table: no adduct; text: names a +46 Da adduct)."),
}


# ===========================================================================
# 1. PROTEIN AND MATRIX LOADINGS
# ===========================================================================
# Every reversible-binding prediction multiplies a per-gram constant by one of
# these. Where the source did not measure it, that is recorded, and where the
# loading is a DECLARED ASSUMPTION it is bounded and the bound is propagated
# into the prediction interval rather than hidden in a point value.

@dataclass(frozen=True)
class MatrixLoading:
    key: str
    display: str
    protein_g_per_l: Optional[float]
    protein_lo_g_per_l: Optional[float]
    protein_hi_g_per_l: Optional[float]
    fat_volume_fraction: Optional[float]
    ph: Optional[float]
    source_anchor: str
    evidence_class: str
    notes: str = ""


MATRIX_LOADING: Mapping[str, MatrixLoading] = {
    "water": MatrixLoading(
        "water", "distilled / ultrapure water", 0.0, 0.0, 0.0, 0.0, None,
        "definition", "measured_ratio", "The reference medium of every ratio."),
    "skim_milk": MatrixLoading(
        "skim_milk", "UHT skim milk (Lactel)", 33.9, 33.9, 33.9, 0.003, None,
        "Meynier2002_extraction.md sec. 2: caseins 29.5 + beta-lg 3.2 + "
        "alpha-la 1.2 g/L, CITED from Swaisgood 1995, NOT measured by Meynier. "
        "Lipid stated only as '< 0.3 %'. pH not stated.",
        "measured_ratio",
        "Protein content is [C] not [M]; carried as a point because all three "
        "cited components are conventional and mutually consistent."),
    "caseinate_1pct": MatrixLoading(
        "caseinate_1pct", "1 % w/v calcium caseinate", 10.0, 10.0, 10.0, 0.0, 7.0,
        "leksrisompong2010_extraction.md sec. 2.1, Table 1.", "measured_ratio",
        "pH 5.5 and 7.0 arms both exist; the pH-7.0 arm is the one used here "
        "because the aldehyde-adduct chemistry is pH-gated (see PH_ADDUCT_GATE)."),
    "gelatin_3pct": MatrixLoading(
        "gelatin_3pct", "3 % w/v gelatin gel", 30.0, 30.0, 30.0, 0.0, None,
        "vega1994 via k2_matrix_and_thresholds.md sec. A.2: 3 % w/v gelatin in "
        "distilled water = 30 g protein/L, LIPID-FREE, dosed at 22 C then held "
        "18 h at 4 C, i.e. NO thermal step after dosing.", "measured_ratio"),
    "soy_paste_hong": MatrixLoading(
        "soy_paste_hong", "autoclaved whole-soybean paste (Hong 2020 model)",
        142.0, 100.0, 200.0, 0.12, None,
        "DECLARED ASSUMPTION. Hong 2020 reports NO composition at all -- pH, "
        "protein, fat, moisture and solids are every one of them unstated "
        "(hong2020_public_manifest.json, matrix_composition: not_stated). The "
        "central value is arithmetic on the printed PREPARATION: 2.5 kg dry "
        "soybeans (conventional ~36-40 % protein, ~20 % oil, dry basis) soaked "
        "in 3 L water for 20 h, autoclaved and ground whole, giving a paste of "
        "roughly 6.0-7.0 kg and therefore ~100-200 g protein/L.",
        "declared_assumption",
        "THE SINGLE LARGEST ASSUMPTION IN THIS LAYER. The composition numbers "
        "are generic food-composition knowledge, NOT corpus values, and the "
        "band is propagated into every soy prediction interval. It is also why "
        "k2 sec. D.4.1's lookup-table key must include matrix LOADING and not "
        "just matrix identity -- and this matrix's loading is unknowable from "
        "its own paper."),
}


# ===========================================================================
# 2. TERM 1 -- REVERSIBLE BINDING, per-gram constants
# ===========================================================================
# `method` and `pH` are first-class because k2 sec. B.3 measured a 35x gap on
# an ALDEHYDE where the ketone gap is 1.3x, between headspace-depletion and
# dialysis determinations, with the mechanism printed by the source: headspace
# depletion counts irreversible cysteine-aldehyde / Schiff chemistry, dialysis
# with 2-mercaptoethanol does not. Pooling them is forbidden by construction:
# `binding_constant_for` refuses to cross the method boundary for aldehydes.

REVERSIBLE_BINDING: Tuple[MatrixParameter, ...] = (
    # ---- from Meynier 2002's within-study RATIOS (Amendment 4: the RATIOS
    #      are FIT; the ABSOLUTE K_aw never is -- static-headspace absolute
    #      scale is suspect by 6.24x, and the offset cancels in a ratio).
    #      K_g = (K_water/K_matrix - 1) / protein_g_per_L, at 33.9 g/L.
    MatrixParameter(
        "kg_isoamyl_acetate_dairy", "per-gram reversible binding constant",
        2.06e-3, "L/g", "isoamyl_acetate", "skim_milk", "static_headspace_partition",
        None, 30.0,
        "Meynier2002_extraction.md sec. 5.1 (water/skim-milk 1.07x) / sec. 2 (33.9 g/L)",
        "measured_ratio", "FIT (Amendment 4: Meynier partition RATIOS)",
        provenance={"absolute_scale_suspect": True, "offset_vs_literature": "6.24x_low",
                    "ratio_cancels_offset": True, "cross_study_cross_method": False}),
    MatrixParameter(
        "kg_amyl_acetate_dairy", "per-gram reversible binding constant",
        5.90e-3, "L/g", "amyl_acetate", "skim_milk", "static_headspace_partition",
        None, 30.0,
        "Meynier2002_extraction.md sec. 5.1 (1.20x) / sec. 2 (33.9 g/L)",
        "measured_ratio", "FIT (Amendment 4)",
        provenance={"absolute_scale_suspect": True, "ratio_cancels_offset": True,
                    "cross_study_cross_method": False}),
    MatrixParameter(
        "kg_ethyl_pentanoate_dairy", "per-gram reversible binding constant",
        9.73e-3, "L/g", "ethyl_pentanoate", "skim_milk", "static_headspace_partition",
        None, 30.0,
        "Meynier2002_extraction.md sec. 5.1 (1.33x) / sec. 2 (33.9 g/L)",
        "measured_ratio", "FIT (Amendment 4)",
        provenance={"absolute_scale_suspect": True, "ratio_cancels_offset": True,
                    "cross_study_cross_method": False,
                    "note": "the closest structural homologue of Hong's "
                            "ethyl-4-methylpentanoate anywhere in the corpus"}),
    MatrixParameter(
        "kg_hexanal_dairy", "per-gram reversible binding constant",
        1.151e-2, "L/g", "hexanal", "skim_milk", "static_headspace_partition",
        None, 30.0,
        "Meynier2002_extraction.md sec. 5.1 (1.39x) / sec. 2 (33.9 g/L)",
        "measured_ratio", "FIT (Amendment 4)",
        provenance={"absolute_scale_suspect": True, "ratio_cancels_offset": True,
                    "cross_study_cross_method": False,
                    "irreversible_share": "<=2% on a headspace-equilibration "
                                          "timescale (Amendment 6 ruling 4)"}),
    MatrixParameter(
        "kg_t_2_hexenal_dairy", "per-gram apparent binding constant",
        1.734e-1, "L/g", "t_2_hexenal", "skim_milk", "static_headspace_partition",
        None, 30.0,
        "Meynier2002_extraction.md sec. 5.1 (6.88x) / sec. 2 (33.9 g/L)",
        "measured_ratio",
        "FIT for the UNSATURATION PENALTY ONLY; QUARANTINED as a binding constant",
        notes="k4b sec. H.3 quarantines Meynier's t-2-hexenal / skim-milk rows: "
              "measured on a visibly browning system at up to 5 800 mg/L. And "
              "Amendment 6 ruling 4 sizes ~22-33 % of this 'partition' as "
              "IRREVERSIBLE CHEMISTRY, so it is not a reversible constant at all. "
              "It is used here only as one of the two FIT-row observations of the "
              "alkenal/alkanal contrast.",
        provenance={"absolute_scale_suspect": True, "quarantined_as_binding": True,
                    "covalent_share_pct": (22.0, 33.0)}),
    # ---- from Leksrisompong 2010's within-study K RATIOS (Amendment 4 FIT).
    #      1 % w/v caseinate = 10 g/L, 40 C. TWO OF THREE ARE NEGATIVE: the
    #      matrix makes the odourant MORE volatile. That is a measurement, and
    #      it is the only thing in the FIT column that lets this layer emit a
    #      shift below 1 at all.
    MatrixParameter(
        "kg_diacetyl_caseinate", "per-gram reversible binding constant",
        7.30e-2, "L/g", "diacetyl", "caseinate_1pct", "static_headspace_partition",
        7.0, 40.0,
        "leksrisompong2010_extraction.md sec. 4.2 (K_water/K_matrix = 1.73x at "
        "0 % fat pH 7.0) / sec. 2.1 (1 % w/v Cas = 10 g/L)",
        "measured_ratio", "FIT (Amendment 4: Leksrisompong K ratios)",
        provenance={"absolute_scale_suspect": True, "offset_vs_literature": "6-17x_low",
                    "ratio_cancels_offset": True, "cross_study_cross_method": False,
                    "ph_dependence_measured": "binds caseinate at pH 7 (P<0.05), "
                                              "NOT at pH 5.5 (N.S.)"}),
    MatrixParameter(
        "kg_delta_decalactone_caseinate", "per-gram reversible binding constant",
        -2.30e-2, "L/g", "delta_decalactone", "caseinate_1pct",
        "static_headspace_partition", 7.0, 40.0,
        "leksrisompong2010_extraction.md sec. 4.2 (K_water/K_matrix = 0.77x)",
        "measured_ratio", "FIT (Amendment 4)",
        notes="NEGATIVE BY MEASUREMENT: 1 % caseinate makes delta-decalactone "
              "MORE volatile, not less. A binding term that can only suppress is "
              "refuted by this cell.",
        provenance={"enhancement_observed": True, "ratio_cancels_offset": True}),
    MatrixParameter(
        "kg_furaneol_caseinate", "per-gram reversible binding constant",
        -6.20e-2, "L/g", "furaneol", "caseinate_1pct", "static_headspace_partition",
        7.0, 40.0,
        "leksrisompong2010_extraction.md sec. 4.2 (K_water/K_matrix = 0.38x)",
        "measured_ratio", "FIT (Amendment 4)",
        notes="NEGATIVE BY MEASUREMENT, and the largest enhancement in the FIT "
              "column. Furaneol is also on the 32-of-47 NO-ADDUCT list, so it "
              "gets no covalent term either.",
        provenance={"enhancement_observed": True, "ratio_cancels_offset": True}),
    # ---- Damodaran 1981 soy (Module 6 FIT; the strongest provenance in the
    #      batch -- the 100 kDa basis is PRINTED, not inferred).
    MatrixParameter(
        "kg_2_heptanone_soy", "per-gram reversible binding constant",
        4.40e-3, "L/g", "2_heptanone", "soy_protein", "equilibrium_dialysis",
        8.0, 25.0,
        "damodaran1981 via k2_matrix_and_thresholds.md sec. (b); n*K/MW, "
        "MW 100 000 g/mol STATED BY SOURCE",
        "measured_ratio", "FIT (declaration D.6 Module 6)",
        provenance={"molar_basis": "stated_by_source", "buffer": "30 mM Tris, 10 mM 2-ME"}),
    MatrixParameter(
        "kg_2_octanone_soy", "per-gram reversible binding constant",
        1.24e-2, "L/g", "2_octanone", "soy_protein", "equilibrium_dialysis",
        8.0, 25.0, "damodaran1981 via k2 sec. (b)", "measured_ratio",
        "FIT (D.6 Module 6)", provenance={"molar_basis": "stated_by_source"}),
    MatrixParameter(
        "kg_2_nonanone_soy", "per-gram reversible binding constant",
        3.72e-2, "L/g", "2_nonanone", "soy_protein", "equilibrium_dialysis",
        8.0, 25.0, "damodaran1981 via k2 sec. (b)", "measured_ratio",
        "FIT (D.6 Module 6)",
        notes="Agrees with Andriot's beta-lactoglobulin value to 1.8x across "
              "three proteins, three methods, three labs, 43 years -- the "
              "strongest cross-validation in the binding batch, and it says the "
              "KETONE per-gram constants are transferable.",
        provenance={"molar_basis": "stated_by_source"}),
    MatrixParameter(
        "kg_nonanal_soy", "per-gram reversible binding constant",
        4.38e-2, "L/g", "nonanal", "soy_protein", "equilibrium_dialysis",
        8.0, 25.0, "damodaran1981 via k2 sec. (b)", "measured_ratio",
        "FIT (D.6 Module 6)",
        notes="DIALYSIS + 2-mercaptoethanol: this constant EXCLUDES the "
              "cysteine-aldehyde chemistry a headspace determination would "
              "count. It must never be pooled with a headspace aldehyde value.",
        provenance={"molar_basis": "stated_by_source",
                    "suppresses_cysteine_aldehyde_chemistry": True}),
    MatrixParameter(
        "kg_hexanal_soy_denatured", "per-gram reversible binding constant",
        1.47e-3, "L/g", "hexanal", "soy_protein", "gel_filtration", None, None,
        "Arai 1970 via damodaran1981, via k2 sec. (b)", "measured_ratio",
        "FIT (D.6 Module 6, via Damodaran)",
        notes="Partially denatured soy. 35x BELOW the headspace-derived lupin "
              "value for the same compound -- the method gap k2 sec. B.3 names.",
        provenance={"basis": "per_gram_no_molar_mass_needed"}),
    # ---- Andriot 2000 beta-lactoglobulin (Module 6 FIT), basis recovered.
    MatrixParameter(
        "kg_2_heptanone_blg", "per-gram reversible binding constant",
        8.97e-3, "L/g", "2_heptanone", "beta_lactoglobulin", "headspace_depletion",
        3.0, 30.0, "andriot2000 via k2 sec. (b) and sec. B.1", "measured_ratio",
        "FIT (D.6 Module 6)",
        provenance={"molar_basis": "recovered_by_arithmetic (36 800 g/mol dimer)"}),
    MatrixParameter(
        "kg_2_octanone_blg", "per-gram reversible binding constant",
        2.58e-2, "L/g", "2_octanone", "beta_lactoglobulin", "headspace_depletion",
        3.0, 30.0, "andriot2000 via k2 sec. (b)", "measured_ratio",
        "FIT (D.6 Module 6)",
        provenance={"molar_basis": "recovered_by_arithmetic"}),
    MatrixParameter(
        "kg_2_nonanone_blg", "per-gram reversible binding constant",
        6.63e-2, "L/g", "2_nonanone", "beta_lactoglobulin", "headspace_depletion",
        3.0, 30.0, "andriot2000 via k2 sec. (b)", "measured_ratio",
        "FIT (D.6 Module 6)",
        provenance={"molar_basis": "recovered_by_arithmetic"}),
)

#: Binding constants that EXIST in the corpus and are deliberately NOT carried
#: here because they are declared HOLD-OUT. The keys are registered so that a
#: caller can see the gap is a governance decision and not an oversight.
HOLDOUT_SEALED_BINDING: Mapping[str, str] = {
    "kg_hexanal_lupin": "Barallat-Perez 2024 lupin -- Module 6 STAR HOLD-OUT (D.6). "
                        "Value not carried in this file.",
    "kg_nonanal_lupin": "Barallat-Perez 2024 lupin -- Module 6 STAR HOLD-OUT.",
    "kg_2_nonanone_lupin": "Barallat-Perez 2024 lupin -- Module 6 STAR HOLD-OUT.",
    "kg_hexanal_mucin": "Barallat-Perez 2024 pig gastric mucin -- Module 6 STAR HOLD-OUT.",
    "kg_nonanal_mucin": "Barallat-Perez 2024 mucin -- Module 6 STAR HOLD-OUT.",
    "kg_2_nonanone_mucin": "Barallat-Perez 2024 mucin -- Module 6 STAR HOLD-OUT.",
}

#: The measured chain-length slope of the per-gram binding constant, per CH2.
#: Two independent proteins, two independent labs, agreeing to 6 %: Andriot
#: beta-lg 2.72x/CH2 and Damodaran soy 2.9x/CH2 (k2 sec. B.6). This is the ONLY
#: licence in this layer to move a constant from one chain length to another,
#: it applies to BINDING CONSTANTS ONLY, and it carries a surrogate flag.
#: NOTE, and it matters: k2 sec. D.3 shows the chain-length structure of the
#: measured THRESHOLD SHIFT does NOT transfer (monotone in gelatin, collapsed in
#: beef). Nothing here extrapolates a threshold shift by chain length.
CHAIN_LENGTH_SLOPE_PER_CH2: float = 2.81  # geometric mean of 2.72 and 2.9
CHAIN_LENGTH_SLOPE_SOURCE = (
    "Andriot 2000 2.72x/CH2 (beta-lactoglobulin, headspace) and Damodaran 1981 "
    "2.9x/CH2 (soy, dialysis), k2_matrix_and_thresholds.md sec. B.6. Both FIT."
)

#: Amendment 6 ruling 2, computed not asserted: reversible binding supplies
#: ~25 % of the hexanal beef/water log-shift and covalent ~0.06 %. This is the
#: CAP on the reversible term's share of any observed log-shift, and the layer
#: refuses to report a reversible explanation above it.
REVERSIBLE_LOG_SHIFT_CEILING: float = 0.25
REVERSIBLE_LOG_SHIFT_CEILING_SOURCE = (
    "Amendment 6 ruling 2 (meynier2004_extraction.md sec. 9): reversible 25.4 % "
    "+ covalent 0.06 % = 25.5 % of the 1 304x hexanal log-shift on a "
    "threshold-panel timescale; a month of storage reaches only ~40 %. "
    "k2_matrix_and_thresholds.md sec. B.4 independently: best-case reversible "
    "explanation 6.2x (log10 0.792) against a measured log10 3.115."
)


# ===========================================================================
# 3. TERM 2 -- THE ALPHA,BETA-UNSATURATION PENALTY
# ===========================================================================
# The only parametric matrix term k2 sec. D.4.4 judged defensible: it
# reproduces in independent matrices, it is chemically motivated (2-alkenals
# are Michael acceptors for lysine/cysteine), and it runs OPPOSITE to
# hydrophobicity (t-2-hexenal has a LOWER log P than hexanal and a larger
# penalty), which is exactly the signature of adduction rather than partition.
#
# FIT ROWS ONLY. The beef observation (2.0x) is EXCLUDED from the fit because
# Brewer 1995 is declared HOLD-OUT and reclassified `dose_added_pre_cook`.

ALPHA_BETA_UNSATURATION_OBSERVATIONS: Tuple[MatrixParameter, ...] = (
    MatrixParameter(
        "unsat_penalty_gelatin", "alkenal / alkanal matrix-shift ratio",
        2.81, "x", "t_2_hexenal / hexanal", "gelatin_3pct", "sensory_BET",
        None, 22.0,
        "vega1994 via k2 sec. A.2: gelatin/water 36.3x (t-2-hexenal) over 12.9x "
        "(hexanal). Equivalently 109/58 ppb in gel times the 4.5/3.0 water ratio.",
        "measured_ratio", "FIT (declaration D.6 Module 7: Vega gelatin ladder)",
        provenance={"thermal_step_after_dosing": False, "concentration_verified": False,
                    "criterion": "75 % detection, UNCORRECTED for chance",
                    "cross_study_cross_method": True,
                    "aqueous_reference_source": "Guadagni 1963/72 via Vega (reference only)"}),
    MatrixParameter(
        "unsat_penalty_dairy_headspace", "alkenal / alkanal matrix-shift ratio",
        4.95, "x", "t_2_hexenal / hexanal", "skim_milk", "static_headspace_partition",
        None, 30.0,
        "Meynier2002_extraction.md sec. 5.1: 6.88x over 1.39x.",
        "measured_ratio", "FIT (Amendment 4: Meynier partition RATIOS)",
        notes="The anhydrous-milk-fat column shows the OPPOSITE ranking with a "
              "mere 1.4x gap (44.0x vs 30.5x): LIPIDS do not discriminate the two "
              "aldehydes, PROTEINS do by 4.9x. That contrast is the cleanest "
              "evidence in the paper that the effect is chemistry, not "
              "hydrophobicity, and it is why the penalty is gated on protein.",
        provenance={"cross_study_cross_method": False,
                    "lipid_control_ratio": 1.44,
                    "quarantine": "measured on a visibly browning system"}),
)

#: The EXCLUDED observation, registered so its exclusion is visible.
UNSATURATION_OBSERVATIONS_EXCLUDED: Mapping[str, str] = {
    "unsat_penalty_beef": "Brewer 1995 gives 2 623x / 1 304x = 2.01x. EXCLUDED: "
                          "Brewer is declaration D.6 HOLD-OUT and reclassified "
                          "`dose_added_pre_cook`, so a large part of both numbers "
                          "is thermal loss before perception, not perception.",
}

#: Ordinal gates the penalty must not contradict (anantharamkrishnan2020b sec. 8).
UNSATURATION_ORDINAL_GATES: Tuple[str, ...] = (
    "G-1 aldehyde reactivity ordering: enals >> hexanal > benzaldehyde ~ furfural "
    "> citral > vanillin (none). A model that ranks any pair the other way is falsified.",
    "G-2 (ordinal form, Amendment 6 ruling 3): enals cross-link the protein at "
    ">=4x lower dose than hexanal; hexanal DOES cross-link at high dose. The "
    "original binary form ('hexanal does not cross-link') is FALSIFIED.",
)


# ===========================================================================
# 4. TERM 3 -- THE COVALENT CEILING (INERT)
# ===========================================================================
# This term contributes EXACTLY ZERO to every point prediction. It exists so
# that the layer can say, in reporting, how large the covalent channel could be
# at a stated temperature and time -- and say honestly that the Ea which decides
# that is unmeasured in every source the corpus holds.

@dataclass(frozen=True)
class CovalentCeiling:
    k2_m_per_s_at_20c: float
    ea_kj_mol: Optional[float]
    reversible_share_headspace_timescale: float
    reversible_share_48h: float
    ambient_half_life_days: Tuple[float, float]
    activation_condition_ea_kj_mol: float
    source_anchor: str
    notes: str
    #: B8: the two fields that turn "unmeasured, could matter" into "measured,
    #: does not". `ea_kj_mol` is now a NUMBER and `ea_status` says so.
    ea_range_kj_mol: Optional[Tuple[float, float]] = None
    ea_status: str = ""
    process_temperature_verdict: str = ""


#: ===================================================================
#: B8 -- THE CEILING IS RETIRED TO "MEASURED NEGLIGIBLE"
#: ===================================================================
#: FIT_HOLDOUT_DECLARATION.md Amendment 17 clause 6, ratified 2026-08-30.
#: Amendment 6 ruling 2 posed the question in binary, quantitative form: this
#: channel matters at process temperature ONLY if its Ea is >= ~70 kJ/mol, and
#: that Ea is unmeasured in every corpus source. WAVE K6b MEASURED IT.
#:
#:   Ea = 15-23 kJ/mol.  The threshold is MISSED BY 3.5-4.7x.
#:
#: Primary: Shepelev & Reineccius 2024, JAFC 72:10579-10583, Fig. 3 -- decanal
#: + beta-lactoglobulin, 14-C scintillation counting of protein-bound label at
#: 25 / 45 / 65 C. Raw day-1 extents give 15.2 kJ/mol (R^2 0.998); the
#: saturation-corrected first-order approach to plateau gives 20.0 (R^2 0.9997);
#: the full envelope over every digitisation choice is 14-23; propagating the
#: printed SDs (n = 4) gives 20 +/- 4 (1 sigma), 95% CI 12-28.
#:
#: FOUR INDEPENDENT CORROBORATIONS, none sharing a method with the primary:
#:   1. Hamzalioglu 2018's HMF + free lysine, aqueous pH 3.5, 5-50 C, HPLC-UV:
#:      10.0 kJ/mol. Different lab, aldehyde, nucleophile presentation, method.
#:   2. Shepelev's OWN PEITC control, same figures and apparatus: 36-43 kJ/mol.
#:      A 2.0-2.5x separation WITHIN one experiment proves the estimator
#:      RESOLVES an Ea -- the low decanal number is not an estimator artefact.
#:   3. Yuan 2023's ordinal panel: hexanal flat at 1,1,1,1 across 22 -> 130 C
#:      in a table that DOES register crossings for compounds that accelerate.
#:   4. Hidalgo 2010's thiol-Michael series, 8 temperatures, 80-180 C: 28-30.
#:
#: WHY THE OLD NUMBER WAS WRONG, in one sentence: every intuition the corpus had
#: about "carbonyl-amine activation energies" was built on BAND C -- Strecker
#: chemistry, browning, fluorescence, all at 50-81 kJ/mol -- and applied to a
#: BAND A reaction. Adduct FORMATION lives at 10-30. The 70 was a band-C number
#: applied to a band-A step.
#:
#: WHAT IT DOES TO THE ARITHMETIC. Acceleration 25 C -> process temperature:
#:   Ea = 20 (measured):  5.1x at 100 C,  9.5x at 140 C,  15.8x at 180 C
#:   Ea = 70 (assumed):   291x,           2 596x,         15 700x
#: Fraction of a hexanal dose consumed covalently, computed from the FAST end
#: of the corpus's own t_1/2 bracket so every figure OVER-states the sink:
#:   UHT 130 C/30 s 0.005% | extrusion 140 C/60 s 0.012% |
#:   retort 100 C/30 min 0.20% | bake 180 C/10 min 0.21%
#: Under the >= 70 assumption those four are 1.0% / 3.3% / 4.6% / 90.7%.
#: THE ASSUMPTION WAS WRONG BY TWO TO THREE ORDERS OF MAGNITUDE IN THE QUANTITY
#: THAT MATTERS. Against the 1 304x (3.115 decades) hexanal cooked-beef/water
#: threshold shift the channel was invoked to explain, 0.21% contributes
#: 0.0009 decades -- 0.03% of it.
#:
#: TWO ARGUMENTS THAT NEED NO ARRHENIUS ARITHMETIC AT ALL, so the verdict does
#: not rest on the Ea value being right:
#:   (i)  THE SINK'S CAPACITY FALLS WITH TEMPERATURE, measured. Day-28 binding
#:        is 25.6 > 20.1 > 17.4 mg/g at 25 / 45 / 65 C, and the 65 C series
#:        PEAKS AT DAY 7 then declines 26%. The authors give the mechanism and
#:        say they favour it: heating alters protein structure so it becomes
#:        less available. Heating does not merely fail to speed the sink up; it
#:        shrinks the sink's ceiling.
#:   (ii) THE EQUILIBRIUM UNWINDS ON HEATING, measured. Zamora 2010's matched
#:        forward/reverse pair -- the corpus's only one -- gives Ea_fwd 44 and
#:        Ea_rev 52, hence dH ~ -8 kJ/mol and K_eq falling 3.0x from 25 to
#:        180 C; demonstrated by acrylamide that was 99.3% gone after 28 d at
#:        60 C being restored to 10.6% by 20 min at 180 C, a 15x recovery of an
#:        analyte that had, to all measurement, disappeared.
#:
#: WHAT THE CHANNEL IS STILL GOOD FOR, stated positively: AMBIENT STORAGE over
#: weeks to months in a high-protein matrix, where it is real and sizeable
#: (decanal loses 14.5-20% of dose over 28 d at 25 C; hexanal's t_1/2 is
#: 18-74 d). That is where the sink belongs, and this model has no storage clock.
#:
#: ⚠ WHAT THE VERDICT DOES NOT SETTLE, carried forward rather than closed:
#: REVERSIBILITY of an aldehyde-protein adduct is STILL never measured. Three
#: papers assert it and none tests it. It is now a well-specified experiment
#: (synthesise the adduct, heat it alone, quantify released aldehyde against a
#: labelled IS) rather than an open-ended gap.
#:
#: ⚠ AND A CONTRADICTION FOUND WHILE DOING THIS, REPORTED NOT RESOLVED: this
#: object carries `ambient_half_life_days = (37, 760)` while the lipid lane's
#: `CovalentSinkCeiling` carries (37, 74). They are two different brackets --
#: 37-760 d is anantharamkrishnan2020b's MS adduct-counting range, 37-74 d is
#: the OVERLAP of that with meynier2004's independent bound. Neither is wrong;
#: they answer different questions. NEITHER IS CHANGED HERE. B8 will not
#: harmonise two measured brackets by picking one after reading a scorecard.
COVALENT_CEILING = CovalentCeiling(
    k2_m_per_s_at_20c=2.5e-5,
    # B8: MEASURED. Was None ("UNMEASURED, everywhere in the corpus").
    ea_kj_mol=20.0,
    ea_range_kj_mol=(15.0, 23.0),
    ea_status=(
        "MEASURED (Shepelev 2024 Fig. 3, 14-C on beta-lactoglobulin, 3 "
        "temperatures): 15-23 kJ/mol, central 20.0. Corroborated at 10.0 by "
        "Hamzalioglu 2018 in a different lab on a different aldehyde by a "
        "different method. The most generous defensible upper bound is 45 "
        "(adding back a hypothetical hydrophobic pre-equilibrium the paper does "
        "not measure and this layer would not want, since the sink acts on "
        "total rather than bound aldehyde) -- and even 45 < 70."
    ),
    process_temperature_verdict=(
        "NEGLIGIBLE AT PROCESS TEMPERATURE, MEASURED. The channel removes "
        "0.005%-0.21% of a hexanal dose in any realistic thermal process "
        "(UHT 130 C/30 s, extrusion 140 C/60 s, retort 100 C/30 min, bake "
        "180 C/10 min), against the 1%-91% the >= 70 kJ/mol assumption implied. "
        "Three independent mechanisms point the same way: the rate barely "
        "rises (Ea 15-23), the CAPACITY FALLS with temperature (day-28 binding "
        "25.6 > 20.1 > 17.4 mg/g at 25/45/65 C), and the EQUILIBRIUM UNWINDS on "
        "heating (Ea_rev 52 > Ea_fwd 44 => K_eq falls 3.0x from 25 to 180 C). "
        "It is an AMBIENT-STORAGE channel -- real and sizeable over weeks in a "
        "high-protein matrix -- and no amount of heat makes it a process "
        "channel. The term stays INERT and contributes exactly 0.0, which it "
        "already did; what changed is that the zero is now MEASURED rather than "
        "ASSUMED. STILL OPEN: adduct reversibility is asserted by three papers "
        "and tested by none."
    ),
    reversible_share_headspace_timescale=0.98,
    reversible_share_48h=0.77,
    ambient_half_life_days=(37.0, 760.0),
    #: KEPT, as the historical decision threshold the measurement is scored
    #: against. It is no longer a condition; it is the number that was missed.
    activation_condition_ea_kj_mol=70.0,
    source_anchor=(
        "Amendment 6 rulings 1, 2 and 4 (meynier2004_extraction.md sec. 9, G-2) "
        "bracketed two-sided against anantharamkrishnan2020b's floor. THE Ea IS "
        "NOW MEASURED: shepelev2024_extraction.md sec. 6 / "
        "k6b_adduct_kinetics_synthesis.md secs. 1b, 2b, 2d, 2e (Amendment 17 "
        "clause 6), with hamzalioglu2018, hidalgo2010 and zamora2010 as the "
        "independent corroborations."),
    notes=(
        "k2(hexanal-Lys) <= 2.5e-5 M^-1 s^-1 at 20 C is a CEILING, measured as "
        "'<= 6 % of a 40 mM dose in 48 h'. The measured analogues are "
        "t-2-hexenal 5.3-7.9e-5 (Lys) and 1.5-2.4e-4 (His). Hexanal binding is "
        ">=98 % reversible on headspace timescales and >=77 % over 48 h; "
        "t-2-hexenal is ~2/3 reversible at 48 h. The covalent channel supplies "
        "~0.06 % of the 1 304x hexanal log-shift on a threshold-panel "
        "timescale, so it does NOT close the matrix gap at ambient. THE ONE "
        "SURVIVING CORNER IS NOW CLOSED TOO (Wave B8): Amendment 6's argument "
        "was that extent is uncapped (lysine is in 360-36 000x excess), so at "
        "process temperature with Ea >= ~70 kJ/mol, 36x would arrive in 36 min "
        "- 2.2 h at 140 C. THE Ea IS MEASURED AT 15-23 kJ/mol, so the "
        "acceleration to 140 C is 9.5x rather than 2 596x and 36x never "
        "arrives. The uncapped-extent argument also fails on its own terms: "
        "the sink's measured CAPACITY FALLS with temperature. The named "
        "wet-lab gap Amendment 6 opened is CLOSED; the one that replaces it is "
        "narrower and different -- nobody has ever tested whether an "
        "aldehyde-protein adduct is reversible."),
)

#: B8: kept as a named constant so the retirement is greppable and so the
#: decision threshold and the measurement that missed it travel together.
COVALENT_CEILING_RETIREMENT = (
    "RETIRED 2026-08-30 (Wave B8, Amendment 17 clause 6). State change: "
    "'UNMEASURED, could matter at process temperature if Ea >= 70 kJ/mol' -> "
    "'MEASURED NEGLIGIBLE, Ea = 15-23 kJ/mol'. NO PREDICTED VALUE MOVES: the "
    "term was already INERT BY RULING and contributed exactly 0.0 to every "
    "point prediction, so what this replaces is an UNMEASURED zero with a "
    "MEASURED one. Do not read a score change into it; if one appears, the "
    "retirement was not the no-op it is declared to be and that is a finding."
)

#: STRUCTURAL GATE 1 -- the 32-of-47 no-adduct negative gate. A compound on
#: this list gets NO covalent term at all, at any temperature. This is the
#: single largest structural constraint in the layer and it falsifies any
#: generic "protein binding removes volatiles" term by 32 counter-examples at a
#: saturating dose (12 ppth, i.e. 10^6-10^7x a food-relevant loading).
NO_ADDUCT_GATE_COMPOUNDS: Tuple[str, ...] = (
    "p-cymene",
    "(l)-menthol", "geraniol", "2-pentanol", "2,3-butanediol",
    "2,6-dimethylphenol", "eugenol", "4-vinylphenol",
    "furaneol", "maltol",
    "butyric acid", "2-hydroxybenzoic acid", "m-toluic acid",
    "methyl anthranilate", "iso-amyl acetate", "ethyl lactate", "methyl salicylate",
    "delta-dodecalactone", "gamma-butyrolactone", "gamma-decalactone",
    "2-heptanone", "2-nonanone", "cyclotene",
    "2,5-dimethylpyrazine", "3-acetylpyridine", "allylamine",
    "vanillin",
    "4-methyl-5-vinylthiazole", "2-methylthiophene", "dimethyl sulfone",
    "dimethyl sulfide", "dimethyl disulfide",
)
assert len(NO_ADDUCT_GATE_COMPOUNDS) == 32, "the 32-of-47 gate must have 32 members"

#: The classes the negative gate covers, which is how it reaches compounds the
#: panel did not itself contain. Every one of these classes is represented on
#: the 32-member list by at least one member, most by two or more.
NO_ADDUCT_GATE_CLASSES: Tuple[str, ...] = (
    "phenol",            # 2,6-dimethylphenol, eugenol, 4-vinylphenol
    "carboxylic_acid",   # butyric acid, 2-hydroxybenzoic acid, m-toluic acid
    "ester",             # methyl anthranilate, iso-amyl acetate, ethyl lactate, methyl salicylate
    "lactone",           # three lactones
    "methyl_ketone",     # 2-heptanone, 2-nonanone
    "alkylpyrazine",     # 2,5-dimethylpyrazine
    "furanone",          # furaneol, maltol
    "alcohol",
    "hydrocarbon",
    "disulfide",         # dimethyl disulfide -- SEE THE CONTRADICTION BELOW
)

#: Classes that DO form adducts, from the same table's positive half.
ADDUCT_POSITIVE_CLASSES: Tuple[str, ...] = (
    "n_alkanal", "branched_alkanal", "alkenal", "thiol", "trisulfide",
    "isothiocyanate", "diketone",
)

#: A CONTRADICTION IN THE SOURCE, carried rather than resolved.
SOURCE_CONTRADICTIONS: Mapping[str, str] = {
    "dimethyl_disulfide_adduct": (
        "anantharamkrishnan2020b Table 2 lists dimethyl disulfide as `no` "
        "(unreactive) and is one of the 32. Its own Results text, pp. 12-13, "
        "says DMDS DOES form a covalent bond with beta-lactoglobulin and names "
        "the adduct (+46 Da, BLG-CysSSMe), while noting the extent is less than "
        "dimethyl trisulfide's. THE TABLE AND THE TEXT CONTRADICT EACH OTHER. "
        "This layer resolves it CONSERVATIVELY -- DMDS gets no covalent term, "
        "matching the tabulated gate -- and reports the contradiction on every "
        "DMDS prediction rather than picking a side silently."),
    "meynier_partition_absolute_scale": (
        "Meynier's own measured air/water K sits 6.24x below the Henry's "
        "constants printed in its own Table I, systematically across five "
        "compounds in two chemical classes with a 1.6x spread. Leksrisompong "
        "lands 6-17x low for the same methodological reason (gas-phase sample "
        "quantified against a liquid-phase calibration). This is NOT a licence "
        "to swap the constant; it is why every absolute here carries +/-0.5 "
        "decades and why only RATIOS are fitted."),
    "unsaturation_penalty_magnitude": (
        "k2 sec. D.3 puts the alpha,beta-unsaturation penalty at ~2-3x from two "
        "sensory matrices; Meynier's headspace contrast is 4.95x. k4b flags the "
        "disagreement explicitly: 'the term's magnitude is matrix-dependent'. "
        "The FIT-only estimate here (3.73x) therefore sits ABOVE k2's stated "
        "2-3x band, and it does so because the one observation that would pull "
        "it down (beef, 2.01x) is HOLD-OUT and excluded."),
}

#: STRUCTURAL GATE 2 -- the pH-3 adduct gate. Carbonyl-lysine adduct formation
#: is ABOLISHED at pH 3 (anantharamkrishnan2020 sec. C1), corroborated in-wave
#: by Leksrisompong's diacetyl/caseinate result (binds at pH 7, N.S. at pH 5.5).
#: Any aldehyde-loss-to-protein term must go to zero at acid pH.
PH_ADDUCT_GATE_BELOW: float = 3.0
PH_ADDUCT_GATE_SOURCE = (
    "anantharamkrishnan2020 sec. C1 (benzaldehyde forms a +88 Da Schiff base at "
    "pH 7-8 and NOTHING at pH 3), declaration Amendment 4 NEITHER row: the paper "
    "is a mechanism_reference with no rates, and its pH-3 adduct gate is a "
    "STRUCTURAL CONSTRAINT ONLY. Corroborated by leksrisompong2010 Table 4/6."
)
#: Between PH_ADDUCT_GATE_BELOW and this value the covalent channel is
#: attenuated but not proven absent; the layer reports `pH_uncertain` rather
#: than interpolating, because no source measures the intermediate.
PH_ADDUCT_GATE_UNCERTAIN_BELOW: float = 5.5


# ===========================================================================
# 5. ABSOLUTE-CONCENTRATION RELIABILITY -- why an absolute is never a point
# ===========================================================================
#: Same-sample HS-SPME dispersion, measured: Xin 2026 and Xin 2026b report the
#: SAME SAMPLES 10-23x apart. The declaration calls this "NOT A DATASET -- a
#: CALIBRATION FACT" and requires it be applied as an `absolute_concentration:
#: false` flag, never scored.
HS_SPME_SAME_SAMPLE_DISPERSION: Tuple[float, float] = (10.0, 23.0)
HS_SPME_DISPERSION_SOURCE = (
    "Declaration D.6 Module 8, 'Xin 2026 vs Xin 2026b same-sample 10-23x "
    "discrepancy -- NOT A DATASET, a CALIBRATION FACT' (sec. E.2.2, sec. B9.1)."
)

#: The air/water partition constant for hexanal spans 9.5x across the whole
#: literature (EPI-Suite 8.7e-3; Meynier measured 1.9e-3; Hall & Anderson
#: 1.8e-2, all at 25-30 C). k4b sec. D.2's ruling is explicit: this is NOT a
#: licence to swap the constant, it is an instruction to carry +/-0.5 decades.
K_AW_UNCERTAINTY_DECADES: float = 0.5
K_AW_UNCERTAINTY_SOURCE = (
    "k4b_paired_thresholds_and_browning.md sec. D.2: total literature spread "
    "9.5x on hexanal K_aw; two independent static-headspace papers land 6-17x "
    "below their own cited references for the same methodological reason."
)

#: Hexanal K_aw in water, the measured surface, dimensionless air/matrix.
#: Guard from k4b hold-out #11: a repo hexanal partition term at 30 C outside
#: 1.9e-3 to 1.8e-2 is outside the ENTIRE measured literature.
K_AW_HEXANAL_WATER_GUARD: Tuple[float, float] = (1.9e-3, 1.8e-2)


# ===========================================================================
# 6. POLICY ASSERTIONS
# ===========================================================================

def assert_no_dft_matrix() -> None:
    """No barrier, constant or bound in this file is computational."""
    for parameter in REVERSIBLE_BINDING + ALPHA_BETA_UNSATURATION_OBSERVATIONS:
        if parameter.evidence_class not in MATRIX_EVIDENCE_CLASSES:
            raise AssertionError(f"{parameter.key}: illegal evidence class")
        if "dft" in parameter.source_anchor.lower():
            raise AssertionError(f"{parameter.key}: DFT provenance is forbidden")


def assert_no_mocked_protein_source() -> None:
    """
    The OLD lane's ``protein_source_registry.json`` is declared
    ``no_verifiable_source`` and raises RuntimeWarnings where it is used. This
    layer must not import it, reproduce it, or carry any protein
    differentiation traceable to it. Every protein loading here comes from a
    named printed composition or from an explicitly DECLARED ASSUMPTION.
    """
    import sys
    forbidden = ("src.matrix_correction", "src.headspace")
    for name in forbidden:
        if name in sys.modules and name in globals().get("__imported_by_b4__", ()):
            raise AssertionError(f"B4 must not import the old lane module {name}")
    for loading in MATRIX_LOADING.values():
        if loading.evidence_class not in ("measured_ratio", "declared_assumption"):
            raise AssertionError(
                f"{loading.key}: protein loading must be printed-composition or a "
                f"declared assumption, never a registry-mocked differentiation")


def binding_constant_for(
    compound: str,
    method_family: str = "headspace",
) -> Optional[MatrixParameter]:
    """
    The per-gram reversible binding constant for ``compound``, or None.

    ``method_family`` is 'headspace' or 'dialysis'. FOR ALDEHYDES THE TWO ARE
    NEVER POOLED (k2 sec. B.3: a 35x gap on an aldehyde where the ketone gap is
    1.3x, with the mechanism printed by the source). For non-aldehydes the
    constraint is not asserted by any source, so both families are searched.
    """
    structure = COMPOUND_STRUCTURE.get(compound)
    is_aldehyde = structure is not None and structure.binding_class in (
        "n_alkanal", "branched_alkanal", "alkenal")
    headspace_methods = ("headspace_depletion", "static_headspace_partition")
    dialysis_methods = ("equilibrium_dialysis", "gel_filtration")
    wanted = headspace_methods if method_family == "headspace" else dialysis_methods
    for parameter in REVERSIBLE_BINDING:
        if parameter.compound != compound:
            continue
        if is_aldehyde and parameter.method not in wanted:
            continue
        if not is_aldehyde and parameter.method not in (wanted + (
                dialysis_methods if method_family == "headspace" else headspace_methods)):
            continue
        return parameter
    return None


def matrix_registry_metadata() -> Dict[str, object]:
    """Everything a report needs to reproduce and defend this layer."""
    return {
        "evidence_classes": list(MATRIX_EVIDENCE_CLASSES),
        "n_binding_constants": len(REVERSIBLE_BINDING),
        "binding_constants": [
            {
                "key": p.key, "compound": p.compound, "medium": p.medium,
                "value_l_per_g": p.value, "method": p.method,
                "pH": p.ph_of_measurement, "temperature_c": p.temperature_c,
                "source": p.source_anchor, "role": p.declaration_role,
                "evidence_class": p.evidence_class,
                "provenance": dict(p.provenance), "notes": p.notes,
            }
            for p in REVERSIBLE_BINDING
        ],
        "holdout_sealed_binding": dict(HOLDOUT_SEALED_BINDING),
        "chain_length_slope_per_ch2": CHAIN_LENGTH_SLOPE_PER_CH2,
        "chain_length_slope_source": CHAIN_LENGTH_SLOPE_SOURCE,
        "reversible_log_shift_ceiling": REVERSIBLE_LOG_SHIFT_CEILING,
        "reversible_log_shift_ceiling_source": REVERSIBLE_LOG_SHIFT_CEILING_SOURCE,
        "unsaturation_observations": [
            {"key": p.key, "value_x": p.value, "medium": p.medium,
             "source": p.source_anchor, "role": p.declaration_role}
            for p in ALPHA_BETA_UNSATURATION_OBSERVATIONS
        ],
        "unsaturation_observations_excluded": dict(UNSATURATION_OBSERVATIONS_EXCLUDED),
        "unsaturation_ordinal_gates": list(UNSATURATION_ORDINAL_GATES),
        "covalent_ceiling": {
            "k2_M_per_s_at_20C": COVALENT_CEILING.k2_m_per_s_at_20c,
            "ea_kJ_per_mol": COVALENT_CEILING.ea_kj_mol,
            "ea_range_kJ_per_mol": list(COVALENT_CEILING.ea_range_kj_mol or ()),
            "ea_status": COVALENT_CEILING.ea_status,
            "process_temperature_verdict": (
                COVALENT_CEILING.process_temperature_verdict),
            "retirement": COVALENT_CEILING_RETIREMENT,
            "reversible_share_headspace": COVALENT_CEILING.reversible_share_headspace_timescale,
            "reversible_share_48h": COVALENT_CEILING.reversible_share_48h,
            "ambient_half_life_days": list(COVALENT_CEILING.ambient_half_life_days),
            "decision_threshold_ea_kj_mol": (
                COVALENT_CEILING.activation_condition_ea_kj_mol),
            "decision_threshold_missed_by_x": round(
                COVALENT_CEILING.activation_condition_ea_kj_mol
                / COVALENT_CEILING.ea_kj_mol, 2),
            # kept under its old key so no downstream reader breaks; the key
            # now names a threshold that was MISSED rather than a live condition
            "activates_if_ea_at_least": COVALENT_CEILING.activation_condition_ea_kj_mol,
            "contribution_to_point_prediction": 0.0,
            "source": COVALENT_CEILING.source_anchor,
            "notes": COVALENT_CEILING.notes,
        },
        "no_adduct_gate": {
            "n_compounds": len(NO_ADDUCT_GATE_COMPOUNDS),
            "compounds": list(NO_ADDUCT_GATE_COMPOUNDS),
            "classes": list(NO_ADDUCT_GATE_CLASSES),
            "source": "anantharamkrishnan2020b Table 2, 32 of 47 unreactive at "
                      "24-48 h at a saturating 12 ppth dose (Amendment 5).",
        },
        "adduct_positive_classes": list(ADDUCT_POSITIVE_CLASSES),
        "ph_adduct_gate": {
            "covalent_zero_at_or_below_pH": PH_ADDUCT_GATE_BELOW,
            "uncertain_below_pH": PH_ADDUCT_GATE_UNCERTAIN_BELOW,
            "source": PH_ADDUCT_GATE_SOURCE,
        },
        "absolute_reliability": {
            "hs_spme_same_sample_dispersion_x": list(HS_SPME_SAME_SAMPLE_DISPERSION),
            "hs_spme_source": HS_SPME_DISPERSION_SOURCE,
            "k_aw_uncertainty_decades": K_AW_UNCERTAINTY_DECADES,
            "k_aw_source": K_AW_UNCERTAINTY_SOURCE,
            "k_aw_hexanal_water_guard": list(K_AW_HEXANAL_WATER_GUARD),
        },
        "matrix_loadings": {
            k: {"protein_g_per_L": v.protein_g_per_l,
                "protein_band_g_per_L": [v.protein_lo_g_per_l, v.protein_hi_g_per_l],
                "fat_volume_fraction": v.fat_volume_fraction, "pH": v.ph,
                "source": v.source_anchor, "evidence_class": v.evidence_class,
                "notes": v.notes}
            for k, v in MATRIX_LOADING.items()
        },
        "source_contradictions": dict(SOURCE_CONTRADICTIONS),
        "refused_by_policy": [
            "a general matrix correction factor (k2 sec. D.1: a uniform 33x "
            "misplaces the extremes by 10x and 28x in opposite directions)",
            "any partition-derived threshold (Amendment 6; refuted 3 ways on "
            "matched samples)",
            "any term monotone in log P (k4b hold-out guard #4; Hong r = -0.05)",
            "protein_source_registry.json's mocked differentiation "
            "(no_verifiable_source)",
            "Barallat-Perez lupin and mucin constants (Module 6 HOLD-OUT)",
            "Leksrisompong's 24 BETs as table entries (Amendment 4 makes only "
            "its K RATIOS fit-eligible)",
            "Hong 2020's soy thresholds as table entries (GATING HOLD-OUT)",
            "any DFT-derived barrier (standing owner policy)",
        ],
    }


assert_no_dft_matrix()
assert_no_mocked_protein_source()
