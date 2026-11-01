"""The reviewed allowlist of CONSTANTS the engine reads at call time (A3 phase 2, 2026-11-01).

tests/unit/test_engine_config.py requires every public UPPERCASE module-level global that a function in
src/kinetic_core reads (by name in its own module, or as ``module.NAME`` through an import of an engine module)
to be either a ConfigKey binding (src/kinetic_core/config.py) or listed here, by category, with the reason it is
not a switch. A new global read at call time fails until it is keyed or listed; a listed name that is no longer
read fails too, so the list cannot rot. (The phase 1 review guessed ~60 such constants; the scan finds these.)

How it was drafted: by the scan, each name sorted by its value at import. Every number, record and fact was
read against its definition: none is set by a script or a test (checked separately, by the assignment scan),
and none selects a structure or a solver path. The one call-time read that was a solver choice,
``engine.RAMP_STEP_K``, became the key ``engine.ramp_step_k``. Tables, provenance strings, labels and paths were
reviewed by category.
"""

from __future__ import annotations

#: category -> why its names are constants, not switches
WHY = {
    'physical': (
        'a physical constant, reference temperature or unit identity'
    ),
    'numeric': (
        'a measured or declared number (a literature value, a reference condition, a declared threshold, or a '
        'fitting / scoring tolerance) that no job varies'
    ),
    'record': (
        'a measured record (a titratable group, a literature parameter, a declared default buffer or '
        'contract)'
    ),
    'table': (
        'a table: species / reaction / parameter tables, key lists, and the arrays and counts built from them '
        'at import'
    ),
    'provenance': (
        'provenance text (a source, anchor, caveat or declared reason), read into reports'
    ),
    'label': (
        'a label: an enum-like string compared or reported, never chosen per job'
    ),
    'path': (
        'a data path (tests redirect these to tmp_path)'
    ),
    'fact': (
        'a fact about the environment (whether the numba kernel compiled), not a choice'
    ),
}

#: category -> module ("" is the package) -> the names, space-separated
CONSTANTS = {
    'physical': {
        'acrylamide': 'T_REF_A_K',
        'engine': 'CELSIUS',
        'matrix_sites': 'CELSIUS R_KJ',
        'observations': 'LABEL_TEMPERATURE_K',
        'parameters': 'R_KJ T_REF_K',
        'parameters_acrylamide': 'R_KJ T_REF_A_K',
        'parameters_furanic': 'M_INV_DAY_TO_L_PER_MMOL_MIN R_KJ T_REF_K',
        'parameters_lipid': 'MEQ_O2_PER_MMOL_LOOH',
        'parameters_sulfur': 'PKA_REFERENCE_T_K R_KJ T_REF_S_K',
        'ph_state': 'LABEL_TEMPERATURE_K PKA_REFERENCE_T_K R_KJ',
        'sulfur': 'LABEL_TEMPERATURE_K T_REF_S_K',
        'uncertainty': 'R_KJ_PER_MOL_K',
        'vessel': 'DISSOLVED_O2_MM_AT_20C O2_MMOL_PER_ML_HEADSPACE',
    },
    'numeric': {
        'acrylamide_conditions': 'AW_MULTIPLIER PH_EXPONENT_ELIMINATION PH_EXPONENT_FORMATION REFERENCE_PH',
        'directional': 'MEASUREMENT_NOISE_FOLD ZERO_FLOOR_UG_PER_L',
        'engine': 'NETWORK_PH OX_AMBIENT_MMOL_L OX_RESERVOIR_DEFAULT_UNITS OX_SAT_MMOL_L',
        'fitting': 'BOUND_FRACTION COST_AGREEMENT DIAGNOSTIC_MEAN_SE DIAGNOSTIC_TREND_P NULL_WEIGHT_TOL',
        'furanic': (
            'BLANK_PH BLANK_STATED_MAX_SD_FRACTION BLANK_SUGAR_LOADING_MMOL_PER_L BLANK_TEMPERATURE_C '
            'BLANK_TIME_MIN'
        ),
        'lipid': 'OLEATE_MOLAR_ANCHOR_CENTRE Y_HEXANAL_PER_LOOH',
        'matrix_oav': (
            'CHAIN_LENGTH_SLOPE_PER_CH2 K_AW_UNCERTAINTY_DECADES PH_ADDUCT_GATE_BELOW '
            'PH_ADDUCT_GATE_UNCERTAIN_BELOW REVERSIBLE_LOG_SHIFT_CEILING'
        ),
        'parameters': 'AW_OF_MEASUREMENT NETWORK_PH',
        'parameters_acrylamide': 'EA_E2_MEASURED_KJ_MOL FITTED_AW_OF_DECLARATION',
        'parameters_dicarbonyl': 'EA_AMA_G_KJ_MOL EA_FED_3DEOXY_KJ_MOL EA_G_GO_AQUEOUS_KJ_MOL',
        'parameters_furanic': 'BLANK_AMINE_LOADING_MMOL_PER_L BLANK_PH BLANK_TEMPERATURE_C FURANONE_PARTITION_EA_KJ_MOL',
        'parameters_glycation': (
            'AW_OF_MEASUREMENT EA_CML_LOSS_KJ_MOL EA_FLP_CEL_KJ_MOL EA_FLP_CML_KJ_MOL EA_FLP_DECAY_KJ_MOL '
            'EA_GLYC_KJ_MOL GLYCATION_FIT_PH'
        ),
        'parameters_matrix': (
            'CHAIN_LENGTH_SLOPE_PER_CH2 K_AW_UNCERTAINTY_DECADES PH_ADDUCT_GATE_BELOW '
            'PH_ADDUCT_GATE_UNCERTAIN_BELOW REVERSIBLE_LOG_SHIFT_CEILING'
        ),
        'parameters_methionine': 'AW_OF_MEASUREMENT EA_MARP_KJ_MOL EA_MSH_DMDS_KJ_MOL METHIONINE_FIT_PH',
        'parameters_proline': 'AW_OF_MEASUREMENT EA_HA_ATHP_KJ_MOL EA_MGO_PRO_KJ_MOL EA_PYRL_AP_KJ_MOL PROLINE_FIT_PH',
        'parameters_pyrazine': 'AW_OF_MEASUREMENT PYRAZINE_FIT_PH PYRAZINE_PH_KNOT PYRAZINE_REFERENCE_PH',
        'parameters_sulfur': (
            'KANG_EA_FREE_CYS_DEPLETION_KJ_MOL KANG_TTCA_FREE_CYS_YIELD_CEILING_MOL_PCT '
            'K_OX_SUPPLY_PER_UNIT_MIN MELE_SITE_DEFAULT_EA_KJ_MOL PH_NORMALISATION '
            'PROTEIN_DISULFIDE_EXCHANGE_USED_L_PER_MMOL_MIN PROTEIN_DISULFIDE_SITES_BLG_10PCT_MMOL_L '
            'STACK_DELTA_H_ADDUCT_KJ_MOL STACK_MEASUREMENT_T_K THIAMINE_H2S_SHARE_FLOOR_PH '
            'THIOL_ADDITION_DEFAULT_EA_KJ_MOL THIOL_QUINONE_EQUILIBRIUM_USED_L_PER_MMOL '
            'ZHANG_EA_THIOL_TO_DISULFIDE_KJ_MOL ZHENG_MECHANISM_BOUNDARY_PH'
        ),
        'ph_state': 'ARP_CARBOXYL_DH ARP_CARBOXYL_PKA',
        'scoring': 'PASS_BAND_LEVEL',
        'thermal_history': 'DEFAULT_START_C',
        'trunk_conditions': 'AW_TABLE_FLOOR PH_EXPONENT_DECADES_PER_UNIT REFERENCE_AW REFERENCE_PH',
        'uncertainty': (
            'B1_T_REF_C B3_T_REF_C B8_T_REF_C CI_LEVEL_PCT MIN_EVALUABLE_CI_WIDTH_LOG10 '
            'PREFACTOR_PRIOR_DECADES'
        ),
        'user_fit': 'FOLD_PASS IDENTIFIABILITY_FLOOR IDENTIFIABILITY_TOL SIGMA_DEFAULT_LOG10 SIGMA_FLOOR_LOG10',
        'user_scoring': 'PASS_BAND_LEVEL',
    },
    'record': {
        'engine': 'DEFAULT_BUFFER',
        'matrix_oav': 'COVALENT_CEILING',
        'parameters_lipid': 'K_HEXANAL_SCHROEN K_LOOH_DECOMP_ANCHOR PROPANAL_HEXANAL_RATIO',
        'parameters_matrix': 'COVALENT_CEILING',
        'parameters_pyrazine': 'K_COND',
        'ph_state': (
            'ACETIC CITRATE CYSTEINE_B9 CYSTEINE_GKL DFG_GROUP FORMIC GLYCINE PHOSPHATE PHTHALATE '
            'PYROPHOSPHATE TTCA_GROUP'
        ),
        'scoring': 'DEFAULT_VALIDATION_CONTRACT',
        'sulfur': 'DEFAULT_BUFFER',
    },
    'table': {
        '': 'FURANIC_PARAMETERS MARTINS_M4 MEASURED_ACRYLAMIDE MEASURED_SULFUR',
        'acrylamide': (
            'ACRYLAMIDE_BY_KEY ACRYLAMIDE_CENTRE_LEDGER ACRYLAMIDE_INDEX ACRYLAMIDE_REACTIONS '
            'ACRYLAMIDE_SINK_REACTIONS ACRYLAMIDE_SOURCE_REACTIONS ACRYLAMIDE_STATE_KEYS '
            'ACRYLAMIDE_STOICHIOMETRY ASN_SCHIFF_AMADORI_SPLIT BALANCED_ELEMENTS FULL_ACRYLAMIDE_REACTIONS '
            'MEASURED_ACRYLAMIDE N_ACRYLAMIDE_STATE OUT_OF_SCOPE SCHIFF_AMADORI_SPLIT TRUNK_REACTIONS'
        ),
        'acrylamide_conditions': (
            'AW_ELIMINATION_SCALE_BAND AW_ELIMINATION_TABLE AW_SCALE_BAND AW_WINDOW ELIMINATION_KEYS '
            'PH_ELIMINATION_KEYS PH_EXPONENT_ELIMINATION_BAND PH_EXPONENT_FORMATION_BAND PH_FORMATION_KEYS '
            'PH_WINDOW SCALED_KEYS'
        ),
        'config': 'BUILD_KEYS KEYS KEY_BY_NAME',
        'ea_class_prior': 'CHAN_1994_TABLE_I_KJ_MOL CLASS_OF_ROUTE MEASURED_EVIDENCE',
        'engine': (
            'ACRYLAMIDE_INDEX AW_TERM_LANES AW_TERM_WINDOWS FURANONE_BANDED_KEYS GLYCATION_TARGET_KEYS LANES '
            'LANE_DEFAULT_TARGETS LIPID_CARRIER_ALIASES MAILLARD_LANES MEASURED_ACRYLAMIDE MEASURED_SULFUR '
            'MOLECULAR_WEIGHT_G_PER_MOL NO_PH_TERM_LANES PH_TERM_WINDOWS PRECURSOR_ALIASES SPECIES_KEYS '
            'SULFUR_INDEX TARGET_ALIASES TRUNK_ARMS UNREPRESENTED_COMPOUNDS'
        ),
        'fit_rows': 'LEDGER_CLASS_TO_METHOD OBSERVABLES POT_FIELDS ROW_FIELDS UNITS',
        'fitting': 'EXTRA_STEPS EXTRA_THRESHOLDS LEAST_SQUARES_SETTINGS SCALE_TAU STEPS THRESHOLDS',
        'furanic': (
            'DECLARED_GAPS EDGE_C_ZERO_BY_DECLARATION FURANIC_PARAMETERS FURANIC_REACTION_KEYS '
            'FURANONE_EA_ASSUMPTION FURANONE_STRUCTURAL_CONSTRAINTS HMF_LIMB_SOURCES '
            'HMF_SOURCE_TOPOLOGY_CORROBORATION NAMING_TRAPS PROHIBITED_DERIVATIONS TRUNK_REACTIONS'
        ),
        'integrate': 'INDEX N_SPECIES REACTIONS SPECIES_KEYS',
        'keyspaces': 'B4_COMPOUND_KEY NO_B4_RECORD',
        'lipid': (
            'CLEAVAGE_MECHANISM FRANKEL1981_OLEATE_SLATE FRANKEL_SLATE FRANKEL_SYSTEM_GEOMETRY '
            'FRANKEL_ZERO_ADDITIVE LIPID_CARRIERS LIPID_KEYS LIPID_OWNED_KEYS LIPID_SPECIES LOOH_POOLS '
            'MOLECULAR_WEIGHT_G_PER_MOL PENTYLFURAN_PER_HEXANAL POSITION_PRODUCTS'
        ),
        'matrix_oav': (
            'ADDUCT_POSITIVE_CLASSES ALPHA_BETA_UNSATURATION_OBSERVATIONS COMPOUND_STRUCTURE '
            'HS_SPME_SAME_SAMPLE_DISPERSION MATRIX_LOADING MATRIX_THRESHOLDS NO_ADDUCT_GATE_CLASSES '
            'NO_ADDUCT_GATE_COMPOUNDS REVERSIBLE_BINDING SEALED_OR_REFUSED_MATRICES '
            'UNIMPLEMENTED_CANDIDATE_TERMS WATER_THRESHOLDS'
        ),
        'matrix_sites': 'BINDING_CLASSES BINDING_OF_SPECIES',
        'network': (
            'BALANCED_ELEMENTS BY_KEY INDEX MARTINS_M4 N_SPECIES REACTIONS SCHIFF_AMADORI_SPLIT SPECIES_KEYS '
            'STOICHIOMETRY TRUNK_CENTRE_LEDGER TRUNK_REACTIONS'
        ),
        'panel': 'RATIO_UNIT_FACTORS UNSCORED_EVIDENCE_CLASSES',
        'parameters': 'CROSS_LAB_COMPARATORS EVIDENCE_CLASSES MARTINS_M4 SCHIFF_AMADORI_SPLIT',
        'parameters_acrylamide': (
            'ASN_SCHIFF_AMADORI_SPLIT CROSS_LAB_CONFLICTS DELIBERATE_UNDERFITS EVIDENCE_CLASSES '
            'HOLDOUT_EXPOSURE_DISCLOSURE MEASURED_ACRYLAMIDE REFUSED_PARAMETERS UNIT_BY_ORDER'
        ),
        'parameters_dicarbonyl': 'DICARBONYL_PARAMETERS DISPUTED_SINK_KEYS FED_3DEOXY_COORDINATES HPD_SINK_KEYS',
        'parameters_furanic': (
            'EDGE_C_ZERO_BY_DECLARATION FITTED_FURANIC FURANIC_PARAMETERS FURANONE_EA_ASSUMPTION '
            'MEASURED_FURANIC PROHIBITED_DERIVATIONS TEMPERATURE_CAP_K'
        ),
        'parameters_lipid': (
            'FRANKEL_ZERO_ADDITIVE LIPID_CARRIERS LIPID_EVIDENCE_CLASSES LIPID_SOURCE_CONTRADICTIONS '
            'PROHIBITED_DERIVATIONS SHIPPED_VALUES_REFUTED'
        ),
        'parameters_matrix': (
            'ADDUCT_POSITIVE_CLASSES ALPHA_BETA_UNSATURATION_OBSERVATIONS COMPOUND_STRUCTURE '
            'HOLDOUT_SEALED_BINDING HS_SPME_SAME_SAMPLE_DISPERSION K_AW_HEXANAL_WATER_GUARD '
            'MATRIX_EVIDENCE_CLASSES MATRIX_LOADING NO_ADDUCT_GATE_CLASSES NO_ADDUCT_GATE_COMPOUNDS '
            'REVERSIBLE_BINDING SOURCE_CONTRADICTIONS UNSATURATION_OBSERVATIONS_EXCLUDED '
            'UNSATURATION_ORDINAL_GATES'
        ),
        'parameters_methionine': 'FROZEN_B18',
        'parameters_pyrazine': 'PYRAZINE_PH_SLOPES',
        'parameters_sulfur': (
            'ALKALINE_PRIORS CONSUMPTION_KEYS CYSTINE_SULFIDE_EXCHANGE_DEFERRED DECAY_EA_BOUNDS '
            'DECAY_EA_BOUNDS_BY_FAMILY DECAY_FAMILIES DELTA_H_IONISATION_KJ_MOL EVIDENCE_CLASSES '
            'FORMATION_EA_PRIOR_SOURCE FORMATION_ROUTE_OF GLUTATHIONE_ARRHENIUS_REFIT '
            'KANG_CYS_CONVERSION_AT_120_MIN KUMAZAWA_FFT_RESIDUAL_FRACTION MEASURED_EA_OVERRIDES '
            'MEASURED_SULFUR MELANOIDIN_SITE_DENSITY_MMOL_PER_G NO_EA_KEYS NO_MEASURED_RATE '
            'PH_TERM_PROVENANCE PKA_25C PKA_GKL_CYSTEINE PROHIBITED_DERIVATIONS '
            'PROTEIN_DISULFIDE_EXCHANGE_BRACKET_L_PER_MMOL_MIN QUINONE_AMINE_K2_CONTEXT_M_INV_S_INV '
            'QUINONE_THIOL_K2_CONTEXT_M_INV_S_INV ROUTE_CLASS_MEDIAN_EA_KJ_MOL THIAMINE_H2S_SHARE_POINTS '
            'THIOL_CHANNELS THIOL_QUINONE_EQUILIBRIUM_M_INV UNIT_BY_ORDER ZHENG_CYSTEINE_THERMOLYSIS '
            'ZHENG_GLUTATHIONE_TABLE_III ZHENG_TABLE_II_PH_LAW'
        ),
        'ph_state': (
            'ACID_YIELD_BOUNDS ARP_AMINE_PKA_BOUNDS CENTRE_KINDS ORGANIC_ACID_MIX PKW_TABLE_C '
            'TITRATABLE_CENTRES'
        ),
        'scoring': 'RATIO_UNIT_FACTORS SHARED_WITH_HOLDOUT_PANEL',
        'species': 'INDEX N_SPECIES SPECIES SPECIES_KEYS',
        'species_acrylamide': (
            'ACRYLAMIDE_INDEX ACRYLAMIDE_MW_G_PER_MOL ACRYLAMIDE_STATE ACRYLAMIDE_STATE_KEYS '
            'N_ACRYLAMIDE_STATE'
        ),
        'species_lipid': 'LIPID_KEYS LIPID_SPECIES MOLECULAR_WEIGHT_G_PER_MOL',
        'species_sulfur': (
            'MOLECULAR_WEIGHT_G_PER_MOL NO_THRESHOLD_AVAILABLE N_SULFUR_STATE ODOUR_THRESHOLD_UG_PER_L '
            'SULFUR_INDEX SULFUR_SPECIES SULFUR_STATE SULFUR_STATE_KEYS'
        ),
        'sulfur': (
            'BALANCED_ELEMENTS CENTRE_LEDGER FULL_REACTIONS MEASURED_SULFUR N_SULFUR_STATE OUT_OF_SCOPE '
            'REACTION_PH_FACTOR SITE_POOLS SULFUR_BY_KEY SULFUR_INDEX SULFUR_REACTIONS SULFUR_STATE_KEYS '
            'SULFUR_STOICHIOMETRY TERMINAL_POOLS TRUNK_REACTIONS'
        ),
        'thermal_history': 'TIME_CONVENTIONS',
        'trunk_conditions': (
            'ATMOSPHERE_FACTORS ATMOSPHERE_VALUES AW_MULTIPLIER_TABLE AW_SCALE_BAND AW_STEPS '
            'OXIDATIVE_ENTRY_STEPS PH_EXPONENT_BAND PH_MEASURED_WINDOW PH_STEPS'
        ),
        'uncertainty': 'CORE_PRIORS DEFINITIONAL_BAND_KINDS RATIO_UNIT_FACTORS SHARED_WITH_HOLDOUT_PANEL',
        'vessel': 'THIOL_PRECURSOR_MARKERS',
    },
    'provenance': {
        'acrylamide_conditions': 'AW_SOURCE PH_SOURCE',
        'engine': 'BUFFER_ABSENT_WARNING',
        'matrix_oav': 'COVALENT_CEILING_RETIREMENT',
        'parameters_furanic': 'MU_M_HAZARD',
        'parameters_lipid': 'DECADIENAL_TO_HEXANAL_GAP FRANKEL_ANCHOR FRANKEL_DOSSIER LOOH_FORMATION_GAP',
        'parameters_matrix': (
            'CHAIN_LENGTH_SLOPE_SOURCE COVALENT_CEILING_RETIREMENT HS_SPME_DISPERSION_SOURCE '
            'K_AW_UNCERTAINTY_SOURCE PH_ADDUCT_GATE_SOURCE REVERSIBLE_LOG_SHIFT_CEILING_SOURCE'
        ),
        'parameters_methionine': 'DENG_SOURCE EA_MARP_DECLARED',
        'parameters_sulfur': (
            'KANG_CYS_ANCHOR KANG_CYS_CAVEATS KANG_TTCA_CEILING_ANCHOR KUMAZAWA_ANCHOR KUMAZAWA_CAVEATS '
            'PROTEIN_DISULFIDE_DERIVATION STACK_TRANSFER_CAVEATS THIOL_PKA_CAVEAT '
            'THIOL_QUINONE_EQUILIBRIUM_ANCHOR ZHANG_DISULFIDE_ANCHOR'
        ),
        'sulfur': 'BUFFER_ABSENT_WARNING',
        'trunk_conditions': 'ATMOSPHERE_TRANSFER AW_SOURCE PH_SOURCE',
        'uncertainty': 'B18_UNSAMPLED_REASON B55_REASON',
    },
    'label': {
        'directional': 'INDEPENDENT NOT_EVALUABLE',
        'ea_class_prior': 'CARBOHYDRATE SULFUR_CLASS',
        'engine': 'ACRYLAMIDE B1_VARIANT HEXOSE_ENTRY_UNIDENTIFIED LIPID SULFUR TRUNK',
        'lipid': 'REFERENCE_AUTOXIDATION_SYSTEM REFERENCE_LINOLEATE_1981_COLUMN REFERENCE_OLEATE_1981_COLUMN',
        'panel': 'PANEL_EXTERNAL_MATRIX PANEL_MAILLARD_PATH_HOLDOUT PANEL_TRUST_LOOP SYNTHETIC_SNAPSHOT_SUFFIX',
        'trunk_conditions': 'ATMOSPHERE_REFERENCE',
        'uncertainty': (
            'ACRYLAMIDE DEFINITIONAL_BAND EMPIRICAL_CLASS FROZEN_NOT_FREE LAPLACE_SAMPLED LIPID '
            'NO_UNCERTAINTY PRINTED_HPD_BAND QUANTIFICATION_EXTRACTION SECOND_LABORATORY_BAND SULFUR TRUNK '
            'UNIDENTIFIED_CAPPED UNIDENTIFIED_FLAT UNIDENTIFIED_NO_BAND UNIDENTIFIED_SAMPLED'
        ),
    },
    'path': {
        'calibration': 'USER_RESULTS_DIR',
        'directional': 'PANEL_PATH',
        'fit_rows': 'FIT_ROWS_DIR LABORATORIES',
        'fit_targets': 'RECORD_GLOB',
        'uncertainty': 'DISCREPANCY_PATH LAPLACE_PATH PROFILE_PATH',
        'user_scoring': 'USER_RESULTS_DIR',
    },
    'fact': {
        'ph_kernel': 'COMPILED',
    },
}


def allowlisted():
    """{(module, name): category} for every listed constant."""
    return {(module, name): category for category, by_module in CONSTANTS.items()
            for module, names in by_module.items() for name in names.split()}
