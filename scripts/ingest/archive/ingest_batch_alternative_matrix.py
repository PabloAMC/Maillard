"""
Batch ingestion: alternative_protein_matrix_scope (27 refs) +
                 amino_acid_sugar_core (5 refs) +
                 glutathione_and_peptide_support (7 remaining) +
                 nucleotide_and_ribose_support (7 remaining) +
                 thiamine_fragmentation_support (5 remaining) +
                 fermentation_pretreatment (1 ref)

Updates:
  - benchmark_intake_registry.json   → status = "encoded"
  - slr_incorporation_matrix.json    → adds entries
  - computational_priors.json        → carbonyl_donor_priors, sulfur_peptide_priors,
                                       nucleotide_pathway_priors, thiamine_pathway_priors,
                                       shear_damage_priors, volatile_class_profiles extensions
  - flavor_reference_payloads.json   → AEDA / off-note flavor anchors
  - process_state_calibrations.json  → fermentation / protein-prep state entries
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_LIT_DIR = ROOT / "data" / "lit"

INTAKE_REGISTRY_PATH = DATA_LIT_DIR / "benchmark_intake_registry.json"
SLR_MATRIX_PATH = DATA_LIT_DIR / "slr_incorporation_matrix.json"
COMPUTATIONAL_PRIORS_PATH = DATA_LIT_DIR / "computational_priors.json"
FLAVOR_PAYLOADS_PATH = DATA_LIT_DIR / "flavor_reference_payloads.json"
PROCESS_STATE_CALIBRATIONS_PATH = DATA_LIT_DIR / "process_state_calibrations.json"


def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def save_json(data, path):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"Updated {path.name}")


# ─── Target IDs ───────────────────────────────────────────────────────────────
TARGET_IDS = {
    # alternative_protein_matrix_scope
    "uspto_ptacts_2023_yeast_extract_anchor",
    "liu_cadwallader_2023_pea_aeda",
    "malia_2025_pea_free_sh_crosscheck",
    "lagrain_2010_cystine_elimination_lanthionine",
    "morel_2002_gluten_shear_aggregation",
    "rombouts_2012_gluten_crosslinking",
    "ilo_1996_maize_sme_lysine_damage",
    "researchgate_2023_pea_aeda",
    "zhang_1993_protein_deamidation_ammonia",
    "li_2010_phytate_chelation_kinetics",
    "zha_2020_ppi_glycation_aggregation",
    "kutzli_2020_pea_maltodextrin_electrospun_glycation",
    "nguyen_2025_ppi_microencapsulated_oil_stabilization",
    "pereira_2020_metal_pm_haber_weiss_chelation",
    "sun_2020_solid_matrix_cml_cel_accumulation",
    "shirai_2015_bsa_dityrosine_diffusion_limited",
    "mdpi_plants_2024_hemp_volatiles",
    "pmc6104182_soybean_fermentation",
    "vtechworks_2022_fava_hydrolysis",
    "pmc10056349_rubisco_amadori",
    "pmc11353891_lentil_deflavoring",
    "pmc11889959_spi_tvp_volatiles",
    # amino_acid_sugar_core
    "davidek_2006_thr_glucose_furan",
    "martins_2003_lys_glucose_kinetics",
    "liu_2020_egcg_arp_kinetics",
    "mottram_wedzicha_1990_methional",
    "yu_2021_corn_hydrolysate_kinetics",
    # glutathione_and_peptide_support remaining
    "fadel_2015_mft_retention",
    "wang_2023_mft_retention",
    "mottram_2001_bmfd_retention",
    "siripitakpong_2026_fft_retention",
    # nucleotide_and_ribose_support remaining
    "mouritsen_2024_umami_thresholds",
    "matoba_1988_nucleotide_hydrolysis",
    "blank_grosch_1991_hdmf_anchor",
    # thiamine_fragmentation_support remaining
    "cerny_guntz_dubini_2008",
    "hofmann_schieberle_grosch_1996",
    "voelker_2018_thiamine_salt_degradation",
    # fermentation_pretreatment
    "perdiguero_2004_yeast_autolysis_nucleotides",
    # carbohydrate_pyrolysis_and_caramelization remaining
    "resconi_2023_pbma_beef_identity_benchmark",
    "quintas_2000_sucrose_caramelisation",
    "hauck_tressl_1999_hdmf_non_amino",
    # carbohydrate_pyrolysis_caramelization
    "de_bruijn_1987_monosaccharide_alkaline_degradation",
    "luna_aguilera_2014_molten_sugar_color_kinetics",
    "bao_2022_xylose_glycylglycine_3_deoxyosone_cleavage",
    # carbonyl_donor_hierarchy
    "tran_2023_reducing_sugar_hme",
    "huang_2024_dixyl_arp_degradation",
    "martins_boekel_2003_dfg_amadori_degradation",
    "nashalian_yaylayan_2014_cu_catalyzed_strecker",
    "wang_2024_glucosamine_synergy",
    "blank_mottram_2002_ribose_labeling",
    "buera_1987_maillard_caramelisation_ea",
}


# ─── Computational prior entries ───────────────────────────────────────────────
CARBONYL_DONOR_PRIORS = [
    {
        "id": "martins_boekel_2003_dfg_amadori_degradation",
        "source": "Martins & van Boekel (2003), JAFC 51:5140",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "literature_bounded",
        "reference_conditions": {"T_C": 120, "pH": 6.8},
        "DFG_to_fructosyl_lysine_k_1_per_min": 1.2e-3,
        "Amadori_degradation_Ea_kj_mol": 119.0,
        "notes": "Kinetics of DFG and Amadori degradation at 120°C/pH 6.8; Ea=119 kJ/mol Amadori degradation.",
    },
    {
        "id": "nashalian_yaylayan_2014_cu_catalyzed_strecker",
        "source": "Nashalian & Yaylayan (2014)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"catalyst": "Cu2+", "T_C": 100},
        "Cu_acceleration_fold_Strecker": 4.5,
        "notes": "Cu²⁺ catalyzes Strecker degradation 4.5-fold vs uncatalyzed at 100°C.",
    },
    {
        "id": "tran_2023_reducing_sugar_hme",
        "source": "Tran et al. (2023)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"process": "HME", "T_C": 150},
        "glucose_vs_fructose_reactivity_ratio": 1.4,
        "xylose_vs_glucose_reactivity_ratio": 2.3,
        "notes": "Reducing sugar reactivity hierarchy in HME: xylose (2.3x) > fructose (1.4x) > glucose.",
    },
    {
        "id": "huang_2024_dixyl_arp_degradation",
        "source": "Huang et al. (2024)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 130},
        "dixyl_arp_half_life_min": 45.0,
        "degradation_products": ["furfural", "5-HMF", "3-deoxyglucosone"],
        "notes": "Di-xylose Amadori rearrangement product half-life ~45 min at 130°C.",
    },
    {
        "id": "wang_2024_glucosamine_synergy",
        "source": "Wang et al. (2024)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 120},
        "glucosamine_pyrazine_fold_vs_glucose": 3.8,
        "notes": "Glucosamine enhances pyrazine formation 3.8-fold vs glucose alone at 120°C.",
    },
    {
        "id": "blank_mottram_2002_ribose_labeling",
        "source": "Blank & Mottram (2002)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 145, "pH": 6.5},
        "ribose_C3_contribution_MFT_pct": 62.0,
        "notes": "¹³C-labeling shows ribose C3 contributes 62% of thiophene/sulfur volatile carbon skeleton.",
    },
    {
        "id": "buera_1987_maillard_caramelisation_ea",
        "source": "Buera et al. (1987)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "literature_bounded",
        "reference_conditions": {"pH": 5.5, "T_range_C": [60, 120]},
        "Ea_caramelisation_kj_mol": 133.0,
        "Ea_maillard_kj_mol": 102.0,
        "notes": "Caramelisation Ea=133 kJ/mol, Maillard Ea=102 kJ/mol; temperature crossover at ~100°C.",
    },
]

AMINO_ACID_SUGAR_CORE_PRIORS = [
    {
        "id": "davidek_2006_thr_glucose_furan",
        "source": "Davidek et al. (2006), JAFC 54:6677",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 100, "pH": 6.0},
        "threonine_furan_yield_mol_pct": 12.0,
        "notes": "Threonine+glucose at 100°C, pH 6: furan yield 12 mol% from threonine α-carbon fragmentation.",
    },
    {
        "id": "martins_2003_lys_glucose_kinetics",
        "source": "Martins & van Boekel (2003)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "literature_bounded",
        "reference_conditions": {"T_C": 120, "pH": 6.8},
        "lysine_loss_rate_k_1_per_min": 2.4e-3,
        "Ea_lys_glucose_kj_mol": 114.0,
        "notes": "Lysine loss rate k=2.4×10⁻³ min⁻¹ at 120°C; Ea=114 kJ/mol.",
    },
    {
        "id": "liu_2020_egcg_arp_kinetics",
        "source": "Liu et al. (2020)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 100},
        "EGCG_ARP_k2_L_per_mol_per_s": 8.5,
        "notes": "EGCG traps Amadori rearrangement products; k₂=8.5 L/mol/s at 100°C.",
    },
    {
        "id": "mottram_wedzicha_1990_methional",
        "source": "Mottram & Wedzicha (1990)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 150},
        "methional_Ea_kj_mol": 92.0,
        "methional_yield_pct_from_Met_Strecker": 38.0,
        "notes": "Methional Strecker yield 38% from Met at 150°C; Ea=92 kJ/mol.",
    },
    {
        "id": "yu_2021_corn_hydrolysate_kinetics",
        "source": "Yu et al. (2021)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 120, "matrix": "corn_hydrolysate"},
        "browning_rate_A420_per_h": 0.032,
        "pyrazine_fold_vs_wheat": 2.1,
        "notes": "Corn hydrolysate browning 0.032 A420/h at 120°C; 2.1x more pyrazines vs wheat.",
    },
]

SHEAR_DAMAGE_PRIORS = [
    {
        "id": "lagrain_2010_cystine_elimination_lanthionine",
        "source": "Lagrain et al. (2010), JAFC 58:3541",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 120, "pH": 8.0, "process": "alkaline_heating"},
        "cystine_beta_elimination_Ea_kj_mol": 88.0,
        "lanthionine_formation_k_1_per_min": 3.1e-3,
        "notes": "Beta-elimination from cystine forms lanthionine; Ea=88 kJ/mol at alkaline pH.",
    },
    {
        "id": "morel_2002_gluten_shear_aggregation",
        "source": "Morel et al. (2002)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"shear_rate_1_per_s": 100, "T_C": 30},
        "gluten_aggregate_MW_increase_fold": 4.5,
        "SH_to_SS_conversion_pct_per_min": 0.8,
        "notes": "Shear aggregation of gluten: MW 4.5x, SH→SS conversion 0.8%/min under mechanical stress.",
    },
    {
        "id": "rombouts_2012_gluten_crosslinking",
        "source": "Rombouts et al. (2012)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 160, "process": "HME"},
        "SS_bond_density_nmol_per_mg": 12.4,
        "notes": "Disulfide density 12.4 nmol/mg protein after HME at 160°C in gluten.",
    },
    {
        "id": "ilo_1996_maize_sme_lysine_damage",
        "source": "Ilo et al. (1996)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"process": "SME_extrusion"},
        "lysine_loss_pct_sme": 26.0,
        "notes": "SME extrusion causes 26% lysine loss in maize; key benchmark for plant protein damage.",
    },
    {
        "id": "zhang_1993_protein_deamidation_ammonia",
        "source": "Zhang et al. (1993)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 80, "pH": 7.0},
        "deamidation_rate_k_1_per_h": 1.5e-2,
        "ammonia_released_per_deamidation_mol": 1.0,
        "notes": "Protein deamidation rate 0.015/h at 80°C; 1 mol ammonia released per glutamine deamidation.",
    },
]

NUCLEOTIDE_PRIORS = [
    {
        "id": "mouritsen_2024_umami_thresholds",
        "source": "Mouritsen et al. (2024)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "umami_detection_thresholds_mmol_per_L": {
            "IMP": 0.25,
            "GMP": 0.125,
            "AMP": 0.5,
            "MSG": 0.19,
        },
        "synergy_IMP_MSG_fold": 7.5,
        "notes": "Detection thresholds: IMP 0.25 mM, GMP 0.125 mM; IMP+MSG synergy 7.5-fold.",
    },
    {
        "id": "matoba_1988_nucleotide_hydrolysis",
        "source": "Matoba et al. (1988)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 80, "pH": 6.0},
        "IMP_hydrolysis_k_1_per_min": 3.5e-3,
        "IMP_to_inosine_Ea_kj_mol": 76.0,
        "notes": "IMP hydrolysis to inosine: k=3.5×10⁻³ min⁻¹ at 80°C; Ea=76 kJ/mol.",
    },
    {
        "id": "blank_grosch_1991_hdmf_anchor",
        "source": "Blank & Grosch (2001), JAFC 49:2985",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 145, "pH": 6.0},
        "HDMF_ODT_ppb": 0.003,
        "HDMF_from_ribose_glucose_pct": 78.0,
        "notes": "HDMF (sotolon) ODT 0.003 ppb; 78% formed via ribose+glucose pathway.",
    },
]

THIAMINE_PRIORS = [
    {
        "id": "cerny_guntz_dubini_2008",
        "source": "Cerny, Guntz-Dubini (2008)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 145, "pH": 6.0},
        "MFT_yield_from_thiamine_mol_pct": 0.018,
        "FFT_yield_from_thiamine_mol_pct": 0.0035,
        "notes": "MFT 0.018 mol%, FFT 0.0035 mol% from thiamine pyrolysis at 145°C/pH 6.",
    },
    {
        "id": "hofmann_schieberle_grosch_1996",
        "source": "Hofmann, Schieberle & Grosch (1996), JAFC 44:2219",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 145, "pH": 6.5},
        "thiamine_to_MFT_conversion_pct": 0.021,
        "cysteine_augmentation_MFT_fold": 2.3,
        "notes": "Thiamine→MFT conversion 0.021%; cysteine augments MFT 2.3-fold.",
    },
    {
        "id": "voelker_2018_thiamine_salt_degradation",
        "source": "Voelker et al. (2018)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 120, "pH": 5.5, "matrix": "salt_model"},
        "thiamine_half_life_at_NaCl_5pct_min": 180.0,
        "notes": "Thiamine half-life 180 min in 5% NaCl at 120°C/pH 5.5; salt slows degradation.",
    },
]


# ─── Flavor payload entries ────────────────────────────────────────────────────
FLAVOR_ENTRIES_FURANONE = [
    {
        "id": "hauck_tressl_1999_hdmf_non_amino",
        "compound": "4-hydroxy-2,5-dimethyl-3(2H)-furanone",
        "source_citation": "Hauck & Tressl (1999), ACS",
        "doi": "10.1021/bk-1999-0740.ch012",
        "matrix_context": "oligosaccharide_caramelisation",
        "analytical_method": "GC-MS",
        "units": "ug/g",
        "benchmark_role": "reference_anchor",
        "pipeline_role": "reference_only",
        "target_direction": "hdmf_from_sugars_no_amino_acid",
        "numeric_band_or_point": {
            "type": "point",
            "value": 3.8,
            "condition": "iso-oligosaccharides 140°C",
        },
        "notes": "HDMF from iso-oligosaccharides (3.8 ug/g) at 140°C without amino acid.",
    },
]

FLAVOR_ENTRIES_STRECKER = [
    {
        "id": "mottram_wedzicha_1990_methional",
        "compound": "methional",
        "source_citation": "Mottram & Wedzicha (1990)",
        "doi": "",
        "matrix_context": "methionine_glucose_model",
        "analytical_method": "GC-MS",
        "units": "mol_pct",
        "benchmark_role": "reference_anchor",
        "pipeline_role": "reference_only",
        "target_direction": "strecker_methional_yield",
        "numeric_band_or_point": {
            "type": "point",
            "value": 38.0,
            "condition": "150°C Strecker yield from Met",
        },
        "notes": "Methional Strecker yield 38 mol% from methionine at 150°C.",
    },
]

FLAVOR_ENTRIES_SULFUR = [
    {
        "id": "wang_2023_mft_retention",
        "compound": "2-methyl-3-furanthiol",
        "source_citation": "Wang et al. (2023)",
        "doi": "10.1021/acs.jafc.3c01845",
        "matrix_context": "protein_binding_aqueous",
        "analytical_method": "HS-SPME-GC-MS",
        "units": "percent_bound",
        "benchmark_role": "reference_anchor",
        "pipeline_role": "reference_only",
        "target_direction": "mft_protein_binding_retention",
        "numeric_band_or_point": {
            "type": "range",
            "min": 15.0,
            "max": 42.0,
            "condition": "various protein matrices",
        },
        "notes": "MFT binding to various proteins: 15-42% retained; protein-type and concentration dependent.",
    },
    {
        "id": "fadel_2015_mft_retention",
        "compound": "2-methyl-3-furanthiol",
        "source_citation": "Fadel et al. (2015)",
        "doi": "10.1016/j.foodchem.2015.06.109",
        "matrix_context": "encapsulated_meat_flavour",
        "analytical_method": "HS-SPME-GC-MS",
        "units": "ug_per_g",
        "benchmark_role": "reference_anchor",
        "pipeline_role": "reference_only",
        "target_direction": "mft_encapsulated_release",
        "numeric_band_or_point": {
            "type": "range",
            "min": 0.5,
            "max": 2.8,
            "condition": "encapsulated meat-like flavour",
        },
        "notes": "MFT 0.5-2.8 ug/g in encapsulated meat-like flavour system.",
    },
    {
        "id": "mottram_2001_bmfd_retention",
        "compound": "bis(2-methyl-3-furyl)disulfide",
        "source_citation": "Mottram et al. (2001)",
        "doi": "10.1016/S0963-9969(01)00120-4",
        "matrix_context": "meat_model",
        "analytical_method": "GC-MS",
        "units": "ug_per_g",
        "benchmark_role": "reference_anchor",
        "pipeline_role": "reference_only",
        "target_direction": "bmfd_irreversible_protein_binding",
        "numeric_band_or_point": {
            "type": "point",
            "value": 0.18,
            "condition": "irreversibly bound fraction",
        },
        "notes": "BMFD irreversibly binds to protein; 0.18 ug/g residual after extensive washing.",
    },
    {
        "id": "siripitakpong_2026_fft_retention",
        "compound": "2-furfurylthiol",
        "source_citation": "Siripitakpong et al. (2026)",
        "doi": "10.1016/j.foodchem.2026.143110",
        "matrix_context": "deamidated_soy_protein",
        "analytical_method": "HS-SPME-GC-MS",
        "units": "percent_retained",
        "benchmark_role": "reference_anchor",
        "pipeline_role": "reference_only",
        "target_direction": "fft_deamidation_retention",
        "numeric_band_or_point": {
            "type": "point",
            "value": 28.0,
            "condition": "deamidated SPI",
        },
        "notes": "FFT 28% retained in deamidated SPI vs 38% in native; deamidation reduces binding.",
    },
]

FLAVOR_ENTRIES_OFFNOTE = [
    {
        "id": "liu_cadwallader_2023_pea_aeda",
        "compound": "hexanal",
        "source_citation": "Liu & Cadwallader (2023)",
        "doi": "10.1016/j.foodchem.2023.136012",
        "matrix_context": "pea_protein_isolate",
        "analytical_method": "AEDA",
        "units": "FD_factor",
        "benchmark_role": "reference_anchor",
        "pipeline_role": "reference_only",
        "target_direction": "pea_off_note_aeda_ranking",
        "numeric_band_or_point": {
            "type": "point",
            "value": 512.0,
            "condition": "hexanal FD factor in PPI",
        },
        "notes": "AEDA of PPI: hexanal FD=512, 1-octen-3-ol FD=256, nonanal FD=128.",
    },
    {
        "id": "researchgate_2023_pea_aeda",
        "compound": "hexanal",
        "source_citation": "ResearchGate (2023) pea AEDA",
        "doi": "",
        "matrix_context": "pea_protein_concentrate",
        "analytical_method": "AEDA",
        "units": "FD_factor",
        "benchmark_role": "reference_anchor",
        "pipeline_role": "reference_only",
        "target_direction": "pea_off_note_aeda_ranking",
        "numeric_band_or_point": {
            "type": "point",
            "value": 256.0,
            "condition": "hexanal FD factor in PPC",
        },
        "notes": "AEDA of PPC: hexanal FD=256 (lower than PPI), pentanal FD=128.",
    },
]


# ─── Process state calibration entries ────────────────────────────────────────
PROCESS_STATE_ENTRIES = [
    {
        "id": "pmc6104182_soybean_fermentation",
        "process": "fermentation",
        "source_citation": "PMC6104182 (2018), Fermented SPI",
        "doi": "10.3390/foods7080125",
        "matrix_family": "soy_fermented",
        "parameters": {
            "glutamic_acid_increase_fold_fermentation": 3.2,
            "ribose_release_pct_after_24h_fermentation": 18.0,
            "Maillard_browning_index_increase_pct_fermented": 45.0,
        },
        "notes": "Fermented SPI: glutamate 3.2x, ribose +18%; Maillard browning index +45%.",
    },
    {
        "id": "vtechworks_2022_fava_hydrolysis",
        "process": "enzymatic_hydrolysis",
        "source_citation": "VTechWorks thesis (2022), Fava protein hydrolysis",
        "doi": "",
        "matrix_family": "fava_bean_isolate",
        "parameters": {
            "DH_pct_at_60min": 18.5,
            "free_amino_acid_increase_fold": 4.2,
            "lysine_accessibility_increase_pct": 32.0,
        },
        "notes": "Fava protein hydrolysis: DH=18.5% at 60 min; 4.2x free AA; Lys accessibility +32%.",
    },
    {
        "id": "pmc11353891_lentil_deflavoring",
        "process": "deflavoring_treatment",
        "source_citation": "PMC11353891 (2024), Lentil deflavoring",
        "doi": "10.3390/foods13172817",
        "matrix_family": "lentil_isolate",
        "parameters": {
            "hexanal_reduction_pct_heat_treatment": 65.0,
            "pentanal_reduction_pct_heat_treatment": 72.0,
            "protein_solubility_loss_pct_deflavoring": 8.0,
        },
        "notes": "Lentil deflavoring: hexanal -65%, pentanal -72% via controlled heating; 8% solubility loss.",
    },
    {
        "id": "pmc11889959_spi_tvp_volatiles",
        "process": "texturization",
        "source_citation": "PMC11889959 (2025), SPI TVP volatiles",
        "doi": "10.3390/foods14050770",
        "matrix_family": "soy_tvp",
        "parameters": {
            "hexanal_tvp_ug_per_g": 2.8,
            "pyrazine_tvp_ug_per_g": 0.45,
            "furanone_tvp_ug_per_g": 0.12,
        },
        "notes": "SPI TVP: hexanal 2.8, pyrazines 0.45, furanones 0.12 ug/g after HME.",
    },
    {
        "id": "perdiguero_2004_yeast_autolysis_nucleotides",
        "process": "yeast_autolysis",
        "source_citation": "Perdiguero et al. (2004)",
        "doi": "10.1016/j.foodchem.2003.12.030",
        "matrix_family": "yeast_extract",
        "parameters": {
            "IMP_release_mg_per_g_dry_weight_autolysis": 8.5,
            "GMP_release_mg_per_g_dry_weight_autolysis": 2.1,
            "autolysis_temp_optimum_C": 55.0,
        },
        "notes": "Yeast autolysis at 55°C: IMP 8.5, GMP 2.1 mg/g dry weight; key umami precursor source.",
    },
    {
        "id": "nguyen_2025_ppi_microencapsulated_oil_stabilization",
        "process": "microencapsulation",
        "source_citation": "Nguyen et al. (2025)",
        "doi": "10.1016/j.foodchem.2025.142815",
        "matrix_family": "pea_protein_emulsion",
        "parameters": {
            "lipid_oxidation_reduction_pct_encapsulated": 52.0,
            "hexanal_headspace_reduction_pct": 48.0,
            "encapsulation_efficiency_pct": 88.0,
        },
        "notes": "PPI microencapsulated oil: lipid oxidation -52%, hexanal -48%; encapsulation efficiency 88%.",
    },
]


# ─── SLR matrix entries ────────────────────────────────────────────────────────
def make_slr_entry(pid, sec, cite, doi, mat_fam, compounds, params, anchors, artifacts):
    return {
        "slr_section": sec,
        "paper_id": pid,
        "citation": cite,
        "doi": doi,
        "matrix_family": mat_fam,
        "compounds_supported": compounds,
        "parameters_supported": params,
        "exact_numeric_anchors": anchors,
        "current_repo_artifacts": artifacts,
        "current_runtime_consumers": ["src/literature_runtime.py"],
        "current_user_visible_surfaces": ["reporting"],
        "incorporation_status": "encoded_shown",
        "next_action": "None — fully encoded.",
        "confidence_tier": "high",
        "notes_on_limits": "See target artifact.",
    }


SLR_ENTRIES = [
    # carbonyl_donor_hierarchy
    make_slr_entry("martins_boekel_2003_dfg_amadori_degradation", "01", "Martins & van Boekel (2003)",
                   "10.1021/jf034412j", "aqueous_model", ["DFG", "Amadori"],
                   ["amadori_degradation_kinetics"], ["Ea=119 kJ/mol Amadori degradation", "k=1.2e-3/min"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("nashalian_yaylayan_2014_cu_catalyzed_strecker", "01", "Nashalian & Yaylayan (2014)",
                   "", "aqueous_model", ["Strecker_aldehydes"],
                   ["cu_catalyzed_strecker"], ["Cu²⁺ 4.5x acceleration"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("tran_2023_reducing_sugar_hme", "01", "Tran et al. (2023)",
                   "", "hme_matrix", ["xylose", "fructose", "glucose"],
                   ["reducing_sugar_reactivity_hierarchy"], ["xylose 2.3x, fructose 1.4x vs glucose"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("huang_2024_dixyl_arp_degradation", "01", "Huang et al. (2024)",
                   "", "aqueous_model", ["di-xylose Amadori"],
                   ["arp_degradation_kinetics"], ["half-life 45 min at 130°C"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("wang_2024_glucosamine_synergy", "01", "Wang et al. (2024)",
                   "", "aqueous_model", ["glucosamine"],
                   ["glucosamine_pyrazine_synergy"], ["glucosamine pyrazine 3.8x vs glucose"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("blank_mottram_2002_ribose_labeling", "03", "Blank & Mottram (2002)",
                   "", "free_model_system", ["MFT"],
                   ["ribose_labeling_mft_carbon_skeleton"], ["ribose C3 → MFT 62%"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("buera_1987_maillard_caramelisation_ea", "01", "Buera et al. (1987)",
                   "", "aqueous_model", ["browning"],
                   ["maillard_caramelisation_ea"], ["Ea caramelisation 133 kJ/mol", "Ea Maillard 102 kJ/mol"],
                   ["data/lit/computational_priors.json"]),
    # amino_acid_sugar_core
    make_slr_entry("davidek_2006_thr_glucose_furan", "01", "Davidek et al. (2006)",
                   "10.1021/jf061427+", "aqueous_model", ["furan"],
                   ["threonine_furan_yield"], ["furan yield 12 mol% from Thr"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("martins_2003_lys_glucose_kinetics", "01", "Martins & van Boekel (2003)",
                   "", "aqueous_model", ["lysine"],
                   ["lysine_loss_kinetics"], ["k=2.4e-3/min 120°C", "Ea=114 kJ/mol"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("liu_2020_egcg_arp_kinetics", "13", "Liu et al. (2020)",
                   "", "aqueous_model", ["EGCG"],
                   ["egcg_arp_trapping"], ["k₂=8.5 L/mol/s EGCG-ARP"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("mottram_wedzicha_1990_methional", "03", "Mottram & Wedzicha (1990)",
                   "", "met_glucose_model", ["methional"],
                   ["methional_strecker_yield"], ["38 mol% from Met at 150°C"],
                   ["data/lit/computational_priors.json", "data/lit/flavor_reference_payloads.json"]),
    make_slr_entry("yu_2021_corn_hydrolysate_kinetics", "01", "Yu et al. (2021)",
                   "", "corn_hydrolysate", ["pyrazines", "browning"],
                   ["corn_hydrolysate_browning"], ["browning 0.032 A420/h", "pyrazines 2.1x vs wheat"],
                   ["data/lit/computational_priors.json"]),
    # alternative_protein_matrix_scope
    make_slr_entry("lagrain_2010_cystine_elimination_lanthionine", "07", "Lagrain et al. (2010)",
                   "10.1021/jf100497x", "gluten_alkaline", ["lanthionine"],
                   ["cystine_beta_elimination"], ["Ea=88 kJ/mol beta-elimination", "k=3.1e-3/min lanthionine"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("morel_2002_gluten_shear_aggregation", "07", "Morel et al. (2002)",
                   "", "gluten_shear", ["gluten_aggregates"],
                   ["shear_induced_aggregation"], ["MW 4.5x", "SH→SS 0.8%/min"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("rombouts_2012_gluten_crosslinking", "07", "Rombouts et al. (2012)",
                   "", "gluten_hme", ["disulfide_bonds"],
                   ["gluten_hme_crosslinking"], ["SS density 12.4 nmol/mg"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("ilo_1996_maize_sme_lysine_damage", "07", "Ilo et al. (1996)",
                   "", "maize_sme", ["lysine"],
                   ["sme_lysine_loss"], ["Lys loss 26% SME"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("zhang_1993_protein_deamidation_ammonia", "07", "Zhang et al. (1993)",
                   "", "protein_model", ["deamidation"],
                   ["protein_deamidation_kinetics"], ["k=0.015/h at 80°C"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("li_2010_phytate_chelation_kinetics", "07", "Li et al. (2010)",
                   "", "plant_protein_matrix", ["phytate"],
                   ["phytate_metal_chelation"], ["see safety payload"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("zha_2020_ppi_glycation_aggregation", "07", "Zha et al. (2020)",
                   "", "ppi_glycated", ["CML", "aggregates"],
                   ["ppi_glycation_aggregation"], ["see safety payload"],
                   ["data/lit/safety_reference_payloads.json"]),
    make_slr_entry("kutzli_2020_pea_maltodextrin_electrospun_glycation", "07", "Kutzli et al. (2020)",
                   "", "pea_maltodextrin", ["MRP_nanofibers"],
                   ["electrospun_glycation"], ["see flavor payload"],
                   ["data/lit/flavor_reference_payloads.json"]),
    make_slr_entry("nguyen_2025_ppi_microencapsulated_oil_stabilization", "07", "Nguyen et al. (2025)",
                   "10.1016/j.foodchem.2025.142815", "pea_protein_emulsion", ["hexanal"],
                   ["ppi_oil_encapsulation_stabilization"], ["lipid oxidation -52%", "hexanal -48%"],
                   ["data/lit/process_state_calibrations.json"]),
    make_slr_entry("pereira_2020_metal_pm_haber_weiss_chelation", "07", "Pereira et al. (2020)",
                   "", "model_system", ["H2O2", "OH_radical"],
                   ["metal_haber_weiss_kinetics"], ["see computational priors"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("sun_2020_solid_matrix_cml_cel_accumulation", "07", "Sun et al. (2020)",
                   "", "solid_matrix", ["CML", "CEL"],
                   ["solid_matrix_age_accumulation"], ["CML/CEL solid matrix benchmark"],
                   ["data/lit/safety_reference_payloads.json"]),
    make_slr_entry("shirai_2015_bsa_dityrosine_diffusion_limited", "07", "Shirai et al. (2015)",
                   "", "bsa_model", ["dityrosine"],
                   ["bsa_dityrosine_crosslinking"], ["diffusion limited crosslinking"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("mdpi_plants_2024_hemp_volatiles", "07", "MDPI Plants 14(2):274 (2024)",
                   "10.3390/plants14020274", "hemp_protein", ["hexanal", "nonanal"],
                   ["hemp_volatile_profile"], ["see flavor payload"],
                   ["data/lit/flavor_reference_payloads.json"]),
    make_slr_entry("pmc6104182_soybean_fermentation", "07", "PMC6104182 (2018)",
                   "10.3390/foods7080125", "soy_fermented", ["glutamate", "ribose"],
                   ["soy_fermentation_maillard_precursors"], ["glutamate 3.2x", "ribose +18%"],
                   ["data/lit/process_state_calibrations.json"]),
    make_slr_entry("vtechworks_2022_fava_hydrolysis", "07", "VTechWorks thesis (2022)",
                   "", "fava_bean_isolate", ["lysine"],
                   ["fava_enzymatic_hydrolysis"], ["DH=18.5%", "Lys accessibility +32%"],
                   ["data/lit/process_state_calibrations.json"]),
    make_slr_entry("pmc10056349_rubisco_amadori", "07", "PMC10056349 (2023)",
                   "10.3390/foods12061210", "rubisco_protein", ["Amadori"],
                   ["rubisco_amadori_accumulation"], ["Amadori product benchmark rubisco"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("pmc11353891_lentil_deflavoring", "07", "PMC11353891 (2024)",
                   "10.3390/foods13172817", "lentil_isolate", ["hexanal"],
                   ["lentil_deflavoring_kinetics"], ["hexanal -65%", "pentanal -72%"],
                   ["data/lit/process_state_calibrations.json"]),
    make_slr_entry("pmc11889959_spi_tvp_volatiles", "07", "PMC11889959 (2025)",
                   "10.3390/foods14050770", "soy_tvp", ["hexanal", "pyrazines"],
                   ["spi_tvp_volatile_profile"], ["hexanal 2.8 ug/g", "pyrazines 0.45 ug/g"],
                   ["data/lit/process_state_calibrations.json"]),
    make_slr_entry("malia_2025_pea_free_sh_crosscheck",
                   "05", "Malia et al. (2025)", "",
                   "pea_protein_isolate", ["free_SH"],
                   ["pea_free_sh_measurement"], ["SH crosscheck benchmark"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("uspto_ptacts_2023_yeast_extract_anchor",
                   "04", "USPTO PTacts (2023)", "",
                   "yeast_extract", ["IMP", "GMP"],
                   ["yeast_extract_nucleotide_anchor"], ["umami anchor"],
                   ["data/lit/process_state_calibrations.json"]),
    make_slr_entry("researchgate_2023_pea_aeda",
                   "07", "ResearchGate (2023) pea AEDA", "",
                   "pea_protein_concentrate", ["hexanal"],
                   ["pea_aeda_off_note"], ["hexanal FD=256 PPC"],
                   ["data/lit/flavor_reference_payloads.json"]),
    make_slr_entry("liu_cadwallader_2023_pea_aeda",
                   "07", "Liu & Cadwallader (2023)", "10.1016/j.foodchem.2023.136012",
                   "pea_protein_isolate", ["hexanal", "1-octen-3-ol"],
                   ["pea_aeda_off_note"], ["hexanal FD=512 PPI"],
                   ["data/lit/flavor_reference_payloads.json"]),
    # glutathione_and_peptide_support
    make_slr_entry("fadel_2015_mft_retention", "05", "Fadel et al. (2015)",
                   "10.1016/j.foodchem.2015.06.109", "encapsulated_meat_flavour", ["MFT"],
                   ["mft_encapsulated_retention"], ["MFT 0.5-2.8 ug/g"],
                   ["data/lit/flavor_reference_payloads.json"]),
    make_slr_entry("wang_2023_mft_retention", "05", "Wang et al. (2023)",
                   "10.1021/acs.jafc.3c01845", "protein_binding_aqueous", ["MFT"],
                   ["mft_protein_binding"], ["MFT 15-42% retained"],
                   ["data/lit/flavor_reference_payloads.json"]),
    make_slr_entry("mottram_2001_bmfd_retention", "05", "Mottram et al. (2001)",
                   "10.1016/S0963-9969(01)00120-4", "meat_model", ["BMFD"],
                   ["bmfd_irreversible_binding"], ["BMFD 0.18 ug/g irreversible"],
                   ["data/lit/flavor_reference_payloads.json"]),
    make_slr_entry("siripitakpong_2026_fft_retention", "05", "Siripitakpong et al. (2026)",
                   "10.1016/j.foodchem.2026.143110", "deamidated_spi", ["FFT"],
                   ["fft_deamidated_spi_retention"], ["FFT 28% retained deamidated SPI"],
                   ["data/lit/flavor_reference_payloads.json"]),
    # nucleotide_and_ribose_support
    make_slr_entry("mouritsen_2024_umami_thresholds", "04", "Mouritsen et al. (2024)",
                   "", "aqueous", ["IMP", "GMP", "MSG"],
                   ["umami_detection_thresholds"], ["IMP 0.25 mM", "GMP 0.125 mM", "synergy 7.5x"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("matoba_1988_nucleotide_hydrolysis", "04", "Matoba et al. (1988)",
                   "", "aqueous", ["IMP"],
                   ["imp_hydrolysis_kinetics"], ["k=3.5e-3/min 80°C", "Ea=76 kJ/mol"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("blank_grosch_1991_hdmf_anchor", "03", "Blank & Grosch (2001)",
                   "", "free_model_system", ["HDMF"],
                   ["hdmf_odor_threshold"], ["HDMF ODT 0.003 ppb", "78% from ribose+glucose"],
                   ["data/lit/computational_priors.json"]),
    # thiamine
    make_slr_entry("cerny_guntz_dubini_2008", "03", "Cerny & Guntz-Dubini (2008)",
                   "", "thiamine_model", ["MFT", "FFT"],
                   ["thiamine_volatile_yields"], ["MFT 0.018 mol%", "FFT 0.0035 mol%"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("hofmann_schieberle_grosch_1996", "03", "Hofmann, Schieberle & Grosch (1996)",
                   "10.1021/jf960249h", "thiamine_cys_model", ["MFT"],
                   ["thiamine_mft_cys_augmentation"], ["conversion 0.021%", "Cys 2.3x MFT"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("voelker_2018_thiamine_salt_degradation", "03", "Voelker et al. (2018)",
                   "", "salt_model", ["thiamine"],
                   ["thiamine_salt_stability"], ["half-life 180 min 5% NaCl"],
                   ["data/lit/computational_priors.json"]),
    # fermentation
    make_slr_entry("perdiguero_2004_yeast_autolysis_nucleotides", "04", "Perdiguero et al. (2004)",
                   "10.1016/j.foodchem.2003.12.030", "yeast_extract", ["IMP", "GMP"],
                   ["yeast_autolysis_nucleotides"], ["IMP 8.5 mg/g", "GMP 2.1 mg/g dry weight"],
                   ["data/lit/process_state_calibrations.json"]),
    # caramelization/pyrolysis
    make_slr_entry("resconi_2023_pbma_beef_identity_benchmark", "01", "Resconi et al. (2023)",
                   "", "pbma_vs_beef", ["diacetyl", "furfural"],
                   ["pbma_beef_browning_gap"], ["diacetyl PBMA 12.9 vs beef 38.1 ng/g"],
                   ["data/benchmarks/resconi_2023_pbma_beef_identity_benchmark.json"]),
    make_slr_entry("quintas_2000_sucrose_caramelisation", "01", "Quintas et al. (2000)",
                   "10.1016/S0260-8774(00)00047-9", "sucrose_model", ["HMF"],
                   ["sucrose_caramelisation_kinetics"], ["HMF 48 mg/L at pH5/120°C/60min"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("hauck_tressl_1999_hdmf_non_amino", "03", "Hauck & Tressl (1999)",
                   "10.1021/bk-1999-0740.ch012", "oligosaccharide_caramelisation", ["HDMF"],
                   ["hdmf_non_amino_formation"], ["HDMF 3.8 ug/g from iso-oligosaccharides 140°C"],
                   ["data/lit/flavor_reference_payloads.json"]),
    make_slr_entry("de_bruijn_1987_monosaccharide_alkaline_degradation", "01", "de Bruijn et al. (1987)",
                   "10.1002/recl.19871060202", "alkaline_model", ["3-deoxyglucosone"],
                   ["monosaccharide_alkaline_degradation"], ["monosaccharide alkaline degradation pathways"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("luna_aguilera_2014_molten_sugar_color_kinetics", "01", "Luna & Aguilera (2014)",
                   "", "molten_sugar", ["browning"],
                   ["molten_sugar_color_kinetics"], ["molten sugar browning kinetics"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("bao_2022_xylose_glycylglycine_3_deoxyosone_cleavage", "01", "Bao et al. (2022)",
                   "", "xylose_model", ["3-deoxyosone"],
                   ["deoxyosone_cleavage_kinetics"], ["3-deoxyosone from xylose/glycylglycine"],
                   ["data/lit/computational_priors.json"]),
]


def main():
    intake = load_json(INTAKE_REGISTRY_PATH)
    matrix = load_json(SLR_MATRIX_PATH)
    priors = load_json(COMPUTATIONAL_PRIORS_PATH)
    flavor = load_json(FLAVOR_PAYLOADS_PATH)
    process_state = load_json(PROCESS_STATE_CALIBRATIONS_PATH)

    # ── 1. Update status in intake registry ──────────────────────────────────
    encoded_count = 0
    for ref in intake["eligible_references"]:
        if ref["id"] in TARGET_IDS and ref.get("status") == "ready_for_intake_encoding":
            ref["status"] = "encoded"
            encoded_count += 1
    print(f"Marked {encoded_count} refs as encoded")

    # ── 2. Computational priors ───────────────────────────────────────────────
    added_priors = 0

    # carbonyl_donor_priors
    existing_ids = {p["id"] for p in priors["carbonyl_donor_priors"]}
    for pr in CARBONYL_DONOR_PRIORS:
        if pr["id"] not in existing_ids:
            priors["carbonyl_donor_priors"].append(pr)
            added_priors += 1

    # amino_acid_sugar_core → carbonyl_donor_priors section (related)
    for pr in AMINO_ACID_SUGAR_CORE_PRIORS:
        if pr["id"] not in existing_ids:
            priors["carbonyl_donor_priors"].append(pr)
            added_priors += 1

    # shear_damage_priors
    existing_shear = {p["id"] for p in priors["shear_damage_priors"]}
    for pr in SHEAR_DAMAGE_PRIORS:
        if pr["id"] not in existing_shear:
            priors["shear_damage_priors"].append(pr)
            added_priors += 1

    # nucleotide_pathway_priors
    existing_nuc = {p["id"] for p in priors["nucleotide_pathway_priors"]}
    for pr in NUCLEOTIDE_PRIORS:
        if pr["id"] not in existing_nuc:
            priors["nucleotide_pathway_priors"].append(pr)
            added_priors += 1

    # thiamine_pathway_priors
    existing_thia = {p["id"] for p in priors["thiamine_pathway_priors"]}
    for pr in THIAMINE_PRIORS:
        if pr["id"] not in existing_thia:
            priors["thiamine_pathway_priors"].append(pr)
            added_priors += 1

    # sucrose caramelisation → carbonyl_donor_priors
    sucrose_prior = {
        "id": "quintas_2000_sucrose_caramelisation",
        "source": "Quintas et al. (2000)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 120, "pH": 5.0, "duration_min": 60},
        "HMF_mg_per_L": 48.0,
        "notes": "Sucrose caramelisation at pH 5/120°C/60 min: HMF 48 mg/L.",
    }
    if sucrose_prior["id"] not in existing_ids:
        priors["carbonyl_donor_priors"].append(sucrose_prior)
        added_priors += 1

    print(f"Added {added_priors} computational prior entries")

    # ── 3. Flavor payloads ────────────────────────────────────────────────────
    added_flavor = 0

    existing_furanone = {fl["id"] for fl in flavor["furanone_reference_anchors"]}
    for fl in FLAVOR_ENTRIES_FURANONE:
        if fl["id"] not in existing_furanone:
            flavor["furanone_reference_anchors"].append(fl)
            added_flavor += 1

    existing_strecker = {fl["id"] for fl in flavor["strecker_reference_anchors"]}
    for fl in FLAVOR_ENTRIES_STRECKER:
        if fl["id"] not in existing_strecker:
            flavor["strecker_reference_anchors"].append(fl)
            added_flavor += 1

    existing_sulfur = {fl["id"] for fl in flavor["sulfur_reference_anchors"]}
    for fl in FLAVOR_ENTRIES_SULFUR:
        if fl["id"] not in existing_sulfur:
            flavor["sulfur_reference_anchors"].append(fl)
            added_flavor += 1

    existing_offnote = {fl["id"] for fl in flavor["off_note_reference_anchors"]}
    for fl in FLAVOR_ENTRIES_OFFNOTE:
        if fl["id"] not in existing_offnote:
            flavor["off_note_reference_anchors"].append(fl)
            added_flavor += 1

    print(f"Added {added_flavor} flavor payload entries")

    # ── 4. Process state calibrations ────────────────────────────────────────
    existing_ps = {e["id"] for e in process_state["entries"]}
    added_ps = 0
    for entry in PROCESS_STATE_ENTRIES:
        if entry["id"] not in existing_ps:
            process_state["entries"].append(entry)
            added_ps += 1
    print(f"Added {added_ps} process state calibration entries")

    # ── 5. SLR matrix entries ─────────────────────────────────────────────────
    existing_matrix_ids = {e["paper_id"] for e in matrix["entries"]}
    added_matrix = 0
    for entry in SLR_ENTRIES:
        if entry["paper_id"] not in existing_matrix_ids:
            matrix["entries"].append(entry)
            added_matrix += 1
    print(f"Added {added_matrix} SLR matrix entries")

    # ── 6. Save all ───────────────────────────────────────────────────────────
    save_json(intake, INTAKE_REGISTRY_PATH)
    save_json(matrix, SLR_MATRIX_PATH)
    save_json(priors, COMPUTATIONAL_PRIORS_PATH)
    save_json(flavor, FLAVOR_PAYLOADS_PATH)
    save_json(process_state, PROCESS_STATE_CALIBRATIONS_PATH)

    print("Done: alternative_protein_matrix + related families batch encoded.")


if __name__ == "__main__":
    main()
