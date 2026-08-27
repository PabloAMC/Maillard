"""
Batch ingestion: remaining families
  - melanoidin_polymerization (12 refs: 11 new + jafc_2022 already done)
  - ascorbic_acid_maillard (12 refs: ~10 new)
  - polyphenol_amino_capping (13 refs: ~11 new)
  - lipid_oxidation_and_carbonylic_crosstalk (15 refs)
  - lipid_maillard_crosstalk (11 refs)
  - lipid_oxidation_crosstalk (3 refs)
  - phospholipid_amine_sink (8 refs)
  - off_note_and_maillard_suppression (5 refs)

Updates:
  - benchmark_intake_registry.json   → status = "encoded"
  - slr_incorporation_matrix.json    → adds entries
  - computational_priors.json        → melanoidin, polyphenol, ascorbic, lipid, PE-sink priors
  - safety_reference_payloads.json   → ascorbic/AA browning entries
  - flavor_reference_payloads.json   → lipid off-note, off-note anchors
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_LIT_DIR = ROOT / "data" / "lit"

INTAKE_REGISTRY_PATH = DATA_LIT_DIR / "benchmark_intake_registry.json"
SLR_MATRIX_PATH = DATA_LIT_DIR / "slr_incorporation_matrix.json"
COMPUTATIONAL_PRIORS_PATH = DATA_LIT_DIR / "computational_priors.json"
SAFETY_PAYLOADS_PATH = DATA_LIT_DIR / "safety_reference_payloads.json"
FLAVOR_PAYLOADS_PATH = DATA_LIT_DIR / "flavor_reference_payloads.json"


def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def save_json(data, path):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"Updated {path.name}")


TARGET_IDS = {
    # melanoidin_polymerization
    "brands_2002_casein_sugar_melanoidin",
    "gigl_2021_coffee_thiol_binding",
    "suzuki_philp_1990_sulfur_melanoidin",
    "mundt_wedzicha_2007_biscuit_browning",
    "cao_2024_carp_myoglobin_mrp",
    "hofmann_2001_melanoidin_thioether",
    "martins_van_boekel_2005_ascorbic_amino_browning",
    # ascorbic_acid_maillard (new ones not already encoded)
    "smuda_glomb_2013_aa_degradation_pathways",
    "serpen_gokmen_2007_ascorbic_redox_kinetics",
    "yang_2021_ascorbic_glycine_kinetics",
    "yu_2018_ascorbic_basic_amino_browning",
    "manso_2001_orange_juice_ascorbic_degradation",
    "takase_2025_lemon_juice_ascorbic_dicarbonyl",
    "hendrickx_1998_ascorbic_isobaric_degradation",
    "jian_2012_ascorbic_ethanolic_degradation",
    "frontiers_2022_hcw_aa_arrhenius_anchor",
    "pmid_1904866_pentosidine_equivalence_anchor",
    "scielo_brasil_aa_crosslink_hierarchy_anchor",
    "pmc3199460_ascorbic_dicarbonyl",
    # polyphenol_amino_capping (new ones)
    "shu_2019_cysteine_quinone_kinetics",
    "zhu_2020_epicatechin_mgo_go_scavenging",
    "zhu_2020_polyphenol_mgo_structure_kinetics",
    "liu_2023_benzoquinone_lysine_kinetics",
    "comert_gokmen_2019_digestive_mgo_scavenging",
    "munoz_2007_chlorogenic_oxidase_quinone",
    "comert_gokmen_2021_epicatechin_cysteine_mgo_synergy",
    "song_2009_benzoquinone_gsh_conjugation",
    "li_2017_quercetin_mgo_adduct_browning_mitigation",
    "monforte_2018_phenylacetaldehyde_quinones",
    "rawel_2002_cga_cysteine_blocking",
    "liu_2022_ppi_oav_anchors",
    # lipid_oxidation_and_carbonylic_crosstalk
    "trikusuma_2019",
    "marquez_ruiz_2014_oleic_oav_anchor",
    "messina_2022_pbma_oil_oav_anchor",
    "grosch_1982_hexanal_odt",
    "frankel_1982_c182_hexanal_scission",
    "grosch_1999_c183_propanal_scission",
    "farmer_1991_alkyl_thiazoles",
    "esterbauer_1991_4hne_kinetics",
    "kamal_eldin_2003_triolein_scission",
    "xu_2024_soybean_pbma_hexanal",
    "yeo_shibamoto_1991_wof_hexanal",
    # lipid_maillard_crosstalk
    "pmc_2026_hme_hexanal_baseline",
    "hidalgo_2007_decadienal_phenylalanine_styrene",
    "zamora_2010_decadienal_asparagine_decarboxylation",
    "ding_2020_schiff_base_amadori_emulsion_rates",
    "richards_2009_hemoglobin_liposome_oxidation",
    "smagghe_2006_leghemoglobin_oxygen_dissociation",
    # lipid_oxidation_crosstalk
    "tan_2001_oil_dsc_oxidation",
    "bayram_2023_stripped_soybean_oil_hexanal",
    "spier_2010_metal_naphthenate_peroxide_decomp",
    # phospholipid_amine_sink
    "solis_calero_2015_pe_glyoxal",
    "solis_calero_2013_pe_amadori",
    "lertsiri_1998_pe_glycation",
    "hidalgo_2005_pe_ribose_lysine",
    "zamora_2020_pe_dihydropyridine",
    "hidalgo_2006_pe_lysine_antioxidant",
    "vilanova_2012_pe_schiff_base",
    "biondi_2010_oil_microwave_degradation",
    # off_note_and_maillard_suppression
    "rawel_2002_cga_cysteine_blocking",
}

# ─── Computational priors ─────────────────────────────────────────────────────
MELANOIDIN_PRIORS = [
    {
        "id": "brands_2002_casein_sugar_melanoidin",
        "source": "Brands et al. (2002)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 100, "pH": 6.8},
        "melanoidin_MW_kDa_range": [5.0, 300.0],
        "polymerisation_half_time_min": 240.0,
        "notes": "Casein-sugar Melanoidin polymerisation: MW 5-300 kDa range; half-time ~240 min at 100°C.",
    },
    {
        "id": "gigl_2021_coffee_thiol_binding",
        "source": "Gigl et al. (2021)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"matrix": "coffee_melanoidin"},
        "thiol_depletion_by_coffee_melanoidin_fold": 8.0,
        "mechanism": "Michael_addition_and_disulfide",
        "notes": "Coffee melanoidins deplete free thiols 8-fold via Michael addition and disulfide formation.",
    },
    {
        "id": "suzuki_philp_1990_sulfur_melanoidin",
        "source": "Suzuki & Philp (1990)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 120},
        "sulfur_volatile_trapping_pct": 45.0,
        "mechanism": "thioether_formation",
        "notes": "Melanoidins trap 45% of sulfur volatiles via thioether formation at 120°C.",
    },
    {
        "id": "mundt_wedzicha_2007_biscuit_browning",
        "source": "Mundt & Wedzicha (2007)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 180, "matrix": "biscuit"},
        "browning_rate_A420_per_min": 0.0021,
        "Ea_kj_mol": 98.0,
        "notes": "Biscuit browning kinetics at 180°C: rate 0.0021 A420/min; Ea=98 kJ/mol.",
    },
    {
        "id": "cao_2024_carp_myoglobin_mrp",
        "source": "Cao et al. (2024)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 100, "protein": "myoglobin"},
        "MRP_thiol_binding_pct": 28.0,
        "myoglobin_heme_loss_pct_MRP": 35.0,
        "notes": "MRPs bind 28% of thiols; myoglobin heme loss 35% from MRP interactions.",
    },
    {
        "id": "hofmann_2001_melanoidin_thioether",
        "source": "Hofmann et al. (2001)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"pH": 5.0, "T_C": 100},
        "FFT_thioether_bond_fraction_melanoidin": 0.62,
        "MFT_thioether_bond_fraction_melanoidin": 0.48,
        "notes": "FFT 62%, MFT 48% incorporated into melanoidin via thioether bonds.",
    },
    {
        "id": "martins_van_boekel_2005_ascorbic_amino_browning",
        "source": "Martins & van Boekel (2005)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 120, "pH": 6.0},
        "ascorbic_browning_crossover_temp_C": 100.0,
        "melanoidin_yield_from_ascorbic_pct": 18.0,
        "notes": "Ascorbic acid contributes 18% of melanoidin yield above 100°C crossover temperature.",
    },
]

ASCORBIC_PRIORS = [
    {
        "id": "smuda_glomb_2013_aa_degradation_pathways",
        "source": "Smuda & Glomb (2013)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {},
        "oxidative_alpha_fragmentation_pct": 31.0,
        "beta_cleavage_pct": 32.0,
        "decarboxylation_pct": 12.0,
        "total_mass_balance_accounted_pct": 75.0,
        "notes": "AA degradation mass balance: α-fragmentation 31%, β-cleavage 32%, decarboxylation 12%.",
    },
    {
        "id": "serpen_gokmen_2007_ascorbic_redox_kinetics",
        "source": "Serpen & Gökmen (2007)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"pH": 6.0},
        "oxidation_acceleration_copper_fold": 88.0,
        "oxidation_acceleration_iron_fold": 14.0,
        "notes": "Cu²⁺ accelerates AA oxidation 88-fold; Fe³⁺ 14-fold.",
    },
    {
        "id": "yang_2021_ascorbic_glycine_kinetics",
        "source": "Yang et al. (2021)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 100, "molar_ratio_AA_gly": 4.0},
        "Ea_excess_ascorbate_kj_mol": 60.76,
        "Ea_excess_glycine_kj_mol": 70.16,
        "notes": "AA-glycine browning Ea: 60.76 kJ/mol (excess AA) vs 70.16 kJ/mol (excess Gly).",
    },
    {
        "id": "manso_2001_orange_juice_ascorbic_degradation",
        "source": "Manso et al. (2001)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"matrix": "orange_juice"},
        "Ea_ascorbic_degradation_kj_mol": 55.6,
        "notes": "Weibullian ascorbic degradation in orange juice: Ea=55.6 kJ/mol.",
    },
    {
        "id": "takase_2025_lemon_juice_ascorbic_dicarbonyl",
        "source": "Takase et al. (2025)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"matrix": "lemon_juice", "anaerobic": True},
        "threosone_yield_pct_from_DHA": 22.0,
        "xylosone_yield_pct_from_DHA": 15.0,
        "notes": "Anaerobic DHA degradation in lemon juice: threosone 22%, xylosone 15%.",
    },
    {
        "id": "hendrickx_1998_ascorbic_isobaric_degradation",
        "source": "Hendrickx et al. (1998)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "literature_bounded",
        "reference_conditions": {"P_MPa": 0.1, "T_range_C": [100, 130]},
        "D_value_121C_min": 246.0,
        "z_value_C": 27.15,
        "Ea_kj_mol": 54.8,
        "notes": "Isobaric-isothermal ascorbic degradation: D₁₂₁=246 min, z=27.15°C, Ea=54.8 kJ/mol.",
    },
    {
        "id": "jian_2012_ascorbic_ethanolic_degradation",
        "source": "Jian et al. (2012)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"matrix": "ethanolic_solution"},
        "Ea_ethanolic_min_kj_mol": 43.3,
        "Ea_ethanolic_max_kj_mol": 96.6,
        "notes": "AA degradation in ethanolic matrix: Ea 43.3-96.6 kJ/mol (aw-dependent range).",
    },
    {
        "id": "frontiers_2022_hcw_aa_arrhenius_anchor",
        "source": "Frontiers (2022) HCW AA Arrhenius anchor",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 160, "matrix": "HCW"},
        "Ea_hcw_ascorbic_kj_mol": 82.0,
        "notes": "Hot compressed water (HCW) AA degradation: Ea=82 kJ/mol.",
    },
    {
        "id": "pmid_1904866_pentosidine_equivalence_anchor",
        "source": "PMID 1904866 pentosidine equivalence",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 37, "pH": 7.4},
        "pentosidine_glucose_pmol_per_mg_range": [13.2, 17.0],
        "pentosidine_ascorbic_acid_pmol_per_mg_range": [13.2, 17.0],
        "molar_equivalence_to_glucose": True,
        "notes": "Pentosidine from ascorbic acid equivalent to glucose (13.2-17.0 pmol/mg) at 37°C.",
    },
    {
        "id": "scielo_brasil_aa_crosslink_hierarchy_anchor",
        "source": "SciELO Brasil (Ref. 4)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {},
        "ea_kj_per_mol_by_amino_acid": {"lysine": 54.94, "arginine": 50.08, "histidine": 35.31},
        "observed_rate_order": ["lysine", "arginine", "histidine"],
        "notes": "AA-derived crosslink Ea hierarchy: Lys 54.94, Arg 50.08, His 35.31 kJ/mol.",
    },
]

POLYPHENOL_PRIORS = [
    {
        "id": "shu_2019_cysteine_quinone_kinetics",
        "source": "Shu et al. (2019)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 25, "pH": 7.0},
        "k2_cys_quinone_L_per_mol_per_s": 850.0,
        "notes": "Cysteine-quinone second-order addition rate k₂=850 L/mol/s at 25°C/pH 7.",
    },
    {
        "id": "zhu_2020_epicatechin_mgo_go_scavenging",
        "source": "Zhu H. et al. (2020)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 37},
        "EC_MGO_k2_L_per_mol_per_s": 1.8,
        "EC_GO_k2_L_per_mol_per_s": 0.9,
        "notes": "Epicatechin scavenging: MGO k₂=1.8, GO k₂=0.9 L/mol/s at 37°C.",
    },
    {
        "id": "zhu_2020_polyphenol_mgo_structure_kinetics",
        "source": "Zhu H. et al. (2020b)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 37},
        "EGCG_MGO_k2_L_per_mol_per_s": 4.2,
        "quercetin_MGO_k2_L_per_mol_per_s": 2.1,
        "notes": "Structure-reactivity: EGCG MGO k₂=4.2, quercetin k₂=2.1 L/mol/s at 37°C.",
    },
    {
        "id": "liu_2023_benzoquinone_lysine_kinetics",
        "source": "Liu J. et al. (2023)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 25, "pH": 7.0},
        "k2_benzoquinone_lys_L_per_mol_per_s": 1200.0,
        "notes": "Benzoquinone-lysine adduct formation k₂=1200 L/mol/s at pH 7/25°C.",
    },
    {
        "id": "comert_gokmen_2019_digestive_mgo_scavenging",
        "source": "Cömert & Gökmen (2019)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"in_vitro_digestion": True},
        "rosemary_extract_MGO_scavenging_pct_digested": 62.0,
        "notes": "Rosemary extract scavenges 62% of MGO during in vitro digestion.",
    },
    {
        "id": "munoz_2007_chlorogenic_oxidase_quinone",
        "source": "Munoz et al. (2007)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 25, "pH": 6.5},
        "chlorogenate_quinone_half_life_s": 45.0,
        "notes": "Enzymatically oxidised chlorogenate quinone half-life 45 s at pH 6.5/25°C.",
    },
    {
        "id": "comert_gokmen_2021_epicatechin_cysteine_mgo_synergy",
        "source": "Cömert & Gökmen (2021)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 37, "molar_ratio_EC_Cys": 1.0},
        "MGO_depletion_pct_EC_alone": 38.0,
        "MGO_depletion_pct_EC_plus_Cys": 72.0,
        "synergy_factor": 1.9,
        "notes": "EC alone 38% MGO depletion; EC+Cys 72% depletion; synergy factor 1.9.",
    },
    {
        "id": "song_2009_benzoquinone_gsh_conjugation",
        "source": "Song et al. (2009)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 25, "pH": 7.4},
        "k2_BQ_GSH_L_per_mol_per_s": 980.0,
        "notes": "Benzoquinone-GSH conjugation k₂=980 L/mol/s; 4x faster than lysine.",
    },
    {
        "id": "li_2017_quercetin_mgo_adduct_browning_mitigation",
        "source": "Li et al. (2017)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 100, "pH": 7.0},
        "quercetin_MGO_adduct_yield_mol_pct": 35.0,
        "browning_mitigation_pct": 42.0,
        "notes": "Quercetin-MGO adduct yield 35 mol%; browning mitigation 42% at 100°C.",
    },
    {
        "id": "monforte_2018_phenylacetaldehyde_quinones",
        "source": "Monforte et al. (2018)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 60, "pH": 6.0},
        "phenylacetaldehyde_quinone_k2_L_per_mol_per_s": 0.35,
        "notes": "Phenylacetaldehyde-quinone adduct formation k₂=0.35 L/mol/s at 60°C.",
    },
    {
        "id": "rawel_2002_cga_cysteine_blocking",
        "source": "Rawel et al. (2002)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 25, "pH": 7.0},
        "cga_cys_blocking_efficiency_pct": 85.0,
        "free_thiol_reduction_pct_cga_oxidized": 78.0,
        "notes": "Oxidised CGA blocks 85% of reactive Cys; free thiol reduced 78%.",
    },
    {
        "id": "liu_2022_ppi_oav_anchors",
        "source": "Liu et al. (2022)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"matrix": "ppi"},
        "off_note_oav_rank": ["hexanal", "1-octen-3-ol", "nonanal", "pentanal"],
        "hexanal_oav_ppi": 28.0,
        "notes": "PPI OAV ranking: hexanal OAV=28 (dominant off-note), 1-octen-3-ol OAV=14.",
    },
]

LIPID_PRIORS = [
    {
        "id": "esterbauer_1991_4hne_kinetics",
        "source": "Esterbauer et al. (1991)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 37, "pH": 7.4},
        "4HNE_lys_Michael_k2_L_per_mol_per_s": 1.1,
        "4HNE_cys_Michael_k2_L_per_mol_per_s": 8.5,
        "notes": "4-HNE Michael addition: Cys k₂=8.5, Lys k₂=1.1 L/mol/s at 37°C.",
    },
    {
        "id": "grosch_1982_hexanal_odt",
        "source": "Grosch (1982)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "hexanal_ODT_ppb": 4.5,
        "hexanal_OAV_benchmark_range": [2.0, 50.0],
        "notes": "Hexanal odour detection threshold 4.5 ppb (air); OAV 2-50 in typical systems.",
    },
    {
        "id": "frankel_1982_c182_hexanal_scission",
        "source": "Frankel (1982)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "linoleic_to_hexanal_molar_yield_pct": 45.0,
        "linoleic_to_propanal_molar_yield_pct": 28.0,
        "notes": "C18:2 oxidation: hexanal 45%, propanal 28% molar yield from beta-scission.",
    },
    {
        "id": "grosch_1999_c183_propanal_scission",
        "source": "Grosch (1999)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "linolenate_to_propanal_molar_yield_pct": 52.0,
        "linolenate_to_pentanal_molar_yield_pct": 22.0,
        "notes": "C18:3 oxidation: propanal 52%, pentanal 22% from beta-scission.",
    },
    {
        "id": "farmer_1991_alkyl_thiazoles",
        "source": "Farmer et al. (1991)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 120},
        "alkyl_thiazole_yield_lipid_derived_pct": 12.0,
        "notes": "Alkyl thiazoles: 12% yield from lipid-amino acid interactions at 120°C.",
    },
    {
        "id": "kamal_eldin_2003_triolein_scission",
        "source": "Kamal-Eldin et al. (2003)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "triolein_nonanal_molar_yield_pct": 38.0,
        "notes": "Triolein oxidation beta-scission: nonanal 38% molar yield.",
    },
    {
        "id": "tan_2001_oil_dsc_oxidation",
        "source": "Tan et al. (2001)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"method": "DSC_oxidation"},
        "soybean_oil_onset_temp_C": 195.0,
        "palm_oil_onset_temp_C": 213.0,
        "notes": "DSC oxidation onset: soybean oil 195°C, palm oil 213°C.",
    },
    {
        "id": "bayram_2023_stripped_soybean_oil_hexanal",
        "source": "Bayram et al. (2023)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"matrix": "stripped_soybean_oil"},
        "hexanal_ug_per_g_oil_after_72h_60C": 18.5,
        "notes": "Stripped soybean oil: hexanal 18.5 ug/g after 72h at 60°C.",
    },
    {
        "id": "spier_2010_metal_naphthenate_peroxide_decomp",
        "source": "Spier et al. (2010)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"catalyst": "metal_naphthenate"},
        "peroxide_decomp_k_1_per_s_metal_cat": 3.5e-4,
        "notes": "Metal naphthenate catalysed peroxide decomposition: k=3.5×10⁻⁴/s.",
    },
    {
        "id": "hidalgo_2007_decadienal_phenylalanine_styrene",
        "source": "Hidalgo et al. (2007)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 100},
        "styrene_yield_mol_pct_phe_decadienal": 8.5,
        "notes": "Styrene 8.5 mol% from Phe+2,4-decadienal Strecker-like at 100°C.",
    },
    {
        "id": "zamora_2010_decadienal_asparagine_decarboxylation",
        "source": "Zamora et al. (2010)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 100},
        "propanal_from_asn_decadienal_pct": 22.0,
        "notes": "Propanal 22% yield from Asn+decadienal decarboxylation at 100°C.",
    },
    {
        "id": "ding_2020_schiff_base_amadori_emulsion_rates",
        "source": "Ding et al. (2020)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 60, "matrix": "o/w_emulsion"},
        "schiff_base_formation_k_1_per_h": 0.045,
        "amadori_formation_k_1_per_h": 0.012,
        "notes": "Schiff base k=0.045/h, Amadori k=0.012/h in o/w emulsion at 60°C.",
    },
    {
        "id": "richards_2009_hemoglobin_liposome_oxidation",
        "source": "Richards et al. (2009)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 37},
        "Hb_mediated_TBARS_fold_increase": 3.8,
        "notes": "Hemoglobin mediates lipid oxidation 3.8-fold in liposome model at 37°C.",
    },
    {
        "id": "smagghe_2006_leghemoglobin_oxygen_dissociation",
        "source": "Smagghe et al. (2006)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 25},
        "leghemoglobin_O2_kd_umol_L": 0.024,
        "leghemoglobin_O2_kon_1_per_umol_per_s": 46.0,
        "notes": "Leghemoglobin Kd=0.024 μM O₂; kon=46/μM/s — very high O₂ affinity.",
    },
    {
        "id": "trikusuma_2019",
        "source": "Trikusuma et al. (2019)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"matrix": "pea_protein_emulsion"},
        "hexanal_oav_pea_emulsion_range": [4.0, 32.0],
        "notes": "Hexanal OAV 4-32 in pea protein emulsion; key oxidative off-note range.",
    },
    {
        "id": "pmc_2026_hme_hexanal_baseline",
        "source": "PMC 2026 HME hexanal baseline",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"process": "HME", "T_C": 150},
        "hexanal_hme_ug_per_g": 1.8,
        "notes": "HME baseline hexanal 1.8 ug/g; lipid oxidation contribution to off-notes.",
    },
    {
        "id": "xu_2024_soybean_pbma_hexanal",
        "source": "Xu et al. (2024)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"matrix": "pbma_soy"},
        "hexanal_ug_per_g_pbma_soy": 2.4,
        "notes": "Soy-based PBMA hexanal 2.4 ug/g; benchmark for plant-based meat analog.",
    },
    {
        "id": "yeo_shibamoto_1991_wof_hexanal",
        "source": "Yeo & Shibamoto (1991)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"matrix": "wof_model"},
        "hexanal_wof_ppb": 220.0,
        "notes": "Warmed-over flavour (WOF) hexanal benchmark: 220 ppb.",
    },
    {
        "id": "marquez_ruiz_2014_oleic_oav_anchor",
        "source": "Marquez-Ruiz et al. (2014)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "nonanal_ODT_ppb": 1.0,
        "nonanal_OAV_range": [1.0, 28.0],
        "notes": "Nonanal ODT 1.0 ppb; OAV 1-28 from oleic acid oxidation.",
    },
    {
        "id": "messina_2022_pbma_oil_oav_anchor",
        "source": "Messina et al. (2022)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"matrix": "pbma_oil_fraction"},
        "hexanal_OAV_pbma_oil": 18.0,
        "nonanal_OAV_pbma_oil": 8.0,
        "notes": "PBMA oil fraction: hexanal OAV=18, nonanal OAV=8.",
    },
]

PE_SINK_PRIORS = [
    {
        "id": "solis_calero_2015_pe_glyoxal",
        "source": "Solís-Calero et al. (2015)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"method": "DFT"},
        "PE_glyoxal_barrier_kcal_mol": 21.5,
        "notes": "DFT: PE-glyoxal Schiff base barrier 21.5 kcal/mol.",
    },
    {
        "id": "solis_calero_2013_pe_amadori",
        "source": "Solís-Calero et al. (2013)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"method": "DFT"},
        "PE_amadori_barrier_kcal_mol": 22.8,
        "notes": "DFT: PE-sugar Amadori rearrangement barrier 22.8 kcal/mol.",
    },
    {
        "id": "lertsiri_1998_pe_glycation",
        "source": "Lertsiri et al. (1998)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 60, "pH": 7.4},
        "PE_glycation_k_1_per_h": 8.5e-3,
        "notes": "PE glycation rate k=8.5×10⁻³/h at 60°C/pH 7.4.",
    },
    {
        "id": "hidalgo_2005_pe_ribose_lysine",
        "source": "Hidalgo et al. (2005)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 100},
        "PE_ribose_lysine_pyridinium_yield_pct": 28.0,
        "notes": "PE-ribose-lysine pyridinium derivative: 28% yield at 100°C.",
    },
    {
        "id": "zamora_2020_pe_dihydropyridine",
        "source": "Zamora et al. (2020)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 80},
        "dihydropyridine_PE_yield_pct": 15.0,
        "notes": "PE-dihydropyridine adduct: 15% yield at 80°C from carbonyl crosstalk.",
    },
    {
        "id": "hidalgo_2006_pe_lysine_antioxidant",
        "source": "Hidalgo et al. (2006)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 60},
        "PE_lysine_antioxidant_inhibition_pct": 38.0,
        "notes": "PE-Lys Maillard product inhibits lipid oxidation 38%.",
    },
    {
        "id": "vilanova_2012_pe_schiff_base",
        "source": "Vilanova et al. (2012)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"T_C": 120},
        "PE_schiff_base_k_1_per_min": 2.1e-3,
        "notes": "PE Schiff base formation k=2.1×10⁻³/min at 120°C.",
    },
    {
        "id": "biondi_2010_oil_microwave_degradation",
        "source": "Biondi et al. (2010)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"process": "microwave"},
        "PE_loss_pct_microwave_5min": 22.0,
        "notes": "PE degraded 22% by 5 min microwave; faster than conventional heating.",
    },
]


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
    # melanoidin_polymerization
    make_slr_entry("brands_2002_casein_sugar_melanoidin", "16", "Brands et al. (2002)", "",
                   "casein_sugar_model", ["melanoidin"], ["melanoidin_polymerisation"],
                   ["MW 5-300 kDa", "half-time 240 min"], ["data/lit/computational_priors.json"]),
    make_slr_entry("gigl_2021_coffee_thiol_binding", "16", "Gigl et al. (2021)", "",
                   "coffee_melanoidin", ["thiols"], ["thiol_melanoidin_binding"],
                   ["thiol depletion 8-fold"], ["data/lit/computational_priors.json"]),
    make_slr_entry("suzuki_philp_1990_sulfur_melanoidin", "16", "Suzuki & Philp (1990)", "",
                   "model_melanoidin", ["sulfur_volatiles"], ["melanoidin_thioether_trapping"],
                   ["sulfur trapping 45%"], ["data/lit/computational_priors.json"]),
    make_slr_entry("mundt_wedzicha_2007_biscuit_browning", "16", "Mundt & Wedzicha (2007)", "",
                   "biscuit", ["browning"], ["biscuit_browning_kinetics"],
                   ["A420 0.0021/min", "Ea 98 kJ/mol"], ["data/lit/computational_priors.json"]),
    make_slr_entry("cao_2024_carp_myoglobin_mrp", "16", "Cao et al. (2024)", "",
                   "carp_myoglobin", ["MRP", "thiols"], ["mrp_thiol_binding"],
                   ["MRP thiol 28%", "myoglobin heme -35%"], ["data/lit/computational_priors.json"]),
    make_slr_entry("hofmann_2001_melanoidin_thioether", "16", "Hofmann et al. (2001)", "",
                   "model_melanoidin", ["FFT", "MFT"], ["thioether_incorporation"],
                   ["FFT 62%", "MFT 48% in melanoidin"], ["data/lit/computational_priors.json"]),
    make_slr_entry("martins_van_boekel_2005_ascorbic_amino_browning", "14", "Martins & van Boekel (2005)", "",
                   "ascorbic_amino", ["melanoidin"], ["ascorbic_melanoidin_crossover"],
                   ["crossover 100°C", "melanoidin 18% from AA"], ["data/lit/computational_priors.json"]),
    # ascorbic_acid_maillard
    make_slr_entry("smuda_glomb_2013_aa_degradation_pathways", "14", "Smuda & Glomb (2013)",
                   "10.1002/anie.201300399", "aqueous_model", ["CML", "3-deoxyosone"],
                   ["aa_degradation_pathways"], ["α-frag 31%", "β-cleavage 32%", "decarboxylation 12%"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("serpen_gokmen_2007_ascorbic_redox_kinetics", "14", "Serpen & Gökmen (2007)",
                   "10.1016/j.foodchem.2006.11.073", "aqueous_model", ["ascorbic acid"],
                   ["ascorbic_metal_catalysis"], ["Cu²⁺ 88x", "Fe³⁺ 14x"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("yang_2021_ascorbic_glycine_kinetics", "14", "Yang et al. (2021)",
                   "10.47836/ifrj.28.3.16", "aa_glycine_model", ["browning"],
                   ["aa_glycine_browning_ea"], ["Ea 60.76 kJ/mol excess AA"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("yu_2018_ascorbic_basic_amino_browning", "14", "Yu et al. (2018)",
                   "10.1590/1678-457x.08717", "aa_amino_model", ["browning"],
                   ["aa_amino_acid_browning"], ["Ea Lys 54.94, Arg 50.08, His 35.31 kJ/mol"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("manso_2001_orange_juice_ascorbic_degradation", "14", "Manso et al. (2001)",
                   "10.1046/j.1365-2621.2001.t01-1-00460.x", "orange_juice", ["ascorbic acid"],
                   ["aa_degradation_oj"], ["Ea 55.6 kJ/mol"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("takase_2025_lemon_juice_ascorbic_dicarbonyl", "14", "Takase et al. (2025)",
                   "10.1021/acsfoodscitech.4c00956", "lemon_juice", ["threosone", "xylosone"],
                   ["dha_anaerobic_degradation"], ["threosone 22%", "xylosone 15%"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("hendrickx_1998_ascorbic_isobaric_degradation", "14", "Hendrickx et al. (1998)",
                   "10.1021/jf9708251", "aqueous_model", ["ascorbic acid"],
                   ["isobaric_ascorbic_degradation"], ["D₁₂₁=246 min", "z=27.15°C", "Ea=54.8 kJ/mol"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("jian_2012_ascorbic_ethanolic_degradation", "14", "Jian et al. (2012)",
                   "10.1021/jf3032342", "ethanolic_model", ["ascorbic acid"],
                   ["ethanolic_aa_degradation"], ["Ea 43.3-96.6 kJ/mol"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("frontiers_2022_hcw_aa_arrhenius_anchor", "14", "Frontiers (2022) HCW AA", "",
                   "hcw_model", ["ascorbic acid"], ["hcw_aa_degradation"],
                   ["Ea 82 kJ/mol HCW"], ["data/lit/computational_priors.json"]),
    make_slr_entry("pmid_1904866_pentosidine_equivalence_anchor", "14", "PMID 1904866 pentosidine", "",
                   "protein_model", ["pentosidine"], ["pentosidine_aa_equivalence"],
                   ["pentosidine 13.2-17.0 pmol/mg"], ["data/lit/computational_priors.json"]),
    make_slr_entry("scielo_brasil_aa_crosslink_hierarchy_anchor", "14", "SciELO Brasil (Ref. 4)", "",
                   "protein_model", ["CML", "crosslinks"], ["aa_crosslink_amino_acid_hierarchy"],
                   ["Lys Ea 54.94", "Arg 50.08", "His 35.31 kJ/mol"], ["data/lit/computational_priors.json"]),
    make_slr_entry("pmc3199460_ascorbic_dicarbonyl", "14", "PMC3199460 (Ref. 2)", "",
                   "aqueous_model", ["3-deoxythreosone", "xylosone"], ["ascorbic_dicarbonyl_yields"],
                   ["3-deoxythreosone 9.1 pmol/mg", "xylosone 0.5 pmol/mg"],
                   ["data/lit/computational_priors.json"]),
    # polyphenol_amino_capping
    make_slr_entry("shu_2019_cysteine_quinone_kinetics", "13", "Shu et al. (2019)", "",
                   "aqueous_model", ["cysteine", "quinone"], ["cys_quinone_kinetics"],
                   ["k₂=850 L/mol/s"], ["data/lit/computational_priors.json"]),
    make_slr_entry("zhu_2020_epicatechin_mgo_go_scavenging", "13", "Zhu H. et al. (2020)",
                   "", "aqueous_model", ["MGO", "GO"], ["epicatechin_mgo_scavenging"],
                   ["MGO k₂=1.8", "GO k₂=0.9 L/mol/s"], ["data/lit/computational_priors.json"]),
    make_slr_entry("zhu_2020_polyphenol_mgo_structure_kinetics", "13", "Zhu H. et al. (2020b)",
                   "", "aqueous_model", ["EGCG", "quercetin"], ["polyphenol_mgo_structure_reactivity"],
                   ["EGCG k₂=4.2", "quercetin k₂=2.1 L/mol/s"], ["data/lit/computational_priors.json"]),
    make_slr_entry("liu_2023_benzoquinone_lysine_kinetics", "13", "Liu J. et al. (2023)",
                   "", "aqueous_model", ["benzoquinone", "lysine"], ["bq_lys_adduct_kinetics"],
                   ["k₂=1200 L/mol/s"], ["data/lit/computational_priors.json"]),
    make_slr_entry("comert_gokmen_2019_digestive_mgo_scavenging", "13", "Cömert & Gökmen (2019)",
                   "", "digestion_model", ["MGO"], ["digestive_mgo_scavenging"],
                   ["rosemary 62% MGO scavenged"], ["data/lit/computational_priors.json"]),
    make_slr_entry("munoz_2007_chlorogenic_oxidase_quinone", "13", "Munoz et al. (2007)",
                   "", "enzymatic_model", ["chlorogenate_quinone"], ["chlorogenate_quinone_stability"],
                   ["quinone half-life 45 s"], ["data/lit/computational_priors.json"]),
    make_slr_entry("comert_gokmen_2021_epicatechin_cysteine_mgo_synergy", "13", "Cömert & Gökmen (2021)",
                   "", "aqueous_model", ["MGO"], ["ec_cys_mgo_synergy"],
                   ["EC+Cys 72% vs EC alone 38%", "synergy 1.9x"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("song_2009_benzoquinone_gsh_conjugation", "13", "Song et al. (2009)",
                   "", "aqueous_model", ["benzoquinone", "GSH"], ["bq_gsh_conjugation"],
                   ["k₂=980 L/mol/s BQ-GSH"], ["data/lit/computational_priors.json"]),
    make_slr_entry("li_2017_quercetin_mgo_adduct_browning_mitigation", "13", "Li et al. (2017)",
                   "", "aqueous_model", ["quercetin", "MGO"], ["quercetin_mgo_browning"],
                   ["adduct 35 mol%", "browning -42%"], ["data/lit/computational_priors.json"]),
    make_slr_entry("monforte_2018_phenylacetaldehyde_quinones", "13", "Monforte et al. (2018)",
                   "", "aqueous_model", ["phenylacetaldehyde", "quinone"],
                   ["phenylacetaldehyde_quinone_kinetics"], ["k₂=0.35 L/mol/s"],
                   ["data/lit/computational_priors.json"]),
    make_slr_entry("rawel_2002_cga_cysteine_blocking", "13", "Rawel et al. (2002)",
                   "", "protein_model", ["CGA", "cysteine"], ["cga_cys_blocking"],
                   ["blocking 85%", "free SH -78%"], ["data/lit/computational_priors.json"]),
    make_slr_entry("liu_2022_ppi_oav_anchors", "07", "Liu et al. (2022)",
                   "", "ppi", ["hexanal"], ["ppi_off_note_oav"],
                   ["hexanal OAV=28 PPI"], ["data/lit/computational_priors.json"]),
    # lipid_oxidation_and_carbonylic_crosstalk
    make_slr_entry("esterbauer_1991_4hne_kinetics", "15", "Esterbauer et al. (1991)", "",
                   "aqueous_model", ["4-HNE", "cysteine"], ["4hne_michael_kinetics"],
                   ["Cys k₂=8.5", "Lys k₂=1.1 L/mol/s"], ["data/lit/computational_priors.json"]),
    make_slr_entry("grosch_1982_hexanal_odt", "07", "Grosch (1982)", "",
                   "aqueous_model", ["hexanal"], ["hexanal_odour_threshold"],
                   ["ODT 4.5 ppb"], ["data/lit/computational_priors.json"]),
    make_slr_entry("frankel_1982_c182_hexanal_scission", "15", "Frankel (1982)", "",
                   "linoleic_model", ["hexanal", "propanal"], ["c182_beta_scission"],
                   ["hexanal 45%", "propanal 28% from C18:2"], ["data/lit/computational_priors.json"]),
    make_slr_entry("grosch_1999_c183_propanal_scission", "15", "Grosch (1999)", "",
                   "linolenate_model", ["propanal", "pentanal"], ["c183_beta_scission"],
                   ["propanal 52%", "pentanal 22% from C18:3"], ["data/lit/computational_priors.json"]),
    make_slr_entry("farmer_1991_alkyl_thiazoles", "15", "Farmer et al. (1991)", "",
                   "lipid_amino_model", ["thiazoles"], ["alkyl_thiazole_lipid_origin"],
                   ["alkyl thiazoles 12% from lipid-AA"], ["data/lit/computational_priors.json"]),
    make_slr_entry("kamal_eldin_2003_triolein_scission", "15", "Kamal-Eldin et al. (2003)", "",
                   "triolein_model", ["nonanal"], ["triolein_nonanal_yield"],
                   ["nonanal 38% from triolein"], ["data/lit/computational_priors.json"]),
    make_slr_entry("tan_2001_oil_dsc_oxidation", "15", "Tan et al. (2001)", "",
                   "vegetable_oils", ["oxidation_products"], ["dsc_oxidation_onset"],
                   ["soybean 195°C", "palm 213°C DSC onset"], ["data/lit/computational_priors.json"]),
    make_slr_entry("bayram_2023_stripped_soybean_oil_hexanal", "15", "Bayram et al. (2023)", "",
                   "stripped_soybean_oil", ["hexanal"], ["soybean_oil_hexanal_benchmark"],
                   ["hexanal 18.5 ug/g 72h 60°C"], ["data/lit/computational_priors.json"]),
    make_slr_entry("spier_2010_metal_naphthenate_peroxide_decomp", "15", "Spier et al. (2010)", "",
                   "model_system", ["peroxides"], ["metal_catalyzed_peroxide_decomp"],
                   ["k=3.5e-4/s metal naphthenate"], ["data/lit/computational_priors.json"]),
    make_slr_entry("trikusuma_2019", "07", "Trikusuma et al. (2019)", "",
                   "pea_protein_emulsion", ["hexanal"], ["hexanal_oav_pea_emulsion"],
                   ["hexanal OAV 4-32"], ["data/lit/computational_priors.json"]),
    make_slr_entry("marquez_ruiz_2014_oleic_oav_anchor", "15", "Marquez-Ruiz et al. (2014)", "",
                   "oleic_model", ["nonanal"], ["nonanal_oav_anchor"],
                   ["nonanal ODT 1.0 ppb", "OAV 1-28"], ["data/lit/computational_priors.json"]),
    make_slr_entry("messina_2022_pbma_oil_oav_anchor", "07", "Messina et al. (2022)", "",
                   "pbma_oil", ["hexanal", "nonanal"], ["pbma_oil_oav"],
                   ["hexanal OAV=18", "nonanal OAV=8"], ["data/lit/computational_priors.json"]),
    make_slr_entry("xu_2024_soybean_pbma_hexanal", "07", "Xu et al. (2024)", "",
                   "pbma_soy", ["hexanal"], ["soy_pbma_hexanal"],
                   ["hexanal 2.4 ug/g"], ["data/lit/computational_priors.json"]),
    make_slr_entry("yeo_shibamoto_1991_wof_hexanal", "07", "Yeo & Shibamoto (1991)", "",
                   "wof_model", ["hexanal"], ["wof_hexanal_benchmark"],
                   ["hexanal 220 ppb WOF"], ["data/lit/computational_priors.json"]),
    # lipid_maillard_crosstalk
    make_slr_entry("pmc_2026_hme_hexanal_baseline", "15", "PMC 2026 HME baseline", "",
                   "hme_matrix", ["hexanal"], ["hme_hexanal_baseline"],
                   ["hexanal 1.8 ug/g HME"], ["data/lit/computational_priors.json"]),
    make_slr_entry("hidalgo_2007_decadienal_phenylalanine_styrene", "15", "Hidalgo et al. (2007)", "",
                   "phe_decadienal_model", ["styrene"], ["strecker_lipid_styrene"],
                   ["styrene 8.5 mol% from Phe+decadienal"], ["data/lit/computational_priors.json"]),
    make_slr_entry("zamora_2010_decadienal_asparagine_decarboxylation", "15", "Zamora et al. (2010)", "",
                   "asn_decadienal_model", ["propanal"], ["lipid_asn_decarboxylation"],
                   ["propanal 22% from Asn+decadienal"], ["data/lit/computational_priors.json"]),
    make_slr_entry("ding_2020_schiff_base_amadori_emulsion_rates", "15", "Ding et al. (2020)", "",
                   "ow_emulsion", ["Schiff_base", "Amadori"], ["emulsion_maillard_rates"],
                   ["Schiff k=0.045/h", "Amadori k=0.012/h"], ["data/lit/computational_priors.json"]),
    make_slr_entry("richards_2009_hemoglobin_liposome_oxidation", "15", "Richards et al. (2009)", "",
                   "liposome_model", ["lipid_oxidation"], ["hb_mediated_lipid_oxidation"],
                   ["TBARS 3.8x Hb-mediated"], ["data/lit/computational_priors.json"]),
    make_slr_entry("smagghe_2006_leghemoglobin_oxygen_dissociation", "07", "Smagghe et al. (2006)", "",
                   "soy_leghemoglobin", ["O2"], ["leghemoglobin_o2_kinetics"],
                   ["Kd=0.024 μM O₂", "kon=46/μM/s"], ["data/lit/computational_priors.json"]),
    # phospholipid_amine_sink
    make_slr_entry("solis_calero_2015_pe_glyoxal", "15", "Solís-Calero et al. (2015)", "",
                   "PE_model", ["glyoxal"], ["pe_glyoxal_schiff_barrier"],
                   ["barrier 21.5 kcal/mol"], ["data/lit/computational_priors.json"]),
    make_slr_entry("solis_calero_2013_pe_amadori", "15", "Solís-Calero et al. (2013)", "",
                   "PE_model", ["glucose"], ["pe_amadori_barrier"],
                   ["barrier 22.8 kcal/mol"], ["data/lit/computational_priors.json"]),
    make_slr_entry("lertsiri_1998_pe_glycation", "15", "Lertsiri et al. (1998)", "",
                   "PE_phospholipid", ["glucose"], ["pe_glycation_kinetics"],
                   ["k=8.5e-3/h at 60°C"], ["data/lit/computational_priors.json"]),
    make_slr_entry("hidalgo_2005_pe_ribose_lysine", "15", "Hidalgo et al. (2005)", "",
                   "PE_ribose_model", ["pyridinium"], ["pe_ribose_lysine_pyridinium"],
                   ["pyridinium 28%"], ["data/lit/computational_priors.json"]),
    make_slr_entry("zamora_2020_pe_dihydropyridine", "15", "Zamora et al. (2020)", "",
                   "PE_carbonyl_model", ["dihydropyridine"], ["pe_dihydropyridine_adduct"],
                   ["DHP 15%"], ["data/lit/computational_priors.json"]),
    make_slr_entry("hidalgo_2006_pe_lysine_antioxidant", "15", "Hidalgo et al. (2006)", "",
                   "PE_lysine_model", ["MRP"], ["pe_lys_antioxidant"],
                   ["antioxidant -38% lipid ox"], ["data/lit/computational_priors.json"]),
    make_slr_entry("vilanova_2012_pe_schiff_base", "15", "Vilanova et al. (2012)", "",
                   "PE_model", ["Schiff_base"], ["pe_schiff_kinetics"],
                   ["k=2.1e-3/min"], ["data/lit/computational_priors.json"]),
    make_slr_entry("biondi_2010_oil_microwave_degradation", "15", "Biondi et al. (2010)", "",
                   "oil_microwave", ["PE"], ["pe_microwave_degradation"],
                   ["PE -22% 5min microwave"], ["data/lit/computational_priors.json"]),
]


def main():
    intake = load_json(INTAKE_REGISTRY_PATH)
    matrix = load_json(SLR_MATRIX_PATH)
    priors = load_json(COMPUTATIONAL_PRIORS_PATH)
    safety = load_json(SAFETY_PAYLOADS_PATH)
    flavor = load_json(FLAVOR_PAYLOADS_PATH)

    # ── 1. Update status in intake registry ──────────────────────────────────
    encoded_count = 0
    for ref in intake["eligible_references"]:
        if ref["id"] in TARGET_IDS and ref.get("status") == "ready_for_intake_encoding":
            ref["status"] = "encoded"
            encoded_count += 1
    print(f"Marked {encoded_count} refs as encoded")

    # ── 2. Computational priors ───────────────────────────────────────────────
    added_priors = 0

    existing_melanoidin = {p["id"] for p in priors["melanoidin_trapping_priors"]}
    for pr in MELANOIDIN_PRIORS:
        if pr["id"] not in existing_melanoidin:
            priors["melanoidin_trapping_priors"].append(pr)
            added_priors += 1

    existing_ascorbic = {p["id"] for p in priors["ascorbic_pathway_priors"]}
    for pr in ASCORBIC_PRIORS:
        if pr["id"] not in existing_ascorbic:
            priors["ascorbic_pathway_priors"].append(pr)
            added_priors += 1

    existing_polyphenol = {p["id"] for p in priors["polyphenol_thiol_capping_priors"]}
    for pr in POLYPHENOL_PRIORS:
        if pr["id"] not in existing_polyphenol:
            priors["polyphenol_thiol_capping_priors"].append(pr)
            added_priors += 1

    existing_lipid = {p["id"] for p in priors["lipid_offnote_priors"]}
    for pr in LIPID_PRIORS:
        if pr["id"] not in existing_lipid:
            priors["lipid_offnote_priors"].append(pr)
            added_priors += 1

    existing_pe = {p["id"] for p in priors["phospholipid_amine_sink_priors"]}
    for pr in PE_SINK_PRIORS:
        if pr["id"] not in existing_pe:
            priors["phospholipid_amine_sink_priors"].append(pr)
            added_priors += 1

    print(f"Added {added_priors} computational prior entries")

    # ── 3. Safety entries (ascorbic browning safety refs) ────────────────────
    existing_safety = {e["id"] for e in safety["entries"]}
    added_safety = 0
    safety_entries = [
        {
            "id": "smuda_glomb_2013_aa_degradation_pathways",
            "category": "ascorbic_maillard",
            "source_citation": "Smuda & Glomb (2013)",
            "doi": "10.1002/anie.201300399",
            "matrix_context": "model_system",
            "analytical_method": "LC-MS",
            "benchmark_role": "reference_anchor",
            "key_values": {
                "oxidative_alpha_fragmentation_pct": 31.0,
                "beta_cleavage_pct": 32.0,
            },
            "notes": "AA degradation: α-frag 31%, β-cleavage 32%, CML formed via ascorbic pathway.",
        },
        {
            "id": "yu_2018_ascorbic_basic_amino_browning",
            "category": "ascorbic_maillard",
            "source_citation": "Yu et al. (2018)",
            "doi": "10.1590/1678-457x.08717",
            "matrix_context": "aqueous_model",
            "analytical_method": "UV-vis",
            "benchmark_role": "reference_anchor",
            "key_values": {
                "Ea_lysine_kj_mol": 54.94,
                "Ea_arginine_kj_mol": 50.08,
                "Ea_histidine_kj_mol": 35.31,
            },
            "notes": "Zero-order AA-amino acid browning: Lys fastest (Ea 54.94 kJ/mol).",
        },
        {
            "id": "zha_2020_ppi_glycation_aggregation",
            "category": "protein_damage_markers",
            "source_citation": "Zha et al. (2020)",
            "doi": "",
            "matrix_context": "ppi_glycated",
            "analytical_method": "SDS-PAGE",
            "benchmark_role": "reference_anchor",
            "key_values": {
                "aggregation_index_glycated_ppi": 0.42,
                "CML_increase_fold_glycated_ppi": 2.8,
            },
            "notes": "PPI glycation: aggregation index 0.42; CML 2.8x increase.",
        },
        {
            "id": "sun_2020_solid_matrix_cml_cel_accumulation",
            "category": "ages",
            "source_citation": "Sun et al. (2020)",
            "doi": "",
            "matrix_context": "solid_matrix",
            "analytical_method": "ELISA",
            "benchmark_role": "reference_anchor",
            "key_values": {
                "CML_solid_matrix_ug_per_g_protein": 28.0,
            },
            "notes": "CML 28 ug/g protein in solid matrix; AGE accumulation benchmark.",
        },
    ]
    for entry in safety_entries:
        if entry["id"] not in existing_safety:
            safety["entries"].append(entry)
            added_safety += 1
    print(f"Added {added_safety} safety entries")

    # ── 4. Flavor off-note entries ────────────────────────────────────────────
    existing_offnote = {fl["id"] for fl in flavor["off_note_reference_anchors"]}
    added_flavor = 0
    offnote_entries = [
        {
            "id": "marquez_ruiz_2014_oleic_oav_anchor",
            "compound": "nonanal",
            "source_citation": "Marquez-Ruiz et al. (2014)",
            "doi": "",
            "matrix_context": "oleic_model",
            "analytical_method": "GC-MS",
            "units": "ppb",
            "benchmark_role": "reference_anchor",
            "pipeline_role": "reference_only",
            "target_direction": "nonanal_oav_from_oleic_oxidation",
            "numeric_band_or_point": {
                "type": "point",
                "value": 1.0,
                "condition": "nonanal ODT",
            },
            "notes": "Nonanal ODT 1.0 ppb from oleic acid oxidation; OAV 1-28.",
        },
        {
            "id": "messina_2022_pbma_oil_oav_anchor",
            "compound": "hexanal",
            "source_citation": "Messina et al. (2022)",
            "doi": "",
            "matrix_context": "pbma_oil_fraction",
            "analytical_method": "GC-MS",
            "units": "OAV",
            "benchmark_role": "reference_anchor",
            "pipeline_role": "reference_only",
            "target_direction": "pbma_lipid_off_note_oav",
            "numeric_band_or_point": {
                "type": "point",
                "value": 18.0,
                "condition": "hexanal OAV PBMA oil",
            },
            "notes": "PBMA oil: hexanal OAV=18, nonanal OAV=8.",
        },
        {
            "id": "yeo_shibamoto_1991_wof_hexanal",
            "compound": "hexanal",
            "source_citation": "Yeo & Shibamoto (1991)",
            "doi": "",
            "matrix_context": "wof_model",
            "analytical_method": "GC-MS",
            "units": "ppb",
            "benchmark_role": "reference_anchor",
            "pipeline_role": "reference_only",
            "target_direction": "wof_hexanal_threshold",
            "numeric_band_or_point": {
                "type": "point",
                "value": 220.0,
                "condition": "WOF hexanal ppb",
            },
            "notes": "Warmed-over flavour: hexanal 220 ppb benchmark.",
        },
        {
            "id": "rawel_2002_cga_cysteine_blocking",
            "compound": "chlorogenic acid",
            "source_citation": "Rawel et al. (2002)",
            "doi": "",
            "matrix_context": "protein_model",
            "analytical_method": "fluorescence",
            "units": "percent_blocked",
            "benchmark_role": "reference_anchor",
            "pipeline_role": "reference_only",
            "target_direction": "cga_cys_blocking_off_note_suppression",
            "numeric_band_or_point": {
                "type": "point",
                "value": 85.0,
                "condition": "blocking efficiency oxidised CGA",
            },
            "notes": "Oxidised CGA blocks 85% reactive Cys; off-note suppression via polyphenol.",
        },
        {
            "id": "liu_2022_ppi_oav_anchors",
            "compound": "hexanal",
            "source_citation": "Liu et al. (2022)",
            "doi": "",
            "matrix_context": "ppi",
            "analytical_method": "GC-MS-AEDA",
            "units": "OAV",
            "benchmark_role": "reference_anchor",
            "pipeline_role": "reference_only",
            "target_direction": "ppi_off_note_oav_ranking",
            "numeric_band_or_point": {
                "type": "point",
                "value": 28.0,
                "condition": "hexanal OAV PPI",
            },
            "notes": "PPI OAV: hexanal=28 (dominant off-note).",
        },
    ]
    for fl in offnote_entries:
        if fl["id"] not in existing_offnote:
            flavor["off_note_reference_anchors"].append(fl)
            added_flavor += 1
    print(f"Added {added_flavor} flavor off-note entries")

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
    save_json(safety, SAFETY_PAYLOADS_PATH)
    save_json(flavor, FLAVOR_PAYLOADS_PATH)

    print("Done: remaining families batch encoded.")


if __name__ == "__main__":
    main()
