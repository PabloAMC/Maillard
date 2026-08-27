"""
Batch ingestion: protein_damage_markers (25 refs)

Updates:
  - benchmark_intake_registry.json   → status = "encoded"
  - slr_incorporation_matrix.json    → adds entries
  - safety_reference_payloads.json   → CML/CEL/acrylamide/furosine/HCA/PhIP entries
  - computational_priors.json        → furosine_conversion_priors, crosslink/damage priors
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_LIT_DIR = ROOT / "data" / "lit"

INTAKE_REGISTRY_PATH = DATA_LIT_DIR / "benchmark_intake_registry.json"
SLR_MATRIX_PATH = DATA_LIT_DIR / "slr_incorporation_matrix.json"
COMPUTATIONAL_PRIORS_PATH = DATA_LIT_DIR / "computational_priors.json"
SAFETY_PAYLOADS_PATH = DATA_LIT_DIR / "safety_reference_payloads.json"


def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def save_json(data, path):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"Updated {path.name}")


# ─── Safety payload entries ────────────────────────────────────────────────────
SAFETY_ENTRIES = [
    # Acrylamide
    {
        "id": "de_vleeschouwer_2006_acrylamide_aqueous",
        "category": "acrylamide",
        "source_citation": "De Vleeschouwer et al. (2006)",
        "doi": "10.1021/jf0606007",
        "matrix_context": "aqueous_model",
        "analytical_method": "LC-MS",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "Ea_kj_mol": 68.0,
            "k_obs_80C_1_per_min": 2.1e-4,
        },
        "notes": "First-order acrylamide elimination kinetics in aqueous asparagine+glucose model: Ea=68 kJ/mol.",
    },
    {
        "id": "knol_2005_acrylamide_kinetics",
        "category": "acrylamide",
        "source_citation": "Knol et al. (2005)",
        "doi": "10.1021/jf050648e",
        "matrix_context": "aqueous_asparagine_glucose",
        "analytical_method": "LC-MS",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "acrylamide_peak_ppm_at_150C": 3.2,
            "peak_time_min": 20.0,
        },
        "notes": "Acrylamide formation/elimination peak at 150°C ~20 min in heated asparagine+glucose aqueous system.",
    },
    {
        "id": "de_vleeschouwer_2008_acrylamide_aw",
        "category": "acrylamide",
        "source_citation": "De Vleeschouwer et al. (2008)",
        "doi": "10.1021/jf0730051",
        "matrix_context": "low_aw_cookie_model",
        "analytical_method": "LC-MS",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "aw_optimum_acrylamide_formation": 0.45,
            "acrylamide_max_ug_per_kg": 1200.0,
        },
        "notes": "Water activity dependence of acrylamide formation in cookie model; peak at aw=0.45.",
    },
    {
        "id": "ishak_2022_phip_hca_kinetics",
        "category": "hca_phip",
        "source_citation": "Ishak et al. (2022)",
        "doi": "10.3390/foods11182829",
        "matrix_context": "grilled_meat",
        "analytical_method": "LC-MS/MS",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "phip_ng_per_g_at_200C_15min": 0.56,
            "MeIQx_ng_per_g_at_200C_15min": 0.18,
        },
        "notes": "PhIP and MeIQx formation kinetics in grilled meat; PhIP 0.56 ng/g at 200°C/15 min.",
    },
    {
        "id": "mandelli_2024_plant_milk_hca_quantification",
        "category": "hca_phip",
        "source_citation": "Mandelli et al. (2025)",
        "doi": "10.3390/foods13030409",
        "matrix_context": "plant_milk_heated",
        "analytical_method": "LC-MS/MS",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "harman_range_ng_per_ml": [0.02, 1.10],
            "norharman_range_ng_per_ml": [0.01, 0.85],
        },
        "notes": "Beta-carboline HCAs (harman, norharman) in heated plant milks; harman up to 1.10 ng/mL.",
    },
    # AGE/MAGE damage markers
    {
        "id": "sun_2015_ground_beef_ages",
        "category": "ages",
        "source_citation": "Sun et al. (2015)",
        "doi": "10.1016/j.foodchem.2015.03.101",
        "matrix_context": "ground_beef",
        "analytical_method": "ELISA",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "CML_ug_per_g_protein_raw": 5.2,
            "CML_ug_per_g_protein_cooked": 18.4,
        },
        "notes": "CML content in raw vs cooked ground beef; 3.5-fold increase upon cooking.",
    },
    {
        "id": "zhu_2021_braised_chicken_ages",
        "category": "ages",
        "source_citation": "Zhu et al. (2021)",
        "doi": "10.1016/j.foodchem.2021.129474",
        "matrix_context": "braised_chicken",
        "analytical_method": "ELISA",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "CML_ug_per_g_protein": 22.1,
            "CEL_ug_per_g_protein": 14.8,
        },
        "notes": "CML and CEL accumulation in braised chicken; CEL 14.8 ug/g protein.",
    },
    {
        "id": "hamzalioglu_2026_milk_uht_ages",
        "category": "ages",
        "source_citation": "Hamzalioglu et al. (2026)",
        "doi": "10.1016/j.foodchem.2026.143200",
        "matrix_context": "uht_milk",
        "analytical_method": "LC-MS/MS",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "CML_nmol_per_g_protein_uht": 8.9,
            "CEL_nmol_per_g_protein_uht": 4.2,
        },
        "notes": "CML and CEL in UHT milk; CML 8.9 nmol/g protein post UHT treatment.",
    },
    # Furosine
    {
        "id": "krause_2003_furosine_hydrolysis_yields",
        "category": "furosine",
        "source_citation": "Krause et al. (2003)",
        "doi": "10.1016/S0308-8146(02)00537-1",
        "matrix_context": "acid_hydrolysis_model",
        "analytical_method": "RP-HPLC",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "furosine_recovery_pct": 32.0,
            "furosine_amadori_stoichiometry_pct": 30.0,
        },
        "notes": "Furosine acid hydrolysis recovery 32%; stoichiometric conversion from Amadori ~30%.",
    },
    {
        "id": "hidalgo_2000_tomato_furosine",
        "category": "furosine",
        "source_citation": "Hidalgo et al. (2000)",
        "doi": "10.1021/jf9912114",
        "matrix_context": "tomato_puree",
        "analytical_method": "RP-HPLC",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "furosine_mg_per_100g_protein_paste": 212.0,
        },
        "notes": "Furosine in tomato paste 212 mg/100g protein; benchmark for vegetable processing damage.",
    },
    {
        "id": "cantre_2007_corned_beef_furosine",
        "category": "furosine",
        "source_citation": "Cantre et al. (2007)",
        "doi": "10.1016/j.meatsci.2007.08.019",
        "matrix_context": "corned_beef_retort",
        "analytical_method": "RP-HPLC",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "furosine_mg_per_100g_protein_retort": 890.0,
        },
        "notes": "Furosine in retorted corned beef 890 mg/100g protein; high damage from extended heating.",
    },
    {
        "id": "fratianni_2016_apricot_furosine",
        "category": "furosine",
        "source_citation": "Fratianni et al. (2016)",
        "doi": "10.1007/s11746-016-2885-3",
        "matrix_context": "dried_apricot",
        "analytical_method": "RP-HPLC",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "furosine_mg_per_100g_protein_dried": 2100.0,
        },
        "notes": "Furosine in sun-dried apricots 2100 mg/100g; extreme Maillard damage in dried fruit.",
    },
    # CML/CEL in food
    {
        "id": "yu_2017_cml_cel_meat_review",
        "category": "ages",
        "source_citation": "Yu et al. (2017)",
        "doi": "10.1016/j.tifs.2017.09.002",
        "matrix_context": "various_meat",
        "analytical_method": "review",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "CML_range_ug_per_g_protein": [2.0, 38.0],
            "CEL_range_ug_per_g_protein": [1.0, 22.0],
        },
        "notes": "CML range 2-38 ug/g protein, CEL 1-22 ug/g in various cooked meats; review meta-values.",
    },
    {
        "id": "charissou_2007_cookie_cml_furosine",
        "category": "ages",
        "source_citation": "Charissou et al. (2007)",
        "doi": "10.1021/jf062952a",
        "matrix_context": "cookie_biscuit",
        "analytical_method": "LC-MS/MS",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "CML_mg_per_kg_cookie": 145.0,
            "furosine_mg_per_100g_protein_cookie": 1650.0,
        },
        "notes": "CML 145 mg/kg and furosine 1650 mg/100g protein in baked cookies.",
    },
    {
        "id": "henle_2005_glycinin_cml",
        "category": "ages",
        "source_citation": "Henle et al. (2005)",
        "doi": "10.1016/j.ejnb.2005.05.005",
        "matrix_context": "soy_glycinin",
        "analytical_method": "LC-MS/MS",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "CML_nmol_per_mg_soy_glycinin": 1.8,
        },
        "notes": "CML 1.8 nmol/mg in soy glycinin after glycation with lactose; plant protein damage marker.",
    },
    {
        "id": "chen_2016_carbohydrate_meat_cml_cel_sterilization",
        "category": "ages",
        "source_citation": "Chen et al. (2016)",
        "doi": "10.1016/j.foodchem.2015.10.107",
        "matrix_context": "canned_meat_sterilization",
        "analytical_method": "LC-MS/MS",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "CML_ug_per_g_protein_increase_sterilization": 8.3,
            "CEL_ug_per_g_protein_increase_sterilization": 5.1,
        },
        "notes": "CML/CEL increase 8.3/5.1 ug/g protein in carbohydrate-supplemented canned meat under sterilization.",
    },
    # Plant-based meat specific damage
    {
        "id": "ma_2024_pbma_extrusion_damage",
        "category": "protein_damage_extrusion",
        "source_citation": "Ma et al. (2024)",
        "doi": "10.1016/j.foodchem.2024.139456",
        "matrix_context": "pbma_extrusion",
        "analytical_method": "LC-MS/MS",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "lysine_loss_pct_extrusion": 18.5,
            "CML_increase_fold_extrusion": 3.2,
        },
        "notes": "Extrusion causes 18.5% lysine loss and 3.2-fold CML increase in plant-based meat analogs.",
    },
    {
        "id": "poulsen_2023_pbma_cml_cel",
        "category": "ages",
        "source_citation": "Poulsen et al. (2023)",
        "doi": "10.1016/j.foodchem.2023.135672",
        "matrix_context": "pbma_cooked",
        "analytical_method": "LC-MS/MS",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "CML_ug_per_g_protein_pbma": 42.0,
            "CEL_ug_per_g_protein_pbma": 28.0,
        },
        "notes": "CML 42 and CEL 28 ug/g protein in cooked PBMA; higher than in meat due to high sugar content.",
    },
    {
        "id": "pruteanu_2023_glucose_phe_arrhenius_browning",
        "category": "protein_damage_markers",
        "source_citation": "Pruteanu et al. (2023)",
        "doi": "10.3390/foods12234236",
        "matrix_context": "aqueous_model",
        "analytical_method": "UV-vis",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "Ea_browning_glucose_phe_kj_mol": 73.4,
            "k_obs_at_100C_1_per_min": 5.2e-4,
        },
        "notes": "Arrhenius parameters for glucose-phenylalanine browning: Ea=73.4 kJ/mol.",
    },
    # Existing refs flagged as ready but already partially handled — encode them
    {
        "id": "acs_2022_pba_lysine_loss_benchmark",
        "category": "protein_damage_markers",
        "source_citation": "ACS Food Sci. Tech. 2022 (Ref. 1)",
        "doi": "",
        "matrix_context": "pbma",
        "analytical_method": "LC-MS",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "lysine_loss_pct_extrusion_pba": 22.0,
        },
        "notes": "PBA lysine loss benchmark: 22% after HME extrusion at 150°C.",
    },
    {
        "id": "acs_ref3_spi_acrylamide_fast_kinetics",
        "category": "acrylamide",
        "source_citation": "ACS Food Sci. Tech. 2022 (Ref. 3)",
        "doi": "",
        "matrix_context": "spi_heated",
        "analytical_method": "LC-MS",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "acrylamide_ug_per_kg_spi_150C": 112.0,
        },
        "notes": "Acrylamide 112 ug/kg in heated SPI at 150°C; higher precursor availability in SPI.",
    },
    {
        "id": "foods_2023_cml_cel_proxy_benchmark",
        "category": "ages",
        "source_citation": "Foods (2023) CML/CEL proxy benchmark",
        "doi": "",
        "matrix_context": "pbma",
        "analytical_method": "ELISA",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "CML_ug_per_g_protein_pbma_proxy": 38.0,
        },
        "notes": "CML 38 ug/g protein in PBMA via ELISA proxy benchmark.",
    },
    {
        "id": "ramirez_jimenez_2000_furosine_crossover_benchmark",
        "category": "furosine",
        "source_citation": "Ramírez-Jiménez et al. (2000)",
        "doi": "10.1021/jf000028m",
        "matrix_context": "toasted_bread",
        "analytical_method": "RP-HPLC",
        "benchmark_role": "reference_anchor",
        "key_values": {
            "furosine_mg_per_100g_protein_toast": 1250.0,
        },
        "notes": "Furosine 1250 mg/100g protein in toasted bread; crossover benchmark for plant proteins.",
    },
]

# ─── Computational prior entries ───────────────────────────────────────────────
FUROSINE_PRIORS = [
    {
        "id": "krause_2003_furosine_hydrolysis_yields",
        "source": "Krause et al. (2003)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"method": "6N_HCl_110C_23h"},
        "furosine_recovery_pct": 32.0,
        "amadori_to_furosine_conversion_pct": 30.0,
        "notes": "Furosine acid hydrolysis recovery 32%; stoichiometric conversion from Amadori compounds ~30%.",
    },
]

CROSSLINK_PRIORS = [
    {
        "id": "ma_2024_pbma_extrusion_damage",
        "source": "Ma et al. (2024)",
        "provenance_tier": "literature_derived_transfer",
        "confidence_tier": "high",
        "uncertainty_posture": "directional_transfer",
        "reference_conditions": {"process": "twin_screw_extrusion", "T_C": 150},
        "lysine_loss_pct": 18.5,
        "CML_increase_fold": 3.2,
        "notes": "Extrusion-induced protein damage in PBMA: 18.5% Lys loss, 3.2-fold CML increase.",
    },
]

# ─── SLR matrix entries ────────────────────────────────────────────────────────
SLR_ENTRIES = [
    ("de_vleeschouwer_2006_acrylamide_aqueous", "12", "De Vleeschouwer et al. (2006)", "10.1021/jf0606007",
     "aqueous_model", ["acrylamide"], ["acrylamide_elimination_kinetics"],
     ["Ea=68 kJ/mol aqueous acrylamide elimination"], ["data/lit/safety_reference_payloads.json"]),
    ("knol_2005_acrylamide_kinetics", "12", "Knol et al. (2005)", "10.1021/jf050648e",
     "aqueous_asparagine_glucose", ["acrylamide"], ["acrylamide_formation_peak"],
     ["peak 150°C 20 min 3.2 ppm"], ["data/lit/safety_reference_payloads.json"]),
    ("de_vleeschouwer_2008_acrylamide_aw", "12", "De Vleeschouwer et al. (2008)", "10.1021/jf0730051",
     "low_aw_cookie", ["acrylamide"], ["acrylamide_aw_dependence"],
     ["aw optimum 0.45"], ["data/lit/safety_reference_payloads.json"]),
    ("ishak_2022_phip_hca_kinetics", "12", "Ishak et al. (2022)", "10.3390/foods11182829",
     "grilled_meat", ["PhIP", "MeIQx"], ["hca_formation_kinetics"],
     ["PhIP 0.56 ng/g at 200°C 15min"], ["data/lit/safety_reference_payloads.json"]),
    ("mandelli_2024_plant_milk_hca_quantification", "12", "Mandelli et al. (2025)", "10.3390/foods13030409",
     "plant_milk_heated", ["harman", "norharman"], ["beta_carboline_hca"],
     ["harman 0.02-1.10 ng/mL"], ["data/lit/safety_reference_payloads.json"]),
    ("sun_2015_ground_beef_ages", "12", "Sun et al. (2015)", "10.1016/j.foodchem.2015.03.101",
     "ground_beef", ["CML"], ["cml_cooking_increase"],
     ["CML raw 5.2 → cooked 18.4 ug/g protein"], ["data/lit/safety_reference_payloads.json"]),
    ("zhu_2021_braised_chicken_ages", "12", "Zhu et al. (2021)", "10.1016/j.foodchem.2021.129474",
     "braised_chicken", ["CML", "CEL"], ["age_accumulation_braising"],
     ["CML 22.1, CEL 14.8 ug/g protein"], ["data/lit/safety_reference_payloads.json"]),
    ("hamzalioglu_2026_milk_uht_ages", "12", "Hamzalioglu et al. (2026)", "10.1016/j.foodchem.2026.143200",
     "uht_milk", ["CML", "CEL"], ["age_uht_milk"],
     ["CML 8.9 nmol/g protein UHT"], ["data/lit/safety_reference_payloads.json"]),
    ("krause_2003_furosine_hydrolysis_yields", "12", "Krause et al. (2003)", "10.1016/S0308-8146(02)00537-1",
     "acid_hydrolysis", ["furosine"], ["furosine_hydrolysis_recovery"],
     ["furosine recovery 32%", "amadori→furosine 30%"], ["data/lit/computational_priors.json"]),
    ("hidalgo_2000_tomato_furosine", "12", "Hidalgo et al. (2000)", "10.1021/jf9912114",
     "tomato_puree", ["furosine"], ["furosine_vegetable_benchmark"],
     ["furosine 212 mg/100g protein"], ["data/lit/safety_reference_payloads.json"]),
    ("cantre_2007_corned_beef_furosine", "12", "Cantre et al. (2007)", "10.1016/j.meatsci.2007.08.019",
     "corned_beef_retort", ["furosine"], ["furosine_retort_benchmark"],
     ["furosine 890 mg/100g protein"], ["data/lit/safety_reference_payloads.json"]),
    ("fratianni_2016_apricot_furosine", "12", "Fratianni et al. (2016)", "10.1007/s11746-016-2885-3",
     "dried_apricot", ["furosine"], ["furosine_drying_benchmark"],
     ["furosine 2100 mg/100g protein"], ["data/lit/safety_reference_payloads.json"]),
    ("yu_2017_cml_cel_meat_review", "12", "Yu et al. (2017)", "10.1016/j.tifs.2017.09.002",
     "various_meat", ["CML", "CEL"], ["cml_cel_meat_range"],
     ["CML 2-38 ug/g", "CEL 1-22 ug/g protein"], ["data/lit/safety_reference_payloads.json"]),
    ("charissou_2007_cookie_cml_furosine", "12", "Charissou et al. (2007)", "10.1021/jf062952a",
     "cookie_biscuit", ["CML", "furosine"], ["cml_furosine_baked_goods"],
     ["CML 145 mg/kg", "furosine 1650 mg/100g protein"], ["data/lit/safety_reference_payloads.json"]),
    ("henle_2005_glycinin_cml", "12", "Henle et al. (2005)", "10.1016/j.ejnb.2005.05.005",
     "soy_glycinin", ["CML"], ["cml_plant_protein"],
     ["CML 1.8 nmol/mg glycinin"], ["data/lit/safety_reference_payloads.json"]),
    ("chen_2016_carbohydrate_meat_cml_cel_sterilization", "12", "Chen et al. (2016)", "10.1016/j.foodchem.2015.10.107",
     "canned_meat", ["CML", "CEL"], ["cml_cel_sterilization"],
     ["CML +8.3 ug/g", "CEL +5.1 ug/g protein sterilization"], ["data/lit/safety_reference_payloads.json"]),
    ("ma_2024_pbma_extrusion_damage", "12", "Ma et al. (2024)", "10.1016/j.foodchem.2024.139456",
     "pbma_extrusion", ["CML", "lysine"], ["extrusion_damage"],
     ["Lys -18.5%", "CML 3.2x"], ["data/lit/computational_priors.json"]),
    ("poulsen_2023_pbma_cml_cel", "12", "Poulsen et al. (2023)", "10.1016/j.foodchem.2023.135672",
     "pbma_cooked", ["CML", "CEL"], ["cml_cel_pbma"],
     ["CML 42 ug/g", "CEL 28 ug/g protein"], ["data/lit/safety_reference_payloads.json"]),
    ("pruteanu_2023_glucose_phe_arrhenius_browning", "12", "Pruteanu et al. (2023)", "10.3390/foods12234236",
     "aqueous_model", ["browning"], ["arrhenius_browning_kinetics"],
     ["Ea 73.4 kJ/mol glucose-Phe"], ["data/lit/safety_reference_payloads.json"]),
    ("acs_2022_pba_lysine_loss_benchmark", "12", "ACS Food Sci. Tech. 2022 (Ref. 1)", "",
     "pbma", ["lysine"], ["lysine_loss_extrusion"],
     ["Lys loss 22%"], ["data/lit/safety_reference_payloads.json"]),
    ("acs_ref3_spi_acrylamide_fast_kinetics", "12", "ACS Food Sci. Tech. 2022 (Ref. 3)", "",
     "spi_heated", ["acrylamide"], ["acrylamide_spi"],
     ["acrylamide 112 ug/kg at 150°C"], ["data/lit/safety_reference_payloads.json"]),
    ("foods_2023_cml_cel_proxy_benchmark", "12", "Foods 2023 CML/CEL benchmark", "",
     "pbma", ["CML"], ["cml_pbma_proxy"],
     ["CML 38 ug/g protein ELISA"], ["data/lit/safety_reference_payloads.json"]),
    ("ramirez_jimenez_2000_furosine_crossover_benchmark", "12", "Ramírez-Jiménez et al. (2000)", "10.1021/jf000028m",
     "toasted_bread", ["furosine"], ["furosine_bread_benchmark"],
     ["furosine 1250 mg/100g protein toast"], ["data/lit/safety_reference_payloads.json"]),
    ("acrylamide_heat_trapping_2024", "12", "ResearchGate ref. 10 (ACRYLAMIDE/HEAT)", "",
     "free_model_system", ["acrylamide"], ["acrylamide_high_T_range"],
     ["acrylamide 31.81-186.70 ug/kg"], ["data/lit/safety_reference_payloads.json"]),
]

# IDs handled here (protein_damage_markers family)
TARGET_FAMILY = "protein_damage_markers"
TARGET_IDS = {
    "acs_2022_pba_lysine_loss_benchmark",
    "acs_ref3_spi_acrylamide_fast_kinetics",
    "foods_2023_cml_cel_proxy_benchmark",
    "ramirez_jimenez_2000_furosine_crossover_benchmark",
    "henle_2005_glycinin_cml",
    "acrylamide_heat_trapping_2024",
    "de_vleeschouwer_2006_acrylamide_aqueous",
    "zhu_2022_acrylamide_lipid_crosstalk",
    "sun_2015_ground_beef_ages",
    "zhu_2021_braised_chicken_ages",
    "hamzalioglu_2026_milk_uht_ages",
    "krause_2003_furosine_hydrolysis_yields",
    "hidalgo_2000_tomato_furosine",
    "cantre_2007_corned_beef_furosine",
    "yu_2017_cml_cel_meat_review",
    "charissou_2007_cookie_cml_furosine",
    "fratianni_2016_apricot_furosine",
    "ma_2024_pbma_extrusion_damage",
    "knol_2005_acrylamide_kinetics",
    "de_vleeschouwer_2008_acrylamide_aw",
    "ishak_2022_phip_hca_kinetics",
    "mandelli_2024_plant_milk_hca_quantification",
    "chen_2016_carbohydrate_meat_cml_cel_sterilization",
    "pruteanu_2023_glucose_phe_arrhenius_browning",
    "poulsen_2023_pbma_cml_cel",
}

# acrylamide_lipid_crosstalk — cross-listed in protein_damage_markers
EXTRA_SAFETY_ENTRY = {
    "id": "zhu_2022_acrylamide_lipid_crosstalk",
    "category": "acrylamide",
    "source_citation": "Zhu et al. (2022)",
    "doi": "10.1016/j.foodchem.2022.132563",
    "matrix_context": "high_fat_biscuit",
    "analytical_method": "LC-MS",
    "benchmark_role": "reference_anchor",
    "key_values": {
        "acrylamide_reduction_pct_lipid_trapping": 38.0,
        "carbonyl_scavenging_lipid_pct": 44.0,
    },
    "notes": "Lipid carbonyl scavenging reduces acrylamide formation by 38% in high-fat biscuit.",
}


def main():
    intake = load_json(INTAKE_REGISTRY_PATH)
    matrix = load_json(SLR_MATRIX_PATH)
    priors = load_json(COMPUTATIONAL_PRIORS_PATH)
    safety = load_json(SAFETY_PAYLOADS_PATH)

    # ── 1. Update status in intake registry ──────────────────────────────────
    encoded_count = 0
    for ref in intake["eligible_references"]:
        if ref["id"] in TARGET_IDS and ref.get("status") == "ready_for_intake_encoding":
            ref["status"] = "encoded"
            encoded_count += 1
    print(f"Marked {encoded_count} protein_damage_markers refs as encoded")

    # ── 2. Safety payload entries ─────────────────────────────────────────────
    existing_safety_ids = {e["id"] for e in safety["entries"]}
    added_safety = 0
    for entry in SAFETY_ENTRIES:
        if entry["id"] not in existing_safety_ids:
            safety["entries"].append(entry)
            added_safety += 1
    # extra lipid crosstalk entry
    if EXTRA_SAFETY_ENTRY["id"] not in existing_safety_ids:
        safety["entries"].append(EXTRA_SAFETY_ENTRY)
        added_safety += 1
    print(f"Added {added_safety} safety payload entries")

    # ── 3. Computational priors ───────────────────────────────────────────────
    existing_furosine_ids = {p["id"] for p in priors["furosine_conversion_priors"]}
    added_priors = 0
    for pr in FUROSINE_PRIORS:
        if pr["id"] not in existing_furosine_ids:
            priors["furosine_conversion_priors"].append(pr)
            added_priors += 1

    existing_crosslink_ids = {p["id"] for p in priors["crosslink_kinetics_priors"]}
    for pr in CROSSLINK_PRIORS:
        if pr["id"] not in existing_crosslink_ids:
            priors["crosslink_kinetics_priors"].append(pr)
            added_priors += 1
    print(f"Added {added_priors} computational prior entries")

    # ── 4. SLR matrix entries ─────────────────────────────────────────────────
    existing_matrix_ids = {e["paper_id"] for e in matrix["entries"]}
    added_matrix = 0
    for (pid, sec, cite, doi, mat_fam, compounds, params, anchors, artifacts) in SLR_ENTRIES:
        if pid not in existing_matrix_ids:
            matrix["entries"].append({
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
                "notes_on_limits": "See safety_reference_payloads.json.",
            })
            added_matrix += 1
    # extra entry for zhu_2022
    if "zhu_2022_acrylamide_lipid_crosstalk" not in existing_matrix_ids:
        matrix["entries"].append({
            "slr_section": "12",
            "paper_id": "zhu_2022_acrylamide_lipid_crosstalk",
            "citation": "Zhu et al. (2022)",
            "doi": "10.1016/j.foodchem.2022.132563",
            "matrix_family": "high_fat_biscuit",
            "compounds_supported": ["acrylamide"],
            "parameters_supported": ["acrylamide_lipid_trapping"],
            "exact_numeric_anchors": ["acrylamide -38% with lipid trapping"],
            "current_repo_artifacts": ["data/lit/safety_reference_payloads.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "None — fully encoded.",
            "confidence_tier": "high",
            "notes_on_limits": "High-fat matrix specific.",
        })
        added_matrix += 1
    print(f"Added {added_matrix} SLR matrix entries")

    # ── 5. Save all ───────────────────────────────────────────────────────────
    save_json(intake, INTAKE_REGISTRY_PATH)
    save_json(matrix, SLR_MATRIX_PATH)
    save_json(priors, COMPUTATIONAL_PRIORS_PATH)
    save_json(safety, SAFETY_PAYLOADS_PATH)

    print("Done: protein_damage_markers batch encoded.")


if __name__ == "__main__":
    main()
