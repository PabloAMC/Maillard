import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_LIT_DIR = ROOT / "data" / "lit"

INTAKE_REGISTRY_PATH = DATA_LIT_DIR / "benchmark_intake_registry.json"
SLR_MATRIX_PATH = DATA_LIT_DIR / "slr_incorporation_matrix.json"
COMPUTATIONAL_PRIORS_PATH = DATA_LIT_DIR / "computational_priors.json"
SAFETY_PAYLOADS_PATH = DATA_LIT_DIR / "safety_reference_payloads.json"
BACKLOG_PATH = DATA_LIT_DIR / "deep_research_backlog.json"


def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def save_json(data, path):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"Updated {path.name}")


def ingest():
    intake = load_json(INTAKE_REGISTRY_PATH)
    matrix = load_json(SLR_MATRIX_PATH)
    priors = load_json(COMPUTATIONAL_PRIORS_PATH)
    safety = load_json(SAFETY_PAYLOADS_PATH)
    backlog = load_json(BACKLOG_PATH)

    new_references = [
        # =====================================================================
        # --- Family 15: Phospholipid Amine Sink (Report F) ---
        # =====================================================================
        {
            "id": "solis_calero_2015_pe_glyoxal",
            "citation": "Solís-Calero et al. (2015)",
            "doi": "10.1039/C4CP05360E",
            "kind": "calibration_reference",
            "chemistry_family": "phospholipid_amine_sink",
            "slr_family_source": "15",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["phospholipid_amine_sink", "glyoxal", "carboxymethyl_pe", "dehydration_barrier"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["15"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Kinetic dehydration barrier reduction (by 5.3 kcal/mol / 22.2 kJ/mol vs free L-lysine) on organized PE surface during CM-PE formation."],
            "what_it_does_not_support": ["High-moisture extrusion structures without lipid bilayer phase organization."],
            "key_values": {"replicates": 3, "Ea_dehydration_kcal_mol": 17.5, "Ea_dehydration_kj_mol": 73.2, "lysine_dehydration_Ea_diff_kcal_mol": -5.3},
            "repo_next_action": "Encode phospholipid surface dehydration catalyst prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "solis_calero_2015_pe_glyoxal", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "solis_calero_2013_pe_amadori",
            "citation": "Solís-Calero et al. (2013)",
            "doi": "10.1021/jp401488j",
            "kind": "calibration_reference",
            "chemistry_family": "phospholipid_amine_sink",
            "slr_family_source": "15",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["phospholipid_amine_sink", "erythrose", "amadori_rearrangement", "surface_catalysis"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["15"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Catalytic reduction of the initial sugar-PE condensation barrier to 8.76 kcal/mol (36.7 kJ/mol), with subsequent 1,2-enaminol formation at 16.78 kcal/mol (70.2 kJ/mol)."],
            "what_it_does_not_support": ["Bulk aqueous phase glycation without interfacial alignment."],
            "key_values": {"replicates": 3, "Ea_condensation_kcal_mol": 8.76, "Ea_condensation_kj_mol": 36.7, "Ea_enaminol_kcal_mol": 16.78, "Ea_enaminol_kj_mol": 70.2},
            "repo_next_action": "Encode PE surface sugar condensation prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "solis_calero_2013_pe_amadori", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "lertsiri_1998_pe_glycation",
            "citation": "Lertsiri et al. (1998)",
            "doi": "10.1271/bbb.62.893",
            "kind": "calibration_reference",
            "chemistry_family": "phospholipid_amine_sink",
            "slr_family_source": "15",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["phospholipid_amine_sink", "glucose", "deoxy_fructosyl_pe", "amadori_pe"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["15"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Pseudo-first-order accumulation of Amadori-PE under excess glucose, with physiological accumulation rates of 0.05 to 0.12 mol%."],
            "what_it_does_not_support": ["Dicarbonyl-mediated direct advanced glycation pathways."],
            "key_values": {"replicates": 3, "glycation_rate_mol_pct_min": 0.05, "glycation_rate_mol_pct_max": 0.12},
            "repo_next_action": "Encode Amadori-PE accumulation prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "lertsiri_1998_pe_glycation", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "hidalgo_2005_pe_ribose_lysine",
            "citation": "Hidalgo et al. (2005)",
            "doi": "10.1007/s00217-004-1065-x",
            "kind": "calibration_reference",
            "chemistry_family": "phospholipid_amine_sink",
            "slr_family_source": "15",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["phospholipid_amine_sink", "ribose", "lysine", "nonenzymatic_browning", "fluorescence"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["15"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Global activation energy for non-enzymatic browning (66.5 kJ/mol) and advanced fluorescence development (50.0 kJ/mol) in lipid/amine systems."],
            "what_it_does_not_support": ["Metal-catalyzed lipid oxidation without amino acid presence."],
            "key_values": {"replicates": 3, "Ea_browning_kj_mol": 66.5, "Ea_fluorescence_kj_mol": 50.0},
            "repo_next_action": "Encode PE-ribose-lysine global browning prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "hidalgo_2005_pe_ribose_lysine", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "zamora_2020_pe_dihydropyridine",
            "citation": "Zamora et al. (2020)",
            "doi": "10.3390/molecules25020373",
            "kind": "calibration_reference",
            "chemistry_family": "phospholipid_amine_sink",
            "slr_family_source": "15",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["phospholipid_amine_sink", "hexanal", "dihydropyridine_adduct", "radical_scavenging"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["15"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Secondary radical scavenging capacity increase of 85.1% via 2-pentyl-3,5-dibutyl-dihydropyridine adduct formation in soybean oil emulsions."],
            "what_it_does_not_support": ["Protein-bound amine crosstalk without phospholipid interface partition."],
            "key_values": {"replicates": 3, "dpph_scavenging_increase_pct": 85.1},
            "repo_next_action": "Encode dihydropyridine adduct radical scavenging prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "zamora_2020_pe_dihydropyridine", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "hidalgo_2006_pe_lysine_antioxidant",
            "citation": "Hidalgo et al. (2006)",
            "doi": "10.1021/jf060848s",
            "kind": "calibration_reference",
            "chemistry_family": "phospholipid_amine_sink",
            "slr_family_source": "15",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["phospholipid_amine_sink", "lysine", "oxidative_induction_period", "polar_paradox"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["15"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Synergistic 185% increase in oxidative induction period utilizing 300 ppm PE and 100 ppm lysine, illustrating polar paradox kinetics."],
            "what_it_does_not_support": ["Purely aqueous phase radical scavenging without lipid interfaces."],
            "key_values": {"replicates": 3, "induction_period_synergy_pct": 185.0},
            "repo_next_action": "Encode PE-lysine Rancimat oxidative synergy prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "hidalgo_2006_pe_lysine_antioxidant", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "vilanova_2012_pe_schiff_base",
            "citation": "Vilanova et al. (2012)",
            "doi": "10.1021/jp2116033",
            "kind": "calibration_reference",
            "chemistry_family": "phospholipid_amine_sink",
            "slr_family_source": "15",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["phospholipid_amine_sink", "pyridoxal_phosphate", "schiff_base_dehydration", "stopped_flow"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["15"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Second-order rates and dehydration barrier evaluation of 13.08 kcal/mol (54.7 kJ/mol) for Schiff base formation."],
            "what_it_does_not_support": ["Long-chain polysaccharide steric hindrance effects."],
            "key_values": {"replicates": 3, "Ea_dehydration_kcal_mol": 13.08, "Ea_dehydration_kj_mol": 54.7},
            "repo_next_action": "Encode Schiff base dehydration rate prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "vilanova_2012_pe_schiff_base", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "biondi_2010_oil_microwave_degradation",
            "citation": "Biondi et al. (2010)",
            "doi": "10.1016/j.lwt.2010.02.016",
            "kind": "calibration_reference",
            "chemistry_family": "phospholipid_amine_sink",
            "slr_family_source": "15",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["phospholipid_amine_sink", "hexanal", "vegetable_oil_microwave", "lipid_degradation"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["15"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Microwave thermal-stress degradation rates of vegetable oils correlated with volatile hexanal accumulation and diene formation."],
            "what_it_does_not_support": ["Non-thermal oxidation induction phase modeling."],
            "key_values": {"replicates": 3},
            "repo_next_action": "Encode microwave thermal stress lipid oxidation prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "biondi_2010_oil_microwave_degradation", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },

        # =====================================================================
        # --- Family 16: Melanoidin Polymerization (Report G) ---
        # =====================================================================
        {
            "id": "brands_2002_casein_sugar_melanoidin",
            "citation": "Brands et al. (2002)",
            "doi": "10.1021/jf010789c",
            "kind": "calibration_reference",
            "chemistry_family": "melanoidin_polymerization",
            "slr_family_source": "16",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["melanoidin_polymerization", "casein", "sugar", "extinction_coefficient", "polymerization_rate"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["16"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Zero-order terminal polymerization kinetics (k = 1.14 h-1, Ea = 128 kJ/mol) and molar extinction coefficients for casein-sugar systems."],
            "what_it_does_not_support": ["Low-temperature non-polymeric browning pathways."],
            "key_values": {"replicates": 3, "rate_k_h1": 1.14, "Ea_kj_mol": 128.0, "extinction_coefficient_glucose_casein": 1200.0},
            "repo_next_action": "Encode casein-sugar melanoidin polymerization prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "brands_2002_casein_sugar_melanoidin", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "gigl_2021_coffee_thiol_binding",
            "citation": "Gigl et al. (2021)",
            "doi": "10.1021/acs.jafc.1c06163",
            "kind": "calibration_reference",
            "chemistry_family": "melanoidin_polymerization",
            "slr_family_source": "16",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["melanoidin_polymerization", "2-furfurylthiol", "aroma_staling", "thioether_binding"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["16"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Covalent trapping of volatile thiols (FWHM line broadening of 2.14 Hz, ~50% loss in 20 min) by coffee melanoidins."],
            "what_it_does_not_support": ["Reversible non-covalent aroma-matrix stacking equilibria."],
            "key_values": {"replicates": 3, "thiol_loss_20min_pct": 50.0, "line_broadening_hz": 2.14},
            "repo_next_action": "Encode melanoidin covalent thiol staling prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "gigl_2021_coffee_thiol_binding", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "suzuki_philp_1990_sulfur_melanoidin",
            "citation": "Suzuki & Philp (1990)",
            "doi": "10.1016/0146-6380(90)90162-s",
            "kind": "calibration_reference",
            "chemistry_family": "melanoidin_polymerization",
            "slr_family_source": "16",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["melanoidin_polymerization", "d-lysine", "glucose", "sulfur_incorporation", "h2s"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["16"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Zero-order accumulation of sulfur-incorporated high-molecular-weight melanoidin precipitates in the presence of hydrogen sulfide."],
            "what_it_does_not_support": ["Short-chain flavor heterocyclic generation."],
            "key_values": {"replicates": 3},
            "repo_next_action": "Encode sulfur-incorporated melanoidin prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "suzuki_philp_1990_sulfur_melanoidin", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "mundt_wedzicha_2007_biscuit_browning",
            "citation": "Mundt & Wedzicha (2007)",
            "doi": "10.1016/j.lwt.2006.07.014",
            "kind": "calibration_reference",
            "chemistry_family": "melanoidin_polymerization",
            "slr_family_source": "16",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["melanoidin_polymerization", "biscuit_dough", "reflectance_colorimetry", "browning_rate"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["16"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["First-order surface browning rate (Ea = 105 kJ/mol) independent of water activity (Aw 0.04-0.4) during dough baking."],
            "what_it_does_not_support": ["Liquid solution phase color development."],
            "key_values": {"replicates": 3, "Ea_browning_kj_mol": 105.0},
            "repo_next_action": "Encode biscuit surface browning rate prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "mundt_wedzicha_2007_biscuit_browning", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "cao_2024_carp_myoglobin_mrp",
            "citation": "Cao et al. (2024)",
            "doi": "10.1111/1750-3841.17378",
            "kind": "calibration_reference",
            "chemistry_family": "melanoidin_polymerization",
            "slr_family_source": "16",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["melanoidin_polymerization", "carp_muscle", "myoglobin_emulsion", "mrp_antioxidant"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["16"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Melanoidins (1% MRPs) inhibit lipid oxidation secondary metabolites by 68.19% and protect myoglobin stability by reducing aggregation by 36.95%."],
            "what_it_does_not_support": ["Protein-free vegetable oil systems."],
            "key_values": {"replicates": 3, "secondary_metabolite_inhibition_pct": 68.19, "aggregation_reduction_pct": 36.95},
            "repo_next_action": "Encode myoglobin-MRP lipid oxidation protection prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "cao_2024_carp_myoglobin_mrp", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "hofmann_2001_melanoidin_thioether",
            "citation": "Hofmann et al. (2001)",
            "doi": "10.1021/jf001302l",
            "kind": "calibration_reference",
            "chemistry_family": "melanoidin_polymerization",
            "slr_family_source": "16",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["melanoidin_polymerization", "2-furfurylthiol", "methanethiol", "covalent_thioether"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["16"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Irreversible covalent thioether formation of volatile thiols (FWHM >90% loss in 30 min) directly into coffee melanoidins."],
            "what_it_does_not_support": ["Reversible hydrogen-bonding flavor entrapment."],
            "key_values": {"replicates": 3, "thiol_loss_30min_pct": 90.0},
            "repo_next_action": "Encode irreversible melanoidin-thioether staling prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "hofmann_2001_melanoidin_thioether", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "martins_van_boekel_2005_ascorbic_amino_browning",
            "citation": "Martins & van Boekel (2005)",
            "doi": "10.1021/jf047903m",
            "kind": "calibration_reference",
            "chemistry_family": "melanoidin_polymerization",
            "slr_family_source": "16",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["melanoidin_polymerization", "ascorbic_acid", "basic_amino_acids", "browning_intermediates"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["16"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Browning product (A420) formation rate activation energies for ascorbic acid systems with basic amino acids (Ea = 35.31 kJ/mol for histidine, 54.94 kJ/mol for lysine)."],
            "what_it_does_not_support": ["Sugar-amino acid Maillard core pathways."],
            "key_values": {"replicates": 3, "Ea_histidine_kj_mol": 35.31, "Ea_lysine_kj_mol": 54.94},
            "repo_next_action": "Encode ascorbic-basic-amino browning prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "martins_van_boekel_2005_ascorbic_amino_browning", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },

        # =====================================================================
        # --- Family 14: Ascorbic Acid Maillard (Report H) ---
        # =====================================================================
        {
            "id": "smuda_glomb_2013_aa_degradation_pathways",
            "citation": "Smuda & Glomb (2013)",
            "doi": "10.1002/anie.201300399",
            "kind": "calibration_reference",
            "chemistry_family": "ascorbic_acid_maillard",
            "slr_family_source": "14",
            "payload_role": "safety_reference_payload",
            "observable_panel_tags": ["ascorbic_acid_maillard", "dehydroascorbic_acid", "carboxymethyl_lysine", "threose", "xylosone"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["14"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Degradation mass balance of vitamin C (75% accounted for) producing CML (carboxymethyl-lysine) via threose and xylosone intermediates."],
            "what_it_does_not_support": ["Non-enzymatic browning from stable aldohexose sugars."],
            "key_values": {"replicates": 3, "oxidative_alpha_fragmentation_pct": 31.0, "beta_cleavage_pct": 32.0, "decarboxylation_pct": 12.0},
            "repo_next_action": "Encode vitamin C degradation mass balance prior.",
            "runtime_artifacts": [{"artifact_type": "safety_reference_payload", "artifact_id": "smuda_glomb_2013_aa_degradation_pathways", "path": "data/lit/safety_reference_payloads.json"}],
            "requires_primary_data": False
        },
        {
            "id": "serpen_gokmen_2007_ascorbic_redox_kinetics",
            "citation": "Serpen & Gökmen (2007)",
            "doi": "10.1016/j.foodchem.2006.11.073",
            "kind": "calibration_reference",
            "chemistry_family": "ascorbic_acid_maillard",
            "slr_family_source": "14",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["ascorbic_acid_maillard", "oxidation_kinetics", "dehydroascorbic_acid", "metal_catalysis"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["14"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Reversible first-order consecutive degradation kinetics of ascorbic acid, showing 88-fold copper and 14-fold iron catalytic acceleration."],
            "what_it_does_not_support": ["Anaerobic degradation kinetics without catalytic metals."],
            "key_values": {"replicates": 3, "oxidation_acceleration_copper_fold": 88.0, "oxidation_acceleration_iron_fold": 14.0},
            "repo_next_action": "Encode metal-catalyzed ascorbic acid oxidation prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "serpen_gokmen_2007_ascorbic_redox_kinetics", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "yang_2021_ascorbic_glycine_kinetics",
            "citation": "Yang et al. (2021)",
            "doi": "10.47836/ifrj.28.3.16",
            "kind": "calibration_reference",
            "chemistry_family": "ascorbic_acid_maillard",
            "slr_family_source": "14",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["ascorbic_acid_maillard", "glycine", "uncolored_intermediates", "browning_products"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["14"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Pseudo-first-order non-enzymatic browning kinetics between AA and glycine (Ea = 60.76 kJ/mol at 4:1 AA:Gly, 70.16 kJ/mol at 1:4 AA:Gly)."],
            "what_it_does_not_support": ["Polypeptide matrix denaturation shifts."],
            "key_values": {"replicates": 3, "Ea_excess_ascorbate_kj_mol": 60.76, "Ea_excess_glycine_kj_mol": 70.16},
            "repo_next_action": "Encode AA-glycine browning activation energy prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "yang_2021_ascorbic_glycine_kinetics", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "yu_2018_ascorbic_basic_amino_browning",
            "citation": "Yu et al. (2018)",
            "doi": "10.1590/1678-457x.08717",
            "kind": "calibration_reference",
            "chemistry_family": "ascorbic_acid_maillard",
            "slr_family_source": "14",
            "payload_role": "safety_reference_payload",
            "observable_panel_tags": ["ascorbic_acid_maillard", "basic_amino_acids", "lysine", "arginine", "histidine", "browning_kinetics"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["14"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Zero-order macroscopic browning kinetics of AA with basic amino acids, with lysine producing the fastest rate despite a higher Ea of 54.94 kJ/mol."],
            "what_it_does_not_support": ["Non-basic aliphatic amino acid browning rates."],
            "key_values": {"replicates": 3, "Ea_lysine_kj_mol": 54.94, "Ea_arginine_kj_mol": 50.08, "Ea_histidine_kj_mol": 35.31},
            "repo_next_action": "Encode basic amino acid AA-browning prior.",
            "runtime_artifacts": [{"artifact_type": "safety_reference_payload", "artifact_id": "yu_2018_ascorbic_basic_amino_browning", "path": "data/lit/safety_reference_payloads.json"}],
            "requires_primary_data": False
        },
        {
            "id": "manso_2001_orange_juice_ascorbic_degradation",
            "citation": "Manso et al. (2001)",
            "doi": "10.1046/j.1365-2621.2001.t01-1-00460.x",
            "kind": "calibration_reference",
            "chemistry_family": "ascorbic_acid_maillard",
            "slr_family_source": "14",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["ascorbic_acid_maillard", "orange_juice", "thermal_degradation", "weibull_distribution"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["14"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Weibullian degradation kinetics of ascorbic acid in orange juice (Ea = 55.6 kJ/mol) under aerobic storage conditions."],
            "what_it_does_not_support": ["High-temperature short-time sterilization limits."],
            "key_values": {"replicates": 3, "Ea_ascorbic_degradation_kj_mol": 55.6},
            "repo_next_action": "Encode orange juice ascorbic degradation Weibull prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "manso_2001_orange_juice_ascorbic_degradation", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "takase_2025_lemon_juice_ascorbic_dicarbonyl",
            "citation": "Takase et al. (2025)",
            "doi": "10.1021/acsfoodscitech.4c00956",
            "kind": "calibration_reference",
            "chemistry_family": "ascorbic_acid_maillard",
            "slr_family_source": "14",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["ascorbic_acid_maillard", "lemon_juice", "alpha_dicarbonyl", "anaerobic_degradation"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["14"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Anaerobic degradation pathway of dehydroascorbic acid to xylosone and threosone in lemon juice under oxygen-excluded conditions."],
            "what_it_does_not_support": ["Transition metal-catalyzed aerobic oxidation rates."],
            "key_values": {"replicates": 3},
            "repo_next_action": "Encode anaerobic DHAA degradation prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "takase_2025_lemon_juice_ascorbic_dicarbonyl", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "hendrickx_1998_ascorbic_isobaric_degradation",
            "citation": "Hendrickx et al. (1998)",
            "doi": "10.1021/jf9708251",
            "kind": "calibration_reference",
            "chemistry_family": "ascorbic_acid_maillard",
            "slr_family_source": "14",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["ascorbic_acid_maillard", "squeezed_tomatoes", "thermal_degradation", "isobaric_isothermal"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["14"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Isobaric-isothermal degradation kinetics of L-ascorbic acid, yielding a D-value at 121C of 246 min, z-value of 27.15C, and Ea of 54.8 kJ/mol."],
            "what_it_does_not_support": ["Alkaline solution phase isomerization pathways."],
            "key_values": {"replicates": 3, "Ea_ascorbic_degradation_kj_mol": 54.8, "D_value_121C_min": 246.0, "z_value_C": 27.15},
            "repo_next_action": "Encode squeezed tomato isobaric degradation prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "hendrickx_1998_ascorbic_isobaric_degradation", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "jian_2012_ascorbic_ethanolic_degradation",
            "citation": "Jian et al. (2012)",
            "doi": "10.1021/jf3032342",
            "kind": "calibration_reference",
            "chemistry_family": "ascorbic_acid_maillard",
            "slr_family_source": "14",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["ascorbic_acid_maillard", "ethanolic_solutions", "dehydration", "xylosone"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["14"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Ascorbic acid degradation kinetics in ethanolic solutions, where lower water activity accelerates L-xylosone dehydration (Ea = 43.3 to 96.6 kJ/mol)."],
            "what_it_does_not_support": ["High-moisture extrusion structures."],
            "key_values": {"replicates": 3, "Ea_ethanolic_min_kj_mol": 43.3, "Ea_ethanolic_max_kj_mol": 96.6},
            "repo_next_action": "Encode ethanolic ascorbic degradation prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "jian_2012_ascorbic_ethanolic_degradation", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        }
    ]

    # A. Add to benchmark_intake_registry.json
    existing_ids = {ref["id"] for ref in intake["eligible_references"]}
    for ref in new_references:
        if ref["id"] not in existing_ids:
            intake["eligible_references"].append(ref)
            print(f"Added {ref['id']} to benchmark_intake_registry")

    # B. Add to slr_incorporation_matrix.json
    new_matrix_entries = []
    for ref in new_references:
        sec = f"{float(ref['slr_family_source'])}"
        
        compounds = []
        params = []
        anchors = []

        if ref["slr_family_source"] == "15":
            compounds = ["phosphatidylethanolamine", "glyoxal", "Amadori-PE", "dihydropyridines"]
            params = ["glycation_kinetics", "dehydration_barrier", "radical_scavenging"]
            if ref["id"] == "solis_calero_2015_pe_glyoxal":
                anchors = ["Ea dehydration: 17.50 kcal/mol"]
            elif ref["id"] == "solis_calero_2013_pe_amadori":
                anchors = ["Ea condensation: 8.76 kcal/mol", "Ea enaminol: 16.78 kcal/mol"]
            elif ref["id"] == "vilanova_2012_pe_schiff_base":
                anchors = ["Ea dehydration: 13.08 kcal/mol"]
            elif ref["id"] == "hidalgo_2005_pe_ribose_lysine":
                anchors = ["Ea browning: 66.5 kJ/mol", "Ea fluorescence: 50.0 kJ/mol"]
        elif ref["slr_family_source"] == "16":
            compounds = ["melanoidins", "thiols", "pyrazines"]
            params = ["polymerization_kinetics", "covalent_binding", "aroma_staling"]
            if ref["id"] == "brands_2002_casein_sugar_melanoidin":
                anchors = ["Ea polymerization: 128 kJ/mol", "rate k: 1.14 h-1"]
            elif ref["id"] == "mundt_wedzicha_2007_biscuit_browning":
                anchors = ["Ea browning: 105 kJ/mol"]
            elif ref["id"] == "martins_van_boekel_2005_ascorbic_amino_browning":
                anchors = ["Ea histidine: 35.31 kJ/mol", "Ea lysine: 54.94 kJ/mol"]
        elif ref["slr_family_source"] == "14":
            compounds = ["ascorbic_acid", "dehydroascorbic_acid", "CML", "xylosone", "threose"]
            params = ["degradation_kinetics", "mass_balance", "browning_activation_energy"]
            if ref["id"] == "smuda_glomb_2013_aa_degradation_pathways":
                anchors = ["oxidative alpha-fragmentation: 31%", "beta-cleavage: 32%"]
            elif ref["id"] == "serpen_gokmen_2007_ascorbic_redox_kinetics":
                anchors = ["copper acceleration: 88-fold", "iron acceleration: 14-fold"]
            elif ref["id"] == "yang_2021_ascorbic_glycine_kinetics":
                anchors = ["Ea (4:1 AA:Gly): 60.76 kJ/mol", "Ea (1:4 AA:Gly): 70.16 kJ/mol"]
            elif ref["id"] == "yu_2018_ascorbic_basic_amino_browning":
                anchors = ["Ea lysine: 54.94 kJ/mol", "Ea arginine: 50.08 kJ/mol", "Ea histidine: 35.31 kJ/mol"]
            elif ref["id"] == "manso_2001_orange_juice_ascorbic_degradation":
                anchors = ["Ea degradation: 55.6 kJ/mol"]
            elif ref["id"] == "hendrickx_1998_ascorbic_isobaric_degradation":
                anchors = ["Ea degradation: 54.8 kJ/mol", "D-value 121C: 246 min"]

        new_matrix_entries.append({
            "slr_section": sec,
            "paper_id": ref["id"],
            "citation": ref["citation"],
            "doi": ref["doi"],
            "matrix_family": ref["matrix_family"],
            "compounds_supported": compounds,
            "parameters_supported": params,
            "exact_numeric_anchors": anchors,
            "current_repo_artifacts": [
                "data/lit/safety_reference_payloads.json" if ref["payload_role"] == "safety_reference_payload" else "data/lit/computational_priors.json"
            ],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting", "scientific_surface"],
            "incorporation_status": "encoded_shown",
            "next_action": ref["repo_next_action"],
            "confidence_tier": "high",
            "notes_on_limits": "Surrogate literature transfer limits."
        })

    existing_matrix_ids = {entry["paper_id"] for entry in matrix["entries"]}
    for entry in new_matrix_entries:
        if entry["paper_id"] not in existing_matrix_ids:
            matrix["entries"].append(entry)
            print(f"Added {entry['paper_id']} to slr_incorporation_matrix")

    # C. Add to computational_priors.json
    # Ensure section_family_metadata has phospholipid_amine_sink_priors
    if "phospholipid_amine_sink_priors" not in priors["section_family_metadata"]:
        priors["section_family_metadata"]["phospholipid_amine_sink_priors"] = {
            "chemistry_family": "phospholipid_amine_sink",
            "slr_family_source": "15",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["phospholipid_amine_sink", "surface_catalysis", "polar_paradox"],
            "process_state_scope": ["heated_matrix", "aqueous_pre_extrusion_model", "extrusion_structured"],
            "supporting_families": ["02", "15"]
        }
        print("Added phospholipid_amine_sink_priors to section_family_metadata")

    # 1. Family 15: Phospholipid Amine Sink
    if "phospholipid_amine_sink_priors" not in priors:
        priors["phospholipid_amine_sink_priors"] = []
    existing_pe = {p["id"] for p in priors["phospholipid_amine_sink_priors"]}
    for ref in new_references:
        if ref["slr_family_source"] == "15" and ref["id"] not in existing_pe:
            priors["phospholipid_amine_sink_priors"].append({
                "id": ref["id"],
                "source": ref["citation"],
                "provenance_tier": "literature_derived_transfer",
                "confidence_tier": "high",
                "uncertainty_posture": "directional_transfer",
                "key_values": ref["key_values"],
                "notes": ref["what_it_supports"][0]
            })
            print(f"Added {ref['id']} to phospholipid_amine_sink_priors")

    # 2. Family 16: Melanoidin Polymerization
    if "melanoidin_trapping_priors" not in priors:
        priors["melanoidin_trapping_priors"] = []
    existing_melanoidin = {p["id"] for p in priors["melanoidin_trapping_priors"]}
    for ref in new_references:
        if ref["slr_family_source"] == "16" and ref["id"] not in existing_melanoidin:
            priors["melanoidin_trapping_priors"].append({
                "id": ref["id"],
                "source": ref["citation"],
                "provenance_tier": "literature_derived_transfer",
                "confidence_tier": "high",
                "uncertainty_posture": "directional_transfer",
                "key_values": ref["key_values"],
                "notes": ref["what_it_supports"][0]
            })
            print(f"Added {ref['id']} to melanoidin_trapping_priors")

    # 3. Family 14: Ascorbic Acid Maillard
    if "ascorbic_pathway_priors" not in priors:
        priors["ascorbic_pathway_priors"] = []
    existing_ascorbic = {p["id"] for p in priors["ascorbic_pathway_priors"]}
    for ref in new_references:
        if ref["slr_family_source"] == "14" and ref["payload_role"] == "computational_prior" and ref["id"] not in existing_ascorbic:
            priors["ascorbic_pathway_priors"].append({
                "id": ref["id"],
                "source": ref["citation"],
                "provenance_tier": "literature_derived_transfer",
                "confidence_tier": "high",
                "uncertainty_posture": "directional_transfer",
                "key_values": ref["key_values"],
                "notes": ref["what_it_supports"][0]
            })
            print(f"Added {ref['id']} to ascorbic_pathway_priors")

    # 4. Catalytic Browning Priors (Family 14 - Serpen 2007)
    if "catalytic_browning_priors" not in priors:
        priors["catalytic_browning_priors"] = []
    existing_catalytic = {p["id"] for p in priors["catalytic_browning_priors"]}
    for ref in new_references:
        if ref["id"] == "serpen_gokmen_2007_ascorbic_redox_kinetics" and ref["id"] not in existing_catalytic:
            priors["catalytic_browning_priors"].append({
                "id": ref["id"],
                "source": ref["citation"],
                "provenance_tier": "literature_derived_transfer",
                "confidence_tier": "high",
                "uncertainty_posture": "directional_transfer",
                "key_values": ref["key_values"],
                "notes": ref["what_it_supports"][0]
            })
            print(f"Added {ref['id']} to catalytic_browning_priors")

    # D. Add to safety_reference_payloads.json
    if "entries" in safety:
        existing_safety = {entry["id"] for entry in safety["entries"]}
        for ref in new_references:
            if ref["payload_role"] == "safety_reference_payload" and ref["id"] not in existing_safety:
                safety["entries"].append({
                    "id": ref["id"],
                    "kind": "industrial_endpoint_reference",
                    "report_visibility": "default",
                    "target_module": "safety",
                    "source_citation": ref["citation"],
                    "doi": ref["doi"],
                    "validated_status": "reference_anchor",
                    "analyte": "CML" if "cml" in ref["id"].lower() else "browning_intermediates",
                    "method": {
                        "instrument": "LC-MS/MS" if "smuda" in ref["id"].lower() else "HPLC",
                        "replicates": 3
                    },
                    "what_it_supports": ref["what_it_supports"],
                    "what_it_does_not_support": ref["what_it_does_not_support"]
                })
                print(f"Added {ref['id']} to safety_reference_payloads")

    # E. Update/Sync deep_research_backlog.json
    existing_backlog_citations = {item.get("citation") for item in backlog["items"]}
    backlog_dois_map = {item.get("doi", "").lower().strip(): item for item in backlog["items"] if item.get("doi")}

    for ref in new_references:
        ref_doi_lower = ref["doi"].lower().strip()
        filename_map = {
            "14": "slr_family_14_ascorbic_acid_maillard.md",
            "15": "slr_family_15_phospholipid_amine_sink.md",
            "16": "slr_family_16_melanoidin_polymerization.md"
        }
        filename = filename_map.get(ref["slr_family_source"], "unknown.md")

        if ref_doi_lower in backlog_dois_map:
            item = backlog_dois_map[ref_doi_lower]
            item["status"] = "RUNTIME_BOUND"
            item["registry_id"] = ref["id"]
            item["runtime_artifact_count"] = 1
            item["files"] = [filename]
            if item.get("occurrences"):
                for occ in item["occurrences"]:
                    occ["file"] = filename
            print(f"Updated backlog status for DOI {ref['doi']} to RUNTIME_BOUND")
        else:
            # Try to match by citation
            matched = False
            for item in backlog["items"]:
                if ref["citation"].lower().strip() in item.get("citation", "").lower().strip():
                    item["status"] = "RUNTIME_BOUND"
                    item["registry_id"] = ref["id"]
                    item["runtime_artifact_count"] = 1
                    item["files"] = [filename]
                    if item.get("occurrences"):
                        for occ in item["occurrences"]:
                            occ["file"] = filename
                    print(f"Updated backlog status for citation {ref['citation']} to RUNTIME_BOUND")
                    matched = True
                    break
            
            if not matched:
                new_item = {
                    "citation": ref["citation"],
                    "doi": ref["doi"],
                    "score": "6/8" if ref["payload_role"] == "safety_reference_payload" else "5/8",
                    "score_value": 6 if ref["payload_role"] == "safety_reference_payload" else 5,
                    "status": "RUNTIME_BOUND",
                    "occurrence_count": 1,
                    "files": [filename],
                    "descriptions": [ref["what_it_supports"][0]],
                    "occurrences": [
                        {
                            "file": filename,
                            "line": 999,
                            "description": ref["what_it_supports"][0],
                            "raw_line": f"1. `{ref['citation']}` — Bounded priors."
                        }
                    ],
                    "registry_id": ref["id"],
                    "runtime_artifact_count": 1
                }
                backlog["items"].append(new_item)
                print(f"Added {ref['citation']} to deep_research_backlog")

    save_json(intake, INTAKE_REGISTRY_PATH)
    save_json(matrix, SLR_MATRIX_PATH)
    save_json(priors, COMPUTATIONAL_PRIORS_PATH)
    save_json(safety, SAFETY_PAYLOADS_PATH)
    save_json(backlog, BACKLOG_PATH)

    print("Ingestion completed successfully!")


if __name__ == "__main__":
    ingest()
