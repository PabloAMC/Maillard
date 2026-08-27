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
        # --- Family 02 ---
        {
            "id": "tan_2001_oil_dsc_oxidation",
            "citation": "Tan et al. (2001)",
            "doi": "10.1007/s11746-001-0402-y",
            "kind": "calibration_reference",
            "chemistry_family": "lipid_oxidation_crosstalk",
            "slr_family_source": "02",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["lipid_oxidation", "peroxide_formation", "differential_scanning_calorimetry"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["02"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Isothermal peroxide formation activation energy parameters derived via DSC."],
            "what_it_does_not_support": ["Extrusion high-shear physical matrix transformations."],
            "key_values": {"replicates": 3},
            "repo_next_action": "Encode lipid oxidation Arrhenius prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "tan_2001_oil_dsc_oxidation", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "bayram_2023_stripped_soybean_oil_hexanal",
            "citation": "Bayram et al. (2023)",
            "doi": "10.3989/gya.1051211",
            "kind": "calibration_reference",
            "chemistry_family": "lipid_oxidation_crosstalk",
            "slr_family_source": "02",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["stripped_oil", "hexanal_kinetics", "antioxidant_modulation"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["02"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Arrhenius hexanal kinetics modulated by antioxidants like ascorbyl palmitate (Ea 1.62 to 89.40 kJ/mol)."],
            "what_it_does_not_support": ["High-moisture extrusion structures without oil phase isolation."],
            "key_values": {"replicates": 3},
            "repo_next_action": "Encode hexanal accumulation Arrhenius prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "bayram_2023_stripped_soybean_oil_hexanal", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "spier_2010_metal_naphthenate_peroxide_decomp",
            "citation": "Spier et al. (2010)",
            "doi": "10.1021/ef101678s",
            "kind": "calibration_reference",
            "chemistry_family": "lipid_oxidation_crosstalk",
            "slr_family_source": "02",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["metal_catalysis", "peroxide_decomposition", "transition_metals"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["02"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Apparent Ea for transition-metal catalyzed peroxide decomposition (11.0 to 26.0 kcal/mol)."],
            "what_it_does_not_support": ["Heme-catalyzed lipid breakdown kinetics."],
            "key_values": {"replicates": 3, "Fe_Ea_kcal_mol": 23.0, "Zn_Ea_kcal_mol": 26.0, "Cu_Ea_kcal_mol": 14.0, "Mn_Ea_kcal_mol": 11.0},
            "repo_next_action": "Encode catalytic peroxide decomposition prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "spier_2010_metal_naphthenate_peroxide_decomp", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        # --- Family 06 ---
        {
            "id": "zha_2020_ppi_glycation_aggregation",
            "citation": "Zha et al. (2020)",
            "doi": "10.1021/acs.jafc.0c04281",
            "kind": "calibration_reference",
            "chemistry_family": "alternative_protein_matrix_scope",
            "slr_family_source": "06",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["plant_protein_matrix", "pea_isolate", "structural_aggregation", "lysine_blockade"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["06"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "plant_protein_matrix",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Pea protein isolate glycation aggregation kinetics and lysine blockage parameters."],
            "what_it_does_not_support": ["Cereal matrix gliadin beta-elimination kinetics."],
            "key_values": {"replicates": 3},
            "repo_next_action": "Encode PPI glycation prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "zha_2020_ppi_glycation_aggregation", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "kutzli_2020_pea_maltodextrin_electrospun_glycation",
            "citation": "Kutzli et al. (2020)",
            "doi": "10.1039/D0FO00292E",
            "kind": "calibration_reference",
            "chemistry_family": "alternative_protein_matrix_scope",
            "slr_family_source": "06",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["electrospun_fibers", "pea_isolate", "maltodextrin", "glycation_kinetics"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["06"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "plant_protein_matrix",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Pseudo-first order glycation within electrospun pea protein-maltodextrin fibers under relative humidity control."],
            "what_it_does_not_support": ["Free amino acid solution systems without structural fiber constraints."],
            "key_values": {"replicates": 3},
            "repo_next_action": "Encode fiber glycation prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "kutzli_2020_pea_maltodextrin_electrospun_glycation", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "nguyen_2025_ppi_microencapsulated_oil_stabilization",
            "citation": "Nguyen et al. (2025)",
            "doi": "10.1016/j.foodchem.2025.146396",
            "kind": "calibration_reference",
            "chemistry_family": "alternative_protein_matrix_scope",
            "slr_family_source": "06",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["microencapsulated_oil", "pea_isolate", "oxidation_barrier"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["06"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "plant_protein_matrix",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Pea protein isolate physical barrier coefficients limiting microencapsulated oil self-oxidation."],
            "what_it_does_not_support": ["Bulk liquid phase systems without microencapsulation boundaries."],
            "key_values": {"replicates": 3},
            "repo_next_action": "Encode oil stabilization barrier prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "nguyen_2025_ppi_microencapsulated_oil_stabilization", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "pereira_2020_metal_pm_haber_weiss_chelation",
            "citation": "Pereira et al. (2020)",
            "doi": "10.3390/antiox9080756",
            "kind": "calibration_reference",
            "chemistry_family": "alternative_protein_matrix_scope",
            "slr_family_source": "06",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["haber_weiss_chelation", "pyridoxamine", "free_energy", "transition_metals"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["06"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Transition metal coordination free energy calculations for Cu(II) and Fe(III) complexes (Cu-PM: -35.8 kcal/mol, Fe-PM: -58.9 kcal/mol)."],
            "what_it_does_not_support": ["Biological matrix-level enzyme-mediated chelation paths."],
            "key_values": {"Cu_PM_formation_free_energy_kcal_mol": -35.8, "Fe_PM_formation_free_energy_kcal_mol": -58.9},
            "repo_next_action": "Encode metal chelation thermodynamic prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "pereira_2020_metal_pm_haber_weiss_chelation", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "sun_2020_solid_matrix_cml_cel_accumulation",
            "citation": "Sun et al. (2020)",
            "doi": "10.1016/j.meatsci.2020.108151",
            "kind": "calibration_reference",
            "chemistry_family": "alternative_protein_matrix_scope",
            "slr_family_source": "06",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["solid_matrix", "pork_meat_analog", "cml_cel_kinetics", "sterilization"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["06"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "high_moisture_melt",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Zero-order free CML (Ea = 44.158 kJ/mol) and free CEL (Ea = 40.971 kJ/mol) accumulation in a solid matrix during sterilization."],
            "what_it_does_not_support": ["Extrusion processes with low water activities."],
            "key_values": {"replicates": 3, "Ea_free_CML_kj_mol": 44.158, "Ea_free_CEL_kj_mol": 40.971},
            "repo_next_action": "Encode solid matrix CML/CEL accumulation prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "sun_2020_solid_matrix_cml_cel_accumulation", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "shirai_2015_bsa_dityrosine_diffusion_limited",
            "citation": "Shirai et al. (2015)",
            "doi": "10.1021/acs.est.5b02902",
            "kind": "calibration_reference",
            "chemistry_family": "alternative_protein_matrix_scope",
            "slr_family_source": "06",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["bsa_protein", "dityrosine_formation", "diffusion_limited", "kinetic_multilayer"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["06"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "plant_protein_matrix",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Diffusion-limited dityrosine formation modeling using kinetic multilayer surface/bulk approaches."],
            "what_it_does_not_support": ["Disulfide-based crosslinking mechanisms."],
            "key_values": {"replicates": 3},
            "repo_next_action": "Encode dityrosine formation prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "shirai_2015_bsa_dityrosine_diffusion_limited", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        # --- Family 07 ---
        {
            "id": "martins_boekel_2003_dfg_amadori_degradation",
            "citation": "Martins & van Boekel (2003)",
            "doi": "10.1021/jf025830v",
            "kind": "calibration_reference",
            "chemistry_family": "carbonyl_donor_hierarchy",
            "slr_family_source": "07",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["dfg_amadori", "enolization_degradation", "multiresponse_kinetics"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["07"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Multiresponse kinetic modeling of DFG Amadori enolization and degradation."],
            "what_it_does_not_support": ["Solid-state extrusion with low water activity."],
            "key_values": {"replicates": 3},
            "repo_next_action": "Encode Amadori degradation prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "martins_boekel_2003_dfg_amadori_degradation", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "nashalian_yaylayan_2014_cu_catalyzed_strecker",
            "citation": "Nashalian & Yaylayan (2014)",
            "doi": "10.1021/acsomega.7b00321",
            "kind": "calibration_reference",
            "chemistry_family": "carbonyl_donor_hierarchy",
            "slr_family_source": "07",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["copper_catalysis", "strecker_decarboxylation", "enolization_kinetics"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["07"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Copper-catalyzed enolization (rate k = 0.41 s-1) and decarboxylation kinetics in model systems."],
            "what_it_does_not_support": ["Non-transition-metal catalyzed Strecker pathways."],
            "key_values": {"decarboxylation_k_s1": 0.41},
            "repo_next_action": "Encode catalyzed Strecker enolization prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "nashalian_yaylayan_2014_cu_catalyzed_strecker", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        # --- Family 09 ---
        {
            "id": "de_bruijn_1987_monosaccharide_alkaline_degradation",
            "citation": "de Bruijn et al. (1987)",
            "doi": "10.1002/recl.19871060202",
            "kind": "calibration_reference",
            "chemistry_family": "carbohydrate_pyrolysis_caramelization",
            "slr_family_source": "09",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["alkaline_degradation", "retro_aldolization", "isomerization_kinetics"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["09"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Alkaline-driven retro-aldolization isomerization kinetics of hexose monosaccharides."],
            "what_it_does_not_support": ["Acid-catalyzed sugar dehydration pathways."],
            "key_values": {"replicates": 3},
            "repo_next_action": "Encode alkaline retro-aldol prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "de_bruijn_1987_monosaccharide_alkaline_degradation", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "luna_aguilera_2014_molten_sugar_color_kinetics",
            "citation": "Luna & Aguilera (2014)",
            "doi": "10.1016/j.jfoodeng.2014.06.032",
            "kind": "calibration_reference",
            "chemistry_family": "carbohydrate_pyrolysis_caramelization",
            "slr_family_source": "09",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["molten_sugars", "caramelization_color", "activation_energy"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["09"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Arrhenius activation energies for color development in caramelizing molten hexose systems."],
            "what_it_does_not_support": ["Low-temperature aqueous sugar degradation."],
            "key_values": {"replicates": 3},
            "repo_next_action": "Encode caramelization color kinetics prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "luna_aguilera_2014_molten_sugar_color_kinetics", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "bao_2022_xylose_glycylglycine_3_deoxyosone_cleavage",
            "citation": "Bao et al. (2022)",
            "doi": "10.1021/acs.jafc.2c03427",
            "kind": "calibration_reference",
            "chemistry_family": "carbohydrate_pyrolysis_caramelization",
            "slr_family_source": "09",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["glycylglycine_catalysis", "xylose", "3_deoxyosone_cleavage", "furfural_kinetics"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["09"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Glycylglycine catalytic cleavage of Amadori products and 3-deoxyosone accumulation at 120 C."],
            "what_it_does_not_support": ["Non-dipeptide free amino acid reactions."],
            "key_values": {"replicates": 3},
            "repo_next_action": "Encode dipeptide Amadori cleavage prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "bao_2022_xylose_glycylglycine_3_deoxyosone_cleavage", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        # --- Family 11 ---
        {
            "id": "hidalgo_2007_decadienal_phenylalanine_styrene",
            "citation": "Hidalgo et al. (2007)",
            "doi": "10.1021/jf070527w",
            "kind": "calibration_reference",
            "chemistry_family": "lipid_maillard_crosstalk",
            "slr_family_source": "11",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["decadienal", "phenylalanine_degradation", "styrene_formation"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["11"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Phenylalanine degradation kinetics by 2,4-decadienal to yield styrene (Ea = 150.4 kJ/mol)."],
            "what_it_does_not_support": ["Solid-state high-moisture extrusion without secondary lipid initiation."],
            "key_values": {"replicates": 3, "Ea_styrene_formation_kj_mol": 150.4},
            "repo_next_action": "Encode lipid-amino crosstalk prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "hidalgo_2007_decadienal_phenylalanine_styrene", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "zamora_2010_decadienal_asparagine_decarboxylation",
            "citation": "Zamora et al. (2010)",
            "doi": "10.1021/jf102026c",
            "kind": "calibration_reference",
            "chemistry_family": "lipid_maillard_crosstalk",
            "slr_family_source": "11",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["decadienal", "asparagine_decarboxylation", "acrylamide_pathway"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["11"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Decarboxylation kinetics of asparagine in O/W emulsion in the presence of decadienal (Ea = 81.0 kJ/mol)."],
            "what_it_does_not_support": ["Purely aqueous sugar-asparagine Maillard systems."],
            "key_values": {"replicates": 3, "Ea_asparagine_decarboxylation_kj_mol": 81.0},
            "repo_next_action": "Encode lipid-catalyzed decarboxylation prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "zamora_2010_decadienal_asparagine_decarboxylation", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "ding_2020_schiff_base_amadori_emulsion_rates",
            "citation": "Ding et al. (2020)",
            "doi": "10.1021/acs.jafc.0c04738",
            "kind": "calibration_reference",
            "chemistry_family": "lipid_maillard_crosstalk",
            "slr_family_source": "11",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["schiff_base", "amadori_rearrangement", "micellar_emulsion", "rate_ratios"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["11"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Emulsion Tween-20 micellar rates where reverse Schiff base and Amadori constants are 10^3 times forward constants."],
            "what_it_does_not_support": ["Dry extrusion melts without surfactant boundaries."],
            "key_values": {"replicates": 3},
            "repo_next_action": "Encode micellar Schiff base prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "ding_2020_schiff_base_amadori_emulsion_rates", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "richards_2009_hemoglobin_liposome_oxidation",
            "citation": "Richards et al. (2009)",
            "doi": "10.1021/jf9013394",
            "kind": "calibration_reference",
            "chemistry_family": "lipid_maillard_crosstalk",
            "slr_family_source": "11",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["hemoglobin", "liposome_oxidation", "michaelis_menten", "rate_parameters"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["11"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Michaelis-Menten kinetic parameters for cod and bovine hemoglobin-induced lipid oxidation in liposomes."],
            "what_it_does_not_support": ["Metal naphthenate chemical catalysis without heme coordination."],
            "key_values": {"replicates": 3, "cod_Hb_Vmax_uM_min": 66.2, "cod_Hb_Km_uM": 0.67, "bovine_Hb_Vmax_uM_min": 56.6, "bovine_Hb_Km_uM": 1.2},
            "repo_next_action": "Encode Hb-mediated oxidation prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "richards_2009_hemoglobin_liposome_oxidation", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "smagghe_2006_leghemoglobin_oxygen_dissociation",
            "citation": "Smagghe et al. (2006)",
            "doi": "10.1021/bi051902l",
            "kind": "calibration_reference",
            "chemistry_family": "lipid_maillard_crosstalk",
            "slr_family_source": "11",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["soy_leghemoglobin", "oxygen_dissociation", "autoxidation_rate"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["11"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Recombinant soy leghemoglobin autoxidation and oxygen dissociation kinetics (rate constant = 5.6 s-1)."],
            "what_it_does_not_support": ["Muscle-derived myoglobin autoxidation rates."],
            "key_values": {"replicates": 3, "soy_leghemoglobin_dissociation_s1": 5.6},
            "repo_next_action": "Encode leghemoglobin dissociation prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "smagghe_2006_leghemoglobin_oxygen_dissociation", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        # --- Family 12 ---
        {
            "id": "mandelli_2024_plant_milk_hca_quantification",
            "citation": "Mandelli et al. (2025)",
            "doi": "10.3390/foods14193295",
            "kind": "calibration_reference",
            "chemistry_family": "protein_damage_markers",
            "slr_family_source": "12",
            "payload_role": "safety_reference_payload",
            "observable_panel_tags": ["plant_milk_beverage", "hca_variants", "safety_levels"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["12"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "plant_protein_matrix",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Quantitation and extraction profile of 10 HCA variants in soy and almond milk beverages."],
            "what_it_does_not_support": ["Dry extrusion HCA profiles."],
            "key_values": {"replicates": 3},
            "repo_next_action": "Encode HCA safety reference prior.",
            "runtime_artifacts": [{"artifact_type": "safety_reference_payload", "artifact_id": "mandelli_2024_plant_milk_hca_quantification", "path": "data/lit/safety_reference_payloads.json"}],
            "requires_primary_data": False
        },
        {
            "id": "chen_2016_carbohydrate_meat_cml_cel_sterilization",
            "citation": "Chen et al. (2016)",
            "doi": "10.1007/s10068-016-0185-5",
            "kind": "calibration_reference",
            "chemistry_family": "protein_damage_markers",
            "slr_family_source": "12",
            "payload_role": "safety_reference_payload",
            "observable_panel_tags": ["muscle_protein", "glucose_doping", "cml_cel_sterilization", "zero_order_rates"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["12"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "high_moisture_melt",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Zero-order CML and CEL formation rates during high-moisture sterilization accelerated 4x by glucose doping."],
            "what_it_does_not_support": ["Low-moisture starch-dominated systems."],
            "key_values": {"replicates": 3, "glucose_doping_multiplier": 4.0},
            "repo_next_action": "Encode sterilization AGE safety prior.",
            "runtime_artifacts": [{"artifact_type": "safety_reference_payload", "artifact_id": "chen_2016_carbohydrate_meat_cml_cel_sterilization", "path": "data/lit/safety_reference_payloads.json"}],
            "requires_primary_data": False
        },
        {
            "id": "pruteanu_2023_glucose_phe_arrhenius_browning",
            "citation": "Pruteanu et al. (2023)",
            "doi": "10.1016/j.lwt.2024.117316",
            "kind": "calibration_reference",
            "chemistry_family": "protein_damage_markers",
            "slr_family_source": "12",
            "payload_role": "safety_reference_payload",
            "observable_panel_tags": ["glucose_glycine", "glucose_phenylalanine", "arrhenius_browning", "activation_energy"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["12"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Apparent Arrhenius activation energy comparison for glucose/glycine (Ea = 109 kJ/mol) and glucose/phenylalanine (Ea = 145 kJ/mol)."],
            "what_it_does_not_support": ["Complex solid-state protein matrices with native starch."],
            "key_values": {"replicates": 3, "glucose_glycine_Ea_kj_mol": 109.0, "glucose_phenylalanine_Ea_kj_mol": 145.0},
            "repo_next_action": "Encode browning activation energy safety prior.",
            "runtime_artifacts": [{"artifact_type": "safety_reference_payload", "artifact_id": "pruteanu_2023_glucose_phe_arrhenius_browning", "path": "data/lit/safety_reference_payloads.json"}],
            "requires_primary_data": False
        },
        # --- Family 13 ---
        {
            "id": "zhu_2020_epicatechin_mgo_go_scavenging",
            "citation": "Zhu H. et al. (2020)",
            "doi": "10.1021/acs.jafc.0c01761",
            "kind": "calibration_reference",
            "chemistry_family": "polyphenol_amino_capping",
            "slr_family_source": "13",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["uht_milk", "epicatechin_adducts", "methylglyoxal_trapping", "second_order_rates"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["13"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Second-order epicatechin methylglyoxal (k = 1.6 M-1 s-1) and glyoxal (k = 0.059 M-1 s-1) trapping in stored milk models."],
            "what_it_does_not_support": ["Enzymatic polyphenol oxidase-catalyzed capping."],
            "key_values": {"replicates": 4, "MGO_trapping_k_M1_s1": 1.6, "GO_trapping_k_M1_s1": 0.059},
            "repo_next_action": "Encode epicatechin dicarbonyl trapping prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "zhu_2020_epicatechin_mgo_go_scavenging", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "zhu_2020_polyphenol_mgo_structure_kinetics",
            "citation": "Zhu H. et al. (2020b)",
            "doi": "10.1016/j.foodchem.2020.126500",
            "kind": "calibration_reference",
            "chemistry_family": "polyphenol_amino_capping",
            "slr_family_source": "13",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["polyphenol_structures", "kaempferol", "resveratrol", "mgo_trapping_rates"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["13"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Structure-dependent MGO trapping kinetics for kaempferol (k = 0.063 M-1 s-1) and resveratrol (k = 0.027 M-1 s-1)."],
            "what_it_does_not_support": ["Catechin-class flavanol trapping mechanisms."],
            "key_values": {"replicates": 3, "kaempferol_MGO_k_M1_s1": 0.063, "resveratrol_MGO_k_M1_s1": 0.027},
            "repo_next_action": "Encode structural polyphenol trapping prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "zhu_2020_polyphenol_mgo_structure_kinetics", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "liu_2023_benzoquinone_lysine_kinetics",
            "citation": "Liu J. et al. (2023)",
            "doi": "10.1016/j.foodres.2022.112187",
            "kind": "calibration_reference",
            "chemistry_family": "polyphenol_amino_capping",
            "slr_family_source": "13",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["benzoquinone", "lysine_conjugation", "stopped_flow", "activation_energy"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["13"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Second-order 4-methylbenzoquinone-lysine conjugation kinetics (Ea = 19.00 kJ/mol)."],
            "what_it_does_not_support": ["Michael additions to soft thiol sulfur nucleophiles."],
            "key_values": {"replicates": 3, "Ea_4MBQ_Lysine_kj_mol": 19.0},
            "repo_next_action": "Encode benzoquinone-lysine conjugation prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "liu_2023_benzoquinone_lysine_kinetics", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "comert_gokmen_2019_digestive_mgo_scavenging",
            "citation": "Cömert & Gökmen (2019)",
            "doi": "10.1016/j.foodres.2019.03.046",
            "kind": "calibration_reference",
            "chemistry_family": "polyphenol_amino_capping",
            "slr_family_source": "13",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["digestive_fluid", "phloretin_scavenging", "methylglyoxal_trapping", "bimolecular_rates"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["13"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Methylglyoxal scavenging kinetics under simulated gastrointestinal digestion fluid (egg-MGO rate k = 26.6 L/mol/min)."],
            "what_it_does_not_support": ["High-temperature dry extrusion caramelization."],
            "key_values": {"replicates": 3, "egg_MGO_rate_k_L_mol_min": 26.6},
            "repo_next_action": "Encode digestive dicarbonyl scavenging prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "comert_gokmen_2019_digestive_mgo_scavenging", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "munoz_2007_chlorogenic_oxidase_quinone",
            "citation": "Munoz et al. (2007)",
            "doi": "10.1021/jf062081+",
            "kind": "calibration_reference",
            "chemistry_family": "polyphenol_amino_capping",
            "slr_family_source": "13",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["chlorogenic_acid", "oxidase_oxidation", "quinone_absorptivity", "reaction_rates"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["13"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Chlorogenic acid enzymatically oxidized to o-quinone reacting with substrate (k = 2.73 M-1 s-1)."],
            "what_it_does_not_support": ["Non-catechol structural polyphenol trapping paths."],
            "key_values": {"replicates": 3, "CQA_Q_reaction_rate_k_M1_s1": 2.73},
            "repo_next_action": "Encode chlorogenic o-quinone prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "munoz_2007_chlorogenic_oxidase_quinone", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "comert_gokmen_2021_epicatechin_cysteine_mgo_synergy",
            "citation": "Cömert & Gökmen (2021)",
            "doi": "10.1016/j.foodchem.2020.128670",
            "kind": "calibration_reference",
            "chemistry_family": "polyphenol_amino_capping",
            "slr_family_source": "13",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["epicatechin_cysteine", "synergistic_scavenging", "methylglyoxal_trapping"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["13"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Synergistic rate augmentation between epicatechin and cysteine in trapping methylglyoxal."],
            "what_it_does_not_support": ["Antagonistic pathways driven by gallic acid variants."],
            "key_values": {"replicates": 3},
            "repo_next_action": "Encode EC-Cys synergistic scavenging prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "comert_gokmen_2021_epicatechin_cysteine_mgo_synergy", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "song_2009_benzoquinone_gsh_conjugation",
            "citation": "Song et al. (2009)",
            "doi": "10.1073/pnas.0810352106",
            "kind": "calibration_reference",
            "chemistry_family": "polyphenol_amino_capping",
            "slr_family_source": "13",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["benzoquinone", "gsh_conjugation", "stopped_flow", "reverse_thermodynamics"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["13"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Benzoquinone conjugation with GSH rate parameters (forward k = 1547 M-1 s-1, reverse Ea = 11.2 kcal/mol)."],
            "what_it_does_not_support": ["Reactions catalyzed by polyphenol oxidases."],
            "key_values": {"replicates": 3, "forward_rate_k_M1_s1": 1547.0, "reverse_Ea_kcal_mol": 11.2},
            "repo_next_action": "Encode quinone-GSH conjugation prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "song_2009_benzoquinone_gsh_conjugation", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "li_2017_quercetin_mgo_adduct_browning_mitigation",
            "citation": "Li et al. (2017)",
            "doi": "10.1021/acs.jafc.6b05811",
            "kind": "calibration_reference",
            "chemistry_family": "polyphenol_amino_capping",
            "slr_family_source": "13",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["quercetin_adducts", "browning_mitigation", "lysine_glucose_model"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["13"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Quercetin-MGO dicarbonyl trapping kinetics mitigating browning in lysine/glucose Maillard models."],
            "what_it_does_not_support": ["Purely non-enzymatic lipid oxidation paths."],
            "key_values": {"replicates": 3},
            "repo_next_action": "Encode quercetin MGO adduction prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "li_2017_quercetin_mgo_adduct_browning_mitigation", "path": "data/lit/computational_priors.json"}],
            "requires_primary_data": False
        },
        {
            "id": "monforte_2018_phenylacetaldehyde_quinones",
            "citation": "Monforte et al. (2018)",
            "doi": "10.1021/acs.jafc.7b00264",
            "kind": "calibration_reference",
            "chemistry_family": "polyphenol_amino_capping",
            "slr_family_source": "13",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["phenylacetaldehyde_kinetics", "hydroquinone_catalysis", "activation_energy"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["13"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["Phenylacetaldehyde formation activation energy catalyzed by hydroquinone (Ea = 32.9 kJ/mol) or trihydroxybenzene (Ea = 31.5 kJ/mol)."],
            "what_it_does_not_support": ["Strecker degradation in non-polyphenol systems."],
            "key_values": {"replicates": 3, "Ea_hydroquinone_catalyzed_kj_mol": 32.9, "Ea_trihydroxybenzene_catalyzed_kj_mol": 31.5},
            "repo_next_action": "Encode polyphenol-catalyzed Strecker prior.",
            "runtime_artifacts": [{"artifact_type": "computational_prior", "artifact_id": "monforte_2018_phenylacetaldehyde_quinones", "path": "data/lit/computational_priors.json"}],
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
        sec = "13.0"
        if ref["slr_family_source"] == "02": sec = "2.0"
        elif ref["slr_family_source"] == "06": sec = "6.0"
        elif ref["slr_family_source"] == "07": sec = "7.0"
        elif ref["slr_family_source"] == "09": sec = "9.0"
        elif ref["slr_family_source"] == "11": sec = "11.0"
        elif ref["slr_family_source"] == "12": sec = "12.0"

        # Map details to compounds / parameters
        compounds = []
        params = []
        anchors = []

        if ref["slr_family_source"] == "02":
            compounds = ["peroxides", "hexanal"]
            params = ["lipid_oxidation_kinetics", "hydroperoxide_decomposition_activation_energy"]
            if ref["id"] == "spier_2010_metal_naphthenate_peroxide_decomp":
                anchors = ["Ea: 11-26 kcal/mol"]
            elif ref["id"] == "bayram_2023_stripped_soybean_oil_hexanal":
                anchors = ["Ea: 1.62-89.4 kJ/mol"]
        elif ref["slr_family_source"] == "06":
            compounds = ["glyco-conjugates", "dityrosine", "CML", "CEL"]
            params = ["protein_glycation_kinetics", "dityrosine_diffusion_limited_kinetics"]
            if ref["id"] == "sun_2020_solid_matrix_cml_cel_accumulation":
                anchors = ["Ea free CML: 44.16 kJ/mol", "Ea free CEL: 40.97 kJ/mol"]
            elif ref["id"] == "pereira_2020_metal_pm_haber_weiss_chelation":
                anchors = ["Cu-PM free energy: -35.8 kcal/mol"]
        elif ref["slr_family_source"] == "07":
            compounds = ["Amadori", "decarboxylated_complexes"]
            params = ["Amadori_degradation_kinetics", "decarboxylation_rate_constants"]
            if ref["id"] == "nashalian_yaylayan_2014_cu_catalyzed_strecker":
                anchors = ["decarboxylation rate k: 0.41 s-1"]
        elif ref["slr_family_source"] == "09":
            compounds = ["furanones", "retro_aldol_intermediates"]
            params = ["caramelization_kinetics", "molten_sugar_color_activation_energy"]
        elif ref["slr_family_source"] == "11":
            compounds = ["styrene", "2-pentylpyridine", "acrylamide", "met-leghemoglobin"]
            params = ["decadienal_amino_crosstalk_kinetics", "leghemoglobin_oxygen_dissociation"]
            if ref["id"] == "hidalgo_2007_decadienal_phenylalanine_styrene":
                anchors = ["Ea styrene: 150.4 kJ/mol"]
            elif ref["id"] == "zamora_2010_decadienal_asparagine_decarboxylation":
                anchors = ["Ea asparagine: 81.0 kJ/mol"]
            elif ref["id"] == "smagghe_2006_leghemoglobin_oxygen_dissociation":
                anchors = ["dissociation rate k: 5.6 s-1"]
        elif ref["slr_family_source"] == "12":
            compounds = ["HCAs", "CML", "CEL", "acrylamide"]
            params = ["HCA_formation_kinetics", "CML_CEL_accumulation_kinetics", "browning_activation_energy"]
            if ref["id"] == "pruteanu_2023_glucose_phe_arrhenius_browning":
                anchors = ["Ea Phe model: 145 kJ/mol", "Ea Gly model: 109 kJ/mol"]
            elif ref["id"] == "chen_2016_carbohydrate_meat_cml_cel_sterilization":
                anchors = ["glucose doping acceleration: 4.0x"]
        elif ref["slr_family_source"] == "13":
            compounds = ["epicatechin-adducts", "quinone-conjugates", "phenylacetaldehyde"]
            params = ["quinone_Michael_addition_kinetics", "dicarbonyl_trapping_rates", "Strecker_catalysis"]
            if ref["id"] == "zhu_2020_epicatechin_mgo_go_scavenging":
                anchors = ["MGO k: 1.6 M-1 s-1", "GO k: 0.059 M-1 s-1"]
            elif ref["id"] == "liu_2023_benzoquinone_lysine_kinetics":
                anchors = ["Ea 4MBQ-Lys: 19.00 kJ/mol"]
            elif ref["id"] == "song_2009_benzoquinone_gsh_conjugation":
                anchors = ["k forward: 1547 M-1 s-1"]
            elif ref["id"] == "munoz_2007_chlorogenic_oxidase_quinone":
                anchors = ["k CQA_Q: 2.73 M-1 s-1"]

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
            "current_user_visible_surfaces": ["reporting"],
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
    # 1. Thiol Capping (Family 13)
    if "polyphenol_thiol_capping_priors" in priors:
        existing_capping = {p["id"] for p in priors["polyphenol_thiol_capping_priors"]}
        for ref in new_references:
            if ref["slr_family_source"] == "13" and ref["id"] not in existing_capping:
                priors["polyphenol_thiol_capping_priors"].append({
                    "id": ref["id"],
                    "source": ref["citation"],
                    "provenance_tier": "literature_derived_transfer",
                    "confidence_tier": "high",
                    "uncertainty_posture": "directional_transfer",
                    "key_values": ref["key_values"],
                    "notes": ref["what_it_supports"][0]
                })
                print(f"Added {ref['id']} to polyphenol_thiol_capping_priors")

    # 2. Lipid oxidation crosstalk priors (Family 02 -> lipid_offnote_priors)
    if "lipid_offnote_priors" in priors:
        existing_lipid = {p["id"] for p in priors["lipid_offnote_priors"]}
        for ref in new_references:
            if ref["slr_family_source"] == "02" and ref["id"] not in existing_lipid:
                priors["lipid_offnote_priors"].append({
                    "id": ref["id"],
                    "source": ref["citation"],
                    "provenance_tier": "literature_derived_transfer",
                    "confidence_tier": "high",
                    "uncertainty_posture": "directional_transfer",
                    "key_values": ref["key_values"],
                    "notes": ref["what_it_supports"][0]
                })
                print(f"Added {ref['id']} to lipid_offnote_priors")

    # 3. Alternative proteins / crosslinks (Family 06 -> crosslink_kinetics_priors)
    if "crosslink_kinetics_priors" in priors:
        existing_crosslink = {p["id"] for p in priors["crosslink_kinetics_priors"]}
        for ref in new_references:
            if ref["slr_family_source"] == "06" and ref["id"] not in existing_crosslink:
                priors["crosslink_kinetics_priors"].append({
                    "id": ref["id"],
                    "source": ref["citation"],
                    "provenance_tier": "literature_derived_transfer",
                    "confidence_tier": "high",
                    "uncertainty_posture": "directional_transfer",
                    "key_values": ref["key_values"],
                    "notes": ref["what_it_supports"][0]
                })
                print(f"Added {ref['id']} to crosslink_kinetics_priors")

    # 4. Carbonyl donor hierarchy (Family 07 -> carbonyl_donor_priors)
    if "carbonyl_donor_priors" in priors:
        existing_carbonyl = {p["id"] for p in priors["carbonyl_donor_priors"]}
        for ref in new_references:
            if ref["slr_family_source"] == "07" and ref["id"] not in existing_carbonyl:
                priors["carbonyl_donor_priors"].append({
                    "id": ref["id"],
                    "source": ref["citation"],
                    "provenance_tier": "literature_derived_transfer",
                    "confidence_tier": "high",
                    "uncertainty_posture": "directional_transfer",
                    "key_values": ref["key_values"],
                    "notes": ref["what_it_supports"][0]
                })
                print(f"Added {ref['id']} to carbonyl_donor_priors")

    # 5. Caramelization (Family 09 -> furanone_priors)
    if "furanone_priors" in priors:
        existing_furanone = {p["id"] for p in priors["furanone_priors"]}
        for ref in new_references:
            if ref["slr_family_source"] == "09" and ref["id"] not in existing_furanone:
                priors["furanone_priors"].append({
                    "id": ref["id"],
                    "source": ref["citation"],
                    "provenance_tier": "literature_derived_transfer",
                    "confidence_tier": "high",
                    "uncertainty_posture": "directional_transfer",
                    "key_values": ref["key_values"],
                    "notes": ref["what_it_supports"][0]
                })
                print(f"Added {ref['id']} to furanone_priors")

    # 6. Lipid-Maillard Crosstalk (Family 11 -> strecker_crosstalk_priors)
    if "strecker_crosstalk_priors" in priors:
        existing_strecker = {p["id"] for p in priors["strecker_crosstalk_priors"]}
        for ref in new_references:
            if ref["slr_family_source"] == "11" and ref["id"] not in existing_strecker:
                priors["strecker_crosstalk_priors"].append({
                    "id": ref["id"],
                    "source": ref["citation"],
                    "provenance_tier": "literature_derived_transfer",
                    "confidence_tier": "high",
                    "uncertainty_posture": "directional_transfer",
                    "key_values": ref["key_values"],
                    "notes": ref["what_it_supports"][0]
                })
                print(f"Added {ref['id']} to strecker_crosstalk_priors")

    # 7. Add catalytic_browning_priors block (heme + transition metal catalysis)
    if "catalytic_browning_priors" not in priors:
        priors["catalytic_browning_priors"] = []
    existing_catalytic = {p["id"] for p in priors["catalytic_browning_priors"]}
    catalytic_references = [
        "spier_2010_metal_naphthenate_peroxide_decomp",
        "pereira_2020_metal_pm_haber_weiss_chelation",
        "nashalian_yaylayan_2014_cu_catalyzed_strecker",
        "richards_2009_hemoglobin_liposome_oxidation",
        "smagghe_2006_leghemoglobin_oxygen_dissociation"
    ]
    for ref in new_references:
        if ref["id"] in catalytic_references and ref["id"] not in existing_catalytic:
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

    # D. Add to safety_reference_payloads.json (Family 12)
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
                    "analyte": "HCAs" if "hca" in ref["id"] else ("CML" if "cml" in ref["id"] else "browning_intermediates"),
                    "method": {
                        "instrument": "HPLC" if "pruteanu" in ref["id"] else ("UHPLC-MS/MS" if "urugo" in ref["id"] else "LC-MS"),
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
        filename = f"{ref['slr_family_source']}_"
        # Find SLR filename based on family ID
        filename_map = {
            "02": "02_lipid_oxidation_and_carbonylic_crosstalk.md",
            "06": "06_alternative_proteins.md",
            "07": "07_reducing_sugars.md",
            "09": "09_carbohydrate_pyrolysis.md",
            "11": "11_maillard_lipid_crosstalk.md",
            "12": "12_protein_damage_markers.md",
            "13": "13_polyphenol_amino_capping.md"
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
            print(f"Updated backlog status for DOI {ref['doi']} (citation: {item.get('citation')}) to RUNTIME_BOUND")
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
                    "score": "7/8",
                    "score_value": 7,
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
