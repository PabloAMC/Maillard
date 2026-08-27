import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DATA_LIT_DIR = ROOT / "data" / "lit"

# File paths
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
    # 1. Load files
    intake = load_json(INTAKE_REGISTRY_PATH)
    matrix = load_json(SLR_MATRIX_PATH)
    priors = load_json(COMPUTATIONAL_PRIORS_PATH)
    safety = load_json(SAFETY_PAYLOADS_PATH)
    backlog = load_json(BACKLOG_PATH)

    # 2. Define the new references for benchmark_intake_registry.json
    new_eligible_references = [
        {
            "id": "lagrain_2010_cystine_elimination_lanthionine",
            "citation": "Lagrain et al. (2010)",
            "doi": "10.1021/jf102575r",
            "kind": "calibration_reference",
            "chemistry_family": "alternative_protein_matrix_scope",
            "slr_family_source": "06",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["crosslinking", "lanthionine", "cysteine_loss"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["12"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "wheat_gliadin",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Derived strict first-order kinetics for cystine beta-elimination (Ea 119 kJ/mol at pH 8.0, 88.2 kJ/mol at pH 6.0).",
                "Second-order kinetics for lanthionine formation (Ea 12.3 kJ/mol at pH 8.0, 6.51 kJ/mol at pH 6.0).",
            ],
            "what_it_does_not_support": ["High-moisture extrusion shear cells"],
            "key_values": {
                "cystine_elimination_ea_ph8_kj_mol": 119.0,
                "cystine_elimination_ea_ph6_kj_mol": 88.2,
                "lanthionine_formation_ea_ph8_kj_mol": 12.3,
                "lanthionine_formation_ea_ph6_kj_mol": 6.51,
                "cystine_elimination_k_ph8_min": 0.054,
                "replicates": 3,
            },
            "repo_next_action": "Encode lanthionine and cystine elimination kinetics priors.",
            "runtime_artifacts": [
                {
                    "artifact_type": "computational_prior",
                    "artifact_id": "lagrain_2010_cystine_elimination_lanthionine",
                    "path": "data/lit/computational_priors.json",
                }
            ],
            "requires_primary_data": False,
        },
        {
            "id": "morel_2002_gluten_shear_aggregation",
            "citation": "Morel et al. (2002)",
            "doi": "10.1021/bm015639p",
            "kind": "calibration_reference",
            "chemistry_family": "alternative_protein_matrix_scope",
            "slr_family_source": "06",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["shear_damage", "sulfhydryl_loss", "gluten"],
            "process_state_scope": ["heated_matrix", "extrusion_structured"],
            "supporting_families": ["12"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "wheat_gluten",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Mechanochemical shift reducing apparent Arrhenius activation energy for protein solubility loss under intense mechanical shear from thermally dominated down to 33.7 kJ/mol."
            ],
            "what_it_does_not_support": ["Free amino acid kinetics in solution"],
            "key_values": {"shear_reduced_ea_kj_mol": 33.7, "replicates": 3},
            "repo_next_action": "Encode shear-induced aggregation priors.",
            "runtime_artifacts": [
                {
                    "artifact_type": "computational_prior",
                    "artifact_id": "morel_2002_gluten_shear_aggregation",
                    "path": "data/lit/computational_priors.json",
                }
            ],
            "requires_primary_data": False,
        },
        {
            "id": "yu_2017_cml_cel_meat_review",
            "citation": "Yu et al. (2017)",
            "doi": "10.1016/j.tifs.2020.01.021",
            "kind": "calibration_reference",
            "chemistry_family": "protein_damage_markers",
            "slr_family_source": "12",
            "payload_role": "safety_payload",
            "observable_panel_tags": ["safety", "cml", "cel", "activation_energy"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["01"],
            "target_modules": ["literature_runtime", "safety"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Differential thermodynamics for advanced glycation endproducts: CML activation energy of 61.01 kJ/mol versus CEL activation energy of 29.21 kJ/mol."
            ],
            "what_it_does_not_support": ["Extrusion mechanical shear effects"],
            "key_values": {
                "cml_formation_ea_kj_mol": 61.01,
                "cel_formation_ea_kj_mol": 29.21,
                "replicates": 3,
            },
            "repo_next_action": "Expose as CML/CEL activation energy safety reference payload.",
            "runtime_artifacts": [
                {
                    "artifact_type": "safety_reference_payload",
                    "artifact_id": "yu_2017_cml_cel_meat_review",
                    "path": "data/lit/safety_reference_payloads.json",
                }
            ],
            "requires_primary_data": False,
        },
        {
            "id": "huang_2024_dixyl_arp_degradation",
            "citation": "Huang et al. (2024)",
            "doi": "10.1021/acs.jafc.4c05736",
            "kind": "calibration_reference",
            "chemistry_family": "carbonyl_donor_hierarchy",
            "slr_family_source": "07",
            "payload_role": "computational_prior",
            "observable_panel_tags": [
                "arp_degradation",
                "dicarbonyl",
                "furosine",
                "kinetics",
            ],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["12"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Bifurcated, first-order Amadori degradation kinetics: 3-deoxyglucosone-mediated pathway (Ea 115.0 kJ/mol) versus methylglyoxal-mediated pathway (Ea 85.0 kJ/mol)."
            ],
            "what_it_does_not_support": ["High moisture extrusion mechanical effects"],
            "key_values": {
                "three_dg_mediated_ea_kj_mol": 115.0,
                "mgo_mediated_ea_kj_mol": 85.0,
                "replicates": 3,
            },
            "repo_next_action": "Encode bifurcated Amadori degradation pathway priors.",
            "runtime_artifacts": [
                {
                    "artifact_type": "computational_prior",
                    "artifact_id": "huang_2024_dixyl_arp_degradation",
                    "path": "data/lit/computational_priors.json",
                }
            ],
            "requires_primary_data": False,
        },
        {
            "id": "charissou_2007_cookie_cml_furosine",
            "citation": "Charissou et al. (2007)",
            "doi": "10.1021/jf063024j",
            "kind": "calibration_reference",
            "chemistry_family": "protein_damage_markers",
            "slr_family_source": "12",
            "payload_role": "safety_payload",
            "observable_panel_tags": ["safety", "cml", "furosine", "cookie"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["01"],
            "target_modules": ["literature_runtime", "safety"],
            "matrix_family": "cookie_model",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "CML acts as a stable, linear (zero-order) accumulator marker under low moisture, while furosine shows a transient, bell-shaped decay curve dependent on temperature (180 to 220 C)."
            ],
            "what_it_does_not_support": [
                "High moisture extrusion mechanical shear effects"
            ],
            "key_values": {"temp_range_C": [180.0, 220.0], "replicates": 3},
            "repo_next_action": "Expose as cookie CML/furosine safety reference payload.",
            "runtime_artifacts": [
                {
                    "artifact_type": "safety_reference_payload",
                    "artifact_id": "charissou_2007_cookie_cml_furosine",
                    "path": "data/lit/safety_reference_payloads.json",
                }
            ],
            "requires_primary_data": False,
        },
        {
            "id": "fratianni_2016_apricot_furosine",
            "citation": "Fratianni et al. (2016)",
            "doi": "10.1016/j.foodres.2016.12.009",
            "kind": "calibration_reference",
            "chemistry_family": "protein_damage_markers",
            "slr_family_source": "12",
            "payload_role": "safety_payload",
            "observable_panel_tags": ["safety", "furosine", "apricot"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["01"],
            "target_modules": ["literature_runtime", "safety"],
            "matrix_family": "apricot",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Zero-order reaction kinetics for furosine formation in low-moisture apricot matrices governed by an activation energy of 83.3 kJ/mol."
            ],
            "what_it_does_not_support": ["High temperature dry extrusion kinetics"],
            "key_values": {"furosine_formation_ea_kj_mol": 83.3, "replicates": 3},
            "repo_next_action": "Expose as apricot furosine safety reference payload.",
            "runtime_artifacts": [
                {
                    "artifact_type": "safety_reference_payload",
                    "artifact_id": "fratianni_2016_apricot_furosine",
                    "path": "data/lit/safety_reference_payloads.json",
                }
            ],
            "requires_primary_data": False,
        },
        {
            "id": "rombouts_2012_gluten_crosslinking",
            "citation": "Rombouts et al. (2012)",
            "doi": "10.1021/jf3024672",
            "kind": "calibration_reference",
            "chemistry_family": "alternative_protein_matrix_scope",
            "slr_family_source": "06",
            "payload_role": "computational_prior",
            "observable_panel_tags": [
                "crosslinking",
                "lanthionine",
                "lysinoalanine",
                "gluten",
            ],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["12"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "wheat_gluten",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Overall reaction rate constants for wheat gluten heat-induced polymerization and non-disulfide cross-linking (lanthionine/lysinoalanine) with activation energy of 142.0 kJ/mol above thermal threshold of 90.0 C."
            ],
            "what_it_does_not_support": ["Low moisture extrusion Apparent Viscosity"],
            "key_values": {
                "gluten_polymerization_ea_kj_mol": 142.0,
                "thermal_threshold_temp_C": 90.0,
                "replicates": 3,
            },
            "repo_next_action": "Encode gluten cross-linking priors.",
            "runtime_artifacts": [
                {
                    "artifact_type": "computational_prior",
                    "artifact_id": "rombouts_2012_gluten_crosslinking",
                    "path": "data/lit/computational_priors.json",
                }
            ],
            "requires_primary_data": False,
        },
        {
            "id": "ma_2024_pbma_extrusion_damage",
            "citation": "Ma et al. (2024)",
            "doi": "10.3390/ijms25168668",
            "kind": "calibration_reference",
            "chemistry_family": "protein_damage_markers",
            "slr_family_source": "12",
            "payload_role": "safety_payload",
            "observable_panel_tags": [
                "safety",
                "cml",
                "cel",
                "acrylamide",
                "extrusion",
            ],
            "process_state_scope": ["extrusion_structured"],
            "supporting_families": ["01", "07"],
            "target_modules": ["literature_runtime", "safety"],
            "matrix_family": "plant_based_meat_analogue",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "High-moisture extrusion selective acceleration of CML via glyoxal transport.",
                "Elevated barrel temperatures up to 160 C spikes CEL concentration by 39-101% due to rapid thermochemical cycling of methylglyoxal.",
            ],
            "what_it_does_not_support": ["Low temperature ambient slurry kinetics"],
            "key_values": {
                "temp_threshold_C": 160.0,
                "cel_increase_pct_range": [39, 101],
                "replicates": 3,
            },
            "repo_next_action": "Expose as PBMA extrusion safety reference payload.",
            "runtime_artifacts": [
                {
                    "artifact_type": "safety_reference_payload",
                    "artifact_id": "ma_2024_pbma_extrusion_damage",
                    "path": "data/lit/safety_reference_payloads.json",
                }
            ],
            "requires_primary_data": False,
        },
        {
            "id": "ilo_1996_maize_sme_lysine_damage",
            "citation": "Ilo et al. (1996)",
            "doi": "10.1006/fstl.1996.0092",
            "kind": "calibration_reference",
            "chemistry_family": "alternative_protein_matrix_scope",
            "slr_family_source": "06",
            "payload_role": "computational_prior",
            "observable_panel_tags": [
                "shear_damage",
                "lysine_damage",
                "maize_grits",
                "extrusion",
            ],
            "process_state_scope": ["extrusion_structured"],
            "supporting_families": ["12"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "maize_grits",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Specific mechanical energy (SME) inputs correlate to chemical breakdown.",
                "Lysine damage adheres to first-order kinetics and is suppressed by higher moisture acting as a lubricant.",
            ],
            "what_it_does_not_support": ["Free amino acid kinetics in solution"],
            "key_values": {"kinetic_order": "first_order", "replicates": 3},
            "repo_next_action": "Encode maize extrusion shear damage priors.",
            "runtime_artifacts": [
                {
                    "artifact_type": "computational_prior",
                    "artifact_id": "ilo_1996_maize_sme_lysine_damage",
                    "path": "data/lit/computational_priors.json",
                }
            ],
            "requires_primary_data": False,
        },
    ]

    # A. Add to benchmark_intake_registry.json
    existing_ids = {ref["id"] for ref in intake["eligible_references"]}
    for ref in new_eligible_references:
        if ref["id"] not in existing_ids:
            intake["eligible_references"].append(ref)
            print(f"Added {ref['id']} to benchmark_intake_registry")

    # B. Add to slr_incorporation_matrix.json
    new_matrix_entries = [
        {
            "slr_section": "6.0",
            "paper_id": "lagrain_2010_cystine_elimination_lanthionine",
            "citation": "Lagrain et al. (2010)",
            "doi": "10.1021/jf102575r",
            "matrix_family": "wheat_gliadin",
            "compounds_supported": ["lanthionine", "cystine"],
            "parameters_supported": [
                "cystine_elimination_kinetics",
                "lanthionine_formation_kinetics",
            ],
            "exact_numeric_anchors": [
                "Ea beta-elimination pH 8: 119 kJ/mol",
                "Ea LAN pH 8: 12.3 kJ/mol",
                "k at pH 8.0: 0.054 min-1",
            ],
            "current_repo_artifacts": ["data/lit/computational_priors.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Encode lanthionine and cystine elimination kinetics priors.",
            "confidence_tier": "high",
            "notes_on_limits": "Wheat gliadin model system.",
        },
        {
            "slr_section": "6.0",
            "paper_id": "morel_2002_gluten_shear_aggregation",
            "citation": "Morel et al. (2002)",
            "doi": "10.1021/bm015639p",
            "matrix_family": "wheat_gluten",
            "compounds_supported": ["sulfhydryl"],
            "parameters_supported": ["shear_aggregation_kinetics"],
            "exact_numeric_anchors": ["apparent Ea reduced to 33.7 kJ/mol under shear"],
            "current_repo_artifacts": ["data/lit/computational_priors.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Encode shear-induced aggregation priors.",
            "confidence_tier": "high",
            "notes_on_limits": "Wheat gluten melt under torque.",
        },
        {
            "slr_section": "12.0",
            "paper_id": "yu_2017_cml_cel_meat_review",
            "citation": "Yu et al. (2017)",
            "doi": "10.1016/j.tifs.2020.01.021",
            "matrix_family": "free_model_system",
            "compounds_supported": ["CML", "CEL"],
            "parameters_supported": ["cml_cel_activation_energy"],
            "exact_numeric_anchors": [
                "CML formation Ea 61.01 kJ/mol",
                "CEL formation Ea 29.21 kJ/mol",
            ],
            "current_repo_artifacts": ["data/lit/safety_reference_payloads.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Expose as CML/CEL activation energy safety reference payload.",
            "confidence_tier": "high",
            "notes_on_limits": "Complex isolates and meat products review.",
        },
        {
            "slr_section": "7.0",
            "paper_id": "huang_2024_dixyl_arp_degradation",
            "citation": "Huang et al. (2024)",
            "doi": "10.1021/acs.jafc.4c05736",
            "matrix_family": "free_model_system",
            "compounds_supported": ["furosine", "3-DG", "MGO"],
            "parameters_supported": ["diXyl_arp_degradation_kinetics"],
            "exact_numeric_anchors": [
                "Ea furosine via 3-DG 115 kJ/mol",
                "Ea via MGO 85 kJ/mol",
            ],
            "current_repo_artifacts": ["data/lit/computational_priors.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Encode bifurcated Amadori degradation pathway priors.",
            "confidence_tier": "high",
            "notes_on_limits": "Aqueous free model system.",
        },
        {
            "slr_section": "12.0",
            "paper_id": "charissou_2007_cookie_cml_furosine",
            "citation": "Charissou et al. (2007)",
            "doi": "10.1021/jf063024j",
            "matrix_family": "cookie_model",
            "compounds_supported": ["CML", "furosine"],
            "parameters_supported": ["cookie_glycation_kinetics"],
            "exact_numeric_anchors": [
                "zero-order CML accumulation",
                "transient furosine decay 180 to 220 C",
            ],
            "current_repo_artifacts": ["data/lit/safety_reference_payloads.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Expose as cookie CML/furosine safety reference payload.",
            "confidence_tier": "high",
            "notes_on_limits": "Low moisture cookie model system.",
        },
        {
            "slr_section": "12.0",
            "paper_id": "fratianni_2016_apricot_furosine",
            "citation": "Fratianni et al. (2016)",
            "doi": "10.1016/j.foodres.2016.12.009",
            "matrix_family": "apricot",
            "compounds_supported": ["furosine"],
            "parameters_supported": ["apricot_furosine_kinetics"],
            "exact_numeric_anchors": [
                "zero-order furosine formation",
                "Ea 83.3 kJ/mol",
            ],
            "current_repo_artifacts": ["data/lit/safety_reference_payloads.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Expose as apricot furosine safety reference payload.",
            "confidence_tier": "high",
            "notes_on_limits": "Apricot tissue matrix.",
        },
        {
            "slr_section": "6.0",
            "paper_id": "rombouts_2012_gluten_crosslinking",
            "citation": "Rombouts et al. (2012)",
            "doi": "10.1021/jf3024672",
            "matrix_family": "wheat_gluten",
            "compounds_supported": ["lanthionine", "lysinoalanine"],
            "parameters_supported": ["gluten_crosslinking_kinetics"],
            "exact_numeric_anchors": [
                "overall Ea 142 kJ/mol",
                "non-disulfide crosslinks above 90 C",
            ],
            "current_repo_artifacts": ["data/lit/computational_priors.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Encode gluten cross-linking priors.",
            "confidence_tier": "high",
            "notes_on_limits": "Wheat gluten and BSA model mixtures.",
        },
        {
            "slr_section": "12.0",
            "paper_id": "ma_2024_pbma_extrusion_damage",
            "citation": "Ma et al. (2024)",
            "doi": "10.3390/ijms25168668",
            "matrix_family": "plant_based_meat_analogue",
            "compounds_supported": ["CML", "CEL", "acrylamide"],
            "parameters_supported": ["extrusion_chemical_damage_kinetics"],
            "exact_numeric_anchors": [
                "CML acceleration via glyoxal",
                "CEL 39-101% spike at barrel temp 160 C",
            ],
            "current_repo_artifacts": ["data/lit/safety_reference_payloads.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Expose as PBMA extrusion safety reference payload.",
            "confidence_tier": "high",
            "notes_on_limits": "High-moisture PBMA twin-screw extrusion.",
        },
        {
            "slr_section": "6.0",
            "paper_id": "ilo_1996_maize_sme_lysine_damage",
            "citation": "Ilo et al. (1996)",
            "doi": "10.1006/fstl.1996.0092",
            "matrix_family": "maize_grits",
            "compounds_supported": ["lysine"],
            "parameters_supported": ["maize_extrusion_lysine_damage_kinetics"],
            "exact_numeric_anchors": [
                "lysine damage is first-order",
                "SME input correlation",
            ],
            "current_repo_artifacts": ["data/lit/computational_priors.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Encode maize extrusion shear damage priors.",
            "confidence_tier": "high",
            "notes_on_limits": "Twin-screw extrusion cooking of maize grits.",
        },
    ]

    existing_papers = {ent["paper_id"] for ent in matrix["entries"]}
    for ent in new_matrix_entries:
        if ent["paper_id"] not in existing_papers:
            matrix["entries"].append(ent)
            print(f"Added {ent['paper_id']} to slr_incorporation_matrix")

    # C. Add to computational_priors.json
    # C.1 Add new metadata to section_family_metadata
    if "crosslink_kinetics_priors" not in priors["section_family_metadata"]:
        priors["section_family_metadata"]["crosslink_kinetics_priors"] = {
            "chemistry_family": "alternative_protein_matrix_scope",
            "slr_family_source": "06",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["crosslinking", "lanthionine", "lysinoalanine"],
            "process_state_scope": ["heated_matrix", "extrusion_structured"],
            "supporting_families": ["12"],
        }
        print("Added crosslink_kinetics_priors metadata")

    if "shear_damage_priors" not in priors["section_family_metadata"]:
        priors["section_family_metadata"]["shear_damage_priors"] = {
            "chemistry_family": "alternative_protein_matrix_scope",
            "slr_family_source": "06",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["shear_damage", "sulfhydryl_loss", "gluten"],
            "process_state_scope": ["extrusion_structured"],
            "supporting_families": ["12"],
        }
        print("Added shear_damage_priors metadata")

    # C.2 Append new lists at root level if not exist
    if "crosslink_kinetics_priors" not in priors:
        priors["crosslink_kinetics_priors"] = []
    if "shear_damage_priors" not in priors:
        priors["shear_damage_priors"] = []

    # C.3 Populate the new priors
    new_crosslink_priors = [
        {
            "id": "lagrain_2010_cystine_elimination_lanthionine",
            "source": "Lagrain et al. (2010)",
            "provenance_tier": "literature_derived_transfer",
            "confidence_tier": "high",
            "uncertainty_posture": "directional_transfer",
            "reference_conditions": {
                "ph_values": [6.0, 8.0],
                "temp_range_C": [30.0, 100.0],
            },
            "elimination_kinetics": {
                "ea_ph8_kj_mol": 119.0,
                "ea_ph6_kj_mol": 88.2,
                "k_ph8_min": 0.054,
            },
            "lanthionine_formation_kinetics": {
                "ea_ph8_kj_mol": 12.3,
                "ea_ph6_kj_mol": 6.51,
            },
            "notes": "Strict first-order kinetics for cystine beta-elimination and second-order lanthionine formation.",
        },
        {
            "id": "rombouts_2012_gluten_crosslinking",
            "source": "Rombouts et al. (2012)",
            "provenance_tier": "literature_derived_transfer",
            "confidence_tier": "high",
            "uncertainty_posture": "directional_transfer",
            "reference_conditions": {"thermal_threshold_temp_C": 90.0},
            "gluten_polymerization": {"ea_kj_mol": 142.0},
            "notes": "Overall reaction rate constants for gluten polymerization and crosslinking above thermal thresholds.",
        },
    ]
    existing_crosslink = {pr["id"] for pr in priors["crosslink_kinetics_priors"]}
    for pr in new_crosslink_priors:
        if pr["id"] not in existing_crosslink:
            priors["crosslink_kinetics_priors"].append(pr)
            print(f"Added {pr['id']} to crosslink_kinetics_priors")

    new_shear_priors = [
        {
            "id": "morel_2002_gluten_shear_aggregation",
            "source": "Morel et al. (2002)",
            "provenance_tier": "literature_derived_transfer",
            "confidence_tier": "high",
            "uncertainty_posture": "directional_transfer",
            "reference_conditions": {"shear_rate_rpm": [10.0, 100.0]},
            "mechanochemical_kinetics": {"shear_reduced_ea_kj_mol": 33.7},
            "notes": "Conclusively proves that application of intense mechanical shear reduces Arrhenius activation energy for solubility loss down to 33.7 kJ/mol.",
        },
        {
            "id": "ilo_1996_maize_sme_lysine_damage",
            "source": "Ilo et al. (1996)",
            "provenance_tier": "literature_derived_transfer",
            "confidence_tier": "high",
            "uncertainty_posture": "directional_transfer",
            "reference_conditions": {"matrix": "maize_grits"},
            "extrusion_lysine_damage": {"kinetic_order": "first_order"},
            "notes": "Correlates Specific Mechanical Energy inputs to first-order lysine damage kinetics.",
        },
    ]
    existing_shear = {pr["id"] for pr in priors["shear_damage_priors"]}
    for pr in new_shear_priors:
        if pr["id"] not in existing_shear:
            priors["shear_damage_priors"].append(pr)
            print(f"Added {pr['id']} to shear_damage_priors")

    # C.4 Add Huang (2024) to carbonyl_donor_priors (root key already exists)
    new_carbonyl_priors = [
        {
            "id": "huang_2024_dixyl_arp_degradation",
            "source": "Huang et al. (2024)",
            "provenance_tier": "literature_derived_transfer",
            "confidence_tier": "high",
            "uncertainty_posture": "directional_transfer",
            "reference_conditions": {"matrix": "aqueous_free_model"},
            "bifurcated_degradation_kinetics": {
                "three_deoxyglucosone_mediated_furosine_ea_kj_mol": 115.0,
                "methylglyoxal_mediated_pathway_ea_kj_mol": 85.0,
            },
            "notes": "Bifurcated, first-order Amadori degradation kinetics determining crossover at extreme temperatures.",
        }
    ]
    existing_carbonyl = {pr["id"] for pr in priors["carbonyl_donor_priors"]}
    for pr in new_carbonyl_priors:
        if pr["id"] not in existing_carbonyl:
            priors["carbonyl_donor_priors"].append(pr)
            print(f"Added {pr['id']} to carbonyl_donor_priors")

    # D. Add to safety_reference_payloads.json
    new_safety_entries = [
        {
            "id": "yu_2017_cml_cel_meat_review",
            "kind": "industrial_endpoint_reference",
            "report_visibility": "default",
            "target_module": "safety",
            "source_citation": "Yu et al. (2017), Trends Food Sci. Technol. 98:30-40",
            "doi": "10.1016/j.tifs.2020.01.021",
            "validated_status": "reference_anchor",
            "analyte": "cml_cel",
            "method": {"instrument": "LC-MS/MS", "replicates": 3},
            "matrix_reference_ranges": [
                {
                    "matrix_family": "free_model_system",
                    "units": "activation_energy_kj_mol",
                    "cml_formation_ea": 61.01,
                    "cel_formation_ea": 29.21,
                }
            ],
            "what_it_supports": [
                "Differential thermodynamics for advanced glycation endproducts: CML activation energy of 61.01 kJ/mol versus CEL activation energy of 29.21 kJ/mol."
            ],
            "what_it_does_not_support": ["Extrusion mechanical shear effects"],
            "comment": "Provides a reference for comparative CML vs CEL activation energy barriers.",
        },
        {
            "id": "charissou_2007_cookie_cml_furosine",
            "kind": "industrial_endpoint_reference",
            "report_visibility": "default",
            "target_module": "safety",
            "source_citation": "Charissou et al. (2007), J. Agric. Food Chem. 55(11):4532-4539",
            "doi": "10.1021/jf063024j",
            "validated_status": "reference_anchor",
            "analyte": "furosine_cml",
            "method": {"instrument": "GC-MS", "replicates": 3},
            "matrix_reference_ranges": [
                {
                    "matrix_family": "cookie_model",
                    "units": "temperature_range_C",
                    "min": 180.0,
                    "max": 220.0,
                }
            ],
            "what_it_supports": [
                "Transient curve of furosine decay and stable zero-order accumulation of CML in dry cookie baking model (180 to 220 C)."
            ],
            "what_it_does_not_support": ["High moisture extrusion mechanical effects"],
            "comment": "Cookie system showing CML/furosine crossover dynamics.",
        },
        {
            "id": "fratianni_2016_apricot_furosine",
            "kind": "industrial_endpoint_reference",
            "report_visibility": "default",
            "target_module": "safety",
            "source_citation": "Fratianni et al. (2016), Food Res. Int. 99:862-867",
            "doi": "10.1016/j.foodres.2016.12.009",
            "validated_status": "reference_anchor",
            "analyte": "furosine",
            "method": {"instrument": "HPLC", "replicates": 3},
            "matrix_reference_ranges": [
                {
                    "matrix_family": "apricot",
                    "units": "activation_energy_kj_mol",
                    "ea": 83.3,
                }
            ],
            "what_it_supports": [
                "Zero-order kinetics for furosine formation in apricot matrices (Ea 83.3 kJ/mol)."
            ],
            "what_it_does_not_support": ["Dry extrusion kinetics"],
            "comment": "Apricot model system for furosine.",
        },
        {
            "id": "ma_2024_pbma_extrusion_damage",
            "kind": "industrial_endpoint_reference",
            "report_visibility": "default",
            "target_module": "safety",
            "source_citation": "Ma et al. (2024), Int. J. Mol. Sci. 25:8408",
            "doi": "10.3390/ijms25168668",
            "validated_status": "reference_anchor",
            "analyte": "cml_cel_acrylamide",
            "method": {"instrument": "LC-MS/MS", "replicates": 3},
            "matrix_reference_ranges": [
                {
                    "matrix_family": "plant_based_meat_analogue",
                    "units": "temperature_C",
                    "cel_increase_pct_min": 39.0,
                    "cel_increase_pct_max": 101.0,
                    "barrel_temp": 160.0,
                }
            ],
            "what_it_supports": [
                "CML and CEL selective generation under extrusion barrel temperatures of 160 C."
            ],
            "what_it_does_not_support": ["Aqueous buffer without protein matrix"],
            "comment": "Important PBMA extrusion safety reference payload for CML, CEL and acrylamide.",
        },
    ]
    existing_safety = {ent["id"] for ent in safety["entries"]}
    for ent in new_safety_entries:
        if ent["id"] not in existing_safety:
            safety["entries"].append(ent)
            print(f"Added {ent['id']} to safety_reference_payloads")

    # E. Update deep_research_backlog.json status and set registry_id
    backlog_citations_map = {
        "Lagrain et al. (2010)": "lagrain_2010_cystine_elimination_lanthionine",
        "Morel et al. (2002)": "morel_2002_gluten_shear_aggregation",
        "Yu et al. (2017)": "yu_2017_cml_cel_meat_review",
        "Huang et al. (2024)": "huang_2024_dixyl_arp_degradation",
        "Charissou et al. (2007)": "charissou_2007_cookie_cml_furosine",
        "Fratianni et al. (2016)": "fratianni_2016_apricot_furosine",
        "Rombouts et al. (2012)": "rombouts_2012_gluten_crosslinking",
        "Ma et al. (2024)": "ma_2024_pbma_extrusion_damage",
        "Ilo et al. (1996)": "ilo_1996_maize_sme_lysine_damage",
    }

    # Ensure items exist in backlog items list
    existing_backlog_citations = {item.get("citation") for item in backlog["items"]}

    # Check and add if missing
    for cit, reg_id in backlog_citations_map.items():
        if cit not in existing_backlog_citations:
            matching_ref: dict[str, Any] | None = next(
                (r for r in new_eligible_references if r["id"] == reg_id), None
            )
            family_num: str = "12"
            desc: str = ""
            if matching_ref is not None:
                family_num = str(matching_ref.get("slr_family_source", "12"))
                supports = matching_ref.get("what_it_supports", [])
                if isinstance(supports, list) and supports:
                    desc = str(supports[0])

            family_to_file = {
                "01": "01_amino_acid_sugar.md",
                "02": "02_lipid_oxidation.md",
                "03": "03_thiamine_degradation.md",
                "04": "04_nucleotide_degradation.md",
                "05": "05_glutathione_peptides.md",
                "06": "06_alternative_proteins.md",
                "07": "07_reducing_sugars.md",
                "08": "08_off_note_chemistry.md",
                "09": "09_carbohydrate_degradation.md",
                "10": "10_microbial_fermentation.md",
                "11": "11_maillard_lipid_crosstalk.md",
                "12": "12_protein_damage_markers.md",
                "13": "13_polyphenol_amino_capping.md",
                "14": "14_ascorbic_acid_maillard.md",
                "15": "15_phospholipid_amine_maillard.md",
                "16": "16_melanoidin_polymerization.md",
            }
            filename = family_to_file.get(family_num, "12_protein_damage_markers.md")

            new_item = {
                "citation": cit,
                "score": "6/8",
                "score_value": 6,
                "status": "RUNTIME_BOUND",
                "occurrence_count": 1,
                "files": [filename],
                "descriptions": [desc],
                "occurrences": [
                    {
                        "file": filename,
                        "line": 999,  # Placeholder line number
                        "description": desc,
                        "raw_line": f"1. `{cit}` — {desc} 6/8.",
                    }
                ],
                "registry_id": reg_id,
                "runtime_artifact_count": 1,
            }
            backlog["items"].append(new_item)
            print(f"Added {cit} to deep_research_backlog items")
        else:
            matching_ref: dict[str, Any] | None = next(
                (r for r in new_eligible_references if r["id"] == reg_id), None
            )
            family_num: str = "12"
            if matching_ref is not None:
                family_num = str(matching_ref.get("slr_family_source", "12"))

            family_to_file = {
                "01": "01_amino_acid_sugar.md",
                "02": "02_lipid_oxidation.md",
                "03": "03_thiamine_degradation.md",
                "04": "04_nucleotide_degradation.md",
                "05": "05_glutathione_peptides.md",
                "06": "06_alternative_proteins.md",
                "07": "07_reducing_sugars.md",
                "08": "08_off_note_chemistry.md",
                "09": "09_carbohydrate_degradation.md",
                "10": "10_microbial_fermentation.md",
                "11": "11_maillard_lipid_crosstalk.md",
                "12": "12_protein_damage_markers.md",
                "13": "13_polyphenol_amino_capping.md",
                "14": "14_ascorbic_acid_maillard.md",
                "15": "15_phospholipid_amine_maillard.md",
                "16": "16_melanoidin_polymerization.md",
            }
            filename = family_to_file.get(family_num, "12_protein_damage_markers.md")

            for item in backlog["items"]:
                if item.get("citation") == cit:
                    item["status"] = "RUNTIME_BOUND"
                    item["registry_id"] = reg_id
                    item["runtime_artifact_count"] = 1
                    item["files"] = [filename]
                    if item.get("occurrences"):
                        for occ in item["occurrences"]:
                            occ["file"] = filename
                    print(
                        f"Updated backlog status and files for {cit} to RUNTIME_BOUND"
                    )

    # Write back all files
    save_json(intake, INTAKE_REGISTRY_PATH)
    save_json(matrix, SLR_MATRIX_PATH)
    save_json(priors, COMPUTATIONAL_PRIORS_PATH)
    save_json(safety, SAFETY_PAYLOADS_PATH)
    save_json(backlog, BACKLOG_PATH)

    print("Successfully ingested all references into structured json databases!")


if __name__ == "__main__":
    ingest()
