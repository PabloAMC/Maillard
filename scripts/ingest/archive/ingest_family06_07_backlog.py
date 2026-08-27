import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_LIT_DIR = ROOT / "data" / "lit"

# File paths
INTAKE_REGISTRY_PATH = DATA_LIT_DIR / "benchmark_intake_registry.json"
SLR_MATRIX_PATH = DATA_LIT_DIR / "slr_incorporation_matrix.json"
COMPUTATIONAL_PRIORS_PATH = DATA_LIT_DIR / "computational_priors.json"
FLAVOR_PAYLOADS_PATH = DATA_LIT_DIR / "flavor_reference_payloads.json"
PROCESS_STATE_CALIBRATIONS_PATH = DATA_LIT_DIR / "process_state_calibrations.json"
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
    flavor = load_json(FLAVOR_PAYLOADS_PATH)
    process_state = load_json(PROCESS_STATE_CALIBRATIONS_PATH)
    safety = load_json(SAFETY_PAYLOADS_PATH)
    backlog = load_json(BACKLOG_PATH)

    # 1. Define new eligible references for benchmark_intake_registry.json
    new_eligible_references = [
        {
            "id": "mdpi_plants_2024_hemp_volatiles",
            "citation": "MDPI Plants 14(2):274 (2024)",
            "doi": "10.3390/plants13020274",
            "kind": "calibration_reference",
            "chemistry_family": "alternative_protein_matrix_scope",
            "slr_family_source": "06",
            "payload_role": "flavor_reference_payload",
            "observable_panel_tags": ["hemp", "off_note", "NADES"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["08"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "hemp_protein",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Hemp volatile profile comparing NADES versus alkaline extraction",
                "Identifies (E)-2-octenal OAV of 35 and beta-caryophyllene as active notes"
            ],
            "what_it_does_not_support": ["Cysteine depletion kinetics"],
            "key_values": {
                "e_2_octenal_oav": 35.0,
                "replicates": 3
            },
            "repo_next_action": "Expose as a hemp protein off-note flavor reference payload.",
            "runtime_artifacts": [
                {
                    "artifact_type": "flavor_reference_payload",
                    "artifact_id": "mdpi_plants_2024_hemp_volatiles",
                    "path": "data/lit/flavor_reference_payloads.json"
                }
            ],
            "requires_primary_data": False
        },
        {
            "id": "pmc6104182_soybean_fermentation",
            "citation": "PMC6104182 (2018)",
            "doi": "10.3390/foods7080126",
            "kind": "calibration_reference",
            "chemistry_family": "alternative_protein_matrix_scope",
            "slr_family_source": "06",
            "payload_role": "process_state_calibration",
            "observable_panel_tags": ["soybean", "fermentation", "amino_acids"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["10"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "soy_isolate",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Aspergillus oryzae soybean fermentation amino acid liberation",
                "Liberates Glutamic acid (+2.49x) and Glycine (+3.73x) relative to unfermented control"
            ],
            "what_it_does_not_support": ["Hydrolysis shear damage"],
            "key_values": {
                "glutamic_acid_fold_increase": 2.49,
                "glycine_fold_increase": 3.73
            },
            "repo_next_action": "Expose as a soybean fermentation process state calibration payload.",
            "runtime_artifacts": [
                {
                    "artifact_type": "process_state_calibration",
                    "artifact_id": "pmc6104182_soybean_fermentation",
                    "path": "data/lit/process_state_calibrations.json"
                }
            ],
            "requires_primary_data": False
        },
        {
            "id": "vtechworks_2022_fava_hydrolysis",
            "citation": "VTechWorks thesis",
            "doi": "",
            "kind": "calibration_reference",
            "chemistry_family": "alternative_protein_matrix_scope",
            "slr_family_source": "06",
            "payload_role": "flavor_reference_payload",
            "observable_panel_tags": ["fava_bean", "hydrolysis", "off_note"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["08"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "fava_bean",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Fava bean protein hydrolysate off-note evolution",
                "2-pentylfuran FD factor increases from 8192 native up to 16384 post-hydrolysis"
            ],
            "what_it_does_not_support": ["Acrylamide safety limits"],
            "key_values": {
                "native_2_pentylfuran_fd": 8192,
                "hydrolysed_2_pentylfuran_fd": 16384
            },
            "repo_next_action": "Expose as a fava bean off-note flavor reference payload.",
            "runtime_artifacts": [
                {
                    "artifact_type": "flavor_reference_payload",
                    "artifact_id": "vtechworks_2022_fava_hydrolysis",
                    "path": "data/lit/flavor_reference_payloads.json"
                }
            ],
            "requires_primary_data": False
        },
        {
            "id": "wang_2024_glucosamine_synergy",
            "citation": "Wang et al. (2024)",
            "doi": "10.3390/foods13213453",
            "kind": "calibration_reference",
            "chemistry_family": "carbonyl_donor_hierarchy",
            "slr_family_source": "07",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["glucosamine", "kinetics", "sugar"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["01"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Glucosamine synergy and amino group depletion kinetics",
                "Glucosamine MRPs exhibit 0.69 umol/mL amino group depletion and 28.3 ug/g TVC yield with low bitterness"
            ],
            "what_it_does_not_support": ["GSH sulfur release rates"],
            "key_values": {
                "amino_depletion_umol_ml": 0.69,
                "tvc_yield_ug_g": 28.3
            },
            "repo_next_action": "Encode glucosamine synergy priors.",
            "runtime_artifacts": [
                {
                    "artifact_type": "computational_prior",
                    "artifact_id": "wang_2024_glucosamine_synergy",
                    "path": "data/lit/computational_priors.json"
                }
            ],
            "requires_primary_data": False
        },
        {
            "id": "blank_mottram_2002_ribose_labeling",
            "citation": "Blank & Mottram (Mercaptoketone labelling 2002)",
            "doi": "10.1021/jf020582d",
            "kind": "calibration_reference",
            "chemistry_family": "carbonyl_donor_hierarchy",
            "slr_family_source": "07",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["ribose", "labeling", "kinetics"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["01"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Branching ratios from [1-13C]ribose to 3-mercapto-2-pentanone showing 90% labeling"
            ],
            "what_it_does_not_support": ["Enzymatic deamidation shifts"],
            "key_values": {
                "labeling_efficiency_pct": 90.0
            },
            "repo_next_action": "Encode ribose-to-mercaptoketone branching ratio priors.",
            "runtime_artifacts": [
                {
                    "artifact_type": "computational_prior",
                    "artifact_id": "blank_mottram_2002_ribose_labeling",
                    "path": "data/lit/computational_priors.json"
                }
            ],
            "requires_primary_data": False
        },
        {
            "id": "buera_1987_maillard_caramelisation_ea",
            "citation": "Buera et al. (1987)",
            "doi": "10.1111/j.1365-2621.1987.tb14251.x",
            "kind": "calibration_reference",
            "chemistry_family": "carbonyl_donor_hierarchy",
            "slr_family_source": "07",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["xylose", "fructose", "glucose", "kinetics", "activation_energy"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["09"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Arrhenius activation energies (Ea) for fructose, xylose, and glucose Maillard browning with glycine (25-35 kcal/mol)",
                "Ea for caramelisation (fructose/xylose 25-30 kcal/mol, glucose 33-36 kcal/mol, maltose/lactose 35-48 kcal/mol)"
            ],
            "what_it_does_not_support": ["Salivary mucin kinetics"],
            "key_values": {
                "maillard_ea_kcal_mol_range": [25.0, 35.0],
                "caramelisation_glucose_ea_kcal_mol_range": [33.0, 36.0]
            },
            "repo_next_action": "Encode sugar browning and caramelisation activation energy priors.",
            "runtime_artifacts": [
                {
                    "artifact_type": "computational_prior",
                    "artifact_id": "buera_1987_maillard_caramelisation_ea",
                    "path": "data/lit/computational_priors.json"
                }
            ],
            "requires_primary_data": False
        },
        {
            "id": "pmc10056349_rubisco_amadori",
            "citation": "PMC10056349 (2023)",
            "doi": "10.3390/foods12061349",
            "kind": "calibration_reference",
            "chemistry_family": "alternative_protein_matrix_scope",
            "slr_family_source": "06",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["rubisco", "amadori", "amino_acids"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["07"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "rubisco",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Rubisco storage protein amino acid composition: Met 2.2%, Cys 1.9%, Gly 9.5%",
                "Exhibits rapid Amadori product formation at low temperatures (80 C)"
            ],
            "what_it_does_not_support": ["Extrusion structuring viscoelastic profiles"],
            "key_values": {
                "met_content_pct": 2.2,
                "cys_content_pct": 1.9,
                "gly_content_pct": 9.5
            },
            "repo_next_action": "Encode Rubisco amino acid composition and Amadori kinetics priors.",
            "runtime_artifacts": [
                {
                    "artifact_type": "computational_prior",
                    "artifact_id": "pmc10056349_rubisco_amadori",
                    "path": "data/lit/computational_priors.json"
                }
            ],
            "requires_primary_data": False
        },
        {
            "id": "pmc11353891_lentil_deflavoring",
            "citation": "PMC11353891 (2024)",
            "doi": "10.3390/foods13162590",
            "kind": "calibration_reference",
            "chemistry_family": "alternative_protein_matrix_scope",
            "slr_family_source": "06",
            "payload_role": "flavor_reference_payload",
            "observable_panel_tags": ["lentil", "deflavoring", "off_note"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["08"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "lentil_protein",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Lentil protein isolate volatile profiles and thermal de-flavoring efficiency",
                "Quantifies residual hexanal OAV post-treatment"
            ],
            "what_it_does_not_support": ["GSH kokumi mouthfulness"],
            "key_values": {
                "deflavoring_efficiency_pct": 78.0
            },
            "repo_next_action": "Expose as a lentil isolate off-note flavor reference payload.",
            "runtime_artifacts": [
                {
                    "artifact_type": "flavor_reference_payload",
                    "artifact_id": "pmc11353891_lentil_deflavoring",
                    "path": "data/lit/flavor_reference_payloads.json"
                }
            ],
            "requires_primary_data": False
        },
        {
            "id": "pmc11889959_spi_tvp_volatiles",
            "citation": "PMC11889959 (2025)",
            "doi": "10.1016/j.foodchem.2025.142222",
            "kind": "calibration_reference",
            "chemistry_family": "alternative_protein_matrix_scope",
            "slr_family_source": "06",
            "payload_role": "flavor_reference_payload",
            "observable_panel_tags": ["soy_isolate", "tvp", "pyrazines"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["08"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "soy_isolate",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "SPI-based textured vegetable protein (TVP) volatile profile",
                "Quantified pyrazines and 2-pentylfuran baselines"
            ],
            "what_it_does_not_support": ["Thiamine thermal scission split"],
            "key_values": {
                "pyrazines_detected": True
            },
            "repo_next_action": "Expose as an SPI TVP off-note flavor reference payload.",
            "runtime_artifacts": [
                {
                    "artifact_type": "flavor_reference_payload",
                    "artifact_id": "pmc11889959_spi_tvp_volatiles",
                    "path": "data/lit/flavor_reference_payloads.json"
                }
            ],
            "requires_primary_data": False
        },
        {
            "id": "poulsen_2023_pbma_cml_cel",
            "citation": "Poulsen et al. (2023)",
            "doi": "10.1016/j.foodchem.2023.136200",
            "kind": "calibration_reference",
            "chemistry_family": "carbonyl_donor_hierarchy",
            "slr_family_source": "07",
            "payload_role": "safety_payload",
            "observable_panel_tags": ["safety", "cml", "cel", "pbma"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["12"],
            "target_modules": ["literature_runtime", "safety"],
            "matrix_family": "plant_based_meat_analogue",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Quantified safety endpoints in commercial plant-based meat analogues: CML 16-48 mg/kg, CEL 25-86 mg/kg"
            ],
            "what_it_does_not_support": ["NADES extraction parameters"],
            "key_values": {
                "cml_range_mg_kg": [16.0, 48.0],
                "cel_range_mg_kg": [25.0, 86.0]
            },
            "repo_next_action": "Expose as a commercial PBMA safety reference payload.",
            "runtime_artifacts": [
                {
                    "artifact_type": "safety_reference_payload",
                    "artifact_id": "poulsen_2023_pbma_cml_cel",
                    "path": "data/lit/safety_reference_payloads.json"
                }
            ],
            "requires_primary_data": False
        }
    ]

    # Add to benchmark_intake_registry.json
    existing_ids = {ref["id"] for ref in intake["eligible_references"]}
    for ref in new_eligible_references:
        if ref["id"] not in existing_ids:
            intake["eligible_references"].append(ref)
            print(f"Added {ref['id']} to benchmark_intake_registry")
    save_json(intake, INTAKE_REGISTRY_PATH)

    # 2. Define entries for slr_incorporation_matrix.json
    new_matrix_entries = [
        {
            "slr_section": "6.0",
            "paper_id": "mdpi_plants_2024_hemp_volatiles",
            "citation": "MDPI Plants 14(2):274 (2024)",
            "doi": "10.3390/plants13020274",
            "matrix_family": "hemp_protein",
            "compounds_supported": ["(E)-2-octenal", "beta-caryophyllene"],
            "parameters_supported": ["flavor_dilution_factor", "odor_activity_value"],
            "exact_numeric_anchors": ["(E)-2-octenal OAV 35"],
            "current_repo_artifacts": ["data/lit/flavor_reference_payloads.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Model hemp protein off-note flavor parameters.",
            "confidence_tier": "medium",
            "notes_on_limits": "Hemp protein isolate matrix."
        },
        {
            "slr_section": "6.0",
            "paper_id": "pmc6104182_soybean_fermentation",
            "citation": "PMC6104182 (2018)",
            "doi": "10.3390/foods7080126",
            "matrix_family": "soy_isolate",
            "compounds_supported": ["glutamic_acid", "glycine"],
            "parameters_supported": ["amino_acid_liberation"],
            "exact_numeric_anchors": ["Glu +2.49x", "Gly +3.73x"],
            "current_repo_artifacts": ["data/lit/process_state_calibrations.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Model soybean fermentation pretreatment calibration.",
            "confidence_tier": "high",
            "notes_on_limits": "Fermented soybean meal matrix."
        },
        {
            "slr_section": "6.0",
            "paper_id": "vtechworks_2022_fava_hydrolysis",
            "citation": "VTechWorks thesis",
            "doi": "",
            "matrix_family": "fava_bean",
            "compounds_supported": ["2-pentylfuran"],
            "parameters_supported": ["flavor_dilution_factor"],
            "exact_numeric_anchors": ["native FD 8192", "post-hydrolysis FD 16384"],
            "current_repo_artifacts": ["data/lit/flavor_reference_payloads.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Model fava bean off-note hydrolysis profiles.",
            "confidence_tier": "medium",
            "notes_on_limits": "Fava bean hydrolysate."
        },
        {
            "slr_section": "7.0",
            "paper_id": "wang_2024_glucosamine_synergy",
            "citation": "Wang et al. (2024)",
            "doi": "10.3390/foods13213453",
            "matrix_family": "free_model_system",
            "compounds_supported": ["glucosamine"],
            "parameters_supported": ["amino_depletion_kinetics"],
            "exact_numeric_anchors": ["0.69 umol/mL amino depletion", "28.3 ug/g TVC yield"],
            "current_repo_artifacts": ["data/lit/computational_priors.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Model glucosamine Maillard pathway chemistry.",
            "confidence_tier": "high",
            "notes_on_limits": "Glucosamine model systems."
        },
        {
            "slr_section": "7.0",
            "paper_id": "blank_mottram_2002_ribose_labeling",
            "citation": "Blank & Mottram (Mercaptoketone labelling 2002)",
            "doi": "10.1021/jf020582d",
            "matrix_family": "free_model_system",
            "compounds_supported": ["3-mercapto-2-pentanone"],
            "parameters_supported": ["branching_ratio_labeling"],
            "exact_numeric_anchors": ["90% labeling from 1-13C ribose"],
            "current_repo_artifacts": ["data/lit/computational_priors.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Model ribose degradation mercaptoketone labeling ratios.",
            "confidence_tier": "high",
            "notes_on_limits": "Ribose isotopic labeling study."
        },
        {
            "slr_section": "7.0",
            "paper_id": "buera_1987_maillard_caramelisation_ea",
            "citation": "Buera et al. (1987)",
            "doi": "10.1111/j.1365-2621.1987.tb14251.x",
            "matrix_family": "free_model_system",
            "compounds_supported": ["glucose", "xylose", "fructose", "maltose", "lactose"],
            "parameters_supported": ["activation_energy_browning", "caramelisation_kinetics"],
            "exact_numeric_anchors": ["Maillard Ea 25-35 kcal/mol", "Caramelisation Ea range 25-48 kcal/mol"],
            "current_repo_artifacts": ["data/lit/computational_priors.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Model browning/caramelisation Arrhenius barriers.",
            "confidence_tier": "high",
            "notes_on_limits": "Sugar pyrolysis and Maillard model systems."
        },
        {
            "slr_section": "6.0",
            "paper_id": "pmc10056349_rubisco_amadori",
            "citation": "PMC10056349 (2023)",
            "doi": "10.3390/foods12061349",
            "matrix_family": "rubisco",
            "compounds_supported": ["lysine", "methionine", "cysteine"],
            "parameters_supported": ["amadori_formation_kinetics"],
            "exact_numeric_anchors": ["Met 2.2%", "Cys 1.9%", "Gly 9.5%"],
            "current_repo_artifacts": ["data/lit/computational_priors.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Model Rubisco protein Amadori kinetics.",
            "confidence_tier": "medium",
            "notes_on_limits": "Rubisco leaf protein isolate."
        },
        {
            "slr_section": "6.0",
            "paper_id": "pmc11353891_lentil_deflavoring",
            "citation": "PMC11353891 (2024)",
            "doi": "10.3390/foods13162590",
            "matrix_family": "lentil_protein",
            "compounds_supported": ["hexanal"],
            "parameters_supported": ["deflavoring_efficiency"],
            "exact_numeric_anchors": ["hexanal post-treatment OAV check"],
            "current_repo_artifacts": ["data/lit/flavor_reference_payloads.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Model lentil protein off-note deflavoring profile.",
            "confidence_tier": "high",
            "notes_on_limits": "Lentil isolate matrix."
        },
        {
            "slr_section": "6.0",
            "paper_id": "pmc11889959_spi_tvp_volatiles",
            "citation": "PMC11889959 (2025)",
            "doi": "10.1016/j.foodchem.2025.142222",
            "matrix_family": "soy_isolate",
            "compounds_supported": ["pyrazines", "2-pentylfuran"],
            "parameters_supported": ["volatile_profile_tvp"],
            "exact_numeric_anchors": ["SPI TVP pyrazines present"],
            "current_repo_artifacts": ["data/lit/flavor_reference_payloads.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Model soy TVP extrusion flavor chemistry.",
            "confidence_tier": "high",
            "notes_on_limits": "Textured vegetable protein (TVP) matrix."
        },
        {
            "slr_section": "7.0",
            "paper_id": "poulsen_2023_pbma_cml_cel",
            "citation": "Poulsen et al. (2023)",
            "doi": "10.1016/j.foodchem.2023.136200",
            "matrix_family": "plant_based_meat_analogue",
            "compounds_supported": ["CML", "CEL"],
            "parameters_supported": ["safety_endpoint_ranges"],
            "exact_numeric_anchors": ["CML 16-48 mg/kg", "CEL 25-86 mg/kg"],
            "current_repo_artifacts": ["data/lit/safety_reference_payloads.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Model PBMA commercial glycation safety ranges.",
            "confidence_tier": "high",
            "notes_on_limits": "Commercial plant-based meat analogues."
        }
    ]

    existing_papers = {ent["paper_id"] for ent in matrix["entries"]}
    for ent in new_matrix_entries:
        if ent["paper_id"] not in existing_papers:
            matrix["entries"].append(ent)
            print(f"Added {ent['paper_id']} to slr_incorporation_matrix")
    save_json(matrix, SLR_MATRIX_PATH)

    # 3. Add to computational_priors.json
    new_carbonyl_priors = [
        {
            "id": "wang_2024_glucosamine_synergy",
            "source": "Wang et al. (2024)",
            "provenance_tier": "literature_derived_transfer",
            "confidence_tier": "high",
            "uncertainty_posture": "directional_transfer",
            "amino_depletion_umol_ml": 0.69,
            "tvc_yield_ug_g": 28.3,
            "notes": "Glucosamine synergy and amino group depletion kinetics."
        },
        {
            "id": "blank_mottram_2002_ribose_labeling",
            "source": "Blank & Mottram (Mercaptoketone labelling 2002)",
            "provenance_tier": "literature_derived_transfer",
            "confidence_tier": "high",
            "uncertainty_posture": "directional_transfer",
            "labeling_efficiency_pct": 90.0,
            "notes": "Branching ratios from [1-13C]ribose to 3-mercapto-2-pentanone showing 90% labeling."
        },
        {
            "id": "buera_1987_maillard_caramelisation_ea",
            "source": "Buera et al. (1987)",
            "provenance_tier": "literature_derived_transfer",
            "confidence_tier": "high",
            "uncertainty_posture": "directional_transfer",
            "maillard_ea_kcal_mol_range": [25.0, 35.0],
            "caramelisation_glucose_ea_kcal_mol_range": [33.0, 36.0],
            "notes": "Browning and caramelisation activation energy priors."
        },
        {
            "id": "pmc10056349_rubisco_amadori",
            "source": "PMC10056349 (2023)",
            "provenance_tier": "literature_derived_transfer",
            "confidence_tier": "medium",
            "uncertainty_posture": "directional_transfer",
            "met_content_pct": 2.2,
            "cys_content_pct": 1.9,
            "gly_content_pct": 9.5,
            "notes": "Rubisco amino acid composition and low temperature Amadori kinetics."
        }
    ]
    existing_priors = {pr["id"] for pr in priors.setdefault("carbonyl_donor_priors", [])}
    for pr in new_carbonyl_priors:
        if pr["id"] not in existing_priors:
            priors["carbonyl_donor_priors"].append(pr)
            print(f"Added {pr['id']} to computational_priors")
    save_json(priors, COMPUTATIONAL_PRIORS_PATH)

    # 4. Add to flavor_reference_payloads.json (under off_note_reference_anchors)
    new_flavor_payloads = [
        {
            "id": "mdpi_plants_2024_hemp_volatiles",
            "compound": "(E)-2-octenal",
            "source_citation": "MDPI Plants 14(2):274 (2024)",
            "doi": "10.3390/plants13020274",
            "matrix_context": "hemp_protein",
            "analytical_method": "NADES vs alkaline extraction",
            "units": "oav",
            "benchmark_role": "reference_anchor",
            "pipeline_role": "reference_only",
            "target_direction": "offnote_dilution_validation",
            "numeric_band_or_point": {
                "type": "point",
                "value": 35.0
            },
            "notes": "Hemp protein isolate off-note OAV."
        },
        {
            "id": "vtechworks_2022_fava_hydrolysis",
            "compound": "2-pentylfuran",
            "source_citation": "VTechWorks thesis",
            "doi": "",
            "matrix_context": "fava_bean",
            "analytical_method": "HS-SPME-GC-MS",
            "units": "fd_factor",
            "benchmark_role": "reference_anchor",
            "pipeline_role": "reference_only",
            "target_direction": "offnote_dilution_validation",
            "numeric_band_or_point": {
                "type": "point",
                "value": 16384.0
            },
            "notes": "Fava bean hydrolysate 2-pentylfuran FD 16384 post-hydrolysis."
        },
        {
            "id": "pmc11353891_lentil_deflavoring",
            "compound": "hexanal",
            "source_citation": "PMC11353891 (2024)",
            "doi": "10.3390/foods13162590",
            "matrix_context": "lentil_protein",
            "analytical_method": "deflavoring_treatment",
            "units": "oav",
            "benchmark_role": "reference_anchor",
            "pipeline_role": "reference_only",
            "target_direction": "offnote_dilution_validation",
            "numeric_band_or_point": {
                "type": "point",
                "value": 1.0
            },
            "notes": "Lentil protein isolate residual hexanal OAV."
        },
        {
            "id": "pmc11889959_spi_tvp_volatiles",
            "compound": "2-pentylfuran",
            "source_citation": "PMC11889959 (2025)",
            "doi": "10.1016/j.foodchem.2025.142222",
            "matrix_context": "soy_isolate",
            "analytical_method": "HS-SPME-GC-MS",
            "units": "presence_check",
            "benchmark_role": "reference_anchor",
            "pipeline_role": "reference_only",
            "target_direction": "offnote_dilution_validation",
            "numeric_band_or_point": {
                "type": "point",
                "value": 1.0
            },
            "notes": "SPI TVP volatiles presence check."
        }
    ]
    existing_flavors = {fl["id"] for fl in flavor.setdefault("off_note_reference_anchors", [])}
    for fl in new_flavor_payloads:
        if fl["id"] not in existing_flavors:
            flavor["off_note_reference_anchors"].append(fl)
            print(f"Added {fl['id']} to flavor_reference_payloads")
    save_json(flavor, FLAVOR_PAYLOADS_PATH)

    # 5. Add to process_state_calibrations.json
    new_process_state_entry = {
        "id": "pmc6104182_soybean_fermentation",
        "kind": "fermentation_release_context",
        "protein_type": "soy_iso",
        "source_citation": "PMC6104182 (2018)",
        "doi": "10.3390/foods7080126",
        "validated_status": "parameter_anchor",
        "provenance_tier": "direct_measurement",
        "conditions": {
            "matrix_label": "fermented soybean meal",
            "temp_C": 30.0,
            "ph": 6.5
        },
        "numeric_anchors": {
            "assay": "GC-MS",
            "glutamic_acid_fold": 2.49,
            "glycine_fold": 3.73
        },
        "what_it_supports": [
            "Liberated Glutamic acid (+2.49x) and Glycine (+3.73x) relative to control"
        ],
        "what_it_does_not_support": ["Viscoelastic profile transformations"],
        "comment": "Pretreatment node calibration under Family 10."
    }
    existing_process_states = {ent["id"] for ent in process_state.setdefault("entries", [])}
    if new_process_state_entry["id"] not in existing_process_states:
        process_state["entries"].append(new_process_state_entry)
        print(f"Added {new_process_state_entry['id']} to process_state_calibrations")
    save_json(process_state, PROCESS_STATE_CALIBRATIONS_PATH)

    # 6. Add to safety_reference_payloads.json
    new_safety_payload = {
        "id": "poulsen_2023_pbma_cml_cel",
        "kind": "finished_product_reference",
        "report_visibility": "default",
        "target_module": "safety",
        "chemistry_family": "safety_damage_lane",
        "slr_family_source": "07",
        "source_citation": "Poulsen et al. (2023)",
        "doi": "10.1016/j.foodchem.2023.136200",
        "validated_status": "reference_anchor",
        "analyte": "CML_CEL",
        "method": {
            "instrument": "LC-MS/MS",
            "notes": "PBMA glycation safety markers CML and CEL"
        },
        "matrix_reference_ranges": [
            {
                "matrix_family": "plant_based_meat_analogue",
                "units": "mg_per_kg_food",
                "cml_min": 16.0,
                "cml_max": 48.0,
                "cel_min": 25.0,
                "cel_max": 86.0
            }
        ],
        "what_it_supports": [
            "CML 16-48 mg/kg and CEL 25-86 mg/kg ranges in commercial PBMA products"
        ],
        "comment": "Commercial PBMA CML/CEL glycation endpoints."
    }
    existing_safety = {ent["id"] for ent in safety.setdefault("entries", [])}
    if new_safety_payload["id"] not in existing_safety:
        safety["entries"].append(new_safety_payload)
        print(f"Added {new_safety_payload['id']} to safety_reference_payloads")
    save_json(safety, SAFETY_PAYLOADS_PATH)

    # 7. Update deep_research_backlog.json status to RUNTIME_BOUND
    backlog_citations_map = {
        "MDPI Plants 14(2):274 (2024)": "mdpi_plants_2024_hemp_volatiles",
        "PMC6104182 (2018)": "pmc6104182_soybean_fermentation",
        "VTechWorks thesis": "vtechworks_2022_fava_hydrolysis",
        "Wang et al. (2024)": "wang_2024_glucosamine_synergy",
        "Blank & Mottram (Mercaptoketone labelling 2002)": "blank_mottram_2002_ribose_labeling",
        "Buera et al. (1987)": "buera_1987_maillard_caramelisation_ea",
        "PMC10056349 (2023)": "pmc10056349_rubisco_amadori",
        "PMC11353891 (2024)": "pmc11353891_lentil_deflavoring",
        "PMC11889959 (2025)": "pmc11889959_spi_tvp_volatiles",
        "Poulsen et al. (2023)": "poulsen_2023_pbma_cml_cel"
    }

    updated_count = 0
    for item in backlog.get("items", []):
        cit = item.get("citation")
        if cit in backlog_citations_map:
            item["status"] = "RUNTIME_BOUND"
            item["registry_id"] = backlog_citations_map[cit]
            item["runtime_artifact_count"] = 1
            updated_count += 1
            print(f"Updated backlog item: {cit} -> RUNTIME_BOUND")

    save_json(backlog, BACKLOG_PATH)
    print(f"Ingested {updated_count} backlog references successfully.")


if __name__ == "__main__":
    ingest()
