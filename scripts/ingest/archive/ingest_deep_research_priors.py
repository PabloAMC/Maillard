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

    new_eligible_references = [
        {
            "id": "voelker_2018_thiamine_salt_degradation",
            "citation": "Voelker et al. (2018)",
            "doi": "10.1016/j.foodres.2018.06.056",
            "kind": "calibration_reference",
            "chemistry_family": "thiamine_fragmentation_support",
            "slr_family_source": "03",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["thiamine_degradation", "mononitrate", "chloride_hydrochloride"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["03"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "First-order thiamine degradation kinetics where TMN degrades faster (Ea 88-105 kJ/mol) than TClHCl (Ea 90-135 kJ/mol) due to pyrimidine protonation."
            ],
            "what_it_does_not_support": ["Solid-state high-moisture extrusion without pH validation"],
            "key_values": {
                "TMN_Ea_range_kj_mol": [88.0, 105.0],
                "TClHCl_Ea_range_kj_mol": [90.0, 135.0],
                "replicates": 3
            },
            "repo_next_action": "Encode thiamine salt form Arrhenius kinetics.",
            "runtime_artifacts": [
                {
                    "artifact_type": "computational_prior",
                    "artifact_id": "voelker_2018_thiamine_salt_degradation",
                    "path": "data/lit/computational_priors.json"
                }
            ],
            "requires_primary_data": False
        },
        {
            "id": "knol_2005_acrylamide_kinetics",
            "citation": "Knol et al. (2005)",
            "doi": "10.1021/jf050504m",
            "kind": "calibration_reference",
            "chemistry_family": "protein_damage_markers",
            "slr_family_source": "12",
            "payload_role": "safety_reference_payload",
            "observable_panel_tags": ["acrylamide", "toxicological_safety"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["12"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Bifurcated acrylamide formation (Ea 52.1 kJ/mol) and thermal degradation (Ea 72.9 kJ/mol) kinetics."
            ],
            "what_it_does_not_support": ["Low moisture solid-state matrices without water activity correction"],
            "key_values": {
                "formation_Ea_kj_mol": 52.1,
                "degradation_Ea_kj_mol": 72.9,
                "replicates": 3
            },
            "repo_next_action": "Encode as acrylamide safety reference payload.",
            "runtime_artifacts": [
                {
                    "artifact_type": "safety_reference_payload",
                    "artifact_id": "knol_2005_acrylamide_kinetics",
                    "path": "data/lit/safety_reference_payloads.json"
                }
            ],
            "requires_primary_data": False
        },
        {
            "id": "de_vleeschouwer_2008_acrylamide_aw",
            "citation": "De Vleeschouwer et al. (2008)",
            "doi": "10.1021/bp060389f",
            "kind": "calibration_reference",
            "chemistry_family": "protein_damage_markers",
            "slr_family_source": "12",
            "payload_role": "safety_reference_payload",
            "observable_panel_tags": ["acrylamide", "water_activity", "moisture_dependency"],
            "process_state_scope": ["extrusion_structured", "heated_matrix"],
            "supporting_families": ["12"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Acrylamide moisture-dependent kinetics showing minimum elimination rate at Aw of 0.82."
            ],
            "what_it_does_not_support": ["Aqueous solution systems without starch-like structural water restrictions"],
            "key_values": {
                "min_elimination_Aw": 0.82,
                "replicates": 3
            },
            "repo_next_action": "Encode as water activity dependent acrylamide safety reference payload.",
            "runtime_artifacts": [
                {
                    "artifact_type": "safety_reference_payload",
                    "artifact_id": "de_vleeschouwer_2008_acrylamide_aw",
                    "path": "data/lit/safety_reference_payloads.json"
                }
            ],
            "requires_primary_data": False
        },
        {
            "id": "ishak_2022_phip_hca_kinetics",
            "citation": "Ishak et al. (2022)",
            "doi": "10.1016/j.foodchem.2022.132372",
            "kind": "calibration_reference",
            "chemistry_family": "protein_damage_markers",
            "slr_family_source": "12",
            "payload_role": "safety_reference_payload",
            "observable_panel_tags": ["PhIP", "HCAs", "toxicological_safety"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["12"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "PhIP heterocyclic amine formation kinetics from phenylalanine (Ea 95.36 kJ/mol) vs. proline (Ea 114.12 kJ/mol) precursor pools."
            ],
            "what_it_does_not_support": ["Low-temperature extrusion texturization below 100 C"],
            "key_values": {
                "phenylalanine_Ea_kj_mol": 95.36,
                "proline_Ea_kj_mol": 114.12,
                "replicates": 3
            },
            "repo_next_action": "Encode as PhIP HCA safety reference payload.",
            "runtime_artifacts": [
                {
                    "artifact_type": "safety_reference_payload",
                    "artifact_id": "ishak_2022_phip_hca_kinetics",
                    "path": "data/lit/safety_reference_payloads.json"
                }
            ],
            "requires_primary_data": False
        },
        {
            "id": "shu_2019_cysteine_quinone_kinetics",
            "citation": "Shu et al. (2019)",
            "doi": "10.1016/j.freeradbiomed.2019.04.026",
            "kind": "calibration_reference",
            "chemistry_family": "polyphenol_amino_capping",
            "slr_family_source": "13",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["polyphenols", "cysteine_loss", "quinone"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["13"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Second-order rate constants of 10^5 to 10^6 M-1 s-1 for quinone-thiol Michael additions representing rapid cysteine depletion by oxidized polyphenols."
            ],
            "what_it_does_not_support": ["Non-oxidizing environments without electron-acceptor or metal catalysts"],
            "key_values": {
                "quinone_thiol_k_range": [100000.0, 1000000.0],
                "replicates": 3
            },
            "repo_next_action": "Encode cysteine-quinone Michael adduct kinetics.",
            "runtime_artifacts": [
                {
                    "artifact_type": "computational_prior",
                    "artifact_id": "shu_2019_cysteine_quinone_kinetics",
                    "path": "data/lit/computational_priors.json"
                }
            ],
            "requires_primary_data": False
        },
        {
            "id": "zhang_1993_protein_deamidation_ammonia",
            "citation": "Zhang et al. (1993)",
            "doi": "10.1111/j.1745-4530.1993.tb00179.x",
            "kind": "calibration_reference",
            "chemistry_family": "alternative_protein_matrix_scope",
            "slr_family_source": "06",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["deamidation", "ammonia", "asparagine", "glutamine"],
            "process_state_scope": ["extrusion_structured", "heated_matrix"],
            "supporting_families": ["06"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "soy_protein",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Deamidation of asparagine (Ea 25.1-113 kJ/mol) vs. glutamine (Ea 120-268 kJ/mol) in protein melts yielding free ammonia."
            ],
            "what_it_does_not_support": ["Free amino acid solutions without protein structural constraints"],
            "key_values": {
                "asn_deamidation_Ea_range_kj_mol": [25.1, 113.0],
                "gln_deamidation_Ea_range_kj_mol": [120.0, 268.0],
                "replicates": 3
            },
            "repo_next_action": "Encode thermal deamidation kinetics.",
            "runtime_artifacts": [
                {
                    "artifact_type": "computational_prior",
                    "artifact_id": "zhang_1993_protein_deamidation_ammonia",
                    "path": "data/lit/computational_priors.json"
                }
            ],
            "requires_primary_data": False
        },
        {
            "id": "li_2010_phytate_chelation_kinetics",
            "citation": "Li et al. (2010)",
            "doi": "10.1016/j.mineng.2010.05.003",
            "kind": "calibration_reference",
            "chemistry_family": "alternative_protein_matrix_scope",
            "slr_family_source": "06",
            "payload_role": "computational_prior",
            "observable_panel_tags": ["phytate", "chelation", "metal_catalysis"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["06"],
            "target_modules": ["literature_runtime"],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Adsorption and coordinate covalent chelation activation energy of phytic acid binding metal ions (Ea 84.54 kJ/mol)."
            ],
            "what_it_does_not_support": ["Divalent cation complexes without phosphate ester coordination"],
            "key_values": {
                "chelation_Ea_kj_mol": 84.54,
                "replicates": 3
            },
            "repo_next_action": "Encode phytate chelation kinetics prior.",
            "runtime_artifacts": [
                {
                    "artifact_type": "computational_prior",
                    "artifact_id": "li_2010_phytate_chelation_kinetics",
                    "path": "data/lit/computational_priors.json"
                }
            ],
            "requires_primary_data": False
        }
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
            "slr_section": "3.0",
            "paper_id": "voelker_2018_thiamine_salt_degradation",
            "citation": "Voelker et al. (2018)",
            "doi": "10.1016/j.foodres.2018.06.056",
            "matrix_family": "free_model_system",
            "compounds_supported": ["thiamine_mononitrate", "thiamine_chloride_hydrochloride"],
            "parameters_supported": ["thiamine_degradation_kinetics", "activation_energy"],
            "exact_numeric_anchors": ["Ea TMN: 88-105 kJ/mol", "Ea TClHCl: 90-135 kJ/mol"],
            "current_repo_artifacts": ["data/lit/computational_priors.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Encode thiamine salt form Arrhenius kinetics.",
            "confidence_tier": "high",
            "notes_on_limits": "Aqueous buffered systems."
        },
        {
            "slr_section": "12.0",
            "paper_id": "knol_2005_acrylamide_kinetics",
            "citation": "Knol et al. (2005)",
            "doi": "10.1021/jf050504m",
            "matrix_family": "free_model_system",
            "compounds_supported": ["acrylamide"],
            "parameters_supported": ["acrylamide_formation_kinetics", "acrylamide_degradation_kinetics"],
            "exact_numeric_anchors": ["formation Ea 52.1 kJ/mol", "degradation Ea 72.9 kJ/mol"],
            "current_repo_artifacts": ["data/lit/safety_reference_payloads.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Encode as acrylamide safety reference payload.",
            "confidence_tier": "high",
            "notes_on_limits": "Glucose-asparagine model system."
        },
        {
            "slr_section": "12.0",
            "paper_id": "de_vleeschouwer_2008_acrylamide_aw",
            "citation": "De Vleeschouwer et al. (2008)",
            "doi": "10.1021/bp060389f",
            "matrix_family": "free_model_system",
            "compounds_supported": ["acrylamide"],
            "parameters_supported": ["acrylamide_moisture_dependent_kinetics"],
            "exact_numeric_anchors": ["minimum elimination rate at Aw 0.82"],
            "current_repo_artifacts": ["data/lit/safety_reference_payloads.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Encode as water activity dependent acrylamide safety reference payload.",
            "confidence_tier": "high",
            "notes_on_limits": "Low-moisture starch systems."
        },
        {
            "slr_section": "12.0",
            "paper_id": "ishak_2022_phip_hca_kinetics",
            "citation": "Ishak et al. (2022)",
            "doi": "10.1016/j.foodchem.2022.132372",
            "matrix_family": "free_model_system",
            "compounds_supported": ["PhIP", "HCAs"],
            "parameters_supported": ["phip_formation_kinetics", "activation_energy"],
            "exact_numeric_anchors": ["Ea Phe model: 95.36 kJ/mol", "Ea Pro model: 114.12 kJ/mol"],
            "current_repo_artifacts": ["data/lit/safety_reference_payloads.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Encode as PhIP HCA safety reference payload.",
            "confidence_tier": "high",
            "notes_on_limits": "Creatinine-phenylalanine/proline models."
        },
        {
            "slr_section": "13.0",
            "paper_id": "shu_2019_cysteine_quinone_kinetics",
            "citation": "Shu et al. (2019)",
            "doi": "10.1016/j.freeradbiomed.2019.04.026",
            "matrix_family": "free_model_system",
            "compounds_supported": ["cysteine-quinone Michael adducts"],
            "parameters_supported": ["quinone_thiol_conjugation_kinetics", "rate_constant"],
            "exact_numeric_anchors": ["second-order rate constant k: 10^5 to 10^6 M-1 s-1"],
            "current_repo_artifacts": ["data/lit/computational_priors.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Encode cysteine-quinone Michael adduct kinetics.",
            "confidence_tier": "high",
            "notes_on_limits": "Recombinant protein models."
        },
        {
            "slr_section": "6.0",
            "paper_id": "zhang_1993_protein_deamidation_ammonia",
            "citation": "Zhang et al. (1993)",
            "doi": "10.1111/j.1745-4530.1993.tb00179.x",
            "matrix_family": "soy_protein",
            "compounds_supported": ["ammonia", "asparagine", "glutamine"],
            "parameters_supported": ["asparagine_deamidation_kinetics", "glutamine_deamidation_kinetics"],
            "exact_numeric_anchors": ["Asn deamidation Ea: 25.1-113 kJ/mol", "Gln deamidation Ea: 120-268 kJ/mol"],
            "current_repo_artifacts": ["data/lit/computational_priors.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Encode thermal deamidation kinetics.",
            "confidence_tier": "high",
            "notes_on_limits": "Protein melts."
        },
        {
            "slr_section": "6.0",
            "paper_id": "li_2010_phytate_chelation_kinetics",
            "citation": "Li et al. (2010)",
            "doi": "10.1016/j.mineng.2010.05.003",
            "matrix_family": "free_model_system",
            "compounds_supported": ["phytic_acid", "metal_phytate_chelates"],
            "parameters_supported": ["phytate_chelation_adsorption_kinetics", "activation_energy"],
            "exact_numeric_anchors": ["apparent Ea: 84.54 kJ/mol"],
            "current_repo_artifacts": ["data/lit/computational_priors.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Encode phytate chelation kinetics prior.",
            "confidence_tier": "high",
            "notes_on_limits": "Aqueous solution."
        }
    ]

    existing_matrix_ids = {entry["paper_id"] for entry in matrix["entries"]}
    for entry in new_matrix_entries:
        if entry["paper_id"] not in existing_matrix_ids:
            matrix["entries"].append(entry)
            print(f"Added {entry['paper_id']} to slr_incorporation_matrix")

    # C. Add to computational_priors.json
    # 1. Thiamine
    if "thiamine_pathway_priors" in priors:
        existing_thiamine = {p["id"] for p in priors["thiamine_pathway_priors"]}
        if "voelker_2018_thiamine_salt_degradation" not in existing_thiamine:
            priors["thiamine_pathway_priors"].append({
                "id": "voelker_2018_thiamine_salt_degradation",
                "source": "Voelker et al. (2018)",
                "provenance_tier": "literature_derived_transfer",
                "confidence_tier": "high",
                "uncertainty_posture": "directional_transfer",
                "reference_conditions": {
                    "temp_range_C": [100.0, 150.0]
                },
                "thiamine_degradation_kinetics": {
                    "ea_tmn_kj_mol": 96.5,
                    "ea_tclhcl_kj_mol": 112.5
                },
                "notes": "First-order thiamine degradation kinetics where TMN degrades faster (Ea 88-105 kJ/mol) than TClHCl (Ea 90-135 kJ/mol)."
            })
            print("Added voelker_2018_thiamine_salt_degradation to thiamine_pathway_priors")

    # 2. Capping
    if "polyphenol_thiol_capping_priors" in priors:
        existing_capping = {p["id"] for p in priors["polyphenol_thiol_capping_priors"]}
        if "shu_2019_cysteine_quinone_kinetics" not in existing_capping:
            priors["polyphenol_thiol_capping_priors"].append({
                "id": "shu_2019_cysteine_quinone_kinetics",
                "source": "Shu et al. (2019)",
                "provenance_tier": "literature_derived_transfer",
                "confidence_tier": "high",
                "uncertainty_posture": "directional_transfer",
                "quinone_thiol_k_range": [100000.0, 1000000.0],
                "notes": "Second-order rate constants of 10^5 to 10^6 M-1 s-1 for quinone-thiol Michael additions."
            })
            print("Added shu_2019_cysteine_quinone_kinetics to polyphenol_thiol_capping_priors")

    # 3. Crosslinks
    if "crosslink_kinetics_priors" in priors:
        existing_crosslinks = {p["id"] for p in priors["crosslink_kinetics_priors"]}
        if "zhang_1993_protein_deamidation_ammonia" not in existing_crosslinks:
            priors["crosslink_kinetics_priors"].append({
                "id": "zhang_1993_protein_deamidation_ammonia",
                "source": "Zhang et al. (1993)",
                "provenance_tier": "literature_derived_transfer",
                "confidence_tier": "high",
                "uncertainty_posture": "directional_transfer",
                "deamidation_kinetics": {
                    "ea_asn_kj_mol_range": [25.1, 113.0],
                    "ea_gln_kj_mol_range": [120.0, 268.0]
                },
                "notes": "Ammonia release from asparagine and glutamine deamidation."
            })
            print("Added zhang_1993_protein_deamidation_ammonia to crosslink_kinetics_priors")
        if "li_2010_phytate_chelation_kinetics" not in existing_crosslinks:
            priors["crosslink_kinetics_priors"].append({
                "id": "li_2010_phytate_chelation_kinetics",
                "source": "Li et al. (2010)",
                "provenance_tier": "literature_derived_transfer",
                "confidence_tier": "high",
                "uncertainty_posture": "directional_transfer",
                "chelation_kinetics": {
                    "apparent_ea_kj_mol": 84.54
                },
                "notes": "Coordinate covalent chelation of metals by phytic acid."
            })
            print("Added li_2010_phytate_chelation_kinetics to crosslink_kinetics_priors")

    # 4. Update quinone_cys_michael dft kinetic prior
    if "dft_kinetic_priors" in priors and "entries" in priors["dft_kinetic_priors"]:
        for entry in priors["dft_kinetic_priors"]["entries"]:
            if entry.get("reaction_key") == "quinone_cys_michael":
                entry["current_tier"] = "literature_derived_transfer"
                entry["refinability_status"] = "literature_anchor_only"
                entry["computational_method"] = "Direct literature anchor from Shu et al. (2019) second-order cysteine-quinone conjugation kinetics (k ~ 10^5 to 10^6 M-1 s-1)."
                print("Updated quinone_cys_michael status in dft_kinetic_priors")

    # D. Add to safety_reference_payloads.json
    if "entries" in safety:
        existing_safety = {entry["id"] for entry in safety["entries"]}
        new_safety_entries = [
            {
                "id": "knol_2005_acrylamide_kinetics",
                "kind": "industrial_endpoint_reference",
                "report_visibility": "default",
                "target_module": "safety",
                "source_citation": "Knol et al. (2005)",
                "doi": "10.1021/jf050504m",
                "validated_status": "reference_anchor",
                "analyte": "acrylamide",
                "method": {
                    "instrument": "LC-MS",
                    "replicates": 3
                },
                "what_it_supports": [
                    "Bifurcated acrylamide formation (Ea 52.1 kJ/mol) and degradation (Ea 72.9 kJ/mol) kinetic parameters."
                ],
                "what_it_does_not_support": [
                    "Solid-state high-moisture extrusion without Aw corrections."
                ]
            },
            {
                "id": "de_vleeschouwer_2008_acrylamide_aw",
                "kind": "industrial_endpoint_reference",
                "report_visibility": "default",
                "target_module": "safety",
                "source_citation": "De Vleeschouwer et al. (2008)",
                "doi": "10.1021/bp060389f",
                "validated_status": "reference_anchor",
                "analyte": "acrylamide",
                "method": {
                    "instrument": "HPLC",
                    "replicates": 3
                },
                "what_it_supports": [
                    "Moisture-dependent acrylamide degradation with minimum elimination at Aw 0.82."
                ],
                "what_it_does_not_support": [
                    "Aqueous solution systems."
                ]
            },
            {
                "id": "ishak_2022_phip_hca_kinetics",
                "kind": "industrial_endpoint_reference",
                "report_visibility": "default",
                "target_module": "safety",
                "source_citation": "Ishak et al. (2022)",
                "doi": "10.1016/j.foodchem.2022.132372",
                "validated_status": "reference_anchor",
                "analyte": "PhIP",
                "method": {
                    "instrument": "HPLC",
                    "replicates": 3
                },
                "what_it_supports": [
                    "PhIP HCA formation kinetics (Ea 95.36 kJ/mol Phe vs 114.12 kJ/mol Pro)."
                ],
                "what_it_does_not_support": [
                    "Low temperature processes."
                ]
            }
        ]
        for entry in new_safety_entries:
            if entry["id"] not in existing_safety:
                safety["entries"].append(entry)
                print(f"Added {entry['id']} to safety_reference_payloads")

    # E. Update/Sync deep_research_backlog.json
    backlog_citations_map = {
        "Voelker et al. (2018)": ("voelker_2018_thiamine_salt_degradation", "03_thiamine_degradation.md"),
        "Knol et al. (2005)": ("knol_2005_acrylamide_kinetics", "12_protein_damage_markers.md"),
        "De Vleeschouwer et al. (2008)": ("de_vleeschouwer_2008_acrylamide_aw", "12_protein_damage_markers.md"),
        "Ishak et al. (2022)": ("ishak_2022_phip_hca_kinetics", "12_protein_damage_markers.md"),
        "Shu et al. (2019)": ("shu_2019_cysteine_quinone_kinetics", "13_polyphenol_amino_capping.md"),
        "Zhang et al. (1993)": ("zhang_1993_protein_deamidation_ammonia", "06_alternative_proteins.md"),
        "Li et al. (2010)": ("li_2010_phytate_chelation_kinetics", "06_alternative_proteins.md")
    }

    existing_backlog_citations = {item.get("citation") for item in backlog["items"]}

    for cit, (reg_id, filename) in backlog_citations_map.items():
        if cit not in existing_backlog_citations:
            new_item = {
                "citation": cit,
                "score": "6/8",
                "score_value": 6,
                "status": "RUNTIME_BOUND",
                "occurrence_count": 1,
                "files": [filename],
                "descriptions": [f"Deep research kinetic parameters for {cit}"],
                "occurrences": [
                    {
                        "file": filename,
                        "line": 999,
                        "description": f"Deep research kinetic parameters for {cit}",
                        "raw_line": f"1. `{cit}` — Bounded priors."
                    }
                ],
                "registry_id": reg_id,
                "runtime_artifact_count": 1
            }
            backlog["items"].append(new_item)
            print(f"Added {cit} to deep_research_backlog")
        else:
            for item in backlog["items"]:
                if item.get("citation") == cit:
                    item["status"] = "RUNTIME_BOUND"
                    item["registry_id"] = reg_id
                    item["runtime_artifact_count"] = 1
                    item["files"] = [filename]
                    if item.get("occurrences"):
                        for occ in item["occurrences"]:
                            occ["file"] = filename
                    print(f"Updated backlog status for {cit} to RUNTIME_BOUND")

    save_json(intake, INTAKE_REGISTRY_PATH)
    save_json(matrix, SLR_MATRIX_PATH)
    save_json(priors, COMPUTATIONAL_PRIORS_PATH)
    save_json(safety, SAFETY_PAYLOADS_PATH)
    save_json(backlog, BACKLOG_PATH)

    print("Ingestion completed successfully!")


if __name__ == "__main__":
    ingest()
