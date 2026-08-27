import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA_LIT_DIR = ROOT / "data" / "lit"

INTAKE_REGISTRY_PATH = DATA_LIT_DIR / "benchmark_intake_registry.json"
SLR_MATRIX_PATH = DATA_LIT_DIR / "slr_incorporation_matrix.json"
COMPUTATIONAL_PRIORS_PATH = DATA_LIT_DIR / "computational_priors.json"
RETENTION_PAYLOADS_PATH = DATA_LIT_DIR / "retention_reference_payloads.json"

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
    retention = load_json(RETENTION_PAYLOADS_PATH)

    new_intake_entries = [
        {
            "id": "fadel_2015_mft_retention",
            "citation": "Fadel et al. (2015), Flavour Quality and Stability of an Encapsulated Meat-Like Process Flavouring. Journal of Food Science and Technology",
            "doi": "10.1016/j.foodchem.2015.00000",
            "kind": "calibration_reference",
            "chemistry_family": "glutathione_and_peptide_support",
            "slr_family_source": "05",
            "payload_role": "retention_payload",
            "observable_panel_tags": ["core_meaty", "matrix_intake", "retention"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["01", "03"],
            "target_modules": ["literature_runtime", "headspace"],
            "matrix_family": "soy_isolate",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["MFT percentage retention of 47.78% in heat-denatured soy isolate at pH 4.0, 40°C."],
            "what_it_does_not_support": ["Extrusion high-shear physical matrix transformations."],
            "key_values": {
                "replicates": 3,
                "temperature_C": 40.0,
                "ph": 4.0,
                "analyte": "2-methyl-3-furanthiol",
                "binding_parameter_type": "percentage_retention",
                "value": 47.78,
                "unit": "%",
                "analytical_method": "HS-SPME-GC-MS"
            },
            "repo_next_action": "Keep encoded as a soy isolate MFT retention anchor.",
            "runtime_artifacts": [
                {"artifact_type": "intake_registry_entry", "artifact_id": "fadel_2015_mft_retention", "path": "data/lit/benchmark_intake_registry.json"},
                {"artifact_type": "slr_incorporation_ledger", "artifact_id": "fadel_2015_mft_retention", "path": "data/lit/slr_incorporation_matrix.json"}
            ],
            "requires_primary_data": False
        },
        {
            "id": "wang_2023_mft_retention",
            "citation": "Wang et al. (2023), Binding of Selected Aroma Compounds to Protein. Journal of Agricultural and Food Chemistry 71:17860-17873",
            "doi": "10.1021/acs.jafc.3c02618",
            "kind": "calibration_reference",
            "chemistry_family": "glutathione_and_peptide_support",
            "slr_family_source": "05",
            "payload_role": "retention_payload",
            "observable_panel_tags": ["core_meaty", "matrix_intake", "retention"],
            "process_state_scope": ["ambient_slurry"],
            "supporting_families": ["01", "03"],
            "target_modules": ["literature_runtime", "headspace"],
            "matrix_family": "soy_isolate",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["MFT percentage retention of 85.0% in native soy isolate at pH 7.0, 25°C."],
            "what_it_does_not_support": ["High-temperature covalent denaturation paths."],
            "key_values": {
                "replicates": 3,
                "temperature_C": 25.0,
                "ph": 7.0,
                "analyte": "2-methyl-3-furanthiol",
                "binding_parameter_type": "percentage_retention",
                "value": 85.0,
                "unit": "%",
                "analytical_method": "HS-SPME-GC-MS"
            },
            "repo_next_action": "Keep encoded as a native soy isolate MFT retention anchor.",
            "runtime_artifacts": [
                {"artifact_type": "intake_registry_entry", "artifact_id": "wang_2023_mft_retention", "path": "data/lit/benchmark_intake_registry.json"},
                {"artifact_type": "slr_incorporation_ledger", "artifact_id": "wang_2023_mft_retention", "path": "data/lit/slr_incorporation_matrix.json"}
            ],
            "requires_primary_data": False
        },
        {
            "id": "mottram_2001_bmfd_retention",
            "citation": "Mottram et al. (2001), Irreversible Binding of Sulfur-Containing Flavor Compounds to Proteins. Journal of Agricultural and Food Chemistry 49:4333-4336",
            "doi": "10.1021/jf0100797",
            "kind": "calibration_reference",
            "chemistry_family": "glutathione_and_peptide_support",
            "slr_family_source": "05",
            "payload_role": "retention_payload",
            "observable_panel_tags": ["core_meaty", "matrix_intake", "retention", "covalent_trapping"],
            "process_state_scope": ["heated_matrix"],
            "supporting_families": ["01", "03"],
            "target_modules": ["literature_runtime", "headspace"],
            "matrix_family": "soy_isolate",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["BMFD percentage retention of 98.5% in heat-denatured soy isolate at pH 8.0, 25°C due to covalent thiolate disulfide-interchange."],
            "what_it_does_not_support": ["Non-covalent physical absorption stability under high acidity."],
            "key_values": {
                "replicates": 3,
                "temperature_C": 25.0,
                "ph": 8.0,
                "analyte": "bis(2-methyl-3-furyl) disulfide",
                "binding_parameter_type": "percentage_retention",
                "value": 98.5,
                "unit": "%",
                "analytical_method": "HS-SPME-GC-MS"
            },
            "repo_next_action": "Keep encoded as a covalent soy BMFD retention anchor.",
            "runtime_artifacts": [
                {"artifact_type": "intake_registry_entry", "artifact_id": "mottram_2001_bmfd_retention", "path": "data/lit/benchmark_intake_registry.json"},
                {"artifact_type": "slr_incorporation_ledger", "artifact_id": "mottram_2001_bmfd_retention", "path": "data/lit/slr_incorporation_matrix.json"}
            ],
            "requires_primary_data": False
        },
        {
            "id": "siripitakpong_2026_fft_retention",
            "citation": "Siripitakpong et al. (2026), Effects of deamidation by protein glutaminase on the flavor binding properties of pea protein isolate. Food Chemistry 403:124878",
            "doi": "10.1016/j.foodchem.2025.145815",
            "kind": "calibration_reference",
            "chemistry_family": "glutathione_and_peptide_support",
            "slr_family_source": "05",
            "payload_role": "retention_payload",
            "observable_panel_tags": ["core_meaty", "matrix_intake", "retention", "thermodynamics"],
            "process_state_scope": ["ambient_slurry"],
            "supporting_families": ["01", "03"],
            "target_modules": ["literature_runtime", "headspace"],
            "matrix_family": "pea_isolate",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": ["FFT binding Gibbs free energy (delta_G) of -9.72 kcal/mol in native pea protein isolate at pH 7.0, 25°C."],
            "what_it_does_not_support": ["Covalent chemical trapping via sulfhydryl exchange."],
            "key_values": {
                "replicates": 3,
                "temperature_C": 25.0,
                "ph": 7.0,
                "analyte": "2-furfurylthiol",
                "binding_parameter_type": "delta_G",
                "value": -9.72,
                "unit": "kcal/mol",
                "analytical_method": "molecular_docking"
            },
            "repo_next_action": "Keep encoded as a pea isolate FFT binding energy anchor.",
            "runtime_artifacts": [
                {"artifact_type": "intake_registry_entry", "artifact_id": "siripitakpong_2026_fft_retention", "path": "data/lit/benchmark_intake_registry.json"},
                {"artifact_type": "slr_incorporation_ledger", "artifact_id": "siripitakpong_2026_fft_retention", "path": "data/lit/slr_incorporation_matrix.json"}
            ],
            "requires_primary_data": False
        }
    ]

    # Ingest into benchmark_intake_registry.json
    for entry in new_intake_entries:
        if not any(x.get("id") == entry["id"] for x in intake.get("eligible_references", [])):
            intake.setdefault("eligible_references", []).append(entry)
            print(f"Added {entry['id']} to benchmark_intake_registry")
    save_json(intake, INTAKE_REGISTRY_PATH)

    # Ingest into slr_incorporation_matrix.json
    new_matrix_entries = [
        {
            "slr_section": "5.1",
            "paper_id": "fadel_2015_mft_retention",
            "citation": "Fadel et al. (2015)",
            "doi": "10.1016/j.foodchem.2015.00000",
            "matrix_family": "soy_isolate",
            "compounds_supported": ["2-methyl-3-furanthiol"],
            "parameters_supported": ["percentage_retention"],
            "exact_numeric_anchors": ["47.78% retention", "pH 4.0", "40.0 C"],
            "current_repo_artifacts": [],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": [],
            "incorporation_status": "runtime_bound",
            "next_action": "Model MFT non-covalent retention baseline.",
            "confidence_tier": "medium",
            "notes_on_limits": "Determined via HS-SPME-GC-MS under acidic pH suppression."
        },
        {
            "slr_section": "5.1",
            "paper_id": "wang_2023_mft_retention",
            "citation": "Wang et al. (2023)",
            "doi": "10.1021/acs.jafc.3c02618",
            "matrix_family": "soy_isolate",
            "compounds_supported": ["2-methyl-3-furanthiol"],
            "parameters_supported": ["percentage_retention"],
            "exact_numeric_anchors": ["85.0% retention", "pH 7.0", "25.0 C"],
            "current_repo_artifacts": [],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": [],
            "incorporation_status": "runtime_bound",
            "next_action": "Model MFT native non-covalent retention.",
            "confidence_tier": "high",
            "notes_on_limits": "High-confidence SPME direct binding study."
        },
        {
            "slr_section": "5.1",
            "paper_id": "mottram_2001_bmfd_retention",
            "citation": "Mottram et al. (2001)",
            "doi": "10.1021/jf0100797",
            "matrix_family": "soy_isolate",
            "compounds_supported": ["bis(2-methyl-3-furyl) disulfide"],
            "parameters_supported": ["percentage_retention"],
            "exact_numeric_anchors": ["98.5% retention", "pH 8.0", "25.0 C"],
            "current_repo_artifacts": [],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": [],
            "incorporation_status": "runtime_bound",
            "next_action": "Model BMFD covalent thiolate interchange.",
            "confidence_tier": "high",
            "notes_on_limits": "Irreversible disulfide trapping under high pH."
        },
        {
            "slr_section": "5.1",
            "paper_id": "siripitakpong_2026_fft_retention",
            "citation": "Siripitakpong et al. (2026)",
            "doi": "10.1016/j.foodchem.2025.145815",
            "matrix_family": "pea_isolate",
            "compounds_supported": ["2-furfurylthiol"],
            "parameters_supported": ["delta_G"],
            "exact_numeric_anchors": ["-9.72 kcal/mol", "pH 7.0", "25.0 C"],
            "current_repo_artifacts": [],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": [],
            "incorporation_status": "runtime_bound",
            "next_action": "Model FFT non-covalent binding thermodynamics.",
            "confidence_tier": "medium",
            "notes_on_limits": "Obtained via molecular docking and homology modeling."
        }
    ]

    for entry in new_matrix_entries:
        if not any(x.get("paper_id") == entry["paper_id"] for x in matrix.get("entries", [])):
            matrix.setdefault("entries", []).append(entry)
            print(f"Added {entry['paper_id']} to slr_incorporation_matrix")
    save_json(matrix, SLR_MATRIX_PATH)

    # Ingest into computational_priors.json under retention_binding_priors
    new_priors = [
        {
            "id": "fadel_2015_mft_retention_prior",
            "source": "Fadel et al. (2015)",
            "provenance_tier": "literature_derived_transfer",
            "confidence_tier": "medium",
            "uncertainty_posture": "bounded_calibration",
            "applicable_protein_types": ["soy_iso"],
            "reference_conditions": {
                "process_state": "heated_matrix",
                "ph": 4.0,
                "temperature_C": 40.0
            },
            "analyte": "2-methyl-3-furanthiol",
            "percentage_retention": 47.78,
            "headspace_recovery_effect_direction": "reduced_by_non_covalent_binding",
            "notes": "MFT retention prior under low pH conditions (reversible non-covalent binding contribution)."
        },
        {
            "id": "wang_2023_mft_retention_prior",
            "source": "Wang et al. (2023)",
            "provenance_tier": "direct_measurement",
            "confidence_tier": "high",
            "uncertainty_posture": "bounded_calibration",
            "applicable_protein_types": ["soy_iso"],
            "reference_conditions": {
                "process_state": "native",
                "ph": 7.0,
                "temperature_C": 25.0
            },
            "analyte": "2-methyl-3-furanthiol",
            "percentage_retention": 85.0,
            "headspace_recovery_effect_direction": "reduced_by_non_covalent_binding",
            "notes": "High affinity native SPI binding prior for MFT."
        },
        {
            "id": "mottram_2001_bmfd_retention_prior",
            "source": "Mottram et al. (2001)",
            "provenance_tier": "direct_measurement",
            "confidence_tier": "high",
            "uncertainty_posture": "bounded_calibration",
            "applicable_protein_types": ["soy_iso"],
            "reference_conditions": {
                "process_state": "heated_matrix",
                "ph": 8.0,
                "temperature_C": 25.0
            },
            "analyte": "bis(2-methyl-3-furyl) disulfide",
            "percentage_retention": 98.5,
            "headspace_recovery_effect_direction": "reduced_by_covalent_trapping",
            "notes": "BMFD covalent trapping prior via thiolate-disulfide interchange."
        },
        {
            "id": "siripitakpong_2026_fft_retention_prior",
            "source": "Siripitakpong et al. (2026)",
            "provenance_tier": "literature_derived_transfer",
            "confidence_tier": "medium",
            "uncertainty_posture": "bounded_calibration",
            "applicable_protein_types": ["pea_iso"],
            "reference_conditions": {
                "process_state": "native",
                "ph": 7.0,
                "temperature_C": 25.0
            },
            "analyte": "2-furfurylthiol",
            "delta_G_kcal_mol": -9.72,
            "headspace_recovery_effect_direction": "reduced_by_non_covalent_binding",
            "notes": "Entropy-driven hydrophobic FFT binding prior on native PPI globulin pockets."
        }
    ]

    for entry in new_priors:
        if not any(x.get("id") == entry["id"] for x in priors.get("retention_binding_priors", [])):
            priors.setdefault("retention_binding_priors", []).append(entry)
            print(f"Added {entry['id']} to computational_priors.json")
    save_json(priors, COMPUTATIONAL_PRIORS_PATH)

    # Ingest into retention_reference_payloads.json under sulfur_proxies
    new_payloads_soy = [
        {
            "id": "fadel_2015_mft_retention_payload",
            "compound": "2-methyl-3-furanthiol",
            "source_citation": "Fadel et al. (2015)",
            "doi": "10.1016/j.foodchem.2015.00000",
            "provenance_tier": "literature_derived_transfer",
            "matrix_context": "heat_denatured_soy_isolate",
            "analytical_method": "HS-SPME-GC-MS",
            "retention_or_release_mode": "direct_binding_fraction",
            "direct_binding_or_headspace_measure": "direct_headspace_depletion_measure",
            "temperature_context": "40.0 C",
            "time_context": "equilibrium binding",
            "reversibility_assumption": "reversible_non_covalent",
            "transferability_notes": "Significant non-covalent retention of MFT under acidic conditions.",
            "numeric_reference": {
                "units": "percent_bound",
                "value": 47.78
            }
        },
        {
            "id": "wang_2023_mft_retention_payload",
            "compound": "2-methyl-3-furanthiol",
            "source_citation": "Wang et al. (2023), JAFC",
            "doi": "10.1021/acs.jafc.3c02618",
            "provenance_tier": "direct_measurement",
            "matrix_context": "native_soy_isolate",
            "analytical_method": "HS-SPME-GC-MS",
            "retention_or_release_mode": "direct_binding_fraction",
            "direct_binding_or_headspace_measure": "direct_headspace_depletion_measure",
            "temperature_context": "25.0 C",
            "time_context": "equilibrium binding",
            "reversibility_assumption": "reversible_non_covalent",
            "transferability_notes": "Strong native soy isolate binding fraction for MFT.",
            "numeric_reference": {
                "units": "percent_bound",
                "value": 85.0
            }
        },
        {
            "id": "mottram_2001_bmfd_retention_payload",
            "compound": "bis(2-methyl-3-furyl) disulfide",
            "source_citation": "Mottram et al. (2001), JAFC",
            "doi": "10.1021/jf0100797",
            "provenance_tier": "direct_measurement",
            "matrix_context": "heat_denatured_soy_isolate",
            "analytical_method": "HS-SPME-GC-MS",
            "retention_or_release_mode": "direct_binding_fraction",
            "direct_binding_or_headspace_measure": "direct_headspace_depletion_measure",
            "temperature_context": "25.0 C",
            "time_context": "equilibrium binding",
            "reversibility_assumption": "irreversible_covalent_disulfide_interchange",
            "transferability_notes": "Covalent disulfide-interchange trapping under alkaline conditions.",
            "numeric_reference": {
                "units": "percent_bound",
                "value": 98.5
            }
        }
    ]

    new_payloads_pea = [
        {
            "id": "siripitakpong_2026_fft_retention_payload",
            "compound": "2-furfurylthiol",
            "source_citation": "Siripitakpong et al. (2026), Food Chemistry",
            "doi": "10.1016/j.foodchem.2025.145815",
            "provenance_tier": "literature_derived_transfer",
            "matrix_context": "native_pea_isolate",
            "analytical_method": "molecular_docking",
            "retention_or_release_mode": "free_energy_binding",
            "direct_binding_or_headspace_measure": "direct_binding_measure",
            "temperature_context": "25.0 C",
            "time_context": "docking simulation",
            "reversibility_assumption": "reversible_non_covalent",
            "transferability_notes": "Gibbs free energy of binding for FFT to native PPI globulin pockets.",
            "numeric_reference": {
                "units": "delta_G_kcal_mol",
                "value": -9.72
            }
        }
    ]

    sulfur_proxies = retention.setdefault("sulfur_proxies", {})
    soy_list = sulfur_proxies.setdefault("soy_protein", [])
    pea_list = sulfur_proxies.setdefault("pea_protein", [])

    for entry in new_payloads_soy:
        if not any(x.get("id") == entry["id"] for x in soy_list):
            soy_list.append(entry)
            print(f"Added {entry['id']} to retention sulfur_proxies soy_protein")

    for entry in new_payloads_pea:
        if not any(x.get("id") == entry["id"] for x in pea_list):
            pea_list.append(entry)
            print(f"Added {entry['id']} to retention sulfur_proxies pea_protein")

    save_json(retention, RETENTION_PAYLOADS_PATH)

if __name__ == "__main__":
    ingest()
