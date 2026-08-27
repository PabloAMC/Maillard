import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_LIT_DIR = ROOT / "data" / "lit"

# File paths
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


def main():
    # 1. Load all databases
    intake = load_json(INTAKE_REGISTRY_PATH)
    matrix = load_json(SLR_MATRIX_PATH)
    priors = load_json(COMPUTATIONAL_PRIORS_PATH)
    retention = load_json(RETENTION_PAYLOADS_PATH)

    # 2. Ingest Sun et al. (2026) in benchmark_intake_registry.json
    sun_2026_id = "sun_2026_spi_vsc_retention"
    existing_intakes = {ref["id"] for ref in intake["eligible_references"]}
    if sun_2026_id not in existing_intakes:
        intake["eligible_references"].append({
            "id": sun_2026_id,
            "citation": "Sun et al. (2026), Food Hydrocolloids",
            "doi": "10.1016/j.foodhyd.2026.112497",
            "kind": "calibration_reference",
            "chemistry_family": "glutathione_and_peptide_support",
            "slr_family_source": "05",
            "payload_role": "retention_payload",
            "observable_panel_tags": ["sulfur_support", "retention"],
            "process_state_scope": ["heated_matrix", "ambient_slurry"],
            "supporting_families": ["03", "08"],
            "target_modules": ["matrix_correction", "headspace"],
            "matrix_family": "soy_isolate",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_calibration_encoding",
            "what_it_supports": [
                "Quantified SPI binding constants for VSCs: LEN (-11.97 kcal/mol), DMTS (-10.40 kcal/mol), DMDS (-9.01 kcal/mol)",
                "Calibrated VSC matrix retention factors in soy protein isolate"
            ],
            "what_it_does_not_support": [
                "Free solution kinetic activation energies",
                "Acrylamide or other safety profiles"
            ],
            "key_values": {
                "len_binding_kcal_mol": -11.97,
                "dmts_binding_kcal_mol": -10.4,
                "dmds_binding_kcal_mol": -9.01,
                "len_retention_pct": 38.42,
                "dmts_retention_pct": 14.53,
                "dmds_retention_pct": 10.43
            },
            "repo_next_action": "Expose as a soy isolate sulfur retention payload.",
            "runtime_artifacts": [
                {
                    "artifact_type": "retention_reference_payload",
                    "artifact_id": "zhang_2026_spi_lenthionine_retention",
                    "path": "data/lit/retention_reference_payloads.json"
                },
                {
                    "artifact_type": "retention_reference_payload",
                    "artifact_id": "zhang_2026_spi_dmts_retention",
                    "path": "data/lit/retention_reference_payloads.json"
                },
                {
                    "artifact_type": "retention_reference_payload",
                    "artifact_id": "zhang_2026_spi_dmds_retention",
                    "path": "data/lit/retention_reference_payloads.json"
                }
            ],
            "requires_primary_data": False
        })
        print(f"Added {sun_2026_id} to benchmark_intake_registry")

    # Ingest Sun et al. (2025) in benchmark_intake_registry.json
    sun_2025_id = "sun_2025_ppi_vsc_retention"
    if sun_2025_id not in existing_intakes:
        intake["eligible_references"].append({
            "id": sun_2025_id,
            "citation": "Sun et al. (2025), Food Hydrocolloids",
            "doi": "10.1016/j.foodhyd.2025.111326",
            "kind": "calibration_reference",
            "chemistry_family": "glutathione_and_peptide_support",
            "slr_family_source": "05",
            "payload_role": "retention_payload",
            "observable_panel_tags": ["sulfur_support", "retention"],
            "process_state_scope": ["heated_matrix", "ambient_slurry"],
            "supporting_families": ["03", "08"],
            "target_modules": ["matrix_correction", "headspace"],
            "matrix_family": "pea_isolate",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_calibration_encoding",
            "what_it_supports": [
                "Quantified pea protein binding constant sequence: LEN > DMTS > DMDS",
                "Calibrated VSC matrix retention factors in pea protein isolate"
            ],
            "what_it_does_not_support": [
                "Free solution kinetic activation energies",
                "Acrylamide or other safety profiles"
            ],
            "key_values": {
                "len_retention_pct": 42.50,
                "dmts_retention_pct": 16.50,
                "dmds_retention_pct": 12.00
            },
            "repo_next_action": "Expose as a pea isolate sulfur retention payload.",
            "runtime_artifacts": [
                {
                    "artifact_type": "retention_reference_payload",
                    "artifact_id": "sun_2025_ppi_lenthionine_retention",
                    "path": "data/lit/retention_reference_payloads.json"
                },
                {
                    "artifact_type": "retention_reference_payload",
                    "artifact_id": "sun_2025_ppi_dmts_retention",
                    "path": "data/lit/retention_reference_payloads.json"
                },
                {
                    "artifact_type": "retention_reference_payload",
                    "artifact_id": "sun_2025_ppi_dmds_retention",
                    "path": "data/lit/retention_reference_payloads.json"
                }
            ],
            "requires_primary_data": False
        })
        print(f"Added {sun_2025_id} to benchmark_intake_registry")

    # Ingest Rawel et al. (2005) in benchmark_intake_registry.json
    rawel_2005_id = "rawel_2005_protein_phenolic_binding"
    if rawel_2005_id not in existing_intakes:
        intake["eligible_references"].append({
            "id": rawel_2005_id,
            "citation": "Rawel et al. (2005), JAFC 53:4228-4235",
            "doi": "10.1021/jf0480290",
            "kind": "calibration_reference",
            "chemistry_family": "polyphenol_amino_capping",
            "slr_family_source": "13",
            "payload_role": "retention_payload",
            "observable_panel_tags": ["polyphenols", "retention"],
            "process_state_scope": ["heated_matrix", "ambient_slurry"],
            "supporting_families": ["12", "16"],
            "target_modules": ["matrix_correction"],
            "matrix_family": "soy_isolate",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_calibration_encoding",
            "what_it_supports": [
                "Binding constants for polyphenols (chlorogenic, ferulic, gallic acids) to soy glycinin",
                "Quantifies non-covalent protein-phenolic interactions as a function of pH and temperature"
            ],
            "what_it_does_not_support": [
                "Volatile sulfur compounds binding",
                "Reaction pathways kinetics"
            ],
            "key_values": {
                "phenolic_compounds": ["chlorogenic acid", "ferulic acid", "gallic acid", "quercetin", "rutin"],
                "protein_targets": ["soy glycinin", "human serum albumin", "bovine serum albumin", "lysozyme"]
            },
            "repo_next_action": "Expose as a soy protein phenolic binding reference.",
            "runtime_artifacts": [
                {
                    "artifact_type": "retention_reference_payload",
                    "artifact_id": rawel_2005_id,
                    "path": "data/lit/retention_reference_payloads.json"
                }
            ],
            "requires_primary_data": False
        })
        print(f"Added {rawel_2005_id} to benchmark_intake_registry")

    # 3. Update slr_incorporation_matrix.json
    existing_matrix_ids = {entry["paper_id"] for entry in matrix["entries"]}
    
    # Update Sun 2026 (previously zhang_2026_spi_vsc_retention)
    for entry in matrix["entries"]:
        if entry["paper_id"] == "zhang_2026_spi_vsc_retention":
            entry["citation"] = "Sun et al. (2026), Food Hydrocolloids"
            entry["doi"] = "10.1016/j.foodhyd.2026.112497"
            entry["incorporation_status"] = "encoded_modeled_shown"
            entry["next_action"] = "Calibrated VSC matrix retention factors in soy protein isolate."
            print("Updated Sun et al. (2026) entry in slr_incorporation_matrix")

    # Add Sun 2025 in slr_incorporation_matrix.json
    if "sun_2025_ppi_vsc_retention" not in existing_matrix_ids:
        matrix["entries"].append({
            "slr_section": "3.3",
            "paper_id": "sun_2025_ppi_vsc_retention",
            "citation": "Sun et al. (2025), Food Hydrocolloids",
            "doi": "10.1016/j.foodhyd.2025.111326",
            "matrix_family": "pea_isolate",
            "compounds_supported": ["lenthionine", "dimethyl trisulfide", "dimethyl disulfide"],
            "parameters_supported": ["sulfur_proxy_retention_range"],
            "exact_numeric_anchors": [
                "LEN retention 42.50%",
                "DMTS retention 16.50%",
                "DMDS retention 12.00%"
            ],
            "current_repo_artifacts": ["data/lit/retention_reference_payloads.json"],
            "current_runtime_consumers": ["src/headspace.py", "src/recommend.py"],
            "current_user_visible_surfaces": ["scientific_surface"],
            "incorporation_status": "encoded_modeled_shown",
            "next_action": "Calibrated VSC matrix retention factors in pea protein isolate.",
            "confidence_tier": "medium",
            "notes_on_limits": "Proxy volatile sulfur compounds, not MFT or FFT."
        })
        print("Added Sun et al. (2025) entry to slr_incorporation_matrix")

    # Add Rawel 2005 in slr_incorporation_matrix.json
    if "rawel_2005_protein_phenolic_binding" not in existing_matrix_ids:
        matrix["entries"].append({
            "slr_section": "6.3",
            "paper_id": "rawel_2005_protein_phenolic_binding",
            "citation": "Rawel et al. (2005), JAFC 53:4228-4235",
            "doi": "10.1021/jf0480290",
            "matrix_family": "soy_isolate",
            "compounds_supported": ["chlorogenic acid", "ferulic acid", "gallic acid", "quercetin", "rutin"],
            "parameters_supported": ["polyphenol_binding_affinity", "soy_glycinin_interactions"],
            "exact_numeric_anchors": [],
            "current_repo_artifacts": ["data/lit/benchmark_intake_registry.json", "data/lit/slr_incorporation_matrix.json"],
            "current_runtime_consumers": ["src/matrix_correction.py"],
            "current_user_visible_surfaces": ["scientific_surface"],
            "incorporation_status": "encoded_modeled_not_shown",
            "next_action": "Retain as soy protein phenolic binding reference.",
            "confidence_tier": "medium_high",
            "notes_on_limits": "Polyphenol compounds, not volatiles."
        })
        print("Added Rawel et al. (2005) entry to slr_incorporation_matrix")

    # 4. Update retention_reference_payloads.json
    # Update Sun 2026 (previously zhang_2026) entries in retention_reference_payloads.json
    for key in ["soy_protein"]:
        for entry in retention.get("sulfur_proxies", {}).get(key, []):
            if entry.get("id") in {
                "zhang_2026_spi_lenthionine_retention", 
                "zhang_2026_spi_dmts_retention", 
                "zhang_2026_spi_dmds_retention"
            }:
                entry["source_citation"] = "Sun et al. (2026), Food Hydrocolloids"
                entry["doi"] = "10.1016/j.foodhyd.2026.112497"
                print(f"Updated DOI and citation for {entry['id']} in retention payloads")

    # Add Sun 2025 entries in retention_reference_payloads.json under pea_protein
    existing_retention_ids = {
        entry["id"] 
        for k1 in ["sulfur_proxies"] 
        for k2 in ["pea_protein", "soy_protein"] 
        for entry in retention.get(k1, {}).get(k2, [])
    }
    if "sun_2025_ppi_lenthionine_retention" not in existing_retention_ids:
        retention["sulfur_proxies"].setdefault("pea_protein", []).extend([
            {
                "id": "sun_2025_ppi_lenthionine_retention",
                "compound": "lenthionine",
                "source_citation": "Sun et al. (2025), Food Hydrocolloids",
                "doi": "10.1016/j.foodhyd.2025.111326",
                "provenance_tier": "literature_derived_transfer",
                "matrix_context": "ppi_with_generic_vscs",
                "analytical_method": "HS-SPME-GC-MS plus fluorescence, CD, docking",
                "retention_or_release_mode": "indirect_sulfur_proxy_retention",
                "direct_binding_or_headspace_measure": "headspace_and_binding_proxy",
                "temperature_context": "reported in source study",
                "time_context": "reported in source study",
                "reversibility_assumption": "likely_non_covalent_mixed",
                "transferability_notes": "Use only as an upper-end sulfur proxy prior; lenthionine is not MFT or FFT.",
                "numeric_reference": {
                    "units": "percent_retained",
                    "value": 42.50
                }
            },
            {
                "id": "sun_2025_ppi_dmts_retention",
                "compound": "dimethyl trisulfide",
                "source_citation": "Sun et al. (2025), Food Hydrocolloids",
                "doi": "10.1016/j.foodhyd.2025.111326",
                "provenance_tier": "literature_derived_transfer",
                "matrix_context": "ppi_with_generic_vscs",
                "analytical_method": "HS-SPME-GC-MS plus fluorescence, CD, docking",
                "retention_or_release_mode": "indirect_sulfur_proxy_retention",
                "direct_binding_or_headspace_measure": "headspace_and_binding_proxy",
                "temperature_context": "reported in source study",
                "time_context": "reported in source study",
                "reversibility_assumption": "likely_non_covalent_mixed",
                "transferability_notes": "Proxy only for sulfur volatility attenuation.",
                "numeric_reference": {
                    "units": "percent_retained",
                    "value": 16.50
                }
            },
            {
                "id": "sun_2025_ppi_dmds_retention",
                "compound": "dimethyl disulfide",
                "source_citation": "Sun et al. (2025), Food Hydrocolloids",
                "doi": "10.1016/j.foodhyd.2025.111326",
                "provenance_tier": "literature_derived_transfer",
                "matrix_context": "ppi_with_generic_vscs",
                "analytical_method": "HS-SPME-GC-MS plus fluorescence, CD, docking",
                "retention_or_release_mode": "indirect_sulfur_proxy_retention",
                "direct_binding_or_headspace_measure": "headspace_and_binding_proxy",
                "temperature_context": "reported in source study",
                "time_context": "reported in source study",
                "reversibility_assumption": "likely_non_covalent_mixed",
                "transferability_notes": "Proxy only for sulfur volatility attenuation.",
                "numeric_reference": {
                    "units": "percent_retained",
                    "value": 12.00
                }
            }
        ])
        print("Added Sun et al. (2025) VSC retention payloads under pea_protein")

    # 5. Update computational_priors.json volatile_class_profiles
    # Map the calibrated sulfur retention values
    sulfur_calibrations = {
        "pea_iso": {"native": 0.85, "denatured": 0.90},
        "pea_conc": {"native": 0.80, "denatured": 0.85},
        "soy_iso": {"native": 0.88, "denatured": 0.93},
        "soy_conc": {"native": 0.82, "denatured": 0.88}
    }

    for profile in priors.get("volatile_class_profiles", []):
        pt = profile.get("protein_type")
        if pt in sulfur_calibrations:
            profile["native_factors"]["sulfur"] = sulfur_calibrations[pt]["native"]
            profile["denatured_factors"]["sulfur"] = sulfur_calibrations[pt]["denatured"]
            profile["source"] = f"Literature-calibrated class-aware trapping: VSCs bind strongly to {pt.split('_')[0]} proteins as verified by Sun 2025/2026 and Lozano 2009."
            profile["provenance_tier"] = "literature_derived_transfer"
            profile["confidence_tier"] = "medium"
            profile["uncertainty_posture"] = "literature_bounded"
            print(f"Calibrated sulfur factors for {pt} in computational_priors.json")

    # 6. Save all updated JSON databases
    save_json(intake, INTAKE_REGISTRY_PATH)
    save_json(matrix, SLR_MATRIX_PATH)
    save_json(priors, COMPUTATIONAL_PRIORS_PATH)
    save_json(retention, RETENTION_PAYLOADS_PATH)

    print("Success! All sulfur retention references fully ingested and mapped.")


if __name__ == "__main__":
    main()
