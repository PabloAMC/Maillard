import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_LIT_DIR = ROOT / "data" / "lit"

# File paths
INTAKE_REGISTRY_PATH = DATA_LIT_DIR / "benchmark_intake_registry.json"
SLR_MATRIX_PATH = DATA_LIT_DIR / "slr_incorporation_matrix.json"
COMPUTATIONAL_PRIORS_PATH = DATA_LIT_DIR / "computational_priors.json"
FLAVOR_PAYLOADS_PATH = DATA_LIT_DIR / "flavor_reference_payloads.json"
SAFETY_PAYLOADS_PATH = DATA_LIT_DIR / "safety_reference_payloads.json"
BACKLOG_PATH = DATA_LIT_DIR / "deep_research_backlog.json"


def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def save_json(data, path):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"Updated {path.name}")


def ingest_family06_07():
    # 1. Load all files
    intake = load_json(INTAKE_REGISTRY_PATH)
    matrix = load_json(SLR_MATRIX_PATH)
    priors = load_json(COMPUTATIONAL_PRIORS_PATH)
    flavor = load_json(FLAVOR_PAYLOADS_PATH)
    safety = load_json(SAFETY_PAYLOADS_PATH)
    backlog = load_json(BACKLOG_PATH)

    # 2. Define the new references for benchmark_intake_registry.json
    new_intake_entries = [
        {
          "id": "researchgate_2023_pea_aeda",
          "citation": "ResearchGate (2023)",
          "doi": "",
          "kind": "calibration_reference",
          "chemistry_family": "alternative_protein_matrix_scope",
          "slr_family_source": "06",
          "payload_role": "benchmark_intake",
          "observable_panel_tags": [
            "pea_protein",
            "aeda",
            "gc-o"
          ],
          "process_state_scope": [
            "heated_matrix"
          ],
          "supporting_families": [
            "08"
          ],
          "target_modules": [
            "matrix_correction",
            "reporting"
          ],
          "matrix_family": "pea_protein_isolate",
          "eligibility": "benchmark_eligible",
          "status": "ready_for_intake_encoding",
          "what_it_supports": [
            "Pea protein heated slurry AEDA and GC-O volatile ranking.",
            "Identifies (E,E)-2,4-heptadienal (FD 4096) and methoxypyrazine (FD 256) as active off-notes."
          ],
          "what_it_does_not_support": [
            "Absolute concentration values for meaty sulfur markers",
            "Extrusion process-shear texturization parameters"
          ],
          "key_values": {
            "ee_2_4_heptadienal_fd_factor": 4096.0,
            "methoxypyrazine_fd_factor": 256.0
          },
          "repo_next_action": "Expose as a pea protein heating volatile AEDA reference payload.",
          "runtime_artifacts": [
            {
              "artifact_type": "flavor_reference_payload",
              "artifact_id": "researchgate_2023_pea_aeda_anchor",
              "path": "data/lit/flavor_reference_payloads.json"
            }
          ],
          "requires_primary_data": False
        },
        {
          "id": "malia_2025_pea_free_sh_crosscheck",
          "citation": "Malia et al. (2025)",
          "doi": "",
          "kind": "calibration_reference",
          "chemistry_family": "alternative_protein_matrix_scope",
          "slr_family_source": "06",
          "payload_role": "benchmark_intake",
          "observable_panel_tags": [
            "pea_protein",
            "free_sh"
          ],
          "process_state_scope": [
            "heated_matrix"
          ],
          "supporting_families": [
            "12"
          ],
          "target_modules": [
            "matrix_correction",
            "literature_runtime"
          ],
          "matrix_family": "pea_protein_isolate",
          "eligibility": "benchmark_eligible",
          "status": "ready_for_intake_encoding",
          "what_it_supports": [
            "Pea protein free-SH levels and DSC thermal denaturation crosscheck.",
            "Cysteine accessibility baseline in heated matrices."
          ],
          "what_it_does_not_support": [
            "Absolute volatile yield curves",
            "Lysine-loss kinetics"
          ],
          "key_values": {
            "assay": "DTNB / Ellman",
            "free_sh_reported": True
          },
          "repo_next_action": "Use as a supportive pea free-SH behavior crosscheck.",
          "runtime_artifacts": [
            {
              "artifact_type": "process_state_calibration",
              "artifact_id": "malia_2025_pea_free_sh_crosscheck",
              "path": "data/lit/process_state_calibrations.json"
            }
          ],
          "requires_primary_data": False
        },
        {
          "id": "tran_2023_reducing_sugar_hme",
          "citation": "Tran et al. (2023 / PMC9846454)",
          "doi": "10.1016/j.carbpol.2022.120468",
          "kind": "calibration_reference",
          "chemistry_family": "carbonyl_donor_hierarchy",
          "slr_family_source": "07",
          "payload_role": "benchmark_intake",
          "observable_panel_tags": [
            "reducing_sugar",
            "starch",
            "hme"
          ],
          "process_state_scope": [
            "extrusion_structured"
          ],
          "supporting_families": [
            "09"
          ],
          "target_modules": [
            "literature_runtime",
            "matrix_correction"
          ],
          "matrix_family": "pea_protein_isolate",
          "eligibility": "benchmark_eligible",
          "status": "ready_for_intake_encoding",
          "what_it_supports": [
            "Starch reducing-sugar liberation during HME (maltose generation)",
            "Compares pea starch (18.4 umol/g) versus waxy maize starch (31.8 umol/g) yield parameters"
          ],
          "what_it_does_not_support": [
            "Free ribose or glucose addition rates",
            "Volatile recovery yields under vacuum packaging"
          ],
          "key_values": {
            "pea_maltose_liberation_umol_g": 18.4,
            "waxy_maize_maltose_liberation_umol_g": 31.8,
            "ti_starch_reduction_pct": 66.0
          },
          "repo_next_action": "Encode starch reducing-sugar liberation priors during HME.",
          "runtime_artifacts": [
            {
              "artifact_type": "computational_prior",
              "artifact_id": "tran_2023_starch_liberation_v1",
              "path": "data/lit/computational_priors.json",
              "sections": [
                "carbonyl_donor_priors"
              ]
            }
          ],
          "requires_primary_data": False
        },
        {
          "id": "henle_2005_glycinin_cml",
          "citation": "Henle (2005)",
          "doi": "10.1007/s00726-005-0187-z",
          "kind": "calibration_reference",
          "chemistry_family": "protein_damage_markers",
          "slr_family_source": "12",
          "payload_role": "benchmark_intake",
          "observable_panel_tags": [
            "safety",
            "cml",
            "glycinin"
          ],
          "process_state_scope": [
            "heated_matrix"
          ],
          "supporting_families": [
            "07"
          ],
          "target_modules": [
            "safety",
            "matrix_correction"
          ],
          "matrix_family": "soy_isolate",
          "eligibility": "benchmark_eligible",
          "status": "ready_for_intake_encoding",
          "what_it_supports": [
            "Quantified CML formation (2.3 umol/g protein) from heated soy glycinin.",
            "Dicarbonyl consumption modeling."
          ],
          "what_it_does_not_support": [
            "Acrylamide kinetics on wheat gluten matrices",
            "Free volatile headspace partition factors"
          ],
          "key_values": {
            "glycinin_cml_umol_g": 2.3
          },
          "repo_next_action": "Encode soy glycinin dicarbonyl consumption safety reference.",
          "runtime_artifacts": [
            {
              "artifact_type": "safety_reference",
              "artifact_id": "henle_2005_glycinin_cml",
              "path": "data/lit/safety_reference_payloads.json"
            }
          ],
          "requires_primary_data": False
        }
    ]

    existing_ids = {ref["id"] for ref in intake["eligible_references"]}
    for ref in new_intake_entries:
        if ref["id"] not in existing_ids:
            intake["eligible_references"].append(ref)
            print(f"Added {ref['id']} to benchmark_intake_registry")

    # 3. Define entries for slr_incorporation_matrix.json
    new_matrix_entries = [
        {
          "slr_section": "3.1",
          "paper_id": "researchgate_2023_pea_aeda",
          "citation": "ResearchGate (2023)",
          "doi": "",
          "matrix_family": "pea_protein_isolate",
          "compounds_supported": [
            "(E,E)-2,4-heptadienal",
            "3-isobutyl-2-methoxypyrazine"
          ],
          "parameters_supported": [
            "flavor_dilution_factor"
          ],
          "exact_numeric_anchors": [
            "heptadienal FD 4096",
            "methoxypyrazine FD 256"
          ],
          "current_repo_artifacts": [
            "data/lit/flavor_reference_payloads.json"
          ],
          "current_runtime_consumers": [
            "src/sensory_scoring.py"
          ],
          "current_user_visible_surfaces": [
            "sensory_reporting"
          ],
          "incorporation_status": "encoded_shown",
          "next_action": "Expose as a pea protein heating volatile AEDA reference payload.",
          "confidence_tier": "high",
          "notes_on_limits": "GC-O/AEDA semi-quantitative values."
        },
        {
          "slr_section": "3.1",
          "paper_id": "malia_2025_pea_free_sh_crosscheck",
          "citation": "Malia et al. (2025), Food Res Int",
          "doi": "",
          "matrix_family": "pea_protein_isolate",
          "compounds_supported": [
            "cysteine"
          ],
          "parameters_supported": [
            "free_sh_accessibility"
          ],
          "exact_numeric_anchors": [
            "free SH in pea protein"
          ],
          "current_repo_artifacts": [
            "data/lit/process_state_calibrations.json"
          ],
          "current_runtime_consumers": [
            "src/literature_runtime.py"
          ],
          "current_user_visible_surfaces": [
            "reporting"
          ],
          "incorporation_status": "encoded_shown",
          "next_action": "Use as a supportive pea free-SH behavior crosscheck.",
          "confidence_tier": "high",
          "notes_on_limits": "Free-SH determination by Ellman assay."
        },
        {
          "slr_section": "3.2",
          "paper_id": "tran_2023_reducing_sugar_hme",
          "citation": "Tran et al. (2023), Carb Polymers",
          "doi": "10.1016/j.carbpol.2022.120468",
          "matrix_family": "pea_protein_isolate",
          "compounds_supported": [
            "maltose"
          ],
          "parameters_supported": [
            "reducing_sugar_liberation"
          ],
          "exact_numeric_anchors": [
            "pea 18.4 umol maltose/g",
            "waxy maize 31.8 umol/g",
            "TI starch 66% reduction"
          ],
          "current_repo_artifacts": [
            "data/lit/computational_priors.json"
          ],
          "current_runtime_consumers": [
            "src/literature_runtime.py"
          ],
          "current_user_visible_surfaces": [
            "reporting"
          ],
          "incorporation_status": "encoded_shown",
          "next_action": "Encode starch reducing-sugar liberation priors during HME.",
          "confidence_tier": "high",
          "notes_on_limits": "HME starch structural breakdown."
        },
        {
          "slr_section": "3.2",
          "paper_id": "henle_2005_glycinin_cml",
          "citation": "Henle (2005), Amino Acids",
          "doi": "10.1007/s00726-005-0187-z",
          "matrix_family": "soy_isolate",
          "compounds_supported": [
            "carboxymethyllysine"
          ],
          "parameters_supported": [
            "cml_formation",
            "dicarbonyl_consumption"
          ],
          "exact_numeric_anchors": [
            "2.3 umol CML/g protein",
            "dicarbonyl consumption estimate"
          ],
          "current_repo_artifacts": [
            "data/lit/safety_reference_payloads.json"
          ],
          "current_runtime_consumers": [
            "src/literature_runtime.py"
          ],
          "current_user_visible_surfaces": [
            "safety_reporting"
          ],
          "incorporation_status": "encoded_shown",
          "next_action": "Encode soy glycinin dicarbonyl consumption safety reference.",
          "confidence_tier": "high",
          "notes_on_limits": "Heated soy glycinin model system."
        }
    ]

    existing_papers = {ent["paper_id"] for ent in matrix["entries"]}
    for ent in new_matrix_entries:
        if ent["paper_id"] not in existing_papers:
            matrix["entries"].append(ent)
            print(f"Added {ent['paper_id']} to slr_incorporation_matrix")

    # 4. Add to computational_priors.json (under carbonyl_donor_priors)
    new_priors = [
        {
          "id": "tran_2023_starch_liberation_v1",
          "source": "Tran et al. (2023 / PMC9846454)",
          "provenance_tier": "literature_derived_transfer",
          "confidence_tier": "high",
          "uncertainty_posture": "directional_transfer",
          "applicable_protein_types": ["pea_iso", "plant_protein_generic"],
          "reference_conditions": {
            "process_state": "extrusion_structured",
            "extrusion_mode": "hme"
          },
          "sugar_liberation_umol_g": {
            "pea_starch_maltose": 18.4,
            "waxy_maize_starch_maltose": 31.8,
            "ti_starch_reduction_pct": 66.0
          },
          "notes": "Starch reducing-sugar liberation during HME. Highlights maltose liberation and TI starch reduction."
        }
    ]

    existing_priors = {pr["id"] for pr in priors["carbonyl_donor_priors"]}
    for pr in new_priors:
        if pr["id"] not in existing_priors:
            priors["carbonyl_donor_priors"].append(pr)
            print(f"Added {pr['id']} to computational_priors")

    # 5. Add to flavor_reference_payloads.json (under off_note_reference_anchors)
    new_flavor_payloads = [
        {
          "id": "researchgate_2023_pea_aeda_anchor",
          "compound": "(E,E)-2,4-heptadienal",
          "source_citation": "ResearchGate (2023)",
          "doi": None,
          "matrix_context": "pea_protein_isolate",
          "analytical_method": "AEDA and GC-O",
          "units": "fd_factor",
          "benchmark_role": "reference_anchor",
          "pipeline_role": "reference_only",
          "target_direction": "offnote_dilution_validation",
          "numeric_band_or_point": {
            "type": "point",
            "value": 4096.0
          },
          "notes": "Pea protein AEDA (E,E)-2,4-heptadienal FD 4096."
        }
    ]

    existing_flavors = {fl["id"] for fl in flavor["off_note_reference_anchors"]}
    for fl in new_flavor_payloads:
        if fl["id"] not in existing_flavors:
            flavor["off_note_reference_anchors"].append(fl)
            print(f"Added {fl['id']} to flavor_reference_payloads")

    # 6. Add to safety_reference_payloads.json
    new_safety_payloads = [
        {
          "id": "henle_2005_glycinin_cml",
          "kind": "finished_product_reference",
          "report_visibility": "default",
          "target_module": "safety",
          "chemistry_family": "safety_damage_lane",
          "slr_family_source": "12",
          "source_citation": "Henle (2005)",
          "doi": "10.1007/s00726-005-0187-z",
          "validated_status": "reference_anchor",
          "analyte": "CML",
          "method": {
            "instrument": "GC-MS",
            "notes": "Heated soy glycinin model system CML quantification"
          },
          "matrix_reference_ranges": [
            {
              "matrix_family": "soy_glycinin",
              "units": "umol_per_g_protein",
              "point": 2.3
            }
          ],
          "what_it_supports": [
            "Upper-bound CML yield (2.3 umol/g protein) in heated soy glycinin matrices",
            "Dicarbonyl consumption estimation baseline for soy glycinin"
          ],
          "comment": "Heated soy glycinin model system."
        }
    ]

    existing_safety = {sf["id"] for sf in safety["entries"]}
    for sf in new_safety_payloads:
        if sf["id"] not in existing_safety:
            safety["entries"].append(sf)
            print(f"Added {sf['id']} to safety_reference_payloads")

    # 7. Update deep_research_backlog.json status to RUNTIME_BOUND and set registry_id
    backlog_citations_map = {
        "ResearchGate (2023)": "researchgate_2023_pea_aeda",
        "Malia et al. (2025)": "malia_2025_pea_free_sh_crosscheck",
        "Tran et al. (2023 / PMC9846454)": "tran_2023_reducing_sugar_hme",
        "Henle (2005)": "henle_2005_glycinin_cml"
    }

    updated_backlog_count = 0
    for item in backlog["items"]:
        cit = item.get("citation")
        if cit in backlog_citations_map:
            item["status"] = "RUNTIME_BOUND"
            item["registry_id"] = backlog_citations_map[cit]
            # Since Malia et al has 1 process calibration, ResearchGate has 1 flavor payload,
            # Tran et al has 1 prior, Henle has 1 safety payload:
            item["runtime_artifact_count"] = 1
            updated_backlog_count += 1
            print(f"Updated backlog status for {cit} to RUNTIME_BOUND")

    # Write back all files
    save_json(intake, INTAKE_REGISTRY_PATH)
    save_json(matrix, SLR_MATRIX_PATH)
    save_json(priors, COMPUTATIONAL_PRIORS_PATH)
    save_json(flavor, FLAVOR_PAYLOADS_PATH)
    save_json(safety, SAFETY_PAYLOADS_PATH)
    save_json(backlog, BACKLOG_PATH)

    print("Successfully ingested Families 06 and 07 references!")


if __name__ == "__main__":
    ingest_family06_07()
