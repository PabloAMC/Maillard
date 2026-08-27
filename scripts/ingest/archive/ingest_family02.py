import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_LIT_DIR = ROOT / "data" / "lit"

# File paths
INTAKE_REGISTRY_PATH = DATA_LIT_DIR / "benchmark_intake_registry.json"
SLR_MATRIX_PATH = DATA_LIT_DIR / "slr_incorporation_matrix.json"
COMPUTATIONAL_PRIORS_PATH = DATA_LIT_DIR / "computational_priors.json"
FLAVOR_PAYLOADS_PATH = DATA_LIT_DIR / "flavor_reference_payloads.json"
BACKLOG_PATH = DATA_LIT_DIR / "deep_research_backlog.json"


def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def save_json(data, path):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"Updated {path.name}")


def ingest_family02():
    # 1. Load all files
    intake = load_json(INTAKE_REGISTRY_PATH)
    matrix = load_json(SLR_MATRIX_PATH)
    priors = load_json(COMPUTATIONAL_PRIORS_PATH)
    flavor = load_json(FLAVOR_PAYLOADS_PATH)
    backlog = load_json(BACKLOG_PATH)

    # Define the 5 new eligible references
    new_eligible_references = [
        {
          "id": "farmer_1991_alkyl_thiazoles",
          "citation": "Farmer & Patterson (1991)",
          "doi": "10.1016/0308-8146(91)90013-I",
          "kind": "calibration_reference",
          "chemistry_family": "lipid_oxidation_and_carbonylic_crosstalk",
          "slr_family_source": "02",
          "payload_role": "benchmark_intake",
          "observable_panel_tags": [
            "lipid",
            "thiazole"
          ],
          "process_state_scope": [
            "heated_matrix"
          ],
          "supporting_families": [
            "01",
            "08"
          ],
          "target_modules": [
            "literature_runtime",
            "lipid_condensation"
          ],
          "matrix_family": "free_model_system",
          "eligibility": "benchmark_eligible",
          "status": "ready_for_intake_encoding",
          "what_it_supports": [
            "Aliphatic aldehyde + H2S + NH3 thiazole cyclisation kinetics and relative yields.",
            "Confirms alkyl thiazole OAV > 1 in lipid-crosstalk modeling."
          ],
          "what_it_does_not_support": [
            "Protein-matrix retention constants",
            "pH-dependent oxidation rate limits"
          ],
          "key_values": {
            "alkyl_thiazole_oav_ceiling": 50.0,
            "aldehyde_precursor_reactivity_ratio": 1.2
          },
          "repo_next_action": "Expose as a lipid-Maillard thiazole cross-coupling mechanistic reference.",
          "runtime_artifacts": [
            {
              "artifact_type": "intake_registry_entry",
              "artifact_id": "farmer_1991_alkyl_thiazoles",
              "path": "data/lit/benchmark_intake_registry.json"
            },
            {
              "artifact_type": "slr_incorporation_ledger",
              "artifact_id": "farmer_1991_alkyl_thiazoles",
              "path": "data/lit/slr_incorporation_matrix.json"
            }
          ],
          "requires_primary_data": False
        },
        {
          "id": "esterbauer_1991_4hne_kinetics",
          "citation": "Esterbauer et al. (1991)",
          "doi": "10.1016/0891-5849(91)90081-A",
          "kind": "calibration_reference",
          "chemistry_family": "lipid_oxidation_and_carbonylic_crosstalk",
          "slr_family_source": "02",
          "payload_role": "benchmark_intake",
          "observable_panel_tags": [
            "lipid",
            "4hne",
            "scavenging"
          ],
          "process_state_scope": [
            "heated_matrix"
          ],
          "supporting_families": [
            "01",
            "08"
          ],
          "target_modules": [
            "literature_runtime",
            "lipid_condensation"
          ],
          "matrix_family": "free_model_system",
          "eligibility": "benchmark_eligible",
          "status": "ready_for_intake_encoding",
          "what_it_supports": [
            "Unsaturated aldehyde (4-HNE) scavenging rate hierarchy: Cys >> His > Lys.",
            "Michael addition rate constants for reactive amino acid sidechains."
          ],
          "what_it_does_not_support": [
            "Saturated aldehyde trapping kinetics",
            "Denaturation-induced accessibility values"
          ],
          "key_values": {
            "cys_Michael_addition_rate_M_1_s_1": 1.33,
            "his_Michael_addition_rate_M_1_s_1": 0.001,
            "lys_Schiff_addition_rate_M_1_s_1": 0.0001
          },
          "repo_next_action": "Encode 4-HNE reactive amino acid scavenging rate priors.",
          "runtime_artifacts": [
            {
              "artifact_type": "computational_prior",
              "artifact_id": "esterbauer_1991_4hne_kinetics_v1",
              "path": "data/lit/computational_priors.json",
              "sections": [
                "lipid_offnote_priors"
              ]
            }
          ],
          "requires_primary_data": False
        },
        {
          "id": "kamal_eldin_2003_triolein_scission",
          "citation": "Kamal-Eldin et al. (2003)",
          "doi": "10.1007/s11745-003-1082-x",
          "kind": "calibration_reference",
          "chemistry_family": "lipid_oxidation_and_carbonylic_crosstalk",
          "slr_family_source": "02",
          "payload_role": "benchmark_intake",
          "observable_panel_tags": [
            "lipid",
            "scission",
            "triolein"
          ],
          "process_state_scope": [
            "heated_matrix"
          ],
          "supporting_families": [
            "11"
          ],
          "target_modules": [
            "literature_runtime",
            "lipid_oxidation"
          ],
          "matrix_family": "free_model_system",
          "eligibility": "benchmark_eligible",
          "status": "ready_for_intake_encoding",
          "what_it_supports": [
            "Triolein oxidation and oleic acid (C18:1) beta-scission yield parameters."
          ],
          "what_it_does_not_support": [
            "Phospholipid fraction contributions",
            "LOX inactivation temperatures"
          ],
          "key_values": {
            "c181_to_octanal_beta_scission_yield_pct": 3.5,
            "c181_to_nonanal_beta_scission_yield_pct": 4.2
          },
          "repo_next_action": "Encode oleic (C18:1) scission volatile split priors.",
          "runtime_artifacts": [
            {
              "artifact_type": "computational_prior",
              "artifact_id": "kamal_eldin_2003_triolein_scission_v1",
              "path": "data/lit/computational_priors.json",
              "sections": [
                "lipid_offnote_priors"
              ]
            }
          ],
          "requires_primary_data": False
        },
        {
          "id": "xu_2024_soybean_pbma_hexanal",
          "citation": "Xu et al. (2024)",
          "doi": "10.3390/foods13091393",
          "kind": "calibration_reference",
          "chemistry_family": "lipid_oxidation_and_carbonylic_crosstalk",
          "slr_family_source": "02",
          "payload_role": "benchmark_intake",
          "observable_panel_tags": [
            "lipid",
            "hexanal",
            "soy_isolate"
          ],
          "process_state_scope": [
            "ambient_slurry"
          ],
          "supporting_families": [
            "06",
            "08"
          ],
          "target_modules": [
            "literature_runtime",
            "headspace"
          ],
          "matrix_family": "soy_isolate",
          "eligibility": "benchmark_eligible",
          "status": "ready_for_intake_encoding",
          "what_it_supports": [
            "Soybean PBMA storage hexanal accumulation trends and baseline concentration values."
          ],
          "what_it_does_not_support": [
            "Free solution kinetic activation energies",
            "Furanone reaction branching ratios"
          ],
          "key_values": {
            "soy_pbma_storage_hexanal_ppb": 185.0,
            "storage_temp_C": 4.0,
            "storage_time_days": 14.0
          },
          "repo_next_action": "Expose as a soybean PBMA storage flavor reference payload.",
          "runtime_artifacts": [
            {
              "artifact_type": "flavor_reference_payload",
              "artifact_id": "xu_2024_soybean_pbma_hexanal_anchor",
              "path": "data/lit/flavor_reference_payloads.json"
            }
          ],
          "requires_primary_data": False
        },
        {
          "id": "yeo_shibamoto_1991_wof_hexanal",
          "citation": "Yeo & Shibamoto (1991)",
          "doi": "10.1021/jf00002a024",
          "kind": "calibration_reference",
          "chemistry_family": "lipid_oxidation_and_carbonylic_crosstalk",
          "slr_family_source": "02",
          "payload_role": "benchmark_intake",
          "observable_panel_tags": [
            "lipid",
            "hexanal",
            "lox"
          ],
          "process_state_scope": [
            "heated_matrix"
          ],
          "supporting_families": [
            "06",
            "08"
          ],
          "target_modules": [
            "literature_runtime",
            "lipid_oxidation"
          ],
          "matrix_family": "free_model_system",
          "eligibility": "benchmark_eligible",
          "status": "ready_for_intake_encoding",
          "what_it_supports": [
            "Phospholipid fraction contribution (82-88%) to warmed-over flavor (WOF) hexanal.",
            "LOX denaturation D-value table (70 C D=28 min, z~10 C)."
          ],
          "what_it_does_not_support": [
            "Direct peptide-cysteine nucleophilicity indices",
            "Absolute volatile release in non-meat systems"
          ],
          "key_values": {
            "phospholipid_hexanal_fraction_pct_range": [
              82.0,
              88.0
            ],
            "lox_d_value_70C_min": 28.0,
            "lox_z_value_C": 10.0
          },
          "repo_next_action": "Encode phospholipid hexanal contribution and LOX denaturation priors.",
          "runtime_artifacts": [
            {
              "artifact_type": "computational_prior",
              "artifact_id": "yeo_shibamoto_1991_wof_hexanal_v1",
              "path": "data/lit/computational_priors.json",
              "sections": [
                "lipid_offnote_priors"
              ]
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
          "slr_section": "2.1",
          "paper_id": "farmer_1991_alkyl_thiazoles",
          "citation": "Farmer & Patterson (1991), Food Chem",
          "doi": "10.1016/0308-8146(91)90013-I",
          "matrix_family": "free_model_system",
          "compounds_supported": ["alkyl thiazoles"],
          "parameters_supported": ["thiazole_cyclisation_selectivity"],
          "exact_numeric_anchors": ["thiazole OAV > 1"],
          "current_repo_artifacts": ["data/lit/benchmark_intake_registry.json"],
          "current_runtime_consumers": ["src/lipid_oxidation.py"],
          "current_user_visible_surfaces": ["reporting"],
          "incorporation_status": "encoded_shown",
          "next_action": "Expose as lipid-Maillard thiazole cross-coupling mechanistic reference.",
          "confidence_tier": "high",
          "notes_on_limits": "Free lipid-amine model system."
        },
        {
          "slr_section": "2.1",
          "paper_id": "esterbauer_1991_4hne_kinetics",
          "citation": "Esterbauer et al. (1991), Free Radic Biol Med",
          "doi": "10.1016/0891-5849(91)90081-A",
          "matrix_family": "free_model_system",
          "compounds_supported": ["4-HNE"],
          "parameters_supported": ["Michael_addition_reactivity"],
          "exact_numeric_anchors": ["Cys >> His > Lys"],
          "current_repo_artifacts": ["data/lit/computational_priors.json"],
          "current_runtime_consumers": ["src/lipid_oxidation.py"],
          "current_user_visible_surfaces": ["reporting"],
          "incorporation_status": "encoded_shown",
          "next_action": "Encode 4-HNE reactive amino acid scavenging rate priors.",
          "confidence_tier": "high",
          "notes_on_limits": "Free lipid-aldehyde-amino-acid kinetics."
        },
        {
          "slr_section": "2.1",
          "paper_id": "kamal_eldin_2003_triolein_scission",
          "citation": "Kamal-Eldin et al. (2003), Lipids",
          "doi": "10.1007/s11745-003-1082-x",
          "matrix_family": "free_model_system",
          "compounds_supported": ["octanal", "nonanal"],
          "parameters_supported": ["c181_scission_yield"],
          "exact_numeric_anchors": ["octanal 3.5%", "nonanal 4.2%"],
          "current_repo_artifacts": ["data/lit/computational_priors.json"],
          "current_runtime_consumers": ["src/lipid_oxidation.py"],
          "current_user_visible_surfaces": ["reporting"],
          "incorporation_status": "encoded_shown",
          "next_action": "Encode oleic (C18:1) scission volatile split priors.",
          "confidence_tier": "high",
          "notes_on_limits": "Free lipid model system."
        },
        {
          "slr_section": "2.2",
          "paper_id": "xu_2024_soybean_pbma_hexanal",
          "citation": "Xu et al. (2024), Foods",
          "doi": "10.3390/foods13091393",
          "matrix_family": "soy_isolate",
          "compounds_supported": ["hexanal"],
          "parameters_supported": ["soy_pbma_storage_hexanal"],
          "exact_numeric_anchors": ["hexanal 185 ppb"],
          "current_repo_artifacts": ["data/lit/flavor_reference_payloads.json"],
          "current_runtime_consumers": ["src/sensory_scoring.py"],
          "current_user_visible_surfaces": ["sensory_reporting"],
          "incorporation_status": "encoded_shown",
          "next_action": "Expose as a soybean PBMA storage flavor reference payload.",
          "confidence_tier": "high",
          "notes_on_limits": "Soy plant-based meat alternative matrix."
        },
        {
          "slr_section": "2.2",
          "paper_id": "yeo_shibamoto_1991_wof_hexanal",
          "citation": "Yeo & Shibamoto (1991), JAFC",
          "doi": "10.1021/jf00002a024",
          "matrix_family": "free_model_system",
          "compounds_supported": ["hexanal"],
          "parameters_supported": ["phospholipid_hexanal_pct", "lox_denaturation_d_value"],
          "exact_numeric_anchors": ["phospholipid fraction 82-88%", "LOX D-value 28 min"],
          "current_repo_artifacts": ["data/lit/computational_priors.json"],
          "current_runtime_consumers": ["src/lipid_oxidation.py"],
          "current_user_visible_surfaces": ["reporting"],
          "incorporation_status": "encoded_shown",
          "next_action": "Encode phospholipid hexanal contribution and LOX denaturation priors.",
          "confidence_tier": "high",
          "notes_on_limits": "LOX denaturation kinetics."
        }
    ]

    existing_papers = {ent["paper_id"] for ent in matrix["entries"]}
    for ent in new_matrix_entries:
        if ent["paper_id"] not in existing_papers:
            matrix["entries"].append(ent)
            print(f"Added {ent['paper_id']} to slr_incorporation_matrix")

    # C. Add to computational_priors.json (under lipid_offnote_priors)
    new_priors = [
        {
          "id": "esterbauer_1991_4hne_kinetics_v1",
          "source": "Esterbauer et al. (1991)",
          "provenance_tier": "literature_derived_transfer",
          "confidence_tier": "high",
          "uncertainty_posture": "directional_transfer",
          "applicable_protein_types": ["free"],
          "reference_conditions": {"ph": 7.4, "temperature_c": 37.0},
          "rate_constants": {
            "cys_Michael_addition_rate_M_1_s_1": 1.33,
            "his_Michael_addition_rate_M_1_s_1": 0.001,
            "lys_Schiff_addition_rate_M_1_s_1": 0.0001
          },
          "notes": "Unsaturated aldehyde (4-HNE) scavenging rate constants."
        },
        {
          "id": "kamal_eldin_2003_triolein_scission_v1",
          "source": "Kamal-Eldin et al. (2003)",
          "provenance_tier": "literature_derived_transfer",
          "confidence_tier": "high",
          "uncertainty_posture": "directional_transfer",
          "applicable_protein_types": ["free"],
          "reference_conditions": {},
          "yields": {
            "c181_to_octanal_beta_scission_yield_pct": 3.5,
            "c181_to_nonanal_beta_scission_yield_pct": 4.2
          },
          "notes": "Oleic acid (C18:1) beta-scission yield splits."
        },
        {
          "id": "yeo_shibamoto_1991_wof_hexanal_v1",
          "source": "Yeo & Shibamoto (1991)",
          "provenance_tier": "literature_derived_transfer",
          "confidence_tier": "high",
          "uncertainty_posture": "directional_transfer",
          "applicable_protein_types": ["free"],
          "reference_conditions": {},
          "yields": {
            "phospholipid_hexanal_fraction_pct_range": [82.0, 88.0],
            "lox_d_value_70C_min": 28.0,
            "lox_z_value_C": 10.0
          },
          "notes": "Phospholipid contribution to hexanal and LOX denaturation kinetics."
        }
    ]

    existing_priors = {pr["id"] for pr in priors["lipid_offnote_priors"]}
    for pr in new_priors:
        if pr["id"] not in existing_priors:
            priors["lipid_offnote_priors"].append(pr)
            print(f"Added {pr['id']} to computational_priors")

    # D. Add to flavor_reference_payloads.json
    new_flavor_payload = {
      "id": "xu_2024_soybean_pbma_hexanal_anchor",
      "compound": "hexanal",
      "source_citation": "Xu et al. (2024)",
      "doi": "10.3390/foods13091393",
      "matrix_context": "soy_isolate",
      "analytical_method": "GC-MS with internal standard",
      "units": "ppb",
      "benchmark_role": "reference_anchor",
      "pipeline_role": "reference_only",
      "target_direction": "hexanal_storage_accumulation",
      "numeric_band_or_point": {
        "type": "point",
        "value": 185.0
      },
      "notes": "Soybean PBMA storage hexanal accumulation baseline of 185 ppb at 4 C after 14 days."
    }

    existing_flavors = {fl["id"] for fl in flavor["off_note_reference_anchors"]}
    if new_flavor_payload["id"] not in existing_flavors:
        flavor["off_note_reference_anchors"].append(new_flavor_payload)
        print(f"Added {new_flavor_payload['id']} to flavor_reference_payloads")

    # E. Update deep_research_backlog.json status to RUNTIME_BOUND and set registry_id
    backlog_citations_map = {
        "Farmer & Patterson (1991)": "farmer_1991_alkyl_thiazoles",
        "Esterbauer et al. (1991)": "esterbauer_1991_4hne_kinetics",
        "Kamal-Eldin et al. (2003)": "kamal_eldin_2003_triolein_scission",
        "Xu et al. (2024)": "xu_2024_soybean_pbma_hexanal",
        "Yeo & Shibamoto (1991)": "yeo_shibamoto_1991_wof_hexanal"
    }

    updated_backlog_count = 0
    for item in backlog["items"]:
        cit = item.get("citation")
        if cit in backlog_citations_map:
            item["status"] = "RUNTIME_BOUND"
            item["registry_id"] = backlog_citations_map[cit]
            item["runtime_artifact_count"] = 1
            updated_backlog_count += 1
            print(f"Updated backlog status for {cit} to RUNTIME_BOUND")

    # Write back all files
    save_json(intake, INTAKE_REGISTRY_PATH)
    save_json(matrix, SLR_MATRIX_PATH)
    save_json(priors, COMPUTATIONAL_PRIORS_PATH)
    save_json(flavor, FLAVOR_PAYLOADS_PATH)
    save_json(backlog, BACKLOG_PATH)
    
    print(f"Successfully ingested all 5 Family 02 references!")


if __name__ == "__main__":
    ingest_family02()
