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
BACKLOG_PATH = DATA_LIT_DIR / "deep_research_backlog.json"


def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def save_json(data, path):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"Updated {path.name}")


def ingest_chunk_2():
    # 1. Load all files
    intake = load_json(INTAKE_REGISTRY_PATH)
    matrix = load_json(SLR_MATRIX_PATH)
    priors = load_json(COMPUTATIONAL_PRIORS_PATH)
    flavor = load_json(FLAVOR_PAYLOADS_PATH)
    process_state = load_json(PROCESS_STATE_CALIBRATIONS_PATH)
    backlog = load_json(BACKLOG_PATH)

    # Define new eligible references for benchmark_intake_registry.json
    new_eligible_references = [
        {
            "id": "liu_2022_ppi_oav_anchors",
            "citation": "Liu (2022), NC State Repository",
            "doi": "",
            "kind": "calibration_reference",
            "chemistry_family": "off_note_and_maillard_suppression",
            "slr_family_source": "08",
            "payload_role": "flavor_reference_payload",
            "observable_panel_tags": [
                "off_note",
                "methoxypyrazine",
                "threshold"
            ],
            "process_state_scope": [
                "heated_matrix"
            ],
            "supporting_families": [
                "02"
            ],
            "target_modules": [
                "literature_runtime",
                "safety"
            ],
            "matrix_family": "pea_isolate",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Baseline OAV table for 12 key pea protein isolate (PPI) volatiles",
                "Methoxypyrazines odor detection threshold (ODT) of 0.002-0.006 ppb"
            ],
            "what_it_does_not_support": [
                "Free solution kinetic activation energies",
                "Browning index evolution rates"
            ],
            "key_values": {
                "methoxypyrazine_odt_ppb_range": [0.002, 0.006],
                "volatiles_count": 12
            },
            "repo_next_action": "Expose as a pea isolate off-note flavor reference payload.",
            "runtime_artifacts": [
                {
                    "artifact_type": "flavor_reference_payload",
                    "artifact_id": "liu_2022_ppi_oav_anchors",
                    "path": "data/lit/flavor_reference_payloads.json"
                }
            ],
            "requires_primary_data": False
        },
        {
            "id": "rawel_2002_cga_cysteine_blocking",
            "citation": "Rawel et al. (2002), JAFC 50:5343-5350",
            "doi": "10.1021/jf020082z",
            "kind": "calibration_reference",
            "chemistry_family": "off_note_and_maillard_suppression",
            "slr_family_source": "08",
            "payload_role": "computational_prior",
            "observable_panel_tags": [
                "guardrail",
                "cysteine_sink",
                "polyphenol_capping"
            ],
            "process_state_scope": [
                "heated_matrix"
            ],
            "supporting_families": [
                "05",
                "13"
            ],
            "target_modules": [
                "literature_runtime"
            ],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Chlorogenic acid (CGA) Cys blocking k_obs is 2.9x faster than Lys",
                "Lysine blocking efficiency of 43% at a 10:1 polyphenol-to-protein ratio"
            ],
            "what_it_does_not_support": [
                "Matrix-assisted cysteine release rates",
                "Radical scavenging mechanism validation"
            ],
            "key_values": {
                "cys_blocking_rate_vs_lys": 2.9,
                "lys_blocked_pct_at_10_1": 43.0
            },
            "repo_next_action": "Encode CGA Cys and Lys blocking rate priors.",
            "runtime_artifacts": [
                {
                    "artifact_type": "computational_prior",
                    "artifact_id": "rawel_2002_cga_cysteine_blocking",
                    "path": "data/lit/computational_priors.json"
                }
            ],
            "requires_primary_data": False
        },
        {
            "id": "quintas_2000_sucrose_caramelisation",
            "citation": "Quintas et al. (2000), J Food Engineering",
            "doi": "10.1016/S0260-8774(00)00047-9",
            "kind": "calibration_reference",
            "chemistry_family": "carbohydrate_pyrolysis_and_caramelization",
            "slr_family_source": "09",
            "payload_role": "computational_prior",
            "observable_panel_tags": [
                "caramelization",
                "hmf",
                "kinetics"
            ],
            "process_state_scope": [
                "heated_matrix"
            ],
            "supporting_families": [
                "07"
            ],
            "target_modules": [
                "literature_runtime"
            ],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Sucrose caramelisation k_obs and Ea versus pH and Temperature",
                "HMF generation of 48 mg/L at pH 5.0, 120 C after 60 min"
            ],
            "what_it_does_not_support": [
                "Free amino acid condensation mechanisms",
                "Extruded texture impact parameters"
            ],
            "key_values": {
                "hmf_mg_per_l_at_ph5_120c_60min": 48.0
            },
            "repo_next_action": "Encode sucrose caramelisation kinetics priors.",
            "runtime_artifacts": [
                {
                    "artifact_type": "computational_prior",
                    "artifact_id": "quintas_2000_sucrose_caramelisation",
                    "path": "data/lit/computational_priors.json"
                }
            ],
            "requires_primary_data": False
        },
        {
            "id": "hauck_tressl_1999_hdmf_non_amino",
            "citation": "Hauck & Tressl (1999), ACS review",
            "doi": "10.1021/bk-1999-0740.ch012",
            "kind": "calibration_reference",
            "chemistry_family": "carbohydrate_pyrolysis_and_caramelization",
            "slr_family_source": "09",
            "payload_role": "flavor_reference_payload",
            "observable_panel_tags": [
                "furanone",
                "hdmf",
                "caramelization"
            ],
            "process_state_scope": [
                "heated_matrix"
            ],
            "supporting_families": [
                "01"
            ],
            "target_modules": [
                "literature_runtime"
            ],
            "matrix_family": "free_model_system",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "HDMF generation from iso-oligosaccharides (3.8 ug/g) and maltose (1.4 ug/g) at 140 C without amino acids",
                "Methylglyoxal (MGO) self-condensation HDMF route yield (2.9 ug/g)"
            ],
            "what_it_does_not_support": [
                "Amadori rearrangement activation barriers",
                "Peptide-bound lysine trapping kinetics"
            ],
            "key_values": {
                "hdmf_from_iso_oligosaccharides_ug_per_g": 3.8,
                "hdmf_from_maltose_ug_per_g": 1.4,
                "hdmf_from_mgo_self_condensation_ug_per_g": 2.9,
                "reaction_temp_C": 140.0
            },
            "repo_next_action": "Expose as a non-amino acid dependent furanone flavor reference payload.",
            "runtime_artifacts": [
                {
                    "artifact_type": "flavor_reference_payload",
                    "artifact_id": "hauck_tressl_1999_hdmf_non_amino",
                    "path": "data/lit/flavor_reference_payloads.json"
                }
            ],
            "requires_primary_data": False
        },
        {
            "id": "perdiguero_2004_yeast_autolysis_nucleotides",
            "citation": "Perdiguero et al. (2004), JAFC 52:7802-7807",
            "doi": "10.1021/jf0494452",
            "kind": "calibration_reference",
            "chemistry_family": "fermentation_pretreatment",
            "slr_family_source": "10",
            "payload_role": "process_state_calibration",
            "observable_panel_tags": [
                "fermentation",
                "nucleotide",
                "yeast"
            ],
            "process_state_scope": [
                "heated_matrix"
            ],
            "supporting_families": [
                "04"
            ],
            "target_modules": [
                "literature_runtime",
                "recommend"
            ],
            "matrix_family": "yeast_ye",
            "eligibility": "benchmark_eligible",
            "status": "ready_for_intake_encoding",
            "what_it_supports": [
                "Yeast autolysis nucleotide release kinetics yielding GMP 2.8 mg/g DW at 55 C after 24h"
            ],
            "what_it_does_not_support": [
                "Free amino acid extraction scaling parameters",
                "Browning intensity at high temperature"
            ],
            "key_values": {
                "gmp_mg_per_g_dw": 2.8,
                "temp_C": 55.0,
                "time_h": 24.0
            },
            "repo_next_action": "Expose as a fermentation process state calibration payload.",
            "runtime_artifacts": [
                {
                    "artifact_type": "process_state_calibration",
                    "artifact_id": "perdiguero_2004_yeast_autolysis_nucleotides",
                    "path": "data/lit/process_state_calibrations.json"
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
            "slr_section": "8",
            "paper_id": "liu_2022_ppi_oav_anchors",
            "citation": "Liu (2022), NC State Repository",
            "doi": "",
            "matrix_family": "pea_isolate",
            "compounds_supported": ["methoxypyrazines"],
            "parameters_supported": ["odor_detection_thresholds"],
            "exact_numeric_anchors": ["methoxypyrazines ODT 0.002–0.006 ppb"],
            "current_repo_artifacts": ["data/lit/flavor_reference_payloads.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Expose as a pea isolate off-note flavor reference payload.",
            "confidence_tier": "high",
            "notes_on_limits": "Pea protein isolate matrix."
        },
        {
            "slr_section": "8.2",
            "paper_id": "rawel_2002_cga_cysteine_blocking",
            "citation": "Rawel et al. (2002), JAFC",
            "doi": "10.1021/jf020082z",
            "matrix_family": "free_model_system",
            "compounds_supported": ["cysteine", "lysine"],
            "parameters_supported": ["cysteine_depletion_kinetics"],
            "exact_numeric_anchors": ["Cys blocking rate 2.9x faster than Lys", "43% Lys blocked at 10:1 ratio"],
            "current_repo_artifacts": ["data/lit/computational_priors.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Encode CGA Cys and Lys blocking rate priors.",
            "confidence_tier": "high",
            "notes_on_limits": "Polyphenol capping kinetics."
        },
        {
            "slr_section": "9.1",
            "paper_id": "quintas_2000_sucrose_caramelisation",
            "citation": "Quintas et al. (2000), J Food Eng",
            "doi": "10.1016/S0260-8774(00)00047-9",
            "matrix_family": "free_model_system",
            "compounds_supported": ["HMF"],
            "parameters_supported": ["sucrose_caramelisation_kinetics"],
            "exact_numeric_anchors": ["HMF 48 mg/L at pH 5.0, 120 C, 60 min"],
            "current_repo_artifacts": ["data/lit/computational_priors.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Encode sucrose caramelisation kinetics priors.",
            "confidence_tier": "high",
            "notes_on_limits": "Sucrose caramelisation system."
        },
        {
            "slr_section": "9.1",
            "paper_id": "hauck_tressl_1999_hdmf_non_amino",
            "citation": "Hauck & Tressl (1999), ACS review",
            "doi": "10.1021/bk-1999-0740.ch012",
            "matrix_family": "free_model_system",
            "compounds_supported": ["HDMF"],
            "parameters_supported": ["furanone_generation_without_aa"],
            "exact_numeric_anchors": ["HDMF from oligosaccharides 3.8 ug/g", "HDMF from maltose 1.4 ug/g"],
            "current_repo_artifacts": ["data/lit/flavor_reference_payloads.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Expose as a non-amino acid dependent furanone flavor reference payload.",
            "confidence_tier": "high",
            "notes_on_limits": "Carbohydrate pyrolysis system."
        },
        {
            "slr_section": "10",
            "paper_id": "perdiguero_2004_yeast_autolysis_nucleotides",
            "citation": "Perdiguero et al. (2004), JAFC",
            "doi": "10.1021/jf0494452",
            "matrix_family": "yeast_ye",
            "compounds_supported": ["GMP"],
            "parameters_supported": ["yeast_autolysis_kinetics"],
            "exact_numeric_anchors": ["GMP 2.8 mg/g DW at 55 C, 24h"],
            "current_repo_artifacts": ["data/lit/process_state_calibrations.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Expose as a fermentation process state calibration payload.",
            "confidence_tier": "high",
            "notes_on_limits": "Yeast cell autolysis matrix."
        },
        {
            "slr_section": "11",
            "paper_id": "acs_2020_raw_pea_hexanal_baseline",
            "citation": "Bi et al. (2020), JAFC",
            "doi": "10.1021/acs.jafc.9b07711",
            "matrix_family": "raw_pea_flour",
            "compounds_supported": ["hexanal"],
            "parameters_supported": ["hexanal_baseline_concentration"],
            "exact_numeric_anchors": ["hexanal 1260 ug/kg"],
            "current_repo_artifacts": ["data/lit/benchmark_intake_registry.json", "data/lit/flavor_reference_payloads.json"],
            "current_runtime_consumers": ["src/literature_runtime.py"],
            "current_user_visible_surfaces": ["reporting"],
            "incorporation_status": "encoded_shown",
            "next_action": "Track raw pea off-note baseline in lipid-Maillard crosstalk model.",
            "confidence_tier": "high",
            "notes_on_limits": "Raw pea flour matrix."
        }
    ]

    existing_papers = {ent["paper_id"] for ent in matrix["entries"]}
    for ent in new_matrix_entries:
        if ent["paper_id"] not in existing_papers:
            matrix["entries"].append(ent)
            print(f"Added {ent['paper_id']} to slr_incorporation_matrix")

    # C. Add to computational_priors.json
    new_polyphenol_priors = [
        {
            "id": "rawel_2002_cga_cysteine_blocking",
            "source": "Rawel et al. (2002)",
            "provenance_tier": "literature_derived_transfer",
            "confidence_tier": "high",
            "uncertainty_posture": "directional_transfer",
            "cys_blocking_rate_vs_lys": 2.9,
            "lys_blocked_percent_at_10_1_ratio": 43.0,
            "effect_direction": "increase_cysteine_sink_and_reduce_free_thiol_budget",
            "notes": "Chlorogenic acid Cys blocking kinetics; Cys blocking rate is 2.9x faster than Lys, 43% Lys blocked at 10:1 ratio."
        }
    ]
    existing_polyphenol_priors = {pr["id"] for pr in priors["polyphenol_thiol_capping_priors"]}
    for pr in new_polyphenol_priors:
        if pr["id"] not in existing_polyphenol_priors:
            priors["polyphenol_thiol_capping_priors"].append(pr)
            print(f"Added {pr['id']} to computational_priors (polyphenol_thiol_capping_priors)")

    new_carbonyl_priors = [
        {
            "id": "quintas_2000_sucrose_caramelisation",
            "source": "Quintas et al. (2000)",
            "provenance_tier": "literature_derived_transfer",
            "confidence_tier": "high",
            "uncertainty_posture": "directional_transfer",
            "reference_conditions": {
                "ph": 5.0,
                "temperature_c": 120.0,
                "time_min": 60.0
            },
            "sucrose_caramelisation": {
                "hmf_mg_per_l_at_ph5_120c_60min": 48.0
            },
            "notes": "Sucrose caramelisation kinetics showing k_obs and Ea versus pH and Temperature. Quantifies HMF 48 mg/L at pH 5.0, 120°C, 60 min."
        }
    ]
    existing_carbonyl_priors = {pr["id"] for pr in priors["carbonyl_donor_priors"]}
    for pr in new_carbonyl_priors:
        if pr["id"] not in existing_carbonyl_priors:
            priors["carbonyl_donor_priors"].append(pr)
            print(f"Added {pr['id']} to computational_priors (carbonyl_donor_priors)")

    # D. Add to flavor_reference_payloads.json
    new_off_note_payload = {
        "id": "liu_2022_ppi_oav_anchors",
        "compound": "methoxypyrazines",
        "source_citation": "Liu (2022)",
        "doi": "",
        "matrix_context": "pea_protein_isolate",
        "analytical_method": "AEDA/GC-MS/O",
        "units": "ppb",
        "benchmark_role": "reference_anchor",
        "pipeline_role": "reference_only",
        "target_direction": "track_pea_off_note_volatiles_and_thresholds",
        "numeric_band_or_point": {
            "type": "band",
            "lower_bound": 0.002,
            "upper_bound": 0.006
        },
        "notes": "Baseline OAV table for 12 key pea protein isolate (PPI) volatiles and methoxypyrazines ODT (0.002-0.006 ppb)."
    }
    existing_off_notes = {fl["id"] for fl in flavor["off_note_reference_anchors"]}
    if new_off_note_payload["id"] not in existing_off_notes:
        flavor["off_note_reference_anchors"].append(new_off_note_payload)
        print(f"Added {new_off_note_payload['id']} to flavor_reference_payloads (off_note_reference_anchors)")

    new_furanone_payload = {
        "id": "hauck_tressl_1999_hdmf_non_amino",
        "compound": "HDMF",
        "source_citation": "Hauck & Tressl (1999)",
        "doi": "10.1021/bk-1999-0740.ch012",
        "matrix_context": "sugar_pyrolysis_without_amino_acids",
        "analytical_method": "SIDA",
        "units": "ug_per_g",
        "benchmark_role": "reference_anchor",
        "pipeline_role": "reference_only",
        "target_direction": "non_amino_acid_dependent_furanone_yield",
        "numeric_band_or_point": {
            "type": "point",
            "value": 3.8
        },
        "key_values": {
            "hdmf_from_iso_oligosaccharides_ug_per_g": 3.8,
            "hdmf_from_maltose_ug_per_g": 1.4,
            "hdmf_from_mgo_self_condensation_ug_per_g": 2.9,
            "reaction_temp_C": 140.0
        },
        "notes": "Quantifies HDMF from iso-oligosaccharides (3.8 ug/g) and maltose (1.4 ug/g) at 140 C without amino acids."
    }
    existing_furanones = {fl["id"] for fl in flavor["furanone_reference_anchors"]}
    if new_furanone_payload["id"] not in existing_furanones:
        flavor["furanone_reference_anchors"].append(new_furanone_payload)
        print(f"Added {new_furanone_payload['id']} to flavor_reference_payloads (furanone_reference_anchors)")

    # E. Add to process_state_calibrations.json
    new_process_state_entry = {
        "id": "perdiguero_2004_yeast_autolysis_nucleotides",
        "kind": "fermentation_release_context",
        "protein_type": "yeast_ye",
        "source_citation": "Perdiguero et al. (2004), J. Agric. Food Chem. 52:7802-7807",
        "doi": "10.1021/jf0494452",
        "validated_status": "parameter_anchor",
        "provenance_tier": "direct_measurement",
        "conditions": {
            "matrix_label": "yeast autolysis",
            "temp_C": 55.0,
            "time_h": 24.0,
            "ph": 5.0
        },
        "numeric_anchors": {
            "assay": "HPLC-UV",
            "gmp_mg_per_g_dw": 2.8,
            "is_amp_d4_used": True,
            "replicates": 3
        },
        "what_it_supports": [
            "Yeast autolysis nucleotide release kinetics yielding GMP 2.8 mg/g DW at 55 C after 24h",
            "Fermentation pretreatment node calibrations in literature_runtime.py"
        ],
        "what_it_does_not_support": [
            "Cysteine or other free amino acid release profiles under same conditions",
            "Alternative enzyme cocktail efficiency metrics"
        ],
        "comment": "Key mechanistic parameter anchor for yeast autolysis nucleotides under Family 10."
    }
    existing_process_states = {ent["id"] for ent in process_state["entries"]}
    if new_process_state_entry["id"] not in existing_process_states:
        process_state["entries"].append(new_process_state_entry)
        print(f"Added {new_process_state_entry['id']} to process_state_calibrations")

    # F. Update deep_research_backlog.json status to RUNTIME_BOUND and set registry_id
    backlog_citations_map = {
        "Liu (2022)": "liu_2022_ppi_oav_anchors",
        "Rawel et al. (2002)": "rawel_2002_cga_cysteine_blocking",
        "Quintas et al. (2000)": "quintas_2000_sucrose_caramelisation",
        "Hauck & Tressl (1999)": "hauck_tressl_1999_hdmf_non_amino",
        "Perdiguero et al. (2004)": "perdiguero_2004_yeast_autolysis_nucleotides",
        "ACS JAFC 2020 (raw peas)": "acs_2020_raw_pea_hexanal_baseline"
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
    save_json(process_state, PROCESS_STATE_CALIBRATIONS_PATH)
    save_json(backlog, BACKLOG_PATH)

    print(f"Successfully ingested all 6 Chunk 2 references!")


if __name__ == "__main__":
    ingest_chunk_2()
