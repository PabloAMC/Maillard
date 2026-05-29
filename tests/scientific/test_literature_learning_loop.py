import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.literature_learning_loop import (
    build_literature_learning_loop_payload,
    render_literature_learning_loop_markdown,
)


def test_learning_loop_payload_links_ready_references_to_runtime_artifacts():
    payload = build_literature_learning_loop_payload(ROOT)
    rows = {row["id"]: row for row in payload["ready_reference_rows"]}

    assert payload["summary"]["ready_reference_count"] == len(rows)
    assert payload["summary"]["ready_reference_count"] >= 5
    assert payload["summary"]["queue_conflict_count"] == 0
    assert rows["trikusuma_2019"]["triage_status"] == "encoded"
    assert rows["trikusuma_2019"]["runtime_artifacts_present"] is True
    for entry_id in {
        "trikusuma_2019",
        "pmc_2026_hme_hexanal_baseline",
        "acs_2022_pba_lysine_loss_benchmark",
        "acs_jafc_3c05991_ppi_spi_partitioning",
        "acs_jafc_3c02618_binding_prior",
        "acs_jafc_0c01925_protein_binding_hierarchy",
        "acs_apts_ref24_3dg_arrhenius_anchor",
        "mottram_nobrega_2002_furanone_bridge",
        "pmc_4419266_pe_interfacial_maillard_kinetics",
        "comunian_2021_thiamine_encapsulation",
        "voelker_2021_thiamine_kinetics",
        "huang_2022_thiamine_metal_catalysis",
        "wang_2012_gsh_xylose_sulfur_uplift",
        "blank_grosch_1991_hdmf_anchor",
        "liu_2023_ppi_offnote_baseline",
        "pmc11049305_spirulina_offnote_anchor",
        "pmc12155365_sunflower_roasted_anchor",
        "pmc_2024_pba_cml_cel_ranges_anchor",
        "pmc_12648097_acrylamide_mitigation_anchor",
        "frontiers_2022_hcw_aa_arrhenius_anchor",
        "pmid_1904866_pentosidine_equivalence_anchor",
        "scielo_brasil_aa_crosslink_hierarchy_anchor",
        "pmc5992167_amadori_pe_burden_anchor",
        "pmc9351765_crosspy_trapping_anchor",
        "uspto_ptacts_2023_yeast_extract_anchor",
        "wageningen_ref9_hme_rework_hydration_anchor",
        "acs_foodscitech_2024_hme_firmness_anchor",
    }:
        assert rows[entry_id]["encoding_status"] == "encoded_runtime_artifact"
    assert rows["trikusuma_2019"]["template_kind"] == "benchmark_payload"
    assert rows["resconi_2023_pbma_beef_identity_benchmark"]["template_kind"] == "benchmark_payload"
    assert rows["resconi_2023_pbma_beef_identity_benchmark"]["runtime_artifacts_present"] is True
    assert rows["lincoln_2025"]["template_kind"] == "computational_prior"
    assert rows["asen_2022"]["source_payload_role"] == "benchmark_intake"
    assert "process_state_calibration" in rows["asen_2022"]["target_payload_types"]
    assert "computational_prior" in rows["asen_2022"]["target_payload_types"]
    assert rows["lincoln_2025"]["encoding_status"] == "encoded_runtime_artifact"
    for entry_id, artifact_id in {
        "trikusuma_2019": "pea_isolate_uht_140C_Trikusuma2019",
        "asen_2022": "asen_2022_pea_denaturation",
        "pmc_2026_hme_hexanal_baseline": "li_2026_spi_wg_hme_hexanal_control_point",
        "acs_2022_pba_lysine_loss_benchmark": "acs_2022_pba_lysine_loss",
        "acs_jafc_3c05991_ppi_spi_partitioning": "acs_jafc_3c05991_ppi_spi_partitioning",
        "resconi_2023_pbma_beef_identity_benchmark": "resconi_2023_pbma_beef_identity_benchmark",
        "acs_apts_ref24_3dg_arrhenius_anchor": "pyrraline_from_3dg",
        "mottram_nobrega_2002_furanone_bridge": "mottram_nobrega_2002_furanone_sulfur_bridge_v1",
        "pmc_4419266_pe_interfacial_maillard_kinetics": "pe_schiff_base",
        "comunian_2021_thiamine_encapsulation": "comunian_2021_thiamine_encapsulation",
        "voelker_2021_thiamine_kinetics": "voelker_2021_thiamine_arrhenius_v1",
        "huang_2022_thiamine_metal_catalysis": "huang_2022_thiamine_metal_catalysis_v1",
        "wang_2012_gsh_xylose_sulfur_uplift": "wang_xu_glutathione_peptide_support_v1",
        "blank_grosch_1991_hdmf_anchor": "blank_grosch_1991_beef_hdmf_band",
        "liu_2023_ppi_offnote_baseline": "liu_2023_ppi_ibmp_band",
        "pmc11049305_spirulina_offnote_anchor": "pmc11049305_spirulina_beta_ionone_oav_floor",
        "pmc12155365_sunflower_roasted_anchor": "pmc12155365_sunflower_4_vinylguaiacol_fd_point",
        "pmc_2024_pba_cml_cel_ranges_anchor": "pmc_2024_pba_cml_cel_ranges",
        "pmc_12648097_acrylamide_mitigation_anchor": "pmc_12648097_acrylamide_mitigation",
        "frontiers_2022_hcw_aa_arrhenius_anchor": "frontiers_2022_hcw_aa_arrhenius_v1",
        "pmid_1904866_pentosidine_equivalence_anchor": "pmid_1904866_aa_pentosidine_equivalence_v1",
        "scielo_brasil_aa_crosslink_hierarchy_anchor": "scielo_brasil_aa_crosslink_hierarchy_v1",
        "pmc5992167_amadori_pe_burden_anchor": "pmc5992167_amadori_pe_food_matrix_burden",
        "pmc9351765_crosspy_trapping_anchor": "pmc9351765_crosspy_mft_scavenging_v1",
        "uspto_ptacts_2023_yeast_extract_anchor": "uspto_ptacts_2023_yeast_extract_mft_oav_band",
        "wageningen_ref9_hme_rework_hydration_anchor": "wageningen_ref9_hme_rework_hydration_collapse",
        "acs_foodscitech_2024_hme_firmness_anchor": "acs_foodscitech_2024_hme_firmness_window",
    }.items():
        assert any(item["artifact_id"] == artifact_id for item in rows[entry_id]["runtime_artifacts"])


def test_learning_loop_reviews_matrix_priors_and_structural_gaps():
    payload = build_literature_learning_loop_payload(ROOT)
    prior_rows = {row["protein_type"]: row for row in payload["matrix_prior_review"]}
    queue_rows = {row["chemistry_family"]: row for row in payload["payload_queue_review"]["queue_by_chemistry_family"]}
    promotion_queue = payload["s11_c_family_promotion_queue"]
    intake_gaps = {row["gap_id"]: row for row in payload["intake_structural_gap_review"]}
    backlog = payload["literature_backlog"]

    assert prior_rows["myco"]["has_accessibility_window"] is True
    assert "directional_only" in prior_rows["myco"]["uncertainty_postures"]
    assert "extrusion_structured" in prior_rows["myco"]["process_state_applicability"]
    assert payload["summary"]["payload_type_queue"]["benchmark_payload"] >= 1
    assert payload["summary"]["payload_type_queue"]["computational_prior"] >= 1
    assert queue_rows["alternative_protein_matrix_scope"]["payload_type_counts"]["process_state_calibration"] >= 1
    assert any(row["gap_id"] == "ppi_meaty_positive_matrix_benchmark" for row in payload["intake_structural_gap_review"])
    assert any(row["gap_id"] == "intact_spi_ppi_quantified_mft_fft" for row in payload["process_gap_review"])
    assert intake_gaps["ppi_meaty_positive_matrix_benchmark"]["closure_outcome"] == "wet_lab_only"
    assert intake_gaps["spi_meaty_positive_matrix_benchmark"]["triage_decision"] == "no_external_candidate_package"
    assert intake_gaps["spi_meaty_positive_matrix_benchmark"]["near_miss_candidates"][0]["entry_id"] == "nishimura_abe_2024"
    assert "same-run meaty sulfur markers" in intake_gaps["meaty_off_flavour_safety_tradeoff_panel"]["benchmark_contract_missing"]
    assert payload["summary"]["families_with_primary_payload_support"] >= 6
    assert payload["summary"]["wet_lab_only_intake_gap_count"] >= 3
    assert promotion_queue["selected_family"]["family_id"] == "carbonyl_donor_hierarchy"
    assert promotion_queue["fallback_family"]["family_id"] in {"thiamine_fragmentation_support", "fermentation_pretreatment"}
    assert promotion_queue["selected_family"]["minimum_runtime_landing"] == "benchmark_payload"
    assert promotion_queue["selected_family"]["reject_narrative_only"] is True
    assert backlog["summary"]["queue_conflict_count"] == 0
    assert any(row["id"] == "resconi_2023_pbma_beef_identity_benchmark" for row in backlog["encoded_reference_rows"])
    assert backlog["minimum_primary_experiment"]["exogenous_precursors"]["D-ribose_mM"] == 1.0

    markdown = render_literature_learning_loop_markdown(payload)
    assert "Literature Learning Loop" in markdown
    assert "Payload Queue Review" in markdown
    assert "Family Payload Coverage" in markdown
    assert "S11.C Family Promotion Queue" in markdown
    assert "Matrix Prior Review" in markdown
    assert "Intake Structural Gaps" in markdown
    assert "Runtime Present" in markdown
    assert "spi_meaty_positive_matrix_benchmark" in markdown
    assert "myco" in markdown