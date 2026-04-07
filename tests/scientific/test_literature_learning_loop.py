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
    assert rows["trikusuma_2019"]["encoding_status"] == "encoded_runtime_artifact"
    assert rows["trikusuma_2019"]["triage_status"] == "ready_benchmark"
    assert rows["trikusuma_2019"]["runtime_artifacts_present"] is True
    assert rows["pmc_2026_hme_hexanal_baseline"]["encoding_status"] == "encoded_runtime_artifact"
    assert rows["acs_2022_pba_lysine_loss_benchmark"]["encoding_status"] == "encoded_runtime_artifact"
    assert rows["trikusuma_2019"]["template_kind"] == "benchmark_payload"
    assert rows["lincoln_2025"]["template_kind"] == "computational_prior"
    assert rows["asen_2022"]["source_payload_role"] == "benchmark_intake"
    assert "process_state_calibration" in rows["asen_2022"]["target_payload_types"]
    assert "computational_prior" in rows["asen_2022"]["target_payload_types"]
    assert rows["lincoln_2025"]["encoding_status"] == "encoded_runtime_artifact"
    assert any(item["artifact_id"] == "pea_isolate_uht_140C_Trikusuma2019" for item in rows["trikusuma_2019"]["runtime_artifacts"])
    assert any(item["artifact_id"] == "asen_2022_pea_denaturation" for item in rows["asen_2022"]["runtime_artifacts"])
    assert any(item["artifact_id"] == "li_2026_spi_wg_hme_hexanal_control_point" for item in rows["pmc_2026_hme_hexanal_baseline"]["runtime_artifacts"])
    assert any(item["artifact_id"] == "acs_2022_pba_lysine_loss" for item in rows["acs_2022_pba_lysine_loss_benchmark"]["runtime_artifacts"])


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