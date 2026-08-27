import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.doe_generator import (  # noqa: E402
    build_extrusion_benchmark_protocol,
    build_extrusion_sop_lock_register,
    render_extrusion_benchmark_protocol_markdown,
    render_extrusion_sop_lock_register_markdown,
)
from src.extrusion_benchmark_landing import (  # noqa: E402
    build_extrusion_disulfide_follow_on_package,
    build_extrusion_disulfide_follow_on_workbook,
    build_extrusion_external_closure_package,
    build_extrusion_external_closure_workbook,
    render_extrusion_disulfide_follow_on_markdown,
    render_extrusion_disulfide_follow_on_workbook_markdown,
    render_extrusion_external_closure_package_markdown,
    render_extrusion_external_closure_workbook_markdown,
)


def test_extrusion_benchmark_protocol_is_repo_backed_but_marks_missing_locked_specs():
    payload = build_extrusion_benchmark_protocol(ROOT)

    assert payload["protocol_id"] == "spi_extrusion_mvp_benchmark_2026"
    assert payload["selected_protein_type"] == "soy_iso"
    # Updated 2026-08-27: the [120, 180] kJ/kg pair sat inside the SME response's
    # old saturated region, so the model predicted both arms identically. The
    # regenerated arms span the real 300-800 kJ/kg twin-screw window.
    assert payload["design_targets"]["sme_levels_kj_per_kg"] == [300.0, 700.0]
    assert payload["repo_backed_lab_specs"]["headspace_method"]["fiber"] == "DVB/CAR/PDMS"
    assert payload["repo_backed_lab_specs"]["headspace_method"]["aldehyde_internal_standard_identity"] == "hexanal-d12"
    assert payload["repo_backed_lab_specs"]["safety_method"]["platform"] == "UHPLC Nexera plus TQMS8040"
    assert payload["repo_backed_lab_specs"]["safety_method"]["internal_standard_identity"] == "13C3-acrylamide"
    assert payload["closest_repo_backed_hme_anchor"]["extrusion_moisture_wt_pct"] == 57.0
    assert payload["closest_repo_backed_hme_anchor"]["screw_speed_rpm"] == 280.0
    assert payload["closest_repo_backed_hme_anchor"]["feed_rate_kg_per_h"] == 4.6
    assert payload["closest_repo_backed_hme_anchor"]["barrel_temp_profile_C"] == [30.0, 90.0, 120.0, 140.0, 150.0, 160.0]
    assert payload["closest_repo_backed_hme_anchor"]["control_off_note_anchors_ug_per_kg"]["hexanal"] == 605.6
    missing_fields = {row["field"] for row in payload["missing_locked_lab_specs"]}
    assert "extruder_model" in missing_fields
    assert "thiol_internal_standard_concentrations" in missing_fields
    assert "hexanal_internal_standard_concentration" in missing_fields
    assert "acrylamide_internal_standard_concentration" in missing_fields

    markdown = render_extrusion_benchmark_protocol_markdown(payload)
    assert "Extrusion Benchmark Protocol" in markdown
    assert "review_ready_missing_locked_lab_specs" in markdown
    assert "DVB/CAR/PDMS" in markdown
    assert "57% moisture" in markdown
    assert "extruder_model" in markdown


def test_extrusion_sop_lock_register_locks_platform_but_keeps_model_and_spikes_unresolved():
    payload = build_extrusion_sop_lock_register(ROOT)

    assert payload["summary"]["status"] == "partially_locked_repo_backed"
    assert payload["locked_fields"]["extruder_platform"] == "twin_screw"
    assert payload["locked_fields"]["screw_speed_rpm"] == 280.0
    unresolved = {row["field"] for row in payload["unresolved_fields"]}
    assert "extruder_model" in unresolved
    assert "thiol_internal_standard_concentrations" in unresolved

    markdown = render_extrusion_sop_lock_register_markdown(payload)
    assert "Extrusion SOP Lock Register" in markdown
    assert "twin_screw" in markdown


def test_extrusion_external_closure_package_materializes_direct_damage_plus_tradeoff_bundle():
    payload = build_extrusion_external_closure_package(ROOT)

    assert payload["summary"]["contract_id"] == "extrusion_external_closure_2026"
    assert payload["summary"]["selected_matrix"] == "soy_iso"
    assert payload["summary"]["status"] == "specified_not_yet_measured"
    assert payload["summary"]["arm_count"] == 2
    direct_markers = {row["marker_id"] for row in payload["selected_matrix_contract"]["direct_marker_panel"]}
    assert {"reactive_lysine_fraction", "furosine", "lysinoalanine"}.issubset(direct_markers)
    assert "2-methyl-3-furanthiol" in payload["same_run_tradeoff_panel"]["same_run_minimum"]
    assert "Hexanal" in payload["same_run_tradeoff_panel"]["same_run_minimum"]

    markdown = render_extrusion_external_closure_package_markdown(payload)
    assert "Extrusion External Closure Package" in markdown
    assert "lysinoalanine" in markdown.lower()
    # Arm id updated 2026-08-27 with the regenerated DoE levels (120/180 ->
    # 300/700 kJ/kg); the arm ids are derived from the protocol's SME levels.
    assert "spi_hme_300_kj_per_kg" in markdown


def test_extrusion_external_closure_workbook_is_per_arm_and_placeholder_based():
    payload = build_extrusion_external_closure_workbook(ROOT)

    assert payload["package_id"] == "spi_extrusion_external_closure_workbook_2026"
    assert payload["selected_matrix"] == "soy_iso"
    assert len(payload["experiments"]) == 2
    assert payload["experiments"][0]["source_kind"] == "external_literature"
    assert payload["experiments"][0]["benchmark_alignment"]["target_benchmark_id"] == "soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026"
    assert payload["experiments"][0]["measured_damage_markers"]["furosine_mg_per_kg"].startswith("REPLACE_WITH_MEASURED_")
    assert "screw speed rpm" in payload["experiments"][0]["required_metadata"]

    markdown = render_extrusion_external_closure_workbook_markdown(payload)
    assert "Extrusion External Closure Workbook" in markdown
    assert "Requires direct markers" in markdown


def test_extrusion_5_8_follow_on_package_is_trigger_ready_and_runtime_aligned():
    payload = build_extrusion_disulfide_follow_on_package(ROOT)

    assert payload["summary"]["follow_on_id"] == "extrusion_5_8_disulfide_retention_follow_on_2026"
    assert payload["summary"]["status"] == "trigger_ready_pending_first_measured_spi_extrusion_panel"
    assert payload["supporting_runtime_sources"]["process_state_calibration_id"] == "raman_sds_extrusion_disulfide_severity"
    assert "acs_jafc_3c02618_mft_disulfide_trapping_v1" in payload["supporting_runtime_sources"]["retention_binding_prior_ids"]
    metric_ids = {row["metric_id"] for row in payload["derived_metrics"]}
    assert "free_sh_retention_fraction" in metric_ids
    assert "sulfur_to_pyrazine_retention_ratio" in metric_ids

    markdown = render_extrusion_disulfide_follow_on_markdown(payload)
    assert "Extrusion 5.8 Follow-On Package" in markdown
    assert "raman_sds_extrusion_disulfide_severity" in markdown


def test_extrusion_5_8_follow_on_workbook_keeps_feed_reference_and_post_extrusion_placeholders():
    payload = build_extrusion_disulfide_follow_on_workbook(ROOT)

    assert payload["package_id"] == "spi_extrusion_5_8_follow_on_workbook_2026"
    assert len(payload["experiments"]) == 2
    assert payload["experiments"][0]["reference_benchmark_id"] == "soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026"
    assert payload["experiments"][0]["feed_reference_assays"]["pre_extrusion_free_sh_umol_per_g"].startswith("REPLACE_WITH_")
    assert payload["experiments"][0]["post_extrusion_process_state_assays"]["furosine_mg_per_kg"].startswith("REPLACE_WITH_")
    assert payload["experiments"][0]["measured_volatiles_ppb"]["2-methyl-3-furanthiol"].startswith("REPLACE_WITH_MEASURED_")

    markdown = render_extrusion_disulfide_follow_on_workbook_markdown(payload)
    assert "Extrusion 5.8 Follow-On Workbook" in markdown
    assert "Feed reference assays" in markdown

def test_extrusion_doe_arms_are_predicted_to_differ():
    """Audit item 3.2: two arms the model scores identically buy no information."""
    payload = build_extrusion_benchmark_protocol(ROOT)
    discrimination = payload["arm_discrimination"]

    assert discrimination["arms_predicted_distinct"] is True
    for field in (
        "effective_temperature_celsius",
        "mean_residence_time_seconds",
        "relative_rtd_spread",
        "die_exit_temperature_celsius",
    ):
        assert field in discrimination["discriminating_fields"]

    signatures = [arm["predicted_signature"] for arm in payload["process_arms"]]
    assert len({tuple(sorted(sig.items())) for sig in signatures}) == len(signatures)

    # And the fields the model genuinely cannot separate are named, not hidden.
    assert "predicted_furosine_mg_per_kg" in discrimination["non_discriminating_fields"]
    markdown = render_extrusion_benchmark_protocol_markdown(payload)
    assert "Arm Discrimination Check" in markdown
