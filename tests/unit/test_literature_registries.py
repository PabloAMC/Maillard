"""Literature-side registries: every builder builds, renders, and pins the decisions it encodes.

Merged 2026-09-03 (test audit) from ten single-test files that each imported one builder and
one renderer and never called the renderer. The value assertions are kept verbatim below,
one block per registry; the renderer smoke test covers all ten at once.
"""
from __future__ import annotations

import pytest

from src.family_ingestion_plan import build_family_ingestion_plan_artifact, render_family_ingestion_plan_markdown
from src.chemistry_family_scope import build_chemistry_family_scope_artifact, render_chemistry_family_scope_markdown
from src.family_strategy_policy import build_family_strategy_policy_artifact, render_family_strategy_policy_markdown
from src.matrix_family_coverage import build_matrix_family_coverage_artifact, render_matrix_family_coverage_markdown
from src.matrix_family_next_action import build_matrix_family_next_action_artifact, render_matrix_family_next_action_markdown
from src.matrix_family_priority_ranking import build_matrix_family_priority_ranking_artifact, render_matrix_family_priority_ranking_markdown
from src.mycoprotein_reference import build_mycoprotein_reference_artifact, render_mycoprotein_reference_markdown
from src.scope_gap_guard import build_scope_gap_guard_artifact, render_scope_gap_guard_markdown
from src.dha_lysinoalanine_external_package import build_dha_lysinoalanine_external_package_artifact, render_dha_lysinoalanine_external_package_markdown
from src.extrusion_external_closure import build_extrusion_external_closure_artifact, render_extrusion_external_closure_markdown


BUILDERS = [
    pytest.param(build_family_ingestion_plan_artifact, render_family_ingestion_plan_markdown, id="family_ingestion_plan"),
    pytest.param(build_chemistry_family_scope_artifact, render_chemistry_family_scope_markdown, id="chemistry_family_scope"),
    pytest.param(build_family_strategy_policy_artifact, render_family_strategy_policy_markdown, id="family_strategy_policy"),
    pytest.param(build_matrix_family_coverage_artifact, render_matrix_family_coverage_markdown, id="matrix_family_coverage"),
    pytest.param(build_matrix_family_next_action_artifact, render_matrix_family_next_action_markdown, id="matrix_family_next_action"),
    pytest.param(build_matrix_family_priority_ranking_artifact, render_matrix_family_priority_ranking_markdown, id="matrix_family_priority_ranking"),
    pytest.param(build_mycoprotein_reference_artifact, render_mycoprotein_reference_markdown, id="mycoprotein_reference"),
    pytest.param(build_scope_gap_guard_artifact, render_scope_gap_guard_markdown, id="scope_gap_guard"),
    pytest.param(build_dha_lysinoalanine_external_package_artifact, render_dha_lysinoalanine_external_package_markdown, id="dha_lysinoalanine_external_package"),
    pytest.param(build_extrusion_external_closure_artifact, render_extrusion_external_closure_markdown, id="extrusion_external_closure"),
]


@pytest.mark.parametrize("build, render", BUILDERS)
def test_every_registry_builds_and_renders_non_empty_markdown(build, render):
    payload = build()
    assert isinstance(payload, dict) and payload
    text = render(payload)
    assert isinstance(text, str) and text.strip().startswith("#")


# ---- from test_family_ingestion_plan.py ----
def test_family_ingestion_plan_prioritizes_first_wave_for_product_value():
    payload = build_family_ingestion_plan_artifact()
    by_family = {row["slr_family"]: row for row in payload["families"]}

    assert payload["summary"]["recommended_first_wave"][:4] == ["02", "07", "10", "08"]
    assert {"11", "12", "13"}.issubset(
        set(payload["summary"]["recommended_first_wave"])
    )
    assert by_family["02"]["runtime_concept"] == "lipid_crosstalk_lane"
    assert by_family["07"]["runtime_concept"] == "carbonyl_donor_hierarchy"
    assert by_family["10"]["runtime_concept"] == "fermentation_pretreatment_node"
    assert by_family["08"]["strategic_posture"] == "guardrail_lane"


def test_family_ingestion_plan_surfaces_deep_research_priority_queue():
    payload = build_family_ingestion_plan_artifact()
    summary = payload["summary"]
    by_family = {row["slr_family"]: row for row in payload["families"]}

    assert summary["backlog_family_count"] >= 5
    assert summary["backlog_citation_count"] >= 20
    assert summary["recommended_runtime_queue"][:3] == ["07", "11", "12"]
    assert summary["recommended_next_family"]["slr_family"] == "07"
    assert (
        "land donor identity"
        in summary["recommended_next_family"]["recommended_slice"]["focus"]
    )
    assert by_family["11"]["deep_research_backlog"]["high_confidence_count"] >= 4
    assert by_family["03"]["deep_research_backlog"]["high_confidence_count"] >= 3


# ---- from test_chemistry_family_scope.py ----
def test_chemistry_family_scope_recommends_lipid_crosstalk_as_next_family():
    payload = build_chemistry_family_scope_artifact()
    by_family = {row["family_id"]: row for row in payload["families"]}

    assert payload["summary"]["recommended_next_family"] == "lipid_oxidation_and_carbonylic_crosstalk"
    assert by_family["amino_acid_sugar_core"]["current_status"] == "first_class_core"
    assert by_family["lipid_oxidation_and_carbonylic_crosstalk"]["current_status"] == "partially_encoded_high_priority"
    assert by_family["nucleotide_and_ribose_support"]["current_status"] == "bounded_lane"


def test_chemistry_family_scope_uses_canonical_numbered_family_ids():
    payload = build_chemistry_family_scope_artifact()
    family_ids = {row["family_id"] for row in payload["families"]}

    assert "carbonyl_donor_hierarchy" in family_ids
    assert "fermentation_pretreatment" in family_ids
    assert "off_note_and_maillard_suppression" in family_ids


# ---- from test_family_strategy_policy.py ----
def test_family_strategy_policy_prioritizes_lipid_crosstalk_and_keeps_core_quantitative():
    payload = build_family_strategy_policy_artifact()
    summary = payload["summary"]

    assert summary["quantitative_trunk_family"] == "amino_acid_sugar_core"
    assert summary["default_next_expansion_family"] == "lipid_oxidation_and_carbonylic_crosstalk"
    assert summary["shared_ingestion_contract"]["machine_readable_only"] is True
    assert summary["lipid_crosstalk_dual_lane_policy"]["observable_lane_payloads"] == ["benchmark_payload", "retention_payload"]
    assert "amino_acid_sugar_core" in summary["family_lane_classification"]["first_class_core"]


# ---- from test_matrix_family_coverage.py ----
def test_matrix_family_coverage_artifact_distinguishes_explicit_from_indirect_support():
    payload = build_matrix_family_coverage_artifact()
    by_family = {row["matrix_family"]: row for row in payload["families"]}

    assert by_family["pea_isolate"]["runtime_posture"] == "directional_matrix"
    assert by_family["pea_isolate"]["support_class"] == "explicit_supported"
    assert by_family["soy_isolate"]["runtime_posture"] == "directional_matrix"
    assert by_family["coconut_oil_co_matrix"]["runtime_posture"] == "indirect_generic_support"
    assert by_family["coconut_oil_co_matrix"]["expansion_status"] == "blocked_on_family_specific_evidence"
    assert "coconut-specific lipid profile" in by_family["coconut_oil_co_matrix"]["what_is_not_supported"]
    assert by_family["mycoprotein"]["expansion_status"] == "bounded_expansion_candidate"
    assert "other_plant_proteins" in payload["summary"]["open_gap_families"]
    assert payload["summary"]["bounded_expansion_candidates"] == ["mycoprotein"]


# ---- from test_matrix_family_next_action.py ----
def test_matrix_family_next_action_advances_only_one_bounded_family():
    payload = build_matrix_family_next_action_artifact()

    assert payload["summary"]["chosen_family"] == "mycoprotein"
    assert payload["summary"]["chosen_action"] == "advance_bounded_next_family"

    by_family = {row["matrix_family"]: row for row in payload["decisions"]}
    assert by_family["mycoprotein"]["decision"] == "advance_now"
    assert by_family["coconut_oil_co_matrix"]["decision"] == "defer_as_scope_gap"
    assert by_family["other_plant_proteins"]["decision"] == "defer_as_scope_gap"
    assert by_family["pea_isolate"]["decision"] == "defer_until_primary_matrix_lane_moves"


# ---- from test_matrix_family_priority_ranking.py ----
def test_matrix_family_priority_ranking_orders_next_scope_choices():
    payload = build_matrix_family_priority_ranking_artifact()

    by_family = {row["matrix_family"]: row for row in payload["families"]}
    assert by_family["pea_isolate"]["scope_priority"] == "active_matrix_priority"
    assert by_family["soy_isolate"]["scope_priority"] == "active_matrix_priority"
    assert by_family["mycoprotein"]["scope_priority"] == "bounded_next_candidate"
    assert by_family["coconut_oil_co_matrix"]["scope_priority"] == "scope_gap_to_rank"
    assert by_family["other_plant_proteins"]["scope_priority"] == "scope_gap_to_rank"

    families_in_order = [row["matrix_family"] for row in payload["families"]]
    assert families_in_order.index("mycoprotein") < families_in_order.index("coconut_oil_co_matrix")
    assert families_in_order.index("coconut_oil_co_matrix") < families_in_order.index("other_plant_proteins")
    assert "mycoprotein" in payload["summary"]["bounded_next_candidates"]
    assert "coconut_oil_co_matrix" in payload["summary"]["scope_gap_priorities"]
    assert by_family["pea_isolate"]["evidence_landing"] == "external_benchmark_plus_calibration_review"
    assert by_family["mycoprotein"]["evidence_landing"] == "bounded_calibration_prior"
    assert by_family["coconut_oil_co_matrix"]["evidence_landing"] == "family_specific_calibration_or_tradeoff_benchmark"


# ---- from test_mycoprotein_reference.py ----
def test_mycoprotein_reference_uses_bounded_priors_and_next_family_decision():
    payload = build_mycoprotein_reference_artifact()

    assert payload["summary"]["matrix_family"] == "mycoprotein"
    assert payload["summary"]["decision"] == "advance_now"
    assert payload["summary"]["evidence_surface"] == "bounded_calibration_prior"
    assert payload["reference_windows"]["denaturation"]["midpoint_celsius"] == 78.0


# ---- from test_scope_gap_guard.py ----
def test_scope_gap_guard_blocks_scope_gap_families():
    payload = build_scope_gap_guard_artifact()

    assert payload["summary"]["blocked_family_count"] == 2
    assert payload["summary"]["blocked_families"] == ["coconut_oil_co_matrix", "other_plant_proteins"]


# ---- from test_dha_lysinoalanine_external_package.py ----
def test_dha_lysinoalanine_external_package_specifies_two_matrix_bundle():
    payload = build_dha_lysinoalanine_external_package_artifact()

    assert payload["summary"]["package_status"] == "specified_not_yet_measured"
    assert payload["summary"]["sprint_name"] == "dha_lysinoalanine_external_package"
    assert payload["summary"]["matrix_count"] == 2
    by_matrix = {row["matrix"]: row for row in payload["matrices"]}
    assert "Lysinoalanine (LAL)" in by_matrix["pea_iso"]["direct_damage_targets"]
    assert "direct_crosslink_marker_external_quantified" in by_matrix["soy_iso"]["would_close_requirements"]
    assert by_matrix["soy_iso"]["supportive_anchor_ids"] == ["de_leyn_2019_thiamine_retention"]


# ---- from test_extrusion_external_closure.py ----
def test_extrusion_external_closure_artifact_tracks_root_blocker_and_next_sprint():
    payload = build_extrusion_external_closure_artifact()

    assert payload["summary"]["matrix_count"] == 2
    assert payload["summary"]["root_blocker"] == "external_quantitative_direct_damage_markers_under_extrusion"
    assert payload["summary"]["selected_next_family_sprint"] == "dha_lysinoalanine_external_package"
    by_matrix = {row["matrix"]: row for row in payload["matrices"]}
    assert by_matrix["pea_iso"]["direct_closure_ready"] is False
    assert "Reactive lysine fraction" in by_matrix["pea_iso"]["direct_markers_missing"]
    assert by_matrix["soy_iso"]["supportive_anchor_ids"] == ["de_leyn_2019_thiamine_retention"]
    assert by_matrix["soy_iso"]["contextual_anchor_ids"] == ["troise_2018_soy_thermal_history"]
