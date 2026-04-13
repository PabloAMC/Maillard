import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.extrusion_benchmark_execution import (  # noqa: E402
    build_example_extrusion_disulfide_follow_on_workbook,
    build_example_extrusion_external_closure_workbook,
    build_extrusion_disulfide_follow_on_execution_artifact,
    build_extrusion_external_closure_execution_artifact,
    export_extrusion_diagnostic_examples_bundle,
    render_extrusion_disulfide_follow_on_execution_markdown,
    render_extrusion_external_closure_execution_markdown,
)



def _completed_closure_workbook() -> dict:
    workbook = build_example_extrusion_external_closure_workbook(ROOT)
    for index, experiment in enumerate(workbook["experiments"]):
        experiment["source_kind"] = "external_literature"
        experiment["provenance"]["origin"] = "contract_lab_external"
        experiment["provenance"]["source_reference"] = f"Contract Lab Lot {index + 1}"
    return workbook


def _completed_follow_on_workbook() -> dict:
    return build_example_extrusion_disulfide_follow_on_workbook(ROOT)


def test_extrusion_external_closure_execution_turns_workbook_into_support_delta_bundle():
    workbook = _completed_closure_workbook()
    payload = build_extrusion_external_closure_execution_artifact(workbook, root=ROOT)

    assert payload["summary"]["experiment_count"] == 2
    assert payload["summary"]["all_direct_damage_markers_complete"] is True
    assert payload["summary"]["status"] == "ready_for_external_landing_review"
    first = payload["experiments"][0]
    assert first["intake_payload"]["source_kind"] == "external_literature"
    assert first["support_delta"]["aligned_benchmark"]["benchmark_id"] == "soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026"
    assert first["support_delta"]["promotion_assessment"]["landing_recommendation"] == "land_in_benchmark_candidate_or_blocker_registry"

    markdown = render_extrusion_external_closure_execution_markdown(payload)
    assert "Extrusion External Closure Execution" in markdown
    assert "ready_for_external_landing_review" in markdown


def test_extrusion_follow_on_execution_computes_5_8_metrics_against_reference_benchmark():
    workbook = _completed_follow_on_workbook()
    payload = build_extrusion_disulfide_follow_on_execution_artifact(workbook, root=ROOT)

    assert payload["summary"]["experiment_count"] == 2
    assert payload["summary"]["furosine_damage_gradient"] > 0.0
    first = payload["experiments"][0]
    second = payload["experiments"][1]
    assert first["metrics"]["free_sh_retention_fraction"] > second["metrics"]["free_sh_retention_fraction"]
    assert second["metrics"]["disulfide_pressure_proxy"] > first["metrics"]["disulfide_pressure_proxy"]
    assert second["metrics"]["furyl_disulfide_to_mft_ratio"] > first["metrics"]["furyl_disulfide_to_mft_ratio"]

    markdown = render_extrusion_disulfide_follow_on_execution_markdown(payload)
    assert "Extrusion 5.8 Follow-On Execution" in markdown
    assert "Furosine damage gradient" in markdown


def test_extrusion_diagnostic_examples_bundle_exports_both_flows_without_claiming_real_data(tmp_path):
    payload = export_extrusion_diagnostic_examples_bundle(tmp_path, root=ROOT)

    assert payload["summary"]["status"] == "diagnostic_examples_generated"
    assert payload["artifacts"]["closure_execution_status"] == "diagnostic_example_only"
    assert payload["artifacts"]["follow_on_execution_status"] == "diagnostic_example_only"
    assert (tmp_path / "extrusion_external_closure_diagnostic_example.yaml").exists()
    assert (tmp_path / "extrusion_disulfide_follow_on_diagnostic_example.yaml").exists()