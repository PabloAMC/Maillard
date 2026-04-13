import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.primary_benchmark_campaign import (  # noqa: E402
    build_matrix_primary_benchmark_campaign,
    build_primary_matrix_external_package,
    build_primary_matrix_external_package_intake_template,
    render_matrix_primary_benchmark_campaign_markdown,
    render_primary_matrix_external_package_markdown,
    render_primary_matrix_external_package_intake_template_markdown,
)


def test_primary_matrix_campaign_materializes_dual_pea_soy_execution_bundle():
    payload = build_matrix_primary_benchmark_campaign(ROOT)

    assert payload["summary"]["protocol_id"] == "ppi_spi_ribose_cysteine_primary_benchmark_2026"
    assert payload["summary"]["selected_matrix"] == "pea_iso"
    assert payload["summary"]["selected_benchmark_id"] == "pea_isolate_ribose_cysteine_100C_45min_Internal2026"
    assert payload["summary"]["fallback_matrix"] == "soy_iso"
    assert payload["summary"]["arm_count"] == 2

    arms = {row["matrix"]: row for row in payload["arms"]}
    pea = arms["pea_iso"]
    soy = arms["soy_iso"]
    assert pea["would_close_requirements"] == [
        "comparator_is_measured_volatiles",
        "external_quantitative_origin",
        "minimum_quantitative_closed_targets",
    ]
    assert pea["remaining_after_protocol"] == ["no_internal_or_directional_dependencies"]
    assert pea["evidence_or_calibration_blockers"] == ["Hexanal", "Nonanal"]
    assert round(float(pea["hexanal_ratio"]), 3) == 1.074
    assert round(float(pea["nonanal_ratio"]), 3) == 0.985
    assert round(float(soy["hexanal_ratio"]), 3) == 1.000
    assert round(float(soy["nonanal_ratio"]), 3) == 1.000

    markdown = render_matrix_primary_benchmark_campaign_markdown(payload)
    assert "Matrix Primary Benchmark Campaign" in markdown
    assert "pea_isolate_ribose_cysteine_100C_45min_Internal2026" in markdown
    assert "Comparator signal is wet-lab measured_volatiles" in markdown


def test_primary_matrix_external_package_stays_pea_first_and_explicitly_not_measured():
    payload = build_primary_matrix_external_package(ROOT)

    assert payload["summary"]["selected_matrix"] == "pea_iso"
    assert payload["summary"]["status"] == "specified_not_yet_measured"
    assert payload["package"]["benchmark_candidate"] == "pea_isolate_ribose_cysteine_100C_45min_Internal2026"
    assert payload["package"]["current_external_anchor"] == "pea_isolate_40C_PratapSingh2021"
    assert payload["package"]["current_external_target_profile"] == "adverse_only"
    assert payload["package"]["would_close_requirements"] == [
        "comparator_is_measured_volatiles",
        "external_quantitative_origin",
        "minimum_quantitative_closed_targets",
    ]
    assert payload["package"]["remaining_after_package"] == ["no_internal_or_directional_dependencies"]

    markdown = render_primary_matrix_external_package_markdown(payload)
    assert "Primary Matrix External Package" in markdown
    assert "pea_isolate_40C_PratapSingh2021" in markdown
    assert "specified_not_yet_measured" in markdown


def test_primary_matrix_external_package_intake_template_is_shareable_but_placeholder_based():
    payload = build_primary_matrix_external_package_intake_template(ROOT)

    assert payload["protein_type"] == "pea_iso"
    assert payload["source_kind"] == "external_literature"
    assert payload["conditions"]["temp_C"] == 95.0
    assert payload["measured_volatiles"]["Hexanal"]["conc_ppb"].startswith("REPLACE_WITH_MEASURED_")
    assert payload["provenance"]["source_reference"] == "REPLACE_WITH_CONTRACT_LAB_OR_PUBLICATION"

    markdown = render_primary_matrix_external_package_intake_template_markdown(payload)
    assert "Primary Matrix External Intake Template" in markdown
    assert "Fill every measured_volatiles entry" in markdown