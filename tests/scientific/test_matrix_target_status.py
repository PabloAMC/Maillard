import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import build_matrix_target_status_artifact
from src.presentation import render_matrix_target_status_markdown


def test_matrix_target_status_distinguishes_quantitative_from_internal_and_directional_support():
    payload = build_matrix_target_status_artifact([
        ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json",
        ROOT / "data" / "benchmarks" / "soy_isolate_40C_PratapSingh2021.json",
        ROOT / "data" / "benchmarks" / "pea_isolate_ribose_cysteine_100C_45min_Internal2026.json",
        ROOT / "data" / "benchmarks" / "soy_isolate_ribose_cysteine_100C_45min_Internal2026.json",
        ROOT / "data" / "benchmarks" / "pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026.json",
        ROOT / "data" / "benchmarks" / "soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026.json",
    ])
    by_id = {row["benchmark_id"]: row for row in payload["benchmarks"]}

    pea_off = by_id["pea_isolate_40C_PratapSingh2021"]
    assert pea_off["support_counts"]["quantitative_closed"] >= 1
    assert pea_off["target_profile"] == "adverse_only"

    pea_meaty = by_id["pea_isolate_ribose_cysteine_100C_45min_Internal2026"]
    assert pea_meaty["support_counts"]["internal_candidate"] >= 1
    assert pea_meaty["support_counts"]["internal_reference_candidate"] >= 1
    assert pea_meaty["target_profile"] == "mixed"
    assert pea_meaty["promotion_ready"] is False
    assert pea_meaty["mechanistic_priority_ready"] is True
    assert pea_meaty["next_best_action"] == "prioritize_mechanistic_refinement"
    assert pea_meaty["promotion_blocker"] == "insufficient externally measured target closure; current comparator is internal reference-only"

    # Audit 2026-08-26: ProtocolPilot payloads are frozen model output; they
    # surface as synthetic and carry reference-grade support, never measured.
    pea_protocol = by_id["pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026"]
    assert pea_protocol["external_data_status"] == "synthetic_diagnostic_only"
    assert pea_protocol["support_counts"]["internal_reference_candidate"] >= 1
    assert pea_protocol["support_counts"]["internal_measured_candidate"] == 0
    assert pea_protocol["mechanistic_priority_ready"] is True
    assert pea_protocol["promotion_blocker"] == "insufficient externally measured target closure; current comparator is synthetic model output (diagnostic only)"

    summary = payload["summary"]
    assert summary["quantitative_closed"] >= 1
    assert summary["internal_candidate"] >= 1
    assert summary["internal_reference_candidate"] >= 1
    assert summary["mechanistic_priority_ready"] >= 1

    markdown = render_matrix_target_status_markdown(payload)
    assert "Matrix Target Status" in markdown
    assert "quantitative_closed" in markdown
    assert "internal_measured_candidate" in markdown
    assert "internal_reference_candidate" in markdown
    assert "synthetic_diagnostic_only" in markdown
    assert "prioritize_mechanistic_refinement" in markdown
