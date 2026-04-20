import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.family_validation_overview import (  # noqa: E402
    build_family_validation_overview_artifact,
    render_family_validation_overview_markdown,
)


def test_family_validation_overview_tracks_quantitative_and_gap_families():
    payload = build_family_validation_overview_artifact([
        ROOT / "data" / "benchmarks" / "cys_glucose_150C_Farmer1999.json",
        ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json",
        ROOT / "data" / "benchmarks" / "soy_isolate_40C_PratapSingh2021.json",
    ])

    assert payload["summary"]["family_count"] >= 10
    assert payload["summary"]["quantitative_family_count"] >= 2
    assert payload["summary"]["quantitative_point_count"] >= 3
    assert any(row["benchmark_count"] > 0 for row in payload["families"])
    assert any(row["has_quantitative_parity"] for row in payload["families"])
    assert any(not row["has_quantitative_parity"] for row in payload["families"])

    markdown = render_family_validation_overview_markdown(payload)
    assert "Family Validation Overview" in markdown
    assert "compound-level quantitative parity" in markdown
    assert "| 01 |" in markdown


def test_family_validation_overview_tracks_integrated_family_runtime_support_and_late_qm_lanes():
    payload = build_family_validation_overview_artifact([
        ROOT / "data" / "benchmarks" / "acrylamide_spi_extrusion_130C_ACSRef3.json",
        ROOT / "data" / "benchmarks" / "cml_cel_commercial_pbma_Foods2023.json",
        ROOT / "data" / "benchmarks" / "furosine_extrusion_crossover_140C_RamirezJimenez2000.json",
    ])

    assert payload["summary"]["integrated_family_count"] >= 14
    by_slr = {row["slr_family"]: row for row in payload["families"]}

    assert by_slr["11"]["has_runtime_support"] is True
    assert by_slr["12"]["has_runtime_support"] is True
    assert by_slr["13"]["has_runtime_support"] is True
    assert by_slr["14"]["has_runtime_support"] is True
    assert by_slr["15"]["has_runtime_support"] is True
    assert by_slr["16"]["has_runtime_support"] is True
    assert by_slr["12"]["benchmark_count"] >= 3
    assert by_slr["15"]["benchmark_count"] == 0
    assert by_slr["16"]["benchmark_count"] == 0