import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.generators.generate_validation_figures import _build_payload


def test_validation_overview_payload_includes_quantitative_matrix_benchmarks():
    payload = _build_payload()

    quantitative = {row["benchmark_id"]: row for row in payload["quantitative_benchmarks"]}
    plotted_ids = {row["benchmark_id"] for row in payload["quantitative_points"]}

    assert "pea_isolate_uht_140C_Trikusuma2019" in quantitative
    assert quantitative["pea_isolate_uht_140C_Trikusuma2019"]["execution_path"] == "matrix_only"
    assert "pea_isolate_uht_140C_Trikusuma2019" in plotted_ids
    assert any(row["execution_path"] == "matrix_only" for row in payload["quantitative_points"])