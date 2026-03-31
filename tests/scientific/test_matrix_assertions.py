import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import build_matrix_benchmark_assertions
from src.presentation import render_matrix_benchmark_assertions_markdown


def test_matrix_benchmark_assertions_cover_external_and_internal_matrix_rows(pea_matrix_assertion_rows):
    rows = pea_matrix_assertion_rows
    markdown = render_matrix_benchmark_assertions_markdown(rows)
    by_id = {row.benchmark_id: row for row in rows}

    pea_off = by_id["pea_isolate_40C_PratapSingh2021"]
    assert pea_off.overall_status == "pass"
    assert pea_off.adverse_order_status == "pass"
    assert pea_off.strict_gate_blocked is True

    pea_meaty = by_id["pea_isolate_ribose_cysteine_100C_45min_Internal2026"]
    assert pea_meaty.top_k_status == "pass"
    assert pea_meaty.ratio_status == "pass"
    assert "mixed" in markdown
    assert "Matrix Benchmark Assertion Report" in markdown
