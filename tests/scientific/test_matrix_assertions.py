import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import build_matrix_benchmark_assertions
from src.presentation import render_matrix_benchmark_assertions_markdown


def test_matrix_benchmark_assertions_cover_external_and_internal_matrix_rows():
    rows = build_matrix_benchmark_assertions([
        ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json",
        ROOT / "data" / "benchmarks" / "pea_isolate_ribose_cysteine_100C_45min_Internal2026.json",
    ])
    markdown = render_matrix_benchmark_assertions_markdown(rows)
    by_id = {row.benchmark_id: row for row in rows}

    pea_off = by_id["pea_isolate_40C_PratapSingh2021"]
    # RE-PINNED 2026-08-27 (Wave M) from "pass"/"pass" to "fail"/"fail". CAUSE: the Wave K/M
    # content correction on this benchmark (see its `content_correction_note`). Molecules
    # 2021, 26, 4104 Table 1 (PMC8271896) reports pea hexanal = 1138.00 ppb (repo had 260,
    # 4.38x low) and hexanol as n.d. (repo had 80 ppb, fabricated). The hexanal observability
    # factor was back-solved from the erroneous 260, so it still predicts 260.6 and the lane
    # now misses by 4.37x against a 2.0x tolerance; and with the paper's real values hexanal
    # outranks 2-pentylfuran, which the model gets the wrong way round. Both failures were
    # INVISIBLE while the reference agreed with the constant fitted to it. Nothing was
    # relaxed: the tolerance and the ranking contract are untouched, and the observability
    # factors were deliberately NOT refitted (owner decision -- see AUDIT.md Round 3).
    assert pea_off.ranking_contract_status == "order_mismatch"
    assert pea_off.max_ratio == pytest.approx(4.366, rel=1e-3)
    assert pea_off.ratio_status == "fail"
    assert pea_off.adverse_order_status == "fail"
    assert pea_off.overall_status == "fail"
    assert pea_off.strict_gate_blocked is True

    pea_meaty = by_id["pea_isolate_ribose_cysteine_100C_45min_Internal2026"]
    assert pea_meaty.top_k_status == "pass"
    assert pea_meaty.ratio_status == "pass"
    assert "mixed" in markdown
    assert "Matrix Benchmark Assertion Report" in markdown
