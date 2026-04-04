import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import build_matrix_benchmark_evidence_audit
from src.presentation import render_matrix_benchmark_evidence_markdown


def test_matrix_evidence_audit_distinguishes_external_off_flavour_from_internal_meaty_candidates():
    rows = build_matrix_benchmark_evidence_audit([
        ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json",
        ROOT / "data" / "benchmarks" / "pea_isolate_ribose_cysteine_100C_45min_Internal2026.json",
        ROOT / "data" / "benchmarks" / "soy_isolate_ribose_cysteine_100C_45min_Internal2026.json",
        ROOT / "data" / "benchmarks" / "pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026.json",
        ROOT / "data" / "benchmarks" / "soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026.json",
    ])
    markdown = render_matrix_benchmark_evidence_markdown(rows)
    by_id = {row.benchmark_id: row for row in rows}

    pea_off = by_id["pea_isolate_40C_PratapSingh2021"]
    assert pea_off.external_data_status == "external_quantitative"
    assert pea_off.target_profile == "adverse_only"
    assert pea_off.promotable is False
    assert "no external meaty-positive targets" in pea_off.promotion_blocker

    pea_meaty = by_id["pea_isolate_ribose_cysteine_100C_45min_Internal2026"]
    assert pea_meaty.external_data_status == "internal_reference_only"
    assert pea_meaty.target_profile == "mixed"
    assert pea_meaty.promotable is False
    assert "internal reference-only" in pea_meaty.promotion_blocker

    soy_meaty = by_id["soy_isolate_ribose_cysteine_100C_45min_Internal2026"]
    assert soy_meaty.external_data_status == "internal_reference_only"
    assert soy_meaty.promotable is False

    pea_protocol = by_id["pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026"]
    assert pea_protocol.external_data_status == "internal_measured_quantitative"
    assert pea_protocol.target_profile == "mixed"
    assert pea_protocol.promotable is False
    assert "internal measured experiment" in pea_protocol.promotion_blocker

    soy_protocol = by_id["soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026"]
    assert soy_protocol.external_data_status == "internal_measured_quantitative"
    assert soy_protocol.promotable is False

    assert "Matrix Benchmark Evidence Audit Report" in markdown
    assert "internal_reference_only" in markdown
    assert "internal_measured_quantitative" in markdown