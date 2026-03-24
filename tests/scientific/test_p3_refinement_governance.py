import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.p3_refinement_governance import (  # noqa: E402
    build_p3_refinement_governance_artifact,
    render_p3_refinement_governance_markdown,
)


BENCHMARKS = [
    ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json",
    ROOT / "data" / "benchmarks" / "pea_isolate_ribose_cysteine_100C_45min_Internal2026.json",
    ROOT / "data" / "benchmarks" / "soy_isolate_ribose_cysteine_100C_45min_Internal2026.json",
]


def test_p3_governance_holds_offline_compute_until_cheap_screening_moves_benchmarks():
    payload = build_p3_refinement_governance_artifact(BENCHMARKS)

    assert payload["summary"]["governing_status"] == "hold_observable_first"
    assert payload["summary"]["mechanistic_priority_benchmark_count"] >= 1
    assert payload["summary"]["approved_offline_job_count"] == 0
    assert any("cheap-first screening produced no benchmark-visible improvement" in blocker for blocker in payload["blockers"])

    benchmark_ids = {row["benchmark_id"] for row in payload["mechanistic_priority_benchmarks"]}
    assert "pea_isolate_ribose_cysteine_100C_45min_Internal2026" in benchmark_ids

    pea_row = next(row for row in payload["mechanistic_priority_benchmarks"] if row["benchmark_id"] == "pea_isolate_ribose_cysteine_100C_45min_Internal2026")
    assert "Hexanal" in pea_row["target_compounds"]
    assert "adverse-marker closure" in pea_row["expected_decision_change"]


def test_p3_governance_markdown_surfaces_family_gate_and_blockers():
    markdown = render_p3_refinement_governance_markdown(build_p3_refinement_governance_artifact(BENCHMARKS))

    assert "P3 Refinement Governance" in markdown
    assert "Family Gate" in markdown
    assert "No selective DFT jobs are currently approved." in markdown