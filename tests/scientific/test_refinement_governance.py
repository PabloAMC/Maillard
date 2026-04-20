import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.selective_refinement_governance import (  # noqa: E402
    build_selective_refinement_governance_artifact,
    render_selective_refinement_governance_markdown,
)


pytestmark = [pytest.mark.slow]


BENCHMARKS = [
    ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json",
    ROOT / "data" / "benchmarks" / "pea_isolate_ribose_cysteine_100C_45min_Internal2026.json",
    ROOT / "data" / "benchmarks" / "soy_isolate_ribose_cysteine_100C_45min_Internal2026.json",
]


def test_refinement_governance_holds_offline_compute_until_cheap_screening_moves_benchmarks():
    payload = build_selective_refinement_governance_artifact(BENCHMARKS)

    assert payload["summary"]["governing_status"] == "hold_observable_first"
    assert payload["summary"]["mechanistic_priority_benchmark_count"] == 2
    assert payload["summary"]["approved_offline_job_count"] == 0
    assert any("cheap-first screening produced no benchmark-visible improvement" in blocker for blocker in payload["blockers"])
    assert len(payload["mechanistic_priority_benchmarks"]) == 2


def test_refinement_governance_markdown_surfaces_family_gate_and_blockers():
    markdown = render_selective_refinement_governance_markdown(build_selective_refinement_governance_artifact(BENCHMARKS))

    assert "Selective Mechanistic Refinement Governance" in markdown
    assert "Family Gate" in markdown
    assert "Geom Preopt" in markdown
    assert "No selective DFT jobs are currently approved." in markdown