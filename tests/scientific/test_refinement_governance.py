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
    # RE-PINNED 2026-08-27 (Wave S1 fix 2): 2 -> 0. CAUSE, and it is a lookup repair rather
    # than a science change: the compound-specific matrix calibration registry is now
    # REACHABLE from the `matrix_precursor_augmented` lane (species used to arrive labelled
    # by canonical SMILES and miss the name-keyed registry entirely), so the Hexanal /
    # Nonanal rows on the two Internal2026 benchmarks stop reading
    # `calibration_evidence_strength = "heuristic"` and stop being scored
    # `mechanistic_blocker` by `_matrix_closure_action`. No mechanistic-priority benchmark
    # survives, so the priority list is empty.
    # THE GOVERNING DECISION IS UNCHANGED and that is the point of pinning 0 rather than
    # relaxing the assertion: `governing_status` is still `hold_observable_first` and
    # `approved_offline_job_count` is still 0, so no DFT compute is unlocked by this. What
    # was lost is a WARNING, not a gate -- see the dated block in
    # tests/scientific/test_matrix_observable_closure_audit.py for why
    # `process_state_mismatch` should probably have its own branch in the rule set, and
    # tasks/audit_remediation.md (Wave S1) for the open item.
    assert payload["summary"]["mechanistic_priority_benchmark_count"] == 0
    assert payload["summary"]["approved_offline_job_count"] == 0
    assert any("cheap-first screening produced no benchmark-visible improvement" in blocker for blocker in payload["blockers"])
    assert payload["mechanistic_priority_benchmarks"] == []


def test_refinement_governance_markdown_surfaces_family_gate_and_blockers():
    markdown = render_selective_refinement_governance_markdown(build_selective_refinement_governance_artifact(BENCHMARKS))

    assert "Selective Mechanistic Refinement Governance" in markdown
    assert "Family Gate" in markdown
    assert "Geom Preopt" in markdown
    assert "No selective DFT jobs are currently approved." in markdown