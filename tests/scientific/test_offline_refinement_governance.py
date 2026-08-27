import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.offline_refinement_governance import (  # noqa: E402
    build_offline_refinement_governance_artifact,
    render_offline_refinement_governance_markdown,
)


BENCHMARKS = [
    ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json",
    ROOT / "data" / "benchmarks" / "pea_isolate_ribose_cysteine_100C_45min_Internal2026.json",
    ROOT / "data" / "benchmarks" / "soy_isolate_ribose_cysteine_100C_45min_Internal2026.json",
]


def test_offline_refinement_governance_holds_compute_until_cheap_screening_moves_benchmarks():
    payload = build_offline_refinement_governance_artifact(BENCHMARKS)

    assert payload["summary"]["governing_status"] == "hold_observable_first"
    # RE-PINNED 2026-08-27 (Wave S1 fix 2): ">= 1" -> exactly 0, and the assertion is now
    # TWO-SIDED rather than a floor, because a floor could not have detected this movement
    # in the other direction either. CAUSE: the compound-specific matrix calibration
    # registry became reachable from the `matrix_precursor_augmented` lane, so the
    # pea/soy Internal2026 Hexanal and Nonanal rows resolve a real record instead of
    # falling through to `heuristic`, and `_matrix_closure_action` stops scoring them
    # `mechanistic_blocker`. Full causal note in
    # tests/scientific/test_matrix_observable_closure_audit.py; the open item (the rule set
    # has no branch for `process_state_mismatch`) is carried in
    # tasks/audit_remediation.md (Wave S1).
    # THE GOVERNING DECISION IS UNCHANGED: still `hold_observable_first`, still zero
    # approved offline jobs, still the same three blockers. No compute was unlocked.
    assert payload["summary"]["mechanistic_priority_benchmark_count"] == 0
    assert payload["summary"]["approved_offline_job_count"] == 0
    assert any("cheap-first screening produced no benchmark-visible improvement" in blocker for blocker in payload["blockers"])

    assert payload["mechanistic_priority_benchmarks"] == [], (
        "a mechanistic-priority benchmark is back. If the matrix calibration registry "
        "lookup regressed to the SMILES-vs-name miss, this is how it shows up."
    )


def test_offline_refinement_governance_markdown_surfaces_family_gate_and_blockers():
    markdown = render_offline_refinement_governance_markdown(build_offline_refinement_governance_artifact(BENCHMARKS))

    assert "Offline Refinement Governance" in markdown
    assert "Family Gate" in markdown
    assert "No selective DFT jobs are currently approved." in markdown