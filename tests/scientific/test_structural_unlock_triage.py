import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.structural_unlock_triage import (  # noqa: E402
    build_structural_unlock_triage,
    render_structural_unlock_triage_markdown,
)


def test_structural_unlock_triage_prefers_primary_matrix_package_and_defers_broader_solver_work():
    payload = build_structural_unlock_triage(ROOT)

    assert payload["summary"]["ready_runtime_backlog_count"] == 0
    assert payload["summary"]["ready_benchmark_backlog_count"] == 0
    assert payload["summary"]["citation_backlog_exhausted"] is True

    ranked = payload["ranked_workstreams"]
    assert ranked[0]["workstream_id"] == "primary_matrix_benchmark_package"
    assert ranked[0]["recommendation"] == "run_now"
    assert ranked[1]["workstream_id"] == "extrusion_benchmark_translation"

    deferred = {row["item_id"]: row for row in payload["deferred_modeling_review"]}
    assert deferred["5.8"]["decision"] == "next_after_s17_baseline"
    assert deferred["5.7"]["decision"] == "defer"
    assert deferred["5.10"]["decision"] == "defer"
    assert deferred["5.11"]["decision"] == "defer"

    markdown = render_structural_unlock_triage_markdown(payload)
    assert "Structural Unlock Triage" in markdown
    assert "primary_matrix_benchmark_package" in markdown
    assert "5.8 Disulfide Bond Evolution / MFT Retention" in markdown
