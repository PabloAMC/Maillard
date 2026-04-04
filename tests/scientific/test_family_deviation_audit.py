import pytest

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.family_deviation_audit import (  # noqa: E402
    build_family_deviation_audit_artifact,
    render_family_deviation_audit_markdown,
)


def test_family_deviation_audit_reports_calibrated_error_tail():
    payload = build_family_deviation_audit_artifact()

    assert payload["summary"]["quantitative_point_count"] > 0
    assert payload["summary"]["families_with_quantitative_points"] >= 1
    assert payload["summary"]["max_observed_ratio"] is not None
    assert payload["summary"]["max_observed_ratio"] == pytest.approx(
        float(payload["summary"]["worst_point"]["compound_ratio"])
    )
    if payload["summary"]["max_observed_ratio"] >= payload["summary"]["high_ratio_threshold"]:
        assert payload["summary"]["high_ratio_point_count"] >= 1
        assert payload["summary"]["root_cause_counts"].get("model_or_mapping_mismatch", 0) >= 1
        assert len(payload["fix_queue"]) >= 1
    else:
        assert payload["summary"]["high_ratio_point_count"] == 0
        assert len(payload["fix_queue"]) == 0
    assert len(payload["top_outliers_by_ratio"]) > 0

    markdown = render_family_deviation_audit_markdown(payload)
    assert "Family Deviation Audit" in markdown
    assert "Top Outliers by Ratio" in markdown
    assert "Root-Cause Distribution" in markdown
    assert "(Cerny, 2008)" in markdown