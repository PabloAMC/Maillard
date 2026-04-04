from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.ts_seed_benchmark_validator import (  # noqa: E402
    _calculate_rmsd,
    build_ts_seed_assessment_artifact,
    render_ts_seed_assessment_markdown,
)


def test_ts_seed_rmsd_is_zero_for_identical_geometries():
    xyz = "3\nHCN\nH 0.0 1.0 0.0\nC -1.0 0.0 0.0\nN 1.0 0.0 0.0\n"

    assert _calculate_rmsd(xyz, xyz) == 0.0


def test_ts_seed_assessment_renders_markdown():
    payload = build_ts_seed_assessment_artifact()
    markdown = render_ts_seed_assessment_markdown(payload)

    assert payload["summary"]["benchmark_id"] == "mlp_ts_seed_benchmark_v1"
    assert "TS Seed Recovery Assessment" in markdown
    assert any(row["candidate_id"] == "aimnet2_shortlist" for row in payload["candidate_assessments"])