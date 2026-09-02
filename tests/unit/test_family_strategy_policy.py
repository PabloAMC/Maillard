import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.family_strategy_policy import build_family_strategy_policy_artifact, render_family_strategy_policy_markdown  # noqa: E402


def test_family_strategy_policy_prioritizes_lipid_crosstalk_and_keeps_core_quantitative():
    payload = build_family_strategy_policy_artifact()
    summary = payload["summary"]

    assert summary["quantitative_trunk_family"] == "amino_acid_sugar_core"
    assert summary["default_next_expansion_family"] == "lipid_oxidation_and_carbonylic_crosstalk"
    assert summary["shared_ingestion_contract"]["machine_readable_only"] is True
    assert summary["lipid_crosstalk_dual_lane_policy"]["observable_lane_payloads"] == ["benchmark_payload", "retention_payload"]
    assert "amino_acid_sugar_core" in summary["family_lane_classification"]["first_class_core"]


