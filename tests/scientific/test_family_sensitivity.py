import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.family_sensitivity import build_family_sensitivity_artifact, render_family_sensitivity_markdown
from src.family_lane_sensitivity import build_family_lane_sensitivity_payload, render_family_lane_sensitivity_markdown


def test_family_sensitivity_builds_barrier_family_impact_artifact_for_benchmark_visible_systems():
    payload = build_family_sensitivity_artifact(
        [
            ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json",
            ROOT / "data" / "benchmarks" / "pea_isolate_ribose_cysteine_100C_45min_Internal2026.json",
            ROOT / "data" / "benchmarks" / "soy_isolate_ribose_cysteine_100C_45min_Internal2026.json",
        ],
        delta_kcal=3.0,
        family_offset_keys={
            "strecker_degradation": "strecker",
            "cysteine_thermolysis": "cys",
        },
    )

    assert payload["summary"]["evaluated_benchmark_count"] == 3
    assert payload["summary"]["family_count"] == 2
    assert len(payload["families"]) == 2
    assert payload["families"][0]["reaction_family"] in {"strecker_degradation", "cysteine_thermolysis"}
    assert len(payload["families"][0]["scenarios"]) == 2
    assert payload["families"][0]["scenarios"][0]["benchmarks"]

    markdown = render_family_sensitivity_markdown(payload)
    assert "Family Sensitivity" in markdown
    assert "Reaction Family" in markdown
    assert "Benchmarks evaluated" in markdown


def test_family_lane_sensitivity_stays_separate_from_barrier_family_sensitivity():
    payload = build_family_lane_sensitivity_payload(
        {
            "family_lane_summary": {
                "02": {
                    "display_name": "Lipid oxidation and carbonylic crosstalk",
                    "active": True,
                    "strategic_posture": "immediate_expansion_lane",
                },
                "10": {
                    "display_name": "Microbial fermentation pretreatment",
                    "active": True,
                    "strategic_posture": "upstream_pretreatment_lane",
                },
            },
            "family_lane_adjustments": {
                "per_lane": {
                    "02": {"target_score_delta": -0.08, "maillard_closure_delta": -0.21, "off_flavour_risk_delta": 0.18},
                    "10": {"target_score_delta": 0.12, "off_flavour_risk_delta": -0.08},
                }
            },
        }
    )

    assert payload["summary"]["family_lane_count"] == 2
    assert payload["summary"]["sensitivity_policy"] == "family_lane_sensitivity_tracks_runtime_toggle_impact_not_barrier_offsets"
    assert payload["family_lanes"][0]["slr_family"] in {"02", "10"}

    markdown = render_family_lane_sensitivity_markdown(payload)
    assert "Family Lane Sensitivity" in markdown
    assert "Toggle Magnitude" in markdown