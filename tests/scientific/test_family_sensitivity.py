import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.family_sensitivity import build_family_sensitivity_artifact, render_family_sensitivity_markdown


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