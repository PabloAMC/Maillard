import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import src.barrier_constants as barrier_constants
from src.refinement_campaign import (
    build_cheap_screening_artifact,
    build_global_sensitivity_artifact,
    build_refinement_impact_artifact,
    build_selective_dft_plan,
)


BENCHMARKS = [
    ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json",
    ROOT / "data" / "benchmarks" / "pea_isolate_ribose_cysteine_100C_45min_Internal2026.json",
    ROOT / "data" / "benchmarks" / "soy_isolate_ribose_cysteine_100C_45min_Internal2026.json",
]


def test_global_sensitivity_spans_barrier_process_and_formulation_axes():
    payload = build_global_sensitivity_artifact(BENCHMARKS)

    assert payload["summary"]["evaluated_benchmark_count"] == 3
    assert payload["summary"]["barrier_family_count"] >= 1
    assert payload["process_axes"]
    assert payload["formulation_axes"]
    assert payload["barrier_families"][0]["reaction_family"]


def test_refinement_campaign_discriminates_candidates_and_quarantines_non_improving_refinements():
    cheap_payload = build_cheap_screening_artifact(BENCHMARKS)
    dft_payload = build_selective_dft_plan(BENCHMARKS)
    impact_payload = build_refinement_impact_artifact(BENCHMARKS)

    assert cheap_payload["summary"]["candidate_count"] >= 1
    assert cheap_payload["summary"]["defer"] + cheap_payload["summary"]["reject"] >= 1
    assert any(row["no_escalation_reason"] != "none" for row in cheap_payload["candidates"])
    assert dft_payload["summary"]["candidate_count"] >= 1
    assert any(row["decision"] in {"defer", "reject"} for row in dft_payload["candidates"])
    assert impact_payload["summary"]["accepted_candidate_count"] == 0
    assert impact_payload["summary"]["patched_total_score"] == impact_payload["summary"]["baseline_total_score"]


def test_get_barrier_consumes_refinement_surrogate_patch(monkeypatch, tmp_path):
    patch_file = tmp_path / "refinement_surrogate_patches.json"
    patch_file.write_text(json.dumps({"accepted_offsets": {"strecker": 1.5}}), encoding="utf-8")

    monkeypatch.setattr(barrier_constants, "REFINEMENT_PATCH_FILE", patch_file)
    monkeypatch.setattr(barrier_constants, "_REFINEMENT_PATCH_CACHE", None)
    monkeypatch.setattr(barrier_constants, "_REFINEMENT_PATCH_MTIME", None)

    barrier, uncertainty = barrier_constants.get_barrier("strecker_degradation")

    assert round(barrier, 2) == 23.5
    assert uncertainty == 3.5