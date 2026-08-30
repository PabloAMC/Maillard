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
    # RE-READ 2026-08-28 (Wave R1) — these two lines were an ACCIDENT on this three-benchmark
    # subset and are now a STRUCTURAL GUARANTEE, and the difference is the whole finding.
    # On the full panel the same code accepted NINE offsets of exactly +/-3.0 kcal/mol
    # (`retro_aldol`, `schiff*`, `thiol_addition*`), each kept because it lowered the score on
    # the benchmark panel the model is then evaluated against, and wrote them into
    # data/lit/refinement_surrogate_patches.json — which barrier_constants.get_barrier() ADDS
    # to the audited FAST_BARRIERS value. Armed, the shipped barriers were
    # schiff_condensation 18.0 (table 15.0), retro_aldol 29.0 (32.0), thiol_addition 31.6
    # (28.6): a ~35x rate factor each at 150 C. That is a fit to the evaluation set, and it
    # was never declared to scripts/ci/fit_target_gate.py.
    # These assertions passed throughout because THIS subset's candidates never advanced. A
    # test that only sees a three-benchmark panel cannot see a panel-wide fit; the guard that
    # can is tests/unit/test_wave_r1_barrier_offset_retirement.py, which pins
    # get_barrier(f) == FAST_BARRIERS[f][0] for every family and fails if accepted_offsets is
    # ever non-empty. Auto-acceptance is now retired, so both lines below hold by construction
    # on ANY panel; candidates are surfaced under `candidate_offsets_not_applied` instead.
    assert impact_payload["summary"]["accepted_candidate_count"] == 0
    assert impact_payload["summary"]["patched_total_score"] == impact_payload["summary"]["baseline_total_score"]


def test_get_barrier_consumes_refinement_surrogate_patch(monkeypatch, tmp_path):
    """The offset ARITHMETIC is deliberately kept; only automatic acceptance was retired.

    ANNOTATED 2026-08-28 (Wave R1). This test proves that an offset in the patch file
    reaches the shipped barrier — which is exactly the property that made the retired
    auto-acceptance dangerous, and exactly the property that a DECLARED barrier fit would
    need in future. It uses a `tmp_path` file, so it says nothing about the tracked one;
    the tracked file must stay empty, and
    tests/unit/test_wave_r1_barrier_offset_retirement.py is what asserts that. Keep both:
    delete this one and a future declared fit has no wiring test; delete that one and the
    2026-08-28 defect can return in silence.
    """
    patch_file = tmp_path / "refinement_surrogate_patches.json"
    patch_file.write_text(json.dumps({"accepted_offsets": {"strecker": 1.5}}), encoding="utf-8")

    monkeypatch.setattr(barrier_constants, "REFINEMENT_PATCH_FILE", patch_file)
    monkeypatch.setattr(barrier_constants, "_REFINEMENT_PATCH_CACHE", None)
    monkeypatch.setattr(barrier_constants, "_REFINEMENT_PATCH_MTIME", None)

    barrier, uncertainty = barrier_constants.get_barrier("strecker_degradation")

    assert round(barrier, 2) == 23.5
    assert uncertainty == 3.5