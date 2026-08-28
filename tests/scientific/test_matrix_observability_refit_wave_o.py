"""Guards on the Wave O matrix-observability refit (2026-08-27, owner-approved).

Why this file exists
--------------------
Wave O refitted the ambient-slurry hexanal observability factors against the
CONTENT-VERIFIED Pratap-Singh anchors (Molecules 2021, 26, 4104 Table 1, Europe PMC
PMC8271896: pea hexanal 1138.00 ppb, soy hexanal 1621.71 ppb), replacing constants that had
been back-solved from the 260 / 380 ppb transcription errors Wave K found.

A refit creates two durable ways to drift, and each gets a test here:

1. **The constants and their record can diverge.** The registry is edited by hand (so the
   provenance block lands with the values); the derivation lives in a generator. If someone
   edits one and not the other, the shipped model and the published derivation stop being
   the same thing, silently. ``test_shipped_constants_match_the_published_refit`` re-derives
   the fit and compares.

2. **A fitted row can quietly re-enter the evidence count.** Both benchmarks are fit targets
   and must stay out of the honest literature-coverage numerator and denominator. That is
   enforced through TWO independent mechanisms -- the registry's ``fitted_to_benchmark``
   label and the fit record's ``fit_leverage`` declaration -- and the second one is new with
   this wave. ``test_the_refit_record_declares_both_benchmarks_as_per_row_fit_targets``
   pins it.

The hold-out consequence of this refit is deliberately NOT tested here; it lives with the
other honest headlines in ``tests/scientific/test_honest_headline_guards.py`` (median fold
error 15.31x -> 42.62x, coverage 5/8 -> 4/8).
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.fit_target_index import fit_leverage_class, load_fit_records
from src.matrix_calibration_registry import _MATRIX_CALIBRATION_RECORDS

REFIT_RECORD = ROOT / "results" / "validation" / "matrix_observability_refit_pratap_singh.json"

FIT_TARGETS = (
    "pea_isolate_40C_PratapSingh2021",
    "soy_isolate_40C_PratapSingh2021",
)


def _record() -> dict:
    # Force-added in .gitignore precisely so this is readable in a fresh clone; a test that
    # skipped itself when the artifact was missing would be the pattern Wave J2 removed.
    assert REFIT_RECORD.exists(), (
        f"{REFIT_RECORD.relative_to(ROOT)} is missing. It is force-added in .gitignore "
        "because the registry cites it from the constants themselves and the fit_target_gate "
        "reads it to know these two benchmarks are fit targets."
    )
    return json.loads(REFIT_RECORD.read_text(encoding="utf-8"))


def test_shipped_constants_match_the_published_refit():
    """Wave O's fit is still shipped -- as a PRODUCT, which is what it always was.

    The record stores, per lane, the pre-Wave-O ``previous_value`` and the ``refit_value``
    it solved for. This originally walked the live registry and compared each shipped
    constant against its ``refit_value`` directly.

    RE-PINNED 2026-08-28 (Wave Y), AND MADE STRICTER RATHER THAN LOOSER. Wave Y relocated
    the shared ambient-hexanal scale out of ``observable_factor`` and into
    ``MATRIX_BENCHMARK_BASE_MARKER_YIELDS['Hexanal']`` (0.205 -> 0.885036), on the unit
    argument that an observability factor cannot exceed 1 while a marker yield multiplying
    an arbitrary ``hydroperoxide_scale`` can. Wave O's fit was NOT withdrawn and no
    parameter was added: the same single constant, the same two verified anchors, the same
    1.0113x residual.

    What Wave O actually solved for was never the factor on its own -- it was the PRODUCT
    ``marker_yield x observable_factor``, because that product is what reaches a prediction
    and the split between its two halves is a convention. So this test now checks the
    product, at the yield of the era each side belongs to. That is strictly stronger than
    the old form: it fails on any edit to EITHER constant that is not compensated in the
    other, which the old form could not see, and it is invariant to a future relocation
    that preserves predictions. Record:
    results/validation/matrix_marker_yield_rederivation.{json,md}.
    """
    from src.benchmark_validation import MATRIX_BENCHMARK_BASE_MARKER_YIELDS

    # The marker yield Wave O's refit_value was expressed against, before Wave Y moved the
    # scale onto it. Hard-coded because it is history, not a live constant.
    PRE_WAVE_Y_HEXANAL_YIELD = 0.205

    live = {
        (r.protein_type, r.process_state, r.compound): float(r.observable_factor)
        for r in _MATRIX_CALIBRATION_RECORDS
    }
    adopted = _record()["adopted"]
    assert adopted, "the refit record adopted nothing"

    for entry in adopted:
        lane = (entry["protein_type"], entry["process_state"], entry["compound"])
        assert lane in live, f"{lane} vanished from the calibration registry"
        assert lane[2] == "hexanal", (
            f"{lane} is not a hexanal lane; the product form below assumes the hexanal yield"
        )
        shipped_product = live[lane] * float(MATRIX_BENCHMARK_BASE_MARKER_YIELDS["Hexanal"])
        published_product = float(entry["refit_value"]) * PRE_WAVE_Y_HEXANAL_YIELD
        delta_dex = abs(math.log10(shipped_product / published_product))
        assert delta_dex <= 5.0e-4, (
            f"{lane} ships yield x factor = {shipped_product:.6g} but the published refit "
            f"says {published_product:.6g} ({delta_dex:.2e} dex apart). The fit and what "
            f"ships have diverged: either re-run "
            f"scripts/generators/refit_matrix_observability_pratap_singh.py and land the "
            f"new record, or restore the constants. NOTE the product, not the factor: "
            f"Wave Y moved the scale between the two halves without changing their product."
        )
        # The record must also still remember where it came from, or the move is
        # unauditable a year from now.
        assert float(entry["previous_value"]) != pytest.approx(float(entry["refit_value"])), (
            f"{lane} records an identical previous and refit value; the record can no "
            f"longer restate its own starting point."
        )


def test_the_one_parameter_fit_left_a_residual_rather_than_an_arithmetic_zero():
    """The 1.0113x residual is the informative part of this fit, so it is pinned.

    Two free factors against two rows would give exactly 1.0000x and would measure nothing.
    ONE shared scale leaves a residual, and its value says the two content-corrected anchors
    agree to 1.1% -- i.e. the transcription error was a common ABSOLUTE-SCALE error and the
    pea-vs-soy release ratio survived it. If this ever reads 1.0000, someone re-introduced
    per-lane freedom and destroyed that measurement.
    """
    record = _record()
    assert record["fit_leverage"]["free_parameters"] == 1, (
        "the refit is published as a ONE-parameter fit; adding parameters makes the residual "
        "below meaningless and must be argued for in AUDIT.md."
    )
    residuals = record["mutual_consistency"]["residual_dex_per_row"]
    assert len(residuals) == 2
    for residual in residuals:
        assert 10.0 ** abs(residual) == pytest.approx(1.0113, abs=5e-4), (
            f"shared-scale residual moved to {10.0 ** abs(residual):.4f}x (published "
            f"1.0113x)."
        )
    # The two per-lane scales bracket the shared one; that is why one row reads under and
    # the other over (asserted directly in test_matrix_only_benchmark.py).
    required = list(record["per_lane_required_scales"].values())
    shared = record["shared_scale_optimum"]
    assert min(required) < shared < max(required)

    assert record["bounds_check"]["hit_a_bound"] is False, (
        "the fitted scale is now within a rounding of a search bound, so the reported "
        "optimum may be an artefact of the range rather than of the data."
    )


def test_the_refit_record_declares_both_benchmarks_as_per_row_fit_targets():
    """Fitted rows stay out of the coverage count -- via the record, not only the registry.

    ``scripts/ci/fit_target_gate.py`` and ``src/fit_target_index.py`` both read this
    declaration. Losing it would not break any prediction; it would quietly let two
    fit-recovery rows count as literature coverage, which is the exact circularity the
    gate exists to block.
    """
    matching = [
        r for r in load_fit_records()
        if r["name"] == "matrix_observability_refit_pratap_singh.json"
    ]
    assert len(matching) == 1, (
        "the Wave O refit record is not being picked up as a fit record; check that it "
        "still matches one of src.fit_target_index.FIT_RECORD_GLOBS and still declares a "
        "fit_target_ids key."
    )
    record = matching[0]
    assert set(record["targets"]) == set(FIT_TARGETS), (
        f"fit targets moved to {sorted(record['targets'])}, published as {sorted(FIT_TARGETS)}"
    )
    assert fit_leverage_class(record) == "per_row_recovery", (
        "the refit is no longer classified per_row_recovery, so its two benchmarks would "
        "re-enter the honest literature-coverage numerator and denominator."
    )
