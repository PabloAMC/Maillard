"""
WAVE B45 (2026-09-11): the two completeness probes. Pre-registration:
results/validation/kinetic_core_b45_prereg.md.

B45 shipped no constant, so neither probe has a fit report or a ship rule to pin its digits.
The pre-registration nonetheless quotes both as evidence, and the argument it builds on them --
that the model is missing two CHANNELS rather than carrying two wrong constants -- fails if either
number drifts. These tests hold the registered predictions and the shape of the artifact.

They assert the PREDICTIONS, not the exact digits, in the two places where the digit is not the
claim: P1's claim is "under 1 %", P2's is "under 5 % lost". The digits are held only to the
precision the pre-registration prints them at, so that a legitimate change to an unrelated part of
the engine is not reported as a falsification of a claim it does not touch.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "results" / "validation" / "kinetic_core_b45_probe.json"
PREREG = ROOT / "results" / "validation" / "kinetic_core_b45_prereg.md"


def _payload():
    return json.loads(ART.read_text())


def test_the_preregistration_exists_and_records_its_outcome_after_its_predictions():
    text = PREREG.read_text()
    assert "## 4. The two probes, and the predictions" in text
    assert "## 6. Outcome (written 2026-09-11, after the probes)" in text
    assert text.index("## 4.") < text.index("## 6."), "the outcome must follow the predictions"


def test_p1_the_binding_block_stays_under_one_percent_on_every_programme():
    """The registered prediction: 'the model binds under 1 % of the hexanal'. Held for BOTH
    isolates and for programmes far harsher than the one Shi ran, because the claim is about the
    channel's size, not about one pot."""
    p1 = _payload()["p1_hexanal_binding"]
    for matrix in ("soy_isolate", "pea_isolate"):
        block = p1[matrix]
        assert block["charged"], f"{matrix} charges no amine pool"
        for label, row in block["programmes"].items():
            assert row is not None, f"{matrix}/{label} has no binding class"
            assert 0.0 < row["bound_fraction"] < 0.01, (matrix, label, row["bound_fraction"])
            lo, hi = row["bound_fraction_corners"]
            assert lo <= row["bound_fraction"] <= hi, (matrix, label)


def test_p1_is_four_orders_below_the_measured_release_and_points_the_other_way():
    """The finding, not the digit: the declared covalent block cannot account for what Shi
    measured, because she measured a RELEASE and the block can only BIND."""
    p1 = _payload()["p1_hexanal_binding"]
    shi = p1["measured_comparator"]
    assert shi["release_factor"] > 2.5, shi
    bound = p1["soy_isolate"]["programmes"]["shi_hold_95C_5min"]["bound_fraction"]
    assert round(bound, 6) == 0.000108, bound          # 0.0108 %, as the prereg prints it
    measured_swing = shi["release_factor"] - 1.0
    assert measured_swing / bound > 1_000.0, (measured_swing, bound)


def test_p2_the_engine_keeps_almost_all_of_a_thiol_only_pot():
    """The registered prediction: 'the model loses under 5 % of the charged cysteine' at 95 C,
    against a measurement in which nearly all of it is gone within 5 min of a 40-60 C ramp."""
    p2 = _payload()["p2_cysteine_removal"]
    five = p2["holds"]["5min"]
    assert five["answerable"]
    assert five["cys_fraction_remaining"] > 0.95, five
    assert round(five["cys_fraction_remaining"], 4) == 0.9949, five   # 99.49 %, as the prereg prints it
    assert p2["holds"]["180min"]["cys_fraction_remaining"] > 0.80, p2["holds"]["180min"]


def test_p2_removal_is_monotone_in_time_so_the_probe_is_reading_a_trajectory():
    p2 = _payload()["p2_cysteine_removal"]["holds"]
    fractions = [p2[k]["cys_fraction_remaining"] for k in ("5min", "30min", "180min")]
    assert fractions == sorted(fractions, reverse=True), fractions


def test_the_probe_declares_that_it_is_not_scored():
    p = _payload()
    what = p["what_this_is"].lower()
    assert "completeness" in what and "neither is scored" in what
    assert p["prereg"] == "results/validation/kinetic_core_b45_prereg.md"
