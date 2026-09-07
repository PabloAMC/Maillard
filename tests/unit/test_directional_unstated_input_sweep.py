"""2026-09-07: an unstated input is SWEPT, not invented -- the verdict stands only when unanimous."""
from __future__ import annotations

from src.kinetic_core import directional


def _claim(values):
    sys = lambda T: {"precursors": {"L-Asparagine": 10.0, "D-Glucose": 10.0}, "temp_C": T, "ph": 6.8, "water_activity": 0.99, "time_min": 20}
    return {
        "claim_id": "T-SWEEP", "statement": "acrylamide rises 140 -> 160 C", "claim_type": "ordering",
        "category": "temperature", "observables": ["Acrylamide"],
        "conditions": {"A": {"label": "160", "system": sys(160)}, "B": {"label": "140", "system": sys(140)}},
        "expected_relation": "A>B", "fit_status": "independent", "unstated_inputs": {"ph": values},
    }


def test_with_input_sets_every_arm_and_drops_the_field():
    c = directional._with_input(_claim([5.0, 7.0]), "ph", 5.0)
    assert c["conditions"]["A"]["system"]["ph"] == 5.0 and c["conditions"]["B"]["system"]["ph"] == 5.0
    assert "unstated_inputs" not in c


def test_a_unanimous_sweep_keeps_the_verdict_and_records_it():
    out = directional.score_claim(_claim([5.0, 6.8, 8.0]), flat_tolerance_pct=5.0)
    assert out["status"] == "agree"
    assert out["unstated_input_sweep"]["input"] == "ph" and set(out["unstated_input_sweep"]["statuses"].values()) == {"agree"}


def test_a_swept_value_outside_the_lane_window_is_held_not_refused():
    # pH 9 lies outside the acrylamide lane's measured window, but BOTH arms share the swept value, so
    # nothing is refused: the factor is held at pH 8 on both and the direction is still scored.
    out = directional.score_claim(_claim([5.0, 9.0]), flat_tolerance_pct=5.0)
    assert out["status"] == "agree"
    assert out["unstated_input_sweep"]["statuses"] == {"5.0": "agree", "9.0": "agree"}


def test_a_verdict_that_depends_on_the_unstated_input_is_not_evaluable(monkeypatch):
    calls = iter(["agree", "disagree"])
    real = directional.score_claim

    def fake(claim, *, flat_tolerance_pct):
        if "unstated_inputs" in claim:
            return real(claim, flat_tolerance_pct=flat_tolerance_pct)
        return {"claim_id": "T-SWEEP", "status": next(calls), "reason": None}

    monkeypatch.setattr(directional, "score_claim", fake)
    out = fake(_claim([5.0, 7.0]), flat_tolerance_pct=5.0)
    assert out["status"] == directional.NOT_EVALUABLE
    assert "depends on the unstated input ph (agree at ph 5, disagree at ph 7)" in out["reason"]
