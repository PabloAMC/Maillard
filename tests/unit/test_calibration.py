"""`maillard calibrate` on a stand-in model: the machinery, not the engine (calibration_prereg.md T1, T3, T4).

The stand-in predicts log10 concentration = base(compound) + slope * (T - 145) / 10 + 0.4 log10(t / 20),
and moves with a calibration's overrides the way the real engine moves with a rate constant: a shift of
one named sulfur coordinate changes MFT by an amount that depends on the conditions. Measurements are
generated from it with a known threefold response factor and a known shift, and the fit must recover
both, with a pure level offset going to the factor and never to the coordinate.
"""
from __future__ import annotations

import ast
import math
from pathlib import Path

import numpy as np
import pytest

from src.kinetic_core import calibration as C
from src.kinetic_core import user_fit as F

ROOT = Path(__file__).resolve().parents[2]
TRUE_SHIFT = -0.40          # the laboratory's pot really runs this coordinate 0.4 dex slower
TRUE_FACTOR_LOG10 = math.log10(3.0)
COORD = "b8.k_nf_mft.log10_k_ref_145C"
MFT = "2-methyl-3-furanthiol"
FFT = "2-furfurylthiol"


def stand_in(true_shift: float):
    """A predictor with the engine's interface: (spec, compounds, calibration) -> Prediction."""

    def predict(spec, compounds, calibration):
        shift = true_shift                 # the laboratory's truth; the shipped model has 0.0 here
        if calibration is not None:
            for o in calibration.overrides:
                if o.coordinate.name == COORD:
                    shift += o.shift
        values = {}
        for c in compounds:
            base = {MFT: 2.0, FFT: 1.8}[c]
            slope = {MFT: 0.30, FFT: 0.10}[c]
            log10 = base + slope * (float(spec["temp_C"]) - 145.0) / 10.0 + 0.4 * math.log10(float(spec["time_min"]) / 20.0)
            if c == MFT:
                # the coordinate the calibration may move: like a rate constant, its effect depends on
                # the conditions (here it grows with time and temperature), so contrasts can see it,
                # while a pure level offset stays with the response factor
                g = 1.0 + 0.5 * (float(spec["temp_C"]) - 145.0) / 10.0 + math.log10(float(spec["time_min"]) / 20.0)
                log10 += shift * g
            values[c] = 10.0 ** log10
            if calibration is not None:
                rf = calibration.response_factors.get(c)
                if rf is not None:
                    values[c] *= rf.factor
        return F.Prediction(values, "sulfur")

    return predict


def measurements(n: int = 6):
    """Records generated from the stand-in with the TRUE shift and the TRUE factor, no noise."""
    truth = stand_in(TRUE_SHIFT)
    systems = []
    for i, (T, t) in enumerate([(120.0, 20.0), (130.0, 20.0), (140.0, 20.0), (145.0, 20.0), (145.0, 40.0), (145.0, 80.0)][:n]):
        p = truth({"temp_C": T, "time_min": t}, [MFT, FFT], None)
        systems.append({
            "name": f"pot_{i}", "precursors": {"L-Cysteine": 10.0, "D-Ribose": 10.0}, "temp_C": T, "time_min": t,
            "ph": 5.0, "aw": 0.98, "matrix": "water",
            "measured": {MFT: {"value": p.values[MFT] * 3.0, "uncertainty_pct": 10}, FFT: {"value": p.values[FFT], "uncertainty_pct": 10}},
            "quantification_class": "stable_isotope_dilution_gcms", "source": {"lab": "stand-in lab"},
        })
    return {"systems": systems}


@pytest.fixture(scope="module")
def result():
    return F.calibrate(measurements(), "stand-in lab", predict_fn=stand_in(0.0), max_coordinates=2)


def test_t1_the_response_factor_and_the_shift_are_recovered(result):
    cal, card = result
    assert abs(cal.response_factors[MFT].log10 - TRUE_FACTOR_LOG10) < 0.05
    assert abs(cal.response_factors[FFT].log10) < 0.05
    moved = {o.coordinate.name: o for o in cal.overrides}
    assert COORD in moved, card["diagnostics"]
    o = moved[COORD]
    assert abs(o.shift - TRUE_SHIFT) < max(2 * o.sigma, 0.05), (o.shift, o.sigma)


def test_t1_no_other_coordinate_moves(result):
    cal, _ = result
    for o in cal.overrides:
        if o.coordinate.name != COORD:
            assert abs(o.shift) < 0.05, (o.coordinate.name, o.shift)


def test_t3_the_holdout_is_never_read_during_the_fit(result):
    cal, card = result
    assert cal.validate_records, "with six records every second one is held out"
    assert not set(card["reads_during_fit"]) & set(cal.validate_records)


def test_the_holdout_improves_on_the_stand_in(result):
    _, card = result
    before, after = card["holdout"]["before"], card["holdout"]["after"]
    assert after["median_fold"] < before["median_fold"]
    assert after["median_fold"] < 1.2


def test_roles_are_assigned_before_the_fit_and_tags_are_kept():
    doc = measurements(6)
    doc["systems"][0]["role"] = "validate"
    records = F.records_from_document(doc)
    fit, validate, note = F.assign_roles(records)
    assert [r.name for r in validate] == ["pot_0"]
    assert len(fit) == 5
    few = measurements(3)
    fit, validate, note = F.assign_roles(F.records_from_document(few))
    assert not validate and "UNVALIDATED" in note


def test_contrasts_cancel_the_factor():
    """The log-ratio between two pots of the same laboratory does not depend on the factor."""
    records = F.records_from_document(measurements(4))
    fit, _, _ = F.assign_roles(records)
    contrasts = F.contrasts_of(fit)
    assert contrasts
    truth = stand_in(TRUE_SHIFT)
    by = {r.name: r for r in fit}
    for c in contrasts:
        pa = truth(by[c.a].spec, [c.compound], None).values[c.compound]
        pb = truth(by[c.b].spec, [c.compound], None).values[c.compound]
        assert abs(c.log10_ratio - math.log10(pb / pa)) < 1e-9


def test_the_calibration_round_trips_through_json(result, tmp_path):
    cal, _ = result
    path = cal.save(tmp_path)
    back = C.Calibration.load(path)
    assert back.as_dict() == cal.as_dict()
    assert path.parent == tmp_path / "stand-in_lab"


def test_overlay_writes_only_the_moved_coordinates():
    from src.kinetic_core.engine import SULFUR, frozen_parameters

    cands = {c.name: (c, v, s, b) for c, v, s, b in C.candidate_coordinates(SULFUR)}
    coord, value, sigma, band = cands[COORD]
    cal = C.Calibration("x", "b9", "2026-09-08", "water", {}, (C.Override(coord, value, sigma, value - 0.3, 0.1, band),), (), ())
    vector = cal.overlay(SULFUR)
    base = frozen_parameters(SULFUR)
    assert abs(vector["log10_k_ref_at_145C"]["k_nf_mft"] - (base["log10_k_ref_at_145C"]["k_nf_mft"] - 0.3)) < 1e-12
    untouched = {k: v for k, v in vector["log10_k_ref_at_145C"].items() if k != "k_nf_mft"}
    assert untouched == {k: v for k, v in base["log10_k_ref_at_145C"].items() if k != "k_nf_mft"}
    assert cal.overlay("trunk") is None and cal.overlay("lipid") is None


def test_candidate_coordinates_have_a_finite_prior_on_every_lane():
    for lane in ("trunk", "sulfur", "acrylamide"):
        cands = C.candidate_coordinates(lane)
        assert cands, lane
        for coord, value, sigma, band in cands:
            assert sigma >= 1e-3 and math.isfinite(value), coord.name
            assert coord.lane == lane


def test_t4_no_generator_or_engine_module_reads_user_results():
    """The shipped artifacts never see a calibration: nothing under scripts/generators or
    src/kinetic_core (other than the calibration and user-scoring modules) names results/user."""
    offenders = []
    for base in (ROOT / "scripts" / "generators", ROOT / "src" / "kinetic_core"):
        for path in base.rglob("*.py"):
            # the two index builders DESCRIBE the directory in their file tables; they read no record
            if path.name in ("calibration.py", "user_fit.py", "user_scoring.py", "build_results_readme.py", "build_data_readme.py"):
                continue
            if "results/user" in path.read_text(encoding="utf-8") or "USER_RESULTS_DIR" in path.read_text(encoding="utf-8"):
                offenders.append(str(path.relative_to(ROOT)))
    assert not offenders, offenders


def test_t4_the_engine_without_a_calibration_is_the_plain_engine():
    from src.comparative_cli import spec_to_core
    from src.kinetic_core.engine import predict

    spec = spec_to_core({"name": "x", "precursors": {"L-Cysteine": 10.0, "D-Ribose": 10.0}, "temp_C": 145.0,
                         "time_min": 20.0, "ph": 5.0, "aw": 0.98})
    a = predict(spec, [MFT])
    b = C.predict_calibrated(spec, [MFT], None)
    assert a.concentrations_ug_per_l == b.concentrations_ug_per_l


def test_a_calibration_from_another_matrix_warns_on_the_answer():
    from src.comparative_cli import spec_to_core

    cal = C.Calibration("x", "b9", "2026-09-08", "pea_isolate", {MFT: C.ResponseFactor(MFT, 0.3, 0.1, 2)}, (), (), ())
    spec = spec_to_core({"name": "x", "precursors": {"L-Cysteine": 10.0, "D-Ribose": 10.0}, "temp_C": 145.0,
                         "time_min": 20.0, "ph": 5.0, "aw": 0.98})
    run = C.predict_calibrated(spec, [MFT], cal)
    assert "warning" in run.run_metadata["calibration"]
    plain = C.predict_calibrated(spec, [MFT], None)
    assert abs(run.concentrations_ug_per_l[MFT] / plain.concentrations_ug_per_l[MFT] - 10**0.3) < 1e-9


def test_a_factors_sigma_widens_the_calibrated_interval():
    from src.comparative_cli import spec_to_core

    spec = spec_to_core({"name": "x", "precursors": {"L-Cysteine": 10.0, "D-Ribose": 10.0}, "temp_C": 145.0,
                         "time_min": 20.0, "ph": 5.0, "aw": 0.98})
    plain = C.predict_calibrated(spec, [MFT], None).absolutes()[MFT]
    cal = C.Calibration("x", "b9", "2026-09-08", "water", {MFT: C.ResponseFactor(MFT, 0.0, 0.3, 4)}, (), (), ())
    wide = C.predict_calibrated(spec, [MFT], cal).absolutes()[MFT]
    plain_w = math.log10(plain.hi_ug_per_l / plain.lo_ug_per_l)
    wide_w = math.log10(wide.hi_ug_per_l / wide.lo_ug_per_l)
    assert wide_w > plain_w
    # 1.645 * 0.3 decades half-width added in quadrature to the band's own half-width
    assert abs(wide_w / 2 - math.hypot(plain_w / 2, 1.645 * 0.3)) < 1e-6
