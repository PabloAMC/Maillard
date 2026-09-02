"""FROZEN-WAVE REGRESSION RECORD (labelled 2026-09-03, test audit). The wave generator this file
tests is frozen (scripts/generators/WAVES.md); these tests fail only if the frozen report, the
network or the parameter tables change. They are the contract of a finished wave, not live checks
of new behaviour.


tests/unit/test_kinetic_core_b2_4.py

BUILD WAVE B2.4 -- the four things this wave must not be able to lose.

B2.4 is not a modelling wave either. It changes no stoichiometry, no species,
no reaction, no benchmark and no pass band. Its whole deliverable is that three
things that were IMPLICIT become DECLARED, and one broken citation becomes a
probe that runs. So these tests are not accuracy tests; they are tests that the
declarations are load-bearing and cannot silently rot:

  1. THE DECLARED WEIGHTING EXISTS, HAS THE DECLARED VALUES, AND THE CONTROL
     VALUE REPRODUCES B2.3's OBJECTIVE EXACTLY. If `shipped` ever stops meaning
     sigma_ph = 0.25, the control stops being a control and every cross-wave
     comparison in the wave report becomes a comparison of two different
     objectives.
  2. THE FOUR SHARED ROWS ARE DECLARED IN BOTH ARTIFACTS, and the two
     declarations name each other. D1 sec. 5 established that four exam points
     and four panel rows are THE SAME MEASUREMENTS; a declaration present in
     only one of the two scorecards would still let a reader treat them as
     independent evidence.
  3. THE ENSEMBLE REPORTS A SPREAD, NOT A POINT -- including the budget
     accounting, because a member that exhausted its evaluation budget is not
     a basin and must not enter a spread statistic.
  4. THE FREE SET IS WHAT THE PRE-REGISTRATION SAYS IT IS: 20 of 48, by four
     disjoint clauses, every key real.

Plus the firewall in the two forms B2.2 established -- a literal grep and a
SYSTEMS walk -- extended to this wave's new fit-side files.

They are unit tests. Anything that reads a generated artifact skips itself if
that artifact has not been generated, so the file runs in seconds either way.
The one test that evaluates the objective does so TWICE and no more.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
B23_FIT_REPORT = ROOT / "results/validation/kinetic_core_b2_3_fit_report.json"
ENSEMBLE = ROOT / "results/validation/kinetic_core_b2_4_ensemble.json"
PREREG = ROOT / "results/validation/kinetic_core_b2_4_prereg.md"
AMINE_PROBE = ROOT / "results/validation/kinetic_core_b2_4_amine_fate_probe.json"


def _executable_text(source: str) -> str:
    """
    Source with comments and string literals removed, so a check on CODE is not
    defeated -- in either direction -- by prose. A docstring that says
    "exam-blind" must not fail a blindness test, and a comment quoting a
    hold-out number must not pass a firewall test by hiding in a docstring.
    """
    import io
    import tokenize
    out = []
    try:
        for tok in tokenize.generate_tokens(io.StringIO(source).readline):
            if tok.type in (tokenize.COMMENT, tokenize.STRING):
                continue
            out.append(tok.string)
    except tokenize.TokenError:
        return source
    return " ".join(out)


def _executable_source(obj) -> str:
    import inspect
    return _executable_text(inspect.getsource(obj))

def _b24():
    import scripts.generators.generate_kinetic_core_b2_4_fit as mod
    return mod


# ===========================================================================
# 1. THE DECLARED WEIGHTING
# ===========================================================================


def test_the_exchange_rate_is_declared_with_all_three_values():
    b24 = _b24()
    assert set(b24.PH_ENDPOINT_WEIGHT) == {"shipped", "half", "measured"}
    assert b24.PH_ENDPOINT_WEIGHT["shipped"] == pytest.approx(1.40)
    assert b24.PH_ENDPOINT_WEIGHT["half"] == pytest.approx(0.70)
    assert b24.PH_ENDPOINT_WEIGHT["measured"] == pytest.approx(0.28)
    assert b24.SIGMA_LOG_REFERENCE == pytest.approx(0.35)
    # every declared value carries a written basis, and it is not a stub
    for tag in b24.PH_ENDPOINT_WEIGHT:
        assert len(b24.PH_ENDPOINT_WEIGHT_BASIS[tag]) > 120


def test_the_shipped_weighting_is_exactly_b2_3s_accidental_sigma():
    """
    THE CONTROL MUST BE A CONTROL. B2.3 hard-codes sigma = 0.25 pH units on the
    three `zhou_final_pH_*` rows. `shipped` is that rate expressed as an
    exchange rate: 0.35 / 0.25 = 1.40. If this drifts, the wave's control
    column silently becomes a fourth weighting.
    """
    b24 = _b24()
    sigma_ph = b24.SIGMA_LOG_REFERENCE / b24.PH_ENDPOINT_WEIGHT["shipped"]
    assert sigma_ph == pytest.approx(0.25, abs=1e-9)

    import scripts.generators.generate_kinetic_core_b2_3_fit as b23
    declared = {row["sigma_log"] for row in b23.ACTIVE_FIT_ROWS
                if row["kind"] == "ph_endpoint"}
    assert declared == {0.25}, (
        "B2.3's pH-endpoint sigma moved; the B2.4 control value must move with "
        f"it or stop being a control. Found {declared}.")


def test_the_weighting_applies_to_the_three_ph_endpoint_rows_and_nothing_else():
    b24 = _b24()
    assert set(b24.PH_ENDPOINT_ROW_IDS) == {
        "zhou_final_pH_from_pH6",
        "zhou_final_pH_from_pH7",
        "zhou_final_pH_from_pH8",
    }


@pytest.mark.skipif(not B23_FIT_REPORT.exists(), reason="B2.3 fit not generated")
def test_the_control_reproduces_b2_3s_cost_and_only_the_ph_rows_move():
    """
    TWO evaluations, and the whole point of both.

    At B2.3's own frozen vector: evaluated under `shipped` the objective must
    return B2.3's published cost EXACTLY, and evaluated under `measured` every
    non-pH residual must be bit-identical while the total cost moves. That is
    the demonstration that the declared scalar touches the three rows it says
    it touches and nothing else.
    """
    import numpy as np
    b24 = _b24()
    x = b24.b23_vector()

    r_shipped = b24.residual_vector(x, b24.PH_ENDPOINT_WEIGHT["shipped"], True)
    cost = 0.5 * float(np.dot(r_shipped, r_shipped))
    published = json.loads(B23_FIT_REPORT.read_text())["objective"]["final_cost"]

    # ==== WAVE B7 DISCLOSURE -- THIS BIT-IDENTITY IS BROKEN, ON PURPOSE ====
    # Through B6 this assertion was `cost == published` to 1e-9, and it was the
    # demonstration that B2.4's declared pH weighting reproduces B2.3's
    # objective exactly at B2.3's own vector.
    #
    # B7 hangs the FURANIC CHANNEL on the TRUNK, and the sulfur network runs
    # the trunk. Four of the eleven new steps drain a species the sulfur lane
    # reads (fructose, 3-deoxyglucosone, 1-deoxyglucosone, methylglyoxal), so
    # B2.3's residuals move even though NOT ONE of its 48 fitted constants has
    # been touched. The objective moves 8.1754 -> 8.3862, i.e. +2.6 %.
    #
    # WHAT THIS WAVE DID **NOT** DO ABOUT IT: refit. B7's licence is the
    # furanic channel's own declared FIT rows and nothing else; refitting the
    # sulfur lane to absorb a trunk change would be exactly the undeclared
    # sigma-exchange move that D1's diagnosis found accounted for 96 % of the
    # B2.2 -> B2.3 exam regression. The consequence is REPORTED instead: the
    # B2.4 vector is now fitted against a slightly different trunk from the one
    # that ships, the discrepancy is 2.6 % of the objective, and the largest
    # single exam row it moves is 1.26x (see the B7 hold-out report, H7).
    # An orchestrator ruling on whether to re-run B2.3/B2.4 against the B7
    # trunk is requested; nothing here presumes one.
    # ==== WAVE Q1: THE TOLERANCE, NOT THE VALUE, WAS THE DEFECT ====
    # B7 re-pinned this literal and kept the pre-B7 `rel=1e-9`. That tolerance
    # was correct when the assertion was `cost == published` -- both sides were
    # then the SAME float, read from the same artifact, so bit-identity was the
    # property under test and 1e-9 asserted it. After B7 the two sides are no
    # longer the same float: the left is a ~50-term least-squares cost RECOMPUTED
    # here, and 1e-9 silently became a demand for bit-reproducible floating-point
    # summation across BLAS/libm builds. That is not a property this repository
    # has ever claimed, and it does not hold: on macOS/arm64 outside the
    # validated container the cost is 8.386176922760855, which differs from the
    # pinned literal by 1.5e-5 absolute / 1.8e-6 relative -- ~1800x the declared
    # tolerance, and RED on a clean tree at 92b746d.
    #
    # The literal is NOT changed: it is the container's value and re-pinning it
    # to this host's would move a scientific number and hide the portability
    # fact. The tolerance is widened to 1e-5, which is ~5x the observed
    # cross-platform spread and still ~2600x tighter than the +2.6 % trunk
    # perturbation this assertion exists to detect. The precedent is B1's own B7
    # re-pins (tests/unit/test_kinetic_core_b1.py), which use rel=2e-3 for the
    # same class of quantity in the same wave.
    #
    # WHAT IS STILL BEING ASSERTED: that B2.3's objective at its own vector sits
    # where B7's trunk perturbation put it, and has not moved again. What is no
    # longer being asserted -- because it was never true -- is that a float sum
    # is bit-identical across platforms.
    # ==== WAVE B8: THE SECOND, MUCH LARGER PROPAGATION -- RE-PINNED ====
    # B7 moved this quantity by 2.6 % because the furanic channel drains four
    # trunk species the sulfur lane reads. B8 moves it by **2.69x**, and for a
    # different and much more direct reason: FIT_HOLDOUT_DECLARATION.md
    # Amendment 17 clause 5 replaces FOUR of the sulfur network's activation
    # energies with MEASUREMENTS the fit cannot move --
    #
    #   k_dimer_mft, k_dimer_fft   no barrier at all  ->  122.2 (Zhang k17)
    #   k_arp_dpo,   k_arp_tdp     the lumped 64.08   ->   85.7 (Zhang k16)
    #
    # -- and re-bands `Ea_decay_thiol_sink` from (20, 250) to Gigl's measured
    # (7, 102), which CLIPS B2.3's own 216.1 down to 102.
    #
    # 8.3862 -> 22.5575 IS THE SIZE OF THE DISAGREEMENT BETWEEN THE CORPUS'S
    # MEASURED BARRIERS AND B2.3'S FITTED ONES, evaluated at B2.3's vector. It
    # is not a defect and it is not absorbed: B8 refits from B2.4-half's
    # incumbent and reports this as its own start cost, so the movement is
    # visible rather than buried inside a converged number.
    #
    # WHAT IS STILL BEING ASSERTED, and it is the same property as before: that
    # B2.3's objective at its own vector sits where the last declared physics
    # change put it and has not moved AGAIN. The tolerance stays at the 1e-5
    # Q1 set, for the reason Q1 gave (cross-platform float summation, ~1.8e-6
    # relative on this objective).
    B8_T_STRUCTURE_SHIFTED_COST = 22.55754870634435
    assert cost == pytest.approx(B8_T_STRUCTURE_SHIFTED_COST, rel=1e-5), (
        f"B2.3's objective at its own vector has moved again: {cost}. It is "
        f"expected to differ from the published {published} by exactly the B7 "
        f"trunk perturbation plus B8's T-structure, and by nothing else.")
    assert cost / published == pytest.approx(2.759, abs=5e-3), (
        "the size of B8's T-structure perturbation of B2.3's objective has "
        "changed; it was 2.759x when Amendment 18 clause 10 recorded it")

    r_measured = b24.residual_vector(x, b24.PH_ENDPOINT_WEIGHT["measured"], True)
    ph_idx = {i for i, row in enumerate(b24.B23.ACTIVE_FIT_ROWS)
              if row["kind"] == "ph_endpoint"}
    for i in range(len(r_shipped)):
        if i in ph_idx:
            continue
        assert r_shipped[i] == pytest.approx(r_measured[i], rel=1e-12), (
            f"row {b24.B23.ACTIVE_FIT_ROWS[i]['id']} moved when only the pH "
            f"exchange rate changed")
    # and the comparator is genuinely weighting-independent
    assert b24.sum_r2_level(r_shipped) == pytest.approx(
        b24.sum_r2_level(r_measured), rel=1e-12)


# ===========================================================================
# 2. THE FOUR SHARED ROWS, DECLARED IN BOTH ARTIFACTS
# ===========================================================================


@pytest.mark.parametrize("tag", ("shipped", "half", "measured"))
def test_the_generated_artifacts_carry_the_shared_declaration(tag):
    exam_path = ROOT / f"results/validation/kinetic_core_b2_4_exam_{tag}.json"
    panel_path = ROOT / f"results/validation/kinetic_core_b2_4_panel_{tag}.json"
    if not exam_path.exists() or not panel_path.exists():
        pytest.skip(f"the {tag} re-sit has not been generated")

    exam = json.loads(exam_path.read_text())
    assert exam["shared_with_holdout_panel"]["n"] == 4
    assert "NOT independent evidence" in exam["shared_with_holdout_panel"][
        "declaration"]
    flagged = [r for r in exam["rows"] if r.get("shared_with_holdout_panel_row")]
    assert len(flagged) == 4

    panel = json.loads(panel_path.read_text())
    shared = [r for r in panel["rows"] if r.get("shared_with_cutover_exam")]
    assert len(shared) == 4


# ===========================================================================
# 3. THE SCORER CONDITIONS -- continuous statistics beside the censored ones
# ===========================================================================


def test_the_panel_reports_median_abs_log10_beside_the_gating_count():
    import scripts.generators.generate_kinetic_core_b2_3_holdout as panel
    rows = [
        {"gating": True, "pass": True, "fold_error": 1.0},
        {"gating": True, "pass": False, "fold_error": 100.0},
        {"gating": False, "pass": False, "fold_error": 0.01},
        {"gating": True, "pass": None, "fold_error": float("nan")},
    ]
    out = panel._continuous_scorecard(rows)
    assert out["median_abs_log10_fold_gating"] == pytest.approx(1.0)
    assert out["geometric_mean_fold_gating"] == pytest.approx(10.0)
    assert out["n_with_a_fold_gating"] == 2
    assert out["n_with_a_fold_all_scored"] == 3
    assert "CENSORED" in out["why_these_are_here"]


@pytest.mark.parametrize("tag", ("shipped", "half", "measured"))
def test_the_generated_panel_carries_the_continuous_statistic(tag):
    path = ROOT / f"results/validation/kinetic_core_b2_4_panel_{tag}.json"
    if not path.exists():
        pytest.skip(f"the {tag} panel has not been generated")
    scorecard = json.loads(path.read_text())["scorecard"]
    for key in ("gating_passed", "gating_rows",
                "median_abs_log10_fold_gating", "geometric_mean_fold_gating"):
        assert key in scorecard, f"{key} missing from the {tag} panel scorecard"
    assert scorecard["median_abs_log10_fold_gating"] is not None


# ===========================================================================
# 4. THE ENSEMBLE -- a spread, not a point
# ===========================================================================


def test_the_free_set_is_twenty_of_forty_eight_by_four_disjoint_clauses():
    b24 = _b24()
    assert len(b24.FREE_KEYS) == 20
    assert len(set(b24.FREE_KEYS)) == 20, "a key is claimed by two clauses"
    assert len(b24.FROZEN_KEYS) == 28
    assert set(b24.FREE_KEYS) | set(b24.FROZEN_KEYS) == set(b24.ALL_KEYS)
    assert not set(b24.FREE_KEYS) & set(b24.FROZEN_KEYS)
    clauses = (b24.FREE_R1_IDENTIFIED, b24.FREE_R2_PH_DRIFT,
               b24.FREE_R3_TOPOLOGY_SWITCHES, b24.FREE_R4_DECAY_BARRIERS)
    assert [len(c) for c in clauses] == [11, 2, 5, 2]
    for clause in clauses:
        for key in clause:
            assert key in b24.ALL_KEYS, f"{key} is not a parameter of this model"
            assert b24.FREE_CLAUSE_OF[key]


@pytest.mark.skipif(not B23_FIT_REPORT.exists(), reason="B2.3 fit not generated")
def test_clause_R1_is_exactly_what_b2_3_calls_identified():
    """
    The free set's first clause is not a hand-picked list: it is whatever the
    B2.3 report's own Gauss-Newton intervals call identified. If a regenerated
    B2.3 report ever changes that set, this test says so rather than letting a
    stale hard-coded tuple pass for a rule.
    """
    b24 = _b24()
    intervals = json.loads(B23_FIT_REPORT.read_text())["parameter_intervals"]
    identified = {k for k, v in intervals.items() if v.get("identified")}
    assert set(b24.FREE_R1_IDENTIFIED) == identified


def test_the_six_starts_are_stratified_and_start_zero_is_the_incumbent():
    import numpy as np
    b24 = _b24()
    base = b24.b23_vector()
    assert b24.start_arm(0) == "incumbent"
    assert np.array_equal(b24.start_vector(0, base), base)
    for start in b24.LOCAL_ARM_STARTS:
        assert b24.start_arm(start) == "local"
    for start in b24.GLOBAL_ARM_STARTS:
        assert b24.start_arm(start) == "global"
    # a local draw never leaves its declared neighbourhood; a global one does
    idx = list(b24.FREE_INDEX)
    local = np.abs(b24.start_vector(1, base) - base)[idx]
    rates = [i for i, key in enumerate(b24.FREE_KEYS) if key.startswith("k_")]
    assert local[rates].max() <= b24.LOCAL_HALF_WIDTH_DECADES + 1e-9
    # every start is reproducible from its seed
    assert np.array_equal(b24.start_vector(3, base), b24.start_vector(3, base))


def test_a_budget_exhausted_member_never_enters_the_spread():
    """
    THE RULE THAT MAKES THE SPREAD STATISTIC MEAN ANYTHING. A member stopped by
    the evaluation budget reports where it GOT TO, not where a basin is. The
    pre-registration excludes it; this pins that the code does.
    """
    b24 = _b24()
    members = [
        {"weight_tag": "shipped", "start": 0, "arm": "incumbent", "cost": 8.0,
         "converged": True, "budget_exhausted": False, "sum_r2_level": 15.0,
         "ph_endpoint_miss_pH_units": {"a": 0.1, "b": -0.2, "c": 0.3}},
        {"weight_tag": "shipped", "start": 1, "arm": "local", "cost": 16.0,
         "converged": True, "budget_exhausted": False, "sum_r2_level": 20.0,
         "ph_endpoint_miss_pH_units": {"a": 0.1, "b": -0.2, "c": 0.3}},
        {"weight_tag": "shipped", "start": 3, "arm": "global", "cost": 900.0,
         "converged": False, "budget_exhausted": True, "sum_r2_level": 800.0,
         "ph_endpoint_miss_pH_units": {"a": 3.0, "b": 3.0, "c": 3.0}},
    ]
    block = b24.ensemble_summary(members)["shipped"]
    assert block["n_members"] == 3
    assert block["n_converged"] == 2
    assert block["n_budget_exhausted"] == 1
    # log10(16/8) = 0.301, i.e. the 900 is excluded
    assert block["spread_S_log10_max_over_min"] == pytest.approx(0.30103, abs=1e-4)
    assert block["gate_G2_at_least_4_of_6_converged"] is False
    assert block["gate_G1_calibration_within_1pH"] is True


def test_the_shipping_criterion_is_mechanical_and_reads_no_score():
    """
    The criterion must be a function of the ensemble and nothing else. If it
    ever grows an argument, this test is where that is noticed.
    """
    import inspect
    b24 = _b24()
    signature = inspect.signature(b24.shipping_choice)
    assert list(signature.parameters) == ["summary"]
    # The check is on EXECUTABLE code, not on prose: the docstring is allowed
    # to say "exam-blind", and it does.
    code = _executable_source(b24.shipping_choice)
    for forbidden in ("exam", "holdout", "hold_out", "panel", "fold",
                      "median", "read_text", "json.load"):
        assert forbidden not in code.lower(), (
            f"the shipping criterion's CODE mentions {forbidden!r}; it is "
            f"supposed to be blind to both scorecards")

    summary = {
        "shipped": {"E": 1.40, "spread_S_log10_max_over_min": 0.90,
                    "gate_G1_calibration_within_1pH": True,
                    "gate_G2_at_least_4_of_6_converged": True},
        "half": {"E": 0.70, "spread_S_log10_max_over_min": 0.20,
                 "gate_G1_calibration_within_1pH": True,
                 "gate_G2_at_least_4_of_6_converged": True},
        "measured": {"E": 0.32, "spread_S_log10_max_over_min": 0.05,
                     "gate_G1_calibration_within_1pH": False,
                     "gate_G2_at_least_4_of_6_converged": True},
    }
    choice = b24.shipping_choice(summary)
    assert choice["shipped"] == "half", (
        "the smallest spread belongs to a weighting that FAILS gate G1 and "
        "must be disqualified")


@pytest.mark.skipif(not ENSEMBLE.exists(), reason="the ensemble has not been run")
def test_the_generated_ensemble_reports_a_spread_and_not_only_a_best():
    payload = json.loads(ENSEMBLE.read_text())
    assert payload["declared_weights"]
    for tag, block in payload["ensemble"].items():
        assert block["n_members"] >= 1
        for key in ("costs", "spread_S_log10_max_over_min", "n_converged",
                    "n_budget_exhausted", "sum_r2_level_per_member",
                    "gate_G1_calibration_within_1pH",
                    "gate_G2_at_least_4_of_6_converged"):
            assert key in block, f"{tag} ensemble is missing {key}"
        assert len(block["costs"]) == block["n_members"]
    choice = payload["shipping_choice"]
    assert choice["shipped"] in set(payload["declared_weights"]) | {None}
    assert "THE EXAM AND THE PANEL DO NOT ENTER IT" in choice["criterion"]


# ===========================================================================
# 5. THE AMINE-FATE PROBE -- it must RUN, which is the whole defect
# ===========================================================================


def test_the_rebuilt_amine_probe_does_not_reference_a_removed_species():
    """
    D1 sec. 7: `scratch/b23_encoding_probe.py` raises KeyError 'AMN' because no
    such species exists. The replacement must not reintroduce the same
    reference, and must derive the released amine from species that DO exist.
    """
    from src.kinetic_core.species_sulfur import SULFUR_INDEX
    text = _executable_text(
        (ROOT / "scripts/generators/probe_amine_fate_b2_4.py").read_text())
    assert not re.search(r"""["']AMN["']""", text), (
        "the rebuilt probe still touches the removed AMN species in "
        "executable code")
    import scripts.generators.probe_amine_fate_b2_4 as probe
    assert probe.AMINE_BEARING, "no amine-bearing species found"
    for key in probe.AMINE_BEARING:
        assert key in SULFUR_INDEX or key == "ACRCYS", (
            f"{key} is declared amine-bearing but is not a species")


@pytest.mark.skipif(not AMINE_PROBE.exists(), reason="the probe has not been run")
def test_the_rebuilt_probe_distinguishes_three_encodings():
    payload = json.loads(AMINE_PROBE.read_text())
    assert payload["verdict"]["encodings_are_distinct"], (
        "the rebuild reproduces the B2.3 defect: two 'encodings' that are the "
        "same code path")


# ===========================================================================
# 6. THE FIREWALL -- literal grep and SYSTEMS walk, extended to this wave
# ===========================================================================

_HOLDOUT_LITERALS = (
    "5.907", "11.439", "60.400", "12.757", "73.157",
    "62.6", "8.195",
    "6.88", "1.28", "3.29", "1.46", "2.4", "1.68", "1.71", "1.62",
    "525.62", "325.22", "50.07", "582.34", "436.63",
    "696.99", "813.65", "59.70", "553.0", "229.0",
)

_B2_4_FIT_SIDE_FILES = (
    "scripts/generators/generate_kinetic_core_b2_4_fit.py",
    "scripts/generators/probe_amine_fate_b2_4.py",
)


def test_no_holdout_literal_appears_on_the_b2_4_fit_side():
    """
    THIS TEST HAS ALREADY EARNED ITS KEEP. On its first run it found the
    literal `1.28` in the B2.4 fit generator: the third declared weighting's
    basis had been derived from D1's Hofmann pH-3-to-pH-7 slope, which is a
    HOLD-OUT quantity. The weighting was re-derived from Kumazawa's declared
    FIT ladder before any member at that weighting ran, and the whole exchange
    is recorded in the pre-registration as Correction C1.

    TWO ESCAPES, both narrow and both auditable:
      * a match preceded by a letter or an underscore is a WAVE NAME
        (`B2.4`, `b2_4`), not a measured value -- a measured value is never
        preceded by a letter, so the exemption costs the test nothing;
      * a line carrying the marker `FIREWALL-OK` is exempted, and the marker
        must be followed by a reason on the same line.
    """
    offences = []
    for relative in _B2_4_FIT_SIDE_FILES:
        text = (ROOT / relative).read_text()
        for literal in _HOLDOUT_LITERALS:
            pattern = re.compile(
                r"(?<![A-Za-z0-9._])" + re.escape(literal) + r"(?![0-9eE])")
            for line_no, line in enumerate(text.splitlines(), 1):
                if not pattern.search(line):
                    continue
                upper = line.upper()
                if "HOLD-OUT" in upper or "IS NOT HERE" in upper:
                    continue
                if "FIREWALL-OK" in upper:
                    continue
                offences.append(
                    f"{relative}:{line_no}: {literal!r} in {line.strip()[:90]}")
    assert not offences, (
        "hold-out literals leaked onto the B2.4 fit side:\n" + "\n".join(offences))


def test_every_firewall_ok_marker_carries_a_reason():
    """An escape hatch with no reason attached is a hole, not an escape."""
    for relative in _B2_4_FIT_SIDE_FILES:
        for line_no, line in enumerate(
            (ROOT / relative).read_text().splitlines(), 1
        ):
            if "FIREWALL-OK" not in line:
                continue
            tail = line.split("FIREWALL-OK", 1)[1].strip(" :-")
            assert len(tail) > 15, (
                f"{relative}:{line_no}: FIREWALL-OK with no stated reason")


def test_the_b2_4_fit_side_never_opens_a_frozen_bundle():
    io_tokens = ("open(", "read_text", "read_bytes", "json.load", "glob(",
                 "iterdir", "rglob")
    for relative in _B2_4_FIT_SIDE_FILES:
        for line_no, line in enumerate(
            (ROOT / relative).read_text().splitlines(), 1
        ):
            if "external_validation" not in line:
                continue
            assert not any(token in line for token in io_tokens), (
                f"{relative}:{line_no} reads a frozen bundle: {line.strip()[:100]}")


def test_the_b2_4_fit_integrates_only_b2_3s_own_fit_systems():
    """
    THE SYSTEMS WALK. B2.4 reuses B2.3's system table wholesale rather than
    declaring its own, so the walk is that it has not grown one -- the fastest
    way for a re-weighting wave to cheat would be to quietly add a pot.
    """
    import scripts.generators.generate_kinetic_core_b2_3_fit as b23
    b24 = _b24()
    assert Path(b24.B23.__file__) == Path(b23.__file__)
    assert not hasattr(b24, "SYSTEMS"), (
        "B2.4 has declared its own system table; it is supposed to reuse "
        "B2.3's unchanged")
    assert len(b24.B23.ACTIVE_FIT_ROWS) == len(b23.ACTIVE_FIT_ROWS) == 58
    for name, spec in b23.SYSTEMS.items():
        t_c, ph = float(spec["t_c"]), float(spec["ph"])
        if name.startswith("kang"):
            assert t_c in (100.0, 120.0), f"{name} at {t_c} C is the hold-out rung"
        if name.startswith("zhou") and not name.startswith("zhou_drift_"):
            assert ph == 7.0, f"{name} at pH {ph} is a held-out column"
        assert ph != 9.0 and ph != 6.5


def test_the_prereg_was_written_and_names_its_own_amendments():
    """
    A pre-registration that can be edited silently is not one. This does not
    prove WHEN it was written -- git does that -- but it pins that the two
    amendments forced by the container measurement are disclosed IN it rather
    than only in the wave report.
    """
    assert PREREG.exists()
    text = PREREG.read_text()
    for token in ("PH_ENDPOINT_WEIGHT", "AMENDMENTS A1 and A2",
                  "W-5", "W-10", "SHIPPING CRITERION",
                  "BEFORE any member finished"):
        assert token in text, f"the pre-registration does not mention {token}"
