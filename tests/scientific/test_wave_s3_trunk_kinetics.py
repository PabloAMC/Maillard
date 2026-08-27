"""
tests/scientific/test_wave_s3_trunk_kinetics.py

Tests for the Wave S3 trunk integrator (`src/trunk_kinetics.py`) and for the
invariants of the trunk rate calibration that produced
`results/validation/trunk_rate_calibration_refit.{json,md}`.

DELIBERATELY NOT PINNED: the fitted rate constants themselves. Pinning the
number a wave has just produced is the circularity earlier rounds of this audit
removed. What IS pinned is everything that makes the number checkable —
stoichiometry, analytic limits, the declared fit corpus, the corpus exclusions,
and the fact that no shipped module imports the calibration lane.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from src.trunk_kinetics import (
    R_KJ,
    SPECIES,
    STEP_FAMILY,
    STEP_KEYS,
    STEP_UNIT,
    T_REF_K,
    arrhenius_k,
    carbon_like_balance,
    derivatives,
    integrate,
    rate_constants_at,
)

ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = ROOT / "results" / "validation" / "trunk_rate_calibration_refit.json"

ZERO = {key: (0.0, 100.0) for key in STEP_KEYS}


def params(**overrides):
    out = dict(ZERO)
    out.update(overrides)
    return out


# ---------------------------------------------------------------------------
# Arrhenius reparameterisation
# ---------------------------------------------------------------------------


def test_arrhenius_returns_k_ref_at_the_reference_temperature():
    assert arrhenius_k(0.037, 97.0, T_REF_K) == pytest.approx(0.037, rel=1e-12)


def test_arrhenius_matches_the_closed_form_off_reference():
    k = arrhenius_k(0.037, 97.0, 393.15)
    expected = 0.037 * math.exp(-(97.0 / R_KJ) * (1.0 / 393.15 - 1.0 / T_REF_K))
    assert k == pytest.approx(expected, rel=1e-12)


def test_arrhenius_is_increasing_in_temperature_for_positive_ea():
    assert arrhenius_k(1.0, 100.0, 353.15) < arrhenius_k(1.0, 100.0, 393.15)


def test_zero_activation_energy_is_temperature_independent():
    assert arrhenius_k(0.5, 0.0, 353.15) == pytest.approx(arrhenius_k(0.5, 0.0, 393.15))


def test_rate_constants_at_covers_every_declared_step():
    assert set(rate_constants_at(ZERO, T_REF_K)) == set(STEP_KEYS)


# ---------------------------------------------------------------------------
# Stoichiometry
# ---------------------------------------------------------------------------


def test_condensation_consumes_one_glucose_and_one_glycine_per_schiff_base():
    k = {key: 0.0 for key in STEP_KEYS}
    k["k_schiff"] = 1e-3
    d = derivatives(np.array([100.0, 100.0, 0.0, 0.0, 0.0, 0.0]), k)
    rate = 1e-3 * 100.0 * 100.0
    assert d[SPECIES.index("Glc")] == pytest.approx(-rate)
    assert d[SPECIES.index("Gly")] == pytest.approx(-rate)
    assert d[SPECIES.index("SB")] == pytest.approx(+rate)


def test_the_two_glycine_releasing_dfg_channels_regenerate_glycine():
    """DFG -> 3DG + Gly and DFG -> 1DG + Gly each return one amine."""
    for step, product in (("k_dfg_3dg", "TDG"), ("k_dfg_1dg", "ODG")):
        k = {key: 0.0 for key in STEP_KEYS}
        k[step] = 0.02
        d = derivatives(np.array([0.0, 0.0, 0.0, 10.0, 0.0, 0.0]), k)
        assert d[SPECIES.index("DFG")] == pytest.approx(-0.2)
        assert d[SPECIES.index("Gly")] == pytest.approx(+0.2)
        assert d[SPECIES.index(product)] == pytest.approx(+0.2)


def test_the_other_dfg_channel_releases_no_glycine():
    """This is the entire reason `k_dfg_other` exists as a separate step."""
    k = {key: 0.0 for key in STEP_KEYS}
    k["k_dfg_other"] = 0.02
    d = derivatives(np.array([0.0, 0.0, 0.0, 10.0, 0.0, 0.0]), k)
    assert d[SPECIES.index("DFG")] == pytest.approx(-0.2)
    assert d[SPECIES.index("Gly")] == pytest.approx(0.0)


def test_sugar_pool_is_non_increasing_and_only_the_declared_sinks_drain_it():
    run = integrate(
        params(k_schiff=(2e-5, 90.0), k_amadori=(0.1, 60.0), k_dfg_3dg=(0.01, 100.0)),
        T_REF_K, {"Glc": 200.0, "Gly": 200.0}, np.linspace(0, 200, 41),
    )
    pool = carbon_like_balance(run)
    # k_3dg_out / k_1dg_out / k_dfg_other / k_glc_other are all zero here, so
    # the pool must be conserved to integrator tolerance.
    assert np.allclose(pool, pool[0], rtol=1e-6)


def test_declared_sinks_do_drain_the_pool():
    run = integrate(
        params(k_glc_other=(0.01, 100.0)),
        T_REF_K, {"Glc": 200.0}, np.array([0.0, 100.0]),
    )
    pool = carbon_like_balance(run)
    assert pool[-1] < pool[0]


# ---------------------------------------------------------------------------
# Analytic limits
# ---------------------------------------------------------------------------


def test_isolated_first_order_decay_matches_the_closed_form_solution():
    """DFG alone with one open channel is a textbook exponential."""
    k = 0.0343
    run = integrate(
        params(k_dfg_3dg=(k, 97.0)),
        T_REF_K, {"DFG": 10.0}, np.array([0.0, 15.0, 30.0, 60.0, 90.0]),
    )
    expected = 10.0 * np.exp(-k * run.times_min)
    assert np.allclose(run.series("DFG"), expected, rtol=1e-6)


def test_glycine_release_integral_tracks_the_two_releasing_channels_only():
    k3, k4, kother = 0.02, 0.01, 0.05
    run = integrate(
        params(k_dfg_3dg=(k3, 90.0), k_dfg_1dg=(k4, 90.0), k_dfg_other=(kother, 90.0)),
        T_REF_K, {"DFG": 10.0}, np.array([0.0, 500.0]),
    )
    # asymptotic release = DFG0 * (k3 + k4) / (k3 + k4 + kother)
    expected = 10.0 * (k3 + k4) / (k3 + k4 + kother)
    assert run.glycine_released[-1] == pytest.approx(expected, rel=1e-4)


def test_fast_amadori_limit_reproduces_the_one_step_condensation():
    """
    With k_amadori enormous the Schiff pool is at steady state and the trunk
    collapses to a single bimolecular Glc + Gly -> DFG step. This is the limit
    Martins' own scheme assumes, and the limit the Schiff/Amadori profile
    likelihood is measured against.
    """
    k1 = 1.8e-5
    grid = np.array([0.0, 30.0, 90.0, 150.0])
    run = integrate(
        params(k_schiff=(k1, 90.0), k_amadori=(1.0e4, 60.0)),
        T_REF_K, {"Glc": 200.0, "Gly": 200.0}, grid,
    )
    # closed form for A + B -> C with equal initial concentrations:
    # 1/A(t) - 1/A0 = k*t
    a0 = 200.0
    expected = 1.0 / (1.0 / a0 + k1 * grid)
    assert np.allclose(run.series("Glc"), expected, rtol=2e-3)
    assert run.series("SB").max() < 1e-3 * a0


def test_slow_amadori_produces_a_visible_schiff_base_lag():
    """The counter-limit: a slow rearrangement makes DFG sigmoidal."""
    grid = np.linspace(0.0, 150.0, 31)
    run = integrate(
        params(k_schiff=(1.8e-5, 90.0), k_amadori=(2e-3, 60.0)),
        T_REF_K, {"Glc": 200.0, "Gly": 200.0}, grid,
    )
    dfg = run.series("DFG")
    # a lag means the curve is convex early: the first increment is smaller
    # than the second.
    assert dfg[1] - dfg[0] < dfg[2] - dfg[1]
    assert run.series("SB").max() > 0.1 * dfg.max()


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------


def test_unknown_species_is_rejected_rather_than_silently_ignored():
    with pytest.raises(KeyError):
        integrate(ZERO, T_REF_K, {"fructose": 1.0}, np.array([0.0, 1.0]))


def test_non_monotonic_time_grid_is_rejected():
    with pytest.raises(ValueError):
        integrate(ZERO, T_REF_K, {"Glc": 1.0}, np.array([0.0, 10.0, 5.0]))


def test_negative_time_is_rejected():
    with pytest.raises(ValueError):
        integrate(ZERO, T_REF_K, {"Glc": 1.0}, np.array([-1.0, 10.0]))


# ---------------------------------------------------------------------------
# Declarations that must not drift
# ---------------------------------------------------------------------------


def test_every_step_declares_a_family_slot_and_a_unit():
    for key in STEP_KEYS:
        assert key in STEP_FAMILY
        assert key in STEP_UNIT


def test_the_two_structural_gap_steps_are_declared_as_having_no_repo_family():
    """
    `k_glc_other` and `k_dfg_other` correspond to NO reaction family in this
    repository. That is a finding, not an oversight, and it must stay visible:
    quietly mapping either onto a neighbouring family would launder a fitted
    nuisance parameter into a chemistry claim.
    """
    assert STEP_FAMILY["k_glc_other"] is None
    assert STEP_FAMILY["k_dfg_other"] is None
    assert all(STEP_FAMILY[k] is not None
               for k in STEP_KEYS if k not in {"k_glc_other", "k_dfg_other"})


def test_the_bimolecular_step_is_the_only_one_with_bimolecular_units():
    bimolecular = [k for k in STEP_KEYS if STEP_UNIT[k] == "L/(mmol*min)"]
    assert bimolecular == ["k_schiff"]


def test_no_shipped_module_imports_the_calibration_lane():
    """
    `src/trunk_kinetics.py` is the calibration lane. If a prediction path ever
    starts importing it, the wave's central claim -- that no shipped number
    moved -- stops being true, and this test is where that must be noticed.
    """
    offenders = []
    for path in (ROOT / "src").rglob("*.py"):
        if path.name == "trunk_kinetics.py":
            continue
        if "trunk_kinetics" in path.read_text(encoding="utf-8"):
            offenders.append(path.relative_to(ROOT).as_posix())
    assert offenders == [], f"trunk_kinetics is imported by shipped code: {offenders}"


# ---------------------------------------------------------------------------
# The fit artifact's declarations
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not ARTIFACT.exists(), reason="calibration artifact not generated")
class TestCalibrationArtifact:
    @staticmethod
    def payload():
        return json.loads(ARTIFACT.read_text(encoding="utf-8"))

    def test_declares_its_fit_corpus_and_leverage(self):
        payload = self.payload()
        assert payload["fit_target_files"], "the fit corpus must be machine-declared"
        leverage = payload["fit_leverage"]
        assert leverage["free_parameters"] == 2 * len(STEP_KEYS)
        assert leverage["fitted_rows"] > 0
        assert leverage["class"] in {"per_row_recovery", "global_low_leverage"}

    def test_every_fit_target_lives_in_the_timeseries_directory(self):
        for target in self.payload()["fit_target_files"]:
            assert target.startswith("data/lit/timeseries/"), (
                f"fit corpus escaped the timeseries directory: {target}"
            )

    def test_the_benchmark_panel_and_holdout_are_declared_excluded(self):
        excluded = " ".join(self.payload()["fit_corpus_excluded"])
        assert "data/benchmarks/**" in excluded
        assert "external_validation" in excluded
        assert "directional_claims_panel" in excluded

    def test_no_holdout_bundle_id_is_named_anywhere_in_the_artifact(self):
        """
        The hold-out guard greps this file's raw text. Naming a bundle id here
        -- even to say it was excluded -- would fail the build, and rightly so.
        """
        assert "mp_holdout_" not in ARTIFACT.read_text(encoding="utf-8")

    def test_every_parameter_carries_a_ci_and_an_identifiability_verdict(self):
        for key, entry in self.payload()["parameters"].items():
            assert entry["identifiability_k"], f"{key} has no rate verdict"
            assert entry["identifiability_ea"], f"{key} has no Ea verdict"
            # a null CI is allowed, but ONLY together with an UNIDENTIFIED /
            # UNCONSTRAINED verdict -- never silently.
            if entry["k_ref_ci95"] is None:
                assert "UNIDENTIFIED" in entry["identifiability_k"]
            if entry["ea_ci95_kj_mol"] is None:
                assert "UNIDENTIFIED" in entry["identifiability_ea"]

    def test_confidence_intervals_were_widened_by_the_reduced_chi_square(self):
        objective = self.payload()["objective"]
        assert objective["confidence_interval_scaling"] == pytest.approx(
            math.sqrt(objective["reduced_chi_square"]), rel=1e-9
        )

    def test_no_parameter_is_resting_on_a_search_bound(self):
        stuck = [k for k, v in self.payload()["parameters"].items() if v["at_bound"]]
        assert stuck == [], f"parameters pinned at a search bound: {stuck}"

    def test_the_brands_cross_check_is_labelled_independent_and_martins_is_not(self):
        payload = self.payload()
        assert "INDEPENDENT" in payload["independent_cross_check_brands"]["independence"]
        assert "NOT INDEPENDENT" in payload["reproducibility_check_martins"]["independence"]
