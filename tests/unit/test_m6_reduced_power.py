"""Case M6's falsification arithmetic and reduced mini-model (src/m6_reduced_power.py): the structural predictions the
case's rivals are defined by, the analogue's carry, the job (a) verdicts the case quotes, and the declared noise in the
adopted judge's terms. Synthetic and printed values only; no engine."""
from __future__ import annotations

import math

import numpy as np
import pytest

from src import m6_reduced_power as M
from src.kinetic_core import predictive_densities as pdn


def _pred(rival, n=400, seed=0):
    th, lr = M.draw(rival, n, np.random.default_rng(seed), defensive=False)
    assert np.all(lr == 0)
    return M.predict_c(th, rival)


def test_elfsson_carry_matches_the_dossier():
    # the dossier: 0.76 M^-1 min^-1 at 100 C, pH 7.5 from the printed Arrhenius line (5.2 M^-1 h^-1 at 50 C); this
    # carry starts from the printed 5.53 at pH 7.42 and adds the pH slope: within 0.1 dex
    assert abs(math.log10(M.elfsson_k2(100.0, 7.5) / 0.76)) < 0.1
    assert M.ELFSSON_PH_SLOPE == pytest.approx(math.log10(80.1) / 2.61, rel=1e-3)


def test_step_units_are_structural_for_c0_and_c3():
    c0, c3 = _pred("C0"), _pred("C3")
    assert np.median(c0["mgo_step"]) == pytest.approx(math.log10(3), abs=0.02)
    assert np.median(c3["mgo_step"]) == pytest.approx(math.log10(9), abs=0.04)
    assert np.median(c0["dilution_step"]) == pytest.approx(math.log10(1 / 16), abs=0.04)
    assert np.median(c3["dilution_step"]) == pytest.approx(math.log10(1 / 256), abs=0.08)
    # no free aminoketone: the trap reads its blank (<= the reporting limit)
    assert np.all(c0["ak_trap"] <= 1e-9) and np.all(c3["ak_trap"] <= 1e-9)


def test_c1_signature_exit_dominated():
    c1 = _pred("C1", n=600)
    # the DMP rate barely follows MGO, dilution acts as second order, the aminoketone is far above the limit
    assert np.median(c1["mgo_step"]) < 0.15
    assert np.median(c1["dilution_step"]) == pytest.approx(math.log10(1 / 16), abs=0.15)
    assert np.median(c1["ak_trap"]) > 0.5


def test_quadratic_onset_has_low_r2_and_linear_has_one():
    t = M.ZHOU_TIMES
    assert M.r2_linear(t, M.lag_shape(t, "linear", 0.0)[None, :])[0] == pytest.approx(1.0)
    assert M.r2_linear(t, (t ** 2)[None, :])[0] == pytest.approx(0.9195, abs=1e-3)


def test_job_a_verdicts_the_case_quotes():
    s = M.strecker_only_amine_exit()
    assert not s["falsified"]
    assert s["hofmann"]["rows"][0]["products_over_consumed"] == pytest.approx(0.181, abs=1e-3)
    k = M.kanzler_one_to_one()
    assert k["falsified"] and k["amine_dependent_min_mM"] > 55
    assert k["consistent_check_piloty_sized_sink"]["ala_consumed_pct_at_lower_bound"] < 5


@pytest.mark.slow
def test_onset_falsification_and_exit_edge():
    a = M.onset_falsification(np.random.default_rng(1), reps=4000)
    assert a["no_exit_falsified"]
    assert 0.01 < a["ks_min_per_min"] < 0.05


def test_designed_units_are_within_run_credit_units():
    m = M.noise_models(M.UNITS_C + M.UNITS_S, M.NoiseSpec(n_rep=2))
    for u, nm in m.items():
        assert nm.role == "credit" and nm.level == "within_run"
        assert {t.kind for t in nm.terms} == {"within_run_sd", "instrument_floor"}
    extra = M.misspecified_extra_sd(("mgo_step",), M.NoiseSpec(n_rep=2), 1.0)
    assert extra["mgo_step"] == pytest.approx(float(m["mgo_step"].terms[0].sd_array[0]))


def test_defensive_proposal_weights_integrate_to_one():
    rng = np.random.default_rng(3)
    x, lw = M._defensive_uniform(0.0, 1.0, 200000, rng)
    assert np.mean(np.exp(lw)) == pytest.approx(1.0, abs=0.02)
    x, lw = M._defensive_normal(0.01, 200000, rng)
    assert np.mean(np.exp(lw)) == pytest.approx(1.0, abs=0.02)


def test_clopper_pearson_bounds():
    lo, hi = M.clopper_pearson(0, 50)
    assert lo == 0.0 and hi == pytest.approx(0.0711, abs=1e-3)
    lo, hi = M.clopper_pearson(50, 50)
    assert hi == 1.0 and lo == pytest.approx(0.9289, abs=1e-3)
    assert pdn.wilson_interval(0, 50)[0] == pytest.approx(0.0, abs=1e-12)
