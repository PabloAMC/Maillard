"""Forced-input Strecker layer (roadmap step 2; scripts/generators/feasibility_strecker_forced_input.py): the module
follows docs/guides/A2_MODULE_CONTRACT.md sec. 3 -- balanced reactions, no colliding keys, exact calibrations, a
quadrature that matches the analytic fed pot, oxygen off = argon, and a noise-free truth that passes itself."""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

from src.kinetic_core import fit_rows

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "generators"))
import feasibility_strecker_forced_input as fs  # noqa: E402


@pytest.fixture(scope="module")
def data():
    return fs.load_data()


def _central(n=1):
    return {c.name.split(".", 1)[1]: np.full(n, float(c.default)) for c in fs.DRAWN}


@pytest.mark.parametrize("reaction", fs.REACTIONS, ids=lambda r: r.key)
def test_every_reaction_is_carbon_nitrogen_and_sulfur_balanced(reaction):
    for element in ("carbon", "nitrogen", "sulfur"):
        left, right = reaction.atom_balance(element, fs.LOOKUP)
        assert left == right, (reaction.key, element, left, right)


def test_new_species_keys_do_not_collide_with_engine_keys():
    engine = fit_rows.species_vocabulary()
    assert not {s.key for s in fs.SPECIES} & engine
    assert {"MGO", "GO", "TDG", "G", "FRAG_C"} <= engine


def test_every_reaction_has_a_rate_constant_and_routes_are_config_not_globals():
    cfg = fs.config("DC_ARP")
    params = fs.calibrate(_central(), fs.load_data())
    k = fs.rate_constants(params, {"cfg": cfg, "o2": np.ones(1)})
    assert set(k) == {r.key for r in fs.REACTIONS}
    off = fs.rate_constants(params, {"cfg": fs.config("DC"), "o2": np.ones(1)})
    assert all(np.all(off[key] == 0) for key in k if key.startswith("r_arp_"))
    assert fs.config("DC_ARP") == cfg          # config() builds a fresh mapping; nothing global changed


def test_fed_pot_quadrature_matches_the_analytic_second_order_solution():
    kappa, t = np.array([7e-4]), 60.0
    d, aa = fs.fed_donor_courses(kappa, np.array([0.0]))       # nu = 0: amine constant
    integ = fs.cum_trapz(d * aa)[0, fs.idx(t)]
    exact = fs.AA0 * math.log(1.0 + kappa[0] * fs.D0 * t) / kappa[0]
    assert integ == pytest.approx(exact, rel=1e-4)


def test_calibration_reproduces_the_fed_argon_rows_exactly(data):
    th = _central()
    cal = fs.calibrate(th, data)
    for d in fs.DONORS:
        dd, aa = fs.fed_donor_courses(cal["kappa"][d], th["nu"])
        integ = fs.cum_trapz(dd * aa)[0, fs.idx(60)]
        lossy = fs.cum_with_loss(dd * aa, cal["k_loss"])[0, fs.idx(60)]
        k, phi = cal["k"][d][0], cal["phi"][d][0]
        assert (1 - phi) * k * lossy == pytest.approx(data.fed[d]["argon"][0], rel=1e-9)
        assert phi * k * integ == pytest.approx(data.fed[d]["argon"][1], rel=1e-9)


def test_the_glyoxal_oxygen_hold_out_is_the_dry_review_number(data):
    s = fs.score(_central(), data, fs.config("DC_ARP"))
    # review (EX-042): omega_C2/C3 from MGO predicts GO's PAA air / argon 6.1 vs 5.5 measured (0.04 dex)
    assert float(s["O2_GO"][0]) == pytest.approx(math.log10(6.1 / 5.5), abs=0.02)


def test_oxygen_off_leaves_only_the_argon_channels(data):
    th = _central()
    cal = fs.calibrate(th, data)
    courses = fs.interpolate(data.forced_inputs)
    g0 = fs.glucose_pot(th, cal, courses, fs.config("DC_ARP"), o2=np.zeros(1))
    g_no = fs.glucose_pot(th, cal, courses, fs.config("DC_ARP", **{"strecker.route.o2_acid": False}), o2=np.zeros(1))
    np.testing.assert_allclose(g0["PAA"], g_no["PAA"])


def test_a_noise_free_truth_passes_its_own_signatures(data):
    th = _central()
    cfg = fs.config("DC_ARP")
    truth = fs.truth_rows(th, data, cfg)
    s = fs.score(th, truth, cfg)
    for key in fs.JOINT:
        assert float(s[key][0]) < 0.02, key


# ---- the blind test on hofmann2000b Fig. 1 (prereg sec. 12) ----

@pytest.fixture(scope="module")
def fig1_meas():
    return fs.fig1_measured(fs.load_fig1())


def test_fig1_table_is_keyed_two_pots_28_rows():
    fig1 = fs.load_fig1()
    assert sum(len(t) for pot in fig1.values() for t in pot.values()) == 28
    assert set(fig1["glc"]) == set(fig1["arp"]) == {"PHEAL", "PHEAC"}


def test_rule_6_drops_exactly_the_hidden_10_min_glucose_points(fig1_meas):
    dropped = {(s, t) for s in ("G5p", "G7") for t, m in fig1_meas[s].items() if m["dropped"]}
    assert dropped == {("G5p", 10.0)}
    assert fig1_meas["G5p"][10.0]["read_error"] > fs.BLIND_TOL["G5p"] / 2


def test_fig1_in_sample_numbers_are_the_keyed_ones(fig1_meas):
    assert 10.0 ** fig1_meas["G6a_p"] == pytest.approx(3.709 / 2.804, rel=2e-3)
    assert 60.0 < fig1_meas["G6b_p"] < 120.0


def _score(ok, early=False):
    return {"blind_points_within_tol": ok, "early_deficit": early}


def test_blind_outcome_rule():
    assert fs.blind_outcome([_score(True)] * 4 + [_score(False)]) == "PASS"           # every 4-subset has >= 3
    assert fs.blind_outcome([_score(True)] * 3 + [_score(False)] * 2) == "INCONCLUSIVE"
    assert fs.blind_outcome([_score(False, True)] * 5) == "EARLY DEFICIT REPRODUCED"
    assert fs.blind_outcome([_score(False, True)] * 4 + [_score(False)]) == "INCONCLUSIVE"


def test_a_model_scored_on_its_own_prediction_passes_every_blind_point(data, fig1_meas):
    theta = {k: float(v[0]) for k, v in _central().items()}
    pred = fs.blind_predict(theta, data)
    tot = pred["glc_PA"] + pred["glc_PAA"]
    own = {**fig1_meas,
           "G5p": {t: dict(m, measured=float(np.log10(tot[fs.idx(t)] / tot[fs.idx(120)]))) for t, m in fig1_meas["G5p"].items()},
           "G7": {t: dict(m, measured=float(np.log10(pred["arp_PA"][fs.idx(t)] / pred["arp_PA"][fs.idx(120)])))
                  for t, m in fig1_meas["G7"].items()}}
    s = fs.blind_score(pred, own)
    assert s["blind_points_within_tol"] and not s["early_deficit"]
    assert max(abs(v) for v in s["G5p"].values()) < 1e-9
