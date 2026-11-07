"""Unit tests for the M4 reduced forced-input mini-model (src/m4_reduced_power.py; case M4 sec. 5.4).

Label conservation, mass balance, each rival's structural prediction, the hand-computed label table, the claim's FM4
rescale, job (a)'s arithmetic, and a synthetic recovery check. Everything here is synthetic: no measured value of a
decisive unit is read.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from src import m4_reduced_power as M

# by hand from the maps (module docstring): (K_u, K_s, L_u, L_s, C2 + C3)
HAND = {
    "L1_M": (1, 0, 0, 1, 0.75), "L1_ring": (0, 0, 0, 1, 0.5), "L1_acyl": (1, 0, 0, 0, 0.5),
    "L6_M": (0, 1, 1, 0, 0.75), "L6_ring": (0, 0, 1, 0, 0.5), "L6_acyl": (0, 1, 0, 0, 0.5),
    "L2_M": (1, 1, 1, 1, 0.75), "L2_ring": (1, 1, 0, 1, 0.75), "L2_acyl": (1, 0, 1, 0, 0.5),
    "L12_M1": (0, 1, 1, 0, 0.0), "L12_M2": (1, 0, 0, 1, 0.5), "CAM": (0, 0, 0, 0, 0.5),
}


def _draws(rival, n=6, seed=3, form=M.PRIMARY):
    th, cnt = M.draw_set(rival, form, n, np.random.default_rng(seed))
    assert M.n_of(th) == n
    return th


def _pot(name):
    return M.POT_BY_NAME[name]


@pytest.mark.parametrize("read", sorted(HAND))
def test_label_table_matches_hand(read):
    tb = M.label_table(read)
    assert (tb["K_u"], tb["K_s"], tb["L_u"], tb["L_s"], tb["C23"]) == HAND[read]


def test_label_maps_are_permutations_and_complementary():
    # each family's C5 keeps five distinct positions; the lost carbon is the sixth; the swap is p -> 7 - p
    for u, s, lost_u in ((M.KEPT_U, M.KEPT_S, 6), (M.LOST_U, M.LOST_S, 1)):
        assert len(set(u)) == 5 and set(range(1, 7)) - set(u) == {lost_u}
        assert tuple(7 - p for p in u) == s


@pytest.mark.parametrize("rival", M.RIVALS)
def test_mass_carbon_cysteine_and_label_balances(rival):
    th = _draws(rival)
    for av in ("a1", "a2", "a2cys", "a3", "a3m"):
        for pot in M.POTS[av][::2]:
            b = M.balances(pot, M.simulate(th, rival, M.PRIMARY, pot, full=True))
            for key, v in b.items():
                assert v < 1e-10, (rival, pot.name, key, v)


def test_r0_mixed_share_is_the_binomial_and_r0_is_below_lod():
    th = _draws("R0")
    p = M.predict(th, "R0", M.PRIMARY, arms=("a1",))
    cam = M.GROUPS["a1:CAM"]
    assert np.allclose(p.share_a[:, cam], 0.5, atol=1e-12) and np.allclose(p.share_b[:, cam], 0.0)
    lev = M.GROUPS["a1:MFT"]
    assert np.all(p.level[:, lev] < 1e-3 * M.LOD_MM["MFT"])        # R0's sink-free bound is > 3 dex under the LOD


def test_d3_arithmetic():
    a = M.d3_arithmetic()
    assert a["mixed_share_c2c3"] == pytest.approx(0.5)
    assert a["max_route_share_f"] == pytest.approx((0.01 + 1.96 * 0.01) / 0.5)
    assert a["z_full_route"] == pytest.approx(49.0)


def _arm(rival, arm, form=M.PRIMARY, n=6, seed=3):
    th = _draws(rival, n, seed, form)
    return {pot.name: M.simulate(th, rival, form, pot) for pot in M.POTS[arm]}


def _mft(r, t):
    return r[("MFTK", t)] + r[("MFTL", t)] + r[("MFTD", t)] + r[("MFTC", t)]


def test_r1_structure():
    a1 = _arm("R1", "a1")
    for name, r in a1.items():
        assert np.all(r[("MFTL", 20.0)] == 0.0) and np.all(r[("MFTD", 20.0)] == 0.0)
        assert np.all(r[("MFTK", 20.0)] > 0.0)                         # MFT under air and N2
    a2 = _arm("R1", "a2")
    for name, r in a2.items():
        if "_GOS_" in name:
            assert np.all(_mft(r, 20.0) == 0.0) and np.all(r[("DPO1", 5.0)] == 0.0)   # fed glucosone: no C5
        else:
            assert np.all(_mft(r, 20.0) > 0.0)                         # fed 1-DG -> 1-DPO -> MFT, no amino acid needed
    p = M.predict(_draws("R1"), "R1", M.PRIMARY, arms=("a1",))
    i = M.GROUPS["a1:L1_M"]
    assert np.allclose(p.share_a[:, i], 1.0) and np.allclose(p.share_b[:, i], -1.0)    # labelled methyl, 1 - e
    i = M.GROUPS["a1:L2_ring"]
    assert np.allclose(p.share_a[:, i], 1.0) and np.allclose(p.share_b[:, i], 0.0)     # C-2 always in the ring


def test_r5_structure():
    a1 = _arm("R5", "a1")
    for name, r in a1.items():
        assert np.all(r[("MFTK", 20.0)] == 0.0)
        if name.endswith("_N2"):
            assert np.all(_mft(r, 20.0) == 0.0)                        # glucosone is oxidative: no MFT under N2
        else:
            assert np.all(r[("MFTL", 20.0)] > 0.0)
    a2 = _arm("R5", "a2")
    for name, r in a2.items():
        if "_GOS_" in name:
            assert np.all(r[("DPO1", 5.0)] > 0.0) and np.all(r[("FA", 20.0)] > 0.0) and np.all(_mft(r, 20.0) > 0.0)
        else:
            assert np.all(_mft(r, 20.0) == 0.0)                        # fed 1-DG has no C5 exit in R5
    p = M.predict(_draws("R5"), "R5", M.PRIMARY, arms=("a1",))
    air = [i for i in M.GROUPS["a1:L2_ring"] if M.OBS[i].pot.endswith("_air")]
    assert np.allclose(p.share_a[:, air], 0.0) and np.allclose(p.share_b[:, air], 1.0)   # C-2 is the methyl (e = 0)


def test_r2_structure():
    a1 = _arm("R2", "a1")
    for r in a1.values():
        assert np.all(r[("NFK", 20.0)] + r[("NFL", 20.0)] + r[("NFC", 20.0)] == 0.0)   # MFT bypasses NF
        assert np.all(r[("MFTD", 20.0)] > 0.0)
    a2 = _arm("R2", "a2")
    for name, r in a2.items():
        assert np.all(r[("DPO1", 5.0)] == 0.0)
        if "_DG1_" in name and "gly33" in name:
            assert np.all(_mft(r, 20.0) > 0.0) and np.all(r[("DDQ", 5.0)] > 0.0)
        else:
            assert np.all(_mft(r, 20.0) == 0.0)                        # Na2S, no amino acid / glucosone: nothing


def test_r4p_structure():
    a3 = _arm("R4P", "a3")
    for name, r in a3.items():
        if "_fed_" in name:
            assert np.all(_mft(r, 60.0) == 0.0)                        # the pool is a sink: no C5 from fed Amadori
            assert np.all(r[("CYS", 60.0)] > 0.0)                      # ... but it returns cysteine
        else:
            assert np.all(_mft(r, 60.0) > 0.0)                         # enolisation x catalysis feeds 1-DPO


def test_claim_rescale_hits_f_and_is_inside_its_bands():
    th, cnt = M.draw_set("CLAIM", M.PRIMARY, 16, np.random.default_rng(9))
    f = M.ref_share(th, M.PRIMARY)
    assert np.max(np.abs(f - th["f_ref"])) < 1e-3
    band = {c.name: (c.lo, c.hi) for c in M.coords("CLAIM")}
    eo = th["log10_odds_odg"] + th["_dlog_odg"]
    assert np.all((eo >= band["log10_odds_odg"][0] - 1e-9) & (eo <= band["log10_odds_odg"][1] + 1e-9))
    # the reference pot's MFT carries the same R5 share as its 1-DPO (one shared downstream)
    r = M.simulate(th, "CLAIM", M.PRIMARY, M.REF_POT)
    share_mft = r[("MFTL", 20.0)] / (r[("MFTK", 20.0)] + r[("MFTL", 20.0)])
    assert np.max(np.abs(share_mft - th["f_ref"])) < 0.05


def test_claim_proposal_weights_are_prior_over_proposal():
    th = M.draw("CLAIM", M.PRIMARY, 4000, np.random.default_rng(1), proposal=True)
    w = np.exp(th["_logw"])
    # E_q[p / q] = 1 and the weighted f is uniform on [0.001, 0.999]: mean 0.5
    assert np.mean(w) == pytest.approx(1.0, abs=0.05)
    assert np.sum(w * th["f_ref"]) / np.sum(w) == pytest.approx(0.5, abs=0.03)


def test_pins_and_form_touches():
    f = M.Form("x", pins=(("log10_cat", 1.0),))
    assert all(c.lo == c.hi == 1.0 for c in M.coords("R4P", f) if c.name == "log10_cat")
    assert M.form_touches(f, "R4P") and not M.form_touches(f, "R1") and not M.form_touches(f, "R0")
    assert M.form_touches(M.Form("o", o2_budget=True), "R5") and not M.form_touches(M.Form("o", o2_budget=True), "R0")


def test_o2_budget_factor():
    pot = _pot("a1_145_pH5_air")
    assert M.o2_factor(pot, M.PRIMARY) == 1.0
    assert M.o2_factor(pot, M.Form("b", o2_budget=True)) == pytest.approx(1 - 33.0 / (4 * M.O2_MM_HALF_FILL))
    assert M.o2_factor(_pot("a1_145_pH5_N2"), M.PRIMARY) == 0.0


def test_censoring_and_share_readability():
    th = _draws("R5", 4)
    p = M.predict(th, "R5", M.PRIMARY, arms=("a1",))
    d = M.synthesize(p, 0, 3, M.PRIMARY, np.random.default_rng(2))
    for i in M.GROUPS["a1:MFT"]:
        if M.OBS[i].pot.endswith("_N2"):
            assert d.censored[i]                                      # R5 makes no MFT under N2
    for i in M.GROUPS["a1:L1_M"]:
        if M.OBS[i].pot.endswith("_N2"):
            assert not d.read[i]


def test_synthetic_recovery_r1_over_r5_and_claim_f():
    rng = np.random.default_rng(4)
    truth = M.predict(_draws("R1", 2, 11), "R1", M.PRIMARY, arms=("a1",))
    sc = {rv: M.predict(M.draw_set(rv, M.PRIMARY, 128, np.random.default_rng(5), proposal=True)[0], rv,
                        M.PRIMARY, arms=("a1",)) for rv in ("R1", "R5")}
    d = M.synthesize(truth, 0, 6, M.PRIMARY, rng)
    lp = {}
    for rv, p in sc.items():
        gl = M.group_loglik(p, d, M.PRIMARY, M.DESIGNS["arm1"])
        lp[rv] = M.lppd(M.design_loglik(gl, M.DESIGNS["arm1"]), p.logw)["lppd"]
    assert lp["R1"] - lp["R5"] > M.THRESHOLD_NATS
    # the claim recovers a mixture's f from arm 1's labels
    th_t, _ = M.draw_set("CLAIM", M.PRIMARY, 40, np.random.default_rng(12))
    j = int(np.argmin(np.abs(th_t["f_ref"] - 0.3)))
    pt = M.predict(M.take(th_t, [j]), "CLAIM", M.PRIMARY, arms=("a1",))
    d = M.synthesize(pt, 0, 6, M.PRIMARY, rng)
    th_s, _ = M.draw_set("CLAIM", M.PRIMARY, 512, np.random.default_rng(13), proposal=True)
    ps = M.predict(th_s, "CLAIM", M.PRIMARY, arms=("a1",))
    gl = M.group_loglik(ps, d, M.PRIMARY, M.DESIGNS["arm1"])
    ll = M.design_loglik(gl, M.DESIGNS["arm1"]) + ps.logw
    w = np.exp(ll - ll.max())
    f_post = float(np.sum(w * th_s["f_ref"]) / np.sum(w))
    assert abs(f_post - th_t["f_ref"][j]) < 0.15


def test_shipped_constants_match_the_engine():
    from src.kinetic_core import engine, sulfur
    p = engine.core_parameters(engine.SULFUR)
    k145 = sulfur.sulfur_rate_constants_at(p, M.T145, 5.0)
    k100 = sulfur.sulfur_rate_constants_at(p, M.T100, 5.0)
    assert k145["r_glc_c2c3"] == pytest.approx(M.K_GLC_HA_145, rel=1e-12)
    assert p["k_glc_ha"].ea_kj_mol == pytest.approx(M.EA_GLC_HA, rel=1e-12)
    assert k100["r_fru_odg"] == pytest.approx(M.K_FRU_ODG_100, rel=1e-12)
    assert k100["r_glc_tdg"] == pytest.approx(M.K_GLC_TDG_100, rel=1e-12)
    assert k100["r_odg_aa"] == pytest.approx(M.K_ODG_AA_100, rel=1e-12)
    ea = M.R_GAS * math.log(k145["r_odg_aa"] / k100["r_odg_aa"]) / (1 / M.T100 - 1 / M.T145)
    assert ea == pytest.approx(M.EA_ODG_AA, abs=0.1)
    assert math.log10(p["k_ama_g"].k_ref) == pytest.approx(M.LOG10_K_AMA_G_100, abs=1e-12)
    assert p["k_ama_g"].ea_kj_mol == pytest.approx(M.EA_AMA_G)
    assert math.log10(p["k_g_go"].k_ref) == pytest.approx(M.LOG10_K_G_GO_100, abs=1e-12)
    assert p["k_g_go"].ea_kj_mol == pytest.approx(M.EA_G_GO)


def test_heatings_count():
    assert M.heatings("arm1_levels_only", 3) == 3 * 8
    assert M.heatings("arm3", 6) == 6 * 24
