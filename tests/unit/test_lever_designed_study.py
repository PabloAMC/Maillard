"""Synthetic tests of the lever table's designed-study path (PREREG A7, R0 - R5), written before any P3 value
exists: scripts/generators/lever_designed_study.py and the R4 term in build_lever_table.py."""
import importlib.util
import itertools
import math
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / "generators" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


LT = _load("build_lever_table")
DS = _load("lever_designed_study")

EFFECT = {"ribose": 0.3, "thiamine": -0.5, "ph": 0.2}     # dex; pH per unit


def _log10_value(sugar, sulfur, ph, buffer="matrix", interaction=None):
    y = 1.0 + EFFECT["ribose"] * (sugar == "ribose") + EFFECT["thiamine"] * (sulfur == "thiamine") \
        + EFFECT["ph"] * (ph - 6.0)
    if buffer == "phosphate":
        y += 0.1
    if interaction:
        y += interaction(sugar, sulfur, ph)
    return y


def pots(lab="L1", cluster="P3", n=4, final_shift=1.0, lod=None, nd=lambda s, u, p, i: False, iso=False,
         interaction=None, with_pair=True, phos_slope_extra=0.0, odorant="methional"):
    out = []
    cells = [(s, u, p, "matrix") for s, u, p in itertools.product(("glucose", "ribose"), ("cysteine", "thiamine"),
                                                                    (6.0, 7.0))]
    if with_pair:
        cells += [("glucose", "cysteine", p, "phosphate") for p in (6.0, 7.0)]
    for s, u, p, b in cells:
        y = _log10_value(s, u, p, b, interaction)
        if b == "phosphate":
            y += phos_slope_extra * (p - 6.0)
        for i in range(n):
            gone = nd(s, u, p, i) and b == "matrix"
            out.append({"lab": lab, "cluster": cluster, "study": "P3", "odorant": odorant, "sugar": s,
                        "sulfur": u, "ph_set": p, "buffer": b, "final_ph": p - 0.4 + (final_shift - 1.0) * (p - 6.0),
                        "value": None if gone else 10 ** y, "lod": lod, "isotope_dilution": iso})
    return out


# ------------------------------------------------------------------------------------------------ R0
def test_pot_sd_from_declared_cv():
    assert DS.sigma_pot(True) == pytest.approx(math.sqrt(0.01 + math.log10(1.1) ** 2))
    assert DS.sigma_pot(False) == pytest.approx(math.sqrt(0.01 + math.log10(1.2) ** 2))


# ------------------------------------------------------------------------------------------------ R1
def test_corner_states():
    det = [{"value": 100.0}, {"value": 1000.0}]
    assert DS.corner_value(det) == pytest.approx({"state": "measured", "y": 2.5, "v": DS.sigma_pot(False) ** 2 / 2,
                                                  "n": 2, "censored_pots": 0})
    b = DS.corner_value([{"value": None, "lod": 10.0}, {"value": None, "lod": 1.0}])
    assert b["state"] == "bounded" and b["y"] == pytest.approx(1.0) and b["v"] == pytest.approx(
        DS.sigma_pot(False) ** 2 / 2)
    assert DS.corner_value([{"value": 5.0}, {"value": None, "lod": None}])["state"] == "excluded"


def test_censored_corner_is_the_tobit_ml_mean_with_observed_information():
    ps = [{"value": 100.0}, {"value": 50.0}, {"value": None, "lod": 30.0}]
    c = DS.corner_value(ps)
    sig = DS.sigma_pot(False)
    det, lim = np.log10([100.0, 50.0]), math.log10(30.0)
    grid = np.linspace(0.5, 3.0, 200001)
    ll = stats.norm.logpdf(det[:, None], grid, sig).sum(0) + stats.norm.logcdf(lim, grid, sig)
    assert c["state"] == "censored" and c["y"] == pytest.approx(grid[np.argmax(ll)], abs=1e-4)
    assert c["y"] < det.mean()            # the not-detected pot pulls the mean down
    assert sig ** 2 / 3 < c["v"] < sig ** 2 / 2   # a censored pot carries less information than a measured one


def test_a4_contrast_is_the_mean_over_corners_of_the_designed_effects():
    ps = [p for p in pots() if p["buffer"] == "matrix"]
    for lever, eff in (("sugar identity", 0.3), ("sulfur source", -0.5), ("pH", 0.2)):
        c = DS.a4_contrast(ps, lever)
        assert c["enters"] and not c["lower_bound"] and c["y"] == pytest.approx(eff)
        assert c["v"] == pytest.approx(8 * DS.sigma_pot(False) ** 2 / 4 / 16)
        assert all(k["y"] == pytest.approx(eff) for k in c["corners"])


def test_bounded_lower_arm_enters_as_a_lower_bound_but_bounded_higher_arm_does_not():
    # thiamine is the lower arm (-0.5 dex); bound 3 of its 4 corners
    lower = [p for p in pots(lod=10 ** 1.05, nd=lambda s, u, p, i: u == "thiamine" and not (s == "glucose" and p == 6.0))
             if p["buffer"] == "matrix"]
    c = DS.a4_contrast(lower, "sulfur source")
    assert c["enters"] and c["lower_bound"] and c["bounded_corners"] == 3
    assert c["y"] == pytest.approx((0.5 + 3 * 1.05) / 4 - 1.25)   # biased toward zero (true -0.5)
    # cysteine is the higher arm; bound 3 of its corners at a limit far above the values
    higher = [p for p in pots(lod=1e4, nd=lambda s, u, p, i: u == "cysteine" and not (s == "glucose" and p == 6.0))
              if p["buffer"] == "matrix"]
    c = DS.a4_contrast(higher, "sulfur source")
    assert not c["enters"] and "higher arm" in c["reason"]


def test_corner_bounded_in_both_arms_needs_one_limit():
    def nd(s, u, p, i):
        return s == "glucose" and p == 6.0
    ps = [dict(p, lod=(1.0 if p["sulfur"] == "cysteine" else 2.0)) for p in pots(nd=nd) if p["buffer"] == "matrix"]
    c = DS.a4_contrast(ps, "sulfur source")
    assert not c["enters"] and "different limits" in c["reason"]
    same = [dict(p, lod=1.0) for p in ps]
    c = DS.a4_contrast(same, "sulfur source")
    assert c["enters"] and sum(k["both_bounded"] for k in c["corners"]) == 1
    assert [k["y"] for k in c["corners"]].count(None) == 1      # no corner contrast for R3


# ------------------------------------------------------------------------------------------------ R2
def test_ph_contrast_enters_only_when_the_set_levels_held():
    rows, rep = DS.designed_rows(pots(), resolve=str)
    ph = [r for r in rows if r["lever_class"] == "pH"]
    assert rep["ph_held"]["L1"]["held"] and len(ph) == 1
    assert ph[0]["band"] == "pH 6-7.5, initial or label" and ph[0]["y"] == pytest.approx(0.2)
    rows, rep = DS.designed_rows(pots(final_shift=0.2), resolve=str)     # final pH 7.0 vs 6.0 differ by 0.2
    assert not rep["ph_held"]["L1"]["held"] and not [r for r in rows if r["lever_class"] == "pH"]
    entry = [c for c in rep["contrasts"] if c["lever"] == "pH"][0]
    assert not entry["enters"] and entry["y"] == pytest.approx(0.2)     # reported, enters no estimate
    missing = pots()
    missing[0]["final_ph"] = None
    assert not DS.ph_held(missing)["held"]


# ------------------------------------------------------------------------------------------------ R3
def test_corner_check_flags_a_reversing_corner_only():
    ps = [p for p in pots() if p["buffer"] == "matrix"]
    assert not DS.corner_check([DS.a4_contrast(ps, "sugar identity")])["flag"]

    def reverse(s, u, p):              # at thiamine, pH 7 ribose loses 2.0 dex: one corner reverses strongly
        return -2.0 if (s == "ribose" and u == "thiamine" and p == 7.0) else 0.0
    rev = [p for p in pots(interaction=reverse, iso=True) if p["buffer"] == "matrix"]

    def big(s, u, p):                  # the main effect stays positive and clear: +1.2 at the three other corners
        return 0.9 if (s == "ribose" and not (u == "thiamine" and p == 7.0)) else 0.0
    both = [p for p in pots(interaction=lambda s, u, p: reverse(s, u, p) + big(s, u, p), iso=True)
            if p["buffer"] == "matrix"]
    c = DS.a4_contrast(both, "sugar identity")
    assert c["y"] > 0
    chk = DS.corner_check([c])
    assert chk["flag"] and chk["label"] == "interaction-dominated"
    # the reversal alone drags the mean to -0.2 dex while three corners stay at +0.3: flagged the other way
    c = DS.a4_contrast(rev, "sugar identity")
    assert c["y"] == pytest.approx(-0.2) and DS.corner_check([c])["flag"]
    # a main effect whose 90 % interval includes zero is never flagged
    weak = [dict(p, isotope_dilution=False) for p in pots(interaction=lambda s, u, p: -0.3 * (s == "ribose")
                                                          + (0.4 if (s == "ribose" and p == 7.0) else 0.0) - 0.4 *
                                                          (s == "ribose" and p == 6.0))
            if p["buffer"] == "matrix"]
    c = DS.a4_contrast(weak, "sugar identity")
    assert c["y"] == pytest.approx(0.0) and not DS.corner_check([c])["flag"]


# ------------------------------------------------------------------------------------------------ R5
def test_buffered_pair_D_and_its_reading_against_the_margin():
    bp = DS.buffered_pair(pots(iso=True))
    assert bp["enters"] and bp["D"] == pytest.approx(0.0, abs=1e-12)
    assert bp["v"] == pytest.approx(4 * DS.sigma_pot(True) ** 2 / 4)
    assert bp["buffer_contrast"]["y"] == pytest.approx(-0.1)          # matrix-buffered minus phosphate
    # one laboratory: at 0.108 dex per pot, D = 0 reads "same"; at 0.128 it cannot (the stated power)
    assert DS.buffered_pair_reading({"L1": bp})["reading"] == "same"
    assert DS.buffered_pair_reading({"L1": DS.buffered_pair(pots())})["reading"] == "unresolved"
    big = DS.buffered_pair(pots(phos_slope_extra=0.8))
    assert big["D"] == pytest.approx(0.8) and DS.buffered_pair_reading({"L1": big})["reading"] == "different"


def test_buffered_pair_enters_the_buffer_cell_against_matrix_buffered():
    rows, rep = DS.designed_rows(pots(), resolve=str)
    buf = [r for r in rows if r["lever_class"] == "buffer"]
    assert len(buf) == 1 and buf[0]["pairs"][0]["level"] == DS.MATRIX_BUFFERED
    assert buf[0]["pairs"][0]["against"] == "phosphate"
    rows, _ = DS.designed_rows(pots(with_pair=False), resolve=str)
    assert not [r for r in rows if r["lever_class"] == "buffer"]


# ------------------------------------------------------------------------------------------------ R4
def _literature_row(lab, y, v=0.01):
    return {"source": f"paper_{lab}", "contrast": "c1", "compound": "methional", "lever_class": "pH", "lab": lab,
            "partition": "test", "focal": True, "peak_area": False, "matrix": ["water"], "temps": [140.0],
            "kind": "numeric", "y": y, "v": v, "band": DS.PH_BAND, "non_monotone": False, "range": [5, 7],
            "censored_arms": 0, "alt": {}, "turning_point": False, "buffer_state": "buffered"}


def test_cluster_counts_as_one_laboratory_and_lolo_leaves_it_out_whole():
    rows = []
    for lab in ("L1", "L2"):
        rows += DS.designed_rows(pots(lab=lab), resolve=str)[0]
    ph = [r for r in rows if r["lever_class"] == "pH"]
    lit = [_literature_row(f"lab{i}", 0.2 + 0.05 * i) for i in range(3)]
    cells = {(c["lever_class"], c["band"]): c for c in LT.cell_table(ph + lit)}
    cell = cells[("pH", DS.PH_BAND)]
    assert cell["labs"] == 4 and cell["design_labs"] == ["L1", "L2"]        # 3 literature + ONE cluster
    assert cell["lolo"]["labs"] == 4
    assert cell["fit"].get("tau_w_fitted") and set(cell["fit"]["tau_w_sensitivity"]) == {"0.1", "0.5"}
    assert cell["pea_answer_shown"] is False          # no outside laboratory has cooked in an isolate
    only = {(c["lever_class"], c["band"]): c for c in LT.cell_table(ph)}[("pH", DS.PH_BAND)]
    assert only["labs"] == 1 and only["reading"] == "too few laboratories"
    sug = {c["lever_class"]: c for c in LT.cell_table(rows)}["sugar identity"]
    assert sug["fit"]["levels"]["ribose"]["labs"] == 1 and sug["fit"]["levels"]["ribose"]["reading"] == \
        "too few laboratories"


def test_new_laboratory_interval_excludes_w():
    rows = []
    for lab in ("L1", "L2"):
        rows += [r for r in DS.designed_rows(pots(lab=lab), resolve=str)[0] if r["lever_class"] == "pH"]
    rows += [_literature_row(f"lab{i}", 0.2 + 0.05 * i) for i in range(3)]
    narrow = LT.numeric_fit(rows, LT.PRIOR_SCALE, tau_w_scale=0.25)
    # adding w as a new-laboratory draw would widen the interval; it is not one
    y = np.array([r["y"] for r in rows]); S = np.diag([r["v"] for r in rows]); X = np.ones((len(rows), 1))
    comps = [LT._incidence([r["lab"] for r in rows]), np.eye(len(rows)), LT._within_cluster(rows)]
    fit = LT.fit_gaussian(y, S, X, comps, [LT.PRIOR_SCALE, LT.PRIOR_SCALE, 0.25])
    assert LT.summarise(fit, 0, [0, 1])["new_lab_pi90"] == pytest.approx(narrow["new_lab_pi90"])
    assert LT.summarise(fit, 0, [0, 1, 2])["new_lab_pi90_width"] > narrow["new_lab_pi90_width"]


def test_protocol_report_has_no_new_laboratory_interval():
    _, rep = DS.designed_rows(pots(), resolve=str)
    pr = rep["protocol"]["methional | sugar identity"]
    assert pr["estimate"] == pytest.approx(0.3) and pr["new_lab_interval"] is None and pr["labs"] == 1


def test_fit_gaussian_scalar_scale_equals_the_repeated_list():
    rng = np.random.default_rng(0)
    y = rng.normal(0.2, 0.2, 6); S = np.diag(np.full(6, 0.01)); X = np.ones((6, 1))
    comps = [LT._incidence(["a", "a", "b", "b", "c", "c"]), np.eye(6)]
    a, b = LT.fit_gaussian(y, S, X, comps, 0.5), LT.fit_gaussian(y, S, X, comps, [0.5, 0.5])
    assert np.allclose(a["w"], b["w"]) and np.allclose(a["means"], b["means"])
