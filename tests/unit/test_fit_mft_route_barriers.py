"""
tests/unit/test_fit_mft_route_barriers.py

E6 + E11 (results/validation/kinetic_core_mft_route_barriers_prereg.md): each pre-registered member's engine state is
ONE Config differing from the shipped one in the B62 hook only, carrying exactly the declared barriers; the free set is
B60's 23 (+ log10 k_dpo_nf and the tied DDP offset for the A members); the coordinates round-trip; ``applied`` puts
the shipped engine back. No integration is run.
"""
from __future__ import annotations

import sys

import numpy as np
import pytest

from src.kinetic_core import config as engine_config
from src.kinetic_core import parameters_sulfur as ps

LUMP = 64.0822869313293
B62 = {"k_tdp_fur": 147.0, "k_fur_decay": 46.5}
EXPECTED = {
    "N1": {**B62},
    "A64": {**B62, "k_dpo_nf": 83.8, "k_nf_mft": LUMP},
    "A100": {**B62, "k_dpo_nf": 83.8, "k_nf_mft": 100.0},
    "A133": {**B62, "k_dpo_nf": 83.8, "k_nf_mft": 133.0},
    "A100-s1": {**B62, "k_dpo_nf": 83.8, "k_nf_mft": 100.0},
    "A133-93": {**B62, "k_dpo_nf": 93.0, "k_nf_mft": 133.0},
    "A64-93": {**B62, "k_dpo_nf": 93.0, "k_nf_mft": LUMP},
}


@pytest.fixture(scope="module")
def gen():
    from tests.support import wave_generator

    sys.modules.pop("fit_mft_route_barriers", None)
    with wave_generator("fit_mft_route_barriers") as module:
        yield module
    sys.modules.pop("fit_mft_route_barriers", None)


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_member_config_is_the_shipped_config_with_the_declared_barriers(gen, name):
    cfg = gen.member_config(name)
    diff = cfg.diff()
    assert set(diff) == {"sulfur.furfural_barriers"}, diff
    hook = cfg["sulfur.furfural_barriers"]
    assert hook.enabled
    assert dict(hook.barriers_kj_mol) == pytest.approx(EXPECTED[name])
    with engine_config.applied(cfg):
        assert ps.measured_ea_overrides() == {**ps.MEASURED_EA_OVERRIDES, **EXPECTED[name]}
    assert ps.FURFURAL_BARRIERS == ps.FurfuralBarriers(), "applied() must restore the shipped hook"
    assert engine_config.Config.from_globals() == engine_config.Config.shipped()


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_member_free_set_and_start(gen, name):
    import generate_kinetic_core_b60_fit as B60

    coords = gen.coordinates(name)
    keys = [c[0] for c in coords]
    extra = [] if name == "N1" else ["k_dpo_nf", "offset_ddp_route"]
    assert keys == list(B60.FREE_KEYS) + extra
    assert len(keys) == 23 + len(extra)
    m = gen.member(name)
    assert m.start == (1 if name.endswith("-s1") else 0)
    theta = gen.start_theta(name)
    lower, upper = np.array([c[2] for c in coords]), np.array([c[3] for c in coords])
    assert np.all(theta >= lower) and np.all(theta <= upper)
    if m.start == 0:
        assert np.array_equal(gen.to_full(theta, name), gen.incumbent_vector())


def test_ddp_offset_moves_both_ddp_keys_together_within_their_bounds(gen):
    import generate_kinetic_core_b60_fit as B60

    base = gen.incumbent_vector()
    lo, hi = B60.full_bounds()
    label, kind, olo, ohi = gen.coordinates("A100")[-1]
    assert label == "offset_ddp_route" and kind == "log10k"
    theta = gen.theta_of(base, "A100")
    assert theta[-1] == 0.0
    for offset in (olo, 2.0, ohi):
        theta[-1] = offset
        x = gen.to_full(theta, "A100")
        for key in ("k_ddp_mft", "k_ddp_mft_hs"):
            i = B60.ALL_KEYS.index(key)
            assert x[i] == pytest.approx(base[i] + offset)
            assert lo[i] - 1e-12 <= x[i] <= hi[i] + 1e-12
        assert x[B60.ALL_KEYS.index("k_dpo_ddp")] == base[B60.ALL_KEYS.index("k_dpo_ddp")]
        assert gen.theta_of(x, "A100") == pytest.approx(theta)


def test_unregistered_members_are_refused(gen):
    for bad in ("N1-s1", "A120", "A100-s2", "B62"):
        with pytest.raises(ValueError):
            gen.member(bad)


def test_ledger_allows_only_named_rows_in_their_shipped_direction(gen):
    rows = {"zhou_pH7_MFT": {"raw": -2.0, "units": "dex"}, "kang_100C_MFT": {"raw": 0.9, "units": "dex"},
            "zhou_pH7_FFT": {"raw": 1.0, "units": "dex"}}
    b60 = {"zhou_pH7_MFT": {"raw": -1.5}, "kang_100C_MFT": {"raw": 0.5}, "zhou_pH7_FFT": {"raw": -0.4}}
    out = gen.ledger(rows, b60)
    assert set(out["allowed_by_ledger"]) == {"zhou_pH7_MFT"}
    assert set(out["violations"]) == {"kang_100C_MFT", "zhou_pH7_FFT"}  # FFT flipped sign: not its shipped direction
    assert out["holds"] is False
