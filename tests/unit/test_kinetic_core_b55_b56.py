"""B55 / B56 (2026-10-28): the single-temperature Ea prior and the structural-discrepancy term
(results/validation/kinetic_core_b55_b56_prereg.md)."""
from __future__ import annotations

import json

import numpy as np
import pytest

from src import data_paths
from src.kinetic_core import ea_class_prior as E
from src.kinetic_core import parameters_sulfur as ps
from src.kinetic_core import uncertainty as U
from src.kinetic_core.engine import SULFUR, core_parameters

T_REF_K = 418.15


def test_chan_values_are_the_ones_the_registry_prints():
    source = ps.FORMATION_EA_PRIOR_SOURCE["thiol_assembly"]
    for _what, value in E.CHAN_1994_TABLE_I_KJ_MOL:
        assert f"{value:.1f}" in source


def test_the_class_prior_is_measured_barriers_only():
    rows = E.measured_activation_energies()
    vals = E.class_values()
    assert len(vals[E.CARBOHYDRATE]) == 20 and len(vals[E.SULFUR_CLASS]) == 15
    assert all(r.ea_kj_mol > 0.0 for r in rows)
    # the sulfur network's own measured barriers are members, read from the module, not re-typed
    assert ps.ZHANG_EA_THIOL_TO_DISULFIDE_KJ_MOL in vals[E.SULFUR_CLASS]
    assert ps.ZHANG_EA_CYS_AMADORI_TO_ALPHA_DC_KJ_MOL in vals[E.CARBOHYDRATE]


def test_targets_are_the_31_constants_with_no_measured_barrier():
    targets = E.single_temperature_targets()
    keys = [k for k, _c, _e in targets]
    assert len(keys) == 31 and len(set(keys)) == 31
    for key, klass, shipped in targets:
        assert key not in ps.MEASURED_EA_OVERRIDES
        assert ps.decay_family_of(key) is None
        assert klass in (E.CARBOHYDRATE, E.SULFUR_CLASS)
        assert (shipped is None) == (key in ps.NO_EA_KEYS)
    assert dict((k, c) for k, c, _ in targets)["k_mmft"] == E.SULFUR_CLASS


def test_empirical_quantile_never_leaves_the_measured_range():
    v = E.class_values()[E.CARBOHYDRATE]
    assert E.empirical_quantile(v, 0.0) == min(v) and E.empirical_quantile(v, 1.0) == max(v)
    assert min(v) <= E.empirical_quantile(v, 0.37) <= max(v)


def test_the_engine_override_holds_k_at_145C_and_moves_the_slope():
    base = core_parameters(SULFUR)
    drawn = core_parameters(SULFUR, frozen={"single_temperature_ea_kj_mol": {"k_pent_dpo": 150.0, "k_mmft": 90.0}})
    for key in ("k_pent_dpo", "k_mmft"):
        assert drawn[key].k_at(T_REF_K) == pytest.approx(base[key].k_at(T_REF_K), rel=1e-12)
        assert drawn[key].k_at(373.15) < base[key].k_at(373.15)
    assert drawn["k_nf_mft"] == base["k_nf_mft"]


def test_the_engine_refuses_to_overwrite_a_measured_barrier():
    with pytest.raises(ValueError, match="measured barrier"):
        core_parameters(SULFUR, frozen={"single_temperature_ea_kj_mol": {"k_arp_dpo": 50.0}})


def _draw(priors, index=0):
    child = np.random.SeedSequence(0).spawn(index + 1)[index]
    return U.draw_from_rng(np.random.default_rng(child), index, priors)


def test_b55_rows_reach_the_draw_and_move_nothing_else():
    off = U.core_priors(single_temperature_ea_prior=False)
    on = U.core_priors(single_temperature_ea_prior=True)
    added = [p for p in on if p.key.startswith("b55.")]
    assert len(added) == 31 and len(on) == len(off) + 31
    a, b = _draw(off), _draw(on)
    block = b.core.maillard["single_temperature_ea_kj_mol"]
    assert set(block) == {k for k, _c, _e in E.single_temperature_targets()}
    lo = min(min(v) for v in E.class_values().values())
    hi = max(max(v) for v in E.class_values().values())
    assert all(lo <= x <= hi for x in block.values())
    assert "single_temperature_ea_kj_mol" not in (a.core.maillard or {})
    shared = [k for k in a.coordinates if k in b.coordinates]
    assert all(a.coordinates[k] == b.coordinates[k] for k in shared)


def test_b56_normals_are_keyed_and_reproducible():
    x = U._discrepancy_normals(0, "bench", "MFT", 5)
    assert np.array_equal(x, U._discrepancy_normals(0, "bench", "MFT", 5))
    assert not np.array_equal(x, U._discrepancy_normals(0, "bench", "FFT", 5))
    # a longer run extends, never reshuffles, the per-draw sequence
    assert np.array_equal(x, U._discrepancy_normals(0, "bench", "MFT", 9)[:5])


def test_b56_sigmas_come_from_the_in_sample_artifact():
    assert U._discrepancy_sigmas(False) == {}
    payload = json.loads((data_paths.VALIDATION_DIR / "kinetic_core_b56_structural_discrepancy.json").read_text())
    assert U._discrepancy_sigmas(True) == payload["sigma_log10_by_lane"]
    assert set(payload["sigma_log10_by_lane"]) == {"sulfur", "trunk", "acrylamide", "lipid"}
    for lane in payload["lanes"].values():
        for fit in lane["fits"]:
            assert fit["fit"].startswith("results/validation/kinetic_core_b"), "in-sample fit reports only"


def test_shipped_switches():
    """B55 SHIPPED, B56 NOT SHIPPED (kinetic_core_b55_b56_ship_rule.md, 2026-10-28)."""
    assert U.B55_SINGLE_TEMPERATURE_EA_PRIOR is True
    assert U.B56_STRUCTURAL_DISCREPANCY is False
    assert sum(1 for p in U.CORE_PRIORS if p.key.startswith("b55.")) == 31
    assert U._discrepancy_sigmas(None) == {}
