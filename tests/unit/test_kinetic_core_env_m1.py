"""ENV-M1 (2026-09-11): one random stream per coordinate (kinetic_core_env_m1_prereg.md)."""
from __future__ import annotations

import numpy as np

from src.kinetic_core import uncertainty as U


def _draw(priors, seed=0, index=0):
    child = np.random.SeedSequence(seed).spawn(index + 1)[index]
    return U.draw_from_rng(np.random.default_rng(child), index, priors)


def test_m2_two_draws_at_the_same_seed_are_identical():
    a, b = _draw(U.CORE_PRIORS), _draw(U.CORE_PRIORS)
    assert a.coordinates == b.coordinates
    assert (a.k_aw_multiplier, a.hs_spme_multiplier) == (b.k_aw_multiplier, b.hs_spme_multiplier)


def test_m1_removing_a_prior_block_moves_no_other_coordinate():
    """The property the whole ENV method note wanted and could not have under one shared stream."""
    full = U.CORE_PRIORS
    for prefix in ("b34.", "b13.", "b18."):
        without = tuple(p for p in full if not p.key.startswith(prefix))
        assert len(without) < len(full)
        a, b = _draw(full), _draw(without)
        shared = [k for k in a.coordinates if k in b.coordinates]
        assert shared and all(a.coordinates[k] == b.coordinates[k] for k in shared), prefix
        assert not any(k.startswith(prefix) for k in b.coordinates)
        # the sulfur lane's joint block and the observable multipliers are untouched as well
        assert (a.k_aw_multiplier, a.hs_spme_multiplier) == (b.k_aw_multiplier, b.hs_spme_multiplier)
        assert a.core.maillard.get("log10_k_ref_at_145C") == b.core.maillard.get("log10_k_ref_at_145C")


def test_a_reordered_prior_table_gives_the_same_draw():
    full = list(U.CORE_PRIORS)
    a = _draw(tuple(full))
    b = _draw(tuple(reversed(full)))
    assert a.coordinates == b.coordinates


def test_the_stream_is_keyed_by_name_not_position():
    s = U._KeyedStreams(np.random.default_rng(3))
    x = s.for_key("b34.k_tdg_ddg.ea_kj_mol").uniform(0, 1)
    y = s.for_key("b34.k_tdg_ddg.ea_kj_mol").uniform(0, 1)
    z = s.for_key("b34.k_glc_tdg.ea_kj_mol").uniform(0, 1)
    assert x == y and x != z
