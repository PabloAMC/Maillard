"""
tests/unit/test_kinetic_core_b51.py

WAVE B51 -- the pH state corrected against its thermodynamic sources (2026-10-27).

The six ``ph_state.B51`` switches must (i) reproduce every earlier wave bit for bit when off, and
(ii) carry exactly the published constants when on. Every test restores the switches it flips.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from src.kinetic_core import parameters_sulfur as ps
from src.kinetic_core import ph_state as PH
from src.kinetic_core import sulfur as S

T25, T100, T145 = 298.15, 373.15, 418.15
_DRIFT = PH.PhDrift(acid_yield=0.9, arp_amine_pka=10.0)
_PHOSPHATE = PH.BufferSpec(kind="phosphate", phosphate_mol_l=0.5, declared=True, source="test")


@pytest.fixture
def switches():
    """Yield a setter for ph_state.B51 and restore the shipped switches afterwards."""
    saved = PH.B51

    def _set(**kw):
        PH.B51 = PH.PhStateCorrections(**kw)

    yield _set
    PH.B51 = saved


def _toy_parameters():
    from tests.unit.test_kinetic_core_b2_2 import _toy_parameters as toy

    return toy()


def test_all_switches_off_is_the_b50_pka_law(switches):
    switches()
    assert ps.pka_at("h2s_1", T145) == pytest.approx(5.944, abs=1e-3)
    assert ps.pka_at("thiol", T145) == pytest.approx(7.023, abs=1e-3)
    assert PH.PHOSPHATE.pkas_at(T145)[1] == pytest.approx(7.017, abs=1e-3)
    assert PH.davies_term(0.5, T145) == 0.0


def test_heat_capacity_raises_every_pk_with_a_negative_dcp(switches):
    switches(heat_capacity=True)
    # phosphate pK2: dH +3.6, dCp -230 J/(K mol) (GKL Table 7.51) -> 7.63 at 145 C
    assert PH.PHOSPHATE.pkas_at(T145)[1] == pytest.approx(7.633, abs=2e-3)
    # the law is exact at the reference temperature
    assert PH.PHOSPHATE.pkas_at(T25) == pytest.approx(PH.PHOSPHATE.pkas, abs=1e-12)


def test_the_heat_capacity_law_reproduces_the_tabulated_water_ion_product(switches):
    """Water is the one acid whose high-T pK is TABULATED; the dCp law must get it right too."""
    switches(heat_capacity=True)
    # water: dH 55.81 kJ/mol, dCp -224 J/(K mol) at 25 C. A CONSTANT dCp is itself an approximation
    # (water's dCp grows more negative with T): 0.04 units off at 100 C, 0.11 at 150 C -- against the
    # 0.27 / 0.75 that constant-dH van 't Hoff misses by.
    for t_c, tol in ((100.0, 0.06), (150.0, 0.15)):
        t = t_c + 273.15
        assert PH.pk_at(13.995, 55.81, -224.0, t) == pytest.approx(PH.pkw_at(t), abs=tol)


def test_cysteine_carries_gkl_table_7_23(switches):
    switches(cysteine_gkl=True)
    assert PH.cysteine_group().pkas == (1.71, 8.36, 10.75)
    assert PH.cysteine_group().delta_h_kj_mol == (-0.6, 36.1, 34.1)
    assert ps.PKA_GKL_CYSTEINE["thiol"][:2] == PH.CYSTEINE_GKL.pkas[1:2] + PH.CYSTEINE_GKL.delta_h_kj_mol[1:2]
    switches()
    assert PH.cysteine_group() is PH.CYSTEINE_B9


def test_h2s_pk1_follows_suleimenov_and_seward(switches):
    assert ps.h2s_pk1_suleimenov_seward(T25) == pytest.approx(6.97, abs=0.02)
    assert ps.h2s_pk1_suleimenov_seward(T145) == pytest.approx(6.49, abs=0.01)
    # the measured correlation has a shallow minimum near 100-130 C; constant-dH has none
    grid = [ps.h2s_pk1_suleimenov_seward(273.15 + c) for c in range(25, 201, 5)]
    assert min(grid) > 6.4 and grid.index(min(grid)) not in (0, len(grid) - 1)
    switches(h2s_suleimenov=True)
    assert ps.pka_at("h2s_1", T145) == pytest.approx(ps.h2s_pk1_suleimenov_seward(T145))


def test_davies_a_is_the_textbook_value_and_rises_with_temperature(switches):
    assert PH.davies_a(T25) == pytest.approx(0.511, abs=0.003)
    assert PH.davies_a(T145) > PH.davies_a(T100) > PH.davies_a(T25)
    switches(ionic_strength=True)
    s = math.sqrt(0.5)
    assert PH.davies_term(0.5, T25) == pytest.approx(PH.davies_a(T25) * (s / (1 + s) - 0.15))


def test_ionic_strength_lowers_phosphate_pk2_by_three_a_f(switches):
    switches(ionic_strength=True)
    af = PH.davies_term(0.5, T25)
    assert PH.PHOSPHATE.pkas_at(T25, 0.5)[1] == pytest.approx(7.198 - 3.0 * af, abs=1e-9)


def test_ionic_strength_of_a_half_molar_phosphate_buffer(switches):
    switches(ionic_strength=True)
    inventory = PH.titratable_inventory({}, _DRIFT, _PHOSPHATE)
    sid = PH.strong_ion_difference(7.0, {}, _DRIFT, _PHOSPHATE)
    ph, ionic = PH.solve_ph_and_ionic_strength(T25, sid, inventory)
    assert ph == pytest.approx(7.0, abs=1e-6)  # the label is recovered at the label temperature
    assert 0.8 < ionic < 1.6                    # Na2HPO4 / NaH2PO4 mix near pH 7


def test_pyrophosphate_and_phthalate_are_chargeable_buffers():
    for kind, attr in (("pyrophosphate", "pyrophosphate_mol_l"), ("phthalate", "phthalate_mol_l")):
        spec = PH.BufferSpec(kind=kind, declared=True, source="test", **{attr: 0.2})
        (group, mmol), = spec.groups()
        assert mmol == pytest.approx(200.0)
        assert group.delta_cp_j_mol_k and len(group.delta_cp_j_mol_k) == len(group.pkas)
        with pytest.raises(ValueError):
            PH.BufferSpec(kind=kind, declared=True, source="test")


def test_the_measured_tables_read_the_cooled_ph_only_when_switched_on(switches):
    p = _toy_parameters()
    switches()
    hot = S.sulfur_rate_constants_at(p, T145, 5.0, ph_label=6.0)
    assert hot["r_cys_h2s"] == pytest.approx(ps.cysteine_thermolysis_k(5.0, T145))
    switches(label_scale_tables=True)
    cool = S.sulfur_rate_constants_at(p, T145, 5.0, ph_label=6.0)
    assert cool["r_cys_h2s"] == pytest.approx(ps.cysteine_thermolysis_k(6.0, T145))


def test_the_sid_is_carried_across_segments_only_when_switched_on(switches):
    p = _toy_parameters()
    grid = np.array([0.0, 30.0])

    def two_segments():
        first = S.integrate_sulfur(p, T145, {"ARP": 20.0, "Cys": 20.0, "OX": 1.0}, grid, ph=6.0,
                                   buffer_spec=PH.DEFAULT_BUFFER, ph_drift=_DRIFT, rtol=1e-8, atol=1e-14)
        state = {k: float(first.concentrations[-1, S.SULFUR_INDEX[k]]) for k in S.SULFUR_STATE_KEYS}
        second = S.integrate_sulfur(p, T145, state, grid, ph=6.0, buffer_spec=PH.DEFAULT_BUFFER,
                                    ph_drift=_DRIFT, rtol=1e-8, atol=1e-14,
                                    sid_mmol_l=first.metadata["strong_ion_difference_mmol_L"])
        return first, second

    switches()
    first, second = two_segments()
    # off: the second segment re-titrates its acidified start back to the label
    assert second.metadata["ph_initial_in_situ"] > first.metadata["ph_final_in_situ"] + 1e-3
    switches(sid_across_segments=True)
    first, second = two_segments()
    assert second.metadata["strong_ion_difference_mmol_L"] == first.metadata["strong_ion_difference_mmol_L"]
    assert second.metadata["ph_initial_in_situ"] == pytest.approx(first.metadata["ph_final_in_situ"], abs=2e-3)
