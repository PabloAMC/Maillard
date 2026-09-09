"""
Wave B10 (2026-09-06): one formation barrier -> two by route; the ambient oxidant
charged consistently; Yiltirak 2026's six within-study folds as fit rows.
Pre-registered in results/validation/kinetic_core_b10_prereg.md.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

from src import data_paths
from src.kinetic_core import engine, parameters_sulfur as ps
from src.kinetic_core.engine import SULFUR, core_parameters, frozen_parameters
from tests.support import wave_generator


# ---------------------------------------------------------------------------
# 1. The route table is total, and one number still means one number
# ---------------------------------------------------------------------------


def test_every_key_that_can_reach_the_formation_barrier_has_exactly_one_route():
    reachable = set(ps.FITTED_SULFUR_KEYS) - set(ps.NO_EA_KEYS) - set(ps.MEASURED_EA_OVERRIDES)
    assert set(ps.FORMATION_ROUTE_OF) == reachable
    assert set(ps.FORMATION_ROUTE_OF.values()) == set(ps.FORMATION_ROUTES)
    by_route = {r: [k for k, v in ps.FORMATION_ROUTE_OF.items() if v == r] for r in ps.FORMATION_ROUTES}
    assert min(len(v) for v in by_route.values()) >= 10
    # the thiol-sink family falls back to the thiol route, the carbonyl-sink family to the sugar route
    for key in ps.DECAY_FAMILY_THIOL_SINK:
        assert ps.FORMATION_ROUTE_OF[key] == "thiol_assembly", key
    for key in ("k_nf_decay", "k_fur_decay", "k_osone_decay"):
        assert ps.FORMATION_ROUTE_OF[key] == "sugar_trunk", key


def test_a_float_and_an_equal_mapping_build_identical_parameters():
    fitted = {k: -2.0 for k in ps.FITTED_SULFUR_KEYS}
    decay = {"thiol_sink": 60.2, "carbonyl_sink": 249.9}
    a = ps.with_fitted_sulfur(fitted, 64.0, decay)
    b = ps.with_fitted_sulfur(fitted, {"sugar_trunk": 64.0, "thiol_assembly": 64.0}, decay)
    for key in a:
        assert a[key].ea_kj_mol == b[key].ea_kj_mol, key
        assert a[key].k_ref == b[key].k_ref, key


def test_a_mapping_routes_each_step_to_its_own_barrier():
    fitted = {k: -2.0 for k in ps.FITTED_SULFUR_KEYS}
    built = ps.with_fitted_sulfur(fitted, {"sugar_trunk": 80.0, "thiol_assembly": 110.0},
                                  {"thiol_sink": 60.2, "carbonyl_sink": 249.9})
    assert built["k_pent_dpo"].ea_kj_mol == pytest.approx(80.0)
    assert built["k_nf_mft"].ea_kj_mol == pytest.approx(110.0)
    assert built["k_fur_fft"].ea_kj_mol == pytest.approx(110.0)
    # family barriers and measured overrides are untouched by the routes
    assert built["k_mft_decay"].ea_kj_mol == pytest.approx(60.2)
    assert built["k_fur_decay"].ea_kj_mol == pytest.approx(249.9)
    for key, value in ps.MEASURED_EA_OVERRIDES.items():
        assert built[key].ea_kj_mol == pytest.approx(value), key
    # without a family barrier the family keys take their route's
    built2 = ps.with_fitted_sulfur(fitted, {"sugar_trunk": 80.0, "thiol_assembly": 110.0})
    assert built2["k_mft_decay"].ea_kj_mol == pytest.approx(110.0)
    assert built2["k_fur_decay"].ea_kj_mol == pytest.approx(80.0)


def test_the_bands_are_the_prefactor_narrowed_bands_around_sourced_centres():
    for route in ps.FORMATION_ROUTES:
        lo, hi = ps.FORMATION_EA_BOUNDS_BY_ROUTE[route]
        centre = ps.FORMATION_EA_PRIOR_CENTRE[route]
        assert lo < centre < hi
        assert hi - lo <= 2 * 48.0 + 1e-9   # 12 decades at 8.0 kJ/mol per decade, both sides
        assert ps.LUMPED_FORMATION_EA_BOUNDS[0] <= lo and hi <= ps.LUMPED_FORMATION_EA_BOUNDS[1]
        assert len(ps.FORMATION_EA_PRIOR_SOURCE[route]) > 60
    assert ps.FORMATION_EA_PRIOR_CENTRE["sugar_trunk"] == pytest.approx(ps.ZHANG_EA_CYS_AMADORI_TO_ALPHA_DC_KJ_MOL)


# ---------------------------------------------------------------------------
# 2. The engine: a pre-B10 report is reproduced; a route block is honoured
# ---------------------------------------------------------------------------


def test_the_shipped_report_without_a_route_block_gives_every_route_the_lumped_barrier():
    frozen = frozen_parameters(SULFUR)
    built = core_parameters(SULFUR)
    lumped = frozen["lumped_formation_Ea_kJ_mol"]
    if "formation_Ea_by_route_kJ_mol" in frozen:
        pytest.skip("the engine reads a B10-or-later report; the route test below covers it")
    for key in ps.FORMATION_ROUTE_OF:
        if ps.decay_family_of(key) is None:
            assert built[key].ea_kj_mol == pytest.approx(lumped), key


def test_a_draw_carrying_a_route_block_moves_only_that_route():
    frozen = frozen_parameters(SULFUR)
    base = core_parameters(SULFUR)
    lumped = frozen.get("formation_Ea_by_route_kJ_mol", {}).get("sugar_trunk", frozen["lumped_formation_Ea_kJ_mol"])
    thiol = frozen.get("formation_Ea_by_route_kJ_mol", {}).get("thiol_assembly", frozen["lumped_formation_Ea_kJ_mol"])
    frozen["formation_Ea_by_route_kJ_mol"] = {"sugar_trunk": lumped, "thiol_assembly": thiol + 20.0}
    moved = core_parameters(SULFUR, frozen=frozen)
    assert moved["k_nf_mft"].ea_kj_mol == pytest.approx(base["k_nf_mft"].ea_kj_mol + 20.0)
    assert moved["k_pent_dpo"].ea_kj_mol == pytest.approx(base["k_pent_dpo"].ea_kj_mol)
    for key in ps.MEASURED_EA_OVERRIDES:
        assert moved[key].ea_kj_mol == pytest.approx(base[key].ea_kj_mol)


# ---------------------------------------------------------------------------
# 3. The ambient oxidant: fit and deployment charge the same number
# ---------------------------------------------------------------------------


def test_engine_and_fit_generator_charge_the_same_ambient_oxidant():
    with wave_generator("generate_kinetic_core_b2_3_fit") as B23:
        assert ps.OX_AMBIENT_MMOL_L == B23.OX_AMBIENT_MMOL_L == 1.0


def test_the_ambient_charge_changes_a_trace_thiol_prediction_by_under_one_percent():
    from src.kinetic_core import panel

    bench = panel.load_bundle(data_paths.BENCHMARKS_DIR / "hofmann1998_ribose_cysteine_145C_20min_pH5.json")
    spec = panel.core_spec(bench, use_buffer=True)
    decl = engine.declare_envelope(spec, list(panel.bundle_targets(bench)))
    operative = core_parameters(SULFUR)
    with_ambient, _ = engine._integrate_program(SULFUR, operative, dict(decl.mapped_precursors), spec.process)
    assert with_ambient["OX"] == pytest.approx(ps.OX_AMBIENT_MMOL_L, rel=0.05)
    none, _ = engine._integrate_program(SULFUR, operative, {**decl.mapped_precursors, "OX": 0.0}, spec.process)
    for key in ("MFT", "FFT"):
        assert abs(with_ambient[key] / none[key] - 1.0) < 0.01, key
    assert none["MFTD"] == 0.0 and with_ambient["MFTD"] > 0.0


# ---------------------------------------------------------------------------
# 4. The generator: rows, vector, bounds, incumbent
# ---------------------------------------------------------------------------


@pytest.fixture
def b10():
    with wave_generator("generate_kinetic_core_b10_fit") as mod:
        mod.configure(True)
        yield mod
        mod.restore()


def _b23():
    import generate_kinetic_core_b2_3_fit as B23
    return B23


def test_b10_installs_sixty_rows_and_twenty_five_free_coordinates(b10):
    B23 = _b23()
    assert len(B23.ACTIVE_FIT_ROWS) == 60
    assert len(b10.FREE_KEYS) == 25 and len(b10.ALL_KEYS) == 49
    assert b10.ALL_KEYS[b10.SUGAR_SLOT] == "Ea_lumped_formation"
    assert b10.ALL_KEYS[b10.THIOL_SLOT] == "Ea_thiol_assembly"
    assert set(b10.FREE_KEYS) >= {"Ea_lumped_formation", "Ea_thiol_assembly"}
    lower, upper = b10.full_bounds()
    assert len(lower) == len(upper) == 49
    assert (lower[b10.SUGAR_SLOT], upper[b10.SUGAR_SLOT]) == ps.FORMATION_EA_BOUNDS_BY_ROUTE["sugar_trunk"]
    assert (lower[b10.THIOL_SLOT], upper[b10.THIOL_SLOT]) == ps.FORMATION_EA_BOUNDS_BY_ROUTE["thiol_assembly"]


def test_the_six_folds_are_table_s3_exactly_and_declare_their_bundles(b10):
    rows = {r["id"]: r for r in b10.B10_FIT_ROWS}
    assert rows["yiltirak_MFT_fold_110C_over_100C"]["target"] == pytest.approx(3.29 / 6.88)
    assert rows["yiltirak_MFT_fold_130C_over_120C"]["target"] == pytest.approx(1.71 / 2.4)
    assert rows["yiltirak_FFT_fold_120C_over_110C"]["target"] == pytest.approx(1.68 / 1.46)
    for r in rows.values():
        assert r["kind"] == "cross_system_ratio" and r["sigma_log"] == 0.10
        assert r["benchmark_id"].endswith("_Yiltirak2026") and r["benchmark_id_b"].endswith("_Yiltirak2026")
        assert "Table S3" in r["anchor"]
    B23 = _b23()
    for name, spec in b10.B10_SYSTEMS.items():
        assert B23.SYSTEMS[name] is spec
        assert spec["initial"] == {"PENT": 25.0, "Cys": 25.0, "OX": ps.OX_AMBIENT_MMOL_L}
        assert spec["ph"] == 5.5 and spec["buffer"].phosphate_mol_l == 0.5
    assert {s["t_c"] * s["minutes"] for s in b10.B10_SYSTEMS.values()} == {24000.0, 13200.0, 7200.0, 3900.0}


def test_the_incumbent_is_b9_with_the_route_barriers_at_their_centres(b10):
    x = b10.incumbent_vector()
    report = json.loads((data_paths.VALIDATION_DIR / "kinetic_core_b9_fit_report.json").read_text())
    fr = report["frozen_parameters"]
    B23 = _b23()
    for i, key in enumerate(B23.PARAM_ORDER):
        assert x[i] == pytest.approx(fr["log10_k_ref_at_145C"][key])
    assert x[b10.SUGAR_SLOT] == pytest.approx(ps.FORMATION_EA_PRIOR_CENTRE["sugar_trunk"])
    assert x[b10.THIOL_SLOT] == pytest.approx(ps.FORMATION_EA_PRIOR_CENTRE["thiol_assembly"])
    assert x[B23.N_K + 1] == pytest.approx(fr["decay_Ea_kJ_mol"]["thiol_sink"])
    # the builder reads both slots
    params = b10.build_parameters(x)
    assert params["k_pent_dpo"].ea_kj_mol == pytest.approx(x[b10.SUGAR_SLOT])
    assert params["k_nf_mft"].ea_kj_mol == pytest.approx(x[b10.THIOL_SLOT])


def test_the_leave_yiltirak_out_variant_is_b9s_objective(b10):
    B23 = _b23()
    b10.configure(False)
    try:
        assert len(B23.ACTIVE_FIT_ROWS) == 54
        assert b10.OUT_FIT_REPORT.name == "kinetic_core_b10_noyil_fit_report.json"
    finally:
        b10.configure(True)
        assert len(B23.ACTIVE_FIT_ROWS) == 60


def test_the_laplace_vector_round_trips_a_route_report(b10):
    from generate_kinetic_core_b8_laplace import coordinate_of, frozen_vector

    x = b10.incumbent_vector()
    fake = {"frozen_parameters": {
        "log10_k_ref_at_145C": {k: float(x[i]) for i, k in enumerate(_b23().PARAM_ORDER)},
        "lumped_formation_Ea_kJ_mol": float(x[b10.SUGAR_SLOT]),
        "decay_Ea_kJ_mol": {"thiol_sink": float(x[b10.SUGAR_SLOT + 1]), "carbonyl_sink": float(x[b10.SUGAR_SLOT + 2])},
        "ph_drift": {"acid_yield_per_sink_event": float(x[b10.SUGAR_SLOT + 3]), "arp_secondary_ammonium_pKa": float(x[b10.SUGAR_SLOT + 4])},
        "formation_Ea_by_route_kJ_mol": b10.route_ea_from_vector(x),
    }}
    np.testing.assert_allclose(frozen_vector(fake), x)
    assert coordinate_of("Ea_thiol_assembly") == {"block": "formation_Ea_by_route_kJ_mol", "key": "thiol_assembly", "kind": "Ea"}
