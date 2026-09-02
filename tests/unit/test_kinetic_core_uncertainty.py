"""
tests/unit/test_kinetic_core_uncertainty.py

Retirement step B2 -- the Monte-Carlo envelope computed ON THE KINETIC CORE.
Contract tests: the priors table's shape and honesty labels, determinism of
the engine's draw API and of the artifact, and the artifact's field
contract. No rate constant is pinned here.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from src import data_paths
from src.kinetic_core import panel
from src.kinetic_core import uncertainty as unc
from src.kinetic_core.engine import (
    ACRYLAMIDE,
    MAILLARD_LANES,
    SULFUR,
    TRUNK,
    CoreDraw,
    FormulationSpec,
    ProcessSpec,
    ThermalProgram,
    core_parameters,
    frozen_parameters,
    predict,
)

ROOT = Path(__file__).resolve().parents[2]

MFT = "2-Methyl-3-furanthiol (MFT)"
FFT = "2-Furfurylthiol (FFT)"


def _sulfur_spec():
    return FormulationSpec(
        "sulfur",
        {"D-Ribose": 100.0, "L-Cysteine": 33.0},
        ProcessSpec(ThermalProgram.isothermal(145.0, 20.0), ph=5.0),
    )


def _acrylamide_spec():
    return FormulationSpec(
        "acrylamide",
        {"D-Glucose": 27.75, "L-Asparagine": 33.3},
        ProcessSpec(ThermalProgram.isothermal(180.0, 30.0), ph=6.0, water_activity=0.99),
    )


def _lipid_spec():
    return FormulationSpec(
        "lipid",
        {"Pea Protein Isolate": 1000.0},
        ProcessSpec(ThermalProgram.isothermal(140.0, 0.1), ph=7.1, matrix="pea_iso"),
    )


def _glucose_spec():
    return FormulationSpec(
        "glucose",
        {"D-Glucose": 555.0},
        ProcessSpec(ThermalProgram.isothermal(121.0, 18.0), ph=4.36),
    )


# ---------------------------------------------------------------------------
# 1. The engine's draw API is a no-op at the centre
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "spec, targets",
    [
        (_sulfur_spec(), [MFT, FFT]),
        (_acrylamide_spec(), ["Acrylamide", "5-Hydroxymethylfurfural (HMF)"]),
        (_glucose_spec(), ["5-HMF", "DMHF"]),
    ],
)
def test_predict_defaults_equal_the_explicit_centre(spec, targets):
    base = predict(spec, targets)
    same = predict(spec, targets, draw=None, size_declared_bands=True)
    assert base.concentrations_ug_per_l == same.concentrations_ug_per_l
    assert base.species_mmol_per_l == same.species_mmol_per_l
    assert base.run_metadata == same.run_metadata
    # a draw whose fields are all None is the centre too
    centre = predict(spec, targets, draw=CoreDraw())
    assert centre.concentrations_ug_per_l == base.concentrations_ug_per_l
    assert CoreDraw().is_centre


@pytest.mark.parametrize("lane", MAILLARD_LANES)
def test_core_parameters_from_the_frozen_vector_are_the_report_parameters(lane):
    direct = core_parameters(lane)
    via = core_parameters(lane, frozen=frozen_parameters(lane))
    assert set(direct) == set(via)
    for key in direct:
        assert direct[key] == via[key], key


def test_frozen_vector_moves_the_fit_coordinates_not_operative_constants():
    frozen = frozen_parameters(ACRYLAMIDE)
    frozen["fitted_Ea_kJ_mol"]["Ea_acr_dp"] += 10.0
    moved = core_parameters(ACRYLAMIDE, frozen=frozen)
    base = core_parameters(ACRYLAMIDE)
    assert moved["k_acr_dp"].ea_kj_mol == pytest.approx(base["k_acr_dp"].ea_kj_mol + 10.0)
    # the shared trunk pairs are untouched
    assert moved["k_mgo_mel"] == base["k_mgo_mel"]


def test_sulfur_frozen_vector_honours_measured_overrides_and_no_ea_keys():
    from src.kinetic_core.parameters_sulfur import MEASURED_EA_OVERRIDES, NO_EA_KEYS

    frozen = frozen_parameters(SULFUR)
    frozen["lumped_formation_Ea_kJ_mol"] += 15.0
    moved = core_parameters(SULFUR, frozen=frozen)
    base = core_parameters(SULFUR)
    for key in MEASURED_EA_OVERRIDES:
        assert moved[key].ea_kj_mol == pytest.approx(MEASURED_EA_OVERRIDES[key])
    for key in NO_EA_KEYS:
        assert moved[key].ea_kj_mol is None
    # a step that shares the lumped Ea moved by exactly the offset
    shared = [
        k for k, p in base.items()
        if getattr(p, "flags", None) and "shares_lumped_formation_Ea" in p.flags
    ]
    assert shared
    assert moved[shared[0]].ea_kj_mol == pytest.approx(base[shared[0]].ea_kj_mol + 15.0)


def test_size_declared_bands_false_equals_the_point_at_the_centre_draw():
    for spec, targets in (
        (_glucose_spec(), ["DMHF", "5-HMF"]),
        (_lipid_spec(), ["hexanal"]),
    ):
        full = predict(spec, targets)
        bare = predict(spec, targets, draw=CoreDraw(), size_declared_bands=False)
        assert bare.concentrations_ug_per_l == full.concentrations_ug_per_l
        assert bare.species_mmol_per_l == full.species_mmol_per_l
        # the bands were priced by re-integration in the full run and not here
        assert (
            full.run_metadata.get("furanic_extra_decades")
            or full.run_metadata.get("lipid_extra_decades")
        )
        assert not bare.run_metadata.get("furanic_extra_decades")
        assert not bare.run_metadata.get("lipid_extra_decades")
        assert bare.run_metadata["declared_bands_sized"] is False


def test_a_lipid_draw_moves_the_lipid_lane_within_the_declared_bands():
    spec = _lipid_spec()
    base = predict(spec, ["hexanal"], draw=CoreDraw(), size_declared_bands=False)
    high = predict(
        spec, ["hexanal"],
        draw=CoreDraw(q10=3.0, lipid_fraction_scale=2.0, peroxide_scale=3.0),
        size_declared_bands=False,
    )
    assert high.concentrations_ug_per_l["hexanal"] > base.concentrations_ug_per_l["hexanal"]
    # a scale far outside the band is clipped to the carrier's declared band
    clipped = predict(
        spec, ["hexanal"], draw=CoreDraw(lipid_fraction_scale=1e6, peroxide_scale=1e6),
        size_declared_bands=False,
    )
    corner = predict(
        spec, ["hexanal"], draw=CoreDraw(lipid_fraction_scale=2.4, peroxide_scale=4.0),
        size_declared_bands=False,
    )
    assert clipped.concentrations_ug_per_l["hexanal"] == pytest.approx(
        corner.concentrations_ug_per_l["hexanal"]
    )


def test_operative_parameters_and_draw_maillard_cannot_both_be_passed():
    with pytest.raises(ValueError):
        predict(
            _sulfur_spec(), [MFT],
            parameters=core_parameters(SULFUR),
            draw=CoreDraw(maillard=frozen_parameters(SULFUR)),
        )


# ---------------------------------------------------------------------------
# 2. The priors table
# ---------------------------------------------------------------------------

REQUIRED_PRIOR_FIELDS = {
    "key", "lane", "kind", "distribution", "centre", "sigma", "band", "unit",
    "source", "sampled", "reason",
}


def test_priors_table_shape():
    priors = unc.CORE_PRIORS
    assert priors == unc.core_priors()
    keys = [p.key for p in priors]
    assert len(keys) == len(set(keys)), "prior keys must be unique"
    for p in priors:
        assert REQUIRED_PRIOR_FIELDS <= set(p.as_dict())
        assert p.sampled is (p.distribution != "fixed")
        if p.sampled and p.distribution in ("normal", "normal_log10"):
            assert p.sigma is not None and p.sigma > 0.0
        if p.sampled and p.distribution in ("uniform", "log_uniform", "log_uniform_dispersion"):
            assert p.band is not None and p.band[0] < p.band[1]


def test_b1_identified_pairs_are_sampled_and_the_unidentified_are_fixed():
    by_key = {p.key: p for p in unc.CORE_PRIORS}
    for key in ("k_mgo_mel", "k_aa_frag"):
        assert by_key[f"b1.{key}.log10_k_ref_100C"].sampled
        assert by_key[f"b1.{key}.ea_kj_mol"].sampled
    for key in ("k_glc_frag", "k_fa_frag"):
        assert not by_key[f"b1.{key}.log10_k_ref_100C"].sampled
        assert not by_key[f"b1.{key}.ea_kj_mol"].sampled


def test_b3_only_the_identified_pair_is_sampled():
    b3 = [p for p in unc.CORE_PRIORS if p.key.startswith("b3.")]
    sampled = sorted(p.key for p in b3 if p.sampled)
    assert sampled == ["b3.Ea_acr_dp", "b3.k_acr_dp.log10_k_ref_160C"]
    report = json.loads(
        (data_paths.VALIDATION_DIR / "kinetic_core_b3_fit_report.json").read_text()
    )
    half = report["parameter_intervals"]["k_acr_dp"]["ci95_halfwidth"]
    by_key = {p.key: p for p in b3}
    assert by_key["b3.k_acr_dp.log10_k_ref_160C"].sigma == pytest.approx(half / 1.96)


def test_every_sulfur_prior_is_unsampled_with_the_declared_reason():
    sulfur = [p for p in unc.CORE_PRIORS if p.lane == SULFUR]
    assert len(sulfur) >= 40  # the 43 log10 k, the lumped Ea, two decay Ea, two drift
    for p in sulfur:
        assert p.sampled is False
        assert p.distribution == "fixed"
        assert p.reason == unc.NO_UNCERTAINTY
    assert SULFUR in {p.lane for p in unc.CORE_PRIORS if not p.sampled}


def test_declared_bands_are_uniform_over_the_declared_band():
    from src.kinetic_core.parameters_furanic import FURANONE_PARTITION_EA_BAND_KJ_MOL
    from src.kinetic_core.parameters_lipid import LIPID_CARRIERS, Q10_ASSUMPTION
    from src.kinetic_core.parameters_matrix import (
        HS_SPME_SAME_SAMPLE_DISPERSION,
        K_AW_UNCERTAINTY_DECADES,
    )

    by_key = {p.key: p for p in unc.CORE_PRIORS}
    q10 = by_key["lipid.q10"]
    assert q10.distribution == "uniform" and q10.band == (Q10_ASSUMPTION.lo, Q10_ASSUMPTION.hi)
    fur = by_key["furanic.partition_ea_offset_kj_mol"]
    assert fur.distribution == "uniform"
    assert fur.band == (-FURANONE_PARTITION_EA_BAND_KJ_MOL, FURANONE_PARTITION_EA_BAND_KJ_MOL)
    kaw = by_key["observable.air_water_partition_constant"]
    assert kaw.band == (-K_AW_UNCERTAINTY_DECADES, K_AW_UNCERTAINTY_DECADES)
    hs = by_key["observable.hs_spme_same_sample_dispersion"]
    assert hs.band == tuple(float(x) for x in HS_SPME_SAME_SAMPLE_DISPERSION)
    pea = LIPID_CARRIERS["pea_protein_isolate"]
    assert by_key["lipid.pea_protein_isolate.lipid_mass_fraction"].band == (pea.lipid_lo, pea.lipid_hi)
    assert by_key["lipid.pea_protein_isolate.peroxide_value_meq_per_kg"].band == (pea.pv_lo, pea.pv_hi)
    # the fed hydroperoxide has a degenerate band and is not sampled
    assert not by_key["lipid.frankel_pure_hydroperoxide.lipid_mass_fraction"].sampled


# ---------------------------------------------------------------------------
# 3. Draws
# ---------------------------------------------------------------------------


def test_sample_draws_is_deterministic_and_seed_sensitive():
    a = unc.sample_draws(6, seed=3)
    b = unc.sample_draws(6, seed=3)
    c = unc.sample_draws(6, seed=4)
    assert [d.coordinates for d in a] == [d.coordinates for d in b]
    assert [d.coordinates for d in a] != [d.coordinates for d in c]
    assert [d.index for d in a] == list(range(6))
    for d in a:
        assert 2.0 <= d.core.q10 <= 3.0
        assert -50.0 <= d.core.furanone_partition_ea_kj_mol <= 50.0
        assert 10 ** -0.5 <= d.k_aw_multiplier <= 10 ** 0.5
        assert 23 ** -0.5 <= d.hs_spme_multiplier <= 23 ** 0.5
        assert d.core.ph_drift is None  # B8: nothing to sample
        assert "log10_k_ref_at_145C" not in (d.core.maillard or {})  # sulfur never moved
        # the B1 Ea draws are clipped to the fit's own search bounds
        assert 20.0 <= d.core.maillard["variant_A_measured_sink"]["k_aa_frag"]["ea_kj_mol"] <= 260.0


def test_a_draw_runs_through_every_lane():
    draw = unc.sample_draws(1, seed=0)[0].core
    for spec, targets in (
        (_sulfur_spec(), [MFT]),
        (_acrylamide_spec(), ["Acrylamide"]),
        (_glucose_spec(), ["DMHF"]),
        (_lipid_spec(), ["hexanal"]),
    ):
        run = predict(spec, targets, draw=draw, size_declared_bands=False)
        assert run.answered
        assert all(v >= 0.0 for v in run.concentrations_ug_per_l.values())


# ---------------------------------------------------------------------------
# 4. The panel and the artifact
# ---------------------------------------------------------------------------


def test_panel_union_tags_and_exclusions():
    scored, skipped = panel.panel_bundles()
    tags = {tag for _, tag in scored}
    assert tags == {
        panel.PANEL_TRUST_LOOP, panel.PANEL_MAILLARD_PATH_HOLDOUT, panel.PANEL_EXTERNAL_MATRIX
    }
    names = {p.name for p, _ in scored}
    assert not any(n.endswith("_Internal2026.json") for n in names)
    assert {s["benchmark_id"] for s in skipped} >= {
        "pea_isolate_ribose_cysteine_100C_45min_Internal2026",
        "soy_isolate_ribose_cysteine_100C_45min_Internal2026",
    }
    for path, _ in scored:
        bench = panel.load_bundle(path)
        assert bench.get("evidence_class") != "diagnostic_only"
    # non-recursive over data/benchmarks: the quarantine never leaks in
    assert not any("quarantined" in str(p) or "step_level_unreachable" in str(p) for p, _ in scored)


def test_exam_generator_reads_the_bundle_through_panel():
    import scripts.generators.generate_cutover_final_exam as exam

    assert exam._core_spec is panel.core_spec
    assert exam._measured_value is panel.measured_value
    assert exam.buffer_from_bundle is panel.buffer_from_bundle
    assert exam.SHARED_WITH_HOLDOUT_PANEL is panel.SHARED_WITH_HOLDOUT_PANEL
    assert len(panel.SHARED_WITH_HOLDOUT_PANEL) == 4


_SMALL_PANEL = (
    (data_paths.BENCHMARKS_DIR / "hofmann1998_ribose_cysteine_145C_20min_pH5.json", "trust_loop"),
    (data_paths.BENCHMARKS_DIR / "hofmann1998_norfuraneol_h2s_145C_20min_pH5.json", "trust_loop"),
    (data_paths.BENCHMARKS_DIR / "pea_isolate_uht_140C_Trikusuma2019.json", "trust_loop"),
    (
        data_paths.MAILLARD_PATH_HOLDOUT_DIR
        / "mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7.json",
        "maillard_path_holdout",
    ),
)

BENCHMARK_FIELDS = {
    "benchmark_id", "bench_file", "panel", "execution_path", "protein_type", "fitted_row",
    "evidence_role", "signal_origin", "fit_target_of", "compounds", "refused_compounds",
    "quantification_family", "quantification_source",
}
COMPOUND_FIELDS = {
    "compound", "measured_ppb", "predicted_p5", "predicted_p50", "predicted_p95",
    "predicted_point", "inside_ci", "ci_width_log10", "lane", "shared_with", "target_unit",
    "observable_multipliers_applied", "in_core_fit",
}
SUMMARY_FIELDS = {
    "benchmark_count", "matched_compound_count", "ci_coverage_hits", "ci_coverage_rate",
    "ci_level_pct", "n_samples", "seed", "honest_literature_coverage", "signal_origin_split",
    "parameter_sources",
}


@pytest.fixture(scope="module")
def small_artifact():
    return unc.propagate_panel(_SMALL_PANEL, n_samples=8, seed=5, workers=1)


def test_artifact_contract_fields(small_artifact):
    payload = small_artifact
    assert SUMMARY_FIELDS <= set(payload["summary"])
    assert payload["summary"]["ci_level_pct"] == 90
    assert payload["summary"]["n_samples"] == 8 and payload["summary"]["seed"] == 5
    for source in payload["summary"]["parameter_sources"]:
        assert {"report", "sha256"} <= set(source) and len(source["sha256"]) == 64
    assert {"hits", "total", "rate", "median_ci_width_log10", "excluded_fitted_rows"} <= set(
        payload["summary"]["honest_literature_coverage"]
    )
    assert set(payload["summary"]["signal_origin_split"]) == {
        "external_literature", "internal_synthetic", "fitted_row", "in_core_fit"
    }
    oos = payload["summary"]["out_of_sample_literature_coverage"]
    lit = payload["summary"]["honest_literature_coverage"]
    assert oos["total"] + oos["rows_the_core_fit_read"]["total"] == lit["total"]
    assert payload["priors"] == [p.as_dict() for p in unc.CORE_PRIORS]
    assert payload["benchmarks"]
    for bench in payload["benchmarks"]:
        assert BENCHMARK_FIELDS <= set(bench)
        assert bench["panel"] in {"trust_loop", "maillard_path_holdout"}
        for row in bench["compounds"]:
            assert COMPOUND_FIELDS <= set(row)
            assert row["predicted_p5"] <= row["predicted_p50"] <= row["predicted_p95"]
            assert row["inside_ci"] == (row["predicted_p5"] <= row["measured_ppb"] <= row["predicted_p95"])
    for item in payload["refused_compounds"]:
        assert {"benchmark_id", "compound", "reason"} <= set(item) and item["reason"]


def test_shared_hofmann_rows_carry_shared_with(small_artifact):
    rows = {
        (b["benchmark_id"], c["compound"]): c
        for b in small_artifact["benchmarks"] for c in b["compounds"]
    }
    key = ("mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7", MFT)
    assert rows[key]["shared_with"] == "hofmann_ribose_pH7_MFT"
    assert rows[("hofmann1998_ribose_cysteine_145C_20min_pH5", MFT)]["shared_with"] is None


def test_predicted_point_is_the_deterministic_prediction(small_artifact):
    for bench in small_artifact["benchmarks"]:
        spec = panel.core_spec(panel.load_bundle(ROOT / bench["bench_file"]))
        for row in bench["compounds"]:
            run = predict(spec, [row["compound"]])
            assert row["predicted_point"] == pytest.approx(
                run.concentrations_ug_per_l[row["compound"]]
            )


def test_refused_rows_use_the_engine_refusal_vocabulary(small_artifact):
    refused = {(r["benchmark_id"], r["compound"]): r for r in small_artifact["refused_compounds"]}
    # H2S + norfuraneol: hydrogen sulfide is not a precursor the core can charge
    key = ("hofmann1998_norfuraneol_h2s_145C_20min_pH5", MFT)
    assert key in refused and "UNMAPPED PRECURSORS" in refused[key]["reason"]
    # 2-pentylfuran is on the named unrepresented list
    key = ("pea_isolate_uht_140C_Trikusuma2019", "2-pentylfuran")
    assert key in refused and "UNREPRESENTED TARGETS" in refused[key]["reason"]


def test_same_seed_gives_an_identical_artifact():
    a = unc.propagate_panel(_SMALL_PANEL[:1], n_samples=6, seed=11, workers=1)
    b = unc.propagate_panel(_SMALL_PANEL[:1], n_samples=6, seed=11, workers=1)
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)
    c = unc.propagate_panel(_SMALL_PANEL[:1], n_samples=6, seed=12, workers=1)
    assert json.dumps(a["benchmarks"], sort_keys=True) != json.dumps(c["benchmarks"], sort_keys=True)


def test_every_trust_loop_benchmark_is_scored_or_refused_with_a_reason():
    """Cheap: the deterministic pass only (n_samples small), over the whole trust loop."""
    scored, _ = panel.panel_bundles()
    trust = [(p, t) for p, t in scored if t == panel.PANEL_TRUST_LOOP]
    payload = unc.propagate_panel(trust, n_samples=5, seed=0, workers=1)
    seen = {b["benchmark_id"] for b in payload["benchmarks"]}
    assert seen == {panel.load_bundle(p)["benchmark_id"] for p, _ in trust}
    for bench in payload["benchmarks"]:
        targets = panel.bundle_targets(panel.load_bundle(ROOT / bench["bench_file"]))
        scored_names = {c["compound"] for c in bench["compounds"]}
        refused_names = {r["compound"] for r in bench["refused_compounds"]}
        assert scored_names | refused_names == set(targets), bench["benchmark_id"]
        assert not (scored_names & refused_names)
        for r in bench["refused_compounds"]:
            assert r["reason"].strip()


def test_markdown_renders(small_artifact):
    text = unc.render_markdown(small_artifact)
    assert "honest literature coverage" in text
    assert "## Priors" in text
    assert unc.NO_UNCERTAINTY not in text.split("## Priors")[0]  # only in the table


# ---------------------------------------------------------------------------
# Observable bands apply to HEADSPACE numbers, not to extraction numbers
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "path, family",
    [
        (data_paths.BENCHMARKS_DIR / "hofmann1998_ribose_cysteine_145C_20min_pH5.json", "extraction"),
        (
            data_paths.MAILLARD_PATH_HOLDOUT_DIR / "mp_holdout_glucose_asparagine_180C_Ye2024.json",
            "extraction",
        ),
        (
            data_paths.EXTERNAL_VALIDATION_DIR / "external_validation_bi_2020_raw_pea_hexanal.json",
            "headspace",
        ),
        (data_paths.BENCHMARKS_DIR / "pea_isolate_uht_140C_Trikusuma2019.json", "undeclared"),
    ],
)
def test_quantification_family_reads_the_bundles_own_class(path, family):
    got, why = panel.quantification_family(panel.load_bundle(path))
    assert got == family, why
    assert why


def test_quantification_family_marker_precedence():
    assert panel.quantification_family({"quantification_class": "stable_isotope_dilution_gcms"})[0] == "extraction"
    assert panel.quantification_family(
        {"content_verification": {"quantification_class": "internal standard, SPME-GC-MS"}}
    )[0] == "headspace"
    assert panel.quantification_family({"quantification_class": "gravimetric"})[0] == "undeclared"
    assert panel.quantification_family({})[0] == "undeclared"


def test_extraction_rows_do_not_carry_the_headspace_bands(small_artifact):
    by_id = {b["benchmark_id"]: b for b in small_artifact["benchmarks"]}
    hofmann = by_id["hofmann1998_ribose_cysteine_145C_20min_pH5"]
    assert hofmann["quantification_family"] == "extraction"
    assert all(not c["observable_multipliers_applied"] for c in hofmann["compounds"])
    # sulfur has no sampled fit uncertainty; without the headspace bands the
    # only movement is the shared B1 trunk, so the width is small (or nil)
    for c in hofmann["compounds"]:
        assert c["ci_width_log10"] is None or c["ci_width_log10"] < 0.5
    trikusuma = by_id["pea_isolate_uht_140C_Trikusuma2019"]
    assert trikusuma["quantification_family"] == "undeclared"
    assert all(c["observable_multipliers_applied"] for c in trikusuma["compounds"])
    policy = small_artifact["summary"]["observable_multiplier_policy"]
    assert policy["rows_by_family"]["extraction"] == len(hofmann["compounds"]) + len(
        by_id["mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7"]["compounds"]
    )
    assert "pea_isolate_uht_140C_Trikusuma2019" in policy["undeclared_bundles"]
