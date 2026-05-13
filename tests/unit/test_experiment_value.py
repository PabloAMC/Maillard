"""Unit tests for src.experiment_value."""

from __future__ import annotations

from pathlib import Path

import pytest

from src.experiment_value import (
    CompoundSpec,
    _decision_relevance,
    _envelope_miss_log10,
    _normalise,
    _suggest_template,
    _voi_score,
    build_ranking_payload,
    filter_by_matrix,
    infer_matrix_family,
    load_compound_specs,
    lookup_spec,
    rank_experiments,
    render_markdown,
    write_artifact,
)


def test_normalise_strips_aliases_and_lowercases():
    assert _normalise("2-Methyl-3-furanthiol (MFT)") == "2 methyl 3 furanthiol"
    assert _normalise("2-furfurylthiol") == "2 furfurylthiol"
    assert _normalise("") == ""


def test_load_compound_specs_aliases_resolve():
    specs = load_compound_specs()
    # Both the long and short alias should resolve to the same MFT spec.
    long_form = lookup_spec("2-Methyl-3-furanthiol (MFT)", specs)
    short_form = lookup_spec("MFT", specs)
    bare = lookup_spec("2-methyl-3-furanthiol", specs)
    assert long_form is not None
    assert short_form is long_form
    assert bare is long_form
    assert long_form.odour_threshold_ug_per_kg == pytest.approx(0.0001)


def test_envelope_miss_inside_ci_is_zero():
    assert _envelope_miss_log10(measured=10.0, p5=1.0, p95=100.0) == 0.0


def test_envelope_miss_above_p95():
    miss = _envelope_miss_log10(measured=1000.0, p5=1.0, p95=100.0)
    assert miss == pytest.approx(1.0)


def test_envelope_miss_below_p5():
    miss = _envelope_miss_log10(measured=0.1, p5=1.0, p95=100.0)
    assert miss == pytest.approx(1.0)


def test_decision_relevance_known_odt_is_log_ratio_clipped():
    assert _decision_relevance(measured=10.0, p50=10.0, odt=1.0) == pytest.approx(1.0)
    assert _decision_relevance(measured=1e6, p50=1e6, odt=1.0) == pytest.approx(5.0)
    # Below ODT clipped to 0.5
    assert _decision_relevance(measured=0.01, p50=0.01, odt=1.0) == pytest.approx(0.5)


def test_decision_relevance_no_odt_returns_one():
    assert _decision_relevance(measured=10.0, p50=10.0, odt=None) == 1.0


def test_voi_inside_ci_uses_width_only():
    s = _voi_score(
        inside_ci=True, envelope_miss_log10=0.0, ci_width_log10=2.0, decision_relevance=2.0
    )
    assert s == pytest.approx(0.3 * 2.0 * 2.0)


def test_voi_outside_ci_includes_miss_term():
    s = _voi_score(
        inside_ci=False, envelope_miss_log10=1.0, ci_width_log10=0.5, decision_relevance=2.0
    )
    # miss_term=2.0, width_term=0.15 -> total=(2.15)*2.0=4.3
    assert s == pytest.approx(4.3, rel=1e-6)


def test_suggest_template_safety_marker():
    template, _ = _suggest_template("acrylamide", ci_width_log10=0.1, inside_ci=False)
    assert template == "missing_absolute_anchor"


def test_suggest_template_meaty():
    template, _ = _suggest_template("2-Methyl-3-furanthiol (MFT)", ci_width_log10=0.1, inside_ci=False)
    assert template == "blocking_benchmark_gap"


def test_rank_experiments_sorted_descending_with_full_payload(tmp_path: Path):
    payload = {
        "benchmarks": [
            {
                "benchmark_id": "bench_a",
                "compounds": [
                    {
                        "compound": "acrylamide",
                        "measured_ppb": 1500.0,
                        "predicted_p5": 1.0,
                        "predicted_p50": 5.0,
                        "predicted_p95": 10.0,
                        "inside_ci": False,
                        "ci_width_log10": 1.0,
                    },
                    {
                        "compound": "hexanal",
                        "measured_ppb": 50.0,
                        "predicted_p5": 10.0,
                        "predicted_p50": 30.0,
                        "predicted_p95": 100.0,
                        "inside_ci": True,
                        "ci_width_log10": 1.0,
                    },
                ],
            }
        ],
    }
    ranked = rank_experiments(payload)
    assert [c.compound for c in ranked][0] == "acrylamide"
    assert ranked[0].rank == 1
    assert ranked[0].voi_score > ranked[1].voi_score
    out = build_ranking_payload(ranked)
    md = render_markdown(out)
    assert "acrylamide" in md
    paths = write_artifact(out, output_dir=tmp_path, basename="evr_smoke")
    assert paths["json"].exists() and paths["md"].exists()


# ---------------------------------------------------------------------------
# S23 — matrix-family attribution and filtering
# ---------------------------------------------------------------------------

def test_infer_matrix_family_explicit_protein_type_wins():
    bench = {"benchmark_id": "some_irrelevant_id", "protein_type": "pea_iso"}
    assert infer_matrix_family(bench["benchmark_id"], bench) == "pea_iso"


def test_infer_matrix_family_id_prefix_for_unattributed_benchmark():
    # wheat_gluten and SPI HVP benchmarks ship without protein_type today;
    # the helper must still resolve them via benchmark_id substring matching.
    assert (
        infer_matrix_family("wheat_gluten_hvp_xylose_120C_PMC9905368", {})
        == "wheat_gluten"
    )
    assert (
        infer_matrix_family("spi_hvp_xylose_120C_PMC9905368", {}) == "soy_iso"
    )


def test_infer_matrix_family_free_precursor_and_unknown():
    assert infer_matrix_family("cys_glucose_150C_Farmer1999", {}) == "free"
    assert infer_matrix_family("thiamine_thermal_120C", {}) == "free"
    assert infer_matrix_family("completely_novel_id", {}) == "unknown"


def _three_matrix_payload():
    return {
        "benchmarks": [
            {
                "benchmark_id": "pea_isolate_uht_140C_demo",
                "protein_type": "pea_iso",
                "compounds": [
                    {
                        "compound": "2-furfurylthiol",
                        "measured_ppb": 50.0,
                        "predicted_p5": 1.0,
                        "predicted_p50": 5.0,
                        "predicted_p95": 10.0,
                        "inside_ci": False,
                        "ci_width_log10": 1.0,
                    }
                ],
            },
            {
                "benchmark_id": "wheat_gluten_hvp_xylose_120C_demo",
                "compounds": [
                    {
                        "compound": "methional",
                        "measured_ppb": 200.0,
                        "predicted_p5": 1.0,
                        "predicted_p50": 5.0,
                        "predicted_p95": 20.0,
                        "inside_ci": False,
                        "ci_width_log10": 1.3,
                    }
                ],
            },
            {
                "benchmark_id": "cys_glucose_150C_demo",
                "compounds": [
                    {
                        "compound": "hexanal",
                        "measured_ppb": 30.0,
                        "predicted_p5": 10.0,
                        "predicted_p50": 25.0,
                        "predicted_p95": 80.0,
                        "inside_ci": True,
                        "ci_width_log10": 0.9,
                    }
                ],
            },
        ],
    }


def test_rank_experiments_attaches_matrix_family_and_payload_exposes_it():
    ranked = rank_experiments(_three_matrix_payload())
    families = {c.matrix_family for c in ranked}
    assert families == {"pea_iso", "wheat_gluten", "free"}
    payload = build_ranking_payload(ranked, matrix_filter=["pea_iso"])
    assert payload["matrix_filter"] == ["pea_iso"]
    assert all("matrix_family" in row for row in payload["candidates"])
    md = render_markdown(payload)
    assert "Matrix filter:" in md
    assert "| Matrix |" in md


def test_filter_by_matrix_keeps_only_requested_and_renumbers():
    ranked = rank_experiments(_three_matrix_payload())
    soy_only = filter_by_matrix(ranked, ["wheat_gluten"])
    assert len(soy_only) == 1
    assert soy_only[0].matrix_family == "wheat_gluten"
    assert soy_only[0].rank == 1


def test_filter_by_matrix_empty_filter_returns_all():
    ranked = rank_experiments(_three_matrix_payload())
    assert len(filter_by_matrix(ranked, None)) == len(ranked)
    assert len(filter_by_matrix(ranked, [])) == len(ranked)
