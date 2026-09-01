import pytest

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.generators.generate_validation_figures import _bench_label, _build_payload


def test_validation_overview_payload_includes_quantitative_matrix_benchmarks():
    payload = _build_payload()

    quantitative = {row["benchmark_id"]: row for row in payload["quantitative_benchmarks"]}
    plotted_ids = {row["benchmark_id"] for row in payload["quantitative_points"]}

    assert "pea_isolate_uht_140C_Trikusuma2019" in quantitative
    assert quantitative["pea_isolate_uht_140C_Trikusuma2019"]["execution_path"] == "matrix_only"
    assert "pea_isolate_uht_140C_Trikusuma2019" in plotted_ids
    assert any(row["execution_path"] == "matrix_only" for row in payload["quantitative_points"])


def test_validation_overview_payload_reports_current_tolerance_closure():
    payload = _build_payload()
    ratios = [float(row["max_ratio"] or 0.0) for row in payload["quantitative_benchmarks"]]
    experimental_ratios = [
        float(row["max_ratio"] or 0.0)
        for row in payload["quantitative_benchmarks"]
        if row["reference_signal_origin"] == "measured_volatiles"
    ]

    assert payload["outside_1_5x_benchmark_count"] == sum(1 for ratio in ratios if ratio > 1.5)
    assert payload["outside_2x_benchmark_count"] == sum(1 for ratio in ratios if ratio > 2.0)
    assert payload["inside_1_5x_benchmark_count"] + payload["outside_1_5x_benchmark_count"] == len(payload["quantitative_benchmarks"])
    assert payload["worst_quantitative_ratio"] == pytest.approx(max(ratios, default=0.0))
    assert payload["experimental_inside_1_5x_benchmark_count"] == sum(1 for ratio in experimental_ratios if ratio <= 1.5)
    assert payload["experimental_outside_1_5x_benchmark_count"] == sum(1 for ratio in experimental_ratios if ratio > 1.5)
    assert payload["experimental_outside_2x_benchmark_count"] == sum(1 for ratio in experimental_ratios if ratio > 2.0)
    assert payload["experimental_worst_quantitative_ratio"] == pytest.approx(max(experimental_ratios, default=0.0))
    # UPDATED 2026-08-27 (cause: reporting-consistency fix, audit remediation Part 2a).
    # `worst_quantitative_ratio` was always a max over ALL quantitative benchmarks, but
    # `worst_quantitative_point` was selected only from the reference_volatiles subset.
    # The two coincided by accident until the Bolton1994 benchmark entered the panel, at
    # which point the artifact reported a worst ratio of 173x next to a worst point of
    # 3.4x. This test pinned the buggy pairing (it asserted the point was Cerny, a
    # reference_volatiles row, while the ratio above it came from the full population).
    # The generator now draws both from the same population and exposes the reference-only
    # view under explicitly named `reference_worst_quantitative_*` keys. The assertions
    # below pin the CONSISTENCY property rather than a specific benchmark id, so the pair
    # can never silently diverge again.
    assert payload["worst_quantitative_population"] == "all_quantitative_benchmarks"
    assert payload["worst_quantitative_point"]["max_ratio"] == pytest.approx(
        payload["worst_quantitative_ratio"]
    ), "worst_quantitative_point and worst_quantitative_ratio must describe the same population"

    reference_ratios = [
        float(row["max_ratio"] or 0.0)
        for row in payload["quantitative_benchmarks"]
        if row["reference_signal_origin"] != "measured_volatiles"
    ]
    assert payload["reference_worst_quantitative_ratio"] == pytest.approx(
        max(reference_ratios, default=0.0)
    )
    assert payload["reference_worst_quantitative_point"]["benchmark_id"] == "thiamine_cys_xylose_145C_Cerny2008"
    assert payload["reference_worst_quantitative_point"]["reference_signal_origin"] == "reference_volatiles"
    assert payload["experimental_worst_quantitative_point"]["benchmark_id"] != "thiamine_cys_xylose_145C_Cerny2008"
    assert payload["experimental_worst_quantitative_point"]["reference_signal_origin"] == "measured_volatiles"


def test_validation_overview_uses_standardized_benchmark_labels():
    assert _bench_label("thiamine_cys_ribose_100C_Hofmann1996") == r"Thiamine + cysteine + ribose, $100\,^{\circ}$C (Hofmann, 1996)"
    assert _bench_label("thiamine_cys_xylose_145C_Cerny2008") == r"Thiamine + cysteine + xylose, $145\,^{\circ}$C (Cerny, 2008) [reference anchor]"
    assert _bench_label("pea_isolate_ribose_cysteine_100C_45min_Internal2026") == r"Pea isolate + ribose + cysteine, $100\,^{\circ}$C (Internal, 2026)"
    assert _bench_label("resconi_2023_pbma_beef_identity_benchmark") == r"PBMA vs beef comparator, $150\,^{\circ}$C (Resconi et al., 2023)"
    # RE-PINNED 2026-08-27 (Wave I): the label now carries a `[QUARANTINED]` suffix.
    # `spi_hvp_xylose_120C_PMC9905368` and its wheat-gluten twin were quarantined as
    # fabricated -- the cited paper reports relative peak areas for glucose/fructose at
    # pH 7.5 and never mentions FFT or MFT. The label entries are KEPT so any surviving
    # reference renders as a name rather than a raw slug, and the suffix is there so a
    # figure or table that somehow still shows one of them cannot look like a live
    # benchmark. See data/benchmarks/quarantined/README.md.
    assert _bench_label("spi_hvp_xylose_120C_PMC9905368") == r"SPI hydrolysate + xylose, $120\,^{\circ}$C (Cho et al., 2023) [QUARANTINED]"
    # RE-PINNED 2026-08-27 (Wave I): see the SPI entry above -- same quarantine.
    assert _bench_label("wheat_gluten_hvp_xylose_120C_PMC9905368") == r"Wheat gluten hydrolysate + xylose, $120\,^{\circ}$C (Cho et al., 2023) [QUARANTINED]"


def test_validation_overview_payload_surfaces_integrated_families_11_to_16():
    payload = _build_payload()

    assert payload["integrated_family_count"] >= 14
    by_slr = {row["slr_family"]: row for row in payload["family_overview"]["families"]}

    assert by_slr["11"]["has_runtime_support"] is True
    assert by_slr["12"]["has_runtime_support"] is True
    assert by_slr["13"]["has_runtime_support"] is True
    assert by_slr["14"]["has_runtime_support"] is True
    assert by_slr["15"]["has_runtime_support"] is True
    assert by_slr["16"]["has_runtime_support"] is True
    assert by_slr["12"]["benchmark_count"] >= 3
    assert by_slr["15"]["benchmark_count"] == 0
    assert by_slr["16"]["benchmark_count"] == 0