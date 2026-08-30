import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import (
    get_matrix_only_target_snapshot_exclusions,
    snapshot_all_benchmark_targets,
    snapshot_benchmark_targets,
)
from src.presentation import render_benchmark_targets_markdown
from src.recommend import _NON_OBSERVABLE_KAW_THRESHOLD


def test_benchmark_targets_snapshot_contains_headspace_metadata():
    rows = snapshot_benchmark_targets(ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json")

    assert rows
    furfural = next(row for row in rows if row.target_name == "Furfural")
    assert furfural.target_type in {"desirable", "competing", "toxic"}
    assert furfural.target_class in {"severity_markers", "furans_furanones", "unknown"}
    assert furfural.evidence_state in {"conditional_calibration", "transferred_prior", "still_missing"}
    assert furfural.headspace_class in {"observable", "low_headspace", "assumed_observable"}
    assert furfural.henry_source_name
    assert furfural.proxy_ppb >= furfural.predicted_ppb
    assert 0.0 <= furfural.observable_ratio <= 1.0


def test_benchmark_targets_markdown_reports_low_headspace_count():
    # Second file retargeted from acrylamide_asparagine_glucose_Parker2012 (quarantined
    # 2026-08-26, dead DOI + tightest-in-collection tolerance) to the rebuilt Bolton 1994
    # thiamine benchmark. This test only needs two renderable free-precursor benchmarks.
    rows = snapshot_all_benchmark_targets([
        ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json",
        ROOT / "data" / "benchmarks" / "thiamine_cys_glucose_120C_Bolton1994.json",
    ])

    markdown = render_benchmark_targets_markdown(rows)

    assert "Benchmark Targets" in markdown
    assert "Evidence State" in markdown
    assert "Panel Class" in markdown
    assert "Headspace" in markdown
    assert "Proxy ppb" in markdown
    assert "Obs/Proxy" in markdown
    assert "Low-headspace rows:" in markdown


def test_benchmark_targets_markdown_reports_matrix_only_exclusions():
    rows = snapshot_all_benchmark_targets([
        ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json",
        ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json",
        ROOT / "data" / "benchmarks" / "soy_isolate_40C_PratapSingh2021.json",
    ])

    markdown = render_benchmark_targets_markdown(
        rows,
        excluded_benchmark_ids=get_matrix_only_target_snapshot_exclusions([
            ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json",
            ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json",
            ROOT / "data" / "benchmarks" / "soy_isolate_40C_PratapSingh2021.json",
        ]),
    )

    assert "Excluded matrix-only benchmarks:" in markdown
    assert "pea_isolate_40C_PratapSingh2021" in markdown
    assert "soy_isolate_40C_PratapSingh2021" in markdown


def test_headspace_classification_and_suppression_follow_the_shipped_kaw_gate(tmp_path):
    """The low-headspace gate is classified from the shipped Kaw table and suppresses output.

    History, because this test has been retargeted twice and then restructured:

    1. It began on ``cys_glucose_150C_Farmer1999`` and was moved to
       ``acrylamide_asparagine_glucose_Parker2012`` when that file was quarantined. Parker was
       itself quarantined on 2026-08-26 (dead DOI, no identifiable source, and the tightest
       validation contract in the collection), so the fixture is gone again.
    2. More importantly, the assertion it made -- ``HMF is low_headspace`` -- is no longer true
       of the shipped data. The Kaw table was rebuilt against Sander 5.0.0 and HMF now carries
       5e-8, above the 1e-8 observability gate; **no benchmark in the panel currently produces
       a single ``low_headspace`` row.** Re-pinning a specific compound to a specific class
       would therefore encode a data decision this test does not own.

    So the test now asserts the *mechanism* rather than one compound's classification: that the
    class is derived from the shipped Kaw against the shipped threshold, and that whenever a
    row does land below the gate its observable output is suppressed well below its proxy. It
    keeps working whichever side of the gate a given compound's Kaw ends up on.

    The synthetic hexose formulation below replaces the quarantined Parker fixture: HMF is only
    reachable from a hexose, so the surviving pentose panel benchmarks (cysteine + ribose /
    xylose) have no HMF target row at all. It carries no reference values, only conditions.
    """
    bench_path = tmp_path / "synthetic_hexose_headspace_probe.json"
    bench_path.write_text(
        json.dumps(
            {
                "benchmark_id": "synthetic_hexose_headspace_probe",
                "source_doi": "synthetic_mechanism_probe",
                "precursors": {
                    "L-Asparagine": {"concentration_mM": 10.0},
                    "D-Glucose": {"concentration_mM": 10.0},
                },
                "conditions": {
                    "temp_C": 180.0,
                    "ph": 6.8,
                    "water_activity": 0.3,
                    "time_min": 20.0,
                },
                "metadata": {
                    "tier": "DIAGNOSTIC",
                    "family": "headspace_mechanism_probe",
                    "execution_path": "free_precursor",
                    "notes": (
                        "Synthetic hexose formulation used only to obtain Furfural and HMF "
                        "target rows. It declares no measured or reference volatiles and is "
                        "never scored against anything."
                    ),
                },
            }
        ),
        encoding="utf-8",
    )

    rows = snapshot_benchmark_targets(bench_path)
    by_name = {row.target_name: row for row in rows}

    # Both markers must be reachable from a hexose at all -- this is the reachability
    # regression the Parker fixture used to carry.
    furfural = by_name["Furfural"]
    hmf = by_name["5-Hydroxymethylfurfural (HMF)"]

    for row in (furfural, hmf):
        assert row.henry_kaw_25c is not None, f"{row.target_name} has no Henry entry"
        expected_class = (
            "observable" if row.henry_kaw_25c >= _NON_OBSERVABLE_KAW_THRESHOLD else "low_headspace"
        )
        assert row.headspace_class == expected_class, (
            f"{row.target_name}: Kaw {row.henry_kaw_25c:g} vs gate "
            f"{_NON_OBSERVABLE_KAW_THRESHOLD:g} implies {expected_class}, got {row.headspace_class}"
        )

    # The suppression invariant, asserted over every row rather than one hand-picked pair:
    # the proxy is the unsuppressed quantity, so it can never sit below the observable output,
    # and anything below the gate must actually be pushed down by a wide margin.
    for row in rows:
        assert row.proxy_ppb >= row.predicted_ppb, (
            f"{row.target_name}: observable output {row.predicted_ppb:g} exceeds its own "
            f"unsuppressed proxy {row.proxy_ppb:g}"
        )
        if row.headspace_class == "low_headspace" and row.proxy_ppb > 0.0:
            assert row.predicted_ppb < row.proxy_ppb / 20.0, (
                f"{row.target_name} is classified low_headspace but its output is barely "
                f"suppressed ({row.predicted_ppb:g} vs proxy {row.proxy_ppb:g})"
            )


def test_matrix_precursor_augmented_targets_are_included_in_snapshot_reports():
    rows = snapshot_all_benchmark_targets([
        ROOT / "data" / "benchmarks" / "pea_isolate_ribose_cysteine_100C_45min_Internal2026.json",
        ROOT / "data" / "benchmarks" / "soy_isolate_ribose_cysteine_100C_45min_Internal2026.json",
    ])

    benchmark_ids = {row.benchmark_id for row in rows}
    assert "pea_isolate_ribose_cysteine_100C_45min_Internal2026" in benchmark_ids
    assert "soy_isolate_ribose_cysteine_100C_45min_Internal2026" in benchmark_ids