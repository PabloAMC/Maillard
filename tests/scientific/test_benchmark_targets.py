import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import render_benchmark_targets_markdown, snapshot_all_benchmark_targets, snapshot_benchmark_targets


def test_benchmark_targets_snapshot_contains_headspace_metadata():
    rows = snapshot_benchmark_targets(ROOT / "data" / "benchmarks" / "cys_glucose_150C_Farmer1999.json")

    assert rows
    furfural = next(row for row in rows if row.target_name == "Furfural")
    assert furfural.target_type in {"desirable", "competing", "toxic"}
    assert furfural.headspace_class in {"observable", "low_headspace", "assumed_observable"}
    assert furfural.henry_source_name


def test_benchmark_targets_markdown_reports_low_headspace_count():
    rows = snapshot_all_benchmark_targets([
        ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json",
        ROOT / "data" / "benchmarks" / "cys_glucose_150C_Farmer1999.json",
    ])

    markdown = render_benchmark_targets_markdown(rows)

    assert "Benchmark Targets" in markdown
    assert "Headspace" in markdown
    assert "Low-headspace rows:" in markdown


def test_benchmark_targets_keep_low_headspace_markers_well_below_observable_outputs():
    rows = snapshot_benchmark_targets(ROOT / "data" / "benchmarks" / "cys_glucose_150C_Farmer1999.json")

    by_name = {row.target_name: row for row in rows}
    furfural = by_name["Furfural"]
    hmf = by_name["5-Hydroxymethylfurfural (HMF)"]

    assert hmf.headspace_class == "low_headspace"
    assert furfural.headspace_class == "observable"
    assert furfural.predicted_ppb > hmf.predicted_ppb * 20.0