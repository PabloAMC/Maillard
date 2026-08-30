import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import build_benchmark_index
from src.presentation import render_benchmark_index_markdown


def test_benchmark_index_marks_matrix_only_scope_gaps_explicitly():
    entries = build_benchmark_index([
        ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json",
        ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json",
        ROOT / "data" / "benchmarks" / "soy_isolate_40C_PratapSingh2021.json",
    ])
    by_id = {entry.benchmark_id: entry for entry in entries}

    assert by_id["cys_ribose_140C_Hofmann1998"].execution_path == "free_precursor"
    assert by_id["cys_ribose_140C_Hofmann1998"].benchmark_engine == "fast_observable"
    assert by_id["cys_ribose_140C_Hofmann1998"].cantera_role == "diagnostic_reference_only"
    assert by_id["cys_ribose_140C_Hofmann1998"].thermodynamic_gating_policy == "diagnostic_only"
    assert by_id["pea_isolate_40C_PratapSingh2021"].execution_path == "matrix_only"
    assert by_id["pea_isolate_40C_PratapSingh2021"].benchmark_engine == "matrix_intake_headspace"
    assert by_id["pea_isolate_40C_PratapSingh2021"].cantera_role == "not_authoritative"
    assert by_id["pea_isolate_40C_PratapSingh2021"].thermodynamic_gating_policy == "not_applicable"
    assert by_id["pea_isolate_40C_PratapSingh2021"].supported is True
    assert by_id["pea_isolate_40C_PratapSingh2021"].strict_ready is False
    assert by_id["pea_isolate_40C_PratapSingh2021"].status in {"pass", "pass-no-ranking", "scale-gap", "ranking-gap"}
    assert by_id["soy_isolate_40C_PratapSingh2021"].execution_path == "matrix_only"
    assert by_id["soy_isolate_40C_PratapSingh2021"].supported is True
    assert by_id["soy_isolate_40C_PratapSingh2021"].strict_ready is False
    assert by_id["soy_isolate_40C_PratapSingh2021"].status in {"pass", "pass-no-ranking", "scale-gap", "ranking-gap"}


def test_benchmark_index_markdown_exposes_execution_path():
    entries = build_benchmark_index([
        ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json",
        ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json",
        ROOT / "data" / "benchmarks" / "soy_isolate_40C_PratapSingh2021.json",
    ])

    markdown = render_benchmark_index_markdown(entries)

    assert "Benchmark Index" in markdown
    assert "Execution Path" in markdown
    assert "Engine" in markdown
    assert "Cantera Role" in markdown
    assert "Thermo Policy" in markdown
    assert "matrix_only" in markdown


def test_benchmark_index_includes_matrix_precursor_augmented_candidates():
    entries = build_benchmark_index([
        ROOT / "data" / "benchmarks" / "pea_isolate_ribose_cysteine_100C_45min_Internal2026.json",
        ROOT / "data" / "benchmarks" / "soy_isolate_ribose_cysteine_100C_45min_Internal2026.json",
    ])
    by_id = {entry.benchmark_id: entry for entry in entries}

    assert by_id["pea_isolate_ribose_cysteine_100C_45min_Internal2026"].execution_path == "matrix_precursor_augmented"
    assert by_id["pea_isolate_ribose_cysteine_100C_45min_Internal2026"].strict_ready is False
    assert by_id["soy_isolate_ribose_cysteine_100C_45min_Internal2026"].execution_path == "matrix_precursor_augmented"
    assert by_id["soy_isolate_ribose_cysteine_100C_45min_Internal2026"].strict_ready is False