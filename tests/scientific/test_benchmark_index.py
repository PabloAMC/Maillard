import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import build_benchmark_index, render_benchmark_index_markdown


def test_benchmark_index_marks_matrix_only_scope_gaps_explicitly():
    entries = build_benchmark_index([
        ROOT / "data" / "benchmarks" / "cys_glucose_150C_Farmer1999.json",
        ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json",
    ])
    by_id = {entry.benchmark_id: entry for entry in entries}

    assert by_id["cys_glucose_150C_Farmer1999"].execution_path == "free_precursor"
    assert by_id["pea_isolate_40C_PratapSingh2021"].execution_path == "matrix_only"
    assert by_id["pea_isolate_40C_PratapSingh2021"].supported is False
    assert by_id["pea_isolate_40C_PratapSingh2021"].status == "unsupported"


def test_benchmark_index_markdown_exposes_execution_path():
    entries = build_benchmark_index([
        ROOT / "data" / "benchmarks" / "cys_glucose_150C_Farmer1999.json",
        ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json",
    ])

    markdown = render_benchmark_index_markdown(entries)

    assert "Benchmark Index" in markdown
    assert "Execution Path" in markdown
    assert "matrix_only" in markdown