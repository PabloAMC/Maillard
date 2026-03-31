import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = ROOT / "scripts" / "generators" / "generate_benchmark_coverage_gaps.py"


def _load_script_module():
    spec = importlib.util.spec_from_file_location("generate_benchmark_coverage_gaps", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_benchmark_coverage_gap_report_flags_missing_external_matrix_meaty_positive():
    module = _load_script_module()
    rows = module._build_rows()

    target = next(
        row for row in rows
        if row["dimension"] == "scientific_gap" and row["category"] == "external_matrix_meaty_positive"
    )

    assert target["status"] == "gap"
    assert target["benchmark_count"] == 0