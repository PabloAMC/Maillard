import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import (
    BenchmarkSummary,
    audit_all_thermodynamic_gating,
    resolve_thermodynamic_gating_mode,
    thermodynamic_gating_materially_improves,
)
from src.presentation import render_thermodynamic_gating_audit_markdown


def _summary(*, status: str, mae: float, max_ratio: float) -> BenchmarkSummary:
    return BenchmarkSummary(
        benchmark_id="synthetic",
        bench_file=ROOT / "data" / "benchmarks" / "cys_glucose_150C_Farmer1999.json",
        tier="PRIMARY",
        family="free_aa_sulfur",
        execution_path="free_precursor",
        benchmark_engine="fast_observable",
        comparator_signal="predicted_ppb",
        cantera_role="diagnostic_reference_only",
        target_snapshot_policy="included",
        thermodynamic_gating_policy="diagnostic_only",
        supported=True,
        reason=None,
        protein_type="free",
        coverage=1.0,
        matched_compounds=3,
        total_compounds=3,
        pearson_r=0.99,
        mae_ppb=mae,
        max_ratio=max_ratio,
        mean_ratio=max_ratio,
        ranking_status="pass",
        scale_status="pass",
        overall_status=status,
        strict_ready=True,
        blocking_issues=[],
        conditions={},
    )


def test_thermodynamic_gating_materiality_requires_real_improvement():
    baseline = _summary(status="pass", mae=100.0, max_ratio=1.40)
    gated = _summary(status="pass", mae=92.0, max_ratio=1.33)

    assert thermodynamic_gating_materially_improves(baseline, gated) is True


def test_thermodynamic_gating_materiality_rejects_worse_status():
    baseline = _summary(status="pass", mae=100.0, max_ratio=1.40)
    gated = _summary(status="scale-gap", mae=90.0, max_ratio=1.60)

    assert thermodynamic_gating_materially_improves(baseline, gated) is False


def test_thermodynamic_gating_audit_reports_policy_for_real_benchmarks():
    rows = audit_all_thermodynamic_gating([
        ROOT / "data" / "benchmarks" / "cys_glucose_150C_Farmer1999.json",
        ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json",
    ])

    by_id = {row.benchmark_id: row for row in rows}
    assert by_id["cys_glucose_150C_Farmer1999"].recommended_policy in {"diagnostic_only", "benchmark_facing_candidate"}
    assert by_id["pea_isolate_40C_PratapSingh2021"].applicable is False

    markdown = render_thermodynamic_gating_audit_markdown(rows)
    assert "Thermodynamic Gating Audit" in markdown
    assert "Recommended Policy" in markdown


def test_auto_thermodynamic_gating_mode_resolves_to_contract_policy():
    free_bench = {
        "benchmark_id": "cys_glucose_150C_Farmer1999",
        "protein_type": "free",
    }
    matrix_bench = {
        "benchmark_id": "pea_isolate_40C_PratapSingh2021",
        "protein_type": "pea_iso",
    }

    assert resolve_thermodynamic_gating_mode(free_bench, "auto") == "off"
    assert resolve_thermodynamic_gating_mode(matrix_bench, "auto") == "off"


def test_auto_thermodynamic_gating_mode_honors_benchmark_facing_override():
    bench = {
        "benchmark_id": "cys_glucose_150C_Farmer1999",
        "protein_type": "free",
        "metadata": {
            "thermodynamic_gating_policy": "benchmark_facing",
        },
    }

    assert resolve_thermodynamic_gating_mode(bench, "auto") == "on"