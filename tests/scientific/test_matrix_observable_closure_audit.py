import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import build_matrix_observable_closure_audit
from src.presentation import render_matrix_observable_closure_audit_markdown


def test_matrix_observable_closure_audit_labels_transfer_vs_mechanistic_actions():
    payload = build_matrix_observable_closure_audit([
        ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json",
        ROOT / "data" / "benchmarks" / "pea_isolate_uht_140C_Trikusuma2019.json",
        ROOT / "data" / "benchmarks" / "soy_isolate_40C_PratapSingh2021.json",
        ROOT / "data" / "benchmarks" / "pea_isolate_ribose_cysteine_100C_45min_Internal2026.json",
        ROOT / "data" / "benchmarks" / "soy_isolate_ribose_cysteine_100C_45min_Internal2026.json",
    ])

    by_id = {row["benchmark_id"]: row for row in payload["benchmarks"]}
    pea_candidate = by_id["pea_isolate_ribose_cysteine_100C_45min_Internal2026"]
    compounds = {row["compound"]: row for row in pea_candidate["compounds"]}
    watchlist = {row["benchmark_id"]: row for row in payload["mechanistic_refinement_watchlist"]}

    assert compounds["2-furfurylthiol"]["closure_action"] == "class_level_transfer_acceptable"
    assert compounds["Hexanal"]["closure_action"] == "mechanistic_blocker"
    assert watchlist["pea_isolate_ribose_cysteine_100C_45min_Internal2026"]["target_compounds"] == ["Hexanal", "Nonanal"]
    assert "named adverse-marker closure" in watchlist["pea_isolate_ribose_cysteine_100C_45min_Internal2026"]["expected_decision_change"]

    markdown = render_matrix_observable_closure_audit_markdown(payload)
    assert "Matrix Observable Closure Audit" in markdown
    assert "class_level_transfer_acceptable" in markdown
    assert "mechanistic_blocker" in markdown
    assert "Mechanistic Refinement Watchlist" in markdown