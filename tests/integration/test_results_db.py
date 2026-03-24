from pathlib import Path
from src.results_db import ResultsDB  # noqa: E402

def test_results_db_basic(tmp_path):
    db_file = tmp_path / "test.db"
    db = ResultsDB(db_path=str(db_file))
    
    reactants = ["C", "O2"]
    products = ["CO2"]
    
    # 1. Add a barrier
    db.add_barrier(reactants, products, 35.5, "wB97M-V", family="combustion")
    
    # 2. Find it back
    res = db.find_barrier(reactants, products)
    assert res is not None
    assert res["delta_g_kcal"] == 35.5
    assert res["method"] == "wB97M-V"
    
    # 3. Test method priority
    db.add_barrier(reactants, products, 40.0, "xtb", family="combustion")
    
    # Should prefer wB97M-V over xtb
    res_best = db.find_barrier(reactants, products, method_preference=["wB97M-V", "xtb"])
    assert res_best["method"] == "wB97M-V"
    
    # Should prefer xtb if requested first
    res_xtb = db.find_barrier(reactants, products, method_preference=["xtb", "wB97M-V"])
    assert res_xtb["method"] == "xtb"

def test_results_db_normalization(tmp_path):
    db_file = tmp_path / "test_norm.db"
    db = ResultsDB(db_path=str(db_file))
    
    # Reactants in different order should map to same reaction
    db.add_barrier(["A", "B"], ["C"], 10.0, "method")
    res = db.find_barrier(["B", "A"], ["C"])
    assert res is not None
    assert res["delta_g_kcal"] == 10.0

def test_results_db_empty(tmp_path):
    db_file = tmp_path / "test_empty.db"
    db = ResultsDB(db_path=str(db_file))
    assert db.find_barrier(["X"], ["Y"]) is None


def test_results_db_persists_ml_adoption_decisions(tmp_path):
    db_file = tmp_path / "test_p4.db"
    db = ResultsDB(db_path=str(db_file))

    db.add_ml_adoption_decision(
        {
            "candidate_id": "mace_off_medium",
            "model_family": "mace_off",
            "model_name": "medium",
            "proposed_role": "local_barrier_surrogate",
            "decision": "quarantine",
            "benchmark_set_id": "maillard_reaction_benchmark_v1",
            "coverage_ratio": 0.5,
            "rank_correlation": -1.0,
            "mean_abs_error_kcal": 50.0,
            "max_abs_error_kcal": 90.0,
            "stop_reasons": ["nonphysical_barrier_predictions"],
            "rationale": "Nonphysical barriers on reaction benchmark.",
            "fallback_comparator": "xTB_plus_selective_DFT",
            "benchmark_visible_gap": "barrier ordering",
            "approved_for_default": False,
        }
    )

    rows = db.list_ml_adoption_decisions("mace_off_medium")
    assert len(rows) == 1
    assert rows[0]["decision"] == "quarantine"
    assert rows[0]["stop_reasons"] == ["nonphysical_barrier_predictions"]

if __name__ == "__main__":
    # For manual debugging
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        test_results_db_basic(Path(td))
