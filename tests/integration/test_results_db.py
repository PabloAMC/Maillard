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

if __name__ == "__main__":
    # For manual debugging
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        test_results_db_basic(Path(td))
