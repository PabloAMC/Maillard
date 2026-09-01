import sqlite3

import pytest

from pathlib import Path
from src.results_db import (  # noqa: E402
    CURRENT_THERMO_VERSION,
    DEFAULT_DB_PATH,
    LEGACY_THERMO_VERSION,
    ResultsDB,
    is_computed_method,
    method_family,
    method_uncertainty_kcal,
)
from src.barrier_constants import get_barrier  # noqa: E402

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


# ---------------------------------------------------------------------------
# Audit remediation 2026-08-27 (items 3.7 / 3.8 / 3.9)
# ---------------------------------------------------------------------------

THIOL_REACTANTS = ["O=Cc1ccco1", "S", "[HH]"]
THIOL_PRODUCTS = ["O", "SCc1ccco1"]


def test_default_db_path_is_repo_anchored_and_cwd_independent(tmp_path, monkeypatch):
    """The answer must not depend on where the process happened to start."""
    from_root = ResultsDB()
    monkeypatch.chdir(tmp_path)
    from_elsewhere = ResultsDB()

    assert from_root.db_path == from_elsewhere.db_path == DEFAULT_DB_PATH
    assert from_root.db_path.is_absolute()
    assert not (tmp_path / "results" / "maillard_results.db").exists()

    # And the served barrier is identical from both working directories.
    assert from_root.get_best_barrier(THIOL_REACTANTS, THIOL_PRODUCTS, "Thiol_Addition_H2") == \
        from_elsewhere.get_best_barrier(THIOL_REACTANTS, THIOL_PRODUCTS, "Thiol_Addition_H2")


def test_curated_registry_wins_over_non_computed_db_rows(tmp_path):
    """A 'literature_heuristic' row must not override the curated constant."""
    db = ResultsDB(db_path=str(tmp_path / "policy.db"))
    db.add_barrier(THIOL_REACTANTS, THIOL_PRODUCTS, 15.0, "literature_heuristic",
                   family="Thiol_Addition")

    curated, _ = get_barrier("Thiol_Addition_H2")
    barrier, source, _ = db.get_best_barrier(THIOL_REACTANTS, THIOL_PRODUCTS, "Thiol_Addition_H2")

    assert barrier == curated
    assert barrier != 15.0
    assert source == "Heuristic"


def test_computed_rows_may_override_the_curated_registry(tmp_path):
    db = ResultsDB(db_path=str(tmp_path / "computed.db"))
    db.add_barrier(THIOL_REACTANTS, THIOL_PRODUCTS, 21.5, "wB97M-V/def2-tzvp",
                   family="Thiol_Addition")

    barrier, source, sigma = db.get_best_barrier(THIOL_REACTANTS, THIOL_PRODUCTS, "Thiol_Addition_H2")
    assert barrier == 21.5
    assert source.startswith("DB:")
    assert sigma == 1.5


def test_shipped_db_serves_the_curated_thiol_barrier():
    """Regression guard on the tracked DB: the 15.0 heuristic row is inert."""
    db = ResultsDB()
    curated, _ = get_barrier("Thiol_Addition_H2")
    barrier, source, _ = db.get_best_barrier(THIOL_REACTANTS, THIOL_PRODUCTS, "Thiol_Addition_H2")
    assert (barrier, source) == (curated, "Heuristic")


def test_method_family_matching_is_prefix_based_not_exact():
    assert method_family("wB97M-V/def2-tzvp") == "wb97m-v"
    assert method_family("wB97M-V//r2SCAN") == "wb97m-v"
    assert method_family("GFN2-xTB") == "xtb"
    assert method_family("mace-off24") == "mace-off"
    assert method_family("literature_heuristic") is None

    # DFT must carry the SMALLER sigma than xTB (this was inverted).
    assert method_uncertainty_kcal("wB97M-V/def2-tzvp") < method_uncertainty_kcal("xtb")
    assert is_computed_method("wB97M-V/def2-tzvp")
    assert not is_computed_method("literature_heuristic")


def test_dft_outranks_xtb_even_with_a_basis_suffix(tmp_path):
    db = ResultsDB(db_path=str(tmp_path / "pref.db"))
    db.add_barrier(["C"], ["O"], 40.0, "xtb", family="test")
    db.add_barrier(["C"], ["O"], 30.0, "wB97M-V/def2-tzvp", family="test")

    res = db.find_barrier(["C"], ["O"])
    assert res["method"] == "wB97M-V/def2-tzvp"
    assert method_uncertainty_kcal(res["method"]) == 1.5


def test_lookup_prefers_the_newest_row_within_a_method_family(tmp_path):
    db = ResultsDB(db_path=str(tmp_path / "newest.db"))
    db.add_barrier(["C"], ["O"], 10.0, "wB97M-V/def2-tzvp", family="test")
    db.add_barrier(["C"], ["O"], 20.0, "wB97M-V/def2-tzvp", family="test")

    # Force distinct timestamps (add_barrier stamps to the second).
    with sqlite3.connect(db.db_path) as conn:
        conn.execute("UPDATE barriers SET timestamp = '2020-01-01 00:00:00' WHERE delta_g_kcal = 10.0")
        conn.execute("UPDATE barriers SET timestamp = '2026-01-01 00:00:00' WHERE delta_g_kcal = 20.0")

    assert db.find_barrier(["C"], ["O"])["delta_g_kcal"] == 20.0


def test_thermo_version_filter_excludes_pre_fix_rows(tmp_path):
    db = ResultsDB(db_path=str(tmp_path / "version.db"))
    db.add_barrier(["C"], ["O"], 12.0, "wB97M-V/def2-tzvp", family="test",
                   thermo_version=LEGACY_THERMO_VERSION)

    assert db.find_barrier(["C"], ["O"]) is None
    assert db.find_barrier(["C"], ["O"], thermo_version=None)["delta_g_kcal"] == 12.0

    db.add_barrier(["C"], ["O"], 18.0, "wB97M-V/def2-tzvp", family="test")
    assert db.find_barrier(["C"], ["O"])["delta_g_kcal"] == 18.0


def test_shipped_db_pre_qrrho_rows_are_excluded_but_preserved():
    db = ResultsDB()
    with sqlite3.connect(db.db_path) as conn:
        rows = conn.execute(
            "SELECT thermo_version, count(*) FROM barriers GROUP BY thermo_version"
        ).fetchall()
    tally = dict(rows)
    # Every shipped row predates the QRRHO fix and is therefore excluded from
    # lookups; none is served, none was destroyed.
    assert tally.get(CURRENT_THERMO_VERSION, 0) == 0
    assert tally.get(LEGACY_THERMO_VERSION, 0) > 0


def test_smiles_are_canonicalized_on_write_and_lookup(tmp_path):
    db = ResultsDB(db_path=str(tmp_path / "canon.db"))
    # Non-canonical spelling on write ...
    db.add_barrier(["OCC(O)C(O)C(O)/C=N/CC(=O)O"], ["OCC(O)C(O)C(=O)CNCC(=O)O"],
                   25.0, "wB97M-V/def2-tzvp", family="Amadori_Rearrangement")
    # ... canonical spelling on lookup, and vice versa.
    res = db.find_barrier(["O=C(O)C/N=C/C(O)C(O)C(O)CO"], ["O=C(O)CNCC(=O)C(O)C(O)CO"])
    assert res is not None and res["delta_g_kcal"] == 25.0

    with sqlite3.connect(db.db_path) as conn:
        count = conn.execute("SELECT count(*) FROM reactions").fetchone()[0]
    assert count == 1


def test_shipped_db_has_no_non_canonical_reactions():
    from src.results_db import canonicalize_smiles
    import json as _json

    db = ResultsDB()
    with sqlite3.connect(db.db_path) as conn:
        rows = conn.execute("SELECT id, reactants_json, products_json FROM reactions").fetchall()
    for rid, r_json, p_json in rows:
        for payload in (r_json, p_json):
            values = _json.loads(payload)
            assert values == sorted(canonicalize_smiles(v) for v in values), f"reaction {rid}"


def test_barrier_rows_record_their_evaluation_temperature(tmp_path):
    db = ResultsDB(db_path=str(tmp_path / "temp.db"))
    db.add_barrier(["C"], ["O"], 22.0, "wB97M-V/def2-tzvp", family="test", temperature_k=423.15)
    assert db.find_barrier(["C"], ["O"])["temperature_k"] == 423.15


def test_non_numeric_delta_g_is_rejected(tmp_path):
    db = ResultsDB(db_path=str(tmp_path / "typing.db"))
    with pytest.raises(ValueError):
        db.add_barrier(["C"], ["O"], "n/a", "wB97M-V/def2-tzvp", family="test")

    # SQLite's dynamic typing lets a bad value in through raw SQL; the reader
    # must skip it rather than crash or coerce it.
    db.add_barrier(["C"], ["O"], 19.0, "wB97M-V/def2-tzvp", family="test")
    with sqlite3.connect(db.db_path) as conn:
        reaction_id = conn.execute("SELECT id FROM reactions").fetchone()[0]
        conn.execute(
            "INSERT INTO barriers (reaction_id, delta_g_kcal, method, converged, timestamp,"
            " temperature_k, thermo_version) VALUES (?, 'not-a-number', 'wB97M-V', 1,"
            " '2027-01-01 00:00:00', 298.15, ?)",
            (reaction_id, CURRENT_THERMO_VERSION),
        )

    res = db.find_barrier(["C"], ["O"])
    assert res is not None and res["delta_g_kcal"] == 19.0


def test_find_barrier_default_preference_is_not_a_mutable_default(tmp_path):
    db = ResultsDB(db_path=str(tmp_path / "mutable.db"))
    db.add_barrier(["C"], ["O"], 30.0, "wB97M-V/def2-tzvp", family="test")

    preference = ["xtb"]
    db.find_barrier(["C"], ["O"], method_preference=preference)
    preference.append("hf")
    # Default lookup must be unaffected by any caller's list mutation.
    assert db.find_barrier(["C"], ["O"])["method"] == "wB97M-V/def2-tzvp"
    import inspect
    default = inspect.signature(ResultsDB.find_barrier).parameters["method_preference"].default
    assert default is None


def test_connections_are_closed(tmp_path):
    db = ResultsDB(db_path=str(tmp_path / "close.db"))
    with db._get_connection() as conn:
        conn.execute("SELECT 1")
    with pytest.raises(sqlite3.ProgrammingError):
        conn.execute("SELECT 1")


if __name__ == "__main__":
    # For manual debugging
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        test_results_db_basic(Path(td))
