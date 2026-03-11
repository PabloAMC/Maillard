import pytest
import sqlite3
import json
from pathlib import Path
from src.cantera_export import CanteraExporter  # noqa: E402

@pytest.mark.skipif(not Path("results/maillard_results.db").exists(), reason="Results DB not found")
def test_database_reaction_balance():
    """Verify that all reactions stored in the SQLite results DB are atom-balanced."""
    db_path = "results/maillard_results.db"
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT id, reactants_json, products_json FROM reactions")
    rows = cursor.fetchall()
    conn.close()

    exporter = CanteraExporter()
    unbalanced = []
    
    for rid, r_json, p_json in rows:
        reactants = json.loads(r_json)
        products = json.loads(p_json)
        
        r_comp = {}
        for r in reactants:
            exporter.add_species(r)
            for k, v in exporter.species[r]["composition"].items():
                r_comp[k] = r_comp.get(k, 0) + v
                
        p_comp = {}
        for p in products:
            exporter.add_species(p)
            for k, v in exporter.species[p]["composition"].items():
                p_comp[k] = p_comp.get(k, 0) + v
                
        if r_comp != p_comp:
            unbalanced.append((rid, r_comp, p_comp))
    
    assert not unbalanced, f"Found {len(unbalanced)} unbalanced reactions in DB: {unbalanced[:3]}..."
