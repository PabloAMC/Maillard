from src.cantera_export import CanteraExporter
import sqlite3
import json

db_path = "results/maillard_results.db"
conn = sqlite3.connect(db_path)
cursor = conn.cursor()
cursor.execute("SELECT id, reactants_json, products_json FROM reactions")
rows = cursor.fetchall()
conn.close()

exporter = CanteraExporter()
for rid, r_json, p_json in rows:
    reactants = json.loads(r_json)
    products = json.loads(p_json)
    
    # Let's inspect BEFORE adding to exporter
    r_comp = {}
    for r in reactants:
        name = exporter.add_species(r)
        for k, v in exporter.species[r]["composition"].items():
            r_comp[k] = r_comp.get(k, 0) + v
            
    p_comp = {}
    for p in products:
        name = exporter.add_species(p)
        for k, v in exporter.species[p]["composition"].items():
            p_comp[k] = p_comp.get(k, 0) + v
            
    if r_comp != p_comp:
        print(f"Reaction ID {rid} IS UNBALANCED:")
        print(f"Reactants: {reactants} -> {r_comp}")
        print(f"Products: {products} -> {p_comp}")
        
    exporter.add_reaction(reactants, products, 10.0)
    
for s in exporter.species.values():
    if s["name"] == "S_13" or s["name"] == "S_0" or s["name"] == "S_14":
        print(f"{s['name']} comp = {s['composition']}")
