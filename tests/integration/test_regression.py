import pytest
import json
import os
import pandas as pd
import numpy as np
from pathlib import Path
from src.smirks_engine import SmirksEngine, Species
from src.cantera_export import CanteraExporter
from src.kinetics import KineticsEngine
from src.results_db import ResultsDB
from src.conditions import ReactionConditions

@pytest.fixture
def regression_data():
    path = Path("data/lit/canonical_systems.json")
    with open(path, "r") as f:
        return json.load(f)

@pytest.fixture
def results_db():
    return ResultsDB(db_path="results/maillard_results.db")

NAME_MAP = {
    "O=CC(O)C(O)C(O)CO": "ribose",
    "O=CC(O)C(O)C(O)C(O)CO": "glucose",
    "NCC(=O)O": "glycine",
    "NCC(O)=O": "glycine",
    "NC(CS)C(=O)O": "cysteine",
    "CC(C)CC(N)C(=O)O": "leucine",
    "O": "water",
    "O=C=O": "CO2",
    "N": "ammonia",
    "S": "H2S",
    "[HH]": "H2",
    "Cc1nccs1": "2-methylthiazole",
    "SCc1ccco1": "2-furfurylthiol",
    "O=Cc1ccco1": "furfural",
    "OCC1=CC=C(C=O)O1": "HMF",
    "CC(C)CC=O": "3-methylbutanal",
    "Cc1cnc(C)cn1": "2,5-dimethylpyrazine",
    "CC=O": "acetaldehyde"
}

LOOKUP = {
    "ribose": "O=CC(O)C(O)C(O)CO",
    "glucose": "O=CC(O)C(O)C(O)C(O)CO",
    "glycine": "NCC(=O)O",
    "cysteine": "NC(CS)C(=O)O",
    "leucine": "CC(C)CC(N)C(=O)O"
}

@pytest.mark.regression
@pytest.mark.parametrize("system_key", ["ribose_cysteine", "glucose_glycine", "ribose_cysteine_leucine"])
def test_canonical_systems(system_key, regression_data, results_db):
    data = regression_data[system_key]
    precursors = data["precursors"]
    targets = data["targets"]
    
    # 1. Setup Engine and Conditions
    cond = ReactionConditions(temperature_celsius=150.0)
    engine_smirks = SmirksEngine(conditions=cond)
    
    # 2. Map Precursors
    precursor_objs = []
    for name, conc in precursors.items():
        smi = LOOKUP.get(name.lower(), name)
        precursor_objs.append(Species(label=name, smiles=smi))
        
    # 3. Discover Network
    steps = engine_smirks.enumerate(precursor_objs, max_generations=3)
    assert len(steps) > 0
    
    # 4. Export to Cantera
    exporter = CanteraExporter()
    for step in steps:
        reactants = [s.smiles for s in step.reactants]
        products = [s.smiles for s in step.products]
        
        # Dual-lookup barrier
        barrier_kcal, _, _ = results_db.get_best_barrier(reactants, products, step.reaction_family)
        
        # Add species with names
        for s in step.reactants + step.products:
            name = NAME_MAP.get(s.smiles, s.label)
            exporter.add_species(s.smiles, name=name)
        
        try:
            exporter.add_reaction(reactants, products, barrier_kcal)
        except ValueError:
            continue
            
    assert len(exporter.reactions) > 0
    
    mech_path = f"tests/temp_{system_key}_mech.yaml"
    exporter.export_yaml(mech_path)
    
    # 5. Simulate
    kinetics = KineticsEngine(temperature_k=423.15) # 150 C
    
    init_state = {}
    for name, conc in precursors.items():
        # Find the name used in the mechanism for this precursor
        found_name = None
        smi = LOOKUP.get(name.lower(), name)
        for _, s_info in exporter.species.items():
            if s_info["smiles"] == smi:
                found_name = s_info["name"]
                break
        if found_name:
            init_state[found_name] = conc
            
    assert len(init_state) == len(precursors)
    
    # Increase time to allow cascades to finish
    simulation_duration = 3600 # 1 hour
    results = kinetics.simulate_network_cantera(mech_path, init_state, (0, simulation_duration))
    
    # 6. Verify Targets
    df = pd.DataFrame(results)
    final_concs = {k: v for k, v in results.items() if k not in ["time", "temperature"] and not k.endswith("_X")}
    
    # Filter out precursors and inert small molecules for "volatiles" ranking
    inerts = ["water", "CO2", "H2", "ammonia", "H2S"]
    inerts_lower = [i.lower() for i in inerts]
    precursor_names_lower = [p.lower() for p in init_state.keys()]
    
    # We also exclude heavy non-volatile intermediates from "aroma" ranking
    # such as Schiff bases, Amadori products, and deoxyosones.
    def is_volatile(name):
        name_lower = name.lower()
        if name_lower in inerts_lower: return False
        if name_lower in precursor_names_lower: return False
        if "base" in name_lower: return False
        if "amadori" in name_lower: return False
        if "osone" in name_lower: return False
        # Also exclude common non-aroma fragments
        if name_lower in ["dehydroalanine", "glyceraldehyde"]: return False
        return True

    volatiles = {k: v[-1] for k, v in final_concs.items() if is_volatile(k)}
    
    # Sort volatiles by final concentration
    sorted_volatiles = sorted(volatiles.items(), key=lambda x: x[1], reverse=True)
    # Characteristic volatiles might not be the absolute most abundant (fragments usually are)
    # but they should be in the top 10.
    top_n_names = [name for name, conc in sorted_volatiles[:10]]
    
    # Log for debugging if test fails
    print(f"\nSystem: {system_key}")
    print(f"Top 10 Volatiles: {sorted_volatiles[:10]}")
    
    # Assert at least one target is in top N
    found_any = any(t in top_n_names for t in targets)
    assert found_any, f"None of targets {targets} found in top {len(top_n_names)} volatiles: {top_n_names}. volatiles map count: {len(volatiles)}"
    
    # Cleanup
    if os.path.exists(mech_path):
        os.remove(mech_path)
