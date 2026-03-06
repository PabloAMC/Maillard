"""
tests/test_pipeline.py — Test suite for Phase 7.1 end-to-end wiring.

Verifies:
1. `precursor_resolver` correctly maps string names to SMILES.
2. `Recommender.predict_from_steps` correctly finds the min-max bottleneck
   barrier via graph traversal on generated steps.
"""

import pytest
from pathlib import Path

from src.precursor_resolver import resolve, resolve_many
from src.pathway_extractor import Species, ElementaryStep
from src.recommend import Recommender

def test_resolver_exact_match():
    sp = resolve("L-Cysteine")
    assert sp.label == "L-Cysteine"
    assert "S" in sp.smiles

def test_resolver_fuzzy_match():
    sp1 = resolve("ribose")
    sp2 = resolve("d-ribose")
    assert sp1.label == "D-Ribose"
    assert sp1.smiles == sp2.smiles

def test_resolver_multiple():
    sps = resolve_many(["glucose", "glycine", "thiamine"])
    assert len(sps) == 3
    assert sps[0].label == "D-Glucose"
    assert sps[1].label == "Glycine"
    assert sps[2].label == "Thiamine (Vitamin B1)"

def test_resolver_unknown():
    with pytest.raises(ValueError, match="Unknown precursor"):
        resolve("nonexistent_sugar")

def test_recommender_predict_from_steps():
    # We mock an empty results file since we rely on dynamic barriers passed in
    mock_results = Path("mock_results.json")
    if not mock_results.exists():
        mock_results.write_text("[]")
        
    rec = Recommender(mock_results)
    
    # We want to form 2-Methyl-3-furanthiol (MFT), canonical SMILES matches desirable targets
    # Let's use the exact SMILES from the desirable_targets.yml: "Cc1occc1S"
    MFT = Species("MFT", "Cc1occc1S")
    FURFURAL = Species("Furfural", "O=Cc1ccco1") # not a target, intermediate
    RIBOSE = Species("Ribose", "OC[C@H]1OC(O)[C@H](O)[C@@H]1O")
    CYS = Species("Cysteine", "N[C@@H](CS)C(=O)O")
    
    # Pathway 1: Ribose -> Furfural -> MFT
    # Barriers: Ribose->Furfural = 30.0, Furfural->MFT = 15.0
    # Min-max barrier required to reach MFT = max(30, 15) = 30.0
    
    step1 = ElementaryStep(reactants=[RIBOSE], products=[FURFURAL])
    step2 = ElementaryStep(reactants=[FURFURAL, CYS], products=[MFT])
    
    # Build dictionary
    s1_key = f"{RIBOSE.smiles}->{FURFURAL.smiles}"
    s2_key = f"{CYS.smiles}+{FURFURAL.smiles}->{MFT.smiles}"
    # Python sorting makes C comes before O, we must sort the smiles:
    s2_key = f"{sorted([CYS.smiles, FURFURAL.smiles])[0]}+{sorted([CYS.smiles, FURFURAL.smiles])[1]}->{MFT.smiles}"
    
    barriers = {
        s1_key: 30.0,
        s2_key: 15.0
    }
    
    # Also add a slower pathway to a competing target like hexanal to test sorting
    HEXANAL = Species("Hexanal", "CCCCCC=O")
    LINOLEIC = Species("Linoleic", "CCCCC=CCC=CCCCCCCCC(=O)O") # Mock pre-cursor
    step3 = ElementaryStep(reactants=[LINOLEIC], products=[HEXANAL])
    s3_key = f"{LINOLEIC.smiles}->{HEXANAL.smiles}"
    barriers[s3_key] = 45.0
    
    steps = [step1, step2, step3]
    initial_pool = {RIBOSE.smiles: 1.0, CYS.smiles: 1.0, LINOLEIC.smiles: 1.0}
    
    predictions = rec.predict_from_steps(steps, barriers, initial_pool)["targets"]
    
    # Assertions
    assert len(predictions) == 2
    
    # Should be sorted by min-max span: MFT (30.0) first, then Hexanal (45.0)
    p1 = predictions[0]
    p2 = predictions[1]
    
    assert p1["name"] == "2-Methyl-3-furanthiol (MFT)"
    assert p1["span"] == 30.0
    assert p1["type"] == "desirable"
    
    assert p2["name"] == "Hexanal"
    assert p2["span"] == 45.0
    assert p2["type"] == "competing"

    # Cleanup mock if needed
    if mock_results.exists():
        mock_results.unlink()
