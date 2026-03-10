from src.smirks_engine import Species, _thiol_addition

def test_thiohemiacetal_pathway():
    # Test that _thiol_addition generates the 2-step H2S pathway
    furfural = Species(label="furfural", smiles="O=Cc1ccco1")
    h2s = Species(label="H2S", smiles="S")
    water = Species(label="water", smiles="O")
    
    steps = _thiol_addition([furfural, h2s, water])
    
    # Verify Step 1: Thiohemiacetal formation
    step1 = next(s for s in steps if s.reaction_family == "Thiohemiacetal_Formation")
    assert "OC(S)c1ccco1" in [p.smiles for p in step1.products]
    
    # Verify Step 2: Dehydration to FFT
    step2 = next(s for s in steps if s.reaction_family == "Thiol_Dehydration")
    assert "SCc1ccco1" in [p.smiles for p in step2.products]
    assert "[S]" in [p.smiles for p in step2.products]

if __name__ == "__main__":
    test_thiohemiacetal_pathway()
    print("Test passed!")
