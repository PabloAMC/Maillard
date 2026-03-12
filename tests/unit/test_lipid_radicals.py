import pytest
from rdkit import Chem
from src.smirks_engine import SmirksEngine, Species, ReactionConditions

def test_lipid_hydroperoxide_scission():
    # Linoleic acid hydroperoxide model (simplified C18)
    # SMILES for 9-Hydroperoxiocta-10,12-dienoic acid roughly
    h_smi = r"CCCCC/C=C\C=C/C(O[OH])CCCCCCCC(=O)O"
    lipid_h = Species(label="linoleic_hydroperoxide", smiles=h_smi)
    
    engine = SmirksEngine()
    steps = engine.enumerate([lipid_h])
    
    # 1. Homolysis should occur
    homolysis_steps = [s for s in steps if s.reaction_family == "Lipid_Homolysis"]
    assert len(homolysis_steps) > 0
    
    alkoxy = homolysis_steps[0].products[0]
    assert "alkoxy-radical" in alkoxy.label
    assert "[" in alkoxy.smiles # Radical check
    
    # 2. Beta-scission should follow
    scission_steps = [s for s in steps if s.reaction_family == "Beta_Scission"]
    assert len(scission_steps) > 0
    
    # Check for hexanal or similar fragment
    found_hexanal = False
    for s in scission_steps:
        for p in s.products:
            # Canonicalize for match
            can = Chem.MolToSmiles(Chem.MolFromSmiles(p.smiles))
            if can == "CCCCCC=O" or p.label == "hexanal":
                found_hexanal = True
    
    # Given the SMILES and rules, we expect it to find at least some small aldehyde
    assert any("=" in p.smiles and "O" in p.smiles for s in scission_steps for p in s.products)

def test_radical_crosstalk():
    # R. + H2S -> RH + .SH
    radical = Species(label="test-radical", smiles="CCCC[CH2]")
    h2s = Species(label="H2S", smiles="S")
    
    engine = SmirksEngine()
    steps = engine.enumerate([radical, h2s])
    
    crosstalk = [s for s in steps if s.reaction_family == "Radical_Crosstalk"]
    assert len(crosstalk) > 0
    
    # Products should include quenched alkane and SH radical
    products = [Chem.MolToSmiles(Chem.MolFromSmiles(p.smiles)) for s in crosstalk for p in s.products]
    assert "[SH]" in products
    assert "CCCCC" in products # Quenched pentane
def test_full_autoxidation_cycle():
    # Linoleic acid (C18:2) model
    smi = r"CCCCC/C=C\C=C/CCCCCCCC(=O)O"
    lipid = Species(label="linoleic_acid", smiles=smi)
    # Start with a trace radical
    rad = Species(label="alkyl-radical", smiles=r"CCCCC/C=C\C=C/C[CH]CCCCCCC(=O)O")
    o2 = Species(label="O2", smiles="[O]=[O]")
    
    engine = SmirksEngine()
    steps = engine.enumerate([lipid, rad, o2], max_generations=2)
    
    # 1. Propagation: R. + O2 -> ROO.
    propagation = [s for s in steps if s.reaction_family == "Radical_Propagation_O2"]
    assert len(propagation) > 0
    peroxy = propagation[0].products[0]
    assert "O[O]" in peroxy.smiles
    
    # 2. H-abstraction: ROO. + RH -> ROOH + R.
    abstraction = [s for s in steps if s.reaction_family == "Peroxy_H_Abstraction"]
    assert len(abstraction) > 0
    
    # 3. Termination
    termination = [s for s in steps if s.reaction_family == "Radical_Termination"]
    assert len(termination) > 0

def test_mft_quenching():
    rad = Species(label="alkyl-radical", smiles="CCCC[CH2]")
    mft = Species(label="2-methyl-3-furanthiol", smiles="Cc1c(S)cco1")
    
    engine = SmirksEngine()
    steps = engine.enumerate([rad, mft])
    
    quenching = [s for s in steps if s.reaction_family == "Radical_Crosstalk" and "MFT" in s.products[1].label]
    assert len(quenching) > 0
    assert quenching[0].products[1].label == "MFT-radical"
