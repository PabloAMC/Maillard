import pytest
from pathlib import Path
from src.pathway_extractor import PathwayExtractor, ElementaryStep, Species

# Define mock test directories
MOCK_DIR = Path(__file__).parent / "mock_rmg_output"
CHEMKIN_DIR = MOCK_DIR / "chemkin"

@pytest.fixture(autouse=True)
def setup_mock_rmg_output(tmp_path):
    """Generates mock RMG output files for testing."""
    chem_dir = tmp_path / "chemkin"
    chem_dir.mkdir(parents=True)
    
    # Mock species dictionary
    dict_content = """
Ribose
// SMILES="OC[C@H]1OC(O)[C@H](O)[C@@H]1O"
1 0
multiplicity 1

Cysteine
// SMILES="N[C@@H](CS)C(=O)O"
1 0
multiplicity 1

Water
// SMILES="O"
1 0
multiplicity 1

Intermediate_1
// SMILES="C(=O)C(O)C"
1 0
multiplicity 1
"""
    (chem_dir / "species_dictionary.txt").write_text(dict_content)
    
    # Mock chem_annotated.inp
    inp_content = """
ELEMENTS
H C N O S
END

SPECIES
Ribose Cysteine Water Intermediate_1
END

REACTIONS
! Reaction family: Condensation
Ribose+Cysteine=Intermediate_1+Water 1.000E+00 0.00 0.00
END
"""
    (chem_dir / "chem_annotated.inp").write_text(inp_content)
    
    return tmp_path

def test_parse_species_dictionary(setup_mock_rmg_output):
    extractor = PathwayExtractor(setup_mock_rmg_output)
    extractor._parse_species_dictionary()
    
    assert "Ribose" in extractor.species_dict
    assert extractor.species_dict["Ribose"].smiles == "OC[C@H]1OC(O)[C@H](O)[C@@H]1O"
    assert "Intermediate_1" in extractor.species_dict
    assert extractor.species_dict["Intermediate_1"].smiles == "C(=O)C(O)C"

def test_parse_chem_inp(setup_mock_rmg_output):
    extractor = PathwayExtractor(setup_mock_rmg_output)
    extractor.run()
    
    assert len(extractor.elementary_steps) == 1
    step = extractor.elementary_steps[0]
    
    r_labels = [r.label for r in step.reactants]
    p_labels = [p.label for p in step.products]
    
    assert "Ribose" in r_labels
    assert "Cysteine" in r_labels
    assert "Intermediate_1" in p_labels
    assert "Water" in p_labels
    
def test_missing_files_raise_error(tmp_path):
    import shutil
    # Clean up the auto-used mock dir
    shutil.rmtree(tmp_path / "chemkin")
    
    empty_extractor = PathwayExtractor(tmp_path)
    (tmp_path / "chemkin").mkdir()
    
    with pytest.raises(FileNotFoundError):
        empty_extractor.run()
