import pytest
from pathlib import Path
from src.pathway_extractor import PathwayExtractor  # noqa: E402

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

def test_extractor_malformed_input(tmp_path):
    """Test robustness against malformed lines in chem_annotated.inp."""
    # Use a unique subdirectory to avoid collision with setup_mock_rmg_output fixture if it ran
    test_dir = tmp_path / "malformed_test"
    chem_dir = test_dir / "chemkin"
    chem_dir.mkdir(parents=True)
    
    # Valid dict
    (chem_dir / "species_dictionary.txt").write_text("A\n// SMILES=\"C\"\n1 0\n")
    
    # Malformed inp (missing '=' or empty species)
    inp_content = """
REACTIONS
A+B  1.0 0 0
A=A  1.0 0 0
END
"""
    (chem_dir / "chem_annotated.inp").write_text(inp_content)
    
    extractor = PathwayExtractor(test_dir)
    steps = extractor.run()
    
    # A+B should be ignored as it has no '=' (it has space)
    # A=A should be kept if A is resolved
    assert len(steps) == 1


# ---------------------------------------------------------------------------
# Audit remediation 2026-08-27 (extractor failure modes)
# ---------------------------------------------------------------------------

def _write_case(tmp_path, name, dict_content, inp_content):
    test_dir = tmp_path / name
    chem_dir = test_dir / "chemkin"
    chem_dir.mkdir(parents=True)
    (chem_dir / "species_dictionary.txt").write_text(dict_content)
    (chem_dir / "chem_annotated.inp").write_text(inp_content)
    return test_dir


def test_unresolved_species_fail_the_step_instead_of_being_dropped(tmp_path):
    test_dir = _write_case(
        tmp_path,
        "unresolved",
        'A\n// SMILES="C"\n\nC\n// SMILES="CO"\n',
        "REACTIONS\nA+Missing=C  1.0 0 0\nEND\n",
    )
    steps = PathwayExtractor(test_dir).run()

    assert len(steps) == 1
    step = steps[0]
    # Stoichiometry preserved: the unknown reactant is NOT silently removed.
    assert [r.label for r in step.reactants] == ["A", "Missing"]
    assert step.source_quality == "unresolved_species"
    assert step.unresolved_species == ["Missing"]


def test_equation_split_tolerates_equals_signs_inside_species_labels(tmp_path):
    test_dir = _write_case(
        tmp_path,
        "equals",
        'C=CC=O\n// SMILES="C=CC=O"\n\nCC=O\n// SMILES="CC=O"\n',
        "REACTIONS\nC=CC=O=CC=O  1.0 0 0\nEND\n",
    )
    steps = PathwayExtractor(test_dir).run()

    assert len(steps) == 1
    assert [r.label for r in steps[0].reactants] == ["C=CC=O"]
    assert [p.label for p in steps[0].products] == ["CC=O"]
    assert steps[0].unresolved_species == []


def test_reversible_arrows_are_flagged_and_can_emit_both_directions(tmp_path):
    dict_content = 'A\n// SMILES="C"\n\nB\n// SMILES="CO"\n'
    inp_content = "REACTIONS\nA<=>B  1.0 0 0\nEND\n"

    flagged = PathwayExtractor(_write_case(tmp_path, "rev_flag", dict_content, inp_content)).run()
    assert len(flagged) == 1
    assert flagged[0].reversible is True
    assert flagged[0].direction == "forward"

    both = PathwayExtractor(
        _write_case(tmp_path, "rev_both", dict_content, inp_content),
        emit_reverse_steps=True,
    ).run()
    assert len(both) == 2
    assert both[1].direction == "reverse"
    assert [r.label for r in both[1].reactants] == ["B"]


def test_irreversible_arrow_is_not_marked_reversible(tmp_path):
    test_dir = _write_case(
        tmp_path,
        "irrev",
        'A\n// SMILES="C"\n\nB\n// SMILES="CO"\n',
        "REACTIONS\nA=>B  1.0 0 0\nEND\n",
    )
    steps = PathwayExtractor(test_dir).run()
    assert len(steps) == 1
    assert steps[0].reversible is False


def test_real_rmg_adjacency_lists_are_parsed(tmp_path):
    from rdkit import Chem

    adjacency = "\n".join([
        "water",
        "1 O u0 p2 c0 {2,S} {3,S}",
        "2 H u0 p0 c0 {1,S}",
        "3 H u0 p0 c0 {1,S}",
        "",
        "acetaldehyde",
        "multiplicity 1",
        "1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}",
        "2 C u0 p0 c0 {1,S} {3,D} {7,S}",
        "3 O u0 p2 c0 {2,D}",
        "4 H u0 p0 c0 {1,S}",
        "5 H u0 p0 c0 {1,S}",
        "6 H u0 p0 c0 {1,S}",
        "7 H u0 p0 c0 {2,S}",
    ])
    test_dir = _write_case(tmp_path, "adjacency", adjacency, "REACTIONS\nEND\n")

    extractor = PathwayExtractor(test_dir)
    extractor.run()

    assert extractor.species_dict["water"].smiles == Chem.CanonSmiles("O")
    assert extractor.species_dict["acetaldehyde"].smiles == Chem.CanonSmiles("CC=O")


def test_unsupported_species_block_raises_instead_of_returning_empty(tmp_path):
    from src.pathway_extractor import UnsupportedSpeciesFormatError

    test_dir = _write_case(
        tmp_path,
        "unsupported",
        "MysterySpecies\nsome free-form text that is neither SMILES nor adjacency\n",
        "REACTIONS\nEND\n",
    )
    with pytest.raises(UnsupportedSpeciesFormatError):
        PathwayExtractor(test_dir).run()
