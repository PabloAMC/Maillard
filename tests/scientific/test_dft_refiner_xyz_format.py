from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.dft_refiner import DFTRefiner, _extract_atom_string, _molecule_to_xyz  # noqa: E402


def test_extract_atom_string_handles_xyz_without_comment_line():
    xyz_content = "2\nH 0.0 0.0 0.0\nO 0.0 0.0 1.0\n"

    atom_string = _extract_atom_string(xyz_content)

    assert atom_string.splitlines() == ["H 0.0 0.0 0.0", "O 0.0 0.0 1.0"]


def test_molecule_to_xyz_round_trips_through_setup_mol():
    refiner = DFTRefiner(solvent_name="water", temp_k=423.15)
    mol = refiner._setup_mol("2\ncomment\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n", basis="sto-3g")

    xyz_content = _molecule_to_xyz(mol, comment="round trip")
    reparsed = refiner._setup_mol(xyz_content, basis="sto-3g")

    assert reparsed.natm == 2
    assert reparsed.atom_symbol(0) == "H"
    assert reparsed.atom_symbol(1) == "H"