"""
tests/test_solvation.py

Phase 9 — Verification suite for SolvationEngine.

Tests the CREST/QCG wrapper, including:
- Correct atom count after solvation.
- Preservation of solute coordinates (< 0.05 Å tolerance) when freeze_core=True.
- Correct structure of the generated .xcontrol file.
- Graceful error handling when CREST is unavailable or times out.
"""

import pytest
import numpy as np
from unittest.mock import patch, MagicMock
from pathlib import Path
from src.solvation import SolvationEngine


# ---------------------------------------------------------------------------
# Test Fixtures
# ---------------------------------------------------------------------------

FORMALDEHYDE_XYZ = """\
4
Formaldehyde
O       0.000000    0.000000    1.200000
C       0.000000    0.000000    0.000000
H       0.000000    0.940000   -0.580000
H       0.000000   -0.940000   -0.580000
"""

N_SOLUTE = 4  # Formaldehyde atom count


def make_solvated_xyz(n_water: int, solute_xyz: str = FORMALDEHYDE_XYZ) -> str:
    """
    Build a mock CREST output: solute atoms in their original positions, then
    n_water * 3 water atoms appended. Simulates a perfect freeze_core output.
    """
    solute_lines = solute_xyz.strip().splitlines()
    n_solute = int(solute_lines[0])
    coord_lines = solute_lines[2:2 + n_solute]

    water_lines = []
    for i in range(n_water):
        oy = 5.0 + i * 3.0  # place waters far from solute
        water_lines += [
            f"O    0.0  {oy:.3f}  0.0",
            f"H    0.96 {oy:.3f}  0.0",
            f"H   -0.24 {oy:.3f}  0.93",
        ]

    total = n_solute + n_water * 3
    lines = [str(total), "Solvated Cluster (mock CREST output)"] + coord_lines + water_lines
    return "\n".join(lines) + "\n"


def extract_coords(xyz_str: str, n_atoms: int) -> np.ndarray:
    """Extract the first n_atoms coordinate rows from an XYZ string."""
    body = xyz_str.strip().splitlines()[2:2 + n_atoms]
    return np.array([[float(v) for v in row.split()[1:4]] for row in body])


# ---------------------------------------------------------------------------
# Unit tests (mock subprocess — CI-safe; do not require crest binary)
# ---------------------------------------------------------------------------

class TestXControlGeneration:
    """Verify the .xcontrol constraint string is correctly constructed."""

    def test_xcontrol_contains_correct_atom_range(self):
        content = SolvationEngine._build_xcontrol(n_solute_atoms=6)
        assert "$fix" in content
        assert "atoms: 1-6" in content
        assert "$end" in content

    def test_xcontrol_single_atom(self):
        content = SolvationEngine._build_xcontrol(n_solute_atoms=1)
        assert "atoms: 1-1" in content

    def test_n_atoms_parsed_correctly(self):
        n = SolvationEngine._parse_n_atoms(FORMALDEHYDE_XYZ)
        assert n == N_SOLUTE


class TestSolvatedClusterAtomCount:
    """After solvation the cluster must contain exactly n_solute + n_water*3 atoms."""

    @pytest.mark.parametrize("n_water", [1, 3, 6])
    def test_correct_atom_count(self, n_water: int):
        mock_xyz = make_solvated_xyz(n_water)
        engine = SolvationEngine()

        with patch("subprocess.run") as mock_run, \
             patch("pathlib.Path.exists", return_value=True), \
             patch("pathlib.Path.write_text"), \
             patch("pathlib.Path.read_text", return_value=mock_xyz):

            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            result = engine.generate_solvated_cluster(FORMALDEHYDE_XYZ, n_water=n_water)

        actual_n = int(result.strip().splitlines()[0])
        assert actual_n == N_SOLUTE + n_water * 3, (
            f"Expected {N_SOLUTE + n_water * 3} atoms, got {actual_n}"
        )


class TestSolutePreservation:
    """
    CRITICAL (Phase 9 gate): solute atoms must not shift by more than 0.05 Å
    when freeze_core=True. This guards TS geometries from being corrupted.
    """

    def test_solute_coords_preserved_within_tolerance(self):
        """With freeze_core=True the first N_SOLUTE atom positions must be unchanged."""
        n_water = 3
        mock_xyz = make_solvated_xyz(n_water)  # solute coords identical to input
        engine = SolvationEngine()

        with patch("subprocess.run") as mock_run, \
             patch("pathlib.Path.exists", return_value=True), \
             patch("pathlib.Path.write_text"), \
             patch("pathlib.Path.read_text", return_value=mock_xyz):

            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            result = engine.generate_solvated_cluster(
                FORMALDEHYDE_XYZ, n_water=n_water, freeze_core=True
            )

        original_coords = extract_coords(FORMALDEHYDE_XYZ, N_SOLUTE)
        output_coords = extract_coords(result, N_SOLUTE)
        displacement = np.linalg.norm(output_coords - original_coords, axis=1)

        assert np.all(displacement < 0.05), (
            f"Solute atoms shifted by {displacement.max():.4f} Å — "
            "freeze_core must prevent any TS geometry distortion."
        )

    def test_freeze_core_passes_xcontrol_flag(self):
        """Confirm -cinp is passed to crest when freeze_core=True."""
        n_water = 2
        mock_xyz = make_solvated_xyz(n_water)
        engine = SolvationEngine()

        with patch("subprocess.run") as mock_run, \
             patch("pathlib.Path.exists", return_value=True), \
             patch("pathlib.Path.write_text"), \
             patch("pathlib.Path.read_text", return_value=mock_xyz):

            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            engine.generate_solvated_cluster(
                FORMALDEHYDE_XYZ, n_water=n_water, freeze_core=True
            )
            call_args = mock_run.call_args[0][0]

        assert "-cinp" in call_args, "-cinp flag must be passed to crest when freeze_core=True"

    def test_no_xcontrol_when_freeze_core_false(self):
        """Confirm -cinp is NOT passed when freeze_core=False."""
        n_water = 2
        mock_xyz = make_solvated_xyz(n_water)
        engine = SolvationEngine()

        with patch("subprocess.run") as mock_run, \
             patch("pathlib.Path.exists", return_value=True), \
             patch("pathlib.Path.write_text"), \
             patch("pathlib.Path.read_text", return_value=mock_xyz):

            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            engine.generate_solvated_cluster(
                FORMALDEHYDE_XYZ, n_water=n_water, freeze_core=False
            )
            call_args = mock_run.call_args[0][0]

        assert "-cinp" not in call_args, "-cinp must NOT be passed when freeze_core=False"


class TestErrorHandling:
    """Engine must gracefully fall back to heuristics when CREST fails or times out."""

    def test_fallback_on_nonzero_exit(self):
        engine = SolvationEngine()
        with patch("subprocess.run") as mock_run, \
             patch("pathlib.Path.write_text"):
            mock_run.return_value = MagicMock(returncode=1, stdout="", stderr="crest died")
            # Should NOT raise, but return heuristic cluster
            result = engine.generate_solvated_cluster(FORMALDEHYDE_XYZ, n_water=1)
            assert "Heuristic Fallback" in result

    def test_fallback_on_slow_crest(self):
        import subprocess as sp
        engine = SolvationEngine(timeout_sec=1)
        with patch("subprocess.run", side_effect=sp.TimeoutExpired("crest", 1)), \
             patch("pathlib.Path.write_text"):
            # Should NOT raise, but return heuristic cluster
            result = engine.generate_solvated_cluster(FORMALDEHYDE_XYZ, n_water=1)
            assert "Heuristic Fallback" in result

    def test_fallback_when_output_missing(self):
        engine = SolvationEngine()
        with patch("subprocess.run") as mock_run, \
             patch("pathlib.Path.write_text"), \
             patch("pathlib.Path.exists", return_value=False):
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            # Should NOT raise, but return heuristic cluster
            result = engine.generate_solvated_cluster(FORMALDEHYDE_XYZ, n_water=1)
            assert "Heuristic Fallback" in result
