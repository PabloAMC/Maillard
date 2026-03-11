"""
src/solvation.py

Phase 9: Explicit Solvation Automation via CREST/QCG.

Uses the CREST Quantum Cluster Growth (QCG) algorithm to generate an explicit
solvation shell of water molecules around a given solute XYZ geometry. When
freeze_core=True (required for transition states), a .xcontrol file is generated
that freezes all solute Cartesian coordinates so only the water molecules are
sampled — preventing CREST from collapsing the TS into a local minimum.
"""

import os
import subprocess
import tempfile
import shutil
import numpy as np
from pathlib import Path
from typing import Optional
from typing import List


class SolvationEngine:
    """
    Handles generation of explicit water clusters natively using CREST.
    The 'freeze_core' option writes a .xcontrol constraint file that pins all
    solute atom coordinates while CREST samples conformers for the water shell.
    This is mandatory for transition state geometries.
    """

    def __init__(self, crest_bin: Optional[str] = None, timeout_sec: int = 300):
        """
        Args:
            crest_bin: Path to the crest binary. Auto-detects from conda_env/ if None.
            timeout_sec: Maximum seconds to wait for crest subprocess (default 5 min).
        """
        if crest_bin is None:
            # Auto-resolve relative to this file's parent (project root)
            project_root: Path = Path(__file__).parent.parent
            home: Path = Path.home()
            conda_prefix = os.environ.get("CONDA_PREFIX")
            
            candidates: List[Path] = []
            if conda_prefix:
                candidates.append(Path(conda_prefix) / "bin" / "crest")
                
            candidates.extend([
                project_root / "conda_env" / "bin" / "crest",
                home / "miniforge3" / "bin" / "crest",
                home / "miniconda3" / "bin" / "crest",
                Path("/usr/local/bin/crest"),
                Path("/opt/homebrew/bin/crest")
            ])
            
            for candidate in candidates:
                if candidate.exists():
                    crest_bin = str(candidate)
                    break
                    
            if not crest_bin:
                # Fall back to PATH lookup
                crest_bin = shutil.which("crest") or "crest"

        self.crest_bin: str = crest_bin
        self.timeout_sec: int = timeout_sec

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def generate_solvated_cluster(
        self,
        xyz_string: str,
        n_water: int = 3,
        freeze_core: bool = True,
    ) -> str:
        """
        Run CREST/QCG to grow a solvation shell of n_water molecules around
        the input geometry and return the best cluster as an XYZ string.

        Args:
            xyz_string:  Input geometry in XYZ format (header + coord lines).
            n_water:     Number of explicit water molecules to add.
            freeze_core: If True, freeze solute Cartesian coordinates so the
                         TS geometry is preserved during conformational sampling.

        Returns:
            The best solvated cluster as an XYZ-format string.

        Raises:
            RuntimeError: If CREST terminates with a non-zero exit code.
            FileNotFoundError: If CREST does not produce an output file.
            TimeoutError: If CREST exceeds self.timeout_sec.
        """
        n_solute_atoms = self._parse_n_atoms(xyz_string)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)

            # 1. Write solute geometry
            input_xyz = tmp / "input.xyz"
            input_xyz.write_text(xyz_string)

            # 2. Write solvent molecule definition (required by CREST -qcg)
            if n_water > 0:
                solvent_xyz = tmp / "water.xyz"
                solvent_xyz.write_text("3\nwater\nO 0.0 0.0 0.0\nH 0.0 0.0 0.9584\nH 0.0 0.9240 -0.2529\n")

            # 3. Write .xcontrol constraint file if freezing solute core
            xcontrol_path: Optional[Path] = None
            if freeze_core:
                xcontrol_content = self._build_xcontrol(n_solute_atoms)
                xcontrol_path = tmp / ".xcontrol"
                xcontrol_path.write_text(xcontrol_content)

            # 4. Build the CREST command
            #    -qcg activates Quantum Cluster Growth mode
            #    -nsolv controls number of solvent molecules
            #    -cinp points to our constraint file
            cmd = [
                self.crest_bin,
                str(input_xyz),
                "-qcg", "water",
                "-nsolv", str(n_water),
            ]
            if xcontrol_path is not None:
                cmd += ["-cinp", str(xcontrol_path)]

            # 5. Run CREST
            try:
                result = subprocess.run(
                    cmd,
                    cwd=str(tmp),
                    capture_output=True,
                    text=True,
                    timeout=self.timeout_sec,
                )
                
                if result.returncode != 0:
                    # Specific check for missing xtb-IFF
                    if "xtb-IFF" in result.stderr or "xtb-IFF" in result.stdout:
                        raise RuntimeError(f"CREST/QCG failed (exit {result.returncode}): xtb-IFF not found. Please ensure xtb-IFF is installed and in your PATH.")
                    
                    raise RuntimeError(f"CREST/QCG failed (exit {result.returncode}): {result.stderr}")

                best_xyz = tmp / "crest_best.xyz"
                if not best_xyz.exists():
                    raise FileNotFoundError("CREST completed but 'crest_best.xyz' was not generated.")

                return best_xyz.read_text()

            except (subprocess.TimeoutExpired, RuntimeError, FileNotFoundError) as e:
                raise RuntimeError(f"Solvation error: {str(e)}") from e


    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_n_atoms(xyz_string: str) -> int:
        """Return the atom count from the first line of an XYZ string."""
        first_line = xyz_string.strip().splitlines()[0].strip()
        return int(first_line)

    @staticmethod
    def _build_xcontrol(n_solute_atoms: int) -> str:
        """
        Generate a CREST .xcontrol constraint block that fixes the Cartesian
        coordinates of atoms 1 through n_solute_atoms (the solute core).

        Syntax: $fix  atoms: 1-N  end

        This prevents CREST from distorting the TS geometry while it samples
        the water-shell conformational space.
        """
        return (
            "$fix\n"
            f"  atoms: 1-{n_solute_atoms}\n"
            "$end\n"
        )
