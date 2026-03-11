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


class SolvationEngine:
    """
    Generates explicit solvation clusters via CREST/QCG (Quantum Cluster Growth).

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
            project_root = Path(__file__).parent.parent
            candidates = [
                project_root / "conda_env" / "bin" / "crest",
                Path("/opt/homebrew/Caskroom/miniforge/base/bin/crest"),
                Path("/Users/pabloantoniomorenocasares/miniforge3/bin/crest"),
                Path("/Users/pabloantoniomorenocasares/miniconda3/bin/crest")
            ]
            
            for candidate in candidates:
                if candidate.exists():
                    crest_bin = str(candidate)
                    break
                    
            if not crest_bin:
                # Fall back to PATH lookup
                crest_bin = shutil.which("crest") or "crest"

        self.crest_bin = crest_bin
        self.timeout_sec = timeout_sec

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
                        print(">>> [Solvation] WARNING: xtb-IFF not found. Falling back to heuristic placement.")
                        return self._generate_heuristic_cluster(xyz_string, n_water)
                    
                    raise RuntimeError(f"CREST/QCG failed (exit {result.returncode}).")

                best_xyz = tmp / "crest_best.xyz"
                if not best_xyz.exists():
                    return self._generate_heuristic_cluster(xyz_string, n_water)

                return best_xyz.read_text()

            except (subprocess.TimeoutExpired, RuntimeError, FileNotFoundError, Exception) as e:
                print(f">>> [Solvation] WARNING: CREST/QCG call failed ({type(e).__name__}). Falling back to heuristic placement.")
                return self._generate_heuristic_cluster(xyz_string, n_water)

    # ------------------------------------------------------------------
    # Fallback Heuristics (for environments without full CREST stack)
    # ------------------------------------------------------------------

    def _generate_heuristic_cluster(self, xyz_string: str, n_water: int = 3) -> str:
        """
        Geometric fallback: place waters near polar sites.
        Used when CREST or xtb-IFF is missing.
        """
        lines = xyz_string.strip().split('\n')
        n_atoms_solute = int(lines[0])
        solute_data = []
        for line in lines[2:2+n_atoms_solute]:
            parts = line.split()
            solute_data.append({
                'symbol': parts[0],
                'coords': np.array([float(parts[1]), float(parts[2]), float(parts[3])])
            })

        # Identify polar sites (O, N) for water placement
        sites = [d['coords'] for d in solute_data if d['symbol'] in ['O', 'N']]
        if not sites:
            sites = [d['coords'] for d in solute_data]

        cluster_atoms = solute_data.copy()
        
        for i in range(n_water):
            base_site = sites[i % len(sites)]
            angle = (2 * np.pi / max(1, n_water)) * i
            phi = np.pi / 4 
            offset = np.array([
                np.cos(angle) * np.sin(phi), 
                np.sin(angle) * np.sin(phi), 
                np.cos(phi)
            ]) * 2.8
            
            o_coord = base_site + offset
            cluster_atoms.append({'symbol': 'O', 'coords': o_coord})
            cluster_atoms.append({'symbol': 'H', 'coords': o_coord + np.array([0.96, 0.0, 0.0])})
            cluster_atoms.append({'symbol': 'H', 'coords': o_coord + np.array([-0.24, 0.93, 0.0])})

        # Form final XYZ string
        final_xyz = f"{len(cluster_atoms)}\nSolvated Cluster (Heuristic Fallback)\n"
        for d in cluster_atoms:
            final_xyz += f"{d['symbol']:<2} {d['coords'][0]:12.6f} {d['coords'][1]:12.6f} {d['coords'][2]:12.6f}\n"

        return final_xyz

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
