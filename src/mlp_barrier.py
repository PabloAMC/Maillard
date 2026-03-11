"""
src/mlp_barrier.py

Phase B: Public ML Potential Integration.
Wraps MACE-OFF24 for rapid activation energy estimation.
"""

import io
from typing import Optional
from contextlib import redirect_stdout, redirect_stderr

try:
    from ase.io import read, write
    from mace.calculators import mace_off
    _MLP_AVAILABLE = True
except ImportError:
    read = write = mace_off = None
    _MLP_AVAILABLE = False


class MLPBarrier:
    """
    Wrapper for MACE-OFF24 interatomic potential.
    Provides energy calculation and barrier estimation.
    """

    def __init__(self, model: str = "medium", device: str = "cpu", default_dtype: str = "float64"):
        """
        Initialize the MACE interatomic potential.
        
        Args:
            model: "small", "medium", or "large". These pull MACE-OFF23 models.
            device: 'cpu' or 'cuda'/'mps'.
            default_dtype: Float precision ('float32' or 'float64').
        """
        if not _MLP_AVAILABLE:
            raise ImportError("MLPBarrier dependencies (mace-torch, ase) are not installed.")

        self.model = model
        self.device = device
        
        # Load the pre-trained foundation model
        assert mace_off is not None
        with io.StringIO() as buf, redirect_stdout(buf), redirect_stderr(buf):
            self.calc = mace_off(
                model=model,
                device=device,
                default_dtype=default_dtype
            )

    def get_energy(self, xyz_string: str) -> float:
        """Calculate potential energy in eV using MACE."""
        assert read is not None
        with io.StringIO(xyz_string.strip()) as f:
            atoms = read(f, format='xyz')
            if isinstance(atoms, list):
                atoms = atoms[-1]
        
        assert atoms is not None
        return atoms.get_potential_energy()

    def estimate_barrier(self, reactant_xyz: str, product_xyz: str) -> Optional[float]:
        """
        Estimates the reaction barrier (kcal/mol) as the difference 
        between reactant and product potential energies.
        
        Enforces strict stoichiometry check (atom count must match).
        """
        # Quick stoichiometry check
        r_lines = [line.strip() for line in reactant_xyz.strip().split('\n') if line.strip()]
        p_lines = [line.strip() for line in product_xyz.strip().split('\n') if line.strip()]
        
        if not r_lines or not p_lines:
            return None
            
        try:
            # First non-empty line of XYZ is the atom count
            r_atoms = int(r_lines[0])
            p_atoms = int(p_lines[0])
        except (ValueError, IndexError):
            return None

        if r_atoms != p_atoms:
            # print(f"[MLPBarrier] Stoichiometry mismatch: {r_atoms} != {p_atoms}. Skipping.")
            return None

        e_reac = self.get_energy(reactant_xyz)
        e_prod = self.get_energy(product_xyz)
        
        # Convert eV to kcal/mol
        ev_to_kcal_mol = 23.0605
        diff = (e_prod - e_reac) * ev_to_kcal_mol
        
        # For a simple ranking, we return max(0, diff) as a floor barrier
        return max(0.0, diff)

    def optimize_geometry(self, xyz_string: str, fmax: float = 0.01) -> str:
        """
        Refine a geometry using MACE-OFF24.
        """
        from ase.optimize import BFGS
        
        assert read is not None
        with io.StringIO(xyz_string) as f:
            atoms = read(f, format='xyz')
            if isinstance(atoms, list):
                atoms = atoms[-1]
        
        assert atoms is not None
        atoms.calc = self.calc
        
        with io.StringIO() as buf, redirect_stdout(buf):
            opt = BFGS(atoms, logfile=None)
            opt.run(fmax=fmax)
            
        assert write is not None
        with io.StringIO() as f:
            write(f, atoms, format='xyz')
            return f.getvalue()
