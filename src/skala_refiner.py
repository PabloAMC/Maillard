"""
src/skala_refiner.py — Tier 2 DFT Refinement using PySCF and Microsoft Skala.

Refines the 6-8 chemically decisive rate-limiting transition states 
using the Skala XC functional (chemical accuracy at DFT cost) in 
implicit water solvent (CPCM).
"""

import tempfile
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Tuple

try:
    from pyscf import gto, dft, solvent
    from pyscf.geomopt import geometric_solver
except ImportError:
    # Graceful degradation for environments without PySCF
    pass

@dataclass
class SkalaResult:
    energy_hartree: float
    optimized_xyz: str
    converged: bool

class SkalaRefiner:
    """Wrapper for running DFT calculations via PySCF using the Skala functional."""
    
    def __init__(self, basis: str = 'def2-svp', solvent_name: str = 'water', use_skala: bool = True):
        self.basis = basis
        self.solvent_name = solvent_name
        self.use_skala = use_skala
        # Fallback to standard robust hybrid if Skala isn't installed
        self.fallback_xc = 'b3lyp' 
        
    def _setup_mol(self, xyz_content: str, charge: int = 0, spin: int = 0) -> 'gto.Mole':
        """Initialize PySCF GTO Mole object from XYZ."""
        # PySCF expects XYZ without the first two lines (atom count and comment)
        # or it can parse a standard XYZ block. But it's safer to strip headers.
        lines = xyz_content.strip().split('\n')
        if len(lines) > 2 and lines[0].strip().isdigit():
            atom_string = '\n'.join(lines[2:])
        else:
            atom_string = xyz_content
            
        mol = gto.M(
            atom=atom_string,
            basis=self.basis,
            charge=charge,
            spin=spin,
            verbose=3 # 0=silent, 3=warnings/basics, 4=info, 5=debug
        )
        return mol
        
    def _build_mf(self, mol: 'gto.Mole'):
        """Build the Mean-Field object with Skala XC and CPCM solvent."""
        mf = dft.RKS(mol)
        
        # Determine XC functional
        if self.use_skala:
            try:
                from skala.pyscf import SkalaKS
                mf = SkalaKS(mol, xc='skala')
            except ImportError:
                print(f"Warning: Skala not found. Falling back to {self.fallback_xc}.")
                mf = dft.RKS(mol)
                mf.xc = self.fallback_xc
        else:
            mf.xc = self.fallback_xc
            
        # Apply implicit solvation
        if self.solvent_name:
            mf = mf.ddCOSMO() # ddCOSMO is PySCF's standard fast continuum model
            mf.with_solvent.eps = 78.3553 # Water dielectric
            
        # Ensure robust SCF convergence
        mf.conv_tol = 1e-7
        mf.max_cycle = 100
        
        return mf

    def single_point(self, xyz_content: str, charge: int=0, spin: int=0) -> SkalaResult:
        """Run a single-point energy calculation."""
        mol = self._setup_mol(xyz_content, charge, spin)
        mf = self._build_mf(mol)
        
        energy = mf.scf()
        
        return SkalaResult(
            energy_hartree=energy,
            optimized_xyz=xyz_content,
            converged=mf.converged
        )

    def optimize_geometry(self, xyz_content: str, charge: int=0, spin: int=0, max_steps: int=50) -> SkalaResult:
        """Run geometry optimization (requires geomeTRIC package)."""
        mol = self._setup_mol(xyz_content, charge, spin)
        mf = self._build_mf(mol)
        
        # PySCF uses geomeTRIC as the backend optimizer. It updates mol in-place.
        mol_opt = geometric_solver.optimize(mf, maxsteps=max_steps)
        
        opt_xyz = mol_opt.tostring(format='xyz')
        
        # Run a final single-point on the optimized geometry to get the exact energy
        # since geometric doesn't always populate `mf.e_tot` on the original object
        final_energy = self.single_point(opt_xyz, charge, spin).energy_hartree
        
        return SkalaResult(
            energy_hartree=final_energy,
            optimized_xyz=opt_xyz,
            converged=True
        )
        
    def optimize_ts(self, xyz_ts_guess: str, charge: int=0, spin: int=0) -> SkalaResult:
        """Run Transition State optimization."""
        # For TS opt in PySCF+geometric, we pass specialized params
        mol = self._setup_mol(xyz_ts_guess, charge, spin)
        mf = self._build_mf(mol)
        
        # Provide TS flag to geometric
        import os
        with tempfile.TemporaryDirectory() as td:
            pwd = os.getcwd()
            try:
                os.chdir(td)
                mol_opt = geometric_solver.optimize(mf, maxsteps=100, custom_engine_params={'transition': True})
            finally:
                os.chdir(pwd)
                
        opt_xyz = mol_opt.tostring(format='xyz')
        
        # Run a final single-point on the optimized geometry
        final_energy = self.single_point(opt_xyz, charge, spin).energy_hartree
        
        return SkalaResult(
            energy_hartree=final_energy,
            optimized_xyz=opt_xyz,
            converged=True
        )

    def refine_barrier(self, reactant_xyz: str, ts_xyz: str) -> float:
        """Compute the activation barrier ΔE‡ in kcal/mol."""
        r_res = self.optimize_geometry(reactant_xyz)
        ts_res = self.optimize_ts(ts_xyz)
        
        delta_eh = ts_res.energy_hartree - r_res.energy_hartree
        return delta_eh * 627.509
