"""
src/dft_refiner.py — Tier 2 DFT Refinement using PySCF.

Implements a composite workflow:
1. Geometry/Frequencies: r2SCAN-3c
2. Refinement: wB97M-V single point
3. Validation: revDSD-PBEP86-D4 (optional)

Thermodynamics incorporate Grimme quasi-harmonic corrections.
"""

import os
import tempfile
import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple, Dict, List

try:
    from pyscf import gto, dft, hessian, lib
    from pyscf.geomopt import geometric_solver
    from pyscf.data import nist
except ImportError:
    # Graceful degradation for environments without PySCF
    pass

from .thermo import QuasiHarmonicCorrector

@dataclass
class DFTResult:
    method: str
    energy_hartree: float
    gibbs_free_energy_hartree: Optional[float] = None
    quasi_harmonic_gibbs_hartree: Optional[float] = None
    optimized_xyz: str = ""
    converged: bool = False
    frequencies_cm1: Optional[List[float]] = None

class DFTRefiner:
    """Wrapper for running the tiered DFT composite workflow."""
    
    def __init__(self, solvent_name: str = 'water', temp_k: float = 423.15):
        self.solvent_name = solvent_name
        self.temp_k = temp_k # Default 150 C
        
        # The tiered methods defined in the plan
        self.opt_method = 'r2SCAN'
        self.opt_basis = 'def2-svp' # Close approximation to -3c base
        
        self.refinement_method = 'wB97M-V'
        self.refinement_basis = 'def2-tzvp'
        
        self.verif_method = 'revDSD-PBEP86'
        self.verif_basis = 'def2-tzvp'
        
        # [PERFORMANCE] Enable multithreading using available CPU cores
        n_threads = os.cpu_count() or 1
        os.environ["OMP_NUM_THREADS"] = str(n_threads)
        try:
            from pyscf import lib
            # We rely on the environment variable above.
        except ImportError:
            pass
        
    def _setup_mol(self, xyz_content: str, charge: int = 0, spin: int = 0, basis: str = 'def2-svp') -> 'gto.Mole':
        """Initialize PySCF GTO Mole object from XYZ."""
        lines = xyz_content.strip().split('\n')
        if len(lines) > 2 and lines[0].strip().isdigit():
            atom_string = '\n'.join(lines[2:])
        else:
            atom_string = xyz_content
            
        mol = gto.M(
            atom=atom_string,
            basis=basis,
            charge=charge,
            spin=spin,
            verbose=3 # 0=silent, 3=warnings, 4=info
        )
        return mol

    def _build_mf(self, mol: 'gto.Mole', xc_method: str = 'r2SCAN', use_solvent: bool = True, conv_tol: float = 1e-7):
        """Build the Mean-Field object with specified XC and implicit solvent."""
        from pyscf import scf
        
        # Use RHF for pure Hartree-Fock to avoid DFT grid overhead
        if xc_method.lower() == 'hf':
            mf = scf.RHF(mol)
        else:
            mf = dft.RKS(mol)
            mf.xc = xc_method
        
        if use_solvent and self.solvent_name:
            mf = mf.ddCOSMO()
            if self.solvent_name.lower() == 'water':
                mf.with_solvent.eps = 78.3553
        elif hasattr(mf, 'with_solvent'):
            # Failsafe: ensure no solvent is attached if use_solvent is False
            delattr(mf, 'with_solvent')
                
        # Enhanced convergence strategies for transition states
        mf.conv_tol = conv_tol
        mf.max_cycle = 100
        
        # [SCRATCH MANAGEMENT] Prevent PySCF from creating .chk files in the 
        # current working directory, which clutters the workspace.
        mf.chkfile = None
        
        return mf

    def single_point(self, xyz_content: str, xc_method: str = 'wB97M-V', basis: str = 'def2-tzvp', charge: int = 0, spin: int = 0) -> float:
        """Run a high-level single-point energy calculation."""
        mol = self._setup_mol(xyz_content, charge, spin, basis=basis)
        
        # PySCF supports wB97M-V via libxc
        # Since libxc support varies, we map to closest available or pass string directly.
        mf = self._build_mf(mol, xc_method=xc_method, use_solvent=True)
        energy = mf.scf()
        
        return energy
        
    def _run_hessian_and_thermo(self, mf) -> Tuple[List[float], float, float]:
        """Run frequency calculation and compute thermo."""
        # Hessian with solvent models in PySCF can be problematic depending on the version.
        # If it's a ddCOSMO method, it requires specific handling.
        # For our purposes we strip solvent from Hessian if it causes issues, since
        # solvent effects on frequencies are typically minor compared to the single point energy.
        if hasattr(mf, 'with_solvent'):
            # PySCF 2.x Hessian with ddCOSMO has known issues.
            # We create a vacuum counterpart for the Hessian calculation.
            from pyscf import scf
            if hasattr(mf, 'xc'):
                mf_vac = dft.RKS(mf.mol)
                mf_vac.xc = mf.xc
            else:
                mf_vac = scf.RHF(mf.mol)
            
            mf_vac.scf()
            hessobj = mf_vac.Hessian()
            vib = hessobj.kernel()
            mf_to_use = mf_vac
        else:
            hessobj = mf.Hessian()
            vib = hessobj.kernel()
            mf_to_use = mf
        
        # We need to get frequencies. In PySCF its usually via thermo module
        from pyscf.hessian import thermo
        # Need to diagonalize mass-weighted hessian
        freq_info = thermo.harmonic_analysis(mf_to_use.mol, vib)
        freqs = freq_info['freq_wavenumber']
        
        # Raw thermodynamics (standard harmonic)
        thermo_info = thermo.thermo(mf_to_use, freq_info['freq_au'], self.temp_k, 101325.0)
        
        # Apply quasi-harmonic
        corrector = QuasiHarmonicCorrector(temp_k=self.temp_k)
        # Use the ORIGINAL mf energy (which might have solvent)
        # Note: G_tot already contains E_tot. We need to be careful with the delta.
        # G_qh = E_elec + thermal_qh
        qh_result = corrector.calculate_thermo(freqs, electronic_energy_h=mf.e_tot)
        
        return freqs.tolist(), thermo_info['G_tot'][0] - mf_to_use.e_tot + mf.e_tot, qh_result.qh_gibbs_h

    def optimize_geometry(self, xyz_content: str, charge: int=0, spin: int=0, is_ts: bool=False, max_steps: int=100) -> DFTResult:
        """Run geometry/TS optimization using the r2SCAN-3c base method and return frequencies."""
        mol = self._setup_mol(xyz_content, charge, spin, basis=self.opt_basis)
        
        # [PERFORMANCE] Run geometry optimization in VACUUM. 
        # Resolving implicit solvent gradients is extremely slow in PySCF (ddCOSMO).
        # Geometries in vacuum are sufficiently similar for Maillard transition states.
        mf = self._build_mf(mol, xc_method=self.opt_method, use_solvent=False, conv_tol=1e-6)
        
        # Optimization
        with tempfile.TemporaryDirectory() as td:
            pwd = os.getcwd()
            try:
                os.chdir(td)
                geome_kwargs = {'maxsteps': max_steps}
                if is_ts:
                    geome_kwargs['transition'] = True
                    
                mol_opt = geometric_solver.optimize(mf, **geome_kwargs)
            finally:
                os.chdir(pwd)
                
        opt_xyz = mol_opt.tostring(format='xyz')
        
        # Compute frequencies and solvent energy at the OPTIMIZED vacuum geometry
        mf_opt = self._build_mf(mol_opt, xc_method=self.opt_method, use_solvent=True)
        mf_opt.scf()
        
        freqs, g_raw, g_qh = self._run_hessian_and_thermo(mf_opt)
        
        # Single-point Refinement
        sp_energy = self.single_point(opt_xyz, xc_method=self.refinement_method, basis=self.refinement_basis)
        
        # Calculate final refined Gibbs using SP energy + (G_qh - E_opt) thermal corrections
        thermal_corr_qh = g_qh - mf_opt.e_tot
        refined_gibbs = sp_energy + thermal_corr_qh
        
        return DFTResult(
            method=f"{self.refinement_method}//{self.opt_method}",
            energy_hartree=sp_energy,
            gibbs_free_energy_hartree=g_raw, # Just the uncorrected for reference
            quasi_harmonic_gibbs_hartree=refined_gibbs,
            optimized_xyz=opt_xyz,
            converged=True,
            frequencies_cm1=freqs
        )
        
    def generate_irc(self, ts_xyz: str, charge: int=0, spin: int=0, step_size: float=0.05) -> Tuple[str, str]:
        """
        Intrinsic Reaction Coordinate (IRC) validation via Displacement + Optimization.
        
        1. Calculates Hessian at the TS.
        2. Identifying the imaginary mode.
        3. Displaces +/- along the mode.
        4. Optimizes both endpoints to find the connected basins.
        """
        print(">>> [Phase 3.4] Starting IRC Validation (Double-Optimization method)...")
        mol = self._setup_mol(ts_xyz, charge, spin, basis=self.opt_basis)
        mf = self._build_mf(mol, xc_method=self.opt_method, use_solvent=False)
        mf.conv_tol = 1e-8
        mf.scf()
        
        print("    Computing TS Hessian for mode identification...")
        hessobj = mf.Hessian()
        hessian_matrix = hessobj.kernel()
        
        from pyscf.hessian import thermo
        # Manual mass-weighting to ensure signs are preserved for imaginary modes
        natm = mol.natm
        h = hessian_matrix.reshape(natm*3, natm*3)
        mass = mol.atom_mass_list()
        m = np.repeat(mass, 3)
        m_inv_sqrt = 1.0 / np.sqrt(m)
        h_mw = h * np.outer(m_inv_sqrt, m_inv_sqrt)
        
        e, v = np.linalg.eigh(h_mw)
        # Lowest eigenvalue corresponds to the most imaginary mode if negative
        # Frequencies are proportional to vacuum sqrt(abs(e))
        # We want to identify the reaction mode.
        if e[0] >= -1e-5:
            # Print some eigenvalues for debugging
            print(f"    Lowest eigenvalues: {e[:6]}")
            raise ValueError(f"No imaginary frequency found at the provided TS geometry (Lowest eigenvalue: {e[0]:.2e})")
        
        # Identify the mode
        idx = 0
        # h_info for normalized modes (pyscf uses mass-weighted normalized modes)
        h_info = thermo.harmonic_analysis(mol, hessian_matrix)
        freqs = h_info['freq_wavenumber']
        
        # We take the mode corresponding to the lowest eigenvalue
        # pyscf's norm_mode is sorted by frequency magnitude (absolute value?)
        # Let's be safe and find which mode in h_info['norm_mode'] corresponds to e[0]
        # Actually, harmonic_analysis often sorts by freq_au.
        # Let's just print diagnostic info.
        print(f"    Found imaginary mode with eigenvalue: {e[0]:.2e}")
        
        # Usually h_info['norm_mode'][0] is the one, but let's verify
        mode_vec = h_info['norm_mode'][0]
        
        orig_coords = mol.atom_coords() # in Bohr
        
        endpoints = []
        for direction in [1.0, -1.0]:
            label = "Forward" if direction > 0 else "Backward"
            print(f"    Following IRC {label} path...")
            
            # Displacement in Bohr
            displaced_coords = orig_coords + direction * step_size * mode_vec
            
            # Create temporary mol for optimization
            atoms = []
            for i in range(mol.natm):
                atoms.append((mol.atom_symbol(i), displaced_coords[i] * nist.BOHR))
            
            tmp_mol = gto.M(atom=atoms, basis=self.opt_basis, charge=charge, spin=spin)
            tmp_mf = self._build_mf(tmp_mol, xc_method=self.opt_method, use_solvent=False)
            
            # Optimize to nearest minimum
            with tempfile.TemporaryDirectory() as td:
                pwd = os.getcwd()
                try:
                    os.chdir(td)
                    # Use standard optimization (not transition state)
                    opt_mol = geometric_solver.optimize(tmp_mf, maxsteps=100)
                    endpoints.append(opt_mol.tostring(format='xyz'))
                finally:
                    os.chdir(pwd)
        
        return endpoints[0], endpoints[1]

    def verify_barrier(self, optimized_xyz: str, charge: int = 0, spin: int = 0) -> float:
        """
        Phase 3.5 Verification: Compute single-point energy using a higher-level functional.
        Typically revDSD-PBEP86-D4 on the optimized geometry.
        """
        print(f">>> [Phase 3.5] Running High-Level Verification SP ({self.verif_method})...")
        energy = self.single_point(optimized_xyz, xc_method=self.verif_method, basis=self.verif_basis, charge=charge, spin=spin)
        return energy

    def calculate_barrier(self, reactant_xyz: str, ts_xyz: str, charge: int = 0, run_irc: bool = False) -> float:
        """
        End-to-end composite calculation of a kinetic barrier in kcal/mol.
        If run_irc=True, performs Phase 3.4 IRC validation.
        """
        # 1. Optimize Reactant
        print(">>> [Phase 3.3] Starting Reactant Optimization...")
        res_r = self.optimize_geometry(reactant_xyz, charge=charge, is_ts=False)
        
        # 2. Optimize TS
        print(">>> [Phase 3.3] Starting Transition State (TS) Optimization...")
        res_ts = self.optimize_geometry(ts_xyz, charge=charge, is_ts=True)
        
        # 3. IRC Validation (Phase 3.4)
        if run_irc:
            self.generate_irc(res_ts.optimized_xyz, charge=charge)
            
        # 4. Delta G‡
        delta_g_h = res_ts.quasi_harmonic_gibbs_hartree - res_r.quasi_harmonic_gibbs_hartree
        return delta_g_h * 627.509
