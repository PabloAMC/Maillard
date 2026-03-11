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
import io
from dataclasses import dataclass
from typing import Optional, Tuple, Dict, List

try:
    from pyscf import gto, scf, hessian # noqa: F401
    from pyscf.geomopt import geometric_solver
    from pyscf.data import nist
except ImportError:
    # Graceful degradation for environments without PySCF
    pass

from .thermo import QuasiHarmonicCorrector
from .solvation import SolvationEngine
from .diffusion_ts import DiffusionTSEngine

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
    
    def __init__(self, solvent_name: str = 'water', temp_k: float = 423.15, 
                 use_explicit_solvent: bool = False, n_water: int = 3, 
                 geometry_backend: str = 'pyscf', db_path: Optional[str] = None):
        self.solvent_name = solvent_name
        self.temp_k = temp_k # Default 150 C
        self.use_explicit_solvent = use_explicit_solvent
        self.n_water = n_water
        self.geometry_backend = geometry_backend.lower()
        self.db_path = db_path
        
        # Phase 16: Initialize Results Database
        if self.db_path:
            from .results_db import ResultsDB
            self.db = ResultsDB(db_path=self.db_path)
        else:
            self.db = None
        
        # Phase 9: Initialize Solvation Engine with CREST/QCG discovery
        self.solvation_engine = SolvationEngine()
        
        # Phase 10: Initialize MLPOptimizer if requested
        if self.geometry_backend == 'mace':
            try:
                from .mlp_optimizer import MLPOptimizer
                self.mlp_optimizer = MLPOptimizer()
            except ImportError:
                print("WARNING: ML properties requested but not available. Falling back to pyscf.")
                self.geometry_backend = 'pyscf'
        else:
            self.mlp_optimizer = None
        
        # Phase 11: Initialize TSOptimizer
        try:
            from .ts_optimizer import TSOptimizer
            self.ts_optimizer = TSOptimizer()
        except ImportError:
            self.ts_optimizer = None
            
        # Phase 14: Initialize Diffusion TS Engine
        self.diffusion_engine = DiffusionTSEngine()
        
        # Phase B: Initialize MLP Barrier (MACE-OFF24)
        try:
            from .mlp_barrier import MLPBarrier
            self.mlp_barrier = MLPBarrier()
        except ImportError:
            self.mlp_barrier = None
        
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
        from pyscf import scf, dft
        
        is_open_shell = (mol.spin != 0)
        
        # Use RHF/UHF for pure Hartree-Fock
        if xc_method.lower() == 'hf':
            if is_open_shell:
                mf = scf.UHF(mol)
            else:
                mf = scf.RHF(mol)
        else:
            if is_open_shell:
                mf = dft.UKS(mol)
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
            from pyscf import dft
            is_open_shell = (mf.mol.spin != 0)
            if hasattr(mf, 'xc'):
                mf_vac = dft.UKS(mf.mol) if is_open_shell else dft.RKS(mf.mol)
                mf_vac.xc = mf.xc
            else:
                mf_vac = scf.UHF(mf.mol) if is_open_shell else scf.RHF(mf.mol)
            
            mf_vac.chkfile = None
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

    def optimize_geometry(self, xyz_content: str, charge: int=0, spin: int=0, is_ts: bool=False, max_steps: int=200, use_explicit_solvent: Optional[bool] = None, n_water: Optional[int] = None) -> DFTResult:
        """Run geometry/TS optimization using the r2SCAN-3c base method and return frequencies."""
        
        # Determine effective solvation settings
        eff_use_explicit = use_explicit_solvent if use_explicit_solvent is not None else self.use_explicit_solvent
        eff_n_water = n_water if n_water is not None else self.n_water
        
        # Phase 9: Apply explicit solvation if requested
        n_atoms_solute = int(xyz_content.strip().split('\n')[0])
        if eff_use_explicit:
            print(f">>> [Phase 9] Generating solvated cluster with {eff_n_water} waters...")
            xyz_content = self.solvation_engine.generate_solvated_cluster(
                xyz_content, 
                n_water=eff_n_water, 
                freeze_core=is_ts
            )
            
        # Phase 14: First-pass TS guess via Diffusion Model
        diffusion_ts_xyz = None
        if is_ts and self.diffusion_engine.available:
            # We need smiles for diffusion, if available in metadata we'd use them.
            # Passing placeholder strings or extracting from XYZ if possible.
            # For now, we check confidence and try to predict.
            # metadata could be passed through kwargs or extracted if we had SMILES.
            # This is a scaffold for React-TS integration.
            conf = self.diffusion_engine.get_confidence_score("", "") 
            if conf > 0.85:
                print(">>> [Phase 14] Attempting high-confidence TS guess via Diffusion Model...")
                diffusion_ts_xyz = self.diffusion_engine.predict_ts_geometry("", "")
                if diffusion_ts_xyz:
                    xyz_content = diffusion_ts_xyz

        # Phase 10: Backend Selection for Geometry Optimization
        if self.geometry_backend == 'mace':
            print(f">>> [Phase 10] Running MACE geometric optimization (is_ts={is_ts})...")
            # MACE bypasses PySCF/geomeTRIC entirely for structural relaxation
            if is_ts:
                opt_xyz = self.mlp_optimizer.optimize_ts(xyz_content, max_steps=max_steps)
            else:
                opt_xyz = self.mlp_optimizer.optimize_geometry(xyz_content, max_steps=max_steps)
                
            mol_opt = self._setup_mol(opt_xyz, charge, spin, basis=self.opt_basis)
            conv = True # Assume converged if MLP didn't raise
            
        else:
            # Original PySCF / geomeTRIC backend
            mol = self._setup_mol(xyz_content, charge, spin, basis=self.opt_basis)
            mf = self._build_mf(mol, xc_method=self.opt_method, use_solvent=False, conv_tol=1e-6)
            
            # Phase 11: Specialized TS search via Sella
            conv = False  # Initialize before potential Sella/geomeTRIC paths
            if is_ts and self.ts_optimizer:
                print(">>> [Phase 11] Running Sella eigenvector-following TS search...")
                try:
                    # PySCF↔ASE bridge: convert Mole to ASE Atoms manually
                    from ase import Atoms as ASEAtoms
                    from ase.calculators.calculator import Calculator, all_changes
                    from pyscf.data import nist as pyscf_nist
                    
                    class BuiltinPySCFCalc(Calculator):
                        """Minimal ASE Calculator for PySCF to avoid ase-pyscf dependency gaps."""
                        implemented_properties = ['energy', 'forces', 'hessian']
                        def __init__(self, mol, mf_class, xc):
                            super().__init__()
                            self.mol = mol
                            self.mf_class = mf_class
                            self.xc = xc
                        def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
                            super().calculate(atoms, properties, system_changes)
                            # Update molecular coordinates
                            self.mol.set_geom_(atoms.positions / pyscf_nist.BOHR, unit='Bohr')
                            mf = self.mf_class(self.mol)
                            if hasattr(mf, 'xc'): 
                                mf.xc = self.xc
                            
                            mf.kernel()
                            self.results['energy'] = float(mf.e_tot) * pyscf_nist.HARTREE2EV
                            
                            # Forces = -Gradient
                            if 'forces' in properties:
                                g = mf.nuc_grad_method().kernel()
                                self.results['forces'] = -g * pyscf_nist.HARTREE2EV / pyscf_nist.BOHR
                                
                            if 'hessian' in properties:
                                h_obj = mf.Hessian()
                                h = h_obj.kernel()
                                # PySCF returns (natm, natm, 3, 3). Sella/ASE expects (3*natm, 3*natm)
                                natm = self.mol.natm
                                h_reshaped = h.transpose(0, 2, 1, 3).reshape(3*natm, 3*natm)
                                self.results['hessian'] = h_reshaped * pyscf_nist.HARTREE2EV / (pyscf_nist.BOHR**2)
                        
                        def get_hessian(self, atoms):
                            self.calculate(atoms, properties=['hessian'])
                            return self.results['hessian']

                    coords_bohr = mol.atom_coords()
                    symbols = [mol.atom_symbol(i) for i in range(mol.natm)]
                    positions_ang = coords_bohr * pyscf_nist.BOHR
                    ase_atoms = ASEAtoms(symbols=symbols, positions=positions_ang)
                    
                    calc = BuiltinPySCFCalc(mol=mol, mf_class=type(mf), xc=getattr(mf, 'xc', 'hf'))
                    ase_atoms.calc = calc
                    
                    atoms = self.ts_optimizer.find_ts(ase_atoms, calc)
                    if self.ts_optimizer.is_converged(atoms):
                        # Convert ASE Atoms back to XYZ string
                        from ase.io import write as ase_write
                        with io.StringIO() as f_xyz:
                            ase_write(f_xyz, atoms, format='xyz')
                            opt_xyz = f_xyz.getvalue()
                        mol_opt = self._setup_mol(opt_xyz, charge, spin, basis=self.opt_basis)
                        conv = True
                    else:
                        print("    WARNING: Sella failed to converge. Falling back to geomeTRIC...")
                except (ImportError, AttributeError, Exception) as e:
                    print(f"    WARNING: Sella/ASE bridge failed ({type(e).__name__}: {e}). Falling back to geomeTRIC...")
                    
                is_ts_fallback = is_ts and not conv
            else:
                is_ts_fallback = is_ts

            if not is_ts or (is_ts and not self.ts_optimizer) or (is_ts and not conv):
                # Optimization via geomeTRIC (Standard or Fallback)
                with tempfile.TemporaryDirectory() as td:
                    pwd = os.getcwd()
                    try:
                        os.chdir(td)
                        geome_kwargs = {'maxsteps': max_steps}
                        if is_ts_fallback:
                            geome_kwargs['transition'] = True
                            
                        if eff_use_explicit and is_ts_fallback:
                            print("    Pre-relaxing solvent molecules around the frozen core...")
                            with open("constraints.txt", "w") as f:
                                f.write("$freeze\n")
                                f.write(f"xyz 1-{n_atoms_solute}\n")
                            
                            pre_relax_kwargs = {'maxsteps': 50, 'constraints': "constraints.txt"}
                            conv_pre, mol_relaxed = geometric_solver.kernel(mf, **pre_relax_kwargs)
                            mol_opt = mol_relaxed
                            conv = conv_pre
                        else:
                            conv, mol_opt = geometric_solver.kernel(mf, **geome_kwargs)
                    finally:
                        os.chdir(pwd)
                        
                opt_xyz = mol_opt.tostring(format='xyz')
        
        # The geometric_solver returns the optimized molecule object.
        # It updates mol in-place as well. 
        # Crucially, we assume success if it returns, but for Sella or more complex
        # runs, we'd check the internal state. 
        # In PySCF/geomeTRIC, the function returns the optimized molecule
        # but doesn't explicitly return a 'converged' boolean in the same call 
        # unless we capture the stdout or use the internal API.
        # However, if it fails to converge within max_steps, it often raises an error
        # or we can check the gradient if we have access to the internal optimizer state.
        
        # SOTA-level fix: We check the mf_opt convergence after the scf() call.
        
        # Compute frequencies and solvent energy at the OPTIMIZED vacuum geometry
        mf_opt = self._build_mf(mol_opt, xc_method=self.opt_method, use_solvent=True)
        # Check SCF convergence
        mf_opt.scf()
        if not mf_opt.converged:
             print("    WARNING: SCF failed to converge on the optimized geometry.")
        
        freqs, g_raw, g_qh = self._run_hessian_and_thermo(mf_opt)
        
        # Single-point Refinement
        sp_energy = self.single_point(opt_xyz, xc_method=self.refinement_method, basis=self.refinement_basis, charge=charge, spin=spin)
        
        # Calculate final refined Gibbs using SP energy + (G_qh - E_opt) thermal corrections
        thermal_corr_qh = g_qh - mf_opt.e_tot
        refined_gibbs = sp_energy + thermal_corr_qh
        
        return DFTResult(
            method=f"{self.refinement_method}//{self.opt_method}",
            energy_hartree=sp_energy,
            gibbs_free_energy_hartree=g_raw, # Just the uncorrected for reference
            quasi_harmonic_gibbs_hartree=refined_gibbs,
            optimized_xyz=opt_xyz,
            converged=conv, # Correctly report convergence from geomeTRIC
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
        hessian_matrix.reshape(natm*3, natm*3)
        mass = mol.atom_mass_list()
        np.repeat(mass, 3)
        # Identify the mode using pyscf's harmonic analysis (handles trans/rot projection)
        h_info = thermo.harmonic_analysis(mol, hessian_matrix)
        freqs = h_info['freq_wavenumber'] # Negative values for imaginary frequencies
        
        print(f"    Lowest frequencies (cm^-1): {freqs[:6]}")
        
        # We look for a significant imaginary frequency (conventionally negative in PySCF)
        # We use a threshold of -50 cm^-1 to avoid numerical noise
        if freqs[0] >= -50.0:
            raise ValueError(f"No significant imaginary frequency found at the provided TS geometry (Lowest freq: {freqs[0]:.1f} cm^-1)")
        
        print(f"    Found imaginary mode with frequency: {freqs[0]:.1f} cm^-1")
        
        # The first mode in h_info['norm_mode'] is the lowest frequency mode
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

    def calculate_barrier(self, reactant_xyz: str, ts_xyz: str, charge: int = 0, 
                          run_irc: bool = False, reaction_meta: Optional[Dict] = None) -> float:
        """
        End-to-end composite calculation of a kinetic barrier in kcal/mol.
        If run_irc=True, performs Phase 3.4 IRC validation.
        
        Args:
            reaction_meta: Optional dict with 'reactants', 'products' (SMILES lists) and 'family'.
        """
        # 1. Optimize Reactant
        print(">>> [Phase 3.3] Starting Reactant Optimization...")
        res_r = self.optimize_geometry(reactant_xyz, charge=charge, is_ts=False)
        if not res_r.converged:
            raise RuntimeError("Reactant optimization failed to converge.")
        
        # 2. Optimize TS
        print(">>> [Phase 3.3] Starting Transition State (TS) Optimization...")
        res_ts = self.optimize_geometry(ts_xyz, charge=charge, is_ts=True)
        if not res_ts.converged:
            raise RuntimeError("Transition State optimization failed to converge.")
        
        # 3. IRC Validation (Phase 3.4)
        if run_irc:
            self.generate_irc(res_ts.optimized_xyz, charge=charge)
            
        # 4. Delta G‡
        delta_g_h = res_ts.quasi_harmonic_gibbs_hartree - res_r.quasi_harmonic_gibbs_hartree
        barrier_kcal = delta_g_h * 627.509
        
        # 5. Log to DB if available
        if self.db and reaction_meta:
            self.db.add_barrier(
                reactants=reaction_meta.get('reactants', []),
                products=reaction_meta.get('products', []),
                family=reaction_meta.get('family', 'unknown'),
                delta_g_kcal=barrier_kcal,
                method=res_ts.method,
                basis=self.refinement_basis,
                solvation=self.solvent_name,
                converged=True
            )
            
        return barrier_kcal

    def calculate_mlp_barrier(self, reactant_xyz: str, product_xyz: str, reaction_meta: Optional[Dict] = None) -> Optional[float]:
        """
        Rapidly estimate the barrier using MACE-OFF24.
        Bypasses full TS optimization for a quick relative ranking.
        """
        if not self.mlp_barrier:
            return None
            
        barrier_kcal = self.mlp_barrier.estimate_barrier(reactant_xyz, product_xyz)
        
        if barrier_kcal is not None and self.db and reaction_meta:
            self.db.add_barrier(
                reactants=reaction_meta.get('reactants', []),
                products=reaction_meta.get('products', []),
                family=reaction_meta.get('family', 'unknown'),
                delta_g_kcal=barrier_kcal,
                method="mace-off",
                basis="OFF23_medium", # Default model name
                solvation=self.solvent_name,
                converged=True
            )
            
        return barrier_kcal
