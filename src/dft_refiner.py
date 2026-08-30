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
import json
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Any, Tuple
from pathlib import Path
import time

from src.dft_geometry_preflight import repair_steric_clash
from src.xyz_common import is_xyz_coordinate_line
from src.xyz_fragment_utils import split_xyz_into_fragments, assign_fragment_charge_spin

class OptimizationTimeoutError(Exception):
    pass


class BadTSGuessRejected(RuntimeError):
    """Raised when the TS guess fails all refinement attempts and is invalid."""
    pass


# Phase-2 hardening: imaginary frequency below this magnitude (cm^-1) is treated
# as numerical noise from residual translation/rotation, not a real reaction
# coordinate mode. A "true" TS must have exactly one imaginary mode whose
# magnitude exceeds this floor.
TS_IMAGINARY_FREQ_FLOOR_CM1 = 50.0


def _classify_ts_frequencies(
    frequencies_cm1: Optional[List[float]],
    *,
    floor_cm1: float = TS_IMAGINARY_FREQ_FLOOR_CM1,
) -> Dict[str, Any]:
    """Classify a TS frequency spectrum into a saddle-point verdict.

    Returns a dict with:
        n_imaginary: count of f < 0
        n_significant_imaginary: count of f < -floor_cm1
        imaginary_freqs_cm1: full list of negative frequencies (sorted ascending)
        most_negative_cm1: smallest (most negative) frequency or None
        is_true_saddle: True iff exactly one significant imaginary mode
        verdict: one of {"true_ts", "minimum", "high_order_saddle",
                          "weak_imaginary", "no_frequencies"}
        floor_cm1: threshold used
    """
    if not frequencies_cm1:
        return {
            "n_imaginary": 0,
            "n_significant_imaginary": 0,
            "imaginary_freqs_cm1": [],
            "most_negative_cm1": None,
            "is_true_saddle": False,
            "verdict": "no_frequencies",
            "floor_cm1": float(floor_cm1),
        }
    imag = sorted([float(f) for f in frequencies_cm1 if f < 0])
    sig_imag = [f for f in imag if f < -floor_cm1]
    n_imag = len(imag)
    n_sig = len(sig_imag)
    most_neg = imag[0] if imag else None
    if n_sig == 1:
        verdict = "true_ts"
    elif n_sig == 0 and n_imag == 0:
        verdict = "minimum"
    elif n_sig == 0 and n_imag > 0:
        verdict = "weak_imaginary"
    else:  # n_sig >= 2
        verdict = "high_order_saddle"
    return {
        "n_imaginary": n_imag,
        "n_significant_imaginary": n_sig,
        "imaginary_freqs_cm1": imag,
        "most_negative_cm1": most_neg,
        "is_true_saddle": (n_sig == 1),
        "verdict": verdict,
        "floor_cm1": float(floor_cm1),
    }
from src.logger import get_logger

logger = get_logger(__name__)

try:
    from pyscf import gto, scf, hessian # noqa: F401
    from pyscf.geomopt import geometric_solver
    from pyscf.data import nist
except ImportError:
    # Graceful degradation for environments without PySCF
    gto = scf = hessian = geometric_solver = nist = None

from .thermo import QuasiHarmonicCorrector
from .solvation import SolvationEngine
from .diffusion_ts import DiffusionTSEngine
from . import xtb_backend as _xtb


def _normalize_reaction_family(value: Any) -> str:
    if isinstance(value, str):
        return value
    if isinstance(value, (list, tuple, set)):
        items = [str(item) for item in value if str(item).strip()]
        if not items:
            return "unknown"
        if len(items) == 1:
            return items[0]
        return ", ".join(items)
    if value is None:
        return "unknown"
    return str(value)


def _normalize_reaction_side(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        return [str(item) for item in value]
    return [str(value)]


# Producers of `reaction_meta` do not agree on key names, and the two that
# drive the DFT campaigns (scripts/run_computational_gap_dft.py and
# src/elementary_step_runner.py) supply neither 'reactants' nor 'products'.
# `.get('reactants', [])` therefore silently filed every DFT barrier under the
# same empty-identity reaction ([] -> []), which is unaddressable by any lookup.
# Accept the known aliases, and refuse to persist a row that has no identity at
# all rather than writing another one into that bucket.
_REACTION_SIDE_ALIASES: Dict[str, Tuple[str, ...]] = {
    "reactants": ("reactants", "reactant_smiles", "reactants_smiles"),
    "products": ("products", "product_smiles", "products_smiles"),
}


def _resolve_reaction_side(reaction_meta: Optional[Dict], side: str) -> List[str]:
    meta = reaction_meta or {}
    for key in _REACTION_SIDE_ALIASES[side]:
        value = meta.get(key)
        if value:
            return _normalize_reaction_side(value)
    return []

@dataclass
class DFTResult:
    method: str
    energy_hartree: float
    gibbs_free_energy_hartree: Optional[float] = None
    quasi_harmonic_gibbs_hartree: Optional[float] = None
    optimized_xyz: str = ""
    converged: bool = False
    frequencies_cm1: Optional[List[float]] = None

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to a JSON-safe dictionary."""
        return {
            "method": self.method,
            "energy_hartree": self.energy_hartree,
            "gibbs_free_energy_hartree": self.gibbs_free_energy_hartree,
            "quasi_harmonic_gibbs_hartree": self.quasi_harmonic_gibbs_hartree,
            "optimized_xyz": self.optimized_xyz,
            "converged": self.converged,
            "frequencies_cm1": self.frequencies_cm1,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "DFTResult":
        """Reconstruct from a dictionary (e.g. loaded from JSON checkpoint)."""
        return cls(
            method=d["method"],
            energy_hartree=d["energy_hartree"],
            gibbs_free_energy_hartree=d.get("gibbs_free_energy_hartree"),
            quasi_harmonic_gibbs_hartree=d.get("quasi_harmonic_gibbs_hartree"),
            optimized_xyz=d.get("optimized_xyz", ""),
            converged=d.get("converged", False),
            frequencies_cm1=d.get("frequencies_cm1"),
        )


def _extract_atom_string(xyz_content: str) -> str:
    lines = [line.rstrip() for line in xyz_content.strip().splitlines() if line.strip()]
    if not lines:
        return ""
    if lines[0].strip().isdigit():
        atom_count = int(lines[0].strip())
        remaining = lines[1:]
        if remaining and not is_xyz_coordinate_line(remaining[0]):
            remaining = remaining[1:]
        atom_lines = remaining[:atom_count]
        return "\n".join(atom_lines)
    return "\n".join(lines)


def _molecule_to_xyz(mol: Any, comment: str = "") -> str:
    coords_angstrom = mol.atom_coords() * nist.BOHR
    lines = [str(mol.natm), comment]
    for atom_index in range(mol.natm):
        symbol = mol.atom_symbol(atom_index)
        x_coord, y_coord, z_coord = coords_angstrom[atom_index]
        lines.append(f"{symbol:<2} {x_coord: .12f} {y_coord: .12f} {z_coord: .12f}")
    return "\n".join(lines) + "\n"

class DFTRefiner:
    """Wrapper for running the tiered DFT composite workflow."""

    def __init__(self, solvent_name: Optional[str] = 'water', temp_k: float = 423.15,
                 use_explicit_solvent: bool = False, n_water: int = 3,
                 geometry_backend: str = 'pyscf', db_path: Optional[str] = None):
        self.solvent_name = solvent_name
        self.temp_k = temp_k # Default 150 C
        self.use_explicit_solvent = use_explicit_solvent
        self.n_water = n_water
        self.geometry_backend = geometry_backend.lower()
        self.db_path = db_path

        if self.db_path:
            from .results_db import ResultsDB
            self.db = ResultsDB(db_path=self.db_path)
        else:
            self.db = None

        self.solvation_engine = SolvationEngine()

        if self.geometry_backend == 'mace':
            try:
                from .mlp_optimizer import MLPOptimizer
                self.mlp_optimizer = MLPOptimizer()
            except ImportError:
                logger.info("WARNING: ML properties requested but not available. Falling back to pyscf.")
                self.geometry_backend = 'pyscf'
        else:
            self.mlp_optimizer = None

        # Lightweight MLP for DFT starting-point pre-relaxation
        self._prerelax_mlp = None
        if self.geometry_backend != 'mace':
            try:
                from .mlp_optimizer import MLPOptimizer
                self._prerelax_mlp = MLPOptimizer(model_name="small")
            except (ImportError, Exception):
                pass
        self._secondary_prerelax_mlp = None
        self._secondary_prerelax_mlp_error: Optional[str] = None

        try:
            from .ts_optimizer import TSOptimizer
            self.ts_optimizer = TSOptimizer()
        except Exception:
            self.ts_optimizer = None

        self.diffusion_engine = DiffusionTSEngine()

        try:
            from .mlp_barrier import MLPBarrier
            self.mlp_barrier = MLPBarrier()
        except ImportError:
            self.mlp_barrier = None

        self.opt_method = 'r2SCAN'
        self.opt_basis = 'def2-svp'
        self.refinement_method = 'wB97M-V'
        self.refinement_basis = 'def2-tzvp'
        self.verif_method = 'revDSD-PBEP86'
        self.verif_basis = 'def2-tzvp'
        self._optimization_checkpoint_context: Optional[Dict[str, Any]] = None

        n_threads = os.cpu_count() or 1
        os.environ["OMP_NUM_THREADS"] = str(n_threads)

    @staticmethod
    def _emit_progress(
        progress_callback: Optional[Callable[[str, Dict[str, Any]], None]],
        phase: str,
        **details: Any,
    ) -> None:
        if progress_callback is None:
            return
        progress_callback(phase, details)

    @staticmethod
    def _checkpoint_slug(value: str) -> str:
        chars: List[str] = []
        for char in value.lower():
            if char.isalnum():
                chars.append(char)
            else:
                chars.append("_")
        slug = "".join(chars)
        while "__" in slug:
            slug = slug.replace("__", "_")
        slug = slug.strip("_")
        return slug or "step"

    def _set_optimization_checkpoint_context(
        self,
        *,
        checkpoint_dir: Optional[Path],
        phase_prefix: str,
        resume_state: Optional[Dict[str, Any]] = None,
    ) -> None:
        if checkpoint_dir is None:
            self._optimization_checkpoint_context = None
            return
        self._optimization_checkpoint_context = {
            "dir": str(checkpoint_dir),
            "phase_prefix": phase_prefix,
            "resume_state": dict(resume_state or {}),
        }

    def _get_optimization_checkpoint_context(self, phase_prefix: str) -> Optional[Dict[str, Any]]:
        ctx = self._optimization_checkpoint_context
        if not isinstance(ctx, dict):
            return None
        if ctx.get("phase_prefix") != phase_prefix:
            return None
        dir_value = ctx.get("dir")
        if not dir_value:
            return None
        return {
            "dir": Path(str(dir_value)),
            "resume_state": dict(ctx.get("resume_state") or {}),
        }

    @staticmethod
    def _is_second_order_scf_wrapper(mf: Any) -> bool:
        class_name = type(mf).__name__.lower()
        if "newton" in class_name or "soscf" in class_name:
            return True
        wrapped = getattr(mf, "_scf", None)
        if wrapped is None:
            return False
        wrapped_name = type(wrapped).__name__.lower()
        return "newton" in wrapped_name or "soscf" in wrapped_name

    @staticmethod
    def _apply_fermi_smearing(mf: Any, scf_module: Any, smearing_sigma: Optional[float]) -> Any:
        if smearing_sigma is None:
            return mf
        if DFTRefiner._is_second_order_scf_wrapper(mf):
            logger.warning(
                "    WARNING: Refusing to apply Fermi smearing to an SOSCF/Newton-wrapped SCF object. "
                "Smearing must be attached to the base engine before second-order SCF is enabled."
            )
            return mf
        try:
            mf = scf_module.addons.smearing_(mf, sigma=smearing_sigma, method='fermi')
            setattr(mf, "_maillard_fermi_smearing_sigma", smearing_sigma)
            logger.info(f"    [SCF] Enabled Fermi smearing with sigma={smearing_sigma:.4f} Ha.")
        except Exception as exc:
            logger.warning(f"    WARNING: Fermi smearing setup failed ({exc}). Continuing without smearing.")
        return mf

    def _best_available_geometry_fallback(
        self,
        xyz_content: str,
        *,
        charge: int,
        spin: int,
        basis: str,
        label: str,
    ) -> Tuple[bool, Any, str]:
        logger.warning(
            f"    [LADDER] {label}: all geometry-optimization strategies exhausted. "
            "Falling back to the best available geometry for single-point evaluation."
        )
        mol = self._setup_mol(xyz_content, charge=charge, spin=spin, basis=basis)
        return False, mol, xyz_content

    @staticmethod
    def _orbital_gap_instability(
        mf: Any,
        *,
        small_gap_ev: float = 0.05,
    ) -> Optional[Dict[str, float | bool]]:
        if getattr(getattr(mf, 'mol', None), 'spin', 0) != 0:
            return None
        mo_energy = getattr(mf, 'mo_energy', None)
        mo_occ = getattr(mf, 'mo_occ', None)
        if mo_energy is None or mo_occ is None:
            return None
        if isinstance(mo_energy, (tuple, list)) or isinstance(mo_occ, (tuple, list)):
            return None

        occ_mask = np.asarray(mo_occ) > 0
        vir_mask = ~occ_mask
        if not occ_mask.any() or not vir_mask.any():
            return None

        occ_energies = np.asarray(mo_energy)[occ_mask]
        vir_energies = np.asarray(mo_energy)[vir_mask]
        homo = float(np.max(occ_energies))
        lumo = float(np.min(vir_energies))
        gap_h = lumo - homo
        gap_ev = gap_h * 27.211386245988
        return {
            'gap_ev': gap_ev,
            'inverted': gap_ev < 0.0,
            'near_degenerate': gap_ev <= small_gap_ev,
        }

    @staticmethod
    def _is_repeated_force_plateau(
        previous_positions_ang: Optional[np.ndarray],
        positions_ang: np.ndarray,
        previous_forces_ev_ang: Optional[np.ndarray],
        forces_ev_ang: np.ndarray,
        previous_energy_ev: Optional[float],
        energy_ev: float,
        *,
        position_tol_ang: float = 1e-6,
        force_tol_ev_ang: float = 5e-3,
        energy_tol_ev: float = 1e-4,
    ) -> bool:
        if previous_positions_ang is None or previous_forces_ev_ang is None or previous_energy_ev is None:
            return False
        same_geometry = float(np.max(np.abs(positions_ang - previous_positions_ang))) <= position_tol_ang
        same_forces = float(np.max(np.abs(forces_ev_ang - previous_forces_ev_ang))) <= force_tol_ev_ang
        same_energy = abs(energy_ev - previous_energy_ev) <= energy_tol_ev
        return same_geometry and same_forces and same_energy

    def _primary_prerelax_backend(self) -> Any:
        if self.geometry_backend == 'mace' and self.mlp_optimizer is not None:
            return self.mlp_optimizer
        return self._prerelax_mlp

    def _secondary_prerelax_backend(self) -> Any:
        if self._secondary_prerelax_mlp is not None:
            return self._secondary_prerelax_mlp
        if self._secondary_prerelax_mlp_error is not None:
            return None
        try:
            from .ubio_optimizer import UBioOptimizer

            self._secondary_prerelax_mlp = UBioOptimizer(device="cpu")
        except Exception as exc:
            self._secondary_prerelax_mlp_error = str(exc)
            logger.warning(
                f"    [MLP FALLBACK] Secondary MLP backend (UBio) unavailable: {exc}. "
                "Only the primary MLP will be tried for pre-relaxation."
            )
            return None
        return self._secondary_prerelax_mlp

    def _iter_prerelax_backends(self):
        primary = self._primary_prerelax_backend()
        if primary is not None:
            yield "primary_mlp", primary
        secondary = self._secondary_prerelax_backend()
        if secondary is not None:
            yield "secondary_mlp", secondary

    def _setup_mol(self, xyz_content: str, charge: int = 0, spin: int = 0, basis: str = 'def2-svp') -> Any:
        """Initialize PySCF GTO Mole object from XYZ."""
        atom_string = _extract_atom_string(xyz_content)

        mol = gto.M(
            atom=atom_string,
            basis=basis,
            charge=charge,
            spin=spin,
            verbose=3 # 0=silent, 3=warnings, 4=info
        )
        return mol

    def _mlp_prerelax(self, xyz_content: str) -> str:
        """Quick MLP geometry pre-relaxation to improve the DFT starting point."""
        try:
            from .mlp_optimizer import compute_geometry_drift_metrics
        except Exception:
            compute_geometry_drift_metrics = None

        for backend_label, backend in self._iter_prerelax_backends():
            try:
                logger.info(f"    Running {backend_label} pre-relaxation...")
                relaxed = backend.optimize_geometry(
                    xyz_content, fmax=0.1, max_steps=200, drift_threshold=1.5
                )
                if compute_geometry_drift_metrics is None:
                    if relaxed.strip() != xyz_content.strip():
                        logger.info(f"    {backend_label} pre-relaxation complete.")
                        return relaxed
                else:
                    metrics = compute_geometry_drift_metrics(xyz_content, relaxed)
                    if float(metrics.get("max_atom_displacement_angstrom") or 0.0) > 1e-4:
                        logger.info(f"    {backend_label} pre-relaxation complete.")
                        return relaxed
                logger.info(f"    {backend_label} produced no usable geometry update.")
            except Exception as exc:
                logger.info(f"    WARNING: {backend_label} pre-relaxation failed ({exc}).")
        return xyz_content

    def _count_mlp_imaginary_modes(self, xyz_content: str) -> Optional[int]:
        mlp_backend = self._primary_prerelax_backend()
        if mlp_backend is None or not hasattr(mlp_backend, 'count_imaginary_modes'):
            return None
        try:
            return int(mlp_backend.count_imaginary_modes(xyz_content))
        except Exception as exc:
            logger.info(f"    WARNING: MLP imaginary-mode check failed ({exc}).")
            return None

    def _build_mf(
        self,
        mol: Any,
        xc_method: str = 'r2SCAN',
        use_solvent: bool = True,
        conv_tol: float = 1e-7,
        harden_scf: bool = False,
        use_soscf: bool = False,
        level_shift: float = 0.0,
        smearing_sigma: Optional[float] = None,
    ):
        """Build the Mean-Field object with specified XC, implicit solvent, and stability settings."""
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

        # [PERFORMANCE] Enable density fitting (RI-J) to accelerate Coulomb
        # integrals.  For hybrid/range-separated functionals (e.g. wB97M-V)
        # this still helps the J part; exact exchange stays exact.
        try:
            mf = mf.density_fit()
        except Exception as exc:
            logger.warning(f"    [RI-J] density_fit() unavailable ({exc}); proceeding without RI acceleration.")

        if use_solvent and self.solvent_name is not None:
            mf = mf.ddCOSMO()
            if self.solvent_name.lower() == 'water':
                mf.with_solvent.eps = 78.3553
        elif hasattr(mf, 'with_solvent'):
            # Failsafe: ensure no solvent is attached if use_solvent is False
            delattr(mf, 'with_solvent')
                
        mf.conv_tol = conv_tol
        mf.max_cycle = 200
        
        if level_shift > 0:
            mf.level_shift = level_shift
            
        # [L2 ROBUSTNESS] Explicitly set memory cap to 8GB to prevent OOM
        # in the container and host system instability on 16GB machines.
        if use_soscf:
            mf.max_memory = 8000
        
        # [SCRATCH MANAGEMENT] Prevent PySCF from creating .chk files in the 
        mf.chkfile = None
        if harden_scf:
            if hasattr(mf, 'level_shift') and level_shift == 0:
                mf.level_shift = 0.5
            if hasattr(mf, 'damp'):
                mf.damp = 0.2
            if hasattr(mf, 'diis_space'):
                mf.diis_space = 12
            if hasattr(mf, 'direct_scf_tol'):
                mf.direct_scf_tol = 1e-13
            mf.init_guess = 'atom'

        # [SCF SAFETY CONTRACT] Smearing + SOSCF/Newton is fundamentally incompatible
        # with PySCF's analytical gradient/Hessian path used by Sella and geomeTRIC:
        # Fermi smearing emits fractional `mo_occ` values, but the downstream
        # derivative code (especially when combined with density_fit + ddCOSMO)
        # assumes integer occupations when building occ/vir boolean masks. This
        # manifests as:
        #   * UKS: `TypeError: 'float' object is not subscriptable` in newton_ah.rotate_mo
        #   * RKS (in full TS pipeline): `ValueError: NumPy boolean array indexing
        #     assignment cannot assign N input values to M output values where the
        #     mask is true` in the Sella/ASE bridge or geomeTRIC.
        # We refuse to combine the two regardless of spin. Smearing is kept as a
        # pure DIIS convergence aid; plateau-breaking with Newton uses level_shift
        # instead (see the TS SCF strategy ladder).
        if use_soscf and smearing_sigma is not None:
            logger.warning(
                "    [SCF] Refusing to combine Fermi smearing with SOSCF/Newton — "
                "PySCF's analytical derivatives do not support fractional occupations. "
                "Dropping smearing; Newton will run with level_shift to break plateaus."
            )
            smearing_sigma = None

        # Apply smearing to the base SCF engine before enabling Newton/SOSCF.
        # PySCF's second-order wrapper is not a safe target for post-hoc smearing.
        mf = self._apply_fermi_smearing(mf, scf, smearing_sigma)

        if use_soscf:
            mf = mf.newton()
            mf.max_cycle = 50
            if hasattr(mf, 'level_shift') and level_shift == 0:
                mf.level_shift = 0.2
        else:
            mf.max_cycle = 200
        
        return mf

    def recover_ts_guess_from_forming_bond(
        self,
        reactant_xyz: str,
        ts_xyz: str,
        *,
        charge: int = 0,
        spin: int = 0,
        reaction_meta: Optional[Dict[str, Any]] = None,
        force_scan: bool = False,
        progress_callback: Optional[Callable[[str, Dict[str, Any]], None]] = None,
    ) -> Dict[str, Any]:
        forming_bond_atoms = list((reaction_meta or {}).get('forming_bond_atoms', []))
        if len(forming_bond_atoms) != 2:
            return {"used": False, "strategy": "unavailable", "ts_xyz": ts_xyz}

        initial_modes = self._count_mlp_imaginary_modes(ts_xyz)
        if initial_modes is not None and initial_modes == 1 and not force_scan:
            return {
                "used": False,
                "strategy": "existing_ts_guess",
                "ts_xyz": ts_xyz,
                "imaginary_mode_count": initial_modes,
            }

        primary_mlp = self._primary_prerelax_backend()
        if primary_mlp is None or not hasattr(primary_mlp, 'find_ts_relaxed_scan'):
            secondary_payload = self._recover_ts_guess_with_secondary_mlp(
                ts_xyz,
                initial_modes=initial_modes,
                progress_callback=progress_callback,
            )
            if secondary_payload is not None:
                return secondary_payload
            return {
                "used": False,
                "strategy": "unavailable",
                "ts_xyz": ts_xyz,
                "imaginary_mode_count": initial_modes,
            }

        try:
            self._emit_progress(progress_callback, "ts_relaxed_scan", label="MLP forming-bond relaxed scan")
            scan_xyz, scan_meta = primary_mlp.find_ts_relaxed_scan(
                reactant_xyz,
                (int(forming_bond_atoms[0]), int(forming_bond_atoms[1])),
            )
            repaired = repair_steric_clash(scan_xyz)
            if repaired.get("repaired") and not repaired.get("repair_complete", True):
                logger.warning("    [TS RECOVERY] Steric repair incomplete after max iterations. Flagging.")
            candidate_xyz = repaired["xyz"] if repaired.get("repaired") else scan_xyz
            return {
                "used": True,
                "strategy": "forming_bond_relaxed_scan",
                "ts_xyz": candidate_xyz,
                "imaginary_mode_count": self._count_mlp_imaginary_modes(candidate_xyz),
                "initial_imaginary_mode_count": initial_modes,
                "scan": scan_meta,
                "steric_repair": repaired if repaired.get("repaired") else None,
            }
        except Exception as exc:
            logger.info(f"    WARNING: forming-bond relaxed scan failed ({exc}).")
            secondary_payload = self._recover_ts_guess_with_secondary_mlp(
                ts_xyz,
                initial_modes=initial_modes,
                progress_callback=progress_callback,
            )
            if secondary_payload is not None:
                return secondary_payload
            return {
                "used": False,
                "strategy": "forming_bond_relaxed_scan_failed",
                "ts_xyz": ts_xyz,
                "imaginary_mode_count": initial_modes,
                "error": str(exc),
            }

    def _recover_ts_guess_with_secondary_mlp(
        self,
        ts_xyz: str,
        *,
        initial_modes: Optional[int],
        progress_callback: Optional[Callable[[str, Dict[str, Any]], None]] = None,
    ) -> Optional[Dict[str, Any]]:
        secondary_backend = self._secondary_prerelax_backend()
        if secondary_backend is None or not hasattr(secondary_backend, 'optimize_ts'):
            return None
        try:
            self._emit_progress(progress_callback, "ts_secondary_mlp_seed", label="Secondary MLP TS seed recovery")
            recovered_xyz = secondary_backend.optimize_ts(ts_xyz, fmax=0.05, max_steps=100)
            repaired = repair_steric_clash(recovered_xyz)
            candidate_xyz = repaired["xyz"] if repaired.get("repaired") else recovered_xyz
            return {
                "used": True,
                "strategy": "secondary_mlp_ts_seed",
                "ts_xyz": candidate_xyz,
                "imaginary_mode_count": self._count_mlp_imaginary_modes(candidate_xyz),
                "initial_imaginary_mode_count": initial_modes,
                "scan": None,
                "steric_repair": repaired if repaired.get("repaired") else None,
            }
        except Exception as exc:
            logger.info(f"    WARNING: secondary MLP TS recovery failed ({exc}).")
            return None

    @staticmethod
    def _is_scf_gradient_failure(exc: Exception) -> bool:
        text = str(exc).lower()
        markers = (
            'nuclear gradients',
            'scf not converged',
            'not converged',
            'rks_scanner',
            # PySCF smearing + Newton + analytical-derivative pathologies:
            # the occ/vir boolean masks in the gradient/Hessian code get
            # inconsistent shapes under fractional occupations. Treat as
            # an escalable SCF failure so the next ladder stage gets a shot.
            'boolean array indexing assignment',
            "'float' object is not subscriptable",
        )
        return any(marker in text for marker in markers)

    def _log_spin_integrity(self, mf: Any) -> None:
        """Calculate and log spin contamination for open-shell systems.
        Theoretical <S^2> for a doublet is 0.75, triplet is 2.0, etc.
        """
        if mf.mol.spin == 0:
            return
            
        try:
            ss, mult = mf.spin_square()
            s_ideal = mf.mol.spin / 2.0
            ss_ideal = s_ideal * (s_ideal + 1.0)
            deviation = (ss - ss_ideal) / ss_ideal if ss_ideal > 0 else ss
            
            message = f"    [SPIN CHECK] <S^2> = {ss:.4f} (Ideal: {ss_ideal:.4f}, Deviation: {deviation*100:.1f}%)"
            if abs(deviation) > 0.1:
                logger.warning(f"{message} -- [SCIENTIFIC CAUTION] High Spin Contamination detected!")
            else:
                logger.info(message)
        except Exception as e:
            logger.warning(f"    [SPIN CHECK] Could not calculate <S^2>: {e}")

    def _run_xtb_optimization(self, xyz_content: str, charge: int, spin: int) -> str:
        """Delegate to :mod:`xtb_backend`."""
        return _xtb.run_xtb_optimization(xyz_content, charge, spin)

    def _run_xtb_ts_safe_prerelax(
        self,
        xyz_content: str,
        charge: int,
        spin: int,
        max_cycles: int = 20,
    ) -> str:
        """Delegate to :mod:`xtb_backend`."""
        return _xtb.run_xtb_ts_safe_prerelax(xyz_content, charge, spin, max_cycles=max_cycles)

    def _optimize_with_geometric_kernel(
        self,
        mf: Any,
        *,
        label: str,
        max_steps: int,
        is_ts_fallback: bool,
        eff_use_explicit: bool,
        n_atoms_solute: int,
        coords: str = "tric",
    ) -> Tuple[bool, Any, Optional[str]]:
        try:
            initial_gradient = mf.nuc_grad_method().kernel()
            max_gradient = float(np.abs(initial_gradient).max())
            logger.info(f"    [GEOMETRIC] Initial max gradient: {max_gradient:.2f} Ha/Bohr")
            if max_gradient > 50000.0:
                raise RuntimeError(f"Initial gradient explosion detected ({max_gradient:.2f} Ha/Bohr)")
            if coords == "tric" and max_gradient > 5000.0:
                logger.warning("    [GEOMETRIC] Gradient too large for TRIC. Forcing immediate Cartesian fallback.")
                return False, mf.mol, None
        except RuntimeError:
            raise
        except Exception as exc:
            logger.info(f"    WARNING: Could not evaluate initial gradient guard ({exc}). Proceeding with geomeTRIC.")
        
        start_time = time.time()
        timeout_seconds = 20 * 3600  # 20 hours soft kill
        
        # [SOTA GUARD] Memory limit detection (cgroup v2)
        mem_limit = 10 * 1024 * 1024 * 1024 # 10GB default
        if os.path.exists("/sys/fs/cgroup/memory.max"):
            try:
                with open("/sys/fs/cgroup/memory.max", "r") as f:
                    val = f.read().strip()
                    if val != "max":
                        mem_limit = int(val)
            except Exception:
                pass

        # [SOTA GUARD] Gradient plateau detection — abort early when the
        # optimizer is making no progress so the ladder can escalate.
        _plateau_energies: list = []
        _PLATEAU_LOOKBACK = 5  # cycles to compare
        _PLATEAU_ABS_TOL = 1e-4  # Ha — if energy spread < this over window, stalled

        def _watchdog_callback(envs):
            # 0. Energy Plateau Detection
            # PySCF's PySCFEngine.calc_new calls callback(locals()) which
            # includes 'energy' (not 'e_tot') as the SCF total energy.
            e_tot = None
            if isinstance(envs, dict):
                e_tot = envs.get("energy", envs.get("e_tot", None))
            if e_tot is not None:
                _plateau_energies.append(float(e_tot))
                if len(_plateau_energies) >= _PLATEAU_LOOKBACK:
                    recent = _plateau_energies[-_PLATEAU_LOOKBACK:]
                    spread = max(recent) - min(recent)
                    if spread < _PLATEAU_ABS_TOL:
                        logger.warning(
                            f"    [WATCHDOG] Energy plateau detected: "
                            f"last {_PLATEAU_LOOKBACK} energies spread={spread:.2e} Ha "
                            f"(< {_PLATEAU_ABS_TOL:.0e} threshold). Aborting to let ladder escalate."
                        )
                        raise OptimizationTimeoutError(
                            f"Energy plateau (spread={spread:.2e} Ha over {_PLATEAU_LOOKBACK} cycles)"
                        )

            # 1. Time Check
            if time.time() - start_time > timeout_seconds:
                raise OptimizationTimeoutError("Internal watchdog time limit reached")
            
            # 2. Memory Check (SOTA Guard)
            if os.path.exists("/sys/fs/cgroup/memory.current"):
                try:
                    with open("/sys/fs/cgroup/memory.current", "r") as f:
                        usage = int(f.read().strip())
                    if usage > 0.9 * mem_limit:
                        logger.warning(f"    [WATCHDOG ERROR] Memory usage ({usage/1e9:.1f}GB) approaching limit ({mem_limit/1e9:.1f}GB). Emergency exit.")
                        raise OptimizationTimeoutError("Internal watchdog memory limit reached")
                except OptimizationTimeoutError:
                    raise
                except Exception:
                    pass

        # [SOTA GUARD] Persistent Workspace
        # Instead of TemporaryDirectory, we use a target-specific folder 
        # so that trajectories (optim.xyz) survive system reboots.
        run_label = f"{label}_{'ts' if is_ts_fallback else 'reactant'}"
        work_dir = Path("data/geometries/dft_runs") / run_label
        work_dir.mkdir(parents=True, exist_ok=True)
        
        pwd = os.getcwd()
        try:
            os.chdir(work_dir)
            geome_kwargs = {
                'maxsteps': max_steps,
                'trajectory': 'optim.xyz',
                'callback': _watchdog_callback,
                'coordsys': coords
            }
            if is_ts_fallback:
                geome_kwargs['transition'] = True

            if eff_use_explicit and is_ts_fallback:
                logger.info("    Pre-relaxing solvent molecules around the frozen core...")
                with open("constraints.txt", "w") as f:
                    f.write("$freeze\n")
                    f.write(f"xyz 1-{n_atoms_solute}\n")

                pre_relax_kwargs = {'maxsteps': 50, 'constraints': "constraints.txt"}
                conv, mol_opt = geometric_solver.kernel(mf, **pre_relax_kwargs)
            else:
                try:
                    logger.info(f"    [GEOMETRIC] System coordinates: {coords.upper()}")
                    conv, mol_opt = geometric_solver.kernel(mf, **geome_kwargs)
                    return conv, mol_opt, None
                except (OptimizationTimeoutError, RuntimeError) as exc:
                    is_scf_fail = isinstance(exc, RuntimeError) and self._is_scf_gradient_failure(exc)
                    if isinstance(exc, OptimizationTimeoutError):
                        logger.warning(f"    [WATCHDOG] Limit reached in {work_dir}. Terminating early to preserve best guess.")
                    elif is_scf_fail:
                        logger.warning(f"    [GEOMETRIC] SCF did not converge during optimization: {exc}")
                        logger.info("    [GEOMETRIC] Recovering best geometry from trajectory...")
                    else:
                        raise
                    opt_xyz = None
                    if os.path.exists("optim.xyz"):
                        with open("optim.xyz", "r") as f:
                            lines = f.readlines()
                        if lines:
                            natoms = int(lines[0].strip())
                            if len(lines) >= natoms + 2:
                                opt_xyz = "".join(lines[-(natoms+2):])
                                logger.info(f"    [GEOMETRIC] Recovered best-guess frame from {work_dir}/optim.xyz")
                                return False, mf.mol, opt_xyz
                    return False, mf.mol, None
        finally:
            os.chdir(pwd)

        return conv, mol_opt, None

    def _optimize_with_pyscf_backend(
        self,
        xyz_content: str,
        charge: int,
        spin: int,
        is_ts: bool,
        max_steps: int,
        eff_use_explicit: bool,
        n_atoms_solute: int,
        *,
        label: str = "default",
        harden_scf: bool = False,
        use_solvent_optimization: bool = False,
        use_soscf: bool = False,
        smearing_sigma: Optional[float] = None,
        level_shift: float = 0.0,
        forming_bond_pair: Optional[Tuple[int, int]] = None,
        ts_v0: Optional[np.ndarray] = None,
    ) -> Tuple[bool, Any, str]:
        mol = self._setup_mol(xyz_content, charge, spin, basis=self.opt_basis)
        mf = self._build_mf(
            mol,
            xc_method=self.opt_method,
            use_solvent=use_solvent_optimization,
            conv_tol=1e-6,
            harden_scf=harden_scf,
            use_soscf=use_soscf,
            smearing_sigma=smearing_sigma,
            level_shift=level_shift,
        )
        # [L3 ROBUSTNESS] During geometry optimization, cap SCF cycles per
        # geometry step at 50 (not 200).  If the electronic structure hasn't
        # converged in 50 cycles it is oscillating; geomeTRIC can still use
        # the approximate gradient to make progress.
        if not use_soscf:
            mf.max_cycle = 50

        conv = False
        if is_ts and self.ts_optimizer:
            logger.info(">>> [Phase 11] Running Sella eigenvector-following TS search...")
            try:
                from ase import Atoms as ASEAtoms
                from ase.calculators.calculator import Calculator, all_changes
                from pyscf.data import nist as pyscf_nist
                from .ts_optimizer import TSOptimizationStall, TSPlateauConverged

                # Capture current SCF settings in a factory so every Sella
                # step builds a fully-configured MF (solvent, smearing, etc.)
                _refiner = self
                _build_kwargs = dict(
                    xc_method=self.opt_method,
                    use_solvent=use_solvent_optimization,
                    conv_tol=1e-6,
                    harden_scf=harden_scf,
                    use_soscf=use_soscf,
                    smearing_sigma=smearing_sigma,
                    level_shift=level_shift,
                )

                class BuiltinPySCFCalc(Calculator):
                    implemented_properties = ['energy', 'forces', 'hessian']

                    def __init__(self, mol):
                        super().__init__()
                        self.mol = mol
                        self.consecutive_scf_failures = 0
                        self.repeated_plateau_evals = 0
                        self._last_positions_ang: Optional[np.ndarray] = None
                        self._last_forces_ev_ang: Optional[np.ndarray] = None
                        self._last_energy_ev: Optional[float] = None
                        self.orbital_instability_evals = 0
                        # Best-so-far tracking + flat-line plateau detector
                        self.best_fmax_ev_ang: float = float("inf")
                        self.best_positions_ang: Optional[np.ndarray] = None
                        self.best_energy_ev: Optional[float] = None
                        self._fmax_history: list[float] = []
                        # If max-min spread of last `_PLATEAU_WINDOW` fmax
                        # values is below `_PLATEAU_REL_TOL` of the mean, the
                        # search is on a numerical plateau and is aborted.
                        self._PLATEAU_WINDOW = 30
                        self._PLATEAU_REL_TOL = 1e-2
                        # P2: SCF density-matrix guess reuse across Sella steps.
                        # Geometry changes are sub-Å between consecutive evals;
                        # using the previous DM as initial guess typically cuts
                        # SCF iterations 2-3x.
                        self._last_dm: Optional[np.ndarray] = None
                        # P3: when an inverted/near-degenerate HOMO-LUMO gap is
                        # observed, surface a hint to the outer ladder so it
                        # can skip ahead to a level-shifted strategy.
                        self.escalate_to_level_shift: bool = False

                    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
                        super().calculate(atoms, properties, system_changes)
                        self.mol.set_geom_(atoms.positions / pyscf_nist.BOHR, unit='Bohr')
                        mf = _refiner._build_mf(self.mol, **_build_kwargs)
                        if not _build_kwargs.get("use_soscf"):
                            mf.max_cycle = 50

                        if self._last_dm is not None:
                            try:
                                mf.kernel(dm0=self._last_dm)
                            except Exception:
                                # Stale DM (e.g. basis change between strategies):
                                # fall back to a fresh guess.
                                self._last_dm = None
                                mf.kernel()
                        else:
                            mf.kernel()

                        if mf.converged:
                            try:
                                self._last_dm = mf.make_rdm1()
                            except Exception:
                                self._last_dm = None
                        else:
                            # Drop stale DM if SCF failed; next step starts fresh
                            self._last_dm = None

                        if not mf.converged:
                            self.consecutive_scf_failures += 1
                            logger.warning(
                                f"    [SELLA-SCF] SCF did not converge "
                                f"({self.consecutive_scf_failures} consecutive failures)"
                            )
                            if self.consecutive_scf_failures >= 5:
                                raise TSOptimizationStall(
                                    f"SCF failed to converge in {self.consecutive_scf_failures} "
                                    f"consecutive Sella steps"
                                )
                        else:
                            self.consecutive_scf_failures = 0

                        instability = _refiner._orbital_gap_instability(mf)
                        if instability is not None:
                            gap_ev = float(instability['gap_ev'])
                            if bool(instability['inverted']) or bool(instability['near_degenerate']):
                                self.orbital_instability_evals += 1
                                logger.warning(
                                    "    [SELLA-SCF] Closed-shell orbital gap instability detected "
                                    f"(gap={gap_ev:.4f} eV, inverted={bool(instability['inverted'])}, "
                                    f"near_degenerate={bool(instability['near_degenerate'])}) "
                                    f"[{self.orbital_instability_evals} consecutive evaluations]"
                                )
                                if self.orbital_instability_evals >= 3:
                                    # P3: hint to outer ladder to skip ahead to
                                    # a level-shifted strategy rather than the
                                    # next sequential one.
                                    self.escalate_to_level_shift = True
                                    raise TSOptimizationStall(
                                        "Closed-shell SCF repeatedly converged to an inverted or near-degenerate "
                                        "HOMO-LUMO gap; escalating optimizer strategy"
                                    )
                            else:
                                self.orbital_instability_evals = 0

                        energy_ev = float(mf.e_tot) * pyscf_nist.HARTREE2EV
                        self.results['energy'] = energy_ev

                        if 'forces' in properties:
                            g = mf.nuc_grad_method().kernel()
                            forces_ev_ang = -g * pyscf_nist.HARTREE2EV / pyscf_nist.BOHR
                            self.results['forces'] = forces_ev_ang

                            positions_ang = np.array(atoms.positions, copy=True)
                            fmax_now = float(
                                np.sqrt((forces_ev_ang ** 2).sum(axis=1).max())
                            )

                            # Track best-so-far so the caller can recover the
                            # lowest-force geometry even if the search stalls.
                            if fmax_now < self.best_fmax_ev_ang:
                                self.best_fmax_ev_ang = fmax_now
                                self.best_positions_ang = positions_ang
                                self.best_energy_ev = energy_ev

                            # P1: bad TS guess gate.  If the FIRST DFT force
                            # evaluation already shows fmax >> typical TS
                            # gradient (~0.1 eV/Å), the geometry is far from
                            # any saddle point.  Warn loudly; if it persists
                            # after a few evals, abort to avoid wasting hours.
                            if len(self._fmax_history) == 0 and fmax_now > 1.0:
                                logger.warning(
                                    f"    [SELLA-GUESS] Initial DFT fmax="
                                    f"{fmax_now:.3f} eV/Å is suspiciously large "
                                    "(target 0.05). TS guess is likely poor; "
                                    "monitoring next 3 evals before deciding."
                                )
                            if (
                                len(self._fmax_history) == 3
                                and min(self._fmax_history) > 2.0
                            ):
                                raise TSPlateauConverged(
                                    f"Bad TS guess: fmax min over first 3 evals "
                                    f"is {min(self._fmax_history):.3f} eV/Å "
                                    "(>>0.1 expected near a saddle). Aborting "
                                    "DFT before more compute is wasted."
                                )

                            # Flat-line plateau detector. When fmax stops
                            # changing meaningfully across a window of evals,
                            # the optimizer is wasting cycles — abort and let
                            # the caller use the best geometry seen so far.
                            self._fmax_history.append(fmax_now)
                            if len(self._fmax_history) >= self._PLATEAU_WINDOW:
                                window = self._fmax_history[-self._PLATEAU_WINDOW:]
                                wmin = min(window)
                                wmax = max(window)
                                wmean = sum(window) / len(window)
                                rel_spread = (wmax - wmin) / max(wmean, 1e-12)
                                if rel_spread <= self._PLATEAU_REL_TOL:
                                    logger.warning(
                                        f"    [SELLA-PLATEAU] fmax flat at "
                                        f"{wmean:.5f} eV/Å over "
                                        f"{self._PLATEAU_WINDOW} evals "
                                        f"(rel spread {rel_spread:.2e} <= "
                                        f"{self._PLATEAU_REL_TOL:.0e}); aborting "
                                        "TS search and using best geometry."
                                    )
                                    raise TSPlateauConverged(
                                        f"Sella TS search plateaued at fmax≈"
                                        f"{wmean:.5f} eV/Å (best ever "
                                        f"{self.best_fmax_ev_ang:.5f}); "
                                        "best geometry will be used."
                                    )

                            if _refiner._is_repeated_force_plateau(
                                self._last_positions_ang,
                                positions_ang,
                                self._last_forces_ev_ang,
                                forces_ev_ang,
                                self._last_energy_ev,
                                energy_ev,
                            ):
                                self.repeated_plateau_evals += 1
                            else:
                                self.repeated_plateau_evals = 0

                            self._last_positions_ang = positions_ang
                            self._last_forces_ev_ang = np.array(forces_ev_ang, copy=True)
                            self._last_energy_ev = energy_ev

                            if self.repeated_plateau_evals >= 6:
                                raise TSPlateauConverged(
                                    "Repeated near-identical DFT force evaluations at unchanged geometry; "
                                    "aborting Sella and using best geometry."
                                )

                        if 'hessian' in properties:
                            h_obj = mf.Hessian()
                            h = h_obj.kernel()
                            natm = self.mol.natm
                            h_reshaped = h.transpose(0, 2, 1, 3).reshape(3 * natm, 3 * natm)
                            self.results['hessian'] = h_reshaped * pyscf_nist.HARTREE2EV / (pyscf_nist.BOHR**2)

                    def get_hessian(self, atoms):
                        self.calculate(atoms, properties=['hessian'])
                        return self.results['hessian']

                coords_bohr = mol.atom_coords()
                symbols = [mol.atom_symbol(i) for i in range(mol.natm)]
                positions_ang = coords_bohr * pyscf_nist.BOHR
                ase_atoms = ASEAtoms(symbols=symbols, positions=positions_ang)

                calc = BuiltinPySCFCalc(mol=mol)
                ase_atoms.calc = calc

                # [P0b] DFT constrained prerelax: relax backbone at DFT
                # level while keeping the forming bond fixed.  This
                # dramatically reduces initial DFT fmax (typically
                # 0.8→0.2 eV/Å) so Sella can make progress.
                if forming_bond_pair is not None and len(forming_bond_pair) == 2:
                    from ase.constraints import FixBondLength
                    from ase.optimize import BFGS as ASE_BFGS
                    fb_i, fb_j = forming_bond_pair
                    d_ij = float(ase_atoms.get_distance(fb_i, fb_j))
                    logger.info(
                        f"    [DFT-PRERELAX] Constrained prerelax: fixing "
                        f"bond ({fb_i},{fb_j}) at d={d_ij:.2f} Å, relaxing "
                        f"backbone at {self.opt_method}..."
                    )
                    prerelax_calc = BuiltinPySCFCalc(mol=mol)
                    prerelax_calc._PLATEAU_WINDOW = 999  # disable plateau
                    ase_atoms.set_constraint(FixBondLength(fb_i, fb_j))
                    ase_atoms.calc = prerelax_calc
                    prerelax_dyn = ASE_BFGS(ase_atoms, logfile=None)
                    prerelax_nsteps = 0
                    try:
                        prerelax_dyn.run(fmax=0.2, steps=25)
                        prerelax_nsteps = prerelax_dyn.nsteps
                    except Exception as pre_exc:
                        logger.warning(
                            f"    [DFT-PRERELAX] Prerelax stopped: {pre_exc}"
                        )
                    forces_now = ase_atoms.get_forces()
                    fmax_after = float(
                        np.sqrt((forces_now ** 2).sum(axis=1).max())
                    )
                    logger.info(
                        f"    [DFT-PRERELAX] Done: fmax={fmax_after:.3f} eV/Å "
                        f"after {prerelax_nsteps} steps"
                    )
                    # Remove constraint, create fresh Sella calculator
                    ase_atoms.set_constraint()
                    mol.set_geom_(
                        ase_atoms.positions / pyscf_nist.BOHR, unit='Bohr'
                    )
                    calc = BuiltinPySCFCalc(mol=mol)
                    ase_atoms.calc = calc

                sella_traj_dir = str(Path("data/geometries/dft_runs") / f"{label}_ts_sella")
                atoms = self.ts_optimizer.find_ts(
                    ase_atoms, calc,
                    trajectory_dir=sella_traj_dir,
                    timeout_seconds=2 * 3600,
                    delta0=0.05,
                    v0=ts_v0,
                    internal=ts_v0 is not None,
                )
                if self.ts_optimizer.is_converged(atoms):
                    from ase.io import write as ase_write
                    import io
                    with io.StringIO() as f_xyz:
                        ase_write(f_xyz, atoms, format='xyz')
                        opt_xyz = f_xyz.getvalue()
                    mol_opt = self._setup_mol(opt_xyz, charge, spin, basis=self.opt_basis)
                    conv = True
                else:
                    logger.info("    WARNING: Sella failed to converge. Falling back to geomeTRIC...")
            except TSPlateauConverged as plateau_exc:
                # Numerical plateau: stop ALL geometry optimization and use the
                # best geometry seen so far. Skip the geomeTRIC fallback (which
                # would also stall on the same near-saddle PES).
                from ase.io import write as ase_write
                import io
                if calc.best_positions_ang is not None:
                    ase_atoms.set_positions(calc.best_positions_ang)
                logger.warning(
                    f"    [SELLA-PLATEAU] {plateau_exc} Best fmax="
                    f"{calc.best_fmax_ev_ang:.5f} eV/Å; using best geometry "
                    "directly for single-point + Hessian (no geomeTRIC fallback)."
                )
                with io.StringIO() as f_xyz:
                    ase_write(f_xyz, ase_atoms, format='xyz')
                    plateau_xyz = f_xyz.getvalue()
                # Attach best geometry to the exception so the outer ladder can
                # use it after the local frame is unwound.
                plateau_exc.best_xyz = plateau_xyz  # type: ignore[attr-defined]
                plateau_exc.best_fmax_ev_ang = calc.best_fmax_ev_ang  # type: ignore[attr-defined]
                # Re-raise so the outer ladder breaks out cleanly to single-point.
                raise
            except (ImportError, AttributeError, Exception) as e:
                # On stall/timeout, update mol with Sella's best geometry
                # so the geomeTRIC fallback starts from partial progress.
                if 'TSOptimizationStall' in type(e).__name__:
                    logger.warning(f"    [SELLA] {e}. Recovering partial geometry for geomeTRIC fallback...")
                    # P3: surface escalation hint to outer ladder via self.
                    if getattr(calc, 'escalate_to_level_shift', False):
                        self._ladder_escalate_to_level_shift = True
                    try:
                        from pyscf.data import nist as _nist
                        mol.set_geom_(ase_atoms.positions / _nist.BOHR, unit='Bohr')
                    except Exception:
                        pass
                else:
                    logger.info(f"    WARNING: Sella/ASE bridge failed ({type(e).__name__}: {e}). Falling back to geomeTRIC...")

            is_ts_fallback = is_ts and not conv
        else:
            is_ts_fallback = is_ts

        if not is_ts or (is_ts and not self.ts_optimizer) or (is_ts and not conv):
            conv, mol_opt, fallback_xyz = self._optimize_with_geometric_kernel(
                mf,
                label=label,
                max_steps=max_steps,
                is_ts_fallback=is_ts_fallback,
                eff_use_explicit=eff_use_explicit,
                n_atoms_solute=n_atoms_solute,
            )
            if not conv and fallback_xyz is None:
                logger.info("    [GEOMETRIC] Retrying failed TRIC optimization in Cartesian coordinates...")
                conv, mol_opt, fallback_xyz = self._optimize_with_geometric_kernel(
                    mf,
                    label=f"{label}_cart",
                    max_steps=max_steps,
                    is_ts_fallback=is_ts_fallback,
                    eff_use_explicit=eff_use_explicit,
                    n_atoms_solute=n_atoms_solute,
                    coords="cart",
                )
            if fallback_xyz:
                opt_xyz = fallback_xyz
            else:
                opt_xyz = _molecule_to_xyz(mol_opt, comment="optimized geometry")

        return conv, mol_opt, opt_xyz

    def single_point(self, xyz_content: str, xc_method: str = 'wB97M-V', basis: str = 'def2-tzvp', charge: int = 0, spin: int = 0) -> float:
        """Run a high-level single-point energy calculation."""
        mol = self._setup_mol(xyz_content, charge, spin, basis=basis)
        
        # PySCF supports wB97M-V via libxc
        # Since libxc support varies, we map to closest available or pass string directly.
        mf = self._build_mf(mol, xc_method=xc_method, use_solvent=True)
        try:
            energy = mf.scf()
        except Exception as exc:
            logger.warning(f"    [SP] Primary single-point crashed: {exc}. Retrying without solvent...")
            mf = self._build_mf(mol, xc_method=xc_method, use_solvent=False, harden_scf=True)
            energy = mf.scf()
        
        if not mf.converged:
            logger.info(f"    WARNING: Single-point SCF did not converge. Retrying with hardened SCF parameters...")
            mf_hard = self._build_mf(mol, xc_method=xc_method, use_solvent=True, harden_scf=True)
            try:
                energy = mf_hard.scf()
            except Exception as exc:
                logger.warning(f"    [SP] Hardened single-point crashed: {exc}. Falling back to vacuum.")
                mf_vac = self._build_mf(mol, xc_method=xc_method, use_solvent=False, harden_scf=True)
                energy = mf_vac.scf()
                mf_hard = mf_vac
            if not mf_hard.converged:
                logger.warning(f"    [CRITICAL] Single-point SCF still failed to converge! Returning unconverged energy as fallback.")
        
        # [SOTA GUARD] Spin Integrity Check
        self._log_spin_integrity(mf)

        return energy
        
    def _run_hessian_and_thermo(self, mf) -> Tuple[List[float], float, float]:
        """
        Run frequency analysis and compute thermochemical corrections. 
        Returns: (frequencies, raw_gibbs_free_energy, quasi_harmonic_gibbs_free_energy)
        """
        natm = mf.mol.natm

        # For molecules with > 20 heavy atoms the full DFT Hessian is
        # prohibitively expensive (∼ 6N gradient evaluations).  Use xTB
        # frequencies for the thermal/entropic correction while keeping the
        # DFT electronic energy.  This mirrors the ORCA --hess-xtb workflow.
        if natm > 20:
            xtb_result = self._run_xtb_hessian_thermo(mf)
            if xtb_result is not None:
                return xtb_result

        try:
            # We must use a vacuum Hessian calculation if the solvent is ddCOSMO because
            # PySCF 2.x Hessian with ddCOSMO has known issues.  Solvent effects on
            # frequencies are typically minor compared to the single-point energy.
            if hasattr(mf, 'with_solvent'):
                xc = getattr(mf, 'xc', None)
                mf_vac = self._build_mf(
                    mf.mol,
                    xc_method=xc or 'hf',
                    use_solvent=False,
                    conv_tol=mf.conv_tol,
                )
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
            thermo_info = thermo.thermo(mf_to_use, freq_info['freq_au'], self.temp_k, 101325)
            
            # Apply quasi-harmonic
            corrector = QuasiHarmonicCorrector(temp_k=self.temp_k)
            # Use the ORIGINAL mf energy (which might have solvent)
            # Note: G_tot already contains E_tot. We need to be careful with the delta.
            # G_qh = E_elec + thermal_qh
            qh_result = corrector.calculate_thermo(freqs, electronic_energy_h=mf.e_tot)
            
            # [SOTA GUARD] Entropy Ceiling
            thermal_corr_kcal = (qh_result.qh_gibbs_h - mf.e_tot) * 627.509
            if abs(thermal_corr_kcal) > 40.0:
                logger.warning(f"    [SCIENTIFIC CAUTION] Unusually large thermal correction: {thermal_corr_kcal:.1f} kcal/mol. "
                               "Check for exploding low-frequency modes.")

            return freqs.tolist(), thermo_info['G_tot'][0] - mf_to_use.e_tot + mf.e_tot, qh_result.qh_gibbs_h
            
        except Exception as exc:
            logger.warning(f"    [WARNING] Hessian/Frequency analysis failed: {exc}")
            logger.info("    Falling back to electronic energy only (zero vibrational correction).")
            # Return electronic energy as a surrogate for Gibbs in the worst-case
            # to avoid throwing away a valid optimization and single-point.
            # Mark with empty freq list so callers can detect missing thermo.
            logger.warning("    [SCIENTIFIC CAUTION] Barrier will lack thermal/entropic corrections. "
                          "Treat as electronic-only estimate.")
            return [], mf.e_tot, mf.e_tot

    def _run_xtb_hessian_thermo(self, mf) -> Optional[Tuple[List[float], float, float]]:
        """Delegate to :mod:`xtb_backend`."""
        mol = mf.mol
        xyz_text = _molecule_to_xyz(mol, comment="xtb hessian input")
        return _xtb.run_xtb_hessian_thermo(
            xyz_text, mol.charge, mol.spin, mf.e_tot, self.temp_k,
        )

    @staticmethod
    def _parse_xtb_frequencies(stdout: str) -> List[float]:
        """Delegate to :mod:`xtb_backend`."""
        return _xtb.parse_xtb_frequencies(stdout)

    def optimize_geometry(self, xyz_content: str, charge: int=0, spin: int=0, is_ts: bool=False, max_steps: int=100, label: str = "default", use_explicit_solvent: Optional[bool] = None, n_water: Optional[int] = None, progress_callback: Optional[Callable[[str, Dict[str, Any]], None]] = None, ts_refine_ctx: Optional[Dict[str, Any]] = None) -> DFTResult:
        """Run geometry/TS optimization using the SOTA Stabilized Ladder:
        Stage 1: MLP (MACE-OFF) -> Stage 2: xTB (Semi-empirical) -> Stage 2b: P0 TS-guess refine -> Stage 3: DFT-Warmup (Level-Shifted) -> Stage 4: DFT-Final.

        Parameters
        ----------
        ts_refine_ctx : dict, optional
            Context for P0 TS-guess refinement.  Keys:
            ``reactant_xyz``, ``product_xyz``, ``path_frames_dir``,
            ``policy`` (auto/always/off), ``threshold_fmax`` (float).
        """
        phase_prefix = "ts" if is_ts else "reactant"
        phase_label = "Transition State" if is_ts else "Reactant"
        ckpt_ctx = self._get_optimization_checkpoint_context(phase_prefix)
        phase_ckpt_dir = ckpt_ctx["dir"] if ckpt_ctx else None
        resume_state = ckpt_ctx.get("resume_state", {}) if ckpt_ctx else {}
        progress_state: Dict[str, Any] = dict(resume_state)
        progress_path = phase_ckpt_dir / f"{phase_prefix}_progress.json" if phase_ckpt_dir else None

        def _save_phase_xyz(file_name: str, xyz: Optional[str]) -> None:
            if phase_ckpt_dir is None or not xyz:
                return
            try:
                phase_ckpt_dir.mkdir(parents=True, exist_ok=True)
                (phase_ckpt_dir / file_name).write_text(xyz, encoding="utf-8")
                logger.info(f"    [CHECKPOINT] Saved geometry {phase_ckpt_dir / file_name}")
            except Exception as exc:
                logger.warning(f"    [CHECKPOINT] WARNING: Failed to save geometry {file_name}: {exc}")

        def _save_phase_progress(*, xyz: Optional[str] = None, **updates: Any) -> None:
            if phase_ckpt_dir is None or progress_path is None:
                return
            try:
                phase_ckpt_dir.mkdir(parents=True, exist_ok=True)
                if xyz:
                    inflight_name = f"{phase_prefix}_inflight.xyz"
                    _save_phase_xyz(inflight_name, xyz)
                    progress_state["current_xyz_file"] = inflight_name
                progress_state.update(updates)
                progress_state["phase_prefix"] = phase_prefix
                progress_state["phase_label"] = phase_label
                progress_state["target_label"] = label
                progress_state["updated_at"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
                progress_path.write_text(json.dumps(progress_state, indent=2), encoding="utf-8")
            except Exception as exc:
                logger.warning(f"    [CHECKPOINT] WARNING: Failed to save progress {progress_path}: {exc}")

        def _save_named_geometry(stage_name: str, xyz: Optional[str]) -> None:
            if not xyz:
                return
            file_name = f"{phase_prefix}_{stage_name}.xyz"
            _save_phase_xyz(file_name, xyz)
            _save_phase_progress(xyz=xyz, last_stage=stage_name, last_geometry_file=file_name)

        resume_geometry_ready = bool(progress_state.get("geometry_ready_for_postprocessing", False))
        resume_mlp_completed = bool(progress_state.get("mlp_completed", False))
        resume_xtb_completed = bool(progress_state.get("xtb_completed", False))
        try:
            resume_strategy_index = int(progress_state.get("next_strategy_index", 0) or 0)
        except (TypeError, ValueError):
            resume_strategy_index = 0
        resume_strategy_index = max(0, resume_strategy_index)
        if progress_state and max_steps > 0:
            logger.info(
                f"    [RESUME] Found partial {phase_label.lower()} checkpoint "
                f"(geometry_ready={resume_geometry_ready}, next_strategy_index={resume_strategy_index})."
            )
        
        # [SOTA GUARD] Check if we can resume from a persistent trajectory
        workspace_dir = Path("data/geometries/dft_runs") / f"{label}_{phase_prefix}"
        traj_path = workspace_dir / "optim.xyz"
        if traj_path.exists() and max_steps > 0 and not progress_state:
            logger.info(f"    [RESUME] Found persistent trajectory in {traj_path}. Checking for usable coordinates...")
            try:
                with open(traj_path, "r") as f:
                    lines = f.readlines()
                if lines:
                    natoms = int(lines[0].strip())
                    if len(lines) >= natoms + 2:
                        xyz_content = "".join(lines[-(natoms+2):])
                        logger.info(f"    [RESUME] Successfully recovered last frame from persistent trajectory.")
            except Exception as e:
                logger.warning(f"    [RESUME] Failed to parse trajectory: {e}")

        current_xyz = xyz_content
        # Determine effective solvation settings from call-time overrides first,
        # then fall back to the instance defaults configured on the refiner.
        effective_n_water = self.n_water if n_water is None else n_water
        if use_explicit_solvent is None:
            eff_use_explicit = bool(self.use_explicit_solvent or (n_water is not None and n_water > 0))
        else:
            eff_use_explicit = bool(use_explicit_solvent)
        if effective_n_water is None:
            effective_n_water = 0
        if effective_n_water <= 0:
            eff_use_explicit = False

        n_atoms_solute = len(_extract_atom_string(xyz_content).splitlines()) if eff_use_explicit else 0
        current_atom_count = len(_extract_atom_string(current_xyz).splitlines())
        if eff_use_explicit and current_atom_count <= n_atoms_solute:
            logger.info(
                "    [SOLVATION] Generating explicit solvent cluster with %d water molecules.",
                effective_n_water,
            )
            current_xyz = self.solvation_engine.generate_solvated_cluster(
                current_xyz,
                n_water=effective_n_water,
                freeze_core=is_ts,
            )

        conv = False
        mol_opt = None
        opt_xyz = None
        mlp_produced_update = False

        # [LADDER] Stage 1: MLP (MACE-OFF)
        if max_steps > 0 and not resume_geometry_ready:
            self._emit_progress(progress_callback, f"{phase_prefix}_prerelaxation", label=f"{phase_label} MLP pre-relaxation")
            if resume_mlp_completed:
                logger.info(f"    [RESUME] Skipping {phase_label.lower()} MLP pre-relaxation (progress checkpoint).")
            else:
                pre_mlp_xyz = current_xyz
                current_xyz = self._mlp_prerelax(current_xyz)
                mlp_produced_update = current_xyz.strip() != pre_mlp_xyz.strip()
                _save_named_geometry("mlp_prerelaxed", current_xyz)
                _save_phase_progress(
                    xyz=current_xyz,
                    current_stage="mlp_prerelaxation",
                    mlp_completed=True,
                    xtb_completed=resume_xtb_completed,
                    next_strategy_index=resume_strategy_index,
                    geometry_ready_for_postprocessing=False,
                )
                resume_mlp_completed = True

        # [LADDER] Stage 2 (TS-safe): xTB crude pre-relax for TS guesses when MLPs reverted.
        # All MLPs returned the input unchanged (drift rejection). The TS guess may have
        # steric clashes that crash SCF in the DFT ladder. Run a minimal xTB GFN2 pass
        # (crude convergence, ~20 cycles, gas-phase) to relieve clashes without driving
        # the geometry into a nearby minimum (which would destroy the saddle character).
        if (
            is_ts
            and max_steps > 0
            and not resume_geometry_ready
            and not resume_xtb_completed
            and not mlp_produced_update
        ):
            try:
                relaxed_xyz = self._run_xtb_ts_safe_prerelax(current_xyz, charge, spin)
                current_xyz = relaxed_xyz
                _save_named_geometry("xtb_ts_safe_prerelaxed", current_xyz)
                _save_phase_progress(
                    xyz=current_xyz,
                    current_stage="xtb_ts_safe_prerelaxation",
                    mlp_completed=resume_mlp_completed,
                    xtb_completed=True,
                    next_strategy_index=resume_strategy_index,
                    geometry_ready_for_postprocessing=False,
                )
                resume_xtb_completed = True
                logger.info("    [LADDER] xTB TS-safe pre-relax produced an updated geometry.")
            except RuntimeError as exc:
                logger.warning(
                    f"    [LADDER] xTB TS-safe pre-relax failed ({exc}). "
                    "Continuing with MLP/input geometry into the DFT ladder."
                )

        # [P0] TS-guess quality gate & refinement
        # Probes fmax and n_imag at xTB level; if the guess is poor, attempts
        # relaxed scan and CI-NEB at xTB to produce a better starting geometry
        # before entering the expensive DFT ladder.
        ts_refinement_result: Optional[Dict[str, Any]] = None
        if (
            is_ts
            and max_steps > 0
            and not resume_geometry_ready
            and ts_refine_ctx is not None
            and ts_refine_ctx.get("policy", "auto") != "off"
        ):
            from .ts_guess_refiner import refine_if_needed, RefinementOutcome
            try:
                outcome: RefinementOutcome = refine_if_needed(
                    ts_xyz=current_xyz,
                    reactant_xyz=ts_refine_ctx.get("reactant_xyz", ""),
                    product_xyz=ts_refine_ctx.get("product_xyz"),
                    charge=charge,
                    spin=spin,
                    threshold_fmax=float(ts_refine_ctx.get("threshold_fmax", 1.0)),
                    path_frames_dir=ts_refine_ctx.get("path_frames_dir"),
                    policy=str(ts_refine_ctx.get("policy", "auto")),
                    enable_gsm_fallback=bool(ts_refine_ctx.get("enable_gsm_fallback", False)),
                    gsm_work_dir=str(phase_ckpt_dir / "gsm_scratch") if phase_ckpt_dir else None,
                    gsm_nodes=int(ts_refine_ctx.get("gsm_nodes", 11)),
                    gsm_max_iters=int(ts_refine_ctx.get("gsm_max_iters", 80)),
                    gsm_timeout_s=float(ts_refine_ctx.get("gsm_timeout_s", 1800.0)),
                )
                # Persist audit trail
                ts_refinement_result = {
                    "initial_probe": outcome.initial_probe,
                    "policy": str(ts_refine_ctx.get("policy", "auto")),
                    "threshold_fmax_ev_ang": float(ts_refine_ctx.get("threshold_fmax", 1.0)),
                    "triggered": outcome.final_method != "original",
                    "attempts": [
                        {"method": a.method, "fmax_ev_ang": a.fmax_ev_ang,
                         "n_imag": a.n_imag, "accepted": a.accepted, "reason": a.reason}
                        for a in outcome.attempts
                    ],
                    "final_method": outcome.final_method,
                    "decision": "accepted" if outcome.accepted else "rejected",
                    "reason": outcome.reason,
                }
                if phase_ckpt_dir:
                    try:
                        (phase_ckpt_dir / "ts_xtb_probe_initial.json").write_text(
                            json.dumps(outcome.initial_probe, indent=2), encoding="utf-8"
                        )
                        (phase_ckpt_dir / "ts_xtb_refinement.json").write_text(
                            json.dumps(ts_refinement_result, indent=2), encoding="utf-8"
                        )
                    except Exception as ckpt_exc:
                        logger.warning(f"    [P0] Checkpoint write failed: {ckpt_exc}")

                if not outcome.accepted:
                    raise BadTSGuessRejected(outcome.reason)

                if outcome.final_method != "original":
                    logger.info(
                        f"    [P0] Replacing TS guess: {outcome.reason}"
                    )
                    current_xyz = outcome.xyz
                    _save_named_geometry("xtb_refined", current_xyz)
                    if phase_ckpt_dir:
                        try:
                            (phase_ckpt_dir / "ts_xtb_refined.xyz").write_text(
                                current_xyz, encoding="utf-8"
                            )
                        except Exception:
                            pass
            except BadTSGuessRejected:
                raise
            except Exception as exc:
                logger.warning(f"    [P0] TS refinement pipeline failed ({exc}); continuing with current guess.")

        # [LADDER] Stage 2: xTB (GFN2-ALPB)
        # xTB is much better at capturing chemical repulsive forces than MLPs for complex crosslinks.
        if not is_ts and max_steps > 0 and not resume_geometry_ready:
            try:
                if resume_xtb_completed:
                    logger.info("    [RESUME] Skipping reactant xTB optimization (progress checkpoint).")
                else:
                    current_xyz = self._run_xtb_optimization(current_xyz, charge, spin)
                    _save_named_geometry("xtb_optimized", current_xyz)
                    _save_phase_progress(
                        xyz=current_xyz,
                        current_stage="xtb_optimization",
                        mlp_completed=resume_mlp_completed,
                        xtb_completed=True,
                        next_strategy_index=resume_strategy_index,
                        geometry_ready_for_postprocessing=False,
                    )
                    resume_xtb_completed = True
            except RuntimeError as exc:
                logger.warning(f"    [LADDER] xTB stage failed ({exc}). Continuing with MLP/input geometry.")
        
        # [P0b] Detect the forming-bond pair so the DFT ladder can perform
        # a constrained pre-relaxation that drains backbone forces before Sella.
        forming_bond_pair: Optional[Tuple[int, int]] = None
        if is_ts and ts_refine_ctx is not None and ts_refine_ctx.get("product_xyz"):
            try:
                from .xtb_backend import detect_forming_bond as _detect_fb
                fb_result = _detect_fb(
                    ts_refine_ctx.get("reactant_xyz", ""),
                    ts_refine_ctx["product_xyz"],
                )
                if fb_result is not None:
                    forming_bond_pair = fb_result
                    logger.info(
                        f"    [P0b] Forming bond for DFT prerelax: "
                        f"atoms {forming_bond_pair}"
                    )
            except Exception as fb_exc:
                logger.warning(
                    f"    [P0b] Could not detect forming bond: {fb_exc}"
                )

        # [P0c] Auto-load pre-computed TS mode (v0) from checkpoint when
        # available.  This is produced by compute_xtb_ts_mode() during a
        # previous GSM fallback run and allows Sella to follow the correct
        # imaginary mode from the very first DFT step.
        if (
            is_ts
            and ts_refine_ctx is not None
            and "ts_v0" not in ts_refine_ctx
            and phase_ckpt_dir is not None
        ):
            _v0_path = phase_ckpt_dir / "ts_v0_xtb.npy"
            if _v0_path.exists():
                try:
                    ts_refine_ctx["ts_v0"] = np.load(str(_v0_path))
                    logger.info(
                        "    [P0c] Loaded pre-computed TS mode from %s",
                        _v0_path.name,
                    )
                except Exception as _v0_exc:
                    logger.warning(
                        "    [P0c] Failed to load ts_v0_xtb.npy: %s", _v0_exc
                    )

        # [LADDER] Stage 3: DFT Geometry Optimization with SCF Strategy Escalation
        #
        # Instead of a monolithic warm-up + refinement split, we try progressively
        # harder SCF strategies.  Each strategy can fail in two ways:
        #   (a) SCF gradient crash  →  recover best geometry, try next strategy
        #   (b) geomeTRIC doesn't converge  →  pass best geometry to next strategy
        # If any strategy converges, we stop.
        _SCF_STRATEGIES = [
            {"name": "standard-DIIS",
             "use_soscf": False, "smearing_sigma": None, "level_shift": 0.0,
             "harden_scf": False, "use_solvent": not is_ts},
            {"name": "hardened-DIIS+smearing",
             "use_soscf": False, "smearing_sigma": 0.005, "level_shift": 0.0,
             "harden_scf": True, "use_solvent": not is_ts},
        ]
        if is_ts:
            # SOTA plateau-breaking ladder for TS searches.
            #
            # Smearing+Newton is *not* included: PySCF's analytical derivatives
            # (needed by Sella/geomeTRIC) break under fractional occupations —
            # this surfaced as `NumPy boolean array indexing assignment` crashes
            # in the lysinoalanine TS run. Instead, we break inverted HOMO-LUMO
            # plateaus with level_shift (raises virtual orbital energies and
            # directly heals the aufbau-violating state) *before* running Newton.
            _SCF_STRATEGIES.extend([
                {"name": "level-shifted-DIIS",
                 "use_soscf": False, "smearing_sigma": None, "level_shift": 0.5,
                 "harden_scf": False, "use_solvent": False},
                {"name": "SOSCF/Newton+level-shift",
                 "use_soscf": True, "smearing_sigma": None, "level_shift": 0.3,
                 "harden_scf": False, "use_solvent": False},
            ])
        else:
            _SCF_STRATEGIES.extend([
                {"name": "SOSCF/Newton",
                 "use_soscf": True, "smearing_sigma": None, "level_shift": 0.0,
                 "harden_scf": False, "use_solvent": True},
                {"name": "vacuum-SOSCF",
                 "use_soscf": True, "smearing_sigma": None, "level_shift": 0.0,
                 "harden_scf": True, "use_solvent": False},
            ])

        if max_steps > 0:
            self._emit_progress(progress_callback, f"{phase_prefix}_geometry_optimization", label=f"{phase_label} geometry optimization")
            last_exc: Optional[Exception] = None
            if resume_geometry_ready:
                logger.info(
                    f"    [RESUME] Found geometry-ready checkpoint for {phase_label.lower()}. "
                    "Skipping geometry stages and resuming post-processing."
                )
                mol_opt = self._setup_mol(current_xyz, charge, spin, basis=self.opt_basis)
                conv = bool(progress_state.get("converged", False))
                opt_xyz = current_xyz
            else:
                start_strategy_index = min(resume_strategy_index, len(_SCF_STRATEGIES))
                if start_strategy_index > 0:
                    logger.info(
                        f"    [RESUME] Resuming ladder from strategy {start_strategy_index + 1}/{len(_SCF_STRATEGIES)}."
                    )

                # P3: index of the first level-shifted strategy in this phase's
                # ladder, used to skip ahead when an inverted HOMO-LUMO gap is
                # detected.  None if no level-shifted strategy exists.
                _level_shift_idx = next(
                    (i for i, s in enumerate(_SCF_STRATEGIES)
                     if s.get("level_shift", 0.0) > 0.0 and not s.get("use_soscf", False)),
                    None,
                )
                self._ladder_escalate_to_level_shift = False

                strat_idx = start_strategy_index
                while strat_idx < len(_SCF_STRATEGIES):
                    strat = _SCF_STRATEGIES[strat_idx]
                    strat_name = strat["name"]
                    strategy_slug = self._checkpoint_slug(strat_name)
                    _save_phase_progress(
                        xyz=current_xyz,
                        current_stage="dft_geometry_optimization",
                        mlp_completed=resume_mlp_completed,
                        xtb_completed=resume_xtb_completed,
                        current_strategy_index=strat_idx,
                        current_strategy_name=strat_name,
                        next_strategy_index=strat_idx,
                        geometry_ready_for_postprocessing=False,
                    )
                    logger.info(f"    [LADDER] Stage 3: DFT optimization — strategy '{strat_name}'")
                    try:
                        conv, mol_opt, opt_xyz = self._optimize_with_pyscf_backend(
                            current_xyz,
                            charge,
                            spin,
                            is_ts,
                            max_steps,
                            eff_use_explicit,
                            n_atoms_solute,
                            label=label,
                            use_soscf=strat["use_soscf"],
                            smearing_sigma=strat["smearing_sigma"],
                            harden_scf=strat["harden_scf"],
                            use_solvent_optimization=strat["use_solvent"],
                            level_shift=strat.get("level_shift", 0.0),
                            forming_bond_pair=forming_bond_pair,
                            ts_v0=ts_refine_ctx.get("ts_v0") if ts_refine_ctx else None,
                        )
                        if opt_xyz:
                            current_xyz = opt_xyz
                            _save_named_geometry(f"ladder_{strat_idx + 1:02d}_{strategy_slug}", current_xyz)
                        if conv:
                            logger.info(f"    [LADDER] Converged with strategy '{strat_name}'.")
                            _save_phase_progress(
                                xyz=current_xyz,
                                current_stage="geometry_complete",
                                current_strategy_index=strat_idx,
                                current_strategy_name=strat_name,
                                last_strategy_index=strat_idx,
                                last_strategy_name=strat_name,
                                last_strategy_status="converged",
                                next_strategy_index=strat_idx,
                                geometry_ready_for_postprocessing=True,
                                converged=True,
                                mlp_completed=resume_mlp_completed,
                                xtb_completed=resume_xtb_completed,
                            )
                            break
                        # Partial progress — feed best geometry to next strategy
                        if opt_xyz:
                            logger.info(f"    [LADDER] Strategy '{strat_name}' did not converge. Using partial geometry for next attempt.")
                        # P3: if Sella surfaced an inverted/near-degenerate
                        # HOMO-LUMO gap, jump straight to the first level-shifted
                        # strategy instead of crawling through smearing first.
                        next_idx = strat_idx + 1
                        if (
                            self._ladder_escalate_to_level_shift
                            and _level_shift_idx is not None
                            and _level_shift_idx > strat_idx
                        ):
                            logger.info(
                                f"    [LADDER] Orbital instability detected — "
                                f"skipping ahead to '{_SCF_STRATEGIES[_level_shift_idx]['name']}'."
                            )
                            next_idx = _level_shift_idx
                            self._ladder_escalate_to_level_shift = False
                        _save_phase_progress(
                            xyz=current_xyz,
                            current_stage="dft_geometry_optimization",
                            current_strategy_index=strat_idx,
                            current_strategy_name=strat_name,
                            last_strategy_index=strat_idx,
                            last_strategy_name=strat_name,
                            last_strategy_status="partial",
                            next_strategy_index=next_idx,
                            geometry_ready_for_postprocessing=False,
                            converged=False,
                            mlp_completed=resume_mlp_completed,
                            xtb_completed=resume_xtb_completed,
                        )
                        strat_idx = next_idx
                        continue
                    except Exception as exc:
                        last_exc = exc
                        # Plateau detected — best geometry was attached to the
                        # exception by _optimize_with_pyscf_backend.  Skip the
                        # remaining ladder strategies and proceed to single-point.
                        if 'TSPlateauConverged' in type(exc).__name__:
                            best_xyz = getattr(exc, 'best_xyz', None)
                            best_fmax = getattr(exc, 'best_fmax_ev_ang', None)
                            logger.warning(
                                f"    [LADDER] Strategy '{strat_name}' detected a "
                                f"numerical plateau (best fmax="
                                f"{best_fmax if best_fmax is not None else 'n/a'}). "
                                "Stopping ladder and using best geometry for "
                                "single-point + Hessian."
                            )
                            if best_xyz:
                                opt_xyz = best_xyz
                                current_xyz = opt_xyz
                                mol_opt = self._setup_mol(opt_xyz, charge, spin, basis=self.opt_basis)
                            _save_named_geometry("sella_plateau_best", current_xyz)
                            _save_phase_progress(
                                xyz=current_xyz,
                                current_stage="geometry_complete",
                                current_strategy_index=strat_idx,
                                current_strategy_name=strat_name,
                                last_strategy_index=strat_idx,
                                last_strategy_name=strat_name,
                                last_strategy_status="plateau_aborted",
                                last_error=str(exc),
                                next_strategy_index=len(_SCF_STRATEGIES),
                                geometry_ready_for_postprocessing=True,
                                converged=False,
                                mlp_completed=resume_mlp_completed,
                                xtb_completed=resume_xtb_completed,
                            )
                            conv = False
                            break
                        if self._is_scf_gradient_failure(exc):
                            logger.warning(f"    [LADDER] Strategy '{strat_name}' hit SCF gradient failure: {exc}")
                            _save_phase_progress(
                                xyz=current_xyz,
                                current_stage="dft_geometry_optimization",
                                current_strategy_index=strat_idx,
                                current_strategy_name=strat_name,
                                last_strategy_index=strat_idx,
                                last_strategy_name=strat_name,
                                last_strategy_status="scf_gradient_failure",
                                last_error=str(exc),
                                next_strategy_index=strat_idx + 1,
                                geometry_ready_for_postprocessing=False,
                                converged=False,
                                mlp_completed=resume_mlp_completed,
                                xtb_completed=resume_xtb_completed,
                            )
                            strat_idx += 1
                            continue
                        if not is_ts:
                            logger.warning(
                                f"    [LADDER] Strategy '{strat_name}' raised a non-recoverable "
                                f"reactant optimisation error: {exc}. Proceeding with best available geometry."
                            )
                            conv, mol_opt, opt_xyz = self._best_available_geometry_fallback(
                                current_xyz,
                                charge=charge,
                                spin=spin,
                                basis=self.opt_basis,
                                label=f"{phase_label} refinement",
                            )
                            if opt_xyz:
                                current_xyz = opt_xyz
                            _save_named_geometry("best_available", current_xyz)
                            _save_phase_progress(
                                xyz=current_xyz,
                                current_stage="geometry_complete",
                                current_strategy_index=strat_idx,
                                current_strategy_name=strat_name,
                                last_strategy_index=strat_idx,
                                last_strategy_name=strat_name,
                                last_strategy_status="exception_fallback",
                                last_error=str(exc),
                                next_strategy_index=len(_SCF_STRATEGIES),
                                geometry_ready_for_postprocessing=True,
                                converged=False,
                                mlp_completed=resume_mlp_completed,
                                xtb_completed=resume_xtb_completed,
                            )
                            break
                        _save_phase_progress(
                            xyz=current_xyz,
                            current_stage="dft_geometry_optimization",
                            current_strategy_index=strat_idx,
                            current_strategy_name=strat_name,
                            last_strategy_index=strat_idx,
                            last_strategy_name=strat_name,
                            last_strategy_status="exception",
                            last_error=str(exc),
                            next_strategy_index=strat_idx,
                            geometry_ready_for_postprocessing=False,
                            converged=False,
                            mlp_completed=resume_mlp_completed,
                            xtb_completed=resume_xtb_completed,
                        )
                        logger.info(f"    [LADDER ERROR] Refinement Ladder failed: {exc}")
                        raise
                    # End of try/except — advance to next strategy.
                    strat_idx += 1
                else:
                    # All strategies exhausted without convergence (while-else)
                    if mol_opt is None and last_exc is not None:
                        logger.info(f"    [LADDER ERROR] All SCF strategies exhausted. Last error: {last_exc}")
                        conv, mol_opt, opt_xyz = self._best_available_geometry_fallback(
                            current_xyz,
                            charge=charge,
                            spin=spin,
                            basis=self.opt_basis,
                            label=f"{'TS' if is_ts else 'Reactant'} refinement",
                        )
                        if opt_xyz:
                            current_xyz = opt_xyz
                    if current_xyz:
                        _save_named_geometry("best_available", current_xyz)
                    _save_phase_progress(
                        xyz=current_xyz,
                        current_stage="geometry_complete",
                        last_strategy_index=len(_SCF_STRATEGIES) - 1,
                        last_strategy_name=_SCF_STRATEGIES[-1]["name"],
                        last_strategy_status="exhausted",
                        next_strategy_index=len(_SCF_STRATEGIES),
                        ladder_exhausted=True,
                        geometry_ready_for_postprocessing=True,
                        converged=False,
                        mlp_completed=resume_mlp_completed,
                        xtb_completed=resume_xtb_completed,
                    )
                    logger.warning("    [LADDER] No strategy converged. Proceeding with best available geometry.")
        else:
            # Checkpoint resume shortcut
            mol_opt = self._setup_mol(current_xyz, charge, spin, basis=self.opt_basis)
            conv = True
            opt_xyz = current_xyz

        # Stage 5: Post-Refinement Single Point & Properties
        if not conv:
            logger.warning(f"    [LADDER WARNING] Final stage did not converge. Using best guess.")
        
        # Use final geometry for high-level single-point and thermo
        final_xyz = opt_xyz if opt_xyz else _molecule_to_xyz(mol_opt)
        _save_named_geometry("geometry_ready", final_xyz)
        _save_phase_progress(
            xyz=final_xyz,
            current_stage="single_point",
            geometry_ready_for_postprocessing=True,
            converged=conv,
            mlp_completed=resume_mlp_completed,
            xtb_completed=resume_xtb_completed,
        )
        
        self._emit_progress(progress_callback, f"{phase_prefix}_single_point",
                            label=f"{phase_label} high-level single-point ({self.refinement_method}/{self.refinement_basis})")
        logger.info(f">>> [Phase 3.3] Running high-level single-point: {self.refinement_method}/{self.refinement_basis}")
        electronic_energy = self.single_point(
            final_xyz,
            xc_method=self.refinement_method,
            basis=self.refinement_basis,
            charge=charge,
            spin=spin
        )
        
        self._emit_progress(progress_callback, f"{phase_prefix}_hessian",
                            label=f"{phase_label} frequency analysis (natoms={getattr(mol_opt, 'natm', '?')})")
        logger.info(f">>> [Phase 3.3] Running frequency analysis & quasi-harmonic thermo...")
        _save_phase_progress(
            xyz=final_xyz,
            current_stage="hessian",
            geometry_ready_for_postprocessing=True,
            converged=conv,
            mlp_completed=resume_mlp_completed,
            xtb_completed=resume_xtb_completed,
        )
        mf_final = self._build_mf(mol_opt, xc_method=self.opt_method, use_solvent=False)
        mf_final.e_tot = electronic_energy # Inject SP energy for thermo output matching
        
        freqs, raw_g, qh_g = self._run_hessian_and_thermo(mf_final)
        _save_phase_progress(
            xyz=final_xyz,
            current_stage="completed",
            geometry_ready_for_postprocessing=True,
            converged=conv,
            completed=True,
            mlp_completed=resume_mlp_completed,
            xtb_completed=resume_xtb_completed,
        )
        
        return DFTResult(
            method=f"{self.refinement_method}/{self.refinement_basis}",
            optimized_xyz=final_xyz,
            energy_hartree=electronic_energy,
            gibbs_free_energy_hartree=raw_g,
            quasi_harmonic_gibbs_hartree=qh_g,
            frequencies_cm1=freqs,
            converged=conv
        )
        
    def generate_irc(self, ts_xyz: str, charge: int=0, spin: int=0, step_size: float=0.05) -> Tuple[str, str]:
        """
        Intrinsic Reaction Coordinate (IRC) validation via Displacement + Optimization.
        
        1. Calculates Hessian at the TS.
        2. Identifying the imaginary mode.
        3. Displaces +/- along the mode.
        4. Optimizes both endpoints to find the connected basins.
        """
        logger.info(">>> [Phase 3.4] Starting IRC Validation (Double-Optimization method)...")
        mol = self._setup_mol(ts_xyz, charge, spin, basis=self.opt_basis)
        mf = self._build_mf(mol, xc_method=self.opt_method, use_solvent=False)
        mf.conv_tol = 1e-8
        mf.scf()
        
        logger.info("    Computing TS Hessian for mode identification...")
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
        
        logger.info(f"    Lowest frequencies (cm^-1): {freqs[:6]}")
        
        # We look for a significant imaginary frequency (conventionally negative in PySCF)
        # We use a threshold of -50 cm^-1 to avoid numerical noise
        if freqs[0] >= -50.0:
            raise ValueError(f"No significant imaginary frequency found at the provided TS geometry (Lowest freq: {freqs[0]:.1f} cm^-1)")
        
        logger.info(f"    Found imaginary mode with frequency: {freqs[0]:.1f} cm^-1")
        
        # The first mode in h_info['norm_mode'] is the lowest frequency mode
        mode_vec = h_info['norm_mode'][0]
        
        orig_coords = mol.atom_coords() # in Bohr
        
        endpoints = []
        for direction in [1.0, -1.0]:
            label = "Forward" if direction > 0 else "Backward"
            logger.info(f"    Following IRC {label} path...")
            
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
                    endpoints.append(_molecule_to_xyz(opt_mol, comment=f"irc endpoint {label.lower()}"))
                finally:
                    os.chdir(pwd)
        
        return endpoints[0], endpoints[1]

    # ------------------------------------------------------------------
    # Robust barrier estimation: xTB geometry + DFT single-point
    # ------------------------------------------------------------------
    def calculate_robust_barrier(
        self,
        reactant_xyz: str,
        ts_guess_xyz: str,
        charge: int = 0,
        spin: int = 0,
        checkpoint_dir: Optional[str] = None,
        progress_callback: Optional[Callable[[str, Dict[str, Any]], None]] = None,
        product_xyz: Optional[str] = None,
        reaction_meta: Optional[Dict] = None,
        enable_gsm_fallback: bool = False,
        gsm_nodes: int = 11,
        gsm_max_iters: int = 80,
        gsm_timeout_s: float = 1800.0,
    ) -> float:
        """Estimate a reaction barrier without DFT geometry optimisation.

        Strategy (≈20-40 min on a laptop for 33-atom systems):
          1. Optimise reactant at xTB level (or reuse cached DFT result).
          2. Refine TS guess at xTB-Sella (cheap, seconds/step) with v0
             from xTB numerical Hessian.
          3. Validate TS via xTB Hessian (require n_imag >= 1).
          4. DFT single-point on reactant & TS (one SCF each).
          5. ZPE / thermal correction from xTB Hessian.
          6. Report barrier with ``"robust_sp"`` quality tag.

        Returns barrier in kcal/mol.
        """
        import pathlib
        _ckpt = pathlib.Path(checkpoint_dir) if checkpoint_dir else None
        if _ckpt:
            _ckpt.mkdir(parents=True, exist_ok=True)

        HA_TO_KCAL = 627.509
        _cb = progress_callback or (lambda *_a, **_kw: None)

        # ── helpers ────────────────────────────────────────────────────
        def _save(name: str, content: str) -> None:
            if _ckpt:
                (_ckpt / name).write_text(content, encoding="utf-8")

        def _load(name: str) -> Optional[str]:
            if _ckpt and (_ckpt / name).exists():
                return (_ckpt / name).read_text(encoding="utf-8")
            return None

        # ── Phase R: Reactant ──────────────────────────────────────────
        # Always xTB-optimize + DFT SP so the geometry level is
        # consistent with the TS (which is also xTB geometry + DFT SP).
        # A cached DFT result is kept only for thermo corrections.
        _cb("robust_sp_reactant", {"label": "Reactant (xTB opt + DFT SP)"})
        logger.info(">>> [robust_sp] Phase R: Reactant")

        reactant_dft_e: Optional[float] = None  # Hartree
        reactant_opt_xyz: Optional[str] = None

        # Load cached DFT result for thermo only
        cached_result = None
        if _ckpt and (_ckpt / "reactant_result.json").exists():
            try:
                data = json.loads((_ckpt / "reactant_result.json").read_text())
                cached_result = DFTResult.from_dict(data)
                logger.info(
                    f"    [robust_sp] Cached DFT reactant available "
                    f"(for thermo): {cached_result.energy_hartree:.6f} Ha"
                )
            except Exception as exc:
                logger.warning(f"    [robust_sp] Cannot load cached reactant: {exc}")

        # Always xTB opt + DFT SP for energy consistency with TS
        cached_sp = _load("robust_reactant_dft_sp.txt")
        if cached_sp is not None:
            reactant_dft_e = float(cached_sp.strip())
            reactant_opt_xyz = _load("robust_reactant_xtb.xyz") or reactant_xyz
            logger.info(f"    [robust_sp] Reusing cached xTB-geom DFT SP: {reactant_dft_e:.6f} Ha")
        else:
            logger.info("    [robust_sp] Optimising reactant at xTB level...")
            reactant_opt_xyz = _xtb.run_xtb_optimization(
                reactant_xyz, charge, spin
            )
            _save("robust_reactant_xtb.xyz", reactant_opt_xyz)
            logger.info("    [robust_sp] DFT single-point on xTB-optimised reactant...")
            reactant_dft_e = self.single_point(
                reactant_opt_xyz,
                xc_method=self.opt_method,
                basis=self.opt_basis,
                charge=charge,
                spin=spin,
            )
            _save("robust_reactant_dft_sp.txt", f"{reactant_dft_e:.10f}")
            logger.info(f"    [robust_sp] Reactant DFT SP energy: {reactant_dft_e:.6f} Ha")

        # ── Phase TS: Transition State ─────────────────────────────────
        _cb("robust_sp_ts", {"label": "TS (xTB-Sella refinement + DFT SP)"})
        logger.info(">>> [robust_sp] Phase TS: xTB-Sella refinement")

        # Extract v0 from xTB Hessian on TS guess
        ts_v0 = None
        v0_path = _ckpt / "ts_v0_xtb.npy" if _ckpt else None
        if v0_path and v0_path.exists():
            try:
                ts_v0 = np.load(str(v0_path))
                logger.info(f"    [robust_sp] Loaded cached v0 ({ts_v0.shape})")
            except Exception:
                ts_v0 = None
        if ts_v0 is None:
            ts_v0 = _xtb.compute_xtb_ts_mode(ts_guess_xyz, charge=charge, spin=spin)
            if ts_v0 is not None and _ckpt:
                np.save(str(_ckpt / "ts_v0_xtb.npy"), ts_v0)

        # Refine TS at xTB-Sella with v0
        sella_result = _xtb.refine_ts_sella_xtb(
            ts_guess_xyz,
            charge=charge,
            spin=spin,
            fmax=0.05,
            max_steps=300,
            timeout_seconds=600.0,
        )
        if sella_result and sella_result.get("xyz"):
            ts_refined_xyz = sella_result["xyz"]
            logger.info(
                f"    [robust_sp] xTB-Sella converged: "
                f"fmax={sella_result['fmax_ev_ang']:.4f} eV/Å, "
                f"E={sella_result['energy_ev']:.4f} eV"
            )
        else:
            logger.warning("    [robust_sp] xTB-Sella failed, using raw TS guess.")
            ts_refined_xyz = ts_guess_xyz
        _save("robust_ts_xtb_sella.xyz", ts_refined_xyz)

        # ── Phase TS-val: xTB Hessian validation ──────────────────────
        _cb("robust_sp_ts_validation", {"label": "TS validation (xTB Hessian)"})
        logger.info(">>> [robust_sp] Phase TS-val: xTB Hessian check")

        ts_probe = _xtb.probe_ts_guess_xtb(ts_refined_xyz, charge=charge, spin=spin)
        ts_n_imag = ts_probe.get("n_imag")
        ts_lowest_freq = ts_probe.get("lowest_freq_cm")
        ts_xtb_energy_eh = ts_probe.get("energy_eh")
        logger.info(
            f"    [robust_sp] TS probe: n_imag={ts_n_imag}, "
            f"lowest_freq={ts_lowest_freq} cm⁻¹, E_xTB={ts_xtb_energy_eh} Ha"
        )

        quality_notes: List[str] = []
        if ts_n_imag is None or ts_n_imag < 1:
            quality_notes.append(
                f"TS has n_imag={ts_n_imag} at xTB level (expected 1). "
                "Barrier is an upper/lower bound estimate."
            )

        # ── GSM fallback: recover real TS when Sella fails ────────────
        if (
            enable_gsm_fallback
            and product_xyz
            and (ts_n_imag is None or ts_n_imag < 1)
        ):
            _cb("robust_sp_gsm_fallback", {"label": "GSM fallback (xTB DE-GSM)"})
            logger.info(
                ">>> [robust_sp] GSM fallback: TS has n_imag=%s, launching DE-GSM",
                ts_n_imag,
            )
            try:
                from .gsm_backend import GSMRunner

                gsm_dir = (_ckpt / "gsm_robust_sp") if _ckpt else pathlib.Path("gsm_robust_sp_scratch")
                runner = GSMRunner(
                    work_dir=gsm_dir,
                    charge=charge,
                    spin=spin,
                    n_nodes=gsm_nodes,
                    max_iters=gsm_max_iters,
                    timeout_s=gsm_timeout_s,
                )
                gsm_reactant = reactant_opt_xyz or reactant_xyz
                # Optimize product at xTB level so both endpoints are
                # at the same level of theory — raw product geometries
                # cause distorted interpolations that break xTB SCF.
                gsm_product = product_xyz
                cached_product = _load("robust_product_xtb.xyz")
                if cached_product:
                    gsm_product = cached_product
                    logger.info("    [robust_sp] Reusing cached xTB-optimised product")
                else:
                    try:
                        gsm_product = _xtb.run_xtb_optimization(
                            product_xyz, charge, spin
                        )
                        _save("robust_product_xtb.xyz", gsm_product)
                        logger.info("    [robust_sp] Product xTB-optimised for GSM")
                    except Exception as opt_exc:
                        logger.warning(
                            "    [robust_sp] Product xTB opt failed (%s), using raw",
                            opt_exc,
                        )
                gsm_result = runner.run_de_gsm(gsm_reactant, gsm_product)

                if gsm_result.converged and gsm_result.ts_xyz:
                    logger.info("    [robust_sp] GSM converged, refining GSM TS with Sella...")
                    # Re-run Sella on GSM TS
                    gsm_sella = _xtb.refine_ts_sella_xtb(
                        gsm_result.ts_xyz,
                        charge=charge,
                        spin=spin,
                        fmax=0.05,
                        max_steps=300,
                        timeout_seconds=600.0,
                    )
                    if gsm_sella and gsm_sella.get("xyz"):
                        gsm_ts_refined = gsm_sella["xyz"]
                        logger.info(
                            "    [robust_sp] GSM+Sella converged: fmax=%.4f eV/Å",
                            gsm_sella["fmax_ev_ang"],
                        )
                    else:
                        gsm_ts_refined = gsm_result.ts_xyz
                        logger.warning("    [robust_sp] GSM+Sella failed, using raw GSM TS node")

                    # Re-validate
                    gsm_probe = _xtb.probe_ts_guess_xtb(gsm_ts_refined, charge=charge, spin=spin)
                    gsm_n_imag = gsm_probe.get("n_imag")
                    logger.info(
                        "    [robust_sp] GSM TS probe: n_imag=%s, lowest_freq=%s cm⁻¹",
                        gsm_n_imag, gsm_probe.get("lowest_freq_cm"),
                    )

                    # Always adopt the GSM TS geometry: it sits at the
                    # energy maximum along the reaction coordinate, which
                    # is more physical than the NEB interpolated node even
                    # when the xTB Hessian shows n_imag=0.
                    ts_refined_xyz = gsm_ts_refined
                    ts_n_imag = gsm_n_imag if gsm_n_imag is not None else 0
                    ts_lowest_freq = gsm_probe.get("lowest_freq_cm")
                    quality_notes = [
                        n for n in quality_notes if "n_imag=" not in n
                    ]
                    if gsm_n_imag is not None and gsm_n_imag >= 1:
                        quality_notes.append(
                            f"GSM fallback recovered TS with n_imag={gsm_n_imag}"
                        )
                        logger.info(
                            "    [robust_sp] GSM recovered true TS (n_imag=%d)",
                            gsm_n_imag,
                        )
                    else:
                        quality_notes.append(
                            f"GSM path-maximum TS adopted (n_imag={gsm_n_imag}, "
                            "upper-bound estimate)"
                        )
                        logger.warning(
                            "    [robust_sp] GSM TS adopted with n_imag=%s "
                            "(path-maximum, upper-bound estimate)",
                            gsm_n_imag,
                        )
                    _save("robust_ts_xtb_sella.xyz", ts_refined_xyz)
                    # Invalidate cached DFT SP — geometry changed
                    if _ckpt and (_ckpt / "robust_ts_dft_sp.txt").exists():
                        (_ckpt / "robust_ts_dft_sp.txt").unlink()
                else:
                    reason = getattr(gsm_result, "reason", "unknown")
                    quality_notes.append(f"GSM fallback did not converge: {reason}")
                    logger.warning("    [robust_sp] GSM did not converge: %s", reason)

                # Persist GSM audit
                if _ckpt:
                    try:
                        gsm_audit = (
                            gsm_result.to_audit_dict()
                            if hasattr(gsm_result, "to_audit_dict")
                            else {"converged": gsm_result.converged}
                        )
                        _save("gsm_robust_sp_audit.json", json.dumps(gsm_audit, indent=2))
                    except Exception:
                        pass
            except Exception as gsm_exc:
                logger.warning("    [robust_sp] GSM fallback failed: %s", gsm_exc)
                quality_notes.append(
                    f"GSM fallback error: {type(gsm_exc).__name__}: {gsm_exc}"
                )

        # ── Phase SP: DFT single-points ───────────────────────────────
        _cb("robust_sp_dft_singlepoints", {"label": "DFT single-points on R & TS"})
        logger.info(">>> [robust_sp] Phase SP: DFT single-point on TS")

        cached_ts_sp = _load("robust_ts_dft_sp.txt")
        if cached_ts_sp is not None:
            ts_dft_e = float(cached_ts_sp.strip())
            logger.info(f"    [robust_sp] Reusing cached TS DFT SP: {ts_dft_e:.6f} Ha")
        else:
            ts_dft_e = self.single_point(
                ts_refined_xyz,
                xc_method=self.opt_method,
                basis=self.opt_basis,
                charge=charge,
                spin=spin,
            )
            logger.info(f"    [robust_sp] TS DFT SP energy: {ts_dft_e:.6f} Ha")
            _save("robust_ts_dft_sp.txt", f"{ts_dft_e:.10f}")

        # ── Phase thermo: ZPE from xTB ────────────────────────────────
        _cb("robust_sp_thermo", {"label": "ZPE/thermal corrections (xTB)"})
        logger.info(">>> [robust_sp] Phase thermo: xTB ZPE corrections")

        # Reactant thermo
        r_zpe_corr = 0.0
        if cached_result and cached_result.quasi_harmonic_gibbs_hartree is not None:
            r_zpe_corr = cached_result.quasi_harmonic_gibbs_hartree - cached_result.energy_hartree
            logger.info(f"    [robust_sp] Reactant thermal from cached DFT: {r_zpe_corr * HA_TO_KCAL:.2f} kcal/mol")
        else:
            r_hess = _xtb.run_xtb_hessian_thermo(
                reactant_opt_xyz or reactant_xyz, charge, spin,
                electronic_energy_h=reactant_dft_e, temp_k=self.temp_k,
            )
            if r_hess:
                r_zpe_corr = r_hess[1] - reactant_dft_e  # qh_gibbs - electronic
                logger.info(f"    [robust_sp] Reactant xTB thermal: {r_zpe_corr * HA_TO_KCAL:.2f} kcal/mol")

        # TS thermo
        ts_hess = _xtb.run_xtb_hessian_thermo(
            ts_refined_xyz, charge, spin,
            electronic_energy_h=ts_dft_e, temp_k=self.temp_k,
        )
        ts_zpe_corr = 0.0
        if ts_hess:
            ts_zpe_corr = ts_hess[1] - ts_dft_e
            logger.info(f"    [robust_sp] TS xTB thermal: {ts_zpe_corr * HA_TO_KCAL:.2f} kcal/mol")

        # ── Barrier ────────────────────────────────────────────────────
        de_electronic = (ts_dft_e - reactant_dft_e) * HA_TO_KCAL
        de_corrected = ((ts_dft_e + ts_zpe_corr) - (reactant_dft_e + r_zpe_corr)) * HA_TO_KCAL

        logger.info(
            f">>> [robust_sp] BARRIER (electronic): {de_electronic:.2f} kcal/mol"
        )
        logger.info(
            f">>> [robust_sp] BARRIER (with xTB thermal): {de_corrected:.2f} kcal/mol"
        )
        if quality_notes:
            for note in quality_notes:
                logger.warning(f"    [robust_sp] QUALITY: {note}")

        # ── Save audit JSON ────────────────────────────────────────────
        audit = {
            "strategy": "robust_sp",
            "reactant_dft_energy_h": reactant_dft_e,
            "ts_dft_energy_h": ts_dft_e,
            "reactant_thermal_correction_h": r_zpe_corr,
            "ts_thermal_correction_h": ts_zpe_corr,
            "barrier_electronic_kcal": round(de_electronic, 2),
            "barrier_corrected_kcal": round(de_corrected, 2),
            "ts_n_imag_xtb": ts_n_imag,
            "ts_lowest_freq_cm": ts_lowest_freq,
            "ts_xtb_sella_fmax": sella_result.get("fmax_ev_ang") if sella_result else None,
            "dft_method": f"{self.opt_method}/{self.opt_basis}",
            "temp_k": self.temp_k,
            "quality_notes": quality_notes,
        }
        _save("robust_sp_audit.json", json.dumps(audit, indent=2))
        _cb("robust_sp_complete", {"label": "Complete", "audit": audit})
        logger.info(f">>> [robust_sp] Done. Returning corrected barrier: {de_corrected:.2f} kcal/mol")

        return de_corrected

    def calculate_barrier(self, reactant_xyz: str, ts_xyz: str, charge: int = 0, spin: int = 0,
                          run_irc: bool = False, reaction_meta: Optional[Dict] = None,
                          checkpoint_dir: Optional[str] = None,
                          progress_callback: Optional[Callable[[str, Dict[str, Any]], None]] = None,
                          product_xyz: Optional[str] = None,
                          ts_refine_policy: str = "auto",
                          ts_refine_fmax_threshold: float = 1.0,
                          path_frames_dir: Optional[str] = None,
                          enable_gsm_fallback: bool = False,
                          gsm_nodes: int = 11,
                          gsm_max_iters: int = 80,
                          gsm_timeout_s: float = 1800.0,
                          fragment_spins: Optional[List[int]] = None) -> float:
        """
        End-to-end composite calculation of a kinetic barrier in kcal/mol.
        If run_irc=True, performs Phase 3.4 IRC validation.

        Checkpoint strategy (3 tiers, checked top-down):
          1. ``{phase}_result.json`` — full DFTResult → zero recomputation
          2. ``{phase}_optimized.xyz`` — geometry only → skip opt, redo SP+thermo
          3. Nothing → full pipeline

        Args:
            reaction_meta: Optional dict with 'reactants', 'products' (SMILES lists) and 'family'.
            checkpoint_dir: If provided, intermediate results are persisted here so partial
                runs survive crashes and restarts.
        """
        import pathlib
        _ckpt = pathlib.Path(checkpoint_dir) if checkpoint_dir else None
        if _ckpt:
            _ckpt.mkdir(parents=True, exist_ok=True)

        # ── Checkpoint helpers ─────────────────────────────────────────────
        def _save_xyz(name: str, xyz: str) -> None:
            if _ckpt is None:
                return
            try:
                (_ckpt / name).write_text(xyz, encoding="utf-8")
                logger.info(f"    [CHECKPOINT] Saved geometry {_ckpt / name}")
            except Exception as exc:
                logger.warning(f"    [CHECKPOINT] WARNING: Failed to save {name}: {exc}")

        def _load_xyz(name: str) -> Optional[str]:
            if _ckpt is None:
                return None
            path = _ckpt / name
            if path.exists():
                try:
                    xyz = path.read_text(encoding="utf-8")
                    logger.info(f"    [RESUME] Found geometry checkpoint: {path}")
                    return xyz
                except Exception as exc:
                    logger.warning(f"    [RESUME] WARNING: Failed to read checkpoint {name}: {exc}")
            return None

        def _save_result(name: str, result: DFTResult) -> None:
            """Persist the *complete* DFTResult as JSON (tier-1 checkpoint)."""
            if _ckpt is None:
                return
            try:
                (_ckpt / name).write_text(
                    json.dumps(result.to_dict(), indent=2), encoding="utf-8"
                )
                logger.info(f"    [CHECKPOINT] Saved full result {_ckpt / name}")
            except Exception as exc:
                logger.warning(f"    [CHECKPOINT] WARNING: Failed to save result {name}: {exc}")

        def _load_result(name: str) -> Optional[DFTResult]:
            """Load a tier-1 result checkpoint. Returns None if absent/corrupt."""
            if _ckpt is None:
                return None
            path = _ckpt / name
            if not path.exists():
                return None
            try:
                data = json.loads(path.read_text(encoding="utf-8"))
                res = DFTResult.from_dict(data)
                logger.info(f"    [RESUME] Restored full DFTResult from {path}")
                return res
            except Exception as exc:
                logger.warning(f"    [RESUME] WARNING: Corrupt result checkpoint {name}: {exc}")
                return None

        def _load_progress(name: str) -> Optional[Dict[str, Any]]:
            if _ckpt is None:
                return None
            path = _ckpt / name
            if not path.exists():
                return None
            try:
                data = json.loads(path.read_text(encoding="utf-8"))
                logger.info(f"    [RESUME] Restored partial progress from {path}")
                return data
            except Exception as exc:
                logger.warning(f"    [RESUME] WARNING: Corrupt progress checkpoint {name}: {exc}")
                return None

        def _run_phase_optimization(
            xyz: str,
            *,
            is_ts: bool,
            max_steps: int,
            resume_state: Optional[Dict[str, Any]] = None,
            ts_refine_ctx: Optional[Dict[str, Any]] = None,
            override_charge: Optional[int] = None,
            override_spin: Optional[int] = None,
            checkpoint_prefix: Optional[str] = None,
        ) -> DFTResult:
            phase = checkpoint_prefix or ("ts" if is_ts else "reactant")
            previous_ctx = self._optimization_checkpoint_context
            self._set_optimization_checkpoint_context(
                checkpoint_dir=_ckpt,
                phase_prefix=phase,
                resume_state=resume_state,
            )
            try:
                return self.optimize_geometry(
                    xyz,
                    charge=override_charge if override_charge is not None else charge,
                    spin=override_spin if override_spin is not None else spin,
                    is_ts=is_ts,
                    max_steps=max_steps,
                    label=target_label,
                    progress_callback=progress_callback,
                    ts_refine_ctx=ts_refine_ctx,
                )
            finally:
                self._optimization_checkpoint_context = previous_ctx

        # ── 0. Atom Consistency Check ──────────────────────────────────────
        def _get_atom_symbols(xyz: str) -> List[str]:
            atoms = []
            for line in xyz.strip().split('\n')[2:]:
                parts = line.split()
                if parts:
                    atoms.append(parts[0])
            return atoms

        r_atoms = _get_atom_symbols(reactant_xyz)
        ts_atoms = _get_atom_symbols(ts_xyz)
        if r_atoms != ts_atoms:
                   raise ValueError(
                       f"Atom sequence mismatch between reactant and TS guess! \n"
                       f"Reactant: {r_atoms}\n"
                       f"TS Guess: {ts_atoms}"
                   )

        target_label = reaction_meta.get('target_id', 'unknown') if reaction_meta else 'unknown'

        # ── 1. Reactant ────────────────────────────────────────────────────
        # Fragment-aware strategy: if the reactant XYZ contains multiple
        # disconnected molecular fragments (bimolecular reaction), each
        # fragment is optimised individually with its own charge/spin.
        # The reference energy is: G(reactant) = Σᵢ G(fragmentᵢ).
        # This avoids pathological SCF convergence for weakly-interacting
        # van der Waals complexes and is thermodynamically rigorous.
        self._emit_progress(progress_callback, "reactant_start", label="Reactant optimization starting")

        reactant_fragments = split_xyz_into_fragments(reactant_xyz)
        n_fragments = len(reactant_fragments)

        if n_fragments > 1:
            logger.info(f">>> [Phase 3.3] Detected {n_fragments} disconnected fragments in reactant — using fragment-based reference energy")
            frag_charge_spin = assign_fragment_charge_spin(
                reactant_fragments, total_charge=charge, total_spin=spin,
                manifest_fragment_spins=fragment_spins,
            )
            logger.info(f">>> [Phase 3.3] Fragment charge/spin assignments: {frag_charge_spin}")

            fragment_results: List[DFTResult] = []
            for frag_idx, (frag_xyz, (frag_charge, frag_spin)) in enumerate(
                zip(reactant_fragments, frag_charge_spin)
            ):
                frag_prefix = f"reactant_frag{frag_idx}"
                frag_label = f"fragment {frag_idx} (charge={frag_charge}, spin={frag_spin})"
                self._emit_progress(progress_callback, f"reactant_frag{frag_idx}_start",
                                    label=f"Optimising {frag_label}")

                res_frag = _load_result(f"{frag_prefix}_result.json")
                if res_frag:
                    logger.info(f">>> [Phase 3.3] Skipping {frag_label} (full result checkpoint)")
                else:
                    ckpt_frag_xyz = _load_xyz(f"{frag_prefix}_optimized.xyz")
                    if ckpt_frag_xyz:
                        logger.info(f">>> [Phase 3.3] {frag_label}: XYZ checkpoint found — SP+thermo only")
                        res_frag = _run_phase_optimization(
                            ckpt_frag_xyz, is_ts=False, max_steps=0,
                            override_charge=frag_charge, override_spin=frag_spin,
                            checkpoint_prefix=frag_prefix,
                        )
                    else:
                        frag_progress = _load_progress(f"{frag_prefix}_progress.json")
                        frag_inflight_xyz = _load_xyz(f"{frag_prefix}_inflight.xyz")
                        if frag_progress and frag_inflight_xyz:
                            if frag_progress.get("geometry_ready_for_postprocessing"):
                                logger.info(f">>> [Phase 3.3] Resuming {frag_label} — geometry ready, SP+thermo only")
                                res_frag = _run_phase_optimization(
                                    frag_inflight_xyz, is_ts=False, max_steps=0,
                                    resume_state=frag_progress,
                                    override_charge=frag_charge, override_spin=frag_spin,
                                    checkpoint_prefix=frag_prefix,
                                )
                            else:
                                logger.info(f">>> [Phase 3.3] Resuming {frag_label} — continuing geometry ladder")
                                res_frag = _run_phase_optimization(
                                    frag_inflight_xyz, is_ts=False, max_steps=100,
                                    resume_state=frag_progress,
                                    override_charge=frag_charge, override_spin=frag_spin,
                                    checkpoint_prefix=frag_prefix,
                                )
                        else:
                            logger.info(f">>> [Phase 3.3] Starting {frag_label} optimisation...")
                            res_frag = _run_phase_optimization(
                                frag_xyz, is_ts=False, max_steps=100,
                                override_charge=frag_charge, override_spin=frag_spin,
                                checkpoint_prefix=frag_prefix,
                            )
                    if not res_frag.converged:
                        logger.warning(f"[WARNING] {frag_label} optimisation hit max cycles. Using best available geometry.")
                    _save_xyz(f"{frag_prefix}_optimized.xyz", res_frag.optimized_xyz)
                    _save_result(f"{frag_prefix}_result.json", res_frag)

                fragment_results.append(res_frag)

            # Synthesise a composite res_r whose Gibbs energy is the sum of fragments
            total_energy = sum(fr.energy_hartree for fr in fragment_results)
            total_gibbs = sum(
                fr.quasi_harmonic_gibbs_hartree for fr in fragment_results
                if fr.quasi_harmonic_gibbs_hartree is not None
            )
            all_converged = all(fr.converged for fr in fragment_results)
            all_freqs: List[float] = []
            for fr in fragment_results:
                if fr.frequencies_cm1:
                    all_freqs.extend(fr.frequencies_cm1)
            combined_xyz = "\n---\n".join(fr.optimized_xyz for fr in fragment_results)

            res_r = DFTResult(
                method=fragment_results[0].method,
                energy_hartree=total_energy,
                gibbs_free_energy_hartree=None,
                quasi_harmonic_gibbs_hartree=total_gibbs if total_gibbs != 0 else None,
                optimized_xyz=combined_xyz,
                converged=all_converged,
                frequencies_cm1=all_freqs if all_freqs else None,
            )
            _save_result("reactant_result.json", res_r)
        else:
            # Single-fragment reactant — original monolithic pathway
            res_r = _load_result("reactant_result.json")   # Tier 1
            if res_r:
                logger.info(">>> [Phase 3.3] Skipping Reactant entirely (full result checkpoint)")
            else:
                ckpt_r_xyz = _load_xyz("reactant_optimized.xyz")   # Tier 2
                if ckpt_r_xyz:
                    logger.info(">>> [Phase 3.3] Skipping Reactant geometry opt (XYZ checkpoint) — SP+thermo still needed")
                    res_r = _run_phase_optimization(ckpt_r_xyz, is_ts=False, max_steps=0)
                else:
                    reactant_progress = _load_progress("reactant_progress.json")
                    reactant_inflight_xyz = _load_xyz("reactant_inflight.xyz")
                    if reactant_progress and reactant_inflight_xyz:
                        if reactant_progress.get("geometry_ready_for_postprocessing"):
                            logger.info(">>> [Phase 3.3] Resuming Reactant from partial checkpoint — geometry ready, SP+thermo only")
                            res_r = _run_phase_optimization(
                                reactant_inflight_xyz,
                                is_ts=False,
                                max_steps=0,
                                resume_state=reactant_progress,
                            )
                        else:
                            logger.info(">>> [Phase 3.3] Resuming Reactant from partial checkpoint — continuing geometry ladder")
                            res_r = _run_phase_optimization(
                                reactant_inflight_xyz,
                                is_ts=False,
                                max_steps=100,
                                resume_state=reactant_progress,
                            )
                    else:
                        logger.info(">>> [Phase 3.3] Starting Reactant Optimization...")
                        res_r = _run_phase_optimization(reactant_xyz, is_ts=False, max_steps=100)
                    if not res_r.converged:
                        logger.warning("[WARNING] Reactant optimization hit max cycles. Using best available geometry.")
                    _save_xyz("reactant_optimized.xyz", res_r.optimized_xyz)
                # After SP+thermo complete, persist the full result
                _save_result("reactant_result.json", res_r)

        logger.info(">>> [Phase 3.3] Reactant refinement complete.")

        # ── 2. Transition State ────────────────────────────────────────────
        self._emit_progress(progress_callback, "ts_start", label="Transition-state optimization starting")

        # Build P0 TS-guess refinement context
        _ts_refine_ctx: Optional[Dict[str, Any]] = None
        if ts_refine_policy != "off":
            _ts_refine_ctx = {
                "reactant_xyz": reactant_xyz,
                "product_xyz": product_xyz,
                "path_frames_dir": path_frames_dir,
                "policy": ts_refine_policy,
                "threshold_fmax": ts_refine_fmax_threshold,
                "enable_gsm_fallback": enable_gsm_fallback,
                "gsm_nodes": gsm_nodes,
                "gsm_max_iters": gsm_max_iters,
                "gsm_timeout_s": gsm_timeout_s,
            }

        ts_seed_xyz = ts_xyz
        if reaction_meta and reaction_meta.get('forming_bond_atoms'):
            recovery = self.recover_ts_guess_from_forming_bond(
                reactant_xyz,
                ts_seed_xyz,
                charge=charge,
                spin=spin,
                reaction_meta=reaction_meta,
                force_scan=False,
                progress_callback=progress_callback,
            )
            ts_seed_xyz = str(recovery.get("ts_xyz", ts_seed_xyz))

        res_ts = _load_result("ts_result.json")   # Tier 1
        if res_ts:
            logger.info(">>> [Phase 3.3] Skipping TS entirely (full result checkpoint)")
        else:
            ckpt_ts_xyz = _load_xyz("ts_optimized.xyz")   # Tier 2
            if ckpt_ts_xyz:
                logger.info(">>> [Phase 3.3] Skipping TS geometry opt (XYZ checkpoint) — SP+thermo still needed")
                res_ts = _run_phase_optimization(ckpt_ts_xyz, is_ts=True, max_steps=0)
            else:
                ts_progress = _load_progress("ts_progress.json")
                ts_inflight_xyz = _load_xyz("ts_inflight.xyz")
                if ts_progress and ts_inflight_xyz:
                    if ts_progress.get("geometry_ready_for_postprocessing"):
                        logger.info(">>> [Phase 3.3] Resuming TS from partial checkpoint — geometry ready, SP+thermo only")
                        res_ts = _run_phase_optimization(
                            ts_inflight_xyz,
                            is_ts=True,
                            max_steps=0,
                            resume_state=ts_progress,
                        )
                    else:
                        logger.info(">>> [Phase 3.3] Resuming TS from partial checkpoint — continuing geometry ladder")
                        res_ts = _run_phase_optimization(
                            ts_inflight_xyz,
                            is_ts=True,
                            max_steps=100,
                            resume_state=ts_progress,
                            ts_refine_ctx=_ts_refine_ctx,
                        )
                else:
                    logger.info(">>> [Phase 3.3] Starting Transition State (TS) Optimization...")
                    res_ts = _run_phase_optimization(ts_seed_xyz, is_ts=True, max_steps=100, ts_refine_ctx=_ts_refine_ctx)
                if not res_ts.converged:
                    logger.warning("[WARNING] TS optimization hit max cycles. Using best available geometry.")
                _save_xyz("ts_optimized.xyz", res_ts.optimized_xyz)
            _save_result("ts_result.json", res_ts)
            
        logger.info(">>> [Phase 3.3] Transition State refinement complete.")

        # 2.5 Validation of TS nature (Phase-2 hardening).
        # Classify the spectrum and persist the verdict so downstream tooling
        # (run_computational_gap_dft.py promotion gate) can refuse to promote
        # a "barrier" that came from a minimum or a higher-order saddle.
        ts_validation = _classify_ts_frequencies(res_ts.frequencies_cm1)
        verdict = ts_validation["verdict"]
        n_sig = ts_validation["n_significant_imaginary"]
        n_imag = ts_validation["n_imaginary"]
        floor = ts_validation["floor_cm1"]
        if verdict == "true_ts":
            logger.info(
                f"    [VALIDATION] True TS verified: 1 significant imaginary mode "
                f"({ts_validation['most_negative_cm1']:.1f} cm^-1, |f| > {floor:.0f})."
            )
        elif verdict == "minimum":
            logger.error(
                "[SOTA WARNING] TS optimization converged to a minimum "
                "(0 imaginary frequencies). Reported barrier will be flagged "
                "as invalid_saddle and will NOT be promoted."
            )
        elif verdict == "weak_imaginary":
            logger.error(
                f"[SOTA WARNING] TS has {n_imag} imaginary mode(s) but none "
                f"exceed the |f| > {floor:.0f} cm^-1 floor "
                "(numerical noise, not a real reaction coordinate). "
                "Will be flagged as invalid_saddle."
            )
        elif verdict == "high_order_saddle":
            logger.warning(
                f"[SOTA WARNING] TS has {n_sig} significant imaginary modes — "
                "higher-order saddle, not a true TS. Will be flagged as invalid_saddle."
            )
        else:  # no_frequencies
            logger.warning(
                "[SOTA WARNING] No frequency analysis available for TS — "
                "cannot validate saddle nature. Will be flagged as invalid_saddle."
            )
        # Attach to result for in-process callers and persist to disk for the
        # promotion-gate script.
        try:
            setattr(res_ts, "ts_validation", ts_validation)
        except Exception:
            pass
        if _ckpt is not None:
            try:
                _ckpt.mkdir(parents=True, exist_ok=True)
                (_ckpt / "ts_validation.json").write_text(
                    json.dumps(ts_validation, indent=2), encoding="utf-8"
                )
            except Exception as exc:
                logger.warning(f"Failed to persist ts_validation.json: {exc}")

        # 2.6 pyGSM hard fallback — opt-in retry when the validated TS is
        # not a true saddle.  We invalidate the TS checkpoint, replace the
        # seed with the GSM-derived geometry, and rerun the TS phase
        # exactly once.  The cascade inside the rerun has its own GSM
        # toggle disabled to guarantee a single attempt.
        gsm_retry_attempted = False
        if (
            enable_gsm_fallback
            and product_xyz
            and not ts_validation.get("is_true_saddle", False)
            and _ts_refine_ctx is not None
        ):
            gsm_retry_attempted = True
            logger.info(
                ">>> [Phase 3.3.5] TS failed saddle gate (verdict=%s); "
                "launching pyGSM hard fallback.", verdict,
            )
            gsm_outcome: Dict[str, Any] = {
                "triggered": True,
                "verdict_pre_gsm": verdict,
            }
            try:
                from .gsm_backend import GSMRunner

                gsm_dir = (_ckpt / "gsm_retry") if _ckpt else Path("gsm_retry")
                runner = GSMRunner(
                    work_dir=gsm_dir,
                    charge=charge,
                    spin=spin,
                    n_nodes=int(_ts_refine_ctx.get("gsm_nodes", gsm_nodes)),
                    max_iters=int(_ts_refine_ctx.get("gsm_max_iters", gsm_max_iters)),
                    timeout_s=float(_ts_refine_ctx.get("gsm_timeout_s", gsm_timeout_s)),
                )
                gsm_result = runner.run_de_gsm(reactant_xyz, product_xyz)
                gsm_outcome["pygsm"] = gsm_result.to_audit_dict()

                if gsm_result.converged and gsm_result.ts_xyz:
                    # Compute xTB Hessian at the GSM TS node to extract the
                    # imaginary mode.  Passing this as v0 to Sella prevents
                    # the optimizer from collapsing to a minimum when the DFT
                    # and xTB PES topologies diverge.
                    ts_v0 = None
                    try:
                        from .xtb_backend import compute_xtb_ts_mode
                        ts_v0 = compute_xtb_ts_mode(
                            gsm_result.ts_xyz, charge=charge, spin=spin,
                        )
                        if ts_v0 is not None:
                            gsm_outcome["ts_v0_extracted"] = True
                            logger.info("    [Phase 3.3.5] xTB TS mode extracted for Sella v0.")
                            # Persist v0 so it survives checkpoint resume.
                            if _ckpt:
                                try:
                                    np.save(str(_ckpt / "ts_v0_xtb.npy"), ts_v0)
                                except Exception:
                                    pass
                        else:
                            gsm_outcome["ts_v0_extracted"] = False
                            logger.warning(
                                "    [Phase 3.3.5] No negative eigenvalue at xTB level; "
                                "Sella will run without v0."
                            )
                    except Exception as v0_exc:
                        logger.warning(f"    [Phase 3.3.5] v0 extraction failed: {v0_exc}")
                        gsm_outcome["ts_v0_error"] = str(v0_exc)

                    # Invalidate stale TS checkpoints so the retry actually
                    # recomputes, but keep the originals for forensics.
                    if _ckpt:
                        for stale in ("ts_result.json", "ts_optimized.xyz",
                                      "ts_progress.json", "ts_inflight.xyz"):
                            stale_path = _ckpt / stale
                            if stale_path.exists():
                                try:
                                    stale_path.rename(stale_path.with_suffix(stale_path.suffix + ".pre_gsm"))
                                except Exception:
                                    pass
                    # Disable GSM inside the retry cascade — single attempt only.
                    retry_ctx = dict(_ts_refine_ctx)
                    retry_ctx["enable_gsm_fallback"] = False
                    if ts_v0 is not None:
                        retry_ctx["ts_v0"] = ts_v0
                    res_ts = _run_phase_optimization(
                        gsm_result.ts_xyz, is_ts=True, max_steps=100,
                        ts_refine_ctx=retry_ctx,
                    )
                    if _ckpt:
                        _save_xyz("ts_optimized.xyz", res_ts.optimized_xyz)
                        _save_result("ts_result.json", res_ts)
                    # Reclassify with the post-GSM frequencies.
                    ts_validation = _classify_ts_frequencies(res_ts.frequencies_cm1)
                    verdict = ts_validation["verdict"]
                    gsm_outcome["verdict_post_gsm"] = verdict
                    gsm_outcome["recovered_true_saddle"] = bool(
                        ts_validation.get("is_true_saddle", False)
                    )
                    if _ckpt:
                        try:
                            (_ckpt / "ts_validation.json").write_text(
                                json.dumps(ts_validation, indent=2), encoding="utf-8"
                            )
                        except Exception as exc:
                            logger.warning(f"Failed to persist post-GSM ts_validation.json: {exc}")
                    if gsm_outcome["recovered_true_saddle"]:
                        logger.info(
                            "    [Phase 3.3.5] pyGSM recovered a true TS "
                            "(verdict=%s).", verdict,
                        )
                    else:
                        logger.warning(
                            "    [Phase 3.3.5] pyGSM TS rerun still failed "
                            "the saddle gate (verdict=%s).", verdict,
                        )
                    try:
                        setattr(res_ts, "ts_validation", ts_validation)
                    except Exception:
                        pass
                else:
                    logger.warning(
                        "    [Phase 3.3.5] pyGSM did not converge: %s",
                        gsm_result.reason,
                    )
                    gsm_outcome["recovered_true_saddle"] = False
            except Exception as gsm_exc:
                logger.warning(f"    [Phase 3.3.5] pyGSM fallback failed: {gsm_exc}")
                gsm_outcome["error"] = f"{type(gsm_exc).__name__}: {gsm_exc}"
                gsm_outcome["recovered_true_saddle"] = False

            if _ckpt:
                try:
                    (_ckpt / "gsm_attempt.json").write_text(
                        json.dumps(gsm_outcome, indent=2), encoding="utf-8"
                    )
                except Exception as exc:
                    logger.warning(f"Failed to persist gsm_attempt.json: {exc}")

        # 3. IRC Validation (Phase 3.4)
        if run_irc:
            self._emit_progress(progress_callback, "irc_validation", label="IRC validation")
            try:
                self.generate_irc(res_ts.optimized_xyz, charge=charge, spin=spin)
            except Exception as irc_exc:
                logger.warning(f"[WARNING] IRC validation failed: {irc_exc}. Proceeding with barrier evaluation.")
            
        # 4. Delta G‡
        self._emit_progress(progress_callback, "barrier_evaluation", label="Barrier evaluation")
        if res_ts.quasi_harmonic_gibbs_hartree is None:
            raise RuntimeError(
                "TS quasi-harmonic Gibbs energy is None — Hessian/thermo step likely failed."
            )
        if res_r.quasi_harmonic_gibbs_hartree is None:
            raise RuntimeError(
                "Reactant quasi-harmonic Gibbs energy is None — Hessian/thermo step likely failed."
            )
        delta_g_h = res_ts.quasi_harmonic_gibbs_hartree - res_r.quasi_harmonic_gibbs_hartree
        barrier_kcal = delta_g_h * 627.509
        logger.info(f">>> [Phase 3.3] Barrier evaluation complete: {barrier_kcal:.2f} kcal/mol")
        self._emit_progress(progress_callback, "barrier_complete", label="Barrier evaluation complete", barrier_kcal=barrier_kcal)
        
        # 5. Log to DB if available
        if self.db and reaction_meta:
            # Accurate data provenance: ensure the database reflects if this was
            # a fully converged optimization or a watchdog-saved "best guess"
            is_fully_converged = res_r.converged and res_ts.converged
            db_reactants = _resolve_reaction_side(reaction_meta, "reactants")
            db_products = _resolve_reaction_side(reaction_meta, "products")
            if not db_reactants and not db_products:
                logger.warning(
                    "[DB] Not persisting barrier for %s: reaction_meta carries no reactant/product "
                    "SMILES, so the row would have no addressable reaction identity.",
                    reaction_meta.get('target_id', 'unknown'),
                )
            else:
                self.db.add_barrier(
                    reactants=db_reactants,
                    products=db_products,
                    family=_normalize_reaction_family(reaction_meta.get('family', 'unknown')),
                    delta_g_kcal=barrier_kcal,
                    method=res_ts.method,
                    basis=self.refinement_basis,
                    solvation=self.solvent_name,
                    converged=is_fully_converged
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
            db_reactants = _resolve_reaction_side(reaction_meta, "reactants")
            db_products = _resolve_reaction_side(reaction_meta, "products")
            if not db_reactants and not db_products:
                logger.warning(
                    "[DB] Not persisting MLP barrier for %s: reaction_meta carries no "
                    "reactant/product SMILES.",
                    reaction_meta.get('target_id', 'unknown'),
                )
            else:
                self.db.add_barrier(
                    reactants=db_reactants,
                    products=db_products,
                    family=_normalize_reaction_family(reaction_meta.get('family', 'unknown')),
                    delta_g_kcal=barrier_kcal,
                    method="mace-off",
                    basis="OFF23_medium", # Default model name
                    solvation=self.solvent_name,
                    converged=True
                )
            
        return barrier_kcal


def compute_irc(
    ts_xyz: str,
    *,
    charge: int = 0,
    spin: int = 0,
    step_size: float = 0.05,
    forward: bool = True,
    backward: bool = True,
    refiner: Optional[DFTRefiner] = None,
) -> Dict[str, Any]:
    if not forward and not backward:
        raise ValueError("At least one IRC direction must be enabled")

    engine = refiner or DFTRefiner()
    backward_xyz, forward_xyz = engine.generate_irc(ts_xyz, charge=charge, spin=spin, step_size=step_size)

    points: List[Tuple[str, float]] = []
    if backward:
        points.append((backward_xyz, float(engine.single_point(backward_xyz, charge=charge, spin=spin))))
    points.append((ts_xyz, float(engine.single_point(ts_xyz, charge=charge, spin=spin))))
    if forward:
        points.append((forward_xyz, float(engine.single_point(forward_xyz, charge=charge, spin=spin))))

    max_point, max_energy = max(points, key=lambda item: item[1])
    return {
        "path_type": "double_optimization_proxy",
        "backward_endpoint": backward_xyz if backward else None,
        "forward_endpoint": forward_xyz if forward else None,
        "energies": [energy for _xyz, energy in points],
        "max_energy": max_energy,
        "max_point": max_point,
    }
