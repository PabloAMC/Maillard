"""xTB (GFN2) subprocess backend — stateless helpers for geometry
optimization, TS pre-relaxation, Hessian/thermo and TS-guess refinement."""

from __future__ import annotations

import io
import os
import subprocess
import tempfile
from functools import lru_cache
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .logger import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# ASE Calculator factory (cached per (charge, spin, solvent))
# ---------------------------------------------------------------------------

@lru_cache(maxsize=8)
def _tblite_available() -> bool:
    try:
        from tblite.ase import TBLite  # noqa: F401
        return True
    except ImportError:
        return False


def get_xtb_ase_calculator(
    charge: int = 0,
    spin: int = 0,
    solvent: Optional[str] = "water",
):
    """Return an ASE Calculator wrapping GFN2-xTB.

    Prefers *tblite-python* (native, fast).  Falls back to a thin
    subprocess wrapper that shells out to the ``xtb`` binary.
    """
    uhf = abs(spin)
    if _tblite_available():
        from tblite.ase import TBLite
        kwargs: Dict[str, Any] = {
            "method": "GFN2-xTB",
            "charge": charge,
            "uhf": uhf,
            "verbosity": 0,
        }
        # tblite solvent support is version-dependent; try ALPB
        if solvent:
            try:
                return TBLite(alpb_solvation=solvent, **kwargs)
            except TypeError:
                pass
        return TBLite(**kwargs)
    # Fallback: subprocess-based calculator (slower but always works)
    return _SubprocessXTBCalc(charge=charge, uhf=uhf, solvent=solvent)


class _SubprocessXTBCalc:
    """Minimal ASE-compatible calculator that delegates to ``xtb`` CLI."""

    implemented_properties = ["energy", "forces"]
    name = "xtb-subprocess"

    def __init__(self, charge: int = 0, uhf: int = 0, solvent: Optional[str] = "water"):
        self.charge = charge
        self.uhf = uhf
        self.solvent = solvent
        self.results: Dict[str, Any] = {}

    # ASE Calculator protocol -------------------------------------------------
    def calculate(self, atoms=None, properties=None, system_changes=None):  # noqa: D401
        if atoms is None:
            raise RuntimeError("_SubprocessXTBCalc requires atoms")
        from ase.io import read as ase_read, write as ase_write
        with tempfile.TemporaryDirectory() as td:
            inp = os.path.join(td, "input.xyz")
            ase_write(inp, atoms, format="xyz")
            cmd = [
                "xtb", "input.xyz", "--grad",
                "--chrg", str(self.charge), "--uhf", str(self.uhf),
            ]
            if self.solvent:
                cmd += ["--alpb", self.solvent]
            res = subprocess.run(cmd, cwd=td, capture_output=True, text=True, timeout=300)
            if res.returncode != 0:
                raise RuntimeError(f"xTB subprocess failed: {res.stderr[:300]}")
            # Parse energy from stdout
            energy_ha = None
            for line in res.stdout.splitlines():
                if "TOTAL ENERGY" in line:
                    parts = line.split()
                    for i, p in enumerate(parts):
                        if p == "TOTAL" and i + 2 < len(parts) and parts[i + 1] == "ENERGY":
                            energy_ha = float(parts[i + 2])
                            break
            if energy_ha is None:
                raise RuntimeError("Could not parse xTB energy")
            # Parse gradient file
            grad_path = os.path.join(td, "gradient")
            if not os.path.exists(grad_path):
                raise RuntimeError("xTB did not produce gradient file")
            grad_lines = open(grad_path).readlines()
            n = len(atoms)
            # xTB gradient file: header, coords, then gradients in Hartree/Bohr
            grad_start = None
            for i, line in enumerate(grad_lines):
                if "cycle" in line.lower():
                    grad_start = i + 1
            if grad_start is None:
                grad_start = 1
            # Skip n coordinate lines, then read n gradient lines
            grad_data = []
            for line in grad_lines[grad_start + n: grad_start + 2 * n]:
                parts = line.split()
                grad_data.append([float(x) for x in parts[:3]])
            if len(grad_data) != n:
                raise RuntimeError(f"Expected {n} gradient lines, got {len(grad_data)}")
            grad_arr = np.array(grad_data)
            # Convert: gradient (Ha/Bohr) → forces (eV/Å) = -grad * 27.2114 / 0.529177
            ha_to_ev = 27.211386245988
            bohr_to_ang = 0.529177210903
            forces_ev_ang = -grad_arr * ha_to_ev / bohr_to_ang
            self.results = {
                "energy": energy_ha * ha_to_ev,
                "forces": forces_ev_ang,
            }

    def get_potential_energy(self, atoms=None, force_consistent=False):
        if "energy" not in self.results:
            self.calculate(atoms)
        return self.results["energy"]

    def get_forces(self, atoms=None):
        if "forces" not in self.results:
            self.calculate(atoms)
        return self.results["forces"]


def _write_xyz(td: str, xyz_content: str) -> str:
    """Write *xyz_content* to ``input.xyz`` inside *td* and return the path."""
    p = os.path.join(td, "input.xyz")
    with open(p, "w") as f:
        f.write(xyz_content)
    return p


def _read_opt_xyz(td: str, tag: str = "xTB") -> str:
    """Read and validate ``xtbopt.xyz`` from *td*."""
    out = os.path.join(td, "xtbopt.xyz")
    if not os.path.exists(out):
        raise RuntimeError(f"{tag} produced no xtbopt.xyz")
    with open(out) as f:
        data = f.read()
    if not data.strip():
        raise RuntimeError(f"{tag} produced an empty xtbopt.xyz")
    return data


def run_xtb_optimization(xyz_content: str, charge: int, spin: int) -> str:
    """GFN2-xTB tight optimization with ALPB Water (reactant ladder Stage 2)."""
    logger.info("    [LADDER] Stage 2: Running xTB (GFN2) ALPB-Water optimization...")
    with tempfile.TemporaryDirectory() as td:
        _write_xyz(td, xyz_content)
        cmd = ["xtb", "input.xyz", "--opt", "tight",
               "--chrg", str(charge), "--uhf", str(spin), "--alpb", "water"]
        try:
            subprocess.run(cmd, cwd=td, capture_output=True, text=True,
                           check=True, errors="replace")
            return _read_opt_xyz(td)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"xTB optimization failed (code {e.returncode}): {e.stderr[:300]}"
            ) from e


def run_xtb_ts_safe_prerelax(
    xyz_content: str, charge: int, spin: int, max_cycles: int = 20,
) -> str:
    """Crude xTB pre-relax safe for TS guesses (gas-phase, cycle-capped)."""
    logger.info(
        f"    [LADDER] Stage 2 (TS-safe): Running xTB (GFN2) crude pre-relax (max {max_cycles} cycles)..."
    )
    with tempfile.TemporaryDirectory() as td:
        _write_xyz(td, xyz_content)
        # Cycle cap via xcontrol (xTB has no --cycles CLI flag).
        xctrl = os.path.join(td, "xcontrol")
        with open(xctrl, "w") as f:
            f.write(f"$opt\n   maxcycle={int(max_cycles)}\n$end\n")
        cmd = ["xtb", "input.xyz", "--opt", "crude", "--input", "xcontrol",
               "--chrg", str(charge), "--uhf", str(spin)]
        try:
            subprocess.run(cmd, cwd=td, capture_output=True, text=True,
                           check=True, errors="replace")
        except subprocess.CalledProcessError as e:
            # Non-zero exit when cycle cap is hit before convergence — expected.
            if not os.path.exists(os.path.join(td, "xtbopt.xyz")):
                raise RuntimeError(
                    f"xTB TS-safe pre-relax failed (code {e.returncode}): {e.stderr[:300]}"
                ) from e
            logger.info("    [LADDER] xTB TS-safe pre-relax hit cycle cap (expected). Accepting last frame.")
        return _read_opt_xyz(td, tag="xTB TS-safe pre-relax")


def run_xtb_hessian_thermo(
    xyz_text: str, charge: int, spin: int,
    electronic_energy_h: float, temp_k: float,
) -> Optional[Tuple[List[float], float, float]]:
    """GFN2-xTB Hessian → quasi-harmonic Gibbs using DFT electronic energy.

    Returns ``(frequencies, qh_gibbs_h, qh_gibbs_h)`` or ``None`` on failure.
    """
    from .thermo import QuasiHarmonicCorrector
    try:
        with tempfile.TemporaryDirectory() as td:
            inp = os.path.join(td, "mol.xyz")
            with open(inp, "w") as f:
                f.write(xyz_text)
            cmd = ["xtb", inp, "--hess", "--chrg", str(charge),
                   "--uhf", str(spin), "--alpb", "water"]
            result = subprocess.run(cmd, cwd=td, capture_output=True, text=True, timeout=600)
            if result.returncode != 0:
                logger.warning(f"    [xTB-HESS] xTB hessian failed (code {result.returncode})")
                return None
            freqs = parse_xtb_frequencies(result.stdout)
            if not freqs:
                logger.warning("    [xTB-HESS] Could not parse xTB frequencies")
                return None
            logger.info(f"    [xTB-HESS] Got {len(freqs)} frequencies (lowest: {freqs[0]:.1f} cm⁻¹)")
            corrector = QuasiHarmonicCorrector(temp_k=temp_k)
            qh = corrector.calculate_thermo(freqs, electronic_energy_h=electronic_energy_h)
            thermal_kcal = (qh.qh_gibbs_h - electronic_energy_h) * 627.509
            if abs(thermal_kcal) > 40.0:
                logger.warning(f"    [SCIENTIFIC CAUTION] Large xTB thermal correction: {thermal_kcal:.1f} kcal/mol.")
            return freqs, qh.qh_gibbs_h, qh.qh_gibbs_h
    except FileNotFoundError:
        logger.warning("    [xTB-HESS] xTB binary not found; skipping.")
        return None
    except subprocess.TimeoutExpired:
        logger.warning("    [xTB-HESS] xTB hessian timed out (>600s); skipping.")
        return None
    except Exception as exc:
        logger.warning(f"    [xTB-HESS] xTB hessian fallback failed ({type(exc).__name__}): {exc}")
        return None


def parse_xtb_frequencies(stdout: str) -> List[float]:
    """Extract vibrational frequencies from xTB stdout."""
    freqs: List[float] = []
    in_freq_block = False
    for line in stdout.splitlines():
        if "Frequency" in line and "Printout" in line:
            in_freq_block = True
            continue
        if in_freq_block:
            if not line.strip():
                break
            parts = line.split()
            if len(parts) >= 2:
                try:
                    freqs.append(float(parts[1]))
                except ValueError:
                    pass
    if not freqs:
        capture = False
        for line in stdout.splitlines():
            if "projected vibrational frequencies" in line.lower():
                capture = True
                continue
            if capture:
                stripped = line.strip()
                if not stripped:
                    if freqs:
                        break
                    continue
                for tok in stripped.split():
                    try:
                        freqs.append(float(tok))
                    except ValueError:
                        pass
    return freqs


# ---------------------------------------------------------------------------
# P0 — TS-guess quality probe & refinement engines
# ---------------------------------------------------------------------------

def _xyz_to_ase_atoms(xyz_text: str):
    """Parse an XYZ string into an ASE Atoms object."""
    from ase.io import read as ase_read
    with tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False) as f:
        f.write(xyz_text)
        f.flush()
        atoms = ase_read(f.name, format="xyz")
    os.unlink(f.name)
    return atoms


def _ase_atoms_to_xyz(atoms) -> str:
    """Serialize ASE Atoms to XYZ string."""
    from ase.io import write as ase_write
    with tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False) as f:
        name = f.name
    ase_write(name, atoms, format="xyz")
    with open(name) as f:
        xyz = f.read()
    os.unlink(name)
    return xyz


def probe_ts_guess_xtb(
    xyz_text: str,
    charge: int = 0,
    spin: int = 0,
) -> Dict[str, Any]:
    """Evaluate TS-guess quality at GFN2-xTB level.

    Returns a dict with keys:
      - ``fmax_ev_ang``: max atomic force norm (eV/Å)
      - ``energy_eh``: total energy (Hartree)
      - ``n_imag``: number of imaginary frequencies
      - ``lowest_freq_cm``: lowest vibrational frequency (cm⁻¹)
      - ``highest_imag_cm``: most negative frequency (cm⁻¹), or None
    """
    logger.info("    [TS-PROBE] Computing xTB forces + Hessian on TS guess...")
    atoms = _xyz_to_ase_atoms(xyz_text)
    calc = get_xtb_ase_calculator(charge=charge, spin=spin, solvent="water")
    atoms.calc = calc

    # Forces → fmax
    try:
        forces = atoms.get_forces()
        fmax = float(np.sqrt((forces ** 2).sum(axis=1).max()))
        energy_ev = float(atoms.get_potential_energy())
    except Exception as exc:
        logger.warning(f"    [TS-PROBE] xTB force evaluation failed: {exc}")
        return {"fmax_ev_ang": float("inf"), "energy_eh": None, "n_imag": None,
                "lowest_freq_cm": None, "highest_imag_cm": None, "error": str(exc)}

    ha_to_ev = 27.211386245988
    energy_eh = energy_ev / ha_to_ev

    # Hessian via subprocess (tblite doesn't expose analytical Hessian easily)
    hess_result = run_xtb_hessian_thermo(
        xyz_text, charge, spin,
        electronic_energy_h=energy_eh, temp_k=298.15,
    )
    n_imag = None
    lowest_freq = None
    highest_imag = None
    if hess_result is not None:
        freqs = hess_result[0]
        imag_freqs = [f for f in freqs if f < -10.0]  # threshold for noise
        n_imag = len(imag_freqs)
        lowest_freq = freqs[0] if freqs else None
        highest_imag = min(imag_freqs) if imag_freqs else None

    logger.info(
        f"    [TS-PROBE] fmax={fmax:.4f} eV/Å, E={energy_eh:.6f} Ha, "
        f"n_imag={n_imag}, lowest_freq={lowest_freq}"
    )
    return {
        "fmax_ev_ang": fmax,
        "energy_eh": energy_eh,
        "n_imag": n_imag,
        "lowest_freq_cm": lowest_freq,
        "highest_imag_cm": highest_imag,
    }


def compute_xtb_ts_mode(
    xyz_text: str,
    charge: int = 0,
    spin: int = 0,
    delta: float = 0.005,
) -> Optional[np.ndarray]:
    """Compute the xTB Hessian numerically and return the TS eigenvector.

    Returns a ``(3*N,)`` array corresponding to the eigenvector of the
    most-negative eigenvalue, or ``None`` if the Hessian has no negative
    eigenvalue (i.e. the geometry is already a minimum at xTB level).

    Parameters
    ----------
    delta : float
        Finite-difference displacement in Å.
    """
    logger.info("    [TS-MODE] Computing xTB numerical Hessian for TS mode extraction...")
    atoms = _xyz_to_ase_atoms(xyz_text)
    calc = get_xtb_ase_calculator(charge=charge, spin=spin, solvent="water")
    atoms.calc = calc

    n = len(atoms)
    ndof = 3 * n
    hessian = np.zeros((ndof, ndof))

    pos0 = atoms.positions.copy()
    try:
        # Central finite differences: H_ij = (f_i(+dj) - f_i(-dj)) / (2*delta)
        for j in range(ndof):
            atom_j, coord_j = divmod(j, 3)
            # +delta
            atoms.positions = pos0.copy()
            atoms.positions[atom_j, coord_j] += delta
            fp = atoms.get_forces().ravel()  # (3N,) in eV/Å
            # -delta
            atoms.positions = pos0.copy()
            atoms.positions[atom_j, coord_j] -= delta
            fm = atoms.get_forces().ravel()
            # Hessian column (force = -gradient, so H = -(f+ - f-)/(2*delta))
            hessian[:, j] = -(fp - fm) / (2.0 * delta)

        # Symmetrise
        hessian = 0.5 * (hessian + hessian.T)

        # Diagonalise
        eigenvalues, eigenvectors = np.linalg.eigh(hessian)
        # eigenvalues are in eV/Å²; most negative first
        n_neg = int(np.sum(eigenvalues < -1.0))  # threshold 1 eV/Å²
        logger.info(
            f"    [TS-MODE] Hessian eigenvalues: lowest={eigenvalues[0]:.2f} eV/Å², "
            f"n_negative(>1)={n_neg}"
        )

        if eigenvalues[0] < -1.0:
            v0 = eigenvectors[:, 0]  # eigenvector of most negative eigenvalue
            logger.info(f"    [TS-MODE] Extracted TS mode (eigenvalue={eigenvalues[0]:.2f} eV/Å²)")
            return v0
        else:
            logger.warning("    [TS-MODE] No significant negative eigenvalue found at xTB level.")
            return None

    except Exception as exc:
        logger.warning(f"    [TS-MODE] xTB Hessian computation failed: {exc}")
        return None
    finally:
        atoms.positions = pos0


def validate_ts_mode(
    imag_vector: np.ndarray,
    expected_atoms: List[int],
    *,
    threshold: float = 0.3,
) -> Dict[str, Any]:
    """Sanity-check that an imaginary mode is concentrated on expected atoms.

    Given the eigenvector of the negative Hessian eigenvalue (shape ``(3*N,)``,
    typically returned by :func:`compute_xtb_ts_mode`), compute the fraction of
    its squared norm carried by the atoms in ``expected_atoms``.

    A mode that truly represents the targeted reaction coordinate (e.g. the
    transferred H plus the donor and acceptor) should concentrate most of its
    motion on those atoms; a "fortuitous" saddle point — for example a methyl
    rotation that happens to be a low-curvature mode — will not.

    Returns a dict with:
      - ``concentration``: float in ``[0, 1]`` — fraction on expected atoms
      - ``pass``: bool — ``concentration >= threshold``
      - ``top_atoms``: list of (atom_idx, fraction) sorted descending,
        truncated to the top 5
    """
    v = np.asarray(imag_vector, dtype=float)
    if v.ndim != 1 or v.size % 3 != 0:
        raise ValueError(f"validate_ts_mode: expected 1-D (3*N,) vector, got shape {v.shape}")
    n_atoms = v.size // 3
    per_atom = (v.reshape(n_atoms, 3) ** 2).sum(axis=1)
    total = float(per_atom.sum())
    if total <= 0.0:
        return {"concentration": 0.0, "pass": False, "top_atoms": []}
    fractions = per_atom / total
    expected = [int(i) for i in expected_atoms if 0 <= int(i) < n_atoms]
    concentration = float(fractions[expected].sum()) if expected else 0.0
    top_idx = np.argsort(fractions)[::-1][:5]
    top_atoms = [(int(i), float(fractions[i])) for i in top_idx]
    return {
        "concentration": concentration,
        "pass": concentration >= threshold,
        "top_atoms": top_atoms,
    }


def detect_forming_bond(
    reactant_xyz: str,
    product_xyz: str,
    *,
    max_product_distance: float = 2.0,
) -> Optional[Tuple[int, int]]:
    """Identify the atom pair forming a new bond between reactant and product.

    Compares pairwise distances; returns the (i, j) with the largest
    shortening whose product distance is < *max_product_distance* Å
    (i.e. a chemical bond in the product).  Returns None if ambiguous.

    Heavy-heavy atom pairs are preferred over X-H pairs to avoid selecting
    proton transfers instead of the actual bond-forming step.
    """
    r_atoms = _xyz_to_ase_atoms(reactant_xyz)
    p_atoms = _xyz_to_ase_atoms(product_xyz)
    n = len(r_atoms)
    if len(p_atoms) != n:
        return None

    best_heavy: Optional[Tuple[int, int]] = None
    best_heavy_delta = 0.0
    best_any: Optional[Tuple[int, int]] = None
    best_any_delta = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            d_r = float(r_atoms.get_distance(i, j))
            d_p = float(p_atoms.get_distance(i, j))
            delta = d_r - d_p  # positive = bond shortened
            if d_p < max_product_distance and delta > 0.3:
                both_heavy = (
                    r_atoms[i].symbol != "H" and r_atoms[j].symbol != "H"
                )
                if both_heavy and delta > best_heavy_delta:
                    best_heavy_delta = delta
                    best_heavy = (i, j)
                if delta > best_any_delta:
                    best_any_delta = delta
                    best_any = (i, j)

    # Prefer heavy-heavy pairs; fall back to any pair (e.g. pure H-transfer).
    best_pair = best_heavy if best_heavy is not None else best_any
    best_delta = best_heavy_delta if best_heavy is not None else best_any_delta

    if best_pair is not None:
        logger.info(
            f"    [FORMING-BOND] Detected pair ({best_pair[0]}, {best_pair[1]}) "
            f"Δd={best_delta:.2f} Å"
        )
        return best_pair
    return None


def run_xtb_relaxed_scan(
    reactant_xyz: str,
    product_xyz: str,
    charge: int = 0,
    spin: int = 0,
    *,
    forming_bond: Optional[Tuple[int, int]] = None,
    n_points: int = 8,
    fmax: float = 0.10,
    max_steps_per_point: int = 50,
) -> Optional[Dict[str, Any]]:
    """Relaxed scan along the forming bond at GFN2-xTB level.

    Returns dict with ``xyz`` (best TS candidate), ``energy_ev``,
    ``fmax_ev_ang``, ``scan_distance``, ``scan_points``, or None on failure.
    """
    from ase.constraints import FixBondLength
    from ase.optimize import BFGS

    if forming_bond is None:
        forming_bond = detect_forming_bond(reactant_xyz, product_xyz)
    if forming_bond is None:
        logger.warning("    [RELAXED-SCAN] Could not detect forming bond; skipping.")
        return None

    i_atom, j_atom = forming_bond
    r_atoms = _xyz_to_ase_atoms(reactant_xyz)
    p_atoms = _xyz_to_ase_atoms(product_xyz)

    d_reactant = float(r_atoms.get_distance(i_atom, j_atom))
    d_product = float(p_atoms.get_distance(i_atom, j_atom))

    # Scan from reactant distance toward product distance,
    # but never go below a physically reasonable minimum bond length.
    _MIN_BOND_DISTANCE = 1.2  # Å — shorter than any covalent single bond
    scan_start = d_reactant
    scan_stop = max(d_product, _MIN_BOND_DISTANCE)
    if abs(scan_start - scan_stop) < 0.2:
        logger.warning("    [RELAXED-SCAN] Start/stop too close; skipping.")
        return None

    step = (scan_stop - scan_start) / n_points
    distances = [scan_start + step * k for k in range(1, n_points + 1)]

    calc = get_xtb_ase_calculator(charge=charge, spin=spin, solvent="water")

    best_atoms = None
    best_energy = None
    best_distance = None
    scan_log: List[Dict[str, float]] = []

    logger.info(
        f"    [RELAXED-SCAN] Scanning bond ({i_atom},{j_atom}) "
        f"from {scan_start:.2f} to {scan_stop:.2f} Å in {n_points} steps"
    )

    for target_d in distances:
        try:
            atoms = r_atoms.copy()
            atoms.calc = calc
            atoms.set_distance(i_atom, j_atom, target_d, fix=0)
            atoms.set_constraint(FixBondLength(i_atom, j_atom))
            with io.StringIO() as buf:
                opt = BFGS(atoms, logfile=buf)
                opt.run(fmax=fmax, steps=max_steps_per_point)
            energy = float(atoms.get_potential_energy())
            f = atoms.get_forces()
            fmax_val = float(np.sqrt((f ** 2).sum(axis=1).max()))
            scan_log.append({"distance": target_d, "energy_ev": energy, "fmax": fmax_val})
            if best_energy is None or energy > best_energy:
                best_atoms = atoms.copy()
                best_energy = energy
                best_distance = target_d
        except Exception as exc:
            logger.warning(f"    [RELAXED-SCAN] Point d={target_d:.2f} failed: {exc}")
            continue

    if best_atoms is None:
        logger.warning("    [RELAXED-SCAN] All scan points failed.")
        return None

    # Prefer an interior maximum (true saddle region) over edge points
    # which may just be compression artifacts.
    if len(scan_log) >= 3:
        for k in range(1, len(scan_log) - 1):
            if (scan_log[k]["energy_ev"] >= scan_log[k - 1]["energy_ev"]
                    and scan_log[k]["energy_ev"] >= scan_log[k + 1]["energy_ev"]):
                # Re-locate best from the interior peak
                if best_distance == distances[-1] or best_distance == distances[0]:
                    # Current best is at an edge — prefer interior peak
                    best_distance = scan_log[k]["distance"]
                    # Re-create atoms at the interior peak distance
                    _peak_atoms = r_atoms.copy()
                    _peak_atoms.calc = get_xtb_ase_calculator(
                        charge=charge, spin=spin, solvent="water"
                    )
                    _peak_atoms.set_distance(i_atom, j_atom, best_distance, fix=0)
                    from ase.constraints import FixBondLength as _FBL
                    _peak_atoms.set_constraint(_FBL(i_atom, j_atom))
                    try:
                        from ase.optimize import BFGS as _BFGS
                        _opt = _BFGS(_peak_atoms, logfile=None)
                        _opt.run(fmax=fmax, steps=max_steps_per_point)
                        _peak_atoms.set_constraint()
                        best_atoms = _peak_atoms.copy()
                        best_energy = float(_peak_atoms.get_potential_energy())
                        logger.info(
                            f"    [RELAXED-SCAN] Preferring interior peak at "
                            f"d={best_distance:.2f} Å over edge maximum"
                        )
                    except Exception:
                        pass  # fall back to original best
                break

    # Remove constraint for final fmax measurement
    best_atoms.set_constraint()
    best_atoms.calc = calc
    forces = best_atoms.get_forces()
    final_fmax = float(np.sqrt((forces ** 2).sum(axis=1).max()))

    logger.info(
        f"    [RELAXED-SCAN] Best at d={best_distance:.2f} Å, "
        f"E={best_energy:.4f} eV, fmax={final_fmax:.4f} eV/Å"
    )
    return {
        "xyz": _ase_atoms_to_xyz(best_atoms),
        "energy_ev": best_energy,
        "fmax_ev_ang": final_fmax,
        "scan_distance": best_distance,
        "scan_points": scan_log,
        "forming_bond": list(forming_bond),
    }


def run_xtb_cineb(
    reactant_xyz: str,
    product_xyz: str,
    charge: int = 0,
    spin: int = 0,
    *,
    path_frames_dir: Optional[str] = None,
    n_images: int = 8,
    fmax: float = 0.10,
    max_steps: int = 200,
    spring_constant: float = 0.1,
) -> Optional[Dict[str, Any]]:
    """CI-NEB at GFN2-xTB level, seeded with existing path frames if available.

    Returns dict with ``xyz`` (climbing-image geometry), ``energy_ev``,
    ``fmax_ev_ang``, or None on failure.
    """
    from ase.io import read as ase_read
    from ase.mep.neb import NEB
    from ase.optimize import FIRE

    r_atoms = _xyz_to_ase_atoms(reactant_xyz)
    p_atoms = _xyz_to_ase_atoms(product_xyz)

    # Try to seed from existing xtbpath frames
    seed_frames: List = []
    if path_frames_dir and os.path.isdir(path_frames_dir):
        frame_files = sorted(
            [f for f in os.listdir(path_frames_dir) if f.startswith("xtbpath_") and f.endswith(".xyz") and f != "xtbpath_ts.xyz"],
            key=lambda f: int(f.replace("xtbpath_", "").replace(".xyz", "")) if f.replace("xtbpath_", "").replace(".xyz", "").isdigit() else 999,
        )
        for ff in frame_files:
            try:
                a = ase_read(os.path.join(path_frames_dir, ff), format="xyz")
                seed_frames.append(a)
            except Exception:
                continue

    # Build images
    def _make_calc():
        return get_xtb_ase_calculator(charge=charge, spin=spin, solvent="water")

    # Endpoints need calculators + cached energies for NEB spring/climb forces
    r_img = r_atoms.copy()
    r_img.calc = _make_calc()
    r_img.get_potential_energy()

    p_img = p_atoms.copy()
    p_img.calc = _make_calc()
    p_img.get_potential_energy()

    if len(seed_frames) >= 3:
        # Subsample to n_images intermediate frames
        indices = np.linspace(0, len(seed_frames) - 1, n_images, dtype=int)
        images = [r_img]
        for idx in indices:
            img = seed_frames[idx].copy()
            img.calc = _make_calc()
            images.append(img)
        images.append(p_img)
        logger.info(f"    [CI-NEB] Seeded with {len(seed_frames)} path frames → {n_images} images")
    else:
        # IDPP interpolation — only for single-fragment systems
        from ase.mep.neb import interpolate as neb_interpolate
        images = [r_img]
        for _ in range(n_images):
            img = r_atoms.copy()
            img.calc = _make_calc()
            images.append(img)
        images.append(p_img)
        neb_interpolate(images)
        logger.info(f"    [CI-NEB] IDPP interpolation with {n_images} images (no path frames)")

    neb = NEB(images, climb=True, k=spring_constant, method="improvedtangent")

    try:
        opt = FIRE(neb, logfile=None)
        opt.run(fmax=fmax, steps=max_steps)
    except Exception as exc:
        logger.warning(f"    [CI-NEB] NEB optimization failed: {exc}")
        return None

    # Find climbing image (highest energy intermediate)
    best_img = None
    best_e = None
    for img in images[1:-1]:
        try:
            e = float(img.get_potential_energy())
            if best_e is None or e > best_e:
                best_e = e
                best_img = img
        except Exception:
            continue

    if best_img is None:
        logger.warning("    [CI-NEB] No valid images after NEB.")
        return None

    forces = best_img.get_forces()
    final_fmax = float(np.sqrt((forces ** 2).sum(axis=1).max()))
    logger.info(f"    [CI-NEB] Best image E={best_e:.4f} eV, fmax={final_fmax:.4f} eV/Å")

    return {
        "xyz": _ase_atoms_to_xyz(best_img),
        "energy_ev": best_e,
        "fmax_ev_ang": final_fmax,
    }


# ---------------------------------------------------------------------------
# xTB-level Sella TS refinement
# ---------------------------------------------------------------------------

def refine_ts_sella_xtb(
    ts_guess_xyz: str,
    charge: int = 0,
    spin: int = 0,
    *,
    fmax: float = 0.05,
    max_steps: int = 150,
    timeout_seconds: float = 300.0,
) -> Optional[Dict[str, Any]]:
    """Run Sella eigenvector-following TS optimisation at GFN2-xTB level.

    Takes a TS guess (e.g. from relaxed scan or NEB) and refines it to an
    actual saddle point at xTB level.  Returns dict with ``xyz``,
    ``energy_ev``, ``fmax_ev_ang``, or *None* on failure.
    """
    import time as _time

    try:
        from sella import Sella as _Sella
    except ImportError:
        logger.warning("    [xTB-SELLA] Sella not installed; skipping xTB TS refinement.")
        return None

    atoms = _xyz_to_ase_atoms(ts_guess_xyz)
    calc = get_xtb_ase_calculator(charge=charge, spin=spin)
    if calc is None:
        logger.warning("    [xTB-SELLA] No xTB calculator available; skipping.")
        return None
    atoms.calc = calc

    logger.info(f"    [xTB-SELLA] Starting Sella TS search (fmax={fmax}, max_steps={max_steps})...")

    try:
        dyn = _Sella(atoms, logfile=None, trajectory=None)
        start = _time.monotonic()
        best_fmax = float("inf")
        best_xyz: Optional[str] = None
        best_energy = 0.0

        def _watchdog():
            nonlocal best_fmax, best_xyz, best_energy
            elapsed = _time.monotonic() - start
            if elapsed > timeout_seconds:
                raise RuntimeError(f"xTB-Sella exceeded {timeout_seconds:.0f}s")
            cached = getattr(atoms.calc, "results", {}).get("forces")
            if cached is not None:
                fmax_now = float(np.sqrt((cached ** 2).sum(axis=1).max()))
                if fmax_now < best_fmax:
                    best_fmax = fmax_now
                    best_xyz = _ase_atoms_to_xyz(atoms)
                    best_energy = float(atoms.get_potential_energy())

        dyn.attach(_watchdog)
        dyn.run(fmax=fmax, steps=max_steps)
    except Exception as exc:
        logger.warning(f"    [xTB-SELLA] Sella failed: {exc}")
        # Fall through — use best geometry found so far if any
        if best_xyz is None:
            return None

    # Final measurement — wrapped because atoms may be in broken SCF state
    try:
        if best_xyz is None:
            best_xyz = _ase_atoms_to_xyz(atoms)
            best_energy = float(atoms.get_potential_energy())
        forces = atoms.get_forces()
        final_fmax = float(np.sqrt((forces ** 2).sum(axis=1).max()))
        if final_fmax < best_fmax:
            best_fmax = final_fmax
            best_xyz = _ase_atoms_to_xyz(atoms)
            best_energy = float(atoms.get_potential_energy())
    except Exception:
        if best_xyz is None:
            return None

    logger.info(
        f"    [xTB-SELLA] Done: fmax={best_fmax:.4f} eV/Å, E={best_energy:.4f} eV"
    )
    return {
        "xyz": best_xyz,
        "energy_ev": best_energy,
        "fmax_ev_ang": best_fmax,
    }
