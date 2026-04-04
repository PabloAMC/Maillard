#!/usr/bin/env python3
"""
scripts/run_dft_c4_c5.py — Step C4 + C5 DFT Runner.

Executes the two-stage DFT refinement for the 6 highest-decision-impact
reaction targets identified in the L4 plan:

  Step C4: r2SCAN-3c geometry optimisation + frequency calculation (via DFTRefiner)
  Step C5: wB97M-V/def2-tzvp single-point on the C4-optimised geometry

Produces per-reaction JSON results in results/dft_c4_c5/<reaction_key>/.

Usage:
  # Run all 6 targets sequentially
  python scripts/run_dft_c4_c5.py

  # Run a single target
  python scripts/run_dft_c4_c5.py --reaction hexanal_radical_quench

  # Dry-run: show job manifest without running DFT
  python scripts/run_dft_c4_c5.py --dry-run

The script is safe to re-run: it skips any target whose output JSON already exists
(unless --force is passed).
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("dft_c4_c5")

# ── Job manifest ──────────────────────────────────────────────────────────────
# 6 highest-decision-impact targets from the L4 plan.
# charge=0, spin=0 for all closed-shell ground states.
# pe_schiff_base and pe_amadori have well-known literature Ea values that serve
# as acceptance anchors for the DFT result.
C4_C5_JOBS: List[Dict[str, Any]] = [
    {
        "reaction_key": "hexanal_radical_quench",
        "slr_family": "11",
        "charge": 0,
        "spin": 0,
        "description": "Cys + hexanal radical quenching — Maillard/lipid competition surface",
        "literature_ea_kcal": None,          # No firm literature prior; use DFT
        "mace_barrier_kcal": 7.89,
        "promoted_tier": "xtb_derived_gfn2", # Will promote to selective_dft_anchor after DFT
        "uncertainty_kj": 42.0,              # ±10 kJ/mol → ±10 kcal/mol × 4.184 floor
        "promotion_ceiling": "ranking_only",
    },
    {
        "reaction_key": "quinone_cys_michael",
        "slr_family": "13",
        "charge": 0,
        "spin": 0,
        "description": "Ortho-quinone + Cys Michael addition — polyphenol precursor sink",
        "literature_ea_kcal": None,
        "mace_barrier_kcal": 20.33,
        "promoted_tier": "xtb_derived_gfn2",
        "uncertainty_kj": 42.0,
        "promotion_ceiling": "bounded_calibration",
    },
    {
        "reaction_key": "aa_ring_open_dicarbonyl",
        "slr_family": "14",
        "charge": 0,
        "spin": 0,
        "description": "Ascorbic acid ring-open → 3-deoxyosone + dicarbonyl source",
        "literature_ea_kcal": 31.7 / 4.184,  # 31.7 kJ/mol from HCW (Ea=7.577 kcal/mol)
        "mace_barrier_kcal": 18.10,
        "promoted_tier": "xtb_derived_gfn2",  # Literature Ea already available; DFT validates
        "uncertainty_kj": 20.0,              # ±5 kJ/mol target for r2SCAN-3c
        "promotion_ceiling": "bounded_calibration",
    },
    {
        "reaction_key": "pe_schiff_base",
        "slr_family": "15",
        "charge": 0,
        "spin": 0,
        "description": "PE headgroup + glucose Schiff base — interfacial sugar diversion",
        "literature_ea_kcal": 92.9 / 4.184,  # 92.9 kJ/mol from SLR 15 (Ea=22.21 kcal/mol)
        "mace_barrier_kcal": 15.20,
        "promoted_tier": "literature_derived_transfer",  # Literature better than MACE here
        "uncertainty_kj": 20.9,              # ±5 kcal/mol on transfer
        "promotion_ceiling": "bounded_calibration",
    },
    {
        "reaction_key": "pe_amadori",
        "slr_family": "15",
        "charge": 0,
        "spin": 0,
        "description": "PE Schiff → Amadori rearrangement — second step of PE-sugar sink",
        "literature_ea_kcal": 82.9 / 4.184,  # 82.9 kJ/mol from SLR 15 (Ea=19.81 kcal/mol)
        "mace_barrier_kcal": 8.94,
        "promoted_tier": "literature_derived_transfer",  # Literature better than MACE here
        "uncertainty_kj": 20.9,
        "promotion_ceiling": "bounded_calibration",
    },
    {
        "reaction_key": "lysinoalanine_crosslink",
        "slr_family": "12",
        "charge": 0,
        "spin": 0,
        "description": "DHA + Lys ε-amine → lysinoalanine crosslink (β-elimination product)",
        "literature_ea_kcal": None,
        "mace_barrier_kcal": 31.75,
        "promoted_tier": "xtb_derived_gfn2",
        "uncertainty_kj": 42.0,
        "promotion_ceiling": "ranking_only",
    },
]

OUTPUT_DIR = ROOT / "results" / "dft_c4_c5"
XTB_INPUTS = ROOT / "data" / "geometries" / "xtb_inputs"


def _load_xyz(reaction_key: str, which: str) -> str:
    """Load reactant or product XYZ from the xtb_inputs directory."""
    xyz_path = XTB_INPUTS / reaction_key / f"{which}.xyz"
    if not xyz_path.exists():
        raise FileNotFoundError(f"XYZ not found: {xyz_path}")
    return xyz_path.read_text()


def _find_ts_xyz(reaction_key: str) -> Optional[str]:
    """Return TS guess XYZ from xTB path, or None if unavailable."""
    rxn_dir = XTB_INPUTS / reaction_key
    # Prefer explicit TS file
    for candidate in ["xtbpath_ts.xyz", "ts.xyz"]:
        p = rxn_dir / candidate
        if p.exists():
            return p.read_text()
    # Fall back: use midpoint of numbered xtbpath_*.xyz files
    path_files = sorted(rxn_dir.glob("xtbpath_*.xyz"))
    if path_files:
        mid = path_files[len(path_files) // 2]
        logger.warning(
            f"  No explicit TS for {reaction_key}; using midpoint path file: {mid.name}"
        )
        return mid.read_text()
    logger.warning(f"  No TS geometry available for {reaction_key}; using reactant as TS guess")
    return None


def run_dft_job(job: Dict[str, Any], *, force: bool = False) -> Dict[str, Any]:
    """
    Execute C4 (r2SCAN-3c geom+freq) and C5 (wB97M-V SP) for one reaction.
    Returns a result dict. Writes JSON to OUTPUT_DIR/<reaction_key>/.
    """
    from src.dft_refiner import DFTRefiner

    key = job["reaction_key"]
    out_dir = OUTPUT_DIR / key
    out_dir.mkdir(parents=True, exist_ok=True)
    result_file = out_dir / "barrier_summary.json"

    if result_file.exists() and not force:
        logger.info(f"  [{key}] Already done — loading cached result")
        return json.loads(result_file.read_text())

    logger.info(f"  [{key}] Loading geometries ...")
    reactant_xyz = _load_xyz(key, "reactant")
    ts_xyz = _find_ts_xyz(key)

    refiner = DFTRefiner(
        solvent_name="water",
        temp_k=423.15,  # 150 °C extrusion temperature
        geometry_backend="pyscf",
    )

    t0 = time.time()

    # ── Step C4: r2SCAN-3c geometry + frequency ───────────────────────────
    logger.info(f"  [{key}] Step C4: r2SCAN-3c geometry optimisation ...")
    res_r = refiner.optimize_geometry(reactant_xyz, charge=job["charge"], spin=job["spin"], is_ts=False)
    if not res_r.converged:
        raise RuntimeError(f"[{key}] Reactant geometry did not converge")

    ts_use = ts_xyz if ts_xyz else reactant_xyz
    res_ts = refiner.optimize_geometry(ts_use, charge=job["charge"], spin=job["spin"], is_ts=True)
    if not res_ts.converged:
        raise RuntimeError(f"[{key}] TS geometry did not converge")

    # ── Step C5: wB97M-V single-point on C4 geometry ─────────────────────
    logger.info(f"  [{key}] Step C5: wB97M-V single-point energy ...")
    sp_r   = refiner.single_point(res_r.optimized_xyz,  xc_method="wB97M-V", basis="def2-tzvp",
                                   charge=job["charge"], spin=job["spin"])
    sp_ts  = refiner.single_point(res_ts.optimized_xyz, xc_method="wB97M-V", basis="def2-tzvp",
                                   charge=job["charge"], spin=job["spin"])

    # ── Compute ΔG‡ ───────────────────────────────────────────────────────
    # Thermal correction: quasi-harmonic G from C4, SP energy from C5
    assert res_ts.quasi_harmonic_gibbs_hartree is not None
    assert res_r.quasi_harmonic_gibbs_hartree is not None
    therm_corr_ts = res_ts.quasi_harmonic_gibbs_hartree - res_ts.energy_hartree  # thermal correction
    therm_corr_r  = res_r.quasi_harmonic_gibbs_hartree  - res_r.energy_hartree

    refined_g_ts = sp_ts + therm_corr_ts
    refined_g_r  = sp_r  + therm_corr_r
    delta_g_h    = refined_g_ts - refined_g_r
    barrier_kcal = delta_g_h * 627.509

    elapsed = time.time() - t0
    logger.info(f"  [{key}] ΔG‡ = {barrier_kcal:.3f} kcal/mol  ({elapsed:.0f}s)")

    result = {
        "reaction_key": key,
        "slr_family": job["slr_family"],
        "method": "wB97M-V//r2SCAN-3c",
        "basis_geom": "def2-svp",
        "basis_sp":   "def2-tzvp",
        "solvent": "water (ddCOSMO)",
        "temp_k": 423.15,
        "barrier_kcal_mol": round(barrier_kcal, 4),
        "barrier_kj_mol":   round(barrier_kcal * 4.184, 4),
        "uncertainty_kj":   job["uncertainty_kj"],
        "ppb_low_factor":   round(1.0 / (2.0 ** (job["uncertainty_kj"] / 4.184 / 1.987e-3 / 423.15 / 2.303)), 4),
        "source_quality":   "selective_dft_anchor",
        "promotion_ceiling": job["promotion_ceiling"],
        "honest_label":     "DFT-based estimate (wB97M-V//r2SCAN-3c), ±factor 2 in absolute concentration, ranking reliable",
        "optimized_reactant_xyz": res_r.optimized_xyz,
        "optimized_ts_xyz":       res_ts.optimized_xyz,
        "quasi_harmonic_gibbs_reactant_h": res_r.quasi_harmonic_gibbs_hartree,
        "quasi_harmonic_gibbs_ts_h":       res_ts.quasi_harmonic_gibbs_hartree,
        "sp_energy_reactant_h": sp_r,
        "sp_energy_ts_h":       sp_ts,
        "elapsed_s": round(elapsed, 1),
        "mace_barrier_kcal_mol": job["mace_barrier_kcal"],
        "literature_ea_kcal_mol": job.get("literature_ea_kcal"),
    }

    # Write optimised geometries
    (out_dir / "optimized_reactant.xyz").write_text(res_r.optimized_xyz)
    (out_dir / "optimized_ts.xyz").write_text(res_ts.optimized_xyz)
    result_file.write_text(json.dumps(result, indent=2))
    logger.info(f"  [{key}] Written to {result_file}")
    return result


def write_job_manifest() -> None:
    """Write the machine-readable C4/C5 job manifest to results/dft_c4_c5/."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest_path = OUTPUT_DIR / "c4_c5_job_manifest.json"
    manifest = {
        "description": "L4 DFT job manifest: Steps C4 (r2SCAN-3c geom+freq) + C5 (wB97M-V SP)",
        "method_chain": "wB97M-V/def2-tzvp // r2SCAN-3c/def2-svp + ddCOSMO(water)",
        "temperature_k": 423.15,
        "jobs": C4_C5_JOBS,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2))
    logger.info(f"Job manifest written to {manifest_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Step C4+C5 DFT runner")
    parser.add_argument("--reaction", default=None, help="Run only this reaction key")
    parser.add_argument("--dry-run", action="store_true", help="Print manifest and exit")
    parser.add_argument("--force", action="store_true", help="Re-run even if output exists")
    args = parser.parse_args()

    write_job_manifest()

    jobs = C4_C5_JOBS
    if args.reaction:
        jobs = [j for j in jobs if j["reaction_key"] == args.reaction]
        if not jobs:
            logger.error(f"Unknown reaction key: {args.reaction}")
            sys.exit(1)

    if args.dry_run:
        print(json.dumps({"jobs": jobs}, indent=2))
        return

    results = []
    failures = []
    for job in jobs:
        try:
            result = run_dft_job(job, force=args.force)
            results.append(result)
        except Exception as exc:
            logger.error(f"  [{job['reaction_key']}] FAILED: {exc}")
            failures.append({"reaction_key": job["reaction_key"], "error": str(exc)})

    # Write summary
    summary_path = OUTPUT_DIR / "c4_c5_run_summary.json"
    summary_path.write_text(json.dumps({"results": results, "failures": failures}, indent=2))
    logger.info(f"Summary: {len(results)} completed, {len(failures)} failed — {summary_path}")

    if failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
