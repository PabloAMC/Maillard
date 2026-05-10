"""Compare React-OT TS seeds against existing xTB TS guesses.

For each CHON-eligible target, this report computes:

* atom-count consistency (React-OT seed vs reactant)
* RMSD (Kabsch-aligned) between the React-OT seed and the reactant
  (sanity: a meaningful TS guess must displace from the reactant by
  >= 0.3 Å — see the TS-guess gates documented in agents.md)
* RMSD between the React-OT seed and the existing xTB TS guess
  (when present)

The output is *evidence only*: it is the artifact required by the
project's observable-first governance rule before React-OT seeds may be
promoted into the Sella DFT seed pipeline. React-OT remains a
geometry-only pathfinder — never a barrier authority.

Outputs:
    results/validation/react_ot_seed_coverage.json
    results/validation/react_ot_seed_coverage.md
"""
from __future__ import annotations

import argparse
import io
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.recover_ts_react_ot_seed import ELIGIBLE_TARGETS  # noqa: E402

DEFAULT_SEED_DIR = REPO_ROOT / "results" / "computational_gap_refinement"
DEFAULT_XTB_INPUTS = REPO_ROOT / "data" / "geometries" / "xtb_inputs"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results" / "validation"

# TS-guess gate from agents.md: minimum pairwise RMS Δ between
# reactant and TS must be at least 0.3 Å.
MIN_REACTANT_TS_RMSD_ANGSTROM = 0.30


def _read_xyz_atoms(path: Path):
    from ase.io import read

    with path.open("r", encoding="utf-8") as handle:
        atoms = read(handle, format="xyz")
    if isinstance(atoms, list):
        atoms = atoms[-1]
    return atoms


def _aligned_rmsd(reference, mobile) -> float:
    import numpy as np
    from ase.build import minimize_rotation_and_translation

    if len(reference) != len(mobile):
        raise ValueError("atom count mismatch")
    minimize_rotation_and_translation(reference, mobile)
    delta = reference.get_positions() - mobile.get_positions()
    return float(np.sqrt(np.mean(np.sum(delta ** 2, axis=1))))


def _symbols(atoms) -> Tuple[str, ...]:
    return tuple(atoms.get_chemical_symbols())


def _evaluate_target(target: str, seed_dir: Path, xtb_inputs: Path) -> Dict[str, Any]:
    seed_xyz = seed_dir / f"{target}_react_ot_seed.xyz"
    reactant_xyz = xtb_inputs / target / "reactant.xyz"
    xtb_ts_xyz = xtb_inputs / target / "xtbpath_ts.xyz"

    record: Dict[str, Any] = {
        "target": target,
        "seed_path": str(seed_xyz.relative_to(REPO_ROOT)),
        "reactant_path": str(reactant_xyz.relative_to(REPO_ROOT)),
        "xtb_ts_path": str(xtb_ts_xyz.relative_to(REPO_ROOT)),
        "seed_present": seed_xyz.exists(),
        "reactant_present": reactant_xyz.exists(),
        "xtb_ts_present": xtb_ts_xyz.exists(),
    }

    if not seed_xyz.exists():
        record["status"] = "missing_seed"
        return record
    if not reactant_xyz.exists():
        record["status"] = "missing_reactant"
        return record

    try:
        seed_atoms = _read_xyz_atoms(seed_xyz)
        reactant_atoms = _read_xyz_atoms(reactant_xyz)
    except Exception as exc:  # pragma: no cover - ASE errors propagated
        record["status"] = "read_error"
        record["error"] = str(exc)
        return record

    record["seed_atom_count"] = len(seed_atoms)
    record["reactant_atom_count"] = len(reactant_atoms)
    record["atom_count_match"] = len(seed_atoms) == len(reactant_atoms)
    record["symbol_sequence_match"] = _symbols(seed_atoms) == _symbols(reactant_atoms)

    if not record["atom_count_match"]:
        record["status"] = "atom_count_mismatch"
        return record

    try:
        rmsd_reactant = _aligned_rmsd(reactant_atoms, seed_atoms)
    except Exception as exc:  # pragma: no cover
        record["status"] = "rmsd_error"
        record["error"] = str(exc)
        return record
    record["rmsd_vs_reactant_angstrom"] = rmsd_reactant
    record["rmsd_vs_reactant_passes_min_gate"] = rmsd_reactant >= MIN_REACTANT_TS_RMSD_ANGSTROM

    rmsd_xtb: Optional[float] = None
    if xtb_ts_xyz.exists():
        try:
            xtb_atoms = _read_xyz_atoms(xtb_ts_xyz)
            if len(xtb_atoms) == len(seed_atoms):
                rmsd_xtb = _aligned_rmsd(xtb_atoms, seed_atoms)
            else:
                record["xtb_atom_count_mismatch"] = True
        except Exception as exc:  # pragma: no cover
            record["xtb_compare_error"] = str(exc)
    record["rmsd_vs_xtb_ts_angstrom"] = rmsd_xtb

    if record["atom_count_match"] and record["rmsd_vs_reactant_passes_min_gate"]:
        record["status"] = "ok"
    else:
        record["status"] = "warn"
    return record


def build_report(seed_dir: Path, xtb_inputs: Path, targets: List[str]) -> Dict[str, Any]:
    records = [_evaluate_target(t, seed_dir, xtb_inputs) for t in targets]
    summary = {
        "n_targets": len(records),
        "n_ok": sum(1 for r in records if r.get("status") == "ok"),
        "n_warn": sum(1 for r in records if r.get("status") == "warn"),
        "n_missing_seed": sum(1 for r in records if r.get("status") == "missing_seed"),
        "n_other": sum(
            1
            for r in records
            if r.get("status") not in {"ok", "warn", "missing_seed"}
        ),
    }
    return {
        "report": "react_ot_seed_coverage",
        "trust_posture": {
            "role": "ts_initialization_geometry_only",
            "is_runtime_authority": False,
            "energy_use_allowed": False,
        },
        "min_reactant_ts_rmsd_angstrom": MIN_REACTANT_TS_RMSD_ANGSTROM,
        "summary": summary,
        "targets": records,
    }


def render_markdown(report: Dict[str, Any]) -> str:
    lines: List[str] = []
    lines.append("# React-OT seed coverage report")
    lines.append("")
    lines.append(
        "_React-OT seeds are geometry-only pathfinders. Sella DFT remains the "
        "single barrier authority._"
    )
    lines.append("")
    s = report["summary"]
    lines.append(
        f"**Summary**: {s['n_ok']} ok / {s['n_warn']} warn / "
        f"{s['n_missing_seed']} missing seeds / {s['n_other']} other "
        f"(of {s['n_targets']} eligible targets)."
    )
    lines.append("")
    lines.append(
        "| Target | Status | Atoms (seed/reactant) | RMSD vs reactant (Å) | "
        "RMSD vs xTB TS (Å) |"
    )
    lines.append("| --- | --- | --- | --- | --- |")
    for r in report["targets"]:
        atoms_cell = (
            f"{r.get('seed_atom_count', '–')}/{r.get('reactant_atom_count', '–')}"
        )
        rmsd_r = r.get("rmsd_vs_reactant_angstrom")
        rmsd_x = r.get("rmsd_vs_xtb_ts_angstrom")
        rmsd_r_cell = f"{rmsd_r:.3f}" if isinstance(rmsd_r, float) else "–"
        rmsd_x_cell = f"{rmsd_x:.3f}" if isinstance(rmsd_x, float) else "–"
        lines.append(
            f"| `{r['target']}` | {r['status']} | {atoms_cell} | "
            f"{rmsd_r_cell} | {rmsd_x_cell} |"
        )
    lines.append("")
    lines.append(
        f"Gate: `rmsd_vs_reactant >= {MIN_REACTANT_TS_RMSD_ANGSTROM:.2f} Å` "
        "(see agents.md TS-guess gates)."
    )
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-dir", default=str(DEFAULT_SEED_DIR))
    parser.add_argument("--xtb-inputs", default=str(DEFAULT_XTB_INPUTS))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument(
        "--target",
        action="append",
        choices=ELIGIBLE_TARGETS,
        help="Restrict to one or more eligible targets (default: all).",
    )
    args = parser.parse_args()

    targets = args.target or list(ELIGIBLE_TARGETS)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    report = build_report(Path(args.seed_dir), Path(args.xtb_inputs), targets)
    json_path = out_dir / "react_ot_seed_coverage.json"
    md_path = out_dir / "react_ot_seed_coverage.md"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True))
    md_path.write_text(render_markdown(report))

    print(json.dumps({"json": str(json_path), "md": str(md_path), "summary": report["summary"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
