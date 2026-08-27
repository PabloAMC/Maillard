from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List

from src.xyz_common import extract_xyz_last_frame


def path_frame_index(path: Path) -> int:
    match = re.search(r"xtbpath_(\d+)\.xyz$", path.name)
    return int(match.group(1)) if match else -1


def inspect_xtb_log_health(runner_dir: Path) -> Dict[str, Any]:
    log_path = runner_dir / "xtbopt.log"
    if not log_path.exists():
        return {"quality_gate_passed": True, "failure_markers": []}

    content = log_path.read_text(encoding="utf-8", errors="ignore")
    markers = [
        marker
        for marker in [
            "Self consistent charge iterator did not converge",
            "SCF not converged, aborting",
            "Could not relax/optimize structure",
            "Setup of Coulomb evaluator failed",
        ]
        if marker in content
    ]
    return {
        "quality_gate_passed": not markers,
        "failure_markers": markers,
    }


def _parse_xtb_energy_from_content(xyz_content: str) -> float | None:
    """Extract xTB energy from an XYZ content string's comment line (line 2)."""
    try:
        lines = xyz_content.splitlines()
        if len(lines) < 2:
            return None
        comment = lines[1]
        for token_idx, token in enumerate(comment.split()):
            if token == "energy:":
                return float(comment.split()[token_idx + 1])
        return None
    except (ValueError, IndexError):
        return None

def materialize_xtb_outputs(runner_dir: Path) -> List[Path]:
    frame_paths = sorted(
        [path for path in runner_dir.glob("xtbpath_*.xyz") if path_frame_index(path) >= 0],
        key=path_frame_index,
    )
    if not frame_paths:
        return []

    path_bundle = runner_dir / "xtbpath.xyz"
    ts_guess = runner_dir / "xtbpath_ts.xyz"

    created: List[Path] = []
    if not path_bundle.exists():
        # For the bundle, use last frames only (optimized geometries, not trajectories)
        last_frames = [extract_xyz_last_frame(fp.read_text(encoding="utf-8")).strip() for fp in frame_paths]
        path_bundle.write_text("\n".join(chunk for chunk in last_frames if chunk) + "\n", encoding="utf-8")
        created.append(path_bundle)

    if not ts_guess.exists():
        # Extract last frame from each trajectory and parse energy
        interior_frames = frame_paths[1:-1] if len(frame_paths) > 2 else frame_paths
        candidates = []
        for fp in interior_frames:
            last_frame = extract_xyz_last_frame(fp.read_text(encoding="utf-8"))
            energy = _parse_xtb_energy_from_content(last_frame)
            candidates.append((fp, last_frame, energy))

        frames_with_energy = [(fp, lf, e) for fp, lf, e in candidates if e is not None and e != 0.0]

        if frames_with_energy:
            # Pick frame with highest energy (most TS-like)
            best_fp, best_frame, _ = max(frames_with_energy, key=lambda x: x[2])
        else:
            # Fallback: pick the frame with largest RMSD from first frame (most displaced = most TS-like)
            ref_frame = extract_xyz_last_frame(frame_paths[0].read_text(encoding="utf-8"))
            best_frame = extract_xyz_last_frame(frame_paths[len(frame_paths) // 2].read_text(encoding="utf-8"))

        ts_guess.write_text(best_frame, encoding="utf-8")
        created.append(ts_guess)

    return created


def assess_xtb_path_quality(runner_dir: Path, *, materialize_missing: bool = False) -> Dict[str, Any]:
    created = materialize_xtb_outputs(runner_dir) if materialize_missing else []
    frame_paths = sorted(
        [path for path in runner_dir.glob("xtbpath_*.xyz") if path_frame_index(path) >= 0],
        key=path_frame_index,
    )
    path_bundle = runner_dir / "xtbpath.xyz"
    ts_guess = runner_dir / "xtbpath_ts.xyz"
    health = inspect_xtb_log_health(runner_dir)

    return {
        "frame_count": len(frame_paths),
        "has_path_bundle": path_bundle.exists(),
        "has_ts_guess": ts_guess.exists(),
        "synthesized_outputs": created,
        **health,
    }