"""One provenance block for every generated artifact.

2026-09-03 (improvement backlog, "shared report helpers"). Before this module the tree carried
six copies of a git-head helper returning three different shapes (``dirty`` as a bool in two,
as the string ``"yes"``/``"no"`` in two, absent in two), three idioms for today's date, four
inline sha256 sites, and five disjoint provenance schemas across ``results/validation``. The
three headline artifacts (scorecard, envelope, directional) carried no git stamp at all.

Everything here is stdlib-only and never raises for a missing repository: a stamp that
cannot be taken is recorded as ``"unknown"`` rather than aborting a 40-minute run.

Read ``git`` as "the checkout the generator ran in": a tracked artifact is regenerated BEFORE
the commit that carries it, so its ``commit`` is that commit's parent and ``dirty`` is
normally true. The ``inputs`` hashes, not the git stamp, are what a freshness check compares.
"""
from __future__ import annotations

import hashlib
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional

from src import data_paths


def git_head(cwd: Path | str = data_paths.REPO_ROOT) -> Dict[str, Any]:
    """``{"commit", "short", "branch", "dirty"}`` for the checkout at ``cwd``.

    ``dirty`` is a bool (uncommitted changes to tracked files). All four fields are
    ``"unknown"`` / ``False`` when git is unavailable, so callers can always format them.
    """

    def _run(*args: str) -> Optional[str]:
        try:
            return subprocess.check_output(
                ["git", *args], cwd=str(cwd), text=True, stderr=subprocess.DEVNULL, timeout=10
            ).strip()
        except (OSError, subprocess.SubprocessError):
            return None

    commit = _run("rev-parse", "HEAD")
    branch = _run("rev-parse", "--abbrev-ref", "HEAD")
    status = _run("status", "--porcelain", "--untracked-files=no")
    return {
        "commit": commit or "unknown",
        "short": (commit or "unknown")[:8],
        "branch": branch or "unknown",
        "dirty": bool(status) if status is not None else False,
    }


def today_iso() -> str:
    """Today's date as ISO-8601, in UTC (the artifacts are compared across machines)."""
    return datetime.now(timezone.utc).date().isoformat()


def sha256_of(path: Path | str) -> Optional[str]:
    """Hex digest of a file's bytes, or ``None`` when the file does not exist."""
    candidate = Path(path)
    if not candidate.is_file():
        return None
    return hashlib.sha256(candidate.read_bytes()).hexdigest()


def input_sources(paths: Iterable[Path | str]) -> List[Dict[str, Optional[str]]]:
    """``[{"path": repo-relative, "sha256": digest-or-None}, ...]`` in the order given.

    Missing inputs are listed with ``sha256: None`` rather than dropped: an artifact that
    silently stopped reading one of its inputs is exactly the drift this exists to expose.
    """
    return [{"path": data_paths.rel(p), "sha256": sha256_of(p)} for p in paths]


def provenance_block(
    artifact: str,
    *,
    generated_by: Optional[str] = None,
    wave: Optional[str] = None,
    inputs: Iterable[Path | str] = (),
) -> Dict[str, Any]:
    """The block every generated artifact carries under ``"provenance"``.

    ``generated_by`` is the repo-relative path of the writer (module or script); ``inputs``
    are the files whose content the artifact depends on (fit reports, panels, registries).
    """
    block: Dict[str, Any] = {
        "artifact": artifact,
        "generated_on": today_iso(),
        "git": git_head(),
    }
    if generated_by is not None:
        block["generated_by"] = generated_by
    if wave is not None:
        block["wave"] = wave
    block["inputs"] = input_sources(inputs)
    return block


#: Keys a freshness comparison must ignore: they change with the clock and the checkout,
#: not with the code or the data.
VOLATILE_KEYS = ("generated_on", "git", "wall_seconds", "workers", "generated")


def strip_volatile(value: Any) -> Any:
    """``value`` with every ``VOLATILE_KEYS`` entry removed, recursively (lists included)."""
    if isinstance(value, dict):
        return {k: strip_volatile(v) for k, v in value.items() if k not in VOLATILE_KEYS}
    if isinstance(value, list):
        return [strip_volatile(v) for v in value]
    return value


#: Cross-architecture floating-point reproducibility measured 2026-09-03: the same scorecard on
#: linux/arm64 vs linux/amd64 differs by at most 7e-8 relative over 49 rows. 1e-6 is a decision-safe
#: tolerance (the 3x band is six orders away) that still catches any real parameter or data change.
FLOAT_REL_TOL = 1e-6


#: Absolute floor for float comparison: differences below it are numerical noise around zero
#: (concentrations of 1e-30 ug/L, log-space steps of 1e-8).
FLOAT_ABS_TOL = 1e-9


def payload_differences(
    tracked: Any, live: Any, *, rel_tol: float = FLOAT_REL_TOL, abs_tol: float = FLOAT_ABS_TOL,
    path: str = "$", limit: int = 20,
) -> List[str]:
    """Where two artifacts differ once volatile keys are stripped; floats compare with ``rel_tol``
    or ``abs_tol`` (either suffices).

    Returns at most ``limit`` human-readable ``"$.summary.hits: 10 != 11"`` lines; an empty
    list means the artifacts match.
    """
    tracked, live = strip_volatile(tracked), strip_volatile(live)
    out: List[str] = []

    def walk(a: Any, b: Any, where: str) -> None:
        if len(out) >= limit:
            return
        if isinstance(a, dict) and isinstance(b, dict):
            for key in sorted(set(a) | set(b)):
                if key not in a or key not in b:
                    out.append(f"{where}.{key}: {'missing in tracked' if key not in a else 'missing in live'}")
                else:
                    walk(a[key], b[key], f"{where}.{key}")
            return
        if isinstance(a, list) and isinstance(b, list):
            if len(a) != len(b):
                out.append(f"{where}: length {len(a)} != {len(b)}")
                return
            for i, (x, y) in enumerate(zip(a, b)):
                walk(x, y, f"{where}[{i}]")
            return
        if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool) and not isinstance(b, bool):
            if a != b and abs(a - b) > max(rel_tol * max(abs(a), abs(b)), abs_tol):
                out.append(f"{where}: {a} != {b}")
            return
        if a != b:
            out.append(f"{where}: {a!r} != {b!r}")

    walk(tracked, live, path)
    return out


def stale_inputs(payload: Mapping[str, Any]) -> List[str]:
    """Inputs recorded in ``payload["provenance"]["inputs"]`` whose file no longer hashes the same.

    The cheap freshness test for an expensive artifact: no regeneration, just the hashes.
    Missing files and files whose hash moved are both reported.
    """
    block = payload.get("provenance") or {}
    out: List[str] = []
    for entry in block.get("inputs", []):
        current = sha256_of(data_paths.REPO_ROOT / entry["path"])
        if current != entry.get("sha256"):
            out.append(f"{entry['path']}: recorded {str(entry.get('sha256'))[:12]} now {str(current)[:12]}")
    return out


__all__ = [
    "FLOAT_ABS_TOL",
    "FLOAT_REL_TOL",
    "VOLATILE_KEYS",
    "git_head",
    "input_sources",
    "payload_differences",
    "provenance_block",
    "sha256_of",
    "stale_inputs",
    "strip_volatile",
    "today_iso",
]
