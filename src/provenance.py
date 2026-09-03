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
from typing import Any, Dict, Iterable, List, Optional

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
VOLATILE_KEYS = ("generated_on", "git", "wall_seconds", "workers")


__all__ = [
    "VOLATILE_KEYS",
    "git_head",
    "input_sources",
    "provenance_block",
    "sha256_of",
    "today_iso",
]
