"""
src/evidence_anchors.py -- every constant's evidence anchor resolves to a file (2026-10-29, guard 4).

A ``KineticParameter`` carries its evidence as free text (``dossier_anchor``): "kinetic_core_b17_prereg.md;
kumazawa2003_extraction.md (the apparent loss ...)". Free-text citations drift (this week: a wrong Yu 2018 volume,
two matrix inputs with no dossier, a "sealed" hold-out the cited paper does not contain). The checkable part of the
text is its file tokens; this module extracts them and resolves each one to a file in the repository: a path, a
dossier name, or a uniquely named tracked file. ``scripts/ci/evidence_anchor_gate.py`` fails on a token that does not
resolve and on a measured constant that cites no dossier; ``build_sulfur_reaction_evidence.py`` reads the same tokens.
"""
from __future__ import annotations

import functools
import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

from src import data_paths

TOKEN = re.compile(r"[A-Za-z0-9_./-]+\.(?:md|json|yml|py)")
MEASURED_CLASSES = ("measured_rate", "measured_activation_energy")


def anchor_tokens(text: Optional[str]) -> List[str]:
    """The file tokens in an anchor, in order, without duplicates."""
    return list(dict.fromkeys(TOKEN.findall(text or "")))


@functools.lru_cache(maxsize=1)
def _tracked() -> frozenset:
    """Every tracked path, as git spells it (empty outside a git checkout)."""
    try:
        return frozenset(subprocess.run(["git", "ls-files"], cwd=data_paths.REPO_ROOT, capture_output=True,
                                        text=True, check=True).stdout.split())
    except Exception:  # noqa: BLE001
        return frozenset()


@functools.lru_cache(maxsize=1)
def _dossier_names() -> frozenset:
    return frozenset(p.name for p in data_paths.EXTRACTION_DOSSIERS_DIR.iterdir())


def resolve(token: str) -> Optional[Path]:
    """The file a token names: a tracked repository path, a dossier, or the one tracked file with that name.
    Names are matched exactly, never through the file system: macOS folds case and CI's Linux does not, which
    let a mis-cased name pass locally and fail in CI on 2026-10-29."""
    if "/" in token:
        rel = token[2:] if token.startswith("./") else token
        return data_paths.REPO_ROOT / rel if rel in _tracked() else None
    if token in _dossier_names():
        return data_paths.EXTRACTION_DOSSIERS_DIR / token
    named = [f for f in _tracked() if f.rsplit("/", 1)[-1] == token]
    return data_paths.REPO_ROOT / named[0] if len(named) == 1 else None


def dossier_stems(text: Optional[str]) -> List[str]:
    """The dossier stems an anchor cites (``kumazawa2003`` for ``kumazawa2003_extraction.md``)."""
    return [Path(t).name[: -len("_extraction.md")] for t in anchor_tokens(text) if t.endswith("_extraction.md")
            and resolve(t) is not None]
