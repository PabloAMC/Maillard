"""
src/target_registry.py — Cached target YAML registry

Centralises all loading of species data YAML files (desirable targets,
off-flavour targets, toxic markers) that were previously loaded ad-hoc
inside the Recommender class on every instantiation.  Using @lru_cache
guarantees that each file is read from disk exactly once per process.
"""
from __future__ import annotations

import yaml
from functools import lru_cache
from pathlib import Path
from typing import Dict, Any

_ROOT = Path(__file__).resolve().parents[1]
_SPECIES_DIR = _ROOT / "data" / "species"


def _load_species_yaml(filename: str) -> Dict[str, Any]:
    """Load a species YAML file and return a name-keyed dict."""
    path = _SPECIES_DIR / filename
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    return {item["name"]: item for item in data.get("compounds", []) if item.get("name")}


@lru_cache(maxsize=None)
def get_toxic_markers() -> Dict[str, Any]:
    """Cached mapping of toxic-marker compound names → YAML entry."""
    return _load_species_yaml("toxic_markers.yml")


@lru_cache(maxsize=None)
def get_desirable_targets() -> Dict[str, Any]:
    """Cached mapping of desirable-target compound names → YAML entry."""
    return _load_species_yaml("desirable_targets.yml")


@lru_cache(maxsize=None)
def get_off_flavour_targets() -> Dict[str, Any]:
    """Cached mapping of off-flavour compound names → YAML entry."""
    return _load_species_yaml("off_flavour_targets.yml")


def get_target_lookup(target_tag: str = "meaty") -> Dict[str, Any]:
    """
    Return the primary target lookup dict keyed by canonical SMILES for the
    given tag.  Currently delegates to ``get_desirable_targets`` for all tags
    because the registry is tag-aware at the recommend layer.  Extend this
    function when separate registries per tag are introduced.
    """
    return get_desirable_targets()
