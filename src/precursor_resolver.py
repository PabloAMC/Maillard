"""
src/precursor_resolver.py — Map user-friendly precursor names to Species objects.

Reads data/species/precursors.yml and provides fuzzy name matching so users
can type "ribose" or "cysteine" instead of full SMILES.
"""

import yaml
from pathlib import Path
from typing import List, Dict, Optional

from src.pathway_extractor import Species  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
_PRECURSORS_PATH = ROOT / "data" / "species" / "precursors.yml"

# Lazily loaded cache
_LOOKUP: Optional[Dict[str, dict]] = None


def _load_precursors() -> Dict[str, dict]:
    """Load precursors.yml and build a case-insensitive lookup dict."""
    global _LOOKUP
    if _LOOKUP is not None:
        return _LOOKUP

    with open(_PRECURSORS_PATH) as f:
        data = yaml.safe_load(f)

    _LOOKUP = {}
    for category in ("amino_acids", "sugars", "exogenous_precursors", "lipids"):
        for entry in data.get(category, []):
            name = entry["name"]
            # Index by multiple keys for fuzzy matching:
            # "L-Cysteine" -> keys: "l-cysteine", "cysteine"
            # "D-Ribose"   -> keys: "d-ribose", "ribose"
            # "Glycine"    -> keys: "glycine"
            keys = [name.lower()]
            # Strip stereochemical prefix (D-, L-)
            for prefix in ("l-", "d-"):
                if name.lower().startswith(prefix):
                    keys.append(name.lower()[len(prefix):])
            # Also allow the name without parenthetical notes
            if "(" in name:
                base = name.split("(")[0].strip()
                keys.append(base.lower())
                # Handle prefixes on the base name too
                for prefix in ("l-", "d-"):
                    if base.lower().startswith(prefix):
                        keys.append(base.lower()[len(prefix):])

            for key in keys:
                _LOOKUP[key] = entry

    return _LOOKUP


def resolve(name: str) -> Species:
    """
    Resolve a precursor name to a Species object.

    Accepts case-insensitive names like "ribose", "D-ribose", "cysteine",
    "L-Cysteine", etc.

    Raises ValueError if the name is not found.
    """
    lookup = _load_precursors()
    key = name.strip().lower()

    if key not in lookup:
        available = sorted(set(e["name"] for e in lookup.values()))
        raise ValueError(
            f"Unknown precursor '{name}'. Available: {', '.join(available)}"
        )

    entry = lookup[key]
    return Species(label=entry["name"], smiles=entry["smiles"])


def resolve_many(names: List[str]) -> List[Species]:
    """Resolve a list of precursor names to Species objects."""
    return [resolve(n) for n in names]


def list_available() -> List[str]:
    """Return all available precursor names."""
    lookup = _load_precursors()
    return sorted(set(e["name"] for e in lookup.values()))
