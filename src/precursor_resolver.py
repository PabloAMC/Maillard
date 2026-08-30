"""
src/precursor_resolver.py — Map user-friendly precursor names to Species objects.

Reads data/species/precursors.yml and provides fuzzy name matching so users
can type "ribose" or "cysteine" instead of full SMILES.
"""

import yaml
from pathlib import Path
from typing import List, Dict, Optional

from src.chem_utils import Species  # noqa: E402

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
    # 2026-08-28 (Wave X): `maillard_intermediates` added.  The tuple below is the
    # ONLY place the resolver learns which YAML blocks are precursor blocks, so a
    # category added to the file but not here is silently invisible -- that is how
    # `lipids` behaved before it was added, and it is why the list is spelled out
    # rather than derived from `data.keys()`: an accidental top-level key in the
    # YAML must not become a resolvable precursor.
    for category in ("amino_acids", "sugars", "exogenous_precursors", "lipids", "maillard_intermediates"):
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

            # 2026-08-28 (Wave X): explicit `aliases`.  The prefix/parenthesis
            # heuristics above cover "D-Ribose" -> "ribose" but cannot cover the
            # synonym pairs the sulfur literature actually uses -- furan-2-aldehyde /
            # furfural, 2-oxopropanal / methylglyoxal / pyruvaldehyde, norfuraneol /
            # NF, hydrogen sulfide / H2S.  Those are DIFFERENT WORDS, not decorations
            # on the same word, so no amount of string surgery reaches them and the
            # synonym has to be declared in the data file next to the SMILES it means.
            for alias in entry.get("aliases", []) or []:
                alias_key = str(alias).strip().lower()
                if alias_key:
                    keys.append(alias_key)

            for key in keys:
                # An alias must never silently steal a name another entry already
                # owns; that would make resolution depend on YAML ordering.
                if key in _LOOKUP and _LOOKUP[key]["name"] != name:
                    raise ValueError(
                        f"precursor key collision: {key!r} is claimed by both "
                        f"{_LOOKUP[key]['name']!r} and {name!r} in "
                        f"{_PRECURSORS_PATH.name}"
                    )
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
