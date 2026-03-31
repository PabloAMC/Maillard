from __future__ import annotations

from typing import Any, Dict, Mapping, Sequence


def normalize_aliases(payload: Mapping[str, Any], aliases: Mapping[str, Sequence[str]]) -> Dict[str, Any]:
    normalized = dict(payload)
    for canonical_key, alias_keys in aliases.items():
        if canonical_key in normalized:
            continue
        for alias_key in alias_keys:
            if alias_key in normalized:
                normalized[canonical_key] = normalized[alias_key]
                break
    return normalized
