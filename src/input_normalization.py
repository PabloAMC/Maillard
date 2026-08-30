from __future__ import annotations

from typing import Any, Dict, Iterable, Mapping, Optional, Sequence


# Every spelling of a process-condition key that occurs somewhere in this repo.
# Formulation grid YAML, campaign `shared_conditions`, `Formulation.to_dict()`
# and hand-written benchmark payloads do not agree on one casing, so any code
# that reads a condition off an untyped mapping must go through
# `resolve_condition_value` rather than guessing a single spelling.
#
# 2026-08-27 (audit remediation): added because the campaign leaderboard roll-up
# read `conditions.get("ph")` / `.get("temp")` straight off the per-run
# formulation dicts and printed pH 0.00 / Temp 0.0 for every row — those dicts
# come from the grid, which carries no conditions at all, while the values the
# pipeline actually used live in the campaign's shared conditions.
CONDITION_KEY_ALIASES: Mapping[str, tuple[str, ...]] = {
    "ph": ("ph", "pH", "PH"),
    "temp": ("temp", "TEMP", "temperature", "temperature_celsius", "temp_C", "temp_c"),
    "aw": ("aw", "AW", "water_activity"),
    "time_minutes": ("time_minutes", "time_min", "TIME_MINUTES"),
    "protein_type": ("protein_type", "protein"),
}


def _lookup(source: Any, key: str) -> Any:
    """Read `key` off a mapping or off an object exposing `.get` (e.g. Formulation)."""
    if source is None:
        return None
    if isinstance(source, Mapping):
        return source.get(key)
    getter = getattr(source, "get", None)
    if callable(getter):
        return getter(key)
    return getattr(source, key, None)


def resolve_condition_value(
    canonical_key: str,
    sources: Iterable[Any],
    default: Any = None,
) -> Any:
    """First non-None value for `canonical_key` across `sources`, honouring aliases.

    `sources` is searched in priority order (most specific first, e.g. the
    per-run formulation, then the campaign's shared conditions). Within one
    source every known spelling of the key is tried before moving on, so a
    caller never has to know which casing a particular file used.
    """
    aliases = CONDITION_KEY_ALIASES.get(canonical_key, (canonical_key,))
    for source in sources:
        for alias in aliases:
            value = _lookup(source, alias)
            if value is not None:
                return value
    return default


def resolve_condition_float(
    canonical_key: str,
    sources: Iterable[Any],
    default: Optional[float] = None,
) -> Optional[float]:
    """`resolve_condition_value` coerced to float, or `default` if absent/uncoercible.

    Returns `None` rather than a fabricated 0.0 when nothing can be resolved:
    a reported pH of 0.00 is a chemistry claim, and an unknown condition must
    read as unknown.
    """
    value = resolve_condition_value(canonical_key, sources, default=None)
    if value is None:
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


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
