"""Which benchmarks were used to FIT the constants that are then scored against them.

2026-08-27 (Wave I).

Motivation
----------
The cold-start red team found that Wave H had re-derived a constant against two
benchmarks and then reported those same two benchmarks' agreement as literature
coverage. That is circular, and nothing in the repo could see it: the fit records and
the scored artifact had no shared vocabulary.

This module gives them one. It reads the fit records under ``results/validation/`` and
answers two questions about any benchmark:

* **Was it a fit target?** -- and of which record.
* **At what leverage?** -- how many free parameters that record fitted against how many
  rows.

Leverage is the distinction that keeps this honest in both directions:

``per_row_recovery`` (parameters per fitted row >= :data:`FIT_LEVERAGE_THRESHOLD`)
    The fit has enough freedom to reproduce its target row by row. Agreement carries no
    information about the model, so such rows are removed from the literature-coverage
    numerator AND denominator and reported separately, with their outcomes. Example: the
    matrix observability factors, one free factor per compound per lane, which is why
    those benchmarks score Pearson exactly 1.000.

``global_low_leverage`` (below the threshold)
    One or two global constants fitted across many rows. Individual-row agreement is
    still informative -- two parameters cannot make twenty-four rows agree -- so these
    rows STAY in the coverage denominator. They are annotated ``fit_target_of`` so a
    reader knows they are not strictly out-of-sample, but they are not excluded, because
    excluding them would quietly delete genuine FAILURES from the count. Example: the
    projection budget's ``reference_conversion_time_min``, a single global scale fitted
    across the whole literature panel.

The rule is deliberately symmetric about which direction it flatters. Excluding a fitted
row removes its misses as well as its hits; a threshold that excluded everything would
make a failing panel look empty rather than failing.

``scripts/ci/fit_target_gate.py`` re-implements this logic with the standard library
only (no ``src`` import), so the gate keeps running in the dependency-free CI tier. The
two must agree; the gate's own checks fail loudly if they drift.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterator, List, Tuple

ROOT = Path(__file__).resolve().parents[1]
VALIDATION_DIR = ROOT / "results" / "validation"

#: Globs whose JSON records are read for fit-target declarations.
FIT_RECORD_GLOBS: Tuple[str, ...] = (
    "*refit*.json",
    "*rederivation*.json",
    "*calibration_refit*.json",
)

#: JSON keys under which a record names what it was fitted against.
FIT_TARGET_KEYS: Tuple[str, ...] = (
    "fit_target",
    "fit_targets",
    "fit_target_files",
    "fit_target_ids",
)

#: Keys marking a record as withdrawn. A retracted record's targets are not live.
RETRACTION_KEYS: Tuple[str, ...] = ("RETRACTED", "retracted")

#: Free parameters per fitted row at or above which a fit can reproduce its own target.
FIT_LEVERAGE_THRESHOLD = 0.5


def _iter_strings(node: Any) -> Iterator[str]:
    if isinstance(node, str):
        yield node
    elif isinstance(node, list):
        for item in node:
            yield from _iter_strings(item)
    elif isinstance(node, dict):
        for item in node.values():
            yield from _iter_strings(item)


def benchmark_id_from_reference(token: str) -> str:
    """Normalise a fit-target reference (id, filename, or path) to a benchmark id."""
    text = str(token).strip()
    if not text:
        return ""
    head = text.split()[0].strip(",;")
    stem = Path(head).name
    if stem.endswith(".json"):
        stem = stem[: -len(".json")]
    return stem


@lru_cache(maxsize=1)
def load_fit_records(validation_dir: str = "") -> Tuple[Dict[str, Any], ...]:
    """Every live (non-retracted) fit record, as ``{name, path, targets, leverage}``."""
    directory = Path(validation_dir) if validation_dir else VALIDATION_DIR
    records: List[Dict[str, Any]] = []
    seen: set[Path] = set()
    for pattern in FIT_RECORD_GLOBS:
        for path in sorted(directory.glob(pattern)):
            if path in seen:
                continue
            seen.add(path)
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            if not isinstance(payload, dict):
                continue
            if any(key in payload for key in RETRACTION_KEYS):
                continue
            targets: List[str] = []
            for key in FIT_TARGET_KEYS:
                if key in payload:
                    targets.extend(
                        benchmark_id_from_reference(raw) for raw in _iter_strings(payload[key])
                    )
            targets = sorted({t for t in targets if t})
            if not targets:
                continue
            records.append(
                {
                    "name": path.name,
                    "path": str(path),
                    "targets": tuple(targets),
                    "leverage": dict(payload.get("fit_leverage", {}) or {}),
                }
            )
    return tuple(records)


def fit_leverage_class(record: Dict[str, Any]) -> str:
    """``"per_row_recovery"`` | ``"global_low_leverage"`` | ``"undeclared"``."""
    leverage = record.get("leverage") or {}
    declared = str(leverage.get("class", "")).strip()
    if declared in {"per_row_recovery", "global_low_leverage"}:
        return declared
    try:
        free = float(leverage["free_parameters"])
        rows = float(leverage["fitted_rows"])
    except (KeyError, TypeError, ValueError):
        return "undeclared"
    if rows <= 0:
        return "undeclared"
    return "per_row_recovery" if (free / rows) >= FIT_LEVERAGE_THRESHOLD else "global_low_leverage"


@lru_cache(maxsize=1)
def fit_target_map() -> Dict[str, Tuple[Dict[str, str], ...]]:
    """benchmark id -> the live fit records naming it, with each one's leverage class."""
    mapping: Dict[str, List[Dict[str, str]]] = {}
    for record in load_fit_records():
        klass = fit_leverage_class(record)
        for benchmark_id in record["targets"]:
            mapping.setdefault(benchmark_id, []).append(
                {"record": record["name"], "leverage_class": klass}
            )
    return {key: tuple(value) for key, value in sorted(mapping.items())}


def is_per_row_fit_target(benchmark_id: str) -> bool:
    """True when some live fit record could reproduce this benchmark row by row.

    Such a benchmark's agreement is algebraic recovery, so its rows are excluded from
    the literature-coverage numerator and denominator and reported separately.
    """
    for entry in fit_target_map().get(str(benchmark_id), ()):
        if entry["leverage_class"] in {"per_row_recovery", "undeclared"}:
            return True
    return False


def fit_target_records_for(benchmark_id: str) -> Tuple[str, ...]:
    """Names of the live fit records that name this benchmark as a target."""
    return tuple(entry["record"] for entry in fit_target_map().get(str(benchmark_id), ()))
