"""Loading curated data: one place, strict, cached.

Companion to :mod:`src.data_paths` (which says WHERE files are; this module says HOW
they are read). Rules:

* **Missing or malformed files raise.** Until 2026-09-01 most loaders returned ``{}``
  when a registry was absent and the model ran on quietly degraded inputs (a CWD-
  relative Henry's-law path silently turned every partition coefficient into the
  hard-coded 0.01). A curated file that cannot be read is a broken checkout, not a
  configuration to tolerate. :class:`DataFileError` names the file.
* **Read once per process.** ``load_json`` / ``load_yaml`` cache by resolved path and
  mtime, so the eleven modules that each opened ``computational_priors.json`` share one
  parse. Tests that rewrite a file in place can call :func:`clear_cache`.
* **Return plain ``dict`` / ``list``.** Validation against a schema is a later phase and
  will live here too; nothing else should grow its own ``json.load`` of a curated file.

Callers that genuinely want "absent is fine" -- the tombstone patch file, an optional
override -- say so explicitly with ``load_json(path, missing_ok=True)``, which returns
``None`` (never ``{}``, so absence cannot be mistaken for an empty registry).
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import yaml

from src import data_paths


class DataFileError(RuntimeError):
    """A curated data file is missing or cannot be parsed."""


class DataFileMissing(DataFileError, FileNotFoundError):
    """The file is not there. Also a ``FileNotFoundError`` so older handlers still catch it."""


class DataFileMalformed(DataFileError, ValueError):
    """The file exists but does not parse. Also a ``ValueError``."""


_CACHE: Dict[Tuple[str, str], Tuple[float, Any]] = {}


def _resolve(path: Path | str) -> Path:
    candidate = Path(path)
    if not candidate.is_absolute():
        candidate = data_paths.REPO_ROOT / candidate
    return candidate.resolve()


def _read(path: Path | str, kind: str, *, missing_ok: bool) -> Optional[Any]:
    resolved = _resolve(path)
    key = (kind, str(resolved))
    try:
        mtime = os.stat(resolved).st_mtime
    except FileNotFoundError:
        if missing_ok:
            return None
        raise DataFileMissing(
            f"curated data file is missing: {data_paths.rel(resolved)} "
            f"(resolved to {resolved}). Curated inputs are read-only and must be "
            f"present in the checkout; see src/data_paths.py."
        ) from None
    cached = _CACHE.get(key)
    if cached is not None and cached[0] == mtime:
        return cached[1]
    try:
        text = resolved.read_text(encoding="utf-8")
        payload = json.loads(text) if kind == "json" else yaml.safe_load(text)
    except (OSError, ValueError, yaml.YAMLError) as exc:
        raise DataFileMalformed(
            f"could not parse {data_paths.rel(resolved)}: {type(exc).__name__}: {exc}"
        ) from exc
    _CACHE[key] = (mtime, payload)
    return payload


def load_json(path: Path | str, *, missing_ok: bool = False) -> Any:
    """Parse a JSON data file. Raises :class:`DataFileError` unless ``missing_ok``."""
    return _read(path, "json", missing_ok=missing_ok)


def load_yaml(path: Path | str, *, missing_ok: bool = False) -> Any:
    """Parse a YAML data file. Raises :class:`DataFileError` unless ``missing_ok``."""
    return _read(path, "yaml", missing_ok=missing_ok)


def load_mapping(path: Path | str, *, missing_ok: bool = False) -> Optional[Dict[str, Any]]:
    """Like ``load_json``/``load_yaml`` (by suffix) but insists on a top-level mapping."""
    suffix = Path(path).suffix.lower()
    payload = load_yaml(path, missing_ok=missing_ok) if suffix in {".yml", ".yaml"} else load_json(
        path, missing_ok=missing_ok
    )
    if payload is None:
        return None
    if not isinstance(payload, dict):
        raise DataFileError(
            f"{data_paths.rel(_resolve(path))}: expected a top-level mapping, got {type(payload).__name__}"
        )
    return payload


def clear_cache() -> None:
    """Forget every cached parse (for tests that rewrite a data file in place)."""
    _CACHE.clear()
