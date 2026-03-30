from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def repo_root() -> Path:
    """Canonical project root (parent of the ``src/`` package directory)."""
    return Path(__file__).resolve().parents[1]


def resolve_optional_path(file_path: Path | str | None, default_path: Path) -> Path:
    return Path(file_path) if file_path is not None else default_path


def load_json_mapping(path: Path | str) -> dict[str, Any]:
    with open(Path(path), "r", encoding="utf-8") as handle:
        return json.load(handle)


def load_optional_json_mapping(path: Path | str) -> dict[str, Any]:
    candidate = Path(path)
    if not candidate.exists():
        return {}
    return load_json_mapping(candidate)

