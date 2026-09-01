from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from src import data_paths


def repo_root() -> Path:
    """Canonical project root (parent of the ``src/`` package directory)."""
    return data_paths.REPO_ROOT


def resolve_optional_path(file_path: Path | str | None, default_path: Path) -> Path:
    return Path(file_path) if file_path is not None else default_path


def load_json_mapping(path: Path | str) -> dict[str, Any]:
    with open(Path(path), "r", encoding="utf-8") as handle:
        return json.load(handle)

