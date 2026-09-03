"""Reading and writing generated artifacts (JSON + markdown twin).

2026-09-03: ``write_artifact`` is the one writer for ``results/`` artifacts. Before, the
scorecard, envelope, directional and experiment-value modules each carried their own with
different ``sort_keys``, trailing-newline and return conventions.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable, Mapping, Optional, Tuple

from src import data_access, data_paths


def repo_root() -> Path:
    """Canonical project root (parent of the ``src/`` package directory)."""
    return data_paths.REPO_ROOT


def resolve_optional_path(file_path: Path | str | None, default_path: Path) -> Path:
    return Path(file_path) if file_path is not None else default_path


def load_json_mapping(path: Path | str) -> dict[str, Any]:
    """A JSON object from disk (strict: a missing file raises)."""
    return dict(data_access.load_json(path))


def dump_json(payload: Mapping[str, Any], *, sort_keys: bool = False) -> str:
    """The one JSON serialisation for artifacts: indent 2, ``default=str``, trailing newline."""
    return json.dumps(payload, indent=2, sort_keys=sort_keys, default=str) + "\n"


def write_artifact(
    payload: Mapping[str, Any],
    json_path: Path | str,
    *,
    render: Optional[Callable[[Mapping[str, Any]], str]] = None,
    md_path: Path | str | None = None,
    sort_keys: bool = False,
) -> Tuple[Path, Optional[Path]]:
    """Write ``payload`` as JSON and, when ``render`` is given, its markdown twin.

    Returns ``(json_path, md_path)``; ``md_path`` is ``None`` when nothing was rendered.
    The markdown twin defaults to the JSON path with a ``.md`` suffix.
    """
    json_path = Path(json_path)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(dump_json(payload, sort_keys=sort_keys), encoding="utf-8")
    if render is None:
        return json_path, None
    md_target = Path(md_path) if md_path is not None else json_path.with_suffix(".md")
    md_target.write_text(render(payload), encoding="utf-8")
    return json_path, md_target


__all__ = ["dump_json", "load_json_mapping", "repo_root", "resolve_optional_path", "write_artifact"]
