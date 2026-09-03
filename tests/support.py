"""Helpers shared by the test suite (2026-09-03, test audit).

Two things live here because several test files had grown their own copy:

* ``executable_code(path)`` -- a file's code with comments, docstrings and string
  literals removed, for firewall tests that must not be fooled by prose that *mentions*
  a forbidden token (a docstring promising "no optimiser" contains the word).
* ``wave_generator(name)`` -- import a frozen wave generator that splices its rows into
  ``generate_kinetic_core_b2_3_fit`` at import time (B8 and B9 do), and undo the splice
  afterwards. Without this, whichever test imported B8 first silently changed the row
  table every later B2.3 / B2.4 test scored against, and the suite passed only in
  alphabetical order.
"""
from __future__ import annotations

import contextlib
import importlib
import io
import sys
import tokenize
from pathlib import Path
from typing import Iterator

ROOT = Path(__file__).resolve().parents[1]
GENERATORS = ROOT / "scripts" / "generators"

#: Wave modules that mutate B2.3's module state when imported.
SPLICING_WAVE_MODULES = ("generate_kinetic_core_b8_fit", "generate_kinetic_core_b9_fit")


def strip_prose(text: str) -> str:
    """Strip comments, docstrings and string literals from Python source; return code tokens."""
    out = []
    for token in tokenize.tokenize(io.BytesIO(text.encode("utf-8")).readline):
        if token.type in (tokenize.COMMENT, tokenize.STRING):
            continue
        out.append(token.string)
    return " ".join(out)


def executable_code(path: Path | str) -> str:
    """``strip_prose`` of a file."""
    return strip_prose(Path(path).read_text(encoding="utf-8"))


def _ensure_generator_path() -> None:
    for entry in (str(ROOT), str(GENERATORS)):
        if entry not in sys.path:
            sys.path.insert(0, entry)


@contextlib.contextmanager
def wave_generator(name: str) -> Iterator[object]:
    """Yield a freshly imported wave generator; restore B2.3's row table on exit."""
    _ensure_generator_path()
    import generate_kinetic_core_b2_3_fit as b23  # noqa: E402

    snapshot = (b23.ACTIVE_FIT_ROWS, b23.FIT_ROWS, dict(b23.SYSTEMS))
    for module in SPLICING_WAVE_MODULES:
        sys.modules.pop(module, None)
    try:
        yield importlib.import_module(name)
    finally:
        b23.ACTIVE_FIT_ROWS, b23.FIT_ROWS = snapshot[0], snapshot[1]
        b23.SYSTEMS.clear()
        b23.SYSTEMS.update(snapshot[2])
        for module in SPLICING_WAVE_MODULES:
            sys.modules.pop(module, None)
