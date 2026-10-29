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
import pkgutil
import sys
import tokenize
import types
from pathlib import Path
from typing import Iterator

ROOT = Path(__file__).resolve().parents[1]
GENERATORS = ROOT / "scripts" / "generators"

#: Wave modules that mutate B2.3's module state when imported.
SPLICING_WAVE_MODULES = (
    "generate_kinetic_core_b8_fit", "generate_kinetic_core_b9_fit", "generate_kinetic_core_b10_fit",
    "generate_kinetic_core_b11_fit", "generate_kinetic_core_b16_fit",
)


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


def _engine_switches() -> dict:
    """Every module-level UPPERCASE value of every ``src.kinetic_core`` module, as first imported.

    The engine's switches are module globals, and the wave generators set theirs AT IMPORT
    (``generate_kinetic_core_b64_fit`` sets ``sulfur.PH_COUPLING = "inline"``; B51 / B58 / B60-B63
    set ``ph_state.B51`` and route switches). A test that imports one, directly or through a gate,
    used to leave the switch set for every later test: 2026-10-29, the B2.4 unit test ran B2.3's
    frozen objective under the inline coupling, crawled for over an hour, and would have checked the
    wrong configuration had it finished. Taken when conftest imports this module, before collection,
    so no test's imports are in it.
    """
    import src.kinetic_core as kinetic_core

    switches = {}
    for info in pkgutil.iter_modules(kinetic_core.__path__):
        module = importlib.import_module(f"src.kinetic_core.{info.name}")
        for name, value in vars(module).items():
            if name.isupper() and not callable(value) and not isinstance(value, types.ModuleType):
                switches[(module, name)] = value
    return switches


_PRISTINE_SWITCHES = _engine_switches()


def reset_engine_switches() -> None:
    """Put every engine switch back to its shipped value. Run before each test (conftest), and by any
    fixture of wider scope, which pytest sets up BEFORE the per-test reset."""
    for (module, name), value in _PRISTINE_SWITCHES.items():
        if getattr(module, name, value) is not value:
            setattr(module, name, value)


def _ensure_generator_path() -> None:
    for entry in (str(ROOT), str(GENERATORS)):
        if entry not in sys.path:
            sys.path.insert(0, entry)


@contextlib.contextmanager
def wave_generator(name: str) -> Iterator[object]:
    """Yield a freshly imported wave generator; restore B2.3's row table on exit."""
    _ensure_generator_path()
    import generate_kinetic_core_b2_3_fit as b23  # noqa: E402

    # B10 also rebinds B2.3's parameter builder (two route barriers); snapshot it too.
    snapshot = (b23.ACTIVE_FIT_ROWS, b23.FIT_ROWS, dict(b23.SYSTEMS), b23.build_parameters)
    for module in SPLICING_WAVE_MODULES:
        sys.modules.pop(module, None)
    # 2026-10-29: the list above is not the whole set. B10, B11 and B17 also rebind B8's keys, and the
    # later waves chain-import them (B64 -> B60 -> ... -> B2.3). So every generator module already
    # loaded gets its namespace back on exit, and every one loaded inside is dropped.
    generators_before = {key: dict(vars(module)) for key, module in sys.modules.items()
                         if key.startswith("generate_")}
    try:
        yield importlib.import_module(name)
    finally:
        for key in [key for key in sys.modules if key.startswith("generate_")]:
            if key in generators_before:
                namespace = vars(sys.modules[key])
                namespace.clear()
                namespace.update(generators_before[key])
            else:
                sys.modules.pop(key)
        b23.ACTIVE_FIT_ROWS, b23.FIT_ROWS = snapshot[0], snapshot[1]
        b23.SYSTEMS.clear()
        b23.SYSTEMS.update(snapshot[2])
        b23.build_parameters = snapshot[3]
        for module in SPLICING_WAVE_MODULES:
            sys.modules.pop(module, None)
