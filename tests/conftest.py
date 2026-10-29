"""Pytest configuration.

Rewritten 2026-09-03 (test audit): the fifteen DFT / Cantera / recommendation-engine
fixtures that lived here had no remaining user (their modules were retired with the legacy
engine), and the auto-marker lists named five test files that no longer exist. What is
left is the one rule still in force: everything under tests/scientific/ carries the
``scientific_regression`` marker. Shared helpers live in tests/support.py.
"""
from __future__ import annotations

import importlib
import pkgutil
import types
from pathlib import Path

import pytest

from tests.support import wave_generator as _wave_generator

TESTS_DIR = Path(__file__).resolve().parent


def pytest_collection_modifyitems(config, items):
    """Everything under tests/scientific/ is a scientific regression test."""
    for item in items:
        path = Path(str(item.fspath)).resolve()
        try:
            relative = path.relative_to(TESTS_DIR)
        except ValueError:
            continue
        if relative.parts[:1] == ("scientific",):
            item.add_marker(pytest.mark.scientific_regression)


def _engine_switches() -> dict:
    """Every module-level UPPERCASE value of every ``src.kinetic_core`` module, as first imported.

    The engine's switches are module globals, and the wave generators set theirs AT IMPORT
    (``generate_kinetic_core_b64_fit`` sets ``sulfur.PH_COUPLING = "inline"``; B51 / B58 / B60-B63
    set ``ph_state.B51`` and route switches). A test that imports one, directly or through a gate,
    used to leave the switch set for every later test: 2026-10-29, the B2.4 unit test ran B2.3's
    frozen objective under the inline coupling, crawled for over an hour, and would have checked the
    wrong configuration had it finished. Taken before collection, so no test's imports are in it.
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


@pytest.fixture(autouse=True)
def _pristine_engine_switches():
    """Each test starts from the engine's shipped switches, whatever an earlier test imported."""
    for (module, name), value in _PRISTINE_SWITCHES.items():
        if getattr(module, name, value) is not value:
            setattr(module, name, value)
    yield


@pytest.fixture
def wave_generator():
    """Factory: ``with wave_generator("generate_kinetic_core_b8_fit") as b8: ...``."""
    return _wave_generator
