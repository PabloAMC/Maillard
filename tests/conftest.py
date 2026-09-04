"""Pytest configuration.

Rewritten 2026-09-03 (test audit): the fifteen DFT / Cantera / recommendation-engine
fixtures that lived here had no remaining user (their modules were retired with the legacy
engine), and the auto-marker lists named five test files that no longer exist. What is
left is the one rule still in force: everything under tests/scientific/ carries the
``scientific_regression`` marker. Shared helpers live in tests/support.py.
"""
from __future__ import annotations

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


@pytest.fixture
def wave_generator():
    """Factory: ``with wave_generator("generate_kinetic_core_b8_fit") as b8: ...``."""
    return _wave_generator
