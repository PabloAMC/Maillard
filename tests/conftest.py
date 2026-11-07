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

from tests.support import reset_engine_switches
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


#: set in the CI job that runs the JAX tests (maillard-jax env): there a skip is a FAILURE, because every JAX test
#: skips itself when jax / diffrax do not import (``pytest.importorskip``), and a job that silently skipped them all
#: would pass while testing nothing
REQUIRE_NO_SKIPS_ENV = "MAILLARD_REQUIRE_NO_SKIPS"


def pytest_sessionfinish(session, exitstatus):
    import os

    if os.environ.get(REQUIRE_NO_SKIPS_ENV, "") != "1":
        return
    reporter = session.config.pluginmanager.get_plugin("terminalreporter")
    skipped = reporter.stats.get("skipped", []) if reporter is not None else []
    if skipped:
        names = sorted({str(getattr(r, "nodeid", "") or getattr(r, "fspath", "")) for r in skipped})
        print(f"\n{REQUIRE_NO_SKIPS_ENV}=1: {len(skipped)} skipped, which fails this job: {names[:20]}")
        session.exitstatus = pytest.ExitCode.TESTS_FAILED


@pytest.fixture(autouse=True)
def _pristine_engine_switches():
    """Each test starts from the engine's shipped switches, whatever an earlier test imported
    (``tests.support.reset_engine_switches``; why: its docstring)."""
    reset_engine_switches()
    yield


@pytest.fixture
def wave_generator():
    """Factory: ``with wave_generator("generate_kinetic_core_b8_fit") as b8: ...``."""
    return _wave_generator
