"""Maillard runtime package.

The public names below are re-exported LAZILY (PEP 562). Until 2026-09-01 this file
imported ``src.pipeline`` eagerly, so ``import src.anything`` -- including
``src.data_paths`` -- pulled in the whole pipeline and read eleven data files off disk
before the caller's first line ran. Import the module you need; the convenience names
still resolve on first access.
"""
from __future__ import annotations

from importlib import import_module
from typing import Any

_LAZY_EXPORTS = {
    "MaillardPipeline": "src.pipeline",
    "FormulationResult": "src.pipeline",
    "ReactionConditions": "src.smirks_engine",
    "SmirksEngine": "src.smirks_engine",
    "FormulationOptimizer": "src.bayesian_optimizer",
    "Formulation": "src.formulation",
}

__all__ = list(_LAZY_EXPORTS)


def __getattr__(name: str) -> Any:
    module_name = _LAZY_EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module 'src' has no attribute {name!r}")
    value = getattr(import_module(module_name), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
