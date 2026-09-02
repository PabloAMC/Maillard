"""Maillard runtime package.

2026-09-03 (retirement step B5b): one engine. The public names below are the kinetic
core's front door and are re-exported LAZILY (PEP 562) so ``import src.data_paths`` does
not integrate anything.
"""
from __future__ import annotations

from importlib import import_module
from typing import Any

_LAZY_EXPORTS = {
    "predict": "src.kinetic_core.engine",
    "FormulationSpec": "src.kinetic_core.engine",
    "ProcessSpec": "src.kinetic_core.engine",
    "ThermalProgram": "src.kinetic_core.engine",
    "EnvelopeDeclaration": "src.kinetic_core.engine",
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
