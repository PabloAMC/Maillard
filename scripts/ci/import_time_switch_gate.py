#!/usr/bin/env python
"""
import_time_switch_gate -- no generator sets an engine switch when it is imported (2026-10-29).

The engine's switches are module globals (``sulfur.PH_COUPLING``, ``ph_state.B51``, ``ps.CYS_AMADORI_ROUTE`` ...).
A generator that assigns one at module level changes the engine for every importer: on 2026-10-29 the
route-calibration test imported B64, inherited the inline pH coupling, and B2.4's frozen-objective test then ran
the wrong configuration for over an hour. ``tests/conftest.py`` now resets the switches per test, but that only
defends; this gate prevents. A switch is set inside a function or a context manager (the route calibration's
``installed()``), never at import.

Checks every ``scripts/generators/*.py``: a top-level statement (outside any function or class) that assigns to
an UPPERCASE attribute of a name (``module.SWITCH = ...``), or calls the engine configuration's setters
(``config.update(...)`` / ``config.install(...)``, A3 phase 2) or ``setattr(module, "SWITCH", ...)``, fails, unless
the file is grandfathered below: the
frozen wave records written before this gate, whose behaviour their artifacts already record. A grandfathered
file that no longer does it must be removed from the list, so the list only shrinks.

    python scripts/ci/import_time_switch_gate.py
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GENERATORS = ROOT / "scripts" / "generators"

#: Frozen wave generators that set switches at import, recorded 2026-10-29. Do not add to this list.
GRANDFATHERED = {
    "generate_kinetic_core_b51_fit.py",
    "generate_kinetic_core_b58_fit.py",
    "generate_kinetic_core_b60_fit.py",
    "generate_kinetic_core_b60_verification.py",
    "generate_kinetic_core_b61_fit.py",
    "generate_kinetic_core_b62_fit.py",
    "generate_kinetic_core_b63_fit.py",
    "generate_kinetic_core_b64_fit.py",
}


#: the engine-configuration calls that set switches (src/kinetic_core/config.py, A3 phase 2)
_CONFIG_MODULE = "src.kinetic_core.config"
_CONFIG_SETTERS = ("install", "update")


def _config_setter_names(tree: ast.Module) -> tuple:
    """(aliases of the config module, bare names bound to its setters) in a file."""
    modules, functions = set(), set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module == "src.kinetic_core":
            modules |= {a.asname or a.name for a in node.names if a.name == "config"}
        elif isinstance(node, ast.ImportFrom) and node.module == _CONFIG_MODULE:
            functions |= {a.asname or a.name for a in node.names if a.name in _CONFIG_SETTERS}
        elif isinstance(node, ast.Import):
            modules |= {a.asname for a in node.names if a.name == _CONFIG_MODULE and a.asname}
    return modules, functions


def _top_level_switch_assignments(tree: ast.Module) -> list:
    hits = []
    config_modules, config_functions = _config_setter_names(tree)

    def setter_calls(node):
        """``config.install(...)`` / ``config.update(...)`` / ``setattr(module, "SWITCH", ...)`` in a statement
        (a top-level call sets switches at import as an assignment does)."""
        out = []
        for call in ast.walk(node):
            if not isinstance(call, ast.Call):
                continue
            f = call.func
            if (isinstance(f, ast.Attribute) and isinstance(f.value, ast.Name) and f.value.id in config_modules
                    and f.attr in _CONFIG_SETTERS) or (isinstance(f, ast.Name) and f.id in config_functions):
                out.append(f"config.{f.attr if isinstance(f, ast.Attribute) else f.id}(...)")
            elif (isinstance(f, ast.Name) and f.id == "setattr" and len(call.args) >= 2
                  and isinstance(call.args[1], ast.Constant) and str(call.args[1].value).isupper()):
                out.append(f"setattr(..., {call.args[1].value!r})")
        return out

    def visit(statements):
        for node in statements:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                continue
            targets = []
            if isinstance(node, ast.Assign):
                targets = node.targets
            elif isinstance(node, (ast.AugAssign, ast.AnnAssign)):
                targets = [node.target]
            for target in targets:
                for t in (target.elts if isinstance(target, ast.Tuple) else [target]):
                    if isinstance(t, ast.Attribute) and isinstance(t.value, ast.Name) and t.attr.isupper():
                        hits.append((node.lineno, f"{t.value.id}.{t.attr}"))
            if isinstance(node, (ast.Expr, ast.Assign, ast.AnnAssign)):
                hits.extend((node.lineno, what) for what in setter_calls(node.value))
            for block in ("body", "orelse", "finalbody", "handlers"):  # if / try / with at top level still run at import
                inner = getattr(node, block, None)
                if isinstance(inner, list) and not isinstance(node, (ast.FunctionDef, ast.ClassDef)):
                    visit([n for n in inner if isinstance(n, ast.stmt)] +
                          [s for h in inner if isinstance(h, ast.ExceptHandler) for s in h.body])

    visit(tree.body)
    return hits


def check() -> list:
    problems = []
    for path in sorted(GENERATORS.glob("*.py")):
        hits = _top_level_switch_assignments(ast.parse(path.read_text(encoding="utf-8"), filename=str(path)))
        if hits and path.name not in GRANDFATHERED:
            problems += [f"{path.name}:{line}: sets {name} at import" for line, name in hits]
        if not hits and path.name in GRANDFATHERED:
            problems.append(f"{path.name}: no longer sets a switch at import; remove it from GRANDFATHERED")
    return problems


def main() -> int:
    problems = check()
    if problems:
        print("import_time_switch_gate: FAIL")
        for p in problems:
            print("  - " + p)
        return 1
    print(f"import_time_switch_gate: PASS ({len(GRANDFATHERED)} frozen generators grandfathered; no new ones)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
