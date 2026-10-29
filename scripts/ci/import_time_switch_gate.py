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
an UPPERCASE attribute of a name (``module.SWITCH = ...``) fails, unless the file is grandfathered below: the
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


def _top_level_switch_assignments(tree: ast.Module) -> list:
    hits = []

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
