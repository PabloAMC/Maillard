#!/usr/bin/env python
"""Hold-out exclusion guard.

The external-validation bundles under ``data/benchmarks/external_validation/`` are
the repo's only genuinely out-of-sample evidence. Their value is destroyed the
moment any of them leaks into calibration -- and the leak would be silent, showing
up only as an implausibly good hold-out score.

Three independent mechanisms keep them out. This guard asserts all three
*statically*, so a refactor cannot quietly remove one:

1. **Data-side label.** Every JSON under ``data/benchmarks/external_validation/``
   declares ``evidence_class: external_validation_only`` at top level.
2. **Code-side refusal.** ``src/matrix_experiment_intake.calibrate_from_intake``
   raises on an ``external_validation_only`` payload rather than calibrating it.
3. **Discovery-side isolation.** ``src.benchmark_validation.get_benchmark_files``
   globs ``*.json`` **non-recursively**, so the ``external_validation/``
   subdirectory is never picked up by default benchmark discovery. A change to
   ``rglob`` or ``**/*.json`` would silently pull the hold-out set into the
   calibration corpus; that is the single highest-consequence one-character
   regression in this repo.

The guard is deliberately static (AST + file reads, no imports of heavy deps), so
it runs in about a second and cannot be defeated by import-time monkeypatching.

Usage:
    python scripts/ci/holdout_guard.py
"""

from __future__ import annotations

import ast
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

HOLDOUT_DIR = ROOT / "data" / "benchmarks" / "external_validation"
INTAKE_MODULE = ROOT / "src" / "matrix_experiment_intake.py"
VALIDATION_MODULE = ROOT / "src" / "benchmark_validation.py"

EVIDENCE_CLASS = "external_validation_only"
GUARDED_FUNCTION = "calibrate_from_intake"
DISCOVERY_FUNCTION = "get_benchmark_files"

# Substrings that must appear in the guard's raise message, so that deleting the
# guard while keeping a same-named function still fails this check.
GUARD_MESSAGE_MARKERS = (
    "Refusing to calibrate",
    "hold-out",
)

# Glob methods / patterns that would defeat the non-recursive discovery invariant.
RECURSIVE_GLOB_METHODS = frozenset({"rglob"})
RECURSIVE_PATTERN_MARKER = "**"


def check_evidence_class(failures: list[str]) -> None:
    if not HOLDOUT_DIR.is_dir():
        failures.append(f"hold-out directory is missing: {HOLDOUT_DIR.relative_to(ROOT)}")
        return

    files = sorted(HOLDOUT_DIR.glob("*.json"))
    if not files:
        failures.append(
            f"no hold-out bundles found in {HOLDOUT_DIR.relative_to(ROOT)} - "
            "the external-validation lane has been emptied"
        )
        return

    for path in files:
        rel = path.relative_to(ROOT).as_posix()
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:
            failures.append(f"{rel}: not parseable JSON ({type(exc).__name__}: {exc})")
            continue

        declared = payload.get("evidence_class")
        if declared != EVIDENCE_CLASS:
            failures.append(
                f"{rel}: evidence_class is {declared!r}, must be {EVIDENCE_CLASS!r}. "
                "Without this label the bundle is eligible for calibration."
            )

    print(f"  [1/3] evidence_class: {len(files)} hold-out bundle(s) checked.")


def _find_function(tree: ast.Module, name: str) -> ast.FunctionDef | ast.AsyncFunctionDef | None:
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name:
            return node
    return None


def check_calibration_guard(failures: list[str]) -> None:
    rel = INTAKE_MODULE.relative_to(ROOT).as_posix()
    if not INTAKE_MODULE.is_file():
        failures.append(f"{rel} is missing")
        return

    source = INTAKE_MODULE.read_text(encoding="utf-8")
    tree = ast.parse(source)
    func = _find_function(tree, GUARDED_FUNCTION)
    if func is None:
        failures.append(f"{rel}: function {GUARDED_FUNCTION}() not found")
        return

    body = ast.get_source_segment(source, func) or ""

    if EVIDENCE_CLASS not in body:
        failures.append(
            f"{rel}:{func.lineno} {GUARDED_FUNCTION}() no longer mentions "
            f"{EVIDENCE_CLASS!r} - the hold-out refusal has been removed."
        )
        return

    # The guard must actually raise, not merely warn or log.
    raises_in_guard = any(
        isinstance(node, ast.Raise) for node in ast.walk(func)
    )
    if not raises_in_guard:
        failures.append(
            f"{rel}:{func.lineno} {GUARDED_FUNCTION}() references {EVIDENCE_CLASS!r} "
            "but contains no raise statement - the guard must refuse, not warn."
        )

    missing = [m for m in GUARD_MESSAGE_MARKERS if m not in body]
    if missing:
        failures.append(
            f"{rel}:{func.lineno} {GUARDED_FUNCTION}() guard message lost expected "
            f"marker(s) {missing!r}. If the wording changed deliberately, update "
            "GUARD_MESSAGE_MARKERS in scripts/ci/holdout_guard.py in the same commit."
        )

    if not failures:
        print(f"  [2/3] calibration guard: {GUARDED_FUNCTION}() raises on {EVIDENCE_CLASS}.")


def check_non_recursive_discovery(failures: list[str]) -> None:
    rel = VALIDATION_MODULE.relative_to(ROOT).as_posix()
    if not VALIDATION_MODULE.is_file():
        failures.append(f"{rel} is missing")
        return

    source = VALIDATION_MODULE.read_text(encoding="utf-8")
    tree = ast.parse(source)
    func = _find_function(tree, DISCOVERY_FUNCTION)
    if func is None:
        failures.append(f"{rel}: function {DISCOVERY_FUNCTION}() not found")
        return

    glob_calls: list[tuple[str, str]] = []  # (method, pattern-or-'<dynamic>')
    for node in ast.walk(func):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            method = node.func.attr
            if method in ("glob", "rglob", "iglob", "walk"):
                if node.args and isinstance(node.args[0], ast.Constant) and isinstance(node.args[0].value, str):
                    pattern = node.args[0].value
                else:
                    pattern = "<dynamic>"
                glob_calls.append((method, pattern))

    if not glob_calls:
        failures.append(
            f"{rel}:{func.lineno} {DISCOVERY_FUNCTION}() performs no recognisable glob. "
            "This guard can no longer verify that hold-out files stay out of default "
            "benchmark discovery; re-express the check or restore a literal glob."
        )
        return

    for method, pattern in glob_calls:
        if method in RECURSIVE_GLOB_METHODS:
            failures.append(
                f"{rel}:{func.lineno} {DISCOVERY_FUNCTION}() uses .{method}() - recursive "
                "discovery pulls data/benchmarks/external_validation/ into the default "
                "benchmark corpus and destroys the hold-out."
            )
        elif RECURSIVE_PATTERN_MARKER in pattern:
            failures.append(
                f"{rel}:{func.lineno} {DISCOVERY_FUNCTION}() globs {pattern!r} - the '**' "
                "wildcard recurses into data/benchmarks/external_validation/ and destroys "
                "the hold-out."
            )
        elif method == "walk" or pattern == "<dynamic>":
            failures.append(
                f"{rel}:{func.lineno} {DISCOVERY_FUNCTION}() uses a non-literal or "
                f"walk-based discovery ({method}, pattern={pattern!r}); this guard cannot "
                "statically prove the hold-out stays excluded."
            )

    if not any(f.startswith(f"{rel}:{func.lineno} {DISCOVERY_FUNCTION}") for f in failures):
        patterns = ", ".join(f"{m}({p!r})" for m, p in glob_calls)
        print(f"  [3/3] discovery isolation: {DISCOVERY_FUNCTION}() uses {patterns} (non-recursive).")


def main() -> int:
    print("holdout_guard: asserting hold-out exclusion invariants")
    failures: list[str] = []

    check_evidence_class(failures)
    check_calibration_guard(failures)
    check_non_recursive_discovery(failures)

    if failures:
        print(f"\nFAIL: {len(failures)} hold-out invariant violation(s):\n", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}\n", file=sys.stderr)
        print(
            "The external-validation set is the repo's only out-of-sample evidence. "
            "If any of these changes was intentional, the hold-out is no longer a "
            "hold-out and every external-validation number must be re-derived.",
            file=sys.stderr,
        )
        return 1

    print("\nholdout_guard: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
