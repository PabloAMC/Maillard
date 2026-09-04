#!/usr/bin/env python
"""Hold-out exclusion guard.

The external-validation bundles under ``data/benchmarks/external_validation/`` are
the repo's only genuinely out-of-sample evidence. Their value is destroyed the
moment any of them leaks into calibration -- and the leak would be silent, showing
up only as an implausibly good hold-out score.

Three independent mechanisms keep them out. This guard asserts all three
*statically*, so a refactor cannot quietly remove one:

1. **Data-side label.** Every JSON under ``data/benchmarks/external_validation/``
   -- *including its subdirectories* -- declares
   ``evidence_class: external_validation_only`` at top level.
2. **Fit-side isolation (2026-09-03, B5b).** No core fit generator
   (``scripts/generators/generate_kinetic_core_*_fit.py``) carries a string that
   names the hold-out directory.
3. **Discovery-side isolation.** ``src.kinetic_core.panel.panel_bundles`` globs
   ``*.json`` **non-recursively** per panel directory and tags each bundle; a
   change to ``rglob`` or ``**/*.json`` would silently pull quarantined bundles
   into the scored panel.

A fourth check was added on 2026-08-27 (Wave U), when the hold-out gained its
first bundles that actually execute the reaction network:

4. **Maillard-path bundles stay out of the fit corpus.** Wave S1 established that
   all four legacy hold-out bundles run ``matrix_only`` and never reach
   ``Recommender.predict_from_steps``; three consecutive waves of network changes
   moved zero hold-out points for that reason. The ``maillard_path/`` bundles are
   the repo's first out-of-sample test of the network itself, and they are the
   pre-registration target of the pending rate-calibration wave -- which makes
   them the single most attractive thing in the tree to quietly fit to. So every
   bundle flagged ``maillard_path_holdout: true`` must (a) live under the hold-out
   directory, (b) carry the evidence_class label, (c) declare
   ``metadata.execution_path: free_precursor`` -- a bundle that silently drifts to
   ``matrix_only`` stops testing the network and the panel would not notice -- and
   (d) have its ``benchmark_id`` absent from the fit/calibration corpus
   (``data/lit/`` registries and any ``results/validation/*refit*`` /
   ``*rederivation*`` record).

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

BENCHMARK_ROOT = ROOT / "data" / "benchmarks"
HOLDOUT_DIR = BENCHMARK_ROOT / "external_validation"
FIT_GENERATOR_DIR = ROOT / "scripts" / "generators"
PANEL_MODULE = ROOT / "src" / "kinetic_core" / "panel.py"
HOLDOUT_PATH_MARKERS = ("external_validation", "maillard_path")

EVIDENCE_CLASS = "external_validation_only"
MAILLARD_PATH_FLAG = "maillard_path_holdout"
REQUIRED_EXECUTION_PATH = "free_precursor"
DISCOVERY_FUNCTION = "panel_bundles"

# Glob methods / patterns that would defeat the non-recursive discovery invariant.
RECURSIVE_GLOB_METHODS = frozenset({"rglob"})
RECURSIVE_PATTERN_MARKER = "**"


def check_evidence_class(failures: list[str]) -> None:
    if not HOLDOUT_DIR.is_dir():
        failures.append(f"hold-out directory is missing: {HOLDOUT_DIR.relative_to(ROOT)}")
        return

    # 2026-08-27 (Wave U): rglob, not glob. The Maillard-path bundles live in the
    # maillard_path/ SUBDIRECTORY (so that they stay out of the matrix-only
    # external_validation_report median, which would otherwise average two
    # different execution paths into one meaningless figure). A non-recursive
    # glob here would have left them entirely unguarded.
    files = sorted(HOLDOUT_DIR.rglob("*.json"))
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

    print(f"  [1/4] evidence_class: {len(files)} hold-out bundle(s) checked.")


def _find_function(tree: ast.Module, name: str) -> ast.FunctionDef | ast.AsyncFunctionDef | None:
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name:
            return node
    return None


def check_fit_generators_never_open_the_holdout(failures: list[str]) -> None:
    """Check 2 (2026-09-03, retirement step B5b): no core FIT generator names the hold-out.

    The legacy intake guard (``matrix_experiment_intake.calibrate_from_intake``) is gone
    with the lane it guarded. The kinetic core is fitted by
    ``scripts/generators/generate_kinetic_core_*_fit.py``; none of them may carry a string
    constant that points into ``data/benchmarks/external_validation/``. (The one hold-out row
    the sulfur fit DID read, Hofmann 1998 xylose pH 5, is declared in
    ``src/kinetic_core/fit_targets.py`` and flagged on every scorecard row; it entered the
    fit as a literal table row, not by opening the bundle -- which is exactly why this check
    exists.)
    """
    generators = sorted(FIT_GENERATOR_DIR.glob("generate_kinetic_core_*_fit.py"))
    if not generators:
        failures.append(f"no core fit generators found under {FIT_GENERATOR_DIR.relative_to(ROOT)}")
        return
    before = len(failures)
    for path in generators:
        rel = path.relative_to(ROOT).as_posix()
        tree = ast.parse(path.read_text(encoding="utf-8"))
        # Docstrings and bare string statements are prose (they may well DISCUSS the
        # hold-out); only strings that feed code count.
        prose = {
            id(node.value)
            for node in ast.walk(tree)
            if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant)
            and isinstance(node.value.value, str)
        }
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, str) and id(node) not in prose:
                value = node.value.strip()
                # A PATH has no spaces (or is a glob / a .json name); a sentence that
                # discloses "... external_validation/ was never opened" is a declaration.
                looks_like_path = (" " not in value) or value.endswith(".json") or "*" in value
                if looks_like_path and any(marker in value for marker in HOLDOUT_PATH_MARKERS):
                    failures.append(
                        f"{rel}:{node.lineno} names the hold-out ({node.value!r}). A fit "
                        "generator must not open, glob or reference the hold-out directory."
                    )
    if len(failures) == before:
        print(f"  [2/4] fit isolation: {len(generators)} core fit generator(s) never name the hold-out.")


def check_non_recursive_discovery(failures: list[str]) -> None:
    """Check 3: the core panel's discovery must stay a literal, non-recursive glob.

    ``src/kinetic_core/panel.py::panel_bundles`` reads three directories on purpose (the
    trust loop, the maillard_path hold-out, the external matrix bundles) and TAGS each
    bundle with its panel. What it must never do is recurse: a ``**`` or ``rglob`` would
    pull quarantined and step-level-unreachable bundles -- or a future subdirectory -- into
    the scored panel silently.
    """
    rel = PANEL_MODULE.relative_to(ROOT).as_posix()
    if not PANEL_MODULE.is_file():
        failures.append(f"{rel} is missing")
        return
    source = PANEL_MODULE.read_text(encoding="utf-8")
    tree = ast.parse(source)
    func = _find_function(tree, DISCOVERY_FUNCTION)
    if func is None:
        failures.append(f"{rel}: function {DISCOVERY_FUNCTION}() not found")
        return
    glob_calls: list[tuple[str, str]] = []
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
            f"{rel}:{func.lineno} {DISCOVERY_FUNCTION}() performs no recognisable glob; "
            "re-express the check or restore a literal glob."
        )
        return
    for method, pattern in glob_calls:
        if method in RECURSIVE_GLOB_METHODS or RECURSIVE_PATTERN_MARKER in pattern or method == "walk" or pattern == "<dynamic>":
            failures.append(
                f"{rel}:{func.lineno} {DISCOVERY_FUNCTION}() uses {method}({pattern!r}): "
                "panel discovery must be a literal non-recursive glob."
            )
    if not any(f.startswith(f"{rel}:{func.lineno} {DISCOVERY_FUNCTION}") for f in failures):
        patterns = ", ".join(f"{m}({p!r})" for m, p in glob_calls)
        print(f"  [3/4] discovery isolation: {DISCOVERY_FUNCTION}() uses {patterns} (non-recursive).")


def check_maillard_path_bundles(failures: list[str]) -> None:
    """Check 4 (2026-08-27, Wave U): the network hold-out must stay a hold-out.

    See the module docstring for why these bundles are the highest-value fit
    target in the repository and therefore need their own invariant.
    """
    flagged: list[tuple[Path, dict]] = []

    # (a) Nothing outside the hold-out directory may claim the flag.
    for path in sorted(BENCHMARK_ROOT.rglob("*.json")):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            continue
        if not isinstance(payload, dict) or payload.get(MAILLARD_PATH_FLAG) is not True:
            continue
        rel = path.relative_to(ROOT).as_posix()
        if HOLDOUT_DIR not in path.parents:
            failures.append(
                f"{rel}: carries {MAILLARD_PATH_FLAG}=true but does not live under "
                f"{HOLDOUT_DIR.relative_to(ROOT).as_posix()}/. A Maillard-path hold-out "
                "point outside the hold-out directory is discoverable by default "
                "benchmark discovery and is no longer out of sample."
            )
            continue
        flagged.append((path, payload))

    if not flagged:
        print(f"  [4/4] maillard-path bundles: none flagged {MAILLARD_PATH_FLAG}=true.")
        return

    # (b)/(c) Label and execution path.
    for path, payload in flagged:
        rel = path.relative_to(ROOT).as_posix()
        if payload.get("evidence_class") != EVIDENCE_CLASS:
            failures.append(
                f"{rel}: {MAILLARD_PATH_FLAG}=true but evidence_class is "
                f"{payload.get('evidence_class')!r}, must be {EVIDENCE_CLASS!r}."
            )
        execution_path = (payload.get("metadata") or {}).get("execution_path")
        if execution_path != REQUIRED_EXECUTION_PATH:
            failures.append(
                f"{rel}: metadata.execution_path is {execution_path!r}, must be "
                f"{REQUIRED_EXECUTION_PATH!r}. The whole point of these bundles is that "
                "they run Recommender.predict_from_steps; on any other path they stop "
                "testing the reaction network and the panel cannot tell."
            )

    # (d) No fit/calibration record may name them.
    ids = {
        str(payload.get("benchmark_id"))
        for _, payload in flagged
        if payload.get("benchmark_id")
    }
    fit_corpus: list[Path] = []
    lit_dir = ROOT / "data" / "lit"
    if lit_dir.is_dir():
        fit_corpus.extend(sorted(lit_dir.rglob("*.json")))
    results_dir = ROOT / "results" / "validation"
    if results_dir.is_dir():
        fit_corpus.extend(
            path
            for path in sorted(results_dir.glob("*.json"))
            if "refit" in path.name or "rederivation" in path.name
        )
    for path in fit_corpus:
        try:
            text = path.read_text(encoding="utf-8")
        except Exception:
            continue
        for benchmark_id in sorted(ids):
            if benchmark_id in text:
                failures.append(
                    f"{path.relative_to(ROOT).as_posix()}: names the Maillard-path hold-out "
                    f"benchmark {benchmark_id!r}. These points are the pre-registration "
                    "target of the pending rate-calibration wave; a calibration record that "
                    "names one has destroyed the only out-of-sample evidence the reaction "
                    "network has."
                )

    if not any("maillard" in f.lower() or any(i in f for i in ids) for f in failures):
        print(
            f"  [4/4] maillard-path bundles: {len(flagged)} flagged, all "
            f"{REQUIRED_EXECUTION_PATH}, none named by any fit record "
            f"({len(fit_corpus)} records scanned)."
        )


def main() -> int:
    print("holdout_guard: asserting hold-out exclusion invariants")
    failures: list[str] = []

    check_evidence_class(failures)
    check_fit_generators_never_open_the_holdout(failures)
    check_non_recursive_discovery(failures)
    check_maillard_path_bundles(failures)

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
