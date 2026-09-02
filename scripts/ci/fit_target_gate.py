#!/usr/bin/env python
"""Fit-target circularity + constant-precision gate.

Why this gate exists
--------------------
Both of the failure modes it checks actually happened in this repository, during the
audit that was supposed to be catching them. This gate is the machinery that makes each
one a build failure rather than a thing a reviewer has to notice.

**Failure mode 1 -- fit, then score, then report as validation.**
Wave H (2026-08-27) re-derived the Methional hydrolysate observability constant against
two benchmarks, and those same two benchmarks were then scored in the coverage headline.
They produced the *only two hits* in "2 of 11 literature rows inside the CI". A constant
solved from a measurement will reproduce that measurement; reporting the reproduction as
agreement is circular. (Both files turned out to be fabricated as well, which is a
separate finding -- but the circularity would have been wrong even if the data had been
real.)

**Failure mode 2 -- fit residue shipped as a constant.**
The matrix calibration registry carried `observable_factor=0.22877612093571738`, an
8-to-17-significant-figure number sitting next to a benchmark row reporting Pearson
R = 1.000. No measurement in this domain supports 17 significant figures; a constant
written that way is the output of a solve, and its precision is the tell. Constants may
still be fitted -- they just have to say so and be rounded to a precision they can
defend, with the prior value retained for audit.

Checks (all offline, all blocking)
----------------------------------
1. **No undisclosed fit-then-score.** Any benchmark named as a fit target by a record
   under ``results/validation/*refit*`` / ``*rederivation*`` that ALSO appears as a
   scored row in ``prediction_uncertainty.json`` must carry ``fitted_row: true`` on that
   scored row.
2. **Fitted rows are excluded from the coverage claim.** The scored artifact must carry
   a ``summary.honest_literature_coverage`` block, its denominator must exclude every
   fitted row, and the fitted rows must be reported separately.
3. (Retired 2026-09-03 with the legacy matrix calibration registry.) The core's own
   fit-row ledger is ``src/kinetic_core/fit_targets.py``, checked by its unit tests; it is
   scored. (Catches the reverse drift: someone marks the registry honestly but the
   reporting layer keeps counting the row.)
4. **Constant precision.** Shipped float constants in the calibration registries must not
   carry more than ``MAX_SIGNIFICANT_FIGURES`` significant figures. ``previous_value=``
   is exempt by design -- it exists precisely to retain the unrounded prior value.

This gate is deliberately dependency-free (stdlib only, ``ast`` for the Python side) so
it runs in the fast CI tier and cannot be taken down by a dependency solve.

Usage
-----
    python scripts/ci/fit_target_gate.py
    python scripts/ci/fit_target_gate.py --verbose
"""

from __future__ import annotations

import argparse
import ast
import json
import re
import sys
from decimal import Decimal
from pathlib import Path
from typing import Any, Dict, Iterator, List, Sequence, Set, Tuple

ROOT = Path(__file__).resolve().parents[2]

VALIDATION_DIR = ROOT / "results" / "validation"
# 2026-09-03 (retirement step B5b): the scored artifact is the KINETIC CORE's envelope.
PREDICTION_UNCERTAINTY = VALIDATION_DIR / "core_prediction_uncertainty.json"

# Records that declare a fit. Any file matching one of these globs is read for fit
# targets; a new refit generator is picked up automatically because it will write a
# record whose name matches.
FIT_RECORD_GLOBS = ("*refit*.json", "*rederivation*.json", "*calibration_refit*.json")

# JSON keys under which a fit record names what it fitted against.
FIT_TARGET_KEYS = ("fit_target", "fit_targets", "fit_target_files", "fit_target_ids")

# Keys marking a record as withdrawn. A retracted record's targets are not treated as
# live fit targets -- but the retraction has to be explicit, which is the point.
RETRACTION_KEYS = ("RETRACTED", "retracted")

MAX_SIGNIFICANT_FIGURES = 6

# Free parameters per fitted row at or above which a fit can reproduce its own target.
# Kept in sync with src.fit_target_index.FIT_LEVERAGE_THRESHOLD (that module cannot be
# imported here: this gate runs in the dependency-free CI tier).
FIT_LEVERAGE_THRESHOLD = 0.5

# Python sources holding shipped calibration constants, and the top-level assignment
# targets within them whose float literals are "shipped constants". Restricting to named
# containers keeps the check meaningful: a stray precision-carrying float in an unrelated
# helper is not a shipped calibration constant.
CONSTANT_SOURCES: Dict[str, Tuple[str, ...]] = {
    # 2026-09-03 (B5b): the legacy containers (matrix_calibration_registry, recommend) are
    # deleted. The kinetic core's fitted constants live in results/validation/
    # kinetic_core_b*_fit_report.json, not in code literals; the declared-assumption
    # tables in src/kinetic_core/parameters_*.py are the next candidates (backlog).
}

# Keyword arguments whose value is deliberately allowed to keep full precision, because
# their entire job is to preserve an unrounded historical value for audit.
PRECISION_EXEMPT_KEYWORDS = frozenset({"previous_value"})


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def _load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:  # pragma: no cover - CI diagnostics
        raise SystemExit(f"fit_target_gate: cannot read {_rel(path)}: {exc}")


def _benchmark_id_from(token: str) -> str:
    """Normalise a fit-target reference to a benchmark id.

    Records name targets either as a benchmark id ("cys_ribose_140C_Hofmann1998") or as a
    file name ("cys_ribose_140C_Hofmann1998.json"), sometimes with a path or a trailing
    parenthetical comment. Take the stem of the first whitespace-delimited token.

    WAVE S3 (2026-08-27): `.yml`/`.yaml` are stripped as well as `.json`. Until this
    wave every fit in the repository was against a benchmark JSON, so stripping `.json`
    alone was sufficient. The trunk rate calibration is fitted against the harvested
    concentration-time series in `data/lit/timeseries/*.yml`, and a declaration of
    `".../martins2005_glucose_glycine_80_100_120C_pH68.yml"` previously normalised to a
    pseudo-id that still carried its suffix -- it could never collide with a scored
    benchmark id, so the declaration passed the gate while buying no circularity
    protection at all. Stripping the suffix makes a YAML corpus a first-class
    declaration on the same footing as a JSON one. The principle the gate enforces is
    unchanged: anything used in a fit must be machine-declared, and nothing declared may
    appear in scored evidence.

    `src/fit_target_index.py::benchmark_id_from_reference` re-implements this and MUST
    be kept in lockstep.
    """
    text = str(token).strip()
    if not text:
        return ""
    head = text.split()[0].strip(",;")
    stem = Path(head).name
    for suffix in (".json", ".yml", ".yaml"):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return stem


def _iter_strings(node: Any) -> Iterator[str]:
    if isinstance(node, str):
        yield node
    elif isinstance(node, list):
        for item in node:
            yield from _iter_strings(item)
    elif isinstance(node, dict):
        for item in node.values():
            yield from _iter_strings(item)


def _significant_figures(literal: str) -> int:
    """Significant figures in a numeric literal, read from the SOURCE TEXT.

    Working from the text rather than the float matters: 0.228776 and
    0.22877612093571738 are different literals, and the second is the one that tells you
    a solver wrote it. Trailing zeros after a decimal point count (they are written
    deliberately); leading zeros never do.
    """
    text = literal.strip().replace("_", "")
    if not text:
        return 0
    mantissa = re.split(r"[eE]", text)[0]
    mantissa = mantissa.lstrip("+-")
    digits = mantissa.replace(".", "")
    digits = digits.lstrip("0")
    if not digits:
        return 0
    if "." not in mantissa:
        digits = digits.rstrip("0") or "0"
    return len(digits)


# ---------------------------------------------------------------------------
# check 1-3: fit-then-score circularity
# ---------------------------------------------------------------------------
def _leverage_class(record: Dict[str, Any]) -> str:
    """Mirror of src.fit_target_index.fit_leverage_class, stdlib only."""
    leverage = record.get("fit_leverage") or {}
    declared = str(leverage.get("class", "")).strip()
    if declared in {"per_row_recovery", "global_low_leverage"}:
        return declared
    try:
        free = float(leverage["free_parameters"])
        rows = float(leverage["fitted_rows"])
    except (KeyError, TypeError, ValueError):
        return "undeclared"
    if rows <= 0:
        return "undeclared"
    return "per_row_recovery" if (free / rows) >= FIT_LEVERAGE_THRESHOLD else "global_low_leverage"


def _collect_declared_fit_targets(verbose: bool) -> Tuple[Dict[str, List[str]], List[str]]:
    """(benchmark id -> records naming it as a target, declaration failures).

    Records are classified by LEVERAGE. Only a ``per_row_recovery`` fit -- one with
    enough free parameters to reproduce its target row by row -- forces its targets to be
    flagged ``fitted_row``. A low-leverage global fit (a couple of constants across many
    rows) is recorded but does not force exclusion, because excluding those rows would
    delete genuine FAILURES from the count rather than expose them. Every live record must
    nevertheless DECLARE its leverage: "undeclared" is a failure, not a default.
    """
    targets: Dict[str, List[str]] = {}
    failures: List[str] = []
    seen: Set[Path] = set()
    for pattern in FIT_RECORD_GLOBS:
        for path in sorted(VALIDATION_DIR.glob(pattern)):
            if path in seen:
                continue
            seen.add(path)
            record = _load_json(path)
            if not isinstance(record, dict):
                continue
            if any(key in record for key in RETRACTION_KEYS):
                if verbose:
                    print(f"    (skipping retracted record {_rel(path)})")
                continue
            declared: List[str] = []
            for key in FIT_TARGET_KEYS:
                if key not in record:
                    continue
                for raw in _iter_strings(record[key]):
                    benchmark_id = _benchmark_id_from(raw)
                    if benchmark_id:
                        declared.append(benchmark_id)
            if not declared:
                continue
            klass = _leverage_class(record)
            if verbose:
                print(f"    {_rel(path)}: leverage={klass}, {len(set(declared))} target(s)")
            if klass == "undeclared":
                failures.append(
                    f"UNDECLARED FIT LEVERAGE: {_rel(path)} names fit targets but has no "
                    "`fit_leverage` block. State `free_parameters` and `fitted_rows` (or an "
                    "explicit `class`). Without it a reader cannot tell one global scale "
                    "fitted across two dozen rows from one free knob per row -- and only "
                    "the second makes agreement meaningless."
                )
                continue
            if klass != "per_row_recovery":
                continue
            for benchmark_id in set(declared):
                targets.setdefault(benchmark_id, []).append(_rel(path))
    return targets, failures


def check_fit_then_score(verbose: bool) -> List[str]:
    declared, failures = _collect_declared_fit_targets(verbose)

    if verbose:
        print(f"    per-row fit targets: {sorted(declared) or 'none'}")

    if not PREDICTION_UNCERTAINTY.exists():
        failures.append(
            f"{_rel(PREDICTION_UNCERTAINTY)} is missing: the scored artifact this gate "
            "checks does not exist, so fit-then-score circularity cannot be ruled out."
        )
        return failures

    payload = _load_json(PREDICTION_UNCERTAINTY)
    rows = payload.get("benchmarks", []) if isinstance(payload, dict) else []
    summary = payload.get("summary", {}) if isinstance(payload, dict) else {}

    scored_flags: Dict[str, bool] = {}
    for row in rows:
        if not isinstance(row, dict):
            continue
        benchmark_id = str(row.get("benchmark_id", ""))
        if not benchmark_id:
            continue
        scored_flags[benchmark_id] = bool(row.get("fitted_row", False))

    # --- check 1: a declared fit target that is also scored must disclose it ----------
    for benchmark_id, records in sorted(declared.items()):
        if benchmark_id not in scored_flags:
            continue
        if not scored_flags[benchmark_id]:
            failures.append(
                f"CIRCULARITY: benchmark '{benchmark_id}' is named as a FIT TARGET by "
                f"{', '.join(sorted(set(records)))} and is ALSO scored as a row in "
                f"{_rel(PREDICTION_UNCERTAINTY)} without `fitted_row: true`. A constant "
                "fitted to a measurement reproduces that measurement; scoring it as "
                "coverage is circular. Either disclose the row (fitted_row: true, which "
                "removes it from the coverage numerator and denominator) or stop fitting "
                "against it."
            )

    # --- check 2: the coverage claim must exclude fitted rows -------------------------
    honest = summary.get("honest_literature_coverage") if isinstance(summary, dict) else None
    if not isinstance(honest, dict):
        failures.append(
            f"{_rel(PREDICTION_UNCERTAINTY)} has no `summary.honest_literature_coverage` "
            "block. The aggregate `ci_coverage_rate` pools literature rows with synthetic "
            "reproducibility rows and with fitted rows, so it cannot be the headline."
        )
    else:
        split = summary.get("signal_origin_split", {}) or {}
        fitted_bucket = split.get("fitted_row", {}) or {}
        fitted_total = int(fitted_bucket.get("total", 0) or 0)
        fitted_not_evaluable = int(fitted_bucket.get("not_evaluable", 0) or 0)
        declared_excluded = int(honest.get("excluded_fitted_rows", 0) or 0)
        if declared_excluded != fitted_total + fitted_not_evaluable:
            failures.append(
                "COVERAGE ACCOUNTING: `honest_literature_coverage.excluded_fitted_rows` is "
                f"{declared_excluded} but the fitted_row bucket holds "
                f"{fitted_total + fitted_not_evaluable} rows. Fitted rows must be listed "
                "separately, not folded in or dropped silently."
            )
        lit_bucket = split.get("external_literature", {}) or {}
        if int(honest.get("total", -1) or 0) != int(lit_bucket.get("total", -2) or 0):
            failures.append(
                "COVERAGE ACCOUNTING: `honest_literature_coverage.total` does not match the "
                "external_literature bucket's evaluable total. The honest denominator must "
                "be exactly the literature rows, with fitted and synthetic rows removed."
            )

    return failures


# ---------------------------------------------------------------------------
# check 4: constant precision
# ---------------------------------------------------------------------------
def _float_literals_with_precision(
    source: str, container_names: Sequence[str]
) -> List[Tuple[str, str, int, int]]:
    """(context, literal_text, lineno, sig_figs) for float literals in the containers."""
    tree = ast.parse(source)
    lines = source.splitlines()
    found: List[Tuple[str, str, int, int]] = []

    def literal_text(node: ast.AST) -> str:
        if not hasattr(node, "lineno"):
            return ""
        segment = ast.get_source_segment(source, node)
        if segment:
            return segment
        return lines[node.lineno - 1] if node.lineno - 1 < len(lines) else ""

    def scan(node: ast.AST, context: str, exempt: bool) -> None:
        if isinstance(node, ast.Call):
            for keyword in node.keywords:
                keyword_exempt = exempt or (keyword.arg in PRECISION_EXEMPT_KEYWORDS)
                scan(keyword.value, f"{context}.{keyword.arg or '*'}", keyword_exempt)
            for arg in node.args:
                scan(arg, context, exempt)
            return
        if isinstance(node, ast.Dict):
            for key, value in zip(node.keys, node.values):
                key_text = literal_text(key) if key is not None else "*"
                key_exempt = exempt or any(
                    token in str(key_text) for token in PRECISION_EXEMPT_KEYWORDS
                )
                scan(value, f"{context}[{key_text}]", key_exempt)
            return
        if isinstance(node, ast.Constant):
            if exempt or not isinstance(node.value, float):
                return
            text = literal_text(node)
            digits = _significant_figures(text)
            found.append((context, text, node.lineno, digits))
            return
        for child in ast.iter_child_nodes(node):
            scan(child, context, exempt)

    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        names = [t.id for t in node.targets if isinstance(t, ast.Name)]
        for name in names:
            if name in container_names:
                scan(node.value, name, False)
    return found


def check_constant_precision(verbose: bool) -> List[str]:
    failures: List[str] = []
    for rel_path, containers in CONSTANT_SOURCES.items():
        path = ROOT / rel_path
        if not path.exists():
            failures.append(f"constant source {rel_path} does not exist")
            continue
        source = path.read_text(encoding="utf-8")
        for context, text, lineno, digits in _float_literals_with_precision(source, containers):
            if digits <= MAX_SIGNIFICANT_FIGURES:
                continue
            failures.append(
                f"PRECISION: {rel_path}:{lineno} ships `{text}` ({digits} significant "
                f"figures, limit {MAX_SIGNIFICANT_FIGURES}) in {context}. No measurement "
                "in this domain supports that precision -- a constant written this way is "
                "the residue of a solve. Round it to a precision it can defend and keep "
                "the exact prior value in `previous_value=` (exempt from this check), "
                "with the fit disclosed."
            )
            if verbose:
                print(f"    flagged {rel_path}:{lineno} {text} ({digits} sf)")
    return failures


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    print("fit_target_gate: checking fit-target circularity and constant precision")

    print("  [1/2] fit-then-score circularity + coverage accounting")
    failures = check_fit_then_score(args.verbose)

    print("  [2/2] shipped-constant precision")
    failures += check_constant_precision(args.verbose)

    if failures:
        print()
        print(f"fit_target_gate: FAIL ({len(failures)} violation(s))")
        for failure in failures:
            print()
            print(f"  - {failure}")
        print()
        return 1

    print()
    print("fit_target_gate: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
