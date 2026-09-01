#!/usr/bin/env python3
"""Score BOTH frozen hold-outs before and after the 2026-08-28 barrier-offset retirement.

WHAT WAS RETIRED. `data/lit/refinement_surrogate_patches.json` carried AUTO-ACCEPTED
barrier offsets of exactly +/-3.0 kcal/mol which `src.barrier_constants.get_barrier()`
applied to EVERY shipped prediction, silently overriding the audited `FAST_BARRIERS`
table.  Measured shipped values before retirement: `Schiff_Base_Formation` 18.0 (table
15.0), `Retro_Aldol_Fragmentation` 29.0 (table 32.0), `Thiol_Addition` 31.6 (table 28.6).
At 150 C, 3.0 kcal/mol is a ~35x rate factor.  The offsets had been accepted BECAUSE they
improved the benchmark panel the model is then scored on -- an undeclared fit to the
evaluation set -- and every one sat exactly at the +/-3.0 search bound, which is a bound
report and not an optimum.  `accepted_offsets` is now permanently empty.

WHY THIS SCRIPT EXISTS.  The two frozen hold-outs are the only surfaces in the repository
that were never available to any fit, so they are the only surfaces on which "the numbers
got worse" is an honest measurement rather than a re-scoring of a fit against itself.

HOW THE "BEFORE" IS RECONSTRUCTED.  NOT by restoring the retired file.  `get_barrier()`
merges a `BARRIER_OFFSETS` environment variable on top of the (now empty) patch file, so
setting that variable to the nine retired offsets -- kept for provenance in the patch
file's own `retirement_note.retired_offsets_kept_for_provenance` -- reproduces the
pre-retirement shipped barriers exactly, in a SUBPROCESS, without any file being touched.
The two arms run in separate interpreters so no module-level cache can leak between them.

NEITHER HOLD-OUT IS REGENERATED OR RE-FROZEN.
  * `results/validation/maillard_path_holdout_frozen_predictions.json` is a
    PRE-REGISTRATION at git 12f43dd.  It is READ.  It is never written by this script.
  * The matrix eight-point hold-out lives under `data/benchmarks/external_validation/`,
    which `scripts/ci/holdout_guard.py` asserts statically is named by no fit record.
Both are scored; nothing about them enters any fit.

Usage:
    python scripts/generators/generate_wave_r1_holdout_prepost.py
    python scripts/generators/generate_wave_r1_holdout_prepost.py --emit   # internal
"""

from __future__ import annotations

import argparse
import json
import math
import os
import statistics
import subprocess
import sys
from datetime import date
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths

VALIDATION = data_paths.VALIDATION_DIR
FROZEN = VALIDATION / "maillard_path_holdout_frozen_predictions.json"
PATCH_FILE = data_paths.REFINEMENT_SURROGATE_PATCHES
FROZEN_COMMIT = "12f43dd"

OUT_JSON = VALIDATION / "holdout_prepost_barrier_offset_retirement.json"
OUT_MD = VALIDATION / "holdout_prepost_barrier_offset_retirement.md"


def retired_offsets() -> Dict[str, float]:
    """The nine offsets that were live before the retirement, read from the patch file's
    own provenance block -- so this script cannot drift from what was actually retired."""
    payload = json.loads(PATCH_FILE.read_text(encoding="utf-8"))
    if payload.get("accepted_offsets"):
        raise SystemExit(
            "accepted_offsets is NON-EMPTY: the retirement has been undone. Refusing to "
            "report a before/after against a live fit."
        )
    note = payload.get("retirement_note", {})
    offsets = note.get("retired_offsets_kept_for_provenance", {})
    if not offsets:
        raise SystemExit("no retired offsets recorded in the patch file's retirement note")
    return {str(k): float(v) for k, v in offsets.items()}


# --------------------------------------------------------------------------------------
# the two hold-out surfaces, measured inside one arm
# --------------------------------------------------------------------------------------


def _fold_stats(folds: List[float]) -> Dict[str, Any]:
    if not folds:
        return {"n": 0}
    errs = [f if f >= 1.0 else 1.0 / f for f in folds]
    dex = [abs(math.log10(f)) for f in folds]
    return {
        "n": len(folds),
        "median_fold_error": statistics.median(errs),
        "median_abs_log10": statistics.median(dex),
        "mean_abs_log10": statistics.fmean(dex),
        "worst_fold_error": max(errs),
        "within_10x": sum(1 for e in errs if e <= 10.0),
    }


def score_maillard_path() -> Dict[str, Any]:
    from scripts.generators.generate_maillard_path_holdout_frozen_predictions import (
        build_frozen_predictions,
    )

    payload = build_frozen_predictions(target_tag="meaty")
    rows: Dict[str, Dict[str, Any]] = {}
    for bench in payload.get("benchmarks", []):
        for compound in bench.get("compounds", []):
            rows[f"{bench['benchmark_id']} / {compound['compound']}"] = {
                "predicted_native_unit": compound.get("predicted_native_unit"),
                "target_value": compound.get("target_value"),
                "fold_error": compound.get("fold_error"),
            }
    return {"summary": dict(payload.get("summary", {})), "rows": rows}


def score_matrix_eight_point() -> Dict[str, Any]:
    from src.external_validation import build_external_validation_report

    payload = build_external_validation_report(n_samples=200, seed=0, target_tag="meaty")
    rows: Dict[str, Dict[str, Any]] = {}
    for bench in payload.get("benchmarks", []):
        bid = bench.get("benchmark_id") or bench.get("benchmark") or "?"
        for comparison in bench.get("compounds", []) or []:
            name = comparison.get("compound")
            if name is None:
                continue
            rows[f"{bid} / {name}"] = {
                "measured_ppb": comparison.get("measured_ppb"),
                # the MC lane reports a posterior; p50 is the point estimate the report's
                # own fold_error is built from
                "predicted_ppb": comparison.get("predicted_p50"),
                "fold_error": comparison.get("fold_error"),
                "within_ci": comparison.get("inside_ci"),
            }
    return {"summary": dict(payload.get("summary", {})), "rows": rows}


def emit_arm() -> Dict[str, Any]:
    from src.barrier_constants import get_barrier

    probes = [
        "Schiff_Base_Formation",
        "Retro_Aldol_Fragmentation",
        "Thiol_Addition",
        "Thiol_Addition_Pentodiulose",
        "Thiol_Addition_Hexose",
        "Amadori_Rearrangement",
    ]
    return {
        "barrier_offsets_env": os.environ.get("BARRIER_OFFSETS"),
        "shipped_barriers": {fam: get_barrier(fam)[0] for fam in probes},
        "maillard_path": score_maillard_path(),
        "matrix_eight_point": score_matrix_eight_point(),
    }


# --------------------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------------------


def run_arm(offsets: Optional[Dict[str, float]]) -> Dict[str, Any]:
    env = dict(os.environ)
    if offsets:
        env["BARRIER_OFFSETS"] = json.dumps(offsets)
    else:
        env.pop("BARRIER_OFFSETS", None)
    proc = subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), "--emit"],
        cwd=str(ROOT),
        env=env,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise SystemExit(f"arm failed (offsets={bool(offsets)})")
    marker = "<<<R1-ARM-JSON>>>"
    _, _, tail = proc.stdout.partition(marker)
    return json.loads(tail.strip())


def _cmp_rows(before: Dict[str, Any], after: Dict[str, Any], key: str) -> List[Dict[str, Any]]:
    out = []
    for name in sorted(set(before) | set(after)):
        b = (before.get(name) or {}).get(key)
        a = (after.get(name) or {}).get(key)
        moved = not (
            b is not None
            and a is not None
            and (b == a or abs(b - a) <= 1e-9 * max(abs(b), abs(a)))
        )
        if b is None and a is None:
            moved = False
        out.append({"target": name, "before": b, "after": a, "moved": moved,
                    "ratio": (a / b) if (b not in (None, 0) and a is not None) else None})
    return out


def build() -> Dict[str, Any]:
    frozen = json.loads(FROZEN.read_text(encoding="utf-8"))
    if not frozen.get("git", {}).get("short", "").startswith(FROZEN_COMMIT):
        raise SystemExit(
            f"the frozen artifact does not name commit {FROZEN_COMMIT}; refusing to score "
            "against a baseline that is not the pre-registration"
        )
    offsets = retired_offsets()
    before = run_arm(offsets)
    after = run_arm(None)

    frozen_rows = {}
    for bench in frozen.get("benchmarks", []):
        for compound in bench.get("compounds", []):
            frozen_rows[f"{bench['benchmark_id']} / {compound['compound']}"] = {
                "predicted_native_unit": compound.get("predicted_native_unit"),
                "target_value": compound.get("target_value"),
                "fold_error": compound.get("fold_error"),
            }

    mp_rows = _cmp_rows(
        before["maillard_path"]["rows"], after["maillard_path"]["rows"], "predicted_native_unit"
    )
    mx_rows = _cmp_rows(
        before["matrix_eight_point"]["rows"], after["matrix_eight_point"]["rows"], "predicted_ppb"
    )
    # against the 12f43dd pre-registration itself, over ITS OWN targets only
    vs_frozen_before = _cmp_rows(
        frozen_rows,
        {k: v for k, v in before["maillard_path"]["rows"].items() if k in frozen_rows},
        "predicted_native_unit",
    )
    vs_frozen_after = _cmp_rows(
        frozen_rows,
        {k: v for k, v in after["maillard_path"]["rows"].items() if k in frozen_rows},
        "predicted_native_unit",
    )

    return {
        "artifact": "holdout_prepost_barrier_offset_retirement",
        "wave": "R1",
        "generated_on": date.today().isoformat(),
        "what_changed": (
            "The auto-accepted +/-3.0 kcal/mol barrier offsets in "
            "data/lit/refinement_surrogate_patches.json were retired on 2026-08-28. They had "
            "been applied by src.barrier_constants.get_barrier() to every shipped prediction, "
            "silently overriding the audited FAST_BARRIERS table. The 'before' arm reconstructs "
            "them through the BARRIER_OFFSETS environment variable in a subprocess; no file was "
            "restored and no hold-out was regenerated or re-frozen."
        ),
        "barriers": {
            "before": before["shipped_barriers"],
            "after": after["shipped_barriers"],
        },
        "maillard_path": {
            "baseline_file": str(FROZEN.relative_to(ROOT)),
            "baseline_commit": frozen.get("git", {}).get("short"),
            "frozen_summary": dict(frozen.get("summary", {})),
            "summary_before": before["maillard_path"]["summary"],
            "summary_after": after["maillard_path"]["summary"],
            "targets_moved_by_the_retirement": sum(1 for r in mp_rows if r["moved"]),
            "targets_total": len(mp_rows),
            "targets_differing_from_the_freeze_before": sum(1 for r in vs_frozen_before if r["moved"]),
            "targets_differing_from_the_freeze_after": sum(1 for r in vs_frozen_after if r["moved"]),
            "frozen_target_count": len(frozen_rows),
            "rows": mp_rows,
        },
        "matrix_eight_point": {
            "summary_before": before["matrix_eight_point"]["summary"],
            "summary_after": after["matrix_eight_point"]["summary"],
            "targets_moved_by_the_retirement": sum(1 for r in mx_rows if r["moved"]),
            "targets_total": len(mx_rows),
            "rows": mx_rows,
        },
    }


def _fmt(value: Any, digits: int = 4) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def render(payload: Dict[str, Any]) -> str:
    out: List[str] = []
    out.append("# Both frozen hold-outs, before and after the barrier-offset retirement (2026-08-28)\n")
    out.append(payload["what_changed"] + "\n")

    out.append("## Shipped barriers (kcal/mol)\n")
    out.append("| family | before (offset applied) | after (= FAST_BARRIERS table) | delta |")
    out.append("|---|---:|---:|---:|")
    for fam, before_value in payload["barriers"]["before"].items():
        after_value = payload["barriers"]["after"][fam]
        out.append(f"| `{fam}` | {before_value} | {after_value} | {after_value - before_value:+.2f} |")
    out.append("")

    mp = payload["maillard_path"]
    out.append("## Hold-out 1 — the free-precursor `maillard_path` pre-registration\n")
    out.append(
        f"Baseline `{mp['baseline_file']}` at commit `{mp['baseline_commit']}`: READ, never "
        "regenerated, never re-frozen.\n"
    )
    out.append(
        f"**{mp['targets_moved_by_the_retirement']} of {mp['targets_total']} targets moved** when "
        "the offsets were retired.\n"
    )
    out.append("| metric | frozen (12f43dd) | before retirement | after retirement |")
    out.append("|---|---:|---:|---:|")
    fs, bs, asum = mp["frozen_summary"], mp["summary_before"], mp["summary_after"]
    for key in (
        "bundle_count",
        "target_count",
        "quantitatively_scored_count",
        "structural_zero_count",
        "median_fold_error",
        "median_abs_log10_error",
        "worst_fold_error",
        "best_fold_error",
        "within_10x",
        "ordinal_pairs_correct",
        "series_directions_correct",
    ):
        out.append(f"| `{key}` | {_fmt(fs.get(key))} | {_fmt(bs.get(key))} | {_fmt(asum.get(key))} |")
    out.append("")
    out.append(
        f"Against the pre-registration's own {mp['frozen_target_count']} targets: "
        f"**{mp['targets_differing_from_the_freeze_before']} differed before the retirement, "
        f"{mp['targets_differing_from_the_freeze_after']} differ after**.\n"
    )
    moved = [row for row in mp["rows"] if row["moved"]]
    if moved:
        out.append("### Targets the retirement moved\n")
        out.append("| target | before | after | after/before |")
        out.append("|---|---:|---:|---:|")
        for row in moved:
            out.append(
                f"| `{row['target']}` | {_fmt(row['before'])} | {_fmt(row['after'])} | "
                f"{_fmt(row['ratio'])}x |"
            )
    else:
        out.append("**No target moved.**\n")
    out.append("")

    mx = payload["matrix_eight_point"]
    out.append("## Hold-out 2 — the eight-point matrix hold-out\n")
    out.append(
        f"**{mx['targets_moved_by_the_retirement']} of {mx['targets_total']} points moved.**\n"
    )
    out.append("| metric | before retirement | after retirement |")
    out.append("|---|---:|---:|")
    for key in (
        "matched_compound_count",
        "ci_coverage_hits",
        "ci_coverage_rate",
        "median_accuracy_fold",
        "median_abs_log10_error",
        "max_fold_error",
    ):
        out.append(
            f"| `{key}` | {_fmt(mx['summary_before'].get(key))} | {_fmt(mx['summary_after'].get(key))} |"
        )
    out.append("")
    moved_mx = [row for row in mx["rows"] if row["moved"]]
    if moved_mx:
        out.append("| point | before ppb | after ppb | after/before |")
        out.append("|---|---:|---:|---:|")
        for row in moved_mx:
            out.append(
                f"| `{row['target']}` | {_fmt(row['before'])} | {_fmt(row['after'])} | "
                f"{_fmt(row['ratio'])}x |"
            )
    else:
        out.append(
            "**No point moved.** The matrix hold-out runs the `matrix_only` execution path, "
            "which never reaches the reaction network (Wave S1), so it is structurally blind to "
            "any barrier. That is a NEGATIVE result about the hold-out, not a clean bill of "
            "health for the model.\n"
        )
    out.append("")
    return "\n".join(out)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--emit", action="store_true", help="internal: run one arm and print JSON")
    args = parser.parse_args()
    if args.emit:
        print("<<<R1-ARM-JSON>>>")
        print(json.dumps(emit_arm()))
        return 0
    payload = build()
    OUT_JSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    OUT_MD.write_text(render(payload), encoding="utf-8")
    print(f"wrote {OUT_JSON.relative_to(ROOT)}")
    print(f"wrote {OUT_MD.relative_to(ROOT)}")
    mp, mx = payload["maillard_path"], payload["matrix_eight_point"]
    print(
        f"maillard_path: {mp['targets_moved_by_the_retirement']}/{mp['targets_total']} moved; "
        f"median fold {mp['summary_before'].get('median_fold_error')} -> "
        f"{mp['summary_after'].get('median_fold_error')}"
    )
    print(f"matrix 8-point: {mx['targets_moved_by_the_retirement']}/{mx['targets_total']} moved")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
