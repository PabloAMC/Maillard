"""
Build Wave B2.4 -- the BLIND RE-SIT, run once per declared weighting.

Two scorers, both EXISTING and both reused rather than forked:

  * the full B2.1/B2.2/B2.3 hold-out panel,
    `generate_kinetic_core_b2_3_holdout`, and
  * the both-ways cutover exam, `generate_cutover_final_exam`.

Neither is re-implemented here. This file rebinds exactly two things before
calling them -- which frozen fit report the parameters come from, and which
basename the artifact is written under -- and does nothing else. A forked
hold-out scorer is a hold-out scorer that can drift, and the whole value of a
blind re-sit is that the scorer did not move.

FIREWALL. `data/benchmarks/external_validation/` is opened by
`generate_cutover_final_exam` and by nothing else, here or anywhere in this
wave. This file reads no measured value: it reads the two scorers' own output.

Usage
-----
    python scripts/generators/generate_kinetic_core_b2_4_scorers.py \
        --weight-tag shipped
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths
if str(ROOT / "scripts" / "generators") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))

WEIGHT_TAGS = ("shipped", "half", "measured")


#: Set by main() when --fit-report overrides the per-weighting default. Used
#: for the per-MEMBER exam re-sits W-2 needs, and for the W-4 control column,
#: which is scored at B2.2's OWN frozen parameters.
_OVERRIDE: Path = None


def fit_report_for(tag: str) -> Path:
    if _OVERRIDE is not None:
        return _OVERRIDE
    return data_paths.VALIDATION_DIR / f"kinetic_core_b2_4_fit_{tag}.json"


def run_panel(tag: str) -> Dict[str, Any]:
    """The full hold-out panel, at this weighting's ensemble-best parameters."""
    import generate_kinetic_core_b2_3_holdout as PANEL

    PANEL.FIT_REPORT = fit_report_for(tag)
    PANEL.OUT_BASENAME = f"kinetic_core_b2_4_panel_{tag}"
    rc = PANEL.main()
    if rc != 0:
        raise SystemExit(f"panel returned {rc} for weighting {tag}")
    return json.loads(
        (data_paths.VALIDATION_DIR / f"{PANEL.OUT_BASENAME}.json").read_text())


def run_exam(tag: str, basename: str = None) -> Dict[str, Any]:
    """
    The both-ways cutover exam, at this weighting's ensemble-best parameters.

    The engine resolves its sulfur parameters through a module-level path,
    `engine._B2_FIT_REPORT`, read at call time by `core_parameters` and
    `core_ph_drift`. Rebinding it is the whole of the intervention: no engine
    behaviour changes, no exam row changes, no band changes.
    """
    import generate_cutover_final_exam as EXAM
    from src.kinetic_core import engine

    previous = engine._B2_FIT_REPORT
    try:
        engine._B2_FIT_REPORT = fit_report_for(tag)
        payload = EXAM.build_exam(target_tag="meaty")
    finally:
        engine._B2_FIT_REPORT = previous

    payload["parameters_from"] = str(fit_report_for(tag).relative_to(ROOT))
    payload["b2_4_weight_tag"] = tag
    base = basename or f"kinetic_core_b2_4_exam_{tag}"
    out_json = data_paths.VALIDATION_DIR / f"{base}.json"
    out_json.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str))
    (data_paths.VALIDATION_DIR / f"{base}.md").write_text(
        EXAM.render_markdown(payload))
    return payload


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--weight-tag", choices=WEIGHT_TAGS, required=True)
    parser.add_argument("--panel-only", action="store_true")
    parser.add_argument("--exam-only", action="store_true")
    parser.add_argument("--fit-report", default=None,
                        help="score at THESE frozen parameters instead of the "
                             "weighting's ensemble-best (W-2's per-member "
                             "re-sits, and W-4's B2.2 control column)")
    parser.add_argument("--basename", default=None)
    args = parser.parse_args(argv)

    global _OVERRIDE
    if args.fit_report:
        _OVERRIDE = Path(args.fit_report)
        if not _OVERRIDE.is_absolute():
            _OVERRIDE = ROOT / args.fit_report
    tag = args.weight_tag
    if not fit_report_for(tag).exists():
        raise SystemExit(
            f"{fit_report_for(tag)} not found -- run the B2.4 fit and its "
            f"--consolidate step first. The scorers never fit anything.")

    if not args.exam_only:
        panel = run_panel(tag)
        sc = panel["scorecard"]
        print(f"[{tag}] PANEL gating {sc['gating_passed']}/{sc['gating_rows']}  "
              f"median|log10| {sc.get('median_abs_log10_fold_gating')}  "
              f"geo {sc.get('geometric_mean_fold_gating')}")
    if not args.panel_only:
        exam = run_exam(tag, args.basename)
        s = exam["summary"]
        print(f"[{tag}] EXAM paired median "
              f"{s['paired_subset']['core']['median_fold_error']}  "
              f"geo-mean {s['core'].get('geometric_mean_fold')}  "
              f"in band {s['core_within_band']}/{s['core_scored']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
