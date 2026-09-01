#!/usr/bin/env python3
"""Generate the report figures and sample output the README shows — from a REAL run.

2026-08-27 (Wave I) — REWRITTEN. Read this before changing it back.
--------------------------------------------------------------------
This generator used to build two `FormulationResult` objects by hand, from hard-coded
numbers, and render them. Nothing in it ever touched the model. The consequences were not
cosmetic:

* The artifact it produced announced **"DECISION CONFIDENCE: HIGH_CONFIDENCE (82/100)"**
  and **"SCIENTIFIC ENVELOPE: TRUSTED"** — both hard-coded — while the repository's own
  README states there is currently *no* "high" tier and that every matrix prediction comes
  out `exploratory`. A reader comparing the two would conclude the README was being
  pessimistic. It was not; the sample was fabricated.
* It labelled 2-methyl-3-furanthiol and 2-furfurylthiol `literature_anchored`, tier
  `high`, score 88, citing *"Pratap-Singh 2021 soy-vs-pea ambient slurry release ratio"* as
  their calibration source. That lane is an aldehyde/furan **matrix observability** anchor
  with no sulfur content whatsoever, and it is now labelled `fitted_to_benchmark` because
  its factors were back-solved from the benchmarks they are scored against. So the
  showcase attributed the flagship compounds to an anchor that is neither about them nor
  independent.
* Because the numbers were invented, the figures could not drift when the model changed.
  The tracked PNGs dated from May and had survived every subsequent recalibration intact —
  the one artifact in the repository guaranteed never to show a regression.

The generator now runs the pipeline. Whatever the model currently predicts is what the
README shows, including tiers that say `exploratory` and confidence lines that say the run
is not trustworthy. If that looks bad, the fix is the model, not the fixture.

Determinism: the two formulations, the conditions and the protein type are pinned below,
and the pipeline is deterministic under a fixed `PYTHONHASHSEED` (see QUICKSTART's
Reproducibility section), so re-running this reproduces the same figures until the science
changes — which is exactly when they should change.

LaTeX: `src/plot_style.py` falls back to mathtext when dvipng is absent (Wave I), so this
runs on the documented conda path as well as in Docker. Set `MAILLARD_STRICT_LATEX=1` to
require the real toolchain.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths
from src.conditions import ReactionConditions
from src.pipeline import MaillardPipeline
from src.reporting import generate_comparison_report, generate_report
from src.usability_reports import prepare_cli_confidence

# Two formulations built exactly the way the documented forward-mode CLI builds them
# (docs/guides/QUICKSTART.md, "First Prediction"), so the report this produces is the same
# artifact a reader gets by running the command in the docs. The pair differs only by
# cysteine, so the intervention waterfall shows a real effect rather than a staged one.
BASELINE_FORMULATION = {
    "name": "Ribose + leucine (no cysteine)",
    "sugars": ["ribose"],
    "amino_acids": ["leucine"],
    "additives": [],
    "lipids": [],
    "molar_ratios": {"ribose": 0.5, "leucine": 0.1},
}
CURRENT_FORMULATION = {
    "name": "Ribose + cysteine + leucine",
    "sugars": ["ribose"],
    "amino_acids": ["cysteine", "leucine"],
    "additives": [],
    "lipids": [],
    "molar_ratios": {"ribose": 0.5, "cysteine": 0.2, "leucine": 0.1},
}

# ReactionConditions has no time field; time lives on the formulation dict the
# pipeline evaluates (`time_minutes`).
CONDITIONS = dict(
    temperature_celsius=105.0,
    pH=5.5,
    water_activity=0.85,
    protein_type="pea_iso",
)
TIME_MINUTES = 45.0


def _run(pipeline: MaillardPipeline, formulation: dict, conditions: ReactionConditions):
    payload = dict(formulation)
    payload.setdefault("time_minutes", TIME_MINUTES)
    # Same entry point the forward-mode CLI uses (scripts/run_pipeline.py:508)...
    result = pipeline.evaluate_single(payload, conditions)
    # ...followed by the same confidence step the CLI applies before rendering
    # (scripts/run_pipeline.py:155). Without it `confidence_metadata` is empty and the
    # report omits the Compound Confidence table entirely -- which is precisely the table
    # the README quotes, so skipping it would put the showcase back out of step with the
    # tool. Doing this here keeps the artifact identical to what the documented command
    # produces.
    precursors = (
        list(payload.get("sugars", []))
        + list(payload.get("amino_acids", []))
        + list(payload.get("additives", []))
        + list(payload.get("lipids", []))
    )
    warnings = prepare_cli_confidence(
        result,
        target_tag=pipeline.target_tag,
        precursor_names=precursors,
        protein_type=str(CONDITIONS["protein_type"]),
        temp_c=float(CONDITIONS["temperature_celsius"]),
        ph=float(CONDITIONS["pH"]),
        aw=float(CONDITIONS["water_activity"]),
        formulation=payload,
        baseline_conditions=conditions,
        designer=pipeline,
    )
    return result, warnings


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=data_paths.rel(data_paths.REPORT_VISUAL_EXAMPLES_DIR))
    parser.add_argument("--docs-asset-dir", default=data_paths.rel(data_paths.DOCS_ASSETS_DIR))
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    docs_asset_dir = Path(args.docs_asset_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    docs_asset_dir.mkdir(parents=True, exist_ok=True)

    pipeline = MaillardPipeline(target_tag="meaty", minimize_tag="beany")
    conditions = ReactionConditions(**CONDITIONS)

    baseline, _baseline_warnings = _run(pipeline, BASELINE_FORMULATION, conditions)
    current, current_warnings = _run(pipeline, CURRENT_FORMULATION, conditions)

    single_dir = generate_report(
        current,
        current_warnings,
        {
            "name": current.name,
            "protein_type": "pea_iso",
            "process_state": "aqueous_pre_extrusion_model",
            "target": "meaty",
            "minimize": "beany",
            "temp": CONDITIONS["temperature_celsius"],
            "ph": CONDITIONS["pH"],
            "time_minutes": TIME_MINUTES,
            "aw": CONDITIONS["water_activity"],
            "generated_by": "scripts/generators/generate_report_visual_examples.py",
            "note": (
                "REAL pipeline output, not a rendering fixture. Regenerate with "
                "scripts/generators/generate_report_visual_examples.py."
            ),
        },
        output_dir=output_dir / "single_run",
        baseline_result=baseline,
    )

    comparison_dir = generate_comparison_report(
        [baseline, current],
        [
            {"name": baseline.name, "protein_type": "pea_iso", "target": "meaty"},
            {"name": current.name, "protein_type": "pea_iso", "target": "meaty"},
        ],
        output_dir=output_dir / "comparison_run",
        campaign_metadata={"campaign_name": "report_visual_examples"},
    )

    copied = []
    for source, destination in (
        (single_dir / "compound_confidence_overlay.png", "report_compound_confidence_overlay.png"),
        (single_dir / "intervention_waterfall.png", "report_intervention_waterfall.png"),
        (comparison_dir / "comparison_intervention_waterfall.png", "report_comparison_intervention_waterfall.png"),
    ):
        if source.exists():
            shutil.copyfile(source, docs_asset_dir / destination)
            copied.append(destination)
            print(f"Wrote {source}")
        else:
            print(f"NOT WRITTEN (figure path unavailable): {source}")

    print(f"Copied {len(copied)} docs asset(s) into {docs_asset_dir}")
    print()
    print(f"Baseline: {baseline.name}  target_score={baseline.target_score:.3f}")
    print(f"Current : {current.name}  target_score={current.target_score:.3f}")
    print(
        "Run tier: "
        f"{(current.confidence_metadata or {}).get('tier', 'unknown')} "
        f"(score {(current.confidence_metadata or {}).get('score', 'n/a')})"
    )
    print("Reminder: the README sample table must be copied from this run, not written by hand.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
