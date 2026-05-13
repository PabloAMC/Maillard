#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.pipeline import FormulationResult, UncertaintyEnvelope
from src.reporting import generate_comparison_report, generate_report


def _example_result(name: str, *, scale: float, risk: float) -> FormulationResult:
    predicted = {
        "Hexanal": (14.0 * scale, 10.0 * scale, 14.0 * scale, 18.0 * scale, "aldehyde", False),
        "Nonanal": (7.0 * scale, 5.5 * scale, 7.0 * scale, 9.5 * scale, "aldehyde", False),
        "2-Methyl-3-furanthiol (MFT)": (6.0 * scale, 3.0 * scale, 6.0 * scale, 9.0 * scale, "sulfur", True),
        "2-furfurylthiol": (4.0 * scale, 2.2 * scale, 4.0 * scale, 6.5 * scale, "sulfur", True),
        "2,5-dimethylpyrazine": (2.2 * scale, 1.5 * scale, 2.2 * scale, 3.0 * scale, "pyrazine", False),
    }
    envelopes = {
        compound: UncertaintyEnvelope(
            compound=compound,
            predicted_ppb=p50,
            predicted_p5=p5,
            predicted_p50=p50,
            predicted_p95=p95,
            ci_level_pct=90,
            support_count=3,
        )
        for compound, (ppb, p5, p50, p95, _klass, _direct_anchor) in predicted.items()
    }
    projection_metadata = {}
    confidence_rows = []
    for compound, (observable_ppb, p5, p50, p95, volatile_class, direct_anchor) in predicted.items():
        projection_metadata[compound.lower()] = {
            "compound": compound,
            "observable_ppb": observable_ppb,
            "volatile_class": volatile_class,
            "calibration_source": "Pratap-Singh 2021 soy-vs-pea ambient slurry release ratio" if direct_anchor else "class_fallback",
            "calibration_evidence_strength": "literature_anchored" if direct_anchor else "heuristic",
            "calibration_fallback_mode": "compound_specific" if direct_anchor else "class_level",
            "direct_anchor": direct_anchor,
        }
        confidence_rows.append(
            {
                "compound": compound,
                "observable_ppb": observable_ppb,
                "tier": "high" if direct_anchor else "medium",
                "score": 88.0 if direct_anchor else 55.0,
                "prediction_mode": "bounded_estimate" if direct_anchor else "directional_hypothesis",
                "reachability_status": "benchmark_aligned" if direct_anchor else "merely_plausible",
                "calibration_source": projection_metadata[compound.lower()]["calibration_source"],
                "observable_assumption_summary": "matrix-calibrated" if direct_anchor else "class-level surrogate",
            }
        )

    return FormulationResult(
        name=name,
        target_score=3.8 * scale,
        off_flavour_risk=risk,
        predicted_ppb={compound: values[2] for compound, values in predicted.items()},
        uncertainty_envelopes=envelopes,
        projection_metadata=projection_metadata,
        confidence_metadata={
            "tier": "high_confidence" if scale >= 1.0 else "mixed_confidence",
            "score": 82.0 if scale >= 1.0 else 61.0,
            "benchmark_neighborhood": "matrix_supported",
            "decision_mode": "bounded_estimate",
            "prediction_mode": "bounded_estimate",
            "recommended_posture": "scientist-facing report",
            "compound_confidence": confidence_rows,
        },
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate reproducible report visual examples for docs assets.")
    parser.add_argument("--output-dir", default="results/validation/report_visual_examples")
    parser.add_argument("--docs-asset-dir", default="docs/assets")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    docs_asset_dir = Path(args.docs_asset_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    docs_asset_dir.mkdir(parents=True, exist_ok=True)

    baseline = _example_result("baseline_meaty_mix", scale=1.0, risk=0.82)
    current = _example_result("glutathione_enriched_mix", scale=1.22, risk=0.58)
    # Force the aldehydes down relative to the sulfur families so the waterfall is legible.
    current.uncertainty_envelopes["Hexanal"].predicted_ppb = 7.0
    current.uncertainty_envelopes["Hexanal"].predicted_p5 = 5.0
    current.uncertainty_envelopes["Hexanal"].predicted_p50 = 7.0
    current.uncertainty_envelopes["Hexanal"].predicted_p95 = 10.0
    current.predicted_ppb["Hexanal"] = 7.0
    current.projection_metadata["hexanal"]["observable_ppb"] = 7.0

    current.uncertainty_envelopes["Nonanal"].predicted_ppb = 4.5
    current.uncertainty_envelopes["Nonanal"].predicted_p5 = 3.2
    current.uncertainty_envelopes["Nonanal"].predicted_p50 = 4.5
    current.uncertainty_envelopes["Nonanal"].predicted_p95 = 6.0
    current.predicted_ppb["Nonanal"] = 4.5
    current.projection_metadata["nonanal"]["observable_ppb"] = 4.5

    single_dir = generate_report(
        current,
        [],
        {
            "name": current.name,
            "protein_type": "pea_iso",
            "process_state": "aqueous_pre_extrusion_model",
            "target": "meaty",
            "minimize": "beany",
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

    shutil.copyfile(single_dir / "compound_confidence_overlay.png", docs_asset_dir / "report_compound_confidence_overlay.png")
    shutil.copyfile(single_dir / "intervention_waterfall.png", docs_asset_dir / "report_intervention_waterfall.png")
    shutil.copyfile(comparison_dir / "comparison_intervention_waterfall.png", docs_asset_dir / "report_comparison_intervention_waterfall.png")

    print(f"Wrote {single_dir / 'compound_confidence_overlay.png'}")
    print(f"Wrote {single_dir / 'intervention_waterfall.png'}")
    print(f"Wrote {comparison_dir / 'comparison_intervention_waterfall.png'}")
    print(f"Copied docs assets into {docs_asset_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())