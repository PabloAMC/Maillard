#!/usr/bin/env python3
"""
Run a shareable multi-formulation campaign from a YAML specification.

This creates:
- per-formulation report bundles
- a side-by-side comparison artifact
- a campaign-level Markdown/JSON summary with provenance
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.conditions import ReactionConditions  # noqa: E402
from src.inverse_design import InverseDesigner  # noqa: E402
from src.reporting import generate_campaign_report, generate_report  # noqa: E402
from src.usability_reports import DomainOfValidityChecker, build_confidence_package  # noqa: E402


def _slugify(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_") or "run"


def _merge_formulation(base: Dict[str, object], overrides: Dict[str, object]) -> Dict[str, object]:
    merged = dict(base)
    for key, value in overrides.items():
        if key == "molar_ratios" and isinstance(value, dict):
            existing = dict(merged.get("molar_ratios", {}))
            existing.update(value)
            merged[key] = existing
        else:
            merged[key] = value
    return merged


def _load_campaign_spec(path: Path) -> Dict[str, object]:
    with open(path, "r") as handle:
        payload = yaml.safe_load(handle) or {}
    if not payload.get("campaign"):
        raise ValueError("Campaign spec must define a 'campaign' section.")
    if not payload.get("formulations"):
        raise ValueError("Campaign spec must define at least one formulation entry.")
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description="Run a shareable Maillard screening campaign.")
    parser.add_argument("--spec", required=True, help="Path to a campaign YAML specification")
    parser.add_argument("--output-dir", default=None, help="Directory where campaign artifacts should be written")
    args = parser.parse_args()

    spec_path = Path(args.spec).resolve()
    campaign_spec = _load_campaign_spec(spec_path)
    campaign_meta = campaign_spec.get("campaign", {})
    shared_conditions = campaign_spec.get("shared_conditions", {}) or {}

    target_tag = str(campaign_meta.get("target_tag", "meaty"))
    minimize_tag = str(campaign_meta.get("minimize_tag", "beany"))
    designer = InverseDesigner(target_tag, minimize_tag)
    checker = DomainOfValidityChecker(target_tag)

    base_conditions = ReactionConditions(
        pH=float(shared_conditions.get("ph", 6.0)),
        temperature_celsius=float(shared_conditions.get("temp", 150.0)),
        water_activity=float(shared_conditions.get("aw", 1.0)),
        protein_type=str(shared_conditions.get("protein_type", "free")),
    )

    requested_formulations: List[Dict[str, object]] = []
    for entry in campaign_spec.get("formulations", []):
        entry = entry or {}
        formulation_name = entry.get("name")
        if not formulation_name:
            raise ValueError("Each campaign formulation must define a name.")
        base = next((item for item in designer.grid if item.get("name") == formulation_name), None)
        if base is None:
            raise ValueError(f"Formulation '{formulation_name}' not found in data/formulation_grid.yml")
        requested_formulations.append(_merge_formulation(base, entry.get("overrides", {})))

    print("======================================================")
    print("      Maillard Shareable Campaign Runner")
    print("======================================================\n")
    print(f"Campaign:  {campaign_meta.get('name', 'Unnamed campaign')}")
    print(f"Target:    {target_tag}")
    print(f"Minimize:  {minimize_tag}")
    print(f"Spec:      {spec_path}")
    print(f"Runs:      {len(requested_formulations)}")
    print("-" * 54)

    results = designer.evaluate_all(base_conditions, grid_override=requested_formulations)
    results_by_name = {item.name: item for item in results}
    ordered_results = [results_by_name[item["name"]] for item in requested_formulations if item["name"] in results_by_name]

    warnings_list = []
    for formulation in requested_formulations:
        result = results_by_name.get(formulation["name"])
        if result is None:
            continue
        precursor_names: List[str] = []
        for key in ("sugars", "amino_acids", "additives", "lipids"):
            precursor_names.extend(formulation.get(key, []))
        protein_type = str(formulation.get("protein_type", shared_conditions.get("protein_type", "free")))
        item_warnings = checker.check(
            precursor_names=precursor_names,
            protein_type=protein_type,
            temp_c=float(formulation.get("temp", shared_conditions.get("temp", base_conditions.temperature_celsius))),
            ph=float(formulation.get("ph", shared_conditions.get("ph", base_conditions.pH))),
        )
        result.confidence_metadata = build_confidence_package(
            result,
            item_warnings,
            precursor_names=precursor_names,
            protein_type=protein_type,
            formulation=formulation,
            baseline_conditions=base_conditions,
            designer=designer,
        )
        warnings_list.append(item_warnings)

    output_dir = Path(args.output_dir) if args.output_dir else Path("reports") / f"campaign_{_slugify(str(campaign_meta.get('name', 'campaign')))}"
    output_dir.mkdir(parents=True, exist_ok=True)

    run_artifacts = []
    for formulation, result, item_warnings in zip(requested_formulations, ordered_results, warnings_list):
        run_dir = output_dir / "runs" / _slugify(result.name)
        generate_report(
            result=result,
            warnings=item_warnings,
            conditions_dict=formulation,
            output_dir=run_dir,
            campaign_metadata=campaign_meta,
        )
        run_artifacts.append(
            {
                "name": result.name,
                "directory": str(run_dir),
                "markdown": str(run_dir / "report.md"),
                "json": str(run_dir / "report.json"),
            }
        )

    campaign_dir = generate_campaign_report(
        campaign_spec=campaign_spec,
        results=ordered_results,
        conditions_list=requested_formulations,
        warnings_list=warnings_list,
        run_artifacts=run_artifacts,
        output_dir=output_dir,
    )

    print(f"\nCampaign artifacts written to: {campaign_dir}")
    print(f"  - {campaign_dir / 'campaign.md'}")
    print(f"  - {campaign_dir / 'campaign.json'}")
    print(f"  - {campaign_dir / 'comparison.md'}")
    print(f"  - {campaign_dir / 'comparison.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())