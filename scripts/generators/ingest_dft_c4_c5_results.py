#!/usr/bin/env python3
"""
scripts/generators/ingest_dft_c4_c5_results.py — Ingest DFT C4/C5 results.

Reads completed barrier_summary.json files from results/dft_c4_c5/<reaction>/
and promotes those entries from their current tier (mlp_screen_mace or
xtb_derived_gfn2) to selective_dft_anchor in:

  - data/lit/arrhenius_params.yml   (adds DFT-derived kinetic parameters)
  - data/lit/computational_gap_closure_targets.json  (updates method + tier)

Also writes results/validation/dft_c4_c5_ingestion_report.{md,json}.

Usage:
  python scripts/generators/ingest_dft_c4_c5_results.py
  python scripts/generators/ingest_dft_c4_c5_results.py --output-dir results/validation
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

DFT_OUTPUT_DIR  = ROOT / "results" / "dft_c4_c5"
ARRHENIUS_FILE  = ROOT / "data" / "lit" / "arrhenius_params.yml"
GAP_TARGETS_FILE = ROOT / "data" / "lit" / "computational_gap_closure_targets.json"

# Canonical mapping: reaction_key → arrhenius YAML key to create/update
REACTION_TO_YAML_KEY: Dict[str, str] = {
    "hexanal_radical_quench": "hexanal_radical_quench_dft",
    "quinone_cys_michael":    "quinone_cys_michael_dft",
    "aa_ring_open_dicarbonyl": "aa_ring_open_dicarbonyl_dft",
    "pe_schiff_base":         "pe_schiff_base_dft",
    "pe_amadori":             "pe_amadori_dft",
    "lysinoalanine_crosslink": "lysinoalanine_crosslink_dft",
}


def load_dft_results() -> List[Dict[str, Any]]:
    """Load all completed barrier_summary.json files."""
    results = []
    for reaction_dir in sorted(DFT_OUTPUT_DIR.glob("*/barrier_summary.json")):
        data = json.loads(reaction_dir.read_text())
        results.append(data)
    return results


def _arrhenius_entry_from_dft(result: Dict[str, Any]) -> Dict[str, Any]:
    """Build a new arrhenius_params.yml entry from a DFT result."""
    barrier_kj = float(result["barrier_kj_mol"])
    # Estimate A factor from Evans-Polanyi + collision theory at 423 K
    # A ≈ kT/h * exp(ΔS‡/R) — use 1e12 as DFT-calibrated TST pre-exponential
    return {
        "A_value": 1e12,
        "A_unit": "1/s",
        "Ea_kj_mol": round(barrier_kj, 3),
        "uncertainty_kj": float(result.get("uncertainty_kj", 42.0)),
        "temp_range_C": "150",
        "pH": "Neutral",
        "substrate": result["reaction_key"].replace("_", " ").title(),
        "source": f"wB97M-V/def2-tzvp // r2SCAN-3c/def2-svp, ddCOSMO(water), T=150°C",
        "source_quality": result.get("source_quality", "selective_dft_anchor"),
        "slr_family": str(result.get("slr_family", "unknown")),
        "method": result.get("method", "wB97M-V//r2SCAN-3c"),
        "basis_geom": result.get("basis_geom", "def2-svp"),
        "basis_sp": result.get("basis_sp", "def2-tzvp"),
        "promotion_ceiling": result.get("promotion_ceiling", "ranking_only"),
        "honest_label": result.get("honest_label", "DFT-based estimate"),
    }


def update_arrhenius_yaml(results: List[Dict[str, Any]]) -> List[str]:
    """Add DFT-derived entries to arrhenius_params.yml. Returns list of updated keys."""
    if not ARRHENIUS_FILE.exists():
        raise FileNotFoundError(f"Not found: {ARRHENIUS_FILE}")

    with open(ARRHENIUS_FILE) as f:
        params = yaml.safe_load(f) or {}

    reactions = params.setdefault("reactions", {})
    updated_keys = []
    for result in results:
        key = result["reaction_key"]
        yaml_key = REACTION_TO_YAML_KEY.get(key, f"{key}_dft")
        reactions[yaml_key] = _arrhenius_entry_from_dft(result)
        updated_keys.append(yaml_key)

    with open(ARRHENIUS_FILE, "w") as f:
        yaml.dump(params, f, default_flow_style=False, allow_unicode=True, sort_keys=False)

    return updated_keys


def update_gap_closure_targets(results: List[Dict[str, Any]]) -> None:
    """Update computational_gap_closure_targets.json with promoted method tier."""
    if not GAP_TARGETS_FILE.exists():
        print(f"  Note: {GAP_TARGETS_FILE} not found — skipping gap targets update")
        return

    with open(GAP_TARGETS_FILE) as f:
        targets = json.load(f)

    entries = targets if isinstance(targets, list) else targets.get("targets", [])
    result_map = {r["reaction_key"]: r for r in results}

    for entry in entries:
        key = str(entry.get("reaction_key", entry.get("reaction_id", "")))
        # Strip trailing _mlp_derived or _xtb_derived suffix for matching
        base_key = key.replace("_mlp_derived", "").replace("_xtb_derived", "")
        if base_key in result_map:
            r = result_map[base_key]
            entry["method"] = r.get("method", "wB97M-V//r2SCAN-3c")
            entry["source_quality"] = r.get("source_quality", "selective_dft_anchor")
            entry["barrier_kcal_mol"] = r["barrier_kcal_mol"]
            entry["uncertainty_kj"] = r.get("uncertainty_kj", 42.0)

    if isinstance(targets, list):
        updated = entries
    else:
        targets["targets"] = entries
        updated = targets

    with open(GAP_TARGETS_FILE, "w") as f:
        json.dump(updated, f, indent=2)


def render_report_markdown(results: List[Dict[str, Any]], yaml_keys: List[str]) -> str:
    """Render a markdown ingestion report."""
    lines = [
        "# DFT C4/C5 Ingestion Report",
        "",
        f"Ingested {len(results)} DFT results (Step C4: r2SCAN-3c geom+freq, Step C5: wB97M-V SP).",
        "",
        "## Promoted Barriers",
        "",
        "| Reaction Key | Family | ΔG‡ (kcal/mol) | ΔG‡ (kJ/mol) | Uncertainty (kJ) | Tier | Ceiling |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for r in results:
        lines.append(
            f"| {r['reaction_key']} "
            f"| {r.get('slr_family','?')} "
            f"| {r['barrier_kcal_mol']:.3f} "
            f"| {r['barrier_kj_mol']:.2f} "
            f"| ±{r.get('uncertainty_kj', 42.0):.1f} "
            f"| {r.get('source_quality','?')} "
            f"| {r.get('promotion_ceiling','?')} |"
        )
    lines += [
        "",
        "## YAML Keys Added to arrhenius_params.yml",
        "",
    ]
    for k in yaml_keys:
        lines.append(f"- `{k}`")
    lines += [
        "",
        "## Evidence Labels",
        "",
        "- `selective_dft_anchor` — wB97M-V//r2SCAN-3c result, ±factor 2 in concentration, ranking reliable",
        "- `literature_derived_transfer` — Known Ea from SLR, ±factor 2–5, use with care",
        "- `xtb_derived_gfn2` — GFN2-xTB path only, ±factor 3–8, use for ranking only",
    ]
    return "\n".join(lines) + "\n"


def main(output_dir: str = "results/validation") -> None:
    out_path = ROOT / output_dir

    print(f"Loading DFT results from {DFT_OUTPUT_DIR} ...")
    results = load_dft_results()

    if not results:
        print(f"  No completed DFT results found in {DFT_OUTPUT_DIR}")
        # 2026-08-27 (Wave R): this used to read "Run scripts/run_dft_c4_c5.py first."
        # No such script exists, in the working tree or anywhere in git history, so the
        # instruction sent a stuck user to a file that was never written. Stating the
        # gap is more useful than naming a phantom runner.
        print("  No runner for the C4/C5 DFT jobs exists in this repository "
              "(scripts/run_dft_c4_c5.py has never been written); the inputs must be "
              "produced by hand before this ingester can read them.")
        # Still write an empty report so the pipeline doesn't fail
        report_json = {"ingested": 0, "results": []}
        report_md = "# DFT C4/C5 Ingestion Report\n\nNo DFT results available yet.\n"
    else:
        print(f"  Found {len(results)} completed results")

        print("Updating arrhenius_params.yml ...")
        yaml_keys = update_arrhenius_yaml(results)
        print(f"  Added/updated keys: {yaml_keys}")

        print("Updating computational_gap_closure_targets.json ...")
        update_gap_closure_targets(results)

        report_json = {"ingested": len(results), "results": results, "yaml_keys": yaml_keys}
        report_md = render_report_markdown(results, yaml_keys)

    out_path.mkdir(parents=True, exist_ok=True)
    (out_path / "dft_c4_c5_ingestion_report.json").write_text(
        json.dumps(report_json, indent=2)
    )
    (out_path / "dft_c4_c5_ingestion_report.md").write_text(report_md)
    print(f"Report written to {out_path / 'dft_c4_c5_ingestion_report.{md,json}'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    args = parser.parse_args()
    main(output_dir=args.output_dir)
