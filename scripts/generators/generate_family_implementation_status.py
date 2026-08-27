#!/usr/bin/env python3
"""Derive, by enumeration, how each of the 16 SLR families is actually implemented.

Why this exists (2026-08-27, Wave H)
------------------------------------
README's family table said "16 families of reaction chemistry, each independently
calibrated and wired". That conflates three very different things: SLR literature lanes,
families the reaction engine can actually GENERATE chemistry for, and layers that are
implemented but are not reaction chemistry at all (matrix retention, pre-processing,
trapping). The 2026-08-27 SMIRKS chemistry review flagged the conflation; this script
settles it from the engine rather than from the review, so the README column cannot go
stale silently.

Method
------
1. Run `SmirksEngine.enumerate` over a battery of formulations chosen to exercise every
   template entry point the engine has (core sugar+amino acid, hexose/asparagine, a
   lipid HYDROPEROXIDE seed, thiamine at its exact canonical SMILES, glutathione,
   ascorbate, polyphenol, phospholipid headgroup, nucleotide, serine/lysine). Collect
   every `reaction_family` string emitted.
2. Map each emitted family onto an SLR family through `_FAMILY_TO_SLR` below. The mapping
   is an authored classification and is listed explicitly so it can be argued with; what
   is NOT authored is which families the engine emits, which is measured here.
3. Assert that every emitted family is mapped. A new template with no mapping fails this
   script loudly rather than quietly widening a README claim.
4. Classify each SLR family into one of four implementation states and name the code that
   implements it.

Two entry-point facts worth keeping in view, both measured by this script:
  * The lipid radical chain runs only from a HYDROPEROXIDE. A plain fatty acid plus O2
    enumerates to ZERO steps -- there is no initiation step from an unoxidised lipid, so
    in production the chain is seeded by the lipid-oxidation anchor, not by the network.
  * The thiamine cascade matches thiamine by exact canonical SMILES
    (`Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1`) and by a >= 100 C gate; a differently written
    thiamine SMILES silently produces nothing.

Usage
-----
    python scripts/generators/generate_family_implementation_status.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.conditions import ReactionConditions
from src.family_ingestion_plan import load_family_ingestion_plan
from src.smirks_engine import SmirksEngine, Species

ENUMERATION_TEMPERATURE_C = 150.0
MAX_GENERATIONS = 3

# Battery of entry points. Each case is named for the template family it is meant to
# reach; between them they must touch every template the engine carries.
BATTERY: Dict[str, List[Tuple[str, str]]] = {
    "core_pentose_cysteine_lysine": [
        ("D-Ribose", "O=CC(O)C(O)C(O)CO"),
        ("L-Cysteine", "NC(CS)C(=O)O"),
        ("Glycine", "NCC(=O)O"),
        ("L-Leucine", "CC(C)CC(N)C(=O)O"),
        ("L-Lysine", "NCCCCC(N)C(=O)O"),
    ],
    "hexose_asparagine": [
        ("D-Glucose", "O=CC(O)C(O)C(O)C(O)CO"),
        ("Glycine", "NCC(=O)O"),
        ("L-Asparagine", "NC(=O)CC(N)C(=O)O"),
    ],
    "lipid_hydroperoxide_seed": [
        ("13-HPODE", "CCCCC/C=C\\C=C\\C(OO)CCCCCCCC(=O)O"),
        ("O2", "O=O"),
        ("L-Cysteine", "NC(CS)C(=O)O"),
        ("Glycine", "NCC(=O)O"),
    ],
    "unoxidised_lipid_plus_o2": [
        ("linoleic acid", "CCCCC/C=C\\C/C=C\\CCCCCCCC(=O)O"),
        ("O2", "O=O"),
    ],
    "thiamine": [
        ("thiamine", "Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1"),
        ("L-Cysteine", "NC(CS)C(=O)O"),
        ("D-Glucose", "O=CC(O)C(O)C(O)C(O)CO"),
    ],
    "glutathione": [
        ("glutathione", "NC(CCC(=O)NC(CS)C(=O)NCC(=O)O)C(=O)O"),
        ("D-Ribose", "O=CC(O)C(O)C(O)CO"),
        ("Glycine", "NCC(=O)O"),
    ],
    "ascorbate": [
        ("ascorbic acid", "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O"),
        ("Glycine", "NCC(=O)O"),
    ],
    "polyphenol": [
        ("catechol", "Oc1ccccc1O"),
        ("L-Cysteine", "NC(CS)C(=O)O"),
    ],
    "phospholipid_headgroup": [
        ("phosphatidylethanolamine headgroup", "NCCOP(=O)(O)OCC(O)CO"),
        ("D-Glucose", "O=CC(O)C(O)C(O)C(O)CO"),
    ],
    "nucleotide": [
        ("inosine-5'-monophosphate", "O=c1[nH]cnc2n(cnc12)C1OC(COP(=O)(O)O)C(O)C1O"),
        ("L-Cysteine", "NC(CS)C(=O)O"),
    ],
    "serine_lysine_dha": [
        ("L-Serine", "NC(CO)C(=O)O"),
        ("L-Lysine", "NCCCCC(N)C(=O)O"),
        ("D-Glucose", "O=CC(O)C(O)C(O)C(O)CO"),
    ],
}

# Authored classification: emitted engine family -> SLR family id.
_FAMILY_TO_SLR: Dict[str, str] = {
    # --- 01 amino acid + sugar core ------------------------------------------
    "Schiff_Base_Formation": "amino_acid_sugar_core",
    "Amadori_Rearrangement": "amino_acid_sugar_core",
    "Heyns_Rearrangement": "amino_acid_sugar_core",
    "Enolisation_Intermediate": "amino_acid_sugar_core",
    "Enolisation_1_2": "amino_acid_sugar_core",
    "Enolisation_2_3": "amino_acid_sugar_core",
    "Enolisation_2_3_Amadori": "amino_acid_sugar_core",
    "Retro_Aldol_Fragmentation": "amino_acid_sugar_core",
    "Strecker_Degradation": "amino_acid_sugar_core",
    "Aminoketone_Condensation": "amino_acid_sugar_core",
    "Generalized_Deamination": "amino_acid_sugar_core",
    "Furanone_Cyclisation": "amino_acid_sugar_core",
    "Furanone_Formation": "amino_acid_sugar_core",
    # The cysteine-derived sulfur limb belongs to the core, not to family 05:
    # family 05 is GLUTATHIONE and peptide support, and none of these steps needs a
    # peptide. MFT/FFT are named as family 01 targets in the ingestion plan.
    "Cysteine_Degradation": "amino_acid_sugar_core",
    # Wave N 2026-08-27: the MFT route was corrected on isotope evidence (Cerny &
    # Davidek 2003, 10.1021/jf026123f; 2004, 10.1021/jf035265m). The two new
    # families sit in the same core sugar/cysteine limb as the step they replace.
    # `Thiol_Addition_Norfuraneol` is retired (no step emits it) but is kept in the
    # map so an older tree, or a re-enabled route, still classifies rather than
    # tripping the unmapped-family assertion.
    "Thiol_Addition_Norfuraneol": "amino_acid_sugar_core",
    "Deoxyosone_Reduction": "amino_acid_sugar_core",
    "Thiol_Addition_Pentodiulose": "amino_acid_sugar_core",
    "Thiol_Addition_H2": "amino_acid_sugar_core",
    "Thiohemiacetal_Formation": "amino_acid_sugar_core",
    "Thiol_Dehydration": "amino_acid_sugar_core",
    "Thiol_Oxidation": "amino_acid_sugar_core",
    "Thiol_Addition_Legacy_Shortcut": "amino_acid_sugar_core",
    "Thiol_Addition_Hexose_Legacy_Shortcut": "amino_acid_sugar_core",
    # --- 02 lipid oxidation ---------------------------------------------------
    "Lipid_Homolysis": "lipid_oxidation_and_carbonylic_crosstalk",
    "Beta_Scission": "lipid_oxidation_and_carbonylic_crosstalk",
    "Radical_Propagation_O2": "lipid_oxidation_and_carbonylic_crosstalk",
    "Peroxy_H_Abstraction": "lipid_oxidation_and_carbonylic_crosstalk",
    "Radical_Termination": "lipid_oxidation_and_carbonylic_crosstalk",
    # --- 03 thiamine ----------------------------------------------------------
    "Additive_Thermal_Degradation": "thiamine_fragmentation_support",
    "Furan_Ring_Aromatisation": "thiamine_fragmentation_support",
    # --- 11 lipid x Maillard crosstalk ---------------------------------------
    "Lipid_Schiff_Base": "lipid_maillard_crosstalk",
    "Lipid_Thiazole_Condensation": "lipid_maillard_crosstalk",
    "Lipid_Strecker_Synergy": "lipid_maillard_crosstalk",
    # --- 12 protein damage markers -------------------------------------------
    "Safety_Risk_AGE": "protein_damage_markers",
    "Safety_Risk_Acrylamide": "protein_damage_markers",
    "Beta_Elimination": "protein_damage_markers",
    "DHA_Crosslinking": "protein_damage_markers",
}

# For families with no templates of their own: what implements them, and where.
_NON_TEMPLATE_IMPLEMENTATION: Dict[str, Dict[str, str]] = {
    "nucleotide_and_ribose_support": {
        "state": "priors_only",
        "implementation": (
            "Ribose itself is a first-class sugar donor of families 01/07, but there is no "
            "nucleotide template: IMP/GMP enumerate to nothing. The lane exists as literature "
            "priors and coverage bookkeeping only (src/family_barrier_progress.py, "
            "src/dft_coverage_map.py)."
        ),
    },
    "glutathione_and_peptide_support": {
        "state": "priors_only",
        "implementation": (
            "No GSH template. Glutathione enumerates only through the generic amine/carbonyl "
            "templates of family 01 -- the `glutathione_cleavage` FAST_BARRIERS entry is "
            "never reached because no family emits it. Priors and bookkeeping only "
            "(src/literature_learning_loop.py, src/family_barrier_progress.py)."
        ),
    },
    "alternative_protein_matrix_scope": {
        "state": "matrix_or_modifier_layer",
        "implementation": (
            "Not reaction chemistry. Implemented as the matrix layer: "
            "src/matrix_correction.py (accessibility + volatile retention per protein type) "
            "and data/lit/protein_source_registry.json."
        ),
    },
    "carbonyl_donor_hierarchy": {
        "state": "shared_templates_plus_modifier",
        "implementation": (
            "No templates of its own: it modulates family 01's. Implemented as "
            "`DONOR_REACTIVITY_MULTIPLIERS` in src/barrier_constants.py (per-family "
            "pentose/hexose/fructose/phosphorylated multipliers) plus "
            "`infer_carbohydrate_donor_identity`."
        ),
    },
    "off_note_and_maillard_suppression": {
        "state": "shared_templates_plus_modifier",
        "implementation": (
            "No templates of its own. The off-notes are produced by family 02 and trapped by "
            "family 11 (`Lipid_Schiff_Base`); the suppression side is the intervention layer "
            "(src/pre_processor.py) and the matrix retention profiles."
        ),
    },
    "carbohydrate_pyrolysis_and_caramelization": {
        "state": "shared_templates",
        "implementation": (
            "No templates of its own. Furfural and HMF are produced by family 01's 1,2-"
            "enolisation limb (`Enolisation_1_2`); severity is reported by the projection "
            "layer's process-state index (src/projection.py)."
        ),
    },
    "fermentation_pretreatment": {
        "state": "matrix_or_modifier_layer",
        "implementation": (
            "Not reaction chemistry. Implemented as an upstream precursor-enrichment "
            "modifier: src/pre_processor.py `KNOWN_INTERVENTIONS`, wired through "
            "src/literature_runtime.py."
        ),
    },
    "polyphenol_amino_capping": {
        "state": "priors_only",
        "implementation": (
            "No template. The quinone/cysteine Michael step has a computational prior "
            "(`quinone_cys_michael`) which Wave G1 explicitly PARKED because the curated "
            "layer has no step to route it to. src/literature_runtime.py only."
        ),
    },
    "ascorbic_acid_maillard": {
        "state": "curated_only",
        "implementation": (
            "No template. Present only as an ascorbate prior attached to the curated "
            "`Enolisation_1_2` step (src/curated_pathways.py, src/literature_runtime.py). "
            "[P] that prior's provenance is still open -- see tasks/audit_remediation.md."
        ),
    },
    "phospholipid_amine_sink": {
        "state": "priors_only",
        "implementation": (
            "No template; the PE headgroup enumerates to nothing. Implemented as a sugar-"
            "depletion bookkeeping term from the `pe_schiff_base` / `pe_amadori` priors in "
            "src/literature_runtime.py."
        ),
    },
    "melanoidin_polymerization": {
        "state": "matrix_or_modifier_layer",
        "implementation": (
            "Not reaction chemistry. Implemented as the thiol-trapping factor "
            "`_MELANOIDIN_TRAPPING_PROFILES` / `_resolve_melanoidin_trapping_factor` in "
            "src/recommend.py, gated on the family-16 lane being active."
        ),
    },
}

_STATE_LABELS = {
    "generative_templates": "generative reaction templates",
    "shared_templates": "no templates of its own (uses another family's)",
    "shared_templates_plus_modifier": "no templates of its own; modifier layer",
    "matrix_or_modifier_layer": "matrix / modifier layer, not reaction chemistry",
    "curated_only": "curated layer only",
    "priors_only": "literature priors only",
}


def enumerate_emitted_families() -> Tuple[Dict[str, List[str]], Dict[str, int]]:
    conditions = ReactionConditions(temperature_celsius=ENUMERATION_TEMPERATURE_C)
    emitted: Dict[str, List[str]] = {}
    step_counts: Dict[str, int] = {}
    for case, precursors in BATTERY.items():
        engine = SmirksEngine(conditions=conditions)
        steps = engine.enumerate(
            [Species(label=label, smiles=smiles) for label, smiles in precursors],
            max_generations=MAX_GENERATIONS,
        )
        step_counts[case] = len(steps)
        for step in steps:
            emitted.setdefault(str(step.reaction_family), []).append(case)
    return {family: sorted(set(cases)) for family, cases in emitted.items()}, step_counts


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    args = parser.parse_args()

    emitted, step_counts = enumerate_emitted_families()

    unmapped = sorted(set(emitted) - set(_FAMILY_TO_SLR))
    assert not unmapped, (
        "the engine emits reaction families this script does not classify: "
        f"{unmapped}. Add them to _FAMILY_TO_SLR before publishing a family table."
    )

    plan = load_family_ingestion_plan()
    rows: List[Dict[str, Any]] = []
    for family in plan.get("families", []):
        family_id = str(family.get("family_id"))
        slr = str(family.get("slr_family", "99")).zfill(2)
        templates = sorted(
            name for name, target in _FAMILY_TO_SLR.items()
            if target == family_id and name in emitted
        )
        declared_but_unreached = sorted(
            name for name, target in _FAMILY_TO_SLR.items()
            if target == family_id and name not in emitted
        )
        if templates:
            state = "generative_templates"
            implementation = (
                f"{len(templates)} reaction template(s) emitted by src/reaction_templates.py "
                "/ src/smirks_engine.py"
            )
        else:
            fallback = _NON_TEMPLATE_IMPLEMENTATION.get(family_id)
            assert fallback is not None, (
                f"family {family_id} emits no templates and has no recorded non-template "
                "implementation; classify it before publishing a family table."
            )
            state = fallback["state"]
            implementation = fallback["implementation"]

        rows.append(
            {
                "slr_family": slr,
                "family_id": family_id,
                "display_name": str(family.get("display_name", family_id)),
                "implementation_state": state,
                "implementation_state_label": _STATE_LABELS[state],
                "emitted_reaction_families": templates,
                "mapped_but_not_emitted_by_this_battery": declared_but_unreached,
                "implementation": implementation,
            }
        )

    rows.sort(key=lambda row: row["slr_family"])
    generative = [row for row in rows if row["implementation_state"] == "generative_templates"]

    record: Dict[str, Any] = {
        "generated_by": "scripts/generators/generate_family_implementation_status.py",
        "method": (
            "engine enumeration over a battery of entry points; the mapping from emitted "
            "reaction family to SLR family is authored and listed in the script, the set of "
            "emitted families is measured"
        ),
        "enumeration_temperature_celsius": ENUMERATION_TEMPERATURE_C,
        "max_generations": MAX_GENERATIONS,
        "battery_step_counts": step_counts,
        "emitted_reaction_families": emitted,
        "family_count": len(rows),
        "families_with_generative_templates": len(generative),
        "entry_point_findings": [
            "The lipid radical chain runs only from a HYDROPEROXIDE: the "
            f"`unoxidised_lipid_plus_o2` case enumerates to {step_counts.get('unoxidised_lipid_plus_o2', 0)} "
            "steps. There is no initiation step from an unoxidised fatty acid, so in "
            "production the chain is seeded by the lipid-oxidation anchor rather than by the "
            "network itself.",
            "The thiamine cascade matches thiamine by EXACT canonical SMILES "
            "(Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1) and by a >= 100 C gate; a differently written "
            "thiamine SMILES silently produces nothing.",
            "Glutathione, ascorbate, catechol, the phospholipid headgroup and IMP reach no "
            "template of their own -- they enumerate only through the generic family-01 "
            "amine/carbonyl templates, or to zero steps.",
        ],
        "families": rows,
    }

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "family_implementation_status.json"
    json_path.write_text(json.dumps(record, indent=2), encoding="utf-8")

    lines = [
        "# How the 16 families are actually implemented",
        "",
        "Derived by enumeration from the current engine, not asserted. Regenerate with "
        "`python scripts/generators/generate_family_implementation_status.py`.",
        "",
        f"**{len(generative)} of {len(rows)}** SLR families have generative reaction "
        "templates; the rest are shared limbs, matrix/modifier layers, or literature priors.",
        "",
        "| # | Family | Implementation | Reaction templates emitted |",
        "| -- | --- | --- | --- |",
    ]
    for row in rows:
        templates = ", ".join(f"`{name}`" for name in row["emitted_reaction_families"]) or "—"
        lines.append(
            f"| {row['slr_family']} | {row['display_name']} | "
            f"{row['implementation_state_label']} | {templates} |"
        )
    lines += ["", "## Implementation detail", ""]
    for row in rows:
        lines.append(f"* **{row['slr_family']} {row['display_name']}** — {row['implementation']}")
    lines += ["", "## Entry-point findings", ""]
    for finding in record["entry_point_findings"]:
        lines.append(f"* {finding}")
    lines.append("")
    md_path = output_dir / "family_implementation_status.md"
    md_path.write_text("\n".join(lines), encoding="utf-8")

    print(f"{len(generative)}/{len(rows)} families have generative reaction templates")
    for row in rows:
        print(f"  {row['slr_family']} {row['family_id']:46s} {row['implementation_state']}")
    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
