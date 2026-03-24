# Family-Aware Calibration Guide

This document explains how the Maillard system handles calibration across chemistry families and how to interpret family-aware results in reports.

## Overview

The Maillard tool models Maillard chemistry through **10 numbered chemistry families** (SLR families 01–10), each with an explicit evidence posture and calibration approach. Rather than using a monolithic "one-size-fits-all" calibration, families are calibrated individually based on available evidence, with clear boundaries between benchmarked, calibration-grade, directional, and gap-only support.

## The 10 Chemistry Families

| Family ID | Name | Evidence Level | Status |
|-----------|------|-----------------|--------|
| 00 | Amino acid–sugar core | Benchmarked | ✓ Core |
| 01 | Strecker pathways | Benchmarked | ✓ Core |
| 02 | Lipid oxidation + crosstalk | Calibration-grade | ✓ Active |
| 03 | Donor hierarchy (sugar donors) | Calibration-grade | ✓ Active |
| 04 | Fermentation pretreatment | Directional | ✓ Active |
| 05 | Thiamine + nucleotide + support | Directional | ✓ Active |
| 06 | Caramelization + degradation | Directional | ⚠ Bounded |
| 07 | Sulfur chemistry (Strecker) | Benchmarked | ✓ Core |
| 08 | Off-notes + adverse markers | Calibration-grade | ✓ Active |
| 09 | Alternative matrices + mycoprotein | Directional | ⚠ Bounded |
| 10 | Safety markers + toxic thresholds | Directional | ⚠ Bounded |

**Note**: The original system (free precursors, pea matrices, soy matrices) is anchored on families 00, 01, and 07 (benchmarked), with 02, 03, 08 providing calibration-grade observable closure.

## Evidence Ladder

Each family occupies one of four rungs on the **evidence ladder**:

### 1. **Benchmarked** (Highest Confidence)
- Prediction is directly compared to external experimental measurements.
- Compounds or variables are part of the quantitative observation panel.
- Calibration is based on fit-to-data rather than prior knowledge or heuristics.
- **Families**: 00 (amino acid–sugar core), 01 (Strecker), 07 (sulfur).

**Examples**:
- MFT (methyl furan-2-carboxylate) concentration in pea/soy matrices.
- Pyrazine ratios (2-methylpyrazine / 2-ethyl-5-methylpyrazine) in controlled reactions.
- Hydrogen sulfide production from cysteine and methionine.

### 2. **Calibration-Grade** (Near-Quantitative)
- Prediction is constrained by observable anchor compounds, but not directly measured.
- Calibration uses directional literature data, literature-derived fold-changes, or matrix-specific transfer rules.
- Confidence is lower than benchmarked work but higher than pure prior.
- **Families**: 02 (lipid oxidation), 03 (donor hierarchy), 08 (off-notes).

**Examples**:
- Lipid-derived aldehydes (hexanal, octanal) retention under cooking, calibrated via headspace vapor data from lipid-only studies.
- Donor-specific baseline Maillard closure pressure, inferred from literature pH and kinetic data without direct pea/soy benchmarks.
- Adverse marker severity (furans, HMF), scaled via literature off-note sensory thresholds.

### 3. **Directional Prior** (Qualitative/Bounded)
- Prediction is supported by literature evidence but not quantitatively anchored.
- Calibration uses bounded ranges, fold-changes, or qualitative plausibility from SLR papers.
- Predictions are useful for ranking but not for absolute scale trust.
- **Families**: 04 (fermentation), 05 (thiamine/nucleotide), 06 (caramelization), 09 (alt matrices), 10 (safety).

**Examples**:
- Fermentation pretreatment effects on precursor pool (e.g., "2–5× nucleotide enrichment after 24 h fermentation" from literature ranges).
- Thiamine availability modulation, set to bounded +5% to +20% improvement based on literature reports without direct calibration.
- Mycoprotein matrix behavior, anchored to vague directional statements ("similar to soy isolate with 10% reduction in volatile yield").

### 4. **Structural Gap** (Unknown/Extrapolation)
- No direct or transferable evidence exists; prediction relies on generic chemistry or analogy.
- Calibration is either absent or purely mechanistic.
- Confidence is lowest; predictions should be treated as exploratory.
- **Part of families**: 06 (caramelization under extreme processing), 09 (coconut oil co-matrices), 10 (rare toxic threshold interactions).

**Examples**:
- Volatile loss under extreme extrusion force (>1000 bar, >180 °C); only mechanistic simulation available.
- Coconut oil lipid oxidation in plant-based matrices; no direct literature data.
- Synergistic toxic threshold interactions between multiple precursor pools.

## Family-Aware Calibration: How It Works

### Ingestion Plan
Each family is mapped to:
- **Payload role**: benchmark, calibration, prior, structural gap, or safety.
- **Target runtime modules**: which code modules apply this family.
- **Observable panel**: which compounds/markers are targets.
- **Benchmarkability status**: yes, partial, or no.
- **Next curation action**: what evidence work would increase confidence.

See `data/lit/family_ingestion_plan.json` for the machine-readable version.

### Data Sources
- **`data/lit/chemistry_family_scope_registry.json`**: Family definitions, strategic classification (core/partial/bounded/gap), and links to SLR sources.
- **`data/lit/benchmark_intake_registry.json`**: Family-tagged benchmark measurements; each entry carries `chemistry_family`, `slr_family_source`, and `payload_role`.
- **`data/lit/computational_priors.json`**: Family-specific priors for mechanistic modifiers (e.g., Maillard closure pressure deltas), calibration ranges, and uncertainty bounds.
- **`data/lit/process_gap_registry.json`**: Structural blockers tied to specific families (e.g., "no measurement of nucleotide accessibility in mycoprotein").
- **`data/lit/matrix_decision_panel.json`**: Observable targets split by family lane, labeled with evidence state and modeling regime.

### Runtime Queries
The runtime uses family-aware query helpers to:
1. Fetch all priors applicable to a compound, family, matrix type, and process state.
2. Determine observable projection factors for a given family lane.
3. Assess open gaps and missing evidence for a given recommendation.

Key functions:
- `query_family_runtime_priors(...)` – returns priors for a family + compound.
- `query_flavor_reference_entries(...)` – scores targets by evidence state.
- `query_retention_reference_entries(...)` – applies family-specific release/trapping rules.
- `build_family_upstream_contract(...)` – applies donor/pretreatment/support overlays before pathway solve.

## Reading Family-Aware Reports

### Recommendation Reports
Every run report includes:

#### Family Evidence Ladder
```json
"family_evidence_ladder": {
  "family_00_amino_acid_sugar_core": {
    "evidence_level": "benchmarked",
    "panel_targets": ["meaty_positives", "furanones"],
    "confidence_summary": "quantitative_core",
    "active_observable_count": 12,
    "open_gaps": []
  },
  "family_02_lipid_oxidation": {
    "evidence_level": "calibration_grade",
    "panel_targets": ["hexanal", "octanal"],
    "confidence_summary": "directional_with_bounds",
    "active_observable_count": 4,
    "open_gaps": ["retention_under_high_temp"]
  },
  ...
}
```

**Interpretation**:
- `benchmarked` → Use this prediction with high confidence.
- `calibration_grade` → Useful for ranking; verify scale experimentally.
- `directional` → Qualitative only; use for priority but not absolute values.
- `open_gaps` → These aspects are not yet modeled.

#### Family Runtime Support Summary
```json
"family_runtime_support_summary": {
  "active_family_count": 7,
  "lane_breakdown": {
    "core_benchmarked": ["family_00", "family_01", "family_07"],
    "active_calibration_grade": ["family_02", "family_03", "family_08"],
    "bounded_directional": ["family_04", "family_05"],
    "structural_gap": ["family_06", "family_09", "family_10"]
  },
  "recommendation_drivers": [
    {"family": "family_02", "effect": "closure_pressure_delta", "magnitude": -0.15, "confidence": "directional"},
    {"family": "family_03", "effect": "donor_limited", "magnitude": "glucose_shortage_10pct", "confidence": "calibration_grade"}
  ]
}
```

**Interpretation**:
- `recommendation_drivers` shows which families actually changed the output.
- `confidence` field indicates trust level for each driver.

#### Family-Specific Open Gaps
```json
"family_specific_open_gaps": {
  "family_02_lipid_oxidation": [
    "lipid_retention_above_170C_incomplete_data",
    "coconut_oil_oxidation_specific_rates_missing"
  ],
  "family_04_fermentation": [
    "nucleotide_enrichment_kinetics_temperature_dependent",
    "mycoprotein_fermentation_interaction_unknown"
  ]
}
```

**What to do**:
- Gaps labeled `incomplete_data` → new experiments could help.
- Gaps labeled `unknown` → structural barriers; may need computational study.
- If a gap affects your formulation decision, consider it a risk and prioritize offline refinement (DFT, MLP, etc.).

### Comparison Reports
When comparing two formulations, family-lane sensitivity shows:

```json
"family_lane_sensitivity": {
  "family_02_lipid_oxidation": {
    "toggled_off_ranking_delta": +2,
    "toggle_impact_signif": "no_change_to_top_three",
    "reason": "both_formulations_low_lipid_fraction"
  },
  "family_03_donor_hierarchy": {
    "toggled_off_ranking_delta": -8,
    "toggle_impact_signif": "changes_winner",
    "reason": "formulation_A_fructose_limited_vs_B_glucose_rich"
  }
}
```

**Interpretation**:
- `toggle_impact_signif: "changes_winner"` → this family is critical for your decision.
- `no_change` → this family does not affect your specific ranking.
- Use this to focus refinement effort on families that matter for your formulation.

## Calibration Across Matrix Families

### Pea Isolate (Matrix 01)
- **Benchmarked families**: 00, 01, 07 (amino acid–sugar, Strecker, sulfur).
- **Calibration-grade**: 02, 03, 08 (calibrated via pea-specific headspace, literature priors).
- **Bounded directional**: 04, 05 (applied with pea-specific pH correction).
- **Structural gap**: 06, 09, 10 (not calibrated; exploratory only).

### Soy Isolate (Matrix 02)
- **Benchmarked**: 00, 01, 07 (core chemistry same across proteins).
- **Calibration-grade**: 02, 03, 08 (recalibrated for soy headspace data; lipid ratios adjust for soy lecithin).
- **Bounded directional**: 04, 05 (applied with soy-specific pH and accessibility correction).
- **Structural gap**: 06, 09, 10 (not calibrated; exploratory only).

### Mycoprotein (Matrix 09)
- **Benchmarked**: 00, 01 (core chemistry assumed to transfer).
- **Calibration-grade**: 02, 07 (partially transferred; no mycoprotein-specific headspace studies).
- **Bounded directional**: 03, 04, 05 (applied with broad ±20% uncertainty).
- **Structural gap**: 06, 08, 09 (off-notes and caramelization unknown; mycoprotein scent threshold unmeasured); 10 (safety interaction unknown).

## How to Interpret "Calibration Across Families"

### Example Scenario
You are comparing Formulation A (soy isolate + ribose) vs. Formulation B (pea isolate + glucose).

**Expected Report Structure**:

1. **Core prediction (Family 00)**: Both formulations rank Formulation B higher (more sugar + pea has better Strecker yield in literature).

2. **Family 03 (Donor hierarchy)**: Formulation B is glucose-rich (high Mailllard potential). Formulation A is ribose-only (lower potential, literature-backed). →Impact: +5% increase in Formulation B ranking.

3. **Family 02 (Lipid oxidation)**: Soy lecithin in Formulation A triggers crosstalk, reducing closure. Pea protein in Formulation B has lower phospholipid load. → Impact: calibration-grade; expect −8% adjustment to Formulation A.

4. **Family 04 (Fermentation)**: Not applicable; no fermentation preprocessing. → Impact: none.

5. **Family 05 (Nucleotide)**: Soy isolate has trace nucleotides; benefit small. Pea isolate similar. → Impact: negligible, directional only.

**Final output**:
- Formulation B wins by ~12% (5% from donor + 8% from lipid avoidance, with core already favoring B).
- **Confidence**: High for families 00, 01, 07; Medium for 02, 03 (calibration-grade); Low for 04, 05 (directional).
- **Open gap**: Nucleotide + fermentation interaction in pea matrix unknown.

## Integration into Your Workflow

### For Scientists Planning Wet-Lab Experiments
1. **Run a prediction** for your candidate formulation.
2. **Check the family evidence ladder** in the report.
3. **Identify which families matter** for your decision (via family lane sensitivity).
4. **Prioritize experiments** on high-impact families with `open_gaps`.
5. **Use structural gaps** as justification for offline DFT or MLP refinement.

### For Developers Extending the System
1. **Add a new family?** Use the `family_ingestion_plan.json` template; don't hardcode family logic.
2. **Calibrate a family?** Add entries to `computational_priors.json` with family tag; do not embed priors in code.
3. **Expose a gap?** Add to `process_gap_registry.json` with family attribution.
4. **Query family logic?** Use `query_family_runtime_priors(...)` etc.; don't write new bespoke getters.

### For Report Readers (Non-Technical)
- **Green checkmark (✓)** = Benchmarked: Trust this prediction.
- **Yellow circle (⚠)** = Calibration-grade: Useful for ranking; verify experimentally.
- **Blue dash (−)** = Directional: Qualitative only; risky to use for absolute values.
- **Red X (✗)** = Structural gap: Not modeled; avoid relying on this unless offline study is available.

## See Also

- [Strategic Roadmap](../../tasks/todo.md) – S1–S8 capture family-aware implementation details.
- [Reporting Module](../src/reporting.py) – Code that generates family-aware report surfaces.
- [Generator Scripts](../../scripts/generators/README.md) – How to regenerate family artifacts.
- [Family Ingestion Plan](../../data/lit/family_ingestion_plan.json) – Machine-readable family definitions.
