# Maillard Formulation Comparison Report

**Date:** 2026-09-01 22:11:32

## 1. Metric Overview
| Metric | Ribose + leucine (no cysteine) | Ribose + cysteine + leucine |
| :--- | :---: | :---: |
| **Target Score** | 0.28 | 0.27 |
| **Off-Flavour Risk** | 0.00 | 0.01 |
| **Safety Score** (2× band risk, higher is worse) | 1.04 | 1.04 |
| **Lysine Budget** | 0.0% | 0.0% |
| **Trapping Eff.** | 0.0% | 0.0% |
| **MFT/Furfural Ratio** | 0.0000 | 0.0023 |
| **Meaty Quality Penalty** | 2.50 | 0.79 |
| **Strecker Balance** | 0.00 | 0.00 |
| **Strecker Penalty** | 0.00 | 0.00 |
| **Pyrazine Burden** | 0.00 | 0.00 |
| **Pyrazine Penalty** | 0.00 | 0.00 |
| **Furanone Penalty** | 0.35 | 0.35 |
| **Confidence** | low (47) | low (47) |
| **Prediction Mode** | directional_only | directional_only |

## 2. Competitive Highlights
- 🏆 **Highest Target Score:** Ribose + leucine (no cysteine) (0.28)
- 🛡️ **Safest Formulation:** Ribose + leucine (no cysteine) (1.04)
- 🍃 **Lowest Off-Flavour Risk:** Ribose + leucine (no cysteine) (0.00)

## Intervention Waterfall
- **baseline:** Ribose + leucine (no cysteine)
- **current:** Ribose + cysteine + leucine
- **Furans:** lowered 0.07 ppb (-41%)
- **Sulfur / Thiols:** raised 0.00 ppb
- **Aldehydes:** lowered 0.00 ppb (-41%)
- **figure:** comparison_intervention_waterfall.png

## 3. Cross-Marker Context
| Formulation | Strecker Balance | Strecker Support | Pyrazine Burden | Sulfur Trapping | Furanone Penalty | Benchmark Neighborhood | Thiamine Pathway | Thiamine Source | Expected Furanones | Missing Furanones |
| :--- | ---: | :---: | ---: | :--- | ---: | :--- | :---: | :--- | :--- | :--- |
| Ribose + leucine (no cysteine) | 0.00 | weak | 0.00 | not_applicable (1.00) | 0.35 | matrix_transfer | False | native_matrix_default_inactive | HEMF, DMHF | HEMF, DMHF |
| Ribose + cysteine + leucine | 0.00 | weak | 0.00 | moderate (0.70) | 0.35 | matrix_transfer | False | native_matrix_default_inactive | HEMF, DMHF | HEMF, DMHF |

## 4. Calibration Contrast
| Formulation | Decision Mode | Benchmark Neighborhood | Top Calibration Source | Observable Assumption | Extrapolation Axes | Missing Data Flags | Benchmark Summary |
| :--- | :--- | :--- | :--- | :--- | :--- | ---: | :--- |
| Ribose + leucine (no cysteine) | directional_hypothesis | matrix_intake_only | Pratap-Singh 2021 pea isolate ambient slurry baseline (generic furan transfer) | static_class_profile | class_level | standard_matrix_support | benchmark_neighborhood, matrix, precursors | 3 | Run uses matrix intake/headspace support and transferred accessibility priors rather than direct ranking benchmarks. |
| Ribose + cysteine + leucine | directional_hypothesis | matrix_intake_only | Pratap-Singh 2021 pea isolate ambient slurry baseline (generic furan transfer) | static_class_profile | class_level | standard_matrix_support | benchmark_neighborhood, matrix, precursors | 3 | Run uses matrix intake/headspace support and transferred accessibility priors rather than direct ranking benchmarks. |

## Trust Surface
| Formulation | Prediction Mode | Decision Mode | Top Reachability | Support Origin | Recommended Use |
| :--- | :--- | :--- | :--- | :--- | :--- |
| Ribose + leucine (no cysteine) | directional_only | directional_hypothesis | chemically_reachable | standard_matrix_support | Use directionally only; absolute concentrations should be treated as provisional. |
| Ribose + cysteine + leucine | directional_only | directional_hypothesis | chemically_reachable | standard_matrix_support | Use directionally only; absolute concentrations should be treated as provisional. |

### Ribose + leucine (no cysteine)
```text

════════════════════════════════════════════════════════════════════════════════
                         📊 MAILLARD DECISION SUMMARY
════════════════════════════════════════════════════════════════════════════════

  🟡 DIRECTIONAL ONLY
      Confidence is low. Use for hypothesis generation, not decision-grade prioritization.
  ────────────────────────────────────────────────────────────────────────────

  [0] DECISION CONFIDENCE: LOW (47/100)
      Benchmark Basis  : matrix_intake_only
      Decision Mode    : directional_hypothesis
      Prediction Mode : directional_only
      Accessibility    : free_like ✅
      Recommended Use  : Use directionally only; absolute concentrations should be treated as provisional.
      Why              : Plant-matrix support is still intake/headspace validated rather than target-ranking validated.
      Why              : Average pathway uncertainty is moderate (3.5 kcal/mol).
      Calibration      : Recommendation extrapolates beyond the strongest support on: benchmark_neighborhood, matrix, precursors.
      Extrapolation    : benchmark_neighborhood, matrix, precursors

  [1] SCIENTIFIC ENVELOPE: TRUSTED ✅

  [2] TOP-LINE METRICS:
      Target Score     : 0.28
      Off-Flavour Risk : 0.00
      Safety Score     : 1.04  (2x band-relative risk; higher is worse; >1.0 = above the action band)

════════════════════════════════════════════════════════════════════════════════

```

### Ribose + cysteine + leucine
```text

════════════════════════════════════════════════════════════════════════════════
                         📊 MAILLARD DECISION SUMMARY
════════════════════════════════════════════════════════════════════════════════

  🟡 DIRECTIONAL ONLY
      Confidence is low. Use for hypothesis generation, not decision-grade prioritization.
  ────────────────────────────────────────────────────────────────────────────

  [0] DECISION CONFIDENCE: LOW (47/100)
      Benchmark Basis  : matrix_intake_only
      Decision Mode    : directional_hypothesis
      Prediction Mode : directional_only
      Accessibility    : free_like ✅
      Recommended Use  : Use directionally only; absolute concentrations should be treated as provisional.
      Why              : Plant-matrix support is still intake/headspace validated rather than target-ranking validated.
      Why              : Average pathway uncertainty is moderate (3.5 kcal/mol).
      Calibration      : Recommendation extrapolates beyond the strongest support on: benchmark_neighborhood, matrix, precursors.
      Extrapolation    : benchmark_neighborhood, matrix, precursors

  [1] SCIENTIFIC ENVELOPE: TRUSTED ✅

  [2] TOP-LINE METRICS:
      Target Score     : 0.27
      Off-Flavour Risk : 0.01
      Safety Score     : 1.04  (2x band-relative risk; higher is worse; >1.0 = above the action band)

════════════════════════════════════════════════════════════════════════════════

```

## 5. Provenance
- **artifact_kind:** formulation_comparison
- **generated_at:** 2026-09-01T22:11:13.014601
- **generator:** scripts/generators/generate_report_visual_examples.py --output-dir results/validation/report_visual_examples --docs-asset-dir docs/assets
- **repository:** workspace | branch cleaning | commit 2208fd0 | dirty True
- **input_fingerprint_sha256:** b1272b62fda83200e95ca6d9e4be97ab29cf3d78f4c9fa684309c1f09aa1bad7
- **scientific_surface:**
  - scientific_reference: docs/reference/SCIENTIFIC_REFERENCE.md
  - benchmark_summary: results/validation/benchmark_summary.md
  - validated_envelope: results/validation/validated_envelope.md
  - validation_overview: results/validation/validation_overview.md
  - matrix_decision_panel: data/lit/matrix_decision_panel.json
  - chemistry_family_scope_registry: data/lit/chemistry_family_scope_registry.json
  - family_ingestion_plan_registry: data/lit/family_ingestion_plan.json
  - family_identifier_contract: results/validation/family_identifier_contract.md
  - family_identifier_contract_json: results/validation/family_identifier_contract.json
  - family_strategy_policy: results/validation/family_strategy_policy.md
  - family_strategy_policy_json: results/validation/family_strategy_policy.json
  - family_payload_coverage: results/validation/family_payload_coverage.md
  - family_payload_coverage_json: results/validation/family_payload_coverage.json
  - matrix_family_coverage_registry: data/lit/matrix_family_coverage_registry.json
  - benchmark_intake_registry: data/lit/benchmark_intake_registry.json
  - computational_priors: data/lit/computational_priors.json
  - slr_incorporation_matrix: results/literature/slr_incorporation_matrix.json
  - flavor_reference_payloads: data/lit/flavor_reference_payloads.json
  - process_state_calibrations: data/lit/process_state_calibrations.json
  - retention_reference_payloads: data/lit/retention_reference_payloads.json
  - process_gap_registry: data/lit/process_gap_registry.json
  - safety_reference_payloads: data/lit/safety_reference_payloads.json
  - primary_benchmark_protocol: docs/protocols/PPI_SPI_PRIMARY_BENCHMARK_PROTOCOL.md
  - primary_benchmark_contract: data/protocols/ppi_spi_primary_benchmark_contract.json
  - literature_learning_loop: results/validation/literature_learning_loop.md
  - literature_learning_loop_json: results/validation/literature_learning_loop.json
  - literature_backlog: results/validation/literature_backlog.md
  - literature_backlog_json: results/validation/literature_backlog.json
  - deep_research_runtime_queue: results/validation/deep_research_runtime_queue.md
  - deep_research_runtime_queue_json: results/validation/deep_research_runtime_queue.json
  - literature_runtime_templates: results/validation/literature_runtime_templates.json
  - family_ingestion_plan: results/validation/family_ingestion_plan.md
  - family_ingestion_plan_json: results/validation/family_ingestion_plan.json
  - matrix_target_status: results/validation/matrix_target_status.md
  - matrix_target_status_json: results/validation/matrix_target_status.json
  - chemistry_family_scope: results/validation/chemistry_family_scope.md
  - chemistry_family_scope_json: results/validation/chemistry_family_scope.json
  - matrix_family_coverage: results/validation/matrix_family_coverage.md
  - matrix_family_coverage_json: results/validation/matrix_family_coverage.json
  - family_sensitivity: results/validation/family_sensitivity.md
  - family_sensitivity_json: results/validation/family_sensitivity.json
  - family_lane_validation: results/validation/family_lane_validation.md
  - family_lane_validation_json: results/validation/family_lane_validation.json
  - refinement_surrogate_patches: data/lit/refinement_surrogate_patches.json
- **campaign_name:** unknown

## 6. Glossary
Plain-language meaning of the labels used above. The model is honest about *how* it knows what it claims; this section names that vocabulary.

**Three different tier vocabularies appear in this report. They are not interchangeable.** Each answers a different question:

**1. `tier` — how well benchmark-supported is *this run*?** Emitted on the run-level *Confidence & Support* block, on every *Compound Confidence* row, and on every *Aggregate Sensory Confidence* row. It is a band on a 0-100 confidence score, and it always travels with a `prediction_mode`:

| `tier` | score band | paired `prediction_mode` | what it licenses |
| :--- | :--- | :--- | :--- |
| `high` | >= 85 | `benchmark_supported_quantitative` | quantitative prioritisation before wet-lab confirmation |
| `medium` | 65-84 | `ranking_supported` | ranking and triage; verify absolute levels experimentally |
| `low` | 45-64 | `directional_only` | direction only; absolute ppb is provisional |
| `exploratory` | < 45 | `hypothesis_only` | hypothesis generation, not decision-grade |

**2. `calibration_evidence_strength` (shown as *Evidence* in the Calibration Summary) — what kind of anchor stands behind a compound's projection?** Strongest to weakest, and only these values are emitted at compound level:
- `literature_anchored` — a published measurement on a directly comparable system backs this compound's retention/partition treatment.
- `conditional_literature_anchored` — literature-anchored, but only under stated conditions (pH / process-state caveats attached).
- `class_anchored` — anchored at compound-class level (e.g. "sulfur volatiles"), not for this molecule specifically.
- `directional_transferred` — a prior transferred from an adjacent matrix or process state; direction is meaningful, magnitude is not.
- `process_state_mismatch` — an anchor exists, but for a different process state than the one you asked about; the nearest state was substituted.
- `heuristic` — no anchor at all; a built-in class default. Ranking use only.

  When the run is out of calibration scope every one of these is demoted one notch (see the scope banner below), and the pre-demotion value is preserved as `scope_demoted_from` in `report.json`.

**3. `confidence_tier` — how strong is the *literature prior* behind a chemistry lane?** This one comes from the curated literature registries (`data/lit/`), not from your run, and uses a five-point scale: `high`, `medium_high`, `medium`, `medium_low`, `low`. It grades the source, not the prediction: a `high` `confidence_tier` prior can still feed an `exploratory` `tier` prediction, because your formulation may sit far from where that prior was measured.

  *Name collision, stated plainly:* the campaign/comparison JSON also carries a key spelled `confidence_tier` that holds the run-level `tier` value (vocabulary 1), kept as a legacy alias. Prefer the `tier` key alongside it; only `confidence_tier` inside `data/lit/` payloads means vocabulary 3.

**Scope banner.** A `⚠️ Out of calibration scope` banner at the top of the report means the matrix or process you asked about lies outside the convex hull of formulations the model has been calibrated against. The predictions are still emitted, but every compound's evidence strength is demoted one notch.

**Reachability** (`reachability_status`). `chemically_reachable` — the compound is produced by an enumerated, barrier-scored pathway from your precursors. `conditionally_reachable` — reachable only under an assumption the run had to make. `merely_plausible` — no enumerated route; the number is a class-level projection.

**Observable assumption** (`observable_assumption_summary`). A pipe-joined triple: retention runtime mode, calibration fallback mode, support origin — e.g. `static_class_profile | class_level | standard_matrix_support` means the volatile's matrix retention came from a static class profile, its calibration fell back to class level, and no special matrix-support route was used.

**Confidence envelope.** `0.038 ppb [0.012-0.089, 90% CI]` means the p50 (median) Monte-Carlo prediction is `0.038 ppb`, with the central 90 % of samples between `0.012` and `0.089 ppb`. A compound printed without an interval had no envelope sampled. Wide envelopes make coverage cheap — read coverage and width together.

**Intervention waterfall.** When two formulations are compared, the per-compound delta is decomposed into class-level (e.g. "thiols", "aldehydes") and per-precursor (e.g. "cysteine", "glutathione") contributions. Per-precursor attribution sums to the compound delta and is explicit about attribution mode.

Full machine-readable trust artifacts: `results/validation/`. Per-compound 90 % envelope: `results/validation/prediction_uncertainty.md`. External hold-out: `results/validation/external_validation_report.md`.

## 7. Recommended next experiment
Ranked by value-of-information (envelope miss × CI width × ODT-anchored decision relevance). Source: `results/validation/experiment_value_ranking.json`.

Scoped to matrix `pea_iso` (filtered from the global ranking).

| Rank | VoI | Benchmark | Matrix | Compound | DoE template | Why this one |
| ---: | ---: | --- | --- | --- | --- | --- |
| 21 | 0.97 | `pea_isolate_ribose_cysteine_100C_45min_Internal2026` | `pea_iso` | 2,5-dimethylpyrazine | `missing_kinetic_dataset` | CI width 6.47 dex; ≈2e-08× ODT (decision_relevance=0.50); wide envelope — time-course narrows the rate-limiting step |
| 22 | 0.97 | `pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026` | `pea_iso` | 2,5-dimethylpyrazine | `missing_kinetic_dataset` | CI width 6.47 dex; ≈2e-08× ODT (decision_relevance=0.50); wide envelope — time-course narrows the rate-limiting step |
| 25 | 0.68 | `pea_isolate_ribose_cysteine_100C_45min_Internal2026` | `pea_iso` | bis(2-methyl-3-furyl) disulfide | `missing_kinetic_dataset` | CI width 4.56 dex; ≈0.06× ODT (decision_relevance=0.50); wide envelope — time-course narrows the rate-limiting step |

_How to use this: run `./scripts/docker_maillard.sh next-experiment --top 3` to materialise pre-filled intake YAMLs and protocol Markdown for each row. Ingest the resulting measurement via `./scripts/docker_maillard.sh ingest --file results.csv ...`._

### Ribose + leucine (no cysteine) — sensory readout
Per-compound odour activity value (OAV = predicted ppb ÷ curated odour threshold). An axis is *above threshold* when at least one of its compounds reaches OAV ≥ 1. Compounds without a curated odour threshold are listed but excluded from axis roll-ups. Source thresholds: `data/species/desirable_targets.yml`, `data/species/off_flavour_targets.yml`.

**Headline:** meaty no ODT; off-notes no ODT; safety no ODT.

### Axis roll-up
| Axis | Compounds (with ODT) | Above threshold | Max OAV | Top contributor |
| --- | ---: | ---: | ---: | --- |
| meaty | 0 | 0 | n/a | _no compound with curated ODT in this run_ |
| off-note | 0 | 0 | n/a | _no compound with curated ODT in this run_ |
| safety | 0 | 0 | n/a | _no compound with curated ODT in this run_ |

### Per-compound OAV (90 % CI)
| Compound | Axis | ODT (μg/kg) | Predicted ppb (p50) | OAV (p50) | OAV p5 | OAV p95 | ≥1? |
| --- | --- | ---: | ---: | ---: | ---: | ---: | :---: |
| 3-methylbutanal | unclassified | 1.5 | 0.00152 | 1.02e-03 | n/a | n/a | · |
| 3-Methylbutanal | unclassified | 1.5 | 0.00152 | 1.02e-03 | n/a | n/a | · |
| furfural | unclassified | 3000 | 0.181 | 6.05e-05 | 5.47e-07 | 6.68e-03 | · |
| Furfural | unclassified | 3000 | 0.181 | 6.05e-05 | 5.47e-07 | 6.68e-03 | · |
| 2,5-dimethylpyrazine | unclassified | 1800 | 7.25e-06 | 4.03e-09 | 1.11e-12 | 3.25e-06 | · |
| 2,5-Dimethylpyrazine | unclassified | 1800 | 7.25e-06 | 4.03e-09 | 1.11e-12 | 3.25e-06 | · |
| O=Cc1ccco1 | unclassified | n/a | 0.181 | n/a | n/a | n/a | — |
| CC(C)CC=O | unclassified | n/a | 0.00152 | n/a | n/a | n/a | — |
| Cc1cnc(C)cn1 | unclassified | n/a | 7.25e-06 | n/a | n/a | n/a | — |

_3/9 predicted compounds have no curated odour threshold; they appear in the per-compound table but do not contribute to axis roll-ups._

### Ribose + cysteine + leucine — sensory readout
Per-compound odour activity value (OAV = predicted ppb ÷ curated odour threshold). An axis is *above threshold* when at least one of its compounds reaches OAV ≥ 1. Compounds without a curated odour threshold are listed but excluded from axis roll-ups. Source thresholds: `data/species/desirable_targets.yml`, `data/species/off_flavour_targets.yml`.

**Headline:** meaty below threshold; off-notes no ODT; safety no ODT.

### Axis roll-up
| Axis | Compounds (with ODT) | Above threshold | Max OAV | Top contributor |
| --- | ---: | ---: | ---: | --- |
| meaty | 4 | 0 | 0.210 | 2-furfurylthiol |
| off-note | 0 | 0 | n/a | _no compound with curated ODT in this run_ |
| safety | 0 | 0 | n/a | _no compound with curated ODT in this run_ |

### Per-compound OAV (90 % CI)
| Compound | Axis | ODT (μg/kg) | Predicted ppb (p50) | OAV (p50) | OAV p5 | OAV p95 | ≥1? |
| --- | --- | ---: | ---: | ---: | ---: | ---: | :---: |
| 2-furfurylthiol | meaty | 0.01 | 0.0021 | 0.210 | 4.60e-04 | 23.2 | · |
| 2-Furfurylthiol (FFT) | meaty | 0.01 | 0.0021 | 0.210 | 4.60e-04 | 23.2 | · |
| 2-methyl-3-furanthiol | meaty | 0.007 | 0.000252 | 0.036 | 3.25e-04 | 3.97 | · |
| 2-Methyl-3-furanthiol (MFT) | meaty | 0.007 | 0.000252 | 0.036 | 3.25e-04 | 3.97 | · |
| bis(2-methyl-3-furyl) disulfide | unclassified | 0.02 | 5.51e-05 | 2.76e-03 | 2.63e-06 | 0.305 | · |
| Bis(2-methyl-3-furyl) disulfide | unclassified | 0.02 | 5.51e-05 | 2.76e-03 | 2.63e-06 | 0.305 | · |
| 3-methylbutanal | unclassified | 1.5 | 0.000903 | 6.02e-04 | n/a | n/a | · |
| 3-Methylbutanal | unclassified | 1.5 | 0.000903 | 6.02e-04 | n/a | n/a | · |
| furfural | unclassified | 3000 | 0.108 | 3.58e-05 | 3.24e-07 | 3.96e-03 | · |
| Furfural | unclassified | 3000 | 0.108 | 3.58e-05 | 3.24e-07 | 3.96e-03 | · |
| 2,5-dimethylpyrazine | unclassified | 1800 | 4.3e-06 | 2.39e-09 | 6.61e-13 | 1.93e-06 | · |
| 2,5-Dimethylpyrazine | unclassified | 1800 | 4.3e-06 | 2.39e-09 | 6.61e-13 | 1.93e-06 | · |
| O=Cc1ccco1 | unclassified | n/a | 0.108 | n/a | n/a | n/a | — |
| CC(C)CC=O | unclassified | n/a | 0.000903 | n/a | n/a | n/a | — |
| Cc1cnc(C)cn1 | unclassified | n/a | 4.3e-06 | n/a | n/a | n/a | — |
| SCc1ccco1 | unclassified | n/a | 0.0021 | n/a | n/a | n/a | — |
| Cc1occc1S | unclassified | n/a | 0.000252 | n/a | n/a | n/a | — |
| Cc1occc1SSc1ccoc1C | unclassified | n/a | 5.51e-05 | n/a | n/a | n/a | — |

_6/18 predicted compounds have no curated odour threshold; they appear in the per-compound table but do not contribute to axis roll-ups._

