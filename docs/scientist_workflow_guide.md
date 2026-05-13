# Scientist Workflow Guide: Flavour Design with Maillard

This guide is the shortest scientist-facing path from formulation idea to report bundle, then from fresh GC-MS data back into the model review loop.

## 1. Refresh The Trust Surface

Start by regenerating the repository-level trust artifacts in Docker:

```bash
./scripts/docker_maillard.sh summary
./scripts/docker_maillard.sh validated-envelope
./scripts/docker_maillard.sh validation-figures
```

Read these first:

- `results/validation/benchmark_summary.md`
- `results/validation/prediction_uncertainty.md`
- `results/validation/external_validation_report.md`

Goal: confirm the current trust headline before you act on any formulation ranking.

## 2. Generate A Scientist-Facing Report

Fastest path: the bundled quickstart demo runs a representative pea-isolate
baseline and a cysteine-enrichment comparison in one command, with no flags
to assemble:

```bash
./scripts/docker_maillard.sh quickstart
```

That writes `results/quickstart/baseline/` and `results/quickstart/comparison/`
and prints where to look. For a custom formulation, run the normal report path
through the Docker wrapper:

```bash
./scripts/docker_maillard.sh run "python scripts/run_pipeline.py \
	--sugars ribose,glucose \
	--amino-acids cysteine,leucine \
	--ratios ribose:0.5,glucose:0.2,cysteine:0.2,leucine:0.1 \
	--ph 5.5 \
	--temp 105 \
	--time-minutes 45 \
	--protein-type pea_iso \
	--target meaty \
	--minimize beany \
	--report \
	--output-dir results/scientist_run"
```

Review these outputs in `results/scientist_run/`:

- `report.md`: narrative summary for scientists and reviewers, with a plain-language glossary at the bottom
- `report.json`: machine-readable payload for downstream analysis
- `compound_confidence_overlay.png`: per-compound p5/p50/p95 whiskers, colored by evidence class

![Compound confidence overlay](assets/report_compound_confidence_overlay.png)

Goal: see both the predicted ppb and the confidence structure that supports each compound.

## 3. Compare Two Candidates And Read The Intervention Waterfall

For a pairwise comparison, run the campaign path so the comparison bundle emits the cross-formulation waterfall:

```bash
./scripts/docker_maillard.sh campaign \
	--names "Cysteine Enrichment (Basic),Premium Meaty Mix" \
	--ph 5.5 \
	--temp 105 \
	--target-tag meaty \
	--minimize-tag beany \
	--campaign-name "Scientist head-to-head" \
	--output-dir results/scientist_compare
```

Review these outputs in `results/scientist_compare/`:

- `comparison.md`: side-by-side report for both candidates
- `comparison.json`: structured comparison payload
- `comparison_intervention_waterfall.png`: class-level delta summary between the baseline and the current candidate

![Intervention waterfall](assets/report_comparison_intervention_waterfall.png)

Goal: answer questions like “did glutathione raise sulfur notes while lowering aldehydes?” without reading raw tables.

## 4. Inspect Trust, Safety, And Process Limits

Use these repository-level artifacts when you need to understand why a result is or is not decision-ready:

```bash
./scripts/docker_maillard.sh explain-formulation MY_FORM_NAME
./scripts/docker_maillard.sh literature-learning-loop
```

Read:

- `results/validation/what_can_i_trust_today.md`
- `results/validation/literature_learning_loop.md`
- `results/validation/calibration_refit_decision.md`

The current 2026-05-10 recalibration decision is explicit: the matrix observable refit was rejected because the panel median residual and the 39/48 trust headline did not improve.

## 5. Ingest New GC-MS Results Before You Commit Them

### 5a. The 30-minute first ingest (template path)

If you have never ingested data before, start from a bundled CSV template
instead of writing the schema from scratch. End-to-end this path takes
≤ 30 minutes and covers ~80 % of plant-protein GC-MS panels.

1. **Copy a template.** Pick the shape that matches your assay:

   - HS-SPME-GC-MS quantitation → [data/ingest_templates/hs_spme_gc_ms_template.csv](../data/ingest_templates/hs_spme_gc_ms_template.csv)
   - Stable-isotope-dilution (SIDA) → [data/ingest_templates/sida_template.csv](../data/ingest_templates/sida_template.csv)

   ```bash
   cp data/ingest_templates/hs_spme_gc_ms_template.csv data/intake/my_run.csv
   ```

2. **Replace the example rows with your numbers.** Keep the canonical
   column names (`compound`, `observed_ppb`, `rsd_pct`, `temperature_c`,
   `time_min`, `ph`, `water_activity`, `source_reference`) — the parser
   also accepts loose aliases like `Analyte Name`, `Aw`, `RSD %`. Mandatory
   fields and aliases are documented in
   [data/ingest_templates/README.md](../data/ingest_templates/README.md).

3. **Run a dry-run preview** (no canonical write):

   ```bash
   ./scripts/docker_maillard.sh ingest \
       --file data/intake/my_run.csv \
       --protein-type pea_iso \
       --process-state heated_matrix \
       --precursor "Pea Protein Isolate=1000.0" \
       --precursor "D-Ribose=1.0" \
       --precursor "L-Cysteine=1.0" \
       --source-reference "My 2026 PPI off-note panel" \
       --output-dir results/validation/ingest_preview
   ```

   Note: matrix conditions (`temperature_c`, `time_min`, `ph`,
   `water_activity`) come from the CSV rows; you do not need to repeat them
   on the CLI when using the templates.

4. **Inspect the preview artifacts** in `results/validation/ingest_preview/`:

   - `*_ingest_preview.md` — does the trust-loop delta and benchmark
     alignment match what you expected?
   - `*_support_delta.md` — does the support-strengthening verdict make
     sense for the matrix you measured?

5. **Confirm** when the preview looks right by re-running the same command
   with `--confirm`. That writes the normalized intake YAML next to the
   preview so the next `summary` / `validation-figures` run picks it up.

> Out-of-scope warning: if your matrix is wheat gluten or mycoprotein, or
> outside the calibrated process envelope, the report will carry a
> `⚠️ Out of calibration scope` banner. That is intentional — the ingest
> still works; the predictions just inherit a one-notch evidence-strength
> demotion until a matrix-specific calibration anchor exists.

### 5b. Direct CLI ingest (advanced)

```bash
./scripts/docker_maillard.sh ingest \
	--file results.csv \
	--protein-type soy_iso \
	--process-state extrusion_structured \
	--precursor "Soy Protein Isolate=1000" \
	--precursor "D-Ribose=1.0" \
	--precursor "L-Cysteine=1.0" \
	--temp-c 140 \
	--ph 6.0 \
	--water-activity 0.45 \
	--time-min 1.0 \
	--source-reference "SPI extrusion GC-MS panel" \
	--output-dir results/validation/ingest_preview
```

This command does **not** write a canonical intake payload until you add `--confirm`.

Preview outputs:

- `*_ingest_preview.md` and `*_ingest_preview.json`: benchmark-added preview, trust-loop delta, top envelope tightening, and top median shifts
- `*_support_delta.md` and `*_support_delta.json`: support-strengthening and promotion-readiness comparison against the closest aligned benchmark

Add `--confirm` when the preview looks correct. That writes the normalized intake YAML alongside the preview artifacts.

## 6. What To Share Externally

If you need to hand the result to collaborators, pair the run-level bundle with these trust artifacts:

- `results/validation/prediction_uncertainty.md`
- `results/validation/external_validation_report.md`
- `results/validation/calibration_refit_decision.md`

That keeps the formulation claim, the current external hold-out state, and the latest accept/reject calibration decision in one packet.
