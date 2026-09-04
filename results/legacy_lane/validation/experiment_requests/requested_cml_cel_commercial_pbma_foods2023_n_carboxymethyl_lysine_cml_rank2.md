# Experiment Request `requested_cml_cel_commercial_pbma_foods2023_n_carboxymethyl_lysine_cml_rank2`

VoI rank **#2**, score **4.08**.
Generated from the current Monte-Carlo benchmark envelope.

## Why this experiment

- Compound: **Nε-(Carboxymethyl)lysine (CML)**
- Envelope miss: 3.08 dex
- 90% CI width: 0.00 dex
- ODT: unknown; decision relevance 1.00
- Rationale: measured outside 90% CI by 3.08 dex; CI width 0.00 dex; safety marker — needs SIDA-grade absolute anchor



## Benchmark context

- Benchmark target: `cml_cel_commercial_pbma_Foods2023`
- Protein type: `free`
- Conditions: {"temp_C": 150, "ph": 6.1, "water_activity": 0.45, "time_min": 20}

## Suggested protocol (`missing_absolute_anchor`)

- Method: **SIDA GC-MS/O**
- Factors:
  - Temperature (95C, 120C)
  - Time (10m, 30m)

### Instructions

Prepare standard aqueous target matrix. Add isotopically labeled internal standards (e.g. 13C-MFT, d3-methional) before extraction. Heat uniformly.

## Analytical context (mirror in intake YAML `analytical_context`)

- `headspace_method`: `HS-SPME-GC-MS`
- `quantification_mode`: `internal_standard_calibrated`
- `replicates`: `3` (minimum)
- `non_detect_policy`: `report_lod_and_do_not_backfill`
- `internal_standards`:
  - `compound_specific_internal_standard_for_n_carboxymethyl_lysine_cml`

## CRO send-to-lab checklist

- [ ] Confirm target compound identity and calibrate over a range spanning BOTH the published value (≈ 32000 ppb) and the model's 90% CI midpoint (≈ 26.5851 ppb). They differ by 1.2e+03× — that disagreement is *why* this experiment is ranked, so a calibration range covering only one of them cannot resolve it.
- [ ] Procure or confirm availability of the suggested isotopically-labeled internal standards above; substitute and **note the swap** in `analytical_context.notes`.
- [ ] Run ≥ 3 biological replicates plus a process blank and a matrix blank.
- [ ] Measure and report LoD and LoQ; do **not** backfill non-detects (`non_detect_policy: report_lod_and_do_not_backfill`).
- [ ] Quantify against the internal standards (NOT semi-quant external calibration); set `quantification_mode: internal_standard_calibrated`.
- [ ] Record the measured `conc_ppb` and `uncertainty_pct` (1σ relative) under `measured_volatiles.Nε-(Carboxymethyl)lysine (CML)` in the pre-filled intake YAML.
- [ ] Fill `provenance.source_doi` (or internal lab batch ID), `provenance.measurement_date`, and `provenance.notes` (instrument, method file).
- [ ] Set `status: measured` and return the YAML; do not edit `request_metadata` (it carries the upstream VoI rationale).

## Data return

- Pre-filled intake YAML: `data/protocols/requested_cml_cel_commercial_pbma_foods2023_n_carboxymethyl_lysine_cml_rank2.yaml`
- Fill in `measured_volatiles.Nε-(Carboxymethyl)lysine (CML).conc_ppb` (and any
  recommended_same_run_extensions) with the lab readout, set `status:
  measured`, and re-run the propagation to update the envelope.
