# Experiment Request `requested_cys_ribose_140c_hofmann1998_2_furfurylthiol_rank3`

VoI rank **#3**, score **8.34**.
Generated from the current Monte-Carlo benchmark envelope.

## Why this experiment

- Compound: **2-furfurylthiol**
- Envelope miss: 0.14 dex
- 90% CI width: 2.68 dex
- ODT: 0.01 µg/kg; decision relevance 4.30
- Rationale: measured outside 90% CI by 0.14 dex; CI width 2.68 dex; ≈2e+04× ODT (decision_relevance=4.30); critical meaty odorant — multi-factor SIDA closes precursor × matrix gap



## Benchmark context

- Benchmark target: `cys_ribose_140C_Hofmann1998`
- Protein type: `free`
- Conditions: {"temp_C": 140, "ph": 5.0, "water_activity": 0.98, "time_min": 30}

## Suggested protocol (`blocking_benchmark_gap`)

- Method: **Multi-factorial Quantitative Headspace SIDA**
- Factors:
  - Precursor Load (1x, 5x)
  - Temperature (90C, 130C)
  - Matrix (SPI, PPI)

### Instructions

Standard PBMA formulation baseline. Use Safe+SPME extraction to capture both highly volatile (H2S) and semi-volatile (pyrazines) simultaneously with SIDA quantitation.

## Analytical context (mirror in intake YAML `analytical_context`)

- `headspace_method`: `HS-SPME-GC-MS`
- `quantification_mode`: `internal_standard_calibrated`
- `replicates`: `3` (minimum)
- `non_detect_policy`: `report_lod_and_do_not_backfill`
- `internal_standards`:
  - `13C-2-furfurylthiol`
  - `hexanal-d12`

## CRO send-to-lab checklist

- [ ] Confirm target compound identity and expected dynamic range (model 90% CI midpoint ≈ 40.8504 ppb).
- [ ] Procure or confirm availability of the suggested isotopically-labeled internal standards above; substitute and **note the swap** in `analytical_context.notes`.
- [ ] Run ≥ 3 biological replicates plus a process blank and a matrix blank.
- [ ] Measure and report LoD and LoQ; do **not** backfill non-detects (`non_detect_policy: report_lod_and_do_not_backfill`).
- [ ] Quantify against the internal standards (NOT semi-quant external calibration); set `quantification_mode: internal_standard_calibrated`.
- [ ] Record the measured `conc_ppb` and `uncertainty_pct` (1σ relative) under `measured_volatiles.2-furfurylthiol` in the pre-filled intake YAML.
- [ ] Fill `provenance.source_doi` (or internal lab batch ID), `provenance.measurement_date`, and `provenance.notes` (instrument, method file).
- [ ] Set `status: measured` and return the YAML; do not edit `request_metadata` (it carries the upstream VoI rationale).

## Data return

- Pre-filled intake YAML: `data/protocols/requested_cys_ribose_140c_hofmann1998_2_furfurylthiol_rank3.yaml`
- Fill in `measured_volatiles.2-furfurylthiol.conc_ppb` (and any
  recommended_same_run_extensions) with the lab readout, set `status:
  measured`, and re-run the propagation to update the envelope.
