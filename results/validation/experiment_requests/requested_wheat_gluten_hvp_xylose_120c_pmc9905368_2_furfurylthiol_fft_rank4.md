# Experiment Request `requested_wheat_gluten_hvp_xylose_120c_pmc9905368_2_furfurylthiol_fft_rank4`

VoI rank **#4**, score **5.06**.
Generated from the current Monte-Carlo benchmark envelope.

## Why this experiment

- Compound: **2-Furfurylthiol (FFT)**
- Envelope miss: 0.83 dex
- 90% CI width: 3.33 dex
- ODT: 0.01 µg/kg; decision relevance 1.79
- Rationale: measured outside 90% CI by 0.83 dex; CI width 3.33 dex; ≈6e+01× ODT (decision_relevance=1.79); critical meaty odorant — multi-factor SIDA closes precursor × matrix gap



## Benchmark context

- Benchmark target: `wheat_gluten_hvp_xylose_120C_PMC9905368`
- Protein type: `free`
- Conditions: {"temp_C": 120, "ph": 6.0, "water_activity": 0.92, "time_min": 30}

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

- [ ] Confirm target compound identity and expected dynamic range (model 90% CI midpoint ≈ 0.00521544 ppb).
- [ ] Procure or confirm availability of the suggested isotopically-labeled internal standards above; substitute and **note the swap** in `analytical_context.notes`.
- [ ] Run ≥ 3 biological replicates plus a process blank and a matrix blank.
- [ ] Measure and report LoD and LoQ; do **not** backfill non-detects (`non_detect_policy: report_lod_and_do_not_backfill`).
- [ ] Quantify against the internal standards (NOT semi-quant external calibration); set `quantification_mode: internal_standard_calibrated`.
- [ ] Record the measured `conc_ppb` and `uncertainty_pct` (1σ relative) under `measured_volatiles.2-Furfurylthiol (FFT)` in the pre-filled intake YAML.
- [ ] Fill `provenance.source_doi` (or internal lab batch ID), `provenance.measurement_date`, and `provenance.notes` (instrument, method file).
- [ ] Set `status: measured` and return the YAML; do not edit `request_metadata` (it carries the upstream VoI rationale).

## Data return

- Pre-filled intake YAML: `data/protocols/requested_wheat_gluten_hvp_xylose_120c_pmc9905368_2_furfurylthiol_fft_rank4.yaml`
- Fill in `measured_volatiles.2-Furfurylthiol (FFT).conc_ppb` (and any
  recommended_same_run_extensions) with the lab readout, set `status:
  measured`, and re-run the propagation to update the envelope.
