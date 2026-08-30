# Experiment Request `requested_thiamine_cys_glucose_120c_bolton1994_2_methyl_3_furanthiol_mft_rank1`

VoI rank **#1**, score **13.36**.
Generated from the current Monte-Carlo benchmark envelope.

## Why this experiment

- Compound: **2-Methyl-3-furanthiol (MFT)**
- Envelope miss: 2.83 dex
- 90% CI width: 0.85 dex
- ODT: 0.007 µg/kg; decision relevance 3.27
- Rationale: measured outside 90% CI by 2.83 dex; CI width 0.85 dex; ≈2e+03× ODT (decision_relevance=3.27); critical meaty odorant in a FREE-PRECURSOR system — the open question is absolute yield vs precursor dose and temperature, not matrix transfer



## Benchmark context

- Benchmark target: `thiamine_cys_glucose_120C_Bolton1994`
- Protein type: `free`
- Conditions: {"temp_C": 120.0, "ph": 5.0, "water_activity": 0.98, "time_min": 60.0}

## Suggested protocol (`free_precursor_sulfur_yield`)

- Method: **Quantitative Headspace SIDA, free-precursor aqueous model system**
- Factors:
  - Precursor Load (1x, 5x of the benchmark's stated molarities)
  - Temperature (100C, 120C, 140C)
  - pH (4.5, 5.0, 6.0)

### Instructions

Aqueous buffered model system ONLY -- do NOT add a protein isolate; this benchmark's system is free precursors in buffer and adding a matrix would answer a different question. Reproduce the benchmark's stated precursor set and molarities as the 1x arm, in sealed vials at the stated headspace ratio. Add deuterated/13C internal standards before heating, not after. Report absolute concentrations against the internal standards, plus LoD/LoQ. Where the benchmark records a value as `assumed_not_from_source` (commonly pH and the molarities), pin it explicitly in the returned YAML so the next comparison is not against an assumption.

## Analytical context (mirror in intake YAML `analytical_context`)

- `headspace_method`: `HS-SPME-GC-MS`
- `quantification_mode`: `internal_standard_calibrated`
- `replicates`: `3` (minimum)
- `non_detect_policy`: `report_lod_and_do_not_backfill`
- `internal_standards`:
  - `13C-2-methyl-3-furanthiol`

## CRO send-to-lab checklist

- [ ] Confirm target compound identity and calibrate over a range spanning BOTH the published value (≈ 13 ppb) and the model's 90% CI midpoint (≈ 0.0178032 ppb). They differ by 730× — that disagreement is *why* this experiment is ranked, so a calibration range covering only one of them cannot resolve it.
- [ ] Procure or confirm availability of the suggested isotopically-labeled internal standards above; substitute and **note the swap** in `analytical_context.notes`.
- [ ] Run ≥ 3 biological replicates plus a process blank and a matrix blank.
- [ ] Measure and report LoD and LoQ; do **not** backfill non-detects (`non_detect_policy: report_lod_and_do_not_backfill`).
- [ ] Quantify against the internal standards (NOT semi-quant external calibration); set `quantification_mode: internal_standard_calibrated`.
- [ ] Record the measured `conc_ppb` and `uncertainty_pct` (1σ relative) under `measured_volatiles.2-Methyl-3-furanthiol (MFT)` in the pre-filled intake YAML.
- [ ] Fill `provenance.source_doi` (or internal lab batch ID), `provenance.measurement_date`, and `provenance.notes` (instrument, method file).
- [ ] Set `status: measured` and return the YAML; do not edit `request_metadata` (it carries the upstream VoI rationale).

## Data return

- Pre-filled intake YAML: `data/protocols/requested_thiamine_cys_glucose_120c_bolton1994_2_methyl_3_furanthiol_mft_rank1.yaml`
- Fill in `measured_volatiles.2-Methyl-3-furanthiol (MFT).conc_ppb` (and any
  recommended_same_run_extensions) with the lab readout, set `status:
  measured`, and re-run the propagation to update the envelope.
