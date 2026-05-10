# Experiment Request `requested_cys_ribose_150c_mottram1994_2_methyl_3_furanthiol_rank5`

VoI rank **#5**, score **5.55**.
Generated from the current Monte-Carlo benchmark envelope.

## Why this experiment

- Compound: **2-methyl-3-furanthiol**
- Envelope miss: inside 90% CI
- 90% CI width: 3.70 dex
- ODT: 0.0001 µg/kg; decision relevance 5.00
- Rationale: CI width 3.70 dex; ≈2e+06× ODT (decision_relevance=5.00); critical meaty odorant — multi-factor SIDA closes precursor × matrix gap
- Goal: expand model coverage
- Budget label: TBD

## Benchmark context

- Benchmark target: `cys_ribose_150C_Mottram1994`
- Protein type: `free`
- Conditions: {"temp_C": 150, "ph": 5.5, "water_activity": 0.95, "time_min": 60}

## Suggested protocol (`blocking_benchmark_gap`)

- Method: **Multi-factorial Quantitative Headspace SIDA**
- Factors:
  - Precursor Load (1x, 5x)
  - Temperature (90C, 130C)
  - Matrix (SPI, PPI)

### Instructions

Standard PBMA formulation baseline. Use Safe+SPME extraction to capture both highly volatile (H2S) and semi-volatile (pyrazines) simultaneously with SIDA quantitation.

## Data return

- Pre-filled intake YAML: `data/protocols/requested_cys_ribose_150c_mottram1994_2_methyl_3_furanthiol_rank5.yaml`
- Fill in `measured_volatiles.2-methyl-3-furanthiol.conc_ppb` (and any
  recommended_same_run_extensions) with the lab readout, set `status:
  measured`, and re-run the propagation to update the envelope.
