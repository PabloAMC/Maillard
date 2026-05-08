# Experiment Request `requested_thiamine_cys_ribose_100c_hofmann1996_2_methyl_3_furanthiol_mft_rank4`

VoI rank **#4**, score **5.85**.
Generated from the current Monte-Carlo benchmark envelope.

## Why this experiment

- Compound: **2-Methyl-3-furanthiol (MFT)**
- Envelope miss: inside 90% CI
- 90% CI width: 4.25 dex
- ODT: 0.0001 µg/kg; decision relevance 4.59
- Rationale: CI width 4.25 dex; ≈4e+04× ODT (decision_relevance=4.59); critical meaty odorant — multi-factor SIDA closes precursor × matrix gap
- Goal: expand model coverage
- Budget label: TBD

## Benchmark context

- Benchmark target: `thiamine_cys_ribose_100C_Hofmann1996`
- Protein type: `free`
- Conditions: {"temp_C": 100, "ph": 5.5, "water_activity": 0.98, "time_min": 30}

## Suggested protocol (`blocking_benchmark_gap`)

- Method: **Multi-factorial Quantitative Headspace SIDA**
- Factors:
  - Precursor Load (1x, 5x)
  - Temperature (90C, 130C)
  - Matrix (SPI, PPI)

### Instructions

Standard PBMA formulation baseline. Use Safe+SPME extraction to capture both highly volatile (H2S) and semi-volatile (pyrazines) simultaneously with SIDA quantitation.

## Data return

- Pre-filled intake YAML: `data/protocols/requested_thiamine_cys_ribose_100c_hofmann1996_2_methyl_3_furanthiol_mft_rank4.yaml`
- Fill in `measured_volatiles.2-Methyl-3-furanthiol (MFT).conc_ppb` (and any
  recommended_same_run_extensions) with the lab readout, set `status:
  measured`, and re-run the propagation to update the envelope.
