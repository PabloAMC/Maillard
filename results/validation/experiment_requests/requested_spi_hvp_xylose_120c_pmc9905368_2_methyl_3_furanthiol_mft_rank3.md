# Experiment Request `requested_spi_hvp_xylose_120c_pmc9905368_2_methyl_3_furanthiol_mft_rank3`

VoI rank **#3**, score **6.87**.
Generated from the current Monte-Carlo benchmark envelope.

## Why this experiment

- Compound: **2-Methyl-3-furanthiol (MFT)**
- Envelope miss: inside 90% CI
- 90% CI width: 6.32 dex
- ODT: 0.0001 µg/kg; decision relevance 3.63
- Rationale: CI width 6.32 dex; ≈4e+03× ODT (decision_relevance=3.63); critical meaty odorant — multi-factor SIDA closes precursor × matrix gap
- Goal: expand model coverage
- Budget label: TBD

## Benchmark context

- Benchmark target: `spi_hvp_xylose_120C_PMC9905368`
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

## Data return

- Pre-filled intake YAML: `data/protocols/requested_spi_hvp_xylose_120c_pmc9905368_2_methyl_3_furanthiol_mft_rank3.yaml`
- Fill in `measured_volatiles.2-Methyl-3-furanthiol (MFT).conc_ppb` (and any
  recommended_same_run_extensions) with the lab readout, set `status:
  measured`, and re-run the propagation to update the envelope.
