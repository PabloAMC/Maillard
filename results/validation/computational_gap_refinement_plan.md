# Computational Gap Refinement Plan

Generated: 2026-05-27T21:12:36+00:00
Targets: 8
xTB-ready targets: 8
DFT-ready targets: 7
Governance: hold_observable_first

| Priority | Target | Family | xTB | DFT | Surrogate | Ceiling | Write-back |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | hexanal_radical_quench | 11 | ready | ready_for_dft | - | ranking_only_support | data/lit/computational_priors.json, results/validation/computational_gap_dft_ingestion_report.md |
| 2 | lysinoalanine_crosslink | 12 | ready | ready_for_dft | DHA+Lys family surrogate | ranking_only_support | data/lit/computational_priors.json, results/validation/computational_gap_dft_ingestion_report.md |
| 3 | aa_ring_open_dicarbonyl | 14 | ready | ready_for_dft | Family 14 HCW surrogate | bounded_calibration | data/lit/computational_priors.json, results/validation/computational_gap_dft_ingestion_report.md |
| 4 | quinone_cys_michael | 13 | ready | ready_for_dft | - | bounded_calibration | data/lit/computational_priors.json, results/validation/computational_gap_dft_ingestion_report.md |
| 5 | pe_schiff_base | 15 | ready | ready_for_dft | - | bounded_calibration | data/lit/computational_priors.json, results/validation/computational_gap_dft_ingestion_report.md |
| 6 | pe_amadori | 15 | ready | ready_for_dft | - | bounded_calibration | data/lit/computational_priors.json, results/validation/computational_gap_dft_ingestion_report.md |
| 7 | pe_amadori_water | 15 | ready | blocked_missing_xtb_ts | - | bounded_calibration | data/lit/computational_priors.json, results/validation/computational_gap_dft_ingestion_report.md |
| 8 | asparagine_sugar_explicit_water_cluster | 13_safety | ready | ready_for_dft | - | uncertainty_narrowing_only | data/lit/computational_priors.json, results/validation/computational_gap_dft_ingestion_report.md |
