# Wave B51 ship rule — the pH state, corrected against its thermodynamic sources

*Pre-registration: `results/validation/kinetic_core_b51_prereg.md`. Rule: SHIP (A) if (A) passes T2, T3, T4; else SHIP (B) if (B) passes them; else DO NOT SHIP.*

## Verdict: **DO NOT SHIP**

| test | baseline (B9, switches off) | (A) B51 refit | pass | (B) B9 vector, switches on | pass |
| --- | --- | --- | --- | --- | --- |
| T2 hold-out sulfur rows within 3x | 1/14 (median 1.87 dex) | 1/14 (median 1.61 dex) | True | 0/14 (median 1.42 dex) | True |
| T3 Hofmann ribose + xylose pH 5 | 3/4 (median 0.37 dex) | 0/4 (median 1.80 dex) | False | 0/4 (median 2.48 dex) | False |
| T4 directional, independent agree | 18/31 | 18/32 | True | 17/31 | False |
| T5 pH-axis claims (recorded) | 6/8 | 5/8 | — | 5/8 | — |
| union panel within 3x (recorded) | 18/50 (median 0.88 dex) | 16/50 (median 0.95 dex) | — | 15/50 (median 0.95 dex) | — |

**T1 (recorded), cost on B47's 52 rows:** baseline 36.46, (A) 20.41, (B) 43.29. Laplace for (A): {'identified': 15, 'free': 23}.

**Per-switch attribution** — B9's frozen vector on its own 54 rows, one switch on at a time:

| switch | cost |
| --- | ---: |
| none | 18.740 |
| heat_capacity | 25.723 |
| ionic_strength | 19.635 |
| cysteine_gkl | 26.094 |
| h2s_suleimenov | 18.882 |
| label_scale_tables | 15.470 |
| sid_across_segments | 18.740 |
| all | 32.463 |

### Rows that moved by more than 0.05 dex

| arm | bundle | compound | fold, baseline | fold, arm | change in abs log10 |
| --- | --- | --- | ---: | ---: | ---: |
| A_refit | hofmann1998_c2c3_recombination_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | 1.49 | 1.83 | +0.090 |
| A_refit | hofmann1998_c2c3_recombination_145C_20min_pH7 | 2-Methyl-3-furanthiol (MFT) | 4.65 | 5.41 | +0.065 |
| A_refit | hofmann1998_furan2aldehyde_h2s_145C_20min_pH5 | 2-Furfurylthiol (FFT) | 1.29 | 1.62 | +0.100 |
| A_refit | hofmann1998_norfuraneol_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | 2.18 | 1.46 | -0.174 |
| A_refit | hofmann1998_norfuraneol_h2s_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | 1.38 | 2.06 | +0.173 |
| A_refit | hofmann1998_ribose_cysteine_145C_20min_pH5 | 2-Furfurylthiol (FFT) | 1.58 | 47.2 | +1.475 |
| A_refit | hofmann1998_ribose_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | 2.76 | 118 | +1.632 |
| A_refit | resconi_2023_pbma_beef_identity_benchmark | furfural | 93.3 | 12.9 | -0.858 |
| A_refit | mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | 2-Furfurylthiol (FFT) | 4.93 | 11.1 | +0.351 |
| A_refit | mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | 2-Methyl-3-furanthiol (MFT) | 23.6 | 13.4 | -0.246 |
| A_refit | mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | 2-Furfurylthiol (FFT) | 134 | 65.6 | -0.309 |
| A_refit | mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | 2-Methyl-3-furanthiol (MFT) | 140 | 4.72e+03 | +1.527 |
| A_refit | mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | 2-Furfurylthiol (FFT) | 1.99 | 37.4 | +1.274 |
| A_refit | mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | 3.82 | 85.4 | +1.349 |
| A_refit | mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 9.31 | 1.42 | -0.816 |
| A_refit | mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 477 | 139 | -0.535 |
| A_refit | mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 29.5 | 3.43 | -0.934 |
| A_refit | mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 367 | 107 | -0.534 |
| A_refit | mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 56 | 5.11 | -1.040 |
| A_refit | mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 216 | 68.4 | -0.500 |
| A_refit | mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 99.7 | 7.35 | -1.133 |
| A_refit | mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Furfurylthiol (FFT) | 131 | 44.1 | -0.472 |
| B_frozen_b9 | hofmann1998_c2c3_recombination_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | 1.49 | 1.8 | +0.082 |
| B_frozen_b9 | hofmann1998_norfuraneol_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | 2.18 | 1.91 | -0.058 |
| B_frozen_b9 | hofmann1998_ribose_cysteine_145C_20min_pH5 | 2-Furfurylthiol (FFT) | 1.58 | 655 | +2.618 |
| B_frozen_b9 | hofmann1998_ribose_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | 2.76 | 175 | +1.801 |
| B_frozen_b9 | resconi_2023_pbma_beef_identity_benchmark | furfural | 93.3 | 148 | +0.201 |
| B_frozen_b9 | mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | 2-Methyl-3-furanthiol (MFT) | 23.6 | 5.73 | -0.614 |
| B_frozen_b9 | mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | 2-Furfurylthiol (FFT) | 134 | 159 | +0.076 |
| B_frozen_b9 | mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | 2-Methyl-3-furanthiol (MFT) | 140 | 674 | +0.682 |
| B_frozen_b9 | mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | 2-Furfurylthiol (FFT) | 1.99 | 519 | +2.417 |
| B_frozen_b9 | mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | 3.82 | 126 | +1.518 |
| B_frozen_b9 | mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 9.31 | 3 | -0.492 |
| B_frozen_b9 | mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 477 | 151 | -0.499 |
| B_frozen_b9 | mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 29.5 | 7.93 | -0.570 |
| B_frozen_b9 | mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 367 | 86.6 | -0.627 |
| B_frozen_b9 | mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 56 | 12.4 | -0.655 |
| B_frozen_b9 | mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 216 | 37.5 | -0.761 |
| B_frozen_b9 | mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 99.7 | 18.2 | -0.739 |
| B_frozen_b9 | mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Furfurylthiol (FFT) | 131 | 16.6 | -0.895 |
