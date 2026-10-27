# Wave B47 ship rule — the sulfur objective, corrected against its PDFs

*Pre-registration: `results/validation/kinetic_core_b47_prereg.md`. Rule: SHIP if T2 (hold-out sulfur rows), T3 (Hofmann pH-5 bundles) and T4 (directional) hold.*

## Verdict: **DO NOT SHIP**

| test | B9 | B47 | pass |
| --- | --- | --- | --- |
| T1 cost on B47's 52 rows (recorded) | 36.46 | 22.90 | — |
| T2 hold-out sulfur rows within 3x | 1/14 (median 1.87 dex) | 2/14 (median 1.68 dex) | True |
| T3 Hofmann pH-5 rows within 3x | 0/2 (median 0.70 dex) | 1/2 (median 0.49 dex) | True |
| T4 directional, independent agree | 25/43 | 24/45 | False |
| T5 Laplace identified (recorded) | 20/23 | 17/23 | — |
| union panel within 3x (recorded) | 12/45 | 14/45 | — |

### T2 rows

| bundle | compound | B9 fold | B47 fold |
| --- | --- | ---: | ---: |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | 2-Furfurylthiol (FFT) | 4.93 | 3.12 |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | 2-Methyl-3-furanthiol (MFT) | 23.6 | 40.8 |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | 2-Furfurylthiol (FFT) | 134 | 108 |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | 2-Methyl-3-furanthiol (MFT) | 140 | 192 |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | 2-Furfurylthiol (FFT) | 1.99 | 1.57 |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | 3.82 | 2.62 |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 9.31 | 5.72 |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 477 | 406 |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 29.5 | 17.5 |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 367 | 310 |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 56 | 32.6 |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 216 | 176 |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 99.7 | 57.2 |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Furfurylthiol (FFT) | 131 | 101 |

### T3 rows

| bundle | compound | B9 fold | B47 fold |
| --- | --- | ---: | ---: |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | 2-Furfurylthiol (FFT) | 6.86 | 4.36 |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | 3.69 | 2.22 |
