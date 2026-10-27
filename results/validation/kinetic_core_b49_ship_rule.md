# Wave B49 ship rule — norfuraneol's other sulfide products

*Pre-registration: `results/validation/kinetic_core_b49_prereg.md`. Rule: SHIP if T1 (pH 6.5 mechanism, held out), T2 (hold-out sulfur rows), T3 (Hofmann pH-5 bundles) and T4 (directional) hold.*

## Verdict: **DO NOT SHIP**

| test | B9 | B49 | pass |
| --- | --- | --- | --- |
| T1 NF + Cys, 140 C: MFT fall pH 4.5 -> 6.5 (measured >= 150x; MFT at 6.5 < 0.001 mol %) | 0.482x (0.16 mol %) | 1.01x (0.0743 mol %) | False |
| T2 hold-out sulfur rows within 3x | 1/14 (median 1.87 dex) | 3/14 (median 1.68 dex) | True |
| T3 Hofmann pH-5 rows within 3x | 0/2 (median 0.70 dex) | 1/2 (median 0.36 dex) | True |
| T4 directional, independent agree | 25/43 | 25/44 | True |
| T5 Laplace identified (recorded) | 20/23 | 17/24 | — |
| union panel within 3x (recorded) | 12/45 | 15/45 | — |

### T2 rows

| bundle | compound | B9 fold | B49 fold |
| --- | --- | ---: | ---: |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | 2-Furfurylthiol (FFT) | 4.93 | 2.35 |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | 2-Methyl-3-furanthiol (MFT) | 23.6 | 55.4 |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | 2-Furfurylthiol (FFT) | 134 | 232 |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | 2-Methyl-3-furanthiol (MFT) | 140 | 549 |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | 2-Furfurylthiol (FFT) | 1.99 | 1.3 |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | 3.82 | 2.04 |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 9.31 | 4.34 |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 477 | 386 |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 29.5 | 13 |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 367 | 300 |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 56 | 24 |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 216 | 168 |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 99.7 | 42 |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Furfurylthiol (FFT) | 131 | 92.1 |

### T3 rows

| bundle | compound | B9 fold | B49 fold |
| --- | --- | ---: | ---: |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | 2-Furfurylthiol (FFT) | 6.86 | 3.17 |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | 3.69 | 1.65 |
