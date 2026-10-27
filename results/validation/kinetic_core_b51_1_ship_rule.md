# Wave B51.1 ship rule — the 25 °C scale for the measured pH tables, and the SID across segments

*Pre-registration: `results/validation/kinetic_core_b51_1_prereg.md`. Rule: SHIP if T2, T3 and T4 hold.*

## Verdict: **DO NOT SHIP**

| test | baseline (B9, all off) | candidate (two switches on) | pass |
| --- | --- | --- | --- |
| T1 B9 cost on its 54 rows (recorded) | 18.740 | 15.470 | — |
| T2 hold-out sulfur rows within 3x | 1/14 (median 1.87 dex) | 1/14 (median 1.87 dex) | True |
| T3 Hofmann ribose + xylose pH 5 | 3/4 (median 0.37 dex) | 3/4 (median 0.37 dex) | True |
| T4 directional, independent agree | 18/31 | 17/31 | False |
| pH-axis claims (recorded) | 6/8 | 6/8 | — |
| union panel within 3x (recorded) | 18/50 (median 0.88 dex) | 18/50 (median 0.88 dex) | — |

### Rows that moved by more than 0.01 dex

| bundle | compound | fold, baseline | fold, candidate | change in abs log10 |
| --- | --- | ---: | ---: | ---: |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | 2-Furfurylthiol (FFT) | 134 | 110 | -0.086 |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | 2-Methyl-3-furanthiol (MFT) | 140 | 115 | -0.086 |
