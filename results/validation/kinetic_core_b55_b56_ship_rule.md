# B55 / B56 ship rules

Pre-registration: `results/validation/kinetic_core_b55_b56_prereg.md`. Four arms of the same code (n = 200, seed 0). Out-of-sample coverage on a FIXED denominator of 44 rows (a not-evaluable interval counts as a miss).

**B55: SHIP.** T1 pass (40 rows checked, 0 moved); T2 pass (median Δwidth 0.697 dex, min 0.013); T3 pass (17 → 20 of 44).

**B56: DO NOT SHIP** (arm D vs B). T4 FAIL: 26/44 = 0.591 against [0.70, 0.95]; median width 2.368 dex (reference 0.455, 1.253 dex). σ_lane: acrylamide 1.188, lipid 0.120, sulfur 0.491, trunk 0.135 dex.

## Out-of-sample coverage by arm and lane (fixed denominator)

| arm | hits / total | rate | median width | acrylamide | lipid | sulfur | trunk |
|---|---|---|---|---|---|---|---|
| A | 17/44 | 0.386 | 1.009 | 4/12 (0.91) | 7/8 (2.19) | 4/18 (1.37) | 2/6 (0.96) |
| B | 20/44 | 0.455 | 1.253 | 4/12 (0.91) | 7/8 (2.19) | 7/18 (1.74) | 2/6 (0.96) |
| C | 25/44 | 0.568 | 2.201 | 8/12 (3.98) | 7/8 (2.20) | 8/18 (2.05) | 2/6 (1.01) |
| D | 26/44 | 0.591 | 2.368 | 8/12 (3.98) | 7/8 (2.20) | 9/18 (2.37) | 2/6 (1.01) |

Artifact's own definitions (hits, total, not evaluable):

- A: honest literature [19, 48, 2] (median width 0.979), out of sample [17, 43, 1]
- B: honest literature [22, 48, 2] (median width 1.142), out of sample [20, 43, 1]
- C: honest literature [30, 50, 0] (median width 2.151), out of sample [25, 44, 0]
- D: honest literature [31, 50, 0] (median width 2.259), out of sample [26, 44, 0]

## T2 rows: sulfur away from 145 °C (A → B)

| benchmark | compound | T (°C) | width A | width B | Δ | p50 A | p50 B | measured | in A | in B |
|---|---|---|---|---|---|---|---|---|---|---|
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 100 | 1.606 | 3.044 | 1.438 | 170.799 | 60.790 | 1.280 | no | yes |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 100 | 1.541 | 2.698 | 1.157 | 8.364 | 9.189 | 6.880 | yes | yes |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 110 | 1.521 | 2.441 | 0.920 | 250.504 | 118.127 | 1.460 | no | no |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 110 | 1.140 | 2.112 | 0.972 | 21.556 | 21.888 | 3.290 | no | yes |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 120 | 1.408 | 1.916 | 0.508 | 256.477 | 162.399 | 1.680 | no | no |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 120 | 0.829 | 1.462 | 0.633 | 48.791 | 46.576 | 2.400 | no | no |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Furfurylthiol (FFT) | 130 | 1.323 | 1.468 | 0.146 | 209.981 | 163.428 | 1.620 | no | no |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 130 | 0.598 | 0.994 | 0.396 | 94.363 | 90.422 | 1.710 | no | no |
| resconi_2023_pbma_beef_identity_benchmark | furfural | 150 | 1.682 | 1.695 | 0.013 | 9.627 | 10.593 | 715.220 | no | no |
| thiamine_cys_glucose_120C_Bolton1994 | 2-Methyl-3-furanthiol (MFT) | 120 | 0.646 | 1.406 | 0.760 | 95.838 | 41.240 | 11.700 | no | yes |
