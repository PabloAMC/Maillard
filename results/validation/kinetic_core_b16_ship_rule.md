# Wave B16 ship rule -- **DO NOT SHIP**

Prereg: `results/validation/kinetic_core_b16_prereg.md` sec. 4. Rule: SHIP if T1 (shape + folds within 0.3 dex), T2 (no B9 row +0.3 dex) and T4 (TTCA rows within 0.1 dex) hold on variant b16.

## Variant `b16` -- thiol-sink Ea 102.0 kJ/mol (ceiling 102); cost 930.98; active bounds ['k_dimer_mft', 'k_dimer_fft', 'k_pent_caramel', 'Ea_decay_thiol_sink']

* **T1 shape:** MFT at 100 C over 30/60/360/720 min = [22.5, 54.8, 109.0, 95.7] ug/L (monotone: False); FFT [2.0, 8.0, 56.3, 59.9] (monotone: True); fold rows worst 1.50 dex -> **FAIL**
  * MFT_fold_60_over_30 -0.17, MFT_fold_360_over_30 -0.79, MFT_fold_720_over_30 -0.86, FFT_fold_60_over_30 +0.43, FFT_fold_360_over_30 -0.05, FFT_fold_720_over_30 -0.02, MFT_145C20min_over_100C360min +1.50
* **T2 in-sample:** worst growth +1.90 dex on `fed_ribose_h2s_MFT`; 16 rows over 0.3, 12 over 0.5 -> **FAIL**
* **T3 Yiltirak 100/110 C levels:** median fold 198.1x (B9) -> 13.6x -> PASS
* **T4 TTCA rows:** 100C -0.68, 120C -1.25, 140C -1.25 dex -> **FAIL**
* **T6 Liu 2023 at 168 C (recorded):** MFT ['12.2', '8.48', '5.86'] ug/L at 20/40/60 min; evaluable True; RIB-T-01 decreasing True; RIB-T-02 flat False

## Variant `b16_lift` -- thiol-sink Ea 160.0 kJ/mol (ceiling 160); cost 802.28; active bounds ['k_dimer_mft', 'k_dimer_fft', 'Ea_decay_thiol_sink']

* **T1 shape:** MFT at 100 C over 30/60/360/720 min = [23.9, 75.3, 357.9, 328.9] ug/L (monotone: False); FFT [7.3, 30.2, 267.5, 337.1] (monotone: True); fold rows worst 0.89 dex -> **FAIL**
  * MFT_fold_60_over_30 -0.09, MFT_fold_360_over_30 -0.47, MFT_fold_720_over_30 -0.54, FFT_fold_60_over_30 +0.40, FFT_fold_360_over_30 -0.11, FFT_fold_720_over_30 -0.01, MFT_145C20min_over_100C360min +0.89
* **T2 in-sample:** worst growth +0.57 dex on `kang_100C_cys_conversion`; 10 rows over 0.3, 3 over 0.5 -> **FAIL**
* **T3 Yiltirak 100/110 C levels:** median fold 198.1x (B9) -> 41.6x -> PASS
* **T4 TTCA rows:** 100C -0.68, 120C -1.25, 140C -1.25 dex -> **FAIL**
* **T6 Liu 2023 at 168 C (recorded):** MFT ['5.59', '3.35', '2.21'] ug/L at 20/40/60 min; evaluable True; RIB-T-01 decreasing True; RIB-T-02 flat False

## T5 -- Laplace (variant b16)

* thiol-sink barrier identified = False, sigma = 52.77702254061807, off its bound = False; 10 of 23 identified -> **FAIL**
