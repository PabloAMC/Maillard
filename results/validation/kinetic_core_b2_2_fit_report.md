# Kinetic core, Build Wave B2.2 -- separated decay barriers + a pH-trajectory state

Modules 1 (sulfur formation) and 2 (thiol consumption) of `docs/reference/FIT_HOLDOUT_DECLARATION.md`.

- network: **47 species, 79 reactions** (15 inherited from B1's trunk, 64 new), carbon, nitrogen AND sulfur balance enforced at import
- objective: **58 declared FIT rows**, **48 free parameters** (FOUR more than B2.1: two decay barriers and two pH-drift constants), final cost **8.991**, reduced chi-square **1.80**
- fitted pH parameters: **0**. The pH shape is structural (acid/base catalysis + measured H2S speciation + Zheng & Ho's measured four-pH thermolysis set).
- branch-fraction constants: **0**. Every split is a ratio of time-integrated mass-action fluxes.
- activation energies in the NAMED CONSUMPTION channels: **0**, by policy (inventory sec. C.1 / B.7) -- unchanged.
- activation energies on the DECAY LUMPS: **two**, one per named family, both fitted. `k_dimer_decay` and `k_h2s_loss` keep B2.1's lumped formation Ea because nothing measures either at two temperatures.

## The two decay barriers

| family | keys | fitted Ea (kJ/mol) | identified by |
|---|---|---:|---|
| `thiol_sink` | `k_mft_decay`, `k_fft_decay`, `k_thiol_decay`, `k_thiolate_loss` | 248.0 | Kumazawa 2003 Fig. 3 (2-furfurylthiol survival at 121 C, formation absent, six declared-FIT pH levels) + Hofma |
| `carbonyl_sink` | `k_fur_decay`, `k_osone_decay`, `k_nf_decay` | 150.2 | Kang 2026 SI T-S4 furfural at the declared 100 and 120 C rungs (the dossier's own best-behaved Arrhenius serie |

The lumped FORMATION barrier is **76.2 kJ/mol** for comparison. A decay barrier BELOW it means the sinks lose ground to formation as temperature rises (which is the direction B2.1's own 140 C diagnosis asked for); ABOVE it means they gain. Neither bound was tightened toward either answer and both barriers were started from independent uniform draws over 20-250 kJ/mol.

## The pH-trajectory state

Two calibrated constants: acid yield per sink event **0.968**, Amadori secondary-ammonium pKa **10.77**. Both are calibrated on Amendment 7's THREE declared drift anchors and on nothing else; every other acid-base constant in `ph_state.py` is textbook and immovable.

| system | buffer | label pH | in-situ start | in-situ end | cooled end | acid, mM | SID, mM | beta0, mM/pH |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| `hofmann_pentose_pH5` | phosphate | 5.00 | 4.94 | 3.83 | 3.59 | 86.7 | 502.4 | 15.7 |
| `hofmann_glucose_pH5` | phosphate | 5.00 | 4.94 | 4.21 | 3.96 | 0.0 | 502.4 | 15.7 |
| `hofmann_fructose_pH5` | phosphate | 5.00 | 4.94 | 3.80 | 3.56 | 0.0 | 502.4 | 15.7 |
| `fed_ribose_h2s` | phosphate | 5.00 | 4.96 | 4.94 | 4.97 | 0.3 | 502.4 | 14.9 |
| `fed_tdp_h2s` | phosphate | 5.00 | 4.96 | 4.34 | 4.09 | 18.4 | 502.4 | 14.9 |
| `fed_furfural_h2s` | phosphate | 5.00 | 4.96 | 4.96 | 5.00 | 0.0 | 502.4 | 14.9 |
| `fed_nf_h2s` | phosphate | 5.00 | 4.96 | 4.96 | 5.00 | 0.0 | 502.4 | 14.9 |
| `fed_nf_cys` | phosphate | 5.00 | 4.95 | 4.95 | 5.00 | 0.0 | 502.4 | 15.4 |
| `fed_c2c3` | phosphate | 5.00 | 4.96 | 4.96 | 5.00 | 0.0 | 502.4 | 14.9 |
| `fed_c2c3_pH3` | phosphate | 3.00 | 3.39 | 3.39 | 3.00 | 0.0 | 437.4 | 223.8 |
| `fed_c2c3_pH7` | phosphate | 7.00 | 6.82 | 6.82 | 7.00 | 0.0 | 694.0 | 287.8 |
| `fed_thiamine` | phosphate | 5.00 | 4.96 | 4.96 | 5.00 | 0.0 | 502.4 | 14.9 |
| `fed_mgo_h2s_1to1` | phosphate | 5.00 | 4.96 | 4.96 | 5.00 | 0.0 | 502.4 | 14.9 |
| `fed_mgo_h2s_1to2` | phosphate | 5.00 | 4.96 | 4.96 | 5.00 | 0.0 | 502.4 | 14.9 |
| `zhou_arp_cys_pH7` | none | 7.00 | 5.90 | 3.50 | 3.50 | 7.7 | 0.9 | 12.8 |
| `zhou_arp_pH7` | none | 7.00 | 5.95 | 3.21 | 3.21 | 5.5 | 0.0 | 0.6 |
| `whitfield_nf_cys` | phosphate ⚠assumed | 4.50 | 4.68 | 4.68 | 4.50 | 0.0 | 498.7 | 15.9 |
| `whitfield_nf_h2s` | phosphate ⚠assumed | 4.50 | 4.68 | 4.68 | 4.50 | 0.0 | 498.8 | 15.7 |
| `vanseeventer_130C` | phosphate ⚠assumed | 5.00 | 4.95 | 4.13 | 3.90 | 31.9 | 502.4 | 14.8 |
| `zhang_fig1_cys` | phosphate ⚠assumed | 4.90 | 4.86 | 4.17 | 3.98 | 24.8 | 501.5 | 13.8 |
| `zhang_fig1_gcys` | phosphate ⚠assumed | 4.90 | 4.86 | 4.17 | 3.98 | 24.8 | 501.5 | 13.8 |
| `cerny_ternary` | phosphate ⚠assumed | 5.00 | 4.91 | 3.50 | 3.28 | 269.2 | 502.4 | 17.2 |
| `kang_ttca_100` | none | 7.00 | 5.59 | 9.69 | 11.42 | 0.8 | 8.5 | 0.2 |
| `kang_ttca_120` | none | 7.00 | 5.31 | 9.40 | 11.42 | 5.1 | 8.5 | 0.1 |
| `kang_freecys_100` | none | 7.00 | 6.08 | 6.16 | 7.08 | 0.0 | 0.4 | 4.9 |
| `kang_freecys_120` | none | 7.00 | 5.90 | 6.12 | 7.22 | 0.0 | 0.4 | 6.1 |
| `kumazawa_pH3_0` | citrate_phosphate | 3.00 | 3.08 | 3.08 | 3.00 | 0.0 | 218.2 | 147.2 |
| `kumazawa_pH4_0` | citrate_phosphate | 4.00 | 3.95 | 3.95 | 4.00 | 0.0 | 300.3 | 66.1 |
| `kumazawa_pH5_0` | citrate_phosphate | 5.00 | 4.93 | 4.93 | 5.00 | 0.0 | 367.0 | 63.6 |
| `kumazawa_pH5_4` | citrate_phosphate | 5.40 | 5.36 | 5.36 | 5.40 | 0.0 | 393.2 | 56.2 |
| `kumazawa_pH6_0` | citrate_phosphate | 6.00 | 6.01 | 6.01 | 6.00 | 0.0 | 435.5 | 85.3 |
| `kumazawa_pH6_4` | citrate_phosphate | 6.40 | 6.39 | 6.39 | 6.40 | 0.0 | 476.0 | 129.8 |
| `zhou_drift_pH6` | none | 6.00 | 5.23 | 3.29 | 3.29 | 7.5 | 0.1 | 2.6 |
| `zhou_drift_pH7` | none | 7.00 | 5.90 | 3.50 | 3.50 | 7.7 | 0.9 | 12.8 |
| `zhou_drift_pH8` | none | 8.00 | 6.86 | 5.07 | 5.07 | 8.3 | 6.5 | 20.2 |

## Read this before reading anything else

**3 of 48 parameters are individually identified.** With 58 rows against 48 free parameters the fit has 10 degrees of freedom, so the row-by-row agreement below is NOT evidence that the model is right. What the corpus determines is a set of RATIOS and SHAPES -- branch shares, MFT/FFT, dimer/monomer, the pH ordering -- and those are what the hold-out report scores. No individual rate constant in this table should be quoted as a measured quantity, cited elsewhere, or carried into another module as if it were one.

## Row-by-row

| row | observed | predicted | fold | source |
|---|---:|---:|---:|---|
| `hofmann_ribose_FFT` | 0.00106 | 0.0009612 | 0.91x | Hofmann 1998 T1, ribose pH 5.0: 12.1 ug/100 mL = 121 ppb |
| `hofmann_ribose_MFT` | 0.001734 | 0.001493 | 0.86x | Hofmann 1998 T1, ribose pH 5.0: 19.8 ug/100 mL = 198 ppb |
| `hofmann_xylose_FFT` | 0.0008409 | 0.0009612 | 1.14x | Hofmann 1998 T1, xylose pH 5.0: 9.6 ug/100 mL = 96 ppb |
| `hofmann_xylose_MFT` | 0.001253 | 0.001493 | 1.19x | Hofmann 1998 T1, xylose pH 5.0: 14.3 ug/100 mL = 143 ppb |
| `hofmann_glucose_FFT` | 0.0002453 | 0.0002044 | 0.83x | Hofmann 1998 T1, glucose pH 5.0: 2.8 ug/100 mL = 28 ppb |
| `hofmann_glucose_MFT` | 0.0001664 | 0.0002828 | 1.70x | Hofmann 1998 T1, glucose pH 5.0: 1.9 ug/100 mL = 19 ppb |
| `hofmann_fructose_FFT` | 0.0002803 | 0.000336 | 1.20x | Hofmann 1998 T1, fructose pH 5.0: 3.2 ug/100 mL = 32 ppb |
| `hofmann_fructose_MFT` | 0.000219 | 0.0001086 | 0.50x | Hofmann 1998 T1, fructose pH 5.0: 2.5 ug/100 mL = 25 ppb |
| `hofmann_ribose_NF_insitu` | 4.779 | 4.825 | 1.01x | Hofmann 1998 T5, ribose pH 5.0: 54 530 ug/100 mL = 4.78 mmol/L |
| `hofmann_ribose_FUR_insitu` | 0.007025 | 0.007121 | 1.01x | Hofmann 1998 T5, ribose pH 5.0: 67.5 ug/100 mL furan-2-aldehyde |
| `fed_ribose_h2s_FFT` | 0.008 | 0.008638 | 1.08x | Hofmann 1998 T3: ribose + H2S -> FFT, 9.2 ug, 0.008 mol% |
| `fed_ribose_h2s_MFT` | 0.01 | 0.01001 | 1.00x | Hofmann 1998 T6: ribose + H2S -> MFT, 15.1 ug, 0.01 mol% |
| `fed_tdp_h2s_FFT` | 0.08 | 0.1183 | 1.48x | Hofmann 1998 T3: 3-deoxyribosulose + H2S -> FFT, 78.6 ug, 0.08 mol% |
| `fed_furfural_h2s_FFT` | 0.48 | 0.4197 | 0.87x | Hofmann 1998 T3: furan-2-aldehyde + H2S -> FFT, 550.8 ug, 0.48 mol% |
| `fed_nf_h2s_MFT` | 0.19 | 0.291 | 1.53x | Hofmann 1998 T4: norfuraneol + H2S -> MFT, 211.2 ug, 0.19 mol% |
| `fed_nf_cys_MFT` | 0.05 | 0.09133 | 1.83x | Hofmann 1998 T4: norfuraneol + cysteine -> MFT, 50.8 ug, 0.05 mol% |
| `fed_c2c3_MFT` | 0.24 | 0.1325 | 0.55x | Hofmann 1998 T10: C2 + C3 -> MFT, 268.1 ug, 0.24 mol% |
| `fed_c2c3_MFT_pH3` | 0.01 | 0.1334 | 13.34x | Hofmann 1998 T10 pH ladder: 15.5 ug at pH 3.0 |
| `fed_c2c3_MFT_pH7` | 0.27 | 0.1014 | 0.38x | Hofmann 1998 T10 pH ladder: 311.5 ug at pH 7.0 |
| `fed_thiamine_MFT` | 0.01 | 0.006426 | 0.64x | Hofmann 1998 T8: thiamin -> MFT, 8.2 ug, 0.01 mol% |
| `fed_mgo_h2s_1to1_MP` | 1.8 | 2.148 | 1.19x | Hofmann 1998 T7: 2-oxopropanal + H2S 1:1 -> 1650 ug, 1.8 mol% |
| `fed_mgo_h2s_1to2_MP` | 4 | 4.239 | 1.06x | Hofmann 1998 T7: 2-oxopropanal + H2S 1:2 -> 3600 ug, 4.0 mol% |
| `zhou_pH7_MFT` | 0.01392 | 0.008565 | 0.62x | Zhou 2023 T1 pH 7: MFT 1588.57 +/- 21.24 ug/L |
| `zhou_pH7_FFT` | 0.006639 | 0.003779 | 0.57x | Zhou 2023 T1 pH 7: FFT 757.965 +/- 13.03 ug/L |
| `zhou_pH7_MFT_over_FFT` | 2.096 | 2.267 | 1.08x | Zhou 2023 T1 pH 7, MFT/FFT = 2.096 [D] |
| `zhou_pH7_dimer_over_MFT` | 0.0323 | 0.008133 | 0.25x | Zhou 2023 T1 pH 7: 102.59 ug/L dimer against 1588.57 MFT (molar ratio  |
| `zhou_pH7_ACTZ` | 9.201e-05 | 9.185e-05 | 1.00x | Zhou 2023 T1 pH 7: 2-acetylthiazole 11.70 +/- 2.14 ug/L |
| `zhou_pH7_FUR_arp_alone` | 0.01394 | 0.03072 | 2.20x | Zhou 2023 T1 pH 7, ARP ALONE: 2-furfural 1339.37 +/- 83.04 ug/L |
| `whitfield_nf_cys_MFT` | 0.15 | 0.1037 | 0.69x | Whitfield & Mottram 1999 T1: NF + cysteine, pH 4.5, 140 C, 0.150 mol%  |
| `whitfield_nf_h2s_MFT` | 0.12 | 0.05621 | 0.47x | Whitfield & Mottram 1999 T1: NF + H2S 1:2, pH 4.5, 0.120 mol% |
| `whitfield_mercaptoketone_over_MFT` | 16.3 | 14.33 | 0.88x | Whitfield 1999: mercaptoketones : MFT = 16.3 : 1 from fed NF |
| `vanseeventer_cys_conversion` | 0.55 | 0.9554 | 1.74x | van Seeventer 2001: cysteine 33 -> 15 mM (55%), 130 C / 20 min |
| `vanseeventer_ribose_conversion` | 0.75 | 0.942 | 1.26x | van Seeventer 2001: ribose 100 -> 25 mM (75%), 130 C / 20 min |
| `zhang_fig1_cys_dimer_over_MFT` | 0.0429 | 0.05601 | 1.31x | Zhang 2024 Fig. 1: 0.115 ng/mL dimer against 1.34 MFT = 8.6% BY MASS = |
| `zhang_fig1_gcys_dimer_over_MFT` | 0.2711 | 1.691 | 6.24x | Zhang 2024 Fig. 1: 1.09 against 2.01 = 54.2% by mass = 0.2711 molar |
| `zhang_fig1_cys_MMFT_over_MFT` | 0.0213 | 0.02137 | 1.00x | Zhang 2024 Fig. 1: MMFT 0.04 against MFT 1.34 ng/mL = 0.0213 molar |
| `cerny_ternary_thiamine_share` | 0.54 | 0.4524 | 0.84x | Cerny 2007b Table 3, full ternary: 54% unlabelled (thiamine) : 46% 13C |
| `cerny_isomer_split` | 1 | 1.003 | 1.00x | Cerny 2007 T4: 2-mercapto-3-pentanone is 94->95% XYLOSE-derived while  |
| `kang_100C_MFT` | 1.084e-05 | 7.585e-06 | 0.70x | Kang 2026 SI T-S4, 100 C: MFT 1.237 ug/L (Tier A, R^2 0.9989) |
| `kang_120C_MFT` | 1.216e-05 | 1.474e-05 | 1.21x | Kang 2026 SI T-S4, 120 C: MFT 1.388 ug/L; the same value appears in th |
| `kang_100C_FFT` | 3.271e-05 | 2.122e-05 | 0.65x | Kang 2026 SI T-S4, 100 C: FFT 3.734 ug/L (Tier A, R^2 0.9992) |
| `kang_120C_FFT` | 3.598e-05 | 4.037e-05 | 1.12x | Kang 2026 SI T-S4, 120 C: FFT 4.107 ug/L |
| `kang_100C_FUR` | 3.519e-05 | 5.544e-05 | 1.58x | Kang 2026 SI T-S4, 100 C: furfural 3.381 ug/L |
| `kang_120C_FUR` | 6.029e-05 | 4.98e-05 | 0.83x | Kang 2026 SI T-S4, 120 C: furfural 5.793 ug/L |
| `kang_100C_MFT_over_FFT` | 0.3311 | 0.3575 | 1.08x | Kang 2026 SI sec. 7b: FFT/MFT = 3.02 at 100 C |
| `kang_100C_cys_conversion` | 0.162 | 0.1618 | 1.00x | Kang 2026 SI Fig. S4 (digitised): free-cysteine conversion 16.2% at 10 |
| `kang_120C_cys_conversion` | 0.387 | 0.3846 | 0.99x | Kang 2026 SI Fig. S4 (digitised): 38.7% at 120 C / 120 min |
| `kang_free_cys_yield_ceiling` | 0.163 | 0.0516 | 0.32x | Kang 2026 SI Fig. S3 (digitised, kang2026_SI_extraction.md sec. 6c): p |
| `kumazawa_FFT_retention_pH3_0` | 0.995 | 0.9641 | 0.97x | Kumazawa 2003 Fig. 3: 99.5% of 1 ppm 2-furfurylthiol survives 121 C /  |
| `kumazawa_FFT_retention_pH4_0` | 0.962 | 0.9484 | 0.99x | Kumazawa 2003 Fig. 3: 96.2% of 1 ppm 2-furfurylthiol survives 121 C /  |
| `kumazawa_FFT_retention_pH5_0` | 0.891 | 0.8165 | 0.92x | Kumazawa 2003 Fig. 3: 89.1% of 1 ppm 2-furfurylthiol survives 121 C /  |
| `kumazawa_FFT_retention_pH5_4` | 0.795 | 0.6457 | 0.81x | Kumazawa 2003 Fig. 3: 79.5% of 1 ppm 2-furfurylthiol survives 121 C /  |
| `kumazawa_FFT_retention_pH6_0` | 0.451 | 0.3113 | 0.69x | Kumazawa 2003 Fig. 3: 45.1% of 1 ppm 2-furfurylthiol survives 121 C /  |
| `kumazawa_FFT_retention_pH6_4` | 0.11 | 0.1699 | 1.54x | Kumazawa 2003 Fig. 3: 11.0% of 1 ppm 2-furfurylthiol survives 121 C /  |
| `zhou_final_pH_from_pH6` | 3.22 | 3.291 | 1.02x | Zhou 2023 Fig. 2 [D, +/-0.06]: ARP + Cys, initial pH 6.0 -> final pH 3 |
| `zhou_final_pH_from_pH7` | 3.42 | 3.496 | 1.02x | Zhou 2023 Fig. 2 [D, +/-0.06]: ARP + Cys, initial pH 7.0 -> final pH 3 |
| `zhou_final_pH_from_pH8` | 5.07 | 5.066 | 1.00x | Zhou 2023 Fig. 2 [D, +/-0.06]: ARP + Cys, initial pH 8.0 -> final pH 5 |
| `yaghmur_fft_share_ceiling` | 0.012 | 0.008754 | 0.73x | Yaghmur et al. 2005 sec. 3.1 p. 226 (<1% FFT conversion of the furfura |

## Parameters

| parameter | log10 k(145 C) | 95% half-width | identified? |
|---|---:|---:|---|
| `k_pent_dpo` | -1.339 | -- | **UNIDENTIFIED** |
| `k_pent_tdp` | -1.547 | -- | **UNIDENTIFIED** |
| `k_dpo_c2c3` | -2.754 | -- | **UNIDENTIFIED** |
| `k_arp_dpo` | -6.069 | -- | **UNIDENTIFIED** |
| `k_arp_tdp` | -1.863 | -- | **UNIDENTIFIED** |
| `k_dpo_nf` | 0.184 | -- | **UNIDENTIFIED** |
| `k_dpo_ptr` | -4.294 | -- | **UNIDENTIFIED** |
| `k_dpo_ddp` | 0.107 | -- | **UNIDENTIFIED** |
| `k_tdp_fur` | -2.453 | -- | **UNIDENTIFIED** |
| `k_ddp_mft` | -5.832 | -- | **UNIDENTIFIED** |
| `k_fur_fft` | -4.257 | -- | **UNIDENTIFIED** |
| `k_nf_mft` | -2.255 | -- | **UNIDENTIFIED** |
| `k_nf_mp3p` | -2.616 | -- | **UNIDENTIFIED** |
| `k_mgo_mp` | -3.912 | -- | **UNIDENTIFIED** |
| `k_ha_mp_mft` | -4.084 | -- | **UNIDENTIFIED** |
| `k_glc_ha` | -0.531 | -- | **UNIDENTIFIED** |
| `k_thi_hmp` | -2.676 | -- | **UNIDENTIFIED** |
| `k_thi_mesh` | -3.248 | -- | **UNIDENTIFIED** |
| `k_hmp_mft` | -2.712 | -- | **UNIDENTIFIED** |
| `k_hmp_mp2p` | -3.440 | -- | **UNIDENTIFIED** |
| `k_cys_actz` | -3.167 | -- | **UNIDENTIFIED** |
| `k_dimer_mft` | 0.220 | -- | **UNIDENTIFIED** |
| `k_dimer_fft` | 0.097 | 1.735 | yes |
| `k_mmft` | -2.332 | -- | **UNIDENTIFIED** |
| `k_mft_decay` | 0.058 | -- | **UNIDENTIFIED** |
| `k_fft_decay` | -0.584 | 3.107 | yes |
| `k_dimer_decay` | -8.086 | -- | **UNIDENTIFIED** |
| `k_nf_decay` | -4.560 | -- | **UNIDENTIFIED** |
| `k_fur_decay` | 0.414 | -- | **UNIDENTIFIED** |
| `k_h2s_loss` | -1.318 | -- | **UNIDENTIFIED** |
| `k_osone_decay` | -0.858 | -- | **UNIDENTIFIED** |
| `k_thiol_decay` | -3.045 | -- | **UNIDENTIFIED** |
| `k_pent_caramel` | -2.922 | -- | **UNIDENTIFIED** |
| `k_pent_thermal` | -3.515 | -- | **UNIDENTIFIED** |
| `k_glc_fur` | -4.283 | -- | **UNIDENTIFIED** |
| `k_arp_tdp_th` | -3.675 | -- | **UNIDENTIFIED** |
| `k_arp_dpo_th` | -0.991 | -- | **UNIDENTIFIED** |
| `k_ddp_mft_hs` | -4.918 | -- | **UNIDENTIFIED** |
| `k_fur_fft_hs` | 0.499 | -- | **UNIDENTIFIED** |
| `k_ttca_cys` | -1.160 | -- | **UNIDENTIFIED** |
| `k_ttca_deg` | -1.987 | -- | **UNIDENTIFIED** |
| `k_cys_thermal` | -2.024 | 0.714 | yes |
| `k_thiolate_loss` | -2.056 | -- | **UNIDENTIFIED** |
| `Ea_lumped_formation` | 76.159 | -- | **UNIDENTIFIED** |
| `Ea_decay_thiol_sink` | 247.975 | -- | **UNIDENTIFIED** |
| `Ea_decay_carbonyl_sink` | 150.221 | -- | **UNIDENTIFIED** |
| `ph_acid_yield_per_sink_event` | 0.968 | -- | **UNIDENTIFIED** |
| `ph_arp_secondary_ammonium_pKa` | 10.767 | -- | **UNIDENTIFIED** |

## Branch shares at the fitted point (computed, never stored)

**hofmann_pentose_pH5**: MFT_share_thiamine_route = 0, MFT_share_sugar_routes = 1, MFT_share_intact_skeleton_route = 0.001029, MFT_share_norfuraneol_route = 0.999, MFT_share_C2_plus_C3_route = 7.063e-11, FFT_share_of_furfural_flux = 0.008754, MFT_consumed_share_dimerisation = 0.0003584, MFT_consumed_share_MMFT = 0, MFT_consumed_share_thioether = 0, MFT_consumed_share_oligomerisation = 0, MFT_consumed_share_unassigned_decay = 0.9432, MFT_consumed_share_thiolate_oxidation = 0.0006706, MFT_consumed_share_protein_disulfide = 0, MFT_share_intact_skeleton_via_neutral_H2S = 0.4505

**cerny_ternary**: MFT_share_thiamine_route = 0.4524, MFT_share_sugar_routes = 0.5476, MFT_share_intact_skeleton_route = 0.0008301, MFT_share_norfuraneol_route = 0.5468, MFT_share_C2_plus_C3_route = 4.175e-11, FFT_share_of_furfural_flux = 0.007554, MFT_consumed_share_dimerisation = 0, MFT_consumed_share_MMFT = 0.002665, MFT_consumed_share_thioether = 0, MFT_consumed_share_oligomerisation = 0, MFT_consumed_share_unassigned_decay = 0.9284, MFT_consumed_share_thiolate_oxidation = 0.0002703, MFT_consumed_share_protein_disulfide = 0, MFT_share_intact_skeleton_via_neutral_H2S = 0.6967

## Contradiction found in the declaration

Cerny 2007 Table 4's MFT column (FIT) and Cerny 2007 Table 5's 1x arm (HOLD-OUT) are the SAME NUMBER: 85% thiamine-derived / 85:15. The declaration's disjointness rule 1 is violated by that pair.

Table 4's MFT column is EXCLUDED from this objective. The fit uses Cerny 2007b's full-ternary 54:46 (a different paper and a different system composition) and Table 4's non-MFT species. Reported, not resolved silently.

## Out of scope for this wave

- **pyrazines** -- strands: Amrani-Hemaimi 1995 Table 2's 40 isotope-fraction cells (a FIT row) and the alanine-vs-glycine ethyl-pyrazine ON/OFF switch (a star HOLD-OUT), plus Zhou 2023 Table 1's pyrazine and methylpyrazine columns. Reason: Pyrazines are a nitrogen-heterocycle lane with its own dicarbonyl and amino-acid bookkeeping; nothing in it is a thiol. Adding it would double the network for no gain on Modules 1-2. Deferred to a later wave, EXPLICITLY -- the on/off switch hold-out is NOT scored below and is not claimed.
- **Sotolon and 2-acetyl-2-thiazoline** -- strands: Hofmann 1996b's Sotolon anchors (13.5 / 764.7 / 273.1 ug, a FIT row) and the oxidant series (a star HOLD-OUT). Reason: Sotolon's precursors (butane-2,3-dione + hydroxyacetaldehyde) and 2-acetyl-2-thiazoline's (HDT/ATD) are a separate sub-network, and the AT numbers are figure-derived LOWER BOUNDS whose Ea inequality (sec. B10.4) is explicitly licensed as a MODEL CHANGE and not a parameter transfer. Deferred.
- **butane-2,3-dione and the Cerny 2004 in-situ ratios** -- strands: Cerny 2004's 54:46 butanedione split, 65:35 thiazole and 87:13 methylthio splits (FIT rows). Reason: Butane-2,3-dione is not in the state vector. The one Cerny 2004 constraint this module DOES carry is structural rather than numeric: 'the C2+C3 route was not relevant at 95 C', which is why that route competes rather than being assigned a share.
- **the Bornhorst norfuraneol Ea and the whole alkaline block** -- strands: Bornhorst 2017/2017b Ea, matrix pair and structural zero. Reason: pH 8.4-9.5 against this module's pH 4.5-7 -- the alkaline-pH wall, sec. C.13, 'neither the rates nor the M-2 values transfer to a pH-5 sulfur benchmark'. Carried as an UNVALIDATED PRIOR in parameters_sulfur.ALKALINE_PRIORS and operative nowhere.
- **ribose vs xylose** -- strands: the 1.38x MFT / 1.26x FFT gap between them at pH 5. Reason: One generic aldopentose. No step in the corpus is measured separately for the two sugars, so a per-sugar constant would be fitted to exactly the two numbers it then reproduces. The gap is reported as a fit residual, not absorbed.
