# Kinetic core, Build Wave B2.1 -- the REVISED sulfur fit

Modules 1 (sulfur formation) and 2 (thiol consumption) of `docs/reference/FIT_HOLDOUT_DECLARATION.md`.

- network: **46 species, 79 reactions** (15 inherited from B1's trunk, 64 new), carbon, nitrogen AND sulfur balance enforced at import
- objective: **55 declared FIT rows**, **44 free parameters**, final cost **11.078**, reduced chi-square **2.01**
- fitted pH parameters: **0**. The pH shape is structural (acid/base catalysis + measured H2S speciation + Zheng & Ho's measured four-pH thermolysis set).
- branch-fraction constants: **0**. Every split is a ratio of time-integrated mass-action fluxes.
- activation energies in the CONSUMPTION lane: **0**, by policy (inventory sec. C.1 / B.7).

## Read this before reading anything else

**7 of 44 parameters are individually identified.** With 55 rows against 44 free parameters the fit has 11 degrees of freedom, so the row-by-row agreement below is NOT evidence that the model is right. What the corpus determines is a set of RATIOS and SHAPES -- branch shares, MFT/FFT, dimer/monomer, the pH ordering -- and those are what the hold-out report scores. No individual rate constant in this table should be quoted as a measured quantity, cited elsewhere, or carried into another module as if it were one.

## Row-by-row

| row | observed | predicted | fold | source |
|---|---:|---:|---:|---|
| `hofmann_ribose_FFT` | 0.00106 | 0.0009159 | 0.86x | Hofmann 1998 T1, ribose pH 5.0: 12.1 ug/100 mL = 121 ppb |
| `hofmann_ribose_MFT` | 0.001734 | 0.001456 | 0.84x | Hofmann 1998 T1, ribose pH 5.0: 19.8 ug/100 mL = 198 ppb |
| `hofmann_xylose_FFT` | 0.0008409 | 0.0009159 | 1.09x | Hofmann 1998 T1, xylose pH 5.0: 9.6 ug/100 mL = 96 ppb |
| `hofmann_xylose_MFT` | 0.001253 | 0.001456 | 1.16x | Hofmann 1998 T1, xylose pH 5.0: 14.3 ug/100 mL = 143 ppb |
| `hofmann_glucose_FFT` | 0.0002453 | 0.0002217 | 0.90x | Hofmann 1998 T1, glucose pH 5.0: 2.8 ug/100 mL = 28 ppb |
| `hofmann_glucose_MFT` | 0.0001664 | 0.0003102 | 1.86x | Hofmann 1998 T1, glucose pH 5.0: 1.9 ug/100 mL = 19 ppb |
| `hofmann_fructose_FFT` | 0.0002803 | 0.0003102 | 1.11x | Hofmann 1998 T1, fructose pH 5.0: 3.2 ug/100 mL = 32 ppb |
| `hofmann_fructose_MFT` | 0.000219 | 0.0001307 | 0.60x | Hofmann 1998 T1, fructose pH 5.0: 2.5 ug/100 mL = 25 ppb |
| `hofmann_ribose_NF_insitu` | 4.779 | 4.115 | 0.86x | Hofmann 1998 T5, ribose pH 5.0: 54 530 ug/100 mL = 4.78 mmol/L |
| `hofmann_ribose_FUR_insitu` | 0.007025 | 0.007161 | 1.02x | Hofmann 1998 T5, ribose pH 5.0: 67.5 ug/100 mL furan-2-aldehyde |
| `fed_ribose_h2s_FFT` | 0.008 | 0.006879 | 0.86x | Hofmann 1998 T3: ribose + H2S -> FFT, 9.2 ug, 0.008 mol% |
| `fed_ribose_h2s_MFT` | 0.01 | 0.01181 | 1.18x | Hofmann 1998 T6: ribose + H2S -> MFT, 15.1 ug, 0.01 mol% |
| `fed_tdp_h2s_FFT` | 0.08 | 0.1363 | 1.70x | Hofmann 1998 T3: 3-deoxyribosulose + H2S -> FFT, 78.6 ug, 0.08 mol% |
| `fed_furfural_h2s_FFT` | 0.48 | 0.9308 | 1.94x | Hofmann 1998 T3: furan-2-aldehyde + H2S -> FFT, 550.8 ug, 0.48 mol% |
| `fed_nf_h2s_MFT` | 0.19 | 0.3554 | 1.87x | Hofmann 1998 T4: norfuraneol + H2S -> MFT, 211.2 ug, 0.19 mol% |
| `fed_nf_cys_MFT` | 0.05 | 0.05413 | 1.08x | Hofmann 1998 T4: norfuraneol + cysteine -> MFT, 50.8 ug, 0.05 mol% |
| `fed_c2c3_MFT` | 0.24 | 0.09747 | 0.41x | Hofmann 1998 T10: C2 + C3 -> MFT, 268.1 ug, 0.24 mol% |
| `fed_c2c3_MFT_pH3` | 0.01 | 0.09747 | 9.75x | Hofmann 1998 T10 pH ladder: 15.5 ug at pH 3.0 |
| `fed_c2c3_MFT_pH7` | 0.27 | 0.09747 | 0.36x | Hofmann 1998 T10 pH ladder: 311.5 ug at pH 7.0 |
| `fed_thiamine_MFT` | 0.01 | 0.009776 | 0.98x | Hofmann 1998 T8: thiamin -> MFT, 8.2 ug, 0.01 mol% |
| `fed_mgo_h2s_1to1_MP` | 1.8 | 2.793 | 1.55x | Hofmann 1998 T7: 2-oxopropanal + H2S 1:1 -> 1650 ug, 1.8 mol% |
| `fed_mgo_h2s_1to2_MP` | 4 | 2.839 | 0.71x | Hofmann 1998 T7: 2-oxopropanal + H2S 1:2 -> 3600 ug, 4.0 mol% |
| `zhou_pH7_MFT` | 0.01392 | 0.006858 | 0.49x | Zhou 2023 T1 pH 7: MFT 1588.57 +/- 21.24 ug/L |
| `zhou_pH7_FFT` | 0.006639 | 0.003953 | 0.60x | Zhou 2023 T1 pH 7: FFT 757.965 +/- 13.03 ug/L |
| `zhou_pH7_MFT_over_FFT` | 2.096 | 1.735 | 0.83x | Zhou 2023 T1 pH 7, MFT/FFT = 2.096 [D] |
| `zhou_pH7_dimer_over_MFT` | 0.0323 | 0.02429 | 0.75x | Zhou 2023 T1 pH 7: 102.59 ug/L dimer against 1588.57 MFT (molar ratio  |
| `zhou_pH7_ACTZ` | 9.201e-05 | 9.209e-05 | 1.00x | Zhou 2023 T1 pH 7: 2-acetylthiazole 11.70 +/- 2.14 ug/L |
| `zhou_pH7_FUR_arp_alone` | 0.01394 | 0.007245 | 0.52x | Zhou 2023 T1 pH 7, ARP ALONE: 2-furfural 1339.37 +/- 83.04 ug/L |
| `whitfield_nf_cys_MFT` | 0.15 | 0.04253 | 0.28x | Whitfield & Mottram 1999 T1: NF + cysteine, pH 4.5, 140 C, 0.150 mol%  |
| `whitfield_nf_h2s_MFT` | 0.12 | 0.08795 | 0.73x | Whitfield & Mottram 1999 T1: NF + H2S 1:2, pH 4.5, 0.120 mol% |
| `whitfield_mercaptoketone_over_MFT` | 16.3 | 12.62 | 0.77x | Whitfield 1999: mercaptoketones : MFT = 16.3 : 1 from fed NF |
| `vanseeventer_cys_conversion` | 0.55 | 0.4542 | 0.83x | van Seeventer 2001: cysteine 33 -> 15 mM (55%), 130 C / 20 min |
| `vanseeventer_ribose_conversion` | 0.75 | 1 | 1.33x | van Seeventer 2001: ribose 100 -> 25 mM (75%), 130 C / 20 min |
| `zhang_fig1_cys_dimer_over_MFT` | 0.0429 | 0.02162 | 0.50x | Zhang 2024 Fig. 1: 0.115 ng/mL dimer against 1.34 MFT = 8.6% BY MASS = |
| `zhang_fig1_gcys_dimer_over_MFT` | 0.2711 | 0.9615 | 3.55x | Zhang 2024 Fig. 1: 1.09 against 2.01 = 54.2% by mass = 0.2711 molar |
| `zhang_fig1_cys_MMFT_over_MFT` | 0.0213 | 0.02128 | 1.00x | Zhang 2024 Fig. 1: MMFT 0.04 against MFT 1.34 ng/mL = 0.0213 molar |
| `cerny_ternary_thiamine_share` | 0.54 | 0.5426 | 1.00x | Cerny 2007b Table 3, full ternary: 54% unlabelled (thiamine) : 46% 13C |
| `cerny_isomer_split` | 1 | 1 | 1.00x | Cerny 2007 T4: 2-mercapto-3-pentanone is 94->95% XYLOSE-derived while  |
| `kang_100C_MFT` | 1.084e-05 | 7.908e-06 | 0.73x | Kang 2026 SI T-S4, 100 C: MFT 1.237 ug/L (Tier A, R^2 0.9989) |
| `kang_120C_MFT` | 1.216e-05 | 2.33e-05 | 1.92x | Kang 2026 SI T-S4, 120 C: MFT 1.388 ug/L; the same value appears in th |
| `kang_100C_FFT` | 3.271e-05 | 1.599e-05 | 0.49x | Kang 2026 SI T-S4, 100 C: FFT 3.734 ug/L (Tier A, R^2 0.9992) |
| `kang_120C_FFT` | 3.598e-05 | 1.701e-05 | 0.47x | Kang 2026 SI T-S4, 120 C: FFT 4.107 ug/L |
| `kang_100C_FUR` | 3.519e-05 | 0.000245 | 6.96x | Kang 2026 SI T-S4, 100 C: furfural 3.381 ug/L |
| `kang_120C_FUR` | 6.029e-05 | 8.109e-05 | 1.34x | Kang 2026 SI T-S4, 120 C: furfural 5.793 ug/L |
| `kang_100C_MFT_over_FFT` | 0.3311 | 0.4944 | 1.49x | Kang 2026 SI sec. 7b: FFT/MFT = 3.02 at 100 C |
| `kang_100C_cys_conversion` | 0.162 | 0.1405 | 0.87x | Kang 2026 SI Fig. S4 (digitised): free-cysteine conversion 16.2% at 10 |
| `kang_120C_cys_conversion` | 0.387 | 0.392 | 1.01x | Kang 2026 SI Fig. S4 (digitised): 38.7% at 120 C / 120 min |
| `kang_free_cys_yield_ceiling` | 0.163 | 0.1032 | 0.63x | Kang 2026 SI Fig. S3 (digitised, kang2026_SI_extraction.md sec. 6c): p |
| `kumazawa_FFT_retention_pH3_0` | 0.995 | 0.7507 | 0.75x | Kumazawa 2003 Fig. 3: 99.5% of 1 ppm 2-furfurylthiol survives 121 C /  |
| `kumazawa_FFT_retention_pH4_0` | 0.962 | 0.7409 | 0.77x | Kumazawa 2003 Fig. 3: 96.2% of 1 ppm 2-furfurylthiol survives 121 C /  |
| `kumazawa_FFT_retention_pH5_0` | 0.891 | 0.6564 | 0.74x | Kumazawa 2003 Fig. 3: 89.1% of 1 ppm 2-furfurylthiol survives 121 C /  |
| `kumazawa_FFT_retention_pH5_4` | 0.795 | 0.5521 | 0.69x | Kumazawa 2003 Fig. 3: 79.5% of 1 ppm 2-furfurylthiol survives 121 C /  |
| `kumazawa_FFT_retention_pH6_0` | 0.451 | 0.3162 | 0.70x | Kumazawa 2003 Fig. 3: 45.1% of 1 ppm 2-furfurylthiol survives 121 C /  |
| `kumazawa_FFT_retention_pH6_4` | 0.11 | 0.1796 | 1.63x | Kumazawa 2003 Fig. 3: 11.0% of 1 ppm 2-furfurylthiol survives 121 C /  |
| `yaghmur_fft_share_ceiling` | 0.012 | 0.01373 | 1.14x | Yaghmur et al. 2005 sec. 3.1 p. 226 (<1% FFT conversion of the furfura |

## Parameters

| parameter | log10 k(145 C) | 95% half-width | identified? |
|---|---:|---:|---|
| `k_pent_dpo` | -5.769 | -- | **UNIDENTIFIED** |
| `k_pent_tdp` | -1.722 | -- | **UNIDENTIFIED** |
| `k_dpo_c2c3` | -1.727 | -- | **UNIDENTIFIED** |
| `k_arp_dpo` | -1.921 | -- | **UNIDENTIFIED** |
| `k_arp_tdp` | -0.297 | -- | **UNIDENTIFIED** |
| `k_dpo_nf` | -2.056 | -- | **UNIDENTIFIED** |
| `k_dpo_ptr` | -5.845 | -- | **UNIDENTIFIED** |
| `k_dpo_ddp` | -2.734 | -- | **UNIDENTIFIED** |
| `k_tdp_fur` | -2.132 | -- | **UNIDENTIFIED** |
| `k_ddp_mft` | -2.349 | -- | **UNIDENTIFIED** |
| `k_fur_fft` | 0.497 | -- | **UNIDENTIFIED** |
| `k_nf_mft` | -2.866 | -- | **UNIDENTIFIED** |
| `k_nf_mp3p` | -1.481 | -- | **UNIDENTIFIED** |
| `k_mgo_mp` | 0.486 | -- | **UNIDENTIFIED** |
| `k_ha_mp_mft` | -4.472 | 1.281 | yes |
| `k_glc_ha` | -2.476 | 2.989 | yes |
| `k_thi_hmp` | -2.902 | -- | **UNIDENTIFIED** |
| `k_thi_mesh` | -2.445 | -- | **UNIDENTIFIED** |
| `k_hmp_mft` | -2.703 | -- | **UNIDENTIFIED** |
| `k_hmp_mp2p` | -3.064 | -- | **UNIDENTIFIED** |
| `k_cys_actz` | -4.630 | -- | **UNIDENTIFIED** |
| `k_dimer_mft` | -0.811 | -- | **UNIDENTIFIED** |
| `k_dimer_fft` | -0.020 | 0.692 | yes |
| `k_mmft` | -2.650 | -- | **UNIDENTIFIED** |
| `k_mft_decay` | -1.021 | 0.647 | yes |
| `k_fft_decay` | -0.637 | 1.852 | yes |
| `k_dimer_decay` | -7.127 | -- | **UNIDENTIFIED** |
| `k_nf_decay` | -8.629 | -- | **UNIDENTIFIED** |
| `k_fur_decay` | 0.500 | -- | **UNIDENTIFIED** |
| `k_h2s_loss` | -0.396 | -- | **UNIDENTIFIED** |
| `k_osone_decay` | -1.034 | -- | **UNIDENTIFIED** |
| `k_thiol_decay` | -0.749 | 0.228 | yes |
| `k_pent_caramel` | -0.599 | -- | **UNIDENTIFIED** |
| `k_pent_thermal` | 0.500 | -- | **UNIDENTIFIED** |
| `k_glc_fur` | -4.471 | -- | **UNIDENTIFIED** |
| `k_arp_tdp_th` | -3.408 | -- | **UNIDENTIFIED** |
| `k_arp_dpo_th` | -7.551 | -- | **UNIDENTIFIED** |
| `k_ddp_mft_hs` | -6.011 | -- | **UNIDENTIFIED** |
| `k_fur_fft_hs` | 0.500 | -- | **UNIDENTIFIED** |
| `k_ttca_cys` | -0.199 | -- | **UNIDENTIFIED** |
| `k_ttca_deg` | 0.500 | -- | **UNIDENTIFIED** |
| `k_cys_thermal` | -2.125 | 0.871 | yes |
| `k_thiolate_loss` | -7.452 | -- | **UNIDENTIFIED** |
| `Ea_lumped_formation` | 119.281 | -- | **UNIDENTIFIED** |

## Branch shares at the fitted point (computed, never stored)

**hofmann_pentose_pH5**: MFT_share_thiamine_route = 0, MFT_share_sugar_routes = 1, MFT_share_intact_skeleton_route = 0.7015, MFT_share_norfuraneol_route = 0.06704, MFT_share_C2_plus_C3_route = 0.2315, FFT_share_of_furfural_flux = 0.01373, MFT_consumed_share_dimerisation = 0.002726, MFT_consumed_share_MMFT = 0, MFT_consumed_share_thioether = 0, MFT_consumed_share_oligomerisation = 0, MFT_consumed_share_unassigned_decay = 0.6137, MFT_consumed_share_thiolate_oxidation = 2.274e-07, MFT_consumed_share_protein_disulfide = 0, MFT_share_intact_skeleton_via_neutral_H2S = 0.9998

**cerny_ternary**: MFT_share_thiamine_route = 0.5426, MFT_share_sugar_routes = 0.4574, MFT_share_intact_skeleton_route = 0.3585, MFT_share_norfuraneol_route = 0.01313, MFT_share_C2_plus_C3_route = 0.08583, FFT_share_of_furfural_flux = 0.01365, MFT_consumed_share_dimerisation = 0, MFT_consumed_share_MMFT = 0.02144, MFT_consumed_share_thioether = 0, MFT_consumed_share_oligomerisation = 0, MFT_consumed_share_unassigned_decay = 0.5754, MFT_consumed_share_thiolate_oxidation = 2.132e-07, MFT_consumed_share_protein_disulfide = 0, MFT_share_intact_skeleton_via_neutral_H2S = 0.9998

## Contradiction found in the declaration

Cerny 2007 Table 4's MFT column (FIT) and Cerny 2007 Table 5's 1x arm (HOLD-OUT) are the SAME NUMBER: 85% thiamine-derived / 85:15. The declaration's disjointness rule 1 is violated by that pair.

Table 4's MFT column is EXCLUDED from this objective. The fit uses Cerny 2007b's full-ternary 54:46 (a different paper and a different system composition) and Table 4's non-MFT species. Reported, not resolved silently.

## Out of scope for this wave

- **pyrazines** -- strands: Amrani-Hemaimi 1995 Table 2's 40 isotope-fraction cells (a FIT row) and the alanine-vs-glycine ethyl-pyrazine ON/OFF switch (a star HOLD-OUT), plus Zhou 2023 Table 1's pyrazine and methylpyrazine columns. Reason: Pyrazines are a nitrogen-heterocycle lane with its own dicarbonyl and amino-acid bookkeeping; nothing in it is a thiol. Adding it would double the network for no gain on Modules 1-2. Deferred to a later wave, EXPLICITLY -- the on/off switch hold-out is NOT scored below and is not claimed.
- **Sotolon and 2-acetyl-2-thiazoline** -- strands: Hofmann 1996b's Sotolon anchors (13.5 / 764.7 / 273.1 ug, a FIT row) and the oxidant series (a star HOLD-OUT). Reason: Sotolon's precursors (butane-2,3-dione + hydroxyacetaldehyde) and 2-acetyl-2-thiazoline's (HDT/ATD) are a separate sub-network, and the AT numbers are figure-derived LOWER BOUNDS whose Ea inequality (sec. B10.4) is explicitly licensed as a MODEL CHANGE and not a parameter transfer. Deferred.
- **butane-2,3-dione and the Cerny 2004 in-situ ratios** -- strands: Cerny 2004's 54:46 butanedione split, 65:35 thiazole and 87:13 methylthio splits (FIT rows). Reason: Butane-2,3-dione is not in the state vector. The one Cerny 2004 constraint this module DOES carry is structural rather than numeric: 'the C2+C3 route was not relevant at 95 C', which is why that route competes rather than being assigned a share.
- **the Bornhorst norfuraneol Ea and the whole alkaline block** -- strands: Bornhorst 2017/2017b Ea, matrix pair and structural zero. Reason: pH 8.4-9.5 against this module's pH 4.5-7 -- the alkaline-pH wall, sec. C.13, 'neither the rates nor the M-2 values transfer to a pH-5 sulfur benchmark'. Carried as an UNVALIDATED PRIOR in parameters_sulfur.ALKALINE_PRIORS and operative nowhere.
- **ribose vs xylose** -- strands: the 1.38x MFT / 1.26x FFT gap between them at pH 5. Reason: One generic aldopentose. No step in the corpus is measured separately for the two sugars, so a per-sugar constant would be fitted to exactly the two numbers it then reproduces. The gap is reported as a fit residual, not absorbed.
