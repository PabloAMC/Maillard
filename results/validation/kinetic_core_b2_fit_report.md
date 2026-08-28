# Kinetic core, Build Wave B2 -- the sulfur fit

Modules 1 (sulfur formation) and 2 (thiol consumption) of `docs/reference/FIT_HOLDOUT_DECLARATION.md`.

- network: **42 species, 67 reactions** (15 inherited from B1's trunk, 52 new), carbon, nitrogen AND sulfur balance enforced at import
- objective: **39 declared FIT rows**, **36 free parameters**, final cost **16.195**, reduced chi-square **10.80**
- fitted pH parameters: **0**. The pH shape is structural (acid/base catalysis + measured H2S speciation + Zheng & Ho's measured four-pH thermolysis set).
- branch-fraction constants: **0**. Every split is a ratio of time-integrated mass-action fluxes.
- activation energies in the CONSUMPTION lane: **0**, by policy (inventory sec. C.1 / B.7).

## Read this before reading anything else

**0 of 36 parameters are individually identified.** With 39 rows against 36 free parameters the fit has 3 degrees of freedom, so the row-by-row agreement below is NOT evidence that the model is right. What the corpus determines is a set of RATIOS and SHAPES -- branch shares, MFT/FFT, dimer/monomer, the pH ordering -- and those are what the hold-out report scores. No individual rate constant in this table should be quoted as a measured quantity, cited elsewhere, or carried into another module as if it were one.

**The single lumped formation activation energy landed ON ITS BOUND at 20.0 kJ/mol.** That is not a barrier estimate, it is the fit saying it wants LESS temperature dependence than the bound allows. The reason is visible in the corpus: the fit panel spans 115-145 C but each temperature comes with a DIFFERENT feedstock, buffer and lab, so temperature is confounded with system and a cross-system Ea is not identifiable even in principle. Reported as a failure to determine, not as a value.

## Row-by-row

| row | observed | predicted | fold | source |
|---|---:|---:|---:|---|
| `hofmann_ribose_FFT` | 0.00106 | 0.001002 | 0.95x | Hofmann 1998 T1, ribose pH 5.0: 12.1 ug/100 mL = 121 ppb |
| `hofmann_ribose_MFT` | 0.001734 | 0.002673 | 1.54x | Hofmann 1998 T1, ribose pH 5.0: 19.8 ug/100 mL = 198 ppb |
| `hofmann_xylose_FFT` | 0.0008409 | 0.001002 | 1.19x | Hofmann 1998 T1, xylose pH 5.0: 9.6 ug/100 mL = 96 ppb |
| `hofmann_xylose_MFT` | 0.001253 | 0.002673 | 2.13x | Hofmann 1998 T1, xylose pH 5.0: 14.3 ug/100 mL = 143 ppb |
| `hofmann_glucose_FFT` | 0.0002453 | 0.0002131 | 0.87x | Hofmann 1998 T1, glucose pH 5.0: 2.8 ug/100 mL = 28 ppb |
| `hofmann_glucose_MFT` | 0.0001664 | 0.000318 | 1.91x | Hofmann 1998 T1, glucose pH 5.0: 1.9 ug/100 mL = 19 ppb |
| `hofmann_fructose_FFT` | 0.0002803 | 0.000321 | 1.15x | Hofmann 1998 T1, fructose pH 5.0: 3.2 ug/100 mL = 32 ppb |
| `hofmann_fructose_MFT` | 0.000219 | 0.0001212 | 0.55x | Hofmann 1998 T1, fructose pH 5.0: 2.5 ug/100 mL = 25 ppb |
| `hofmann_ribose_NF_insitu` | 4.779 | 3.21 | 0.67x | Hofmann 1998 T5, ribose pH 5.0: 54 530 ug/100 mL = 4.78 mmol/L |
| `hofmann_ribose_FUR_insitu` | 0.007025 | 0.006155 | 0.88x | Hofmann 1998 T5, ribose pH 5.0: 67.5 ug/100 mL furan-2-aldehyde |
| `fed_ribose_h2s_FFT` | 0.008 | 0.007774 | 0.97x | Hofmann 1998 T3: ribose + H2S -> FFT, 9.2 ug, 0.008 mol% |
| `fed_ribose_h2s_MFT` | 0.01 | 0.009428 | 0.94x | Hofmann 1998 T6: ribose + H2S -> MFT, 15.1 ug, 0.01 mol% |
| `fed_tdp_h2s_FFT` | 0.08 | 0.5766 | 7.21x | Hofmann 1998 T3: 3-deoxyribosulose + H2S -> FFT, 78.6 ug, 0.08 mol% |
| `fed_furfural_h2s_FFT` | 0.48 | 0.5 | 1.04x | Hofmann 1998 T3: furan-2-aldehyde + H2S -> FFT, 550.8 ug, 0.48 mol% |
| `fed_nf_h2s_MFT` | 0.19 | 0.3543 | 1.86x | Hofmann 1998 T4: norfuraneol + H2S -> MFT, 211.2 ug, 0.19 mol% |
| `fed_nf_cys_MFT` | 0.05 | 0.01685 | 0.34x | Hofmann 1998 T4: norfuraneol + cysteine -> MFT, 50.8 ug, 0.05 mol% |
| `fed_c2c3_MFT` | 0.24 | 0.1022 | 0.43x | Hofmann 1998 T10: C2 + C3 -> MFT, 268.1 ug, 0.24 mol% |
| `fed_c2c3_MFT_pH3` | 0.01 | 0.1022 | 10.22x | Hofmann 1998 T10 pH ladder: 15.5 ug at pH 3.0 |
| `fed_c2c3_MFT_pH7` | 0.27 | 0.1022 | 0.38x | Hofmann 1998 T10 pH ladder: 311.5 ug at pH 7.0 |
| `fed_thiamine_MFT` | 0.01 | 0.0165 | 1.65x | Hofmann 1998 T8: thiamin -> MFT, 8.2 ug, 0.01 mol% |
| `fed_mgo_h2s_1to1_MP` | 1.8 | 1.964 | 1.09x | Hofmann 1998 T7: 2-oxopropanal + H2S 1:1 -> 1650 ug, 1.8 mol% |
| `fed_mgo_h2s_1to2_MP` | 4 | 3.422 | 0.86x | Hofmann 1998 T7: 2-oxopropanal + H2S 1:2 -> 3600 ug, 4.0 mol% |
| `zhou_pH7_MFT` | 0.01392 | 0.001818 | 0.13x | Zhou 2023 T1 pH 7: MFT 1588.57 +/- 21.24 ug/L |
| `zhou_pH7_FFT` | 0.006639 | 0.0007669 | 0.12x | Zhou 2023 T1 pH 7: FFT 757.965 +/- 13.03 ug/L |
| `zhou_pH7_MFT_over_FFT` | 2.096 | 2.37 | 1.13x | Zhou 2023 T1 pH 7, MFT/FFT = 2.096 [D] |
| `zhou_pH7_dimer_over_MFT` | 0.0323 | 0.006889 | 0.21x | Zhou 2023 T1 pH 7: 102.59 ug/L dimer against 1588.57 MFT (molar ratio  |
| `zhou_pH7_ACTZ` | 9.201e-05 | 9.308e-05 | 1.01x | Zhou 2023 T1 pH 7: 2-acetylthiazole 11.70 +/- 2.14 ug/L |
| `zhou_pH7_FUR_arp_alone` | 0.01394 | 0.01607 | 1.15x | Zhou 2023 T1 pH 7, ARP ALONE: 2-furfural 1339.37 +/- 83.04 ug/L |
| `whitfield_nf_cys_MFT` | 0.15 | 0.03634 | 0.24x | Whitfield & Mottram 1999 T1: NF + cysteine, pH 4.5, 140 C, 0.150 mol%  |
| `whitfield_nf_h2s_MFT` | 0.12 | 0.468 | 3.90x | Whitfield & Mottram 1999 T1: NF + H2S 1:2, pH 4.5, 0.120 mol% |
| `whitfield_mercaptoketone_over_MFT` | 16.3 | 18.01 | 1.10x | Whitfield 1999: mercaptoketones : MFT = 16.3 : 1 from fed NF |
| `vanseeventer_cys_conversion` | 0.55 | 0.2123 | 0.39x | van Seeventer 2001: cysteine 33 -> 15 mM (55%), 130 C / 20 min |
| `vanseeventer_ribose_conversion` | 0.75 | 0.5286 | 0.70x | van Seeventer 2001: ribose 100 -> 25 mM (75%), 130 C / 20 min |
| `zhang_fig1_cys_dimer_over_MFT` | 0.0429 | 0.06301 | 1.47x | Zhang 2024 Fig. 1: 0.115 ng/mL dimer against 1.34 MFT = 8.6% BY MASS = |
| `zhang_fig1_gcys_dimer_over_MFT` | 0.2711 | 1.763 | 6.50x | Zhang 2024 Fig. 1: 1.09 against 2.01 = 54.2% by mass = 0.2711 molar |
| `zhang_fig1_cys_MMFT_over_MFT` | 0.0213 | 0.02302 | 1.08x | Zhang 2024 Fig. 1: MMFT 0.04 against MFT 1.34 ng/mL = 0.0213 molar |
| `cerny_ternary_thiamine_share` | 0.54 | 0.3203 | 0.59x | Cerny 2007b Table 3, full ternary: 54% unlabelled (thiamine) : 46% 13C |
| `cerny_isomer_split` | 1 | 1.266 | 1.27x | Cerny 2007 T4: 2-mercapto-3-pentanone is 94->95% XYLOSE-derived while  |
| `yaghmur_fft_share_ceiling` | 0.012 | 0.0143 | 1.19x | Yaghmur et al. 2005 sec. 3.1 p. 226 (<1% FFT conversion of the furfura |

## Parameters

| parameter | log10 k(145 C) | 95% half-width | identified? |
|---|---:|---:|---|
| `k_pent_dpo` | -5.368 | -- | **UNIDENTIFIED** |
| `k_pent_tdp` | -4.590 | -- | **UNIDENTIFIED** |
| `k_dpo_c2c3` | -4.519 | -- | **UNIDENTIFIED** |
| `k_arp_dpo` | -2.920 | -- | **UNIDENTIFIED** |
| `k_arp_tdp` | -0.899 | -- | **UNIDENTIFIED** |
| `k_dpo_nf` | -2.231 | -- | **UNIDENTIFIED** |
| `k_dpo_ptr` | -3.417 | -- | **UNIDENTIFIED** |
| `k_dpo_ddp` | -3.254 | -- | **UNIDENTIFIED** |
| `k_tdp_fur` | -2.250 | -- | **UNIDENTIFIED** |
| `k_ddp_mft` | -4.621 | -- | **UNIDENTIFIED** |
| `k_fur_fft` | -1.525 | -- | **UNIDENTIFIED** |
| `k_nf_mft` | -4.548 | -- | **UNIDENTIFIED** |
| `k_nf_mp3p` | -2.288 | -- | **UNIDENTIFIED** |
| `k_mgo_mp` | -3.270 | -- | **UNIDENTIFIED** |
| `k_ha_mp_mft` | -4.942 | -- | **UNIDENTIFIED** |
| `k_glc_ha` | -0.750 | -- | **UNIDENTIFIED** |
| `k_thi_hmp` | -1.504 | -- | **UNIDENTIFIED** |
| `k_thi_mesh` | -0.598 | -- | **UNIDENTIFIED** |
| `k_hmp_mft` | -3.300 | -- | **UNIDENTIFIED** |
| `k_hmp_mp2p` | -1.008 | -- | **UNIDENTIFIED** |
| `k_cys_actz` | -2.227 | -- | **UNIDENTIFIED** |
| `k_dimer_mft` | -0.594 | -- | **UNIDENTIFIED** |
| `k_dimer_fft` | -3.853 | -- | **UNIDENTIFIED** |
| `k_mmft` | -3.281 | -- | **UNIDENTIFIED** |
| `k_mft_decay` | -2.540 | -- | **UNIDENTIFIED** |
| `k_fft_decay` | -0.726 | -- | **UNIDENTIFIED** |
| `k_dimer_decay` | -6.518 | -- | **UNIDENTIFIED** |
| `k_nf_decay` | -3.586 | -- | **UNIDENTIFIED** |
| `k_fur_decay` | 0.332 | -- | **UNIDENTIFIED** |
| `k_h2s_loss` | -6.621 | -- | **UNIDENTIFIED** |
| `k_osone_decay` | -9.029 | -- | **UNIDENTIFIED** |
| `k_thiol_decay` | -0.677 | -- | **UNIDENTIFIED** |
| `k_pent_caramel` | -2.909 | -- | **UNIDENTIFIED** |
| `k_pent_thermal` | -1.352 | -- | **UNIDENTIFIED** |
| `k_glc_fur` | -3.029 | -- | **UNIDENTIFIED** |
| `Ea_lumped_formation` | 20.000 | -- | **UNIDENTIFIED** |

## Branch shares at the fitted point (computed, never stored)

**cerny_ternary**: MFT_share_thiamine_route = 0.3203, MFT_share_sugar_routes = 0.6797, MFT_share_intact_skeleton_route = 0.5896, MFT_share_norfuraneol_route = 0.09006, MFT_share_C2_plus_C3_route = 2.108e-09, FFT_share_of_furfural_flux = 0.03494, MFT_consumed_share_dimerisation = 0, MFT_consumed_share_MMFT = 0.05112, MFT_consumed_share_thioether = 0, MFT_consumed_share_oligomerisation = 0, MFT_consumed_share_unassigned_decay = 0.02095

**hofmann_pentose_pH5**: MFT_share_thiamine_route = 0, MFT_share_sugar_routes = 1, MFT_share_intact_skeleton_route = 0.7049, MFT_share_norfuraneol_route = 0.2951, MFT_share_C2_plus_C3_route = 4.388e-09, FFT_share_of_furfural_flux = 0.0143, MFT_consumed_share_dimerisation = 0.003543, MFT_consumed_share_MMFT = 0, MFT_consumed_share_thioether = 0, MFT_consumed_share_oligomerisation = 0, MFT_consumed_share_unassigned_decay = 0.01296

## Contradiction found in the declaration

Cerny 2007 Table 4's MFT column (FIT) and Cerny 2007 Table 5's 1x arm (HOLD-OUT) are the SAME NUMBER: 85% thiamine-derived / 85:15. The declaration's disjointness rule 1 is violated by that pair.

Table 4's MFT column is EXCLUDED from this objective. The fit uses Cerny 2007b's full-ternary 54:46 (a different paper and a different system composition) and Table 4's non-MFT species. Reported, not resolved silently.

## Out of scope for this wave

- **pyrazines** -- strands: Amrani-Hemaimi 1995 Table 2's 40 isotope-fraction cells (a FIT row) and the alanine-vs-glycine ethyl-pyrazine ON/OFF switch (a star HOLD-OUT), plus Zhou 2023 Table 1's pyrazine and methylpyrazine columns. Reason: Pyrazines are a nitrogen-heterocycle lane with its own dicarbonyl and amino-acid bookkeeping; nothing in it is a thiol. Adding it would double the network for no gain on Modules 1-2. Deferred to a later wave, EXPLICITLY -- the on/off switch hold-out is NOT scored below and is not claimed.
- **Sotolon and 2-acetyl-2-thiazoline** -- strands: Hofmann 1996b's Sotolon anchors (13.5 / 764.7 / 273.1 ug, a FIT row) and the oxidant series (a star HOLD-OUT). Reason: Sotolon's precursors (butane-2,3-dione + hydroxyacetaldehyde) and 2-acetyl-2-thiazoline's (HDT/ATD) are a separate sub-network, and the AT numbers are figure-derived LOWER BOUNDS whose Ea inequality (sec. B10.4) is explicitly licensed as a MODEL CHANGE and not a parameter transfer. Deferred.
- **butane-2,3-dione and the Cerny 2004 in-situ ratios** -- strands: Cerny 2004's 54:46 butanedione split, 65:35 thiazole and 87:13 methylthio splits (FIT rows). Reason: Butane-2,3-dione is not in the state vector. The one Cerny 2004 constraint this module DOES carry is structural rather than numeric: 'the C2+C3 route was not relevant at 95 C', which is why that route competes rather than being assigned a share.
- **the Bornhorst norfuraneol Ea and the whole alkaline block** -- strands: Bornhorst 2017/2017b Ea, matrix pair and structural zero. Reason: pH 8.4-9.5 against this module's pH 4.5-7 -- the alkaline-pH wall, sec. C.13, 'neither the rates nor the M-2 values transfer to a pH-5 sulfur benchmark'. Carried as an UNVALIDATED PRIOR in parameters_sulfur.ALKALINE_PRIORS and operative nowhere.
- **ribose vs xylose** -- strands: the 1.38x MFT / 1.26x FFT gap between them at pH 5. Reason: One generic aldopentose. No step in the corpus is measured separately for the two sugars, so a per-sugar constant would be fitted to exactly the two numbers it then reproduces. The gap is reported as a fit residual, not absorbed.
