# Core prediction uncertainty (Monte-Carlo envelope on the kinetic core)

n_samples = 200, seed = 0, CI level = 90 %.

* benchmarks with an envelope: **32** of 40 on the panel; matched rows **49**; refused rows 18
* mixed-population coverage: 11/49 (0.224)
* **honest literature coverage: 11/44 (0.250)**, median CI width 1.069 log10; 5 not evaluable; 0 fitted rows excluded
* out-of-sample literature coverage: 8/35 (5 not evaluable); rows the core fit read: {'hits': 3, 'total': 9, 'not_evaluable': 0}
* sampled priors 33, fixed 45; lanes with NO sampled fit uncertainty: none
* observable bands (K_aw, HS-SPME) applied by quantification family -- rows: headspace 2, extraction 37, undeclared 10; undeclared bundles: acrylamide_spi_extrusion_130C_ACSRef3, cys_ribose_140C_Hofmann1998, external_validation_li_2026_spi_wg_hme_control, external_validation_liu_2023_ppi_offnote_baseline, pea_isolate_40C_PratapSingh2021, pea_isolate_uht_140C_Trikusuma2019, resconi_2023_pbma_beef_identity_benchmark, soy_isolate_40C_PratapSingh2021, thiamine_cys_glucose_120C_Bolton1994

## Per panel

| panel | hits | total | rate | median width (log10) | not evaluable |
|---|---|---|---|---|---|
| external_matrix | 2 | 4 | 0.500 | 2.096 | 0 |
| maillard_path_holdout | 4 | 25 | 0.160 | 0.930 | 5 |
| trust_loop | 5 | 15 | 0.333 | 1.639 | 0 |

## Rows

| benchmark | panel | compound | unit | measured | point | p5 | p50 | p95 | inside | width | obs bands | lane | role |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| acrylamide_spi_extrusion_130C_ACSRef3 | trust_loop | acrylamide | ppb | 150.000 | 0.035 | 0.00726 | 0.038 | 0.234 | no | 1.508 | yes (undeclared) | acrylamide | predictive |
| cys_ribose_140C_Hofmann1998 | trust_loop | 2-methyl-3-furanthiol | ppb | 342.000 | 72.778 | 11.604 | 76.404 | 505.099 | yes | 1.639 | yes (undeclared) | sulfur | predictive |
| cys_ribose_140C_Hofmann1998 | trust_loop | 2-furfurylthiol | ppb | 200.000 | 78.932 | 10.021 | 84.954 | 924.890 | yes | 1.965 | yes (undeclared) | sulfur | predictive |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | ppb | 32.000 | 338.275 | 32.464 | 492.714 | 3.31e+03 | no | 2.008 | no (extraction) | sulfur | predictive [in core fit] |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 25.000 | 10.463 | 1.359 | 6.539 | 16.633 | no | 1.088 | no (extraction) | sulfur | predictive [in core fit] |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | ppb | 28.000 | 236.804 | 8.446 | 430.649 | 3.49e+03 | yes | 2.617 | no (extraction) | sulfur | predictive [in core fit] |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 19.000 | 33.088 | 4.517 | 19.537 | 50.626 | yes | 1.050 | no (extraction) | sulfur | predictive [in core fit] |
| hofmann1998_norfuraneol_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 1.02e+03 | 2.32e+03 | 1.48e+03 | 2.25e+03 | 2.73e+03 | no | 0.267 | no (extraction) | sulfur | predictive [in core fit] |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | ppb | 121.000 | 1.21e+03 | 233.237 | 1.13e+03 | 6.52e+03 | no | 1.447 | no (extraction) | sulfur | predictive [in core fit] |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 198.000 | 697.088 | 385.854 | 668.033 | 1.08e+03 | no | 0.447 | no (extraction) | sulfur | predictive [in core fit] |
| pea_isolate_40C_PratapSingh2021 | trust_loop | hexanal | ppb | 1.14e+03 | 0.339 | 0.036 | 0.369 | 4.381 | no | 2.089 | yes (undeclared) | lipid | predictive |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | hexanal | ppb | 782.000 | 22.839 | 1.223 | 25.195 | 459.646 | no | 2.575 | yes (undeclared) | lipid | predictive |
| resconi_2023_pbma_beef_identity_benchmark | trust_loop | furfural | ppb | 715.220 | 9.46e-05 | 2.74e-15 | 0.000578 | 2.058 | no | 14.876 | yes (undeclared) | sulfur | predictive |
| soy_isolate_40C_PratapSingh2021 | trust_loop | hexanal | ppb | 1.62e+03 | 0.267 | 0.028 | 0.290 | 3.449 | no | 2.089 | yes (undeclared) | lipid | predictive |
| thiamine_cys_glucose_120C_Bolton1994 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 13.000 | 17.383 | 3.569 | 18.677 | 115.067 | yes | 1.508 | yes (undeclared) | sulfur | predictive |
| mp_holdout_fructose_asparagine_180C_Lin2022 | maillard_path_holdout | Acrylamide | ppb | 1.86e+03 | 225.440 | 79.043 | 207.025 | 689.414 | no | 0.941 | no (extraction) | acrylamide | predictive |
| mp_holdout_fructose_asparagine_180C_Lin2022 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | ppb | 1.23e+04 | 1.94e+03 | 1.94e+03 | 1.94e+03 | 1.94e+03 | no | 1.14e-11 | no (extraction) | acrylamide | predictive |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | maillard_path_holdout | DMHF | ppb | 1.15e+03 | 21.841 | 7.588 | 21.565 | 64.513 | no | 0.930 | no (extraction) | acrylamide | predictive |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | ppb | 5.73e+04 | 2.8e+04 | 2.8e+04 | 2.8e+04 | 2.8e+04 | no | 2.8e-10 | no (extraction) | acrylamide | predictive |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | maillard_path_holdout | DMHF | ppb | 5.89e+03 | 21.841 | 7.588 | 21.565 | 64.513 | no | 0.930 | no (extraction) | acrylamide | predictive |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | ppb | 1.01e+05 | 2.8e+04 | 2.8e+04 | 2.8e+04 | 2.8e+04 | no | 2.8e-10 | no (extraction) | acrylamide | predictive |
| mp_holdout_glucose_asparagine_180C_10min_Chang2021 | maillard_path_holdout | Acrylamide | ppb | 28.000 | 6.75e+03 | 2.37e+03 | 6.2e+03 | 1.79e+04 | no | 0.878 | no (extraction) | acrylamide | predictive |
| mp_holdout_glucose_asparagine_180C_30min_Chang2021 | maillard_path_holdout | Acrylamide | ppb | 1.46e+03 | 4.03e+03 | 1.41e+03 | 3.7e+03 | 1.23e+04 | yes | 0.941 | no (extraction) | acrylamide | predictive |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | maillard_path_holdout | Acrylamide | ppb | 832.000 | 4.03e+03 | 1.41e+03 | 3.7e+03 | 1.23e+04 | no | 0.941 | no (extraction) | acrylamide | predictive |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | ppb | 7e+03 | 2.11e+03 | 2.11e+03 | 2.11e+03 | 2.11e+03 | no | 3.98e-10 | no (extraction) | acrylamide | predictive |
| mp_holdout_glucose_asparagine_180C_Ye2024 | maillard_path_holdout | Acrylamide | umol_per_mol_limiting_precursor | 140.580 | 7.04e+03 | 2.47e+03 | 6.46e+03 | 2.15e+04 | no | 0.941 | no (extraction) | acrylamide | predictive |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | ppb | 1.74e+04 | 1.46e+03 | 1.46e+03 | 1.46e+03 | 1.46e+03 | no | 8.25e-09 | no (extraction) | trunk | predictive |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 7.000 | 93.626 | 5.114 | 187.290 | 1.26e+03 | yes | 2.390 | no (extraction) | sulfur | predictive |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 3.000 | 33.217 | 4.540 | 19.784 | 50.780 | no | 1.049 | no (extraction) | sulfur | predictive |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 6.000 | 0.117 | 0.000138 | 0.147 | 2.441 | no | 4.247 | no (extraction) | sulfur | predictive |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 4.000 | 22.940 | 2.179 | 12.052 | 34.407 | yes | 1.198 | no (extraction) | sulfur | predictive |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 229.000 | 1.63e+03 | 297.819 | 1.55e+03 | 8.11e+03 | no | 1.435 | no (extraction) | sulfur | predictive (shared: hofmann_ribose_pH3_FFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 553.000 | 22.314 | 11.967 | 21.364 | 36.427 | no | 0.483 | no (extraction) | sulfur | predictive (shared: hofmann_ribose_pH3_MFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 12.000 | 0.095 | 0.00509 | 0.094 | 0.889 | no | 2.242 | no (extraction) | sulfur | predictive (shared: hofmann_ribose_pH7_FFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 25.000 | 0.174 | 0.029 | 0.149 | 0.248 | no | 0.929 | no (extraction) | sulfur | predictive (shared: hofmann_ribose_pH7_MFT) |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 96.000 | 268.846 | 70.428 | 289.514 | 1.25e+03 | yes | 1.250 | no (extraction) | sulfur | predictive [in core fit] |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 143.000 | 526.834 | 302.821 | 497.524 | 720.876 | no | 0.377 | no (extraction) | sulfur | predictive [in core fit] |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 6.880 | 60.743 | 23.585 | 50.708 | 95.985 | no | 0.610 | no (extraction) | sulfur | predictive |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 1.280 | 634.638 | 309.860 | 713.210 | 1.39e+03 | no | 0.652 | no (extraction) | sulfur | predictive |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 3.290 | 92.511 | 44.499 | 81.401 | 134.034 | no | 0.479 | no (extraction) | sulfur | predictive |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 1.460 | 557.702 | 279.335 | 691.546 | 1.54e+03 | no | 0.740 | no (extraction) | sulfur | predictive |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 2.400 | 128.587 | 68.877 | 119.775 | 191.160 | no | 0.443 | no (extraction) | sulfur | predictive |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 1.680 | 392.614 | 192.731 | 517.118 | 1.38e+03 | no | 0.856 | no (extraction) | sulfur | predictive |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 1.710 | 163.367 | 88.788 | 160.725 | 266.976 | no | 0.478 | no (extraction) | sulfur | predictive |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 1.620 | 242.131 | 123.975 | 330.230 | 920.255 | no | 0.871 | no (extraction) | sulfur | predictive |
| external_validation_bi_2020_raw_pea_hexanal | external_matrix | hexanal | ppb | 1.26e+03 | 0.339 | 0.036 | 0.369 | 4.381 | no | 2.089 | yes (headspace) | lipid | predictive |
| external_validation_bi_2020_roasted_pea_hexanal | external_matrix | hexanal | ppb | 324.000 | 88.598 | 8.016 | 101.125 | 1.23e+03 | yes | 2.184 | yes (headspace) | lipid | predictive |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | hexanal | ppb | 605.600 | 69.695 | 5.987 | 64.082 | 759.830 | yes | 2.104 | yes (undeclared) | lipid | predictive |
| external_validation_liu_2023_ppi_offnote_baseline | external_matrix | hexanal | ppb | 1.13e+04 | 0.339 | 0.036 | 0.369 | 4.381 | no | 2.089 | yes (undeclared) | lipid | predictive |

## Refused rows

| benchmark | panel | compound | reason |
|---|---|---|---|
| cml_cel_commercial_pbma_Foods2023 | trust_loop | Nε-(Carboxymethyl)lysine (CML) | UNREPRESENTED TARGETS: Nε-(Carboxymethyl)lysine (CML) -- not a species in any core lane, and not on the named unrepresented-compound list either: the engine has no vocabulary entry for it. |
| cml_cel_commercial_pbma_Foods2023 | trust_loop | Nε-(Carboxyethyl)lysine (CEL) | UNREPRESENTED TARGETS: Nε-(Carboxyethyl)lysine (CEL) -- not a species in any core lane, and not on the named unrepresented-compound list either: the engine has no vocabulary entry for it. |
| furosine_extrusion_crossover_140C_RamirezJimenez2000 | trust_loop | furosine | UNREPRESENTED TARGETS: furosine -- not a species in any core lane, and not on the named unrepresented-compound list either: the engine has no vocabulary entry for it. |
| hofmann1998_c2c3_recombination_145C_20min_pH3 | trust_loop | 2-Methyl-3-furanthiol (MFT) | UNMAPPED PRECURSORS 'Hydroxyacetaldehyde', 'Mercapto-2-propanone': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge. |
| hofmann1998_c2c3_recombination_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | UNMAPPED PRECURSORS 'Hydroxyacetaldehyde', 'Mercapto-2-propanone': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge. |
| hofmann1998_c2c3_recombination_145C_20min_pH7 | trust_loop | 2-Methyl-3-furanthiol (MFT) | UNMAPPED PRECURSORS 'Hydroxyacetaldehyde', 'Mercapto-2-propanone': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge. |
| hofmann1998_furan2aldehyde_h2s_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | UNMAPPED PRECURSORS 'Furan-2-aldehyde', 'Hydrogen sulfide': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge. |
| hofmann1998_norfuraneol_h2s_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | UNMAPPED PRECURSORS 'Hydrogen sulfide': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge. |
| pea_isolate_40C_PratapSingh2021 | trust_loop | 2-pentylfuran | UNREPRESENTED TARGETS: 2-pentylfuran -- The B6 lipid lane exists, but 2-pentylfuran is NOT in Frankel 1989's six-product slate and no branch fraction for the linoleate -> alkylfuran route is measured anywhere in the fit  |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | 2-pentylfuran | UNREPRESENTED TARGETS: 2-pentylfuran -- The B6 lipid lane exists, but 2-pentylfuran is NOT in Frankel 1989's six-product slate and no branch fraction for the linoleate -> alkylfuran route is measured anywhere in the fit  |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | nonanal | UNREPRESENTED TARGETS: nonanal -- the lipid lane exists and nonanal is a species in it, but its ONLY parent is the OLEATE hydroperoxide pool and the oleate -> nonanal branch fraction is measured NOWHERE in the fit corpus |
| soy_isolate_40C_PratapSingh2021 | trust_loop | 2-pentylfuran | UNREPRESENTED TARGETS: 2-pentylfuran -- The B6 lipid lane exists, but 2-pentylfuran is NOT in Frankel 1989's six-product slate and no branch fraction for the linoleate -> alkylfuran route is measured anywhere in the fit  |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | maillard_path_holdout | Furfural | LANE CONFLICT: this request needs both the acrylamide and sulfur lanes at once. They do not compose -- the acrylamide network deliberately omits every sulfur step (acrylamide.OUT_OF_SCOPE), because composing them would s |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | maillard_path_holdout | Furfural | LANE CONFLICT: this request needs both the acrylamide and sulfur lanes at once. They do not compose -- the acrylamide network deliberately omits every sulfur step (acrylamide.OUT_OF_SCOPE), because composing them would s |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | 1-hexanol | UNREPRESENTED TARGETS: 1-hexanol -- The B6 lipid lane exists and forms the SIX products Frankel 1989 measured, but 1-hexanol is not one of them and NO aldehyde-reduction step is measured anywhere in the corpus -- in a th |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | 2-pentylfuran | UNREPRESENTED TARGETS: 2-pentylfuran -- The B6 lipid lane exists, but 2-pentylfuran is NOT in Frankel 1989's six-product slate and no branch fraction for the linoleate -> alkylfuran route is measured anywhere in the fit  |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | nonanal | UNREPRESENTED TARGETS: nonanal -- the lipid lane exists and nonanal is a species in it, but its ONLY parent is the OLEATE hydroperoxide pool and the oleate -> nonanal branch fraction is measured NOWHERE in the fit corpus |
| external_validation_liu_2023_ppi_offnote_baseline | external_matrix | nonanal | UNREPRESENTED TARGETS: nonanal -- the lipid lane exists and nonanal is a species in it, but its ONLY parent is the OLEATE hydroperoxide pool and the oleate -> nonanal branch fraction is measured NOWHERE in the fit corpus |

## Bundles kept off the scored panel

- pea_isolate_ribose_cysteine_100C_45min_Internal2026 (trust_loop): synthetic snapshot (_Internal2026): legacy-model output, not a measurement
- soy_isolate_ribose_cysteine_100C_45min_Internal2026 (trust_loop): synthetic snapshot (_Internal2026): legacy-model output, not a measurement

## Priors

| key | lane | kind | distribution | centre | sigma | band | sampled | reason |
|---|---|---|---|---|---|---|---|---|
| b1.k_glc_frag.log10_k_ref_100C | trunk | fitted_rate | fixed | -8.000 | - | - | no | unidentified in the fit report (log10_k_ref_stderr is null) |
| b1.k_glc_frag.ea_kj_mol | trunk | fitted_ea | fixed | 180.695 | - | - | no | unidentified in the fit report (ea_stderr_kj_mol is null) |
| b1.k_mgo_mel.log10_k_ref_100C | trunk | fitted_rate | normal_log10 | -1.643 | 0.114 | - | yes | identified in the fit report (stderr reported) |
| b1.k_mgo_mel.ea_kj_mol | trunk | fitted_ea | normal | 20.043 | 31.414 | [20.000, 260.000] | yes | identified in the fit report (stderr reported); draws are CLIPPED to the fit's own search bounds FITTED_EA_BOUNDS, which |
| b1.k_fa_frag.log10_k_ref_100C | trunk | fitted_rate | fixed | -7.460 | - | - | no | unidentified in the fit report (log10_k_ref_stderr is null) |
| b1.k_fa_frag.ea_kj_mol | trunk | fitted_ea | fixed | 20.531 | - | - | no | unidentified in the fit report (ea_stderr_kj_mol is null) |
| b1.k_aa_frag.log10_k_ref_100C | trunk | fitted_rate | normal_log10 | -1.928 | 0.082 | - | yes | identified in the fit report (stderr reported) |
| b1.k_aa_frag.ea_kj_mol | trunk | fitted_ea | normal | 20.000 | 20.865 | [20.000, 260.000] | yes | identified in the fit report (stderr reported); draws are CLIPPED to the fit's own search bounds FITTED_EA_BOUNDS, which |
| b3.k_int1_mel.log10_k_ref_160C | acrylamide | fitted_rate | fixed | -4.496 | - | - | no | unidentified in the fit report (ci95_halfwidth 1316.322310548673 above identified_threshold 1.0) |
| b3.k_acr_dp.log10_k_ref_160C | acrylamide | fitted_rate | normal_log10 | -0.898 | 0.227 | - | yes | identified in the fit report (ci95_halfwidth below the identified_threshold) |
| b3.k_gln_glc.log10_k_ref_160C | acrylamide | fitted_rate | fixed | -7.836 | - | - | no | unidentified in the fit report (ci95_halfwidth 2157390.4724324015 above identified_threshold 1.0) |
| b3.k_lys_glc.log10_k_ref_160C | acrylamide | fitted_rate | fixed | -2.596 | - | - | no | unidentified in the fit report (ci95_halfwidth 27.163628456550253 above identified_threshold 1.0) |
| b3.k_ala_glc.log10_k_ref_160C | acrylamide | fitted_rate | fixed | -6.872 | - | - | no | unidentified in the fit report (ci95_halfwidth 446550.9087512145 above identified_threshold 1.0) |
| b3.k_acr_gln.log10_k_ref_160C | acrylamide | fitted_rate | fixed | -1.895 | - | - | no | unidentified in the fit report (ci95_halfwidth 3.125748413772875 above identified_threshold 1.0) |
| b3.k_acr_lys.log10_k_ref_160C | acrylamide | fitted_rate | fixed | -2.492 | - | - | no | unidentified in the fit report (ci95_halfwidth 5.46047262651422 above identified_threshold 1.0) |
| b3.k_acr_ala.log10_k_ref_160C | acrylamide | fitted_rate | fixed | -8.825 | - | - | no | unidentified in the fit report (ci95_halfwidth 2305878.1472097607 above identified_threshold 1.0) |
| b3.Ea_int1_mel | acrylamide | fitted_ea | fixed | 260.000 | - | - | no | unidentified in the fit report (ci95_halfwidth 122057.37351828751 above identified_threshold 60.0) |
| b3.Ea_acr_dp | acrylamide | fitted_ea | normal | 136.109 | 24.543 | - | yes | identified in the fit report (ci95_halfwidth below the identified_threshold) |
| b3.Ea_competitor_sugar | acrylamide | fitted_ea | fixed | 20.000 | - | - | no | unidentified in the fit report (ci95_halfwidth 4383.491567802573 above identified_threshold 60.0) |
| b7.k_dpo_af.log10_k | trunk | fitted_rate | normal_log10 | -5.395 | 0.041 | - | yes | the B7 report carries no parameter stderr; its residual sigma_log10 is used as the log10-k stderr proxy, per the plan. k |
| b8.k_pent_dpo.log10_k_ref_145C | sulfur | fitted_rate | fixed | -0.545 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_pent_tdp.log10_k_ref_145C | sulfur | fitted_rate | fixed | -0.502 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_dpo_c2c3.log10_k_ref_145C | sulfur | fitted_rate | fixed | -3.016 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_arp_dpo.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -2.325 | 2.198 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_arp_tdp.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -1.790 | 0.422 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_dpo_nf.log10_k_ref_145C | sulfur | fitted_rate | fixed | 0.434 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_dpo_ptr.log10_k_ref_145C | sulfur | fitted_rate | fixed | -4.347 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_dpo_ddp.log10_k_ref_145C | sulfur | fitted_rate | fixed | 0.042 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_tdp_fur.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -3.005 | 0.317 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_ddp_mft.log10_k_ref_145C | sulfur | fitted_rate | fixed | -6.541 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_fur_fft.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -0.550 | 0.668 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_nf_mft.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -2.042 | 0.151 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_nf_mp3p.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -2.859 | 0.323 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_mgo_mp.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -3.891 | 0.279 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_ha_mp_mft.log10_k_ref_145C | sulfur | fitted_rate | fixed | -3.870 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_glc_ha.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -0.757 | 0.495 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_thi_hmp.log10_k_ref_145C | sulfur | fitted_rate | fixed | -2.588 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_thi_mesh.log10_k_ref_145C | sulfur | fitted_rate | fixed | -3.273 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_hmp_mft.log10_k_ref_145C | sulfur | fitted_rate | fixed | -2.609 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_hmp_mp2p.log10_k_ref_145C | sulfur | fitted_rate | fixed | -3.437 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_cys_actz.log10_k_ref_145C | sulfur | fitted_rate | fixed | -2.974 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_dimer_mft.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | 0.500 | 0.327 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_dimer_fft.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | 0.500 | 0.505 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_mmft.log10_k_ref_145C | sulfur | fitted_rate | fixed | -1.984 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_mft_decay.log10_k_ref_145C | sulfur | fitted_rate | fixed | 0.134 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_fft_decay.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -0.657 | 0.099 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_dimer_decay.log10_k_ref_145C | sulfur | fitted_rate | fixed | -9.563 | 5.34e+05 | [-10.000, 0.500] | no | unidentified_direction_in_laplace_covariance |
| b8.k_nf_decay.log10_k_ref_145C | sulfur | fitted_rate | fixed | -2.182 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_fur_decay.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | 0.436 | 0.542 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_h2s_loss.log10_k_ref_145C | sulfur | fitted_rate | fixed | -1.351 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_osone_decay.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -1.213 | 0.475 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_thiol_decay.log10_k_ref_145C | sulfur | fitted_rate | fixed | -9.483 | 4.06e+04 | [-10.000, 0.500] | no | unidentified_direction_in_laplace_covariance |
| b8.k_pent_caramel.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -2.597 | 0.573 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_pent_thermal.log10_k_ref_145C | sulfur | fitted_rate | fixed | -3.437 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_glc_fur.log10_k_ref_145C | sulfur | fitted_rate | fixed | -4.923 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_arp_tdp_th.log10_k_ref_145C | sulfur | fitted_rate | fixed | -4.354 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_arp_dpo_th.log10_k_ref_145C | sulfur | fitted_rate | fixed | -1.290 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_ddp_mft_hs.log10_k_ref_145C | sulfur | fitted_rate | fixed | -7.693 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_fur_fft_hs.log10_k_ref_145C | sulfur | fitted_rate | fixed | -0.00955 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_ttca_cys.log10_k_ref_145C | sulfur | fitted_rate | fixed | -0.579 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_ttca_deg.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -1.694 | 0.375 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_cys_thermal.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -2.052 | 0.253 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_thiolate_loss.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -1.834 | 0.614 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.lumped_formation_Ea_kJ_mol | sulfur | fitted_ea | fixed | 64.082 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.decay_Ea_kJ_mol.thiol_sink | sulfur | fitted_ea | fixed | 102.000 | 27.224 | [7.000, 102.000] | no | unidentified_direction_in_laplace_covariance |
| b8.decay_Ea_kJ_mol.carbonyl_sink | sulfur | fitted_ea | fixed | 166.675 | 63.584 | [20.000, 250.000] | no | unidentified_direction_in_laplace_covariance |
| b8.ph_drift.acid_yield_per_sink_event | sulfur | fitted_ph_drift | fixed | 3.76e-06 | 0.289 | [0.000, 1.000] | no | unidentified_direction_in_laplace_covariance |
| b8.ph_drift.arp_secondary_ammonium_pKa | sulfur | fitted_ph_drift | normal | 7.092 | 0.438 | [5.000, 11.000] | yes | laplace_covariance_at_b8_optimum |
| lipid.q10 | lipid | declared_band | uniform | 2.449 | - | [2.000, 3.000] | yes | declared corner band, sampled uniform over it |
| lipid.pea_protein_isolate.lipid_mass_fraction | lipid | declared_band | log_uniform | 0.025 | - | [0.010, 0.060] | yes | declared corner band, sampled log-uniform over it as ONE scale shared by every carrier in a draw (CoreDraw.lipid_fractio |
| lipid.pea_protein_isolate.peroxide_value_meq_per_kg | lipid | declared_band | log_uniform | 10.000 | - | [2.000, 40.000] | yes | declared corner band, sampled log-uniform over it as ONE scale shared by every carrier in a draw (CoreDraw.peroxide_scal |
| lipid.soy_protein_isolate.lipid_mass_fraction | lipid | declared_band | log_uniform | 0.020 | - | [0.008, 0.050] | yes | declared corner band, sampled log-uniform over it as ONE scale shared by every carrier in a draw (CoreDraw.lipid_fractio |
| lipid.soy_protein_isolate.peroxide_value_meq_per_kg | lipid | declared_band | log_uniform | 10.000 | - | [2.000, 40.000] | yes | declared corner band, sampled log-uniform over it as ONE scale shared by every carrier in a draw (CoreDraw.peroxide_scal |
| lipid.frankel_pure_hydroperoxide.lipid_mass_fraction | lipid | declared_band | fixed | 1.000 | - | [1.000, 1.000] | no | degenerate band (fed hydroperoxide: the fraction is the definition) |
| lipid.frankel_pure_hydroperoxide.peroxide_value_meq_per_kg | lipid | declared_band | fixed | 2e+03 | - | [2e+03, 2e+03] | no | degenerate band (fed hydroperoxide: PV is the definition) |
| furanic.partition_ea_offset_kj_mol | trunk | declared_band | uniform | 0.000 | - | [-50.000, 50.000] | yes | declared corner band on the furanone PARTITION barrier, sampled uniform |
| observable.air_water_partition_constant | observable | observable | log_uniform | 0.000 | - | [-0.500, 0.500] | yes | declared +/-0.5 dex band on K_aw, sampled uniform in log10 |
| observable.hs_spme_same_sample_dispersion | observable | observable | log_uniform_dispersion | 15.166 | - | [10.000, 23.000] | yes | measured 10-23x same-sample dispersion: D is drawn log-uniform over [10, 23] and the multiplier is D^u with u uniform on |
