# Core prediction uncertainty (Monte-Carlo envelope on the kinetic core)

n_samples = 200, seed = 0, CI level = 90 %.

* benchmarks with an envelope: **32** of 40 on the panel; matched rows **49**; refused rows 18
* mixed-population coverage: 10/49 (0.204)
* **honest literature coverage: 10/44 (0.227)**, median CI width 1.170 log10; 5 not evaluable; 0 fitted rows excluded
* out-of-sample literature coverage: 10/43 (5 not evaluable); rows the core fit read: {'hits': 0, 'total': 1, 'not_evaluable': 0}
* sampled priors 35, fixed 43; lanes with NO sampled fit uncertainty: none
* observable bands (K_aw, HS-SPME) applied by quantification family -- rows: headspace 2, extraction 37, undeclared 10; undeclared bundles: acrylamide_spi_extrusion_130C_ACSRef3, cys_ribose_140C_Hofmann1998, external_validation_li_2026_spi_wg_hme_control, external_validation_liu_2023_ppi_offnote_baseline, pea_isolate_40C_PratapSingh2021, pea_isolate_uht_140C_Trikusuma2019, resconi_2023_pbma_beef_identity_benchmark, soy_isolate_40C_PratapSingh2021, thiamine_cys_glucose_120C_Bolton1994

## Per panel

| panel | hits | total | rate | median width (log10) | not evaluable |
|---|---|---|---|---|---|
| external_matrix | 2 | 4 | 0.500 | 2.096 | 0 |
| maillard_path_holdout | 4 | 25 | 0.160 | 0.941 | 5 |
| trust_loop | 4 | 15 | 0.267 | 1.652 | 0 |

## Rows

| benchmark | panel | compound | unit | measured | point | p5 | p50 | p95 | inside | width | obs bands | lane | role |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| acrylamide_spi_extrusion_130C_ACSRef3 | trust_loop | acrylamide | ppb | 150.000 | 0.035 | 0.00726 | 0.038 | 0.234 | no | 1.508 | yes (undeclared) | acrylamide | predictive |
| cys_ribose_140C_Hofmann1998 | trust_loop | 2-methyl-3-furanthiol | ppb | 342.000 | 76.277 | 11.623 | 79.309 | 521.636 | yes | 1.652 | yes (undeclared) | sulfur | predictive |
| cys_ribose_140C_Hofmann1998 | trust_loop | 2-furfurylthiol | ppb | 200.000 | 57.694 | 4.357 | 63.675 | 1.03e+03 | yes | 2.373 | yes (undeclared) | sulfur | predictive |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | ppb | 32.000 | 1.67e+03 | 183.821 | 2.41e+03 | 1.07e+04 | no | 1.766 | no (extraction) | sulfur | predictive |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 25.000 | 2.98e-10 | 7.9e-11 | 2.59e-10 | 8.52e-10 | no | 1.032 | no (extraction) | sulfur | predictive |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | ppb | 28.000 | 2.3e+03 | 275.074 | 3.25e+03 | 1.41e+04 | no | 1.710 | no (extraction) | sulfur | predictive |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 19.000 | 1.37e-09 | 3.48e-10 | 1.16e-09 | 3.92e-09 | no | 1.051 | no (extraction) | sulfur | predictive |
| hofmann1998_norfuraneol_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 1.02e+03 | 2.36e+03 | 1.43e+03 | 2.29e+03 | 2.81e+03 | no | 0.293 | no (extraction) | sulfur | predictive [in core fit] |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | ppb | 121.000 | 834.391 | 55.756 | 798.227 | 7.92e+03 | yes | 2.153 | no (extraction) | sulfur | predictive |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 198.000 | 731.065 | 348.934 | 682.743 | 1.16e+03 | no | 0.521 | no (extraction) | sulfur | predictive |
| pea_isolate_40C_PratapSingh2021 | trust_loop | hexanal | ppb | 1.14e+03 | 0.339 | 0.036 | 0.369 | 4.381 | no | 2.089 | yes (undeclared) | lipid | predictive |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | hexanal | ppb | 782.000 | 22.839 | 1.223 | 25.195 | 459.646 | no | 2.575 | yes (undeclared) | lipid | predictive |
| resconi_2023_pbma_beef_identity_benchmark | trust_loop | furfural | ppb | 715.220 | 7.665 | 1.685 | 11.247 | 75.170 | no | 1.650 | yes (undeclared) | sulfur | predictive |
| soy_isolate_40C_PratapSingh2021 | trust_loop | hexanal | ppb | 1.62e+03 | 0.267 | 0.028 | 0.290 | 3.449 | no | 2.089 | yes (undeclared) | lipid | predictive |
| thiamine_cys_glucose_120C_Bolton1994 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 13.000 | 17.314 | 3.560 | 18.602 | 114.536 | yes | 1.508 | yes (undeclared) | sulfur | predictive |
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
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 7.000 | 837.456 | 156.497 | 1.29e+03 | 5.47e+03 | no | 1.544 | no (extraction) | sulfur | predictive |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 3.000 | 1.4e-09 | 3.66e-10 | 1.23e-09 | 3.98e-09 | no | 1.037 | no (extraction) | sulfur | predictive |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 6.000 | 3.027 | 0.310 | 2.776 | 6.705 | yes | 1.335 | no (extraction) | sulfur | predictive |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 4.000 | 1.23e-09 | 1.53e-10 | 1e-09 | 2.97e-09 | no | 1.288 | no (extraction) | sulfur | predictive |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 229.000 | 1.14e+03 | 82.609 | 1.06e+03 | 9.94e+03 | yes | 2.080 | no (extraction) | sulfur | predictive (shared: hofmann_ribose_pH3_FFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 553.000 | 23.470 | 10.725 | 21.966 | 39.503 | no | 0.566 | no (extraction) | sulfur | predictive (shared: hofmann_ribose_pH3_MFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 12.000 | 0.090 | 0.00373 | 0.085 | 0.760 | no | 2.309 | no (extraction) | sulfur | predictive (shared: hofmann_ribose_pH7_FFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 25.000 | 0.178 | 0.025 | 0.149 | 0.256 | no | 1.020 | no (extraction) | sulfur | predictive (shared: hofmann_ribose_pH7_MFT) |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 96.000 | 192.701 | 21.700 | 211.607 | 1.61e+03 | yes | 1.871 | no (extraction) | sulfur | predictive |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 143.000 | 548.334 | 288.962 | 514.430 | 761.050 | no | 0.421 | no (extraction) | sulfur | predictive |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 6.880 | 64.116 | 14.121 | 51.237 | 101.783 | no | 0.858 | no (extraction) | sulfur | predictive |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 1.280 | 615.517 | 228.932 | 697.859 | 1.75e+03 | no | 0.883 | no (extraction) | sulfur | predictive |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 3.290 | 96.974 | 30.077 | 81.614 | 143.827 | no | 0.680 | no (extraction) | sulfur | predictive |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 1.460 | 541.772 | 200.517 | 682.398 | 2.07e+03 | no | 1.014 | no (extraction) | sulfur | predictive |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 2.400 | 134.711 | 57.742 | 122.305 | 212.969 | no | 0.567 | no (extraction) | sulfur | predictive |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 1.680 | 367.708 | 134.388 | 484.971 | 1.87e+03 | no | 1.142 | no (extraction) | sulfur | predictive |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 1.710 | 170.934 | 85.775 | 165.530 | 294.625 | no | 0.536 | no (extraction) | sulfur | predictive |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 1.620 | 213.745 | 79.207 | 294.094 | 1.25e+03 | no | 1.197 | no (extraction) | sulfur | predictive |
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
| b8.k_arp_dpo.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -2.310 | 1.805 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_arp_tdp.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -1.727 | 0.607 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_dpo_nf.log10_k_ref_145C | sulfur | fitted_rate | fixed | 0.434 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_dpo_ptr.log10_k_ref_145C | sulfur | fitted_rate | fixed | -4.347 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_dpo_ddp.log10_k_ref_145C | sulfur | fitted_rate | fixed | 0.042 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_tdp_fur.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -3.037 | 0.340 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_ddp_mft.log10_k_ref_145C | sulfur | fitted_rate | fixed | -6.541 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_fur_fft.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -0.631 | 1.028 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_nf_mft.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -2.025 | 0.165 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_nf_mp3p.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -2.868 | 0.348 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_mgo_mp.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -3.981 | 0.320 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_ha_mp_mft.log10_k_ref_145C | sulfur | fitted_rate | fixed | -3.870 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_glc_ha.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -8.589 | 1.37e-11 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_thi_hmp.log10_k_ref_145C | sulfur | fitted_rate | fixed | -2.588 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_thi_mesh.log10_k_ref_145C | sulfur | fitted_rate | fixed | -3.273 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_hmp_mft.log10_k_ref_145C | sulfur | fitted_rate | fixed | -2.609 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_hmp_mp2p.log10_k_ref_145C | sulfur | fitted_rate | fixed | -3.437 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_cys_actz.log10_k_ref_145C | sulfur | fitted_rate | fixed | -2.974 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_dimer_mft.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | 0.500 | 0.354 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_dimer_fft.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | 0.500 | 0.549 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_mmft.log10_k_ref_145C | sulfur | fitted_rate | fixed | -1.984 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_mft_decay.log10_k_ref_145C | sulfur | fitted_rate | fixed | 0.134 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_fft_decay.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -0.651 | 0.107 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_dimer_decay.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -9.710 | 1.2e-05 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_nf_decay.log10_k_ref_145C | sulfur | fitted_rate | fixed | -2.182 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_fur_decay.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | 0.470 | 0.437 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_h2s_loss.log10_k_ref_145C | sulfur | fitted_rate | fixed | -1.351 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_osone_decay.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -1.187 | 0.449 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_thiol_decay.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -9.166 | 9.89e-05 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_pent_caramel.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -2.535 | 0.653 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_pent_thermal.log10_k_ref_145C | sulfur | fitted_rate | fixed | -3.437 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_glc_fur.log10_k_ref_145C | sulfur | fitted_rate | fixed | -4.923 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_arp_tdp_th.log10_k_ref_145C | sulfur | fitted_rate | fixed | -4.354 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_arp_dpo_th.log10_k_ref_145C | sulfur | fitted_rate | fixed | -1.290 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_ddp_mft_hs.log10_k_ref_145C | sulfur | fitted_rate | fixed | -7.693 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_fur_fft_hs.log10_k_ref_145C | sulfur | fitted_rate | fixed | -0.00955 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_ttca_cys.log10_k_ref_145C | sulfur | fitted_rate | fixed | -0.579 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.k_ttca_deg.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -1.698 | 0.408 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_cys_thermal.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -2.066 | 0.275 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_thiolate_loss.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -1.832 | 0.661 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.lumped_formation_Ea_kJ_mol | sulfur | fitted_ea | fixed | 64.082 | - | - | no | frozen in the B8 fit (not a free coordinate) |
| b8.decay_Ea_kJ_mol.thiol_sink | sulfur | fitted_ea | fixed | 102.000 | 29.240 | [7.000, 102.000] | no | unidentified_direction_in_laplace_covariance |
| b8.decay_Ea_kJ_mol.carbonyl_sink | sulfur | fitted_ea | fixed | 174.922 | 60.566 | [20.000, 250.000] | no | unidentified_direction_in_laplace_covariance |
| b8.ph_drift.acid_yield_per_sink_event | sulfur | fitted_ph_drift | fixed | 0.000359 | 0.340 | [0.000, 1.000] | no | unidentified_direction_in_laplace_covariance |
| b8.ph_drift.arp_secondary_ammonium_pKa | sulfur | fitted_ph_drift | normal | 7.062 | 0.482 | [5.000, 11.000] | yes | laplace_covariance_at_b8_optimum |
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
