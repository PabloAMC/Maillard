# Core prediction uncertainty (Monte-Carlo envelope on the kinetic core)

n_samples = 200, seed = 0, CI level = 90 %.

* benchmarks with an envelope: **27** of 37 on the panel; matched rows **39**; refused rows 25
* mixed-population coverage: 5/39 (0.128)
* **honest literature coverage: 5/33 (0.152)**, median CI width 1.369 log10; 6 not evaluable; 0 fitted rows excluded
* out-of-sample literature coverage: 5/32 (6 not evaluable); rows the core fit read: {'hits': 0, 'total': 1, 'not_evaluable': 0}
* sampled priors 41, fixed 37; lanes with NO sampled fit uncertainty: none
* observable bands (K_aw, HS-SPME) applied by quantification family -- rows: headspace 8, extraction 31, undeclared 0

## Per panel

| panel | hits | total | rate | median width (log10) | not evaluable |
|---|---|---|---|---|---|
| external_matrix | 1 | 4 | 0.250 | 2.049 | 0 |
| maillard_path_holdout | 3 | 21 | 0.143 | 0.905 | 5 |
| trust_loop | 1 | 8 | 0.125 | 1.858 | 1 |

## Rows

| benchmark | panel | compound | unit | measured | point | p5 | p50 | p95 | inside | width | obs bands | lane | role |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| acrylamide_spi_extrusion_130C_ACSRef3 | trust_loop | acrylamide | ppb | 150.000 | 0.035 | 0.035 | 0.035 | 0.035 | no | 0.00106 | no (extraction) | acrylamide | predictive |
| hofmann1998_norfuraneol_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 1.02e+03 | 2.35e+03 | 1.44e+03 | 2.25e+03 | 2.73e+03 | no | 0.279 | no (extraction) | sulfur | predictive [in core fit] |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | ppb | 121.000 | 830.643 | 45.420 | 1.01e+03 | 8.09e+03 | yes | 2.250 | no (extraction) | sulfur | predictive |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 198.000 | 730.638 | 383.849 | 654.728 | 1.06e+03 | no | 0.443 | no (extraction) | sulfur | predictive |
| pea_isolate_40C_PratapSingh2021 | trust_loop | hexanal | ppb | 1.14e+03 | 0.339 | 0.030 | 0.289 | 3.410 | no | 2.051 | yes (headspace) | lipid | predictive |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | hexanal | ppb | 782.000 | 22.839 | 0.717 | 19.721 | 508.110 | no | 2.851 | yes (headspace) | lipid | predictive |
| resconi_2023_pbma_beef_identity_benchmark | trust_loop | furfural | ppb | 715.220 | 7.665 | 1.365 | 9.930 | 63.129 | no | 1.665 | yes (headspace) | sulfur | predictive |
| soy_isolate_40C_PratapSingh2021 | trust_loop | hexanal | ppb | 1.62e+03 | 0.267 | 0.024 | 0.227 | 2.684 | no | 2.051 | yes (headspace) | lipid | predictive |
| thiamine_cys_glucose_120C_Bolton1994 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 11.700 | 235.973 | 48.710 | 103.791 | 219.453 | no | 0.654 | no (extraction) | sulfur | predictive |
| mp_holdout_fructose_asparagine_180C_Lin2022 | maillard_path_holdout | Acrylamide | ppb | 1.86e+03 | 225.440 | 89.867 | 229.650 | 712.828 | no | 0.899 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_fructose_asparagine_180C_Lin2022 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | ppb | 1.23e+04 | 1.94e+03 | 1.94e+03 | 1.94e+03 | 1.94e+03 | no | 3.36e-05 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | maillard_path_holdout | DMHF | ppb | 1.15e+03 | 21.841 | 7.780 | 24.604 | 62.585 | no | 0.905 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | ppb | 5.73e+04 | 2.8e+04 | 2.8e+04 | 2.8e+04 | 2.8e+04 | no | 0.000205 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | maillard_path_holdout | DMHF | ppb | 5.89e+03 | 21.841 | 7.780 | 24.604 | 62.585 | no | 0.905 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | ppb | 1.01e+05 | 2.8e+04 | 2.8e+04 | 2.8e+04 | 2.8e+04 | no | 0.000205 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_asparagine_180C_10min_Chang2021 | maillard_path_holdout | Acrylamide | ppb | 28.000 | 6.75e+03 | 2.69e+03 | 6.84e+03 | 1.83e+04 | no | 0.833 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_asparagine_180C_30min_Chang2021 | maillard_path_holdout | Acrylamide | ppb | 1.46e+03 | 4.03e+03 | 1.61e+03 | 4.11e+03 | 1.27e+04 | no | 0.899 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | maillard_path_holdout | Acrylamide | ppb | 832.000 | 4.03e+03 | 1.61e+03 | 4.11e+03 | 1.27e+04 | no | 0.899 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | ppb | 7e+03 | 2.11e+03 | 2.11e+03 | 2.11e+03 | 2.11e+03 | no | 0.000586 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_asparagine_180C_Ye2024 | maillard_path_holdout | Acrylamide | umol_per_mol_limiting_precursor | 140.580 | 7.04e+03 | 2.81e+03 | 7.17e+03 | 2.23e+04 | no | 0.899 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | ppb | 1.74e+04 | 1.46e+03 | 1.46e+03 | 1.46e+03 | 1.46e+03 | no | 7.04e-07 | no (extraction) | trunk | external_holdout |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 229.000 | 1.13e+03 | 64.571 | 1.38e+03 | 1.04e+04 | yes | 2.208 | no (extraction) | sulfur | external_holdout (shared: hofmann_ribose_pH3_FFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 553.000 | 23.470 | 11.858 | 21.144 | 35.332 | no | 0.474 | no (extraction) | sulfur | external_holdout (shared: hofmann_ribose_pH3_MFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 12.000 | 0.090 | 0.00478 | 0.085 | 0.597 | no | 2.097 | no (extraction) | sulfur | external_holdout (shared: hofmann_ribose_pH7_FFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 25.000 | 0.178 | 0.032 | 0.153 | 0.255 | no | 0.902 | no (extraction) | sulfur | external_holdout (shared: hofmann_ribose_pH7_MFT) |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 96.000 | 191.058 | 28.100 | 275.662 | 1.61e+03 | yes | 1.757 | no (extraction) | sulfur | external_holdout |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 143.000 | 546.398 | 312.866 | 484.552 | 708.893 | no | 0.355 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 6.880 | 64.086 | 1.920 | 8.886 | 45.547 | yes | 1.375 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 1.280 | 610.255 | 18.475 | 179.284 | 869.569 | no | 1.673 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 3.290 | 96.892 | 6.508 | 21.702 | 73.254 | no | 1.051 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 1.460 | 535.496 | 28.900 | 238.136 | 1.05e+03 | no | 1.561 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 2.400 | 134.513 | 19.271 | 48.146 | 116.080 | no | 0.780 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 1.680 | 363.445 | 37.181 | 252.246 | 987.961 | no | 1.424 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 1.710 | 170.534 | 47.418 | 96.612 | 181.001 | no | 0.582 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 1.620 | 211.740 | 37.318 | 195.546 | 872.925 | no | 1.369 | no (extraction) | sulfur | external_holdout |
| external_validation_bi_2020_raw_pea_hexanal | external_matrix | hexanal | ppb | 1.26e+03 | 0.339 | 0.030 | 0.289 | 3.410 | no | 2.051 | yes (headspace) | lipid | external_holdout |
| external_validation_bi_2020_roasted_pea_hexanal | external_matrix | hexanal | ppb | 324.000 | 88.598 | 6.700 | 80.020 | 746.145 | yes | 2.047 | yes (headspace) | lipid | external_holdout |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | hexanal | ppb | 605.600 | 69.695 | 5.170 | 55.355 | 570.820 | no | 2.043 | yes (headspace) | lipid | external_holdout |
| external_validation_liu_2023_ppi_offnote_baseline | external_matrix | hexanal | ppb | 1.13e+04 | 0.339 | 0.030 | 0.289 | 3.410 | no | 2.051 | yes (headspace) | lipid | external_holdout |

## Refused rows

| benchmark | panel | compound | reason |
|---|---|---|---|
| cml_cel_commercial_pbma_Foods2023 | trust_loop | Nε-(Carboxymethyl)lysine (CML) | UNREPRESENTED TARGETS: Nε-(Carboxymethyl)lysine (CML) -- not a species in any core lane, and not on the named unrepresented-compound list either: the engine has no vocabulary entry for it. |
| cml_cel_commercial_pbma_Foods2023 | trust_loop | Nε-(Carboxyethyl)lysine (CEL) | UNREPRESENTED TARGETS: Nε-(Carboxyethyl)lysine (CEL) -- not a species in any core lane, and not on the named unrepresented-compound list either: the engine has no vocabulary entry for it. |
| hofmann1998_c2c3_recombination_145C_20min_pH3 | trust_loop | 2-Methyl-3-furanthiol (MFT) | UNMAPPED PRECURSORS 'Hydroxyacetaldehyde', 'Mercapto-2-propanone': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge. |
| hofmann1998_c2c3_recombination_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | UNMAPPED PRECURSORS 'Hydroxyacetaldehyde', 'Mercapto-2-propanone': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge. |
| hofmann1998_c2c3_recombination_145C_20min_pH7 | trust_loop | 2-Methyl-3-furanthiol (MFT) | UNMAPPED PRECURSORS 'Hydroxyacetaldehyde', 'Mercapto-2-propanone': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge. |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_furan2aldehyde_h2s_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | UNMAPPED PRECURSORS 'Furan-2-aldehyde', 'Hydrogen sulfide': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge. |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_norfuraneol_h2s_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | UNMAPPED PRECURSORS 'Hydrogen sulfide': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge. |
| pea_isolate_40C_PratapSingh2021 | trust_loop | 2-pentylfuran | UNREPRESENTED TARGETS: 2-pentylfuran -- The lipid lane exists, but 2-pentylfuran is NOT in Frankel 1989's six-product slate and no branch fraction for the linoleate -> alkylfuran route is measured anywhere in the fit cor |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | 2-pentylfuran | UNREPRESENTED TARGETS: 2-pentylfuran -- The lipid lane exists, but 2-pentylfuran is NOT in Frankel 1989's six-product slate and no branch fraction for the linoleate -> alkylfuran route is measured anywhere in the fit cor |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | nonanal | UNREPRESENTED TARGETS: nonanal -- the lipid lane exists and nonanal is a species in it, but its ONLY parent is the OLEATE hydroperoxide pool and the oleate -> nonanal branch fraction is measured NOWHERE in the fit corpus |
| soy_isolate_40C_PratapSingh2021 | trust_loop | 2-pentylfuran | UNREPRESENTED TARGETS: 2-pentylfuran -- The lipid lane exists, but 2-pentylfuran is NOT in Frankel 1989's six-product slate and no branch fraction for the linoleate -> alkylfuran route is measured anywhere in the fit cor |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | maillard_path_holdout | Furfural | LANE CONFLICT: this request needs both the acrylamide and sulfur lanes at once. They do not compose -- the acrylamide network deliberately omits every sulfur step (acrylamide.OUT_OF_SCOPE), because composing them would s |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | maillard_path_holdout | Furfural | LANE CONFLICT: this request needs both the acrylamide and sulfur lanes at once. They do not compose -- the acrylamide network deliberately omits every sulfur step (acrylamide.OUT_OF_SCOPE), because composing them would s |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | 1-hexanol | UNREPRESENTED TARGETS: 1-hexanol -- The lipid lane exists and forms the SIX products Frankel 1989 measured, but 1-hexanol is not one of them and NO aldehyde-reduction step is measured anywhere in the corpus -- in a therm |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | 2-pentylfuran | UNREPRESENTED TARGETS: 2-pentylfuran -- The lipid lane exists, but 2-pentylfuran is NOT in Frankel 1989's six-product slate and no branch fraction for the linoleate -> alkylfuran route is measured anywhere in the fit cor |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | nonanal | UNREPRESENTED TARGETS: nonanal -- the lipid lane exists and nonanal is a species in it, but its ONLY parent is the OLEATE hydroperoxide pool and the oleate -> nonanal branch fraction is measured NOWHERE in the fit corpus |
| external_validation_liu_2023_ppi_offnote_baseline | external_matrix | nonanal | UNREPRESENTED TARGETS: nonanal -- the lipid lane exists and nonanal is a species in it, but its ONLY parent is the OLEATE hydroperoxide pool and the oleate -> nonanal branch fraction is measured NOWHERE in the fit corpus |

## Bundles kept off the scored panel

- pea_isolate_ribose_cysteine_100C_45min_Internal2026 (trust_loop): synthetic snapshot (_Internal2026): legacy-model output, not a measurement
- soy_isolate_ribose_cysteine_100C_45min_Internal2026 (trust_loop): synthetic snapshot (_Internal2026): legacy-model output, not a measurement

## Priors

| key | lane | kind | distribution | centre | sigma | band | sampled | reason |
|---|---|---|---|---|---|---|---|---|
| b1.k_glc_frag.log10_k_ref_100C | trunk | fitted_rate | fixed | -8.000 | - | - | no | unidentified_in_the_fit: no band is declared, still fixed (a recorded gap) |
| b1.k_glc_frag.ea_kj_mol | trunk | fitted_ea | uniform_band | 180.695 | - | [137.832, 223.559] | yes | unidentified_in_the_fit: declared band capped by the prefactor prior |
| b1.k_mgo_mel.log10_k_ref_100C | trunk | fitted_rate | normal_log10 | -1.643 | 0.114 | - | yes | identified in the fit report (stderr reported) |
| b1.k_mgo_mel.ea_kj_mol | trunk | fitted_ea | normal | 20.043 | 31.414 | [20.000, 260.000] | yes | identified in the fit report (stderr reported); draws are CLIPPED to the fit's own search bounds FITTED_EA_BOUNDS, which |
| b1.k_fa_frag.log10_k_ref_100C | trunk | fitted_rate | fixed | -7.460 | - | - | no | unidentified_in_the_fit: no band is declared, still fixed (a recorded gap) |
| b1.k_fa_frag.ea_kj_mol | trunk | fitted_ea | uniform_band | 20.531 | - | [20.000, 63.394] | yes | unidentified_in_the_fit: declared band capped by the prefactor prior |
| b1.k_aa_frag.log10_k_ref_100C | trunk | fitted_rate | normal_log10 | -1.928 | 0.082 | - | yes | identified in the fit report (stderr reported) |
| b1.k_aa_frag.ea_kj_mol | trunk | fitted_ea | normal | 20.000 | 20.865 | [20.000, 260.000] | yes | identified in the fit report (stderr reported); draws are CLIPPED to the fit's own search bounds FITTED_EA_BOUNDS, which |
| b3.k_int1_mel.log10_k_ref_160C | acrylamide | fitted_rate | fixed | -4.496 | - | - | no | unidentified_in_the_fit: no band is declared, still fixed (a recorded gap) |
| b3.k_acr_dp.log10_k_ref_160C | acrylamide | fitted_rate | normal_log10 | -0.898 | 0.227 | - | yes | identified in the fit report (ci95_halfwidth below the identified_threshold) |
| b3.k_gln_glc.log10_k_ref_160C | acrylamide | fitted_rate | fixed | -7.836 | - | - | no | unidentified_in_the_fit: no band is declared, still fixed (a recorded gap) |
| b3.k_lys_glc.log10_k_ref_160C | acrylamide | fitted_rate | fixed | -2.596 | - | - | no | unidentified_in_the_fit: no band is declared, still fixed (a recorded gap) |
| b3.k_ala_glc.log10_k_ref_160C | acrylamide | fitted_rate | fixed | -6.872 | - | - | no | unidentified_in_the_fit: no band is declared, still fixed (a recorded gap) |
| b3.k_acr_gln.log10_k_ref_160C | acrylamide | fitted_rate | fixed | -1.895 | - | - | no | unidentified_in_the_fit: no band is declared, still fixed (a recorded gap) |
| b3.k_acr_lys.log10_k_ref_160C | acrylamide | fitted_rate | fixed | -2.492 | - | - | no | unidentified_in_the_fit: no band is declared, still fixed (a recorded gap) |
| b3.k_acr_ala.log10_k_ref_160C | acrylamide | fitted_rate | fixed | -8.825 | - | - | no | unidentified_in_the_fit: no band is declared, still fixed (a recorded gap) |
| b3.Ea_int1_mel | acrylamide | fitted_ea | uniform_band | 260.000 | - | [210.244, 260.000] | yes | unidentified_in_the_fit: declared band capped by the prefactor prior |
| b3.Ea_acr_dp | acrylamide | fitted_ea | normal | 136.109 | 24.543 | - | yes | identified in the fit report (ci95_halfwidth below the identified_threshold) |
| b3.Ea_competitor_sugar | acrylamide | fitted_ea | uniform_band | 20.000 | - | [20.000, 69.755] | yes | unidentified_in_the_fit: declared band capped by the prefactor prior |
| b7.k_dpo_af.log10_k | trunk | fitted_rate | normal_log10 | -5.395 | 0.041 | - | yes | the furanic-channel fit report carries no parameter stderr; its residual sigma_log10 is used as the log10-k stderr proxy |
| b8.k_pent_dpo.log10_k_ref_145C | sulfur | fitted_rate | fixed | -0.545 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_pent_tdp.log10_k_ref_145C | sulfur | fitted_rate | fixed | -0.502 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_dpo_c2c3.log10_k_ref_145C | sulfur | fitted_rate | fixed | -3.016 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_arp_dpo.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -2.310 | 1.805 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_arp_tdp.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -1.727 | 0.607 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_dpo_nf.log10_k_ref_145C | sulfur | fitted_rate | fixed | 0.434 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_dpo_ptr.log10_k_ref_145C | sulfur | fitted_rate | fixed | -4.347 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_dpo_ddp.log10_k_ref_145C | sulfur | fitted_rate | fixed | 0.042 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_tdp_fur.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -3.037 | 0.340 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_ddp_mft.log10_k_ref_145C | sulfur | fitted_rate | fixed | -6.541 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_fur_fft.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -0.631 | 1.028 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_nf_mft.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -2.025 | 0.165 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_nf_mp3p.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -2.868 | 0.348 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_mgo_mp.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -3.981 | 0.320 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_ha_mp_mft.log10_k_ref_145C | sulfur | fitted_rate | fixed | -3.870 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_glc_ha.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -8.589 | 1.37e-11 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_thi_hmp.log10_k_ref_145C | sulfur | fitted_rate | fixed | -2.588 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_thi_mesh.log10_k_ref_145C | sulfur | fitted_rate | fixed | -3.273 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_hmp_mft.log10_k_ref_145C | sulfur | fitted_rate | fixed | -2.609 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_hmp_mp2p.log10_k_ref_145C | sulfur | fitted_rate | fixed | -3.437 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_cys_actz.log10_k_ref_145C | sulfur | fitted_rate | fixed | -2.974 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_dimer_mft.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | 0.500 | 0.354 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_dimer_fft.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | 0.500 | 0.549 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_mmft.log10_k_ref_145C | sulfur | fitted_rate | fixed | -1.984 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_mft_decay.log10_k_ref_145C | sulfur | fitted_rate | fixed | 0.134 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_fft_decay.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -0.651 | 0.107 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_dimer_decay.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -9.710 | 1.2e-05 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_nf_decay.log10_k_ref_145C | sulfur | fitted_rate | fixed | -2.182 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_fur_decay.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | 0.470 | 0.437 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_h2s_loss.log10_k_ref_145C | sulfur | fitted_rate | fixed | -1.351 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_osone_decay.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -1.187 | 0.449 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_thiol_decay.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -9.166 | 9.89e-05 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_pent_caramel.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -2.535 | 0.653 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_pent_thermal.log10_k_ref_145C | sulfur | fitted_rate | fixed | -3.437 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_glc_fur.log10_k_ref_145C | sulfur | fitted_rate | fixed | -4.923 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_arp_tdp_th.log10_k_ref_145C | sulfur | fitted_rate | fixed | -4.354 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_arp_dpo_th.log10_k_ref_145C | sulfur | fitted_rate | fixed | -1.290 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_ddp_mft_hs.log10_k_ref_145C | sulfur | fitted_rate | fixed | -7.693 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_fur_fft_hs.log10_k_ref_145C | sulfur | fitted_rate | fixed | -0.00955 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_ttca_cys.log10_k_ref_145C | sulfur | fitted_rate | fixed | -0.579 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_ttca_deg.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -1.698 | 0.408 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_cys_thermal.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -2.066 | 0.275 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_thiolate_loss.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -1.832 | 0.661 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.lumped_formation_Ea_kJ_mol | sulfur | fitted_ea | fixed | 64.082 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.decay_Ea_kJ_mol.thiol_sink | sulfur | fitted_ea | uniform_band | 102.000 | 29.240 | [7.000, 102.000] | yes | unidentified_in_the_fit: sampled across its declared band |
| b8.decay_Ea_kJ_mol.carbonyl_sink | sulfur | fitted_ea | uniform_band | 174.922 | 60.566 | [126.889, 222.954] | yes | unidentified_in_the_fit: declared band capped by the prefactor prior |
| b8.ph_drift.acid_yield_per_sink_event | sulfur | fitted_ph_drift | fixed | 0.000359 | 0.340 | [0.000, 1.000] | no | unidentified_in_the_fit: its bound is DEFINITIONAL (the quantity is a fraction), not a measurement of where the value li |
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
