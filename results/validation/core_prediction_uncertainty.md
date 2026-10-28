# Core prediction uncertainty (Monte-Carlo envelope on the kinetic core)

n_samples = 200, seed = 0, CI level = 90 %.

* benchmarks with an envelope: **33** of 42 on the panel; matched rows **66**; refused rows 27
* mixed-population coverage: 36/66 (0.545)
* **honest literature coverage: 36/65 (0.554)**, median CI width 1.869 log10; 1 not evaluable; 0 fitted rows excluded
* out-of-sample literature coverage: 33/59 (1 not evaluable); rows the core fit read: {'hits': 3, 'total': 6, 'not_evaluable': 0}
* B55 single-temperature Ea prior: on (31 rows); B56 structural discrepancy: off
* sampled priors 97, fixed 45; lanes with NO sampled fit uncertainty: none
* observable bands (K_aw, HS-SPME) applied by quantification family -- rows: headspace 25, extraction 41, undeclared 0

## Per panel

| panel | hits | total | rate | median width (log10) | not evaluable |
|---|---|---|---|---|---|
| external_matrix | 3 | 4 | 0.750 | 2.186 | 0 |
| maillard_path_holdout | 23 | 46 | 0.500 | 2.006 | 1 |
| trust_loop | 10 | 15 | 0.667 | 1.297 | 0 |

## Rows

| benchmark | panel | compound | unit | measured | point | p5 | p50 | p95 | inside | width | obs bands | lane | role |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| acrylamide_spi_extrusion_130C_ACSRef3 | trust_loop | acrylamide | ppb | 150.000 | 0.016 | 0.00741 | 0.016 | 0.025 | no | 0.536 | no (extraction) | acrylamide | predictive |
| hofmann1998_c2c3_recombination_145C_20min_pH3 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 310.000 | 4.03e+03 | 3.05e+03 | 4.03e+03 | 4.07e+03 | no | 0.125 | no (extraction) | sulfur | predictive [in core fit] |
| hofmann1998_c2c3_recombination_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 5.36e+03 | 2.92e+03 | 2.11e+03 | 2.96e+03 | 3.45e+03 | no | 0.213 | no (extraction) | sulfur | predictive [in core fit] |
| hofmann1998_c2c3_recombination_145C_20min_pH7 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 6.23e+03 | 1.15e+03 | 537.211 | 1.18e+03 | 1.74e+03 | no | 0.510 | no (extraction) | sulfur | predictive [in core fit] |
| hofmann1998_furan2aldehyde_h2s_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | ppb | 1.1e+04 | 1.7e+04 | 9.261 | 7.07e+03 | 5.48e+04 | yes | 3.772 | no (extraction) | sulfur | predictive [in core fit] |
| hofmann1998_norfuraneol_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 1.02e+03 | 1.07e+03 | 471.481 | 1.07e+03 | 1.87e+03 | yes | 0.598 | no (extraction) | sulfur | predictive [in core fit] |
| hofmann1998_norfuraneol_h2s_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 4.22e+03 | 8.27e+03 | 3.73e+03 | 7.07e+03 | 9.07e+03 | yes | 0.386 | no (extraction) | sulfur | predictive [in core fit] |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | ppb | 121.000 | 346.617 | 16.447 | 177.976 | 1.22e+03 | yes | 1.869 | no (extraction) | sulfur | predictive |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 198.000 | 93.633 | 37.127 | 95.921 | 219.557 | yes | 0.772 | no (extraction) | sulfur | predictive |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | hexanal | ppb | 782.000 | 353.858 | 67.930 | 363.927 | 2.25e+03 | yes | 1.521 | yes (headspace) | lipid | predictive |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | 2-pentylfuran | ppb | 163.000 | 64.442 | 12.355 | 69.268 | 414.631 | yes | 1.526 | yes (headspace) | lipid | predictive |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | nonanal | ppb | 24.000 | 15.606 | 2.463 | 17.466 | 166.040 | yes | 1.829 | yes (headspace) | lipid | predictive |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | (E,E)-2,4-decadienal | ppb | 46.900 | 75.702 | 3.163 | 69.137 | 1.24e+03 | yes | 2.593 | yes (headspace) | lipid | predictive |
| resconi_2023_pbma_beef_identity_benchmark | trust_loop | furfural | ppb | 715.220 | 46.453 | 1.67e-08 | 24.152 | 1.07e+04 | yes | 11.809 | yes (headspace) | sulfur | predictive |
| thiamine_cys_glucose_120C_Bolton1994 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 11.700 | 232.916 | 19.570 | 82.882 | 388.216 | no | 1.297 | no (extraction) | sulfur | predictive |
| mp_holdout_arp_glutathione_120C_120min_pH75_Zhou2023c | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 8.100 | 1.31e+03 | 18.899 | 322.541 | 3.84e+03 | no | 2.308 | yes (headspace) | sulfur | external_holdout |
| mp_holdout_arp_glutathione_120C_120min_pH75_Zhou2023c | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 116.000 | 129.540 | 0.019 | 5.193 | 622.566 | yes | 4.517 | yes (headspace) | sulfur | external_holdout |
| mp_holdout_arp_glutathione_120C_120min_pH75_Zhou2023c | maillard_path_holdout | Furfural | ppb | 716.000 | 24.138 | 0.011 | 4.736 | 386.463 | no | 4.562 | yes (headspace) | sulfur | external_holdout |
| mp_holdout_arp_glutathione_120C_30min_pH75_Zhou2023c | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 59.100 | 844.862 | 12.449 | 188.362 | 1.71e+03 | yes | 2.137 | yes (headspace) | sulfur | external_holdout |
| mp_holdout_arp_glutathione_120C_30min_pH75_Zhou2023c | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 124.000 | 25.671 | 0.024 | 1.927 | 185.895 | yes | 3.888 | yes (headspace) | sulfur | external_holdout |
| mp_holdout_arp_glutathione_120C_30min_pH75_Zhou2023c | maillard_path_holdout | Furfural | ppb | 113.000 | 4.595 | 0.028 | 1.727 | 69.759 | no | 3.400 | yes (headspace) | sulfur | external_holdout |
| mp_holdout_arp_glutathione_120C_60min_pH75_Zhou2023c | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 64.100 | 1.75e+03 | 18.278 | 324.855 | 3.93e+03 | yes | 2.333 | yes (headspace) | sulfur | external_holdout |
| mp_holdout_arp_glutathione_120C_60min_pH75_Zhou2023c | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 383.000 | 73.488 | 0.025 | 4.560 | 445.495 | yes | 4.258 | yes (headspace) | sulfur | external_holdout |
| mp_holdout_arp_glutathione_120C_60min_pH75_Zhou2023c | maillard_path_holdout | Furfural | ppb | 397.000 | 8.835 | 0.030 | 3.386 | 194.757 | no | 3.814 | yes (headspace) | sulfur | external_holdout |
| mp_holdout_arp_glutathione_120C_90min_pH75_Zhou2023c | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 57.100 | 1.65e+03 | 20.615 | 348.926 | 3.95e+03 | yes | 2.282 | yes (headspace) | sulfur | external_holdout |
| mp_holdout_arp_glutathione_120C_90min_pH75_Zhou2023c | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 270.000 | 110.837 | 0.024 | 4.762 | 576.440 | yes | 4.387 | yes (headspace) | sulfur | external_holdout |
| mp_holdout_arp_glutathione_120C_90min_pH75_Zhou2023c | maillard_path_holdout | Furfural | ppb | 575.000 | 16.172 | 0.021 | 4.112 | 325.617 | no | 4.190 | yes (headspace) | sulfur | external_holdout |
| mp_holdout_fructose_asparagine_180C_Lin2022 | maillard_path_holdout | Acrylamide | ppb | 1.86e+03 | 203.194 | 68.996 | 199.516 | 565.424 | no | 0.914 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_fructose_asparagine_180C_Lin2022 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | ppb | 1.23e+04 | 1.97e+03 | 1.1e+03 | 1.88e+03 | 2.73e+03 | no | 0.395 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | maillard_path_holdout | DMHF | ppb | 1.15e+03 | 21.841 | 6.595 | 21.187 | 60.187 | no | 0.960 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | ppb | 5.73e+04 | 5.54e+04 | 8.91e+03 | 3.3e+04 | 1.59e+05 | yes | 1.253 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | maillard_path_holdout | DMHF | ppb | 5.89e+03 | 21.841 | 6.595 | 21.187 | 60.187 | no | 0.960 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | ppb | 1.01e+05 | 5.54e+04 | 8.91e+03 | 3.3e+04 | 1.59e+05 | yes | 1.253 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_asparagine_180C_10min_Chang2021 | maillard_path_holdout | Acrylamide | ppb | 28.000 | 6.11e+03 | 1.94e+03 | 5.68e+03 | 1.42e+04 | no | 0.865 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_asparagine_180C_30min_Chang2021 | maillard_path_holdout | Acrylamide | ppb | 1.46e+03 | 3.68e+03 | 1.25e+03 | 3.62e+03 | 1.03e+04 | yes | 0.915 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | maillard_path_holdout | Acrylamide | ppb | 832.000 | 3.68e+03 | 1.25e+03 | 3.62e+03 | 1.03e+04 | no | 0.915 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | ppb | 7e+03 | 2.55e+03 | 1.12e+03 | 2.18e+03 | 1.21e+04 | yes | 1.031 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_asparagine_180C_Ye2024 | maillard_path_holdout | Acrylamide | umol_per_mol_limiting_precursor | 140.580 | 7e+03 | 2.43e+03 | 6.84e+03 | 1.94e+04 | no | 0.902 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | ppb | 1.74e+04 | 1.87e+03 | 471.658 | 1.32e+03 | 4.41e+03 | no | 0.970 | no (extraction) | trunk | external_holdout |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | maillard_path_holdout | 3-deoxyglucosone | ppb | 5.22e+04 | 6.64e+04 | 1.68e+04 | 4.85e+04 | 1.57e+05 | yes | 0.969 | no (extraction) | trunk | external_holdout |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | maillard_path_holdout | 3,4-dideoxyglucosone | ppb | 5.55e+04 | 8.37e+03 | 2.1e+03 | 5.79e+03 | 2.04e+04 | no | 0.988 | no (extraction) | trunk | external_holdout |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | maillard_path_holdout | glucosone | ppb | 7.5e+03 | 119.059 | 118.920 | 119.086 | 119.141 | no | 0.000807 | no (extraction) | trunk | external_holdout |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | maillard_path_holdout | glyoxal | ppb | 5.6e+03 | 160.794 | 168.424 | 200.190 | 205.037 | no | 0.085 | no (extraction) | trunk | external_holdout |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | maillard_path_holdout | methylglyoxal | ppb | 2.6e+03 | 2.57e+03 | 622.049 | 1.72e+03 | 5.59e+03 | yes | 0.953 | no (extraction) | trunk | external_holdout |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 229.000 | 1.52e+03 | 56.235 | 710.590 | 6.78e+03 | yes | 2.081 | no (extraction) | sulfur | external_holdout (shared: hofmann_ribose_pH3_FFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 553.000 | 24.707 | 8.920 | 25.378 | 70.653 | no | 0.899 | no (extraction) | sulfur | external_holdout (shared: hofmann_ribose_pH3_MFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 12.000 | 0.170 | 0.00319 | 0.075 | 0.328 | no | 2.012 | no (extraction) | sulfur | external_holdout (shared: hofmann_ribose_pH7_FFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 25.000 | 0.0051 | 3.28e-05 | 0.00429 | 0.021 | no | 2.809 | no (extraction) | sulfur | external_holdout (shared: hofmann_ribose_pH7_MFT) |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 96.000 | 346.617 | 16.447 | 177.976 | 1.22e+03 | yes | 1.869 | no (extraction) | sulfur | external_holdout |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 143.000 | 93.633 | 37.127 | 95.921 | 219.557 | yes | 0.772 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 6.880 | 10.334 | 0.194 | 3.569 | 83.929 | yes | 2.636 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 1.280 | 299.499 | 0.271 | 12.621 | 422.813 | yes | 3.194 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 3.290 | 12.015 | 0.596 | 5.519 | 76.404 | yes | 2.108 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 1.460 | 246.542 | 0.621 | 19.017 | 367.092 | yes | 2.772 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 2.400 | 12.859 | 1.389 | 7.672 | 52.086 | yes | 1.574 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 1.680 | 167.266 | 1.658 | 23.995 | 285.632 | yes | 2.236 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 1.710 | 12.681 | 2.527 | 9.590 | 35.304 | no | 1.145 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 1.620 | 96.098 | 2.626 | 23.394 | 166.040 | no | 1.801 | no (extraction) | sulfur | external_holdout |
| mp_holdout_xylose_cysteine_phosphate_120C_50min_pH65_Zhang2025b | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 22.630 | 0.510 | 0.023 | 0.336 | 3.353 | no | 2.173 | yes (headspace) | sulfur | external_holdout |
| mp_holdout_xylose_cysteine_phosphate_120C_50min_pH65_Zhang2025b | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 0.290 | 2.414 | 0.014 | 0.406 | 6.299 | yes | 2.654 | yes (headspace) | sulfur | external_holdout |
| mp_holdout_xylose_cysteine_phosphate_120C_50min_pH65_Zhang2025b | maillard_path_holdout | Bis(2-methyl-3-furyl) disulfide | ppb | 18.250 | 0.000636 | 2.91e-06 | 0.000195 | 0.011 | no | 3.577 | yes (headspace) | sulfur | external_holdout |
| mp_holdout_xylose_cysteine_phosphate_120C_50min_pH65_Zhang2025b | maillard_path_holdout | 2-Acetylthiazole | ppb | 114.130 | 0.065 | 0.00455 | 0.049 | 0.455 | no | 1.999 | yes (headspace) | sulfur | external_holdout |
| external_validation_bi_2020_roasted_pea_hexanal | external_matrix | hexanal | ppb | 324.000 | 109.490 | 9.240 | 109.689 | 1.6e+03 | yes | 2.237 | yes (headspace) | lipid | external_holdout |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | 2-pentylfuran | ppb | 5.63e+03 | 15.387 | 1.115 | 13.019 | 171.111 | no | 2.186 | yes (headspace) | lipid | external_holdout |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | hexanal | ppb | 605.600 | 70.264 | 5.084 | 59.299 | 779.558 | yes | 2.186 | yes (headspace) | lipid | external_holdout |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | nonanal | ppb | 74.370 | 23.885 | 1.732 | 20.209 | 265.608 | yes | 2.186 | yes (headspace) | lipid | external_holdout |

## Refused rows

| benchmark | panel | compound | reason |
|---|---|---|---|
| cml_cel_commercial_pbma_Foods2023 | trust_loop | Nε-(Carboxymethyl)lysine (CML) | GLYCATION TARGETS 'Nε-(Carboxymethyl)lysine (CML)' (wave B20) run on the trunk lane only: the acrylamide lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + ami |
| cml_cel_commercial_pbma_Foods2023 | trust_loop | Nε-(Carboxyethyl)lysine (CEL) | GLYCATION TARGETS 'Nε-(Carboxyethyl)lysine (CEL)' (wave B20) run on the trunk lane only: the acrylamide lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + amin |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| pea_isolate_40C_PratapSingh2021 | trust_loop | hexanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's own vessel says so (closure = 'no cook'), and the physics agrees: over this thermal progr |
| pea_isolate_40C_PratapSingh2021 | trust_loop | 2-pentylfuran | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's own vessel says so (closure = 'no cook'), and the physics agrees: over this thermal progr |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | 2,5-dimethylpyrazine | PYRAZINE TARGETS '2,5-dimethylpyrazine' (wave B18) run on the trunk lane only: the lipid lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + amine pot that reso |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | methional | METHIONINE CHAIN TARGETS 'methional' (wave B22) run on the trunk lane only: the lipid lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + amine pot that resolve |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | 2-acetyl-1-pyrroline | 2-ACETYL-1-PYRROLINE TARGETS '2-acetyl-1-pyrroline' (wave B24) run on the trunk lane only: the lipid lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + amine p |
| soy_isolate_40C_PratapSingh2021 | trust_loop | hexanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's own vessel says so (closure = 'no cook'), and the physics agrees: over this thermal progr |
| soy_isolate_40C_PratapSingh2021 | trust_loop | 2-pentylfuran | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's own vessel says so (closure = 'no cook'), and the physics agrees: over this thermal progr |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | maillard_path_holdout | Furfural | LANE CONFLICT: this request needs both the acrylamide and sulfur lanes at once. They do not compose -- the acrylamide network deliberately omits every sulfur step (acrylamide.OUT_OF_SCOPE), because composing them would s |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | maillard_path_holdout | Furfural | LANE CONFLICT: this request needs both the acrylamide and sulfur lanes at once. They do not compose -- the acrylamide network deliberately omits every sulfur step (acrylamide.OUT_OF_SCOPE), because composing them would s |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| external_validation_bi_2020_raw_pea_hexanal | external_matrix | hexanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's own vessel says so (closure = 'no cook'), and the physics agrees: over this thermal progr |
| external_validation_bi_2020_raw_pea_hexanal | external_matrix | nonanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's own vessel says so (closure = 'no cook'), and the physics agrees: over this thermal progr |
| external_validation_bi_2020_roasted_pea_hexanal | external_matrix | furaneol | THIS POT CHARGES NO PRECURSOR THAT COULD MAKE 'furaneol'. The charge declares only a matrix/lipid carrier, which is not a precursor: it resolves to a hydroperoxide pool for the lipid lane and charges NOTHING on the trunk |
| external_validation_bi_2020_roasted_pea_hexanal | external_matrix | 2,5-dimethylpyrazine | PYRAZINE TARGETS '2,5-dimethylpyrazine' (wave B18) run on the trunk lane only: the lipid lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + amine pot that reso |
| external_validation_bi_2020_roasted_pea_hexanal | external_matrix | furfural | THIS POT CHARGES NO PRECURSOR THAT COULD MAKE 'furfural'. The charge declares only a matrix/lipid carrier, which is not a precursor: it resolves to a hydroperoxide pool for the lipid lane and charges NOTHING on the trunk |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | 1-hexanol | UNREPRESENTED TARGETS: 1-hexanol -- The lipid lane exists and forms the SIX products Frankel 1989 measured, but 1-hexanol is not one of them and NO aldehyde-reduction step is measured anywhere in the corpus -- in a therm |
| external_validation_liu_2023_ppi_offnote_baseline | external_matrix | hexanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's own vessel says so (closure = 'no cook'), and the physics agrees: over this thermal progr |
| external_validation_liu_2023_ppi_offnote_baseline | external_matrix | nonanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's own vessel says so (closure = 'no cook'), and the physics agrees: over this thermal progr |

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
| b8.k_arp_dpo.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -2.634 | 2.500 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_arp_tdp.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -1.965 | 0.939 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_dpo_nf.log10_k_ref_145C | sulfur | fitted_rate | fixed | 0.434 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_dpo_ptr.log10_k_ref_145C | sulfur | fitted_rate | fixed | -4.347 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_dpo_ddp.log10_k_ref_145C | sulfur | fitted_rate | fixed | 0.042 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_tdp_fur.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -3.491 | 0.376 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_ddp_mft.log10_k_ref_145C | sulfur | fitted_rate | fixed | -6.541 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_fur_fft.log10_k_ref_145C | sulfur | fitted_rate | uniform_band | -1.643 | 15.718 | [-10.000, 0.500] | yes | unidentified_in_the_fit: sampled across its declared band |
| b8.k_nf_mft.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -2.666 | 0.274 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_nf_mp3p.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -3.793 | 0.578 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_mgo_mp.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -3.946 | 0.398 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_ha_mp_mft.log10_k_ref_145C | sulfur | fitted_rate | fixed | -3.870 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_glc_ha.log10_k_ref_145C | sulfur | fitted_rate | uniform_band | -8.589 | - | [-10.000, 0.500] | yes | unidentified_in_the_fit: sampled across its declared band |
| b8.k_thi_hmp.log10_k_ref_145C | sulfur | fitted_rate | fixed | -2.588 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_thi_mesh.log10_k_ref_145C | sulfur | fitted_rate | fixed | -3.273 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_hmp_mft.log10_k_ref_145C | sulfur | fitted_rate | fixed | -2.609 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_hmp_mp2p.log10_k_ref_145C | sulfur | fitted_rate | fixed | -3.437 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_cys_actz.log10_k_ref_145C | sulfur | fitted_rate | fixed | -2.974 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_dimer_mft.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | 0.500 | 0.404 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_dimer_fft.log10_k_ref_145C | sulfur | fitted_rate | uniform_band | -1.460 | 67.997 | [-10.000, 0.441] | yes | unidentified_in_the_fit: sampled across its declared band |
| b8.k_mmft.log10_k_ref_145C | sulfur | fitted_rate | fixed | -1.984 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_mft_decay.log10_k_ref_145C | sulfur | fitted_rate | fixed | 0.134 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_fft_decay.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -0.828 | 0.315 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_dimer_decay.log10_k_ref_145C | sulfur | fitted_rate | uniform_band | -9.139 | - | [-10.000, -0.868] | yes | unidentified_in_the_fit: sampled across its declared band |
| b8.k_nf_decay.log10_k_ref_145C | sulfur | fitted_rate | fixed | -2.182 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_fur_decay.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -0.417 | 0.884 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_h2s_loss.log10_k_ref_145C | sulfur | fitted_rate | fixed | -1.351 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_osone_decay.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -1.162 | 0.604 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_thiol_decay.log10_k_ref_145C | sulfur | fitted_rate | uniform_band | -7.614 | - | [-10.000, -1.353] | yes | unidentified_in_the_fit: sampled across its declared band |
| b8.k_pent_caramel.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -2.249 | 0.833 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_pent_thermal.log10_k_ref_145C | sulfur | fitted_rate | fixed | -3.437 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_glc_fur.log10_k_ref_145C | sulfur | fitted_rate | fixed | -4.923 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_arp_tdp_th.log10_k_ref_145C | sulfur | fitted_rate | fixed | -4.354 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_arp_dpo_th.log10_k_ref_145C | sulfur | fitted_rate | fixed | -1.290 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_ddp_mft_hs.log10_k_ref_145C | sulfur | fitted_rate | fixed | -7.693 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_fur_fft_hs.log10_k_ref_145C | sulfur | fitted_rate | fixed | -0.00955 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_ttca_cys.log10_k_ref_145C | sulfur | fitted_rate | fixed | -0.579 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.k_ttca_deg.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -2.212 | 1.892 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_cys_thermal.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -2.124 | 0.353 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.k_thiolate_loss.log10_k_ref_145C | sulfur | fitted_rate | normal_log10 | -1.547 | 0.471 | [-10.000, 0.500] | yes | laplace_covariance_at_b8_optimum |
| b8.lumped_formation_Ea_kJ_mol | sulfur | fitted_ea | fixed | 64.082 | - | - | no | frozen in the sulfur fit (not a free coordinate; source uncertainty unrecorded) |
| b8.decay_Ea_kJ_mol.thiol_sink | sulfur | fitted_ea | normal | 101.254 | 32.822 | [7.000, 102.000] | yes | laplace_covariance_at_b8_optimum |
| b8.decay_Ea_kJ_mol.carbonyl_sink | sulfur | fitted_ea | uniform_band | 143.538 | 140.102 | [95.506, 191.570] | yes | unidentified_in_the_fit: declared band capped by the prefactor prior |
| b8.ph_drift.acid_yield_per_sink_event | sulfur | fitted_ph_drift | fixed | 7.5e-05 | 0.535 | [0.000, 1.000] | no | unidentified_in_the_fit: its bound is DEFINITIONAL (the quantity is a fraction), not a measurement of where the value li |
| b8.ph_drift.arp_secondary_ammonium_pKa | sulfur | fitted_ph_drift | normal | 7.068 | 0.599 | [5.000, 11.000] | yes | laplace_covariance_at_b8_optimum |
| b13.k_da_sink.log10_k_100C | trunk | fitted_rate | uniform_band | - | - | [-12.000, -1.940] | yes | two_laboratories_disagree: flat across the disagreement |
| b13.k_go_sink.ea_kj_mol | trunk | fitted_ea | uniform_band | 0.000 | - | [0.000, 150.800] | yes | two_laboratories_disagree: flat across the disagreement |
| b13.k_odg_da.log10_k_100C | trunk | fitted_rate | uniform_band | -5.640 | - | [-4.630, -1.960] | yes | two_laboratories_disagree: flat across the disagreement |
| b13.k_hmf_self.log10_k_100C | trunk | fitted_rate | uniform_band | -6.047 | - | [-6.050, -0.960] | yes | two_laboratories_disagree: flat across the disagreement |
| b13.k_glc_g.log10_k_100C | trunk | fitted_rate | fixed | -7.272 | - | - | no | NOT SAMPLED, and not by oversight. k_ddg_hmf (1.13x) and k_go_sink's RATE (1.87x) are cross-laboratory agreements inside |
| b13.k_g_go.log10_k_100C | trunk | fitted_rate | fixed | -2.451 | - | - | no | NOT SAMPLED, and not by oversight. k_ddg_hmf (1.13x) and k_go_sink's RATE (1.87x) are cross-laboratory agreements inside |
| b13.k_ddg_hmf.log10_k_100C | trunk | fitted_rate | fixed | -0.924 | - | - | no | NOT SAMPLED, and not by oversight. k_ddg_hmf (1.13x) and k_go_sink's RATE (1.87x) are cross-laboratory agreements inside |
| b34.k_glc_tdg.log10_k_100C | trunk | fitted_rate | uniform_band | -5.027 | - | [-5.406, -4.828] | yes | SAMPLED over the source's own printed 95 % HPD (Kocadagli & Gokmen 2016 Table 2). The interval is the measurement's, not |
| b34.k_glc_tdg.ea_kj_mol | trunk | fitted_ea | uniform_band | 107.200 | - | [54.500, 159.900] | yes | SAMPLED over the source's own printed 95 % HPD (Kocadagli & Gokmen 2016 Table 2). The interval is the measurement's, not |
| b34.k_fru_int.log10_k_100C | trunk | fitted_rate | uniform_band | -2.963 | - | [-2.994, -2.934] | yes | SAMPLED over the source's own printed 95 % HPD (Kocadagli & Gokmen 2016 Table 2). The interval is the measurement's, not |
| b34.k_fru_int.ea_kj_mol | trunk | fitted_ea | uniform_band | 100.400 | - | [93.800, 107.000] | yes | SAMPLED over the source's own printed 95 % HPD (Kocadagli & Gokmen 2016 Table 2). The interval is the measurement's, not |
| b34.k_fru_odg.log10_k_100C | trunk | fitted_rate | uniform_band | -5.130 | - | [-5.221, -5.054] | yes | SAMPLED over the source's own printed 95 % HPD (Kocadagli & Gokmen 2016 Table 2). The interval is the measurement's, not |
| b34.k_fru_odg.ea_kj_mol | trunk | fitted_ea | uniform_band | 99.300 | - | [77.500, 121.100] | yes | SAMPLED over the source's own printed 95 % HPD (Kocadagli & Gokmen 2016 Table 2). The interval is the measurement's, not |
| b18.log10_k_go_ak_100C | trunk | fitted_rate | normal_log10 | -6.542 | 0.082 | [-9.191, -5.191] | yes | laplace_covariance_at_b8_optimum |
| b18.ea_go_ak_kj_mol | trunk | fitted_ea | uniform_band | 103.100 | 16.158 | [100.590, 103.100] | yes | unidentified_direction_flat_across_its_declared_band |
| b18.log10_k_mgo_ak_100C | trunk | fitted_rate | normal_log10 | -7.530 | 0.075 | [-10.093, -6.093] | yes | laplace_covariance_at_b8_optimum |
| b18.ea_mgo_ak_kj_mol | trunk | fitted_ea | uniform_band | 114.900 | 16.097 | [111.660, 114.900] | yes | unidentified_direction_flat_across_its_declared_band |
| b18.ph_slope_above_7_decades_per_unit | trunk | fitted_rate | fixed | 0.197 | 0.039 | [0.000, 1.500] | no | IDENTIFIED AND STILL NOT SAMPLED, and the interval is narrower for it. This coordinate is a module-level constant in par |
| b18.ph_slope_below_7_decades_per_unit | trunk | fitted_rate | fixed | 0.580 | 0.046 | [0.000, 1.500] | no | IDENTIFIED AND STILL NOT SAMPLED, and the interval is narrower for it. This coordinate is a module-level constant in par |
| b39.log10_k_tdg_ddg_100C | trunk | fitted_rate | normal_log10 | -1.933 | 0.062 | [-4.428, -0.428] | yes | laplace_covariance_at_b8_optimum |
| b39.log10_k_ddg_tdg_100C | trunk | fitted_rate | normal_log10 | -1.965 | 0.179 | [-4.428, -0.428] | yes | laplace_covariance_at_b8_optimum |
| b39.log10_k_ddg_dgal_100C | trunk | fitted_rate | normal_log10 | -1.812 | 0.243 | [-4.428, -0.428] | yes | laplace_covariance_at_b8_optimum |
| b39.log10_k_dgal_ddg_100C | trunk | fitted_rate | normal_log10 | -1.698 | 0.258 | [-4.428, -0.428] | yes | laplace_covariance_at_b8_optimum |
| b39.log10_k_ddg_hmf_100C | trunk | fitted_rate | normal_log10 | -1.450 | 0.064 | [-2.924, 1.076] | yes | laplace_covariance_at_b8_optimum |
| lipid.q10 | lipid | declared_band | uniform | 2.449 | - | [2.000, 3.000] | yes | declared corner band, sampled uniform over it |
| lipid.pea_protein_isolate.lipid_mass_fraction | lipid | declared_band | log_uniform | 0.025 | - | [0.010, 0.060] | yes | declared corner band, sampled log-uniform over it as ONE scale shared by every carrier in a draw (CoreDraw.lipid_fractio |
| lipid.pea_protein_isolate.peroxide_value_meq_per_kg | lipid | declared_band | log_uniform | 10.000 | - | [2.000, 40.000] | yes | declared corner band, sampled log-uniform over it as ONE scale shared by every carrier in a draw (CoreDraw.peroxide_scal |
| lipid.soy_protein_isolate.lipid_mass_fraction | lipid | declared_band | log_uniform | 0.020 | - | [0.008, 0.050] | yes | declared corner band, sampled log-uniform over it as ONE scale shared by every carrier in a draw (CoreDraw.lipid_fractio |
| lipid.soy_protein_isolate.peroxide_value_meq_per_kg | lipid | declared_band | log_uniform | 10.000 | - | [2.000, 40.000] | yes | declared corner band, sampled log-uniform over it as ONE scale shared by every carrier in a draw (CoreDraw.peroxide_scal |
| lipid.frankel_pure_hydroperoxide.lipid_mass_fraction | lipid | declared_band | fixed | 1.000 | - | [1.000, 1.000] | no | degenerate band (fed hydroperoxide: the fraction is the definition) |
| lipid.frankel_pure_hydroperoxide.peroxide_value_meq_per_kg | lipid | declared_band | fixed | 2e+03 | - | [2e+03, 2e+03] | no | degenerate band (fed hydroperoxide: PV is the definition) |
| trunk.aw_multiplier_scale | trunk | declared_band | uniform | 1.000 | - | [0.000, 1.200] | yes | declared band: 0 = Bell 1995's fixed-molality plateau, 1.2 = the source's 95 % CI |
| trunk.amadori_ph_exponent_decades_per_unit | trunk | declared_band | uniform | 0.690 | - | [0.370, 0.920] | yes | declared band: the six Martins 2003 per-step ratios span it |
| acrylamide.aw_multiplier | acrylamide | declared_band | uniform | 1.000 | - | [0.410, 1.390] | yes | declared band: the source's a_w point estimates and the 0.92 column's 95 % HPD, relative to the shipped constant |
| acrylamide.aw_elimination_deficit_scale | acrylamide | declared_band | uniform | 1.000 | - | [0.000, 1.200] | yes | declared band: 0 = no elimination a_w effect, 1.2 = 1.2x De Vleeschouwer 2007's k_E shape (SEs up to 90 %) |
| acrylamide.ph_exponent_formation_decades_per_unit | acrylamide | declared_band | uniform | 0.235 | - | [0.114, 0.281] | yes | declared band: the potato-matrix slope minus its SE to the phosphate slope plus its SE |
| acrylamide.ph_exponent_elimination_decades_per_unit | acrylamide | declared_band | uniform | 0.149 | - | [0.116, 0.155] | yes | declared band: the two measured slopes and their SEs |
| sulfur.oxygen.k_cys_ox.log10_k | sulfur | declared_band | fixed | 0.000 | - | [-5.000, -1.000] | no | not a free coordinate of the shipped report (B11 not shipped): the consumer is zero by declaration and the structure ine |
| sulfur.oxygen.k_red_ox.log10_k | sulfur | declared_band | fixed | 0.000 | - | [-5.000, -1.000] | no | not a free coordinate of the shipped report (B11 not shipped): the consumer is zero by declaration and the structure ine |
| sulfur.oxygen.reservoir_scale | sulfur | declared_band | fixed | 1.000 | - | [0.100, 10.000] | no | not a free coordinate of the shipped report (B11 not shipped): the reservoir is inert while the consumers are zero |
| furanic.partition_ea_offset_kj_mol | trunk | declared_band | uniform | 0.000 | - | [-50.000, 50.000] | yes | declared corner band on the furanone PARTITION barrier, sampled uniform |
| observable.air_water_partition_constant | observable | observable | log_uniform | 0.000 | - | [-0.500, 0.500] | yes | declared +/-0.5 dex band on K_aw, sampled uniform in log10 |
| observable.hs_spme_same_sample_dispersion | observable | observable | log_uniform_dispersion | 15.166 | - | [10.000, 23.000] | yes | measured 10-23x same-sample dispersion: D is drawn log-uniform over [10, 23] and the multiplier is D^u with u uniform on |
| b55.k_pent_dpo.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [30.000, 237.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_pent_tdp.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [30.000, 237.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_dpo_c2c3.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [30.000, 237.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_dpo_nf.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [30.000, 237.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_dpo_ptr.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [30.000, 237.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_dpo_ddp.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [30.000, 237.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_tdp_fur.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [30.000, 237.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_ddp_mft.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [29.675, 133.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_fur_fft.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [29.675, 133.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_nf_mft.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [29.675, 133.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_nf_mp3p.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [29.675, 133.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_mgo_mp.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [29.675, 133.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_ha_mp_mft.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [29.675, 133.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_glc_ha.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [30.000, 237.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_thi_hmp.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [29.675, 133.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_thi_mesh.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [29.675, 133.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_hmp_mft.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [29.675, 133.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_hmp_mp2p.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [29.675, 133.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_cys_actz.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [29.675, 133.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_mmft.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | - | - | [29.675, 133.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_dimer_decay.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [29.675, 133.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_h2s_loss.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [29.675, 133.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_pent_caramel.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [30.000, 237.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_pent_thermal.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [30.000, 237.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_glc_fur.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [30.000, 237.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_arp_tdp_th.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [30.000, 237.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_arp_dpo_th.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [30.000, 237.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_ddp_mft_hs.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [29.675, 133.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_fur_fft_hs.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [29.675, 133.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_ttca_cys.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [30.000, 237.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
| b55.k_ttca_deg.ea_kj_mol | sulfur | fitted_ea | empirical_measured_class | 64.082 | - | [30.000, 237.000] | yes | B55: the barrier was never measured (lumped formation Ea, unidentified since B2.3, or none at all); drawn from the EMPIR |
