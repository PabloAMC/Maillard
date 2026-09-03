# Core panel scorecard (kinetic core on the union panel)

pass band = 3.0x; contracts from each bundle's `validation_contract.scale_thresholds`, else the global default.

* panel: **40** benchmarks, 32 scored, 8 fully refused; rows **49**, refused rows 18
* within 3x: 6/49 (0.122); median fold 50.058, geometric mean 173.403, worst 8.39e+10
* evidence roles (core): {'external_holdout': 21, 'predictive': 19}
* predictive benchmarks passing their contract: thiamine_cys_glucose_120C_Bolton1994; strict-ready: thiamine_cys_glucose_120C_Bolton1994
* **honest literature: 6/49 within band** (0.122), 32 benchmarks, median fold 50.058, geometric mean 173.403
* **out-of-sample: 5/48 within band** (0.104), 31 benchmarks, median fold 51.180, geometric mean 189.706
* **rows the sulfur fit read: 1/1 within band** (1.000), 1 benchmarks, median fold 2.323, geometric mean 2.323

## Per panel / role / lane

| split | key | benchmarks | rows | within band | rate | contract passes | strict-ready | median fold | geo-mean fold |
|---|---|---|---|---|---|---|---|---|---|
| panel | external_matrix | 4 | 4 | 0 | 0.000 | 0 | 0 | 1.86e+03 | 250.597 |
| panel | maillard_path_holdout | 17 | 30 | 4 | 0.133 | 0 | 0 | 39.767 | 85.044 |
| panel | trust_loop | 19 | 15 | 2 | 0.133 | 1 | 1 | 52.302 | 653.494 |
| evidence_role | external_holdout | 21 | 34 | 4 | 0.118 | 0 | 0 | 39.767 | 96.574 |
| evidence_role | predictive | 19 | 15 | 2 | 0.133 | 1 | 1 | 52.302 | 653.494 |
| lane | acrylamide | - | 12 | 2 | 0.167 | 0 | 0 | 7.285 | 21.945 |
| lane | lipid | - | 7 | 0 | 0.000 | 0 | 0 | 3.36e+03 | 430.826 |
| lane | sulfur | - | 29 | 4 | 0.138 | 0 | 0 | 56.130 | 359.098 |
| lane | trunk | - | 1 | 0 | 0.000 | 0 | 0 | 11.930 | 11.930 |

## Benchmarks

| benchmark | panel | tier | role | rows | coverage | max ratio | mean log10 | contract (ratio / log10) | status | strict | in core fit |
|---|---|---|---|---|---|---|---|---|---|---|---|
| acrylamide_spi_extrusion_130C_ACSRef3 | trust_loop | PRIMARY | predictive | 1/1 | 1.000 | 4.25e+03 | 3.628 | 1.50 / 0.200 | scale-gap | no | - |
| cml_cel_commercial_pbma_Foods2023 | trust_loop | PRIMARY | predictive | 0/2 | 0.000 | - | - | 1.80 / 0.250 | refused | no | - |
| cys_ribose_140C_Hofmann1998 | trust_loop | REFERENCE | predictive | 2/2 | 1.000 | 4.484 | 0.596 | 1.50 / 0.100 | scale-gap | no | - |
| furosine_extrusion_crossover_140C_RamirezJimenez2000 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 2.00 / 0.300 | refused | no | - |
| hofmann1998_c2c3_recombination_145C_20min_pH3 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - |
| hofmann1998_c2c3_recombination_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - |
| hofmann1998_c2c3_recombination_145C_20min_pH7 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 2/2 | 1.000 | 8.39e+10 | 6.321 | 1.10 / 0.041 | scale-gap | no | - |
| hofmann1998_furan2aldehyde_h2s_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 2/2 | 1.000 | 1.39e+10 | 6.028 | 1.10 / 0.041 | scale-gap | no | - |
| hofmann1998_norfuraneol_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 1/1 | 1.000 | 2.323 | 0.366 | 1.10 / 0.041 | scale-gap | no | 1 |
| hofmann1998_norfuraneol_h2s_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 2/2 | 1.000 | 6.896 | 0.703 | 1.10 / 0.041 | scale-gap | no | - |
| pea_isolate_40C_PratapSingh2021 | trust_loop | PRIMARY | predictive | 1/2 | 0.500 | 3.36e+03 | 3.526 | 2.00 / 0.120 | coverage-gap | no | - |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | PRIMARY | predictive | 1/3 | 0.333 | 34.240 | 1.535 | 2.00 / 0.120 | coverage-gap | no | - |
| resconi_2023_pbma_beef_identity_benchmark | trust_loop | SECONDARY | predictive | 1/1 | 1.000 | 93.308 | 1.970 | 1.50 / 0.100 | scale-gap | no | - |
| soy_isolate_40C_PratapSingh2021 | trust_loop | PRIMARY | predictive | 1/2 | 0.500 | 6.08e+03 | 3.784 | 2.00 / 0.120 | coverage-gap | no | - |
| thiamine_cys_glucose_120C_Bolton1994 | trust_loop | PRIMARY | predictive | 1/1 | 1.000 | 1.332 | 0.124 | 3.00 / 0.480 | pass-no-ranking | yes | - |
| thiamine_cys_xylose_145C_Cerny2008 | trust_loop | PRIMARY | predictive | 0/0 | 0.000 | - | - | 1.50 / 0.100 | refused | no | - |
| mp_holdout_fructose_asparagine_180C_Lin2022 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 8.246 | 0.859 | 1.50 / 0.100 | scale-gap | no | - |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | maillard_path_holdout | REFERENCE | external_holdout | 2/3 | 0.667 | 52.799 | 1.017 | 1.50 / 0.100 | coverage-gap | no | - |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | maillard_path_holdout | REFERENCE | external_holdout | 2/3 | 0.667 | 269.863 | 1.494 | 1.50 / 0.100 | coverage-gap | no | - |
| mp_holdout_glucose_asparagine_180C_10min_Chang2021 | maillard_path_holdout | REFERENCE | external_holdout | 1/1 | 1.000 | 241.056 | 2.382 | 1.50 / 0.100 | scale-gap | no | - |
| mp_holdout_glucose_asparagine_180C_30min_Chang2021 | maillard_path_holdout | REFERENCE | external_holdout | 1/1 | 1.000 | 2.763 | 0.441 | 1.50 / 0.100 | scale-gap | no | - |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 4.845 | 0.603 | 1.50 / 0.100 | scale-gap | no | - |
| mp_holdout_glucose_asparagine_180C_Ye2024 | maillard_path_holdout | REFERENCE | external_holdout | 1/1 | 1.000 | 50.058 | 1.699 | 1.50 / 0.100 | scale-gap | no | - |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | maillard_path_holdout | REFERENCE | external_holdout | 1/1 | 1.000 | 11.930 | 1.077 | 1.50 / 0.100 | scale-gap | no | - |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 2.15e+09 | 5.705 | 1.10 / 0.041 | scale-gap | no | - |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 3.24e+09 | 4.904 | 1.10 / 0.041 | scale-gap | no | - |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 23.562 | 1.034 | 1.10 / 0.041 | scale-gap | no | - |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 140.309 | 2.136 | 1.10 / 0.041 | scale-gap | no | - |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 3.835 | 0.443 | 1.10 / 0.041 | scale-gap | no | - |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 480.873 | 1.826 | 1.50 / 0.100 | scale-gap | no | - |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 371.077 | 2.019 | 1.50 / 0.100 | scale-gap | no | - |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 218.874 | 2.045 | 1.50 / 0.100 | scale-gap | no | - |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 131.941 | 2.060 | 1.50 / 0.100 | scale-gap | no | - |
| external_validation_bi_2020_raw_pea_hexanal | external_matrix | PRIMARY | external_holdout | 1/1 | 1.000 | 3.72e+03 | 3.570 | 2.00 / 0.120 | scale-gap | no | - |
| external_validation_bi_2020_roasted_pea_hexanal | external_matrix | PRIMARY | external_holdout | 1/1 | 1.000 | 3.657 | 0.563 | 2.00 / 0.120 | scale-gap | no | - |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | PRIMARY | external_holdout | 1/4 | 0.250 | 8.689 | 0.939 | 2.00 / 0.120 | coverage-gap | no | - |
| external_validation_liu_2023_ppi_offnote_baseline | external_matrix | PRIMARY | external_holdout | 1/2 | 0.500 | 3.34e+04 | 4.524 | 2.00 / 0.120 | coverage-gap | no | - |

## Rows

| benchmark | compound | unit | measured | predicted | fold | within band | within contract | interval (ug/L) | measured inside | lane | in core fit |
|---|---|---|---|---|---|---|---|---|---|---|---|
| acrylamide_spi_extrusion_130C_ACSRef3 | acrylamide | ppb | 150.000 | 0.035 | 4.25e+03 | no | no | [0.00505, 0.247] | no | acrylamide | no |
| cys_ribose_140C_Hofmann1998 | 2-methyl-3-furanthiol | ppb | 342.000 | 76.277 | 4.484 | no | no | [10.906, 533.489] | yes | sulfur | no |
| cys_ribose_140C_Hofmann1998 | 2-furfurylthiol | ppb | 200.000 | 57.694 | 3.467 | no | no | [8.249, 403.517] | yes | sulfur | no |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | 2-Furfurylthiol (FFT) | ppb | 32.000 | 1.67e+03 | 52.302 | no | no | [239.297, 1.17e+04] | no | sulfur | no |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | ppb | 25.000 | 2.98e-10 | 8.39e+10 | no | no | [4.26e-11, 2.08e-09] | no | sulfur | no |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | 2-Furfurylthiol (FFT) | ppb | 28.000 | 2.3e+03 | 82.016 | no | no | [328.338, 1.61e+04] | no | sulfur | no |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | ppb | 19.000 | 1.37e-09 | 1.39e+10 | no | no | [1.96e-10, 9.57e-09] | no | sulfur | no |
| hofmann1998_norfuraneol_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | ppb | 1.02e+03 | 2.36e+03 | 2.323 | yes | no | [337.400, 1.65e+04] | yes | sulfur | yes |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | 2-Furfurylthiol (FFT) | ppb | 121.000 | 834.391 | 6.896 | no | no | [119.299, 5.84e+03] | yes | sulfur | no |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | ppb | 198.000 | 731.065 | 3.692 | no | no | [104.525, 5.11e+03] | yes | sulfur | no |
| pea_isolate_40C_PratapSingh2021 | hexanal | ppb | 1.14e+03 | 0.339 | 3.36e+03 | no | no | [0.012, 9.428] | no | lipid | no |
| pea_isolate_uht_140C_Trikusuma2019 | hexanal | ppb | 782.000 | 22.839 | 34.240 | no | no | [0.234, 2.23e+03] | yes | lipid | no |
| resconi_2023_pbma_beef_identity_benchmark | furfural | ppb | 715.220 | 7.665 | 93.308 | no | no | [1.096, 53.611] | no | sulfur | no |
| soy_isolate_40C_PratapSingh2021 | hexanal | ppb | 1.62e+03 | 0.267 | 6.08e+03 | no | no | [0.00944, 7.545] | no | lipid | no |
| thiamine_cys_glucose_120C_Bolton1994 | 2-Methyl-3-furanthiol (MFT) | ppb | 13.000 | 17.314 | 1.332 | yes | yes | [2.475, 121.095] | yes | sulfur | no |
| mp_holdout_fructose_asparagine_180C_Lin2022 | Acrylamide | ppb | 1.86e+03 | 225.440 | 8.246 | no | no | [32.233, 1.58e+03] | no | acrylamide | no |
| mp_holdout_fructose_asparagine_180C_Lin2022 | 5-Hydroxymethylfurfural (HMF) | ppb | 1.23e+04 | 1.94e+03 | 6.323 | no | no | [277.676, 1.36e+04] | yes | acrylamide | no |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | DMHF | ppb | 1.15e+03 | 21.841 | 52.799 | no | no | [2.223, 214.536] | no | acrylamide | no |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | 5-Hydroxymethylfurfural (HMF) | ppb | 5.73e+04 | 2.8e+04 | 2.048 | yes | no | [4e+03, 1.96e+05] | yes | acrylamide | no |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | DMHF | ppb | 5.89e+03 | 21.841 | 269.863 | no | no | [2.223, 214.536] | no | acrylamide | no |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | 5-Hydroxymethylfurfural (HMF) | ppb | 1.01e+05 | 2.8e+04 | 3.599 | no | no | [4e+03, 1.96e+05] | yes | acrylamide | no |
| mp_holdout_glucose_asparagine_180C_10min_Chang2021 | Acrylamide | ppb | 28.000 | 6.75e+03 | 241.056 | no | no | [965.032, 4.72e+04] | no | acrylamide | no |
| mp_holdout_glucose_asparagine_180C_30min_Chang2021 | Acrylamide | ppb | 1.46e+03 | 4.03e+03 | 2.763 | yes | no | [576.354, 2.82e+04] | yes | acrylamide | no |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | Acrylamide | ppb | 832.000 | 4.03e+03 | 4.845 | no | no | [576.354, 2.82e+04] | yes | acrylamide | no |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | 5-Hydroxymethylfurfural (HMF) | ppb | 7e+03 | 2.11e+03 | 3.320 | no | no | [301.447, 1.47e+04] | yes | acrylamide | no |
| mp_holdout_glucose_asparagine_180C_Ye2024 | Acrylamide | umol_per_mol_limiting_precursor | 140.580 | 7.04e+03 | 50.058 | no | no | - | - | acrylamide | no |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 5-Hydroxymethylfurfural (HMF) | ppb | 1.74e+04 | 1.46e+03 | 11.930 | no | no | [208.533, 1.02e+04] | no | trunk | no |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | 2-Furfurylthiol (FFT) | ppb | 7.000 | 837.456 | 119.637 | no | no | [119.737, 5.86e+03] | no | sulfur | no |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | 2-Methyl-3-furanthiol (MFT) | ppb | 3.000 | 1.4e-09 | 2.15e+09 | no | no | [2e-10, 9.78e-09] | no | sulfur | no |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | 2-Furfurylthiol (FFT) | ppb | 6.000 | 3.027 | 1.982 | yes | no | [0.433, 21.173] | yes | sulfur | no |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | 2-Methyl-3-furanthiol (MFT) | ppb | 4.000 | 1.23e-09 | 3.24e+09 | no | no | [1.76e-10, 8.62e-09] | no | sulfur | no |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | 2-Furfurylthiol (FFT) | ppb | 229.000 | 1.14e+03 | 4.965 | no | no | [162.556, 7.95e+03] | yes | sulfur | no (shared: hofmann_ribose_pH3_FFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | 2-Methyl-3-furanthiol (MFT) | ppb | 553.000 | 23.470 | 23.562 | no | no | [3.356, 164.153] | no | sulfur | no (shared: hofmann_ribose_pH3_MFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | 2-Furfurylthiol (FFT) | ppb | 12.000 | 0.090 | 133.596 | no | no | [0.013, 0.628] | no | sulfur | no (shared: hofmann_ribose_pH7_FFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | 2-Methyl-3-furanthiol (MFT) | ppb | 25.000 | 0.178 | 140.309 | no | no | [0.025, 1.246] | no | sulfur | no (shared: hofmann_ribose_pH7_MFT) |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | 2-Furfurylthiol (FFT) | ppb | 96.000 | 192.701 | 2.007 | yes | no | [27.552, 1.35e+03] | yes | sulfur | no |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | ppb | 143.000 | 548.334 | 3.835 | no | no | [78.399, 3.84e+03] | yes | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | ppb | 6.880 | 64.116 | 9.319 | no | no | [9.167, 448.435] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Furfurylthiol (FFT) | ppb | 1.280 | 615.517 | 480.873 | no | no | [88.005, 4.31e+03] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | ppb | 3.290 | 96.974 | 29.475 | no | no | [13.865, 678.249] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Furfurylthiol (FFT) | ppb | 1.460 | 541.772 | 371.077 | no | no | [77.461, 3.79e+03] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | ppb | 2.400 | 134.711 | 56.130 | no | no | [19.261, 942.187] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Furfurylthiol (FFT) | ppb | 1.680 | 367.708 | 218.874 | no | no | [52.574, 2.57e+03] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | ppb | 1.710 | 170.934 | 99.961 | no | no | [24.440, 1.2e+03] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Furfurylthiol (FFT) | ppb | 1.620 | 213.745 | 131.941 | no | no | [30.561, 1.49e+03] | no | sulfur | no |
| external_validation_bi_2020_raw_pea_hexanal | hexanal | ppb | 1.26e+03 | 0.339 | 3.72e+03 | no | no | [0.012, 9.428] | no | lipid | no |
| external_validation_bi_2020_roasted_pea_hexanal | hexanal | ppb | 324.000 | 88.598 | 3.657 | no | no | [4.054, 1.94e+03] | yes | lipid | no |
| external_validation_li_2026_spi_wg_hme_control | hexanal | ppb | 605.600 | 69.695 | 8.689 | no | no | [2.132, 2.28e+03] | yes | lipid | no |
| external_validation_liu_2023_ppi_offnote_baseline | hexanal | ppb | 1.13e+04 | 0.339 | 3.34e+04 | no | no | [0.012, 9.428] | no | lipid | no |

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
