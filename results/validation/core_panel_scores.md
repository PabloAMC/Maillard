# Core panel scorecard (kinetic core on the union panel)

pass band = 3.0x; contracts from each bundle's `validation_contract.scale_thresholds`, else the global default.

* panel: **40** benchmarks, 32 scored, 8 fully refused; rows **49**, refused rows 18
* within 3x: 8/49 (0.163); median fold 11.072, geometric mean 35.986, worst 7.56e+06
* evidence roles (core): {'predictive': 40}
* predictive benchmarks passing their contract: thiamine_cys_glucose_120C_Bolton1994; strict-ready: thiamine_cys_glucose_120C_Bolton1994
* **honest literature: 8/49 within band** (0.163), 32 benchmarks, median fold 11.072, geometric mean 35.986
* **out-of-sample: 4/40 within band** (0.100), 27 benchmarks, median fold 31.179, geometric mean 58.757
* **rows the sulfur fit read: 4/9 within band** (0.444), 5 benchmarks, median fold 3.521, geometric mean 4.072

## Per panel / role / lane

| split | key | benchmarks | rows | within band | rate | contract passes | strict-ready | median fold | geo-mean fold |
|---|---|---|---|---|---|---|---|---|---|
| panel | external_matrix | 4 | 4 | 0 | 0.000 | 0 | 0 | 1.86e+03 | 250.597 |
| panel | maillard_path_holdout | 16 | 28 | 2 | 0.071 | 0 | 0 | 26.451 | 28.216 |
| panel | trust_loop | 20 | 17 | 6 | 0.353 | 1 | 1 | 4.699 | 34.025 |
| evidence_role | predictive | 40 | 49 | 8 | 0.163 | 1 | 1 | 11.072 | 35.986 |
| lane | acrylamide | - | 12 | 2 | 0.167 | 0 | 0 | 7.285 | 21.945 |
| lane | lipid | - | 7 | 0 | 0.000 | 0 | 0 | 3.36e+03 | 430.826 |
| lane | sulfur | - | 29 | 6 | 0.207 | 0 | 0 | 10.571 | 25.194 |
| lane | trunk | - | 1 | 0 | 0.000 | 0 | 0 | 11.930 | 11.930 |

## Benchmarks

| benchmark | panel | tier | role | rows | coverage | max ratio | mean log10 | contract (ratio / log10) | status | strict | in core fit |
|---|---|---|---|---|---|---|---|---|---|---|---|
| acrylamide_spi_extrusion_130C_ACSRef3 | trust_loop | PRIMARY | predictive | 1/1 | 1.000 | 4.25e+03 | 3.628 | 1.50 / 0.200 | scale-gap | no | - |
| cml_cel_commercial_pbma_Foods2023 | trust_loop | PRIMARY | predictive | 0/2 | 0.000 | - | - | 1.80 / 0.250 | refused | no | - |
| cys_ribose_140C_Hofmann1998 | trust_loop | REFERENCE | predictive | 2/2 | 1.000 | 4.699 | 0.538 | 1.50 / 0.100 | scale-gap | no | - |
| furosine_extrusion_crossover_140C_RamirezJimenez2000 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 2.00 / 0.300 | refused | no | - |
| hofmann1998_c2c3_recombination_145C_20min_pH3 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - |
| hofmann1998_c2c3_recombination_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - |
| hofmann1998_c2c3_recombination_145C_20min_pH7 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 2/2 | 1.000 | 10.571 | 0.701 | 1.10 / 0.041 | scale-gap | no | 2 |
| hofmann1998_furan2aldehyde_h2s_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 2/2 | 1.000 | 8.457 | 0.584 | 1.10 / 0.041 | scale-gap | no | 2 |
| hofmann1998_norfuraneol_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 1/1 | 1.000 | 2.280 | 0.358 | 1.10 / 0.041 | scale-gap | no | 1 |
| hofmann1998_norfuraneol_h2s_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 2/2 | 1.000 | 9.981 | 0.773 | 1.10 / 0.041 | scale-gap | no | 2 |
| hofmann1998_xylose_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 2/2 | 1.000 | 3.684 | 0.507 | 1.10 / 0.041 | scale-gap | no | 2 |
| pea_isolate_40C_PratapSingh2021 | trust_loop | PRIMARY | predictive | 1/2 | 0.500 | 3.36e+03 | 3.526 | 2.00 / 0.120 | coverage-gap | no | - |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | PRIMARY | predictive | 1/3 | 0.333 | 34.240 | 1.535 | 2.00 / 0.120 | coverage-gap | no | - |
| resconi_2023_pbma_beef_identity_benchmark | trust_loop | SECONDARY | predictive | 1/1 | 1.000 | 7.56e+06 | 6.878 | 1.50 / 0.100 | scale-gap | no | - |
| soy_isolate_40C_PratapSingh2021 | trust_loop | PRIMARY | predictive | 1/2 | 0.500 | 6.08e+03 | 3.784 | 2.00 / 0.120 | coverage-gap | no | - |
| thiamine_cys_glucose_120C_Bolton1994 | trust_loop | PRIMARY | predictive | 1/1 | 1.000 | 1.337 | 0.126 | 3.00 / 0.480 | pass-no-ranking | yes | - |
| thiamine_cys_xylose_145C_Cerny2008 | trust_loop | PRIMARY | predictive | 0/0 | 0.000 | - | - | 1.50 / 0.100 | refused | no | - |
| mp_holdout_fructose_asparagine_180C_Lin2022 | maillard_path_holdout | REFERENCE | predictive | 2/2 | 1.000 | 8.246 | 0.859 | 1.50 / 0.100 | scale-gap | no | - |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | maillard_path_holdout | REFERENCE | predictive | 2/3 | 0.667 | 52.799 | 1.017 | 1.50 / 0.100 | coverage-gap | no | - |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | maillard_path_holdout | REFERENCE | predictive | 2/3 | 0.667 | 269.863 | 1.494 | 1.50 / 0.100 | coverage-gap | no | - |
| mp_holdout_glucose_asparagine_180C_10min_Chang2021 | maillard_path_holdout | REFERENCE | predictive | 1/1 | 1.000 | 241.056 | 2.382 | 1.50 / 0.100 | scale-gap | no | - |
| mp_holdout_glucose_asparagine_180C_30min_Chang2021 | maillard_path_holdout | REFERENCE | predictive | 1/1 | 1.000 | 2.763 | 0.441 | 1.50 / 0.100 | scale-gap | no | - |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | maillard_path_holdout | REFERENCE | predictive | 2/2 | 1.000 | 4.845 | 0.603 | 1.50 / 0.100 | scale-gap | no | - |
| mp_holdout_glucose_asparagine_180C_Ye2024 | maillard_path_holdout | REFERENCE | predictive | 1/1 | 1.000 | 50.058 | 1.699 | 1.50 / 0.100 | scale-gap | no | - |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | maillard_path_holdout | REFERENCE | predictive | 1/1 | 1.000 | 11.930 | 1.077 | 1.50 / 0.100 | scale-gap | no | - |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | maillard_path_holdout | REFERENCE | predictive | 2/2 | 1.000 | 13.375 | 1.085 | 1.10 / 0.041 | scale-gap | no | - |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | maillard_path_holdout | REFERENCE | predictive | 2/2 | 1.000 | 51.361 | 1.235 | 1.10 / 0.041 | scale-gap | no | - |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | maillard_path_holdout | REFERENCE | predictive | 2/2 | 1.000 | 24.783 | 1.123 | 1.10 / 0.041 | scale-gap | no | - |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | maillard_path_holdout | REFERENCE | predictive | 2/2 | 1.000 | 143.814 | 2.129 | 1.10 / 0.041 | scale-gap | no | - |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | maillard_path_holdout | REFERENCE | predictive | 2/2 | 1.000 | 495.811 | 1.821 | 1.50 / 0.100 | scale-gap | no | - |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | maillard_path_holdout | REFERENCE | predictive | 2/2 | 1.000 | 381.988 | 2.016 | 1.50 / 0.100 | scale-gap | no | - |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | maillard_path_holdout | REFERENCE | predictive | 2/2 | 1.000 | 233.699 | 2.049 | 1.50 / 0.100 | scale-gap | no | - |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | maillard_path_holdout | REFERENCE | predictive | 2/2 | 1.000 | 149.464 | 2.077 | 1.50 / 0.100 | scale-gap | no | - |
| external_validation_bi_2020_raw_pea_hexanal | external_matrix | PRIMARY | predictive | 1/1 | 1.000 | 3.72e+03 | 3.570 | 2.00 / 0.120 | scale-gap | no | - |
| external_validation_bi_2020_roasted_pea_hexanal | external_matrix | PRIMARY | predictive | 1/1 | 1.000 | 3.657 | 0.563 | 2.00 / 0.120 | scale-gap | no | - |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | PRIMARY | predictive | 1/4 | 0.250 | 8.689 | 0.939 | 2.00 / 0.120 | coverage-gap | no | - |
| external_validation_liu_2023_ppi_offnote_baseline | external_matrix | PRIMARY | predictive | 1/2 | 0.500 | 3.34e+04 | 4.524 | 2.00 / 0.120 | coverage-gap | no | - |

## Rows

| benchmark | compound | unit | measured | predicted | fold | within band | within contract | interval (ug/L) | measured inside | lane | in core fit |
|---|---|---|---|---|---|---|---|---|---|---|---|
| acrylamide_spi_extrusion_130C_ACSRef3 | acrylamide | ppb | 150.000 | 0.035 | 4.25e+03 | no | no | [0.00505, 0.247] | no | acrylamide | no |
| cys_ribose_140C_Hofmann1998 | 2-methyl-3-furanthiol | ppb | 342.000 | 72.778 | 4.699 | no | no | [10.406, 509.016] | yes | sulfur | no |
| cys_ribose_140C_Hofmann1998 | 2-furfurylthiol | ppb | 200.000 | 78.932 | 2.534 | yes | no | [11.285, 552.062] | yes | sulfur | no |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | 2-Furfurylthiol (FFT) | ppb | 32.000 | 338.275 | 10.571 | no | no | [48.366, 2.37e+03] | no | sulfur | yes |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | ppb | 25.000 | 10.463 | 2.389 | yes | no | [1.496, 73.180] | yes | sulfur | yes |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | 2-Furfurylthiol (FFT) | ppb | 28.000 | 236.804 | 8.457 | no | no | [33.857, 1.66e+03] | no | sulfur | yes |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | ppb | 19.000 | 33.088 | 1.741 | yes | no | [4.731, 231.421] | yes | sulfur | yes |
| hofmann1998_norfuraneol_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | ppb | 1.02e+03 | 2.32e+03 | 2.280 | yes | no | [331.245, 1.62e+04] | yes | sulfur | yes |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | 2-Furfurylthiol (FFT) | ppb | 121.000 | 1.21e+03 | 9.981 | no | no | [172.681, 8.45e+03] | no | sulfur | yes |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | ppb | 198.000 | 697.088 | 3.521 | no | no | [99.667, 4.88e+03] | yes | sulfur | yes |
| hofmann1998_xylose_cysteine_145C_20min_pH5 | 2-Furfurylthiol (FFT) | ppb | 96.000 | 268.846 | 2.800 | yes | no | [38.439, 1.88e+03] | yes | sulfur | yes |
| hofmann1998_xylose_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | ppb | 143.000 | 526.834 | 3.684 | no | no | [75.325, 3.68e+03] | yes | sulfur | yes |
| pea_isolate_40C_PratapSingh2021 | hexanal | ppb | 1.14e+03 | 0.339 | 3.36e+03 | no | no | [0.012, 9.428] | no | lipid | no |
| pea_isolate_uht_140C_Trikusuma2019 | hexanal | ppb | 782.000 | 22.839 | 34.240 | no | no | [0.234, 2.23e+03] | yes | lipid | no |
| resconi_2023_pbma_beef_identity_benchmark | furfural | ppb | 715.220 | 9.46e-05 | 7.56e+06 | no | no | [1.35e-05, 0.000662] | no | sulfur | no |
| soy_isolate_40C_PratapSingh2021 | hexanal | ppb | 1.62e+03 | 0.267 | 6.08e+03 | no | no | [0.00944, 7.545] | no | lipid | no |
| thiamine_cys_glucose_120C_Bolton1994 | 2-Methyl-3-furanthiol (MFT) | ppb | 13.000 | 17.383 | 1.337 | yes | yes | [2.485, 121.582] | yes | sulfur | no |
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
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | 2-Furfurylthiol (FFT) | ppb | 7.000 | 93.626 | 13.375 | no | no | [13.386, 654.830] | no | sulfur | no |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | 2-Methyl-3-furanthiol (MFT) | ppb | 3.000 | 33.217 | 11.072 | no | no | [4.749, 232.324] | no | sulfur | no |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | 2-Furfurylthiol (FFT) | ppb | 6.000 | 0.117 | 51.361 | no | no | [0.017, 0.817] | no | sulfur | no |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | 2-Methyl-3-furanthiol (MFT) | ppb | 4.000 | 22.940 | 5.735 | no | no | [3.280, 160.448] | yes | sulfur | no |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | 2-Furfurylthiol (FFT) | ppb | 229.000 | 1.63e+03 | 7.124 | no | no | [233.256, 1.14e+04] | no | sulfur | no (shared: hofmann_ribose_pH3_FFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | 2-Methyl-3-furanthiol (MFT) | ppb | 553.000 | 22.314 | 24.783 | no | no | [3.190, 156.067] | no | sulfur | no (shared: hofmann_ribose_pH3_MFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | 2-Furfurylthiol (FFT) | ppb | 12.000 | 0.095 | 125.955 | no | no | [0.014, 0.666] | no | sulfur | no (shared: hofmann_ribose_pH7_FFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | 2-Methyl-3-furanthiol (MFT) | ppb | 25.000 | 0.174 | 143.814 | no | no | [0.025, 1.216] | no | sulfur | no (shared: hofmann_ribose_pH7_MFT) |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | ppb | 6.880 | 60.743 | 8.829 | no | no | [8.685, 424.846] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Furfurylthiol (FFT) | ppb | 1.280 | 634.638 | 495.811 | no | no | [90.739, 4.44e+03] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | ppb | 3.290 | 92.511 | 28.119 | no | no | [13.227, 647.032] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Furfurylthiol (FFT) | ppb | 1.460 | 557.702 | 381.988 | no | no | [79.739, 3.9e+03] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | ppb | 2.400 | 128.587 | 53.578 | no | no | [18.385, 899.352] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Furfurylthiol (FFT) | ppb | 1.680 | 392.614 | 233.699 | no | no | [56.135, 2.75e+03] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | ppb | 1.710 | 163.367 | 95.536 | no | no | [23.358, 1.14e+03] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Furfurylthiol (FFT) | ppb | 1.620 | 242.131 | 149.464 | no | no | [34.619, 1.69e+03] | no | sulfur | no |
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
