# Core panel scorecard (kinetic core on the union panel)

pass band = 3.0x; contracts from each bundle's `validation_contract.scale_thresholds`, else the global default.

* panel: **37** benchmarks, 23 scored, 14 fully refused; rows **44**, refused rows 25
* within 3x: 9/44 (0.205); median fold 10.622, geometric mean 18.798, worst 9.63e+03
* evidence roles (core): {'external_holdout': 21, 'predictive': 16}
* predictive benchmarks passing their contract: NONE; strict-ready: NONE
* **honest literature: 9/44 within band** (0.205), 23 benchmarks, median fold 10.622, geometric mean 18.798
* **out-of-sample: 8/43 within band** (0.186), 22 benchmarks, median fold 11.930, geometric mean 19.736
* **rows the sulfur fit read: 1/1 within band** (1.000), 1 benchmarks, median fold 2.316, geometric mean 2.316

## Per panel / role / lane

| split | key | benchmarks | rows | within band | rate | contract passes | strict-ready | median fold | geo-mean fold |
|---|---|---|---|---|---|---|---|---|---|
| panel | external_matrix | 4 | 4 | 0 | 0.000 | 0 | 0 | 6.173 | 13.711 |
| panel | maillard_path_holdout | 17 | 31 | 5 | 0.161 | 0 | 0 | 29.450 | 21.897 |
| panel | trust_loop | 16 | 9 | 4 | 0.444 | 0 | 0 | 3.690 | 12.787 |
| evidence_role | external_holdout | 21 | 35 | 5 | 0.143 | 0 | 0 | 23.562 | 20.756 |
| evidence_role | predictive | 16 | 9 | 4 | 0.444 | 0 | 0 | 3.690 | 12.787 |
| lane | acrylamide | - | 12 | 2 | 0.167 | 0 | 0 | 7.728 | 23.060 |
| lane | lipid | - | 7 | 3 | 0.429 | 0 | 0 | 3.042 | 6.071 |
| lane | sulfur | - | 19 | 2 | 0.105 | 0 | 0 | 29.450 | 30.281 |
| lane | trunk | - | 6 | 2 | 0.333 | 0 | 0 | 22.182 | 10.320 |

## Benchmarks

| benchmark | panel | tier | role | rows | coverage | max ratio | mean log10 | contract (ratio / log10) | status | strict | in core fit | O2 : thiol |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| acrylamide_spi_extrusion_130C_ACSRef3 | trust_loop | PRIMARY | predictive | 1/1 | 1.000 | 9.63e+03 | 3.984 | 1.50 / 0.200 | scale-gap | no | - | continuous |
| cml_cel_commercial_pbma_Foods2023 | trust_loop | PRIMARY | predictive | 0/2 | 0.000 | - | - | 1.80 / 0.250 | refused | no | - | not_applicable |
| hofmann1998_c2c3_recombination_145C_20min_pH3 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 1.32 |
| hofmann1998_c2c3_recombination_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 1.32 |
| hofmann1998_c2c3_recombination_145C_20min_pH7 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 1.32 |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 0/2 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 0.27 |
| hofmann1998_furan2aldehyde_h2s_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 1.32 |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 0/2 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 0.27 |
| hofmann1998_norfuraneol_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 1/1 | 1.000 | 2.316 | 0.365 | 1.10 / 0.041 | scale-gap | no | 1 | 1.32 |
| hofmann1998_norfuraneol_h2s_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 1.32 |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 2/2 | 1.000 | 6.865 | 0.702 | 1.10 / 0.041 | scale-gap | no | - | 0.27 |
| pea_isolate_40C_PratapSingh2021 | trust_loop | PRIMARY | predictive | 0/2 | 0.000 | - | - | 2.00 / 0.120 | refused | no | - | not_applicable |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | PRIMARY | predictive | 3/3 | 1.000 | 2.529 | 0.311 | 2.00 / 0.120 | scale-gap | no | - | continuous |
| resconi_2023_pbma_beef_identity_benchmark | trust_loop | SECONDARY | predictive | 1/1 | 1.000 | 93.308 | 1.970 | 1.50 / 0.100 | scale-gap | no | - | not_applicable |
| soy_isolate_40C_PratapSingh2021 | trust_loop | PRIMARY | predictive | 0/2 | 0.000 | - | - | 2.00 / 0.120 | refused | no | - | not_applicable |
| thiamine_cys_glucose_120C_Bolton1994 | trust_loop | PRIMARY | predictive | 1/1 | 1.000 | 20.168 | 1.305 | 3.00 / 0.480 | scale-gap | no | - | 2.31 |
| mp_holdout_fructose_asparagine_180C_Lin2022 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 9.149 | 0.881 | 1.50 / 0.100 | scale-gap | no | - | ambiguous |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | maillard_path_holdout | REFERENCE | external_holdout | 2/3 | 0.667 | 52.799 | 1.017 | 1.50 / 0.100 | coverage-gap | no | - | ambiguous |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | maillard_path_holdout | REFERENCE | external_holdout | 2/3 | 0.667 | 269.863 | 1.494 | 1.50 / 0.100 | coverage-gap | no | - | ambiguous |
| mp_holdout_glucose_asparagine_180C_10min_Chang2021 | maillard_path_holdout | REFERENCE | external_holdout | 1/1 | 1.000 | 218.216 | 2.339 | 1.50 / 0.100 | scale-gap | no | - | ambiguous |
| mp_holdout_glucose_asparagine_180C_30min_Chang2021 | maillard_path_holdout | REFERENCE | external_holdout | 1/1 | 1.000 | 2.522 | 0.402 | 1.50 / 0.100 | scale-gap | no | - | ambiguous |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 4.423 | 0.575 | 1.50 / 0.100 | scale-gap | no | - | ambiguous |
| mp_holdout_glucose_asparagine_180C_Ye2024 | maillard_path_holdout | REFERENCE | external_holdout | 1/1 | 1.000 | 49.824 | 1.697 | 1.50 / 0.100 | scale-gap | no | - | 0.184 mmol O2, no thiol |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | maillard_path_holdout | REFERENCE | external_holdout | 6/6 | 1.000 | 62.994 | 1.014 | 1.50 / 0.100 | ranking-gap | no | - | ambiguous |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | maillard_path_holdout | REFERENCE | external_holdout | 0/2 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 0.27 |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | maillard_path_holdout | REFERENCE | external_holdout | 0/2 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 0.27 |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 23.562 | 1.033 | 1.10 / 0.041 | scale-gap | no | - | 0.27 |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 140.337 | 2.137 | 1.10 / 0.041 | scale-gap | no | - | 0.27 |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 3.821 | 0.441 | 1.10 / 0.041 | scale-gap | no | - | 0.27 |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 476.762 | 1.824 | 1.50 / 0.100 | scale-gap | no | - | 1.98 |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 366.778 | 2.017 | 1.50 / 0.100 | scale-gap | no | - | 1.98 |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 216.336 | 2.042 | 1.50 / 0.100 | scale-gap | no | - | 1.98 |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 130.704 | 2.058 | 1.50 / 0.100 | scale-gap | no | - | 1.98 |
| external_validation_bi_2020_raw_pea_hexanal | external_matrix | PRIMARY | external_holdout | 0/1 | 0.000 | - | - | 2.00 / 0.120 | refused | no | - | not_applicable |
| external_validation_bi_2020_roasted_pea_hexanal | external_matrix | PRIMARY | external_holdout | 1/1 | 1.000 | 3.657 | 0.563 | 2.00 / 0.120 | scale-gap | no | - | open |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | PRIMARY | external_holdout | 3/4 | 0.750 | 365.609 | 1.328 | 2.00 / 0.120 | coverage-gap | no | - | continuous |
| external_validation_liu_2023_ppi_offnote_baseline | external_matrix | PRIMARY | external_holdout | 0/2 | 0.000 | - | - | 2.00 / 0.120 | refused | no | - | not_applicable |

## Rows

| benchmark | compound | unit | measured | predicted | fold | within band | within contract | interval (ug/L) | measured inside | lane | in core fit |
|---|---|---|---|---|---|---|---|---|---|---|---|
| acrylamide_spi_extrusion_130C_ACSRef3 | acrylamide | ppb | 150.000 | 0.016 | 9.63e+03 | no | no | [0.00223, 0.109] | no | acrylamide | no |
| hofmann1998_norfuraneol_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | ppb | 1.02e+03 | 2.35e+03 | 2.316 | yes | no | [336.411, 1.65e+04] | yes | sulfur | yes |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | 2-Furfurylthiol (FFT) | ppb | 121.000 | 830.643 | 6.865 | no | no | [118.763, 5.81e+03] | yes | sulfur | no |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | ppb | 198.000 | 730.638 | 3.690 | no | no | [104.464, 5.11e+03] | yes | sulfur | no |
| pea_isolate_uht_140C_Trikusuma2019 | hexanal | ppb | 782.000 | 353.839 | 2.210 | yes | no | [3.623, 3.46e+04] | yes | lipid | no |
| pea_isolate_uht_140C_Trikusuma2019 | 2-pentylfuran | ppb | 163.000 | 64.442 | 2.529 | yes | no | [0.660, 6.29e+03] | yes | lipid | no |
| pea_isolate_uht_140C_Trikusuma2019 | nonanal | ppb | 24.000 | 15.606 | 1.538 | yes | yes | [0.160, 1.52e+03] | yes | lipid | no |
| resconi_2023_pbma_beef_identity_benchmark | furfural | ppb | 715.220 | 7.665 | 93.308 | no | no | [1.096, 53.611] | no | sulfur | no |
| thiamine_cys_glucose_120C_Bolton1994 | 2-Methyl-3-furanthiol (MFT) | ppb | 11.700 | 235.970 | 20.168 | no | no | [33.738, 1.65e+03] | no | sulfur | no |
| mp_holdout_fructose_asparagine_180C_Lin2022 | Acrylamide | ppb | 1.86e+03 | 203.194 | 9.149 | no | no | [29.052, 1.42e+03] | no | acrylamide | no |
| mp_holdout_fructose_asparagine_180C_Lin2022 | 5-Hydroxymethylfurfural (HMF) | ppb | 1.23e+04 | 1.95e+03 | 6.308 | no | no | [278.358, 1.36e+04] | yes | acrylamide | no |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | DMHF | ppb | 1.15e+03 | 21.841 | 52.799 | no | no | [2.223, 214.536] | no | acrylamide | no |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | 5-Hydroxymethylfurfural (HMF) | ppb | 5.73e+04 | 2.8e+04 | 2.048 | yes | no | [4e+03, 1.96e+05] | yes | acrylamide | no |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | DMHF | ppb | 5.89e+03 | 21.841 | 269.863 | no | no | [2.223, 214.536] | no | acrylamide | no |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | 5-Hydroxymethylfurfural (HMF) | ppb | 1.01e+05 | 2.8e+04 | 3.599 | no | no | [4e+03, 1.96e+05] | yes | acrylamide | no |
| mp_holdout_glucose_asparagine_180C_10min_Chang2021 | Acrylamide | ppb | 28.000 | 6.11e+03 | 218.216 | no | no | [873.594, 4.27e+04] | no | acrylamide | no |
| mp_holdout_glucose_asparagine_180C_30min_Chang2021 | Acrylamide | ppb | 1.46e+03 | 3.68e+03 | 2.522 | yes | no | [526.139, 2.57e+04] | yes | acrylamide | no |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | Acrylamide | ppb | 832.000 | 3.68e+03 | 4.423 | no | no | [526.139, 2.57e+04] | yes | acrylamide | no |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | 5-Hydroxymethylfurfural (HMF) | ppb | 7e+03 | 2.19e+03 | 3.194 | no | no | [313.351, 1.53e+04] | yes | acrylamide | no |
| mp_holdout_glucose_asparagine_180C_Ye2024 | Acrylamide | umol_per_mol_limiting_precursor | 140.580 | 7e+03 | 49.824 | no | no | - | - | acrylamide | no |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 5-Hydroxymethylfurfural (HMF) | ppb | 1.74e+04 | 1.46e+03 | 11.930 | no | no | [208.533, 1.02e+04] | no | trunk | no |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 3-deoxyglucosone | ppb | 5.22e+04 | 4.7e+04 | 1.111 | yes | yes | [6.72e+03, 3.29e+05] | yes | trunk | no |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 3,4-dideoxyglucosone | ppb | 5.55e+04 | 1.71e+03 | 32.435 | no | no | [244.652, 1.2e+04] | no | trunk | no |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | glucosone | ppb | 7.5e+03 | 119.059 | 62.994 | no | no | [17.023, 832.715] | no | trunk | no |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | glyoxal | ppb | 5.6e+03 | 160.794 | 34.827 | no | no | [22.990, 1.12e+03] | no | trunk | no |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | methylglyoxal | ppb | 2.6e+03 | 2.03e+03 | 1.281 | yes | yes | [290.205, 1.42e+04] | yes | trunk | no |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | 2-Furfurylthiol (FFT) | ppb | 229.000 | 1.13e+03 | 4.932 | no | no | [161.492, 7.9e+03] | yes | sulfur | no (shared: hofmann_ribose_pH3_FFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | 2-Methyl-3-furanthiol (MFT) | ppb | 553.000 | 23.470 | 23.562 | no | no | [3.356, 164.149] | no | sulfur | no (shared: hofmann_ribose_pH3_MFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | 2-Furfurylthiol (FFT) | ppb | 12.000 | 0.090 | 133.631 | no | no | [0.013, 0.628] | no | sulfur | no (shared: hofmann_ribose_pH7_FFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | 2-Methyl-3-furanthiol (MFT) | ppb | 25.000 | 0.178 | 140.337 | no | no | [0.025, 1.246] | no | sulfur | no (shared: hofmann_ribose_pH7_MFT) |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | 2-Furfurylthiol (FFT) | ppb | 96.000 | 191.058 | 1.990 | yes | no | [27.317, 1.34e+03] | yes | sulfur | no |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | ppb | 143.000 | 546.398 | 3.821 | no | no | [78.122, 3.82e+03] | yes | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | ppb | 6.880 | 64.086 | 9.315 | no | no | [9.163, 448.229] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Furfurylthiol (FFT) | ppb | 1.280 | 610.255 | 476.762 | no | no | [87.252, 4.27e+03] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | ppb | 3.290 | 96.892 | 29.450 | no | no | [13.853, 677.673] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Furfurylthiol (FFT) | ppb | 1.460 | 535.496 | 366.778 | no | no | [76.564, 3.75e+03] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | ppb | 2.400 | 134.513 | 56.047 | no | no | [19.232, 940.802] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Furfurylthiol (FFT) | ppb | 1.680 | 363.445 | 216.336 | no | no | [51.964, 2.54e+03] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | ppb | 1.710 | 170.534 | 99.727 | no | no | [24.382, 1.19e+03] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Furfurylthiol (FFT) | ppb | 1.620 | 211.740 | 130.704 | no | no | [30.274, 1.48e+03] | no | sulfur | no |
| external_validation_bi_2020_roasted_pea_hexanal | hexanal | ppb | 324.000 | 88.598 | 3.657 | no | no | [4.054, 1.94e+03] | yes | lipid | no |
| external_validation_li_2026_spi_wg_hme_control | 2-pentylfuran | ppb | 5.63e+03 | 15.387 | 365.609 | no | no | [0.471, 502.987] | no | lipid | no |
| external_validation_li_2026_spi_wg_hme_control | hexanal | ppb | 605.600 | 69.695 | 8.689 | no | no | [2.132, 2.28e+03] | yes | lipid | no |
| external_validation_li_2026_spi_wg_hme_control | nonanal | ppb | 72.660 | 23.885 | 3.042 | no | no | [0.731, 780.763] | yes | lipid | no |

## Refused rows

| benchmark | panel | compound | reason |
|---|---|---|---|
| cml_cel_commercial_pbma_Foods2023 | trust_loop | Nε-(Carboxymethyl)lysine (CML) | GLYCATION TARGETS 'Nε-(Carboxymethyl)lysine (CML)' (wave B20) run on the trunk lane only: the acrylamide lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + ami |
| cml_cel_commercial_pbma_Foods2023 | trust_loop | Nε-(Carboxyethyl)lysine (CEL) | GLYCATION TARGETS 'Nε-(Carboxyethyl)lysine (CEL)' (wave B20) run on the trunk lane only: the acrylamide lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + amin |
| hofmann1998_c2c3_recombination_145C_20min_pH3 | trust_loop | 2-Methyl-3-furanthiol (MFT) | UNMAPPED PRECURSORS 'Hydroxyacetaldehyde', 'Mercapto-2-propanone': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge. |
| hofmann1998_c2c3_recombination_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | UNMAPPED PRECURSORS 'Hydroxyacetaldehyde', 'Mercapto-2-propanone': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge. |
| hofmann1998_c2c3_recombination_145C_20min_pH7 | trust_loop | 2-Methyl-3-furanthiol (MFT) | UNMAPPED PRECURSORS 'Hydroxyacetaldehyde', 'Mercapto-2-propanone': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge. |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_furan2aldehyde_h2s_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | UNMAPPED PRECURSORS 'Furan-2-aldehyde', 'Hydrogen sulfide': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge. |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_norfuraneol_h2s_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | UNMAPPED PRECURSORS 'Hydrogen sulfide': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge. |
| pea_isolate_40C_PratapSingh2021 | trust_loop | hexanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's own vessel says so (closure = 'no cook'), and the physics agrees: over this thermal progr |
| pea_isolate_40C_PratapSingh2021 | trust_loop | 2-pentylfuran | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's own vessel says so (closure = 'no cook'), and the physics agrees: over this thermal progr |
| soy_isolate_40C_PratapSingh2021 | trust_loop | hexanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's own vessel says so (closure = 'no cook'), and the physics agrees: over this thermal progr |
| soy_isolate_40C_PratapSingh2021 | trust_loop | 2-pentylfuran | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's own vessel says so (closure = 'no cook'), and the physics agrees: over this thermal progr |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | maillard_path_holdout | Furfural | LANE CONFLICT: this request needs both the acrylamide and sulfur lanes at once. They do not compose -- the acrylamide network deliberately omits every sulfur step (acrylamide.OUT_OF_SCOPE), because composing them would s |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | maillard_path_holdout | Furfural | LANE CONFLICT: this request needs both the acrylamide and sulfur lanes at once. They do not compose -- the acrylamide network deliberately omits every sulfur step (acrylamide.OUT_OF_SCOPE), because composing them would s |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| external_validation_bi_2020_raw_pea_hexanal | external_matrix | hexanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's own vessel says so (closure = 'no cook'), and the physics agrees: over this thermal progr |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | 1-hexanol | UNREPRESENTED TARGETS: 1-hexanol -- The lipid lane exists and forms the SIX products Frankel 1989 measured, but 1-hexanol is not one of them and NO aldehyde-reduction step is measured anywhere in the corpus -- in a therm |
| external_validation_liu_2023_ppi_offnote_baseline | external_matrix | hexanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's own vessel says so (closure = 'no cook'), and the physics agrees: over this thermal progr |
| external_validation_liu_2023_ppi_offnote_baseline | external_matrix | nonanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's own vessel says so (closure = 'no cook'), and the physics agrees: over this thermal progr |

## Bundles kept off the scored panel

- pea_isolate_ribose_cysteine_100C_45min_Internal2026 (trust_loop): synthetic snapshot (_Internal2026): legacy-model output, not a measurement
- soy_isolate_ribose_cysteine_100C_45min_Internal2026 (trust_loop): synthetic snapshot (_Internal2026): legacy-model output, not a measurement
