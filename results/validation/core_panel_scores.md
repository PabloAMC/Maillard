# Core panel scorecard (kinetic core on the union panel)

pass band = 3.0x; contracts from each bundle's `validation_contract.scale_thresholds`, else the global default.

* panel: **42** benchmarks, 33 scored, 9 fully refused; rows **66**, refused rows 27
* within 3x: 21/66 (0.318); median fold 7.871, geometric mean 14.181, worst 2.87e+04
* evidence roles (core): {'external_holdout': 26, 'predictive': 16}
* predictive benchmarks passing their contract: hofmann1998_norfuraneol_cysteine_145C_20min_pH5; strict-ready: hofmann1998_norfuraneol_cysteine_145C_20min_pH5
* **honest literature: 21/66 within band** (0.318), 33 benchmarks, median fold 7.871, geometric mean 14.181
* **out-of-sample: 14/56 within band** (0.250), 25 benchmarks, median fold 11.802, geometric mean 19.205
* **rows the sulfur fit read: 7/10 within band** (0.700), 8 benchmarks, median fold 2.036, geometric mean 2.595

## Per panel / role / lane

| split | key | benchmarks | rows | within band | rate | contract passes | strict-ready | median fold | geo-mean fold |
|---|---|---|---|---|---|---|---|---|---|
| panel | external_matrix | 4 | 4 | 1 | 0.250 | 0 | 0 | 5.866 | 13.054 |
| panel | maillard_path_holdout | 22 | 47 | 10 | 0.213 | 0 | 0 | 14.295 | 19.411 |
| panel | trust_loop | 16 | 15 | 10 | 0.667 | 1 | 1 | 2.210 | 5.422 |
| evidence_role | external_holdout | 26 | 51 | 11 | 0.216 | 0 | 0 | 9.308 | 18.816 |
| evidence_role | predictive | 16 | 15 | 10 | 0.667 | 1 | 1 | 2.210 | 5.422 |
| lane | acrylamide | - | 12 | 4 | 0.333 | 0 | 0 | 7.697 | 20.300 |
| lane | lipid | - | 8 | 5 | 0.625 | 0 | 0 | 2.744 | 5.019 |
| lane | sulfur | - | 40 | 10 | 0.250 | 0 | 0 | 13.648 | 17.256 |
| lane | trunk | - | 6 | 2 | 0.333 | 0 | 0 | 7.968 | 7.470 |

## Benchmarks

| benchmark | panel | tier | role | rows | coverage | max ratio | mean log10 | contract (ratio / log10) | status | strict | in core fit | O2 : thiol |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| acrylamide_spi_extrusion_130C_ACSRef3 | trust_loop | PRIMARY | predictive | 1/1 | 1.000 | 9.63e+03 | 3.984 | 1.50 / 0.200 | scale-gap | no | - | continuous |
| cml_cel_commercial_pbma_Foods2023 | trust_loop | PRIMARY | predictive | 0/2 | 0.000 | - | - | 1.80 / 0.250 | refused | no | - | not_applicable |
| hofmann1998_c2c3_recombination_145C_20min_pH3 | trust_loop | PRIMARY | predictive | 1/1 | 1.000 | 13.001 | 1.114 | 1.10 / 0.041 | scale-gap | no | 1 | 1.32 |
| hofmann1998_c2c3_recombination_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 1/1 | 1.000 | 1.837 | 0.264 | 1.10 / 0.041 | scale-gap | no | 1 | 1.32 |
| hofmann1998_c2c3_recombination_145C_20min_pH7 | trust_loop | PRIMARY | predictive | 1/1 | 1.000 | 5.434 | 0.735 | 1.10 / 0.041 | scale-gap | no | 1 | 1.32 |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 0/2 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 0.27 |
| hofmann1998_furan2aldehyde_h2s_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 1/1 | 1.000 | 1.544 | 0.189 | 1.10 / 0.041 | scale-gap | no | 1 | 1.32 |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 0/2 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 0.27 |
| hofmann1998_norfuraneol_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 1/1 | 1.000 | 1.057 | 0.024 | 1.10 / 0.041 | pass-no-ranking | yes | 1 | 1.32 |
| hofmann1998_norfuraneol_h2s_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 1/1 | 1.000 | 1.958 | 0.292 | 1.10 / 0.041 | scale-gap | no | 1 | 1.32 |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 2/2 | 1.000 | 2.865 | 0.391 | 1.10 / 0.041 | scale-gap | no | 2 | 0.27 |
| pea_isolate_40C_PratapSingh2021 | trust_loop | PRIMARY | predictive | 0/2 | 0.000 | - | - | 2.00 / 0.120 | refused | no | - | not_applicable |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | PRIMARY | predictive | 4/7 | 0.571 | 2.529 | 0.286 | 2.00 / 0.120 | coverage-gap | no | - | continuous |
| resconi_2023_pbma_beef_identity_benchmark | trust_loop | SECONDARY | predictive | 1/1 | 1.000 | 15.397 | 1.187 | 1.50 / 0.100 | scale-gap | no | - | not_applicable |
| soy_isolate_40C_PratapSingh2021 | trust_loop | PRIMARY | predictive | 0/2 | 0.000 | - | - | 2.00 / 0.120 | refused | no | - | not_applicable |
| thiamine_cys_glucose_120C_Bolton1994 | trust_loop | PRIMARY | predictive | 1/1 | 1.000 | 19.907 | 1.299 | 3.00 / 0.480 | scale-gap | no | - | 2.31 |
| mp_holdout_arp_glutathione_120C_120min_pH75_Zhou2023c | maillard_path_holdout | REFERENCE | external_holdout | 3/3 | 1.000 | 161.229 | 1.243 | 1.50 / 0.100 | ranking-gap | no | - | ambiguous |
| mp_holdout_arp_glutathione_120C_30min_pH75_Zhou2023c | maillard_path_holdout | REFERENCE | external_holdout | 3/3 | 1.000 | 24.592 | 1.077 | 1.50 / 0.100 | ranking-gap | no | - | ambiguous |
| mp_holdout_arp_glutathione_120C_60min_pH75_Zhou2023c | maillard_path_holdout | REFERENCE | external_holdout | 3/3 | 1.000 | 44.933 | 1.268 | 1.50 / 0.100 | ranking-gap | no | - | ambiguous |
| mp_holdout_arp_glutathione_120C_90min_pH75_Zhou2023c | maillard_path_holdout | REFERENCE | external_holdout | 3/3 | 1.000 | 35.556 | 1.132 | 1.50 / 0.100 | ranking-gap | no | - | ambiguous |
| mp_holdout_fructose_asparagine_180C_Lin2022 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 9.149 | 0.878 | 1.50 / 0.100 | scale-gap | no | - | ambiguous |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | maillard_path_holdout | REFERENCE | external_holdout | 2/3 | 0.667 | 52.799 | 0.868 | 1.50 / 0.100 | coverage-gap | no | - | ambiguous |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | maillard_path_holdout | REFERENCE | external_holdout | 2/3 | 0.667 | 269.863 | 1.345 | 1.50 / 0.100 | coverage-gap | no | - | ambiguous |
| mp_holdout_glucose_asparagine_180C_10min_Chang2021 | maillard_path_holdout | REFERENCE | external_holdout | 1/1 | 1.000 | 218.216 | 2.339 | 1.50 / 0.100 | scale-gap | no | - | ambiguous |
| mp_holdout_glucose_asparagine_180C_30min_Chang2021 | maillard_path_holdout | REFERENCE | external_holdout | 1/1 | 1.000 | 2.522 | 0.402 | 1.50 / 0.100 | scale-gap | no | - | ambiguous |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 4.423 | 0.542 | 1.50 / 0.100 | scale-gap | no | - | ambiguous |
| mp_holdout_glucose_asparagine_180C_Ye2024 | maillard_path_holdout | REFERENCE | external_holdout | 1/1 | 1.000 | 49.824 | 1.697 | 1.50 / 0.100 | scale-gap | no | - | 0.184 mmol O2, no thiol |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | maillard_path_holdout | REFERENCE | external_holdout | 6/6 | 1.000 | 62.994 | 0.873 | 1.50 / 0.100 | ranking-gap | no | - | ambiguous |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | maillard_path_holdout | REFERENCE | external_holdout | 0/2 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 0.27 |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | maillard_path_holdout | REFERENCE | external_holdout | 0/2 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 0.27 |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 22.382 | 1.086 | 1.10 / 0.041 | scale-gap | no | - | 0.27 |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 4.91e+03 | 2.770 | 1.10 / 0.041 | scale-gap | no | - | 0.27 |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 3.611 | 0.371 | 1.10 / 0.041 | scale-gap | no | 2 | 0.27 |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 233.983 | 1.273 | 1.50 / 0.100 | scale-gap | no | - | 1.98 |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 168.865 | 1.395 | 1.50 / 0.100 | scale-gap | no | - | 1.98 |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 99.563 | 1.364 | 1.50 / 0.100 | scale-gap | no | - | 1.98 |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 59.320 | 1.322 | 1.50 / 0.100 | scale-gap | no | - | 1.98 |
| mp_holdout_xylose_cysteine_phosphate_120C_50min_pH65_Zhang2025b | maillard_path_holdout | REFERENCE | external_holdout | 4/4 | 1.000 | 2.87e+04 | 2.568 | 1.50 / 0.100 | ranking-gap | no | - | 0.44 |
| external_validation_bi_2020_raw_pea_hexanal | external_matrix | PRIMARY | external_holdout | 0/2 | 0.000 | - | - | 2.00 / 0.120 | refused | no | - | not_applicable |
| external_validation_bi_2020_roasted_pea_hexanal | external_matrix | PRIMARY | external_holdout | 1/4 | 0.250 | 2.959 | 0.471 | 2.00 / 0.120 | coverage-gap | no | - | open |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | PRIMARY | external_holdout | 3/4 | 0.750 | 365.609 | 1.331 | 2.00 / 0.120 | coverage-gap | no | - | continuous |
| external_validation_liu_2023_ppi_offnote_baseline | external_matrix | PRIMARY | external_holdout | 0/2 | 0.000 | - | - | 2.00 / 0.120 | refused | no | - | not_applicable |

## Rows

| benchmark | compound | unit | measured | predicted | fold | within band | within contract | interval (ug/L) | measured inside | lane | in core fit |
|---|---|---|---|---|---|---|---|---|---|---|---|
| acrylamide_spi_extrusion_130C_ACSRef3 | acrylamide | ppb | 150.000 | 0.016 | 9.63e+03 | no | no | [0.00223, 0.109] | no | acrylamide | no |
| hofmann1998_c2c3_recombination_145C_20min_pH3 | 2-Methyl-3-furanthiol (MFT) | ppb | 310.000 | 4.03e+03 | 13.001 | no | no | [576.226, 2.82e+04] | no | sulfur | yes |
| hofmann1998_c2c3_recombination_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | ppb | 5.36e+03 | 2.92e+03 | 1.837 | yes | no | [417.443, 2.04e+04] | yes | sulfur | yes |
| hofmann1998_c2c3_recombination_145C_20min_pH7 | 2-Methyl-3-furanthiol (MFT) | ppb | 6.23e+03 | 1.15e+03 | 5.434 | no | no | [163.920, 8.02e+03] | yes | sulfur | yes |
| hofmann1998_furan2aldehyde_h2s_145C_20min_pH5 | 2-Furfurylthiol (FFT) | ppb | 1.1e+04 | 1.7e+04 | 1.544 | yes | no | [2.43e+03, 1.19e+05] | yes | sulfur | yes |
| hofmann1998_norfuraneol_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | ppb | 1.02e+03 | 1.07e+03 | 1.057 | yes | yes | [153.577, 7.51e+03] | yes | sulfur | yes |
| hofmann1998_norfuraneol_h2s_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | ppb | 4.22e+03 | 8.27e+03 | 1.958 | yes | no | [1.18e+03, 5.78e+04] | yes | sulfur | yes |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | 2-Furfurylthiol (FFT) | ppb | 121.000 | 346.617 | 2.865 | yes | no | [49.558, 2.42e+03] | yes | sulfur | yes |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | ppb | 198.000 | 93.633 | 2.115 | yes | no | [13.387, 654.885] | yes | sulfur | yes |
| pea_isolate_uht_140C_Trikusuma2019 | hexanal | ppb | 782.000 | 353.858 | 2.210 | yes | no | [1.087, 1.15e+05] | yes | lipid | no |
| pea_isolate_uht_140C_Trikusuma2019 | 2-pentylfuran | ppb | 163.000 | 64.442 | 2.529 | yes | no | [0.198, 2.1e+04] | yes | lipid | no |
| pea_isolate_uht_140C_Trikusuma2019 | nonanal | ppb | 24.000 | 15.606 | 1.538 | yes | yes | [0.048, 5.08e+03] | yes | lipid | no |
| pea_isolate_uht_140C_Trikusuma2019 | (E,E)-2,4-decadienal | ppb | 46.900 | 75.702 | 1.614 | yes | yes | [0.233, 2.46e+04] | yes | lipid | no |
| resconi_2023_pbma_beef_identity_benchmark | furfural | ppb | 715.220 | 46.453 | 15.397 | no | no | [6.642, 324.896] | no | sulfur | no |
| thiamine_cys_glucose_120C_Bolton1994 | 2-Methyl-3-furanthiol (MFT) | ppb | 11.700 | 232.916 | 19.907 | no | no | [33.302, 1.63e+03] | no | sulfur | no |
| mp_holdout_arp_glutathione_120C_120min_pH75_Zhou2023c | 2-Methyl-3-furanthiol (MFT) | ppb | 8.100 | 1.31e+03 | 161.229 | no | no | [186.721, 9.13e+03] | no | sulfur | no |
| mp_holdout_arp_glutathione_120C_120min_pH75_Zhou2023c | 2-Furfurylthiol (FFT) | ppb | 116.000 | 129.540 | 1.117 | yes | yes | [18.521, 906.021] | yes | sulfur | no |
| mp_holdout_arp_glutathione_120C_120min_pH75_Zhou2023c | Furfural | ppb | 716.000 | 24.138 | 29.663 | no | no | [3.451, 168.822] | no | sulfur | no |
| mp_holdout_arp_glutathione_120C_30min_pH75_Zhou2023c | 2-Methyl-3-furanthiol (MFT) | ppb | 59.100 | 844.862 | 14.295 | no | no | [120.796, 5.91e+03] | no | sulfur | no |
| mp_holdout_arp_glutathione_120C_30min_pH75_Zhou2023c | 2-Furfurylthiol (FFT) | ppb | 124.000 | 25.671 | 4.830 | no | no | [3.670, 179.545] | yes | sulfur | no |
| mp_holdout_arp_glutathione_120C_30min_pH75_Zhou2023c | Furfural | ppb | 113.000 | 4.595 | 24.592 | no | no | [0.657, 32.138] | no | sulfur | no |
| mp_holdout_arp_glutathione_120C_60min_pH75_Zhou2023c | 2-Methyl-3-furanthiol (MFT) | ppb | 64.100 | 1.75e+03 | 27.283 | no | no | [250.041, 1.22e+04] | no | sulfur | no |
| mp_holdout_arp_glutathione_120C_60min_pH75_Zhou2023c | 2-Furfurylthiol (FFT) | ppb | 383.000 | 73.488 | 5.212 | no | no | [10.507, 513.985] | yes | sulfur | no |
| mp_holdout_arp_glutathione_120C_60min_pH75_Zhou2023c | Furfural | ppb | 397.000 | 8.835 | 44.933 | no | no | [1.263, 61.796] | no | sulfur | no |
| mp_holdout_arp_glutathione_120C_90min_pH75_Zhou2023c | 2-Methyl-3-furanthiol (MFT) | ppb | 57.100 | 1.65e+03 | 28.811 | no | no | [235.215, 1.15e+04] | no | sulfur | no |
| mp_holdout_arp_glutathione_120C_90min_pH75_Zhou2023c | 2-Furfurylthiol (FFT) | ppb | 270.000 | 110.837 | 2.436 | yes | no | [15.847, 775.207] | yes | sulfur | no |
| mp_holdout_arp_glutathione_120C_90min_pH75_Zhou2023c | Furfural | ppb | 575.000 | 16.172 | 35.556 | no | no | [2.312, 113.106] | no | sulfur | no |
| mp_holdout_fructose_asparagine_180C_Lin2022 | Acrylamide | ppb | 1.86e+03 | 203.194 | 9.149 | no | no | [29.052, 1.42e+03] | no | acrylamide | no |
| mp_holdout_fructose_asparagine_180C_Lin2022 | 5-Hydroxymethylfurfural (HMF) | ppb | 1.23e+04 | 1.97e+03 | 6.246 | no | no | [281.109, 1.38e+04] | yes | acrylamide | no |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | DMHF | ppb | 1.15e+03 | 21.841 | 52.799 | no | no | [2.223, 214.535] | no | acrylamide | no |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | 5-Hydroxymethylfurfural (HMF) | ppb | 5.73e+04 | 5.54e+04 | 1.033 | yes | yes | [7.92e+03, 3.88e+05] | yes | acrylamide | no |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | DMHF | ppb | 5.89e+03 | 21.841 | 269.863 | no | no | [2.223, 214.535] | no | acrylamide | no |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | 5-Hydroxymethylfurfural (HMF) | ppb | 1.01e+05 | 5.54e+04 | 1.816 | yes | no | [7.92e+03, 3.88e+05] | yes | acrylamide | no |
| mp_holdout_glucose_asparagine_180C_10min_Chang2021 | Acrylamide | ppb | 28.000 | 6.11e+03 | 218.216 | no | no | [873.594, 4.27e+04] | no | acrylamide | no |
| mp_holdout_glucose_asparagine_180C_30min_Chang2021 | Acrylamide | ppb | 1.46e+03 | 3.68e+03 | 2.522 | yes | no | [526.139, 2.57e+04] | yes | acrylamide | no |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | Acrylamide | ppb | 832.000 | 3.68e+03 | 4.423 | no | no | [526.139, 2.57e+04] | yes | acrylamide | no |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | 5-Hydroxymethylfurfural (HMF) | ppb | 7e+03 | 2.55e+03 | 2.744 | yes | no | [364.770, 1.78e+04] | yes | acrylamide | no |
| mp_holdout_glucose_asparagine_180C_Ye2024 | Acrylamide | umol_per_mol_limiting_precursor | 140.580 | 7e+03 | 49.824 | no | no | - | - | acrylamide | no |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 5-Hydroxymethylfurfural (HMF) | ppb | 1.74e+04 | 1.87e+03 | 9.308 | no | no | [267.288, 1.31e+04] | no | trunk | no |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 3-deoxyglucosone | ppb | 5.22e+04 | 6.64e+04 | 1.271 | yes | yes | [9.49e+03, 4.64e+05] | yes | trunk | no |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 3,4-dideoxyglucosone | ppb | 5.55e+04 | 8.37e+03 | 6.628 | no | no | [1.2e+03, 5.86e+04] | yes | trunk | no |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | glucosone | ppb | 7.5e+03 | 119.059 | 62.994 | no | no | [17.023, 832.715] | no | trunk | no |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | glyoxal | ppb | 5.6e+03 | 160.794 | 34.827 | no | no | [22.990, 1.12e+03] | no | trunk | no |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | methylglyoxal | ppb | 2.6e+03 | 2.57e+03 | 1.010 | yes | yes | [367.997, 1.8e+04] | yes | trunk | no |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | 2-Furfurylthiol (FFT) | ppb | 229.000 | 1.52e+03 | 6.652 | no | no | [217.793, 1.07e+04] | yes | sulfur | no (shared: hofmann_ribose_pH3_FFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | 2-Methyl-3-furanthiol (MFT) | ppb | 553.000 | 24.707 | 22.382 | no | no | [3.533, 172.806] | no | sulfur | no (shared: hofmann_ribose_pH3_MFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | 2-Furfurylthiol (FFT) | ppb | 12.000 | 0.170 | 70.756 | no | no | [0.024, 1.186] | no | sulfur | no (shared: hofmann_ribose_pH7_FFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | 2-Methyl-3-furanthiol (MFT) | ppb | 25.000 | 0.0051 | 4.91e+03 | no | no | [0.000729, 0.036] | no | sulfur | no (shared: hofmann_ribose_pH7_MFT) |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | 2-Furfurylthiol (FFT) | ppb | 96.000 | 346.617 | 3.611 | no | no | [49.558, 2.42e+03] | yes | sulfur | yes |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | 2-Methyl-3-furanthiol (MFT) | ppb | 143.000 | 93.633 | 1.527 | yes | no | [13.387, 654.885] | yes | sulfur | yes |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | ppb | 6.880 | 10.334 | 1.502 | yes | no | [1.478, 72.281] | yes | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Furfurylthiol (FFT) | ppb | 1.280 | 299.499 | 233.983 | no | no | [42.821, 2.09e+03] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | ppb | 3.290 | 12.015 | 3.652 | no | no | [1.718, 84.036] | yes | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Furfurylthiol (FFT) | ppb | 1.460 | 246.542 | 168.865 | no | no | [35.250, 1.72e+03] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | ppb | 2.400 | 12.859 | 5.358 | no | no | [1.839, 89.941] | yes | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Furfurylthiol (FFT) | ppb | 1.680 | 167.266 | 99.563 | no | no | [23.915, 1.17e+03] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | ppb | 1.710 | 12.681 | 7.416 | no | no | [1.813, 88.695] | no | sulfur | no |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Furfurylthiol (FFT) | ppb | 1.620 | 96.098 | 59.320 | no | no | [13.740, 672.122] | no | sulfur | no |
| mp_holdout_xylose_cysteine_phosphate_120C_50min_pH65_Zhang2025b | 2-Methyl-3-furanthiol (MFT) | ppb | 22.630 | 0.510 | 44.392 | no | no | [0.073, 3.565] | no | sulfur | no |
| mp_holdout_xylose_cysteine_phosphate_120C_50min_pH65_Zhang2025b | 2-Furfurylthiol (FFT) | ppb | 0.290 | 2.414 | 8.326 | no | no | [0.345, 16.887] | no | sulfur | no |
| mp_holdout_xylose_cysteine_phosphate_120C_50min_pH65_Zhang2025b | Bis(2-methyl-3-furyl) disulfide | ppb | 18.250 | 0.000636 | 2.87e+04 | no | no | [9.1e-05, 0.00445] | no | sulfur | no |
| mp_holdout_xylose_cysteine_phosphate_120C_50min_pH65_Zhang2025b | 2-Acetylthiazole | ppb | 114.130 | 0.065 | 1.76e+03 | no | no | [0.00929, 0.454] | no | sulfur | no |
| external_validation_bi_2020_roasted_pea_hexanal | hexanal | ppb | 324.000 | 109.490 | 2.959 | yes | no | [4.569, 2.62e+03] | yes | lipid | no |
| external_validation_li_2026_spi_wg_hme_control | 2-pentylfuran | ppb | 5.63e+03 | 15.387 | 365.609 | no | no | [0.138, 1.71e+03] | no | lipid | no |
| external_validation_li_2026_spi_wg_hme_control | hexanal | ppb | 605.600 | 70.264 | 8.619 | no | no | [0.631, 7.82e+03] | yes | lipid | no |
| external_validation_li_2026_spi_wg_hme_control | nonanal | ppb | 74.370 | 23.885 | 3.114 | no | no | [0.215, 2.65e+03] | yes | lipid | no |

## Refused rows

| benchmark | panel | compound | reason |
|---|---|---|---|
| cml_cel_commercial_pbma_Foods2023 | trust_loop | Nε-(Carboxymethyl)lysine (CML) | GLYCATION TARGETS 'Nε-(Carboxymethyl)lysine (CML)' (wave B20) run on the trunk lane only: the acrylamide lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + ami |
| cml_cel_commercial_pbma_Foods2023 | trust_loop | Nε-(Carboxyethyl)lysine (CEL) | GLYCATION TARGETS 'Nε-(Carboxyethyl)lysine (CEL)' (wave B20) run on the trunk lane only: the acrylamide lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + amin |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| pea_isolate_40C_PratapSingh2021 | trust_loop | hexanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's source says so, recorded as closure = 'never heated' (its timed block is the headspace in |
| pea_isolate_40C_PratapSingh2021 | trust_loop | 2-pentylfuran | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's source says so, recorded as closure = 'never heated' (its timed block is the headspace in |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | 2,5-dimethylpyrazine | PYRAZINE TARGETS '2,5-dimethylpyrazine' (wave B18) run on the trunk lane only: the lipid lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + amine pot that reso |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | methional | METHIONINE CHAIN TARGETS 'methional' (wave B22) run on the trunk lane only: the lipid lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + amine pot that resolve |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | 2-acetyl-1-pyrroline | 2-ACETYL-1-PYRROLINE TARGETS '2-acetyl-1-pyrroline' (wave B24) run on the trunk lane only: the lipid lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + amine p |
| soy_isolate_40C_PratapSingh2021 | trust_loop | hexanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's source says so, recorded as closure = 'never heated' (its timed block is the headspace in |
| soy_isolate_40C_PratapSingh2021 | trust_loop | 2-pentylfuran | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's source says so, recorded as closure = 'never heated' (its timed block is the headspace in |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | maillard_path_holdout | Furfural | LANE CONFLICT: this request needs both the acrylamide and sulfur lanes at once. They do not compose -- the acrylamide network deliberately omits every sulfur step (acrylamide.OUT_OF_SCOPE), because composing them would s |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | maillard_path_holdout | Furfural | LANE CONFLICT: this request needs both the acrylamide and sulfur lanes at once. They do not compose -- the acrylamide network deliberately omits every sulfur step (acrylamide.OUT_OF_SCOPE), because composing them would s |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| external_validation_bi_2020_raw_pea_hexanal | external_matrix | hexanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's source says so, recorded as closure = 'never heated' (its timed block is the headspace in |
| external_validation_bi_2020_raw_pea_hexanal | external_matrix | nonanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's source says so, recorded as closure = 'never heated' (its timed block is the headspace in |
| external_validation_bi_2020_roasted_pea_hexanal | external_matrix | furaneol | THIS POT CHARGES NO PRECURSOR THAT COULD MAKE 'furaneol'. The charge declares only a matrix/lipid carrier, which is not a precursor: it resolves to a hydroperoxide pool for the lipid lane and charges NOTHING on the trunk |
| external_validation_bi_2020_roasted_pea_hexanal | external_matrix | 2,5-dimethylpyrazine | PYRAZINE TARGETS '2,5-dimethylpyrazine' (wave B18) run on the trunk lane only: the lipid lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + amine pot that reso |
| external_validation_bi_2020_roasted_pea_hexanal | external_matrix | furfural | THIS POT CHARGES NO PRECURSOR THAT COULD MAKE 'furfural'. The charge declares only a matrix/lipid carrier, which is not a precursor: it resolves to a hydroperoxide pool for the lipid lane and charges NOTHING on the trunk |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | 1-hexanol | UNREPRESENTED TARGETS: 1-hexanol -- The lipid lane exists and forms the SIX products Frankel 1989 measured, but 1-hexanol is not one of them and NO aldehyde-reduction step is measured anywhere in the corpus -- in a therm |
| external_validation_liu_2023_ppi_offnote_baseline | external_matrix | hexanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's source says so, recorded as closure = 'never heated' (its timed block is the headspace in |
| external_validation_liu_2023_ppi_offnote_baseline | external_matrix | nonanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's source says so, recorded as closure = 'never heated' (its timed block is the headspace in |

## Bundles kept off the scored panel

- pea_isolate_ribose_cysteine_100C_45min_Internal2026 (trust_loop): synthetic snapshot (_Internal2026): legacy-model output, not a measurement
- soy_isolate_ribose_cysteine_100C_45min_Internal2026 (trust_loop): synthetic snapshot (_Internal2026): legacy-model output, not a measurement
