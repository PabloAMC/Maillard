# Core panel scorecard (kinetic core on the union panel)

pass band = 3.0x; contracts from each bundle's `validation_contract.scale_thresholds`, else the global default.

* panel: **42** benchmarks, 10 scored, 32 fully refused; rows **23**, refused rows 70
* within 3x: 10/23 (0.435); median fold 4.423, geometric mean 7.600, worst 365.609
* **the denominator fell by refusal, not by any prediction improving**: since 2026-11-05 the 27 rows whose absolute level the thiol lane does not ship (H0 / EX-059, D2 branch (i)) are refused and out of every count here
* evidence roles (core): {'external_holdout': 26, 'predictive': 16}
* predictive benchmarks passing their contract: NONE; strict-ready: NONE
* **honest literature: 10/23 within band** (0.435), 10 benchmarks, median fold 4.423, geometric mean 7.600
* **out-of-sample: 10/23 within band** (0.435), 10 benchmarks, median fold 4.423, geometric mean 7.600
* **rows the sulfur fit read: 0/0 within band** (-), 0 benchmarks, median fold -, geometric mean -

## Per panel / role / lane

| split | key | benchmarks | rows | within band | rate | contract passes | strict-ready | median fold | geo-mean fold |
|---|---|---|---|---|---|---|---|---|---|
| panel | external_matrix | 4 | 3 | 0 | 0.000 | 0 | 0 | 8.619 | 21.408 |
| panel | maillard_path_holdout | 22 | 16 | 6 | 0.375 | 0 | 0 | 5.757 | 8.816 |
| panel | trust_loop | 16 | 4 | 4 | 1.000 | 0 | 0 | 1.912 | 1.930 |
| evidence_role | external_holdout | 26 | 19 | 6 | 0.316 | 0 | 0 | 6.690 | 10.142 |
| evidence_role | predictive | 16 | 4 | 4 | 1.000 | 0 | 0 | 1.912 | 1.930 |
| lane | acrylamide | - | 10 | 4 | 0.400 | 0 | 0 | 6.786 | 12.333 |
| lane | lipid | - | 7 | 4 | 0.571 | 0 | 0 | 2.529 | 5.413 |
| lane | trunk | - | 6 | 2 | 0.333 | 0 | 0 | 5.757 | 5.039 |

## Benchmarks

| benchmark | panel | tier | role | rows | coverage | max ratio | mean log10 | contract (ratio / log10) | status | strict | in core fit | O2 : thiol |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| acrylamide_spi_extrusion_130C_ACSRef3 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.50 / 0.200 | refused | no | - | continuous |
| cml_cel_commercial_pbma_Foods2023 | trust_loop | PRIMARY | predictive | 0/2 | 0.000 | - | - | 1.80 / 0.250 | refused | no | - | not_applicable |
| hofmann1998_c2c3_recombination_145C_20min_pH3 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 1.32 |
| hofmann1998_c2c3_recombination_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 1.32 |
| hofmann1998_c2c3_recombination_145C_20min_pH7 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 1.32 |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 0/2 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 0.27 |
| hofmann1998_furan2aldehyde_h2s_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 1.32 |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 0/2 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 0.27 |
| hofmann1998_norfuraneol_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 1.32 |
| hofmann1998_norfuraneol_h2s_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 1.32 |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | trust_loop | PRIMARY | predictive | 0/2 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 0.27 |
| pea_isolate_40C_PratapSingh2021 | trust_loop | PRIMARY | predictive | 0/2 | 0.000 | - | - | 2.00 / 0.120 | refused | no | - | not_applicable |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | PRIMARY | predictive | 4/7 | 0.571 | 2.529 | 0.286 | 2.00 / 0.120 | coverage-gap | no | - | continuous |
| resconi_2023_pbma_beef_identity_benchmark | trust_loop | SECONDARY | predictive | 0/1 | 0.000 | - | - | 1.50 / 0.100 | refused | no | - | not_applicable |
| soy_isolate_40C_PratapSingh2021 | trust_loop | PRIMARY | predictive | 0/2 | 0.000 | - | - | 2.00 / 0.120 | refused | no | - | not_applicable |
| thiamine_cys_glucose_120C_Bolton1994 | trust_loop | PRIMARY | predictive | 0/1 | 0.000 | - | - | 3.00 / 0.480 | refused | no | - | 2.31 |
| mp_holdout_arp_glutathione_120C_120min_pH75_Zhou2023c | maillard_path_holdout | REFERENCE | external_holdout | 0/3 | 0.000 | - | - | 1.50 / 0.100 | refused | no | - | ambiguous |
| mp_holdout_arp_glutathione_120C_30min_pH75_Zhou2023c | maillard_path_holdout | REFERENCE | external_holdout | 0/3 | 0.000 | - | - | 1.50 / 0.100 | refused | no | - | ambiguous |
| mp_holdout_arp_glutathione_120C_60min_pH75_Zhou2023c | maillard_path_holdout | REFERENCE | external_holdout | 0/3 | 0.000 | - | - | 1.50 / 0.100 | refused | no | - | ambiguous |
| mp_holdout_arp_glutathione_120C_90min_pH75_Zhou2023c | maillard_path_holdout | REFERENCE | external_holdout | 0/3 | 0.000 | - | - | 1.50 / 0.100 | refused | no | - | ambiguous |
| mp_holdout_fructose_asparagine_180C_Lin2022 | maillard_path_holdout | REFERENCE | external_holdout | 1/2 | 0.500 | 9.149 | 0.961 | 1.50 / 0.100 | coverage-gap | no | - | ambiguous |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | maillard_path_holdout | REFERENCE | external_holdout | 2/3 | 0.667 | 52.799 | 0.868 | 1.50 / 0.100 | coverage-gap | no | - | ambiguous |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | maillard_path_holdout | REFERENCE | external_holdout | 2/3 | 0.667 | 269.863 | 1.345 | 1.50 / 0.100 | coverage-gap | no | - | ambiguous |
| mp_holdout_glucose_asparagine_180C_10min_Chang2021 | maillard_path_holdout | REFERENCE | external_holdout | 1/1 | 1.000 | 218.216 | 2.339 | 1.50 / 0.100 | scale-gap | no | - | ambiguous |
| mp_holdout_glucose_asparagine_180C_30min_Chang2021 | maillard_path_holdout | REFERENCE | external_holdout | 1/1 | 1.000 | 2.522 | 0.402 | 1.50 / 0.100 | scale-gap | no | - | ambiguous |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | maillard_path_holdout | REFERENCE | external_holdout | 2/2 | 1.000 | 4.423 | 0.542 | 1.50 / 0.100 | scale-gap | no | - | ambiguous |
| mp_holdout_glucose_asparagine_180C_Ye2024 | maillard_path_holdout | REFERENCE | external_holdout | 1/1 | 1.000 | 49.824 | 1.697 | 1.50 / 0.100 | scale-gap | no | - | 0.184 mmol O2, no thiol |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | maillard_path_holdout | REFERENCE | external_holdout | 6/6 | 1.000 | 34.237 | 0.702 | 1.50 / 0.100 | ranking-gap | no | - | 0.272 mmol O2, no thiol |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | maillard_path_holdout | REFERENCE | external_holdout | 0/2 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 0.27 |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | maillard_path_holdout | REFERENCE | external_holdout | 0/2 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 0.27 |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | maillard_path_holdout | REFERENCE | external_holdout | 0/2 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 0.27 |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | maillard_path_holdout | REFERENCE | external_holdout | 0/2 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 0.27 |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | maillard_path_holdout | REFERENCE | external_holdout | 0/2 | 0.000 | - | - | 1.10 / 0.041 | refused | no | - | 0.27 |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | maillard_path_holdout | REFERENCE | external_holdout | 0/2 | 0.000 | - | - | 1.50 / 0.100 | refused | no | - | 1.98 |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | maillard_path_holdout | REFERENCE | external_holdout | 0/2 | 0.000 | - | - | 1.50 / 0.100 | refused | no | - | 1.98 |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | maillard_path_holdout | REFERENCE | external_holdout | 0/2 | 0.000 | - | - | 1.50 / 0.100 | refused | no | - | 1.98 |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | maillard_path_holdout | REFERENCE | external_holdout | 0/2 | 0.000 | - | - | 1.50 / 0.100 | refused | no | - | 1.98 |
| mp_holdout_xylose_cysteine_phosphate_120C_50min_pH65_Zhang2025b | maillard_path_holdout | REFERENCE | external_holdout | 0/4 | 0.000 | - | - | 1.50 / 0.100 | refused | no | - | 0.44 |
| external_validation_bi_2020_raw_pea_hexanal | external_matrix | PRIMARY | external_holdout | 0/2 | 0.000 | - | - | 2.00 / 0.120 | refused | no | - | not_applicable |
| external_validation_bi_2020_roasted_pea_hexanal | external_matrix | PRIMARY | external_holdout | 0/4 | 0.000 | - | - | 2.00 / 0.120 | refused | no | - | open |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | PRIMARY | external_holdout | 3/4 | 0.750 | 365.609 | 1.331 | 2.00 / 0.120 | coverage-gap | no | - | continuous |
| external_validation_liu_2023_ppi_offnote_baseline | external_matrix | PRIMARY | external_holdout | 0/2 | 0.000 | - | - | 2.00 / 0.120 | refused | no | - | not_applicable |

## Rows

| benchmark | compound | unit | measured | predicted | fold | within band | within contract | interval (ug/L) | measured inside | lane | in core fit |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pea_isolate_uht_140C_Trikusuma2019 | hexanal | ppb | 782.000 | 353.858 | 2.210 | yes | no | [1.087, 1.15e+05] | yes | lipid | no |
| pea_isolate_uht_140C_Trikusuma2019 | 2-pentylfuran | ppb | 163.000 | 64.442 | 2.529 | yes | no | [0.198, 2.1e+04] | yes | lipid | no |
| pea_isolate_uht_140C_Trikusuma2019 | nonanal | ppb | 24.000 | 15.606 | 1.538 | yes | yes | [0.048, 5.08e+03] | yes | lipid | no |
| pea_isolate_uht_140C_Trikusuma2019 | (E,E)-2,4-decadienal | ppb | 46.900 | 75.702 | 1.614 | yes | yes | [0.233, 2.46e+04] | yes | lipid | no |
| mp_holdout_fructose_asparagine_180C_Lin2022 | Acrylamide | ppb | 1.86e+03 | 203.194 | 9.149 | no | no | [29.052, 1.42e+03] | no | acrylamide | no |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | DMHF | ppb | 1.15e+03 | 21.841 | 52.799 | no | no | [2.223, 214.535] | no | acrylamide | no |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | 5-Hydroxymethylfurfural (HMF) | ppb | 5.73e+04 | 5.54e+04 | 1.033 | yes | yes | [7.92e+03, 3.88e+05] | yes | acrylamide | no |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | DMHF | ppb | 5.89e+03 | 21.841 | 269.863 | no | no | [2.223, 214.535] | no | acrylamide | no |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | 5-Hydroxymethylfurfural (HMF) | ppb | 1.01e+05 | 5.54e+04 | 1.816 | yes | no | [7.92e+03, 3.88e+05] | yes | acrylamide | no |
| mp_holdout_glucose_asparagine_180C_10min_Chang2021 | Acrylamide | ppb | 28.000 | 6.11e+03 | 218.216 | no | no | [873.594, 4.27e+04] | no | acrylamide | no |
| mp_holdout_glucose_asparagine_180C_30min_Chang2021 | Acrylamide | ppb | 1.46e+03 | 3.68e+03 | 2.522 | yes | no | [526.139, 2.57e+04] | yes | acrylamide | no |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | Acrylamide | ppb | 832.000 | 3.68e+03 | 4.423 | no | no | [526.139, 2.57e+04] | yes | acrylamide | no |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | 5-Hydroxymethylfurfural (HMF) | ppb | 7e+03 | 2.55e+03 | 2.744 | yes | no | [364.770, 1.78e+04] | yes | acrylamide | no |
| mp_holdout_glucose_asparagine_180C_Ye2024 | Acrylamide | umol_per_mol_limiting_precursor | 140.580 | 7e+03 | 49.824 | no | no | - | - | acrylamide | no |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 5-Hydroxymethylfurfural (HMF) | ppb | 1.74e+04 | 1.96e+03 | 8.891 | no | no | [279.807, 1.37e+04] | no | trunk | no |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 3-deoxyglucosone | ppb | 5.22e+04 | 6.55e+04 | 1.255 | yes | yes | [9.37e+03, 4.58e+05] | yes | trunk | no |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 3,4-dideoxyglucosone | ppb | 5.55e+04 | 8.3e+03 | 6.690 | no | no | [1.19e+03, 5.8e+04] | yes | trunk | no |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | glucosone | ppb | 7.5e+03 | 219.059 | 34.237 | no | no | [31.320, 1.53e+03] | no | trunk | no |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | glyoxal | ppb | 5.6e+03 | 1.16e+03 | 4.824 | no | no | [165.967, 8.12e+03] | yes | trunk | no |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | methylglyoxal | ppb | 2.6e+03 | 3.45e+03 | 1.327 | yes | yes | [493.442, 2.41e+04] | yes | trunk | no |
| external_validation_li_2026_spi_wg_hme_control | 2-pentylfuran | ppb | 5.63e+03 | 15.387 | 365.609 | no | no | [0.138, 1.71e+03] | no | lipid | no |
| external_validation_li_2026_spi_wg_hme_control | hexanal | ppb | 605.600 | 70.264 | 8.619 | no | no | [0.631, 7.82e+03] | yes | lipid | no |
| external_validation_li_2026_spi_wg_hme_control | nonanal | ppb | 74.370 | 23.885 | 3.114 | no | no | [0.215, 2.65e+03] | yes | lipid | no |

## Refused rows

| benchmark | panel | compound | reason |
|---|---|---|---|
| acrylamide_spi_extrusion_130C_ACSRef3 | trust_loop | acrylamide | REFUSED BY THE DATA: precursor charge not printed. Ma 2024 (secs. 3.1-3.2, Table 1) adds no asparagine or glucose and prints no free asparagine, reducing sugar, pH or water activity, so the pot the row was integrated at  |
| cml_cel_commercial_pbma_Foods2023 | trust_loop | Nε-(Carboxymethyl)lysine (CML) | GLYCATION TARGETS 'Nε-(Carboxymethyl)lysine (CML)' (wave B20) run on the trunk lane only: the acrylamide lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + ami |
| cml_cel_commercial_pbma_Foods2023 | trust_loop | Nε-(Carboxyethyl)lysine (CEL) | GLYCATION TARGETS 'Nε-(Carboxyethyl)lysine (CEL)' (wave B20) run on the trunk lane only: the acrylamide lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + amin |
| hofmann1998_c2c3_recombination_145C_20min_pH3 | trust_loop | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| hofmann1998_c2c3_recombination_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| hofmann1998_c2c3_recombination_145C_20min_pH7 | trust_loop | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_furan2aldehyde_h2s_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_norfuraneol_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| hofmann1998_norfuraneol_h2s_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | NOT EVALUABLE: REFUSED, FALSIFIED pH BALANCE (S3; M1b part (i), CF-084, CF-134): computed on the sulfur lane's pH solved from the shipped charge balance, which S3 falsified on CF-084's controls with no rival credited (M1 |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: REFUSED, FALSIFIED pH BALANCE (S3; M1b part (i), CF-084, CF-134): computed on the sulfur lane's pH solved from the shipped charge balance, which S3 falsified on CF-084's controls with no rival credited (M1 |
| pea_isolate_40C_PratapSingh2021 | trust_loop | hexanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's source says so, recorded as closure = 'never heated' (its timed block is the headspace in |
| pea_isolate_40C_PratapSingh2021 | trust_loop | 2-pentylfuran | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's source says so, recorded as closure = 'never heated' (its timed block is the headspace in |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | 2,5-dimethylpyrazine | PYRAZINE TARGETS '2,5-dimethylpyrazine' (wave B18) run on the trunk lane only: the lipid lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + amine pot that reso |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | methional | METHIONINE CHAIN TARGETS 'methional' (wave B22) run on the trunk lane only: the lipid lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + amine pot that resolve |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | 2-acetyl-1-pyrroline | 2-ACETYL-1-PYRROLINE TARGETS '2-acetyl-1-pyrroline' (wave B24) run on the trunk lane only: the lipid lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + amine p |
| resconi_2023_pbma_beef_identity_benchmark | trust_loop | furfural | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| soy_isolate_40C_PratapSingh2021 | trust_loop | hexanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's source says so, recorded as closure = 'never heated' (its timed block is the headspace in |
| soy_isolate_40C_PratapSingh2021 | trust_loop | 2-pentylfuran | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's source says so, recorded as closure = 'never heated' (its timed block is the headspace in |
| thiamine_cys_glucose_120C_Bolton1994 | trust_loop | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| mp_holdout_arp_glutathione_120C_120min_pH75_Zhou2023c | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| mp_holdout_arp_glutathione_120C_120min_pH75_Zhou2023c | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| mp_holdout_arp_glutathione_120C_120min_pH75_Zhou2023c | maillard_path_holdout | Furfural | NOT EVALUABLE: REFUSED, FALSIFIED pH BALANCE (S3; M1b part (i), CF-084, CF-134): computed on the sulfur lane's pH solved from the shipped charge balance, which S3 falsified on CF-084's controls with no rival credited (M1 |
| mp_holdout_arp_glutathione_120C_30min_pH75_Zhou2023c | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| mp_holdout_arp_glutathione_120C_30min_pH75_Zhou2023c | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| mp_holdout_arp_glutathione_120C_30min_pH75_Zhou2023c | maillard_path_holdout | Furfural | NOT EVALUABLE: REFUSED, FALSIFIED pH BALANCE (S3; M1b part (i), CF-084, CF-134): computed on the sulfur lane's pH solved from the shipped charge balance, which S3 falsified on CF-084's controls with no rival credited (M1 |
| mp_holdout_arp_glutathione_120C_60min_pH75_Zhou2023c | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| mp_holdout_arp_glutathione_120C_60min_pH75_Zhou2023c | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| mp_holdout_arp_glutathione_120C_60min_pH75_Zhou2023c | maillard_path_holdout | Furfural | NOT EVALUABLE: REFUSED, FALSIFIED pH BALANCE (S3; M1b part (i), CF-084, CF-134): computed on the sulfur lane's pH solved from the shipped charge balance, which S3 falsified on CF-084's controls with no rival credited (M1 |
| mp_holdout_arp_glutathione_120C_90min_pH75_Zhou2023c | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| mp_holdout_arp_glutathione_120C_90min_pH75_Zhou2023c | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| mp_holdout_arp_glutathione_120C_90min_pH75_Zhou2023c | maillard_path_holdout | Furfural | NOT EVALUABLE: REFUSED, FALSIFIED pH BALANCE (S3; M1b part (i), CF-084, CF-134): computed on the sulfur lane's pH solved from the shipped charge balance, which S3 falsified on CF-084's controls with no rival credited (M1 |
| mp_holdout_fructose_asparagine_180C_Lin2022 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | REFUSED BY THE DATA: within-source conflict (CF-127): the prose prints 12.28 ppm citing Fig. 4A, but Fig. 4A's asparagine + fructose 30-min bar (water) reads ~5.0 ppm (SD to ~6) on an axis checked against the printed 91. |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | maillard_path_holdout | Furfural | LANE CONFLICT: this request needs both the acrylamide and sulfur lanes at once. They do not compose -- the acrylamide network deliberately omits every sulfur step (acrylamide.OUT_OF_SCOPE), because composing them would s |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | maillard_path_holdout | Furfural | LANE CONFLICT: this request needs both the acrylamide and sulfur lanes at once. They do not compose -- the acrylamide network deliberately omits every sulfur step (acrylamide.OUT_OF_SCOPE), because composing them would s |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: REFUSED, FALSIFIED pH BALANCE (S3; M1b part (i), CF-084, CF-134): computed on the sulfur lane's pH solved from the shipped charge balance, which S3 falsified on CF-084's controls with no rival credited (M1 |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: REFUSED, FALSIFIED pH BALANCE (S3; M1b part (i), CF-084, CF-134): computed on the sulfur lane's pH solved from the shipped charge balance, which S3 falsified on CF-084's controls with no rival credited (M1 |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: REFUSED, FALSIFIED pH BALANCE (S3; M1b part (i), CF-084, CF-134): computed on the sulfur lane's pH solved from the shipped charge balance, which S3 falsified on CF-084's controls with no rival credited (M1 |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: REFUSED, FALSIFIED pH BALANCE (S3; M1b part (i), CF-084, CF-134): computed on the sulfur lane's pH solved from the shipped charge balance, which S3 falsified on CF-084's controls with no rival credited (M1 |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| mp_holdout_xylose_cysteine_phosphate_120C_50min_pH65_Zhang2025b | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: REFUSED, FALSIFIED pH BALANCE (S3; M1b part (i), CF-084, CF-134): computed on the sulfur lane's pH solved from the shipped charge balance, which S3 falsified on CF-084's controls with no rival credited (M1 |
| mp_holdout_xylose_cysteine_phosphate_120C_50min_pH65_Zhang2025b | maillard_path_holdout | 2-Furfurylthiol (FFT) | NOT EVALUABLE: ABSOLUTE LEVEL REFUSED (H0 / EX-059): the thiol lane ships relative predictions only (D2 branch (i): tau_ABS 0.964 dex, U >= 1.60; absolute thiol evidence from one laboratory). Lifted by absolute measureme |
| mp_holdout_xylose_cysteine_phosphate_120C_50min_pH65_Zhang2025b | maillard_path_holdout | Bis(2-methyl-3-furyl) disulfide | NOT EVALUABLE: REFUSED, FALSIFIED pH BALANCE (S3; M1b part (i), CF-084, CF-134): computed on the sulfur lane's pH solved from the shipped charge balance, which S3 falsified on CF-084's controls with no rival credited (M1 |
| mp_holdout_xylose_cysteine_phosphate_120C_50min_pH65_Zhang2025b | maillard_path_holdout | 2-Acetylthiazole | NOT EVALUABLE: REFUSED, FALSIFIED pH BALANCE (S3; M1b part (i), CF-084, CF-134): computed on the sulfur lane's pH solved from the shipped charge balance, which S3 falsified on CF-084's controls with no rival credited (M1 |
| external_validation_bi_2020_raw_pea_hexanal | external_matrix | hexanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's source says so, recorded as closure = 'never heated' (its timed block is the headspace in |
| external_validation_bi_2020_raw_pea_hexanal | external_matrix | nonanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's source says so, recorded as closure = 'never heated' (its timed block is the headspace in |
| external_validation_bi_2020_roasted_pea_hexanal | external_matrix | hexanal | REFUSED BY THE DATA: net loss of a carried level (CF-129): raw peas carry 1260 +/- 114 ug/kg hexanal and the roasted peas measure 324 +/- 34.3 (Bi 2020 Table 4), a 3.9x loss on roasting. The lipid lane adds a carried lev |
| external_validation_bi_2020_roasted_pea_hexanal | external_matrix | furaneol | THIS POT CHARGES NO PRECURSOR THAT COULD MAKE 'furaneol'. The charge declares only a matrix/lipid carrier, which is not a precursor: it resolves to a hydroperoxide pool for the lipid lane and charges NOTHING on the trunk |
| external_validation_bi_2020_roasted_pea_hexanal | external_matrix | 2,5-dimethylpyrazine | PYRAZINE TARGETS '2,5-dimethylpyrazine' (wave B18) run on the trunk lane only: the lipid lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + amine pot that reso |
| external_validation_bi_2020_roasted_pea_hexanal | external_matrix | furfural | THIS POT CHARGES NO PRECURSOR THAT COULD MAKE 'furfural'. The charge declares only a matrix/lipid carrier, which is not a precursor: it resolves to a hydroperoxide pool for the lipid lane and charges NOTHING on the trunk |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | 1-hexanol | UNREPRESENTED TARGETS: 1-hexanol -- The lipid lane exists and forms the SIX products Frankel 1989 measured, but 1-hexanol is not one of them and NO aldehyde-reduction step is measured anywhere in the corpus -- in a therm |
| external_validation_liu_2023_ppi_offnote_baseline | external_matrix | hexanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's source says so, recorded as closure = 'never heated' (its timed block is the headspace in |
| external_validation_liu_2023_ppi_offnote_baseline | external_matrix | nonanal | THIS POT WAS NEVER COOKED, so what it measures is what the raw material ARRIVED WITH, and this lane models FORMATION. The bundle's source says so, recorded as closure = 'never heated' (its timed block is the headspace in |

## Bundles kept off the scored panel

- pea_isolate_ribose_cysteine_100C_45min_Internal2026 (trust_loop): synthetic snapshot (_Internal2026): legacy-model output, not a measurement
- soy_isolate_ribose_cysteine_100C_45min_Internal2026 (trust_loop): synthetic snapshot (_Internal2026): legacy-model output, not a measurement
