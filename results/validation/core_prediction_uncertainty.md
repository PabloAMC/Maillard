# Core prediction uncertainty (Monte-Carlo envelope on the kinetic core)

n_samples = 200, seed = 0, CI level = 90 %.

* benchmarks with an envelope: **23** of 37 on the panel; matched rows **45**; refused rows 32
* mixed-population coverage: 12/45 (0.267)
* **honest literature coverage: 12/42 (0.286)**, median CI width 0.919 log10; 3 not evaluable; 0 fitted rows excluded
* out-of-sample literature coverage: 12/41 (3 not evaluable); rows the core fit read: {'hits': 0, 'total': 1, 'not_evaluable': 0}
* sampled priors 55, fixed 46; lanes with NO sampled fit uncertainty: none
* observable bands (K_aw, HS-SPME) applied by quantification family -- rows: headspace 9, extraction 36, undeclared 0

## Per panel

| panel | hits | total | rate | median width (log10) | not evaluable |
|---|---|---|---|---|---|
| external_matrix | 3 | 4 | 0.750 | 2.132 | 0 |
| maillard_path_holdout | 4 | 28 | 0.143 | 0.914 | 3 |
| trust_loop | 5 | 10 | 0.500 | 1.438 | 0 |

## Rows

| benchmark | panel | compound | unit | measured | point | p5 | p50 | p95 | inside | width | obs bands | lane | role |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| acrylamide_spi_extrusion_130C_ACSRef3 | trust_loop | acrylamide | ppb | 150.000 | 0.016 | 0.008 | 0.016 | 0.026 | no | 0.512 | no (extraction) | acrylamide | predictive |
| hofmann1998_norfuraneol_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 1.02e+03 | 2.35e+03 | 1.62e+03 | 2.28e+03 | 2.86e+03 | no | 0.247 | no (extraction) | sulfur | predictive [in core fit] |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | ppb | 121.000 | 830.643 | 37.394 | 652.651 | 8.57e+03 | yes | 2.360 | no (extraction) | sulfur | predictive |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 198.000 | 730.638 | 411.166 | 672.826 | 1.23e+03 | no | 0.475 | no (extraction) | sulfur | predictive |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | hexanal | ppb | 782.000 | 353.839 | 73.661 | 402.251 | 2.03e+03 | yes | 1.441 | yes (headspace) | lipid | predictive |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | 2-pentylfuran | ppb | 163.000 | 64.442 | 13.522 | 74.365 | 368.263 | yes | 1.435 | yes (headspace) | lipid | predictive |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | nonanal | ppb | 24.000 | 15.606 | 2.362 | 17.719 | 141.735 | yes | 1.778 | yes (headspace) | lipid | predictive |
| pea_isolate_uht_140C_Trikusuma2019 | trust_loop | (E,E)-2,4-decadienal | ppb | 46.900 | 75.955 | 3.279 | 70.934 | 1.13e+03 | yes | 2.538 | yes (headspace) | lipid | predictive |
| resconi_2023_pbma_beef_identity_benchmark | trust_loop | furfural | ppb | 715.220 | 7.665 | 1.405 | 11.014 | 90.298 | no | 1.808 | yes (headspace) | sulfur | predictive |
| thiamine_cys_glucose_120C_Bolton1994 | trust_loop | 2-Methyl-3-furanthiol (MFT) | ppb | 11.700 | 235.970 | 48.719 | 104.031 | 216.682 | no | 0.648 | no (extraction) | sulfur | predictive |
| mp_holdout_fructose_asparagine_180C_Lin2022 | maillard_path_holdout | Acrylamide | ppb | 1.86e+03 | 203.194 | 69.866 | 186.040 | 573.774 | no | 0.914 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_fructose_asparagine_180C_Lin2022 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | ppb | 1.23e+04 | 1.95e+03 | 859.244 | 1.93e+03 | 1.95e+03 | no | 0.355 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | maillard_path_holdout | DMHF | ppb | 1.15e+03 | 21.841 | 7.414 | 21.229 | 62.215 | no | 0.924 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | ppb | 5.73e+04 | 2.8e+04 | 3.99e+03 | 2.7e+04 | 2.8e+04 | no | 0.846 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | maillard_path_holdout | DMHF | ppb | 5.89e+03 | 21.841 | 7.414 | 21.229 | 62.215 | no | 0.924 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | ppb | 1.01e+05 | 2.8e+04 | 3.99e+03 | 2.7e+04 | 2.8e+04 | no | 0.846 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_asparagine_180C_10min_Chang2021 | maillard_path_holdout | Acrylamide | ppb | 28.000 | 6.11e+03 | 1.88e+03 | 5.38e+03 | 1.54e+04 | no | 0.913 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_asparagine_180C_30min_Chang2021 | maillard_path_holdout | Acrylamide | ppb | 1.46e+03 | 3.68e+03 | 1.27e+03 | 3.36e+03 | 1.04e+04 | yes | 0.914 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | maillard_path_holdout | Acrylamide | ppb | 832.000 | 3.68e+03 | 1.27e+03 | 3.36e+03 | 1.04e+04 | no | 0.914 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | ppb | 7e+03 | 2.19e+03 | 917.342 | 2.15e+03 | 2.2e+03 | no | 0.380 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_asparagine_180C_Ye2024 | maillard_path_holdout | Acrylamide | umol_per_mol_limiting_precursor | 140.580 | 7e+03 | 2.51e+03 | 6.39e+03 | 2.06e+04 | no | 0.914 | no (extraction) | acrylamide | external_holdout |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | maillard_path_holdout | 5-Hydroxymethylfurfural (HMF) | ppb | 1.74e+04 | 1.46e+03 | 1.07e+03 | 1.45e+03 | 1.46e+03 | no | 0.135 | no (extraction) | trunk | external_holdout |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | maillard_path_holdout | 3-deoxyglucosone | ppb | 5.22e+04 | 4.7e+04 | 4.7e+04 | 4.7e+04 | 4.7e+04 | no | 1.47e-06 | no (extraction) | trunk | external_holdout |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | maillard_path_holdout | 3,4-dideoxyglucosone | ppb | 5.55e+04 | 1.71e+03 | 1.71e+03 | 1.71e+03 | 1.71e+03 | no | 1.05e-06 | no (extraction) | trunk | external_holdout |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | maillard_path_holdout | glucosone | ppb | 7.5e+03 | 119.059 | 119.059 | 119.059 | 119.059 | no | 2.05e-06 | no (extraction) | trunk | external_holdout |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | maillard_path_holdout | glyoxal | ppb | 5.6e+03 | 160.794 | 167.306 | 197.878 | 205.081 | no | 0.088 | no (extraction) | trunk | external_holdout |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | maillard_path_holdout | methylglyoxal | ppb | 2.6e+03 | 2.03e+03 | 1.55e+03 | 1.97e+03 | 2.16e+03 | no | 0.145 | no (extraction) | trunk | external_holdout |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 229.000 | 1.13e+03 | 67.468 | 875.815 | 1.13e+04 | yes | 2.223 | no (extraction) | sulfur | external_holdout (shared: hofmann_ribose_pH3_FFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 553.000 | 23.470 | 12.689 | 21.612 | 41.846 | no | 0.518 | no (extraction) | sulfur | external_holdout (shared: hofmann_ribose_pH3_MFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 12.000 | 0.090 | 0.00441 | 0.075 | 0.495 | no | 2.050 | no (extraction) | sulfur | external_holdout (shared: hofmann_ribose_pH7_FFT) |
| mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 25.000 | 0.178 | 0.042 | 0.153 | 0.262 | no | 0.801 | no (extraction) | sulfur | external_holdout (shared: hofmann_ribose_pH7_MFT) |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 96.000 | 191.058 | 23.187 | 191.119 | 1.63e+03 | yes | 1.847 | no (extraction) | sulfur | external_holdout |
| mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 143.000 | 546.398 | 334.422 | 499.113 | 785.835 | no | 0.371 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 6.880 | 64.086 | 1.661 | 9.566 | 53.145 | yes | 1.505 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 1.280 | 610.255 | 20.981 | 184.346 | 1.12e+03 | no | 1.729 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 3.290 | 96.892 | 5.920 | 23.853 | 85.259 | no | 1.158 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 1.460 | 535.496 | 29.321 | 251.659 | 1.28e+03 | no | 1.640 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 2.400 | 134.513 | 17.946 | 53.943 | 144.072 | no | 0.905 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 1.680 | 363.445 | 34.506 | 257.905 | 1.12e+03 | no | 1.511 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | maillard_path_holdout | 2-Methyl-3-furanthiol (MFT) | ppb | 1.710 | 170.534 | 45.250 | 102.179 | 228.443 | no | 0.703 | no (extraction) | sulfur | external_holdout |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | maillard_path_holdout | 2-Furfurylthiol (FFT) | ppb | 1.620 | 211.740 | 34.512 | 194.988 | 944.426 | no | 1.437 | no (extraction) | sulfur | external_holdout |
| external_validation_bi_2020_roasted_pea_hexanal | external_matrix | hexanal | ppb | 324.000 | 88.598 | 7.395 | 69.413 | 1.06e+03 | yes | 2.157 | yes (headspace) | lipid | external_holdout |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | 2-pentylfuran | ppb | 5.63e+03 | 15.387 | 1.170 | 10.861 | 158.531 | no | 2.132 | yes (headspace) | lipid | external_holdout |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | hexanal | ppb | 605.600 | 69.695 | 5.301 | 49.195 | 718.039 | yes | 2.132 | yes (headspace) | lipid | external_holdout |
| external_validation_li_2026_spi_wg_hme_control | external_matrix | nonanal | ppb | 72.660 | 23.885 | 1.817 | 16.860 | 246.079 | yes | 2.132 | yes (headspace) | lipid | external_holdout |

## Refused rows

| benchmark | panel | compound | reason |
|---|---|---|---|
| cml_cel_commercial_pbma_Foods2023 | trust_loop | Nε-(Carboxymethyl)lysine (CML) | GLYCATION TARGETS 'Nε-(Carboxymethyl)lysine (CML)' (wave B20) run on the trunk lane only: the acrylamide lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + ami |
| cml_cel_commercial_pbma_Foods2023 | trust_loop | Nε-(Carboxyethyl)lysine (CEL) | GLYCATION TARGETS 'Nε-(Carboxyethyl)lysine (CEL)' (wave B20) run on the trunk lane only: the acrylamide lane's network keeps the topology its fit was run on and carries these species inert. Ask for them in a sugar + amin |
| hofmann1998_c2c3_recombination_145C_20min_pH3 | trust_loop | 2-Methyl-3-furanthiol (MFT) | UNMAPPED PRECURSORS 'Hydroxyacetaldehyde', 'Mercapto-2-propanone': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge. /  |
| hofmann1998_c2c3_recombination_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | UNMAPPED PRECURSORS 'Hydroxyacetaldehyde', 'Mercapto-2-propanone': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge. /  |
| hofmann1998_c2c3_recombination_145C_20min_pH7 | trust_loop | 2-Methyl-3-furanthiol (MFT) | UNMAPPED PRECURSORS 'Hydroxyacetaldehyde', 'Mercapto-2-propanone': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge. /  |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_furan2aldehyde_h2s_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | UNMAPPED PRECURSORS 'Furan-2-aldehyde', 'Hydrogen sulfide': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge. / THIS PO |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | 2-Furfurylthiol (FFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (FFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | NOT EVALUABLE: HEXOSE ENTRY UNIDENTIFIED (MFT): the only route from a hexose to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary measurement identifies (the primary-evidence refit left them  |
| hofmann1998_norfuraneol_h2s_145C_20min_pH5 | trust_loop | 2-Methyl-3-furanthiol (MFT) | UNMAPPED PRECURSORS 'Hydrogen sulfide': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge. |
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
| b13.k_da_sink.log10_k_100C | trunk | fitted_rate | uniform_band | - | - | [-12.000, -1.940] | yes | two_laboratories_disagree: flat across the disagreement |
| b13.k_go_sink.ea_kj_mol | trunk | fitted_ea | uniform_band | 0.000 | - | [0.000, 150.800] | yes | two_laboratories_disagree: flat across the disagreement |
| b13.k_odg_da.log10_k_100C | trunk | fitted_rate | uniform_band | -5.640 | - | [-4.630, -1.960] | yes | two_laboratories_disagree: flat across the disagreement |
| b13.k_hmf_self.log10_k_100C | trunk | fitted_rate | uniform_band | -6.047 | - | [-6.050, -0.960] | yes | two_laboratories_disagree: flat across the disagreement |
| b13.k_glc_g.log10_k_100C | trunk | fitted_rate | fixed | -7.272 | - | - | no | NOT SAMPLED, and not by oversight. k_tdg_ddg (1.5x), k_ddg_hmf (1.13x) and k_go_sink's RATE (1.87x) are the first cross- |
| b13.k_g_go.log10_k_100C | trunk | fitted_rate | fixed | -2.451 | - | - | no | NOT SAMPLED, and not by oversight. k_tdg_ddg (1.5x), k_ddg_hmf (1.13x) and k_go_sink's RATE (1.87x) are the first cross- |
| b13.k_tdg_ddg.log10_k_100C | trunk | fitted_rate | fixed | -2.428 | - | - | no | NOT SAMPLED, and not by oversight. k_tdg_ddg (1.5x), k_ddg_hmf (1.13x) and k_go_sink's RATE (1.87x) are the first cross- |
| b13.k_ddg_hmf.log10_k_100C | trunk | fitted_rate | fixed | -0.924 | - | - | no | NOT SAMPLED, and not by oversight. k_tdg_ddg (1.5x), k_ddg_hmf (1.13x) and k_go_sink's RATE (1.87x) are the first cross- |
| b18.log10_k_go_ak_100C | trunk | fitted_rate | normal_log10 | -6.542 | 0.082 | [-9.191, -5.191] | yes | laplace_covariance_at_b8_optimum |
| b18.ea_go_ak_kj_mol | trunk | fitted_ea | uniform_band | 103.100 | 16.158 | [100.590, 103.100] | yes | unidentified_direction_flat_across_its_declared_band |
| b18.log10_k_mgo_ak_100C | trunk | fitted_rate | normal_log10 | -7.530 | 0.075 | [-10.093, -6.093] | yes | laplace_covariance_at_b8_optimum |
| b18.ea_mgo_ak_kj_mol | trunk | fitted_ea | uniform_band | 114.900 | 16.097 | [111.660, 114.900] | yes | unidentified_direction_flat_across_its_declared_band |
| b18.ph_slope_above_7_decades_per_unit | trunk | fitted_rate | fixed | 0.197 | 0.039 | [0.000, 1.500] | no | IDENTIFIED AND STILL NOT SAMPLED, and the interval is narrower for it. This coordinate is a module-level constant in par |
| b18.ph_slope_below_7_decades_per_unit | trunk | fitted_rate | fixed | 0.580 | 0.046 | [0.000, 1.500] | no | IDENTIFIED AND STILL NOT SAMPLED, and the interval is narrower for it. This coordinate is a module-level constant in par |
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
