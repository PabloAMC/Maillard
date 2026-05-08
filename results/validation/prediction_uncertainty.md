# Prediction Uncertainty Envelope

_Monte Carlo propagation of barrier-family offset priors (additive Gaussian, kcal/mol) through the benchmark evaluator. CI = 90% (P5–P95)._

**Headline trust metric**: measured value lies inside 90% CI for **39 / 48** matched compounds (**81.2%**).

Samples per benchmark: 200; seed 0; benchmarks evaluated: 16.

## Priors

| Family offset key | σ (kcal/mol) | Source |
| --- | --- | --- |
| `schiff_condensation` | 1.50 | barrier_constants_default |
| `amadori_rearrangement` | 1.50 | barrier_constants_default |
| `1,2-enolisation` | 1.50 | barrier_constants_default |
| `2,3-enolisation` | 1.50 | barrier_constants_default |
| `strecker_degradation` | 1.50 | barrier_constants_default |
| `cysteine_thermolysis` | 1.50 | barrier_constants_default |
| `thiol_addition` | 1.50 | barrier_constants_default |
| `thiol_addition_hexose` | 1.50 | barrier_constants_default |
| `thiol_oxidation` | 1.50 | barrier_constants_default |
| `aminoketone_condensation` | 1.50 | barrier_constants_default |
| `retro_aldol` | 1.50 | barrier_constants_default |
| `dehydration` | 1.50 | barrier_constants_default |
| `beta_elimination` | 1.50 | barrier_constants_default |
| `lipid_thiazole` | 1.50 | barrier_constants_default |

## Per-benchmark envelopes

### `acrylamide_asparagine_glucose_Parker2012`

- Execution path: `free_precursor`
- Matched compounds with envelope: 1

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| acrylamide | 1.5e+03 | 1.52e+03 | 1.52e+03 | 1.52e+03 | 0.00 | ✗ |

### `acrylamide_spi_extrusion_130C_ACSRef3`

- Execution path: `free_precursor`
- Matched compounds with envelope: 1

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| acrylamide | 62.6 | 60.4 | 60.4 | 60.4 | 0.00 | ✗ |

### `cml_cel_commercial_pbma_Foods2023`

- Execution path: `free_precursor`
- Matched compounds with envelope: 2

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| Nε-(Carboxymethyl)lysine (CML) | 32 | 26.6 | 26.6 | 26.6 | 0.00 | ✗ |
| Nε-(Carboxyethyl)lysine (CEL) | 55 | 51.7 | 51.7 | 51.7 | 0.00 | ✗ |

### `cys_glucose_150C_Farmer1999`

- Execution path: `free_precursor`
- Matched compounds with envelope: 3

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| 2-methyl-3-furanthiol | 15 | 0.0059 | 35 | 785 | 5.12 | ✓ |
| furfural | 450 | 171 | 408 | 478 | 0.45 | ✓ |
| pyrazine | 30 | 0.0274 | 23.7 | 196 | 3.85 | ✓ |

### `cys_ribose_140C_Hofmann1998`

- Execution path: `free_precursor`
- Matched compounds with envelope: 2

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| 2-methyl-3-furanthiol | 342 | 0.411 | 282 | 742 | 3.26 | ✓ |
| 2-furfurylthiol | 200 | 54.9 | 173 | 271 | 0.69 | ✓ |

### `cys_ribose_150C_Mottram1994`

- Execution path: `free_precursor`
- Matched compounds with envelope: 3

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| 2-methyl-3-furanthiol | 120 | 0.223 | 243 | 1.12e+03 | 3.70 | ✓ |
| bis(2-methyl-3-furyl) disulfide | 45 | 0.0172 | 17.6 | 494 | 4.46 | ✓ |
| furfural | 890 | 94 | 526 | 802 | 0.93 | ✗ |

### `furosine_extrusion_crossover_140C_RamirezJimenez2000`

- Execution path: `free_precursor`
- Matched compounds with envelope: 1

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| furosine | 87 | 86.5 | 86.5 | 86.5 | 0.00 | ✗ |

### `pea_isolate_ribose_cysteine_100C_45min_Internal2026`

- Execution path: `matrix_precursor_augmented`
- Matched compounds with envelope: 7

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| furfural | 0.0317 | 0.000742 | 0.0272 | 0.0353 | 1.68 | ✓ |
| 2-furfurylthiol | 0.00502 | 0.000117 | 0.0044 | 0.00649 | 1.74 | ✓ |
| 2-methyl-3-furanthiol | 0.00342 | 2.04e-06 | 0.00616 | 0.0393 | 4.28 | ✓ |
| bis(2-methyl-3-furyl) disulfide | 0.000324 | 3.24e-08 | 7.02e-05 | 0.00842 | 5.41 | ✓ |
| 2,5-dimethylpyrazine | 0.000234 | 1.5e-08 | 0.000168 | 0.00421 | 5.45 | ✓ |
| Hexanal | 0.000196 | 0.000196 | 0.000196 | 0.000196 | 0.00 | ✓ |
| Nonanal | 0.000102 | 0.000102 | 0.000102 | 0.000102 | 0.00 | ✓ |

### `pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026`

- Execution path: `matrix_precursor_augmented`
- Matched compounds with envelope: 6

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| 2-furfurylthiol | 0.0051 | 0.000117 | 0.0044 | 0.00649 | 1.74 | ✓ |
| 2-methyl-3-furanthiol | 0.0036 | 2.04e-06 | 0.00616 | 0.0393 | 4.28 | ✓ |
| bis(2-methyl-3-furyl) disulfide | 0.00031 | 3.24e-08 | 7.02e-05 | 0.00842 | 5.41 | ✓ |
| 2,5-dimethylpyrazine | 0.00024 | 1.5e-08 | 0.000168 | 0.00421 | 5.45 | ✓ |
| Hexanal | 0.00021 | 0.000196 | 0.000196 | 0.000196 | 0.00 | ✗ |
| Nonanal | 0.0001 | 0.000102 | 0.000102 | 0.000102 | 0.00 | ✗ |

### `resconi_2023_pbma_beef_identity_benchmark`

- Execution path: `free_precursor`
- Matched compounds with envelope: 1

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| furfural | 1.04e+03 | 342 | 817 | 956 | 0.45 | ✗ |

### `soy_isolate_ribose_cysteine_100C_45min_Internal2026`

- Execution path: `matrix_precursor_augmented`
- Matched compounds with envelope: 7

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| furfural | 0.276 | 0.00647 | 0.238 | 0.308 | 1.68 | ✓ |
| 2-furfurylthiol | 0.00707 | 0.000165 | 0.00619 | 0.00913 | 1.74 | ✓ |
| 2-methyl-3-furanthiol | 0.00482 | 2.87e-06 | 0.00868 | 0.0554 | 4.28 | ✓ |
| bis(2-methyl-3-furyl) disulfide | 0.000456 | 4.57e-08 | 9.88e-05 | 0.0119 | 5.41 | ✓ |
| 2,5-dimethylpyrazine | 0.000332 | 2.13e-08 | 0.000238 | 0.00598 | 5.45 | ✓ |
| Hexanal | 0.000203 | 0.000203 | 0.000203 | 0.000203 | 0.00 | ✓ |
| Nonanal | 0.000105 | 0.000105 | 0.000105 | 0.000105 | 0.00 | ✓ |

### `soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026`

- Execution path: `matrix_precursor_augmented`
- Matched compounds with envelope: 6

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| 2-furfurylthiol | 0.00707 | 0.000165 | 0.00619 | 0.00913 | 1.74 | ✓ |
| 2-methyl-3-furanthiol | 0.00482 | 2.87e-06 | 0.00868 | 0.0554 | 4.28 | ✓ |
| bis(2-methyl-3-furyl) disulfide | 0.000456 | 4.57e-08 | 9.88e-05 | 0.0119 | 5.41 | ✓ |
| 2,5-dimethylpyrazine | 0.000332 | 2.13e-08 | 0.000238 | 0.00598 | 5.45 | ✓ |
| Hexanal | 0.000203 | 0.000203 | 0.000203 | 0.000203 | 0.00 | ✓ |
| Nonanal | 0.000105 | 0.000105 | 0.000105 | 0.000105 | 0.00 | ✓ |

### `spi_hvp_xylose_120C_PMC9905368`

- Execution path: `free_precursor`
- Matched compounds with envelope: 3

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| 2-Furfurylthiol (FFT) | 0.42 | 0.0023 | 0.703 | 4.61 | 3.30 | ✓ |
| 2-Methyl-3-furanthiol (MFT) | 0.18 | 2.44e-05 | 0.424 | 50.6 | 6.32 | ✓ |
| Methional | 1.88 | 0.111 | 1.37 | 2.09 | 1.28 | ✓ |

### `thiamine_cys_ribose_100C_Hofmann1996`

- Execution path: `free_precursor`
- Matched compounds with envelope: 1

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| 2-Methyl-3-furanthiol (MFT) | 3.14 | 0.00119 | 3.92 | 21.1 | 4.25 | ✓ |

### `thiamine_cys_xylose_145C_Cerny2008`

- Execution path: `free_precursor`
- Matched compounds with envelope: 1

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| 2-Methyl-3-furanthiol (MFT) | 2.47 | 0.00299 | 3.43 | 22.3 | 3.87 | ✓ |

### `wheat_gluten_hvp_xylose_120C_PMC9905368`

- Execution path: `free_precursor`
- Matched compounds with envelope: 3

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| 2-Furfurylthiol (FFT) | 0.61 | 0.00375 | 0.949 | 7.52 | 3.30 | ✓ |
| 2-Methyl-3-furanthiol (MFT) | 0.34 | 4e-05 | 0.69 | 50 | 6.10 | ✓ |
| Methional | 3.44 | 0.138 | 2.48 | 3.91 | 1.45 | ✓ |
