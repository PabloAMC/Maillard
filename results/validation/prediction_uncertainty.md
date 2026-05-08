# Prediction Uncertainty Envelope

_Monte Carlo propagation of barrier-family offset priors (additive Gaussian, kcal/mol) through the benchmark evaluator. CI = 90% (P5–P95)._

**Headline trust metric**: measured value lies inside 90% CI for **39 / 48** matched compounds (**81.2%**).

Samples per benchmark: 200; seed 0; benchmarks evaluated: 16.

## Priors

| Key | Kind | σ (kcal/mol or log) | Source |
| --- | --- | --- | --- |
| `schiff_condensation` | barrier | 1.50 | barrier_constants_default |
| `amadori_rearrangement` | barrier | 1.50 | barrier_constants_default |
| `1,2-enolisation` | barrier | 1.50 | barrier_constants_default |
| `2,3-enolisation` | barrier | 1.50 | barrier_constants_default |
| `strecker_degradation` | barrier | 1.50 | barrier_constants_default |
| `cysteine_thermolysis` | barrier | 1.50 | barrier_constants_default |
| `thiol_addition` | barrier | 1.50 | barrier_constants_default |
| `thiol_addition_hexose` | barrier | 1.50 | barrier_constants_default |
| `thiol_oxidation` | barrier | 1.50 | barrier_constants_default |
| `aminoketone_condensation` | barrier | 1.50 | barrier_constants_default |
| `retro_aldol` | barrier | 1.50 | barrier_constants_default |
| `dehydration` | barrier | 1.50 | barrier_constants_default |
| `beta_elimination` | barrier | 1.50 | barrier_constants_default |
| `lipid_thiazole` | barrier | 1.50 | barrier_constants_default |
| `matrix_headspace` | observable | 0.30 | observable_multiplier_default |
| `henry_kaw` | observable | 0.30 | observable_multiplier_default |
| `matrix_retention` | observable | 0.20 | observable_multiplier_default |

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
| 2-methyl-3-furanthiol | 15 | 0.0043 | 31.7 | 587 | 5.14 | ✓ |
| furfural | 450 | 227 | 412 | 480 | 0.32 | ✓ |
| pyrazine | 30 | 0.0278 | 19.3 | 225 | 3.91 | ✓ |

### `cys_ribose_140C_Hofmann1998`

- Execution path: `free_precursor`
- Matched compounds with envelope: 2

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| 2-methyl-3-furanthiol | 342 | 1.16 | 257 | 742 | 2.81 | ✓ |
| 2-furfurylthiol | 200 | 58.9 | 174 | 273 | 0.67 | ✓ |

### `cys_ribose_150C_Mottram1994`

- Execution path: `free_precursor`
- Matched compounds with envelope: 3

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| 2-methyl-3-furanthiol | 120 | 0.457 | 206 | 1.09e+03 | 3.38 | ✓ |
| bis(2-methyl-3-furyl) disulfide | 45 | 0.0163 | 15.9 | 555 | 4.53 | ✓ |
| furfural | 890 | 124 | 539 | 806 | 0.81 | ✗ |

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
| furfural | 0.0317 | 0.00123 | 0.0277 | 0.0353 | 1.46 | ✓ |
| 2-furfurylthiol | 0.00502 | 0.000172 | 0.00443 | 0.00608 | 1.55 | ✓ |
| 2-methyl-3-furanthiol | 0.00342 | 2.13e-06 | 0.0058 | 0.0378 | 4.25 | ✓ |
| bis(2-methyl-3-furyl) disulfide | 0.000324 | 2.82e-08 | 7.91e-05 | 0.0114 | 5.61 | ✓ |
| 2,5-dimethylpyrazine | 0.000234 | 7.15e-08 | 0.000114 | 0.0058 | 4.91 | ✓ |
| Hexanal | 0.000196 | 0.000196 | 0.000196 | 0.000196 | 0.00 | ✓ |
| Nonanal | 0.000102 | 0.000102 | 0.000102 | 0.000102 | 0.00 | ✓ |

### `pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026`

- Execution path: `matrix_precursor_augmented`
- Matched compounds with envelope: 6

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| 2-furfurylthiol | 0.0051 | 0.000172 | 0.00443 | 0.00608 | 1.55 | ✓ |
| 2-methyl-3-furanthiol | 0.0036 | 2.13e-06 | 0.0058 | 0.0378 | 4.25 | ✓ |
| bis(2-methyl-3-furyl) disulfide | 0.00031 | 2.82e-08 | 7.91e-05 | 0.0114 | 5.61 | ✓ |
| 2,5-dimethylpyrazine | 0.00024 | 7.15e-08 | 0.000114 | 0.0058 | 4.91 | ✓ |
| Hexanal | 0.00021 | 0.000196 | 0.000196 | 0.000196 | 0.00 | ✗ |
| Nonanal | 0.0001 | 0.000102 | 0.000102 | 0.000102 | 0.00 | ✗ |

### `resconi_2023_pbma_beef_identity_benchmark`

- Execution path: `free_precursor`
- Matched compounds with envelope: 1

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| furfural | 1.04e+03 | 454 | 824 | 959 | 0.32 | ✗ |

### `soy_isolate_ribose_cysteine_100C_45min_Internal2026`

- Execution path: `matrix_precursor_augmented`
- Matched compounds with envelope: 7

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| furfural | 0.276 | 0.0108 | 0.242 | 0.308 | 1.46 | ✓ |
| 2-furfurylthiol | 0.00707 | 0.000242 | 0.00624 | 0.00855 | 1.55 | ✓ |
| 2-methyl-3-furanthiol | 0.00482 | 3e-06 | 0.00816 | 0.0532 | 4.25 | ✓ |
| bis(2-methyl-3-furyl) disulfide | 0.000456 | 3.96e-08 | 0.000111 | 0.0161 | 5.61 | ✓ |
| 2,5-dimethylpyrazine | 0.000332 | 1.02e-07 | 0.000163 | 0.00824 | 4.91 | ✓ |
| Hexanal | 0.000203 | 0.000203 | 0.000203 | 0.000203 | 0.00 | ✓ |
| Nonanal | 0.000105 | 0.000105 | 0.000105 | 0.000105 | 0.00 | ✓ |

### `soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026`

- Execution path: `matrix_precursor_augmented`
- Matched compounds with envelope: 6

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| 2-furfurylthiol | 0.00707 | 0.000242 | 0.00624 | 0.00855 | 1.55 | ✓ |
| 2-methyl-3-furanthiol | 0.00482 | 3e-06 | 0.00816 | 0.0532 | 4.25 | ✓ |
| bis(2-methyl-3-furyl) disulfide | 0.000456 | 3.96e-08 | 0.000111 | 0.0161 | 5.61 | ✓ |
| 2,5-dimethylpyrazine | 0.000332 | 1.02e-07 | 0.000163 | 0.00824 | 4.91 | ✓ |
| Hexanal | 0.000203 | 0.000203 | 0.000203 | 0.000203 | 0.00 | ✓ |
| Nonanal | 0.000105 | 0.000105 | 0.000105 | 0.000105 | 0.00 | ✓ |

### `spi_hvp_xylose_120C_PMC9905368`

- Execution path: `free_precursor`
- Matched compounds with envelope: 3

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| 2-Furfurylthiol (FFT) | 0.42 | 0.00173 | 0.672 | 4.91 | 3.45 | ✓ |
| 2-Methyl-3-furanthiol (MFT) | 0.18 | 3.01e-05 | 0.247 | 40.6 | 6.13 | ✓ |
| Methional | 1.88 | 0.0908 | 1.57 | 2.09 | 1.36 | ✓ |

### `thiamine_cys_ribose_100C_Hofmann1996`

- Execution path: `free_precursor`
- Matched compounds with envelope: 1

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| 2-Methyl-3-furanthiol (MFT) | 3.14 | 0.0017 | 3.69 | 20.2 | 4.08 | ✓ |

### `thiamine_cys_xylose_145C_Cerny2008`

- Execution path: `free_precursor`
- Matched compounds with envelope: 1

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| 2-Methyl-3-furanthiol (MFT) | 2.47 | 0.0028 | 2.98 | 21.4 | 3.88 | ✓ |

### `wheat_gluten_hvp_xylose_120C_PMC9905368`

- Execution path: `free_precursor`
- Matched compounds with envelope: 3

| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |
| --- | --- | --- | --- | --- | --- | --- |
| 2-Furfurylthiol (FFT) | 0.61 | 0.00283 | 0.856 | 8.09 | 3.46 | ✓ |
| 2-Methyl-3-furanthiol (MFT) | 0.34 | 5.07e-05 | 0.407 | 42.5 | 5.92 | ✓ |
| Methional | 3.44 | 0.127 | 2.83 | 3.91 | 1.49 | ✓ |
