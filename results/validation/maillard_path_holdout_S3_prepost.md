# Wave S3 — hold-out pre/post against the pre-registered baseline

Generated 2026-08-27. Baseline: `results/validation/maillard_path_holdout_frozen_predictions.json` at commit `12f43dd`. READ, never regenerated. The pre-registration is the file, not a re-run.

## AS SHIPPED — what this wave actually changes

**0/22 targets moved** (0 improved, 0 worsened).

| measure | frozen | shipped |
|---|---:|---:|
| `median_fold_error` | 6.0388 | 6.0388 |
| `median_abs_log10_error` | 0.7809 | 0.7809 |
| `within_10x` | 12 | 12 |
| `worst_fold_error` | 52.5864 | 52.5864 |
| `best_fold_error` | 1.5168 | 1.5168 |
| `structural_zero_count` | 1 | 1 |
| `ordinal_pairs_correct` | 8 | 8 |
| `series_directions_correct` | 3 | 3 |

### Series directions (the three Wave U structural errors live here)

| series | compound | measured ratio | frozen predicted | shipped predicted | frozen ok? | shipped ok? |
|---|---|---:|---:|---:|---|---|
| acrylamide_time_series_Chang2021 | Acrylamide | 52.11 | 1.24 | 1.24 | True | True |
| ph_pair_Schibilsky2019 | 5-Hydroxymethylfurfural (HMF) | 1.758 | 0.7743 | 0.7743 | False | False |
| ph_pair_Schibilsky2019 | DMHF | 5.111 | 2.479 | 2.479 | True | True |
| ph_pair_Schibilsky2019 | Furfural | 1.471 | 0.7743 | 0.7743 | False | False |
| sulfur_temperature_series_Yiltirak2026 | 2-Furfurylthiol (FFT) | 1.266 | 5.475 | 5.475 | True | True |
| sulfur_temperature_series_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 0.2485 | 4.551 | 4.551 | False | False |

**No target moved. Every prediction is bit-identical to the pre-registration.**

## COUNTERFACTUAL — what applying the derived barriers WOULD do (measured; not shipped)

Barrier overrides applied in memory only:

* `amadori_rearrangement` → 23.20 kcal/mol
* `enolisation_intermediate` → 25.55 kcal/mol
* `2,3-enolisation` → 24.40 kcal/mol
* `1,2-enolisation` → 23.49 kcal/mol
* `furanone_amino_acid_reduction` → 20.99 kcal/mol

**21/22 targets moved** (11 improved, 10 worsened).

| measure | frozen | counterfactual |
|---|---:|---:|
| `median_fold_error` | 6.0388 | 7.6110 |
| `median_abs_log10_error` | 0.7809 | 0.8814 |
| `within_10x` | 12 | 13 |
| `worst_fold_error` | 52.5864 | 55.4654 |
| `best_fold_error` | 1.5168 | 1.1215 |
| `structural_zero_count` | 1 | 1 |
| `ordinal_pairs_correct` | 8 | 6 |
| `series_directions_correct` | 3 | 3 |

### Series directions (the three Wave U structural errors live here)

| series | compound | measured ratio | frozen predicted | counterfactual predicted | frozen ok? | counterfactual ok? |
|---|---|---:|---:|---:|---|---|
| acrylamide_time_series_Chang2021 | Acrylamide | 52.11 | 1.24 | 1.493 | True | True |
| ph_pair_Schibilsky2019 | 5-Hydroxymethylfurfural (HMF) | 1.758 | 0.7743 | 0.3008 | False | False |
| ph_pair_Schibilsky2019 | DMHF | 5.111 | 2.479 | 2.427 | True | True |
| ph_pair_Schibilsky2019 | Furfural | 1.471 | 0.7743 | 0.3008 | False | False |
| sulfur_temperature_series_Yiltirak2026 | 2-Furfurylthiol (FFT) | 1.266 | 5.475 | 8.715 | True | True |
| sulfur_temperature_series_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 0.2485 | 4.551 | 7.245 | False | False |

### Targets that moved

| benchmark | compound | target | frozen fold | counterfactual fold |
|---|---|---:|---:|---:|
| mp_holdout_fructose_asparagine_180C_Lin2022 | 5-Hydroxymethylfurfural (HMF) | 1.228e+04 ppb | 4.984 | 3.833 |
| mp_holdout_fructose_asparagine_180C_Lin2022 | Acrylamide | 1859 ppb | 18.49 | 24.04 |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | 5-Hydroxymethylfurfural (HMF) | 5.725e+04 ppb | 2.363 | 3.084 |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | DMHF | 1153 ppb | 6.039 | 16.03 |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | Furfural | 1633 ppb | 11.3 | 8.66 |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | 5-Hydroxymethylfurfural (HMF) | 1.006e+05 ppb | 5.363 | 18.02 |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | DMHF | 5894 ppb | 2.929 | 7.611 |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | Furfural | 2402 ppb | 5.952 | 1.771 |
| mp_holdout_glucose_asparagine_180C_10min_Chang2021 | Acrylamide | 28 ppb | 2.19 | 1.589 |
| mp_holdout_glucose_asparagine_180C_30min_Chang2021 | Acrylamide | 1459 ppb | 19.19 | 55.47 |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | 5-Hydroxymethylfurfural (HMF) | 7000 ppb | 6.161 | 6.167 |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | Acrylamide | 832 ppb | 10.94 | 31.63 |
| mp_holdout_glucose_asparagine_180C_Ye2024 | Acrylamide | 140.6 umol_per_mol_limiting_precursor | 1.749 | 40.78 |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 1.28 ppb | 12.16 | 4.41 |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 6.88 ppb | 1.517 | 4.181 |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 1.46 ppb | 22.96 | 9.977 |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 3.29 ppb | 2.581 | 1.122 |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 1.68 ppb | 36.95 | 18.78 |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 2.4 ppb | 5.989 | 3.045 |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Furfurylthiol (FFT) | 1.62 ppb | 52.59 | 30.37 |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 1.71 ppb | 12.07 | 6.972 |

