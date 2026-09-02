# Both frozen hold-outs, before and after the barrier-offset retirement (2026-08-28)

The auto-accepted +/-3.0 kcal/mol barrier offsets in data/lit/refinement_surrogate_patches.json were retired on 2026-08-28. They had been applied by src.barrier_constants.get_barrier() to every shipped prediction, silently overriding the audited FAST_BARRIERS table. The 'before' arm reconstructs them through the BARRIER_OFFSETS environment variable in a subprocess; no file was restored and no hold-out was regenerated or re-frozen.

## Shipped barriers (kcal/mol)

| family | before (offset applied) | after (= FAST_BARRIERS table) | delta |
|---|---:|---:|---:|
| `Schiff_Base_Formation` | 18.0 | 15.0 | -3.00 |
| `Retro_Aldol_Fragmentation` | 29.0 | 32.0 | +3.00 |
| `Thiol_Addition` | 31.6 | 28.6 | -3.00 |
| `Thiol_Addition_Pentodiulose` | 28.6 | 28.6 | +0.00 |
| `Thiol_Addition_Hexose` | 29.65 | 29.65 | +0.00 |
| `Amadori_Rearrangement` | 23.0 | 23.0 | +0.00 |

## Hold-out 1 — the free-precursor `maillard_path` pre-registration

Baseline `results/validation/maillard_path_holdout_frozen_predictions.json` at commit `12f43dd`: READ, never regenerated, never re-frozen.

**31 of 32 targets moved** when the offsets were retired.

| metric | frozen (12f43dd) | before retirement | after retirement |
|---|---:|---:|---:|
| `bundle_count` | 12 | 17 | 17 |
| `target_count` | 22 | 32 | 32 |
| `quantitatively_scored_count` | 21 | 31 | 31 |
| `structural_zero_count` | 1 | 1 | 1 |
| `median_fold_error` | 6.0388 | 10.0462 | 10.8638 |
| `median_abs_log10_error` | 0.7809 | 1.0020 | 1.0360 |
| `worst_fold_error` | 52.5864 | 565.2249 | 506.4244 |
| `best_fold_error` | 1.5168 | 1.2516 | 1.2435 |
| `within_10x` | 12 | 15 | 15 |
| `ordinal_pairs_correct` | 8 | 10 | 10 |
| `series_directions_correct` | 3 | 3 | 3 |

Against the pre-registration's own 22 targets: **21 differed before the retirement, 8 differ after**.

### Targets the retirement moved

| target | before | after | after/before |
|---|---:|---:|---:|
| `mp_holdout_fructose_asparagine_180C_Lin2022 / 5-Hydroxymethylfurfural (HMF)` | 61256.5885 | 61206.6652 | 0.9992x |
| `mp_holdout_fructose_asparagine_180C_Lin2022 / Acrylamide` | 100.6405 | 100.5585 | 0.9992x |
| `mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 / 5-Hydroxymethylfurfural (HMF)` | 24224.7473 | 24232.9737 | 1.0003x |
| `mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 / DMHF` | 6980.4149 | 6963.7087 | 0.9976x |
| `mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 / Furfural` | 18457.0327 | 18463.3004 | 1.0003x |
| `mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 / 5-Hydroxymethylfurfural (HMF)` | 18763.6634 | 18763.6690 | 1.0000x |
| `mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 / DMHF` | 17263.1400 | 17263.1374 | 1.0000x |
| `mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 / Furfural` | 14296.1883 | 14296.1925 | 1.0000x |
| `mp_holdout_glucose_asparagine_180C_10min_Chang2021 / Acrylamide` | 61.6326 | 61.3208 | 0.9949x |
| `mp_holdout_glucose_asparagine_180C_30min_Chang2021 / Acrylamide` | 76.3805 | 76.0366 | 0.9955x |
| `mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 / 5-Hydroxymethylfurfural (HMF)` | 43127.9604 | 43128.2655 | 1.0000x |
| `mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 / Acrylamide` | 76.3805 | 76.0366 | 0.9955x |
| `mp_holdout_glucose_asparagine_180C_Ye2024 / Acrylamide` | 80.4224 | 80.3778 | 0.9994x |
| `mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 / 2-Furfurylthiol (FFT)` | 668.9628 | 951.2337 | 1.4220x |
| `mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 / 2-Methyl-3-furanthiol (MFT)` | 1695.6748 | 1519.2732 | 0.8960x |
| `mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 / 2-Furfurylthiol (FFT)` | 286.0801 | 343.3847 | 1.2003x |
| `mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 / 2-Methyl-3-furanthiol (MFT)` | 15.7394 | 22.8502 | 1.4518x |
| `mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 / 2-Furfurylthiol (FFT)` | 1021.2252 | 1408.6073 | 1.3793x |
| `mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 / 2-Methyl-3-furanthiol (MFT)` | 3079.7911 | 2877.6790 | 0.9344x |
| `mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 / 2-Furfurylthiol (FFT)` | 523.2364 | 624.7049 | 1.1939x |
| `mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 / 2-Methyl-3-furanthiol (MFT)` | 321.3562 | 316.2148 | 0.9840x |
| `mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 / 2-Furfurylthiol (FFT)` | 1168.3002 | 1428.7610 | 1.2229x |
| `mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 / 2-Methyl-3-furanthiol (MFT)` | 750.8980 | 713.7397 | 0.9505x |
| `mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 / 2-Furfurylthiol (FFT)` | 12.8592 | 15.4393 | 1.2006x |
| `mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 / 2-Methyl-3-furanthiol (MFT)` | 8.6111 | 8.5556 | 0.9936x |
| `mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 / 2-Furfurylthiol (FFT)` | 27.7062 | 33.1748 | 1.1974x |
| `mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 / 2-Methyl-3-furanthiol (MFT)` | 15.8614 | 15.7090 | 0.9904x |
| `mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 / 2-Furfurylthiol (FFT)` | 51.3362 | 61.2810 | 1.1937x |
| `mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 / 2-Methyl-3-furanthiol (MFT)` | 26.4237 | 26.0732 | 0.9867x |
| `mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 / 2-Furfurylthiol (FFT)` | 70.4328 | 83.9927 | 1.1925x |
| `mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 / 2-Methyl-3-furanthiol (MFT)` | 37.2963 | 36.7348 | 0.9849x |

## Hold-out 2 — the eight-point matrix hold-out

**0 of 8 points moved.**

| metric | before retirement | after retirement |
|---|---:|---:|
| `matched_compound_count` | 8 | 8 |
| `ci_coverage_hits` | 3 | 3 |
| `ci_coverage_rate` | 0.3750 | 0.3750 |
| `median_accuracy_fold` | 93.6837 | 93.6837 |
| `median_abs_log10_error` | 1.9717 | 1.9717 |
| `max_fold_error` | 2474.3860 | 2474.3860 |

**No point moved.** The matrix hold-out runs the `matrix_only` execution path, which never reaches the reaction network (Wave S1), so it is structurally blind to any barrier. That is a NEGATIVE result about the hold-out, not a clean bill of health for the model.

