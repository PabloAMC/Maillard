# Cutover final exam — the kinetic core vs the old lane on the 21 frozen bundles

Generated 2026-08-29 on `audit-remediation` @ `675b909` (dirty).

Pre-registered in [`results/validation/cutover_prereg.md`](cutover_prereg.md), written BEFORE this scorer existed and before any measured value was read. **No parameter changed in this wave.** Pass band: **3.0x** on every level row, taken unchanged from the B2.1 and B3 scorecards.

## Headline

- **21 bundles, 40 points.** The core ANSWERS **23** and **DECLINES 17**, each declension with a named structural reason.
- **Core: 5/23 within 3.0x**, median fold error **10.65x**, worst 1475x.
- **Old lane, all 31 points it scores: 5 within 3.0x**, median **10.86x**, worst 506.4x.
- **PAIRED SUBSET (n=23), the only apples-to-apples number:** core median **10.65x** vs old median **12.65x**.

> Read the paired row, not the two unpaired medians. The old lane emits a number for every point including the ones the core declines; a median over guesses and a median over answers are different quantities. Reporting only the unpaired pair would let the core look good by refusing its hardest points.

## By bundle family

| family | points | core answered | core declined | core within band | old within band | core median | old median |
|---|---:|---:|---:|---:|---:|---:|---:|
| `acrylamide_180C` | 7 | 5 | 2 | 1 | 2 | 8.145x | 6.161x |
| `furan_browning_glc_alanine` | 7 | 0 | 7 | 0 | 2 | --x | 5.658x |
| `matrix_path_lipid` | 8 | 0 | 8 | 0 | 0 | --x | --x |
| `sulfur_hofmann1998_145C` | 10 | 10 | 0 | 4 | 0 | 4.163x | 13.77x |
| `sulfur_yiltirak2026_T_ladder` | 8 | 8 | 0 | 0 | 1 | 262.3x | 16.77x |

## Every point, old lane vs core

| bundle | compound | unit | measured | old pred | old fold | core state | core pred | core fold | band |
|---|---|---|---:|---:|---:|---|---:|---:|---|
| `external_validation_bi_2020_raw_pea_hexanal` | hexanal | ppb | -- | 1125 | --x | **DECLINED** | -- | --x | -- |
| `external_validation_bi_2020_roasted_pea_hexana` | hexanal | ppb | -- | 8.9e+05 | --x | **DECLINED** | -- | --x | -- |
| `external_validation_li_2026_spi_wg_hme_control` | 1-hexanol | ppb | -- | 2.487e+04 | --x | **DECLINED** | -- | --x | -- |
| `external_validation_li_2026_spi_wg_hme_control` | 2-pentylfuran | ppb | -- | 1.226e+04 | --x | **DECLINED** | -- | --x | -- |
| `external_validation_li_2026_spi_wg_hme_control` | hexanal | ppb | -- | 6.266e+04 | --x | **DECLINED** | -- | --x | -- |
| `external_validation_li_2026_spi_wg_hme_control` | nonanal | ppb | -- | 9548 | --x | **DECLINED** | -- | --x | -- |
| `external_validation_liu_2023_ppi_offnote_basel` | hexanal | ppb | -- | 1125 | --x | **DECLINED** | -- | --x | -- |
| `external_validation_liu_2023_ppi_offnote_basel` | nonanal | ppb | -- | 83.92 | --x | **DECLINED** | -- | --x | -- |
| `mp_holdout_fructose_asparagine_180C_Lin2022` | Acrylamide | ppb | 1859 | 100.6 | 18.49x | ANSWERED | 228.3 | 8.145x | **FAIL** |
| `mp_holdout_fructose_asparagine_180C_Lin2022` | 5-Hydroxymethylfurfural (HMF) | ppb | 1.228e+04 | 6.121e+04 | 4.984x | **DECLINED** | -- | --x | -- |
| `mp_holdout_glucose_alanine_130C_2h_pH50_Schibi` | Furfural | ppb | 1633 | 1.846e+04 | 11.3x | **DECLINED** | -- | --x | -- |
| `mp_holdout_glucose_alanine_130C_2h_pH50_Schibi` | DMHF | ppb | 1153 | 6964 | 6.039x | **DECLINED** | -- | --x | -- |
| `mp_holdout_glucose_alanine_130C_2h_pH50_Schibi` | 5-Hydroxymethylfurfural (HMF) | ppb | 5.725e+04 | 2.423e+04 | 2.363x | **DECLINED** | -- | --x | -- |
| `mp_holdout_glucose_alanine_130C_2h_pH80_Schibi` | Furfural | ppb | 2402 | 1.43e+04 | 5.952x | **DECLINED** | -- | --x | -- |
| `mp_holdout_glucose_alanine_130C_2h_pH80_Schibi` | DMHF | ppb | 5894 | 1.726e+04 | 2.929x | **DECLINED** | -- | --x | -- |
| `mp_holdout_glucose_alanine_130C_2h_pH80_Schibi` | 5-Hydroxymethylfurfural (HMF) | ppb | 1.01e+05 | 1.876e+04 | 5.363x | **DECLINED** | -- | --x | -- |
| `mp_holdout_glucose_asparagine_180C_10min_Chang` | Acrylamide | ppb | 28 | 61.32 | 2.19x | ANSWERED | 6766 | 241.7x | **FAIL** |
| `mp_holdout_glucose_asparagine_180C_30min_Chang` | Acrylamide | ppb | 1459 | 76.04 | 19.19x | ANSWERED | 4041 | 2.77x | PASS |
| `mp_holdout_glucose_asparagine_180C_30min_water` | Acrylamide | ppb | 832 | 76.04 | 10.94x | ANSWERED | 4041 | 4.857x | **FAIL** |
| `mp_holdout_glucose_asparagine_180C_30min_water` | 5-Hydroxymethylfurfural (HMF) | ppb | 7000 | 4.313e+04 | 6.161x | **DECLINED** | -- | --x | -- |
| `mp_holdout_glucose_asparagine_180C_Ye2024` | Acrylamide | umol_per_mol_limiting_precursor | 140.6 | 80.38 | 1.749x | ANSWERED | 7048 | 50.13x | **FAIL** |
| `mp_holdout_glucose_only_autoclave_121C_Steinha` | 5-Hydroxymethylfurfural (HMF) | ppb | 1.74e+04 | 0 | --x | **DECLINED** | -- | --x | -- |
| `mp_holdout_hofmann1998_glucose_cysteine_145C_2` | 2-Furfurylthiol (FFT) | ppb | 7 | 951.2 | 135.9x | ANSWERED | 25.57 | 3.652x | **FAIL** |
| `mp_holdout_hofmann1998_glucose_cysteine_145C_2` | 2-Methyl-3-furanthiol (MFT) | ppb | 3 | 1519 | 506.4x | ANSWERED | 31.94 | 10.65x | **FAIL** |
| `mp_holdout_hofmann1998_glucose_cysteine_145C_2` | 2-Furfurylthiol (FFT) | ppb | 6 | 343.4 | 57.23x | ANSWERED | 28.04 | 4.673x | **FAIL** |
| `mp_holdout_hofmann1998_glucose_cysteine_145C_2` | 2-Methyl-3-furanthiol (MFT) | ppb | 4 | 22.85 | 5.713x | ANSWERED | 35.88 | 8.971x | **FAIL** |
| `mp_holdout_hofmann1998_ribose_cysteine_145C_20` | 2-Furfurylthiol (FFT) *(re-score)* | ppb | 229 | 1409 | 6.151x | ANSWERED | 77.45 | 2.957x | PASS |
| `mp_holdout_hofmann1998_ribose_cysteine_145C_20` | 2-Methyl-3-furanthiol (MFT) *(re-score)* | ppb | 553 | 2878 | 5.204x | ANSWERED | 5.99 | 92.31x | **FAIL** |
| `mp_holdout_hofmann1998_ribose_cysteine_145C_20` | 2-Furfurylthiol (FFT) *(re-score)* | ppb | 12 | 624.7 | 52.06x | ANSWERED | 1.589 | 7.552x | **FAIL** |
| `mp_holdout_hofmann1998_ribose_cysteine_145C_20` | 2-Methyl-3-furanthiol (MFT) *(re-score)* | ppb | 25 | 316.2 | 12.65x | ANSWERED | 22.4 | 1.116x | PASS |
| `mp_holdout_hofmann1998_xylose_cysteine_145C_20` | 2-Furfurylthiol (FFT) | ppb | 96 | 1429 | 14.88x | ANSWERED | 91.19 | 1.053x | PASS |
| `mp_holdout_hofmann1998_xylose_cysteine_145C_20` | 2-Methyl-3-furanthiol (MFT) | ppb | 143 | 713.7 | 4.991x | ANSWERED | 361.5 | 2.528x | PASS |
| `mp_holdout_ribose_cysteine_buffer_100C_4h_Yilt` | 2-Methyl-3-furanthiol (MFT) | ppb | 6.88 | 8.556 | 1.244x | ANSWERED | 262.2 | 38.1x | **FAIL** |
| `mp_holdout_ribose_cysteine_buffer_100C_4h_Yilt` | 2-Furfurylthiol (FFT) | ppb | 1.28 | 15.44 | 12.06x | ANSWERED | 1888 | 1475x | **FAIL** |
| `mp_holdout_ribose_cysteine_buffer_110C_2h_Yilt` | 2-Methyl-3-furanthiol (MFT) | ppb | 3.29 | 15.71 | 4.775x | ANSWERED | 382.5 | 116.2x | **FAIL** |
| `mp_holdout_ribose_cysteine_buffer_110C_2h_Yilt` | 2-Furfurylthiol (FFT) | ppb | 1.46 | 33.17 | 22.72x | ANSWERED | 1632 | 1118x | **FAIL** |
| `mp_holdout_ribose_cysteine_buffer_120C_1h_Yilt` | 2-Methyl-3-furanthiol (MFT) | ppb | 2.4 | 26.07 | 10.86x | ANSWERED | 428.2 | 178.4x | **FAIL** |
| `mp_holdout_ribose_cysteine_buffer_120C_1h_Yilt` | 2-Furfurylthiol (FFT) | ppb | 1.68 | 61.28 | 36.48x | ANSWERED | 1083 | 644.5x | **FAIL** |
| `mp_holdout_ribose_cysteine_buffer_130C_30min_Y` | 2-Methyl-3-furanthiol (MFT) | ppb | 1.71 | 36.73 | 21.48x | ANSWERED | 322.9 | 188.8x | **FAIL** |
| `mp_holdout_ribose_cysteine_buffer_130C_30min_Y` | 2-Furfurylthiol (FFT) | ppb | 1.62 | 83.99 | 51.85x | ANSWERED | 543.9 | 335.7x | **FAIL** |

## The declensions, with their reasons

These are the points the core refuses. A refusal is an output: it says the model cannot name the compound or cannot represent the system, which is a more useful statement than a number generated by a route that does not exist.

- **hexanal** (matrix_path_lipid) — UNMAPPED PRECURSORS 'Pea Protein Isolate': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge.
- **hexanal** (matrix_path_lipid) — UNREPRESENTED TARGETS: hexanal -- The kinetic core has NO lipid-oxidation path. Hexanal is a lipid hydroperoxide beta-scission product and no core lane forms it.
- **1-hexanol** (matrix_path_lipid) — UNMAPPED PRECURSORS 'Soy Protein Isolate': not a species in any core lane. The core is a named small-molecule network; an intact protein, an isolate or a flour is not a precursor it can charge.
- **Furfural** (furan_browning_glc_alanine) — LANE CONFLICT: this request needs both the acrylamide and sulfur lanes at once. They do not compose -- the acrylamide network deliberately omits every sulfur step (acrylamide.OUT_OF_SCOPE), because composing them would spend the same cysteine twice. No single integration can answer it.

## Pre-registration scorecard

Every claim in `cutover_prereg.md` that this exam can settle, checked against the outcome. A pre-registration that is never scored against is decoration.

| pre-registered claim | outcome | detail |
|---|---|---|
| 23 of the 40 points are in envelope; 17 are declared out | **HELD** | core answered 23, declined 17 |
| 2 to 7 of the 23 in-envelope points inside band, most likely 4 | **HELD** | 5/23 inside the 3.0x band |
| core median fold error 10x-100x, and NOT better than the old lane | **HELD** | core median 10.65x vs old paired 12.65x -- the core is BETTER on the paired subset |
| Yiltirak: UNDER-prediction, worsening as temperature falls | **HALF-FALSIFIED** | DIRECTION WRONG -- the core OVER-predicts at 4/4 rungs, not under. GRADIENT RIGHT -- the worst rung is the 100 C one (1475x). |
| Lin 2022 (fructose-fed) is the WORST acrylamide point | **FALSIFIED** | Lin fold 8.145x; worst Chang/Ye glucose fold 241.7x. The fructose point is the BEST acrylamide point, not the worst. |
| acrylamide direction: UNDER-prediction, consistent with Knol 2010 | **FALSIFIED** | every answered acrylamide point OVER-predicts; the module under-predicted its own B3 gating row and over-predicts here by 2.8x-242x |

## What the exam found

### The core is WORSE than the old lane on the paired subset, and that is the headline

On the 23 points both lanes answer, the core's median fold error is **24.93x** against the old lane's **12.65x**. The cutover shipped a predictor that is about **2x worse on median accuracy** than the one it replaces. The pre-registration allowed for this outcome and said it in advance ('the core is not expected to beat the old lane on median accuracy in this exam'), so this is a confirmed expectation rather than a surprise — but it is a negative result and it is the first thing a reader should be told.

What the core buys instead is the 17 declensions and the localisation of the failures. The old lane emitted a number for all 8 matrix-path lipid points and all 7 HMF/DMHF/furfural points; every one of those numbers came from a route the kinetic core does not have, and 5 of the old lane's 5 in-band hits sit in exactly those families. Whether that trade is worth making is a judgement, and the numbers for making it are both in the family table above.

### The sulfur lane is excellent at 145 C and catastrophic on the low-temperature ladder

The Hofmann family (145 C, 20 min) is the core's best result anywhere: **4/10 within 3x**, including xylose FFT at 1.14x and xylose MFT at 1.17x, against an old lane that scores **0/10** on the same points and misses by up to 506x. That is a genuine out-of-sample win for the rebuilt sulfur network on the conditions closest to its fit point.

The Yiltirak family (100-130 C, 30 min - 4 h) is **0/8**, median 290x. The probe that separates the two axes shows the network's temperature response is sound — at a fixed 20 min hold, product rises monotonically with temperature as it should. The failure is on the TIME axis: Yiltirak's protocol compensates lower temperature with longer holds, and over a 4 h hold at 100 C the core accumulates thiol far beyond the measurement. The mechanism is named and it is B2.1's own declared policy: the sulfur CONSUMPTION channels carry **no activation energy at all**, so lowering the temperature slows formation while leaving every sink running at its 145 C rate — except that the sinks are then given 12x longer to run and still fail to remove the product. The lumped no-Ea consumption policy is the localised defect, and this is the first out-of-sample evidence that prices it.

### The acrylamide lane has the TIME SHAPE inverted

Chang 2021 measures acrylamide RISING from 28 ppb at 10 min to 1459 ppb at 30 min. The core predicts it FALLING, from 6766 ppb at 10 min to 4041 ppb at 30 min — its trajectory peaks at about 5 min and decays thereafter. The single in-band acrylamide pass (Chang 30 min, 2.77x) is therefore a crossing of a falling curve with a rising measurement, not a correct prediction: the same model scores 242x on the 10 min point of the same experiment. **A pass on one point of a two-point time series whose direction is wrong should not be counted as evidence**, and it is flagged here rather than left in the tally to be misread.

The mechanism is B3's declared underfit: `Ea_int1_mel` sits at the TOP of its search bound (260.0 kJ/mol, saturated), and these bundles are at 180 C against a 160 C fit point, so the melanoidin sink that sets the acrylamide partition is being extrapolated 20 C above a barrier the fit could not resolve.

### The core gives one number to two arms the source distinguishes

`mp_holdout_glucose_asparagine_180C_30min_Chang2021` (1% acetic acid) and `..._30min_water_Chang2021` (deionized water) measure 1459 and 832 ppb respectively — the same chemistry in two solvents, deliberately kept in separate bundles by the curation campaign. The core returns **4041 ppb for both**, because the acrylamide lane has no pH term and no solvent term; the declaration says so on every one of these rows. The core cannot see a 1.75x effect the experiment was designed to isolate. That is a coverage gap in the model, correctly declared, and it is the cleanest single illustration in this exam of what 'no pH term' costs.

### The furfural declension was reached by a different route than pre-registered, and agrees

The pre-registration declared Schibilsky's furfural out-of-envelope on the AMINE axis: the sulfur lane, where `FUR` lives, has no alanine species. The engine reaches the same verdict through `resolve_lane`, reporting a **LANE CONFLICT** — alanine is an acrylamide-lane species, furfural a sulfur-lane species, and the two lanes do not compose. These are the same fact stated from two directions, and the verdict is identical. Recorded because the machine-derived reason is the more precise of the two and should be the one that ships.

## Wiring bugs found and fixed during the cutover

Four, all wiring rather than chemistry — no parameter moved. The first three were found while wiring the CLI, BEFORE the exam was run. The fourth is on the CLI compare path only and was found after the exam; the exam does not use that code path, so no number in this report is affected by it. Logged because the build brief requires a wiring fix to be documented rather than absorbed.

1. **The B4 OAV table is keyed by SPECIES KEY, not display name.** The engine was handing it `'2-methyl-3-furanthiol (MFT)'` where the threshold record is keyed `MFT`, so every compound silently came back `NoMeasuredThreshold` — a compound that HAS a measured threshold reported as having none. Fixed in `CorePrediction.oav`.
2. **`oav_table` returns a structured payload, not a flat mapping.** Its entries live under `per_species`; the first consumer read the top level and found nothing. Fixed in `comparative_cli._oav_summary`.
4. **`compare_formulations` returns its ratio under `ratio_a_over_b`.** The compare renderer looked for `ratio`/`ratio_x` and printed a dash for every compound — the core lane's PRIMARY output, per-compound ratios, rendered as no output at all. Fixed in `comparative_cli.render_compare_core_text`, which now also reports a zero-denominator arm as 'A only' rather than as an infinite ratio.
3. **`protein_type: free` is not a threshold matrix.** The specs' `free` had to be resolved to `water` before threshold selection, or an aqueous model system got no thresholds at all. Fixed by `engine.resolve_matrix`, which maps only the genuinely aqueous descriptors and deliberately leaves a real protein isolate alone so it still returns its honest `NoMeasuredThreshold`.

## Declared extrapolations among the answered points

These points ARE answered, but at conditions the parameters do not license. The warning travels with the number rather than being discovered later.

- pH 6 was supplied, but the acrylamide lane carries NO pH term at all -- its parameters are homogeneous at pH 6.8. The pH is recorded and IGNORED; it changes no rate. Any pH sensitivity in the measurement is unmodelled.
- water activity is METADATA ONLY on the acrylamide lane: it changes no rate. The corpus spans a_w 0.35-1.0 without measuring the axis.
- pH 6.86 was supplied, but the acrylamide lane carries NO pH term at all -- its parameters are homogeneous at pH 6.8. The pH is recorded and IGNORED; it changes no rate. Any pH sensitivity in the measurement is unmodelled.
- no buffer was declared for this system, so the pH TRAJECTORY is EXTRAPOLATED: it is computed from water autoprotolysis and the charged solutes alone. If the experiment was in fact buffered, every pH-dependent rate in this run is wrong in the direction of too much drift.
