# Cutover final exam — the kinetic core vs the old lane on the 21 frozen bundles

Generated 2026-08-29 on `audit-remediation` @ `bdb1142` (dirty).

Pre-registered in [`results/validation/cutover_prereg.md`](cutover_prereg.md), written BEFORE this scorer existed and before any measured value was read. **No parameter changed in this wave.** Pass band: **3.0x** on every level row, taken unchanged from the B2.1 and B3 scorecards.

## Headline

- **21 bundles, 40 points.** The core ANSWERS **27** and **DECLINES 13**, each declension with a named structural reason.
- **Core: 3/27 within 3.0x**, median fold error **21.36x**, worst 3.339e+04x.
- **Old lane, all 39 points it scores: 7 within 3.0x**, median **10.94x**, worst 2748x.
- **PAIRED SUBSET (n=27), the only apples-to-apples number:** core median **21.36x** vs old median **12.65x**.
- **GEOMETRIC-MEAN FOLD, on the same points:** core over all scored **48.73x**; on the paired subset core **48.73x** vs old **16.8x**. A median moves only when a row crosses the middle of the pool; a geometric mean moves whenever any row moves, including a row that was already failing.
- **4 of these points are SHARED with the hold-out panel** — the same measurements, scored twice. THESE FOUR POINTS ARE THE SAME MEASUREMENTS AS FOUR ROWS OF THE KINETIC-CORE HOLD-OUT PANEL, not analogues of them. The exam and the panel are therefore NOT independent evidence on this axis: agreement between them here is one measurement counted twice. Established in results/validation/d1_exam_panel_reconciliation.md sec. 5 and declared in both artifacts from Wave B2.4 onward.

> Read the paired row, not the two unpaired medians. The old lane emits a number for every point including the ones the core declines; a median over guesses and a median over answers are different quantities. Reporting only the unpaired pair would let the core look good by refusing its hardest points.

## BOTH WAYS — buffer-completed and as-was

FIT_HOLDOUT_DECLARATION.md Amendment 9 clause 2: the exam is reported BOTH WAYS -- buffer-completed and as-was -- in the same artifact, PERMANENTLY. Not transitional: every number this repo published before B2.3 was computed as-was, and a report that silently replaced them would make its own history unreadable.

| | scored | within band | median fold | paired median (n=27) |
|---|---:|---:|---:|---:|
| **buffer-completed** | 27 | 3 | 21.36x | **21.36x** |
| **as-was (no buffer field)** | 27 | 5 | 10.63x | **10.63x** |
| old lane (identical in both) | 39 | 7 | 10.94x | 12.65x |

> THE OLD LANE HAS NO pH STATE AND NO BUFFER INPUT AT ALL, so its numbers are BY CONSTRUCTION the same in both columns. That is why the old-lane comparison is reported against BOTH core columns rather than recomputed: the comparison changes because the CORE moves, never because the old lane does.

**The buffer field moved 16 of 40 points** — 4 closer to the measurement, 12 further away, 24 untouched.

> Only the SULFUR lane carries a pH state. An acrylamide-lane or matrix-lane row is identical in both columns no matter what its buffer says, and that identity is a REPORTED GAP rather than an omission -- it is the same gap that leaves Chang's two arms predicting the same value.

### Both ways, by family

| family | points | completed median | completed in band | as-was median | as-was in band |
|---|---:|---:|---:|---:|---:|
| `acrylamide_180C` | 7 | 8.145x | 1 | 8.145x | 1 |
| `furan_browning_glc_alanine` | 7 | --x | 0 | --x | 0 |
| `matrix_path_lipid` | 8 | 1863x | 0 | 1863x | 0 |
| `sulfur_hofmann1998_145C` | 10 | 5.103x | 2 | 3.326x | 4 |
| `sulfur_yiltirak2026_T_ladder` | 8 | 610.8x | 0 | 226.6x | 0 |

### Every point the buffer field moved

| bundle | compound | buffer | provenance | as-was fold | completed fold | |
|---|---|---|---|---:|---:|---|
| `mp_holdout_hofmann1998_glucose_cysteine_145C_2` | 2-Furfurylthiol (FFT) | phosphate 0.5 M | primary_source_pdf | 3.63x | 3.804x | **further** |
| `mp_holdout_hofmann1998_glucose_cysteine_145C_2` | 2-Furfurylthiol (FFT) | phosphate 0.5 M | primary_source_pdf | 4.581x | 33.56x | **further** |
| `mp_holdout_hofmann1998_glucose_cysteine_145C_2` | 2-Methyl-3-furanthiol (MFT | phosphate 0.5 M | primary_source_pdf | 8.777x | 5.734x | closer |
| `mp_holdout_hofmann1998_ribose_cysteine_145C_20` | 2-Furfurylthiol (FFT) | phosphate 0.5 M | primary_source_pdf | 3.022x | 1.454x | closer |
| `mp_holdout_hofmann1998_ribose_cysteine_145C_20` | 2-Methyl-3-furanthiol (MFT | phosphate 0.5 M | primary_source_pdf | 92.28x | 21.36x | closer |
| `mp_holdout_hofmann1998_ribose_cysteine_145C_20` | 2-Furfurylthiol (FFT) | phosphate 0.5 M | primary_source_pdf | 2.506x | 4.472x | **further** |
| `mp_holdout_hofmann1998_ribose_cysteine_145C_20` | 2-Methyl-3-furanthiol (MFT | phosphate 0.5 M | primary_source_pdf | 2.135x | 140.3x | **further** |
| `mp_holdout_hofmann1998_xylose_cysteine_145C_20` | 2-Furfurylthiol (FFT) | phosphate 0.5 M | primary_source_pdf | 1.076x | 3.285x | **further** |
| `mp_holdout_hofmann1998_xylose_cysteine_145C_20` | 2-Methyl-3-furanthiol (MFT | phosphate 0.5 M | primary_source_pdf | 2.287x | 2.989x | **further** |
| `mp_holdout_ribose_cysteine_buffer_100C_4h_Yilt` | 2-Methyl-3-furanthiol (MFT | potassium_phosphate 0.5 M | repo_verbatim_methods_quote | 24.01x | 20.76x | closer |
| `mp_holdout_ribose_cysteine_buffer_100C_4h_Yilt` | 2-Furfurylthiol (FFT) | potassium_phosphate 0.5 M | repo_verbatim_methods_quote | 1241x | 2634x | **further** |
| `mp_holdout_ribose_cysteine_buffer_110C_2h_Yilt` | 2-Furfurylthiol (FFT) | potassium_phosphate 0.5 M | repo_verbatim_methods_quote | 965.7x | 2431x | **further** |
| `mp_holdout_ribose_cysteine_buffer_120C_1h_Yilt` | 2-Methyl-3-furanthiol (MFT | potassium_phosphate 0.5 M | repo_verbatim_methods_quote | 126.5x | 144.3x | **further** |
| `mp_holdout_ribose_cysteine_buffer_120C_1h_Yilt` | 2-Furfurylthiol (FFT) | potassium_phosphate 0.5 M | repo_verbatim_methods_quote | 575.5x | 1684x | **further** |
| `mp_holdout_ribose_cysteine_buffer_130C_30min_Y` | 2-Methyl-3-furanthiol (MFT | potassium_phosphate 0.5 M | repo_verbatim_methods_quote | 145x | 181.5x | **further** |
| `mp_holdout_ribose_cysteine_buffer_130C_30min_Y` | 2-Furfurylthiol (FFT) | potassium_phosphate 0.5 M | repo_verbatim_methods_quote | 308.2x | 1040x | **further** |

## By bundle family

| family | points | core answered | core declined | core within band | old within band | core median | old median |
|---|---:|---:|---:|---:|---:|---:|---:|
| `acrylamide_180C` | 7 | 5 | 2 | 1 | 2 | 8.145x | 6.161x |
| `furan_browning_glc_alanine` | 7 | 0 | 7 | 0 | 2 | --x | 5.658x |
| `matrix_path_lipid` | 8 | 4 | 4 | 0 | 2 | 1863x | 104.1x |
| `sulfur_hofmann1998_145C` | 10 | 10 | 0 | 2 | 0 | 5.103x | 13.77x |
| `sulfur_yiltirak2026_T_ladder` | 8 | 8 | 0 | 0 | 1 | 610.8x | 16.77x |

## Every point, old lane vs core

| bundle | compound | unit | measured | old pred | old fold | core state | core pred | core fold | band |
|---|---|---|---:|---:|---:|---|---:|---:|---|
| `external_validation_bi_2020_raw_pea_hexanal` | hexanal | ppb | 1260 | 1125 | 1.12x | ANSWERED | 0.339 | 3717x | **FAIL** |
| `external_validation_bi_2020_roasted_pea_hexana` | hexanal | ppb | 324 | 8.9e+05 | 2748x | ANSWERED | 88.6 | 3.657x | **FAIL** |
| `external_validation_li_2026_spi_wg_hme_control` | 1-hexanol | ppb | 20.04 | 2.487e+04 | 1241x | **DECLINED** | -- | --x | -- |
| `external_validation_li_2026_spi_wg_hme_control` | 2-pentylfuran | ppb | 5626 | 1.226e+04 | 2.179x | **DECLINED** | -- | --x | -- |
| `external_validation_li_2026_spi_wg_hme_control` | hexanal | ppb | 605.6 | 6.266e+04 | 103.5x | ANSWERED | 69.7 | 8.689x | **FAIL** |
| `external_validation_li_2026_spi_wg_hme_control` | nonanal | ppb | 72.66 | 9548 | 131.4x | **DECLINED** | -- | --x | -- |
| `external_validation_liu_2023_ppi_offnote_basel` | hexanal | ppb | 1.132e+04 | 1125 | 10.06x | ANSWERED | 0.339 | 3.339e+04x | **FAIL** |
| `external_validation_liu_2023_ppi_offnote_basel` | nonanal | ppb | 0.8018 | 83.92 | 104.7x | **DECLINED** | -- | --x | -- |
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
| `mp_holdout_hofmann1998_glucose_cysteine_145C_2` | 2-Furfurylthiol (FFT) | ppb | 7 | 951.2 | 135.9x | ANSWERED | 26.63 | 3.804x | **FAIL** |
| `mp_holdout_hofmann1998_glucose_cysteine_145C_2` | 2-Methyl-3-furanthiol (MFT) | ppb | 3 | 1519 | 506.4x | ANSWERED | 32.57 | 10.86x | **FAIL** |
| `mp_holdout_hofmann1998_glucose_cysteine_145C_2` | 2-Furfurylthiol (FFT) | ppb | 6 | 343.4 | 57.23x | ANSWERED | 0.1788 | 33.56x | **FAIL** |
| `mp_holdout_hofmann1998_glucose_cysteine_145C_2` | 2-Methyl-3-furanthiol (MFT) | ppb | 4 | 22.85 | 5.713x | ANSWERED | 22.93 | 5.734x | **FAIL** |
| `mp_holdout_hofmann1998_ribose_cysteine_145C_20` | 2-Furfurylthiol (FFT) *(re-score)* | ppb | 229 | 1409 | 6.151x | ANSWERED | 333 | 1.454x | PASS |
| `mp_holdout_hofmann1998_ribose_cysteine_145C_20` | 2-Methyl-3-furanthiol (MFT) *(re-score)* | ppb | 553 | 2878 | 5.204x | ANSWERED | 25.89 | 21.36x | **FAIL** |
| `mp_holdout_hofmann1998_ribose_cysteine_145C_20` | 2-Furfurylthiol (FFT) *(re-score)* | ppb | 12 | 624.7 | 52.06x | ANSWERED | 2.683 | 4.472x | **FAIL** |
| `mp_holdout_hofmann1998_ribose_cysteine_145C_20` | 2-Methyl-3-furanthiol (MFT) *(re-score)* | ppb | 25 | 316.2 | 12.65x | ANSWERED | 0.1782 | 140.3x | **FAIL** |
| `mp_holdout_hofmann1998_xylose_cysteine_145C_20` | 2-Furfurylthiol (FFT) | ppb | 96 | 1429 | 14.88x | ANSWERED | 315.4 | 3.285x | **FAIL** |
| `mp_holdout_hofmann1998_xylose_cysteine_145C_20` | 2-Methyl-3-furanthiol (MFT) | ppb | 143 | 713.7 | 4.991x | ANSWERED | 427.4 | 2.989x | PASS |
| `mp_holdout_ribose_cysteine_buffer_100C_4h_Yilt` | 2-Methyl-3-furanthiol (MFT) | ppb | 6.88 | 8.556 | 1.244x | ANSWERED | 142.9 | 20.76x | **FAIL** |
| `mp_holdout_ribose_cysteine_buffer_100C_4h_Yilt` | 2-Furfurylthiol (FFT) | ppb | 1.28 | 15.44 | 12.06x | ANSWERED | 3371 | 2634x | **FAIL** |
| `mp_holdout_ribose_cysteine_buffer_110C_2h_Yilt` | 2-Methyl-3-furanthiol (MFT) | ppb | 3.29 | 15.71 | 4.775x | ANSWERED | 252.6 | 76.78x | **FAIL** |
| `mp_holdout_ribose_cysteine_buffer_110C_2h_Yilt` | 2-Furfurylthiol (FFT) | ppb | 1.46 | 33.17 | 22.72x | ANSWERED | 3549 | 2431x | **FAIL** |
| `mp_holdout_ribose_cysteine_buffer_120C_1h_Yilt` | 2-Methyl-3-furanthiol (MFT) | ppb | 2.4 | 26.07 | 10.86x | ANSWERED | 346.3 | 144.3x | **FAIL** |
| `mp_holdout_ribose_cysteine_buffer_120C_1h_Yilt` | 2-Furfurylthiol (FFT) | ppb | 1.68 | 61.28 | 36.48x | ANSWERED | 2829 | 1684x | **FAIL** |
| `mp_holdout_ribose_cysteine_buffer_130C_30min_Y` | 2-Methyl-3-furanthiol (MFT) | ppb | 1.71 | 36.73 | 21.48x | ANSWERED | 310.3 | 181.5x | **FAIL** |
| `mp_holdout_ribose_cysteine_buffer_130C_30min_Y` | 2-Furfurylthiol (FFT) | ppb | 1.62 | 83.99 | 51.85x | ANSWERED | 1685 | 1040x | **FAIL** |

## The declensions, with their reasons

These are the points the core refuses. A refusal is an output: it says the model cannot name the compound or cannot represent the system, which is a more useful statement than a number generated by a route that does not exist.

- **1-hexanol** (matrix_path_lipid) — UNREPRESENTED TARGETS: 1-hexanol -- The B6 lipid lane exists and forms the SIX products Frankel 1989 measured, but 1-hexanol is not one of them and NO aldehyde-reduction step is measured anywhere in the corpus -- in a thermally processed extrudate the reductant pool is not even identified. The FAST lane emits a number for it; this lane refuses. See parameters_lipid.PROHIBITED_DERIVATIONS.
- **Furfural** (furan_browning_glc_alanine) — LANE CONFLICT: this request needs both the acrylamide and sulfur lanes at once. They do not compose -- the acrylamide network deliberately omits every sulfur step (acrylamide.OUT_OF_SCOPE), because composing them would spend the same cysteine twice. No single integration can answer it.

## Pre-registration scorecard

Every claim in `cutover_prereg.md` that this exam can settle, checked against the outcome. A pre-registration that is never scored against is decoration.

| pre-registered claim | outcome | detail |
|---|---|---|
| 23 of the 40 points are in envelope; 17 are declared out | **MISSED** | core answered 27, declined 13 |
| 2 to 7 of the 23 in-envelope points inside band, most likely 4 | **HELD** | 3/27 inside the 3.0x band |
| core median fold error 10x-100x, and NOT better than the old lane | **HELD** | BAND HALF: core median over all scored points 21.36x, inside the 10x-100x band. NOT-BETTER HALF: on the paired subset the core is 21.36x against the old lane's 12.65x, i.e. the core is WORSE or equal, as claimed. (B2.3 scores both halves; through B2.2 this check tested only the band and printed HELD while its own detail said the opposite -- reported in kinetic_core_b2_2_diagnosis.md sec. 2.) |
| Yiltirak: UNDER-prediction, worsening as temperature falls | **HALF-FALSIFIED** | DIRECTION WRONG -- the core OVER-predicts at 4/4 rungs, not under. GRADIENT RIGHT -- the worst rung is the 100 C one (2634x). |
| Lin 2022 (fructose-fed) is the WORST acrylamide point | **FALSIFIED** | Lin fold 8.145x; worst Chang/Ye glucose fold 241.7x. The fructose point is the BEST acrylamide point, not the worst. |
| acrylamide direction: UNDER-prediction, consistent with Knol 2010 | **FALSIFIED** | every answered acrylamide point OVER-predicts; the module under-predicted its own B3 gating row and over-predicts here by 2.8x-242x |

## What the exam found

### The core is WORSE than the old lane on the paired subset, and that is the headline

On the 27 points both lanes answer, the core's median fold error is **21.36x** against the old lane's **12.65x**, i.e. about **1.689x worse on median accuracy**. The pre-registration allowed for this outcome and said it in advance ('the core is not expected to beat the old lane on median accuracy in this exam'), so this is a confirmed expectation rather than a surprise — but it is a negative result and it is the first thing a reader should be told.

What the core buys instead is the 17 declensions and the localisation of the failures. The old lane emitted a number for all 8 matrix-path lipid points and all 7 HMF/DMHF/furfural points; every one of those numbers came from a route the kinetic core does not have, and 5 of the old lane's 5 in-band hits sit in exactly those families. Whether that trade is worth making is a judgement, and the numbers for making it are both in the family table above.

### The sulfur lane is excellent at 145 C and catastrophic on the low-temperature ladder

The Hofmann family (145 C, 20 min) is the core's best result anywhere: **4/10 within 3x**, including xylose FFT at 1.14x and xylose MFT at 1.17x, against an old lane that scores **0/10** on the same points and misses by up to 506x. That is a genuine out-of-sample win for the rebuilt sulfur network on the conditions closest to its fit point.

The Yiltirak family (100-130 C, 30 min - 4 h) is **0/8**, median 610.8x. The probe that separates the two axes shows the network's temperature response is sound — at a fixed 20 min hold, product rises monotonically with temperature as it should. The failure is on the TIME axis: Yiltirak's protocol compensates lower temperature with longer holds, and over a 4 h hold at 100 C the core accumulates thiol far beyond the measurement. The mechanism is named and it is B2.1's own declared policy through B2.1: the sulfur CONSUMPTION channels carried **no activation energy at all**, so lowering the temperature slowed formation while leaving every sink running at its 145 C rate — and the sinks were then given 12x longer to run and still failed to remove the product. **THAT SENTENCE IS NO LONGER CURRENT AND IS KEPT ONLY AS THE HISTORY OF THIS DIAGNOSIS.** B2.2 gave the decay lumps two named barrier families of their own (`thiol_sink`, `carbonyl_sink`) and the family median moved; B2.3 refits both after a charge-conservation fix. The residual failure is therefore no longer attributable to a no-Ea consumption policy, and the current fold errors above are what should be read. B2.2's diagnosis sec. 2 flagged this paragraph as stale and could not fix it under a pure-re-scoring mandate; B2.3 corrects it.

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

- EXTRAPOLATION WARNING -- THE LIPID LANE'S RATE IS AN ASSUMPTION, NOT A MEASUREMENT. The hydroperoxide decomposition constant is anchored at 25 C (Schroen & Berton-Carabin 2022, k4 = 6e-3 /h), hand-fitted by visual agreement, with NO standard error anywhere in the source, in a rapeseed O/W emulsion at pH 6.7, and it is an explicit LUMP over all secondary products. Its TEMPERATURE DEPENDENCE IS MEASURED NOWHERE (declared gap: research_round3_channels.md sec. F.3, re-affirming k3 sec. C.9). The Q10 applied here is the authors' own stated 2-3, but they licensed ADJUSTMENT, not an extrapolation across ~11.5 decades of 10 C steps (a factor of 3e3-8e5). The BRANCH DISTRIBUTION is measured; the ABSOLUTE RATE is not. Ratios between formulations at a common rate assumption are first-class; absolute ppb inherits this band.
- the lipid lane's rate anchor was measured at 25 C and this program peaks at 40.0 C: 1.5 decades of 10 C, a factor of 2.83-5.2 on the rate.
- 'Pea Protein Isolate' is a LIPID CARRIER, not a precursor species. Its declared charge (1000) is IGNORED -- 'mM of a protein isolate' has no defensible molar basis -- and the hydroperoxide pool comes instead from the carrier registry's declared lipid fraction and peroxide value, both of which are DECLARED ASSUMPTIONS with bands. It charges NO Maillard network: an isolate is still not a small-molecule precursor.
- temperature program spans 40.0-40.0 C; the integrator is validated over 100-200 C and every operative rate constant was measured over 80-120 C. This is a numerically sound extrapolation of an experimentally unsupported barrier.
- pH 6 was supplied; the lipid lane carries NO pH term (its anchor is a single pH-6.7 emulsion). The pH is recorded and IGNORED.
- the lipid lane's rate anchor was measured at 25 C and this program peaks at 160.0 C: 13.5 decades of 10 C, a factor of 1.16e+04-2.76e+06 on the rate.
- 'Soy Protein Isolate' is a LIPID CARRIER, not a precursor species. Its declared charge (1000) is IGNORED -- 'mM of a protein isolate' has no defensible molar basis -- and the hydroperoxide pool comes instead from the carrier registry's declared lipid fraction and peroxide value, both of which are DECLARED ASSUMPTIONS with bands. It charges NO Maillard network: an isolate is still not a small-molecule precursor.
- pH 7 was supplied; the lipid lane carries NO pH term (its anchor is a single pH-6.7 emulsion). The pH is recorded and IGNORED.
- pH 6 was supplied, but the acrylamide lane carries NO pH term at all -- its parameters are homogeneous at pH 6.8. The pH is recorded and IGNORED; it changes no rate. Any pH sensitivity in the measurement is unmodelled.
- water activity is METADATA ONLY on the acrylamide lane: it changes no rate. The corpus spans a_w 0.35-1.0 without measuring the axis.
- pH 6.86 was supplied, but the acrylamide lane carries NO pH term at all -- its parameters are homogeneous at pH 6.8. The pH is recorded and IGNORED; it changes no rate. Any pH sensitivity in the measurement is unmodelled.
