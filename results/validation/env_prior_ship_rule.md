# ENV-B18 and ENV-B13 ship rule: INSTALL

*Rule: INSTALL if T1 (the prior rows exist), T2 (every row the priors reach widens, and no other row moves beyond the MEASURED noise floor) and T3 (no median moves) hold; T4 reported.*

| test | result | pass |
|---|---|---|
| T1 prior rows | B18 6 rows, 4 sampled; B13 8 rows, 4 sampled | True |
| T2 widths | 5 rows the priors reach; 34 they do not; violations 0 | True |
| T3 medians | floor 0.121 dex from 44 seed pairs; 0 unreached medians beyond it | True |
| T4 | 5 widened; newly inside 0 | reported |

## The measured Monte-Carlo noise floor

Two runs of the SAME priors at different seeds, 44 rows compared. Relative difference in interval width: median 3.22%, 95th percentile 13.37%, worst 16.99%.

This is why the pre-registrations' "not one row may narrow" could not be tested as written. The sampler draws every coordinate from one stream, so adding coordinates re-shuffles every later draw. Measuring the floor turns an untestable rule into a testable one.

## What actually widened

| compound | benchmark | before | after | change |
|---|---|---:|---:|---:|
| 5-Hydroxymethylfurfural (HMF) | mp_holdout_glucose_alanine_130C_2h_pH50_Schi | 0.000 | 0.846 | +0.846 dex |
| 5-Hydroxymethylfurfural (HMF) | mp_holdout_glucose_alanine_130C_2h_pH80_Schi | 0.000 | 0.846 | +0.846 dex |
| 5-Hydroxymethylfurfural (HMF) | mp_holdout_glucose_asparagine_180C_30min_wat | 0.010 | 0.380 | +0.370 dex |
| 5-Hydroxymethylfurfural (HMF) | mp_holdout_fructose_asparagine_180C_Lin2022 | 0.001 | 0.355 | +0.355 dex |
| 5-Hydroxymethylfurfural (HMF) | mp_holdout_glucose_only_autoclave_121C_Stein | 0.000 | 0.135 | +0.135 dex |

Widths are the 90 % interval in decades. A row at 0.000 before was being published with NO interval at all: the model was asserting it exactly.
