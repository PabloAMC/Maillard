# ENV-B34 ship rule: INSTALL

*Rule: INSTALL if T1 (the prior rows exist), T2 (every reached row widens; no other row moves beyond the MEASURED noise floor) and T3 (no unreached median moves) hold; T4 reported. Pre-registration `results/validation/kinetic_core_env_b34_prereg.md`.*

| test | result | pass |
|---|---|---|
| T1 prior rows | 8 rows, 8 sampled | True |
| T2 widths | 9 named rows must widen; 11 other Maillard-lane rows moved (reported); lipid rows bit-identical 8/8; violations 0 | True |
| T3 medians | lipid (unreached) medians moved: 0 (must be 0, exact under ENV-M1); Maillard-lane medians moved through shared glucose: 11, largest 0.157 dex; the old two-seed floor for comparison 0.069 dex | True |
| T4 | 9 widened; newly inside 4, newly outside 0; coverage [12, 42, 3] -> [16, 44, 1]; median width 0.9148 -> 0.9874 dex | reported |

## The measured Monte-Carlo noise floor

Two runs of the SAME priors at different seeds, 45 rows compared. Relative difference in interval width: median 5.63%, worst 20.98%. The floor is the observed MAXIMUM, not a quantile.

## The rows the priors reach

| compound | benchmark | width before | width after | delta |
|---|---|---:|---:|---:|
| 5-Hydroxymethylfurfural (HMF) | mp_holdout_fructose_asparagine_180C_Lin2022 | 0.247 | 0.383 | +0.136 dex |
| DMHF | mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | 0.883 | 0.960 | +0.077 dex |
| 5-Hydroxymethylfurfural (HMF) | mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | 0.662 | 1.122 | +0.460 dex |
| DMHF | mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | 0.883 | 0.960 | +0.077 dex |
| 5-Hydroxymethylfurfural (HMF) | mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | 0.662 | 1.122 | +0.460 dex |
| 5-Hydroxymethylfurfural (HMF) | mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | 0.268 | 0.824 | +0.556 dex |
| 5-Hydroxymethylfurfural (HMF) | mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 0.091 | 0.969 | +0.878 dex |
| 3-deoxyglucosone | mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 0.000 | 0.991 | +0.991 dex |
| 3,4-dideoxyglucosone | mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 0.000 | 0.984 | +0.984 dex |

Newly inside their interval: 5-Hydroxymethylfurfural (HMF) in `mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019`, 5-Hydroxymethylfurfural (HMF) in `mp_holdout_glucose_asparagine_180C_30min_water_Chang2021`, 3-deoxyglucosone in `mp_holdout_glucose_only_autoclave_121C_Steinhagen2021`, methylglyoxal in `mp_holdout_glucose_only_autoclave_121C_Steinhagen2021`
