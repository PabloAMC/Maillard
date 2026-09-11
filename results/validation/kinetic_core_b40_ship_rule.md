# Wave B40 ship rule: DO NOT SHIP

_SHIP if P1 and P4 and P6 and Leitzen 3-DG <= 3x (P3 decisive half); P2, P5 reported (prereg sec. 5). Judged on frozen artifacts: the fit report and the panel before/after._

- **P1** rows fit: **HELD** — maxima within 0.3 dex True, t_max in brackets True, χ²_red 0.15
- **P2** k_tdg_ddg moves less than in B39: HELD ({'b39_delta_dex': 0.5656201217641403, 'b40_delta_dex': 0.5342132713252836, 'held': True}) — reported
- **P3** Leitzen 3-DG at or under 3× (decisive): **HELD** — (1.1110425180991021, 1.3618958740086156)
- **P4** the Leitzen hold-out, never read by the fit: **HELD**

| Leitzen 2021 row | fold error before | after |
|---|---:|---:|
| 3,4-dideoxyglucosone | 32.43x | 6.01x |
| 3-deoxyglucosone | 1.11x | 1.36x |
| 5-Hydroxymethylfurfural (HMF) | 11.93x | 7.19x |

- **P5** ≥ 3 pinned: REFUTED (2 pinned; collinear [{'a': 'log10_k_tdg_ddg_100C', 'b': 'log10_k_dgal_ddg_100C', 'corr': 0.9551700902642359}, {'a': 'log10_k_dgal_ddg_100C', 'b': 'log10_k_ddg_hmf_100C', 'corr': 0.9614178830956657}]) — reported
- **P6** no row within 3× leaves the band: **REFUTED** — left [('mp_holdout_glucose_only_autoclave_121C_Steinhagen2021', 'methylglyoxal')]; entered [('mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019', '5-Hydroxymethylfurfural (HMF)'), ('mp_holdout_glucose_asparagine_180C_30min_water_Chang2021', '5-Hydroxymethylfurfural (HMF)')]; within 10/45 → 11/45

## Every other row that moved by more than 1 %

| benchmark | compound | before | after |
|---|---|---:|---:|
| mp_holdout_fructose_asparagine_180C_Lin2022 | 5-Hydroxymethylfurfural (HMF) | 6.31x | 6.22x |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | 5-Hydroxymethylfurfural (HMF) | 2.05x | 1.07x |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | 5-Hydroxymethylfurfural (HMF) | 3.60x | 1.64x |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | 5-Hydroxymethylfurfural (HMF) | 3.19x | 2.59x |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 3,4-dideoxyglucosone | 32.43x | 6.01x |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 3-deoxyglucosone | 1.11x | 1.36x |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 5-Hydroxymethylfurfural (HMF) | 11.93x | 7.19x |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | methylglyoxal | 1.28x | 33.18x |
