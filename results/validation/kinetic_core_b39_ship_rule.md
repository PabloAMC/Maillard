# Wave B39 ship rule: DO NOT SHIP

_SHIP if P1 and P4 and P6; P2, P3, P5 reported (prereg sec. 6). Judged on frozen artifacts: the fit report and the panel before/after._

- **P1** rows fit: **REFUTED** — maxima within 0.3 dex True, t_max in brackets False, χ²_red 0.58
- **P2** k_tdg_ddg up 0.5–1.2 dex: HELD (Δ +0.57 dex) — reported
- **P3** k_ddg_hmf down ≥ 0.5 dex: HELD (Δ -0.60 dex) — reported
- **P4** the Leitzen hold-out, never read by the fit: **REFUTED**

| Leitzen 2021 row | fold error before | after |
|---|---:|---:|
| 3,4-dideoxyglucosone | 32.43x | 8.09x |
| 3-deoxyglucosone | 1.11x | 1.23x |
| 5-Hydroxymethylfurfural (HMF) | 11.93x | 12.21x |

- **P5** ≥ 3 pinned: HELD (5 pinned; collinear []) — reported
- **P6** no row within 3× leaves the band: **HELD** — left []; entered [('mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019', '5-Hydroxymethylfurfural (HMF)'), ('mp_holdout_glucose_asparagine_180C_30min_water_Chang2021', '5-Hydroxymethylfurfural (HMF)')]; within 10/45 → 12/45

## Every other row that moved by more than 1 %

| benchmark | compound | before | after |
|---|---|---:|---:|
| mp_holdout_fructose_asparagine_180C_Lin2022 | 5-Hydroxymethylfurfural (HMF) | 6.31x | 6.24x |
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | 5-Hydroxymethylfurfural (HMF) | 2.05x | 1.03x |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | 5-Hydroxymethylfurfural (HMF) | 3.60x | 1.70x |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | 5-Hydroxymethylfurfural (HMF) | 3.19x | 2.72x |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 3,4-dideoxyglucosone | 32.43x | 8.09x |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 3-deoxyglucosone | 1.11x | 1.23x |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 5-Hydroxymethylfurfural (HMF) | 11.93x | 12.21x |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | methylglyoxal | 1.28x | 1.38x |
