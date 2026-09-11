# Wave B41 ship rule: SHIP

_SHIP if P1, P2 (Leitzen: 3,4-DGE < 10x, 3-DG <= 3x, HMF within 3 %, METHYLGLYOXAL within 3 % of 1.28x) and P3 (no row leaves the band); P4 reported (prereg sec. 5). Judged on frozen artifacts: the fit report and the panel before/after._

- **P1** rows fit: **HELD** — maxima within 0.3 dex True, t_max in brackets True, χ²_red 0.12
- **P2** k_tdg_ddg moves less than in B39: HELD ({'b39_delta_dex': 0.5656201217641403, 'b41_delta_dex': 0.49508700356741375, 'held': True}) — reported
- **P2 decisive rows**: 3-DG ≤ 3× **HELD** (1.1110425180991021, 1.2712448571778727); methylglyoxal within 3 % **HELD** (1.2809566500613123, 1.010171781531992)
- **P4** the Leitzen hold-out, never read by the fit: **HELD**

| Leitzen 2021 row | fold error before | after |
|---|---:|---:|
| 3,4-dideoxyglucosone | 32.43x | 6.63x |
| 3-deoxyglucosone | 1.11x | 1.27x |
| 5-Hydroxymethylfurfural (HMF) | 11.93x | 9.31x |
| methylglyoxal | 1.28x | 1.01x |

- **P5** ≥ 3 pinned: HELD (4 pinned; collinear []) — reported
- **P6** no row within 3× leaves the band: **HELD** — left []; entered [('mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019', '5-Hydroxymethylfurfural (HMF)'), ('mp_holdout_glucose_asparagine_180C_30min_water_Chang2021', '5-Hydroxymethylfurfural (HMF)')]; within 10/45 → 12/45

## Every other row that moved by more than 1 %

| benchmark | compound | before | after |
|---|---|---:|---:|
| mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019 | 5-Hydroxymethylfurfural (HMF) | 2.05x | 1.03x |
| mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019 | 5-Hydroxymethylfurfural (HMF) | 3.60x | 1.82x |
| mp_holdout_glucose_asparagine_180C_30min_water_Chang2021 | 5-Hydroxymethylfurfural (HMF) | 3.19x | 2.74x |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 3,4-dideoxyglucosone | 32.43x | 6.63x |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 3-deoxyglucosone | 1.11x | 1.27x |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | 5-Hydroxymethylfurfural (HMF) | 11.93x | 9.31x |
| mp_holdout_glucose_only_autoclave_121C_Steinhagen2021 | methylglyoxal | 1.28x | 1.01x |
