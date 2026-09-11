# B22 -- the methionine chain on the sugar path

*Generated 2026-09-09 by `scripts/generators/generate_kinetic_core_b22_fit.py`; pre-registration `results/validation/kinetic_core_b22_prereg.md`; docs/reference/FIT_HOLDOUT_DECLARATION.md Amendment 32.*

Objective: nine zero-order rates (Pan 2025, 100 / 120 / 140 C) as mean formation rates over 30-600 s in Pan's pot; four free. Cost 7373.94 on 9 rows, reduced chi-square 1474.8.

## The fitted coordinates

| coordinate | value | sigma | identified | on bound | band |
|---|---|---|---|---|---|
| log10_identity_ratio_met_over_gly | 2.000 | 2.5 | False | True | [-2.0, 2.0] |
| log10_k_mtal_msh_100C | -0.281 | 6.16 | False | False | [-5.0, 1.0] |
| ea_mtal_msh_kj_mol | 20.000 | 643 | False | True | [20.0, 150.0] |
| log10_k_msh_dmds_100C | 2.000 | 8.77 | False | True | [-4.0, 2.0] |

## The rows (mean formation rate, umol L-1 min-1)

| row | printed | model | residual (dex) | decisive |
|---|---|---|---|---|
| pan_MTAL_rate_100C | 0.0109 | 2.97e-08 | -5.57 | True |
| pan_MTAL_rate_120C | 0.101 | 3e-06 | -4.53 | True |
| pan_MTAL_rate_140C | 0.539 | 0.000124 | -3.64 | True |
| pan_MSH_rate_100C | 0.00947 | 4.52e-08 | -5.32 | True |
| pan_MSH_rate_120C | 0.0141 | 6e-06 | -3.37 | True |
| pan_MSH_rate_140C | 0.038 | 0.000393 | -1.99 | True |
| pan_DMDS_rate_100C | 9e-05 | 2.45e-15 | -10.56 | False |
| pan_DMDS_rate_120C | 0.00042 | 1.54e-10 | -6.43 | False |
| pan_DMDS_rate_140C | 0.00112 | 2.54e-06 | -2.64 | False |

## Diagnostics

- Deng 2022, methional at 120 C, 30min: model 1.58e+04 ug/L vs printed 18.52 (+2.93 dex)
- Deng 2022, methional at 120 C, 60min: model 1.32e+04 ug/L vs printed 25.27 (+2.72 dex)
- Deng 2022, methional at 120 C, 120min: model 6.45e+03 ug/L vs printed 75.76 (+1.93 dex)
- Deng 2022, methional at 120 C, 180min: model 3.11e+03 ug/L vs printed 87.77 (+1.55 dex)
- Deng 2022, rising 30 -> 120 min: False
- Chin & Lindsay 1994, methanethiol half-life at 30 C: model 69.2 min vs 17.0 min with copper
- Pan's pot at 140 C / 10 min, mmol/L: {"GO": 0.00944, "MGO": 0.0424, "MTAL": 1.18e-06, "MSH": 3.74e-06, "DMDS": 2.42e-08, "MET": 0.268}
