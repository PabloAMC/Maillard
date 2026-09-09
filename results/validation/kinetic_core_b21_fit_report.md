# B21 -- the aqueous glucosone route to glyoxal

*Generated 2026-09-09 by `scripts/generators/generate_kinetic_core_b21_fit.py`; pre-registration `results/validation/kinetic_core_b21_prereg.md`; docs/reference/FIT_HOLDOUT_DECLARATION.md Amendment 31.*

Objective: six printed first-order constants (Hamzalioglu 2026, 110-140 C), log10 residuals; two log10 constants at 100 C free, barriers declared. Cost 3.543 on 6 rows, reduced chi-square 0.89.

## The fitted constants (log10 at 100 C; barriers declared)

| coordinate | value | sigma (dex) | identified |
|---|---|---|---|
| log10_k_ama_g_100C | -2.0339 | 0.085 | True |
| log10_k_g_go_aqueous_100C | -0.5104 | 0.144 | True |

The glass value of k_g_go at 100 C was 0.00354 per minute (Ea 93.8); the aqueous value is 0.309.

## The rows

| row | printed | model | residual (dex) | decisive |
|---|---|---|---|---|
| ham_k_ama_g_110C | 0.019 | 0.0175 | -0.04 | True |
| ham_k_ama_g_120C | 0.021 | 0.0321 | +0.18 | True |
| ham_k_ama_g_130C | 0.025 | 0.0571 | +0.36 | False |
| ham_k_ama_g_140C | 0.15 | 0.0988 | -0.18 | True |
| ham_k_g_go_120C | 0.33 | 0.331 | +0.00 | True |
| ham_k_g_go_130C | 0.35 | 0.341 | -0.01 | False |

## Diagnostics, before glass

- Quan 2020 glyoxal at 100C, 21 min: model 0.0000 mmol/L, printed [0.052, 0.127] (-4.21 dex from the range)
- Quan 2020 glyoxal at 130C, 21 min: model 0.0004 mmol/L, printed [0.144, 0.605] (-2.53 dex from the range)
- Xia 2022 at 130 C, 80 min: glyoxal 0.001 vs methylglyoxal 8.138 mmol/L; glyoxal above: False
- Leahy 1989 total pyrazine, 95 C, 2 h: model 17.4 ug/L vs 13100 (-2.88 dex; B18 recorded -2.88)
- B1 browning hold-out: median fold 1.43, max 2.85, within 2x 0.95, within 3x 1.00

## Diagnostics, after candidate

- Quan 2020 glyoxal at 100C, 21 min: model 0.0234 mmol/L, printed [0.052, 0.127] (-0.35 dex from the range)
- Quan 2020 glyoxal at 130C, 21 min: model 0.3871 mmol/L, printed [0.144, 0.605] (+0.00 dex from the range)
- Xia 2022 at 130 C, 80 min: glyoxal 1.227 vs methylglyoxal 7.823 mmol/L; glyoxal above: False
- Leahy 1989 total pyrazine, 95 C, 2 h: model 21.9 ug/L vs 13100 (-2.78 dex; B18 recorded -2.88)
- B1 browning hold-out: median fold 1.31, max 2.79, within 2x 0.95, within 3x 1.00
