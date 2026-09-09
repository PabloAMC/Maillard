# B20 -- the glycation arm: protein-bound lysine as a reactant on the trunk lane

*Generated 2026-09-09 by `scripts/generators/generate_kinetic_core_b20_fit.py`; pre-registration `results/validation/kinetic_core_b20_prereg.md`; docs/reference/FIT_HOLDOUT_DECLARATION.md Amendment 30.*

Objective: ten printed rate constants (Nguyen 2016 M1, 120 / 130 C), log10 residuals over the printed interval; five log10 constants at 100 C free, five barriers declared. Cost 17.206 on 10 rows, reduced chi-square 3.44; best start 1.

## The fitted constants (log10 at 100 C; barriers declared)

| coordinate | value | sigma (dex) | identified | band |
|---|---|---|---|---|
| log10_k_glyc_100C | -4.7073 | 0.085 | True | [-6.51, -2.51] |
| log10_k_flp_cml_100C | -3.3174 | 0.175 | True | [-4.86, -0.86] |
| log10_k_flp_cel_100C | -3.6521 | 0.112 | True | [-5.29, -1.29] |
| log10_k_flp_decay_100C | -1.8560 | 0.152 | True | [-3.91, 0.09] |
| log10_k_cml_loss_100C | -0.9827 | 0.253 | True | [-2.83, 1.17] |

## The rows

| row | printed | model | residual (dex) | decisive |
|---|---|---|---|---|
| nguyen_k_glyc_120C | 0.00015 +/- 3e-05 | 9.59e-05 | -0.19 | True |
| nguyen_k_glyc_130C | 0.00016 +/- 2.2e-05 | 0.0002 | +0.10 | True |
| nguyen_k_flp_cml_120C | 0.0088 +/- 0.0066 | 0.00307 | -0.46 | True |
| nguyen_k_flp_cml_130C | 0.006 +/- 0.0016 | 0.00724 | +0.08 | True |
| nguyen_k_flp_cel_120C | 0.0023 +/- 0.0027 | 0.00101 | -0.36 | False |
| nguyen_k_flp_cel_130C | 0.002 +/- 0.0003 | 0.00202 | +0.01 | True |
| nguyen_k_flp_decay_120C | 0.052 +/- 0.033 | 0.0685 | +0.12 | True |
| nguyen_k_flp_decay_130C | 0.15 +/- 0.034 | 0.143 | -0.02 | True |
| nguyen_k_cml_loss_120C | 0.29 +/- 0.27 | 0.104 | -0.45 | True |
| nguyen_k_cml_loss_130C | 0.077 +/- 0.033 | 0.104 | +0.13 | True |

## Nguyen's pot, integrated (LYSP 16, glucose 150 mmol/L, pH 6.8, 30 min), mmol/L

- 120C: CML 0.0529 (printed range [0.025, 0.135]), CEL 0.0462, fructosyl-lysine 2.061, lysine lost 14.1 %
- 130C: CML 0.1230 (printed range [0.025, 0.135]), CEL 0.1038, fructosyl-lysine 1.798, lysine lost 14.2 %

- Berk 2021, fructosyl-lysine -> CML at 180 C (dry sesame): model 0.299 vs 0.00554 per minute (+1.73 dex)
- Hamzalioglu 2026, lactulosyl-lysine -> CML at 110C (milk): model 0.00125 vs 0.00017 (+0.87 dex)
- Hamzalioglu 2026, lactulosyl-lysine -> CML at 120C (milk): model 0.00307 vs 0.0036 (-0.07 dex)
- Hamzalioglu 2026, lactulosyl-lysine -> CML at 130C (milk): model 0.00724 vs 0.0014 (+0.71 dex)
- Hamzalioglu 2026, lactulosyl-lysine -> CML at 140C (milk): model 0.0164 vs 0.0011 (+1.17 dex)
- expanded soybean, 110 C, 60 min: about 25 % of the bound lysine lost (a moist solid with sucrose, which the engine cannot charge; direction only)
