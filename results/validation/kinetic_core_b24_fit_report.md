# B24 -- 2-acetyl-1-pyrroline from proline

*Generated 2026-09-09 by `scripts/generators/generate_kinetic_core_b24_fit.py`; pre-registration `results/validation/kinetic_core_b24_prereg.md`; docs/reference/FIT_HOLDOUT_DECLARATION.md Amendment 33.*

Objective: five 30-minute yields (Hofmann 1998b Tables 7 and 9, 100 C, pH 7), log10 residuals; two log10 constants at 100 C free, barriers declared. Cost 53.693 on 5 rows, reduced chi-square 17.90.

## The fitted constants (log10 at 100 C, L/(mmol min); barriers declared)

| coordinate | value | sigma (dex) | identified |
|---|---|---|---|
| log10_k_pyrl_ap_100C | -2.5618 | 0.532 | True |
| log10_k_mgo_pro_100C | -6.6208 | 0.664 | True |

## The rows (yield at 30 min, mol % of the basis species)

| row | printed | model | residual (dex) |
|---|---|---|---|
| hof_t7_e1_pyrl2_mgo10 | 28.7 | 43.02 | +0.18 |
| hof_t7_e2_pyrl2_mgo2 | 5.3 | 10.54 | +0.30 |
| hof_t9_pro400_mgo4 | 0.0058 | 0.0002991 | -1.29 |
| hof_t9_pro400_mgo40 | 0.0125 | 0.01678 | +0.13 |
| hof_t9_pro400_mgo400 | 0.0179 | 0.2594 | +1.16 |

## Diagnostics

- Hofmann experiment 3 (1-pyrroline 10 + methylglyoxal 2 mmol/L): model 41.2 mol % of the methylglyoxal vs printed 0.33 (+2.10 dex): the source's suppression by excess pyrroline is not written
- apparent barrier of the whole cascade, glucose 100 + proline 100 mmol/L, pH 7, 75-115 C: model 591 kJ/mol vs Chan & Reineccius 1994's 60.2
- the B18 pot at 120 C / 60 min, mmol/L: without proline {"PZ": 0.00394, "MPZ": 0.00131, "DMP": 0.00044, "AP": 0.0}; with 10 mmol/L proline {"PZ": 0.00472, "MPZ": 0.00157, "DMP": 0.000525, "AP": 0.000752}
