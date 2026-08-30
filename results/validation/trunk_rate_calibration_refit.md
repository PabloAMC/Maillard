# Trunk rate calibration — the first rate-level fit of the sugar trunk

Wave S3 · generated 2026-08-27 · git `b2a1b20` (audit-remediation)

## (a) Lane

Fitted in a **dedicated integrator**, `src/trunk_kinetics.py`, not in either shipped lane. The FAST screening lane performs no time integration at all; the Cantera export lane integrates but discards the initial molarity (`src/kinetics.py:373` sets `phase.X`, which Cantera normalises), lacks four of the ten measured species, and lumps three separately-measured enolisation steps onto one A/Ea pair. Both blockers were measured, not assumed — see the module docstring. **No shipped prediction imports the fitted module.**

## (b) Fitted parameters

176 values, 16 free parameters, reduced chi-square 14.692. 1/60 random starts reached the reported optimum.

The sigmas are the REPLICATE scatter the data report about themselves. A reduced chi-square above 1 therefore measures how far the model's systematic error exceeds that scatter: here residuals run 3.8x the replicate sigma. Every confidence interval below has been widened by that same factor; unscaled J^T J intervals would have been 3.8x too tight.

| step | repo family | k(100 °C) | 95% CI | unit | rate identifiability |
|---|---|---:|---|---|---|
| `k_schiff` | schiff_condensation | 1.71e-05 | 1.48e-05 – 1.97e-05 | L/(mmol*min) | well constrained (95% CI inside +/-26%) |
| `k_amadori` | amadori_rearrangement | 0.1536 | 0.104 – 0.226 | 1/min | constrained (95% CI inside a factor of 2) |
| `k_dfg_3dg` | enolisation_intermediate | 0.006462 | 0.00396 – 0.0105 | 1/min | constrained (95% CI inside a factor of 2) |
| `k_dfg_1dg` | 2,3-enolisation | 0.03048 | 0.00558 – 0.167 | 1/min | weakly constrained (95% CI up to a factor of 10) |
| `k_dfg_other` | *(none — structural gap)* | 0.004312 | 2.95e-08 – 630 | 1/min | SLOPPY (order-of-magnitude only) |
| `k_3dg_out` | 1,2-enolisation | 0.1047 | 0.0568 – 0.193 | 1/min | constrained (95% CI inside a factor of 2) |
| `k_1dg_out` | furanone_amino_acid_reduction | 3.031 | 0.534 – 17.2 | 1/min | weakly constrained (95% CI up to a factor of 10) |
| `k_glc_other` | *(none — structural gap)* | 0.000271 | 4.31e-05 – 0.0017 | 1/min | weakly constrained (95% CI up to a factor of 10) |

| step | Ea kJ/mol | 95% CI | Ea identifiability |
|---|---:|---|---|
| `k_schiff` | 94.4 | 86 – 103 | well constrained |
| `k_amadori` | 34.7 | 13 – 57 | constrained |
| `k_dfg_3dg` | 144.8 | 116 – 174 | constrained |
| `k_dfg_1dg` | 73.0 | 15 – 131 | weakly constrained |
| `k_dfg_other` | 89.5 | -314 – 493 | UNCONSTRAINED SIGN (95% CI reaches -314 kJ/mol; an Ea <= 0 is unphysical, so the data do not determine this step's T-dependence) |
| `k_3dg_out` | 135.6 | 96 – 176 | weakly constrained |
| `k_1dg_out` | 40.1 | -21 – 101 | UNCONSTRAINED SIGN (95% CI reaches -21 kJ/mol; an Ea <= 0 is unphysical, so the data do not determine this step's T-dependence) |
| `k_glc_other` | 146.1 | 26 – 266 | SLOPPY (the 80-120 C window is too narrow for this step) |

### Sloppy directions

Hessian condition number 5.58e+08. Three sloppiest eigendirections:

* eigenvalue `2.28e-05` — Ea k_dfg_other (-0.98), Ea k_1dg_out (+0.11), Ea k_dfg_1dg (+0.11)
* eigenvalue `0.000277` — Ea k_glc_other (+0.98), Ea k_1dg_out (-0.10), Ea k_dfg_1dg (-0.10)
* eigenvalue `0.00131` — Ea k_1dg_out (+0.68), Ea k_dfg_1dg (+0.61), Ea k_3dg_out (-0.27)

## (c) Cross-checks

### Brands (2002) — INDEPENDENT

INDEPENDENT of this fit: different amine (protein-bound lysine epsilon-amino, not free glycine), different sugar:amine ratio (10:1, not 1:1), different author, different data. None of it entered the objective.

| step | fitted @120 °C | Brands @120 °C | ratio | fitted Ea | Brands Ea |
|---|---:|---:|---:|---:|---:|
| condensation (sugar + amine -> Amadori) | 8.037e-05 | 0.00024 | **0.33×** | 94 | 114 |
| TOTAL Amadori degradation | 0.189 | 0.2805 | **0.67×** | — | 126 |

### Martins (2003) M4 — REPRODUCIBILITY, NOT INDEPENDENCE

NOT INDEPENDENT. Martins fitted these constants to the same measurements this script fits. Agreement is a reproducibility result -- a different implementation, scheme and objective recovering the same rates -- and is NOT evidence from a second experiment.

Bimolecular comparators are converted to pseudo-first-order at the experiment's own 200 mmol/L glycine before comparison; comparing them raw understates them 200-fold.

| fitted step | Martins step(s) | fitted k(100 °C) | Martins comparator | ratio | ΔEa kJ/mol |
|---|---|---:|---:|---:|---:|
| `k_schiff` | 1 | 1.71e-05 | 1.61e-05 | **1.06×** | -3 |
| `k_dfg_3dg` | 4 | 0.006462 | 0.0111 | **0.58×** | +48 |
| `k_dfg_1dg` | 7 | 0.03048 | 0.0157 | **1.94×** | -34 |
| `k_dfg_other` | 6 | 0.004312 | 0.00708 | **0.61×** | -35 |
| `k_3dg_out` | 5 + 9 | 0.1047 | 0.1969 | **0.53×** | +52 |
| `k_1dg_out` | 8 | 3.031 | 1.45 | **2.09×** | -36 |
| `k_glc_other` | 2 | 0.000271 | 0.00164 | **0.17×** | +23 |

## (d) Schiff / Amadori

Pseudo-first-order at 200 mmol/L glycine, 100 °C: Schiff `0.003421` /min versus Amadori `0.1536` /min — ratio `44.9`.

| Amadori/Schiff ratio | Δcost | χ² statistic | rejected at 95%? |
|---:|---:|---:|---|
| 0.1 | 37441.125 | 74882.25 | **yes** |
| 0.3 | 18454.383 | 36908.77 | **yes** |
| 1 | 7121.461 | 14242.92 | **yes** |
| 3 | 2514.520 | 5029.04 | **yes** |
| 10 | 563.349 | 1126.70 | **yes** |
| 20 | 115.413 | 230.83 | **yes** |
| 30 | 22.317 | 44.63 | **yes** |
| 35 | 2.274 | 4.55 | **yes** |
| 40 | 1.754 | 3.51 | no |
| 45 | 0.059 | 0.12 | no |
| 50 | 7.493 | 14.99 | **yes** |
| 55 | 16.251 | 32.50 | **yes** |
| 60 | 36.458 | 72.92 | **yes** |
| 80 | 46.722 | 93.44 | **yes** |
| 100 | 101.556 | 203.11 | **yes** |
| 300 | 312.215 | 624.43 | **yes** |
| 1e+03 | 421.597 | 843.19 | **yes** |
| 1e+04 | 455.653 | 911.31 | **yes** |
| 1e+05 | 457.224 | 914.45 | **yes** |
| 1e+06 | 462.716 | 925.43 | **yes** |

**Not rejected at 95%: ratio 40 – 45.** FAST_BARRIERS implies 5.0e-05 (Schiff faster); `arrhenius_params.yml` implies 3.3e+04 (Amadori faster). **Both are rejected — but only one of them has the sign right.**

## (f) Propagation to the screening lane

**Applied: False.** The mapping is principled but its dominant uncertainty is an ASSUMED prefactor, not the fitted rate. Setting A = 1e11 s^-1 for every step folds the real prefactor, the activation entropy and the buffer catalysis into Ea at 1.71 kcal/mol per decade -- wider than every fitted confidence interval in this report. Applying it would move shipped predictions on the strength of an assumption rather than on the strength of the measurement. Published with the arithmetic so the owner can decide; FAST_BARRIERS is UNCHANGED by this wave. The counterfactual effect of applying it was measured anyway and is recorded in the wave ledger.

Prefactor used: `1e+11` — recommend.py:1216 calls get_reference_pre_exponential() with no argument, so the screening lane uses 1e11 for every family; the per-family YAML A values cancel exactly in benchmark_validation.py:662-673 and are dead in this lane.

| step | family | derived kcal/mol | shipped kcal/mol | Δ |
|---|---|---:|---:|---:|
| `k_schiff` | schiff_condensation | *excluded* | 15.0 | — |
| `k_amadori` | amadori_rearrangement | 23.20 | 23.00 | +0.20 |
| `k_dfg_3dg` | enolisation_intermediate | 25.55 | 21.00 | +4.55 |
| `k_dfg_1dg` | 2,3-enolisation | 24.40 | 28.00 | -3.60 |
| `k_3dg_out` | 1,2-enolisation | 23.49 | 28.00 | -4.51 |
| `k_1dg_out` | furanone_amino_acid_reduction | 20.99 | 28.00 | -7.01 |

## Per-response fit

| experiment | species | T °C | n | median fold error |
|---|---|---:|---:|---:|
| dfgdeg | DFG | 100 | 9 | 1.381× |
| dfgdeg | DFG | 120 | 9 | 1.238× |
| gg | 1-DG | 80 | 11 | 1.065× |
| gg | 1-DG | 100 | 14 | 1.157× |
| gg | 1-DG | 120 | 9 | 1.125× |
| gg | 3-DG | 80 | 16 | 1.071× |
| gg | 3-DG | 100 | 12 | 1.079× |
| gg | 3-DG | 120 | 8 | 1.063× |
| gg | DFG | 80 | 16 | 1.060× |
| gg | DFG | 100 | 15 | 1.044× |
| gg | DFG | 120 | 13 | 1.024× |
| gg | glucose | 80 | 17 | 1.012× |
| gg | glucose | 100 | 15 | 1.025× |
| gg | glucose | 120 | 12 | 1.020× |

## Out-of-objective diagnostics

```json
{
  "glycine_release_yield": [
    {
      "T_C": 100.0,
      "measured_final_glycine_mmol_L": 8.34,
      "predicted_final_glycine_mmol_L": 8.949383742481544,
      "measured_yield_pct": 83.4,
      "predicted_yield_pct": 89.49383742481544
    },
    {
      "T_C": 120.0,
      "measured_final_glycine_mmol_L": 7.81,
      "predicted_final_glycine_mmol_L": 9.009259592427,
      "measured_yield_pct": 78.1,
      "predicted_yield_pct": 90.09259592427
    }
  ],
  "ph55_out_of_objective": [
    {
      "T_C": 100.0,
      "pH": 5.5,
      "n": 10,
      "median_fold_error": 6.475857851815408
    },
    {
      "T_C": 120.0,
      "pH": 5.5,
      "n": 9,
      "median_fold_error": 14.499999486021016
    }
  ],
  "glucose_budget_100C_150min": {
    "total_glucose_lost_mmol_L": 80.25551604441982,
    "via_condensation_mmol_L": 73.91516475724251,
    "via_declared_nuisance_channel_mmol_L": 6.342280552596105,
    "nuisance_share": 0.07902420178104792
  }
}
```

