# Kinetic core B1 -- HOLD-OUT report (predicted vs measured)

Generated 2026-08-28 on `audit-remediation` @ `05a5fae8` (dirty).

Parameters were frozen by `results/validation/kinetic_core_b1_fit_report.json` BEFORE any hold-out value was read. **Nothing was changed after reading them.**

## Hold-out 1 -- Martins 2005 step-9 browning (epsilon 0.64)

- Held-out response: the melanoidin trajectory at 80/100/120 C, 38 points.
- Held-out constant: epsilon = 0.64 +/- 0.03 L/(mmol*cm) at 470 nm.

### Headline (variant B -- genuinely out-of-sample)

The melanoidin sink constant was estimated from the reactant side only (3-DG and glycine); the browning response never entered the objective.

- **median fold error 1.45x**, max 2.73x
- median signed log10 error +0.161 (the model over-predicts browning)
- within 2x: 95% of points; within 3x: 100%
- **verdict: PASS -- median 1.45x, inside the 1.5x band**

| variant | T (C) | n | median fold error | max fold error | direction | implied epsilon L/(mmol*cm) |
|---|---:|---:|---:|---:|---|---:|
| B (reactant-side sink, out-of-sample) | 80 | 10 | 1.21x | 1.74x | over-predicts | 0.288 |
| B (reactant-side sink, out-of-sample) | 100 | 14 | 1.46x | 1.61x | over-predicts | 0.435 |
| B (reactant-side sink, out-of-sample) | 120 | 14 | 1.47x | 2.73x | over-predicts | 0.432 |
| A (Martins' measured sink, reproducibility) | 80 | 10 | 1.22x | 1.80x | over-predicts | 0.276 |
| A (Martins' measured sink, reproducibility) | 100 | 14 | 1.49x | 1.66x | over-predicts | 0.427 |
| A (Martins' measured sink, reproducibility) | 120 | 14 | 1.49x | 2.81x | over-predicts | 0.429 |

The **implied epsilon** column is the extinction coefficient that would make the model's predicted melanoidin concentration reproduce the measured absorbance. The held-out value is 0.64 +/- 0.03.

### Secondary (variant A -- Martins' own step-9 constant, NOT out-of-sample)

- median fold error 1.47x, max 2.81x, signed log10 +0.169
- **This is a reproducibility check, not a prediction.** Martins fitted this constant to the very response being scored, so a good number here carries much less evidential weight than the same number in variant B.

### Predicted vs measured, point by point (variant B)

| T (C) | t (min) | measured mmol/L | predicted mmol/L | measured A470 | predicted A470 | fold error |
|---:|---:|---:|---:|---:|---:|---:|
| 80 | 30 | 0.000 | 0.003 | 0.000 | 0.002 | 1.01x |
| 80 | 45 | 0.000 | 0.012 | 0.000 | 0.008 | 1.06x |
| 80 | 45 | 0.050 | 0.012 | 0.032 | 0.008 | 1.18x |
| 80 | 60 | 0.000 | 0.031 | 0.000 | 0.020 | 1.15x |
| 80 | 60 | 0.050 | 0.031 | 0.032 | 0.020 | 1.08x |
| 80 | 90 | 0.050 | 0.108 | 0.032 | 0.069 | 1.23x |
| 80 | 120 | 0.100 | 0.245 | 0.064 | 0.157 | 1.48x |
| 80 | 150 | 0.210 | 0.440 | 0.134 | 0.282 | 1.56x |
| 80 | 150 | 0.260 | 0.440 | 0.166 | 0.282 | 1.39x |
| 80 | 180 | 0.310 | 0.689 | 0.198 | 0.441 | 1.74x |
| 100 | 30 | 0.420 | 0.799 | 0.269 | 0.511 | 1.61x |
| 100 | 30 | 0.570 | 0.799 | 0.365 | 0.511 | 1.30x |
| 100 | 30 | 0.680 | 0.799 | 0.435 | 0.511 | 1.14x |
| 100 | 45 | 1.310 | 2.022 | 0.838 | 1.294 | 1.47x |
| 100 | 45 | 1.360 | 2.022 | 0.870 | 1.294 | 1.42x |
| 100 | 60 | 2.210 | 3.548 | 1.414 | 2.271 | 1.56x |
| 100 | 60 | 2.370 | 3.548 | 1.517 | 2.271 | 1.46x |
| 100 | 60 | 2.420 | 3.548 | 1.549 | 2.271 | 1.43x |
| 100 | 90 | 4.590 | 6.889 | 2.938 | 4.409 | 1.48x |
| 100 | 90 | 4.730 | 6.889 | 3.027 | 4.409 | 1.44x |
| 100 | 90 | 4.850 | 6.889 | 3.104 | 4.409 | 1.40x |
| 100 | 120 | 6.320 | 10.155 | 4.045 | 6.499 | 1.59x |
| 100 | 120 | 6.420 | 10.155 | 4.109 | 6.499 | 1.56x |
| 100 | 120 | 6.900 | 10.155 | 4.416 | 6.499 | 1.46x |
| 120 | 5 | 0.050 | 0.483 | 0.032 | 0.309 | 2.73x |
| 120 | 5 | 0.100 | 0.483 | 0.064 | 0.309 | 2.28x |
| 120 | 10 | 1.360 | 2.453 | 0.870 | 1.570 | 1.70x |
| 120 | 10 | 1.690 | 2.453 | 1.082 | 1.570 | 1.40x |
| 120 | 10 | 1.990 | 2.453 | 1.274 | 1.570 | 1.21x |
| 120 | 15 | 3.680 | 5.054 | 2.355 | 3.234 | 1.35x |
| 120 | 15 | 3.740 | 5.054 | 2.394 | 3.234 | 1.33x |
| 120 | 15 | 3.890 | 5.054 | 2.490 | 3.234 | 1.28x |
| 120 | 20 | 4.630 | 7.682 | 2.963 | 4.916 | 1.63x |
| 120 | 20 | 5.110 | 7.682 | 3.270 | 4.916 | 1.48x |
| 120 | 20 | 5.370 | 7.682 | 3.437 | 4.916 | 1.42x |
| 120 | 30 | 8.050 | 12.402 | 5.152 | 7.937 | 1.53x |
| 120 | 30 | 8.310 | 12.402 | 5.318 | 7.937 | 1.48x |
| 120 | 30 | 8.380 | 12.402 | 5.363 | 7.937 | 1.47x |

## Hold-out 2 -- Knol 2005 epsilon = 282 L/(mol*cm), glucose/asparagine

- **NOT PREDICTED -- structural gap, declared**
- the kinetic core carries the melanoidin pool ELEMENTALLY (carbon and nitrogen) and has no representation of the polymer's optical properties. An amine-specific extinction coefficient is not derivable from any mass-action rate, and the corpus contains no measured relation between amine identity and epsilon from which one could be built. Inventing a functional form to score this row would defeat the purpose of holding it out.
- Martins (glycine) 640 vs Knol (asparagine) 282 L/(mol*cm) = 2.27x.
- any model that carries a single epsilon across amino acids is wrong by 2.27x on browning whenever the amine changes between glycine and asparagine.
- the implied epsilons above come from hold-out 1 and are a restatement of that hold-out's fold errors in optical units. They are NOT a prediction of Knol's value and must not be read as one.

## What this does and does not establish

- The reactant-side estimate of the melanoidin sink constant is 0.927x Martins' own value, recovered WITHOUT seeing the browning response. That agreement is the substantive result of this wave: the trunk's own reactant balance determines the browning flux.
- epsilon CANCELS in the concentration comparison, because Martins' melanoidin response is itself A470/epsilon. Hold-out 1 therefore tests the predicted FLUX into the sink, and only the implied-epsilon column restates it in optical units. This is stated so the pass is not over-read as a test of epsilon itself.
- Hold-out 2 (Knol's epsilon) is NOT predicted and cannot be: the module has no optical model. That is a declared structural gap, not a near-miss.
- The whole evaluation sits at pH 6.8 and dilute aqueous conditions. The module has no pH term and no a_w term, so none of these numbers transfer off that point.
