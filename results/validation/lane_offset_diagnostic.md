# Per-lane offset diagnostic (wave B46)

Signed offset `dex = log10(predicted / measured)`. Positive means the model reads high.
A diagnostic, not a fit: nothing here is scored and no constant moves.

| lane | n | pots | median signed dex | reads | sign consistency | systematic | median abs dex |
|---|---:|---:|---:|---|---:|---|---:|
| sulfur | 19 | 11 | +0.84 | high | 79% | YES | 1.47 |
| trunk | 6 | 1 | -0.90 | low | 83% | YES | 0.90 |
| acrylamide | 11 | 7 | -0.44 | low | 73% | no | 0.80 |
| lipid | 8 | 3 | -0.45 | low | 88% | YES | 0.45 |

## What each lane's offset tracks

Spearman rank correlation of the signed offset against each stated condition. |rho| >= 0.6 is the pre-registered threshold for 'tracks'. **n is small everywhere; read the correlations with it.**

| lane | covariate | rho | n | distinct levels | tracks (pre-registered) | trend supported by the design |
|---|---|---:|---:|---:|---|---|
| acrylamide | temp_C | +0.52 | 11 | 2 |  |  |
| acrylamide | time_min | -0.13 | 11 | 4 |  |  |
| acrylamide | ph | -0.10 | 11 | 3 |  |  |
| acrylamide | water_activity | +0.59 | 11 | 3 |  |  |
| lipid | temp_C | -0.87 | 8 | 2 | YES | NO -- only 2 distinct value(s) of temp_C in this lane: a rank correlation here is a 2-group comparison, not a trend |
| lipid | time_min | -0.81 | 8 | 3 | YES | YES |
| lipid | ph | +0.81 | 8 | 3 | YES | YES |
| lipid | water_activity | n/a | 4 | 1 |  |  |
| sulfur | temp_C | -0.82 | 19 | 6 | YES | YES |
| sulfur | time_min | +0.70 | 19 | 5 | YES | YES |
| sulfur | ph | +0.13 | 19 | 5 |  |  |
| sulfur | water_activity | +0.14 | 19 | 3 |  |  |
| trunk | temp_C | n/a | 6 | 1 |  |  |
| trunk | time_min | n/a | 6 | 1 |  |  |
| trunk | ph | n/a | 6 | 1 |  |  |
| trunk | water_activity | n/a | 6 | 1 |  |  |

## How to read a hit

A lane that is systematic AND tracks a covariate is where a missing process is most likely to live, because parameter uncertainty has already been ruled out globally: `core_prediction_uncertainty.json` reports that uncapping every prior moves coverage from 19 % to 21 % against a nominal 90 %.

It does NOT say which process. An offset that grows with temperature fits a wrong barrier, a wrong Q10, a missing temperature-dependent channel, or a measurement whose efficiency changes with temperature. Wave B45 met the last of those.

