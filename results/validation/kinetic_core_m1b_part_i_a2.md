# M1b part (i): Amendment A2 on the stored draws (prereg v3.8 secs. 21.1, 21.4)

*`scripts/generators/fit_m1b_ph_balance.py --a2-noise`: the buffered-pot noise term, the verdict per choice, per-signature false-kill rates. One engine run (beta_ref). A1's records stand beside it.*

## sigma_max from Meynier's 20-pot statement

| edge | sigma_max (MC) | P(all 20 within) there | closed form |
|---|---|---|---|
| 0.2 | 0.135 | 0.0503 | 0.1352 |
| 0.25 | 0.168 | 0.0524 | 0.1690 |

## Buffer capacities (cooled, mmol/L per pH)

| pot | pH | I (M) | beta, I = 0 | beta, Davies |
|---|---|---|---|---|
| Meynier Cys label 4.5, t = 0 | 4.50 | 0.615 | 5.6 | 13.2 |
| Meynier Cys label 5, t = 0 | 5.00 | 0.646 | 9.4 | 33.4 |
| Meynier Cys label 5.5, t = 0 | 5.50 | 0.721 | 25.1 | 74.5 |
| Meynier Cys label 6, t = 0 | 6.00 | 0.851 | 62.6 | 113.6 |
| Meynier Cys label 6.5, t = 0 | 6.50 | 1.001 | 110.1 | 107.8 |
| Zhou ARP label 8, B0-pinned 60 min (beta_ref) | 7.10 | 0.018 | 2.42 | 2.38 |

## The choices

| choice | sd per label (4.5 .. 6.5) | edge | in D | basis |
|---|---|---|---|---|
| sd_0.071 | 0.071 / 0.071 / 0.071 / 0.071 / 0.071 | 0.2 | True | meter only (A1) |
| sd_0.100 | 0.100 / 0.100 / 0.100 / 0.100 / 0.100 | 0.2 | True | a point inside D |
| sigma_max | 0.135 / 0.135 / 0.135 / 0.135 / 0.135 | 0.2 | True | Meynier's 20-pot statement, edge 0.2 |
| sigma_max_edge0.25 | 0.168 / 0.168 / 0.168 / 0.168 / 0.168 | 0.25 | True | Meynier's statement read at 0.1 resolution, edge 0.25 |
| beta_s_x0.5 | 0.074 / 0.071 / 0.071 / 0.071 / 0.071 | 0.2 | True | A1's shape term through the Davies buffer-capacity ratio |
| beta_s_x1 | 0.082 / 0.073 / 0.071 / 0.071 / 0.071 | 0.2 | True | A1's shape term through the Davies buffer-capacity ratio |
| beta_s_x2 | 0.108 / 0.078 / 0.072 / 0.071 / 0.071 | 0.2 | True | A1's shape term through the Davies buffer-capacity ratio |
| shape_unscaled_s_x1 | 0.237 / 0.237 / 0.237 / 0.237 / 0.237 | 0.2 | False | the review's form: A1's shape term put on Meynier unscaled |
| shape_unscaled_s_x2 | 0.458 / 0.458 / 0.458 / 0.458 / 0.458 | 0.2 | False | the review's form: A1's shape term put on Meynier unscaled |

## On the measured data: draws passing all four (minimum over A1's seven forms)

| rival | sd_0.071 | sd_0.100 | sigma_max | sigma_max_edge0.25 | beta_s_x0.5 | beta_s_x1 | beta_s_x2 | shape_unscaled_s_x1 | shape_unscaled_s_x2 | verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| B0_pinned | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | falsified on the screen's draws (A2: killed at every choice in D) |
| B0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | falsified on the screen's draws (A2: killed at every choice in D) |
| B1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | not falsified (A2: passes at every choice in D) |
| B2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | falsified on the screen's draws (A2: killed at every choice in D) |
| B3 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | falsified on the screen's draws (A2: killed at every choice in D) |
| B13 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | not falsified (A2: passes at every choice in D) |

## Per-signature false-kill rate (sets from the rival's own best draw)

| generation \| family | rival | Hou level | Hou shape | Mey | Zhou | all four | Mey all within |
|---|---|---|---|---|---|---|---|
| sd_0.071|primary | B0 | 0.000 | 0.000 | 0.000 | 0.215 | 0.215 | 0.000 |
| sd_0.071|primary | B1 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.330 |
| sd_0.071|primary | B2 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| sd_0.071|primary | B3 | 0.000 | 0.000 | 0.045 | 0.000 | 0.045 | 0.000 |
| sd_0.071|primary | B13 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.225 |
| sd_0.071|offset_shape | B0 | 0.000 | 0.000 | 0.000 | 0.440 | 0.440 | 0.000 |
| sd_0.071|offset_shape | B1 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.335 |
| sd_0.071|offset_shape | B2 | 0.000 | 0.000 | 0.005 | 0.000 | 0.005 | 0.000 |
| sd_0.071|offset_shape | B3 | 0.000 | 0.000 | 0.065 | 0.010 | 0.075 | 0.000 |
| sd_0.071|offset_shape | B13 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.260 |
| sigma_max|primary | B0 | 0.000 | 0.000 | 0.000 | 0.210 | 0.210 | 0.000 |
| sigma_max|primary | B1 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.135 |
| sigma_max|primary | B2 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| sigma_max|primary | B3 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| sigma_max|primary | B13 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.170 |
| sigma_max|offset_shape | B0 | 0.000 | 0.000 | 0.000 | 0.450 | 0.450 | 0.000 |
| sigma_max|offset_shape | B1 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.135 |
| sigma_max|offset_shape | B2 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| sigma_max|offset_shape | B3 | 0.000 | 0.000 | 0.015 | 0.010 | 0.025 | 0.000 |
| sigma_max|offset_shape | B13 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.185 |

## Per-signature kill power (max over the other truths)

| generation \| family | rival | Hou level | Hou shape | Mey | Zhou | all four |
|---|---|---|---|---|---|---|
| sd_0.071|primary | B0 | 0.000 | 0.000 | 1.000 | 1.000 | 1.000 |
| sd_0.071|primary | B1 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| sd_0.071|primary | B2 | 0.000 | 0.000 | 1.000 | 0.000 | 1.000 |
| sd_0.071|primary | B3 | 0.000 | 0.000 | 1.000 | 1.000 | 1.000 |
| sd_0.071|primary | B13 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| sd_0.071|offset_shape | B0 | 0.000 | 0.005 | 1.000 | 0.990 | 1.000 |
| sd_0.071|offset_shape | B1 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| sd_0.071|offset_shape | B2 | 0.000 | 0.000 | 1.000 | 0.000 | 1.000 |
| sd_0.071|offset_shape | B3 | 0.000 | 0.005 | 1.000 | 0.990 | 1.000 |
| sd_0.071|offset_shape | B13 | 0.000 | 0.000 | 0.000 | 0.000 | 0.005 |
| sigma_max|primary | B0 | 0.000 | 0.000 | 0.980 | 1.000 | 1.000 |
| sigma_max|primary | B1 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| sigma_max|primary | B2 | 0.000 | 0.000 | 0.980 | 0.000 | 0.980 |
| sigma_max|primary | B3 | 0.000 | 0.000 | 0.980 | 1.000 | 1.000 |
| sigma_max|primary | B13 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| sigma_max|offset_shape | B0 | 0.000 | 0.005 | 0.990 | 0.995 | 0.995 |
| sigma_max|offset_shape | B1 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| sigma_max|offset_shape | B2 | 0.000 | 0.000 | 0.990 | 0.000 | 0.990 |
| sigma_max|offset_shape | B3 | 0.000 | 0.005 | 0.990 | 0.995 | 0.995 |
| sigma_max|offset_shape | B13 | 0.000 | 0.000 | 0.000 | 0.000 | 0.010 |
