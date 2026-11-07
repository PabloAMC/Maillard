# M1b part (i): S3's posterior-predictive test (prereg v3.8 sec. 21.3)

*p per signature at the converged fit (Laplace predictive), Holm across the four (alpha 0.05). Decisive: B0 at F-A. Shown: the primary A1 form and the offset / shape form at s, per A2 choice.*

## B0_attribution_fit: falsified under S3 (inadequate after Holm in every cell in D)

Laplace chi2_red 2.548; unidentified (infinite sigma, not propagated): none; failed pot runs 0.

| cell | in D | Hou level | Hou shape | Mey | Zhou | rejected after Holm |
|---|---|---|---|---|---|---|
| primary|sd_0.071 | True | 0.868 | 4.93e-09 | 0.0001 | 0.0506 | sigma_hou_shape, sigma_mey |
| primary|sd_0.100 | True | 0.868 | 4.93e-09 | 0.000105 | 0.0506 | sigma_hou_shape, sigma_mey |
| primary|sigma_max | True | 0.868 | 4.93e-09 | 0.00012 | 0.0506 | sigma_hou_shape, sigma_mey |
| primary|sigma_max_edge0.25 | True | 0.868 | 4.93e-09 | 0.000905 | 0.0506 | sigma_hou_shape, sigma_mey |
| primary|beta_s_x0.5 | True | 0.868 | 4.93e-09 | 7.5e-05 | 0.0506 | sigma_hou_shape, sigma_mey |
| primary|beta_s_x1 | True | 0.868 | 4.93e-09 | 6e-05 | 0.0506 | sigma_hou_shape, sigma_mey |
| primary|beta_s_x2 | True | 0.868 | 4.93e-09 | 9e-05 | 0.0506 | sigma_hou_shape, sigma_mey |
| primary|shape_unscaled_s_x1 | False | 0.868 | 4.93e-09 | 0.000295 | 0.0506 | sigma_hou_shape, sigma_mey |
| primary|shape_unscaled_s_x2 | False | 0.868 | 4.93e-09 | 0.000405 | 0.0506 | sigma_hou_shape, sigma_mey |
| offset_shape_s|sd_0.071 | True | 0.871 | 0.0162 | 5e-05 | 0.0744 | sigma_mey, sigma_hou_shape |
| offset_shape_s|sd_0.100 | True | 0.871 | 0.0162 | 7e-05 | 0.0744 | sigma_mey, sigma_hou_shape |
| offset_shape_s|sigma_max | True | 0.871 | 0.0162 | 8.5e-05 | 0.0744 | sigma_mey, sigma_hou_shape |
| offset_shape_s|sigma_max_edge0.25 | True | 0.871 | 0.0162 | 0.0007 | 0.0744 | sigma_mey, sigma_hou_shape |
| offset_shape_s|beta_s_x0.5 | True | 0.871 | 0.0162 | 7.5e-05 | 0.0744 | sigma_mey, sigma_hou_shape |
| offset_shape_s|beta_s_x1 | True | 0.871 | 0.0162 | 6.5e-05 | 0.0744 | sigma_mey, sigma_hou_shape |
| offset_shape_s|beta_s_x2 | True | 0.871 | 0.0162 | 5e-05 | 0.0744 | sigma_mey, sigma_hou_shape |
| offset_shape_s|shape_unscaled_s_x1 | False | 0.871 | 0.0162 | 0.00033 | 0.0744 | sigma_mey, sigma_hou_shape |
| offset_shape_s|shape_unscaled_s_x2 | False | 0.871 | 0.0162 | 0.00047 | 0.0744 | sigma_mey, sigma_hou_shape |

## B0_confirmation_fit: falsified under S3 (inadequate after Holm in every cell in D)

Laplace chi2_red 1.250; unidentified (infinite sigma, not propagated): none; failed pot runs 0.

| cell | in D | Hou level | Hou shape | Mey | Zhou | rejected after Holm |
|---|---|---|---|---|---|---|
| primary|sd_0.071 | True | 0.835 | 2.08e-11 | 0.0008 | 2.06e-05 | sigma_hou_shape, sigma_zhou, sigma_mey |
| primary|sd_0.100 | True | 0.835 | 2.08e-11 | 0.00095 | 2.06e-05 | sigma_hou_shape, sigma_zhou, sigma_mey |
| primary|sigma_max | True | 0.835 | 2.08e-11 | 0.00105 | 2.06e-05 | sigma_hou_shape, sigma_zhou, sigma_mey |
| primary|sigma_max_edge0.25 | True | 0.835 | 2.08e-11 | 0.00283 | 2.06e-05 | sigma_hou_shape, sigma_zhou, sigma_mey |
| primary|beta_s_x0.5 | True | 0.835 | 2.08e-11 | 0.000755 | 2.06e-05 | sigma_hou_shape, sigma_zhou, sigma_mey |
| primary|beta_s_x1 | True | 0.835 | 2.08e-11 | 0.000975 | 2.06e-05 | sigma_hou_shape, sigma_zhou, sigma_mey |
| primary|beta_s_x2 | True | 0.835 | 2.08e-11 | 0.00123 | 2.06e-05 | sigma_hou_shape, sigma_zhou, sigma_mey |
| primary|shape_unscaled_s_x1 | False | 0.835 | 2.08e-11 | 0.0008 | 2.06e-05 | sigma_hou_shape, sigma_zhou, sigma_mey |
| primary|shape_unscaled_s_x2 | False | 0.835 | 2.08e-11 | 0.00026 | 2.06e-05 | sigma_hou_shape, sigma_zhou, sigma_mey |
| offset_shape_s|sd_0.071 | True | 0.842 | 0.014 | 0.000755 | 0.00141 | sigma_mey, sigma_zhou, sigma_hou_shape |
| offset_shape_s|sd_0.100 | True | 0.842 | 0.014 | 0.00105 | 0.00141 | sigma_mey, sigma_zhou, sigma_hou_shape |
| offset_shape_s|sigma_max | True | 0.842 | 0.014 | 0.00108 | 0.00141 | sigma_mey, sigma_zhou, sigma_hou_shape |
| offset_shape_s|sigma_max_edge0.25 | True | 0.842 | 0.014 | 0.00296 | 0.00141 | sigma_zhou, sigma_mey, sigma_hou_shape |
| offset_shape_s|beta_s_x0.5 | True | 0.842 | 0.014 | 0.00079 | 0.00141 | sigma_mey, sigma_zhou, sigma_hou_shape |
| offset_shape_s|beta_s_x1 | True | 0.842 | 0.014 | 0.001 | 0.00141 | sigma_mey, sigma_zhou, sigma_hou_shape |
| offset_shape_s|beta_s_x2 | True | 0.842 | 0.014 | 0.00126 | 0.00141 | sigma_mey, sigma_zhou, sigma_hou_shape |
| offset_shape_s|shape_unscaled_s_x1 | False | 0.842 | 0.014 | 0.0007 | 0.00141 | sigma_mey, sigma_zhou, sigma_hou_shape |
| offset_shape_s|shape_unscaled_s_x2 | False | 0.842 | 0.014 | 0.000235 | 0.00141 | sigma_mey, sigma_zhou, sigma_hou_shape |

## B2_confirmation_fit: falsified under S3 (inadequate after Holm in every cell in D)

Laplace chi2_red 1.339; unidentified (infinite sigma, not propagated): none; failed pot runs 0.

| cell | in D | Hou level | Hou shape | Mey | Zhou | rejected after Holm |
|---|---|---|---|---|---|---|
| primary|sd_0.071 | True | 0.899 | 0.263 | 0 | 0.127 | sigma_mey |
| primary|sd_0.100 | True | 0.899 | 0.263 | 0 | 0.127 | sigma_mey |
| primary|sigma_max | True | 0.899 | 0.263 | 1.5e-05 | 0.127 | sigma_mey |
| primary|sigma_max_edge0.25 | True | 0.899 | 0.263 | 0.000255 | 0.127 | sigma_mey |
| primary|beta_s_x0.5 | True | 0.899 | 0.263 | 0 | 0.127 | sigma_mey |
| primary|beta_s_x1 | True | 0.899 | 0.263 | 0 | 0.127 | sigma_mey |
| primary|beta_s_x2 | True | 0.899 | 0.263 | 5e-06 | 0.127 | sigma_mey |
| primary|shape_unscaled_s_x1 | False | 0.899 | 0.263 | 0.00018 | 0.127 | sigma_mey |
| primary|shape_unscaled_s_x2 | False | 0.899 | 0.263 | 0.000195 | 0.127 | sigma_mey |
| offset_shape_s|sd_0.071 | True | 0.899 | 0.308 | 0 | 0.145 | sigma_mey |
| offset_shape_s|sd_0.100 | True | 0.899 | 0.308 | 0 | 0.145 | sigma_mey |
| offset_shape_s|sigma_max | True | 0.899 | 0.308 | 1.5e-05 | 0.145 | sigma_mey |
| offset_shape_s|sigma_max_edge0.25 | True | 0.899 | 0.308 | 0.00021 | 0.145 | sigma_mey |
| offset_shape_s|beta_s_x0.5 | True | 0.899 | 0.308 | 0 | 0.145 | sigma_mey |
| offset_shape_s|beta_s_x1 | True | 0.899 | 0.308 | 0 | 0.145 | sigma_mey |
| offset_shape_s|beta_s_x2 | True | 0.899 | 0.308 | 0 | 0.145 | sigma_mey |
| offset_shape_s|shape_unscaled_s_x1 | False | 0.899 | 0.308 | 0.00016 | 0.145 | sigma_mey |
| offset_shape_s|shape_unscaled_s_x2 | False | 0.899 | 0.308 | 0.00017 | 0.145 | sigma_mey |

## B3_confirmation_fit: falsified under S3 (inadequate after Holm in every cell in D)

Laplace chi2_red 1.083; unidentified (infinite sigma, not propagated): none; failed pot runs 0.

| cell | in D | Hou level | Hou shape | Mey | Zhou | rejected after Holm |
|---|---|---|---|---|---|---|
| primary|sd_0.071 | True | 0.955 | 1.51e-13 | 0 | 0.275 | sigma_mey, sigma_hou_shape |
| primary|sd_0.100 | True | 0.955 | 1.51e-13 | 0 | 0.275 | sigma_mey, sigma_hou_shape |
| primary|sigma_max | True | 0.955 | 1.51e-13 | 0 | 0.275 | sigma_mey, sigma_hou_shape |
| primary|sigma_max_edge0.25 | True | 0.955 | 1.51e-13 | 1e-05 | 0.275 | sigma_hou_shape, sigma_mey |
| primary|beta_s_x0.5 | True | 0.955 | 1.51e-13 | 0 | 0.275 | sigma_mey, sigma_hou_shape |
| primary|beta_s_x1 | True | 0.955 | 1.51e-13 | 0 | 0.275 | sigma_mey, sigma_hou_shape |
| primary|beta_s_x2 | True | 0.955 | 1.51e-13 | 0 | 0.275 | sigma_mey, sigma_hou_shape |
| primary|shape_unscaled_s_x1 | False | 0.955 | 1.51e-13 | 0 | 0.275 | sigma_mey, sigma_hou_shape |
| primary|shape_unscaled_s_x2 | False | 0.955 | 1.51e-13 | 1e-05 | 0.275 | sigma_hou_shape, sigma_mey |
| offset_shape_s|sd_0.071 | True | 0.958 | 0.0116 | 0 | 0.325 | sigma_mey, sigma_hou_shape |
| offset_shape_s|sd_0.100 | True | 0.958 | 0.0116 | 0 | 0.325 | sigma_mey, sigma_hou_shape |
| offset_shape_s|sigma_max | True | 0.958 | 0.0116 | 0 | 0.325 | sigma_mey, sigma_hou_shape |
| offset_shape_s|sigma_max_edge0.25 | True | 0.958 | 0.0116 | 0 | 0.325 | sigma_mey, sigma_hou_shape |
| offset_shape_s|beta_s_x0.5 | True | 0.958 | 0.0116 | 0 | 0.325 | sigma_mey, sigma_hou_shape |
| offset_shape_s|beta_s_x1 | True | 0.958 | 0.0116 | 0 | 0.325 | sigma_mey, sigma_hou_shape |
| offset_shape_s|beta_s_x2 | True | 0.958 | 0.0116 | 0 | 0.325 | sigma_mey, sigma_hou_shape |
| offset_shape_s|shape_unscaled_s_x1 | False | 0.958 | 0.0116 | 0 | 0.325 | sigma_mey, sigma_hou_shape |
| offset_shape_s|shape_unscaled_s_x2 | False | 0.958 | 0.0116 | 3.5e-05 | 0.325 | sigma_mey, sigma_hou_shape |

