# B56 structural-discrepancy term (in-sample fit residuals only)

Pre-registration: `results/validation/kinetic_core_b55_b56_prereg.md`. sigma_L^2 = sum_f c_f sum_i r_i^2 / sum_f n_used_f, c_f = n_f/(n_f - p_f); r_i the log10 residual at the frozen optimum of a two-sided row the fit read. A DECLARED MODEL-ERROR ASSUMPTION estimated in sample only, not a measurement; the parameter share is the dof correction and is not subtracted again.

| lane | sigma (log10) | RMS uncorrected | rows used | fits |
|---|---:|---:|---:|---|
| sulfur | 0.491 | 0.372 | 49 | kinetic_core_b9_laplace_covariance.json (n=54, p=23, used 49, c=1.74) |
| trunk | 0.135 | 0.132 | 379 | kinetic_core_b1_fit_report.json (n=354, p=8, used 354, c=1.02); kinetic_core_b7_fit_report.json (n=3, p=1, used 3, c=1.50); kinetic_core_b18_fit_report.json (n=10, p=6, used 10, c=2.50); kinetic_core_b41_fit_report.json (n=12, p=5, used 12, c=1.71) |
| acrylamide | 1.188 | 0.945 | 17 | kinetic_core_b3_fit_report.json (n=30, p=11, used 17, c=1.58) |
| lipid | 0.120 | 0.069 | 18 | kinetic_core_b6_fit_report.json (n=18, p=12, used 18, c=3.00) |
