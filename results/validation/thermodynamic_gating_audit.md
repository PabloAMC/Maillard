# Thermodynamic Gating Audit

| Benchmark | Path | Applicable | Baseline Status | Gated Status | Δ MAE ppb | Δ Max Ratio | Material | Recommended Policy | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| acrylamide_asparagine_glucose_Parker2012 | free_precursor | no | n/a | n/a | n/a | n/a | no | diagnostic_only | Thermodynamic gating audit is currently only meaningful for non-safety free-precursor FAST benchmarks. |
| cys_glucose_150C_Farmer1999 | free_precursor | yes | pass | pass | -33.82 | -0.015 | no | diagnostic_only | Thermodynamic gating does not materially improve benchmark error under the current threshold and remains diagnostic-only. |
| cys_ribose_140C_Hofmann1998 | free_precursor | yes | partial-pass | scale-gap | -59.97 | -5.597 | no | diagnostic_only | Thermodynamic gating does not materially improve benchmark error under the current threshold and remains diagnostic-only. |
| cys_ribose_150C_Mottram1994 | free_precursor | yes | pass | scale-gap | 55.18 | -0.434 | no | diagnostic_only | Thermodynamic gating does not materially improve benchmark error under the current threshold and remains diagnostic-only. |
| pea_isolate_40C_PratapSingh2021 | matrix_only | no | n/a | n/a | n/a | n/a | no | diagnostic_only | Thermodynamic gating audit is currently only meaningful for non-safety free-precursor FAST benchmarks. |
| soy_isolate_40C_PratapSingh2021 | matrix_only | no | n/a | n/a | n/a | n/a | no | diagnostic_only | Thermodynamic gating audit is currently only meaningful for non-safety free-precursor FAST benchmarks. |

Audited benchmarks: 6
Material improvements: 0
