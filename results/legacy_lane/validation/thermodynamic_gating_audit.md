# Thermodynamic Gating Audit

| Benchmark | Path | Applicable | Baseline Status | Gated Status | Δ MAE ppb | Δ Max Ratio | Material | Recommended Policy | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| acrylamide_spi_extrusion_130C_ACSRef3 | free_precursor | no | n/a | n/a | n/a | n/a | no | diagnostic_only | Thermodynamic gating audit is currently only meaningful for non-safety free-precursor FAST benchmarks. |
| cml_cel_commercial_pbma_Foods2023 | free_precursor | no | n/a | n/a | n/a | n/a | no | diagnostic_only | Thermodynamic gating audit is currently only meaningful for non-safety free-precursor FAST benchmarks. |
| cys_ribose_140C_Hofmann1998 | free_precursor | yes | scale-gap | scale-gap | 0.15 | 0.497 | yes | benchmark_facing_candidate | Thermodynamic gating materially improves benchmark error without degrading supported coverage/status. |
| furosine_extrusion_crossover_140C_RamirezJimenez2000 | free_precursor | no | n/a | n/a | n/a | n/a | no | diagnostic_only | Thermodynamic gating audit is currently only meaningful for non-safety free-precursor FAST benchmarks. |
| hofmann1998_c2c3_recombination_145C_20min_pH3 | free_precursor | yes | scale-gap | scale-gap | 0.00 | 0.000 | no | diagnostic_only | Thermodynamic gating does not materially improve benchmark error under the current threshold and remains diagnostic-only. |
| hofmann1998_c2c3_recombination_145C_20min_pH5 | free_precursor | yes | scale-gap | scale-gap | 0.00 | 0.000 | no | diagnostic_only | Thermodynamic gating does not materially improve benchmark error under the current threshold and remains diagnostic-only. |
| hofmann1998_c2c3_recombination_145C_20min_pH7 | free_precursor | yes | scale-gap | scale-gap | 0.00 | 0.000 | no | diagnostic_only | Thermodynamic gating does not materially improve benchmark error under the current threshold and remains diagnostic-only. |
| hofmann1998_fructose_cysteine_145C_20min_pH5 | free_precursor | yes | scale-gap | scale-gap | -1.24 | -0.459 | no | diagnostic_only | Thermodynamic gating does not materially improve benchmark error under the current threshold and remains diagnostic-only. |
| hofmann1998_furan2aldehyde_h2s_145C_20min_pH5 | free_precursor | yes | scale-gap | scale-gap | 0.00 | 0.000 | no | diagnostic_only | Thermodynamic gating does not materially improve benchmark error under the current threshold and remains diagnostic-only. |
| hofmann1998_glucose_cysteine_145C_20min_pH5 | free_precursor | yes | scale-gap | scale-gap | 261.68 | 18.900 | yes | benchmark_facing_candidate | Thermodynamic gating materially improves benchmark error without degrading supported coverage/status. |
| hofmann1998_norfuraneol_cysteine_145C_20min_pH5 | free_precursor | yes | scale-gap | scale-gap | 0.00 | 0.000 | no | diagnostic_only | Thermodynamic gating does not materially improve benchmark error under the current threshold and remains diagnostic-only. |
| hofmann1998_norfuraneol_h2s_145C_20min_pH5 | free_precursor | yes | scale-gap | scale-gap | -106.50 | -0.206 | no | diagnostic_only | Thermodynamic gating does not materially improve benchmark error under the current threshold and remains diagnostic-only. |
| hofmann1998_ribose_cysteine_145C_20min_pH5 | free_precursor | yes | scale-gap | scale-gap | 391.81 | 7.441 | yes | benchmark_facing_candidate | Thermodynamic gating materially improves benchmark error without degrading supported coverage/status. |
| pea_isolate_40C_PratapSingh2021 | matrix_only | no | n/a | n/a | n/a | n/a | no | diagnostic_only | Thermodynamic gating audit is currently only meaningful for non-safety free-precursor FAST benchmarks. |
| pea_isolate_ribose_cysteine_100C_45min_Internal2026 | matrix_precursor_augmented | no | n/a | n/a | n/a | n/a | no | diagnostic_only | Thermodynamic gating audit is currently only meaningful for non-safety free-precursor FAST benchmarks. |
| pea_isolate_uht_140C_Trikusuma2019 | matrix_only | no | n/a | n/a | n/a | n/a | no | diagnostic_only | Thermodynamic gating audit is currently only meaningful for non-safety free-precursor FAST benchmarks. |
| resconi_2023_pbma_beef_identity_benchmark | free_precursor | yes | scale-gap | scale-gap | -328.34 | -0.459 | no | diagnostic_only | Thermodynamic gating does not materially improve benchmark error under the current threshold and remains diagnostic-only. |
| soy_isolate_40C_PratapSingh2021 | matrix_only | no | n/a | n/a | n/a | n/a | no | diagnostic_only | Thermodynamic gating audit is currently only meaningful for non-safety free-precursor FAST benchmarks. |
| soy_isolate_ribose_cysteine_100C_45min_Internal2026 | matrix_precursor_augmented | no | n/a | n/a | n/a | n/a | no | diagnostic_only | Thermodynamic gating audit is currently only meaningful for non-safety free-precursor FAST benchmarks. |
| thiamine_cys_glucose_120C_Bolton1994 | free_precursor | yes | scale-gap | scale-gap | 0.01 | 5481.235 | yes | benchmark_facing_candidate | Thermodynamic gating materially improves benchmark error without degrading supported coverage/status. |
| thiamine_cys_xylose_145C_Cerny2008 | free_precursor | yes | scale-gap | scale-gap | 0.65 | 20.765 | yes | benchmark_facing_candidate | Thermodynamic gating materially improves benchmark error without degrading supported coverage/status. |

Audited benchmarks: 21
Material improvements: 5
