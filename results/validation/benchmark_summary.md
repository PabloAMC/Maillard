# Benchmark Summary

| Benchmark | Tier | Family | Protein | Process State | Execution Path | Engine | Cantera Role | Thermo Policy | Ranking Contract | Status | Strict Ready | Coverage | Pearson R | Max Ratio | MAE ppb | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| acrylamide_asparagine_glucose_Parker2012 | PRIMARY | safety | free | n/a | free_precursor | fast_observable | diagnostic_reference_only | diagnostic_only | n/a | pass-no-ranking | yes | 100.0% | n/a | 1.033 | 48.59 | validated |
| cys_glucose_150C_Farmer1999 | PRIMARY | free_aa_sugar_discrimination | free | n/a | free_precursor | fast_observable | diagnostic_reference_only | diagnostic_only | n/a | pass | yes | 100.0% | 1.000 | 1.222 | 2.53 | validated |
| cys_ribose_140C_Hofmann1998 | PRIMARY | free_aa_sulfur | free | n/a | free_precursor | fast_observable | diagnostic_reference_only | diagnostic_only | n/a | pass-no-ranking | yes | 100.0% | n/a | 1.442 | 53.39 | validated |
| cys_ribose_150C_Mottram1994 | PRIMARY | free_aa_sulfur | free | n/a | free_precursor | fast_observable | diagnostic_reference_only | diagnostic_only | n/a | pass | yes | 100.0% | 0.995 | 1.380 | 94.57 | validated |
| pea_isolate_40C_PratapSingh2021 | PRIMARY | matrix_headspace | pea_iso | ambient_slurry | matrix_only | matrix_intake_headspace | not_authoritative | not_applicable | pass | pass | no | 100.0% | 1.000 | 1.002 | 0.34 | matrix-only intake path is executable but not yet in the strict release gate |
| pea_isolate_ribose_cysteine_100C_45min_Internal2026 | SECONDARY | matrix_precursor_augmented | pea_iso | aqueous_pre_extrusion_model | matrix_precursor_augmented | fast_observable_matrix | not_authoritative | not_applicable | pass | pass | no | 100.0% | 1.000 | 1.000 | 0.00 | matrix-only intake path is executable but not yet in the strict release gate |
| soy_isolate_40C_PratapSingh2021 | PRIMARY | matrix_headspace | soy_iso | ambient_slurry | matrix_only | matrix_intake_headspace | not_authoritative | not_applicable | pass | pass | no | 100.0% | 1.000 | 1.001 | 0.16 | matrix-only intake path is executable but not yet in the strict release gate |
| soy_isolate_ribose_cysteine_100C_45min_Internal2026 | SECONDARY | matrix_precursor_augmented | soy_iso | aqueous_pre_extrusion_model | matrix_precursor_augmented | fast_observable_matrix | not_authoritative | not_applicable | pass | pass | no | 100.0% | 1.000 | 1.000 | 0.00 | matrix-only intake path is executable but not yet in the strict release gate |

Supported benchmarks: 8/8
Benchmarks without blocking coverage/ranking gaps: 8/8
Strict-ready benchmarks: 4/8
