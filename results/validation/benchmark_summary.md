# Benchmark Summary

| Benchmark | Tier | Family | Protein | Execution Path | Status | Strict Ready | Coverage | Pearson R | Max Ratio | MAE ppb | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| acrylamide_asparagine_glucose_Parker2012 | PRIMARY | safety | free | free_precursor | partial-pass | yes | 100.0% | n/a | 1.033 | 48.59 | validated |
| cys_glucose_150C_Farmer1999 | PRIMARY | free_aa_sugar_discrimination | free | free_precursor | pass | yes | 100.0% | 1.000 | 1.222 | 2.53 | validated |
| cys_ribose_140C_Hofmann1998 | PRIMARY | free_aa_sulfur | free | free_precursor | partial-pass | yes | 100.0% | n/a | 1.442 | 53.38 | validated |
| cys_ribose_150C_Mottram1994 | PRIMARY | free_aa_sulfur | free | free_precursor | pass | yes | 100.0% | 0.995 | 1.381 | 94.61 | validated |
| pea_isolate_40C_PratapSingh2021 | PRIMARY | matrix_headspace | pea_iso | matrix_only | pass | no | 100.0% | 1.000 | 1.002 | 0.34 | matrix-only intake path is executable but not yet in the strict release gate |
| soy_isolate_40C_PratapSingh2021 | PRIMARY | matrix_headspace | soy_iso | matrix_only | pass | no | 100.0% | 1.000 | 1.001 | 0.16 | matrix-only intake path is executable but not yet in the strict release gate |

Supported benchmarks: 6/6
Benchmarks without blocking coverage/ranking gaps: 6/6
Strict-ready benchmarks: 4/6
