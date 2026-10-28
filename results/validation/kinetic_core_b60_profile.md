# B8 sulfur objective: slice profiles around the frozen optimum

cost at optimum 27.132; 329 evaluations; verdicts {'asymmetric': 6, 'bound_limited': 9, 'flat': 1, 'quadratic': 7}

| coordinate | kind | sigma used | rise -2s | rise +2s | quadratic +2s | asymmetry | slice min k | verdict |
|---|---|---|---|---|---|---|---|---|
| k_tdp_fur | log10k | 0.376 | 21.785 | 19.729 | 20.874 | 1.10 | 0 | **quadratic** |
| k_nf_mft | log10k | 0.274 | 9.709 | 15.742 | 10.434 | 1.62 | 0 | **quadratic** |
| k_nf_mp3p | log10k | 0.578 | 4.287 | 6.802 | 4.320 | 1.59 | 0 | **quadratic** |
| k_mgo_mp | log10k | 0.398 | 3.863 | 3.387 | 3.742 | 1.14 | 0 | **quadratic** |
| k_glc_ha | log10k | 1.05 | 0.000 | 0.000 | 0.000 | - | -3 | **bound_limited** |
| k_dimer_fft | log10k | 1.05 | -0.026 | 4.142 | 0.012 | - | -3 | **bound_limited** |
| k_fft_decay | log10k | 0.315 | 3.866 | 57.556 | 12.713 | 14.89 | 0 | **asymmetric** |
| k_fur_decay | log10k | 0.884 | 3.714 | 6.414 | 13.456 | 1.73 | 0 | **bound_limited** |
| k_pent_caramel | log10k | 0.833 | 7.430 | 18.405 | 7.553 | 2.48 | 0 | **asymmetric** |
| k_cys_thermal | log10k | 0.353 | 3.737 | 5.855 | 4.808 | 1.57 | 0 | **quadratic** |
| k_thiolate_loss | log10k | 0.471 | 12.496 | 617.181 | 65.470 | 49.39 | 0 | **asymmetric** |
| ph_acid_yield_per_sink_event | acid_yield | 0.1 | -0.001 | 2.042 | 1.406 | - | -3 | **bound_limited** |
| ph_arp_secondary_ammonium_pKa | pka | 0.599 | 5.745 | 8.182 | 9.543 | 1.42 | 0 | **quadratic** |
| k_fur_fft | log10k | 1.05 | 0.005 | 2.895 | 0.187 | 532.31 | 0 | **asymmetric** |
| k_osone_decay | log10k | 0.604 | 2.146 | 165.554 | 16.549 | 77.14 | 0 | **asymmetric** |
| k_dimer_decay | log10k | 1.05 | 0.000 | 0.000 | 0.000 | 1.48 | 0 | **bound_limited** |
| k_ttca_deg | log10k | 1.89 | -0.029 | 71.548 | 11.683 | - | -1 | **bound_limited** |
| k_thiol_decay | log10k | 1.05 | -0.000 | 0.000 | 0.000 | - | 3 | **flat** |
| Ea_decay_thiol_sink | Ea | 32.8 | 103.988 | -0.065 | 33.019 | - | 1 | **bound_limited** |
| Ea_decay_carbonyl_sink | Ea | 23 | 2.580 | 0.791 | 1.290 | 3.26 | 0 | **asymmetric** |
| k_dimer_mft | log10k | 0.404 | 14.722 | 0.000 | 5.813 | 30081125.28 | 0 | **bound_limited** |
| k_arp_dpo | log10k | 2.5 | 0.136 | 16.802 | 5.061 | 123.78 | 0 | **bound_limited** |
| k_arp_tdp | log10k | 0.939 | 40.699 | 65.225 | 42.486 | 1.60 | 0 | **quadratic** |

## Data-admissible bands of the unidentified coordinates (slice cost rise <= 3.59)

| coordinate | declared band | admissible band |
|---|---|---|
| k_glc_ha | [-10, 0.5] | [-10, 0.5] |
| k_dimer_fft | [-10, 0.5] | [-10, 0.4408] |
| ph_acid_yield_per_sink_event | [0, 1] | [0, 0.4189] |
| k_fur_fft | [-10, 0.5] | [-10, 0.5] |
| k_dimer_decay | [-10, 0.5] | [-10, -0.8679] |
| k_thiol_decay | [-10, 0.5] | [-10, -1.353] |
| Ea_decay_carbonyl_sink | [20, 250] | [91.9, 250] |
