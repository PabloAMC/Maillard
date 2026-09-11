# Wave B29 fit report

Cost 20.804 on 4 ratios, 2 free. Lever on k_ama_g, k_glc_g; observable AKG + AKM (the total Strecker flux, not one dicarbonyl's).

Fitted factors: {"air": 1.0, "argon": 0.001, "air_cu": 99.944}

| ratio | printed | model | residual (dex) |
|---|---:|---:|---:|
| arp_air_over_argon | 9.2 | 2.188 | -0.624 |
| glc_air_over_argon | 3.5 | 2.756 | -0.104 |
| arp_aircu_over_air | 2.5 | 1.387 | -0.256 |
| glc_aircu_over_air | 1.9 | 2.142 | +0.052 |

Laplace: sigma {"log10_f_argon": 42.085, "log10_f_air_cu": 18.967}; identified {'log10_f_argon': False, 'log10_f_air_cu': False}; on bound {'log10_f_argon': True, 'log10_f_air_cu': True}.

Air is exactly 1: parameters identical True, observable identical True.
