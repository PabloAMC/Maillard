# Wave B29 fit report

Cost 5.597 on 4 ratios, 2 free. Lever on k_ama_g, k_glc_g; observable AKG.

Fitted factors: {"air": 1.0, "argon": 0.0494, "air_cu": 14.6668}

| ratio | printed | model | residual (dex) |
|---|---:|---:|---:|
| arp_air_over_argon | 9.2 | 5.63 | -0.213 |
| glc_air_over_argon | 3.5 | 5.804 | +0.220 |
| arp_aircu_over_air | 2.5 | 1.729 | -0.160 |
| glc_aircu_over_air | 1.9 | 2.291 | +0.081 |

Laplace: sigma {"log10_f_argon": 0.248, "log10_f_air_cu": 1.615}; identified {'log10_f_argon': True, 'log10_f_air_cu': False}; on bound {'log10_f_argon': False, 'log10_f_air_cu': False}.

Air is exactly 1: parameters identical True, observable identical True.
