# Wave B22b fit report

Cost 5.671 on 5 rows, 2 free. Model peaks at 30 min; the source peaks at 120.

| minutes | printed umol/L | model | residual (dex) |
|---:|---:|---:|---:|
| 30 | 0.244 | 0.7811 | +0.505 |
| 60 | 0.52 | 0.7811 | +0.177 |
| 120 | 1.887 | 0.7811 | -0.383 |
| 180 | 1.477 | 0.7811 | -0.277 |
| 240 | 0.822 | 0.7811 | -0.022 |

Laplace: sigma {"log10_k_marp_mtal_120C": 0.092, "log10_k_marp_loss_120C": 0.092}; identified {'log10_k_marp_mtal_120C': True, 'log10_k_marp_loss_120C': True}; on bound {'log10_k_marp_mtal_120C': False, 'log10_k_marp_loss_120C': False}.

## The two-arm check, not fitted

| minutes | printed ARP / binary | model ARP / binary |
|---:|---:|---:|
| 30 | 1.37 | 7.81e+29 |
| 60 | 2.14 | 7.81e+29 |
| 120 | 2.6 | 7.81e+29 |
| 180 | 1.75 | 7.81e+29 |
| 240 | 1.1 | 7.81e+29 |

> The binary arm rests on charging methionine as GLYCINE for the Amadori chemistry, a declared substitution B22 already found overshoots Deng's pot by 1.5 to 2.9 decades. This check is where that shows.
