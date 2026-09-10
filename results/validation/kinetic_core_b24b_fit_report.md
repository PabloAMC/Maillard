# Wave B24b fit report

Cost 79.463 on 4 rows, 2 free coordinates. B24's two constants held, not refitted.

| row | target | predicted | residual (dex) |
|---|---:|---:|---:|
| switch_mgo4 | 0.16 | 0.5964 | +0.571 |
| switch_mgo40 | 0.51 | 0.9159 | +0.254 |
| switch_mgo400 | 12.8 | 1.248 | -1.011 |
| pyrl_excess_ap_molpct | 0.33 | 32.46 | +1.993 |

Laplace: sigma {"log10_k_ha_athp_100C": 0.944, "log10_k_pyrl_loss_100C": 7.951}; identified {'log10_k_ha_athp_100C': True, 'log10_k_pyrl_loss_100C': False}; on bound {'log10_k_ha_athp_100C': False, 'log10_k_pyrl_loss_100C': True}.

## Checks, not fitted

| check | printed | model | dex |
|---|---:|---:|---:|
| hof_t7_e1_pyrl2_mgo10 | 28.7 mol % | 32.26 | +0.051 |
| hof_t7_e2_pyrl2_mgo2 | 5.3 mol % | 7.726 | +0.164 |
| ATHP ratio to pH 7 at pH 5.0 | 0.08333 | 0.8198 | +0.993 |
| ATHP ratio to pH 7 at pH 7.0 | 1 | 1 | +0.000 |
| ATHP ratio to pH 7 at pH 9.0 | 3.556 | 1.015 | -0.545 |

> Neither block is in the objective. The fed rows check that adding the branch did not break what B24 already fitted; the ladder checks B18's pH slopes on a chemistry they were not fitted on.
