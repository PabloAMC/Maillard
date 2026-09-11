# Wave B24b ship rule: DO NOT SHIP

*Rule: SHIP if T1 (the switch), T2 (the pyrroline excess), T3 (B24's fed rows) and T5 hold. Pre-registration `results/validation/kinetic_core_b24b_prereg.md`.*

| test | result | pass |
|---|---|---|
| T1 the switch | ordering increasing True; worst switch_mgo400 +1.01 dex; model spans 2.1x against a printed 80x | False |
| T2 pyrroline excess | 32.5 mol % against 0.33, +1.99 dex | False |
| T3 B24's fed rows | worst +0.164 dex | True |
| T4 identification | sigma {"log10_k_ha_athp_100C": 0.94, "log10_k_pyrl_loss_100C": 7.95}; on bound {'log10_k_ha_athp_100C': False, 'log10_k_pyrl_loss_100C': True} | reported |
| T4b the pH transfer | worst +0.99 dex | reported |
| T5 nothing else moves | panel identical True | True |
