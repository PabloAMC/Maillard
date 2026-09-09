# Wave B24 ship rule: DO NOT SHIP

*Rule: SHIP if T1 (fed rows within 0.3 dex, proline rows within 0.5 and in order), T2 (no panel row moves 0.05 dex) and T5 (both identified, off their bounds) hold; T3 and T4 reported. Pre-registration `results/validation/kinetic_core_b24_prereg.md`.*

Cost 53.69 on 5 rows; optimum {"log10_k_pyrl_ap_100C": -2.562, "log10_k_mgo_pro_100C": -6.621}

| test | result | pass |
|---|---|---|
| T1 rows | fed 1-pyrroline {"hof_t7_e1_pyrl2_mgo10": 0.18, "hof_t7_e2_pyrl2_mgo2": 0.3}; proline {"hof_t9_pro400_mgo4": -1.29, "hof_t9_pro400_mgo40": 0.13, "hof_t9_pro400_mgo400": 1.16}; order True | False |
| T2 panel | 39 numbers, worst change 0.00e+00 dex | True |
| T3 apparent barrier | model 591 kJ/mol vs Chan & Reineccius 1994's 60.2 | reported |
| T4 excess pyrroline | model 41.2 mol % of methylglyoxal vs printed 0.33 (+2.10 dex) | reported |
| T5 identification | sigma {"log10_k_pyrl_ap_100C": 0.53, "log10_k_mgo_pro_100C": 0.66}; on bound False | True |
