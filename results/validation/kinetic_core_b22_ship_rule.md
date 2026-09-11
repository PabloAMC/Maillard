# Wave B22 ship rule: DO NOT SHIP

*Rule: SHIP if T1 (methional and methanethiol rows within 0.3 dex, disulfide within 0.5, release barrier off its bounds), T2 (no panel row moves 0.05 dex) and T4 (every coordinate identified, the ratio inside its band) hold; T3 and T5 reported. Pre-registration `results/validation/kinetic_core_b22_prereg.md`.*

Cost 7374 on 9 rows; optimum {"log10_identity_ratio_met_over_gly": 2.0, "log10_k_mtal_msh_100C": -0.281, "ea_mtal_msh_kj_mol": 20.0, "log10_k_msh_dmds_100C": 2.0}

| test | result | pass |
|---|---|---|
| T1 rows | methional and methanethiol {"pan_MTAL_rate_100C": -5.57, "pan_MTAL_rate_120C": -4.53, "pan_MTAL_rate_140C": -3.64, "pan_MSH_rate_100C": -5.32, "pan_MSH_rate_120C": -3.37, "pan_MSH_rate_140C": -1.99}; disulfide {"pan_DMDS_rate_100C": -10.56, "pan_DMDS_rate_120C": -6.43, "pan_DMDS_rate_140C": -2.64}; release barrier on bound True | False |
| T2 panel | 39 numbers, worst change 0.00e+00 dex | True |
| T3 Deng 2022 | methional at 120 min +1.93 dex; rising 30 -> 120 min False; all {"30min": 2.93, "60min": 2.72, "120min": 1.93, "180min": 1.55} | reported (False) |
| T4 identification | sigma {"log10_identity_ratio_met_over_gly": 2.5, "log10_k_mtal_msh_100C": 6.16, "ea_mtal_msh_kj_mol": 643.44, "log10_k_msh_dmds_100C": 8.77}; on bound {"log10_identity_ratio_met_over_gly": true, "log10_k_mtal_msh_100C": false, "ea_mtal_msh_kj_mol": true, "log10_k_msh_dmds_100C": true}; ratio inside band False | False |
| T5 Chin & Lindsay 1994 | methanethiol half-life at 30 C: model 69.2 min vs 17.0 with copper | reported |
