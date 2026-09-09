# Wave B20 ship rule: SHIP

*Rule: SHIP if T1 (the prereg's decisive rows within 0.3 dex), T2 (no scored panel row moves 0.05 dex) and T5 (every coordinate identified, off its bound) hold; T3 and T4 reported. Pre-registration `results/validation/kinetic_core_b20_prereg.md`.*

Cost 17.21 on 10 rows; fitted constants (log10 at 100 C): {"log10_k_glyc_100C": -4.707, "log10_k_flp_cml_100C": -3.317, "log10_k_flp_cel_100C": -3.652, "log10_k_flp_decay_100C": -1.856, "log10_k_cml_loss_100C": -0.983}

| test | result | pass |
|---|---|---|
| T1 decisive rows | worst nguyen_k_glyc_120C -0.19 dex; all decisive {"nguyen_k_glyc_120C": -0.19, "nguyen_k_glyc_130C": 0.1, "nguyen_k_flp_cml_130C": 0.08, "nguyen_k_flp_decay_120C": 0.12, "nguyen_k_flp_decay_130C": -0.02, "nguyen_k_flp_cel_130C": 0.01, "nguyen_k_cml_loss_130C": 0.13}; reported {"nguyen_k_flp_cml_120C": -0.46, "nguyen_k_flp_cel_120C": -0.36, "nguyen_k_cml_loss_120C": -0.45} | True |
| T2 panel untouched | 39 predicted numbers compared; worst change 0.00e+00 dex; refused rows [25, 25] | True |
| T3 Nguyen's pot 120 C / 30 min | CML 0.0529 mmol/L (printed [0.025, 0.135]), CEL 0.0462, lysine lost 14.1 % | True |
| T4 other laboratories | Berk 2021 180 C (dry sesame) +1.73 dex; Hamzalioglu 2026 milk {"110C": 0.87, "120C": -0.07, "130C": 0.71, "140C": 1.17} dex; expanded soybean, 110 C, 60 min: about 25 % of the bound lysine lost (a moist solid with sucrose, which the engine cannot charge; direction only) | reported |
| T5 identification | sigma {"log10_k_glyc_100C": 0.085, "log10_k_flp_cml_100C": 0.175, "log10_k_flp_cel_100C": 0.112, "log10_k_flp_decay_100C": 0.152, "log10_k_cml_loss_100C": 0.253}; on bound False | True |
