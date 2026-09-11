# Wave B21 ship rule: SHIP

*Rule: SHIP if T1 (decisive rows within 0.3 dex), T2 (browning hold-out median below 2 and within 3x, no panel row moves 0.3 dex), T3 (Quan's glyoxal within the widened range) and T6 (both coordinates identified) hold; T4 and T5 reported. Pre-registration `results/validation/kinetic_core_b21_prereg.md`.*

Cost 3.54 on 6 rows; fitted (log10 at 100 C): {"log10_k_ama_g_100C": -2.034, "log10_k_g_go_aqueous_100C": -0.51}

| test | result | pass |
|---|---|---|
| T1 decisive rows | worst ham_k_ama_g_120C +0.18 dex; {"ham_k_ama_g_110C": -0.04, "ham_k_ama_g_120C": 0.18, "ham_k_ama_g_140C": -0.18, "ham_k_g_go_120C": 0.0}; reported {"ham_k_ama_g_130C": 0.36, "ham_k_g_go_130C": -0.01} | True |
| T2a browning hold-out | median fold 1.43 -> 1.31; within 3x 1.00 -> 1.00 | True |
| T2b panel | 39 predicted numbers; worst change 0.000 dex; moved over 0.01 dex: 0 | True |
| T3 Quan 2020 glyoxal | {"100C": {"model_mmol_l": 0.0234, "printed": [0.052, 0.127], "dex_from_range": -0.35}, "130C": {"model_mmol_l": 0.3871, "printed": [0.144, 0.605], "dex_from_range": 0.0}} | True |
| T4 Xia 2022 ordering | glyoxal 1.23 vs methylglyoxal 7.82 mmol/L at 130 C / 80 min; glyoxal above: False | reported |
| T5 Leahy total pyrazine | model 21.9 ug/L vs 13100 (-2.78 dex; B18 -2.88) | reported |
| T6 identification | sigma {"log10_k_ama_g_100C": 0.085, "log10_k_g_go_aqueous_100C": 0.144}; on bound False | True |
