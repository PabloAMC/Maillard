# Selective DFT Plan

| Reaction Family | DFT Decision | Tier | Cheap Improvement | Offsets |
| --- | --- | --- | ---: | --- |
| thiol_oxidation | defer | cheap_first_hold | 1650.41 | {"thiol_oxidation": -1.0} |
| amadori_rearrangement | defer | cheap_first_hold | 359.76 | {"amadori": 1.0, "amadori_rearrangement": 1.0} |
| 1,2-enolisation | defer | cheap_first_hold | 194.89 | {"1,2-enolisation": 1.0, "1,2_enolisation": 1.0, "enol": 1.0, "enolisation": 1.0, "enolisation_1_2": 1.0} |
| aminoketone_condensation | defer | cheap_first_hold | -36.76 | {"aminoketone_condensation": 1.0} |
| retro_aldol | reject | not_applicable | 2599.78 | {"retro_aldol": -3.0, "retro_aldol_fragmentation": -3.0} |
| schiff_condensation | reject | not_applicable | 708.02 | {"lipid_schiff_base": 3.0, "schiff": 3.0, "schiff_base_formation": 3.0, "schiff_condensation": 3.0} |
| thiol_addition | reject | not_applicable | 6.12 | {"thiol_addition": 3.0, "thiol_addition_h2": 3.0, "thiol_addition_legacy_shortcut": 3.0} |
| strecker_degradation | reject | not_applicable | 0.32 | {"strecker": -3.0, "strecker_degradation": -3.0} |
| cysteine_thermolysis | reject | not_applicable | 0.00 | {"cys": -1.0, "cysteine_degradation": -1.0, "cysteine_thermolysis": -1.0} |

Run-now jobs: 0
Deferred jobs: 4
Rejected jobs: 5
