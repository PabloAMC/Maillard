# Cheap Refinement Screening

| Reaction Family | Screen Decision | Improvement | Direction | Offsets | Mechanistic Benchmarks |
| --- | --- | ---: | --- | --- | --- |
| retro_aldol | advance | 2599.78 | down_3.0 | {"retro_aldol": -3.0, "retro_aldol_fragmentation": -3.0} | none |
| schiff_condensation | advance | 708.02 | up_3.0 | {"lipid_schiff_base": 3.0, "schiff": 3.0, "schiff_base_formation": 3.0, "schiff_condensation": 3.0} | none |
| thiol_addition | advance | 6.12 | up_3.0 | {"thiol_addition": 3.0, "thiol_addition_h2": 3.0, "thiol_addition_legacy_shortcut": 3.0} | none |
| strecker_degradation | advance | 0.32 | down_3.0 | {"strecker": -3.0, "strecker_degradation": -3.0} | none |
| thiol_oxidation | defer | 1650.41 | down_1.0 | {"thiol_oxidation": -1.0} | pea_isolate_ribose_cysteine_100C_45min_Internal2026, pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026, soy_isolate_ribose_cysteine_100C_45min_Internal2026, soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026 |
| amadori_rearrangement | defer | 359.76 | up_1.0 | {"amadori": 1.0, "amadori_rearrangement": 1.0} | pea_isolate_ribose_cysteine_100C_45min_Internal2026, pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026, soy_isolate_ribose_cysteine_100C_45min_Internal2026, soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026 |
| 1,2-enolisation | defer | 194.89 | up_1.0 | {"1,2-enolisation": 1.0, "1,2_enolisation": 1.0, "enol": 1.0, "enolisation": 1.0, "enolisation_1_2": 1.0} | pea_isolate_ribose_cysteine_100C_45min_Internal2026, pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026, soy_isolate_ribose_cysteine_100C_45min_Internal2026, soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026 |
| aminoketone_condensation | defer | -36.76 | up_1.0 | {"aminoketone_condensation": 1.0} | pea_isolate_ribose_cysteine_100C_45min_Internal2026, pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026, soy_isolate_ribose_cysteine_100C_45min_Internal2026, soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026 |
| cysteine_thermolysis | reject | 0.00 | down_1.0 | {"cys": -1.0, "cysteine_degradation": -1.0, "cysteine_thermolysis": -1.0} | none |

Advance candidates: 4
Deferred candidates: 4
Rejected candidates: 1
