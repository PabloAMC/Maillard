# Hexanal Nonanal Calibration Closure

| Lane | Protein | Primary Benchmark | Reference Benchmark | Closure Action | Hexanal Ratio | Hexanal In Bounds | Nonanal Ratio | Nonanal In Bounds | Next Best Action |
| --- | --- | --- | --- | --- | ---: | --- | ---: | --- | --- |
| pea_protocol_pilot_vs_internal2026 | pea_iso | pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026 | pea_isolate_ribose_cysteine_100C_45min_Internal2026 | calibration_closed | 1.000 | yes | 1.000 | yes | retain_internal_calibration_route_and_seek_external_quantitative_mixed_matrix_evidence |
| soy_protocol_pilot_vs_internal2026 | soy_iso | soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026 | soy_isolate_ribose_cysteine_100C_45min_Internal2026 | calibration_closed | 1.000 | yes | 1.000 | yes | retain_internal_calibration_route_and_seek_external_quantitative_mixed_matrix_evidence |

Lanes checked: 2
Closed lanes: 2
Hazard lanes: 0
Closed markers: 4 / 4
Default ratio bounds: [0.5, 2.0]
Policy: hexanal_nonanal_protocol_pilot_to_internal2026_ratio_check_closes_internal_calibration_routes_but_not_external_promotion_claims

## Prediction Validation Chain

| Lane | Protein | Compound | Internal2026 ppb | ProtocolPilot2026 ppb | Ratio | Closure State | Next Decision Gate |
| --- | --- | --- | ---: | ---: | ---: | --- | --- |
| pea_protocol_pilot_vs_internal2026 | pea_iso | Hexanal | 0.742533 | 0.742533 | 1.000 | calibration_closed | retain_internal_calibration_route_and_seek_external_quantitative_mixed_matrix_evidence |
| pea_protocol_pilot_vs_internal2026 | pea_iso | Nonanal | 0.0389915 | 0.0389915 | 1.000 | calibration_closed | retain_internal_calibration_route_and_seek_external_quantitative_mixed_matrix_evidence |
| soy_protocol_pilot_vs_internal2026 | soy_iso | Hexanal | 1.70062 | 1.70062 | 1.000 | calibration_closed | retain_internal_calibration_route_and_seek_external_quantitative_mixed_matrix_evidence |
| soy_protocol_pilot_vs_internal2026 | soy_iso | Nonanal | 0.0425151 | 0.0425151 | 1.000 | calibration_closed | retain_internal_calibration_route_and_seek_external_quantitative_mixed_matrix_evidence |
