# Objective Progress

| Objective | Closed / Target | Status | Resolved In Last Step | Prediction Effect |
| --- | ---: | --- | --- | --- |
| External mixed meaty-positive package | 0 / 2 | blocked_on_external_data | explicit_required_external_package_for_pea_and_soy, promotion_delta_if_package_lands_today | No promotion-ready lane is unlocked yet; the repo now states exactly which package would move the decision gate. |
| Hexanal/Nonanal ambiguity | 4 / 4 | closed_internal_calibration_route | prediction_validation_chain_exposed, closed_marker_counts_visible | All tracked Hexanal/Nonanal markers are within the accepted internal ratio band, so adverse-marker drift is calibration-closed without upgrading promotion posture. |
| Extrusion direct damage closure package | 0 / 3 | blocked_on_external_data | shared_extrusion_blocker_encoded_as_contract_artifact, dha_lysinoalanine_external_package_specified | Extrusion reporting now states the exact direct-marker package still missing for DHA/LAL and lysine-damage claims, and names the external measurement bundle required to move closure while keeping thiamine retention and soy thermal-history anchors explicitly partial. |

## Hexanal Nonanal Prediction Change

| Protein | Compound | Internal2026 ppb | ProtocolPilot2026 ppb | Ratio | Closure State |
| --- | --- | ---: | ---: | ---: | --- |
| pea_iso | Hexanal | 0.742533 | 0.742533 | 1.000 | calibration_closed |
| pea_iso | Nonanal | 0.0389915 | 0.0389915 | 1.000 | calibration_closed |
| soy_iso | Hexanal | 1.70062 | 1.70062 | 1.000 | calibration_closed |
| soy_iso | Nonanal | 0.0425151 | 0.0425151 | 1.000 | calibration_closed |

## Pea Soy Promotion Delta

| Protein | Mixed External Present Today | Would Close If Package Lands | Still Remaining After Package |
| --- | --- | --- | --- |
| pea_iso | no | comparator_is_measured_volatiles, external_quantitative_origin, minimum_quantitative_closed_targets | no_internal_or_directional_dependencies |
| soy_iso | no | comparator_is_measured_volatiles, external_quantitative_origin, minimum_quantitative_closed_targets | no_internal_or_directional_dependencies |

## Pea Soy Mixed External Package Delta

| Matrix | Status | Benchmark Candidate | Current Anchor | Desirable Targets | Adverse Targets | Would Close | Remaining After Package |
| --- | --- | --- | --- | --- | --- | --- | --- |
| pea_iso | specified_not_yet_measured | pea_isolate_ribose_cysteine_100C_45min_Internal2026 | pea_isolate_40C_PratapSingh2021 | 2-methyl-3-furanthiol, 2-furfurylthiol, bis(2-methyl-3-furyl) disulfide, 2,5-dimethylpyrazine | Hexanal, 2-pentylfuran, furfural, 1-octen-3-ol | external_quantitative_origin, comparator_is_measured_volatiles, minimum_quantitative_closed_targets | no_internal_or_directional_dependencies |
| soy_iso | specified_not_yet_measured | soy_isolate_ribose_cysteine_100C_45min_Internal2026 | soy_isolate_40C_PratapSingh2021 | 2-methyl-3-furanthiol, 2-furfurylthiol, bis(2-methyl-3-furyl) disulfide, 2,5-dimethylpyrazine | Hexanal, 2-pentylfuran, furfural, 1-octen-3-ol | external_quantitative_origin, comparator_is_measured_volatiles, minimum_quantitative_closed_targets | no_internal_or_directional_dependencies |

## Extrusion Direct Closure Delta

| Matrix | Direct Closure Ready | Missing Direct Markers | Supportive Anchors | Contextual Anchors | Remaining Requirements |
| --- | --- | --- | --- | --- | --- |
| pea_iso | no | Reactive lysine fraction, Lysinoalanine (LAL) | none | none | direct_crosslink_marker_external_quantified, lysine_damage_marker_external_quantified, traceable_source_metadata |
| soy_iso | no | Reactive lysine fraction, Lysinoalanine (LAL) | de_leyn_2019_thiamine_retention | troise_2018_soy_thermal_history | direct_crosslink_marker_external_quantified, lysine_damage_marker_external_quantified, traceable_source_metadata |

## DHA Lysinoalanine External Package Delta

| Matrix | Status | Direct Damage Targets | Paired Meaty Targets | Metadata Required | Would Close | Remaining After Package |
| --- | --- | --- | --- | --- | --- | --- |
| pea_iso | specified_not_yet_measured | Reactive lysine fraction, Lysinoalanine (LAL) | 2-methyl-3-furanthiol, 2-furfurylthiol | extruder_configuration, specific_mechanical_energy_kj_per_kg, feed_moisture_pct, die_exit_temperature_celsius, residence_time_seconds, post-process_sampling_delay_minutes | direct_crosslink_marker_external_quantified, lysine_damage_marker_external_quantified | promotion_review_with_same_run_tradeoff_panel |
| soy_iso | specified_not_yet_measured | Reactive lysine fraction, Lysinoalanine (LAL) | 2-methyl-3-furanthiol, 2-furfurylthiol | extruder_configuration, specific_mechanical_energy_kj_per_kg, feed_moisture_pct, die_exit_temperature_celsius, residence_time_seconds, post-process_sampling_delay_minutes | direct_crosslink_marker_external_quantified, lysine_damage_marker_external_quantified | promotion_review_with_same_run_tradeoff_panel |

Selected next family sprint: dha_lysinoalanine_external_package
Policy: objective_progress_must_state_what_changed_without_turning_internal_closure_into_external_promotion
