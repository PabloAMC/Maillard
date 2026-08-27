# Glossary for Scientists

This document translates the computational framework's internal terminology into standard
food science concepts.

It is split in two. **Part 1** defines the labels you will literally see in `report.md`,
`report.json` and the validation artifacts. **Part 2** defines concepts used in this
documentation and in conversation, which are *not* emitted as labels anywhere.

---

## Part 1 — Labels emitted in output

### The three tier vocabularies

Three different tier scales coexist in one report. They are not interchangeable, and
confusing them is the single most common misreading of the output.

| Field | Values | Grades what? |
|---|---|---|
| `tier` | `high`, `medium`, `low`, `exploratory` | How well **this run** is benchmark-supported |
| `calibration_evidence_strength` (shown as *Evidence*) | `literature_anchored`, `conditional_literature_anchored`, `class_anchored`, `directional_transferred`, `process_state_mismatch`, `heuristic` | What kind of **anchor** stands behind one compound's projection |
| `confidence_tier` | `high`, `medium_high`, `medium`, `medium_low`, `low` | How strong a curated **literature prior** in `data/lit/` is — a property of the source, not of your prediction |

A `high` `confidence_tier` prior can perfectly well feed an `exploratory` `tier` prediction,
because your formulation may sit far from the conditions that prior was measured at.

| Term | What it means in practice |
|---|---|
| **`tier`** | Band on a 0-100 confidence score for the current run: `high` (>= 85), `medium` (65-84), `low` (45-64), `exploratory` (< 45). Appears at run level, per compound, and per sensory tag. |
| **`prediction_mode`** | The action each `tier` licenses, one-to-one: `benchmark_supported_quantitative` (high), `ranking_supported` (medium), `directional_only` (low), `hypothesis_only` (exploratory). |
| **`decision_mode`** | Coarser roll-up of `prediction_mode`: `quantitative_recommendation` only when the run is `benchmark_supported_quantitative`, otherwise `directional_hypothesis`. |
| **`calibration_evidence_strength`** | Strength of the anchor behind a compound: `literature_anchored` (published measurement on a comparable system) > `conditional_literature_anchored` (literature-anchored under stated conditions) > `class_anchored` (anchored at compound-class level) > `directional_transferred` (prior transferred from an adjacent matrix/process) > `process_state_mismatch` (anchor exists but for a different process state; nearest was substituted) > `heuristic` (no anchor; built-in class default). |
| **`scope_demoted_from`** | Present in `report.json` when a run is out of calibration scope. Every evidence strength is demoted one notch and this records the pre-demotion value. |
| **`evidence_state`** | Where a compound sits on the evidence ladder: `externally_benchmarked`, `conditional_calibration`, `transferred_prior`, `still_missing`. |
| **`reachability_status`** | `chemically_reachable` — an enumerated, barrier-scored pathway produces the compound from your precursors. `conditionally_reachable` — reachable only under an assumption the run had to make. `merely_plausible` — no enumerated route; the number is a class-level projection. |
| **`reachability_basis`** | Why reachability was granted: `direct_anchor`, `transferred_or_refined_support`, `mechanistic_surrogate_only`. |
| **`observable_assumption_summary`** | Pipe-joined triple `retention_runtime_mode \| calibration_fallback_mode \| support_origin`, e.g. `static_class_profile \| class_level \| standard_matrix_support`. |
| **`retention_runtime_mode`** | How the volatile's matrix retention was computed: `static_class_profile` (fixed per compound class), `direct_binding_plus_ph_release_reference` (measured binding plus pH-dependent release), `sulfur_binding_prior` / `pyrazine_binding_prior` (family-specific priors). |
| **`calibration_fallback_mode`** | `class_level` (fell back to the compound class) or `nearest_process_state` (used the closest calibrated process state rather than yours). |
| **`support_origin`** | Which support route the matrix correction took: `standard_matrix_support`, `lower_regime_transfer`, `extrusion_specific_support`, `extrusion_extrapolation`. |
| **`process_state`** | Coarse process bucket the run was assigned (e.g. `aqueous_pre_extrusion_model`, `extrusion_structured`), used to select calibrations. |
| **`process_regime`** | Severity bucket derived from temperature and water activity: `free_aqueous`, `matrix_hydrated`, `extrusion_like`, `extrusion_heavy`. The last two force a `prediction_mode` downgrade. |
| **`benchmark_neighborhood`** | Which benchmark family the run is nearest to: `primary_free_precursor`, `free_precursor_partial_analogy`, `matrix_intake_only`, `exploratory_matrix`, `sparse_precursor_analogy`. |
| **`target_class`** | What role a compound plays in the panel: `sulfur_positives`, `adverse_lipid_markers`, `pyrazines`, `severity_markers`, `pretreatment_state_markers`. |
| **`panel_role`** | How a row is used in validation: `scored` (counts toward the headline), `constrained` (bounded but not scored), `diagnostic`, `report_only`. |
| **`skipped_formulations`** | Formulations dropped before scoring because a precursor name could not be resolved. Listed with the offending name; they are absent from every ranking in that run. |
| **`input_fingerprint_sha256`** | SHA-256 over the *scientific* inputs of a run. Deliberately excludes `--output-dir`, `--report` and argv, so the same science fingerprints identically wherever it is written. |
| **Benchmark-validated** (or "strict-ready") | The model prediction is anchored to a specific, quantified, peer-reviewed literature experiment. Emitted as `strict_ready` in the validation artifacts. |
| **FAST proxy / FAST mode** | The primary, laptop-speed prediction mode. It uses heuristically calibrated energy barriers based on literature rather than running heavy, hours-long Quantum Mechanics calculations for every molecule. |
| **Validated Envelope** (or "Validated Domain") | The physical boundary (pH, temp, time, ingredients) where the model has literature backing. Going outside this envelope means the output is a hypothesis, not a proven fact. Emitted as `validated_envelope.{md,json}`. |
| **VoI (Value of Information)** | A priority score ranking which wet-lab experiment is most valuable to run next. It targets compounds where the model is highly uncertain (wide prediction bounds) *and* which have a high sensory impact (far above odour detection threshold). |
| **ODT (Odour Detection Threshold)** | The lowest concentration of a volatile compound that can be detected by the human nose (typically µg/kg of water/matrix). Used to calculate how decision-relevant a compound's yield is. |
| **Order of Magnitude (or Dex)** | Logarithmic steps representing a ten-fold ($10\times$) difference. A change of 1 dex is a $10\times$ difference in concentration; 2 dex is $100\times$. |
| **p5 / p95 (Confidence Whiskers)** | The statistical bounds of a prediction. There is a 90% probability that the actual concentration falls between the lower bound ($p_5$) and the upper bound ($p_{95}$). Rendered inline as `12.4 ppb [3.1-49.7, 90% CI]`. |
| **Melanoidin Trapping & Retention** | Physical or chemical trapping of volatile aroma molecules (especially sulfur thiols) by complex brown polymers (melanoidins) formed during late-stage Maillard browning. |
| **SIDA (Stable Isotope Dilution Assay)** | A highly precise mass spectrometry method using isotopically labeled standards to measure absolute concentration (ppb) in complex plant matrices. The preferred method for generating new benchmark data; appears in benchmark source metadata. |

---

## Part 2 — Concepts (not emitted as labels)

These are used in the documentation and in discussion, but you will not find them as field
values in any artifact.

| Term | What it means in practice |
|---|---|
| **Observable Projection** | The conversion layer that takes theoretical chemical yields in a beaker and corrects them for what a sensory panel would actually experience (pH-dependent headspace release, matrix trapping, degradation). In output this layer's fingerprints are the `projection_metadata` block and the `observable_*` fields above. |
| **Target-ranking Validation** | A benchmark where the model must predict rank order, not just presence (e.g. hexanal far more abundant than an oxidized lipid note). Its results surface in the benchmark artifacts as ranking assertions, not as a label. |
| **Intake Check** | A basic test that the software runs and produces the right compounds, without guaranteeing absolute predicted concentrations. Corresponds to the matrix intake contract, whose runs appear under `benchmark_neighborhood: matrix_intake_only`. |
| **Matrix Retargeting / Matrix Calibration** | Adjusting theoretical free-system predictions (water/buffer) downward to account for physical trapping inside structured plant proteins such as pea or soy isolates. The mechanism is visible through `calibration_source`, `support_origin` and `calibration_fallback_mode`, not under this name. |
