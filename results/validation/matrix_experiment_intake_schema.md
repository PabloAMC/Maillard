# Matrix Experiment Intake Schema

Contract id: matrix_experiment_intake_v1
Schema version: 1.0
Description: Machine-readable intake schema for comparing a new matrix experiment payload against the current benchmark surface without synthetic closure.

## Required Fields

Top-level: experiment_id, source_kind, protein_type, process_state, conditions, formulation, measured_volatiles, provenance
Conditions: temp_C, ph, water_activity, time_min
Formulation: precursors
Provenance: origin, source_reference, source_doi, measurement_date, notes

## Allowed Values

Source kinds: external_literature, internal_experiment, synthetic_diagnostic
Protein types: pea_iso, soy_iso, myco, free

## Policies

- Partial compound panels are allowed. Do not backfill unmeasured compounds.
- The comparison artifact is not itself a promotion decision.
- External decision readiness still depends on the matrix promotion contract and closure audit.
- evidence_class=external_validation_only means the payload is a hold-out surface: materialize it, evaluate it, but never absorb it into calibration or promotion counts.
