# Matrix Experiment Ingestion

This guide explains how to compare a new matrix measurement set against the current model without turning the result into synthetic closure.

## What This Workflow Is For

Use this when you have a new matrix measurement set from:

- external literature
- an internal wet-lab experiment
- a synthetic diagnostic payload that should be compared, but not promoted as evidence

The workflow answers five practical questions:

1. Which compounds matched the current model?
2. Which compounds are externally closed, internally supported, only transferred, or still open?
3. Did the new measurement strengthen benchmark readiness?
4. What blocker still prevents external decision readiness?
5. Should the new data land in a benchmark payload, calibration payload, blocker registry, or diagnostic-only surface?

## Scope Reminder

Soy and pea are not the only matrix families tracked by the repository.

Current matrix-family scope includes:

- pea isolate: directional matrix support with executable benchmarks
- soy isolate: directional matrix support with executable benchmarks
- soy hydrolysate: qualitative intake support
- mycoprotein: bounded directional priors
- extrusion-heavy systems: exploratory process regime
- coconut oil co-matrix: indirect generic support only
- other plant proteins: explicit open gap until runtime-facing evidence exists

See results/validation/matrix_family_coverage.md for the current scope boundary.

## Intake Contract

The canonical machine-readable schema is:

- data/protocols/matrix_experiment_intake_schema.json

An example payload is provided at:

- data/protocols/example_matrix_experiment_intake.yaml

Required top-level fields:

- experiment_id
- source_kind
- protein_type
- process_state
- conditions
- formulation
- measured_volatiles
- provenance

Partial compound panels are allowed. Do not backfill unmeasured compounds.

## Run The Comparison

Generate the schema artifact:

```bash
./scripts/docker_maillard.sh experiment-intake-schema
```

Compare an intake payload against the model and the closest aligned matrix benchmark:

```bash
./scripts/docker_maillard.sh compare-experiment data/protocols/example_matrix_experiment_intake.yaml
```

This generates support-delta artifacts in results/validation.

## Key Outputs

- results/validation/matrix_experiment_intake_schema.md
- results/validation/example_matrix_experiment_intake_support_delta.md
- results/validation/example_matrix_experiment_intake_support_delta.json

## How To Interpret The Result

Support statuses:

- quantitative_closed: externally measured and benchmark-usable for that compound
- internal_candidate: useful for reproducibility or mechanistic triage, not release-grade promotion
- directional_support: transferred or class-level support, useful for prioritization but still bounded
- open_gap: no usable support yet

Support deltas:

- strengthened: the new data improved support relative to the aligned baseline
- unchanged: no support-level change
- weakened: the new data reduced support or exposed a contradiction
- new_compound: the new data added a previously untracked measured target

Readiness changes:

- promoted_to_external_decision_ready
- evidence_strengthened_not_yet_promotable
- promotion_blocker_shifted
- no_material_change

## Promotion Discipline

The comparison artifact is not itself the promotion decision.

Promotion still depends on the matrix promotion contract in:

- results/validation/matrix_promotion_contract.md

And the compound-by-compound closure audit in:

- results/validation/matrix_observable_closure_audit.md

If a new measurement strengthens support but still leaves transferred or internal-only dependencies, keep the lane below external decision readiness and encode the blocker explicitly.