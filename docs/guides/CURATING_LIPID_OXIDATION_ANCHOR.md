# Curating a lipid-oxidation external hold-out anchor

> **Gate (hard rule): do not invent numbers.** This checklist describes *how* to add a
> real, literature-sourced lipid-oxidation anchor to the frozen external hold-out. Do not
> add an entry until a qualifying paper is in hand. The 4 existing bundles in
> `data/benchmarks/external_validation/` remain the frozen test set.

## Why this matters

S27 established that the matrix lipid-aldehyde over-prediction on the external hold-out is a
**per-(matrix, process_state) calibration gap**, not a kinetic-shape problem
(see [src/lipid_oxidation.py](../../src/lipid_oxidation.py) and
[data/lit/lipid_oxidation_calibration.json](../../data/lit/lipid_oxidation_calibration.json)).
Workstream B made the model *honestly uncertain* on uncalibrated process-states (external
hold-out 0/8 → 3/8 inside 90% CI), but the extreme-processing cases (HME extrusion, roasting)
remain genuine misses. **Closing them requires real measured anchors at those process-states**,
not a wider prior or a tuned cap. That is the curation task below.

## What qualifies as an anchor (intake criteria)

A candidate paper must report **all** of the following for a well-defined plant-protein matrix:

1. **Quantitative** headspace concentrations (ppb or convertible) for at least one lipid-oxidation
   marker: hexanal, nonanal, 2-pentylfuran, 1-hexanol, or (E,E)-2,4-decadienal.
2. **SIDA-grade or internal-standard-calibrated** quantification (HS-SPME-GC-MS/SIM acceptable;
   semi-quant peak areas are **not** acceptable for a hold-out anchor).
3. A **defined matrix** (protein type + isolate/concentrate) and **process-state** parameters:
   temperature, time, and (ideally) water activity / moisture regime / SME for extrusion.
4. A process-state that **fills a calibration gap** — i.e. roasting, HME extrusion, or another
   high-severity regime not already pinned by the calibration registry. (Ambient slurry is
   already represented.)
5. Enough method detail to set `analytical_context` (method, quantification_mode, replicates).

If any of 1–3 is missing, the paper is **not** hold-out-grade — log it in the backlog instead.

## Steps to land an anchor (once a qualifying paper is in hand)

1. **Extract** the measured marker concentrations and the full process-state, in the paper's own
   units. Record the DOI. Convert to ppb explicitly; show the conversion.
2. **Add a bundle spec** to `_HOLDOUT_BUNDLE_SPECS` in
   [src/external_validation.py](../../src/external_validation.py), mirroring the existing entries
   (`bundle_id`, `anchor_ids`, `matrix_context`, `protein_type`, `process_state`, `conditions`,
   `precursors`, `benchmark_alignment`, `analytical_context`). Keep
   `evidence_class = external_validation_only` so it stays out of calibration.
3. **Generate** the benchmark JSON under `data/benchmarks/external_validation/` via the existing
   payload generator (`scripts/generators/generate_external_validation_payloads.py`).
4. **Regenerate** the report: `./scripts/docker_maillard.sh run "python scripts/generators/generate_external_validation_report.py"`.
   The new anchor will appear with its honest fold-error and CI status — **report it as-is**.
5. **Do not** retune `max_conversion_fraction`, the observable sigmas, or any prior to make the new
   anchor pass. Anchors are a test set, not a fit target. If the new point is a large miss, that is
   a true signal that the process-state needs *dedicated calibration data*, which is a separate,
   in-panel ingestion task (`./scripts/docker_maillard.sh ingest`).

## If the paper is good but not hold-out-grade

Log it in the literature backlog (`results/validation/literature_backlog.*`) with the reason it
was rejected (e.g. semi-quant only, undefined matrix). It may still inform the in-panel
calibration once a quantitative companion is found.

## Related

- Mechanism + S27 findings: [docs/architecture.md](../architecture.md) (matrix & headspace layer)
- Validation contract: [docs/reference/VALIDATION_CONTRACT.md](../reference/VALIDATION_CONTRACT.md)
- Diagnostic: `scripts/diagnose_lipid_bias.py` (measured vs predicted lipid markers, in-panel + hold-out)
