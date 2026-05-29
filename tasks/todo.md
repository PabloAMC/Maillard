# Maillard Strategic Roadmap

> **Backup of the previous long-form roadmap**: [tasks/.todo.md.bak](.todo.md.bak). Closed sprints have full review blocks preserved there and in git history; this file is the live working surface only.

## Priority order (set 2026-05-23 after literature-first realignment)

> **Operating constraint (user, 2026-05-13):** *"Do not complicate the tool more than it can do really really well."* Every new lane below has a stop-the-line gate that prevents adding a feature that the kernel cannot back with calibrated evidence.

1. **ACTIVE SPRINT: S25 — Phased & QC-Gated Literature Ingestion** — Systematically ingest the 99 backlog references from `deep_research_backlog.json` to ground our kinetics and matrices in empirical data. Ingested in 3 family-based chunks with strict QC validation steps to avoid poor quality data.
2. **PARKED — Targeted Literature Search for Unresolved Gaps** — Cast a targeted net for missing matrix benchmarks (`ppi_meaty_positive_matrix_benchmark`, `spi_meaty_positive_matrix_benchmark`, and `meaty_off_flavour_safety_tradeoff_panel`) after Phase 1 backlog ingestion completes.
3. **PARKED — Triangulated Barrier Strategy for Families 11–15** (gated by observables, parked pending literature ingestion).
4. **PARKED — pyGSM integration & global barrier run** (gated by observables, parked).
5. **PARKED — `pe_amadori` multi-step DFT closure** (parked, heavy compute).
6. **MAINTENANCE BACKLOG** — `2026-04-17 DFT Ladder Bulletproofing` (3/4 items remaining) + `2026-04-13 Architecture Hardening`.
7. **STRATEGIC BACKLOG** — `S14 monolith decomposition`, `S19 web interface`. Defer.
8. **DONE / SUPERSEDED** — see [Recent closures](#recent-closures).

---

## Active Sprint — 2026-05-23 S25 — Phased & QC-Gated Literature Ingestion

**Branch (proposed)**: `s25-literature-ingestion`. **TL;DR**: Ground the model in empirical data by ingesting the 99 backlog references from `deep_research_backlog.json` in three QC-gated chunks.

### QC Validation checklist per reference (The "Gate")
1. **Extraction Validation**: Extract temperature, duration, pH, and matrix details. Confirm quantification method (GC-MS, SAFE, SIDA) and unit conversions (ppb/ppm). Verify physical constraints.
2. **Registry Mapping**: Encode under [data/lit/benchmark_intake_registry.json](data/lit/benchmark_intake_registry.json). Ensure unique IDs and non-colliding matrix families.
3. **Payload Generation**: Create/update corresponding computational priors or payloads under `data/lit/`.
4. **Verification**: Run `pytest tests/unit -q` and `python -m src.literature_learning_loop` inside Docker. Ensure no calibration inflation.

### Lane S25-A: Phase 1 — Chunk 1 Ingestion: Core & Matrix Basics (27 items)
- [x] A.1 Families `01` (Core - completed), `02` (Lipid Oxidation - completed), `06` (Alternative Proteins - completed), `07` (Reducing Sugars - completed).
- [x] A.2 Extract, QC-validate, and encode Chunk 1 references into `benchmark_intake_registry.json` (27/27 references completed).
- [x] A.3 Generate computational priors and payloads for Chunk 1 (27/27 references completed).
- [x] A.4 Validate Chunk 1 via Docker tests and literature learning loop (Families 01, 02, 06, 07 validated).

### Lane S25-B: Phase 2 — Chunk 2 Ingestion: Flavor, Degradation & Fermentation (33 items)
- [x] B.1 Families `08` (Off-notes), `09` (Carbohydrate degradation), `10` (Fermentation), `11` (Lipid-Maillard crosstalk).
- [x] B.2 Extract, QC-validate, and encode Chunk 2 references into `benchmark_intake_registry.json`.
- [x] B.3 Generate computational priors and payloads for Chunk 2.
- [x] B.4 Validate Chunk 2 via Docker tests and literature learning loop.

### Lane S25-C: Phase 3 — Chunk 3 Ingestion: Advanced Caps, Damage & Polymers (41 items)
- [x] C.1 Families `03` (Thiamine), `04` (Nucleotides), `05` (Glutathione), `12` (Protein damage), `13` (Polyphenols), `14` (Ascorbic acid), `15` (Phospholipids), `16` (Melanoidins).
- [x] C.2 Extract, QC-validate, and encode Chunk 3 references into `benchmark_intake_registry.json`.
- [x] C.3 Generate computational priors and payloads for Chunk 3.
- [x] C.4 Validate Chunk 3 via Docker tests and literature learning loop.

### Lane S25-D: Phase 4 — Targeted Literature Search for Key Matrix Gaps
- [x] D.1 Search literature databases (OpenAlex, Europe PMC, and web engines) for `ppi_meaty_positive_matrix_benchmark`, `spi_meaty_positive_matrix_benchmark`, and `meaty_off_flavour_safety_tradeoff_panel` (completed; confirmed no single literature package closes these combined matrix-aroma-safety contracts).
- [x] D.2 QC-validate and ingest matching references (none found to ingest; confirmed as genuine structural literature gaps requiring wet-lab measurements).

### Lane S25-E: Phase 5 — Re-evaluation & Roadmap Review
- [x] E.1 Review updated `literature_learning_loop.md` and `family_deviation_audit.md` (completed; verified error tails are within expected tolerances, and literature ingestion backlog has been successfully processed).
- [x] E.2 Decide next strategic steps (completed; determined that literature grounding is complete, and the next priority is transitioning to primary matrix and extrusion wet-lab data collection as outlined in `structural_unlock_triage.py`).

---

## Parked — 2026-05-12 `pe_amadori` multi-step DFT closure

> **STATUS: PARKED behind S22.** The multi-step pipeline ([src/elementary_step_runner.py](src/elementary_step_runner.py) + [src/microsolvation.py](src/microsolvation.py)) is in place; the only configured target (`pe_amadori_water_shuttle`) crashes at the GSM call with `pyGSM exception: IndexError: index 34 is out of bounds for axis 0 with size 33` (see [results/computational_gap_refinement/multistep/pe_amadori_water_shuttle/target_result.json](results/computational_gap_refinement/multistep/pe_amadori_water_shuttle/target_result.json)). When unparked: M1 fixes the IndexError; M2 encodes the real two-step mechanism (O0→C2 then C1→N3); M3 writes back uncertainty narrowing only (no anchor promotion). Hard timebox on M2 (~24h DFT wall per step).

---

## Parked — 2026-04-22 Triangulated Barrier Strategy for Families 11–15

> **STATUS: PARKED pending literature ingestion.** Reactivation criterion: completion of Sprint 25 literature grounding and re-evaluation.

**Targets (7):** `hexanal_radical_quench`, `quinone_cys_michael`, `lysinoalanine_crosslink`, `aa_ring_open_dicarbonyl`, `pe_schiff_base`, `pe_amadori`, `asparagine_sugar_explicit_water_cluster`.
**React-OT eligible (CHON):** `lysinoalanine_crosslink`, `aa_ring_open_dicarbonyl`, `pe_schiff_base`, `asparagine_sugar_explicit_water_cluster`.
**Provenance matrix already published**: [results/validation/qm_barrier_provenance.md](results/validation/qm_barrier_provenance.md).

---

## Parked — 2026-04-19 pyGSM integration & global barrier run

> **STATUS: PARKED pending literature ingestion.**

---

## Maintenance Backlog (cheap, behavior-preserving)

### DFT Ladder Bulletproofing (2026-04-17)
- [x] Preserve best-available geometry fallback for non-TS optimizations instead of raising the last DFT exception. Verified by a new regression in [tests/unit/test_dft_refiner_fallback.py](tests/unit/test_dft_refiner_fallback.py) and a full Docker unit pass (**633 passed, 1 skipped**).
- [ ] Make trajectory resume authoritative: skip MLP/xTB re-relaxation when a DFT-derived checkpoint is recovered.
- [ ] Harden the high-level single-point stage against runtime exceptions, not only SCF non-convergence.
- [x] Add focused unit coverage for non-TS best-geometry SP fallback, authoritative resume, and SP exception fallback. Three new tests in [tests/unit/test_dft_refiner_fallback.py](../tests/unit/test_dft_refiner_fallback.py): `test_non_ts_sp_uses_best_available_fallback_geometry` (proves SP+Hessian consume the fallback xyz, not the original input), `test_geometry_ready_resume_skips_mlp_xtb_and_dft_ladder` (authoritative resume), `test_single_point_recovers_from_solvent_scf_exception` (SP runtime-exception fallback).
- [ ] Revalidate in Docker on the computational-gap lane before any heavy DFT restart.

### Architecture Hardening (2026-04-13)
- [ ] Move LaTeX-backed figure config into one shared helper; fail explicitly if the toolchain is absent. **Note**: also closes the LaTeX-fallback contract referenced by the (now-closed) UX sprint Lane B.
- [ ] Reduce duplicated plot-style bootstrapping across validation generators.
- [ ] Split computational-gap execution state from scientific promotion state more cleanly (a synthesized xTB path with SCF warnings must not look equivalent to a clean path downstream).
- [ ] Promote proxy-readiness diagnostics into machine-readable artifacts.
- [ ] Continue shrinking duplicated script-local helper logic when a shared `src/` helper makes the contract clearer.
- [ ] Keep public docs aligned with machine-readable artifacts.

---

## Strategic Backlog (defer until next sprint closes)

### S14 — Monolith decomposition
- [ ] `src/literature_runtime.py` (~4.3k LoC) → loader / query / calibration / governance.
- [ ] `src/benchmark_validation.py` (~2.6k LoC) → registry / evaluation / reporting / assertion.
- [ ] `src/recommend.py` (~1.3k LoC) → projection / observable mapping / scoring.
- [ ] Triage slowest 10 tests; introduce fast/slow/full marks.

### S19 — Minimal web interface
- [ ] FastAPI form → `run_pipeline` → embedded HTML report. Out of scope until adoption pull from a real scientist.

### Reproducibility & provenance
- [ ] Pin Docker image tag and embed it in every report artifact.

---

## Mission

Build the most useful computational tool for scientists who want to imitate meat-like Maillard chemistry in plant-based systems.

## Product Thesis

This is not mainly a problem of making one chemistry engine more exact. It is a problem of combining: (1) a quantitatively credible free-precursor core, (2) matrix-aware observability and accessibility, (3) process-aware confidence boundaries, (4) scientist-facing reporting that states what is benchmarked, what is transferred, and what is blocked.

Product question:
> *Which formulation and process changes are most worth testing next if the goal is meat-like aroma under plant-matrix constraints?*

## Success Criteria

- [ ] A scientist can use the tool to narrow a wet-lab or literature-backed matrix campaign, not just inspect simulated chemistry.
- [ ] At least one matrix lane becomes meaningfully closer to external decision readiness. **PARTIAL** — `## 7.` recommendation is now matrix-scoped, so soy/pea/wheat-gluten/myco formulations get a campaign focused on their own family rather than the global top miss (sprint 2026-05-13 S23).
- [ ] Every prediction carries explicit uncertainty bounds (90% CI). **DONE** (compound-confidence overlay + intervention waterfall, sprint 2026-05-10).
- [ ] A scientist can ingest new experimental data with a single command and see its impact before committing. **DONE** (`maillard ingest`, sprint 2026-05-10 Lane C).
- [x] The tool can rank candidate experiments by expected calibration gain and auto-generate protocols. **DONE** (`experiment-value-ranking`, `next-experiment`, and `## 7. Recommended next experiment`, sprint 2026-05-12).
- [ ] Transfer assumptions are validated via leave-one-out cross-validation. **DONE** ([src/cross_validation.py](src/cross_validation.py), sprint 2026-05-10 Lane D).
- [ ] Out-of-calibration predictions carry an explicit scope banner and demoted tiers. **DONE** (sprint 2026-05-10b Lane E).
- [ ] First 10 minutes of repo use produce a real report bundle without hand-assembling CLI flags. **DONE** (`./scripts/docker_maillard.sh quickstart`, sprint 2026-05-11).
- [x] All benchmarks pass inside the 1.5× acceptance band on the calibration set.
- [x] Trust headline 39/48 (81.2 %) holds under defensive scope.

---

## Recent closures

| Sprint | Date | Headline |
|---|---|---|
| 2026-05-13 DFT Ladder Bulletproofing — SP fallback / resume / SP-exception coverage | 2026-05-13 | Three new tests in [tests/unit/test_dft_refiner_fallback.py](../tests/unit/test_dft_refiner_fallback.py): non-TS best-geometry SP fallback (verifies SP+Hessian consume fallback xyz), authoritative geometry-ready resume (MLP+xTB+DFT-ladder all skipped), and `single_point` runtime-exception fallback (vacuum + hardened SCF retry). Closes item 4 of the DFT Ladder Bulletproofing backlog. Unit suite **657 passed / 1 skipped**. |
| 2026-05-13 CRO send-to-lab checklist + `analytical_context` | 2026-05-13 | Each `next-experiment` per-target protocol Markdown now ships an `## Analytical context` block (mirrors intake-schema field names) and an 8-step `## CRO send-to-lab checklist`; the pre-filled intake YAML carries default `analytical_context` keyed off the matched DoE template + compound-specific isotopically-labeled internal-standard hints (with a slug-keyed placeholder fallback). Returned YAML lands cleanly via `scripts/ingest_results.py` with no field translation. 3 new unit tests; unit suite **654 passed / 1 skipped**. |
| 2026-05-13 S24 Tier A.2 Sensory readout panel | 2026-05-13 | New `src/sensory_readout.py` (`compute_oav`, `classify_axis`, `roll_up_axes`, `build_sensory_readout`); `## 8. Sensory readout` block in `report.md` and per-formulation in `comparison.md` (axis roll-up + per-compound 90 % CI OAV table); 10 new unit tests; A.3 (tweak recommender) and A.1 (lipid-oxidation anchor curation) explicitly gated, not shipped. Unit suite **651 passed / 1 skipped**. |
| 2026-05-13 S23 Matrix-aware experiment recommendations | 2026-05-13 | `infer_matrix_family()` attributes every benchmark to a canonical family; `--matrix` flag on `experiment-value-ranking` + `next-experiment`; `## 7. Recommended next experiment` filters by the formulation's `protein_type` (with explicit fallback when no candidate matches). Unit suite 641/1. |
| 2026-05-12 S22 Most-Valuable-Experiment surface | 2026-05-12 | `experiment-value-ranking` + `next-experiment` Docker aliases, `## 7. Recommended next experiment` in report.md/comparison.md, scientist-facing VoI surface shipped. |
| 2026-05-11 Scientist-grade UX hardening | 2026-05-11 | Audience-grouped help, `quickstart` alias, docker-only docs, `## 6. Glossary` footer in report.md/comparison.md. Unit suite 632/1. |
| 2026-05-10b Defensive scope & honest priors | 2026-05-10 | `--scope-check` gate, lipid-oxidation `external_failing` posture, wheat_gluten/myco surrogates, ingest CSV templates, per-precursor intervention view. |
| 2026-05-10 Trust loop → scientist's hands | 2026-05-10 | External validation report (`0/8` honest baseline), confidence overlay + intervention waterfall in `--report`, `maillard ingest` CLI, observable-projection refit (rejected with explicit decision). |
| 2026-04-21 Family 14 TS Recovery | 2026-04-22 | Both recovery strategies failed; bounded_calibration HCW surrogate (7.58 kcal/mol) kept; `refinability_status=stack_limit_reached`. |
| 2026-04-21 SOTA TS pipeline (method B core) | 2026-04-2x | Multi-step pipeline implemented ([src/elementary_step_runner.py](src/elementary_step_runner.py), [src/microsolvation.py](src/microsolvation.py)). First multi-step target (`pe_amadori_water_shuttle`) crashes at GSM call — see proposed sprint above. |
| 2026-04-09 Maintainability Audit II / Surrogate Closure / Balanced Geometry / Git Tracking | 2026-04-09–10 | Helper extraction, surrogate posture cleanup, balanced geometry recovery, generated-artifact untracking. |
| 2026-04-13 QM Queue Realignment / Scientific Delivery Sequence | 2026-04-13 | Folded into `2026-04-22 Triangulated Barrier Strategy`. |

Detail for any closure is in [tasks/.todo.md.bak](.todo.md.bak), git history of `tasks/todo.md`, and `/memories/repo/`.
