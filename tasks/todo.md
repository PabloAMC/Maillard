# Maillard Strategic Roadmap

> **Backup of the previous long-form roadmap**: [tasks/.todo.md.bak](.todo.md.bak). Closed sprints have full review blocks preserved there and in git history; this file is the live working surface only.

## Priority order (set 2026-05-13 after the Tier A audit)

> **Operating constraint (user, 2026-05-13):** *"Do not complicate the tool more than it can do really really well."* Every new lane below has a stop-the-line gate that prevents adding a feature that the kernel cannot back with calibrated evidence.

> **PROPOSED (2026-09-01) — `cleaning` Phase 2: data restructure.** Planning only; see [tasks/data_restructure_plan.md](data_restructure_plan.md). Prerequisite: the `cleaning` branch does not pass CI (~25 test files + `scripts/run_pipeline.py` import modules deleted in Phase 1, `2dbe6a9`). Eight owner decisions listed in §5 of the plan.

1. **CLOSED (2026-06-22) — S27 Matrix lipid-oxidation recalibration (Workstream A) — landed as a refactor + a KEY NEGATIVE RESULT.** Goal was to close the matrix lipid-aldehyde over-prediction driving the external hold-out (0/8, median 36x). Root cause traced: hydroperoxide kinetics in `src/lipid_oxidation.py` are an unbounded Arrhenius × time model; the per-process_state matrix calibration registry pins in-panel to measured (in-panel reads 1.0x), but external-only hold-out hits unpinned process_states and the raw load leaks through (roasted-pea 160C 2748x, SPI/WG HME 24–748x).
   - **Finding**: A1 saturation CANNOT fix the hold-out without regressing the headline. In-panel (progress ~0.06 @100C/45min) and hold-out (HME ~0.07) share the same kinetic regime; only roasted-160C is true extrapolation. A cap strong enough to bend the hold-out perturbs in-panel trace lipid points (Hexanal/Nonanal ~2e-4 ppb, near-zero-width MC CIs) → 37/48 dropped to 31/48 even at c=12. ⇒ the gap is a per-(matrix, process_state) **calibration** problem, not a kinetic-shape problem.
   - **Shipped**: A1 saturation mechanism implemented but **DISABLED by default** (`max_conversion_fraction: null`) — default path is byte-identical to pre-S27, headline stays **37/48**, hold-out unchanged at 36.02x. A2 (externalize scale/Ea/branching/profiles to `data/lit/lipid_oxidation_calibration.json`, per-matrix `LipidProfile`, shared kinetics core) + A3 (shared core, per-stage branching retained) landed as a behavior-preserving refactor. Harnesses `scripts/diagnose_lipid_bias.py` + `scripts/calibrate_lipid_oxidation.py`; tests `tests/scientific/test_lipid_oxidation_saturation.py` (48 passing in the lipid/matrix/projection set). Gate honored: no invented numbers; 4 `external_validation/` bundles stayed frozen.
   - **Next (Workstream B, NOT YET DONE)**: process-state-aware uncertainty — detect uncalibrated process_states and widen CI / lower trust there instead of emitting confident wrong point predictions; size matrix observable sigma from residuals (DONE 2026-08-26: leave-lane-out derivation gives RMS ln-sigma 2.86, 90% CI [1.98, 5.48]; shipped 2.0 sits at the lower edge — see results/validation/matrix_sigma_residual_derivation.md); tier barrier sigma by provenance. This is the real lever for hold-out coverage. Plus curation (the original #4 lane) for high-T/extrusion anchors.
2. **NO OTHER ACTIVE SPRINT** — last closures: `2026-05-13 CRO send-to-lab checklist + analytical_context` (unit suite 654/1) and `2026-05-13 DFT Ladder Bulletproofing` (unit suite 657/1).
3. **DEFERRED (gate before doing) — S25 Tier A.3 — Formulation-tweak recommender.** Promising but risks overpromising: a one-shot grid sweep over `±cysteine`, `±ribose`, `±lipid`, `±pH`, `±temp` would feel like a recipe optimizer when the underlying kinetics are still calibrated against ≤48 benchmark targets. **Gate**: ship only after a scientist asks for it AND we can pin a reproducibility test that shows the top-3 tweaks come back with non-overlapping 90 % CIs on the meaty axis. Until then the existing `--report` already tells them which precursor is rate-limiting (`bottleneck_precursor`).
4. **SUPERSEDED by S27 (2026-06-22) — S26 Tier A.1 — Lipid-oxidation external anchor.** Was PARKED as a curation-only task; the model-side recalibration is now in progress under S27 above. Curation gate still stands for any NEW hold-out anchor: do not invent numbers; the existing 4 `external_validation/` bundles remain the frozen test set.
5. **PARKED — `pe_amadori` multi-step DFT closure**. Real chemistry win but heavy compute; pick up once a scientist actually asks for that barrier (or once `next-experiment` ranks it into the top-3).
6. **MAINTENANCE BACKLOG** — `2026-04-17 DFT Ladder Bulletproofing` (3/4 items remaining) + `2026-04-13 Architecture Hardening`. Cheap, behavior-preserving; pick up between sprint lanes.
7. **GATED — observable-first** — `2026-04-22 Triangulated Barrier Strategy` and `2026-04-19 pyGSM`. PAUSED until ≥1 wet-lab observable for Families 11–15 lands via `maillard ingest`. Note: the live top recommendation is still `cys_glucose_150C_Farmer1999 :: 2-methyl-3-furanthiol` (VoI 7.70).
8. **STRATEGIC BACKLOG** — `S14 monolith decomposition`, `S19 web interface`. Defer.
9. **DONE / SUPERSEDED** — see [Recent closures](#recent-closures).

---

## Closed Sprint — 2026-05-13 S24 Tier A.2 — Sensory readout panel

**Branch (proposed)**: `s24-sensory-readout`. **TL;DR**: Translate the kernel's per-compound `predicted_p5/p50/p95` ppb into Odour Activity Values (OAV = predicted / odour_threshold) and roll them into three axes (`meaty`, `off_note`, `safety`) using the existing keyword vocabulary in [src/experiment_value.py](src/experiment_value.py). Add a `## 8. Sensory readout` block to `report.md` and per-formulation in `comparison.md`. Additive only: zero changes to the kinetic kernel, the recommender, or the projection layer.

### Why this slice (and why nothing else from Tier A this sprint)
A scientist looks at the report and sees `2-methyl-3-furanthiol predicted_p50 = 0.42 ppb` and has to mentally divide by the odour threshold to know whether it's perceptible. We already curate ODT for every relevant compound (`data/species/desirable_targets.yml`, `off_flavour_targets.yml`) and `load_compound_specs()` in [src/experiment_value.py](src/experiment_value.py) already exposes them. This is a 1-day translation layer with high scientist value and zero new science. A.3 (tweak recommender) and A.1 (lipid-oxidation anchor) were audited but are **explicitly out of scope** — see Priority order #2 and #3 above for the gates.

### Lane S24-A — `src/sensory_readout.py`
- [x] A.1 `compute_oav(predicted_ppb, odour_threshold_ug_per_kg)` returns `None` when ODT is missing/non-positive/NaN; clamps negative predictions to 0.
- [x] A.2 `roll_up_axes(rows)` keyed by `meaty | off_note | safety` using the same keyword classifiers as `_suggest_template`. Excludes compounds without ODT from `max_oav` / `above_threshold_count` but counts them in `compounds_without_odt`.
- [x] A.3 `build_sensory_readout(result)` returns per-compound OAV (with p5 / p50 / p95) + axis roll-ups + a one-line headline.

### Lane S24-B — Report wiring
- [x] B.1 `_render_sensory_readout_markdown(result, *, heading="## 8. Sensory readout")` in [src/reporting.py](src/reporting.py) emits an axis roll-up table and a per-compound OAV table sorted by descending OAV (compounds without ODT pushed to the bottom).
- [x] B.2 Called from `generate_report` immediately after `_render_next_experiment_markdown`.
- [x] B.3 Called inside the per-formulation loop in `generate_comparison_report` with `heading=f"### {res.name} — sensory readout"` so the heading hierarchy stays sane.
- [x] B.4 Graceful degradation: empty `predicted_ppb` renders a one-line stub ("nothing to score"); a missing-ODT inventory note is emitted when at least one compound has no curated threshold.

### Lane S24-C — Tests
- [x] C.1 New [tests/unit/test_sensory_readout.py](tests/unit/test_sensory_readout.py) covers `compute_oav` math + `None` paths, axis classifier (incl. safety > meaty > off-note precedence), axis roll-up with mixed ODT presence, `build_sensory_readout` against a fixture `FormulationResult` with curated specs (MFT/Hexanal/Acrylamide), and the empty-prediction path.
- [x] C.2 [tests/unit/test_usability_reports.py](tests/unit/test_usability_reports.py) asserts `## 8. Sensory readout` appears after `## 7.` in `report.md` and that `comparison.md` carries a per-formulation `— sensory readout` heading.
- [x] C.3 `pytest tests/unit -q` in Docker → **651 passed, 1 skipped, 2 warnings** in 285.89s (was 641/1, +10 new).

### Out of scope this sprint (held to)
- A.3 (formulation-tweak recommender). A.1 (lipid-oxidation anchor curation).
- New keywords / new axis classifiers (e.g. green, fatty, roasted). Punt to a real scientist request.
- Per-axis weighting / sensory PCA. Same gate.
- Any change to the kinetic kernel, projection layer, or VoI ranker.

### Review
- **Shipped**: [src/sensory_readout.py](src/sensory_readout.py) (`compute_oav`, `classify_axis`, `roll_up_axes`, `build_sensory_readout`); `## 8. Sensory readout` in `generate_report` and a per-formulation block in `generate_comparison_report`; 10 new unit tests; full unit suite green at 651/1.
- **Scientist value**: every report now states OAV with explicit p5/p95 bounds for every predicted compound that has a curated odour threshold, plus a three-line axis roll-up so a scientist sees "meaty above threshold; off-notes below threshold; safety clear" without doing arithmetic.
- **Held the line**: Tier A.3 (tweak recommender) and A.1 (lipid-oxidation anchor) were not shipped despite being on the audit list — both have explicit gates in Priority order #2 and #3 of this file.
- **Follow-up trigger**: open S25 (A.3) only when a scientist asks for it AND a reproducibility test for non-overlapping CIs across top-3 tweaks is in place. Open S26 (A.1) only when a real SIDA-grade lipid-oxidation paper is in hand.

---

## Closed Sprint — 2026-05-12 S22 Most-Valuable-Experiment surface

> Closed; full review collapsed into the [Recent closures](#recent-closures) table. Long-form review is preserved in git history and `tasks/.todo.md.bak`.

---

## Parked — 2026-05-12 `pe_amadori` multi-step DFT closure

> **STATUS: PARKED behind S22.** The multi-step pipeline ([src/elementary_step_runner.py](src/elementary_step_runner.py) + [src/microsolvation.py](src/microsolvation.py)) is in place; the only configured target (`pe_amadori_water_shuttle`) crashes at the GSM call with `pyGSM exception: IndexError: index 34 is out of bounds for axis 0 with size 33` (see [results/computational_gap_refinement/multistep/pe_amadori_water_shuttle/target_result.json](results/computational_gap_refinement/multistep/pe_amadori_water_shuttle/target_result.json)). When unparked: M1 fixes the IndexError; M2 encodes the real two-step mechanism (O0→C2 then C1→N3); M3 writes back uncertainty narrowing only (no anchor promotion). Hard timebox on M2 (~24h DFT wall per step).

---

## Active Plan — 2026-04-22 Triangulated Barrier Strategy for Families 11–15

> **STATUS: GATED — observable-first.** Reactivation criterion: ≥1 wet-lab observable for any of Families 11–15 lands via `maillard ingest`. Order then is **D → C(inventory only) → B → A**. Hard stop: if Vía B (React-OT) yields seeds for ≥2/4 eligible targets and the rest are classified by provenance, Vía A (OA-ReactDiff) is skipped.

**Targets (7):** `hexanal_radical_quench`, `quinone_cys_michael`, `lysinoalanine_crosslink`, `aa_ring_open_dicarbonyl`, `pe_schiff_base`, `pe_amadori`, `asparagine_sugar_explicit_water_cluster`.
**React-OT eligible (CHON):** `lysinoalanine_crosslink`, `aa_ring_open_dicarbonyl`, `pe_schiff_base`, `asparagine_sugar_explicit_water_cluster`.
**Provenance matrix already published**: [results/validation/qm_barrier_provenance.md](results/validation/qm_barrier_provenance.md).

### Lanes (collapsed; details in git history of this file)
- **B (React-OT pilot)** — wrapper [scripts/recover_ts_react_ot_seed.py](scripts/recover_ts_react_ot_seed.py) is in repo; `models/external/react_ot/` has `provenance.json`. Pending: run on the 4 eligible targets when the gate opens.
- **A (OA-ReactDiff fallback)** — only if B fails on ≥3/4 eligible targets. Not yet installed.
- **C (wet-lab backlog)** — paralelo, no bloquea. Spec writing only; no experiments without trigger.
- **D (governance)** — `qm_barrier_provenance.{md,json}` exists; `aa_ring_open_dicarbonyl.refinability_status="stack_limit_reached"` already committed.

---

## Active Plan — 2026-04-19 pyGSM integration & global barrier run

> **STATUS: ⏸ PAUSED — gated by observables.** Same gate as above. The Sella + GFN2 + r2SCAN/def2-svp stack already failed twice on Family 14; without a benchmark observable a successful pyGSM run still moves nothing in `--report`. Do not start Phase 1.1 until the gate opens.

7 targets, 4 phases (backend → 3-target validation → 7-target diagnosis → global run). Detail preserved in git history; reactivate only after the gate opens.

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
