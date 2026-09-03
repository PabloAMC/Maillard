# Data restructure plan — `cleaning` Phase 2

> Planning round, 2026-09-01. Nothing here has been executed. Companion to the Phase 1 prune
> (`2dbe6a9`). Evidence comes from a full read of `data/`, `src/`, `scripts/`, `tests/`, the CI
> gates and `AUDIT.md`; every claim below cites a path so it can be re-checked.

## 0. TL;DR

`data/` is not one database. It is five overlapping file stores (literature registries, benchmarks,
species, protocols, QM fixtures) plus 400 MB of untracked scratch, with **no shared keys, no
schemas, no single access path, and no read-only boundary**. The runtime reaches it through
95 modules that each re-derive `repo_root / "data" / ...`, joins records by free-text names
through five alias tables and a fuzzy matcher, silently substitutes `{}` when a file is missing,
and writes generated artifacts back into the "curated inputs" tree.

The plan is six phases. Phases 0–2 are behaviour-preserving housekeeping that should land
regardless (fix the broken branch, purge, one access layer). Phase 3 (keys) is the real fix and
the only phase that can change model outputs. Phases 4–5 (schemas, docs) are what keeps it
clean afterwards.

**Blocking prerequisite:** the `cleaning` branch does not currently pass CI. Phase 1 deleted
36 `src/` modules and 5 `scripts/`, but ~25 test files and `scripts/run_pipeline.py` still
import them (§4, Phase 0).

---

## 1. Diagnosis — what "frankenstein" means, concretely

### D1. No access layer: 95 modules hard-code paths, some relative to CWD
- 50 modules in `src/` and 45 in `scripts/` each contain `Path(__file__).resolve().parents[1] / "data" / ...`
  (54 occurrences) or a literal `"data/..."` string (58 occurrences). No `src/paths.py`, no
  env var for the data root. Six modules derive the same path twice in one file.
- `src/headspace.py:59` defaults to the CWD-relative string `"data/lit/henry_constants.yml"`;
  `HeadspaceModel()` is called with no argument at `src/benchmark_validation.py:844` and ~20 test
  sites. If CWD ≠ repo root, the constants load as `{}` and every Kaw becomes the hard-coded
  `0.01` (`headspace.py:69`) with no error.
- `scripts/calibrate_barriers.py:263` writes CWD-relative. `tests/integration/test_regression.py:16`
  reads CWD-relative.
- Path strings are baked into report output (`src/reporting.py:1106-1169`, 18 literals) and four
  tests assert those exact strings, so the physical layout is part of the public output contract.

### D2. `data/` is read-write: 11 write paths, generated artifacts committed as "curated"
| Writer | Writes into `data/` | Readers |
|---|---|---|
| `src/matrix_calibration_optimizer.py:415` | `calibration_history/calib_*.json` (102 files, all pytest output from `tests/integration/test_matrix_calibration_loop.py:16`, which does not monkeypatch the dir) | none |
| `src/experiment_request.py:419` | `protocols/requested_*.yaml` (5) | none; MD twins live in `results/validation/experiment_requests/` |
| `src/external_validation.py:1016-1029` | `benchmarks/external_validation/*.json` + `protocols/external_validation/*.yaml` (same 4 names, two formats) | hold-out lane |
| `scripts/deep_research_tracker.py:256`, `scripts/sync_backlog.py:110`, `scripts/ingest_deep_research_markdown.py:182` | `lit/deep_research_backlog.json` (430 KB, three independent writers) | 2 src, 5 scripts |
| `scripts/calibrate_barriers.py:263` | `lit/calibration_offsets.json` | none (README says unwired) |
| `scripts/generators/generate_slr_family_reports.py:559` | `Gemini_Deep_Research/slr_family_*.md` (16) | none |
| `scripts/generators/refresh_internal_reproducibility_snapshots.py` | `benchmarks/*_{Internal2026,ProtocolPilot2026}.json`, `protocols/*_protocol_pilot_intake.yaml` | tests |
| `src/matrix_recalibration.py:192` | `lit/matrix_calibration_offsets.json` (absent on disk; read back by `matrix_calibration_registry.py:12`) | runtime feedback loop through `data/` |
- Status ledgers with counts live in `data/lit`: `dft_coverage_map.json` (0 readers),
  `slr_incorporation_matrix.json` (350 KB; hand-maintains `current_runtime_consumers`, i.e. the
  import graph, by hand). `agents.md:20` forbids exactly this.
- Mirror image: `results/validation/qm_barrier_provenance.json` is a hand-curated input with no
  generator, read at runtime by `src/uncertainty_propagation.py:513`. `docs/validation/directional_claims_panel.yml`
  and `directional_accuracy_report.md` are curated inputs in `docs/`, parsed by `src/directional_reliability.py:57-62`.

### D3. No keys: every join is a string match
- **Compounds.** Keyed by display name everywhere. In `data/benchmarks/**` alone: `Hexanal`/`hexanal`,
  `Acrylamide`/`acrylamide`, `2-Methyl-3-furanthiol (MFT)`/`2-methyl-3-furanthiol`, and
  `Furan-2-aldehyde (FA)` = furfural; 25 keys for ~19 compounds. `flavor_reference_payloads.json`
  uses `compound`, `safety_reference_payloads.json` uses `analyte` (7 rows `null`), species YAML
  uses `name`. The runtime papers over this with **five** independent normaliser/alias tables
  (`benchmark_validation.py:90-140`, `external_validation.py:61-67`, `experiment_value.py:146-198`,
  `sensory.py:170-212`, `matrix_calibration_registry.py:773`) and a `difflib.SequenceMatcher`
  fallback at ratio 0.75 (`benchmark_validation.py:551-582`). `test_regression.py:53-55` keeps its
  own synonym table. SMILES/InChI exist in the YAMLs but are never used as keys.
- **Papers.** 171 of 227 distinct DOIs appear in more than one file; no shared id. DOI
  `10.3390/molecules28073151` has 13 ids across 4 files; `10.1590/1678-457x.08717` has 3 ids in
  `benchmark_intake_registry.json` alone with wrong years. PDF stems, dossier stems, benchmark
  ids, intake ids and SLR ids are five namespaces (14-row mismatch table in Appendix B), including
  two documented file-swap traps (`Kocada2016.pdf` vs `Kocadagli2016.pdf`, the two `Zhai2023` files).
- **Reaction families.** Four keyspaces: `arrhenius_params.yml` (`enolisation`), `FAST_BARRIERS`
  in `src/barrier_constants.py` (`enolisation_intermediate`), `calibration_offsets.json` (`enol`),
  `reaction_benchmark_set.json` (`1,2-enolisation`). Bridged by hard-coded dicts at
  `barrier_constants.py:740-748, 807-813`. The YAML header itself records a kJ/kcal collision.
- **Chemistry families.** 16 canonical slugs; `benchmark_intake_registry.json` uses 5 aliases
  across 11 rows; `slr_section` mixes `01`, `5`, `11` and `11.0`, `1.1`, `"1.3,6.1,7.4"`.
- **Matrix families.** 8 canonical values in `matrix_family_coverage_registry.json`; the intake
  registry uses 86. Two vocabularies (`pea_isolate` vs `pea_iso`) bridged by one alias map
  (`literature_runtime.py:103-110`) and a substring rule over benchmark ids (`experiment_value.py:52-62`).

### D4. No schemas, one validator, and it drifted from its own schema file
- Benchmark JSON has 4 top-level variants and 9 distinct `measured_volatiles` key sets; `ph` is
  bench-initial in 20 bundles and post-autoclave in one (disclosed only in prose). No validator:
  `src/benchmark_validation.py` validates model vs. measurement, not structure.
- `computational_priors.json` (217 KB) is 22 unrelated tables with 22 record shapes in one file.
- The only validator (`src/matrix_experiment_intake.py:64`) hard-codes the enums in Python
  (`:32-34`) instead of reading `matrix_experiment_intake_schema.json`, and they already disagree
  (`wheat_gluten`). The `requested_*.yaml` it supposedly emits fail it three ways.
- No pydantic/jsonschema anywhere. `tests/unit/test_data_integrity.py` reads no data file. No
  referential-integrity test between any two data files.
- Confidence vocabulary: `agents.md` mandates `bounded_calibration / transferred_literature /
  surrogate_family / xtb_derived`; three of the four occur **zero** times in `data/`. The real
  fields are `confidence_tier` (with synonyms `low_medium` and `medium_low` both live),
  `provenance_tier`, `uncertainty_posture`; 19 of 27 JSON files carry none.

### D5. Failure is silent
Most loaders return `{}` on a missing or malformed file (`literature_runtime.py:40`,
`barrier_constants.py:530`, `lipid_oxidation.py:76`, `recommend.py:74`); two substitute
hard-coded fallback constants that then look data-derived. 11 files are read at import time
into module globals across 7 modules. A stale path after a move degrades the model quietly
rather than failing.

### D6. Dead, orphaned and ghost files
- **Zero consumers:** `data/rmg_extensions/` (13 families; no `rmgpy` in the env),
  `data/reactions/curated_pathways.py` (imports the deleted `src.pathway_extractor`; stale copy of
  `src/curated_pathways.py`), `data/reactions/reaction_families.yml` (existence-checked only),
  `lit/dft_coverage_map.json`, `lit/ts_seed_benchmark_set.json`, `lit/computational_gap_closure_targets.json`,
  `lit/computational_gap_multistep_targets.json`, `lit/calibration_offsets.json`, all 33 tracked
  `data/geometries/` files (readers deleted in Phase 1; 15 are tracked while gitignored),
  `results_db.ml_adoption_decisions` table + 2 methods, `docs/notebooks/results/maillard_results.db`.
- **Test-only:** `lit/canonical_systems.json`.
- **Ghosts (referenced, absent):** `lit/ribose_glycine_2021.json` (README), `lit/matrix_calibration_offsets.json`
  (read at `matrix_calibration_registry.py:12`), `lit/deep_research_candidate_registry.json` (written
  by `ingest_deep_research_markdown.py`, cited by the intake registry's notes), `src/dft_coverage_map.py`
  (cited by `generate_family_implementation_status.py:199`), `data/benchmarks/quarantined/{cys_ribose_150C_Mottram1994,cys_glucose_150C_Farmer1999}.json`
  (`calibrate_barriers.py:70-71`), `Gemini_Deep_Research/raw/11_maillard_lipid_crosstalk.md`
  (`computational_priors.json:2118`), `xtb_inputs/pe_amadori/{reactant,product}.xyz`, `data/reactions/rmg_validation_cases/`.

### D7. Untracked bulk hides inside `data/` — the `data/qm` lesson, still live
`.gitignore` blanket-ignores `data/*` then whitelists. On disk but invisible to git and every
audit: 818 geometry files (31 MB, generated QM scratch), `data/logs/` (193 MB of dead stdout,
zero writers, zero readers), `data/calibration_history/` (102 files), `data/articles/` (134 PDFs,
162 MB, genuine sources, never cited by filename from any registry). AUDIT.md §"Open items 3"
documents what this pattern already cost once.

### D8. Documentation describes a different repository
`data/lit/README.md` documents one ghost file, omits five real ones (`binding_constants.yml`,
`chemistry_family_scope_registry.json`, `extrusion_damage_reference_payloads.json`,
`lipid_oxidation_calibration.json`, `process_gap_registry.json`), and calls a run ledger
"curated". `timeseries/*.yml` headers say "NOTHING IN THIS FILE IS WIRED" while the README marks
two as fit corpus. `src/external_validation.py:7` says 19 benchmarks; there are 23.

### D9. Regeneration hazards
- `scripts/generators/complete_benchmark_buffer_fields.py` patched `conditions.buffer` into all
  21 hold-out bundles, but `_HOLDOUT_BUNDLE_SPECS` (`src/external_validation.py:71-84`) has no
  `buffer` key. Re-running `generate_external_validation_payloads.py` deletes the provenance
  blocks. This already happened once (`external_validation.py:113-116`).
- `generate_slr_family_reports.py:44-56` deletes nine legacy filenames on every run.

### D10. The `Internal2026` / `ProtocolPilot2026` pairs are one synthetic snapshot under two labels
Byte-identical `conc_ppb` for all 7 compounds; one is filed as `measured_volatiles`, the other as
`reference_volatiles`. Calling model output "measured" is the confusion the rest of the repo
works hardest to avoid. Heavily test-referenced (17 + 15 + 3 + 3 files), so retire with care.

---

## 2. Design principles for the target state

1. **`data/` is read-only at runtime.** Anything a script or test writes goes to `results/`
   (or `tests/` tmp). Enforced by a CI gate that diffs `data/` after the test run.
2. **One key per entity.** Compounds by InChIKey, papers by `paper_id` (one per DOI), chemistry /
   matrix / reaction families and process states by closed enumerations in `data/keys/`. Every
   other file references keys; no file re-spells a name.
3. **One record envelope.** Adopt the shape `binding_constants.yml` and `timeseries/*.yml`
   already use: `sources[]` + `records[]`, each record with `id`, `paper_id`, `quote`,
   `verification_status`, one tier field. Flat `entries[]`, never nested class→protein→list.
4. **One access path.** `src/data_paths.py` (named path constants, one data root, env override)
   and `src/data_access.py` (typed loaders; raise on missing; `lru_cache`; no import-time I/O).
   `grep '"data/' src` must return only `data_paths.py`.
5. **Schemas in the repo, enforced in CI.** JSON Schema per record type under `data/schemas/`,
   a stdlib+pyyaml `scripts/ci/schema_gate.py` next to the three existing gates, and one
   referential-integrity test.
6. **Generated ≠ curated ≠ forensic.** Three trees, three rules: `data/` (hand-edited, cited),
   `results/` (regenerable, never hand-edited), `data/research_corpus/` (LLM dumps kept only so
   `citation_gate.py` can detect laundered citations; never a provenance source).
7. **Names say what things are.** No wave shorthand (`k4b_`, `phase33_`), per `agents.md`.

---

## 3. Target layout

```
data/
  README.md                         one line per file; CI checks no ghosts, no omissions
  keys/                             NEW — the join layer
    compounds.yml                   id=InChIKey; name, cas, smiles, aliases[], class,
                                    odour_threshold{value,unit,paper_id}, sensory_tags[]
                                    ← species/{desirable,off_flavour,toxic}_targets.yml,
                                      sensory_tags.yml, lit/matrix_decision_panel.json,
                                      all in-code alias tables
    precursors.yml                  ← species/precursors.yml, one section schema, inchikey
    papers.yml                      one row per DOI: paper_id, citation, doi, pdf, dossier
    chemistry_families.json         ← chemistry_family_scope_registry + family_ingestion_plan
                                      (+ the 10 rows of dft_coverage_map)
    matrix_families.json            ← matrix_family_coverage_registry; long+short ids in one row
    reaction_families.json          NEW: arrhenius key ↔ FAST_BARRIERS key ↔ benchmark key
    process_states.json
  schemas/                          NEW — JSON Schema per record type
    benchmark.schema.json  reference_record.schema.json  prior_record.schema.json
    intake.schema.json (← protocols/matrix_experiment_intake_schema.json)
  constants/                        measured physical constants
    arrhenius_params.yml  henry_constants.yml  binding_constants.yml
    lipid_oxidation_calibration.json  interventions.yml
  priors/                           ← computational_priors.json split by section (22 → ~10 files)
    matrix_accessibility.json  thiamine.json  nucleotide.json  carbonyl_donor.json
    sulfur_peptide.json  lipid_offnote.json  ascorbic.json  melanoidin.json  dft_kinetic.json …
  references/                       observable anchors, one flat entries[] envelope each
    flavor.json  safety.json (+extrusion_damage)  retention.json  process_state.json
  benchmarks/
    calibration/*.json              today's top level (minus the synthetic pair, see decision 6)
    holdout/matrix/*.json           ← external_validation/
    holdout/free_precursor/*.json   ← external_validation/maillard_path/
    retired/quarantined/            retired/unreachable/
  timeseries/                       unchanged (already the model)
  dossiers/<paper_id>.md            ← lit/extraction_dossiers/; synthesis/ for k1…k6b, round3/4
  intake/
    benchmark_intake_registry.json  slimmed: references paper_id; derived columns dropped
    gaps.json                       ← process_gap_registry + intake.structural_gaps
  protocols/
    contracts/*.json                the 4 *_contract.json, one shape
    example_matrix_experiment_intake.yaml
  qm/
    reaction_barriers.json          ← lit/reaction_benchmark_set.json
    geometries.json                 ← lit/geometry_benchmark_set + ts_seed_benchmark_set
    mlp_registry.json               ← mlp_candidate_registry + mlp_external_benchmark_evidence
    qm_barrier_provenance.json      ← results/validation/ (curated input, misfiled)
    unverified/phase33_*.json, phase35_*.json   (or deleted — decision 2)
  formulation_grid.yml
  articles/                         untracked PDFs + tracked manifest.json {file→doi,paper_id,dossier}
  research_corpus/                  ← Gemini_Deep_Research/raw + 13 dumps; forensic only

docs/examples/                      ← cli_examples/, campaigns/
docs/templates/                     ← ingest_templates/
docs/reference/reaction_families.yml  ← data/reactions/ (a design doc)
docs/retired/maillard_validation_benchmarks.md  (self-declared fabricated sections)

results/                            generated, gitignored unless a test/gate needs it
  validation/experiment_requests/*.{yaml,md}   ← protocols/requested_*.yaml
  validation/reproducibility/                  ← protocols/*_protocol_pilot_intake.yaml
  validation/holdout_bundles/*.yaml            ← protocols/external_validation/
  literature/deep_research_backlog.json  slr_incorporation_status.json  dft_coverage_map.json
  literature/slr_family_reports/*.md
  calibration/calibration_offsets.json  calibration_history/

tests/fixtures/canonical_systems.json

src/data_paths.py   src/data_access.py   scripts/ci/schema_gate.py   scripts/ci/data_readonly_gate.py
tests/unit/test_data_referential_integrity.py

DELETE
  data/rmg_extensions/  data/reactions/curated_pathways.py  data/logs/
  data/geometries/ (33 tracked + 818 untracked)  data/calibration_history/
  docs/notebooks/results/maillard_results.db  results_db.ml_adoption_decisions (+2 methods)
  lit/{dft_coverage_map,ts_seed_benchmark_set,computational_gap_*}.json (after merge/move)
```

Whether the top-level stays `data/lit/…` or flattens as above is decision 7; the plan works either way.

---

## 4. Phased execution

Each phase is one PR, behaviour-preserving unless stated, with the unit + scientific suites run
in Docker before and after and the headline numbers (`tests/scientific/test_honest_headline_guards.py`)
pinned unchanged. Moves use `git mv` so history survives. Effort is in focused sessions.

### Phase 0 — Stabilise the `cleaning` branch (prerequisite, ~1 session) — IN PROGRESS 2026-09-01
Owner decisions received 2026-09-01: (i) full tree backup before deletions —
`../Maillard_backup_2026-09-01_pre-data-restructure.tar.gz` (298 MB, includes `.git` and all
untracked data) + git tag `pre-data-restructure-2026-09-01`; (ii) the QM/DFT lane is gone for
good, so its data, loaders, tests and docs are to be removed, not preserved.
- [x] 16 test files: `from src.pathway_extractor import Species/ElementaryStep` → `src.chem_utils`.
- [x] Deleted QM-only tests (owner-authorised): `tests/scripts/{test_run_computational_gap_dft,
      test_run_computational_gap_xtb,test_react_ot_seed_coverage,test_import_react_ot_colab_artifacts,
      test_open_react_ot_colab}.py`, `tests/scientific/{test_mlp_assessment,test_refinement_campaign,
      test_refinement_governance,test_refinement_watchlist,test_offline_refinement_governance}.py`,
      the `PathwayExtractor` test in `tests/unit/test_error_handling.py`, empty `tests/qm/`.
- [x] `tests/unit/test_wave_r1_barrier_offset_retirement.py`: dropped the `src.refinement_campaign`
      import; the retirement-note test now guards the tracked file directly.
- [x] `scripts/run_pipeline.py`: removed `XTBScreener` import and the `--xtb` flag.
- [x] `.github/workflows/ci.yml`: removed the `tests-qm` job and the `qm` dispatch input.
- [x] `scripts/docker_maillard.sh`: removed dead `run_computational_gap_job`.
- [x] D9 fixed properly: the four matrix hold-out bundles are frozen evidence. `write_holdout_bundles`
      is write-once by default (`overwrite=True` carries curated-only keys forward); the diff that
      forced this is bigger than `conditions.buffer` — Wave W primary-source corrections to the Bi 2020
      values/uncertainties live only in the committed JSON. Two new tests pin it.
- [x] `data/calibration_history` leak: writer moved to `results/calibration_history/`; the
      integration test monkeypatches the dir.
- [x] `data/reactions/curated_pathways.py` deleted (broken import, zero consumers).
- [x] Baseline in Docker (commit `8c35ec8`): unit+scripts **1347 passed**; integration+scientific **428 passed, 2 xfailed**.

### Phase 1a — QM/DFT lane purge (owner-authorised 2026-09-01) — IN PROGRESS
Everything the xTB → DFT → MLP refinement lane left behind, deleted rather than re-homed:
- [x] `data/qm/`, `data/geometries/` (33 tracked + 818 untracked), `data/rmg_extensions/`, `data/logs/`,
      `data/calibration_history/`, `models/external/`, `docs/guides/COMPUTATIONAL_GAP_RUNBOOK.md`,
      `docs/notebooks/results/maillard_results.db`, `environment.react_ot.yml`.
- [x] `data/lit/{dft_coverage_map, computational_gap_closure_targets, computational_gap_multistep_targets,
      ts_seed_benchmark_set, geometry_benchmark_set, mlp_candidate_registry, mlp_external_benchmark_evidence,
      reaction_benchmark_set, calibration_offsets}.json`.
- [x] `results/validation/`: `pygsm_scratch/` (106 files) + 31 QM-lane artifacts (refinement_*, computational_gap_*,
      selective_dft_plan, qm_barrier_provenance, react_ot_*, xtb_*, offline_dft_jobs, cheap_refinement_screening).
- [x] `src/{authority_benchmark_data, reaction_benchmark, solvation, skip_policy_registry}.py` + their tests;
      `scripts/{calibrate_barriers.py, generators/generate_rmg_inputs.py, generators/generate_skip_registry.py}`.
- [x] `src/results_db.py`: `ml_adoption_decisions` table + 2 methods dropped (table also dropped from the tracked DB).
- [x] `src/uncertainty_propagation.py`: `_apply_qm_provenance_narrowing` (a documented no-op) removed.
- [x] `src/reporting.py`: 21 QM-artifact keys dropped from the scientific surface.
- [x] `scripts/ci/citation_gate.py`: check 5 (barrier-source disclosure over `data/qm/`) retired; `data/qm` glob removed.
- [x] Census re-pinned 120/98/80 → 102/80/80 in `test_honest_headline_guards.py` and README; the data/qm split test removed.
- [x] Docs: README (QM section → "No quantum chemistry in the loop", QM-operator row), `agents.md`, `CONTRIBUTING.md`,
      `VALIDATION_CONTRACT.md`, CI header comments, `.gitignore` (all xtb/geometry rules).
- [x] Environment: `environment.yml` reduced to imported deps only; `Dockerfile`, `scripts/ci/setup_env.sh`,
      `docker_maillard.sh bootstrap/status`, `scripts/check_env.py` stripped of torch/jax/xtbiff/pyGSM/MACE. **The running
      container still has the old env; a rebuild is not required for tests to pass.**
- [x] Gates: citation_gate / holdout_guard / fit_target_gate PASS; collection 1768 tests, 0 errors.
- [x] Full suites in Docker (commit `22b5dee`): unit+scripts **1340 passed**; integration+scientific **426 passed, 2 xfailed**.

### Phase 1b — Re-home what is not curated input — IN PROGRESS
- [x] `data/cli_examples/*`, `data/campaigns/*` → `docs/examples/`; `data/ingest_templates/*` → `docs/templates/`;
      `data/lit/canonical_systems.json` → `tests/fixtures/`; `data/reactions/reaction_families.yml` → `docs/reference/`
      (so `data/reactions/` is gone).
- [x] `data/protocols/requested_*.yaml` → `results/validation/experiment_requests/` (next to their MD twins);
      `src/experiment_request.py` writes there now.
- [x] `data/protocols/external_validation/*.yaml` → `data/benchmarks/external_validation/intake/` (frozen evidence next
      to the JSON it materialized; `EXTERNAL_VALIDATION_PROTOCOL_DIR` repointed).
- [x] `data/Gemini_Deep_Research/slr_family_*.md` (16, generated) → `results/literature/slr_family_reports/`, untracked;
      generator repointed and creates the dir.
- [x] `.gitignore` data rules rewritten as an explicit whitelist; `git ls-files -i -c --exclude-standard` is empty.
- [x] `data/lit/README.md`: QM category removed, ghosts removed, five undocumented files added.
- [ ] Deferred to Phase 2 (many consumers, do after `data_paths.py`): `deep_research_backlog.json`,
      `slr_incorporation_matrix.json` → `results/literature/`.
- [ ] Not moved, on purpose: `data/benchmarks/maillard_validation_benchmarks.md` (cited as a fixed historical path in
      ~35 places; the audit narrative depends on the path); `data/protocols/*_protocol_pilot_intake.yaml` (the
      `comparison_contract` half is hand-written; only `conc_ppb` is refreshed) — decision 6 in §5 still open;
      `data/Gemini_Deep_Research/` rename → Phase 5.
- [x] Full suites in Docker (commit `623b850`): unit+scripts **1340 passed**; integration+scientific **426 passed, 2 xfailed**.

### Phase 1 — Purge and re-home (original checklist; SUPERSEDED by 1a/1b above, kept for the exit criteria)
Pure deletes and moves of files with 0–2 consumers. Expected model-output change: none.
- [ ] Delete dead: `data/rmg_extensions/`, `data/reactions/curated_pathways.py`, `data/logs/`,
      `data/geometries/**` (tracked and untracked), `data/calibration_history/`,
      `docs/notebooks/results/maillard_results.db`, `lit/ts_seed_benchmark_set.json` (fold its 2 rows
      into `geometry_benchmark_set.json`), `lit/computational_gap_*.json`, `lit/dft_coverage_map.json`,
      `lit/calibration_offsets.json` (+ `scripts/calibrate_barriers.py`, whose two fit targets no longer exist).
      `ResultsDB`: drop `ml_adoption_decisions` and its two methods.
- [ ] Move generated → `results/`: `protocols/requested_*.yaml`, `protocols/*_protocol_pilot_intake.yaml`,
      `protocols/external_validation/*.yaml`, `lit/deep_research_backlog.json`, `lit/slr_incorporation_matrix.json`,
      `Gemini_Deep_Research/slr_family_*.md`. Repoint the writers (§1 D2 table) and the 2+5 readers of the backlog.
- [ ] Move examples/docs → `docs/`: `cli_examples/`, `campaigns/`, `ingest_templates/`,
      `reactions/reaction_families.yml`, `benchmarks/maillard_validation_benchmarks.md`.
      `lit/canonical_systems.json` → `tests/fixtures/`. `results/validation/qm_barrier_provenance.json` → `data/qm/`.
- [ ] Rename `Gemini_Deep_Research/` → `research_corpus/` (raw + 13 dumps only); update
      `citation_gate.py` regexes (`:216` and the `raw/\d\d_` alternative).
- [ ] `data/articles/manifest.json` generated from the dossier §0 IDENTITY tables (`scripts/generators/build_article_manifest.py`);
      the file-swap traps in `src/kinetic_core/parameters_furanic.py:124-126` become data.
- [ ] Rewrite `.gitignore` `# Data rules` as an explicit whitelist with no tracked-but-ignored files
      (`git ls-files -i -c --exclude-standard` must be empty). Remove the geometry/xtb rules.
- [ ] Add `scripts/ci/data_readonly_gate.py`: fail if `git status --porcelain data/` is non-empty after the test tiers.
- [ ] Update the three CI gates' globs and the count pins in `test_honest_headline_guards.py:959-1016`.
- **Exit:** `data/` holds only hand-curated, cited inputs; on-disk size drops from ~400 MB to ~170 MB (articles) / ~8 MB tracked.

### Phase 2 — One access layer (~2–3 sessions) — IN PROGRESS 2026-09-01
Do this **before** any renames inside `data/lit`, so later moves are one-line edits.
- [x] `src/data_paths.py` (65+ constants, `MAILLARD_DATA_ROOT`/`MAILLARD_RESULTS_ROOT`, `rel()`, `benchmark_path()`).
- [x] `src/__init__.py` re-exports are lazy (PEP 562): `import src.x` no longer imports the pipeline and reads 11 files.
- [x] `src/data_access.py`: `load_json` / `load_yaml` / `load_mapping` raise `DataFileError` on missing/malformed,
      mtime-keyed cache, `missing_ok=True` returns `None` (never `{}`), `clear_cache()`. Not adopted yet (2b).
- [x] 2a (commit `97069f1`): ~110 hard-coded path sites in `src/` and `scripts/` migrated onto `data_paths`; 41 private
      `ROOT` definitions removed; CWD-relative defaults in `headspace.py` / `sensory.py` removed (probe: `HeadspaceModel()`
      from `/tmp` now loads 31 constants, previously 0); `matrix_calibration_offsets.json` feedback loop routed to
      `results/calibration/`; four tests with CWD-relative paths fixed (`test_temporal_fast` had pointed at a nonexistent
      ramp CSV and silently tested the non-ramp path); latent `NameError` in `report_html.residual_section` fixed.
      The three stdlib-only CI gates keep literal globs by design.
- [x] 2b: `data_access` adopted by every curated-file loader — `literature_runtime`, `matrix_correction`, `barrier_constants`
      (`_load_dft_anchor_metadata` no longer substitutes the inline fallback table; `get_arrhenius_params` no longer returns
      `None` for every family when the YAML is absent), `recommend` (`_load_henry_lookup`, `_load_yaml_db`; `_load_ramp` now
      raises on a missing/malformed user ramp instead of logging and running isothermally), `lipid_oxidation`,
      `protein_binding`, `headspace`, `literature_family_registry`, the six `_load_json` copies, `curated_pathways`,
      `safety`, `matrix_prior_registry`, `matrix_targets`, `benchmark_validation.load_benchmark` (now cached),
      `precursor_resolver`, `experiment_value`, `sensory`, `pipeline`. `DataFileMissing` subclasses `FileNotFoundError`
      and `DataFileMalformed` subclasses `ValueError`, so existing handlers keep working.
- [x] 2c decided AGAINST for now: (i) the 17 remaining import-time loads cost one cached parse each and nothing outside
      their modules imports the globals — making them lazy means rewriting internal references in `literature_runtime`
      (4.5k LoC), which belongs with the S14 monolith decomposition; (ii) moving `deep_research_backlog.json` /
      `slr_incorporation_matrix.json` to `results/` would silently drop 97 `no_verifiable_source` records and 331 DOI
      fields out of the census and the citation gate — do it in Phase 4 with an explicit re-pin, not as a file move.
- **Exit:** `grep -rn '"data/' src scripts | grep -v data_paths.py` is empty; `pytest` from a
      non-root CWD passes; model outputs byte-identical.

### Phase 3 — Keys (the real fix, ~3–5 sessions; the only phase that can change outputs) — IN PROGRESS 2026-09-01
Executed ADDITIVELY: the key files are generated from the existing data and every record keeps its own id;
integrity tests make every spelling / DOI resolve. Rewriting record ids is a later, separate step.
- [x] 3a `data/keys/compounds.yml` (74 entries: 61 molecules keyed by RDKit InChIKey, 2 classes, 6 marker sets,
      5 process levers) generated by `scripts/generators/build_compound_registry.py` from the species YAML + the
      decision panel + 27 hand-entered structures (labelled as such). `src/compound_keys.py` is the one resolver
      (exact, normalised; no similarity). `tests/unit/test_compound_keys.py`: every compound spelling under `data/`
      resolves, no spelling claimed twice, InChIKeys unique, the retired `BENCHMARK_NAME_ALIASES` literal is a
      subset of the registry-derived match sets, generator output current. `benchmark_validation._best_prediction_match`
      now uses `compound_keys.match_norms` (the literal table is deleted; the `SequenceMatcher` fallback remains and
      is the next thing to retire once a test shows it never fires on the panel).
- [x] 3b `data/keys/papers.yml` (246 DOIs) generated by `scripts/generators/build_paper_registry.py`; `paper_id` =
      intake id (161) > SLR id (24) > dossier stem > DOI slug (61); `record_ids` lists every citing record per file.
      `src/paper_keys.py` (`normalise_doi` handles Wiley DOIs with parentheses and markdown emphasis; `for_doi`).
      `tests/unit/test_paper_keys.py`: DOI set under `data/` == registry, ids unique, generator current.
      Labelled `autogenerated_unverified`; it asserts which records share a paper, not that a DOI is right.
- [ ] 3b follow-ups: only 5 of the 26 dossiers that print a DOI match a cited DOI (44 print none) — see the
      MISS list in the session log; `data/articles/manifest.json` still to do (needs the dossier identity tables).
- [ ] 3c `keys/chemistry_families.json`, `keys/matrix_families.json`, `keys/reaction_families.json` as closed enums
      (not started; the intake registry's 86 `matrix_family` values and the 5 chemistry-family aliases are unchanged).
- [ ] 3d rewrite record references to `paper_id` / compound `id`, retire the remaining alias tables
      (`external_validation._COMPOUND_ALIASES`, `experiment_value._alias_keys`, `sensory.alias_map`) — all of their
      spellings are already covered by the registry (tested), so this is now mechanical.
- [ ] `keys/papers.yml`: one `paper_id` per normalised DOI (227), built by script from the intake
      registry + SLR matrix + dossier IDENTITY tables; case-fold DOIs. Rewrite every record's paper
      reference to `paper_id`; keep `citation_repair`/`doi_repair` history in one `repairs.yml`,
      not inline in six files. Rename dossiers to `<paper_id>.md`.
- [ ] `keys/compounds.yml`: seeded from `matrix_decision_panel.json` (already has `aliases`) +
      species YAMLs + `henry_constants` names + every alias table in code + `test_regression.py:53-55`.
      InChIKey from RDKit on the existing SMILES; fail the build on collisions. Then: benchmark
      `measured_volatiles` keyed by compound id (name kept as `display_name`); the five alias tables
      and the `SequenceMatcher` fallback replaced by one `resolve_compound()` over this file.
      **Expect a handful of benchmark rows to change match status** where the fuzzy matcher was
      silently merging or splitting compounds. Record every such change in the PR.
- [ ] `keys/reaction_families.json`: the explicit crosswalk that `barrier_constants.py:740-748,807-813`
      currently hard-codes; the kJ/kcal collision in `arrhenius_params.yml`'s header gets resolved
      or explicitly quarantined here.
- [ ] `keys/chemistry_families.json`, `keys/matrix_families.json`, `keys/process_states.json` as closed
      enums; `family_identifier_contract.py` goes from reporting unknown ids to rejecting them.
      Collapse the intake registry's 86 `matrix_family` values onto the 8 canonical + a free-text `matrix_note`.
- [ ] `tests/unit/test_data_referential_integrity.py`: every `paper_id`, compound id, family id,
      matrix id in every file resolves; every odour threshold has a `paper_id`.
- **Exit:** no free-text join anywhere in `src/`; headline guards re-pinned with the diff explained.

### Phase 4 — Schemas and merges (~3–4 sessions) — STARTED 2026-09-01
- [x] `data/schemas/benchmark.schema.json` (JSON Schema 2020-12) written from the 56 committed payloads: required core
      (`benchmark_id`, `precursors{concentration_mM}`, `conditions{temp_C,ph,water_activity,time_min}`, `metadata{tier,family,
      execution_path,notes}`), closed enums for tier / execution_path / evidence_class / protein_type, `measured_volatiles`
      XOR `reference_volatiles`, `holdout_targets` only with `reference_volatiles`, every volatile entry carries `conc_ppb`
      or a self-explaining `conc_ppb_deliberately_absent`; optional `ph_basis`. Provenance blocks stay free-form.
- [x] `scripts/ci/schema_gate.py` (pyyaml + jsonschema, no science stack) — schema + `benchmark_id == stem` + every
      volatile spelling resolves in `compounds.yml` + hold-out bundles labelled. Wired into the CI gates tier. All 56 files pass.
- [ ] `data/schemas/benchmark.schema.json`: one core, `conditions` with typed optional `buffer`,
      `measured_volatiles` **xor** `reference_volatiles` with the semantics documented, compound
      keys from `keys/compounds.yml`, `ph_basis: initial | final`. `scripts/ci/schema_gate.py` (stdlib + pyyaml).
- [ ] One record envelope (§2 P3) for `references/*`, `priors/*`, `constants/*`; split
      `computational_priors.json` by section; merge `extrusion_damage` → `safety`, the two MLP files,
      `family_ingestion_plan` → `chemistry_families`, `process_gap_registry` + `structural_gaps` → `gaps.json`,
      measured rows of `retention_reference_payloads` → `binding_constants.yml`.
- [ ] Unify the four `*_contract.json` shapes (`protocol_id` vs `package_id`). Make
      `matrix_experiment_intake.py` read its enums from the schema file.
- [ ] Confidence vocabulary: pick one (decision 4); dedupe `low_medium`/`medium_low`; add the tier
      field to the 19 files lacking it, defaulting honestly (`unlabelled`).
- [ ] Collapse `process_state_calibrations.kind` (20 values / 26 rows) to a real vocabulary.
- **Exit:** every file in `data/` validates against a schema in CI; `data/lit/README.md`'s file map is generated.

### Phase 5 — Docs and contracts (~1 session) — PARTIAL 2026-09-01
- [x] `data/README.md` generated by `scripts/generators/build_data_readme.py` from `git ls-files data` + a hand-written
      description map; `tests/unit/test_data_readme.py` fails on an undescribed file or a ghost. `data/lit/README.md`
      now points at it and keeps only the ingestion workflow.
- [x] `agents.md` / `CONTRIBUTING.md`: layout, data-access rules (`data_paths`, `data_access`, keys, schemas, gates), pitfalls.
- [x] `tasks/lessons.md`: six patterns from this restructure.
- [ ] Remaining: stale prose (`timeseries/*.yml` headers, `external_validation.py:7` "19 benchmarks"); confidence-vocabulary
      decision (§5.4); `data/Gemini_Deep_Research` rename; the seven owner decisions in §5 that are still open.
- [ ] `data/README.md` generated from `data_paths.py` + schema titles; CI fails on ghosts/omissions.
- [ ] Correct `agents.md` / `CONTRIBUTING.md`: `data/` rules, layout map, tier vocabulary, the
      "curated pathways live in data/" exception (goes away with Phase 1).
- [ ] Fix stale prose: `timeseries/*.yml` "NOTHING IS WIRED" headers, `external_validation.py:7`
      "19 benchmarks", `COMPUTATIONAL_GAP_RUNBOOK.md` (lane deleted in Phase 1).
- [ ] `tasks/lessons.md`: the pattern "gitignored dirs inside `data/` are invisible to audits".

Total: roughly 11–16 sessions. Recommended cut if time is short: Phases 0–2 in full (high
value, ~zero risk), Phase 3 compounds + papers only, Phase 4 the benchmark schema only.

---

## 5. Decisions needed from the owner — RESOLVED 2026-09-01 (owner round 2)

1. **Mocked protein-source registry → withdrawn.** Trace: on a default run none of its numbers apply (only
   `--protein-source X`, `--protein-type myco`, or a quarantined benchmark carrying `protein_source` reach them); the
   kinetic-core B4 lane already refuses the file by policy; four of five fields are unit-less invented scores and
   `aa_composition` rescales the user's own input rather than seeding an endogenous pool. Deleted with its code
   paths. If plant-source differentiation is wanted later, it must come from measured inputs: a per-source amino-acid
   composition with units and a citation that seeds an endogenous precursor pool, and lipoxygenase presence reported
   as a cited fact, not a multiplier.
2. **Synthetic benchmark pairs → one per matrix.** All four files were frozen model output (labelled as such); the two
   `ProtocolPilot2026` files were byte-identical twins of `Internal2026` under a `measured_volatiles` label, and the
   "protocol pilot vs internal" hexanal/nonanal closure lane compared the two (always 1.0). Twins and lane deleted;
   `Internal2026` kept, labelled `diagnostic_only`. Still inside the scored MC panel (the README discloses the
   synthetic share); moving them out of the panel entirely is the next honest step and is listed under follow-ups.
3. **Confidence vocabulary → keep the load-bearing three, drop the rest.** `source_status` and `provenance_tier` feed
   blocking citation-gate checks and the census; `confidence_tier` moves one number (furanone lane). `uncertainty_posture`
   (175 records), `validated_status` (58) and `value_basis` (15) are read by nothing but report echoes: deleted, no
   output changes. `low_medium` → `medium_low` (undocumented synonym). Furanone tier map extended to every documented tier.
4. **Keys → continue (safer long term):** rewrite record references onto registry ids and retire the remaining alias
   tables and the fuzzy fallback.
5. **Ledgers → `results/literature/`** with the citation gate and census extended to scan them there.
6. **Research corpus → keep, rename** to signal it is not provenance.

### Round-2 execution log (2026-09-02)
- [x] Mocked protein-source registry withdrawn with its code paths (`literature_runtime`, `matrix_correction`, `recommend`,
      `pipeline`, `run_pipeline --protein-source`, `comparative_cli`, `benchmark_validation`); B4 lane policy reworded to
      "withdrawn"; census 102/80/80 → 87/65/65 re-pinned in the guards, README and the model card. Output-neutral on default runs.
- [x] Synthetic pairs collapsed: `*_ProtocolPilot2026.json` deleted (byte-identical twins), `*_Internal2026.json` kept and
      labelled `diagnostic_only`; the "protocol pilot vs internal" hexanal/nonanal ratio lane (always 1.0) deleted with its
      artifacts; the campaign arm no longer claims `calibration_passed`. Panel 23 → 21, MC panel 20 → 18 benchmarks /
      47 → 35 matched rows (the 12 rows that left were the model scored against a copy of its own output);
      `honest_literature_coverage` 4/13 unchanged. Evidence-role split 14/5/4 → 14/5/2.
- [x] Vocabulary: `uncertainty_posture` (175), `validated_status` (58) deleted from data and their echo sites; `low_medium` →
      `medium_low` (7); furanone tier map covers all five documented tiers and fails on an unknown one.
- [x] Ledgers moved to `results/literature/` (tracked); citation gate, census and the paper registry scan them there.
- [x] `data/Gemini_Deep_Research` → `data/research_corpus`; the citation gate's digest regex catches both spellings.
- [x] Artifacts regenerated (MC panel, summaries, matrix artifacts, figures, ranking, model card).
- [x] 3d: the five alias sites resolve through `compound_keys`; the token-overlap and `difflib` fallbacks in the
      benchmark matcher and the substring fallback in `experiment_value.lookup_spec` are gone. The 35-row MC pin held,
      proving no scored row depended on them.
- [~] 3c (STAGED, NOT COMMITTED, end of 2026-09-02): the 11 intake records using the five chemistry-family aliases were
      rewritten to canonical ids; `literature_family_registry._CANONICAL_FAMILY_ALIASES` deleted;
      `tests/unit/test_chemistry_family_keys.py` closes the axis (16 ids; intake + process-gap registries must use them).
      Targeted tests: 77 passed. FULL SUITES NOT RUN. Before committing: run both tiers.
      Open finding: regenerating `results/validation/literature_backlog.*` drops 6 intake rows that are NOT touched by
      this change (`fadel_2015_mft_retention`, `farmer_1991_alkyl_thiazoles`, `mottram_2001_bmfd_retention`,
      `nishimura_abe_2024`, `siripitakpong_2026_fft_retention`, `wang_2023_mft_retention`; 203 → 197 encoded
      references) — i.e. the tracked artifact is stale relative to the registry for a reason that predates 3c
      (`generate_literature_backlog` was not in the round-2 regeneration list). Diagnose why those six leave before
      regenerating; the four regenerated files were reverted to HEAD for now.
      Matrix-family axis: 86 free-text `matrix_family` values in the intake registry describe paper matrices, not the
      8 product matrices — needs a per-record mapping (rename to `matrix_context` + optional canonical `matrix_family`),
      not a mechanical fix. Reaction-family crosswalk deferred: the legacy engine that owns those keys is being retired.
- [ ] Next session: (1) finish/commit 3c as above; (2) roadmap item 1 — retire the SMIRKS/Hammond engine and the
      volatile-budget layer behind one typed pipeline contract. A read-only footprint map was requested from a subagent
      at the end of 2026-09-02; if its report is not in the transcript, redo the map (entry points, keep/retire/seam per
      module, data files, tests, what the MC envelope samples).

### Product decision (owner, 2026-09-02)
One tool, one engine. Retire the legacy SMIRKS/Hammond path (`src/smirks_engine.py`, `src/reaction_templates.py`
enumeration, `FAST_BARRIERS`, the fixed volatile-budget renormalisation) and keep the kinetic core as the single
source of concentrations; ranking/directional output by default, absolute concentrations only after the user
calibrates on their own data (`maillard ingest`). Architecture follow-ups are listed in the session record.

### Product roadmap (owner-approved 2026-09-02; item 6 reserved)
1. **One pipeline with a typed contract.** formulation + process recipe → concentrations (kinetic core) → matrix and
   headspace physics → observables → sensory axes → ranking + next experiment, as one callable with typed input/output
   and a stable JSON result; CLI, HTML report and notebooks consume it. Retire the legacy SMIRKS/Hammond engine and the
   volatile-budget renormalisation as part of this step.
2. **Process as a trajectory.** A first-class process-recipe schema (time, temperature, moisture, pH along the path;
   extrusion zones, frying, retorting) replacing scalar temp/time; Cantera becomes the one integrator if a parity test
   against the current propagator passes.
3. **Matrix by measurable quantities.** Protein source/grade, moisture and a_w, pH, lipid content, free amino-acid pool,
   thiol content; cited defaults labelled as defaults. No lookup multipliers.
4. **Calibrate-on-your-data as the core loop.** `maillard ingest` fits the few matrix parameters on the user's GC-MS table
   with the hold-out firewall kept; calibrations versioned per matrix; every prediction names the calibration it rests on.
5. **Sensory targets as the headline output.** Source the odour thresholds (1 of 26 carries a source today); add reference
   profiles of cooked meat and commercial analogues; report distance to target on the meaty / off-note axes.
6. *(Owner reservation: keep the absolute-ppb Monte-Carlo intervals as they are.)* Add rank-stability / P(A beats B)
   alongside them rather than instead of them.
7. **Next-experiment recommendation as the second default output**, with the protocol sheet generated.
8. **Cut the surface.** One CLI (`predict`, `compare`, `rank`, `ingest`, `next-experiment`), one report, one `validate`
   command that regenerates the evidence dashboard; guards pin contracts, not numbers. Notebook-first.
Not to build: quantum chemistry, ML potentials, generic reaction enumeration, absolute claims without user calibration.
Sequence: 1 → 8 → 4 (then 2, 3, 5, 7). Data lever: make the PPI/SPI protocol trivially reproducible.

### Engine retirement — footprint map (2026-09-03, read-only survey)
The legacy path is not a side lane; it is the whole validation harness:
- **Two front doors, two engines.** `scripts/run_pipeline.py` → `MaillardPipeline` → `SmirksEngine.enumerate` →
  `ResultsDB.get_best_barrier` (FAST_BARRIERS + refinement patches) → `Recommender.predict_from_steps` →
  `_apply_output_projection` (the fixed volatile budget, `recommend.py:947-950`) → `predicted_ppb` → `reporting`.
  `scripts/maillard.py --lane core` → `comparative_cli.predict_core` → `kinetic_core.engine.predict` (own 15-step
  mass-action network in `kinetic_core/network.py`, parameters from `results/validation/kinetic_core_b*_fit_report.json`,
  imports nothing from `src/` but `data_paths`) → `report_html` / `explain`. The core never sees SMIRKS steps.
- **Everything that validates runs on the legacy engine:** `benchmark_validation.py` (2852 lines, 36 test files)
  builds `SmirksEngine` + `predict_from_steps` directly; the Monte-Carlo envelope (`uncertainty_propagation.py`)
  samples 14 FAST family barrier offsets via `BARRIER_OFFSETS` and its artifact `prediction_uncertainty.json` feeds
  `pipeline` envelopes, `experiment_value.rank` (used by the CORE `rank` verb), external/cross-validation, ingest,
  matrix recalibration, and every headline-guard pin (panel 21, MC 18/35, 4/13, 1.44 dex). The core's own scoring is
  `scripts/generators/generate_cutover_final_exam.py` only.
- **Legacy-only modules (retire):** smirks_engine, reaction_templates, curated_pathways, barrier_constants
  (2 seams: HEME_*, arrhenius_rate_constant → conditions/recommend, both legacy), kinetics, results_db,
  family_sensitivity, cantera_export, projection, projection_utils, bayesian_optimizer, pre_processor,
  trunk_kinetics, and the harness modules pipeline, recommend, conditions, benchmark_validation, external_validation,
  cross_validation, reporting, presentation, usability_reports, sensory, safety, lipid_oxidation, headspace,
  matrix_correction (none imported by kinetic_core / report_html).
- **Seams:** projection_metadata (a type, pulled through literature_runtime/benchmark_types), matrix_targets (used by
  the literature/family registries, not by the core, which has `matrix_oav.py`), uncertainty_propagation (rewrite).
- **Legacy-only data:** arrhenius_params.yml, refinement_surrogate_patches.json, four computational_priors sections,
  interventions.yml, formulation_grid.yml, results/maillard_results.db, precursors.yml (core takes names raw).
- **MC envelope after retirement:** sample the core's fitted (k, Ea) pairs and its declared assumption bands
  (Q10, lipid fraction, PV, furanic partition, pH drift) by re-integrating `predict()` per draw; drop the
  BARRIER_OFFSETS routing; re-point every consumer of `prediction_uncertainty.json`.
- **Three riskiest couplings:** (1) `prediction_uncertainty.json` shared by both worlds; (2) `benchmark_validation`
  is the legacy harness, not a neutral scorer — the panel must be re-based on the core before deletion;
  (3) `literature_runtime` + family/matrix registries sit on legacy types.
Consequence: "retire the second engine" = rebuild the validation harness, the envelope, the matrix/headspace
layer and the report on the core. Owner decision on staging requested 2026-09-03.

### Retirement, staged (option B) — execution log
- [x] B1a (2026-09-03, commit `0112107`; unit+scripts 1338, integration+scientific 418 green): `src/curated_pathways.py` (the hand-maintained mirror of the SMIRKS chemistry),
      `scripts/generate_reaction_network.py` (drew it), `scripts/reproduce_use_cases.py` (legacy demo),
      `tests/unit/test_data_integrity.py` (atom balance of the mirror) deleted; two soundness tests trimmed. The
      Bayesian optimiser is deferred: its helpers are used by the front-door UX tests and it dies with `pipeline.py`.
      The Cantera lane stays (roadmap item 2 may make it the integrator); `results_db`, `pre_processor`,
      `family_sensitivity`, `trunk_kinetics` stay until the harness moves.
      Also in this step: quarantined benchmarks referenced from the intake registry are labelled
      `artifact_status: quarantined` and `deep_research_tracker` no longer counts them as runtime-bound (the corrected
      path had made them look live); the `reaction_templates` ↔ `smirks_engine` import cycle is broken with a lazy
      accessor (it only resolved when `smirks_engine` happened to be imported first).
- [x] B2 core Monte-Carlo envelope — DESIGN SETTLED 2026-09-03; BUILT AND VERIFIED 2026-09-03 (execution log below the design):
      * `predict(spec, targets, parameters=...)` already accepts an operative-parameter override; add a `CoreDraw`
        (maillard overrides, Q10, lipid-fraction and PV scales, furanone partition Ea, pH drift) and
        `size_declared_bands=False` so a draw does not double-count the corner re-integrations; add
        `core_parameters(lane, frozen=...)` so draws are made in FIT-REPORT space (shared Ea scalars stay shared,
        `MEASURED_EA_OVERRIDES` / `NO_EA_KEYS` honoured).
      * Sampling: B1 identified pairs from their stderr (`k_mgo_mel`, `k_aa_frag`; the other two are unidentified →
        fixed); B3 `k_acr_dp`/`Ea_acr_dp` from `ci95_halfwidth`; B7 `sigma_log10` 0.0414; declared bands (Q10,
        lipid fraction, PV, furanone partition, K_aw, HS-SPME) sampled UNIFORM / log-uniform over the declared band,
        not as lognormal 90 % (they were declared as corner bands). Independent draws, `SeedSequence(0)`, n=200.
      * B8 (sulfur, most of the panel) carries NO uncertainty in its fit report. Decision: compute a Laplace
        (finite-difference Jacobian) covariance at the B8 optimum rather than invent a 0.4-dex prior; until it
        exists the artifact must mark sulfur priors `sampled: false, reason: no_uncertainty_in_fit_report`.
      * Panel on the core: union of the trust-loop panel, the 17 maillard_path hold-outs and the 4 matrix bundles,
        each row tagged `panel:`; `honest_literature_coverage` computed identically across the union; the 4 shared
        Hofmann rows carry `shared_with`. The two `*_Internal2026` snapshots are legacy-model output and leave the
        scored panel (regenerate later as core frozen predictions). Vocabulary gaps that refuse external rows
        (hydroxyacetaldehyde, mercapto-2-propanone, hydrogen sulfide as precursors) are recorded, not patched
        silently.
      * Artifact keeps every field consumers read (benchmarks[].{benchmark_id, bench_file, execution_path,
        protein_type, fitted_row, compounds[].{compound, measured_ppb, predicted_p5/p50/p95, inside_ci,
        ci_width_log10}}, summary.honest_literature_coverage.*, summary.signal_origin_split.*, summary counts);
        adds `predicted_point`, `lane`, `refused_compounds[]`, `priors[]` with `sampled`/`source`.
      * Delete rather than port: `uncertainty_propagation.py` sampling, `matrix_recalibration.py`, the MC sections of
        `external_validation.py` and `data_ingest.py`, `pipeline.build_formulation_uncertainty_envelopes`.
      * Effort: ~2 sessions (+1 for the B8 Laplace covariance and the precursor aliases).
      * EXECUTION LOG (2026-09-03). Files: `src/kinetic_core/engine.py` (`CoreDraw`, `frozen_parameters(lane)`,
        `core_parameters(lane, frozen=)`, `predict(..., draw=, size_declared_bands=)`, lipid lane takes q10/scales
        as arguments, `corners=False`); `src/kinetic_core/panel.py` (THE bundle→spec mapping, lifted out of the exam
        generator, plus panel membership and `quantification_family`); `src/kinetic_core/uncertainty.py` (priors
        table, sampler, `propagate_panel`, artifact + markdown); `scripts/generators/generate_core_prediction_uncertainty.py`;
        `tests/unit/test_kinetic_core_uncertainty.py` (33 tests). The exam generator now imports from `panel.py`
        and its JSON is identical to HEAD's code output (checked in a detached worktree).
      * DEVIATION FROM THE DESIGN, deliberate: the K_aw (±0.5 dex) and HS-SPME same-sample dispersion (10–23×)
        bands are facts about a HEADSPACE number. The first build multiplied every ppb row by them, which put a
        near-constant 1.33-dex width on SIDA-extracted thiols and HPLC/LC-MS acrylamide and HMF (nothing in those
        numbers passed through an air/water partition or a fibre). Now gated by the bundle's own
        `content_verification.quantification_class`: headspace → applied; extraction (SIDA, HPLC-UV, LC-MS/MS,
        internal-standard GC/MS of an extract) → not applied; undeclared → applied, and the row/summary say so.
        Rows: headspace 2, extraction 37, undeclared 10 (bundles: acrylamide_spi_extrusion_130C_ACSRef3,
        cys_ribose_140C_Hofmann1998, li_2026, liu_2023, pea/soy PratapSingh2021, Trikusuma2019, resconi_2023,
        Bolton1994, Cerny2008 — declaring their `quantification_class` is a curated-data task for B4).
      * RESULT `results/validation/core_prediction_uncertainty.{json,md}` (n=200, seed 0): 32/40 panel benchmarks
        answered, 49 rows, 18 refused; honest literature coverage **5/20 (0.25), median width 0.94 dex**,
        24 rows NOT EVALUABLE (sulfur rows quantified by extraction: the sulfur lane has no sampled uncertainty and
        does not respond to the sampled B1 pairs — width ~1e-6 dex), 5 fitted rows excluded. Per panel:
        external_matrix 2/4 (width 2.10), maillard_path_holdout 1/11 (19 n.e.), trust_loop 4/10 (5 n.e.).
        Mixed-population 7/49. Legacy artifact for comparison: 4/13 (unchanged, untouched until B4). The honest
        reading: the core's envelope is EVALUABLE on the acrylamide, trunk and lipid lanes only; sulfur needs the
        B8 Laplace covariance before any sulfur coverage number means anything.
      * Timing: n=50 single worker 471 s; n=200 with `--workers 4` 2290 s wall while both test tiers ran alongside
        (user CPU 29 min). Result is independent of the worker count (per-draw `SeedSequence.spawn`).
      * Verification: py_compile + pyflakes clean; five gates PASS; unit+scripts 1371, integration+scientific
        418 + 2 xfail. README, headline guards, legacy `prediction_uncertainty.json` and the model card untouched.
      * FOUND IN PASSING: the tracked `results/validation/cutover_final_exam.json` was already stale before B2 —
        HEAD's code gives 3/34 core rows within band, the artifact says 4/34 (four MFT/FFT sulfur rows moved with an
        earlier fit-report/data change). Regenerate at B4 with the re-pin.
      * B5 NOTE: `uncertainty.py` still imports `benchmark_evidence_role`, `benchmark_signal_origin`,
        `get_benchmark_metadata` from `src/benchmark_validation.py`; those three helpers must move to a core-side
        metadata module before the legacy harness is deleted.
- [x] B3 core scoring of the union panel — BUILT AND VERIFIED 2026-09-03.
      * Files: `src/kinetic_core/scoring.py` (`score_benchmark`, `score_panel`, markdown, artifact),
        `scripts/generators/generate_core_panel_scores.py` → `results/validation/core_panel_scores.{json,md}`;
        `src/kinetic_core/fit_targets.py` (which panel rows the CORE's fits read); `src/benchmark_metadata.py`
        (engine-neutral `get_benchmark_metadata`, `benchmark_signal_origin`, `benchmark_evidence_role`,
        `matrix_source_anchor`, `resolve_scale_thresholds`, moved out of the legacy harness, which re-exports them);
        `tests/unit/test_kinetic_core_scoring.py` (19 tests). Legacy vocabulary kept: per-bundle scale contract
        (`validation_contract.scale_thresholds`, else the global default), `scale_status`, `ranking_status`,
        `overall_status`, `strict_ready` with the same tier/execution-path eligibility, `blocking_issues`.
      * FOUND: the sulfur lane's fit (B2–B8, 62 rows / 23 free) read eight panel rows the fit-target index could not
        see — the Hofmann 1998 Table 1 pH-5 rows (ribose, glucose, fructose) AND the xylose pH-5 row that sits on the
        maillard_path HOLD-OUT panel, plus five fed-intermediate step rows. Declared in `fit_targets.CORE_SULFUR_FIT_ROWS`
        (test-checked against the bundles). Leverage 0.37/row < 0.5 → `global_low_leverage` by the repo's rule: rows
        stay in, annotated `in_core_fit`, and an `out_of_sample` split is printed beside `honest_literature`.
        The CORE's evidence role ignores the legacy engine's fit-recovery declarations (matrix observability factors
        fitted to Pratap-Singh/Trikusuma, the projection budget) and prints them as `legacy_evidence_role`: on the
        union panel the core has predictive 40 / fit_recovery 0 / internal_synthetic 0 (legacy 21-file panel: 14/5/2).
      * RESULT (`core_panel_scores.md`): 40 bundles, 32 scored, 49 rows, 18 refused. Within 3x: 8/49; honest
        literature 8/49 (median fold 11.1x, geometric mean 36x); **out-of-sample 4/40 (median 31x)**; rows the
        sulfur fit read 4/9. Per panel: trust loop 5/15, hold-out 3/30, external matrix 0/4. Per lane: acrylamide
        2/12, sulfur 6/29, lipid 0/7, trunk 0/1. **Contract passes / strict-ready: 1 — `thiamine_cys_glucose_120C_Bolton1994`**
        (MFT 13 → 17.4 ppb, 1.34x, under the bundle's own 3x/0.48-dex contract; PRIMARY, free_precursor, not a fit
        row). Legacy: 0/23 strict-ready. This number moves at B4; the guards are untouched here.
      * Envelope re-run with the B3 annotations (n=200, 6 workers, 2562 s under load): honest literature coverage
        5/20 → **7/25** because the five legacy fit-recovery rows are predictive on the core; out-of-sample 6/23
        (17 not evaluable); rows the core fit read 1/2 (+7 not evaluable). `results/validation/cutover_final_exam.*`
        regenerated (was stale; 3/34 within band; the README's "42.23x paired" prose is B4's).
      * Verification: py_compile + pyflakes clean; five gates PASS; integration+scientific 418 + 2 xfail;
        unit+scripts 1389.
- [x] B4 publish core numbers next to legacy for one release — DONE 2026-09-03.
      * README: the exam table/prose re-pinned to the regenerated artifact (3/34 within 3x; paired 24.78x vs 10.86x,
        "2.3x worse", was 42.23x / "3.9x worse" from the B7 scoring; trajectory row added for the B8 re-score;
        sulfur family 1/10 median 12.22x, ladder median 122.5x); NEW subsection "The core scored on its own" with the
        legacy-vs-core table (panel 21 vs 40; roles 14/5/2 vs 40/0/0; strict-ready 0/23 vs 1 of 40 = Bolton 1994;
        90 % CI coverage 4/13 vs 7 of 25 evaluable, 6 of 23 out-of-sample, 24 sulfur rows not evaluable) and the
        three disclosures (sulfur fit read 8 panel rows incl. the xylose hold-out; sulfur has no uncertainty;
        headspace bands gated by quantification method). AUDIT.md: one dated row.
      * Model card: `collect_core_panel()` (recomputed live, ~15 s) and a "kinetic core" validity-domain cell beside
        the legacy strict-ready cell; README block regenerated.
      * Guards: `tests/scientific/test_core_headline_guards.py` pins the core artifacts (scorecard recomputed live
        and diffed against the tracked file — the freshness check the exam lacked; envelope and exam read tracked)
        and the README/AUDIT sentences that quote them. Legacy guards untouched: both sets hold for this release.
      * Deliberately NOT done here: the calibration badge (line 6), the "When to trust" table and the legacy headline
        sections still describe the legacy lane — they go with it at B5.
      * Verification: five gates PASS; `generate_model_card.py --check` clean; unit+scripts 1389; integration+scientific
        426 + 2 xfail (the 8 core guards included).
- [ ] B5 delete the legacy engine and its harness — OWNER GO-AHEAD 2026-09-03 ("let's make sure we are moving
      towards a better tool"). Staged:
  - [x] B5a seams (2026-09-03): the surviving tool's import closure (kinetic_core, `scripts/maillard.py`,
        `comparative_cli`, `model_card`, `report_html`, `explain_compound`, `experiment_value`, the CI gates, the
        kept generators) no longer imports any legacy module — verified by a static import-graph check.
        `comparative_cli`: FAST/screening lane, `--lane`, `--absolute`, `--target-tag/--minimize-tag` deleted;
        `model_card`: legacy collectors (benchmark panel, free-precursor / matrix hold-outs, observability modes,
        directional) replaced by the core scorecard + envelope; the card now SAYS directional skill is not yet
        measured on the core (the retired lane's 24/36 does not transfer); `matrix_calibration_registry` no longer
        imports `matrix_correction` (classifier inlined); `benchmark_metadata` drops the legacy matrix-factor
        fit-recovery declarations; `experiment_value` (the `rank` verb) reads the core envelope
        (`data_paths.CORE_PREDICTION_UNCERTAINTY`); `legacy_evidence_role` removed from the scorecard/envelope.
        Tests pruned accordingly (screening-lane tests). Gates PASS, model-card check clean, unit+scripts 1374,
        integration+scientific 426 + 2 xfail.
  - [x] B5b delete (2026-09-03). Rule: a module goes if it is not reachable from the kept scripts and the front
        door (static import graph). DELETED: 62 src modules (the SMIRKS engine `smirks_engine`/`reaction_templates`/
        `barrier_constants`/`kinetics`/`results_db`/`cantera_export`/`trunk_kinetics`; the harness
        `benchmark_validation`/`pipeline`/`recommend`/`conditions`/`projection*`/`presentation`/`reporting`/
        `usability_reports`/`sensory`/`headspace`/`matrix_correction`/`uncertainty_propagation`/`external_validation`/
        `cross_validation`/`matrix_recalibration`/`data_ingest`/`matrix_experiment_intake`/`matrix_calibration_*`/
        `bayesian_optimizer`/`pre_processor`; the family/matrix promotion reporting layer; `literature_runtime`
        (4.5k lines, orphaned), `safety` (the second acrylamide/CML model), `precursor_resolver`, `formulation`,
        `thermo`, `protein_binding`, `experiment_request`, and other orphans), 68 scripts (run_pipeline,
        optimize_formulation, run_campaign, run_cantera_kinetics, ingest_results, explain_formulation,
        request_experiment, the cutover exam generator, every legacy refit/report generator), 124 test files.
        ARCHIVED (git mv, nothing reads them): 108 artifacts + `maillard_results.db` → `results/legacy_lane/` with a
        README. Kept in place: `trunk_rate_calibration_refit.*` (the core's B1 reads it), `holdout_frozen/`, the
        frozen cutover exam/prereg, every `kinetic_core_*` report. DATA: `formulation_grid.yml`, `sensory_tags.yml`,
        `refinement_surrogate_patches.json`, the three matrix-intake protocol files deleted; `data/README.md`
        regenerated; `data/keys/papers.yml` rebuilt (278 papers, was 285: the DOIs cited only from
        `barrier_constants.py`); registries' `target_runtime_modules` / `current_runtime_assets` / `artifacts`
        entries that named deleted modules rewritten as "retired 2026-09-03 (B5b): …" (36 strings), archived
        artifact pointers repointed to `results/legacy_lane/`.
        GATES re-based: `holdout_guard` check 2 = no core fit generator names the hold-out directory (path-like
        string constants, prose excluded), check 3 = `panel.panel_bundles` is a literal non-recursive glob;
        `fit_target_gate` reads `core_prediction_uncertainty.json`, registry check retired, constant-precision
        containers emptied (backlog: point it at `parameters_*.py`). `src/__init__.py` exports the core's front
        door. `tests/unit/test_cli_scripts.py` rewritten for the surviving entry points; the data census guard moved
        to `tests/scientific/test_data_headline_guards.py`. Verification: five gates PASS, model-card check clean,
        unit+scripts 569, integration+scientific 77 (646 tests survive of 1,800).
  - [x] B5c docs (2026-09-03): README rewritten around the one engine (296 lines, was 1,600; the old README, the old
        QUICKSTART, the two notebooks and the campaign spec are archived verbatim under `docs/history/`); new
        QUICKSTART; `scripts/docker_maillard.sh` rewritten (core-scores / core-envelope / model-card / gates /
        data-readme / keys / the four verbs; the 40-odd legacy lanes gone); AUDIT.md and VALIDATION_CONTRACT.md carry
        a dated status banner (the contract's body is history until its backlog rewrite); agents.md layout updated;
        engine refusal strings no longer name "the FAST lane". Verification: five gates PASS, model-card check
        clean, every README/QUICKSTART/agents.md link resolves, unit+scripts 569, integration+scientific 77.
- [x] B7 (2026-09-03, the product claim): `src/kinetic_core/directional.py` scores all 69 claims of
      `docs/validation/directional_claims_panel.yml` on the core through the front door; artifact
      `results/validation/core_directional_scores.{json,md}`. Rules stated in the module: refused arms and
      unrepresented observables are NOT EVALUABLE; identical predictions across a moved axis the lane has no term for
      are MISSES flagged `mechanism_absent`; the sulfur lane refusing for "no sulfur source" is a structural ZERO for
      a sulfur observable (CYS-01 agrees). RESULT: strictly independent **18 of 30** evaluable (60 %), 22 independent
      claims not evaluable (16 prose-only, 2,5-dimethylpyrazine / 2-pentylfuran / nonanal-from-oleate / H2S and
      hydroxyacetaldehyde precursors refused); excluding pH and a_w **15 of 22**; pH + a_w **3 of 8**; per axis:
      sugar 5/9, temperature 5/7, time 2/2, cysteine 3/3, pH 3/5 (caution, on the floor), a_w 0/3, ranking 0/1.
      All claims: 22/42. The retired lane's 24/36 is a different engine on a different evaluable subset.
      `directional_reliability` now reads the artifact (independent split) and derives its axis notes from it;
      `compare_core` carries a `reliability` block (axes moved, per-axis verdicts, governing = weakest) and the
      renderer prints it; the model card's second sentence and per-axis cells are back, from the artifact; guards pin
      the numbers and README quotes them; `docker_maillard.sh core-directional`. The legacy report is bannered.
- [x] B8 (2026-09-03): `scripts/generators/generate_kinetic_core_b8_laplace.py` → `results/validation/
      kinetic_core_b8_laplace_covariance.json`: Gauss-Newton covariance at the frozen B8 optimum (J by finite
      differences of B8's own 62-row residual vector, 23 free coordinates, ~80 s; reduced chi-square 1.03, i.e. the
      declared row sigmas are consistent with the fit). Identification rule: sigma below the kind's threshold AND
      not (optimum on a declared bound with sigma > 10 % of the band) → **18 of 23** identified; not sampled:
      `k_dimer_decay`, `k_thiol_decay` (flat), `Ea_decay_carbonyl_sink` (sigma 64 > 60), `Ea_decay_thiol_sink`
      (on its 102 kJ/mol ceiling, sigma 27), `ph_acid_yield_per_sink_event` (at 0, sigma 0.29). `uncertainty.py`
      draws the identified sub-space JOINTLY (Cholesky) per draw, clips to bounds, routes into the sulfur blocks and
      `CoreDraw.ph_drift`. ENVELOPE (n=200): honest literature coverage **11 of 44** (was 7/25; 5 not evaluable, was
      24), median width 1.07 dex (was 0.94), out-of-sample **8 of 35** (was 6/23); every lane sampled. Restored the
      two test files B5b had deleted for importing the exam generator (`test_kinetic_core_uncertainty`,
      `test_kinetic_core_b2_4`) minus their exam tests. CAVEAT recorded: the covariance is a local Gaussian at an
      optimum that sits on two bounds; a profile-likelihood or MCMC pass is the honest next step (backlog).
### Post-retirement steps (owner-approved 2026-09-03: "1 go ahead, 2/3 whatever is best for the tool, 4 go ahead, 5 agree, 6 yes, 7 go ahead")
- [x] Step 1 — fit rows declare their bundles: the sulfur fit's level rows in `generate_kinetic_core_b2_3_fit.py`
      carry `benchmark_id` / `benchmark_compound`; `generate_core_fit_targets.py` writes
      `results/validation/kinetic_core_b8_fit_targets.json` (`fit_target_ids`, `fit_leverage` 23/62, rows);
      `fit_target_index` and the fit-target gate glob `*_fit_targets.json`; `src/kinetic_core/fit_targets.py` reads
      the record and the hand table is gone; hold-out guard check 4 sees core fit records (14 records scanned).
- [x] Step 2 — the xylose row: `mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5` → `data/benchmarks/
      hofmann1998_xylose_cysteine_145C_20min_pH5.json` (trust loop, PRIMARY, `hold_out_history` block, declared
      in-fit). Chosen over dropping it from the fit because a refit is a new wave (B9, backlog). Panel: trust loop
      20, hold-out 16 (Wave U count re-pinned 17 → 16 with the reason).
- [x] Step 3 — the wave generators are frozen by hash manifest, not by moving them (ten test files import them by
      path): `build_wave_manifest.py` → `results/validation/wave_generators_manifest.json`;
      `tests/scientific/test_wave_generators_frozen.py`; the rule and the wave log in `scripts/generators/WAVES.md`.
- [x] Step 5 — axis refusal: `engine.axis_refusal` refuses a comparison whose arms differ in water activity (no lane
      has a term) or in pH on trunk / acrylamide / lipid (no pH term by declaration); `compare` returns
      `comparable: False` with the reason; `directional.py` scores such claims NOT EVALUABLE. Directional headline
      18/30 → **18/27** (25 independent claims not evaluable, 4 of them refused; pH on the sulfur lane 3/5; a_w 0/0).
- [x] Step 6 — bring your own measurement: `src/kinetic_core/user_scoring.py` + `maillard score DOC` (`--template`):
      scores measured ppb like the panel scorecard (fold, 3x band, interval, refusals) and writes a bundle-shaped
      record with provenance under `results/user/` (never `data/`); explicitly scoring, not calibration.
- [x] Step 4 — slice profiles of the B8 objective (`generate_kinetic_core_b8_profile.py` → `results/validation/
      kinetic_core_b8_profile.{json,md}`, 161 evaluations, 340 s): verdicts {'asymmetric': 6, 'bound_limited': 9, 'quadratic': 8}.
      quadratic: k_tdp_fur, k_nf_mft, k_mgo_mp, k_glc_ha, k_fft_decay, k_pent_caramel, ph_arp_secondary_ammonium_pKa, k_arp_tdp. asymmetric: k_nf_mp3p, k_cys_thermal, k_thiolate_loss, k_osone_decay, k_ttca_deg, Ea_decay_carbonyl_sink.
      bound_limited (the declared band ends within 2 sigma): k_dimer_fft, k_fur_decay, ph_acid_yield_per_sink_event, k_fur_fft, k_dimer_decay, k_thiol_decay, Ea_decay_thiol_sink, k_dimer_mft, k_arp_dpo.
      READING: the Laplace picture holds locally on 8 coordinates and is one-sided on 6; on 9 the DECLARED BANDS,
      not the data, end the slice — the search ranges are active constraints of the B8 optimum. The slice rise at
      2 marginal sigma is far above 2 on the quadratic coordinates (e.g. k_tdp_fur 17.8), i.e. the marginal sigmas
      are inflated by strong correlations; the joint draw in the envelope is the right object, the marginal table
      is not. The thiol-sink ceiling (102 kJ/mol) is one of the active bounds: whether it is Gigl 2021's measured
      range or a search convenience is the first question for B9, together with the fit/validate split.
- [x] Step 7 — B6 test audit: `tasks/test_audit.md` (49 files, per-file verdicts). Executed: the public-API smoke test
      deleted; ten "markdown mentions phrase" prose tests deleted; the seven frozen-wave contract files labelled as
      regression records. Left in the audit's own backlog: merge the four generator-text checks into the freeze guard,
      rewrite `test_v1_reports`' five substring checks, collapse the fourteen literature-side files into one
      parametrised builder test, fix the coverage/fork configuration.
- [x] **B9 RUN 2026-09-03, committed `1067358`** (`generate_kinetic_core_b9_fit.py`, prereg `kinetic_core_b9_prereg.md`): B8 minus the
      eight Hofmann T1 level rows, 54 rows, same free set / bands / weighting / protocol. Both starts converge to cost
      18.74 (B8's vector scores 18.87 on the 54 rows); active bounds: k_dimer_mft, k_dimer_fft, Ea_decay_thiol_sink
      (ceiling), acid yield (floor). Laplace at B9: 20/23 identified, chi2_red 1.21. PREREG CHECK: hold-out 4/30
      (not worse), returned bundles **1 of 8** within 3x (expected 2-4; falsifier was 0/8 AND hold-out -2) → B9 SHIPS.
      THE FINDING: glucose and fructose MFT predict ZERO without the level rows — the hexose→thiol route had no
      step-level support and was carried by the very rows it was scored on. Scorecard 8/49 → **6/49**, out-of-sample
      4/40 → **5/48** (only the C2+C3 pH-5 row is a fit row now), median fold 11x → 50x; directional 18/27 → **17/27**.
      Xylose returned to the hold-out (17 bundles). B8's fit-target record removed (superseded). Envelope on B9 (n=200):
      **10/44** evaluable literature rows (was 11/44; median width 1.17 dex, was 1.07), out of sample **10/43** (was 8/35).
      B9 slice profile: 4 quadratic / 9 asymmetric / 3 flat / 7 bound-limited — the declared bands are still active
      constraints and the hexose coordinates are now flat (unconstrained by step-level data). The bands question and a
      hexose→thiol step-level measurement are the next science items.
      Original definition: the sulfur refit that drops the eight Hofmann 1998 Table 1 LEVEL rows (ribose / glucose /
      fructose / xylose pH 5, FFT and MFT) from the objective and keeps the step-level rows (fed intermediates T3 /
      T4 / T10, conversions, ratios, Kang / Zhou / Whitfield / Cerny). Rule going forward: rate constants,
      activation energies, fed-intermediate yields, conversions and within-study ratios are FIT evidence; end-to-end
      concentrations in full precursor systems are VALIDATION evidence. B9 returns all four Hofmann pH-5 bundles to
      validation (xylose to the hold-out). Check before shipping: absolute scale must stay identified from the fed
      rows (yields per mmol fed) — the profile of the new optimum decides.

- [ ] B6 (owner request 2026-09-03, not urgent; AFTER B5 so the audit covers the tests that survive): test-suite
      design audit — are the tests well designed, which files should be restructured for coverage, and which tests
      assert nothing (tautologies, pinned-artifact echoes, contract tests without a failure mode). Deliverable: a
      per-file verdict (keep / merge / rewrite / delete) with the coverage gaps named, then execute it.

### Original decision list (kept for the record)

1. **`protein_source_registry.json`** is self-labelled `value_basis: mocked_placeholder` and
   `source_status: no_verifiable_source`, yet feeds four `kinetic_core` modules. Options: (a) gate
   it — multipliers default to 1.0 and the report says "matrix multipliers unsourced" (recommended;
   honest, cheap); (b) keep as-is with the warning; (c) source the 14 profiles (a research task,
   not a cleaning task). *(Withdrawn 2026-09-01 — see resolved item 1 above.)*
2. **`data/qm/phase33_*` / `phase35_*`** (27 values, no verifiable source, zero runtime readers,
   pinned by `test_honest_headline_guards` and `citation_gate` check 5). Keep under
   `qm/unverified/` with the label, or delete and drop the pins. Recommended: delete; the
   AUDIT.md narrative already preserves the finding.
3. **Research corpus** (`Gemini_Deep_Research/raw` + dumps, ~2 MB tracked): keep in-repo as the
   forensic corpus `citation_gate` scans (recommended), or move off-repo and relax the gate.
4. **Confidence vocabulary**: migrate data to the four labels `agents.md` mandates, or amend
   `agents.md` to the three fields actually in use. Recommended: amend the docs; the data's
   `provenance_tier` + `uncertainty_posture` carry more information than the four labels.
5. **Scope**: stop after Phase 2, or commit to Phase 3. Phase 3 is where maintainability actually
   comes from; without it the alias tables regrow.
6. **The synthetic `Internal2026` / `ProtocolPilot2026` pairs**: keep one per matrix, relabel
   `evidence_class: synthetic_snapshot` and `reference_volatiles`, retire the other and repoint the
   ~38 test references. Or keep both and only fix the `measured_volatiles` label.
7. **Layout depth**: flatten `data/lit/*` into `data/{keys,constants,priors,references,…}` as in §3,
   or keep `lit/` as a parent. Flattening is recommended: "lit" currently means "everything".
8. **Script deletions** (Appendix C) need explicit confirmation per `agents.md`.

---

## 7. Improvement backlog (standing owner request 2026-09-03: "suggest ideas that improve the codebase and add them here")


Not scheduled; each is a self-contained change with its motivation. Strike or promote as decided.

**Done 2026-09-03 (backlog pass, commits after `66c6fe3`):** fit-target declarations (step 1), the hold-out the fit read
(B9), README numbers via the model card + guards (B4/B5c), the exam scorer (frozen history), shared helpers
(`src/provenance.py`, `src/report_format.py`, `src/artifact_io.write_artifact`; the three headline artifacts now
carry a `provenance` block with git + input hashes), `data_paths.bundle_path`, CLI headline numbers read from the
artifacts (`comparative_cli.core_caveat`, `maillard.py` help), dead code (`benchmark_evidence_role`,
`SCREENING_CAVEAT`, the model card's unproduced hold-out branch, five `_load_json` wrappers, two dead
`data_paths` constants), stale pointers to `benchmark_validation` / `matrix_calibration_registry` /
`get_benchmark_files`.

**Done 2026-09-03 (backlog pass 7, the last):** one evidence-role vocabulary in code (`fit_targets.EVIDENCE_ROLES`:
fit_recovery / internal_synthetic / external_holdout / predictive; the bundle-level `evidence_class` stays as the data-side
marker the gates check, so no hold-out bundle was rewritten); `trust` now also requires the 95 % Wilson lower bound above
0.50 (pH 4/5 → caution; thresholds not moved); the deep-research tracker writes `deep_research_gap_analysis.json` beside the
ledger instead of overwriting it; five caller-less scripts deleted (`calibrate_lipid_oxidation`, `check_env`,
`check_references`, `find_all_doi_gaps`, `sync_backlog`); five more bundles verified from PubMed Central full text
(Pratap-Singh ×2, Li 2026, Hernandez/Resconi 2023, Ma 2024): classes declared, `content_verification` blocks with the
quoted method sentences and table values. Envelope after the five verifications: ACSRef3 (isotope-dilution LC-MS/MS) loses its band-only interval, 9/43 → **9/42** literature rows,
out of sample **9/41**, rows by family headspace 8 / extraction 39 / undeclared 2. **Thiol-sink ceiling: KEPT.** The (7, 102) kJ/mol band is Gigl 2021's covalent
capture range (`parameters_sulfur.py`, B8 prereg T1), a measurement, not a search bound; that B9 presses it (bound-limited
in the profile) is the finding, recorded, and the Laplace leaves the coordinate at the bound rather than sampling across a
measured limit. A profile likelihood over 23 coordinates (~20 h) is not worth the hold-out it would not change.

**Done 2026-09-03 (backlog pass 3):** `results/README.md` generated by `build_results_readme.py` (exact / glob / directory
entries, ghost + orphan checks, `--check`, `tests/unit/test_results_readme.py`); orphans archived under
`results/legacy_lane/` (`experiment_requests/`, `active_learning_requests.json`, `family_coverage.png`; the archive
README is now tracked); `scripts/ci/artifact_freshness_gate.py` (sixth gate: regenerate-and-compare for the
scorecard / directional / fit-target record modulo volatile keys, every `--check` runner, and provenance input
hashes for every artifact that carries a block, the envelope included).

- [x] **Envelope cost — RESOLVED 2026-09-03 (backlog pass 6): native arm64 container.** The container now follows the host
      architecture (`scripts/docker_maillard.sh`; the amd64 one stays reachable via `MAILLARD_PLATFORM=linux/amd64`).
      Natively: tier 1 196 s → 83 s, 12 envelope draws 127 s → 46 s with six workers (serial pass of point predictions
      remains), full envelope ~5 min instead of 35-43. Cross-architecture reproducibility measured at ≤ 7e-8 relative on the
      scorecard; `FLOAT_REL_TOL = 1e-6` (+ `FLOAT_ABS_TOL = 1e-9`) in `src/provenance.py` and the scorecard guard. Two
      reproducibility leaks fixed on the way: the scorecard's blocking string printed a 1e10 ratio to three decimals (now 4
      significant figures), and the directional scorer took 1e-31 µg/L as a number (now floored to zero below 1 pg/L; one
      noise-ratio MISS became NOT EVALUABLE: 17/27 → 17/26). Original finding follows.
- [ ] (history) **Envelope cost: the Docker container serialises worker processes (found 2026-09-03).** Measured in the amd64
      container (Rosetta on Apple Silicon): one sulfur predict is 0.17-0.33 s (LSODA, Python RHS, no file I/O once the
      engine caches its fit reports); `propagate_panel` at 12 draws takes 127 s with 1 worker, 127 s with 6 forked
      workers, 134 s with 6 spawned workers. Per worker: 64 s wall for 10.8 s CPU; two workers take 2x, three take 3x.
      A probe shows the same for pure-Python loops once numpy is imported in the parent, while a numpy-free pool
      parallelises (12 tasks: 1.4 s -> 0.4 s). Affinity is full, BLAS threads pinned to 1 change nothing. The cause is
      below the code (emulation layer); the fix is a native arm64 container: `MAILLARD_PLATFORM=linux/arm64` with its
      own container name and conda volume (script override added). Validate: probe parallelism, run both tiers, and the
      freshness gate must still match the tracked scorecard (rel 1e-9) before the default flips.
- [x] **`quantification_class` declared on every panel bundle and enforced (2026-09-03, backlog pass 5).** Three declared
      from text already in the repo (Bolton 1994 → `internal_standard_gcms_sim`, Liu 2023 → `hs_spme_gcmsms_external_standard_curve`,
      Trikusuma 2019 → `dynamic_headspace_gcmsms`; the marker list learned `dhs`); eleven set to an explicit `undeclared`
      with a `quantification_note` saying why (source not on disk, figure read-off, repo-internal derivation, model
      snapshot). `benchmark.schema.json` carries the enum; `schema_gate` fails a panel bundle whose class resolves to no
      family unless it says `undeclared`. Still open: fetch PMC8271896 (Pratap-Singh ×2) and the Li 2026 / Resconi
      method sentences to declare four more. Envelope effect: Bolton's row had only band width and is now NOT EVALUABLE:
      10/44 → **9/43** literature, 10/43 → **9/42** out of sample, rows by family headspace 4 / extraction 38 / undeclared 7.

- [x] **Test coverage (closed 2026-09-03):** wave generators are frozen by policy and omitted from coverage on purpose
      (`WAVES.md`, hash manifest); the fork/multiprocessing coverage configuration is fixed (`uncertainty.py` 92 %).
- [x] **Shape-only tests (closed 2026-09-03 by the two test-audit passes):** prose tests deleted, wave contract files
      labelled frozen records, artifact-reading tests either recompute live (guards) or are contract tests by design.
- [ ] **Bundle verification coverage: 16 panel bundles have no `content_verification` block and only 31 of 278 papers
      carry an extraction dossier** (the quantification-class half of this item closed on 2026-09-03). The trust-loop
      bundles (Pratap-Singh, Trikusuma, Resconi, Bolton, Cerny, ACSRef3) and the matrix bundles are the ones without a
      primary-source pass; one pass per bundle, recorded in the bundle, closes it. Needs the PDFs (`data/articles/`).

## 6. Risks and guardrails

| Risk | Guardrail |
|---|---|
| A stale path after a move silently degrades the model (D5) | Phase 2 makes missing files raise; do Phase 2 before any `data/lit` rename |
| The three CI gates and `test_honest_headline_guards` hard-code globs and counts | Update in the same PR; the gates are stdlib-only and fast to re-run |
| Report output embeds paths and four tests assert them | Route through `data_paths`; treat the output strings as a versioned contract |
| Hold-out leakage: `holdout_guard` relies on `get_benchmark_files` being non-recursive | `benchmarks/calibration/` vs `benchmarks/holdout/` makes the separation physical; keep the guard's rglob check |
| Phase 3 changes compound matching | Diff `benchmark_summary.json` before/after; explain each changed row in the PR; re-pin headline guards last |
| Regeneration deletes provenance (D9) | Phase 0 idempotency test |
| Pytest writes into `data/` again | `data_readonly_gate.py` |
| Losing git history on moves | `git mv`, one directory per commit |

---

## Appendix A — Disposition of every tracked path under `data/` (current names)

| Path | Consumers (src/scripts/tests) | Disposition |
|---|---|---|
| `lit/arrhenius_params.yml` | 5/4/4 | keep → `constants/` |
| `lit/henry_constants.yml` | 2/0/2 | keep → `constants/` |
| `lit/binding_constants.yml` | 3/1/1 | keep → `constants/` (envelope template for the rest) |
| `lit/lipid_oxidation_calibration.json` | 4/2/1 | keep → `constants/` |
| `lit/calibration_offsets.json` | 0 readers | delete (generated, unwired, fit targets gone) |
| `lit/computational_priors.json` | 10/2/4 | split → `priors/*.json` |
| `lit/flavor_reference_payloads.json` | 5/2/2 | keep → `references/flavor.json`, one envelope |
| `lit/safety_reference_payloads.json` | 5/3/1 | keep → `references/safety.json` (+ extrusion_damage) |
| `lit/extrusion_damage_reference_payloads.json` | 1/0/0 | merge into safety |
| `lit/retention_reference_payloads.json` | 5/2/3 | measured rows → `binding_constants.yml`; withdrawn rows → `references/retention.json` (quarantine list) |
| `lit/process_state_calibrations.json` | 8/1/1 | keep → `references/process_state.json`; collapse `kind` |
| `lit/protein_source_registry.json` | 4/1/1 | decision 1 — withdrawn 2026-09-01 |
| `lit/matrix_decision_panel.json` | 3/0/2 | seed for `keys/compounds.yml` |
| `lit/matrix_family_coverage_registry.json` | 4/0/1 | → `keys/matrix_families.json` |
| `lit/chemistry_family_scope_registry.json` | 4/1/1 | → `keys/chemistry_families.json` |
| `lit/family_ingestion_plan.json` | 4/4/1 | merge into `keys/chemistry_families.json` |
| `lit/process_gap_registry.json` | 4/1/1 | merge → `intake/gaps.json` |
| `lit/benchmark_intake_registry.json` | 9/8/4 | keep → `intake/`; reference `paper_id`; drop derived columns |
| `lit/slr_incorporation_matrix.json` | 1/3/1 | generated → `results/literature/` |
| `lit/deep_research_backlog.json` | 2/5/1 | generated → `results/literature/` |
| `lit/dft_coverage_map.json` | 0 | delete (10 rows fold into chemistry_families) |
| `lit/computational_gap_closure_targets.json`, `_multistep_targets.json` | 0 | delete |
| `lit/geometry_benchmark_set.json` | 1/0/2 | → `qm/geometries.json` (+ ts_seed rows); rename misused `chemistry_family` field |
| `lit/ts_seed_benchmark_set.json` | 0 | merge then delete |
| `lit/reaction_benchmark_set.json` | 2/0/2 | → `qm/reaction_barriers.json` |
| `lit/mlp_candidate_registry.json`, `mlp_external_benchmark_evidence.json` | 1/0/2 each | merge → `qm/mlp_registry.json` |
| `lit/refinement_surrogate_patches.json` | 3/1/2 | keep (tombstone read by a retirement guard) → `constants/` |
| `lit/canonical_systems.json` | 0/0/1 | → `tests/fixtures/` |
| `lit/README.md` | — | regenerate |
| `lit/extraction_dossiers/*` (75) | cited as strings by ~15 files | → `dossiers/<paper_id>.md`; `k*`/`round*` → `dossiers/synthesis/` |
| `lit/timeseries/*` (4 + README) | 1/4/0 + gates | keep as-is; fix stale headers |
| `benchmarks/*.json` top-level (23) | glob | → `benchmarks/calibration/`; schema; decision 6 |
| `benchmarks/external_validation/*.json` (4) | generated | → `benchmarks/holdout/matrix/`; fix generator (D9) |
| `benchmarks/external_validation/maillard_path/*.json` (17) | hand-authored, frozen | → `benchmarks/holdout/free_precursor/` |
| `benchmarks/quarantined/`, `step_level_unreachable/` | READMEs cited | → `benchmarks/retired/{quarantined,unreachable}/` |
| `benchmarks/maillard_validation_benchmarks.md` | quoted by 5 src | → `docs/retired/` |
| `protocols/*_contract.json` (4) | 1–4 each | → `protocols/contracts/`, one shape |
| `protocols/matrix_experiment_intake_schema.json` | 2/1/0 | → `schemas/intake.schema.json`; enums move into it |
| `protocols/example_matrix_experiment_intake.yaml` | 2/0/2 | keep in `protocols/` |
| `protocols/{pea,soy}_iso_protocol_pilot_intake.yaml` | generated | → `results/validation/reproducibility/` |
| `protocols/external_validation/*.yaml` (4) | generated twins | → `results/validation/holdout_bundles/` |
| `protocols/requested_*.yaml` (5) | 0 | → `results/validation/experiment_requests/` |
| `species/{desirable,off_flavour,toxic}_targets.yml`, `sensory_tags.yml` | 2–4 src each | → `keys/compounds.yml` (Phase 3); until then unchanged |
| `species/precursors.yml` | 2/2/0 | → `keys/precursors.yml`, one section schema |
| `reactions/curated_pathways.py` | 0 (broken import) | delete |
| `reactions/reaction_families.yml` | existence check only | → `docs/reference/` |
| `rmg_extensions/**` (26) | 0 | delete |
| `qm/phase33_*`, `phase35_*` | loader + parse test | decision 2 |
| `qm/irc_validation_cases/**` | parse test only | with decision 2 |
| `geometries/**` (33 tracked) | 0 after Phase 1 | delete |
| `Gemini_Deep_Research/slr_family_*.md` (16) | 0 | → `results/literature/slr_family_reports/` |
| `Gemini_Deep_Research/{raw/*, 13 dumps, README}` | `citation_gate`, 3 scripts | → `research_corpus/` (decision 3) |
| `cli_examples/*`, `campaigns/*`, `ingest_templates/*` | docs + 1–2 tests | → `docs/examples/`, `docs/templates/` |
| `formulation_grid.yml`, `interventions.yml` | 2–3 src | keep (`interventions` → `constants/`) |
| untracked: `articles/` | prose citations | keep untracked + tracked `manifest.json` |
| untracked: `logs/`, `calibration_history/`, `geometries/{dft_runs,dft_checkpoints,multistep,xtb_inputs}` | 0 | delete |
| `results/maillard_results.db` | `pipeline.py:222,554`, `benchmark_validation.py:709` | keep; drop `ml_adoption_decisions` |
| `docs/notebooks/results/maillard_results.db` | 0 | delete |
| `results/validation/qm_barrier_provenance.json` | `uncertainty_propagation.py:513` | curated input → `data/qm/` |
| `docs/validation/directional_claims_panel.yml`, `directional_accuracy_report.md` | `directional_reliability.py:57,62` | leave in `docs/` (the report is the record by design); register in `data_paths.py` |

## Appendix B — Paper-identifier mismatches (sample of 14)

| PDF | Dossier | Benchmark id | Intake id | SLR id | Defect |
|---|---|---|---|---|---|
| `Kocada2016.pdf` | `kocadagli2016jafc` | — | — | — | shorter stem is the JAFC paper; documented swap trap |
| `Kocadagli2016.pdf` | `kocadagli2016foodchem` | — | — | — | text layer garbled |
| `Zhai2023.pdf` / `Zhai2023b.pdf` | `zhai2023foodchem` / `zhai2023jafc` | — | — | — | swapped relative to the K6a brief |
| `Goncouglu2016.pdf` | `goncuoglutas2017` | — | — | — | surname, year and no registry entry |
| `Goncouglu2026.pdf` | none | — | — | — | different paper, one digit away |
| `Kang2026(_supplementary).pdf` | `kang2026(_SI)` | — | — | — | case + suffix |
| `frankel1989.pdf` | none | — | `frankel_1982_…` | same | citation says 1983; three years |
| `Hamzalioglu2018.pdf` | `hamzalioglu2018` | — | `hamzalioglu_2026_…` | same | registry id is a different 2026 paper |
| `Zhou2023b.pdf` | `zhou2023` | — | — | — | dossier drops the `b` |
| `Xin2026b.pdf` | `Xin2026b` (capitalised) | — | — | — | dossier casing inconsistent |
| none | none | `pea_isolate_uht_140C_Trikusuma2019` | `trikusuma_2019` | two ids | benchmark with no PDF; SLR double id |
| `hofmann1998.pdf` | `hofmann1998_reconciliation` | 9 benchmarks, 2 casings | — | — | most-cited source has no registry id |
| `1-s2.0-S0308814622010068-main.pdf`, `jf0480290.pdf` | none | — | — | — | publisher-default names |

## Appendix C — Scripts proposed for deletion (owner confirmation required)

- `scripts/calibrate_barriers.py` — both fit targets quarantined/deleted 2026-08-26; output unwired.
- `scripts/generators/generate_rmg_inputs.py` — targets nonexistent, gitignored `data/reactions/rmg_validation_cases/`.
- `scripts/sync_backlog.py` **or** `scripts/ingest_deep_research_markdown.py` — three writers of one file; keep `deep_research_tracker.py`.
- `tests/scripts/{test_run_computational_gap_dft,test_run_computational_gap_xtb,test_react_ot_seed_coverage,test_import_react_ot_colab_artifacts,test_open_react_ot_colab}.py` — import scripts deleted in Phase 1.
- `tests/scientific/{test_refinement_campaign,test_refinement_watchlist,test_refinement_governance,test_offline_refinement_governance,test_mlp_assessment}.py`, `tests/unit/test_wave_r1_barrier_offset_retirement.py` — import `src` modules deleted in Phase 1 (check whether any surviving assertion is worth porting first).
