# Maillard Strategic Roadmap

## Active Plan — 2026-04-19 pyGSM integration & global barrier run

### Objetivos del usuario
1. Integrar pyGSM y conectarlo con el resto del pipeline (TS-guess → DFT).
2. Probar utilidad sobre los 7 targets de computational gap, no sólo los 3 ya validados.
3. Pulir todo y lanzar la corrida global para obtener barreras finales.

### Contexto duro (ya verificado)
- `pyGSM` no está en PyPI; instalable sólo desde `git+https://github.com/ZimmermanGroup/pyGSM.git`. Su `xtb` (xtb-python) transitivo no compila en el container actual → instalar con `--no-deps` y reutilizar `get_xtb_ase_calculator()` (tblite).
- Validación dura ya en producción (`_classify_ts_frequencies`, floor 50 cm⁻¹, `ts_validation.json` consumido por la promotion gate). Todo lo nuevo debe seguir pasando esa puerta.
- xTB sirve sólo como pathfinder; las magnitudes son inflated. DFT sigue siendo la única fuente de la altura final.
- 7 targets en scope: `hexanal_radical_quench`, `lysinoalanine_crosslink`, `aa_ring_open_dicarbonyl`, `quinone_cys_michael`, `pe_schiff_base`, `pe_amadori`, `asparagine_sugar_explicit_water_cluster`.
- Estado actual: hexanal y lysino flagueados `invalid_saddle`; aa_ring interrumpido; los otros 4 con seeds xTB pero sin DFT corrido.

### Fase 1 — pyGSM backend (núcleo)
- [ ] **1.1** Instalar pyGSM con `--no-deps` en el container y persistirlo en `environment.yml`/`Dockerfile`. Smoke test: import + ejemplo Müller-Brown.
- [ ] **1.2** Crear `src/gsm_backend.py` con `GSMRunner.run_de_gsm(reactant_xyz, product_xyz, *, charge, spin, n_nodes, max_iters, work_dir, timeout_s) -> GSMResult`. LOT delgado sobre `get_xtb_ase_calculator()`. Mover `_suppress_native_output` desde `scripts/diagnose_xtb_gsm_viability.py` a `src/_native_io.py`.
- [ ] **1.3** Conectar como cuarto escalón en el cascade de `src/ts_guess_refiner.py` (gobernado por `enable_gsm_fallback: bool = False`).
- [ ] **1.4** Conectar como fallback duro en `src/dft_refiner.py`: si `ts_validation.is_true_saddle == False` → relanzar Phase 3.2+3.3 con guess GSM. Persistir `gsm_attempt.json` en el checkpoint dir. Nunca usar la energía xTB como barrera final.
- [ ] **1.5** Exponer flags en `scripts/run_computational_gap_dft.py`: `--enable-gsm-fallback`, `--gsm-nodes`, `--gsm-max-iters`, `--gsm-timeout-s`. Counters en summary.

### Fase 2 — Validación contra los 3 targets ya estudiados
- [ ] **2.1** Tests unitarios en `tests/unit/test_gsm_backend.py` con potencial Müller-Brown analítico (saddle conocido).
- [ ] **2.2** Smoke real con xTB sobre `hexanal_radical_quench`, `lysinoalanine_crosslink`, `aa_ring_open_dicarbonyl`. Criterio Fase 2: ≥ 2/3 generan guess que pasa el gate xTB. Si no, **STOP & re-plan**.
- [ ] **2.3** Validación DFT light (`--fast`, sin IRC) sobre los que pasen 2.2: comprobar `is_true_saddle == True` y barrera plausible (5–80 kcal/mol).

### Fase 3 — Extensión a los 7 targets y diagnóstico de utilidad
- [ ] **3.1** Extender `scripts/diagnose_xtb_gsm_viability.py` a los 4 restantes; agregar a panel `xtb_gsm_viability_2026_04_19_full.json` y actualizar memory.
- [ ] **3.2** Smoke pyGSM sobre los 4 nuevos. Documentar honestamente cuáles fallan (Schiff base/Amadori probablemente requieren proton-transfer concertado).
- [ ] **3.3** Asignar tier por target: `tier_A` (entra en corrida global), `tier_B` (necesita ajuste), `tier_C` (intractable, etiquetado).

### Fase 4 — Pulido y corrida global
- [ ] **4.1** Limpieza de logging/dead code; docstrings sólo en superficie pública nueva.
- [ ] **4.2** `pytest tests/` completo en Docker, 0 fallos.
- [ ] **4.3** Borrar `_dft_execution.json` stale de hexanal/lysino. Lanzar `run_computational_gap_dft.py --execute --enable-gsm-fallback` sobre `tier_A` (+ `tier_B` con seguimiento). Recolectar en `results/computational_gap_refinement/computational_gap_dft_global_2026_04_19.json`.
- [ ] **4.4** Reporte: tabla por target (barrier, verdict, GSM sí/no, tier). Memory `/memories/repo/pygsm_global_run_2026_04_19.md`.

### Riesgos
- pyGSM puede tropezar con `asparagine_sugar_explicit_water_cluster` (muchos átomos).
- xTB con tblite es inestable en multiplicidades altas; saltar si aplica.
- Mapping de `forming_bond_pair` al sistema de constraints de pyGSM (DLC vs Cartesian).
- `--no-deps` requiere verificar que numpy/scipy/ase/networkx estén ya cubiertos por `environment.yml`.

### Política de ejecución
- Validar Fase 1 + Fase 2 antes de tocar la corrida global.
- Si Fase 2 < 2/3, parar y re-planear.
- Correcciones del usuario → `tasks/lessons.md`.

### Review
- _pendiente_

---

## Active Follow-up — 2026-04-17 DFT Ladder Bulletproofing Plan

- [ ] Preserve best-available geometry fallback for non-TS optimizations instead of raising the last DFT exception after all strategies exhaust.
- [ ] Make trajectory resume authoritative: if a persistent DFT trajectory or checkpoint is recovered, skip MLP/xTB re-relaxation unless explicitly requested.
- [ ] Harden the high-level single-point stage against hard runtime exceptions, not only SCF non-convergence, so recovered geometries still yield an energy artifact when scientifically defensible.
- [ ] Add focused unit coverage for: non-TS best-geometry single-point fallback, authoritative resume semantics, and single-point exception fallback behavior.
- [ ] Revalidate in Docker on the computational-gap lane before any new long DFT restart is considered complete.

Implementation notes 2026-04-17:

- Change order must be conservative:
	1. encode the intended contracts in tests first,
	2. unify fallback semantics for TS and non-TS,
	3. make resume skip upstream relaxers only when the recovered geometry is DFT-derived,
	4. harden single-point error handling without silently changing scientific meaning,
	5. rerun the targeted Docker/unit validation lane,
	6. only then trust the next heavy execution result.
- Scientific guardrail: fallback energies are acceptable only when surfaced as best-effort artifacts; they must not be mislabeled as fully converged solvent-corrected thermochemistry.

## Active Follow-up — 2026-04-16 SCF Engine Stabilization

- [x] Diagnose the SOSCF / Fermi-smearing engine-construction failure in the DFT builder.
- [x] Enforce a hard compatibility rule that disables Newton / SOSCF whenever Fermi smearing is active.
- [x] Split the ladder into a Stage 3 warm-up with DIIS + Fermi smearing and a Stage 4 refinement with Newton and no smearing.
- [x] Add focused unit coverage for the builder compatibility rule, the Newton-only refinement path, and the staged ladder behavior.
- [x] **V6 architecture overhaul**: replace the monolithic warm-up + refinement split with SCF strategy escalation.

Review 2026-04-16 (V6 — final):

- Root cause (v1–v3): smearing + Newton/SOSCF are fundamentally incompatible in PySCF for meta-GGAs like r2SCAN.
- Root cause (v4–v5): even after separating them, the warm-up stage used extreme settings (level_shift=1.0, conv_tol=1e-6, 200 cycles) that still failed to converge — and `_optimize_with_geometric_kernel` crashed fatally on unconverged SCF instead of recovering partial geometry.
- Architecture fix (two layers):
  1. `_optimize_with_geometric_kernel` now catches SCF gradient RuntimeErrors and recovers the best geometry from the geomeTRIC trajectory instead of crashing.
  2. `optimize_geometry` replaced the fixed Stage 3/4 warm-up+refinement with a **strategy escalation loop**: standard-DIIS → hardened-DIIS+smearing → SOSCF/Newton → vacuum-SOSCF. Each failed strategy passes its best partial geometry to the next.
- Also fixed a `DFTResult` dataclass construction bug (`electronic_energy_hartree` → `energy_hartree`).
- Focused Docker validation: 20 tests passed across `test_dft_refiner_scf_builder.py`, `test_dft_refiner_fallback.py`, `test_dft_refiner_irc_api.py`, and `test_run_computational_gap_dft.py`.

## Active Follow-up — 2026-04-16 Barrier Readiness Validation & Low-Risk Runtime Cleanup

- [x] Refresh the live DFT preflight state in Docker instead of trusting stale execution JSON.
- [x] Re-run the targeted computational-gap Docker tests before pruning runtime helpers.
- [x] Confirm whether the barrier workflow is actually launch-ready at the preflight layer.
- [x] Remove only low-risk duplicated XYZ helper logic in the runtime and support modules.
- [x] Re-run the targeted Docker validation lane after the helper cleanup.

Review 2026-04-16:

- Docker preflight now refreshes all seven computational-gap targets to `ready_for_dft`, replacing stale `running` statuses in per-target execution JSON artifacts.
- The targeted computational-gap validation lane is green in Docker: 19 tests passed across the DFT runner, xTB runner, forming-bond inference, xTB quality gates, proxy readiness, and the scientific refinement-plan artifact.
- The runtime remains only preflight-ready, not universally execution-proven: `lysinoalanine_crosslink` still has a recorded full DFT gradient-convergence failure, so the cleanup pass must stay behavior-preserving.
- The first cleanup pass is limited to low-risk deduplication of XYZ helper logic and similar non-scientific surface area.
- Post-cleanup validation stayed green: the same 19 Docker tests still pass, and representative preflight reruns for `lysinoalanine_crosslink`, `pe_schiff_base`, and `quinone_cys_michael` still return `ready_for_dft`.

## Active Follow-up — 2026-04-16 Technical-Debt Reduction Pass

- [x] Audit the newest standalone scripts and probes for real references before keeping them.
- [x] Delete clearly orphaned code and redundant wrappers when the core runtime already exposes the same behavior.
- [x] Re-run targeted reference checks to make sure the cleanup did not break active imports or documented lanes.
- [x] Record which recent additions remain intentionally active so the next cleanup pass does not revisit settled files.

Review 2026-04-16:

- The active `src/` hardening modules remain justified: they are still wired into runtime or tests.
- Removed the root-level callback probe and the manual `fix_lysinoalanine_docking.py` helper because neither had live references.
- Kept `generate_explicit_water_seed.py` because the Family 13 explicit-water lane is still an open roadmap item, and kept `run_dft_queue.sh` because it is an operational wrapper rather than a proven orphan.
- Targeted post-cleanup checks found no remaining references to the deleted files, and the touched explicit-water generator is error-free.

## Active Follow-up — 2026-04-16 Import Recovery And MLP Fallback Audit

- [x] Restore compatibility imports required by the remaining collection failures.
- [x] Repair the partially applied `dft_refiner.py` / `mlp_optimizer.py` backend-fallback changes and revalidate static analysis.
- [ ] Re-run the previously failing unit/scientific tests and confirm collection is clean.
- [ ] Decide the safe runtime posture for non-MACE MLP fallback candidates when MACE and xTB fail.
- [ ] Remove or keep explicitly dead code paths only after test confirmation.

Review 2026-04-16:
- Restored the missing compatibility surfaces for matrix benchmark materialization, DFT coverage metadata, projection canonicalization, and MLP backend helper APIs.
- Repaired the broken intermediate state in `src/dft_refiner.py` so the pre-relax/TS-recovery path is coherent again and a structured `compute_irc(...)` API is exported for tests.
- Static validation is clean on the touched files; test reruns and the final code-pruning decision are the remaining gates.

## Active Follow-up — 2026-04-13 Steric Clash Fix & DFT Queue Update

- [x] Diagnose `lysinoalanine_crosslink` and `pe_schiff_base` DFT preflight steric clash blockers.
- [x] Root-cause identified: multi-fragment RDKit embeddings produced overlapping fragments (0.47 Å H–C clash in lysinoalanine, 0.30 Å O–H clash in pe_schiff_base). The xTB path search returned identity paths (reactant ≈ TS), making the previously reported ~3 kcal/mol barrier for lysinoalanine invalid.
- [x] Sanity-checked the three passing targets (`hexanal_radical_quench`, `quinone_cys_michael`, `pe_amadori`). Result: `quinone_cys_michael` also has an identity path (pairwise RMS delta = 0.0000 Å) despite passing the 0.6 Å gate — its xTB path is fake too.
- [x] Regenerated reactant/product geometries for all three broken targets with proper fragment separation (`scripts/fix_steric_clashes.py`). All three now pass the 0.6 Å minimum distance gate with balanced atom counts preserved.
- [ ] Re-run xTB path search for `lysinoalanine_crosslink`, `pe_schiff_base`, and `quinone_cys_michael` after the geometry fix (requires Docker):
	```
	./scripts/docker_maillard.sh computational-gap-xtb lysinoalanine_crosslink
	./scripts/docker_maillard.sh computational-gap-xtb pe_schiff_base
	./scripts/docker_maillard.sh computational-gap-xtb quinone_cys_michael
	```
- [ ] After xTB completes, verify the new TS guesses pass the DFT preflight gate (pairwise RMS delta >> 0.01 Å, confirming a real reaction path).
- [ ] Re-run DFT preflight to confirm DFT-readiness:
	```
	./scripts/docker_maillard.sh computational-gap-dft-preflight lysinoalanine_crosslink
	./scripts/docker_maillard.sh computational-gap-dft-preflight pe_schiff_base
	./scripts/docker_maillard.sh computational-gap-dft-preflight quinone_cys_michael
	```

Review 2026-04-13 (steric clash):
- The stale xTB path outputs for all three targets have been removed. The reactant and product XYZ files have been regenerated from corrected SMILES with proper multi-fragment spacing.
- `aa_ring_open_dicarbonyl` remains actively running in DFT (TS geometry optimization phase, 4+ hours elapsed).
- `hexanal_radical_quench` and `pe_amadori` remain the only two targets with both clean preflight AND real xTB paths. They should be the first post-Family-14 DFT launches.
- The updated DFT queue priority is: (1) wait for `aa_ring_open_dicarbonyl` to complete, (2) launch `hexanal_radical_quench`, (3) launch `pe_amadori`, (4) after xTB reruns: `lysinoalanine_crosslink`, `pe_schiff_base`, `quinone_cys_michael`.


## Active Follow-up — 2026-04-13 QM Queue Realignment

- [x] Confirm whether the rerun of `pe_schiff_base` produces a usable xTB TS guess or remains proxy-only.
- [x] Confirm whether `aa_ring_open_dicarbonyl` is still actively running in DFT before launching any additional heavy job.
- [x] Decide whether Family 15 PE lanes should enter the formal DFT manifest or stay as proxy lanes.
- [x] Align README and the computational-gap runbook with the real selective-QM queue and with the semantics of the embedded validation plot.

Review 2026-04-13:
- `aa_ring_open_dicarbonyl` remains actively running in DFT, so `hexanal_radical_quench` should not be launched yet.
- Managed Docker xTB artifacts now exist for both `pe_schiff_base` and `pe_amadori`, and the Family 15 pair gate passes cleanly in `results/validation/computational_gap_proxy_readiness.{md,json}`.
- The formal queue is therefore `hexanal_radical_quench`, `lysinoalanine_crosslink`, `aa_ring_open_dicarbonyl`, `quinone_cys_michael`, `pe_schiff_base`, and `pe_amadori`, while `asparagine_sugar_explicit_water_cluster` and `melanoidin_radical_trapping` remain outside the copy-paste DFT-ready set.
- The README plot `docs/assets/validation_overview.png` is intentionally a benchmark-parity plot, not a direct visualization of every computational-gap target, so the README and runbook were updated to state that explicitly.

## Active Follow-up — 2026-04-13 Scientific Delivery Sequence

- [ ] Finish the restarted controlled `aa_ring_open_dicarbonyl` DFT run, ingest the result, and decide whether the Family 14 surrogate can be narrowed or promoted.
- [ ] Family 14 completion contract:
	- preflight rule: require `results/computational_gap_refinement/aa_ring_open_dicarbonyl_dft_execution.json` to report a clean geometry gate before any heavy DFT restart
	- current restart state: the legacy uncontrolled run was stopped and the target was relaunched through the phase-tracked runner using repaired seed `data/geometries/xtb_inputs/aa_ring_open_dicarbonyl/xtbpath_0.xyz`
	- execution artifact to inspect: `results/computational_gap_refinement/aa_ring_open_dicarbonyl_dft_execution.json`
	- post-run actions: ingest DFT result, compare against the HCW surrogate, and decide whether the live prior stays bounded-calibration or gets narrowed
	- acceptance rule: no write-back if the result only reproduces the current surrogate uncertainty without tightening it
- [ ] Launch `hexanal_radical_quench` DFT immediately after Family 14 finishes; keep it ranking-only unless the ingestion artifacts remain benchmark-honest.
- [ ] Family 11 execution contract:
	- only start after Family 14 frees the machine
	- ingest as ranking-only support even if the barrier is chemically plausible
	- require the write-back artifact to keep observable-first governance explicit
- [ ] Re-run or stabilize `lysinoalanine_crosslink` only after the geometry-refresh pass so we do not waste DFT cycles on fragile surrogate coordinates.
- [ ] Family 12 stabilization contract:
	- complete the public-geometry refresh first
	- compare SCF stability against the current surrogate seeds before spending new DFT cycles
	- only then rerun xTB and DFT
- [x] Decide the Family 15 promotion rule from evidence, not convenience: the managed pair gate now passes, so `pe_schiff_base` and `pe_amadori` enter the official DFT-ready manifest as bounded-calibration candidates.
- [ ] Family 15 promotion contract:
	- artifact of record: `results/validation/computational_gap_proxy_readiness.{md,json}`
	- pair gate: both `pe_schiff_base` and `pe_amadori` need `xtbpath.xyz`, `xtbpath_ts.xyz`, and a clean xTB quality gate
	- managed rule: `completed` or `completed_cached` execution artifacts are acceptable only when the managed quality gate also passes
	- proxy rule: synthesized outputs with SCF warnings remain proxy-only and must not be treated as formal DFT-ready evidence
- [ ] Generate the explicit-water seed for `asparagine_sugar_explicit_water_cluster` and keep its ceiling at uncertainty narrowing only.
- [ ] Family 13 safety-lane contract:
	- generate the explicit-water seed as a separate preparatory job
	- do not merge it into the standard xTB path-search queue
	- keep all downstream claims capped at uncertainty narrowing only
- [ ] Keep Family 16 (`melanoidin_radical_trapping`) MLP-first unless a smaller fragment proxy yields a cleaner xTB path than the current polymer proxy.
- [ ] Family 16 triage contract:
	- shortlist smaller fragment proxies first
	- only escalate to xTB if a fragment produces a cleaner path than the current polymer proxy
	- keep the lane out of the formal DFT queue until that fragment screen exists
- [ ] After the QM queue stabilizes, write back only the deltas that genuinely tighten priors or uncertainty bands in the live runtime.

Review 2026-04-13:
- The scientific order is now explicit: Family 14 completion, then Family 11, then Family 13/15 expansion, then a refreshed Family 12 attempt, while Family 16 stays outside the formal queue until a smaller-fragment screen exists.
- Family 14 is no longer running as an opaque legacy job. The old host-side DFT process was stopped after the new geometry gate flagged the original TS guess as unsafe, the runner auto-repaired the seed to `xtbpath_0.xyz`, and the target was relaunched through the controlled Docker flow with phase tracking.
- Family 13 to 16 are already present in the runtime surface; the remaining scientific work is promotion quality and uncertainty narrowing, not basic model inclusion.

## Active Follow-up — 2026-04-13 Architecture Hardening

- [ ] Remove hidden environment assumptions from scientific generators: LaTeX-backed figure generation should be configured once and fail explicitly if the required toolchain is absent.
- [ ] Reduce duplicated plot-style bootstrapping across validation generators by using one shared helper and one environment contract.
- [ ] Split computational-gap execution state from scientific promotion state more cleanly: a synthesized xTB path with SCF warnings must not look equivalent to a clean path in downstream tooling.
- [ ] Promote proxy-readiness diagnostics into machine-readable artifacts so Family 15 and future lanes are not tracked only in chat or ad hoc shell output.
- [ ] Continue shrinking duplicated script-local helper logic when a shared `src/` helper makes the contract clearer without expanding the public surface.
- [ ] Keep public docs aligned with machine-readable artifacts; README should distinguish benchmark parity, runtime family coverage, and selective-QM queue state.

Detailed implementation sequence:

- [ ] Plotting contract hardening:
	- shared helper in `src/plot_style.py`
	- explicit failure if the LaTeX toolchain is incomplete
	- no script-local style bootstrapping or silent fallback
- [ ] xTB quality-state hardening:
	- shared helper for path materialization and log-health inspection
	- explicit `completed_with_quality_warnings` state in xTB execution outputs
	- downstream promotion logic must read quality, not just file presence
- [ ] Proxy artifact generation:
	- write `computational_gap_proxy_readiness.{md,json}` from source
	- encode pair-level promotion rules for Family 15
	- keep public docs and runbooks tied to that artifact instead of shell transcripts

Review 2026-04-13:
- The next architecture wins are not broad rewrites. They are contract hardening and duplication removal around plotting, computational-gap status encoding, and public artifact generation.
- The maintainability bar for new scientific scripts is now stricter: less local setup, fewer environment assumptions, and clearer machine-readable status transitions.

## Active Follow-up — 2026-04-09 Maintainability Audit II

- [x] Inventory structural duplication, stale wrappers, and scripts that still expose machine-specific or ad hoc execution patterns.
- [x] Identify outdated code paths that are now bypassed by the current runtime, Docker entrypoints, or generator layout.
- [x] Propose the smallest safe structural reductions that improve discoverability without breaking active scientific workflows.
- [x] Validate any recommended removals or rewrites against the current test and execution surface before landing changes.

Review 2026-04-09:
- In progress. The first cleanup wave removed obvious dead files. This pass is focused on higher-value maintainability problems: duplicated entrypoints, hidden execution paths, legacy helper patterns, and modules whose structure no longer matches the current runtime surface.
- Confirmed debt clusters after code archaeology: duplicated scientist-facing CLI/report wiring across `run_pipeline.py`, `optimize_formulation.py`, `compare_formulations.py`, and `run_campaign.py`; package-internal `sys.path` mutation inside `src/` and even `data/reactions/curated_pathways.py`; documentation that still promotes `python scripts/...` while the validated Docker wrapper exposes only part of that surface; and cross-generator imports such as `generate_validation_figures.py` importing `_build_rows` from `generate_benchmark_coverage_gaps.py`.
- The next safe refactor target is not broad deletion. It is extracting a small shared scientist-facing execution helper, moving curated pathway Python definitions out of `data/`, and breaking generator-to-generator library imports by promoting shared builders into `src/`.
- Landed 2026-04-10: `src/scientist_cli.py` now centralizes CLI-side confidence assessment, presentation rendering, and optional report generation so the four scientist-facing entrypoints no longer duplicate that wiring. `src/benchmark_coverage_gaps.py` now owns the reusable coverage-gap builder/markdown logic, while `generate_validation_figures.py` no longer imports private helpers from another generator script.
- Follow-up 2026-04-10: that first helper extraction was then compacted to reduce total code size. The CLI confidence helper was absorbed into `src/usability_reports.py` and the extra `src/scientist_cli.py` module was removed. The current maintainability batch is net negative in size (`130` insertions, `230` deletions across the touched files).
- Also landed 2026-04-10: removed unnecessary `sys.path` mutation from `src/recommend.py`, `src/smirks_engine.py`, and `data/reactions/curated_pathways.py`. Those modules now rely on the normal repo import surface instead of self-mutating package state.
- Landed 2026-04-10: the hand-curated pathway module was moved from `data/reactions/curated_pathways.py` into `src/curated_pathways.py`, its ad hoc `__main__` block and decorative spacing were removed, and consumers were rewired to import from `src`. That step alone is net negative in size (`3` insertions, `266` deletions across the touched files) while removing Python code from the data tree.
- Landed 2026-04-10: `scripts/docker_maillard.sh` was compacted without changing its public command surface. Generator dispatch now uses shared helpers and the `validation-figures` bundle is centralized in one lane. The wrapper refactor is net negative in size (`99` insertions, `109` deletions) and was validated with `bash -n`, `./scripts/docker_maillard.sh summary`, `./scripts/docker_maillard.sh coverage-gaps`, and `./scripts/docker_maillard.sh validation-figures`.
- Focused Docker validation passed after the refactor: `tests/unit/test_cli_scripts.py`, `tests/unit/test_benchmark_coverage_gaps.py`, `tests/scientific/test_validation_figures.py`, `tests/unit/test_data_integrity.py`, `tests/unit/test_target_observability_metadata.py`, and `tests/unit/test_smirks_engine.py` (`68 passed`).

## Active Follow-up — 2026-04-09 Public 3D Geometry Refresh

- [ ] Download public 3D templates for the three active computational-gap targets and record provenance.
- [ ] Preserve the mapped atom ordering required by xTB path search while replacing local surrogate coordinates with more realistic 3D seeds.
- [ ] Re-run xTB for `hexanal_radical_quench`, `lysinoalanine_crosslink`, and `aa_ring_open_dicarbonyl` after refreshing the seeds.
- [ ] Re-run per-target DFT jobs and capture whether the refreshed seeds reduce immediate SCF failures.

Review 2026-04-09:
- In progress. The current local XYZ files are balanced and xTB-ready, but they were built from simplified mapped SMILES. The refresh pass is using public PubChem 3D records where available, while keeping the local atom ordering invariant so `xtb --path` still interpolates the intended atoms.

## Active Follow-up — 2026-04-09 Git Tracking Audit

- [x] Inventory the files currently tracked even though they now match ignore rules.
- [x] Quantify the largest tracked files and the heaviest historical blobs.
- [x] Classify candidates into keep in git, move to releases/LFS, or stop tracking entirely.
- [x] Document the exact untracking command set and any history-rewrite follow-up if the repo should be slimmed down on GitHub.

Review 2026-04-09:
- The immediate untracking surface is dominated by generated xTB path-search XYZ files under `data/geometries/xtb_inputs/`, solver CSV/PNG outputs under the RMG validation case, one notebook SQLite artifact, and the currently tracked `results/` outputs already covered by `.gitignore`.
- The largest live tracked file is `data/lit/foods-15-00912.pdf` at about `2.0 MB`; the largest historical blob is a removed `Miniforge3-Darwin-arm64.sh` installer at about `65.6 MB`, which means the GitHub size concern is partly historical, not just the current tree.
- Current tracked content is roughly `12.5 MB` under `data/`, `1.47 MB` under `docs/`, `1.45 MB` under `src/`, and `0.49 MB` under `results/`, so the highest-return cleanup target is generated scientific artifacts rather than source code.
- Direct `git rm --cached` candidates already matched by ignore rules: `102` xTB path-search XYZ files (`~7.44 MB`), `24` solver CSV/PNG outputs (`~1.32 MB`), `7` `results/` artifacts (`~0.49 MB`), and `docs/notebooks/results/maillard_results.db` (`~28.7 KB`).
- Recommended policy: keep source code, curated machine-readable inputs, and hand-maintained docs in git; stop tracking generated scientific artifacts and local notebook databases; consider moving literature PDFs to external storage or LFS only if offline archival is a real product requirement.
- Minimal cleanup command for the current tree: `git ls-files -ci --exclude-standard -z | xargs -0 git rm --cached --`.
- Executed cleanup on 2026-04-09: the working tree now has `135` staged untracking changes (`102` xTB XYZ artifacts, `24` solver outputs, `7` `results/` artifacts, `1` notebook DB, `1` literature PDF), and `git ls-files -ci --exclude-standard` now returns `0`.
- If GitHub repo size still matters after that, the historical cleanup target is the removed `Miniforge3-Darwin-arm64.sh` blob (`b538bbc0986fd68c0c9e483670b3a529d96c7383`), which should be purged with `git filter-repo` or BFG rather than by changing the working tree alone.

## Mission

Build the most useful computational tool for scientists who want to imitate meat-like Maillard chemistry in plant-based systems.

## Product Thesis

This is not mainly a problem of making one chemistry engine more exact.

It is a problem of combining:

- a quantitatively credible free-precursor core
- matrix-aware observability and accessibility
- process-aware confidence boundaries
- scientist-facing reporting that states what is benchmarked, what is transferred, and what remains blocked by missing external evidence

The product question remains:

> Which formulation and process changes are most worth testing next if the goal is meat-like aroma under plant-matrix constraints?

## Active Follow-up — 2026-04-09 Surrogate Closure Pass

- [x] Audit the blocked computational-gap targets against the live prior registry and runtime consumers.
- [x] Replace the invalid lysinoalanine xTB-derived posture with an explicit DHA-plus-lysine mechanistic surrogate anchored in the repo's DHA_Crosslinking family.
- [x] Make the ascorbic-acid ring-opening prior explicitly describe its bounded Family 14 surrogate posture rather than implying a blocked xTB path is the active authority.
- [x] Regenerate the computational-gap refinement plan/manifests and add focused tests that lock the new surrogate metadata in place.

Review 2026-04-09:
- `lysinoalanine_crosslink` no longer pretends to be backed by a valid xTB elementary path. Its live prior is now an explicit DHA-plus-lysine family surrogate anchored to the repo's `DHA_Crosslinking` rule (`E0 = 16 kcal/mol`) and kept at `ranking_only` posture.
- `aa_ring_open_dicarbonyl` remains a bounded-calibration prior, but the plan and prior registry now say explicitly that the active authority is the Family 14 HCW surrogate, not the blocked placeholder xTB pair.
- The refinement plan and xTB/DFT manifests now surface surrogate metadata directly for blocked targets so the user sees both truths at once: the current geometries are unusable, and a bounded mechanistic fallback exists.
- Focused validation passed in the `maillard` conda environment: `tests/scientific/test_computational_gap_refinement_plan.py` and the Family 14 runtime check in `tests/unit/test_literature_runtime.py`.

## Active Follow-up — 2026-04-09 Balanced Geometry Recovery

- [x] Replace the invalid placeholder XYZ pairs for `lysinoalanine_crosslink` and `aa_ring_open_dicarbonyl` with balanced mapped surrogate reactions.
- [x] Regenerate xTB input geometries for both targets and verify reactant/product atom counts match.
- [x] Execute the xTB wrapper for both targets and materialize `xtbpath.xyz` plus `xtbpath_ts.xyz` so the DFT manifest can advance.
- [x] Regenerate the computational-gap refinement plan/manifests and align focused tests with the new ready/ready_for_dft state.

Review 2026-04-09:
- `lysinoalanine_crosslink` now uses a balanced dehydroalanine-plus-lysine mapped pair (`35/35` atoms) instead of the previous invalid placeholder.
- `aa_ring_open_dicarbonyl` now uses a balanced dehydroascorbic-acid-plus-water to 2,3-diketogulonic-acid mapped pair (`21/21` atoms).
- The xTB execution wrapper completed successfully for both regenerated targets and materialized `xtbpath.xyz` and `xtbpath_ts.xyz`, which moves both lanes from xTB-blocked to `ready_for_dft` in the live plan surface.
- The published computational-gap plan/manifests were regenerated from source and focused validation now passes again in the `maillard` conda environment, including the stale Family 11 runtime expectation that had drifted behind the promoted hexanal prior.

---

## P1 — Immediate Execution Slice (2026-04-05)

### S15. Literature Extraction Sprint — Top-20 Backlog Papers

**Rationale:** 164 of 172 scored citations (95%) remain in BACKLOG. Among those are ~30 papers scored 8/8 containing directly encodable constants (Ea values, partition coefficients, stoichiometries). The current runtime surface uses only 21 + 8 staged references. A focused extraction sprint on the highest-scored papers would approximately double the parametric surface with no wet-lab work required.

**Status 2026-04-08:** the compact runtime-first backlog is effectively exhausted (`ready_runtime=0`, `ready_benchmark=0`). Remaining S15 work is now structural: only keep extracting literature when it unlocks a blocked benchmark surface, a cross-family transfer, or a still-missing executable benchmark.

**Highest Priority Now (2026-04-08):**

- S20. No-Wet-Lab Computational Closure: formal UQ (S20.1) and leave-one-out cross-validation (S20.2) are the highest-leverage items because they transform every existing prediction from a point estimate into an honest uncertainty range, and validate transfer assumptions without needing new data.
- S13. Scientist-Facing Visual Output: model credibility has overtaken report usability, so the next product bottleneck is turning validated predictions into charts and confidence views a food scientist can act on immediately.
- S21/S22. Data Ingestion + Experiment Recommender: close the scientist-model feedback loop so that when data does arrive, it can be ingested frictionlessly and the most valuable next experiment is always visible.
- S17. Extrusion Benchmark Validation: the extrusion lane is architecturally useful but still lacks an external-decision-ready benchmark surface.
- S15 structural unlock triage: the compact queue is exhausted, so further literature work should only continue when it opens a blocked benchmark, matrix lane, or cross-family transfer.

- [x] Re-rank the remaining backlog by structural unlock potential, not just raw score, to isolate the few citations that can still move calibration or executable benchmark coverage.
- [ ] Extract only the constants that open a blocked runtime or benchmark surface; do not keep broadening compact transfer payloads with redundant priors.
- [ ] Prefer landings that create a new executable benchmark or tighten a cross-family structural transfer instead of another single-benchmark local correction.
- [ ] Update the Deep Research backlog state only when the landing is wired into runtime, benchmark, or governance surfaces that survive Docker validation.
- [ ] Validate that no existing benchmark regresses after each structural batch of new payloads or contract changes.

Review 2026-04-08:
`results/validation/structural_unlock_triage.{md,json}` now makes the post-backlog posture explicit: the curated citation queue is exhausted (`ready_runtime=0`, `ready_benchmark=0`), so the next leverage is not another compact extraction batch. The ranked structural sequence is now primary PPI/SPI benchmark package first, extrusion benchmark translation second, and retention follow-on third. Deferred modeling items `5.7`, `5.10`, and `5.11` stay deferred; `5.8` is the only adjacent follow-on worth revisiting immediately after the first extrusion benchmark lands.

### S15.1. Active Entry Criteria After Runtime-First Exhaustion

- [ ] Open a new S15 landing only when it unlocks a blocked benchmark, executable workbook, or cross-family transfer named in `results/validation/structural_unlock_triage.{md,json}` or the no-wet-lab computational gap plan.
- [ ] Land the citation as a compact runtime, benchmark, or governance surface that survives focused Docker validation; do not keep expanding narrative-only payloads.
- [ ] Stop a candidate immediately if it does not change benchmark coverage, benchmark accuracy, or an executable scientist workflow.

Completed 2026-04-05 to 2026-04-08:

- Runtime-first queue publishing, auto-advance, and registry landings are complete.
- The Cerny 2008 failure repair and the SLR-identified runtime corrections are complete.
- The structural continuation slices through Families 02-16 are complete, including the Family 03 dilute-loading closure and the final ingredient-source EUC slice.
- The curated literature backlog now sits at `ready_runtime=0` and `ready_benchmark=0`, so S15 is no longer a broad sprint lane.

Verification 2026-04-08:
Focused Docker subsets covering runtime landing, backlog governance, sulfur/thiamine closure, payload wiring, and benchmark guards passed across the final continuation batches. The roadmap no longer needs the historical per-slice checklist here; the active rule is exception-only literature landing.

---

### S16. README Rewrite for Food Scientists

**Rationale:** The current README (359 lines) is written for computational chemists. The target user — a food scientist in an alternative protein company — needs to see what goes in, what comes out, and how much to trust it, within 2 minutes.

- [x] Rewrite README as a ~100-line landing page: "What this does → Install → Run → Example output → Where to learn more."
- [x] Add a sample output snapshot (radar chart, ranked formulations table, confidence annotations) so scientists see value before installing.
- [x] Replace internal vocabulary ("bench-neighborhood," "family-lane transparency," "validated envelope") with plain-language equivalents.
- [x] Simplify the architecture diagram to a 3-box version ("Ingredients + Process → Maillard Engine → Predictions + Confidence") as the first visual. Move the detailed Mermaid diagram to `docs/architecture.md`.
- [x] Consolidate the 16 validation artifact links into a single "click here to see what's validated today" reference.
- [x] Move detailed content to dedicated documents: `docs/VALIDATION.md`, `docs/PHILOSOPHY.md`.
- [x] Add explicit positioning statement against existing tools (NIST WebBook, RMG, manual literature review) so scientists understand the unique value.

---

### S13. Scientist-Facing Visual Output

**Rationale:** The pipeline currently outputs raw JSON which masks the deterministic and capability-bounded nature of the framework. A food scientist doesn't need generic sensory radar charts (which hide physics behind arbitrary scores). They need a Decision-Risk Dashboard that shows *why* a suggestion is made and *how much they should trust it*. This is the single biggest reporting gap.

- [ ] Add a **Relative Intervention Waterfall Chart** to `--report` output showing what changed: e.g., "Adding glutathione raised thiols by +X% but reduced total aldehydes by -Y%."
- [ ] Add **Confidence Bounds Overlay** plotting predictions with an error bar representing the Evidence Lane (e.g., narrow bounds for strict benchmarks vs. wide bounds for directional extrapolations).
- [ ] Add a **Capabilities Heatmap** showing Matrix Types (columns) against Process States (rows) to visually translate the 150+ reference papers into proof of capability.
- [ ] Add a parity plot export to the validation-figures command that a scientist can embed in presentations.
- [ ] Refactor the Pareto frontier visualization for multi-objective optimization (meaty-positive vs. safety vs. off-note suppression) to overlay extrapolation risk metrics.

---

## P2 — High-Impact Medium-Effort

### S17. Extrusion Benchmark Validation

**Rationale:** The extrusion lane now has executable protocol, landing, follow-on, and diagnostic surfaces. The remaining gap is no longer architecture; it is the absence of real same-run SPI or PPI extrusion measurements. Without wet-lab access, the correct goal is to narrow mechanism and process uncertainty computationally while keeping benchmark-closure claims withheld.

#### S17.1. Extrusion Benchmark Landing Infrastructure

- [x] Specify the minimum viable extrusion benchmark: one protein type (PPI or SPI), two SME levels, one barrel temperature, measuring MFT + hexanal + furosine simultaneously.
- [x] Publish the protocol, SOP lock register, closure package, and intake workbooks as shareable artifacts in `results/validation/`.
- [x] Expose executable workbook processors and diagnostic example bundles in Docker for both external closure and item `5.8` follow-on analysis.
- [ ] Land the first non-placeholder SPI or PPI extrusion workbook with real measurements; until then keep the lane explicitly non-closure.
- [ ] Leave unresolved lab-SOP fields explicit until an external SOP fixes the exact extruder brand/model and isotope spike concentrations.

Review 2026-04-08:
`results/validation/extrusion_benchmark_protocol.{md,json}` now publishes a review-ready S17.1 artifact generated from `src/doe_generator.py`. The current minimum viable design selects `soy_iso` with a fixed `145 C` barrel temperature and `120/180 kJ/kg` SME arms, measuring `MFT`, `hexanal`, and `furosine` with same-run extensions for `FFT`, `2-pentylfuran`, `2,5-dimethylpyrazine`, and `acrylamide`. The protocol is now tightened with the closest structured HME translation anchor already present in the repo: `57%` moisture, `280 rpm`, `4.6 kg/h`, and a `30/90/120/140/150/160 C` barrel profile plus control off-note anchors from `Li et al. (2026)`. The remaining blocker is now narrower: the repo still does not canonize a specific extruder model or exact spike concentrations for `[2H2]-MFT`, `[2H2]-FFT`, `hexanal-d12`, and `13C3-acrylamide`.

Review 2026-04-08, primary matrix package:
`results/validation/matrix_primary_benchmark_campaign.{md,json}` and `results/validation/primary_matrix_external_package.{md,json}` are now regenerated from live source in `src/primary_benchmark_campaign.py` instead of surviving as orphaned output-only artifacts. The repo now explicitly materializes the dual-arm pea/soy benchmark campaign, the pea-first external package, the current adverse-only external anchor, and the exact promotion delta that the wet-lab package is expected to close.

Review 2026-04-08, SOP and landing execution:
`results/validation/extrusion_sop_lock_register.{md,json}` now separates what the repo truly locks from what remains unresolved in S17.1: the process envelope is canonized as a twin-screw HME translation with fixed moisture, rpm, feed rate, barrel profile, die-exit temperature, and analytical platform identities, while exact brand/model and exact isotope spike concentrations remain explicitly unresolved. `results/validation/extrusion_external_closure_package.{md,json}` and `results/validation/extrusion_external_closure_workbook.{md,yaml}` now convert the existing external-closure contract into an exact per-arm landing bundle for the first SPI extrusion panel, including direct damage markers, same-run sulfur/off-note tradeoff targets, and required extrusion metadata. `src/extrusion_benchmark_execution.py` plus `scripts/generators/process_extrusion_external_closure_workbook.py` now turn a filled workbook into per-arm intake payloads and support-delta artifacts, and `./scripts/docker_maillard.sh extrusion-closure-workbook ...` exposes that path directly inside Docker. `results/validation/primary_matrix_external_package_intake_template.{md,yaml}` now turns the pea-first package into a shareable landing template with placeholder measurement fields, so the next wet-lab package can be ingested without inventing values in-source.

Review 2026-04-08, 5.8 landing package:
`results/validation/extrusion_disulfide_follow_on_package.{md,json}` and `results/validation/extrusion_disulfide_follow_on_workbook.{md,yaml}` now pre-wire item `5.8` to the first measured SPI extrusion panel instead of leaving it as a narrative next step. The follow-on package is explicitly keyed to the runtime's existing disulfide context and retention priors (`raman_sds_extrusion_disulfide_severity`, `acs_jafc_3c02618_mft_disulfide_trapping_v1`, `acs_jafc_0c01925_protein_binding_hierarchy_v1`) and defines the exact same-run observables and derived ratios needed before any benchmark-fitted sulfur-retention coefficient is claimed. `src/extrusion_benchmark_execution.py` plus `scripts/generators/process_extrusion_disulfide_follow_on_workbook.py` now evaluate a filled workbook against the measured soy protocol-pilot reference, compute the 5.8 derived metrics, and expose the lane via `./scripts/docker_maillard.sh extrusion-follow-on-workbook ...`.

Review 2026-04-08, diagnostic rehearsal:
`results/validation/extrusion_diagnostic_examples_bundle.json` plus the paired diagnostic workbooks and execution reports now exercise both extrusion workbook lanes end-to-end without inventing wet-lab truth. This is the right stopping point until real measurements exist.

#### S17.2. Extrusion Model Extensions

- [x] Add volatile stripping correction at the die (flash-vaporization loss based on die temperature and compound vapor pressure).
- [x] Add shear-volatile coupling beyond the simple linear SME→ΔT slope (cell-wall rupture → precursor release, protein aggregation → trapping landscape).
- [x] Evaluate whether a simple RTD (residence time distribution) model is needed or if the sequential-zone model is sufficient for the target use case.

Review 2026-04-08, S17.2:
`src/extrusion.py` now extends the sequential-zone extrusion profile with die-exit temperature, mean residence-time estimate, RTD spread, explicit RTD sufficiency evaluation, and compound-class transport modifiers that combine shear-driven precursor release with aggregation/trapping and die-loss pressure. `src/headspace.py`, `src/sensory.py`, and `src/pipeline.py` now consume that extrusion transport layer so extrusion-aware runs do not just report higher effective temperature; they also propagate bounded die stripping and class-specific release penalties into headspace-facing sensory predictions. For the current scientist-facing HME use case, the new RTD assessment keeps the sequential-zone model as sufficient until screw-geometry-specific benchmark data exists.

#### S17.3. Computational Bridge Without Wet Lab

- [ ] Keep extrusion status at `computationally narrowed` rather than `benchmark closed` until same-run extrusion observables exist.
- [ ] Use xTB to generate motif-level path-search seeds for sulfur retention, hexanal suppression, lysinoalanine formation, and acrylamide explicit-water clusters that could modify extrusion-facing priors.
- [ ] Use MLP only as an offline accelerator for geometry preoptimization where benchmarked (`mace_mp_small` today); do not use MLP barrier predictions as authority.
- [ ] Run selective r2SCAN-3c plus the existing single-point refinement chain only for motifs that survive xTB screening and map cleanly to a runtime patch or uncertainty-band tightening.
- [ ] Write back only cached correction deltas, uncertainty narrowing, and provenance artifacts; do not invent synthetic extrusion benchmark measurements.

Review 2026-04-08:
`results/validation/computational_gap_plan_no_wet_lab.md` now captures the no-wet-lab closure policy, the remaining irreducible wet-lab gaps, and the exact xTB, MLP, and DFT roles that are still honest for extrusion-adjacent uncertainty reduction.

### S18. Selective xTB/DFT/MLP Computational Closure

**Rationale:** Current refinement governance still says `hold_observable_first`, so selective QM should run only as a bounded prior-refinement lane, not as a substitute for missing matrix benchmarks. The useful no-wet-lab targets remain narrow and motif-specific.

- [ ] Family 11: run `hexanal_radical_quench` as xTB path search → optional MLP geometry preopt → r2SCAN-3c refinement, then ingest only as a ranking-only off-note suppression anchor.
- [ ] Family 12: run `lysinoalanine_crosslink` through xTB seed generation plus r2SCAN-3c refinement, then ingest only as a bounded AGE/ALE damage prior.
- [ ] Family 14: generate and refine `aa_ring_open_dicarbonyl`, then ingest only as a bounded stealth-browning prior.
- [ ] Safety lane: evaluate the asparagine-sugar transition state in a minimal explicit-water cluster and use the result only to tighten the acrylamide uncertainty band, not to claim matrix closure.
- [ ] Expand geometry benchmark coverage before adopting any MLP `ts_initialization` role beyond the current `mace_mp_small` geometry-preopt allowance.

Review 2026-04-08:
`results/validation/refinement_governance.md` still reports `hold_observable_first` with `approved_offline_jobs=0`, `results/validation/selective_dft_plan.md` still reports `run_now=0`, and `results/validation/mlp_adoption_notes.md` adopts only `mace_mp_small` for offline geometry preoptimization while quarantining `mace_off_medium` barrier predictions. The computational closure lane should therefore update priors and uncertainty bands only after backtesting, never benchmark status.

### S19. Web Interface (Minimal)

**Rationale:** Every interaction currently requires Docker + command-line. For food scientists, this is a major adoption barrier. A minimal web interface would increase adoption by an order of magnitude.

- [ ] Build a minimal Flask/FastAPI web interface with a formulation input form, "run prediction" button, and visual report output.
- [ ] Serve the radar chart, kinetic traces, and safety dashboard from S13 in the web response.
- [ ] Include a "download report" button for shareable PDF/HTML export.

### S20. No-Wet-Lab Computational Closure

**Rationale:** External-decision-ready matrix benchmarks cannot be closed computationally. But the framework can still make significant, honest progress by: (1) replacing heuristic trust tiers with formal uncertainty quantification, (2) validating transfer assumptions via cross-validation, (3) mining existing literature for external validation points, and (4) identifying which parameters actually drive formulation decisions vs. which are noise. The gap between "directionally useful" and "quantitatively calibrated with explicit uncertainty" is closable without a single new experiment.

**Strategic posture:** These items produce the maximum trust improvement per unit of effort without wet-lab access. They transform the tool from "the model says X" to "the model says X ± Y, and here's exactly what would need to be measured to shrink Y."

#### S20.1. Formal Uncertainty Quantification

- [ ] Build `src/uncertainty_propagation.py`: define prior distributions for each uncertain parameter (barrier heights: ±3 kcal/mol, matrix corrections: ±30%, retention factors: ±50%).
- [ ] Implement Monte Carlo propagation (N=500–1000) through `MaillardPipeline.evaluate()` to produce per-compound percentile predictions (5th/25th/50th/75th/95th).
- [ ] Add parametric sensitivity indices (Sobol or Morris method) to rank which parameters drive decision-relevant variance.
- [ ] Store as `UncertaintyEnvelope` dataclass in `FormulationResult` and wire into `--report` output as `predicted: 0.038 ppb [0.012—0.089, 90% CI]`.

#### S20.2. Leave-One-Out Cross-Validation

- [ ] Build `src/cross_validation.py`: for each of the 19 benchmarks, exclude it from calibration and re-evaluate all remaining benchmarks.
- [ ] Compare LOO residuals to full-model residuals per benchmark and per compound class.
- [ ] Identify "load-bearing" benchmarks (removing them degrades generalization significantly).
- [ ] Wire results into the validation surface so scientists see transfer robustness, not just accuracy.

#### S20.3. Observable Projection Re-calibration

- [ ] Systematically fit matrix correction offsets by minimizing prediction error across the 19 benchmarks per protein type.
- [ ] Compare fitted offsets against current heuristic values in `calibration_offsets.json`.
- [ ] Adopt improved offsets only if LOO (S20.2) confirms they generalize; otherwise document that current heuristics are already near-optimal.

#### S20.4. Extended Sensitivity Analysis

- [ ] Extend `family_sensitivity.py` barrier offset sweeps to also perturb matrix correction factors (per protein type), headspace partition coefficients (Henry constants), and retention/release factors.
- [ ] Rank all parameters by decision impact: "which parameter change would flip a formulation recommendation?"
- [ ] Output as a ranked importance chart for S13 visual dashboard.

#### S20.5. Literature Mining for External Validation

- [ ] Scan the 150-paper literature backlog for quantitative ppb measurements in plant protein matrices that can serve as hold-out external validation points.
- [ ] For each usable measurement: create an external validation payload via `matrix_experiment_intake.py` WITHOUT fitting to it.
- [ ] Report standalone external validation accuracy: "on N external measurements not used in calibration, the framework achieves X× median accuracy."

---

### S21. Data Ingestion Pipeline

**Rationale:** The repo already has strong intake infrastructure (`matrix_experiment_intake.py`, extrusion workbook processors, intake templates with YAML placeholders). But a scientist with GC-MS results doesn't know which module to use, can't ingest raw instrument files, and can't see how new data changes the validation surface before committing. The ingestion path should be as frictionless as the prediction path.

#### S21.1. Universal Intake CLI (`maillard ingest`)

- [ ] Build `src/data_ingest.py` with a CLI entry point: `maillard ingest --file results.csv --protein-type soy_iso --process-state extrusion_structured`.
- [ ] Auto-detect file format (CSV/Excel/YAML/JSON) and map column headers to the existing intake schema using fuzzy matching.
- [ ] Present the scientist with a confirmation summary before generating a normalized YAML payload via `matrix_experiment_intake.py`.
- [ ] Run `build_matrix_experiment_support_delta_artifact()` automatically to compute the ingestion impact.

#### S21.2. Impact Preview

- [ ] Before committing ingested data, display a clear summary: benchmarks added, validation surface change, promotion readiness change, compounds strengthened/weakened, median ratio shift.
- [ ] Require explicit confirmation before persisting.
- [ ] Log the full before/after state for audit.

#### S21.3. Benchmark Versioning

- [ ] Maintain `data/benchmarks/CHANGELOG.yaml` recording: when each benchmark was added, active calibration offsets at ingestion time, validation surface state before/after, and source provenance (DOI, lab, date).
- [ ] Enable reproducibility: "this prediction was made with benchmark set v12, which included the following 20 external validations."

---

### S22. Most-Valuable-Experiment Recommender

**Rationale:** The repo has static DoE templates (`doe_generator.py`) and specific benchmark campaign protocols (`primary_benchmark_campaign.py`). But there's no system that ranks experiments by expected calibration gain. A scientist with limited lab time needs the answer to: "Which single experiment will improve the model's predictions the most?" This is the feature that closes the loop between computational modeling and wet-lab prioritization.

**Dependency:** Requires S20.1 (UQ) and S20.4 (extended sensitivity) to be implemented first.

#### S22.1. Value-of-Information Engine

- [ ] Build `src/experiment_value.py`: for each potential experiment in the gap registry, estimate Expected Calibration Gain (ECG) = Σ [ current_uncertainty × decision_relevance × probability_of_improvement ] per compound.
- [ ] `current_uncertainty` from S20.1 UQ propagation.
- [ ] `decision_relevance` from S20.4 sensitivity analysis (does this compound flip formulation recommendations?).
- [ ] `probability_of_improvement` from how close the current evidence is to closing the gap.
- [ ] Output a ranked experiment list with ECG scores and estimated lab time.

#### S22.2. Experiment Protocol Generator

- [ ] Wrap existing `doe_generator.py` + `primary_benchmark_campaign.py` with ECG ranking, budget-aware experiment selection (combine multiple cheap experiments vs. one expensive one), and protocol export as PDF/Markdown.
- [ ] Auto-generate a pre-filled intake template (S21) for each recommended experiment so results can be directly ingested.
- [ ] CLI entry point: `maillard request-experiment --budget "3 lab days" --protein-type soy_iso --goal "meaty aroma optimization"`.

#### S22.3. Gap Visualization

- [ ] Add to the `--report` dashboard: a Gap Heatmap showing Matrix types × Process states × Compound classes, color-coded by uncertainty level from S20.1.
- [ ] Embed experiment suggestions in the report: "To improve this prediction from ±300% to ±50%, measure MFT headspace in SPI extrudate at 145°C."

---

### S14. Codebase Health & Maintainability

**Rationale:** `benchmark_validation.py` (2,611 lines) and `recommend.py` (1,338 lines) are monoliths that impede contribution and debugging. `literature_runtime.py` (4,324 lines) is the largest source file and the hidden monolith. Test suite runtime (~1h40m) blocks iteration speed.

- [ ] Decompose `literature_runtime.py` into modular components: loader, query, calibration, governance.
- [ ] Decompose `benchmark_validation.py` into modular components: registry, evaluation, reporting, assertion.
- [ ] Decompose `recommend.py` into modular components: concentration projection, observable mapping, scoring.
- [ ] Triage the test suite for performance: identify and optimize the slowest 10 tests, introduce pytest marks for fast/slow/full lanes.

---

## P3 — Strategic / Deferred

### S12. Scaling the Literature Pipeline & Uncertainty

#### S12.1: Formal Uncertainty Quantification (UQ)

- [ ] Replace narrative "trust heuristics" (e.g., Extrusion Exploratory Mode) with explicit mathematical confidence intervals (e.g., via parametric variance or Gaussian Processes) for out-of-domain predictions.
- [ ] Propagate UQ bounds into the predicted volatile headspace (ppb) figures so scientists know the exact variance of un-benchmarked estimates.

#### S12.2: Automated LLM-Assisted Payload Extraction

- [ ] Build an automated ingestion pipeline that parses eligible Deep Research summaries into canonical `benchmark_payload` JSONs to accelerate closing the ~150-paper backlog.
- [ ] Include a strict human-in-the-loop review interface to guarantee the 8-point SLR criteria are strictly maintained before merging into the main pipeline.

#### S12.3: Model-Guided Active Learning (DoE Feedback Loop)

- [ ] Formalize the "Structural Gaps" into explicit Design of Experiments (DoE) workflows.
- [ ] Implement an API so that when the system identifies a critical gap (e.g., lack of MFT/FFT data in SPI extrudates), it auto-generates a precise wet-lab protocol optimised for maximal model calibration gain.

### Deferred Scientific Modeling Backlog

Status 2026-04-08: keep `5.7`, `5.10`, and `5.11` deferred until the first primary matrix and extrusion benchmark panels exist. Treat `5.8` as the only immediate follow-on once the SPI extrusion benchmark produces same-run SH, furosine, and volatile data.

#### 5.7 Bidirectional Lipid-Maillard Crosstalk

- [ ] Add dicarbonyl-lipid oxidation promotion pathway in `lipid_oxidation.py`.
- [ ] Add melanoidin antioxidant capacity as a time-dependent LOPs suppressant.
- [ ] Validate against Report 11 crosstalk heuristics.

#### 5.8 Disulfide Bond Evolution / MFT Retention

- [ ] Model free-SH to disulfide kinetics as a function of SME and temperature.
- [ ] Link that state variable to MFT headspace recovery in the volatile retention model.

#### 5.10 Sunflower Chlorogenic Acid Off-Note

- [ ] Add temperature-gated 4-vinylguaiacol penalty for sunflower-containing formulations.
- [ ] Include chlorogenic acid to lysine covalent adduct formation as a lysine-accessibility sink.

#### 5.11 Transport / Diffusion Model for Volatile Release

- [ ] Design a 1D Fickian diffusion slab model for volatile release during cooling or serving.
- [ ] Integrate it with volatile retention factors as a compound-class-specific alternative to scalar correction.

### S9. Skipped Test Triage & QM Optionality

Status: supporting infrastructure, not current product bottleneck.

- [ ] Resume only after S14 if skipped-test cleanup blocks deterministic confidence in the active scientist workflow.

#### S9.1: Inventory and classify skipped tests

- [ ] Build a machine-readable skip registry from `pytest -rs` grouped by reason, file, and dependency class.
- [ ] Classify skips into: `not_implemented_module`, `missing_external_dataset`, `missing_optional_backend`, `long_running_campaign`.
- [ ] Add owner and unblock criteria for each skip cluster.

#### S9.2: Quasi-harmonic correction decision and implementation path

- [ ] Implement `src/quasi_harmonic_correction.py` with pure-Python deterministic numerics and no heavy backend coupling.
- [ ] Replace unconditional skips in `tests/benchmarks/test_quasi_harmonic_correction.py` with executable deterministic tests plus `skipif` only for optional integration points.

#### S9.3: Barrier and IRC benchmark skip policy

- [ ] Keep Phase 3.3/3.4 benchmark tests gated by explicit dataset/backend markers rather than unconditional `pytest.skip` inside the test body.
- [ ] Encode a run contract: default CI lane runs deterministic/unit and mock-backed checks; optional QM lane runs benchmark suites when datasets are mounted.

#### S9.4: DFT and ML-potential complement policy for skipped-test conversion

- [ ] Tie each formerly skipped QM test cluster to one of three execution lanes: `deterministic_helper_lane`, `optional_mlp_acceleration_lane`, `optional_dft_authority_lane`.
- [ ] Ensure report and validation artifacts expose lane provenance so users can distinguish deterministic numerics, MLP-assisted predictions, and DFT-confirmed evidence.

### P3–Refinement. Selective Mechanistic Refinement

Active only for matrix benchmarks that remain `mechanistic_priority` after observable closure review. Do not expand broad xTB or DFT activity beyond the specific targets in S18.

### P4. MLP Adoption

Offline accelerator lane only. External MLP evaluation must not substitute for missing matrix benchmarks or observable anchors.

### P6. Matrix-Family Expansion Beyond Pea and Soy

Keep matrix-family coverage explicit in artifacts. Do not broaden family-level scope faster than the evidence surface can support.

### Reproducibility & Provenance

- [ ] Add version-pinned reproducibility snapshots (`pip freeze` equivalent or Docker image tag) to every report artifact so predictions are exactly reproducible for peer review and regulatory contexts.
- [ ] Document the exact mapping between report provenance metadata and the reproducibility snapshot.

---

## Current Product Status

### Strong today

- Quantitative screening is inside the 1.5x band across the current executable surface: 13/13 experimental and 14/14 total quantitative benchmarks, with 11 strict-ready benchmarks and worst experimental ratio 1.442x.
- 16 chemistry families are tracked; 7 are benchmark-linked and 4 currently expose compound-level quantitative parity.
- 54 quantitative compound points are now summarized in the validation surface; the core Maillard family sits at median ratio 1.025x and mean |log10 error| 0.026.
- Family-aware ingestion, runtime, validation, and reporting are operational.
- Pea and soy matrix paths are executable and useful for directional prioritisation.
- README-facing validation imagery is regenerated from the official validation workflow and copied into `docs/assets/`.
- Family 03 dilute mixed-loading closure has been re-verified in Docker after landing.
- Trust language, evidence posture, and family-lane transparency are already visible in reports.
- Intake-registry state model is normalised with explicit `triage_status`, `encoding_status`, and `runtime_artifacts_present` fields.
- Canonical literature backlog artifact exposes three exclusive queues: `ready_runtime`, `ready_benchmark`, `wet_lab_blocked`.
- Extrusion process modeling with SME coupling, moisture-regime bifurcation, and pre-extrusion damage baselines is operational.
- Safety markers (acrylamide, furosine, CML/CEL, LAL) are integrated into the pipeline.
- Cheap-first refinement screening pipeline (barrier offset sweeps with benchmark-visible decision gating) is operational.
- DoE generator templates exist for 5 gap types (blocking benchmark, missing anchor, missing kinetic, missing process-state, missing flavor anchor).
- SLR covers 990 lines with 16+ benchmark-eligible references and 5 structural gaps formally identified.

### Still blocking scientist value

- **No formal uncertainty quantification**: predictions carry heuristic trust tiers instead of explicit ±bounds, so scientists cannot distinguish "I can act on this" from "this is a guess" (S20.1).
- **No matrix benchmark** is yet external-decision-ready; status is `computationally_narrowed` at best.
- **No data ingestion path for non-experts**: a scientist with GC-MS results has no single-command way to feed data back into the model (S21).
- **No experiment prioritization**: the repo knows its gaps but doesn't tell scientists which measurement would improve predictions the most (S22).
- **Scientist-facing visual output** is still too thin in the main report surface: no intervention waterfall, confidence overlay, capabilities heatmap, or presentation-grade parity export (S13).
- **No web interface** — all interaction requires Docker + command-line, a major adoption barrier for food scientists (S19).
- **Extrusion modeling is not yet benchmark-calibrated enough for release-facing use** even though the process-state lane exists (S17).
- **No cross-validation infrastructure**: transfer assumptions have not been validated via leave-one-out or holdout schemes (S20.2).
- Mixed pea and soy meaty-positive targets still rely on transferred or internal-candidate observable support.
- Six chemistry families still have zero benchmark-linked closure.
- `literature_runtime.py` (4,324 lines), `benchmark_validation.py` (2,611 lines), and `recommend.py` (1,338 lines) are monoliths that impede contribution and debugging (S14).
- Test suite runtime (~1h40m) impedes development iteration (S14).
- No version-pinned reproducibility or benchmark versioning for peer review (S21.3).

---

## Success Criteria

- [ ] A scientist can use the tool to narrow a wet-lab or literature-backed matrix campaign, not just inspect simulated chemistry.
- [ ] Free-precursor predictions remain quantitatively stable while matrix usefulness improves.
- [ ] At least one matrix lane becomes meaningfully closer to external decision readiness.
- [ ] Reports and artifacts make promotion blockers explicit enough that synthetic closure is difficult.
- [ ] Expensive compute stays offline, sparse, and justified by benchmark-visible decisions.
- [ ] Visual output makes predictions actionable without requiring the scientist to interpret raw JSON.
- [x] All benchmarks pass inside the 1.5× acceptance band.
- [x] A food scientist can understand what the tool does and run a first prediction within 10 minutes of reading the README.
- [ ] The runtime parametric surface grows by ≥2× through literature extraction before any wet-lab work.
- [ ] Multi-objective tradeoffs (meaty vs. safe vs. off-note) are visible as Pareto frontiers, not collapsed scores.
- [ ] Every prediction carries explicit uncertainty bounds (90% CI) from Monte Carlo propagation (S20.1).
- [ ] A scientist can ingest new experimental data with a single command and see its impact before committing (S21).
- [ ] The tool can rank candidate experiments by expected calibration gain and auto-generate protocols (S22).
- [ ] Transfer assumptions are validated via leave-one-out cross-validation across the entire benchmark set (S20.2).

---

## Completed Foundations

Sprints S0–S11 are complete as of 2026-04-05. All foundational architecture, family-aware ingestion and runtime, matrix observable closure, scientist experiment intake loop, family promotion contracts, intake-registry normalisation, Deep Research runtime queue publishing, extrusion process modeling (SME/moisture), safety marker integration (acrylamide, furosine, CML/CEL, LAL), cheap-first refinement screening, and DoE generator templates have been implemented and validated.

---

## Active Plan — 2026-04-21 SOTA TS pipeline + elementary-step decomposition (method B)

### Contexto
- `pe_amadori` con un solo TS dry: 73.6 kcal/mol (gap ~54 vs lit 19.8).
- Mecanismo real: dos transferencias 1,3-H acopladas (O0→C2, C1→N3). DE-GSM no encuentra TS razonable porque busca un único saddle para dos eventos simultáneos en un anillo de 3.
- Intentos con agua explícita (4 iteraciones, ver `memories/repo/pe_amadori_water_shuttle_attempt_2026_04_21.md`): chocan por geometría (waters de R chocan con N3 desplazado en P) y por pyGSM TRIC en sistemas con fragmentos desconectados.
- **Decisión estratégica del usuario**: implementar opción B (decomposición secuencial) y hacerlo *general*, no ad-hoc. Además: revisar si el pipeline de búsqueda de TS está al estado del arte.

### Auditoría del estado del arte (¿estamos SOTA?)

**Lo que tenemos**:
| Etapa | Implementación | Estado SOTA |
|---|---|---|
| Endpoints | RDKit ETKDG (1 confórmero) + xTB opt | ⚠️ No ensemble. Riesgo de basin-hopping. |
| Path search | xTB `--path` + DE-GSM (fallback) | ✅ DE-GSM es SOTA para concerted; ⚠️ no CI-NEB ni AFIR |
| TS guess | Pico de energía del string | ⚠️ Sin validación n_imag intermedia |
| TS refine | Sella (eigenvector following) + v0 de xTB Hessian | ✅ Cerca del SOTA |
| Validación | n_imag, magnitud, IRC opcional (proxy ±desplazamiento) | ⚠️ IRC no analítico; sin visualización del modo |
| Solvación | ddCOSMO/ALPB implícito; CREST-QCG existe pero **no integrado** | ⚠️ Para protones móviles, cluster-continuum es gold standard |
| Multi-step | **No existe**; cada `target` = 1 TS | ❌ Bloqueante para reacciones acopladas |

**Gaps relevantes para Maillard** (proton transfers, reacciones acopladas, agua catalítica):
1. **Decomposición de pasos elementales** — bloqueador #1 para `pe_amadori`. Sin esto, ninguna mejora del refinamiento ayuda.
2. **Cluster-continuum solvation generalizable** — la lógica ad-hoc de `generate_pe_amadori_water_seed.py` debe promoverse a módulo.
3. **Mode-coordinate sanity check** — comparar el vector imaginario del TS con los átomos esperados (donor, acceptor, H transferido). Cheap, pero hoy no existe.
4. **CI-NEB como segundo path-method** — útil cuando GSM falla. Disponible en ASE (`ase.neb.NEB` con `climb=True`).
5. **GP-NEB / ML surrogate** — no necesario hoy (los pasos elementales son baratos en xTB).

### Plan de implementación

#### Fase A — Limpieza (preparar terreno)
- [ ] **A.1** Identificar y borrar dead code:
  - `scripts/_tmp_trim_recommend.py`, `scripts/get_benchmark_lipids_tmp.py` (one-offs evidentes por nombre)
  - `scripts/generate_explicit_water_seed.py` (será reemplazado por módulo general)
  - `scripts/generate_pe_amadori_water_seed.py` (idem; conservar como referencia hasta que módulo nuevo esté validado)
  - Verificar `src/diffusion_ts.py`, `src/mlp_barrier.py` — si son scaffolding sin uso real, borrar
  - Auditar `calculate_barrier()` vs `calculate_robust_barrier()` en `src/dft_refiner.py`: si el primero ya no se llama desde nada activo, eliminarlo (solo quedará robust_sp como camino oficial)
- [ ] **A.2** Tests siguen verdes (pytest unit lane).

#### Fase B — Microsolvation generalizable
- [ ] **B.1** Crear `src/microsolvation.py`. API:
  ```python
  def place_proton_shuttle_water(coords, syms, donor_idx, acceptor_idx, h_idx, *, oh_bond=0.97, hoh_deg=104.5, offset=1.6) -> tuple[np.ndarray, list[str]]:
      """Place 1 explicit water bridging donor→H→acceptor in a 5-membered TS topology."""
  
  def place_water_shell(coords, syms, hot_atoms, n_waters, *, frozen_solute=True) -> ...:
      """Generic water placement: drop n waters near `hot_atoms`, seeded by molecular normal + perturbation."""
  
  def kabsch_align(P, Q) -> tuple[rot, t_p, t_q]:
      """Heavy-atom RMSD alignment (move into shared utils)."""
  
  def relax_solvent(xyz_path, solute_atoms, *, level='xtb-alpb', work_dir) -> Path:
      """Constrained xTB opt: solute frozen, waters relaxed."""
  ```
- [ ] **B.2** Tests unitarios: water geometry sanity (OH=0.97±0.02, HOH≈104.5±2°, no clashes).

#### Fase C — Decomposición de pasos elementales (corazón de método B)
- [ ] **C.1** Schema en `data/lit/computational_gap_closure_targets.json`:
  ```json
  {
    "reaction_key": "pe_amadori_decomposed",
    "type": "multi_step",
    "elementary_steps": [
      {
        "step_id": "step1_oh_to_c2",
        "name": "1,3-H shift O0→C2 (water-shuttled)",
        "reactant_xyz": "data/geometries/.../step1/reactant.xyz",
        "product_xyz":  "data/geometries/.../step1/product.xyz",
        "donor_atom": 0, "acceptor_atom": 2, "h_atom": 17,
        "microsolvation": {"n_waters": 1, "topology": "proton_shuttle"}
      },
      {
        "step_id": "step2_ch_to_n3",
        "name": "1,3-H shift C1→N3 (water-shuttled)",
        ...
      }
    ],
    "literature": {"barrier_kcal_mol": 19.81, "source": "..."},
    "aggregation": "max_step"
  }
  ```
- [ ] **C.2** `src/elementary_step_runner.py`:
  ```python
  class ElementaryStepRunner:
      def run_step(self, step_spec) -> StepResult:
          # 1. apply microsolvation
          # 2. xTB --path or DE-GSM
          # 3. xTB-Sella refine
          # 4. validate n_imag + mode-coordinate sanity
          # 5. DFT SP (robust_sp)
          # 6. return ΔG‡_step
      
      def run_target(self, target_spec) -> TargetResult:
          # iterate elementary_steps; aggregate
  ```
- [ ] **C.3** Wire en `scripts/run_computational_gap_dft.py`: detectar `type: multi_step` y dispatch a `ElementaryStepRunner.run_target()` en lugar del flujo single-step.
- [ ] **C.4** Mantener compatibilidad: targets sin `type: multi_step` siguen al pipeline existente.

#### Fase D — Mode-coordinate sanity check (cheap SOTA win)
- [ ] **D.1** En `src/xtb_backend.py`, función nueva:
  ```python
  def validate_ts_mode(coords, syms, imag_vector, expected_atoms, *, threshold=0.3) -> dict:
      """Check that |imag_vector[expected_atoms]| / |imag_vector|_total > threshold.
      
      Returns {'pass': bool, 'concentration': float, 'top_atoms': [...]}.
      """
  ```
- [ ] **D.2** Llamar en `ElementaryStepRunner.run_step()` después de validar n_imag. Si falla → marcar el TS como `mode_mismatch` (no es el TS esperado aunque sea saddle válido).

#### Fase E — Aplicar a pe_amadori
- [ ] **E.1** Definir `pe_amadori_decomposed` con dos pasos:
  - Step 1: Schiff base hidratado → enol intermediario (H17 transfer O→C2)
  - Step 2: enol → Amadori (H18 transfer C1→N3)
- [ ] **E.2** Ejecutar `python scripts/run_computational_gap_dft.py --target pe_amadori_decomposed --strategy robust_sp --execute`.
- [ ] **E.3** Verificar: cada paso n_imag=1, mode-coordinate-sanity OK, ΔG‡_max comparado con lit 19.8 kcal/mol.

#### Fase F — Cierre
- [ ] **F.1** Tests unitarios: mock target multi-step + step runner.
- [ ] **F.2** Update README "validation" section con resultado pe_amadori antes/después.
- [ ] **F.3** Memoria de repo: `pe_amadori_decomposed_validation_2026_04_2X.md`.

### Criterios de éxito
1. ✅ `pe_amadori_decomposed` ΔG‡_max dentro de 5 kcal/mol de lit 19.8 (vs gap 54 actual).
2. ✅ Mismo schema funciona para registrar futuras reacciones acopladas (ej. `aa_ring_open_dicarbonyl` si tiene mecanismo multipaso).
3. ✅ Suite unit tests sigue verde (≥531 pass).
4. ✅ ≥3 archivos/scripts dead code eliminados.
5. ✅ Validación n_imag + mode-coordinate sanity en cada TS.

### Lo que NO hacemos en este sprint (parking lot)
- CI-NEB como segundo path-method (pyGSM ya cubre bien si endpoints están bien planteados).
- ML-NEB / GP-NEB (overkill para pasos elementales en xTB).
- Ensemble de confórmeros (CREST conformational sampling) — defer hasta que sea bottleneck demostrado.
- Refactor de `dft_refiner.py` (2980 líneas) — bigger fish a aguas fuera del foco.

