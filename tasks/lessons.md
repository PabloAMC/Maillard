# Lessons Learned

## Kinetic Modeling

- **Concentration Scaling**: Always standardize units (e.g., ppb) across the entire pipeline. Mixing molar ratios, ppm, and ppb causes detection thresholds to fail and scoring to zero out.
- **Boolean Guarding**: When calculating ratios (like lysine budget or masking factors), always guard against zero-sum denominators (`if sum > 0: ... else: 0.0`).
- **Detection Thresholds**: In scientific benchmarking, ensure that safety scoring and off-flavor detection use sufficiently sensitive thresholds (e.g., 1e-20 ppb) to capture the effects of barrier-shifting interventions (like calcium carbonate) even at high temperatures where formation is rapid.

## Architecture & Data Structures

- **Dataclass Defaults**: Adding required fields to shared dataclasses like `FormulationResult` breaks existing unit tests that use positional arguments or incomplete mocks. Always provide sensible defaults for new fields.
- **Redundant Scaling**: Avoid "stacking" Boltzmann factors. If the recommendation engine already accounts for kinetic probability in its concentrations, the scoring layer should use them directly rather than re-applying the exponential.
- **SMILES vs Common Names**: Scientific benchmarks and sensory modules often use fuzzy name matching. Using common names as keys in `predicted_ppb` dictionaries is often more robust than raw SMILES, provided the pipeline's canonicalizer (e.g., `_canon`) handles non-SMILES strings gracefully (by returning the original string with a warning).

## Sensory Modeling

- **Unit Propagation**: Psychophysical models (Stevens' Law) are extremely sensitive to input units. A 1000x error (ppm vs ppb) can drop concentrations below the odor detection threshold (ODT), causing `KeyError` or zero scores.
- **Mocking Robustness**: When mocking databases for sensory tests, ensure all internal maps (`smiles_map`, `chemical_to_descriptor`) are populated consistently with the implementation's lookup strategy.
- **Intervention Transparency**: Ensure that high-level evaluation methods (like `evaluate_single`) correctly propagate all formulation attributes (especially `interventions`) to the underlying evaluation loops, otherwise the resulting scores will lack scientific sensitivity.

## Rule-Based Engines

- **Label Robustness**: Exact label matching for chemicals (e.g., "leucine" vs. "l-leucine") can cause silent failures in rule-based engines. Always normalize labels to a canonical, case-insensitive form before dictionary lookups to ensure all intended pathways fire correctly.

## Verification & Benchmarking

- **Soundness over Shortcuts**: When a regression test fails, avoid "fixing" it with ad-hoc string exclusions (e.g., `if "=n" in name`). Instead, identify the underlying physical property (e.g., volatility defined by MW and H-bond capacity) that justifies the exclusion. This ensures the test remains a valid scientific benchmark rather than a self-fulfilling prophecy.
- **Equilibrium Dynamics**: In closed-system kinetic models (like Cantera batch reactors), irreversibility is a dangerous simplification for Schiff bases. Always model these as reversible equilibria to avoid "trapping" volatiles at near-zero concentrations, which causes false negatives in aroma prediction.
- **Metric Validity**: Never report correlation metrics like Pearson R for fewer than 3 matched compounds. With two points they are mathematically misleading and can mask obvious calibration failures.
- **Workflow Discipline**: After a user correction about process, update the working plan and the lessons file immediately before continuing implementation.
- **Environment Assumptions**: If the user states that an external tool is a hard requirement, do not silently add behavioral fallbacks that change outputs. First verify the tool path, then fix the environment or fail explicitly with a clear contract.
- **LaTeX Contract**: For Matplotlib science plots, “latex available” means both the binary and the required TeX packages are usable by a real compile. Check the full toolchain via one shared preflight helper; do not trust `which latex` alone.
- **Priority Reassessment**: If the user says another model or chat already advanced the work, re-audit priorities from repo state (`git status`/diff, backlog, repo memory) before answering. Do not trust the previous summary alone.
- **Artifact-Driven Replanning**: Before keeping an old validation-closure sprint active, re-read the current quantitative artifacts. If high-ratio tails and broad residual errors are already gone, shift the roadmap from error triage to the next real product bottleneck, usually matrix promotion logic, observable closure, or evidence ingestion.
- **Execution Environment**: In this repository, verification must run inside Docker using the `maillard` conda environment on Python 3.12 when the user requests it. Do not fall back to host Python tools for scientific validation.
- **Step Framing**: When the user asks for it, state the difficulty of the next substantial step and whether GPT-5.4 or GPT-5 mini should handle it before proceeding.
- **Task Annotation**: If the user asks for model recommendations, write the GPT-5.4 vs GPT 5 mini split directly into the active plan or backlog artifact, not only in chat updates.
- **Docker-First Triage**: Before changing scoring or QM code, reproduce the exact failing subset inside the Docker `maillard` environment. If the subset is already green there, do not edit production logic just to match a stale host-side failure report.
- **Scientific Inspection Ergonomics**: If a benchmark diagnostic requires long `python -c` or heredoc commands through the Docker wrapper, promote it into a named script or wrapper subcommand. Quote-heavy ad hoc commands are not a reproducible interface.
- **Documentation Scope**: When the user asks for a structural repository review for newcomers, do not answer with a short doc touch-up. Rework the entry layer explicitly: README, quickstart, trust-and-limitations guidance, command map, and a clear grouping of research versus development notes.
- **Backlog Encoding**: When the user asks to turn a strategic recommendation into action, do not leave it in chat. Update the canonical docs and encode the actual execution sequence, constraints, and decision rules in `tasks/todo.md`.
- **Reference Naming Consistency**: When validation figures, README visuals, or reporting artifacts surface literature labels, normalize them to one citation style consistently (for example, `Author et al. (Year)`). Do not mix author-year labels with raw PMC IDs, ACS shorthand, or tracker-specific aliases in the same user-facing surface.
- **Descriptive File Naming**: Do not introduce or preserve opaque shorthand like `p3` or `family03` in new filenames when a descriptive scientific scope is available. Prefer names that explain the benchmark or workflow under test.
- **Execution Lane Naming**: Do not invent new user-facing script, artifact, or command names with opaque phase labels like `wave1`, `p3`, `p4`, or chemistry shorthand like `c4_c5`. Name them after the scientific job they perform, for example `computational_gap_dft_ingestion_report` or `refinement-governance`.
- **Roadmap Pruning Discipline**: When cleaning `tasks/todo.md`, compress completed slices aggressively but do not delete deferred strategic backlog that still encodes product direction. Preserve long-horizon modeling lanes under an explicit deferred section instead of dropping them.
- **No-Wet-Lab Constraint**: If the user clarifies that no wet-lab loop is available, rebuild the plan around literature-closable surfaces, machine-readable evidence maps, and selective offline refinement. Separate structural gaps that truly require new primary data from gaps that computation or curation can still reduce.
- **Active Job Artifact Safety**: Before deleting ignored outputs under `results/` or other runtime artifact directories, check for active long-running processes that use those paths for heartbeats or incremental persistence. If a job is live, preserve its output directory or coordinate a restart window first.
- **QM Execution Ownership**: If a legacy DFT process still owns a target execution JSON, stop it before writing a new preflight or restart artifact for that same target. Otherwise the old heartbeat can overwrite the repaired state and make a blocked or ready seed look ambiguous.
- **SLR-to-Repo Reaudit**: If the user asks to step back after literature updates, do not continue from the previous roadmap blindly. Re-read the active SLR, compare it against the machine-readable calibration artifacts and user-facing reports, and rebuild the backlog around the delta between literature data, encoded priors, and exposed outputs.
- **Optional Backend Imports**: For QM helpers that depend on optional backends like Sella/JAX, catch runtime import failures as well as missing-package errors. Unsupported CPU features must degrade to an explicit skip or fallback path, not a pytest collection error.
- **Projection Stability**: Do not use raw FAST pathway weights to allocate absolute volatile budget. They are useful for ranking but can collapse the projection onto sulfur donors or other non-benchmark carriers. Keep the budget projection based on stable endpoint filters and use weights only if they are explicitly re-scaled and benchmark-calibrated.
- **Temporal Consistency**: Do not calibrate benchmark branches until the main evaluation path and the diagnostic path propagate the same `time_minutes` value into FAST prediction. Otherwise the apparent branch errors can be artifacts of inconsistent temporal semantics rather than chemistry.
- **Stop Local Recalibration Loops**: If Docker sweeps show that no single local barrier tweak improves Hofmann, Mottram, and Farmer together, stop tuning sulfur barriers and move the work to the concentration/headspace projection layer. That is a structural blocker, not a chemistry-family constant problem.
- **Quantitative Visibility**: If a validation overview claims to summarize quantitative performance, do not filter out matrix-only or matrix-augmented numeric benchmarks. Show the full quantitative surface and differentiate trust tiers with visual encoding, not silent omission.
- **Calibration Single-Application Rule**: When `HeadspaceModel.get_matrix_benchmark_headspace_factor()` already includes the matrix observable calibration, downstream benchmark code must not multiply `calibration_observable_factor` a second time. Keep `proxy`, dynamic release, and total observable scaling as separate quantities in metadata.
- **No Synthetic Closure**: If external mixed matrix benchmarks do not exist, do not treat internally constructed mixed benchmarks or directional support as promotion closure. Re-route those cases into mechanistic triage or explicit external-data blockers, and encode that distinction in artifacts and roadmap text.
- **Report Variable Scope**: In `generate_report(...)`, keep shared markdown locals (for example `calibration`) initialized outside conditional blocks; regression tests often use sparse `FormulationResult` fixtures without `confidence_metadata`, and uninitialized locals can crash report rendering.
- **Figure Review Discipline**: Before considering a validation figure done, inspect visible outliers and legend labels together. A quantitatively correct artifact can still be misleading if a known reference-only outlier dominates the plot or if benchmark labels mix naming conventions across papers, internal harnesses, and protocol pilots.
- **Script Deletion Safety**: A grep scan finding zero in-code references to a generator script does not prove it is unused. Users invoke scripts via ad-hoc CLI commands (`docker_maillard.sh run "python scripts/generators/..."`) that leave no trace in the codebase. Before deleting any script, confirm with the user or wire them into `docker_maillard.sh` so their usage is discoverable.
- **Multi-Fragment Geometry Embedding**: When generating 3D coordinates from SMILES containing `.` (disconnected fragments), RDKit's `EmbedMolecule(ETKDG())` can place fragments on top of each other, producing sub-angstrom interatomic distances. Always embed and MMFF-optimize each fragment independently, then combine with explicit spacing. Do not optimize after combining — MMFF can collapse fragments back together.
- **xTB Path Validation**: A completed xTB path search (returncode 0) does not guarantee a real reaction path was found. Always check that the pairwise RMS delta between reactant and TS geometries is >> 0.01 Å. If it is ~1e-7 Å, the path is an identity transformation and the reported barrier is numerical noise, not chemistry. The DFT preflight gate catches steric clashes but not degenerate paths where the reactant geometry passes but xTB returns it unchanged.
- **Multi-Frame Trajectory Files**: xTB `--path` writes `xtbpath_N.xyz` files as multi-frame trajectories (each containing the full optimization history for path point N). The FIRST frame is always the initial/reactant geometry; the LAST frame is the optimized geometry for that path point. Any code that reads these files must extract the last frame, not just call `_parse_xyz()` which reads the first frame. Failing to do this causes the DFT repair logic to treat every candidate as identical to the reactant, silently selecting a worthless TS guess.
- **Minimum Displacement Gate**: A TS guess that passes all "not too distorted" checks can still be useless if it is identical to the reactant. Always enforce a minimum pairwise RMS delta (≥ 0.3 Å) between the TS guess and the reactant. Without this, the system can spend 10+ hours computing a meaningless barrier of ~1 kcal/mol that is just SCF numerical noise.
- **Defense-in-Depth for Expensive Compute**: Before committing to multi-hour DFT calculations, run every cheap sanity check possible: (1) minimum displacement in preflight, (2) energy-based TS selection from xTB path, (3) pre-execution RMSD sanity check (< 0.5 Å = abort). Each gate costs milliseconds but can save hours.

### Intermolecular Pre-Docking Constraints
- **Pattern**: When generating multi-fragment coordinate seeds for xTB, simple translation often leads to huge physical gaps between reactive centers (e.g., 9.95 Å) because internal geometry dictates fragment facing. This breaks xTB's linear topological interpolation.
- **Correction**: Do not enforce temporary dummy topological `Chem.Bond` features to pull fragments together, as this breaks RDKit's implicit valence logic and fails embeddings silently. 
- **The Correct Way**: Use proper forcefield properties natively via `AllChem.MMFFGetMoleculeForceField`. Apply `ff.MMFFAddDistanceConstraint()` to softly pull specific reactant atoms to a close Van der Waals distance (~3.0 Å), and sequentially call `ff.Minimize()`. This organically rotates entire molecule fragments like rigid bodies into their correct pre-reactive topology.

### ML Force Field Version Pinning
- **Pattern**: Installing optional ML libraries (e.g., UBio-MolFM, nequip-allegro) via `pip install` can silently upgrade shared deep learning dependencies like `e3nn` even though they appear unrelated.
- **Correction**: After any experimental `pip install` in the Docker container, verify the MACE calculator loads (`from mace.calculators import mace_off`) before considering the environment stable. If it breaks, downgrade `e3nn` back to `==0.5.1`.
- **The Correct Way**: Pin exact compatible versions for critical shared deps directly in `environment.yml`. For MACE 0.3.x, this means `e3nn==0.5.1`. Never leave it unpinned in an environment that also installs nequip or other e3nn-dependent libs.

### DFT Geometry Checkpointing
- **Pattern**: Long DFT geometry optimizations (8–12 hours) were being lost silently when the downstream solvent SCF or single-point step failed – the optimized XYZ geometry existed only in `tempfile.TemporaryDirectory` which was wiped on process exit.
- **Correction**: Checkpoint the optimized reactant/TS geometry to a deterministic path on disk immediately after `optimize_geometry()` returns, before the solvent SCF is attempted.
- **The Correct Way**: Pass a `checkpoint_dir` to `calculate_barrier()` and call `_save_checkpoint()` after every major optimization phase. This guarantees that even if the final energy evaluation fails, the expensive geometry is recoverable and can be used as a warm start on restart.

### Solvent SCF Convergence Robustness
- **Pattern**: A valid optimized TS geometry (19h of compute) was discarded because the final `ddCOSMO` solvent SCF step failed to converge – a one-shot attempt with a hard `RuntimeError`.
- **Correction**: Implement a 3-strategy retry ladder: (1) standard solvent SCF → (2) SOSCF + solvent → (3) hardened vacuum SCF + level-shift + SOSCF. Only raise if all 3 fail.
- **The Correct Way**: The vacuum fallback loses solvent correction energy (~0.5–2 kcal/mol accuracy hit) but preserves the most expensive result – a valid TS geometry and its energy – rather than discarding all compute.

### Pre-Launch Code Verification
- **Pattern**: An edit to fix a non-fatal pre-optimization fallback introduced broken indentation that made the entire SOSCF retry block syntactically wrong, which would have crashed every single job at that recovery branch.
- **Correction**: After any multi-chunk edit, always run `python -m py_compile` on every changed file before launching long-running jobs.
- **The Correct Way**: Treat expensive multi-hour compute jobs like production deployments. Run a syntax check + smoke import test as the final gate before launch. This takes 30 seconds and can save days of re-running.

### React-OT Inputs Need Atom-Mapped Reactant/Product XYZs
- **Pattern**: Running React-OT (deepprinciple/react-ot) on our `data/geometries/xtb_inputs/<target>/{reactant,product}.xyz` pairs produced the warning `RMSD is too large, please double check the input xyz files` followed by a NaN-position cascade for every CHON-eligible target (`lysinoalanine_crosslink`, `aa_ring_open_dicarbonyl`, `pe_schiff_base`, `asparagine_sugar_explicit_water_cluster`). All four returned `missing_ts_output` despite a clean Colab GPU run.
- **Root cause**: React-OT requires reactant/product XYZ files where atom *i* in R corresponds to atom *i* in P. Our `reactant.xyz` and `product.xyz` are produced by independent docking + relaxation (RDKit + xTB), so atoms are never in matched order. Same atom count is necessary but not sufficient.
- **The Correct Way**: Before sending a target to React-OT, generate atom-mapped R/P pairs (e.g. via an RDKit reaction template `RunReactants` → align indices, or RXNMapper, or read mapped SMILES from curated_pathways and re-embed both endpoints from the mapped molecule). Do not feed independently-relaxed geometries to React-OT.
- **Operational corollary**: React-OT remains a pathfinder; if atom-mapping cannot be produced for a target, `react_ot_seed_coverage.md` must record `missing_seed` and the xTB → DFT lane proceeds unchanged. Never silently downgrade a barrier on a NaN'd React-OT output.

### Colab Notebook Hardening for External Repos
- **Pattern**: A single Colab notebook for an external repo (React-OT) hit five distinct import failures in sequence: missing module after `pip install -e .` (broken upstream packaging), `numpy.rec` missing (forced numpy<2 broke Colab's scipy), `torch_scatter` ABI mismatch (forced torch reinstall on top of Colab's preinstalled torch), `ase.neb` missing (ASE moved NEB to `ase.mep`), `ipdb`/`wandb` missing (upstream imports them at module load).
- **The Correct Way**: For any third-party repo on Colab — (a) NEVER reinstall torch; derive the PyG wheel index from the preinstalled `torch.__version__`. (b) NEVER pin numpy/pandas; Colab's scipy is built against numpy 2.x. (c) Do the editable install with a fallback to `sys.path` injection because some upstream packages have broken `setup.py`/`pyproject.toml`. (d) Install the obvious dev-dep imports the upstream module-loads (here `ipdb`, `wandb`). (e) Add compat shims (`sys.modules['ase.neb'] = shim`) for any deprecated import paths the upstream still uses. (f) Always print a diagnostic dict (`torch`, `cuda`, `numpy`, `scipy`, `torch_scatter` versions) so failures get attributed correctly.
- **Workflow corollary**: Colab pulls notebooks live from GitHub — local-only patches are invisible. Push the fix to the branch the user opens.


## 2026-09-01 — data restructure (cleaning branch)

- **A gitignored directory under `data/` is invisible to every audit, gate and reviewer.** `data/qm/` shipped 18
  unsourced barrier windows for four months because `data/*` was ignored wholesale. Rule: `.gitignore` whitelists
  data files explicitly; `git ls-files -i -c --exclude-standard` must be empty.
- **Silent fallbacks hide broken installs as bad science.** `HeadspaceModel()` resolved its constants relative to the
  CWD and fell back to `Kaw = 0.01` for every compound with no message; `_load_ramp` logged a warning and ran
  isothermally; a dozen loaders returned `{}` for a missing registry. Rule: curated files are loaded through
  `src/data_access.py`, which raises. A test that points a loader at a nonexistent path must expect an exception.
- **Tests that write into the repo tree accumulate.** 102 `calibration_history` files were pytest droppings. Rule:
  `scripts/ci/data_readonly_gate.py` runs after every tier; output dirs are monkeypatched to `tmp_path`.
- **Regeneration deletes hand corrections.** The hold-out bundles carried primary-source corrections that the
  generator's spec never learned; a regeneration would have reverted them (and once did, for a citation note).
  Rule: generators of frozen evidence are write-once and carry curated-only keys forward on `--overwrite`.
- **Free-text joins need a registry, not a fifth alias table.** 148 spellings for ~85 compounds, five alias tables,
  one `difflib` fallback at 0.75. Rule: `data/keys/compounds.yml` + `compound_keys.resolve`; every spelling under
  `data/` must resolve (tested); new spellings are added to the generator seeds, not to code.
- **Scripted edits need an assertion per replacement.** A `str.find` that returned `-1` silently duplicated a test
  file; a whitespace mismatch silently skipped an edit. Rule: assert `old in text` before every replace, compile
  after every batch, and re-read the tail of a test log before committing (a Phase 1a commit claimed a pass with two
  failures in the log).
- **A rule's positive control is a data audit.** Writing the 10-hydroperoxide-to-1-octen-3-ol rule (2026-09-08)
  failed its own control because the registry's 1-octen-3-ol SMILES had seven carbons while its InChI had eight
  (`off_flavour_targets.yml`, copied into the generated registry). Rule: every structure a rule names is checked by
  RDKit against an independent identifier (InChI or atom count) when it enters; a curated SMILES is data, not truth.

- **A pre-registered ranking comparison should ask about a direction, not a pair.** The cultivated-tissue
  invariance test (2026-09-13) declared T1 as "one reversal carries ≥ 50 % of the disagreements". The engine
  disagreed with the naive rule in 40 % of draws, robustly (100 % envelope survival), and in one consistent
  direction (glucose demoted every time), but the direction surfaced in whichever pair happened to be restorable
  in the draw, so no pair reached half and the verdict fell to indeterminate. Rule: when a comparison is between
  rankings, the concentration clause names a *precursor* that moves consistently in one direction across draws,
  never a recurring *pair*. Also: before declaring the candidate set, ask the engine which names it takes and
  which lanes compose; three of eight declared candidates and the whole trunk arm were refused by construction,
  and the thresholds had been sized for a set the engine could not rank (section 6 of the pre-registration).

- **An automated read of a table is a claim about the table, not the table.** The 2026-09-14 composition read wrote
  seven dossiers from fetched HTML and abstracts. Checked against the PDFs the same day: Kim 2024b matched; Bischof
  2023 did not on a single number (the fetch returned the α-glucose row for one breed and the β-glucose row for the
  other under "glucose", and IMP, inosine, leucine and methionine values that appear nowhere in Table 1, while the
  dossier's own flag named a milder fault, repeated columns); and the review used for two beef sugar corners cites a
  paper whose title says CHICKEN (Aliani 2013), which the fetch never showed because the reference list did not come
  through. Two of three box corners that rested on unchecked reads were wrong in kind, not in digit. Rule: a range
  enters a box, a dossier's section 2, or a note as `sourced` only from the measuring paper's table read by eye from
  the PDF on disk; an automated read is `secondary` at best and its dossier says so in the first line. A value taken
  from a review is traced to the primary's TITLE (species, matrix) before it is used even as a corner. And when the
  PDF arrives, check its first author against the dossier stem: `hwang2026` is Dashmaa et al.

- **Quote the value the engine runs, not the prior it was searched around; and search for a paper by its stem.**
  On 2026-10-09 the thiamine route's barrier was described to the owner, and to four reading agents, as "about
  100 kJ/mol": that is `FORMATION_EA_PRIOR_CENTRE["thiol_assembly"]`, the centre of a search band. The engine runs
  the frozen lumped barrier, 64.08 kJ/mol, because B10's split was re-merged; three dossiers had to be corrected.
  The same afternoon `schieberle2000` was reported as "not in the model" from a grep for author-name patterns,
  when it has a dossier, a registry entry and a role in B16. Rule: before stating a parameter, read it from the
  frozen fit the engine loads (`core_prediction_uncertainty.json` lists every one with its source); before
  saying a paper is absent, grep its dossier stem (`<author><year>`) in `data/keys/papers.yml` and
  `data/lit/extraction_dossiers/`.

- **The project's purpose frames every brief, whatever branch is checked out.** On 2026-10-09 a decision brief
  on "which experiments are worth running" was written from the `cultivated-tissue` branch and framed the whole
  model's value around cultivated meat. Maillard exists to help build meat aroma into **plant-based** products
  (README, first paragraph); cultivated tissue is one side study. Rule: before writing anything that argues for
  the model's value, re-read the README's opening paragraph and state the purpose in the first lines; a side
  branch's question is an example, never the frame.

- **Deliverables live in the repository, and outward emails point to it.** On 2026-10-09 an approved
  "write up the map" was put in a private claude.ai doc, outside the repo, then linked from a draft email the
  owner had not seen it in. It duplicated `docs/guides/EXPERIMENTS_BRIEF.md`. Rule: a write-up about this
  model goes into the repo beside what it summarises (extend an existing doc before creating one), passes
  through a PR like everything else, and an email to a third party is short and links only to published
  repository files.

- **A headline number is only as good as its definition; read the test before repeating it.** On
  2026-10-27 the owner was told, for a meeting with WUR, that "63 unsourced records are used when the model
  runs", repeated from README. The test behind it counted any record under `data/lit/` as "reaching runtime";
  a trace showed the front door opens none of those files. Rule: before quoting a repository statistic to
  anyone outside, open the test or generator that pins it and check that its definition measures what the
  sentence claims; where it does not, measure (trace, recount) before speaking.

- **A verification claim names its scope, or it reads as total.** The owner believed the literature was sane
  because the record said "every defective anchor was repaired or explicitly labelled" and "the citation gate
  blocks a confabulated DOI". Both were true of a narrower scope than they sounded: labelled records stayed in
  the tree, and the gate never sees a live DOI that points at the wrong paper, nor any value. Rule: a sentence
  that says something is verified states what was checked (which files, against what) and what was not; an
  unsourced record is deleted, not labelled and kept.

## Textbook constants carry citations that were never checked (2026-10-27, B51)

- **What happened.** The cysteine thiol and ammonium ionisation enthalpies (26 and 43 kJ/mol) had been
  cited to Goldberg, Kishore & Lennen 2002 for two months; GKL's Table 7.23 prints 36.1 and 34.1. The H2S
  pK1 was extrapolated by constant-dH van 't Hoff to 5.94 at 145 °C; Suleimenov & Seward measured 6.49.
  No ionisation heat capacity was carried anywhere, though the module itself tabulates pKw because
  constant-dH van 't Hoff is 0.49 units wrong for water at 145 °C. Every sulfur fit since B2.1 was
  calibrated on these speciation factors.
- **Rule.** "Textbook" and "immovable" are claims, not exemptions: an equilibrium constant gets the
  same check against its cited table as a rate constant gets against its paper. Any temperature
  extrapolation over more than ~50 K carries the source's dCp, or a measured high-temperature
  correlation, or says why not.

## A misfit was blamed on a mechanism before the paths and the numerics were checked (2026-10-28, B60)

- **What happened.** After B51, the Hofmann ribose/xylose pH 5 rows sat ~100x low and the record blamed
  "the thiolate-acting thiol sink". A flux attribution two days later showed that sink takes 0–4 % of MFT's
  loss in those pots. The real causes were (1) a pH trajectory solved on 41 uniform nodes that missed a
  transient lasting seconds, and a pH fixed point that oscillated ~1 pH unit and was stopped after 2–4
  sweeps without checking its residual (MFT moved 39x when both were fixed); and (2) a trunk that destroys
  cysteine in ~1 s where three measured pots say hours, because endpoint fit rows constrain branch ratios
  but not timescales. The residual `ph_fixed_point_residual` had been written to every run's metadata and
  read by nobody.
- **Rule.** Before attributing a misfit to a mechanism: (a) show the numerics converged at the speed of
  the fastest transient (refine the grid, iterate the fixed point to tolerance); (b) run
  `build_sulfur_flux_attribution.py` and confirm the blamed channel actually carries the flux; (c) compare
  the intermediates' lifetimes with any time course in the corpus. A diagnostic a code writes must be
  asserted somewhere, or it is decoration.
- **Confirmed by the result.** On B60's converged numerics those four rows went from 37–118x off to
  1.5–3.6x with no change to the sink, and the B51 model-card entry was withdrawn. Two more things showed
  up. First, B51's reported fit cost (20.4) was itself flattered by the unconverged numerics: the same
  vector costs 33.2 integrated correctly. A fit's cost is only as good as the integrator inside its
  objective. Second, a legacy switch must cover EVERY numerical change: the closed flux budget also moved
  a frozen wave's share rows (4.6e-5 relative). Run the frozen-wave tests before declaring a numerics
  correction backward-compatible.

## A pseudo-inverse reported "identified" for directions the data never saw (2026-10-28, B60)

- **What happened.** The Laplace covariance drops JᵀJ's near-null eigen-directions (to keep it PSD) and
  then reads σ off the diagonal. A dropped direction contributes ZERO variance, so a constant with an
  exactly zero Jacobian column (`k_glc_ha`: no fit row contains glucose) came out σ = 1e-9 and
  "IDENTIFIED". `k_thiol_decay` and `k_dimer_decay` did the same, in B51 and in B60. The envelope sampled
  "identified" constants from the covariance, so it pinned all three at the optimum.
- **Rule.** Truncating a pseudo-inverse is a statement that those directions have INFINITE variance, and
  the output must say so: record each coordinate's weight in the dropped subspace and mark it
  unidentified. More generally, a sigma far below the step used to differentiate (1e-9 against a 0.05
  step) is a red flag, not a precision.


## The flat-draw box and the flux-less barrier (2026-10-28, B60 envelope and the B62 design)

- **What happened (1).** After the Laplace fix, three unidentified coordinates were drawn flat across the
  optimiser's box [-10, 0.5]. The fit rows reject the top 1.5-2 decades of that box (slice cost +714 at
  the top of `k_thiol_decay`), and those draws put 6-9 dex intervals on rows the coordinates were fitted
  on. Coverage rose from 20/50 to 38/66, mostly bought with width.
- **Rule (1).** An optimiser bound is not a prior. Draw an unidentified coordinate over the range the data
  ADMIT: the profile's walk to each bound, cut at the same Δχ² the Laplace uses. Before reporting a
  coverage gain, look at the widths that bought it.
- **What happened (2).** I told Pablo that three sources put furfural FORMATION at ~128 kJ/mol. Two of
  them (Garrett 1969, Danon 2014) measure amine-free, acid-catalysed pentose -> furfural, a route that
  carries ~0 flux in every pot we score. The flux attribution showed that `r_pent_tdp` (amine-catalysed)
  makes 95-100 % of the 3-deoxyosone. Pablo approved the wave on that claim; I caught it afterwards and
  corrected it.
- **Rule (2).** Before mapping a literature barrier or rate onto an engine step, check that the step
  carries the flux in the pots it will be scored on (`sulfur_flux_attribution.md`). Also check that the
  measured system's catalysis regime is the operative one. A barrier for an idle route says nothing about
  the busy one.
- **And (3), compute.** I ran 12+ processes on 10 cores (load ~170). Every job ran at ~50 % and the user
  had no CPU left. Cap my parallel work at about half the cores, nice it, and stop superseded probes.
  `scipy.optimize.least_squares`' `max_nfev` does not count finite-difference Jacobian evaluations, so
  budget a k-parameter refit at ~(k+1) x max_nfev integrations.

## 2026-10-29 — an independent review failed B59a and corrected B64 (lessons from my own errors)
- **Check a new reaction's REDOX, not only C/N/S.** `r_dpo_ddp_h2s` balanced carbon, nitrogen and sulfur but emitted
  the 2e-reduced DDP together with H2S, double-counting the enaminol's reducing equivalents. For any step that moves
  sulfur between oxidation states, write the half-reactions and check H and electrons before running anything.
- **A yield is not a bound on a branching share** unless every unit of product is known to come from the branch.
  Ask what else makes the product (a blank?) and whether the yield scales with a co-reactant: if it does, a fixed
  share is the wrong functional form.
- **Read the source's numbers, not its prose.** Schieberle 2000's "a factor of 13" is 7.9x on its own tables. A
  discriminator taken from a sentence must be recomputed from the tables, with the discrepancy as its uncertainty.
- **Amendments after results are selection.** Refitting a variant because a sweep lowered its cost, or carrying "the
  best-fitting" baselines forward, is choosing on scores. Run every pre-registered alternative, or declare the
  selection as post hoc in the record.
- **"V1-V5 PASS" must name what it covers.** A verification run on two readings does not certify a third.
- **Pre-register the mechanical ship rule with the audit**, and make the engine check (item 4) a non-tautological
  reproduction of the source pot.

## 2026-10-29 (later) — 15 refits launched before the cheap question was asked (owner: "get clarity before expensive runs")
- **Before any refit, run the change at the incumbent's constants with NO refit** and score it row by row and on
  every hold-out. Five minutes on one core showed what hours of refits would have buried: which rows the change
  breaks, why (a flux budget), and that the change already moves the hold-outs the right way. A refit is for
  rebalancing a change you understand, never for discovering what it does.
- **"Refit only the coordinates the data identify" (WAVES.md) applies to every kind-2 refit**, and the identified set
  is a property of the rows AND the structure: after a structural change it must be re-derived, not inherited.
- **List every direct measurement of each NEW species before writing audit item 5.** B59b's audit missed Zhai 2019
  Fig. 3A (the route's own intermediates at four temperatures), which contradicts the transferred formation barrier.
  Grep the dossiers for the species' names (TTCA, Cys-Amadori, ARP), not only for the step's name.
- **When a fit moves a constant a long way, ask which measurement in the objective sees that constant.** If none
  does, the move is compensation, and the fix is data that see it, not a bound.
- **Check the record for the same idea before proposing it.** "Add Schieberle's 100 °C series as fit rows" was B16
  (2026-09-07), and it failed for reasons that changed the plan.
- **Copy the codebase's optimiser settings, and restart-check every fit.** The route calibration's first fits left
  `least_squares` at its default finite-difference step (~1e-8, below the ODE noise at rtol 1e-7) and reported a
  convergence that a restart from the other arm's optimum disproved (8.27 → 4.45). Every fit here uses
  `diff_step = 3e-2`. "Converged" means two starts agree, not that the optimiser's status flag says so.
- **Background agents must save as they go.** Two agents (a digitisation, a literature search) were lost when the
  session ended mid-task with nothing written. Every long agent gets an output file it appends to after each unit
  of work.
- **Question inherited rules on their merits** (owner, 2026-10-29). "Rates, fed yields and ratios fit; levels
  validate" encoded a real worry (calibration leaking into constants) with the wrong remedy (excluding data by
  type); the right tool was a measurement model. A written rule is evidence of someone's reasoning, not proof of it.
- **scipy's `diff_step` is RELATIVE** (3 % of the coordinate's value). At a coordinate near zero the step collapses
  to machine epsilon and the Jacobian column is the ODE integrator's noise; across units it is inconsistent (3 kJ/mol
  on an Ea of 100, 0.018 dex on a log10 k of -0.6). Every fit copying B8's setting shares this for coordinates near
  zero (B60's acid yield sits at 7.5e-5). `src/kinetic_core/fitting.py` uses an absolute step per coordinate kind.
  Found by writing a test with a known answer on a rough objective, which is the point of shared, tested tooling.
- **A parameter shift measured in the full fit's sigmas can hide a tension**, because the tension inflates the full
  fit's chi-square and so its sigmas. Leave-one-group-out also needs the predictive check: the held-out group's
  residuals under the model fitted without it, against that fit's own scatter (an F test).
- **A dossier's summary label is not the source.** Zhai 2021's kinetic rows were run at pH 5.5 but labelled pH 7 in
  the dossier heading; B16 and stage 1 both inherited it. Before a row enters a fit, its conditions are checked
  against the paper's own text (the measurement index should carry the quoted condition, not the summary).
- **Generators set engine switches AT IMPORT, and a test that imports one leaks them to every later test**
  (2026-10-29). `generate_kinetic_core_b64_fit` sets `sulfur.PH_COUPLING = "inline"` at import; the route-calibration
  test imported it, and B2.4's frozen-objective test then ran under the inline coupling: over an hour of a unit tier
  that takes minutes, and a check of the wrong configuration had it finished. The B2.3 row splicing leaked the same way
  (cost 22.56 -> 23.44). Guards now: `tests/conftest.py` resets every `src.kinetic_core` UPPERCASE global to its
  pre-collection value before each test, and `tests/support.wave_generator` restores every generator module's
  namespace on exit. Rules: a test imports a wave generator only through `wave_generator`; a new generator sets its
  switches inside a context (as the route calibration's `installed()` does), not at import. Found by pairing each
  earlier test file with the slow one (a bisection that took one minute); "probably slow" was the wrong first guess.
- **A test that is slow in the suite and fast alone is test pollution, not slowness.** Time it alone first.
- **A phase that regenerates headline artifacts is not closed on the unit tier** (2026-10-29). Phase 0 was declared
  closed on 1091 unit tests and the gates; CI's scientific tier then failed five headline guards (pinned 17/60, 33/59,
  21/38, 17 bundles, "no Hofmann pot in the fit") that the regeneration had legitimately moved. The guards were right
  to fail, but the prose in README had also gone stale, and a merge was delayed. Rule: after regenerating any
  scorecard, envelope or baseline, run `tests/scientific/test_core_headline_guards.py` (about 10 min) before calling it
  done, and re-pin with the reason in the same commit.
- **Local checks must be case-exact** (2026-10-29). macOS folds case and CI's Linux does not: `xin2026b` and
  `zhang2024_extraction.md` passed every local gate and failed CI. Gates match names against directory listings or
  git's tracked paths, never `Path.exists()`.
- **Never brief an agent with `python -I -c 'from src ...'`** (2026-10-29). `-I` drops the working directory from
  `sys.path`, so every such check fails with `No module named src` before reading anything; three adjudicators hit it.
  `-I` is for scripts that read untrusted downloads, not for repo imports.
- **Verification is scoped to the table that reads the source** (2026-10-29). Keyed by ledger id alone, one
  sub-model's verified reading cleared another sub-model's pending rows under the same id. A permission check keyed
  coarser than the thing it protects leaks across users of the key.
- **Read a function's sign convention from its code, not from memory** (2026-10-29). The sub-model 1 probe report
  and prereg said "residual = measured minus predicted"; `fit_rows.residual` returns predicted minus measured. The
  findings were stated from the numbers and stood, but a reader acting on the label would have read every bias
  backwards. Labels that describe a computation are written next to the line that computes it.
- **Pin BLAS / numba threads to one per worker** when forking `MAILLARD_WORKERS` processes beside another parallel
  job (2026-10-29): the probe batch next to the unit tier drove the load average to ~250 and the probes to ~7 % CPU
  each. `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=NUMBA_NUM_THREADS=1`.
- **A wait loop must not match itself** (2026-10-30). `until ! pgrep -f "pytest tests/unit"; do ...` never ends: the
  loop's own shell command line contains the pattern. Three such loops ran ~8 h. Wait on a PID (`wait`, `kill -0
  $PID`), on a log line, or on a pattern the waiter cannot contain (`pgrep -f "[p]ytest tests/unit"`).
- **Before a push, run every CI gate, not only the one that failed last** (2026-10-30). A fix for the freshness gate
  was pushed after the unit tier and that gate passed locally; the data-integrity tier then failed on a register entry
  added the same morning (CF-068: a kind outside the enum, a document name where ledger ids belong). Run
  `./scripts/docker_maillard.sh gates` (or every step of `.github/workflows/ci.yml`'s data-integrity job) before a
  push that touches `docs/reference/*_register.yml`, `data/`, or the generators.
- **A wait loop needs a way out** (2026-10-31). `until grep -q "<line>" log; do sleep 30; done` waited ~12 h for a
  freshness-gate line that the producing run (a detached `&` job whose chain stopped after pytest) never printed. Wait
  on the producing process, not on a string (`run_in_background` and its completion notice), or give the loop a
  deadline and print what it last saw. Never start a detached `&` job: its output and exit are invisible.

**Generated READMEs count TRACKED files: regenerate them after `git add`** (2026-10-31). `build_data_readme.py` and
`build_results_readme.py` list what git tracks, so running them before staging new files leaves both stale in the
commit (unit tests `test_*_readme_is_current` fail). Order: stage everything, regenerate the READMEs, stage them,
commit. And a verification chain must not stop at the first failing step: the env wrapper runs `set -eo pipefail`,
so `pytest ... | tail` kills the rest of the chain -- redirect each step to its own log instead.

**Even text munging goes through the env wrapper** (2026-10-31, again). I parsed `gh pr view` JSON and edited the PR
body with system `python3` -- the rule has no exception for "not repo code". Use
`./scripts/docker_maillard.sh run "python ..."`, or `jq` / `gh -q` for JSON.

**Backticks inside `docker_maillard.sh run "..."` are command substitution** (2026-10-31). A Markdown `code` span in
a heredoc passed in double quotes ran the script name as a command and silently dropped it from the text. Write such
edits with the Edit tool, or single-quote the wrapper argument.

**I started a detached `&` job again** (2026-10-31 05:08, the B7a arms queue) -- against the wait-loop lesson above.
The harness cannot notify on it. Remedy used: a background waiter on the job's PID with a deadline. Rule: a queue is
always launched as the background command itself (run_in_background), never with a trailing `&`.

**The paper registry records citation LINE numbers** (2026-10-31): any edit above a cited line in a register or
ledger moves them, and `build_paper_registry.py --check` (unit test + freshness gate) fails. Regenerate it with the
READMEs after staging, as part of every commit that touches docs/reference or data/lit.

**Describe a paper only after its dossier exists** (2026-10-31). I gave the owner one-line descriptions of the Tai &
Ho papers from a literature agent's summary ("cysteine vs cystine"); the PDFs compare cysteine with its sulfinic /
cysteic acids and GSH with its sulfonic acid. The DOIs were right, the contents were not. When listing papers to
download, say "reported to" and keep the claim to what a verified abstract states.

**An agent's factual claim about a paper is a lead, not a fact** (2026-10-31, again). "Zhai 2021's pots are 0.2 M
phosphate" came from a dossier agent's report; I wrote it into a pre-registration; the dossier and the fit rows say
unbuffered water, and the cysteine review caught it. Before a claim about conditions enters a prereg, grep the
dossier / fit table for it (one command).

**An engine edit stales every signature record, even when bit-identical** (2026-10-31, twice). The scorers record
the sha256 of sulfur.py / parameters_sulfur.py, so B7a's and B7a-2's switches (defaults bit-identical) left ~30
signature and chance records stale for the freshness gate. Remedy now tracked:
`scripts/generators/regenerate_signature_records.py` re-scores every stale signature / chance record (and lists the
other generators to re-run: data_wishlist, sulfur_flux_attribution). Run it after any src/kinetic_core edit, before
the full verification.

## 2026-10-31 — Side problems go to subagents, not the main thread
- **Pattern:** the feasibility tester's validation (v1–v5, power test, oracle check) ate the main thread for hours while the plan's active step (commit, verify, push, split) waited. Owner: "leave rabbit hole problems to subagents, so you do not get off track and cover more ground."
- **Rule:** when a sub-question appears mid-step, brief a background subagent (files it may touch, what not to touch, compute cap, no commits) and continue the active step; fold the result in when it reports.
- 2026-10-31 (again): a text-patch script was run with system python3 -I. Even pure text patches go through the wrapper: ./scripts/docker_maillard.sh run "python <patch.py>".

**An overnight run on a laptop must hold the machine awake, and a per-pot wall-clock limit is not sleep-safe**
(2026-11-01, the thiol leave-one-lab-out). The laptop hibernated on 1 % battery from 03:29 to 09:33 in the middle of
a fold's fit; `run_pot_limited` uses SIGALRM, so pots in flight would time out on wake and enter the optimiser as
failed residuals. Remedy used: stop, resume from the per-segment checkpoints, relaunch under `caffeinate -i` on AC.
Rule: launch long jobs as `caffeinate -i nice -n 10 ...` and check `pmset -g batt` first; checkpoint every fit segment
so a stop costs one segment. And a summary averaged over folds must skip a fold with no scored row (Davidek's species
have no null answer): counting it as a difference of 0 pulls the mean to zero and shrinks its se.
## 2026-10-31 — A rival's band derived from the test it is scored on, and data left unscored
- **Pattern:** an investigation agent admitted "slow cysteine exchange". Its rate ceiling had been derived from the same xylose contrast it was then scored on, and two concentrated-pot datasets that contradict it (same lab and another lab) were left unscored under a reason that did not apply. It also counted 50 C storage data as a second lab for a processing loss law. The independent red team caught all three.
- **Rule:** a rival's bands come from data OTHER than its scored tests; every on-disk dataset that bears on the question is scored or its exclusion justified in the prereg; a "second lab" must measure the same channel. Every structural verdict gets a red team before it moves the engine.

## 2026-10-31 — An edit to the conflicts register or the measurement ledger stales sulfur_reaction_evidence
- **Pattern:** twice in one day a branch failed the freshness gate after only register / ledger edits, because `results/validation/sulfur_reaction_evidence.json` hashes `docs/reference/conflicts_register.yml` and `data/lit/measurement_ledger/trunk.yml`.
- **Rule:** after any edit to either file, run `python scripts/generators/build_sulfur_reaction_evidence.py` before staging (it takes seconds), and `regenerate_signature_records.py --dry-run` lists anything else stale.

## 2026-11-01 — Dossier agents skip the DOI line; check it after every ingestion batch
- **Pattern:** 9 of 36 new dossiers had no DOI line, although the brief asked for one; the paper registry's "dossiers without a printed DOI" count (30 -> 41) was the only signal.
- **Rule:** after an ingestion batch, compare that count before and after, and grep each new dossier for a `DOI:` line in its header (before the first `##` / `###`); add missing ones only after a Crossref check.

## 2026-11-01 — `git add -A` swept an agent's half-written files into a commit
- **Pattern:** a docs-only commit (the queue made private) used `git add -A` while a leave-one-lab-out agent was writing in the same worktree; its unfinished prereg, generator and test were committed and pushed, and the next full suite failed on the undescribed prereg.
- **Rule:** while any agent may write in a worktree, stage by explicit path (`git add <files>`), never `-A` / `.`; check `git status --short` against the intended file list before committing.
- **Recurred 2026-11-03:** `git add data/lit` (a directory, not `-A`) swept two dossiers another agent was finishing into an ingest commit whose message did not name them (amended before push). A directory is not an explicit path: name the files.

## 2026-11-01 — A red team's number is a lead too
- **Pattern:** the engine-error red team ranked E4 (the pyrazine pH law) as SUPPORTED on "methylpyrazine pH 7 / 5 = 170 vs Leahy's 14". 14 was pyrazine's ratio; methylpyrazine's is 215, so the shipped engine already matched. The premortem inherited the error, and only the pre-registration agent, reading the table for its test, caught it.
- **Rule:** before a red team's or reviewer's number orders work or a change, the coordinator checks it against the dossier table (one grep). Reviews reduce error; they do not certify numbers.

## 2026-11-01 — Check the alias table before calling a PDF "wrong"
- **Pattern:** an agent reported Zhou2023.pdf as "the wrong paper", and I passed it on to the owner as "please replace". data/keys/pdf_dossier_aliases.yml had recorded since 2026-08 that Zhou2023.pdf is Zhu et al. (kafirin) and that the wanted paper is on disk as Zhou2023b.pdf.
- **Rule:** before telling the owner a file is missing or wrong, grep pdf_dossier_aliases.yml and the dossier's "Source on disk" line, and check +/- 1 year and the b / c suffixes.

## 2026-11-01 — A frozen artifact must be checked against its inputs at freeze time
- **Pattern:** the benchmark's frozen v1 split was checked by regenerating it from the LIVE cited documents. A legitimate later edit (the owner's decision note appended to the E4 prereg) made the check fail although the split was unchanged.
- **Rule:** a frozen record checks against `git show <freeze commit>:<input>`; live drift is reported separately and never silently fails or silently passes.

## 2026-11-01 — Check "held-out" against the frozen split files; cumulative yields are not pool sizes
- **Pattern:** the buffer / pyrazine resolution called Larralde 1995 and Lee 2024 "TEST" for E5 / E13. Neither sits in any frozen split: both are keyed held-out tables, candidate TEST for a later version (contrast-benchmark README sec. 10). In the same record a pool-invariance test read sqrt(cumulative yield) as the instantaneous aminoketone pool and "rejected" the lower-order extra MPZ route at |z| 3-6. A cumulative yield grows even at a constant pool; with rate-based proxies the pools barely moved and the test had no power. The red team caught both, after the claims had reached the draft register text.
- **Rule:** before writing "TEST", "DEV" or "held-out" for a source, grep `docs/reference/contrast_benchmark/split.yml` and `split_v2_per_change.yml` for its cell or lab and say which version; a keyed table outside them is "keyed held-out (not in a frozen split)". Before using a time course as a pool, check whether the printed quantity is cumulative (yield, amount formed) or standing (concentration at t); a pool test needs the standing quantity or a rate-based proxy, and points in net loss carry no pool information.

## 2026-11-01 — A worker's switches are not the parent's: post-processing in the parent read a different state
- **Pattern:** the sub-model fits installed each pot's switches in `run_task`, then predicted the rows in the PARENT. Serially the parent held the last task's state; with `MAILLARD_WORKERS > 1` it held its own, so an entry arm's pH rows were predicted without the trunk solutes the integration had in its charge balance (Cao 2017's hold-out pH, up to 0.14 pH; found by A3 phase 2's serial / fork / persistent comparison). The fitted rows were not pH rows, but every hold-out pH a forked run reported was computed in the wrong state.
- **Rule:** a result is read in the configuration it was computed in: build each task's `Config` in the caller (`parallel.configured`) and read its result under it (`parallel.run_under`). A check of a parallel path compares serial, forked and persistent runs, not one of them against a record.

## 2026-11-01 — A limit argument in a pre-registration must be checked on the rows' own regime
- **Pattern:** CF-077's prereg predicted the B18 fit unaffected by the cross-pair factor because the fit-time Leahy pots ran in the minor-pool limit of the FAST-condensation partition. At the fit commit the Leahy rows moved up to +0.51 dex (cost 2.64 -> 15.4): in the dilute acid arm the aminoketone lifetime (1 / sqrt(2 k F) ~ 400 min) exceeds the window, so the condensation is not fast there. "Declared fast" was true of Zhou's 20 mM pots and silently assumed for every pot.
- **Rule:** before a prediction leans on a limit (fast step, steady state, minor pool), compute the dimensionless number that licenses it (here k_cond x pool x window) for each scored row, and state it in the prereg. Running the check at the fit's own commit (a detached worktree) is cheap and decided it.

## 2026-11-04 — Never pause a fit that runs under a wall-clock limit; size concurrency to free memory at launch
- **Pattern:** H0's pair fits thrashed on a 16 GB machine (a VM held 10 GB). I first raised concurrency (2 -> 5 fits), then paused two with SIGSTOP to relieve memory. Each pot integration runs under a 120 s wall-clock alarm (`run_pot_limited`), so paused or memory-starved integrations "timed out" and entered the objective as failed pot runs: three of ten records reported 55 / 117 / 50 failures where every other fit had 0. They had to be re-run blind (amendment A2).
- **Rule:** before a long run, check free memory (`vm.swapusage`, `top -o mem`), not only cores; set concurrency so the fits fit in RAM. Never SIGSTOP a running fit; reduce concurrency only at launch. After any run, read `failed_pot_runs` before using a record: a non-zero count means the record depends on machine state.
## 2026-11-03 — The plan may be newer on an unmerged branch; quote the one plan, and delete superseded rules where they live
- **Pattern:** asked for the goal, I quoted the 2026-10-28 "usefulness bar" (>= 40 % within x3, >= 65 % directional) from main's `tasks/data_restructure_plan.md` and README model card. The owner had removed round-number thresholds in the 11-02 rewrite of `STATUS_AND_PLAN.md` (sec. 7), which sat on a cloud branch with no PR, so main still printed the bar, and the model card still computed MET / NOT MET against it.
- **Rule:** before stating the plan or how success is judged, run `git log --all --since=<last week> -- docs/guides/STATUS_AND_PLAN.md` and read any commit not on main. A rule the owner withdraws is deleted where it lives (code, generated card, backlog, register), not only marked, in the same PR that lands the plan.

## 2026-11-04 — Read a scored item's sign convention before stating which way a measurement moves
- **Pattern:** reading the W1.1 score, I told the owner that Meynier's norfuraneol "rises 0.4 - 1.2 dex" with pH and built a plan argument on it. The items are `m = log10(a / b)` with `a` the LOWER pH, so a positive m means norfuraneol FALLS as pH rises; the fit-row header says so ("fall with pH in every table"). The direction decides which rivals are plausible (a rising norfuraneol favours a split moving to 2,3-enolisation; a falling one equally fits base-catalysed decay).
- **Rule:** before saying which way a measured contrast goes, read the scorer's definition of the item's sign (`p_of`, `m`) and check it against the source table or its dossier note. Say "X at pH 6.5 is N dex below pH 4.5", not "rises" or "falls" alone.

## 2026-11-04 — Check which steps a switch touches, and what earlier waves found, before naming it as a plan step
- **Pattern:** I added M1b to the plan as "the existing `ph.trunk_terms` switch", to give Meynier's ribose pots a pH response. The switch slopes only the three hexose Amadori-decay steps (`trunk_conditions.PH_STEPS`), so it is inert in every pentose pot, the prompting ones included. The same plan's M1 restated an oxidant structure (B27, gated out before its fit on 2026-09-11), a headspace O2 reservoir (B11, DO NOT SHIP) and a fitted disulfide release (B17b, DO NOT SHIP) without their outcomes. The case drafts caught both, before any build.
- **Rule:** before a plan row names a switch or structure, read the code path it changes (which parameter keys, which pots reach them) and grep `tasks/data_restructure_plan.md`, `docs/reference/experiment_register.yml` and `results/validation/` for earlier waves on it; the row cites their outcomes.

## 2026-11-05 — A measurement model is part of the design: check floors, offsets and densities against the data's own form
- **Pattern:** M1's frozen prereg passed two independent design reviews, yet its verdict ("not separated") rested on a linearisation artefact. Three measurement-model choices did it: (1) detection floors set at measured values (the 30- and 50-min data points) flattened one rival to a zero-width line at a datum; (2) H0's between-laboratory tau (~1 dex, from absolute thiol levels) was applied to a C/C0-normalised series where a calibration offset cancels, hiding the incumbent's main failure; (3) the floor sensitivities used a linearised density where the rival's prediction is strongly non-linear (exact lppd -746 vs linearised -175), and the MC check ran only on the primary setting.
- **Rule:** before freezing, for each decisive series state its value basis (absolute, ratio, C/C0) and apply a between-lab offset only where the basis keeps one; set floors from the assay (LOD / LOQ, calibration range), never from a measured point; and run the exact (or adequately sampled MC) density on EVERY form that can flip a sign, not only the primary. A reviewer should check the deciding number's density form, not only the rule's logic.

## 2026-11-05 — A brief must not ask for a band that contains the scored datum; screen structures before designing a fit around them
- **Pattern:** briefing the M1 + M2 prereg, I asked for the reservoir's release-rate band to be "wide enough to include" Zhang 2025b's measured timescale. Zhang's pot was a decisive unit, so the band would have been set by the test it is scored on (the 2026-10-31 lesson). The drafter caught it. The same draft built a whole design around a pH-dependent exchange that a 1 h plug-in screen then showed cannot reproduce Zhai 2021's pH series or Zhang's mass balance; the design's supply axis emptied.
- **Rule:** in a brief, name where a band comes from (a source that is not scored), never which value it must reach. Before drafting a fit design around a structure, run the mechanism-signature screen (plug-in runs on the signatures the evidence names) and design only around the structures that pass; record which pots the screen has seen, since they cannot be decisive afterwards.

## 2026-11-05 — Compare a pH on the scale it was measured on (cooled vs in situ)
- **Pattern:** a screen reported that the engine's pH in Cao 2017's pot "starts at 3.86 where Cao measures 4.67" and barely drifts. I carried it into the plan, an agent brief and the owner's summary as an upstream fault that justified a re-sequence. The 3.86 was the engine's in-situ pH at 90 C; Cao read pH after cooling. On the cooled scale (`observations.ph_cooled`) the engine starts at 5.50 like Cao and misses by ~0.2, all of it one known allocation (CF-084). A ~1-unit in-situ / cooled split is expected chemistry here, not a fault.
- **Rule:** every pH number in a record, plan line or brief names its scale (cooled / label 25 C, or in situ at T). Compare the engine to a measurement only on the measurement's scale, using the observable the probe scores. Before acting on a surprising residual, reproduce it on the scored observable.

## 2026-11-05 — Never wait on a job with `until ! pgrep -f "<pattern>"`
- **Pattern:** subagents waited for the gates with `until ! pgrep -f "docker_maillard.sh gates" >/dev/null; do sleep 10; done`. `pgrep -f` matches full command lines, and the waiting shell's own line contains the pattern, so the loop could never end. Six such waiters (one on a log line that a finished job never wrote) hung for ~11 h, and the H0-application agent sat "waiting on its own background work" the whole time.
- **Rule:** run long commands in the background and wait for the completion notification; or wait on a known PID (`while kill -0 $PID`) or a file the job is sure to write, with a timeout. If a pattern match is unavoidable, use the bracket trick (`pgrep -f "[d]ocker_maillard.sh gates"`). Every subagent brief that may wait on a job says so.

## 2026-11-05 — Blind keying may run early, but its commit lands on top of the frozen design
- **Pattern:** to save time I had two blind readers digitise Hofmann 2000b Figs. 2–3 while M3's prereg was still a draft, and told the reconciling agent to commit and push the keyed rows. The M3 review then required that the frozen design be an ANCESTOR of the keyed-data commit (otherwise git cannot show the design was fixed before the data existed in the repo). I stopped the push in time.
- **Rule:** digitising early is fine if the readers are blind and the output stays outside the repo (scratchpad, hash recorded). The keyed file is committed only on a branch descended from the design's freeze commit, and the scorer checks that ancestry plus the file hash. Design authors are told not to open the scratch readings.

## 2026-11-05 — `bootstrap` upgrades every unpinned package; a unit-tuned proposal at a kinked peak is too narrow
- **Pattern (env):** adding arviz to `environment.yml` and running `./scripts/docker_maillard.sh bootstrap` (`conda env update`) also upgraded numpy 2.5.3 -> 2.5.4, numba, scipy's build, networkx, pytest, pip, setuptools, pyyaml's build and coverage: the numbers in `results/validation/` come from the old set. Reverted by reinstalling the exact previous builds (`conda list --export` before and after, diffed).
- **Rule:** before `bootstrap`, export the env (`conda list --export > before.txt`); install the new package first with `conda install -n maillard <pkg>=<ver>` and check the plan lists only NEW packages; after `bootstrap`, diff the export and restore any upgraded build. Only add, never upgrade. Since M1.1 the packages that moved are pinned in `environment.yml` (versions; scipy and pyyaml also their build numbers), and a bootstrap leaves the env unchanged (export diff empty); unpinned ones can still move.
- **Pattern (densities):** M1's R0 Cys predictive under the alternative floor (M1.1's regression, from a 25-point engine scan) has a kink at its peak, where the legacy-floored 50-min prediction crosses 0.22 mM. The local inverse Hessian there gives sd 0.0016 against the integrand's ~0.03; a Student-t proposal built on it has k-hat 0.92 and ESS 45 (refused), a proposal as wide as the integrand has k-hat < 0 and agrees with a fine trapezoid (-747.4 vs -747.5).
- **Rule:** scale a unit-tuned proposal from the integrand's width over a step comparable to the posterior (or a profile), not from a curvature at a point; any floor, bound or switch can put a kink at the mode. PSIS's k-hat is the check that catches it; never skip it.

## 2026-11-05 — A reference generator moves the engine only through Config or the engine's own arguments
- **Pattern:** the JAX lipid module's reference generator patched `parameters_lipid.K_LOOH_DECOMP_ANCHOR` (restored in a `finally`) to get engine finite differences for a rate-multiplier coordinate. `test_engine_config.py::test_every_switch_a_script_or_test_sets_has_a_config_key` failed the unit tier on it: any script assignment of an engine UPPERCASE global counts as a switch without a key, restored or not.
- **Rule:** an engine reference moves a coordinate only through `config.applied(Config...replace(...))` or a function argument (`integrate_lipid(q10=, donor=)`). A coordinate the engine has no knob for is checked against the port's own central differences, and the reference says so; never add a patch to get an engine number.
- **Also:** do not carry "the faster solver" from one lane to another. Kvaerno4 beat Kvaerno5 on the sulfur bake-off pot, but took ~4x more steps on the lipid pots, where pools fall far below atol. The solver is a per-network choice, measured on that network's pots.

## 2026-11-05 — Generated READMEs go stale when main moves under an open PR
- **Pattern:** two PRs (#60, #63) failed their unit tests on `data/README.md` / `results/README.md` being stale. CI tests the PR merged into main; other PRs had merged new dossiers and results files meanwhile, and `data/README.md` counts committed files, so regenerating it before the commit that adds them leaves it one step behind.
- **Rule:** last step before pushing any branch that touches `data/`, `results/` or the registries: merge `origin/main`, regenerate `build_paper_registry.py`, `build_data_readme.py` and `build_results_readme.py`, commit, then run each with `--check` (they must print "current"). Re-do it whenever a sibling PR merges first.

## 2026-11-05 — Never write an author or file key that was not read from the record
- **Pattern:** a side session named Lopez et al. 2019 "Zhou 2019" (and queued `zhou2019_mushroom.pdf`) from a title-only OpenAlex row; the PDF's first page shows Lopez, Kerley, Jenkinson, Luckett & Munafo.
- **Rule:** author names and paper keys come from the authorship record (OpenAlex / Crossref authorships) or the PDF's first page, never from a title or memory. `openalex.py verify` (W1.4) checks every `papers.yml` entry.

## 2026-11-05 — Run a literature search against the question, not the topic
- **Pattern:** Round 3 read Vauthey 2000 for mechanism, filed it as "microemulsion", and skipped it. A forward-citation sweep keyed to H0's question (a second laboratory's absolute MFT / FFT) found it as the best lead on disk for lifting H0's refusal.
- **Rule:** a case's sweep is keyed to its own decisive criterion (sec. 7 item 2). A source skipped for one question gets re-screened for another; record the skip reason with the question it answered.
