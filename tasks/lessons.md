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
