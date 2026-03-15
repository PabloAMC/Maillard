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
- **Execution Environment**: In this repository, verification must run inside Docker using the `maillard` conda environment on Python 3.12 when the user requests it. Do not fall back to host Python tools for scientific validation.
- **Docker-First Triage**: Before changing scoring or QM code, reproduce the exact failing subset inside the Docker `maillard` environment. If the subset is already green there, do not edit production logic just to match a stale host-side failure report.
- **Scientific Inspection Ergonomics**: If a benchmark diagnostic requires long `python -c` or heredoc commands through the Docker wrapper, promote it into a named script or wrapper subcommand. Quote-heavy ad hoc commands are not a reproducible interface.
- **Optional Backend Imports**: For QM helpers that depend on optional backends like Sella/JAX, catch runtime import failures as well as missing-package errors. Unsupported CPU features must degrade to an explicit skip or fallback path, not a pytest collection error.
- **Projection Stability**: Do not use raw FAST pathway weights to allocate absolute volatile budget. They are useful for ranking but can collapse the projection onto sulfur donors or other non-benchmark carriers. Keep the budget projection based on stable endpoint filters and use weights only if they are explicitly re-scaled and benchmark-calibrated.
- **Temporal Consistency**: Do not calibrate benchmark branches until the main evaluation path and the diagnostic path propagate the same `time_minutes` value into FAST prediction. Otherwise the apparent branch errors can be artifacts of inconsistent temporal semantics rather than chemistry.
- **Stop Local Recalibration Loops**: If Docker sweeps show that no single local barrier tweak improves Hofmann, Mottram, and Farmer together, stop tuning sulfur barriers and move the work to the concentration/headspace projection layer. That is a structural blocker, not a chemistry-family constant problem.
