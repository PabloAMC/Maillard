# Scientific Validation Guide

## 1. How we validate "The Numbers"
The Maillard framework uses a **layered validation strategy** to ensure that predicted concentrations and sensory profiles align with reality.

### A. Kinetic Reliability (Pearson R)
We use `scripts/compare_sim_to_lit.py` to correlate our microkinetic yields against peer-reviewed experimental data curated in the **[Literature Benchmark Reference](file:///Users/pabloantoniomorenocasares/Developer/Maillard/data/benchmarks/maillard_validation_benchmarks.md)**.
- **Benchmark**: Ribose-Glycine model systems and other PRIMARY-tier systems.
- **Goal**: Maintain a Pearson R > 0.85 for major volatiles (MFT, 2,5-DMP, Methional).
- **Status**: Checked every release against `data/lit/ribose_glycine_2021.json` and associated benchmark data.

### B. Sensory Fidelity (Stevens' Law)
Because the human nose is logarithmic, absolute ppm values are less important than **perceived intensity (OAV)**.
- We use a square-root scaling factor ($Intensity = \text{Conc}^{0.5}$) for sensory radar charts.
- This ensures that a 10x increase in precursor concentration doesn't result in an unrealistic 10x flavor intensity, matching psychophysical reality.

### C. Safety Conservatism
Safety scores are designed to be **conservative** (penalizing risk even at low concentrations).
- We use Pareto-ranking to ensure that "Meaty" flavor is never optimized at the expense of crossing an "Acrylamide" threshold.

## 2. Accurately Describing Blind Spots
It is critical to be honest about what the tool *does not* yet know. These are documented in [docs/use_cases/pea_protein_report.md](file:///Users/pabloantoniomorenocasares/Developer/Maillard/docs/use_cases/pea_protein_report.md) and verified by `tests/scientific/test_blind_spots.py`:

- **Peptide Accessibility**: Most amino acids in plant isolates are buried inside globulin proteins. Yields will be lower for intact proteins vs. hydrolysates.
- **Volatile Partitioning**: High fiber/fat matrices trap flavors. 10 ppm of MFT in water smells much stronger than 10 ppm in a high-fiber pea patty.
- **Non-Heme Catalysis**: Transition metals (Iron) accelerate pyrazine formation faster than simple thermal kinetics suggest.

## 3. How to verify for your system
1. Run `python scripts/reproduce_use_cases.py` to see the baseline comparisons.
2. Check the `tests/scientific/` directory for regression tests against literature.
3. Compare the `Sensory Radar` against your own descriptive sensory panel data.
