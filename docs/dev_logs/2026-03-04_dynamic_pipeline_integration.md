# Dev Log: Phase 7.1 End-to-End Pipeline Integration

## Date: 2026-03-04

### Summary
Today we completed the architectural transition from static, curated pathways (Phases 2-5) to a fully dynamic, combinatorial generative engine (Phase 7.1). By connecting the `SmirksEngine` directly to the `Recommender`, the pipeline can now evaluate non-deterministic, highly branched reaction networks on the fly.

### Architectural Decisions

#### 1. Shift from Static Pathways to Combinatorial Hypergraphs
- **Context:** The MVP relied on hardcoded sequences (`curated_pathways.py`) to bypass RMG limitations. While informative for tuning the physics engine (xTB), it could not discover novel cross-reactivity between diverse precursors.
- **Decision:** The `Recommender` now accepts a dynamically generated `List[ElementaryStep]` and constructs a reaction hypergraph to identify targets. This transforms the tool from an analysis script into a true exploratory engine.

#### 2. Bottleneck Ranking via Hypergraph Relaxation
- **Context:** Standard graph traversal (like Dijkstra) fails for multi-reactant steps (e.g., Strecker, Aminoketone Condensation). A shortest-path algorithm might find a low-barrier path to `Reactant A`, but if `Reactant B` requires a massive 60 kcal/mol barrier to form, the true kinetic bottleneck for the step must be the maximum of both.
- **Decision:** We implemented a modified Bellman-Ford topological relaxation in `predict_from_steps()`. The "distance" (kinetic bottleneck) to a product node is defined as `max(max(reactant_distances), step_intrinsic_barrier)`. This ensures that computationally predicted pathways are physically realistic and correctly bottlenecked by their most restricted precursor branch.

#### 3. Dual-Mode Evaluation (Fast vs Rigorous)
- **Context:** Full GFN2-xTB geometric relaxation and NEB calculations take minutes to hours representing a severe UX bottleneck for formulation scientists needing rapid iteration.
- **Decision:** We architected a dual-mode `run_pipeline.py`. 
    - **Fast Mode (Default):** Applies rule-based Hammond heuristic fallbacks (~15-40 kcal/mol depending on family) modified by pH conditions, providing instant sensory feedback (fractions of a second).
    - **Rigorous Mode (`--xtb`):** Triggers the complete semi-empirical physics engine for publishable or novel results.

### Current Status
The generative architecture successfully supports arbitrary precursor pools. Running "Ribose + Cysteine" via the FAST heuristic mode instantly enumerates steps and ranks 2-Furfurylthiol (FFT) as the dominant meaty aroma target.

### Update: Phase 7 Completion (Inverse Design & PBMA Metrics)
Following the core pipeline integration, we extended the engine to support **Plant-Based Meat Alternative (PBMA) industrial workflows**:
1. **PBMA Additives & Trapping (`7.2` & `7.4`)**: Integrated Thiamine degradation, Glutathione cleavage, and Hexanal/Nonanal lipid oxidation products. Implemented a **Lipid Trapping Efficiency** metric to quantify how well amino acids scavenge off-flavor aldehydes into Schiff bases.
2. **Lysine Budget / DHA Competition (`7.6`)**: Added a stoichiometric tracker that compares the kinetic flux of Lysine towards Maillard (flavor) vs. Dehydroalanine crosslinking (texture), exposing the physical competition inherent in plant protein extrusion.
3. **Inverse Design Engine (`7.5`)**: Built an optimization layer (`src/inverse_design.py`) that evaluates a matrix of predefined industrial formulations (`formulation_grid.yml`) against desired sensory profiles (`sensory_tags.yml`). By running the FAST engine across the grid and scoring the generated targets (e.g., maximizing `meaty` while minimizing `beany`), the tool successfully acts as a recommendation engine for food scientists. 

Phase 7 is officially complete. The pipeline now functions as a complete end-to-end sandbox for both forward simulation and inverse commercial design.
