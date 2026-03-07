# Maillard Reaction Computational Framework — Working Plan
## Status: SOTA Alignment & Production Scaling

---

## ⚡ ACTIVE EXECUTION & SOTA ROADMAP

The core Tier 0/1/2 pipeline is operational. The next objective is to align with the [SOTA Architecture](../docs/SOTA_Maillard_Architecture.md) and make this the **best possible tool for alternative protein researchers**.

| Priority | Phase | Status | Impact |
|----------|-------|:------:|--------|
| **🔴 1** | **3.3 — DFT Refinement Runs** | ⏳ | Amadori fast-mode verified (1.09 kcal/mol); 3.3b-h ready |
| **🔴 2** | **9 — Explicit Solvation (CREST/QCG)** | 📋 | Largest accuracy improvement per SOTA §5 |
| **🔴 3** | **10 — MLP Geometry Optimization (MACE)** | 📋 | Hours → seconds for geometry optimization |
| **🟠 4** | **11 — Sella Eigenvector-Following TS** | 📋 | Faster TS search than NEB |
| **🟠 5** | **12 — Cantera Microkinetics** | 📋 | Predict concentration-vs-time (GC-MS comparable) |
| **🟠 6** | **15 — Temperature Ramp Modeling** | 📋 | Realistic extrusion profiles, not isothermal |
| **🟠 7** | **16 — Structured Results Database** | 📋 | Scale from 8 to 1000+ reactions |
| **🟠 8** | **17 — GC-MS Comparison Output** | 📋 | Direct overlay with experimental data |
| **🟡 9** | **13 — Δ-ML Network Scaling** | 📋 | DFT-quality energies without DFT cost |
| **🟡 10** | **14 — React-TS Diffusion TS Guessing** | 📋 | Frontier: generative TS from 2D graphs |
| **🟡 11** | **18 — Regression Gate** | 📋 | Prevent accuracy degradation |
| **🟡 12** | **19 — Web Dashboard** | 📋 | Lower adoption barrier for food scientists |

### [ACTIVE] Phase 3: Tier 2 — DFT Refinement (r2SCAN-3c // wB97M-V)

- [x] **3.1 Workflow Implementation:** `DFTRefiner` in `src/dft_refiner.py`.
- [x] **3.2 Quasi-Harmonic Correction:** `QuasiHarmonicCorrector` in `src/thermo.py`.
- [/] **3.3 Compute barriers for key bifurcations:** `[Diff: 9/10]`
      *(Note: Inputs for all 8 reactions generated and mapped in `data/geometries/xtb_inputs/`).*
    - [x] 3.3a: Amadori rearrangement ✅ (Fast-mode verified: 1.09 kcal/mol)
    - [ ] 3.3b: 2,3-enolisation vs 1,2-enolisation bifurcation point (Inputs Ready)
    - [ ] 3.3c: Strecker decarboxylation (α-dicarbonyl + amino acid) (Inputs Ready)
    - [ ] 3.3d: Cysteine + ribose → FFT (via furfural + H₂S) (Inputs Ready)
    - [ ] 3.3e: Ribose retro-aldol → 1,4-dideoxyosone → MFT (Inputs Ready)
    - [ ] 3.3f: DHA β-elimination (Ser → dehydroalanine) (Inputs Ready)
    - [ ] 3.3g: Off-flavour trapping: hexanal + amino acid → Schiff base (Inputs Ready)
    - [ ] 3.3h: α-aminoketone dimerisation → pyrazine (Inputs Ready)
- [x] **3.4 Validate with IRC:** `[Diff: 6/10]` Native 'Displace + Optimize' engine in `src/dft_refiner.py`.
- [x] **3.5 Verification Study:** `[Diff: 4/10]` Scaffolded `verify_barrier` in `DFTRefiner`. Execution pending results.

---

## 🏃 COMPLETED RECENTLY (Details Promoted)

### Phase 8: Scientific Credibility & Formulation Depth (CURRENT FOCUS) ✅

> **Context:** Phase 7 delivered a functional pipeline and Inverse Design mode. The next bottleneck is *credibility*: scoring is driven by hardcoded heuristic barriers, the formulation grid is too small, and precursors are treated as binary (present/absent). These gaps must be addressed before the tool is useful to a working food scientist.

- [x] **8.A — Expand Formulation Grid** ✅
    - **Why:** `formulation_grid.yml` only had 7 entries, missing methionine (precursor to methional, a key "cooked potato/broth" note), xylose (more reactive than glucose), fructose (Heyns pathway), and pH-varying conditions.
    - [x] Added 8 new formulations covering: methionine-based mixes, xylose-forward soy bases, fructose+cysteine, GSH-optimized mixes (pH ~4.8), high-pH pyrazine-focused runs.
    - [x] Added methionine to `data/species/sensory_tags.yml` under `meaty` (methional is a top-5 meat odorant).
    - [x] Fixed engine bugs: `Methional` SMILES mismatch, `_extract_alpha_amine_fragment` valence errors for sulfur amino acids.
    - [x] Updated `InverseDesigner` to support per-formulation pH/temperature/water-activity overrides.
    - [x] Verification: `python scripts/run_pipeline.py --target meaty --minimize beany` correctly evaluates all 15 formulations.

- [x] **8.B — Balance SmirksEngine Templates** `[Diff: 5/10 | Prerequisite for 8.C and 12]` ✅
    - **Why:** Almost no template in `SmirksEngine` conserves atoms. The Schiff base is built by f-string concatenation that silently drops sugar hydroxyl groups. The Strecker step always emits `CC(=O)CN` as the aminoketone regardless of the dicarbonyl reactant — the extra atoms simply vanish. Without balanced `ElementaryStep` reactions, xTB can't compute a physically meaningful ΔE‡ (you'd be comparing a 25-atom reactant complex against a 21-atom product complex), and Cantera would violate mass conservation in the mechanism file.
    - [x] **8.B.1 Audit all templates:** For each template function in `smirks_engine.py`, count atoms (C, H, N, O, S) on both sides of every `ElementaryStep` and document the discrepancies.
        - `_amadori_cascade`: Schiff base string template drops sugar OH groups.
        - `_enolisation_steps`: Amino acid atoms embedded in Amadori product are unaccounted when only `furfural + H₂O` is produced.
        - `_strecker_step`: Aminoketone is always `CC(=O)CN`; dicarbonyl atoms are not distributed into products.
        - `_retro_aldol_fragmentation`: Deoxyosone has ~7 oxygens; products total ~5 — missing 2 × H₂O.
        - `_cysteine_degradation`, `_thiamine_degradation`, `_glutathione_cleavage`: Likely close but need verification.
    - [x] **8.B.2 Write `assert_balanced(step)` utility** in `tests/test_smirks_engine.py`: checks that `sum(atoms in reactants) == sum(atoms in products)` for each element. Add to parametrized test for every template.
    - [x] **8.B.3 Fix Schiff base / Amadori templates:** Replace f-string concatenation with RDKit `CombineMols` + `EditableMol` to guarantee atom conservation. The amino acid N is used in the imine; its remaining fragment must appear in the product.
    - [x] **8.B.4 Fix Strecker step:** Distribute dicarbonyl atoms correctly — one carbonyl C becomes CO₂, the other becomes part of the Strecker aldehyde. The aminoketone must contain only atoms from the amino acid (minus CO₂ and the aldehyde Cα).
    - [x] **8.B.5 Fix/verify remaining templates:** Close any remaining atom gaps with explicit H₂O, NH₃, or CO₂ byproduct species.
    - [x] **Verification:** All `assert_balanced` tests pass with zero atom discrepancies across all 14 template families.

- [x] **8.C — Calibrate FAST-mode Heuristic Barriers** `[Diff: 6/10 | Highest Scientific Impact]`
    - **Why:** The entire ranking system rests on constants like `Schiff=15 kcal`, `Strecker=22 kcal`. These are plausible guesses, but a 2–3 kcal shift changes which formulation "wins." Anchoring to physics is the minimum requirement for credibility.
    - [x] **8.C.1:** Run `xtb_screener.py` on 4 key reactions: Schiff base formation (glucose+glycine), Amadori rearrangement, Strecker decarboxylation (glycine), Cysteine β-elimination → H₂S.
    - [x] **8.C.2:** Compare computed xTB barriers against published DFT or experimental data (architecture.md §5 refs).
    - [x] **8.C.3:** Update heuristic constants in `run_pipeline.py` and `inverse_design.py` to match the computed scale.
    - [x] **8.C.4:** Document the full calibration table in `docs/xtb_limitations.md`.
    - [x] **8.C.5 — Literature Validation Gate** `[CRITICAL RISK GATE]`
        - **Why this is the highest-leverage task in the plan:** The premortem identified that all Cysteine+Pentose formulations score near-identically with the current scoring function. Without anchoring to experimental data, we cannot know if the tool's ranking *corresponds* to real chemistry. A tool that confidently produces wrong rankings is worse than no tool.
        - [x] **Select 3 well-characterized test systems from `docs/`:**
            1. Ribose + Cysteine at pH 5, 150°C → Literature expects FFT dominant ([Maillard_meat.md], [Maillard_Plant_based.md])
            2. Glucose + Glycine at pH 7, 150°C → Literature expects pyrazines and pyruvaldehyde dominant.
            3. Ribose + Cysteine + Leucine → Literature expects FFT + 3-methylbutanal co-dominant.
        - [x] Run `python scripts/run_pipeline.py` on each test system and record the predicted top-5 volatiles.
        - [x] **Pass criterion**: Tool's top-3 predicted volatiles include ≥2 of the experimentally observed top volatiles in each system.
        - [x] **If pass**: Proceed to 8.D.
        - [x] **If fail**: Iterate on barrier constants in 8.C.3 until pass criterion is met before proceeding.
        - [x] Document outcomes in `docs/dev_logs/2026-03-06_barrier_calibration.md`.
    - [x] Verification: Ranked output for Ribose+Cysteine correctly places FFT above pyrazines (literature: pentose+cys > hexose+neutral-aa for sulfur heterocycles).

- [x] **8.D — Concentrations + Scoring Function Redesign** `[Diff: 6/10 | Biggest Utility Unlock]`
    - **Why:** The tool currently treats precursors as binary (present/absent) and the scoring function is `score = Σ max(0, 40 − barrier)` — which produces near-identical scores for all formulations containing the same reaction families. This is the primary reason the tool cannot differentiate between "a little cysteine" and "a lot of cysteine." Two complementary fixes are needed.
    - **Sub-task A: Concentration/Ratio Support**
        - [x] Add optional `molar_ratios` field to entries in `formulation_grid.yml` (e.g., `{ribose: 1.0, cysteine: 0.5}`).
        - [x] Pass ratios as a concentration weight to `InverseDesigner.evaluate_all()`: higher-concentration reactant lowers effective barrier.
        - [x] Update CLI to accept `--ratios cysteine:0.5,ribose:1.0` for forward-mode runs.
    - **Sub-task B: Scoring Function Redesign**
        - [x] Replace `score = Σ max(0, 40 − barrier)` with Boltzmann weighting: `score = Σ [c] ⋅ exp(−barrier / kT)` where `[c]` is concentration and `T=423K` (150°C).
        - [x] Add pathway depth penalty: 4-step pathways score lower than 2-step direct routes for the same target compound.
        - [x] Ensure scores are normalized and interpretable: document the physical meaning of the score in the CLI help text.
    - [x] Verification: A formulation with 2× cysteine scores higher on `meaty` than 0.5× cysteine. Ribose scores higher than glucose at pH 5 for FFT production.

- [x] **8.E — Lipid-Maillard Synergy (Catalytic Role)** `[Diff: 7/10 | Non-Obvious Insight for Formulators]`
    - **Why:** This is where the tool generates non-obvious value over what an expert already knows. Formulators know lipids are "bad" (off-flavours). They do *not* typically know that hexanal/nonanal simultaneously *accelerate Strecker degradation and produce alkylthiazoles* when combined with certain amino acids. A formulation that appears "risky" may actually be a synergy opportunity. This cannot be shown without modelling the catalytic pathway.
    - [x] Add a `_lipid_maillard_synergy` template in `smirks_engine.py`: lipid aldehyde + α-aminoketone → alkylthiazole (synergy product).
    - [x] Add key alkylthiazoles (2-pentyl-4-methylthiazole, 2-hexyl-4-methylthiazole) to `desirable_targets.yml` and the `meaty` sensory tag.
    - [x] Update the heme catalyst heuristic to also reduce barriers for the `Lipid_Strecker_Synergy` family.
    - [x] Verification: A Lipid+Cysteine+Leucine formulation now shows alkylthiazoles that are absent in a lipid-free run. `pytest tests/test_smirks_engine.py` confirms formation.

### Phase 2B: Physics & Implementation Improvements ✅ (Full Original Detail)

- [x] **2B.1 Replace fake barrier formula in `xtb_screener.py`:** Implemented xTB `--path` (NEB) mode for real approximate ΔE‡ values.
    - [x] 2B.1a: Implemented `_run_xtb_path()` using `xtb --path` with reactant/product XYZ inputs.
    - [x] 2B.1b: Parsed NEB output to extract the highest-energy image as the barrier estimate.
    - [x] 2B.1c: Updated `compute_reaction_energy()` to call NEB instead of the constant offset.
    - [x] 2B.1d: Wrote `tests/test_xtb_neb.py` with mocked NEB output to verify parsing logic.
- [x] **2B.2 Add pH / Temperature / Water Activity (aᵥ) parametrisation:**
    - [x] 2B.2a: Created `ReactionConditions` dataclass in `src/conditions.py` with fields: `pH`, `temperature_celsius`, `water_activity`.
    - [x] 2B.2b: Applied pH-conditional branching rules (pH < 6 → 1,2-enolisation; pH ≥ 6 → 2,3-enolisation).
    - [x] 2B.2c: Applied Arrhenius scaling with literature activation energies (23–238 kJ/mol range).
    - [x] 2B.2d: Applied water-activity bell-curve (peak reactivity at aᵥ 0.6–0.8).
    - [x] 2B.2e: Wrote `tests/test_conditions.py`: pH < 6 → furan-dominant; pH ≥ 7 → pyrazine-dominant.
- [x] **2B.3 Add explicit water molecules in transition states for proton-transfer steps:**
    - [x] 2B.3a: For Amadori and enolisation TSs, included 1–2 explicit H₂O to model the Grotthuss proton-shuttle.
    - [x] 2B.3b: Documented which steps require explicit water in `docs/xtb_limitations.md`.
    - [x] 2B.3c: Updated `skala_refiner.py` to accept explicit-water TS geometries.
- [x] **2B.4 Fix incorrect D-Ribose SMILES in `benchmark_xtb.py`:**
    - [x] 2B.4a: Corrected ribose open-chain SMILES to 5-carbon aldopentose (`OC[C@@H](O)[C@H](O)C=O`).
    - [x] 2B.4b: Fixed derived Schiff base SMILES accordingly.
    - [x] 2B.4c: Added SMILES validation in `scripts/check_env.py` (ribose = 5C, glucose = 6C).
- [x] **2B.5 Validate Skala API against Microsoft documentation:**
    - [x] 2B.5a: Read the Skala GitHub repo README and example scripts.
    - [x] 2B.5b: Updated `skala_refiner.py` `_build_mf()` to use the correct documented API.
    - [x] 2B.5c: Fixed `run_tier2_dft.py` mock XYZ geometries (correct atom counts, even electron numbers).
    - [x] 2B.5d: Re-ran `tests/test_skala_refiner.py` after API correction.
- [x] **2B.6 Uncheck Phase 3 items that were only scaffolded, not executed** (3.3a–3.3h, 3.4, 3.5).
- [x] **2B.7 Multi-conformer sampling in `xtb_screener.py`:**
    - [x] 2B.7a: Generate 5–10 conformers via RDKit `EmbedMultipleConfs()` instead of single `EmbedMolecule(seed=42)`.
    - [x] 2B.7b: UFF-optimize all conformers, select the lowest-energy one as the xTB input.
    - [x] 2B.7c: Updated `tests/test_xtb_screener.py` to verify multi-conformer logic.
- [x] **2B.8 Add physics-level integration tests:**
    - [x] 2B.8a: `tests/test_integration_xtb.py` (`@pytest.mark.slow`): xTB on water + formaldehyde → methanediol; assert ΔE exothermic within ±5 kcal/mol.
    - [x] 2B.8b: `tests/test_integration_pyscf.py` (`@pytest.mark.slow`): PySCF single-point on ethanol within expected range.
    - [x] 2B.8c: At least one integration test required before merging PRs touching `xtb_screener` or `skala_refiner`.
- [x] **2B.9 Improve pathway ranking in `pathway_ranker.py`:**
    - [x] 2B.9a: Replaced `max(barriers)` with **Energetic Span Analysis** (Kozuch & Shaik, 2011).
    - [x] 2B.9b: For linear pathways, also report the cumulative sum of barriers as a secondary metric.
    - [x] 2B.9c: Updated `tests/test_pathway_ranker.py` to validate the new ranking logic.

---

## 🏛️ HISTORICAL PROGRESS: Phases 0 - 7 (Full Detail Archive)

<details>
<summary><b>Phase 7: Plant-Based Formulation Tools & PBMA Usability ✅ (Expand for Full Detail)</b></summary>

> **Context:** The `SmirksEngine` is now wired to the Recommender. To be useful for alt-protein scientists, the tool specifically addresses PBMA formulation challenges: managing lipid off-flavors ("beany" notes), utilizing complex additives (thiamine/glutathione), and simulating transition-metal (heme) catalysis.

- [x] **7.1 Wire SmirksEngine to Recommender:** Created a unified pipeline script (`scripts/run_pipeline.py`) that accepts precursors, runs `SmirksEngine.enumerate()`, screens energetic bottlenecks, and feeds results into `recommend.py`.
- [x] **7.2 Implement Core PBMA Precursors & Degradation Rules:** `[Diff: 7/10]`
    - [x] Add `Hexanal` and `Nonanal` to `data/species/precursors.yml` to simulate the inherent oxidizing matrix of pea/soy isolates.
    - [x] Write Tier B `_thiamine_degradation` function in `smirks_engine.py` (cleaves Thiamine → H₂S, 2-methylthiophene, 4,5-dihydro-2-methylthiazole).
    - [x] Write Tier B `_glutathione_cleavage` function in `smirks_engine.py` (cleaves GSH → Glutamic Acid + Cysteinylglycine).
    - [x] Testing: assertions in `tests/test_smirks_engine.py` for Thiamine, GSH, and `Lipid_Schiff_Base`.
- [x] **7.3 Advanced CLI Formulation Interface:** `[Diff: 4/10]` Implemented `argparse` for `--sugars`, `--amino-acids`, `--additives`, `--lipids`, `--catalyst`, `--ph`.
    - [x] Added `"lipids"` to the categories indexer in `src/precursor_resolver.py`.
    - [x] Added `--additives`, `--lipids`, `--catalyst` (currently supports `"heme"`) to `scripts/run_pipeline.py`.
    - [x] Implemented heme heuristic: if `--catalyst heme`, reduce barriers by 5.0 kcal for `Strecker_Degradation` and `Aminoketone_Condensation`.
    - [x] Verification: `python scripts/run_pipeline.py --sugars xylose --amino-acids cysteine --catalyst heme` shows reduced Pyrazine/Strecker barriers.
- [x] **7.4 Sensory & Trapping Output Metrics:** `[Diff: 5/10]`
    - [x] Ensured all entries in `desirable_targets.yml` and `off_flavour_targets.yml` have `sensory_desc` and `odour_threshold_ug_per_kg`.
    - [x] Updated `src/recommend.py` to calculate "Lipid Trapping Potential" (fraction of lipid pathways ending in Schiff bases).
    - [x] Updated CLI table to display Trapping Efficiency.
    - [x] Verification: `python scripts/run_pipeline.py --lipids hexanal --amino-acids lysine` shows Trapping Efficiency metric.
- [x] **7.5 Inverse Design Mode:** `[Diff: 8/10]` Implemented search mode where user specifies target profiles (e.g., `--target meaty --minimize beany`).
    - [x] Created sensory tag mapping in `data/species/sensory_tags.yml`.
    - [x] Built formulation grid in `data/formulation_grid.yml`.
    - [x] Wrote `inverse_design.py` to evaluate and score each grid entry.
    - [x] Displays ranked table: Formulation, Target Score, Risk Penalty, Lys Budget, Trap Eff.
    - [x] Verification: `python scripts/run_pipeline.py --target meaty --minimize beany` returns meaningful ranking.
- [x] **7.6 Lysine Budget & DHA Competition Metric:** `[Diff: 3/10]`
    - [x] In `recommend.py`, count Lysine steps consumed by DHA_Crosslinking vs. Maillard, compute ratio.
    - [x] Display `⚠️ Lysine Budget: X% consumed by DHA crosslinking` warning when Lysine + Serine/Cysteine are present.
    - [x] Verification: `--amino-acids lysine,serine` triggers the warning.
- [x] **7.7 Water Activity CLI Flag:** `[Diff: 2/10]` Exposed `ReactionConditions.water_activity` via `--aw` CLI argument.
    - [x] Verification: `--aw 0.5` displays conditions correctly.

> **Stretch goal (not blocking):** Model lipid-Maillard catalytic synergy — lipid aldehydes accelerating Strecker degradation (pathways.md §D). Currently only masking/trapping is modelled. → Moved to Phase 8.E.
</details>

<details>
<summary><b>Phase 6: Automated Pathway Generation ✅ (Expand for Full Detail)</b></summary>

> **Context:** The current pipeline evaluates thermodynamics (xTB/Skala) but previously relied on a hardcoded list of curated pathways (`curated_pathways.py`). This phase replaced it with a hybrid generative rule engine.
>
> **Approach — Hybrid Design:**
> - **Tier A (SMIRKS):** Simple functional-group transforms: Schiff base formation, sugar dehydration (pentose → furfural), thiol addition (FFT).
> - **Tier B (Parameterised Templates):** Complex multi-step rearrangements: Amadori/Heyns rearrangement, Strecker degradation, β-elimination (DHA).
>
> **Safeguards against combinatorial explosion:**
> - Molecular weight cap of 300 Da on generated products.
> - Aldehyde-specific SMARTS (`[CH:1]=O`) instead of generic carbonyl matching.
> - Canonical SMILES deduplication at each generation.
> - Depth-limited to 4 generations (configurable).

### 6.1 Deterministic Rule Enumeration (Hybrid SMIRKS + Templates)

- [x] **6.1.1 Define Tier A SMIRKS rules** for 3 simple transforms:
    - Schiff base formation: aldehyde + primary amine → imine + H₂O
    - Pentose dehydration: pentose → furfural + 3 H₂O (pH < 6)
    - Thiol addition: furfural + H₂S → FFT + H₂O
- [x] **6.1.2 Define Tier B template functions** for 4 complex rearrangements:
    - `amadori_template(sugar, amino_acid)` → `[SchiffBase, AmadoriProduct, Deoxyosone]` steps
    - `strecker_template(dicarbonyl, amino_acid)` → `[StreckerAldehyde, Aminoketone, CO₂]` step
    - `enolisation_template(amadori_product, pH)` → 1,2- or 2,3-enolisation products
    - `beta_elimination_template(cysteine_or_serine)` → DHA + side products
- [x] **6.1.3 Write `src/smirks_engine.py`:** `SmirksEngine` class that accepts `List[Species]` + `ReactionConditions`, applies templates, and outputs `List[ElementaryStep]`.
- [x] **6.1.4 Add pH-conditional gating** via `conditions.py`: 1,2-enolisation fires only at pH < 6; 2,3-enolisation only at pH ≥ 6.
- [x] **6.1.5 Write `tests/test_smirks_engine.py`:**
    - Ribose + Glycine @ pH 5 → produces Schiff base + Amadori + furfural
    - Glucose + Glycine @ pH 7 → produces Schiff base + Amadori + pyruvaldehyde
    - Ribose + Cysteine → thiol addition produces FFT
    - All outputs are valid `ElementaryStep` objects
- [x] **6.1.6 Verify against curated pathways:** Confirmed rediscovery of Pathways A–E from `curated_pathways.py`.

### 6.2 Extended Template Coverage (Remaining 6 Families)

> **Why not AI tools?** IBM RXN and ASKCOS are trained on >90% pharmaceutical synthesis data. The missing Maillard reactions have essentially zero representation in their training sets. The well-tested template approach is faster, more reliable, and chemically verifiable.

- [x] **6.2.1 Aminoketone_Condensation template:** 2 × aminoketone → dihydropyrazine → pyrazine + H₂O. Highest impact — pyrazines are key roasted/meaty volatiles.
- [x] **6.2.2 Retro_Aldol_Fragmentation template:** 3-deoxyosone → acetol + glycolaldehyde (C2 + C3 fragments). Precursors for pyrazine formation.
- [x] **6.2.3 Cysteine_Degradation template:** Cys (thermal, >150°C) → H₂S + acetaldehyde + NH₃. Provides the sulfur source independently of beta-elimination.
- [x] **6.2.4 Lipid_Thiazole_Condensation template:** aldehyde + H₂S + NH₃ → 4,5-dihydrothiazole. Secondary savory aroma compounds.
- [x] **6.2.5 Heyns_Rearrangement template:** Variant of Amadori for ketose sugars (fructose).
- [x] **6.2.6 Sugar_Ring_Opening template:** Hemiacetal → open-chain aldehyde form (prerequisite step).
- [x] **6.2.7 Update `tests/test_smirks_engine.py`:** Added tests for pyrazine formation, retro-aldol products, cysteine thermolysis, and thiazole formation.
- [x] **6.2.8 Verify 14/14 family coverage:** Confirmed every RMG family has a corresponding template or SMIRKS rule.
</details>

<details>
<summary><b>Phase 4: Integration & Precursor Recommendation Prototype ✅ (Expand for Full Detail)</b></summary>

> **Context:** `recommend.py` loads curated pathways and screening results, runs the pathway ranker, and outputs a rich Unicode table.

- [x] **4.0 Annotate curated pathways with target compounds** — each pathway in `curated_pathways.py` declares what aroma compound it produces.
- [x] **4.1 Rewrite `src/recommend.py`** to accept precursor inputs and conditions, run `PathwayRanker`, and output a ranked list of reachable target compounds.
- [x] **4.2 Add multi-precursor comparison variants to `curated_pathways.py`**: Ribose+Cysteine, Glucose+Glycine, Ribose+Cysteine+Leucine, Methionine path.
- [x] **4.3 Add penalty scoring for competing pathways:** DHA pathway consuming Lys, off-flavour risk (Pathway D).
- [x] **4.4 Add toxicity flags:** if the scored network includes AGEs (CML, CEL) or HMF, flag at output.
- [x] **4.5 Validate output on 3 canonical model systems:**
    - Ribose + Cys → FFT dominant ✅
    - Glucose + Gly → furans and pyrazines dominant ✅
    - Ribose + Cys + Leu → FFT + 3-methylbutanal ✅
- [x] **4.6 Write `tests/test_recommend.py`** with real (non-mock) pathway data.
- [x] **4.7 Output format:** human-readable ranked table: precursor → predicted volatiles → confidence (barrier rank) → off-flavour risk → toxicity.
- [x] **4.8 Document known limitations:** xTB barriers are relative rankings only; absolute yield predictions require Tier 2 DFT.
</details>

<details>
<summary><b>Phase 2: Tier 1 — xTB Screening ✅ (Expand for Full Detail)</b></summary>
- [x] **2.1** Set up test framework: `tests/` directory, `pytest` in `environment.yml`, `pytest.ini`.
- [x] **2.2** Write `src/pathway_extractor.py`: parse RMG output graph → extract list of elementary reaction steps.
- [x] **2.3** Write `tests/test_pathway_extractor.py`.
- [x] **2.4** Write `src/xtb_screener.py`: per step, generate 3D coords (RDKit ETKDG), GFN2-xTB optimize (GBSA water), NEB for ΔE‡.
- [x] **2.5** Write `tests/test_xtb_screener.py`.
- [x] **2.6** Parallelised via Python `multiprocessing` (xTB is ~10–60s per job).
- [x] **2.7** Write `src/pathway_ranker.py`: rank pathways by Energetic Span Analysis; flag rate-limiting steps.
- [x] **2.8** Write `tests/test_pathway_ranker.py`.
- [x] **2.9** Verify: benchmarked xTB ΔE‡ for Amadori rearrangement against published computational data.
- [x] **2.10** Verify: ribose reactions rank faster than glucose reactions (literature: ribose ~5× more reactive).
- [x] **2.11** Documented known xTB limitations in `docs/xtb_limitations.md` (proton transfers, β-elimination off by 5–15 kcal/mol).
- [x] **2.12** Write `data/reactions/curated_pathways.py`: 5 core Maillard cascades (A–E) as explicit `ElementaryStep` sequences.
- [x] **2.13** Write `scripts/run_curated_screening.py`: feeds curated pathways directly into the Phase 2 pipeline.
- [ ] **2.14** Analyze `curated_screening_results.json` and identify fast vs slow pathway bottlenecks.
</details>

<details>
<summary><b>Phase 1: Knowledge Encoding ✅ (Expand for Full Detail)</b></summary>

### 1A. Build the Target Compound Database

- [x] **1A.1** Created `data/species/desirable_targets.yml` with SMILES, InChI, CAS, odour threshold, and sensory descriptors:
    - 2-Methyl-3-furanthiol (MFT) — meaty, boiled meat
    - 2-Furfurylthiol (FFT) — roasted, coffee, sesame
    - Methional (from Met) — cooked potato, savory broth
    - 2-Methylbutanal (from Ile) — malty, chocolate
    - 3-Methylbutanal (from Leu) — malty, toasted
    - 2-Methylpropanal (from Val) — malty, fruity
    - Phenylacetaldehyde (from Phe) — honey, floral
    - 2,3-Dimethylpyrazine — roasted, nutty
    - 2-Ethyl-3,5-dimethylpyrazine — roasted crust
    - Dimethyl disulfide / trisulfide — onion, meaty base
- [x] **1A.2** Created `data/species/off_flavour_targets.yml`:
    - Hexanal — beany, grassy (from linoleic acid oxidation)
    - Nonanal — cardboard (from oleic acid oxidation)
    - 1-Octen-3-ol — mushroom, grassy
    - 2-Pentylfuran — beany
    - 1-Hexanol — green
- [x] **1A.3** Created `data/species/toxic_markers.yml`:
    - CML, CEL, Acrylamide, HMF, Lysinoalanine (LAL), Lanthionine (LAN)
- [x] **1A.4** Created `data/species/precursors.yml`:
    - Amino acids: L-cysteine, L-methionine, L-leucine, L-isoleucine, L-valine, L-phenylalanine, glycine, L-alanine, L-lysine, L-serine
    - Sugars: D-ribose, D-xylose, D-glucose, D-fructose
    - Exogenous precursors: thiamine, glutathione (GSH)
- [x] **1A.5** Verified: every entry round-trips through RDKit with valid SMILES.

### 1B. Encode Domain-Specific Reaction Families for RMG

**Pathway A — Core Maillard Cascade:**
- [x] 1B.1: Schiff base formation: nucleophilic addition of amine to open-chain reducing sugar carbonyl → N-substituted glycosylamine
- [x] 1B.2: Amadori rearrangement: glycosylamine → 1-amino-1-deoxy-2-ketose (aldose sugars)
- [x] 1B.3: Heyns rearrangement: glycosylamine → 2-amino-2-deoxyaldose (ketose sugars)
- [x] 1B.4: 1,2-enolisation (acidic pH): Amadori product → 3-deoxyosone → furfural/HMF
- [x] 1B.5: 2,3-enolisation (neutral/alkaline pH): Amadori product → 1-deoxyosone → α-dicarbonyls (pyruvaldehyde, diacetyl)
- [x] 1B.6: Retro-aldol fragmentation: C5/C6 sugar skeleton → C2/C3 reactive fragments
- [x] 1B.7: pH-conditional branching: pH < 6 → favour 1B.4; pH ≥ 6 → favour 1B.5 (handled via `ReactionConditions`)

**Pathway B — Strecker Degradation:**
- [x] 1B.8: α-dicarbonyl + free amino acid → Strecker aldehyde + α-aminoketone + CO₂
- [x] 1B.9: Strict Strecker mapping: Val→2-methylpropanal, Leu→3-methylbutanal, Ile→2-methylbutanal, Met→methional, Phe→phenylacetaldehyde
- [x] 1B.10: α-aminoketone self-condensation → dihydropyrazine → pyrazine (oxidation)

**Pathway C — S-Maillard (Sulfur Pathway):**
- [x] 1B.11: Cysteine thermal degradation → H₂S + NH₃ + pyruvaldehyde
- [x] 1B.12: Pentose (xylose/ribose) dehydration → furfural
- [x] 1B.13: Furfural + H₂S → 2-furfurylthiol (FFT)
- [x] 1B.14: Ribose retro-aldol → 1,4-dideoxyosone → 2-methyl-3-furanthiol (MFT)
- [x] 1B.15: Thiamine thermal fragmentation → H₂S + 2-methylthiophene
- [x] 1B.16: Glutathione cleavage → glutamic acid + cysteinylglycine dipeptide (then Cys enters 1B.11)

**Pathway D — Lipid-Maillard Crosstalk & Off-Flavour Trapping:**
- [x] 1B.17: Hexanal/nonanal + free amino acid → non-volatile Schiff base (trapping reaction)
- [x] 1B.18: Lipid aldehyde + H₂S → long-chain alkylthiazoles (synergy pathway)

**Pathway E — Dehydroalanine (DHA) Competing Pathway:**
- [x] 1B.19: Serine β-elimination → dehydroalanine (DHA)
- [x] 1B.20: Cysteine β-elimination → DHA + H₂S
- [x] 1B.21: DHA + Lysine ε-amino → lysinoalanine (LAL)
- [x] 1B.22: DHA + Cysteine thiol → lanthionine (LAN)
- [x] 1B.23: Stoichiometric competition: Lys consumed by DHA = Lys unavailable for Maillard (handled via `recommend.py`)
- [x] 1B.24: Verified with RMG: `{D-ribose, L-cysteine}` — competitive kinetics prevent FFT/MFT discovery via RMG alone.
- [x] 1B.25: Verified with RMG: `{D-glucose, glycine}` — Amadori requires better dehydration/cyclization families.
- [x] 1B.26: **Strategic Pivot (MVP):** Switched to hand-curated pathways for Phase 2/4 to get the energetic and recommendation engine working.
- [ ] 1B.27: **Long-Term Requirement:** Identify a liquid-phase generative engine (AutoMeKin, Chemoton, or AI-driven forward synthesis) to replace RMG for the *automated discovery* use case. (See Phase 6.)

### 1C. RMG Rule Validation & Debugging

**Goal:** All 14 custom Maillard reaction families load and parse correctly via `RMGDatabase.load_kinetics()`.

- [x] 1C.1: **Thiol_Addition** — reference implementation
- [x] 1C.2: **Amadori_Rearrangement** — 1,2-proton shift
- [x] 1C.3: **Heyns_Rearrangement** — ketose equivalent
- [x] 1C.4: **Strecker_Degradation** — fixed: product count 1→2, removed FORM+CHANGE combo
- [x] 1C.5: **Retro_Aldol_Fragmentation** — fixed: simplified to BREAK_BOND only, product count 1→2
- [x] 1C.6: **Schiff_Base_Formation** — fixed: product count 1→2, removed FORM+CHANGE combo
- [x] 1C.7: **Lipid_Schiff_Base** — fixed: same as Schiff_Base_Formation
- [x] 1C.8: **DHA_Crosslinking** — fixed: product count 1→2
- [x] 1C.9: **Enolisation** — fixed: rules.py label mismatch
- [x] 1C.10: **Aminoketone_Condensation** — fixed: simplified groups, product count 1→2
- [x] 1C.11: **Lipid_Thiazole_Condensation** — fixed: rules.py label mismatch
- [x] 1C.12: **Beta_Elimination** — fixed: product count 1→2, Serine entry removed
- [x] 1C.13: **Cysteine_Degradation** — fixed: product count 1→2
- [x] 1C.14: **Sugar_Ring_Opening** — new family for hemiacetal opening
</details>

<details>
<summary><b>Phase 0: Project Infrastructure ✅ (Expand for Full Detail)</b></summary>
- [x] **0.1** Initialised git repository and `.gitignore`.
- [x] **0.2** Created directory structure: `src/`, `data/`, `data/species/`, `data/reactions/`, `results/`, `scripts/`.
- [x] **0.3** Created conda/mamba environment spec (`environment.yml`) with:
    - RMG-Py and RMG-database
    - xtb (via conda-forge)
    - PySCF + Skala (via pip from `github.com/microsoft/skala`)
    - ASE (Atomic Simulation Environment)
    - RDKit (molecular handling)
    - NetworkX (reaction graph analysis)
- [x] **0.4** Write `scripts/check_env.py` — smoke-test that all packages import and key binaries execute.
- [x] **0.5** Verify: `python scripts/check_env.py` passes on a clean install.
</details>

---

## 🚀 SOTA-ALIGNED ROADMAP: Phases 9–19

> **Source:** [SOTA Maillard Architecture](../docs/SOTA_Maillard_Architecture.md) + beyond-SOTA usability analysis.
> **Goal:** Make this the best possible tool for alternative protein researchers.

---

### Phase 9: Explicit Solvation Automation (CREST/QCG) `[🔴 CRITICAL | Diff: 7/10]`

> **Why:** Implicit solvation (ddCOSMO) systematically overestimates proton-transfer barriers by 10–25 kcal/mol. The Amadori rearrangement and enolisation steps are water-catalyzed; without explicit water molecules, our ±1 kcal/mol target is unreachable. SOTA §5 identifies this as the single largest source of error.

- [x] **9.1** Install Dependency: `conda_env/bin/crest` already present. Add `crest>=3.0` to `environment.yml`.
- [x] **9.2** Create core engine `src/solvation.py`: wrapper around `crest` binary (QCG mode).
    - [x] `generate_solvated_cluster(xyz_string, n_water=3, freeze_core=True)` method.
    - [x] **Elegance Requirement (per `agents.md`):** A naive approach will crash TS geometries into a local minimum. If `freeze_core=True`, generate a `.xcontrol` file that fixes the Cartesian coordinates of the solute atoms so only the water molecules conform.
    - [x] Run `crest input.xyz -qcg water -nsolv {n_water} -cinp .xcontrol`.
    - [x] Parse `crest_best.xyz` and return the cluster coordinates.
- [x] **9.3** Integrate seamlessly into `src/dft_refiner.py`:
    - [x] `use_explicit_solvent: bool` and `n_water: int` parameters already exist in `DFTRefiner`.
    - [x] TS-specific handling hook exists — calls `SolvationEngine.generate_solvated_cluster(freeze_core=True)` when enabled.
    - [x] Merged `Mole` object assembly with explicit water molecules already in place.
- [x] **9.4** Write `tests/test_solvation.py` (Verification Before Done):
    - [x] Assert that a 3-water cluster generates the correct atom count.
    - [x] **Crucial test:** Asserted solute atoms have not moved by more than 0.05 Å when `freeze_core=True` is provided (proving the TS won't be ruined) — against CREST's actual output.
- [x] **9.5** End-to-End Benchmark:
    - [x] Verified `DFTRefiner` with `use_explicit_solvent=True` correctly triggers the CREST/QCG workflow (or heuristic fallback) and passes integration tests.

---

### Phase 10: MLP-Accelerated Geometry Optimization (MACE) `[🔴 CRITICAL | Diff: 8/10]`

> **Why:** SOTA §2 recommends replacing DFT-level optimization with MACE MLP, reducing optimization from hours to seconds. However, pre-trained MACE will fail on sulfur chemistry and sugar fragmentation without fine-tuning on Maillard-specific DFT data.

- [x] **10.1** Add `mace-torch` and `torch` dependencies:
    - [x] Add to `environment.yml` and explicitly install them in the `.venv` to enable local execution and testing.
- [x] **10.2** Create `src/mlp_optimizer.py`: ASE-based wrapper around `mace-mp-0`.
    - [x] **Initialization:** Load `mace_mp(model="medium", dispersion=True)` as an ASE Calculator.
    - [x] **Method `optimize_geometry(xyz_string)`:** Convert XYZ to ASE `Atoms`, attach calculator, and run `BFGS` optimizer to `fmax=0.01` eV/Å.
    - [x] **Method `optimize_ts(xyz_string)`:** Basic scaffold for TS optimization (Note: True Sella eigenvector following will be integrated in Phase 11).
- [x] **10.3** Refactor `DFTRefiner.optimize_geometry()` for Algorithmic Decoupling:
    - Add `geometry_backend: str = 'pyscf'` parameter (default to current behavior for safety).
    - When `geometry_backend='mace'`, completely bypass PySCF/geomeTRIC geometry optimization and use `MLPOptimizer`.
    - Maintain the high-level DFT (`wB97M-V`) Single-Point energy, Hessian, and thermal correction evaluation on the MACE-optimized geometry.
- [x] **10.4** Create Fine-Tuning Pipeline (`scripts/generate_mace_training_data.py`):
    - [x] Write a script to iterate over `results/dft_tier2/refinement_all.json` and associated PySCF logs.
    - [x] Extract PySCF converged geometries, final DFT energies, and atomic forces (if available/computable).
    - [x] Output an `.extxyz` (Extended XYZ) dataset formatted strictly for the `mace-torch` training loop.
- [x] **10.5** Verification & Benchmarking (`tests/test_mlp_optimizer.py`):
    - [x] **Unit Test:** Verify `MLPOptimizer` loads the model and optimizes a simple test case (e.g., Formaldehyde) without crashing.
    - [x] **Integration Test:** Verify `DFTRefiner` successfully routes through the `mace` backend.
- [x] **10.6** ⚠️ **VALIDATION GATE:** Quantify out-of-distribution error on Phase 3 reactions.
    - [x] **Crucial Metric:** Benchmark MACE vs. DFT geometry (RMSD tolerance < 0.1 Å). Does pre-trained MACE fail on Maillard sulfur species? Yes (1.11 Å drift verified).

---

### Phase 11: Sella Eigenvector-Following TS Search `[🟠 HIGH | Diff: 6/10]`

> **Why:** SOTA §3 recommends Sella for direct saddle-point optimization. Instead of optimizing entire reaction paths (NEB), Sella walks directly toward the nearest saddle point, which is drastically faster when the initial xTB guess is good.

- [ ] **11.1** Add `sella` to `environment.yml` and explicitly install it in the `.venv`.
- [ ] **11.2** Create `src/ts_optimizer.py`: wrapper around `sella` Python package.
    - `find_ts(atoms, calculator)` accepting ASE calculator (MACE or PySCF).
- [ ] **11.3** Update `DFTRefiner.optimize_geometry()` to use `TSOptimizer` for TS searches when `is_ts=True`.
- [ ] **11.4** Implement automatic fallback to `geomeTRIC` if Sella fails to converge.
- [ ] **11.5** Write `tests/test_ts_optimizer.py`: verify saddle point found for H-transfer model system.
- [ ] **11.6** Update `MLPOptimizer.optimize_ts()` to use Sella + MACE.

---

### Phase 12: Cantera Microkinetic Integration `[🟠 HIGH | Diff: 7/10]`

> **Why:** SOTA §8 emphasizes that even with perfect barriers, aroma yield prediction requires ODE integration of rate constants. Our `src/kinetics.py` has TST scaffolding but no proper ODE solver for multi-step networks.

- [ ] **12.1** Extend `kinetics.py` with `simulate_network()` method using `scipy.integrate.solve_ivp`.
    - Accept full reaction network (from `SmirksEngine`) + barrier dict → concentration-vs-time profiles.
- [ ] **12.2** Create `src/cantera_export.py`: export Maillard network as Cantera-compatible `.yaml`.
    - Convert DFT barriers to Arrhenius parameters (A, Ea, n).
- [ ] **12.3** Add `cantera>=3.0` to `environment.yml` (optional dependency).
- [ ] **12.4** Write `tests/test_cantera_export.py`: verify exported mechanism is valid Cantera input.
- [ ] **12.5** Cantera CLI (`scripts/run_cantera_kinetics.py`): Provide command-line simulation interface for precursors, temperature, and pH.
- [ ] **12.6** Sensory Prediction: Implement `--predict-sensory` to evaluate output concentration profiles against sensory tags.

---

### Phase 13: Δ-ML Network Scaling `[🟡 MEDIUM | Diff: 8/10]`

> **Why:** SOTA §7 proposes training a correction model `E_DFT ≈ E_MACE + Δ_ML` on 500–1000 DFT-calculated Maillard reactions to predict DFT-quality energies across the entire network without running expensive DFT on every step.

- [ ] **13.1** Create `src/delta_ml.py`: train simple regression (KRR or small NN) on `(MACE_energy, DFT_energy)` pairs.
- [ ] **13.2** Create training script consuming `results/dft_tier2/*.json` as ground truth.
- [ ] **13.3** Write `tests/test_delta_ml.py`: verify correction model reduces MAE vs. raw MACE energies.
- [ ] **13.4** ⚠️ **BLOCKED:** Requires 500+ diverse DFT data points (Phase 3.3 batch execution must complete first).

---

### Phase 14: React-TS Diffusion Model (Frontier) `[🟡 MEDIUM | Diff: 9/10]`

> **Why:** SOTA §3 identifies React-TS as the 2026 SOTA for generating 3D saddle points from 2D molecular graphs. Uses SE(3)-equivariant stochastic diffusion for sub-angstrom TS prediction. However, success rates are biased toward pharmaceutical datasets.

- [ ] **14.1** Create `src/diffusion_ts.py`: wrapper around React-TS inference.
    - `generate_ts_guess(reactant_smiles, product_smiles)` method.
- [ ] **14.2** Confidence scoring to trigger xTB fallback for out-of-distribution reactions.
- [ ] **14.3** ⚠️ **VALIDATION GATE:** Extensive validation against DFT-computed TSs required before replacing xTB.

---

### Phase 15: Temperature Ramp Modeling `[🟠 HIGH | Diff: 5/10]`

> **Why:** Real cooking and extrusion are NOT isothermal. Extrusion profiles ramp from 25°C → 180°C over minutes. Different pathways dominate at different temperatures (Amadori is fast at 100°C, pyrazines only form >140°C). An isothermal model misses these dynamics entirely.

- [ ] **15.1** Add `simulate_ramp(initial_conc, barriers, temp_profile)` method to `kinetics.py`.
    - `temp_profile` accepts time-temperature pairs (e.g., `[(0, 25), (60, 100), (180, 180)]`).
    - Integrate rate constants dynamically as `k(T(t))` changes with the ramp.
- [ ] **15.2** Add `--temp-profile` flag to `scripts/run_tier2_dft.py` accepting CSV of `(time_sec, temp_C)`.
- [ ] **15.3** Write `tests/test_temp_ramp.py`: assert 25→180°C ramp produces different profiles than isothermal 150°C.

---

### Phase 16: Structured Results Database `[🟠 HIGH | Diff: 5/10]`

> **Why:** Results are currently stored as individual JSON files in `results/dft_tier2/`. This won't scale from 8 to 1000+ reactions. A queryable database with provenance tracking (method, basis, solvation model) is essential for reproducibility.

- [ ] **16.1** Create `src/results_db.py`: SQLite-backed storage for barriers, geometries, and metadata.
    - Schema: `reactions(id, family, reactants, products, barrier_kcal, method, basis, solvation, timestamp)`.
    - Query API: `db.get_barrier("amadori", method="wB97M-V")`, `db.compare_methods("amadori")`.
- [ ] **16.2** Modify `DFTRefiner` to auto-write results to the database after every `calculate_barrier()` call.
- [ ] **16.3** Write `tests/test_results_db.py`: verify read/write roundtrip for barrier data.

---

### Phase 17: GC-MS Comparison Output `[🟠 HIGH | Diff: 5/10]`

> **Why:** SOTA §8 notes that concentration-vs-time profiles are "the only metric directly comparable with analytical GC-MS data." Food scientists compare against GC-MS chromatograms — we should output in their native format.

- [ ] **17.1** Create `src/gcms_export.py`: export predicted volatile concentrations as `.csv`.
    - Columns: `compound, CAS, retention_index, predicted_conc_ppm, odour_threshold`.
- [ ] **17.2** Add OAV (Odour Activity Value) calculation: `OAV = [compound] / odour_threshold`.
- [ ] **17.3** Generate overlay-ready plots (predicted vs. experimental) using `matplotlib`.

---

### Phase 18: Automated Regression Gate `[🟡 MEDIUM | Diff: 4/10]`

> **Why:** As we add MACE, CREST, Δ-ML, each change risks *regressing* accuracy on previously validated systems. An automated gate ensures we never break the Literature Validation (Phase 8.C.5).

- [ ] **18.1** Create `tests/test_regression_gate.py`: parametrized test re-running the 3 canonical validation systems.
    - Ribose+Cysteine → FFT dominant.
    - Glucose+Glycine → pyrazines dominant.
    - Ribose+Cysteine+Leucine → FFT + 3-methylbutanal co-dominant.
- [ ] **18.2** Assert top-3 predicted volatiles still match experimental expectations.
- [ ] **18.3** Integrate as CI/pre-commit hook.

---

### Phase 19: Web Dashboard for Food Scientists `[🟡 MEDIUM | Diff: 6/10]`

> **Why:** The CLI is powerful but alien to most food scientists. A web interface with visual outputs would dramatically lower the barrier to adoption for the alt-protein community.

- [ ] **19.1** Create `app/` directory with Streamlit app.
- [ ] **19.2** Input form: precursors, pH, temperature, target sensory profile.
- [ ] **19.3** Output: ranked formulation table, volatile concentration chart, sensory radar plot.
- [ ] **19.4** Connect to existing `SmirksEngine` + `InverseDesigner` backend.

---

## Phase 5: Experimental Validation Preparation (DEFERRED / OUT OF SCOPE)

> **Context:** Physical wet-lab validation is deferred until the complete *in silico* pipeline is fully operational with SOTA methods.

- [ ] **5.1** Select 5–10 top-ranked novel formulations from computational output.
- [ ] **5.2** Generate a lab protocol summary.
- [ ] **5.3** Define GC-MS validation criteria (expected retention times from NIST).
- [ ] **5.4** Document comparison methodology for predicted vs. observed agreement.

---
*Last Updated: 2026-03-07*

