# Maillard Reaction Computational Framework — Working Plan
## Status: Active Development

---

## ⚡ Prioritized Action Plan (as of 2026-03-05)

**Phase 7 & 8.A complete.** The tool now has a working forward simulation pipeline, Inverse Design mode, and an expanded 15-entry formulation grid. The next challenge is scientific credibility: the compound rankings are driven by hardcoded heuristic barriers, not real physics. Without anchoring these to validated data, scientists cannot trust the recommendations over their own intuition.

Ordered by **impact on alt-protein formulation scientists** (the end-user):

| Priority | Task | Status | Why |
|----------|------|:------:|-----|
| ~~1–7~~ | Phase 7 (all sub-tasks) | ✅ | Full forward pipeline + Inverse Design mode operational. |
| ~~8.A~~ | **8.A — Expand formulation grid** | ✅ | Grid expanded to 15 entries; covers methionine, xylose, fructose, pH variants. Engine fixed for AA fragments. |
| **🔴 Next** | **8.0 — Balance SmirksEngine templates** | ☐ | **Prerequisite for 8.B and 8.2.** Templates don't conserve atoms — Schiff base uses f-string concatenation that silently drops sugar hydroxyls; Strecker always emits `CC(=O)CN` regardless of the dicarbonyl. Without balanced reactions, xTB ΔE‡ values are physically meaningless. |
| **🔴 High** | **8.B — Calibrate FAST-mode barriers** | ☐ | Rankings are hardcoded guesses (e.g. Schiff=15, Strecker=22 kcal). Needs 8.0 first. |
| **🟡 Med** | **8.C — Add concentration/ratio support** | ☐ | Tool treats precursors as binary (present/absent). Formulators need "how much?", not just "which?". |
| **🟡 Med** | **8.D — Model lipid-Maillard synergy** | ☐ | Lipid aldehydes *catalyze* Strecker (pathways.md §D) — not just a trapping target. Currently only the masking side is modelled. |
| 🟢 Low | 3.3 — Full DFT barriers (Skala) | ☐ | Replaces all heuristics with chemically accurate kinetics. High effort; needed before publication-quality output. |
| 🟢 Low | 8.2 — Cantera Kinetic Simulation | ☐ | Time/temperature cooking profiles. Needs 8.0 and 8.B first. |

> **Deferred:** Phase 5 (Experimental validation), Phase 8.1/8.2 (full kinetic modeling), AI-models (IBM RXN).

---

## 🏃 ACTIVE & UPCOMING PHASES

## Phase 8: Scientific Credibility & Formulation Depth (CURRENT FOCUS)

> **Context:** Phase 7 delivered a functional pipeline and Inverse Design mode. The next bottleneck is *credibility*: scoring is driven by hardcoded heuristic barriers, the formulation grid is too small, and precursors are treated as binary (present/absent). These gaps must be addressed before the tool is useful to a working food scientist.

- [x] **8.A — Expand Formulation Grid** ✅
    - **Why:** `formulation_grid.yml` only had 7 entries, missing methionine (precursor to methional, a key "cooked potato/broth" note), xylose (more reactive than glucose), fructose (Heyns pathway), and pH-varying conditions.
    - [x] Added 8 new formulations covering: methionine-based mixes, xylose-forward soy bases, fructose+cysteine, GSH-optimized mixes (pH ~4.8), high-pH pyrazine-focused runs.
    - [x] Added methionine to `data/species/sensory_tags.yml` under `meaty` (methional is a top-5 meat odorant).
    - [x] Fixed engine bugs: `Methional` SMILES mismatch, `_extract_alpha_amine_fragment` valence errors for sulfur amino acids.
    - [x] Updated `InverseDesigner` to support per-formulation pH/temperature/water-activity overrides.
    - [x] Verification: `python scripts/run_pipeline.py --target meaty --minimize beany` correctly evaluates all 15 formulations.

- [ ] **8.0 — Balance SmirksEngine Templates** `[Diff: 5/10 | Prerequisite for 8.B and 8.2]`
    - **Why:** Almost no template in `SmirksEngine` conserves atoms. The Schiff base is built by f-string concatenation that silently drops sugar hydroxyl groups. The Strecker step always emits `CC(=O)CN` as the aminoketone regardless of the dicarbonyl reactant — the extra atoms simply vanish. Without balanced `ElementaryStep` reactions, xTB can't compute a physically meaningful ΔE‡ (you'd be comparing a 25-atom reactant complex against a 21-atom product complex), and Cantera would violate mass conservation in the mechanism file.
    - [ ] **8.0.1 Audit all templates:** For each template function in `smirks_engine.py`, count atoms (C, H, N, O, S) on both sides of every `ElementaryStep` and document the discrepancies.
        - `_amadori_cascade`: Schiff base string template drops sugar OH groups.
        - `_enolisation_steps`: Amino acid atoms embedded in Amadori product are unaccounted when only `furfural + H₂O` is produced.
        - `_strecker_step`: Aminoketone is always `CC(=O)CN`; dicarbonyl atoms are not distributed into products.
        - `_retro_aldol_fragmentation`: Deoxyosone has ~7 oxygens; products total ~5 — missing 2 × H₂O.
        - `_cysteine_degradation`, `_thiamine_degradation`, `_glutathione_cleavage`: Likely close but need verification.
    - [ ] **8.0.2 Write `assert_balanced(step)` utility** in `tests/test_smirks_engine.py`: checks that `sum(atoms in reactants) == sum(atoms in products)` for each element. Add to parametrized test for every template.
    - [ ] **8.0.3 Fix Schiff base / Amadori templates:** Replace f-string concatenation with RDKit `CombineMols` + `EditableMol` to guarantee atom conservation. The amino acid N is used in the imine; its remaining fragment must appear in the product.
    - [ ] **8.0.4 Fix Strecker step:** Distribute dicarbonyl atoms correctly — one carbonyl C becomes CO₂, the other becomes part of the Strecker aldehyde. The aminoketone must contain only atoms from the amino acid (minus CO₂ and the aldehyde Cα).
    - [ ] **8.0.5 Fix/verify remaining templates:** Close any remaining atom gaps with explicit H₂O, NH₃, or CO₂ byproduct species.
    - [ ] **Verification:** All `assert_balanced` tests pass with zero atom discrepancies across all 14 template families.

- [ ] **8.B — Calibrate FAST-mode Heuristic Barriers** `[Diff: 6/10 | Highest Scientific Impact]`
    - **Why:** The entire ranking system (Inverse Design scores, Lysine Budget, Trapping Efficiency) rests on constants like `Schiff=15 kcal`, `Strecker=22 kcal`. These are plausible guesses, but a 2–3 kcal shift changes which formulation "wins." Anchoring to real data is the single most impactful improvement.
    - [ ] 8.B.1: Run `xtb_screener.py` on 4 key reactions: Schiff base formation (glucose+glycine), Amadori rearrangement, Strecker decarboxylation (glycine), Cysteine β-elimination → H₂S.
    - [ ] 8.B.2: Compare computed xTB barriers against published DFT or experimental data (see `architecture.md` §5 refs).
    - [ ] 8.B.3: Update heuristic constants in `run_pipeline.py` and `inverse_design.py` to match the computed scale.
    - [ ] 8.B.4: Document the calibration in `docs/xtb_limitations.md`.
    - [ ] Verification: Ranked output for Ribose+Cysteine correctly places FFT above pyrazines, consistent with literature (pentose+cys > hexose+neutral-aa for sulfur heterocycles).

- [ ] **8.C — Add Concentration/Ratio Support to Formulation Grid** `[Diff: 5/10]`
    - **Why:** Scientists don't just choose ingredients — they choose ratios. The literature (`Maillard_Plant_based.md` §GSH) specifies optimal GSH conditions at pH 4.79, 131°C, 95 min. The tool currently can't distinguish a high-cysteine from a low-cysteine formulation.
    - [ ] Add optional `molar_ratios` field to entries in `formulation_grid.yml` (e.g., `{ribose: 1.0, cysteine: 0.5}`).
    - [ ] Update `InverseDesigner.evaluate_all()` to pass relative concentrations as a weight on barrier scores (higher-concentration reactant → lower effective barrier).
    - [ ] Update CLI to accept `--ratios cysteine:0.5,ribose:1.0` for forward-mode runs.
    - [ ] Verification: A formulation with 2× cysteine should score higher on `meaty` targets than the same formulation with 0.5× cysteine.

- [ ] **8.D — Model Lipid-Maillard Synergy (Pathway D Catalysis)** `[Diff: 7/10]`
    - **Why:** `pathways.md §D` explicitly describes lipid aldehydes (hexanal, nonanal) *catalyzing* Strecker degradation and reacting with H₂S to form alkylthiazoles — a "species-specific" flavor lever. Currently the tool only models hexanal as a nuisance to be trapped, missing this positive role entirely.
    - [ ] Add a `Lipid_Strecker_Catalysis` template in `smirks_engine.py`: lipid aldehyde + α-aminoketone → alkylthiazole / pyrazine (synergy product).
    - [ ] Add entries for key alkylthiazoles (e.g., 4-methyl-5-thiazoleethanol, benzothiazole) to `desirable_targets.yml` and the `meaty` tag in `sensory_tags.yml`.
    - [ ] Update the heme catalyst heuristic to also reduce barriers for the new synergy family.
    - [ ] Verification: A formulation with both lipids and cysteine/thiamine should yield alkylthiazoles not seen in a lipid-free or cysteine-free run.

---

## ✅ COMPLETED PHASES

## Phase 7: Plant-Based Formulation Tools & PBMA Usability ✅

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

> **Stretch goal (not blocking):** Model lipid-Maillard catalytic synergy — lipid aldehydes accelerating Strecker degradation (pathways.md §D). Currently only masking/trapping is modelled. → Moved to Phase 8.D.

---

## Phase 6: Automated Pathway Generation (Tier 0 Replacement) ✅

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

---

## Phase 4: Integration & Precursor Recommendation Prototype ✅

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

---

## Phase 3: Tier 2 — DFT Refinement with Skala ⏳

Focus on the 6–8 chemically decisive reactions only:

- [x] **3.1** Write `src/skala_refiner.py`: PySCF + Skala XC functional template for TS optimisation with CPCM water solvent.
- [x] **3.2** Write `tests/test_skala_refiner.py`: run a fast, low-basis single-point energy calculation on a small molecule to ensure PySCF+Skala integration works.
- [ ] **3.3 Compute barriers for key bifurcations:** `[Diff: 9/10]`
    - [ ] 3.3a: Amadori rearrangement (Schiff base → 1-amino-1-deoxy-2-ketose)
    - [ ] 3.3b: 2,3-enolisation vs 1,2-enolisation bifurcation point
    - [ ] 3.3c: Strecker decarboxylation (α-dicarbonyl + amino acid)
    - [ ] 3.3d: Cysteine + ribose → FFT (via furfural + H₂S)
    - [ ] 3.3e: Ribose retro-aldol → 1,4-dideoxyosone → MFT
    - [ ] 3.3f: DHA β-elimination (Ser → dehydroalanine)
    - [ ] 3.3g: Off-flavour trapping: hexanal + amino acid → Schiff base
    - [ ] 3.3h: α-aminoketone dimerisation → pyrazine
- [ ] **3.4 Validate with IRC:** `[Diff: 6/10]` Confirm each TS connects to correct reactant/product via intrinsic reaction coordinate.
- [ ] **3.5 Comparison Study:** `[Diff: 3/10]` Compare Skala barriers against published DFT or CCSD(T) data where available.

---

## Phase 2B: Physics & Implementation Improvements ✅

> Identified during cross-referencing the codebase against `Maillard_meat.md`, `Maillard_Plant_based.md`, `pathways.md`, and `architecture.md`.

### 🔴 Critical — Barriers, Conditions, and Solvation

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

### 🟡 Moderate — Data Correctness and API Validation

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

### 🟢 Minor — Robustness and Accuracy

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

## Phase 2: Tier 1 — xTB Energy Screening Pipeline ✅

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

---

## Phase 1C: RMG Rule Validation & Debugging ✅

**Goal:** All 14 custom Maillard reaction families load and parse correctly via `RMGDatabase.load_kinetics()`.

### Key Lessons Learned

- `groups.py` requires: `template(reactants=[...])` + `recipe(actions=[...])` + `entry()` + `tree()`
- Recipe atom labels (`*1`, `*2`...) must map to starred atoms in the adjacency list entries
- Adjacency list describes pre-reaction state — do NOT include bonds the recipe will form
- All bonds must be reciprocal (atom A lists B → B must list A)
- Use `db.local_context` (not `{}`) when loading to avoid `ArrheniusEP` errors in rules.py
- **Valence Violations:** Fixed 5-valent carbon/nitrogen errors in Enolisation, Amadori, Heyns, and Retro_Aldol families.
- **Product Count:** RMG `apply_recipe` requires template product count to match fragment count.

### Progress (14 families total — ALL PASSING ✅)

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

---

## Phase 1: Target Compound Library & Domain Knowledge Encoding ✅

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

---

## Phase 0: Project Infrastructure ✅

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

---

## ⏸️ DEFERRED PHASES

## Phase 8.1/8.2: Quantitative Kinetic Modeling (DEFERRED)

> **Context:** Precise temporal yield predictions (e.g., "how much FFT at t=12m") require robust DFT barriers. This is scientifically valuable but deferred until the core categorical recommendation tool (Phase 7) is usable.

- [ ] **8.1 Tier 2 DFT Refinement (from Phase 3.3):** Compute exact barriers using PySCF/Skala for rate-limiting steps. xTB barriers are too noisy for strict kinetic solvers.
- [ ] **8.2 Cantera Kinetic Simulation:**
    - Write a Cantera YAML mechanism generator that converts DFT barriers to Arrhenius rate constants.
    - Run microkinetic time-temperature profiles to plot concentration-vs-time curves for target volatiles.

## Phase 5: Experimental Validation Preparation (DEFERRED / OUT OF SCOPE)

> **Context:** As a purely computational framework, physical wet-lab validation is currently deferred until the complete *in silico* pipeline (including automated generation and kinetic filtering) is fully operational.

- [ ] **5.1** Select 5–10 top-ranked novel formulations from computational output.
- [ ] **5.2** Generate a lab protocol summary.
- [ ] **5.3** Define GC-MS validation criteria (expected retention times from NIST).
- [ ] **5.4** Document comparison methodology for predicted vs. observed agreement.
