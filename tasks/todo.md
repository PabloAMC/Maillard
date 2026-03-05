# Maillard Reaction Computational Framework — Working Plan
## Status: Planning — Awaiting approval before execution

---

## ⚡ Prioritized Action Plan (as of 2026-03-04)

Ordered by **impact on alt-protein formulation scientists** (the end-user):

| Priority | Task | Status | Why |
|----------|------|:------:|-----|
| ~~1~~ | 7.1 — Wire SmirksEngine → Recommender | ✅ | Core value: combinatorial screening of novel formulations. |
| ~~2~~ | 7.2 — PBMA Precursors & Degradation Rules | ✅ | Thiamine, Glutathione, Lipid trapping — the PBMA-specific chemistry. |
| ~~3~~ | 7.3 — Advanced CLI Formulation Interface | ✅ | Scientists can run realistic formulations without editing Python. |
| ~~4~~ | 7.4 — Sensory & Trapping Output Metrics | ✅ | Translates SMILES into "Meaty", "Beany" — the language formulators speak. |
| **5** | 7.6 — Lysine Budget & DHA Competition | ☐ | Without this, the tool **overestimates** aroma yields (pathways.md §E). |
| **6** | 7.7 — Water Activity (`--aw`) Flag | ☐ | Critical processing variable for extrusion/spray-drying; trivial to add. |
| **7** | 7.5 — Inverse Design Mode | ☐ | Answers the scientist's real question: "What do I add to get meaty?" |
| 8 | 3.3 — DFT barriers (Skala) | ☐ | Replaces noisy heuristics with accurate kinetics. |
| 9 | 8.2 — Cantera Kinetics | ☐ | Time/temperature cooking simulations (needs DFT barriers first). |

> **Deferred:** Phase 5 (Experimental validation), Phase 8 (Kinetic modeling), AI-models (IBM RXN).

---

## 🏃‍♂️ ACTIVE & UPCOMING PHASES

## Phase 7: Plant-Based Formulation Tools & PBMA Usability (NEW CORE FOCUS)

> **Context:** The `SmirksEngine` is now wired to the Recommender, but the tool currently acts as a *generic* Maillard simulator. To be useful for alt-protein scientists, it must specifically address PBMA formulation challenges highlighted in the literature: managing lipid off-flavors ("beany" notes), utilizing complex additives (thiamine/glutathione), and simulating transition-metal (heme) catalysis.

- [x] 7.1 **Wire SmirksEngine to Recommender:** Create a unified pipeline script (`scripts/run_pipeline.py`) that accepts simple precursors, runs `SmirksEngine.enumerate()`, screens energetic bottlenecks, and feeds results directly into `recommend.py`.
- [x] 7.2 **Implement Core PBMA Precursors & Degradation Rules:** `[Diff: 7/10 | Model: Claude Opus / Gemini Pro]` Extrude the `SmirksEngine` to handle the heavy lifters of PBMA formulation:
    - **Add Lipids:** [x] Add `Hexanal` and `Nonanal` to `data/species/precursors.yml` to simulate the inherent oxidizing matrix of pea/soy isolates.
    - **Add Thiamine Template:** [x] Write a Tier B `_thiamine_degradation` function in `smirks_engine.py` that cleaves Thiamine into H₂S, 2-methylthiophene, and 4,5-dihydro-2-methylthiazole.
    - **Add Glutathione Template:** [x] Write a Tier B `_glutathione_cleavage` function in `smirks_engine.py` that cleaves GSH into Glutamic Acid and the highly reactive Cysteinylglycine dipeptide.
    - **Testing:** [x] Add assertions in `tests/test_smirks_engine.py` to ensure Thiamine and GSH break down correctly, and ensure the existing `Lipid_Schiff_Base` rule correctly consumes Hexanal + an amino acid.
- [x] 7.3 **Advanced CLI Formulation Interface:** [x] `[Diff: 4/10 | Model: Gemini Flash]` Implement complex `argparse` allowing users to run realistic formulations: e.g., `--sugars xylose --amino-acids cysteine --additives thiamine --lipids hexanal --catalyst heme --ph 6.0`.
    - **Update Resolver:** [x] Add `"lipids"` to the categories indexer in `src/precursor_resolver.py`.
    - **Update CLI Arguments:** [x] In `scripts/run_pipeline.py`, add `--additives`, `--lipids`, and `--catalyst` (currently supports `"heme"`).
    - **Implement Heme Heuristic:** [x] In `FAST` mode, if `--catalyst heme` is set, reduce barriers by 5.0 kcal for `Strecker_Degradation` and `Aminoketone_Condensation` families.
    - **Verification:** [x] Run `python scripts/run_pipeline.py --sugars xylose --amino-acids cysteine --catalyst heme` and verify Pyrazine/Strecker barriers are reduced compared to a control run.
- [x] 7.4 **Sensory & Trapping Output Metrics:** [x] `[Diff: 5/10 | Model: Claude Sonnet]`
    - **OAV & Sensory Lookup:** [x] Ensure all entries in `desirable_targets.yml` and `off_flavour_targets.yml` have `sensory_desc` and `odour_threshold_ug_per_kg`.
    - **Trapping Score Logic:** [x] Update `src/recommend.py` to calculate a "Lipid Trapping Potential" based on the fraction of lipid pathways that end in Schiff bases vs. remaining free.
    - **Output Enhancement:** [x] Update the CLI table to display sensory descriptors and the calculated Trapping Efficiency.
    - **Verification:** [x] Run `python scripts/run_pipeline.py --lipids hexanal --amino-acids lysine` and verify the "Trapping Efficiency: XX%" metric is displayed.
- [x] 7.5 **Inverse Design Mode:** [x] `[Diff: 8/10 | Model: Claude Opus]` Implement a search mode where the user specifies target profiles (e.g., `--target meaty --minimize beany`), and the tool evaluates a predefined grid of precursor formulations to recommend the optimal industrial blend.
    - **Define Sensory Tags:** [x] Create a mapping from target names to human-readable tags (`meaty`, `roasted`, `beany`, `grassy`) in `data/species/sensory_tags.yml`.
    - **Build Formulation Grid:** [x] Define a grid of realistic precursor combinations spanning pea/soy bases with various sugar+additive combos in `data/formulation_grid.yml`.
    - **Grid Evaluation Engine:** [x] Write `inverse_design.py` that runs the pipeline over each grid entry, collects the predicted targets, and scores each formulation against the user's `--target` / `--minimize` profile.
    - **Ranking & Output:** [x] Display a ranked table of top-N formulations sorted by match score, with columns: Formulation, Target Score, Risk Penalty, Lys Budget, Trap Eff.
    - **Verification:** [x] Run `python scripts/run_pipeline.py --target meaty --minimize beany` and confirm meaningful ranking.
- [x] 7.6 **Lysine Budget & DHA Competition Metric:** [x] `[Diff: 3/10 | Model: Gemini Flash]` Surface the stoichiometric competition between the DHA crosslinking pathway and the Maillard cascade to the user.
    - **Logic:** [x] In `recommend.py`, count how many generated steps consume Lysine (via DHA_Crosslinking) vs. Maillard (Schiff base), and compute a ratio.
    - **Output:** [x] Display a `⚠️ Lysine Budget: X% consumed by DHA crosslinking` warning in the PBMA metrics block when Lysine + Serine/Cysteine are in the pool.
    - **Verification:** [x] Run with `--amino-acids lysine,serine` and confirm the warning appears.
- [x] 7.7 **Water Activity CLI Flag:** [x] `[Diff: 2/10 | Model: Gemini Flash]` Expose the existing `ReactionConditions.water_activity` field via a `--aw` CLI argument in `run_pipeline.py`.
    - **Verification:** [x] Run with `--aw 0.5` and confirm conditions are displayed correctly.

> **Stretch goal (not blocking):** Model lipid-Maillard catalytic synergy — lipid aldehydes accelerating Strecker degradation (pathways.md §D). Currently only masking/trapping is modelled.

## Phase 3: Tier 2 — DFT Refinement with Skala

Focus on the 6–8 chemically decisive reactions only:

- [x] 3.1 Write `src/skala_refiner.py`: PySCF + Skala XC functional template for TS optimisation with CPCM water solvent
- [x] 3.2 Write `tests/test_skala_refiner.py`: run a fast, low-basis single-point energy calculation on a small molecule to ensure PySCF+Skala integration works
- [ ] 3.3 **Compute barriers for key bifurcations:** `[Diff: 9/10 | Model: Gemini Pro (Long Context for Logs)]`
    - [ ] 3.3a Amadori rearrangement (Schiff base → 1-amino-1-deoxy-2-ketose)
    - [ ] 3.3b 2,3-enolisation vs 1,2-enolisation bifurcation point
    - [ ] 3.3c Strecker decarboxylation (α-dicarbonyl + amino acid)
    - [ ] 3.3d Cysteine + ribose → FFT (via furfural + H₂S)
    - [ ] 3.3e Ribose retro-aldol → 1,4-dideoxyosone → MFT
    - [ ] 3.3f DHA β-elimination (Ser → dehydroalanine)
    - [ ] 3.3g Off-flavour trapping: hexanal + amino acid → Schiff base
    - [ ] 3.3h α-aminoketone dimerisation → pyrazine
- [ ] 3.4 **Validate with IRC:** `[Diff: 6/10 | Model: Claude Sonnet]` (intrinsic reaction coordinate) to confirm each TS connects to correct reactant/product
- [ ] 3.5 **Comparison Study:** `[Diff: 3/10 | Model: Gemini Flash]` Compare Skala barriers against published DFT or CCSD(T) data where available

---

## ⏸️ DEFERRED PHASES

## Phase 8: Quantitative Kinetic Modeling (DEFERRED)

> **Context:** Precise temporal yield predictions (e.g., "how much FFT at t=12m") require robust DFT barriers. This is scientifically valuable but deferred until the core categorical recommendation tool (Phase 7) is usable.

- [ ] 8.1 **Tier 2 DFT Refinement (from Phase 3.3):** Compute exact barriers using PySCF/Skala for rate-limiting steps. xTB barriers are too noisy for strict kinetic solvers.
- [ ] 8.2 **Cantera Kinetic Simulation:** 
    - Write a Cantera YAML mechanism generator that converts DFT barriers to Arrhenius rate constants.
    - Run microkinetic time-temperature profiles to plot concentration-vs-time curves for target volatiles.

## Phase 5: Experimental Validation Preparation (DEFERRED / OUT OF SCOPE)

> **Context:** As a purely computational framework, physical wet-lab validation is currently deferred until the complete *in silico* pipeline (including automated generation and kinetic filtering) is fully operational.

- [ ] 5.1 Select 5–10 top-ranked novel formulations from computational output
- [ ] 5.2 Generate a lab protocol summary
- [ ] 5.3 Define GC-MS validation criteria (expected retention times from NIST)
- [ ] 5.4 Document comparison methodology for predicted vs. observed agreement

---

## ✅ COMPLETED PHASES

## Phase 0: Project Infrastructure

- [x] 0.1 Initialise git repository and `.gitignore`
- [x] 0.2 Create directory structure: `src/`, `data/`, `data/species/`, `data/reactions/`, `results/`, `scripts/`
- [x] 0.3 Create conda/mamba environment spec (`environment.yml`) with:
    - RMG-Py and RMG-database
    - xtb (via conda-forge)
    - PySCF + Skala (via pip from `github.com/microsoft/skala`)
    - ASE (Atomic Simulation Environment)
    - RDKit (molecular handling)
    - NetworkX (reaction graph analysis)
- [x] 0.4 Write `scripts/check_env.py` — smoke-test that all packages import and key binaries execute
- [x] 0.5 Verify: `python scripts/check_env.py` passes on a clean install

---

## Phase 1: Target Compound Library & Domain Knowledge Encoding

### 1A. Build the target compound database

- [x] 1A.1 Create `data/species/desirable_targets.yml` with SMILES, InChI, CAS, odour threshold, and sensory descriptor for each:
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
- [x] 1A.2 Create `data/species/off_flavour_targets.yml`:
    - Hexanal — beany, grassy (from linoleic acid oxidation)
    - Nonanal — cardboard (from oleic acid oxidation)
    - 1-Octen-3-ol — mushroom, grassy
    - 2-Pentylfuran — beany
    - 1-Hexanol — green
- [x] 1A.3 Create `data/species/toxic_markers.yml`:
    - CML (Nε-carboxymethyllysine)
    - CEL (Nε-carboxyethyllysine)
    - Acrylamide
    - HMF (5-hydroxymethylfurfural)
    - Lysinoalanine (LAL) — DHA pathway marker
    - Lanthionine (LAN) — DHA pathway marker
- [x] 1A.4 Create `data/species/precursors.yml` for the input molecules:
    - Amino acids: L-cysteine, L-methionine, L-leucine, L-isoleucine, L-valine, L-phenylalanine, glycine, L-alanine, L-lysine, L-serine
    - Sugars: D-ribose, D-xylose, D-glucose, D-fructose
    - Exogenous precursors: thiamine, glutathione (GSH)
- [x] 1A.5 Verify: every entry has a valid SMILES that round-trips through RDKit

### 1B. Encode domain-specific reaction families for RMG

Each family below corresponds to a pathway identified in `pathways.md` (A–E) and the literature.

**Pathway A — Core Maillard Cascade:**
- [x] 1B.1 Schiff base formation: nucleophilic addition of amine to open-chain reducing sugar carbonyl → N-substituted glycosylamine
- [x] 1B.2 Amadori rearrangement: glycosylamine → 1-amino-1-deoxy-2-ketose (for aldose sugars)
- [x] 1B.3 Heyns rearrangement: glycosylamine → 2-amino-2-deoxyaldose (for ketose sugars)
- [x] 1B.4 1,2-enolisation (acidic pH): Amadori product → 3-deoxyosone → furfural/HMF
- [x] 1B.5 2,3-enolisation (neutral/alkaline pH): Amadori product → 1-deoxyosone → α-dicarbonyls (pyruvaldehyde, diacetyl)
- [x] 1B.6 Retro-aldol fragmentation: C5/C6 sugar skeleton → C2/C3 reactive fragments
- [x] 1B.7 Encode pH-conditional branching rule: if pH < 6 → favour 1B.4; if pH ≥ 6 → favour 1B.5 (Handled via ReactionConditions during screening)

**Pathway B — Strecker Degradation:**
- [x] 1B.8 α-dicarbonyl + free amino acid → Strecker aldehyde + α-aminoketone + CO₂
- [x] 1B.9 Strict Strecker mapping: Val→2-methylpropanal, Leu→3-methylbutanal, Ile→2-methylbutanal, Met→methional, Phe→phenylacetaldehyde
- [x] 1B.10 α-aminoketone self-condensation → dihydropyrazine → pyrazine (oxidation)

**Pathway C — S-Maillard (Sulfur Pathway):**
- [x] 1B.11 Cysteine thermal degradation → H₂S + NH₃ + pyruvaldehyde
- [x] 1B.12 Pentose (xylose/ribose) dehydration → furfural
- [x] 1B.13 Furfural + H₂S → 2-furfurylthiol (FFT)
- [x] 1B.14 Ribose retro-aldol → 1,4-dideoxyosone → 2-methyl-3-furanthiol (MFT)
- [x] 1B.15 Thiamine thermal fragmentation → H₂S + 2-methylthiophene
- [x] 1B.16 Glutathione cleavage → glutamic acid + cysteinylglycine dipeptide (then Cys enters 1B.11)

**Pathway D — Lipid-Maillard Crosstalk & Off-Flavour Trapping:**
- [x] 1B.17 Hexanal/nonanal + free amino acid → non-volatile Schiff base (trapping reaction)
- [x] 1B.18 Lipid aldehyde + H₂S → long-chain alkylthiazoles (synergy pathway)

**Pathway E — Dehydroalanine (DHA) Competing Pathway:**
- [x] 1B.19 Serine β-elimination → dehydroalanine (DHA)
- [x] 1B.20 Cysteine β-elimination → DHA + H₂S
- [x] 1B.21 DHA + Lysine ε-amino → lysinoalanine (LAL)
- [x] 1B.22 DHA + Cysteine thiol → lanthionine (LAN)
- [x] 1B.23 Encode stoichiometric competition: Lys consumed by DHA = Lys unavailable for Maillard (Handled via recommend.py)

- [x] 1B.24 Verify: run RMG with `{D-ribose, L-cysteine}` input → **Result:** Ribose reacts via Thiol_Addition, but competitive kinetics prevent FFT/MFT discovery.
- [x] 1B.25 Verify: run RMG with `{D-glucose, glycine}` input → **Result:** Amadori product requires better dehydration/cyclization families.
- [x] 1B.26 **Strategic Pivot (MVP):** Switched to hand-curated pathways for Phase 2/4 to get the energetic and recommendation engine working. 
- [ ] 1B.27 **Long-Term Requirement:** RMG is optimized for gas-phase combustion. We must identify a better liquid-phase generative engine (e.g., AutoMeKin, Chemoton, or AI-driven forward synthesis models) to fulfill the tool's true purpose of *automated discovery*. (See Phase 6).

---

## Phase 1C: RMG Rule Validation & Debugging

**Goal:** All 14 custom Maillard reaction families must load and parse correctly via `RMGDatabase.load_kinetics()`.

### Key lessons learned:
- `groups.py` needs: `template(reactants=[...])` + `recipe(actions=[...])` + `entry()` + `tree()`
- Recipe atom labels (`*1`, `*2`...) must map to starred atoms in the adjacency list entries
- Adjacency list describes pre-reaction state — do NOT include bonds the recipe will form
- All bonds must be reciprocal (atom A lists B → B must list A)
- Use `db.local_context` (not `{}`) when loading to avoid `ArrheniusEP` errors in rules.py
- **Valence Violations:** Fixed 5-valent carbon/nitrogen errors in Enolisation, Amadori, Heyns, and Retro_Aldol families by correcting proton transfers in recipes.
- **Product Count:** RMG `apply_recipe` requires template product count to match fragment count. Fixed with `products=["Product1", "Product2"]`.

### Progress (14 families total — ALL PASSING ✅):
- [x] 1C.1 **Thiol_Addition** — reference implementation ✅
- [x] 1C.2 **Amadori_Rearrangement** — 1,2-proton shift ✅
- [x] 1C.3 **Heyns_Rearrangement** — ketose equivalent ✅
- [x] 1C.4 **Strecker_Degradation** — fixed: product count 1→2, removed FORM+CHANGE combo ✅
- [x] 1C.5 **Retro_Aldol_Fragmentation** — fixed: simplified to BREAK_BOND only, product count 1→2 ✅
- [x] 1C.6 **Schiff_Base_Formation** — fixed: product count 1→2, removed FORM+CHANGE combo ✅
- [x] 1C.7 **Lipid_Schiff_Base** — fixed: same as Schiff_Base_Formation ✅
- [x] 1C.8 **DHA_Crosslinking** — fixed: product count 1→2 ✅
- [x] 1C.9 **Enolisation** — fixed: rules.py label mismatch ✅
- [x] 1C.10 **Aminoketone_Condensation** — fixed: simplified groups, product count 1→2 ✅
- [x] 1C.11 **Lipid_Thiazole_Condensation** — fixed: rules.py label mismatch ✅
- [x] 1C.12 **Beta_Elimination** — fixed: product count 1→2, rules.py Serine entry removed ✅
- [x] 1C.13 **Cysteine_Degradation** — fixed: product count 1→2 ✅
- [x] 1C.14 **Sugar_Ring_Opening** — New family for hemiacetal opening ✅

## Phase 2: Tier 1 — xTB Energy Screening Pipeline

- [x] 2.1 Set up test framework: create `tests/` directory, add `pytest` to `environment.yml`, and configure `pytest.ini`
- [x] 2.2 Write `src/pathway_extractor.py`: parse RMG output graph → extract list of elementary reaction steps (reactant, product, reaction family)
- [x] 2.3 Write `tests/test_pathway_extractor.py`: parse a known small RMG graph (e.g., Ribose+Cys output) and assert correct elementary steps are extracted
- [x] 2.4 Write `src/xtb_screener.py`: for each elementary step:
    - Generate 3D coordinates via RDKit ETKDG
    - Optimise reactant and product geometries with GFN2-xTB (GBSA water solvent)
    - Run relaxed coordinate scan or NEB for crude TS estimate
    - Record ΔErxn and ΔE‡
- [x] 2.5 Write `tests/test_xtb_screener.py`: test coordinate generation and a single fast optimization on a trivial molecule (e.g., water or ammonia)
- [x] 2.6 Parallelise via Python `multiprocessing` (xTB is ~10–60s per job)
- [x] 2.7 Write `src/pathway_ranker.py`: rank pathways by cumulative barrier height; flag rate-limiting steps
- [x] 2.8 Write `tests/test_pathway_ranker.py`: provide mock energetic data and assert correct ordering and rate-limiting step identification
- [x] 2.9 Verify: benchmark xTB ΔE‡ for Amadori rearrangement against any published computational data
- [x] 2.10 Verify: confirm relative ordering — ribose reactions should rank faster than glucose reactions (literature: ribose ~5× more reactive)
- [x] 2.11 Document known xTB limitations: proton transfers in water, ionic intermediates, β-elimination barriers may be off by 5–15 kcal/mol
- [x] 2.12 Write `data/reactions/curated_pathways.py`: Define 5 core Maillard cascades (A-E) as explicit `ElementaryStep` sequences to bypass RMG limitations
- [x] 2.13 Write `scripts/run_curated_screening.py`: Runner script that feeds curated pathways directly into the Phase 2 pipeline for screening
- [ ] 2.14 Analyze `curated_screening_results.json` and identify the fast vs slow pathway bottlenecks

---

## Phase 2B: Physics \u0026 Implementation Improvements (Post-Audit)

> Identified during cross-referencing the codebase against `Maillard_meat.md`, `Maillard_Plant_based.md`, `pathways.md`, and `architecture.md`.

### 🔴 Critical — Barriers, Conditions, and Solvation

- [x] 2B.1 **Replace fake barrier formula in `xtb_screener.py`:** the current `max(0, ΔE) + 15.0` constant is physically meaningless. Implement xTB `--path` (NEB) mode to scan the reaction coordinate between reactant and product complexes, yielding real approximate ΔE‡ values.
    - [x] 2B.1a Implement `_run_xtb_path()` method using `xtb --path` with reactant/product XYZ inputs
    - [x] 2B.1b Parse the NEB output to extract the highest-energy image as the barrier estimate
    - [x] 2B.1c Update `compute_reaction_energy()` to call NEB instead of the constant offset
    - [x] 2B.1d Write `tests/test_xtb_neb.py` with mocked NEB output to verify parsing logic

- [x] 2B.2 **Add pH / Temperature / Water Activity (aᵥ) parametrisation across the pipeline:**
    - [x] 2B.2a Create a `ReactionConditions` dataclass in `src/conditions.py` with fields: `pH`, `temperature_celsius`, `water_activity`
    - [x] 2B.2b In `pathway_extractor.py` or a new `src/condition_filter.py`: apply pH-conditional branching rules — if pH < 6, upweight 1,2-enolisation pathways; if pH ≥ 6, upweight 2,3-enolisation (per `Maillard_meat.md` §2.2 and `pathways.md` §A)
    - [x] 2B.2c Apply Arrhenius scaling: use literature activation energies (23–238 kJ/mol range) to modulate relative pathway rates with temperature
    - [x] 2B.2d Apply water-activity bell-curve: peak reactivity at aᵥ 0.6–0.8, suppressed at < 0.3 (glass transition) and > 0.9 (dilution), per `Maillard_meat.md` §6.2
    - [x] 2B.2e Write `tests/test_conditions.py`: assert that pH < 6 returns furan-dominant pathways and pH ≥ 7 returns pyrazine-dominant pathways

- [x] 2B.3 **Add explicit water molecules in transition states for proton-transfer steps:**
    - [x] 2B.3a For Amadori rearrangement and enolisation TS geometries, include 1–2 explicit H₂O in the XYZ input to model the Grotthuss proton-shuttle mechanism
    - [x] 2B.3b Document which steps require explicit water vs. implicit-only in `docs/xtb_limitations.md`
    - [x] 2B.3c Update `skala_refiner.py` to accept explicit-water TS geometries for Tier 2 refinement

### 🟡 Moderate — Data Correctness and API Validation

- [x] 2B.4 **Fix incorrect D-Ribose SMILES in `benchmark_xtb.py`:**
    - [x] 2B.4a Correct ribose open-chain SMILES from current 6-carbon aldohexose (`OC[C@H](O)[C@H](O)[C@@H](O)C=O`) to correct 5-carbon aldopentose (`OC[C@@H](O)[C@H](O)C=O`)
    - [x] 2B.4b Fix derived Schiff base SMILES accordingly
    - [x] 2B.4c Add a SMILES validation check in `scripts/check_env.py` that round-trips all precursor SMILES through RDKit and asserts correct heavy-atom counts (ribose = 5C, glucose = 6C)

- [x] 2B.5 **Validate Skala API against actual Microsoft documentation:**
    - [x] 2B.5a Read the Skala GitHub repo (`github.com/microsoft/skala`) README and example scripts
    - [x] 2B.5b Update `skala_refiner.py` `_build_mf()` to use the correct, documented API (currently speculative `skala.apply_skala(mf)`)
    - [x] 2B.5c Fix `run_tier2_dft.py` mock XYZ geometries to have correct atom counts and even electron numbers (current ones cause PySCF odd-electron errors)
    - [x] 2B.5d Re-run `tests/test_skala_refiner.py` after API correction

- [x] 2B.6 **Uncheck Phase 3 items that were only scaffolded, not executed:**
    - [x] 2B.6a Items 3.3a–3.3h, 3.4, 3.5 were marked complete but no actual DFT barriers were computed — revert to `[ ]` and re-execute once Skala API is validated and real TS geometries are available

### 🟢 Minor — Robustness and Accuracy

- [x] 2B.7 **Multi-conformer sampling in `xtb_screener.py`:**
    - [x] 2B.7a Generate N conformers (e.g., 5–10) via RDKit `EmbedMultipleConfs()` instead of single `EmbedMolecule(seed=42)`
    - [x] 2B.7b UFF-optimize all conformers, select the lowest-energy one as the xTB input
    - [x] 2B.7c Update `tests/test_xtb_screener.py` to verify multi-conformer logic

- [x] 2B.8 **Add physics-level integration tests (not just mocks):**
    - [x] 2B.8a Write `tests/test_integration_xtb.py` (marked `@pytest.mark.slow`): run actual xTB on water + formaldehyde → methanediol and assert ΔE is exothermic and within ±5 kcal/mol of literature
    - [x] 2B.8b Write `tests/test_integration_pyscf.py` (marked `@pytest.mark.slow`): run PySCF single-point on ethanol and assert energy is within expected range
    - [x] 2B.8c Require at least one integration test to pass before merging any PR that touches `xtb_screener` or `skala_refiner`

- [x] 2B.9 **Improve pathway ranking method in `pathway_ranker.py`:**
    - [x] 2B.9a Replace `max(barriers)` with **Energetic Span Analysis** (Kozuch \u0026 Shaik, 2011): rank by the energy difference between the TOF-determining TS and the TOF-determining intermediate
    - [x] 2B.9b For linear (non-branching) pathways, also report the cumulative sum of barriers as a secondary metric
    - [x] 2B.9c Update `tests/test_pathway_ranker.py` to validate the new ranking logic against known test cases

---

## Phase 4: Integration & Precursor Recommendation Prototype

> **Current status:** Phase 4 is complete. `recommend.py` loads curated pathways and screening results, runs the pathway ranker, and outputs a rich Unicode table.

**🔴 Highest Priority — Wire up real recommendations:**
- [x] 4.0 **Annotate curated pathways with target compounds** — each pathway in `curated_pathways.py` must declare what aroma compound it produces (e.g., `C_S_Maillard_FFT → produces: FFT`).
- [x] 4.1 **Rewrite `src/recommend.py`** to accept precursor inputs and conditions, load relevant curated pathways, run `PathwayRanker`, and output a ranked list of reachable target compounds. Remove all mocked logic.
- [x] 4.2 **Add multi-precursor comparison variants to `curated_pathways.py`:** at minimum:
    - Ribose + Cysteine (current)
    - Glucose + Glycine (Core Maillard, furans/pyrazines)
    - Ribose + Cysteine + Leucine (adds Strecker 3-methylbutanal)
    - Methionine path (methional / cooked-potato note)
- [x] 4.3 **Add penalty scoring for competing pathways:** DHA pathway consuming Lys, off-flavour risk (Pathway D).
- [x] 4.4 **Add toxicity flags:** if the scored network includes AGEs (CML, CEL) or HMF, flag at output.
- [x] 4.5 **Validate output on 3 canonical model systems** (per `architecture.md` §5 Phase 2):
    - Ribose + Cys → expect FFT dominant
    - Glucose + Gly → expect furans and pyrazines dominant
    - Ribose + Cys + Leu → expect FFT + 3-methylbutanal
- [x] 4.6 **Write `tests/test_recommend.py`** with real (non-mock) pathway data to confirm the output is physically reasonable.
- [x] 4.7 **Output format:** produce a human-readable ranked table: precursor → predicted volatiles → confidence (barrier rank) → off-flavour risk → toxicity
- [x] 4.8 **Document known limitations:** xTB barriers are relative rankings only; absolute yield predictions require Tier 2 DFT.

---

## Phase 6: Automated Pathway Generation (Tier 0 Replacement)

> **Context:** The current pipeline securely evaluates thermodynamics (xTB/Skala) but relies on a hardcoded list of curated pathways (`curated_pathways.py`). To truly explore the combinatorial ingredient space, we must build a hybrid generative pipeline. This phase is the **core intellectual contribution** of the project.

### 6.1 Deterministic Rule Enumeration (Hybrid SMIRKS + Templates)

> **Goal:** Replace the hand-curated `PATHWAYS` dict with a rule engine that automatically enumerates Maillard reaction pathways for *any* sugar + amino acid combination.
>
> **Approach — Hybrid Design:**
> - **Tier A (SMIRKS):** Simple functional-group transforms where atom mapping is unambiguous: Schiff base formation (aldehyde + amine → imine), sugar dehydration (pentose → furfural, hexose → HMF), thiol addition (furfural + H₂S → FFT). These are robust, well-defined SMIRKS.
> - **Tier B (Parameterised Templates):** Complex multi-step rearrangements where SMIRKS would be chemically misleading: Amadori/Heyns rearrangement, Strecker degradation, β-elimination (DHA). These are encoded as *template functions* that accept any sugar/amino acid pair and return the correct `ElementaryStep` sequence by SMILES substitution.
>
> **Safeguards against combinatorial explosion:**
> - Molecular weight cap of 300 Da on generated products.
> - Aldehyde-specific SMARTS (`[CH:1]=O`) instead of generic carbonyl matching.
> - Canonical SMILES deduplication at each generation.
> - Depth-limited to 3 generations (configurable).

- [x] 6.1.1 **Define Tier A SMIRKS rules** for 3 simple transforms:
    - Schiff base formation: aldehyde + primary amine → imine + H₂O
    - Pentose dehydration: pentose → furfural + 3 H₂O (pH < 6)
    - Thiol addition: furfural + H₂S → FFT + H₂O
- [x] 6.1.2 **Define Tier B template functions** for 4 complex rearrangements:
    - `amadori_template(sugar, amino_acid)` → `[SchiffBase, AmadoriProduct, Deoxyosone]` steps
    - `strecker_template(dicarbonyl, amino_acid)` → `[StreckerAldehyde, Aminoketone, CO₂]` step
    - `enolisation_template(amadori_product, pH)` → 1,2- or 2,3-enolisation products
    - `beta_elimination_template(cysteine_or_serine)` → DHA + side products
- [x] 6.1.3 **Write `src/smirks_engine.py`:** `SmirksEngine` class that:
    - Accepts `List[Species]` + `ReactionConditions`
    - Applies Tier B templates first (to generate core intermediates)
    - Then applies Tier A SMIRKS iteratively on the growing pool (with pruning)
    - Outputs `List[ElementaryStep]` compatible with `xtb_screener.py`
- [x] 6.1.4 **Add pH-conditional gating** via `conditions.py`: 1,2-enolisation template fires only at pH < 6; 2,3-enolisation only at pH ≥ 6.
- [x] 6.1.5 **Write `tests/test_smirks_engine.py`:**
    - Ribose + Glycine @ pH 5 → produces Schiff base + Amadori + furfural
    - Glucose + Glycine @ pH 7 → produces Schiff base + Amadori + pyruvaldehyde
    - Ribose + Cysteine → thiol addition produces FFT
    - All outputs are valid `ElementaryStep` objects
- [x] 6.1.6 **Verify against curated pathways:** Run on 4 canonical systems, confirm it rediscovers Pathways A–E from `curated_pathways.py`.

### 6.2 Extended Template Coverage (Remaining 6 Families)

> **Goal:** Extend `smirks_engine.py` to cover *all* 14 known Maillard reaction families, completing the deterministic reaction enumeration layer.
>
> **Why not AI tools?** IBM RXN and ASKCOS are trained on >90% pharmaceutical synthesis data. The missing Maillard reactions (aminoketone condensation, retro-aldol fragmentation, cysteine thermolysis) have essentially zero representation in their training sets. Extending the well-tested template approach is faster, more reliable, and chemically verifiable.

- [x] 6.2.1 **Aminoketone_Condensation template:** 2 × aminoketone → dihydropyrazine → pyrazine + H₂O. Highest impact — pyrazines are key roasted/meaty volatiles.
- [x] 6.2.2 **Retro_Aldol_Fragmentation template:** 3-deoxyosone → acetol + glycolaldehyde (C2 + C3 fragments). These are precursors for pyrazine formation.
- [x] 6.2.3 **Cysteine_Degradation template:** Cys (thermal, >150°C) → H₂S + acetaldehyde + NH₃. Provides the sulfur source independently of beta-elimination.
- [x] 6.2.4 **Lipid_Thiazole_Condensation template:** aldehyde + H₂S + NH₃ → 4,5-dihydrothiazole. Secondary savory aroma compounds.
- [x] 6.2.5 **Heyns_Rearrangement template:** Variant of Amadori for ketose sugars (fructose). Minor modification to existing template.
- [x] 6.2.6 **Sugar_Ring_Opening template:** Hemiacetal → open-chain aldehyde form (trivial prerequisite step).
- [x] 6.2.7 **Update `tests/test_smirks_engine.py`:** Add tests for pyrazine formation, retro-aldol products, cysteine thermolysis, and thiazole formation.
- [x] 6.2.8 **Verify 14/14 family coverage:** Run on all canonical systems and confirm every RMG family has a corresponding template or SMIRKS rule.

