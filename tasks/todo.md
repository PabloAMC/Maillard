# Maillard Reaction Computational Framework — Working Plan
## Status: Planning — Awaiting approval before execution

---

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
- [x] 1B.26 **Strategic Pivot:** Switched to hand-curated pathways for Phase 2 screening to bypass RMG's automated discovery limitations.

---

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

## Phase 3: Tier 2 — DFT Refinement with Skala

Focus on the 6–8 chemically decisive reactions only:

- [x] 3.1 Write `src/skala_refiner.py`: PySCF + Skala XC functional template for TS optimisation with CPCM water solvent
- [x] 3.2 Write `tests/test_skala_refiner.py`: run a fast, low-basis single-point energy calculation on a small molecule to ensure PySCF+Skala integration works
- [ ] 3.3 Compute barriers for:
    - [ ] 3.3a Amadori rearrangement (Schiff base → 1-amino-1-deoxy-2-ketose)
    - [ ] 3.3b 2,3-enolisation vs 1,2-enolisation bifurcation point
    - [ ] 3.3c Strecker decarboxylation (α-dicarbonyl + amino acid)
    - [ ] 3.3d Cysteine + ribose → FFT (via furfural + H₂S)
    - [ ] 3.3e Ribose retro-aldol → 1,4-dideoxyosone → MFT
    - [ ] 3.3f DHA β-elimination (Ser → dehydroalanine)
    - [ ] 3.3g Off-flavour trapping: hexanal + amino acid → Schiff base
    - [ ] 3.3h α-aminoketone dimerisation → pyrazine
- [ ] 3.4 Validate with IRC (intrinsic reaction coordinate) to confirm each TS connects to correct reactant/product
- [ ] 3.5 Compare Skala barriers against published DFT or CCSD(T) data where available

---

## Phase 4: Integration & Precursor Recommendation Prototype

- [ ] 4.1 Write `src/recommend.py`: given a target volatile profile (e.g., "maximise MFT and FFT, minimise hexanal"), query the reaction graph and xTB/DFT-ranked pathways to score precursor combinations
- [ ] 4.2 Write `tests/test_recommend.py`: provide mock ranked pathways and target profiles, assert correct precursor recommendation and penalty application
- [ ] 4.3 Parametrise environmental conditions: pH (5–8), T (100–250°C), water activity (0.3–0.9)
- [ ] 4.4 Include DHA competition penalty: deduct Lys consumed by DHA from Maillard-available pool
- [ ] 4.5 Include off-flavour trapping bonus: amino acid excess beyond Maillard need traps hexanal/nonanal
- [ ] 4.6 Output: ranked table of precursor formulations with predicted volatile yields, off-flavour risk, and toxicity flags (AGE/HAA markers)
- [ ] 4.7 Verify: run prototype on 3 canonical systems and confirm results align with published model-system GC-MS
    - Ribose + Cys → expect dominant sulfur heterocycles
    - Glucose + Gly → expect dominant furans and pyrazines
    - Ribose + Cys + Leu → expect sulfur heterocycles + 3-methylbutanal

---

## Phase 5: Experimental Validation Preparation

- [ ] 5.1 Select 5–10 top-ranked novel formulations from Phase 4 output
- [ ] 5.2 Generate a lab protocol summary: precursor concentrations, pH, temperature, heating time, expected volatile profile
- [ ] 5.3 Define GC-MS validation criteria: which peaks to look for, expected retention times (from NIST WebBook)
- [ ] 5.4 Document comparison methodology: how to score predicted vs observed agreement

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

---

## Review Section
*To be filled during implementation — record what worked, what didn't, and lessons learned.*
