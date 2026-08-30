# QM Barrier Provenance — Computational-Gap Targets (Families 11-15 + 13-safety)

_Generated: 2026-04-21T12:21:11.943468+00:00_  

## Scope
7 computational-gap targets (Families 11-15) + acrylamide safety target — literature-anchored barrier provenance.

## Global decision
Scientific judgment on 2026-04-21: xTB/DFT values obtained for computational-gap targets were never accurate enough to carry as runtime anchors. Retired the xTB-derived numbers for hexanal_radical_quench and quinone_cys_michael; marked aa_ring_open_dicarbonyl as stack_limit_reached after two failed Sella-DFT recovery strategies; registered asparagine Ea from Knol 2009 directly. Diffusion-based TS generators (React-OT, OA-ReactDiff) could not be installed due to unresolved Python/torch-geometric/yarp dependency conflicts in maillard_validation; wrappers and tests were retired. Further refinement escalated to wet-lab backlog (Vía C).

## Per-target barrier table

| Target | F | Ea (kcal/mol) | Ea (kJ/mol) | Tier | Ceiling | Status | Literature |
|---|---|---|---|---|---|---|---|
| `hexanal_radical_quench` | 11 | — | — | no_literature_anchor | wet_lab_required | wet_lab_pending | _none_ |
| `quinone_cys_michael` | 13 | 6.93 | 29.0 | literature_family_surrogate | bounded_calibration | literature_anchor_only | Thiol Michael addition family, ResearchGate acrylamide-removal review, ~29 kJ/mol |
| `aa_ring_open_dicarbonyl` | 14 | 7.58 | 31.7 | literature_derived_transfer | bounded_calibration | stack_limit_reached | HCW Family 14 surrogate (no DOI in repo) |
| `pe_schiff_base` | 15 | 22.21 | 92.9 | literature_derived_transfer | bounded_calibration | literature_derived_transfer | PMC4419266 — phospholipid-amine Maillard, Ea=92.9 kJ/mol |
| `pe_amadori` | 15 | 19.81 | 82.9 | literature_derived_transfer | bounded_calibration | literature_derived_transfer | PMC4419266; PMC5992167 — Ea=82.9 kJ/mol consensus |
| `lysinoalanine_crosslink` | 12 | 16.00 | 67.0 | family_rule_surrogate | ranking_only | family_rule_surrogate | DHA family ArrheniusEP rule, E0=16 kcal/mol |
| `asparagine_sugar_explicit_water_cluster` | 13_safety | 30.83 | 129.0 | literature_derived_transfer | — | literature_derived_transfer | Knol et al. 2009 Food Chem, DOI 10.1016/j.foodchem.2009.11.049 — Ea=129 kJ/mol; Parker 2012 elimination leg Ea~105 kJ/mol |

## Detailed entries

### `hexanal_radical_quench` (SLR Family 11)
- **Label:** off-note lipid radical quench (hexanal H-abstraction)
- **Mechanism:** H-abstraction by phenoxyl/thiyl radical from hexanal
- **Active barrier:** _no anchor_  (runtime falls back to literature_runtime default)
- **Current tier → Target tier:** no_literature_anchor → wet_lab_anchor
- **Promotion ceiling:** wet_lab_required
- **Computational method:** No reliable anchor available. xTB-derived value (17.88 kcal/mol) retired 2026-04-21 because (i) GFN2 H-abstraction barriers for radical quench reactions are systematically biased and (ii) no direct literature Ea exists in repo for hexanal-quench by phenoxyl/thiyl radicals.
- **Active Arrhenius key:** `hexanal_radical_quench_no_anchor`
- **Honest label:** Hexanal radical quench has no reliable kinetic anchor. Runtime falls back to literature_runtime barrier_kj_mol=31.72 default; treat any prediction qualitatively only.
- **Refinability status:** wet_lab_pending
- **Refinability note:** Requires either: (a) hexanal H-abstraction Ea from radical-trapping literature ingestion, or (b) wet-lab measurement of hexanal suppression vs cysteine/sulfur partner concentration in oxidized-oil pre-cooking. See data/lit/deep_research_backlog.json for ResearchGate Ref 11 (lipid off-notes) ingestion.
- **Evidence paths:**
  - `data/geometries/dft_checkpoints/hexanal_radical_quench`
  - `results/computational_gap_refinement/hexanal_radical_quench_dft_execution.json`
  - `results/computational_gap_refinement/hexanal_radical_quench_xtb_execution.json`

### `quinone_cys_michael` (SLR Family 13)
- **Label:** polyphenol o-quinone + Cys Michael addition
- **Mechanism:** thiol Michael addition
- **Active barrier:** 6.93 kcal/mol (29.0 kJ/mol)  ±15.0 kJ/mol
- **Current tier → Target tier:** literature_family_surrogate → literature_derived_transfer
- **Promotion ceiling:** bounded_calibration
- **Computational method:** Literature family surrogate: thiol Michael addition to alpha,beta-unsaturated carbonyl (acrylamide-removal review, ResearchGate, ~29 kJ/mol). xTB-derived value was unreliable and has been retired.
- **Active Arrhenius key:** `quinone_cys_michael_thiol_addition_family`
- **Honest label:** Family surrogate (Michael thiol addition, ~29 kJ/mol literature). Same mechanism class as o-quinone + Cys; bounded calibration, +/- factor 3-8.
- **Refinability status:** literature_anchor_only
- **Refinability note:** Direct ingestion of 13_polyphenol_amino_capping deep-research source (currently in deep_research_backlog) would close this gap with a primary measurement.
- **Literature refs:**
  - Thiol Michael addition family, ResearchGate acrylamide-removal review, ~29 kJ/mol
- **Evidence paths:**
  - `data/geometries/dft_checkpoints/quinone_cys_michael`
  - `results/computational_gap_refinement/quinone_cys_michael_dft_execution.json`
  - `results/computational_gap_refinement/quinone_cys_michael_xtb_execution.json`

### `aa_ring_open_dicarbonyl` (SLR Family 14)
- **Label:** ascorbic-acid ring opening (stealth browning)
- **Mechanism:** ring-opening dicarbonyl source
- **Active barrier:** 7.58 kcal/mol (31.7 kJ/mol)  ±20.0 kJ/mol
- **Current tier → Target tier:** literature_derived_transfer → selective_dft_anchor
- **Promotion ceiling:** bounded_calibration
- **Computational method:** Family 14 HCW literature surrogate (bounded upstream dicarbonyl source)
- **Active Arrhenius key:** `aa_ring_open_dicarbonyl_hcw_family14_surrogate`
- **Honest label:** HCW Family 14 surrogate for an upstream ascorbic-acid dicarbonyl source, DFT validation pending, ±factor 2-5; bounded calibration
- **Refinability status:** stack_limit_reached
- **Refinability note:** Family 14 selective-DFT closure attempted with two strategies (interpolation+xTB constrained relax, xTB-CI-NEB+TBLite). Sella DFT verdict: minimum on both seeds; mode-match score 0.78 / |omega|=193 cm-1 below acceptance gate. Stack limit reached for concerted multi-coordinate TS at GFN2-xTB seed quality. Surrogate retained as bounded_calibration; do not re-attempt without substantially improved seed (e.g., diffusion-based geometry generator with public CHON checkpoint that supports this composition).
- **Literature refs:**
  - HCW Family 14 surrogate (no DOI in repo)
- **Evidence paths:**
  - `data/geometries/dft_checkpoints/aa_ring_open_dicarbonyl`
  - `results/computational_gap_refinement/aa_ring_open_dicarbonyl_dft_execution.json`
  - `results/computational_gap_refinement/aa_ring_open_dicarbonyl_interp_recovery.json`
  - `results/computational_gap_refinement/aa_ring_open_dicarbonyl_ts_dft_validation.json`
  - `results/computational_gap_refinement/aa_ring_open_dicarbonyl_xtb_execution.json`

### `pe_schiff_base` (SLR Family 15)
- **Label:** phosphatidylethanolamine Schiff base
- **Mechanism:** PE + reducing sugar Schiff base
- **Active barrier:** 22.21 kcal/mol (92.9 kJ/mol)  ±20.9 kJ/mol
- **Current tier → Target tier:** literature_derived_transfer → selective_dft_anchor
- **Promotion ceiling:** bounded_calibration
- **Computational method:** Literature transfer (SLR 15, phospholipid-amine Maillard)
- **Active Arrhenius key:** `pe_schiff_base_lit_derived`
- **Honest label:** SLR 15 literature Ea=92.9 kJ/mol, ±factor 2-5; bounded calibration
- **Refinability status:** literature_derived_transfer
- **Literature refs:**
  - PMC4419266 — phospholipid-amine Maillard, Ea=92.9 kJ/mol
- **Evidence paths:**
  - `data/geometries/dft_checkpoints/pe_schiff_base`
  - `results/computational_gap_refinement/pe_schiff_base_dft_execution.json`
  - `results/computational_gap_refinement/pe_schiff_base_xtb_execution.json`

### `pe_amadori` (SLR Family 15)
- **Label:** phosphatidylethanolamine Amadori
- **Mechanism:** PE-Schiff → PE-Amadori
- **Active barrier:** 19.81 kcal/mol (82.9 kJ/mol)  ±20.9 kJ/mol
- **Current tier → Target tier:** literature_derived_transfer → selective_dft_anchor
- **Promotion ceiling:** bounded_calibration
- **Computational method:** Literature transfer (SLR 15, phospholipid-amine Maillard)
- **Active Arrhenius key:** `pe_amadori_lit_derived`
- **Honest label:** SLR 15 literature Ea=82.9 kJ/mol, ±factor 2-5; bounded calibration
- **Refinability status:** literature_derived_transfer
- **Literature refs:**
  - PMC4419266
  - PMC5992167 — Ea=82.9 kJ/mol consensus
- **Evidence paths:**
  - `data/geometries/dft_checkpoints/pe_amadori`
  - `results/computational_gap_refinement/pe_amadori_dft_execution.json`
  - `results/computational_gap_refinement/pe_amadori_water_xtb_execution.json`
  - `results/computational_gap_refinement/pe_amadori_xtb_execution.json`

### `lysinoalanine_crosslink` (SLR Family 12)
- **Label:** DHA-Lys crosslink (safety)
- **Mechanism:** DHA + Lys Michael addition
- **Active barrier:** 16.00 kcal/mol (67.0 kJ/mol)  ±42.0 kJ/mol
- **Current tier → Target tier:** family_rule_surrogate → selective_dft_anchor
- **Promotion ceiling:** ranking_only
- **Computational method:** DHA_Crosslinking family-rule surrogate (ArrheniusEP E0=16 kcal/mol)
- **Active Arrhenius key:** `lysinoalanine_crosslink_dha_family_surrogate`
- **Honest label:** DHA-plus-lysine family surrogate, assumes prior dehydroalanine formation, ±factor 5-15; ranking only
- **Refinability status:** family_rule_surrogate
- **Literature refs:**
  - DHA family ArrheniusEP rule, E0=16 kcal/mol
- **Evidence paths:**
  - `data/geometries/dft_checkpoints/lysinoalanine_crosslink`
  - `results/computational_gap_refinement/lysinoalanine_crosslink_dft_execution.json`
  - `results/computational_gap_refinement/lysinoalanine_crosslink_xtb_execution.json`

### `asparagine_sugar_explicit_water_cluster` (SLR Family 13_safety)
- **Label:** acrylamide pathway (Asn + reducing sugar)
- **Mechanism:** lumped Asn + reducing-sugar -> acrylamide
- **Active barrier:** 30.83 kcal/mol (129.0 kJ/mol)  ±10.0 kJ/mol
- **Current tier → Target tier:** literature_derived_transfer → literature_derived_transfer
- **Promotion ceiling:** —
- **Computational method:** —
- **Active Arrhenius key:** `—`
- **Honest label:** —
- **Refinability status:** literature_derived_transfer
- **Literature refs:**
  - Knol et al. 2009 Food Chem, DOI 10.1016/j.foodchem.2009.11.049 — Ea=129 kJ/mol
  - Parker 2012 elimination leg Ea~105 kJ/mol
- **Evidence paths:**
  - `results/computational_gap_refinement/asparagine_sugar_explicit_water_cluster_dft_execution.json`
  - `results/computational_gap_refinement/asparagine_sugar_explicit_water_cluster_xtb_execution.json`

## Decisions captured (2026-04-21)
1. **hexanal_radical_quench**: xTB-derived 17.88 kcal/mol retired. No reliable literature anchor exists in repo. Runtime falls back to literature_runtime barrier_kj_mol=31.72 default. Wet-lab measurement pending (Vía C).
2. **quinone_cys_michael**: xTB-derived 17.84 kcal/mol retired. Replaced with thiol-Michael addition family surrogate ~29 kJ/mol (6.93 kcal/mol). Same mechanism class; bounded calibration, ±factor 3-8.
3. **aa_ring_open_dicarbonyl**: Marked stack_limit_reached. Sella+GFN2-xTB seed quality insufficient for this concerted TS. Surrogate 7.58 kcal/mol retained as bounded_calibration.
4. **asparagine_sugar_explicit_water_cluster**: Uses Knol 2009 literature Ea=30.83 kcal/mol (129 kJ/mol) through src/safety.py. DFT explicit-water cluster work retired; not needed.
5. **pe_schiff_base / pe_amadori / lysinoalanine_crosslink**: unchanged — already literature-derived-transfer / family-rule surrogates.
6. **Diffusion-based TS generators (React-OT, OA-ReactDiff)**: installation blocked in maillard_validation by upstream dependency conflicts; scaffolding retired.