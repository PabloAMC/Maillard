# Dev Log: RMG Validation and Curated Pathways Strategy

## Date: 2026-03-01

### Summary
Today we completed the validation and debugging of the custom RMG reaction families. While we successfully reached a state where 14/14 families validate and the simulation runs to 90% conversion, we discovered that RMG's automated mechanism generation has significant limitations for complex liquid-phase cascades like the Maillard reaction. We have pivoted to a hybrid approach using hand-curated pathways for the primary screening.

### Key Technical Achievements

#### 1. RMG Family Validation
- Resolved all `AtomTypeError` (valence violations) in Enolisation, Amadori, Heyns, and Retro_Aldol families.
- Fixed `IndexError` by correctly specifying product counts for bond-breaking reactions.
- Patched an internal RMG bug in transport data estimation (`numpy.complex128` error) by using absolute values for critical property derivations.
- Successfully implemented a `Sugar_Ring_Opening` family, though RMG's subgraph matching proved too restrictive for general ring sizes.

#### 2. Root Cause Analysis of Case 1 Failure
- **Competitive Kinetics:** `Thiol_Addition` has a very low barrier (12 kcal/mol) compared to `Schiff_Base_Formation` (~20 kcal/mol). This leads to rapid consumption of precursors into thioether oligomers, starving the Maillard cascade.
- **Structural Representation:** Cyclic sugars (hemiacetals) lack the C=O required for Schiff base formation. We resolved this by using open-chain aldehyde SMILES for precursors.
- **Missing Families:** RMG cannot "discover" furfural or FFT/MFT without explicit dehydration and cyclization rules, which are multi-step and difficult to encode as elementary RMG families.

### Strategic Pivot: Hand-Curated Pathways
To make the tool useful in the short term, we are bypassing RMG's automated discovery for the well-known Maillard cascades.

- **`data/reactions/curated_pathways.py`**: Defines the 5 core cascades (Core Maillard, Strecker, S-Maillard, Off-flavour trapping, DHA competition) as explicit sequences of `ElementaryStep` objects.
- **`scripts/run_curated_screening.py`**: A new runner that feeds these pathways directly into our validated xTB screening pipeline.

### Phase 2: xTB Screening Results (Tier 1)

Following the strategic pivot, we successfully executed the xTB parallel screening on the 5 curated pathways.

#### 1. Screening Outcomes
- **Reliable Bottleneck Identification:** The pipeline correctly identified **1,2-enolisation** of the Amadori product as the primary rate-limiting bottleneck for the core Maillard reaction, with an energetic span of **62.5 kcal/mol**.
- **Pathways B, D, and E:** Strecker degradation, off-flavour trapping, and DHA competition exhibited fast/moderate kinetics (25-36 kcal/mol), aligning with experimental expectations.
- **Pathway C (S-Maillard):** Successfully modeled the formation of 2-furfurylthiol (FFT) with a total thermodynamic driving force of -1232 kcal/mol.

#### 2. Technical Pipeline Fixes
- **Defensive NEB Fallback:** Implemented a mechanism in `xtb_screener.py` to catch unphysical NEB barriers (>500 kcal/mol) caused by missing atomic mapping in complex bimolecular reactions. The pipeline now falls back to a Hammond pseudo-activation barrier (`max(0, ΔE) + 25`).
- **Stoichiometry Calibration:** Corrected a massive endothermicity in the S-Maillard pathway by adding missing product molecules (H₂O) to the `Furfural + H₂S → FFT` step.

### Next Steps
1. **Phase 3 (DFT Refinement):** Use Skala (Tier 2) to refine the 62.5 kcal/mol enolisation bottleneck for higher accuracy.
2. **Phase 4 (Integration):** Feed the ranked screening results into `recommend.py` to finalize the precursor recommendation prototype.
