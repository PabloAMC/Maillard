# SMIRKS Template Provenance & Audit

This document clarifies the origin, chemical justification, and validation status of the reaction templates used in `src/smirks_engine.py`.

## 1. Governance Model
The Maillard Framework uses a **Hybrid Generator**:
- **Tier A (SMIRKS)**: Generic functional group transforms (Schiff base, Thiol addition).
- **Tier B (Templates)**: Complex rearrangements where atom-mapping is non-trivial (Strecker, Amadori).

## 2. Rule Registry

### Tier A: Generic Transforms
| Rule Name | Reaction Family | Chemical Basis | Provenance |
| :--- | :--- | :--- | :--- |
| `schiff_base_lipid` | Lipid_Schiff_Base | Nucleophilic attack of R-NH2 on lipid aldehydes | *Martins et al. (2001)* |
| `beta_scission_alkoxy` | Beta_Scission | Radical fragmentation of alkoxy lipids | *Frankel (1985)* |
| `radical_propagation_o2` | Radical_Propagation | O2 trap of alkyl radicals | *Frankel (1985)* |

### Tier B: Curated Templates
These are implemented as Python functions to handle specific stoichiometry and H2O/CO2 release.

| Template | Mechanism | Key Intermediate | Validation Status |
| :--- | :--- | :--- | :--- |
| `_amadori_cascade` | Glycation | N-substituted-1-amino-1-deoxy-2-ketose | Validated vs. Glucose/Glycine benchmarks |
| `_strecker_step` | Oxidative Decarboxylation | Strecker Aldehyde + CO2 | Validated vs. Methionine/Thiamine systems |
| `_beta_elimination_steps` | DHA Pathway | Dehydroalanine (DHA) | Functional but requires rheology link |
| `_thiazole_condensation` | Heterocycle synthesis | 2-alkylthiazoles | Heuristic (Yaylayan 1990) |

## 3. Physical Constraints (The "Black Box" Guard)
To prevent the "combinatorial explosion" common in unsupervised rule engines, the following filters are applied in `smirks_engine.py`:
1. **MW Capping**: Volatile products must be < 300 Da.
2. **Canonical Deduplication**: Every intermediate is canonicalized via RDKit before recruitment.
3. **Atomic Balance**: All templates are strictly atom-balanced, including H2O and CO2 co-products.

## 4. Current Audit Gaps
- **Nitrate/Nitrite Interventions**: Not yet modeled.
- **Metal Catalysis**: SMIRKS rules for iron coordination are present as placeholders.
