# Dev Log: 2026-03-05 — SmirksEngine Mass Conservation Audit

## Context
As the project moves into Phase 8.C (Calibrate FAST-mode Barriers), the physical validity of the `SmirksEngine` output became a critical bottleneck. Tier 1 screening (xTB) requires atom-balanced elementary steps to provide meaningful $\Delta E^\ddagger$ and $\Delta E_{rxn}$ values. 

## Key Decisions

### 1. Enforced Mass Conservation Policy
We have adopted a strict mass conservation policy. Every reaction template in `smirks_engine.py` must now strictly conserve atoms. 
- Discrepancies identified in the audit (dropped hydroxyls, missing $H_2O$, missing $NH_3$, or unaccounted $H_2$) have been resolved.
- A new utility `assert_balanced(step)` in `tests/test_smirks_engine.py` now enforces this at the unit-test level.

### 2. Hybrid Modeling Architecture (SMARTS + Handcrafted)
We have refined the engine to use a hybrid approach:
- **Tier A (SMARTS Rules)**: Used for simple, low-order stoichiometry (e.g., Schiff Base formation via `Lipid_Schiff_Base`).
- **Tier B (Handcrafted Functions)**: Used for complex rearrangements (Amadori, Strecker) and any reaction involving 3+ reactants (Thiol Addition, Thiazole Condensation).
- **Rationale**: Generic SMARTS matching loops for $N > 2$ reactants are combinatorially expensive and prone to false-positives (misidentifying organic thiols as $H_2S$, etc.). Handcrafted logic with explicit SMILES guards for small species ($H_2$, $H_2S$, $NH_3$) ensures scientific accuracy.

### 3. Canonicalizing $H_2$ as a Reductant
The formation of several aromatic volatiles (pyrazines, thiazoles, furfurylthiol) from their precursors is an oxidative process that releases hydrogen. To maintain mass balance without modeling complex radical redox chains, we have canonicalized molecular hydrogen ($H_2$) as an explicit byproduct or reactant depending on the direction.

## Technical Implementation Details
- **RDKit over Strings**: Replaced f-string/regex SMILES manipulation in the Schiff base and Amadori templates with `AllChem.ReactionFromSmarts` to preserve atom mapping and molecular fragments.
- **Dynamic Partitioning**: The Strecker step now dynamically partitions atoms from the dicarbonyl reactant into the resulting products using a fallback RDKit reaction if the dicarbonyl is unmapped.
- **Engine Hardening**: The `_apply_smirks` loop now discards any `ElementaryStep` where one or more products fail RDKit sanitization or validation, preventing "atom leaks" from invalid SMARTS applications.

## Status
- Core Maillard templates are 100% atom-balanced.
- Full integration suite passes (79 tests).
- Ready for Phase 8.C (xTB-based barrier calibration).
