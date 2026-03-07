# Development Log: Transition State Architecture Decision
**Date:** 2026-03-07
**Context:** Phase 3.3 (Tier 2 DFT Refinement of Key Bifurcations)

## Problem Statement
While implementing the Tier 2 DFT Refinement pipeline (`src/dft_refiner.py`), we encountered a bottleneck in generating the initial 3D geometries for the Transition States (TS) of the 8 target Maillard reactions. The sheer expense and fragility of *ab-initio* TS scans (e.g., using String method or NEB) led to considering an alternative architecture borrowed from the **AltZyme** project.

## The AltZyme "Core Template" Proposal
AltZyme accelerates enzyme mutant evaluations by taking a known, pre-calculated True Transition State for the Wild Type (the "core template"), swapping the surrounding amino acid sidechains (mutagenesis), applying fractional bond restraints to freeze the reaction center, and running classical constrained relaxation to resolve steric clashes. 

We briefly considered adapting this for the Maillard framework: generating a generic TS core for each reaction family (e.g., Amadori rearrangement) and simply grafting different sugar and amino acid R-groups onto it.

## Why it Fails for Maillard Chemistry
The proposal was rejected after careful analysis for the following chemical and physical reasons:

1. **Fundamental Substituent Variability:** In enzymology, the actual reacting substrate is usually identical, only the non-covalent protein scaffold changes. In Maillard chemistry, the **reactants themselves** are highly variable. A Ribose (pentose) Schiff base has a fundamentally different steric bulk, hydrogen bonding network, and optimal folding path compared to a Glucose (hexose) Schiff base. 
2. **Backbone Conformational Shifts:** The Amadori arrangement relies on proton transfers and specific dihedral angles that are highly sensitive to the identity of the amino acid. A Cysteine (polar, bulky, sulfur-containing) vs. Glycine (no sidechain) will dramatically alter the true lowest-energy TS conformation. 
3. **Artificial Strain Generation:** If we enforce a strict 3D constraint on the core reacting atoms (derived from a Glycine+Ribose template) and graft on larger substituents (like Leucine+Glucose), the constrained relaxation will introduce massive, artificial strain energies into the surrounding bonds. The subsequent `wB97M-V` single point calculation would incorrectly report an astronomically high barrier merely as an artifact of the frozen core constraint.

## Selected Architecture
To maintain scientific credibility and ensure that our computed barriers accurately reflect the differences between specific formulations (e.g., matching experimental literature that shows Pentose is more reactive than Hexose), we **must** perform unique structural minimizations for every sugar/AA pairing.

### The Pipeline Workflow:
1. **Tier 1 Initial Guess (xTB):** We will use semi-empirical methods (`xtb`) or specialized reaction path generators (`RMG-Py`) to compute the true minimum-energy TS geometry for the exact pair in question.
2. **Tier 2 Refinement (PySCF):** The resulting `.xyz` structure from Tier 1 is fed directly into `scripts/run_tier2_dft.py` to undergo high-level `r2SCAN-3c` optimization and `wB97M-V` electronic refinement.

We are establishing a firm rule: **No generic 3D TS templates will be used in the Maillard pipeline.** Every variant requires an explicitly computed Tier 1 structural minimum.
