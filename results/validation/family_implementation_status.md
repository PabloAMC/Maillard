# How the 16 families are actually implemented

Derived by enumeration from the current engine, not asserted. Regenerate with `python scripts/generators/generate_family_implementation_status.py`.

**5 of 16** SLR families have generative reaction templates; the rest are shared limbs, matrix/modifier layers, or literature priors.

| # | Family | Implementation | Reaction templates emitted |
| -- | --- | --- | --- |
| 01 | Amino acid-sugar Maillard core | generative reaction templates | `Amadori_Rearrangement`, `Aminoketone_Condensation`, `Cysteine_Degradation`, `Deoxyosone_Reduction`, `Enolisation_1_2`, `Enolisation_2_3`, `Enolisation_2_3_Amadori`, `Enolisation_Intermediate`, `Furanone_Cyclisation`, `Furanone_Formation`, `Generalized_Deamination`, `Retro_Aldol_Fragmentation`, `Schiff_Base_Formation`, `Strecker_Degradation`, `Thiohemiacetal_Formation`, `Thiol_Addition_H2`, `Thiol_Addition_Hexose_Legacy_Shortcut`, `Thiol_Addition_Pentodiulose`, `Thiol_Dehydration`, `Thiol_Oxidation` |
| 02 | Lipid oxidation and carbonylic crosstalk | generative reaction templates | `Beta_Scission`, `Lipid_Homolysis`, `Peroxy_H_Abstraction`, `Radical_Propagation_O2`, `Radical_Termination` |
| 03 | Thiamine degradation and sulfur support | generative reaction templates | `Additive_Thermal_Degradation`, `Furan_Ring_Aromatisation` |
| 04 | Nucleotide degradation and ribose support | literature priors only | — |
| 05 | Glutathione and peptide support | literature priors only | — |
| 06 | Alternative protein matrix scope | matrix / modifier layer, not reaction chemistry | — |
| 07 | Reducing sugar and carbonyl donor hierarchy | no templates of its own; modifier layer | — |
| 08 | Plant off-notes and Maillard suppression | no templates of its own; modifier layer | — |
| 09 | Carbohydrate pyrolysis and caramelization | no templates of its own (uses another family's) | — |
| 10 | Microbial fermentation pretreatment | matrix / modifier layer, not reaction chemistry | — |
| 11 | Maillard/Lipid Crosstalk | generative reaction templates | `Lipid_Schiff_Base`, `Lipid_Strecker_Synergy`, `Lipid_Thiazole_Condensation` |
| 12 | Protein Damage Markers | generative reaction templates | `Beta_Elimination`, `DHA_Crosslinking`, `Safety_Risk_AGE`, `Safety_Risk_Acrylamide` |
| 13 | Polyphenol-Amino Capping | literature priors only | — |
| 14 | Ascorbic Acid Maillard | curated layer only | — |
| 15 | PE Stealth Sugar Sink | literature priors only | — |
| 16 | Melanoidin Polymerization | matrix / modifier layer, not reaction chemistry | — |

## Implementation detail

* **01 Amino acid-sugar Maillard core** — 20 reaction template(s) emitted by src/reaction_templates.py / src/smirks_engine.py
* **02 Lipid oxidation and carbonylic crosstalk** — 5 reaction template(s) emitted by src/reaction_templates.py / src/smirks_engine.py
* **03 Thiamine degradation and sulfur support** — 2 reaction template(s) emitted by src/reaction_templates.py / src/smirks_engine.py
* **04 Nucleotide degradation and ribose support** — Ribose itself is a first-class sugar donor of families 01/07, but there is no nucleotide template: IMP/GMP enumerate to nothing. The lane exists as literature priors and coverage bookkeeping only (src/family_barrier_progress.py, src/dft_coverage_map.py).
* **05 Glutathione and peptide support** — No GSH template. Glutathione enumerates only through the generic amine/carbonyl templates of family 01 -- the `glutathione_cleavage` FAST_BARRIERS entry is never reached because no family emits it. Priors and bookkeeping only (src/literature_learning_loop.py, src/family_barrier_progress.py).
* **06 Alternative protein matrix scope** — Not reaction chemistry. Implemented as the matrix layer: src/matrix_correction.py (accessibility + volatile retention per protein type) and data/lit/protein_source_registry.json.
* **07 Reducing sugar and carbonyl donor hierarchy** — No templates of its own: it modulates family 01's. Implemented as `DONOR_REACTIVITY_MULTIPLIERS` in src/barrier_constants.py (per-family pentose/hexose/fructose/phosphorylated multipliers) plus `infer_carbohydrate_donor_identity`.
* **08 Plant off-notes and Maillard suppression** — No templates of its own. The off-notes are produced by family 02 and trapped by family 11 (`Lipid_Schiff_Base`); the suppression side is the intervention layer (src/pre_processor.py) and the matrix retention profiles.
* **09 Carbohydrate pyrolysis and caramelization** — No templates of its own. Furfural and HMF are produced by family 01's 1,2-enolisation limb (`Enolisation_1_2`); severity is reported by the projection layer's process-state index (src/projection.py).
* **10 Microbial fermentation pretreatment** — Not reaction chemistry. Implemented as an upstream precursor-enrichment modifier: src/pre_processor.py `KNOWN_INTERVENTIONS`, wired through src/literature_runtime.py.
* **11 Maillard/Lipid Crosstalk** — 3 reaction template(s) emitted by src/reaction_templates.py / src/smirks_engine.py
* **12 Protein Damage Markers** — 4 reaction template(s) emitted by src/reaction_templates.py / src/smirks_engine.py
* **13 Polyphenol-Amino Capping** — No template. The quinone/cysteine Michael step has a computational prior (`quinone_cys_michael`) which Wave G1 explicitly PARKED because the curated layer has no step to route it to. src/literature_runtime.py only.
* **14 Ascorbic Acid Maillard** — No template. Present only as an ascorbate prior attached to the curated `Enolisation_1_2` step (src/curated_pathways.py, src/literature_runtime.py). [P] that prior's provenance is still open -- see tasks/audit_remediation.md.
* **15 PE Stealth Sugar Sink** — No template; the PE headgroup enumerates to nothing. Implemented as a sugar-depletion bookkeeping term from the `pe_schiff_base` / `pe_amadori` priors in src/literature_runtime.py.
* **16 Melanoidin Polymerization** — Not reaction chemistry. Implemented as the thiol-trapping factor `_MELANOIDIN_TRAPPING_PROFILES` / `_resolve_melanoidin_trapping_factor` in src/recommend.py, gated on the family-16 lane being active.

## Entry-point findings

* The lipid radical chain runs only from a HYDROPEROXIDE: the `unoxidised_lipid_plus_o2` case enumerates to 0 steps. There is no initiation step from an unoxidised fatty acid, so in production the chain is seeded by the lipid-oxidation anchor rather than by the network itself.
* The thiamine cascade matches thiamine by EXACT canonical SMILES (Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1) and by a >= 100 C gate; a differently written thiamine SMILES silently produces nothing.
* Glutathione, ascorbate, catechol, the phospholipid headgroup and IMP reach no template of their own -- they enumerate only through the generic family-01 amine/carbonyl templates, or to zero steps.
