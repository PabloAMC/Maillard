# Maillard Strategic Roadmap

## Focus

Build the most useful computational tool for scientists designing meat-like Maillard outcomes in plant-based matrices without overstating weak evidence.

## Current State

- [x] Free-precursor chemistry is quantitatively strong.
- [x] Matrix evidence, blocker classification, and scientist reporting surfaces are operational.
- [x] S9 authority infrastructure is executable and wired into the scientific lane.
- [x] P6 family scope is explicit, ranked, and now has a single next-family decision artifact.
- [x] A bounded runtime-executable mycoprotein reference is landed.
- [x] The smallest external evidence package for pea and soy promotion is explicit.
- [x] One primary external matrix package is now selected and specified for the next closure move.
- [x] Hexanal and Nonanal now have an explicit internal-resolution artifact that does not overclaim promotion.
- [x] Scope-gap families are explicitly guarded out of active expansion.
- [x] The full scientific Docker lane is green again after fixing the benchmark-target generator import and the refinement-campaign script wrapper.
- [x] Lipid-Maillard crosstalk is now benchmark-visible in the reaction-governance surface through the trapping_hexanal authority fixture bridge.
- [x] DHA/Lysinoalanine competition is now benchmark-visible in the reaction-governance surface through the dha_elimination authority fixture bridge.
- [x] Thiamine and glutathione support families are now explicitly benchmark-closed in the reaction-governance surface instead of relying on generic sulfur-family aliases.
- [x] Extrusion confidence/reporting now exposes a DHA/LAL closure panel that distinguishes direct closure markers from competition-only inference.
- [ ] No mixed matrix lane is externally decision-ready yet.
- [ ] Mixed matrix promotion still lacks externally measured closure.
- [ ] Selective refinement still has no benchmark-visible win.

## Task: Verify README images

- [x] Check `README.md` for duplicate images
- [x] Verify image paths in `README.md` exist and load correctly

## Task: Restructure README

- [x] Audit `README.md` for redundant sections (Accuracy/Validation/Regimes)
- [x] Consolidate "How to Use" and "Key Commands" into a logical flow
- [x] Improve "Quick Start" for new users
- [x] Ensure scientific posture (honesty about gaps) is preserved but clearly bucketed
- [x] Group validation artifacts and diagnostic commands to reduce early-document noise

## Active Priorities

### Matrix Closure

- [ ] Close one pea or soy mixed lane enough to materially improve external credibility.
- [x] Keep adverse-marker work evidence-first and calibration-first unless a new audit proves mechanistic leverage.
- [x] Make the smallest external evidence package explicit for the next matrix promotion step.

### S9 Authority Lane

- [x] Quasi-harmonic helper coverage is deterministic and executable.
- [x] Phase 3.3 and 3.5 mounted benchmark fixtures execute in the authority lane.
- [x] IRC proxy fixtures exercise the stable API contract.
- [x] The standard scientific lane runs S9 authority checks directly.
- [x] The generated skip registry currently reports zero skips.

### P6 Matrix Expansion

- [x] Matrix-family coverage distinguishes active lanes, bounded next candidates, and scope gaps.
- [x] Matrix-family ranking orders scope by impact and closability.
- [x] A single next-family decision is now encoded in the reporting surface.
- [x] Mycoprotein is the bounded next-family candidate for this cycle.
- [x] Coconut-oil and other plant-protein families remain explicit deferred scope gaps.

## Immediate Next Steps

1. [x] Land one externally meaningful matrix evidence package for pea or soy.
2. [x] Reduce ambiguity around Hexanal and Nonanal without inflating promotion posture.
3. [x] Keep mycoprotein bounded to calibration-grade expansion until a runtime-executable reference is landed.
4. [x] Hold scope-gap families out of active expansion until family-specific evidence exists.

## Success Criteria

- [x] A scientist can see exactly what evidence would change the next matrix decision.
- [x] Reports clearly separate internal calibration support from external promotion-ready support.
- [x] Scope does not expand faster than the evidence surface.

---

## Strategic Review: 26 March 2026

### Q1 — README Accessibility for Non-Computational Scientists

The README accessibility audit is now closed. The main scientist-facing onboarding gaps identified in the review have been addressed:

- [x] **Add a "What is the Maillard reaction?" primer** before the main positioning section.
- [x] **Add a "What does this tool output?" concrete example** so a scientist can see the report shape immediately.
- [x] **Add a minimal glossary** for FAST, SMIRKS, Strecker degradation, Amadori product, headspace partitioning, DFT, MLP, and adjacent terms.
- [x] **Make the Scientist Workflow concrete** with worked commands and explicit decision flow.
- [x] **Add a "Limitations at a glance" table** so trust boundaries are visible without reading the full deep dive.

### Q2 — Chemistry Coverage Gaps (Tool Completeness)

The amino acid + sugar core is strong. Re-audit of the codebase shows several previously suspected gaps are already implemented at some level, but not all of them are benchmark-closed or prominently surfaced. The remaining priority gaps are:

- [ ] **Lipid-Maillard crosstalk as a full kinetic lane**: hexanal/nonanal trapping (Schiff base formation with amino acids) is partially modelled but the aldehyde-catalysed Strecker acceleration and alkylthiazole/alkylpyrazine formation from lipid-Maillard adducts are not benchmark-closed lanes. This is the next highest-impact family sprint per `chemistry_family_scope_2026_03_20.md`.
- [x] **Lipid-Maillard crosstalk benchmark visibility**: the reaction benchmark surface now includes `trapping_hexanal -> lipid_strecker_synergy`, so family 02 is no longer invisible to DFT/MLP governance and coverage reporting.
- [x] **DHA/Lysinoalanine benchmark visibility**: the reaction benchmark surface now includes `dha_elimination -> beta_elimination`, so the lysine-loss competition lane is no longer invisible to DFT/MLP governance and coverage reporting.
- [x] **DHA/Lysinoalanine extrusion closure surface**: extrusion-facing confidence and reporting now show lysine-budget severity plus whether direct DHA/LAL closure markers are present, so scientist-facing outputs no longer hide this lane behind a generic extrusion warning.
- [ ] **DHA/Lysinoalanine external extrusion closure**: DHA and lysinoalanine reaction templates, competition metrics, benchmark visibility, and an explicit extrusion closure surface now exist, but the lane still lacks an extrusion-conditioned external benchmark package for promotion-grade claims.
- [x] **Thiamine thermal degradation closure**: the reaction benchmark surface now includes `thiamine_fragmentation -> thiamine_degradation`, so the family is promoted as an explicit benchmark-visible FAST lane for sulfur aroma prediction and extrusion-retention governance.
- [x] **Glutathione/cysteinyl glycine closure**: the reaction benchmark surface now includes `glutathione_release -> glutathione_cleavage`, so additive cleavage and downstream sulfur release are benchmark-closed for scientist-facing governance.
- [ ] **Heyns pathway closure**: Heyns support is present in priors, but fructose-dominant validation remains materially weaker than the Amadori/aldose trunk.
- [x] **Time-series prediction exposure**: the Docker workflow now exposes concentration-versus-time and temperature-profile generation through `kinetics-timeseries`.
- [x] **Wet-lab data ingestion loop**: the Docker workflow already exposes experiment comparison and benchmark materialization through `compare-experiment` and `materialize-experiment-benchmark`.
- [x] **Sensitivity analysis command**: the Docker workflow now surfaces family sensitivity generation as a primary scientist-facing command.

### Q3 — Realistic Path to DFT/MLP Coverage of All Families

Conclusions from repo memory, architecture.md, and P4 governance:

**What is realistic:**
- ~80% of important flavour families can be covered using curated literature Arrhenius parameters through the SLR ingestion workflow. The machine-readable SLR contract is already in place.
- Selective xTB screening can provide first-pass barrier estimates for families where no literature Ea exists. This is cheap enough to run for every new family before committing to DFT.
- DFT (at the B3LYP/6-31G* or ωB97X-D3/def2-TZVP level) is realistic for 3–5 decisive transition states per new family. NOT for exhaustive coverage of all pathways.
- MLP (mace_mp_small) is validated for geometry pre-optimisation only (mean RMSD ≈ 0.165 Å). It is not a barrier surrogate.

**What is NOT realistic without external resources:**
- Full DFT coverage of all Maillard families in a single sprint. The conformational search + TS optimisation cost is prohibitive across the full reaction graph.
- Using MLPs as barrier predictors — mace_off_medium is quarantined as a barrier surrogate per P4 governance.
- Replacing literature-calibrated Arrhenius priors with ab initio barriers everywhere — the kinetic priors are more calibrated to experimental concentrations than raw DFT Ea values would be.

**Concrete tiered plan:**
- [ ] **Tier 0 — Literature first**: for each new family sprint, run the SLR ingestion workflow to extract all available Arrhenius parameters before touching xTB or DFT.
- [ ] **Tier 1 — xTB screening gate**: if a key pathway step has no literature Ea, run xTB at GFN2-xTB level as a cheap screening gate. Only escalate to DFT if the xTB barrier is within the sensitive range (30–100 kJ/mol) or clearly non-physical.
- [ ] **Tier 2 — Selective DFT**: commission DFT only for steps that (a) have no literature anchor, (b) failed the xTB plausibility gate, and (c) are benchmark-visible — i.e., the step directly governs a compound with a measured literature concentration. Priority candidates from current P4 analysis: AIMNet2 or OrbMol for TS initialisation on lipid-Maillard crosstalk and DHA crosslinking TS.
- [ ] **Tier 3 — MLP geometry assist**: use mace_mp_small for conformer pre-optimisation to reduce DFT wall time. Do not use for energy or barrier ranking.
- [x] **Document the coverage map**: a machine-readable `data/lit/dft_coverage_map.json` generator now exists, along with scientist-facing validation artifacts in `results/validation/dft_coverage_map.{md,json}`.

## Root-Cause Next Plan — 26 March 2026

- [x] **Unify the external extrusion blocker as a contract artifact**: `data/protocols/extrusion_external_closure_contract.json` now encodes direct DHA/LAL, reactive-lysine/furosine, and promotion-grade extrusion requirements in one machine-readable place.
- [x] **Wire the contract into objective progress and reporting**: the new artifact now feeds `results/validation/extrusion_external_closure.{md,json}`, `objective_progress`, and the scientific reporting surface so the repo states which direct markers are missing versus merely contextual.
- [x] **Backfill literature-closable anchors into that contract**: `process_state_calibrations.json` now includes the De Leyn extrusion thiamine-retention anchor, and the contract explicitly reuses the Troise soy lysine-damage context without inflating either one into direct closure.
- [x] **Only then choose the next chemistry-family sprint**: the selected next family sprint is now **DHA/Lysinoalanine external package work**, because it is the highest-impact lane still blocked mainly by missing direct extrusion marker packages rather than by internal visibility or governance gaps.
- [x] **Specify the DHA/Lysinoalanine external package itself**: `data/protocols/dha_lysinoalanine_external_package_contract.json` and `results/validation/dha_lysinoalanine_external_package.{md,json}` now state the exact direct-damage markers, paired meaty-positive panel, extrusion metadata, and remaining post-package requirement without pretending that closure already exists.

## Dual External Package Follow-Through — 26 March 2026

- [x] **Materialize the mixed meaty-positive package for both pea and soy**: a dedicated dual-package artifact now names the exact external benchmark bundles for both matrices instead of only surfacing a single selected primary lane.
- [x] **Keep the single-lane primary package as a prioritization view**: the current selected-matrix artifact is preserved, but its policy now states that it does not replace the dual pea/soy package.
- [x] **Wire the dual package into objective progress and reporting**: the new artifact is exposed alongside the DHA/LAL package so both open external-data blockers are visible as explicit measurement bundles.
- [x] **Regenerate validation artifacts and rerun directed Docker tests**: the dual-package reporting surface and objective-progress deltas now render correctly, and the directed Docker suite passed.
