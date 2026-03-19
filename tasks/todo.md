# Maillard SLR-To-Repository Execution Plan

## Mission

Rebuild the scientific backlog around the actual delta between:

- what the SLR now contains;
- what the repository has already encoded into machine-readable artifacts and production logic;
- what the user can actually see in reports and campaign outputs.

This plan is no longer a generic calibration roadmap. It is a step-by-step SLR incorporation program plus an architecture for the remaining gaps that literature cannot close.

## Non-Negotiable Rules

- [x] Docker plus conda remains the only authoritative validation path.
- [ ] Do not treat SLR prose as incorporated science until it is encoded in a machine-readable artifact or enforced in code.
- [ ] Do not treat encoded science as user-visible until it appears in report JSON, Markdown, and comparison outputs.
- [ ] Prefer literature anchor -> calibration payload -> computational prior -> cheap runtime surrogate -> offline semiempirical/QM -> selective DFT -> localized ML potential.
- [ ] Do not push DFT or ML potentials into the daily runtime path.
- [ ] Every expensive computation must write back a reusable correction artifact, not a one-off notebook result.
- [ ] Do not overclaim matrix quantitation where the SLR only provides directional or transferred support.

## What "Best Possible Tool" Means Here

### Scientific quality

- [ ] Every material claim in the recommendation surface is traceable to either a quantitative paper, a transferred prior, a censoring surrogate, or a computational artifact.
- [ ] Every key compound shown to a scientist carries an evidence tier and calibration source.
- [ ] Important PBMA-vs-meat flavor markers are not silently omitted just because the current model grew around sulfur alone.

### Product quality

- [ ] A scientist can see not only the predicted ranking, but which literature anchors actually drove it.
- [ ] Comparison reports reveal when a formulation wins for the wrong reason, such as furfural inflation or pyrazine overshoot.
- [ ] The repo exposes what is encoded, what is transferred, what is still missing, and what is being approximated computationally.

### Architecture quality

- [ ] Cheap-first calibration remains the default path.
- [ ] DFT is selective and benchmark-driven.
- [ ] ML potentials are only used where repeated motif classes justify a trained local surrogate.

## Audit Result — 2026-03-19

### Already encoded in the repo with real machine-readable artifacts

- [x] Squeo 2023 acrylamide endpoint reference in data/lit/safety_reference_payloads.json.
- [x] Asen 2022 pea denaturation anchor in data/lit/process_state_calibrations.json.
- [x] Li 2025 pea free-SH accessibility anchor in data/lit/process_state_calibrations.json.
- [x] Pratap-Singh 2021 ambient pea/soy observable calibration in src/matrix_calibration_registry.py.
- [x] Heated soy hexanal and censored 2-pentylfuran carryover from Shu 2024 in src/matrix_calibration_registry.py.
- [x] Process-state and matrix priors for pea and soy in data/lit/computational_priors.json.
- [x] Meaty-quality ratio and penalty in src/pipeline.py, src/bayesian_optimizer.py, and src/reporting.py.
- [x] Melanoidin sulfur trapping surrogate in src/recommend.py.
- [x] DFT triage contract in src/dft_refinement_contract.py.

### Present in the SLR but not yet encoded as machine-readable science

- [x] Foods 12(10):1967 PBMA acrylamide plus AGE reference band.
- [x] Choi 2024 acrylamide plus asparagine correlation payload.
- [x] MDPI Foods 11(21):3505 SPI free-SH anchor for native versus heated soy isolate.
- [x] JAFC 10.1021/acs.jafc.3c05991 PPI/SPI hexanal binding payload.
- [x] Xu 2023 SPI temporal off-flavor decline versus time at 95 C.
- [x] Karolkowski 2021 PPI pH-dependent off-flavor release payload.
- [x] Hernandez 2023 Strecker aldehyde sub-benchmarks.
- [x] Hernandez 2023 pyrazine and furfural sub-benchmarks as PBMA reference outputs.
- [x] Trikusuma 2019 pea-protein aroma interaction payload.
- [x] Troise 2018 soy Amadori and furosine calibration payload.
- [x] Lincoln 2025 SPI/PPI polyphenol-modulated Strecker/lipid cross-talk payload.
- [x] Laemont 2023 pyrazine pH/sugar control prior.
- [x] Arsa 2022 pyrazine pH optimum prior.
- [x] Hao 2025 sugar-type to pyrazine profile prior for soy hydrolysate.
- [x] Cerny 2007 thiamine-versus-xylose MFT contribution split.
- [x] Blank and Fay 1996 HEMF/DMHF mechanistic prior.
- [x] Watanabe 2015 DMHF meat-side reference anchor if full quantitative values can be recovered.
- [x] Zhang 2026 indirect sulfur volatile retention prior for SPI.

### Encoded partially but not modeled at the full level implied by the SLR

- [x] Ince 2024 is acknowledged directionally, but hexanal release is not yet modeled as a temperature-dependent reversible retention surface.
- [x] Hofmann 2001 is reflected by a cheap trapping surrogate, but not yet tied to an explicit browning index shown to the scientist.
- [ ] Nishimura 2024 is recognized as qualitative hydrolysate evidence, but its peptide-bound cysteine and protease-specific pyrazine findings are not encoded as structured priors.
- [x] Hernandez 2023 drives the MFT/furfural quality idea, but its Strecker and pyrazine subdata are not yet first-class output references.

### Encoded internally but not shown well enough to the user

- [x] Comparison reports now expose MFT/furfural ratio.
- [x] Comparison reports now expose meaty_quality_penalty.
- [x] Comparison reports now expose sulfur trapping behavior.
- [ ] Comparison reports still do not expose compact key projection-factor deltas beyond sulfur trapping.
- [x] Reports now show a compact evidence ladder per key compound.
- [x] Reports now distinguish direct anchor, transferred prior, mechanistic surrogate, and computational refinement at the compound row level.

### Structural gaps confirmed by the SLR and still unresolved by public data

- [ ] Quantitative aqueous PPI plus ribose plus cysteine benchmark with MFT, FFT, hexanal, and adverse markers.
- [ ] Quantitative aqueous SPI plus ribose plus cysteine benchmark with the same panel.
- [ ] Direct matrix-resolved MFT/FFT retention and trapping curves versus browning severity.
- [ ] Simultaneous Ellman plus OPA plus DSC benchmark in PPI/SPI thermal experiments.
- [ ] Aqueous SPI/PPI acrylamide kinetics across controlled T/t in wet matrices.
- [ ] Absolute pyrazine kinetics in PPI/SPI with internal standards.
- [ ] Absolute HEMF and DMHF in PPI/SPI ribose systems.

## Track 0: Build The SLR Incorporation Matrix

### Goal

- [x] Make it impossible to confuse cited literature with encoded calibration.

### New artifacts

- [x] data/lit/slr_incorporation_matrix.json
- [x] data/lit/flavor_reference_payloads.json
- [x] data/lit/retention_reference_payloads.json
- [x] data/lit/process_gap_registry.json

### Track 0 deliverables

- [x] One canonical incorporation ledger covering every SLR paper and every quantitative sub-payload.
- [x] One flavor-reference payload file for non-safety compounds and meat/PBMA anchor bands.
- [x] One retention-reference payload file for headspace and binding data.
- [x] One explicit structural-gap registry separated from codable literature.

### Track 0 schema contract

#### slr_incorporation_matrix.json

- [x] One object per paper or per extracted sub-payload when one paper supports multiple independent model surfaces.
- [x] Required fields:
  - [x] slr_section
  - [x] paper_id
  - [x] citation
  - [x] doi
  - [x] matrix_family
  - [x] compounds_supported
  - [x] parameters_supported
  - [x] exact_numeric_anchors
  - [x] current_repo_artifacts
  - [x] current_runtime_consumers
  - [x] current_user_visible_surfaces
  - [x] incorporation_status
  - [x] next_action
  - [x] confidence_tier
  - [x] notes_on_limits

#### flavor_reference_payloads.json

- [x] Organize by:
  - [x] sulfur_reference_anchors
  - [x] strecker_reference_anchors
  - [x] pyrazine_reference_anchors
  - [x] carbonyl_reference_anchors
  - [x] furanone_reference_anchors
- [x] Every entry must include:
  - [x] matrix_context
  - [x] analytical_method
  - [x] units
  - [x] benchmark_role
  - [x] target_direction
  - [x] numeric_band_or_point

#### retention_reference_payloads.json

- [x] Organize by compound class and then by matrix family.
- [x] Every entry must include:
  - [x] retention_or_release_mode
  - [x] direct_binding_or_headspace_measure
  - [x] temperature_context
  - [x] time_context
  - [x] reversibility_assumption
  - [x] transferability_notes

#### process_gap_registry.json

- [x] One object per unresolved structural gap.
- [x] Required fields:
  - [x] gap_id
  - [x] gap_type
  - [x] blocks_modules
  - [x] why_literature_cannot_close_it
  - [x] cheapest_next_step
  - [x] computational_fallback
  - [x] wet_lab_requirement

### Exact tasks

1. [x] Create one row per paper or quantitative sub-payload in the SLR.
2. [x] For each row, record:
  - [x] paper id
  - [x] compounds or parameters supported
  - [x] exact numeric anchors
  - [x] current repo location if already encoded
  - [x] current user-visible surface if already shown
  - [x] evidence status: encoded, partially encoded, not encoded
  - [x] action type: encode now, show now, compute later, requires wet lab
3. [x] Split compound-level flavor references away from safety references so reporting can consume them directly.
4. [x] Record structural gaps separately so they stop being mixed with codable literature.

### Track 0 execution order

1. [x] Populate the matrix first for all papers already present in benchmark_intake_registry.json, process_state_calibrations.json, computational_priors.json, and safety_reference_payloads.json.
2. [x] Add all SLR papers currently cited only in docs/slr_benchmark_evaluation.md.
3. [x] Promote papers with direct numeric anchors into flavor_reference_payloads.json or retention_reference_payloads.json.
4. [x] Leave non-codable papers in the matrix with explicit next_action values rather than inventing payloads.
5. [x] Validate that every Track 1 ingestion item is traceable back to one matrix row.

### Acceptance criteria

- [x] Every paper in docs/slr_benchmark_evaluation.md is mapped to a concrete repo status.
- [x] We can answer, from one JSON, whether a datum is only cited, encoded, modeled, and shown.
- [x] We can generate a machine-readable diff between SLR content and repo content without rereading prose.

## Track 1: Encode All SLR Data That We Already Have But Have Not Yet Ingested

### 1A. Safety references

#### Goal

- [x] Expand the safety surface beyond Squeo 2023 so the repo reflects the full SLR safety layer.

#### Files

- [x] data/lit/safety_reference_payloads.json
- [x] src/reporting.py
- [x] src/safety.py

#### Track 1A deliverables

- [x] Safety payloads cover industrial endpoints, PBMA endpoints, precursor correlation anchors, and free-system kinetics provenance.
- [x] Reports show which safety statements are endpoint references versus kinetic estimates.

#### Exact tasks

1. [x] Add Foods 12(10):1967 PBMA acrylamide range and AGE context.
2. [x] Add Choi 2024 acrylamide plus asparagine correlation payload.
3. [x] Add Knol 2009, Claeys 2005, and Sen and Gokmen 2023 as kinetic model provenance payloads.
4. [x] Add Zilic 2014 only as contextual heat-transfer modifier evidence, not as a benchmark.
5. [x] Add explicit payload tags for:
  - [x] industrial_endpoint_reference
  - [x] finished_product_reference
  - [x] precursor_correlation_reference
  - [x] kinetic_model_reference
  - [x] contextual_process_modifier
6. [x] Define which safety payloads should surface in default user reports and which should stay in extended provenance only.

#### Acceptance criteria

- [x] Reports can show both industrial endpoint bands and free-system kinetic provenance.
- [x] Safety output distinguishes endpoint references from kinetic models.
- [x] No safety reference is misrepresented as if it were a controlled benchmark kinetics dataset.

### 1B. Process-state and retention references

#### Goal

- [x] Fill the current soy and retention blind spots that the SLR already covers, except for the still-unextracted Trikusuma 2019 quantitative table.

#### Files

- [x] data/lit/process_state_calibrations.json
- [x] data/lit/retention_reference_payloads.json
- [x] data/lit/computational_priors.json
- [x] src/matrix_calibration_registry.py
- [x] src/headspace.py
- [x] src/matrix_correction.py

#### Track 1B deliverables

- [x] Soy process-state anchors stop depending only on repo synthesis where direct literature exists.
- [x] Retention/reference payloads exist for aldehydes and sulfur proxy volatiles.
- [x] Computational priors become explicit about which entries are direct anchor, transferred anchor, or surrogate interpolation.

#### Exact tasks

1. [x] Add MDPI Foods 11(21):3505 as a soy free-SH anchor with native and heated states.
2. [x] Add JAFC 3c05991 as a PPI/SPI hexanal binding payload.
3. [x] Add Karolkowski 2021 as a pH-dependent PPI release payload.
4. [x] Add Xu 2023 as a temporal SPI off-flavor trend payload.
5. [x] Add Zhang 2026 as an indirect sulfur-volatile retention prior for soy.
6. [x] Add Trikusuma 2019 as a pea mixed-pathway aroma interaction payload if the quantitative table is recoverable enough for intake encoding.
7. [x] Add Troise 2018 as a lysine-loss/furosine calibration payload for soy thermal history.
8. [x] Add explicit provenance_tier normalization across these payloads:
  - [x] direct_measurement
  - [x] literature_derived_transfer
  - [x] mechanistic_surrogate
  - [x] repo_literature_synthesis
  - [x] manual_fallback
9. [x] Define an explicit censoring strategy for below-detection literature statements so they are numerically stable and visibly labeled.
10. [x] Add field-level notes for transferability limits whenever the source experiment is not the benchmark matrix or not the benchmark process window.

#### Acceptance criteria

- [x] Soy process-state priors are no longer based only on repo literature synthesis when the SLR already supplies a direct anchor.
- [x] PPI/SPI aldehyde retention is grounded in explicit payloads rather than only family-level static factors.
- [x] Every retention payload states whether the chemistry is reversible, likely non-covalent, or only indirectly inferred.

### 1C. Flavor reference payloads now missing from the repo

#### Goal

- [x] Promote non-sulfur flavor markers and cross-pathway markers into first-class reference artifacts.

#### Files

- [x] data/lit/flavor_reference_payloads.json
- [x] src/reporting.py
- [x] src/recommend.py

#### Track 1C deliverables

- [x] A flavor reference file that no longer collapses meat-like quality to sulfur anchors only.
- [x] Explicit PBMA-vs-meat reference bands for Strecker, pyrazine, furfural, and key carbonyl markers.

#### Exact tasks

1. [x] Add Hernandez 2023 sub-payloads for:
  - [x] 2-methylbutanal
  - [x] 3-methylbutanal
  - [x] methional
  - [x] phenylacetaldehyde
  - [x] benzaldehyde
  - [x] pyrazines
  - [x] furfural
  - [x] 2,3-butanedione
  - [x] acetoin
2. [x] Add Hofmann 1997 and Hernandez 2023 as explicit MFT/FFT target bands in a flavor reference payload, not only in prose.
3. [x] Add Watanabe 2015 DMHF reference if quantitative meat values can be extracted robustly enough.
4. [x] Add benchmark_role tags such as:
  - [x] reference_anchor
  - [x] pbma_counterexample
  - [x] directional_comparison_anchor
  - [x] low_confidence_mechanistic_anchor
5. [x] Define which reference compounds should affect optimization, which should affect reporting only, and which should remain future-facing until chemistry support improves.

#### Acceptance criteria

- [x] The reference surface is no longer sulfur-only.
- [x] PBMA-vs-meat sub-markers can be surfaced in reports even if their predictive chemistry is still partial.
- [x] The user can see when a formulation is close to PBMA practice but still far from meat-like multivariate quality.

### Track 1 execution order

1. [x] Fill safety_reference_payloads.json first because those payloads are easiest to encode and lowest risk.
2. [x] Fill process_state_calibrations.json and retention_reference_payloads.json second because Track 2 depends on them.
3. [x] Fill flavor_reference_payloads.json third so reporting and optimization can consume the new reference surfaces together.
4. [x] Update slr_incorporation_matrix.json after each ingestion batch so the audit ledger stays authoritative.

### Track 1 done definition

- [x] All Track 1 papers listed in the audit result are represented in one machine-readable artifact.
- [x] No Track 1 paper remains only in docs/slr_benchmark_evaluation.md.
- [x] Each ingested payload has a declared consumer path: reporting, matrix calibration, headspace, safety, scoring, or future-only.

## Track 2: Convert SLR Knowledge Into Production Logic

### 2A. Dynamic retention and release

#### Goal

- [x] Replace static retention assumptions where the SLR already says the process is temperature-dependent, reversible, or browning-dependent.

#### Files

- [x] src/headspace.py
- [x] src/matrix_correction.py
- [x] src/recommend.py
- [x] tests/unit/test_headspace.py

#### Exact tasks

1. [x] Introduce a retention/release model that can depend on:
  - [x] temperature
  - [x] process_state
  - [x] matrix family
  - [x] browning index or projection severity
2. [x] Use Ince 2024 and JAFC 3c05991 to drive aldehyde reversibility rather than a single static factor.
3. [x] Use Xu 2023 to create a cheap temporal attenuation prior for SPI off-flavors.
4. [x] Keep the melanoidin sulfur surrogate, but tie it explicitly to a browning or severity narrative in reporting.

#### Acceptance criteria

- [x] Hexanal release differs meaningfully between low and high temperature soy conditions.
- [x] Projection metadata shows which part came from retention, release, calibration, and browning-linked sulfur trapping.

### 2B. Strecker flavor axis

#### Goal

- [x] Stop treating the flavor surface as if sulfur alone explained meatiness.

#### Files

- [x] src/smirks_engine.py
- [x] src/pipeline.py
- [x] src/sensory.py
- [x] data/lit/flavor_reference_payloads.json
- [x] tests/unit/

#### Exact tasks

1. [x] Promote Strecker aldehydes to explicit target compounds in scoring and reporting.
2. [x] Use Hernandez 2023 values as PBMA-vs-meat reference bands.
3. [x] Make 2-methylbutanal a secondary meaty-quality marker because the SLR identifies it as one of the clearest beef discriminators.
4. [x] Add Lincoln 2025 as a qualitative cross-talk prior for glucose plus lipid plus polyphenol systems.

#### Acceptance criteria

- [x] Reports can explain when a formulation is sulfur-positive but Strecker-poor.
- [x] Comparison outputs no longer collapse everything meaty into one scalar without cross-marker context.

### 2C. Pyrazine control layer

#### Goal

- [x] Encode the SLR's pH-, sugar-, and peptide-dependent pyrazine findings as a real control surface.

#### Files

- [x] data/lit/computational_priors.json
- [x] src/conditions.py or new src/pyrazine_control.py
- [x] src/pipeline.py
- [x] src/reporting.py

#### Exact tasks

1. [x] Add pH-linked pyrazine propensity priors from Laemont 2023 and Arsa 2022.
2. [x] Add sugar-type pyrazine-vs-furfural priors from Laemont 2023 and Hao 2025.
3. [x] Add peptide-sequence and peptide-size priors from Wang 2021 papers.
4. [x] Introduce a pyrazine burden or overshoot term in optimization and comparison reporting.

#### Acceptance criteria

- [x] The repo can show why high-pH or fructose-heavy systems become pyrazine-heavy.
- [x] Scientists can see when a candidate wins on sulfur but loses by roasted or earthy overshoot.

### 2D. Furanone layer

#### Goal

- [x] Add positive furanone chemistry to the roadmap as a first-class axis rather than letting furfural stand in for all oxygenated flavor chemistry.

#### Files

- [x] src/smirks_engine.py
- [x] src/sensory.py
- [x] data/lit/flavor_reference_payloads.json
- [x] data/lit/computational_priors.json

#### Exact tasks

1. [x] Encode HEMF and DMHF as candidate target compounds in the species and sensory layer.
2. [x] Use Blank and Fay 1996 as the mechanistic prior for pentose plus alanine and glycine routes.
3. [x] Use Watanabe 2015 as a meat-side reference if quantitative extraction is robust enough.
4. [x] Keep their confidence low until matrix-specific measurements exist.

#### Acceptance criteria

- [x] The chemistry layer can represent HEMF and DMHF.
- [x] Reports show them as mechanistically expected but low-confidence where appropriate.

### 2E. Thiamine provenance layer

#### Goal

- [x] Distinguish native SPI/PPI systems from additive-rich PBMA systems that can source MFT from thiamine pathways.

#### Files

- [x] src/smirks_engine.py
- [x] src/recommend.py
- [x] data/lit/computational_priors.json
- [x] src/reporting.py

#### Exact tasks

1. [x] Add an explicit formulation metadata field for thiamine availability.
2. [x] Use Cerny 2007 as a mechanistic split prior for MFT attribution.
3. [x] Mark thiamine-assisted sulfur formation as inactive by default for native SPI/PPI benchmark conditions.
4. [x] Surface this as a provenance flag in scientist-facing outputs.

#### Acceptance criteria

- [x] A commercial PBMA-like formulation and a native SPI/PPI benchmark are no longer interpreted as if they shared the same sulfur source mix.

## Track 3: Show Calibration To The User Properly

### Goal

- [x] Make calibration visible, legible, and decision-useful.

### Files

- [x] src/reporting.py
- [x] tests/unit/test_usability_reports.py
- [x] Single-run and comparison report outputs generated by src/reporting.py

### Exact tasks

1. [x] Add a compound-level evidence ladder to single-run reports with direct_anchor, transferred_prior, mechanistic_surrogate, and computational_refinement fields.
2. [x] Extend comparison reports to include MFT/furfural ratio, meaty_quality_penalty, pyrazine burden, Strecker support marker, sulfur trapping summary, and benchmark neighborhood context.
3. [x] Add a calibration summary block per run listing the top literature anchors actually used.
4. [x] Add a missing-data block listing key target compounds that are still hypothesis-only or structurally unsupported.
5. [x] Add a benchmark-neighborhood summary that distinguishes free-system anchor, matrix transfer, hydrolysate proxy, and structural gap.

### Acceptance criteria

- [x] A scientist can tell, from one report, which literature anchors drove the prediction.
- [x] A comparison artifact can reveal when the winning formulation is chemically suspect even if its target score is high.

## Track 4: Immediate Incorporation Sequence From Highest To Lowest Value

### Phase 1: Must encode now because the SLR already contains the data

1. [x] Add Foods 12(10):1967 and Choi 2024 to safety payloads.
2. [x] Add MDPI Foods 11(21):3505 to soy process-state calibrations.
3. [x] Add JAFC 3c05991 to retention payloads.
4. [x] Add Hernandez 2023 Strecker and pyrazine sub-benchmarks.
5. [x] Add Karolkowski 2021 and Xu 2023 to pH/time release priors.

### Phase 2: Must convert to production logic now

1. [x] Replace static aldehyde retention with a temperature-aware reversible model.
2. [x] Add pyrazine burden logic.
3. [x] Add Strecker markers to the score and reporting surface.
4. [x] Expose calibration evidence in comparison reports.

### Phase 3: Important but lower-confidence literature-derived additions

1. [x] Add Troise 2018 soy thermal-history payload.
2. [x] Add Lincoln 2025 polyphenol-modulated cross-talk prior.
3. [x] Add HEMF/DMHF mechanistic priors.
4. [x] Add thiamine attribution layer.

## Track 5: Computational Alternatives For The Remaining True Gaps

### Principle

- [ ] Only unresolved gaps that remain benchmark-relevant after Tracks 1 to 4 enter this lane.

### Gap map to method class

#### 5A. Sulfur selectivity and sulfur-loss chemistry

- [ ] Best cheap method: local sensitivity ranking plus explicit meaty-quality and trapping constraints.
- [ ] Best intermediate method: GFN2-xTB plus CREST conformer search and microhydrated motif screening.
- [ ] Best DFT method: r2SCAN-3c geometry plus frequency workflow, with single-point omegaB97M-V or wB97X-D/def2-TZVPP and cluster-continuum water correction on the top motif set.
- [ ] Best focal check: DLPNO-CCSD(T1) or similarly cheap coupled-cluster correction only for 1 to 3 decisive barriers.
- [ ] Best external-MLP option: evaluate published universal or reactive foundation potentials only as an offline prescreen for conformers or motif stability, never as the sole source of barrier values without chemistry-family validation.

#### 5B. Peptide-bound cysteine reactivity

- [ ] Best cheap method: peptide-state priors driven by hydrolysis level, peptide length, and terminal-residue rules.
- [ ] Best intermediate method: capped dipeptide and tripeptide motif panel screened with GFN2-xTB plus CREST in water.
- [ ] Best DFT method: selective DFT on the top motif families only after xTB ranking stabilizes.
- [ ] Best ML option: local delta model trained on motif-level DFT corrections only if many related peptide motifs must be screened repeatedly.
- [ ] Best external-MLP option: test whether a published foundation MLP reproduces relative peptide conformer ordering well enough to accelerate prescreening before xTB, but require a benchmark against xTB or DFT on a held-out motif set.

#### 5C. Retention and release of aldehydes and sulfur volatiles

- [ ] Best cheap method: temperature-aware reversible retention surrogates fitted to literature payloads.
- [ ] Best intermediate method: representative residue-cluster binding calculations on glycinin or beta-conglycinin motif clusters with xTB.
- [ ] Best DFT method: r2SCAN-3c or similar for motif clusters, but only for calibration deltas, not whole-protein simulation.
- [ ] Best ML option: not first-line; only justified if many repeated motif classes require rapid rescoring after a sufficient DFT seed set exists.
- [ ] Best external-MLP option: use published local-structure potentials only for geometry relaxation or rapid clustering of binding motifs; do not trust internet MLP energies as final retention calibrations without a DFT cross-check set.

#### 5D. Wet-matrix acrylamide kinetics

- [ ] Best cheap method: endpoint safety bands plus free-system kinetics plus explicit uncertainty inflation.
- [ ] Best intermediate method: none unless new internal wet-matrix data exists.
- [ ] Best DFT method: very low priority, because aqueous acrylamide prediction is dominated by matrix-scale transport and availability rather than a single missing gas-phase barrier.
- [ ] Best ML option: not justified before wet-lab kinetics exist.

#### 5E. Furanones and pyrazines in matrix systems

- [ ] Best cheap method: literature-derived propensity priors plus explicit confidence demotion.
- [ ] Best intermediate method: xTB motif ranking for sugar plus amino-acid or peptide branching.
- [ ] Best DFT method: only for the handful of branching steps that materially move furanone versus furfural or pyrazine versus sulfur competition.
- [ ] Best ML option: local reaction-family surrogate only after a curated DFT seed set exists.

## Track 6: Bleeding-Edge Architecture For Computation Without Breaking Practicality

### 6A. Recommended layered architecture

#### Layer 0. Literature and calibration layer

- [ ] Machine-readable payloads for every usable SLR datum.
- [ ] Evidence-tiered priors and censoring surrogates.
- [ ] User-facing calibration ledger.

#### Layer 1. Cheap mechanistic surrogate layer

- [ ] Projection-severity surrogates.
- [ ] Reversible retention models.
- [ ] Pyrazine and Strecker propensity priors.
- [ ] Thiamine attribution flags.

#### Layer 2. Semiempirical screening layer

- [ ] GFN2-xTB plus CREST plus ALPB or GBSA water.
- [ ] Motif library for sulfur, peptide, and retention clusters.
- [ ] Automatic uncertainty tagging based on motif novelty and conformational spread.

#### Layer 3. Selective DFT refinement layer

- [ ] r2SCAN-3c or B97-3c for affordable geometry and thermochemistry.
- [ ] Cluster-continuum water treatment for microhydrated proton transfer and sulfur steps.
- [ ] Higher-level single points only on the smallest decisive subset.
- [ ] Write-back as cached surrogate_patch JSON artifacts.

#### Layer 4. Local ML potential layer

- [ ] Only for repeated motif classes that survive DFT triage and need many related rescoring calls.
- [ ] Use active learning with committee uncertainty.
- [ ] Train only local family-specific models, not a fantasy global Maillard potential.
- [ ] Restrict use to offline acceleration of motif scans, never as the sole source of a new production correction.

#### Layer 4B. External foundation MLP option

- [ ] Offer the user a second lane based on published or internet-available SOTA MLPs.
- [ ] Restrict this lane to offline exploration, prescreening, geometry relaxation, or uncertainty-guided candidate narrowing.
- [ ] Require a calibration card before use:
  - [ ] model name and source
  - [ ] domain the model was trained on
  - [ ] whether sulfur chemistry, charged states, and aqueous motifs are in-domain
  - [ ] validation subset against xTB or DFT on Maillard-relevant motifs
  - [ ] failure criteria that disable the model for a motif family
- [ ] Surface both options to the user explicitly:
  - [ ] local in-house delta MLPs built after DFT seeding
  - [ ] external foundation MLPs used as prescreening accelerators

### 6B. Recommended methods by realism and value

- [ ] Default geometry and conformer workflow: CREST plus GFN2-xTB.
- [ ] Default affordable DFT refinement: r2SCAN-3c.
- [ ] Default decisive single-point refinement: omegaB97M-V or wB97X-D class with def2-TZVPP and solvation.
- [ ] Default highest-accuracy spot check: DLPNO-CCSD(T1) on the smallest decisive motif subset.
- [ ] Default local MLP strategy if justified: MACE-style local potential or equivalent committee-based local force field, trained only on motif families with enough DFT labels.

### 6C. Hard gate for ML potentials

- [ ] Do not start with MLPs to compensate for missing literature.
- [ ] Do not train an MLP until a motif family has at least dozens to hundreds of internally consistent DFT-labeled structures or pathway points.
- [ ] Do not let an MLP create a new production correction without a benchmark-visible validation step.
- [ ] Do not use an external internet MLP as authoritative chemistry unless it has passed a Maillard-specific validation card on our motif families.
- [ ] Always let the user choose between the external-MLP lane and the in-house-local-MLP lane when both are technically viable.

## Track 7: DFT And MLP Approval Gates

### Approve selective DFT only if all are true

- [ ] The gap remains benchmark-visible after literature encoding and cheap surrogates are in place.
- [ ] The missing parameter affects ranking, observability, or safety in a visible way.
- [ ] The target can be compressed into a reusable surrogate artifact.
- [ ] The validation target is explicit in Docker.

### Approve a local ML potential only if all are true

- [ ] A DFT seed set already exists for one repeated motif class.
- [ ] The same family will otherwise require many near-redundant DFT calculations.
- [ ] Uncertainty can be estimated and surfaced.
- [ ] The model is used only offline to accelerate motif exploration or delta prediction.

### Approve an external SOTA MLP only if all are true

- [ ] The model can be run locally or reproducibly through the project workflow.
- [ ] Its training domain plausibly covers the motif class being screened.
- [ ] It passes a small Maillard-specific validation subset against xTB or DFT.
- [ ] It is used as an accelerator or prescreener, not as unverified final truth.

## Track 8: User-Facing Success Criteria

- [x] Single-run report shows the top literature anchors used.
- [x] Single-run report shows which key compounds are anchor-driven versus prior-driven versus surrogate-driven.
- [x] Comparison report shows MFT/furfural ratio, pyrazine burden, Strecker support, and sulfur trapping summary.
- [ ] The repo can show, for each missing gap, whether the next action is literature encoding, wet-lab acquisition, xTB screening, DFT refinement, or no action.

## Review Note — 2026-03-19

- The previous backlog was still too chemistry-general. The next correct move is to finish the SLR-to-repo audit as a first-class artifact, ingest the literature data that already exists but is not yet encoded, expose calibration provenance to the scientist, and reserve DFT or local ML potentials for the much smaller set of gaps that remain genuinely unresolved after that work.
- [ ] current cheap model and why it failed
- [ ] chosen method
- [ ] solvent model
- [ ] charge and spin
- [ ] expected outputs
- [ ] surrogate write-back rule
- [ ] validation target in Docker

### Default cheap-first method by target class

- [ ] sulfur selectivity: barrier sensitivity plus quality-ratio penalty before DFT
- [ ] peptide accessibility: priors plus semiempirical conformer screening before DFT
- [ ] matrix release: calibration registry plus reversible release surrogates before DFT
- [ ] browning-dependent sulfur trapping: severity surrogate before DFT

### Current candidate list to maintain

- [ ] Keep the active target manifest in data/lit/computational_gap_closure_targets.json.
- [ ] Update it whenever a gap is closed cheaply so DFT does not expand by inertia.

## Verification Plan

### Unit and focused integration lanes

- [ ] tests/unit/test_headspace.py
- [ ] tests/unit/test_budget_projection.py
- [ ] tests/unit/test_pipeline_projection_contract.py
- [ ] tests/unit/test_bayesian_optimizer.py
- [ ] tests/unit/test_usability_reports.py
- [ ] focused ranking and contract tests in tests/integration/test_recommendation_engine.py

### Docker commands

- [ ] ./scripts/docker_maillard.sh pytest tests/unit/test_headspace.py tests/unit/test_budget_projection.py tests/unit/test_pipeline_projection_contract.py tests/unit/test_bayesian_optimizer.py
- [ ] ./scripts/docker_maillard.sh pytest tests/unit/test_usability_reports.py
- [ ] ./scripts/docker_maillard.sh pytest tests/integration/test_recommendation_engine.py::TestMaillardPipelineEvaluation::test_evaluate_all_returns_results tests/integration/test_recommendation_engine.py::TestMaillardPipelineEvaluation::test_evaluate_all_result_attributes tests/integration/test_recommendation_engine.py::TestMaillardPipelineEvaluation::test_evaluate_all_results_sorted
- [ ] ./scripts/docker_maillard.sh run python scripts/run_campaign.py --spec data/campaigns/heated_soy_tradeoff_screen.yml --output-dir results/heated_soy_tradeoff_screen

### Promotion rule

- [ ] No new correction is considered complete until both focused tests and one Docker campaign or benchmark-facing artifact show the expected direction.

## Immediate Next Sequence

### Phase 0: Already implemented or in motion

- [x] Make compound-specific matrix calibration affect the observable projection.
- [x] Add heated soy 2-pentylfuran suppression surrogate.
- [x] Add meaty-quality penalty from MFT/furfural ratio.
- [x] Add a first dynamic melanoidin trapping surrogate in the projection layer.
- [x] Create a heated-soy campaign spec.
- [x] Create a cheap-first computational gap target manifest.

### Phase 1: Short-term completion

1. [ ] Run the heated soy campaign in Docker and inspect rank shifts and sulfur penalties.
2. [ ] Add report-facing visibility for melanoidin_trapping_factor if scientist readability needs it beyond raw projection metadata.
3. [ ] Re-audit whether the new sulfur surrogate improves or harms free-system benchmark behavior.

### Phase 2: Cheap-first gap closure

1. [ ] Implement peptide-bound cysteine priors rather than new reaction families first.
2. [ ] Expand heated-state adverse-marker calibration only where a real anchor exists.
3. [ ] Add uncertainty labeling to computational_gap_closure_targets.json if it becomes part of scientist-facing reporting.

### Phase 3: DFT triage only after cheap closure

1. [ ] Run sensitivity analysis to rank which sulfur reactions truly deserve DFT.
2. [ ] Approve only Tier 1A or Tier 1B targets for the first offline batch.
3. [ ] Emit cached correction artifacts, not ad hoc notebook outputs.

## Review Note — 2026-03-19

- The repository no longer needs more abstract strategy before acting. The next useful work is to tighten calibration in the observable layer, keep matrix realism explicit, and use DFT only as a narrow offline correction tool once the cheap surrogates stop buying meaningful accuracy.

## Validation Follow-Up — 2026-03-19

1. [x] Repair broken validation/report-generation imports so Docker validation lanes complete end-to-end again.
2. [x] Expand validation_overview to include quantitative matrix-only benchmarks without confusing them with strict-gate free-precursor trust.
3. [x] Audit all numeric benchmarks in Docker, including Trikusuma 2019, for current coverage and ratio error.
4. [x] Tighten matrix-only calibration until numeric benchmarks are quantitatively defensible.
5. [x] Re-run summary/index/figures and focused scientific tests after calibration updates.

## Review Note — 2026-03-19

- The validation surface now shows all quantitative benchmarks, including matrix-only and matrix-augmented lanes, while keeping strict-ready trust visually distinct from broader executable evidence.
- Matrix-only observable calibration must be applied exactly once. In this repo, `HeadspaceModel.get_matrix_benchmark_headspace_factor()` already carries the matrix observable factor, so `benchmark_validation` must separate proxy release from total observable scaling rather than multiplying the registry factor twice.

## Validation Contract V2 — 2026-03-19

1. [x] Add a second scale metric to the validation contract so benchmark pass/fail no longer depends only on worst-case max ratio.
2. [x] Support explicit per-benchmark scale thresholds in benchmark payloads for the authoritative free-precursor set.
3. [x] Surface the new validation metric in benchmark summary reporting.
4. [x] Document the tool's precision tiers and what they depend on in README.md, including the distinction between kinetic barriers and observable-layer calibration.
5. [x] Re-run focused validation tests and regenerate benchmark summary artifacts after the contract update.

## Track 9: Codebase Restructuring (Clean Code)

### Goal
- Repartition massive logic clusters into focused modules without altering science or validation boundaries.
- Separate reaction templates from the SMIRKS engine.
- Separate projection math from recommendation logic.
- Fix naming misalignment in pipeline orchestration.

### Exact tasks
1. [x] Move structural template logic (`_amadori_cascade`, etc.) out of `src/smirks_engine.py` into a new `src/reaction_templates.py`.
2. [x] Move projection logic (`ProjectionBudget`, etc.) out of `src/recommend.py` into a new `src/projection.py`.
3. [x] Rename `src/inverse_design.py` to `src/pipeline.py` and `InverseDesigner` to `MaillardPipeline`. Update imports in `scripts/run_pipeline.py` and tests.
4. [x] Run full unit/integration test suite (`./scripts/docker_maillard.sh pytest`) to verify absolute parity.
5. [x] Legacy Code Cleanup (Dead code removal)
    - [x] Remove `src/pyrazine_control.py` (Dead code)
    - [x] Remove `src/skala_refiner.py` & `tests/qm/test_skala_refiner.py` (Superseded by `dft_refiner.py`)
    - [x] Remove `scripts/migrate_results_to_db.py` (One-off script)
    - [x] Remove `tmp_mech.yaml` (Root leftover)
    - [x] Update `test_contracts.py` & `test_backend_plumbing.py` to use `DFTRefiner`.
6. [/] Phase 2: Presentation & Logic Separation
    - [ ] Create `src/presentation.py` and move Markdown/CLI rendering logic.
    - [ ] Create `src/projection_utils.py` and move data helpers like `build_projection_rows`.
    - [ ] Deduplicate `reporting.py` and `usability_reports.py`.
    - [ ] Refactor `benchmark_validation.py` to reuse main projection paths.
    - [ ] Run full test suite and verify report outputs.
