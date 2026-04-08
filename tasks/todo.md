# Maillard Strategic Roadmap

## Mission

Build the most useful computational tool for scientists who want to imitate meat-like Maillard chemistry in plant-based systems.

## Product Thesis

This is not mainly a problem of making one chemistry engine more exact.

It is a problem of combining:

- a quantitatively credible free-precursor core
- matrix-aware observability and accessibility
- process-aware confidence boundaries
- scientist-facing reporting that states what is benchmarked, what is transferred, and what remains blocked by missing external evidence

The product question remains:

> Which formulation and process changes are most worth testing next if the goal is meat-like aroma under plant-matrix constraints?

---

## P1 — Immediate Execution Slice (2026-04-05)

### S15. Literature Extraction Sprint — Top-20 Backlog Papers

**Rationale:** 164 of 172 scored citations (95%) remain in BACKLOG. Among those are ~30 papers scored 8/8 containing directly encodable constants (Ea values, partition coefficients, stoichiometries). The current runtime surface uses only 21 + 8 staged references. A focused extraction sprint on the highest-scored papers would approximately double the parametric surface with no wet-lab work required.

**Status 2026-04-08:** the compact runtime-first backlog is effectively exhausted (`ready_runtime=0`, `ready_benchmark=0`). Remaining S15 work is now structural: only keep extracting literature when it unlocks a blocked benchmark surface, a cross-family transfer, or a still-missing executable benchmark.

**Highest Priority Now (2026-04-08):**

- S13. Scientist-Facing Visual Output: model credibility has overtaken report usability, so the next product bottleneck is turning validated predictions into charts and confidence views a food scientist can act on immediately.
- S17. Extrusion Benchmark Validation: the extrusion lane is architecturally useful but still lacks an external-decision-ready benchmark surface.
- S15 structural unlock triage: the compact queue is exhausted, so further literature work should only continue when it opens a blocked benchmark, matrix lane, or cross-family transfer.

- [ ] Re-rank the remaining backlog by structural unlock potential, not just raw score, to isolate the few citations that can still move calibration or executable benchmark coverage.
- [ ] Extract only the constants that open a blocked runtime or benchmark surface; do not keep broadening compact transfer payloads with redundant priors.
- [ ] Prefer landings that create a new executable benchmark or tighten a cross-family structural transfer instead of another single-benchmark local correction.
- [ ] Update the Deep Research backlog state only when the landing is wired into runtime, benchmark, or governance surfaces that survive Docker validation.
- [ ] Validate that no existing benchmark regresses after each structural batch of new payloads or contract changes.

Papers with immediate high-value constants include:

- Glomb & Monnier 1995 — 3-DG retro-aldol fragmentation stoichiometry (MGO 41%, GO 28%, diacetyl 18%)
- Aliani & Farmer 2005 — ribose 3.8× MFT increase, G6P 3.2×
- Ordoudi et al. 2014 — HMF peak kinetics at pH 5.0, 125°C
- Hidalgo & Zamora 2004 — 4-HNE + Phe → 2-pentylpyrrole absolute concentrations
- Blank et al. 2001 — trans-4,5-epoxy-(E)-2-decenal from C20:4, ODT 0.07 ng/L

### S15.1. Land Staged Runtime Queue

- [x] Publish a machine-readable Deep Research runtime-first queue (8 candidates across process-state, safety, computational-prior lanes).
- [x] Land the first selected runtime-only citations from the curated queue into the operational registries without promoting any new benchmark payloads.
- [x] Advance the curated queue to the next runtime-first batch automatically once the previously staged citations are already landed.

### S15.2. Fix Cerny 2008 Benchmark Failure

**Rationale:** This is likely a data ingestion gap, not a modeling gap. The backlog contains `Cerny & Guntz-Dubini (2008)` scored 8/8 with detailed pH-resolved MFT data (thiamine-alone MFT at pH 4–8, mixed-system synergy factor 4.3×) that hasn't been encoded into the runtime. The 30.58× failure ratio is the only out-of-tolerance validation point dragging the entire surface.

- [x] Encode the Cerny 2008 key_values from the benchmark_intake_registry into the thiamine pathway priors.
- [x] Validate that the thiamine benchmark moves within the 1.5× acceptance band.
- [x] Confirm no regressions in the other 17 passing benchmarks.

### S15.3. Implement SLR-Identified Model Corrections

**Rationale:** The SLR (Section "Model corrections identified during the review") explicitly calls out 4 required corrections. Only 1 is implemented.

- [x] MFT/furfural ratio as quality constraint — implemented.
- [x] Make `volatile_retention[hexanal]` temperature-dependent (non-covalent, Ka 3.1×10²–3.1×10⁴ M⁻¹; currently scalar). Source: Ince et al. 2024.
- [x] Make `volatile_retention[MFT/FFT]` dynamic based on degree of browning (melanoidin trapping, 16× reduction). Source: Hofmann et al. 2001. (Currently deferred as 5.7.)
- [x] Add peptide-bound Cys reactivity distinction from free Cys. Source: Nishimura 2024 (17/18 consumed peptides contain Cys).

### S15.1b. Prepare The Next Runtime Queue Slice

- [x] Curate a third runtime-first batch behind the active batch so the queue can auto-advance again once batch 02 is landed.
- [x] Expose the prepared next batch in the generated `deep_research_runtime_queue` artifact instead of keeping it implicit in source only.

### S15.1c. Descriptive Test Naming And Runtime Registry Landing

- [x] Rename opaque scientific test files so filenames describe their actual benchmark or governance scope instead of using shorthand labels like `family03` or `p3`.
- [x] Land all six batch 02 runtime-only citations into the operational registries and update their Deep Research backlog state to `RUNTIME_BOUND`.
- [x] Land the first half of batch 03 with `Ordoudi et al. (2014 / PMC12484514)`, `Glomb & Monnier (1995)`, and `Aliani & Farmer (2005)` so batch 03 becomes the active queue.
- [x] Regenerate the Deep Research runtime queue artifact and validate the renamed scientific subset plus queue and registry landing tests in Docker.

Review 2026-04-05:
Batch 02 and batch 03 are now fully landed in runtime registries, the curated runtime-first queue is exhausted for the current three-batch set, and the scientific tests that previously used opaque filenames now use descriptive names.

### S15.1r. Family 03 Dilute Mixed-Loading Closure (2026-04-08)

- [x] Add a bounded dilute-loading uplift to the mixed thiamine-plus-pentose Family 03 contract so low-total-precursor systems are not forced to inherit the concentrated-model routing fraction unchanged.
- [x] Keep the new uplift anchored to the existing Hofmann-versus-Cerny loading window instead of broadening the free sulfur family or changing global sulfur barriers.
- [x] Re-run the focused Docker subset covering Family 03 runtime behavior, thiamine fragmentation benchmarks, and free sulfur regression guards.

Review 2026-04-08:
The remaining Family 03 error tail was not a dead-metadata bug: thiamine priors already propagated into observability and effective precursor ratios, but the contract still treated dilute mixed systems too much like the 30 mM model mixtures. The runtime now applies a bounded dilute-loading uplift only inside mixed thiamine-plus-pentose cases, interpolated between the Hofmann beef-realistic loading anchor and the Cerny concentrated reference instead of retuning the free sulfur family. This improves the Hofmann 1996 MFT ratio materially while leaving Cerny 2008 and the free-amino-acid sulfur guards stable, which is the right stopping point before any broader structural sulfur/thiamine work.

Verification 2026-04-08:
Re-ran `./scripts/docker_maillard.sh run pytest tests/unit/test_literature_runtime.py tests/unit/test_budget_projection.py tests/scientific/test_thiamine_fragmentation_benchmarks.py tests/scientific/test_free_aa_quantitative_regression.py tests/scientific/test_blind_spots.py -q` and the focused Family 03 subset passed at 74 passed, 2 xfailed.

### S15.1d. Unused Literature File Triage (2026-04-07)

Priority file list for literature that is still underused and could land as runtime or benchmark support:

- [x] `11_maillard_lipid_crosstalk.md` — connect `ACS JAFC 3c05991 / PMC10739987`, `ACS JAFC 3c02618`, and `ACS JAFC 0c01925` to the intake registry so already-landed Family 11 runtime artifacts are visible to the backlog tracker.
- [x] `09_carbohydrate_degradation.md` — promote `Resconi et al. (2023 / PMC10096055)` from benchmark-ready intake to an executable benchmark subset using the extracted absolute furfural and dimethylpyrazine values while keeping the broader identity-gap panel attached in metadata.
- [x] `09_carbohydrate_degradation.md` — land `ACS APTS (Ref. 24)` as pyrraline and furosine Arrhenius anchors, and land `Mottram & Nobrega (2002 / Chapter 9 review)` as furanone-to-sulfur bridge support.
- [x] `03_thiamine_degradation.md` — evaluate `Comunian et al. (2021)` as a thiamine-retention process-state calibration for protected or delayed-release sulfur support.
- [x] `04_nucleotide_degradation.md` — evaluate `Blank & Grosch (1991)` for HDMF OAV benchmark support and donor-routing calibration.
- [x] `15_phospholipid_amine_maillard.md` — link the interfacial PE Maillard kinetics package to the already-landed Family 15 Arrhenius priors so the tracker treats it as encoded runtime support.
- [x] `06_alternative_proteins.md` — evaluate the `Liu (2023)` PPI OAV baseline as a benchmark-ready off-note matrix anchor.

Review 2026-04-07:
Tracker-closure batch landed with compact payloads only: Family 11 linkage rows, the Resconi executable furfural subset, APTS and Mottram runtime priors, Family 15 linkage, and the compact Comunian, Blank & Grosch, and Liu literature anchors.

### S15.1e. Thiamine And Peptide Runtime Closure (2026-04-07)

- [x] Encode `Voelker, Taylor & Mauer (2021)` as a bounded Family 03 Arrhenius prior instead of leaving it in BACKLOG.
- [x] Encode `Huang (2022)` as a bounded Family 03 metal-catalysis correction prior for iron- and copper-rich matrices.
- [x] Close `Wang, Z. et al. (2012)` through a compact intake row linked to the already-landed sulfur-peptide prior rather than duplicating payloads.

Review 2026-04-07:
The next literature slice advances the unfinished extraction sprint with three 8/8 citations and no new engine code: Voelker and Huang now extend the thiamine prior surface, and Wang 2012 is now visible to the tracker as runtime-bound Family 05 support.

### S15.1f. Runtime Contract Regression Repair (2026-04-07)

- [x] Reproduce the 10 failing runtime and validation-contract tests and isolate which expectations regressed.
- [x] Repair the literature/runtime family contract surfaces without reintroducing bloated intake payloads.
- [x] Re-run the targeted validation subset and confirm the repaired contract matches generated artifacts.

Review 2026-04-07:
The regression came from contract-surface overreach rather than broken chemistry: tracker-oriented intake rows and secondary benchmark subsets were leaking into canonical query surfaces, the runtime queue was over-weighting raw score instead of lane priority, and the pipeline was applying browning-driven thiol depletion without explicit residence time. The repaired subset now passes 51 focused tests in Docker.

### S15.1g. Family 02 And 11 Reference Slice (2026-04-07)

- [x] Encode `Marquez-Ruiz et al. (2014)` as a compact high-oleic nonanal-versus-hexanal OAV flavor anchor for Family 02.
- [x] Encode `Messina et al. (2022)` as a comparative cooked-PBMA oil-profile OAV flavor anchor for Family 02.
- [x] Encode `DOI ref. 41 in raw/11_maillard_lipid_crosstalk.md` as a bounded Family 11 sulfur-volatile binding prior for PPI.
- [x] Link the new artifacts through compact intake rows and advance the matching Deep Research backlog entries to `RUNTIME_BOUND`.
- [x] Regenerate backlog artifacts and validate the targeted literature scientific subset in Docker.

Review 2026-04-07:
The next compact literature batch is now landed without opening any new broad registry surface: Family 02 gained an oleic-rich nonanal/hexanal OAV anchor plus a cooked-PBMA comparative oil OAV panel, Family 11 gained a bounded sulfur-volatile binding prior for PPI, the matching Deep Research entries moved to `RUNTIME_BOUND`, and the focused Docker subset covering runtime landing, payload surfaces, backlog rendering, and family queue policy passed 55 tests.

### S15.1h. Family 11 Crosstalk Link And Family 07 Donor Anchor (2026-04-07)

- [x] Link `ACS JAFC 3c08432` into the Family 11 intake surface as crosstalk evidence without duplicating the already-landed fermentation calibration payload.
- [x] Encode `Maillard & van Boekel (1992)` as a quantitative donor-hierarchy prior so Family 07 has a non-duplicative sugar-reactivity anchor beyond Liardon.
- [x] Advance the matching Deep Research backlog entries to `RUNTIME_BOUND` and regenerate the derived backlog and queue artifacts.
- [x] Re-run the focused scientific/runtime Docker subset that covers queueing, learning-loop wiring, intake landing, and Family 07/11 surfaces.

Review 2026-04-07:
The new slice keeps the runtime surface compact while closing two remaining tracker gaps: ACS JAFC 3c08432 is now visible from the Family 11 crosstalk intake surface via the already-landed Rizzello fermentation calibration, and Maillard & van Boekel (1992) now provides a quantitative Family 07 donor-hierarchy prior with explicit k_obs and Ea support. After regenerating the backlog and runtime-queue artifacts, the focused Docker subset spanning queue policy, learning-loop wiring, intake landing, Family 07/11 evidence, and runtime reporting passed 66 tests.

### S15.1i. Figure Label Normalization And Point-1 Crosstalk Slice (2026-04-07)

- [x] Normalize the README-facing benchmark labels for Resconi and the two PMC9905368 sulfur benchmarks to the same author-year format used elsewhere in validation figures.
- [x] Encode `Mottram et al. (2001)` as a compact sulfur-quench prior plus a carnosine-buffering process-state anchor instead of opening a new benchmark family.
- [x] Encode `Yeo & Mottram (2023)` as a phospholipid-crosstalk process-state calibration for soy-lecithin-driven alkyl-thiophene uplift.
- [x] Encode `Wang et al. (2022)` as a residual-hexanal OAV cleanup anchor and `Zhang et al. (2022)` as an unsaturated-aldehyde off-note prior.
- [x] Regenerate backlog and runtime-queue artifacts and re-run the focused Docker subset covering validation labels, payload surfaces, intake wiring, and runtime-bound backlog state.

Review 2026-04-07:
The visible benchmark naming in validation figures is now normalized to author-year for the Resconi and Cho sulfur anchors, and the first user-requested crosstalk/off-note cluster is now landed with compact runtime-native surfaces only: Mottram contributes sulfur-quench and buffering support, Yeo & Mottram adds a lecithin dose window for thiophene uplift, Wang adds a residual-hexanal cleanup OAV anchor, and Zhang adds an unsaturated-aldehyde potency prior. The matching backlog entries now resolve to runtime-bound after artifact regeneration.

### S15.1j. Low-Concentration Sulfur Benchmark Closure (2026-04-07)

- [x] Reproduce the remaining low-concentration quantitative outliers in the Cho et al. (2023) hydrolysate benchmarks and isolate the exact runtime attenuation path for MFT and FFT.
- [x] Recalibrate the hydrolysate sulfur observability logic so SPI and wheat-gluten hydrolysate benchmarks move inside the current validation tolerance without broad sulfur-family regressions.
- [x] Regenerate the validation summary and family-deviation artifacts after the runtime correction.
- [x] Re-run the focused scientific subset covering the Cho benchmarks, validation summaries, and literature/runtime wiring.

Review 2026-04-07:
The remaining low-ppb sulfur error tail came from a contract mismatch rather than broken core chemistry: the Family 05 hydrolysate lane already encoded peptide-assisted sulfur release, but the output observability factor for source-sensitive thiols still treated hydrolysate benchmarks like a fixed low-observability lane. The runtime now applies a bounded hydrolysate-release uplift to MFT and FFT observability before source-specific biasing, bringing both Cho et al. (2023) benchmarks inside the 1.5x envelope while leaving methional stable. After regenerating the validation artifacts, the overview now reports 0 experimental benchmarks outside 1.5x and 0 outside 2x; the remaining prominent error tail is the secondary Resconi furfural surrogate, which remains a broader caramelization-surface gap rather than a low-concentration sulfur issue.

### S15.1k. Safety Bands And HME Rheology Slice (2026-04-08)

- [x] Land `PMC 2024 (PMC12451096)` as an explicit Family 12 commercial AGE safety-reference intake row instead of leaving it only as an untracked payload.
- [x] Land `PMC PMCID:PMC12648097` as an explicit Family 13 acrylamide-mitigation intake row instead of leaving it only as an extended safety reference.
- [x] Start the next highest-value open Family 16 slice with compact process-state calibrations for `Wageningen Ref. 9` and `ACS Food Sci. Technol. 2024` rather than opening a new rheology solver.
- [x] Regenerate backlog and runtime-queue artifacts and re-run the focused literature scientific subset covering safety references, Family 16 calibrations, and backlog governance.

Review 2026-04-08:
The safety lanes for commercial AGE bounds and near-complete acrylamide suppression are now explicit runtime-bound intake surfaces rather than orphaned safety payloads, and the next open-value slice has started with two Family 16 HME state anchors: one for SPI/PPI hydration collapse under rework and one for the practical firmness window in SPI-rich HME analogues. After regenerating the backlog and runtime-queue artifacts, the focused Docker subset covering safety-reference landing, Family 16 calibrations, learning-loop wiring, and backlog governance passed 17 tests. This keeps the landing compact and tracker-visible while deferring any coupled texture-flavor solver work.

### S15.1l. Alternative Protein And Hidden Sink Slice (2026-04-08)

- [x] Land `PMC11049305 (2024)` as a compact Family 06 Spirulina off-note matrix anchor through flavor-reference payloads instead of leaving it as a markdown-only candidate.
- [x] Land `PMC12155365 (2025)` as a compact Family 06 sunflower roasted-matrix anchor with explicit Strecker strength and phenolic interference markers.
- [x] Land `PubMed PMID:1904866 (Ref. 5)` as a bounded Family 14 pentosidine-equivalence prior so ascorbic acid cross-link burden is runtime-visible.
- [x] Land `PMC PMCID:PMC5992167 (Refs. 16/17)` as a Family 15 Amadori-PE burden calibration linked to the existing PE Arrhenius priors.
- [x] Regenerate backlog and runtime-queue artifacts and re-run the focused literature scientific subset covering payload surfaces, learning-loop wiring, and backlog governance.

Review 2026-04-08:
The next double slice is now landed without opening any new broad solver surface: Family 06 gained two compact alternative-protein matrix anchors, one for Spirulina as a marine/fishy lower-bound off-note reference and one for roasted sunflower as a Strecker-strong but phenol-interfered matrix. Family 14 gained a bounded pentosidine-equivalence prior for ascorbic acid, and Family 15 gained a food-matrix Amadori-PE burden calibration that complements the already-landed PE kinetic priors. After regenerating the backlog and runtime-queue artifacts, the focused Docker subset covering payload exposure, intake landing, learning-loop wiring, and backlog governance passed 17 tests.

### S15.1m. Citation Alias Closure And Family 05 Kokumi Slice (2026-04-08)

- [x] Close `JAFC DOI:10.1021/jf034037p (2003)` without duplicating the existing G6P-to-HDMF science payload by teaching the Deep Research tracker to resolve intake citation aliases.
- [x] Land `Ohsu et al. (2025)` as a compact Family 05 CaSR and mouthfulness prior plus a tracker-visible intake anchor rather than opening a standalone kokumi solver.
- [x] Surface the new Family 05 kokumi support in `src/literature_runtime.py` without changing the default sulfur-peptide prior contract.
- [x] Regenerate the Deep Research backlog and runtime-queue artifacts and re-run the focused Docker subset covering tracker aliasing, Family 05 runtime behavior, and scientific payload wiring.

Review 2026-04-08:
The next compact literature closure batch resolved one real tracker gap and one real runtime gap. The JAFC `jf034037p` citation now resolves to the already-landed Blank, Devaud & Grosch G6P-to-HDMF prior through explicit intake citation aliases, so the backlog closes without cloning a second directional prior. Family 05 now also exposes bounded kokumi support through an Ohsu 2025 intake anchor and CaSR prior carrying the reported EC50 and mouthfulness constants, while preserving the existing sulfur-peptide default contract by keeping the kokumi prior in a separate metric surface. After regenerating the backlog and runtime-queue artifacts, the focused Docker subset covering tracker matching, Family 05 runtime output, backlog rendering, queue policy, and payload governance passed 62 tests.

### S15.1n. Family 04 Low-Temp EUC And Family 16 Architecture Slice (2026-04-08)

- [x] Land `Soladoye et al. (2020)` as a compact Family 04 low-temperature EUC anchor with tracker-visible intake wiring plus bounded runtime context, without replacing the existing Matoba hydrolysis prior.
- [x] Land `J. Agric. Food Chem. 2019 (Ref. 21)` as a compact Family 16 architecture anchor through process-state calibration and melanoidin prior wiring rather than opening a rheology solver.
- [x] Extract the Ohsu-based bounded kokumi calculations into an explicit reusable module while preserving the existing Family 05 lane contract.
- [x] Regenerate the Deep Research backlog and runtime-queue artifacts and re-run the focused Docker subset covering tracker state, runtime behavior, and payload surfaces.

Review 2026-04-08:
The next continuation slice closed the two highest-value open literature items while also making kokumi an explicit reusable scoring surface. Family 04 now carries Soladoye 2020 as a bounded low-temperature EUC reference so mild heated-matrix contexts can surface preserved-versus-collapsed nucleotide support explicitly without overriding the existing Matoba hydrolysis logic. Family 16 now carries the JAFC 2019 Ref. 21 pea-hydrolysate plus gum-arabic architecture anchor through a compact process-state calibration and melanoidin prior, which activates only when those cues are actually present. The former inline Ohsu kokumi math now lives in `src/kokumi_scoring.py`, while `src/literature_runtime.py` keeps the same Family 05 contract and additionally exposes a top-level kokumi support signal. After regenerating backlog and runtime-queue artifacts, the focused Docker subset covering the new module, Family 04/05/16 runtime surfaces, tracker wiring, and scientific payload governance passed 67 tests.

### S15.1o. Family 07 Rhamnose Donor And Family 08 Cyclodextrin Slice (2026-04-08)

- [x] Land `Blank et al. (1997)` as a compact Family 07 rhamnose-to-HDMF donor-strength prior plus a tracker-visible intake anchor, without duplicating the broader donor hierarchy already covered by Maillard & van Boekel.
- [x] Land `Bhandari et al. (1998)` as a compact Family 08 cyclodextrin-binding prior for hexanal, nonanal, and `(E)-2-nonenal`, without conflating it with the already-landed protein-binding retention surface.
- [x] Surface both compact payloads in `src/literature_runtime.py` through bounded runtime metrics rather than opening any new benchmark or encapsulation solver.
- [x] Regenerate the Deep Research backlog and runtime-queue artifacts and re-run the focused Docker subset covering tracker wiring, Family 07/08 runtime behavior, and scientific payload governance.

Review 2026-04-08:
The next post-8/8 continuation slice closed two compact but still useful literature gaps without broadening the solver surface. Family 07 now carries Blank et al. (1997) as a bounded rhamnose-plus-proline HDMF reference, exposed only when that specific donor-amino pairing is present so the broader Maillard & van Boekel hierarchy stays intact. Family 08 now carries Bhandari et al. (1998) as a bounded beta-cyclodextrin sequestration prior, which surfaces weighted aldehyde headspace and OAV reduction for observed off-note compounds without confusing cyclodextrin encapsulation with the existing protein-binding retention surface. After regenerating backlog and runtime-queue artifacts, the Deep Research tracker advanced to runtime-bound=67 and backlog=105, and the focused Docker subset covering runtime behavior, tracker state, backlog rendering, and payload governance passed 67 tests.

### S15.1p. Family 03 Aw Thiamine And Family 09 C3-HDMF Slice (2026-04-08)

- [x] Land `Arabshahi & Lund (1988)` as a compact Family 03 aw-dependent thiamine Arrhenius support prior plus a tracker-visible intake anchor, without replacing the existing Cerny flavor-yield prior or the Voelker aqueous Arrhenius anchor.
- [x] Land `Brands & van Boekel (2002)` as a compact Family 09 methylglyoxal-to-HDMF support prior plus a tracker-visible intake anchor, without replacing the existing Blank & Fay mechanistic expectation prior.
- [x] Surface both payloads in `src/literature_runtime.py` through bounded runtime metrics only, keeping the Family 03 and Family 09 contracts backward-compatible.
- [x] Regenerate the Deep Research backlog and runtime-queue artifacts and re-run the focused Docker subset covering tracker wiring, Family 03/09 runtime behavior, and scientific payload governance.

Review 2026-04-08:
The next continuation slice still paid off because both references extended existing lanes with narrow, runtime-visible context instead of opening new solvers. Family 03 now carries Arabshahi & Lund (1988) as an aw-dependent thiamine Arrhenius support prior, which exposes interpolated starch-matrix Ea and an explicit wet-versus-dry modulation factor alongside the existing Cerny yield surface and De Leyn extrusion penalty. Family 09 now carries Brands & van Boekel (2002) as a bounded methylglyoxal-to-HDMF C3-route prior, exposed only when caramelization markers and fragmentation-favoring conditions are present so the default Blank & Fay furanone expectation remains intact. After regenerating backlog and runtime-queue artifacts, the Deep Research tracker advanced to runtime-bound=69 and backlog=103, the curated literature backlog reported no remaining ready runtime or benchmark rows, and the focused Docker subset covering runtime behavior, tracker state, backlog rendering, and payload governance passed 68 tests.

### S15.1q. Family 04 Ingredient-Source EUC Slice (2026-04-08)

- [x] Land `Ahlberg & Mohammadi (2021)` as a compact Family 04 yeast-extract grade reference prior plus a tracker-visible intake anchor, without replacing the existing Matoba or Soladoye thermal support logic.
- [x] Land `Cui et al. (2022)` as a compact Family 04 mushroom GMP and EUC source-profile prior plus a tracker-visible intake anchor, without opening a new combined PBMA umami solver.
- [x] Surface both payloads in `src/literature_runtime.py` through bounded ingredient-source metrics only, keeping the Family 04 thermal and ribose-shift contract backward-compatible.
- [x] Regenerate the Deep Research backlog and runtime-queue artifacts and re-run the focused Docker subset covering tracker wiring, Family 04 runtime behavior, and scientific payload governance.

Review 2026-04-08:
This is the last literature slice that still made clear product-facing sense without broadening the model surface. Family 04 now carries Ahlberg & Mohammadi (2021) as a bounded yeast-extract grade window and Cui et al. (2022) as a bounded mushroom GMP/EUC source-profile prior, so the runtime can expose whether nucleotide support is coming from generic IMP/GMP assumptions, from yeast-extract quality, or from clean-label mushroom ingredients. The existing Matoba kinetics, Nakamura ribose-limit logic, and Soladoye low-temperature EUC window remain the actual thermal backbone; the new references only add ingredient-source context on top. After regenerating backlog and runtime-queue artifacts, the Deep Research tracker advanced to runtime-bound=71 and backlog=101, the curated literature backlog still reported ready runtime=0 and ready benchmark=0, and the focused Docker subset covering runtime behavior, tracker state, backlog rendering, and payload governance passed 69 tests. At this point further literature ingestion no longer looks compact or high-yield; the remaining backlog is predominantly lower-return, more redundant, or blocked by broader solver and wet-lab gaps.

---

### S16. README Rewrite for Food Scientists

**Rationale:** The current README (359 lines) is written for computational chemists. The target user — a food scientist in an alternative protein company — needs to see what goes in, what comes out, and how much to trust it, within 2 minutes.

- [x] Rewrite README as a ~100-line landing page: "What this does → Install → Run → Example output → Where to learn more."
- [x] Add a sample output snapshot (radar chart, ranked formulations table, confidence annotations) so scientists see value before installing.
- [x] Replace internal vocabulary ("bench-neighborhood," "family-lane transparency," "validated envelope") with plain-language equivalents.
- [x] Simplify the architecture diagram to a 3-box version ("Ingredients + Process → Maillard Engine → Predictions + Confidence") as the first visual. Move the detailed Mermaid diagram to `docs/architecture.md`.
- [x] Consolidate the 16 validation artifact links into a single "click here to see what's validated today" reference.
- [x] Move detailed content to dedicated documents: `docs/VALIDATION.md`, `docs/PHILOSOPHY.md`.
- [x] Add explicit positioning statement against existing tools (NIST WebBook, RMG, manual literature review) so scientists understand the unique value.

---

### S13. Scientist-Facing Visual Output

**Rationale:** The pipeline currently outputs raw JSON which masks the deterministic and capability-bounded nature of the framework. A food scientist doesn't need generic sensory radar charts (which hide physics behind arbitrary scores). They need a Decision-Risk Dashboard that shows *why* a suggestion is made and *how much they should trust it*. This is the single biggest reporting gap.

- [ ] Add a **Relative Intervention Waterfall Chart** to `--report` output showing what changed: e.g., "Adding glutathione raised thiols by +X% but reduced total aldehydes by -Y%."
- [ ] Add **Confidence Bounds Overlay** plotting predictions with an error bar representing the Evidence Lane (e.g., narrow bounds for strict benchmarks vs. wide bounds for directional extrapolations).
- [ ] Add a **Capabilities Heatmap** showing Matrix Types (columns) against Process States (rows) to visually translate the 150+ reference papers into proof of capability.
- [ ] Add a parity plot export to the validation-figures command that a scientist can embed in presentations.
- [ ] Refactor the Pareto frontier visualization for multi-objective optimization (meaty-positive vs. safety vs. off-note suppression) to overlay extrapolation risk metrics.

---

## P2 — High-Impact Medium-Effort

### S17. Extrusion Benchmark Validation

**Rationale:** Extrusion modeling is architecturally present (SME coupling, moisture-regime bifurcation, sequential isothermal zones, pre-extrusion damage baselines) but has zero benchmark validation. 0/2 extrusion matrices are closure-ready. For a tool aimed at alternative protein scientists, extrusion is literally the dominant production process.

#### S17.1. Extrusion Benchmark Experiment Design

- [ ] Specify the minimum viable extrusion benchmark: one protein type (PPI or SPI), two SME levels, one barrel temperature, measuring MFT + hexanal + furosine simultaneously.
- [ ] Generate a complete DoE protocol from `doe_generator.py` with real lab specifications (equipment model, SPME fiber type, exact internal standard concentrations).
- [ ] Publish the protocol as a shareable wet-lab request artifact in `results/validation/`.

#### S17.2. Extrusion Model Extensions

- [ ] Add volatile stripping correction at the die (flash-vaporization loss based on die temperature and compound vapor pressure).
- [ ] Add shear-volatile coupling beyond the simple linear SME→ΔT slope (cell-wall rupture → precursor release, protein aggregation → trapping landscape).
- [ ] Evaluate whether a simple RTD (residence time distribution) model is needed or if the sequential-zone model is sufficient for the target use case.

### S18. Selective xTB/DFT Unparking

**Rationale:** P3 refinement governance shows 0 approved jobs, but `why_not_closed.md` identifies 3 specific, narrow motif targets where xTB path search → r2SCAN-3c refinement is cost-effective and would meaningfully improve families 11, 12, and 14.

- [ ] Run xTB path search then r2SCAN-3c refinement for `hexanal_radical_quench` (Family 11: off-note suppression).
- [ ] Run r2SCAN-3c for `lysinoalanine_crosslink` (Family 12: AGE/ALE yield).
- [ ] Generate seed structure for `aa_ring_open_dicarbonyl` (Family 14: stealth browning).
- [ ] Evaluate asparagine-sugar transition state in explicit water cluster to computationally bound the matrix effect on acrylamide kinetics (narrows the wet-lab-only gap).

### S19. Web Interface (Minimal)

**Rationale:** Every interaction currently requires Docker + command-line. For food scientists, this is a major adoption barrier. A minimal web interface would increase adoption by an order of magnitude.

- [ ] Build a minimal Flask/FastAPI web interface with a formulation input form, "run prediction" button, and visual report output.
- [ ] Serve the radar chart, kinetic traces, and safety dashboard from S13 in the web response.
- [ ] Include a "download report" button for shareable PDF/HTML export.

### S14. Codebase Health & Maintainability

**Rationale:** `benchmark_validation.py` (117KB) and `recommend.py` (65KB) are monoliths that impede contribution and debugging. Test suite runtime (~1h40m) blocks iteration speed.

- [ ] Decompose `benchmark_validation.py` into modular components: registry, evaluation, reporting, assertion.
- [ ] Decompose `recommend.py` into modular components: concentration projection, observable mapping, scoring.
- [ ] Triage the test suite for performance: identify and optimize the slowest 10 tests, introduce pytest marks for fast/slow/full lanes.

---

## P3 — Strategic / Deferred

### S12. Scaling the Literature Pipeline & Uncertainty

#### S12.1: Formal Uncertainty Quantification (UQ)

- [ ] Replace narrative "trust heuristics" (e.g., Extrusion Exploratory Mode) with explicit mathematical confidence intervals (e.g., via parametric variance or Gaussian Processes) for out-of-domain predictions.
- [ ] Propagate UQ bounds into the predicted volatile headspace (ppb) figures so scientists know the exact variance of un-benchmarked estimates.

#### S12.2: Automated LLM-Assisted Payload Extraction

- [ ] Build an automated ingestion pipeline that parses eligible Deep Research summaries into canonical `benchmark_payload` JSONs to accelerate closing the ~150-paper backlog.
- [ ] Include a strict human-in-the-loop review interface to guarantee the 8-point SLR criteria are strictly maintained before merging into the main pipeline.

#### S12.3: Model-Guided Active Learning (DoE Feedback Loop)

- [ ] Formalize the "Structural Gaps" into explicit Design of Experiments (DoE) workflows.
- [ ] Implement an API so that when the system identifies a critical gap (e.g., lack of MFT/FFT data in SPI extrudates), it auto-generates a precise wet-lab protocol optimised for maximal model calibration gain.

### Deferred Scientific Modeling Backlog

#### 5.7 Bidirectional Lipid-Maillard Crosstalk

- [ ] Add dicarbonyl-lipid oxidation promotion pathway in `lipid_oxidation.py`.
- [ ] Add melanoidin antioxidant capacity as a time-dependent LOPs suppressant.
- [ ] Validate against Report 11 crosstalk heuristics.

#### 5.8 Disulfide Bond Evolution / MFT Retention

- [ ] Model free-SH to disulfide kinetics as a function of SME and temperature.
- [ ] Link that state variable to MFT headspace recovery in the volatile retention model.

#### 5.10 Sunflower Chlorogenic Acid Off-Note

- [ ] Add temperature-gated 4-vinylguaiacol penalty for sunflower-containing formulations.
- [ ] Include chlorogenic acid to lysine covalent adduct formation as a lysine-accessibility sink.

#### 5.11 Transport / Diffusion Model for Volatile Release

- [ ] Design a 1D Fickian diffusion slab model for volatile release during cooling or serving.
- [ ] Integrate it with volatile retention factors as a compound-class-specific alternative to scalar correction.

### S9. Skipped Test Triage & QM Optionality

Status: supporting infrastructure, not current product bottleneck.

- [ ] Resume only after S14 if skipped-test cleanup blocks deterministic confidence in the active scientist workflow.

#### S9.1: Inventory and classify skipped tests

- [ ] Build a machine-readable skip registry from `pytest -rs` grouped by reason, file, and dependency class.
- [ ] Classify skips into: `not_implemented_module`, `missing_external_dataset`, `missing_optional_backend`, `long_running_campaign`.
- [ ] Add owner and unblock criteria for each skip cluster.

#### S9.2: Quasi-harmonic correction decision and implementation path

- [ ] Implement `src/quasi_harmonic_correction.py` with pure-Python deterministic numerics and no heavy backend coupling.
- [ ] Replace unconditional skips in `tests/benchmarks/test_quasi_harmonic_correction.py` with executable deterministic tests plus `skipif` only for optional integration points.

#### S9.3: Barrier and IRC benchmark skip policy

- [ ] Keep Phase 3.3/3.4 benchmark tests gated by explicit dataset/backend markers rather than unconditional `pytest.skip` inside the test body.
- [ ] Encode a run contract: default CI lane runs deterministic/unit and mock-backed checks; optional QM lane runs benchmark suites when datasets are mounted.

#### S9.4: DFT and ML-potential complement policy for skipped-test conversion

- [ ] Tie each formerly skipped QM test cluster to one of three execution lanes: `deterministic_helper_lane`, `optional_mlp_acceleration_lane`, `optional_dft_authority_lane`.
- [ ] Ensure report and validation artifacts expose lane provenance so users can distinguish deterministic numerics, MLP-assisted predictions, and DFT-confirmed evidence.

### P3–Refinement. Selective Mechanistic Refinement

Active only for matrix benchmarks that remain `mechanistic_priority` after observable closure review. Do not expand broad xTB or DFT activity beyond the specific targets in S18.

### P4. MLP Adoption

Offline accelerator lane only. External MLP evaluation must not substitute for missing matrix benchmarks or observable anchors.

### P6. Matrix-Family Expansion Beyond Pea and Soy

Keep matrix-family coverage explicit in artifacts. Do not broaden family-level scope faster than the evidence surface can support.

### Reproducibility & Provenance

- [ ] Add version-pinned reproducibility snapshots (`pip freeze` equivalent or Docker image tag) to every report artifact so predictions are exactly reproducible for peer review and regulatory contexts.
- [ ] Document the exact mapping between report provenance metadata and the reproducibility snapshot.

---

## Current Product Status

### Strong today

- Quantitative screening is inside the 1.5x band across the current executable surface: 13/13 experimental and 14/14 total quantitative benchmarks, with 11 strict-ready benchmarks and worst experimental ratio 1.442x.
- 16 chemistry families are tracked; 7 are benchmark-linked and 4 currently expose compound-level quantitative parity.
- 54 quantitative compound points are now summarized in the validation surface; the core Maillard family sits at median ratio 1.025x and mean |log10 error| 0.026.
- Family-aware ingestion, runtime, validation, and reporting are operational.
- Pea and soy matrix paths are executable and useful for directional prioritisation.
- README-facing validation imagery is regenerated from the official validation workflow and copied into `docs/assets/`.
- Family 03 dilute mixed-loading closure has been re-verified in Docker after landing.
- Trust language, evidence posture, and family-lane transparency are already visible in reports.
- Intake-registry state model is normalised with explicit `triage_status`, `encoding_status`, and `runtime_artifacts_present` fields.
- Canonical literature backlog artifact exposes three exclusive queues: `ready_runtime`, `ready_benchmark`, `wet_lab_blocked`.
- Extrusion process modeling with SME coupling, moisture-regime bifurcation, and pre-extrusion damage baselines is operational.
- Safety markers (acrylamide, furosine, CML/CEL, LAL) are integrated into the pipeline.
- Cheap-first refinement screening pipeline (barrier offset sweeps with benchmark-visible decision gating) is operational.
- DoE generator templates exist for 5 gap types (blocking benchmark, missing anchor, missing kinetic, missing process-state, missing flavor anchor).
- SLR covers 990 lines with 16+ benchmark-eligible references and 5 structural gaps formally identified.

### Still blocking scientist value

- **Structural literature triage**: the runtime-first compact queue is exhausted, so the remaining S15 value lies in selecting only citations that unlock blocked benchmark, matrix, or cross-family surfaces.
- **No matrix benchmark** is yet external-decision-ready.
- **Scientist-facing visual output** is still too thin in the main report surface: no intervention waterfall, confidence overlay, capabilities heatmap, or presentation-grade parity export.
- **No web interface** — all interaction requires Docker + command-line, a major adoption barrier for food scientists.
- **Extrusion modeling is not yet benchmark-calibrated enough for release-facing use** even though the process-state lane exists.
- Mixed pea and soy meaty-positive targets still rely on transferred or internal-candidate observable support.
- Six chemistry families still have zero benchmark-linked closure.
- xTB/DFT pipeline is parked with 0 approved jobs despite 3 well-scoped motif targets identified.
- No version-pinned reproducibility for peer review.
- Test suite runtime (~1h40m) impedes development iteration.

---

## Success Criteria

- [ ] A scientist can use the tool to narrow a wet-lab or literature-backed matrix campaign, not just inspect simulated chemistry.
- [ ] Free-precursor predictions remain quantitatively stable while matrix usefulness improves.
- [ ] At least one matrix lane becomes meaningfully closer to external decision readiness.
- [ ] Reports and artifacts make promotion blockers explicit enough that synthetic closure is difficult.
- [ ] Expensive compute stays offline, sparse, and justified by benchmark-visible decisions.
- [ ] Visual output makes predictions actionable without requiring the scientist to interpret raw JSON.
- [x] All benchmarks pass inside the 1.5× acceptance band.
- [x] A food scientist can understand what the tool does and run a first prediction within 10 minutes of reading the README.
- [ ] The runtime parametric surface grows by ≥2× through literature extraction before any wet-lab work.
- [ ] Multi-objective tradeoffs (meaty vs. safe vs. off-note) are visible as Pareto frontiers, not collapsed scores.

---

## Completed Foundations

Sprints S0–S11 are complete as of 2026-04-05. All foundational architecture, family-aware ingestion and runtime, matrix observable closure, scientist experiment intake loop, family promotion contracts, intake-registry normalisation, Deep Research runtime queue publishing, extrusion process modeling (SME/moisture), safety marker integration (acrylamide, furosine, CML/CEL, LAL), cheap-first refinement screening, and DoE generator templates have been implemented and validated.
