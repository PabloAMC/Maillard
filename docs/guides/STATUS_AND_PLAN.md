# Status and plan (the one visible plan; updated 2026-10-31)

*Supersedes the plan scattered across STRATEGY_2026-10-29.md secs. 7-9, the pre-registrations and the run logs;
those stay as the record. Sources for every number: two read-only reviews of 2026-10-31 (model state; engineering)
and the artifacts they cite in `results/validation/`.*

## 1. The goal

Given a formulation (amino acids incl. cysteine, sugars, thiamine, buffer, fat) and a process (time, temperature,
pH, vessel), predict the key odorants of alternative-protein flavours -- MFT, FFT, H2S, Strecker aldehydes,
pyrazines, lipid aldehydes -- with honest uncertainty that holds in a laboratory the model never saw, built on the
RIGHT MECHANISMS rather than tuned numbers (owner, 2026-10-30). Measures: leave-one-laboratory-out error per lane,
90 % interval coverage, evidence grounding, blind tests (WUR); judged per STRATEGY sec. 9 (mechanism signatures
lead, calibrated out-of-sample prediction decides adequacy, scores diagnose).

## 2. Where the shipped model stands (B60 calibration; nothing on the branch is shipped)

- Out of sample: 14 / 56 rows within x3 (median 11.8x); the "typical yield" baseline does better (15 rows, 8.1x).
  Within-study contrasts: 0.31 dex against 0.28 for "predict no change". 90 % intervals cover 53 % out of sample.
- Sulfur: 10 / 40 within x3; zero out-of-sample anchors; 0 - 2 % of flux on measured constants; cysteine and
  pentose consumed in seconds where pots measure hours. Acrylamide and lipid: the null model wins. Pyrazines 2.9
  decades low on their own ship rule. Strecker aldehydes, methional, 2-acetyl-1-pyrroline: absent or refused.
  Lanes do not compose (a recipe with both acrylamide and sulfur chemistry is refused).
- What it is good for today: the value-of-information experiment ranking, and ratios inside the Hofmann 145 C family.
  Every lane's absolute prediction is graded do-not-use in the model card. The usefulness bar is not met.

## 3. What the branch (PR #30) has built -- all behind off switches

- Foundations: ~12 two-reader-verified fit tables (~1,800 rows); a fitting module with a real measurement model
  (per-figure scales, series noise, a between-paper error per method class); pre-registration + independent review
  before every fit; the acceptance protocol and a mechanism-signature scorer (15 within-paper tests + isotope routes).
- Findings from evidence (they do not depend on a fit converging): acetic acid transfers across laboratories,
  formic acid does not; the pentose + amino-acid reaction was missing and is now built (Cao hold-out pH RMS
  0.08 at arm E's optimum; the recorded 0.06 was predicted without the trunk solutes, the forked-pool pH bug, re-scored in
  `kinetic_core_joint_fit_entry_ph_rescore.md`); a fitted acetic-acid sink is ruled out by isotope data (CF-069); MFT's route in the lane contradicts
  isotope labelling (Cerny 2003), and the route mix changes with temperature (Cerny 95 C vs Liu 168 C); two
  between-lab conflicts recorded (CF-068: x120 for the same pot).
- Findings from fits, PROVISIONAL (premortem, 2026-10-31): every joint-fit arm so far stopped on its evaluation
  budget (200 evaluations = 11 optimiser iterations for 15 - 16 coordinates), not at an optimum. "MFT's loss ~8x
  slower than shipped" (also CF-013) and "2-acetylthiazole's constant 1000x high" are fit results; at R3En's optimum
  the acetic sink (-0.76 dex) and the shipped MFT decay (-0.02) come back. They stand only once converged fits agree.
  **Converged E (2026-10-31, EX-024): MFT decay at the shipped value (+0.01), the acetic sink partly back (-1.37),
  H2S loss at its lower bound -- "MFT loss ~8x slower" does NOT survive convergence.** The isotope evidence against
  an acetic sink (CF-069) stands, so the fit's sink is absorbing another error (CF-069's first explanation).
- Mechanism signatures: the round-2 model 1 / 15 (despite reduced chi-square ~5); with the pentose entry 5 / 15;
  R3E 5, R3En 4 -- a reshuffle, read off unconverged fits.
  **Chance baseline (EX-025): constants drawn uniformly within the box pass a median 1 of 11 groups (P(>= 3) 0.22),
  so 2 - 3 groups is not yet evidence. Credit goes to the groups that pass <= 5 % by chance: FFT / furfural pH
  direction, temperature dependence, MFT pH direction, Cerny's routes, the entry pH law -- only the entry pH law passes (R3E, R3En).**

## 4. What a well-architected Maillard model needs, and how far we are (owner, 2026-10-31: "keep the big picture")

The chemistry is ONE cascade; the code is several lanes that each carry their own copy of part of it. What a sound
model has, layer by layer, with today's state:

| layer | what it must be | today |
|---|---|---|
| 1. Shared trunk | sugars (pentose, hexose, reducing disaccharides) + ANY amine (free amino acids by class, Cys special, Asn, protein-bound Lys) -> glycosylamine -> Amadori / Heyns -> 1,2- and 2,3-enolisation -> 3- and 1-deoxyosones, plus amine-free caramelisation | split: hexose + glycine in the trunk lane, pentose + Cys / ARP in the sulfur lane, asparagine in the acrylamide lane; the pentose + amine entry exists (off); no amine classes |
| 2. Dicarbonyl fragmentation | the C-C cleavage families (retro-aldol, hydrolytic beta, oxidative alpha) giving C2 / C3 fragments (MGO, glyoxal, diacetyl, glycolaldehyde, hydroxyacetone), acids (formic, acetic), furans (furfural, HMF, furanones) -- each family with class constants | piecemeal per lane; 1-deoxypentosone lacks exits (B7); the formic route unmeasured |
| 3. Strecker | dicarbonyl DONOR (any of layer 2's) + amino acid -> Strecker aldehyde + aminoketone + CO2, identity factors per amino acid | only MGO / GO + glycine (a pyrazine-yield constant, ~2.3 dex low); cysteine's own Strecker in the sulfur lane |
| 4. Heterocycles and thiols | aminoketones -> pyrazines; H2S (cysteine, Strecker, thiamine) + carbonyls -> MFT / FFT / thiazoles / thiophenes by named routes with their own barriers and pH laws | sulfur lane rich but its route mix contradicts labelling (round 3 / B6); pyrazines 2.9 decades low |
| 5. Losses | thiol oxidation / disulfides / adducts, aldehyde aldol / imines, melanoidins | lumped sinks, some fitted to absorb other errors (CF-069) |
| 6. Conditions as physics | pH by a charge balance over every titratable species with the acids the network makes; route-specific Arrhenius; buffer (phosphate) catalysis; water activity / phase; oxygen and headspace; thermal history | pH balance real but blind to the trunk's solutes by default; one lumped barrier for most formation steps; aw and phosphate partly |
| 7. Lipid and matrix | lipid oxidation as a lane that FEEDS carbonyls into layers 3-5 (lipid carbonyls are Strecker reagents; thiol + lipid adducts); protein binding / release and headspace partition | lipid a declared-Q10 lane that does not talk to the others; three capped matrix terms |
| 8. Observation and inference | observation operators per method; a measurement model with between-lab error; posterior ensembles for intervals; leave-one-lab-out, coverage, mechanism signatures; experiment design by value of information | built this week (between-paper term, signatures, register); not yet run on any lane's leave-one-lab-out |
| 9. Software | ONE network assembled from modules (each contributes species, reactions, rate laws, charge ledger, pH law), one explicit configuration object, one integrator, one fit harness | lanes are separate networks that do not compose; ~23 switches are mutable globals set in 87 places |

**What it means for the order of work.** Layers 3 - 5 for any new odorant need the donors of BOTH sugar classes and
every amine, so a new lane built on today's split would duplicate the trunk a third time. The architectural step --
one composable network with an explicit configuration (layers 1 and 9) -- therefore comes BEFORE the Strecker and
pyrazine lanes, right after the sulfur round now in flight. It is a refactor first (bit-identical to today's lanes,
pinned by tests), then the chemistry moves onto it (the pentose and hexose entries as one trunk; amine classes;
donor-generic Strecker). Estimate: 3 - 5 days of work for the refactor and the trunk merge, ~1 day of regeneration and
CI; P(bit-identical on the first pass) ~0.6, the rest found by the pinning tests.

**The mechanisms most likely missing (2026-10-31, after the night's B-section work; P = probability that adding it
changes the alt-protein odorant predictions materially).** The chemistry is qualitatively known (Hodge's scheme,
Yaylayan's "chemical pools"); what the model lacks is the HUB -- the deoxyosones' formation split and fates -- and
the agents that steer it:
1. The amine as a variable: protein-bound lysine / peptide N-termini as the dominant amine in a real alt-protein
   matrix, and amine catalysis of enolisation and of dicarbonyl cleavage (the B7a diagnosis: DPO's beta-cleavage
   fitted on glycine pots and applied amine-free). P 0.85.
2. Oxygen as a reagent: autoxidative (Wolff) glucosone / glyoxal, oxidative alpha-dicarbonyl cleavage, the Namiki
   C2 route from the Schiff base; Hofmann 2000 Table 4 measures x4 more phenylacetic acid in air than under argon
   from the same MGO pot. The engine's oxygen pool is thin. P 0.7.
3. Buffer (phosphate) catalysis of the early steps (the Cao-vs-Zhai ring-opening conflict; well known qualitatively).
   P 0.6.
4. Redox coupling: reductones, the reduction that forms 1,4-dideoxyosones coupled to Strecker degradation
   (Nedvidek 1992), H2S / thiol redox. P 0.5.
5. Donor-generic Strecker degradation (every dicarbonyl x every amino acid) -> Strecker aldehydes, aminoketones ->
   pyrazines, cysteine's H2S. Absent today. P 0.9 for the odorant set (the aldehydes and pyrazines are half of it).
6. Lipid carbonyls crossing into the Maillard network (Strecker donors, thiol traps). P 0.6 for plant-based meat.
7. Matrix and phase: low water activity (real cooking is not a dilute pot), protein binding of thiols (> 98 % of MFT
   bound to pea protein, Xu 2022), melanoidins as sinks / binders. P 0.7 for product-level numbers, ~0 for pot data.

## 5. The plan

**Roadmap (2026-10-31 morning, revised after its premortem the same morning; owner: "do whatever is best", no
designed experiment for now). This ordering governs; the numbered steps below are its detail.** Hours are work;
compute is unattended.

| # | step | addresses (sec. 4 list) | work / compute | decides |
|---|---|---|---|---|
| 1 | **B-close, capped at 24 h (to 11-01 ~08:00):** dossiers of the new papers (Roberts, Parker, Tai & Ho x2, Davidek 2006a/b, Novotny 2008); B7a-2 (DPO's beta-cleavage amine-assisted, with CYSTEINE's amine order declared before the fit -- Zhai's TTCA pots are not amine-free); the stop rule: freeze the thiols at the admitted structure with the measured between-lab spread, leave-one-lab-out overnight. B7b moves into step 4. | 2 | ~4 h / ~10 h | the thiol lane's honest interval |
| 2 | **Odorant layers with forced inputs, now:** Strecker aldehydes scored from MEASURED dicarbonyls (Hofmann 2000's 111 rows, Cremer 2000, Balagiannis 2009), pyrazines from measured aminoketone precursors where printed -- no hub, no refactor needed; each layer judged on its own, so errors stop propagating down the chain. The yardstick extended with the keyed Strecker / pyrazine within-study contrasts, its power reported, paired differences with study bootstraps. | 1 | 17 - 21 h | the first scored odorant lane (~11-07) |
| 3 | **Identifiability audit before the hub:** a within-study factor matrix (for amine, phosphate, O2, pH, temperature: which studies vary the factor INSIDE one study; < ~2 such studies -> a declared band, not a fitted coordinate) and a synthetic model-recovery / power test of the signature machinery (data generated from two rival hub structures with the fitted between-lab noise; does the scorer pick the true one; a third of the signatures held out). | 2, 3, 5 | 10 - 12 h + overnight | which hub factors can be fitted at all |
| 4 | **A3 explicit engine configuration + persistent pool; A4 one fit harness.** | -- | 11 - 13 h | safety of what follows |
| 5 | **C1 the deoxyosone HUB as a standalone module written to the composable-network interface,** seeded from the published multiresponse models (Martins & van Boekel 2005, glucose + glycine; Brands & van Boekel 2002, casein sugars) as null / rivals; B7b (buffer catalysis of the ring opening) folded in; only the factors step 3 admits. | 2, 3, 5, 7 | 2 - 3 days | everything downstream |
| 6 | **A2 the composable network:** the lanes refactored onto the hub module (bit-identical box 3 days), amine classes, oxygen as a species. | enabling | 3 - 5 days | composition of real recipes |
| 7 | **C2 - C5 on the network:** donor-generic Strecker (step 2's forced-input layer, now fed by the hub), pyrazines, thiols re-fitted on the hub (Cerny, Schieberle as tests), lipid carbonyls into the network. | 1, 6 | 7 - 9 days | the odorant set |
| 8 | **C6 cut to a one-page sensitivity on declared bands** (protein-bound lysine, thiol binding, low a_w); product-level numbers out of scope until designed data exist. | 4 | ~4 h | -- |
| 9 | **D adequacy per lane:** leave-one-lab-out with the between-paper error, coverage WITH interval width; Phase 3 regeneration of what passes its audit. | -- | ~1 day / overnight each | what ships |

**Mechanisms before constants (owner, 2026-10-31).** A candidate structure is first tested for FEASIBILITY: can any constant set inside the measured / declared bands reproduce the mechanism signatures (isotope routes, branching ratios, scale-free ratios, presence / absence, the direction of a contrast)? -- by sampling or dry grids, before any fit. A structure that cannot pass is rejected without fitting; fits follow only for structures that pass, and their cost is read after the signatures. Evidence: three B7 arms gained 20 - 70 cost units while failing every scale-free test. Literature mechanism evidence (labelling, fed intermediates, trapping) outranks more fitting.

**Course correction (2026-10-31 ~12:00, after the owner asked to step back and an independent red-team review).** The hub feasibility line (dpo prereg secs. 16 - 21) drifted off this roadmap: it skipped row 3 (identifiability audit + power test) and row 2 (the forced-input Strecker layer), and row 1 never ran its closing leave-one-lab-out. Its signatures became design targets; the trap it found has P ~0.2 and two counter-indications. Decisions: no trap fit; the feasibility method gets negative / positive controls before any more structures; no new global engine switches until A3's configuration object exists; the branch is pushed and split into one PR per sub-model; then rows 2 and 3 in order, and the row-1 leave-one-lab-out overnight.

**How a mechanism is identified (owner, 2026-10-31: "how can we find the true mechanisms?").** Fits to concentration time courses select constants, not structures. Evidence ranks, strongest first:
1. isotope labelling;
2. fed intermediates;
3. trapping;
4. directly measured intermediates and lifetimes;
5. mass balance closure;
6. one factor varied inside one study (order, pH, temperature, O2, buffer);
7. fitted time courses.

Per open mechanism:
- (a) at least two rival structures plus the published null are written BEFORE any result is read;
- (b) the literature experiments of kinds 1-6 that discriminate them are listed, with each rival's prediction;
- (c) the discrimination's POWER is measured on synthetic data from each rival at realistic noise, and an underpowered test is not scored;
- (d) scoring uses the between-lab error, and a single-lab result is PROVISIONAL until one cross-lab check;
- (e) constants are pinned by physical chemistry (measured analogues, classes, pKa) so a structure cannot hide behind free constants;
- (f) where no published experiment discriminates, the item is a declared band and the cheapest discriminating experiment goes to the designed-experiment brief (`rank` verb).

**Mechanism triage (2026-10-31, after its premortem; frozen until the 11-07 checkpoint).** Ranking:
1. Strecker / pyrazine flux with forced inputs (single-lab Hofmann, with a cross-lab check);
2. CF-072 (its own session / PR);
3. thiol losses and cysteine fate, anchored on Cerny;
4. the power / identifiability of the signature machinery, with an exit rule: below 0.6 power, use paired held-out log density or a declared band, 12 h time-box;
5. pyrazines (pyrrole sink);
6. the deoxyosone hub, reclassified from known-wrong to uncertain (single lab, underpowered);
7. oxygen;
8. thermal history and observation operators (OPD, SAFE / SPME, steam distillation), H2S volatility, pH drift in unbuffered pots;
9. phosphate, protein-bound amine, matrix and lipid on the sensitivity page until designed data exist.

Investigation runs in parallel: <= 3 agents, read-only on `src/`, writing standalone modules to [A2_MODULE_CONTRACT.md](A2_MODULE_CONTRACT.md) sec. 3. Each has neutral briefs, a prereg first, positive and negative controls, and "inconclusive" allowed. Implementation is serial: A3, then A2, in a writer session after the B7 branch merges. Forced-input layers are admitted on their own scores, so admission does not wait for A2. P(>= 2 mechanism changes admitted by evidence by 11-21) ~0.5.

**Plan from 2026-11-01 (after the night's atlas and its premortem; supersedes the order above where they differ).**

**Basis.**
- The pathway-condition atlas: [pathway_atlas/README.md](../reference/pathway_atlas/README.md); 32 families, 562 sourced cells, ~9 % wrong per cell; 25 regimes settled, 233 unsettled, 1790 never studied.
- [ENGINE_VS_LITERATURE.md](../reference/pathway_atlas/ENGINE_VS_LITERATURE.md): 20 ranked engine mechanism errors; its sec. 5 red-team amendment is binding.
- A3 phase 1 (Config) merged in PR #34, after an independent review.

**Rules added by the premortem:**
- **(i) No circular scoring.** Before anything is scored, a frozen, hashed split:
  - DEVELOPMENT = every cell cited in ENGINE_VS_LITERATURE secs. 2 - 5 or in any prereg;
  - TEST = whole laboratories not cited there, plus newly keyed rows;
  - lab ids normalised;
  - held-out numbers re-read in the PDF by two readers.
  
  A change's tests come from TEST. The cells that diagnosed it are development checks.
- **(ii) Power before score,** with lab-cluster bootstraps. Route-share changes (E6 / E11) are judged on held-out LABELLING cells, not on contrasts.
- **(iii) Incumbent first.** Each change is scored at the incumbent constants before any refit. A compensation ledger lets a B60 row get worse only if a named open error explains it, declared before the run. One joint refit after E1 + E6 + E11, not three.
- **(iv) A time-box per error item,** and an atlas cell moves the model only after its PDF is re-read.

**Order** (amended by the plan review adopted 2026-11-01 evening, below: an admission gate of >= 2 held-out labs
per error, track H gated on H0, the Strecker layer before E6 + E11; see "Active order from here"):
1. **Land the night:** the ingest PR; EX-056 (the thiol leave-one-lab-out) recorded BEFORE anything it imports changes.
2. **The contrast benchmark, built on the split:** sign and log-ratio per contrast, study / lab bootstraps, power per change. The shipped engine and both null baselines are scored on DEVELOPMENT and TEST separately. A "no change" null has no sign, so the sign score is against the measured direction only.
3. **E4 (the pyrazine pH law off the Strecker step) + CF-077 (the methylpyrazine cross-pair factor):** one owner-approved regeneration, E4 capped at 3 h.
4. **A3 phase 2:**
   - the owner decides review item 1 (refuse non-default build keys) and item 6 (field-level vs dataclass keys);
   - script statics move into jobs; scripts migrate; a persistent pool;
   - a stable digest; an allowlist-based coverage test.
   
   It is on the critical path of E1 / E2 / E10 / E5, which need new switches.
5. **In parallel with 4:** E1 (a Glc + Cys -> 1-DG step) and E2 (per-step enolisation exponents) written as standalone modules to the A2 contract. The forced-input Strecker layer RESTORED (triage item 1): its pre-registered next step is one Amadori decay shared by the fed-ARP and glucose pots (EX-051 follow-up).
6. **After 4:** E6 + E11 (per-route MFT barriers) and E10 via B59b. Then the joint refit; then E5 (buffer catalysis), E3 (the Strecker-aldehyde lane), E7 (lipid coupling). Each: prereg -> feasibility by search -> benchmark on TEST -> red team -> owner ship decision.

**Parallel tracks:**
- **(P1) Keying:** only rows with paired arms, from papers outside DEVELOPMENT, each with a named consumer (a benchmark contrast or an open mechanism), into TEST. Stop if fewer than half of the first 30 qualify.
- **(P2) Owner:** one decision memo with a default for each item (E4 + CF-077, A3 items 1 and 6, the author emails drafted and ready, the E3 + E1 designed-experiment brief).
- **(P3) A literature census** of the most-unknown regimes, low priority.

**Inference track H: hierarchical Bayesian multiresponse and full posterior sampling (owner, 2026-11-01: "add the
full hierarchical sample", with paid Colab compute available; revised after its premortem the same day).**

*Why.* The thiol leave-one-lab-out (EX-056) is NOT ADEQUATE: 6 labs (5 scorable, one holding 341 of 383 held-out
rows), 16 free constants, in-sample between-lab tau 0 - 0.17 dex against ~1 dex held out, 90 % coverage 0.34, and
folds moving constants by 3 - 4.5 dex. In sample the fitted constants act as lab effects. A hierarchical model gives
intervals for a NEW lab; it calibrates uncertainty, it does not fix a wrong structure (P(it alone makes the thiol
lane beat both nulls out of sample) ~0.15), and it can HIDE a wrong structure in lab effects, which the tests below
must catch. Today one full-cost evaluation is 2.5 - 21 s on CPU with finite-difference Jacobians, so the bottleneck
is the code, not the chip.

| # | stage | needs | work / compute | acceptance |
|---|---|---|---|---|
| H0 | **Out-of-sample tau from EX-056, post-processing only.** Nested leave-one-lab-out (tau calibrated without the scored lab), weighted by LAB not row, declared a lower bound (the structure, series noise and warm starts were set on all labs). Second split: another lane's leave-one-lab-out or labs newly keyed by P1; repeated grouped lab partitions only once >= 10 labs exist. | nothing (no A3 dependency) | 4 - 6 h / 1 h | a proper score (held-out log density + interval score) against both nulls, with coverage AND width reported; pre-declared: absolute-class tau > 0.8 dex -> the lane ships relative within-lab predictions only |
| H1 | **Hierarchical specification, pre-registered.** Lab effects on OBSERVATION levels and on bounded latent CONDITIONS (heat-up lag, pH drift, buffer); on a rate constant only where >= 3 labs inform it. Hyperpriors from measured bands and Ea classes; non-centred. Masking tests: (a) lab effects regressed on recorded pH / T / buffer -- a trend names a missing mechanism; (b) lab effects must not change within-lab frozen signature verdicts; (c) a lab spread > 0.3 dex on a physical constant is structural. | H0, identifiability audit | 4 - 6 h | independent review |
| H2 | **Tool bake-off on a forced-input layer** (the Strecker layer with measured dicarbonyls, written to the A2 contract): JAX + diffrax + NumPyro vs a compiled CPU right-hand side, float64. Positive control: simulation-based calibration on synthetic truth. Negative control: a wrong structure WITH lab effects on k, fitted to data simulated on EX-056's lopsided design, must be detected by the masking tests. GPU benchmarked here (T4 / L4 / A100 vs the 8-core CPU) before any queue is sent to Colab. | H1 | 8 - 12 h / hours | trajectories within < 1 % of the observation sigma of the numpy layer; gradients vs central differences; SBC passes; the negative control is caught |
| H3 | **Fast differentiable engine.** First profile; then a 1-day spike on the charge balance alone: pH as an algebraic state of an index-1 DAE in CasADi / IDAS (exact forward and adjoint sensitivities, float64, CPU), against JAX (diffrax has no BDF or DAE; the nested brentq / ionic-strength / warm-start pH solve is the risk; only the inline B64 coupling can be ported). The winner compiles A2's arrays. | A2a merged, H2, a lane passing its structure tests | 3 - 5 days | panel and signature records within the stated tolerance (not bit-identical); gradient checks |
| H4 | **Full hierarchical posterior on each lane that passes its structure tests.** At 16 - 30 constants plus lab effects: SMC or ensemble samplers, or a Gaussian-process emulator with history matching (an explicit model-discrepancy term names structural error); NUTS last, started from Laplace / Pathfinder with a dense mass matrix and a sampler-health stop rule (divergences, tree depth, solver failures). Exact leave-one-lab-out by refitting, compared with Laplace + H0's tau by paired held-out log density. GPU batching (pots x chains x folds) only if H2's benchmark shows a gain. | H3 | 1 - 2 days / 6 - 24 h per run | R-hat < 1.01, bulk ESS > 400, no unresolved divergences; held-out proper score beats Laplace + tau |

*GPU facts (premortem, March 2026 sources).* Stiff kinetics needs float64. FP64 TFLOPS: T4 ~0.25, L4 ~0.47 (both
1/32 - 1/64 of FP32), A100 9.7 (19.5 tensor), H100 26 - 34. Small stiff solves are bound by kernel launches, not
FLOPs, so P(the GPU beats the 8-core CPU by < 3x) ~0.6. Colab Pro / Pro+ buy compute units; an A100 costs ~5 - 8
units/h, premium GPUs are "subject to availability", and 24 h background runs are Pro+ only. Money is not the
limit; disconnects are. Exact gradients cut a 16-constant least-squares fit by ~3 - 8x (mostly compilation), and a
~100-dim gradient by ~100x.

*Rules.*
- A separate `maillard-jax` env (jax / jaxlib, diffrax, numpyro, optimistix, casadi); a gate checks that the core
  `environment.yml` lock is unchanged, so frozen evidence does not drift through a numpy / scipy bump.
- Colab runs are proposed to the owner first; the notebook clones the repo at a recorded sha with a pinned
  requirements file; chains checkpoint to Drive; Colab returns posterior DRAWS only. Scoring, coverage and every
  `results/validation/` artifact are made from the draws locally through `./scripts/docker_maillard.sh`. Nothing
  private (PDFs, the download queue) goes to Colab.

*Calibration (premortem):* H0 by 11-07 ~0.7; H2 passes SBC by 11-14 ~0.35; H3 by end of November ~0.2; a converged
hierarchical posterior on the thiol lane by mid-December ~0.15; given calibration, intervals widen > 2x ~0.85.

**Consolidation wave V (owner, 2026-11-01: "take a step back, send many agents to verify the current information,
solve disagreements, and identify citation chains that may help us finalize the literature").** After the
resolution batch in flight is red-teamed and applied; ~2 days; <= 8 agents at a time; agents write patches to the
scratchpad and ONE writer applies them in order.
- V1 Verify, ranked by influence: every number a prereg, a TEST contrast or a fit row cites, plus a random sample
  of the rest (to re-measure the atlas's ~9 % error rate); figure-read numbers flagged. Verifiers are BLIND to the
  engine's residuals. Corrections to TEST rows land with a new split hash before anything is scored on them. Every
  dossier's DOI line checked against Crossref.
- V2 Disagreements: each open CF entry triaged into resolvable-from-corpus / needs-a-paper / needs-an-experiment;
  one agent per mechanism cluster resolves the first kind.
- V3 Citation chains: best-first search per open mechanism (backward and forward via OpenAlex / Crossref, starting
  from the dossiers with labelling, fed-intermediate or rate-constant evidence), de-duplicated against the corpus;
  stop at < 1 useful candidate per 10 abstracts or a cap of 30 abstracts and 5 full texts per mechanism (pre-1990
  reference lists are patchy in OpenAlex, so read those PDFs' reference lists directly). The ranked download list
  goes to chat with verified DOIs, never into the repo.
- Output: one consolidation record (verified / changed / still open and why); after it every literature item is
  "settled", "needs paper X", or "needs an experiment".

**Missing pieces carried as declared limitations until built:**
- in-situ pH: the probes clamp pH while the literature pots drift (CF-084);
- observation operators;
- intermediate water activity and protein-bound amine, which is where alt-protein recipes actually sit;
- the WUR scoring hook.

**Checkpoint 2026-11-07:**
- (a) the benchmark with its frozen split exists, and the shipped engine is scored;
- (b) E4 + CF-077 have verdicts;
- (c) the A3 phase 2 PR is green;
- (d) E1 is pre-registered as a standalone module.
- (e) H0 (out-of-sample between-lab tau from EX-056) is scored against both nulls.

P(met) ~0.5 (premortem). P(>= 2 evidence-admitted mechanism changes shipped by 11-14) ~0.45.

Rules throughout: premortem of each plan change; pre-registration + independent review before every fit; converged
fits only; chance baselines AND (from step 3) a measured power before credit; one PR per sub-model; a mechanism
without data -> DOIs to the owner.
**Calibration (revised by the premortem):** P(ranked relative predictions with honest intervals in ~4 weeks) ~0.45
with this ordering (~0.30 with the first draft); P(the hub's partition signatures above chance after step 5) ~0.4;
P(a scored Strecker lane by ~11-07) ~0.6 (forced inputs remove the hub dependency); P(absolute odorants within x3
across laboratories for most rows) ~0.15 without designed data.

**A. Consolidate and de-risk (this week; ~1.5 days of work).**
1. Merge PR #30 once CI is green (owner's call; nothing shipped moves). From then on, one PR per sub-model, code
   separate from regenerated results.
2. Rescue the stranded fix in the main checkout (uncommitted since 2026-10-28: the `k_odg_da` band excluding the
   shipped value, a pH fixed-point warning) onto a branch off main (owner's OK; it touches the user's working copy).
3. Engine configuration made explicit: one context manager for every fit's switches (restores on exit), then a
   frozen config passed per job, and a persistent worker pool instead of forking per evaluation (the macOS
   fork-safety hazard and the per-evaluation cost go away). 8 - 10 h.
4. One fit harness in `src/` for the sub-model fits (three drifted copies today). 3 h.

**B. The sulfur lane's mechanisms (time-boxed: R3T + B7, 48 h of wall clock from 2026-10-31; premortem).**
0. FIRST, convergence (premortem #2): no arm verdict is read off a fit that stopped on its budget. Converge E
   (warm start, larger budget, its series noise re-estimated), then fit every round-3 rival from its own optimum
   with E's measurement model FIXED, so the costs are one objective and a nested arm can only match or beat E.
   EX-015 / EX-016 / EX-017 are inconclusive until then. ~1 h of work + ~7 h of unattended compute (night of 10-31).
   **Done 2026-10-31 02:57 (EX-024):** E converged; under E's fixed noise E 1141, R3E 1140 (route off), R3En 1076,
   R3Es 1069 (slope ties zero).
5. Round 3: the pH law of the intact-C5 MFT route. **Verdict (converged): R3En admitted** -- the data accept DDP
   formation without the shipped free-amine law; Cerny pot D fails in every arm (0.98 - 0.99), so the admitted
   structure still makes MFT through fed norfuraneol (R4f, B7a next).
6. The route mix vs temperature (R3T): BEFORE freeing the norfuraneol -> MFT barrier, re-score the optimum with a
   heat-up ramp on the Hofmann 145 C autoclave pots (premortem #3; ~1 h, no fit): a barrier fitted against a misread
   ramp is a tuned number. Then separate barriers for the DDP, norfuraneol and C2 + C3 routes, tested on Cerny 95 C,
   Hofmann 145 C, Liu 168 C.
7. The 1-deoxypentosone partition and the cysteine sink (2026-10-31, R4a review: these OUTRANK further MFT-route arms):
   the engine puts Zhai 2019's fed-TTCA 3-DX / 1-DX at 20 - 330 (measured 2.1 - 2.6), the fitted C2 + C3 drain on DPO
   at +2.4 dex, free cysteine 2 - 40x too high in the fed-TTCA pots (a missing sink), and in glycine pots DPO
   over-feeds norfuraneol (25 - 36 mol% vs <= 5 %). Rival DPO sinks and a cysteine sink, tested on Zhai 2019 / 2021,
   Cao 2017 and Davidek. Then the pH model on unbuffered pots with the pentose acids in. Scoped after its review
   (kinetic_core_dpo_partition_prereg.md sec. 7): B7a the partition only (two-level measurement model; split-law
   and TDP-exit arms); B7b the cysteine fate (buffer catalysis of the ring-chain step) only after a cost-by-block check.
   **Night of 10-31 (EX-030, EX-033 - EX-036):** R4f rejected (TTCA's direct step stays; no in-class barrier passes
   Cerny); the two-level measurement model built (bit-identical off); the admitted structure fails the partition's
   scale-free signatures in time and pH (wrong sign) under either measurement model; NO B7a arm passes (split law,
   TDP exit, both, second MGO parent). Diagnosis: DPO's base-favoured beta-cleavage (fitted on Davidek's glycine
   pots) cancels any split law in amine-free TTCA pots -> next candidate (to design and review): the beta-cleavage
   amine-assisted; TDP's caramelisation source. **Time-box: 48 h from 10-31 00:30 -> B closes 11-02 00:30.**
   Meanwhile (queued after step 0): R4f -- TTCA's direct breakdown to DPO as a free offset (its lower bound = the
   measured Amadori route only), the NF -> MFT barrier profiled at 64 / 100 / 133 kJ/mol (thiol prereg secs. 14 - 16).
8. Adequacy: leave-one-laboratory-out with the between-paper error on the finalist (6 - 10 h compute; overnight or
   Colab), then Phase 3: one regeneration shipping every switch that passes its evidence audit (B59b's de Roos route,
   the entry, the corrections).
   *Success:* the signature groups (below) and Cerny's route tests; out-of-sample coverage 0.80 - 0.97 REPORTED
   WITH the interval width (wide intervals buy coverage); beat both null baselines out of sample. P ~ 0.4.
   *Stop rule (defined, premortem #5):* a plateau is no gain of >= 2 mechanism groups over two consecutive arms, or
   the 48 h time-box; then ship the thiols with the measured between-lab spread as their uncertainty and go to A2 / C.
   The thiol pH arm moves to the backlog.

**Rules adopted from the premortem (2026-10-31).** (a) No verdict from a budget-exhausted fit. (b) Rival arms are
compared under one fixed measurement model (and, for adequacy, by paired held-out log density). (c) An arm may not
lose a signature an earlier arm passed without a stated mechanistic reason (warning sign: a constant moving > 1 dex
between arms that differ by one coordinate). (d) The signature list is frozen as v1 (hashed), counted by mechanism
group, each labelled in-sample or held out, against a chance pass rate from constants drawn within their bands (~4 h).
(e) The scorer reads an arm's configuration from the fit's JSON, not from its file name (~1 - 2 h). (f) A3's config
context manager comes before A2; A2's bit-identical refactor is time-boxed to 3 days.

**A2. One composable network (sec. 4, layers 1 and 9; 3 - 5 days).** After B's round in flight: modules that
contribute species / reactions / rate laws / charge ledger / pH law to one network, one configuration object, one
integrator -- bit-identical to today's lanes first, then the pentose and hexose entries as one trunk, amine
classes, donor-generic Strecker. Prerequisite for C, so new lanes do not copy the trunk a third time.

**C. The lanes alt-protein flavour needs most (each ~2 - 3 days, on the composable network).**
9. Strecker aldehydes (methylbutanals, phenylacetaldehyde, methional): universal in plant proteins; kinetic data in
   the corpus (Balagiannis 2009, Cremer 2000, Chan 1994b). 10. Pyrazines: close the 2.9-decade deficit with the
   aminoketone route's measured barriers. 11. Lipid aldehydes: a measured rate, a pH term, the lipid x thiol
   interaction. 12. Lanes that compose, so a real recipe runs as one system.

**D. Usefulness (in parallel, owner-facing).**
13. Designed data beats more literature: the between-lab spread (10 - 1000x) is the binding limit, so one controlled
    experiment per open mechanism (chosen by the experiment ranking) is worth more than any fit; the WUR blind test
    is scored the day its data arrive.
14. The deliverable for formulators: RELATIVE predictions (rank formulations / process changes by odorant profile,
    with intervals) before absolutes -- within-study contrasts are where a mechanistic model can be trusted first.
    First step (premortem #8, 4 - 6 h): a within-study contrast benchmark (phosphate +/-, pH, hold time, cysteine
    +/-), scored shipped vs branch at every sub-model change. Real products also need protein binding (Xu 2022:
    > 98 % of MFT bound to pea protein) -- layer 7, not yet scheduled.
15. Owner action (premortem #9, ~1 h): send the one-page designed-experiment brief; designed data is the only thing
    that breaks the x10 - x1000 between-lab spread.

**Calibration.** P(the sulfur lane beats both null baselines out of sample within ~2 weeks) ~0.4; P(a useful
Strecker-aldehyde lane within ~1 week of starting it) ~0.5; P(absolute odorant predictions within x3 across
laboratories for most rows) ~0.15 -- the realistic near-term product is ranked, relative predictions with honest
intervals (~0.55).

**Owner decisions, 2026-11-01 (~10:30; defaults approved):**
1. CF-077 (the methylpyrazine cross-pair factor) fixed in its own session, as a PR for the owner's approval.
2. E4: no change; a declared band and a designed experiment (the prereg's decision note).
3. A3 phase 2:
   - refuse non-default build-phase keys;
   - keep dataclass-level keys and add field-level bands for the fields feasibility samples.
4. The E6 + E11 refit runs after the thiol leave-one-lab-out; it does not ship.
5. The coordinator drafts the author data requests and the designed-experiment brief for the owner to review and send.

**Plan review, 2026-11-01 afternoon (cloud session). ADOPTED by the owner, 2026-11-01 evening ("Yes, PR and adopt"):
all six items, as written. Where they differ from the "Plan from 2026-11-01" order above, they govern.** Basis: the
contrast benchmark's power tables (`docs/reference/contrast_benchmark/README.md` secs. 5, 8, 11), the 2026-10-09
brief ([EXPERIMENTS_BRIEF.md](EXPERIMENTS_BRIEF.md): "stop scraping literature"), the plan above.
1. **An admission gate before work starts, not only before scoring.** Under split v2, 8 of 20 errors have < 2
   held-out laboratories (E1, E2, E4, E8, E10, E15, E16, E17); the keyed queue lifts E4 to 3 labs and E10 to 2 once a
   new split version admits it. Rule: no prereg, feasibility search or fit for an error until its TEST(E) has >= 2
   labs in a frozen split version; otherwise it is a declared band plus a designed-experiment line. Default: adopt.
   Effect: step 5's E1 and E2 modules wait; E4 is already a declared band (owner decision 2).
2. **Order after the benchmark scorer and H0:** the forced-input Strecker layer (E3, testable, triage item 1,
   P(a scored lane by 11-07) ~0.6), then E6 + E11 on >= 8 cores, then the pyrazine items E13 / E14 (testable).
   A2 starts when a second lane needs the shared trunk, not before.
3. **Track H gated on H0.** H1 - H4 start only after H0 is scored. If tau_ABS > 0.8 dex the thiol lane ships
   relative predictions only, and H4's gain is then small. Default: adopt.
4. **Wave V cut to what scores.** V1 on the 11 single-reader v2 TEST contrasts (finishes v2) and on the TEST rows a
   scorer reads; V2 for CF entries that name a TEST cell; V3 (citation chains) only for an error that rule 1 blocks
   for lack of labs. Default: adopt.
5. **Owner actions with the largest effect:** send the designed-experiment brief and the 5 author data requests
   (`data/articles/OUTREACH_DRAFTS.md`, laptop only). The between-lab spread (x10 - x1000) is the binding limit; no
   step above removes it.
6. **Profile before the next long refit:** 2 h on one joint-fit evaluation (finite-difference Jacobian, pH solve,
   integrator) before E6 + E11's ~1 day queue; exact Jacobians cut a 16-constant fit ~3 - 8x (H3's first step).

*Premortem of these changes.* (a) Rule 1 starves the mechanism work if P1 keying stalls; P1 needs the PDFs, so it
runs on the laptop (P(>= 2 more errors testable by 11-14) ~0.5). (b) The Strecker layer may fail its blind test
again (EX-051 INCONCLUSIVE, one decay constant cannot serve both pots); its exit is a declared band with the
designed experiment named. (c) Cutting V3 may miss a paper that would settle a CF; rule 1 re-opens V3 per error.
P(checkpoint 11-07 met) with these changes ~0.55 (from ~0.5): fewer parallel tracks, the same dependencies.

**Done 2026-11-01 (cloud session):** PR #36 (A3 phase 2) merged after its full suite; the forked-pool pH bug
re-scored (`results/validation/kinetic_core_joint_fit_entry_ph_rescore.md`: 21 records reproduce their recorded
values to 1e-6; corrected Cao RMS 0.06 - 0.22; the three R4f records flip from fail to within; Zhou unchanged).

**Done 2026-11-01 evening (laptop session):** PR #38 (the ingest) merged; PR #39 merged: four papers keyed (Ruan 2025,
Feng 2024 JSFA, de Sousa Fontes 2024, Feng 2023 Front Nutr), all MARGINAL as data, and the repo's records citing
Feng 2023 corrected against the paper (an unsupported 82 kJ/mol record deleted). No shipped number moved.

**Active order from here (the adopted review, items 1 - 3):**
1. the contrast benchmark scorer on split v2 (its prereg v2 is in `results/validation/`);
2. H0, the out-of-sample between-lab tau from EX-056;
3. the forced-input Strecker layer (E3);
4. E6 + E11 on >= 8 cores, after the 2 h profile (item 6);
5. the pyrazine items E13 / E14.

E1 / E2 wait for >= 2 held-out labs (item 1); H1 - H4 wait for H0 (item 3); A2 waits for a second lane that needs
the trunk. Owner: item 5 (the brief and the 5 author requests).
