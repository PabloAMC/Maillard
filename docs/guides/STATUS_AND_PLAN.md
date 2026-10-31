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
  formic acid does not; the pentose + amino-acid reaction was missing and is now built (Cao hold-out pH within
  0.06); a fitted acetic-acid sink is ruled out by isotope data (CF-069); MFT's route in the lane contradicts
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

## 5. The plan

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
