# Overnight log, 2026-10-29 / 30

Owner asleep; hourly checkpoints against the north star (strategy sec. 1: leave-one-laboratory-out error per lane,
interval coverage, evidence grounding, blind tests) and the purpose (predict MFT, FFT, H2S, Strecker aldehydes,
pyrazines, lipid aldehydes for alternative-protein flavour design). Each entry: what was done, what it showed, whether
the next step still serves the goal, what changed.

## Checkpoint 1 -- 23:15

**Done since the evening.** All four sub-model tables READY (1,587 verified rows; sub-model 1 now includes Martins
2005's acids, two readers). No-refit probes of all four (`results/validation/kinetic_core_*_submodel_probe*.md`),
including arms with the four pre-registered route switches on (B58 TTCA equilibrium, B59b cysteine Amadori route,
B59a Strecker sulfur, B61 norfuraneol hydroxide decay). Inline pH coupling checked against the shipped fixed point
on sub-model 1's table: within 7e-4 (pH or dex) on all 451 rows, 3.6x faster. Fit script for sub-model 1 written and
smoke-tested; an independent design review is running before the real fit (procedure step 5).

**What the probes say, by size of the miss (this reorders the priorities):**

| misfit | size | data | kind | serves the goal via |
|---|---|---|---|---|
| methylglyoxal in pentose systems | 10^2-10^4 x low | 7 Jiangnan papers, calibrated MGO | branching constant: r_dpo_c2c3 exists at ~3e-4 of 1-DX's exits, never fitted to MGO (corrected at checkpoint 2) | C2+C3 route to MFT (Hofmann's labelling), cysteine Strecker -> H2S; pyrazines later |
| MFT / norfuraneol from hexose + cysteine | predicted 0 vs 19-25 ug/L | Hofmann 1998 T1/T2/T5 | missing hexose -> C5 entry (labelling: intact C5; corrected at checkpoint 2) | glucose-based reaction flavours |
| shipped route without B58/B59b | 1-10 dex | 5 papers | switches off | everything downstream; with them on, 0.1-0.8 dex, no refit |
| acid supply's pH dependence | 3 vs 79 mmol/L measured, flat predicted | Davidek, Martins | code path (trunk pH terms not applied in the sulfur lane) | every pH-dependent rate |
| pentose acid supply at 0.4 M | 0.8-0.9 pH units too acid | Cao, Zhao | constant or form | pH of high-concentration reaction flavours |
| thiol pH direction (ribose + Cys) | MFT 46x low, FFT 46x high at pH 7 | Hofmann T2 | CF-012 / CF-017 | MFT / FFT ratio, the core flavour target |
| MFT at 100 C, early | 10-60x low | Schieberle 2000 | rate or route | low-temperature processes |

**Does the plan still serve the goal?** Partly. Sub-model 1's fit is still right to run (it is cheap now, ~7 s per
evaluation, and pH feeds everything). But the two largest misses are MISSING ROUTES, not constants, and both sit
on the path to MFT. Fitting constants around a missing route would bias them (the k_aa_frag lesson). Change: after
sub-model 1's fit, the next work is the literature and the route definitions for (1) pentose fragmentation to MGO
and hydroxyacetaldehyde and (2) the hexose + cysteine route, probed against the Jiangnan MGO rows and Hofmann's
hexose pots before any thiol-constant fit. Sub-model 2's fit is deprioritised: with the pre-registered switches on it
is already within ~x3 on most rows; its open items are TTCA at 120-140 C pH 5.5 (-1.1 dex) and cysteine loss in
Zhai 2021 Fig. 6 (+1.1 dex).

## Checkpoint 2 -- 00:05

**Done.** The design review of sub-model 1 found seven blocking flaws (double-counted priors and drift, glycine
missing from the charge balance, noise that hid misfit, a pH-stat on the wrong pH scale ...); all fixed and the
fixes, tolerances and hold-outs pre-registered (prereg sec. 5) before the fit started at 23:33 (arm B running, 28
min in; arm A after; ~1.5 h each). Hold-out for sub-model 1 (Martins 2005b pH-stat and drift acids) read by two
readers and verified; the frozen-prediction script is ready. Hofmann 1998 Table 9 (mercaptopropanone from glucose +
cysteine) read by two readers into the thiol table. Literature agent: the pentose C2 + C3 step EXISTS (branching
constant never fitted to MGO); hexose MFT keeps an intact C5 (labelling); ranked download list for the morning.
The sub-model 3 design review overturned my "DPO -> acid" joint hypothesis (the acid yield is 7.5e-5; the acid in
pentose pots is cysteine's carboxyl carried out by the Strecker step r_dpo_ddp) and found that per-pot scales hid a
3-4 decade, pH-dependent 1-deoxypentosone level misfit; revised design recorded (prereg sec. 5).

**Does it serve the goal?** Yes, and the review changed what "next" means. The constant that couples sub-models 2-4
is `k_dpo_ddp` (cysteine Strecker by 1-deoxypentosone): it sets cysteine loss (route), the 1-DX level and MGO
(deoxyosones), the intact-C5 MFT route (thiols), and the pentose pots' acidity (sub-model 1's pentose block, which
its own parameters cannot fix). That is the single highest-leverage constant for MFT, the purpose's lead odorant.
Sub-model 3's step 0 (shared per-figure scales; an HA sink; Mikami's NF -> MGO share; the engine state named) comes
next, then its arms with `k_dpo_ddp` free -- after sub-model 1 has fixed the pH machinery they run under.

**Lesson recorded for myself:** two of tonight's design drafts (sub-model 1's fit, sub-model 3's hypothesis) were
wrong in ways an independent reviewer found in ~15 minutes each; the review step paid for itself both times. Keep it
before every fit.

## Checkpoint 3 -- 01:10

**What happened.** Sub-model 1 arm B finished (47 min) and FAILED in a telling way: all three pH slopes on their
1.5 bound, the acid yield on its bound, 6 / 10 data sets outside tolerance. Its diagnostics traced it to two flaws in
MY implementation of the review's fixes: DFG titrated with the calibrated pentose-ARP pKa (coupling the two blocks;
Davidek's NaOH diagnostic at -22 dex) and the trunk slopes on the label scale (importing an unbuffered DFG pH-stat's
own thermal shift into slopes measured in phosphate). Fixed with DFG's measured pKa (Mossine 1994, 8.18) and in-situ
slopes; the start's Davidek acids back to +0.16 dex. Recorded as Amendment 1 (prereg), the failed fit kept as
superseded, arm A stopped, both arms re-run (night queue: sub-model 1 B, A; then sub-model 3 A2, A3, A4, A1; ETA
~06:30-07:00). Sub-model 3 per-figure probe: the 1-deoxypentosone misfit across Zhai 2021's ladder is 2.05 dex RMS
once pots share a scale (0.41 when each pot had its own).

**Does it serve the goal?** Yes. The arm-B failure cost ~50 min of compute and found a real coupling error before it
could bias a shipped constant; this is the misfit protocol working. Risk to watch: if the re-run still puts the
slopes on a bound, the pH dependence of the trunk's acid supply is not a slope on these three steps, and the
mechanism step applies (e.g. phosphate catalysis, which Davidek's pH-5 water / phosphate contrast shows and the lane
lacks) -- pre-registered as an expected flag.

**Process note for the morning.** Three design flaws tonight (two reviewer-found, one fit-found), all mine, all
caught before shipping. The pattern: every "fix" I write gets the same scrutiny as the original design.

## Checkpoint 4 -- 03:10

**Results.** Sub-model 1 arm B (re-run, after Amendment 1): two starts converged; in-sample the hexose acids improve
(Martins 2005 FA / AA within tolerance), but two Amadori slopes sit on their 1.5 bound, a pH trend remains
(rho +0.40), leave-one-laboratory-out fails for every laboratory, and the pentose block is unchanged (its cause is
sub-model 3's k_dpo_ddp, per that review). FROZEN HOLD-OUT (Martins 2005b, committed before scoring): acetic acid
within tolerance in all five runs (incl. the constant-pH pH-stat, a mode the fit never saw); formic acid
under-predicted in all four (x10 at pH 5.5, x3 at 6.0, x2 at 7.5, x3.4 pH-stat) -- the over-steep formic slope.
Martins 2005b's own Table 3 puts the ladder's pH dependence on the glucose-consuming steps (condensation 0.50,
Glc -> Fru 0.88), which the lane leaves flat: Amendment 2 pre-registers arm C (those steps sloped), which predicts
the formic slope relaxes.

**Lost time.** Arm A hung for ~1 h on one pot (one worker at 100 % CPU, eight idle): a stiff corner of parameter
space. Fixed with a 120 s wall-clock limit per pot (a counted failure, tested); queue restarted: arm C, arm A,
then sub-model 3 A2, A3, A4, A1. ETA: sub-model 1 ~05:30, sub-model 3 A2 / A3 by ~07:30.

**Does it serve the goal?** Yes, and the hold-out is the first lab-independent evidence of the night: the acid
supply's pH dependence transfers to a mode the fit never saw for acetic acid and fails for formic acid, which names
the step to look at. Watch item: if arm C also leaves the formic slope on its bound, the 3-DG -> formic exit is not
a single pH-sloped step (Martins 2005b's Table 3 gives it pD 0.09 +/- 1.10 -- unidentified there too), and that
step goes to the literature.

## Checkpoint 5 -- 04:05

Arm C running for 1 h (12 constants: ~13 Jacobian points per iteration, all nine workers busy, no straggler since
the time limit). No new result; the plan stands (arm C's frozen hold-out decides whether the formic-slope problem is
the missing glucose-step slopes or the formic exit itself). Arm A and sub-model 3 follow; ETA for sub-model 3's
first two arms slips to ~08:00.

## Checkpoint 6 -- 05:05

**Results since 4.** Arm C (the glucose-consuming steps sloped): in-sample it did what Amendment 2 predicted
(Amadori slopes inside Martins' band, the pH trend gone, Davidek within tolerance, leave-one-out better) -- by putting
the condensation slope on its bound (3x Martins' 0.50). The frozen hold-out rejected it: the low-pH drifting runs lose
10-200x of their acid. Arm B transfers acetic acid to all five hold-out runs. Sub-model 1's done-test fails for both;
nothing ships; the results and the next round's requirements (acid amounts across pH as fit data, a fresh hold-out,
the formic exit to the literature) are in its prereg.

**Does the plan still serve the goal?** Yes, with one adjustment. The hold-out did exactly its job -- it stopped an
in-sample improvement that would have damaged predictions at low pH, the regime of acidic meat-flavour processing.
Arm A (the acetic-acid sink, 50 min in) finishes ~05:30; then sub-model 3's A2 and A3 (the k_dpo_ddp parent arm and
the 3-deoxypentosone route), which matter more for MFT than A4 / A1. Adjustment: after A3, run the scientific tier
and the gates (PR prerequisite) instead of A4 / A1, so the branch is reviewable when the owner wakes.

## Checkpoint 7 -- 09:10 (owner awake: "you decide on the priorities, be thoughtful")

**Results since 6.** The first joint deoxyosone + thiol fit (s042) closes the dicarbonyls. Its report printed raw
residuals before each semi-quantitative scale; read after the scale, the Jiangnan thiol figures fit in shape and the
misses are absolute: Schieberle under, Yiltirak over (CF-068, x120 between laboratories), the Hofmann pH ladder
inverted, 2-acetylthiazole x1000. Flux diagnosis: H2S only from (measured) cysteine thermolysis, the Strecker H2S
gone. The round-2 design review then found that round 1 keyed scales by ledger instead of by figure (an
implementation deviation from the reviewed design: all three round-1 runs repeated, queue running), that
2-acetylthiazole is not a leaf, and -- the big one -- that the lane's final pH in Zhou 2023's unbuffered ARP pots is
right at pH 6 and 1 / 2.2 units too high at pH 7 / 8. Probe: the lane makes the same ~3.5 mM net acid at every
initial pH and has NO pentose route to formic or acetic acid; the real pots make more acid at higher pH.

**Does the plan still serve the goal?** Re-ordered. The pH of near-neutral unbuffered pots is upstream of every
pH-dependent odorant step (thiols, thiazoles, pyrazines), and product formulations run at pH 5.5 - 7: the pentose acid
supply is now the gating item for the thiol pH arm, so it goes to the literature (search running) while the cheap,
independent thiol arms (R2b, then R2a + R2b) run on the corrected round 1. Duplicate supplements removed (verified
identical but for the download-date watermark).

## Checkpoint 8 -- 12:00 (owner: "take a step back and review where we are")

**Against the north star** (leave-one-lab-out error, interval coverage, evidence grounding, blind tests): nothing
has shipped and no lane has a passing leave-one-lab-out; coverage is a Phase 3 number; grounding improved (Desclaux's
measured enolisations, Davidek's isotope-dilution acids, the frozen-factor tests); WUR still has no data. The gains
of the last day are diagnostic and structural, not scored accuracy.

**What moved:** the measurement model (scales per figure; residuals read after the scale) took the joint fit from
chi2 21.8 to 5.4 and turned two "failures" into fits; round 1 is robust across three runs; MFT decay is x8 slower
than shipped (CF-013 independently); 2-acetylthiazole's constant is fixed by its own data (x10^-3). Four structural
gaps found, each by evidence, not by score: no pentose + amino-acid entry (built, switch off); a fitted acetic-acid
sink that the isotope-dilution data rule out (CF-069); the pH of unbuffered near-neutral pots (no pentose acid
route); and a missing EARLY MFT source (Schieberle's MFT is late, not low).

**What worries me:** (1) the sub-model sequence has eroded into a growing joint fit (8 -> 15 coordinates), the
pattern the strategy's diagnosis 3 warns about; the coupling is real (the 1-deoxypentosone pool feeds norfuraneol,
MGO and acetic acid), so the guard is not to un-join it but to require leave-one-paper-out and profile checks before
any acceptance. (2) Between-laboratory spread for absolute thiols is 10-100x (CF-014, CF-068): that bounds what
"lab-independent" can mean for MFT / FFT, and the honest product may be the spread itself. (3) Sub-model 1 failed and
the sequence went on (owner's call to run the joint fits); today's findings are exactly the foundation it was meant
to give, so the pentose entry is the right next step. (4) The other odorant lanes (Strecker aldehydes, pyrazines,
lipid aldehydes) have had no Phase 2 time.

## Checkpoint 9 -- night of 2026-10-30 / 31, ~00:40 (owner: "keep the big picture; every hour remind yourself")

**Goal (unchanged):** the right mechanisms for the alternative-protein odorants, with honest lab-independent
uncertainty; the shipped model loses to a typical-yield baseline out of sample (STATUS_AND_PLAN sec. 2), and the
architecture is split into lanes that do not compose (sec. 4) -- the composable network is now plan step A2, before
new lanes. **State:** branch 4810355 + uncommitted R3T arms / prereg sec. 12; PR #30 push waits on the freshness gate;
R3Es fitting; R3T under design review; two readers keying Hofmann 2000 (Strecker anchor); experiment register 21
entries, gated. **Decided this hour:** round 3's failure is a temperature conflict, not a pH one (EX-016); R3T frees
the norfuraneol -> MFT barrier (k_nf_mft is referenced at 145 C, so the barrier sets its 95 C value); the dry check
lowers its expected pass to 0.2. **Is this the most valuable next step?** Yes: MFT is the alt-protein meaty odorant
and its route is the sulfur lane's biggest mechanistic contradiction; the Strecker keying is data work that survives
any architecture.

## Checkpoint 10 -- 2026-10-31 ~00:45 (after a context compaction; owner: "remember the state of the plan, the priorities")

**Goal (unchanged):** the right mechanisms for the alt-protein odorants with honest lab-independent uncertainty
(STATUS_AND_PLAN sec. 1); the one plan is STATUS_AND_PLAN sec. 5. **Two independent reads changed tonight's order:**
- *Premortem* (STATUS_AND_PLAN sec. 5 "Rules adopted"): every joint-fit arm (E, R3E, R3En; R3Es too) stopped on its
  budget -- 200 evaluations = 11 optimiser iterations -- so EX-015 / EX-016 are re-marked inconclusive; the "MFT
  loss x8 slower" and "actz x10^-3" findings are provisional (R3En's optimum restores the acetic sink and the shipped
  decay). New rules: no verdict from a budget-exhausted fit; rivals under ONE fixed measurement model; no silently
  lost signature; signatures frozen v1 by mechanism group with a chance baseline; section B time-boxed to 48 h.
- *R3T design review* (thiol prereg sec. 13): R3T withdrawn. The binding conflict is Cerny 95 C vs Schieberle 100 C
  (same chemistry, 5 K apart): the DDP route kills Schieberle's MFT lag, the NF route makes fed NF dominate Cerny's pot.
  The question is now mechanistic -- what makes MFT late at 100 C -- and a literature pass is running (DOIs to the owner).
**Running:** R3Es (budget-limited, scored for the record); then the convergence queue (scratchpad converge.sh: E
converged, E under its own fixed noise, then R3En / R3Es / R3E from their optima under E's noise; ftol 1e-4,
9 workers, ~5 - 8 h); Hofmann 2000 adjudication (agent). **Priorities after the queue:** (1) read the converged
rivals -> round-3 verdict; (2) pre-register the lag rivals from the literature pass; (3) signatures v1 frozen with
groups and chance baseline; (4) leave-one-lab-out on the converged finalist; (5) A3 config context manager, then A2.
**Is this the most valuable next step?** Yes: every verdict so far rests on unconverged fits; converging them is
cheap compute and protects every later decision.

## Checkpoint 11 -- 2026-10-31 ~01:00

**Goal:** unchanged (sec. 1). **State:** commits 5096e72 (premortem folded in) and 8a25ae7 (signatures v2, R3Es,
Hofmann 2000, T1) on submodel1-ph, verifying on a clean checkout before the PR #30 push. Convergence queue running
since 00:55 (E conv -> E fixed -> R3En / R3Es / R3E fixed). Chance baseline (v2, R3En box, 200 draws) running.
**Learned this hour:** the lag literature pass turns the MFT question into a supply-rate one: Schieberle's MFT is
near-linear after a ~9 min dead time and stops at ~6 h (a reservoir), and T1 shows start conditions explain ~0.5 dex
of the early excess at 30 min but not at 60 min -- the intact-C5 route's early supply is ~10x too fast. Rivals R4a /
R4b drafted (prereg sec. 14), to be reviewed in the day. R3Es: budget-exhausted, same picture as R3En.
**Is this the most valuable next step?** Yes: the converged rivals decide whether round 3 said anything; R4 is the
mechanistic answer either way, so its design waits for the owner's review, not for compute.

## Checkpoint 12 -- 2026-10-31 01:40

**Goal:** unchanged. **State:** branch 9 commits ahead of origin (last 1d201d5); full verification of 1d201d5 running
on a clean checkout (the first run of 8a25ae7 found the READMEs stale -- regenerated before staging; lesson added).
E convergence fit running (since 00:55); R3En / R3Es / R3E under E's noise follow; R4f (3 profile fits) queued
after; chance baseline (v2) finishing. Register: 30 entries, every finished run in. **Finished / decided since 01:00:**
R3Es scored (inconclusive, budget-exhausted); the early intact-C5 supply traced to r_ttca_deg (EX-028); R4a reviewed
and replaced by R4f (EX-029 / EX-030); B7 redefined as the 1-deoxypentosone partition + the cysteine sink, ahead of
any further MFT-route arm. The cron's step list is stale (R3T withdrawn, Hofmann done): the live list is
STATUS_AND_PLAN sec. 5. **Most valuable next step:** the B7 design (no compute needed to draft): the cysteine sink and
the DPO partition are upstream of every thiol, so getting them right is worth more than any route arm.

## Checkpoint 13 -- 2026-10-31 02:40

**Goal:** unchanged. **Finished since 01:40:** converged E (status 2) and E under fixed noise (the reference, cost
1141); R3En fixed (1076: the data accept DDP formation, 65 units for one coordinate -- EX-016's reading survives
convergence); R3Es fixed (1069: the pH slope ties zero -> R3En admitted); Cerny pot D fails in every arm (0.98 - 0.99).
Converged E reverses "MFT loss ~8x slower" (decay back at shipped) -- STATUS sec. 3 corrected. Chance baseline: a
median random draw passes 1 of 11 groups, P(>= 3) 0.22, so 2 - 3 groups is not evidence; credit goes to the <= 5 %
groups. B7 groundwork + review: B7a (partition: two-level measurement model, split-law and TDP-exit arms) inside the
box, B7b (cysteine fate, buffer catalysis) after a cost-by-block check; CF-071 registered. Verification of 1d201d5:
gates PASS, unit 1150 PASS, scientific 108 PASS, freshness running. **Next:** R3E fixed -> EX-024 verdict; push when
freshness passes; R4f queued; B7a code after the queue (fitting.py is imported by every queued fit).
**Most valuable?** Yes. The cheapest high-value work left tonight is reading results as they land and keeping the
record exact; new code waits for the queue.

## Checkpoint 14 -- 2026-10-31 03:41

**Goal:** unchanged. **Finished since 02:40:** R3E fixed (1140 = E: the shipped law switches the route off) closes
EX-024 -- R3En admitted, the round-3 readings survive convergence. 1d201d5 pushed to PR #30 after the full clean
verification; PR description updated. R4f (EX-030): TTCA's direct step stays at its value (+0.2 +- 0.14), the NF ->
MFT barrier is unidentified (cost flat within 4.8), Cerny fails at every in-class barrier -- rejected. B7a: the two-level
measurement model coded, reviewed (algebra checked), 48 tests, bit-identical when off; cherry-picked; its scorer shows
the admitted structure failing the scale-free partition signatures in time (1 / 8) and in pH (wrong sign, 1 / 6)
(EX-033). **Running:** the B7a null under the two-level model (w 0.2 then 0.4); the R-S / R-T / R-ST arms being coded
in a separate worktree, their design (prereg sec. 8) under independent review; df6b80f verifying (6 commits unpushed).
**Most valuable?** Yes: the partition is now shown structurally wrong by tests a scale cannot absorb, and it is
upstream of every thiol. Morning: a summary for the owner (decisions: PR #30 merge; the stranded main-checkout fix;
the designed-experiment brief; the four DOIs).

## Checkpoint 15 -- 2026-10-31 04:40

**Goal:** unchanged. **Since 03:41:** df6b80f verified clean (gates, unit 1150, scientific 108, freshness) and pushed;
B7a arms implemented and code-reviewed (nothing blocking), starts moved near the null and declared (prereg sec. 10),
the scorer gained a chance mode; the B7a arms' design review (sec. 9) replaced R-S's linear law with a curved one,
R-T's constant with the trunk's own, and added B7D (a long-lived DPO with a second MGO parent) after the observation
that TDP and DPO decline at the same rate in Zhai 2023 JAFC. B7a null (w 0.2) converged: partition signatures
unchanged (1 / 8, 1 / 6, 1 / 2) -- the structure, not the measurement model, is wrong (EX-035). **Next:** w 0.4 null
(running) -> cherry-pick the arm commits -> regenerate engine-dependent artifacts -> chance baseline -> the four arms.
**Most valuable?** Yes. The morning summary is the other must-do before ~08:00.

## Checkpoint 16 -- 2026-10-31 05:40

**Goal:** unchanged. **Since 04:40:** B7a null w 0.4 = w 0.2 (the conclusion does not depend on the response-factor
width); the B7a arm commits cherry-picked; every engine-dependent artifact re-hashed (all signature counts unchanged:
the defaults are bit-identical); B7D fitted -- r_tdp_mgo identified but small (-6.75 +- 0.19), c2c3 not released
(+2.37), partition signatures 2 / 8, 0 / 6, 1 / 2: the long-lived-DPO rival is not chosen through this MGO parent.
**Running:** B7S (since 05:23), then B7T, B7ST; then the null under its own fixed noise and both chance baselines.
10 commits ahead of PR #30; one full verification after the arms. **Most valuable?** Yes; the next judgement is
whether ANY B7a structure moves the partition signatures -- if none does, the partition's fault is in a step none of
the arms touches (the deoxyosones' detection / the quinoxaline response, or a TTCA -> deoxyosone route with its own
split), and that becomes the owner's morning question.

## Checkpoint 17 -- 2026-10-31 06:40

**Goal:** unchanged. **Since 05:40:** the four B7a arms converged and none passes the partition signatures (EX-036);
the null under the same fixed noise costs 1093: B7S -54 / 2 coordinates and B7D -19 / 1 are real gains without the
signature shape, B7T +53 rejected. Flux attribution at B7S's optimum: DPO's base-favoured beta-cleavage (fitted on
Davidek's glycine pots) takes 51 - 61 % of DPO at pH 8 and cancels the split law in amine-free TTCA pots; 39 - 45 % of
TDP is pH-flat caramelisation. **Running:** the B7a chance baseline, then the thiol v2 chance baseline. **Next:** full
verification of HEAD, push PR #30, morning summary. **Most valuable?** Yes: close the night with a verified, pushed
record and a clear next question (the amine-assisted beta-cleavage), not another unreviewed arm.

## Checkpoint 18 -- 2026-10-31 07:40

**Goal:** unchanged. **Since 06:40:** both chance baselines regenerated (thiol v2 identical; B7a partition thresholds
strict, chance <= 2 %, the fitted arms at or below chance -- EX-037). Verification of 8339590: gates, scientific
108 / 108 PASS; one unit test and the freshness gate failed on a stale paper registry (citation line numbers moved by
the register edits; lesson added); regenerated in 3d606a3, re-verifying (gates + freshness, unit in parallel). PR #30's
pushed head df6b80f is CI-green. **Next:** push 3d606a3; the morning summary below; stop the hourly cron when the
owner is back. **Most valuable?** Yes -- the night's record verified and pushed before the owner reads it.

## Morning summary -- 2026-10-31 ~07:45 (for the owner)

**The night in one paragraph.** The premortem was right about the biggest risk: every joint-fit arm had stopped on its
evaluation budget, so round 3 was re-run converged and under one fixed measurement model. The converged round-3
verdict stands (the data accept DDP formation; R3En admitted), but one headline finding did not survive ("MFT loss ~8x
slower" -- converged E puts MFT decay back at the shipped value). Every attempt to reconcile Cerny 2003's isotope
labelling failed (R3T withdrawn on review; R4f rejected; no in-class barrier passes), and the trail led upstream: the
1- / 3-deoxypentosone partition is structurally wrong in time and in pH, by scale-free tests a calibration scale
cannot absorb -- the fitted structure scores AT OR BELOW random constants on them (EX-033, EX-037). Four reviewed B7a
arms (split pH law, TDP exit, both, a second MGO parent) improve the cost but none reproduces the partition; the flux
check says why: DPO's base-favoured beta-cleavage to acetic acid, fitted on Davidek's GLYCINE pots, drains DPO exactly
where the split would make more of it in amine-free TTCA pots.

**What is now true (with register ids).** Converged round 3: R3En admitted (EX-024). R4f rejected (EX-030). Chance
baselines: thiol signatures -- a random draw passes a median 1 of 11 groups, so 2 - 3 is not evidence (EX-025);
partition signatures -- strict tests, chance <= 2 % (EX-037). Two-level measurement model in `fitting.py` (a shared
figure scale + a response-factor term; bit-identical when off; 48 tests). Hofmann 2000 (Strecker anchor) keyed by two
readers + adjudication, 111 rows (EX-026; CF-070). Reviews this night: premortem, R3T, R4a, B7 draft, B7a arms, B7a
code -- each changed the design before any compute was spent on it.

**Decisions for you.**
1. Merge PR #30 (head 3d606a3 verified locally: 12 gates, unit 1179, scientific 108, freshness; CI re-runs on push; nothing shipped moves). P(regret) low.
2. The stranded fix in the main checkout (uncommitted since 10-28: the `k_odg_da` band, a pH fixed-point warning):
   rescue onto a branch off main? It touches your working copy, so it waits for your OK.
3. Send the one-page designed-experiment brief (premortem #9): designed data is the only lever on the x10 - x1000
   between-lab spread. A time-resolved pot-D experiment (fed NF + [13C5]ribose, 0.5 - 6 h at 100 C) would settle the
   MFT route question directly.
4. Papers (none blocks today's work): Parker et al. 2014 https://doi.org/10.1016/B978-0-12-398549-1.00040-4 (ribose
   thiazolidine limits H2S); Roberts et al. 1987 https://doi.org/10.1021/jm00393a034 (ribose-cysteine thiazolidine
   ring opening); Tai & Ho 1997 https://doi.org/10.1021/jf970297t and 1998 https://doi.org/10.1021/jf971111t
   (cysteine / cystine, thiol oxidation). For the amine-assisted beta-cleavage question, Davidek's own mechanism papers
   are on disk; I will list any gap with DOIs once the design is drafted.

**Proposed next (section B time-box closes 11-02 00:30).** B7a-2: DPO's beta-cleavage amine-assisted (its base factor
times a free-amine term; one structural change, no new coordinate) -- design, independent review, then fit; if the
partition still does not move, section B stops by its rule and A2 (the composable network) starts, with the thiols
shipped on the measured between-lab spread. P(B7a-2 passes (ii) by segment) ~0.3.
