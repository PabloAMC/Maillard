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
