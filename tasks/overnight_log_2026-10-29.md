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
