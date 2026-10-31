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
- Findings that hold: acetic acid transfers across laboratories, formic acid does not; methylglyoxal fits; MFT's
  loss is ~8x slower than shipped (CF-013 independently); 2-acetylthiazole's constant was 1000x high; the
  pentose + amino-acid reaction was missing and is now built and identified (Cao hold-out pH within 0.06); the
  fitted acetic-acid sink is ruled out by isotope data (CF-069); MFT's route in the lane contradicts isotope
  labelling (Cerny 2003), and the route mix changes with temperature (Cerny 95 C vs Liu 168 C); two between-lab
  conflicts recorded (CF-068: x120 for the same pot).
- Mechanism signatures: the round-2 model 1 / 15 (despite reduced chi-square ~5); with the pentose entry 5 / 15.

## 4. The plan

**A. Consolidate and de-risk (this week; ~1.5 days of work).**
1. Merge PR #30 once CI is green (owner's call; nothing shipped moves). From then on, one PR per sub-model, code
   separate from regenerated results.
2. Rescue the stranded fix in the main checkout (uncommitted since 2026-10-28: the `k_odg_da` band excluding the
   shipped value, a pH fixed-point warning) onto a branch off main (owner's OK; it touches the user's working copy).
3. Engine configuration made explicit: one context manager for every fit's switches (restores on exit), then a
   frozen config passed per job, and a persistent worker pool instead of forking per evaluation (the macOS
   fork-safety hazard and the per-evaluation cost go away). 8 - 10 h.
4. One fit harness in `src/` for the sub-model fits (three drifted copies today). 3 h.

**B. The sulfur lane's mechanisms (1 - 2 days of work + overnight compute).**
5. Round 3: the pH law of the intact-C5 MFT route (rivals running; the shipped free-amine law already failed).
6. The route mix vs temperature: separate barriers for the DDP, norfuraneol and C2 + C3 routes, tested on Cerny
   95 C, Hofmann 145 C, Liu 168 C.
7. The glycine-system partition (1-deoxypentosone over-feeds norfuraneol, 25 - 36 mol% vs <= 5 %); the pH model on
   unbuffered pots with the pentose acids in; then the thiol pH arm (HS- / neutral ratio, furfural's pH direction).
8. Adequacy: leave-one-laboratory-out with the between-paper error on the finalist (6 - 10 h compute; overnight or
   Colab), then Phase 3: one regeneration shipping every switch that passes its evidence audit (B59b's de Roos route,
   the entry, the corrections).
   *Success:* >= 11 / 15 signatures and Cerny's route tests; out-of-sample coverage 0.80 - 0.97; beat both null
   baselines out of sample. P ~ 0.4.
   *Stop rule:* if signatures and coverage plateau after 5-7, ship the thiols with the measured between-lab spread as
   their uncertainty and go to C.

**C. The lanes alt-protein flavour needs most (each ~2 - 3 days).**
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

**Calibration.** P(the sulfur lane beats both null baselines out of sample within ~2 weeks) ~0.4; P(a useful
Strecker-aldehyde lane within ~1 week of starting it) ~0.5; P(absolute odorant predictions within x3 across
laboratories for most rows) ~0.15 -- the realistic near-term product is ranked, relative predictions with honest
intervals (~0.55).
