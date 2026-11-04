# H0 (EX-059): independent results review

- Date: 2026-11-04. Reviewer: an independent agent with no part in the design, the scorer or the run. Read-only.
- Scope: the report `kinetic_core_thiol_lolo_nested_tau.{json,md}` against the prereg (v3 + sec. 14), amendments A1
  and A2, and the controls.

## Verdict on the report: no blocking problem

- The report applies the prereg's rules as written (D1 floors, D2 branch logic and test inversion, D3 criteria,
  precedence, set combination, FRAGILE, owner notes).
- Every key number reproduces to < 1e-9 from the fold and pair records with an estimator written from sec. 2's text
  (not the scorer's functions): tau_ABS(all) 0.756183 (warm) / 0.963970 (reset); omega(all) 0.4187 / 0.3901;
  tau_FIG(all) 1.2337 / 1.2354; every outer lab's nested omega(S), tau_ABS(S), tau_FIG; (b)'s 90 % coverage
  0.842 / 0.906 (warm, lab-weighted / pooled) and 0.806 / 0.728 (reset).
- A1 and A2 changed no verdict rule; the commit order is freeze 1fb63f5a, A1 f0352833, A2 9513f252, records
  e41c6995 and fb4ba613, scoring at fb4ba613 on a clean tree.

## How the result must be read (MAJOR, interpretation)

1. **The branch that ships (relative, within-lab) is not better than a flat line.** Within-group shape LPD, engine
   minus null, is negative at 4 of 5 labs in both sets; in the reset set the null is better by the D3 rule (-0.204
   nats/row, SE 0.095, p = 2/32). The two nulls are identical on this branch ("no change within a group").
2. **"Relative" protects only within one level group** (paper x species for absolute rows; paper x figure x species
   for figure series). Cross-species or cross-figure ratios carry the level error (~1 - 1.7 dex). A two-row ratio in
   one group still carries 2 omega^2: ~0.6 dex SD, a 90 % interval of about x/÷ 10.
3. **tau_ABS pools species; the thiol evidence is one laboratory.** DFA holds the only held-out absolute MFT / FFT
   rows (4 groups); Jiangnan's and Ochanomizu's absolute rows are mostly MGO, erring in opposite directions:
   species x laboratory structure, not one laboratory effect (for H1).
4. **tau_FIG(all) = 1.23 dex** (drop-one 0.96 - 1.38), about 4x the shipped declared figure width (0.30): figure-
   series levels are also ~1 dex off for a new laboratory. D2 governs absolute rows only.
5. **"CALIBRATED" means "cleared the 0.70 floor".** In the reset set three of five labs cover 0.63 - 0.72 at nominal
   0.90; read (b)'s 90 % intervals as ~75 - 80 %. The Student-t sensitivity comes closer to nominal.

## MINOR and NOTE

- The warm set alone gives D2 branch (ii) (0.756 <= 0.8); the larger set (reset, 0.964) gives (i) as the prereg
  says. Disclosed here. U is censored at the grid ceiling (U > 1.60); irrelevant to branch (i).
- UNSTABLE fired only in warm, by crossing 0.8 (drop Jiangnan: 0.756 -> 0.851, ratio 1.2); uninformative (P5 0.63).
- The refusal is robust: without Jiangnan, DFA + Ochanomizu give 0.85 / 0.96; Paule-Mandel 0.97; row-weighted 1.09;
  L = 0.60. Reviewer's estimate: P(true tau_ABS > 0.8 for the shipped model on a new lab) ~0.75. For the relative
  branch the declared leaks bias omega low, so its intervals are not conservative.
- The engine's best-supported skill is on the refused branch: on absolute rows it beats typical yield at all three
  absolute labs (p = 1/8, the minimum) and no change at all three in the reset set (nulls' tau_ABS 1.05 / 1.35).
- All D3 comparisons are "driven by few rows"; DFA's interval-score differences dominate the IS means.
- Report cosmetics to fix in the generator later: print "FRAGILE: no"; print the per-set D2 branch; "U > 1.60
  (grid ceiling)"; mark Reading's tau_FIG 0.300 as a fallback in the md.
- A2 auditability: the replaced records' sha256 are now in A2; the re-runs started 2 minutes before A2's rule was
  committed (17:31 vs 17:33); the seven clean records were committed in one batch.

## Plain-language summary

H0 measured, without letting each lab see its own data, how far off the thiol lane is for a laboratory it has never
seen. Absolute concentrations are off by about 1 dex (a factor ~10) typically, and the data rule out less than ~0.6
dex, so absolute predictions for a new lab are rightly refused. The engine does beat chemistry-free guesses on absolute
values at the three labs that measured them, but modestly, and for the thiols themselves that rests on one lab. What
the rule lets ship is relative prediction within one study series, and there the engine was not better than assuming
nothing changes; in one set it was slightly worse. A 90 % interval on a within-series ratio spans about x/÷ 10, and
nominal 90 % intervals cover ~75 - 80 %. With 5 labs (3 for absolute values) nothing here can show the engine beats a
null. Provisional until H0' brings two or more new laboratories.
