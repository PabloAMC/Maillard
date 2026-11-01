# Red team: thiol loss (a) and cysteine fate (b), 2026-10-31

Read-only review of scratchpad/thiol/{prereg.yaml, mft_loss.py, cys_fate.py, potD_sensitivity.py, *_results.json, *_run.log}
against the dossiers. Re-ran potD_sensitivity.py (reproduces: loss law moves pot D share p05-p95 -0.39..+0.26 dex).
Ran one extra check inline (no file): patent Fig. 1 shape under pre-added IS, below.

## 1. Transcription (all match)
- us10076126 sec. 3a: "0.5 | 15", "1 | 124", "2 | 575", "4 | 1855" ng/g.
- vanseeventer2001 Table 1 MFT row: "59 | > 90 | > 90 | 69 | 42". Code: entry 4 lam = -ln(0.31)/1440 = 8.13e-4 OK.
- schieberle2000 Table IV: "4.5 | 13.8 | 156 | 179". liu2023b Table 2: "1252.67 / 760.33 / 628.00".
- cerny2003: "the norfuraneol route supplies 7 % of the MFT". Zhai 2021 Fig. 4b / 5a / 5b / 5c arrays match the dossier
  tables to the digit; Fig. 6 Cys/PENT from route_submodel.yml.

## 2. The spike bound
- Patent text (pdftotext, col. 8): precursors weighed "together with an internal reference of the labeled [2H3]-MFT";
  then heated. Pre-heating addition is real. That the plotted quantities are IS ratios, and that the IS was still
  detectable at 4 h, is NOT printed ("The quantities of MFT and FFT are plotted"); IS amounts, LOD not printed.
  The threshold band [1e-3, 0.1] is a judgement. CF-013 already asserted "at least 5x too fast" before this prereg.
- Stronger, unused evidence: with pre-added IS the reported curve is R(t) = int s(t') exp(Lambda(t')) dt'.
  Fitting the 4 points (shared offset, max |resid|) with the prereg's supply family (lag 0-30, power 0-3,
  exhaustion 2-30 h): lam 0 -> 0.006 dex; 0.01 -> 0.036; 0.015 -> 0.061; 0.02 -> 0.12; 0.03 -> 0.38; 0.095 -> 3.9.
  With exhaustion down to 20 min: lam 0.03 -> 0.012, 0.05 -> 0.037, 0.095 -> 1.6. So A0 is rejected by the shape
  under any supply, the bound is lam(110 C) <~ 0.02 (exhaustion >= 2 h) or <~ 0.05 (fast exhaustion).
  Conditional on IS-ratio quantification (P ~0.85) and one patent with no replicates.
- Matrix: 0.25 M Cys. van Seeventer entries 2 vs 4 say Cys slows loss >= 2x, so the bound is for a Cys-rich pot.

## 3. van Seeventer 50 C as a test at 95-145 C
- Different system: 50 uM thiol in clean 0.5 M phosphate (entries 2-5), or a 5 uM process flavour (entry 1), 24 h,
  zero-order regressions with no R or SE. The repo's own PROHIBITED_DERIVATIONS: "the processing sink and the storage
  sink are DIFFERENT CHANNELS"; the dossier: "its own sink is expressly not the one operating at 115-120 C".
- Ta2 is not a fair test: A3's lam0 and K_c are set from entries 2 and 4 and then scored on entries 2-5 (double use);
  A1 is anchored on entries 1/4 and is pot-independent, so it fails entries 2/3 by construction; A2 is zero in buffer
  by definition. "Power 0.95-1.0" measures the construction. Entry 1 is not scored at all (A3 would miss it).
- A0 fails Ta2 for being too SLOW at 50 C (median 31 %/day vs 59-69). That is the opposite direction from Ta1, so
  "A0 infeasible on two labs" does not support "too fast"; together they only say Ea_eff < 101 IF one channel spans
  50-110 C, which is the prohibited assumption.
- Ea-dependence: lam(110) bound is Ea-free. Extrapolated bounds are not: lam(145) <= 0.02 x exp(Ea/R dT) = 0.13
  (Ea 57) ... 0.29 (Ea 101) ... 0.66 (Ea 133); vs engine 1.38 the factor is 2-10x. ">= 3x at 145 C" holds only if
  Ea <~ 105-115.

## 4. Lab independence
- Grouping right. Liu 2023b is CAAS (Zhang Dequan), not BTBU. van Seeventer (Quest) uses Hofmann & Schieberle's recipe:
  independent lab, shared pot design. "Too fast" rests on ONE lab (Nestle patent); Ta2 points the other way.
- (b): all tests Zhai 2021 (one lab, one method, HPLC-ELSD); K from Cao 2017 (BTBU). Cao itself contradicts B1 (below).

## 5. Power
- Synthetic truths are drawn from the same bands and scored with the same likelihood: power to separate bands,
  not robustness to misspecification (in-situ pH, ELSD non-linearity, Cys-dependent law).
- Supply power band in code (0.0, 1.5) vs prereg 0.5-1.5; Ta5 nuisance (lag 2-8, exhaustion 5-120 min) not in prereg.
  Ta5 feasibility <= 0.10 for every rival incl. A0: the 168 C supply family is too narrow to fit a decelerating decline.
- Ta1 power: A1|A3 0.55, A0|A1 0.59: only A0 vs NC-a / A2 is powered.
- (b) Tb1 noise 0.07 dex is reasonable (Cys CV 1-17 %). Tb2 changed from the prereg's yield ratio to a level with a
  shared 0.2 dex offset; B1 over-predicts by ~0.45 dex at every pH and passes only via that offset (feasibility 0.15).
- B1's k_f upper edge 0.17 L/mol/min was derived from the SAME xylose contrast (cysteine review: "So at 100 C,
  k_f <~ 0.17 L/mol/min (mine)"). Tb1 is partly double use; the decade between B1's edge and B0 (1.7) was never sampled.

## 6. Unlisted structures
- (a) non-monotone Cys law (van Seeventer Fig. 1: "a clear minimum at ... around 50 mM", figure-only, never digitised)
  plus sugar protection (entry 5 slower than 4; A3 predicts faster). Reversible sequestration (MFT-S-S-Cys, Guth's
  disulfide reverting 62 %) not listed: it makes "loss" pot- and analysis-dependent.
- Shared (a)+(b): a thiol/disulfide redox buffer with sugar-derived reductones: Cys and reductones protect MFT
  (entries 4, 5; spike survival at 0.25 M Cys) and keep Cys reduced (+Xyl raises free Cys, Tb1). Prereg splits it into
  A3 and B4 with unshared constants, and A3 has no reductone term.
- (b) Tb3 xylose sink: at pH 8 +Xyl, xylose falls 6.39 -> 2.99 mM over 80-120 min while TTCA is < LOD, so this xylose
  is not consumed via TTCA (it would need TTCA turnover > 8/min). It is a non-Cys sink, or pH confounding: the unbuffered
  xylose control loses the same at every label pH (9.53 / 9.40 / 9.40 / 9.36), consistent with its pH collapsing,
  while TTCA/Cys buffer the test pot. Code uses the control's k_x and the label pH (q. 8). So Tb3 neither supports
  B2/B0 nor refutes slow exchange, but "fails all" should be read as a measurement-model problem.

## 7. Do the scores follow?
- (a) Yes for A0 on Ta1 (the shape argument is stronger still). Ta2 "no structure passes" is by construction.
- (b) Tb1+Tb2: slow family +14.3..+14.9 vs B0 -6.1: correct. Tb1+Tb2+Tb3 joint: B0 -32.08, B1 -32.37, B4 -31.93: a
  tie. "B0 rejected" depends on dropping Tb3 (defensible per sec. 6, but post hoc), plus the TTCA gate (0/300).
- Zhai 2019 (same lab), Tb4 unscored: TTCA 14.5 mM already at t = 0 after ~10 min at RT (82.6 Cys + 165 Xyl, pH 7.4),
  so k_f(RT) >= ~1e-4 L/mmol/min; 70 C, 0-20 min: 1.33 mM/min / (~55 x ~140) ~ 1.7e-4. Cao 2017 relaxation:
  k_f >= ~4.4e-4 L/mmol/min at 90 C. Both exceed B1's 100 C ceiling 1.7e-4. B1 is plain mass action, so its
  "unspecified concentration law" (the prereg's reason for not scoring Tb4) does not exist. Slow, uncatalysed,
  pH-flat exchange is contradicted by two labs. What survives is "the engine's single pH-flat k_r/K can't hold
  in both the dilute and the concentrated pots".

## 8. Code
- Units consistent (min, mM, L/mmol; %/day over 1440 min; log10 vs ln is handled: Ta1 divides by ln 10).
- Label pH used for k_cd, Kallen K, k_x in unbuffered pots (declared, but it drives Tb2 pH slope and Tb3).
- k_cd lower edge 3.5e-3 vs the prereg's 4.5e-3; A3 entry-5 free-Cys nuisance [20, 150] mM can't go below 20 (at
  K 60 with 250 mM ribose, equilibrium free Cys is ~9 mM), same direction, so no change to the verdict.
- No bug that flips a verdict.
