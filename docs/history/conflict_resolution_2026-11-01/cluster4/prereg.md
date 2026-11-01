# Cluster 4 pre-registration: observation operators and vessel O2 in the thiol and Strecker levels (2026-11-01)

Written before `cluster4_tests.py` was run. Honesty note: the source numbers were READ during scoping (dossiers,
register, fit tables, the s042 / R3En residual tables), and some back-of-envelope ratios were done by hand while
reading. The rules and expectations below are fixed before the scripted computation; nothing is tuned after it.

Conflicts: CF-068, CF-073, CF-031, CF-016, CF-015, the cause half of CF-014, the ingest findings (Grosch 1992 SIDA
timing, Mottram 1996b SDE reduction, MFT loss pot dependence, Weenen 1992 donor localisation, pyrazine R across
Shibamoto / Leahy / Ghent), and EX-050's verdict.

Bands come from data other than the scored contrast; every on-disk data set bearing on a question is scored or its
exclusion justified (sec. 9).

## 1. CF-068 (Schieberle 2000 vs Yiltirak 2026 MFT/FFT at 100 C, ~x130 in MFT)

Rivals (chemistry and operator, stated neutrally):
- R68a **Vessel O2.** 17 mL air over 3 mL (O2 : Cys 2.0; O2 : (Cys + pentose) ~1.0) keeps Yiltirak's pot oxidising:
  cysteine to cystine (metal-catalysed, Ehrenberg 1989 form), MFT to disulfides / oxidative adducts. DFA's autoclave
  (100 mL in 200 mL, O2 : reductant ~0.065) is O2-starved after the reductones consume it.
- R68b **Free-thiol operator.** Both are SIDA with the standard added after heating, so both count free thiol only;
  differences in the bound (disulfide / adduct) share then read as level differences. (Does not separate the two by
  itself; needs R68a or chemistry to set the bound share.)
- R68c **DFA high.** Garching's 100 C series over-reads (CF-055 internal inconsistency; no SD, no n).
- R68d **Yiltirak work-up low.** PFBBr derivatisation of a 30 uL aliquot in 1 M NaOH; matrix-matched calibration
  from the unheated EMULSION used for buffer samples.
- R68e **Mass action + pH alone** (null): the 4x pentose, 1.3x Cys, 0.5 pH unit difference explain the gap.

Test T68 (between-lab level residuals, engine-free): residual of each absolute-quantified (isotope standard)
pentose + cysteine MFT level against the DFA 100 C anchor (Schieberle 2000), adjusted for (i) concentration, three
variants: first order in cysteine, first order in pentose, product; (ii) pH, 0.45 dex / pH unit above pH 5 (Hofmann
1998 T2, pH 5 -> 7; a DFA within-study ladder, not the scored pots); (iii) temperature, formation Ea band 80 - 130
kJ/mol; (iv) time, t^1.35 (Schieberle 60 -> 360 min; the anchor's own shape, used only for the anchor's 4 h
interpolation); (v) xylose / ribose 0.72 (Hofmann 1998 T1, 145 C). Pots: Yiltirak 100 C / 4 h, US 8,920,862 buffer
(95 C / 6 h), US 10,076,126 buffer (110 C / 4 h, standards co-heated: residual bracketed for lambda 0 - 0.01 /min).
Covariates: O2 : Cys, O2 : (Cys + pentose), IS timing, enrichment method, concentration, lab.
Decision rule (triage): one covariate explains >= 70 % of the between-paper log variance AND every within-study
contrast agrees in sign -> EXPLAINED by that covariate; otherwise the spread is a declared between-lab band. With
n <= 5 every regression is EXPLORATORY (reported with df and a permutation p, not used to move a constant).
Negative control: the same normalisation applied to Hofmann 1998 145 C ribose vs xylose (same lab, same vessel):
must give |residual| <= 0.2 dex after the xylose factor (it is how the factor was set, so a pass is a check of the
arithmetic, not evidence).

Expectations: mass action + pH explain <= 1.1 dex of the ~2.1 dex gap (P 0.8); residual Yiltirak <= -1.0 dex under
every variant (P 0.8); the Nestle pots sit between DFA and Reading (P 0.6); O2 : reductant explains >= 70 % of the
variance (P 0.2); verdict OPEN / declared band (P 0.75), EXPLAINED (P 0.15), RESOLVED (P 0.1).

## 2. MFT in-pot loss: one constant or a pot-dependent law (EX-050; ingest finding 3)

Rivals: L0 pot-independent first order (the engine's k_mft_decay family, any Ea in [40, 133]); L1 cysteine-suppressed
(lambda0 / (1 + Cys/K_c), the EX-050 A3 form); L2 electrophile capture competed by cysteine (k_E E / (1 + K_cys
Cys), reversible per Esterbauer 1975); L3 O2-driven (k_ox O2); L4 apparatus loss (Grosch's open reflux condenser).
Test T-L: Grosch 1992 Table VIII lambda (re-derived: SDE / reflux = lambda T / (1 - e^(-lambda T)), +-1 SD) transported
to 110 C over Ea in {40, 60, 90, 101, 133}, against the patent bound lambda(110 C) <= 0.02 (slow exhaustion) / 0.05
(fast exhaustion) from the EX-050 red team (bands from the patent, not from Grosch). L0 survives only if a single
(lambda, Ea) satisfies both. Expectation: L0 survives only for Ea <= ~60 with fast exhaustion (P 0.7 that the
feasible L0 region is that small); K_c bounded to [~1, ~250] mM (P 0.6).

## 3. CF-015 (MFT/FFT 3.6 vs 0.50 between the two Nestle patents)

Rivals: R15a sugar : cysteine ratio (3:1 favours FFT via furfural); R15b standard timing (co-heated standards inflate
MFT/FFT by exp((lambda_M - lambda_F) t)); R15c temperature (95 vs 110 C); R15d concentration.
Tests: (i) within-study sign, Xu 2008 (ribose x3 at fixed Cys, 140 C, one lab, RF-free double ratio): R15a predicts
MFT/FFT falls with more sugar; (ii) cross-pot: MFT/FFT of every isotope-quantified pentose + Cys pot by ratio class and
IS timing (Hofmann 145 C, Schieberle 100 C, Yiltirak 100 - 130 C, both patents): R15b alone predicts the pre-IS pot
above ALL post-IS pots; R15a predicts the 1:1 pots above the 3:1 pots whatever the IS timing; (iii) the lambda_M -
lambda_F that R15b needs for the full 0.86 dex, against Kumazawa 2003's FFT survival (pH 5, 121 C) and the EX-050 bound.
Expectations: (i) sign agrees (P 0.8); (ii) 1:1 > 3:1 across labs (P 0.7); R15b needs lambda_M ~0.01 - 0.02 /min, i.e.
not excluded (P 0.6). Verdict EXPLAINED by R15a with R15b a co-term (P 0.5).

## 4. CF-073 (Weenen 2001 vs Hofmann 2000, Strecker phenylacetaldehyde from fed donors at pH ~3)

Rivals: R73a donor-generic operator (SDE stripping avoids a generic PA sink; FID / undecane quantification; a scale);
R73b donor-specific in-pot PA sink in the closed vial (PA + donor co-products) avoided by stripping; R73c
concentration-dependent share (a higher-order donor sink; 26 vs 100 mM); R73d donor speciation / handling (commercial
aqueous MGO / GO); R73e typo / basis.
Test T73: per-donor log ratio Weenen / Hofmann for MGO, GO, 3-DG (and glucose, two Weenen prints) at pH 3 - 3.2, with
noise 0.2 dex per Weenen value (its own 5x glucose irreproducibility is the upper end; 0.3 used as sensitivity) and
0.1 dex per Hofmann digitised value. R73a predicts equal ratios (chi-square on the spread, df = 2 for the three
dicarbonyls); a second-order competitor bounds R73c at <= log10(100/26) = 0.59 dex.
Expectations: R73a rejected at p < 0.01 (P 0.8); residual MGO gap after the maximal R73c term >= 0.6 dex (P 0.75).
Verdict OPEN, localised (P 0.75).

## 5. CF-031 (lipid lowers vs raises meaty thiols)

Rivals: R31a operator (headspace without isotope standards under-reads when an oil phase holds the thiol); R31b lipid
kind / oxidation load (electrophiles capture thiols; Doorn, Esterbauer); R31c structure (interfaces / mesophase raise
yields); R31d register reading error.
Tests: within-study signs per study (Whitfield 1988, Elmore 2002, Xu 2011, Yiltirak 2026 all arms, US 8,920,862,
du2025); operator size: partition factor 1 / (1 + phi K_ow) for MFT (log K_ow 1.5 - 2.0) at each study's lipid volume
fraction phi, against the observed fall.
Expectations: the operator explains <= 0.3 dex of the >= 1 dex falls (P 0.8); the register's "Yiltirak emulsion raises"
is relative to buffer + oil, not to buffer (P 0.9; checkable at once). Verdict EXPLAINED in part (P 0.5).

## 6. CF-016 (furfural 67.5 vs 5.2 ug, same Munich pot)

Test T16: operator correction of the two outside ratios the 10-29 analysis used (Meynier 1995 SDE / octadecane; Cerny
2007 SPME TIC): SPME fibre / water partition scaled as K_ow^0.8 (log K_ow furfural 0.41 measured, MFT 1.5 - 2.0
estimated). Expectation: the corrected Cerny ratio moves toward Hofmann's 3.4 (P 0.65). Verdict stays OPEN (P 0.85),
probabilities shift toward Schieberle under-recovery.

## 7. CF-014 (cause half)

Test: the covariate table (vessel, fill, IS, storage) for Zhou 2022 / 2023 / 2023c. Expectation: no recorded covariate
differs (P 0.85) -> OPEN, author request.

## 8. Pyrazine R (Shibamoto ~440 vs Leahy 7 - 70 vs Ghent <= 3.5)

Test: the 2,5- / 2,6-DMP co-elution bound for Ghent from the 2,6 / 2,5 ratios Leahy (resolved) prints per amine class
(lysine: no 2,6-DMP; asparagine: ~1.1x) and Shibamoto (NH3: 2.7x); the residual lysine-glucose gap Leahy / Ghent
after the bound. Expectation: co-elution explains <= 0.3 dex for Ghent's lysine pots (P 0.75); Shibamoto's excess is a
condition (8 M NH3, no amino acid), not an operator (P 0.8).

## 9. Data scored and exclusions

Scored: every isotope-quantified pentose + cysteine MFT level with stated conditions (Hofmann 1998 T1 / T2, Schieberle
2000 T4, Yiltirak 2026 S3 buffer, US 10,076,126 Fig. 1 - 2 buffer, US 8,920,862 Ex. 5 buffer, Grosch 1992 T5 no. 2 and
no. 1). Excluded from the level regression, with reasons: Bolton 1993 (thiamine supplies 92 % of MFT S, solids matrix,
dodecane added to the extract: no in-pot standard); Elmore 2002, Whitfield 1988, Xu 2011 / 2008, Meynier 1995, the
Jiangnan / BTBU pots (headspace or SPME areas without isotope standards: scale-free series, used only for within-study
signs); Grosch 1992 T5 no. 2 enters as a same-institute (DFA) point, not a second lab; Farmer 1989 / 1992 (no MFT
numbers).
