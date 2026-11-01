# Cluster 3 pre-registration: in-situ pH and speciation (2026-11-01, written before any test script ran)

Honesty note: while reading the sources I did back-of-envelope speciation arithmetic for CF-065 and CF-054 in my head
(no script, no engine run). The expectations below are stated after that arithmetic and are therefore not naive priors;
the decision rules were fixed before any script was written. Existing sub-model 1 artifacts (fit arms A/B/C, hold-out
score) were READ before this prereg; they are cited as evidence already on disk, not as a new test.

Shared bands (declared from data OTHER than the scored quantities):
- H2S pK1(T): Suleimenov & Seward 1997 correlation as shipped (`parameters_sulfur.h2s_pk1_suleimenov_seward`);
  Davies shift via `pka_at(..., ionic_strength)`.
- H2S Kaw(T): `data/lit/henry_constants.yml` (Sander: H_cp 1.0e-3 mol/(m3 Pa), B = 2100 K). Band: van 't Hoff to the
  reaction T (upper) vs clamped at 100 C as `src/headspace.py` does (lower).
- NH4+ pKa 9.245, dH 51.95 kJ/mol (CRC / GKL textbook), dCp 0 (band +-40 J/(K mol)); NH3 headspace ignored (Kaw ~1e-3).
- Ea of a sulfide + carbonyl step: 40-100 kJ/mol (class band; no measured value on disk), central 70.
- Vessel: Hofmann 50 mL in a 200 mL autoclave (Vg/Vl = 3, printed). Zhai 2023 JAFC: not stated; band Vg/Vl 0.5-3
  (same lab's stated vessels: Wang 2026b 10 mL in 15 mL = 0.5; Hou 2026 10 mL in "high-pressure bottles").
- Effective time at temperature: Hofmann 10-20 of 20 min (autoclave come-up not stated); Zhai 6-10 of the first 10 min.
- Replicate CV of total organic acid in one recipe: ~10 % (Davidek 2002 Table 5 printed SD / mean 0.10, 0.10).

## CF-065 (FFT from furfural + sulfide; MFT from norfuraneol + sulfide; Hofmann vs Zhai 2023 JAFC, 50-180x)

Rivals:
- R1 COND, speciation + headspace + temperature: the nucleophile is dissolved HS-; Hofmann's pH-5, 0.5 M phosphate,
  Vg/Vl = 3, 145 C pot has far less dissolved HS- per mol of sulfide charged than Zhai's (NH4)2S pot at pH 7, 120 C.
  Predicts: implied second-order k on dissolved HS- agrees within 0.5 dex.
- R2 CHEM, NH3 co-reactant: (NH4)2S brings 20 mM ammonia; imine / iminium activation or a separate N route consumes
  furfural and norfuraneol fast. Predicts: a residual > 0.5 dex after R1, and furanoid loss much faster than the
  sulfide product formation accounts for.
- R3 OP, Zhai semi-quantitative thiol scale (RF = 1): predicts a residual of similar size for FFT and MFT (a common
  factor), not removable by any condition correction.
- R4 TU: basis slip (mol % of precursor vs of H2S; both 1 mmol so no effect expected).
Test: compute k_HS (and k on total dissolved sulfide) for both labs with the bands above; Monte Carlo / corner scan.
Decision: corrected median gap <= 0.5 dex and band includes 0 -> RESOLVED (COND); 0.5-1.0 dex -> EXPLAINED in part,
residual to R2 / R3; > 1.0 dex after correction -> speciation is a minor part; residual is declared.
Negative control: Hofmann Table 7 (same lab, pH 5, 1 vs 2 mmol H2S): the operator must predict ~2x (measured 2.18x);
and the operator applied to a matched-condition pair (Hofmann vs itself at swapped vessel fill) must leave k unchanged.
Expectations: P(gap <= 0.5 dex after R1) 0.2; P(R1 removes >= 1 dex but leaves > 0.5) 0.55; P(R1 removes < 1 dex) 0.25.
Note: the "Hofmann Table 8" check in the triage does not apply: Table 8 feeds mercapto-2-propanone, not H2S.

## CF-054 (Chen 2000 vs Chen & Ho 2002 final pH 6.05 vs 4.96 after an alkaline start; 3.86 vs 4.01 after pH 5)

Rivals: E1 COND run-to-run acid-yield spread amplified by the alkaline start; E2 TU / OP reading or adjustment.
Test: charge balance of the 180 C pot (ribose 100 + Cys 100 mM, NaOH to the label pH): the acid needed to reach each
endpoint, from both starting pH; ratio required vs the ~10 % replicate CV.
Decision: if the 1.1-unit gap needs <= +-20 % (2 sigma) acid-yield spread from pH 8.5 AND the pH-5 pair's 0.15 gap
needs a comparable or smaller spread -> EXPLAINED (E1, P >= 0.6). If the pH 8.5 gap needs > 50 % -> E2 favoured.
Expectation: P(E1 passes) 0.65.

## CF-069 (acetic acid loss r_aa_frag 0.0164 /min vs Davidek 2008 plateau)

Rivals: M1 MOD, r_aa_frag absorbs the pH-driven slowing of acid production in a fit that holds pH or under-couples it;
M2 OP, late analytical loss of AA (volatility, extraction); M3 CHEM, a real glucose-system AA consumer.
pH angle tested here: (a) AA volatility at in-situ pH from Kaw (AA H_cp ~ 4e1 mol/(m3 Pa)) -- predicts negligible loss
in sealed tubes; (b) existing sub-model 1 arms A (k_aa_frag free) vs B (off), both with the inline pH coupling: if B
fits Martins 2005 AA and the held-out Martins 2005b AA as well as A, the loss is not needed once pH is coupled (M1).
Decision: B within tolerance on both and within 0.05 dex RMS of A -> EXPLAINED (M1); otherwise OPEN.
Expectation: P(M1) 0.55. (Arm E on Davidek 2008 remains the register's test; not run here.)

## CF-012 residual (Zhou endpoints; the pH-3 norfuraneol misfit)

Test: (i) out-of-sample positive control of the pentose drift model on Hou 2026's MEASURED trajectory (GX-ARP 5 mM,
water, 7.0 -> 5.05 at 20 min -> 3.55 at 120 min), shipped drift constants and sub-model 1's arm-B values, no refit;
(ii) in-situ pH of Hofmann's buffered pH 3 / 5 / 7 pots at 145 C (does the pH-3 pot run nearer 3.5?).
Decision for the systematic fix: the pH model must reproduce Hou within 0.3 pH RMS (positive control) and keep a
0.2 M pyrophosphate / 0.5 M phosphate pot within 0.2 units (negative control, Meynier's measured <= 0.2 drift).
Expectation: P(Hou within 0.3 RMS) 0.25 (sub-model 1 already fails Cao, Zhao, Zhou by 0.7-0.8).
