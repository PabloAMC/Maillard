# Pyrazine rivals: red-team notes (2026-10-31)

Scope: prereg.yaml, pyrazine_rivals.py, pyrazine_rivals_result.json, posthoc_summary.json, the dossiers, the engine.
View formed before reading README.md; comparison at the end. Nothing in src/, data/ or the pyrazine/ folder was edited.
Scratch scripts: pyrazine_redteam_work/label_check.py, pyrazine_redteam_work/p2b_ridge.py (run through docker_maillard.sh;
a first pass of label_check.py was run with system python3, contrary to the repo rule, then re-run in the env).

## 0. Timestamps
`ls -lT`: prereg.yaml 13:51:49 < pyrazine_rivals.py 13:54:26 < result.json / run1.log 13:54:38 < posthoc 13:57. The
prereg's text says "written 13:48-14:20" and "amendment_1 14:30, BEFORE the script existed", but its mtime is 13:51:49
and the script's 13:54:26. The order on disk is right, but the stated clock times are wrong. Only the last write is
visible, so earlier runs cannot be excluded.

## 1. Transcription (checked cell by cell)
- Leahy Table I (leahy1989_extraction.md:99-156): all 42 (k, r2, n) triples used match. Table III (:208-211) matches.
- Ho rows: zhou2023 T1 (:157-160), zhou2022 T3 (:111-114), zhou2023 Fig 3 (:284-298, :369-383), zhou2023c T1 (:119-122) match.
- zhou2024 Table 1 (:97-108) matches. Adams Table 3 (:120-127) matches. Amrani Table 2 (:176-177) is 0 %.
- Low 2007 (:101-102): methylpyrazine is 4.8 % labelled at 15 min and 0 % at 30 min; the script uses 0.07 (conservative).
  2,5-DMP is 6.7 % (SL 6.3 + DL 0.4); the script uses 0.063. Both are trivial.
- Not a transcription error but an omission: zhou2023c ARP + Gly (PZ 5.97, MPZ ND, DMP25 1.16 ug/L) and ARP + Cys
  (DMP ND) are listed in the prereg S2 data but never scored (section 5 below).

## 2. P0 is NOT the engine. The engine is tighter: R = 1, not <= 4
network.py:361-372 sets three condensations on ONE k_cond with stoichiometry {AKG: 2}, {AKM: 2}, {AKG: 1, AKM: 1}.
`reaction_rates` evaluates rate = k * prod(c^nu). So PZ = k A^2, DMP = k M^2, MPZ = k A M, and
R = MPZ^2/(PZ DMP) = 1. Verified in the env: AKG 0.3, AKM 0.7, k = 1 gives R = 1.0.
The docstrings (parameters_pyrazine.py:27-28, network.py:372, b18 fit reports :53) claim "rate_MPZ = 2 sqrt(rate_PZ rate_DMP)",
which would give R = 4. **That is an engine bug: the cross pair lacks its statistical factor 2.** The prereg copied the
docstring, not the code. Consequences:
- Every calibrated Ho row with R >= 1.86 rejects the ENGINE (e.g. zhou2023 ARP+Cys pH 8: R 2.43, z ~ 3.0; Cys+MGO pH 7
  90 min: R 3.73, z ~ 3.8). So P0's S1/S2 rejection is two-lab, stronger than claimed.
- Methylpyrazine is a pure prediction (B18 fits only the PZ and DMP rates). Against observed R it is low by sqrt(R_obs):
  x2 against the engine's own design, x3.5 for Leahy Lys-Glc, x7 for Lys-Rib.

## 3. The mixing bound: validity and artefacts
The bound R <= 4 needs:
1. each pyrazine from one pairing only (extra PZ or DMP routes lower R, so only extra MPZ routes break it);
2. k_AM^2 = 4 k_AA k_MM (random pairing);
3. EQUAL DHP -> pyrazine yields (no pair-specific diversion);
4. no differential post-formation loss;
5. the same weight measure for all three products (same window and volume), so that Cauchy-Schwarz gives
   (int A M)^2 <= int A^2 int M^2.
Under these, heterogeneity only lowers R, and the bound is valid.
- Leahy's slopes come from different windows: DMP regressions drop points below detection (n 14-16 against 20-22),
  and the intercepts are positive (dossier flag 6), so condition 5 is only approximate for the rate rows. The 2 h
  cumulative end points (Table III), the clean C-S statistic, give the same picture (R 7.5-52).
- RF error cannot do it. Leahy's RFs were purged from borate (in matrix), and NPD response per mole is roughly N-count
  based. A shared RF bias is ONE number for every Leahy pot (same compounds, same method), yet R spans 5-70 across
  pots, and Lys-Rib needs > 1 dex.
- Artefacts the power check does not contain:
  - **NPD sees every N compound.** MPZ identity was MS-confirmed on one Lys-Glc 6 h sample only, so Lys-Rib (the
    strongest row) has no MS check against co-eluters on DB-225.
  - Capped tubes for up to 24 h at 95 C: alkylpyrazines are more volatile, so DMP loss raises R.
  - Rough probability that the bound violation is artefactual: ~0.15 for Lys-Rib, ~0.3 for the 2-4x rows.
- **A missing rival, P1y: pair-specific DHP fate (condition 3).** The AA-DHP is unstable, and the MM-DHP is diverted
  to aldehyde adducts (Adams Scheme 1A). On R and labels, P1y is observationally identical to P1x. Unlike P1x it is
  pot-dependent, so the post hoc cross-sugar argument against P1x does not touch it. Leahy saw no TrMP or
  ethyl-DMP, which excludes volatile diversion but not polymeric diversion.

## 4. Adams pyrrole mass balance (S5a)
- Both products use the same IS with "response factors calculated". The pyrrole reference is 78 % pure (<= 0.1 dex if
  uncorrected).
- DMP (MW 108, volatile, water-soluble) went through CHCl3 extraction and a 10x evaporative concentration. Its
  recovery can be well below the pyrrole's (MW 165) unless the RFs were run through the work-up, which is not stated.
  Even a 10x DMP loss leaves B > 0 in the pH 4 rows, but the 1:1 rows then drop below 2 sigma.
- **The interpretive assumption is the weak point.** Pyrrole 34 is C9 = 3 x C3. Adams' "aminoacetone + 2 MGO" is a
  proposed scheme with no label, and NH3 or amine + 3 MGO gives the same molecule. If the pyrrole does not pass
  through FREE aminoacetone, S5a says nothing about P0's flux booking.
- S5b (the ratio rises with MGO) is predicted by both readings.
- Conditions: ~850 mM, pH 4 acetate.

## 5. Circularity and unscored data
Band and key issues:
- op.mh_leak_mpz top 0.30: D_obs ~ 0.13 was hand-computed BEFORE the band was set (prereg disclosure), and
  D_obs ~ 0.13 needs r ~ 0.26 to explain it. The band top decides S3's primary outcome. The basis ("methylpyridines
  ~0.3") is not from the scored data, but this is the riskiest band.
- op.near_loq_dex was set from the very rows it demotes (zhou2022 / zhou2023 ARP pH 8, R 7.7 and 14.2). Those rows
  are then dropped as unpowered. The direction is conservative.
- pyr.fM_sd [0, 0.2] is in the code but not in the prereg. Its top is about the DMP overdispersion of zhou2024
  (s ~ 0.19), but its effect on D is ~0.02.
- S6's GO lifetime offsets (-0.5, +0.4) are undeclared.
- S6 and S7 bands come from other tables of the same paper or table. These are minor.

Unscored data:
- **zhou2023c ARP + Gly**: MPZ ND while PZ and DMP are present. R < 0.19 / 0.11 / 0.035 for an MPZ limit of
  1.15 / 0.88 / 0.5 ug/L. Every pairing rival has R >= 4h >= 0.4 (h_any >= 0.1). This row is FORBIDDEN by P1, P1x
  and P2, unless h < 0.05 or the Ho low-level data are unreliable. It is in the prereg data list but was neither
  scored nor excused. README: "Every Cys row has R < 4, which no rival forbids". Rows near zero ARE forbidden.
- Ho non-Cys calibrated rows: all 4 (ARP pH 8 x2, GSH, Glu) sit above R_max (sign p ~ 1/16). Every Cys row sits
  below. This is direction-only second-lab support, unreported.

## 6. Lab independence
- Reineccius: one lab (Leahy).
- Ho network: one lab. Zhou 2022/2023/2023c/2024 share the same ARP and SPME method; Hwang 1995 is Rutgers/Ho;
  Deng is Ho network.
- Ghent: one lab (Adams; Van Lancker is the same lab).
- Nestle and Reading are independent of each other.
- DFA: one lab.
- P0's rejection: Ghent (S5a, S7) + Reineccius (S1) + Ho (S2, given R = 1). Robust.
- P1's rejection: Reineccius only on prereg'd tests.

## 7. Power check
- The synthetic noise equals the scoring noise, so power and size are tautological against the declared noise. The
  rival family is only the nested P1x / P2, with no systematic-artefact rivals (co-elution, loss, window mismatch)
  and no P1y.
- The S3 null lets f_M float over [0.1, 0.7], although DMP pins it at 0.43-0.46. This throws away power.
- Size at the cap (0.10, two-sided).
- S6 at 0.38 was correctly dropped.

## 8. P1x
- It is pre-registered, not post hoc. Physically it is the EXPECTED default: imine formation with aldehydes is
  10-1e3x faster than with methyl ketones. So R = 4 random pairing is closer to a strawman than P1x is.
- It explains R > 4 in one pot. It does not explain the 10x spread across sugars at fixed T with h >= 0.5. The
  posthoc argument is post hoc, and the h_linear basis is weakened by Leahy flag 6.
- **A test the analysis did not run (post hoc, mine): cross-compound 15N mean label.** Pair constants do not move
  labels, so under P1, P1x and P1y: label(MPZ) = (label(PZ) + label(DMP)) / 2.
  - Zhou 2024 at 30 min: 0.371 observed against 0.489 predicted, a deficit of 11.7 pt (~9 sigma at 1 pt per fraction).
  - At 120 min the deficit is 8.2 pt.
  - Inverting [M-H]+ leakage for MPZ and DMP alike: the deficit stays 9-10 pt at r = 0.3 and 7 pt at r = 0.5.
    Leakage cannot rescue it, because DMP leaks at least as much as MPZ.
  - The only escape is M-label heterogeneity CORRELATED with A flux, which DMP's overdispersion (s ~ 0.19) makes
    conceivable.
  - So README's "the [M-H]+ question alone decides P1x vs P2" is too strong. The label data already disfavour P1x/P1y.
- P2b with a low-label (ARP-N) extra route fits the MPZ isotopomers for w >= 0.49, i.e. R/h >= 15.
  zhou2022 / 2023 ARP pH 8 give R 7.7 / 14.2 (near LOQ). Consistent, not tight.

## 9. Code
- lab_level: "weakest_rejecting_row" returns the MAX-z row (mislabel). The prereg's per-system min-z rule is only in
  posthoc_summary.py.
- S4: P2a_aa "passes_obs": False is hard-coded. At phi = 0.5 the allowed w is 0.30, i.e. exactly the material
  threshold (a boundary, not a rejection). P2a_aa (Gly/Ala HCHO, CH3CHO) is irrelevant to Leahy's Lys/Asn pots
  anyway.
- Fig 3 rows add half the calibration term inside sig_pot and the full shared term again (double count, conservative).
- `rd` is unused.
- No bug changes a verdict except the P0 reading (section 2).

## 10. Comparison with README
- P0 REJECTED: agree. The reason should be corrected (R = 1 engine bug), which makes it two-lab on S1/S2.
- P1 INCONCLUSIVE by rule: agree by rule; P(wrong) ~0.9 with the post hoc label test.
- P1x "ADMITTED, PROVISIONAL": disagree. Inconclusive, leaning rejected as the sole explanation (cross-sugar spread
  plus the label deficit).
- P2b admitted provisional: agree, but it is admitted largely through flexibility (free e per pot).
- P2a_aa rejected in 180 C matrices: agree, but it is irrelevant to the S1 violation.
- Missed: P1y, the Gly row, the label identity, the engine factor-2 bug.
