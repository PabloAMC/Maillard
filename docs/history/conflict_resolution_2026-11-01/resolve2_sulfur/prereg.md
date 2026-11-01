# resolve2_sulfur pre-registration (2026-11-01, evening)

Prior knowledge declared before writing: resolution_brief.md; ingest_findings.md Batch 2 (it already states the
dossier writers' readings: Karpinski phi ~1, Radomski pKa 5.2-6.1, van den Ouweland "4-mercapto-2-methylfuran",
Wei M 253, Zhang 2022 disulfide); dossiers karpinski1990, radomski1989, bustillotrueba2018; cluster5 prereg / verdict /
search_c5.json (E4c S1 best: k_r45(100) 4.3e-4 /min, Ea 41.5, s_lo 0.06, pKa_T 5.61 at the band floor);
split_v2_per_change.yml E6 / E7 / E10 blocks. Sec. 1 is written before any run; secs. 2-4 before their checks.

## Split status (decides what a result may be called)

- E10 (TTCA): DEV labs btbu_beijing (Xie: Cao, Zhang 2025b), ho_network (Jiangnan: Zhai 2019 / 2021, Kang, Wei),
  kallen_pennsylvania; TEST = F04#26 Mizutani 1960 only (UNTESTABLE on held-out, 1 lab). **Every CF-076 battery test
  (T1-T8) is a DEV-lab test.** Q1 is therefore a feasibility re-check of a provisional DEV admission, never a held-out
  score; nothing here can promote E4c.
- E6 (NF route): TEST labs btbu_beijing, caas_beijing (F23#5-8). E7: TEST labs ho_network, nestle/cerny, seville
  (F25#17, F30#11/15/16); reading_mottram_parker, caas_beijing, huazhong_pan are DEV(E7). Zhang 2021 (Hefei) and
  Zhang 2022 (Hefei): lab id to be checked; heldout_zhang2021_table3.yml is refused by fit_rows.load (TEST data).

## 1. CF-076 / E10: does the provisional E4c (E2c) admission survive phi ~1 and pKa_T 5.2-6.1?

**Rivals (neutral).**
- R-a: release is slow at 100 C because the C-N release barrier is LOW (Ea ~20-45 kJ/mol) and the pentose ring
  releases like TC at 25 C (phi ~1): the admission survives with a declared low Ea.
- R-b: release Ea is class-typical (50-80), phi ~1: the 25 C anchor then predicts k_r(100 C) 10-50x above what the
  DEV pots allow, so in-situ pH + Kallen shape cannot be the cause: CF-076 returns to OPEN.
- R-c: the 25 C anchor does not transfer (2-alkyl/aryl vs polyhydroxyalkyl acid limb; phosphate in Karpinski's
  pH 6 buffer; L-arabino vs D-xylo; ionic strength): phi(100 C) != phi(25 C). Not testable on disk; stated.

**Bands (none from T1-T8):**
- k25_6 = uncatalysed k_r(25 C, pH 6.0) of a pentose-Cys ring: [3e-5, 2.5e-4] /min, log-uniform. Upper: Karpinski
  AThz 1.8e-4 x 1.3 (RSD 30 %). Lower: AThz / Thz = 1.0 (0.55-1.9) x Kallen-corrected TC 5.4e-5 (= 3.0e-5), equal to
  1.8e-4 / 3.3 x 0.55 (Karpinski's TC is x3.3 Kallen's, read as phosphate / ionic-strength inflation of his pH 6 runs).
- Coupling: k_r45(100) = k25_6 x 10^(-1.1 s_lo) x exp(Ea/R (1/298.15 - 1/373.15)) (the E2c shape is flat from 5.6
  to 6.6, so pH 6.0 sits on the plateau). s_lo stays [0, 0.85] (Kallen shape). Bustillo's maximum at 4.4-5.2 would
  make s_lo < 0, which raises k_r45 and only makes feasibility harder: the Kallen shape is the lenient case.
- pKa_T: [5.2, 6.1] (Radomski charge balance) replaces [5.6, 6.6]. dpKa_CT unchanged.
- Ea: class band [40, 100] (primary); profiles at fixed Ea 30 / 50 / 64 (64 = cluster-5 transfer value).
- Everything else: cluster-5 bands unchanged.

**Tests and sets:** cluster-5 S1 (T1, T2, T4, T5, T6, T8a, T8b; T3 a shared-skeleton failure) governing; S3 core
(T4-T6, T8a, T8b) sensitivity. Search: differential evolution, maxiter 30, popsize 12, seed 3 (as cluster 5).
**Controls:** (i) positive: uncoupled E4c with the new pKa band must stay feasible (<= 2) on S1, else the wrapper is
broken; (ii) negative: NC1 / NC2 at mid-band must still fail S1 (> 3).

**Decision rule (cluster 5 / EX-054):** best max |z| <= 2 feasible; > 3 rejected; else inconclusive.
- coupled E2c or E4c feasible on S1 with Ea in [40, 100] -> admission SURVIVES (provisional, DEV only); report the Ea
  and k25_6 corner it needs;
- feasible only at Ea < 40 (profile) -> WEAKENED: "EXPLAINED (provisional) conditional on release Ea < 40", P down;
- infeasible (> 2) at every Ea profiled -> CF-076 returns to OPEN.

**Expectations:** positive control passes 0.9; NC fail 0.95; coupled family feasible in [40, 100]: 0.45; feasible
only below 40: 0.3; infeasible everywhere: 0.25. At Ea 64 fixed: infeasible (> 2) 0.7.

## 2. E6: van den Ouweland 1975 "4-mercapto-2-methylfuran" (written before re-reading the dossier and PDF)

Rivals: (a) misnaming of MFT (2-methyl-3-furanthiol = "3-mercapto-2-methylfuran"); (b) a real 4-thiol regio-isomer
(2-methyl-4-furanthiol); (c) transcription error in our dossier. Checks: drawn structure, MS (M 114 either way, so MS
cannot separate a vs b unless fragments / RI given), the paper's own numbering, what Whitfield 1999 and later work cite
it for. Rule: structure drawn with SH adjacent to CH3 (C-3) -> (a) or (c); SH at C-4 drawn -> (b); no structure and no
discriminating data -> UNDECIDED with P. Expectation: P(a or c) 0.6.
Second question: does "NF + H2S -> MFT" stand without van den Ouweland? List every fed-NF pot on disk with MFT
measured (Cerny pot D, Hofmann, Givaudan, Mottram / Whitfield) and its label evidence. Rule: >= 1 independent lab
with fed NF (or labelled NF) producing MFT -> the route does not depend on van den Ouweland. Expectation 0.85.

## 3. Wei 2022 mass (arithmetic, before re-reading the dossier text)

TTCA C8H15NO6S? Check: Cys C3H7NO2S (121.16) + xylose C5H10O5 (150.13) - H2O = C8H15NO6S 253.27. Deoxyosone
(3-deoxyxylosone C5H8O4, 132.11) + Cys - H2O = C8H13NO5S 235.26; + Cys without water loss = C8H15NO6S 253.27
(the carbinolamine / hemithioacetal). Prediction: the dossier's "2 H short" compares TTCA with a wrong M (237) and is
wrong; the corrected point is that m/z 254 cannot separate TTCA from ARP-reversion or deoxyosone capture by formula.

## 4. E7 / CF-031: which MFT species did each lab measure? (before re-reading Xu 2011 / Zhang 2021 methods)

Rivals: (a) a measurement-species artefact: labs that report free MFT only see oil "suppress" or "raise" depending on
whether disulfide formation is counted; total MFT (free + 2 x bis(2-methyl-3-furyl) disulfide) behaves the same in
all labs; (b) real capture by lipid oxidation products (aldehyde / hydroperoxide trapping of thiols), lab-independent,
dose-dependent in unsaturation / oxidation; (c) phase partition (oil takes MFT out of the water phase). Rule: if every
lab that measured the disulfide shows total-MFT suppression of comparable size (within x3) and labs that measured free
MFT only diverge, (a) explains the sign conflict and (b) the residual suppression; (c) is ruled out where suppression
exceeds the partition bound. Minimal E7 structure: only one supported by >= 2 labs on the same channel. Expectation:
P(species artefact explains the SIGN disagreement) 0.6; P(a cross-lab-supported E7 structure exists) 0.35.

**Amendment 1 (10:55, before any search result was read; the first run was killed at 6 min by my error, no output
read):** compute budget. maxiter 20 instead of 30 (popsize 12, seed 3, 2 workers). Job order: E4c coupled (S1), E2c
coupled (S1), E4c coupled Ea 64 (S1), positive control E4c UNcoupled with the new pKa band (S1, same maxiter, so it
also measures whether maxiter 20 reaches the cluster-5 optimum ~1.45), E4c coupled Ea 30 (S1). Further profiles only if
time allows. Negative controls already run: NC1 min max|z| 58.2, NC2 33.7 over 20 draws (both fail S1 as designed).

**Amendment 2 (POST HOC, after reading the first two results: both best points sit at the k25_6 floor, 3.0-3.2e-5,
with Ea 40.6-42.9 near the class floor).** Labelled diagnostic, not a verdict input: profile k25_6 FIXED at the
face-value Karpinski AThz 1.8e-4 /min (no phosphate deflation) and at the geometric mid-band 8.7e-5, E2c only (the better-scoring member; compute budget), Ea
free in [40, 100], S1. Expectation: face value infeasible (> 2) 0.75.
Amendment 2b (POST HOC, after the Ea 30 profile came back at 1.10 with k25_6 again at its floor): in the face-value
profiles Ea is widened to [10, 100] (below the class band on purpose), so the question they answer is "what Ea would
the printed Karpinski rate need", not "is it feasible in the class band".
