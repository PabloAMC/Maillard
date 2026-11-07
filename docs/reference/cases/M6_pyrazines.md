# Case M6: pyrazines

**STATUS (2026-11-07): NOT FROZEN. Items 1 - 4 complete; item 5 sized but with one pair undecidable by any computable
design (the claim vs the incumbent on part C; "Freeze record" at the end says exactly what blocks the freeze).
Development grade on disk; decider W6 (decisive-measurements entry `pyrazine_ring_aminoketone`, revised, and a new
entry `pyrazine_amine_budget`, proposed in the PR body).**

- **Falsified on disk** (spent or seen data, constant-free structures, falsify-only; sec. 4.2):
  - **FX-C1:** a free-aminoketone pool with a finite condensation at Elfsson 1999's class and **no fast first-order
    aminoketone exit** (MF-35 candidate (a), CK-52's form) is falsified by Zhou 2024's six zero-order R2 (Holm p <=
    0.026 at the least favourable noise; the analogue would need a constant 3.2 - 5.1 dex faster). The claim R1
    survives only with an exit of **ks >= 0.028 /min** at 20 mM MGO (its pyrrole, or any first-order aminoketone loss).
  - **FX-S2:** "the amine-dependent MGO loss is a 1:1 amine-consuming adduct" is falsified within Kanzler 2026 (>= 61 %
    of the Ala must go; 100 % is left; p 0.0024).
- **Not falsified:** the incumbent's "Strecker is the only amine exit" (FX-S1: Holm p 0.09 / 0.11 under the largest
  listed noise); R3; R1 with a fast exit; the adduct sink; every pairing rival EX-052 left; the supply rival.
- **Credits: none.** Every unit on disk is spent, seen, cross-design or a TEST cell (sec. 2.10).
- **The minimum deciding design** (sec. 4.4, power in sec. 5.4), within-run units only, three fed Ala + MGO pots:
  - **arm M** (one pot, MGO raised x3 at 40 min; n = 2 analyses per aliquot) decides **C3 against C0 and C1**
    (P 0.94 - 0.98 each way; false credit <= 0.06; the M1.1 gate passes);
  - **arm S** (two pots, no buffer and 0.4 M phosphate; n = 3) decides **S0 / S1 / S5** pairwise (P 0.86 - 0.98);
    both pots are necessary;
  - **C1 (the claim) vs C0 (the incumbent) is NOT decided** by any computable design: arm M gives 0.62 - 0.70 at n = 3,
    and the arm that separates them structurally (the aminoketone trap K) leaves the null 23 % unresolved, so under
    the adopted bound C1 can never be credited. This blocks the freeze.

**Item 5 follows sec. 7 item 5 as ADOPTED** (PR #102 / #105; `src/kinetic_core/predictive_densities.py`). The measurement
model of sec. 5.1 is declared in that module's terms (`UnitNoiseModel`, `NoiseTerm`), and the power check (sec. 5.3 - 5.4)
runs the adopted judge (`build_null_pool`, `judge_from_draws`, the misspecified false-credit gate). Items 6 - 9 are
appended later and never edited into this part. **No engine was run for this case, on any pot**; no engine constant
was fitted. The changelog is at the end.

*Plan step M6 ([STATUS_AND_PLAN sec. 5, W2](../../guides/STATUS_AND_PLAN.md); errors E13 (F19-a, CF-077) and E14
(F19-b, F21-a, F20-a) of [ENGINE_VS_LITERATURE](../pathway_atlas/ENGINE_VS_LITERATURE.md) sec. 3, amended in secs. 5.1
and 7 (E13 / E14 and the pyrazine rivals), with CF-095 (the dicarbonyl-amine adduct sink) and E4 (decided: a declared
band, owner 2026-11-01) as neighbours; M0's pyrazine lane MF-34, MF-35, MF-36 and groups G6, G8, G2
([misfit inventory](../misfit_inventory.md)). Items 1 - 5 of sec. 7, written 2026-11-04 and revised 2026-11-07, before
any result of the change exists. M6 shares the Strecker steps `r_go_ak` / `r_mgo_ak` with case
[M3](M3_strecker_aldehydes.md) phase 2.*

**Rules this case applies (sec. 7 as adopted 2026-11-06 / 07).**
- **Spent data** (sec. 7 item 2, M4's FM2): data an M0 check or an earlier diagnostic used to nominate a rival are
  development evidence; they may FALSIFY a rival whose prediction there is constant-free and never CREDIT one. The rule
  applies to R0 exactly as to the new rivals (S1 - S4).
- **A falsification** is a posterior-predictive p after Holm across the datasets the rival was tested on (S3, M1.1).
  On values with no run-to-run noise it stands only under the largest noise any listed source supports. The pyrazine
  lane has **no H0** and no pooled reaction-class between-laboratory spread (`src/evidence_tiers.py` carries none),
  so: within-study values use their printed SD (n - 1 df, an unstated n taken as 2) with reading errors as floors;
  cross-laboratory arithmetic uses H0's omega (0.419 dex per value) and tau_ABS (0.964 dex per laboratory), the only
  measured terms in the programme, as "the largest listed", and each falsification prints its break-even SD.
- **No cross-laboratory comparison can credit** (no measured between-laboratory spread of any pyrazine-lane
  observable); no noise term is transferred across observables; no bound is used as an SD.
- **B18's six constants are fitted-unidentified** (CF-136: B18 does not reproduce on today's network; S 2.6 -> 295) and
  `k_cond` is structural (`incumbent_evidence_tiers.json`). R0 is therefore tested as a structure with its constants
  in their S4 bands, never at B18's point alone.
- **pH scale.** Unless stated, a pH is the label pH (prepared, room temperature). Zhou 2024 and Kanzler 2026 are
  unbuffered (initial pH only); Adams, Piloty and Hofmann are buffered (label).

**Seen by the author.**
- 2026-11-04 text: engine outputs on development or calibration rows only (the B18 outcome; EX-057's moves of Leahy's
  pH-ratio rows and PH-06 / TEMP-04; CF-095's probe of Piloty's pot; EX-052's scores S1 - S7). Measured values only
  for the v2 E13 / E14 TEST cells (F19#0, #1, #13, #17; F20#1, #2, #3, #13; F21#8, #10), for the keyed held-out Lee
  2024 and Sedehi 2013 tables as quoted in CF-097 / CF-095, and for the UC Davis papers quoted in CF-087 / CF-098.
  No engine output on any of them. The four E13 TEST studies were neither scored nor excused by EX-052 or EX-058; that
  debt stands.
- 2026-11-07 revision: the dossiers' printed values of Zhou 2024 (Tables 1, 2, SI), Adams 2008 (Tables 1 - 3), Piloty
  1979 (Tables 1 - 7 as transcribed), Hofmann 2000 (Tables 3, 4), Elfsson 1999, Kanzler 2026 (Figs. 1D - F, 2B);
  CF-077, CF-087, CF-095, CF-135, CF-136; EX-021, EX-052, EX-057, EX-058; the misfit inventory; the S1 record's
  pyrazine rows. The citation sweep (sec. 2.9) read titles and abstracts of unseen works only. **[dry]** marks
  arithmetic on printed values (job (a) of `src/m6_reduced_power.py`); it is not a result of the change.

**Split notation.** "v1" = [split.yml](../contrast_benchmark/split.yml); "v2" =
[split_v2_per_change.yml](../contrast_benchmark/split_v2_per_change.yml).
- v2 E13: DEV labs Ho network, Minnesota (Reineccius). TEST: F19#0 / F20#1@amranihemaimi1995 (Nestle), F19#1 /
  F20#2 / F20#3@low2007 (Reading), F19#13@vanlancker2012 (Ghent), F19#17 / F20#13@shibamoto1977 (UC Davis).
- v2 E14: DEV labs DFA, Ghent, Ho network, Karolinska, Nestle. TEST: F19#1, F19#17, F20#2, F20#3, F20#13,
  F21#8@zamora2020, F21#10@elmore1997. Excluded as calibration input: the Leahy cells F19#8 - #11, F20#12.
- So Ghent is DEV for E14 and a TEST laboratory for E13; anything this case reads from Adams 2008 about the mixing
  ratio R is E13's held-out evidence and stays in item 7.
- v1 TEST cells named below: F19#17, F20#13 (UC Davis), F19#19 (Yu 2018, NUS), F21#11 - #13 (De Haan / Sedehi),
  F21#5 (Zhou 2000, Kentucky). Keyed held-out, outside the frozen splits: `heldout_lee2024_fig2.yml`,
  `heldout_sedehi2013_tab1.yml`.
- Piloty 1979, Koehler 1969 / 1970, Jousse 2002, Kanzler 2026 and Bemis-Young 1993 have no atlas cell. Piloty,
  Koehler, Jousse and Kanzler sit in no frozen split; Bemis-Young is UC Davis, an E13 TEST laboratory.

## 1. The claim

Two changes downstream of the Strecker step, plus the condition the plan row does not name.

| part | step(s) | today (B60 + CF-077, read from the code) | M6 |
|---|---|---|---|
| S. dicarbonyl-amine adduct sink | new: MGO (and GO) + amino acid -> non-Strecker adduct, irreversible, pH-flat | the only MGO + amine step is the Strecker `r_mgo_ak` (k ~3e-8 L/(mmol min) at 100 C); `r_mgo_mel` consumes no amine ([network.py](../../../src/kinetic_core/network.py) ~180, ~359 - 365) | on, its consumption law from Piloty 1979 (CF-095's candidate): a rate constant >= 5.6e-3 M^-1 min^-1 at 100 C [dry, a lower bound], first order in each |
| C. finite aminoketone fate | `r_akg_pz`, `r_akm_dmp`, `r_ak_mpz` on `k_cond`; new: aminoacetone + MGO -> pyrrole | `k_cond` = 1e4 L/(mmol min), Ea 0, no pH ("declared fast" from Jousse 2002; [parameters_pyrazine.py](../../../src/kinetic_core/parameters_pyrazine.py) ~53, ~122 - 124); cross pair x2 since CF-077 (PR #37); no pyrrole, no dihydropyrazine | `k_cond` from a measured class (Elfsson's pH law and barrier); **the pyrrole as a first-order aminoketone exit that job (a) requires to be fast: ks >= 0.028 /min at 20 mM MGO** (sec. 4.2) |
| not claimed | the pyrazine supply (the trunk's GO / MGO in water), the amine identity, the mixing ratio R | B18: glyoxal reaches a sugar pot only through B13's dry-glass glucosone; glycine stands for every amine; R = 4 (random pairing) since CF-077 | named as rivals (part L, part P) because the level gap and the pairing sit there (sec. 3) |

**What the 2026-11-07 revision changes in the claim (before any result of the change).** The 2026-11-04 text read D1
as "the measured analogue rejects the finite condensation as built, unless the ring bypasses free aminoketone". That
holds only for a pool with no fast exit. The claim's own pyrrole is a first-order exit at fixed MGO, and with it the
pool relaxes at ks, not at the condensation's tau. Job (a) puts the edge at ks >= 0.028 /min, and at that edge the
claim needs a Strecker constant within 0.22 dex of the rate Hofmann 2000 implies (EX-021), a DMP yield of ~0.005 per
Strecker event, and ~170 uM of exit product by 120 min in Zhou's DMP pot (sec. 4.2). This is stated as the claim's
signature, not as evidence for it (Zhou and Hofmann are spent).

**Earlier waves on these structures** (cited, none re-run):

| structure | wave | outcome | kind |
|---|---|---|---|
| Strecker x2 fitted, condensation declared fast | B18 (SHIPS, two caveats) | Zhou 2024's six rates within 0.07 dex, Leahy's four pH ratios within 0.13; **T3 failed** (Leahy pyrazine : 2,5-DMP 23 : 1, model 6e-8 : 1), **T4 failed** (total 2.9 dex low, glycine for lysine); apparent barriers 300 - 460 kJ/mol vs measured 150 - 177: "the trunk's dicarbonyl SUPPLY" | fit, on fed-dicarbonyl rows |
| B18 at its frozen point today | PR #107, CF-136 | **DRIFTED-MATERIAL**: S 2.64 -> 295 (B21's glucosone step, CF-077's factor); `k_go_ak` / `k_mgo_ak` log k and both pH slopes outside today's local 95 % region; the six coordinates fitted-unidentified (S4 bands) | reproducibility audit |
| the glucosone route to glyoxal | B21 | Leahy total unchanged (-2.8 dex); the trunk's MGO 7.8 mM where Xia's panel tops at 0.5 mM | fit |
| P0 (engine) vs P1, P1x, P2a / b, NC | EX-052 | P0 rejected (Ghent S5a, S7; Reineccius S1; Ho S2); P1 inconclusive; P1x provisional; P2b provisional; P2a_aa rejected by 13C | constant-free signatures |
| the cross-pair factor | EX-057 / CF-077 | fixed (R = 4); the frozen B18 pH slopes are conditional on it (Leahy rows +0.51 dex at the fit commit); a B18 refit awaits the owner | diagnosis |
| Piloty, R pool invariance, form of E5 | EX-058 | Piloty not discriminating for E4 (P 0.85); exposes an adduct sink >= ~1e3x the engine's Strecker consumption (a lower bound); P2b back to ~0.45 | arithmetic + probe |
| the Strecker step's pH law (E4) | E4 prereg | expected REFUSE; declared band | prereg |
| B18's Strecker constant | EX-021 | a pyrazine-yield constant ~2.3 dex below Hofmann's fed MGO + Phe rate | arithmetic |

**Errors stack (CF-107).** B18's T3 / T4 failures are whole-model outcomes on the trunk. Candidate compensating errors:
- the trunk's GO supply (dry-glass glucosone only) and B13's glyoxal sink (98 % of a fed 20 mM GO in 2 h);
- the trunk's MGO level (B21) and E2; MGO's own loss law (CF-122, CF-135: B1's `k_mgo_mel` drifted with B21);
- glycine for lysine (Leahy's amine; Lys >> Asn x35, F19#10);
- the Leahy pH law on the Strecker step (E4), conditional on CF-077 and drifted (CF-136);
- k_go_ak's 70x conflict between Zhou 2024 and 2025b (B18 post-run note);
- instant condensation and the missing adduct sink (this case).

## 2. All the evidence, by kind

### 2.1 Isotope labelling

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Zhou 2024 (F19#3, F20#10, F03#15) | Ala-xylose ARP + 15N-Ala, 120 C: pyrazine 15N share flat (54.5 -> 54.3 %) over 30 - 120 min; free Ala supplies ~55 % of ring N; MPZ 37 vs 49 % (the label-average deficit, EX-052 red team) | aminoketones from free amino acid on ARP-derived dicarbonyls; pool mixing; pairing (part P) | DEV(E13, E14); calibration input; **spent** (EX-052 S3) |
| Zhou 2025b (F19#2, F20#11) | pyrazine isotopomers close to one well-mixed C2 pool; MPZ M+3 absent at 30 / 60 min, ~4.5 % at 90 / 120 | a slow minor MPZ channel | DEV; seen |
| Hwang 1995b (F21#3) | > 90 % of 1-methylpyrrole N is glycine's | pyrrole N from the amino acid (a different pyrrole) | DEV (Ho network) |
| Amrani-Hemaimi 1995, Low 2007 (F19#0, #1; F20#1 - #3) | MPZ / DMP ring carbons unlabelled from labelled amino acids; DMP 25 - 30 % from labelled Gly | aminoketone pairing; aldehyde addition | **v2 TEST(E13 / E14)**: item 7 only |
| Fang 2009, Mundt 2004 (F31#0 - #3) | 27 - 33 % vs ~2/3 of glycine decarboxylated in melanoidins | total Strecker-type flux into polymer | Fang **v1 TEST**: listed, not used; Mundt DEV |

### 2.2 Fed intermediates

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Zhou 2024 (F19#4; zhou2024 dossier) | Ala + GO or + MGO 20 + 20 mM, pH 8 (initial, unbuffered), 100 / 110 / 120 C: pyrazine 0.0279 / 0.0791 / 0.1507, 2,5-DMP 0.0035 / 0.0100 / 0.0230 umol/(L min); **zero-order c = kt + b over 0 - 120 min, R2 0.963 - 0.995** (five time points, triplicates, no SD printed); no MGO, GO, Ala or pyrrole measured in these pots (dossier flag 8) | the onset shape (D1); B18's fit rows | DEV; calibration input; **spent** |
| Adams 2008 (F21#0 - #2, F19#6, #7, #12) | fed MGO + Ala (~0.85 M), 100 - 130 C: pyrrole 0.50 / 1.1 / 3.7 mol % of Ala at Ala : MGO 2:1 / 1:1 / 1:3; pyrrole : DMP 38 / 55 / 231; pH 4 acetate vs pH 7 phosphate: pyrrole x80, **DMP x2.1 up at acid**; phosphate vs bicarbonate DMP x1.8, pyrrole x0.36; pH 4 time series 10 / 30 / 60 min (DMP nd / 0.020 / 0.026) | the aminoketone sink and its pH (D2) | DEV(E14); **spent** (EX-052 S5a, S5b, S7; M5 part B) |
| Piloty 1979 (piloty1979 dossier; CF-095) | MGO 312 mM + 13 amino acids 16 mM each, 0.4 M phosphate, 100 C, 3 h: loss Gly 85 / 86 / 82 %, Ala 32 / 27 / 35 %, Phe 45 / 46 / 49 % at pH 4.52 / 6.97 / 9.50; amino-acid analyser (irreversible loss), no products measured | the adduct sink's existence and pH-flat law (D3) | no frozen split; **spent** (CF-095, EX-058) |
| Hofmann 2000 Table 3 (hofmann2000 dossier) | MGO + Phe 0.1 + 0.1 M, 0.5 M phosphate pH 7, 98 C, air: (PA + PAA) / Phe consumed **0.18 at 10 min -> 0.058 at 500 min**, in ONE pot. Integrity flag: the Phe column coincides with Hofmann 2000b's glucose / Phe column | consumption vs Strecker products in one pot (D3) | DEV; calibration input of M3; **seen** (M0 MF-29 / MF-30 read Hofmann 2000) |
| Hofmann 2000 Table 4 (F16#0 - #1) | 0.1 M MGO + Phe, pH 7, 98 C, 60 min: PA 11.1 umol/mmol (argon) | the Strecker product's yield on the same donor (D3, D4) | DEV; calibration input; spent (EX-021) |
| Granvogl 2006 (F16#7) | MGO + Leu / Ile / Phe / Met 5 mM, 170 C, pH 7: aldehyde 2.25 - 12.2 mol %, Strecker amine <= 1.1 % | the same | DEV (DFA); seen |
| **Kanzler 2026** (kanzler2026 dossier; added 2026-11-07) | MGO 100 mM +- one of ten amino acids 100 mM, unbuffered pH 5.0 (initial), 100 C, triplicates: amino acids speed MGO loss x2 - 60 at 30 min (MGO left 32 % with Ala vs 65 % alone) **while Ala is not consumed (100 % left at 120 min, SD 1)**, Lys ~2 %, Phe 37 %; MGO alone second order (k2 1.85e-4 mM^-1 min^-1) | part S: whether the amine-dependent MGO loss is an amine sink (D3'); the MGO supply in fed pots | no frozen split; **seen** (CK-46 lists it SEEN; M4 sweep, first reader) |
| Elfsson 1999 (F19#22) | 5-aminolevulinic acid (an alpha-aminoketone), 50 C: second order, k2 0.069 / 5.53 M^-1 h^-1 at pH 4.81 / 7.42 (x80); none at pH 2.35; Ea 43.7 kJ/mol (39.4 - 48.0); EDTA no effect | the only measured aminoketone condensation (an analogue): a band source, never a decisive unit | DEV(E14) |
| Van Lancker 2010 (F20#7 - #8) | glyoxal-only pots: MPZ 1 - 2 % of pyrazine; free Ala vs Lys-Ala | pair supply | DEV |
| Cerny 1994, Weenen 1992 (F20#4, F19#18) | fed MGO + Ala ~5x the hexose level; O2 no effect on methylpyrazines | supply; no O2 term | DEV |
| Sedehi 2013 (F21#12 - #13; keyed table) | MGO + amino acids at 4 - 25 C: methylimidazole ~1 % of MG loss; MG loss > 20x amine loss | the adduct sink's chemistry | **v1 TEST**; keyed held-out: item 7 only |
| Yu 2018 (F19#19) | Glc + Gly 70 - 90 C, one multiresponse fit | barriers | **v1 TEST**: listed, not used |

### 2.3 Trapping

None on disk for an aminoketone. Zhou 2024 traps GO / MGO with OPD **in the ARP runs only** (Fig. 3, figure-only; the
2026-11-04 text said "Fig. 4" and implied the fed pots; corrected: no dicarbonyl was measured in the fed-dicarbonyl
pots). No aminoacetone or aminoacetaldehyde has been measured in any pot (`parameters_pyrazine` names this as
k_cond's missing measurement); the sweep (sec. 2.9) found none.

### 2.4 Measured intermediates and lifetimes

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Zhou 2024 Fig. 3 | GO, MGO, 1-DX, 3-DX, ARP, pH in the **ARP** pots (figure-only) | the donor's lifetime in ARP pots | DEV |
| Kanzler 2026 Fig. 2B | MGO with Ala 100 / 31.8 / 14.0 / 4.3 % at 0 / 30 / 60 / 120 min (100 mM each, pH 5) | MGO's lifetime with an amine (the forced input a fed pot needs) | seen |
| Deng 2022 (deng2022 dossier) | Met + Glc 200 + 200 mM, pH 7.5, 120 C: parent pyrazine 0.09 umol/L in 120 min, ~200x below Zhou's fed-GO pot | sugar pots are supply-limited | DEV (Ho network); M3 phase-2 reserved |
| Zhou 2023 (F19#5, F20#9) | ARP + Cys, pH 6 / 7 / 8: DMP 66 / 106 / 237 ug/L rising while PZ (47 -> 15) and MPZ fall | pair-specific pH response | DEV; calibration input |
| Zhou 2022 (F19#20) | Cys 0 -> 20 mM at pH 8: pyrazine x34; pH 7 -> 8 x14 | cysteine on the supply | DEV |

### 2.5 Mass balance

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Adams 2008 (EX-052 S5a) | pyrrole / (2 DMP) 19 - 116 in all five pH 4 rows; pH 7 rows do not reject | the pyrrole carries most aminoacetone at pH 4, if it forms from FREE aminoacetone (red team: unshown) | spent |
| Hofmann 2000 Table 3 | Strecker products / Phe consumed 0.18 -> 0.058 in one pot | most consumed amine does not give the Strecker products (one laboratory, one pot) | seen; integrity flag |
| Piloty 1979 vs fed-MGO Strecker yields (D3) | 82 - 86 % glycine and 33 % Phe consumed in 3 h, against aldehyde yields of 1 - 12 % per amino acid in fed-MGO pots of two other designs | the same, cross-design | spent |
| Kanzler 2026 | MGO lost with Ala 95.7 mM by 120 min, at most 35.1 mM by the MGO-alone path [dry], yet 0 % Ala lost | the amine-dependent MGO loss consumes no amine | seen |

### 2.6 One factor varied inside one study

| source (record) | factor -> result | bears on | split |
|---|---|---|---|
| Leahy 1989 / 1989a (F19#8 - #11, F20#12) | pH 5 / 7 / 9: k(9) / k(5) 21 - 60; Lys vs Asn x34.9; a_w maximum ~0.75; R 7.5 - 52 | the cascade's pH and amine identity | DEV; calibration input (the four pH ratios are B18 fit rows); spent |
| Zhou 2023c (F19#21; CF-087) | ARP + Gly: PZ 5.97, DMP 1.16, MPZ ND: R < 0.19 | pairing (E13) | DEV; owed to every pairing score |
| Hwang 1995 (F19#23, F20#6) | second amino acid: totals Lys x4.34 ... Arg x0.84 | amine identity | DEV |
| Scalone 2015 (F19#15) | peptides give pyrazines, free amino acids alone none | peptide class | DEV |
| Kanzler 2026 Fig. 2F | MGO + Ala 10 - 100 mM: excess MGO loss ~[Ala]^0.75 - 0.84 | the amine order of MGO's amine-catalysed loss (CF-122) | seen |
| Piloty 1979 Table 4 | pH 4.52 / 6.97 / 9.50 at 100 C: consumption flat (rate ratio 1.0 - 1.2) | the sink's pH law | spent |
| Van Lancker 2012 (F19#13), Shibamoto 1977 (F19#17, F20#13), Low 2007 (F20#2) | dipeptide x20 - 34; O2, metals; glycine flooding | R and supply | **v2 TEST(E13 / E14)**: item 7 only |
| Whitfield 1988 (F30#3) | lecithin lowers methylpyrazine x2.1 | fat (M5) | DEV; spent by M5 |

### 2.7 Fitted time courses and fits

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Jousse 2002 (jousse2002 dossier) | 11-step network; I + I -> pyrazines "Fast" (not fitted); Strecker "fast" | the only source of the instant condensation | no frozen split |
| engine records | B18, B21, EX-021, EX-052, EX-057, EX-058, E4 prereg, CF-136 | sec. 1 table | -- |

### 2.8 Left out, with reasons

Aldehyde-addition routes to ethyl / extra-methyl pyrazines (Amrani F20#0, Cerny 1994, Low 2006; P9): P2a, rejected for
MPZ by EX-052, not M6. Koehler 1969 / 1970, Wong 1988, Bemis-Young 1993, Shibamoto 1976 / 1977 and Lee 2024: the
mixing ratio R (part P), held out or seen-by-designer; item 7. Lipid pyridines (F21#5 - #10): E20. Methoxypyrazines:
not Maillard. Huang 1989 (pyrazine rates above 120 C): glucose-fed (sweep exclusion), not on disk.

**Counts by kind (sources).** Labelling 5 (two held out); fed intermediates 13 (two held out); trapping 0; measured
intermediates 5; mass balance 4; one factor within a study 9 (three held out); fitted 1 plus records.

### 2.9 The citation sweep (sec. 7 item 2)

**The criterion it was keyed to** (the case's decisive question, not its topic). Does the work hold
- (a) an alpha-aminoketone measured or trapped in a heated sugar or dicarbonyl + amino-acid pot while pyrazines form;
- (b) amino-acid consumption AND Strecker products in the SAME fed-dicarbonyl pot;
- (c) a fed GO / MGO + amino-acid pyrazine course with points <= 30 min, or at >= 2 total concentrations at fixed ratio;
- (d) a measured rate constant for aminoketone self-condensation or aminoacetone + MGO (the pyrrole);
- (e) a label that separates a free-aminoketone route from an imine route to the ring;
- (f) a second laboratory's MGO + amino-acid pyrrole yield, or adduct yields with amine consumption?

**Run (2026-11-07).** `scripts/tools/openalex.py sweep` through the wrapper, anonymous, ~1 request per 1.1 s.
- Seeds, forward and backward: Zhou 2024 (10.1021/acs.jafc.4c03706), Adams 2008 (10.1021/jf0726785), Elfsson 1999
  (10.1016/S0928-0987(98)00009-8), Piloty 1979 (10.1007/BF01186494). Forward only: Hofmann 2000 (10.1021/jf990954c),
  Jousse 2002 (10.1111/j.1365-2621.2002.tb08772.x), Kanzler 2026 (10.1016/j.foodchem.2026.151175; no citing work yet).
  All seven resolved; Hofmann returned 170 of the 174 citations OpenAlex counts and Jousse 152 of 155 (no page failed).
- Result: 616 distinct works, 41 already held, **575 unseen**, screened on titles and abstracts (246 have no abstract
  in OpenAlex and were screened on title alone; abstracts read for ~95 chemistry-relevant works).

| decision | n |
|---|---|
| food or matrix product study, no defined fed pot | 213 |
| read individually: does not meet (a) - (f) | 158 |
| biomedical or other field (mostly 5-ALA photodynamic therapy via the Elfsson seed; CO2 capture, atmospheric, tracers) | 108 |
| review, book or chapter without new pot data | 75 |
| computed (DFT, QSPR): excluded by the no-computed-barrier rule | 5 |
| candidates | 16 |

Notable exclusions: Shu 1998 (acetoin, an acyloin, not a fed dicarbonyl); Xia 2022 (concentration dependence, but an
Amadori feed); Liu 2026 and Zhao 2026 (sugar-carbon labels: cannot separate the routes); Hawkins 2018, Gen 2018
(dicarbonyl + ammonium, atmospheric); Marchand 2000 (thiazoles); Huang 1989, Kristianto 2021 / 2025 (glucose-fed,
220 - 280 C); Mendoza-Sarmiento 2023, Zhan 2025 (DFT). A reserve pool of 5-ALA stability papers without abstracts
(McCarron 2004 / 2005, de Blois 2002, Kaliszewski 2007, Donnelly 2005) is listed, not a candidate.

**Candidates, roles fixed before any reading** (first author from the OpenAlex record; P = the chance it holds what
its role needs). None is read; each enters only by amendment, declared before its data are read.

| # | candidate | crit. | role (fixed 2026-11-07) | P |
|---|---|---|---|---|
| 1 | Chuyen 1972, Strecker degradation of alanine with glyoxal (10.1080/00021369.1972.10860387) | b | decisive for part S (D3 within one pot) only if Ala loss AND acetaldehyde / CO2 are tabulated for the same glyoxal pot with initial concentrations | 0.55 |
| 2 | Fujimaki 1971, decarboxylation of amino acids with glyoxal (10.1080/00021369.1971.10860189) | b | band source (Strecker rate across ~20 amino acids), never decisive; same laboratory as #1 | 0.7 |
| 3 | Piloty 1979 part II (10.1007/bf01186495) | b, c | decisive for (b) or (c) only if it prints amino-acid loss with Strecker aldehyde or pyrazines, or an early pyrazine course, in a fed-dicarbonyl pot; same laboratory as the seed (not a second laboratory) | 0.4 |
| 4 | Gonzalez-Rodriguez 2025, phosphate-catalysed non-symmetric pyrazines in water (10.1039/d5gc04772b) | d | band source for k_cond only if it reports a measured aminoacetone dimerisation constant or course with pH and T (its DFT parts excluded); never decisive | 0.35 |
| 5 | Rudroff 2025 (preprint of #4; 10.21203/rs.3.rs-6759688/v1) | d | evidence only; one laboratory with #4 | 0.3 |
| 6 | Butler 1992, non-enzymatic cyclic dimerisation of 5-ALA (10.1016/s0040-4020(01)80465-8) | d | band source (a second laboratory for the ALA analogue), never decisive | 0.5 |
| 7 | Bunke 2000, degradation of 5-ALA (10.1002/1520-6017(200010)89:10<1335::aid-jps11>3.0.co;2-#) | d | band source, never decisive | 0.6 |
| 8 | Novo 1996, chemical instability of 5-ALA (10.1016/1011-1344(96)07285-5) | d | band source, never decisive | 0.4 |
| 9 | Hunter 2005, 5-ALA dimerises via a protonated dihydropyrazine (10.1016/j.abb.2005.03.007) | d | band source, never decisive | 0.3 |
| 10 | Chandrashaker 2012, Knorr vs Fischer-Fink pyrroles in water (10.1016/j.tet.2012.05.080) | d | evidence only (the pyrrole's mechanism) | 0.2 |
| 11 | Baldensperger 2018, alpha-oxoamide AGEs from GO / MGO + Boc-lysine (10.1021/acs.jafc.7b05813) | f | decisive for part S only if adduct yields and lysine loss are in the same incubation (Glomb laboratory) | 0.35 |
| 12 | Henle 1994, imidazolone from Arg + MGO (10.1007/bf01192954) | f | evidence only (protein-bound) | 0.2 |
| 13 | Xing 2023, mechanochemical Strecker of glyoxal + amino acids (10.1016/j.foodchem.2023.138071) | b | evidence only; decisive for (b) only if a heated aqueous control has amino-acid loss plus aldehyde | 0.25 |
| 14 | Guerra 2012, double Schiff bases of 2,3-butanedione with glycine (10.1021/jf303658m) | e | decisive for R1 vs R3 only if its 120 C aqueous pot gives label fractions separating the routes (butanedione, not MGO; same laboratory as #13) | 0.25 |
| 15 | Akagawa 2003, oxidative deamination of benzylamine by glycoxidation (10.1016/s0968-0896(02)00584-9) | b | evidence only (not an amino acid) | 0.3 |
| 16 | Xu 2026, "building-block coupling" in dimethylpyrazine formation (10.1016/j.foodres.2026.118976) | c | decisive for C (onset or order) only if it has a fed MGO / GO + amino-acid course with points <= 30 min or >= 2 concentrations | 0.15 |

**No unseen work points to (a)**, an aminoketone measured in a pot. The strongest leads are #1 for part S and #4 for
k_cond's class. The sweep JSON (scratch, not committed: the download list stays private) is named in the PR.
**P(a candidate, once read, decides any carried pair) ~0.2** (#1 for S0 vs S1 is the main route; nothing decides R1
vs R3 on the candidates as listed).

### 2.10 Spent and seen data (sec. 7 item 2)

**Spent: these data nominated a rival or a part, or were fitted, and are development evidence only (falsify-only).**

| data | spent by | spent for | consequence |
|---|---|---|---|
| Zhou 2024 Table 2 (rates; the zero-order R2) | B18 fit rows; the 2026-11-04 D1 arithmetic (nominated R3 and the k_cond question) | every C-part rival and R0's Strecker constant | falsify-only: FX-C1 (sec. 4.2) |
| Zhou 2024 Table 1, SI Table S1 (labels) | EX-052 S3 (label average) | the pairing rivals (part P) | falsify-only |
| Leahy 1989 / 1989a (rates, pH ratios, distributions) | B18 fit rows (pH slopes); EX-052 S1; EX-057 | E4, R0's pH law, pairing | falsify-only; B18 refit owner-owed |
| Adams 2008 Tables 2 - 3 | EX-052 S5a, S5b, S7; M5 part B | C-part pyrrole sink; R0's pH law on the fed step; buffer | falsify-only |
| Piloty 1979 | CF-095, EX-058 (nominated the adduct sink) | part S | falsify-only (FX-S1 (ii)) |
| Hofmann 2000 Table 4; EX-021 | EX-021 (nominated D4's yield < 1) | R0's yield 1; the Strecker constant (M3) | falsify-only |
| Zhou 2023, Zhou 2023c (CF-087), the UC Davis papers (CF-098) | B18 calibration (Zhou 2023); CF-087 / CF-098 (low R with glycine) | part P | falsify-only |
| B18 T3 / T4 (Leahy totals and ratio), B21 | nominated part L (supply, amine identity) | part L | falsify-only; part L is decided outside M6 |

**Seen but not spent** (engine numbers or M0 readings on record, but they nominated no M6 rival): Hofmann 2000 Table 3
(M0 MF-29 / MF-30 read Hofmann 2000; this case first uses Table 3's consumption column), Kanzler 2026 (CK-46 lists it
SEEN, M4 read it), Granvogl 2006, Zhou 2025b, Zhou 2022, Hwang 1995 / 1995b, Van Lancker 2010, Cerny 1994, Weenen 1992,
Scalone 2015, Deng 2022 (also M3 phase-2 reserved), Jousse 2002, Elfsson 1999 (a band source). Seen data obey the same
rule here: they may falsify a constant-free structure and never credit.

**Held out, item 7 only:** the v2 E13 / E14 TEST cells (Amrani-Hemaimi 1995, Low 2007, Van Lancker 2012, Shibamoto
1977, Zamora 2020, Elmore 1997), the v1 TEST cells (Fang 2009, Yu 2018, De Haan / Sedehi, Zhou 2000), the keyed
held-out Lee 2024 and Sedehi 2013 tables.

**Unspent data on disk that a rival pair could be decided on: none.** Every pot that bears on a carried pair is spent,
seen under the structure it would test, cross-design, or held out. So nothing on disk can credit (sec. 4.3).

### 2.11 M0's groups: what an M6 rival predicts for their members on seen data

Directions are structural or [dry]. No engine was run.

**G6 (dicarbonyl-amine adducts: MF-24, -25, -29; sign check on -34).**

| member | what an M6 rival predicts there | reading |
|---|---|---|
| MF-25 (Piloty: the engine consumes >= 1e3x too little amine) | S1 closes it by construction (its constant is Piloty's lower bound); S0 cannot; S5 closes it only in phosphate / mixtures | M6's own misfit; spent; FX-S1 (ii) tests S0 across laboratories (not falsified under the largest listed noise) |
| MF-24 (MGO self-loss: first order on its Ea floor vs Kanzler's second order, amine-accelerated) | **FX-S2:** the amine-accelerated MGO loss consumes no amine, so it is not S1's adduct; it is S6 (an amine-catalysed, amine-neutral MGO loss), which moves the MGO supply, not the amine budget | S6 belongs to MF-24's decider (`methylglyoxal_self_loss`, M4); M6's design measures MGO in the pot as a forced input instead of predicting it |
| MF-29 (argon PA + PAA: every M3 rival ~5x high) | S1 lowers PA through the amine budget only where the amine is limiting; at Hofmann 2000b's 0.1 M ARP-derived donor it moves PA by kS / (kS + kad) of the amine consumed, i.e. not the aldehyde per Strecker event | no M6 rival predicts MF-29's x5 at a size cited constants fix; MF-29 stays with M3 (CK-01, CK-50) |
| MF-34 sign check (total pyrazines 2.9 dex low) | every S and C rival other than S0 / C0 LOWERS pyrazines (S1 takes amine and dicarbonyl; C1's yield ~0.006 per event) | confirmed: M6 cannot close MF-34; part L carries it |

**G8 (hexose dicarbonyl / osone supply: MF-12, -23, -30, -34).** Part L (L4: supply in water) is the G8 piece for
MF-34. **MF-30** (B18's Strecker constant ~2.3 dex below Hofmann's): **C1 predicts it** in the right direction and
size: at the edge job (a) leaves, C1 needs kS within 0.22 dex of Hofmann's implied rate, i.e. a yield ~10^-2.3 per
Strecker event, which is what MF-30's candidate (a) ("the constant absorbs the supply deficit") would otherwise
attribute to the fit. On spent data this is a consistency, never a credit. MF-12, -23: no M6 rival touches them.

**G2 (the enolisation / pH law: MF-36).** MF-36 (the pyrazine pH law double counts on the Strecker step) is E4's,
decided as a declared band (owner, 2026-11-01). M6's C1 adds a pH law of its own (Elfsson's x80 between pH 7.4 and
4.8 on k_cond); job (a) on Adams' pH 4 / 7 pair (sec. 4.2) shows that C1 then needs the Strecker flux x28 at pH 4
(acetate) over pH 7 (phosphate), against the direction of Leahy's slope: a prediction the designed pot can test at
two held pH values, not a falsification on disk (buffer and pH are confounded in Adams).

**MF-35 (the mixing ratio R).** Its candidate (a), "k_cond at the analogue's measured constant", is **C1n, falsified**
by FX-C1 unless a fast exit is added (then it is C1). Candidates (b) segregated formation and (c) low-level data are
part P's Pseg and Plow (sec. 3), decided by the E13 TEST data in item 7.

## 3. The rivals

The case's rivals R0 - R4 are kept as ids; the comparison runs per part, as in M5, because the parts are decided by
different measurements. Every rival is named from the case (2026-11-04) or from M0's groups, with its mechanism
signature: what it predicts that the others do not, independent of its constants.

**Part C: the aminoketone fate (E14; MF-35 (a)).**

| id | rival (case id) | structure | signature (constant-free) |
|---|---|---|---|
| C0 | R0, R2 (incumbent; "sink only") | instant condensation (k_cond 1e4 L/(mmol min)), 1/2 DMP per aminoketone, no exit | no lag; DMP rate prop. to [A][M]: x3 MGO step -> x3 (+0.48 dex), x4 dilution -> x1/16 (-1.20); free AK ~20 nM (below any trap); DMP per Strecker event 1/2 |
| C1 | R1 (the claim, as it survives job (a)) | a free AK pool; k_cond at Elfsson's class; a first-order exit (the pyrrole: AK + MGO) with **ks >= 0.028 /min** at 20 mM MGO | a lag of ~1.5 / ks (<= ~55 min); in the exit-dominated regime the DMP rate does NOT follow MGO (x3 step -> ~0 dex) and falls x1/16 on x4 dilution; free AK at uM levels; DMP per event ~10^-2 |
| C1n | R1 without a fast exit (MF-35 (a), CK-52) | the same, exit slower than ks_min | **FALSIFIED** (FX-C1): a t^2 onset over hours |
| C3 | R3 (bypass) | the ring forms from two dicarbonyl - amino-acid imines (or the Strecker intermediate) without free AK; the pyrrole a separate MGO-excess product | no free AK; DMP prop. to ([A][M])^2 in the dilute imine limit: x3 MGO -> x9 (+0.95), x4 dilution -> x1/256 (-2.41); no lag beyond imine equilibration |

**Part S: the amine budget (CF-095; G6; MF-25).**

| id | rival (case id) | structure | signature |
|---|---|---|---|
| S0 | R0, R3, R4 as written | the Strecker step is the only amine exit | Strecker events / amine consumed = 1 (less product loss), in every pot |
| S1 | R1, R2 (the claim's sink) | MGO + amino acid -> adduct, irreversible, pH-flat, buffer-flat, k >= 5.6e-3 M^-1 min^-1 at 100 C [dry] | events / consumed = kS / (kS + kad) < 1, the same with and without phosphate, flat in pH |
| S1K | the sink as CK-46 sizes it from Kanzler's MGO loss | the amine-dependent MGO loss is that 1:1 adduct | **FALSIFIED** (FX-S2) |
| S5 | MF-25 (b), CF-095 explanation 2 | the consumption is a phosphate / mixture artefact (Lys / Arg cross-links, 0.4 M phosphate) | events / consumed = 1 in a single-amine, phosphate-free pot; < 1 in 0.4 M phosphate |
| S6 | MF-24 (a), CF-122 (G6) | an amine-catalysed MGO loss that consumes no amine (Kanzler: aldol / melanoidin) | moves MGO, not the amine; with MGO measured as the design's forced input it predicts nothing on part S's units. Carried, decided under `methylglyoxal_self_loss` (M4), not here |

**Part P: pairing and the mixing ratio R (E13; MF-35 (b), (c); EX-052).** P1 random pairing with one encounter constant
(R = 4; the engine since CF-077), P1x pair-specific constants, P1y pair-specific dihydropyrazine fates, P2b an extra
MPZ route (provisional), Pseg segregated formation (h < 0.05; CF-087), Plow low-level data (a measurement rival).
Falsified earlier: P0 (R = 1, retired by CF-077's fix), P2a_aa (13C, two laboratories). Signatures as EX-052's prereg
and CF-087. **Their decisive data are the E13 TEST studies (item 7) and Zhou 2024's raw ion intensities (W6);** M6's
design does not decide part P, and the case says so rather than sizing a design for it.

**Part L: the level (R4; MF-34; G8).** L4 the trunk's GO / MGO supply in water; L4b the amine identity (Lys vs Gly).
Signature: they RAISE sugar-pot pyrazines and leave fed-dicarbonyl pots unchanged; every M6 sink lowers both. Decided
in M2 / M4 / W3 (supply, amine classes); M6 records that its sinks deepen MF-34 (sign check above).

**Why S6 and part P are carried but not sized here.** S6's question is MGO's own loss law (M4's entry); part P's decisive
data exist and are held out. Sizing a design for either here would duplicate a decider already named.

## 4. The decisive data

### 4.1 What cannot decide here

- **Leahy's pH ratios**: B18 fit rows, conditional on CF-077 and drifted (CF-136); sugar + lysine cascades where the
  supply dominates. Not decisive.
- **The four E13 TEST studies**: item 7 only; they measure R (part P), not the sinks. Their scores or written excuses
  are owed before any E13 verdict (EX-058).
- **Cross-design comparisons (D3 across laboratories, D4)**: no measured between-laboratory spread exists for any
  pyrazine-lane observable, so under item 5 they are not decisive; they may falsify only under the largest listed
  noise (sec. 4.2).
- **Absolute pyrazine levels in sugar pots (D6)**: other errors move them (part L, E4, CF-136); not decisive.

### 4.2 The falsification pass (data on disk; falsify-only)

Every number below is from job (a) of [`src/m6_reduced_power.py`](../../../src/m6_reduced_power.py), recorded in
[`kinetic_core_m6_reduced_power.md`](../../../results/validation/kinetic_core_m6_reduced_power.md). [dry] throughout.

**FX-C1: Zhou 2024's onset against a free-aminoketone pool (Part C).** The datum: six zero-order R2 (0.963 - 0.995) of
c = kt + b at t = 0, 30, 60, 90, 120 min. The test statistic is R2 itself (within one run, so run-level terms cancel);
its reference is P(R2 >= printed) at the truth's shape, **maximised over the noise SD (additive and proportional,
0.3 % - 150 % of the final value)**: the least favourable noise for the falsification, since no SD is printed.
- **C1n (no fast exit).** The pool relaxes at tau = 1 / (2 sqrt(F k2)), F = 2 x the measured rate (yield 1 at steady
  state), k2 at Elfsson's fastest corner (Ea 48 kJ/mol, pH 8). tau = 6 - 61 h across the six units, so the 2-h onset
  is ~t^3; noiseless R2 0.821 - 0.822. p max 0.0015 - 0.026, **Holm 0.009 - 0.026 on every unit: falsified.** A lag
  under 10 min needs k2 3.2 - 5.1 dex above Elfsson's fastest corner.
- **C1 (with a first-order exit).** The six R2 survive (every Holm p >= 0.05) for **ks >= 0.028 /min** (lag ~55 min at
  the edge; 20 000 replicates per point; a 8 000-replicate pilot gave 0.023, the lower edge of C1's band in sec. 5.1,
  so the band is 0.08 dex wider than job (a) allows, on the side favourable to C1). At that edge (100 C, pH 8, DMP):
  free AK 51 uM; Strecker flux 1.44 uM/min; kS 3.6e-3 M^-1 min^-1 = B18 + 2.08 dex = Hofmann-implied -0.22 dex; DMP
  per Strecker event 0.0049; exit product ~173 uM by 120 min.
- **C0, C3:** linear onsets; not tested by R2 (they predict what is seen).
- **Caveat:** the heat-up is unprinted; a heat-up adds lag to every rival, so ignoring it is the favourable side for
  C1 (the falsification is conservative).

**FX-S1: the incumbent's "Strecker is the only amine exit" (S0; Part S).** Constant-free: amine consumed = Strecker events.
- (i) **Within one pot**, Hofmann 2000 Table 3: (PA + PAA) / Phe consumed 0.181 at 10 min (-0.74 dex) falling to 0.058
  at 500 min. No replicate is printed; under omega on each of the two assays (sqrt 2 x 0.419 = 0.59 dex) the one-sided
  p is 0.105. Break-even SD 0.45 dex. Integrity flag: the Phe column coincides with Hofmann 2000b's glucose / Phe
  column (dossier caveat), so even this reading is conditional.
- (ii) **Across laboratories**, Piloty's Ala loss (27 % at pH 6.97; 32 / 35 % at 4.52 / 9.50) against the Strecker
  constant Zhou's DMP rate implies under yield 1 (kS 1.75e-5 M^-1 min^-1): even with no MGO loss in Piloty's pot (an
  upper bound on S0's exposure) S0 consumes <= 0.1 %; the gap is 2.51 dex (2.59 / 2.64). Under tau_ABS per laboratory
  and omega per value (SD 1.49 dex) p = 0.046 with no buffer term and 0.196 with phosphate at M5's kappa_P edge
  (x17 at 0.4 M). Break-even between-laboratory SD 0.99 / 0.35 dex.
- **Holm across (i) and (ii): 0.105 / 0.092. S0 is NOT falsified.** It would be if a within-pot replicate SD below
  0.45 dex were printed for Hofmann's ratio, or a measured pyrazine-lane between-laboratory spread below ~1 dex existed.

**FX-S2: "the amine-dependent MGO loss is a 1:1 amine-consuming adduct" (S1K; Part S).** Kanzler 2026, one study,
triplicate pots: MGO lost 95.7 mM by 120 min with Ala; the MGO-alone path (second order, k2 1.85e-4 mM^-1 min^-1) can
take at most 35.1 mM on the with-Ala MGO course; so >= 60.6 mM is amine-dependent and a 1:1 adduct leaves <= 39.4 % of
the Ala. Observed 100 % (SD 1, read +- 3; the read error used as the per-point floor): t(2 df) 20.2, **p 0.0024:
falsified.** The Piloty-sized S1 sink on the same MGO course consumes 1.8 % of the Ala, consistent with Kanzler's ~0 %:
**S1 and Kanzler agree; S1K does not.** The consequence for the network: the amine-dependent MGO loss (S6) and the
amine sink (S1) are two pieces, not one (CK-46's single adduct fails).

**D2, Adams 2008 pH 4 vs 7 (130 C) under C1 (not a falsification).** C1's single pool with Elfsson's pH law (k2 x 1/80 -
1/154 at pH 4) holds Adams' pair only if the aminoketone flux is x28 at pH 4 (acetate) over pH 7 (phosphate) and the
pyrrole's constant x4.4 - 6.2. Buffer and pH are confounded; no constant-free contradiction.

**R0's pH law on the fed Strecker step (EX-052 S7):** at B18's shipped slope (0.58 dex/pH below 7) R0 predicts DMP(4) /
DMP(7) = -1.74 dex against +0.32 (sigma 0.37): falsified at that point. But the slope is fitted-unidentified (CF-136)
and E4's band includes 0, where R0 predicts 0 dex (p 0.38). Under S4 R0 is carried at the band, and E4's outputs keep
their declared band (owner, 2026-11-01).

### 4.3 What survives, and that nothing on disk can credit

| part | falsified on disk | carried | can data on disk credit? |
|---|---|---|---|
| C | C1n | C0 (R0 / R2), C1 (with ks >= 0.028 /min), C3 | no: Zhou is spent, Adams is spent, no aminoketone measured, no dilution or step series |
| S | S1K | S0 (fragile: break-even 0.45 dex within pot), S1, S5, S6 | no: Piloty spent; Hofmann Table 3 seen and flagged; Kanzler measures the MGO side |
| P | P0, P2a_aa (EX-052) | P1, P1x, P1y, P2b, Pseg, Plow | item 7 (E13 TEST) and W6 (Zhou's raw ions) |
| L | -- | L4, L4b | outside M6 |

**S3 consequences (declared, not applied).** No incumbent structure is falsified on disk (S0 carried; C0 carried;
R0's pH law carried at its band). So M6 refuses no shipped output today. What the carried rivals mean for shipped
pyrazines once S3 runs (an S3 case of its own): C1 and C3 lower fed-dicarbonyl pyrazines relative to C0 only if C0's
Strecker constant is held; with each rival's constants from the same spent rows the shipped Zhou-fit levels agree by
construction, and the rivals disagree on the MGO-excess and dilution levers (C1 ~0, C0 x3, C3 x9 on a x3 MGO
step), so **the MGO-excess lever on pyrazines is shown only as "not negative", its size refused** (from no effect to
x9), and the dilution lever's size spans x1/16 - x1/256. S1 vs S0: amine-limited pots (low-amine foods, fed MGO) differ in the amine
budget; S3 widens those by the S1 band (direction shown: S1 lowers both amine and pyrazines).

### 4.4 The designed decider (W6)

One pot class, built so every decisive unit is a **within-run** value (a run-level term cancels, so the absence of a
pyrazine-lane reproducibility term does not block credit): Ala (or one other amino acid) + MGO 20 + 20 mM, 100 C, pH
held by pH-stat at 7.0, a microvessel with a measured heat-up (<= 1 min, recorded), MGO measured in every aliquot
(OPD) as the forced input.

| arm | what is done in the pot | the unit | separates |
|---|---|---|---|
| T (onset) | aliquots at 10 and 40 min (more points optional) | log10 DMP(10) / DMP(40) | C1's lag vs C0 / C3 |
| M (MGO step) | at 40 min MGO raised x3 (dose weighed, +-2 %) | log10 DMP rate 45 - 60 min / 20 - 35 min | C1 (~0) vs C0 (+0.48) vs C3 (+0.95) |
| D (dilution step) | at 40 min the pot diluted x4 with held buffer | log10 rate ratio, as M | C3 (-2.41) vs C0 / C1 (-1.20) |
| K (aminoketone trap) | an aliquot at 40 min into a trap for aminoacetone (e.g. derivatisation + LC-MS; reporting limit <= 0.5 uM; recovery measured if possible, else 0.1 - 1 assumed) | log10 (trapped AK + blank) / reporting limit | C1 (uM) vs C0, C3 (none) |
| S (amine budget) | two pots, one pH-stat without buffer, one in 0.4 M phosphate; at 60 min Ala consumed (amino-acid analysis) and Strecker events (13CO2 from [1-13C]Ala, or acetaldehyde + acetic acid trapped; recovery 0.9 - 1.0) | log10 events / consumed, per pot | S0 (0) vs S1 (< 0 both) vs S5 (0 / < 0) |

Readouts per aliquot: n replicate analyses (n = 2 or 3); SIDA for DMP and the pyrrole; per-reading CV <= 20 %
(0.08 dex) declared as the design target (conservative for SIDA; see sec. 5.3 for why 10 % is not the primary). The pyrrole (2,5-diacetyl-3-methylpyrrole) is read in every aliquot as
evidence (its level is C1's exit product, ~10^2 uM), not as a unit (no rival but C1 predicts it constant-free).

**Minimum deciding design and its power: sec. 5.4** (arm M decides C3; arm S decides part S; C1 vs C0 is not
decided by any design the shared estimator can read).

## 5. Tolerances, measurement model and power

### 5.1 The measurement model per unit (sec. 7 item 5 as adopted; `predictive_densities`)

**On-disk units (falsify-only, sec. 4.2).**

| unit | value level | noise as declared | role |
|---|---|---|---|
| FX-C1 Zhou R2 (x6) | within run (a shape) | no SD printed: the reference maximised over the per-point SD (0.3 - 150 %), additive and proportional | falsification only (spent) |
| FX-S1 (i) Hofmann ratio | within one pot, two assays | no replicate printed: omega 0.419 dex on each assay as the largest listed term; break-even printed | falsification only (seen; flagged) |
| FX-S1 (ii) Piloty vs Zhou | cross-laboratory absolute | tau_ABS 0.964 per laboratory + omega per value (largest listed), buffer factor as a sensitivity | not decisive for credit (no measured pyrazine-lane between-laboratory spread); falsification only |
| FX-S2 Kanzler Ala | within run (relative to t0 on the plate) | run SD of triplicate pots (1 %, n - 1 = 2 df, a Student t), the figure read error (3 %) as the floor | falsification only (seen) |

**Designed units (credit-capable).** Each is a `UnitNoiseModel(level="within_run")` with two terms, declared in
`src/m6_reduced_power.noise_models`:
- `within_run_sd`: the printed within-run replicate SD of the design (n analyses per aliquot), per reading CV 0.08 dex;
  a ratio of two readings carries sqrt(2 / n) x 0.08, a rate ratio (four readings) 2 / sqrt(n) x 0.08; df = values x
  readings x (n - 1); integrated over its inverse-chi-square posterior and redrawn in the null (`replicate_sds`);
- `instrument_floor`: 0.01 dex (a floor, never the whole noise).
- No `run_offset` term (it cancels within run; the module refuses it there), no discrepancy term (none is needed for a
  within-run value), no between-laboratory term, no term transferred from another observable.
- **Forms.** No experimental-fact form is declared: the step doses (+- 0.01 dex), the heat-up (0.8 - 1.2 min), the trap
  blank (0.2 - 1 x the reporting limit), the trap recovery (0.1 - 1) and the event readout's recovery (0.9 - 1.0) are
  measured or bounded design nuisances, drawn as parameters of EVERY rival alike. One modelling choice is declared:
  k2's ALA -> aminoacetone transfer, primary 0 dex, sensitivity -1 / +1 dex (required agreement). **The sensitivity
  was NOT run** (about 3 h of compute per pair; listed in the freeze record).
- **Bands (S4-style, declared with items 1 - 5).** C1: log kS within +-0.8 dex of Hofmann's implied rate (job (a));
  ks in [0.023, 1] /min (lower edge the job (a) pilot; the full job (a) gives 0.028, see sec. 4.2; upper edge
  unsourced). S1 / S5: kad / kS in [1, 30] (Hofmann 0.06 - 0.18,
  Piloty). C0 and C3 predict every C unit structurally (their constants cancel in within-run values).
- **Precision.** The density estimator's precision requirement is the shared floor (ESS >= 100), declared here; the
  scoring ensembles are drawn from a defensive proposal (0.9 prior + 0.1 edge-weighted, weights p / q passed to the
  judge) so that a wrong rival's posterior, pinned at a band edge or a nuisance tail, keeps an effective sample.

### 5.2 Floors

Instrument floor 0.01 dex on every designed reading (declared). The trap's reporting limit 0.5 uM (a design
requirement, not a data point). On disk: Kanzler's 3 % read error; no other floor enters.

### 5.3 The reduced-model power check (BUILT AND RUN, 2026-11-07)

Record: [`kinetic_core_m6_reduced_power.md`](../../../results/validation/kinetic_core_m6_reduced_power.md) (numbers in
the `.json`), written by [`scripts/generators/m6_reduced_power_check.py`](../../../scripts/generators/m6_reduced_power_check.py)
from [`src/m6_reduced_power.py`](../../../src/m6_reduced_power.py). SYNTHETIC data only; no engine; niced, one thread,
<= 2 processes.

How it runs the adopted rule:
- per pair (a, b), scoring ensembles of each rival (defensive proposal), held-out truths from the prior, a parametric
  bootstrap null of >= 1000 replicates per rival (`build_null_pool`: decisive data and replicate SDs redrawn; no spent
  rows and no forced inputs in the reduced model, declared), then `judge_from_draws` per synthetic set;
- **power** = P(credit the generating rival), each direction; **false credit** = P(credit the other); both with exact
  (Clopper-Pearson) 95 % bounds;
- **the misspecified false-credit gate** (M1.1): the generator adds, beyond the scoring noise, an extra per-point term
  equal to the declared per-point SD (the pyrazine lane has no within-laboratory term; a run offset cancels in every
  unit by construction); its false credits k / n feed `judge_pair` as the gate, which needs the Wilson lower end <= 0.05
  and n large enough to fail it;
- the arms' necessity by leave-one-arm-out and single-arm variants.

Declared limits: the rivals' structures are exact in the mini-model (powers are upper bounds); MGO is held at its
charge (the design measures it: with a measured course the null must redraw it); the spent rows enter only through
C1's band, not as conditioning rows; the heat-up is measured; the trap reading is modelled with a blank and an
uncalibrated recovery.

### 5.4 Results: power, false credit, the arms that are necessary

**Run 2026-11-07** (two processes, niced, one thread; 50 synthetic sets per truth, 20 000 scoring draws per rival,
1000 null replicates per rival with re-estimation on 60 000-draw ensembles, per-reading CV 0.08 dex). P = P(credit the
generating rival), 95 % Clopper-Pearson in brackets. "Gate" = the M1.1 misspecified false credit (k / 100, Wilson lower
end must not exceed 0.05): it passes in **every** cell, at x0, x0.5 and x1 of the per-point SD (worst Wilson lower end
0.041). "Unresolved" = the larger null's unresolved fraction after re-estimation; at >= 0.05 the bound makes the rival
tested against that null uncreditable.

| pair | design | n | P(first) | P(second) | false credit (CP upper) | unresolved | reading |
|---|---|---:|---|---|---|---:|---|
| C3 vs C0 | M | 2 | 0.94 [0.83, 0.99] | 0.98 [0.89, 1.00] | 0.06 (0.17) | 0.003 | **decides** |
| C3 vs C0 | M | 3 | 0.96 [0.86, 1.00] | 0.96 [0.86, 1.00] | 0.02 (0.11) | 0.000 | decides |
| C3 vs C0 | D | 2 | 0.92 [0.81, 0.98] | 0.94 [0.83, 0.99] | 0.00 (0.07) | 0.030 | decides (alternative to M) |
| C3 vs C0 | M+D / M+D+K | 2 | 0.00 | 0.88 / 1.00 | 0.00 (0.07) | 0.053 / 0.114 | **adding arms kills C3's credit** (unresolved >= 0.05) |
| C3 vs C1 | M | 2 | 0.94 [0.83, 0.99] | 0.98 [0.89, 1.00] | 0.02 (0.11) | 0.006 | **decides** |
| C3 vs C1 | M | 3 | 0.96 [0.86, 1.00] | 0.96 [0.86, 1.00] | 0.00 (0.07) | 0.005 | decides |
| C3 vs C1 | D | 2 | 0.96 [0.86, 1.00] | 0.90 [0.78, 0.97] | 0.02 (0.11) | 0.033 | decides (alternative) |
| C3 vs C1 | M+D / M+D+K | 2 | 0.00 / 0.00 | 0.98 / 0.00 | 0.00 | 0.057 / 0.455 | not creditable / not computable |
| C1 vs C0 | M | 2 | 0.74 [0.60, 0.85] | 0.42 [0.28, 0.57] | 0.10 (0.22) | 0.000 | ranks, does not decide |
| C1 vs C0 | M | 3 | 0.70 [0.55, 0.82] | 0.62 [0.47, 0.75] | 0.04 (0.14) | 0.000 | ranks, does not decide |
| C1 vs C0 | M+D | 2 | 0.62 [0.47, 0.75] | 0.50 [0.36, 0.64] | 0.04 (0.14) | 0.001 | ranks |
| C1 vs C0 | D | 2 | 0.16 | 0.00 | 0.06 | 0.000 | blind (both -1.20 dex) |
| C1 vs C0 | M+D+K | 2 | **0.00** | 0.66 [0.51, 0.79] | 0.00 (0.07) | **0.232** | **C1 uncreditable: the null under C0 is 23 % unresolved** |
| S1 vs S0 | plain + phosphate | 3 | 0.96 [0.86, 1.00] | 0.92 [0.81, 0.98] | 0.00 (0.07) | 0.044 | **decides** |
| S1 vs S0 | plain + phosphate | 2 | 0.34 [0.21, 0.49] | 0.90 [0.78, 0.97] | 0.04 (0.14) | 0.049 | S1 nearly uncreditable at n = 2 |
| S1 vs S0 | plain / phosphate only | 2 | 0.90 / 0.94 | 0.92 / 0.90 | 0.06 / 0.04 | 0.031 / 0.020 | decides this pair alone |
| S1 vs S5 | plain + phosphate | 3 | 0.98 [0.89, 1.00] | 0.98 [0.89, 1.00] | 0.00 (0.07) | 0.023 | **decides** |
| S1 vs S5 | plain only / phosphate only | 2 | 0.96 / 0.00 | 0.96 / 0.00 | 0.04 / 0.00 | 0.030 / 0.000 | the plain pot is necessary |
| S5 vs S0 | plain + phosphate | 3 | 0.98 [0.89, 1.00] | 0.86 [0.73, 0.94] | 0.00 (0.07) | 0.017 | **decides** |
| S5 vs S0 | plain only / phosphate only | 2 | 0.00 / 0.94 | 0.00 / 0.96 | 0.00 / 0.04 | 0.000 / 0.022 | the phosphate pot is necessary |

**The minimum deciding design (for what is decidable).**
- **Arm M alone** (one pot: Ala + MGO 20 + 20 mM, 100 C, pH-stat 7.0, MGO raised x3 at 40 min, DMP by SIDA in
  aliquots over 20 - 35 and 45 - 60 min, n = 2 analyses of each aliquot, the within-run SD estimated from replicate
  aliquots with independent quench and work-up): **C3 vs C0 and C3 vs C1, P >= 0.94 each way**. D is an equivalent
  alternative (>= 0.90), not an addition: combining arms raises the nulls' unresolved fraction past 0.05 and blocks
  C3's credit.
- **Arm S** (two pots, no buffer and 0.4 M phosphate, n = 3): **S0 / S1 / S5 pairwise, P 0.86 - 0.98**; both pots are
  necessary (each single pot is blind to one pair), and n = 2 is not enough for S1 vs S0 (0.34).
- Three pots in all, each with its MGO course read. Arm K (the aminoketone trap) and arm T (onset) stay in the design
  as evidence (K is the direct measurement the decisive-measurement entry names), but under the shared estimator they
  cannot credit.

**What the check found about the rule (for the judge's writer).** Every failure of power above is a density the shared
estimator does not resolve, not a lack of information: a wrong rival pinned at a band corner (C1 on C0-like data
needs the condensation-dominated corner; C3 on C0 data needs a ~10 SD dose excursion) has PSIS k-hat >= 0.7 even at
60 000 draws, and the bound then sets those replicates at +- infinity. Two consequences: (i) adding a decisive arm can
LOWER power, which the rule did not anticipate; (ii) the estimator also flags a rival with deterministic predictions
(identical log-weights) as k-hat = inf (`importance_lppd(np.zeros(n))` is "unaccepted"), which is why every rival here
carries its measured design nuisances. Making K computable needs a quadrature or SMC density path in
`build_null_pool` / `score_from_draws` for rivals with <= 3 effective coordinates (C0, C3 and C1 here): **~4 - 8 h of
work in a reviewed shared module, plus ~3 h of compute to re-run the C pairs.**


## 6. What M6 needs from the one network (W3)

- **Species:** AKG, AKM (aminoketones) as real pools, the pyrazines (PZ, MPZ, 2,5- and 2,6-DMP), the pyrrole, and
  MGO- / GO-amino-acid adduct pools (free amino acid and protein-bound: B20's CEL is the protein analogue); for C3 a
  dihydropyrazine from imines; **for S6 an amine-catalysed MGO loss that consumes no amine, kept apart from S1's
  amine-consuming adduct** (FX-S2).
- **One Strecker step with two products:** `r_go_ak` / `r_mgo_ak` emit the Strecker aldehyde (M3) and the
  aminoketone from one event, per amino-acid class (M4's classes), so M3 and M6 cannot fit the same step twice.
- **Supply from the composed trunk:** GO and MGO from the hexose and pentose modules in water, not only B13's glass.
- **Per-step pH laws** that keep supply and step apart (E4), and hooks for M5's buffer factor (Adams: phosphate vs
  bicarbonate x1.8) and a_w term (Leahy's maximum).
- **The forced-input driver** for the fed-dicarbonyl pots (MGO's measured course clamped; the designed pot's steps).
- **Isotopologue pools** (15N amino acid, 13C sugar, [1-13C] amino acid for 13CO2) so the labelling rows are observables.
- **Cost of running the carried rivals (S3):** C1, C3, S1, S5 as switches in the NumPy engine; not built; ~1 day.

## Open questions for the build

1. **Decide M3 phase 2 and M6 together.** D4 ties the Strecker constant to the pyrazine yield; C1 now predicts EX-021's
   gap. The owner sees any move of shipped pyrazines first (M3 case, open question 5).
2. **Which k_cond class?** FX-C1 rules out Elfsson's constant without a fast exit; with the exit, k_cond barely moves
   the designed units. Do not fit k_cond on pyrazine levels. Candidate #4 of the sweep is the only lead for a measured
   aminoacetone constant (a band source, by amendment).
3. **The adduct sink's products and stoichiometry.** Piloty measures loss only; Sedehi (held out) is at 4 - 25 C.
   FX-S2 separates S1 from S6; the build carries both.
4. **Part L belongs elsewhere.** Re-word the plan's M6 row so M6 is not expected to fix the level gap, and so it no
   longer says the analogue "rejects the finite condensation as built" without the exit (sec. 1; for the plan's writer).
5. **Re-fit B18 under CF-077 and CF-136** (owner decision, EX-057) before any M6 comparison uses the Leahy rows.
6. **Score or excuse the four E13 TEST studies** in item 7, with and without UC Davis (CF-098), reporting Lee 2024 at
   face value and with MPZ / 10 (CF-097).
7. **W6 requests:** the designed pot (sec. 4.4); Zhou 2024's supplementary ion data (EX-052); any aminoketone
   measurement; Hofmann 2000 Table 3's Phe column re-read against the PDF (the FX-S1 (i) integrity flag; the owner's
   PDF).

## Freeze record

**NOT FROZEN (2026-11-07).** Items 1 - 4 are complete (rivals named with signatures from the case and M0's groups;
every source classified spent / seen / held out / decisive; the falsification pass with numbers; the measurement model
in the adopted judge's terms). Item 5 is sized. What blocks the freeze, exactly:

1. **The claim vs the incumbent on part C (C1 vs C0) has no deciding design that is computable under the adopted
   judge.** The best computable design (arm M, n = 3) separates them at 0.62 [0.47, 0.75]. The arm that separates them
   structurally (K, the aminoketone trap: uM vs none) leaves the null under C0 23 % unresolved, so C1 can never be
   credited. A W6 request sized now would either omit the claim's own decider or ask for an arm the judge cannot read.
   **Unblocks:** a quadrature / SMC density path for low-dimensional rivals in `predictive_densities` (its own
   reviewed PR; ~4 - 8 h), then re-running the C pairs (~3 h of compute).
2. **The k2 transfer sensitivity (-1 / +1 dex) is declared but not run** (~3 h of compute). It moves the condensation
   share, which is exactly what lets C1 mimic C0 on arm M, so it may change item 1's numbers.
3. **No independent premortem / confirmation has read items 1 - 5** (the programme's practice before M4's and M5's
   freezes). Not run here, by scope.
4. **C1's band lower edge (0.023 /min) is the job (a) pilot's, 0.08 dex below the full job (a)'s 0.028.** Declared,
   conservative for C1; re-run with 0.028 when item 1 is re-run.

What would NOT change on unblocking: the falsifications (FX-C1, FX-S2), S0 not falsified (FX-S1), the rivals, the
spent / seen classification, arm M's and arm S's sizes for the pairs they decide.

When frozen, the freeze sha goes in the PR body (a commit cannot hold its own sha).


## Revision log

### Changelog: 2026-11-07, items 1 - 5 brought to a freezable state (branch `m6-case-ready`)

- Header, rules and seen-by lines added (the adopted item 5; no H0 in the lane; CF-136's tiers).
- Sec. 1: the claim's exit stated as job (a) requires it (ks >= 0.028 /min); the 2026-11-04 reading of D1 corrected.
- Sec. 2: Kanzler 2026 and Hofmann 2000 Table 3 added (parts S); Zhou 2024's "Fig. 4 GO / MGO in the fed pots"
  corrected (ARP pots, Fig. 3); sec. 2.9 (the sweep), 2.10 (spent / seen) and 2.11 (M0's groups) added.
- Sec. 3: rivals per part (C, S, P, L) with signatures; C1n, S1K, S5, S6 and part P named from M0 and the records.
- Sec. 4: the falsification pass (FX-C1, FX-S1, FX-S2, D2, R0's pH law); the designed decider written.
- Sec. 5: the measurement model in the adopted judge's terms; the reduced power check built and run; NOT FROZEN
  (freeze record).
