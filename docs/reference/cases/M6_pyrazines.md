# Case M6: pyrazines

*Plan step M6 ([STATUS_AND_PLAN sec. 5, W2](../../guides/STATUS_AND_PLAN.md); errors E13 (F19-a, CF-077) and E14
(F19-b, F21-a, F20-a) of [ENGINE_VS_LITERATURE](../pathway_atlas/ENGINE_VS_LITERATURE.md) sec. 3, amended in secs. 5.1
and 7 (E13 / E14 and the pyrazine rivals), with CF-095 (the dicarbonyl-amine adduct sink) and E4 (decided: a declared
band, owner 2026-11-01) as neighbours. Items 1 - 5 of sec. 7, written 2026-11-04 before any result of the change
exists. No engine was run for this case. Items 6 - 9 are appended later and never edited into this part. M6 shares the
Strecker steps `r_go_ak` / `r_mgo_ak` with case [M3](M3_strecker_aldehydes.md) phase 2.*

**Seen by the author.** Engine outputs read, all on development or calibration rows: the B18 outcome
(`results/validation/kinetic_core_b18_prereg.md` sec. 6: Leahy's 95 C distribution and total, apparent barriers);
EX-057's moves of Leahy's pH-ratio rows and of PH-06 / TEMP-04 under CF-077; CF-095's probe of Piloty's pot (the engine
consumes 1.9e-5 to 1.6e-3 of the glycine); EX-052's scores (S1 - S7). Measured values only for the v2 E13 / E14 TEST
cells (F19#0, #1, #13, #17; F20#1, #2, #3, #13; F21#8, #10), for the keyed held-out Lee 2024 and Sedehi 2013 tables
as quoted in CF-097 / CF-095, and for the UC Davis papers (Shibamoto 1976, Wong 1988, Bemis-Young 1993) quoted in
CF-087 / CF-098. No engine output on any of them. The four E13 TEST studies were neither scored nor excused by EX-052
or EX-058; that debt stands.

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
- Piloty 1979, Koehler 1969 / 1970, Jousse 2002 and Bemis-Young 1993 have no atlas cell. Piloty, Koehler and Jousse sit
  in no frozen split; Bemis-Young is UC Davis, an E13 TEST laboratory.

## 1. The claim

Two changes downstream of the Strecker step, plus the condition the plan row does not name.

| part | step(s) | today (B60 + CF-077, read from the code) | M6 |
|---|---|---|---|
| S. dicarbonyl-amine adduct sink | new: MGO (and GO) + amino acid -> non-Strecker adduct, irreversible, pH-flat | the only MGO + amine step is the Strecker `r_mgo_ak` (k ~3e-8 L/(mmol min) at 100 C); `r_mgo_mel` consumes no amine ([network.py](../../../src/kinetic_core/network.py) ~180, ~359 - 365) | on, its consumption law from Piloty 1979 (CF-095's candidate) |
| C. finite aminoketone fate | `r_akg_pz`, `r_akm_dmp`, `r_ak_mpz` on `k_cond`; new: aminoacetone + MGO -> pyrrole | `k_cond` = 1e4 L/(mmol min), Ea 0, no pH ("declared fast" from Jousse 2002; [parameters_pyrazine.py](../../../src/kinetic_core/parameters_pyrazine.py) ~53, ~122 - 124); cross pair x2 since CF-077 (PR #37); no pyrrole, no dihydropyrazine | `k_cond` from a measured class with its pH law and barrier; the pyrrole as a competing aminoketone sink |
| not claimed | the pyrazine supply (the trunk's GO / MGO in water) and the amine identity | B18: glyoxal reaches a sugar pot only through B13's dry-glass glucosone; glycine stands for every amine | named as rival R4 because the level gap sits there (sec. 4) |

**Earlier waves on these structures** (cited, none re-run):

| structure | wave | outcome | kind |
|---|---|---|---|
| Strecker x2 fitted, condensation declared fast | B18 (SHIPS, two caveats) | Zhou 2024's six rates within 0.07 dex, Leahy's four pH ratios within 0.13; **T3 failed** (Leahy pyrazine : 2,5-DMP 23 : 1, model 6e-8 : 1), **T4 failed** (total 2.9 dex low, glycine for lysine); apparent barriers 300 - 460 kJ/mol vs measured 150 - 177: "the trunk's dicarbonyl SUPPLY" | fit, on fed-dicarbonyl rows |
| the glucosone route to glyoxal | B21 | Leahy total unchanged (-2.8 dex); the trunk's MGO 7.8 mM where Xia's panel tops at 0.5 mM | fit |
| P0 (engine) vs P1, P1x, P2a / b, NC | EX-052 | P0 rejected (Ghent S5a, S7; Reineccius S1; Ho S2); P1 inconclusive; P1x provisional; P2b provisional; P2a_aa rejected by 13C | constant-free signatures |
| the cross-pair factor | EX-057 / CF-077 | fixed (R = 4); the frozen B18 pH slopes are conditional on it (Leahy rows +0.51 dex at the fit commit); a B18 refit awaits the owner | diagnosis |
| Piloty, R pool invariance, form of E5 | EX-058 | Piloty not discriminating for E4 (P 0.85); exposes an adduct sink >= ~1e3x the engine's Strecker consumption (a lower bound); P2b back to ~0.45 | arithmetic + probe |
| the Strecker step's pH law (E4) | E4 prereg | expected REFUSE; declared band | prereg |
| B18's Strecker constant | EX-021 | a pyrazine-yield constant ~2.3 dex below Hofmann's fed MGO + Phe rate | arithmetic |

**Errors stack (CF-107).** B18's T3 / T4 failures are whole-model outcomes on the trunk. Candidate compensating errors:
- the trunk's GO supply (dry-glass glucosone only) and B13's glyoxal sink (98 % of a fed 20 mM GO in 2 h);
- the trunk's MGO level (B21) and E2;
- glycine for lysine (Leahy's amine; Lys >> Asn x35, F19#10);
- the Leahy pH law on the Strecker step (E4), conditional on CF-077;
- k_go_ak's 70x conflict between Zhou 2024 and 2025b (B18 post-run note);
- instant condensation and the missing adduct sink (this case).

The S5a / S7 rejections of P0 in EX-052 are constant-free signatures, not fits, and are stronger evidence.

## 2. All the evidence, by kind

### 2.1 Isotope labelling

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Zhou 2024 (F19#3, F20#10, F03#15) | Ala-xylose ARP + 15N-Ala, 120 C: pyrazine 15N share flat (54.5 -> 54.3 %) over 30 - 120 min; free Ala supplies ~55 % of ring N; MPZ 37 vs 49 % | aminoketones from free amino acid on ARP-derived dicarbonyls; pool mixing | DEV(E13, E14); calibration input |
| Zhou 2025b (F19#2, F20#11) | pyrazine isotopomers close to one well-mixed C2 pool; MPZ M+3 absent at 30 / 60 min, ~4.5 % at 90 / 120 | a slow minor MPZ channel | DEV |
| Hwang 1995b (F21#3) | > 90 % of 1-methylpyrrole N is glycine's | pyrrole N from the amino acid | DEV (Ho network) |
| Amrani-Hemaimi 1995, Low 2007 (F19#0, #1; F20#1 - #3) | MPZ / DMP ring carbons unlabelled from labelled amino acids; DMP 25 - 30 % from labelled Gly | aminoketone pairing; aldehyde addition | **v2 TEST(E13 / E14)**: item 7 only |
| Fang 2009, Mundt 2004 (F31#0 - #3) | 27 - 33 % vs ~2/3 of glycine decarboxylated in melanoidins | total Strecker-type flux into polymer | Fang **v1 TEST**: listed, not used; Mundt DEV |

### 2.2 Fed intermediates

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Zhou 2024 (F19#4; zhou2024 dossier) | Ala + GO or + MGO 20 + 20 mM, pH 8, 100 / 110 / 120 C: pyrazine 0.0279 / 0.0791 / 0.1507, 2,5-DMP 0.0035 / 0.0100 / 0.0230 umol/(L min); **zero-order c = kt + b over 0 - 120 min, R2 > 0.96** | the onset shape (D1); B18's fit rows | DEV; calibration input |
| Adams 2008 (F21#0 - #2, F19#6, #7, #12) | fed MGO + Ala (~0.85 M), 100 - 130 C: pyrrole 0.50 / 1.1 / 3.7 mol % of Ala at Ala : MGO 2:1 / 1:1 / 1:3; pyrrole : DMP 38 / 55 / 231; pH 4 acetate vs pH 7 phosphate: pyrrole x80, **DMP x2.1 up at acid**; phosphate vs bicarbonate DMP x1.8, pyrrole x0.36; Pro lowest | the aminoacetone sink and its pH (D2) | DEV(E14) |
| Piloty 1979 (piloty1979 dossier; CF-095) | MGO 312 mM + 13 amino acids 16 mM each, 0.4 M phosphate, 100 C, 3 h: glycine loss 85 / 86 / 82 % at pH 4.5 / 7.0 / 9.5; Phe 33 %, Leu 26 %; amino-acid analyser (irreversible loss), no products measured | the adduct sink's existence and pH-flat law (D3) | no frozen split |
| Hofmann 2000 Table 4 (F16#0 - #1) | 0.1 M MGO + Phe, pH 7, 98 C, 60 min: PA 11.1 umol/mmol (argon) | the Strecker product's yield on the same donor (D3, D7) | DEV; calibration input |
| Granvogl 2006 (F16#7) | MGO + Leu / Ile / Phe / Met 5 mM, 170 C, pH 7: aldehyde 2.25 - 12.2 mol %, Strecker amine <= 1.1 % | the same | DEV (DFA) |
| Elfsson 1999 (F19#22) | 5-aminolevulinic acid (an alpha-aminoketone), 50 C: second order, k2 0.069 / 5.53 M^-1 h^-1 at pH 4.81 / 7.42 (x80); none at pH 2.35; Ea 43.7 kJ/mol; EDTA no effect | the only measured aminoketone condensation (an analogue) | DEV(E14) |
| Van Lancker 2010 (F20#7 - #8) | glyoxal-only pots: MPZ 1 - 2 % of pyrazine; free Ala vs Lys-Ala | pair supply | DEV |
| Cerny 1994, Weenen 1992 (F20#4, F19#18) | fed MGO + Ala ~5x the hexose level; O2 no effect on methylpyrazines | supply; no O2 term | DEV |
| Sedehi 2013 (F21#12 - #13; keyed table) | MGO + amino acids at 4 - 25 C: methylimidazole ~1 % of MG loss; MG loss > 20x amine loss | the adduct sink's chemistry | **v1 TEST**; keyed held-out: item 7 only |
| Yu 2018 (F19#19) | Glc + Gly 70 - 90 C, one multiresponse fit | barriers | **v1 TEST**: listed, not used |

### 2.3 Trapping

None on disk for an aminoketone. Zhou 2024 traps GO / MGO with OPD (Fig. 4, figure-only); no aminoacetone or
aminoacetaldehyde has been measured in any pot (`parameters_pyrazine` names this as k_cond's missing measurement).

### 2.4 Measured intermediates and lifetimes

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Zhou 2024 Fig. 4 | GO and MGO time courses in the fed pots (figure-only) | the donor's lifetime against B13's sink | DEV |
| Deng 2022 (deng2022 dossier) | Met + Glc 200 + 200 mM, pH 7.5, 120 C: parent pyrazine 0.09 umol/L in 120 min, ~200x below Zhou's fed-GO pot | sugar pots are supply-limited | DEV (Ho network) |
| Zhou 2023 (F19#5, F20#9) | ARP + Cys, pH 6 / 7 / 8: DMP 66 / 106 / 237 ug/L rising while PZ (47 -> 15) and MPZ fall | pair-specific pH response | DEV; calibration input |
| Zhou 2022 (F19#20) | Cys 0 -> 20 mM at pH 8: pyrazine x34; pH 7 -> 8 x14 | cysteine on the supply | DEV |

### 2.5 Mass balance

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Adams 2008 (EX-052 S5a) | pyrrole / (2 DMP) 19 - 116 in all five pH 4 rows; pH 7 rows do not reject | the pyrrole carries most aminoacetone at pH 4, if it forms from FREE aminoacetone (red team: unshown) | DEV(E14) |
| Piloty 1979 vs fed-MGO Strecker yields (D3) | [dry] 82 - 86 % glycine and 33 % Phe consumed in 3 h, against aldehyde yields of 1 - 12 % per amino acid in fed-MGO pots of two other designs | most amine consumed by MGO does not give the Strecker products | cross-design |

### 2.6 One factor varied inside one study

| source (record) | factor -> result | bears on | split |
|---|---|---|---|
| Leahy 1989 / 1989a (F19#8 - #11, F20#12) | pH 5 / 7 / 9: k(9) / k(5) 21 - 60; Lys vs Asn x34.9; a_w maximum ~0.75; R 7.5 - 52 | the cascade's pH and amine identity | DEV; calibration input (the four pH ratios are B18 fit rows) |
| Zhou 2023c (F19#21; CF-087) | ARP + Gly: PZ 5.97, DMP 1.16, MPZ ND: R < 0.19 | pairing (E13) | DEV; owed to every pairing score |
| Hwang 1995 (F19#23, F20#6) | second amino acid: totals Lys x4.34 ... Arg x0.84 | amine identity | DEV |
| Scalone 2015 (F19#15) | peptides give pyrazines, free amino acids alone none | peptide class | DEV |
| Van Lancker 2012 (F19#13), Shibamoto 1977 (F19#17, F20#13), Low 2007 (F20#2) | dipeptide x20 - 34; O2, metals; glycine flooding | R and supply | **v2 TEST(E13 / E14)**: item 7 only |
| Whitfield 1988 (F30#3) | lecithin lowers methylpyrazine x2.1 | fat (M5) | DEV |

### 2.7 Fitted time courses and fits

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Jousse 2002 (jousse2002 dossier) | 11-step network; I + I -> pyrazines "Fast" (not fitted); Strecker "fast" | the only source of the instant condensation | no frozen split |
| engine records | B18, B21, EX-021, EX-052, EX-057, EX-058, E4 prereg | sec. 1 table | -- |

**Left out, with reasons.** Aldehyde-addition routes to ethyl / extra-methyl pyrazines (Amrani F20#0, Cerny 1994, Low
2006; P9): P2a, rejected for MPZ by EX-052, not M6. Koehler 1969 / 1970, Wong 1988, Bemis-Young 1993, Shibamoto 1976 /
1977 and Lee 2024: the mixing ratio R (E13), held out or seen-by-designer. Lipid pyridines (F21#5 - #10): E20.
Methoxypyrazines: not Maillard. Huang 1989 (pyrazine rates above 120 C): on B18's fetch list, not on disk.

**Counts by kind (sources).** Labelling 5 (two held out); fed intermediates 11 (two held out); trapping 0; measured
intermediates 4; mass balance 2; one factor within a study 7 (three held out); fitted 1 plus records.

## 3. The rivals

| id | explanation |
|---|---|
| R0 | **current engine**: Strecker rate-determining, condensation instant (yield 1 per aminoketone pair), no adduct sink, no pyrrole, Leahy's pH law on the Strecker step (E4 band) |
| R1 | **the claim**: + the adduct sink (MGO / GO + amino acid, pH-flat) + a finite k_cond from a measured class (Elfsson's pH law and barrier) + aminoacetone + MGO -> pyrrole competing for the aminoketone |
| R2 | the adduct sink only; condensation stays instant |
| R3 | the ring does not pass through free aminoketone: a dihydropyrazine forms from the dicarbonyl - amino-acid imine (or the Strecker intermediate), so neither k_cond nor a free-aminoketone pyrrole sink sets the pyrazine yield; the pyrrole is a separate MGO-excess product |
| R4 | **supply**: the 2.9-dex level gap is the trunk's GO / MGO supply in water and the amine identity (B18's own diagnosis); M6's sinks are second-order corrections |

## 4. The decisive data

| # | measurement (record) | where the rivals differ | separates |
|---|---|---|---|
| D1 | **onset of fed-dicarbonyl pyrazines** (Zhou 2024: zero-order rise, R2 > 0.96 over 0 - 120 min) | R1 with Elfsson's class: the aminoketone pool takes tau = 1 / (2 sqrt(F k2)) to reach steady state, so DMP rises quadratically after a lag; R0, R2, R3: no lag | **[dry] the ALA class cannot transfer:** at 100 C the aminoacetone flux is F = 2 x 0.0035 = 0.007 umol/(L min) and Elfsson's k2 extrapolates to 0.76 M^-1 min^-1 (pH 7.5), so tau ~ 115 h (the elfsson dossier's own figure); a lag under 30 min needs k2 >= 1 / (4 tau^2 F) ~ 4e4 M^-1 min^-1, ~5e4x the analogue. A measured constant then rejects R1-as-built unless the ring bypasses free aminoketone (R3). Zhou's linearity is partly in-sample (B18 fits its rates) |
| D2 | **pyrrole vs DMP at one supply** (Adams 2008: pyrrole : DMP 38 / 55 / 231 as MGO rises; pH 4 vs 7 pyrrole x80 and DMP x2.1 up) | R1: pyrrole and DMP compete for one free aminoacetone pool, so at fixed supply DMP must fall when the pyrrole rises; at pH 4 DMP rises only if the supply rises more than x2.1 times the pyrrole's diversion. R3: the two products move independently. R0: no pyrrole; DMP(4) / DMP(7) set by the Leahy slope (EX-052 S7: P0 -1.74 dex, observed +0.32) | R0 from R1 / R3 (done, EX-052); R1 vs R3 needs aminoacetone measured in the same pot, which no source on disk has |
| D3 | **amine consumption vs Strecker products on fed MGO** (Piloty pH-flat 82 - 86 %; Hofmann PA 1.1 % per Phe; Granvogl 2 - 12 mol %) | R0, R4: the Strecker step is the only amine exit, so consumption ~ products; R1, R2: most consumption is adduct, pH-flat | R1 / R2 from R0 on the amine budget. Cross-design (MGO 312 mM vs 0.1 M vs 5 mM; 100 vs 98 vs 170 C), so tolerance H0; no single pot measures both |
| D4 | **pyrazine yield per Strecker event** (EX-021: B18's constant ~2.3 dex below Hofmann's fed MGO + Phe rate) | R0: yield 1 per aminoketone pair, so B18's constant is the Strecker rate; R1, R2, R3: the yield is < 1, and if M3 moves `r_mgo_ak` to Hofmann's rate the yield must be ~10^-2.3 to keep Zhou's pyrazine rates | R0 from the rest, jointly with M3. Phe vs Ala, pH 7 vs 8, 0.1 M vs 20 mM: tolerance H0 |
| D5 | **dilution** of a fed-dicarbonyl pot at fixed ratio | R1 (second-order condensation against first-order sinks): pyrazine per Strecker event falls on dilution; R0, R2, R3: unchanged | R1 from the rest. **No such series on disk** (the designed E13 - E14 pot, sec. 3.1, asks for it) |
| D6 | **the level in sugar pots** (B18 T3 / T4: Leahy pyrazine : DMP 23 : 1 and total; Deng 2022 ~200x below fed GO) | R1, R2 lower pyrazines further; only R4 raises them | the plan row's premise (sec. 4 note) |

**The plan's decisive data, checked.**
- "Leahy's pH ratios": Leahy 1989 / 1989a are calibration inputs; the four pH ratios are B18's fit rows, exact in
  every rival after a refit, and conditional on CF-077 (EX-057). They cannot decide. They are also sugar + lysine
  cascades, where the supply dominates (B18's apparent barriers 300 - 460 kJ/mol).
- "The four held-out E13 studies": v2 TEST(E13) (Amrani-Hemaimi 1995, Low 2007, Van Lancker 2012, Shibamoto 1977),
  item 7 only. They measure R (pairing), not the adduct sink or k_cond. Their scores or written excuses are owed before
  any E13 verdict (EX-058).

**What the data cannot decide.** R1 against R3 (does the ring pass through free aminoketone?) needs an aminoketone
measured in a pot; none exists. D3 and D4 are cross-design. No dilution series (D5). The pyrrole sink rests on one
laboratory (Ghent) at ~0.85 M.

## 5. Tolerances from the data

| comparison | noise as the source states it | where stated |
|---|---|---|
| D1 Zhou 2024 | zero-order R2 per temperature (pyrazine 0.979 / 0.992 / 0.963); pH 8.0 +- 0.02; SPME calibration R2 >= 0.992; replicates not stated in the dossier | zhou2024 dossier |
| D1 Elfsson 1999 | k2 SE 6.3e-3 (pH 4.81, n 9) and 0.58 M^-1 h^-1 (pH 7.42, n 7); Ea slope CI; extrapolation beyond 85 C flagged | elfsson1999 dossier |
| D2 Adams 2008 | +- replicate spread, n not stated; pyrrole reference 78 % pure (<= 0.1 dex) | adams2008 dossier; EX-052 red team sec. 4 |
| D3 Piloty 1979 | amino-acid analysis mean error +-3 %; replicates not stated | piloty1979 dossier |
| D3, D4 Hofmann 2000 Table 4 | replicates not stated; no SD; ratios below ~1.3x directional | hofmann2000 dossier |
| D3 Granvogl 2006 | as printed (granvogl2006 dossier; not re-read here) | -- |
| D6 Leahy 1989 | duplicates early, triplicates later | leahy1989a dossier |
| across laboratories and designs (D3, D4) | between-lab spread from H0 (EX-059, pending) | -- |

## 6. What M6 needs from the one network (W3)

- **Species:** AKG, AKM (aminoketones) as real pools, the pyrazines (PZ, MPZ, 2,5- and 2,6-DMP), the pyrrole, and
  MGO- / GO-amino-acid adduct pools (free amino acid and protein-bound: B20's CEL is the protein analogue); for R3 a
  dihydropyrazine.
- **One Strecker step with two products:** `r_go_ak` / `r_mgo_ak` emit the Strecker aldehyde (M3) and the
  aminoketone from one event, per amino-acid class (M4's classes), so M3 and M6 cannot fit the same step twice.
- **Supply from the composed trunk:** GO and MGO from the hexose and pentose modules in water, not only B13's glass.
- **Per-step pH laws** that keep supply and step apart (E4), and hooks for M5's buffer factor (Adams: phosphate vs
  bicarbonate x1.8) and a_w term (Leahy's maximum).
- **The forced-input driver** for the fed-dicarbonyl pots (Zhou's GO / MGO courses clamped).
- **Isotopologue pools** (15N amino acid, 13C sugar) so the labelling rows are observables.

## Open questions for the build

1. **Decide M3 phase 2 and M6 together.** D4 ties the Strecker constant to the pyrazine yield; the owner sees any move
   of shipped pyrazines first (M3 case, open question 5).
2. **Which k_cond class?** D1 rules out Elfsson's constant as transferable. Declare a band from D1's lower bound, or
   build R3; do not fit k_cond on pyrazine levels.
3. **The adduct sink's products and stoichiometry.** Piloty measures loss only; Sedehi (held out) is at 4 - 25 C.
   Declare whether the adduct returns amine (reversible imines were hydrolysed in Piloty's analysis).
4. **R4 belongs elsewhere.** The level gap is the trunk's supply in water (B18 T3 / T4, B21) and the amine class (M4).
   Re-word the plan's R2 row so M6 is not expected to fix it.
5. **Re-fit B18's pH slopes** under CF-077 (owner decision, EX-057) before any M6 comparison uses the Leahy rows.
6. **Score or excuse the four E13 TEST studies** in item 7, with and without UC Davis (CF-098), reporting Lee 2024 at
   face value and with MPZ / 10 (CF-097).
7. **W6 request:** Zhou 2024's supplementary ion data ([M-H]+ correction; EX-052) and any aminoketone measurement.

Seen by the author: as stated at the top. No engine output on any row named in sec. 4 was computed for this case.
