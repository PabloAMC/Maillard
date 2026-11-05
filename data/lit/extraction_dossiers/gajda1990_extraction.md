# Gajda, Nagy & Burger 1990 — EXTRACTION
DOI: 10.1039/dt9900003155 (not printed in the PDF; Crossref 2026-11-05: "Proton and zinc(II) complexes of 2-(polyhydroxyalkyl)thiazolidine-4-carboxylic acid derivatives", Gajda, Nagy, Burger, J. Chem. Soc. Dalton Trans. 1990, 3155. **Matches the brief and the PDF.**)
### Potentiometric protonation constants, 25.0 °C, I = 0.1 M NaClO4, of ten cysteine + aldose thiazolidine-4-carboxylic acids, **including the D-xylose ring (ligand L3 = TTCA): ring-N log K = 5.32 ± 0.02, carboxyl log K ≈ 1.46**, plus the 2-propyl ring (6.25) and calorimetric ΔH of ring-N protonation for five ligands (−20.0 to −32.3 kJ/mol). Room-temperature NMR shows no acyclic (imine) or free-aldose signal at pH 1.5-8. Zinc(II) complexes (not used here). **No formation K, no rates, one temperature.**

**Read 2026-11-05 from `data/articles/gajda1990.pdf`** (6 pp., Dalton Trans. pp. 3155-3160; RSC scan with a text
layer). Prose and Tables 2-4 from `pdftotext -layout`; **Figure 1 (ligand structures) and Table 1 read from a 200 dpi
render of p. 3156** (scratchpad `img/gaj-2.png`), because the text layer loses every structure and so the ligand key.
The ligand key below is MINE, from the Fischer projections of Fig. 1 against the sugar list in the Experimental and
the elemental analyses of Table 1 (see Flag 1); the paper never prints "L3 = D-xylose" in words.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Department of Inorganic and Analytical Chemistry, A. József University, Szeged (Burger; coordination
chemistry). Ligands made by the method of Bognár et al. (Debrecen; ref 3 includes Györgydeák and Szilágyi, see
`szilagyi1979_extraction.md`). Not in `lab_ids.yml`; proposed id `szeged_burger`. Independent of every Maillard group
in the corpus and of the Warsaw group (Radomski, Karpinski).

## 1. What this paper answers for the model

The question (case M1c sec. 5.2, row "R-K: pKa_T"; CF-076): the ring-nitrogen pKa of TTCA, which sets where a
Kallen-shaped K_app(pH) turns over. Until now the corpus had only an inference from Radomski 1989's solution pH
(`radomski1989_extraction.md` sec. 1 answer 3: 5.2-6.1) and the formaldehyde ring's 6.24 (Kallen 1971).

| model need | what Gajda 1990 gives | tag |
|---|---|---|
| TTCA ring-N pKa | **5.32 ± 0.02** (L3, the D-xylose ring; 25.0 °C, I 0.1 M NaClO4, concentration scale) | [M] (identity [Z], Flag 1) |
| TTCA carboxyl pKa | log β012 − log β011 = 6.78 − 5.32 = **1.46** (authors: "log K_CO2H ≈ 1.5" for all) | [M] / [Z] |
| other pentose rings, same method | D-lyxose 5.10, D-arabinose 5.37, D-ribose 5.43, **L-arabinose 5.50** (= Karpinski's AThz) | [M] |
| an aliphatic, non-sugar ring | 2-propyl-TCA (L11) **6.25 ± 0.02** (= Karpinski's PrThz) | [M] |
| T-dependence of the ring-N pKa | ΔH° of ring-N protonation −20.0 to −32.3 kJ/mol (five ligands; **not measured for L3**); group II (L3's group) member L5 −25.2 ± 0.5 | [M] |
| open-chain / free-aldose fraction | none seen by 1H / 13C NMR at room temperature over the titration range (pH 1.5-8) | [Q] |
| C-2 epimer ratio | recrystallised ligands are 1 : 1 2(S) : 2(R) by NMR; the ratio "depends on the pH" (ref 16, unpublished) | [Q] / [C] |
| formation K, release / opening rates, any T ≠ 25 °C | **none** | — |

## 2. Conditions

| item | value |
|---|---|
| Ligands | L1-L10 from L-cysteine·HCl + aldose (pentoses: equimolar, water 20-40 cm³, pyridine, ~12 h, ethanol precipitation, recrystallised from water-ethanol; hexoses: methanol, reflux 1-2 h); L11 2-propyl-TCA; L12 thiaproline (literature values only). Purity by elemental analysis and acid-base titration (Table 1) |
| Titration | potentiometric, computer-controlled, **25.0 ± 0.1 °C**, N2 atmosphere, **I = 0.1 mol dm⁻³ NaClO4**, Radiometer G222B glass electrode vs Ag/AgCl; e.m.f. converted to **[H⁺] (concentration scale)** by their eqn. (2) with junction corrections; pKw 13.75; **six independent titrations** per ligand; range pH 1.5-8; ligand concentration not printed for the protonation runs (zinc runs: ~1.5 × 10⁻³ M ligand, Fig. 3 caption, exponent lost in the text layer) |
| Evaluation | PSEQUAD; overall constants β_pqr for pM + qL + rH ⇌ MpLqHr; β011, β012 are the cumulative protonation constants |
| Calorimetry | LKB-2107 batch calorimeter, three runs per ligand, NaOH dilution heat measured, water ionisation heat −54.9 kJ/mol (0.1 M KCl) |
| NMR | Varian VFT-100, room temperature, D2O, acetonitrile internal standard |
| Buffer | none (titration) |

## 3. Tables re-typed

### Ligand key [Z] (from Fig. 1, p. 3156; Flag 1)

Fig. 1 draws R (the C-2 substituent) as a Fischer chain from C-1' (attached to ring C-2) downward; C-1' = sugar C-2.

| ligand | OH sides C-1' / C-2' / C-3' (C-4') | chain end | sugar | Table 1 M (found / calc) |
|---|---|---|---|---|
| L1 | L / L / R | CH2OH | **D-lyxose** | 253.0 / 253.3 |
| L2 | L / L / R / R | CH2OH | D-mannose | 302.0 / 301.3 (monohydrate) |
| **L3** | **R / L / R** | CH2OH | **D-xylose (TTCA)** | 250.7 / 253.3 |
| L4 | R / L / R / R | CH2OH | D-glucose | 280.6 / 283.3 |
| L5 | L / R / R | CH2OH | D-arabinose | 253.4 / 253.3 |
| L6 | R / R / R | CH2OH | D-ribose | 260.2 / 263.3 [!] |
| L7 | R / L / L | CH2OH | **L-arabinose** | 260.8 / 263.3 [!] |
| L8 | R / L / L / R | CH2OH | D-galactose | 303.2 / 301.3 (monohydrate) |
| L9 | R / R / L / L | CH3 | L-rhamnose | 267.0 / 268.3 |
| L10 | (C-1' = CH2) L / R / R | CH2OH | 2-deoxy-D-glucose | 267.9 / 268.3 |
| L11 | — | propyl | butanal | 173.5 / 175.3 |
| L12 | — | H | formaldehyde (thiaproline, literature only) | — |

Check: the authors pair "L1 and L2, L3 and L4, L7 and L8" as sharing the C-1'-C-3' configuration (p. 3157); the key
gives lyxose/mannose, xylose/glucose, L-arabinose/galactose, which do share sugar C-2 to C-4. ✔

### Table 2 (p. 3157). Protonation constants (log values), 25 °C, I 0.1 M NaClO4 [M]

β011: ring-N (amino) protonation; β012: cumulative, so log K(COOH) = log β012 − log β011 [Z].

| group | ligand | log β011 (this work) | log β012 (this work) | log K(COOH) [Z] | literature β011 / β012 (ref) |
|---|---|---|---|---|---|
| I | L1 (D-lyxose) | 5.10 ± 0.02 | 6.56 ± 0.05 | 1.46 | |
| I | L2 (D-mannose) | 5.14 ± 0.02 | 6.64 ± 0.06 | 1.50 | 5.23 / 7.38 (6) |
| II | **L3 (D-xylose)** | **5.32 ± 0.02** | **6.78 ± 0.05** | **1.46** | **5.43 / 7.66 (6)** |
| II | L4 (D-glucose) | 5.31 ± 0.02 | 6.77 ± 0.06 | 1.46 | 5.37 / 7.52 (6) |
| II | L5 (D-arabinose) | 5.37 ± 0.02 | 6.84 ± 0.05 | 1.47 | |
| II | L6 (D-ribose) | 5.43 ± 0.03 | 6.94 ± 0.05 | 1.51 | 5.51 / 7.57 (6) |
| III | L7 (L-arabinose) | 5.50 ± 0.02 | 6.93 ± 0.06 | 1.43 | 5.66 / 7.92 (6) |
| III | L8 (D-galactose) | 5.53 ± 0.02 | 6.96 ± 0.05 | 1.43 | 5.67 / 7.85 (6) |
| III | L9 (L-rhamnose) | 5.51 ± 0.02 | 6.99 ± 0.06 | 1.48 | |
| IV | L10 (2-deoxy-D-glucose) | 5.53 ± 0.02 | 7.08 ± 0.06 | 1.55 | |
| — | L11 (2-propyl) | 6.25 ± 0.02 | 7.84 ± 0.05 | 1.59 | |
| — | L12 (thiaproline) | — | — | | 6.109 / 7.616 (13, 298 K); 6.104 / 7.829 (10, 310 K) |
| — | cysteine | — | — | | 10.29 / 18.445 / 20.36 (22b) |
| — | glycine | — | — | | 9.56 / 11.92 (22b) |

Ref 6 = Weitzel et al. 1959 (pH 2-6 only); the authors attribute the β012 differences to "acidity errors" in it [Q].
Errors are as printed (the authors do not say whether ± is a standard deviation; six titrations).

### Table 3 (p. 3158). Thermodynamics of the first (ring-N) protonation, 25 °C [M]

| ligand (group) | log β011 | ΔG° / kJ mol⁻¹ | ΔH° / kJ mol⁻¹ | ΔS° / J K⁻¹ mol⁻¹ |
|---|---|---|---|---|
| L1 (I) | 5.10 | −29.14 | −20.0 ± 0.4 | 31.1 ± 1.2 |
| L5 (II) | 5.37 | −30.65 | −25.2 ± 0.5 | 18.7 ± 1.6 |
| L7 (III) | 5.50 | −31.41 | −27.8 ± 0.4 | 12.4 ± 1.4 |
| L9 (III) | 5.51 | −31.49 | −27.5 ± 0.5 | 13.6 ± 1.6 |
| L10 (IV) | 5.53 | −31.53 | −32.3 ± 0.5 | −2.3 ± 2.0 |

L3 (TTCA) was not measured calorimetrically. The authors explain the spread by an intramolecular N-H···O(C-3') hydrogen
bond that opens on protonation [P].

### Table 4 (p. 3159). Zinc(II) formation constants — not transcribed (no zinc in any pot). L3's log β110 = 4.05 ± 0.05,
log β120 = 6.86 ± 0.11 [M], for completeness.

### NMR (p. 3156) [Q]

The acyclic form would show a CH=N 1H signal and hydrolysis an aldose carbonyl 13C signal; "the latter signals did not
appear" at any pH of the study, with or without zinc. Sugar rings are "stable and do not open under mild conditions"
[Q, citing ref 6]. Detection limit not stated.

## 4. Arithmetic (mine) [Z]

**1. TTCA ring-N pKa at process temperature (van 't Hoff, ΔCp = 0).** pKa(T) = 5.32 + (ΔH/(R ln 10))(1/T − 1/298.15),
with ΔH the protonation enthalpy taken from the measured spread (L3's own not measured):

| ΔH (kJ/mol) | 90 °C | 100 °C | 120 °C |
|---|---|---|---|
| −20.0 (L1, group I) | 4.69 | 4.62 | 4.47 |
| **−25.2 (L5, L3's group II)** | **4.53** | **4.43** | **4.25** |
| −27.8 (L7) | 4.45 | 4.34 | 4.14 |
| −32.3 (L10) | 4.31 | 4.18 | 3.95 |

Caveat: the two literature values for thiaproline (L12) at 298 and 310 K differ by only 0.005 in log β011, which taken
at face value would imply ΔH ≈ 0. They come from two laboratories and media (refs 10, 13), so they are not a test of
ΔH, but they are a reason to carry ΔH = 0 (pKa 5.32 at every T) as the other end of a sensitivity.

**2. The Kallen-shaped K ratio R-K needs.** With the case's own form K_app ∝ (1 + 10^(pH − pKa)) (cysteine
zwitterion fraction ≈ 1 below pH 7), K(6.5) / K(4.4) = 18.1 / 14.4 / 10.2 / 3.4 at pKa 5.2 / **5.32** / 5.50 / 6.1
(the case's 3.4 - 18x reproduced at its band ends). At the in-situ pKa of arithmetic 1 (4.53 / 4.43 / 4.25 / 3.95):
54 / 61 / 74 / 93x. **For ANY ring-N pKa the ratio cannot exceed 10^(6.5 − 4.4) = 126x.** So a Kallen-shaped K
driven by the ring-N pKa alone cannot give the ~300x the case says (c) needs, whatever Gajda's number is.

**3. Against Karpinski 1990 (same rings).** AThz (L7, pKa 5.50): protonated-ring fraction 0.24 at pH 6.0 and 0.012 at
pH 7.4, a ratio of 19; Karpinski's release constants fall only 3.8x (1.8 × 10⁻⁴ → 4.8 × 10⁻⁵ min⁻¹). Release at
25 °C is therefore NOT proportional to the protonated ring; a neutral-ring term carries most of it at pH 7.4.

**4. Against Kallen 1971.** TC ring-N pKa 6.24 (I 1.0 M KCl) vs thiaproline 6.10 (literature here, other media):
consistent within medium effects. The sugar chain lowers the ring-N pKa by 0.7-1.1 units relative to thiaproline,
and by 0.7-1.2 relative to the 2-propyl ring (6.25), which the authors attribute to the inductive effect of the polyol
chain [Q].

## 5. What it can and cannot test in the model

- **Can:** fix TTCA's ring-N pKa at 25 °C (5.32 ± 0.02), on a concentration pH scale at I 0.1 M; give its
  temperature shift through a same-class ΔH band (−20 to −32 kJ/mol); give the carboxyl pKa (~1.46) for charge
  balance (`ph_state`); support, qualitatively, a small open-chain fraction at room temperature (P1a).
- **Cannot:** give K for Cys + Xyl ⇌ TTCA, or any rate; test anything at process temperature directly; separate
  in-situ from cooled-scale pH effects; size P1a or P1b numerically (no NMR detection limit).
- **Scale caveat:** a concentration-scale pKa at I 0.1 M is ~0.1 unit above the activity-scale value a pH meter in a
  0.2 M phosphate pot would correspond to (γH ≈ 0.8) [Z, approximate]; smaller than the T effect.

## 6. Conflicts with existing records

| record | what it says | what Gajda says | verdict |
|---|---|---|---|
| `radomski1989_extraction.md` (pKa inferred 5.2-6.1 from solution pH 3.7-4.0) | inference, no titration | 5.32 ± 0.02, titrated | **consistent**; Gajda replaces the inference (it sits in the lower part of the band) |
| `radomski1989_extraction.md` (TTCA 2R : 2S = 4.06 at equilibrium, 37 °C) | epimers 4 : 1 | recrystallised ligands 1 : 1, ratio pH-dependent | **apparent conflict**; probably freshly dissolved crystals vs equilibrated solution; Gajda gives no equilibrium ratio. Radomski's 4.06 stands for P2 |
| `kallen1971_extraction.md` (TC pKa 6.24, I 1.0 M) | formaldehyde ring | thiaproline 6.10 (lit.), sugar rings 5.1-5.5 | consistent; the sugar rings are ~0.9 units more acidic, so Kallen's 6.24 is not TTCA's |
| `karpinski1990_extraction.md` (AThz, PrThz release vs pH) | no ring pKa | AThz 5.50, PrThz 6.25 | adds the missing pKa (arithmetic 3) |
| `pesek1975_extraction.md`, `nagasawa1981_extraction.md` | sugar/aliphatic rings stable in neutral water | no open chain by NMR at pH 1.5-8 | consistent |
| Weitzel 1959 (ref 6, not on disk) | L3 5.43 / 7.66 | 5.32 / 6.78 | β011 within 0.11; the β012 difference is attributed to acidity errors in ref 6 [Q] |

## 7. Flags

1. **[!] The ligand key is mine.** Fig. 1 shows structures only; L3 = D-xylose is read from the Fischer projection
   (C-1' OH right, C-2' left, C-3' right, CH2OH) and confirmed by the authors' pairing L3/L4 (xylose/glucose share
   C-2 to C-4) and by Table 1 (pentose analysis). P(L3 is the D-xylose ring) ~0.95.
2. **[!] Table 1 calc. M for L6 and L7 (263.3) does not match a pentose ring (253.3)** although their C/H/N match it;
   a typesetting slip in the M column (the found values 260 sit between). L6/L7 are pentose rings by Fig. 1.
3. **Concentration scale, I 0.1 M NaClO4, 25 °C only**; process pots are 0.2 M phosphate or unbuffered water at
   90-140 °C.
4. **L3's ΔH is not measured**; the T-shift in arithmetic 1 borrows group II's L5.
5. The 1 : 1 epimer statement is for recrystallised material (no equilibration time stated).
6. The paper's title in print says "Proton and Zinc(II) Complexes ..."; the brief's "protonation constants" is the
   content, not the title.

## M1c role (frozen bd641ae8)

**Role fixed by the case (sec. 4.3):** band source for TTCA's ring-N pKa (R-K). Sec. 5.2 names this paper as the
replacement for Radomski's inferred pKa_T band, by amendment.

**Band it can supply:**
- **Coordinate:** R-K's pKa_T.
- **Numbers:** 5.32 ± 0.02 at 25.0 °C (concentration scale, I 0.1 M NaClO4). At process temperature, by the
  same-class ΔH band −20 to −32 kJ/mol (ΔCp = 0): 4.31-4.69 at 90 °C, 4.18-4.62 at 100 °C, 3.95-4.47 at 120 °C;
  ΔH = 0 (no shift) as the opposite sensitivity.
- **Transfer:** direct for the compound (it IS the D-xylose ring) at 25 °C; the T-shift is a same-class transfer.
  On the cooled scale (case open question 1) the 25 °C value is the relevant one; in situ, the shifted one.
- **What it does to R-K (arithmetic, no engine):** K(6.5)/K(4.4) = 14.4x at 25 °C and 54-93x in situ, against
  the ~300x needed; the ceiling for any ring-N pKa is 126x. R-K's arithmetic failure on (c) stands under this band.
- **Not decisive:** 25 °C, no K, no rate. It measures nothing the case would call decisive (sec. 4.3 requires K or
  a release rate at ≥ 90 °C).
- **Secondary (qualitative):** P1a f_oc ~0 at room temperature by NMR (consistent with the 0-0.05 band, no
  detection limit, so it cannot narrow it).
