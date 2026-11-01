# Scalone, Cucu, De Kimpe & De Meulenaer 2015 — EXTRACTION
### *"Influence of Free Amino Acids, Oligopeptides and Polypeptides on the Formation of Pyrazines in Maillard Model Systems"*
### Tryptic whey-protein hydrolysate (or native whey, or whey + the hydrolysate's free amino acids, or free lysine) with glucose; "roasting" = freeze-dried on sand, 180 °C / 90 min; "sterilization" = 0.1 M K-phosphate pH 7.8 solution, 120 °C (inside) / 2 h; a water-activity series 0.16-0.75 (roasting). HS-SPME-GC-MS **peak areas only**; 2,5- and 2,6-DMP co-elute; **unsubstituted pyrazine is not reported in any table.**

**Read 2026-10-31 from `data/articles/scalone2015.pdf`** (26 pp.: the ACS "Just Accepted" manuscript, pp. 1-25 plus a
cover page; the version of record is JAFC 63:5364-5372). All numbers from `pdftotext -layout` of Tables 1-5 (pp.
21-24); every row checked for column count (8 / 9 / 6 / 6 columns). Figs. 1-3 (schemes, a gel-permeation trace) not
used for numbers. Derived numbers by hand.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value |
|---|---|
| Authors | Gustavo Luis Leonardo Scalone, Tatiana Cucu, Norbert De Kimpe, Bruno De Meulenaer (Ghent University; nutriFOODchem and Sustainable Organic Chemistry) |
| Venue | *J. Agric. Food Chem.* **2015**, 63, 5364-5372 (manuscript posted 14 May 2015; the PDF on disk is the accepted manuscript, not the typeset article [!]) |
| DOI | **10.1021/acs.jafc.5b01129** (printed; Crossref 2026-10-31: same title, vol. 63, pp. 5364-5372, 2015) |
| Lab | Ghent (De Kimpe co-author): same group as Adams 2008 and Van Lancker 2010 / 2012 |

## 1. What this paper answers for the model

| need | what Scalone 2015 gives | tag |
|---|---|---|
| pyrazine ratios for R | **No.** Pyrazine (99 %) was bought as a standard but PZ appears in no table; R cannot be computed | [M] / [!] |
| MPZ / 2,5(6)-DMP (areas) | roasting, hydrolysate 15 mg + glucose 7.5 mg: 23.10 / 42.39 (0.54); a_w 0.33: 54.45 / 419.99 (0.13); free Lys roasting 0.11-0.16; free Lys sterilization 0.76-1.6; hydrolysate sterilization 0.14-0.65 [Z] | [M] / [Z] |
| response factors | absolute peak area of each pyrazine, no IS, no calibration; dry vs wet headspace partition differs **16x** for a 2 µg 2-ethyl-3-methylpyrazine spike (107.89 vs 6.68 x 10⁶) [M]: roasting and sterilization areas are not comparable | [M] |
| 13C / 15N | none | — |
| condensation k, pH, Ea | none (two temperature-time-water combinations, not a T series) | — |
| pyrrole vs pyrazine | none | — |
| **peptide vs free amino acid** | hydrolysate > native whey ≈ native whey + the hydrolysate's free amino acids (Table 4); free amino acids alone (at the hydrolysate's level) + glucose: no pyrazines; free lysine (14.61 mg, 0.1 mmol) gives ~10x the 2,5(6)-DMP area of 15 mg hydrolysate | [M] |
| a_w | roasting pyrazines peak at a_w 0.33 and fall to 0.75 (Table 5); 2,5(6)-DMP 42 → 420 → 212; 3-ethyl-2,5-DMP 5.1 → 263 → 163 | [M] |

**Answers.**

1. **No R, so no test of the mixing bound; does not distinguish P1x / P1y from P2b.**
2. **Peptides vs free amino acids, per amine** [Z]: the authors conclude the role of peptides is "dramatically
   underestimated" [Q]. Per mole of alpha-amine the picture is less dramatic: a tryptic digest of 15 mg whey protein
   carries very roughly 10-15 µmol of new N-termini (one cleavage per ~1.1 kDa; my estimate, not measured), and gives
   2,5(6)-DMP 42 area units; 0.1 mmol free lysine gives ~430. Per alpha-amine that is ~3-4 vs ~4 — comparable.
   What the paper does show cleanly is that the hydrolysate's FREE amino acids (~1.5 % of its mass, Table 1) cannot
   explain its pyrazines (Table 4). Confidence in the per-amine comparison: low (~0.4), since the N-terminus count
   is my estimate and lysine's epsilon-amine adds to both.
3. **The a_w series is confounded** [!]: the a_w 0.16 column of Table 5 is the Table 2 roasting column for 7.5 mg
   glucose (identical values: 23.10 ± 0.36, 42.39 ± 3.04, 11.00 ± 0.51, ...), i.e. freeze-dried and capped at once;
   the a_w 0.33-0.75 samples were held 7 days at 30 °C over salt solutions before roasting. The 10-50x jump from
   0.16 to 0.33 therefore mixes a_w with a 7-day 30 °C pre-incubation (early Maillard / Amadori formation). The
   0.33 → 0.75 decline (all pre-incubated) is the cleaner a_w effect: 2,5(6)-DMP −0.3 dex.

## 2. Conditions [M]

| item | value |
|---|---|
| Protein | whey protein isolate LACPRODAN DI-9224 (84 % protein), 5 mg/mL in 0.1 M K-phosphate pH 7.8, 95 °C / 5 min; trypsin 1:20 w/w, 37 °C, 20 h; 95 °C / 5 min to inactivate. Native WPI: no detectable free amino acids |
| Free amino acids in the hydrolysate (Table 1, mg/100 mg) | Lys 0.9916, Leu 0.1340, Val 0.0960, Ile 0.0728, Arg 0.0546, Tyr 0.0421, Thr 0.0389, Trp 0.0359, Phe 0.0297, Met 0.0218, Ala 0.0104 |
| Sugar | D-glucose 1.25-20 mg per 15 mg hydrolysate (Table 2; protein : glucose 1:0.08-1:1.33); 2.5-37 mg per 14.61 mg lysine (Table 3; text says 1:2.53 to 1:0.15 and "8 different ratios", the table has 9 columns and 2.5 / 14.61 = 0.17) [!] |
| Buffer / pH | all samples made up to 5 mL in 0.1 M K-phosphate pH 7.8 (roasting samples then freeze-dried with 1.5 g sand); pH not re-measured [!] |
| "Sterilization" | 5 mL solution in a capped SPME vial, oil bath 130 °C to reach 120 °C inside, 2 h, ice |
| "Roasting" | freeze-dried on 1.5 g sand, capped, oven 180 °C, 90 min |
| a_w | freeze-dried a_w 0.16 (dew-point meter); equilibrated 7 days at 30 °C over MgCl2 (32 % RH), K2CO3 (43), NaBr (56), KCl (67), NaCl (75); measured a_w 0.33 / 0.38 / 0.53 / 0.62 / 0.75 |
| Analysis | HS-SPME DVB/CAR/PDMS, 30 min at 35 °C; DB-5; EI scan; identification by LRI + MS (α: authentic standard; β: tentative); triplicates (3 independent determinations), Games-Howell letters |

## 3. Quantitative results [M] (peak area x 10⁶, mean ± SD of 3)

### 3.1 Table 2 (hydrolysate 15 mg + glucose), MPZ and 2,5(6)-DMP

| glucose (mg) | 1.25 | 2.5 | 5 | 7.5 | 10 | 15 | 17.5 | 20 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| MPZ roasting | 4.15 | 9.46 | 16.21 | 23.10 | 22.54 | 20.97 | 18.86 | 20.87 |
| MPZ sterilization | nd | nd | 0.94 | 1.75 | 3.46 | 5.45 | 5.35 | 6.32 |
| 2,5(6)-DMP roasting | 24.45 | 35.22 | 37.81 | 42.39 | 35.25 | 25.41 | 22.75 | 20.22 |
| 2,5(6)-DMP sterilization | nd | 5.65 | 6.75 | 9.64 | 7.49 | 11.16 | 8.32 | 9.66 |

Also roasting only: 2,3-DMP 1.13-2.54; 2-ethylpyrazine 3.45-4.96; 2-ethyl-5-methyl 3.79-11.00; 3-ethyl-2,5-DMP
2.42-6.78; 2-ethyl-3,5-DMP 1.26-2.28 (β); 5-ethyl-2,3-DMP 1.03-1.23 (β).

### 3.2 Table 3 (lysine 14.61 mg + glucose), MPZ and 2,5(6)-DMP

| glucose (mg) | 2.5 | 4.5 | 9 | 13.5 | 18 | 22.5 | 27 | 31.5 | 37 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| MPZ roasting | 49.05 | 36.80 | 35.96 | 42.74 | 48.49 | 45.16 | 40.69 | 43.36 | 40.57 |
| MPZ sterilization | nd | 6.75 | 13.28 | 15.11 | 16.97 | 25.45 | 30.33 | 34.82 | 38.87 |
| 2,5(6)-DMP roasting | 430.58 | 426.55 | 402.27 | 450.48 | 448.55 | 433.84 | 392.45 | 293.06 | 251.48 |
| 2,5(6)-DMP sterilization | nd | 8.92 | 19.07 | 25.24 | 20.95 | 27.52 | 29.69 | 25.57 | 24.37 |

Roasting also: 2-ethyl-5-methyl up to 278.17, 3-ethyl-2,5-DMP up to 135.47, 5-ethyl-2,3-DMP (β) up to 222.19,
2,5-diethyl (β) up to 94.04, among 13 pyrazines.

### 3.3 Table 4 (roasting): hydrolysate vs native whey vs whey + free AA (glucose 7.5 / 15 mg)

MPZ: 23.10 / 20.97 vs 8.99 / 8.53 vs 9.66 / 8.06; 2,5(6)-DMP: 42.39 / 25.41 vs 9.02 / 9.29 vs 10.99 / 9.37. 2,3-DMP,
ethyl-, 3-ethyl-2,5-, 2-ethyl-3,5-, 5-ethyl-2,3-DMP only with the hydrolysate. Free AA (Table 1 amounts) + glucose:
"no pyrazines were detected (data not shown)".

### 3.4 Table 5 (roasting; hydrolysate 15 mg + glucose 7.5 mg) vs a_w

| a_w | 0.16 | 0.33 | 0.38 | 0.53 | 0.62 | 0.75 |
|---|---:|---:|---:|---:|---:|---:|
| MPZ | 23.10 | 54.45 | 37.99 | 30.40 | 22.39 | 19.71 |
| 2,5(6)-DMP | 42.39 | 419.99 | 343.52 | 327.52 | 287.13 | 211.5 |
| 3-ethyl-2,5-DMP | 5.14 | 262.84 | 249.39 | 247.43 | 254.57 | 163.00 |
| 2-ethyl-5-methyl | 11.00 | 122.23 | 101.74 | 97.00 | 107.68 | 85.2 |

## 4. Mechanism evidence and its kind

- **Within-study factors** (kind 6), Ghent: amine form (hydrolysate / native / native + free AA / free Lys), sugar
  ratio, a_w, dry vs wet heating. Peak areas only, triplicates.
- The peptide aminoketone mechanism (Fig. 2) is Van Lancker 2010's Scheme 1 [C] / [P]; nothing here tests it.
- The fall of pyrazines at high glucose is attributed to acid formation lowering pH [P] (pH not measured).

## 5. Bearing on open questions

| item | effect | why |
|---|---|---|
| EX-052 mixing test / CF-077 | LEAVES UNCHANGED | no PZ |
| P1x / P1y vs P2b | **does not distinguish** | |
| evidence map 1.1 (peptides, protein-bound amine) | SUPPORTS (direction): a protein hydrolysate makes pyrazines that its free amino acids cannot | the engine's Strecker step needs a free alpha-amino acid; a hydrolysate's peptides are outside it |
| evidence map 1.6 (water activity) | adds a confounded series | the 0.16 point lacks the 7-day pre-incubation |
| triage 9 (protein-bound amine on the sensitivity page) | SUPPORTS keeping it a band | per-amine reactivity of peptides vs Lys is not resolvable from areas |

## 6. Disagreements found

1. **Self-reuse within the paper** [!]: Table 5 a_w 0.16 = Table 2 roasting at 7.5 mg glucose (all eight shared rows
   identical) = Table 4 hydrolysate / 7.5 mg, except Table 4 prints 2-ethylpyrazine 4.23 (vs 4.24), 2-ethyl-5-methyl
   10.70 ± 0.11 (vs 11.00 ± 0.51) and 5-ethyl-2,3-DMP 1.03 ± 0.01 (vs 1.06 ± 0.06). Same for Table 4 / 15 mg vs Table
   2 / 15 mg (2,3-DMP 2.17 vs 2.16; 2-ethyl-5-methyl 8.81 vs 8.80). One data set printed three times with small
   transcription drift; it is one measurement, not three.
2. **Lysine-ratio text vs Table 3**: "8 different ratios ... 1:2.53 to 1:0.15" vs 9 columns starting at 2.5 / 14.61 =
   1:0.17. Cosmetic.
3. **With Van Lancker 2010 (same lab)**: there, peptides beat free amino acids per mole of peptide (1 mmol dipeptide
   vs 1 mmol each free amino acid); here, free lysine beats the hydrolysate per mass by ~10x. Not a contradiction:
   different basis (mole of peptide vs mass of hydrolysate), dry vs aqueous, 180 vs 130 °C.

## 7. Ledger-ready measurements

1. Table 4 (hydrolysate vs native vs native + free AA, 2 glucose levels): a within-study peptide-vs-free-amino-acid
   factor for whey (peak_area_only).
2. Table 5 a_w 0.33-0.75 (pre-incubated arm only): an a_w factor for roasting pyrazines; exclude the 0.16 point or key
   it with the confound.
3. Table 1: the free amino acid composition of a tryptic whey hydrolysate (mg/100 mg).
