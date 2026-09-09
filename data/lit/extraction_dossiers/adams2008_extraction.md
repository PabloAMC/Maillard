# Adams 2008 — EXTRACTION (twenty amino acids + 1,3-dihydroxyacetone neat at 90 C / 30 min, SPME peak areas; alanine + methylglyoxal ~0.85 mol/L in buffer, calibrated yields of 2,5(6)-dimethylpyrazine and 2,5-diacetyl-3-methylpyrrole)
### The amino-acid screen for the aldehyde-addition step: 2,5(6)-dimethylpyrazine from every amino acid, the extra substituent from the Strecker aldehyde (Gly -> trimethyl, Ala -> ethyldimethyl, Val -> isobutyldimethyl) or from sugar-fragment aldehydes; and, in water with excess MGO, a pyrrole outcompetes the pyrazine.

**Source on disk:** `data/articles/adams2008.pdf` (owner's download, 2026-09-08). Read-only extraction from
the pypdf text layer in the scratchpad; Table 1's sparse rows lose their columns in the text layer, so
every Table 1 cell below was assigned from `pdftotext -bbox` word coordinates against the amino-acid
header positions (all offsets within 5 pt of a 22 pt column pitch; the assignments agree with every
product the Results text attributes to an amino acid). Tables 2 and 3 are clean. Schemes 1-3 are
mechanisms, not data. Repo status before this dossier: roadmap §5b lists the ethyl- and
trimethylpyrazines with "the aldehyde-addition step; amino-acid identity" missing.

## 0. Identity

| field | value |
|---|---|
| Title | "Formation of Pyrazines and a Novel Pyrrole in Maillard Model Systems of 1,3-Dihydroxyacetone and 2-Oxopropanal" |
| Authors | An Adams, Viviana Polizzi, Norbert De Kimpe (Ghent); Martinus van Boekel (Wageningen) |
| Venue | J. Agric. Food Chem. 2008, 56, 2147-2153 |
| DOI | 10.1021/jf0726785 |
| Naming | 2-oxopropanal = methylglyoxal (MGO); DHA = 1,3-dihydroxyacetone (the dimer, 98 %), used as an in-situ MGO source by dehydration; 2,5(6)-dimethylpyrazine = the pair unresolved on HP5-MS; **"3-ethyl-2,5-dimethylpyrazine" (LRI 1082) is the same molecule as 2-ethyl-3,6-dimethylpyrazine (Cerny 1994's isomer I); "2-ethyl-3,5-dimethylpyrazine" (LRI 1088) is Cerny's II**; compound 34 = 2,5-diacetyl-3-methyl-1H-pyrrole; 35 = 2,5(6)-dimethylpyrazine; 21 = 1-amino-2-propanone (aminoacetone) |
| Companions | Adams 2004 (proline + DHA -> ATHP, ref 15); Van Lancker 2012 (`vanlancker2012_extraction.md`, same laboratory, peptides); Shu 1999 (`shu1999_extraction.md`, ref 12); Amrani-Hemaimi 1995 (ref 19, 13C-alanine: the ethyl of 3-ethyl-2,5-dimethylpyrazine is 100 % from alanine) |

## 1. Why it matters

B18 gives the trunk pyrazine, methylpyrazine and 2,5-dimethylpyrazine from a single aminoketone per
dicarbonyl. The next step is the extra substituent. This paper runs the one dicarbonyl the trunk has
most of (MGO, supplied in situ from DHA) against all twenty protein amino acids and reads off which
extra pyrazines each makes: the amino acid's Strecker aldehyde is added to the 3,6-dimethyl-2,5-
dihydropyrazine (Scheme 1A) — formaldehyde from glycine gives trimethylpyrazine, acetaldehyde from
alanine gives 3-ethyl-2,5-dimethylpyrazine, 2-methylpropanal from valine gives the isobutyl analogue —
but the same adducts also appear from arginine, lysine and serine, which make no such aldehyde, so the
aldehyde pool is partly the sugar's own. In water at ~1 mol/L with MGO in excess the main product is
not a pyrazine at all but 2,5-diacetyl-3-methylpyrrole (one aminoacetone + two MGO), and the paper's
Table 3 gives calibrated yields of both against pH, buffer, time, temperature and the Ala : MGO ratio —
the only calibrated 2,5(6)-dimethylpyrazine yields from fed MGO + amino acid in the corpus besides
Zhou 2024, at forty times Zhou's concentration.

## 2. Methods as they matter to a model

- **Screen (Table 1):** 5 mmol amino acid + 5 mmol DHA ground together, **no solvent**, 20 mL headspace
  vial, **90 C oil bath, 30 min**, ice. Headspace SPME 30 min at 30 C (50/30 µm DVB/CAR/PDMS), 2 min
  desorption at 250 C; HP5-MS 30 m x 0.25 mm x 0.25 µm; 35 C (5 min) -> 80 C at 2 C/min -> 250 C at
  20 C/min. **Quantity: GC-MS peak area x 1e-6, no internal standard** -> peak_area_only.
  Identification by MS + LRI vs Wagner 1999 (ref 4); 2-ethyl-6-methyl- and the acyl-pyrazines
  tentative (footnote a).
- **Aqueous pots (Tables 2, 3):** 10 mmol amino acid in 10 mL **1 M phosphate pH 7.0** + equimolar
  dicarbonyl (**MGO as 1800 µL of 40 % aqueous**, ~10.5 mmol), pH re-set to 7 with 1 N NaOH, **100 C,
  30 min**, oil bath, ice. Nominal "1 M"; with the added MGO solution the volume is ~11.8 mL so
  **~850 mmol/L each**. Table 3 varies: 130 C; NaH2PO4/Na2HPO4 or NaHCO3/H2CO3 at pH 7; acetate at pH 4;
  10 / 30 / 60 min; Ala : MGO 1:1, 2:1, 1:3 (the 1:3 arm adds ~5.4 mL MGO solution, so its
  concentrations are more dilute still). Table 2: seven amino acids at 1 M (0.5 M for Asp and Phe,
  solubility). Text says the alanine pot was "1 M, 90 C, 30 min, pH 7" where Methods say 100 C (flag 4).
- **Work-up:** pH to 9.0 (2 N NaOH), CHCl3 3 x 10 mL, MgSO4, concentrated to ~1/10, **internal standard
  2-acetyl-1-methylpyrrole** (150 µL of 1 % in CH2Cl2), GC-MS with **response factors calculated** ->
  calibrated. Yields printed as "%" **with respect to the amino acid**; the basis (molar) is not stated
  but an internal-standard yield is molar by construction; read as **mol %** (flag 3). Table 3 carries
  a +/- (replicate spread; n not stated).
- **Pyrrole 34 identity:** isolated by silica chromatography, 1H/13C NMR, MS (m/z 150 base, 165 M+),
  LRI 1495; confirmed by synthesis from 3-methylpyrrole (78 % pure reference). Intermediate
  2,5-diacetyl-3-hydroxy-4-methylpyrrole (32) tentatively identified (MS, LRI 1644).
- **Reporting units summary:** Table 1 peak areas (x 1e6); Tables 2-3 % of amino acid.

## 3. Tables re-typed

### Table 1. "Pyrazines (GC-MS Peak Area x 1e-6) Detected in the Headspace of Model Reactions of 1,3-Dihydroxyacetone with Various Amino Acids (90 C, 30 min)"

Column order in the paper: Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val.
Cells not listed = not detected. LRI on DB-5 / HP5-MS.

| compound | LRI | cells (amino acid: area x 1e-6) |
|---|---:|---|
| methylpyrazine | 818 | Ala 1.20; Arg 14.18; Asn 0.60; Gly 3.35; Lys 10.34; Ser 38.48 |
| 2,5(6)-dimethylpyrazine | 912 | **all twenty** — see the full row below |
| 2-ethyl-6-methylpyrazine (tent.) | 1000 | Ala 5.70; Gly 11.58; Ser 6.14 |
| trimethylpyrazine | 1003 | Arg 148.63; Gly 18.65; Thr 1.55 |
| 3-ethyl-2,5-dimethylpyrazine (= 2-ethyl-3,6-dimethyl-) | 1082 | Ala 3.63; Arg 44.97; Cys 0.40; Gly 3.55; Lys 33.24; Ser 7.67 |
| 2-ethyl-3,5-dimethylpyrazine | 1088 | Ala 5.29 |
| 2-methyl-propylpyrazine (b) | 1092 | Ala 5.82 |
| 3-ethenyl-2,5-dimethylpyrazine | 1103 | Arg 11.94 |
| 2-acetyl-5-methylpyrazine (tent.) | 1126 | Gly 1.46 |
| 2-acetyl-6-methylpyrazine (tent.) | 1132 | Ala 0.29; Gln 0.60; Gly 0.46; Val 2.58 |
| 2-methyl-(E-1-propenyl)pyrazine (b) | 1163 | Ala 6.89 |
| 2-acetyl-3,5-dimethylpyrazine (tent.) | 1173 | Ala 0.93; Gly 4.19 |
| 2,5-dimethyl-3-(2-methylpropyl)pyrazine | 1187 | Val 0.63 |
| 2-propanoyl-5-methylpyrazine (tent.) | 1194 | Ala 2.45; Asp 1.78 |
| 2-propanoyl-6-methylpyrazine (tent.) | 1198 | Ala 12.3; Asp 2.00; Glu 0.68; Thr 0.57 |

(b) "Elution order of the isomers is not determined in ref 4."

2,5(6)-dimethylpyrazine row (area x 1e-6):

| Ala | Arg | Asn | Asp | Cys | Gln | Glu | Gly | His | Ile | Leu | Lys | Met | Phe | Pro | Ser | Thr | Trp | Tyr | Val |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 84.8 | 493.9 | 164.3 | 15.6 | 36.5 | 29.5 | 9.5 | 330.5 | 25.4 | 16.6 | 11.2 | 385.5 | 15.3 | 27.9 | 1.6 | 145.8 | 22.1 | 18.4 | 8.0 | 11.1 |

Rank: Arg > Lys > Gly > Asn > Ser > Ala > Cys > Gln > Phe > His > Thr > Trp > Ile > Asp > Met > Leu >
Val > Glu > Tyr > Pro. Text: only the more stable Strecker aldehydes were seen in the headspace (Ile,
Leu, Met, Phe, Val); proline gives "very little pyrazines"; the largest pyrazine variety from Ala, Gly,
Arg, Ser.

### Table 2. "Yields of 2,5-Diacetyl-3-methyl-1H-pyrrole (34) from the Model Reaction of 2-Oxopropanal (24) with Various Amino Acids (1 M, 100 C, 30 min, Phosphate Buffer)"

| amino acid | yield of 34 (%) |
|---|---:|
| alanine | 0.02 |
| asparagine | 0.36 |
| aspartic acid (0.5 M, solubility) | 1.56 |
| glutamic acid | 1.81 |
| leucine | 0.16 |
| phenylalanine (0.5 M) | 0.04 |
| tryptophan | 0.02 |

Amino acids not listed gave no detectable 34. The authors attribute the Asp and Glu values to the pH
falling below 7 (the 1 M buffer could not hold it).

### Table 3. "Yields of 2,5-Diacetyl-3-methyl-1H-pyrrole (34) and 2,5(6)-Dimethylpyrazine (35) from the Model Reaction of Alanine (36) with 2-Oxopropanal (24)" (yields with respect to alanine; nominal 1 M; +/- as printed)

| time (min) | temp (C) | Ala : MGO | buffer | initial pH | yield of 34 (%) | yield of 35 = 2,5(6)-DMP (%) |
|---:|---:|---|---|---:|---:|---:|
| 30 | 130 | 1:1 | NaH2PO4/Na2HPO4 | 7 | 0.020 +/- 0.0073 | 0.02 |
| 30 | 130 | 1:1 | NaHCO3/H2CO3 | 7 | 0.055 +/- 0.0094 | 0.011 +/- 0.00064 |
| 30 | 130 | 1:1 | CH3COOH/CH3COONa | 4 | 1.6 +/- 0.040 | 0.042 +/- 0.005 |
| 10 | 100 | 1:1 | acetate | 4 | 0.080 +/- 0.039 | nd |
| 30 | 100 | 1:1 | acetate | 4 | 1.1 +/- 0.087 | 0.020 +/- 0.00025 |
| 60 | 100 | 1:1 | acetate | 4 | 1.3 +/- 0.045 | 0.026 +/- 0.0030 |
| 30 | 100 | 2:1 | acetate | 4 | 0.50 +/- 0.019 | 0.013 +/- 0.0029 |
| 30 | 100 | 1:3 | acetate | 4 | 3.7 +/- 0.074 | 0.016 +/- 0.0040 |

Text on the aqueous alanine pot (pH 7, phosphate, extract): "2,5-dimethylpyrazine, 3-ethyl-2,5-
dimethylpyrazine, and some furanones ... with an unknown compound [34] as the main reaction product".
No similar pyrrole from alanine with glyoxal, phenylglyoxal or 2,3-pentanedione.

## 4. Numbers and steps the repository can use

Registry keys: `methylpyrazine`, `2_5_dimethylpyrazine` / `2_6_dimethylpyrazine` (one unresolved number),
`trimethylpyrazine`, `2_ethyl_3_5_dimethylpyrazine` (whose stored SMILES is in fact the 3-ethyl-2,5-
isomer; `cerny1994_extraction.md` flag 8) exist. 3-Ethyl-2,5-dimethylpyrazine as its own key,
2-ethyl-6-methylpyrazine, the acyl-, ethenyl-, propenyl- and isobutyl-pyrazines, 2,5-diacetyl-3-
methylpyrrole, DHA, methylglyoxal, aminoacetone, the amino acids: **not in registry**.

| quantity or step | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| Ala + MGO -> 2,5(6)-DMP | 0.02 (phosphate) / 0.011 (bicarbonate) | mol % of Ala | 1:1, ~850 mmol/L each, pH 7, 130 C, 30 min | Table 3 rows 1-2 | fed_intermediate_yield (calibrated; the phosphate row has no +/-) |
| Ala + MGO -> 2,5(6)-DMP, pH 4 | 0.042 (130 C); nd / 0.020 / 0.026 at 10 / 30 / 60 min (100 C) | mol % of Ala | 1:1, acetate | Table 3 | fed_intermediate_yield, short time series |
| Ala + MGO -> 2,5(6)-DMP vs ratio | 0.013 (2:1) / 0.020 (1:1) / 0.016 (1:3) | mol % of Ala | 100 C, 30 min, pH 4 | Table 3 | fed_intermediate_yield (per MGO the 2:1 and 1:1 arms are equal: 0.026 vs 0.020) |
| Ala + MGO -> pyrrole 34 | 0.020 / 0.055 / 1.6 (130 C; phosphate 7 / bicarbonate 7 / acetate 4); 0.080 / 1.1 / 1.3 (100 C, pH 4, 10 / 30 / 60 min); 0.50 (2:1); 3.7 (1:3) | mol % of Ala | as Table 3 | Table 3 | fed_intermediate_yield (the competing sink for aminoacetone AND for two MGO) |
| pyrrole 34 by amino acid | Ala 0.02; Asn 0.36; Asp 1.56; Glu 1.81; Leu 0.16; Phe 0.04; Trp 0.02; others nd | mol % | 1 M (0.5 M Asp, Phe), pH 7 phosphate, 100 C, 30 min | Table 2 | fed_intermediate_yield (pH-confounded for Asp, Glu) |
| pH 4 / pH 7 (phosphate), 130 C, 30 min | 2,5(6)-DMP 2.1x; pyrrole 80x | — | same pot | Table 3 | within_study_ratio (acid helps the pyrazine a little and the pyrrole a lot — opposite in sign to B18's pH slope below 7, flag 6) |
| phosphate / bicarbonate at pH 7, 130 C | 2,5(6)-DMP 1.8x; pyrrole 0.36x | — | | Table 3 | within_study_ratio (buffer-anion catalysis of the pyrazine; cf. Rizzi 2004) |
| 130 C / 100 C at pH 4, 30 min | 2,5(6)-DMP 2.1x; pyrrole 1.5x | — | | Table 3 | within_study_ratio (a two-point temperature ratio; B18's 115 kJ/mol barrier would give ~16x) |
| pyrrole : 2,5(6)-DMP (molar) | 1 (phosphate 7, 130 C); 5 (bicarbonate 7); 38 (acetate 4, 130 C); 55 (acetate 4, 100 C, 30 min); 38 (2:1); 231 (1:3) | — | | Table 3 | within_study_ratio (at ~1 M the pyrrole is the main fate of aminoacetone whenever the pot is acid or MGO-rich) |
| 2,5(6)-DMP screen by amino acid | the 20-cell row above (Arg 493.9 ... Pro 1.6) | peak area x 1e-6 | neat DHA + amino acid 1:1, 90 C, 30 min, SPME | Table 1 | peak_area_only (the amino-acid reactivity order in the Strecker/aminoketone step, headspace-weighted) |
| aldehyde adducts of the dimethyldihydropyrazine | trimethyl (R = H): Arg 148.63, Gly 18.65, Thr 1.55; 3-ethyl-2,5-dimethyl (R = CH3): Arg 44.97, Lys 33.24, Ser 7.67, Ala 3.63, Gly 3.55, Cys 0.40; 2,5-dimethyl-3-isobutyl (R = iPr): Val 0.63; 2-ethyl-3,5-dimethyl: Ala 5.29 only | peak area x 1e-6 | same | Table 1 | peak_area_only (the amino-acid-specific adduct fires for Gly, Ala, Val; the same adducts from Arg, Lys, Ser come from sugar-derived formaldehyde and acetaldehyde) |
| adduct : 2,5(6)-DMP, same amino acid | trimethyl / DMP: Gly 0.056, Arg 0.30, Thr 0.07; 3-ethyl-2,5-DMP / DMP: Ala 0.043, Arg 0.091, Lys 0.086, Ser 0.053, Gly 0.011; 2-ethyl-3,5-DMP / DMP: Ala 0.062 | ratio of peak areas of different compounds | | Table 1 | peak_area_only; orientation only (the extra-substituent branch is a few % of the parent at 90 C, dry) |
| 2-ethyl-3,5-dimethyl : 3-ethyl-2,5-dimethyl from alanine | 1.46 (5.29 / 3.63) | ratio of peak areas | neat DHA, 90 C | Table 1 | peak_area_only (OPPOSITE to Cerny 1994's I : II = 10 from MGO + Ala at 180 C; flag 7) |
| non-detects that bound the rule | Leu, Ile, Met, Phe: no 3-alkyl-2,5-dimethylpyrazine from their Strecker aldehydes although the aldehydes were in the headspace; Pro: only 2,5(6)-DMP 1.6 | — | 90 C, 30 min, dry | Table 1 | peak_area_only (non-detects) |
| apparent B18-style constant (mine, orientation only) | 0.020 mol % DMP of 850 mmol/L in 30 min = 5.6e-3 mmol/(L·min); k_Strecker = 2 x rate / ([Ala][MGO]) ~ 1.6e-8 L/(mmol·min) at BOTH 100 C / pH 4 (acetate) and 130 C / pH 7 (phosphate) | L/(mmol·min) | assumes constant reactants; MGO at 1 M self-condenses (authors), so this is a lower bound | derived from Table 3 | derived; B18 stores 3.0e-8 at 100 C / pH 6.8 and would give ~5e-7 at 130 C / pH 7 and ~7e-10 at 100 C / pH 4: the 130 C pot is ~30x BELOW B18, the pH 4 pot ~20x ABOVE. Do not fit; the concentrations are 40x Zhou's and the MGO pool is not constant |

**The step in the authors' words (Scheme 1A).** Two aminoacetones (from the Strecker of MGO)
condense to 3,6-dimethyl-2,5-dihydropyrazine (1). Either it oxidises to 2,5-dimethylpyrazine, or a
ring proton is abstracted, the carbanion adds to an aldehyde R-CHO (the Strecker aldehyde, or a
sugar-fragment aldehyde), and water is eliminated to give 3-(CH2R)-2,5-dimethylpyrazine (5) with no
oxidation step. With an alpha-dicarbonyl in place of the aldehyde an acylpyrazine forms (Scheme 1B),
but the acylpyrazines seen are mostly disubstituted, so the authors prefer Scheme 2 for them
(DHA + aldehyde -> aminoketone 20 by an Amadori-type amination, + aminoacetone -> 2-acyl-5-methyl-
pyrazine). 2-Ethyl-6-methylpyrazine from glycine: MGO + formaldehyde aldol -> 2-oxo-3-butenal, its
aminoketone + aminoacetone -> the vinyl dihydropyrazine -> rearrangement (no oxidation). The
Scheme 3 pyrrole: aminoacetone + 2 MGO with one reduction step (reductones), competing with the
pyrazine for aminoacetone; favoured by acid, by MGO : amino acid >= 3, and by concentration.

## 5. Flags

1. **Table 1 is headspace SPME peak area, no internal standard, no replicates**: the DMP row orders the
   amino acids by (Strecker reactivity x headspace partition); cross-compound ratios inherit fibre
   selectivity. Record, never convert.
2. **Table 1's pot is a neat solid melt at 90 C** (5 + 5 mmol ground together); DHA dehydrates to MGO in
   situ; the "MGO" concentration is undefined and the temperature is 30-40 C below any cook the
   repository runs.
3. **Yield basis in Tables 2-3 is "% with respect to alanine", molar or mass unstated**; read as mol %
   (an internal-standard method with response factors is molar). Molar vs mass differ by MW ratio
   ~1.2 (DMP 108 / Ala 89), i.e. within the other uncertainties.
4. **Text/Methods temperature conflict for the aqueous alanine pot** (90 C in Results, 100 C in
   Methods and Table 2's caption). Table 3 states its temperatures explicitly.
5. **Concentrations are nominal**: "1 M" is 10 mmol in 10 mL buffer plus 1.8 mL (1:1) or ~5.4 mL (1:3)
   of 40 % MGO solution; ~850 mmol/L for 1:1. The 1 M phosphate buffer is itself a strong catalyst
   (Table 3 phosphate vs bicarbonate 1.8x on DMP) and did not hold pH 7 against Asp or Glu.
6. **The pH direction disagrees with B18's**: here acetate pH 4 gives 2.1x MORE 2,5(6)-DMP than
   phosphate pH 7 at 130 C, where B18's fitted slope below 7 (0.58 dec/pH, from Leahy's lysine +
   glucose whole-cascade ratios) predicts 55x LESS. A fed-MGO step and a sugar cascade need not share
   a pH law (the dicarbonyl supply is pH-driven in the cascade), and the buffers differ (acetate vs
   phosphate); but this is a directional conflict to carry.
7. **Isomer split inverts between laboratories**: from alanine at 90 C (dry DHA) the 3,5-isomer (II)
   area exceeds the 3,6-isomer (I) 1.5x; Cerny 1994 (MGO + Ala, 180 C, kieselguhr/oil) has I : II =
   10-19 by SIDA. Peak areas vs isotope dilution, dry vs roast, 90 vs 180 C — the reason is not
   identifiable from either paper.
8. **2,5- and 2,6-dimethylpyrazine co-elute** (DB-5); the repo's two keys get one number.
9. **Several identifications are tentative** (2-ethyl-6-methyl-, all acyl-pyrazines) and the
   "2-methyl-propylpyrazine" / "2-methyl-(E-1-propenyl)pyrazine" isomer positions are undetermined.
10. **No time series beyond 10 / 30 / 60 min at one condition**, and the 10 min DMP point is nd, so the
   DMP curve shape (induction vs linear) cannot be read.
11. **Registry naming trap**: Adams' "3-ethyl-2,5-dimethylpyrazine" and Cerny's "2-ethyl-3,6-
   dimethylpyrazine" are the same compound, and it is the structure the registry stores under the
   key for the OTHER isomer (`2_ethyl_3_5_dimethylpyrazine`). Any keying of this paper's two
   ethyldimethylpyrazine rows waits on the registry fix. Not edited here.

## 6. What the aldehyde-addition rule would look like (mine, on Adams' evidence)

**In words.** After B18's condensation of two aminoacetones the model holds 3,6-dimethyl-2,5-
dihydropyrazine (implicitly, since the condensation is declared fast and the product is booked as
2,5-dimethylpyrazine). The rule splits that intermediate: a fraction oxidises to 2,5-dimethylpyrazine
(the B18 product) and a fraction adds an aldehyde R-CHO present in the pot and dehydrates to
3-(CH2R)-2,5-dimethylpyrazine — formaldehyde -> trimethylpyrazine (`trimethylpyrazine`), acetaldehyde
-> 3-ethyl-2,5-dimethylpyrazine (= 2-ethyl-3,6-dimethyl-, Cerny's I, the weak odorant), 2-methylpropanal
-> 2,5-dimethyl-3-(2-methylpropyl)pyrazine. The aldehyde pool is the union of the amino acids' Strecker
aldehydes (glycine -> HCHO; alanine -> CH3CHO; valine -> iPrCHO; the trunk's Strecker step already
books formaldehyde from glycine) and the sugar-fragment aldehydes (formaldehyde and glycolaldehyde
from the retro-aldol of trioses; acetaldehyde), which is why arginine, lysine and serine give the
same adducts. The mixed dihydropyrazine (aminoacetone + 2-aminopropanal, or the glycine-derived
aminoacetaldehyde + aminoacetone) is the route to the potent 2-ethyl-3,5-dimethylpyrazine (Cerny II):
it needs the minor aminoketone and is 5-20x below the 3,6-isomer at 120-180 C (Shu, Cerny), though
not in Adams' 90 C dry pot. Branch magnitude Adams offers: adduct : parent DMP a few % at 90 C dry
(trimethyl / DMP 0.056 from glycine; 3-ethyl-2,5-DMP / DMP 0.043 from alanine), i.e. oxidation wins
unless the aldehyde is abundant. The rule should NOT fire for a peptide amine (Van Lancker 2012: no
Strecker aldehyde, adducts absent), and its Leu/Ile/Met/Phe adducts were absent at 90 C but present
at 130 C in water (Van Lancker Table 4), so it needs a temperature or solvent condition, not a flat
"any aldehyde".

**Positive control (should fire):** `CC1=NCC(C)=NC1` (3,6-dimethyl-2,5-dihydropyrazine) + `CC(C)C=O`
(2-methylpropanal, valine's Strecker aldehyde) -> `CC(C)Cc1nc(C)cnc1C` (2,5-dimethyl-3-(2-methylpropyl)-
pyrazine, Adams Table 1: Val 0.63) + H2O; and `CC1=NCC(C)=NC1` + `C=O` -> `Cc1cnc(C)c(C)n1`
(trimethylpyrazine, Adams: Gly 18.65) + H2O.

**Negative control (must not fire):** the aromatic pyrazine has no acidic ring CH2: `Cc1cnc(C)cn1`
(2,5-dimethylpyrazine) + `CC(C)C=O` -> no 2,5-dimethyl-3-(2-methylpropyl)pyrazine; and a secondary
amine gives no aminoketone to start from: proline + DHA -> no 3-substituted pyrazine at all (Adams
Table 1: proline's only pyrazine is 2,5(6)-DMP at 1.6, the lowest of the twenty), the proline route
going instead to ATHP (`hofmann1998b_extraction.md`).
