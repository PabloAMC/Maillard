# Yu & Ho 1995 — EXTRACTION (methionine and methionine sulfoxide, with and without glucose, 180 °C / 1 h in sealed water: methional against the dimethyl polysulfides)
### Single end point, four systems, one relative-quantified table. The methional : dimethyl-polysulfide split, and the sulfoxide switch that sends the sulfur to DMDS/DMTS instead.

**Source on disk:** `data/articles/yu1995.pdf` (owner's download, 2026-09-08). The `pdftotext` layer of
Table 2 was column-scrambled (compound names, MW, RI and the four value columns printed as separate
blocks); page 3 was re-extracted with `pypdf` in layout mode, which restored row alignment, and every
row below is from that pass. Figure 1 (bar chart of the same yields) not read. Table 3 (mass spectra of
tentatively identified compounds) not carried.

## 0. Identity

| field | value |
|---|---|
| Title | "Volatile Compounds Generated from Thermal Reaction of Methionine and Methionine Sulfoxide with or without Glucose" |
| Authors | Tung-Hsi Yu (Da-Yeh Institute of Technology, Taiwan), Chi-Tang Ho (Rutgers) |
| Venue | J. Agric. Food Chem. 1995, 43, 1641–1646; received 3 April 1995, accepted 12 April 1995 |
| DOI | 10.1021/jf00054a040 (not printed in the PDF; JF940522W is the manuscript code) |
| Systems | M = DL-methionine; M+G = + glucose; MS = DL-methionine sulfoxide; MS+G = + glucose; 0.01 mol each in 200 mL water, sealed, 180 °C, 1 h |
| Analytes | 62 volatiles, GC-FID relative to hexadecane, mg per mol amino acid; GC-MS identification |

## 1. Why it matters

Programme 6 wants methional, MeSH and the two disulfides from methionine. This is the one paper on disk
that quantifies methional, DMDS and DMTS from methionine in the same run, with and without a sugar, at
a cook temperature (180 °C), and repeats it with the sulfoxide. The pattern: from methionine the main
volatile is methional (M+G: 1127 mg/mol = 1.08 mol %) and DMDS is a third of it molar; from methionine
sulfoxide the order inverts — DMDS 2529 mg/mol (2.7 mol %) and DMTS 748 mg/mol (0.59 mol %) with
methional at 0.08 mol %. Glucose lifts every channel 6–26×. S-methyl methanethiosulfonate and dimethyl
tetrasulfide appear only from the sulfoxide, which is the sulfenic-acid route of Chin 1994's reactions
(g)–(h) seen at cook temperature. Also a MeSH sink that is not a disulfide: 1,1,3-tris(methylthio)propane
(methional + 2 MeSH) at 44 mg/mol in M+G. One time point; no MeSH; relative quantification.

## 2. Methods as they matter to a model

- **Charge:** 0.01 mol DL-methionine (Sigma > 99 %) or 0.01 mol DL-methionine sulfoxide (Aldrich 98 %),
  with or without 0.01 mol α-D-glucose (96 %), in **200 mL distilled water** → **50 mmol/L each**, 1:1.
  **No buffer; pH not stated or measured.**
- **Vessel / heating:** 0.3 L Hoke SS-DOT stainless-steel sample cylinders, sealed → ~100 mL headspace
  (air, presumably; not stated: ≈ 0.9 mmol O2 against 10 mmol methionine, my estimate). **180 °C in a GC
  oven, 1 h**; heat-up time not stated; cooled to room temperature.
- **Isolation:** whole reaction mass + 10 mL internal-standard solution (0.077 g hexadecane / 200 mL
  CH2Cl2 → 3.85 mg hexadecane per run), extracted with 200 mL redistilled CH2Cl2 for 3 h; dried (Na2SO4),
  Kuderna-Danish to minimum volume, N2 to 0.2 mL. ⚠ This work-up loses methanethiol (bp 6 °C), H2S and
  DMS; **MeSH is not reported anywhere in the paper.**
- **GC-FID:** Varian 3400, DB-1 60 m × 0.25 mm × 1 µm, 40 °C (5 min) → 260 °C at 2 °C/min (60 min);
  split 50:1. **Quantities as mg per mol of methionine (or sulfoxide) by peak area against the internal
  standard** — i.e. FID response factors taken as 1 (not stated otherwise), no calibration per compound.
- **GC-MS:** Finnigan MAT 8230 magnetic sector, EI 70 eV, same column; ID by library MS + Kovats RI;
  some compounds tentative (Table 3).
- **Replicates:** none stated.
- **Odor:** sniffed descriptions only (Table 1).
- **Conversions used below:** mol % = (mg/mol) ÷ MW ÷ 10. MW: DMDS 94.20; methional 104.15; DMTS 126.26;
  dimethyl tetrasulfide 158.32; S-methyl methanethiosulfonate 126.19; (methylthio)cyclopentane 116.22;
  2-methyl-5-(methylthio)furan 128.19; 2,6-dithianonane 164.32; 1,1,3-tris(methylthio)propane 182.35;
  4-(methylthio)-2-butanone 118.20; 1,1-bis(methylthio)ethane 122.25.

## 3. Tables re-typed

### Table 1. Odor descriptions of the reaction solutions

| sample | odor |
|---|---|
| methionine | fermented radish or cabbage note with baked potato undertone |
| methionine + glucose | fermented radish or cabbage note with burned, caramellic and baked potato undertone |
| methionine sulfoxide | sulfuryl, fermented radish-like, with black mushroom undertone |
| methionine sulfoxide + glucose | sulfuryl, fermented radish-like, burned, caramellic with black mushroom and baked potato undertone |

### Table 2. "Volatile Compounds Identified from the Thermal Reaction of Methionine or Methionine Sulfoxide with or without Glucose" — mg per mol of methionine or methionine sulfoxide; nd = not detected

**Compounds generated from the thermal degradation of methionine or methionine sulfoxide**

| compound | MW | RI | M | M+G | MS | MS+G |
|---|---:|---:|---:|---:|---:|---:|
| **dimethyl disulfide** | 94 | 731 | **19.1** | **332.2** | **96.1** | **2528.7** |
| **methional** | 104 | 875 | **90.5** | **1126.9** | **14.8** | **87.6** |
| (methylthio)cyclopentane | 116 | 943 | 17.1 | nd | 23.9 | nd |
| **dimethyl trisulfide** | 126 | 957 | **2.5** | **nd** | **39.1** | **748.4** |
| S-methyl methylthiosulfonate | 126 | 1092 | nd | nd | 3.1 | 15.6 |
| 2-methyl-5-(methylthio)furan | 128 | 1196 | 1.5 | nd | 3.4 | 15.8 |
| dimethyl tetrasulfide | 158 | 1206 | nd | nd | 4.5 | 24.8 |
| 2,6-dithianonane | 164 | 1316 | 2.2 | nd | 0.9 | nd |
| 1,1,3-tris(methylthio)propane | 182 | 1510 | nd | 44.1 | nd | nd |
| 2-[(methylthio)methyl]-5-(methylthio)-2-pentenal | 190 | 1567 | nd | 72.5 | nd | nd |

**Compounds generated from the thermal degradation of glucose**

| compound | MW | RI | M | M+G | MS | MS+G |
|---|---:|---:|---:|---:|---:|---:|
| 2-butanone | 72 | 564 | nd | 18.6 | nd | 8.9 |
| 1-hydroxy-2-propanone | 74 | 625 | nd | 76.8 | nd | 61.6 |
| acetoin | 88 | 690 | nd | 2.8 | nd | 27.7 |
| 3-hydroxy-2-pentanone | 102 | 783 | nd | 13.7 | nd | 14.1 |
| 2,3-pentanedione | 100 | 790 | nd | 8.1 | nd | 5.8 |
| furfural | 96 | 829 | nd | 14.5 | nd | 17.6 |
| furfuryl alcohol | 98 | 849 | nd | 9.1 | nd | 41.9 |
| 5-methyl-2-furfural | 110 | 950 | nd | 3.6 | nd | 21.7 |
| 2-hydroxycyclohexanone | 114 | 971 | nd | 13.3 | nd | 16.7 |
| phenol | 94 | 981 | nd | 15.2 | nd | 24.1 |
| cyclotene | 112 | 1005 | nd | 56.6 | nd | 117.1 |
| 2-hydroxy-3-methylbenzaldehyde | 136 | 1125 | nd | 9.9 | nd | 7.6 |
| 2,3-dihydro-3,5-dihydroxy-6-methyl-4H-pyran-4-one | 144 | 1144 | nd | nd | nd | 9.1 |
| 5-(hydroxymethyl)furfural | 126 | 1228 | nd | nd | nd | 13.3 |

**Compounds generated from the thermal interactions of glucose and methionine or methionine sulfoxide**

| compound | MW | RI | M | M+G | MS | MS+G |
|---|---:|---:|---:|---:|---:|---:|
| ethanethioic acid S-methyl ester | 90 | 671 | nd | 3.4 | nd | 1.5 |
| pyrazine | 80 | 716 | nd | 97.7 | nd | 80.4 |
| pyridine | 79 | 753 | nd | nd | nd | 13.7 |
| methylpyrazine | 94 | 802 | nd | 94.1 | nd | 170.7 |
| 3-methylpyridine | 93 | 845 | nd | 23.9 | 6.5 | 86.3 |
| 2,5-dimethylpyrazine | 108 | 891 | nd | 165.1 | nd | 192.5 |
| ethylpyrazine | 108 | 895 | nd | 1.38 | nd | 55.6 |
| 2,3-dimethylpyrazine | 108 | 899 | nd | 66.2 | nd | 63.2 |
| vinylpyrazine | 106 | 922 | nd | 0.3 | nd | 18.1 |
| 1,1-bis(methylthio)ethane | 122 | 935 | nd | 2.4 | nd | nd |
| 3-hydroxy-2-thiabutane | 92 | 941 | nd | 1.3 | nd | nd |
| 4-(methylthio)-2-butanone | 118 | 962 | nd | 64.5 | nd | nd |
| 2-ethyl-5-methylpyrazine | 122 | 985 | nd | 75.7 | nd | 160.8 |
| 2-formylthiophene | 112 | 987 | nd | nd | 1.3 | nd |
| trimethylpyrazine | 122 | 992 | nd | 128.9 | nd | 91.7 |
| 2-vinyl-5-methylpyrazine | 120 | 999 | nd | 0.6 | nd | 13.2 |
| 2-acetylpyrrole | 109 | 1051 | nd | 20.1 | nd | 12.3 |
| 2-ethyl-3,5-dimethylpyrazine | 126 | 1065 | nd | 17.8 | nd | 6.6 |
| 3-ethyl-2,5-dimethylpyrazine | 126 | 1071 | nd | 12.3 | nd | 35.1 |
| 2-methyl-5-[(methylthio)methyl]furan | 142 | 1080 | nd | 32.5 | nd | nd |
| 2-methyl-5-(1-propenyl)pyrazine | 134 | 1110 | nd | 7.2 | nd | 6.8 |
| 2-formyl-5-methylthiophene | 126 | 1129 | nd | nd | 0.8 | 2.3 |
| 3-(methylthio)propanoic acid methyl ester | 134 | 1142 | nd | 119.6 | nd | nd |
| 3,5-diethyl-2-methylpyrazine | 150 | 1152 | nd | nd | nd | 4.4 |
| 3,4-dihydro-2H-thiopyran-3-one | 114 | 1155 | nd | nd | 0.7 | nd |
| 3-(methylthio)propanoic acid ethyl ester | 148 | 1167 | nd | 23.1 | nd | nd |
| 2-formylpyrrole | 95 | 1234 | nd | nd | nd | 20.3 |
| 3-[(methylthio)methyl]pyridine | 139 | 1243 | nd | 37.3 | nd | nd |
| 5,6,7,8-tetrahydroquinoline | 133 | 1257 | nd | nd | nd | 31.5 |
| 2-pyridinemethanol | 109 | 1305 | nd | nd | nd | 24.4 |
| 5-[(methylthio)methyl]furfuryl alcohol | 158 | 1335 | nd | nd | nd | 6.5 |
| 2-[(methylthio)methyl]-3-(2-furyl)acrolein (isomer 1) | 182 | 1456 | nd | 25.1 | nd | nd |
| 2-[(methylthio)methyl]-3-(2-furyl)acrolein (isomer 2) | 182 | 1468 | nd | 17.1 | nd | nd |
| 2-[(methylthio)methyl]-3-(3-furyl)acrolein (isomer 1) | 182 | 1480 | nd | 44.1 | nd | nd |
| 2-[(methylthio)methyl]-3-(3-furyl)acrolein (isomer 2) | 182 | 1488 | nd | 11.3 | nd | nd |
| 1-[3-(methylthio)propyl]-2-formylpyrrole | 183 | 1504 | nd | 103.5 | nd | nd |
| 2,3-dimethyl-5-[(methylthio)propyl]pyrazine | 196 | 1556 | nd | 10.3 | nd | nd |
| 2,5-dimethyl-3-[(methylthio)propyl]pyrazine | 196 | 1601 | nd | 31.3 | nd | nd |

Footnotes: MW molecular weight; RI Kovats retention index (DB-1); nd not detected. "ethylpyrazine 1.38" is
printed with two decimals where the rest of the table has one.

### Derived from Table 2 (my arithmetic)

| quantity | M | M+G | MS | MS+G |
|---|---:|---:|---:|---:|
| DMDS, mmol/mol (mol %) | 0.203 (0.020 %) | 3.53 (0.35 %) | 1.02 (0.10 %) | 26.8 (2.7 %) |
| methional, mmol/mol (mol %) | 0.869 (0.087 %) | 10.8 (1.08 %) | 0.142 (0.014 %) | 0.841 (0.084 %) |
| DMTS, mmol/mol (mol %) | 0.020 (0.0020 %) | nd | 0.310 (0.031 %) | 5.93 (0.59 %) |
| dimethyl tetrasulfide, mmol/mol | nd | nd | 0.028 | 0.157 |
| S-methyl methanethiosulfonate, mmol/mol | nd | nd | 0.025 | 0.124 |
| 1,1,3-tris(methylthio)propane, mmol/mol | nd | 0.242 (= 0.48 mmol MeSH-eq) | nd | nd |
| methional / DMDS (molar) | 4.3 | 3.1 | 0.14 | 0.031 |
| DMTS / DMDS (molar) | 0.098 | < (nd) | 0.30 | 0.22 |
| glucose effect, DMDS | — | × 17.4 | — | × 26.3 |
| glucose effect, methional | — | × 12.5 | — | × 5.9 |
| glucose effect, DMTS | — | nd from 2.5 | — | × 19.1 |
| sulfoxide effect (no glucose), DMDS / DMTS / methional | — | — | × 5.0 / × 15.6 / × 0.16 | — |
| S in methional + 2·DMDS + 2·DMTS + 2·S4, % of amino-acid S | 0.13 % | 1.8 % | 0.29 % | 6.7 % |

## 4. Numbers the repository can use

Registry keys: `methional`, `dimethyl_disulfide`, `dimethyl_trisulfide`, `furfural`, `methylpyrazine`,
`2_5_dimethylpyrazine`, `2_3_dimethylpyrazine`, `2_ethylpyrazine`, `trimethylpyrazine`,
`2_ethyl_3_5_dimethylpyrazine`, `pyrazines` exist; `methanethiol` exists but is not measured here.
Not in registry: methionine, methionine sulfoxide, dimethyl tetrasulfide, S-methyl methanethiosulfonate,
1,1,3-tris(methylthio)propane, pyrazine (parent).

| quantity | value | unit | conditions | source | evidence class | registry key |
|---|---|---|---|---|---|---|
| methional from Met alone | 90.5 (0.087 mol %) | mg/mol Met | 50 mM Met, water, unbuffered, sealed steel, 180 °C, 1 h | Table 2 | level_only | methional |
| methional from Met + glucose 1:1 | 1126.9 (1.08 mol %) | mg/mol Met | same + 50 mM glucose | Table 2 | level_only | methional |
| DMDS from Met alone / Met + glucose | 19.1 / 332.2 (0.020 / 0.35 mol %) | mg/mol Met | same | Table 2 | level_only | dimethyl_disulfide |
| DMTS from Met alone / Met + glucose | 2.5 / nd (0.0020 mol % / nd) | mg/mol Met | same | Table 2 | level_only | dimethyl_trisulfide |
| DMDS from MetSO alone / MetSO + glucose | 96.1 / 2528.7 (0.10 / 2.7 mol %) | mg/mol MetSO | same, sulfoxide | Table 2 | level_only | dimethyl_disulfide |
| DMTS from MetSO alone / MetSO + glucose | 39.1 / 748.4 (0.031 / 0.59 mol %) | mg/mol MetSO | same | Table 2 | level_only | dimethyl_trisulfide |
| methional from MetSO alone / + glucose | 14.8 / 87.6 | mg/mol MetSO | same | Table 2 | level_only | methional |
| methional : DMDS molar, Met systems | 4.3 (M), 3.1 (M+G) | — | 180 °C, 1 h | derived | within_study_ratio | methional, dimethyl_disulfide |
| DMTS : DMDS molar | 0.10 (M); 0.30 (MS); 0.22 (MS+G); nd (M+G) | — | same | derived | within_study_ratio | both |
| glucose multiplier on Met → methional / DMDS | 12.5 / 17.4 | — | 1:1 glucose, 180 °C, 1 h | derived | within_study_ratio | methional, dimethyl_disulfide |
| MetSO vs Met, DMDS / DMTS / methional | × 5.0 / × 15.6 / × 0.16 (no glucose); × 7.6 / — / × 0.078 (with glucose) | — | same | derived | within_study_ratio | as above |
| 1,1,3-tris(methylthio)propane (methional·2 MeSH dithioacetal) | 44.1 (0.24 mmol/mol) | mg/mol Met | M+G only | Table 2 | level_only | not in registry |
| S-methyl methanethiosulfonate | 3.1 / 15.6 (MS / MS+G); nd from Met | mg/mol | same | Table 2 | level_only | not in registry |
| pyrazine ladder, M+G | pyrazine 97.7, methyl- 94.1, 2,5-dimethyl- 165.1, 2,3-dimethyl- 66.2, trimethyl- 128.9, 2-ethyl-5-methyl- 75.7, 2-ethyl-3,5-dimethyl- 17.8 | mg/mol Met | same | Table 2 | level_only | pyrazines etc. |

## 5. Flags

1. **One temperature (180 °C), one time (1 h), no replicates.** Nothing kinetic; the glucose and
   sulfoxide multipliers are single-run ratios.
2. **Relative FID quantification** against hexadecane with response factors implicitly 1, after a
   CH2Cl2 / Kuderna-Danish work-up that discriminates against volatile and polar compounds. Levels of
   DMDS (bp 110 °C) are certainly low-biased by the concentration step; ratios between compounds of
   similar volatility (methional 165 °C, DMTS 170 °C) are safer than DMDS ratios.
3. **Methanethiol is not measured** (lost in the work-up), nor H2S or DMS.
4. **Unbuffered water, pH unknown**; 50 mM methionine alone sits near its isoelectric pH (~5.7, mine),
   glucose acids will lower it during the hour. The methional → MeSH elimination is base-sensitive
   (Schutte 1972), so the split reported here is at an undefined, probably falling, pH.
5. **Stainless-steel vessel**: Fe/Cr/Ni surfaces at 180 °C are a possible oxidation catalyst for MeSH →
   DMDS (Chin 1994 shows metal dependence at 30 °C); the oxidant inventory (air headspace ~100 mL) is
   not stated.
6. **DMTS "nd" with glucose but 2.5 mg/mol without** while DMDS rises 17× — inconsistent with DMTS
   forming from DMDS or MeSH alone. Consistent with Chin 1994's finding that DMTS needs a second sulfur
   (H2S / sulfane), which methionine does not supply and which glucose would not add; the sulfoxide,
   via methanesulfenic acid and its thiosulfinate/thiosulfonate, does.
7. **The sulfoxide result is a warning for the isolate programme**: protein-bound methionine in a
   plant isolate is partly oxidised (the paper cites MetSO in soy concentrate); its heating gives 26×
   more DMDS and DMTS-forming sulfur than reduced methionine, with almost no methional. A model charging
   "methionine" should know the oxidation state.
8. **A MeSH sink that is not a disulfide**: 1,1,3-tris(methylthio)propane at 0.24 mmol/mol (0.48 mmol
   MeSH-equivalents per mol Met) in M+G is 7 % of the DMDS-bound MeSH there; 1,1-bis(methylthio)ethane
   (ethanal dithioacetal, cf. Schutte 1972) and 3-hydroxy-2-thiabutane are the same chemistry with
   acetaldehyde. These need a carbonyl partner, which a Maillard pot has in excess.
