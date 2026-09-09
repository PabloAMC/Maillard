# Van Lancker 2012 — EXTRACTION (nine dipeptides and three tripeptides against their free amino acids with glucose, methylglyoxal or glyoxal; water, pH 8, 130 C / 2 h; pyrazines as SBSE-GC-MS peak areas)
### Peptides as the amine: the N-terminal residue decides; Gly-, Ala-, Ser- and Lys-terminated dipeptides make MORE 2,5(6)-dimethyl- and trimethylpyrazine than their free amino acids, Pro-, Val-, Leu-terminated ones far less; unsubstituted pyrazine and the amino-acid-specific pyrazines always favour the free amino acids.

**Source on disk:** `data/articles/vanlancker2012.pdf` (owner's download, 2026-09-08). Read-only
extraction from the pypdf text layer in the scratchpad; the twelve-column Tables 4-6 were re-checked
row by row against `pdftotext -layout` (the layer merges adjacent dashes, which the layout view
resolves). All nine pyrazine tables and Table 10 are re-typed below; the amino-acid-specific rows of the
wide glucose Table 4 are given as non-zero cells. EVERYTHING in this paper is a GC-MS peak area (SBSE,
or SPME for Table 10): peak_area_only throughout, recorded as printed and never converted. Repo status
before this dossier: roadmap §5d asks what a protein isolate's peptides offer as the amine.

## 0. Identity

| field | value |
|---|---|
| Title | "Impact of the N-Terminal Amino Acid on the Formation of Pyrazines from Peptides in Maillard Model Systems" |
| Authors | Fien Van Lancker, An Adams, Norbert De Kimpe (Ghent University) |
| Venue | J. Agric. Food Chem. 2012, 60, 4697-4708 |
| DOI | 10.1021/jf301315b |
| Naming | 2,5(6)-dimethylpyrazine = the unresolved 2,5-/2,6- pair on HP5-MS; "X + Y" = the equimolar mixture of free amino acids; "XY" = the dipeptide; "theor recovery (%)" = the theoretical SBSE (PDMS Twister) recovery of that compound from water, printed by the authors |
| Companion | Van Lancker, Adams, De Kimpe 2010, JAFC 58, 2470 (Lys-X dipeptides; the Methods this paper cites for amounts and volumes; not on disk) |

## 1. Why it matters

A protein isolate brings almost no free amino acid and a great deal of N-terminus and lysine side
chain (roadmap §5d). B18's Strecker step needs a free alpha-amino acid (decarboxylation). This paper
shows, with the same dicarbonyls the trunk carries (MGO, glyoxal) and with glucose, that peptides
make pyrazines WITHOUT Strecker decarboxylation — the aminoketone forms by an imine / 1,5-H-shift /
2-azadiene hydrolysis route (their Scheme 1) that leaves the peptide's carbon skeleton intact and
therefore gives no Strecker aldehyde. Consequences the data bear out: (i) with MGO, Gly-, Ala-, Ser-,
Lys-terminated dipeptides give 4-34x the 2,5(6)-dimethylpyrazine area of the free amino acids;
(ii) the amino-acid-specific alkylpyrazines (3-ethyl-2,5-dimethyl from Ala, isobutyl from Val,
isopentyl from Leu) are largely absent with peptides, so the aldehyde-addition step is starved when the
amine is a peptide; (iii) unsubstituted pyrazine (from glyoxal) is always LOWER with peptides;
(iv) N-terminal proline gives nothing; N-terminal Val/Leu little; tripeptides less than dipeptides;
(v) the peptide advantage inverts above 150 C (Table 10). Directional, peak-area evidence only.

## 2. Methods as they matter to a model

- **Charges (from the table headers and Table 10):** dipeptide **1 mmol**, or **1 mmol of each** free amino
  acid (the "Gly (2 mmol)" arm of Tables 4-6 is the free-glycine control for GlyGly; "Gly (3 mmol)" for
  GlyGlyGly). So the free-amino-acid arm carries **twice** (dipeptides) or three times (tripeptides) the
  alpha-amino groups of the peptide arm. **The glucose amount, the dicarbonyl amount and the water
  volume are given only by reference to Van Lancker 2010 (not on disk); this paper states only that the
  dicarbonyls were used at "a 10-fold lower concentration" than glucose "to avoid too many
  self-condensation reactions".** No mmol/L can be assigned from this paper.
- **Conditions:** "unbuffered aqueous conditions at pH 8 at 130 C for 2 h" (initial pH; no buffer, on
  purpose, because phosphate and carboxylate catalyse). Vessel not described here.
- **Sampling:** Stir Bar Sorptive Extraction (SBSE), 30 min at 35 C, thermal desorption (Gerstel TDS2);
  GC-MS Agilent 6890 / 5973, HP5-MS 30 m x 0.25 mm x 0.25 µm. Identification by LRI (experimental vs
  literature) and MS; "tentatively identified" where marked.
- **Quantity reported:** **GC-MS peak area x 1e8** (Tables 1-9), **x 1e6** (Table 10, SPME). No internal
  standard, no calibration. The "theor recovery" column (pyrazine 0.3 %, methylpyrazine 0.8,
  dimethylpyrazine 2.0/1.6, 2-ethylpyrazine 2.3, trimethylpyrazine 4.1, tetramethylpyrazine 8.4,
  2,3-diethyl-5-methylpyrazine 30.0, 2-acetylpyrazine 0.8) says that areas of DIFFERENT compounds are
  not comparable (a pyrazine area under-represents its amount ~14x relative to a trimethylpyrazine
  area); areas of ONE compound across arms are comparable if the matrices are alike.
- **Table 10 (temperature series):** glycine 2 mmol or diglycine 1 mmol + glucose, 2 h, pH 8, at 100 /
  130 / 150 / 180 C; volatiles by **SPME** (PDMS/CAR/DVB, 30 min, 35 C) — the authors note its 130 C
  values therefore differ from Table 4's.
- **Replicates:** none stated.

## 3. Tables re-typed

Blank cell = not detected ("−"); tr = trace; areas x 1e8 unless noted. "rec" = theoretical SBSE recovery (%).

### Table 1. Glucose + X-Lys dipeptides vs free amino acids (2 h, 130 C)

| compound (LRI exptl) | Gly + Lys | GlyLys | Ala + Lys | AlaLys | Val + Lys | ValLys | rec |
|---|---:|---:|---:|---:|---:|---:|---:|
| pyrazine (759) | 1.46 | 0.01 | 1.70 | 0.01 | 0.80 | 0.02 | 0.3 |
| methylpyrazine (822) | 1.30 | 0.03 | 1.00 | 0.06 | 0.32 | 0.01 | 0.8 |
| 2,5(6)-dimethylpyrazine (908) | 1.53 | 4.79 | 0.92 | 4.94 | 0.41 | 0.07 | 2.0/1.6 |
| 2-ethylpyrazine (911) | 0.01 | | 0.09 | | | | 2.3 |
| 2,3-dimethylpyrazine (913) | 0.19 | 0.32 | 0.14 | 0.07 | 0.07 | 0.01 | 1.6 |
| 2-ethyl-6-methylpyrazine (995) | | 0.01 | 0.04 | 0.01 | | | |
| 2-ethyl-5-methylpyrazine (998) | 0.06 | 3.03 | 0.17 | 0.99 | 0.10 | | |
| trimethylpyrazine (998) | 0.58 | 0.51 | 0.16 | 0.21 | 0.04 | | 4.1 |
| 2-ethenyl-6-methylpyrazine (1015) | 0.01 | 0.01 | | 0.01 | | | |
| 2-ethenyl-5-methylpyrazine (1018) | 0.01 | 0.36 | 0.01 | 0.30 | 0.01 | | |
| 2-(2-methylpropyl)pyrazine (1064, tent.) | | | | | 0.12 | | |
| 3-ethyl-2,5-dimethylpyrazine (1072) | 0.09 | 0.13 | 1.73 | 0.25 | 0.04 | | |
| 2-ethyl-3,5-dimethylpyrazine (1080) | 0.01 | 0.10 | 0.01 | 0.01 | | | |
| tetramethylpyrazine (1082) | 0.02 | | | | | | 8.4 |
| 5-ethyl-2,3-dimethylpyrazine (1082) | 0.02 | 0.41 | 0.02 | 0.02 | 0.01 | | |
| 2,5-diethylpyrazine (1087) | | 0.13 | | 0.01 | | | |
| 3-ethenyl-2,5-dimethylpyrazine (1093) | 0.02 | 0.42 | 0.03 | 0.45 | 0.01 | | |
| methyl-(2-methylpropyl)pyrazine (1142) | | | | | 0.01 | | |
| 2,3-diethyl-5-methylpyrazine (1165) | | 0.01 | 0.06 | | | | 30.0 |
| 3,5-diethyl-2-methylpyrazine (1167) | tr | 0.02 | 0.12 | 0.01 | 0.01 | | |
| 2,3,5-trimethyl-6-ethylpyrazine (1169) | 0.03 | 0.02 | 0.27 | | | | |
| acetylethylpyrazine (1178, tent.) | | 0.01 | | | | | |
| 2,5-dimethyl-3-(2-methylpropyl)pyrazine (1197) | | | | | 0.85 | | |
| 2-(2'-furyl)pyrazine (1255, tent.) | 0.05 | | 0.06 | | 0.02 | | |
| (2-methylpropyl)trimethylpyrazine (1276, tent.) | | | | | 0.12 | | |
| 1,4-dimethylpyrrolo(1,2a)pyrazine (1382, tent.) | 0.10 | | 0.18 | | | | |
| **total pyrazines** | **5.48** | **10.33** | **6.71** | **7.36** | **2.96** | **0.11** | |
| pyrazines, % of total peak area | 40.1 | 86.7 | 26.3 | 79.8 | 18.1 | 3.2 | |

### Table 2. Methylglyoxal + X-Lys dipeptides vs free amino acids (2 h, 130 C)

| compound (LRI) | Gly + Lys | GlyLys | Ala + Lys | AlaLys | Val + Lys | ValLys | rec |
|---|---:|---:|---:|---:|---:|---:|---:|
| methylpyrazine (822) | 0.03 | | 0.05 | 0.04 | 0.03 | | 0.8 |
| 2,5(6)-dimethylpyrazine (908) | 2.21 | 8.31 | 3.32 | 16.81 | 2.79 | 3.86 | 2.0/1.6 |
| 2,3-dimethylpyrazine (913) | 0.75 | tr | 0.07 | 0.33 | 0.04 | tr | 1.6 |
| 2-ethyl-5-methylpyrazine (998) | 0.03 | 0.13 | 0.01 | 0.21 | 0.01 | 0.03 | |
| trimethylpyrazine (998) | 3.40 | 5.19 | 0.66 | 2.90 | 0.51 | 0.18 | 4.1 |
| 3-ethyl-2,5-dimethylpyrazine (1072) | 0.14 | 0.29 | 3.34 | 2.01 | 0.19 | 0.09 | |
| 2-ethyl-3,5-dimethylpyrazine (1080) | 0.03 | 0.07 | | | 0.05 | | |
| tetramethylpyrazine (1082) | 0.10 | 0.23 | | | | | 8.4 |
| 5-ethyl-2,3-dimethylpyrazine (1082) | 0.28 | 0.05 | tr | 0.03 | | | |
| 3-ethenyl-2,5-dimethylpyrazine (1093) | | | | 0.05 | | | |
| 2-acetyl-5-methylpyrazine (1123) | 0.02 | 0.10 | 0.03 | 0.07 | 0.04 | | |
| 2-acetyl-6-methylpyrazine (1133) | 0.03 | tr | 0.04 | 0.07 | 0.02 | | |
| 2,3-diethyl-5-methylpyrazine (1165) | | tr | 0.02 | tr | | | 30 |
| 3,5-diethyl-2-methylpyrazine (1167) | | | 0.02 | tr | | | |
| 2,5-dimethyl-3-propylpyrazine (1169) | | | | | 0.08 | | |
| 2,3,5-trimethyl-6-ethylpyrazine (1169) | 0.02 | | 0.08 | | | | |
| acetyldimethylpyrazine (1180, tent.) | | 0.01 | 0.01 | 0.08 | 0.01 | | |
| 3,5-dimethyl-2-(2-methylpropyl)pyrazine (1211) | | | | | 1.46 | 0.01 | |
| acetyldimethylpyrazine (1219, tent.) | 0.05 | 0.02 | 0.04 | 0.14 | 0.01 | | |
| 2,3-dimethyl-5-(2-methylpropyl)pyrazine (1229) | | 0.06 | 0.02 | 0.21 | 0.02 | | |
| 2,5-dimethyl-3-(E-1-propenyl)pyrazine (1241) | | 0.08 | | 1.13 | | | |
| trimethyl-(2-methylpropyl)pyrazine (1276, tent.) | | | | | 0.03 | | |
| **total pyrazines** | **7.13** | **14.53** | **7.75** | **24.09** | **5.31** | **4.17** | |
| pyrazines, % of total peak area | 62.5 | 58.9 | 53.8 | 53.1 | 25.9 | 45.3 | |

### Table 3. Glyoxal + X-Lys dipeptides vs free amino acids (2 h, 130 C)

| compound (LRI) | Gly + Lys | GlyLys | Ala + Lys | AlaLys | Val + Lys | ValLys | rec |
|---|---:|---:|---:|---:|---:|---:|---:|
| pyrazine (759) | 3.64 | 1.50 | 3.58 | 2.63 | 3.59 | 0.36 | 0.3 |
| methylpyrazine (822) | 0.06 | 0.04 | 0.04 | 0.09 | 0.11 | | 0.8 |
| 2-ethylpyrazine (911) | | | 0.30 | | | | 2.3 |
| 2-(2-methylpropyl)pyrazine (1064, tent.) | | | | | 0.37 | | |
| **total pyrazines** | **3.70** | **1.54** | **3.92** | **2.71** | **4.07** | **0.36** | |
| pyrazines, % of total peak area | 75.7 | 93.9 | 66.6 | 89.1 | 26.6 | 54.4 | |

### Table 4. Glucose + X-Gly dipeptides vs free amino acids (2 h, 130 C) — registry-relevant rows and totals

| compound (LRI) | Gly (2 mmol) | GlyGly | Ala + Gly | AlaGly | Val + Gly | ValGly | Leu + Gly | LeuGly | Ser + Gly | SerGly | Pro + Gly | ProGly | rec |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| pyrazine (759) | 0.02 | | 0.05 | | 0.01 | | 0.02 | | 0.08 | | 0.03 | | 0.3 |
| methylpyrazine (822) | 0.02 | 0.02 | 0.05 | 0.03 | 0.02 | | 0.04 | 0.03 | 0.09 | 0.09 | 0.06 | | 0.8 |
| 2,5(6)-dimethylpyrazine (908) | 0.11 | 2.38 | 0.57 | 3.29 | 0.16 | | 0.51 | 0.57 | 0.18 | 7.19 | 0.03 | | 2.0/1.6 |
| 2,3-dimethylpyrazine (913) | 0.04 | 0.20 | 0.04 | 0.06 | 0.05 | | 0.18 | 0.13 | 0.09 | 0.38 | | | 1.6 |
| ethenylpyrazine (926) | 0.02 | | 0.02 | 0.01 | | | 0.02 | tr | 0.01 | | 0.02 | | |
| 2-ethyl-6-methylpyrazine (995) | tr | 0.01 | 0.03 | 0.01 | tr | | 0.01 | tr | tr | 0.01 | | | |
| trimethylpyrazine (998) | 0.19 | 0.85 | 0.29 | 0.15 | 0.14 | | 0.27 | | 0.21 | 1.06 | 0.03 | | 4.1 |
| 2-ethyl-5-methylpyrazine (998) | 0.02 | 1.68 | 0.13 | 1.63 | 0.05 | | 0.17 | 0.23 | 0.03 | 2.55 | | | |
| 2-ethenyl-6-methylpyrazine (1015) | | | | 0.01 | | | tr | | | 0.03 | | | |
| 2-acetylpyrazine (1018) | | | | | | | | | 0.02 | | | | 0.8 |
| 2-ethenyl-5-methylpyrazine (1018) | 0.01 | 0.24 | 0.01 | 0.71 | | | 0.04 | 0.16 | 0.03 | 0.76 | | | |
| 3-ethyl-2,5-dimethylpyrazine (1072) | tr | 0.04 | 0.93 | 0.05 | 0.02 | | 0.01 | | 0.02 | 0.04 | | | |
| 2-ethyl-3,5-dimethylpyrazine (1080) | | 0.08 | 0.02 | 0.03 | tr | | 0.02 | | 0.01 | 0.12 | | | |
| tetramethylpyrazine (1082) | 0.12 | | 0.04 | | 0.01 | | 0.04 | | 0.03 | | 0.01 | | 8.4 |
| 5-ethyl-2,3-dimethylpyrazine (1082) | 0.01 | 0.24 | 0.06 | 0.05 | 0.03 | | 0.05 | 0.01 | 0.03 | 0.24 | | | |
| 2,5-diethylpyrazine (1087) | | 0.23 | 0.01 | 0.10 | | | tr | 0.01 | | 0.08 | | | |
| 3-ethenyl-2,5-dimethylpyrazine (1093) | tr | 0.43 | 0.03 | 0.40 | | | 0.01 | 0.05 | 0.01 | 1.79 | | | |
| 5-ethenyl-2,3-dimethylpyrazine (1115) | | 0.19 | | 0.17 | | | | 0.06 | | 0.15 | | | |
| 2-acetyl-5-methylpyrazine (1123) | | | | | | | | | 0.01 | 0.03 | | | |
| 2,3-diethyl-5-methylpyrazine (1165) | | tr | 0.04 | tr | | | 0.01 | | | | | | 30.0 |
| 3,5-diethyl-2-methylpyrazine (1167) | | 0.01 | 0.10 | | | | 0.01 | | | | | | |
| 2,3,5-trimethyl-6-ethylpyrazine (1169) | 0.02 | | 0.31 | | 0.01 | | | | 0.01 | | | | |
| **total pyrazines** | **0.57** | **6.59** | **2.72** | **6.70** | **1.48** | **0.00** | **8.83** | **1.27** | **0.84** | **14.52** | **0.18** | **0.00** | |
| pyrazines, % of total peak area | 60.7 | 87.2 | 50.5 | 84.1 | 12.0 | 0.0 | 16.4 | 13.7 | 54.5 | 90.3 | 1.1 | 0.0 | |

Amino-acid-specific rows of Table 4 (non-zero cells only; all tentatively identified unless an LRI
reference is given): 2-(3-methylbutyl)pyrazine (1184) Leu + Gly 0.17; 3-isobutyl-2,5-dimethylpyrazine
(1197) Val + Gly 0.72; methyl-(3-methylbutyl)pyrazines (1250 / 1253 / 1265) Leu + Gly 0.08 / 0.28 / 0.22,
LeuGly tr at 1265; trimethyl-(2-methylpropyl)pyrazine (1276) Val + Gly 0.26; dimethyl-(3-methylbutyl)-
pyrazines (1315 / 1327 / 1337) Leu + Gly **4.40** / 0.03 / 0.05, LeuGly 0.02 at 1315;
trimethyl-(2-methylbutyl)pyrazine (1383) Leu + Gly 0.70; trimethyl-(3-methylbutyl)pyrazine (1387)
Leu + Gly 1.49. The leucine arm's total (8.83) is mostly these.

### Table 5. Methylglyoxal + X-Gly dipeptides vs free amino acids (2 h, 130 C)

| compound (LRI) | Gly (2 mmol) | GlyGly | Ala + Gly | AlaGly | Val + Gly | ValGly | Leu + Gly | LeuGly | Ser + Gly | SerGly | Pro + Gly | ProGly | rec |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| pyrazine (759) | | | | | | | | | 0.06 | | | | 0.3 |
| methylpyrazine (822) | | 0.02 | 0.01 | 0.01 | | | | | 0.09 | 0.02 | | | 0.8 |
| 2,5(6)-dimethylpyrazine (908) | 0.67 | 22.72 | 1.03 | 20.91 | 1.03 | 1.98 | 4.66 | 37.09 | 1.30 | 35.30 | 0.90 | | 2.0/1.6 |
| 2,3-dimethylpyrazine (913) | 0.05 | | 0.05 | | | | 0.08 | | 0.06 | | | | 1.6 |
| trimethylpyrazine (998) | 2.94 | 23.01 | 4.12 | 6.25 | 5.20 | 0.23 | 15.13 | 7.13 | 2.83 | 51.02 | 4.72 | | 4.1 |
| methyl-(1-methylethyl)pyrazine (1050) | | | | | | | 0.19 | | | | | | |
| 3-ethyl-2,5-dimethylpyrazine (1072) | 0.01 | 0.43 | 0.34 | 0.34 | 0.03 | | 0.54 | 0.50 | 0.15 | 0.17 | 0.03 | | |
| 2-ethyl-3,5-dimethylpyrazine (1080) | 0.09 | 0.23 | 0.09 | 0.04 | 0.14 | | 0.38 | 0.05 | 0.09 | 0.45 | 0.09 | | |
| tetramethylpyrazine (1082) | 0.38 | 0.66 | 0.10 | | 0.17 | | 1.29 | 0.01 | 0.33 | 0.43 | 0.39 | | 8.4 |
| 5-ethyl-2,3-dimethylpyrazine (1082) | 0.10 | | 0.37 | 0.06 | 0.72 | | 0.22 | 0.08 | | 0.66 | | | |
| 3-ethenyl-2,5-dimethylpyrazine (1093) | 0.03 | 0.04 | 0.01 | 0.02 | | | | 0.07 | 0.08 | 0.05 | 0.02 | | |
| 2-acetyl-5-methylpyrazine (1123) | 0.03 | 0.06 | 0.04 | 0.04 | 0.03 | | 0.15 | 0.09 | 0.08 | 0.11 | 0.09 | | |
| 2-acetyl-6-methylpyrazine (1133) | | | 0.01 | 0.03 | 0.02 | | 0.08 | 0.06 | 0.01 | | 0.03 | | |
| 5-isopropyl-2,3-dimethylpyrazine (1140) | | | | | 0.03 | | 1.00 | | | | | | |
| 2,3-diethyl-5-methylpyrazine (1165) | | | 0.01 | | | | | | 0.01 | | | | 30.0 |
| 3,5-diethyl-2-methylpyrazine (1167) | | | 0.03 | | | | | | 0.01 | | | | |
| 2,3,5-trimethyl-6-ethylpyrazine (1169) | | | 0.10 | | 0.07 | | 0.12 | | | | | | |
| 2,5-diethyl-3-methylpyrazine (1176) | 0.03 | | 0.01 | | 0.04 | | 0.04 | | | | | | |
| 2,5-dimethyl-3-(2-methylpropyl)pyrazine (1197) | | | | | 0.04 | | 0.04 | 0.15 | | | | | |
| 3,5-dimethyl-2-(2-methylpropyl)pyrazine (1211) | | | | | | | 0.16 | | | | | | |
| acetyldimethylpyrazine (1219, tent.) | 0.03 | 0.03 | 0.07 | 0.03 | 0.05 | | 0.24 | 0.05 | 0.08 | 0.25 | 0.04 | | |
| 2,3-dimethyl-5-(2-methylpropyl)pyrazine (1220) | | | | | | | 0.33 | | | | | | |
| 2,5-dimethyl-3-(E-1-propenyl)pyrazine (1229) | | 0.10 | | 0.12 | | | | 0.12 | | 0.08 | 0.01 | | |
| 2,3-dimethyl-5-(E-1-propenyl)pyrazine (1241) | | 0.04 | | 0.10 | | | 0.02 | 0.13 | | 0.03 | | | |
| 2-isopropenyl-dimethylpyrazine (1244, tent.) | | | | 0.13 | | | 0.02 | 0.26 | | 0.05 | | | |
| methyl-(3-methylbutyl)pyrazine (1250, tent.) | | | | | | | 0.03 | | | | | | |
| methyl-(3-methylbutyl)pyrazine (1253, tent.) | | | | | | | 0.10 | 0.01 | | | | | |
| (2-methylpropyl)trimethylpyrazine (1276) | | | | | 0.02 | | 0.04 | | | | | | |
| dimethyl-(3-methylbutyl)pyrazine (1315, tent.) | | | | | | | 5.15 | 1.46 | | | | | |
| dimethyl-(3-methylbutyl)pyrazine (1327, tent.) | | | | | | | 0.07 | | | | | | |
| trimethyl-(2-methylbutyl)pyrazine (1383, tent.) | | | | | | | 0.24 | | | | | | |
| trimethyl-(3-methylbutyl)pyrazine (1387, tent.) | | | | | | | 1.39 | | | | | | |
| **total pyrazines** | **4.36** | **47.33** | **6.39** | **28.09** | **7.61** | **2.22** | **31.72** | **47.26** | **5.18** | **88.63** | **6.32** | **0.00** | |
| pyrazines, % of total peak area | 58.2 | 77.8 | 63.5 | 61.8 | 45.0 | 31.2 | 21.7 | 63.6 | 44.8 | 89.3 | 35.2 | 0.0 | |

### Table 6. Glyoxal + X-Gly dipeptides vs free amino acids (2 h, 130 C)

| compound (LRI) | Gly (2 mmol) | GlyGly | Ala + Gly | AlaGly | Val + Gly | ValGly | Leu + Gly | LeuGly | Ser + Gly | SerGly | Pro + Gly | ProGly | rec |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| pyrazine (759) | 1.80 | 0.83 | 3.06 | 0.54 | 1.87 | | 3.39 | 0.07 | 2.68 | 0.14 | 1.62 | | 0.30 |
| methylpyrazine (822) | 0.16 | 0.16 | 0.11 | 0.10 | 0.06 | | 0.13 | 0.04 | 0.12 | 0.10 | 0.08 | | 0.80 |
| 2,5(6)-dimethylpyrazine (908) | 0.05 | 0.04 | 0.03 | 0.05 | 0.04 | | | 0.06 | 0.04 | | | | 2.0/1.6 |
| 2-ethylpyrazine (911) | | | 0.21 | | | | | | | | | | 2.30 |
| 2,3-dimethylpyrazine (913) | 0.10 | 0.02 | 0.10 | 0.03 | 0.03 | | | 0.01 | 0.06 | | | | 1.60 |
| 2-(1-methylethyl)pyrazine (977, tent.) | | | | | 0.03 | | 0.03 | | | | | | |
| 2-ethyl-6-methylpyrazine (995) | | | 0.06 | | | | | | | | | | |
| trimethylpyrazine (998) | 0.18 | | | | 0.09 | | 0.04 | | 0.09 | | | | 4.10 |
| 2-ethyl-5-methylpyrazine (998) | | | 0.14 | 0.01 | | | | | | | | | |
| 2-acetylpyrazine (1018) | | | 0.12 | | 0.02 | | 0.02 | | 0.02 | | 0.02 | | 0.80 |
| methyl-(1-methylethyl)pyrazine (1050 / 1057) | | | | | 0.06 / 0.03 | | | | | | | | |
| 2-(2-methylpropyl)pyrazine (1064, tent.) | | | | | 0.28 | | 0.04 | | | | | | |
| 3-ethyl-2,5-dimethylpyrazine (1072) | | | 0.05 | | | | | | 0.02 | | | | |
| 2-ethyl-3,5-dimethylpyrazine (1080) | | | 0.04 | | | | | | | | | | |
| 5-ethyl-2,3-dimethylpyrazine (1082) | | | 0.14 | | | | | | 0.02 | | | | |
| 2-acetyl-5-methylpyrazine (1123) | 0.04 | | | | | | | | | | | | |
| 2,3-dimethyl-5-(1-methylethyl)pyrazine (1144) | | | | | 0.05 | | 0.12 | | | | | | |
| methyl-(2-methylpropyl)pyrazines (1144 / 1146 / 1153) | | | | | 0.17 / 0.25 / 0.32 | | — / 0.05 / 0.02 | | | | | | |
| 2-(3-methylbutyl)pyrazine (1184, tent.) | | | | | | | 3.02 | 0.14 | | | | | |
| 2,5-dimethyl-3-(2-methylpropyl)pyrazine (1197) | | | | | 0.04 | | | | | | | | |
| 3,5-dimethyl-2-(2-methylpropyl)pyrazine (1211) | | | | | 0.16 | | | | | | | | |
| 2,3-dimethyl-5-(2-methylpropyl)pyrazine (1220) | | | | | 0.27 | | 0.01 | | | | | | |
| methyl-(3-methylbutyl)pyrazines (1253 / 1265, tent.) | | | | | | | 1.60 / 1.18 | 0.02 / 0.09 | | | | | |
| 1-(2-pyrazinyl)-3-methyl-1-butanone (1266, tent.) | | | | | 0.12 | | 0.19 | | | | | | |
| dimethyl-(3-methylbutyl)pyrazines (1315 / 1327 / 1337, tent.) | | | | | | | 0.34 / 0.39 / 0.83 | | | | | | |
| trimethyl-(2-methylbutyl)pyrazine (1383, tent.) | | | | | | | 0.25 | | | | | | |
| **total pyrazines** | **2.32** | **1.06** | **4.09** | **0.73** | **3.97** | **0.00** | **11.67** | **0.43** | **3.04** | **0.24** | **1.72** | **0.00** | |
| pyrazines, % of total peak area | 96.9 | 87.9 | 73.4 | 81.1 | 14.0 | 0.0 | 9.9 | 8.9 | 82.2 | 71.8 | 45.3 | 0.0 | |

### Table 7. Glucose + tripeptides vs free amino acids (2 h, 130 C)

| compound (LRI) | Gly (3 mmol) | GlyGlyGly | Lys + Gly + Gly | LysGlyGly | Lys + Ala + Pro | LysAlaPro | rec |
|---|---:|---:|---:|---:|---:|---:|---:|
| pyrazine (759) | 0.02 | | 0.47 | | 1.18 | | 0.3 |
| methylpyrazine (822) | 0.02 | 0.01 | 0.22 | 0.01 | 0.32 | 0.04 | 0.8 |
| 2,5(6)-dimethylpyrazine (908) | 0.06 | 0.74 | 0.31 | 0.10 | 0.20 | 0.80 | 2.0/1.6 |
| ethylpyrazine (911) | | | | | 0.13 | | 2.3 |
| 2,3-dimethylpyrazine (913) | 0.04 | 0.08 | 0.06 | | | 0.20 | 1.6 |
| 2-ethyl-6-methylpyrazine (995) | | 0.42 | 0.03 | | 0.01 | 0.05 | |
| trimethylpyrazine (998) | 0.07 | 0.10 | 0.16 | | 0.04 | 0.07 | 4.1 |
| 2-ethyl-5-methylpyrazine (998) | | | | 0.03 | | | |
| 2-ethenyl-5-methylpyrazine (1018) | | 0.09 | | | | 0.06 | |
| 3-ethyl-2,5-dimethylpyrazine (1072) | 0.01 | 0.01 | 0.04 | | 0.56 | 0.03 | |
| 2-ethyl-3,5-dimethylpyrazine (1080) | | 0.02 | | | | | |
| tetramethylpyrazine (1082) | 0.14 | | 0.04 | | | | 8.4 |
| 5-ethyl-2,3-dimethylpyrazine (1082) | | 0.02 | 0.04 | | | | |
| 2,3-diethyl-5-methylpyrazine (1165) | | | | | 0.02 | | 30.0 |
| 2,3,5-trimethyl-6-ethylpyrazine (1169) | 0.02 | | 0.03 | | 0.12 | | |
| **total pyrazines** | **0.38** | **1.48** | **1.39** | **0.13** | **2.58** | **1.24** | |
| pyrazines, % of total peak area | 30.8 | 60.2 | 38.8 | 7.5 | 11.8 | 37.2 | |

### Table 8. Methylglyoxal + tripeptides vs free amino acids (2 h, 130 C)

| compound (LRI) | Gly (3 mmol) | GlyGlyGly | Lys + Gly + Gly | LysGlyGly | Lys + Ala + Pro | LysAlaPro | rec |
|---|---:|---:|---:|---:|---:|---:|---:|
| methylpyrazine (822) | | | tr | | tr | tr | 0.8 |
| 2,5(6)-dimethylpyrazine (908) | 0.83 | 8.75 | 1.42 | 2.60 | 2.07 | 5.38 | 2.0/1.6 |
| 2-ethyl-5-methylpyrazine (998) | | 0.11 | tr | tr | | 0.02 | |
| trimethylpyrazine (998) | 6.18 | 5.78 | 3.45 | 3.40 | 0.46 | 3.64 | 4.1 |
| 3-ethyl-2,5-dimethylpyrazine (1072) | 0.03 | 0.25 | 0.34 | 0.04 | 1.36 | 0.30 | |
| 2-ethyl-3,5-dimethylpyrazine (1080) | 0.20 | 0.04 | 0.12 | 0.02 | tr | 0.04 | |
| tetramethylpyrazine (1082) | 1.01 | 0.12 | 0.22 | 0.04 | | 0.05 | 8.4 |
| 2-acetyl-5-methylpyrazine (1123) | 0.09 | tr | 0.02 | 0.33 | tr | 0.24 | |
| 2-acetyl-6-methylpyrazine (1133) | 0.05 | | 0.02 | | | 0.06 | |
| acetyldimethylpyrazine (1219, tent.) | 0.06 | 0.01 | 0.05 | 0.23 | 0.01 | 0.37 | |
| 2-acetyl-3,5,6-trimethylpyrazine (1269, tent.) | 0.02 | | tr | tr | | | |
| **total pyrazines** | **8.47** | **15.05** | **5.64** | **6.67** | **3.90** | **10.09** | |
| pyrazines, % of total peak area | 66.9 | 87.2 | 54.8 | 62.8 | 27.9 | 51.3 | |

### Table 9. Glyoxal + tripeptides vs free amino acids (2 h, 130 C)

| compound (LRI) | Gly (3 mmol) | GlyGlyGly | Lys + Gly + Gly | LysGlyGly | Lys + Ala + Pro | LysAlaPro | rec |
|---|---:|---:|---:|---:|---:|---:|---:|
| pyrazine (759) | 2.39 | 0.05 | 2.76 | 0.13 | 4.88 | 1.33 | 0.3 |
| methylpyrazine (822) | 0.11 | 0.02 | 0.13 | | 0.08 | 0.03 | 0.8 |
| 2,5(6)-dimethylpyrazine (908) | 0.10 | | 0.02 | | | | 2.0/1.6 |
| 2-ethylpyrazine (911) | | | | | 0.75 | | 2.3 |
| 2,3-dimethylpyrazine (913) | 0.10 | | 0.01 | | | | 1.6 |
| trimethylpyrazine (998) | 0.83 | | 0.01 | | | | 4.1 |
| acetylpyrazine (1018) | | | | | 0.03 | | |
| 2,6-diethylpyrazine (1072) | | | | | 0.03 | | |
| 2-ethyl-3,5-dimethylpyrazine (1080) | 0.01 | | 0.01 | | | | |
| tetramethylpyrazine (1082) | 0.04 | | | | | | 8.4 |
| 5-ethyl-2,3-dimethylpyrazine (1082) | 0.08 | | | | | | |
| 3-ethenyl-2,5-dimethylpyrazine (1093) | 0.01 | | | | | | |
| 2-acetyl-5-methylpyrazine (1123) | 0.02 | | | | | | |
| **total pyrazines** | **3.70** | **0.07** | **2.95** | **0.13** | **5.77** | **1.36** | |
| pyrazines, % of total peak area | 98.4 | 38.7 | 84.7 | 17.0 | 68.5 | 92.1 | |

### Table 10. Volatiles (SPME peak area x 1e6) from glucose + glycine (2 mmol) or diglycine (1 mmol), 2 h, pH 8, at four temperatures

| compound (LRI) | 100 C Gly | 100 C GlyGly | 130 C Gly | 130 C GlyGly | 150 C Gly | 150 C GlyGly | 180 C Gly | 180 C GlyGly |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 2-methylfuran | | | | | 6.8 | 29.2 | 10.0 | 99.8 |
| 2-ethylfuran | | | | | | 8.5 | 1.3 | 4.1 |
| 2,5-dimethylfuran | | | | | 2.0 | 10.1 | 8.8 | 26.1 |
| 2-vinylfuran | | | | | 3.2 | 1.8 | 9.9 | 6.3 |
| pyrazine (754) | | | 0.9 | | 8.3 | tr | 40.9 | 2.3 |
| dihydro-2-methyl-3(2H)-furanone (798) | | | | | 4.0 | 3.1 | 3.0 | 7.7 |
| 2-methylpyrazine (817) | | | | | 2.7 | tr | 60.8 | 8.5 |
| furfural (822) | | | | | 20.6 | 19.6 | 25.4 | 102.9 |
| 2-acetylfuran (904) | | | 0.1 | | 7.7 | 3.4 | 32.4 | 18.0 |
| 2,5(6)-dimethylpyrazine (906) | | | 0.7 | 15.2 | 14.4 | 16.6 | 210.6 | 54.2 |
| 2,3-dimethylpyrazine (911) | | | | | | | 53.2 | |
| 5-methylfurfural (961) | | | | | 5.5 | 6.3 | 161.3 | 290.9 |
| 2-ethyl-6-methylpyrazine (992) | | | | | | | 20.5 | |
| trimethylpyrazine (998) | | | 0.4 | 1.3 | 10.2 | 1.8 | 405.9 | 9.9 |
| 2-ethyl-5-methylpyrazine (998) | | | | 6.9 | 1.9 | 4.7 | 113.4 | 17.4 |
| 3-ethyl-2,5-dimethylpyrazine (1075) | | | | | | | 22.4 | tr |
| 2-ethyl-3,5-dimethylpyrazine (1080) | | | | | 0.8 | | 73.8 | tr |
| 5-ethyl-2,3-dimethylpyrazine (1082) | | | | | | | 157.5 | tr |
| 2,5-diethylpyrazine (1089) | | | | 0.5 | 1.8 | tr | 33.6 | tr |
| 2-(2-furanylmethyl)-5-methylfuran (1188) | | | | | | | | 17.4 |

At 100 C nothing was detected in either arm.

## 4. Numbers and steps the repository can use

Registry keys: `methylpyrazine`, `2_5_dimethylpyrazine` (and `2_6_dimethylpyrazine`; the paper cannot
separate them), `2_3_dimethylpyrazine`, `2_ethylpyrazine`, `trimethylpyrazine`, `tetramethylpyrazine`,
`2_ethyl_3_5_dimethylpyrazine` (SMILES caveat: see `cerny1994_extraction.md` flag 8), `furfural`,
`2_acetylfuran` exist. Pyrazine itself: family class `pyrazines` only, no molecule id. 2-Ethyl-5(6)-
methyl-, 3-ethyl-2,5-dimethyl-, 2,3-diethyl-5-methylpyrazine, the peptides, glyoxal, methylglyoxal:
**not in registry**. Every value below is a peak area; the ratios are ratios of peak areas of the SAME
compound between arms and are not converted to amounts.

| quantity or step | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| 2,5(6)-DMP from MGO: dipeptide / free amino acids | GlyGly / Gly 34; AlaGly / (Ala + Gly) 20; SerGly / (Ser + Gly) 27; LeuGly / (Leu + Gly) 8.0; ValGly / (Val + Gly) 1.9; ProGly 0 / 0.90; GlyLys / (Gly + Lys) 3.8; AlaLys / (Ala + Lys) 5.1; ValLys / (Val + Lys) 1.4 | ratio of peak areas | 1 mmol peptide vs 1 mmol each amino acid; MGO; water, pH 8, 130 C, 2 h | Tables 2, 5 | peak_area_only (within-study, same compound) |
| trimethylpyrazine from MGO: dipeptide / free | GlyGly 7.8; AlaGly 1.5; SerGly 18; LeuGly 0.47; ValGly 0.044; ProGly 0; GlyLys 1.5; AlaLys 4.4; ValLys 0.35 | ratio of peak areas | same | Tables 2, 5 | peak_area_only |
| 2,5(6)-DMP from glucose: dipeptide / free | GlyGly 22; AlaGly 5.8; SerGly 40; LeuGly 1.1; ValGly 0 / 0.16; ProGly 0 / 0.03; GlyLys 3.1; AlaLys 5.4; ValLys 0.17 | ratio of peak areas | glucose, same | Tables 1, 4 | peak_area_only |
| tripeptides with MGO, 2,5(6)-DMP peptide / free | GlyGlyGly / Gly(3) 10.5; LysGlyGly 1.8; LysAlaPro 2.6 | ratio | | Table 8 | peak_area_only (smaller advantage than dipeptides) |
| tripeptides with glucose, total pyrazines peptide / free | GlyGlyGly 3.9; LysGlyGly 0.09; LysAlaPro 0.48 | ratio | | Table 7 | peak_area_only (no consistent direction) |
| unsubstituted pyrazine from glyoxal: peptide / free | GlyGly 0.46; AlaGly 0.18; LeuGly 0.02; SerGly 0.05; ValGly 0; ProGly 0; GlyLys 0.41; AlaLys 0.73; ValLys 0.10; GlyGlyGly 0.02; LysGlyGly 0.05; LysAlaPro 0.27 | ratio of peak areas | glyoxal | Tables 3, 6, 9 | peak_area_only (always < 1: the parent pyrazine wants the free amino acid) |
| amino-acid-specific pyrazines, peptide / free | 3-ethyl-2,5-DMP from alanine: glucose AlaGly / (Ala + Gly) 0.05 (0.05 vs 0.93), AlaLys / (Ala + Lys) 0.14; MGO AlaGly 1.0 (0.34 vs 0.34), AlaLys 0.60; 2,5-dimethyl-3-isobutyl from valine: glucose 0 vs 0.72; 3,5-dimethyl-2-isobutyl, MGO: ValLys 0.01 vs 1.46; dimethyl-(3-methylbutyl) from leucine: glucose 0.02 vs 4.40, MGO 1.46 vs 5.15 | areas | | Tables 1, 2, 4, 5 | peak_area_only (Strecker-aldehyde adducts need the free amino acid; the one exception is Ala + Gly / AlaGly with MGO, equal) |
| pyrazines' share of all volatiles | dipeptides 80-90 % (GlyGly 87, AlaGly 84, SerGly 90, GlyLys 87, AlaLys 80) vs free amino acids 26-61 % | % of total peak area | glucose | Tables 1, 4 | peak_area_only (the peptide pot is a "cleaner" pyrazine pot) |
| N-terminal proline | ProGly: 0.00 total pyrazines with glucose, MGO and glyoxal; Pro + Gly gives 0.18 / 6.32 / 1.72 and proline-specific pyrrolizines | areas | | Tables 4-6 | peak_area_only (non-detect; no peptide hydrolysis) |
| temperature crossover, glucose + Gly (2 mmol) vs GlyGly (1 mmol) | 2,5(6)-DMP 0.7 vs 15.2 (130 C); 14.4 vs 16.6 (150 C); 210.6 vs 54.2 (180 C); trimethyl 0.4 vs 1.3; 10.2 vs 1.8; 405.9 vs 9.9; nothing at 100 C | SPME area x 1e6 | 2 h, pH 8 | Table 10 | peak_area_only ladder (the dipeptide's advantage holds at 130 C, is gone by 150 C, and inverts 4-40x at 180 C; furans and furfurals take over for the dipeptide) |
| trimethyl : 2,5(6)-dimethyl with MGO, free glycine | 4.4 (2 mmol Gly); 7.4 (3 mmol Gly) vs 1.0 (GlyGly), 0.66 (GlyGlyGly) | ratio of areas of DIFFERENT compounds — recovery 4.1 vs 2.0 % | | Tables 5, 8 | peak_area_only; cross-compound, orientation only: free glycine + MGO makes proportionally more of the formaldehyde adduct than diglycine does |

**The mechanism the authors give for peptides (Scheme 1).** Peptide N-terminus + alpha-dicarbonyl ->
imine; deprotonation and a 1,5-H shift -> 4-hydroxy-2-azadiene (the carbonyl of the future
aminoketone enolised); hydrolysis of the imino group -> the alpha-aminoketone + an alpha-keto-acyl
peptide. No decarboxylation, no Strecker aldehyde. The aminoketone is the SAME one an amino acid would
give (aminoacetone from MGO, aminoacetaldehyde from glyoxal), so the pyrazine set is the same minus the
Strecker-aldehyde adducts. Why dipeptides beat free amino acids with glucose: intramolecular
protonation of the glycosylamine imine by the C-terminal carboxylate (de Kok & Rosing 1994) catalyses
the Amadori rearrangement; tripeptides and free amino acids lack the geometry. Why Val/Leu/Pro
N-termini are poor: steric hindrance of the alpha-deprotonation (Val, Leu); secondary amine and no
Amadori catalysis (Pro).

**What this means for a protein-isolate amine (roadmap §5d).** A rule "peptide N-terminus + dicarbonyl ->
aminoketone" with the same products as R07/R28 but NO Strecker aldehyde, weighted by the N-terminal
residue (Gly/Ala/Ser/Lys high; Val/Leu low; Pro zero), would reproduce this paper's directions. The
absolute rate is unmeasurable here (peak areas, unstated concentrations). The unsubstituted-pyrazine
deficit with peptides is unexplained by the authors and would be a hold-out shape.

## 5. Flags

1. **Peak areas only, no internal standard, no calibration, no replicates.** The recovery column shows
   cross-compound areas differ by up to 100x in what they represent; only same-compound comparisons
   across arms are meaningful, and even those assume the matrix (peptide vs amino acid solution, pH
   drift) does not change the SBSE partition.
2. **Concentrations are not in this paper**: glucose and dicarbonyl amounts and the water volume are by
   reference to Van Lancker 2010. Do not assign mmol/L until that paper is on disk.
3. **Unequal amine charge**: 1 mmol dipeptide vs 1 + 1 mmol free amino acids (2 mmol alpha-amino
   groups); the dipeptide "advantage" ratios above are therefore conservative by up to 2x on a
   per-amino-group basis.
4. **Unbuffered, initial pH 8**; pH after 2 h at 130 C not reported; peptide and amino-acid arms may
   drift differently.
5. **2,5- and 2,6-dimethylpyrazine co-elute** on HP5-MS; the registry's two keys share one number.
6. **Glucose and dicarbonyl arms are not comparable** (10x lower dicarbonyl concentration).
7. **Table 10 used SPME**, the rest SBSE; the authors say the 130 C values differ between Tables 4 and 10.
8. **Many identifications are tentative** (marked); the isobutyl/isopentyl pyrazines rest on MS only.
9. The "theor recovery" values are the authors' theoretical PDMS partition numbers, not measured
   recoveries; treat them as a warning, not a correction factor.
