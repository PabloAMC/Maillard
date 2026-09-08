# Zhang et al. 2020b — EXTRACTION (odour-active volatiles of raw pea milk vs soy milk, HS-SPME GC-O-MS with OAVs; LOX / HPL / ADH / lipase activities; free fatty acids)
### The "beany note before any heat" measured at its source: soaked seeds blended in water at room temperature, no heating, LOX pathway enzymes assayed in the same milk.

**Source on disk:** `data/articles/zhang2020b.pdf` (owner's download, 2026-09-08; Journal Pre-proof,
Food Chemistry, FOCH 127469). Read from the scratchpad text layer (`zhang2020b.txt`); all six
tables (1A, 1B, 2A, 2B, 2C, 3, 4A, 4B) are clean in the text layer and are re-typed below; Figs 1-2
are FIGURE-ONLY. NOT the same paper as `zhang2020_extraction.md` (Zhang et al., alpha-dicarbonyls in
glucose-glutamate; different authors, different journal); the "b" suffix separates them. Repo status
before this dossier: the roadmap (`tasks/roadmap_for_scientists.md` §5) names "hexanal from
lipoxygenase during processing" as needing "its own module and its own data programme";
`results/validation/matrix_sites_prereg.md` §4 says the pea/soy hexanal misses "are a storage and
lipoxygenase question". No LOX activity, no pea free-fatty-acid datum and no raw-pea-milk hexanal
level existed in `data/lit/` before this dossier.

## 0. Identity

| field | value as printed |
|---|---|
| Title | "Key volatile off-flavor compounds in peas (Pisum sativum L.) and their relations with the endogenous precursors and enzymes using soybean (Glycine max) as a reference" |
| Authors | Caimeng Zhang, Yufei Hua (corresponding), Xingfei Li, Xiangzhen Kong, Yeming Chen — Jiangnan University, Wuxi |
| Venue | Food Chemistry, journal pre-proof; PII S0308-8146(20)31331-5; received 30 Mar 2020, revised 15 Jun 2020, accepted 29 Jun 2020. Volume/pages of the version of record are not printed in this PDF. |
| DOI | 10.1016/j.foodchem.2020.127469 |
| Naming | "pea milk" / "soy milk" = the 3000 g supernatant of soaked seeds blended in water at room temperature, diluted to 2 % protein; "cooked" milk = the same from seeds boiled 50 min before blending. LOX-1/-2/-3 are the soybean-style isoform names assigned by assay pH (9.0 / 6.8 / 7.1). "Relative quantities" (Table 2B) and "content" (Table 2C) are two different quantifications of the same peaks (§2). |
| Compound registry | hexanal -> `hexanal`; nonanal -> `nonanal`; 1-hexanol -> `1_hexanol`; 1-octen-3-ol -> `1_octen_3_ol`; 2-pentylfuran -> `2_pentylfuran`; (E)-2-octenal -> `e_2_octenal`; (E,E)-2,4-decadienal -> **not in registry** (registry has only `e_e_2_4_heptadienal`); (E,E)-2,4-nonadienal, (E)-2-hexenal, octanal, decanal, 1-octen-3-one, 1-octanol, 1-nonanol, (E)-2-decenal, dihydro-5-pentyl-2(3H)-furanone -> not in registry; 2-methoxy-3-isopropyl-(5 or 6)-methylpyrazine -> not in registry as a species (registry has `3_isobutyl_2_methoxypyrazine` and the group `methoxypyrazines`). |

## 1. Why it matters

Need (b), directly. This is a measurement of the LOX-pathway volatiles formed when raw yellow pea is
wet-processed with no heat at all (soak 12 h at 4 C, blend 3 min at room temperature, filter,
centrifuge), quantified against authentic standards, alongside the free linoleic acid available and
the specific activities of the enzymes that make hexanal from it (LOX-2, LOX-3, HPL) and that turn
hexanal into 1-hexanol (ADH), with soybean as the reference run through the identical protocol. It
gives the repository: (i) an end-of-process hexanal level for raw pea milk (164.18 ug/L at 2 %
protein) and for soy milk (437.39 ug/L), validation-only; (ii) the within-study soy/pea ratios the
authors use to argue that enzyme activity, not lipid content, sets the hexanal level (hexanal +
1-hexanol soy/pea = 1.44 while free linoleic acid soy/pea = 4.1); (iii) a molar yield of C6 products
on free linoleic acid (~7-11 % in pea milk, §4), which is the kind of number a LOX module needs to
bound its conversion; (iv) the precursor pool (free linoleic acid ~14,300 ug/L in pea milk, fatty
acid profile of the FFA, 38.94 % C18:2); (v) LOX-2 / LOX-3 / HPL / ADH specific activities in the
pea milk in stated units; (vi) OAVs and thresholds for the compounds the roadmap lists (hexanal,
nonanal, 2-pentylfuran, 1-hexanol). It gives NO time course: one process, one time point. Need (a)
is touched only in that the milk is a 2 % pea-protein dispersion, so every level is already net of
whatever binding the pea protein exerts at 50 C in the SPME vial.

## 2. Methods as they matter to a model

- **Seeds.** Yellow pea, Saskatchewan 2018 harvest; soybean, Heilongjiang 2018; stored 4 C; work
  done 2019. Seed composition Table 1A (pea: 19.82 % protein, 2.17 % lipid, 88.62 % dry matter).
- **Milk preparation (verbatim core):** "One hundred grams of pea and soybean seeds were rinsed and
  soaked in 300 mL distilled water at 4 C for 12 h. ... Distilled water was added in the soaked or
  cooked seeds to a total weight of 800 g and the whole mixture was ground at room temperature using
  a blender ... at medium speed for 3 min. The resulting slurries were filtered through two layers
  of muslin, and the filtrates were centrifuged at 3,000 g for 15 min at 4 C. The obtained
  supernatants were designated as pea and soy milk". So **100 g seed per 800 g slurry (12.5 % w/w),
  no heat, ~3 min of active blending plus filtration and 15 min centrifugation before the
  sample exists**; no enzyme-inactivation step. Cooked control: 100 g unsoaked seeds in 2000 mL
  boiling water 50 min, then the same. "Before analysis, the pea and soy milk were uniformly diluted
  to 2% protein content (Table 1A)." Pea milk as analysed: 2.00 % protein, 3.51 % dry matter, 0.19 %
  lipid; soy milk: 2.00 % protein, 4.45 % dry matter, 1.06 % lipid.
- **HS-SPME.** 5 mL milk + 0.5 g NaCl in a 15-mL vial; internal standard 1 uL of 2-methyl-3-heptanone
  at 0.25 ug/uL in methanol (= 0.25 ug per 5 mL = **50 ug/L**); DVB/CAR/PDMS 50/30 um fibre;
  **50 C, 10 min equilibration, 30 min extraction**; desorption 250 C, 3 min. GC-O-MS: Bruker QP 460,
  DB-WAX 30 m x 0.25 mm x 0.25 um; 40 C -> 85 C at 45 C/min -> 200 C at 9 C/min -> 250 C at
  45 C/min, hold 3 min; He 2 mL/min; 1:1 split to MS and sniff port; EI 70 eV, m/z 40-350.
  Identification by NIST 2.0, LRI (C7-C27 alkanes), odour, and authentic standards (column "S").
- **Two quantifications.** (1) Table 2A/2B, verbatim footnote: "Relative quantities were calculated
  using the internal standard 2-methyl-3-heptanone." (single-IS semi-quantification, response
  factor 1). (2) Table 2C, verbatim: "The standard curves of volatile compounds identified by
  GC-O-MS were established using 2% fat cow's milk as the matrix (Yuan & Chang, 2007), and were used
  for identification and quantification of the key volatile odor compounds. A series of
  concentrations of the authentic flavor standards with added 2-methyl-3-heptanone were prepared.
  The extraction and analysis conditions were same to those of volatiles." Footnote b of Table 2C:
  "x is the peak area relative to that of the internal standard 2-methyl-3-heptanone, and y is the
  concentration (ug/L) in the sample relative to that of the internal standard." The calibration
  matrix is 2 % fat cow's milk, not pea milk (0.19 % lipid) — flag 2. The two hexanal numbers for
  pea milk are 140.41 (2B) and 164.18 ug/L (2C).
- **OAV** = content (Table 2C) / odour threshold in water (van Gemert 2011); contribution rate =
  OAV / sum of OAVs of the 13 off-flavour compounds. Sensory: 10 trained panellists, 0-5 scale, six
  references (beany, oil-oxidation, sweet, earthy, grain husk, mushroom); Fig. 1 spider plot,
  FIGURE-ONLY.
- **Free fatty acids.** TLC isolation, FAME, GC on DB-WAX (Yoshida 2007). FFA = 1.75 g/100 g lipid
  in pea seed, 1.56 in soybean; composition Table 1B.
- **LOX assay (verbatim core):** "the protein contents of the pea and soy milk were adjusted to
  0.5% with M/15 phosphate buffer (pH 6.8, 7.1) or M/15 borate buffer (pH 9.0). Twenty milliliter of
  0.5% protein content sample was taken and heated at 100 C for 20 min as blank. The substrate
  solution was prepared by mixing 157.2 uL of linoleic acid, 157.2 uL of Tween-20, and 10 mL of
  deionized water. The solution was then clarified by adding of 1 mL 1 M NaOH, and diluted to
  200 mL M/15 phosphate buffer (pH 6.8, 7.1) and 200 mL M/15 borate buffer (pH 9.0), respectively.
  The substrate solution (29 mL) was stirred for 2 min at 25 C, and the reaction was started by
  adding 1 mL of milk. One unit of LOX activity was defined as an increase in absorbance of 0.001 at
  234 nm (LOX-1, pH 9.0; LOX-2, pH 6.8) or 280 nm (LOX-3, pH 7.1) per minute per milligram of
  protein." Arithmetic: 157.2 uL linoleic acid x 0.902 g/mL (handbook density) = 141.8 mg =
  0.506 mmol; if the whole stock is brought to 200 mL, substrate = **2.5 mM**, and in the 30-mL
  reaction (29 + 1 mL) **~2.4 mM linoleic acid with 5 mg milk protein at 25 C**. Whether one 200-mL
  batch per buffer was made from one stock each, or one stock was split, is not stated (flag 5).
  The enzyme source is the **milk itself** (not a separate extract) — "The crude enzyme samples
  were all extracted from soaked pea and soybean seeds." 1 U = 0.001 A/min/mg protein; with the
  conjugated-diene epsilon(234) of ~25,000 M^-1 cm^-1 (not printed) 0.001 A/min = 0.04 uM/min in a
  1-cm cell; how the 30-mL reaction was read (aliquots? path length?) is not stated, so a conversion
  to umol/min needs the reaction volume assumption (flag 6).
- **HPL assay.** Separate extract: 20 g soaked seeds in 100 mL Bis-Tris 0.1 M pH 6.8 with 0.5 %
  Triton X-100, 4 mM DTT, 0.5 % PVP, 25,000 g 45 min. Substrate 13-hydroperoxy-linoleic acid made
  with soybean LOX-1 (pH 9, 3 h, 4 C), 10 uL of 10 mM in 2.99 mL 0.1 M phosphate pH 6.5
  (= 33 uM); 10 uL extract; 25 C, 2 min; unit = loss of "1 uM" (sic; umol intended) of substrate at
  234 nm per min per mg protein.
- **ADH assay.** 20 g soaked seeds / 100 mL 0.1 M Na-phosphate pH 8.5; medium 0.1 M phosphate
  pH 8.5, 10 mM mercaptoethanol, 100 mM ethanol, 860 nM NAD (sic), 25 C, 2 min; unit = 0.001 A(340)
  per min per mg protein (ethanol -> acetaldehyde direction, i.e. the reverse of hexanal reduction).
- **Lipase assay.** pH-stat at pH 8.0 on 1.5 mM TAG emulsion; IU = 1 umol FFA/min; per mg protein.
  Not detectable in either seed.
- **LC-MS/MS proteomics.** FASP digestion of milk and oil-body proteins, Q Exactive, MaxQuant
  1.3.0.5 label-free (LFQ) intensities, Table 4A.
- **Pyrazine in dry seed.** Freeze-dried seed milled under liquid N2, HS-SPME at 50 C for
  increasing times; plateau after 2 h at 0.039 +/- 0.003 ug/g (text; the curve is Fig. 2,
  FIGURE-ONLY).
- **Statistics.** Triplicates, mean +/- SD, p <= 0.05 (SPSS 22).

## 3. Tables re-typed

### Table 1A. "Characterization of pea/soybean seeds and milk."

| Samples | Pea seeds | Pea milk | Soybean seeds | Soy milk |
|---|---|---|---|---|
| Dry matter (%) | 88.62 +/- 0.05 b | 3.51 +/- 0.02 a | 88.81 +/- 0.06 b | 4.45 +/- 0.02 a |
| Crude protein content (%) | 19.82 +/- 0.08 b | 2.00 +/- 0.01 a | 34.70 +/- 0.09 c | 2.00 +/- 0.03 a |
| Crude lipid content (%) | 2.17 +/- 0.05 c | 0.19 +/- 0.02 a | 16.75 +/- 0.17 d | 1.06 +/- 0.03 b |
| Total chlorophyll (%) | 0.0079 +/- 0.0005 a | - | 0.0028 +/- 0.0002 a | - |
| FFA (g/100 g lipid) | 1.75 +/- 0.11 a | - | 1.56 +/- 0.08 a | - |

"-" = not tested. Letters: row-wise significance (p < 0.05).

### Table 1B. "Free fatty acids compositions of pea and soybean seeds." (% of FFA)

| Fatty acid | Pea (%) | Soybean (%) |
|---|---|---|
| C16:0 | 19.56 +/- 1.06 b | 12.13 +/- 0.57 a |
| C18:0 | 8.22 +/- 0.37 a | 9.03 +/- 0.26 a |
| C18:1 | 11.37 +/- 0.78 a | 17.90 +/- 0.55 b |
| C18:2 | 38.94 +/- 1.12 a | 35.07 +/- 0.91 a |
| C18:3 | 16.87 +/- 0.61 a | 19.20 +/- 0.62 a |
| C20:0 | 2.64 +/- 0.11 a | 4.08 +/- 0.22 b |
| C20:1 | 1.76 +/- 0.09 a | 1.66 +/- 0.08 a |
| C22:0 | 0.65 +/- 0.03 a | 0.92 +/- 0.05 b |

Column sums: pea 100.01, soybean 99.99.

### Table 2A. "Quantitative comparison of different chemical classes in pea and soy milk." (relative quantities, ug/L, by IS)

| Class | Pea milk (ug/L) | Soy milk (ug/L) | n compounds pea | n compounds soy |
|---|---|---|---|---|
| Alcohols | 431.37 +/- 52.60 b | 269.40 +/- 20.49 a | 36 | 24 |
| Aldehydes | 197.38 +/- 20.85 a | 490.05 +/- 35.29 b | 22 | 21 |
| Ketones | 17.28 +/- 2.23 a | 27.37 +/- 3.26 b | 12 | 12 |
| Esters | 27.15 +/- 2.53 b | 9.77 +/- 0.85 a | 6 | 7 |
| Pyrazines | 1.21 +/- 0.19 | ND | 1 | ND |
| Furans | 7.03 +/- 0.70 a | 6.37 +/- 0.98 a | 1 | 2 |
| Others | 2.10 +/- 0.35 a | 1.53 +/- 0.31 a | 2 | 2 |
| Total | 683.52 +/- 79.45 | 804.49 +/- 61.18 | 80 | 68 |

Re-summed: pea classes total 683.52 (exact); soy 804.49 (exact). The header count "68" for soy
contradicts the text's "60 aroma compounds ... from soy milk" (flag 9). The full 80/68-compound list
is not in the paper (only classes and the 17 odour-active ones).

### Table 2B. "Odorants identified by GC-O-MS and relative quantities calculated for the aromas of pea and soy milk." (ug/L, by IS; odour class p/n/u; ID method)

| Compound | LRI (DB-WAX) | Pea milk | Soy milk | Odour (class) | ID |
|---|---|---|---|---|---|
| Hexanal | 1060 | 140.41 +/- 10.20 a | 374.06 +/- 23.19 b | Grassy, green (n) | MS, O, S |
| (E)-2-Hexenal | 1198 | 3.21 +/- 0.67 a | 17.22 +/- 1.59 b | Leaf (u) | MS, O, S |
| 2-Pentylfuran | 1206 | 7.03 +/- 0.70 b | 3.61 +/- 0.57 a | Green bean (n) | MS, O, S |
| Octanal | 1266 | 3.40 +/- 0.61 a | 4.20 +/- 0.33 b | Lemon, fruity (p) | MS, O |
| 1-Octen-3-one | 1276 | ND | 2.21 +/- 0.29 | Mushroom (u) | MS, O, S |
| 1-Hexanol | 1338 | 239.60 +/- 23.58 b | 175.37 +/- 11.32 a | Lemon, grass, green (n) | MS, O, S |
| Nonanal | 1362 | 14.00 +/- 3.74 a | 15.39 +/- 2.43 a | Plastic, citrus (u) | MS, O, S |
| (E)-2-Octenal | 1406 | 7.96 +/- 1.76 a | 9.62 +/- 0.98 b | Cucumber, vegetable (n) | MS, O, S |
| 1-Octen-3-ol | 1413 | 87.22 +/- 13.82 b | 32.14 +/- 2.34 a | Mushroom (u) | MS, O, S |
| 2-Methoxy-3-isopropyl-(5 or 6)-methyl pyrazine | 1417 | 1.21 +/- 0.19 | ND | Earthy, spicy, plastic, hay (u) | MS, O, S |
| Decanal | 1471 | 1.40 +/- 0.25 a | 1.36 +/- 0.24 a | Earthy, mushroom (u) | MS, O, S |
| 1-Octanol | 1528 | 24.62 +/- 5.06 b | 7.22 +/- 0.55 a | Oily, aldehydic (u) | MS, O, S |
| (E)-2-Decenal | 1612 | 5.75 +/- 1.12 a | 11.05 +/- 1.38 b | Orange (p) | MS, O |
| 1-Nonanol | 1627 | 18.29 +/- 3.05 b | 7.38 +/- 1.16 a | Rose-orange (p) | MS, O |
| (E,E)-2,4-Nonadienal | 1676 | 0.79 +/- 0.19 a | 1.41 +/- 0.19 b | Cucumber, green (n) | MS, O, S |
| (E,E)-2,4-Decadienal | 1746 | 0.29 +/- 0.06 a | 0.72 +/- 0.19 b | Fatty (u) | MS, O, S |
| Dihydro-5-pentyl-2(3H)-furanone | 2021 | ND | 0.44 +/- 0.12 | Peach, coconut (p) | MS, O |

Footnotes: LRI on C7-C27 alkanes; quantities by IS 2-methyl-3-heptanone; p/n/u = pleasant / neutral /
unpleasant at the sniff port; S = confirmed with authentic standard; ND = not detected.

### Table 2C. "Contribution rates of odor-active compounds with OAVs greater than 1 in pea and soy milk." (content by calibration curve, ug/L)

| Compound | Calibration equation (y = conc. ug/L rel. IS; x = rel. peak area) | R^2 | Threshold (ug/L in water) | Pea content (ug/L) | Pea OAV | Pea contribution (%) | Soy content (ug/L) | Soy OAV | Soy contribution (%) |
|---|---|---|---|---|---|---|---|---|---|
| Hexanal | y = (0.0000004) x + 18.60 | 0.9994 | 4.5 | 164.18 +/- 11.93 a | 36.48 +/- 2.65 | 10.73 +/- 0.78 | 437.39 +/- 27.11 b | 97.20 +/- 6.03 | 19.86 +/- 1.23 |
| (E)-2-Hexenal | y = (0.0000008) x + 1.05 | 0.9991 | 40 | 8.16 +/- 1.71 a | < 1 | - | 43.82 +/- 4.04 b | 1.10 +/- 0.11 | 0.22 +/- 0.02 |
| 2-Pentylfuran | y = (0.0000012) x + 3.62 | 0.9986 | 4.8 | 31.66 +/- 3.13 b | 6.60 +/- 0.65 | 1.94 +/- 0.19 | 16.26 +/- 2.57 a | 3.39 +/- 0.54 | 0.69 +/- 0.11 |
| 1-Octen-3-one | y = (0.0000014) x + 1.21 | 0.9968 | 0.007 | ND | -* | - | 1.64 +/- 0.21 | 233.71 +/- 30.43 | 47.76 +/- 6.22 |
| 1-Hexanol | y = (0.0000005) x + 21.35 | 0.9961 | 200 | 387.39 +/- 38.12 b | 1.94 +/- 0.19 | 0.57 +/- 0.06 | 283.54 +/- 18.30 a | 1.42 +/- 0.09 | 0.29 +/- 0.02 |
| Nonanal | y = (0.0000005) x + 2.25 | 0.9952 | 3.5 | 7.88 +/- 2.11 a | 2.25 +/- 0.60 | 0.66 +/- 0.18 | 8.66 +/- 1.37 a | 2.48 +/- 0.39 | 0.51 +/- 0.08 |
| (E)-2-Octenal | y = (0.0000007) x - 1.36 | 0.9992 | 4 | 9.84 +/- 2.17 a | 2.46 +/- 0.54 | 0.72 +/- 0.16 | 11.88 +/- 1.20 b | 2.97 +/- 0.30 | 0.61 +/- 0.06 |
| 1-Octen-3-ol | y = (0.0000004) x - 12.26 | 0.9988 | 7 | 105.10 +/- 16.65 b | 15.02 +/- 2.38 | 4.42 +/- 0.70 | 38.73 +/- 2.82 a | 5.53 +/- 0.40 | 1.13 +/- 0.08 |
| 2-Methoxy-3-isopropyl-(5 or 6)-methyl pyrazine | y = (0.0000007) x + 1.32 | 0.9985 | 0.02 | 4.13 +/- 0.65 | 205.65 +/- 32.65 | 60.46 +/- 9.60 | ND | - | - |
| Decanal | y = (0.0000006) x + 0.66 | 0.9896 | 0.9 | 1.82 +/- 0.32 a | 2.02 +/- 0.36 | 0.59 +/- 0.11 | 1.76 +/- 0.31 a | 1.96 +/- 0.68 | 0.40 +/- 0.14 |
| 1-Octanol | y = (0.0000008) x + 8.57 | 0.9927 | 54 | 61.16 +/- 12.57 b | 1.13 +/- 0.23 | 0.33 +/- 0.07 | 7.93 +/- 1.36 a | < 1 | - |
| (E,E)-2,4-Nonadienal | y = (0.000001) x + 1.88 | 0.9938 | 0.06 | 2.18 +/- 0.53 a | 36.35 +/- 8.85 | 10.69 +/- 2.60 | 3.87 +/- 0.52 b | 64.53 +/- 8.65 | 13.19 +/- 1.77 |
| (E,E)-2,4-Decadienal | y = (0.0000016) x + 1.69 | 0.9957 | 0.05 | 1.51 +/- 0.31 a | 30.22 +/- 6.26 | 8.89 +/- 1.84 | 3.75 +/- 0.99 b | 75.02 +/- 19.80 | 15.33 +/- 4.05 |

Footnotes: thresholds from van Gemert 2011; contents from calibration curves of authentic
standards; "* Not calculated"; ND = not detected. Arithmetic check: hexanal pea OAV 164.18 / 4.5 =
36.48 (exact); soy 437.39 / 4.5 = 97.20 (exact). Note the calibration intercepts are large for
hexanal (18.60) and 1-hexanol (21.35) relative to the pea contents.

### Table 3. "Effects of heat treatment on the contents of 1-octen-3-ol and 1-octen-3-one." (ug/L)

| Compound | Cooked pea milk | Pea milk | Cooked soy milk | Soy milk |
|---|---|---|---|---|
| 1-Octen-3-ol | 38.08 +/- 2.29 b | 105.10 +/- 16.65 c | 12.32 +/- 0.95 a | 38.73 +/- 2.82 b |
| 1-Octen-3-one | ND | ND | 0.17 +/- 0.01 a | 1.64 +/- 0.21 b |

Ratios: pea 1-octen-3-ol cooked/raw = 0.362 (text: "about 64%" decrease); soy 0.318 (text 68 %);
soy 1-octen-3-one 0.104 (text 90 %). **Hexanal in the cooked milks is not reported** (flag 1).

### Table 4A. "Identification of the endogenous enzymes using LC-MS/MS analysis."

| No. | Protein | Accession | Score | Peptides | Unique | Coverage (%) | LFQ intensity |
|---|---|---|---|---|---|---|---|
| 1 (pea milk) | Lipoxygenase-3 | P09918.1 | 323.31 | 35 | 32 | 39.7 | 72,801,000,000 |
| 2 (pea milk) | Lipoxygenase-2 | P14856.1 | 323.31 | 34 | 2 | 44.9 | 56,738,000,000 |
| 3 (pea milk) | Hydroperoxide lyase | AGH32771.1 | 133.65 | 2 | 2 | 16.5 | 1,833,000 |
| 4 (pea milk) | Alcohol dehydrogenase | P12886.1 | 148.49 | 14 | 14 | 35 | 21,132,000,000 |
| 5 (pea milk) | Phospholipase C | CAA75546.2 | 101.28 | 6 | 4 | 9.1 | 424,150,000 |
| 6 (soy milk) | Lipoxygenase-3 | NP_001235383.2 | 246.36 | 35 | 25 | 43.2 | 56,395,000,000 |
| 7 (soy milk) | Lipoxygenase-2 | NP_001237685.2 | 323.31 | 43 | 30 | 53.9 | 36,226,000,000 |
| 8 (soy milk) | Lipoxygenase-1 | NP_001236153.2 | 323.31 | 39 | 3 | 55.4 | 24,545,000,000 |
| 9 (soy milk) | Hydroperoxide lyase | AGH32771.1 | 126.36 | 2 | 2 | 16.1 | 4,003,000 |
| 10 (soy milk) | Alcohol dehydrogenase | XP_003526673.1 | 106.75 | 10 | 4 | 33.2 | 3,609,850,000 |
| 11 (soy milk) | Phospholipase A-2-activating protein | XP_003537897.1 | 89.43 | 7 | 7 | 9.6 | 468,020,000 |
| 12 (pea oil body) | Lipase | BAA85654.1 | 53.35 | 3 | 1 | 7.2 | 58,979,810 |
| 13 (soy oil body) | Lipase 3-like isoform X2 | XP_014631379.1 | 57.52 | 1 | 1 | 7.4 | 13,599,000 |

Note: pea LOX-2 rests on 2 unique peptides (vs 32 for LOX-3); HPL on 2 peptides in both species;
the pea HPL is matched to a soybean-database accession (AGH32771.1 appears for both).

### Table 4B. "Determination of LOX pathway enzyme activities under each optimal pH."

| Specific activity | Pea | Soybean |
|---|---|---|
| LOX-1 (U/mg protein), pH 9.0, A234 | ND | 14010 +/- 563 |
| LOX-2 (U/mg protein), pH 6.8, A234 | 2160 +/- 38 b | 1445 +/- 30 a |
| LOX-3 (U/mg protein), pH 7.1, A280 | 150 +/- 6 b | 122 +/- 4 a |
| HPL (umol/(min mg)), pH 6.5 | 0.15 +/- 0.03 a | 0.33 +/- 0.05 b |
| ADH (U/mg protein), pH 8.5, A340 | 178 +/- 12 b | 31 +/- 3 a |
| Lipase (IU/mg protein), pH 8.0 | ND | ND |

1 U (LOX, ADH) = 0.001 absorbance units per min per mg protein; assays at 25 C; "protein" = protein
of the milk (LOX) or of the seed extract (HPL, ADH). ND = not detected.

## 4. Numbers the repository can use

| quantity | value | unit | conditions | source location | evidence class | binding record fit |
|---|---|---|---|---|---|---|
| Hexanal, raw pea milk | 164.18 +/- 11.93 | ug/L | 2 % protein pea milk (12.5 % seed w/w slurry, soak 12 h / 4 C, blend 3 min RT, no heat); HS-SPME 50 C 30 min; calibrated in 2 % cow's milk | Table 2C | measured, **level_only** (end of process) | none (not a binding number). Validation target for a future LOX module; the repository must not fit it. |
| Hexanal, raw pea milk, IS-relative | 140.41 +/- 10.20 | ug/L | same | Table 2B | measured, level_only | second quantification of the same peak; 0.855x the calibrated value |
| Hexanal, raw soy milk | 437.39 +/- 27.11 (calibrated); 374.06 (relative) | ug/L | same protocol, soybean | Tables 2C / 2B | measured, level_only | none |
| 1-Hexanol, pea / soy milk | 387.39 +/- 38.12 / 283.54 +/- 18.30 (calibrated); 239.60 / 175.37 (relative) | ug/L | same | Tables 2C / 2B | measured, level_only | none; the ADH product of hexanal |
| Nonanal, pea / soy milk | 7.88 +/- 2.11 / 8.66 +/- 1.37 | ug/L | same | Table 2C | measured, level_only | none |
| 2-Pentylfuran, pea / soy milk | 31.66 +/- 3.13 / 16.26 +/- 2.57 | ug/L | same | Table 2C | measured, level_only | none |
| 1-Octen-3-ol, pea / soy milk (raw) | 105.10 +/- 16.65 / 38.73 +/- 2.82 | ug/L | same | Tables 2C, 3 | measured, level_only | none |
| 1-Octen-3-ol, cooked-seed pea / soy milk | 38.08 +/- 2.29 / 12.32 +/- 0.95 | ug/L | seeds boiled 50 min before blending | Table 3 | measured, level_only | none; enzymatic share of 1-octen-3-ol = 64 % (pea), 68 % (soy) by the authors' subtraction |
| (E,E)-2,4-Decadienal, pea / soy | 1.51 +/- 0.31 / 3.75 +/- 0.99 | ug/L | same | Table 2C | measured, level_only | none |
| (E,E)-2,4-Nonadienal, pea / soy | 2.18 +/- 0.53 / 3.87 +/- 0.52 | ug/L | same | Table 2C | measured, level_only | none |
| (E)-2-Octenal, pea / soy | 9.84 +/- 2.17 / 11.88 +/- 1.20 | ug/L | same | Table 2C | measured, level_only | none |
| Methoxypyrazine, pea milk | 4.13 +/- 0.65 | ug/L | same | Table 2C | measured, level_only | none; pre-formed in the seed (0.039 +/- 0.003 ug/g dry seed, text), not made in processing |
| Odour thresholds used | hexanal 4.5; nonanal 3.5; 2-pentylfuran 4.8; 1-hexanol 200; 1-octen-3-ol 7; (E,E)-2,4-decadienal 0.05; (E,E)-2,4-nonadienal 0.06; 1-octen-3-one 0.007; methoxypyrazine 0.02 | ug/L in water | van Gemert 2011 | Table 2C | secondary (compilation) | none; comparable with the repo's threshold files if wanted |
| Free linoleic acid, pea / soy milk | "about 14,300 and 58,000" | ug/L | milk at 2 % protein | §3.6 text | measured, level_only (precision as printed: two significant figures) | none; **the precursor pool** for a LOX module: 14,300 ug/L / 280.45 = 51 uM (pea), 207 uM (soy) |
| FFA of seed lipid | 1.75 +/- 0.11 (pea), 1.56 +/- 0.08 (soy) | g / 100 g lipid | seed | Table 1A | measured | none |
| FFA composition, pea | C18:2 38.94 %, C18:3 16.87 %, C18:1 11.37 %, C16:0 19.56 %, C18:0 8.22 % | % of FFA | seed | Table 1B | measured | none |
| Seed lipid, pea / soy | 2.17 / 16.75 | % | seed | Table 1A | measured | none |
| LOX-2 specific activity, pea / soy milk | 2160 +/- 38 / 1445 +/- 30 | U (0.001 A234/min) per mg milk protein | pH 6.8, 25 C, ~2.4 mM linoleic acid + Tween-20, 5 mg protein / 30 mL | Table 4B | measured | none; **the enzyme datum** for a LOX module, in an absorbance unit |
| LOX-3 specific activity, pea / soy | 150 +/- 6 / 122 +/- 4 | U (0.001 A280/min) per mg | pH 7.1, 25 C | Table 4B | measured | none |
| LOX-1, pea / soy | ND / 14010 +/- 563 | U (0.001 A234/min) per mg | pH 9.0 | Table 4B | measured | none |
| HPL, pea / soy | 0.15 +/- 0.03 / 0.33 +/- 0.05 | umol/(min mg protein) | 13-HPOD 33 uM, pH 6.5, 25 C, seed extract | Table 4B | measured | none |
| ADH, pea / soy | 178 +/- 12 / 31 +/- 3 | U (0.001 A340/min) per mg | ethanol -> acetaldehyde direction, pH 8.5 | Table 4B | measured | none |
| Lipase | ND | IU/mg | pH-stat, TAG | Table 4B | measured (null) | none |
| Within-study ratio: (hexanal + 1-hexanol) soy / pea | 1.44 (relative-quantity basis: 549.43 / 380.01 = 1.446); 1.31 on the calibrated basis (720.93 / 551.57) | x | same milks | §3.7 text; Tables 2B/2C | measured ratio | none; the authors' key ratio; note it moves with the quantification basis (flag 3) |
| Within-study ratio: LOX-2 x HPL, soy / pea | 1.46 (recomputed 1.47) | x | Table 4B | §3.7 | derived from measured | none; the authors match it to 1.44 above |
| Within-study ratio: free linoleic acid soy / pea | 4.06 | x | text | §3.6 | derived | none; the argument that lipid is not limiting |
| Molar yield of C6 products on free linoleic acid, pea milk | 10.6 % (calibrated: 1.64 uM hexanal + 3.79 uM 1-hexanol = 5.43 uM on 51.0 uM); 7.3 % (relative basis) | mol % | end of process, no heat | derived from Tables 2B/2C and §3.6 | derived, level_only | none; a **conversion bound** for a LOX module (pea); soy: 3.5 % (calibrated) / 2.6 % (relative) on 207 uM |
| Hexanol : hexanal molar ratio, pea / soy | 2.31 / 0.64 (calibrated); 1.67 / 0.46 (relative) | mol/mol | same | derived from 2B/2C | derived | none; ADH 5.74x higher in pea (Table 4B) is the authors' explanation |
| Enzymatic share of 1-octen-3-ol | 64 % (pea), 68 % (soy) | % | raw vs boiled-seed milk | Table 3 / §3.4 | derived from measured | none |

## 5. Flags

1. **The number the LOX programme most wants is missing: hexanal in the cooked-seed milk.** Table 3
   reports only 1-octen-3-ol and 1-octen-3-one for the enzyme-inactivated control, so the
   enzymatic share of hexanal (the analogue of the 64 % for 1-octen-3-ol) cannot be taken from this
   paper.
2. **Two quantifications, two hexanal numbers** (140.41 relative vs 164.18 ug/L calibrated for pea;
   374.06 vs 437.39 for soy). The calibration was run in 2 % fat cow's milk, whose lipid (2 %)
   differs from pea milk (0.19 %) and soy milk (1.06 %); the partition of hexanal into 2 % milk fat
   at 50 C makes the calibrated pea number, if anything, an overestimate. Use the calibrated value
   as the level and carry both.
3. **The authors' headline ratios are computed on the relative-quantity basis** ("379 and 549 ug/L"
   = Table 2B hexanal + 1-hexanol), not on the calibrated contents they present as the key numbers;
   on the calibrated basis soy/pea is 1.31, not 1.44, and the match to the LOX-2 x HPL ratio (1.46)
   loosens. The direction of the argument survives; the "quite close" does not.
4. **"Free linoleic acid ... about 14,300 and 58,000 ug/L"** is printed in the text only, with no
   table, SD or method line for how the FFA of the milk (as opposed to the seed) was measured.
   Cross-check from Table 1A/1B: pea seed FFA = 2.17 x 1.75 / 100 = 0.038 g / 100 g seed, x 38.94 %
   C18:2 = 14.8 mg free linoleic acid per 100 g seed, into ~0.8 L slurry = 18.5 mg/L if fully
   recovered; 14,300 ug/L is 77 % of that — consistent. For soy the same arithmetic gives
   115 mg/L vs 58 mg/L stated (50 %), also plausible after filtration.
5. **LOX substrate concentration is ambiguous** (one 157.2-uL stock "diluted to 200 mL ... and
   200 mL ..., respectively"): 2.5 mM if one stock per buffer, half that if split. Either way the
   assay is at substrate saturation for LOX (Km typically tens of uM), so the specific activities
   are Vmax-like.
6. **LOX unit is an absorbance rate, not a molar rate,** and the optical geometry of a 30-mL
   reaction is not described. A conversion needs epsilon(234) (~25,000 M^-1 cm^-1, not printed) and
   the reaction volume: 2160 U/mg = 2.16 A/min/mg -> 86 uM/min per mg protein -> ~2.6 umol/min per
   mg milk protein if the whole 30 mL is the reaction volume. Record the printed unit; do not
   promote the conversion.
7. **LOX-3 read at 280 nm** (oxodiene / ketodiene product) is a different observable from LOX-1/-2 at
   234 nm; the LOX-2 and LOX-3 numbers are not additive.
8. **ADH was assayed in the ethanol-oxidising direction** (NAD+, pH 8.5), not hexanal reduction; the
   5.74x pea/soy ratio is transferred to hexanal -> 1-hexanol by the authors' assumption.
9. **Internal inconsistencies:** text says 60 soy compounds, Table 2A says 68; the text's "4.04 +/-
   0.31 ug/L if it is assumed that all pyrazine was transferred" does not follow from 0.039 ug/g x
   100 g into 800 g slurry (4.9 ug/L) and implies ~0.97 L of milk; the HPL unit prints "1 uM of
   substrate" for what must be 1 umol.
10. **Pea LOX-2 identification rests on 2 unique peptides** (Table 4A; LOX-3 has 32). The
    "LOX-2" activity at pH 6.8 is an assay-pH definition and may include LOX-3; the proteomics does
    not isolate which isoform carries the pH-6.8 activity.
11. **No time course anywhere in the volatile data.** Processing is ~3 min blending + filtration +
    15 min centrifugation at 4 C, plus 10 + 30 min at 50 C in the SPME vial with active enzyme
    present (no inactivation step); the hexanal level therefore integrates formation during the
    SPME incubation too. A LOX module cannot extract a rate from this paper; it can only be checked
    against the end level.
12. **Pea protein is present at 2 % during SPME** (the milk itself), so every level is net of
    protein binding at 50 C; this is exactly the situation the matrix layer describes, and it means
    the hexanal here is a headspace-available level, not a formed amount.
13. **Journal pre-proof.** Page and volume numbers absent; values could change in the version of
    record (not on disk).
14. **Registry gaps** the roadmap list implies: (E,E)-2,4-decadienal and octanal are not in
    `data/keys/compounds.yml` although the task brief assumed they were; the pea methoxypyrazine
    here is the isopropyl-methyl congener, not the registered isobutyl one.
