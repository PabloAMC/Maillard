# Whitfield & Mottram 2001 — EXTRACTION (fed norfuraneol + cysteine or H2S at pH 6.5, 140 C / 60 min)
### The pH-6.5 twin of Whitfield & Mottram 1999 (pH 4.5): same lab, same charges, same instrument.

**Source on disk:** `data/articles/whitfield2001.pdf` (owner's download, 2026-09-07). Read-only extraction
from `pdftotext -layout`; the text layer of both tables is clean (all rows re-typed below and the class
totals re-summed). Repo status before this dossier: FIT_HOLDOUT_DECLARATION lists this paper as
**HOLD-OUT** (the pH-6.5 collapse of MFT) with the note that its H2S column should not carry a standalone
row; k3 §A.3.1 already carries "NF + cysteine, pH 6.5: MFT nd (< 0.1 µg/10 mg) = < 0.0010 mol %".

## 0. Identity

| field | value |
|---|---|
| Title | "Heterocyclic Volatiles Formed by Heating Cysteine or Hydrogen Sulfide with 4-Hydroxy-5-methyl-3(2H)-furanone at pH 6.5" |
| Authors | Frank B. Whitfield (Food Science Australia, North Ryde) and Donald S. Mottram (Reading) |
| Venue | J. Agric. Food Chem. 2001, 49, 816-822 |
| DOI | 10.1021/jf0008644 |
| Naming | the paper calls norfuraneol "HMF" (= 4-hydroxy-5-methyl-3(2H)-furanone, NF in the repo). NOT hydroxymethylfurfural. |
| Companion | ref 8 = Whitfield & Mottram 1999, JAFC 47, 1626-1634 (same experiment at pH 4.5), the repo's FIT source for the NF channel |

## 1. Why it matters

The model's MFT-formation step (norfuraneol + H2S -> MFT) is fitted to Whitfield 1999 at pH 4.5 / 140 C
and Hofmann 1998 at pH 5 / 145 C and has no measured pH dependence. This paper is the SAME fed-norfuraneol
experiment (same charges, vessel, oven, trap, GC-MS, internal standard, reporting unit) two pH units
higher. It is therefore the only like-for-like pH pair on the NF -> MFT step in the corpus. The result
is a non-detect: at pH 6.5 neither system produced any furanthiol, thiophenethiol, or furyl disulfide
above the detection limit, and neither smelled meaty.

## 2. Methods as they matter to a model

- **Buffer:** 0.5 M phosphate, pH 6.5, distilled water.
- **Cysteine system:** 11.4 mg HMF (norfuraneol, MW 114.10 = **0.0999 mmol**) + 12.1 mg cysteine
  (MW 121.16 = **0.0999 mmol**) in 2 mL buffer -> **50 mmol/L each, 1:1**.
- **H2S system:** "1 mL of a solution containing 12.1 mg of cysteine, in buffer, and 1 mL of the
  saturated solution of hydrogen sulfide (~6.6 mg)". ⚠ As printed this sentence omits HMF and includes
  cysteine, which contradicts the title, Table 1's column headings ("cysteine" vs "H2S"), the whole
  Results section ("reaction between HMF and hydrogen sulfide") and the H2S column's product spectrum
  (furanones, thiophenones and 2-acetylthiophene need the furanone carbon skeleton). Read as a typo for
  "11.4 mg of HMF"; the intended charge is presumably 0.1 mmol HMF + ~0.19 mmol H2S (6.6 mg / 34.08 =
  0.194 mmol; a 0.2 M saturated solution x 1 mL = 0.2 mmol, consistent) in 2 mL -> **~50 mmol/L NF,
  ~97 mmol/L H2S, H2S:NF ~2:1**, the same 1:2 ratio the repo records for the 1999 H2S run. The HMF
  charge in the H2S ampule is nevertheless UNSTATED in this paper.
- **H2S solution:** saturated (~0.2 M) by passing H2S through pH 6.5 phosphate buffer at 0 C.
- **Vessel / heating:** 5 mL glass ampules, flame-sealed (so the H2S and volatiles stay in), oven at
  **140 C for 60 min**. No stirring mentioned. Atmosphere: sealed air headspace, not stated further.
- **Replicates:** each reaction in duplicate; both values printed in Table 1; Table 2 gives means.
- **Isolation:** after cooling, diluted with 20 mL 0.5 M phosphate pH 6.5; internal standard methyl
  decanoate 100 µg in 0.1 mL ethanol added; dynamic headspace onto Tenax GC while stirring at 60 C.
- **GC-MS:** Varian 1440 GC, 50 m x 0.32 mm BP5, thermal desorption 260 C / 5 min, cryofocused; Varian-MAT
  311A double-focusing MS. Identification: MS + LRI vs authentics where available (footnote b).
- **Quantification: RELATIVE.** "approximate concentrations ... determined by comparing their GC-MS
  chromatogram peak areas with the area of the internal standard, methyl decanoate, which was taken as
  100 µg, and assuming all response factors were 1." Reported as **µg per 10 mg of HMF used**.
  Detection limit **0.1 µg/10 mg HMF** (3x noise); "trace" (tr) = 0.1-1 µg/10 mg HMF; nd = not
  detected. Dynamic headspace at 60 C means the numbers are headspace-recoverable amounts, not solution
  concentrations. No stable-isotope dilution, no calibration.
- **Conversion the repo uses:** 10 mg HMF = 87.64 µmol, so mol % (product/NF fed) = µg / MW_product /
  0.8764. For MFT (MW 114.17): 0.1 µg/10 mg = **0.0010 mol %**; 1 µg = 0.010 mol %.
- **Sensory:** three assessors. HMF + cysteine: "caramel, roasted, and fried bacon rind"; HMF + H2S:
  "caramel, metallic, and medicinal". "Neither system had the meatlike aromas previously observed when
  these mixtures were reacted at pH 4.5".

## 3. Tables re-typed

### Table 1. "Volatile Compounds Obtained from Reactions between 4-Hydroxy-5-methyl-3(2H)-furanone and Cysteine or Hydrogen Sulfide"

Columns: approximate concn (µg/10 mg of HMF), two duplicate values per system; method of ID; LRI.
Footnote a: "Concentrations (µg/10 mg of HMF obtained by comparing GC-MS peak area with that from 100 µg
of methyl decanoate internal standard added to the HMF solution before volatile collection); duplicate
analyses are shown; nd, not detected (limit of detection ~0.1 µg/10 mg of HMF); tr, between 0.1 and 1
µg/10 mg of HMF." Footnote b: MS + LRI = authentic compound; MS = literature spectrum; ms = interpreted
spectrum. Footnote e: "Previously incorrectly reported as the 2,4-isomer (8)." Mass-spectral columns
omitted here.

| no. | compound | cysteine (dup 1, dup 2) | H2S (dup 1, dup 2) | ID | LRI |
|---:|---|---:|---:|---|---:|
| 1 | 2,3-pentanedione | 7, 11 | 16, nd | MS+LRI | 680 |
| 2 | 3-penten-2-one | tr, nd | 8, nd | MS+LRI | 755 |
| 3 | 4,5-dimethyloxazole | 1, 2 | nd, nd | MS+LRI | 771 |
| 4 | 2-hexanone | 1, 1 | 4, 29 | MS+LRI | 792 |
| 5 | 3,4-hexanedione | 3, 3 | 6, nd | MS+LRI | 800 |
| 6 | 4,5-dihydro-2-methyl-3(2H)-furanone | 6, 8 | 15, 8 | MS | 806 |
| 7 | **3-mercaptobutan-2-one** | **2, 5** | **1, nd** | MS+LRI | 816 |
| 8 | methylpyrazine | 3, 2 | nd, nd | MS+LRI | 825 |
| 9 | 2,4,5-trimethyloxazole | 20, 14 | nd, nd | MS+LRI | 846 |
| 10 | 2,5-dimethylthiophene | 6, 22 | 24, 4 | MS+LRI | 867 |
| 11 | 2,4-dimethylthiophene | nd, nd | 9, 11 | MS+LRI | 875 |
| 12 | 2-methyl-2-thiazoline | 6, 4 | nd, nd | MS+LRI | 876 |
| 13 | 2,4-hexanedione | tr, nd | 9, 7 | MS+LRI | 880 |
| 14 | 2,4-dimethylthiazole | tr, nd | nd, nd | MS+LRI | 883 |
| 15 | 3,4-dimethylthiophene | 1, 4 | 1, nd | MS+LRI | 887 |
| 16 | **3-mercaptopentan-2-one** | **2, 4** | **nd, nd** | MS+LRI | 902 |
| 17 | 2-methyl-2-cyclopenten-1-one | 4, 3 | 6, nd | MS+LRI | 905 |
| 18 | 2,5(or 2,6)-dimethylpyrazine | 16, 11 | nd, nd | MS+LRI | 910 |
| 19 | 2,5-dimethylthiazole | 1, nd | nd, nd | MS+LRI | 915 |
| 20 | 2,3-dimethylpyrazine | 2, nd | nd, nd | MS+LRI | 919 |
| 21 | 4-ethyl-2,5-dimethyloxazole | 1, nd | nd, nd | MS | 919 |
| 22 | 2-ethyl-4,5-dimethyloxazole | 2, 1 | nd, nd | MS | 925 |
| 23 | 4,5-dimethylthiazole | 13, 11 | nd, nd | MS+LRI | 933 |
| 24 | 2-ethyl-1H-pyrrole | 1, 2 | nd, nd | MS+LRI | 944 |
| 25 | 1-(2-furyl)-2-propanone | 3, 4 | nd, nd | MS | 952 |
| 26 | 2-ethyl-5-methylthiophene | 1, 4 | 9, 17 | MS+LRI | 959 |
| 27 | 4,5-dihydro-5-methylthiophen-3(2H)-one | 6, 6 | 15, 20 | MS | 982 |
| 28 | 4,5-dihydro-2-methylthiophen-3(2H)-one | 39, 40 | 48, 75 | MS+LRI | 990 |
| 29 | 2,4,5-trimethylthiazole | 7, 7 | nd, nd | MS+LRI | 997 |
| 30 | 2-ethyl-(5 or 6)-methylpyrazine | 8, 7 | nd, nd | MS+LRI | 997 |
| 31 | 2,3,5-trimethylpyrazine | 10, 7 | nd, nd | MS+LRI | 1002 |
| 32 | 1-(2-furyl)-1-propanone | 3, 3 | 4, nd | MS+LRI | 1008 |
| 33 | 4,5-dihydro-2,5-dimethylthiophen-3(2H)-one (E or Z) | 31, 28 | 17, 29 | MS | 1016 |
| 34 | 4,5-dihydro-2,5-dimethylthiophen-3(2H)-one (E or Z) | 5, 8 | 4, 6 | MS | 1027 |
| 35 | 2,3-dimethyl-2-cyclopenten-1-one | 6, 4 | 7, 10 | MS+LRI | 1040 |
| 36 | 1-(5-methyl-2-furyl)-2-propanone | 2, 2 | 6, 7 | MS | 1047 |
| 37 | 2-ethyl-3,6-dimethylpyrazine | 3, 2 | nd, nd | MS+LRI | 1078 |
| 38 | 3-methyl-1,2-dithiolan-4-one | nd, nd | 6, nd | MS | 1071 |
| 39 | 4,5-dihydro-2-ethylthiophen-3(2H)-one | 9, 8 | 17, 19 | MS | 1082 |
| 40 | nitrogen compound A | 34, 31 | nd, nd | ms | 1089 |
| 41 | 2-acetylthiophene | nd, nd | 55, 76 | MS+LRI | 1092 |
| 42 | 3,5-dimethyl-1,2-dithiolan-4-one (E or Z) | 1, nd | 7, 25 | MS | 1098 |
| 43 | 2-ethyl-5-methyl-4,5-dihydrothiophen-3(2H)-one | 4, 4 | 7, nd | ms | 1104 |
| 44 | 2-thiazolyl-1-propanone | 2, 2 | nd, nd | MS | 1119 |
| 45 | 2-formyl-5-methylthiophene | 3, 4 | nd, nd | MS+LRI | 1124 |
| 46 | 3,4,5-trimethyl-2-furfural | 2, 2 | nd, nd | MS | 1130 |
| 47 | 1-(3-thienyl)-2-propanone | 1, nd | nd, nd | MS | 1134 |
| 48 | 3,5-dimethyl-1,2,4-trithiolane (E or Z) | 7, 8 | nd, nd | MS+LRI | 1138 |
| 49 | 3,5-dimethyl-1,2,4-trithiolane (E or Z) | 7, 14 | nd, nd | MS+LRI | 1144 |
| 50 | 1-(dimethyl-2-furyl)-2-propanone | 5, 5 | 6, 15 | MS | 1151 |
| 51 | 2-acetyl-5-methylthiophene | 1, nd | 6, nd | MS+LRI | 1157 |
| 52 | nitrogen compound B | 3, 2 | nd, nd | ms | 1167 |
| 53 | 3-ethyl-1,2-dithiolan-4-one | nd, nd | 1, nd | MS | 1169 |
| 54 | sulfur compound MW 142 | 1, nd | 9, 19 | ms | 1171 |
| 55 | 1-(3-thienyl)-1-propanone | 2, nd | nd, nd | MS+LRI | 1183 |
| 56 | 2,3-dihydro-6-methylthieno[2,3c]furan | 1, nd | nd, nd | MS+LRI | 1199 |
| 57 | 3-ethyl-2-formylthiophene | nd, nd | 22, 21 | MS | 1206 |
| 58 | nitrogen compound C | 30, 23 | nd, nd | ms | 1219 |
| 59 | 3-ethyl-5-methyl-1,2,4-trithiolane (E or Z) | 1, nd | 4, 7 | MS | 1242 |
| 60 | 3-ethyl-5-methyl-1,2,4-trithiolane (E or Z) | 1, 2 | 5, 7 | MS | 1250 |
| 61 | 3-methyl-1,2,4-trithiane | 3, 3 | nd, nd | MS+LRI | 1254 |
| 62 | nitrogen compound D | 4, 2 | nd, nd | ms | 1314 |
| 63 | a dihydrothienothiophene | 8, 5 | nd, nd | MS | 1319 |
| 64 | a methylthienothiophene | 4, 3 | 9, 5 | MS | 1357 |
| 65 | a dihydromethylthienothiophene | 1, nd | 4, nd | MS | 1378 |
| 66 | a dihydromethylthienothiophene | 4, 2 | nd, nd | MS | 1409 |
| 67 | a dihydromethylthienothiophene | 7, 3 | nd, nd | MS | 1418 |
| 68 | **3-(2-methyl-3-furyldithio)-2-butanone** | **1, 1** | **nd, nd** | MS+LRI | 1501 |
| | **total** | **369, 354** | **367, 417** | | |

**Species the model carries that do NOT appear in Table 1 (i.e. nd, < 0.1 µg/10 mg HMF, in both
systems):** 2-methyl-3-furanthiol (MFT), 2-furfurylthiol (FFT), 2-mercapto-3-pentanone,
2-methyl-3-thiophenethiol, bis(2-methyl-3-furyl) disulfide, 2-furfural, hydrogen sulfide (not a
Tenax-trappable analyte here in any case). Norfuraneol itself (the reactant) is not reported.
Compound 68 is the one MFT-containing product: a mixed disulfide of MFT with 3-mercaptobutan-2-one,
1 µg/10 mg HMF in both cysteine duplicates, nd in H2S.

### Table 2. "Comparison of Classes of Compounds Found in Reactions between HMF and Cysteine or Hydrogen Sulfide at pH 4.5 and 6.5"

Footnote a: "Mean concentrations (µg/10 mg of HMF) from duplicate analyses; nd, not detected."
Footnote b: "From ref 8" (= Whitfield & Mottram 1999). Column sums re-checked: 580, 345, 616, 372 all
reproduce.

| compound class | cysteine pH 4.5 (b) | cysteine pH 6.5 | H2S pH 4.5 (b) | H2S pH 6.5 |
|---|---:|---:|---:|---:|
| dithiolanones and dithianones | 41 | 1 | 443 | 20 |
| trithiolanes and trithianes | nd | 23 | nd | 12 |
| thiophenes | 94 | 25 | 21 | 128 |
| thiophenones | 94 | 94 | 39 | 129 |
| aliphatic and alicyclic ketones | 28 | 22 | 39 | 51 |
| **thiols and mercaptoketones** | **294** | **7** | **44** | **1** |
| disulfides | 24 | 1 | 18 | nd |
| furan derivatives | 5 | 23 | 12 | 31 |
| thiazoles and oxazoles | nd | 42 | nd | nd |
| pyrazines | nd | 36 | nd | nd |
| other nitrogen compounds | nd | 71 | nd | nd |
| **total** | **580** | **345** | **616** | **372** |

Note: Table 2's pH 6.5 totals (345, 372) are lower than the Table 1 duplicate means (361.5, 392); the
class table evidently omits a few unclassified entries (e.g. compound 54). Use Table 1 for compound
levels and Table 2 only for class ratios.

Other numbers in the text: 3-ethyl-2-formylthiophene (57) at pH 4.5 was 17-57 µg/10 mg HMF in the
cysteine system and nd with H2S; at pH 6.5 the reverse (21-22 with H2S, nd with cysteine). Class shares
at pH 6.5: thiophenones 25-27 % (cys) / 29-36 % (H2S) of total volatiles; thiophenes 4-10 % / 31-34 %;
trithiolanes 4-7 % / 2-3 %; dithiolanones trace / 4-6 %; nitrogen compounds 41 % of the cysteine total
(pyrazines 8-11 %, thiazoles 6 %, oxazoles 5-7 %, unidentified N compounds 16-19 %).

## 4. What the repo could take

### 4.1 Fed-intermediate rows at pH 6.5, 140 C, 60 min (the pH twin of the FIT rows)

Conversion: 10 mg HMF = 87.64 µmol; mol % = µg / MW / 0.8764.

| system | product | µg/10 mg HMF (duplicates) | mol % of NF fed | role |
|---|---|---:|---:|---|
| NF + cysteine 1:1, pH 6.5 | MFT (free) | nd, nd (< 0.1) | **< 0.0010** | the HOLD-OUT row already declared |
| NF + H2S ~1:2, pH 6.5 | MFT (free) | nd, nd (< 0.1) | **< 0.0010** | same, but the HMF charge in this ampule is unstated (§2) — do not carry as a standalone row, per the existing declaration; usable as a second non-detect |
| NF + cysteine 1:1, pH 6.5 | 3-(2-methyl-3-furyldithio)-2-butanone (68) | 1, 1 | 0.0053 (MW 216.3) | MFT-equivalent bound in a mixed disulfide: ~0.53 µg MFT-eq/10 mg -> ~0.005 mol % |
| NF + cysteine 1:1, pH 6.5 | 3-mercaptopentan-2-one | 2, 4 (mean 3) | **0.029** (MW 118.2) | the only mercaptoketone that survives at pH 6.5 |
| NF + H2S, pH 6.5 | 3-mercaptopentan-2-one | nd, nd | < 0.0010 | |
| NF + cysteine 1:1, pH 6.5 | 3-mercaptobutan-2-one | 2, 5 (mean 3.5) | 0.038 (MW 104.2) | |
| NF + H2S, pH 6.5 | 3-mercaptobutan-2-one | 1, nd | ~0.005 | |
| both | 2-mercapto-3-pentanone, FFT, 2-methyl-3-thiophenethiol, furyl disulfides | nd | < 0.0010 each | |
| NF + cysteine 1:1, pH 6.5 | total volatiles | 369, 354 | — | |
| NF + H2S, pH 6.5 | total volatiles | 367, 417 | — | |

### 4.2 Within-study pH ratios, pH 4.5 -> 6.5 (same lab, method, charges; 140 C / 60 min)

The pH-4.5 compound levels come from Whitfield 1999 via the repo's k3 §A.3.1 / §A.10 (MFT 15, 15 µg/10 mg
NF in the cysteine system = 0.150 mol %; 0.120 mol % in the H2S system; 3-mercapto-2-pentanone 74.5 and
2-mercapto-3-pentanone 77.5 µg/10 mg) — not re-verified against `whitfield1999.pdf` in this pass. The
class-level pH-4.5 numbers are printed in this paper's Table 2.

| quantity | pH 4.5 | pH 6.5 | ratio 4.5 / 6.5 |
|---|---:|---:|---:|
| MFT, cysteine system (free MFT) | 15 µg/10 mg (0.150 mol %) | < 0.1 (< 0.0010 mol %) | **>= 150** |
| MFT, cysteine system (free + disulfide-bound MFT-eq) | >= 15 | <= ~0.63 | >= ~24 (lower bound if 1999's MFT-bearing disulfides are not added back on the pH-4.5 side; they exist there — disulfide class 24 µg) |
| MFT, H2S system | 0.120 mol % | < 0.0010 mol % | >= 120 (charge caveat) |
| 3-mercapto-2-pentanone, cysteine | 74.5 | 3 | **~25** |
| 2-mercapto-3-pentanone, cysteine | 77.5 | < 0.1 | **>= 775** |
| thiols + mercaptoketones class, cysteine | 294 | 7 | **42** |
| thiols + mercaptoketones class, H2S | 44 | 1 | 44 |
| dithiolanones + dithianones, cysteine / H2S | 41 / 443 | 1 / 20 | 41 / 22 |
| disulfides, cysteine / H2S | 24 / 18 | 1 / nd | 24 / > 180 |
| thiophenones, cysteine / H2S | 94 / 39 | 94 / 129 | 1.0 / 0.30 (thiophenones do NOT fall) |
| thiophenes, cysteine / H2S | 94 / 21 | 25 / 128 | 3.8 / 0.16 |
| total volatiles, cysteine / H2S | 580 / 616 | 345 / 372 | **1.7 / 1.7** |
| trithiolanes + trithianes, cysteine / H2S | nd / nd | 23 / 12 | new at 6.5 |
| N-heterocycles (thiazoles+oxazoles, pyrazines, other N), cysteine | nd | 149 (41 % of total) | new at 6.5 |

### 4.3 Directional claims with numbers

- **MFT from fed NF collapses >= 150x between pH 4.5 and 6.5 while total volatile output falls only
  1.7x.** The NF carbon is still converted (thiophenones unchanged at 94; furan derivatives up 5 -> 23);
  what changes is where the sulfur goes. This is a pH effect on the thiol-forming steps, not on NF
  consumption. Authors: "lower pH values are essential for the reaction of hydrogen sulfide with such
  HMF products to form mercaptoketones, furan- and thiophenethiols, and dithiolanones."
- **The pH-6.5 sulfur sinks are trithiolanes/trithianes (23, 12 µg), thiophenes (25, 128) and
  thiophenones (94, 129), not thiols or disulfides.** For a model that needs to know where H2S goes when
  MFT is not made, this is the measured partition (relative units).
- **Ammonia becomes available at pH 6.5**: 149 µg of N-heterocycles in the cysteine system vs nd at pH
  4.5 — cysteine's degradation route is itself pH-switched.
- **The two mercaptopentanone isomers respond differently to pH**: 3-mercapto-2-pentanone survives at
  6.5 (falls ~25x), 2-mercapto-3-pentanone (the repo's NF-diagnostic isomer) is gone (>= 775x). If the
  model gives both isomers the same pH law it will miss this.
- Both duplicates of the MFT-bearing mixed disulfide (68) are 1 µg: MFT is FORMED in trace at pH 6.5
  and is captured as a disulfide with 3-mercaptobutan-2-one. So the true pH-6.5 MFT yield is not zero
  but ~0.005 mol % as MFT-equivalent — still ~30x below 0.150.

### 4.4 Role already declared

FIT_HOLDOUT_DECLARATION: **HOLD-OUT** (the >= 150x collapse). This dossier changes nothing in that
declaration; it adds the isomer-specific ratios (§4.2) and the sink partition (§4.3) as further
directional hold-out shapes, and the disulfide-bound MFT bound as a caveat on the 150x number.

## 5. Caveats

1. **Relative quantification, response factors assumed 1**, dynamic headspace at 60 C, Tenax. Levels
   are µg-equivalents of methyl decanoate; mol % values above inherit that. Ratios within the same
   compound across pH are far safer than any level.
2. **The H2S-system charge is unstated** (Methods sentence omits HMF; §2). The H2S column's product
   spectrum proves NF was present, but its amount is an inference (11.4 mg by analogy with the cysteine
   ampule and with 1999). Existing repo rule stands: no standalone H2S-column row.
3. **Single time point (60 min)**, so nothing on the removal of thiols; the absence of furyl disulfides
   AND of free MFT at pH 6.5 together argue that MFT was not made and then oxidised — it was mostly not
   made (except the 1 µg captured in 68).
4. **Duplicate spread is wide for some entries** (e.g. 2,5-dimethylthiophene 6 vs 22; 2-hexanone 4 vs 29;
   3,5-dimethyl-1,2-dithiolan-4-one 7 vs 25); treat anything below ~5 µg as order-of-magnitude.
5. **Table 2 class totals (345, 372) do not equal Table 1 totals (361.5, 392)** — a few compounds are
   unclassified; use Table 1 for levels.
6. 0.5 M phosphate: the pH-4.5 comparator (1999) used the same buffer strength per this paper's
   "as described previously (8)"; buffer catalysis is therefore common to both arms.
7. The pH-4.5 compound levels quoted in §4.2 are taken from the repo's k3 inventory, not from this PDF;
   only the class totals (Table 2) and the 17-57 µg figure for compound 57 are printed here.
