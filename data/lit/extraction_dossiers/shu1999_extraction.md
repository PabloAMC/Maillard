# Shu 1999 — EXTRACTION (neat serine, threonine and their 1:1 mixture at 12 % moisture, 120 C / 4 h and 300 C / 7 min; pyrazines in ppm of amino acid by GC-FID)
### Pyrazines without a sugar: serine gives pyrazine, ethyl- and diethylpyrazines; threonine gives 2,5-dimethyl-, trimethyl- and ethyldimethylpyrazines.

**Source on disk:** `data/articles/shu1999.pdf` (owner's download, 2026-09-08). Read-only extraction from
the pypdf text layer in the scratchpad. Table 1's sparse rows lose their column positions in the text
layer; the assignment below was fixed by the printed column totals (all six re-sum exactly) and by the
Results text, which names the major and minor pyrazines of each arm. Not the same paper as
`shu1988_extraction.md`. Repo status before this dossier: roadmap §5b lists the ethyl- and
trimethylpyrazines with "the aldehyde-addition step; amino-acid identity" missing.

## 0. Identity

| field | value |
|---|---|
| Title | "Pyrazine Formation from Serine and Threonine" |
| Author | Chi-Kuen Shu (R. J. Reynolds Tobacco, Winston-Salem) |
| Venue | J. Agric. Food Chem. 1999, 47, 4332-4335 |
| DOI | 10.1021/jf9813687 |
| Naming | "ppm" = µg of pyrazine per g of amino acid charged (mass basis); "2,5-dimethyl-" and "2,6-dimethyl-" are resolved (DB-Wax), unlike the DB-5 work of Adams 2008 |
| Cited by | Adams 2008 (ref 12) as the amino-acid-only route |

## 1. Why it matters

For Programme 6 the question is which aminoketone plus which aldehyde makes which alkylpyrazine.
Shu's system has no carbohydrate at all: the aminoketones and the aldehydes both come from the
hydroxy-amino acid, so the product pattern is a clean test of the combinatorics. Threonine, whose
decarbonylation/dehydration gives aminoacetone (the same aminoketone the trunk's MGO + glycine
Strecker step makes), gives 2,5-dimethylpyrazine as the main product with trimethyl- and
2-ethyl-3,6-dimethylpyrazine beside it; serine, which gives aminoacetaldehyde, gives pyrazine with
ethylpyrazine beside it. The 120 C / 4 h arms are the only cooking-relevant ones and show 0.002-0.016
mol % conversion; the 300 C arms are pyrolysis. The mixture arm shows that cross-condensation
(methylpyrazine, 2-ethyl-6-methylpyrazine) dominates when both aminoketone pools are present, and that
the mixture makes 1.6-2.8x more pyrazine than the average of its parts.

## 2. Methods as they matter to a model

- **Charge:** 300 mg amino acid + 36 µL water ("12 % moisture level") in an enclosed Parr vessel. No
  buffer, no sugar, no pH. Arms: L-serine; L-threonine; 50/50 w/w serine/threonine (150 + 150 mg).
  Moles: serine 300 mg / 105.09 = **2.855 mmol**; threonine 300 / 119.12 = **2.518 mmol**; mixture
  1.427 + 1.259 = **2.686 mmol** (mean molar mass 111.7 g/mol). The system is a wet solid, not a
  solution; no mmol/L applies.
- **Heating:** oven, **120 C for 4 h** or **300 C for 7 min**. Cooled to room temperature.
- **Isolation:** CH2Cl2 4 x 5 mL; 300 C extracts injected directly; 120 C extracts concentrated to
  0.5 mL under N2 (an extra step for the low-temperature arms only).
- **Identification:** GC/MS, DB-Wax 60 m x 0.32 mm x 0.15 µm, 50-200 C at 6 C/min, EI 70 eV.
- **Quantification: RELATIVE, GC-FID**, n-hexadecane internal standard, **"it was assumed that the FID
  response factor of n-hexadecane was the same as those of the alkylpyrazines"**. Reported as
  **ppm of serine and/or threonine used**. No replicates stated; no detection limit stated; "trace"
  pyrazines (e.g. 2,5- and 2,6-dimethylpyrazine from serine) not quantified.
- **Conversion the repo uses:** mol % of amino acid = ppm x 1e-4 x MW(amino acid) / MW(pyrazine).
  For the mixture use 111.7 for the amino acid. MW: pyrazine 80.09; methylpyrazine 94.11;
  dimethyl- and ethylpyrazine 108.14; trimethyl- and ethylmethylpyrazine 122.17; ethyldimethyl- and
  diethylpyrazine 136.19. Example: 97.6 ppm pyrazine from serine = 97.6e-4 x 105.09 / 80.09 =
  0.0128 mol %.

## 3. Tables re-typed

### Table 1. "Pyrazines Identified and the Quantitative Data, at Parts per Million" (ppm of amino acid charged; blank = not reported)

| pyrazine | 120 C/4 h Ser | 120 C/4 h Thr | 120 C/4 h Ser/Thr | 300 C/7 min Ser | 300 C/7 min Thr | 300 C/7 min Ser/Thr |
|---|---:|---:|---:|---:|---:|---:|
| pyrazine | 97.6 | | 32.3 | 1477 | | 354 |
| methyl- | 6.1 | | 32.0 | 245 | | 880 |
| 2,5-dimethyl- | | 11.4 | 9.9 | | 1100 | 898 |
| 2,6-dimethyl- | | | 3.1 | | 390 | 501 |
| ethyl- | 21.3 | | 8.2 | 1025 | | 348 |
| 2-ethyl-6-methyl- | | | 12.6 | 131 | | 1077 |
| 2-ethyl-5-methyl- | | | | | | 281 |
| trimethyl- | | 10.0 | 7.4 | | 271 | 705 |
| 2,6-diethyl- | | | | 391 | | 247 |
| 2-ethyl-3,6-dimethyl- | | 0.4 | 11.4 | | 632 | 2291 |
| 2-ethyl-3,5-dimethyl- | | | | | 83 | 458 |
| **total** | **125.0** | **21.8** | **116.9** | **3269** | **2476** | **8040** |

Column re-sums: 97.6 + 6.1 + 21.3 = 125.0; 11.4 + 10.0 + 0.4 = 21.8; 32.3 + 32.0 + 9.9 + 3.1 + 8.2 +
12.6 + 7.4 + 11.4 = 116.9; 1477 + 245 + 1025 + 131 + 391 = 3269; 1100 + 390 + 271 + 632 + 83 = 2476;
354 + 880 + 898 + 501 + 348 + 1077 + 281 + 705 + 247 + 2291 + 458 = 8040. All six reproduce, which
pins the sparse-row assignment. The text confirms: serine 120 C "pyrazine, ethylpyrazine (major),
methylpyrazine (minor)"; threonine 120 C "2,5-dimethylpyrazine, trimethylpyrazine (major),
2-ethyl-3,6-dimethylpyrazine (minor)"; serine 300 C adds 2,6-diethyl- (major) and 2-ethyl-6-methyl-
(minor); threonine 300 C adds 2,6-dimethyl- and 2-ethyl-3,5-dimethyl- (minor).

### Table 1 converted to mol % of amino acid charged (mine; FID response factors assumed 1 carry through)

| pyrazine | Ser 120 | Thr 120 | Ser/Thr 120 | Ser 300 | Thr 300 | Ser/Thr 300 |
|---|---:|---:|---:|---:|---:|---:|
| pyrazine | 0.0128 | | 0.0045 | 0.194 | | 0.049 |
| methyl- | 0.00068 | | 0.0038 | 0.027 | | 0.104 |
| 2,5-dimethyl- | | 0.00126 | 0.00102 | | 0.121 | 0.093 |
| 2,6-dimethyl- | | | 0.00032 | | 0.043 | 0.052 |
| ethyl- | 0.00207 | | 0.00085 | 0.100 | | 0.036 |
| 2-ethyl-6-methyl- | | | 0.00115 | 0.0113 | | 0.098 |
| 2-ethyl-5-methyl- | | | | | | 0.026 |
| trimethyl- | | 0.00098 | 0.00068 | | 0.026 | 0.064 |
| 2,6-diethyl- | | | | 0.030 | | 0.020 |
| 2-ethyl-3,6-dimethyl- | | 0.000035 | 0.00094 | | 0.055 | 0.188 |
| 2-ethyl-3,5-dimethyl- | | | | | 0.0073 | 0.038 |
| **sum** | **0.0156** | **0.0023** | **0.0133** | **0.362** | **0.253** | **0.769** |

### Mechanism as printed (Figures 1 and 2 are schemes, no data)

- Serine: decarbonylation + dehydration -> aminoacetaldehyde (1) -> dimer -> pyrazine. Decarbonylation
  + deamination -> glycolaldehyde; decarboxylation + deamination -> acetaldehyde. Aldol of
  glycolaldehyde + acetaldehyde -> 1,2-butanedione -> Strecker with serine -> 1-amino-2-butanone (2),
  2-aminobutanal (3). 1 + 2 or 1 + 3 -> ethylpyrazine; 2 + 3 -> 2,6-diethylpyrazine. Glycolaldehyde
  aldol/retro-aldol -> aminoacetone (4), 2-aminopropanal (5); 2 + 5 or 3 + 4 -> 2-ethyl-6-methyl-;
  1 + 4 or 1 + 5 -> methylpyrazine.
- Threonine: decarbonylation + dehydration -> aminoacetone (4) -> dimer -> 2,5-dimethylpyrazine;
  4 + 5 -> 2,6-dimethyl-. Aldol of hydroxyacetone, retro-aldol, Strecker -> 3-amino-2-pentanone (6),
  2-amino-3-pentanone (7); 6 + 4 -> 2-ethyl-3,6-dimethyl-; 7 + 4 -> 2-ethyl-3,5-dimethyl-.
  Hydroxyacetone + formaldehyde -> 3-amino-2-butanone (8); 8 + 4 -> trimethylpyrazine.

## 4. Numbers and steps the repository can use

Registry keys: `methylpyrazine`, `2_5_dimethylpyrazine`, `2_6_dimethylpyrazine`, `2_ethylpyrazine`,
`trimethylpyrazine`, `2_ethyl_3_5_dimethylpyrazine` (see flag 7 on its SMILES) exist. Pyrazine itself
has only the family class `pyrazines` (the bare alias "pyrazine" means the family) — **no molecule id**.
2-Ethyl-6-methyl-, 2-ethyl-5-methyl-, 2,6-diethyl- and 2-ethyl-3,6-dimethylpyrazine, serine, threonine:
**not in registry**.

| quantity or step | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| pyrazine from serine alone | 97.6 ppm = 0.0128 mol % | of Ser | 120 C, 4 h, 12 % moisture, neat | Table 1 | level_only (FID, RF = 1; end of cook) |
| ethylpyrazine from serine alone | 21.3 ppm = 0.0021 mol % | of Ser | same | Table 1 | level_only |
| methylpyrazine from serine alone | 6.1 ppm = 0.0007 mol % | of Ser | same | Table 1 | level_only |
| 2,5-dimethylpyrazine from threonine alone | 11.4 ppm = 0.0013 mol % | of Thr | same | Table 1 | level_only |
| trimethylpyrazine from threonine alone | 10.0 ppm = 0.0010 mol % | of Thr | same | Table 1 | level_only |
| 2-ethyl-3,6-dimethylpyrazine from threonine alone | 0.4 ppm = 0.00004 mol % | of Thr | same | Table 1 | level_only |
| total pyrazines, 120 C | 125.0 / 21.8 / 116.9 ppm | of amino acid | Ser / Thr / mixture | Table 1 | level_only |
| total pyrazines, 300 C | 3269 / 2476 / 8040 ppm (0.36 / 0.25 / 0.77 mol %) | of amino acid | 7 min pyrolysis | Table 1 | level_only (outside any cook) |
| ethyl : parent pyrazine, serine | 0.22 (120 C); 0.69 (300 C) | — | | Table 1 | within_study_ratio (the aldehyde-carrying aminoketone gains with temperature) |
| trimethyl : 2,5-dimethyl, threonine | 0.88 (120 C); 0.25 (300 C) | — | | Table 1 | within_study_ratio |
| 2-ethyl-3,6-dimethyl : 2,5-dimethyl, threonine | 0.035 (120 C); 0.57 (300 C) | — | | Table 1 | within_study_ratio |
| 2-ethyl-3,6-dimethyl : 2-ethyl-3,5-dimethyl | 7.6 (Thr 300 C); 5.0 (mixture 300 C); 3,5-isomer absent at 120 C | — | | Table 1 | within_study_ratio (same direction as Cerny 1994's I : II ~10 from MGO + alanine) |
| 2,5- : 2,6-dimethyl, threonine 300 C | 2.8 | — | | Table 1 | within_study_ratio (homo-condensation of aminoacetone over the cross with 2-aminopropanal) |
| mixture synergy | 116.9 vs (125.0 + 21.8)/2 = 73.4 -> 1.6x (120 C); 8040 vs 2872 -> 2.8x (300 C) | — | 50/50 w/w | Table 1 | within_study_ratio |
| cross products only in the mixture | methylpyrazine 32.0 (vs 6.1 Ser, 0 Thr); 2-ethyl-6-methyl 12.6; 2,6-dimethyl 3.1; 2-ethyl-3,6-dimethyl 11.4 (vs 0.4 from Thr) at 120 C | ppm | | Table 1 | level_only (the aminoketone pools cross-condense statistically, and the mixed pool is more productive than either alone) |
| serine vs threonine at 120 C | Ser total 5.7x Thr total | — | | Table 1 | within_study_ratio (authors ascribe it to melting points 222 vs 256 C) |

## 5. Flags

1. **Relative quantification:** FID with n-hexadecane, response factors assumed 1; no replicates or
   detection limits printed. Order-of-magnitude levels; within-arm ratios are safer.
2. **No sugar, no water phase, no pH.** The pot is a wet crystalline solid at 12 % moisture in a
   sealed bomb; "120 C / 4 h" is the oven, the sample's own temperature history is unknown; the
   300 C arm is pyrolysis and belongs to no cook the repository runs.
3. **Two time-temperature points only**, not a series; no rate is extractable. The 120 C arms are the
   only ones with any bearing on the roadmap.
4. **The mechanism is inferred, not labelled.** Which aminoketone pairs actually formed each pyrazine
   is Shu's proposal; the paper has no isotope work. In particular the trialkylpyrazines are explained
   by two aminoketones (one carrying the extra carbon from an aldol-derived dicarbonyl), NOT by the
   dihydropyrazine + aldehyde addition of Cerny 1994 / Adams 2008; both give the same products and
   this paper cannot tell them apart.
5. **Trace products not quantified**: 2,5- and 2,6-dimethylpyrazine from serine are stated as trace
   but absent from Table 1, so the serine column is a floor on pyrazine diversity.
6. **Extract concentration only for the 120 C arms** (0.5 mL under N2) — pyrazine (bp 115 C) and
   methylpyrazine losses are plausible, which would bias the 120 C parent-pyrazine numbers low relative
   to the 300 C arms.
7. **Registry SMILES check (mine, RDKit on the host):** `compounds.yml` gives `2_ethyl_3_5_dimethylpyrazine`
   the SMILES `CCc1nc(C)cnc1C`, which is 2-ethyl-3,6-dimethylpyrazine (InChIKey WHMWOHBXYIZFPF);
   2-ethyl-3,5-dimethylpyrazine is `CCc1ncc(C)nc1C` (JZBCTZLGKSYRSF), and the entry's InChI string is a
   C7 species. Shu's two isomers are distinct rows in Table 1; whoever keys them must fix the registry
   first. Not edited here.

## 6. What the aldehyde-addition rule would look like (mine, on Shu's evidence)

**In words.** A dihydropyrazine formed from two aminoketones (the trunk's aminoketone pool: aminoacetone
from MGO + amine, aminoacetaldehyde from glyoxal + amine) can, instead of oxidising to the pyrazine,
add an aldehyde RCHO at a ring CH and dehydrate to a pyrazine carrying an extra CH2R group. From
aminoacetone's own dihydropyrazine (3,6-dimethyl-2,5-dihydropyrazine) formaldehyde gives
trimethylpyrazine, acetaldehyde gives 2-ethyl-3,6-dimethylpyrazine (major) or, from the cross
dihydropyrazine with 2-aminopropanal, 2-ethyl-3,5-dimethylpyrazine (minor, ~1/5 to 1/10 of the
3,6-isomer here and in Cerny 1994); from aminoacetaldehyde's dihydropyrazine acetaldehyde gives
ethylpyrazine and formaldehyde methylpyrazine. Shu's serine arm supplies acetaldehyde and
glycolaldehyde from the amino acid itself, so ethylpyrazine appears without any sugar; his threonine
arm supplies formaldehyde and hydroxyacetone, so trimethyl- appears. The rule needs an aldehyde pool
that includes sugar-fragment aldehydes (formaldehyde, acetaldehyde, glycolaldehyde), not only the
Strecker aldehyde of the amino acid in the pot.

**Positive control (should fire):** `CC1=NCC(C)=NC1` (3,6-dimethyl-2,5-dihydropyrazine) + `CC=O`
(acetaldehyde) -> `CCc1nc(C)cnc1C` (2-ethyl-3,6-dimethylpyrazine) + H2O; and `CC1=NCC(C)=NC1` + `C=O`
-> `Cc1cnc(C)c(C)n1` (trimethylpyrazine) + H2O.

**Negative control (must not fire):** the aromatic pyrazine has no enolisable ring CH2, so
`Cc1cnc(C)cn1` (2,5-dimethylpyrazine) + `CC=O` -> no ethyldimethylpyrazine; and an aldehyde with no
alpha-hydrogen requirement is irrelevant here, but a KETONE must not add: `CC1=NCC(C)=NC1` + `CC(C)=O`
(acetone) -> no isopropyl-dimethylpyrazine (none is reported by Shu, Cerny or Adams).

**Calibration Shu offers the rule:** at 120 C the extra-alkyl products are 0.03-0.9 of the parent
dimethylpyrazine (trimethyl 0.88; 2-ethyl-3,6-dimethyl 0.035 from threonine), i.e. the addition
competes with oxidation only when the aldehyde is abundant; the isomer split 3,6 : 3,5 >= 5 : 1 says
the homo-dihydropyrazine (two aminoacetones) dominates the cross (aminoacetone + 2-aminopropanal),
consistent with B18 carrying a single aminoketone per dicarbonyl.
