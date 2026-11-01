# Ruan, Hu, Su, Li, Yin, Han, Liu & Cai 2025 — EXTRACTION
DOI: 10.1016/j.foodres.2025.117512 (printed p. 1; Crossref 2026-11-01: "Interactive effect of various fatty acids and the Maillard reaction of sesame meal enzymatic hydrolysates on aroma generation", Food Research International 221, article 117512, print 2025-12. **Matches the brief.**) Not to be confused with `ruan2014_extraction.md` (unrelated).
### Sesame-meal enzymatic hydrolysate (SMH, 1.5 g in 10 mL 0.2 M phosphate, pH 6.5) heated alone or with 2.0 g of ONE free fatty acid (C16:0, C18:0, C20:0, C18:1, C18:2, C18:3), 50 mL sealed "high-pressure bottle", oil bath **160 C, 40 min, one time point**; aqueous phase separated from the fatty-acid phase by centrifugation and analysed alone: total monosaccharides (after acid hydrolysis), 17 free amino acids, A294 / A420, **final pH (measured)**, DPPH, fatty-acid loss (GC-FID), 65 aroma-active volatiles (SAFE, GC-O-MS, single internal standard 4-nonanol). n = 3. **No thiols, no time or temperature series, no oil-phase volatiles, no volatiles for the fatty-acid-only pots.**

**Read 2026-11-01 from `data/articles/ruan2025.pdf`** (15 pp., born-digital). Prose from `pdftotext -layout`;
**Table 1 (pp. 8-10) checked cell by cell against 170 dpi renders** (scratchpad `ruan2025_work/p8-08.png`, `p9-09.png`,
`p10-10.png`; the text layer is correct but scrambled in column order); **Fig. 1A/B read from 400 dpi crops of p. 5**
(`fig1a-05.png`, `fig1b-05.png`; read error +-0.3 mg/g, +-0.2 mg/mL); **Fig. 2 and Fig. 4A-D read from the 170 dpi
render of pp. 5 and 7** (read error +-0.01 pH, +-0.005 absorbance, +-0.5 % DPPH, +-0.7 % loss). Fig. 3 (correlation
matrix), Fig. 5 (class % and heatmap) and Fig. 6 (sensory) not transcribed (not model quantities).

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** School of Food Science and Technology, Henan University of Technology, Zhengzhou (corresponding
author Wen-ting Yin, yin.wenting@haut.edu.cn; Hu B.-b. first co-writer). The same group's sesame-oil / lipid x Maillard
series (Hu et al. 2024 FRI 186:114397, lysine-glucose + sesame oil; Hu et al. 2025 Food Chem 492:144985; Yin et al.
2021, 2023, 2024 [C]) are cited but **none is keyed in `data/lit/extraction_dossiers/`** (grep "sesame", "Henan
University of Technology", 2026-11-01). **Not in `data/keys/laboratories.yml`: NEW**, proposed id `haut_yin_w`
(strings e.g. "Henan Univ. of Technology (Yin W-t)"). [!] The ledgers already carry "Zhengzhou Univ. of Light
Industry / Henan Univ. of Technology (An, Bai)" (li2026, lipid.yml / aroma_other.yml): a different group (extrusion,
An / Bai), same city and one shared institution name; keep apart (laboratories.yml rule: group, not institution).

## 1. What this paper answers for the model

| model need | what Ruan 2025 gives | tag |
|---|---|---|
| CF-031: lipid effect on MFT / FFT | **nothing**: no thiol reported (SAFE + full-scan GC-MS; sulfur volatiles are 0.22 % of the SMH total); only thiazoles / thiophenes at 0.01-0.15 mg/kg | [M] |
| CF-033: lipid-type volatiles in lipid-free pots | defatted-SMH control: hexanal, nonanal, (E)-2-heptenal, 2-pentylfuran, 1-penten-3-ol all **ND**; they appear only in unsaturated-FA arms | [M] |
| CF-088: aldehyde fall when Maillard reactants meet lipid | **no**: the fatty-acid-only pots were heated but their volatiles are not reported (only FA loss, Fig. 2) | — |
| lipid-vs-no-lipid pyrazines (cross-products / competition) | total aroma-active pyrazines 27.74 (SMH) -> 10.40-13.86 mg/kg in **every** FA arm, saturated included | [M] |
| sugar / amino-acid consumption with lipid | total monosaccharides 68.70 -> 28.93 (SMH) -> 28.6-21.7 mg/g; free AA 35.69 -> 28.74 -> 28.3-26.0 mg/mL | [M] / [F] |
| pH | **measured** after heating: 5.93 (SMH) vs 5.51-5.67 (FA arms); initial 6.5 (buffer label) | [F] / [M] |
| fatty-acid oxidation rate | **not usable**: "loss" 12-45 % in 40 min, saturated ≈ unsaturated; O2-impossible in a closed bottle (sec. 6) | [F] / [Z] |
| browning | A294 (1000x dilution), A420 (200x) per arm, single time | [F] |
| absolute concentrations | volatiles: single IS (4-nonanol), no response factors stated = semi-quantitative; sugars HPAEC-PAD and amino acids by analyser (external standards implied, not stated) | [M] / [!] |

**Answers.**
1. **CF-031 (thiols):** the paper gives nothing on MFT, FFT or any thiol (confidence ~98 %). The detected S-volatiles
   (thiazole, 2-methylthiazole, 2-acetylthiazole, dihydro-3(2H)-thiophenone, 2-thiophenemethanol) sum to 0.05-0.15
   mg/kg, near the table's 0.01 mg/kg resolution. Precursor cysteine is tiny (raw ~0.4 mg/mL [F, coarse], -75 % on
   heating [M]); methionine ~1.6 mg/mL raw [F, coarse].
2. **CF-033:** a weak negative control on the "lipid volatiles in lipid-free pots" side (confidence 70 % that it is
   a real absence rather than an LOD artefact): Soxhlet-defatted SMH heated at 160 C shows none of hexanal, nonanal,
   (E)-2-heptenal, 2-pentylfuran; the saturated-FA arms show none either. Detection limit not stated (smallest printed
   cell 0.01 mg/kg), so it cannot exclude the µg/kg levels Wang 2012 / Xu 2013 report.
3. **CF-088:** cannot be addressed (confidence ~95 %): no lipid-alone volatile slate.
4. **Pyrazine fall with lipid** (confidence 75 % that it is mostly a measurement-phase effect, not chemistry): the fall
   is as large with stearic / arachidic acid (which barely oxidise) as with linolenic acid, and its size grows with
   pyrazine lipophilicity, which is the signature of partition into the 2 g fatty-acid phase that was discarded before
   SAFE (sec. 5, [Z]). The authors read it as carbonyl competition for Maillard precursors [Q].
5. **Sugar / amino-acid consumption:** fatty acids add 0.5-10 percentage points of sugar loss and 1.4-7.7 points of
   amino-acid loss to the 58 % / 19 % from heating SMH alone [M]. The increments are small and confounded (sec. 6:
   lower pH, mass-basis dilution, ionised fatty acid in the aqueous phase). Confidence the increment is a real
   lipid-driven Maillard acceleration: ~50 %.

## 2. Conditions

- **SMH preparation** [M]: Zhuzhi No. 22 sesame (2024 harvest), ground, Soxhlet-defatted (petroleum ether), 60 mesh;
  10 g / 100 mL water at 50 C; sequential cellulase (12 U/mg, pH 4.8, 50 C, 2 h), xylanase (150 U/mg, pH 5.0, 50 C,
  2 h), alcalase (10 U/mg, pH 9.0, 50 C, 2 h), flavorzyme (0.8 U/mg, pH 6.0, 60 C, 2 h); 90 C / 10 min inactivation
  after each step; supernatant (4850 g) freeze-dried. **Composition of SMH** (raw, unheated) [M]: total
  monosaccharides 68.70 mg/g (glucose 49.5 %, galactose 26.5 %, mannose 11.2 %, xylose 6.0 %, arabinose 5.8 %,
  rhamnose 1.0 %); 17 free amino acids 35.69 mg/mL (no proline listed [!]; arginine largest ~5.8 mg/mL [F]).
- **Pot** [M]: 1.5 g SMH in 10 mL 0.2 M sodium phosphate pH 6.5 +- 0.01 (= 150 g/L SMH); + 2.0 g fatty acid (equal
  MASS, not equimolar: 6.4 mmol C20:0 to 7.8 mmol C16:0, ~7.1 mmol C18 acids [Z]); fatty-acid-only controls: 2.0 g FA in
  10 mL buffer. Oil : aqueous ≈ 0.22 v/v [Z] (density ~0.9). Chosen to "mimic the lipid-to-protein ratio in sesame
  seeds" [V].
- **Heating** [M]: 50 mL high-pressure bottles, DF-101S oil bath, 160 +- 1 C, 40 min, ice-water quench; heat-up
  included in the 40 min [!]. Headspace ~38 mL, gas not stated (air presumed [Z]) -> ~0.33 mmol O2 [Z].
- **pH** [M]: initial 6.5 (buffer set before SMH and FA addition [!]: the mixture's starting pH is not reported);
  final pH of the aqueous phase measured (precision meter; temperature not stated).
- **Work-up** [M]: centrifuge 3100 g, 10 min; **aqueous phase only** for sugars, amino acids, browning, pH, DPPH and
  volatiles; FA phase for GC-FID of fatty acids.
- **Replicates** [M]: triplicate (n = 3), mean +- SD, ANOVA + Duncan (p < 0.05). Not stated whether n = 3 are
  independent heatings or replicate injections [!].

## 3. Analytical methods

- **Monosaccharides** [M]: aqueous phase freeze-dried; 5 mg in 72 % H2SO4 (25 C, 45 min), diluted, 105 C / 2.5 h,
  neutralised; HPAEC-PAD (CarboPac PA-20). [!] This is **total hydrolysable carbohydrate** (free + oligo / poly),
  not free reducing sugar, and it is **per g of freeze-dried aqueous phase** (includes ~0.27 g phosphate salts per pot
  [Z] and any dissolved fatty acid / soap). Calibration not described.
- **Free amino acids** [M]: 10 % sulfosalicylic acid 4 : 1, 4 C 60 min, 15,450 g, 5-fold dilution; Sykam S-433D
  cation exchange (post-column derivatisation implied). Units mg/mL; whether referred to the reaction solution after
  the dilutions is not stated [!].
- **Fatty acids** [M]: BF3-methanol methylation of the FA phase, BHT as "internal standard" [!] (BHT is not a FAME; an
  unusual IS), GC-FID DB-FastFAME. "Thermal loss rate" formula not given.
- **Browning** [M]: A420 at 200-fold and A294 at 1000-fold dilution of the aqueous phase.
- **Volatiles** [M]: 20 g aqueous phase + 75 mL CH2Cl2 + 150 µL 4-nonanol (1 mg/mL = 150 µg IS), shaken 10 h at 25 C,
  SAFE 50 C, Vigreux to 5 mL, Bemelmans to 1 mL; HP-5MS, full scan 33-350 m/z; GC-O by three panellists;
  identification by NIST17, RI and odour vs standards. **Quantitation by "the internal standard method"** [V]
  (citing Yin et al. 2023): no response factors or calibration curves reported -> **semi-quantitative, single IS**.
  Units mg/kg of aqueous phase. LOD not stated.
- [!] **Retention indices on HP-5MS are wrong for many identifications**: e.g. 2-ethyl-6-methylpyrazine 761
  (literature ~1000), 2,3-dimethylpyrazine 824 (~920), pyridine 1211 (~750), benzeneacetaldehyde 1569 (~1045),
  dihydro-3(2H)-thiophenone 1521 (~950), 2-methylthiazole 1118 (~800), nonanal 1194 (~1104) (literature values from
  general knowledge, not from a file on disk [Z]). Either a column / RI-series error or misidentification; identities
  rest on MS + GC-O + standards. Treat individual identities, especially the low-level ones, as uncertain.

## 4. Quantitative results

### 4.1 Sugars and amino acids (Fig. 1, p. 5)

Totals; printed values [M], the rest [F] (+-0.3 mg/g, +-0.2 mg/mL). Consumption % is relative to SMH (raw).

| arm | total monosaccharides, mg/g | consumed, % | glucose, mg/g | free AA, mg/mL | consumed, % |
|---|---|---|---|---|---|
| SMH raw (unheated) | **68.70** [M] | — | 33.9 (49.5 % [M]) | **35.69** [M] | — |
| SMH | **28.93** [M] | **57.89** [M] | 11.9 (-64.94 % [M]) | **28.74** [M] | 19.43 [M] (19.5 % in text) |
| SMH-C20:0 | 28.6 | **58.41** [M] | 11.9 | 28.3 | **20.86** [M] |
| SMH-C18:0 | 28.0 | 59.2 [Z] | 12.2 | 28.1 | 21.3 [Z] |
| SMH-C16:0 | 28.0 | 59.2 [Z] | 12.0 | 27.1 | 24.1 [Z] |
| SMH-C18:1 | 26.2 | 61.9 [Z] | 11.3 | 27.2 | 23.8 [Z] |
| SMH-C18:2 | 25.9 | 62.3 [Z] | 11.0 | 26.6 | 25.5 [Z] |
| SMH-C18:3 | 21.7 | **68.38** [M] | 9.1 | 26.0 | **27.16** [M] |

Duncan letters [F]: sugars a / b / b / c / c / d / d / e (raw, SMH, C20:0, C18:0, C16:0, C18:1, C18:2, C18:3); amino
acids a / b / c / c / d / d / e / f. Individual sugars after heating SMH [F]: galactose 18.3 -> 5.9 (-70.54 % [M]),
mannose 7.7 -> 4.5, arabinose 3.9 -> 3.0, xylose 4.2 -> 2.9 mg/g. Amino acids [F, coarse]: cysteine raw ~0.4 mg/mL,
-75.0 % on heating [M]; lysine raw ~1.2 -> ~0.8 mg/mL (-37.9 % [M]); methionine raw ~1.6, SMH ~1.1, SMH-C18:3 ~1.2
mg/mL. Molar order of magnitude [Z, basis caveats sec. 3]: free AA ~0.27 M (mean MW ~130); carbohydrate ~12 mg/mL
≈ 70 mM hexose-equivalents if the dry aqueous phase is SMH + buffer salts.

### 4.2 Fatty-acid "thermal loss rate", % (Fig. 2, p. 5)

| | C20:0 | C18:0 | C16:0 | C18:1 | C18:2 | C18:3 |
|---|---|---|---|---|---|---|
| text, "fatty acid-only" [M] | 11.80 | 29.55 | 33.66 | 30.77 | 32.64 | 34.22 |
| text, "fatty acid-SMH" [M] | 16.98 | 33.43 | 43.84 | 36.60 | 42.92 | 45.28 |
| Fig. 2 blue = "Fatty acid-only" [F] | 17.0 | 33.4 | 43.8 | 36.6 | 42.9 | 45.3 |
| Fig. 2 red = "Fatty acid-SMH" [F] | 11.8 | 29.5 | 33.7 | 30.8 | 32.6 | 34.2 |

[!] **Text and figure swap the two series**: the text's "fatty acid-only" numbers are the figure's red (FA-SMH) curve
and vice versa. The text's words ("significantly decreased" with SMH, "suppressed fatty acid thermal loss") agree
with the **figure** (FA-SMH lower). Read: SMH lowers the loss by 5-11 points [F].

### 4.3 Final pH, browning, DPPH (Fig. 4, p. 7) [F]

| arm | final pH | A294 (x1000 dil.) | A420 (x200 dil.) | DPPH scavenging, % |
|---|---|---|---|---|
| SMH | 5.93 (a) | 0.307 (e) | 0.338 (e) | 87.3 (e) |
| SMH-C20:0 | 5.67 (b) | 0.324 (d) | 0.344 (e) | 90.5 (d) |
| SMH-C18:0 | 5.64 (c) | 0.347 (c) | 0.356 (d) | 91.4 (cd) |
| SMH-C16:0 | 5.53 (d) | 0.340 (c) | 0.373 (c) | 92.5 (bc) |
| SMH-C18:1 | 5.56 (d) | 0.338 (c) | 0.391 (b) | 92.2 (bc) |
| SMH-C18:2 | 5.55 (d) | 0.360 (b) | 0.404 (b) | 93.2 (ab) |
| SMH-C18:3 | 5.51 (e) | 0.404 (a) | 0.418 (a) | 94.1 (a) |

Initial pH 6.5 (buffer) [M]. DPPH increase 3.70-7.68 % [M] is relative (3.2-6.8 points [Z]). Undiluted A420 of SMH
≈ 68 [Z] (very dark).

### 4.4 Aroma-active volatiles, mg/kg aqueous phase, mean +- SD, n = 3 (Table 1, pp. 8-10) [M]

Column order SMH / C20:0 / C18:0 / C16:0 / C18:1 / C18:2 / C18:3. ND = not detected. Duncan letters omitted.

**Pyrazines (lipid-vs-no-lipid and cross-products)**

| compound (printed RI) | SMH | C20:0 | C18:0 | C16:0 | C18:1 | C18:2 | C18:3 |
|---|---|---|---|---|---|---|---|
| pyrazine (863) | 1.48 ± 0.04 | 1.57 ± 0.09 | 1.11 ± 0.13 | 1.32 ± 0.14 | 1.44 ± 0.05 | 1.45 ± 0.02 | 1.44 ± 0.08 |
| methylpyrazine (915) | 5.42 ± 0.22 | 3.99 ± 0.16 | 2.97 ± 0.34 | 3.52 ± 0.04 | 3.86 ± 0.18 | 3.87 ± 0.01 | 3.76 ± 0.17 |
| 2,5-dimethylpyrazine (991) | 13.64 ± 0.35 | 5.49 ± 0.22 | 4.10 ± 0.51 | 4.44 ± 0.44 | 5.01 ± 0.25 | 4.81 ± 0.03 | 4.99 ± 0.02 |
| 2,3-dimethylpyrazine (824 [!]) | 0.61 ± 0.02 | 0.14 ± 0.02 | 0.10 ± 0.02 | 0.14 ± 0.03 | 0.15 ± 0.01 | 0.16 ± 0.04 | 0.17 ± 0.01 |
| ethenylpyrazine (1112) | ND | 0.02 ± 0.00 | 0.01 ± 0.00 | ND | ND | 0.02 ± 0.00 | 0.02 ± 0.00 |
| 2-ethyl-6-methylpyrazine (761 [!]) | 0.57 ± 0.01 | 0.75 ± 0.02 | 0.64 ± 0.10 | 0.17 ± 0.00 | 0.69 ± 0.02 | 0.56 ± 0.14 | 0.64 ± 0.14 |
| 2-ethyl-5-methylpyrazine (1065) | 2.50 ± 0.07 | 0.77 ± 0.03 | 0.60 ± 0.09 | 0.71 ± 0.07 | 0.70 ± 0.03 | 0.64 ± 0.01 | ND |
| 2,3,5-trimethylpyrazine (1070) | 2.44 ± 0.08 | 0.59 ± 0.01 | 0.42 ± 0.04 | 0.33 ± 0.05 | 0.53 ± 0.03 | 0.45 ± 0.05 | 0.46 ± 0.04 |
| 2-ethenyl-6-methylpyrazine (1134) | 0.04 ± 0.01 | 0.03 ± 0.00 | 0.02 ± 0.00 | ND | 0.05 ± 0.00 | 0.05 ± 0.00 | 0.05 ± 0.00 |
| acetylpyrazine (1175) | ND | ND | ND | ND | 0.13 ± 0.00 | ND | ND |
| 3-ethyl-2,5-dimethylpyrazine (1128) | 0.07 ± 0.00 | 0.08 ± 0.01 | 0.06 ± 0.00 | 0.07 ± 0.01 | 0.08 ± 0.00 | 0.07 ± 0.00 | 0.08 ± 0.01 |
| 2,3-dimethyl-5-ethylpyrazine (1144) | 0.28 ± 0.00 | ND | ND | ND | ND | ND | ND |
| 2,6-diethylpyrazine (1145) | 0.06 ± 0.01 | ND | 0.05 ± 0.00 | ND | ND | 0.08 ± 0.02 | 0.08 ± 0.00 |
| 2-butyl-3-methylpyrazine (1440) | 0.63 ± 0.15 | 0.44 ± 0.00 | 0.32 ± 0.00 | 0.22 ± 0.00 | 0.39 ± 0.00 | 0.40 ± 0.01 | ND |
| **sum of pyrazines** | **27.74 ± 0.74** | 13.86 ± 0.5 | 10.40 ± 0.15 | 10.92 ± 1.06 | 13.03 ± 0.31 | 12.56 ± 0.63 | 11.69 ± 0.55 |

Pyridines: 2-butylpyridine 0.10 ± 0.00 in C18:2 only (a lipid x Maillard cross-product, from 2,4-decadienal + NH3
[C]); pyridine 0.03-0.04 everywhere; 2-methylpyridine 0.01 (C20:0, C18:1); 1-(2-pyridinyl)ethanone 0.02 (SMH), 0.01
(C18:3). Pyrroles sum 0.81-1.14 (no lipid trend).

**Lipid-oxidation volatiles**

| compound (printed RI) | SMH | C20:0 | C18:0 | C16:0 | C18:1 | C18:2 | C18:3 |
|---|---|---|---|---|---|---|---|
| hexanal (842) | ND | ND | ND | ND | 0.04 ± 0.01 | 0.26 ± 0.02 | 0.07 ± 0.00 |
| nonanal (1194) | ND | ND | ND | ND | 0.14 ± 0.03 | ND | ND |
| (E)-2-heptenal (996) | ND | ND | ND | ND | ND | 0.04 ± 0.00 | ND |
| 2-pentylfuran (1042) | ND | ND | ND | ND | 0.05 ± 0.00 | 0.33 ± 0.01 | 0.07 ± 0.00 |
| 1-penten-3-ol (751) | ND | ND | ND | ND | ND | ND | 1.49 ± 0.07 |
| pentanoic acid (1149) | ND | ND | ND | ND | ND | 0.93 ± 0.05 | 0.31 ± 0.00 |
| 2-heptanone (895) | 0.05 ± 0.00 | ND | ND | ND | ND | ND | 0.05 ± 0.00 |
| 2-ethyl-5-methylfuran (1001) | ND | ND | ND | ND | 0.12 ± 0.00 | ND | 0.13 ± 0.02 |

Missing from the slate: (E,E)-2,4-decadienal, (E)-2-octenal, (E)-2-nonenal, 1-octen-3-ol, octanal, propanal,
(E,E)-2,4-heptadienal (the linolenate markers). Only GC-O-active compounds are listed ("aroma-active"), so absence
from Table 1 is not proof of absence [!].

**Sulfur heterocycles** (no thiols)

| compound (printed RI) | SMH | C20:0 | C18:0 | C16:0 | C18:1 | C18:2 | C18:3 |
|---|---|---|---|---|---|---|---|
| thiazole (879) | 0.03 ± 0.00 | 0.02 ± 0.00 | 0.03 ± 0.00 | 0.03 ± 0.00 | 0.06 ± 0.00 | 0.06 ± 0.00 | 0.06 ± 0.01 |
| 2-methylthiazole (1118 [!]) | ND | ND | 0.02 ± 0.00 | ND | ND | ND | 0.03 ± 0.00 |
| 2-acetylthiazole (1149) | 0.02 ± 0.01 | ND | 0.02 ± 0.00 | ND | ND | ND | ND |
| dihydro-3(2H)-thiophenone (1521 [!]) | 0.09 ± 0.00 | 0.07 ± 0.00 | 0.05 ± 0.01 | 0.05 ± 0.01 | 0.07 ± 0.01 | 0.06 ± 0.00 | 0.08 ± 0.00 |
| 2-thiophenemethanol (1406) | 0.01 ± 0.00 | ND | 0.05 ± 0.00 | ND | 0.04 ± 0.01 | 0.05 ± 0.00 | 0.07 ± 0.00 |

**Furans and Strecker aldehydes (trunk)**

| compound (printed RI) | SMH | C20:0 | C18:0 | C16:0 | C18:1 | C18:2 | C18:3 |
|---|---|---|---|---|---|---|---|
| furfural (984) | 0.27 ± 0.01 | 0.40 ± 0.02 | 0.35 ± 0.07 | 0.35 ± 0.04 | 0.37 ± 0.01 | 0.33 ± 0.01 | 0.46 ± 0.03 |
| 5-methylfurfural (1001) | 0.14 ± 0.02 | 0.30 ± 0.01 | 0.24 ± 0.04 | 0.36 ± 0.05 | 0.27 ± 0.01 | 0.24 ± 0.01 | 0.26 ± 0.02 |
| 5-HMF (1294) | ND | 0.65 ± 0.11 | 0.64 ± 0.02 | 0.19 ± 0.00 | 0.68 ± 0.10 | 0.68 ± 0.05 | 0.67 ± 0.15 |
| 2-furanmethanol (905) | 19.34 ± 0.00 | 18.48 ± 0.66 | 13.21 ± 0.96 | 12.97 ± 0.56 | 17.85 ± 0.25 | 18.27 ± 0.21 | 19.70 ± 0.97 |
| 2-acetylfuran (998) | ND | 0.12 ± 0.01 | 0.12 ± 0.00 | ND | ND | 0.10 ± 0.00 | 0.13 ± 0.00 |
| "furanone" (1150) | ND | ND | 1.76 ± 0.44 | 2.36 ± 0.40 | ND | ND | ND |
| 4-hydroxy-5-methyl-3(2H)-furanone (1122) | 0.64 ± 0.13 | 1.18 ± 0.17 | 1.08 ± 0.26 | 0.26 ± 0.00 | 1.73 ± 0.36 | 2.20 ± 0.01 | 0.80 ± 0.00 |
| maltol (1278) | 4.05 ± 0.16 | 4.86 ± 0.29 | 3.45 ± 0.04 | 2.94 ± 0.19 | 4.12 ± 0.75 | 4.34 ± 0.03 | ND |
| 2-methylpropanal (<700) | ND | 0.72 ± 0.18 | 0.09 ± 0.00 | 0.21 ± 0.00 | 1.69 ± 0.00 | 0.34 ± 0.00 | 0.43 ± 0.03 |
| 2-methylbutanal (700) | 1.19 ± 0.07 | 0.58 ± 0.04 | 0.54 ± 0.09 | 0.48 ± 0.07 | 0.47 ± 0.04 | 0.46 ± 0.01 | 0.58 ± 0.01 |
| benzaldehyde (997) | 0.19 ± 0.01 | 0.12 ± 0.01 | 0.10 ± 0.02 | 0.09 ± 0.01 | 0.10 ± 0.01 | 0.09 ± 0.00 | 0.11 ± 0.01 |
| benzeneacetaldehyde (1569 [!]) | 1.92 ± 0.11 | 3.05 ± 0.13 | 2.32 ± 0.29 | 2.32 ± 0.01 | 2.81 ± 0.10 | 2.44 ± 0.07 | 4.25 ± 0.13 |
| 3-hydroxy-2-butanone (acetoin; "3-hydroxy-2-oxobutane", 775) | 5.91 ± 0.21 | 8.27 ± 0.33 | 5.92 ± 0.73 | 6.74 ± 0.73 | 8.43 ± 0.04 | 8.51 ± 0.06 | 9.32 ± 0.05 |
| 2,3-pentanedione (739) | 0.23 ± 0.03 | 0.16 ± 0.01 | 0.12 ± 0.02 | 0.11 ± 0.03 | 0.13 ± 0.02 | 0.17 ± 0.02 | 0.15 ± 0.02 |
| 3-methylbutanoic acid (937) | 2.22 ± 0.35 | 2.58 ± 0.61 | 1.68 ± 0.21 | 1.49 ± 0.23 | 2.43 ± 0.32 | 2.10 ± 0.13 | 2.80 ± 0.32 |

Not transcribed (in Table 1, minor for the model): pyrroles, phenolics (2-methoxyphenol 0.84 SMH -> 0.34-0.46),
dihydro-2-methyl-3(2H)-furanone, 2(5H)-furanone, 2-ethyl-5-methylfuran, 2-methylfuran, 1,3-benzodioxol-5-ol,
apocynin, 1-hydroxy-2-butanone, 2-butanone, cyclotene, 2,3-butanediol, esters, benzeneacetic acid. Table 2 (ROAVs) not
transcribed (computed from Table 1 and Van Gemert thresholds; several thresholds implausible, e.g. 2-furanmethanol
0.0125 mg/kg [!]).

## 5. Mechanism evidence and its kind

- **Within-study factor, single endpoint**: lipid yes / no and lipid type (chain length, unsaturation) at one T, one t.
  No labelling, no time course, no dose series.
- **The pyrazine fall is partition-shaped** [Z]. Ratio (FA arm / SMH), range over the six arms: pyrazine 0.75-1.06,
  methylpyrazine 0.55-0.74, 2,5-DMP 0.30-0.40, 2-ethyl-5-methylpyrazine 0.24-0.31 (ND in C18:3), 2,3,5-TMP 0.14-0.24,
  2,3-DMP 0.16-0.28. A water / fatty-acid partition with phase ratio phi ≈ 0.22 and K ≈ K_ow gives a water-phase
  retention 1 / (1 + phi K): ~0.9 (pyrazine, log K_ow ≈ -0.3), ~0.75 (methylpyrazine, ≈ 0.2), ~0.5 (2,5-DMP, ≈ 0.6),
  ~0.2-0.3 (trimethyl / ethylmethyl, ≈ 1.0-1.5) — log K_ow values approximate, from general knowledge, not from a file
  on disk. Same order and similar size as observed, and **independent of unsaturation** (stearic and arachidic acid
  cut 2,5-DMP as much as linolenic). Exceptions: 2,3-DMP falls more than 2,5-DMP, and 2-ethyl-6-methylpyrazine (RI
  761, identity doubtful) does not fall. Probability the bulk of the pyrazine fall is phase partition plus the lower
  pH, not carbonyl competition: ~0.7.
- **The authors' "fatty acids raise alkylpyrazines"** (2-ethenyl-6-methyl-, 2,6-diethyl-, 3-ethyl-2,5-dimethyl-) [Q]
  rests on differences of 0.01-0.03 mg/kg at the table's resolution; not a usable cross-product signal.
- **Lipid-specific products** (present only with the matching acid): hexanal / 2-pentylfuran / (E)-2-heptenal /
  pentanoic acid / 2-butylpyridine maximal with C18:2; 1-penten-3-ol only with C18:3; nonanal only with C18:1. These
  are the expected parent-acid assignments [C] and so are qualitative identity checks of the lipid slates, at very low
  levels (0.04-1.5 mg/kg) because the aqueous phase alone was sampled.
- **Strecker aldehydes**: 2-methylbutanal halves with any FA (partition-compatible: log K_ow ~1.2 [Z]), while
  2-methylpropanal appears only with FA (0.09-1.69) and benzeneacetaldehyde rises (x1.2-2.2). 2-methylpropanal ND in
  SMH but 2-methylbutanal 1.19 is odd for a hydrolysate rich in Val and Ile [!].

## 6. Flags

1. **Aqueous phase only** [!]: every volatile, sugar and amino-acid number is for the aqueous phase after the ~2.2 mL
   fatty-acid phase was removed by centrifugation. Lipophilic volatiles (lipid aldehydes, alkylpyrazines, any thiol
   or disulfide) are under-reported by an unknown, compound-specific factor; lipid-vs-no-lipid contrasts of these are
   confounded by partition (the same confound as CF-031's oil-phase explanation and Zhang 2022's water-phase data).
2. **Fatty-acid "loss" cannot be oxidation** [Z]: ~38 mL headspace of air holds ~0.33 mmol O2; one O2 per hydroperoxide
   limits oxidation to ≤ ~5 mol % of the ~7 mmol FA, yet 12-45 % "loss" is reported, and palmitic (43.8 %) ≈ linoleic
   (42.9 %), against the known oxidisability order (linolenic >> linoleic > oleic >> saturated; the engine's
   `parameters_lipid.oleate_fraction` docstring takes oleate ~10x below linoleate). The "loss" is a recovery /
   methylation / partition artefact (e.g. ionised FA staying in the aqueous phase at pH 6.5 in 0.2 M phosphate).
   **Never key it as an oxidation rate.**
3. **Fig. 2 text vs figure swap** of the FA-only and FA-SMH series (sec. 4.2).
4. **Monosaccharides = total hydrolysable carbohydrate per g dry aqueous phase**, not free reducing sugar (H2SO4
   hydrolysis); the dry-mass basis shifts between arms (dissolved soap adds mass and dilutes mg/g), which by itself
   could produce part of the 0.5-10 point "extra consumption" [Z].
5. **pH confound**: final pH 5.93 (SMH) vs 5.51-5.67 (FA arms), and the starting pH of the FA + buffer mixtures is not
   reported. Part of the drop is plausibly plain acid-base (7 mmol carboxylic acid vs 2 mmol phosphate) rather than
   oxidation acids [Z]. pH measurement temperature not stated.
6. **Equal mass, not equimolar** fatty-acid charges (6.4-7.8 mmol).
7. **Single internal standard, no response factors** (4-nonanol): mg/kg are semi-quantitative; between-compound
   comparisons are uncalibrated; within-compound arm contrasts are usable if partition is accounted for.
8. **RI column inconsistent with HP-5MS** for at least seven compounds (sec. 3); identities of low-level compounds
   uncertain.
9. Heat-up inside the 40 min in a 50 mL bottle in an oil bath (no product T trace); headspace gas not stated.
10. Text says "2,6-dimethylpyrazine (ROAV = 2.07)" where Table 2 has 2,6-diethylpyrazine [!] (typo).
11. Replicate structure (independent heatings or injections) not stated; many SDs are 0.00 at two decimals.

## 7. Bearing on open questions

- **CF-031** (lipid effect on MFT / FFT) — LEAVES UNCHANGED: no thiols. Indirect: the pyrazine data show how large a
  pure aqueous-sampling partition effect can be with a 0.22 v/v free-fatty-acid phase (x0.3-0.4 for 2,5-DMP even with
  stearic acid), which supports keeping the "physical partition" explanation alive for every water-phase thiol datum
  (P that this changes CF-031's status: ~5 %).
- **CF-033** (lipid volatiles in lipid-free pots) — weak counter-datum: defatted SMH at 160 C, ND for hexanal,
  nonanal, (E)-2-heptenal, 2-pentylfuran, LOD unstated (P of changing the entry: ~10 %).
- **CF-088** (aldehyde fall with Maillard reactants) — LEAVES UNCHANGED: no lipid-alone volatile slate.
- **CF-011 / phosphate**: 0.2 M phosphate in every arm; not varied.
- **Lipid x Maillard pyrazine coupling (lipid lane cross-products)**: does not support a "lipid lowers pyrazines by
  carbonyl competition" term; if anything it argues for a phase-partition measurement operator before any chemical
  term (P ~0.7, sec. 5).
- **Trunk**: one-point (160 C, 40 min) sugar / amino-acid consumption in a complex hydrolysate with measured final pH;
  usable at most as a coarse matrix check, not a fit row (carbohydrate is total hydrolysable, composition undefined).

## 8. Ledger-ready measurements

- Table 1 lipid-oxidation volatiles x 7 arms (lipid.yml; aqueous phase only, single IS).
- Table 1 pyrazines x 7 arms, lipid vs no lipid (lipid.yml; partition confound).
- Table 1 sulfur heterocycles x 7 arms (sulfur.yml; no thiols).
- Table 1 furans / Strecker aldehydes / acetoin x 7 arms (trunk.yml).
- Fig. 1 total carbohydrate and free amino acids x 8 bars (trunk.yml).
- Fig. 4 final pH, A294, A420 x 7 arms (trunk.yml).
- Fig. 2 fatty-acid "loss" with / without SMH (lipid.yml) — key only with the flag that it is not oxidation.
- Not keyable: any rate constant, Ea, thiol, oil-phase concentration, or the FA-only pots' volatiles (not reported).
