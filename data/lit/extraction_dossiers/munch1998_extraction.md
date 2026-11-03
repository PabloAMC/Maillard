# Münch & Schieberle 1998 — EXTRACTION (file stem `munch1998`)
DOI: 10.1021/jf980511t (printed on p. 4695 and as the page-4701 footer "JF980511T"). Crossref 2026-11-03: "Quantitative
Studies on the Formation of Key Odorants in Thermally Treated Yeast Extracts Using Stable Isotope Dilution Assays",
J. Agric. Food Chem. **46**(11):4695-4701, authors Petra Münch, Peter Schieberle, published 1998-10-07. **Matches the
PDF and the brief.** Received 15 May 1998, accepted 20 Aug 1998.
**Companion:** `munch1997_extraction.md` (the AEDA paper on the same extracts; amino-acid tables reused here).
[!] **Abbreviations differ from Hofmann 2000:** here **PA = phenylalanine, PAA = phenylacetaldehyde**; in
`hofmann2000_extraction.md` PA = phenylacetaldehyde, PAA = phenylacetic acid. Key by name, never by these letters.

### 18 odorants by SIDA in heated commercial yeast extract (CYE), self-prepared autolysate (SPYE) and two LMW (< 1000 Da) baker's-yeast fractions (I fresh, II NaCl-fermented), all 1 + 3 w/w in water, 145 °C / 20 min, laboratory autoclave; free amino acids of the fermentation variants; CYE + glucose or 2-oxopropanal (amounts not printed); and **defined Strecker pots at 145 °C in 0.1 M phosphate pH 7.0: Leu / Phe / Met 1 mmol + glucose or 2-oxopropanal 3 mmol (20 min), and a Phe + 2-oxopropanal time course 20-150 min with Phe loss.** No sugars, thiamine, nucleotides or pH of any extract; no MFT routes; no replicates of the heatings.

**Read 2026-11-03 from `data/articles/munch1998.pdf`** (7 pp., ACS typeset, text layer present). Prose from
`pdftotext -layout`; **Tables 1-10 read by eye from 150 dpi page renders** (pp. 3-6, scratchpad
`ingest/munch1998_work/p-3..6.png`; text layer agreed with the renders). Figs. 1-3 are structures / synthesis /
mechanism, no data. Every OAV and mol % re-computed by `awk` [Z].

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Institut für Lebensmittelchemie der TU München **and** Deutsche Forschungsanstalt für
Lebensmittelchemie, Lichtenbergstrasse 4, Garching; corresponding author **Peter Schieberle**; DFG grant Schi 399/2-2.
Existing id **`dfa_hofmann_schieberle`**; proposed lab string **"DFA Garching (Hofmann / Schieberle)"** (listed; trunk's
Hofmann 2000 entries use the equivalent "DFA Munich (Hofmann / Schieberle)"). Same lab, same SIDA standards and
autoclave as `hofmann2000_extraction.md` (Hofmann, **Münch** & Schieberle 2000) — **a second system and a second
temperature, not a second laboratory.** Other same-lab dossiers: as listed in `munch1997_extraction.md`.

## 1. What this paper answers for the model

| model need | what Münch 1998 gives | tag |
|---|---|---|
| Strecker aldehyde yield from fed dicarbonyl, defined pot | Table 9: Leu / Phe / Met + MGO or glucose, 145 °C / 20 min, µg and mol % | [M] |
| Strecker time course + amino-acid loss | Table 10: Phe + MGO, 20 / 30 / 60 / 90 / 120 / 150 min, PAA µg, µmol, % Phe degraded | [M] |
| measured dicarbonyl trajectory (forced input) | **none** — MGO is dosed, never measured | [!] |
| odorant amounts in heated yeast extracts | Tables 3-5 (18 / 18 / 15 odorants, µg/kg **dry** material, SIDA) | [M] |
| MFT, FFT in yeast extract | CYE 530 / 580, SPYE 22 / 29 µg/kg DW; LMW-YE FFT < 0.1; **MFT not reported for LMW-YE** | [M] |
| MFT route shares (thiamine vs ribose + Cys vs C2 + C3) | **nothing** (no thiamine, ribose, nucleotide, labelling) | [!] |
| precursor composition | 6 amino acids × 4 extracts (Table 7), fermentation variants (Table 6); Cys in `munch1997` only | [M] |
| dicarbonyl / sugar effect in a real matrix | Table 8: CYE + glucose or + MGO (doses **not printed**) | [M] / [!] |
| pyrazines | 2,3-diethyl-5-methylpyrazine only (SIDA) | [M] |
| furanones | HDMF, sotolon (SIDA) | [M] |
| kinetics / T series | one temperature (145 °C); time only in Table 10 | [M] |

**Answers.**
1. **System** (confidence ~90 %). Same materials as `munch1997`: CYE (German manufacturer; baker's yeast from Fa.
   Wieninger, Passau), SPYE (90 g DW yeast autolysed 24 h at 50 °C with NaCl 5 g + ethyl acetate 4.5 g per L, pH
   5.6 → 5.2, freeze-dried), LMW-YE (18 g DW yeast, cell mill 75 s under liquid CO2, 0.1 M phosphate pH 7.0,
   ultrafiltrate < 1000 Da, freeze-dried) from fresh yeast or after 12 h at 32 °C under pure O2 in a nutrient medium
   with/without glucose 100 g/L and/or NaCl 35 g/L (four fermentations, Table 6). **Heating: 5-50 g freeze-dried
   material in water 1 + 3 by weight, laboratory autoclave (Roth type II), 20 min, 145 °C**; closed; heat-up,
   headspace and pH not reported. **Composition printed: six free amino acids only** (Orn, Pro, Leu, Ile, Met, Phe;
   Tables 6-7); no sugars, ribose, nucleotides, thiamine, cysteine (Cys is in `munch1997` Table 3), water activity or
   pH. **No temperature series; one time series (Table 10, model pot only).**
2. **SIDA odorants** (confidence ~97 %). 18 labelled-standard assays (Table 1): MFT and FFT (as 4-vinylpyridine
   thioethers, [2H] standards), 3-methylbutanal, 2-methylbutanal (via [2H3]-3-methylbutanal, response factor 0.54 —
   an isomer standard, not an isotopologue), methional, phenylacetaldehyde, sotolon, HDMF, 2,3-butanedione,
   3-mercapto-2-pentanone, 2,3-diethyl-5-methylpyrazine, 2-AP, 2-acetyl-2-thiazoline, 2-/3-methylbutanoic acid,
   phenylacetic acid, 2-phenylethanol, 3-methylbutanol. **Not quantified: bis(2-methyl-3-furyl) disulfide,
   methanethiol, H2S, other pyrazines, furfural, 2-methylpropanal.** Amounts in sec. 4 (µg/kg dry material; µg per
   30 mL pot in Tables 9-10).
3. **Precursor / spiking experiments** (confidence ~80 %).
   (a) **Strecker from MGO vs glucose, defined pots (Table 9)**: at equal molar dose (3 mmol), MGO raises the
   aldehyde ×10.8 (Leu → 3-methylbutanal, 0.4 → 4.3 mol %), ×3.7 (Phe, 0.3 → 1.0) and ×3.3 (Met, 299 → 989 µg) [M/Z];
   donor order Leu ≫ Phe ≈ Met on MGO. (b) **Table 10**: Phe loss decelerates (average pseudo-first-order k
   0.0186 min⁻¹ at 20 min → 0.0080 at 150 min [Z]) while the PAA yield on Phe consumed rises from 3.3 to 7.6 mol % at
   90 min, then falls to 6.2 at 150 min; PAA amount plateaus at 41.6-43.7 µmol from 90 min [M]. (c) **In CYE (Table 8)**
   glucose raises 3-methylbutanal ×2.35, methional ×1.56, PAA ×2.44, 2-AP ×5.6, and MGO raises them ×6.45 / ×26.5 /
   ×14.8 [Z]; **both additions lower MFT (×4.3, ×3.4) and FFT (×3.6, ×4.8)** [Z] — directional only, doses unknown.
   (d) **Fermentation**: NaCl raises the six amino acids ×4 (Met) to ×25 (Pro); glucose + NaCl raises Orn ×14; Orn ×14
   → 2-AP ×6.6, Leu ×5 → 3-methylbutanal ×8.0 [M/Z]. **Nothing on thiamine vs ribose + Cys for MFT, no labelling.**
4. **Usability** (calibrated).
   - **Defined Strecker pots (Tables 9, 10): yes, a same-lab cross-SYSTEM / cross-TEMPERATURE check for E3**
     (P(useful as a check row set) ~0.6; P(changes a fitted constant) ~0.15). Hofmann 2000 (`strecker_hofmann2000*`)
     is at 98 °C, 0.5 M phosphate, 100 + 100 mM; Münch 1998 is at **145 °C, 0.1 M phosphate, 33 mM amino acid + 100 mM
     dicarbonyl or glucose, sealed autoclave**. Contrasts [Z]: Phe + MGO, Phe degraded 18 % at 30 min (98 °C) vs 31 %
     at 20 min (145 °C): an apparent Ea of only 19-28 kJ/mol from two points, which says the 145 °C pot is MGO- or
     heat-up-limited, not a clean Arrhenius pair; phenylacetaldehyde at the same nominal [MGO]0 (0.1 M) is 9.5 µmol per
     mmol Phe at 98 °C / 30 min vs 10.1 at 145 °C / 20 min. Glucose pots: 3-methylbutanal 0.089 mol % (98 °C, 1:1, 30 min,
     Hofmann) vs 0.40 mol % (145 °C, 3:1, 20 min); phenylacetaldehyde 0.042 vs 0.27 mol %. **Missing to simulate**: MGO
     trajectory (not measured), heat-up profile of the autoclave, final pH, Leu / Phe / Ile species in the engine
     (only Met → methional exists, B22), the unidentified Phe-loss method (as in Hofmann 2000 Table 3). Forced-input
     in the strict sense (measured dicarbonyl) **no** (~0.95); a dosed-dicarbonyl check like Hofmann 2000 Table 3,
     **yes**.
   - **Yeast-extract pots (Tables 3-5, 7-8): not a forced-input pot** (P(simulable) ~0.05): no reducing sugars,
     ribose / nucleotides, thiamine, dicarbonyls, pH, buffer; amino acids for six species only; amounts on a dry-weight
     basis of different materials. Usable as **apparent molar conversions** of amino acid → Strecker aldehyde in a real
     matrix (CYE 0.016-0.037 %, LMW-YE 0.06-0.20 %; sec. 4i) — an order-of-magnitude plausibility band for the plant-
     based target, nothing tighter (P(helps calibrate) ~0.1).
   - **MFT route shares (E6): no** (~0.97). The CYE MFT (530 µg/kg DW) is the strongest yeast-extract MFT level in the
     corpus by SIDA, but with no thiamine or ribose data it cannot apportion routes. Directional: added glucose / MGO
     lowers both thiols in CYE.
   - **Pyrazines (E13/E14): no** — 2,3-diethyl-5-methylpyrazine is not an engine species and its precursors are not
     given.
5. **Ledger entries**: sec. 7 — six (two trunk, two sulfur, two aroma_other).

## 2. Conditions

- **Extract heating** [M]: freeze-dried material 5-50 g "suspended in water (1 + 3 by weight)" [V] (25 % w/w solids;
  `munch1997` used tap water 1:5 w/v); laboratory autoclave type II (Roth), **20 min, 145 °C**. Labelled standards
  5-10 µg each added after cooling to aliquots; ether extraction, high-vacuum sublimation. Water type (tap or
  distilled) not stated here.
- **Defined Strecker pots** [M]: amino acid 1 mmol + glucose 3 mmol (or 2-oxopropanal 3 mmol) in **30 mL phosphate
  buffer 0.1 mol/L, pH 7.0**, autoclave, 145 °C, 20 min (Table 9 fn). Table 10: L-Phe 1 mmol + 2-oxopropanal 3 mmol,
  same buffer, 20-150 min (fn says "20 min" [!]). Concentrations [Z]: amino acid 33 mM, dicarbonyl / glucose 100 mM.
- **CYE + glucose / 2-oxopropanal (Table 8)** [M]: doses, pot size **not printed** [!].
- **Fermentation (Table 6)** [M]: 30 g DW (100 g wet) yeast per L tap water with glucose 100 g, NaCl 35 g, ZnSO4
  28.2 g, (NH4)H2PO4 10 g, KH2PO4 4 g, CuSO4 1.35 g, Na citrate 1.14 g, MgSO4 0.4 g, pyridoxal 12 mg, meso-inositol
  100 mg, Ca pantothenate 5 mg; pH 5.0; 12 h at 32 °C under pure O2; four variants omitting glucose, NaCl or both
  (control A = neither). [!] Differs from the 1997 medium (MgSO4 0.8 g; thiamine 44 mg, biotin 200 mg and the Fe salt
  listed there are absent here).
- **pH**: initial pH 7.0 for the model pots; extracts not measured; no final pH anywhere.
- **Replicates** [M]: heatings apparently single (no ±). Method repeatability: a model mixture of all odorants +
  standards worked up 3×, values within ±5 % [V/Q]. Table 6 amino acids "mean values of triplicates" (no SD printed).

## 3. Analytical methods

- **SIDA / mass chromatography** [M]: HRGC (Fisons 5160; FFAP, SE-54, DB-Wax) to an ion-trap ITD 800 in CI
  (methanol); selected ions and response factors in Table 1 (0.54-1.00). Absolute quantification against labelled
  standards ([2H] or [13C2]), except **2-methylbutanal (quantified on [2H3]-3-methylbutanal, RF 0.54)** and
  **2-methylbutanoic acid (on [2H2-3]-3-methylbutanoic acid)** — isomer standards [!].
- **Thiols** [M]: extract spiked with [2H3]-MFT and [2H2]-FFT (10 µg; MFT stock reduced with LiAlH4 to break the
  disulfide, then calibrated on unlabelled FFT), pH 7.5, 4-vinylpyridine (2 g), 2 h dark, excess quenched with Cys
  3 g; pyridylethyl thioethers isolated by acid / base partitioning; m/z 220 / 222 / 223.
- **OAV** = concentration / odour threshold in water (Table 2, from Rychlik et al. 1998 [C]); [!] the OAVs divide
  **µg per kg dry material** by µg/L of water, i.e. they are not OAVs of the heated solution (÷ 4 for the 1 + 3 pot).
- **Amino acids**: as in `munch1997` (ninhydrin analyser). **Phe loss in Table 10: method not stated** [!].

## 4. Data

### 4a. Table 1 (p. 4697): SIDA ions, standards, response factors [M]

| odorant | m/z | standard | m/z | RF | phase |
|---|---|---|---|---|---|
| 2-acetyl-1-pyrroline | 112 | [2H2-6]-2-AP | 114-118 (sum) | 0.98 | SE-54 |
| 2-acetyl-2-thiazoline | 130 | [2H4] | 134 | 1.00 | SE-54 |
| 2,3-butanedione | 101 | [13C4] | 105 | 1.00 | DB-Wax |
| 2,3-diethyl-5-methylpyrazine | 151 | [2H3] | 154 | 0.96 | SE-54 |
| 2-furfurylthiol | 81 | [2H2] | 83 | 0.89 | SE-54 |
| sotolon | 129 | [13C2] | 131 | 1.00 | FFAP |
| HDMF | 129 | [13C2] | 131 | 1.00 | FFAP |
| methional | 105 | [2H3] | 108 | 0.98 | SE-54 |
| 3-mercapto-2-pentanone | 119 | [2H2] | 121 | 0.88 | SE-54 |
| 2-methylbutanal | 69 | [2H3]-3-methylbutanal | 71-72 | 0.54 | DB-Wax |
| 3-methylbutanal | 69 | [2H3]-3-methylbutanal | 71-72 | 0.94 | DB-Wax |
| 2-/3-methylbutanoic acid | 85 | [2H2-3]-3-methylbutanoic acid | 87-88 (sum) | 0.85 | FFAP |
| 3-methylbutanol | 71 | [2H2] | 73 | 0.95 | FFAP |
| phenylacetaldehyde | 121 | [13C2] | 123 | 1.00 | SE-54 |
| phenylacetic acid | 137 | [13C2] | 139 | 1.00 | FFAP |
| phenylethanol | 105 | [13C2] | 107 | 1.00 | FFAP |
| FFT pyridylethyl thioether | 220 | [2H2] | 222 | 1.00 | SE-54 |
| MFT pyridylethyl thioether | 220 | [2H2] (sic; d-1 is [2H3]) | 223 | 1.00 | SE-54 |

### 4b. Table 2 (p. 4697): odour thresholds in water, µg/L [M, from Rychlik 1998 C]

MFT 0.007; FFT 0.01; 3-methylbutanal 0.4; methional 1.8; sotolon 0.3; phenylacetaldehyde 4; 3-mercapto-2-pentanone
0.7; 2-methylbutanal 4; 2,3-butanedione 15; HDMF 25; 2,3-diethyl-5-methylpyrazine 0.09; 2-AP 0.1; 3-methylbutanoic
acid 750; 2-methylbutanoic acid 540; 2-acetyl-2-thiazoline 1; phenylacetic acid 1000 (triangle test, Schieberle &
Hofmann 1997b [C]); 2-phenylethanol 1000; 3-methylbutanol 1000.

### 4c. Tables 3 and 4 (p. 4698): SIDA concentrations (µg/kg) and OAVs, heated CYE and SPYE [M]

Basis: CYE per kg dry yeast extract; SPYE per kg **dry weight of the fresh baker's yeast**.

| no. | odorant | CYE µg/kg | CYE OAV | SPYE µg/kg | SPYE OAV |
|---|---|---|---|---|---|
| 1 | 2-methyl-3-furanthiol (MFT) | 530 | 75710 | 22 | 3140 |
| 2 | 2-furfurylthiol (FFT) | 580 | 58000 | 29 | 2900 |
| 3 | 3-methylbutanal | 1456 | 3640 | 617 | 1543 |
| 4 | methional | 1095 | 608 | 287 | 159 |
| 5 | sotolon | 117 | 390 | 64 | 213 |
| 6 | phenylacetaldehyde | 1250 | 313 | 1 | 0.3 |
| 7 | 3-mercapto-2-pentanone | 170 | 243 | nd (not determined) | nd |
| 8 | 2-methylbutanal | 965 | 241 | 370 | 93 |
| 9 | 2,3-butanedione | 1932 | 129 | 190 | 13 |
| 10 | HDMF | 2471 | 99 | 1352 | 54 |
| 11 | 2,3-diethyl-5-methylpyrazine | 2.7 | 30 | 53 | 589 |
| 12 | 2-acetyl-1-pyrroline | 2.2 | 22 | 6 | 60 |
| 13 | 3-methylbutanoic acid | 3714 | 5 | 3043 | 4 |
| 14 | 2-methylbutanoic acid | 1486 | 3 | 1217 | 2 |
| 15 | 2-acetyl-2-thiazoline | 2.8 | 3 | 8 | 8 |
| 16 | phenylacetic acid | 2660 | 3 | 932 | 1 |
| 17 | 2-phenylethanol | 476 | 0.5 | 43500 | 43 |
| 18 | 3-methylbutanol | 55 | 0.1 | 5535 | 6 |

Every OAV reproduces as concentration / Table 2 threshold to the printed rounding [Z].

### 4d. Table 5 (p. 4698): LMW-YE from fresh yeast (I) and NaCl + nutrients 12 h yeast (II), µg/kg (OAV) [M]

| odorant | I | II | OAV I | OAV II |
|---|---|---|---|---|
| FFT | <0.1 | <0.1 | <10 | <10 |
| 3-methylbutanal | 100 | 801 | 250 | 2000 |
| methional | 92 | 225 | 51 | 125 |
| phenylacetaldehyde | 81 | 351 | 20 | 88 |
| 2-methylbutanal | 84 | 527 | 21 | 132 |
| 2,3-butanedione | 467 | 1314 | 31 | 87 |
| HDMF | 25 | 316 | 1 | 13 |
| 2,3-diethyl-5-methylpyrazine | 1 | 2 | 11 | 24 |
| 2-acetyl-1-pyrroline | 5 | 33 | 50 | 330 |
| 3-methylbutanoic acid | 5915 | 5714 | 8 | 8 |
| 2-methylbutanoic acid | 2366 | 2286 | 4 | 4 |
| 2-acetyl-2-thiazoline | 13 | 2 | 13 | 2 |
| phenylacetic acid | 180 | 165 | 0.2 | 0.2 |
| 2-phenylethanol | 2470 | 3737 | 3 | 4 |
| 3-methylbutanol | 4556 | 1921 | 5 | 2 |

"15 of the 18": **MFT, sotolon and 3-mercapto-2-pentanone not reported** for LMW-YE [!]. OAV rounding: II butanedione
87 (88 [Z]), II DEMP 24 (22 [Z]), I phenylethanol 3 (2.5 [Z]).

### 4e. Table 6 (p. 4699): free amino acids (mg/g DW) of the LMW fraction after four fermentations, n = 3 means [M]

| amino acid | A control | B (A + glucose) | C (A + NaCl) | D (A + glucose + NaCl) |
|---|---|---|---|---|
| proline | 0.04 | 0.05 | 1.00 | 1.52 |
| ornithine | 0.21 | 0.21 | 0.80 | 2.98 |
| leucine | 0.25 | 0.29 | 1.31 | 1.24 |
| isoleucine | 0.20 | 0.19 | 1.24 | 0.90 |
| methionine | 0.07 | 0.07 | 0.27 | 0.16 |
| phenylalanine | 0.19 | 0.19 | 1.02 | 0.82 |

### 4f. Table 7 (p. 4700): precursor amino acids (mg/kg DW) and, in parentheses, the odorant formed (mg/kg DW) [M]

| amino acid (odorant) | CYE | SPYE | LMW-YE-I | LMW-YE-II |
|---|---|---|---|---|
| ornithine | 2654 | 2719 | 210 | 2980 |
| proline | 3679 | 4004 | 40 | 1520 |
| (2-AP, from Orn + Pro) | (2 × 10⁻³) | (6 × 10⁻³) | (5 × 10⁻³) | (33 × 10⁻³) |
| leucine (3-methylbutanal) | 14003 (1.46) | 10280 (0.62) | 250 (0.1) | 1240 (0.80) |
| isoleucine (2-methylbutanal) | 8748 (0.97) | 6466 (0.37) | 200 (0.08) | 900 (0.53) |
| methionine (methional) | 4280 (1.1) | 2534 (0.29) | 70 (0.09) | 160 (0.23) |
| phenylalanine (phenylacetaldehyde) | 9364 (1.25) | 5447 (1 × 10⁻³) | 190 (0.08) | 820 (0.35) |

Parenthesised values = Tables 3-5 odorant µg/kg ÷ 1000 [Z, every cell checked]; footnote b "for explanation see
text" — they are **mass**, not molar, amounts. LMW-YE-I / II = Table 6 A / D × 1000 [Z].

### 4g. Table 8 (p. 4700): CYE heated alone, + glucose (A), + 2-oxopropanal (B); µg/kg (OAV) [M]

| odorant | control | A | B |
|---|---|---|---|
| 3-methylbutanal | 1456 (3640) | 3418 (8545) | 9391 (23477) |
| methional | 1095 (608) | 1707 (948) | 29068 (16149) |
| phenylacetaldehyde | 1250 (313) | 3056 (764) | 18510 (4627) |
| 2-acetyl-1-pyrroline | 2.2 (22) | 12.3 (123) | nd |
| MFT | 530 (75714) | 122 (17428) | 158 (22571) |
| FFT | 580 (58000) | 163 (16300) | 122 (12200) |

Control = Table 3 (MFT OAV printed 75710 there, 75714 here; 530 / 0.007 = 75714 [Z]). Doses of glucose and
2-oxopropanal: **not printed**.

### 4h. Table 9 (p. 4700): Strecker aldehyde from amino acid 1 mmol + glucose 3 mmol, or + 2-oxopropanal 3 mmol (b), 30 mL 0.1 M phosphate pH 7.0, 145 °C, 20 min [M]

| precursor | aldehyde | µg | mol % printed | mol % [Z] |
|---|---|---|---|---|
| L-leucine | 3-methylbutanal | 343 | 0.4 | 0.398 |
| L-leucine (b) | 3-methylbutanal | 3693 | 4.3 | 4.29 |
| L-phenylalanine | phenylacetaldehyde | 325 | 0.3 | 0.270 |
| L-phenylalanine (b) | phenylacetaldehyde | 1217 | 1.0 | 1.01 |
| L-methionine | methional | 299 | 0.3 | 0.287 |
| L-methionine (b) | methional | 989 | **0.7** | **0.949** |

[!] **Met + MGO row: 989 µg is 0.95 mol %, not 0.7** (MW 104.17); every other row reproduces. One of the two printed
numbers is wrong (P(µg right) ~0.5). mol % is on the amino acid charged.

### 4i. Table 10 (p. 4700): phenylacetaldehyde (PAA) from L-Phe 1 mmol + 2-oxopropanal 3 mmol, 30 mL 0.1 M phosphate pH 7.0, 145 °C [M]

| time (min) | PAA µg | PAA µmol | Phe degraded % | PAA mol % of Phe degraded | [Z] PAA mol % of Phe charged | [Z] mean k, −ln(1−f)/t (min⁻¹) |
|---|---|---|---|---|---|---|
| 20 | 1217 | 10.1 | 31 | 3.3 | 1.01 | 0.0186 |
| 30 | 2216 | 18.5 | 40 | 4.6 | 1.85 | 0.0170 |
| 60 | 2875 | 24.0 | 47 | 5.1 | 2.40 | 0.0106 |
| 90 | 4992 | 41.6 | 55 | 7.6 | 4.16 | 0.0089 |
| 120 | 5215 | 43.5 | 62 | 7.0 | 4.35 | 0.0081 |
| 150 | 5241 | 43.7 | 70 | 6.2 | 4.37 | 0.0080 |

The µmol column uses MW ≈ 120 (5215 / 120.15 = 43.40) [Z]. The 20 min row is the Table 9 Phe + MGO pot (same 1217 µg).
Pot concentrations [Z]: Phe0 33.3 mM, MGO0 100 mM; PAA 0.34-1.46 mM.

**Same-lab comparison, Hofmann 2000 Table 3** (98 °C, 0.5 M phosphate, Phe 0.1 M + MGO 0.1 M, closed vial, from
`hofmann2000_extraction.md`): phenylacetaldehyde 6.2 / 9.5 / 10.6 / 12.3 µmol per mmol Phe at 10 / 30 / 60 / 120 min;
Phe degraded 9 / 18 / 36 / 49 % [C]. At 145 °C the aldehyde reaches ×3.2 the 98 °C maximum (43.7 vs 13.5 µmol/mmol)
[Z], and does not fall within 150 min.

### 4j. Derived apparent conversions in the extracts (mol odorant / mol free amino acid, %) [Z]

| precursor → odorant | CYE | SPYE | LMW-YE-I | LMW-YE-II |
|---|---|---|---|---|
| Leu → 3-methylbutanal | 0.016 | 0.009 | 0.061 | 0.098 |
| Ile → 2-methylbutanal | 0.017 | 0.009 | 0.064 | 0.089 |
| Met → methional | 0.037 | 0.016 | 0.188 | 0.201 |
| Phe → phenylacetaldehyde | 0.018 | 0.00003 | 0.059 | 0.059 |
| Orn + Pro → 2-AP | 0.00004 | 0.0001 | 0.0023 | 0.0008 |
| Cys → MFT / FFT (Cys from `munch1997` Table 3) | 0.026 / 0.029 | 0.011 / 0.015 (SPYA batch) | — | — |

Both sides per kg dry material as printed; assumes the free amino acid is the only precursor (peptides ignored).
Compare the defined glucose pots at 145 °C (Table 9): 0.27-0.40 %, i.e. ×7-25 above CYE and ×1.3-7 above the LMW-YEs.

## 5. Flags and internal consistency

1. [!] **Table 9 Met + MGO: 989 µg vs 0.7 mol %** (computes 0.95). Text "methional being increased 28-fold" (Table 8 B)
   computes ×26.5.
2. [!] **Table 10 footnote says "reacted for 20 min"** for a 20-150 min series; Phe-loss method never described.
3. [!] **LMW-YE-I is "fresh bakers' yeast" (Table 5 title, text) but "identical with A ... in Table 6"** [V], and A is the
   12 h control fermentation (nutrients without glucose or NaCl). And `munch1997` Table 7 "fresh" LMW-YE has different
   values (e.g. Pro 0.54, Orn 0.74, Leu 0.17, Met 0.76 mg/g vs 0.04 / 0.21 / 0.25 / 0.07 here). P(I = control
   fermentation, title wrong) ~0.6.
4. [!] **Reused composition, changed basis**: Table 7 CYE amino acids = `munch1997` Table 9 "mg reacted" per 100 g
   × 10 (e.g. Leu 1400 mg → 14003 mg/kg; 0.980-0.983 × the 1997 mg/g), while SPYE values (Orn 2719, Leu 10280,
   Met 2534) are ×1.7-3.5 the 1997 SPYA mg/g — a different batch or basis, not stated.
5. [!] **Water ratio 1 + 3 w/w** here vs **1:5 w/v** in `munch1997`: the AEDA and SIDA pots are not the same pot.
6. [!] OAVs are computed per kg **dry** material against water thresholds; the heated solution is ×4 more dilute.
   Table 3 / 4 footnote b points to "Table 1" for thresholds (they are in Table 2).
7. [!] 2-methylbutanal and 2-methylbutanoic acid quantified on 3-methyl isomer standards (RF 0.54 for the aldehyde).
8. [!] Table 1 lists the MFT thioether standard as [2H2] (m/z 223, one more than the FFT [2H2] at 222); d-1 is
   [2H3]-MFT per the methods — likely a typo in the table.
9. [!] No heating replicates; the ±5 % is work-up repeatability of a model mixture.
10. Text checks that pass [Z]: MFT / FFT / butanedione SPYE vs CYE "20- or 10-fold lower" (24 / 20 / 10.2); Orn ×14,
    2-AP ×6.6, Leu ×5.0, 3-methylbutanal ×8.0; NaCl ×3.9 (Met) to ×25 (Pro); 2-AP / 3-methylbutanal / HDMF ×6.6 / 8.0 /
    12.6 ("6-12"); Leu MGO / glucose ×10.8 ("11-fold"); Table 10 mol % on Phe degraded all reproduce.

## 6. Existing records in the repo that cite this paper

`grep -rIn` of `data/ docs/ src/ scripts/ tasks/` for `jf980511t`, `munch`/`münch`/`muench`, title words and the
distinctive numbers (1217, 3693, 5241, 29068, 75710) on 2026-11-03: **no ledger, fit, registry or code record cites or
carries this paper's data.** Two dossiers cite its threshold table second-hand:

| record | claim | paper | verdict |
|---|---|---|---|
| `dreher2003_extraction.md` l. 82 | MFT threshold 6.14 × 10⁻⁸ mM "[C: Munch & Schieberle 1998]" | Table 2: 0.007 µg/L = 6.13 × 10⁻⁸ mM (MW 114.17) [Z] | **match** (Münch takes it from Rychlik 1998) |
| `gao2020b_extraction.md` l. 200 | "Münch & Schieberle 1998 [MFT 0.007, FFT 0.01 µg/L]" | Table 2: 0.007 and 0.01 | **match** |

Adjacent, not a citation: `data/lit/benchmark_intake_registry.json` `uspto_ptacts_2023_yeast_extract_anchor` carries a
yeast-extract `mft_oav_band` [2000, 8000] with `analytical_method: sida`, sourced to US 9,943,096 B2. This paper's
SIDA OAVs (per kg dry, flag 6) are 3140 (SPYE, inside the band) and 75710 (CYE, ×9.5 above its top): the band is not
a ceiling for yeast extracts. The engine's thresholds (`src/kinetic_core/species_sulfur.py`: MFT 0.005, FFT 0.006 µg/L)
are not attributed to this paper and differ from its 0.007 / 0.01.

## 7. Proposed ledger entries (YAML, NOT written into the ledgers)

```yaml
  # trunk.yml (Strecker aldehydes from fed dicarbonyls)
  - id: munch1998_tab9_strecker_glucose_vs_mgo_145c
    dossier: munch1998_extraction.md
    paper: "Muench & Schieberle 1998 J Agric Food Chem 46:4695"
    lab: "DFA Garching (Hofmann / Schieberle)"
    system: "L-leucine, L-phenylalanine or L-methionine 1 mmol + glucose 3 mmol or 2-oxopropanal 3 mmol in 30 mL (33 + 100 mM)"
    medium: "phosphate buffer 0.1 mol/L pH 7.0, autoclave (closed), air not controlled"
    temperature_c: [145]
    ph: {initial: [7.0], measured_during: false, final: null}
    time_min: [20]
    species: [3-methylbutanal, phenylacetaldehyde, methional]
    kind: endpoint
    location: "Table 9 (dossier sec. 4h)"
    status: digitised_in_dossier
    quant_class: sida
    simulable: true
    condition_quote: "The amino acid (1 mmol) and glucose (3 mmol) were dissolved in phosphate buffer (30 mL; 0.1 mol/L; pH 7.0) and reacted for 20 min at 145 °C in an autoclave."
    flags: [ph_start_only, no_replicates, heat_up_not_reported, "met_mgo_row_ug_vs_molpct_inconsistent", engine_species_missing]
    note: "6 cells; donor ratios on MGO (Leu 4.3 >> Phe 1.0 ~ Met 0.7-0.95 mol %) and MGO/glucose ratios (x10.8 / x3.7 / x3.3) are within-study contrasts. Same lab as hofmann2000 (98 C, 0.5 M phosphate, 1:1): a second temperature and buffer, not a second lab. Leu / Phe not engine species; Met -> methional is (B22). Phe + MGO row = Table 10 at 20 min."
  - id: munch1998_tab10_phe_mgo_timecourse_145c
    dossier: munch1998_extraction.md
    paper: "Muench & Schieberle 1998 J Agric Food Chem 46:4695"
    lab: "DFA Garching (Hofmann / Schieberle)"
    system: "L-phenylalanine 1 mmol + 2-oxopropanal 3 mmol in 30 mL (33 + 100 mM)"
    medium: "phosphate buffer 0.1 mol/L pH 7.0, autoclave (closed), air not controlled"
    temperature_c: [145]
    ph: {initial: [7.0], measured_during: false, final: null}
    time_min: [20, 30, 60, 90, 120, 150]
    species: [phenylacetaldehyde, phenylalanine]
    kind: time_series
    location: "Table 10 (dossier sec. 4i)"
    status: digitised_in_dossier
    quant_class: other
    simulable: true
    condition_quote: "L-Phenylalanine (1 mmol) and 2-oxopropanal (3 mmol) were dissolved in phosphate buffer (30 mL; pH 7.0; 0.1 mol/L)"
    flags: [ph_start_only, no_replicates, phe_method_not_stated, footnote_time_contradicts_table, heat_up_not_reported, mgo_not_measured]
    note: "6 times; PAA by SIDA ([13C2]), Phe loss method unstated (hence other). PAA plateaus ~43.5 umol from 90 min while Phe keeps falling (55 -> 70 %). 145 C partner of hofmann2000_table3_mgo_phe_timecourse (98 C, 0.5 M phosphate, 1:1); confounded contrast (T, buffer, Phe:MGO). PA = phenylalanine here (opposite of Hofmann 2000's PA)."
  # sulfur.yml
  - id: munch1998_tab3_4_5_yeast_extract_thiols
    dossier: munch1998_extraction.md
    paper: "Muench & Schieberle 1998 J Agric Food Chem 46:4695"
    lab: "DFA Garching (Hofmann / Schieberle)"
    system: "commercial yeast extract, self-prepared autolysate, LMW (<1000 Da) baker's-yeast fractions I / II; 1 + 3 w/w in water"
    medium: "aqueous slurry, unbuffered, pH not stated; laboratory autoclave"
    temperature_c: [145]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [20]
    species: ["MFT", "FFT", "3-mercapto-2-pentanone"]
    kind: endpoint
    location: "Tables 3, 4, 5 (dossier secs. 4c, 4d)"
    status: digitised_in_dossier
    quant_class: sida
    simulable: false
    condition_quote: "The freeze-dried material (5-50 g) was suspended in water (1 + 3 by weight) and thermally treated for 20 min at 145 °C"
    flags: [matrix_food, ph_not_stated, composition_incomplete, dry_weight_basis, no_replicates]
    note: "ug/kg dry material (CYE: per kg extract; SPYE / LMW: per kg DW yeast). CYE MFT 530 / FFT 580; SPYE 22 / 29; LMW FFT < 0.1, MFT not reported. No thiamine, ribose, nucleotides: no route shares. Cys only in munch1997 Table 3."
  - id: munch1998_tab8_cye_glucose_mgo_addition
    dossier: munch1998_extraction.md
    paper: "Muench & Schieberle 1998 J Agric Food Chem 46:4695"
    lab: "DFA Garching (Hofmann / Schieberle)"
    system: "commercial yeast extract alone, + glucose, or + 2-oxopropanal (doses not printed)"
    medium: "aqueous slurry (1 + 3 w/w presumably), pH not stated; laboratory autoclave"
    temperature_c: [145]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [20]
    species: ["MFT", "FFT", 3-methylbutanal, methional, phenylacetaldehyde, 2-acetyl-1-pyrroline]
    kind: ratio
    location: "Table 8 (dossier sec. 4g)"
    status: digitised_in_dossier
    quant_class: sida
    simulable: false
    condition_quote: ""
    flags: [matrix_food, dose_not_stated, ph_not_stated, "also_relevant: aroma_other"]
    note: "Directional only: glucose and MGO both lower MFT (x4.3, x3.4) and FFT (x3.6, x4.8) and raise the Strecker aldehydes (MGO: methional x26.5). One table, sulfur and Strecker rows together."
  # aroma_other.yml
  - id: munch1998_tab3_4_5_yeast_extract_odorants
    dossier: munch1998_extraction.md
    paper: "Muench & Schieberle 1998 J Agric Food Chem 46:4695"
    lab: "DFA Garching (Hofmann / Schieberle)"
    system: "commercial yeast extract, self-prepared autolysate, LMW (<1000 Da) baker's-yeast fractions I / II; 1 + 3 w/w in water"
    medium: "aqueous slurry, unbuffered, pH not stated; laboratory autoclave"
    temperature_c: [145]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [20]
    species: [3-methylbutanal, 2-methylbutanal, methional, phenylacetaldehyde, 2-acetyl-1-pyrroline, "2,3-diethyl-5-methylpyrazine", "2,3-butanedione", furaneol, sotolon, 2-acetyl-2-thiazoline, 3-methylbutanoic acid, 2-methylbutanoic acid, phenylacetic acid, 2-phenylethanol, 3-methylbutanol]
    kind: endpoint
    location: "Tables 3, 4, 5 non-sulfur rows; Table 7 (dossier secs. 4c, 4d, 4f, 4j)"
    status: digitised_in_dossier
    quant_class: sida
    simulable: false
    condition_quote: "The freeze-dried material (5-50 g) was suspended in water (1 + 3 by weight) and thermally treated for 20 min at 145 °C"
    flags: ["complements: munch1998_tab3_4_5_yeast_extract_thiols", matrix_food, composition_incomplete, dry_weight_basis, isomer_standard_for_2_methylbutanal, "also_relevant: trunk"]
    note: "Apparent amino-acid -> Strecker aldehyde conversions 0.016-0.037 % (CYE), 0.06-0.20 % (LMW-YE) [dossier sec. 4j]: a plausibility band for matrix Strecker yields, not a simulable pot (no sugars, dicarbonyls, pH)."
  - id: munch1998_tab6_7_yeast_amino_acids_fermentation
    dossier: munch1998_extraction.md
    paper: "Muench & Schieberle 1998 J Agric Food Chem 46:4695"
    lab: "DFA Garching (Hofmann / Schieberle)"
    system: "LMW fraction of baker's yeast after 12 h at 32 C under O2 in nutrient medium +- glucose 100 g/L +- NaCl 35 g/L; plus CYE and SPYE"
    medium: "dry materials, before heating"
    temperature_c: []
    ph: {initial: null, measured_during: false, final: null}
    time_min: []
    species: [proline, ornithine, leucine, isoleucine, methionine, phenylalanine]
    kind: other
    location: "Tables 6, 7 (dossier secs. 4e, 4f)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "Based on dry weight. Data are mean values of triplicates."
    flags: [matrix_food, "lmw_ye_i_fresh_vs_control_fermentation_conflict", "cye_values_reused_from_munch1997_basis_changed"]
    note: "Precursor composition for the Table 5 / 7 correlation (Orn x14 -> 2-AP x6.6; Leu x5 -> 3-methylbutanal x8). Six amino acids only; no SD printed."
```

## 8. Supplementary material

None referenced (no Supporting Information paragraph; 1998 JAFC). Two items are unpublished: the synthesis of
[13C2]-phenylacetaldehyde and [13C2]-2-phenylethanol ("Schnermann and Schieberle, unpublished results", Fig. 2) and the
phenylacetic acid threshold (Schieberle & Hofmann 1997b, an ACS meeting abstract [C]). Neither carries model data.
The thresholds come from the Rychlik, Schieberle & Grosch 1998 DFA compilation (ISBN 3-9803426-5-4) [C], not on disk.
