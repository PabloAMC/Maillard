# Zheng, Brown, Ledig, Mussinan & Ho 1997 — EXTRACTION
### *"Formation of Sulfur-Containing Flavor Compounds from Reactions of Furaneol and Cysteine, Glutathione, Hydrogen Sulfide, and Alanine/Hydrogen Sulfide"*
### Furaneol (DMHF, 0.1 M) with an equimolar cysteine, glutathione, Na2S or alanine + Na2S, in 0.1 M phosphate set to pH 7.0 before the reagents went in. Parr bomb, 130 °C, 3 h, then SDE and GC-FID. One table of 90 entries in "mg/mol", all quantified against Furaneol with a response factor of 1. **No 2-methyl-3-furanthiol, no 2-furfurylthiol and no 2,5-dimethyl-3-furanthiol in any pot.** The ring-S swap product 2,5-dimethyl-4-hydroxy-3(2H)-thiophenone forms in all four (0.15–1.26 mg/mol). Glutathione gives a subset of the cysteine products, at about a quarter of the sulfur-volatile total. Read 2026-10-27.

**Full read of `data/articles/zheng1997.pdf`** (4 pages, the ACS print PDF with a clean text layer).
Table 1 is set in two side-by-side halves with blank cells, so the text layer cannot be trusted to
assign values to columns. **Every Table 1 cell below was assigned by reading a 300-dpi raster of
p. 896**, not the text layer. No figures. Companions: `zheng1994_extraction.md` (the same first and
last author: H2S release kinetics from cysteine and glutathione, cited here),
`shu1988_extraction.md` (DMHF + cysteine across pH 2.2/5.1/7.1, the only other fed DMHF + cysteine
experiment), and `whitfield1999_extraction.md` / `whitfield2001_extraction.md` (the C5 homologue
norfuraneol + H2S or cysteine).

**Provenance:** **[M]** measured, printed as a number · **[Q]** qualitative · **[Z]** derived by me
from the paper's own numbers, never printed · **[!]** flag.

---

## 0. IDENTITY

| field | value as printed |
|---|---|
| Authors | Yan Zheng (Rutgers and IFF), Sharon Brown, Walter O. Ledig, Cynthia Mussinan (IFF), Chi-Tang Ho* (Rutgers) |
| Title | "Formation of Sulfur-Containing Flavor Compounds from Reactions of Furaneol and Cysteine, Glutathione, Hydrogen Sulfide, and Alanine/Hydrogen Sulfide" |
| Venue | *J. Agric. Food Chem.* **1997**, 45 (3), **894–897** |
| DOI | **not printed on the PDF.** The PDF prints the ACS manuscript codes "S0021-8561(96)00624-3" (p. 894 footer) and "JF960624H" (p. 897). By the ACS convention these give 10.1021/jf960624h. **Confirmed on Crossref (api.crossref.org, 2026-10-27):** that DOI resolves to this title, JAFC 45, 894–897, March 1997 |
| Affiliations | Dept. of Food Science, Cook College, NJ Agricultural Experiment Station, Rutgers, New Brunswick NJ; R&D Center, International Flavors and Fragrances, Union Beach NJ |
| Dates | received 15 Aug 1996, accepted 30 Dec 1996; ASAP abstract 15 Feb 1997 |
| PDF character | one table (Table 1, p. 896), no figures, no structures, no replicates, no error bars. GC-FID quantification against Furaneol as the only reference |

---

## 1. THE ONE-PARAGRAPH ANSWER

Equimolar 0.1 M Furaneol (2,5-dimethyl-4-hydroxy-3(2H)-furanone, the model's `DMHF`) was heated
with 0.1 M cysteine (F/C), glutathione (F/G), alanine + Na2S (F/A/S) or Na2S alone (F/S). Each pot
was 100 mL of 0.1 M sodium phosphate pH 7.0 in a 500 mL Parr bomb, 130 °C for 3 h. Volatiles were
taken by Likens–Nickerson SDE and quantified by GC-FID against Furaneol with every relative response
factor set to 1. Sulfur-containing volatiles summed [Z] to **23.1 (F/C, 15 compounds), 6.1 (F/G, 7),
39.6 (F/A/S, 22) and 64.2 (F/S, 19) mg/mol**. The authors read this as **H2S availability being
limiting**. **No 2-methyl-3-furanthiol (MFT), no 2-furfurylthiol (FFT), no disulfide of either and
no 2,5-dimethyl-3-furanthiol (the C6 analogue of MFT) appears in any pot.** The ring-oxygen-to-sulfur
swap product **2,5-dimethyl-4-hydroxy-3(2H)-thiophenone (the model's `DMHFS`) forms in all four:
0.22 / 0.15 / 1.26 / 1.16 mg/mol**. So do 2,5-dimethylthiophene (0.31 / 0.95 / 0.70 / 0.56) and
3,5-dimethyl-1,2,4-trithiolane (1.31 / 1.88 / 2.47 / 3.74). **Thiirane (11.48) and 2-methylthiophene
(0.25) appear only with cysteine.** The authors assign them to cysteine's Strecker degradation
(mercaptoacetaldehyde), which glutathione cannot undergo because its cysteine amine is in a peptide
bond. Both mercaptopentanone isomers form **only in the two Na2S pots** (2-mercapto-3-pentanone
0.39 / 0.46; 3-mercapto-2-pentanone 0.49 / 0.31), and 2-mercapto-3-butanone peaks at 22.10 mg/mol
in F/S. **[!] The four pots almost certainly did not share pH 7.** 0.1 M of Na2S or glutathione
exceeds the capacity of 0.1 M phosphate, and post-mixing pH was not measured (§2.1). The
glutathione-versus-cysteine and sulfide-versus-amino-acid contrasts are confounded with pH.

---

## 2. METHODS [M] (p. 895)

| item | as printed | [Z] / note |
|---|---|---|
| Reactants | "Equimolar (0.1 M) Furaneol and cysteine, glutathione, sodium sulfide, and alanine or sodium sulfide" | the "alanine **or** sodium sulfide" is a typo for "and"; Table 1's F/A/S has both. Salt forms (cysteine free base or HCl; Na2S anhydrous or nonahydrate) **not stated** [!] |
| Amounts | 0.1 M each in **100 mL** | 0.01 mol (1.28 g) Furaneol per pot [Z] |
| Buffer | "**0.1 M sodium phosphate buffer solution at pH 7.0**"; reagents "mixed and dissolved in" it | pH set **before** the reagents were added; **no post-mixing, hot or post-reaction pH** [!], see §2.1 |
| Vessel | **500 mL Parr bomb** (Parr Instrument Co.), "placed in a preheated reaction vessel" | ~400 mL headspace over 100 mL liquid [Z]; atmosphere not stated (presumably air) |
| Heating | glycerol bath, **130 °C, 3 h**; cooled to room temperature under running tap water | come-up time not stated. Bath temperature, not liquid temperature |
| Storage | "stored in a refrigerator until the isolation" | duration not stated |
| Isolation | reaction mass + 100 mL water rinse into a 1000 mL flask; **Likens–Nickerson SDE into 100 mL CH2Cl2**; dried over Na2SO4; Kuderna–Danish to ~4 mL; N2 stream to **200 µL** | SDE time not stated. **Volatile H2S and part of the low boilers (thiirane, acetaldehyde) are lost or under-recovered by this workup** |
| GC | HP 5890, **FID**; 50 m × 0.32 mm i.d. OV-1 (methyl silicone), 0.52 µm; 40 → 225 °C at 2 °C/min, hold 30 min; He 27.0 cm/s; 1 µL, split 400 mL/min | |
| Identification | mass spectrum + retention index (I_E, ethyl-ester series C1–C18) against the IFF MS library and RI data bank "obtained using the authentic compounds" | GC-MS conditions not given |
| Quantification | Furaneol as reference: response factor **8.74 × 10⁷ mg⁻¹** (1500 ppm Furaneol in CH2Cl2); "Assuming the relative response factors of the other volatiles are the same as Furaneol"; mg/mol = (area × dilution factor) / (molarity of sample × response factor of Furaneol) | see §2.2 for what "mg/mol" can mean |
| Replicates | **none stated** | |

### 2.1 [Z] The pH of each pot was probably not 7 [!]

This is my estimate, not a measurement. It uses 25 °C pKa values and assumes the free-acid or
free-base reagent forms; the paper states neither. 0.1 M phosphate at pH 7.0 holds ~0.042 M H2PO4⁻
and ~0.058 M HPO4²⁻ (pKa2 ≈ 6.86). Each reagent was 0.1 M, about the size of the whole buffer:

| pot | added acid/base load | estimated starting pH |
|---|---|---|
| F/C | cysteine zwitterion (pKa 1.96 / 8.18 / 10.28) is net neutral at pH 7 | **~7** if free base; **~2** if cysteine·HCl (0.1 M HCl exhausts the 0.058 M HPO4²⁻) |
| F/G | glutathione (pKa 2.12 / 3.59 / 8.75 / 9.65) releases its Gly-COOH proton (~0.1 M) on going to pH 7 | **~3.5–4**: the HPO4²⁻ is exhausted and the remaining protons sit on the Gly carboxyl |
| F/A/S | Na2S 0.1 M (two proton acceptors) plus alanine (amine pKa 9.7) | **~10**: S²⁻/HS⁻ uses up the H2PO4⁻, and alanine's ammonium then buffers |
| F/S | Na2S 0.1 M | **~12.5–13**: S²⁻ uses up the 0.042 M H2PO4⁻ and the rest hydrolyses to HS⁻ + OH⁻ |

Under these assumptions **the glutathione pot is ~3 pH units more acidic than the cysteine pot** (if
cysteine was the free base), **and the sulfide pots are strongly alkaline**. So "sulfide gives more
S-volatiles than cysteine" and "glutathione gives fewer than cysteine" each mix a donor effect with a
pH effect. At high pH, H2S is fully ionised (HS⁻) and DMHF (enol pKa ≈ 8.5, literature value, not
from this paper) is partly ionised. Those are the opposite of the conditions Shu & Ho 1988 found
favourable for DMHFS.

### 2.2 What "mg/mol" means [!]

The table footnote says "Calculations based on Furaneol". The natural reading is **mg of compound
per mol of Furaneol charged** (0.01 mol per pot). The printed equation divides by the **molarity**
(0.1 M), not by the moles (0.01 mol), so the unit could instead be mg per (mol/L), numerically 10×
smaller than per mol charged. Which one applies depends on the undefined "dilution factor". **Ratios
between cells are unaffected. Absolute molar yields are uncertain by up to ×10**, on top of the
untested RRF = 1 assumption and the unknown SDE recovery.

---

## 3. TABLE 1 AS PRINTED [M]

"Table 1. Volatile Compounds Identified from the Reaction of Furaneol with Cysteine (F/C), Glutathione
(F/G), Alanine and Sodium Sulfide (F/A/S), and Sodium Sulfide (F/S)". Unit: amount^a (mg/mol).
Footnotes: ^a "Calculations based on Furaneol." ^b "Starting material." ^c "Coeluted with
2-ethyl-5-methylpyrazine." **Blank = not listed for that reaction** (the paper never says whether
blank means not detected or below a limit). The table's own class headings and alphabetical order
are kept; the two printed halves are merged.

### Aldehydes/Ketones

| compound | F/C | F/G | F/A/S | F/S |
|---|---:|---:|---:|---:|
| (E)-3-penten-2-one | 1.07 | 0.50 | | |
| (E)-4-hexen-3-one | | 0.12 | | |
| 1-hydroxy-2-butanone | | 0.18 | | |
| 1-hydroxy-2-propanone | | 1.48 | | |
| 2,3,4-trimethyl-2-cyclopentenone | | | 2.12 | |
| **2,3-butanedione** | **34.49** | **10.14** | **14.97** | **14.47** |
| 2,3-dimethyl-2-cyclopentenone | 0.16 | | | |
| 2,3-hexanedione | 0.17 | 0.20 | 0.35 | 0.36 |
| **2,3-pentanedione** | **0.19** | **0.27** | **0.10** | **0.20** |
| 2,4-hexanedione | | 2.07 | 0.73 | 1.47 |
| 2-butenal | | 0.49 | <0.02 | |
| 2-hexanone | 0.17 | 0.44 | | |
| 2-hydroxy-3-butanone | 0.10 | 0.12 | | |
| 2-hydroxy-3-pentanone | 1.82 | 0.68 | 0.78 | |
| 2-methylpropanal | | 0.04 | <0.02 | 0.04 |
| 2-penten-3-one | | 2.36 | 0.17 | |
| 3,4-hexanedione | 0.54 | 1.02 | | |
| 3-hexanone | 0.49 | 0.93 | | |
| 3-hydroxy-2-pentanone | 0.30 | 0.74 | | 0.17 |
| 3-methyl-2-butanone | 0.10 | 0.11 | 0.06 | 0.16 |
| 3-methyl-3-hexen-2-one | | | 0.51 | 1.30 |
| 3-pentanone | 2.44 | 2.38 | 2.52 | 1.45 |
| 4-hexen-3-one | 0.11 | | | |
| 4-penten-2-one | 1.19 | | | |
| acetaldehyde | 0.35 | 5.42 | | |
| acetone | | 0.07 | | |

### Furans

| compound | F/C | F/G | F/A/S | F/S |
|---|---:|---:|---:|---:|
| 2,4,5-trimethyl-3(2H)-furanone | 0.13 | 0.71 | | 0.62 |
| 2,5-dimethyl-3(2H)-furanone | 0.14 | 0.28 | 0.07 | 0.04 |
| **2,5-dimethyl-4-hydroxy-3(2H)-furanone^b (residual Furaneol)** | **4.36** | **7.64** | **0.23** | **3.51** |
| 4,5-dimethyl-2(3H)-dihydrofuranone | | | | 1.24 |

### Pyrazines

| compound | F/C | F/G | F/A/S | F/S |
|---|---:|---:|---:|---:|
| 2-ethyl-5-methylpyrazine | 0.40 | | | |
| trimethylpyrazine | c | | | |
| 2,5-dimethylpyrazine | 0.37 | 0.14 | | |

### Sulfides (the table's heading; it includes the thiols)

| compound | F/C | F/G | F/A/S | F/S |
|---|---:|---:|---:|---:|
| 2,5-dimethyl-4-mercapto-3(2H)-furanone | 0.11 | | 1.89 | |
| **2-mercapto-3-butanone** (≡ 3-mercapto-2-butanone) | 0.21 | | 4.63 | **22.10** |
| **2-mercapto-3-pentanone** (model `MP3P`) | | | **0.39** | **0.46** |
| **3-mercapto-2-pentanone** (model `MP2P`) | | | **0.49** | **0.31** |
| 2-methyl-2-[(1-methylethyl)-thiol]propane | | | | 0.93 |
| 2-methyl-2-pentanethiol | | | 0.11 | 0.10 |
| 2-methyl-3-(thioacetyl)furan | | | | 0.36 |
| **3,5-dimethyl-1,2,4-trithiolane** | **1.31** | **1.88** | **2.47** | **3.74** |
| 3-methyl-2-butanethiol | | | 5.44 | 15.58 |
| 4-(ethylthiol)butan-2-one | | | 0.38 | |
| 4-(methylthio)phenol | | | 2.20 | |
| dibutyl disulfide | | | 1.88 | |
| ethanethioic acid, S-methyl ester | | | 3.42 | |
| **thiirane** | **11.48** | | | |

### Thiazoles/Thiazolines

| compound | F/C | F/G | F/A/S | F/S |
|---|---:|---:|---:|---:|
| 2,4,5-trimethylthiazole | 1.64 | 0.17 | | |
| 2,4,5-trimethyl-3-thiazoline | 0.67 | | | |
| 2,4-dimethylthiazole | 2.03 | 2.48 | | |
| 2-acetylthiazole | | 0.22 | | |
| 2-thiazolyl ethyl ketone | 1.07 | | | |
| 4,5-dimethylisothiazole | | 0.24 | | |

### Thiophenes

| compound | F/C | F/G | F/A/S | F/S |
|---|---:|---:|---:|---:|
| 2,3,4-trimethylthiophene | | | | 0.47 |
| 2,4-dimethylthiophene | | | 0.11 | 0.17 |
| **2,5-dimethylthiophene** | **0.31** | **0.95** | **0.70** | **0.56** |
| 2,5-dimethyl-2-hydroxy-3(2H)-thiophenone | | | 3.25 | 6.28 |
| **2,5-dimethyl-4-hydroxy-3(2H)-thiophenone** (model `DMHFS`) | **0.22** | **0.15** | **1.26** | **1.16** |
| 2-acetylthiophene | | | 2.14 | 3.13 |
| 2-acetyl-3-methylthiophene | | | 0.62 | 0.64 |
| 2-acetyl-5-methylthiophene | | | 3.90 | |
| 2-ethyltetrahydrothiophen-3-one (left half of the table) | | | 2.46 | |
| 2-ethyl-5-methylthiophene | | | 0.16 | 1.92 |
| 2-ethyltetrahydrothiophen-3-one (right half of the table) | | | | 4.60 |
| 2-methylthiophene | 0.25 | | | |
| 2-methyltetrahydrothiophen-3-one | 0.82 | | | |
| 3-acetylthiophene | | | 0.34 | 0.48 |
| 3-acetyl-2,5-dimethylthiophene | | | 1.33 | 1.18 |
| 5-methyltetrahydrothiophen-3-one | 0.33 | | | |
| tetrahydrothiophen-3-one | 0.09 | | | |
| thiophene | 2.60 | | | |

**[!] Printing issues in Table 1:**
1. **2-ethyltetrahydrothiophen-3-one is listed twice**, once in each half (F/A/S 2.46 left; F/S 4.60
   right). This is either one compound split across two rows or a misnamed isomer. Not resolvable.
2. **"2-methyl-3-(thioacetyl)furan"** is ambiguous. It could be S-(2-methyl-3-furyl) thioacetate
   (the acetylated MFT, CH3C(O)S–) or a true thioacetyl (CH3C(S)–) furan. Either way it carries a
   2-methyl-3-furyl skeleton with sulfur at C-3, in the F/S pot only, and it is **not** free MFT.
3. The text (p. 895) calls the trithiolane "3,5-dimethyl-1,**2,3**-trithiolane"; the abstract and
   Table 1 say 1,**2,4**. The 1,2,4 isomer is the one that forms from acetaldehyde + H2S.
4. **2-mercapto-3-butanone and 3-mercapto-2-butanone are the same molecule**; the model has no species
   for it.

---

## 4. TEXT RESULTS [Q] (pp. 895–896)

- "carbonyls are the major products in the glutathione reaction, while sulfur-containing compounds
  predominate in the volatiles of the reactions of alanine/sodium sulfide and sodium sulfide". The
  authors attribute this to "the higher concentration of hydrogen sulfide as S²⁻" than cysteine or
  glutathione generate.
- "cysteine is more reactive than glutathione". Pyrazines and thiazoles/thiazolines are "more, in
  number and amount" in F/C than in F/G. This agrees with Zhang & Ho 1989/1991 (glucose or
  2,4-decadienal with Cys vs GSH, pH 7.5, 180 °C, 1 h).
- Mechanism (authors): "cysteine ... will participate in flavor reaction via the Maillard reaction
  and Strecker degradation more readily than via thermal degradation. Glutathione, however, will
  degrade to smaller intermediates such as hydrogen sulfide."
- Thiirane and 2-methylthiophene come from Strecker-derived mercaptoacetaldehyde: cyclisation gives
  thiirane, and condensation with acetone gives 2-methylthiophene. Hence they appear only in F/C.
- "2,5-dimethylthiophene is formed from the exchange of the oxygen atom in 2,5-dimethylfuran with the
  sulfur atom in hydrogen sulfide". [!] 2,5-Dimethylfuran is not in Table 1.
- "the majority of the thiols identified in all four reactions results from the reaction of hydrogen
  sulfide with aliphatic aldehydes or alcohols, the thermal degradation products of Furaneol"
  (citing Rizzi 1995, allylic alcohol substitution).
- "A higher reaction temperature will produce more intermediates and subsequently more volatile
  compounds" (citing Zheng & Ho 1994 for the temperature dependence of glutathione degradation).

---

## 5. DERIVED QUANTITIES [Z]

### 5.1 Class totals per pot (mg/mol, Table 1 cells summed; "<0.02" counted as 0.01)

| class | F/C | F/G | F/A/S | F/S |
|---|---:|---:|---:|---:|
| sulfur-containing (Sulfides + Thiazoles + Thiophenes) | **23.14** (15 cpds) | **6.09** (7) | **39.57** (22) | **64.17** (19) |
| of which thiirane | 11.48 | — | — | — |
| sulfur-containing, excluding the two cysteine-Strecker markers (thiirane, 2-methylthiophene) | 11.41 | 6.09 | 39.57 | 64.17 |
| aldehydes/ketones | 43.69 (16) | 29.76 (21) | 22.33 (12) | 19.62 (9) |
| of which 2,3-butanedione | 34.49 | 10.14 | 14.97 | 14.47 |

"Carbonyls are the major products in the glutathione reaction" holds as a share of that pot
(29.8 carbonyl vs 6.1 sulfur). In absolute terms F/C carries more carbonyl (43.7), mostly
2,3-butanedione.

### 5.2 Molar yields, IF mg/mol is mg per mol of Furaneol charged (§2.2; could be ×10 low)

| species | MW | F/C | F/G | F/A/S | F/S | unit |
|---|---:|---:|---:|---:|---:|---|
| DMHFS (the thiophenone) | 144.19 | 1.5 × 10⁻³ | 1.0 × 10⁻³ | 8.7 × 10⁻³ | 8.0 × 10⁻³ | mmol/mol DMHF |
| residual DMHF (recovered) | 128.13 | 0.034 | 0.060 | 0.0018 | 0.027 | mmol/mol DMHF |
| 2,3-butanedione | 86.09 | 0.40 | 0.12 | 0.17 | 0.17 | mmol/mol |
| 2,3-pentanedione | 100.12 | 1.9 × 10⁻³ | 2.7 × 10⁻³ | 1.0 × 10⁻³ | 2.0 × 10⁻³ | mmol/mol |
| 2-mercapto-3-butanone | 104.17 | 2.0 × 10⁻³ | — | 0.044 | 0.21 | mmol/mol |
| thiirane | 60.12 | 0.19 | — | — | — | mmol/mol |
| 3,5-dimethyl-1,2,4-trithiolane | 152.3 | 8.6 × 10⁻³ | 0.012 | 0.016 | 0.025 | mmol/mol |
| MP3P / MP2P | 118.20 | — / — | — / — | 3.3 / 4.1 × 10⁻³ | 3.9 / 2.6 × 10⁻³ | mmol/mol |

**[!] The "residual Furaneol" is 3.4 × 10⁻⁵ of the charge in F/C**, i.e. 0.003 %. Either >99.99 %
of the DMHF was consumed in 3 h at 130 °C, or, much more likely, SDE recovers very little of this
highly water-soluble furanone. **The residual cannot be read as conversion**, and the four residuals
(4.36 / 7.64 / 0.23 / 3.51) cannot rank DMHF consumption between pots: recovery depends on pH, and
the pots' pH differed (§2.1). Every yield in this table is therefore an extract amount, a **floor**
on formation and not a value.

### 5.3 Glutathione vs cysteine, compound by compound

| behaviour | compounds (F/C vs F/G, mg/mol) |
|---|---|
| **similar, within ~×3, both directions** (the products that need only H2S + Furaneol fragments) | DMHFS 0.22 vs 0.15; trithiolane 1.31 vs 1.88; 2,5-dimethylthiophene 0.31 vs 0.95; 2,4-dimethylthiazole 2.03 vs 2.48; 2,3-pentanedione 0.19 vs 0.27 |
| **cysteine only** (Strecker or amine-dependent) | thiirane 11.48; thiophene 2.60; 2-methylthiophene 0.25; three tetrahydrothiophen-3-ones (0.82, 0.33, 0.09); 2,4,5-trimethyl-3-thiazoline 0.67; 2-thiazolyl ethyl ketone 1.07; 2,5-dimethyl-4-mercapto-3(2H)-furanone 0.11; 2-mercapto-3-butanone 0.21; 2-ethyl-5-methylpyrazine 0.40 |
| **glutathione only** | 2-acetylthiazole 0.22; 4,5-dimethylisothiazole 0.24 |
| **lower with glutathione** | 2,4,5-trimethylthiazole 1.64 → 0.17 (×0.10); 2,5-dimethylpyrazine 0.37 → 0.14 |

**Verdict on "does glutathione behave like cysteine?"** Only as an **H2S donor**. Its products are
almost a **subset** of the H2S-only products and appear at similar levels. Cysteine adds a family that
needs its free α-amino group (thiirane, thiophene, 2-methylthiophene, the tetrahydrothiophenones,
the thiazoline), and that family is half of the F/C sulfur total. By Zheng & Ho 1994 the two donors
release H2S at nearly the same first-order rate at 100–110 °C, pH 7 (ratio GSH/Cys 0.98–1.16;
`zheng1994_extraction.md` §4). So the **sulfide supply should be similar; the difference is the
amine chemistry, plus the pH confound of §2.1.**

---

## 6. WHAT THIS CAN AND CANNOT ANCHOR IN THE MODEL

**The model's relevant pieces:** `DMHF` (trunk, B7), `DMHFS` made by `r_dmhf_h2s` (DMHF + H2S →
DMHFS, "Edge C", **rate exactly zero by declaration**: `parameters_furanic.EDGE_C_ZERO_BY_DECLARATION`,
which also calls Shu & Ho 1988 "the only fed-precursor DMHF + cysteine experiment in the corpus");
`GSH` as "a sulfide source only" (Zheng & Ho 1994); the fed-norfuraneol branch `NF` + H2S → `MFT` /
`MP3P` / `SHET` (Whitfield & Mottram 1999/2001); `MP3P` labelled "NORFURANEOL-diagnostic", `MP2P`
"THIAMINE-diagnostic".

1. **Edge C (`r_dmhf_h2s`): a second fed DMHF + sulfur experiment exists, and it supplies a product
   amount, but still no rate.** This paper reports DMHFS in mg per mol of DMHF charged, and a
   (recovered) residual DMHF, which Shu & Ho do not. That is more than Shu & Ho give. It is still a
   **single time point, single temperature, unknown recovery, RRF assumed 1, unit ambiguous ×10, pH
   unknown**. The only defensible use is a **floor**: DMHFS/DMHF₀ ≥ ~1 × 10⁻⁶ mol/mol (F/C) and
   ≥ ~8 × 10⁻⁶ (F/S) after 3 h at 130 °C. A floor on second-order k from F/S [Z]:
   8 × 10⁻⁶ / (0.1 M × 180 min) ≈ 4 × 10⁻⁷ M⁻¹ min⁻¹, i.e. **≥ ~4 × 10⁻¹⁰ L mmol⁻¹ min⁻¹**. The
   assumptions behind that are no secondary consumption, recovery = 1 and [H2S] constant. **It is
   orders of magnitude too weak to move anything**, and it says only that the true rate is not
   exactly zero. **Recommendation: keep Edge C at zero.** Cite this paper as the second fed
   experiment, so the "only fed-precursor DMHF + cysteine experiment" wording in
   `EDGE_C_ZERO_BY_DECLARATION`, `engine.py` and `species_sulfur.py` is out of date. Do not fit it.
2. **A sign / presence HOLD-OUT that the shipped zero already fails:** DMHFS is **present** in all four
   pots at nominal pH 7.0 (buffered), including F/C (0.22). Shu & Ho 1988 found it not detected at
   pH 7.1 (unbuffered, area %). The two are not in conflict, because detection limits and workups
   differ. But "DMHFS > 0 when DMHF meets a sulfide at ~pH 7, 130 °C" is a qualitative fact a zero
   rate cannot reproduce. Since the zero is declared, this is a documented, pre-registered miss.
3. **MFT is absent when DMHF is the only furanone fed (all four pots).** The model has no DMHF → MFT
   route, so this agrees. **It does NOT test the norfuraneol branch**: Furaneol is the C6 homologue,
   not norfuraneol (C5), and the paper never feeds NF. Do not use it as an NF hold-out.
4. **Both mercaptopentanone isomers form from DMHF + Na2S with no pentose, thiamine or norfuraneol**
   (MP3P 0.39 / 0.46, MP2P 0.49 / 0.31; F/A/S, F/S). 2,3-Pentanedione is present in every pot
   (0.10–0.27), which is the Mottram 1995 dione + H2S route. **The model's "NORFURANEOL-diagnostic"
   (`MP3P`) and "THIAMINE-diagnostic" (`MP2P`) labels therefore hold only where no DMHF-derived C5
   dione is present.** In a hexose pot that makes DMHF, both isomers have a third, sugar-side
   source. Qualitative. Not convertible to a rate (no time course, extract amounts only), and it was
   seen only in the high-pH sulfide pots, not with cysteine.
5. **Glutathione as "a sulfide source only" (`GSH`): supported in kind.** F/G gives no product that
   needs a free cysteine amine, and the H2S-type products it does give match F/C within ~×3 (§5.3).
   It also refutes treating **cysteine** as only a sulfide source. Half of the F/C sulfur-volatile
   total (thiirane 11.48 of 23.14) comes from Strecker chemistry of cysteine's own carbons
   (mercaptoacetaldehyde). The model does not carry thiirane, thiophene or 2-methylthiophene, so
   there is nothing to test now. The ~3.8× GSH < Cys sulfur-volatile ratio is **not usable as a
   magnitude**: it is confounded by the probable ~3-unit pH difference (§2.1).
6. **The product split resembles the NF + H2S result (`SHET`).** Fed furanone + sulfide gives mostly
   thiophenes, thiophenones, trithiolanes and mercaptoketones, and no furanthiol. In F/S the
   thiophene/thiophenone class is 20.6 mg/mol [Z] of 64.2. That is the same qualitative message as
   the Whitfield & Mottram rows behind `SHET`, for the C6 homologue. The model has no DMHF-side
   `SHET` pool, so it is nothing to test now.

**Cannot anchor:** any rate constant, Ea, pH dependence or time dependence (one time, one
temperature, one nominal pH); any absolute concentration (RRF = 1, unknown SDE recovery, unit
ambiguity); DMHF conversion; H2S (not measured, lost in SDE); anything about MFT or FFT except that
DMHF does not make them.

---

## 7. CAVEATS

1. **pH**: set before the 0.1 M reagents went into a 0.1 M buffer; never measured afterwards. My §2.1
   estimates put the four pots at roughly pH 7 (or 2), 3.5–4, ~10 and ~12.5–13. Every cross-pot
   comparison is donor × pH.
2. **Quantification** against Furaneol with RRF = 1 for compounds of very different FID response
   (C2 thiirane to C8 acetylthiophenes); SDE recovery uncorrected and certainly compound-specific.
3. **Unit**: "mg/mol" with an equation that divides by molarity; possibly ×10 off "per mol charged".
4. **No replicates, no blanks** (Furaneol heated alone, or cysteine alone, was not run), so the
   Furaneol-only background of carbonyls is unknown.
5. **Parr bomb with a 400 mL headspace**: volatile H2S and low boilers partition out of the liquid;
   not quantified.
6. **Identification** by library MS + RI; no per-compound confirmation flags. Two naming issues
   (§3: duplicated row, ambiguous "thioacetyl") and one text/table isomer inconsistency.
7. **Na2S as an "H2S" surrogate** delivers sulfide as HS⁻/S²⁻ at high pH. That is not the same
   reagent as H2S released gradually from cysteine at pH 7.
