# Wang, Lin, Jia, Sun, Wang & Chen 2026 — EXTRACTION (xylose–cysteine at 120 °C over 150 minutes: the three canonical thiol dimers, as odour-activity values only)

**Source on disk:** `data/articles/Wang2026b.pdf` (10.1 MB; downloaded 2026-09-11). Read 2026-09-11
via `pdftotext -layout` with Table 1 verified against a page render. Wave B37.

| field | value |
|---|---|
| Title | "Characterization and intrinsic correlations of xylose-cysteine Maillard reaction products across different reaction processes: insights into reaction degree, functional group variations, antioxidant potential and aroma formation" |
| Venue | Food Chemistry 521 (2026) 149913 |
| DOI | 10.1016/j.foodchem.2026.149913 |
| System | **D-xylose 0.300 g + L-cysteine 0.242 g, 1:1 molar, in 10 mL of 0.2 mol/L NaH₂PO₄ buffer at pH 5.5**, 15 mL pressure-resistant vessel, oil bath, **120 °C**, 60/90/120/150 min, plus an unreacted 0-min control |
| Quantification | HS-SPME (50/30 µm DVB/CAR/PDMS), 3.0 mL, 60 °C; cyclohexanone internal standard, **"under the assumption of equal response factors among analytes"** |

## 1. What this paper has that no other in the corpus has

It is the only source on disk that reports **all three canonical thiol dimers together** — bis(2-methyl-3-furyl)
disulfide, difurfuryl disulfide and the mixed 2-furfuryl 2-methyl-3-furyl disulfide — alongside their
two parent thiols, in one xylose–cysteine pot, against time. The repository's sulfur lane makes
exactly these species and has had no external time course for them.

## 2. And the reason it cannot be scored: the concentrations are not printed

The 22-compound quantification is in **Table S1, which is supplementary and not on disk**.
Figure 4A ("Concentration of sulfur-containing volatile compounds with OAV > 100, ng/mL" against
reaction time) and Figure 4E (class totals) are **figure-only, with a broken axis and no printed
values**. What Table 1 prints is **odour-activity values**, which are concentration divided by an
odour threshold taken from a compilation. Back-multiplying OAV by threshold would reconstruct a
concentration, and reconstruction is not transcription: it is not done here and the repository must
not do it.

## 3. Table 1 — odour-activity values (no standard deviations printed)

| compound | threshold (mg/kg) | 0 min | 60 min | 90 min | 120 min | 150 min |
|---|---:|---:|---:|---:|---:|---:|
| **2-Methyl-3-furanthiol** | 0.00001 | – | 54 373 000 | 183 061 000 | 187 738 000 | 240 732 000 |
| **2-Furfurylthiol** | 0.00013 | – | 2 112 153 | 10 515 000 | 12 503 307 | 16 309 538 |
| **Bis(2-methyl-3-furyl) disulfide** (printed as 3,3′-dithiobis[2-methyl-furan]) | 0.00000002 | – | 1.865e11 | 2.928e11 | 3.636e11 | 4.361e11 |
| **Difurfuryl disulfide** | 0.00015 | – | 506 067 | 4 585 333 | 5 571 533 | 6 821 067 |
| **2-Furfuryl 2-methyl-3-furyl disulfide** | 0.5 | – | 145.98 | 1375.60 | 1671.46 | 2046.32 |
| 3-Mercapto-2-pentanone | 0.0007 | – | 34 985.71 | 274 871.43 | 373 300.00 | 445 171.43 |
| 5-Methyl-2-thiophenecarboxaldehyde | 0.001 | – | – | – | 16 610 | 19 780 |
| Furfural | 9.562 | – | 26.69 | 40.94 | 43.20 | 44.27 |
| 2-Methylfuran | 4 | – | 4.92 | 7.87 | 16.98 | 36.26 |
| 2-Thiophenecarboxaldehyde | 5 | – | – | – | 4.87 | 3.67 |
| 2-Methylthiophene | 3 | – | – | – | – | 8.19 |
| Decane | 3.6 | 17.01 | – | – | – | – |

## 4. What survives as usable evidence: the ORDER of appearance

Verbatim: *"At 60 min, five key meaty aroma compounds (2-furfurylthiol, 2-methyl-3-furanthiol,
2-furfuryl 2-methyl-3-furyl disulfide, 3,3′-dithiobis[2-methyl-furan], and 3-mercapto-2-pentanone)
rapidly emerged. As the reaction progressed to 90 min, a new disulfide (difurfuryldisulfide) was
generated. Subsequently, at 120 min, 5-methyl-2-thiophenecarboxaldehyde appeared."*

And the mechanism the repository's dimer step assumes, stated by the authors: *"Disulfides were
typically formed through the oxidative dimerization of primary monomeric thiols. Therefore,
difurfuryldisulfide could only be synthesized after sufficient accumulation of its precursor,
2-furfurylthiol, in the reaction system."*

Compound counts across the series: 2, 8, 15, 20 and 21 volatiles at 0/60/90/120/150 min. Class
shares: thioethers fell from 71 % to 54 % of total volatiles while thiophenes rose from 16 % to 34 %.

**Within-OAV ratios are transferable where the two compounds share a pot and a time**, because the
thresholds are constants: for example the mixed disulfide's OAV rises 14.0-fold from 60 to 150 min
while its parent 2-furfurylthiol rises 7.7-fold, which says the dimer outruns the monomer. That is
a within-study ratio and therefore fit-class evidence under the owner's rule — but it is a ratio of
two different compounds, so it also carries the equal-response-factor assumption, and it is recorded
here rather than used.

## 5. Verdict

The right experiment, the right compounds, and the numbers are in a supplementary file that did not
come with the PDF. The one action it supports is a request: **Table S1 of this paper would give the
sulfur lane its first external thiol-and-dimer time course.** Added to the reading list.

## Supplementary material (read 2026-10-29)

**Source:** `data/articles/wang2026b_supplementary.docx` (owner download 2026-10-29, not committed). Read
with macOS `textutil` and, because empty cells drop out of the text conversion, cell by cell from the
docx XML (scratchpad script `docx_table_dump.py`). No title page; identity from content (22 compounds,
MRPs-0 … MRPs-150 of xylose and cysteine, matching sec. 0 and Table 1).

**What it holds.**
- **Table S1** — "Identification of volatile compounds and their relative abundances, semi-quantitative
  concentrations in the Maillard reaction products (MRPs) determined by gas chromatography-mass
  spectrometry." 22 compounds × 5 times; relative abundance (%) and **semi-quantitative concentration
  (ng/mL)**, mean ± SD, n = 3, Duncan letters.
- **Table S2** — Pearson correlation matrix of 16 indicators (pH, A294, A420, fluorescence, DPPH, ABTS,
  reducing power, three e-nose sensors W5S/W1W/W2W, and six volatiles). Correlations only; **no model
  input; not transcribed**.
- No figures.

**Conditions** (main text, sec. 0): D-xylose 0.300 g + L-cysteine 0.242 g (2 + 2 mmol, 1:1) in 10 mL
0.2 mol/L NaH2PO4, pH 5.5, 15 mL pressure vessel, oil bath, 120 °C; sampled at 0, 60, 90, 120, 150 min.
The SI footnote restates only the times ("MRPs of xylose and cysteine obtained at 0, 60, 90, 120, and
150 min").

**Quantification basis: semi-quantitative, response factor 1** — cyclohexanone IS "under the assumption
of equal response factors among analytes" (main text), and the SI column header itself says
"Semi-quantitative concentration (ng/mL)". HS-SPME of 3.0 mL at 60 °C.

### Table S1 (as printed). Relative abundance (%) | semi-quantitative concentration (ng/mL, mean ± SD, n = 3; letters = Duncan p < 0.05). "-" / ND = not detected; one empty cell left empty.

| No. | compound | RI | CAS | formula | % 0 | % 60 | % 90 | % 120 | % 150 | ng/mL 0 | ng/mL 60 | ng/mL 90 | ng/mL 120 | ng/mL 150 |
|---:|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| | *Furan compounds* | | | | | | | | | | | | | |
| 1 | 2-Methylfuran | 848 | 534-22-5 | C5H6O | - | 0.39 | 0.85 | 0.64 | 0.84 | ND | 19.68 ± 0.35c | 61.47 ± 15.81b | 67.93 ± 11.81b | 145.02 ± 18.38a |
| 2 | Furfural | 1459 | 98-01-1 | C5H4O2 | - | 5.07 | 2.41 | 4.76 | 2.23 | ND | 255.22 ± 5.35c | 391.50 ± 27.60b. | 413.05 ± 6.92ab | 423.35 ± 19.61a |
| 3 | 2,5-Dimethylfuran | 940 | 625-86-5 | C6H8O | - | - | 0.08 | 1.08 | 0.38 | ND | ND | 24.65 ± 13.03b | 31.22 ± 1.33b | 75.64 ± 24.95a |
| 4 | **2-Methyl-3-furanthiol** | 1305 | 28588-74-1 | C5H6OS | - | 11.75 | 10.64 | 6.88 | 16.43 | ND | **543.73 ± 124.61c** | **1830.61 ± 191.56b** | **1877.38 ± 47.22b** | **2407.32 ± 420.60a** |
| 5 | **2-Furfurylthiol** | 1430 | 98-02-2 | C5H6OS | - | 5.71 | 9.85 | 17.70 | 10.91 | ND | **274.58 ± 7.72d** | **1366.95 ± 89.99c** | **1625.43 ± 209.03b** | **2120.24 ± 93.03a** |
| | *Alkane compounds* | | | | - | | | | | | | | | |
| 6 | Decane | 996 | 124-18-5 | C10H22 | 78.24 | - | - | - | - | 61.25 ± 1.66 | ND | ND | ND | ND |
| | *Thiophene compounds* | | | | | | | | | | | | | |
| 7 | 2-Methylthiophene | 1080 | 554-14-3 | C5H6S | - | - | - | - | 0.12 | ND | ND | ND | ND | 24.58 ± 5.38 |
| 8 | 2-Thiophenecarboxaldehyde | 1682 | 98-03-3 | C5H4OS | - | - | - | 0.34 | 0.09 | ND | ND | ND | 24.37 ± 3.76a | 18.33 ± 0.31b |
| 9 | 2-Acetyl-3-methylthiophene | 1750 | 13679-72-6 | C7H8OS | - | - | - | 1.01 | 0.28 | ND | ND | ND | 23.06 ± 8.88b | 44.90 ± 1.93a |
| 10 | 5-Methyl-2-thiophenecarboxaldehyde | 1797 | 13679-70-4 | C6H6OS | - | - | - | 1.14 | 0.13 | ND | ND | ND | 16.61 ± 1.70b | 19.78 ± 2.22a |
| 11 | Thieno[3,2-b] thiophene | 1863 | 251-41-2 | C6H4S2 | - | - | 0.04 | 0.60 | 0.72 | ND | ND | 19.41 ± 4.28c | 48.05 ± 7.67b | 131.59 ± 16.80a |
| 12 | 2,5-Thiophenedicarboxaldehyde | 1899 | 932-95-6 | C6H4O2S | - | - | 0.85 | 2.82 | 2.77 | ND | ND | 169.28 ± 8.14c | 261.46 ± 20.75b | 458.52 ± 34.61a |
| 13 | 2-Methylthieno[3,2-b] thiophene | 1929 | 13393-75-4 | C7H6S2 | - | - | 0.06 | 0.96 | 1.43 | ND | ND | 32.61 ± 3.66c | 69.95 ± 15.23b | 242.92 ± 36.87a |
| | *Ketone compounds* | | | | | | | | | | | | | |
| 14 | **3-Mercapto-2-pentanone** | 1354 | 67633-97-0 | C5H10S | - | 0.48 | 1.06 | 2.99 | 1.56 | ND | 24.49 ± 2.33d | 192.41 ± 13.85c | 261.31 ± 14.91b | 311.62 ± 5.29a |
| | *Phenolic compounds* | | | | | | | | | | | | | |
| 15 | 4-(Methylthio)thiophenol | 1951 | 1122-97-0 | C7H8S2 | - | - | 0.73 | 1.84 | 2.95 | ND | ND | 111.96 ± 20.47c | 142.48 ± 13.26b | 444.21 ± 10.87a |
| 16 | 2,4-Di-tert-butylphenol | 2309 | 96-76-4 | C14H22O | 21.76 | 0.34 | 0.65 | 1.68 | 1.35 | 16.73 ± 0.21c | 15.15 ± 1.36c | 150.87 ± 32.43b | 170.50 ± 32.43b | 237.34 ± 33.43a |
| | *Thioether compounds* | | | | | | | | | | | | | |
| 17 | **3,3´-Dithiobis[2-methyl-furan]** (bis(2-methyl-3-furyl) disulfide) | 2145 | 28588-75-2 | C10H10O2S2 | - | 74.67 | 66.36 | 43.83 | 45.40 | ND | **3730.50 ± 142.70d** | **5856.47 ± 117.54c** | **7271.36 ± 585.83b** | **8722.94 ± 244.80a** |
| 18 | **2-Furfuryl 2-methyl-3-furyl disulfide** | 2389 | 109537-55-5 | C10H10O2S2 | - | 1.59 | 5.10 | 5.98 | 5.97 | ND | **75.91 ± 3.98d** | **687.80 ± 49.75c** | **835.73 ± 67.61b** | **1023.16 ± 140.27a** |
| 19 | **Difurfuryldisulfide** | 2558 | 4437-20-1 | C10H10O2S2 | - | *(empty)* | 0.76 | 2.98 | 2.93 | ND | ND | **142.97 ± 10.24c** | **266.02 ± 59.27b** | **480.59 ± 85.01a** |
| | *Ester Compounds* | | | | | | | | | | | | | |
| 20 | 1,2-Benzenedicarboxylic acid, bis(2-methylpropyl) ester | 2531 | 84-69-5 | C16H22O4 | - | - | - | 0.72 | 0.53 | ND | ND | ND | 72.61 ± 11.75b | 112.79 ± 11.71a |
| | *Others* | | | | | | | | | | | | | |
| 21 | Hexahydro-1,2,4,5-tetrazine-3,6-dithione | 1571 | 36239-33-5 | C2H4N4S2 | - | - | 0.55 | 0.97 | 1.40 | ND | ND | 100.71 ± 12.07b | 82.09 ± 5.94b | 174.64 ± 39.80a |
| 22 | 2-Pyrimidinamine | 2288 | 109-12-6 | C4H5N3 | - | - | - | 1.09 | 1.57 | ND | ND | ND | 108.01 ± 18.50b | 320.97 ± 30.05a |

Footnote as printed: "RI, retention index; CAS, chemical abstracts service; ND, not detected. All data
were expressed as mean ± standard deviation (SD) based on three independent measurements. Different
letters in the same line indicated significant differences (p < 0.05)." (The stray "." after furfural's
"b" at 60 min is printed.)

**Checks and findings (mine).**
1. **★ Table 1's OAVs are inflated 10^6 by a unit error.** Every Table 1 OAV equals the S1 concentration in
   **ng/mL** divided by the threshold in **mg/kg**, with no conversion: MFT 60 min 543.73 / 0.00001 =
   54 373 000 (printed 54 373 000); FFT 274.58 / 0.00013 = 2 112 154 (2 112 153); bis(2-methyl-3-furyl)
   disulfide 3730.50 / 2e-8 = 1.865e11 (1.865e11); furfural 255.22 / 9.562 = 26.69 (26.69); decane
   61.25 / 3.6 = 17.01. Since 1 ng/mL ≈ 10^-6 mg/kg, **the true OAVs are 10^6 smaller**: MFT 54-241,
   FFT 2.1-16, bis(2-methyl-3-furyl) disulfide 1.9e5-4.4e5, difurfuryl disulfide 0.95-3.2 (above 1 only at
   120-150 min), 2-furfuryl 2-methyl-3-furyl disulfide ~0.0001-0.002, 3-mercapto-2-pentanone 0.035-0.45,
   5-methyl-2-thiophenecarboxaldehyde 0.017-0.020, furfural/2-methylfuran/decane ≪ 1. The paper's
   "OAV > 100" selection (Fig. 4A) and its odour-activity ranking are wrong in magnitude; within-compound
   time ratios (sec. 4) are unaffected.
2. **Sec. 2 is now closed: the concentrations are printed** (no back-reconstruction from OAV needed). The
   sec. 4 ratio "mixed disulfide rises 14.0-fold 60 → 150 min while FFT rises 7.7-fold" reproduces from
   S1 (1023.16/75.91 = 13.5; 2120.24/274.58 = 7.72; the 14.0 came from rounded OAVs).
3. **Compound counts per time** (2, 8, 15, 20, 21) and the order of appearance quoted in sec. 4 reproduce
   from S1. The relative-abundance columns sum to 100 % (e.g. 90 min 99.99 %).
4. **Dimer vs parent, semi-quantitative (RF = 1), ng/mL:** bis(2-methyl-3-furyl) disulfide is **6.9x MFT at
   60 min and 3.6x at 150 min** by mass; difurfuryl disulfide is 0.10-0.23x FFT (90-150 min). In moles
   (MW MFT 114.17, dimer 226.31) that is ~3.5 dimer molecules per free MFT at 60 min, i.e. ~6.9x the free
   MFT in MFT-equivalents. With RF = 1 these
   cross-compound ratios carry the response-factor error, so they are order-of-magnitude only.
5. **Contaminant-looking peaks:** decane (only at 0 min), 2,4-di-tert-butylphenol (present at 0 min) and
   diisobutyl phthalate are common SPME/plastic artefacts; they enter the relative-abundance denominators.
6. Non-monotone: 2-thiophenecarboxaldehyde 24.37a → 18.33b (120 → 150 min); hexahydro-1,2,4,5-tetrazine-
   3,6-dithione 100.71 → 82.09 (90 → 120 min, same letter).
7. **No pH, cysteine, xylose or dicarbonyl values** are in the SI (Table S2 correlates pH and A294/A420 with
   volatiles but prints no pH values).

**Model-usable:** yes (semi-quantitative). The corpus's first printed thiol-and-dimer time course in one
xylose-cysteine pot: MFT, FFT, bis(2-methyl-3-furyl) disulfide, difurfuryl disulfide, the mixed disulfide
and 3-mercapto-2-pentanone at 60/90/120/150 min, 120 °C, pH 5.5, with SDs. Within-compound time shapes
are usable; absolute ng/mL and cross-compound ratios carry RF = 1.
