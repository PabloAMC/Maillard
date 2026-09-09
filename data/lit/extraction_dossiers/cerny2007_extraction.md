# Cerny & Briffod 2007 — EXTRACTION (pH 4.0-7.0 ladder for sulfur volatiles from [13C5]xylose + cysteine + thiamin, 145 C / 20 min)
### The only five-point pH series on MFT / FFT / furfural / mercaptopentanones in the corpus — but peak areas, and a ternary system with thiamin.

**Source on disk:** `data/articles/cerny2007.pdf` (owner's download, 2026-09-07). Read-only extraction from
`pdftotext -layout`; Tables 1-5 all have clean text layers and were re-typed in full. Repo status before
this dossier (FIT_HOLDOUT_DECLARATION / k3 §B10): **Table 2 pH ladder = HOLD-OUT (directional only)**;
**Table 4 isotope splits across pH = FIT**; **Table 5 concentration pair (85:15 vs 54:46) = HOLD-OUT**.
Not to be confused with `cerny2007b.pdf` (Cerny, LWT 2007, the pH-5 origin-of-carbons paper this one
cites as ref 9).

## 0. Identity

| field | value |
|---|---|
| Title | "Effect of pH on the Maillard Reaction of [13C5]Xylose, Cysteine, and Thiamin" |
| Authors | Christoph Cerny and Matthieu Briffod (Firmenich SA, Geneva) |
| Venue | J. Agric. Food Chem. 2007, 55, 1552-1556 |
| DOI | 10.1021/jf062874w |
| Companion | ref 9 = Cerny, LWT 2007 (in press at the time): same system at pH 5.0, 2x the concentration |

## 1. Why it matters

The model's thiol-forming steps (NF + H2S -> MFT; furfural + H2S -> FFT) are fitted at pH 4.5-5 only and
carry no pH dependence. This paper heats the same pentose-cysteine chemistry (xylose instead of ribose)
at pH 4.0 / 5.0 / 5.5 / 6.0 / 7.0, 145 C / 20 min — the Hofmann 1998 temperature and time — and reads
MFT, FFT, furfural, both mercaptopentanones and eight other sulfur volatiles at every pH. The 13C5-xylose
label separates, at each pH, the fraction of MFT and 3-mercapto-2-pentanone that came from the pentose
from the fraction that came from thiamin, so a pentose-only pH shape can be reconstructed. The price:
GC-TIC peak areas (no response factors, no calibration), SPME (pH-dependent partition), and 50 mM
thiamin in the pot.

## 2. Methods as they matter to a model

- **Charges (Table 1, all in 1.00 mL of 0.5 mol/L potassium phosphate buffer):** cysteine 6.25 mg
  (= 51.6 µmol; text says 50 µmol) -> **~50 mmol/L**; xylose 22.50 mg (149.9 µmol unlabelled; 145 µmol
  as [13C5]xylose, MW 155.1) -> **150 mmol/L**; thiamin hydrochloride 16.85 mg (49.96 µmol) ->
  **50 mmol/L**. Molar ratio **xylose : cysteine : thiamin = 3 : 1 : 1**. Runs A-D use [13C5]xylose
  (99 % enrichment) at pH 4.00 / 5.00 / 6.00 / 7.00; runs E-I use unlabelled xylose at pH 4.00 / 5.00 /
  5.50 / 6.00 / 7.00.
- **Buffer / pH:** "Different buffers with pH values from 4.0 to 7.0 were used" — potassium phosphate
  0.5 mol/L throughout (K2HPO4 / KH2PO4). Post-reaction pH not reported.
- **Vessel / heating:** Teflon vials (Infochroma), stirred, heated metal block (Reacti-Therm), **145 C,
  20 min**. Headspace and atmosphere not stated (closed vials).
- **Replicates:** unlabelled runs E-I in **triplicate**; Table 2 = mean of triplicates, "The standard
  deviation did not exceed 23%." Labelled runs A-D: replication not stated (isotope ratios only).
- **Analysis:** reacted solution transferred to 20 mL headspace vials; HS-SPME (DVB/CAR/PDMS fibre per
  the abbreviations list; method of ref 10) 35 min at 40 C; GC-MS, HP-5MS, 40 C (5 min) -> 5 C/min ->
  260 C -> 15 C/min -> 280 C.
- **Quantification: RELATIVE ONLY.** "Peak areas were integrated using the TIC signals from the
  experiments with unlabeled xylose and expressed as mean values from triplicates." "The data in Table 2
  represent integrated peak areas and are not corrected by MS response factors or taking into account
  the different partition coefficients of volatiles between the sample matrix and the SPME fiber. Also,
  for certain compounds, the pH might have an effect on the partition coefficient." No internal
  standard is mentioned. Nothing in the paper can be converted to mol %.
- **Isotope ratios:** from relative intensities of the molecular-ion cluster (Table 3) and from the two
  analysed ions M and M+5 (Table 4); pH 5.5 was not run with label.
- Identification: authentics for all but 6, 8, 11, 12 (tentative, literature RI/MS).

## 3. Tables re-typed

### Table 1. "Model Reactions" (footnote a: "Reaction in phosphate buffer (0.5 mol/L; 1.00 mL) at 145 °C (20 min)"; b: "Amount (mg)")

| run | A | B | C | D | E | F | G | H | I |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| pH | 4.00 | 5.00 | 6.00 | 7.00 | 4.00 | 5.00 | 5.50 | 6.00 | 7.00 |
| cysteine (mg) | 6.25 | 6.25 | 6.25 | 6.25 | 6.25 | 6.25 | 6.25 | 6.25 | 6.25 |
| [13C5]xylose (mg) | 22.50 | 22.50 | 22.50 | 22.50 | | | | | |
| xylose (mg) | | | | | 22.50 | 22.50 | 22.50 | 22.50 | 22.50 |
| thiamin hydrochloride (mg) | 16.85 | 16.85 | 16.85 | 16.85 | 16.85 | 16.85 | 16.85 | 16.85 | 16.85 |

### Table 2. "Influence of the pH on the Formation of Sulfur Volatiles from Xylose, Cysteine, and Thiamin (GC-TIC Peak Areas × 10^6)"

Footnote a: identified vs authentic reference compounds unless marked c; b: RI on HP-5MS; c: tentative
(literature MS and RI, refs 11-14). Means of triplicates, SD <= 23 %.

| no. | compound | RI | pH 4.0 | pH 5.0 | pH 5.5 | pH 6.0 | pH 7.0 |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | 4,5-dihydro-2-methyl-3(2H)-furanone | 800 | 17 | 18 | 18 | 94 | 104 |
| 2 | 3-mercapto-2-butanone | 807 | 69 | 81 | 83 | 59 | 73 |
| 3 | **2-furaldehyde** | 831 | **208** | **158** | **165** | **0** | **0** |
| 4 | **2-methyl-3-furanthiol** | 871 | **415** | **336** | **371** | **539** | **391** |
| 5 | **3-mercapto-2-pentanone** | 901 | **145** | **141** | **147** | **62** | **12** |
| 6 | **2-mercapto-3-pentanone** (c) | 911 | **0** | **0** | **0** | **19** | **5** |
| 7 | **2-furfurylthiol** | 915 | **431** | **368** | **364** | **185** | **0** |
| 8 | 4,5-dihydro-2-methyl-3-furanthiol (c) | 939 | 226 | 43 | 41 | 18 | 7 |
| 9 | 2-methyl-3-(methylthio)furan | 950 | 14 | 12 | 11 | 11 | 0 |
| 10 | 4,5-dihydro-2-methyl-3(2H)-thiophenone | 988 | 0 | 0 | 0 | 372 | 470 |
| 11 | 3-methyl-1,2-dithian-4-one (c) | 1162 | 0 | 0 | 0 | 75 | 82 |
| 12 | 3-acetyl-1,2-dithiolane (c) | 1201 | 0 | 0 | 0 | 117 | 240 |
| 13 | 5-(2-hydroxyethyl)-4-methylthiazole | 1261 | 525 | 54 | 43 | 161 | 80 |

### Table 3. "Isotope Ratios (Percent) of Sulfur Volatiles Formed from [13C5]Xylose, Cysteine, and Thiamin, Reacted at Different pH Values"

Footnote a: "[12C5]" column = the unlabelled compound from the unlabelled-xylose reaction at pH 4;
b: na = not analyzed.

2-methyl-3-furanthiol (4):

| m/z | [12C5] | pH 4 | pH 5 | pH 6 | pH 7 |
|---:|---:|---:|---:|---:|---:|
| 112 | 0.6 | 0.5 | 0.6 | 0.5 | <0.1 |
| 113 | 17.8 | 13.3 | 15.1 | 16.1 | 14.2 |
| 114 | 72.4 | 53.6 | 61.0 | 65.3 | 57.7 |
| 115 | 5.4 | 3.9 | 4.5 | 4.7 | 4.2 |
| 116 | 3.6 | 2.7 | 3.1 | 3.2 | 2.9 |
| 117 | 0.2 | 0.7 | 0.5 | 0.4 | 0.5 |
| 118 | <0.1 | 5.8 | 3.5 | 2.3 | 4.7 |
| 119 | <0.1 | 18.2 | 10.8 | 7.0 | 14.7 |
| 120 | <0.1 | 0.5 | 0.3 | 0.2 | 0.4 |
| 121 | <0.1 | 0.8 | 0.6 | 0.3 | 0.7 |

3-mercapto-2-pentanone (5):

| m/z | [12C5] | pH 4 | pH 5 | pH 6 | pH 7 |
|---:|---:|---:|---:|---:|---:|
| 116 | 0.5 | 1.0 | 1.3 | 3.1 | 2.1 |
| 117 | <0.1 | <0.1 | <0.1 | <0.1 | <0.1 |
| 118 | 89.1 | 67.1 | 75.9 | 75.5 | 76.0 |
| 119 | 5.6 | 4.6 | 4.9 | 6.5 | 4.6 |
| 120 | 4.3 | 3.4 | 3.9 | 4.0 | 4.3 |
| 121 | 0.3 | 0.6 | 0.4 | 1.0 | 0.5 |
| 122 | <0.1 | 1.4 | 0.8 | 0.7 | 1.3 |
| 123 | <0.1 | 20.4 | 12.1 | 8.2 | 10.2 |
| 125 | <0.1 | 0.3 | 0.2 | 0.4 | 0.2 |
| 126 | <0.1 | 1.2 | 0.5 | 0.6 | 0.8 |

4,5-dihydro-2-methyl-3-furanthiol (8):

| m/z | [12C5] | pH 4 | pH 5 | pH 6 | pH 7 |
|---:|---:|---:|---:|---:|---:|
| 114 | 0.7 | 1.4 | 1.3 | 1.3 | na |
| 115 | 6.1 | 6.1 | 6.0 | 6.3 | na |
| 116 | 83.3 | 82.4 | 82.4 | 83.0 | na |
| 117 | 5.6 | 5.8 | 5.9 | 5.3 | na |
| 118 | 4.1 | 3.8 | 3.9 | 3.7 | na |
| 119 | 0.3 | 0.5 | 0.5 | 0.4 | na |
| 120-123 | <0.1 | <0.1 | <0.1 | <0.1 | na |

5-(2-hydroxyethyl)-4-methylthiazole (13):

| m/z | [12C5] | pH 4 | pH 5 | pH 6 | pH 7 |
|---:|---:|---:|---:|---:|---:|
| 141 | 0.2 | 0.4 | 0.3 | 0.2 | <0.1 |
| 142 | <0.1 | <0.1 | <0.1 | <0.1 | <0.1 |
| 143 | 91.3 | 88.3 | 88.5 | 87.6 | 88.4 |
| 144 | 5.0 | 6.6 | 6.7 | 6.9 | 6.7 |
| 145 | 3.3 | 4.2 | 4.3 | 4.6 | 4.4 |
| 146 | 0.2 | 0.3 | 0.2 | 0.5 | 0.3 |
| 147 | <0.1 | 0.2 | <0.1 | <0.1 | <0.1 |
| 148 | <0.1 | <0.1 | <0.1 | 0.2 | 0.2 |
| 149-151 | <0.1 | <0.1 | <0.1 | <0.1 | <0.1 |

(2-furfurylthiol is discussed in the text as "exclusively 5 times labeled (m/z 119)" but its Table 3
block is not printed in the text layer; the pdftotext output lists only the four blocks above. Table 4
carries FFT.)

### Table 4. "Proportion of Isotopomers Formed from [13C5]Xylose, Cysteine, and Thiamin at pH 4-7"

Footnote a: unlabelled (%); b: labelled (%). Blank cells = compound not detected at that pH.

| no. | compound | m/z analysed | pH 4.0 12C / 13C | pH 5.0 12C / 13C | pH 6.0 12C / 13C | pH 7.0 12C / 13C | labelled C atoms |
|---:|---|---|---|---|---|---|---:|
| 1 | 4,5-dihydro-2-methyl-3(2H)-furanone | 100; 105 | >99 / <1 | >99 / <1 | >99 / <1 | >99 / <1 | 5 |
| 2 | 3-mercapto-2-butanone | 104; 108 | >95 / <5 | >98 / <2 | >98 / <2 | >95 / <5 | 4 |
| 3 | 2-furaldehyde | 96; 101 | <1 / >99 | <1 / >99 | — | — | 5 |
| 4 | **2-methyl-3-furanthiol** | 114; 119 | **75 / 25** | **85 / 15** | **90 / 10** | **80 / 20** | 5 |
| 5 | **3-mercapto-2-pentanone** | 118; 123 | **77 / 23** | **86 / 14** | **90 / 10** | **88 / 12** | 5 |
| 6 | **2-mercapto-3-pentanone** | 118; 123 | — | — | **6 / 94** | **<5 / >95** | 5 |
| 7 | **2-furfurylthiol** | 114; 119 | **<1 / >99** | **<1 / >99** | **<1 / >99** | — | 5 |
| 8 | 4,5-dihydro-2-methyl-3-furanthiol | 116; 121 | >99 / <1 | >99 / <1 | >99 / <1 | >99 / <1 | 5 |
| 9 | 2-methyl-3-(methylthio)furan | 128; 129; 133; 134 | 35 / 3; 8; 54 | 39 / 3; 8; 50 | 40 / 2; 9; 49 | — | 1; 5; 6 |
| 10 | 4,5-dihydro-2-methyl-3(2H)-thiophenone | 116; 121 | — | — | >99 / <1 | >98 / <2 | 5 |
| 11 | 3-methyl-1,2-dithian-4-one | 148; 153 | — | — | 91 / 9 | 92 / 8 | 5 |
| 12 | 3-acetyl-1,2-dithiolane | 148; 153 | — | — | >99 / <1 | >99 / <1 | 5 |
| 13 | 5-(2-hydroxyethyl)-4-methylthiazole | 143; 149 | >99 / <1 | >99 / <1 | >99 / <1 | >99 / <1 | 6 |

Text: "[13C5]4, 10-25%; [13C5]5, 10-23%". Sum of the labelled cluster in Table 3 (m/z 118-121 for MFT:
25.3 / 15.2 / 9.8 / 20.5 %; m/z 122-126 for 3-MP: 23.3 / 13.6 / 9.9 / 12.5 %) reproduces Table 4.

### Table 5. "Proportion of Isotopomers from the Reaction of [13C5]Xylose, Cysteine, and Thiamin (pH 5) at Different Precursor Concentrations"

Footnote a: "Precursor concentrations: A, xylose (0.3 m), cysteine (0.1 m), thiamin (0.1 m); B, xylose
(0.15 m), cysteine (0.05 m), thiamin (0.05 m)." (m = mol/L; B = this paper, A = ref 9.)

| no. | compound | m/z analysed | A unlabelled : labelled (%) | B unlabelled : labelled (%) |
|---:|---|---|---|---|
| 4 | 2-methyl-3-furanthiol | 114; 119 | 54 : 46 | 85 : 15 |
| 5 | 3-mercapto-2-pentanone | 118; 123 | 60 : 40 | 86 : 14 |

Text on Zeiler-Hilgart (ref 26), quoted numbers only: thiamin 120 mmol/L at pH 5.7 produced 13x more
MFT (45 vs 3.4 µg/L) than ribose + cysteine (both 120 mmol/L) under the same conditions; 1500x more at
meat-level concentrations (ref 27).

## 4. What the repo could take

Nothing here is absolute; there are no FIT rows in the mol % sense. Everything below is a within-study
ratio or a sign, at 145 C / 20 min, 0.5 M phosphate, xylose 150 / cysteine 50 / thiamin 50 mmol/L.

### 4.1 The pH shape of each thiol (Table 2, all sources, normalised to pH 5.0 = 1.00)

| compound | pH 4.0 | pH 5.0 | pH 5.5 | pH 6.0 | pH 7.0 |
|---|---:|---:|---:|---:|---:|
| MFT (total) | 1.24 | 1.00 | 1.10 | **1.60** | 1.16 |
| FFT | 1.17 | 1.00 | 0.99 | **0.50** | **0** |
| furfural | 1.32 | 1.00 | 1.04 | **0** | **0** |
| 3-mercapto-2-pentanone | 1.03 | 1.00 | 1.04 | **0.44** | **0.085** |
| 2-mercapto-3-pentanone | 0 | 0 | 0 | 19 (abs.) | 5 (abs.) |
| 4,5-dihydro-2-methyl-3-furanthiol | 5.3 | 1.00 | 0.95 | 0.42 | 0.16 |
| 4,5-dihydro-2-methyl-3(2H)-thiophenone | 0 | 0 | 0 | 372 (abs.) | 470 (abs.) |
| MFT / FFT (area ratio) | 0.96 | 0.91 | 1.02 | 2.9 | inf |

### 4.2 The XYLOSE-derived share — the quantity comparable to the model's pentose lane

Area x labelled fraction (Table 2 x Table 4). pH 5.5 has no label run. Uncertainty: >= 23 % on the area,
unstated on the fraction; treat as +-30-40 %.

| compound | pH 4.0 | pH 5.0 | pH 6.0 | pH 7.0 | note |
|---|---:|---:|---:|---:|---|
| MFT from xylose (area x 10^6) | 415 x 0.25 = **104** | 336 x 0.15 = **50** | 539 x 0.10 = **54** | 391 x 0.20 = **78** | rel. to pH 5: 2.1 / 1.0 / 1.1 / 1.6 — a shallow U, NOT the rise the total MFT shows at pH 6 |
| MFT from thiamin | 311 | 286 | 485 | 313 | the pH-6 peak of total MFT is a THIAMIN peak |
| 3-MP from xylose | 145 x 0.23 = **33** | 141 x 0.14 = **20** | 62 x 0.10 = **6** | 12 x 0.12 = **1.4** | rel.: 1.7 / 1.0 / 0.31 / 0.07 — monotone fall, ~24x from pH 4 to 7 |
| 3-MP from thiamin | 112 | 121 | 56 | 11 | |
| 2-MP (94-95 % xylose) | 0 | 0 | 18 | 5 | appears only at pH >= 6, entirely pentose-derived |
| FFT (>99 % xylose) | 431 | 368 | 185 | 0 | fall 2.3x from pH 4 to 6, gone at 7 |
| furfural (>99 % xylose) | 208 | 158 | 0 | 0 | gone at pH 6 while FFT is still 185 — furfural is consumed faster than it is made at pH 6, or not made |
| xylose-MFT / FFT | 0.24 | 0.14 | 0.29 | inf | pentose-lane MFT:FFT rises with pH |

### 4.3 Directional claims with numbers (the authors' words where possible)

- "2-furaldehyde (3), 2-furfurylthiol (7), and 2-methyl-3(methylthio)furan (9) were formed only at
  acidic pH, and no traceable amounts were detected at pH 7.0": FFT 431 -> 368 -> 364 -> 185 -> 0;
  furfural 208 -> 158 -> 165 -> 0 -> 0. **The furfural + H2S -> FFT lane is switched off between pH
  5.5 and 7.** The model's FFT step must carry a pH factor that reaches ~0 by pH 7 at 145 C.
- "Peak areas decreased with increasing pH value for 8 and 3-mercapto-2-pentanone (5)": 145 -> 12
  (12x, total) and 33 -> 1.4 (24x, xylose-derived) from pH 4 to 7.
- "The pH had no obvious influence on the formation of ... 2-methyl-3-furanthiol (4)": total MFT stays
  within 336-539 over pH 4-7 — but this is thiamin-buffered; the xylose share varies 2x (§4.2).
- "Compounds 6, 11, 12, and ... (10) were detected only when the reaction was carried out at pH 6.0 and
  7.0": 2-mercapto-3-pentanone, the dithianone, the acetyldithiolane and the thiophenone are pH >= 6
  species, which the authors attribute to "the stronger involvement of hydrogen sulfide at higher pH".
  For the model this is a measured H2S partition shift with pH — the same direction as Whitfield 2001's
  thiophenone/trithiolane sinks at pH 6.5.
- "The pH during the reaction had no major influence on the isotopomer distribution": MFT 25/15/10/20 %
  xylose, 3-MP 23/14/10/12 % xylose — route mix moves at most 2.5x across three pH units, far less than
  the compound levels.
- Meynier & Mottram 1995 comparison (text): there, MFT, 3-MP and FFT "were highest at pH 4.5 and
  decreased with increasing pH"; here only 5 and 7 follow, and "the peak areas for 4 remained
  relatively constant" — the difference being the thiamin.
- Table 5: doubling all three precursors at pH 5 moves the xylose share of MFT from 15 % to 46 % and of
  3-MP from 14 % to 40 % (already B10.1 in the repo).

### 4.4 Roles already declared

Table 2 HOLD-OUT (directional); Table 4 FIT; Table 5 HOLD-OUT. This dossier adds the derived
xylose-only pH shape (§4.2) as the row that should be compared to the model's pentose lane rather than
the raw MFT total, and the furfural/FFT switch-off between pH 5.5 and 7 as the sharpest claim for the
FFT step.

## 5. Caveats

1. **Peak areas, TIC, SPME, no internal standard, no response factors**: never levels, and the authors
   themselves warn that "the pH might have an effect on the partition coefficient". Thiols (pKa ~ 10)
   and furfural are neutral over pH 4-7, so their partition is probably pH-insensitive; the caveat
   bites for anything ionisable.
2. **Thiamin at 50 mM** is a second H2S and MFT source. 75-90 % of the MFT and 77-90 % of 3-MP are
   thiamin-derived at this concentration; only the xylose-derived shares (§4.2) speak to the model's
   pentose lane, and they carry two stacked uncertainties. Thiamin also consumes / supplies H2S in
   ways the pentose-cysteine model does not represent.
3. **Xylose, not ribose**; 3:1:1 ratio; 145 C / 20 min — the Hofmann 1998 frame, so the temperature and
   time match the fit panel but the sugar does not.
4. **pH 5.5 has no isotope run**; the xylose share at 5.5 is interpolated at best.
5. **Zeros are "no traceable amounts"** with no stated LOD; a 0 in Table 2 should be read as
   "< a few x 10^6 area units", not as exactly zero.
6. **Table 3's FFT block is not in the text layer** (only the four blocks above). The text and Table 4
   carry FFT as >99 % labelled at pH 4-6.
7. Post-reaction pH is not reported; 0.5 M phosphate should hold it, but at 145 C in Teflon this is
   assumed, not measured.
