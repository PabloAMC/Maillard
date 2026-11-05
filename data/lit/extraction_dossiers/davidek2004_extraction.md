# Davidek, Kraehenbuehl, Devaud, Robert & Blank 2004 — EXTRACTION
DOI: 10.1021/ac048925a (printed on p. 140; Crossref 2026-11-05: same title, five authors, Analytical Chemistry 77(1), 140-147, issued 2004-11-30. **Matches the brief.** The file stem says 2004, the issue says 2005.)
### *"Analysis of Amadori Compounds by High-Performance Cation Exchange Chromatography Coupled to Tandem Mass Spectrometry"* (Anal. Chem. 77:140, issue dated 1 Jan 2005, published on the web 30 Nov 2004)
### A METHODS paper: HPCEC (IonPac CS17) with MS/MS or ECD for hexose- and pentose-derived Amadori compounds; MS/MS fragments of nine glucose- and eight xylose-derived Amadori compounds. **One kinetic series: xylulosylglycine (Xyl-Gly) from D-xylose + glycine 1.0 + 1.0 mol/L in 0.2 M phosphate, label pH 6, 90 °C, 0-120 min (Fig. 8B), duplicate.** Food levels: Fru-Glu in dried tomato, Fru-Pro in dried apricot. **No cysteine.**

**Read 2026-11-05 from `data/articles/davidek2004.pdf`** (8 pp., Acrobat Distiller, born-digital, clean text layer).
Prose and tables from `pdftotext -layout`; **Fig. 8B digitised from a 400 dpi render** (scratchpad `m1c/f8.png`,
`m1c/f8b.png`). Read error on Fig. 8B: **+-0.0015 mol/mol** (y), +-1 min (x). Reader: one reader (Claude, first reader).
No value here is proposed for a ledger, band or fit row; a second reader is required before any is.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure by
me · **[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[Z]** derived by me, never
printed · **[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Tomas Davidek, Karin Kraehenbuehl, Stéphanie Devaud, Fabien Robert, Imre Blank* (Nestlé Research Center, Lausanne) |
| Venue | Analytical Chemistry **77**(1), 140-147 (1 Jan 2005); received 23 Jul 2004, accepted 12 Oct 2004, web 30 Nov 2004 |
| DOI | **10.1021/ac048925a** (printed) |
| Same laboratory | `davidek2003_extraction.md` (the anion-exchange predecessor, ref. 20 here), `davidek2002_extraction.md` (ref. 19) |

## 1. What this paper answers for the model

| model need | what Davidek 2004 gives | tag |
|---|---|---|
| a cysteine Amadori (CAMA) assay or level | **none**; glycine and seven other amino acids, no cysteine | -- |
| pentose Amadori formation and loss, one pot | **Xyl-Gly: 0.104 mol/mol xylose at 15 min, falling to 0.055 at 120 min** (90 °C, label pH 6, 1.0 M each) | [F] |
| whether a pentose Amadori survives the analytical chain | on HPAEC (the 2003 method): "failed due to its rapid on-column decomposition" [V]; on HPCEC: samples analysed within 8 h or frozen at -20 °C, "no degradation of Xyl-Gly was observed" [Q] (no numbers) | [Q] |
| tautomers of a pentose Amadori compound | synthetic Xyl-Gly in D2O: **alpha / beta / open-chain keto = 30 / 40 / 30** (1H NMR, H-C1 and H-C3 integrals); elutes as ONE symmetric peak on CS17, one MRM transition | [M] / [Q] |
| sugar quantification | not possible on this column: "parent sugars ... are not retained on cation exchange columns" [Q] | [Q] |

**Answer.** Evidence only, as the M1c case says. Two items touch the case without deciding anything: (i) a pentose
Amadori compound's tautomers (30 % open chain) co-elute and share one MRM transition here, which bears on, but does
not size, the P4a idea of a CAMA isomer outside an MRM window; (ii) a handling-stability statement for a pentose
Amadori compound ("no degradation" within 8 h), which is the kind of check R-cal asks for, made on glycine's Amadori
compound, not on cysteine or xylose, and printed without numbers.

## 2. Conditions [M]

| item | as printed |
|---|---|
| Xyl-Gly kinetics | D-xylose 3.75 g (25 mmol) + glycine 1.88 g (25 mmol) in 25 mL 0.2 mol/L phosphate, **pH 6**, **90 °C, 2 h**; 1 mL aliquots at intervals, PVDF 0.45 um, diluted into the calibration range, HPCEC-ECD; **"All samples were prepared in duplicate"**. [Z] 1.0 mol/L each (solid volume ignored; the true volume is larger, ~27 mL, so ~0.9 mol/L) |
| pH scale | label pH of the buffer; no measured pH; the drift at 1 M reactants in 0.2 M phosphate is not reported |
| Vessel, atmosphere | not printed |
| HPCEC-ECD | Dionex DX500, AS-50 (10 uL loop), GP-50; IonPac CS17 250 x 2 mm + CG17 50 x 2 mm; water 0.25 mL/min, 20 min; 150 mmol/L formic acid wash 10 min; post-column NaOH 300 mmol/L at 0.5 mL/min; gold electrode, triple or quadruple waveform (ref. 19) |
| HPCEC-MS/MS | Agilent 1100 + Q-Trap; CS17, isocratic 0.1 % formic acid, 0.25 mL/min, 30 °C; ESI+, TIS 450 °C, 5000 V; MRM [M+H]+ -> [M+H-H2O]+ at CV 20 V |
| Xyl-Gly standard | xylose 270 g + water 45 mL, 90 °C; glycine 33.8 g + sodium disulfite 40.5 g, 90 °C 20 min; Dowex 50WX4; 19 g, 22 % yield, **98 % purity by NMR**; hygroscopic, stored under argon |
| Calibration | external, fresh solutions; ranges for Xyl-Gly **not printed**. Fru-Glu linear to 300 mg/L (r^2 0.9999), Fru-Pro (r^2 0.9979), 1-300 ng injected |
| LOD / LOQ | **not printed** for any compound [!] |
| Retention reproducibility | < 5 % (Fru-Phe), < 3 % (others) |
| Recoveries | Fru-Glu 83 % (dried tomato), Fru-Pro 93 % (dried apricot) |

## 3. Quantitative results

### Fig. 8B (p. 146): Xyl-Gly, mol per mol xylose, 90 °C, pH 6 [F]

Axis: 0 -> 0.15 mol/mol, 4420 px per mol/mol at 400 dpi; x 6.77 px/min. Error bars (duplicates) visible at 60 and
120 min only; elsewhere within the marker.

| t, min | 0 | 15 | 30 | 45 | 60 | 90 | 120 |
|---|---|---|---|---|---|---|---|
| Xyl-Gly, mol/mol | 0 | **0.104** | 0.092 | 0.082 | 0.073 (+-0.002) | 0.062 | **0.055 (+-0.004)** |
| Xyl-Gly, mol/L [Z] (x 1.0 M) | 0 | 0.104 | 0.092 | 0.082 | 0.073 | 0.062 | 0.055 |

Text: "reaching ~12 mol %, followed by a slow decline" [V, abridged]; abstract: "~12 mol %, 15 min". **The plotted
15-min point is 0.104 (10.4 %), not 12 %** [!]. The text's 12 % may refer to glycine rather than xylose, a different
basis, or a rounding; the paper does not say. The axis label says per mol xylose.

Derived [Z]:
- Formation: >= 0.104 / 15 = **6.9e-3 mol/mol/min** average over 0-15 min (a floor; the peak may lie before 15 min).
- Net loss after the peak: ln(0.104 / 0.055) / 105 min = **6.1e-3 /min** (15-120 min); 4.5e-3 /min over 60-120 min.
  These are net (formation continues), so the true first-order loss of Xyl-Gly is faster.

### Food levels [M]
- Fru-Glu in dried tomatoes: **1536 mg/100 g** (recovery 83 %); Eichner et al. found 676 and 3788 mg/100 g in tomato
  powder [C].
- Fru-Pro in dried apricots: **206 mg/100 g** (recovery 93 %).

### Tables [M]
- Table 1: CS17 retention times (min) of nine Fru-AAs in water / 15 mmol/L formic acid / 90 % formic + 10 % ACN, e.g.
  Fru-Gly 4.28 / 4.15 / 3.91, Fru-Phe 8.22 / 7.92 / 5.34; dead volume 2.50 min.
- Tables 2-3: MS/MS fragments of nine glucose- and eight xylose-derived Amadori compounds. Xyl-Gly [M+H]+ 208 -> 190,
  172, 154, 144, 126, 97, 88, 76; Xyl-Gly co-elutes with Fru-Gly and is told apart by the 30-unit shift. The
  [M+H-3H2O-H2CO]+ fragment is characteristic of hexose Amadori compounds and absent for xylose ones [Q].

## 4. What the paper can and cannot test in the model

| model part | can test? | why |
|---|---|---|
| CAMA formation (`r_cama_form`) or decay | **no** | glycine, not cysteine |
| a pentose Amadori pool's time scale (glycine) | as context only | Xyl-Gly peaks <= 15 min at 90 °C and loses ~half by 120 min; a glycine analogue |
| P4a (a second CAMA isomer outside Zhang 2025b's MRM) | no | it shows tautomers of Xyl-Gly co-eluting under one transition on a different column; it does not size or exclude a CAMA isomer on Zhang's method |
| R-cal (handling loss) | no | a qualitative "no degradation" for Xyl-Gly within 8 h; no numbers; not cysteine or xylose |
| free sugar | no | not retained on the column |

## 5. Conflicts with existing records

None in kind. Internal: 10.4 % (Fig. 8B) vs "~12 mol %" (text, abstract) for the Xyl-Gly peak [!]. No other
Xyl-Gly time course is on disk to compare.

## M1c case criterion (frozen bd641ae8)

The frozen case (`docs/reference/cases/M1c_thiol_supply_structure.md` on `origin/m1c-case`, sec. 4.3) lists davidek2004
as **evidence only (glycine, not cysteine)**, with no decisive measurement named and no band role.

**Verdict: NOT MET (no criterion applies; evidence only, as declared).**
- No cysteine, no free pentose, no TTCA or CAMA is measured.
- The one kinetic series is Xyl-Gly from 1.0 + 1.0 mol/L xylose + glycine, 0.2 M phosphate, label pH 6, 90 °C:
  0.104 mol/mol at 15 min falling to 0.055 at 120 min (duplicates).
- What it adds as evidence: the predecessor HPAEC method destroys a pentose Amadori compound on-column; on HPCEC,
  Xyl-Gly's three forms (alpha / beta / open = 30 / 40 / 30) elute as one peak under one MRM transition; and Xyl-Gly
  was stable for 8 h of sample handling (no numbers).
