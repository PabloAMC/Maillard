# Yiltirak, Balagiannis, Koek, Koch & Elmore 2026 (10.1016/j.foodres.2026.118600) — single-paper extraction, programme step R0, 2026-09-06

**Source files:** `data/articles/Yiltirak2026.pdf` (main text, 11 pp., born-digital Elsevier PDF, CC BY) and
`data/articles/Yiltirak2026-supplementary.docx` (the publisher's `mmc1.docx`, 27.5 MB, five embedded PNGs).
Read method: the PDF text layer read directly (pp. 1-9); the docx read out of `word/document.xml` with a
zipfile + ElementTree pass (all three tables re-typed below from that XML, not from a rendering).
Both files are gitignored; this dossier is the tracked record.

---

## 0. PAPER IDENTITY — MATCHES THE BUNDLES' EXPECTED IDENTITY

| field | value |
|---|---|
| Authors | Suleyman Yiltirak (a), Dimitris P. Balagiannis (a), Jan Koek (b), Jens Koch (c), J. Stephen Elmore (a, *) |
| Affiliations | (a) Dept. of Food and Nutritional Sciences, University of Reading; (b) Foods Innovation Centre Unilever, Wageningen; (c) Symrise AG, Holzminden |
| Title | **"Formation of aroma-potent thiols in heat-stable vegetable oil-in-water emulsion models containing ribose and cysteine"** |
| Venue | *Food Research International* **231** (2026) 118600 |
| DOI | **10.1016/j.foodres.2026.118600** |
| Dates | Received 2 Oct 2025; revised 29 Jan 2026; accepted 2 Feb 2026; online 3 Feb 2026 |
| Licence | CC BY 4.0 |
| Data availability | "Data will be made available on request." |

Correct file for the four hold-out bundles `mp_holdout_ribose_cysteine_buffer_{100C_4h,110C_2h,120C_1h,130C_30min}_Yiltirak2026`.

---

## 1. ONE-PARAGRAPH VERDICT

The paper is what the bundles said it is: a **stable-isotope-dilution** (MFT-d3, FFT-d2, PFBBr derivatisation,
HS-SPME-GC-MS in SIM, matrix-matched five-level calibration 0.2-10 ug/L) measurement of MFT, FFT and
3-mercapto-2-butanone in four model systems (buffer; buffer + 10 % canola oil; buffer + 3 % sucrose ester
PS750; the full emulsion) at **four time-temperature combinations chosen as "equivalent" cooks**
(100 C / 4 h, 110 C / 2 h, 120 C / 1 h, 130 C / 0.5 h), n = 3. **All eight scored hold-out values are
Table S3 exactly (sec. 6).** What the PDF adds that the bundles did not carry, and what programme step R1
records as a physical-state input: **3 mL of sample in a 20 mL PTFE-capped Duran tube, heated under the air
it was closed with, flushed with argon only AFTER heating, buffer made in tap water.** The headspace holds
about twice as much O2 as there is cysteine (sec. 5). The other 40 numbers in Table S3 are within-lab
contrasts (matrix at fixed cook; cook at fixed matrix) that no scored row uses today (sec. 4).

---

## 2. METHODS, VERBATIM (the passages the bundles and the R1 vessel block rest on)

**2.1 Chemicals** (p. 2): "L-Cysteine (>=97%), D-(-)-ribose (>=98%), 2-methyl-3-furanthiol (99%),
2,3,4,5,6-pentafluorobenzyl bromide (99%) and Tween 20 ... were supplied by Sigma-Aldrich Co. (Gillingham,
UK). 2-Furfurylthiol (99%) and 3-mercapto-2-butanone (99%) were purchased from Tokyo Chemical Industry UK
Ltd. 2-(Methyl-d3)-3-furanthiol (MFT-d3) and 2-furfuryl-alpha,alpha-d2-thiol (FFT-d3) were purchased from
aromaLAB GmbH (Martinsried, Germany)." ... "**Potassium phosphate buffer (0.5 M, pH 5.5) was made using
potassium dihydrogen phosphate and dipotassium hydrogen phosphate in tap water.**"
(The SI Figure S1 caption and sec. 2.6 call the FFT standard FFT-d2; the chemicals list prints "FFT-d3"
once. Two deuteriums on the CH2 is the compound; the ions listed in 2.6, m/z 83/181/296 for FFT-d2 against
81/181/294 for FFT, are consistent with d2.)

**2.3 Preparation of reaction model systems** (pp. 2-3): "Four different models were prepared to understand
the effect of interface with and without oil on the formation of potent thiol compounds. The models were:
i) emulsion system (87% buffer, 3% sucrose ester, 10% canola oil), ii) buffer/oil system (90% buffer, 10%
canola oil), iii) buffer/emulsifier (97% buffer, 3% sucrose ester), iv) buffer system (100% buffer). All
buffer phases were prepared by dissolving ribose (25 mM) and cysteine (25 mM) in potassium phosphate buffer
(0.5 M, pH 5.5). The emulsification procedure was applied to all systems."

**2.4 Heating of samples** (p. 3): "After emulsification, model systems (3 mL) with 3% sucrose ester (PS750)
were transferred to 20-mL Duran tubes with screw-caps lined with PTFE gaskets (Duran Wheaton Kimble,
Wertheim, Germany) and they were heated at 100 C/4 h, 110 C/2 h, 120 C/1 h and 130 C/0.5 h in a heating
block equipped with a magnetic stirrer (RS900 reaction station; Electrochemical, Cambridge, UK). These
temperature combinations were selected as practical experimental conditions to represent slow and fast
cooking techniques for creating meat aroma (Miao et al., 2024). PTFE-covered magnetic stir bars (12 mm x
4.5 mm) were used to help maintain the one-phase structure during heating. After heating, the tubes were
cooled immediately to room temperature in an ice bath to stop the thermal reaction. Subsequently, they were
flushed with Pureshield Argon (B.O.C., UK) to prevent further oxidation and stored at -20 C until analysis.
Each heating condition was conducted three times. All systems without precursors were prepared and heated
under the same conditions and used as blanks for aroma analysis."

**2.5 Derivatization** (p. 3): "An aliquot of sample (30 uL) was mixed and vortexed with 10 uL of 20 ug/L
MFT-d3 and FFT-d2, and extracted with 6 mL ice-cold NaOH (1 M) in a 20-mL vial with metal screw-cap and
septum, suitable for solid-phase micro-extraction (SPME). Then, derivatization was performed by adding 50 uL
of pentafluorobenzyl bromide (PFBBr) solution (20 uL of PFBBr in 50 mL redistilled ethanol), vortex mixing
and resting at room temperature for 10-15 min to allow the reaction to proceed. Tartaric acid (0.5 g) was
added to reduce the pH to between 4 and 5, to stop the derivatization reaction as well as to optimize SPME
extraction and MS detection."

**2.6 Analysis** (p. 3): SPME 50/30 um DVB/CAR/PDMS, equilibration 80 C / 30 min, exposure 80 C / 20 min,
DB-5MS UI 30 m, EI 70 eV, **SIM**: "m/z 113, 181, 296 for MFT, m/z 116, 181, 297 for MFT-d3, m/z 81, 181,
294 for FFT, m/z 83, 181, 296 for FFT-d2 and m/z 72, 241, 181, 284 for 3M2B ... Analysis was performed in
triplicate for each sample. Concentrations were determined by using a matrix-matched calibration curve
prepared from an unheated emulsion containing equal amounts of internal standards as the samples, where
the concentration of the spiked MFT and FFT ranged from 0.2 ug/L to 10 ug/L. Limit of detection (LOD) and
limit of quantification (LOQ) were calculated based on the standard deviation of a linear response and a
slope (Fig. S3 and Table S1)."

**Method S1 (SI), hexanal**: dynamic headspace (2.5 mL sample, 20 uL of 5 ppm 3-furaldehyde, 60 C / 1 h
onto Tenax TA), TD-GC-MS in scan mode. "**Semi-quantification** was achieved by comparing the peak area of
the 3-furaldehyde internal standard with that of hexanal, assuming that both compounds have similar
response factors." Hexanal values live only in Figure S5 (an image in the docx; not digitised here).

---

## 3. TABLES, VERBATIM (re-typed from `word/document.xml`)

**Table S1.** "Linear regression equation, correlation coefficient (R2) limit of detection (LOD), and limit of
quantification (LOQ) values of derivatized MFT, FFT and 3M2B." (y = analyte response / internal-standard
response; x = ug/L)

| | Equation | LOD (ug/L) | LOQ (ug/L) | R2 |
|---|---|---:|---:|---:|
| MFT | y = 0.1336x - 0.0035 | 0.10 | 0.32 | 0.9919 |
| FFT | y = 0.8986x - 0.1472 | 0.28 | 0.85 | 0.9947 |
| 3M2B | y = 0.0071x - 0.0013 | 0.08 | 0.24 | 0.9959 |

**Table S3.** "The concentrations of MFT, FFT and 3M2B in model systems after heating for 4 h, 2 h, 1 h and
0.5 h at 100 C, 110 C, 120 C and 130 C, respectively." Mean +/- SD, ug/L. Footnote: "Different lowercase
letters indicate significant differences within the same column between temperatures for each model
separately (p < 0.05, Tukey HSD). Different uppercase letters indicate significant differences within the
same column between models for each temperature separately (p < 0.05, Tukey HSD)."

| system | cook | MFT (ug/L) | FFT (ug/L) | 3M2B (ug/L) |
|---|---|---|---|---|
| Emulsion | 100 C / 4 h | 3.73 +/- 0.31 aB | 2.74 +/- 0.30 aA | 0.66 +/- 0.11 aB |
| Emulsion | 110 C / 2 h | 2.83 +/- 0.08 bB | 2.81 +/- 0.21 aA | 0.71 +/- 0.08 aB |
| Emulsion | 120 C / 1 h | 1.28 +/- 0.07 cC | 1.28 +/- 0.09 cB | 0.62 +/- 0.02 aB |
| Emulsion | 130 C / 0.5 h | 1.29 +/- 0.08 cC | 2.17 +/- 0.15 bA | 0.24 +/- 0.05 bB |
| Buffer | 100 C / 4 h | **6.88 +/- 0.98 aA** | **1.28 +/- 0.08 cB** | 1.54 +/- 0.03 aA |
| Buffer | 110 C / 2 h | **3.29 +/- 0.06 bB** | **1.46 +/- 0.19 abB** | 1.16 +/- 0.08 bA |
| Buffer | 120 C / 1 h | **2.4 +/- 0.16 bcB** | **1.68 +/- 0.21 aA** | 1.21 +/- 0.21 bA |
| Buffer | 130 C / 0.5 h | **1.71 +/- 0.12 cB** | **1.62 +/- 0.11 abB** | 0.69 +/- 0.06 cA |
| Buffer+Oil | 100 C / 4 h | 0.72 +/- 0.10 aC | 1.19 +/- 0.06 abB | 0.46 +/- 0.03 bcC |
| Buffer+Oil | 110 C / 2 h | 0.75 +/- 0.01 aC | 1.32 +/- 0.04 aB | 0.78 +/- 0.07 aAB |
| Buffer+Oil | 120 C / 1 h | 0.62 +/- 0.13 abD | 0.96 +/- 0.15 bB | 0.65 +/- 0.10 abB |
| Buffer+Oil | 130 C / 0.5 h | 0.45 +/- 0.04 bD | 0.93 +/- 0.17 bC | 0.29 +/- 0.11 cB |
| Buffer+Emulsifier | 100 C / 4 h | 7.57 +/- 0.75 aA | 1.32 +/- 0.14 bB | 0.62 +/- 0.03 cB |
| Buffer+Emulsifier | 110 C / 2 h | 4.15 +/- 0.52 bA | 1.35 +/- 0.06 bB | 1.05 +/- 0.26 abAB |
| Buffer+Emulsifier | 120 C / 1 h | 4.15 +/- 0.24 bA | 1.87 +/- 0.06 aA | 1.39 +/- 0.19 aA |
| Buffer+Emulsifier | 130 C / 0.5 h | 2.81 +/- 0.28 cA | 1.80 +/- 0.290 aAB | 0.73 +/- 0.02 bcA |

Bold = the eight values the four hold-out bundles score. "1.80 +/- 0.290" is printed with three decimals in
the source; carried as printed. Table S2 (creaming index) is not chemistry and is not transcribed.

**Main-text corroboration** (sec. 3.3, p. 8; Fig. 6 is a bar chart of the same data): "The formation of MFT
was highest during heating at 100 C for 4 h, reaching 6.88 +/- 0.98 ug/L and 7.57 +/- 0.75 ug/L in buffer
and buffer/emulsifier system, respectively, and decreased by 75% and 63% after heating at 130 C for 0.5 h
(p < 0.05)." ... "the amount of MFT in the emulsion was 5.1-, 3.8-, 3.5- and 2.9-fold higher than in the
buffer/oil system after heating at 100 C, 110 C, 120 C and 130 C for 4 h, 2 h, 1 h and 0.5 h,
respectively (p < 0.05)." ... (p. 9) "the maximum formation of FFT was found in emulsion models as 2.74 +/-
0.3 ug/L and 2.81 +/- 0.21 ug/L after heating at 100 C for 4 h and at 110 C for 2 h" ... "Comparing heating
at 100 C for 4 h to 130 C for 0.5 h, FFT generation was reduced by 21% in both emulsion and buffer/oil
systems whereas 27% more FFT was formed in both buffer and buffer/emulsifier systems (p < 0.05)."
Checks: 6.88 x 0.25 = 1.72 (table 1.71); 3.73/0.72 = 5.2, 2.83/0.75 = 3.8, 1.28/0.62 = 2.1 (the text says
3.5 -- the 120 C emulsion/oil ratio does NOT reproduce from Table S3; the other three do), 1.29/0.45 = 2.9;
1.62/1.28 = 1.27 (buffer FFT, "27% more").

---

## 4. WITHIN-LAB CONTRASTS THE TABLE SUPPORTS (none scored today)

Temperature-time folds, same matrix (130 C / 0.5 h over 100 C / 4 h):

| matrix | MFT fold | FFT fold | 3M2B fold |
|---|---:|---:|---:|
| Buffer | 0.25 (falls, a>c) | 1.27 (rises, c<ab) | 0.45 |
| Buffer+Emulsifier | 0.37 | 1.36 | 1.18 |
| Emulsion | 0.35 | 0.79 | 0.36 |
| Buffer+Oil | 0.63 | 0.78 | 0.63 |

Matrix ratios at fixed cook (emulsion / buffer+oil; buffer+emulsifier / buffer):

| cook | MFT emul/oil | FFT emul/oil | MFT emulsifier/buffer | FFT emulsifier/buffer |
|---|---:|---:|---:|---:|
| 100 C / 4 h | 5.2 | 2.3 | 1.10 (ns, A/A) | 1.03 |
| 110 C / 2 h | 3.8 | 2.1 | 1.26 (A/B) | 0.92 |
| 120 C / 1 h | 2.1 | 1.3 | 1.73 (A/B) | 1.11 |
| 130 C / 0.5 h | 2.9 | 2.3 | 1.64 (A/B) | 1.11 |

Reading for the model: (i) under a cook chosen to be "equivalent" (Q10 ~ 2, i.e. an effective 57 kJ/mol
between 100 and 130 C), MFT in buffer falls 4x while FFT rises 1.27x -- MFT's net accumulation has an
effective barrier well BELOW the cook's, FFT's is at it. The core (B9) predicts the opposite signs on both
(MFT 64 -> 171, FFT 615 -> 214 ug/L across the same four bundles). (ii) Oil alone SUPPRESSES both thiols
(buffer+oil is the lowest arm everywhere: MFT 0.45-0.75 against 1.71-6.88 in buffer); the interface restores
and exceeds it. The authors read this as partitioning of MFT and intermediates to the interface, away from
the aqueous degradation. (iii) The emulsifier alone raises MFT modestly (1.1-1.7x) and leaves FFT flat.

---

## 5. THE PHYSICAL STATE THE BUNDLES DID NOT CARRY (R1 vessel block)

From sec. 2.4 and 2.1: **3 mL** of model system in a **20 mL** Duran tube, PTFE-gasket screw cap, magnetic
stirring, closed under laboratory air (nothing in the methods displaces it before heating; the argon flush
is explicitly "after heating ... to prevent further oxidation"), buffer in **tap water**.

Headspace air at closure: 20 - 3 = 17 mL. At 1 atm and 293 K, n = PV/RT = 1.013e5 x 1.7e-5 / (8.314 x 293)
= 7.07e-4 mol of air; x 0.209 = **1.48e-4 mol O2 = 0.148 mmol**. Dissolved O2 in 3 mL of air-saturated
water at 20 C (~0.27 mM) is 8e-4 mmol, negligible. Cysteine charged: 25 mM x 3 mL = **0.075 mmol**.
**O2 : cysteine = 2.0 mol/mol** (each O2 oxidises two thiols to one disulfide, so ~4 thiol equivalents per
cysteine). For comparison Hofmann & Schieberle 1998 heated 100 mL in a 200 mL autoclave with 3.3 mmol
cysteine: 0.87 mmol O2, ratio **0.26**; Bolton 1994 heated 33.3 g in a 125 mL vial with ~0.39 mmol cysteine:
~0.8 mmol O2, ratio **~2.0**. The core is near-unbiased on Hofmann (O2-poor) and over-predicts MFT ~20x on
Bolton and 9-100x on Yiltirak (both O2-rich). The oxidant is a state variable the network carries (`OX`)
but charges only from cystine; step R2(a) charges it from this block.

---

## 6. VERIFICATION OF THE FOUR HOLD-OUT BUNDLES (2026-09-06)

| bundle | MFT target | Table S3 | FFT target | Table S3 | verdict |
|---|---:|---:|---:|---:|---|
| ..._100C_4h_Yiltirak2026 | 6.88 +/- 0.98 | 6.88 +/- 0.98 | 1.28 +/- 0.08 | 1.28 +/- 0.08 | exact |
| ..._110C_2h_Yiltirak2026 | 3.29 +/- 0.06 | 3.29 +/- 0.06 | 1.46 +/- 0.19 | 1.46 +/- 0.19 | exact |
| ..._120C_1h_Yiltirak2026 | 2.40 +/- 0.16 | 2.4 +/- 0.16 | 1.68 +/- 0.21 | 1.68 +/- 0.21 | exact |
| ..._130C_30min_Yiltirak2026 | 1.71 +/- 0.12 | 1.71 +/- 0.12 | 1.62 +/- 0.11 | 1.62 +/- 0.11 | exact |

Precursors (25 mM + 25 mM), buffer (0.5 M potassium phosphate, pH 5.5, tap water), cook (four pairs) and
quantification class (stable-isotope dilution, matrix-matched calibration) all confirmed from the PDF. The
bundles' `conditions.buffer.provenance_class` moves from `repo_verbatim_methods_quote` to
`primary_source_pdf` in the same change as this dossier. Unit basis: ug/L of a 100 % aqueous 0.5 M buffer,
read as ppb at ~1 kg/L, no conversion -- unchanged.

---

## 7. WHAT IS AND IS NOT INGESTED, AND WHY

- **Ingested (already):** the 8 buffer-arm MFT/FFT values as hold-out targets.
- **Now recorded, not scored:** the other 40 Table S3 values (sec. 3) and the contrasts of sec. 4. The three
  lipid/emulsifier arms need an interface or partition term no lane carries; 3M2B is not a core species and
  was not spiked with a deuterated analogue (quantified against the MFT/FFT internal standards only).
- **Not ingested:** hexanal (Fig. S5; semi-quantitative by the authors' own statement); Tables 1-3 and S2
  (particle size, viscosity, creaming, surface tension -- physical characterisation of the emulsions).
- **Caveats carried:** tap-water buffer (trace metals uncontrolled); pH measured before heating only; the
  120 C emulsion/oil MFT ratio in the text (3.5) does not reproduce from the table (2.1); "FFT-d3" vs
  "FFT-d2" naming inconsistency in the source; the SI's Table S3 is the only place the numbers are printed.
