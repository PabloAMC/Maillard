# Feng et al. 2022 — extra-added GSH / Cys / Gly / Glu on the thermal degradation of the ribose–glutathione Amadori product (RG-ARP)

### Wave K6a per-paper extraction. 2026-08-29. **Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration was touched.**

### ★ THE CORPUS'S SECOND SULFUR TEMPERATURE LADDER, AND A DIRECT TEST OF KANG 2026's LOW LEG — SAME SPAN (100→120 °C), SAME pH (7.0), SAME HOLD (120 min), SAME CLASS OF CHEMISTRY (Cys released *in situ* from a pentose-derived MRI).

---

## ⚠️ HEADLINE — THE ANSWER TO THE KANG SWITCH-ON QUESTION, UP FRONT

> **Kang's "flat low leg" for MFT and FFT does NOT replicate at face value, and its most
> distinctive feature — free thiols being flatter than their own sulfur class — is
> CONTRADICTED. But once Feng's own measured precursor depletion is divided out, the flatness
> comes back and lands almost exactly on Kang's numbers.**

| quantity, **100 → 120 °C, pH 7, 120 min** | **Kang 2026** (TTCA/Cys) | **Feng 2022, ARP alone** | **Feng 2022, ARP + GSH** |
|---|---|---|---|
| **MFT fold-change** | **× 1.12** | **× 1.645** | **× 3.381** |
| **FFT fold-change** | **× 1.10** | **× 2.071** | **× 4.087** |
| **sulfur class fold-change** | **× 2.566** | **× 2.020** | **× 3.972** |
| **MFT ÷ its own sulfur class** | **0.44** ⭠ *the anomaly* | **0.81** | **0.85** |
| **FFT ÷ its own sulfur class** | **0.43** ⭠ *the anomaly* | **1.03** | **1.03** |
| apparent Eₐ(MFT), 2-point `[D]` | 6.9 kJ mol⁻¹ | **30.3 kJ mol⁻¹** | 74.3 kJ mol⁻¹ |
| apparent Eₐ(FFT), 2-point `[D]` | 5.8 kJ mol⁻¹ | **44.4 kJ mol⁻¹** | 85.9 kJ mol⁻¹ |
| apparent Eₐ(sulfur class), 2-point `[D]` | 57.5 kJ mol⁻¹ | **42.9 kJ mol⁻¹** | 84.1 kJ mol⁻¹ |

**Three findings, in order of confidence:**

1. **✅ THE CLASS-LEVEL LOW LEG REPLICATES WELL.** Kang's sulfur-class 100→120 response is
   ×2.566 (Eₐ 57.5); Feng's ARP-alone arm gives ×2.020 (Eₐ 42.9). Two independent
   Cys-releasing pentose systems agree on the total-sulfur temperature response of the low leg
   to within 27 % in fold-change and 15 kJ mol⁻¹ in apparent Eₐ. **The repo's low-leg sulfur
   slope is corroborated.**
2. **❌ KANG'S SPECIFIC THIOL ANOMALY DOES NOT REPLICATE.** In Kang, MFT and FFT run at
   **0.43–0.44 ×** their own class rate — they are the flattest things in the table. In Feng's
   closest-matched arm they run at **0.81 and 1.03 ×** their class rate; FFT simply *tracks*
   the class. **Feng does not see free thiols decoupling from the sulfur class over 100→120 °C.**
3. **⭐ BUT NORMALISED FOR PRECURSOR CONSUMPTION, FENG'S ARP-ALONE THIOLS GO FLAT AND MATCH
   KANG ALMOST EXACTLY.** Feng measures ARP depletion directly (Fig 2a): **93.50 % of the ARP
   is consumed by 120 min at 120 °C**, against ≈ 48 % at 100 °C — a **× 1.95 ± 0.16** ratio of
   moles-consumed `[D]`. Dividing that out:

   | | raw × | ÷ 1.95 consumed | **Kang raw ×** |
   |---|---|---|---|
   | **MFT per mole ARP consumed** | 1.645 | **× 0.84** | 1.12 |
   | **FFT per mole ARP consumed** | 2.071 | **× 1.06** | 1.10 |
   | sulfur class per mole consumed | 2.020 | **× 1.04** | 2.57 |

   **⇒ Per mole of precursor actually destroyed, Feng's MFT and FFT are FLAT (× 0.84 / × 1.06)
   over 100→120 °C — the same conclusion Kang reached, reached by a different route.** The raw
   × 1.6 / × 2.1 in Feng is almost entirely *extent*, not *selectivity*. **Kang's low leg is
   substantively vindicated; the disagreement is in whether the extent term was already divided
   out.** (This correction rests on a digitised 100 °C ARP endpoint and a nominal 20 mmol L⁻¹
   t = 0 — see §6c for the full uncertainty.)

**⚠️ AND THE THING THIS PAPER CANNOT DO:** Feng has **exactly two temperatures, both inside
Kang's low leg.** It can corroborate the low leg's *magnitude*. **It can never test for the
slope BREAK at 120→140 °C** — that is the whole of Kang's switch-on claim, and Feng is silent
on it. Two points determine a line with zero degrees of freedom.

**Provenance codes:** **[M]** measured and printed · **[C]** cited from another work ·
**[F]** fitted by the authors · **[D]** derived or digitised by this wave, never printed ·
**[NEG]** verified negative · **[!]** integrity flag.

---

## §0. IDENTITY — **CONFIRMED, NO MISMATCH**

| field | value | how verified |
|---|---|---|
| file on disk | `data/articles/feng2022.pdf`, **3,602,080 bytes**, 11 pages, PDF 1.4 | `ls`, `pdfinfo` |
| SHA-256 | `92b2380426c6776497a4b060c21ee319dbd50f44ee972e1b9d22fe13110a7afa` | `shasum -a 256` |
| **title** | ***"Promoted Formation of Pyrazines and Sulfur-Containing Volatile Compounds through Interaction of Extra-Added Glutathione or Its Constituent Amino Acids and Secondary Products of Thermally Degraded N-(1-Deoxy-D-ribulos-1-yl)-Glutathione"*** | p. 9095 title block + PDF `Title` metadata |
| **authors** | **Linhui Feng, Heping Cui, Pusen Chen, Khizar Hayat, Xiaoming Zhang\*, Chi-Tang Ho\*** | p. 9095 |
| affiliations | Feng, Cui, Chen, Zhang — **State Key Laboratory of Food Science and Technology, School of Food Science and Technology, Collaborative Innovation Center of Food Safety and Quality Control in Jiangsu Province, Jiangnan University, Wuxi 214122, P. R. China** · Hayat — **Dept. of Kinesiology, Nutrition, and Health, Miami University, Oxford, Ohio 45056, USA** · Ho — **Dept. of Food Science, Rutgers University, New Brunswick, New Jersey 08901, USA** | p. 9104 |
| ORCIDs | Zhang 0000-0002-1673-6026 · Ho 0000-0001-8273-2085 · Cui 0000-0002-0183-0095 | p. 9104 |
| **DOI** | **`10.1021/acs.jafc.2c02949`** | footer, every page |
| **journal / vol / pages / year** | ***J. Agric. Food Chem.* 2022, 70, 9095−9105** | running footer, every page |
| **dates** | **Received April 28, 2022 · Revised July 2, 2022 · Accepted July 6, 2022 · Published July 15, 2022** | p. 9095 |
| version | **published version of record**, © 2022 American Chemical Society; downloaded via Tsinghua Univ. 2022-11-05 | p. 9095 stamp |
| producer | Arbortext Publishing Engine → PDFlib+PDI 9.1.2p4 (C++/Win64), modified with iText 4.2.0 | `pdfinfo` |
| **funding** | **National Natural Science Foundation of China 32172330; National First-Class Discipline Program of Food Science and Technology JUFSTR20180204; Collaborative Innovation Center of Food Safety and Quality Control in Jiangsu Province** | p. 9104 |
| **competing interests** | *"The authors declare no competing financial interest."* | p. 9104, Notes |
| keywords | Amadori rearrangement product; glutathione; cysteine; pyrazines; sulfur-containing volatile compounds; furans | p. 9095 |

**⇒ Title, authors, journal, volume, pages, year and the two-temperature (100 / 120 °C) design
all match the expected identity exactly. No mismatch to flag.**

**⚠️ Supporting Information is NOT on disk.** The article's SI (`Table S1` calibration curves +
LOD + LOQ; `Tables S2, S3` OAVs; `Figure S1` NMR; `Figure S2` Rib+GSH totals; `Figure S3`
thiazoles) is a separate PDF that is not in `data/articles/`. See §11.

---

## §0b. ⚠️ THE μ → m RASTER CHECK (Amendment 4) — **DISCHARGED. FENG 2022 IS CLEAN.**

**Binding requirement satisfied. Every unit that carries a number used anywhere in this dossier
was verified against a raster render, not the text layer.** Method:
`pdftoppm -r 600/700/900 -png -f N -l N feng2022.pdf out`, cropped to the artefact and read as
an image, glyph by glyph. Table 3 is typeset **rotated 90°** on p. 9102 and was rotated back
with `sips -r 90` before reading.

| artefact | text-layer says | **raster says** | verdict |
|---|---|---|---|
| **Table 1 caption** (p. 9098) | `Volatile Compounds (μg/L)` | **`Volatile Compounds (𝛍g/L)` — bold-italic Greek mu, unmistakable, 900 dpi** | ✅ **NOT corrupted** |
| **Table 2 caption** (p. 9101) | `Volatile Compounds (μg/L)` | **`Volatile Compounds (𝛍g/L)`** | ✅ **NOT corrupted** |
| **Table 3 caption** (p. 9102, rotated) | `Pyrazines (μg/L)` | **`Pyrazines (𝛍g/L)`** | ✅ **NOT corrupted** |
| **internal standard** (p. 9097) | `2 μL of 0.018 μg/μL` | **`2 μL of 0.018 μg/μL 1,2-dichlorobenzene`** — three true mus in one clause | ✅ **NOT corrupted** |
| **column film thickness** (p. 9097) | `30 m × 0.25 mm × 0.25 μm` | **`30 m × 0.25 mm × 0.25 μm`** | ✅ **NOT corrupted** |
| **Fig. 2g y-axis** (p. 9099) | *(vector image, no text layer)* | **`Concentration (μg/L)`** | ✅ raster-only, mu present |
| **Fig. 2a–f, h y-axes** | *(image)* | **`Concentration of ARP / GSH / PCA / Cys-Gly / Cys / Gly (mmol/L)`** — genuine `mmol/L`, **not** a corrupted `μmol/L` (values 0–30 with 20 mmol/L loading confirm the scale) | ✅ raster-only |
| **Fig. 3a–f y-axes** (p. 9100) | *(image)* | **`Concentration (mmol/L)`** — genuine, values 0–0.5 against a 20 mmol L⁻¹ substrate | ✅ raster-only |
| **body text, p. 9102** | `0.04 μg/L`, `0.43 μg/L`, `0.08 μg/L` | **`0.04 μg/L`, `0.43 μg/L`, `0.08 μg/L`** | ✅ **NOT corrupted** |
| **body text, p. 9103** | `36.41 μg/L`, `2.39 μg/L`, `70.84 μg/L` | **`36.41 μg/L`, `2.39 μg/L`, `70.84 μg/L`** | ✅ **NOT corrupted** |
| **quantification sentence** (p. 9098) | `the concentration (μg/L)` | **`x and y represent the concentration (μg/L) and the peak area ratio…`** | ✅ **NOT corrupted** |

**Per-table basis of verification.**
- **Table 1** — raster, p. 9098 at 600 dpi: caption (900 dpi), both header rows, **all 6 furan
  rows, the furan subtotal, all 7 sulfur rows, the sulfur subtotal and the total**, plus both
  footnotes. Every raster reading matched the text layer to the last decimal. **Table 1 is
  verified on a raster basis end-to-end.**
- **Table 2** — raster, p. 9101 at 600 dpi: caption (900 dpi), both header rows, **all 7 furan
  rows + subtotal and all 9 sulfur rows** through 3-(methylthio)-2-butanone. Matched the text
  layer exactly. The sulfur subtotal / pyrazine block / total were read from the text layer and
  are **independently validated by the arithmetic audit in §3c, which closes to ±0.01 on every
  one of the 18 aggregates** — a self-check that fails on any digit error.
- **Table 3** — raster, p. 9102, rotated 90° and de-rotated: caption, header row (all 8 column
  labels), **both compound rows and the totals row, and both footnotes**. This was essential:
  the text layer serialises this rotated table into unusable interleaved fragments (see §3d).
- **Figures 2 and 3** — raster only; there is no text layer inside the plot areas.

**⇒ Unlike `Kang2026.pdf`, this ACS/Arbortext file embeds a real Greek `μ` throughout. No
μ→m correction is required anywhere in Feng 2022. The hazard was real and was checked; it did
not fire here.** ⚠️ Note the asymmetry: `Kang2026.pdf` (RSC-typeset) **is** corrupted, so any
cross-comparison of Feng's numbers with Kang's must use `kang2026_SI_extraction.md`'s corrected
units, not `Kang2026.pdf`'s text layer.

---

## §1. SYSTEM AND CONDITIONS `[M]`

### 1a. Preparation and purification of the RG-ARP

> *"Rib (0.2 mol/L) and GSH (0.1 mol/L) were mixed (3.04 g) in 50 mL of deionized water. Then,
> the mixture was adjusted to pH 7.5 using NaOH (6 mol/L) and heated at 80 °C under atmospheric
> pressure for 40 min in a magnetically stirred water bath. Further vacuum dehydration was
> performed at 80 °C for 10 min. The dehydrated solid was redissolved in 50 mL of deionized
> water to obtain the RG-ARP with a 58.00% yield."* — p. 9097

| item | value | anchor |
|---|---|---|
| precursors | **ribose 0.2 mol L⁻¹ + glutathione 0.1 mol L⁻¹** (2 : 1 sugar : peptide), 3.04 g total in 50 mL deionised water | p. 9097 |
| **preparation pH** | **7.5**, set with NaOH 6 mol L⁻¹ | p. 9097 |
| preparation heat | **80 °C, atmospheric pressure, 40 min**, magnetically stirred water bath | p. 9097 |
| dehydration | **vacuum, 80 °C, 10 min**; solid redissolved in 50 mL water | p. 9097 |
| **preparation yield** | **58.00 %** — defined as *"the ratio of the measured molar concentration of the RG-ARP to the initial molar concentration of GSH"* | p. 9097 |
| purification stage 1 | **Dowex 50WX8 ion-exchange resin, 26 × 300 mm, H⁺ form**; distilled water at 2 mL/min to elute unreacted Rib; then **ammonium hydroxide 0.10 mol L⁻¹ at 1 mL/min** to isolate the GSH + RG-ARP mixture | p. 9097 |
| purification stage 2 | **Waters 1525 HPLC + Alltech 3300 ELSD**, amide column **4.6 × 250 mm, 5 μm**, 1 mL/min, **40 °C**; mobile phase 0.10 % formic acid in water (A) / acetonitrile (B), **85 → 80 % B over 20 min**; ELSD evaporation tube 55 °C, N₂ 1.5 L/min; eluate evaporated then **lyophilised** | p. 9097 |
| **purity of the RG-ARP** | ★ **98.05 %** | p. 9097, Chemicals |
| identity | **N-(1-deoxy-D-ribulos-1-yl)-glutathione**, C₁₅H₂₅O₁₀N₃S, M = 439 Da; [M+H]⁺ calculated **440.1333**, measured **440.1296** (Δ ≈ 8.4 ppm) | p. 9097 |
| MS/MS fragments | m/z **422.1237** [M−H₂O+H]⁺ · **404.1131** [M−2H₂O+H]⁺ · **244.0829** [M−2H₂O−Cys-Gly+H]⁺ · **308.0929** [M−Rib+H]⁺ · **386.1018** [M−3H₂O+H]⁺ · **226.0738** [M−3H₂O−Cys-Gly+H]⁺ · **179.0514** [Cys-Gly+H]⁺ · **130.0534** [PCA+H]⁺ | p. 9097 |
| NMR | Bruker DRX 400 MHz, PABBO 5 mm probe, 50 mg in D₂O, D₂O signals as internal standard. **Five isomers** (1, 3′-α, 3′-β, 3-α, 3-β); δ_C **205.21 ppm** (C=O) proves the chain-keto form; δ_C 102.92–99.62 (C-2) / δ_H 3.31–3.21 for 3′-α/3′-β; δ_C 76.21 (C-1) / δ_H 4.26 and 2.96–2.87 for 3-α/3-β. **Spectra are in Figure S1 (SI, not on disk).** | p. 9097 |

### 1b. The model reaction — **the arms, the loadings, the vessel**

> *"The solutions of the RG-ARP (20 mmol/L) with or without GSH, Cys, Glu, or Gly (20 mmol/L),
> Rib with GSH (20 mmol/L), and GSH (20 mmol/L) were prepared and adjusted to pH 7.0 with NaOH
> (6 mol/L). Then, the solutions were injected into the pressure-resistant glass bottle and
> thermally stirred in an oil bath at 100 and 120 °C for different times (30, 60, 90, and 120
> min). All reaction products were immediately transferred to ice water to stop the reaction."*
> — p. 9097

| item | value | anchor |
|---|---|---|
| **substrate** | **RG-ARP, 20 mmol L⁻¹**, in every ARP arm | p. 9097 |
| **arms (7 distinct systems)** | ① **ARP alone** · ② **ARP + GSH** · ③ **ARP + Cys** · ④ **ARP + Glu** · ⑤ **ARP + Gly** · ⑥ **Rib + GSH** (comparator) · ⑦ **GSH alone** (comparator) | p. 9097 |
| **extra-added co-substrate loading** | ★ **20 mmol L⁻¹ in every case — equimolar with the ARP, 1 : 1** | p. 9097 |
| solvent | **water** (deionised); NaOH 6 mol L⁻¹ used only to set pH | p. 9097 |
| **⚠️ reaction volume** | **NOT STATED anywhere.** Only the 3 mL aliquot taken for GC and the 0.5 mL taken for α-dicarbonyl derivatisation are given. `[NEG]` | — |
| **vessel** | **"pressure-resistant glass bottle"** — sealed and pressure-rated, so the 120 °C runs are above the atmospheric boiling point and are genuinely pressurised | p. 9097 |
| **heating apparatus** | **oil bath, thermally stirred** | p. 9097 |
| **pH** | ★ **initial 7.0 only**, set with NaOH 6 mol L⁻¹ | p. 9097 |
| **⚠️ buffer** | ★ **NONE. No buffer of any kind is used in the reaction.** HEPES 1 mol L⁻¹ appears *only* in the off-line OPD derivatisation reagent, never in the reaction vessel. `[NEG]` | p. 9097–9098 |
| **⚠️ pH drift** | ★ **NEVER MEASURED. No final pH is reported for any arm, any temperature, any time.** `[NEG]` — see §8 | — |
| **temperatures** | ★ **100 and 120 °C. Two levels, no third.** | p. 9097 |
| **times** | **30, 60, 90, 120 min** at 120 °C in Tables 1–2; ★ **only 120 min at 100 °C in Tables 1–2** (the 100 °C time course exists in Figs 2g and 3a/3b — see §4) | Tables 1–3 |
| **replication** | **n = 3**, *"All experiments were executed in triplicate"* | p. 9098 |
| **quench** | ★ *"All reaction products were immediately transferred to ice water to stop the reaction."* | p. 9097 |
| statistics | mean ± SD, **SPSS 22.0**, p < 0.05, letters a–d in a row | p. 9098 |
| reagents | Rib, GSH, PCA, L-Glu, Cys-Gly, Gly, Cys, GO, MGO, 1,2-dichlorobenzene, 3,4-hexanedione, HEPES, formic acid, DTPA, methanol, OPD, acetonitrile — Sigma-Aldrich (Shanghai). 16 authentic flavour standards from Aladdin (Shanghai). **⚠️ No purities are stated for any reagent except the RG-ARP (98.05 %).** | p. 9097 |

**⚠️ The most important condition note for the repo.** The 100 and 120 °C columns of Tables 1 and
2 are **yields at a fixed 120 min hold in an unbuffered, sealed, pH-7-initial aqueous system with
no measured pH drift and no measured final pH.** They are *not* rate measurements, and the
system's pH at the moment of measurement is unknown at both temperatures. Every derived Eₐ in
§6 inherits both defects.

---

## §2. ANALYTICAL METHOD AND QUANTIFICATION BASIS — **REQUIRED SECTION**

### 2a. Volatiles: HS-SPME-GC-MS

| item | value | anchor |
|---|---|---|
| **internal standard** | ★ **1,2-dichlorobenzene**, **2 μL of a 0.018 μg μL⁻¹ solution in methanol** added to a **3 mL** sample ⇒ **36 ng into 3 mL = 12 μg L⁻¹ nominal in the liquid** `[D]` | p. 9097, raster-verified |
| vial | **15 mL**, PTFE/BYTL septum | p. 9097 |
| **SPME fibre** | **divinylbenzene/Carboxen/polydimethylsiloxane (DVB/CAR/PDMS), Supelco USA.** ⚠️ **film thickness NOT stated** `[NEG]` | p. 9097 |
| **extraction** | **60 °C for 20 min**, magnetic stirring, fibre in the headspace | p. 9097 |
| desorption | GC injector **250 °C, 10 min, splitless** | p. 9097 |
| **column** | **DB-Wax, 30 m × 0.25 mm × 0.25 μm**, J&W Scientific | p. 9097, raster-verified |
| oven programme | 40 °C (2 min) → 80 °C at 9 °C/min → hold 2 min → 150 °C at 6 °C/min → hold 3 min → 250 °C at 8 °C/min → hold 3 min | p. 9097 |
| carrier | helium, **1.5 mL/min** | p. 9097 |
| MS | Agilent 7890B GC / 5977B MSD, scan **m/z 35–450**, ion source **230 °C**, EI **70 eV** | p. 9097 |
| identification | RI vs **C₇–C₃₀ n-alkanes** on DB-WAX + **NIST 17** library + **authentic standards**. **Every single compound row in Tables 1, 2 and 3 carries the `S` flag** = *"volatile compounds authenticated by available standard compounds"* | Tables 1–3 footnote b |

### 2b. **THE QUANTIFICATION BASIS — one bolded sentence, as required**

> *"The quantification of volatile compounds was achieved by the limit of detection (LOD), the
> limit of quantitation (LOQ), and the calibration curves established through the available
> standard compounds (Table S1). x and y represent the concentration (μg/L) and the peak area
> ratio of the volatile standard compound to the internal standard, respectively."*
> — p. 9098, **raster-verified verbatim**

**⇒ FENG 2022's VOLATILE QUANTIFICATION IS ABSOLUTE, NOT SINGLE-IS SEMI-QUANT: every reported
compound has its own authentic external standard and its own measured calibration curve of
peak-area-ratio versus concentration, with the 1,2-dichlorobenzene internal standard serving
only as the ratio denominator, so no response factor is assumed equal to 1 anywhere.**

**⚠️ But it is absolute-in-principle and unauditable-in-practice.** The curves, slopes,
intercepts, R², LODs, LOQs, calibrated ranges and number of calibration levels are **all in
Table S1, which is in the Supporting Information and not on disk.** Consequences:
- I cannot check whether the reported 0.03–74.25 μg L⁻¹ values fall inside or outside the
  calibrated span — the same unquantified uncertainty that `kang2026_SI_extraction.md` §2c
  identified as the largest single risk in Kang's Table S4.
- I cannot run the internal-plausibility check that exposed Kang's impossible 2-methylthiazole
  curve. **Nothing here has been checked for that failure mode.**

### 2c. **THIS WAVE'S RULE, APPLIED EXPLICITLY**

The wave rule states: *a single-IS semi-quant ladder is still usable for activation energies as
within-study SHAPE, because a constant unknown response factor cancels in a ratio and therefore
in an Arrhenius slope.* **Feng 2022 clears that bar a fortiori — it is not semi-quant at all.**
Formally: for any compound *i*, the reported concentration is
`C_i = (A_i / A_IS - b_i) / m_i`, and the 120 °C / 100 °C fold-change is a ratio in which the
per-compound slope `m_i` cancels exactly, so **the Arrhenius slope is invariant to the
calibration slope even if that slope is wrong**, provided the intercept `b_i` is small relative
to the signal (unverifiable here — Table S1 is missing). **Fold-changes and apparent Eₐ are the
robust quantities. Use them.**

**What that does NOT license — spelled out:**
1. **Not absolute yields.** These are **headspace** SPME concentrations back-calculated through
   liquid-phase standards. Feng applies **no partition correction** and reports **no recovery**.
   The number 7.22 μg L⁻¹ for MFT is an operational quantity tied to *this* fibre, *this*
   60 °C / 20 min extraction and *this* 3 mL-in-15 mL headspace ratio. `[NEG]`
2. **Not cross-paper magnitude comparison.** Feng's MFT at 120 °C (7.22 μg L⁻¹, ARP alone) and
   Kang's (1.388 μg L⁻¹) differ 5.2-fold, but they were extracted on different fibres against
   different internal standards from different substrates at different loadings. **That ratio
   is meaningless and must never enter a benchmark.** Only the *within-paper* 120/100 ratios are
   comparable across the two papers.
3. **Not any compound whose response factor could plausibly change with matrix.** The
   ARP + GSH arm carries ~20 mmol L⁻¹ more solute, and thiols in particular are
   matrix-sensitive (oxidation on the fibre, metal catalysis, disulfide exchange). **Fold-changes
   should be compared *within* an arm — ARP-alone 120/100, or ARP+GSH 120/100 — never
   ARP+GSH-120 ÷ ARP-alone-100.**
4. **Not a rate.** See §6a.

### 2d. Non-volatiles: two separate, differently-grounded methods

| analyte set | method | **basis** | verdict |
|---|---|---|---|
| **RG-ARP, GSH, PCA, Cys-Gly, Cys, Gly** | 0.22 μm filtration → **UPLC−ESI−MS/MS**, Waters Synapt MALDI Q-TOF, **BEH C18 2.1 × 50 mm, 1.7 μm**, 100 % MeCN (A) / 0.10 % formic acid in water (B), B 98 → 60 % over 10 min, 0.3 mL/min, 45 °C, m/z 50–1000, scan 0.5 s, interscan 0.02 s, extractor 7 V, cone 20 V, capillary 3.5 kV, resolution ≥ 15 000, desolvation 400 °C, source 100 °C | ★ **ABSOLUTE — six printed calibration curves** (§2e) | **USE** |
| **GO, MGO** | OPD derivatisation → **HPLC-DAD 315 nm**, SunFire C18 4.6 × 150 mm 5 μm | ★ **ABSOLUTE — external calibration curves printed** (§2e) | **USE** |
| **1-DO, 3-DO** | same OPD/HPLC-DAD run | ⚠️ ***"quantitated using the internal calibration method due to the lack of standards"*** — **no authentic deoxyosone standards exist; the response factor is assumed** | **RATIO-ONLY** |

**⇒ ⚠️ THE DEOXYOSONE NUMBERS ARE THE ONE SEMI-QUANT DATASET IN THIS PAPER.** 1-DO and 3-DO are
quantified against the 3,4-hexanedione internal standard with an **assumed** response factor.
Their absolute mmol L⁻¹ values must not be used; their **ratios** (within a panel, across time,
across temperature) are sound, because the assumed factor cancels. This is exactly the case the
wave rule was written for, and it applies to 1-DO/3-DO **and to nothing else in this paper**.

**α-dicarbonyl derivatisation detail** `[M]`: OPD **5 g L⁻¹** + DTPA **11 mmol L⁻¹** dissolved in
**HEPES buffer 1 mol L⁻¹**; **0.5 mL derivatising agent + 0.5 mL reaction solution + 20 μL
internal standard (3,4-hexanedione, 0.96 mg mL⁻¹)**, incubated **25 °C, 12 h, in darkness**;
detection **315 nm**; mobile phase 0.10 % formic acid/water (A) + methanol (B) at 1.0 mL/min,
35 °C; gradient 0–10 min 5 % B, 10–20 min 5→30 % B, 20–40 min 30→40 % B, 40–60 min 100 % B.

### 2e. Every printed calibration curve `[F]` — **the six non-volatile curves + two dicarbonyl curves**

| # | analyte | calibration curve as printed | R² | anchor |
|---|---|---|---|---|
| 1 | **RG-ARP** | y = 2.00 × 10⁶ x + 9.65 × 10⁴ | **0.9926** | p. 9097 |
| 2 | **GSH** | y = 6.54 × 10⁵ x + 5.47 × 10⁴ | **0.9988** | p. 9097 |
| 3 | **PCA** (pyroglutamic acid) | y = 1.00 × 10⁶ x + 8.94 × 10⁵ | **0.9978** | p. 9097 |
| 4 | **Cys-Gly** | y = 2.00 × 10⁶ x + 1.10 × 10⁴ | **0.9992** | p. 9097 |
| 5 | **Cys** | y = 3.00 × 10⁵ x + 1.78 × 10³ | **0.9940** | p. 9097 |
| 6 | **Gly** | y = 8.45 × 10⁵ x + 9.45 × 10⁵ | **0.9983** | p. 9097 |
| 7 | **Glu** | y = 1.94 × 10⁵ x + 5.26 × 10⁴ | **0.9999** | p. 9097 |
| 8 | **GO** (glyoxal) | y = 4.00 × 10⁷ x − 2.91 × 10³ | **0.9998** | p. 9098 |
| 9 | **MGO** (methylglyoxal) | y = 3.00 × 10⁷ x − 4.27 × 10³ | **0.9997** | p. 9098 |

**⚠️ [!] Every slope is printed to at most three significant figures and five of the nine are
printed as a bare `n.00 × 10^k`** — i.e. rounded to 1–3 s.f. **Curves 3 and 6 have intercepts
comparable to their slopes** (PCA: intercept 8.94 × 10⁵ against slope 1.00 × 10⁶ ⇒ x-intercept
−0.894 in the concentration unit; Gly: 9.45 × 10⁵ / 8.45 × 10⁵ ⇒ x-intercept **−1.118**). If the
concentration unit is mmol L⁻¹, those intercepts are **≈ 10 % of the reported PCA values and
≈ 6 % of the Gly values** — a systematic offset that is not negligible. `[!]` **⚠️ No units are
stated for x on any of the nine curves, and no calibration range, no LOD, no LOQ and no number
of levels are given for any of them.** `[NEG]`

---

## §3. EVERY TABLE, CELL BY CELL

### 3a. **TABLE 1 — the ARP-alone ladder** `[M]`

**Caption, raster-verified verbatim:** *"Table 1. Volatile Compounds (μg/L) Identified in the
RG-ARP Model Reactions at Different Temperatures (100 and 120 °C)ᵃ"* — p. 9098.
**Footnotes, raster-verified verbatim:** *"ᵃResults are expressed as mean ± standard deviation
(n = 3). Different superscript letters (a−d) in the same row show significant differences
(p < 0.05). "-", not detected. ᵇRI, detected by the DB-WAX column with n-alkanes; MS, mass
spectra identified by the NIST 17 database; S, volatile compounds authenticated by available
standard compounds."*
**Conditions:** RG-ARP 20 mmol L⁻¹, initial pH 7.0, unbuffered, sealed pressure-resistant glass
bottle, oil bath, ice-water quench, HS-SPME-GC-MS, n = 3. **Units μg L⁻¹ throughout.**

| compound | RI | ID | **120 °C / 30 min** | **120 °C / 60 min** | **120 °C / 90 min** | **120 °C / 120 min** | **100 °C / 120 min** |
|---|---|---|---|---|---|---|---|
| **Furans / 6** | | | | | | | |
| furan | 798 | RI, MS, S | – | – | – | **0.43 ± 0.07ᵃ** | – |
| **furfural** | 1493 | RI, MS, S | 27.20 ± 0.59ᶜ | 52.96 ± 3.07ᵇ | 70.99 ± 3.93ᵃ | **74.25 ± 1.14ᵃ** | **9.28 ± 1.16ᵈ** |
| 2-methylfuran | 839 | RI, MS, S | 0.66 ± 0.06ᵈ | 1.00 ± 0.14ᶜ | 1.41 ± 0.08ᵇ | **2.00 ± 0.15ᵃ** | **0.64 ± 0.03ᵈ** |
| **4-hydroxy-5-methyl-3(2H)-furanone** | 2208 | RI, MS, S | 22.45 ± 0.69ᵈ | 31.32 ± 1.23ᶜ | 38.78 ± 1.22ᵇ | **41.69 ± 0.59ᵃ** | **8.82 ± 0.09ᵉ** |
| 1-(2-furanyl)-ethanone | 1534 | RI, MS, S | 0.03 ± 0.01ᵈ | 0.18 ± 0.03ᶜ | 0.50 ± 0.05ᵇ | **0.71 ± 0.05ᵃ** | – |
| 1-(2-furanyl)-1-propanone | 1569 | RI, MS, S | – | 0.16 ± 0.02ᶜ | 0.38 ± 0.03ᵇ | **0.66 ± 0.13ᵃ** | – |
| **subtotal** | | | **50.34** | **85.63** | **112.05** | **119.74** | **18.73** |
| **Sulfur-Containing Volatile Compounds / 7** | | | | | | | |
| 2-methylthiophene | 1093 | RI, MS, S | 0.09 ± 0.04ᵈ | 0.43 ± 0.03ᶜ | 1.00 ± 0.10ᵇ | **1.32 ± 0.13ᵃ** | **1.27 ± 0.14ᵃ** |
| 2-thiophenecarboxaldehyde | 1739 | RI, MS, S | – | 0.21 ± 0.01ᶜ | 0.28 ± 0.02ᵇ | **0.37 ± 0.01ᵃ** | – |
| thiazole | 1183 | RI, MS, S | – | 0.14 ± 0.01ᶜ | 0.21 ± 0.02ᵇ | **0.29 ± 0.01ᵃ** | – |
| 4-methyl-5-ethylthiazole | 1467 | RI, MS, S | – | 0.28 ± 0.04ᶜ | 0.41 ± 0.02ᵇ | **0.98 ± 0.05ᵃ** | **0.26 ± 0.02ᶜ** |
| ★ **2-methyl-3-furanthiol (MFT)** | 1305 | RI, MS, S | 4.45 ± 0.38ᶜ | 5.95 ± 0.04ᵇ | 6.62 ± 0.64ᵃ˒ᵇ | ★ **7.22 ± 0.44ᵃ** | ★ **4.39 ± 0.22ᶜ** |
| ★ **2-furfurylthiol (FFT)** | 1434 | RI, MS, S | 6.04 ± 0.33ᵈ | 8.27 ± 0.52ᶜ | 9.66 ± 0.52ᵇ | ★ **10.46 ± 0.41ᵃ** | ★ **5.05 ± 0.02ᵉ** |
| 3-(methylthio)-2-butanone | 1378 | RI, MS, S | 5.92 ± 0.20ᵇ | 9.31 ± 0.49ᵃ | 10.16 ± 0.26ᵃ | **10.10 ± 1.01ᵃ** | **4.25 ± 0.30ᶜ** |
| **subtotal** | | | **16.51** | **24.60** | **28.34** | **30.73** | **15.21** |
| **total compounds** | | | **66.85** | **110.23** | **140.39** | **150.47** | **33.94** |

**⚠️ Note the missing cells:** at **100 °C only the 120-min column exists in this table**. Four
of the seven sulfur compounds and three of the six furans are **"not detected" at 100 °C/120 min**
— so seven of the thirteen rows cannot yield a fold-change at all. `[NEG]`
**⚠️ Note the RI drift vs Kang 2026** `[!]`: Feng gives MFT **RI 1305** and FFT **RI 1434**;
Kang's Table S4 gives MFT 1302 and FFT 1426 on the same DB-WAX phase. Within-tolerance, and both
carry `S` / external-standard confirmation, so identity is not in doubt. Feng's furfural RI is
**1493** against Kang's **1457** — a 36-unit gap that is larger than comfortable but again both
are standard-confirmed. `[!]` Noted, not fatal.

### 3b. **TABLE 2 — the ARP + GSH ladder** `[M]`

**Caption, raster-verified verbatim:** *"Table 2. Volatile Compounds (μg/L) Identified in the
RG-ARP and GSH Model Reactions at Different Temperatures (100 and 120 °C)ᵃ"* — p. 9101.
Footnotes identical to Table 1. **Conditions as Table 1, plus GSH 20 mmol L⁻¹ (equimolar).**

| compound | RI | ID | **120 °C / 30 min** | **120 °C / 60 min** | **120 °C / 90 min** | **120 °C / 120 min** | **100 °C / 120 min** |
|---|---|---|---|---|---|---|---|
| **Furans / 7** | | | | | | | |
| furan | 798 | RI, MS, S | – | – | **0.14 ± 0.03ᵃ** | – | – |
| **furfural** | 1493 | RI, MS, S | 26.37 ± 1.99ᶜ | 37.80 ± 2.15ᵇ | 39.90 ± 1.60ᵇ | **44.06 ± 1.98ᵃ** | **6.97 ± 0.75ᵈ** |
| 2-methylfuran | 839 | RI, MS, S | 0.70 ± 0.06ᵈ | 1.50 ± 0.17ᶜ | 2.11 ± 0.12ᵇ | **4.61 ± 0.42ᵃ** | **1.24 ± 0.32ᶜ** |
| **4-hydroxy-5-methyl-3(2H)-furanone** | 2208 | RI, MS, S | 9.21 ± 0.21ᵈ | 11.47 ± 0.42ᶜ | 18.25 ± 0.55ᵇ | **21.89 ± 1.19ᵃ** | **6.95 ± 0.25ᵉ** |
| 1-(2-furanyl)-ethanone | 1534 | RI, MS, S | – | 0.10 ± 0.02ᶜ | 0.17 ± 0.03ᵇ | **0.32 ± 0.03ᵃ** | – |
| 1-(2-furanyl)-1-propanone | 1569 | RI, MS, S | – | 0.11 ± 0.02ᶜ | 0.30 ± 0.03ᵇ | **0.46 ± 0.01ᵃ** | – |
| 2(5H)-furanone | 1704 | RI, MS, S | – | – | – | **2.69 ± 0.15ᵃ** | – |
| **subtotal** | | | **36.28** | **50.98** | **60.88** | **74.04** | **15.16** |
| **Sulfur-Containing Volatile Compounds / 9** | | | | | | | |
| thiophene | 1045 | RI, MS, S | – | – | – | **0.04 ± 0.01ᵃ** | – |
| 2-methylthiophene | 1093 | RI, MS, S | 0.31 ± 0.03ᵈ | 1.66 ± 0.12ᶜ | 2.20 ± 0.36ᵇ | **3.60 ± 0.38ᵃ** | **1.39 ± 0.03ᶜ** |
| 2-thiophenecarboxaldehyde | 1739 | RI, MS, S | 0.12 ± 0.02ᵈ | 0.45 ± 0.02ᶜ | 0.84 ± 0.09ᵇ | **0.98 ± 0.07ᵃ** | **0.05 ± 0.01ᵈ** |
| thiazole | 1183 | RI, MS, S | – | 0.18 ± 0.01ᶜ | 0.26 ± 0.03ᵇ | **0.37 ± 0.03ᵃ** | – |
| 4-methyl-5-ethylthiazole | 1467 | RI, MS, S | – | 0.53 ± 0.04ᶜ | 0.90 ± 0.03ᵇ | **1.30 ± 0.26ᵃ** | **0.29 ± 0.03ᵈ** |
| 2-methyl-2-thiazoline | 2125 | RI, MS, S | – | – | – | **0.72 ± 0.27ᵃ** | – |
| ★ **2-methyl-3-furanthiol (MFT)** | 1305 | RI, MS, S | 11.27 ± 0.28ᵈ | 21.88 ± 1.01ᶜ | 28.63 ± 0.76ᵇ | ★ **32.15 ± 0.49ᵃ** | ★ **9.51 ± 0.18ᵉ** |
| ★ **2-furfurylthiol (FFT)** | 1434 | RI, MS, S | 14.31 ± 0.41ᵈ | 30.45 ± 0.36ᶜ | 39.35 ± 0.13ᵇ | ★ **44.63 ± 0.34ᵃ** | ★ **10.92 ± 2.04ᵉ** |
| 3-(methylthio)-2-butanone | 1378 | RI, MS, S | 11.69 ± 0.23ᵈ | 27.89 ± 0.75ᶜ | 30.04 ± 0.69ᵇ | **35.53 ± 0.31ᵃ** | **7.88 ± 0.13ᵉ** |
| **subtotal** | | | **37.70** | **83.04** | **102.22** | **119.32** | **30.04** |
| **Pyrazines / 1** | | | | | | | |
| methylpyrazine | 1194 | RI, MS, S | – | – | – | **0.04 ± 0.01ᵃ** | – |
| **subtotal** | | | **0.00** | **0.00** | **0.00** | **0.04** | **0.00** |
| **total compounds** | | | **73.98** | **134.02** | **165.10** | **193.40** | **45.20** |

**⚠️ [!] The `furan` row is non-monotone in time:** detected at 90 min (0.14) but **not** at 30,
60 or **120** min. Every other row in both tables rises monotonically. Either furan is being
consumed faster than it forms above 90 min, or that cell is an artefact. It contributes 0.19 %
of the 90-min furan subtotal, so nothing downstream is affected — but flag it.

### 3c. **INTERNAL-CONSISTENCY AUDIT — Tables 1 and 2 pass completely** `[D]`

Every one of the 18 printed subtotals and totals recomputed from its own member rows:

| aggregate | recomputed | printed | Δ |
|---|---|---|---|
| T1 furan subtotal, 120 °C, 30/60/90/120 min | 50.34 / 85.62 / 112.06 / 119.74 | 50.34 / 85.63 / 112.05 / 119.74 | ≤ 0.01 |
| T1 furan subtotal, 100 °C | **18.74** | 18.73 | 0.01 (rounding) |
| T1 sulfur subtotal, 120 °C, 30/60/90/120 min | 16.50 / 24.59 / 28.34 / 30.74 | 16.51 / 24.60 / 28.34 / 30.73 | ≤ 0.01 |
| T1 sulfur subtotal, 100 °C | **15.22** | 15.21 | 0.01 |
| T1 total, 120 °C / 100 °C | 150.47 / 33.94 | 150.47 / 33.94 | **0.000** |
| T2 furan subtotal, 120 °C, 30/60/90/120 min | 36.28 / 50.98 / 60.87 / 74.03 | 36.28 / 50.98 / 60.88 / 74.04 | ≤ 0.01 |
| T2 furan subtotal, 100 °C | **15.16** | 15.16 | **0.000** |
| T2 sulfur subtotal, 120 °C, 30/60/90/120 min | 37.70 / 83.04 / 102.22 / 119.32 | same | **0.000** |
| T2 sulfur subtotal, 100 °C | **30.04** | 30.04 | **0.000** |
| T2 total, 120 °C / 100 °C | 193.40 / 45.20 | 193.40 / 45.20 | **0.000** |

**⇒ Every aggregate closes to ≤ 0.01 μg L⁻¹ (pure decimal rounding). Tables 1 and 2 are
arithmetically sound in all 18 aggregates.** Combined with the raster verification of §0b, and
with the independent external anchor of §4b (Fig. 2g reproduces both sulfur subtotals to the
digitisation limit), **Tables 1 and 2 are high-confidence artefacts.**
**⇒ Note the contrast with `kang2026_SI_extraction.md` §4c.** Kang's Tables S4/S5 pass their
arithmetic audit but fail an SD-consistency audit (identical means, different SDs). **Feng has
no duplicated condition anywhere, so that particular test cannot be run** — the SDs here are
neither confirmed nor impeached. Treat them as reported, not as validated.

### 3d. **TABLE 3 — pyrazines across all four extra-added arms, both temperatures** `[M]`

**⚠️ Table 3 is typeset ROTATED 90° on p. 9102 and its text layer is unusable** — `pdftotext`
interleaves the caption, the eight column headers and the three data rows with unrelated body
text from the facing column, producing fragments that cannot be reassembled without ambiguity.
**This transcription is taken entirely from a de-rotated 600 dpi raster.**

**Caption, raster-verified verbatim:** *"Table 3. Pyrazines (μg/L) Identified in the RG-ARP and
GSH or Its Constituent Amino Acid Model Reactions at Different Temperatures (100 and 120 °C)ᵃ"*
**Column super-header, raster-verified: `120 min`** — ★ **Table 3 is a single-time-point table.
There is no time course for pyrazines at either temperature.**
**Footnotes:** identical wording to Tables 1 and 2.
**Conditions:** RG-ARP 20 mmol L⁻¹ + the named co-substrate at 20 mmol L⁻¹, pH 7.0 initial,
120 min, n = 3.

| compound | RI | ID | **ARP+Glu 100 °C** | **ARP+Cys 100 °C** | **ARP+Gly 100 °C** | **ARP+GSH 100 °C** | **ARP+Glu 120 °C** | **ARP+Cys 120 °C** | **ARP+Gly 120 °C** | **ARP+GSH 120 °C** |
|---|---|---|---|---|---|---|---|---|---|---|
| **pyrazine** | 1160 | RI, MS, S | – | **0.43 ± 0.02ᵇ** | – | – | – | **0.53 ± 0.17ᵃ** | – | – |
| **methylpyrazine** | 1194 | RI, MS, S | **0.06 ± 0.00ᶜ** | **0.08 ± 0.01ᵇ** | **0.04 ± 0.01ᵈ** | – | **0.08 ± 0.01ᵇ** | **0.12 ± 0.01ᵃ** | **0.06 ± 0.00ᶜ** | **0.04 ± 0.01ᵈ** |
| **total compounds** | | | **0.06** | **0.51** | **0.04** | **0.00** | **0.08** | **0.65** | **0.06** | **0.04** |

**Arithmetic audit `[D]`:** 0.43 + 0.08 = **0.51** ✔ · 0.53 + 0.12 = **0.65** ✔ · all six
single-compound totals ✔. **Table 3 closes exactly.**

**Structural readings `[D]`:**
1. ★ **Pyrazine (the unsubstituted parent) is formed ONLY in the ARP + Cys arm**, at both
   temperatures, and by nobody else. **Cys is the only co-substrate that makes it.** The paper's
   own explanation: *"Cys was more likely to react with GO to specifically generate pyrazine
   that was not detected in other reaction models"* (p. 9102).
2. **Rank order of total pyrazines is identical at both temperatures:
   Cys ≫ Glu > Gly > GSH.** At 100 °C: 0.51 / 0.06 / 0.04 / 0.00. At 120 °C: 0.65 / 0.08 /
   0.06 / 0.04. **The selectivity ordering is temperature-invariant.**
3. **The extra-added free amino acids all beat the intact peptide.** GSH gives **zero** pyrazine
   at 100 °C and 0.04 μg L⁻¹ at 120 °C — the paper's headline mechanism, that GSH must first be
   hydrolysed to Cys-Gly then to Cys before it can act as a Strecker nucleophile, and that
   hydrolysis is slow.
4. Author's own arithmetic, verified: *"The concentration of the pyrazines formed by the
   ARP + Cys-120 °C model was 16.25 times greater than that formed by the ARP + GSH-120 °C
   model"* — **0.65 / 0.04 = 16.25** ✔ exactly.
5. **⚠️ 100 → 120 °C fold-changes for pyrazines are tiny and near the reported SD:**
   ARP+Cys **× 1.27**, ARP+Glu **× 1.33**, ARP+Gly **× 1.50**, ARP+GSH **undefined (0 → 0.04)**.
   Apparent Eₐ 14.6 / 17.4 / 24.7 kJ mol⁻¹ `[D]`. **All three are inside 2 SD of the printed
   uncertainties (pyrazine at 120 °C is 0.53 ± 0.17, a 32 % relative SD) and none should be
   treated as a resolved temperature response.**
6. ★ **[!] THERE IS NO ARP-ALONE COLUMN IN TABLE 3.** The paper's central claim — that adding a
   nucleophile *promotes* pyrazine formation — has **no no-addition control at either
   temperature**. Table 1 (ARP alone) contains no pyrazine row at all, so we cannot tell whether
   ARP alone gives 0.00 or 0.05 μg L⁻¹. **The "promotion" is unquantified against a baseline.**
   `[NEG]` — this is the single largest logical gap in the paper.

### 3e. **Numbers stated only in the body text `[M]`** — several of them "(data not shown)"

| quantity | value | arm / condition | anchor | status |
|---|---|---|---|---|
| ARP degradation | **93.50 %** total decrease at 120 min | ARP alone, 120 °C | p. 9100 | `[M]` |
| ARP concentration | **15.86 → 9.86 mmol L⁻¹** | ARP + GSH, 100 °C, 30 → 120 min | p. 9101 | `[M]` |
| ARP concentration | **2.97 → 0.83 mmol L⁻¹** | ARP + GSH, 120 °C, 30 → 120 min | p. 9102 | `[M]` |
| PCA released | **26.83 mmol L⁻¹** at 120 min | ARP + GSH, 120 °C | p. 9102 | `[M]` |
| Cys-Gly | **4.01 mmol L⁻¹** at 120 min | ARP + GSH, 120 °C | p. 9102 | `[M]` |
| residual GSH | **11.84 mmol L⁻¹** | ARP + Cys, 100 °C, 120 min | p. 9102 | `[M]` |
| residual Cys | **12.68 mmol L⁻¹** | ARP + Cys, 100 °C, 120 min | p. 9102 | `[M]` |
| **reacted Cys** | **9.15 mmol L⁻¹** | ARP + Cys, 100 °C, 120 min | p. 9102 | ⚠️ **"(data not shown)"** |
| **reacted GSH** | **11.79 mmol L⁻¹** | ARP + GSH, 100 °C, 120 min | p. 9102 | `[M]` in text |
| **reacted Glu** | **2.27 mmol L⁻¹** | ARP + Glu, 100 °C, 120 min | p. 9102 | ⚠️ **"(data not shown)"** |
| **reacted Gly** | **2.17 mmol L⁻¹** | ARP + Gly, 100 °C, 120 min | p. 9102 | ⚠️ **"(data not shown)"** |
| **reacted GSH / reacted Cys** | **33.32 / 25.84 mmol L⁻¹** | 120 °C, 120 min | p. 9103 | ⚠️ **"(data not shown)"** |
| **3-(methylthio)-2-butanone** | ★ **36.41 μg L⁻¹** | **ARP + Cys, 100 °C, 120 min** | p. 9103 | ⚠️ **"(data not shown)"** |
| 3-DO / 1-DO | **0.003 / 0.012 mmol L⁻¹** | ARP + Cys, 100 °C, 120 min | p. 9102 | `[M]`, = Fig 3e |
| 3-DO / 1-DO / GO / MGO | **0.038 / 0.109 / 0.057 / 0.009 mmol L⁻¹** | ARP + GSH, 100 °C, 120 min | p. 9102 | `[M]`, = Fig 3b |
| GO / MGO | **0.103 / 0.010 mmol L⁻¹** | ARP + Glu, 100 °C, 120 min | p. 9102 | `[M]`, = Fig 3e |
| GO | **0.215 / 0.165 mmol L⁻¹** | ARP + Glu / ARP + Gly, 120 °C | p. 9102 | `[M]`, = Fig 3f |
| thiazoles total | **2.39 μg L⁻¹** | ARP + GSH, 120 °C, 120 min | p. 9103 | `[M]`, = T2 (0.37+1.30+0.72) ✔ |
| thiazole ratio | **× 2.48** ARP+Cys-120 ÷ ARP+GSH-120 | 120 min | p. 9103 | `[M]`, Fig S3 (SI, not on disk) |
| Rib + GSH total volatiles | **23.06 → 70.84 μg L⁻¹**, 30 → 120 min, 120 °C | comparator arm | p. 9100/9103 | `[M]`, Fig S2 (SI) |
| furfural share of furans | **49.55 % (100 °C) / 62.01 % (120 °C)** at 120 min | ARP alone | p. 9100 | `[F]`-style; **audit `[D]`: 9.28/18.73 = 49.55 % ✔ · 74.25/119.74 = 62.01 % ✔ exactly** |
| MFT + 3-(methylthio)-2-butanone | **× 3.89** ARP+GSH-120 ÷ ARP+GSH-100, 120 min | | p. 9103 | **audit `[D]`: (32.15+35.53)/(9.51+7.88) = 67.68/17.39 = 3.892 ✔** |
| sulfur : furans in ARP+GSH-120 | **× 1.98** | 120 min | p. 9103 | **audit `[D]`: 119.32/74.04 = 1.612 ✗ — the paper says 1.98** `[!]` |

**⚠️ [!] One text-vs-table arithmetic failure.** p. 9103 states *"the concentration of
sulfur-containing volatile compounds was 1.98 times higher than that of furans in the ARP + GSH
model for 120 min."* From Table 2 the ratio is **119.32 / 74.04 = 1.612**, not 1.98. **No
combination of Table 2's printed cells gives 1.98.** Every other numeric claim in the text
audits correctly (six checked above, all exact), so this is most likely a stale number from an
earlier draft. **Use 1.612 from the table, not 1.98 from the text.** ★ The paper's *qualitative*
claim — that adding GSH inverts the furan : sulfur balance — survives: ARP alone gives
sulfur/furan = 30.73/119.74 = **0.257**; ARP + GSH gives **1.612**. **A 6.3-fold swing in
branching.**

**⚠️ [!] `2-methyl-2-thiazoline` is printed with RI 2125 in Table 2.** Its literature DB-WAX RI
is ≈ 1310–1360. An RI of 2125 on a wax column is far too late for a C₄H₇NS thiazoline. Either
the RI is a typo or the identification is wrong despite the `S` flag. It contributes 0.60 % of
the ARP+GSH-120 sulfur subtotal. **Exclude 2-methyl-2-thiazoline from any structural claim; it
does not affect any subtotal materially.**

---

## §4. EVERY FIGURE — what was digitised, how, and with what error

**Digitisation method, stated once.** Pages 9099 (Fig 2) and 9100 (Fig 3) were rendered with
`pdftoppm -r 600 -png` and read as images; the axis tick labels supply the pixel-to-value
calibration linearly between the printed gridline values (e.g. Fig 2a: 0 and 20 mmol L⁻¹ on the
left axis). **No automated tracing was used** — readings are visual against the rendered tick
grid. **Every digitised value is marked `[D]`.**

**⭐ The digitisation is externally validated.** Fig 2g's four 120-min endpoints and Fig 3b's
four 120-min endpoints are independently pinned by Tables 1/2 and by the body text
respectively. My readings reproduce all eight to within **± 0.5 μg L⁻¹** and **± 0.005 mmol L⁻¹**.
**Those are the honest error bars for every other point on those panels.**

### 4a. **Figure 1 — schemes only, NOT digitisable** `[NEG]`

**Caption verbatim:** *"Figure 1. Total-ion chromatogram (a), UPLC−MS/MS spectrum (b), and the
proposed fragmentation mechanism of N-(1-deoxy-D-ribulos-1-yl)-glutathione (c); the possible
formation pathway of furans, furanthiols, and thiophenes (d), and the possible formation pathway
of pyrazines and thiazoles (e)."* — p. 9096.
Panels (a) and (b) are a TIC and an MS/MS spectrum with **no legible axis values and no peak
table**; the m/z assignments are given in the body text instead (§1a) and are the only usable
content. Panels (c), (d), (e) are **hand-drawn mechanism schemes with no quantities**.
**⇒ Nothing in Figure 1 is digitisable, and nothing numeric is lost by not digitising it.**

### 4b. **Figure 2 — the precursor and product time courses. DIGITISED.**

**Caption verbatim:** *"Figure 2. Time-course levels of ARP (a), GSH (b), PCA (c), Cys-Gly (d),
Cys (e), Gly (f), and sulfur-containing volatile compounds (g) in the ARP and ARP + GSH models;
the concentrations of GSH, PCA, Cys-Gly, Cys, and Gly (h) in the GSH models for 120 min at an
initial pH value of 7.0 at different temperatures (100 and 120 °C)."* — p. 9099.
**Series legend, identical on panels a–g:** ▲ ARP-120 °C · ■ ARP-100 °C · ● ARP+GSH-120 °C ·
— (no marker) ARP+GSH-100 °C. **Panel h is a bar pair:** ▫ GSH-100 °C, ▪ GSH-120 °C.
Error bars are drawn on most points and are small (≲ 3 % of value).

#### Panel (a) — **ARP concentration, mmol L⁻¹** `[D]`, ± 0.4

| min | **ARP-100 °C** | **ARP-120 °C** | **ARP+GSH-100 °C** | **ARP+GSH-120 °C** |
|---|---|---|---|---|
| 30 | **17.8** | **4.55** | **15.85** *(text: 15.86 ✔)* | **2.95** *(text: 2.97 ✔)* |
| 60 | **16.7** | **2.90** | 14.6 | 1.45 |
| 90 | **14.05** | **1.70** | 13.1 | 1.05 |
| 120 | **9.9** | **1.28** | 10.9 | **0.85** *(text: 0.83 ✔)* |

**⚠️ [!] The two 100 °C curves cross between 90 and 120 min, and the ARP-100 vs ARP+GSH-100
assignment at 120 min cannot be resolved from the raster.** The text pins ARP+GSH-100 as
**15.86 → 9.86 mmol L⁻¹**, and 15.86 matches the *lower* 30-min point while 9.86 matches the
*lower* 120-min point — but the two curves visibly swap order in between. Both endpoints are
text-anchored; only the middle assignment is ambiguous. **For §6c I use ARP-alone at
100 °C / 120 min = 10.4 ± 0.8 mmol L⁻¹, the midpoint of the two candidate readings, and carry
that uncertainty through.**

**⭐ Panel (a) is the paper's most valuable single artefact for the repo**, because it converts
a yield ladder into a yield-per-precursor-consumed ladder. See §6c.
**⚠️ Note there is no t = 0 point on any panel; the earliest reading is 30 min.** At 120 °C the
ARP has already fallen from a nominal 20 to **4.55** mmol L⁻¹ by the first sample — **77 % of
the substrate is gone before the first observation.** Any kinetic reading of the 120 °C series
is therefore reading the tail of a reaction whose fast phase was never observed. `[!]`

#### Panels (b)–(f) — GSH, PCA, Cys-Gly, Cys, Gly, mmol L⁻¹ `[D]`, ± 0.3 (b, c, f), ± 0.05 (d), ± 0.02 (e)

Read at **120 min** only (the shape is monotone in every panel and the endpoint carries the
information):

| panel | analyte | **ARP-100 °C** | **ARP-120 °C** | **ARP+GSH-100 °C** | **ARP+GSH-120 °C** |
|---|---|---|---|---|---|
| b | **GSH** | ~0.6 | ~0.4 | **18.3** | **5.8** |
| c | **PCA** | 4.0 | 8.6 | 6.2 | **26.8** *(text: 26.83 ✔)* |
| d | **Cys-Gly** | ~0.06 | ~0.10 | 0.62 | **4.0** *(text: 4.01 ✔)* |
| e | **Cys** | ~0.02 | ~0.06 | **0.64** | **0.79** |
| f | **Gly** | 2.85 | 4.6 | 4.7 | **17.9** |

**⭐ Structural readings `[D]`:**
1. **Cys is the scarcest species in the system by two orders of magnitude.** At the most
   productive condition (ARP + GSH, 120 °C, 120 min) Gly reaches **17.9 mmol L⁻¹** while Cys
   reaches only **0.79** — a **22.6 : 1 imbalance** on two amino acids that Cys-Gly hydrolysis
   must release **equimolar**. **⇒ ≈ 96 % of the Cys released is consumed as fast as it appears.**
   This is the paper's central mechanistic claim and the figure supports it quantitatively.
2. **PCA ≫ Glu.** PCA reaches 26.8 mmol L⁻¹ while free Glu never appears as its own panel — the
   glutamyl residue cyclises to pyroglutamate rather than hydrolysing. Explains why extra-added
   *Glu* outperforms extra-added *GSH* for pyrazines (Table 3).
3. **GSH survives far better at 100 °C than at 120 °C**: 18.3 vs 5.8 mmol L⁻¹ residual at
   120 min in the ARP + GSH arm — a **× 3.16 difference in survival**, apparent Eₐ of GSH
   consumption **≈ 30 kJ mol⁻¹** `[D]` on a crude (20−C)/… basis (see §6d).

#### Panel (g) — ★ **sulfur-containing volatiles, μg L⁻¹ — THE 100 °C TIME COURSES TABLES 1 AND 2 OMIT** `[D]`, ± 0.8

| min | **ARP-100 °C** ■ | **ARP-120 °C** ▲ | **ARP+GSH-100 °C** — | **ARP+GSH-120 °C** ● |
|---|---|---|---|---|
| 30 | ★ **4.6** `[D]` | 16.5 *(T1: 16.51 ✔)* | ★ **15.3** `[D]` | 37.8 *(T2: 37.70 ✔)* |
| 60 | ★ **9.0** `[D]` | 24.3 *(T1: 24.60 ✔)* | ★ **23.8** `[D]` | 82.8 *(T2: 83.04 ✔)* |
| 90 | ★ **11.4** `[D]` | 28.0 *(T1: 28.34 ✔)* | ★ **27.6** `[D]` | 102.2 *(T2: 102.22 ✔)* |
| 120 | 15.1 *(T1: 15.21 ✔)* | 30.8 *(T1: 30.73 ✔)* | 30.0 *(T2: 30.04 ✔)* | 119.3 *(T2: 119.32 ✔)* |

**⭐ This is the most useful `[D]` product of the whole extraction. Eight of the sixteen cells
are already in Tables 1 and 2 and my digitisation reproduces all eight to ≤ 0.4 μg L⁻¹ — so the
six starred cells are new, calibrated data with a ± 0.8 μg L⁻¹ bar.** They give the repo the
**100 °C sulfur-class time course, which the tables do not print**, and with it the ability to
compare *initial rates* rather than only 120-min yields:

| slope over 30 → 60 min, μg L⁻¹ min⁻¹ `[D]` | value | ratio 120/100 |
|---|---|---|
| ARP alone, 100 °C | **0.147** | |
| ARP alone, 120 °C | **0.260** | **× 1.77** |
| ARP + GSH, 100 °C | **0.283** | |
| ARP + GSH, 120 °C | **1.500** | **× 5.30** |

**⇒ Early-rate ratios (× 1.77 ARP alone) are LOWER than the 120-min yield ratios (× 2.02), and
in the ARP + GSH arm they are HIGHER (× 5.30 vs × 3.97) — because the 120 °C ARP + GSH run is
saturating hard between 90 and 120 min while the 100 °C run is not.** Apparent Eₐ from the
30–60 min slopes: **ARP alone 34.8 kJ mol⁻¹, ARP + GSH 101.7 kJ mol⁻¹** `[D]`. Compare the
yield-based 42.9 and 84.1. **The two definitions disagree by 8 and 18 kJ mol⁻¹ respectively.
That spread is an honest floor on the systematic uncertainty of *any* Eₐ taken off this ladder.**

#### Panel (h) — GSH-alone comparator, 120 min, bar chart, mmol L⁻¹ `[D]`, ± 0.15

| species | **GSH-100 °C** | **GSH-120 °C** | significance as printed |
|---|---|---|---|
| **GSH** (residual) | **13.1** | **4.7** | ** |
| **PCA** | **4.5** | **10.6** | ** |
| **Cys-Gly** | **1.85** | **2.85** | * |
| **Cys** | **1.80** | **5.35** | ** |
| **Gly** | **3.0** | **9.15** | ** |

**⇒ In the absence of the ARP, Cys survives to 5.35 mmol L⁻¹ at 120 °C — 6.8 × its level in the
ARP + GSH arm (0.79).** Direct confirmation that the ARP's carbonyl pool is what consumes the
Cys. `[D]` **Apparent Eₐ of GSH hydrolysis to Cys, 100→120 °C: `ln(5.35/1.80) × 60.98` =
66.4 kJ mol⁻¹** `[D]` — a genuinely useful number for the sulfur lane's H₂S-source term, and the
only clean amino-acid-release Eₐ in the paper. **⚠️ It is a yield-at-120-min Eₐ, not a rate Eₐ,
and the GSH-alone system is unbuffered like everything else.**
**⚠️ [!] Mass-balance check `[D]`:** at 120 °C, PCA (10.6) + Cys-Gly (2.85) + residual GSH (4.7)
= 18.15 against a 20 mmol L⁻¹ nominal load — closes to 91 %. At 100 °C: 4.5 + 1.85 + 13.1 =
19.45 = 97 %. **Both close acceptably.** But Cys (5.35) and Gly (9.15) at 120 °C are **not
equimolar** (they should be, from Cys-Gly hydrolysis) — a 1.71 : 1 excess of Gly, consistent
with Cys being lost to H₂S/volatiles even without an ARP present.

### 4c. **Figure 3 — the α-dicarbonyl time courses and the amino-acid comparison. DIGITISED.**

**Caption verbatim:** *"Figure 3. Time-course levels of 1-DO, 3-DO, GO, and MGO in the
ARP-100 °C (a), ARP + GSH-100 °C (b), ARP-120 °C (c), and ARP + GSH-120 °C (d) models; the
concentrations of 1-DO, 3-DO, GO, and MGO in the thermal reaction of RG-ARP and Glu, Cys, or Gly
for 120 min at 100 (e) and 120 °C (f) at an initial pH value of 7.0."* — p. 9100.

**⚠️ [!] LEGEND / CAPTION NOMENCLATURE MISMATCH — RESOLVED, NOT AN ERROR.** The plot legends
read **`3-DP`** and **`1-DP`** while the caption and the entire body text read **`3-DO`** and
**`1-DO`**. `DP` = **deoxy*pentos*one**, which is the chemically correct name here (the sugar is
ribose, a pentose); `DO` = deoxyosone, the generic. **They are the same two species.** I
initially suspected the ▲/■ series were swapped relative to the text, and **checked it at
600 dpi: the UPPER curve carries ■ = 1-DP and the LOWER carries ▲ = 3-DP, which is exactly what
the text asserts** (*"a lower concentration of 3-DO compared to that of 1-DO was observed in the
ARP-only models"*, p. 9100). **No swap. The figure and text agree.**

**Series:** ■ solid = **1-DP (1-DO)** · ▲ solid = **3-DP (3-DO)** · ● long-dashed = **GO** ·
fine-dotted, no marker = **MGO**. Y-axis **Concentration (mmol/L)**, 0–0.5 on panels a–d.

#### Panels (a)–(d), all four series, mmol L⁻¹ `[D]`, ± 0.008 on (a)/(b), ± 0.012 on (c)/(d)

| panel / min | **1-DO** | **3-DO** | **GO** | **MGO** |
|---|---|---|---|---|
| **(a) ARP alone, 100 °C** | | | | |
| 30 | **0.400** | 0.161 | 0.048 | ~0.003 |
| 60 | 0.365 | 0.121 | 0.060 | ~0.005 |
| 90 | 0.283 | 0.095 | 0.065 | ~0.007 |
| 120 | **0.212** | **0.072** | **0.072** | **~0.010** |
| **(b) ARP + GSH, 100 °C** | | | | |
| 30 | **0.285** | 0.135 | 0.041 | ~0.002 |
| 60 | 0.237 | 0.093 | 0.045 | ~0.004 |
| 90 | 0.170 | 0.058 | 0.052 | ~0.007 |
| 120 | **0.109** ✔*(text 0.109)* | **0.038** ✔*(text 0.038)* | **0.057** ✔*(text 0.057)* | **0.009** ✔*(text 0.009)* |
| **(c) ARP alone, 120 °C** | | | | |
| 30 | **0.248** | 0.121 | 0.076 | ~0.005 |
| 60 | 0.083 | 0.082 | 0.090 | ~0.005 |
| 90 | 0.038 | 0.030 | 0.115 | ~0.006 |
| 120 | **0.022** | **0.021** | **0.145** | **~0.007** |
| **(d) ARP + GSH, 120 °C** | | | | |
| 30 | **0.119** | ~0.045 | 0.043 | ~0.005 |
| 60 | 0.023 | 0.020 | 0.032 | ~0.005 |
| 90 | 0.017 | 0.014 | 0.045 | ~0.005 |
| 120 | **0.014** | **0.012** | **0.058** | **~0.005** |

**⭐ Panel (b)'s four 120-min values are all four independently pinned by the body text and my
readings match all four exactly to the printed 3 d.p. — this is the calibration anchor that
sets the ± 0.008 bar on the whole of panels (a) and (b).**

**⭐ Structural readings `[D]`:**
1. **1-DO exceeds 3-DO at every point in every panel** — ratio 2.5–2.9 at 100 °C, falling to
   ≈ 1.0–1.2 by 120 min at 120 °C. **The 1-DO : 3-DO branching is itself temperature-dependent
   and converges at high temperature.** The paper reads this as 3-DO being the less stable of the
   two at pH 7.
2. ★ **GO is the only α-dicarbonyl that RISES with time, and it rises hardest at 120 °C**
   (0.076 → 0.145 mmol L⁻¹ in the ARP-alone arm, **× 1.91 over 90 min**) while both deoxyosones
   collapse (1-DO × 0.089). **Retro-aldol cleavage of the deoxyosones to GO outruns GO's own
   consumption once the temperature is high enough.** In the ARP + GSH arm GO is *suppressed*
   (0.058 vs 0.145 at 120 °C, **× 0.40**) — the added nucleophile traps it. This is the paper's
   headline trapping mechanism, and it is quantified here.
3. **MGO is ≈ 0.005–0.010 mmol L⁻¹ everywhere and never varies by more than 2 ×.** With
   the digitisation bar at ± 0.008 that is **≈ 100 % relative uncertainty. Do not use any MGO
   value from panels (a)–(d).** `[!]` (The panel (e)/(f) bars are better resolved.)
4. **100 → 120 °C fold-changes at 120 min, ARP alone `[D]`:** 1-DO **× 0.104** (Eₐ **−138**
   kJ mol⁻¹), 3-DO **× 0.292** (Eₐ **−75**), **GO × 2.014 (Eₐ +42.7)**, MGO × 0.7 (unusable).
   **Negative apparent Eₐ for the deoxyosones is exactly right and exactly expected — they are
   *intermediates*, and their steady-state level falls when their sink speeds up faster than
   their source. ⇒ A pool concentration is not a rate. This is the cleanest demonstration in the
   paper of why yield-at-fixed-time Eₐ must never be read as a barrier.**

#### Panels (e) and (f) — the amino-acid comparison bars, 120 min, mmol L⁻¹ `[D]`, ± 0.006

Bars in legend order: ▨ 3-DP · ▫ 1-DP · ▪ GO · ▫(dotted) MGO.

| arm | **3-DO** | **1-DO** | **GO** | **MGO** |
|---|---|---|---|---|
| **(e) 100 °C** | | | | |
| ARP + Glu | 0.011ᵃ˒ᵇ | 0.028ᵇ | **0.103ᵃ** ✔*(text 0.103)* | 0.011ᵃ ✔*(text 0.010)* |
| ARP + Cys | **0.003ᵇ** ✔*(text 0.003)* | **0.012ᵇ** ✔*(text 0.012)* | 0.029ᶜ | 0.006 |
| ARP + Gly | 0.018ᵃ | 0.046ᵃ | **0.077ᵇ** | 0.009ᵃ˒ᵇ |
| **(f) 120 °C** | | | | |
| ARP + Glu | 0.011ᵇ | 0.013ᵃ | **0.215ᵃ** ✔*(text 0.215)* | 0.021ᵃ |
| ARP + Cys | 0.002ᶜ | 0.003ᵇ | 0.031ᶜ | 0.008ᶜ |
| ARP + Gly | 0.013ᵃ | 0.015ᵃ | **0.165ᵇ** ✔*(text 0.165)* | 0.017ᵇ |

**⭐ Five of the twelve `[D]` values here are text-pinned and all five match. Calibration confirmed.**

**⭐ The single sharpest structural result in the paper `[D]`:**

| arm | **GO at 100 °C** | **GO at 120 °C** | **fold** | **apparent Eₐ** |
|---|---|---|---|---|
| ARP + Glu | 0.103 | **0.215** | **× 2.09** | **44.8 kJ mol⁻¹** |
| ARP + Gly | 0.077 | **0.165** | **× 2.14** | **46.5 kJ mol⁻¹** |
| ★ **ARP + Cys** | **0.029** | **0.031** | ★ **× 1.07** | ★ **4.1 kJ mol⁻¹** |

**⇒ Glu and Gly let GO accumulate and its accumulation doubles with temperature. Cys clamps GO
at ≈ 0.03 mmol L⁻¹ and the clamp is TEMPERATURE-INVARIANT.** Cys's trapping of glyoxal is fast
enough at 100 °C that raising to 120 °C changes nothing — **the trap is not rate-limiting at
either temperature. This is a clean, model-testable statement about the Cys–GO channel** and it
does not depend on the semi-quant deoxyosone calibration (GO has a printed external curve,
R² = 0.9998, §2e). **RATIO-ONLY on absolute magnitude (digitised), STRUCTURAL on the invariance.**

---

## §5. PUBLISHED KINETICS `[F]` — **THERE ARE NONE. `[NEG]`**

**⚠️ Feng 2022 contains no rate constant, no reaction order, no activation energy, no half-life,
no rate law, no kinetic model, and no fitted temperature dependence of any kind.** The word
"rate" appears only qualitatively (*"the thermal degradation rate of the RG-ARP was
accelerated"*). There is nothing here to audit and nothing to refit.

**The only quantitative temperature statements the authors make are bare fold-changes:**

| author claim | printed value | **my audit `[D]`** | verdict |
|---|---|---|---|
| *"The concentration of furans formed at 120 °C was more than 6.39 times that at 100 °C"* (Abstract) | **× 6.39** | 119.74 / 18.73 = **× 6.393** | ✔ **exact** |
| *"the total decrease was 93.50 % at 120 min"* (ARP, 120 °C) | **93.50 %** | implies ARP = 1.30 mmol L⁻¹ from 20; Fig 2a reads **1.28 ± 0.4** | ✔ |
| *"1.98 times higher than that of furans"* | **× 1.98** | 119.32 / 74.04 = **× 1.612** | ✗ **fails — see §3e** `[!]` |
| *"16.25 times greater"* (pyrazines, Cys vs GSH at 120 °C) | **× 16.25** | 0.65 / 0.04 = **× 16.25** | ✔ **exact** |
| *"3.89 times higher"* (MFT + methylthiobutanone, 120 vs 100 °C) | **× 3.89** | 67.68 / 17.39 = **× 3.892** | ✔ |
| *"2.48 times higher"* (thiazoles, Cys vs GSH at 120 °C) | **× 2.48** | Fig S3, **not on disk** | untestable `[NEG]` |
| *"49.55 and 62.01 % of the total furans"* (furfural share) | | 9.28/18.73 = **49.55 %**; 74.25/119.74 = **62.01 %** | ✔ **exact ×2** |

**⇒ Six of the seven auditable author claims reproduce exactly from the printed tables. The
seventh (× 1.98) does not and is almost certainly stale. Overall the paper's internal arithmetic
is strong.**

**⇒ Every activation energy anywhere in this dossier is `[D]` — derived by this wave, never
printed by the authors. Nothing here should ever look measured.**

---

## §6. DERIVED `[D]` — the two-point apparent activation energies

### 6a. **⚠️ THE CAVEAT CHAIN — read this before any number in §6b**

A two-point Eₐ from `Eₐ = R · ln(C₂/C₁) / (1/T₁ − 1/T₂)`. For T₁ = 373.15 K and
T₂ = 393.15 K, `1/T₁ − 1/T₂ = 1.36333 × 10⁻⁴ K⁻¹`, so

> **Eₐ [kJ mol⁻¹] = 60.98 × ln(fold-change)**

*(Arithmetic validated against the reference finding: Kang's published 120→140 °C legs of
97.8 kJ mol⁻¹ for MFT (× 4.26) and 69.2 for FFT (× 2.78) are reproduced by this formula as
**97.86** and **69.04**. The method is sound; only its inputs are in question.)*

1. **Two points determine a line with zero degrees of freedom.** There is no residual, so
   **R² is undefined and meaningless**, no goodness-of-fit statistic exists, and **no curvature
   can be detected — by construction.** A 2-point Eₐ is a re-expression of a ratio, not a fit.
2. ★ **A 2-point Eₐ can CORROBORATE a slope magnitude. It can NEVER test for a slope BREAK.**
   Feng's two temperatures both sit *inside* Kang's low leg. Nothing in this paper bears on
   whether the slope changes at 140 °C. **Any use of Feng as evidence for or against the
   switch-on is a category error.**
3. **Yield-at-fixed-time is not a rate constant.** Both columns are 120-min yields. If the
   product is at plateau at either temperature the ratio compresses; if the product is itself
   consumed the ratio can go *negative* (demonstrated directly for 1-DO and 3-DO in §4c).
4. ★ **Precursor depletion depresses apparent Eₐ at the top of the ladder — and Feng measures
   exactly how much.** 93.50 % of the ARP is consumed at 120 °C vs ≈ 48 % at 100 °C. **Every
   raw Eₐ below is therefore a LOWER BOUND on any true kinetic barrier**, and the correction is
   quantified in §6c.
5. **Semi-quant scale** — applies to 1-DO and 3-DO only (§2d); everything else is
   external-standard quantified but with **unavailable curves** (§2b).
6. **Unbuffered pH drift.** No buffer, initial pH 7.0, **final pH never measured at either
   temperature**. Maillard systems acidify, and the acidification is temperature-dependent, so
   the two columns are almost certainly at **different final pH**. The pH sensitivity is known
   to be large for exactly these analytes — `kang2026_SI_extraction.md` §5c measures MFT falling
   **× 0.57** and FFT **× 0.53** from pH 5.5 to 8.0. **⇒ A drift of even one pH unit between the
   two temperatures could account for a factor of ≈ 1.3 in the fold-change, i.e. ± 16 kJ mol⁻¹
   of the apparent Eₐ. This is an unbounded, undeclared systematic on every number in §6b.**
7. **Headspace, not liquid.** SPME at 60 °C — the partition coefficients of MFT and FFT are
   themselves temperature-dependent, but the *extraction* temperature is constant at 60 °C for
   all samples, so this cancels in the ratio. **This one is fine.**

### 6b. **Two-point apparent Eₐ, every compound with non-zero values at both temperatures**

**ARM ① — ARP ALONE (Table 1, 120 min).** *The arm most comparable to Kang 2026.*

| # | compound | **100 °C** | **120 °C** | **fold-change** | **apparent Eₐ, kJ mol⁻¹** `[D]` |
|---|---|---|---|---|---|
| 1 | **furfural** | 9.28 | 74.25 | **× 8.001** | **126.8** |
| 2 | **4-hydroxy-5-methyl-3(2H)-furanone** | 8.82 | 41.69 | **× 4.727** | **94.7** |
| 3 | 2-methylfuran | 0.64 | 2.00 | × 3.125 | **69.5** |
| | **FURAN SUBTOTAL** | 18.73 | 119.74 | **× 6.393** | **113.1** |
| 4 | 4-methyl-5-ethylthiazole | 0.26 | 0.98 | × 3.769 | **80.9** |
| 5 | 3-(methylthio)-2-butanone | 4.25 | 10.10 | × 2.376 | **52.8** |
| 6 | ★ **2-furfurylthiol (FFT)** | 5.05 | 10.46 | ★ **× 2.071** | ★ **44.4** |
| 7 | ★ **2-methyl-3-furanthiol (MFT)** | 4.39 | 7.22 | ★ **× 1.645** | ★ **30.3** |
| 8 | **2-methylthiophene** | 1.27 | 1.32 | ★ **× 1.039** | ★ **2.4** |
| | **SULFUR SUBTOTAL** | 15.21 | 30.73 | **× 2.020** | **42.9** |
| | **GRAND TOTAL** | 33.94 | 150.47 | × 4.433 | **90.8** |
| | *free thiols (MFT + FFT)* | 9.44 | 17.68 | **× 1.873** | **38.3** |
| | *thiazoles (all)* | 0.26 | 1.27 | × 4.885 | **96.7** |
| | *thiophenes (all)* | 1.27 | 1.69 | **× 1.331** | **17.4** |
| | not computable (absent at 100 °C) | \[furan, 1-(2-furanyl)-ethanone, 1-(2-furanyl)-1-propanone, 2-thiophenecarboxaldehyde, thiazole\] | | | `[NEG]` |

**ARM ② — ARP + GSH (Table 2, 120 min).**

| # | compound | **100 °C** | **120 °C** | **fold-change** | **apparent Eₐ, kJ mol⁻¹** `[D]` |
|---|---|---|---|---|---|
| 1 | ⚠️ **2-thiophenecarboxaldehyde** | **0.05** | 0.98 | × 19.600 | **181.5** ⚠️ *0.05 is at or below LOQ; ± 0.01 on the denominator swings Eₐ by ± 40. **Do not use.*** |
| 2 | **furfural** | 6.97 | 44.06 | × 6.321 | **112.5** |
| 3 | 4-methyl-5-ethylthiazole | 0.29 | 1.30 | × 4.483 | **91.5** |
| 4 | 3-(methylthio)-2-butanone | 7.88 | 35.53 | × 4.509 | **91.8** |
| 5 | ★ **2-furfurylthiol (FFT)** | 10.92 | 44.63 | ★ **× 4.087** | ★ **85.9** |
| 6 | 2-methylfuran | 1.24 | 4.61 | × 3.718 | **80.1** |
| 7 | ★ **2-methyl-3-furanthiol (MFT)** | 9.51 | 32.15 | ★ **× 3.381** | ★ **74.3** |
| 8 | **4-hydroxy-5-methyl-3(2H)-furanone** | 6.95 | 21.89 | × 3.150 | **70.0** |
| 9 | 2-methylthiophene | 1.39 | 3.60 | × 2.590 | **58.0** |
| | **FURAN SUBTOTAL** | 15.16 | 74.04 | **× 4.884** | **96.7** |
| | **SULFUR SUBTOTAL** | 30.04 | 119.32 | **× 3.972** | **84.1** |
| | **GRAND TOTAL** | 45.20 | 193.40 | × 4.279 | **88.7** |
| | *free thiols (MFT + FFT)* | 20.43 | 76.78 | **× 3.758** | **80.7** |
| | *thiazoles (all)* | 0.29 | 2.39 | × 8.241 | **128.6** |
| | *thiophenes (all)* | 1.44 | 4.62 | × 3.208 | **71.1** |
| | not computable (absent at 100 °C) | \[thiophene, thiazole, 2-methyl-2-thiazoline, 2(5H)-furanone, 1-(2-furanyl)-ethanone, 1-(2-furanyl)-1-propanone, methylpyrazine\] | | | `[NEG]` |

**ARM ③ — pyrazines (Table 3, 120 min).** All three computable fold-changes are inside 2 SD;
see §3d item 5. ARP+Cys × 1.27 (14.6 kJ mol⁻¹), ARP+Glu × 1.33 (17.4), ARP+Gly × 1.50 (24.7).
**None is a resolved measurement.**

**⭐ THE ORDERING, WHICH IS THE ROBUST PART.** In the ARP-alone arm the apparent Eₐ **ranks
cleanly by chemical class**:

> **thiophenes (2.4–17.4) < free thiols (30.3–44.4) < methylthio-ketone (52.8) <
> thiazoles (80.9) < furans (69.5–126.8)**

**The oxygen heterocycles that FEED the thiols carry apparent barriers 2–4 × larger than the
thiols themselves.** Since MFT comes from HMMF + H₂S and FFT from furfural + H₂S, and the furan
precursors have Eₐ 94.7 and 126.8 while their thiol products have 30.3 and 44.4, **the
thiolation step cannot be rate-limiting on this leg — the sulfur supply is.** That is a
structural claim the repo can test, and it is independent of every calibration concern.

### 6c. ⭐ **THE PRECURSOR-NORMALISED LADDER — the one correction Feng uniquely enables**

Feng is the only paper in this wave that measures its own substrate depletion alongside its
product yields. That converts a yield ratio into a **selectivity** ratio.

**Assumptions, stated:** t = 0 concentration is the nominal **20 mmol L⁻¹** loading (the earliest
observation is 30 min, so t = 0 is not observed). ARP-alone at 120 °C / 120 min = **1.30**
mmol L⁻¹ (from the text's 93.50 % decrease). ARP-alone at 100 °C / 120 min = **10.4 ± 0.8**
mmol L⁻¹ (`[D]`, Fig 2a, midpoint of the two ambiguous candidate curves — §4b). ARP + GSH from
the text-anchored 9.86 and 0.83 mmol L⁻¹.

| | **ARP alone** | **ARP + GSH** |
|---|---|---|
| ARP consumed, 100 °C | 20 − 10.4 = **9.6 ± 0.8** | 20 − 9.86 = **10.14** |
| ARP consumed, 120 °C | 20 − 1.30 = **18.70** | 20 − 0.83 = **19.17** |
| **consumption ratio, 120/100** | ★ **× 1.95 ± 0.16** | ★ **× 1.89** |

| product | raw × | **÷ consumption** | Eₐ of the normalised ratio `[D]` |
|---|---|---|---|
| **MFT, ARP alone** | 1.645 | ★ **× 0.84 (+0.08/−0.06)** | **−10.6 kJ mol⁻¹** |
| **FFT, ARP alone** | 2.071 | ★ **× 1.06 (+0.10/−0.08)** | **+3.6 kJ mol⁻¹** |
| sulfur class, ARP alone | 2.020 | **× 1.04** | +2.4 |
| furan class, ARP alone | 6.393 | **× 3.28** | +72.4 |
| 2-methylthiophene, ARP alone | 1.039 | **× 0.53** | −38.6 |
| **MFT, ARP + GSH** | 3.381 | **× 1.79** | +35.4 |
| **FFT, ARP + GSH** | 4.087 | **× 2.16** | +47.0 |
| sulfur class, ARP + GSH | 3.972 | **× 2.10** | +45.2 |

**⭐ Read this carefully — it is the most important derived result in the dossier.**
1. **In the ARP-alone arm, MFT and FFT yield PER MOLE OF PRECURSOR DESTROYED is flat over
   100→120 °C: × 0.84 and × 1.06.** Kang's raw × 1.12 and × 1.10 sit **inside** those bands.
   **⇒ Once extent is removed, Feng and Kang agree on the free thiols to within the
   digitisation error.**
2. **The furan class does NOT go flat under the same normalisation (× 3.28).** So the correction
   is not a trivial rescaling that flattens everything — **it flattens exactly the two compounds
   Kang said were flat, and leaves the furans steep.** That is a non-trivial, discriminating
   result.
3. **The ARP + GSH arm stays steep even after normalisation (× 1.8–2.2)** — because the added
   GSH supplies sulfur *independently* of ARP degradation, so thiol yield is no longer
   ARP-limited and the extent correction no longer applies to the sulfur term. **Consistent, and
   a good internal control on the whole argument.**
4. **⚠️ The correction rests on a digitised, partly ambiguous 100 °C ARP endpoint and an
   unobserved t = 0.** If the true 100 °C endpoint is 9.9 rather than 10.4 the consumption ratio
   becomes × 1.85 and MFT normalises to × 0.89; if it is 10.9, × 2.05 and × 0.80. **The
   conclusion "flat" is robust across that whole range; the exact value is not.** `[!]`

### 6d. Other derived quantities `[D]`

| # | quantity | value | basis | note |
|---|---|---|---|---|
| 1 | **Eₐ of GSH → Cys release, GSH-alone system** | **66.4 kJ mol⁻¹** | Fig 2h, Cys 1.80 → 5.35 mmol L⁻¹ | yield-at-120-min, unbuffered |
| 2 | **Eₐ of GSH consumption, GSH-alone system** | **≈ 34 kJ mol⁻¹** | Fig 2h, residual 13.1 → 4.7 ⇒ consumed 6.9 → 15.3, × 2.22 | ditto |
| 3 | **Eₐ of GO accumulation, ARP alone** | **+42.7 kJ mol⁻¹** | Fig 3a/3c, 0.072 → 0.145 | a *pool*, not a rate |
| 4 | **Eₐ of GO accumulation, ARP + Cys** | ★ **+4.1 kJ mol⁻¹** | Fig 3e/3f, 0.029 → 0.031 | ★ **the Cys clamp** |
| 5 | **Eₐ of 1-DO pool** | **−138 kJ mol⁻¹** | Fig 3a/3c, 0.212 → 0.022 | ★ **negative — proof a pool ≠ a rate** |
| 6 | **Eₐ from 30–60 min initial slopes, sulfur class, ARP alone** | **34.8 kJ mol⁻¹** | Fig 2g `[D]` | vs 42.9 from yields |
| 7 | **Eₐ from 30–60 min initial slopes, sulfur class, ARP + GSH** | **101.7 kJ mol⁻¹** | Fig 2g `[D]` | vs 84.1 from yields |
| 8 | **systematic spread between yield-Eₐ and slope-Eₐ** | ★ **8–18 kJ mol⁻¹** | items 6, 7 | ★ **an honest floor on the uncertainty of any Eₐ from this ladder** |
| 9 | Cys deficit vs Gly, ARP + GSH 120 °C | **22.6 : 1** (Gly 17.9 / Cys 0.79) | Fig 2e,f | ⇒ ≈ 96 % of released Cys consumed |
| 10 | sulfur : furan branching, 120 °C 120 min | ARP alone **0.257** → ARP + GSH **1.612** | T1, T2 | **× 6.27 swing on GSH addition** |
| 11 | FFT : MFT ratio, ARP alone | 100 °C **1.150** · 120 °C **1.449** | T1 | rises with T |
| 12 | FFT : MFT ratio, ARP + GSH | 100 °C **1.148** · 120 °C **1.388** | T2 | ★ **identical to ARP alone at 100 °C, and nearly so at 120 °C** |

**⭐ Item 12 is a clean cross-paper test.** `kang2026_SI_extraction.md` §5c found the FFT : MFT
branching ratio **invariant to pH (2.91 / 3.02 / 2.73)** but spanning **0.15–9.2** across
nitrogen co-substrates, concluding *"the FFT : MFT branching ratio is set by the nitrogen
co-substrate, NOT by pH."* **Feng adds a third axis: at fixed co-substrate, the ratio is also
nearly invariant to the extra-added nucleophile (1.150 vs 1.148 at 100 °C — a 0.2 % difference
across a 20 mmol L⁻¹ GSH addition) and only weakly sensitive to temperature (+26 % over 20 K).**
**⇒ Kang's claim generalises: FFT : MFT is a robust branching constant of the sugar–thiol
chemistry, ≈ 1.15 in the ribose–GSH system and ≈ 3.0 in the xylose–Cys/TTCA system.** The
2.6-fold gap between the two systems is the *sugar/precursor* effect, and it is now bracketed by
two independent papers. **STRUCTURAL, testable, and the second-strongest transferable result
here.**

---

## §7. ⭐ **THE KANG SWITCH-ON CROSS-CHECK — REQUIRED SECTION**

### 7a. The reference finding, restated

Kang et al. 2026 (TTCA 10 mmol L⁻¹ / Cys, pH 7.0, 120 min, sealed pressure vessel, oil bath,
DB-WAX, external calibration, n = 3 — `kang2026_SI_extraction.md` Table S4) reports that MFT and
FFT are **strongly non-Arrhenius**:

| | 100 °C | 120 °C | 140 °C | **low leg 100→120** | **high leg 120→140** |
|---|---|---|---|---|---|
| MFT, μg L⁻¹ | 1.237 | 1.388 | 5.907 | **× 1.12, Eₐ 6.9** | **× 4.26, Eₐ 97.9** |
| FFT, μg L⁻¹ | 3.734 | 4.107 | 11.439 | **× 1.10, Eₐ 5.8** | **× 2.78, Eₐ 69.0** |
| **sulfur class** | 13.978 | 35.866 | 60.400 | **× 2.566, Eₐ 57.5** | **× 1.684, Eₐ 35.2** |
| 2-acetylthiazole | 3.079 | 8.795 | 9.858 | × 2.856, Eₐ 64.0 | × 1.121, Eₐ 7.7 |

**Two distinct claims are bundled here and they must be separated:**
- **(K1) the free thiols SWITCH ON:** their Eₐ climbs 6.9 → 97.9 (MFT) and 5.8 → 69.0 (FFT).
- **(K2) the sulfur class as a whole DECELERATES:** Eₐ falls 57.5 → 35.2, and 2-acetylthiazole
  collapses 64.0 → 7.7.

### 7b. **Why Feng is a valid test of the LOW leg**

| design element | **Kang 2026** | **Feng 2022** | comparable? |
|---|---|---|---|
| temperature span | 100 / 120 / 140 °C | **100 / 120 °C** | ✅ **exactly the low leg** |
| hold time | **120 min** | **120 min** | ✅ identical |
| initial pH | **7.0** | **7.0** | ✅ identical |
| buffer | **none** (NaOH-adjusted) | **none** (NaOH-adjusted) | ✅ identical, and identically unmonitored |
| vessel | sealed pressure-rated glass, oil bath | **sealed pressure-resistant glass, oil bath** | ✅ identical |
| quench | ice bath | **ice water** | ✅ identical |
| precursor | **TTCA** (xylose–Cys thiazolidine MRI), 10 mmol L⁻¹ | **RG-ARP** (ribose–GSH Amadori), 20 mmol L⁻¹ | ⚠️ **different MRI, different pentose, 2 × loading** |
| sulfur source | Cys, covalently bound in TTCA | **Cys, released from GSH → Cys-Gly → Cys in situ** | ⚠️ **Feng's release is slower and rate-limiting** |
| analysis | HS-SPME-GC-MS, DB-WAX, external standards, n = 3 | **HS-SPME-GC-MS, DB-WAX, external standards, n = 3** | ✅ same class of method |
| internal standard | 1,2-dichlorobenzene | **1,2-dichlorobenzene** | ✅ **identical** |

**⇒ The two designs are matched on nine of eleven axes, including the internal standard. The two
differences are the identity of the pentose-derived MRI and the fact that Feng's cysteine must
first be liberated by peptide hydrolysis.** That second difference matters and is discussed in
§7e. **The comparison is legitimate; it is the closest-matched replication available in the
corpus.**

### 7c. **THE FULL 100 → 120 °C FOLD-CHANGE TABLE, every compound and every class subtotal**

| class | compound | **ARP alone** | **ARP + GSH** | **Kang 2026** |
|---|---|---|---|---|
| **free thiols** | **MFT** | **× 1.645** | × 3.381 | **× 1.12** |
| | **FFT** | **× 2.071** | × 4.087 | **× 1.10** |
| | *thiol subtotal* | **× 1.873** | × 3.758 | × 1.11 |
| | 3-mercapto-2-butanone | — | — | × 1.380 |
| | 3-mercapto-2-pentanone | — | — | × 1.334 |
| **thiophenes** | 2-methylthiophene | ★ **× 1.039** | × 2.590 | 0 → 0 (n/a) |
| | 2-thiophenecarboxaldehyde | 0 → 0.37 | × 19.600 ⚠️ | × 24.4 |
| | thiophene | 0 → 0 | 0 → 0.04 | 0 → 0 |
| | *thiophene subtotal* | **× 1.331** | **× 3.208** | × 2.61 |
| **thiazoles** | thiazole | 0 → 0.29 | 0 → 0.37 | × 2.005 |
| | 4-methyl-5-ethylthiazole | × 3.769 | × 4.483 | — |
| | 2-methyl-2-thiazoline | — | 0 → 0.72 | — |
| | 2-acetylthiazole | — | — | × 2.856 |
| | *thiazole subtotal* | **× 4.885** | **× 8.241** | ≈ × 3.4 |
| **thio-ketone** | 3-(methylthio)-2-butanone | × 2.376 | × 4.509 | — |
| **★ SULFUR CLASS** | **subtotal** | ★ **× 2.020** | ★ **× 3.972** | ★ **× 2.566** |
| **furans / O-het** | furfural | × 8.001 | × 6.321 | × 1.713 |
| | 4-hydroxy-5-methyl-3(2H)-furanone | × 4.727 | × 3.150 | 0 → 0 |
| | 2-methylfuran | × 3.125 | × 3.718 | 0 → 0.019 |
| | furan | 0 → 0.43 | 0 → 0 | 0 → 0.138 |
| | 1-(2-furanyl)-ethanone | 0 → 0.71 | 0 → 0.32 | 0 → 0.136 |
| | 2(5H)-furanone | — | 0 → 2.69 | 0 → 0.072 |
| **O-HETEROCYCLE CLASS** | **subtotal** | **× 6.393** | **× 4.884** | **× 1.821** |
| **pyrazines** | methylpyrazine | — | 0 → 0.04 | 0 → 0 |
| | (across Table 3 arms) | Cys × 1.27, Glu × 1.33, Gly × 1.50, GSH 0 → 0.04 | | 0 → 0 |
| **GRAND TOTAL** | | **× 4.433** | **× 4.279** | **× 2.421** |

### 7d. **⚠️ THE VERDICT, IN THE REQUIRED VOCABULARY**

**On (K2), the class-level low leg: ✅ REPLICATES.**
Kang's sulfur class runs × 2.566 (Eₐ 57.5). Feng's ARP-alone arm runs **× 2.020 (Eₐ 42.9)**;
the ARP + GSH arm runs × 3.972 (84.1). **Kang's value sits between Feng's two arms, and Feng's
closest-matched arm agrees to within 27 % in fold-change / 15 kJ mol⁻¹ in Eₐ.** Given that the
two studies use different MRIs and different pentoses, that is a **good** replication. **The
repo's low-leg class Eₐ of ≈ 43–58 kJ mol⁻¹ is now supported by two independent papers.**

**On (K1), the free-thiol flatness at face value: ❌ CONTRADICTED.**
Kang has MFT at **0.44 ×** its own class rate and FFT at **0.43 ×** — the two thiols decouple
sharply downward from everything else in the table. **Feng sees no such decoupling: MFT runs at
0.81 × its class and FFT at 1.03 × (ARP alone), 0.85 × and 1.03 × (ARP + GSH).** In Feng, FFT
simply *is* the sulfur class. **The distinctive signature of Kang's low leg — free thiols much
flatter than their own class — does not appear in Feng in raw form.**

**On (K1) after the extent correction: ✅ REPLICATES, and closely.**
Feng uniquely measures its own precursor depletion. **Per mole of ARP destroyed, Feng's MFT runs
× 0.84 (+0.08/−0.06) and FFT × 1.06 (+0.10/−0.08) over 100 → 120 °C — flat, and bracketing
Kang's raw × 1.12 and × 1.10.** The furan class does *not* flatten under the same correction
(× 3.28), so this is not an artefact of the normalisation. **⇒ The physical content of Kang's
flat low leg is corroborated. The discrepancy in raw numbers is an extent effect, and Feng
supplies the measurement that removes it.** `[!]` This depends on a digitised endpoint; see §6c.

**On the SLOPE BREAK: ⛔ UNTESTABLE — and this is the required, load-bearing statement.**
**Feng has exactly two temperatures and both lie inside Kang's low leg.** A two-point Eₐ has zero
degrees of freedom, so **it cannot detect curvature by construction**. Feng contains **no
evidence whatsoever, for or against, the 120 → 140 °C switch-on that is the whole of Kang's
non-Arrhenius claim.** ⚠️ **Any use of Feng 2022 to support, refute, or calibrate a slope break
is invalid and must be refused.**

**Additional corroborations Feng supplies, unprompted:**
- **2-methylthiophene at × 1.039 (Eₐ 2.4 kJ mol⁻¹) is the flattest compound in either paper** —
  a genuinely near-zero apparent barrier on the low leg, independently confirming that
  **flat low-leg behaviour is real in this chemistry** even though it lands on a different
  compound than in Kang.
- **The thiophene class runs × 1.33 in Feng and × 2.61 in Kang; the thiazole class × 4.89 vs
  ≈ × 3.4.** Same rank order (thiophenes flattest, thiazoles steepest), different magnitudes.
- ⚠️ **The O-heterocycle class does NOT agree.** Feng gives **× 6.39** and Kang **× 1.82** — a
  3.5-fold disagreement, the largest in the table. Feng's furan yields are dominated by furfural
  and HMMF from a *ribose*-derived Amadori compound; Kang's come from a *xylose*-derived
  thiazolidine in which the sugar is already committed to the Cys adduct. **The furan channel is
  not transferable between these two systems. Do not cross-calibrate it.**

### 7e. **Why Feng's raw thiol fold-changes are steeper than Kang's — the mechanistic reading**

In Kang, the cysteine is **already bound in TTCA** and is released by a single ring-opening. In
Feng, cysteine must be liberated by **two sequential peptide hydrolyses**
(GSH → PCA + Cys-Gly → Cys + Gly), and Fig 2h shows that release is itself strongly
temperature-dependent: **Cys goes 1.80 → 5.35 mmol L⁻¹ from 100 to 120 °C, Eₐ ≈ 66 kJ mol⁻¹.**
**⇒ Feng's thiol ladder carries an extra, ~66 kJ mol⁻¹ sulfur-supply term in series that Kang's
does not.** That alone predicts a steeper raw thiol response in Feng, and is fully consistent
with the extent-corrected agreement of §6c. **The two papers are not in conflict; they are
measuring the same thiolation chemistry through different sulfur-delivery bottlenecks.**

---

## §8. VERIFIED NEGATIVES `[NEG]` — what a reader might hope is here, and is not

| # | what is missing | why someone would look for it |
|---|---|---|
| 1 | ★ **Any third temperature.** 100 and 120 °C only. | **No curvature test, no slope break, no 3-point Arrhenius. The single largest limitation.** |
| 2 | ★ **Any rate constant, order, Eₐ, half-life or kinetic model.** | The paper is descriptive throughout; §5. |
| 3 | ★ **Any measured final pH, at any temperature, in any arm.** Only "adjusted to pH 7.0" initially, and **no buffer**. | pH drift is the largest undeclared systematic on every Eₐ here (§6a item 6). |
| 4 | ★ **Any browning, absorbance, colour or melanoidin measurement.** | Nothing at 420 nm, 294 nm or any wavelength. The paper measures no browning at all. |
| 5 | ★ **Any H₂S measurement.** H₂S is invoked as the central intermediate on pp. 9101, 9103, 9104 and **is never measured**. | Same gap `kang2026_SI_extraction.md` §6a found in Kang. **Neither sulfur paper in the corpus measures H₂S.** |
| 6 | **Any ammonia measurement.** Invoked for thiazole formation (p. 9103), never measured. | |
| 7 | **The calibration curves, LODs, LOQs and calibrated ranges for the volatiles** (Table S1, SI). | Cannot check whether reported values are inside the calibrated span, nor screen for a Kang-style impossible curve. |
| 8 | **The OAVs** (Tables S2, S3, SI). The text asserts MFT and FFT have OAV > 1 but prints **no OAV and no odour threshold**. | The repo's threshold lane gets nothing from this paper. |
| 9 | **Any 100 °C time course in Tables 1 or 2** — only the 120-min column. | *Partly recoverable*: Fig 2g gives the sulfur-class 100 °C course and Fig 3a/3b the dicarbonyls (§4b, §4c) `[D]`. |
| 10 | ★ **Any ARP-alone column in Table 3.** No pyrazine baseline at either temperature. | The paper's "promotion" claim has **no control** (§3d item 6). |
| 11 | ★ **Any volatile data for the ARP + Cys / ARP + Glu / ARP + Gly arms except pyrazines.** | The most interesting arm chemically (ARP + Cys) has **no MFT, no FFT, no thiophene, no thiazole number** — except one "(data not shown)" value, 36.41 μg L⁻¹ of 3-(methylthio)-2-butanone. |
| 12 | **The Rib + GSH comparator table** (Fig S2, SI). Only two totals appear in the text (23.06 → 70.84 μg L⁻¹). | |
| 13 | **The thiazole comparison figure** (Fig S3, SI) backing the "× 2.48" claim. | |
| 14 | **The NMR spectra** (Fig S1, SI). | ARP structural proof rests on chemical shifts quoted in text only. |
| 15 | **The reaction volume.** Never stated. | Headspace ratio and pressure in the sealed vessel are therefore unknown. |
| 16 | **The SPME fibre film thickness.** Only "DVB/CAR/PDMS". | |
| 17 | **Reagent purities** for anything except the RG-ARP (98.05 %). | |
| 18 | **Any recovery, spike or matrix-effect assessment**; any headspace-partition correction. | The μg L⁻¹ values are operational, not thermodynamic (§2c item 1). |
| 19 | **Any mass or sulfur balance.** Sulfur in = 20 mmol L⁻¹ (ARP) or 40 (ARP + GSH); sulfur out in volatiles ≈ 0.5 μmol L⁻¹. **Six orders of magnitude unaccounted, unremarked.** | |
| 20 | **Any n or replicate structure specific to Figures 2 and 3** beyond the blanket "all experiments in triplicate". | |
| 21 | **Any 1-DO / 3-DO authentic standard** — explicitly absent, hence the internal-calibration semi-quant. | §2d. |
| 22 | **A t = 0 point on any time course.** Earliest observation is 30 min, by which time 77 % of the ARP is gone at 120 °C. | §4b. |
| 23 | **Any statistical comparison BETWEEN the 100 and 120 °C columns.** The Duncan letters run *within a row across the four 120 °C times*; the 100 °C column carries its own letter but the design of the test across temperatures is not stated. | ⚠️ **The significance of every fold-change in §6b is therefore untested by the authors.** |

---

## §9. CONSOLIDATED PARAMETER TABLE

| # | parameter | value | units | condition | provenance | source anchor |
|---|---|---|---|---|---|---|
| 1 | RG-ARP loading | 20 | mmol L⁻¹ | all ARP arms | `[M]` | p. 9097 |
| 2 | co-substrate loading (GSH/Cys/Gly/Glu) | 20 | mmol L⁻¹ | equimolar | `[M]` | p. 9097 |
| 3 | initial pH | 7.0 | — | all arms, NaOH-set, **unbuffered** | `[M]` | p. 9097 |
| 4 | final pH | **not measured** | — | — | `[NEG]` | — |
| 5 | temperatures | 100, 120 | °C | oil bath, sealed vessel | `[M]` | p. 9097 |
| 6 | times | 30, 60, 90, 120 | min | 120 °C; **120 only at 100 °C** in tables | `[M]` | Tables 1–3 |
| 7 | replication | 3 | — | all | `[M]` | p. 9098 |
| 8 | RG-ARP purity | 98.05 | % | as prepared | `[M]` | p. 9097 |
| 9 | RG-ARP preparation yield | 58.00 | % of GSH | 80 °C, 40 min, pH 7.5 | `[M]` | p. 9097 |
| 10 | internal standard | 1,2-dichlorobenzene, 36 ng in 3 mL | ≈12 μg L⁻¹ | HS-SPME | `[M]`/`[D]` | p. 9097 |
| 11 | SPME extraction | 60 °C / 20 min, DVB/CAR/PDMS | — | 3 mL in 15 mL vial | `[M]` | p. 9097 |
| 12 | column | DB-Wax 30 m × 0.25 mm × 0.25 μm | — | | `[M]` | p. 9097 |
| 13 | **MFT, ARP alone** | **4.39 ± 0.22 / 7.22 ± 0.44** | μg L⁻¹ | 100 / 120 °C, 120 min | `[M]` | Table 1 |
| 14 | **FFT, ARP alone** | **5.05 ± 0.02 / 10.46 ± 0.41** | μg L⁻¹ | 100 / 120 °C, 120 min | `[M]` | Table 1 |
| 15 | **sulfur subtotal, ARP alone** | **15.21 / 30.73** | μg L⁻¹ | 100 / 120 °C, 120 min | `[M]` | Table 1 |
| 16 | **furan subtotal, ARP alone** | **18.73 / 119.74** | μg L⁻¹ | 100 / 120 °C, 120 min | `[M]` | Table 1 |
| 17 | **MFT, ARP + GSH** | **9.51 ± 0.18 / 32.15 ± 0.49** | μg L⁻¹ | 100 / 120 °C, 120 min | `[M]` | Table 2 |
| 18 | **FFT, ARP + GSH** | **10.92 ± 2.04 / 44.63 ± 0.34** | μg L⁻¹ | 100 / 120 °C, 120 min | `[M]` | Table 2 |
| 19 | **sulfur subtotal, ARP + GSH** | **30.04 / 119.32** | μg L⁻¹ | 100 / 120 °C, 120 min | `[M]` | Table 2 |
| 20 | ★ **MFT fold-change** | **× 1.645** (ARP) / **× 3.381** (+GSH) | — | 100 → 120 °C, 120 min | `[D]` | §6b |
| 21 | ★ **FFT fold-change** | **× 2.071** (ARP) / **× 4.087** (+GSH) | — | 100 → 120 °C, 120 min | `[D]` | §6b |
| 22 | ★ **sulfur-class fold-change** | **× 2.020** (ARP) / **× 3.972** (+GSH) | — | 100 → 120 °C, 120 min | `[D]` | §6b |
| 23 | ★ **apparent Eₐ, MFT** | **30.3** (ARP) / 74.3 (+GSH) | kJ mol⁻¹ | 2-point, 100→120 °C | `[D]` | §6b |
| 24 | ★ **apparent Eₐ, FFT** | **44.4** (ARP) / 85.9 (+GSH) | kJ mol⁻¹ | 2-point, 100→120 °C | `[D]` | §6b |
| 25 | ★ **apparent Eₐ, sulfur class** | **42.9** (ARP) / 84.1 (+GSH) | kJ mol⁻¹ | 2-point, 100→120 °C | `[D]` | §6b |
| 26 | apparent Eₐ, furan class | 113.1 (ARP) / 96.7 (+GSH) | kJ mol⁻¹ | 2-point | `[D]` | §6b |
| 27 | apparent Eₐ, furfural | 126.8 (ARP) / 112.5 (+GSH) | kJ mol⁻¹ | 2-point | `[D]` | §6b |
| 28 | apparent Eₐ, HMMF | 94.7 (ARP) / 70.0 (+GSH) | kJ mol⁻¹ | 2-point | `[D]` | §6b |
| 29 | apparent Eₐ, 2-methylthiophene | ★ **2.4** (ARP) | kJ mol⁻¹ | 2-point — flattest in the corpus | `[D]` | §6b |
| 30 | ★ **ARP consumption ratio** | **× 1.95 ± 0.16** (ARP) / × 1.89 (+GSH) | — | 100 → 120 °C, 120 min | `[D]` | §6c |
| 31 | ★ **MFT per mole ARP consumed** | ★ **× 0.84 (+0.08/−0.06)** | — | 100 → 120 °C | `[D]` | §6c |
| 32 | ★ **FFT per mole ARP consumed** | ★ **× 1.06 (+0.10/−0.08)** | — | 100 → 120 °C | `[D]` | §6c |
| 33 | ARP degradation | 93.50 | % | ARP alone, 120 °C, 120 min | `[M]` | p. 9100 |
| 34 | ARP residual | 9.86 / 0.83 | mmol L⁻¹ | ARP + GSH, 100 / 120 °C, 120 min | `[M]` | pp. 9101–9102 |
| 35 | PCA released | 26.83 | mmol L⁻¹ | ARP + GSH, 120 °C, 120 min | `[M]` | p. 9102 |
| 36 | Cys-Gly | 4.01 | mmol L⁻¹ | ARP + GSH, 120 °C, 120 min | `[M]` | p. 9102 |
| 37 | Cys, GSH-alone system | 1.80 / 5.35 | mmol L⁻¹ | 100 / 120 °C, 120 min | `[D]` Fig 2h | §4b |
| 38 | ★ **Eₐ, GSH → Cys release** | ★ **66.4** | kJ mol⁻¹ | 2-point, GSH alone, 120 min | `[D]` | §6d |
| 39 | GO, ARP alone | 0.072 / 0.145 | mmol L⁻¹ | 100 / 120 °C, 120 min | `[D]` Fig 3 | §4c |
| 40 | GO, ARP + Cys | 0.029 / 0.031 | mmol L⁻¹ | 100 / 120 °C, 120 min | `[D]` Fig 3e,f | §4c |
| 41 | ★ **Eₐ, GO accumulation, ARP + Cys** | ★ **4.1** | kJ mol⁻¹ | the Cys clamp | `[D]` | §6d |
| 42 | GO, ARP + Glu / + Gly at 120 °C | 0.215 / 0.165 | mmol L⁻¹ | 120 min | `[M]` | p. 9102 |
| 43 | MGO, ARP + Glu | 0.010 / 0.021 | mmol L⁻¹ | 100 / 120 °C | `[M]`/`[D]` | p. 9102, Fig 3f |
| 44 | 1-DO / 3-DO, ARP + GSH 100 °C | 0.109 / 0.038 | mmol L⁻¹ | 120 min | `[M]` **semi-quant** | p. 9102 |
| 45 | 1-DO / 3-DO, ARP + Cys 100 °C | 0.012 / 0.003 | mmol L⁻¹ | 120 min | `[M]` **semi-quant** | p. 9102 |
| 46 | pyrazine, ARP + Cys | 0.43 ± 0.02 / 0.53 ± 0.17 | μg L⁻¹ | 100 / 120 °C, 120 min | `[M]` | Table 3 |
| 47 | methylpyrazine, ARP + Cys | 0.08 ± 0.01 / 0.12 ± 0.01 | μg L⁻¹ | 100 / 120 °C | `[M]` | Table 3 |
| 48 | pyrazine total, ARP + GSH | 0.00 / 0.04 | μg L⁻¹ | 100 / 120 °C | `[M]` | Table 3 |
| 49 | ★ **FFT : MFT branching** | **1.150 / 1.449** (ARP), **1.148 / 1.388** (+GSH) | — | 100 / 120 °C, 120 min | `[D]` | §6d item 12 |
| 50 | sulfur : furan branching | 0.257 (ARP) → 1.612 (+GSH) | — | 120 °C, 120 min | `[D]` | §6d item 10 |
| 51 | Cys : Gly imbalance | 1 : 22.6 | — | ARP + GSH, 120 °C, 120 min | `[D]` Fig 2e,f | §4b |
| 52 | GO calibration | y = 4.00 × 10⁷ x − 2.91 × 10³, R² 0.9998 | — | HPLC-DAD 315 nm | `[F]` | p. 9098 |
| 53 | MGO calibration | y = 3.00 × 10⁷ x − 4.27 × 10³, R² 0.9997 | — | ditto | `[F]` | p. 9098 |
| 54 | Cys calibration | y = 3.00 × 10⁵ x + 1.78 × 10³, R² 0.9940 | — | UPLC-MS/MS | `[F]` | p. 9097 |
| 55 | RG-ARP calibration | y = 2.00 × 10⁶ x + 9.65 × 10⁴, R² 0.9926 | — | UPLC-MS/MS | `[F]` | p. 9097 |
| 56 | volatile calibration curves | **Table S1 — NOT ON DISK** | — | — | `[NEG]` | SI |
| 57 | odour thresholds / OAVs | **Tables S2, S3 — NOT ON DISK** | — | — | `[NEG]` | SI |

---

## §10. USABILITY VERDICTS

| artefact | **verdict** | reason |
|---|---|---|
| **Table 1 & Table 2 fold-changes (100 → 120 °C), all compounds and class subtotals** | ★ **USE-Q** | External-standard quantification, raster-verified, arithmetic closes on all 18 aggregates, n = 3, matched design to Kang. Qualified by: unbuffered/unmeasured pH drift; yield-not-rate; the authors ran no cross-temperature significance test. |
| **Table 1 & Table 2 absolute μg L⁻¹ values** | **RATIO-ONLY** | Headspace SPME with no partition correction, no recovery, and calibration curves unavailable (Table S1 missing). The magnitudes are operational to this method. |
| **Absolute magnitude compared with Kang 2026 or any other paper** | ⛔ **REFUSE** | Different fibre, different substrate, different loading, no cross-lab standard. See §2c item 2. |
| **Two-point apparent Eₐ (§6b), any compound** | **PRIOR-ONLY** | Zero degrees of freedom; corroborates a magnitude, tests nothing. Lower bound because of precursor depletion. Systematic floor of 8–18 kJ mol⁻¹ from the yield-vs-slope comparison (§4b). |
| ★ **Precursor-normalised thiol response (§6c): MFT × 0.84, FFT × 1.06** | ★ **STRUCTURAL** | The sharpest transferable result. Rests on a digitised, partly ambiguous 100 °C ARP endpoint and an unobserved t = 0 — but the *conclusion* ("flat") is robust across the whole candidate range. |
| **Kang low-leg class Eₐ corroboration (× 2.02 vs × 2.57)** | ★ **USE-Q** | Two independent, closely matched designs agreeing within 27 %. |
| **Kang slope-BREAK evidence** | ⛔ **REFUSE** | Only two temperatures, both on the low leg. Untestable by construction. **Never cite Feng on the break.** |
| **Fig 2g digitised 100 °C sulfur time courses** (6 new cells) | **USE-Q `[D]`** | Validated against 8 table-anchored points to ≤ 0.4 μg L⁻¹; ± 0.8 bar. Gives initial rates as well as yields. |
| **Fig 2h GSH-alone amino-acid release; Eₐ(Cys release) = 66.4 kJ mol⁻¹** | **USE-Q `[D]`** | Bar chart, well resolved, mass balance closes to 91–97 %. Still a yield-at-120-min Eₐ. |
| **Fig 3 GO and MGO absolute values** | **USE-Q `[D]`** | GO/MGO have printed external curves (R² 0.9998/0.9997) and five of twelve panel (e)/(f) readings are text-pinned. MGO on panels (a)–(d) is **REFUSE** — ≈ 100 % relative digitisation error. |
| **Fig 3 1-DO and 3-DO values** | **RATIO-ONLY** | *"Quantitated using the internal calibration method due to the lack of standards"* — response factor assumed. Ratios cancel it; absolutes do not. |
| ★ **The Cys–GO clamp: GO invariant at 0.029 → 0.031 mmol L⁻¹, Eₐ 4.1 kJ mol⁻¹** | ★ **STRUCTURAL** | Clean, model-testable, rests on the *externally calibrated* GO channel not the semi-quant deoxyosones. |
| ★ **FFT : MFT branching ≈ 1.15, invariant to added nucleophile** | ★ **STRUCTURAL** | Extends `kang2026_SI_extraction.md` §5c's pH-invariance to a third axis. Second-strongest transferable result. |
| **Sulfur : furan branching swing, × 6.27 on GSH addition** | **STRUCTURAL** | Large, unambiguous, direction-only. (The paper's own "1.98 ×" statement of this is wrong — use 1.612; §3e.) |
| **Table 3 pyrazine ordering (Cys ≫ Glu > Gly > GSH, at both temperatures)** | **STRUCTURAL** | Robust ordering, temperature-invariant. |
| **Table 3 pyrazine fold-changes (× 1.27–1.50)** | ⛔ **REFUSE** | All inside 2 SD of the printed uncertainty; the 120 °C pyrazine point has a 32 % relative SD. |
| **The "promotion by extra-added nucleophile" claim** | **PRIOR-ONLY** | **No ARP-alone control column exists in Table 3.** The effect is unquantified against a baseline. |
| **2-thiophenecarboxaldehyde Eₐ = 181.5 kJ mol⁻¹ (ARP + GSH)** | ⛔ **REFUSE** | Denominator is 0.05 μg L⁻¹, at or below LOQ; ± 0.01 swings Eₐ by ± 40 kJ mol⁻¹. |
| **2-methyl-2-thiazoline (RI 2125)** | ⛔ **REFUSE** | RI implausible by ~800 units for this structure; identification suspect despite the `S` flag. |
| **Any browning, H₂S, ammonia, OAV, threshold, LOD/LOQ, or calibration-range parameter** | ⛔ **REFUSE** | Not in this document. §8. |
| **The paper's "1.98 ×" sulfur : furan ratio** | ⛔ **REFUSE** | Fails the arithmetic audit (true value 1.612). §3e. |
| **Any pH-dependence claim** | ⛔ **REFUSE** | Single initial pH, no drift measurement, no buffer. Use `kang2026_SI_extraction.md` Table S5 for pH. |

---

## §11. DECLARED GAPS

### 11a. **The Supporting Information IS needed, and specifically for four things**

`Feng2022_supplementary.pdf` (available free at `https://pubs.acs.org/doi/10.1021/acs.jafc.2c02949`)
is **not in `data/articles/`.** It would close, in order of value to the repo:

1. ★ **Table S1 — "Calibration curves, LOD, and LOQ of volatile compounds."** This is the single
   most valuable missing artefact. It would (i) let me check whether MFT at 4.39 μg L⁻¹ and the
   0.03–0.05 μg L⁻¹ trace values are inside the calibrated span, (ii) enable the
   impossible-curve screen that caught Kang's 2-methylthiazole, (iii) supply LOQs so that the
   ⚠️-flagged 2-thiophenecarboxaldehyde Eₐ can be either rescued or definitively refused, and
   (iv) let the seven "not detected at 100 °C" cells be converted from blanks into **bounded
   upper limits**, which would turn seven refused fold-changes into seven one-sided bounds.
2. ★ **Tables S2 and S3 — "OAVs of volatile compounds… under different reaction models and
   conditions."** These carry the **odour thresholds** the paper used, for MFT, FFT and the
   thiophenes/thiazoles, in this matrix. The repo's threshold lane currently gets **nothing**
   from Feng 2022. Kang's dossier has the same hole. **Two OAV tables would be a genuine
   addition to `k4b_paired_thresholds_and_browning.md`.**
3. **Figure S2 — Rib + GSH total volatiles.** Would supply a full non-ARP comparator ladder at
   both temperatures, currently reduced to two numbers in the text (23.06 → 70.84 μg L⁻¹).
   **This is the direct-Maillard control against the MRI route and it is the natural bridge to
   the rest of the corpus.**
4. **Figure S3 — thiazole concentrations across models.** Would supply the ARP + Cys arm's
   thiazole ladder, backing the untestable "× 2.48" claim (§5).
5. Figure S1 — NMR spectra of the purified RG-ARP. Structural confirmation only; low priority.

### 11b. Gaps the SI will **not** close — permanent limitations of this study

| # | gap | consequence |
|---|---|---|
| 1 | ★ **Only two temperatures.** | **No curvature, no slope break, ever.** The switch-on question cannot be answered from this paper no matter what the SI contains. |
| 2 | ★ **No final pH, no buffer.** | An unbounded systematic on every Eₐ. Would require re-running the experiment. |
| 3 | ★ **No ARP-alone pyrazine control.** | The paper's own central claim stays unquantified. |
| 4 | **No H₂S, no ammonia.** | The mechanism's key intermediates are inferred throughout, measured never. |
| 5 | **No volatile data for the ARP + Cys / Glu / Gly arms beyond pyrazines.** | The chemically most interesting arm is a near-total blank; "(data not shown)" appears six times. |
| 6 | **No t = 0; 77 % of the ARP already gone at the first 120 °C sample.** | The fast phase of the 120 °C reaction is unobserved and unrecoverable. |
| 7 | **No cross-temperature significance test.** | Every fold-change in §6b is statistically untested by the authors. |
| 8 | **No mass or sulfur balance.** | Six orders of magnitude of sulfur unaccounted. |

### 11c. What this dossier hands the next wave

- ★ **A second, closely-matched low-leg sulfur ladder that corroborates Kang's class-level
  Eₐ (42.9 vs 57.5 kJ mol⁻¹).**
- ★ **The extent correction (§6c) that reconciles Feng's raw thiol fold-changes with Kang's
  flat low leg — and the digitised ARP-depletion data that makes it possible.**
- ★ **A hard refusal: neither paper can test the 120 → 140 °C break; only Kang's own 140 °C
  column speaks to it, and it is a single unreplicated point.** ⚠️ **The corpus currently has
  exactly ONE measurement of the high leg. That is the top declared gap this wave leaves open.**
- **Six new `[D]` cells of 100 °C sulfur-class time course, and 40 `[D]` cells of α-dicarbonyl
  time course, none of it printed in any table.**
- **Two transferable structural constants: the Cys–GO clamp (Eₐ ≈ 4 kJ mol⁻¹) and the
  FFT : MFT branching ratio (≈ 1.15 in ribose–GSH, ≈ 3.0 in xylose–Cys).**
- **A clean demonstration, from the 1-DO pool's Eₐ of −138 kJ mol⁻¹, that a fixed-time pool
  concentration is not a rate — worth citing whenever a yield ladder is read as a barrier.**
