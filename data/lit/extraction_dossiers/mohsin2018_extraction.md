# Mohsin, Schmitt, Kanzler, Epping, Flemig & Hornemann 2018 — EXTRACTION (D-glucose + L-alanine 1:1, SOLVENT-FREE, 10 min in an aluminium tray at twelve temperatures from 130 to 200 C, dialysed at 12-14 kDa; FTIR, 13C CP/MAS NMR, EPR, MALDI-ToF-MS, SEC and CHN elemental analysis)

### THE C/N SPREAD GETS WIDER, NOT NARROWER: **C/N ≈ 13:1 at 130 C rising to 21:1 at 200 C** in glucose/alanine — every value far ABOVE the trunk's predicted 8.42 to 9.94, in the opposite direction from Mundt & Wedzicha's 7.64, and driven by a mechanism decarboxylation cannot produce (carbohydrate:amine incorporation above 1:1). **And with ALANINE the structural landmarks all move up by exactly one carbon per nitrogen: the paper itself prints 9:1 for a 1:1 unit and 15:1 for a 2:1 unit, where glycine would give 8 and 14.**

**Source on disk:** `data/articles/mohsin2018.pdf` (33 pp.). **This is the ACCEPTED MANUSCRIPT, not
the typeset article**: the cover page prints "To appear in: Food Chemistry", "PII:
S0308-8146(17)31934-9", "Reference: FOCH 22092", received 27 June 2017 / revised 3 November 2017 /
accepted 30 November 2017, and every body page carries an "ACCEPTED MANUSCRIPT" watermark. **No
volume number and no page numbers are printed anywhere in the file.** The `pdftotext -layout` text
layer is a clean digital layer (double-spaced Word manuscript) and came through intact. **The
C/N sentence, which is the whole reason this paper is in the wave, was verified against a 200-dpi
raster of manuscript page 14** (`scratchpad/img/moh16top.png`) and matches the text layer word for
word. Figures 1-5 (SEC chromatograms, FTIR spectra, NMR spectra, EPR, and the proposed glyoxal
polymerisation mechanism) are **FIGURE-ONLY**. **THE SUPPLEMENTARY DATA IS NOT ON DISK** — and that
is where the elemental analysis lives: **Table S-1 (the full C/N series), Fig. S-5 and S-6 (the C,
H, N, O contents), Fig. S-2, S-3 (UV/Vis and SEC), Fig. S-4 (reactant IR), Fig. S-7 and S-8 (MALDI
spectra) and Fig. S-1 (proposed structures) are all absent.** Only the two C/N endpoints quoted in
the running text survive.

## 0. Identity

| field | value |
|---|---|
| Title | "Structural characterization of melanoidin formed from D-glucose and L-alanine at different temperatures applying FTIR, NMR, EPR and MALDI-ToF-MS" |
| Authors | **Ghassan Faisal Mohsin** [a], **Franz-Josef Schmitt** [b], **Clemens Kanzler** [a] (corresponding, clemens.kanzler@tu-berlin.de), **Jan Dirk Epping** [c], **Sabine Flemig** [d], **Andrea Hornemann** [e]. [a] Institut für Lebensmitteltechnologie und Lebensmittelchemie, TU Berlin; [b] Max-Volmer-Labor für Biophysikalische Chemie, TU Berlin; [c] Anorganische und Analytische Chemie, TU Berlin; [d] BAM Bundesanstalt für Materialforschung und -prüfung; [e] Physikalisch-Technische Bundesanstalt |
| Venue | *Food Chemistry* — **the accepted-manuscript PDF prints no volume and no pages**. Cite by DOI |
| **DOI, exactly as printed in the PDF** | **`https://doi.org/10.1016/j.foodchem.2017.11.115`** (printed on the cover page as `DOI: https://doi.org/10.1016/j.foodchem.2017.11.115`, and again in the "Please cite this article as" block, where line-wrapping breaks it as `doi: https://doi.org/10.1016/j.foodchem.` / `2017.11.115`). **Note the DOI year is 2017 while the file is named `mohsin2018`** — the article was accepted in 2017 and appeared in 2018 |
| Funding | "the Ministry of High Education and Scientific Research in Iraq … ministerial order no. 4846 in 23/09/2013". Acknowledgement thanks **Prof. Lothar W. Kroh** |
| Abbreviations as the paper uses them | **Glc** = D-glucose; **Ala** = L-alanine; **MR** = Maillard reaction; **EA** = elemental analysis; **ARP** = Amadori rearrangement product; **SEC** = size-exclusion chromatography |
| The lineage that matters | The preparation follows **Cämmerer & Kroh (1995)**, cited twice in Methods (§2.2 and §2.3) and once in the Conclusions. **`cammerer1994_extraction.md` is that paper and is already on disk**, and `species.py` already quotes its glucose/glycine pair (C/N 7 at 60 C, 9 at 100 C). Kroh is thanked in the acknowledgement. **So this is the same protocol as the corpus's existing glycine anchor, run with alanine.** See section 1 |
| Repo status before this dossier | **Not cited anywhere in `src/kinetic_core/`.** No extraction dossier existed |

## 1. Why it matters

**`src/kinetic_core/species.py` sets `MELANOIDIN_REPEAT_UNIT_CARBON = 8`,
`MELANOIDIN_REPEAT_UNIT_NITROGEN = 1`** — six carbons from 3-deoxyglucosone plus two from an intact
glycine, per Martins & van Boekel's step 9 — and the module already carries two long notes saying
that unit is falsified. The first, `MELANOIDIN_REPEAT_UNIT_FALSIFYING_MEASUREMENT`, records Mundt &
Wedzicha 2004's **7.64 ± 0.21 at 70 C**, *below* the structural floor of 8.0, with a ¹⁴C
reconstruction attributing it to about two thirds of the glycine arriving **decarboxylated**. The
second, `MELANOIDIN_REPEAT_UNIT_SAME_SYSTEM_AT_COOKING_TEMPERATURE`, records Martins & van Boekel
2003's **11 at 120 C, 15→11 at 100 C pH 6.8, 19→16 at 100 C pH 5.5**, and closes with the
instruction that "**the C/N diagnostic should be read as a spread rather than as a bound until it
runs**."

**This paper widens the spread and, more usefully, splits it into two mechanisms that push in
opposite directions.**

| source | system | T (C) | measured C/N |
|---|---|---:|---:|
| Cämmerer & Kroh 1995 (`cammerer1994_extraction.md`) | Glc/**Gly** | 60 | 7 |
| Mundt & Wedzicha 2004 (`mundt2004_extraction.md`) | Glc/**Gly**, aqueous pH 5.5, dialysed | 70 | **7.64 ± 0.21** |
| Cämmerer & Kroh 1995 | Glc/**Gly** | 100 | 9 |
| Martins & van Boekel 2003 (`martins2003c_extraction.md`) | Glc/**Gly**, pH 5.5 | 100 | 19 → 16 |
| Martins & van Boekel 2003 | Glc/**Gly**, pH 6.8 | 100 | 15 → 11 |
| Martins & van Boekel 2003 | Glc/**Gly**, pH 6.8 | 120 | 11 |
| **Mohsin 2018 (this paper)** | **Glc/Ala**, solvent-free | **130** | **≈ 13** |
| **Mohsin 2018** | **Glc/Ala**, solvent-free | **200** | **21** |
| **the trunk predicts** | Glc/Gly | — | **8.42 to 9.94** |

**Every Mohsin value is above every trunk prediction.** The lowest one, 13 at 130 C, is **1.31x the
trunk's upper prediction of 9.94 (mine)**; the highest, 21 at 200 C, is **2.11x (mine)**.

**What alanine changes — the brief's question, answered arithmetically and from the paper's own
sentence.** Alanine is C3H7NO2; glycine is C2H5NO2. **One extra carbon per nitrogen.** The paper
states its own landmarks: "In case of an incorporation of **one mole carbohydrate per mole amino
acid, the ratio should be 9:1 (three carbons of Ala and six carbons of Glc per nitrogen)** and in
case of **two mole carbohydrate per mole amino acid 15:1**." So every structural landmark moves up
by exactly 1 (mine, applying the paper's own convention to both amines):

| incorporation | Glc/**Gly** intact | Glc/**Gly** decarboxylated | Glc/**Ala** intact | Glc/**Ala** decarboxylated |
|---|---:|---:|---:|---:|
| 1 carbohydrate : 1 amino acid | **8** (the trunk's floor) | **7** | **9** (the paper prints this) | **8** |
| 2 : 1 | 14 | 13 | **15** (the paper prints this) | 14 |
| 3 : 1 | 20 | 19 | **21** (mine) | 20 |

**Two consequences the repository should record.**

1. **The trunk's floor of 8.0 is the glycine floor. In alanine it is 9.0.** A C/N diagnostic that
   compares a glycine-parameterised trunk against an alanine measurement is off by one carbon per
   nitrogen before any chemistry happens. **Subtracting 1 puts Mohsin on a glycine basis at ≈ 12 at
   130 C and 20 at 200 C (mine)** — and 12 at 130 C sits directly on top of Martins' 11 at 120 C.
   **On a glycine basis the five laboratories form a single monotonic rising series: 7 at 60 C,
   7.64 at 70 C, 9-19 at 100 C, 11 at 120 C, 12 at 130 C, 20 at 200 C.** The trunk's 8.42-9.94
   corresponds to roughly **60-100 C**, not to a 120-180 C cook. That is the sharpest thing this
   dossier can say.
2. **The high values cannot be explained by decarboxylation, because decarboxylation lowers C/N.**
   Mundt's mechanism — two thirds of the amine entering without its carboxyl — moves C/N *down*, by
   1 per decarboxylated nitrogen. Mohsin's 21:1 is *above* even the intact 2:1 unit of 15:1, and
   lands exactly on the intact **3:1** value of 21 (mine). The paper's own reading is the same:
   "**The ratio of carbohydrate to amino component in the examined melanoidin samples is closer to
   2:1 than 1:1**", and in the Conclusions, "with increasing temperature the ratio of nitrogen to
   carbon decreases indicating that **less amino acid is incorporated in the melanoidin backbone per
   mole carbohydrate**". **So the repeat unit is falsified from below at 70 C by decarboxylation and
   from above at 130-200 C by carbohydrate-rich stoichiometry, and a single fixed unit cannot do
   both.** Any wave that revisits `MELANOIDIN_REPEAT_UNIT_CARBON` needs a *temperature-dependent*
   carbohydrate:amine ratio, not just a second nitrogen pool.

**Does the paper measure how much of the amine enters DECARBOXYLATED? No.** There is no isotope
label, no ¹⁴C, no CO2 measurement, no carboxyl titration and no quantification of any carboxyl
signal. What it has is indirect and one-directional, and section 4 carries it as such: **at 130 C
the ¹³C CP/MAS spectrum shows a sharp signal at 175 ppm assigned to a carboxyl carbon "probably
belonging to the Ala residues in the melanoidin", and the whole spectrum "exhibits strong
similarities to fructosyl alanine — the Amadori rearrangement product"**, which is by construction
non-decarboxylated. **That is positive evidence for intact incorporation at 130 C and nothing
more.** Above 150 C the paper's own attribution of the growing 1717 cm⁻¹ COOH/C=O band is to
**integrated glyoxylic acid**, not to alanine carboxyl, so the FTIR cannot be read as a
decarboxylation tracker either. See flag 5.

## 2. Methods as they matter to a model

- **The pot — solvent-free, and this is not a detail.** "The melanoidins were obtained by **mixing
  Glc and Ala in a molar ratio of 1:1** and **heating for 10 min in an aluminium tray** at
  temperatures of **130, 140, 150, 152, 154, 156, 158, 160, 170, 180, 190, and 200 C** (Cämmerer &
  Kroh, 1995)." The Conclusions call it "a **solvent free system**". **There is no water, no buffer,
  no pH, no water activity, no atmosphere control, no stirring and no headspace statement anywhere
  in the paper.** Compare Mundt & Wedzicha (0.25 M each in 0.2 M acetate, pH 5.5, 70 C) and Martins
  & van Boekel (aqueous, pH 5.5 and 6.8) — **the two sources the trunk's notes already carry are
  aqueous and buffered; this one is a dry melt.**
- **Twelve temperatures, ONE time (10 min).** The fine spacing at 150-160 C (152, 154, 156, 158)
  exists because that is where colour, solubility and the spectra change. **There is no time series
  at any temperature**, so nothing here tests whether C/N is time-invariant — which matters, because
  Martins found it flat over 15-60 min at 120 C and *falling* from 15 to 11 over 30-180 min at
  100 C.
- **Dialysis, and what it means for the number.** "The solid products were ground in a mortar and
  dialysed afterwards. Batch dialysis was performed by **dissolving 5 g of melanoidin in 300 ml of
  distilled water** in dialysis tubes. **The distilled water was exchanged every 10 h, until a total
  dialysis time of 136 h** was reached. After dialysis, all samples were **freeze-dried**." The
  tubing has a **molecular weight cut-off of 12 000-14 000 Da and a pore size of 1.5-3.0 nm**. So
  the elemental analysis is of the **retained > ~12 kDa fraction**, which is the same basis as
  Mundt & Wedzicha's "MW > 12 500" — **the one axis on which these two are directly comparable.**
- **THE AUTHORS' OWN WARNING ABOUT THE LOW-TEMPERATURE SAMPLES.** From §3.1: "Domain (D) … consists
  of low molecular weight substances, **including the reactants Glc and Ala that could not be
  removed completely after dialysis. The high amounts of low molecular weight compounds in the
  samples at 130 and 140 C should be considered for the interpretation of following data.**"
  **The 13:1 at 130 C is one of exactly those two samples.** See flag 2.
- **Elemental analysis (§2.11).** "The elemental analysis of **carbon, hydrogen, and nitrogen** was
  performed on a **FlashEA 1112 Organic Elemental Analyzer** (Thermo Fisher Scientific, Dreieich).
  **The amount of oxygen was calculated on the basis that the sample contains only carbon, nitrogen,
  hydrogen, and oxygen (O = 100 % − C − H − N).** For each determination **1-3 mg** of sample were
  used and **every sample was measured at least in duplicate**. The results are presented as **means
  ± standard deviation (S.D.)**." **The SDs exist and are in Fig. S-5/S-6, which are not on disk.**
- **The C/N ratio's basis is MOLAR (atomic), not mass.** The paper puts "the **molar ratio** of the
  lost hydrogen and oxygen is around 2:1" and "the ratio of C to N is around 13:1" in the same
  paragraph, and immediately explains 9:1 as "**three carbons of Ala and six carbons of Glc per
  nitrogen**" — an atom count. **Same basis as `melanoidin_c_over_n()`.** No conversion is needed.
- **SEC (§2.10).** 1 mg in 1 mL water, 0.2 µm nylon syringe filter, method of Wegener, Kaufmann &
  Kroh 2017. **Four molecular-weight domains A, B, C, D**; 130 C is "mainly … the low molecular
  weight fraction D"; B appears ≥ 150 C and A at 200 C. **D is present in every sample except
  200 C.**
- **UV/Vis (§2.4).** 0.1 mg in 1 mL PBS pH 7.4, 280-800 nm, 1 cm quartz. **A420 rises linearly
  130-150 C then stays constant above 150 C "because the melanoidins formed under these conditions
  are partly insoluble in PBS buffer".**
- **¹³C CP/MAS NMR (§2.8).** Bruker Avance 400, 100.56 MHz ¹³C, MAS 10 kHz, 4 mm HX probe, ¹H π/2 =
  3.1 µs, TPPM decoupling, contact time 2.0 ms, recycle delay 2 s, referenced to external TMS via
  adamantane. **Spectra shown for 130, 160, 170 and 180 C only.**
- **FTIR (§2.5-2.6).** 1 mg melanoidin + 200 mg KBr pressed at ~7 t to a pellet ~0.5 mm x 1.3 cm;
  Bruker Vertex-80v + Hyperion 3000 microscope, 128x2 FPA, 15x objective; **3900-900 cm⁻¹,
  transmission, 128 co-added scans, 4 cm⁻¹ resolution**; 50 pixel spectra averaged; Table 1 band
  positions from a **Lorentz profile fitted to the 130 C melanoidin**.
- **EPR (§2.9).** Magnetech MiniScope MS 100; **339.0 mT flux density, 15.0 mT sweep, 30 s sweep
  time, 0.150 mT modulation, 10 dB attenuation**; 2.2-3.6 mg in an NMR tube; duplicate, normalised
  to sample weight, means ± S.D.
- **MALDI-ToF-MS (§2.7).** Bruker Autoflex III, linear and reflex, ~19 kV, **m/z 1-20 kDa**, DCTB
  matrix at 20 mg/mL in THF, two spotting preparations (ground steel and anchor chip), "**The target
  and matrix delivering the highest signal amplitude were chosen**" — see flag 8.

## 3. Tables re-typed

**THE PAPER HAS EXACTLY ONE NUMBERED TABLE, AND IT IS AN FTIR BAND-ASSIGNMENT TABLE.** The only
elemental table, **Table S-1**, is supplementary and **not on disk**. Marks: `[M]` measured, `[C]`
cited, `[F]` fitted.

### Table 1 (manuscript p. 30-32). "MIR modes of Glc, Ala, and the melanoidin samples prepared at 130 and 160 °C. ν: stretching; δ: deformation/bending; existent: "+"; non-existent: "-". Inhomogeneously broadened: "broad". Band positions were determined from a Lorentz profile fitted to all modes found in the Glc/Ala melanoidin sample prepared at 130 °C."

| position / cm⁻¹ | tentative assignment | Glc | Ala | 130 °C | 160 °C | Reference (as printed) |
|---|---|:-:|:-:|:-:|:-:|---|
| 3600-3200 | ν (O-H) | + | – | + | + | (Victorio, Buquiran, & del Rosario, 2007) |
| **3080** | **ν (NH⁺)** | – | **+** | **+** | **–** | (Rubinsztain, Yariv, Ioselis, Aizenshtat, & Ikan, 1986) |
| 2943 | ν (CH3) | + | + | + | + | 2970-2950 cm⁻¹ (Coates, 2000) |
| 2914 | ν (CH2) | + | + | + | + | 2935-2915 cm⁻¹ (Coates, 2000) |
| 2891-2603 | ν (C-H) | – | + | + | broad | 2900-2880 cm⁻¹, 2820-2780 cm⁻¹ (Coates, 2000) |
| **1717** | **ν (COOH) or ν (C=O)** | – | – | **–** | **+** | ν (COOH): 1725-1700 cm⁻¹ (Coates, 2000); ν (C=O): 1720-1700 cm⁻¹ (Tipson, 1968); ν (C=O): 1710 cm⁻¹ (Cämmerer & Kroh 1995) |
| **1645** | **ν (C=O) or amide I** | – | + | **+** | **–** | ν (C=O): 1670-1620 cm⁻¹ (Tipson, 1968); ν (C=O), Amide I: 1670-1620 cm⁻¹ (Mecozzi, Acquistucci, Nisini, & Conti, 2014) |
| 1622 | ν (C=O) or ν (C=C), ν (C=N) | + | + | + | + | ν (C=O) or ν (C=C): 1615 cm⁻¹ (Victorio, Buquiran, & del Rosario, 2007); ν (C=C) or ν (C=N): 1630 cm⁻¹ (Cämmerer & Kroh, 1995) |
| 1593 | δ (N-H) | – | + | + | broad | 1650-1590 cm⁻¹ (Coates, 2000), 1620-1590 cm⁻¹ (Tipson, 1968) |
| 1458 | δ (C-H) | + | + | + | + | 1460 cm⁻¹ (Tipson, 1968) |
| 1438 | δ (CH3) | + | + | + | broad | 1470-1430 cm⁻¹ (Coates, 2000) |
| 1379 | δ (CH3) | + | – | + | + | 1380-1370 cm⁻¹ (Coates, 2000) |
| 1362 | ν (C-N) or δ (O-H) | + | + | + | broad | ν (C-N): 1360-1310 cm⁻¹ (Coates, 2000) |
| 1340 | ν (C-N) or δ (O-H) | + | – | + | – | δ (O-H): 1410-1310 cm⁻¹ (Coates, 2000) |
| 1306 | ν (C-N) or δ (O-H) | + | – | + | – | — |
| 1147 | δ (C-O) | + | + | + | – | 1150 cm⁻¹, 1100 cm⁻¹ (Coates, 2000) |
| 1114 | δ (C-O) | + | + | + | – | 1150-1100 cm⁻¹ (Tipson, 1968) |
| 1104-1014 | δ (C-O) or δ (C-H) | + | + | + | + | δ (C-H): 1005-990 cm⁻¹ (Tipson, 1968) |
| 995 | δ (C-H) | + | + | + | + | 995-985 cm⁻¹ (Coates, 2000) |

All band positions `[M]`; all "+"/"–"/"broad" `[M]`; the Reference column is `[C]` throughout; the
positions themselves are `[F]` (Lorentz fits to the 130 C spectrum).

### The elemental-analysis statements, transcribed verbatim (manuscript p. 14; verified against the raster)

> "The loss of water is also supported by results of EA (Fig. S-5 and S-6, supplementary data),
> which shows an **increasing carbon content, but decreasing content of hydrogen and oxygen with
> increasing temperature**. More importantly, **the molar ratio of the lost hydrogen and oxygen is
> around 2:1**. In addition, **the nearly constant content of nitrogen** indicates that at higher
> temperatures less Ala is integrated in the corresponding melanoidins. **The ratio of carbohydrate
> to amino component in the examined melanoidin samples is closer to 2:1 than 1:1. In melanoidins
> prepared at 130 °C the ratio of C to N is around 13:1 and increases to 21:1 at 200 °C (Table S-1,
> supplementary data).** In case of an incorporation of one mole carbohydrate per mole amino acid,
> **the ratio should be 9:1 (three carbons of Ala and six carbons of Glc per nitrogen)** and in case
> of two mole carbohydrate per mole amino acid **15:1**. Altogether, the elemental composition of
> the melanoidins prepared **between 160 and 180 °C is close to the proposed melanoidin structure
> reported by Cämmerer & Kroh (1995)** (Fig. S-6, supplementary data)."

**That paragraph is the entire elemental content of the paper as it exists on disk.** Two C/N
values, one qualitative direction for C, H, N and O, one H:O loss ratio, and two structural
landmarks. **The other ten temperatures' C/N values are in Table S-1 and are unavailable.**

### Other quantitative statements in the running text

| statement, verbatim | value | mark |
|---|---|---|
| "melanoidins … formed from D-glucose/L-alanine (Glc/Ala) at varying temperatures **between 130 and 200 °C**" | 130-200 C | `[M]` |
| SEC: "The sample prepared at 130 °C **mainly consists of the low molecular weight fraction D**"; domains "**B (≥ 150 °C)** and **A (200 °C)**" are formed with rising temperature | — | `[M]` |
| UV/Vis: "Between **130 and 150 °C** the browning of the water soluble fraction **increases linearly**, but **above 150 °C the absorbance of the extracts at 420 nm stays constant**" | — | `[M]` |
| SEC/colour: cites Borrelli et al. 2002 that "the high molecular weight fraction from **Glc/Gly** melanoidin was responsible for **80 %** of the total brown colour" | 80 % | `[C]` |
| NMR, 130 C: "two sharp peaks at around **15 ppm** and **175 ppm**. The signal at 15 ppm can be assigned to a **methyl carbon** and the signal at 175 ppm to a **carboxyl carbon, both probably belonging to the Ala residues**"; carbohydrate moiety "broad signals at around **45 ppm** and **70 ppm**. **C-2 of Ala** is to be found at around 70 ppm, too" | — | `[M]` |
| NMR: "the spectrum exhibits **strong similarities to fructosyl alanine – the Amadori rearrangement product (ARP) of Glc and Ala**" | — | `[M]` |
| NMR, rising T: new signals "**between 100 ppm and 160 ppm**", "characteristic for **sp² hybridized carbons** as found in conjugated double bond systems" | — | `[M]` |
| EPR: "The radical signal … **increases from 140 °C to 200 °C** (Fig. 4)" | — | `[M]`, magnitudes figure-only |
| MALDI: "Melanoidin samples prepared at **150 °C and 180 °C** exhibit polymer patterns … with **mass differences of m/z 74.11 (150 °C) and m/z 58.13 (180 °C)**. At 150 °C the peaks occur **between m/z 2500-5500** and at 180 °C between **m/z 1500-2500**" | 74.11 / 58.13 | `[M]` |
| MALDI assignment: **glyoxal (m/z 58)** and **glyoxylic acid (m/z 74)** | — | `[C]` (Thornalley 1999; Novotný 2008) |
| Conclusions: "these polymers with molecular weights **between m/z 1800-4000** are suspected to be **only by-products**" | — | `[M]`/`[F]` — note this range disagrees with the 1500-5500 span in §3.6, see flag 8 |

### Derived numbers (mine, arithmetic on the printed values — NOT the paper's)

**The structural landmark table in section 1 is mine** except for the two entries the paper prints
(Ala 1:1 = 9, Ala 2:1 = 15). Everything else in it is one carbon added or removed per nitrogen.

- **21:1 is exactly the intact 3:1 carbohydrate:alanine value** (3 x 6 + 3 = 21, mine). The paper
  stops at 2:1 = 15 and does not extend the series; the 200 C measurement lands on the next rung.
- **13:1 sits between the intact 1:1 (9) and 2:1 (15) landmarks**, at a nominal
  **1.67 carbohydrates per alanine (mine, linear interpolation: (13 − 3)/6 = 1.67)**. The paper's
  "closer to 2:1 than 1:1" is consistent with this.
- **21:1 corresponds to 3.0 carbohydrates per alanine (mine, (21 − 3)/6 = 3.00)**.
- **Mohsin on a glycine basis (mine, subtracting 1 C per N)**: **≈ 12 at 130 C, 20 at 200 C.**
  This conversion is valid **only** if the carbohydrate:amine stoichiometry is unchanged between the
  two amines, which nothing measures. Flagged.
- **Against the trunk's 8.42-9.94**: 13/9.94 = **1.31x**, 21/9.94 = **2.11x**, 13/8.42 = **1.54x**,
  21/8.42 = **2.49x** (all mine).
- **Against the alanine structural floor of 9.0** (mine, from the paper's own 9:1): 13/9 = **1.44x**,
  21/9 = **2.33x**. **Unlike Mundt's 7.64, neither Mohsin value is below its floor** — this paper
  does not falsify a floor, it falsifies a *ceiling* assumption.

## 4. Numbers the repository can use

**All rows: D-glucose + L-alanine, 1:1 molar, SOLVENT-FREE, 10 min in an open aluminium tray, then
batch-dialysed 136 h against distilled water through 12-14 kDa tubing and freeze-dried. CHN by
FlashEA 1112, O by difference, 1-3 mg per determination, at least duplicate, means ± SD (the SDs
are in Fig. S-5/S-6 and are NOT on disk). No pH, no water activity, no atmosphere, no time series.**

### The C/N series — the brief's first question

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| **melanoidin C/N at 130 C** | **"around 13:1"** | mol C / mol N | Glc/Ala 1:1, solvent-free, 10 min, dialysed > ~12 kDa | manuscript p. 14, running text (Table S-1 not on disk) | **measured_ratio** — but see flags 1 and 2 |
| **melanoidin C/N at 200 C** | **"21:1"** | mol C / mol N | as above | manuscript p. 14 | **measured_ratio** |
| **melanoidin C/N at 140, 150, 152, 154, 156, 158, 160, 170, 180, 190 C** | **NOT AVAILABLE** — measured, tabulated in **Table S-1**, which is supplementary and **not on disk** | — | — | — | — |
| C, H, N, O weight percentages at any temperature | **NOT AVAILABLE** — Fig. S-5 and S-6, not on disk | — | — | — | — |
| direction of C with temperature | **increasing** | — | 130-200 C | p. 14 | **level_only** |
| direction of H and O with temperature | **decreasing** | — | 130-200 C | p. 14 | **level_only** |
| direction of N with temperature | "**nearly constant** content of nitrogen" | — | 130-200 C | p. 14 | **level_only** — note this is a weight fraction, and it is the reason C/N rises |
| molar ratio of lost H to lost O | **≈ 2:1** | mol/mol | 130-200 C | p. 14 | **measured_ratio** — i.e. the mass loss is water, as the paper argues |
| carbohydrate : amino incorporation | "**closer to 2:1 than 1:1**" | mol/mol | 130-200 C, whole sample set | p. 14 | **measured_ratio**, qualitative |
| implied carbohydrates per alanine | **1.67 at 130 C, 3.00 at 200 C** (mine, from the paper's own 9:1 / 15:1 convention) | mol/mol | as above | derived | **derived_assumption (mine)** |
| best structural match | "the elemental composition of the melanoidins prepared **between 160 and 180 C** is close to the proposed melanoidin structure reported by **Cämmerer & Kroh (1995)**" | — | 160-180 C | p. 14 | **level_only** — the Cämmerer structure is two 3-deoxyglucosones plus one amino acid, per the Conclusions |

### The structural landmarks — what alanine changes

| quantity | value | source | evidence class |
|---|---|---|---|
| **C/N of an intact 1 Glc : 1 Ala unit** | **9:1** | **printed by the paper**, p. 14 | **derived_assumption** (a stoichiometric identity, stated by the authors) |
| **C/N of an intact 2 Glc : 1 Ala unit** | **15:1** | **printed by the paper**, p. 14 | derived_assumption |
| C/N of an intact 3 Glc : 1 Ala unit | 21:1 (mine) | derived | derived_assumption (mine) |
| C/N of a **decarboxylated** 1:1 Ala unit | 8:1 (mine) | derived | derived_assumption (mine) — **numerically equal to the trunk's glycine floor, for an unrelated reason** |
| **the glycine offset** | **−1 C per N at every landmark** (mine) | derived from Ala C3 vs Gly C2 | **derived_assumption (mine)** — the single most transferable statement in this dossier |
| Mohsin's series on a glycine basis | ≈ 12 at 130 C, 20 at 200 C (mine) | derived | **derived_assumption (mine)** — valid only at constant carbohydrate:amine stoichiometry |

### Decarboxylation — the brief's third question

| quantity | value | anchor | evidence class |
|---|---|---|---|
| **fraction of alanine incorporated DECARBOXYLATED** | **NOT MEASURED.** No isotope label, no ¹⁴C, no CO2 determination, no carboxyl titration, no quantified carboxyl signal at any temperature | — | — |
| carboxyl carbon present in the 130 C melanoidin | **YES** — "a sharp peak at around 175 ppm … assigned … to a **carboxyl carbon**, probably belonging to the **Ala residues** in the melanoidin" | ¹³C CP/MAS, manuscript p. 13; Fig. 3 | **level_only** — presence, not amount, and "probably" is the paper's own hedge |
| the 130 C polymer resembles the intact Amadori product | **YES** — the spectrum "exhibits **strong similarities to fructosyl alanine – the Amadori rearrangement product (ARP) of Glc and Ala**", and the paper concludes the melanoidin "most likely includes its precursors in form of Glc, Ala, the corresponding ARP, and Maillard intermediates with **intact carbohydrate backbone**, derivatized with Ala" | p. 13 | **level_only** — fructosyl alanine is by definition non-decarboxylated |
| COOH/C=O band at 1717 cm⁻¹ | **absent at 130 C, present at 160 C**, and "**grows in strength with rising temperature**" | Table 1; §3.3 | **level_only** — **and the paper attributes it to integrated glyoxal/glyoxylic acid, NOT to alanine carboxyl** (Conclusions), so it is NOT a decarboxylation readout |
| NH⁺ band at 3080 cm⁻¹ | present at 130, 140, 150 C; **absent above** | Table 1; §3.3 | level_only — bears on the amine's protonation state, not on its carboxyl |
| amide I / C=O at 1645 cm⁻¹ | present at 130 C, **vanishes above 150 C** | Table 1; §3.3 | level_only |
| **the mechanism the paper does invoke for falling N** | "at lower temperatures the linkage of the carbohydrates is realized by **integration of amino acids**, whereas at high temperatures **aldol (condensation) reactions without any involvement of amino compounds are favoured**"; alternatively "the **polymerization of low molecular weight MR products, such as glyoxal and glyoxylic acid**, and their integration in the melanoidins" | Conclusions | **derived_assumption** — two competing hypotheses, neither tested |

**Stated plainly, because the brief asks it directly: this paper measures NOTHING that bears
quantitatively on how much of the amine enters decarboxylated.** It has one qualitative datum
pointing the other way (an intact carboxyl at 130 C), and its explanation for the rising C/N is
amine-free carbon addition, which is a different mechanism from decarboxylation and moves C/N in the
opposite direction from Mundt's.

### Everything else the paper measures

| quantity | value | conditions | anchor | evidence class |
|---|---|---|---|---|
| A420 of the water-soluble fraction | rises **linearly 130-150 C**, then **constant above 150 C** (attributed to insolubility, not to saturation of colour formation) | 0.1 mg/mL PBS pH 7.4 | §3.1; Fig. S-2 not on disk | **level_only** — **a browning ceiling that is an artefact of solubility is directly relevant to any A470-based browning readout** |
| molecular-weight shift | four SEC domains A-D; **D dominant at 130 C**, **B forms ≥ 150 C**, **A at 200 C**; D absent only at 200 C | 1 mg/mL water, 0.2 µm filtered | §3.1; Fig. 1 | **level_only** |
| radical concentration | **increases monotonically from 140 C to 200 C** | 2.2-3.6 mg solid | §3.5; Fig. 4 | **level_only** — magnitudes are figure-only |
| MALDI repeat masses | **m/z 74.11 at 150 C**, **m/z 58.13 at 180 C** | DCTB, 19 kV | §3.6 | **measured_ratio** (mass differences), assigned to **glyoxylic acid** and **glyoxal** respectively `[C]` |
| MALDI mass windows | m/z **2500-5500** (150 C), **1500-2500** (180 C) | as above | §3.6 | level_only — **contradicted by the Conclusions' "1800-4000"**, flag 8 |
| **any rate constant, half-life or activation energy** | **NOT PRESENT** — one time point (10 min) at every temperature | — | — | — |
| **any pH, water activity, or moisture** | **NOT PRESENT** | — | — | — |
| **any yield, concentration or mass balance** | **NOT PRESENT** — the paper is entirely structural | — | — | — |

## 5. Flags

1. **"Around 13:1" is an approximation, and the exact value is in a supplementary table that is not
   on disk.** The paper writes "around 13:1" for 130 C and a bare "21:1" for 200 C, and points at
   Table S-1 for both. **The ten intermediate temperatures — including 150, 160, 170 and 180 C,
   which are the ones a cooking model actually needs — are unavailable.** Nothing here can be
   interpolated. **This is the single highest-value retrieval in the wave: Table S-1 is one
   supplementary file away.**
2. **THE AUTHORS THEMSELVES WARN AGAINST THE 130 C DATUM.** §3.1: "the reactants Glc and Ala …
   could not be removed completely after dialysis. **The high amounts of low molecular weight
   compounds in the samples at 130 and 140 C should be considered for the interpretation of following
   data.**" The 13:1 is a 130 C sample. **The contamination pushes C/N in both directions at once** —
   residual glucose is C6 with no nitrogen (raises C/N without bound), residual alanine is C3 with
   one nitrogen (C/N = 3, lowers it) — and the paper quantifies neither. **The low-temperature end of
   this series, which is the end the trunk cares about, is the untrustworthy end.**
3. **Solvent-free is a different chemistry from every other C/N source in the corpus.** Mundt &
   Wedzicha: 0.25 M each in 0.2 M acetate at pH 5.5. Martins & van Boekel: aqueous at pH 5.5 and
   6.8. Cämmerer & Kroh: the protocol this paper copies. **This paper has no water at all, no pH,
   and no buffer**, and its own introduction lists "temperature, reaction time, **water content and
   pH value**" as having "a strong influence on the structure of the resulting melanoidins". A dry
   melt at 200 C dehydrates far harder than an aqueous system — which is exactly what the H:O = 2:1
   loss ratio says — so **the rising C/N is partly a dehydration effect that an aqueous model will
   not reproduce.**
4. **One time point.** 10 min, at every temperature. Martins measured C/N *falling* from 15 to 11
   over 30-180 min at 100 C, so C/N is not necessarily a state function of temperature alone.
   **Nothing here constrains the time axis, and 10 min is short compared with most of the corpus's
   cooking windows.**
5. **The one carboxyl observation is a "probably" and it is not quantified.** The 175 ppm assignment
   reads "**probably** belonging to the Ala residues in the melanoidin", and a carboxyl carbon at
   175 ppm in a Glc/Ala melt could equally be glyoxylic acid, a sugar acid, or a formate — the paper
   itself proposes glyoxylic acid incorporation at higher temperatures. **The signal is never
   integrated, never compared between temperatures, and never converted to a per-nitrogen carboxyl
   count.** It is evidence that *some* intact carboxyl exists at 130 C, and no more.
6. **The temperature range does not overlap the trunk's operating window at its lower end.** The
   series starts at 130 C. The trunk's browning and its C/N prediction are exercised at 100-150 C;
   Mundt is at 70 C and Martins at 100-120 C. **Between 70 C and 130 C this paper says nothing**, and
   between 130 C and 200 C nothing else in the corpus says anything. The two datasets abut rather
   than overlap.
7. **The +1 carbon-per-nitrogen glycine correction is mine and it is only exactly right at fixed
   stoichiometry.** Alanine's extra methyl also changes its Strecker chemistry (acetaldehyde rather
   than formaldehyde), its solubility in a melt, and its reactivity; nothing here measures whether
   Glc:amine incorporation is the same for the two amines. **Use the correction to compare trends,
   not to transfer a number.**
8. **Two internal inconsistencies in the MALDI account.** §3.6 gives the polymer windows as m/z
   2500-5500 (150 C) and 1500-2500 (180 C); the Conclusions give "molecular weights between m/z
   1800-4000". And the method states the target and matrix were chosen because they gave "**the
   highest signal amplitude**", i.e. optimised for signal rather than for representativeness — the
   Conclusions duly concede these polymers "are suspected to be **only by-products** that probably
   cross-link heterogeneous melanoidin moieties that **could not be characterized using MALDI-ToF-MS,
   yet**". **Nothing in the MALDI section is a composition of the melanoidin.**
9. **Oxygen is by difference, not measured.** "O = 100 % − C − H − N". Any sulfur, ash, residual
   sodium or adsorbed water is silently booked as oxygen. **This does not affect C/N**, which uses
   only the two directly measured elements, but it does affect the H:O = 2:1 argument.
10. **The A420 plateau above 150 C is a solubility artefact, and the paper says so.** "the
    absorbance of the extracts at 420 nm stays constant, **because the melanoidins formed under
    these conditions are partly insoluble in PBS buffer and these fractions cannot contribute to the
    browning of the corresponding extracts**." **Any model calibrated on an A420 or A470 readout in a
    high-temperature dry system is measuring the soluble fraction only.** That is a general warning
    for the trunk's browning response, independent of C/N.
11. **This is the accepted manuscript, not the version of record.** The publisher's own boilerplate
    on the cover page: "**during the production process errors may be discovered which could affect
    the content**". No volume or page numbers exist to cite. Any number here should be re-checked
    against the typeset article before it is installed.
12. **What to request.** (i) **Table S-1** — the full twelve-temperature C, H, N, O and C/N table
    with its standard deviations. This one file would turn two endpoints into a twelve-point series
    across 130-200 C and would be, by a wide margin, the best-resolved C/N temperature dependence in
    the corpus. Also Fig. S-5 and S-6. (ii) The same protocol run with **glycine** instead of
    alanine, which would measure the +1 offset rather than assuming it. (iii) An isotope-labelled or
    CO2-trapping version, which is the only thing that would let a decarboxylated fraction be
    compared with Mundt's 0.289 : 0.662 : 1.
