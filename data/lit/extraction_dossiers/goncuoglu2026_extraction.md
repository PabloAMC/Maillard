# Göncüoğlu Taş, Balagiannis, Ghawi, Gökmen & Parker 2026 — EXTRACTION (twin-screw extrusion of cornmeal + pea protein isolate, 0-70 % PPI, 15-17 % moisture, 400-500 rpm, barrel to 145-160 °C; acrylamide, 42 volatiles including 15 pyrazines, and physical properties)

### NOT A KINETIC PAPER: the same first author and the same Hacettepe laboratory as the 2016 hazelnut study, but ten years later, in Reading, on an extruder — there is not one rate constant, not one activation energy and not one reaction network in it, and its whole quantitative contribution to this repository is one acrylamide level, two raw-material precursor concentrations, and a formation-versus-composition ordering for pyrazines that the engine cannot yet be asked to reproduce.

**Source on disk:** `data/articles/Goncouglu2026.pdf` (Food Chemistry: X **36** (2026) 104019, open
access CC BY, 12 pp. of body plus references).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/Goncouglu2026.txt`, 1222 lines). **Tables 2 and 4 came through clean and are
re-typed in full below. Tables 1 and 3 are GARBLED** — the text layer transposes them into a
vertical ribbon that interleaves row labels with values from several columns at once, so column
membership cannot be recovered without rendering the page. What is legible from them is transcribed
as isolated statements, marked as such, and **no number from Table 1 or Table 3 is assigned to a
named sample by me.** **Supplementary Material Tables S1-S4 are NOT on disk** — S1 (retention
indices and quantifier ions), **S2 (the free amino acid and sugar profile of the raw materials)**,
**S3 (the semi-quantified volatiles of the raw materials and the composition series)** and **S4 (the
volatiles of the screw-speed and temperature series)**. Since every volatile number in this paper
lives in S3 and S4, **the entire aroma dataset is off disk**, and Figure 1 (acrylamide) and Figure 2
(the two PCA plots) are images.

## 0. Identity

| field | value |
|---|---|
| Title | "Effect of extrusion parameters and feed composition on physical characteristics, aroma profile and acrylamide content in pea protein-enriched corn extrudates" |
| Authors | Neslihan Göncüoğlu Taş ᵃ, Dimitris P. Balagiannis ᵇ, Sameer Khalil Ghawi ᵇ, Vural Gökmen ᵃ, **Jane K. Parker ᵇ (corresponding, j.k.parker@reading.ac.uk)**. ᵃ Food Quality and Safety (FoQuS) Research Group, Department of Food Engineering, Hacettepe University, 06800 Beytepe, Ankara, Türkiye. ᵇ Department of Food and Nutritional Sciences, University of Reading, Whiteknights, Reading RG6 6DZ, UK |
| Venue | **Food Chemistry: X 36 (2026) 104019**, part of the special issue "Acrylamide Research". Received 10 February 2026, revised 1 May 2026, accepted 21 May 2026, online 22 May 2026 |
| DOI | `10.1016/j.fochx.2026.104019` |
| Licence | CC BY 4.0, open access |
| Naming | CM = cornmeal; PPI = pea protein isolate; M15/M16/M17 = 15/16/17 % feed moisture; SS400/450/500 = screw speed in rpm; T1-T4 = the four barrel temperature profiles; E1-E19 = the nineteen extrudates; SME = specific mechanical energy |
| Companions on disk | `goncuoglu2016_extraction.md`, `goncuoglutas2016_extraction.md`, `goncuoglutas2017_extraction.md` (the hazelnut paper, all three on `Goncouglu2016.pdf`); `balagiannis2009_extraction.md`, `balagiannis2010_extraction.md`, `balagiannis2015_extraction.md` (the second author's own kinetic-modelling work, which is where a reader looking for constants from this group of people should go); `parker2013_extraction.md`; `hamzalioglu2026_extraction.md` and `kocadagli2021_extraction.md` (the Hacettepe lane) |

## 1. Why it matters

**The brief was to say what is new against the 2016 paper. The honest answer is: everything and
nothing.** Everything, because it is a different food (an extruded corn-pea snack, not a roasted
hazelnut), a different process (a twin-screw extruder at 3.4 kg/h, not an oven), a different
question (product quality and acrylamide safety, not mechanism), a different set of responses
(acrylamide, 42 volatiles, bulk density, hardness, colour — not sugars and dicarbonyls), a
different second institution (Reading joins Hacettepe), and a different statistical apparatus
(ANOVA with Tukey and PCA, not multiresponse regression). Nothing, because **the one thing the 2016
paper gave this repository — a table of rate constants from a second laboratory's multiresponse fit
— has no counterpart here at all.** There is no reaction network, no differential equation, no rate
constant, no activation energy, no reference temperature, no Athena Visual Studio, no determinant
criterion and no highest-posterior-density interval anywhere in the paper. The authors say so
themselves when the composition effects get complicated: "Explaining acrylamide formation in the
context of changing feed composition, including variations in asparagine, reducing sugars, and lipid
oxidation products, is challenging. **Therefore, model system studies may be more appropriate for
elucidating the underlying mechanism.**"

So this dossier is short by the paper's own nature, and it exists for three reasons.

1. **It closes a filename ambiguity.** `data/articles/` holds `Goncouglu2016.pdf` and
   `Goncouglu2026.pdf`, and three dossiers already point at the first. Anyone reaching for
   "Göncüoğlu 2026" expecting a follow-up multiresponse fit should be told in one line that it is
   not one. **The 2016 hazelnut paper remains the only Göncüoğlu Taş multiresponse fit on disk.**
2. **It supplies precursor levels and an acrylamide level in a real low-moisture food at a
   temperature the acrylamide lane cares about.** The barrel runs to 145-160 °C at 15-17 % moisture,
   which is the low-moisture, high-temperature corner that `parameters_acrylamide.py`'s conditions
   module is built around, and the raw-material asparagine and reducing-sugar contents are printed
   in the body text (not only in the absent Table S2). See section 4.
3. **It is a directional test the pyrazine lane could one day be scored against.** Wave B18
   (`network.PYRAZINE_REACTIONS`) makes pyrazine, methylpyrazine and 2,5-dimethylpyrazine from
   glyoxal and methylglyoxal by Strecker deamination and aminoketone condensation, and B21's
   outcome section records that after the glyoxal fix the predicted ordering turns to
   **parent > methyl > 2,5-dimethyl**. This paper identifies **fifteen pyrazines**, finds **none in
   either raw material**, finds **nine in the corn-only extrudate**, and reports that total pyrazine
   peaks at 30 % PPI and falls at 70 % because the sugar side is diluted. **Every one of those
   numbers is in the absent Table S3/S4**, so the claim is available only as an ordering in prose —
   which is why section 4 classes it `level_only` and `within_study_ratio` and not as a benchmark.

What this paper does **not** give the repository: any rate constant; any barrier; any reference
temperature; any concentration-time course (there is no time axis anywhere — residence time is never
reported as a number, only discussed); any sugar or amino acid measurement of the *extrudates*
(only of the two raw materials, and only in the absent Table S2, quoted selectively in prose); any
dicarbonyl measurement of any kind; any browning measurement beyond L\*a\*b\*; and any absolute
volatile concentration (the volatiles are **semi-quantified** against three internal standards, so
they are areas normalised to a standard, not concentrations — section 2).

## 2. Methods as they matter to a model

- **Raw materials.** **Pea protein isolate**, 80 % protein, Nutraceuticals Group Europe (Merstham,
  UK). **Cornmeal medium**, "67 % carbohydrate, 13 % fat, 7 % protein, 0.4 % fibre", Natco Foods
  (Buckingham, UK). Acrylamide-d3 from Toronto Research Chemicals; everything else Sigma-Aldrich or
  Merck.
- **Feed preparation.** 1 kg batches of CM:PPI at **100:0, 85:15, 70:30, 50:50 and 30:70 (w/w)**,
  dry-mixed by whisk, then 300 mL drinking water added in a Kenwood mixer, sieved to a uniform
  particle size, **dried in a tray oven at 70 °C for 24 h**. That drying step is itself a thermal
  history the model would have to carry, and the authors attribute the acrylamide already present in
  the raw materials to "conditions to which they were exposed during drying process".
- **Extruder.** Thermo Fisher **Process 16 Hygienic** twin-screw, **3 mm die**, eight
  temperature-controlled zones, screw **L/D = 40**, pressure probe at the exit zone. Dry feed into
  zone 1 by gravimetric feeder; liquid into zone 2 by peristaltic pump.
- **The constant conditions of the composition series.** Barrel profile
  **40-60-80-100-120-140-145-145 °C** (this is T1), **screw speed 400 rpm**, **dry feed rate
  3.4 kg/h**. Varied: CM:PPI ratio (five levels) and feed moisture (**15, 16, 17 %**), with the
  liquid feed rate adjusted to hit the target moisture.
- **The process series.** Screw speed **400, 450, 500 rpm** and four barrel profiles —
  **T1 up to 145 °C, T2 up to 150 °C, T3 up to 155 °C, T4 up to 160 °C** — on the 50:50 and 30:70
  feeds at 15 % moisture. T4 is printed in full in the discussion: **40-60-80-110-140-155-160-160 °C**.
- **THE TEMPERATURE IS A BARREL SET-POINT, NOT THE MATERIAL TEMPERATURE.** The methods say melt
  temperature — "the actual temperature of the material inside the barrel" — was recorded for each
  run, and it appears in Tables 1 and 3, **both of which are garbled in the text layer**. So the
  numbers that would let a kinetic model be run against this paper (the melt temperature and the
  residence time) are on the page but not legible from the text layer, and residence time is never
  reported at all. **This is the single reason no benchmark row can be built from this paper**
  (Flags 1).
- **Moisture.** Sartorius MA150 infrared analyser at 105 °C, triplicate, on both feed and
  extrudate. Feed moisture **15-17 %**; extrudate moisture appears in the garbled Table 1 (values of
  the form "9.4 ± 0.4", "9.7 ± 0.7", "10.3 ± 0.6" are legible in the ribbon and are of the right
  size to be extrudate moisture in %, but I cannot say which sample each belongs to). **No water
  activity is measured or quoted anywhere in the paper.**
- **pH.** **Not measured.** There is no pH number in this paper.
- **Sugars (raw materials only).** HPAEC-PAD: Thermo Dionex ICS-6000, gold working electrode, PdH
  reference, gold-carbo-quad waveform, IntAmp mode; CarboPac PA-1 (4 × 250 mm) with PA1 guard at
  20 °C; 1 mL/min; multi-step NaOH gradient (A and D water, B 16 mM NaOH, C 250 mM NaOH; 35 % A /
  10 % C / 55 % D for 25 min, then 35 % A / 50 % C / 15 % D for 5 min, held 10 min, back to initial
  in 5 min, held 5 min; **70 min total**); 20 µL injection; triplicate. Calibration: glucose and
  fructose **0.1-2 mg/L**, sucrose and maltose **1-10 mg/L**. Extraction: 0.25 g in 5 mL water with
  125 µL each Carrez I and II, 9500 g for 5 min, 0.2 µm nylon.
- **Free amino acids (raw materials only).** Three extractions (10 + 5 + 5 mL water), vortex 5 min
  and 6654 g for 5 min each; 0.2 mL of pooled extract + 0.8 mL acetonitrile, centrifuged, 0.2 µm
  PTFE. Agilent 1200 HPLC + 6410 triple quadrupole, ESI positive, MRM, conditions from Kocadağlı
  2013. External calibration **1 to 50 µM** in water:acetonitrile (20:80).
- **Acrylamide.** Triple-stage extraction of **100 mg ground sample** with 10 mM formic acid in
  water containing **10 µg/L acrylamide-d3**, 2 mL of extraction solution in total, with a Carrez
  clarification in the first step; clean-up on an **Oasis MCX** cartridge; **Waters Xevo TQ-S cronos
  triple quadrupole with Acquity H Class Plus UPLC**. MRM for acrylamide-d3 m/z **75 → 58**
  (fragmentor 70 V, collision energy 9 V); Q1 unit resolution, Q2 wide. Linear calibration
  **0.1-10 µg/L** with the d3 internal standard. **LOD 0.08 µg/kg, LOQ 0.26 µg/kg**
  (S/N of 3 and 10). This is a proper isotope-dilution LC-MS/MS method — far better than the
  HPLC-UV of Knol 2005 — and its detection limit is two orders of magnitude below the values
  measured.
- **Volatiles: SEMI-QUANTIFIED, NOT QUANTIFIED.** 1 g ground extrudate in a 20 mL headspace vial +
  2 mL saturated NaCl (35 %) containing **three internal standards: 3-methyl-2-butanone and
  3-furaldehyde at 0.5 mg/L and isopropylpyrazine at 0.05 mg/L**. CTC PAL autosampler, SPME, Agilent
  7890A GC, helium at 1.2 mL/min, **DB-5MS 30 m × 0.25 mm × 0.25 µm**, oven 40 °C for 2 min then to
  **280 °C at 4 °C/min**. Assignment of standards: "**The area of 3-furaldehyde for furans,
  isopropylpyrazine for pyrazines, 3-methyl-2-butanone for Strecker aldehydes and the rest of the
  compounds were used for calculations.**" So every reported volatile is a **peak-area ratio to one
  of three standards**, not a concentration, and no response factors were determined. **No volatile
  number from this paper may be compared with a modelled concentration** (Flags 3).
- **Physical measurements.** Bulk density of feed by 250 mL cylinder, triplicate; bulk density of
  extrudate by glass-bead volumetric displacement, five replicates. Texture (hardness in N,
  fracturability in mm) and colour L\*a\*b\* as tabulated. Statistics: ANOVA with Tukey at
  p = 0.05, Student's t-test at p ≤ 0.05 for the raw-material comparison, PCA in XLStat.
- **Reference temperature.** **There is none, because there is no fitted constant.** The house rule
  that a constant fitted at a stated reference temperature must carry that temperature is satisfied
  here vacuously: nothing in this paper is fitted.

## 3. Tables re-typed

### Table 2 (clean). "Physical properties of cornmeal-pea protein isolate extrudates with different compositions at varying moisture contents"

Superscript letters are Tukey groupings within a column; footnote ᵃ (printed under Table 4 and
applying to both): "The values within the same column followed by the same letters are not
significantly different (p = 0.05) according to Tukey's test." Footnote ᵇ: "E: extrudate, CM:
cornmeal, PPI: pea protein isolate, M: moisture content." The diameter column is printed with a
"+" where a "±" is clearly meant, and is transcribed as printed.

| | sample | Bulk Density (g/mL) | Diameter (mm) | Hardness (N) | Fracturability (mm) | L* | a* | b* |
|---|---|---|---|---|---|---|---|---|
| E1 | CM:PPI (100:0)-M15 | 0.16 ± 0.001 ᶜ | 10.3 + 1.0 ᵇ | 7.0 ± 1.9 ᶜ | 28 ± 5 ᵇᶜ | 52 ± 5 ᵇᶜ | − 3.1 ± 0.6 ᵈ | 19 ± 2 ᵃ |
| E2 | CM:PPI (85:15)-M15 | 0.13 ± 0.01 ᶜ | 11.5 + 1.2 ᵃ | 8.9 ± 1.8 ᵃᵇᶜ | 31 ± 2 ᵃ | 47 ± 5 ᶜ | − 0.9 ± 0.3 ᶜ | 14 ± 2 ᵇᶜ |
| E3 | CM:PPI (70:30)-M15 | 0.28 ± 0.01 ᵇᶜ | 7.6 + 1.2 ᶜ | 9.4 ± 2.2 ᵃ | 31 ± 2 ᵃᵇ | 48 ± 4 ᶜ | 1.0 ± 0.4 ᵇ | 14 ± 3 ᵇᶜ |
| E4 | CM:PPI (50:50)-M17 | 0.37 ± 0.02 ᵇᶜ | 5.5 + 0.7 ᵈ | 9.0 ± 1.3 ᵃᵇ | 26 ± 2 ᶜᵈ | 48 ± 3 ᶜ | 2.3 ± 0.4 ᵃ | 14 ± 2 ᵇᶜ |
| E5 | CM:PPI (50:50)-M16 | 0.35 ± 0.02 ᵇᶜ | 5.7 + 0.6 ᵈ | 7.1 ± 0.8 ᵇᶜ | 26 ± 3 ᶜᵈ | 51 ± 2 ᵇᶜ | 2.7 ± 0.6 ᵃ | 15 ± 2 ᵇ |
| E6 | CM:PPI (50:50)-M15 | 0.44 ± 0.02 ᵇ | 5.2 + 0.6 ᵈ | 9.0 ± 1.2 ᵃᵇ | 27 ± 2 ᶜᵈ | 54 ± 3 ᵃᵇ | 3.1 ± 1.0 ᵃ | 14 ± 3 ᵇᶜ |
| E7 | CM:PPI (30:70)-M17 | 1.28 ± 0.12 ᵃ | 3.5 + 0.1 ᵉ | 2.8 ± 0.6 ᵈ | 24 ± 2 ᵈ | 59 ± 3 ᵃ | 3.0 ± 0.9 ᵃ | 13 ± 1 ᵇᶜ |
| E8 | CM:PPI (30:70)-M16 | 1.45 ± 0.12 ᵃ | 3.6 + 0.2 ᵉ | 2.9 ± 0.4 ᵈ | 24 ± 2 ᶜᵈ | 59 ± 4 ᵃ | 2.6 ± 1.1 ᵃ | 12 ± 2 ᶜ |
| E9 | CM:PPI (30:70)-M15 | 1.32 ± 0.26 ᵃ | 3.6 + 0.4 ᵉ | 2.7 ± 0.8 ᵈ | 24 ± 2 ᶜᵈ | 58 ± 4 ᵃ | 3.2 ± 1.5 ᵃ | 13 ± 2 ᵇᶜ |

### Table 4 (clean). "Physical properties of cornmeal-pea protein isolate extrudates (50:50 and 30:70 blends) at varying screw speeds and barrel temperature profiles"

| | sample | Bulk Density (g/mL) | Diameter (mm) | Hardness (N) | Fracturability (mm) | L* | a* | b* |
|---|---|---|---|---|---|---|---|---|
| E10 | CM:PPI (50:50)-SS500 | 0.42 ± 0.02 ᶜ | 5.0 + 0.6 ᵃ | 7.7 ± 1.0 ᵃᵇ | 25 ± 3 ᵃᵇ | 51 ± 1 ᵈ | 3.0 ± 0.5 ᵃᵇ | 15 ± 1 ᵃᵇᶜ |
| E11 | CM:PPI (50:50)-SS450 | 0.39 ± 0.02 ᶜ | 5.1 + 0.4 ᵃ | 8.0 ± 1.3 ᵃᵇ | 26 ± 3 ᵃᵇ | 52 ± 2 ᶜᵈ | 3.1 ± 1.2 ᵃᵇ | 14 ± 3 ᵃᵇᶜᵈ |
| E6 | CM:PPI (50:50)-SS400 | 0.44 ± 0.02 ᶜ | 5.2 + 0.6 ᵃ | 9.0 ± 1.2 ᵃ | 27 ± 2 ᵃᵇ | 54 ± 3 ᵇᶜᵈ | 3.1 ± 1.0 ᵃᵇ | 14 ± 3 ᵃᵇᶜᵈ |
| E12 | CM:PPI (50:50)-T2 | 0.43 ± 0.02 ᶜ | 5.1 + 0.6 ᵃ | 7.5 ± 0.8 ᵇ | 27 ± 1 ᵃ | 53 ± 2 ᵇᶜᵈ | 3.5 ± 1.2 ᵃ | 15 ± 2 ᵃ |
| E13 | CM:PPI (50:50)-T3 | 0.43 ± 0.02 ᶜ | 5.3 + 0.7 ᵃ | 7.8 ± 1.4 ᵃᵇ | 26 ± 2 ᵃᵇ | 53 ± 1 ᵇᶜᵈ | 2.9 ± 1.0 ᵃᵇ | 14 ± 3 ᵃᵇᶜᵈᵉ |
| E14 | CM:PPI (50:50)-T4 | 0.44 ± 0.03 ᶜ | 5.1 + 0.6 ᵃ | 7.0 ± 1.0 ᵇ | 26 ± 3 ᵃᵇ | 52 ± 3 ᵇᶜᵈ | 3.3 ± 1.0 ᵃᵇ | 15 ± 2 ᵃᵇ |
| E15 | CM:PPI (30:70)-SS500 | 1.53 ± 0.16 ᵃ | 3.6 + 0.2 ᵇ | 2.5 ± 0.4 ᶜ | 24 ± 2 ᵇ | 56 ± 1 ᵃᵇᶜ | 3.7 ± 0.4 ᵃ | 13 ± 1 ᵃᵇᶜᵈᵉ |
| E16 | CM:PPI (30:70)-SS450 | 1.10 ± 0.06 ᵇ | 3.8 + 0.3 ᵇ | 2.7 ± 0.6 ᶜ | 24 ± 2 ᵃᵇ | 56 ± 4 ᵃᵇ | 3.0 ± 0.7 ᵃᵇ | 12 ± 1 ᵇᶜᵈᵉ |
| E9 | CM:PPI (30:70)-SS400 | 1.32 ± 0.26 ᵃᵇ | 3.6 + 0.4 ᵇ | 2.7 ± 0.8 ᶜ | 24 ± 2 ᵃᵇ | 58 ± 4 ᵃ | 3.2 ± 0.5 ᵃᵇ | 13 ± 2 ᵃᵇᶜᵈᵉ |
| E17 | CM:PPI (30:70)-T2 | 1.48 ± 0.10 ᵃ | 3.4 + 0.3 ᵇ | 2.2 ± 0.5 ᶜ | 24 ± 2 ᵃᵇ | 59 ± 3 ᵃ | 2.5 ± 0.9 ᵃᵇ | 12 ± 1 ᶜᵈᵉ |
| E18 | CM:PPI (30:70)-T3 | 1.40 ± 0.11 ᵃ | 3.7 + 0.5 ᵇ | 2.2 ± 0.7 ᶜ | 25 ± 1 ᵃᵇ | 60 ± 2 ᵃ | 2.4 ± 0.6 ᵃᵇ | 12 ± 1 ᵈᵉ |
| E19 | CM:PPI (30:70)-T4 | 1.27 ± 0.10 ᵃᵇ | 3.6 + 0.4 ᵇ | 2.5 ± 0.7 ᶜ | 24 ± 2 ᵃᵇ | 60 ± 4 ᵃ | 2.0 ± 0.8 ᵇ | 11 ± 2 ᵉ |

**The colour columns are the only browning-adjacent numbers in this paper.** They run the wrong way
for a browning reading: L\* *rises* from 52 (corn only) through 47-48 (15-30 % PPI) to 58-60 (70 %
PPI), because pea protein isolate is a pale powder that dilutes the corn's own pigments, and b\*
falls from 19 to 11-13 for the same reason. **L\* here is a composition effect, not a Maillard
extent**, and the paper never uses it as one. It cannot be set against the trunk's browning
response.

### Tables 1 and 3 — GARBLED, and what is legible

Table 1: "Process and material characteristics during extrusion of cornmeal-pea protein isolate feed
compositions at varying moisture contents." Table 3: the same for the screw-speed and
temperature-profile series. Both are wide tables whose columns are the samples (E1-E9 and
E10-E19), and the text layer has flattened them into a vertical ribbon of numbers with the row
labels stacked separately. The row labels *are* legible from Table 1, and they are: **Bulk Density
of Dry Feed (g/mL); Moisture Content of Dry Feed (%); Moisture Content of Feed (%); Moisture
Content of Extrudate (%); Flow Rate of Liquid Feed (mL/min)** — plus, from the ribbon, values that
must belong to torque, pressure, melt temperature, screw speed and SME.

**No cell from either table is assigned to a sample here.** Isolated values that are legible and
whose identity is fixed by the body text rather than by the ribbon:

| statement | value | where |
|---|---|---|
| SME range across the composition series at 15 % moisture | **540 to 684 kJ/kg** | Results 3.1, prose |
| torque range across the same series | **13.0 to 18.0 N·m** | Results 3.1, prose |
| barrel profile T1 (composition series, constant) | 40-60-80-100-120-140-145-145 °C | Methods 2.3 and Results 3.1 |
| barrel profile T4 | 40-60-80-110-140-155-160-160 °C | Results 3.3.2, prose |
| barrel profiles T2 and T3 | "up to 150 °C" and "up to 155 °C" (the eight zone set-points are **not** printed) | Results 3.1 |
| feed rate (constant) | 3.4 kg/h | Methods 2.3 |
| screw speed (composition series, constant) | 400 rpm | Methods 2.3 |
| direction of the SME and torque effect | increasing PPI **lowers** both, because starch is what resists flow | Results 3.1 |
| direction of the screw-speed effect | lower screw speed gives **higher** torque and **lower** SME | Results 3.1 |

### Every number printed in the running text

| quantity | value | where |
|---|---|---|
| pea protein isolate, protein content | 80 % | Methods 2.1 |
| cornmeal composition | 67 % carbohydrate, 13 % fat, 7 % protein, 0.4 % fibre | Methods 2.1 |
| **acrylamide in raw cornmeal** | **2.1 ± 0.2 µg/kg** | Results 3.2 |
| **acrylamide in raw pea protein isolate** | **12.2 ± 1.5 µg/kg** | Results 3.2 |
| **acrylamide, corn-only extrudate (E1)** | **19.8 ± 2.9 µg/kg** | Results 3.2 |
| **acrylamide, 15 % PPI extrudate (E2)** | **37.4 ± 1.3 µg/kg** (a significant rise, p < 0.05) | Results 3.2 |
| acrylamide at 30 / 50 / 70 % PPI | "did not cause additional increases"; **decreased** at 50 % and 70 % relative to 15 %. **The values are in Figure 1a only** | Results 3.2 |
| **asparagine in cornmeal** | **128 mg/kg** | Results 3.2, quoting the absent Table S2 |
| **asparagine in pea protein isolate** | **15 mg/kg** | Results 3.2, quoting Table S2 |
| **reducing sugar, pea protein isolate** | **0.011 %** | Results 3.2 |
| **total sugar, pea protein isolate** | **0.031 %** | Results 3.2 |
| **reducing sugar, cornmeal** | **0.4 %** | Results 3.2 |
| **total sugar, cornmeal** | **1.0 %** | Results 3.2 |
| total free amino acids, cornmeal vs PPI | "not significantly different" (no number printed) | Results 3.2 |
| EC benchmark for non-whole-grain cereals | 150 µg/kg (Reg. (EU) 2017/2158) | Results 3.2 |
| acrylamide LOD / LOQ | 0.08 / 0.26 µg/kg | Methods 2.6 |
| comparator: extruded breakfast cereal, corn semolina | < LOQ to 210 µg/kg at 96-243 Wh/kg SME, 11.7-20 % moisture, barrel 180 °C (Lipinski 2025) | Results 3.2 |
| comparator: corn grits + defatted press cakes | acrylamide below LOD (< 3.79 µg/kg) (Jozinović 2024) | Results 3.2 |
| volatiles semi-quantified in pea protein isolate | **42 compounds** | Results 3.3.1 |
| volatiles semi-quantified in cornmeal | **26 compounds** | Results 3.3.1 |
| pyrazines in the raw materials | **none detected, in either** | Results 3.3.1 |
| pyrazines in the corn-only extrudate | **9** | Results 3.3.1 |
| pyrazines identified across all extrudates | **15** (P1-P15, listed in the Fig. 2 caption) | Results 3.3.1 |
| 2-furfural maximum | at 50 % PPI, "approximately **five times** higher than those in the cornmeal-only extrudate" | Results 3.3.1 |
| dimethyl disulfide at 70 % PPI | "approximately doubled", **not significant** (p > 0.05) | Results 3.3.1 |
| pyrazine response to barrel temperature | going from T1 to T4 the listed pyrazines "**nearly doubled**"; T2 and T3 gave no significant increase | Results 3.3.2 |
| odour threshold, 3-ethyl-2,5-dimethylpyrazine | 0.4 µg/L in water (van Gemert 2011) | Results 3.3.1 |
| odour threshold, 2-ethyl-3-methylpyrazine | 0.4 µg/L in water | Results 3.3.1 |
| odour threshold, 2,5-dimethylpyrazine | 800 µg/L in water | Results 3.3.1 |
| odour threshold, 2,6-dimethylpyrazine | 200 µg/L in water | Results 3.3.1 |
| odour threshold, trimethylpyrazine | 400 µg/L in water | Results 3.3.1 |
| cited flavour-dilution factor for six PPI off-flavours | FD = 243, with (E,E)-2,4-decadienal OAV 26 431 (Li 2024b) | Results 3.3.1 |

**Everything else is figure-only or off disk.** Figure 1a and 1b are the acrylamide bar charts for
all nineteen extrudates — per house rule, not typed. Figure 2a and 2b are PCA biplots. Tables S1-S4
are not in the PDF, and **S3 and S4 hold every volatile number in the study**.

### The compound list in the Figure 2 caption (typed because it is the paper's own compound registry)

Strecker aldehydes: 2-methylpropanal, 3-methylbutanal, 2-methylbutanal, phenylacetaldehyde.
Lipid-derived aldehydes: butanal, pentanal, hexanal, 2-hexenal, heptanal, (E)-2-heptenal, octanal,
(E)-2-octenal, nonanal, (E,E)-2,4-decadienal. Other aldehydes: benzaldehyde, 2-furfural.
Diketones: **2,3-butanedione**, 2,3-pentanedione, 2,3-octanedione. Ketones: 2-butanone,
2-hexanone, 2-heptanone, 6-methyl-5-hepten-2-one, 2-octanone, 2,3-dimethyl-2-cyclopenten-1-one,
(E)-3-octen-2-one, (E,E)-3,5-octadien-2-one, 2-nonanone, (E,Z)-3,5-octadien-2-one, 2-decanone.
Alcohols: 1-penten-3-ol, 3-methylbutanol, 2-methylbutanol, 1-pentanol, 1-hexanol, 1-heptanol,
1-octen-3-ol, 1-octanol, 1-nonanol. Furans: 2-ethylfuran, 2-pentylfuran. Sulfur: dimethyl
disulfide. Acids: butanoic, 3-methylbutanoic, 2-methylbutanoic, pentanoic, hexanoic.
**Pyrazines (P1-P15): pyrazine, 2-methylpyrazine, 2,5(6)-dimethylpyrazine, 2-ethylpyrazine,
2,3-dimethylpyrazine, 2-ethyl-5-methylpyrazine, 2-ethyl-6-methylpyrazine, trimethylpyrazine,
2-ethyl-3-methylpyrazine, 2-ethenyl-6-methylpyrazine, 3-ethyl-2,5-dimethylpyrazine,
2-methyl-6-(1-propenyl)pyrazine, 2,3-diethyl-5(or 6)-methylpyrazine, 2-methyl-3,5-diethylpyrazine,
3-propyl-2,5-dimethylpyrazine.**

### Arithmetic (mine)

**1. The acrylamide in the corn-only extrudate is nine times what the feed brought in.**
19.8 µg/kg out of a feed carrying 2.1 µg/kg of cornmeal-derived acrylamide is a **9.4-fold rise
(mine)**, so ~17.7 µg/kg was made during extrusion. For E2 the feed is 85 % cornmeal + 15 % PPI, so
its inherited acrylamide is 0.85 × 2.1 + 0.15 × 12.2 = **3.6 µg/kg (mine)**, and the measured
37.4 µg/kg is a **10.3-fold rise**. So both extrudates make roughly ten times what they inherit,
and the *difference* between them — 37.4 against 19.8, nearly a doubling — is not explained by the
0.7 µg/kg difference in inherited acrylamide. That is the paper's own puzzle, and its own answer is
lipid-derived carbonyls from the pea protein, not asparagine.

**2. The asparagine arithmetic confirms the puzzle.** With cornmeal at 128 mg/kg and PPI at
15 mg/kg, the feed asparagine is 128 mg/kg at 0 % PPI and **111 mg/kg at 15 % PPI (mine)** — it
*falls* by 13 % while the acrylamide nearly doubles. At 50 % PPI it is 71.5 mg/kg and at 70 %
49.9 mg/kg (both mine). **Acrylamide and asparagine move in opposite directions over the first
step of the series and in the same direction thereafter**, which is why the authors reach for lipid
oxidation and then recommend model systems.

**3. Reducing sugar is not limiting at any composition.** Cornmeal 0.4 % reducing sugar against
asparagine 128 mg/kg = 0.0128 %: **sugar is in ~31-fold mass excess over asparagine in cornmeal
(mine)**, and in ~7-fold excess in PPI (0.011 % against 15 mg/kg). On a molar basis with glucose at
180.16 and asparagine at 132.12 the excesses are ~23-fold and ~5-fold (mine). So in every
formulation the acrylamide-forming condensation is asparagine-limited, which supports the authors'
reading and is the same regime the acrylamide lane's `k_asn_glc` operates in.

## 4. Kinetic numbers the repository can use

**There are none.** No rate constant, no barrier, no reference temperature. What the paper offers is
levels and orderings, and the table below is honest about that.

**Registry mapping (`data/keys/compounds.yml`).** Present and keyed among this paper's compounds:
`acrylamide`, `2_3_butanedione`, `2_methylbutanal`, `2_methylpropanal`, `3_methylbutanal`,
`phenylacetaldehyde`, `benzaldehyde`, `furfural`, `hexanal`, `heptanal`, `nonanal`, `e_2_octenal`,
`1_hexanol`, `1_octen_3_ol`, `2_pentylfuran`, `dimethyl_disulfide`, `methylpyrazine`,
`2_5_dimethylpyrazine`, `2_3_dimethylpyrazine`, `2_6_dimethylpyrazine`, `2_ethylpyrazine`,
`trimethylpyrazine`, `tetramethylpyrazine` (not found here), `pyrazines` (the family id). **Absent
and used by this paper:** 2-ethyl-3-methylpyrazine, 2-ethyl-5-methylpyrazine,
2-ethyl-6-methylpyrazine, 3-ethyl-2,5-dimethylpyrazine, 2-ethenyl-6-methylpyrazine,
2-methyl-6-(1-propenyl)pyrazine, 2,3-diethyl-5(or 6)-methylpyrazine, 2-methyl-3,5-diethylpyrazine,
3-propyl-2,5-dimethylpyrazine, bare pyrazine itself (the parent, which wave B18's `PZ` species
models and which the registry does **not** key), 2,3-pentanedione, 2,3-octanedione, 2-ethylfuran,
the C4-C6 acids, and every remaining alcohol and ketone.

Every row below shares: **cornmeal + pea protein isolate feed dried 24 h at 70 °C, extruded on a
Thermo Process 16 twin-screw at 3.4 kg/h through a 3 mm die, barrel set-points to 145 °C (T1) or
160 °C (T4), 15-17 % feed moisture, 400-500 rpm, no measured pH, no measured water activity, melt
temperature and residence time not legible on disk.**

| step | quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|---|
| — | acrylamide, raw cornmeal | 2.1 ± 0.2 | µg/kg | unextruded | — | Results 3.2 | **level_only** |
| — | acrylamide, raw pea protein isolate | 12.2 ± 1.5 | µg/kg | unextruded | — | Results 3.2 | **level_only** |
| Asn + sugar → acrylamide (net) | **acrylamide, corn-only extrudate E1** | **19.8 ± 2.9** | µg/kg | 100:0 CM:PPI, 15 % moisture, T1 (to 145 °C), 400 rpm | — | Results 3.2 | **level_only** — the one extrudate acrylamide value printed as a number |
| " | acrylamide, 15 % PPI extrudate E2 | 37.4 ± 1.3 | µg/kg | 85:15, 15 % moisture, T1, 400 rpm | — | Results 3.2 | **level_only** |
| " | acrylamide, E3 to E19 | — | µg/kg | 17 further formulations | — | **Figure 1a and 1b** | **figure_only** |
| " | net acrylamide formed during extrusion, E1 | 17.7 (**mine**, 19.8 − 2.1) | µg/kg | as E1 | — | derived (mine) | derived_assumption |
| " | net formed, E2 | 33.8 (**mine**, 37.4 − 0.85×2.1 − 0.15×12.2) | µg/kg | as E2 | — | derived (mine) | derived_assumption |
| precursor | **asparagine, cornmeal** | **128** | mg/kg | raw material | — | Results 3.2, quoting the **absent** Table S2 | **level_only** |
| precursor | **asparagine, pea protein isolate** | **15** | mg/kg | raw material | — | Results 3.2, quoting Table S2 | **level_only** |
| precursor | reducing sugar / total sugar, cornmeal | 0.4 / 1.0 | % (w/w, basis not stated) | raw material | — | Results 3.2 | **level_only** |
| precursor | reducing sugar / total sugar, pea protein isolate | 0.011 / 0.031 | % (w/w) | raw material | — | Results 3.2 | **level_only** |
| precursor ratio | sugar over asparagine, cornmeal | ~31 by mass, ~23 by mole (**mine**) | — | raw material | — | derived (mine) | within_study_ratio |
| acrylamide vs composition | ordering | rises 0 → 15 % PPI, flat 15 → 30 %, **falls** at 50 % and 70 % | — | 15 % moisture, T1, 400 rpm | — | Results 3.2 + Fig. 1a | **within_study_ratio** (the shape is prose; the values are figure-only) |
| acrylamide vs moisture | **no significant effect** (p > 0.05) over 15-17 % | — | — | — | — | Results 3.2 | measured null |
| acrylamide vs screw speed and barrel profile | **no significant effect** (p > 0.05) | — | — | 50:50 and 30:70 feeds | — | Results 3.2 | measured null |
| pyrazines, raw materials | **none detected**, in cornmeal or in pea protein isolate | — | unextruded | — | Results 3.3.1 | **measured null** |
| pyrazines, corn-only extrudate | 9 of the 15 present | count | E1 | — | Results 3.3.1 | measured (a count, not a level) |
| total pyrazine vs composition | maximum at **30 % PPI**; no further rise to 50 %; **falls** at 70 % because the sugar side is diluted | — | 15 % moisture, T1, 400 rpm | — | Results 3.3.1 | **within_study_ratio** (the levels are in the absent Table S3) |
| total pyrazine vs barrel temperature | **"nearly doubled"** from T1 (to 145 °C) to T4 (to 160 °C), for eight named pyrazines; T2 and T3 gave no significant rise | — | 50:50 and 30:70 | — | Results 3.3.2 | **within_study_ratio** |
| 2-furfural | maximum at 50 % PPI, ~**5x** the corn-only extrudate | — | — | — | Results 3.3.1 | within_study_ratio |
| 2,3-butanedione | "primarily formed through Maillard reaction during extrusion"; **not affected** by PPI level | — | — | — | Results 3.3.1 | measured null on the composition axis |
| Strecker aldehydes (2-methylpropanal, 2- and 3-methylbutanal) | rise to 30 % PPI, flat to 50 %, **fall** at 70 % | — | — | — | Results 3.3.1 | within_study_ratio |
| phenylacetaldehyde | rises on extrusion, **no** significant composition effect | — | — | — | Results 3.3.1 | measured null |
| 2-methylpropanal in raw materials | **not detected**; formed during extrusion of cornmeal | — | — | — | Results 3.3.1 | measured null |
| every volatile level | — | peak-area ratio to one of three internal standards | — | — | **Tables S3 and S4, absent from disk** | not on disk |
| melt temperature and residence time | — | °C, s | — | — | Tables 1 and 3, **garbled in the text layer** | not legible |
| L*, a*, b* of all 19 extrudates | see Tables 2 and 4 | — | — | — | Tables 2 and 4 | measured — but a **composition** effect, not a browning extent (section 3) |

### What is new against the 2016 paper, point by point

| axis | Göncüoğlu Taş & Gökmen 2016 (hazelnut) | this paper (2026, extrudates) |
|---|---|---|
| question | mechanism | product quality and food safety |
| method | multiresponse kinetic modelling, Athena Visual Studio, determinant criterion, 95 % HPD | ANOVA + Tukey + PCA |
| output | **26 rate constants at three temperatures** and (in the absent Table S4) 26 barriers at T_b = 160 °C | **no constant of any kind** |
| matrix | whole roasted hazelnut, 2.5-5 % moisture, >56 % oil | extruded corn + pea protein, 15-17 % feed moisture, ~9-10 % extrudate moisture |
| thermal treatment | static oven, 150-170 °C, 15-120 min, isothermal by assumption | twin-screw extruder, eight zones 40 → 145-160 °C, residence time not reported |
| responses | sucrose, glucose, fructose, total amino acids, 3-DG, 1-DG, 3,4-DG, glyoxal, methylglyoxal, dimethylglyoxal, HMF | acrylamide, 42 semi-quantified volatiles, bulk density, diameter, hardness, fracturability, L\*a\*b\* |
| dicarbonyls | seven, quantified by o-PDA/quinoxaline HPLC-MS with calibration curves | **none** (2,3-butanedione and 2,3-pentanedione appear only as headspace volatiles) |
| acrylamide | **not measured** | measured by isotope-dilution UPLC-MS/MS, LOQ 0.26 µg/kg |
| pyrazines | **not measured** | 15 identified, semi-quantified, none in the raw materials |
| institutions | Hacettepe alone | Hacettepe **and Reading** (Parker, Balagiannis, Ghawi) |
| use to this repository | a second laboratory's constants for six trunk steps | precursor levels, one acrylamide level, and orderings |

**The one methodological continuity worth recording** is that both papers stop at the same wall.
The 2016 paper says its barriers "could not be explained by the Arrhenius equation" and asks for a
wider temperature range; this one says explaining acrylamide against composition "is challenging"
and asks for model systems. Both are real-food studies whose authors conclude that the real food is
where the kinetics stops being identifiable. That is a datum about the corpus, and it is the reason
the trunk is fitted on Martins' aqueous pot and not on a food.

## 5. Flags

1. **The two tables that would make this paper usable are the two that are garbled.** Tables 1 and 3
   carry the **melt temperature** ("the actual temperature of the material inside the barrel") and
   the process variables; the text layer flattens them into an unassignable ribbon. Without a melt
   temperature and without a residence time — which is **never** reported, in any table or in the
   text — no thermal history exists, so no acrylamide prediction can be made for any of the nineteen
   extrudates, however good the measurement is. **To use this paper the two tables must be read from
   the page images and the residence time requested from the authors.**
2. **The acrylamide dataset is a figure.** Only E1 (19.8 ± 2.9) and E2 (37.4 ± 1.3 µg/kg) are
   printed as numbers; the other seventeen are bars in Figure 1a and 1b. Per house rule they are not
   typed. Requesting the underlying values is item two on the list below.
3. **The volatiles are SEMI-QUANTIFIED against three internal standards and are not concentrations.**
   3-Furaldehyde for furans, isopropylpyrazine for pyrazines, 3-methyl-2-butanone for the Strecker
   aldehydes and "the rest of the compounds"; no response factors were measured. So a "five times
   higher" for 2-furfural or a "nearly doubled" for pyrazines is a ratio of area ratios within one
   method, which is fine as an ordering and inadmissible as a level. **No volatile number from this
   paper may enter a benchmark against a modelled µg/kg.** And they are all in the absent
   Tables S3/S4 in any case.
4. **The feed was dried 24 h at 70 °C before extrusion, and the raw materials already contain
   acrylamide.** Cornmeal 2.1 and pea protein isolate 12.2 µg/kg, which the authors attribute to the
   suppliers' own drying. So each run starts with a non-zero acrylamide and a partly reacted
   precursor pool of unknown extent — an initial condition the model would have to be given, and one
   that the 24 h at 70 °C in the authors' own tray oven adds to.
5. **The lipid route is not in the repository and this paper says it is what drives the effect.**
   The authors' explanation for the acrylamide rise at 15 % PPI, having ruled out asparagine and
   sugar by their own numbers, is "lipid oxidation products with reactive carbonyl groups from pea
   protein isolate, which can promote Maillard reaction". The engine's acrylamide lane has no lipid
   carbonyls at all. Any attempt to score this paper would be scoring a mechanism the model does not
   contain.
6. **The "total sugar" and "reducing sugar" percentages have no stated basis.** 0.4 % and 1.0 % for
   cornmeal, 0.011 % and 0.031 % for PPI — w/w presumably, but wet or dry basis is not said, and the
   sugars were measured by HPAEC-PAD with calibration for glucose, fructose, sucrose and maltose
   individually, so a per-sugar breakdown exists in the absent Table S2 and only the totals reach
   the body text.
7. **Colour is a dilution effect here.** L\* rises and b\* falls as pea protein replaces corn, and
   the paper does not read L\* as browning. Do not treat Table 2 or Table 4's L\* as a browning
   response; the trunk's browning hold-out gets nothing from this paper.
8. **The composition and process series are confounded by which feeds were used.** The screw-speed
   and temperature series were run **only on 50:50 and 30:70** feeds, chosen "because they more
   clearly reflect protein-driven effects" — which are exactly the two compositions where the
   acrylamide had already fallen back and where the authors ascribe the null result to "the limited
   asparagine content". So the finding that temperature and screw speed do not affect acrylamide is
   established on the two lowest-asparagine feeds only, and does not generalise to E1 or E2.
9. **What this paper does not contain**: any rate constant; any activation energy; any reference
   temperature; any residence time; any pH; any water activity; any dicarbonyl measurement; any
   sugar or amino acid measurement of the extrudates (only of the raw materials); any melanoidin or
   browning measurement; any replicate count for the extrusion runs themselves (the analytical
   replicates are stated, n = 3 or 5, but each formulation appears to have been extruded once).
10. **What to request from the authors** (data availability says "Data will be made available on
    request"): (i) **Supplementary Tables S2, S3 and S4** — the raw-material amino acid and sugar
    profile and the two volatile tables, which is the whole aroma dataset; (ii) the numeric
    acrylamide values behind Figure 1a and 1b for E1-E19; (iii) **Tables 1 and 3 in a legible form,
    especially the melt temperature column**; (iv) **the residence time distribution**, without
    which no thermal history exists; (v) the zone-by-zone set-points of profiles T2 and T3, only T1
    and T4 being printed.
11. **Registry gaps against `data/keys/compounds.yml`.** The registry keys `acrylamide` and most of
    the Strecker aldehydes and lipid aldehydes this paper reports, but **not the parent pyrazine
    itself** — which wave B18 models as the species `PZ` and which this paper lists as P1 — nor nine
    of the fifteen pyrazines it identifies, nor 2,3-pentanedione, 2,3-octanedione or 2-ethylfuran.
    If the pyrazine lane is ever to be scored against a real extruded food, the parent pyrazine
    needs a registry id.
