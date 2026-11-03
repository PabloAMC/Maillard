# Palazoğlu & Gökmen 2008 — EXTRACTION (file stem `palazoglu2008`)
DOI: 10.1021/jf073046l (printed p. 6162 footer). Crossref 2026-11-03: "Reduction of Acrylamide Level in French Fries
by Employing a Temperature Program during Frying", J. Agric. Food Chem. **56**(15):6162-6166, authors Palazoğlu,
Gökmen, issued 2008-07-15. **Matches the brief and the PDF.** Not to be confused with Palazoğlu & Gökmen 2008 J Food
Sci 73:E109 (10.1111/j.1750-3841.2008.00676.x; not on disk), which Chan 2020 also cites.

### Single Agria potato strips (8.5 × 8.5 × 70 mm; reducing sugar 1.4 ± 0.2 g/kg, free Asn 7.2 ± 0.6 g/kg fresh) fried in sunflower oil with an oil-temperature step (170 or 190 °C for 1-2 min, then 150 °C) for 3 or 4 min total; acrylamide by GC-MS (bromination, ¹³C₃-acrylamide IS). Plus a 170 °C time series at three sugar levels (0.6 / 1.4 / 2.4 g/kg). Measured values are compared with a 3-D finite-difference heat/mass model using cited Gökmen & Şenyuva 2006 fructose-Asn constants. **No measured product T, no moisture, no kinetics fitted here.**

**Read 2026-11-03 from `data/articles/palazoglu2008.pdf`** (5 pp., ACS, scanned-quality figures). Prose from
`pdftotext -layout`. Fig. 1 and Fig. 2 (p. 2 = 6163) read from 300 dpi crops (scratchpad
`ingest/palazoglu2008_work/fig1hi-2.png`, `fig2hi-2.png`); read error ±3 ng/g (Fig. 2a), ±5 (Fig. 2b), ±20 (Fig. 1).

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Food Engineering, **Hacettepe University** (Gökmen, corresponding) and University of Mersin
(Palazoğlu). **Existing id `hacettepe_gokmen`** (`data/keys/laboratories.yml` l. 60-62, string "Hacettepe (Gokmen)";
its evidence line names this paper). Same lab and same strip geometry / fryer: `gokmen2006_extraction.md`; also
`gokmen2012_extraction.md`. Cited as refs. 15 / 25 by `chan2020_extraction.md`.

## 1. What this paper answers for the model

| model need | what Palazoğlu 2008 gives | tag |
|---|---|---|
| precursor charge of the fried potato | reducing sugar (Glc + Fru) 1.4 ± 0.2 g/kg, free Asn 7.2 ± 0.6 g/kg (fresh, HPLC) | [M] |
| acrylamide vs time at constant oil T | 170 °C, 3 / 4 / 6 min at 0.6 / 1.4 / 2.4 g/kg sugar (Fig. 1; 1.4 has 3 and 4 min only) | [F] |
| acrylamide after temperature-programmed frying | 8 schedules, 3 and 4 min total, duplicates with error bars (Fig. 2) | [F] |
| measured product T or moisture | **none** (moisture "on the basis of our previous experience", ref. 17 = Gökmen 2006) | [!] |
| acrylamide kinetics | Eqs. 7-8, **cited** from Gökmen & Şenyuva 2006 Food Addit Contam 23:348 (not on disk): k′1 = exp(15.438 − 11284/T), k2 = exp(7.4728 − 5950/T), units not stated | [C] |
| Ea of formation / elimination | [Z] 93.8 and 49.5 kJ/mol from Eqs. 7-8 (fructose-Asn model system, 120-200 °C) | [Z] from [C] |

**Answers.**
1. **Does it test Chan 2020's step-frying predictions? No** (confidence ~95 %). Chan's schedule is 170 °C for 3 or
   4 min then 150 °C to 9 min (predicted 178 / 197 µg/kg vs 418 measured for 170 °C × 9 min). Palazoğlu's longest
   run is 4 min total with the step at 1 or 2 min. No schedule overlaps. The **direction** agrees (step reduces
   acrylamide by 35-63 % here vs Chan's simulated ~53-57 %), but that is not a quantitative test, and the tubers
   (1.4 g/kg) are not the 2006 tubers Chan simulated (Chan nevertheless took its sugar charge from this paper).
2. **Do the measured step-frying values carry information for the lane? Yes, a little** (confidence ~70 %): 8
   whole-strip values with a measured precursor charge, SIDA-quantified, one strip per fry so oil T is steady. But
   3-min values (70-116 ng/g) sit near the stated LOQ of 50 ng/g, no product T is measured, and the model needs a
   heat-transfer path to use them (the lane has none for frying; cf. `knol2009_*`).
3. **Is the model–experiment agreement a validation? Weak** (confidence ~75 %). The heat/mass parameters are from
   Yıldız et al. 2007 (ref. 25) for the same strip, the kinetics from a model system, and a hard 120 °C switch; no
   fitted constant is printed, so the agreement (within ±20 %, Fig. 2) is an out-of-sample comparison only if
   nothing was tuned — not stated. The 2.4 g/kg curve under-predicts 4 min by ~25 % and the 0.6 g/kg curve sits
   ~12 % above 6 min [F].
4. **The cited constants (Eqs. 7-8) matter for CF-039**: a same-lab apparent Ea of formation of 93.8 kJ/mol [Z],
   between Gökmen 2012's silica 83.94 and the 159-168 kJ/mol of Claeys / De Vleeschouwer. Medium of the cited
   experiment not given here [!] (the source is not on disk).

## 2. Conditions
- Potatoes cv. **Agria**, local market; strips 8.5 × 8.5 mm, 70 mm long [M]; sugar and Asn as above (fresh basis,
  ± presumably SD of duplicates) [M]. How the 0.6 and 2.4 g/kg lots of Fig. 1 were obtained is **not described** [!].
- Sunflower oil. First stage (170 / 190 °C) in a 5 L bath (Precisterm), then a quick transfer to a 1 L beaker on an
  IKA hot plate at 150 °C. **One strip per fry**, fully immersed on a steel wire assembly; oil T monitored, "no
  changes" observed [Q]. Duplicate fries [M].
- Schedules [M]: 3 min total: 170 × 3; 170 × 2 + 150 × 1; 170 × 1 + 150 × 2; 190 × 1 + 150 × 2. 4 min total: 170 × 4;
  170 × 2 + 150 × 2; 170 × 1 + 150 × 3; 190 × 1 + 150 × 3. (The abstract lists only the 4-min set.)
- No pH, a_w, moisture, oil uptake, colour, product T.

## 3. Analytical methods
Acrylamide: GC-MS (Agilent 5973, InnoWax), brominated derivative, SIM; **¹³C₃-acrylamide internal standard**
(500 ng added to 2 g homogenate, methanol extraction) → isotope dilution. Linear 0.01-1.0 µg/mL; **LOD 15, LOQ
50 ng/g** in French fries; CV ≤ 10 % (three repeats) [M]. ng/g of homogenised fried strip (as-is basis) [Z, ~85 %].
Sugars: HPLC-RI (Shodex SH-1011), LOD 0.01 g/100 g glucose. Asn: FMOC derivative, HPLC-FLD, LOD 0.1 mg/100 g [M].
Error bars in Fig. 2: from duplicate fries, statistic not stated [!].

## 4. Data
### 4a. Fig. 2, acrylamide (ng/g) after programmed frying, sugar 1.4 g/kg [F]
| schedule | model | experimental ± bar |
|---|---|---|
| 170 °C × 3 min | 125 | 116 ± 3.5 |
| 170 × 2 + 150 × 1 | 87.5 | 93 ± 11 |
| 170 × 1 + 150 × 2 | 56.5 | 70 ± 7 |
| 190 × 1 + 150 × 2 | 86 | 102.5 ± 6.5 |
| 170 °C × 4 min | 323 | 313 ± 18 |
| 170 × 2 + 150 × 2 | 180 | 158 ± 3 |
| 170 × 1 + 150 × 3 | 134 | 116 ± 4 |
| 190 × 1 + 150 × 3 | 186 | 205 ± 10 |

Reductions vs the 170 °C control [Z]: 4 min 50 / 63 / 35 % (exp.), 44 / 59 / 43 % (model); 3 min 20 / 40 / 12 % (exp.).
Text: "reduced by 58 %" for 170 × 1 + 150 × 3 [M] — equals the **model** reduction (59 %), not the measured one
(63 %) [Z]; abstract: "reduced by half" [Q].

### 4b. Fig. 1, acrylamide (ng/g) vs time at 170 °C, three sugar levels; markers measured, lines model [F]
| reducing sugar | 3 min | 4 min | 6 min | model at 6 min |
|---|---|---|---|---|
| 0.6 g/kg | 45 | 121 | 444 (bar ±60) | ~505 |
| 1.4 g/kg | 112 | 314 | — | ~1190 |
| 2.4 g/kg | 229 | 723 | 2045 (bar ±140) | ~2000 |
The 1.4 g/kg markers equal Fig. 2's 170 °C controls (112 vs 116; 314 vs 313) [Z]: same data, not independent.

### 4c. Cited kinetic constants (Eqs. 7-8, from ref. 23) [C], derived values [Z]
Scheme A + B → C (acrylamide) → D, reduced to A → C → D with k′1 = k1·C_B (zero order in Asn); first order in
fructose; 120-200 °C; C_A0 = 5 µmol in the source. Ea′1 = 11284 R = **93.8 kJ/mol**, Ea2 = 5950 R = **49.5 kJ/mol**
[Z]. k′1 = 1.3e-5 / 4.4e-5 / 1.3e-4 and k2 = 1.4e-3 / 2.6e-3 / 4.6e-3 at 150 / 170 / 190 °C (units not stated;
min⁻¹ or s⁻¹) [Z]. In the frying model, reaction runs only in nodes ≥ 120 °C; no transport of sugar or acrylamide.

## 5. Flags and internal consistency
1. [!] "58 %" is the model's reduction, not the measured 63 % (sec. 4a).
2. [!] Units of k′1 and k2 not stated; constants are from a fructose-Asn model system applied to potato with Glc +
   Fru treated as fructose [Q]; no Asn dependence.
3. [!] 3-min values 56-125 ng/g lie within ~1-2.5 × LOQ (50 ng/g).
4. [!] Fig. 1 sugar levels (0.6, 2.4 g/kg) have no described source; the Agria lot is 1.4.
5. [!] Moisture of the programmed fries asserted from Gökmen 2006, not measured; no product T measured in this paper.
6. Same-lab cross-check [Z]: 170 °C × 4 min here 313 ng/g vs Gökmen 2006 interpolated ~170 ng/g (six strips per
   fry, LC-MS external standard, unknown tubers) — a ×1.8 difference that tuber sugar or analytics could explain.

## 6. Existing records in the repo that cite this paper
`grep -rIn` (2026-11-03) of `data/ docs/ src/ scripts/ tasks/` for `jf073046l`, `palazoglu`, "Palazoglu",
"Palazoğlu", "temperature program": hits only in `chan2020_extraction.md` (refs. 15 / 25, correctly identified),
`data/keys/laboratories.yml` l. 62 (evidence line: correct, Gökmen corresponding, Hacettepe), and indirectly the
`chan2020_figs2_4_gokmen2006_fries_170c_replot` entry (`acrylamide.yml` l. 1226), whose `system` "reducing sugar
6.8 mol/m3 (as computed by Chan from ref. 15)" derives from this paper's 1.4 g/kg: 6.8 mol/m³ × 0.180 kg/mol =
1.22 kg/m³ ≈ 1.1-1.2 g/kg at ρ ~1.05-1.1 [Z] — consistent with 1.4 ± 0.2 at the low edge, but it is **this paper's
Agria lot, not the 2006 tubers** (see `gokmen2006_extraction.md` sec. 6). The cited Gökmen & Şenyuva 2006 kinetics
paper (Food Addit Contam 23:348) appears only as a deep-research bibliography line
(`data/research_corpus/Deep_Research_Protein_Matrix_Kinetics_Report.md` l. 186) and in two furfural notes
(`trunk.yml` l. 1815, `agcam2022`, `k5a_hmf_synthesis`) that concern a different finding; no constant from it is in
the lane (`grep 11284\|15.438\|5950` in `src/`: none).

## 7. Proposed ledger entries (acrylamide.yml; NOT written)
```yaml
  # ---- palazoglu2008 (ingested 2026-11-03) ----
  - id: palazoglu2008_fig2_temperature_programmed_frying
    dossier: palazoglu2008_extraction.md
    paper: "Palazoglu & Gokmen 2008 J Agric Food Chem 56:6162 (10.1021/jf073046l)"
    lab: "Hacettepe (Gokmen)"
    system: "potato strips cv. Agria 8.5 x 8.5 x 70 mm; reducing sugar 1.4 +- 0.2 g/kg, free Asn 7.2 +- 0.6 g/kg fresh"
    medium: "deep-fat frying in sunflower oil, one strip per fry; 170 or 190 C stage (5 L bath) then 150 C (1 L beaker), 3 or 4 min total"
    temperature_c: [150, 170, 190]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [3, 4]
    species: [acrylamide]
    kind: endpoint
    location: "Fig. 2a/2b, 8 schedules, experimental (duplicates, bars) and model bars (dossier sec. 4a)"
    status: digitised_in_dossier
    quant_class: sida
    simulable: false
    condition_quote: "the potato strip was quickly transferred into a 1 L oil bath"
    flags: [matrix_food, low_moisture_no_aw_term, heating_up_included_in_time, near_loq]
    note: "Measured precursor charge plus an oil-T step: a frying hold-out for a measured-T(t) path once one exists; no product T or moisture measured. Does not overlap Chan 2020's 9-min step schedules. Text '58 %' is the model's reduction (measured 63 %)."
  - id: palazoglu2008_fig1_acrylamide_sugar_level_170c
    dossier: palazoglu2008_extraction.md
    paper: "Palazoglu & Gokmen 2008 J Agric Food Chem 56:6162 (10.1021/jf073046l)"
    lab: "Hacettepe (Gokmen)"
    system: "potato strips 8.5 x 8.5 x 70 mm at three reducing-sugar levels 0.6 / 1.4 / 2.4 g/kg (source of the 0.6 and 2.4 lots not described)"
    medium: "deep-fat frying in sunflower oil, 170 C, one strip per fry"
    temperature_c: [170]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [3, 6]
    species: [acrylamide]
    kind: time_series
    location: "Fig. 1 markers (3/4/6 min; 1.4 g/kg at 3 and 4 min only = Fig. 2 controls) and model lines (dossier sec. 4b)"
    status: digitised_in_dossier
    quant_class: sida
    simulable: false
    condition_quote: "Effect of initial total reducing sugar content and frying time on acrylamide level"
    flags: [matrix_food, low_moisture_no_aw_term, heating_up_included_in_time]
    note: "Sugar-dose response of fried-potato acrylamide at fixed oil T; Asn given only for the 1.4 g/kg lot. 1.4 g/kg points duplicate Fig. 2 controls."
  - id: palazoglu2008_eq7_8_cited_fructose_asn_constants
    dossier: palazoglu2008_extraction.md
    paper: "Palazoglu & Gokmen 2008 J Agric Food Chem 56:6162, citing Gokmen & Senyuva 2006 Food Addit Contam 23:348"
    lab: "Hacettepe (Gokmen)"
    system: "fructose + asparagine model system (C_A0 5 umol), sugar -> acrylamide -> D, zero order in Asn"
    medium: "not stated in this paper (source not on disk)"
    temperature_c: [120, 200]
    ph: {initial: null, measured_during: false, final: null}
    time_min: []
    species: [acrylamide]
    kind: rate_constant
    location: "Eqs. 7-8 p. 6164: ln k'1 = 15.438 - 11284/T, ln k2 = 7.4728 - 5950/T (Ea 93.8 / 49.5 kJ/mol derived, dossier sec. 4c)"
    status: text_only
    quant_class: unknown
    simulable: false
    condition_quote: "Because acrylamide formation was first order with respect to fructose and zeroth order"
    flags: [secondary_source, units_not_stated]
    note: "Same-lab apparent formation Ea 93.8 kJ/mol (CF-039 evidence). Re-key to Gokmen & Senyuva 2006 if downloaded; medium and units must be read there."
```

## 8. Supplementary material
None referenced. The cited kinetic source (Gökmen & Şenyuva 2006, Food Addit Contam 23:348, 10.1080/02652030500482355 per Crossref 2026-11-03; the
DOI is not printed in this paper; `Gemini_deep_research_A.md` l. 136 mislabels it as Knol 2005) and the heat/mass parameters (Yıldız, Palazoğlu & Erdoğdu 2007 J Food Eng 79:11)
are outside this paper; neither is on disk.
