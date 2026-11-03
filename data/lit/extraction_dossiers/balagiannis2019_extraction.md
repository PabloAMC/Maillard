# Balagiannis, Mottram, Higley, Smith, Wedzicha & Parker 2019 — EXTRACTION (file stem `balagiannis2019`)
DOI: 10.1016/j.foodchem.2019.01.075 (printed p. 236). Food Chem. **284**:236-244; received 20 Aug 2018, accepted 11 Jan
2019, online 19 Jan 2019. Crossref 2026-11-03: "Kinetic modelling of acrylamide formation during the finish-frying of
french fries with variable maltose content", Food Chemistry 284, 236-244, issued 2019-06, authors Balagiannis,
Mottram, Higley, Smith, Wedzicha, Parker. **Matches the PDF and the brief.** No naming caveat; the other
`balagiannis2009/2010/2015` dossiers are different papers by the same first author.
**Companion:** `parker2012b_extraction.md` (Parker et al. 2012 JAFC 60:9321): same team, same pilot line; its
glucose- and fructose-dip data are refitted here together with the new maltose arms.

### Ranger Russet shoestring fries (2013 crop), pilot line as Parker 2012 but with a **maltose** dip (0, 1, 2, 3, 4, 14 % w/v), finish-fried in canola oil at 165 °C (6 arms × 10 times, 0.5-5 min) plus a second lot (2 % maltose) at 175 and 185 °C. Maltose, glucose, acrylamide (Fig. 1) and fructose, sucrose, free amino acids, Asn (Fig. S1, **SI not on disk**), defatted dry weight. Two multiresponse models fitted jointly to these and the Parker 2012 data: **Model 1** (generic pathway only; maltose via Int1mal with glucose "peel-off") and **Model 2** (+ direct sugar–Asn glycoconjugate routes kg, kf, km). Constants printed **at 165 °C only; no Ea, no T treatment described** (heat and mass transfer deferred to the Balagiannis 2016 3-D model, not on disk).

**Read 2026-11-03 from `data/articles/balagiannis2019.pdf`** (9 pp., Elsevier born-digital). Prose from
`pdftotext -layout`; **pp. 238-242 rendered at 150 dpi** (`ingest/balagiannis2019_work/p-3.png` … `p-7.png`);
**Table 1 and Fig. 4 from 200 dpi crops** (`t1-4.png`, `f4-7.png`), **Fig. 2 parameter table read by eye on the
150 dpi render** (`p-5.png`), Fig. 3 slopes from a 110 dpi crop (`f3-6.png`, printed on the figure). Fig. 1 digitised
at t = 0 and t = 300 s: read error ±2 µmol/kg (acrylamide, 0-75 axis), ±3 mmol/kg (maltose, 0-180), ±1.5 mmol/kg
(glucose, 0-50). Derived numbers [Z] by hand.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Department of Food and Nutritional Sciences, **University of Reading**; ConAgra Foods (Higley, Smith;
both since moved, per footnotes) and Leeds (Wedzicha). Corresponding author **J. K. Parker**. Existing lab id
**`reading_mottram`** ("Reading (Mottram / Parker / Elmore)"). Same-lab dossiers: `parker2012b` (direct predecessor),
`parker2013`, `parker2014`, `balagiannis2009`, `balagiannis2010`, `balagiannis2015`, `goncuoglu2026`, `kocadagli2021`.

## 1. What this paper answers for the model

| model need | what Balagiannis 2019 gives | tag |
|---|---|---|
| acrylamide formation constants incl. maltose | Fig. 2 table: Model 1 (8 parameters) and Model 2 (11), 95 % CI, **at 165 °C** | [M] |
| direct sugar–Asn (Schiff base / glycoconjugate) rate | Model 2: kg (Glc), kf (Fru) ~5-6 × 10⁻⁷ kg mmol⁻¹ s⁻¹; km (maltose) at lower bound 0 | [M] |
| Ea / T dependence | **none printed**; 175/185 °C data used in the fit but how T enters is not stated | [!] |
| maltose, glucose, acrylamide vs time | Fig. 1: 8 arms × 11 points (0-300 s); Table 1 endpoints | [F]/[M] |
| fructose, sucrose, FAA, Asn vs time (maltose arms) | **Fig. S1 only** (SI, absent); sucrose "did not alter" [Q] | [!] |
| fit quality | Fig. 3 regression-through-origin slopes and R² per response | [M] |
| surface / core T, moisture | **not measured or shown here**; "described by Parker et al. (2012)" / Balagiannis 2016 | [C] |
| colour | Table 1 Agtron after 5 min | [M] |
| pH, a_w | not measured | [!] |

**Answers.**
1. **Scheme** (Fig. 2, Kinetic Mechanism 1; confidence 98 %): glucose + FAA →(k1) (1 + p) Int2 — **one step**, the
   Parker 2012 Int1 is gone; fructose →(k6) (1 + p) Int2; glucose →(k8) fructose (**no k9** reverse step);
   maltose + FAA →(k1mal) Int1mal →(k2mal) glucose + (1 + p) Int2 (the Mundt & Wedzicha 2005 "peel-off" [C]);
   Int2 + FAA →(k3) Maillard products, Int2 + Asn →(k3) acrylamide with yield F_Asn. p = extra moles of Int2 per mole
   of sugar backbone (fragmentation). Asn is now an explicit species (Parker 2012 used R_Asn = 0.37). Mechanism 2 adds
   Glc + Asn →(kg), Fru + Asn →(kf), Mal + Asn →(km) → Int → acrylamide (fast).
2. **Constants** (Fig. 2, 165 °C; "× 10⁴" columns): see sec. 4a. Model 1 k6 = 0.0165 s⁻¹ and k8 = 6.8 × 10⁻³ s⁻¹ agree
   with Parker 2012 Table 2 (0.018, 7.0 × 10⁻³) within CIs; **k1 is ×3.5 lower** (2.2 vs 7.6 × 10⁻⁴ kg mmol⁻¹ s⁻¹)
   because the scheme and the transfer model changed (confidence ~85 %); F_Asn halves (3.2 vs 5.9 × 10⁻³) and is
   defined on a different rate (k3[Int2][Asn] vs Parker's Int2 formation × R_Asn). **Constants from the two papers must
   not be mixed.**
3. **T range / Ea: not given** (confidence 95 %). The text has no "activation", "kJ" or "Arrhenius"; the constants are
   labelled "correspond to a temperature of 165 °C" [V]. Most likely the 2016 3-D model with the Parker 2012 fixed
   Ea = 100 kJ/mol (~70 %), but not shown.
4. **Maltose contributes < 10 %** of acrylamide at equimolar sugar (Fig. 4 simulation, 40 mmol/kg each, 165 °C: maltose
   ~8 µmol/kg vs glucose ~33 and fructose ~34 at 300 s [F]; < 5 % before 200 s [M]). Model 2 returns km = 0: no direct
   maltose–Asn route detectable [Q]. These are model outputs, not measurements.
5. **Measured maltose effect is saturating**: Δacrylamide at 5 min 4.5 (0 %), 10.9 (1 %), 20.1 (2 %), then 17.5-20.8
   µmol/kg for 3-14 % dips (Table 1) while maltose loss rises 10.7 → 52.4 mmol/kg [M]. Molar yield per added-sugar
   lost, net of the M0 background: maltose 0.15 % (2 % dip) → 0.03 % (14 %) vs glucose 0.20 % and fructose 0.34 %
   (2 % dips, Parker data, gross) [Z]. The model does **not** reproduce the saturation (sec. 5, flag 2).
6. **Does the lane use any of this? No** (confidence ~97 %): no hit for the DOI, Balagiannis, maltose, F_Asn, k1mal,
   kg/kf/km or fries in `src/kinetic_core/` (acrylamide modules, `parameters_acrylamide.py`). The lane has no maltose
   species and no frying geometry. Transferable as constants: **no** (same reasons as Parker 2012, plus no Ea). Usable
   as a directional claim (maltose ≪ glucose ≈ fructose for acrylamide per mole; colour per sugar lost also lower) and,
   with Fig. S1, a fry hold-out.

## 2. Conditions

- **Potato** [M]: Ranger Russet, 2013 crop, stored 7-9 °C. Strips 7.5 × 7.5 mm, 7.5-12.5 cm.
- **Pre-processing** [M]: wash, peel, preheat 54 °C 30 min, cut; **blanch 80 °C, 10 min** (Parker 2012: 74 °C, though
  "as described by Parker" [!]); dip 30 s, 70-75 °C, maltose 0, 1, 2, 3, 4, 14 % w/v (14 % to reach ~1.4 % wet-weight
  maltose, like blanched sweet potato [Q]); dry 60 °C, 10 min (removes the ~10 % water uptake); par-fry 190 °C 50 s
  canola; pack, freeze −20 °C.
- **Finish-fry** [M]: ~220 g per batch, canola oil, 165 °C at 0.5-5.0 min (10 times); lot 2 (2 % maltose) at 175 and
  185 °C. Quench in liquid N2, store −20 °C. Fryer not restated (Parker 2012: 20.5 kW Frymaster, ~30 L).
- **Two potato lots** [M]: lot 1 raw fries Glc 12.7, total FAA 131 mmol/kg dry; lot 2 Glc 6.1, FAA 238.
- Replicates, pH, a_w: not stated / not measured [!].

## 3. Analytical methods

As Parker 2012 [C]: acrylamide by SPE + LC-MS/MS with **13C3-acrylamide IS (sida)**; sugars (Glc, Fru, maltose,
sucrose) HPLC with external calibration; 18 free amino acids OPA/Fmoc HPLC; moisture oven 125 °C; fat as FAMEs (here
"by GC–MS"; Parker 2012 says capillary GC vs C11:0 IS); Agtron E-30 colour. All on defatted dry weight. Modelling in
Athena Visual Studio; ± are 95 % confidence intervals.

## 4. Data

### 4a. Fig. 2 parameter table (p. 240), 165 °C, printed "× 10⁴" [M]; lane units [Z]

| parameter (unit) | Model 1 est. ×10⁴ | 95 % CI ×10⁴ | Model 2 est. ×10⁴ | 95 % CI ×10⁴ | Model 1 in lane units [Z] |
|---|---|---|---|---|---|
| k1 (mmol⁻¹ kg s⁻¹) | 2.2 | ±0.03 (1 %) | 2.1 | ±0.3 (16 %) | 0.0132 kg mmol⁻¹ min⁻¹ |
| k3 (mmol⁻¹ kg s⁻¹) | 72 | ±25 (35 %) | 220 | ±62 (28 %) | 0.432 kg mmol⁻¹ min⁻¹ |
| k6 (s⁻¹) | 165 | ±9 (5 %) | 161 | ±13 (8 %) | 0.99 min⁻¹ |
| k8 (s⁻¹) | 68 | ±8 (12 %) | 73 | ±15 (20 %) | 0.41 min⁻¹ |
| p (—) | 8550 | ±2510 (29 %) | 6640 | ±2713 (41 %) | 0.855 |
| F_Asn (—) | 32 | ±4 (14 %) | 18 | ±7 (38 %) | 3.2 × 10⁻³ |
| k1mal (mmol⁻¹ kg s⁻¹) | 2.3 | ±0.1 (5 %) | 2.3 | ±0.1 (5 %) | 0.0138 kg mmol⁻¹ min⁻¹ |
| k2mal (s⁻¹) | 162 | ±25 (15 %) | 245 | ±60 (25 %) | 0.97 min⁻¹ |
| kg (mmol⁻¹ kg s⁻¹) | — | — | 0.006 | ±0.003 (45 %) | (M2) 3.6 × 10⁻⁵ kg mmol⁻¹ min⁻¹ |
| kf (mmol⁻¹ kg s⁻¹) | — | — | 0.005 | ±0.002 (40 %) | (M2) 3.0 × 10⁻⁵ kg mmol⁻¹ min⁻¹ |
| km (mmol⁻¹ kg s⁻¹) | — | — | 0 (lower bound) | — | — |

Text [M]: p = 0.85 (Model 1) and 0.66 (Model 2); F_Asn 32 × 10⁻⁴ = 0.3 % of Asn consumed → acrylamide. kg ±0.003 /
0.006 is 50 %, printed 45 % (rounded estimate) [Z].
Comparison with Parker 2012 Table 2 [Z]: k6 18 ± 2 vs 16.5 ± 0.9 (×10⁻³ s⁻¹), k8 7.0 ± 1.7 vs 6.8 ± 0.8: overlap; k1
7.6 ± 1.0 vs 2.2 ± 0.03 (×10⁻⁴): do not overlap (different scheme).

### 4b. Table 1 (p. 239): 5 min at 165 °C, defatted dry weight [M]

| | Mal 0 % | 1 % | 2 % | 3 % | 4 % | 14 % | Glc 2 % (Parker 2012) | Fru 2 % (Parker 2012) |
|---|---|---|---|---|---|---|---|---|
| Sug_ini (mmol/kg) | – | 16.4 | 24.2 | 41.6 | 47.3 | 164 | 63.0 | 72.8 |
| Sug_fin (mmol/kg) | – | 7.1 | 13.5 | 20.7 | 26.3 | 111 | 27.3 | 50.4 |
| Sug_ini − Sug_fin | – | 9.2 | 10.7 | 20.9 | 21.0 | 52.4 | 35.8 | 22.4 |
| Acr_fin − Acr_ini (µmol/kg) | 4.5 | 10.9 | 20.1 | 17.5 | 19.4 | 20.8 | 72.7 | 75.4 |
| Colour (Agtron) | 70.8 | 63.2 | 56.1 | 52.9 | 53.4 | 34.6 | 26.3 | 34.8 |

"Sug" = the added sugar of that column. Footnote d calls Acr_fin − Acr_ini "concentration … after finish frying"
(it is a difference) [!]. Differences recompute within 0.1 except 1 % (9.3) and 14 % (53, from rounded 164/111) [Z].

### 4c. Fig. 1 (p. 238), t = 0 → t = 300 s [F]; Model 1 line at 300 s

| arm | acrylamide µmol/kg | Model 1 at 300 s | maltose mmol/kg | glucose mmol/kg |
|---|---|---|---|---|
| M0 (165 °C) | 1 → 5 | ~18 | — | 12 → 5 |
| M1 | 2 → 13 | ~18 | 16 → 7 | 10 → 7 |
| M2 | 3 → 23 | ~23 | 24 → 13 | 13 → 7 |
| M3 | 3 → 20 | ~28 | 41 → 21 | 16 → 10 |
| M4 | 4 → 23 | ~30 | 47 → 26 | 13 → 12 (scatter to 20) |
| M14 | 4 → 24 | ~43 | 164 → 111 | 11 → 15 (rises) |
| M2′175 (lot 2, 175 °C) | 4 → 49 | ~38 | 25 → 9 | 3 → 3 |
| M2′185 (lot 2, 185 °C) | 4 → 72 | ~72 | 25 → 8 | 3 → 6 |

Printed maxima [M]: M0 5.4, M1 13.4, M2-M14 20-24 µmol/kg. Glucose stays flat (M4) or rises (M14), which the authors
attribute to glucose released from maltose [Q]. Fructose, FAA, Asn for these arms: Fig. S1 only.

### 4d. Fig. 3 (p. 241), observed vs predicted, Model 1, all data sets; y = bx through origin [M]

glucose b = 1.0045 (R² 0.9583); fructose 1.0351 (0.9806); FAA 1.0225 (0.6558); acrylamide 0.9646 (0.9064); maltose
0.9964 (0.9914); asparagine 1.0067 (0.4536).

## 5. Flags and internal consistency checks

1. [!] **No Ea and no temperature or moisture treatment in this paper**, yet 175/185 °C arms are in the fit. The
   constants are 165 °C values inside an unshown 3-D transfer model (Balagiannis et al. 2016, ACS Symp. Ser.
   "Browned flavors", p. 67 [C], not on disk).
2. [!] **Model 1 misses the measured acrylamide in several arms** [F]: M0 ×3.5 over (18 vs 5), M14 ×1.8 over (43 vs 24),
   M3/M4 ×1.3-1.4 over, M2′175 ×0.8 under; the measured plateau over 2-14 % dips is not reproduced (model rises
   monotonically with maltose). The text calls the fit "very good" [V]. Same low-sugar over-prediction as Parker
   2012's G1/F1 arms.
3. [!] Text says R² "generally > 0.95"; Fig. 3 prints 0.906 (acrylamide), 0.656 (FAA), 0.454 (Asn). Slope range printed
   "0.96–0.103" — typo for 1.035. Regression through the origin inflates R².
4. [!] **Units slip**: maltose losses "10.7 and 52.4 μmol/kg" [V] in text; Table 1 has mmol/kg.
5. [!] Fig. 4 caption sets [FAA]₀ = 146 mmol/kg; text says 164 ("the same as in M0").
6. [!] Model 1 rate law printed as d[Acr]/dt = k3[Int2][Asn] ∗ F_Asn ∗ Asn: a stray "∗ Asn"; Model 2 version has none.
7. [!] Blanch 80 °C here vs 74 °C in Parker 2012 though "as described by Parker".
8. [!] The text introduces "Kinetic Mechanism 2, Fig. 1" [V]; it is drawn in Fig. 2 (Fig. 1 is data).
9. [!] Replicates and error bars absent; Asn and fructose data for the maltose arms are SI-only.
10. Checks that pass [Z]: Table 1 maltose ini/fin agree with Fig. 1 endpoints; Glc/Fru 2 % rows agree with Parker 2012
    Figs. 1/2 G5/F5 (63 → 27, 72 → 50; Δacr 72, 75); "8.2 to 86 µmol/kg" = Parker's F1 and F5 finals; F_Asn 0.3 %.

## 6. Existing records in the repo that cite this paper

`grep -rIn` (2026-11-03) of `data docs src scripts tasks` for `foodchem.2019.01.075`, `balagiannis2019`,
`Balagiannis.{0,30}2019`, maltose in `src/kinetic_core/`; also the Balagiannis records in
`data/lit/measurement_ledger/{acrylamide,trunk,sulfur}.yml` (none cite this paper).

| record | claim | verdict |
|---|---|---|
| `data/keys/pdf_dossier_aliases.yml` l. 107 (pending) | Balagiannis et al. 2019 Food Chem 284:236, DOI, kinetic model of acrylamide in finish-fried fries vs maltose (Reading) | **matches** |
| `src/kinetic_core/` acrylamide lane | — | **no use** |
| `src/kinetic_core/species.py` l. 357-431 | cites Mundt & Wedzicha **2004** (glucose + glycine, C/N of melanoidin) | different paper from this one's ref. Mundt & Wedzicha **2005** (maltose + glycine, JAFC 53:6798) [C]; no conflict |

No other record cites the paper.

## 7. Proposed ledger entries (`data/lit/measurement_ledger/acrylamide.yml`; NOT written)

```yaml
  - id: balagiannis2019_fig2_models1_2_constants_165c
    dossier: balagiannis2019_extraction.md
    paper: "Balagiannis, Mottram, Higley, Smith, Wedzicha & Parker 2019 Food Chem 284:236"
    lab: "Reading (Mottram / Parker / Elmore)"
    system: "par-fried Ranger Russet fries, maltose dip 0-14 % (this paper) + glucose / fructose dip 0-2 % (Parker 2012); Mal 16-164, Glc 3-63, Fru 6-72, FAA 130-240 mmol/kg defatted dry wt"
    medium: "finish-fried in canola oil, pilot fryer, ~220 g lots"
    temperature_c: [165, 175, 185]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [0, 5]
    species: [acrylamide, asparagine, glucose, fructose, maltose, total_free_amino_acids]
    kind: rate_constant
    location: "Fig. 2 table (Model 1: k1, k3, k6, k8, p, F_Asn, k1mal, k2mal; Model 2 adds kg, kf, km) with 95 % CI, all at 165 C (dossier sec. 4a)"
    status: digitised_in_dossier
    quant_class: sida
    simulable: false
    condition_quote: "The rate constants correspond to a temperature of 165 C."
    flags: [matrix_food, low_moisture_no_aw_term, ph_not_measured, no_replicates, "no_ea_printed", "local_constants_inside_slab_transfer_model", "complements: parker2012b_tab2_finish_fry_multiresponse_constants"]
    note: "Refit of Parker 2012 data + maltose arms with a changed scheme (Glc -> Int2 in one step, no Fru->Glc, Asn explicit): k1 3.5x below Parker's, F_Asn redefined; do not mix with parker2012b constants. Transfer model (Balagiannis 2016 3-D) and Ea not in paper. km = 0 (no direct maltose-Asn route)."
  - id: balagiannis2019_fig1_maltose_glucose_acrylamide_timecourses
    dossier: balagiannis2019_extraction.md
    paper: "Balagiannis, Mottram, Higley, Smith, Wedzicha & Parker 2019 Food Chem 284:236"
    lab: "Reading (Mottram / Parker / Elmore)"
    system: "par-fried fries, maltose dip 0/1/2/3/4/14 % (lot 1) and 2 % (lot 2)"
    medium: "finish-fried in canola oil at 165 C (lot 1) and 175 / 185 C (lot 2), quenched in liquid N2"
    temperature_c: [165, 175, 185]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [0, 5]
    species: [acrylamide, maltose, glucose]
    kind: time_series
    location: "Fig. 1 (8 arms x 11 times 0-300 s, with Model 1 lines); t = 0 and 300 s digitised (dossier sec. 4c)"
    status: digitised_in_dossier
    quant_class: sida
    simulable: false
    condition_quote: "Each batch (~220 g) was finish-fried at 165 C for 10 different time points"
    flags: [matrix_food, low_moisture_no_aw_term, heating_up_included_in_time, no_replicates, ph_not_measured]
    note: "Defatted dry-weight basis. Acrylamide plateaus at 20-24 umol/kg for 2-14 % maltose dips; Model 1 over-predicts M0 x3.5 and M14 x1.8. Fructose, FAA and Asn for these arms only in Fig. S1 (absent)."
  - id: balagiannis2019_tab1_maltose_glucose_fructose_dip_endpoints
    dossier: balagiannis2019_extraction.md
    paper: "Balagiannis, Mottram, Higley, Smith, Wedzicha & Parker 2019 Food Chem 284:236"
    lab: "Reading (Mottram / Parker / Elmore)"
    system: "par-fried fries, maltose dip 0-14 % and (from Parker 2012) 2 % glucose and 2 % fructose dips"
    medium: "finish-fried in canola oil, 165 C, 5 min"
    temperature_c: [165]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [5]
    species: [acrylamide, maltose, glucose, fructose, colour_agtron]
    kind: endpoint
    location: "Table 1 (added-sugar ini/fin/loss, acrylamide increase, Agtron colour) (dossier sec. 4b)"
    status: digitised_in_dossier
    quant_class: sida
    simulable: false
    condition_quote: "after finish frying at 165 C for 5 min"
    flags: [matrix_food, low_moisture_no_aw_term, no_replicates, "complements: balagiannis2019_fig1_maltose_glucose_acrylamide_timecourses"]
    note: "Directional: per mole of added sugar lost, acrylamide yield maltose 0.03-0.15 % vs glucose 0.20 %, fructose 0.34 % (dossier Z); colour per sugar lost also lower for maltose. Glc/Fru columns are Parker 2012 data (also parker2012b_figs1_2_finish_fry_timecourses_165c)."
  - id: balagiannis2019_figS1_all_responses_maltose_arms
    dossier: balagiannis2019_extraction.md
    paper: "Balagiannis, Mottram, Higley, Smith, Wedzicha & Parker 2019 Food Chem 284:236"
    lab: "Reading (Mottram / Parker / Elmore)"
    system: "par-fried fries, maltose dip 0-14 % (lot 1, 165 C) and 2 % (lot 2, 175 / 185 C)"
    medium: "finish-fried in canola oil"
    temperature_c: [165, 175, 185]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [0, 5]
    species: [fructose, sucrose, total_free_amino_acids, asparagine, acrylamide, maltose, glucose]
    kind: time_series
    location: "Appendix A, supplementary Fig. S1 (all responses); Figs. S2-S3 Model 2 fits"
    status: supplement_missing
    quant_class: sida
    simulable: false
    condition_quote: "The data are all displayed graphically in supplementary Fig. S1."
    flags: [matrix_food, low_moisture_no_aw_term, supplement_missing_s1]
    note: "Only source of Asn and fructose time courses for the maltose arms (Asn enters Model 1/2 explicitly). Download from the publisher's Appendix A if a fry hold-out is ever built."
```

Fig. 4 (relative sugar contributions) is a model simulation, not a data set: no entry.

## 8. Supplementary material

**Appendix A** (online at the article DOI, "Supplementary data"): **Fig. S1** (all measured responses: maltose,
glucose, fructose, sucrose, free amino acids, Asn, acrylamide vs time for every maltose arm) and **Figs. S2-S3**
(Model 2 fits). **Not on disk.** Model-relevant content only there: the Asn and fructose time courses of the maltose
arms. Needed only if the fries become a hold-out; not needed to judge the constants (confidence ~85 %).
