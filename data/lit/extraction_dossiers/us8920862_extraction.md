# US 8,920,862 B2, Nestec, 2014 — EXTRACTION (cysteine 6.46 g + xylose 24 g in 79.54 g 0.2 mol/L sodium phosphate pH 5.5, 5 mL in silanised crimped 20 mL vials, silicone bath 95 °C; MFT 9.4 and FFT 18.9 µmol/mol xylose after 6 h in buffer, 94.7 and 139 in a structured-lipid arm; SIDA HS-SPME-GC-MS; ONE time point printed)

### ONE ENDPOINT, NOT A TIME COURSE. The text says vials were withdrawn "at defined time intervals", but only the 6 h result is printed (no figure, no table for Example 5). The buffer arm is a clean SIDA endpoint (standards added after heating, two vials per time) at a very concentrated charge (xylose 1.45, cysteine 0.48 mol/kg solution; 3:1). Also printed: xylose loss with glycine in buffer, 120 °C, 30 min (22.3 %).

**Source on disk:** `data/articles/US8920862.pdf` (19 pp.: front page + page 2, three drawing sheets of emulsion cartoons,
specification cols 1-28). Text read with `pdftotext -layout -q`; Example 5 (cols 24-25, PDF pp. 17-18) re-checked on 300-dpi
renders (scratchpad `img/n862_pre.png`, `img/n862_left.png`). The three figures are schematic drawings of emulsion
structures — no data.

**Status of the evidence.** A patent, not peer reviewed; written to show the structured-lipid phase beats water.

## 0. Identity

| field | value |
|---|---|
| Title | "Maillard flavor compositions and methods for making such compostions" [sic] ✔ |
| Patent number | **US 8,920,862 B2** (5 claims, 3 drawing sheets) |
| Assignee | Nestec SA, Vevey (CH) |
| Inventors | Laurent Sagalowicz (Blonay, CH); Tomas Davidek (Correvon, CH); Florian Viton (Lausanne, CH); Haiqing Yu (Ellisville, MO, US); Martin Leser (Dublin, OH, US) |
| Appl. No. | 12/737,247; PCT/US2009/003711, PCT filed **Jun. 22, 2009**; §371 date Apr. 15, 2011 |
| PCT publication | WO2010/008452, Jan. 21, 2010 |
| Priority (as printed) | **US provisional 61/132,976, filed Jun. 24, 2008** |
| Prior publication | US 2011/0189367 A1, Aug. 4, 2011 |
| Date of patent (grant) | **Dec. 30, 2014** |
| Sibling | US 10,076,126 B2 (`us10076126_extraction.md`): same Nestlé group (Davidek, Viton), same SIDA-SPME method family |

## 1. Why it matters

Brief: a second cysteine + pentose MFT/FFT time course in water. **This patent does not deliver one** — it delivers a single
6 h endpoint at 95 °C for the buffer arm. Its value is (i) a clean level at a long time and a temperature below Schieberle's
100 °C, (ii) the only SIDA measurement in the family where the labelled standards were added after heating, and (iii) an
MFT/FFT ratio in buffer (0.50) that contradicts the sibling patent's (3.6 at 110 °C, 4 h).

## 2. Methods as they matter to a model

### 2a. Example 5, buffer arm ("Pre-reaction in phosphate buffer", col. 24 l. 59 – col. 25 l. 2) ✔ verified on the page image

- **Charge.** "A solution of cysteine (6.46 g) and xylose (24 g) in sodium phosphate buffer (79.54 g; 0.2 mol/L; pH 5.5)".
  (The brief's "79.5 g" is 79.54 g.) L-cysteine (Fluka), D-xylose (Fluka); buffer from NaH2PO4·H2O (Merck), Milli-Q water.
  Arithmetic (mine; MW 121.16, 150.13):
  - cysteine **53.32 mmol**, xylose **159.86 mmol** → **xylose : cysteine = 3.00 : 1**.
  - total solution mass 6.46 + 24 + 79.54 = **110.00 g**.
  - per kg solution: cysteine **0.485 mol/kg**, xylose **1.453 mol/kg**.
  - per kg buffer: cysteine **0.670 mol/kg**, xylose **2.010 mol/kg**.
  - per litre: density not printed; at an assumed ≈1.08 g/mL (≈22 wt % sugar) ≈ 0.52 mol/L Cys and 1.57 mol/L xylose
    (my estimate, ±5 %).
- **pH.** 0.2 mol/L phosphate set to pH 5.5 **before** dissolving 0.67 mol/kg cysteine — the precursor outnumbers the buffer
  ~3:1 (0.2 mol/L × ≈0.08 L buffer ≈ 16 mmol phosphate, ≈0.145 mol/kg solution, vs 0.485 mol/kg Cys), so the starting pH of the reaction solution may differ from 5.5.
  **pH of the solution not measured / not reported, before or after heating.** The results line says "at pH=5.5" (nominal).
- **Vessel.** "dispatched (à 5 ml) into silanised 20 mL crimp cap vials (Chromacol)" — 5 mL liquid, ≈15 mL headspace, sealed.
- **Heating.** "heated in a silicone bath at 95° C." Come-up not stated (5 mL in glass, plausibly 3-6 min; my estimate).
- **Sampling.** "At defined time intervals two vials were taken out of the silicone bath and internal standards were added
  (5.53 µg [2H2]-FFT in 50 µL pentane and 14.6 µg [2H3]-MFT in 50 µL pentane) into each vial. The vials were agitated with a
  vortex and after cooling down to room temperature they were stored over-night in a refrigerator." Analysed next day.
  → **Standards added at the end of heating (hot vial), not before.** Duplicate vials per time point.
  Standard per vial: [2H3]-MFT 14.6 µg / 117.2 g/mol = 0.125 µmol; [2H2]-FFT 5.53 µg / 116.2 = 0.048 µmol.
- **Analysis.** HS-SPME, PDMS/DVB 65 µm (Supelco); ≥ 1 h equilibration at 20 °C, 30 min exposure at 20 °C, no agitation;
  desorption 5 min at 250 °C (0.75 mm liner); GC 6890A / MSD 5973N, HP-PONA 50 m × 0.20 mm × 0.50 µm; oven 25 °C (2 min)
  → 40 °C/min → 50 °C (1 min) → 6 °C/min → 240 °C (10 min); EI 70 eV; **isotope dilution in scan mode, MFT m/z 114 vs
  [2H3]-MFT 117, FFT m/z 114 vs [2H2]-FFT 116.** Response factors not printed.

### 2b. Example 5, structured-lipid arm ("Pre-reaction in structured lipid", col. 25 ll. 3-20)

Cysteine 1.62 g + xylose 6.0 g in 17.38 g of the same buffer (25.00 g; Cys 13.37 mmol, xylose 39.97 mmol; ratio 2.99:1;
**0.535 mol/kg Cys and 1.599 mol/kg xylose per kg aqueous solution** — ~10 % more concentrated than the buffer arm).
0.75 g of it into 20 mL crimp vials containing Dimodan U (distilled monoglycerides, 2.82 g) and soybean oil (1.41 g)
(aqueous phase 15 wt % of 4.98 g). Heated at 95 °C until liquid, vortexed to form the mesophase, heated at 95 °C "for a
defined period of time", vortexed every 60 min. Standards (printed "5.53 g ... 50 L" — a µg/µL typesetting loss) added, then
**heated another 5 min** and vortexed; cooled, refrigerated overnight, same analysis.

### 2c. Example 1 (cols 20-22) — xylose consumption with glycine, 120 °C, 30 min

Third composition (the aqueous control): glycine 0.011 g (0.14 mmol) + xylose 0.021 g (0.14 mmol) in **10 mL phosphate buffer
0.2 mol/L pH 6.0** (**14 mmol/L each**), 30 mL Pyrex tube, 120 °C, 30 min (heating medium not stated for this one; the
other two compositions used a silicone bath). Residual xylose by HPAEC-PAD (CarboPac PA-1, duplicate injections,
external calibration).

## 3. Numbers, re-typed

### 3a. Example 5 results (col. 25 ll. 45-49), verbatim values

"After 6 hours of reaction at pH=5.5, the amount of volatile in micromol/mol xylose was 9.4 for MFT in buffer, 94.7 for MFT in
structured lipid, 18.9 of FFT in buffer and 139 of FFT in structured lipid phase." ✔ (page image)

| arm | compound | value as printed (µmol/mol xylose) | µmol/mol Cys (mine, ×3.0) | µmol/kg solution (mine) | ng/g solution (mine) |
|---|---|---|---|---|---|
| buffer | **MFT** | **9.4** | 28.2 | 13.7 | **1560** |
| buffer | **FFT** | **18.9** | 56.7 | 27.5 | **3140** |
| structured lipid | MFT | 94.7 | 283 | 151 per kg aqueous phase; 22.8 per kg vial | 2600 per g vial |
| structured lipid | FFT | 139 | 416 | 222 per kg aqueous phase; 33.5 per kg vial | 3820 per g vial |

Conversions: µmol/mol xylose × 1.453 mol xylose/kg solution (buffer arm) or × 1.599 mol/kg aqueous phase (lipid arm);
× 114.17 g/mol. Ratios: **MFT/FFT = 0.50 in buffer, 0.68 in structured lipid**; lipid/buffer = 10.1× (MFT), 7.4× (FFT).
No ± printed; "two vials" per time — whether the value is a mean of two is implied, not stated.

### 3b. Example 1 (col. 22 ll. 15-27) — residual xylose after 30 min at 120 °C, xylose + glycine 0.2 % / 0.1 % by weight

| composition | residual xylose | xylose degraded |
|---|---|---|
| water-in-oil microemulsion (Method 1) | 1.5 % | 98.5 % |
| structured oil-in-water emulsion (WO2007/060177) | 75.8 % | 24.2 % |
| **aqueous, 0.2 M phosphate pH 6.0, 14 mM each** | **77.7 %** | **22.3 %** |

Arithmetic: apparent first-order k(xylose) = −ln 0.777 / 0.5 h = **0.50 h⁻¹** at 120 °C, pH 6.0, 14 mM xylose + 14 mM glycine,
0.2 M phosphate (single determination, duplicate injection; includes come-up of 10 mL in a 30 mL tube).

### 3c. Everything else

Examples 2-4 and 6-8: palatability trials (dogs, 20; cats, 402) and sensory descriptions; formulations in Tables 1-5
(% w/w). Example 4's aqueous control (glucose 3.84 / xylose 6.42 / glycine 8.13 / cysteine 2.60 / proline 6.65 % in water,
105 °C, 60 min) gives only "light-orange coloration". **No H2S, no other volatile, no other time point, no pH reading
anywhere.**

## 4. Kinetic numbers the repository can use

Grades: **M** = printed number; **F** = figure-digitised; **Z** = not reported / zero information; **!** = flagged caveat.

| # | observable | value | unit | conditions | grade | source |
|---|---|---|---|---|---|---|
| 1 | charge, Cys / xylose | 0.485 / 1.453 | mol/kg solution | 6.46 g + 24 g in 79.54 g 0.2 M NaPi pH 5.5 | M (masses) | col. 24 |
| 2 | MFT, buffer, 6 h | 9.4 | µmol/mol xylose | 95 °C, sealed 20 mL vial, 5 mL, silicone bath; SIDA after heating | M ! (pH drift) | col. 25 |
| 3 | FFT, buffer, 6 h | 18.9 | µmol/mol xylose | same | M ! | col. 25 |
| 4 | MFT/FFT, buffer, 6 h | 0.50 | — | same | M (derived) | col. 25 |
| 5 | MFT / FFT, structured lipid, 6 h | 94.7 / 139 | µmol/mol xylose | W/O mesophase, 15 % aqueous | M (not an aqueous pot) | col. 25 |
| 6 | xylose loss with glycine | 22.3 | % in 30 min | 120 °C, 0.2 M NaPi pH 6.0, 14 mM each | M | col. 22 |
| 7 | MFT/FFT at other times | — | — | "defined time intervals" withdrawn, not printed | Z | — |
| 8 | final pH, H2S, cysteine / xylose consumption, error | — | — | — | Z | — |

## 5. What this can and cannot anchor

**Can:**
- **A 6 h endpoint level** for MFT and FFT in a concentrated cysteine + xylose phosphate pot at 95 °C, by SIDA with standards
  added after heating (the cleaner design of the two Nestlé patents).
- **An MFT/FFT ratio < 1** in buffer under excess xylose (3:1) at 95 °C.
- **A xylose-loss rate** (row 6) for an aqueous xylose + glycine pot at 120 °C, pH 6.0 — a core-lane check, not sulfur.

**Cannot:**
- **No time course** — the brief's hope for a second hour-scale series is not met by this document.
- Starting and final pH unknown; at 0.67 mol/kg cysteine the 0.2 M buffer cannot hold pH 5.5 against cysteine's own
  acid-base chemistry and the acids formed (my inference).
- Very high sugar (≈22 wt %): water activity and viscosity are not those of the dilute model pots.

## 6. Flags

1. **! Only the 6 h point is printed** though a time series was sampled.
2. **! pH 5.5 is the buffer's nominal pH, not the reaction's.** Cys : phosphate ≈ 3 : 1 molar.
3. **Typesetting loss** in the lipid arm: "5.53 g ... in 50 L" = 5.53 µg in 50 µL (the buffer arm prints µg/µL correctly).
4. **The buffer-arm MFT/FFT (0.50) contradicts the sibling patent's buffer arm (3.6 at 110 °C, 4 h, 1:1, 0.25 mol/kg)** by 7×.
   Differences: temperature, sugar excess, concentration, time, and standard-addition timing (sibling adds before heating).
   No single cause can be assigned from the patents.
5. Standards added to a **hot** vial and vortexed: losses of the native thiols during cool-down are corrected, but losses
   during the 6 h are not (correct for "amount present at 6 h").
6. Registry: patent not checked against `data/keys/papers.yml` (no edit made).

## Use in the model

- **Not a hold-out time course.** Use as a **single-point validation level** (grade M !): cysteine 0.485 + xylose 1.453 mol/kg,
  0.2 M phosphate nominal pH 5.5, 95 °C, 6 h, sealed → MFT 13.7 µmol/kg (1.56 µg/g), FFT 27.5 µmol/kg (3.14 µg/g), each with a
  declared band for the unknown pH (I would not tighten below ×2 without a pH).
- **Directional claim:** FFT > MFT in buffer at 95 °C, 6 h, xylose:Cys 3:1 (MFT/FFT = 0.50). Pair with the sibling's MFT > FFT at
  110 °C, 4 h, 1:1 as a joint test of what controls the ratio (sugar excess favours furfural → FFT? temperature favours the
  MFT branch?) — report, do not fit, until a pot with a measured pH reproduces either.
- **Core lane:** row 6 (xylose + glycine, 120 °C, pH 6.0, 22 % loss in 30 min) can be scored as a one-point xylose-loss check.
