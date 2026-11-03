# Gökmen, Palazoğlu & Şenyuva 2006 — EXTRACTION (file stem `gokmen2006`)
DOI: 10.1016/j.jfoodeng.2005.08.030 (printed p. 972 footer and in the PDF title field). Crossref 2026-11-03:
"Relation between the acrylamide formation and time–temperature history of surface and core regions of French
fries", J. Food Eng. **77**(4):972-976, authors Gökmen, Palazoğlu, Şenyuva, issued 2006-12 (online 26 Sep 2005).
**Matches the brief and the PDF.**

### French-fry strips (8.5 × 8.5 × 70 mm, variety not named) fried in sunflower oil at 150 / 170 / 190 °C for 3, 6, 9 min; surface and core T logged every 1 s by thermocouple; acrylamide by LC-APCI-MS (external standard) in the whole strip, a "2 mm surface layer" and the remaining core; whole-strip moisture. **No sugar, no asparagine, no a_w, no moisture method, no error bars, no kinetics.** This is the primary of the data Chan 2020 re-plotted.

**Read 2026-11-03 from `data/articles/gokmen2006.pdf`** (5 pp., Elsevier, born-digital). Prose from `pdftotext
-layout`. Figures read from renders in scratchpad `ingest/gokmen2006_work/`: Fig. 2 (p. 3, 200 dpi `fig2b-3.png`),
Fig. 3 (p. 4, 300 and 600 dpi `fig3hi-4.png`, `fig3zoom-4.png`, `fig3right-4.png`), Fig. 4 and Fig. 5 (p. 4, 300 dpi
`fig4hi-4.png`, `fig5hi-4.png`; Fig. 5 0-230 s at 600 dpi `fig5early-4.png`), Fig. 6 (p. 5, 300 dpi `fig6hi-5.png`).
Read errors stated per table.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Department of Food Engineering, **Hacettepe University**, Ankara; corresponding author **Vural
Gökmen**; Palazoğlu (Mersin, heat-transfer side), Şenyuva (TÜBİTAK ATAL, which ran the LC-MS). **Existing id
`hacettepe_gokmen`** in `data/keys/laboratories.yml` (string "Hacettepe (Gokmen)", added 2026-11-03; its evidence
line already names this paper). Same lab: `gokmen2012_extraction.md`, `palazoglu2008_extraction.md` (same strips,
same fryer, same two authors), `hamzalioglu2018`, `kocadagli2016*`, `goncuoglu2016` dossiers. Re-plotted by
`chan2020_extraction.md` (different lab).

## 1. What this paper answers for the model

| model need | what Gökmen 2006 gives | tag |
|---|---|---|
| acrylamide vs time in a fried potato, whole strip | 3 T × 3 times (Fig. 3 = Fig. 2), 9 min values printed (10 / 418 / 1422 ng/g) | [F] / [M] |
| acrylamide surface vs core | Figs. 4-6, 3 T × 3 times; 9 min surface 72 / 2747 / 6476, core 0 / 0 / 376 ng/g printed | [F] / [M] |
| measured product temperature history | surface and core T at 1 s, oil T, 0-540 s, 3 oil T (Figs. 4-6) | [F] |
| moisture loss | whole-strip moisture (%), 4 times × 3 T (Fig. 3); method not described | [F] [!] |
| precursor charge (Asn, glucose, fructose) | **none**; variety not named | [!] |
| a_w, pH, oil uptake | **none** | [!] |
| below-120 °C formation | 115 °C oil, 40 min: 320 ng/g (basis whole / surface unstated) | [M] |
| rate constants / Ea | **none** | — |

**Answers.**
1. **Is Chan 2020's re-plot faithful to the primary? Yes** (confidence ~90 %). Chan's 170 °C points (sec. 4b of
   `chan2020_extraction.md`) against my reading: water 76.5 / 42.0 / 32.4 / 24.9 vs 77 / 42.4 / 32.7 / 25.2 %; whole-
   strip acrylamide 103 / 250 / 417 vs 110 / 245 / 418 (printed) ng/g; surface T 77 / 104 / 117 / 124 / 125 at 10 / 30
   / 50 / 70 / 90 s vs 79 / 105 / 118 / 124 / 125; core T 25 / 65 / 91 / 101 / 104 vs 25 / 64 / 91 / 101 / 103; surface
   T at 3 / 6 / 9 min 128 / 136 / 141.5 vs ~130 / ~137.5 / ~142.5 (Chan ~1 °C low late). All within the joint read
   error. **Chan used the whole-strip acrylamide**, not the surface layer, and only the 170 °C arm of three.
2. **What the primary adds over Chan: the 150 and 190 °C arms, the surface / core split and the oil-T trace**
   (confidence ~95 %). Three oil temperatures with measured surface T histories make this a small T-series frying
   hold-out for any measured-T(t) acrylamide path, *if* a precursor charge can be assumed (none printed; see 3).
3. **Can the pot be simulated? Not as published** (confidence ~85 %): no Asn, no sugar, no variety, no a_w. The
   same group's Palazoğlu & Gökmen 2008 (Agria, 1.4 ± 0.2 g/kg reducing sugar, 7.2 ± 0.6 g/kg Asn) is a different
   tuber lot two years later; using it here is an assumption (Chan did exactly that, via its ref. 15). P(that the
   2006 tubers had 1.4 g/kg within ±50 %) ~45 % [Z]: Palazoğlu's 4-min whole strip at 170 °C (313 ng/g) is ~1.8 ×
   the 2006 value interpolated at 4 min (~170 ng/g), different analytics notwithstanding.
4. **Is the surface/core split usable quantitatively? Only the whole-strip numbers** (confidence ~75 %). The
   surface + core values do not reconstruct the whole-strip values for any plausible mass split of a 2 mm layer
   (flag 1).

## 2. Conditions
- Potato strips 8.5 × 8.5 mm (French-fry cutter), length 70 mm; variety, supplier, storage **not stated**.
- Sunflower oil, 5 L bath (Precisterm, J.P. Selecta); oil 150, 170, 190 °C [M]; oil T trace in Figs. 4-6 holds
  ~150 / ~170 / ~190 °C within ±2 °C [F] (no drop on loading visible).
- **Six strips per fry**, "sampled in duplicates at 3, 6, and 9 min" [V]; plus one 115 °C fry, 40 min (text).
- Thermocouples: 36-gauge type T, one inserted until flush with the surface, one at the core; logged every 1 s
  (Keithley 2700 + 7700 multiplexer) [M]. Thermocouple placement in Fig. 1 schematic only.
- Moisture: plotted (Fig. 3, "Moisture, %"), **no method, no basis** given (wet basis assumed, ~80 %) [!].
- No pH, a_w, oil content, colour.

## 3. Analytical methods
LC-APCI-MS, Agilent 1100, Inertsil ODS-3, isocratic; SIM m/z 72 and 55. **External calibration** (working standards
0.1-2.0 µg/mL); **no isotope-labelled internal standard** stated [M]. Sample: 1 g finely ground fried sample + Carrez,
made to 10 mL; so ng/g is per g of fried product as analysed (wet, oil-containing basis) [Z, ~85 %]. Method validation
cited to Şenyuva & Gökmen 2005 [C]; **no LOD/LOQ, recovery, or ± printed**. Surface sample = "2 mm surface layer"
cut off with a knife; core = the rest; whole strips analysed separately [M].

## 4. Data
### 4a. Whole-strip acrylamide (ng/g) and moisture (%), Fig. 3 (= Fig. 2 for acrylamide) [F], printed [M] in bold
Read error ±10 ng/g (±30 at 190 °C upper range), ±0.5 % moisture. Times are 180 / 360 / 540 s.

| oil T | AA 3 min | AA 6 min | AA 9 min | moisture 0 | 3 min | 6 min | 9 min |
|---|---|---|---|---|---|---|---|
| 150 °C | ~0 | ~3 | **10** | 77 | 48.5 | 38.1 | 33.0 |
| 170 °C | 110 | 245 | **418** | 77 | 42.4 | 32.7 | 25.2 |
| 190 °C | 441 | 896 | **1422** (fig. ~1420) | 77 | 35.2 | 24.0 | 15.6 |

Fig. 2 (AA vs oil T at 3 / 6 / 9 min) re-plots the same points: 170 °C 110 / 251 / 423, 190 °C 442 / 902 / 1425 [F].
Text: moisture 38-45 % reached at ~6, 3 and < 3 min at 150, 170, 190 °C [M]; 190 °C × 3 min ≈ 40 × 150 °C × 9 min
[M] (441 / 10 = 44 [Z]).

### 4b. Surface / core acrylamide (ng/g), Figs. 4-6 [F], printed [M] in bold
Read error ±2 ng/g (Fig. 4), ±30 (Fig. 5), ±50 (Fig. 6 surface), ±20 (Fig. 6 core).

| oil T | surface 3 / 6 / 9 min | core 3 / 6 / 9 min |
|---|---|---|
| 150 °C | 0 / 51.8 / **72** | 0 / 0 / 0 |
| 170 °C | 907 / 1894 / **2747** | 0 / 0 / 0 |
| 190 °C | 2635 / 4781 / **6476** | 156 / 281 / **376** |

### 4c. Temperature histories (°C), Figs. 4-6 [F]; read error ±1.5 °C for 0-230 s (600 dpi crops of all three figures, `fig4early`, `fig5early`, `fig6early`), ±2.5 °C later
| t (s) | 150: surface / core | 170: surface / core | 190: surface / core |
|---|---|---|---|
| 0 | 21 / 21 | 21 / 20 | 22 / 22 |
| 10 | ~60 / ~25 | 79 / 25 | ~98 / ~27 |
| 30 | ~91 / ~55 | 105 / 64 | ~122 / ~70 |
| 50 | ~104 / ~86 | 118 / 91 | ~132 / ~95 |
| 70 | ~110 / ~100 | 124 / 101 | ~139 / ~103 |
| 100 | ~113 / ~104 | 125 / 103.5 | ~145 / ~104 |
| 180 | ~115 / ~104 | ~130 / 103 | ~148 / ~104 |
| 360 | ~117 / ~104 | ~137.5 / 102.5 | ~156 / ~103 |
| 540 | ~117 (max ~119.5 at ~520 s) / ~104 | ~142.5 / 102.5 | ~164 / ~104 |

Oil: 150.5 / 170 / 190.5 °C, flat [F]. Text: core never above 103-104 °C; surface ≤ 120 °C at 150 °C [M]. The
"surface" thermocouple is flush with the surface, i.e. it reads a near-surface point inside the crust, not the oil-crust
interface [Z].

## 5. Flags and internal consistency
1. [!] **Surface + core do not add up to the whole strip.** If all AA at 150 / 170 °C is in the surface layer, that
   layer must be (whole) / (surface) = 10/72 = 0.14, 418/2747 = 0.15 of the analysed mass; at 190 °C (1422 − 376) /
   (6476 − 376) = 0.17 [Z]. A 2 mm layer off four faces is 72 % of the raw cross-section (47 % for two faces); after
   drying (surface ~3 %, core ~65 % water) its mass share is still ≥ ~0.4 [Z]. The constant ratio 0.14-0.17 says the
   effective surface sample was thin (~0.3-0.5 mm equivalent) or the bases differ. P(the "2 mm" layer is not what
   produced the surface numbers, or the whole-strip and surface values are on different bases) ~75 %. Use the
   whole-strip values for any model check; treat the surface values as a profile index.
2. [!] No precursor data, no variety, no moisture method; ± not shown though duplicates were run.
3. [!] External calibration, no internal standard, no LOD; the 150 °C whole-strip values (3-10 ng/g) are almost
   certainly below a typical LC-MS LOQ of that era (~20-50 ng/g; Palazoğlu 2008 gives LOQ 50 ng/g for its GC-MS) [Z].
4. [!] 115 °C / 40 min / 320 ng/g: basis (whole vs surface) and moisture not given; the "below 120 °C" claim rests on
   a near-surface thermocouple that may read below the crust surface T (flag 4c note).
5. Text and figures agree on every printed value (418, 1422, 2747, 6476, 376, 72, 10) within read error.

## 6. Existing records in the repo that cite this paper
`grep -rIn` (2026-11-03) of `data/ docs/ src/ scripts/ tasks/` for `jfoodeng.2005.08.030`, `gokmen2006`, "Gokmen 2006",
"Gökmen 2006", "Palazoglu": hits only in `chan2020_extraction.md`, the two `chan2020_*` ledger entries
(`data/lit/measurement_ledger/acrylamide.yml` l. 1226-1261), and `data/keys/laboratories.yml` l. 62 (evidence line).

**`chan2020_figs2_4_gokmen2006_fries_170c_replot`**, field by field:
- `paper` / `lab`: match (authors, journal, volume, page, DOI; Hacettepe).
- `system` "77.2 % water": t0 moisture ~77 % [F] — **match**. "reducing sugar 6.8 mol/m3 (as computed by Chan from
  ref. 15)": **not from this paper**; Chan's ref. 15 is Palazoğlu & Gökmen 2008 (other tubers). Mismatch of
  attribution: the sugar belongs to the 2008 pot.
- `medium` "single strip per fry": **mismatch** — six strips per fry in a 5 L bath, duplicates per time. (Single strip
  per fry is Palazoğlu 2008.)
- `temperature_c: [170]`: correct for Chan's subset; **the primary has 150 / 170 / 190 (+115) °C**.
- `species`, `time_min [0, 9]`, `kind: time_series`: match.
- `quant_class: unknown`: primary is LC-MS with external standards → `calibrated` (no IS).
- `note` "basis of ug/kg not stated": still true in the primary (ng/g of ground fried sample, implicitly as-is).
  "Asn not given": true in the primary too.
- Chan dossier sec. 4b numbers: all **match** (answer 1).
**Recommendation**: supersede this entry by the two below (keep it as a pointer, flag `secondary_source`,
`superseded_by_primary`), and correct `medium` and the sugar attribution.

## 7. Proposed ledger entries (acrylamide.yml; NOT written)
```yaml
  # ---- gokmen2006 (ingested 2026-11-03) ----
  - id: gokmen2006_fig3_whole_strip_acrylamide_moisture_150_170_190c
    dossier: gokmen2006_extraction.md
    paper: "Gokmen, Palazoglu & Senyuva 2006 J Food Eng 77:972 (10.1016/j.jfoodeng.2005.08.030)"
    lab: "Hacettepe (Gokmen)"
    system: "potato strips 8.5 x 8.5 x 70 mm, variety not named, ~77 % moisture at t0; Asn and sugars not measured"
    medium: "deep-fat frying in sunflower oil, 5 L bath, six strips per fry, duplicates per time"
    temperature_c: [150, 170, 190]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [0, 9]
    species: [acrylamide, water content]
    kind: time_series
    location: "Fig. 3 (= Fig. 2) whole-strip acrylamide and moisture at 3/6/9 min x 3 oil T; 9 min values printed p. 975 (dossier sec. 4a); 115 C x 40 min = 320 ng/g in text"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "Six potato strips were placed in the fryer and sampled in duplicates"
    flags: [matrix_food, low_moisture_no_aw_term, heating_up_included_in_time, concentrations_not_stated]
    note: "Primary of chan2020_figs2_4 (Chan used the 170 C arm). No precursor charge: simulable only with an assumed Asn/sugar (Palazoglu 2008 Agria lot is a different tuber batch). External-standard LC-MS, no IS, no LOD; 150 C values near zero. Moisture method not stated."
  - id: gokmen2006_figs4_6_surface_core_temperature_acrylamide
    dossier: gokmen2006_extraction.md
    paper: "Gokmen, Palazoglu & Senyuva 2006 J Food Eng 77:972 (10.1016/j.jfoodeng.2005.08.030)"
    lab: "Hacettepe (Gokmen)"
    system: "potato strips 8.5 x 8.5 x 70 mm; 2 mm surface layer vs remaining core"
    medium: "deep-fat frying in sunflower oil; type-T thermocouples flush with surface and at core, 1 s logging"
    temperature_c: [150, 170, 190]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [0, 9]
    species: [acrylamide, surface temperature, core temperature]
    kind: time_series
    location: "Figs. 4-6 (surface and core T every ~10 s marker, oil T; surface and core acrylamide at 3/6/9 min); 9 min values printed (dossier sec. 4b-4c)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "A 2 mm surface layer (the surface sample) was removed from each strip"
    flags: [matrix_food, low_moisture_no_aw_term, heating_up_included_in_time, concentrations_not_stated, mass_balance_inconsistent]
    note: "Measured surface/core T(t) at three oil T: the thermal input for a measured-T(t) path. Surface + core acrylamide do not reconstruct the whole strip (surface share needed 0.14-0.17 vs >= ~0.4 for a 2 mm layer); use surface values as a profile index only."
```

## 8. Supplementary material
None referenced; no SI, appendix or data statement. Everything is in Figs. 2-6 and the text.
