# Kondjoyan, Chevolleau, Portanguen, Molina, Ikonic, Clerjon & Debrauwer 2016 — EXTRACTION (file stem `kondjoyan2016`)
Food Chem. **213**:641-646, DOI 10.1016/j.foodchem.2016.06.118 (printed p. 1). Received 17 Jul 2015, online 1 Jul 2016.
Crossref 2026-11-03: "Relation between crust development and heterocyclic aromatic amine formation when air-roasting a
meat cylinder", Food Chemistry 213, 641-646, authors Kondjoyan, Chevolleau, Portanguen, Molina, Ikonic, Clerjon,
Debrauwer, issued 2016-12. **Matches the PDF and the brief.** No naming caveat.

### Beef *Longissimus thoracis* cylinders (48 mm Ø × 50-55 mm) roasted under a hot-air jet at four jet temperatures (124 / 158 / 192 / 210 °C) for 20-90 min; 7 HAAs (ng/g freeze-dried crust) and coloured-crust thickness (mm) per condition (Table 3), creatine vs depth at 192 °C / 60 min (Fig. 2). **No temperature-time curve of surface, crust or core is printed; no water content; no kinetic constant.**

**Read 2026-11-03 from `data/articles/kondjoyan2016.pdf`** (6 pp., Elsevier, born-digital). Prose from
`pdftotext -layout`; **Table 1 (p. 642), Table 3 and Fig. 2 (p. 644) read by eye from 150 dpi page renders**
(scratchpad `ingest/kondjoyan2016_work/p-2.png`, `p-4.png`); Fig. 3 (p. 645) from a 150 dpi crop (`fig3-5.png`).
Table 2 (literature compilation) from the text layer, checked against the render of p. 643 only for the "present
study" rows. Fig. 2 bars read ±0.03 g/100 g DM.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** UR370 Qualité des Produits Animaux (QuaPA), **INRA Clermont-Ferrand-Theix**, Saint-Genès-Champanelle
(heat transfer, cooking device, crust); HAA analytics at INRA UMR1331 TOXALIM, Toulouse (Chevolleau, Molina,
Debrauwer). Corresponding author **Alain Kondjoyan**; senior author for the analytics Laurent Debrauwer.
**Not in `data/keys/laboratories.yml`: NEW**, proposed id **`inra_qpa_kondjoyan`**, string
"INRA Clermont QuaPA (Kondjoyan)". Same INRA unit (UR370 QuaPA) but a different team from `thomas2014_extraction.md`
(Berdagué, aroma/GC-O); do not merge for leave-one-out on that ground alone. Same group's cited precursors, not on
disk: Kondjoyan 2010a (Food Chem 119:19, superheated steam slices), 2010b (Food Chem 123:659, modelling HAAs in slices),
2014 (Meat Sci 97:323), Portanguen 2014 (Food Bioprocess Technol 7:3308, crust mechanisms, the source of Fig. 3).

## 1. What this paper answers for the model

| model need (mechanism item 8, thermal history) | what Kondjoyan 2016 gives | tag |
|---|---|---|
| surface / crust T(t) | **not printed.** IR pyrometer measured surface T "continuously" [Q], but only Table 1's jet-T ranges appear | [!] |
| core T(t) | **none** | — |
| air T, velocity, humidity, h | pipe-outlet T 160/190/225/260 °C, jet T 2 mm above surface 120-220 °C (ranges), nozzle-surface distance 36 mm [M]; **no jet velocity, no humidity, no h** | [M] / [!] |
| depth temperature structure | Fig. 3: thickness of the ≥200 / 150 / 120 / 110 / 100 °C zones after 210 °C / 70 min (single condition, cited from Portanguen 2014) | [C] |
| water loss / crust moisture vs time | **none** (no weight loss, no water content, no aw) | — |
| crust thickness vs time | Table 3: coloured-crust thickness (mm) at 11 conditions | [M] |
| HAA amounts vs time | Table 3: 7 HAAs, ng/g FDP, 4 jet T × 2-3 times, ± of triplicate analyses | [M] |
| precursor gradient | Fig. 2: creatine (g/100 g DM) at 5 depths, 192 °C / 60 min; Table 2: creatine mg/100 g raw / 192 °C 20 min / 210 °C 90 min | [F] / [M] |
| HAA formation k, Ea | **none fitted** (authors' ANOVA only) | — |

**Answers.**
1. **Measured T-time profiles? No** (confidence ~97 %). The IR surface record exists in the lab but is not shown;
   the paper prints jet-temperature ranges (Table 1) and uses a nominal "Roasting Temp" (124/158/192/210 °C, Table 3)
   that falls inside the jet ranges, i.e. it is the **air-jet temperature, not the product temperature** (~85 %;
   the methods sentence says Table 1 gives "actual temperatures reached on the surface" [V] while the column is headed
   Tjet [!]). Fig. 3 (cited) is the only in-product temperature information: ≥200 °C in the top 0.5 ± 0.5 mm,
   150 °C at 1.8, 120 °C at 3.5, 110 °C at 4.0, 100 °C (evaporation front) at 5.0 ± 0.5 mm after 70 min at 210 °C.
2. **Water loss / crust thickness vs time? Crust thickness only** (confidence ~95 %). Thickness rises with time at
   every T (Table 3); the "linear with time" claim is Portanguen 2014's [C]. At 210 °C it accelerates
   (0.033 mm/min over 20-40 min, 0.065 mm/min over 40-90 min [Z]), so linearity is not shown here [!].
3. **HAA vs time and units:** ng compound per g freeze-dried product (FDP) of the crust plus a < 0.5 mm non-brown
   margin; MeIQx 1.3 → 417 ng/g FDP across the grid; PhIP 0.3 → 256 (sec. 4a). **No per-fresh-weight or per-area
   values, no crust mass**, so exposure or a mass balance cannot be built [!].
4. **Fitted constants? None.** An indicative apparent Ea from MeIQx at 20 min vs jet T (124 → 210 °C) is
   ~92 kJ/mol [Z: ln(181.1/1.3) / (1/397.15 − 1/483.15 K) × R], but it uses air T, not reacting-layer T, and crust
   thickness changes 11-fold over the same span; it should not be carried (P(meaningful) ~15 %).
5. **Usable as a ready heating programme for a Maillard pot? No** (confidence ~90 %). Nothing gives T(t) at the
   surface or crust. What it does give: (a) a qualitative constraint that the reacting layer of a roast is much
   colder on average than the air (the authors' explanation of half the HAAs of iso-T slices [Q]); (b) a
   crust-thickness-vs-time series that a future transfer model could be checked against; (c) one cited depth
   isotherm map (Fig. 3). For an actual surface T(t) under the same jet, the source is Portanguen 2014 or Kondjoyan
   2010b (not on disk; P(they print surface T curves) ~70 %).

## 2. Conditions

- **Meat** [M]: *Longissimus thoracis*, 18-month heifers, aged 12 d vacuum-packed, frozen at −20 °C; cylinder 48 mm Ø ×
  50-55 mm, in a PEEK holder, spring-pushed against a top grid (partial shrinkage compensation). Composition, DM, fat,
  pH: **not printed** [!]. DM ≈ 27-29 % of cooked centre meat [Z, from 501 mg/100 g ≙ 1.85 g/100 g FDP and the Fig. 2
  caption's 550 mg/100 g ≙ 1.9 g/100 g DM].
- **Heating** [M]: open vertical air jet onto the top face; pipe outlet 160 / 190 / 225 / 260 °C → jet T (thermocouple
  2 mm above the surface) 120-130 / 150-160 / 190-200 / 200-220 °C; nozzle 36 mm from surface. Times (Table 1): 20, 40,
  90; 20, 40; 20, 30, 60, 90; 20, 40, 60, 90 min. **The 60 min points of the two hottest arms are not in Table 3** [!].
  Mimics a forced-convection oven [Q]. Humidity of the jet: not stated (dry air implied).
- **Stop** [M]: slide under a 45-55 m/s jet of 3-5 °C air (Ranque-Hilsch vortex tube).
- **Sampling** [M]: one half of the top cut for crust thickness (binocular images, ≥ 5 locations, GIMP + Matlab);
  HAAs on the coloured crust plus "less than 0.5 mm" of non-brown area (Fig. 1) — the thickness slice is described as
  crust + 2-3 mm [!, two different cuts, not explained]. Creatine: slices cut from the top at 192 °C / 20 and 60 min
  ("four 5 mm-thick slices" in the text, five depths in Fig. 2 [!]).
- **Replicates** [M]: each cooked sample extracted and analysed in triplicate; the number of cylinders per condition
  is **not stated** (± is probably analytical SD on one cylinder, ~70 %) [!].

## 3. Analytical methods

- **HAAs** [M]: 13 HAAs by UHPLC-APCI-MS/MS (Accela 600 + TSQ Vantage, MRM, 2 transitions each), Hypersil GOLD C8;
  1 g lyophilised sample, 1 M NaOH, IS **TriMeIQx** (50 pg/µL; a structural analogue, not isotope-labelled), LLE
  (Extrelut, CH2Cl2) + Oasis MCX SPE. Linear r² > 0.99 over 1-1000 pg/µL; LoD 0.015-0.05, LoQ 0.05-0.2 ng/g; area-ratio
  RSD 5-10 %; no matrix effect by post-column infusion [Q]. Recovery: not printed [!]. Trp-P-1, Trp-P-2, AαC, MeAαC
  never above traces; Glu-P-1 and IQ near LoQ only at 210 °C / 90 min.
- **Creatine** [M]: ¹H NMR (3.03 and 3.92 ppm singlets) in pH 5.0 phosphate extract vs a TSP-d4 capillary calibrated to
  a 1.76 mmol/L creatine standard. Creatinine not measured.
- **Crust thickness** [M]: visual edge of the brown-coloured zone on calibrated binocular images; mean ± SD over
  ≥ 5 locations (SD not printed in Table 3).
- **Statistics**: ANOVA + Tukey HSD (α 0.05); grey cells in Table 3 = no significant time effect.

## 4. Data

### 4a. Table 3 (p. 644): HAAs in the developing crust, ng/g FDP, ± as printed (n = 3 analyses) [M]

| line | t (min) | jet T (°C) | crust (mm) | MeIQx | IQx | 4,8-DiMeIQx | PhIP | MeIQ | DMIP | Glu-P-2 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 20 | 124 | 0.17 | 1.3 ± 0.2 | nd | nd | 0.3 ± 0.1 | nd | nd | nd |
| 2 | 40 | 124 | 0.19 | 19.4 ± 0.6 | 1.7 ± 0.2 | nd | 5.1 ± 0.2ᵍ | traces | nd | nd |
| 3 | 90 | 124 | 0.42 | 41.1 ± 7.0 | 2.3 ± 0.1 | 2.2 ± 0.1 | 5.8 ± 0.7ᵍ | traces | 8.7 ± 0.9 | 0.5 |
| 4 | 20 | 158 | 0.38 | 10.9 ± 5.7ᵍ | 0.51 | traces | 3.6 ± 1.4ᵍ | traces | nd | nd |
| 5 | 40 | 158 | 0.72 | 13.7 ± 1.3ᵍ | nd | traces | 1.2 ± 0.5ᵍ | traces | 2.9 | trace |
| 6 | 20 | 192 | 0.86 | 35.5 ± 0.5 | 2.2 ± 0.3 | 2.0 ± 0.1 | 14.5 ± 1.5 | nd | 6.3 ± 0.7 | traces |
| 7 | 30 | 192 | 1.42 | 88.5 ± 8.8 | 5.1 ± 0.6 | 5.2 ± 0.2 | 21.0 ± 1.6 | nd | 13.1 ± 1.9 | traces |
| 8 | 90 | 192 | 3.85 | 217.0 ± 1.6 | 9.8 ± 0.3 | 11.9 ± 0.2 | 98.0 ± 2.2 | 1.1 ± 0.2 | 41.0 ± 1.5 | 21.6 ± 0.7 |
| 9 | 20 | 210 | 1.94 | 181.1 ± 2.1 | 10.3 ± 0.3ᵍ | 11.4 ± 0.3 | 95.0 ± 1.7 | 0.5 ± 0.2 | 27.7 ± 0.4 | 9.4 ± 0.5 |
| 10 | 40 | 210 | 2.59 | 220.4 ± 8.6 | 11.7 ± 1.0ᵍ | 14.0 ± 0.5 | 156.9 ± 4.8 | 1.8 ± 0.2 | 48.8 ± 2.2 | 24.9 ± 1.4 |
| 11* | 90 | 210 | 5.82 | 417.0 ± 14.8 | 21.6 ± 1.0 | 26.9 ± 0.9 | 255.9 ± 10.5 | 2.7 ± 0.2 | 99.7 ± 6.1 | 48.0 ± 5.9 |

ᵍ grey in the PDF (no significant time effect). * "Crust too thick, not accepted by consumers" [V]. nd = not
detected. Line 5 is the authors' recommended low-HAA roast condition [Q].

### 4b. Fig. 2 (p. 644): creatine vs depth, 192 °C jet, 60 min, g/100 g DM, ±0.03 [F]

| depth (mm) | 2.5 | 5 | 10 | 15 | 20 | raw (line) |
|---|---|---|---|---|---|---|
| creatine | 0.70 | 1.55 | 1.65 | 1.73 | 1.85 | 1.9 ± 0.2 [M] |

Text [M]: centre 501 mg/100 g cooked meat = 1.85 g/100 g FDP; 20 min profile "ca. 20 % less" (not shown) [Q].
Table 2 present-study rows [M]: raw 470.7, 192 °C / 20 min 416.9, 210 °C / 90 min 101.0 mg/100 g (basis not stated;
the 210 °C / 90 min creatine run is not described in the methods [!]).

### 4c. Fig. 3 (p. 645): zone thickness after 210 °C / 70 min [C, Portanguen 2014; printed values]

T ≥ 200 °C: 0.5 ± 0.5 mm; 150 °C: 1.8 ± 0.5; 120 °C: 3.5 ± 0.5; 110 °C: 4.0 ± 0.5; 100 °C (evaporation front):
5.0 ± 0.5 mm; below, the liquid-dominant zone ≤ 100 °C. (Read as the depth of each isotherm, ~80 %.)

## 5. Flags and internal consistency

1. [!] "Roasting Temp" = jet (air) T, not product T; the methods call Table 1 surface temperatures while its column
   is Tjet. The IR surface record is never shown.
2. [!] Text p. 645: MeIQx and PhIP at 158 °C / 20 min "approximately 35 times lower" [V] than at 158 °C / 40 min.
   Table 3 gives 10.9 vs 13.7 (MeIQx, ×1.3) and 3.6 vs 1.2 (PhIP, **decrease**). The ×35 is not reproducible from any
   pair of rows (124 °C 20 → 40 min gives ×15 / ×17).
3. [!] 158 °C is non-monotonic vs 124 °C: at 40 min MeIQx 13.7 < 19.4 and PhIP 1.2 < 5.1 although the crust is
   3.8 × thicker; with SD 5.7 on 10.9 (line 4), line 4/5 may be one cylinder each with large variability.
4. [!] "Crust thicknesses greater than 3 mm were voluntarily discarded" [V], yet lines 8 (3.85) and 11 (5.82 mm) are
   in Table 3 (line 11 starred only as unacceptable to consumers).
5. [!] Table 1 lists 60 min runs at the two hottest arms; Table 3 has none.
6. [!] Methods: "four 5 mm-thick slices"; Fig. 2 shows five depths (2.5-20 mm).
7. [!] HAA basis is FDP of the crust sample; crust DM and sample mass are not printed, so no fresh-weight or
   per-cylinder amounts can be derived.
8. Text ranges check: category 1 (lines 1, 2, 4, 5) MeIQx < 20, PhIP ≤ 5.1 (text "< 5" [V]; line 2 is 5.1); category 2
   MeIQx 35-41, PhIP 5.8-14.5 (text 35-40 / 5-15): **match** within rounding.

## 6. Existing records in the repo that cite this paper

`grep -rIn -i` for `kondjoyan`, `portanguen`, `foodchem.2016.06.118` over `data/ docs/ src/ scripts/ tasks/`
(2026-11-03). No ledger, registry, prior or code record cites it. Hits:

| record | claim | paper | verdict |
|---|---|---|---|
| `data/keys/pdf_dossier_aliases.yml` l. 101 (pending) | Food Chem 213:641, DOI as above; crust vs HAAs, air-roasting a meat cylinder; "measured T / crust" | crust thickness measured; T measured but **not reported** (jet-T ranges only) | citation **matches**; "measured T" **overstated** |
| `docs/guides/LITERATURE_SEARCH_2026-10_VARIABLES.md` l. 59 | class B; "crust dehydration, surface T and a Maillard-type product (HAA) measured together" | HAAs and crust thickness measured; **no dehydration / water data; surface T not shown** | **partly**: only HAA + crust thickness |
| same, l. 59 context quote | "the cooking process involves chemical transformations, such as Maillard ..." | a citing paper's sentence (column is "context (citing paper)"); not in this PDF | n/a (not a claim about the text) |
| same, l. 43 (V1 verdict) | "coupling to a reaction exists for HAAs in beef (Kondjoyan, Sprague)" | this paper couples **no** transport model to a reaction; the authors say it "would require accurate heat and mass transfer modelling" [V] | **mismatch for Kondjoyan 2016** (true of Kondjoyan 2010b / 2014, not on disk) |
| same, l. 224 | "Next in line ... Kondjoyan 2016 (crust + HAA)" | — | matches |
| same, l. 63 | "Kondjoyan 2009", DOI 10.1016/j.foodchem.2009.02.081 | Crossref 2026-11-03: Food Chem **119**:19-26, issued **2010-03** (this paper's "Kondjoyan et al., 2010a") | **year: 2010** (2009 is the online date); title and DOI match |
| `data/lit/extraction_dossiers/mandin1999_extraction.md` l. 106 | "Kondjoyan & Berdagué 1996" | LRI compilation (QuaPA aroma team), a different work | unrelated to this paper |

## 7. Proposed ledger entries (not written)

The ledgers have no thermal-history lane. HAAs are no lane's species either; `aroma_other.yml`'s scope includes "meat
... odorant and precursor-composition data", the nearest fit. **The two HAA / creatine entries go to
`aroma_other.yml`** (simulable: false; no HAA chemistry in the engine). The crust-thickness series is a transport
observable: no ledger fits; it is proposed under the minimal `thermal_history` field set defined in
`sprague2011_extraction.md` sec. 7 (shared proposal).

```yaml
  # ---- kondjoyan2016 : Kondjoyan et al. 2016 (air-jet roasted beef cylinder, crust HAAs) (ingested 2026-11-03) ----
  - id: kondjoyan2016_tab3_crust_haa_jet_temperature_time
    dossier: kondjoyan2016_extraction.md
    paper: "Kondjoyan, Chevolleau, Portanguen, Molina, Ikonic, Clerjon & Debrauwer 2016 Food Chem 213:641"
    lab: "INRA Clermont QuaPA (Kondjoyan)"
    system: "beef Longissimus thoracis cylinder 48 mm x 50-55 mm, top face under a hot-air jet; crust + <0.5 mm non-brown margin analysed"
    medium: "whole muscle (composition not printed), open air jet, nozzle 36 mm, humidity and velocity not stated"
    temperature_c: [124, 210]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [20, 90]
    species: [MeIQx, IQx, "4,8-DiMeIQx", PhIP, MeIQ, DMIP, Glu-P-2]
    kind: time_series
    location: "Table 3 (11 conditions, 4 jet T x 2-3 times), transcribed in dossier sec. 4a"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "air temperature at the pipe outlet was either 160, 190, 225, or 260 °C"
    flags: [temperature_is_air_jet_not_product, product_T_not_reported, basis_ng_per_g_freeze_dried_crust, crust_mass_not_reported, analogue_internal_standard_TriMeIQx, replicates_analytical_only, text_35x_claim_contradicts_table, sixty_min_runs_missing, "also_relevant: thermal_history"]
    note: "HAAs, not Maillard aroma; the engine has no HAA species. Use only as a thermal-history sanity check (reacting layer far colder than air; HAAs track crust thickness). Apparent Ea vs jet T ~92 kJ/mol [Z], not to be carried."
  - id: kondjoyan2016_fig2_creatine_depth_profile
    dossier: kondjoyan2016_extraction.md
    paper: "Kondjoyan, Chevolleau, Portanguen, Molina, Ikonic, Clerjon & Debrauwer 2016 Food Chem 213:641"
    lab: "INRA Clermont QuaPA (Kondjoyan)"
    system: "same beef cylinder, 192 C jet, 60 min; slices at 2.5-20 mm depth"
    medium: "whole muscle"
    temperature_c: [192]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [60]
    species: [creatine]
    kind: endpoint
    location: "Fig. 2 (5 bars, digitised +-0.03 g/100 g DM, dossier sec. 4b); Table 2 present-study rows (raw, 192 C/20 min, 210 C/90 min)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "Creatine content was measured on meat cylinders exposed to the hot air jet at 192 °C"
    flags: [figure_only_digitised, creatinine_not_measured, slice_count_text_vs_figure, table2_basis_not_stated]
    note: "Precursor depletion toward the surface (1.85 -> 0.70 g/100 g DM in the top 2.5 mm). A depth-resolved precursor gradient, the only one on disk for a roast; meat-only (creatine has no plant analogue)."
```

Thermal-history candidate (field set of `sprague2011_extraction.md` sec. 7):

```yaml
  - id: kondjoyan2016_tab3_crust_thickness_vs_time
    dossier: kondjoyan2016_extraction.md
    paper: "Kondjoyan et al. 2016 Food Chem 213:641"
    lab: "INRA Clermont QuaPA (Kondjoyan)"
    product: "beef LT cylinder 48 mm x 50-55 mm, one face heated"
    heating_mode: air_jet_impingement
    setpoint: "pipe outlet 160/190/225/260 C; jet T 2 mm above surface 120-220 C"
    boundary: {air_t_c: [124, 210], air_velocity_m_s: null, humidity: null, h_w_m2_k: null}
    observable: crust_thickness_mm
    positions: [coloured_crust]
    time_min: [20, 90]
    representation: data_points
    location: "Table 3 crust column (dossier sec. 4a); Fig. 3 isotherm depths at 210 C/70 min [cited, Portanguen 2014]"
    usable_as_programme: false
    status: digitised_in_dossier
    flags: [no_temperature_time_curve, air_T_only, crust_SD_not_printed, linearity_claim_cited_not_shown]
```

## 8. Supplementary material

None referenced (no "Supplementary data" / Appendix in the PDF; searched the text layer for "upplementary",
"Appendix"). The unprinted IR surface-temperature records are the model-relevant missing item; they would be in
Portanguen 2014 / Kondjoyan 2010b, not in SI.
