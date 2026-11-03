# Kovácsné Oroszvári, Bayod, Sjöholm & Tornberg 2006 — EXTRACTION (file stem `oroszvari2006`)
**J. Food Eng. 76(2):169-178**, DOI **10.1016/j.jfoodeng.2005.05.018** (printed p. 169 footer). Crossref 2026-11-03:
"The mechanisms controlling heat and mass transfer on frying of beefburgers. III. Mass transfer evolution during
frying", Bea Kovácsné Oroszvári, Elena Bayod, Ingegerd Sjöholm, Eva Tornberg, issued 2006-09, J. Food Eng. 76(2),
169-178 (online 3 Aug 2005). **Matches the PDF and the brief.** Part 3 of the series whose Part 2 is
`oroszvari2005_extraction.md` (cited here as "2005a"); Part 1 is cited here as **2005b**, J. Food Eng. 67(4):499-506
(the Part 2 dossier calls it "2004"; the issue is 2005) — not on disk, and it holds the patty size, frying rig,
thermocouples and moisture / fat methods this paper only cites.

### Frozen (-20 °C) minced-beef patties from five raw materials (shank 6.7 % fat, rib + added water 12.6 %, lean brisket 15.3 %, rib 16.1 %, fat brisket 39 %; 1.5 % NaCl), D = 10 cm, double-sided contact frying at pan 175 °C, stopped at 30 / 60 / 90 / 120 s or at centre 72 °C (121-151 s); n = 2 per point. T(t) at the centre (5 mm) and 2 mm below the surface for all five meats (Fig. 1); water and fat loss (as % of the initial water / fat) vs time (Figs. 2-3); shrinkage, porosity and fluxes vs time (Figs. 4-9). **No total-loss numbers, no surface or crust temperature, no fitted transfer constant (h = 260 W m⁻² K⁻¹ is cited), no chemistry, no tables of results.**

**Read 2026-11-03 from `data/articles/oroszvari2006.pdf`** (10 pp., Elsevier, born-digital, all plots are vector
graphics). Prose from `pdftotext`; **pp. 4-8 rendered at 110 dpi** (scratchpad `ingest/oroszvari2006_work/p-0N.png`)
to identify the figures; **Figs. 1-3 digitised from the PDF vector paths** (`pdftocairo -svg`, vertices converted
with the axis tick positions; scratch `fig1_pts.txt`, `p5_fig2_raw.txt`; scripts `scratchpad/scripts/svg_*.py`).
Read error: ±0.3 s / ±0.3 °C (Fig. 1 polyline vertices), ±0.2 % on loss values and ±0.3 % on SE bars (Figs. 2-3).
The x-position of every Fig. 2-3 marker equals that meat's own Tm from Fig. 1 to ±0.3 °C, which is how the
(unlabelled, all-diamond) Fig. 2 markers were assigned to meats [Z]. Figs. 4-9 seen, not digitised. Tables 1-2 read
from the render and the text layer.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Department of Food Technology, Engineering and Nutrition, **Lund Institute of Technology** (Lund
University); corresponding author **Bea Kovácsné Oroszvári**, senior author **Eva Tornberg**. **Existing id
`lund_tornberg`** ("Lund (Tornberg)", `data/keys/laboratories.yml`, whose evidence already names the 2006 paper).
Same lab as `oroszvari2005_extraction.md` (Part 2); no other dossier from this group.

## 1. What this paper answers for the model

| model need (mechanism item 8, thermal history) | what Part 3 gives | tag |
|---|---|---|
| T(t) at positions | Fig. 1: centre (5 mm) and 2 mm, five meats, pan 175 °C, 6 time points, mean of n = 2 (no ± shown) | [F] |
| surface / crust T | **not measured**; Fig. 11 schematic prints ~100 °C at 2 mm at the end, crust seen from ~90 s | [M]/[Q] |
| water loss vs time | Fig. 2: WL / initial water (%) per meat at 30-151 s, mean ± SE | [F] |
| fat loss vs time | Fig. 3: FL / initial fat (%) per meat, mean ± SE | [F] |
| total (weighed) loss vs time | **not reported** although Eq. 1 defines it; a + b reconstructed from Figs. 2-3 × Table 1 | [Z] |
| fitted constants | slopes only: water 0.4 % °C⁻¹ of Tm (R² = 0.97); fat 0.6 then 0.2 % °C⁻¹ (two fattest meats) | [M] |
| heat-transfer coefficient | h = 260 W m⁻² K⁻¹ **cited** (Dagerskog & Sörenfors 1978) for a conduction simulation | [C] |
| patty thickness / mass | **not printed** (Part 1); height "no changes" during frying [Q] | [!] |
| shrinkage, porosity, fluxes | Figs. 4-9 (functional; not digitised) | [F] |
| colour / chemistry | none | [!] |

**Answers.**
1. **Measured water / fat / total loss vs frying time: water and fat yes, total no** (confidence ~95 %). Water and
   fat loss are given only normalised to the initial water / fat content, per meat, at 30, 60, 90, 120 s (not for
   rib + water and lean brisket at 120 s) and at t72 (sec. 4b). Only one pan temperature (175 °C) and one diameter
   (10 cm). The weighed total loss (Eq. 1) is never shown; the a + b I reconstruct (sec. 4c) assumes every lost gram is
   water or fat, which is the authors' own Eq. 2-3 assumption.
2. **T(t): yes, better than Part 2** (~95 %): five meats, two depths, means of two patties, 6 points each; the
   meats do not differ significantly [Q]. The hottest measured layer (2 mm) ends at 97-101 °C; nothing at the surface.
3. **Fitted transfer or kinetic constants: none of the paper's own** (~95 %). The 0.4 % °C⁻¹ is a linear slope of
   WL/IW against the mean of two thermocouples, not a transfer coefficient. The only h is cited, and the authors'
   conduction model with it overshoots the measured centre time by 87-144 % (64-106 % with T-dependent k) [M] — their
   argument that contraction-driven flow speeds heating, not a calibrated parameter.
4. **Usable as a heating / drying programme? Partly** (P ≈ 40 % that it enters a programme): Fig. 1 is a
   sub-surface (2 mm) and core programme of a frozen-start beef patty, and Figs. 2-3 give whole-patty water and fat
   loss on the same time grid — a whole-patty drying programme, not a crust moisture history. For a crust pot (where
   the Maillard chemistry is) P ≈ 5 %. Transfer to plant-based patties is weak (frozen start, beef-fat drip,
   protein contraction specific to meat).

## 2. Conditions

- **Meats** [M] (Table 1, n = 3): fat / water % (S.D.): shank 6.7 (0.8) / 71.1 (1.3); rib + added water 12.6 (0.6) /
  68.1 (1.4); lean brisket 15.3 (0.6) / 64.4 (1.3); rib 16.1 (0.5) / 64.3 (1.2); fat brisket (lean brisket minced with
  tallow) 39.0 (1.6) / 46.4 (1.5). Ground 8 then 3 mm plate, 1.5 % NaCl, vacuum packed, frozen. Water added to rib:
  amount not printed.
- **Patties** [M]: D = 10 cm, prepared "as described in" Part 1 [V]; thickness and mass not printed [!]; fried from
  the frozen state (Fig. 1 starts at -20 °C).
- **Design** (Table 2) [M]: 175 °C pan; stop at 30, 60, 90, 120 s and at centre 72 °C; 120 s omitted for rib + water and
  lean brisket. **t72 = 146 / 121 / 124 / 146 / 151 s** (shank / rib + water / lean brisket / rib / fat brisket). Each
  point "repeated two times" [V].
- **Heating** [M]: double-sided pan fryer, set-up of Part 1 (not described here); no oil stated; no contact pressure.
- **Measurements** [M]: weight frozen before and after frying; diameter at 8 points, 5 min after frying; moisture and
  fat after frying (Part 1 methods); density (Part 2 method). pH, aw, O2: not measured.

## 3. Analytical methods

- Temperatures: thermocouples, acquisition as in Part 1; centre "(5 mm)" and 2 mm below the surface;
  Tm = (T_centre + T_2mm)/2 [M].
- Losses (Eqs. 1-3) [M]: TL = (Wb − Wa)/Wb × 100; water loss a and fat loss b (% of raw patty) from the fried water and
  fat contents via Wf = (Wi − a)/(100 − a − b) × 100 and Ff = (Fi − b)/(100 − a − b) × 100. Plotted as a/Wi (WL/IW)
  and b/Fi (FL/IF).
- Fluxes (Eqs. 4-8): total and water flux per total area A_T = 2πr² + 2πrh, fat flux per lateral area 2πrh, divided
  by frying time [M]. Porosity from shrinkage (Eq. 10) and from densities (Eqs. 17-18, Choi & Okos), pure-component
  densities 995.4 / 916.4 / 1318.5 kg m⁻³ (water / fat / protein, 22 °C) [M].
- ± in Figs. 2-9: standard errors of n = 2 [M].

## 4. Data

### 4a. Fig. 1, T (°C) at pan 175 °C, mean of 2 patties [F, vector vertices, ±0.3 °C]

| t (s) | S 5 mm | RW 5 mm | LB 5 mm | R 5 mm | FB 5 mm | S 2 mm | RW 2 mm | LB 2 mm | R 2 mm | FB 2 mm |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | -20.0 | -20.0 | -20.0 | -20.0 | -20.0 | -20.0 | -20.0 | -20.0 | -20.0 | -20.0 |
| 30 | -11.1 | -10.7 | -12.6 | -12.1 | -15.4 | 29.6 | 34.7 | 28.6 | 29.1 | 19.4 |
| 60 | -5.2 | -4.3 | -5.6 | -6.1 | -7.9 | 56.9 | 68.1 | 71.4 | 62.6 | 55.6 |
| 90 | 3.7 | 15.7 | 12.9 | 11.1 | 3.1 | 74.1 | 92.7 | 91.7 | 82.9 | 78.3 |
| 120 | 47.6 | — | — | 47.2 | 40.8 | 86.6 | — | — | 93.6 | 91.3 |
| t72 | 72.3 (146) | 71.8 (121) | 72.3 (124) | 72.7 (146) | 72.3 (151) | 97.3 | 101.0 | 100.5 | 100.1 | 99.1 |

S shank, RW rib + water, LB lean brisket, R rib, FB fat brisket. Points are joined by straight dashed lines (no
intermediate readings). The Tm line is the **five-meat mean**: -20 / 7.9 / 29.1 / 47.2 / 75.0 / 86.2 °C at
0 / 30 / 60 / 90 / 120 / 137 s (137 = mean t72) [F]; it equals the mean of the ten series above to within 0.6 °C,
with the 121 / 124 s end points of rib + water and lean brisket pooled into the 120 s mean (75.3 vs 75.0) [Z]. Two
bold guide lines on the 2 mm traces (slopes ≈ 1.38 and 0.27 °C s⁻¹ [Z]) cross at the printed **70 s / 79 °C**
("increases steeply for about 70 s" [Q]). Fig. 11 schematic prints rounded mean temperatures (2 mm / centre) [M]:
30 s ~34 / -12; 60 s ~57 / -5; 90 s ~82 / 8; 120 s ~90 / 48; end ~100 / 72 °C; crust drawn from 90 s; whole patty
thawed between 60 and 90 s [Q].

### 4b. Figs. 2-3, water loss / initial water and fat loss / initial fat (%), mean ± SE (n = 2) [F, vector, ±0.2 %]

| meat | t (s) | Tm (°C) | WL/IW | FL/IF |
|---|---|---|---|---|
| shank | 30 / 60 / 90 / 120 / 146 | 9.1 / 26.0 / 39.0 / 67.1 / 84.6 | 12.9 ± 2.1 / 18.0 ± 3.2 / 22.3 ± 0.8 / 30.7 ± 3.6 / 43.7 ± 0.9 | 2.1 ± 0.9 / 6.5 ± 4.0 / 7.0 ± 4.9 / 12.9 ± 4.1 / 12.4 ± 0.5 |
| rib + water | 30 / 60 / 90 / 121 | 12.1 / 32.0 / 54.4 / 86.5 | 5.4 ± 1.4 / 22.3 ± 0.5 / 29.9 ± 3.1 / 42.9 ± 0.5 | 0.9 ± 0.9 / 11.0 ± 0.9 / 29.0 ± 6.0 / 36.3 ± 2.0 |
| lean brisket | 30 / 60 / 90 / 124 | 8.2 / 32.6 / 52.3 / 86.5 | 11.8 ± 0.2 / 18.7 ± 1.4 / 25.5 ± 2.6 / 40.2 ± 1.0 | 2.1 ± 2.2 / 10.8 ± 1.4 / 38.4 ± 10.0 / 41.4 ± 0.9 |
| rib | 30 / 60 / 90 / 120 / 146 | 8.7 / 28.1 / 47.1 / 70.4 / 86.5 | 9.0 ± 3.6 / 14.4 ± 0.2 / 31.3 ± 1.5 / 35.9 ± 2.7 / 46.0 ± 0.3 | 18.7 ± 4.6 / 24.1 ± 1.6 / 46.1 ± 3.3 / 47.1 ± 5.1 / 50.8 ± 1.6 |
| fat brisket | 30 / 60 / 90 / 120 / 151 | 2.1 / 23.9 / 41.1 / 66.2 / 85.8 | 9.3 ± 4.9 / 17.7 ± 4.1 / 21.0 ± 2.1 / 32.3 ± 5.0 / 44.4 ± 2.1 | 21.1 ± 8.3 / 33.5 ± 6.1 / 51.1 ± 5.0 / 58.8 ± 7.0 / 62.1 ± 8.9 |

Both figures also draw a (0 %, Tm = -20 °C, t = 0) origin point, which is assumed, not measured [!]. Fig. 2 fitted line
[F]: through (−20, 0) and (87, 42.1), slope 0.39 % °C⁻¹, matching the printed **0.4 % °C⁻¹, R² = 0.97** [M]. Text:
"approximately 40-45 %" of initial water lost at the end [V/M] (Fig. 2: 40.2-46.0 ✓); shank fat loss "5-10 %"
with almost no rise [Q] (Fig. 3: 2-13 %).

### 4c. Reconstructed loss as % of the raw patty [Z]: a = WL/IW × Wi/100, b = FL/IF × Fi/100 (Table 1)

| meat | 30 s a / b / a+b | 60 s | 90 s | 120 s | t72 |
|---|---|---|---|---|---|
| shank | 9.2 / 0.1 / 9.3 | 12.8 / 0.4 / 13.2 | 15.9 / 0.5 / 16.3 | 21.8 / 0.9 / 22.7 | 31.1 / 0.8 / 31.9 |
| rib + water | 3.7 / 0.1 / 3.8 | 15.2 / 1.4 / 16.6 | 20.4 / 3.7 / 24.0 | — | 29.2 / 4.6 / 33.8 |
| lean brisket | 7.6 / 0.3 / 7.9 | 12.0 / 1.7 / 13.7 | 16.4 / 5.9 / 22.3 | — | 25.9 / 6.3 / 32.2 |
| rib | 5.8 / 3.0 / 8.8 | 9.3 / 3.9 / 13.1 | 20.1 / 7.4 / 27.6 | 23.1 / 7.6 / 30.7 | 29.6 / 8.2 / 37.8 |
| fat brisket | 4.3 / 8.2 / 12.5 | 8.2 / 13.1 / 21.3 | 9.7 / 19.9 / 29.7 | 15.0 / 22.9 / 37.9 | 20.6 / 24.2 / 44.8 |

Text check [M]: fat brisket fat loss "about 24 %" of the whole burger ✓ (24.2), water "about 18 %" (20.6, 2.6 points
higher); "next fattiest" water 24 % / fat 8 % — rib gives 29.6 / 8.2, lean brisket 25.9 / 6.3 [!] (sec. 5).

### 4d. Figs. 4-9 (not digitised)
Fig. 4 diameter shrinkage vs Tm (fat brisket ~22 %, others ~12 % at the end [M]); Figs. 5-6 density- and
shrinkage-based porosity vs 2 mm T (from 1-4 % to ~35 %, steep above ~90 °C [M]); Fig. 7 total flux vs centre T
(with Pan et al. 2000 for comparison); Figs. 8-9 water and fat flux (g cm⁻² s⁻¹) vs time. Fluxes are the same losses
divided by area and time, so they add nothing independent without the (unprinted) patty height.

## 5. Flags and internal consistency

1. [!] **Part 2 vs Part 3 at the same nominal condition** (fat brisket, D 10 cm, pan 175 °C, frozen start): t72 = 99 s
   in Part 2 Fig. 1 vs **151 s** here; 2 mm T at 30 / 60 / 90 s = 36 / 75 / 93 °C (Part 2) vs 19 / 56 / 78 °C (here).
   Water loss/IW at the end agrees (45.5 vs 44.4 %). Part 2's own text gives 132 s as the shortest frying time, so the
   Part 2 Fig. 1 trace is the outlier (P ≈ 70 %); a thinner patty in that run is a candidate cause.
2. [!] **Total loss (Eq. 1) defined, measured (weights) but never reported**; a + b from Eqs. 2-3 is not compared with
   TL, so mass closure (protein / salt in the drip) is unchecked.
3. [!] Text "next fattiest" water / fat losses 24 / 8 % do not match the figures for either candidate meat (rib 29.6 /
   8.2, lean brisket 25.9 / 6.3) [Z]; fat brisket water "about 18 %" vs 20.6 [Z].
4. [!] Figs. 2-9 x-axes are each meat's own Tm (matches Fig. 1 per meat), but the top "frying time" / Tm labels
   (30 … 137 s; 8, 29, 47, 75, 86 °C) are five-meat means; points of one time are spread over ±6 °C of Tm.
5. [!] SE of n = 2 (half the range); several bars are 0.2-0.5 % (shank 146 s, rib 60 s) next to ±8-10 % elsewhere.
6. [!] Patty thickness, mass, thermocouple type, logging interval and fryer details all deferred to Part 1 (not on
   disk); the centre probe at "5 mm" implies ~10 mm thickness [Z, unverified].
7. [!] Fig. 9 legend labels are truncated / duplicated ("Shan", "Brisket" twice) in the PDF.

## 6. Existing records in the repo that cite this paper

`grep -rIn` 2026-11-03 for `jfoodeng.2005.05.018`, `oroszvari2006`, "Mass transfer evolution", "Kovácsné … 2006"
across `data/ docs/ src/ scripts/ tasks/`: four hits, no code, no ledger entry of its own.

| record | field / claim | repo value | paper | verdict |
|---|---|---|---|---|
| `docs/guides/LITERATURE_SEARCH_2026-10_VARIABLES.md` l. 86 | title, journal, DOI | Part III title, J Food Eng, 10.1016/j.jfoodeng.2005.05.018 | same | **match** |
| same | what it gives | "water and fat loss vs time during frying" | Figs. 2-3, one pan T, five meats | **match** |
| same | quote | "Water losing of beef burger was about 80% from during pan-frying" | no 80 % statement anywhere in Part 3; the "about 80 %" drip share is Part 2's (`oroszvari2005_extraction.md` sec. 1) | **mismatch: misattributed quote** |
| `data/lit/measurement_ledger/thermal_history.yml` l. 134 (note on `oroszvari2005_tab3_tab4_figs2_4_endpoint_mass_loss`) | "Time-resolved water/fat loss is in Part 3 (…), not on disk" | time-resolved water/fat loss: yes, Figs. 2-3; DOI matches | **content match; "not on disk" now stale** |
| `oroszvari2005_extraction.md` ll. 7-10, 269 | Part 3 title, DOI, verdict | as above | match; its Part 1 year "2004" vs 2005b, 67(4):499-506 cited here | **match** (Part 1 year minor) |
| `data/keys/laboratories.yml` l. 98 (`lund_tornberg`) | evidence "Kovacsne Oroszvari et al. 2005 / 2006: Eva Tornberg, Food Technology, Lund University" | Dept. Food Technology, Engineering and Nutrition, Lund Institute of Technology; Tornberg senior | **match** |

## 7. Proposed ledger entries (`data/lit/measurement_ledger/thermal_history.yml`; not written)

```yaml
  # ---- oroszvari2006 ----
  - id: oroszvari2006_fig1_subsurface_centre_t_vs_time_five_meats
    dossier: oroszvari2006_extraction.md
    paper: "Kovacsne Oroszvari, Bayod, Sjoholm & Tornberg 2006 J Food Eng 76:169"
    lab: "Lund (Tornberg)"
    product: "frozen minced-beef patties, 5 raw materials (fat 6.7-39 %, water 46-71 %), 1.5 % NaCl, D = 10 cm (thickness not printed)"
    heating_mode: double_sided_contact_frying
    setpoint: "pan 175 C, stop at 30/60/90/120 s or centre 72 C (t72 121-151 s)"
    boundary: {pan_t_c: [175], contact_pressure: null, h_w_m2_k: null, initial_t_c: -20}
    observable: temperature_c
    positions: [centre_5mm, subsurface_2mm]
    time_min: [0, 2.5]
    representation: data_points
    location: "Fig. 1 (10 series x 6 points, mean of n = 2, no SE), digitised from vector paths in dossier sec. 4a; Fig. 11 rounded means"
    usable_as_programme: partial
    status: digitised_in_dossier
    flags: [frozen_start, no_surface_or_crust_temperature, thickness_not_printed, points_joined_by_straight_lines, contradicts_part2_fig1_same_condition]
    note: "2 mm layer ends at 97-101 C; meats not significantly different. Fat brisket t72 151 s here vs 99 s in Part 2 Fig. 1 (same nominal condition)."
  - id: oroszvari2006_figs2_3_water_fat_loss_vs_time
    dossier: oroszvari2006_extraction.md
    paper: "Kovacsne Oroszvari, Bayod, Sjoholm & Tornberg 2006 J Food Eng 76:169"
    lab: "Lund (Tornberg)"
    product: "as oroszvari2006_fig1_subsurface_centre_t_vs_time_five_meats"
    heating_mode: double_sided_contact_frying
    setpoint: "pan 175 C, stop at 30/60/90/120 s or centre 72 C"
    boundary: {pan_t_c: [175], contact_pressure: null, h_w_m2_k: null, initial_t_c: -20}
    observable: [water_loss_pct_of_initial_water, fat_loss_pct_of_initial_fat]
    positions: [whole_patty]
    time_min: [0.5, 2.5]
    representation: data_points
    location: "Figs. 2-3 (5 meats x 4-5 times, mean +- SE, n = 2), digitised in dossier sec. 4b; % of raw patty reconstructed in sec. 4c [derived]"
    usable_as_programme: partial
    status: digitised_in_dossier
    flags: [frozen_start, total_loss_not_reported, normalised_losses_only, origin_point_assumed, x_axis_is_mean_of_two_thermocouples, single_pan_temperature]
    note: "Whole-patty drying programme: WL/IW rises ~linearly to 40-46 % at t72 (0.4 % per C of Tm, R2 0.97); fat loss saturates for the two fattest meats. No local (crust) moisture."
  - id: oroszvari2006_figs4_9_shrinkage_porosity_flux_vs_time
    dossier: oroszvari2006_extraction.md
    paper: "Kovacsne Oroszvari, Bayod, Sjoholm & Tornberg 2006 J Food Eng 76:169"
    lab: "Lund (Tornberg)"
    product: "as oroszvari2006_fig1_subsurface_centre_t_vs_time_five_meats"
    heating_mode: double_sided_contact_frying
    setpoint: "pan 175 C"
    boundary: {pan_t_c: [175], contact_pressure: null, h_w_m2_k: null, initial_t_c: -20}
    observable: [shrinkage_pct, porosity, total_flux_g_cm2_s, water_flux_g_cm2_s, fat_flux_g_cm2_s]
    positions: [whole_patty]
    time_min: [0.5, 2.5]
    representation: data_points
    location: "Figs. 4-9 (not digitised)"
    usable_as_programme: false
    status: figure_only
    flags: [functional_property, derived_from_same_losses, height_not_printed]
```

## 8. Supplementary material
None referenced (no Appendix, no supplementary data; Elsevier 2006). Everything model-relevant is in Figs. 1-3 and
Tables 1-2. The missing pieces (patty thickness / mass, rig, thermocouple placement) are in **Part 1**, Kovácsné
Oroszvári, Sjöholm & Tornberg 2005, J. Food Eng. 67(4):499-506 — not on disk; worth fetching only if a
heat-transfer model of these patties is built.
