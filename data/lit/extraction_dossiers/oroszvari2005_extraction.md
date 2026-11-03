# Kovácsné Oroszvári, Bayod, Sjöholm & Tornberg 2005 — EXTRACTION (file stem `oroszvari2005`)
DOI: 10.1016/j.jfoodeng.2004.10.013 (printed p. 18 footer). Crossref 2026-11-03: "The mechanisms controlling heat and
mass transfer on frying of beefburgers. Part 2: The influence of the pan temperature and patty diameter", J. Food
Eng. **71**(1):18-27, authors Bea Kovácsné Oroszvári, Elena Bayod, Ingegerd Sjöholm, Eva Tornberg, issued 2005-11
(online 22 Dec 2004). **Matches the PDF and the brief.**
[!] **This is Part 2 of a three-part series.** `docs/guides/LITERATURE_SEARCH_2026-10_VARIABLES.md` does **not** list
this paper; it lists **Part 3** (Kovácsné Oroszvári 2006, "III. Mass transfer evolution during frying",
10.1016/j.jfoodeng.2005.05.018, verdict MAYBE, "water and fat loss vs time"). The time-resolved mass loss the search
wanted is in Part 3, not here. Part 1 (Kovácsné, Sjöholm & Tornberg 2004 J. Food Eng. 67:499, composition and
comminution) holds the frying set-up, thermocouple placement and analytical methods this paper only cites. Neither is
on disk.

### Frozen (-20 °C) minced-beef patties (fat 6.7 / 16.1 / 39 %, 1.5 % NaCl), D = 3-12.7 cm, double-sided contact frying at pan 100 / 150 / 175 °C until the centre reaches 72 °C (90-190 s). Measured T(t) at the centre (5 mm) and 2 mm below the surface for ONE patty type (Fig. 1); end-of-frying total / water / fat loss, fluxes, porosity, shrinkage (Tables 3-4, Figs. 2-6). **No surface or crust temperature, no colour, no mass loss vs time, no heat-transfer coefficient, no chemistry.**

**Read 2026-11-03 from `data/articles/oroszvari2005.pdf`** (10 pp., Elsevier, born-digital, figures are vector
graphics). Prose from `pdftotext -layout`; **pp. 4-9 rendered at 150 dpi** and Figs. 1-3 at 300 dpi (scratchpad
`ingest/oroszvari2005_work/p-0N.png`, `fig1hi-04.png`, `fig2hi-07.png`, `fig3hi-07.png`). **Table 3 (rotated, p. 6) read
from the render** and checked cell by cell against the text layer; Table 4 (p. 9) text layer checked against render.
**Figs. 1-4 digitised from the PDF vector paths** (`pdftocairo -svg`, path coordinates converted with the axis tick
positions; scratch `p4_shapes.txt`, `p7_shapes.txt`): read error is the line width, about ±0.3 °C / ±0.3 s for Fig. 1
lines, ±1.5 °C for the 175 °C diamond series (10 s bin means of the marker vertices), ±0.1 % / ±0.05 °C /
±0.00001 g cm⁻² s⁻¹ for the bar heights of Figs. 2-4. Figs. 5-6 (porosity) not digitised.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Department of Food Technology, Engineering and Nutrition, **Lund Institute of Technology** (Lund
University), Sweden; corresponding author **Bea Kovácsné Oroszvári**, senior author **Eva Tornberg** (Tornberg's
meat-science group; Sjöholm co-author). **Not in `data/keys/laboratories.yml`: NEW**, proposed id **`lund_tornberg`**,
lab string "Lund (Tornberg)". No existing dossier from this group (grep of `extraction_dossiers/` for Tornberg,
Kovácsné, Oroszvári, Lund: none).

## 1. What this paper answers for the model

| model need (mechanism item 8, thermal history) | what Kovácsné Oroszvári 2005 gives | tag |
|---|---|---|
| T(t) at positions in a frying patty | Fig. 1: centre (5 mm) and 2 mm sub-surface, pan 100 / 150 / 175 °C, fat brisket, D = 10 cm, ONE run each (frozen start) | [F] |
| surface / crust T(t) | **not measured**; only pan setpoint and 2 mm T. 2 mm end T 84-102 °C (Fig. 3) | [!] |
| end-of-frying 2 mm T by meat × pan T × D | Fig. 3 (18 bars ± SE), Tables 3-4 main effects | [F]/[M] |
| weight / water / fat loss vs time | **not given**: only at the end of frying (centre 72 °C); time courses are in Part 3 (2006) | [!] |
| end-point total / water / fat loss, fluxes | Table 3 (pan T, D, meat main effects), Table 4 (rib, 5 diameters), Figs. 2, 4 | [M]/[F] |
| evaporation vs drip split | "about 80 %" of water loss is drip even at 175 °C [Q], from 100 °C-pan water loss / 175 °C-pan water loss | [Q]/[Z] |
| heat-transfer coefficient h, contact resistance, evaporation front | **none** (no model in this paper); crust forms above ~102 °C [Q] | [!] |
| surface colour / browning, L*a*b* | **none** (crust judged by eye only) | [!] |
| composition | Table 1 fat / water of four raw materials | [M] |

**Answers.**
1. **Measured T(t) at positions: yes, but thin** (confidence ~95 %). Fig. 1 is one patty type (fat brisket, 39 % fat,
   D = 10 cm) at three pan temperatures, two depths, from -20 °C to the 72 °C centre endpoint, one trace per condition
   (no replicate traces, no ± ). Other cells have only the final 2 mm temperature (Fig. 3).
2. **Is any of it a Maillard-relevant temperature? No** (~90 %). The hottest measured point is 2 mm below the surface
   and ends at 84-102 °C (Fig. 3); the 2 mm trace at 175 °C passes 90 °C only after ~85 s. Maillard chemistry in a
   fried patty runs in the crust (top < 1 mm, above 100 °C, between ~100 °C and the 150-175 °C pan); this paper
   measures neither the crust nor the surface. Pan T is an upper bound only.
3. **Weight / water / fat loss vs time: no** (~99 %). Every loss is at t72 (centre 72 °C). Fluxes (Eqs. 4-5) are
   losses divided by area and by t72, i.e. time-averaged.
4. **Fitted constants (h, k, Ea, evaporation front): none** (~99 %). The paper is an ANOVA of a factorial design.
5. **Usable as a heating-programme input?** Only for an interior (≤ 100 °C) layer of a frozen-start patty
   (P ≈ 30 % that anything here enters a programme; it would be the 2 mm trace as a "sub-crust" example). Not for a
   crust pot (P ≈ 5 %). For plant-based patties the frozen start (latent heat of thawing dominates 0-70 s at the centre)
   and beef fat drip make the transfer weak.
6. **Usable as a browning check? No** (~99 %): no colour or chemistry measured.

## 2. Conditions

- **Meat** [M]: beef fat brisket (FB; lean brisket minced with tallow), lean brisket (LB), rib (R), shank (S); ground 8 mm
  then 3 mm plate; **1.5 % NaCl**, hand-mixed 5 min; vacuum packed and frozen. Composition Table 1 (sec. 4a).
- **Patties** [M]: preparation, thickness and mass "as described in Kovácsné et al. (2004)" — **not printed here** [!].
  Centre thermocouple at "5 mm" implies ~10 mm thickness [Z, unverified]. Fried **from the frozen state**; Fig. 1 starts
  at about -20 °C [F].
- **Design** (Table 2) [M]: Exp. Ia FB/LB/S × pan 100/150/175 °C × D 10 cm; Ib the same at D 3 cm; Exp. II rib × pan
  150/175 °C × D 3.0/5.3/8.1/10.0/12.7 cm. "Full experimental design with two replicates" [V].
- **Heating** [M]: double-sided contact frying (set-up of Part 1, not described); stop when centre = 72 °C. Pan surface
  drops when a large patty is laid on [Q]. No oil stated, no surface T, no contact pressure, no h.
- **pH, aw, O2**: not measured / not relevant to this paper.

## 3. Analytical methods

- Temperatures: thermocouples, data acquisition "as described in Kovácsné et al. (2004)" [M]; positions centre (5 mm)
  and 2 mm below the surface. Thermocouple type, logging interval and calibration not printed [!].
- Weights before frying (frozen) and after; total loss TL = (Wb − Wa)/Wb (Eq. 1). Water and fat after frying by the
  Part 1 methods; water and fat loss from mass balance (Eqs. 2-3). Fluxes Wr, Fr in g cm⁻² s⁻¹ (Eqs. 4-5; water over
  total area AT, fat over lateral area AL only). Shrinkage from 8 diameters; porosity εs (volume balance) and εd
  (density by rapeseed-oil displacement, ρ 914 kg m⁻³; component densities water 995.4, fat 916.4, protein
  1318.5 kg m⁻³ at 22 °C, Choi & Okos) [M].
- Statistics: Minitab 12.2; tables are **least-square means** of 3-way (Table 3) / 2-way (Table 4) ANOVA; figures
  "mean values ... with standard errors" [V] (n = 2 per cell).

## 4. Data

### 4a. Table 1, raw-material composition (% w/w, ± "Std") [M]

| sample | fat | water |
|---|---|---|
| fat brisket (FB) | 39.0 ± 1.6 | 46.4 ± 1.5 |
| rib (R) | 16.1 ± 0.5 | 64.3 ± 1.2 |
| lean brisket (LB) | 15.3 ± 0.6 | 64.4 ± 1.3 |
| shank (S) | 6.7 ± 0.8 | 71.1 ± 1.3 |

### 4b. Fig. 1, measured T(t) (°C), fat brisket, D = 10 cm, one trace each [F, vector paths]

| t (s) | 100 °C, 2 mm | 150 °C, 2 mm | 175 °C, 2 mm | 100 °C, centre | 150 °C, centre | 175 °C, centre |
|---|---|---|---|---|---|---|
| 0 | -16.5 | -19.6 | -18 | -19.6 | -20.7 | -20.1 |
| 10 | -13.5 | -11.3 | -9.4 | -19.1 | -19.6 | -19.6 |
| 20 | -6.7 | 14.9 | 13.7 | -18.1 | -18.1 | -17.6 |
| 30 | 5.7 | 38.2 | 36.1 | -16.0 | -15.5 | -14.9 |
| 40 | 21.7 | 49.6 | 50.2 | -14.0 | -12.4 | -11.3 |
| 50 | 35.1 | 57.8 | 60.1 | -11.3 | -9.3 | -8.3 |
| 60 | 46.4 | 64.0 | 74.9 | -9.3 | -6.2 | -5.7 |
| 70 | 52.1 | 68.7 | 84.6 | -7.2 | -3.6 | -1.6 |
| 80 | 55.7 | 74.3 | 88.6 | -5.7 | 4.7 | 32.5 |
| 90 | 58.8 | 82.5 | 93.3 | -4.1 | 22.7 | 55.2 |
| 100 | 61.9 | 89.8 | 96.2 (end ≈ 100 s) | -1.6 | 37.7 | 72.8 (end 99 s) |
| 110 | 65.1 | 94.4 | — | 7.7 | 52.7 | — |
| 120 | 68.1 | 95.5 (97.6 at end, 124 s) | — | 19.6 | 67.6 (73.3 at end, 124 s) | — |
| 130 / 140 / 150 | 71.2 / 74.8 / 78.4 | — | — | 28.4 / 36.7 / 45.4 | — | — |
| 160 / 170 / 180 / 190 | 81.6 / 84.7 / 88.3 / 90.8 | — | — | 53.7 / 60.4 / 67.6 / 72.3 | — | — |

Centre crosses 0 °C at ~70 / 78 / 102 s (175 / 150 / 100 °C) [F]; no flat thawing plateau is visible at the centre
(the trace rises steadily from -20 to ~-2 °C, then steeply) [F]. The 175 °C 2 mm series is markers only (± 1.5 °C).

### 4c. Fig. 3, final T 2 mm below the surface (°C), mean of 2, ± SE not transcribed [F, vector]

| pan T | D | shank | lean brisket | fat brisket |
|---|---|---|---|---|
| 100 | 3 cm | 87.1 | 88.2 | 90.8 |
| 100 | 10 cm | 91.1 | 83.7 | 89.0 |
| 150 | 3 cm | 99.6 | 100.2 | 99.5 |
| 150 | 10 cm | 94.8 | 93.3 | 99.8 |
| 175 | 3 cm | 101.6 | 101.1 | 100.6 |
| 175 | 10 cm | 91.1 | 96.2 | 97.9 |

Text: lowest 83.7 °C (LB, D 10, 100 °C), highest 101.6 °C (shank, D 3, 175 °C) [M] — both match the bars.

### 4d. Fig. 2, water loss / initial water content (%), and Fig. 4 fluxes (g cm⁻² s⁻¹), mean of 2 [F, vector]

| pan T | D | water loss/initial water: S / LB / FB | water flux ×10⁴: S / LB / FB | fat flux ×10⁴: S / LB / FB |
|---|---|---|---|---|
| 100 | 3 | 40.0 / 41.3 / 37.0 | 6.5 / 6.1 / 4.1 | 0.7 / 4.4 / 17.5 |
| 100 | 10 | 33.7 / 34.6 / 31.7 | 5.6 / 5.5 / 3.6 | 2.0 / 7.2 / 38.5 |
| 150 | 3 | 48.1 / 53.8 / 49.5 | 8.7 / 11.4 / 7.4 | 0.7 / 7.4 / 24.4 |
| 150 | 10 | 38.7 / 34.1 / 41.8 | 7.9 / 7.2 / 6.2 | 3.0 / 7.2 / 50.9 |
| 175 | 3 | 50.9 / 57.0 / 54.3 | 11.9 / 13.1 / 8.1 | 1.5 / 6.9 / 24.4 |
| 175 | 10 | 39.6 / 38.8 / 45.5 | 9.0 / 8.8 / 8.0 | 1.5 / 9.6 / 60.2 |

Text check: at 100 °C "about 33 %" (D 10) and "about 39 %" (D 3) [M] → bar means 33.3 / 39.4 [Z] ✓.

### 4e. Table 3, main effects (least-square means; 3-way ANOVA; Exp. Ia + Ib) [M]

| quantity | pan 100 | 150 | 175 | p | D 3.0 | D 10.0 | p | shank | LB | FB | p |
|---|---|---|---|---|---|---|---|---|---|---|---|
| thawing time, centre (s) | 84.3 | 72.4 | 67.2 | <0.001 | 63.4 | 85.9 | <0.001 | 70.6 | 77.5 | 75.8 | 0.041 |
| total frying time (s) | 167.6 | 135.7 | 121.0 | <0.001 | 118.0 | 164.9 | <0.001 | 138.4 | 141.6 | 144.3 | 0.560 |
| final T 2 mm (°C) | 88.5 | 97.8 | 98.3 | <0.001 | 96.4 | 93.3 | 0.003 | 94.5 | 94.0 | 96.1 | 0.179 |
| total loss (%) | 33.6 | 38.1 | 39.2 | 0.001 | 41.0 | 32.9 | <0.001 | 28.3 | 34.2 | 48.3 | <0.001 |
| fat loss (%) | 11.9 | 12.2 | 13.0 | 0.694 | 13.3 | 11.4 | 0.092 | 2.5 | 6.3 | 28.3 | <0.001 |
| water loss (%) | 22.0 | 26.7 | 27.5 | 0.001 | 28.9 | 21.8 | <0.001 | 28.5 | 27.9 | 19.7 | <0.001 |
| fat loss / initial fat (%) | 45.6 | 46.2 | 59.6 | 0.464 | 51.3 | 49.6 | 0.871 | 37.4 | 41.1 | 72.9 | 0.012 |
| water loss / initial water (%) | 36.0 | 44.2 | 46.1 | <0.001 | 47.7 | 36.4 | <0.001 | 40.1 | 43.3 | 42.9 | 0.316 |
| fat flux (g cm⁻² s⁻¹ ×10⁴) | 11.2 | 15.2 | 17.6 | 0.280 | 8.8 | 20.5 | 0.001 | 2.5 | 7.1 | 34.3 | <0.001 |
| water flux (g cm⁻² s⁻¹ ×10⁴) | 5.2 | 8.1 | 9.8 | <0.001 | 8.5 | 6.9 | 0.001 | 8.3 | 8.7 | 6.1 | <0.001 |
| porosity εs | 0.188 | 0.248 | 0.228 | 0.021 | 0.297 | 0.145 | <0.001 | 0.147 | 0.218 | 0.299 | <0.001 |
| diameter shrinkage (%) | 13.6 | 13.9 | 15.0 | 0.560 | 13.2 | 15.1 | 0.090 | 11.2 | 12.5 | 18.8 | <0.001 |

### 4f. Table 4, rib meat, main effects (2-way ANOVA; Exp. II) [M]

| quantity | pan 150 | 175 | p | D 3.0 | 5.3 | 8.1 | 10.0 | 12.7 | p |
|---|---|---|---|---|---|---|---|---|---|
| thawing time, centre (s) | 73.1 | 67.7 | 0.278 | 57.5 | 60.1 | 70.6 | 82.8 | 80.9 | 0.011 |
| total frying time (s) | 134.6 | 118.9 | 0.146 | 98.5 | 101.1 | 120.5 | 148.3 | 165.2 | 0.002 |
| Tf 2 mm (°C) | 93.8 | 95.9 | 0.463 | 99.4 | 89.8 | 93.2 | 94.1 | 97.9 | 0.233 |
| total loss (%) | 34.1 | 38.6 | 0.030 | 37.6 | 41.0 | 38.7 | 33.4 | 31.0 | 0.021 |
| fat loss (%) | 7.4 | 7.8 | 0.688 | 7.3 | 8.7 | 8.8 | 7.5 | 5.7 | 0.233 |
| water loss (%) | 25.3 | 28.0 | 0.056 | 28.9 | 26.5 | 27.7 | 26.3 | 23.9 | 0.202 |
| fat loss / initial fat (%) | 45.8 | 48.2 | 0.688 | 45.0 | 53.8 | 54.5 | 46.3 | 35.4 | 0.233 |
| water loss / initial water (%) | 39.3 | 43.6 | 0.056 | 44.9 | 41.2 | 43.1 | 40.9 | 37.2 | 0.202 |
| fat flux (×10⁴) | 12.3 | 15.0 | 0.249 | 7.1 | 14.3 | 17.8 | 15.1 | 13.8 | 0.092 |
| water flux (×10⁴) | 8.9 | 11.2 | 0.010 | 10.8 | 11.7 | 11.2 | 8.7 | 7.7 | 0.019 |
| porosity εd | 0.071 | 0.147 | 0.028 | – | – | 0.115 | 0.143 | 0.068 | 0.148 |
| porosity εs | 0.244 | 0.289 | 0.059 | 0.348 | 0.275 | 0.300 | 0.296 | 0.115 | <0.001 |
| diameter shrinkage (%) | 11.2 | 9.7 | 0.185 | 6.5 | 10.0 | 10.3 | 12.0 | 13.4 | 0.013 |

Figs. 5 (εs by meat × pan T × D) and 6 (εs, εd vs D at 150 / 175 °C) — porosity, functional; not digitised.

### 4g. Derived [Z]

- Drip share: water loss/initial at pan 100 / 175 = 36.0 / 46.1 = **0.78** (Table 3) → the authors' "about 80 %"
  reproduces, *given* their premise that a 100 °C pan evaporates nothing.
- Evaporative water flux attributable to the hotter pan: (9.8 − 5.2)×10⁻⁴ = 4.6×10⁻⁴ g cm⁻² s⁻¹ ≈ 4.6 g m⁻² s⁻¹
  averaged over the whole patty area and t72 (Table 3).
- Mass balance checks: fat + water loss vs total loss (Table 3) 33.9 / 38.9 / 40.5 vs 33.6 / 38.1 / 39.2 (pan T);
  FB 48.0 vs 48.3; LB 34.2 vs 34.2; **shank 31.0 vs 28.3** (gap 2.7 points). Ratios to Table 1 composition reproduce
  (FB fat 28.3/39.0 = 72.6 vs 72.9; shank water 28.5/71.1 = 40.1 vs 40.1; rib 25.3/64.3 = 39.3 vs 39.3).
- LSM vs plain bar means (Fig. 2): pan 150 °C 44.3 vs Table 44.2; **pan 175 °C 47.7 vs Table 46.1**; pan 100 36.4 vs
  36.0. Water flux 100 / 175: 5.23 / 9.82 vs 5.2 / 9.8 ✓.

## 5. Flags and internal consistency

1. [!] **Fig. 1 vs text on frying times**: text (Exp. Ia) gives the shortest frying time as 132 s (leanest meat,
   175 °C); Fig. 1 shows the **fat brisket** patty at 175 °C reaching 72.8 °C at the centre at **~99 s** (150 °C:
   ~124 s). Either Fig. 1 is a single unrepresentative run or the 132 s is a cell mean; unresolved. The 100 °C end
   (190 s) matches the text.
2. [!] **Text says cooking temperature "had no influence on the total loss"** [V] while Table 3 gives pan-T p = 0.001
   (33.6 → 39.2 %) and Table 4 p = 0.030.
3. [!] Shrinkage direction: text says shrinkage decreased with pan T (Table 4: 11.2 → 9.7, p = 0.185, ns); Table 3
   shows 13.6 → 15.0 (p = 0.560, ns). Neither is significant.
4. [!] Density-porosity diameters: text "8.3, 10 and 12.7 cm"; Table 4 column 8.1.
5. [!] The 80 %-drip conclusion rests on the assumption that a 100 °C pan gives no evaporation, while the 2 mm layer
   reaches 84-91 °C and the contact face ~100 °C; some evaporation at 100 °C is likely. P(drip share really ≥ 70 % at
   175 °C) ~55 %.
6. [!] Patty thickness, mass, thermocouple type and logging interval, pan surface temperature recovery: all deferred
   to Part 1; one trace per condition in Fig. 1; n = 2 replicates per cell.
7. [!] Frozen start (-20 °C) and beef-fat drip (fat flux up to 6×10⁻³ g cm⁻² s⁻¹) make the thermal history specific to
   frozen meat patties.

## 6. Existing records in the repo that cite this paper

`grep -rIn -i` of `data/ docs/ src/ scripts/ tasks/` for `oroszv`, `kovácsn`/`kovacsn`, `jfoodeng.2004.10.013`,
`beefburger` (2026-11-03): **only** `data/keys/pdf_dossier_aliases.yml` l. 103 (the pending triage line). Its claims
— J Food Eng 2005, DOI, "centre / sub-surface T, water and fat loss" of double-sided-fried beefburgers — **match**;
"water and fat loss" are end-of-frying values only (not stated in the line). No ledger, registry or code record.
The VARIABLES search doc cites **Part 3** (2006), a different paper (sec. header).

## 7. Proposed ledger entries (not written)

**No entry for `trunk.yml` or `aroma_other.yml`**: the paper measures no chemical species, colour, odorant or
precursor; nothing a lane simulates. `aroma_other.yml`'s "matrix / meat ... composition data" scope would admit
Table 1 alone, but four fat / water percentages of beef cuts help no plant-based pot (P(useful) < 5 %); not proposed.
The T(t) and mass-loss data are transport observables; they are proposed under the `thermal_history` field set of
`sprague2011_extraction.md` sec. 7 (as used by `kondjoyan2016_extraction.md` sec. 7), a ledger that does not yet exist.

```yaml
  # ---- oroszvari2005 : Kovacsne Oroszvari et al. 2005 J Food Eng 71:18 (Part 2) (ingested 2026-11-03) ----
  - id: oroszvari2005_fig1_subsurface_centre_t_vs_time
    dossier: oroszvari2005_extraction.md
    paper: "Kovacsne Oroszvari, Bayod, Sjoholm & Tornberg 2005 J Food Eng 71:18"
    lab: "Lund (Tornberg)"
    product: "frozen minced-beef patty, fat brisket 39 % fat, 1.5 % NaCl, D = 10 cm (thickness ~10 mm, not printed)"
    heating_mode: double_sided_contact_frying
    setpoint: "pan 100 / 150 / 175 C, stop at centre 72 C"
    boundary: {pan_t_c: [100, 150, 175], contact_pressure: null, h_w_m2_k: null, initial_t_c: -20}
    observable: temperature_c
    positions: [centre_5mm, subsurface_2mm]
    time_min: [0, 3.2]
    representation: data_points
    location: "Fig. 1 (6 traces, one run each), digitised from vector paths in dossier sec. 4b"
    usable_as_programme: partial
    status: digitised_in_dossier
    flags: [frozen_start, no_surface_or_crust_temperature, single_trace_per_condition, thickness_not_printed, text_fig1_frying_time_mismatch]
    note: "Hottest measured point is 2 mm deep and ends at <= 100 C; a sub-crust programme at most, not a crust (Maillard) programme."
  - id: oroszvari2005_tab3_tab4_figs2_4_endpoint_mass_loss
    dossier: oroszvari2005_extraction.md
    paper: "Kovacsne Oroszvari, Bayod, Sjoholm & Tornberg 2005 J Food Eng 71:18"
    lab: "Lund (Tornberg)"
    product: "frozen minced-beef patties, 4 cuts (fat 6.7-39 %), D = 3-12.7 cm"
    heating_mode: double_sided_contact_frying
    setpoint: "pan 100 / 150 / 175 C, stop at centre 72 C (t72 90-190 s)"
    boundary: {pan_t_c: [100, 150, 175], contact_pressure: null, h_w_m2_k: null, initial_t_c: -20}
    observable: [total_loss_pct, water_loss_pct, fat_loss_pct, water_flux_g_cm2_s, fat_flux_g_cm2_s, final_t_2mm_c, porosity, shrinkage_pct]
    positions: [whole_patty, subsurface_2mm]
    time_min: [1.5, 3.2]
    representation: endpoint_lsmeans_and_cell_means
    location: "Tables 3-4 (main-effect LS means, p-values) and Figs. 2-4 (18 cells, n = 2, SE bars), dossier sec. 4c-4f"
    usable_as_programme: false
    status: digitised_in_dossier
    flags: [endpoint_only_no_time_course, lsmeans_not_cell_means, text_contradicts_table_total_loss_vs_pan_t, drip_share_rests_on_no_evaporation_assumption]
    note: "Time-resolved water/fat loss is in Part 3 (Kovacsne Oroszvari 2006, 10.1016/j.jfoodeng.2005.05.018), not on disk."
```

## 8. Supplementary material

None referenced (text layer searched for "upplementary", "Appendix": no hits). Nothing model-relevant is SI-only. The
missing pieces (patty geometry, thermocouple method, set-up) are in **Part 1** (Kovácsné, Sjöholm & Tornberg 2004,
J. Food Eng. 67:499), and mass loss vs time in **Part 3** (2006); neither on disk.
