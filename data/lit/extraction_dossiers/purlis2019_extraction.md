# Purlis & Salvadori 2009 — EXTRACTION (file stem `purlis2019`)
DOI: 10.1016/j.foodres.2009.03.007 (printed p. 865 footer). Crossref 2026-11-03: "Modelling the browning of bread
during baking", Food Res. Int. **42**(7):865-870, authors Emmanuel Purlis, Viviana O. Salvadori, issued 2009-08
(received 8 Oct 2008, accepted 17 Mar 2009). **Matches the PDF and the brief.**
[!] **Year: the file stem says 2019, the paper is 2009.** The stem is kept so the coverage gate matches
(`data/keys/pdf_dossier_aliases.yml` l. 104 says the same).

### French-bread dough cylinders (100 g, 0.04 m × 0.1 m; flour 100, water 54.1, salt / sucrose / margarine 1.6, dry yeast 1.2 parts) baked in an electric static oven at 180 / 200 / 220 °C, natural (0 m/s) and forced (0.9 m/s) convection, 40 min (220 °C: 30 min). **Measured: surface L* (computer vision, every 5 min).** Surface T(t) and a_w(t) are **simulated** by the authors' earlier model (Purlis & Salvadori 2009b), not measured. Fitted: first-order L* decay with an Arrhenius-like k(T, a_w), Eq. 13, four constants, no uncertainties.

**Read 2026-11-03 from `data/articles/purlis2019.pdf`** (6 pp., Elsevier, born-digital). Prose from
`pdftotext -layout`; **pp. 3-5 rendered at 150 dpi** (scratchpad `ingest/purlis2019_work/p-3..5.png`); **Eq. 9 at
300 dpi** (`eq9-3.png`), **Eq. 13 at 400 dpi** (`eq13-5.png`); **Fig. 3 and Fig. 4 at 300 dpi** (`fig3hi-4.png`,
`fig4hi-4.png`), read error ±0.7 L* and ±2-3 °C; Figs. 5, 6, 8a from the 150 dpi page renders, read error ±0.5 L*
(Fig. 5), ±0.02 a_w (Fig. 6), ±1 L* (Fig. 8a). Fig. 2 (photographs) and Fig. 7 (3-D surface of Eq. 13) not digitised.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** **CIDCA** (Centro de Investigación y Desarrollo en Criotecnología de Alimentos, CONICET-UNLP, La Plata)
and MODIAL, Facultad de Ingeniería, UNLP, Argentina; corresponding author **Emmanuel Purlis**, senior author
**Viviana O. Salvadori**. **Not in `data/keys/laboratories.yml`: NEW**, proposed id **`cidca_salvadori`**, lab string
"CIDCA La Plata (Purlis / Salvadori)". No existing dossier from this group. Same-group companions cited, none on disk:
Purlis & Salvadori 2007 J. Food Eng. 80:1107 (bread browning kinetics vs weight loss); 2009a/2009b J. Food Eng.
91:428 / 434 (moving-boundary baking model and its validation — the source of the simulated T and a_w used here).

## 1. What this paper answers for the model

| model need (mechanism item 8 / trunk browning) | what Purlis & Salvadori 2009 gives | tag |
|---|---|---|
| measured T(t) at positions | **none measured**; surface T(t) simulated (Fig. 4), 6 conditions | [F], simulated |
| weight / water loss vs time | **none**; surface a_w(t) simulated for 10-40 min only (Fig. 6) | [F], simulated |
| surface colour vs time | **L* measured** every 5 min, 6 conditions (Fig. 3; Fig. 5 first 10 min; Fig. 8a forced convection with fits) | [F] |
| browning kinetics: order, k, T and a_w dependence | Eq. 9 first-order in L* (to L* = 0); Eq. 13 k(T, a_w) with 4 printed constants, no SE | [M] |
| Ea | **not given** (authors: "fit parameters without physical meaning" [V]); apparent A·R ≈ 73-76 kJ/mol | [Z] |
| onset conditions | browning starts when surface T > 120 °C and a_w < 0.6 | [Q] |
| h, evaporation front, crust thickness | **none** (in 2009a/b) | [!] |
| chemistry (sugars, HMF, melanoidin, A420) | **none** | [!] |

**Answers.**
1. **Is the browning model a Maillard rate law in surface T and water content with printed constants? No, not a
   Maillard rate law** (confidence ~90 %). It is an empirical first-order decay of CIE lightness, dL*/dt = −k L*
   (Eq. 9), with k = k0 exp(−A/T), k0 = k1 + k2/a_w, A = k3 + k4/a_w (Eqs. 10-12). The constants **are printed**
   (Eq. 13): **k1 = 7.9233×10⁶, k2 = 2.7397×10⁶ min⁻¹; k3 = 8.7015×10³, k4 = 49.4738 K; k in min⁻¹, T in K** [M]
   (the × sign is lost in the typesetting, "7.923310⁶" [!]; the reading is confirmed by reproducing Fig. 7's scale and
   the Fig. 8a slopes, sec. 4e). No chemical species, no sugar or amine dependence, caramelisation and Maillard
   lumped; the L* → 0 asymptote (black) is unphysical for a real crust.
2. **Water dependence**: through surface a_w only (not moisture), both in the pre-factor and in A; k has a maximum in
   a_w 0.1-0.2 at fixed T [Q] — reproduces [Z]. The a_w values are simulated, so the a_w term is fitted against model
   output, not data (P(the a_w term would survive measured a_w) ~40 %).
3. **Usable as a heating-programme input?** As a *simulated* crust programme only: Fig. 4 gives surface T(t) rising
   25 → 101-119 °C in 5 min, then 140-192 °C by 30-40 min. It is a plausible shape for a crust pot (P ≈ 25 % that
   it is the right thing to feed; measured surface T would be preferred, e.g. Kondjoyan, Feyissa). Bread, not a
   plant-protein matrix.
4. **Usable as a browning check?** Quantitatively no (P ≈ 5 %): the engine predicts melanoidin / A420 in a well-mixed
   pot and has no L* observation operator; bread browning includes sucrose caramelisation. Qualitatively, the onset
   (T > 120 °C, a_w < 0.6) and the ordering of the six conditions are directional claims (P ≈ 40 % worth logging).

## 2. Conditions

- **Dough** [M]: wheat flour 100 %, water 54.1 %, salt 1.6 %, **sugar 1.6 %** (sucrose presumably, not specified),
  margarine 1.6 %, dry yeast 1.2 % (flour basis); 10 min mixing; 100 g cylinders ~0.1 m × 0.04 m; proving 1.5 h at
  ambient T under film (volume doubled). Flour protein / reducing sugar not printed.
- **Oven** [M]: Ariston FM87-FC electric static oven, **180, 200, 220 °C** (±3.3 °C, T-type thermocouples), natural
  convection (v = 0) or forced (v = 0.9 m/s); perforated tray, two samples per test, **two tests per condition**.
- **Duration** [F]: 40 min (180, 200 °C); 30 min (220 °C) (Figs. 3, 4).
- **pH, O2, humidity**: not reported.

## 3. Analytical methods

- **Colour** [M]: computer vision, IP camera 0.4 m above the samples, 640 × 480 jpg every 5 min (~5 s acquisition,
  oven next to camera), "common (non-ideal) light conditions" [V]; MATLAB RGB → XYZ (sRGB matrix, D65 white) → CIE
  L*a*b* (Eqs. 1-8); calibrated against Minolta CR300 plates before each run. Each point is the mean of four images
  (2 samples × 2 tests), ~10⁵ pixels each. **No ± shown on any L* point** [!]. a* and b* computed, not reported.
- **Surface T and a_w** [M]: numerical simulation with the model of Purlis & Salvadori 2009b at the same operating
  conditions; the model's h, isotherm and evaporation treatment are not printed here [!].
- **Fit** [M]: darkening-phase points only (surface T > 120 °C); Levenberg-Marquardt (lsqcurvefit) on Eq. 9 integrated
  by ode45; mean absolute relative error 3.6095 % (n = 34) (Eq. 14). No parameter SE, no residual plot beyond Fig. 8b.

## 4. Data

### 4a. Fig. 3, measured surface L* vs time (mean of 4 images, no ±) [F, 300 dpi, ±0.7]

nat = natural convection (filled symbols), forc = forced (open symbols).

| t (min) | 180 nat | 180 forc | 200 nat | 200 forc | 220 nat | 220 forc |
|---|---|---|---|---|---|---|
| 0 | 88.7 | 87.6 | 87.6 | 92.8 | 87.3 | 89.7 |
| 5 | 91.2 | 88.8 | 90.6 | ~89 (hidden) | 85.9 | 86.3 |
| 10 | 90.2 | 85.4 | 85.4 | 79.8 | 70.2 | 70.8 |
| 15 | 85.2 | 80.5 | 75.7 | 71.8 | 59.1 | 59.1 |
| 20 | 79.2 | 75.8 | 69.8 | 65.3 | 50.0 | 51.8 |
| 25 | 76.7 | 71.3 | 64.5 | 59.7 | 43.2 | 43.5 |
| 30 | 73.2 | 68.9 | 60.6 | 55.8 | 37.4 | 37.7 |
| 35 | 70.4 | 65.8 | 56.7 | 52.5 | — | — |
| 40 | 67.7 | 63.6 | 53.0 | 48.7 | — | — |

t = 0 and 5 min values from the finer-scale Fig. 5 [F, ±0.5] (the "enlightenment" lag: L* rises 0-5 min at 180 nat,
200 nat, 180 forc). The 200 forc 5 min point is hidden behind other symbols in both figures [!]. Fig. 8a (forced
convection only) agrees with this table within 0.5 L*.

### 4b. Fig. 4, simulated surface temperature (°C) [F, 300 dpi, ±2-3; simulated, not measured]

| t (min) | 180 nat | 180 forc | 200 nat | 200 forc | 220 nat | 220 forc |
|---|---|---|---|---|---|---|
| 0 | 25 | 25 | 25 | 25 | 25 | 25 |
| 5 | 102 | 105 | 107 | 110 | 112 | 119 |
| 10 | 111 | 116 | 120 | 134 | 140 | 150 |
| 15 | 122 | 130 | 136 | 149 | 158 | 168 |
| 20 | 131 | 139 | 147 | 158 | 170 | 178 |
| 25 | 137 | 145 | 155 | 166 | 179 | 186 |
| 30 | 142 | 150 | 162 | 171 | 186 | 192 |
| 35 | 147 | 154 | 167 | 175 | — | — |
| 40 | 151 | 157 | 171 | 178 | — | — |

### 4c. Fig. 6, simulated surface water activity, darkening phase [F, 150 dpi, ±0.02; simulated]

| t (min) | 180 nat | 180 forc | 200 nat | 200 forc | 220 nat | 220 forc |
|---|---|---|---|---|---|---|
| 10 | 0.58 | 0.55 | 0.69 | 0.31 | 0.40 | 0.22 |
| 15 | 0.51 | 0.40 | 0.56 | 0.16 | 0.19 | 0.11 |
| 20 | 0.39 | 0.28 | 0.40 | 0.10 | 0.11 | 0.065 |
| 25 | 0.29 | 0.20 | 0.29 | ~0.07 | 0.075 | ~0.045 |
| 30 | 0.22 | 0.15 | 0.21 | ~0.05 | 0.055 | ~0.04 |
| 35 | 0.13 | 0.10 | 0.12 | 0.03 | — | — |
| 40 | 0.14 | 0.10 | 0.12 | 0.03 | — | — |

### 4d. Eq. 13 (p. 869), the fitted browning rate constant [M]

k (min⁻¹) = (7.9233×10⁶ + 2.7397×10⁶ / a_w) · exp[ −(8.7015×10³ + 49.4738 / a_w) / T ], T in K, a_w surface water
activity; with dL*/dt = −k L* (Eq. 9). Fit on n = 34 darkening-phase points (T_surface > 120 °C), all six conditions
pooled; e_abs = 3.6095 %. No standard errors, no correlation matrix.

### 4e. Derived from Eq. 13 [Z]

| T (°C) | a_w | k (min⁻¹) | A·R (kJ/mol) |
|---|---|---|---|
| 120 | 0.6 | 0.0025 | 73.0 |
| 160 | 0.2 | 0.023 | 74.4 |
| 180 | 0.15 | 0.058 | 75.1 |
| 200 | 0.1 | 0.128 | 76.5 |

Checks: (i) Fig. 7 tops out at ~0.12 min⁻¹ at 200 °C, a_w ≈ 0.1 ✓ and is ~0 at 120 °C, a_w 0.6 ✓ — confirms the ×10ⁿ
reading. (ii) Fig. 8a slopes: 220 forc, ln(70.9/37.7)/20 min = 0.032 min⁻¹ vs Eq. 13 at the run's mean simulated
surface (~175 °C, a_w ~0.08) ≈ 0.039; 180 forc, ln(80.5/63.6)/25 = 0.0094 vs Eq. 13 at ~145 °C, a_w ~0.2 ≈ 0.011 —
consistent to ~20 %. (iii) A·R, an Arrhenius-equivalent, is 73-77 kJ/mol over the window; k1 and k3 are strongly
correlated in any such fit, so A·R is the only combination worth quoting, and the authors disclaim it.

## 5. Flags and internal consistency

1. [!] **Only L* is measured.** T and a_w are model outputs of a separate paper; the kinetic fit's inputs therefore
   carry that model's (unreported) errors. "Thermal history" here means *simulated* history.
2. [!] **First-order decay to L* = 0** has no plateau (L*∞); over 40 min it fits, beyond that it will over-darken.
3. [!] Eq. 13 typesetting drops the multiplication signs ("7.923310⁶"); the reading 7.9233×10⁶ etc. is confirmed by
   sec. 4e. P(misreading) < 5 %.
4. [!] **No error bars** on L*; replicates 2 tests × 2 samples; "non-ideal" lighting; the 0-5 min "enlightenment"
   (L* rising) is attributed to surface smoothing, a measurement artefact [Q].
5. [!] Fig. 6 non-monotonic in oven T: at 10 min, simulated a_w is higher at 200 °C natural (0.69) than at 180 °C
   natural (0.58); 200 forced drops far faster than 200 natural. Plausible only if the simulated crust front behaves
   very differently with v; unexplained.
6. [!] n = 34 fit points vs my count of L* points with simulated surface T > 120 °C (Fig. 4): 6 + 6 + 6-7 + 7 + 5 + 5
   = 35-36 (180 nat, 180 forc, 200 nat, 200 forc, 220 nat, 220 forc) [Z] — consistent within the borderline points.
7. [!] Sucrose 1.6 % plus flour sugars; the authors note caramelisation contributes at high T and low a_w [Q]; Eq. 13
   cannot separate it.

## 6. Existing records in the repo that cite this paper

`grep -rIn -i` of `data/ docs/ src/ scripts/ tasks/` for `purlis`, `salvadori`, `foodres.2009.03.007`
(2026-11-03): two records.
- `data/keys/pdf_dossier_aliases.yml` l. 104: "Purlis & Salvadori 2009 (NOT 2019 ...) Food Res Int 42:865 ...
  modelling bread browning during baking" — **matches** (authors, year, journal, volume, page, DOI).
- `docs/guides/LITERATURE_SEARCH_2026-10_VARIABLES.md` l. 66 (V1 table): title, venue, DOI **match**; class Q,
  verdict MAYBE. Field by field: "non-isothermal browning kinetics driven by local T and a_w" — **matches**, but the
  local T and a_w are **simulated**, which the row does not say; "crust 140-160 °C" — **not in the paper as such**:
  simulated surface T spans 102-192 °C over the runs (Fig. 4), 140-160 °C is only the 180-200 °C endpoints; the
  paper's thresholds are T > 120 °C and a_w < 0.6. The quoted context ("accounting for local temperature and water
  activity variations") is from a citing paper; the abstract here has a close variant ("the variation of local
  temperature and water activity is taken into account" [V]) — **matches in substance**.
No ledger, registry or code record.

## 7. Proposed ledger entries (not written)

`trunk.yml` carries time-resolved browning; the L* series and Eq. 13 go there, `simulable: false` (no L* operator;
bread crust with simulated T / a_w). The simulated T(t) / a_w(t) are a `thermal_history` candidate (field set of
`sprague2011_extraction.md` sec. 7, as in `kondjoyan2016_extraction.md`); not `aroma_other.yml` (no odorant,
precursor or matrix-composition data).

```yaml
  # ---- purlis2019 : Purlis & Salvadori 2009 Food Res Int 42:865 (stem says 2019) (ingested 2026-11-03) ----
  - id: purlis2019_fig3_surface_lightness_vs_time
    dossier: purlis2019_extraction.md
    paper: "Purlis & Salvadori 2009 Food Res Int 42:865 (repo key purlis2019)"
    lab: "CIDCA La Plata (Purlis / Salvadori)"
    system: "French bread dough (flour 100, water 54.1, salt/sugar/margarine 1.6, yeast 1.2 parts), 100 g cylinders"
    medium: "bread crust, electric static oven, natural or forced (0.9 m/s) convection; surface T and a_w simulated, not measured"
    temperature_c: [180, 200, 220]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [0, 40]
    species: ["CIELAB L* (surface, computer vision)"]
    kind: time_series
    location: "Fig. 3 (6 conditions x 5-min points), Fig. 5 (0-10 min), Fig. 8a (forced, with fits); dossier sec. 4a"
    status: digitised_in_dossier
    quant_class: other
    simulable: false
    condition_quote: "three constant oven temperatures: 180, 200 and 220 °C"
    flags: [matrix_food, oven_air_temperature_not_product, product_temperature_simulated, no_error_bars, colour_not_chemistry, caramelisation_confounded]
    note: "Lag 0-10 min (L* rises: surface smoothing); darkening once simulated surface T > 120 C and a_w < 0.6. Directional checks only; no L*-to-melanoidin operator."
  - id: purlis2019_eq13_lightness_rate_constant_t_aw
    dossier: purlis2019_extraction.md
    paper: "Purlis & Salvadori 2009 Food Res Int 42:865 (repo key purlis2019)"
    lab: "CIDCA La Plata (Purlis / Salvadori)"
    system: "French bread crust (as above)"
    medium: "bread crust; k(T, a_w) fitted against simulated surface T and a_w"
    temperature_c: [120, 192]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [10, 40]
    species: ["first-order L* decay constant k (min^-1)"]
    kind: rate_constant
    location: "Eq. 13 p. 869 (k1..k4 printed, no SE), Fig. 7 surface; dossier sec. 4d-4e"
    status: digitised_in_dossier
    quant_class: other
    simulable: false
    condition_quote: "only data from the darkening phase was considered"
    flags: [matrix_food, empirical_colour_kinetics, inputs_simulated, no_parameter_uncertainty, l_star_decays_to_zero, multiplication_signs_lost_in_print]
    note: "Not an Ea (authors); A*R = 73-77 kJ/mol over 120-200 C [Z]. k max at a_w 0.1-0.2. n = 34, e_abs 3.6 %."
```

Thermal-history candidate:

```yaml
  - id: purlis2019_fig4_fig6_simulated_surface_t_aw
    dossier: purlis2019_extraction.md
    paper: "Purlis & Salvadori 2009 Food Res Int 42:865 (repo key purlis2019)"
    lab: "CIDCA La Plata (Purlis / Salvadori)"
    product: "French bread, 100 g cylinder 0.04 m x 0.1 m"
    heating_mode: oven_natural_or_forced_convection
    setpoint: "oven air 180 / 200 / 220 C (+-3.3 C); v = 0 or 0.9 m/s"
    boundary: {air_t_c: [180, 200, 220], air_velocity_m_s: [0, 0.9], humidity: null, h_w_m2_k: null}
    observable: [surface_temperature_c, surface_water_activity]
    positions: [surface]
    time_min: [0, 40]
    representation: simulation
    location: "Fig. 4 (T, 0-40 min) and Fig. 6 (a_w, 10-40 min), model of Purlis & Salvadori 2009b; dossier sec. 4b-4c"
    usable_as_programme: partial
    status: digitised_in_dossier
    flags: [simulated_not_measured, model_parameters_in_other_paper, bread_not_plant_protein, aw_non_monotonic_in_oven_t]
    note: "A crust-shaped T(t) with a_w(t): 25 -> ~105-119 C in 5 min, 140-192 C by 30-40 min."
```

## 8. Supplementary material

None referenced (text layer searched for "upplementary", "Appendix": no hits). The model-relevant missing items are in
the companion papers, not SI: the baking model's h, isotherm and evaporation scheme (Purlis & Salvadori 2009a/b,
J. Food Eng. 91:428-442) and the earlier weight-loss browning model (2007, J. Food Eng. 80:1107); none on disk.
