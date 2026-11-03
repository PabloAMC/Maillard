# Sprague & Colvin 2011 — EXTRACTION (file stem `sprague2011`)
Food Res. Int. **44**:789-797, DOI 10.1016/j.foodres.2011.01.011 (printed p. 1). Received 20 Oct 2010, accepted 8 Jan 2011.
Crossref 2026-11-03: "A mixture-enthalpy fixed-grid model for temperature evolution and heterocyclic-amine formation in
a frying beef patty", Food Research International 44, 789-797, authors Sprague, Colvin, issued 2011-04. **Matches the
PDF and the brief.** No naming caveat.

### Modelling paper: 2-D axisymmetric heat + liquid-water transport (mixture enthalpy, vaporisation at 102 °C, shrinkage) for a pan-fried beef patty, Legendre spectral finite elements; two free parameters (h_pan, τv) tuned on Tran et al. 2002, tested on Ikediala et al. 1996, then an Arrhenius HAA law fitted to Salmon et al. 2000. **No new measurement: every measured curve is re-plotted from those three papers (one, the Tran pan-T record, otherwise unpublished). No food-surface or crust temperature is measured; A is not reported.**

**Read 2026-11-03 from `data/articles/sprague2011.pdf`** (9 pp., Elsevier, born-digital). Prose and equations from
`pdftotext -layout`, equations and Table 1 (p. 793) checked against a 150 dpi render (`ingest/sprague2011_work/p-5.png`).
**Fig. 4a-c (p. 794) read from 300 dpi crops** (`fig4a-6.png`, `fig4b-6.png`, `fig4c-6.png`, read error ±2 °C,
±25 s); Fig. 3 (p. 793), Figs. 5-9 (pp. 794-795) from 150 dpi renders (±3 °C); **Fig. 10 (p. 796) from a 300 dpi crop**
(`fig10-8.png`, log axis, ±0.05 decade ≈ ±12 %). Equation checks by `awk` [Z].

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited / re-plotted from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

**Laboratory.** Center for Computational Biology, **University of California, Merced**; corresponding author
**Michael A. Sprague** (now NREL, Golden), senior **Michael E. Colvin**. **Not in `data/keys/laboratories.yml`: NEW**,
proposed id **`ucmerced_colvin`**, string "UC Merced (Sprague / Colvin)". The data-generating groups are different:
Tran 2002 and Salmon 2000 are Lawrence Livermore (Knize / Felton; Colvin co-authored Tran 2002) → proposed
**`llnl_felton_knize`** for any entry keyed on those data; Ikediala 1996 (Correia, Nova Scotia) → no id proposed.
No same-lab dossier in `data/lit/extraction_dossiers/` (grep merced / livermore / felton / knize: no hit).

## 1. What this paper answers for the model

| model need (mechanism item 8, thermal history) | what Sprague 2011 gives | tag |
|---|---|---|
| pan / contact-surface T(t) | Fig. 3 + Eq. 20 (Tran, T⁰pan 146.9 °C, 0-2500 s); Fig. 9 + Eq. 21 (Ikediala, 180 °C, 0-480 s): **measured pan T with a printed analytical fit** | [C]/[M] |
| food surface / crust T(t) | **not measured**; model field at 2500 s only (Fig. 5a colour map, 31.5-140 °C) | [P] |
| interior T(t) at positions | Fig. 4: 3 heights (6, 9, 12 mm from the pan, r ≈ 0), 2-3 thermocouples each, 0-2500 s; Fig. 8: 3.5, 5.5, 9.5 mm, 0-480 s | [C]/[F] |
| air T, h | T_air 20 °C, h_air 40 W/(m² °C) (chosen from literature 7-60); h_pan 150 W/(m² °C) **fitted** (literature 208-8268) | [M] |
| water transport / latent heat | L = 2.260 × 10⁶ J/kg; T_boil 102 °C; vapour rate (Xw − Xw,eq)/τv, τv 0.1 s fitted; De (Shilton 2002) Dd 1.43 × 10⁻⁶ m²/s, DT 1.58 × 10³ K, DX 6.72 × 10⁻² kg/kg, zero below Td 30 °C; Xw,eq = Xw⁰ exp(−0.0132 (T − 30)) | [M]/[C] |
| water loss vs time | **no measured weight loss**; model water-fraction map at 2500 s only (Fig. 5b, 0.147-0.622) | [P] |
| HAA amounts | Fig. 10: total [HA] (ng/g) vs T⁰pan 160-250 °C, single- vs multi-flip, 4-5 patties each (Salmon 2000) | [C]/[F] |
| HAA kinetics | dHA/dt = A exp(−Ea/RgT); **Ea 25.3 kcal/mol = 105.9 kJ/mol** fitted to the single-flip slope (Tran's FE model: 18 kcal/mol = 75.3 kJ/mol); **A never reported** (model normalised at 160 °C) | [M]/[Z] |
| thermophysical properties | Choi & Okos polynomials (cp, k, ρ of fat, protein, water) and vapour properties (Eqs. 16-17) | [C] |

**Answers.**
1. **Measured T-time profiles: yes, but interior and pan only** (confidence ~95 %). Pan T is the only "surface"
   record; interior thermocouples at 3.5-12 mm never exceed ~98 °C (sec. 4b). The Maillard-relevant layer (pan face,
   modelled at ~140 °C at 2500 s with a pan at ~162 °C, Fig. 5a [P]) is **not measured**.
2. **Water loss / crust thickness vs time: neither measured.** Fig. 5b shows modelled liquid-water volume fraction
   down to 0.147 at the pan face at 2500 s (initial 0.623). No crust.
3. **HAA amounts vs time: no time course**, endpoints only (cooked to centre 70.6 °C, 360-532 s [Z from Eqs. 22-23]);
   units ng/g; the paper's nomenclature says per g **uncooked** meat [M], Salmon's basis and which HAAs are summed are
   not stated here [!].
4. **Fitted constants with units:** h_pan = 150 W/(m² °C) (±20 % moves T by ≤ 3.5 °C, Fig. 6 [M]); τv = 0.1 s (results
   insensitive over 0.05-0.2 s, final T 95.7-98.7 °C at z = 6 mm [M]); Ea = 25.3 × 10³ cal/mol; A in
   ng-HA g⁻¹ s⁻¹ but **not given**. The Ea is an apparent, transport-convolved value; Eq. 15 is called
   "first-order" [V] but has no concentration term (zero order, no precursor depletion) [!].
5. **Usable as a ready heating programme for a Maillard pot? Partly** (confidence ~75 % for "partly, not directly").
   Eq. 20 and Eq. 21 are clean, printed, measured **pan-contact** programmes (a boundary condition), usable directly
   only if the pot is the pan-face layer and one accepts T_layer ≈ T_pan (biased high; the model's own surface is
   ~20 °C below the pan at 2500 s [P]). The interior curves are a heat-penetration programme for the core (≤ 100 °C:
   slow Maillard, relevant to Strecker / lipid lanes in the bulk). A crust T(t) would have to come from running this
   (or a simpler 1-D) transfer model with h_pan 150 and the Eq. 20 boundary [P]. **The most transferable lesson**:
   the same Salmon data give Ea 75 or 106 kJ/mol depending on the thermal model behind them; an apparent Ea from
   cooked food is a property of the thermal model as much as of the chemistry (bears on B64).

## 2. Conditions (of the three re-used experiments; Table 1, p. 793) [M]

| | Tran 2002 | Ikediala 1996 | Salmon 2000 |
|---|---|---|---|
| patty R⁰ × Z⁰ (m) | 0.0475 × 0.0250 | 0.0450 × 0.0150 | 0.0450 × 0.0150 |
| fat / protein / water (vol. fraction) | 0.238 / 0.139 / 0.623 | 0.121 / 0.169 / 0.710 | 0.238 / 0.139 / 0.623 |
| T⁰ (°C) | 13.5 | 7.0 | 10.6 |
| T⁰pan (°C) | 146.9 (Eq. 20) | 180.0 (Eq. 21) | 160-250 ("Variable"; Eq. 20 shape offset) |
| ρ⁰ / k⁰ / c⁰p | 1024 / 0.4390 / 3330 | 1044 / 0.4672 / 3513 | 1024 / 0.4358 / 3329 |
| frying | single-sided, 2500 s | single-sided, 500 s | flipped once at 300 s or every 60 s, to centre 70.6 °C |
| h_air / h_pan / τv | 40 / 150 / 0.1 (all three) | | |

Thawed 20 % fat-by-weight beef patties (Tran, Salmon) [M]. T_air 20 °C. Heating mode: electric pan / griddle
(not described). No humidity, no pH, no composition beyond the fractions. Numerics: NR × NZ = 14 × 9 or 14 × 5, N = 6,
Δt 0.046-0.051 s, explicit RK3.

## 3. Analytical / modelling methods

- **Model** (short): mixture enthalpy over fat, protein, liquid water and vapour (Eqs. 1-5), conductive only (no
  advection), convective BCs h_air (air faces) and h_pan (pan face, Eq. 8); liquid water by diffusion with De(T, Xw)
  (Eq. 10) and boundary Xw = Xw,eq(T) (Eqs. 11-12); vapour generated at T ≥ 102 °C at rate (Xw − Xw,eq)/τv and assumed
  to leave immediately; protein and fat immobile; uniform shrinkage α = (V(t+Δt)/V(t))^⅓ (Eq. 19). HAAs by integrating
  Eq. 15 at every node, volume-averaged.
- **Measurements re-used**: thermocouples (Tran: 2-3 per height; Ikediala: one series per height), pan surface T
  records; Salmon's HAA analysis is not described [!]. No ± or replicate statistics are printed anywhere.

## 4. Data

### 4a. Pan-temperature programmes (printed fits; checked) [M], t̂ = t/(1 s), T in °C

- **Eq. 20 (Tran, T⁰pan 146.9)**: T = 135.2 − 0.01036 t̂ + 11.74 exp(−0.06445 t̂) for t̂ < 216;
  T = 133.0 + 0.02166 (t̂ − 216) − 3.93 × 10⁻⁶ (t̂ − 216)² for t̂ ≥ 216.
  Checks [Z]: T(0) = 146.94; branches meet (132.96 vs 133.0); minimum ~133 °C at 216 s; T(2500) = 161.97 (Fig. 3 data
  ~161-162 [F]); crosses 146.9 at ~957 s (text: "exceeds ... for t > 975 s" [V]).
- **Eq. 21 (Ikediala, 180)**: T = 180 − 0.11117 t̂ + 0.00044487 t̂² − 4.7894 × 10⁻⁷ t̂³.
  Checks [Z]: T(0) = 180 (first data point ~178.7 [F]); minimum 171.6 °C at 174 s; T(480) = 176.2 (data ~176.2 [F]).
  Valid only 0-480 s (cubic diverges beyond).

### 4b. Fig. 4 (Tran 2002, r ≈ 0, T⁰pan 146.9): measured interior T, °C, ±2 °C [F]

| height above pan | series | 250 s | 500 s | 1000 s | 1500 s | 2000 s | 2500 s |
|---|---|---|---|---|---|---|---|
| 6 mm | circles | ~54 | 79 | 93 | 97 | 97 | 98 |
| 6 mm | triangles | ~44 | 65 | 83 | 91 | 93.5 | 96.5 |
| 9 mm | circles | — | 57 | 75 | 85 | 90.5 | 96 |
| 9 mm | triangles | — | 53.5 | 75 | 85.5 | 90.5 | 95 |
| 9 mm | squares | — | 49 | 70 | 78.5 | 83 | 87 |
| 12 mm | triangles | — | 43 | 63 | 75 | 80 | 84 |
| 12 mm | circles | — | 38 | 59 | 70 | 76.5 | 81 |
| 12 mm | squares | — | 35.5 | 57 | 67 | 73 | 77 |

All start at 13.5 °C. The 6 mm pair spreads 10-14 °C at 500 s (probe position uncertainty); a plateau at 96-98 °C
from ~1500 s (6 mm) marks the evaporation front, not 100-102 °C [!]. Model (LSFE, τv 0.1) tracks the lower 6 mm
series and the lower 9 / 12 mm series; τv → ∞ (no vaporisation) overshoots to ~115 °C at 6 mm.

### 4c. Fig. 8 (Ikediala 1996, T⁰pan 180, 15 mm patty): measured interior T, °C, ±3 °C [F]

| height | 100 s | 200 s | 300 s | 400 s | 480 s |
|---|---|---|---|---|---|
| 3.5 mm | 56 | 73 | 83 | 90 | 94 |
| 5.5 mm | 31 | 54 | 67 | 77.5 | 82.5 |
| 9.5 mm | 13 | 25 | 37 | 49.5 | 60 |

(Initial 7 °C.) Independent test of h_pan / τv: model agrees within ~3 °C except 3.5 mm at 50-120 s (data up to ~8 °C
above model) [F].

### 4d. Fig. 10 (Salmon 2000): total [HA] (ng/g) vs T⁰pan, individual patties, ±12 % [F]

| T⁰pan (°C) | single flip (exp.) | multi flip (exp.) |
|---|---|---|
| 160 | 2.1, 1.36, 1.27, 0.88, 0.79 | 0.24, 0.21, 0.12, 0.090, 0.082 |
| 180 | 6.2, 3.4, 2.5, 2.3 | 0.55, 0.15, 0.13, 0.11 |
| 200 | 8.8, 7.6, 6.2, 3.2 | 0.75, 0.71, 0.52, 0.48, 0.092 |
| 250 | 21.1, 19.8, 15.6, 13.0 | 5.2, 4.7, 4.2, 3.6 |

Overlapping markers may hide points (n not printed). LSFE single-flip curve passes ~1.3 (160) → ~18 (250); LSFE
multi-flip ~0.42 (160) → ~6.8 (250), i.e. still 2-5× above the multi-flip data at 160-200 °C (authors: drippings loss
not modelled [Q]). Model curves are **normalised to the single-flip mean at 160 °C**, so only the slope carries
information.

### 4e. Cooking times (Eqs. 22-23, from Salmon 2000) [M], checked [Z]

t_cook,sing = 3600 / (18.21 − 0.1289 T + 3.630 × 10⁻⁴ T²) s → 523 (160), 532 (180), 518 (200), 415 (250) s;
t_cook,mult = 3600 / (18.59 − 0.1144 T + 3.130 × 10⁻⁴ T²) s → 434, 442, 437, 377 s; 360 s at 260 °C.

## 5. Flags and internal consistency

1. [!] **No original data**; measured curves are re-plots, the Tran pan-T record "(unpublished)" [V]. Cite the
   primary papers for any datum; this paper is the carrier of Eq. 20 only.
2. [!] h_pan and τv were fitted on the Tran data that Fig. 4 then "validates"; only Fig. 8 (Ikediala) is out-of-sample,
   and there the authors themselves used h_pan = 250 in their own model.
3. [!] HAA law: zero order in precursors (called first order); A unreported; normalisation at 160 °C leaves only Ea;
   HAAs summed over unnamed species; basis "uncooked meat" vs Salmon's (not stated here).
4. [!] Table 1 gives t_max 360-523 s for Salmon; Eq. 22 gives up to 532 s at 170-190 °C (×1.02; minor).
5. [!] Nomenclature gives De units as "kg-water kg-solid⁻¹" [V]; De is a diffusivity (m²/s, as Dd). Typo.
6. Mixture k⁰ (0.4390) reproduces from Eq. 16 and the Tran fractions at 13.5 °C [Z: 0.4390]: **match**.
7. The 6 mm plateau at 96-98 °C sits below the model's T_boil = 102 °C; the model reaches the same plateau through
   τv, not through T_boil (both fine within ±2 °C read error).

## 6. Existing records in the repo that cite this paper

`grep -rIn -i` for `sprague`, `colvin`, `foodres.2011.01.011` over `data/ docs/ src/ scripts/ tasks/` (2026-11-03).
No ledger, registry, prior or code record cites it. Hits:

| record | claim | paper | verdict |
|---|---|---|---|
| `data/keys/pdf_dossier_aliases.yml` l. 102 (pending) | Food Res Int 44:789, DOI as above; mixture-enthalpy model of T and HAA formation, frying beef patty | same | **matches** |
| `docs/guides/LITERATURE_SEARCH_2026-10_VARIABLES.md` l. 60 | class M; "template for coupling phase change, heat transfer and reaction kinetics in a frying patty"; MAYBE | yes (Eqs. 5, 13, 15) | **matches**; add: no original data, A not given |
| same, l. 43 | "coupling to a reaction exists for HAAs in beef (Kondjoyan, Sprague)"; transport models give "surface temperature and moisture vs time with measured validation" | coupling: **yes** (Sprague). Validation: **interior T only**, no surface T, no moisture data | **partly** for Sprague |

## 7. Proposed ledger entries (not written)

**HAA data set → `aroma_other.yml`** (nearest scope; simulable: false). **Thermal programmes: no existing ledger fits**
(every ledger keys a chemical species and pot conditions; none has positions, boundary conditions or a T(t)
representation). Proposed minimal field set for a new `thermal_history` lane (or a section of
`docs/validation/thermal_history_facts.yml`, which today records only the pots' own come-up): `id, dossier, paper, lab,
product, heating_mode, setpoint, boundary {air_t_c, air_velocity_m_s, humidity, h_w_m2_k}, observable, positions,
time_min, representation (data_points | printed_fit | model_output), location, usable_as_programme, status, flags`.
`representation: printed_fit` is the only kind a pot can take as-is.

```yaml
  # ---- sprague2011 : Sprague & Colvin 2011 (frying-patty transport model; Salmon 2000 HAAs re-plotted) (ingested 2026-11-03) ----
  - id: sprague2011_fig10_salmon2000_total_haa_vs_pan_temperature
    dossier: sprague2011_extraction.md
    paper: "Sprague & Colvin 2011 Food Res Int 44:789 (data re-plotted from Salmon et al. 2000 J Natl Cancer Inst 92:1773)"
    lab: "LLNL (Felton / Knize)"
    system: "thawed 20 % fat ground-beef patty, R 45 mm x 15 mm, pan-fried to centre 70.6 C, flipped once (300 s) or every 60 s"
    medium: "whole patty, open pan, air 20 C"
    temperature_c: [160, 250]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [6.0, 8.9]
    species: [total_HAA]
    kind: endpoint
    location: "Fig. 10 (individual patties, log axis, digitised +-12 %, dossier sec. 4d); Ea 25.3 kcal/mol in sec. 4.3 text"
    status: digitised_in_dossier
    quant_class: unknown
    simulable: false
    condition_quote: "pan-fried with various initial pan temperatures until the center temperature reached 70.6 °C"
    flags: [secondary_replot, temperature_is_initial_pan_T, haa_species_not_named, basis_uncooked_meat_per_nomenclature, model_normalised_at_160C, A_not_reported, zero_order_rate_law, "also_relevant: thermal_history"]
    note: "Apparent Ea of total HAAs: 105.9 kJ/mol with this transport model vs 75.3 kJ/mol with Tran 2002's (same data). An Ea from cooked food depends on the thermal model behind it. No HAA chemistry in the engine."
```

Thermal-history candidates (proposed lane):

```yaml
  - id: sprague2011_eq20_tran2002_pan_temperature_programme
    dossier: sprague2011_extraction.md
    paper: "Sprague & Colvin 2011 Food Res Int 44:789 (pan record of Tran et al. 2002, unpublished there)"
    lab: "LLNL (Felton / Knize)"
    product: "beef patty R 47.5 mm x 25 mm, 20 % fat, T0 13.5 C"
    heating_mode: pan_contact_single_sided
    setpoint: "pan preheated to 146.9 C"
    boundary: {air_t_c: 20, air_velocity_m_s: null, humidity: null, h_w_m2_k: {pan: 150, air: 40, pan_fitted: true}}
    observable: pan_surface_temperature
    positions: [pan_surface]
    time_min: [0, 41.7]
    representation: printed_fit
    location: "Eq. 20 and Fig. 3 (dossier sec. 4a)"
    usable_as_programme: boundary_only
    status: digitised_in_dossier
    flags: [pan_not_food_surface, secondary_replot]
  - id: sprague2011_eq21_ikediala1996_pan_temperature_programme
    dossier: sprague2011_extraction.md
    paper: "Sprague & Colvin 2011 Food Res Int 44:789 (pan record of Ikediala et al. 1996 J Food Sci 61:796)"
    lab: "Ikediala / Correia (no id)"
    product: "beef patty R 45 mm x 15 mm, X_fat 0.121, T0 7.0 C"
    heating_mode: pan_contact_single_sided
    setpoint: "pan 180 C"
    boundary: {air_t_c: 20, air_velocity_m_s: null, humidity: null, h_w_m2_k: {pan: 150, air: 40, pan_fitted: true}}
    observable: pan_surface_temperature
    positions: [pan_surface]
    time_min: [0, 8.0]
    representation: printed_fit
    location: "Eq. 21 and Fig. 9 (dossier sec. 4a); cubic, valid 0-480 s only"
    usable_as_programme: boundary_only
    status: digitised_in_dossier
    flags: [pan_not_food_surface, secondary_replot, fit_diverges_beyond_480s]
  - id: sprague2011_fig4_tran2002_interior_temperatures
    dossier: sprague2011_extraction.md
    paper: "Sprague & Colvin 2011 Food Res Int 44:789 (thermocouples of Tran et al. 2002 Food Chem Toxicol 40:673)"
    lab: "LLNL (Felton / Knize)"
    product: "as sprague2011_eq20_tran2002_pan_temperature_programme"
    heating_mode: pan_contact_single_sided
    setpoint: "pan Eq. 20"
    boundary: {air_t_c: 20, air_velocity_m_s: null, humidity: null, h_w_m2_k: {pan: 150, air: 40, pan_fitted: true}}
    observable: product_temperature
    positions: ["z 6 mm (2 probes)", "z 9 mm (3)", "z 12 mm (3)", "r ~ 0"]
    time_min: [0, 41.7]
    representation: data_points
    location: "Fig. 4a-c (300 dpi, +-2 C, dossier sec. 4b)"
    usable_as_programme: core_only
    status: digitised_in_dossier
    flags: [interior_below_100C, probe_spread_10C, used_to_fit_h_pan]
  - id: sprague2011_fig8_ikediala1996_interior_temperatures
    dossier: sprague2011_extraction.md
    paper: "Sprague & Colvin 2011 Food Res Int 44:789 (thermocouples of Ikediala et al. 1996)"
    lab: "Ikediala / Correia (no id)"
    product: "as sprague2011_eq21_ikediala1996_pan_temperature_programme"
    heating_mode: pan_contact_single_sided
    setpoint: "pan Eq. 21"
    boundary: {air_t_c: 20, air_velocity_m_s: null, humidity: null, h_w_m2_k: {pan: 150, air: 40, pan_fitted: true}}
    observable: product_temperature
    positions: ["z 3.5 mm", "z 5.5 mm", "z 9.5 mm", "r ~ 0"]
    time_min: [0, 8.0]
    representation: data_points
    location: "Fig. 8 (150 dpi, +-3 C, dossier sec. 4c)"
    usable_as_programme: core_only
    status: digitised_in_dossier
    flags: [out_of_sample_test_of_h_pan, secondary_replot]
```

Transport constants (h_air 40, h_pan 150 W/(m² °C); τv 0.1 s; L 2.260 × 10⁶ J/kg; T_boil 102 °C; De constants;
δwe 0.0132 °C⁻¹) are recorded in sec. 1; no ledger entry until a transfer model exists to read them.

## 8. Supplementary material

None referenced (text layer searched for "upplementa", "Appendix": no hit). The model-relevant missing items are in
the primary papers, not SI: Tran 2002 (Food Chem Toxicol 40:673; its pan record unpublished), Ikediala 1996
(J Food Sci 61:796), Salmon 2000 (J Natl Cancer Inst 92:1773, HAA species and basis). None is on disk.
