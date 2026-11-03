# Bassama, Brat, Bohuon, Hocine, Boulanger & Günata 2011 — EXTRACTION (file stem `bassama2011`)
DOI: 10.1016/j.foodres.2011.03.018 (printed p. 1452). Crossref 2026-11-03: "Acrylamide kinetic in plantain during
heating process: Precursors and effect of water activity", Food Res. Int. **44**:1452-1458, issued 2011-06, authors
J. Bassama, P. Brat, P. Bohuon, **B. Hocine**, R. Boulanger, Z. Günata. **Matches the PDF.**
[!] The brief's author list omits **B. Hocine** (fourth author, CIRAD); the paper and Crossref both carry six authors.

### Freeze-dried, milled green-purchased plantain (Musa AAB "barraganete") rehydrated to initial a_w 0.972 / 0.904 / 0.430 (25 °C; 70 / 41 / 8 % water wb), ~10 g paste in a sealed steel test cell, two-bath oil heating (250 °C boost bath, then 140 / 160 / 180 / 200 °C), up to 90 min; acrylamide by isotope-dilution UPLC-MS/MS. Single-response fit of Asn → AA → DP (first order, Arrhenius, T_ref 170 °C) on the measured centre T(t): Table 2 gives k_F, k_E, Ea_F, Ea_E per a_w with 95 % CI. **Asparagine and sugars measured at t = 0 only; no Asn / sugar time course, no pH, no T(t) curves printed, a_w at 25 °C before heating only.**

**Read 2026-11-03 from `data/articles/bassama2011.pdf`** (7 pp., Elsevier, born-digital). Prose from
`pdftotext -layout`; pp. 1455-1458 rendered at 150 dpi (scratchpad `ingest/bassama2011_work/p-4..7.png`) for Table 1,
Table 2, Figs. 1, 2, 4; **Fig. 3 time-course panels read from 400 dpi crops** (`f3a-5.png`, `f3b-5.png`, `f3c-5.png`),
read error ±0.04 mg/kg db (±0.05 in the a_w 0.430 panel, where symbols overlap below 400 s) and ±30 s.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** UMR QualiSud, Montpellier: CIRAD (Bassama, Brat, Hocine, Boulanger), Montpellier SupAgro (Bohuon),
Université Montpellier II (Günata). Corresponding author **Philippe Bohuon** (SupAgro). **Not in
`data/keys/laboratories.yml`: NEW**, proposed id **`montpellier_qualisud_bohuon`**, lab string
"CIRAD / SupAgro QualiSud (Bohuon)". No same-lab dossier in `data/lit/extraction_dossiers/` (grep bohuon / qualisud /
cirad: none). The same test cell is described in Bassama et al. 2010 Food Chem 123:558 (phenolics, model system) [C],
not on disk. **Cited as ref. 33 (kinetic scheme source) by Chan 2020** (`chan2020_extraction.md`).

## 1. What this paper answers for the model

| model need | what Bassama 2011 gives | tag |
|---|---|---|
| acrylamide formation / elimination k vs T in a food matrix | Table 2: k_Fref, k_Eref (s⁻¹, 170 °C), Ea_F, Ea_E, ± 95 % CI, 3 a_w | [M] |
| a_w dependence of k | three initial a_w (0.972 / 0.904 / 0.430); k_F ×8.1, k_E ×14 from 0.972 to 0.430 | [M] / [Z] |
| acrylamide vs time and T | Fig. 3: 3 a_w × 4 T × 7-9 times, mg/kg "dry basis" | [F] |
| precursor charge of the heated flour | Asn 32.3 ± 3.2, Glc 8610 ± 700, Fru 7420 ± 600, Suc 20400 ± 1600 mg/100 g db | [M] |
| Asn and sugar loss during heating | **not measured** (single response) | [!] |
| product temperature history | measured (centre thermocouple) and fed to the ODE solver, **not printed** | [!] |
| sorption isotherm of the matrix | Fig. 2 (25 °C, water addition), not digitised beyond the three working points | [F] |
| pH | **not measured** | [!] |
| Asn in 10 Musa cultivars; Asn vs ripening | Table 1; Fig. 1 | [M] / [F] |

**Answers.**
1. **Is k_F really a formation rate constant of the lane's kind? No** (confidence ~90 %). Eq. 1 sends every consumed
   Asn to acrylamide (dc_Asn/dt = −k_F c_Asn, dc_AA/dt = k_F c_Asn − k_E c_AA) and only acrylamide was fitted, so k_F
   is a lumped *yield-weighted* Asn → AA rate (pseudo-first-order in Asn, sugars ~600× in molar excess [Z]), not an
   Asn loss rate. At 170 °C, a_w 0.972, k_F = 0.96e-5 s⁻¹ implies 3 % Asn loss per hour; real Asn loss in such
   matrices is far faster at ~1 mol % yield. Comparable only to other single-response (Claeys-type) k_F, after
   conversion: the lane's Leuven k_F are second-order (M⁻¹ min⁻¹) in equimolar powders.
2. **Does a_w change the rates? Yes, strongly, and in the opposite direction from what the lane declares**
   (confidence ~75 % that the effect is real for this matrix; ~95 % that the printed numbers say so). k_F rises
   ×5.1 (0.972 → 0.904) and ×1.6 (0.904 → 0.430); k_E ×6.1 and ×2.4 [Z from Table 2]. Ea does not move within the
   CI (EaF 168-178, EaE 146-175 kJ/mol). The lane (`src/kinetic_core/acrylamide_conditions.py`) carries formation
   FLAT in a_w over 0.34-0.99 (band 0.41-1.39 × the 0.92 value) and elimination *falling* on the dry side
   (De Vleeschouwer 2007 k_E ratio 0.33-0.76 below 0.92). Bassama has both rising on the dry side, 0.972 inside the
   lane's window. Caveats: (a) only three a_w, one matrix, one flour; (b) a_w measured at 25 °C before heating;
   (c) 1:1000 sugar excess, which the authors argue removes the mobility limitation that governs equimolar powders
   [Q]; (d) the 0.972 arm has the worst fit at 140 °C (sec. 4c). This is a candidate B54 hold-out and a direct
   challenge to the declared flat formation multiplier, not a constant to adopt.
3. **Temperature dependence:** Ea_F 168-178 kJ/mol and Ea_E 146-175 kJ/mol [M] — in line with the Leuven
   single-response formation barriers (De Vleeschouwer 2007 Ea_F 154-182 kJ/mol, `devleeschouwer2007_extraction.md`;
   Claeys 2005 168 kJ/mol, `claeys2005_extraction.md`), and at the top of the Leuven Ea_E range (44-152). Fitted on
   non-isothermal T(t), so heat-up is in the fit (good) but the T(t) data are not available to re-run it.
4. **Are the constants internally consistent with Fig. 3? Yes, within ~20 %** (confidence ~75 %) [Z]: the
   quasi-steady plateau k_F/k_E × Asn0 × (71.08 / 132.12) with Asn0 = 323 mg/kg db gives 1.17 (a_w 0.972, 160 °C),
   1.01 (0.904, 180 °C) and 0.43 (0.430, 140 °C) mg/kg db against fitted-curve plateaus of ~1.15, ~1.24 and 0.48.
   This supports a molar conversion and "dry basis" = per kg dry matter (see flag 1).
5. **Asparagine in plantain** [M]: cooking bananas / plantains 0-73 mg/100 g db (Harton highest), dessert bananas
   272-321; ripening of barraganete drops Asn from ~86 to ~15-26 mg/100 g db in 5-11 days [F].

## 2. Conditions

- **Matrix** [M]: plantain "barraganete" (AAB), retail Montpellier, maturity stage 2 (green); 10 ± 0.2 mm slices,
  liquid-N2 frozen, freeze-dried (CRYO), ball-milled, stored in vacuum jar. **Which aliquot (green or ripened) made
  the kinetics flour is not stated** [!]; its Asn (32.3) sits between green (~86) and ripe (~15-26) in Fig. 1, and
  its sugars (36 % of dry matter) are high for green plantain [Z, ~60 % that it was partly ripened].
- **Precursors of the kinetics flour** [M]: Asn 32.3 ± 3.2; glucose 8610 ± 700; fructose 7420 ± 600; sucrose
  20400 ± 1600 mg/100 g db. Mass ratio sugars : Asn ≈ 1128 [Z]; molar ratio ≈ 610 (sucrose as one molecule) or ≈ 850
  (as two hexoses) [Z]; the paper's "molecular ratio … 1:1000" [V] is the mass ratio.
- **a_w** [M]: water added directly to the flour, 48 h equilibration in a vacuum jar, a_w (Aqualab, 25 °C) then
  heat. Targets 0.972 / 0.904 / 0.430 ↔ 70 / 41 / 8 % water wb. 12 h more equilibration shifts a_w by ~0.05 [M].
- **Vessel / heating** [M]: hermetically sealed test cell (design ref. Bassama 2010), ~10 g paste; 1.5 mm Type J
  thermocouple at the geometric centre of one compartment, logged; first bath silicone oil at 250 °C to speed
  heat-up, transfer < 3 s to a 5 L fryer at 140-200 °C, both stirred at 800 rpm; local T within 3 °C below T∞ [M].
  Come-up time and the T(t) curves are **not printed** [!]. Cooling in ice water, neglected in the model [M].
- **Times** [F]: ~3-45 min at 180 / 200 °C, ~5-60 min at 160 °C, ~5-90 min at 140 °C (Fig. 3; exact schedule not
  printed).
- **Headspace / O2 / pH**: not stated / not measured.
- **Replicates** [M]: five replicates only at 15 min (140, 160 °C) and 7 min (180, 200 °C); CV < 8 %. Other points
  appear single.

## 3. Analytical methods

- **Acrylamide** [M]: 1 g sample + 1 µg labelled acrylamide (printed "3H2-acrylamide"; m/z 75 → 58, i.e. d3
  [Z]), aqueous 0.2 % formic extraction, Carrez, two-step SPE (Isolute Multimode, ENV+), UPLC (Atlantis dC18) -
  ESI-MS/MS (TSQ Quantum Ultra) SRM 72 → 55 / 75 → 58. **Stable-isotope dilution**: calibrated. No LOD printed.
- **Sugars** [M]: 80 % / 50 % ethanol extraction, HPAE-PAD (Dionex DX600, CarboPac MA-1).
- **Asparagine** [M]: NaHCO3 pH 9 extraction, Dabsyl derivatisation, RP-HPLC (Krause 1995 method). ± = SE (n = 3).
- **Kinetics** [M]: Eq. 1 ODEs, Arrhenius Eq. 2 re-parameterised at T_ref 170 °C, Matlab ode15s fed with measured
  T(t); four parameters per a_w by nlinfit on acrylamide only; ± = 95 % asymptotic CI (nlparci).

## 4. Data

### 4a. Table 2 (p. 1457): fitted constants, T_ref = 170 °C, 140-200 °C, ± 95 % CI [M]

| a_w (25 °C) | k_Fref (×10⁻⁵ s⁻¹) | k_Eref (×10⁻² s⁻¹) | Ea_F (kJ/mol) | Ea_E (kJ/mol) |
|---|---|---|---|---|
| 0.972 | 0.96 ± 0.16 | 0.14 ± 0.04 | 177.6 ± 18.6 | 174.7 ± 21.4 |
| 0.904 | 4.90 ± 1.33 | 0.85 ± 0.26 | 168.5 ± 15.8 | 146.2 ± 16.8 |
| 0.430 | 7.74 ± 2.36 | 2.02 ± 0.64 | 177.8 ± 14.4 | 155.3 ± 15.1 |

Derived [Z]: in min⁻¹ at 170 °C, k_F = 5.8e-4 / 2.9e-3 / 4.6e-3 and k_E = 0.084 / 0.51 / 1.21. At 160 °C (the
Leuven T_ref), k_E = 0.028 / 0.20 / 0.53 min⁻¹ (a_w 0.972 / 0.904 / 0.430), vs Leuven 2008 k_Eref 0.05-0.11 min⁻¹ at
a_w 0.88-0.99 (`acrylamide_conditions.KE_TABLE`) and Leuven 2007 0.16-0.49 min⁻¹ at a_w 0.34-0.92 (`KE_2007`):
same order of magnitude; Bassama slower at 0.972, comparable at 0.904, faster at 0.430.
**Fig. 4** (Arrhenius plot, kF axis ×10⁻⁴ s⁻¹, kE axis ×10⁻² s⁻¹) agrees with Table 2: interpolated to 170 °C the
a_w 0.972 points read 0.098 ×10⁻⁴ and 0.14 ×10⁻² s⁻¹ [F, ±10 %]; the 200 °C / 140 °C kF ratio read off the a_w 0.972
line, ~670, returns Ea ≈ 176 kJ/mol [Z].

### 4b. Precursors

**Table 1, Asn (mg/100 g db)** [M] (± SE, n = 3 where printed): Africa (AAB plantain) n.d.; Dominico 13.32 ± 0.22;
Dominico harton 40.06 ± 1.20; Harton 73.32 ± 4.33; Cubano blanco 40.41 ± 0.66; Cavendish (AAA) 282.47 ± 12.68; Gros
michel 317.02; Rollizo 321.47; Fhia17 (AAAA) 271.41; Cachaco (ABB) 101.97. Flours oven-dried 40 °C / 48 h (Colombia).
**Fig. 1, barraganete Asn vs ripening (mg/100 g db)** [F, ±2]: green 86 (SE ~9), 4 d 62, 5 d 22, 6 d 26, 7 d 24, 8 d
15, 11 d 26.

### 4c. Fig. 3: acrylamide (mg/kg "dry basis") vs time (s) [F]

| a_w | T (°C) | (t s, AA) |
|---|---|---|
| 0.972 | 140 | (300, 0.00) (590, 0.04) (900, 0.07) (1810, 0.17) (3630, 0.29) (4540, 0.58) (5440, 0.37) |
| 0.972 | 160 | (300, 0.15) (590, 0.26) (900, 0.60) (1200, 0.55) (1810, 0.70) (2710, 0.96) (3630, 0.90) |
| 0.972 | 180 | (165, 0.41) (290, 0.95) (405, 0.87) (530, 1.16) (600, 1.21) (745, 1.34) (1200, 1.45) (1810, 1.61) (2710, 1.34) |
| 0.972 | 200 | (165, 0.99) (290, 1.36) (405, 1.45) (530, 1.38) (650, 1.27) (890, 1.20) (1200, 1.15) (1810, 1.08) (2710, 0.65) |
| 0.904 | 140 | (300, 0.11) (600, 0.19) (900, 0.27) (1820, 0.39) (3640, 0.56) (4540, 0.64) (5450, 0.71) |
| 0.904 | 160 | (300, 0.65) (600, 0.73) (900, 0.81) (1210, 0.90) (1820, 1.04) (2720, 1.03) (3640, 0.90) |
| 0.904 | 180 | (185, 0.55) (290, 0.94) (455, 1.24) (530, 1.41) (735, 1.32) (900, 1.47) (1820, 1.17) |
| 0.904 | 200 | (195, 1.30) (300, 1.48) (455, 1.19) (530, 1.10) (660, 0.87) (900, 0.59) (1210, 0.43) (1820, 0.23) |
| 0.430 | 140 | (290, 0.09) (600, 0.22) (900, 0.22) (1820, 0.44) (3650, 0.50) (4540, 0.50) (5450, 0.50) |
| 0.430 | 160 | (290, 0.49) (600, 0.65) (900, 0.72) (1220, 0.71) (1820, 0.60) (2730, 0.43) (3650, 0.35) |
| 0.430 | 180 | (185, 0.82) (290, 0.92) (410, 0.88) (530, 0.76) (720, 0.66) (900, 0.56) (1820, 0.15) |
| 0.430 | 200 | (185, 0.62) (290, ~0.89) (410, 0.60) (530, 0.35) (660, 0.31) (900, 0.16) (1210, 0.06) (1820, 0.04) |

All arms start at 0 at t = 0. Symbol assignment in the a_w 0.430 panel below 400 s is uncertain (□ and △ overlap;
the 200 °C fitted curve peaks at ~1.0 near 300 s). Fitted-curve features [F]: a_w 0.972 140 °C fit reaches only 0.26
at 5400 s (data 0.37-0.58: under-fit); a_w 0.430 160 °C fit stays ~0.6 to 3600 s (data fall to 0.35); 180 °C
elimination is poorly fitted at every a_w, as the authors say [Q].
Text maxima [M]: at 180 °C, maximum reached after 5 / 15 / 30 min for a_w 0.430 / 0.904 / 0.972 (Fig. 3 agrees:
~290-410 s, 900 s, 1810 s). Conversion to wet basis [Z]: maxima 0.48 / 0.87 / 0.85 mg/kg wb (0.972 / 0.904 / 0.430),
matching the abstract's "max ~ 0.9 ppm, wb" [V]. Molar yield at the a_w 0.972 maximum ≈ 0.93 mol % of Asn0 [Z].

### 4d. Fig. 2 isotherm (25 °C) [F]
Water addition (this work, to a_w 0.7 only with error bars) vs Johnson & Brennan 2000 (sulfuric-acid method, to
a_w 0.99); close to a_w 0.7; at a_w 0.43 ≈ 8 g/100 g wb, 0.90 ≈ 41, 0.97 ≈ 68-70 — consistent with sec. 2.4. Not
digitised further (the plant-protein isotherm B54 needs is a different matrix).

## 5. Flags and internal consistency

1. [!] **"Dry basis" is defined oddly**: the text says concentrations were divided by "the amount of water initially
   added" [V paraphrase], which would make the a_w 0.430 values per kg of 8 % water. The wb conversion of the abstract
   and the Table 2 ↔ Fig. 3 plateau check (sec. 1 ans. 4) both fit "per kg dry matter" instead (~75 % confidence).
2. [!] **Parity plots (Fig. 3b) run to CAA_exp ≈ 3.2 / 2.7 / 2.0** whereas no time-course point exceeds 1.61 / 1.48 /
   0.92: the parity panels carry points or a scaling absent from the time courses (unexplained).
3. [!] Text p. 1457: at a_w 0.904 maximum "after 5 min at 200 °C instead of 9 min at 180 °C" [V], but the next
   paragraph (and Fig. 3) give 15 min at 180 °C for a_w 0.904.
4. [!] Text "×10 and ×15" (k_F, k_E, 0.972 → 0.430): Table 2 gives ×8.1 and ×14.4 [Z]. "×1000 with T (140 → 200 °C)":
   Ea 177 kJ/mol gives ×690, Ea_E 146-175 gives ×300-650 [Z]; "×10,000 global" is the product of rounded factors.
5. [!] T(t) measured but not printed; the fit cannot be re-run isothermally without it. Come-up time unknown
   (250 °C boost bath; probably < 1-2 min but unverified).
6. [!] Initial a_w at 25 °C only; the sealed cell keeps water but a_w at 140-200 °C differs; whether water is lost to
   cell headspace is not stated.
7. [!] Kinetics flour's maturity not stated (sec. 2); precursor values are for that flour only.
8. [!] Most kinetic points single; replicate CV < 8 % at four (T, t) cells only; residuals larger at t < 20 min [Q].
9. [!] Model ignores sugars (justified by the 1:1000 excess) — the opposite reduction to Chan 2020 (first order in
   reducing sugar, Asn ignored), which cites this paper as a kinetic source.

## 6. Existing records in the repo that cite this paper

`grep -rIn` of `data/ docs/ src/ scripts/ tasks/` (2026-11-03) for `foodres.2011.03.018`, `bassama`, `plantain`,
`bohuon`, `gunata`/`günata`: **no hits** outside the triage line in `data/keys/pdf_dossier_aliases.yml` (C9). The
acrylamide lane does not use this paper; nothing to check field by field. The a_w code it bears on is
`src/kinetic_core/acrylamide_conditions.py` (AW_WINDOW 0.34-0.99, AW_MULTIPLIER 1.0, AW_ELIMINATION_TABLE from
De Vleeschouwer 2007) — see sec. 1 ans. 2.

## 7. Proposed ledger entries (acrylamide.yml; NOT written)

```yaml
  - id: bassama2011_tab2_acrylamide_kf_ke_ea_aw_plantain
    dossier: bassama2011_extraction.md
    paper: "Bassama, Brat, Bohuon, Hocine, Boulanger & Gunata 2011 Food Res Int 44:1452"
    lab: "CIRAD / SupAgro QualiSud (Bohuon)"
    system: "freeze-dried plantain flour (Musa AAB barraganete; Asn 32.3, Glc 8610, Fru 7420, Suc 20400 mg/100 g db), rehydrated to a_w 0.972 / 0.904 / 0.430"
    medium: "~10 g paste in a hermetically sealed steel test cell, 250 C boost bath then stirred oil bath; measured centre T(t)"
    temperature_c: [140, 160, 180, 200]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [0, 90]
    species: [acrylamide]
    kind: rate_constant
    location: "Table 2: k_Fref (1e-5 /s), k_Eref (1e-2 /s) at T_ref 170 C, Ea_F, Ea_E, +- 95% asymptotic CI, per initial a_w (dossier sec. 4a)"
    status: digitised_in_dossier
    quant_class: sida
    simulable: false
    condition_quote: "Three initial water activities were chosen in order to perform experimentations: 0.972, 0.904, and 0.430."
    flags: [matrix_food, aw_measured_25c_before_heating, single_response_fit, asn_and_sugars_t0_only, temperature_history_not_printed, "complements: bassama2011_fig3_acrylamide_timecourses_aw"]
    note: "Lumped Asn->AA (k_F, first order in Asn, sugars ~600x molar excess) and AA->DP (k_E). k_F x8, k_E x14 from a_w 0.972 to 0.430, Ea flat: contradicts the lane's flat formation / dry-side-falling elimination (acrylamide_conditions.py, Leuven powders). B54 hold-out candidate, not a constant to adopt."
  - id: bassama2011_fig3_acrylamide_timecourses_aw
    dossier: bassama2011_extraction.md
    paper: "Bassama, Brat, Bohuon, Hocine, Boulanger & Gunata 2011 Food Res Int 44:1452"
    lab: "CIRAD / SupAgro QualiSud (Bohuon)"
    system: "freeze-dried plantain flour rehydrated to a_w 0.972 / 0.904 / 0.430 (70 / 41 / 8 % water wb)"
    medium: "sealed steel test cell in oil baths; measured centre T(t), not printed"
    temperature_c: [140, 160, 180, 200]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [0, 90]
    species: [acrylamide]
    kind: time_series
    location: "Fig. 3 left panels, 3 a_w x 4 T x 7-9 times, mg/kg dry basis, read at 400 dpi +-0.04 (dossier sec. 4c)"
    status: digitised_in_dossier
    quant_class: sida
    simulable: false
    condition_quote: "The samples (plantain paste at different moisture content, weight ~10 g) were introduced into the hermetically sealed test cell"
    flags: [matrix_food, aw_measured_25c_before_heating, heating_up_included_in_time, temperature_history_not_printed, no_replicates, dry_basis_definition_ambiguous]
    note: "Rise-then-fall at 180/200 C in all a_w. Peak 1.61 / 1.48 / 0.92 mg/kg db (0.972 / 0.904 / 0.430) = ~0.5-0.9 mg/kg wb, ~0.9 mol% of Asn0. Not simulable until the lane has a plantain/starch matrix and the non-isothermal T(t); usable as a directional a_w check (peak time 30 / 15 / 5 min at 180 C)."
  - id: bassama2011_tab1_fig1_asparagine_musa_ripening
    dossier: bassama2011_extraction.md
    paper: "Bassama, Brat, Bohuon, Hocine, Boulanger & Gunata 2011 Food Res Int 44:1452"
    lab: "CIRAD / SupAgro QualiSud (Bohuon)"
    system: "10 Musa cultivars (Colombia, 40 C oven-dried flours) and barraganete plantain ripened 0-11 d at 17 C after ethylene"
    medium: "unheated flours"
    temperature_c: []
    ph: {initial: null, measured_during: false, final: null}
    time_min: []
    species: [asparagine, glucose, fructose, sucrose]
    kind: endpoint
    location: "Table 1 (Asn, 10 cultivars), Fig. 1 (Asn vs ripening day, read +-2), sec. 3.3 text (Asn, Glc, Fru, Suc of the kinetics flour)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "the results were 32.3 +- 3.2, 8610 +- 700, 7420 +- 600 and 20400 +- 1600 mg/100 g dry basis"
    flags: [matrix_food, not_a_kinetic_dataset]
    note: "Precursor composition only; the sec. 3.3 values are the initial charge of the kinetics pots (which maturity aliquot is not stated)."
```

## 8. Supplementary material
None referenced (no SI, appendix or supplementary table in the paper). The only model-relevant missing item is the
measured T(t) record of each run, which the paper does not publish.
