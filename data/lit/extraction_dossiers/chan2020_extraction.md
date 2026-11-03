# Chan 2020 — EXTRACTION (file stem `chan2020`)
DOI: 10.3390/foods9020200 (printed p. 1 footer). Crossref 2026-11-03: "Computer Simulation with a Temperature-Step
Frying Approach to Mitigate Acrylamide Formation in French Fries", Foods **9**(2):200, issued 2020-02-16, **single
author Der-Sheng Chan**. **Matches the PDF.**
[!] **Naming caveat: the brief and the triage line say "Chan et al."; the paper has one author.** Not to be confused
with `chan1994_extraction.md` / `chan1994b_extraction.md` (Chan F. & Reineccius, sulfur kinetics).

### A 2-D finite-element (COMSOL 5.1) heat + moisture transfer model of one 8.5 × 8.5 mm potato strip frying at 170 °C, coupled to a lumped two-step acrylamide scheme (reducing sugar → AA → D, first order, switched on only at T ≥ 120 °C), then used to simulate a temperature-step (170 → 150 °C after 3 or 4 min) process. **No experiment of its own**: the "validation" re-plots Gökmen, Palazoğlu & Şenyuva 2006 (J Food Eng 77:972) — water, surface / core T and acrylamide at 170 °C — against parameters the author marks "estimated". **No asparagine, no fitting method, no uncertainties; the step-frying predictions are not tested.**

**Read 2026-11-03 from `data/articles/chan2020.pdf`** (15 pp., MDPI, pdfTeX; the text layer interleaves a
"Foods 2019, 8, x FOR PEER REVIEW" draft with the final text, so the text dump is unreliable). Prose from
`pdftotext -layout` / `-raw`; **Table 1 and Figs. 2-4 read from 150 dpi renders of pp. 7-8** (scratchpad
`ingest/chan2020_work/p-07.png`, `p-08.png`), read error ±1 % water, ±2 °C, ±10 µg/kg, ±0.05 min. Figs. 5-10
(simulation only) not rendered; their printed text values are quoted.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Department of Information Technology, **Lee-Ming Institute of Technology**, New Taipei City, Taiwan;
sole and corresponding author **Der-Sheng Chan**. **Not in `data/keys/laboratories.yml`: NEW**, proposed id
**`leeming_chan_ds`** ("Lee-Ming Inst. Tech. (Chan DS)") — but it is a modelling-only source; **every measured number
belongs to the Hacettepe group (Gökmen)**, ledger string "Hacettepe (Gokmen)" (used by `gokmen2012_*` entries; also
not yet in `laboratories.yml`). **Same person, other group:** Der-Sheng Chan is second author of Lin, Chan, Kao &
Sung 2021 Polymers 13:1901 (`lin2021_extraction.md`, ledger lab "NTOU (Sung)", glucose/fructose + Asn acrylamide
with chitosan at 180 °C); since Chan 2020 contributes no measurement, no lab merge is needed.

## 1. What this paper answers for the model

| model need | what Chan 2020 gives | tag |
|---|---|---|
| acrylamide kinetics in a fried potato matrix | k0 = 22 s⁻¹ (shared by both steps), Ea,AA 44.5, Ea,D 21.3 kJ/mol, all "estimated" | [M] |
| which precursor is rate-limiting | first order in reducing sugar (C_Ai = 6.8 mol/m³), asparagine absent | [M] |
| measured acrylamide vs time | 4 points at 170 °C (0, 3, 6, 9 min), re-plotted from Gökmen 2006 | [F, secondary] |
| measured product temperature | surface and core T, 9 points each, 0-9 min, 170 °C (Gökmen 2006) | [F, secondary] |
| measured water loss | average water content, 4 points (Gökmen 2006) | [F, secondary] |
| heat / mass transport parameters | h_t 105 W/m²K, k_m 5.2e-4 m/s, D0 5.0e-6 m²/s, Ea,dif 27.6 kJ/mol, k_evp 8e-3 s⁻¹ (all "estimated") | [M] |
| validation of the step-frying prediction | **none** | [!] |
| a_w, pH, Asn, sugar loss | **none** | [!] |

**Answers.**
1. **Which kinetic scheme and constants come from the literature? The scheme only; no constant.** (confidence ~90 %).
   Eq. 17 (A + B → AA → D) is attributed to refs 25, 29, 33, 34 [C] = Palazoğlu & Gökmen 2008 J Food Sci 73:E109
   (10.1111/j.1750-3841.2008.00676.x), Carrieri et al. 2009 J Food Eng 95:90, **Bassama et al. 2011**
   (`bassama2011_extraction.md`) and Corradini & Peleg 2006 Crit Rev Food Sci Nutr 46:489. First order in reducing
   sugar is credited to refs 15, 25 (Palazoğlu & Gökmen 2008 JAFC 56:6162, 10.1021/jf073046l, and J Food Sci). The
   120 °C onset (S_AA) is credited to ref. 23 (Keramat review). **Parker 2012, Knol, Wedzicha, Balagiannis are not
   cited**: the Reading multi-step French-fry scheme (`parker2012b`, `balagiannis2019`, sibling dossiers) is not
   used. Of Table 1's 16 inputs, only C_Ai and H ([15]), C_wi ([30], Pedreschi 2008), H_evp ([36]) and T_i (set)
   are sourced; Ea,AA, Ea,D, k0 and every transport parameter are footnote "e", "estimated" [M], with no
   estimation procedure, objective or uncertainty described [!].
2. **Is the model validated on measured T and acrylamide? Compared, not validated** (confidence ~80 %). The
   comparison is one condition (170 °C) of one external data set (Gökmen 2006, ref. 20), and the 12 "estimated"
   parameters were almost certainly tuned to that same data set, so the agreement is a calibration. Quality
   of the agreement [Z from sec. 4]: water at 3 min 42 % measured vs ~50 % simulated (worst point); surface T at
   9 min 141.5 measured vs ~148 simulated; core T well reproduced (~103 °C plateau); acrylamide 103 / 250 / 417 µg/kg
   measured vs ~55 / 250 / ~452 simulated at 3 / 6 / 9 min. The 150 °C step predictions (178 and 197 µg/kg) have no
   measurement, although ref. 15 (Palazoğlu & Gökmen 2008, temperature-programmed frying) would supply one [C,
   content not verified].
3. **Are the constants usable by the lane? No** (confidence ~95 %). (a) Single shared k0 for formation and
   degradation with Ea,D < Ea,AA makes k2 = 0.068 s⁻¹ at 170 °C (AA lifetime ~15 s) and k1 = 1.25e-4 s⁻¹ [Z], so AA
   is always at quasi-steady state k1/k2 · C_A, set only by ΔEa = 23.2 kJ/mol; (b) no Asn, no intermediate; (c) a hard
   120 °C threshold; (d) the constants are tuned, undocumented, and the units of the AA output (mol/m³ → µg/kg) are
   not given. Their only use is as a record that a sugar-limited lumped model with a 120 °C switch can match the
   Gökmen 170 °C curve.
4. **The useful content is the re-plotted Gökmen 2006 data**, which belong in the ledger under the primary paper
   (10.1016/j.jfoodeng.2005.08.030, J Food Eng 77:972; not on disk) once it is downloaded: one fry, one T, with surface
   and core T histories — a candidate frying hold-out for the measured-T(t) path, like `knol2009_*`.

## 2. Conditions (of the re-plotted experiment; as described by Chan, not by the primary)

- **Product** [M]: potato strips 8.5 × 8.5 × 70 mm; initial composition 77.2 % water, 18.8 % carbohydrate, 2.7 %
  protein (ref. 30); C_wi = 35,778 mol/m³; reducing sugar C_Ai = 6.8 mol/m³ "calculated from [15]" (cultivar,
  Asn, pre-treatment not given). Oil type ignored by assumption [M].
- **Frying** [M]: 170 °C, 0-9 min; single strip per fry in the source studies, so no oil T drop [Q]. T_i 293.15 K.
- **Model domain** [M]: half cross-section (x, y), called "axisymmetrically two-dimensional" [V] though the
  equations are Cartesian [!]; volume change neglected; evaporation switch S_W at 103 °C (Eq. 1), though Eq. 4's text
  says ≥ 373.15 K [!]; AA switch S_AA at 120 °C.
- **Step frying (simulated only)** [M]: oil 170 °C for t_H = 3 or 4 min, then 150 °C (f = 0.95489 on the kelvin
  temperature).
- **pH, a_w, headspace, replicates**: not applicable / not given.

## 3. Analytical methods
None in this paper. The acrylamide, temperature and water methods of Gökmen 2006 are not described here (quant
class unknown from this source). Model: COMSOL Multiphysics 5.1, FEM; mixing-rule properties (density, c_p,
conductivity) from Choi-Okos-type polynomials (ref. 35) [M].

## 4. Data

### 4a. Table 1 (p. 7), input parameters [M]; footnotes a calculated from [15], b from [30], d from [36], e estimated, s set

| par. | value | src | par. | value | src |
|---|---|---|---|---|---|
| C_Ai | 6.8 mol/m³ | a | H_evp | 2.3 × 10⁶ J/kg | d |
| C_eq | 9265 mol/m³ | e | h_t | 105 W/(m² K) | e |
| C_wi | 35,778 mol/m³ | b | k_evp | 8.0 × 10⁻³ s⁻¹ | e |
| D_0 | 5.0 × 10⁻⁶ m²/s | e | k_m | 5.2 × 10⁻⁴ m/s | e |
| Ea,AA | 44.5 kJ/mol | e | k0 | 22 s⁻¹ | e |
| Ea,D | 21.3 kJ/mol | e | M_A | 0.172 kg/mol | e |
| Ea,dif | 27.6 kJ/mol | e | W | 8.5 × 10⁻³ m | e |
| H | 8.5 × 10⁻³ m | a | T_i | 293.15 K | s |

Derived [Z]: k1(170 °C) = 22·exp(−44500/(8.314·443.15)) = 1.25 × 10⁻⁴ s⁻¹; k2(170 °C) = 0.068 s⁻¹; k1/k2 = 1.8e-3
at 170, 8.3e-4 at 120 °C. M_A = 0.172 kg/mol is labelled molecular weight of A (reducing sugar per the
nomenclature) but is neither glucose / fructose (0.180) nor acrylamide (0.071) [!].

### 4b. Re-plotted Gökmen 2006 measurements at 170 °C, with Chan's simulated curve [F]

| t (min) | water (% wb) meas. / sim. | surface T (°C) meas. | core T (°C) meas. | AA (µg/kg) meas. / sim. |
|---|---|---|---|---|
| 0 | 76.5 / 76 | 20 | 20 | 0 / 0 |
| 0.16 | — | 77 | 25 | — |
| 0.5 | — | 104 | 65 | — |
| 0.84 | — | 117 | 91 | — |
| 1.18 | — | 124 | 101 | — |
| 1.52 | — | 125 | 104 | — |
| 3 | 42.0 / ~50 | 128 | 103 | 103 / ~55 |
| 6 | 32.4 / ~32 | 136 | 103 | 250 / ~250 |
| 9 | 24.9 / ~20 | 141.5 | 104 | 417 / ~452 |

Text values [M]: water "32.4 and 25.1 %" at 6 and 9 min (= the measured points, though worded as model output);
final surface T "148 °C" (= the simulated curve, measured ~141.5); core ~103 °C after 1.5 min; acrylamide "418 µg/kg"
at 9 min (= the measured point; the curve ends ~452); lag 1.5 min. Basis of µg/kg (wet / dry) not stated [!].

### 4c. Simulation-only outputs (Figs. 5-10; text values) [M]
Step frying: final average water 36.8 % (t_H 3 min) and 35.0 % (4 min); average strip T plateau 113 °C; AA 178 and
197 µg/kg at 9 min. 9 min contour ranges: T 103-160 °C (no step) vs 103-136 °C (step); AA 0-2950 vs 0-1950 µg/kg,
highest at the outer corner crust.

## 5. Flags and internal consistency

1. [!] **Reduction stated three ways**: abstract "decreased … to 57 % of" [V]; sec. 3.2 "57 % (418→178,
   (418 − 178)/418)"; conclusion "a decrease of 43 % (418→178, 178/418)". Arithmetic: 178/418 = 42.6 % remaining,
   a 57.4 % reduction [Z]. The abstract and conclusion are both wrong in wording.
2. [!] **Mixed bases in that comparison**: 418 µg/kg is the measured 9 min point (simulated ~452) while 178 is
   simulated; model-vs-model gives 178/452 = 39 % remaining (~65 % confidence that 418 is the measured value).
3. [!] Step frying: average strip T "constant at 113 °C" below the 120 °C onset, yet 178-197 µg/kg forms (crust above
   120 °C); and water: "only 3 %-4 % higher" with the step vs 36.8 / 35.0 against 25.1 % (9 min, no step) = 10-12
   points; a later sentence gives "27.6 % … higher … (19.6 %)" [V], unclear.
4. [!] Shared k0 for k1 and k2 (Eqs. 21-22, Table 1); Ea,D (21.3) below Ea,AA (44.5) — the opposite of Bassama's
   ordering and far below every measured elimination Ea in the corpus (44-175 kJ/mol) [Z].
5. [!] [Z] Quasi-steady AA = (k1/k2)·C_A·M_AA/ρ ≈ 0.4-0.8 mg/kg at 120-160 °C (ρ ~1100 kg/m³) — enough for the
   418 µg/kg average but ×4 below the 2950 µg/kg corner maximum of Fig. 10; the mol/m³ → µg/kg conversion (M_A = 0.172?
   dry basis?) is unstated, so the contours cannot be reproduced from Table 1 (~60 % that an unstated factor enters).
6. [!] W is "estimated" while H, the same 8.5 mm, is "calculated from [15]".
7. [!] The draft text layer is embedded beneath the final text (MDPI PDF); quote from page renders only.

## 6. Existing records in the repo that cite this paper

`grep -rIn` of `data/ docs/ src/ scripts/ tasks/` (2026-11-03) for `foods9020200`, `chan2020`, "Chan et al." /
"Chan, D", `Der-Sheng`, `temperature-step`: **no hits** outside the triage line in `data/keys/pdf_dossier_aliases.yml`
(C6) and the author string of `lin2021_extraction.md` (a different paper, sec. Laboratory). Its data sources (Gökmen 2006 `jfoodeng.2005.08.030`, Palazoğlu & Gökmen 2008 `jf073046l` /
`1750-3841.2008.00676.x`) are also uncited in the repo and not on disk. The acrylamide lane does not use this paper;
nothing to check field by field.

## 7. Proposed ledger entries (acrylamide.yml; NOT written)

```yaml
  - id: chan2020_figs2_4_gokmen2006_fries_170c_replot
    dossier: chan2020_extraction.md
    paper: "Chan 2020 Foods 9:200 (10.3390/foods9020200), re-plotting Gokmen, Palazoglu & Senyuva 2006 J Food Eng 77:972 (10.1016/j.jfoodeng.2005.08.030)"
    lab: "Hacettepe (Gokmen)"
    system: "potato strips 8.5 x 8.5 x 70 mm, 77.2 % water, reducing sugar 6.8 mol/m3 (as computed by Chan from ref. 15)"
    medium: "deep-fat frying, 170 C, single strip per fry"
    temperature_c: [170]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [0, 9]
    species: [acrylamide, water content, surface temperature, core temperature]
    kind: time_series
    location: "Fig. 2 (water, 4 pts), Fig. 3 (surface and core T, 9 pts each), Fig. 4 (acrylamide, 4 pts), 150 dpi (dossier sec. 4b)"
    status: digitised_in_dossier
    quant_class: unknown
    simulable: false
    condition_quote: "The result of mathematical modeling are compared with experimental data obtained from the study of Gokmen et al."
    flags: [secondary_source, matrix_food, heating_up_included_in_time, low_moisture_no_aw_term, concentrations_not_stated]
    note: "Measured T(t) at surface and core plus acrylamide and water at one frying T; re-key to the primary (Gokmen 2006) when on disk. Asn not given; basis of ug/kg not stated. Candidate frying hold-out on the measured-T(t) path (cf. knol2009_fig2_tab1_crisp_acrylamide_frying)."
  - id: chan2020_tab1_lumped_sugar_acrylamide_constants
    dossier: chan2020_extraction.md
    paper: "Chan 2020 Foods 9:200 (10.3390/foods9020200)"
    lab: "Lee-Ming Inst. Tech. (Chan DS)"
    system: "French-fry strip model: reducing sugar -> acrylamide -> degraded products, first order, active at T >= 120 C"
    medium: "2-D FEM heat and moisture transfer model (COMSOL 5.1); no experiment"
    temperature_c: [170]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [0, 9]
    species: [acrylamide, reducing sugar]
    kind: rate_constant
    location: "Table 1: k0 22 /s (shared), Ea,AA 44.5, Ea,D 21.3 kJ/mol, all footnote e 'estimated' (dossier sec. 4a)"
    status: text_only
    quant_class: other
    simulable: false
    condition_quote: "The acrylamide formation can be described using the first-order reaction with respect to reducing sugar"
    flags: [model_estimated_not_measured, no_uncertainty, secondary_source]
    note: "Not a measurement: tuned to the Gokmen 2006 170 C curve, estimation method not described, one shared k0, no Asn. Recorded so it is not mistaken for a measured French-fry constant; do not fit or check against it."
```

## 8. Supplementary material
None referenced (no SI, appendix or data availability statement beyond the article). The model-relevant data are in
the primary sources (Gökmen et al. 2006, ref. 20; Palazoğlu & Gökmen 2008, refs. 15 and 25), not in any SI.
