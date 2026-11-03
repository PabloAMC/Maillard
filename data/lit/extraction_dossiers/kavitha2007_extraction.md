# Kavitha & Modi 2007 — EXTRACTION (file stem `kavitha2007`)
DOI: 10.1016/j.lwt.2006.07.014 (printed p. 1280). Crossref 2026-11-03: "Effect of water activity and temperature on
degradation of 5′-inosine monophosphate in a meat model system", authors Kavitha, Modi, *LWT - Food Science and
Technology* **40**(7):1280-1286, issued 2007-09. **Matches the PDF and the brief.** No naming caveat.

### IMP (120 mg/100 g added, ~126 mg/100 g total) in washed, freeze-dried chicken breast or leg "meat fibre", brought to nominal a_w 0.70 / 0.80 / 0.90 with 0.2 M pyrophosphate pH 5.8, 2 g in 5 mL screw-cap glass tubes, heated **1 h** at 80 °C (water bath) or 120 °C (bench-top autoclave). IMP, inosine, hypoxanthine by HPLC-UV 254 nm, mg/100 g, n = 4. **One time point, two temperatures; no rate constant, no Ea, no ribose / R5P, no recovery control.**

**Read 2026-11-03 from `data/articles/kavitha2007.pdf`** (7 pp., Elsevier "Article in Press" proof layout,
born-digital). Prose from `pdftotext -layout`; Tables 1-3 (p. 1283) from a 150 dpi render; the a_w recipe (p. 1282)
from a 300 dpi crop (the text layer drops "µ": the volumes are **µl**, not ml); **Fig. 1 (p. 1284) and Fig. 2
(p. 1285) read from 200 dpi page renders and 200 dpi crops of the A panels** (scratchpad
`ingest/kavitha2007_work/pg-5.png`, `pg-6.png`, `f1a-5.png`, `f2a-6.png`). Bar read error ±0.1 mg/100 g on the IMP
axes (0-8, 0-10), ±0.02 on the 0-2.5 axes, ±0.01 on inosine axes (0-0.6, 0-0.25).

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Department of Meat, Fish and Poultry Technology, **Central Food Technological Research Institute
(CFTRI), Mysore**, India; corresponding author **V.K. Modi**. **Not in `data/keys/laboratories.yml`: NEW**, proposed id
**`cftri_mysore_modi`**, strings "CFTRI Mysore (Modi)" and the ledger's existing "CFTRI Mysore / Nottingham (Taylor)"
(Modi 2008, `modi2008_extraction.md`: Modi corresponding, same washed-meat-fibre + 0.2 M pyrophosphate pH 5.8 a_w
recipe; the Nottingham co-authors were Modi's host lab). Same-lab cited work not on disk: Vani, Modi, Kavitha et al.
2005 LWT 39:627 (IMP degradation vs pH and T, aqueous and fibre) — the method source for the fibre and the HPLC.

## 1. What this paper answers for the model

| model need | what Kavitha & Modi 2007 give | tag |
|---|---|---|
| IMP loss rate k(T, a_w) | **none printed**; one 1-h endpoint per cell; ~93-98 % lost in every spiked cell | [F] |
| order of IMP loss | not testable (one time point) | [!] |
| Ea | **none printed**; two-point "Ea" from my first-order k is 7-10 kJ/mol | [Z] |
| products | inosine and hypoxanthine (mg/100 g, Figs. 1-2); **ribose and ribose-5-phosphate not measured** | [F] |
| mass balance | inosine + hypoxanthine recover **2-4 mol %** of the IMP lost | [Z] |
| a_w as a variable | nominal a_w 0.70 / 0.80 / 0.90 by buffer volume; measured a_w of the mixes not reported | [M] / [!] |
| moisture basis | derivable: ~0.24 / 0.26 / 0.33 g water per g dry fibre | [Z] |
| pH | buffer pH 5.8; fibre pH 6.21 (breast) / 6.43 (leg) in 1:9 water; no pH after heating | [M] |
| chicken precursor composition | Table 3: IMP, inosine, hypoxanthine in fresh muscle and washed fibre (mg/100 g) | [M] |

**Answers.**
1. **k vs a_w and T: the paper gives no k** (confidence ~99 %). Taking the residual IMP at 1 h against the nominal
   start (126 / 127 mg/100 g) and assuming first order, k = 2.7-3.1 h⁻¹ at 80 °C and 3.7-4.3 h⁻¹ at 120 °C [Z]
   (sec. 4c). The a_w effect within a temperature is ≤ ±7 % in k (spiked arms); the 80 → 120 °C effect is ×1.27-1.41.
   These are not rate constants of a chemical step: at 93-98 % conversion a single point cannot separate first-order
   decay from a fast loss plus a floor, and the un-spiked controls give **lower** apparent k than the spiked arms
   (breast 80 °C: 1.8-2.2 vs 3.0-3.1 h⁻¹), which first-order kinetics forbids [Z].
2. **Ea: none printed; the derivable 7-10 kJ/mol is not an activation energy** (confidence ~90 %). Aqueous IMP
   hydrolysis (Matoba et al. 1988, 0.4 mM, sealed tubes; `matoba1988_extraction.md` sec. 4a [C]): k = 0.010 h⁻¹ at
   80 °C pH 4 (measured), ~0.005 h⁻¹ at 80 °C pH 7 and 0.43-0.50 h⁻¹ at 120 °C (extrapolated with Ea 114-128 kJ/mol).
   At the fibre's pH ~5.8-6.4 that predicts ~0.5-1 % IMP loss in 1 h at 80 °C [Z]; this paper reports ~95 %, **k
   300-600 × faster** at 80 °C and ~8-10 × at 120 °C. Water-limited hydrolysis at
   a_w 0.7 should be slower, not faster. Most likely causes (no control in the paper tests them): incomplete
   perchloric-acid recovery of IMP from heat-denatured fibre (no spiked-unheated recovery is reported; the "initial"
   125.99 mg/100 g is nominal, natural + added) ~55 %; residual muscle phosphatase / nucleosidase activity during
   preparation or heat-up ~25 %; genuine fast matrix-catalysed hydrolysis ~10 %; other ~10 %.
3. **Products: inosine (≤ 0.44 mg/100 g) and hypoxanthine (0.8-1.5 mg/100 g) only** (confidence ~95 %). Hypoxanthine >
   inosine in every spiked arm, consistent with direct N-glycosidic cleavage (IMP → Hx + R5P) rather than
   dephosphorylation first, but the 2-4 mol % balance means the products say little about where the IMP went.
4. **Usable for the a_w dependence the model lacks (B54)? No** (confidence ~85 %). Reasons: IMP is not a trunk
   species; a_w is nominal (recipe, not measured on the mix); one endpoint near full conversion saturates any a_w
   signal; recovery is unverified; and the inosine / hypoxanthine panels of Fig. 1A and Fig. 2A are the same data
   (sec. 5, flag 1). At most a qualitative sign: "a_w 0.7-0.9 changes 1-h IMP loss in chicken fibre by < 10 %".
   P(changes any model parameter) < 3 %. For B54's concentration-basis question the derivable moisture (0.24-0.33 g/g
   d.b. for a_w 0.7-0.9) is one weak point of a chicken-myofibril sorption curve at 25 °C, nothing more.
5. **The registry entry that cites this paper carries values the paper does not contain** (sec. 6; confidence ~98 %).

## 2. Conditions

- **Matrix** [M]: Indian country hen (*Gallus gallus*), 8 birds per replicate, breast and leg minced (4.5 mm), washed
  **8 ×** in 0.1 M phosphate pH 6.0 (OD525 of the washings, Table 1), centrifuged 6000 g, freeze-dried (−40 °C,
  0.3 mbar), powdered, stored −20 °C (text "20 1C", sign lost). Four separate fibre batches = the 4 replicates.
- **Fibre composition** (Table 2, g/100 g, mean ± SD n = 4) [M]: breast moisture 3.25 ± 0.10, fat 4.83 ± 0.05,
  protein 86.5 ± 0.48, ash 2.80 ± 0.22; leg 3.50 ± 0.08, 5.70 ± 0.08, 85.63 ± 0.50, 2.43 ± 0.05. Fibre a_w 0.24 ± 0.03
  (breast), 0.29 ± 0.04 (leg); pH (10 g in 90 mL water) 6.21 ± 0.06 / 6.43 ± 0.05; muscle pH 5.56 / 5.49.
- **Model system** [M]: 7 g fibre + **1400 µl** IMP solution (6 mg/mL in deionised water, pH 5.5; salt form not
  stated) → 120 mg IMP/100 g; + **20 / 110 / 620 µl** 0.2 M pyrophosphate pH 5.8 for a_w 0.7 / 0.8 / 0.9. Controls:
  1400 µl water in place of IMP solution. Mixed under liquid N2 in a mortar, 2 g per 5 mL screw-cap tube, −20 °C.
  [Z] water ≈ 1.65 / 1.74 / 2.25 g per ~6.8 g dry solids = 0.24 / 0.26 / 0.33 g/g d.b. (19.6 / 20.5 / 24.9 % w.b.).
- **Heating** [M]: 1 h, 80 °C water bath or 120 °C bench-top autoclave (J.P. Selecta); come-up time not reported;
  "no loss in weight" after heating. Headspace / O2: air in the tube (not stated otherwise).
- **Replicates** [M]: 4 fibre batches × each cell; ANOVA + Duncan (STATISTICA).

## 3. Analytical method

- Extraction [M] (Ozogul et al. 2000 [C]): 2 g in 15 mL 0.6 M HClO4, Polytron, 3000 g, Whatman 41; filtrate to pH 5.5
  with 6 M KOH, KClO4 removed; purine internal standard 0.1 mL (0.004 mM) added **after** extraction (so it corrects
  injection, not recovery).
- HPLC [M]: Shimadzu LC-10A, UV 254 nm, Supelcosil LC-18-T 15 cm × 4.6 mm 3 µm, 0.1 M phosphate pH 6.0 / methanol
  gradient (18 min), 1.5 mL/min, 20 µl. External standards 0.20 mM IMP, inosine, hypoxanthine; "concentration was
  calculated using the area for each peak" [V]. **Calibrated (external standard)**; no LOD, no recovery, no
  spiked-unheated sample [!]. ± = SD, n = 4.

## 4. Data

### 4a. Table 3 (p. 1283), mg/100 g, mean ± SD (n = 4) [M]

| | breast fresh | breast fibre | leg fresh | leg fibre |
|---|---|---|---|---|
| IMP | 166.25 ± 1.26 | 5.99 ± 0.24 | 125.00 ± 3.74 | 7.31 ± 0.60 |
| inosine | 90.50 ± 1.30 | 0.34 ± 0.05 | 28.45 ± 0.62 | 1.18 ± 0.18 |
| hypoxanthine | 24.50 ± 1.00 | 1.24 ± 0.21 | 5.61 ± 0.49 | 2.71 ± 0.34 |

Basis (wet muscle vs dry fibre) not stated in the table [!]; fibre values are presumably per 100 g fibre powder.

### 4b. Figs. 1-2: residual IMP and products after 1 h, mg/100 g (n = 4) [F]

"Spiked" = fibre-IMP system (Figs. 1A, 2A; nominal start 125.99 breast [M], 127.31 leg [Z]); "control" = fibre
without added IMP (Figs. 1B, 2B; start 5.99 / 7.31 natural IMP [M], Table 3).

| muscle | system | species | 80 °C a_w 0.7 / 0.8 / 0.9 | 120 °C a_w 0.7 / 0.8 / 0.9 |
|---|---|---|---|---|
| breast | spiked | IMP | 6.0 / 5.9 / 6.5 | 2.0 / 1.7 / 2.1 |
| breast | spiked | inosine | 0.41 / 0.43 / 0.35 | 0.25 / 0.14 / 0.19 |
| breast | spiked | hypoxanthine | 1.00 / 1.28 / 1.53 | 0.79 / 0.88 / 0.90 |
| breast | control | IMP | 0.69 / 0.95 / 0.71 | 0.40 / 0.38 / 0.51 |
| breast | control | inosine | 0.17 / n.d. / 0.07 | n.d. / n.d. / n.d. |
| breast | control | hypoxanthine | 0.19 / 0.19 / 1.38 | 0.10 / 0.13 / 0.12 |
| leg | spiked | IMP | 8.5 / 7.6 / 7.0 | 3.0 / 2.65 / 3.2 |
| leg | spiked | inosine | 0.41 / 0.43 / 0.35 | 0.25 / 0.14 / 0.19 |
| leg | spiked | hypoxanthine | 1.00 / 1.28 / 1.53 | 0.80 / 0.89 / 0.91 |
| leg | control | IMP | 2.07 / 2.11 / 1.55 | 1.00 / 0.70 / 0.66 |
| leg | control | inosine | 0.177 / 0.138 / 0.138 | 0.077 / 0.155 / n.d. |
| leg | control | hypoxanthine | 0.77 / 0.62 / 0.66 | 0.20 / 0.16 / 0.25 |

Printed reductions [M] (text pp. 1283-1285): breast spiked 80 °C a_w 0.8 / 0.9: 95.33 / 94.90 %; 120 °C: 98.65 /
98.37 % (the text attaches these to "a_w 0.80 and 0.90"); breast control 80 °C: 88.44 / 84.02 / 88.14 % (a_w 0.7 /
0.8 / 0.9); leg 120 °C: 97.67 / 97.95 / 97.51 % (cited as Fig. 2B, but these reproduce Fig. 2A, the spiked arm: 3.0 /
2.65 / 3.2 of 127.3 → 97.6 / 97.9 / 97.5 % [Z]). My bar reads reproduce every printed percentage within 0.2 points.

### 4c. Derived apparent first-order k = ln(C0/C1h) / 1 h, h⁻¹ [Z] (read error ±2 %)

| muscle, system | 80 °C a_w 0.7 / 0.8 / 0.9 | 120 °C a_w 0.7 / 0.8 / 0.9 | two-point "Ea" (kJ/mol) |
|---|---|---|---|
| breast spiked | 3.04 / 3.06 / 2.96 | 4.14 / 4.31 / 4.09 | 8.9 / 9.9 / 9.3 |
| leg spiked | 2.71 / 2.82 / 2.90 | 3.75 / 3.87 / 3.68 | 9.4 / 9.1 / 6.9 |
| breast control | 2.16 / 1.84 / 2.13 | 2.71 / 2.76 / 2.46 | — |
| leg control | 1.26 / 1.24 / 1.55 | 1.99 / 2.35 / 2.41 | — |

Molar balance, breast spiked 80 °C a_w 0.9 [Z]: IMP lost 119.5 mg/100 g = 0.343 mmol (as free acid, 348.2 g/mol);
hypoxanthine 1.53 mg = 0.011 mmol, inosine 0.35 mg = 0.0013 mmol → **3.6 mol %** recovered as measured products
(2-4 % across arms).

## 5. Flags and internal consistency

1. [!] **Fig. 1A and Fig. 2A inosine and hypoxanthine panels carry the same bars** (heights and error bars agree within
   0.02 mg/100 g in all 12 bars, both with the red-hatched 120 °C bars) but **different significance letters**
   (inosine 80 °C: a/a/a in Fig. 1A, a/b/c in Fig. 2A). Confidence ~90 % that one panel set was duplicated; which
   muscle the numbers belong to is unknowable. Neither panel set should be used.
2. [!] **No recovery control**: start values are nominal; IS added after extraction. The ~95 % loss at 80 °C in 1 h is
   300-600 × the aqueous hydrolysis rate (Matoba 1988; sec. 1, ans. 2) and the products close 2-4 mol %.
3. [!] **Abstract vs figures**: breast-fibre natural IMP "5.1 ± 1.20" in the abstract vs 5.99 ± 0.24 (Table 3) and
   5.99 ± 0.25 (text p. 1283). The abstract's a_w statements (breast: a_w 0.8 lowest loss at 80 °C; 0.7 = 0.9 at
   120 °C; leg 120 °C: 0.8 and 0.9 > 0.7) match the **control** panels (1B, 2B), not the spiked panels the abstract
   describes (Fig. 1A has a/a/a at 80 °C and p/q/q at 120 °C; Fig. 2A has p/q/p at 120 °C).
4. [!] Leg 120 °C reductions are cited to Fig. 2B but are Fig. 2A numbers (sec. 4b).
5. [!] "P < 0.5" printed three times for P < 0.05 (pp. 1284-1285); "1200 1C" for 120 °C (p. 1283).
6. [!] Hypoxanthine "in decreasing trend ... with a_w" (p. 1285) means it rises with a_w (the conclusion says
   decreasing "as the a_w is decreased"); true only for the 80 °C spiked arms, flat at 120 °C.
7. [!] First-order inconsistency: controls (6-7 mg/100 g start) lose a smaller fraction than spiked arms (126 mg/100 g).
8. [!] a_w nominal; measured only for dry fibre; a_w at 80 / 120 °C not considered. Buffer volumes printed as "ml" in
   the text layer, µl in the typeset page (checked at 300 dpi).
9. [!] Layer hen "Indian domesticated layer" (abstract) vs "Indian country hen" (methods); IMP salt form not stated.

## 6. Existing records in the repo that cite this paper

grep 2026-11-03 (`kavitha`, `lwt.2006.07.014`, `Modi (2007)` across `data/ docs/ src/ scripts/ tasks/`): the
registry record below, its DOI-repair note, a second DOI-repair note (Mundt & Wedzicha), `data/keys/papers.yml`
l. 3155, two raw Gemini corpus files that mis-assign this DOI to Mundt & Wedzicha (already repaired), and
`results/validation/key_value_trace_report.md` l. 133. **No reader in `src/`, `scripts/` or tests** of the record's
key values.

### 6a. `data/lit/benchmark_intake_registry.json` → `nakamura_1988_imp_ribose_release` (l. 608)

| field | repo value | paper | verdict |
|---|---|---|---|
| id | nakamura_1988_imp_ribose_release | Kavitha & Modi 2007; no Nakamura, no 1988 | **misleading id** (the doi_repair note itself says Nakamura is unsupported) |
| citation | Kavitha & Modi (2007), LWT 40:1280 | same | **matches** |
| doi | 10.1016/j.lwt.2006.07.014 | same (Crossref) | **matches** |
| kind | calibration_reference | single 1-h endpoints, no constant | **overstated** |
| chemistry_family | nucleotide_and_ribose_support | IMP, inosine, hypoxanthine; **ribose not measured** | **partly** |
| observable_panel_tags | nucleotide_support, ribose_delivery, rate_limiting_step | nucleotide data yes; no ribose, no rates | **ribose_delivery, rate_limiting_step: not in paper** |
| process_state_scope | heated_matrix | heated chicken fibre | **matches** |
| matrix_family | aqueous_and_food_matrix_comparison | fibre only; the aqueous arm is Vani 2005 [C] | **mismatch** |
| eligibility / status | benchmark_eligible_parameter / ready_for_calibration_encoding | no parameter in the paper; flags sec. 5 | **not justified** |
| what_it_supports [0] | "bounded ribose-release ceiling" (modest pentose delivery) | no ribose measured; IMP 93-98 % gone in 1 h at 80-120 °C | **not in paper** (if anything the opposite) |
| what_it_supports [1] | "inosine cleavage remains the slow step under consumer-cook conditions" | inosine ≤ 0.44 mg/100 g, hypoxanthine > inosine in every spiked arm; no rates | **not in paper**; the data do not show inosine accumulating |
| what_it_supports [2] | pool vs ribose-delivery distinction | — | **not in paper** |
| what_it_does_not_support | extrusion; sensory synergy; "Matrix-specific nucleotide retention" | it is a matrix-specific study | first two fine; **third inverted** |
| key_values.imp_to_inosine_k_h_inv_at_100c_ph6_5 | 0.053 | no 100 °C, no pH 6.5, no k (apparent 2.7-4.3 h⁻¹ [Z]) | **not in paper**; 0.053 h⁻¹ is Matoba 1988's measured IMP k at **pH 7**, 100 °C (5.3 × 10⁻² h⁻¹, t½ 13.1 h; `matoba1988_extraction.md` [C]), carried in `matoba_1988_nucleotide_hydrolysis_v1` |
| key_values.inosine_to_ribose_k_h_inv_at_100c_ph6_5 | 0.018 | no inosine rate, no ribose | **not in paper**; origin unknown |
| key_values.effective_ribose_release_fraction_at_100c_30min | 0.06 | no 100 °C, no 30 min, no ribose | **not in paper** |
| key_values.dominant_constraint | n_glycosidic_cleavage_rate_limiting | no rate comparison | **not in paper** |
| runtime_artifacts | matoba_1988_nucleotide_hydrolysis_v1 | that record's source names Matoba / Mouritsen / Yamaguchi / "Nakamura 1988", not Kavitha | **link to other papers** |
| doi_repair.basis | Kavitha "measures 5'-IMP degradation (and hence ribose release)" vs T and a_w | IMP loss yes; ribose release **not measured** | **half right**: the repair re-pointed the DOI to a paper that does not hold the key values |
| citation_repair | → Kavitha & Modi (2007) | — | **matches** |

### 6b. `benchmark_intake_registry.json` → `mundt_wedzicha_2007_biscuit_browning` doi_repair (l. 7657)

"10.1016/j.lwt.2006.07.014 resolves to Kavitha & Modi (2007) ... a different article in the same LWT issue":
identity **matches**; "same issue" **wrong** (Crossref: Kavitha 40(7):1280, Sept 2007; Mundt & Wedzicha 40(6):1078).

### 6c. `data/keys/papers.yml` l. 3155 and `results/validation/key_value_trace_report.md` l. 133

papers.yml: doi and citation **match**; paper_id inherits the misleading Nakamura id; `dossier: null` (this dossier
now exists; not edited). Trace report attributes **0.018, 0.053** to Kavitha & Modi 2007: **not in paper**.

## 7. Proposed ledger entries (aroma_other.yml; not written)

```yaml
  - id: kavitha2007_tab3_chicken_muscle_fibre_imp_inosine_hx
    dossier: kavitha2007_extraction.md
    paper: "Kavitha & Modi 2007 LWT 40:1280"
    lab: "CFTRI Mysore (Modi)"
    system: "Indian country hen breast and leg: fresh muscle and 8x-phosphate-washed freeze-dried fibre"
    medium: "raw muscle / dry fibre powder (moisture 3.25-3.50 g/100 g); no heating"
    temperature_c: []
    ph: {initial: [5.56, 5.49, 6.21, 6.43], measured_during: false, final: null}
    time_min: []
    species: [IMP, inosine, hypoxanthine]
    kind: endpoint
    location: "Table 3 (mg/100 g, mean +- SD n = 4); dossier sec. 4a"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "Values are mean ± S.D. (n = 4)"
    flags: [matrix_precursor_composition, matrix_food, basis_wet_vs_dry_not_stated, no_recovery_reported]
    note: "Chicken nucleotide pools (layer hen). Abstract prints breast-fibre IMP 5.1 +- 1.20 vs Table 3 5.99 +- 0.24. Engine has no IMP species."
  - id: kavitha2007_fig1_fig2_imp_inosine_hx_aw_temperature_1h
    dossier: kavitha2007_extraction.md
    paper: "Kavitha & Modi 2007 LWT 40:1280"
    lab: "CFTRI Mysore (Modi)"
    system: "washed freeze-dried chicken breast or leg fibre 7 g + IMP 120 mg/100 g (or none, control), 2 g in 5 mL screw-cap tube"
    medium: "0.2 M pyrophosphate pH 5.8, 20 / 110 / 620 uL per 7 g fibre for nominal a_w 0.70 / 0.80 / 0.90 (~0.24 / 0.26 / 0.33 g water/g dry)"
    temperature_c: [80, 120]
    ph: {initial: [5.8], measured_during: false, final: null}
    time_min: [60]
    species: [IMP, inosine, hypoxanthine]
    kind: endpoint
    location: "Figs. 1A/1B (breast), 2A/2B (leg), bars read at 200 dpi; text reduction percentages pp. 1283-1285; dossier sec. 4b-4c"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "subjected to heat treatment for 1 h of 80 °C in a water bath"
    flags: [single_time_point, no_recovery_control, figure_panels_duplicated_1A_2A_products, aw_nominal_not_measured, abstract_describes_control_panels, mass_balance_2_to_4_mol_pct, "also_relevant: trunk"]
    note: "IMP 93-98 % lost in 1 h at 80 and 120 C in every cell; a_w effect < 10 % in apparent k; apparent Ea 7-10 kJ/mol; 300-600x the aqueous hydrolysis rate (Matoba 1988) at 80 C. Not a rate constant: recovery artefact or enzymatic loss likely. Not usable for B54. No ribose / R5P measured."
```

## 8. Supplementary material

None referenced (no SI, appendix or supplementary table in the paper). Nothing model-relevant is missing.
