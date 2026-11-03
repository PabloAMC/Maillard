# Kuchiba, Kaizaki, Matoba & Hasegawa 1990 — EXTRACTION (file stem `kuchiba1990`)
DOI: 10.1021/jf00093a001 (not printed; Crossref 2026-11-03: "Depressing effect of salts on thermal degradation of
inosine 5'-monophosphate and guanosine 5'-monophosphate in aqueous solution", J. Agric. Food Chem. **38**(3):
**593-598**, authors Kuchiba, Kaizaki, Matoba, Hasegawa, issued 1990-03. **Matches the PDF and the brief.**)
The PDF has 6 pages: pp. 593-597 plus p. 598, whose lower part is the start of an unrelated paper (Alhir, Markakis
& Chandan, *Penicillium caseicolum* lipase) — ignore it.
**Companion paper, same group:** `matoba1988_extraction.md` (pH × T kinetics, same buffers, same HPLC; cited here
as Matoba et al. 1988).

### IMP or GMP 0.27 mM in 0.01 M acetate (pH 4.0) or 0.01 M MOPS (pH 7.0) with 0.25 M added salt (NaCl, KCl, CaCl2, MnCl2, MgCl2, CH3COONa, Na2SO4, NaH2PO4), pH readjusted, sealed Pyrex tubes, 100 °C, 0-7.5 h. First-order k ± per salt (Table I); ribose 5-phosphate (100 °C) and inosine (115 °C) k per salt at pH 4 (Table II); IMP / inosine / hypoxanthine time courses with 5 salts (Fig. 2); a CaCl2 dose pair 0.05 / 0.25 M (Fig. 1); Ni²⁺ spectra and IR of metal-nucleotide crystals (mechanism only). **Monovalent salts (NaCl, KCl) at 0.25 M change k by ≤ ~10-17 %; divalent cations cut k 1.6-2× at pH 4 and 4-20× at pH 7. One temperature; no salt dose series beyond CaCl2 0.05 / 0.25 M; no ribose / R5P from nucleotides.**

**Read 2026-11-03 from `data/articles/kuchiba1990.pdf`.** Prose from `pdftotext -layout`; **pp. 594-595 rendered at
200 dpi** (scratchpad `ingest/kuchiba1990_work/p-2.png`, `p-3.png`), Tables I-II read from them; **Fig. 1 from a
400 dpi crop** (`fig1b-2.png`, read error ±0.01 ln units) and **Fig. 2 from a 300 dpi crop** (`fig2-3.png`, ±4 µM).
Figs. 3-7 (spectra, IR, crystal structure) not digitised. Ratios [Z] by hand.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Dept. of Food Science and Nutrition, **Nara Women's University** (Kuchiba, Matoba, Hasegawa) and
Dept. of Chemistry, Osaka University (Kaizaki, the IR / coordination chemistry). **Corresponding author Teruyoshi
Matoba** (asterisk); Hasegawa senior. **Not in `data/keys/laboratories.yml`: NEW**, proposed id
**`nara_womens_matoba`**, lab string "Nara Women's Univ. (Matoba / Hasegawa)" — same group as `matoba1988`
(identical buffers, column, mobile phase, retention times). For leave-one-laboratory-out the two papers are one lab.

## 1. What this paper answers for the model

| model need | what Kuchiba 1990 gives | tag |
|---|---|---|
| salt (NaCl) correction to IMP / GMP loss | Table I: 0.25 M NaCl k vs control, pH 4 and 7, 100 °C | [M] |
| divalent-cation (Ca, Mg, Mn) effect | Table I, 0.25 M; CaCl2 0.05 M in Fig. 1 only | [M] / [F] |
| salt effect on R5P dephosphorylation and inosine N-glycosidic cleavage | Table II (pH 4; R5P 100 °C, inosine 115 °C) | [M] |
| product split under salt | Fig. 2 (IMP, inosine, Hx; 5 salts; pH 4) | [F] |
| T dependence of salt effect | **none** (100 °C only, despite "90-115 °C" in methods) | [!] |
| dose-response | CaCl2 0 / 0.05 / 0.25 M only, figure-only for 0.05 M | [!] |
| ribose / R5P release from nucleotides | **not measured** | [!] |
| mechanism | Ni²⁺ 1:1 complex (continuous variation), IR shifts → N(7)-metal-water-phosphate ring proposed | [P] |

**Answers.**
1. **Salt effect size** [Z from Table I, ratios k(salt)/k(control); ~85 % confidence the cells are read right]:

   | salt (0.25 M) | IMP pH 4 | GMP pH 4 | IMP pH 7 | GMP pH 7 |
   |---|---|---|---|---|
   | NaCl | 0.93 | 0.91 | 0.92 | 0.83 |
   | KCl | 0.94 | 0.90 | 0.85 | 0.88 |
   | CH3COONa / NaH2PO4 | 0.90 / 0.92 | 0.94 / 0.83 | NA | NA |
   | Na2SO4 | 0.78 | 0.75 | 0.67 | 0.86 |
   | CaCl2 | 0.61 | 0.58 | 0.27 | 0.17 |
   | MgCl2 | 0.53 | 0.46 | 0.10 | 0.05 |
   | MnCl2 | 0.50 | 0.32 | ND | ND |

   For the model: **NaCl at ~1.5 % w/v (0.25 M) is a ≤ 10-17 % slowing, within ~2-3 × the printed ± — set the salt
   factor to 1.0** (~80 % confidence that a real NaCl effect, if any, is < 20 %). Divalent cations at 0.25 M are
   ~100 × food levels (meat Mg²⁺ / Ca²⁺ ~1-3 mM, my rough figure, ~60 % confidence); the paper gives no low-dose
   point except CaCl2 0.05 M (Fig. 1: k ≈ 5.6-5.9 × 10⁻² h⁻¹, ratio ~0.8 [F]), so the effect at mM levels is
   unknown and plausibly small (~70 %).
2. **Which bond** [Q, ~65 % confidence in the authors' reading]: divalent cations slow the phosphoester cleavage
   of the nucleotide (less inosine formed), while they **accelerate** inosine N-glycosidic cleavage (Table II B:
   CaCl2 1.31 ×, MnCl2 2.35 ×) and free R5P dephosphorylation (Table II A: MnCl2 1.71 ×). Sulfate does the
   opposite on the glycosidic bond (inosine 0.42 ×, less Hx in Fig. 2). Reading Table II as ×10⁻² h⁻¹ (flag 1).
3. **Cross-check with Matoba 1988** [Z]: controls 7.2 / 10.9 (pH 4) and 4.8 / 7.7 (pH 7) × 10⁻² h⁻¹ vs Matoba's
   8.0 / 11.0 and 5.3 / 8.5: IMP ~10 % slower here, GMP pH 4 equal, GMP pH 7 9 % slower (different supplier, 0.27 vs
   0.4 mM). R5P control 9.1 vs Matoba 8.3 (+10 %). Between-paper reproducibility within one lab ≈ ±10 %.
4. **Inosine Ea** [Z]: inosine control 5.5 × 10⁻² h⁻¹ at 115 °C (Table II) vs 1.8 × 10⁻² at 100 °C (Matoba text)
   gives Ea ≈ 90 kJ/mol (two points from two papers; ±25 kJ/mol) — the only T information on the N-glycosidic step,
   i.e. on free-ribose release from inosine.

## 2. Conditions

- **Reactants** [M]: IMP, GMP from Kyowa Hakko Kogyo (Tokyo); other reagents Nacalai Tesque / Wako. Salt form of
  IMP / GMP not stated. **0.27 × 10⁻³ mol dm⁻³** nucleotide (Table I a; Fig. 1 starts at ln ≈ 5.60 = 270 µM [F]);
  R5P, inosine 2-3 × 10⁻³ mol dm⁻³ (Table II a). Salts **0.25 mol dm⁻³** (Table I c); CaCl2 also 0.05 M (Fig. 1).
- **Buffer / pH** [M]: 0.01 M sodium acetate pH 4.0, 0.01 M MOPS pH 7.0; solution readjusted to pH 4.0 / 7.0 after
  salt addition; "insignificant change in pH" during heating observed [Q] (no numbers). Ionic strength not stated:
  0.25 M NaCl → I ≈ 0.26; CaCl2, MgCl2, MnCl2, Na2SO4 → I ≈ 0.76 [Z] — the divalent arms confound specific
  coordination with ionic strength.
- **Vessel / heating** [M]: sealed Pyrex tube; "90-115 °C" stated, data at **100 °C** (IMP, GMP, R5P) and **115 °C**
  (inosine) only; ice-bath quench. Heating medium, headspace not stated.
- **Time** [M]: 0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 7.5 h (last points plot at ~7.9 h in Figs. 1-2 [F], flag 4).
- **Replicates** [M]: "All the experiments were repeated once" [V] — i.e. probably duplicates; ± undefined [!].

## 3. Analytical methods

HPLC-UV as Matoba 1988 (Shimadzu LC-6A, 254 nm, Cosmosil 5C18-P, methanol-0.02 M KH2PO4 15:85, 0.5 mL/min; same
retention times). Pi by Bartlett's method (printed "Bartrett (1959)" in text, 1958 in references [!]). First-order
k from semi-log slopes. Ni²⁺ visible spectra (continuous variation, total 4 × 10⁻² M), elemental analysis and KBr
FT-IR of metal-nucleotide crystals: mechanism only, no rate data.

## 4. Data

### 4a. Table I: k_obsd (×10⁻² h⁻¹) ± (undefined), 100 °C, salt 0.25 M [M]

| salt | IMP pH 4 | GMP pH 4 | IMP pH 7 | GMP pH 7 |
|---|---|---|---|---|
| control (no salt) | 7.2 ± 0.07 | 10.9 ± 0.12 | 4.8 ± 0.17 | 7.7 ± 0.87 |
| NaCl | 6.7 ± 0.21 | 9.9 ± 0.34 | 4.4 ± 0.20 | 6.4 ± 0.18 |
| KCl | 6.8 ± 0.21 | 9.8 ± 0.31 | 4.1 ± 0.08 | 6.8 ± 0.20 |
| CaCl2 | 4.4 ± 0.16 | 6.3 ± 0.21 | 1.3 ± 0.05 | 1.3 ± 0.01 |
| MnCl2 | 3.6 ± 0.21 | 3.5 ± 0.14 | ND | ND |
| MgCl2 | 3.8 ± 0.13 | 5.0 ± 0.75 | 0.5 ± 0.15 | 0.4 ± 0.08 |
| CH3COONa | 6.5 ± 0.23 | 10.2 ± 0.20 | NA | NA |
| Na2SO4 | 5.6 ± 0.13 | 8.2 ± 0.58 | 3.2 ± 0.13 | 6.6 ± 0.21 |
| NaH2PO4 | 6.6 ± 0.32 | 9.1 ± 0.07 | NA | NA |

ND = not detected, NA = not analysed (footnotes e, f). Text: k with 0.25 M CaCl2 is 61 % of control [M] (Table I:
4.4/7.2 = 0.61 [Z], consistent).

### 4b. Table II: k_obsd, pH 4, salt 0.25 M; A = ribose 5-phosphate (100 °C, Pi release), B = inosine (115 °C) [M]

| salt | A R5P | B inosine |
|---|---|---|
| control | 9.1 ± 0.20 | 5.5 ± 0.16 |
| NaCl | 10.9 ± 0.07 | 6.5 ± 0.12 |
| CaCl2 | ND | 7.2 ± 0.07 |
| MnCl2 | 15.6 ± 0.05 | 12.9 ± 0.27 |
| Na2SO4 | 10.4 ± 0.16 | 2.3 ± 0.58 |

Printed unit "h⁻¹" with no ×10⁻² (flag 1).

### 4c. Fig. 1 (IMP, pH 4, 100 °C, CaCl2), ln[IMP (µM)] at the last point (~7.5-8 h), ±0.01 [F]

control 5.07; 0.05 M CaCl2 5.16; 0.25 M CaCl2 5.24 (all from 5.60). Slopes over 7.5 h: 0.071 / 0.059 / 0.048 h⁻¹
(over 8 h: 0.066 / 0.056 / 0.045) [Z]; the 0.05 M arm has no printed k.

### 4d. Fig. 2 (IMP 270 µM, pH 4, 100 °C, salts 0.25 M), µM at the last point (~7.5 h), ±4 µM [F]

| arm | IMP | inosine | hypoxanthine | sum |
|---|---|---|---|---|
| control | 160 | 95 | 28 | 283 |
| NaCl | 165 | 87 | 26 | 278 |
| Na2SO4 | 175 | 85 | 11 | 271 |
| CaCl2 | 188 | 60 | 26 | 274 |
| MnCl2 | 201 | 43 | 21 | 265 |

[Z] ln(270/IMP)/7.5 h reproduces Table I within ~0.5 × 10⁻² h⁻¹ for every arm (control 0.070, NaCl 0.065, Na2SO4
0.057, CaCl2 0.048, MnCl2 0.039). Sums 98-105 % of 270 µM. Control proportions match Matoba 1988 Fig. 3a scaled
(inosine 86, Hx 27 per 270 µM).

## 5. Flags and internal consistency

1. [!] **Table II unit**: header "k_obsd, h⁻¹" lacks the ×10⁻² of Table I. 9.1 h⁻¹ for R5P would be ~100 × Matoba's
   8.3 × 10⁻² h⁻¹ and a 5-min half-life inconsistent with 0-7.5 h sampling: **read as ×10⁻² h⁻¹** (~95 %).
2. [!] **CaCl2 / R5P "ND"** while the text says each salt accelerated R5P hydrolysis: plausibly Ca phosphate
   precipitation masks the Pi assay (my guess, ~50 %). MnCl2 at pH 7 "ND" for both nucleotides is unexplained
   (possible Mn hydroxide / phosphate precipitation, or no measurable loss) [!].
3. [!] Methods say 90-115 °C; all nucleotide data are at 100 °C. Abstract-level conclusions are single-T.
4. [!] Last plotted times in Figs. 1-2 read at ~7.9 h vs the stated 7.5 h (axis registration or a different last
   time); slopes in sec. 4c depend on it by ~6 %.
5. [!] Divalent-salt arms carry ~3 × the ionic strength of NaCl / KCl arms; Na2SO4 (I ≈ 0.75) also slows, so part
   of the "coordination" effect may be ionic strength [Z].
6. Checked, consistent: "Mg2+ and Mn2+ were more effective than Ca2+" [Q] holds in every Table I cell; the Fig. 2
   order of remaining IMP (MnCl2 > CaCl2 > Na2SO4 > NaCl ≈ control) matches the text.
7. [!] The salt-free control k are ~10 % below Matoba 1988 for IMP (answer 3); use ±10 % as between-paper scatter.
8. [!] Reference "Bartrett" / 1958-1959 and "Kuchiba … 1988 … 1989" garbled in the reference list.

## 6. Existing records in the repo that cite this paper

Grep 2026-11-03 (`grep -rIn -iE 'kuchiba|00093a001'` over `data docs src scripts tasks`): **no data record**; none
in `computational_priors.json`, `benchmark_intake_registry.json`, `papers.yml` or any ledger. Mentions only:
- `docs/guides/LITERATURE_SEARCH_2026-10_ROUND3.md` l. 357: authors, title, DOI 10.1021/jf00093a001 — **match**
  Crossref; "same group as Matoba 1988" — **match**; "salt effect on the rate" — **match**. l. 433 lists the PDF.
- `data/keys/pdf_dossier_aliases.yml` l. 92 (triage line, JAFC 38(3)) — **match**.

## 7. Proposed ledger entries (trunk.yml; NOT written)

Same ledger as `matoba1988` (nucleotide / R5P as ribose donor). Not simulable (no IMP species).

```yaml
  # ---- kuchiba1990 ----
  - id: kuchiba1990_tab1_salt_effect_imp_gmp_k
    dossier: kuchiba1990_extraction.md
    paper: "Kuchiba, Kaizaki, Matoba & Hasegawa 1990 J Agric Food Chem 38:593"
    lab: "Nara Women's Univ. (Matoba / Hasegawa)"
    system: "IMP or GMP 0.27 mM + salt 0.25 M (NaCl, KCl, CaCl2, MnCl2, MgCl2, CH3COONa, Na2SO4, NaH2PO4)"
    medium: "0.01 M sodium acetate pH 4.0 or 0.01 M MOPS pH 7.0, pH readjusted after salt; sealed Pyrex tubes"
    temperature_c: [100]
    ph: {initial: [4.0, 7.0], measured_during: true, final: null}
    time_min: [0, 450]
    species: [IMP, GMP]
    kind: rate_constant
    location: "Table I (k_obsd x10^-2 h^-1 +/-; dossier sec. 4a, ratios sec. 1)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "Concentration of added salt was 0.25 mol dm-3"
    flags: [single_temperature, single_salt_dose, ionic_strength_confounded, mn_ph7_not_detected, plus_minus_undefined]
    note: "NaCl/KCl 0.83-0.94x control; Ca/Mg/Mn 0.32-0.61x at pH 4, 0.05-0.27x at pH 7. pH 'insignificant change' stated without numbers. Same lab as matoba1988 (controls ~10 % lower)."
  - id: kuchiba1990_tab2_salt_effect_r5p_inosine_k
    dossier: kuchiba1990_extraction.md
    paper: "Kuchiba, Kaizaki, Matoba & Hasegawa 1990 J Agric Food Chem 38:593"
    lab: "Nara Women's Univ. (Matoba / Hasegawa)"
    system: "ribose 5-phosphate or inosine 2-3 mM + salt 0.25 M (NaCl, CaCl2, MnCl2, Na2SO4)"
    medium: "0.01 M sodium acetate pH 4.0; sealed Pyrex tubes"
    temperature_c: [100, 115]
    ph: {initial: [4.0], measured_during: true, final: null}
    time_min: [0, 450]
    species: [ribose 5-phosphate, phosphate, inosine]
    kind: rate_constant
    location: "Table II (R5P at 100 C by Pi release; inosine at 115 C by inosine loss; dossier sec. 4b)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: ""
    flags: [unit_scale_missing_read_as_1e-2_per_h, ca_r5p_not_detected, single_temperature_per_species]
    note: "Inosine control 5.5e-2 h^-1 at 115 C + matoba1988 1.8e-2 at 100 C -> N-glycosidic Ea ~90 kJ/mol (dossier estimate, two papers)."
  - id: kuchiba1990_fig1_fig2_imp_inosine_hx_time_courses_salts
    dossier: kuchiba1990_extraction.md
    paper: "Kuchiba, Kaizaki, Matoba & Hasegawa 1990 J Agric Food Chem 38:593"
    lab: "Nara Women's Univ. (Matoba / Hasegawa)"
    system: "IMP 0.27 mM + CaCl2 0 / 0.05 / 0.25 M (Fig. 1); + NaCl, CaCl2, MnCl2, Na2SO4 0.25 M (Fig. 2)"
    medium: "0.01 M sodium acetate pH 4.0; sealed Pyrex tubes"
    temperature_c: [100]
    ph: {initial: [4.0], measured_during: true, final: null}
    time_min: [0, 450]
    species: [IMP, inosine, hypoxanthine]
    kind: time_series
    location: "Fig. 1 (ln IMP) and Fig. 2 (IMP, inosine, Hx); endpoints digitised dossier sec. 4c-4d"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: ""
    flags: [figure_read_error_4uM, endpoints_only_digitised, last_time_point_registration]
    note: "Only source of a CaCl2 dose point (0.05 M, ~0.8x control) and of the product split under salt."
```

## 8. Supplementary material

None referenced (1990 JAFC, no Supporting Information mentioned). Nothing model-relevant is missing. The
authors' own earlier paper on IMP in umami seasoning powder vs water activity (Kuchiba et al., J. Home Econ. Jpn.
1989, 40:677) is cited and would be the only low-moisture IMP data from this lab (not on disk; Japanese journal).
