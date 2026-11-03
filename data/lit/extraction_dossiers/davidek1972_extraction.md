# Davídek, Velíšek & Janíček 1972 — EXTRACTION (file stem `davidek1972`)
J. Davídek, J. Velíšek & G. Janíček, "Stability of inosinic acid, inosine and hypoxanthine in aqueous solutions" (A
Research Note), *J. Food Sci.* **37**(5):789-790 (1972). DOI: 10.1111/j.1365-2621.1972.tb02753.x (not printed on the
scan). **Crossref 2026-11-03**: same title, authors J. Davídek, J. Velíšek, G. Janíček, J. Food Sci. 37(5) 789-790,
issued 1972-09: **matches the PDF.**
[!] **Not the Nestlé Davidek.** `davidek2002/2005/2006a/2006b/2008_extraction.md` are by T. Davidek (Nestlé Research
Center, Lausanne; lab `nestle_blank`). This note is from Prague, 30 years earlier; first names are not printed
("J. Davídek"). Do not map it to `nestle_blank`.

### IMP (1 mg/mL), inosine or hypoxanthine (2 mg/mL), each alone, in buffers pH 2-12 (I ≤ 0.20, buffer species not given except acetate at pH 4), **90 °C**; IMP also at pH 3 / 4 / 5 over 60-98 °C (Ea only). IMP, inosine, hypoxanthine by ion-exchange columns + UV 248 nm. Scheme IMP → inosine (+ Pi) → hypoxanthine (+ ribose) → products, consecutive first order. Fig. 2 log k vs pH (three curves), Fig. 3 time courses at pH 4 / 6 / 8. **No table; units of k and of time not printed (inferred below); ribose and Pi drawn in the scheme but not measured.**

**Read 2026-11-03 from `data/articles/davidek1972.pdf`** (2 pp. scan with OCR layer). Prose from `pdftotext -layout`
(equations garbled; read from the render). **Fig. 1-3 read from a 300 dpi page render** (scratchpad
`ingest/davidek1972_work/p-1.png`), **Fig. 2 from a 400 dpi crop with a 0.1-log-unit grid overlaid**
(`f2grid.png`, script `scripts/grid_davidek.py`), **Fig. 3 from a 600 dpi crop** (`f3L-1.png`). Read error: Fig. 2
±0.07 on log k (curve width plus calibration); Fig. 3 ±0.03 on relative concentration and ±0.03 × 10⁴ on time.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Dept. of Food Chemistry & Analysis, **Chemical Technical University (VŠCHT / ICT), Prague 6**; no
corresponding author marked; Davídek first, Velíšek second. **Not in `data/keys/laboratories.yml`: NEW there**;
proposed id **`ict_prague_velisek`** (the id already used in `docs/reference/contrast_benchmark/lab_ids.yml` for
"ICT Prague (Novotny / Velisek)"), with the extra lab string "ICT Prague (Davidek / Velisek)". Evidence: same
department and Velíšek as co-author in both; same-lab dossier `novotny2008_extraction.md` (Novotný, Cejpek & Velíšek
2008). Merging across 36 years is the conservative choice for leave-one-laboratory-out (~60 % that methods overlap
little; the merge still errs on the safe side). Same-lab method refs [C]: Davídek & Janíček 1971 and Davídek et al.
1971 (Z. Lebensm. Unters. Forsch. 146), not on disk.

## 1. What this paper answers for the model

| model need | what Davidek 1972 gives | tag |
|---|---|---|
| IMP → inosine + Pi, k1 vs pH | Fig. 2 solid curve, 90 °C, pH ~1-13 | [F] |
| inosine → hypoxanthine + ribose, k2 vs pH (the ribose-releasing step) | Fig. 2 dashed curve, 90 °C | [F] |
| hypoxanthine → products, k3 vs pH | Fig. 2 dotted curve, 90 °C | [F] |
| Ea of k1 | 34.0 / 30.4 / 28.1 kcal/mol at pH 3 / 4 / 5, 60-98 °C (= 142.3 / 127.2 / 117.6 kJ/mol [Z]) | [M] |
| Ea of k2, k3 | **none** | [!] |
| time courses IMP / inosine / hypoxanthine | Fig. 3, pH 4 / 6 / 8, 90 °C, relative to [IMP]0 | [F] |
| ribose, R5P, Pi | **not measured** (scheme only) | [!] |
| buffer catalysis | general acid/base catalysis, H3O+ > OH⁻ > HA > A⁻ > H2O, H3O+ ~100× OH⁻; no I effect in acetate pH 4 up to 0.20 | [Q] |
| GMP | "practically the same" rate as IMP | [Q] |

**Answers.**
1. **Rate constants vs pH and T, with units** (confidence ~80 % on the unit reading). Fig. 2's axis is labelled
   "log k-7" with no unit; Fig. 3's time axis "τ" (0-5 × 10⁴) with no unit. **The only self-consistent reading is
   y = log10(k / s⁻¹) + 7 and τ in minutes** [Z]: k1 fitted to the Fig. 3 IMP decays (4.1 × 10⁻⁴, 1.2 × 10⁻⁴,
   7 × 10⁻⁵ τ⁻¹ at pH 4 / 6 / 8) gives 1.83 / 1.30 / 1.07 under that reading against 1.83 / 1.31 / 1.02 read from
   the Fig. 2 curve, and the inosine tail at pH 6 (2.9 × 10⁻⁵ min⁻¹) gives 0.68 against the dashed curve's 0.72.
   No other unit pair (s/s, min/min, h/h) comes within 1.5 log units. Independent support: Shaoul & Sporns 1987,
   extrapolated to 90 °C, agrees within 0.1 log at pH 4 and 8 (`shaoul1987_extraction.md` answer 4). At 90 °C
   [F, Z]: IMP t½ ≈ 28 h (pH 4), 55 h (pH 5), 95 h (pH 6), 180 h (pH 8); minimum near pH 10-11.
2. **Products and buffer** (confidence ~70 %). Measured: IMP, inosine, hypoxanthine (UV 248 nm after ion-exchange
   separation). Fig. 3 is the evidence for the consecutive scheme (inosine rises then falls; hypoxanthine lags).
   Ribose and Pi are written into Fig. 1 by stoichiometry only. Buffers: "buffer solutions" pH 2-12, I ≤ 0.20,
   species not named (acetate at pH 4 is mentioned); buffer catalysis was seen, so k depends on buffer identity and
   concentration that the paper does not give [!].
3. **Usable as an IMP → ribose(-5-P) source term? As the second step of IMP → inosine → ribose, at 90 °C only;
   not for R5P** (confidence ~75 %). The scheme is sequential: phosphate leaves first, so the sugar released is
   **free ribose from inosine (k2)**, not R5P in this paper's scheme (no direct IMP → hypoxanthine arrow, and at pH 4-8 Fig. 3 shows
   hypoxanthine only after inosine builds up). k2 at 90 °C [F, Z]: log10(k2/s⁻¹) + 7 = 2.23 (pH 3), 1.52 (pH 4), 1.06
   (pH 5), 0.72 (pH 6), 0.61 (pH 7-8) → t½(inosine) ≈ 3.1 h, 58 h, 167 h, 370 h, 470 h. **No Ea for k2**, so
   extrapolating to cooking temperatures needs a declared assumption. Illustration only [Z, P]: with an assumed
   Ea(k2) = 120 kJ/mol, pH 5.5, 140 °C: t½(inosine) ≈ 2 h, i.e. ~25 % of an IMP dose would reach free ribose in a
   1 h cook (k1 is ~10× faster there). The model-relevant result: **at pH 5-8 the ribose-releasing step is 3-4×
   slower than dephosphorylation** (at 90 °C; k2 < k1 above pH ~3 in Fig. 2), so IMP acts as a slow, delayed ribose source.
3b. **Against the sibling `matoba1988_extraction.md`** [Z]: (i) Matoba estimates ~18 % of IMP going directly to base + R5P at
   pH 4, 100 °C, a route this scheme omits; Fig. 3 here (hypoxanthine lagging inosine) does not exclude a branch of
   that size (~50 %). (ii) Inosine glycosidic k at pH 4, 100 °C: Matoba 1.8 × 10⁻² h⁻¹ vs this paper's 90 °C k2
   scaled with an assumed 118 kJ/mol, ~3.4 × 10⁻² h⁻¹: ×1.9. (iii) IMP k1 at 100 °C, pH 4: 9.2 h here (Ea 127.2) vs
   Matoba 8.7 h and Shaoul 8.3 h; at pH 7 ~48 h here vs 13.1 and 16.7 h. This paper is slow at pH 5-7 only.
4. **What it does NOT give**: no printed rate constants (figure only), no units, no Ea for inosine or hypoxanthine,
   no ribose/Pi data, no buffer identities or concentrations, no pH after heating, no replicates, no temperature
   above 98 °C, no identification of the hypoxanthine (k3) or ribose degradation products.

## 2. Conditions

- **Reactants** [M]: IMP 1 mg/mL (~1.9-2.9 mM, salt and hydrate form not stated) [Z]; inosine 2 mg/mL (~7.5 mM) or hypoxanthine
  2 mg/mL (~15 mM; above its room-temperature solubility of ~0.7 mg/mL, fine at 90 °C) [Z, !]. Suppliers, purity not
  given. Each compound alone (k2, k3 from separate inosine and hypoxanthine runs presumably).
- **Buffer / pH** [M]: buffers pH 2-12 (Fig. 2 curves span pH ~1-13 [F, !]), ionic strength ≤ 0.20; acetate (pH 4)
  for the ionic-strength test; buffer-ratio series at constant pH for the catalysis claim. pH reading temperature
  and post-heating pH not given.
- **Temperature / time** [M]: 90 °C for all pH series; IMP Ea from pH 3, 4, 5 over 60-98 °C. Time to τ = 5 × 10⁴
  (min, inferred: ~35 d) [Z]. Vessel, heating mode, headspace: **not stated** [!].
- **Replicates**: not stated.

## 3. Analytical methods

IMP on Dowex 1 × 8 (formate) 1 × 4 cm, 10 mL sample: water elutes the degradation products, 4 N HCl elutes IMP
(Davídek & Janíček 1971 [C]). Inosine and hypoxanthine on Dowex 50W × 4 (Na⁺) 1 × 20 cm, formic acid/formate pH
3.0. All quantified by UV at 248 nm (absolute, against the compounds' absorbance presumably; no standard curve or
LOD given). Rate constants from Eqs. 1-3 (consecutive first-order integrated forms, after Zahradník 1959 [C]); fit
method not given. Fig. 3 y axis is relative concentration ([X]/[IMP]0, "1" at top) [F].

## 4. Data

### 4.1 Fig. 2 — log10(k / s⁻¹) + 7 vs pH at 90 °C [F, ±0.07; unit reading per answer 1]

| pH | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| IMP (k1, solid) | 3.36 | 2.74 | 2.22 | 1.83 | 1.54 | 1.31 | 1.15 | 1.02 | 0.93 | 0.87 | 0.86 | 0.91 | 1.15 |
| inosine (k2, dashed) | — | 3.19 | 2.23 | 1.52 | 1.06 | 0.72 | 0.61 | 0.61 | 0.64 | 0.68 | 0.72 | 0.81 | 1.11 |
| hypoxanthine (k3, dotted) | — | 2.95 | 1.96 | 1.27 | 1.04 | 0.87 | 0.76 | 0.65 | 0.60 | 0.60 | 0.65 | 0.74 | 0.89 |

Curve ends: IMP starts at pH ~0.9 (3.40); inosine at pH ~1.7 (3.71); hypoxanthine at pH ~1.5 (3.44). The three
cross near pH 3 (~2.2). These are smooth curves; data points are not shown.

### 4.2 Fig. 3 — relative concentration vs τ (× 10⁴, min inferred) at 90 °C [F, ±0.03]

| panel | IMP | inosine | hypoxanthine |
|---|---|---|---|
| pH 4 | t½ ≈ 0.17; 0.09 at 0.59; ~0.01 at 1.05 | max 0.46 at 0.46; 0.21 at 1.03; 0.10 at 1.40; 0.05 at 1.98 | 0.41 at 0.71; 0.36 at 0.96; 0.24 at 1.94; 0.07 at 3.48; 0.05 at 4.96 |
| pH 6 | t½ ≈ 0.55; 0.29 at 0.98; 0.11 at 1.91; 0.056 at 2.49; 0.023 at 3.23 | 0.20 at 0.21; 0.40 at 0.47; 0.57 at 0.98; max 0.62 at 1.48; 0.48 at 2.97; 0.37 at 3.99; 0.27 at 5.0 | 0.025 at 0.49; 0.08 at 0.98; 0.14 at 1.48; 0.22 at 1.82; 0.27 at 2.99; 0.29 at 3.5-4.0; 0.28 at 4.50 |
| pH 8 | 0.78 at 0.30; 0.70 at 0.55; 0.49 at 1.04; 0.28 at 1.82; 0.19 at 2.42; 0.084 at 3.49; 0.02 at 4.86 | 0.22 at 0.40; 0.59 at 1.01; 0.54 at 1.99; 0.54 at 3.0; 0.48 at 4.01; 0.40 at 5.0 | 0.04 at 0.49; 0.06 at 0.99; 0.13 at 1.96; 0.22 at 3.0; 0.27 at 4.0; 0.30 at 5.0 |

First-order k1 from these points [Z]: 4.1 × 10⁻⁴ (pH 4), 1.2 × 10⁻⁴ (pH 6), 7.0 × 10⁻⁵ (pH 8) per min.

### 4.3 Text [M]

IMP Ea (60-98 °C): 34.0 (pH 3), 30.4 (pH 4), 28.1 (pH 5) kcal/mol. H3O⁺ catalysis ~100× OH⁻. No ionic-strength
effect in acetate pH 4 up to I = 0.20. GMP ≈ IMP.

## 5. Flags and internal consistency checks

- [!] **Units not printed** on either figure axis; reading derived from Fig. 2 ↔ Fig. 3 consistency (answer 1).
  Confidence ~80 %; the cross-paper agreement with Shaoul 1987 at pH 4 and 8 raises it from ~70 %.
- **Mass balance** (Fig. 3) [Z]: IMP + inosine + hypoxanthine at τ = 3 × 10⁴: 0.29 (pH 4), 0.76 (pH 6), 0.88 (pH 8);
  the deficit is hypoxanthine destroyed (k3), consistent with k3 > 0 in Fig. 2. pH 4 hypoxanthine tail
  (k ≈ 7.9 × 10⁻⁵ min⁻¹ → 1.13 on the Fig. 2 scale) vs dotted curve 1.27: within 0.15 log.
- pH 4 inosine late decay (1.5 × 10⁻⁴ min⁻¹ → 1.40) vs dashed 1.52: 0.12 log; acceptable.
- Ea trend (142 → 118 kJ/mol, pH 3 → 5) is consistent with a growing glycosidic share at low pH; Shaoul 1987's pH 5
  Ea (125.5) lies inside [Z].
- **Against Shaoul 1987** [Z]: Davidek's k1 is 2-3× lower at pH 5-7 (0.34-0.51 log) though equal at pH 4 and 8; no
  plateau at pH 4-6. Different buffers (catalysis) are the likely cause; not resolvable.
- Fig. 2 curves extend to pH ~1 and 13, outside the stated pH 2-12 [!]. Vessel, heating, headspace not stated.

## 6. Existing records in the repo that cite this paper

`grep -rIn -i -E "tb02753|davidek.{0,20}1972|Janicek|Janíček|Velíšek"` over `data/ docs/ src/ scripts/ tasks/`: **no
record cites this paper** (Velíšek hits are novotny2008 / contrast-benchmark lab ids only). The alias line
`data/keys/pdf_dossier_aliases.yml:95` ("research note, 2 pp.: stability of inosinic acid, inosine and hypoxanthine
in water") **matches**. The IMP half-life prior `matoba_1988_nucleotide_hydrolysis_v1`
(`data/lit/computational_priors.json:956`) does not cite it; for comparison, this paper at 100 °C, pH 5 (Ea 117.6
kJ/mol from 90 °C) gives t½ ≈ 19.5 h [Z] vs the prior's 12.0 h and Shaoul's 8.4 h.

## 7. Proposed ledger entries (trunk.yml; not written)

```yaml
  - id: davidek1972_fig2_imp_inosine_hx_logk_vs_ph_90c
    dossier: davidek1972_extraction.md
    paper: "Davidek, Velisek & Janicek 1972 J Food Sci 37:789"
    lab: "ICT Prague (Davidek / Velisek)"
    system: "IMP 1 mg/mL, inosine 2 mg/mL or hypoxanthine 2 mg/mL, each alone"
    medium: "aqueous buffers pH 2-12, I <= 0.20, species not stated (acetate at pH 4)"
    temperature_c: [90]
    ph: {initial: [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], measured_during: false, final: null}
    time_min: [0, 50000]
    species: [IMP, inosine, hypoxanthine]
    kind: rate_constant
    location: "Fig. 2 (log k vs pH, three curves; dossier sec. 4.1, 13 pH read)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "All studies were made at 90°C."
    flags: [units_inferred, coarse_read, ph_start_only, condition_uncertain, vessel_not_stated, no_replicates]
    note: "k1 IMP -> inosine + Pi, k2 inosine -> hypoxanthine + ribose (the ribose-releasing step), k3 hypoxanthine -> products. Axis read as log10(k/s-1)+7 with Fig. 3 time in min (dossier answer 1, ~80 %). Buffer catalysis reported, buffers not given. Not the Nestle Davidek."
  - id: davidek1972_fig3_imp_inosine_hx_time_courses_90c
    dossier: davidek1972_extraction.md
    paper: "Davidek, Velisek & Janicek 1972 J Food Sci 37:789"
    lab: "ICT Prague (Davidek / Velisek)"
    system: "IMP 1 mg/mL alone"
    medium: "aqueous buffers, species not stated"
    temperature_c: [90]
    ph: {initial: [4, 6, 8], measured_during: false, final: null}
    time_min: [0, 50000]
    species: [IMP, inosine, hypoxanthine]
    kind: time_series
    location: "Fig. 3 (relative to [IMP]0; dossier sec. 4.2, 6-10 points per curve)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: ""
    flags: [units_inferred, ph_start_only, condition_uncertain, vessel_not_stated, no_replicates]
    note: "Consecutive IMP -> inosine -> hypoxanthine; mass balance 0.29 / 0.76 / 0.88 at 3e4 min (pH 4 / 6 / 8). Ribose not measured. Time axis unit inferred (min)."
  - id: davidek1972_text_imp_ea_ph3_5
    dossier: davidek1972_extraction.md
    paper: "Davidek, Velisek & Janicek 1972 J Food Sci 37:789"
    lab: "ICT Prague (Davidek / Velisek)"
    system: "IMP 1 mg/mL alone"
    medium: "aqueous buffers pH 3, 4, 5, species not stated"
    temperature_c: [60, 98]
    ph: {initial: [3, 4, 5], measured_during: false, final: null}
    time_min: []
    species: [IMP]
    kind: rate_constant
    location: "p. 790 text: Ea 34.0 / 30.4 / 28.1 kcal/mol (no k, no temperatures listed)"
    status: text_only
    quant_class: calibrated
    simulable: false
    condition_quote: "in buffer solutions pH 3, 4 and 5 in the temperature region"
    flags: [ph_start_only, condition_not_stated:time_min, no_error_bars]
    note: "Ea of IMP loss only (142.3 / 127.2 / 117.6 kJ/mol); temperatures within 60-98 C not listed. No Ea for inosine or hypoxanthine."
```

## 8. Supplementary material

None (2-page research note). Nothing model-relevant elsewhere; the analytical method is in Davídek & Janíček 1971
and Davídek et al. 1971 (Z. Lebensm. Unters. Forsch. 146:1 and 146:83 [C]), not on disk.
