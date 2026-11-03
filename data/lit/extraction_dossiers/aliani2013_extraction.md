# Aliani, Farmer, Kennedy, Moss & Gordon 2013 — EXTRACTION (file stem `aliani2013`)
Meat Science **94**(1):55-62, DOI 10.1016/j.meatsci.2012.11.032. Crossref 2026-11-03: "Post-slaughter changes in ATP
metabolites, reducing and phosphorylated sugars in chicken meat", authors Michel Aliani, Linda J. Farmer, James T.
Kennedy, Bruce W. Moss, Alan Gordon, vol. 94 issue 1 pp. 55-62, issued 2013-05. **Matches the PDF and the brief.**
No naming caveat (online 2012, print 2013; the stem uses the print year, as the paper's citation line does).

### CHICKEN breast (M. pectoralis major), commercial Ross 308, raw; ATP metabolites, free ribose / glucose and their phosphates by HPLC, µmol/g WET weight. Exp. 1: 15-123 min post-slaughter on the line (3 stages × 6 birds). Exp. 2: 4-200 h at 4 °C (6 birds, repeated sampling). **No fructose (not detected), no amino acids, no cysteine, no thiamine, no moisture, no cooking, no volatiles, no temperature series.**

**Read 2026-11-03 from `data/articles/aliani2013.pdf`** (8 pp., Elsevier, born-digital). Prose from
`pdftotext -layout`. **Tables 1 (p. 4) and 2 (p. 5) read from 300 dpi crops** (scratchpad
`ingest/aliani2013_work/t1.png`, `t2.png`); Table 3 and Figs. 3-5 from 200 dpi page renders (`p-5.png` to
`p-7.png`). Figs. 4 and 5 plot the Table 2 means (checked by eye, every plotted point agrees with the table to the
read error of ±0.1 µmol/g); they were not digitised separately. Sums and linear / log-linear fits by `awk` on
`ingest/aliani2013_work/t2.tsv` [Z].

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Agri-Food and Biosciences Institute (AFBI), Newforge Lane, Belfast, with Queen's University Belfast;
corresponding author **Linda J. Farmer** (AFBI); first author Aliani (since moved to Univ. of Manitoba). **Not in
`data/keys/laboratories.yml`: NEW**, proposed id **`afbi_belfast_farmer`**, lab string "AFBI Belfast (Farmer)".
No existing dossier from this group. Same group, not on disk as a dossier: Aliani & Farmer 2005a JAFC 53:6067
(repo record `aliani_2005_donor_potency_nucleotide_context` in `data/lit/benchmark_intake_registry.json` /
`data/keys/papers.yml`) and 2005b JAFC 53:6455 (ribose addition to chicken), both cited here as methods / motivation.

## 1. What this paper answers for the model

| model need | what Aliani 2013 gives | tag |
|---|---|---|
| free ribose in raw chicken muscle | 0.01-0.14 (Exp. 1), 0.02-0.38 µmol/g wet (Exp. 2, 4-200 h) | [M] |
| ribose phosphate (all isomers; R5P dominant) | 0.04-0.09 (Exp. 1), 0.26-0.44 µmol/g (Exp. 2) | [M] |
| free glucose, glucose phosphate (G6P dominant) | glucose 0.61-1.94, glucose phosphate 0.22-0.73 µmol/g | [M] |
| fructose, mannose and their phosphates | **not detected** (no LOD given for them) | [Q] |
| ATP, ADP, AMP, IMP, inosine, hypoxanthine | Tables 1, 2, recovery-corrected | [M] |
| pH | Exp. 2 per time point, 6.19-6.33 (high pHu: stressed birds) | [M] |
| moisture / dry matter | **not measured**; wet-weight basis only | [!] |
| IMP → inosine → Hx + ribose rate | 6 time points at 4 °C, one temperature; enzymatic, not thermal | [M] / [Z] |
| per-bird spread | only a pooled ANOVA SEM per row; ATP at 15 min 0.9-4.5 µmol/g (text) | [M] |
| thermal (cooking) chemistry of any kind | **none** | — |

**Answers.**
1. **Composition vs time, units, n, method** (confidence ~97 % on the transcription). All values µmol/g **wet
   weight**, recovery-corrected. Exp. 1: n = 6 different birds per stage (18 birds), stages frozen ca. 15, 33 and
   123 min post-slaughter. Exp. 2: the same 6 birds sampled at 4, 28, 55, 100, 150, 200 h, held at 4 ± 1 °C.
   Sugars by HPAEC-PAD (Dionex CarboPac PA1, Exp. 2) or post-column-derivatisation HPLC (Exp. 1); phosphates
   by alkaline-phosphatase conversion to the parent sugar (so "ribose phosphate" and "glucose phosphate" are
   sums over isomers; GC-MS says R5P and G6P dominate); nucleotides by RP-HPLC, UV 254 nm. Internal standards
   rhamnose (sugars) and cytidine (nucleotides). **Fructose was not detected** and no fructose-6-phosphate
   value is printed; the brief's list of F6P is not answered by this paper.
2. **A rate for an IMP-as-ribose-source branch?** Only a 4 °C enzymatic (post-mortem) one, and it is not a
   ribose-release rate that a cooking model could use (confidence ~85 %). [Z] from Table 2: IMP falls
   0.0138 µmol/g/h (linear, 4-200 h; 0.0172 over 4-150 h as the authors describe) or first-order
   k ≈ 3.7 × 10⁻³ h⁻¹ (t½ ≈ 190 h); Hx rises 0.0109 µmol/g/h (linear, Sig Quad NS); free ribose rises only
   0.0018 µmol/g/h, decelerating. Free ribose accumulates at ~13 % of the IMP loss rate and the total
   ribose-containing pool falls 9.03 → 7.25 µmol/g, i.e. ribose is consumed by other routes [Q]. One
   temperature, no Arrhenius, no thermal hydrolysis. A grep of `src/kinetic_core/` found no IMP / inosine /
   R5P species, so there is no branch for it to calibrate today (~85 %).
3. **Usable as precursor composition for a chicken box?** Partly (confidence ~65 % it is fit for a stub-replacing
   range, ~20 % it is representative of retail chicken). `src/cultivated_tissue_invariance.py` has only
   `beef` and `cultivated_muscle` in `BOX` (the pig construct of kim2024b feeds `cultivated_muscle`); there is
   no chicken entry. If one is added, this paper gives ribose, R5P, glucose, G6P and IMP for breast over the
   retail window. At an **assumed** 75 % moisture [!] (not in the paper): ribose 0.027-0.51 mM (retail
   28-100 h: 0.12-0.37), R5P 0.35-0.59, glucose 0.81-2.6, G6P 0.29-0.97, IMP 3.6-12.7 mM [Z]. Caveats: one
   supplier, 6 birds per experiment, high pHu (6.2-6.3, low glycogen: glucose likely at the low end of
   chicken), breast only (no thigh), no cysteine / thiamine / leucine.
4. **Ledger entries**: two (sec. 7), Table 1 and Table 2, in `aroma_other.yml`.

## 2. Conditions

- **Birds:** standard Ross 308 broilers, one commercial supplier; within an experiment same batch and farm.
  Exp. 1 in summer, Exp. 2 in winter (different occasions, different birds). Electrical stunning, commercial line.
- **Exp. 1 line temperatures:** air 15 °C for 13 min, 7-8 °C "pre chill" 30 min, 1-2 °C "chill" 2 h. Stages
  A/B/C sampled ~7, 25, 115 min, frozen ~15, 33, 123 min (dry ice, then −80 °C). Strip of breast at the keel
  edge, ~3 g portions, duplicates.
- **Exp. 2:** whole retail overwrapped birds taken at stage C, iced (0 °C, ~40 min), then **4 ± 1 °C to 200 h**.
  3 g samples from a 12-site grid (A-L) on the mid-breast, balanced across birds and times; skin / cling film
  folded back between samplings (so the cut surface is exposed to air between samplings) [!].
- **pH:** LanceFET probe at each Exp. 2 sampling; pHu = mean of 4-55 h values. Separate occasion: pH and
  internal T at 6 times 2-115 min on another 6-bird group per time (values not printed in the text; Fig. 1 only).
- **Matrix:** raw muscle; no heating; no moisture or fat reported.

## 3. Analytical methods

- **Exp. 2 extraction:** 3 g + 0.5 mL rhamnose 20 mM and cytidine 10 mM (IS) + 6 mL 0.6 M HClO₄ + 1 mL water,
  homogenise, centrifuge, KOH-neutralise, filter; one extract for sugars and nucleotides. Ion-exchange clean-up
  (Dowex 50WX4 + Marathon WBA) before sugars. Phosphates: alkaline phosphatase in 50 mM carbonate buffer pH 9.3.
- **Sugars:** HPAEC-PAD, CarboPac PA1, 180 mM NaOH isocratic-equivalent (10:90 water : 200 mM NaOH). Identity
  confirmed by spiking and GC-MS of TMS derivatives.
- **Nucleotides:** RP-HPLC, UV 254 nm (Exp. 1 Prodigy ODS-3 gradient; Exp. 2 HyPurity Aquastar, 50 mM KH₂PO₄
  pH 2.5 isocratic).
- **Recoveries [M]:** Exp. 1 glucose / ribose / G6P / R5P 93 / 90 / 64 / 93 % (CV 5 / 1 / 10 / 4 %, n = 6);
  Exp. 2 114 / 109 / 107 / 113 % (CV 3 / 3 / 39 / 36 %, n = 5). Nucleotides (both experiments) IMP / inosine /
  Hx / ATP / ADP / AMP 87 / 104 / 94 / 75 / 63 / 99 % (CV 14 / 4 / 14 / 8 / 21 / 3 %). Results are corrected
  for recovery [Q].
- **LOD:** ribose nd = below 0.012 µg (0.08 nmol) per 20 µL injection, stated as 0.01 µmol/g wet [M]. No other LOD.
- **±:** tables give one pooled **SEM** per row from ANOVA (n = 6), not per-cell SDs. Fig. 5 error bars are not
  defined in the caption [!].

## 4. Data

### 4a. Table 1 — Exp. 1, processing line, µmol/g wet weight, n = 6 birds per stage [M]

| metabolite | Stage A (~15 min) | Stage B (~33 min) | Stage C (~123 min) | SEM | Sig time | Sig lin |
|---|---|---|---|---|---|---|
| ATP | 2.42 | 1.78 | 1.33 | 0.394 | NS | NS |
| ADP | 1.01 | 0.67 | 0.80 | 0.167 | NS | NS |
| AMP | 0.98 | 0.58 | 0.83 | 0.258 | NS | NS |
| IMP | 7.78 | 9.52 | 8.29 | 0.788 | NS | NS |
| inosine | 2.30 | 2.52 | 1.55 | 0.385 | NS | NS |
| hypoxanthine | 0.14 | 0.10 | 0.29 | 0.055 | NS | P<0.05 |
| ribose | 0.01 | nd | 0.14 | 0.007 | P<0.001 | P<0.001 |
| ribose phosphate | 0.07 | 0.04 | 0.09 | 0.015 | NS | NS |
| glucose | 1.05 | 1.94 | 1.93 | 0.392 | NS | NS |
| glucose phosphate | 0.35 | 0.36 | 0.33 | 0.130 | NS | NS |
| total ribose-containing | 14.6 | 15.1 | 13.0 | 0.762 | NS | NS |
| total base-containing | 14.6 | 15.2 | 13.1 | 0.754 | NS | NS |

Times: freezing times from sec. 2.2.1; the table header gives stages only.

### 4b. Table 2 — Exp. 2, chilled storage at 4 °C, µmol/g wet weight, same 6 birds [M]

| metabolite | 4 h | 28 h | 55 h | 100 h | 150 h | 200 h | SEM | Sig time | lin | quad |
|---|---|---|---|---|---|---|---|---|---|---|
| ATP | 0.27 | 0.29 | 0.15 | 0.15 | 0.10 | 0.12 | 0.022 | <0.001 | <0.001 | <0.01 |
| ADP | 0.29 | 0.20 | 0.17 | 0.17 | 0.16 | 0.16 | 0.005 | <0.001 | <0.001 | <0.001 |
| AMP | 0.78 | 0.28 | 0.33 | 0.27 | 0.29 | 0.30 | 0.019 | <0.001 | <0.001 | <0.001 |
| IMP | 5.31 | 4.88 | 4.21 | 3.49 | 2.82 | 2.72 | 0.120 | <0.001 | <0.001 | <0.001 |
| inosine | 2.03 | 2.35 | 2.79 | 3.16 | 3.21 | 3.29 | 0.121 | <0.001 | <0.001 | <0.001 |
| hypoxanthine | 0.14 | 0.40 | 0.75 | 1.26 | 1.76 | 2.27 | 0.093 | <0.001 | <0.001 | NS |
| ribose | 0.02 | 0.09 | 0.18 | 0.28 | 0.34 | 0.38 | 0.016 | <0.001 | <0.001 | <0.001 |
| ribose phosphate | 0.32 | 0.44 | 0.37 | 0.30 | 0.26 | 0.28 | 0.028 | <0.001 | <0.001 | NS |
| glucose | 0.83 | 0.97 | 1.04 | 0.87 | 0.73 | 0.61 | 0.064 | <0.01 | <0.001 | <0.05 |
| glucose phosphate | 0.73 | 0.61 | 0.52 | 0.39 | 0.22 | 0.39 | 0.072 | <0.001 | <0.001 | <0.05 |
| pH | 6.19 | 6.30 | 6.25 | 6.28 | 6.28 | 6.33 | 0.024 | <0.05 | <0.01 | NS |
| total ribose-containing | 9.03 | 8.53 | 8.19 | 7.82 | 7.17 | 7.25 | 0.159 | <0.001 | <0.001 | <0.05 |
| total base-containing | 8.84 | 8.40 | 8.39 | 8.50 | 8.33 | 8.86 | 0.138 | <0.05 | NS | <0.01 |

"P<" prefixes dropped in the significance columns. pH at 28 h: SD 0.13 (text) [M].

### 4c. Derived [Z]

- mg/g wet (MW ribose 150.13, glucose 180.16): ribose max 0.38 µmol/g = **0.057 mg/g**; glucose 0.61-1.94 µmol/g
  = **0.11-0.35 mg/g**.
- mM in tissue water at an assumed 75 % moisture (µmol/g ÷ 0.75): see sec. 1, answer 3.
- 4 °C rates (sec. 1, answer 2). Ribose log-linear fit is meaningless (growth, not decay); use the linear slope.

### 4d. Not transcribed

Table 3: Pearson r of glucose and glucose phosphate against pHu at each Exp. 2 time (r −0.73 to −0.96, n = 6,
df = 4); Fig. 3: per-bird glucose / G6P vs pHu at 4 h and 200 h (6 points each). Relevant only as "lower pHu,
more residual sugar" [Q]; not model inputs.

## 5. Flags and internal consistency

- **Sums check** [Z]: Table 2 base-containing 8.82 / 8.40 / 8.40 / 8.50 / 8.34 / 8.86 vs printed 8.84 / 8.40 /
  8.39 / 8.50 / 8.33 / 8.86; ribose-containing 9.02 / 8.53 / 8.20 / 7.82 / 7.18 / 7.25 vs 9.03 / 8.53 / 8.19 /
  7.82 / 7.17 / 7.25. Table 1 likewise within 0.03. All within rounding: **the tables are internally consistent.**
- Text claims checked [Z]: "5-6 % of the potential ribose" free (0.38/7.25 = 5.2 %) or as phosphate (0.44/8.53
  = 5.2 %): match. Exp. 2 sum 8.6 is "60 %" of Exp. 1's 14.3: match. Ribose "20-fold" 4 → 200 h: 19×, match.
- [!] **Results vs Table 2 on glucose:** Results say only glucose was unaffected by storage; Table 2 gives
  Sig time P<0.01 (lin P<0.001) and the Discussion says glucose declines after 55 h. The table governs.
- [!] **Time labels in Exp. 1 are inconsistent**: sampled ~7 / 25 / 115 min, frozen ~15 / 33 / 123 min, and the
  Discussion speaks of 10 / 20 / 115 min. Use freezing times ±8 min.
- [!] Text says IMP fell "5.3 to 2.7 between 4 and 150 h"; Table 2 has 2.82 at 150 h and 2.72 at 200 h.
- [!] ATP "2.4 µmol/g (1.3 mg/g)": 2.42 × 507 = 1.23 mg/g as free acid; the 1.3 fits the disodium salt. Trivial.
- [!] Exp. 2 sugar-phosphate CV 36-39 % (n = 5) and recoveries > 100 % (107-114 %): phosphate values are
  imprecise; the free-sugar values (CV 3 %) are the reliable ones.
- [!] Exp. 1 and Exp. 2 used different extraction methods and differ by 40 % in total nucleotide pool (season,
  birds or method; the authors cannot separate them). Do not splice the two into one time course.
- [!] Birds were stressed (high pHu, very rapid post-mortem glycolysis [Q]); low glucose and low ATP at
  15 min are likely atypical for unstressed broilers.

## 6. Existing records in the repo that cite this paper

Grep (`aliani`, `meatsci.2012.11.032`) over `data/`, `docs/`, `src/`, `scripts/`, `tasks/`: seven places.
**Headline: the review's triple "1.48, 0.33, 0.26 mg/g of glucose, fructose, ribose" is not in Aliani 2013 in
any unit, tissue or time** (confidence ~97 %). The paper prints µmol/g only; its maximum glucose is 0.35 mg/g, its
maximum ribose 0.057 mg/g, and fructose was not detected. 0.26 does appear once in Table 2 (ribose PHOSPHATE,
150 h, µmol/g), a coincidence. Where the review's numbers come from is unresolved: glucose 1.48 and fructose
0.33 mg/g are beef-sized (Koutsidis 2008a AA-silage glucose 8.24 mmol/kg = 1.48 mg/g; Koutsidis 2008b day-1
fructose 1.81 mmol/kg = 0.33 mg/g; day-21 ribose 1.67 mmol/kg = 0.25 mg/g) [Z], so the review may have cited the
wrong paper rather than the wrong species (~40 %; no single Koutsidis column gives all three).

| record | field / claim | repo value | paper | verdict |
|---|---|---|---|---|
| `src/cultivated_tissue_invariance.py` l. 161 (beef ribose note) | "the cited 0.26 mg/g point (Aliani 2013 via Hwang 2026, 2.3 mM)" | 0.26 mg/g ribose, 2.3 mM | no 0.26 mg/g ribose; chicken free ribose ≤ 0.057 mg/g (≤ 0.51 mM at 75 %) | **mismatch**: not in the paper. Arithmetic itself right (0.26/150.13/0.75 = 2.31 mM) but the premise is false; the note should say the point is unsourced (review value not found in its cited primary) |
| `aroma_other.yml` `hwang2026_text_meat_sugar_contents_cited` | system "cited from Aliani 2013, which measured chicken" | — | Aliani measured chicken, but these values are not Aliani's | **mismatch** (implies the values are chicken values) |
| same | note "no dossier, no PDF" | — | PDF on disk, this dossier | **stale** |
| same | species [glucose, fructose, ribose, ATP] | — | fructose not detected in Aliani 2013 | **mismatch** for the Aliani attribution |
| `hwang2026_extraction.md` sec. 2 end | "the chicken numbers happen to fall inside the beef ranges ... coincidence" | — | the triple is not chicken data; real chicken glucose 0.8-2.6 mM and ribose ≤ 0.51 mM are **below** the beef ranges (2.4-15, 0.33-2.2 mM) | **mismatch** |
| same, sec. 1 | "species attribution is the review's error" | — | true that Aliani is chicken; but the numbers are not Aliani's either | **incomplete**: citation error, not only species |
| same | "Aliani 2013 itself has no PDF on disk and no dossier" | — | both now exist | **stale** |
| `bischof2023_extraction.md` l. 63 | "Aliani 2013's cited point (11 mM)" glucose | 11 mM (1.48 mg/g) | Aliani max glucose 1.94 µmol/g = 2.6 mM at 75 % | **mismatch** (label it "Hwang 2026's cited point") |
| `koutsidis2008b_extraction.md` l. 87 | "a cited 0.26 mg/g point from Aliani 2013 via Hwang 2026" | 0.26 mg/g | not in paper | **mismatch** (as above) |
| `tasks/lessons.md` l. 144 | review cites a paper whose title says CHICKEN | — | correct | **matches**; the deeper fault (values absent from the primary) strengthens the rule |
| `tasks/data_restructure_plan.md` l. 2888 | "chicken paper ... cannot source a beef range" | — | correct | **matches** |

`aliani_2005_donor_potency_nucleotide_context` (JAFC 2005) is a different paper by the same group; not checked here.

## 7. Proposed ledger entries (`data/lit/measurement_ledger/aroma_other.yml`; NOT written)

```yaml
  - id: aliani2013_tab2_chicken_breast_nucleotides_sugars_chilled_storage
    dossier: aliani2013_extraction.md
    paper: "Aliani, Farmer, Kennedy, Moss & Gordon 2013 Meat Sci 94:55"
    lab: "AFBI Belfast (Farmer)"
    system: "raw chicken breast (M. pectoralis major), Ross 308, 6 birds sampled repeatedly"
    medium: "raw muscle, whole overwrapped carcass held at 4 +/- 1 C; no heating"
    temperature_c: [4]
    ph: {initial: [6.19], measured_during: true, final: [6.33]}
    time_min: [240, 12000]
    species: [ATP, ADP, AMP, IMP, inosine, hypoxanthine, ribose, ribose phosphate, glucose, glucose phosphate]
    kind: time_series
    location: "Table 2 (umol/g wet weight, means of 6, pooled SEM; transcribed in dossier sec. 4b); Figs. 4-5 replot it"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "held at 4 +/- 1 C until the end of the experiment (200 h)"
    flags: [matrix_precursor_composition, matrix_food, moisture_not_reported, enzymatic_not_thermal]
    note: "CHICKEN breast precursor pools over the retail window (4-200 h post-slaughter). Fructose not detected. Phosphates are isomer sums (R5P, G6P dominant), CV 36-39 %. Stressed birds (pHu 6.2-6.3). Could set a chicken column of the composition BOX (none exists); IMP loss / Hx gain give 4 C enzymatic rates only."
  - id: aliani2013_tab1_chicken_breast_nucleotides_sugars_processing_line
    dossier: aliani2013_extraction.md
    paper: "Aliani, Farmer, Kennedy, Moss & Gordon 2013 Meat Sci 94:55"
    lab: "AFBI Belfast (Farmer)"
    system: "raw chicken breast, Ross 308, 3 processing stages x 6 different birds"
    medium: "raw muscle on a commercial line (air 15 C, then 7-8 C, then 1-2 C); frozen ~15, 33, 123 min"
    temperature_c: [15, 8, 2]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [15, 33, 123]
    species: [ATP, ADP, AMP, IMP, inosine, hypoxanthine, ribose, ribose phosphate, glucose, glucose phosphate]
    kind: endpoint
    location: "Table 1 (umol/g wet weight, means of 6, pooled SEM; transcribed in dossier sec. 4a)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "sample freezing times post-slaughter of ca. 15, 33 and 123 min"
    flags: [matrix_precursor_composition, matrix_food, moisture_not_reported, condition_uncertain]
    note: "Pre-rigor chicken breast; different extraction method and birds from Table 2 (total nucleotides 14.3 vs 8.6 umol/g), do not splice. Stage time labels inconsistent in the text (7/25/115 sampled vs 10/20/115 in discussion)."
```

Suggested revision of the existing entry's note (for the owner; not written): `hwang2026_text_meat_sugar_contents_cited`
→ "The review's 1.48 / 0.33 / 0.26 mg/g are not in its cited primary Aliani 2013 (chicken; µmol/g; ribose ≤ 0.057
mg/g, fructose nd; aliani2013_extraction.md). Origin unknown; unused."

## 8. Supplementary material

None referenced (no Appendix, no Supplementary Table). The per-bird data the text mentions ("not shown") are not
published. Nothing model-relevant is missing from the PDF.
