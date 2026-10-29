# Hofmann & Schieberle 1995b — EXTRACTION (cysteamine + methylglyoxal 10 mmol/L each, 0.5 M Na/K phosphate pH 7.0, 6 °C, 18 h → 2-(1-hydroxyethyl)-4,5-dihydrothiazole (HDT); isolated HDT 0.114 mmol/L heated in water 50-100 °C, 1-60 min → 2-acetyl-2-thiazoline (AT) by SIDA; AT stability at 100 °C (water), 145 °C (phosphate pH 5) and in oil; pH 3-9 at reflux)

### THE ONLY KINETIC DATA IN THE CORPUS ON THE CYSTEAMINE ROUTE — but it is the route to **2-acetyl-2-thiazoline (AT)**, not to 2-acetylthiazole, and it starts from **cysteamine**, not cysteine. It gives the HDT → AT step (max ~11 mol % at 100 °C / 10 min, apparent Ea 57.4 kJ/mol), a steep pH law (×13 from pH 3 to 9), and AT's own fast aqueous decay (96 % in 60 min at 100 °C; > 99 % in 20 min at 145 °C, pH 5). **It measures nothing about how cysteamine is made from cysteine**: that step is only proposed in one sentence.

**Source on disk:** `data/articles/hofmann1995b.pdf` (5 pp., J. Agric. Food Chem. 1995, 43 (11),
2946-2950). Read from the `pdftotext -layout` text layer (OCR, legible); **Tables 3, 4 and 5 were
checked cell by cell against a 150 dpi render of p. 2947** (scratchpad `fork_nf/h95b_p-2.png`). The
text layer reproduces every cell correctly. Fig. 5 (Arrhenius plot) was not digitised; its printed
slope is used. Figs. 1-3 are chromatograms/mass spectra and Figs. 4, 6 structures — no data. No SI.
Not to be confused with the corpus's `hofmann1995.pdf` (`hofmann1995_extraction.md`, the AEDA
paper, JAFC 43:2187), which this paper cites as its motivation.

## 0. Identity

| field | value |
|---|---|
| Title | "Studies on the Formation and Stability of the Roast-Flavor Compound 2-Acetyl-2-thiazoline" |
| Authors | Thomas Hofmann, Peter Schieberle (corresponding) — Deutsche Forschungsanstalt für Lebensmittelchemie / Institut für Lebensmittelchemie der TU München, Garching |
| Venue | J. Agric. Food Chem. 1995, **43** (11), 2946-2950. Received 6 June 1995, accepted 23 August 1995; Advance ACS Abstracts 1 October 1995. Manuscript code printed: **JF950347Z** |
| DOI | **not printed in the PDF**. Brief's 10.1021/jf00059a032 is consistent (vol. 43 issue 11 = jf00059; p. 2946) and matches `LITERATURE_SEARCH_2026-10_ROUND3.md` row G7. **Identity verified** on title, authors, journal, volume, pages. |
| Abbreviations | HDT = 2-(1-hydroxyethyl)-4,5-dihydrothiazole (registry name in `sulfur.py` comment: "HDT/ATD"); AT = 2-acetyl-2-thiazoline; SIDA = stable isotope dilution assay ([²H]₄-AT) |

## 1. Why it matters

Gap (f): the engine predicts 2-acetylthiazole **1 750× low** in Zhang 2025b (xylose + Cys 40/40 mM,
pH 6.5 phosphate, 120 °C, 50 min; measured 114 ng/mL, semi-quantitative). `SULFUR_FLUX_MAP.md`
§5b attributes this to the missing cysteamine / thiazole route "(Hofmann & Schieberle 1995)". The
engine's only 2-acetylthiazole channel is `r_cys_actz` (Cys + MGO → ACTZ + CBX, `k_cys_actz`,
log10 k(145 °C) = −2.97 L mmol⁻¹ min⁻¹ in B60). 2-Acetyl-2-thiazoline is listed in `sulfur.py` as
a **deferred lane** ("2-acetyl-2-thiazoline's (HDT/ATD) are a separate sub-network").

What this paper actually provides:
- **Mechanism** (Fig. 6, proposed, not proven): cysteamine + MGO → aminoacetal → thioacetal → HDT
  (a *dihydro*thiazole carrying the reduced 1-hydroxyethyl side chain); HDT → AT by oxidation with
  a second MGO, or by disproportionation of HDT tautomers (giving AT + 2-(1-hydroxyethyl)thiazolidine,
  "not yet been identified"). HDT forms readily at 6 °C.
- **HDT → AT kinetics in water** (Table 3): five temperatures × six times.
- **pH law of HDT → AT** at reflux (Table 5): 1.8 → 22.8 µg from pH 3 to pH 9.
- **AT decay** (Table 4): fast in water and in hot phosphate, slow in oil.

What it does NOT provide: any experiment starting from **cysteine**; any **cysteine → cysteamine**
rate or yield (the Conclusions only propose cysteamine forms "by an oxidative decarboxylation
(Strecker reaction)" — no data); any **2-acetylthiazole** (the aromatic thiazole) measurement or the
AT → 2-acetylthiazole oxidation; any measurement of HDT itself (only AT is quantified); the
yield of HDT from cysteamine + MGO (the synthesis route via 2-hydroxypropanenitrile gives 61 %, a
different reaction).

## 2. Methods as they matter to a model

- **HDT formation pot.** "Cysteamine hydrochloride (20 mmol) and methylglyoxal (20 mmol)" in Na/K
  phosphate "2 L; 0.5 mol/L; pH 7.0", "allowed to react for 18 h at 6 °C" → **10 mmol/L each**.
  Products extracted, HDT isolated (flash chromatography + HPLC), identified by ¹H/¹³C NMR (Tables
  1, 2) and synthesis. **No yield of HDT is given.** RI DB-5: HDT 1124, AT 1111, 2-acetylthiazolidine 1152.
- **HDT → AT (Table 3).** "HDT (150 µg) was dissolved in water (10 mL; pH 7.0) and heated" →
  **0.114 mmol/L HDT** (C5H9NOS, M 131.2). Unbuffered water, pH set at start only. Vessel and
  heating mode not stated for 50-90 °C. Triplicate means, no SD printed.
- **pH series (Table 5).** "HDT (150 µg) was refluxed for 20 min in water (10 mL)", pH "adjusted with
  either HCl (1 mol/L) or NaOH (1 mol/L) prior to boiling" — **reflux ≈ 100 °C**, pH not held.
  Expt 6: + Na/K phosphate 5 mmol (0.5 mol/L) at pH 7.0; expt 7: in sunflower oil (10 mL).
- **AT decay (Table 4).** "2-Acetyl-2-thiazoline (120 µg) was refluxed in tap water (10 mL)" →
  **0.093 mmol/L AT** (C5H7NOS, M 129.2). Footnote b: phosphate 50 mL, pH 5.0, 0.5 mol/L, **145 °C,
  laboratory autoclave** (i.e. the Hofmann ribose/cysteine pot conditions); footnote c: sunflower oil,
  100 °C.
- **Quantification.** SIDA with [²H]₄-AT, ether extraction (aqueous) or high-vacuum distillation
  (oil), HRGC-MS/CI (methanol), ITD 800. Triplicates.
- **Arrhenius.** "log d[concn HDT]/dt vs 1/temp" over "10 min; temperature range of 50-90 °C";
  slope "−3.0 × 10³" → Ea **57.4 kJ/mol** (2.303 × 8.314 × 3000 = 57.4 ✔).

## 3. Tables re-typed (verified against the page render)

### Table 3 (verbatim). AT generated from HDT (150 µg in 10 mL water, pH 7.0), µg, mean of triplicates

| time (min) | 50 °C | 60 °C | 75 °C | 90 °C | 100 °C |
|---|---|---|---|---|---|
| 1 | < 0.1 | 0.2 | 0.6 | 1.8 | 7.3 |
| 5 | 0.3 | 1.0 | 2.2 | 7.2 | 15.3 |
| 10 | 0.8 | 1.9 | 3.8 | 8.9 | 16.1 |
| 20 | 1.5 | 2.6 | 5.6 | 10.0 | 10.3 |
| 30 | 1.8 | 3.2 | 7.0 | 10.3 | 6.7 |
| 60 | 2.4 | 4.1 | 10.0 | 7.6 | 1.8 |

**Same table as molar yield, mol % of HDT charged (mine; × 0.677):**

| time (min) | 50 °C | 60 °C | 75 °C | 90 °C | 100 °C |
|---|---|---|---|---|---|
| 1 | < 0.07 | 0.14 | 0.41 | 1.22 | 4.94 |
| 5 | 0.20 | 0.68 | 1.49 | 4.87 | 10.4 |
| 10 | 0.54 | 1.29 | 2.57 | 6.02 | **10.9** |
| 20 | 1.02 | 1.76 | 3.79 | 6.77 | 6.97 |
| 30 | 1.22 | 2.17 | 4.74 | 6.97 | 4.54 |
| 60 | 1.62 | 2.78 | 6.77 | 5.14 | 1.22 |

The paper's "highest yield of AT (11%)" = 16.1/150 by mass = 10.9 mol % ✔.

### Table 4 (verbatim). Degradation of AT during heating, %

| heating time (min) | degradation (%) | condition |
|---|---|---|
| 0 | 0 | 120 µg AT refluxed in 10 mL tap water |
| 5 | 7 | same |
| 10 | 31 | same |
| 20 | 59 | same |
| 60 | 96 | same |
| 20 | **> 99** | 0.5 M Na/K phosphate pH 5.0, 50 mL, **145 °C** autoclave (b) |
| 20 | **2.5** | sunflower oil, 100 °C (c) |

### Table 5 (verbatim). AT from HDT (150 µg, 10 mL water, reflux 20 min) vs pH

| expt | pH | AT (printed "mg") | mol % (mine) |
|---|---|---|---|
| 1 | 3.0 | 1.8 | 1.2 |
| 2 | 5.0 | 4.0 | 2.7 |
| 3 | 5.6 | 5.6 | 3.8 |
| 4 | 7.0 | 10.2 | 6.9 |
| 5 | 9.0 | 22.8 | **15.4** |
| 6 | 7.0 + phosphate 5 mmol | 10.3 | 7.0 |
| 7 | 7.0, in sunflower oil (10 mL) | 13.8 | 9.3 |

**The unit "mg" in the column header is a misprint for µg** (150 µg of HDT cannot give 22.8 mg; the
text's "yields of AT at pH 9.0 amounted to 15%" is 22.8/150 µg). Expt 4 (10.2 µg, 20 min) reproduces
Table 3's 100 °C / 20 min cell (10.3 µg), so reflux ≈ 100 °C and the two tables are the same protocol.

### Arithmetic (mine)

1. **Apparent Ea depends on the time window.** Linear fit of ln(AT amount) vs 1/T: 10-min amounts,
   50-90 °C → **56.8 kJ/mol** (reproduces the printed 57.4); 5-min amounts, 50-100 °C → **75.5 kJ/mol**;
   1-min amounts, 60-100 °C → **89 kJ/mol**. At 90 °C the 10-min amount is already ~87 % of that
   series' maximum, so the 10-min window compresses the high-T end. **Read 57.4 kJ/mol as a lower
   estimate; the initial-rate barrier of HDT → AT is plausibly 75-90 kJ/mol.**
2. **AT aqueous decay at 100 °C, first order.** −ln(1−d)/t = 0.015 (5 min), 0.037 (10), 0.045 (20),
   0.054 min⁻¹ (60). The early points carry the reflux heat-up lag; slope of ln(1−d) over 10-60 min
   = **0.057 min⁻¹ (t½ ≈ 12 min)**; carry **k_AT,dec(100 °C, water, pH unbuffered) ≈ 0.045-0.057 min⁻¹**.
3. **AT decay at 145 °C, pH 5 phosphate:** > 99 % in 20 min → **k ≥ 0.23 min⁻¹** (t½ ≤ 3 min). With
   (2), a two-point apparent Ea ≥ ~44 kJ/mol, confounded by medium (tap water vs 0.5 M phosphate pH 5).
4. **Consecutive A → B → C fits per temperature** (HDT → AT → products) are **not identifiable** from
   Table 3 (k₁/k₂ swap and the AT branching fraction trade off; fits tried, not reported).
   The AT maximum (≤ 11 mol %) says most HDT does **not** become AT; the paper's
   disproportionation scheme caps AT at 50 %, and the rest of the HDT fate is unmeasured.
5. **pH law (100 °C, 20 min):** AT(pH 9)/AT(pH 3) = 12.7 ("nearly 13"); log-slope ≈ 0.18 dex per pH
   unit (pH 3-9), steeper between pH 5.6 and 9 (0.22 dex/pH). Because it is a 20-min accumulated
   amount of an unstable product, it mixes formation and decay pH dependences.

## 4. Kinetic numbers the repository can use

| step / observable | quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|---|
| HDT → AT | molar yield vs t, 5 T | Table 3 molar table | mol % | 0.114 mM HDT, water pH 7.0 (unbuffered), 50-100 °C, 1-60 min | Table 3 | **measured time series** (SIDA, n = 3) |
| HDT → AT | max yield | **10.9** | mol % | 100 °C, 10 min | Table 3 / text | measured |
| HDT → AT | apparent Ea | **57.4** (printed); 75-89 from shorter windows | kJ/mol | 50-90 °C, 10-min amounts | Fig. 5 / arithmetic 1 | measured (lumped, window-dependent) |
| HDT → AT | pH dependence | 1.2 / 2.7 / 3.8 / 6.9 / 15.4 at pH 3 / 5 / 5.6 / 7 / 9 | mol % | reflux 20 min, HCl/NaOH-set pH | Table 5 | measured |
| HDT → AT | phosphate effect | none (10.2 vs 10.3 µg) | — | pH 7, 0.5 M Na/K phosphate | Table 5 | measured |
| AT decay | first-order k | **≈ 0.045-0.057** | min⁻¹ | 0.093 mM AT, tap water reflux (~100 °C) | Table 4 / arithmetic 2 | derived_from_measured |
| AT decay | first-order k | **≥ 0.23** | min⁻¹ | 0.5 M phosphate pH 5.0, 145 °C | Table 4 b | derived bound |
| AT decay | 20-min loss | 2.5 | % | sunflower oil, 100 °C | Table 4 c | measured |
| cysteine → cysteamine | rate / yield | **NOT REPORTED** | — | — | — | — |
| AT → 2-acetylthiazole | rate / yield | **NOT REPORTED** | — | — | — | — |

## 5. What this can and cannot anchor

**Can:**
1. A two-step sub-network **cysteamine + MGO → HDT → AT → fragments** with a measured HDT → AT time
   course at five temperatures and an AT decay constant at 100 °C — the formation-and-decay pair
   the deferred "2-acetyl-2-thiazoline" lane needs.
2. A **stability prior for AT**: in water and in hot phosphate AT has t½ of minutes, so any 145 °C
   or 120 °C pot measures AT only as a quasi-steady state between formation and decay.
3. A **shape hold-out**: AT from HDT rising ~13× from pH 3 to 9 at 100 °C.

**Cannot:**
- Anchor anything from cysteine. The cysteine → cysteamine decarboxylation is unmeasured here and,
  per `LITERATURE_SEARCH_2026-10_ROUND3.md` G7, nowhere in the corpus.
- Say anything about **2-acetylthiazole** (Zhang 2025b's analyte). AT is the 4,5-dihydro compound;
  its dehydrogenation to the aromatic thiazole is neither shown nor rated.
- Give HDT formation kinetics (one 6 °C / 18 h pot, no yield).

## 6. Flags

1. **Table 5's "AT (mg)" is a misprint for µg** (text and mass balance). Do not transcribe as mg.
2. **The 57.4 kJ/mol is window-dependent** (arithmetic 1); shorter windows give 75-89 kJ/mol.
3. **pH not buffered** in Tables 3 and 5 (water, HCl/NaOH-set); the HDT/AT pair is basic and the
   pH will drift. "pH 7.0" in Table 3 is a starting value.
4. **Very dilute**: 0.114 mM HDT, 0.093 mM AT; the oxidation step HDT + MGO → AT cannot occur in
   Table 3 (no MGO present), so Table 3 measures the **non-oxidative** (disproportionation/air)
   route only — a lower bound on HDT → AT in a real pot that contains MGO.
5. Heat-up for "reflux" not described; the 1- and 5-min points carry an unknown lag.
6. Registry: AT and HDT are absent from `data/keys/compounds.yml` (grep "thiazoline" finds only the
   `sulfur.py` deferral note); `ACTZ` = 2-acetylthiazole exists.

## 7. Use in the model

- **Gap (f), the 2-acetylthiazole miss (Zhang 2025b, 1 750× low): this paper does NOT close it.**
  It supports the *existence* of a cysteamine route, but to 2-acetyl-2-**thiazoline**, and it supplies
  no number on either end that Zhang's analyte needs (cysteine → cysteamine; thiazoline → thiazole).
  A model change that adds "Cys → cysteamine → … → ACTZ" would have every new rate unmeasured except
  HDT → AT and AT decay. That makes it a **declared-assumption** extension, not a measured one;
  pre-register it with Zhang's 2-acetylthiazole as the hold-out and say so. Note also the
  alternative the engine already has: `r_cys_actz` (Cys + MGO) at log10 k(145 °C) = −2.97; the
  1 750× could equally be a too-low MGO pool at pH 6.5 — Mikami 2017 (`mikami2017_extraction.md`)
  shows norfuraneol fragmenting to MGO and diacetyl at pH 8, a source the engine lacks.
- **If the AT lane is opened** (the `sulfur.py` deferral): use Table 3 (molar table) as fit rows for
  HDT → AT, Table 4's 100 °C series for k_AT,dec, the 145 °C / pH 5 bound as a constraint, and
  Table 5 as a pH-shape hold-out. Carry Ea(HDT → AT) as a band **57-90 kJ/mol**, not 57.4 ± small.
- **Caveats:** different amine (cysteamine), water not buffer, sub-mM concentrations, AT not
  2-acetylthiazole, no cysteine anywhere in the paper.
