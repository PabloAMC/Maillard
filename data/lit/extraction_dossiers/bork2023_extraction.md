# Bork, Baumann, Stobernack, Rohn & Kanzler 2023 (Antioxidants 12:1788) — EXTRACTION (file stem `bork2023`)
### *"Colorants and Antioxidants Deriving from Methylglyoxal and Heterocyclic Maillard Reaction Intermediates"*
### Norfuraneol (NF) or methylglyoxal (MGO), 20 mM, alone or 1:1 with furfural (FF), maltol (MAL) or pyrrole-2-carbaldehyde (PA), in **unbuffered water adjusted to pH 5.0 ± 0.1**, sealed ampoules, **130 °C, 120 min**; reactant conversion by HPLC-DAD **relative to t = 0** (Fig. 4, 120 min only); colour (A420, CIELab), SEC, TEAC / DPPH, HRMS. **NF alone: 48 ± 4 % left after 120 min** (one endpoint; the methods name 0 / 30 / 60 / 120 min samples, but only 120 min is reported). A same-lab repeat of `bork2022` Fig. 1D, not a new laboratory and not a time course.

**Read 2026-11-05 (first reader) from `data/articles/bork2023.pdf`** (20 pp., MDPI open access). [!] The PDF carries a
second, overlaid "FOR PEER REVIEW" text layer on several pages, so `pdftotext -layout` interleaves two versions; the
prose was taken from `pdftotext -raw` and checked against 130 dpi page renders. **Fig. 4 read from a 400 dpi crop of
p. 9** (scratchpad `img/bork23_fig4b-09.png`): 7.70 px per %, read error ±1 % of scale; the text prints the same values.
Supplementary material (Figs. S1-S7, one PDF inside `bork2023_supplementary.zip`) read 2026-11-05 into **sec. 8**: no NF
time points beyond 120 min's, no vessel / O2 / final-pH detail; Fig. S4 gives the partners alone.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative statement by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me ·
**[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | **Leon Valentin Bork** (co-corresponding, l.bork@tu-berlin.de), Maximilian Baumann, Tobias Stobernack (BfR Berlin), Sascha Rohn, **Clemens Kanzler** (co-corresponding) |
| Affiliation | TU Berlin, Institute of Food Technology and Food Chemistry, Dept. of Food Chemistry and Analysis; BfR, Dept. of Chemical and Product Safety |
| Venue | *Antioxidants* **2023**, 12(9), 1788 |
| DOI | **10.3390/antiox12091788** (printed) — matches the brief |
| Dates | received 7 Sep 2023; revised 16 Sep 2023; accepted 19 Sep 2023; published 21 Sep 2023 |

**Laboratory.** Kanzler group, TU Berlin — **the same laboratory as `bork2022` / `bork2022b`** (proposed
`tu_berlin_kanzler`; not yet in `data/keys/laboratories.yml`). Ref. [41] of this paper is `bork2022`. The BfR co-author
ran HRMS only [Z, from the methods]. **For leave-one-lab-out, one laboratory with the two 2022 papers.**

## 1. What this paper answers for the model

| model need | what Bork 2023 gives | tag |
|---|---|---|
| NF decay vs time at pH 5, 130 °C, a lab not yet in the repo | **One endpoint, not a time course: NF alone 48 % (±4 % SD, n = 3) of the initial after 120 min.** The lab is new to the repo, but `bork2022` (same lab) has the time course | [M] / [F] |
| first-order k | if first order from t = 0: **k(0-120 min) = ln(1/0.48)/120 = 6.1 × 10⁻³ min⁻¹** (5.5-6.8 × 10⁻³ for ±1 SD); a two-point estimate that ASSUMES the order and ignores come-up | [Z] |
| buffered or drifting | **unbuffered**; pH adjusted to 5.0 ± 0.1 with NaOH / HCl (the paper notes the resulting "high concentration of Na⁺ and Cl⁻" in the TEAC discussion [Q]); **not re-measured** | [M] / [!] |
| O2 | sealed ampoules; headspace, degassing, chelator **not stated** | [!] |
| MGO alone | **16 % left at 120 min** ("conversion of 84 %") | [M] |
| NF + MGO | NF **7 %**, MGO **1 %** left at 120 min | [M] |
| NF with FF / PA / MAL | NF left **33 %** (FF), **55 %** (PA), **57 %** (MAL); partner left 67 % (FF), 81 % (PA), 97 % (MAL, n.s.) | [M] / [F] |
| MGO with FF / PA / MAL | MGO left 28 % (FF), 17 % (PA), 17 % (MAL); partner left 96 / 91 / 77 % (n.s. vs 100) | [M] / [F] |

## 2. Conditions [M]

| item | as printed |
|---|---|
| Pots (sec. 2.2) | equimolar **20 mmol/L**: NF 57.0 mg, MGO 40 wt% 90.1 mg (0.50 mmol [Z]), FF 48.0 mg, PA 47.6 mg, MAL 63.1 mg, each in **25 mL water**, alone or as binary mixtures |
| pH | "after adjusting the pH to 5.0 ± 0.1" [V]; initial, at room temperature; no buffer |
| Vessel / heating | sealed ampoules at 130 °C; samples at 0, 30, 60, 120 min; freezer quench −20 °C; centrifuged 8175 g |
| Headspace / O2 / metals | not stated |
| Replicates | triplicate pots; mean ± SD; two-way ANOVA + Tukey |

## 3. Assays [M]

- **NF, FF, MAL, PA:** HPLC-DAD 285 nm (Shimadzu LC-20AD, Nucleosil 120-5 C18, 0.075 % acetic acid / MeOH gradient),
  sample 1:20; **"quantified relative to the unheated samples (t = 0 min)"** [V] — a relative assay, no calibration,
  no LOD / LOQ. **Different from `bork2022`'s ion-pair method with external calibration.**
- **MGO:** OPD (50 mM, 1:1, 24 h dark RT), quinoxaline HPLC-DAD 318 nm, relative to t = 0; no LOD.
- Colour: A420 (diluted if > 0.8), CIELab; TEAC (ABTS), DPPH of ethyl acetate extracts; SEC; ESI(+)-Orbitrap HRMSⁿ.

## 4. Data

### 4.1 Fig. 4B — NF systems, 120 min, 130 °C, pH 5 (initial) [F ±1 %, matching text [M]]

| pot | NF left (%) | partner left (%) |
|---|---|---|
| **NF alone** | **47.8 (SD ≈ 3.9)** — text "48 %" | — |
| NF / FF | 32.5 (SD ≈ 1.5) — text "about 33 %" | 66.8 (SD ≈ 3) — text "67 %" |
| NF / PA | 54.6 (SD ≈ 0.5) | 81.0 (SD ≈ 0.5) — text "81 %" |
| NF / MAL | 57.0 (SD ≈ 1) | 97.3 (SD ≈ 3) |

Text: NF / PA and NF / MAL "approximately 45 %" conversion, end concentration 55 %, "significantly higher" than the 48 %
of NF alone [M]. So **PA and MAL slow NF loss** (letters c vs f), FF speeds it.

### 4.2 Fig. 4A — MGO systems, 120 min [M from text; F for SD]

MGO alone 16 % left (84 % converted); MGO / FF 28 %; MGO / PA and MGO / MAL ~17 %; MGO / NF 1 % (MGO) and 7 % (NF).
FF, PA, MAL with MGO: 96 / 91 / 77 % left (not significant).

### 4.3 Arithmetic [Z]

1. **NF alone, k(0-120 min) = 6.1 × 10⁻³ min⁻¹** if first order (5.5-6.8 × 10⁻³ at ±1 SD). Window: one interval, 0-120
   min, no come-up correction (ampoules, 25 mL fill, heat-up not stated).
2. **Same lab, same nominal pot, `bork2022` Fig. 1D: 59.2 % at 120 min (k(0-120) 4.4 × 10⁻³; full-series k 4.8-5.1 ×
   10⁻³).** The two papers differ by 11 % of scale at 120 min — about 3× the larger replicate SD — so **within-lab,
   between-paper repeatability is ~±0.1 in ln (≈ ±15 % in k)**, worse than the within-paper SD suggests. Different HPLC
   method (relative vs external calibration) is one candidate cause; unstated headspace / ampoule size is another.
3. MGO alone 16 % at 120 min vs `bork2022` 13.7 %: agree to ~2 % of scale.

## 5. What it can and cannot test in the model

**Can:** confirm the order of magnitude of NF loss in unbuffered water at pH 5 / 130 °C (k ≈ 5-6 × 10⁻³ min⁻¹, i.e. ~4×
the engine's shipped pH-free `r_nf_decay` at 130 °C, 1.42 × 10⁻³ min⁻¹, and ~19× B61's hydroxide term at pH 5, 3.2 ×
10⁻⁴ min⁻¹ — numbers and caveats in `bork2022_extraction.md` sec. 5); give the lab's own repeatability; show that FF
accelerates and PA / MAL retard NF loss at 20 mM.

**Cannot:** give a time course (only 120 min reported); give an Ea (one temperature); separate O2 (not stated); give
absolute concentrations (relative assay); the antioxidant / colour results have no engine observable.

**Comparisons asked for in the brief [Z]:**
- vs **Mikami 2017** (fed NF 1.7-1.8 mM, 0.2 M phosphate pH 8.0, 100 °C, k = 0.032 min⁻¹): at a declared 103.7 kJ/mol
  Bork 2023's 6.1 × 10⁻³ at 130 °C becomes 5.1 × 10⁻⁴ at 100 °C → Mikami is ~60× faster over 3 pH units (0.6 dex / pH);
  with 143.5 kJ/mol, ~160× (0.73 dex / pH). Sub-decade either way.
- vs **Lau 2003** (k_d ≲ 6 × 10⁻³ min⁻¹ at 2σ, whey gel, 116 °C, pH not reported): Bork 2023 moved to 116 °C is 1.3-2.0 ×
  10⁻³ min⁻¹ (143.5-103.7 kJ/mol) — **inside Lau's bound**; no tension.
- vs **Nakamura 2020 Fig. 7** (DTPA and N2 suppress NF loss at pH 8): not testable here (no chelator, no degassing arm).

## 6. Conflicts and agreements with the repo's records

1. **[!] Same-lab repeat disagrees with `bork2022` beyond the printed SDs** (48 ± 4 vs 59 ± 3 % at 120 min). Treat the
   TU Berlin NF-alone decay as **k(130 °C, pH 5) ≈ 4.4-6.1 × 10⁻³ min⁻¹ (lab range)**, and enter the two papers as one
   laboratory with a between-paper error term, not as two independent confirmations.
2. **CONTRADICTS the engine's acid-side NF sink** (same finding as `bork2022`; not independent).
3. MGO / NF 7 % / 1 % at 120 min agrees with `bork2022` Fig. 1C (NF 7.5 %, MGO 2.3 %).

## 7. Proposed ledger entry (trunk.yml; NOT written — first reader only, second reader pending)

```yaml
  - id: bork2023_fig4_norfuraneol_methylglyoxal_left_120min_ph5_130c
    dossier: bork2023_extraction.md
    paper: "Bork, Baumann, Stobernack, Rohn & Kanzler 2023 Antioxidants 12:1788"
    lab: "TU Berlin (Kanzler / Rohn)"
    system: "norfuraneol or methylglyoxal 20 mM, alone or 1:1 with furfural, maltol or pyrrole-2-carbaldehyde 20 mM"
    medium: "water, no buffer, pH 5.0 +- 0.1 set with NaOH/HCl before heating; sealed ampoules; headspace not stated"
    temperature_c: [130]
    ph: {initial: [5.0], measured_during: false, final: null}
    time_min: [120]
    species: [norfuraneol, methylglyoxal, furfural]
    kind: endpoint
    location: "Fig. 4 p. 9 and text sec. 3.2 (dossier sec. 4)"
    status: digitised_in_dossier
    quant_class: semi_quantitative
    simulable: true
    read_by: "first reader 2026-11-05; second reader pending"
    condition_quote: "heated separately or as mixtures in sealed ampoules at 130 C after adjusting the pH to 5.0"
    flags: [ph_start_only, unbuffered, o2_not_stated, single_time_point, relative_to_t0, same_lab_as_bork2022]
    note: "NF alone 48 +- 4 % at 120 min (k(0-120) 6.1e-3 if first order); bork2022 Fig. 1D same pot 59 % -> lab range 4.4-6.1e-3. MGO alone 16 %."
```

`quant_class: semi_quantitative` because the assay is relative to t = 0 (no calibration); for a fraction-remaining
observable that is adequate, but it is not an absolute concentration.

## 8. Supplementary material (read 2026-11-05, first reader)

**Source.** `data/articles/bork2023_supplementary.zip`, one member, `antioxidants-2627216-supplementary.pdf` (8 pp.,
Figs. S1-S7, no tables); archive listed before extraction (one file, no path components), extracted into a scratch
directory. Identified by content: title, authors and pots (MGO or NF with FF / MAL / PA) are this paper's; the member's
name carries the journal (Antioxidants). No name crossing for this paper. Text by `pdftotext`; figures from 80 and 300 dpi
renders; **Fig. S4 digitised at 300 dpi** (4.74 px per %, axis from its ticks; bar tops at the outline centre; read error
±1 %). Other bar charts read by eye (read error stated per row). Script: scratchpad `bork_tools/bars.py`.

### 8.1 Items asked for

| item | found? | value |
|---|---|---|
| ampoule / vial volume, fill, headspace, degassing, atmosphere, chelator / metals | **NOT found** | no methods text in the SI |
| pH after heating | **NOT found** | captions: "pH 5" (initial) only |
| **NF alone vs time; the 0 / 30 / 60 min NF samples of the methods** | **NOT found** | the SI reports nothing at 30 or 60 min for any species; NF alone appears only as 0 min colour (Figs. S1, S3) and in an HRMS spectrum (Fig. S7). The 120 min point (48 ± 4 %) remains the only NF number |
| NF temperature or pH series | **NOT found** | 130 °C, pH 5 only |

### 8.2 Every figure

| item | caption (paraphrased) | content / values | tag |
|---|---|---|---|
| Fig. S1 | browning (A420 × dilution) of (A) the unheated starting solutions, binary and single, aqueous and ethyl acetate (EE) extracts, pH 5; (B) FF, MAL, PA heated alone, 130 °C, 120 min | read ±0.02: (A) 0 min aqueous: MGO / NF 0.49, NF / FF 0.25, NF / PA 0.18, NF / MAL 0.13, **NF alone 0.20 (SD ≈ 0.06)**, MGO / FF 0.07, FF 0.06, others ≈ 0; EE ≤ ±0.05 with large SD. (B) 120 min aqueous FF 0.10, MAL 0.005, PA 0.007; EE −0.02 to −0.04 (SD ≈ 0.08-0.1, i.e. zero) | [F] |
| Fig. S2 | lightness L* of the starting solutions (0 min) and after 120 min, aqueous and EE | **[!] the graphic is a copy of Fig. S1 (A420 axis, same bars); no L* data are in the SI.** Figure / caption mismatch | [!] |
| Fig. S3 | a* and b* of the starting solutions (A, B) and of MAL and PA heated alone, 130 °C, 120 min (C, D) | read ±0.3 (A, B): b* at 0 min MGO / NF 31.5, NF alone 14 (SD ≈ 3), NF / FF 14, NF / MAL 12.5, NF / PA 10, MGO / FF 4, FF 3.3; a* MGO / NF −6.8, NF alone −2.8, NF / FF −3.0, NF / MAL −2.3, NF / PA −1.9. (C, D) MAL / PA after heating: a*, b* within ±0.1 / ±0.4 of zero | [F] |
| **Fig. S4** | conversion of FF, PA and MAL heated alone, 130 °C, pH (5), 120 min | relative concentration (vs t = 0): **FF 103 % (SD ≈ 6), PA 98 % (SD ≈ 6), MAL 111 % (SD ≈ 12)** — no loss of any heterocycle alone within SD [F ±1 %] | [F] |
| Fig. S5 | SEC (RI and 420 nm) of FF, MAL, PA heated alone | traces at half the scale of the main Fig. 4; weak RI signals ~16-18 min and ~22 min; no colour above 10 kDa | [Q] |
| Fig. S6 | antioxidant activity of MGO, FF, PA at 0 and 120 min (130 °C, pH 5): (A) TEAC aqueous, (B) DPPH of EE extracts | read ±0.03 mmol TE/L: TEAC MGO 0.13 → **1.40 (SD ≈ 0.18)**, FF 0.13 → 0.17, PA 0.07 → 0.03; DPPH MGO 0.30 → 0.38, FF 0.32 → 0.30, PA 0.10 → 0.41 (SDs 0.1-0.27) | [F] |
| Fig. S7 | HR-ESI-MS of the EE extracts of all binary pots and MGO, NF alone (solvent dimer removed in three spectra) | spectra, m/z 100-1000; not digitised | [Q] |

**For the model.** (i) Fig. S4 closes the partner baselines of sec. 4.1: FF, PA and MAL do not decay alone at
pH 5 / 130 °C / 120 min, so **FF's 33 % loss and PA's 19 % loss in the NF pots are cross-reaction with NF** [Z].
(ii) The unheated NF solution is already yellow (b* ≈ 14, A420 ≈ 0.2): the colour baseline is not zero at t = 0 — a
colour-assay caveat only. (iii) **Nothing on O2, final pH, NF at 30 / 60 min, or temperature.** The k of sec. 4.3
(6.1 × 10⁻³ min⁻¹, two-point, 0-120 min) stands unchanged, and the within-lab conflict with `bork2022` (sec. 6 item 1)
is neither resolved nor sharpened by the SI.
